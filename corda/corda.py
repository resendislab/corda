#  corda.py
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

from cobra.solvers import solver_dict, get_solver_name
from cobra.manipulation.modify import convert_to_irreversible
from cobra import Model, Metabolite, Reaction
import numpy as np
from collections import Counter
from .util import safe_revert_reversible

TOL = 1e-6  # Tolerance to judge whether a flux is close to zero

class CORDA(object):

    def __init__(self, model, confidence, met_prod=None, flux_target=1, n=10, 
        noise=1, penalty=1e4, support=2, solver=None, **solver_kwargs):
        """Initialize parameters and model"""
        self.model = model.copy()
        self.objective = model.objective.copy()
        convert_to_irreversible(self.model)
        
        # Map confidences from forward to backward reactions
        self.conf = {}
        for r in self.model.reactions:
            r.objective_coefficient = 0
            r.upper_bound = len(self.model.reactions)*penalty + 1.0
            if r.id in confidence: self.conf[r.id] = confidence[r.id]
            elif "reflection" in r.notes: 
                rev = self.model.reactions.get_by_id(r.notes["reflection"])
                if rev.id not in confidence: 
                    ValueError("{} missing from confidence!".format(rev.id))
                self.conf[r.id] = confidence[rev.id]
            else:
                ValueError("{} missing from confidence!".format(r.id))
        
        self.tflux = 1
        self.n = n
        self.noise = noise
        self.support = support
        self.penalty = penalty
        self.solver = solver_dict[get_solver_name() if solver is None else solver]
        self.corda_solver = solver_dict[get_solver_name() if solver is None else solver]
        self.sargs = solver_kwargs
        
        self.m_pen = Metabolite("penalty")
        self.r_pen = Reaction("EX_penalty")
        self.r_pen.notes["mock"] = self.m_pen
        self.r_pen.add_metabolites({self.m_pen: -1})
        self.r_pen.objective_coefficient = 1.0
        self.r_pen.upper_bound = 1e16
        
        if met_prod:
            for mid in met_prod:
                r = Reaction("EX_CORDA_" + mid)
                r.notes["mock"] = mid
                r.add_metabolites({self.model.metabolites.get_by_id(mid): -1})
                self.model.add_reaction(r)
                self.conf[r.id] = 3

    def __perturb(self, lp, m):
        noise = np.random.uniform(high=self.noise, size=len(m.reactions))
        pen_id = m.metabolites.index("penalty")
        
        for i, r in enumerate(m.reactions):
            if self.m_pen in r.metabolites and r.id != "EX_penalty":
                self.corda_solver.change_coefficient(lp, pen_id, 
                    i, r.metabolites[self.m_pen] + noise[i])

    def associated(self, targets, conf=None):
        """Gets the associated reactions for the target reactions"""
        
        if conf == None: conf = self.conf
        
        tidx = ((i, r.id) for i, r in enumerate(self.model.reactions) \
            if r.id in targets)
        m = self.model.copy()
        m.add_reaction(self.r_pen)
        
        for r in m.reactions:
            if r.id == "EX_penalty": continue
            if conf[r.id] == 2 or conf[r.id] == 1: 
                r.add_metabolites({self.m_pen: np.sqrt(self.penalty)})
            elif conf[r.id] == -1: 
                r.add_metabolites({self.m_pen: self.penalty})
        
        lp = self.solver.create_problem(self.model)
        corda_lp = self.corda_solver.create_problem(m)
        
        needed = {}
        for ti, rid in tidx:
            self.solver.change_variable_objective(lp, ti, 1.0)
            sol = self.solver.solve_problem(lp, objective_sense="maximize", **self.sargs)
            sol = self.solver.format_solution(lp, self.model)
            if(sol.f < TOL): 
                raise ValueError("Reaction {} can never carry flux!".format(rid))
            self.solver.change_variable_objective(lp, ti, 0.0)
            
            upper = m.reactions.get_by_id(rid).upper_bound
            self.corda_solver.change_variable_bounds(corda_lp, ti, 
                self.tflux, upper)
            needed[rid] = np.array([])
            for i in range(self.n):
                self.__perturb(corda_lp, m)
                sol = self.corda_solver.solve_problem(corda_lp, 
                    objective_sense="minimize", **self.sargs)
                sol = self.corda_solver.format_solution(corda_lp, m)
                need = [r for r in sol.x_dict if sol.x_dict[r] > TOL \
                    and r != "EX_penalty" and conf[r] in [-1, 1, 2]]
                need = np.hstack([needed[rid], need])
                needed[rid] = np.unique(need)
            
            self.corda_solver.change_variable_bounds(corda_lp, ti, 0.0, upper)
            
        return needed

    def build(self):
        """Constructs a tissue-specific model"""
        
        # First iteration - find reactions required for high confidence
        include = [r.id for r in self.model.reactions if self.conf[r.id] == 3]
        need = self.associated(include)
        add = np.unique([x for v in need.values() for x in v])
        for a in add: self.conf[a] = 3
        
        # Second iteration - add the best no confidence and independent medium confidence
        include = [r.id for r in self.model.reactions if self.conf[r.id] == 1 \
            or self.conf[r.id] == 2]
        need = self.associated(include)
        add = [x for v in need.values() for x in v if self.conf[x] == -1]
        count = Counter(add)
        add = [k for k in count.keys() if count[k] >= self.support]
        for a in add: self.conf[a] = 3
        
        not_included = [rid for rid in self.conf if self.conf[rid] == -1]
        for rid in not_included: 
            self.model.reactions.get_by_id(rid).upper_bound = 0.0
        lp = self.solver.create_problem(self.model)
        for i, r in enumerate(self.model.reactions):
            if self.conf[r.id] == 1 or self.conf[r.id] == 2:
                self.solver.change_variable_objective(lp, i, 1.0)
                sol = self.solver.solve_problem(lp, objective_sense="maximize", 
                    **self.sargs)
                sol = self.solver.format_solution(lp, self.model)
                if sol.f > TOL: self.conf[r.id] == 3
            self.solver.change_variable_objective(lp, i, 0.0)
        
        # Third iteration block all non-included N+M add free reactions
        for rid, co in self.conf.items():
            if co == 1 or co == 2:
                self.model.reactions.get_by_id(rid).upper_bound = 0.0
            if co == 0:
                self.conf[rid] = -1
        need = self.associated([k for k in self.conf if self.conf[k] == 3])
        add = np.unique([x for v in need.values() for x in v])
        for a in add: self.conf[a] = 3

    def cobra_model(self, name, reversible=True, bound = 1000):
        new_mod = Model(name)
        for rid in self.conf:
            r = self.model.reactions.get_by_id(rid)
            if self.conf[rid] == 3 and "mock" not in r.notes:
                r.upper_bound = bound
                new_mod.add_reaction(r)
        
        if reversible: safe_revert_reversible(new_mod)
        new_mod.change_objective(self.objective)
        return new_mod
