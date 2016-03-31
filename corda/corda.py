#  corda.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

from cobra.solvers import solver_dict, get_solver_name
from cobra.manipulation.modify import convert_to_irreversible, \
    revert_to_reversible
from cobra import Model, Reaction
import numpy as np
from collections import Counter
import sys
import re
from os import devnull

TOL = 1e-6  # Tolerance to judge whether a flux is non-zero
UPPER = 1e6 # default upper bound

class CORDA(object):

    def __init__(self, model, confidence, met_prod=None, n=5,
        penalty_factor=100, support=5, solver=None, **solver_kwargs):
        """Initialize parameters and model"""
        self.model = model.copy()
        self.objective = model.objective.copy()

        # Add metabolic targets as mock reactions
        arrow_re = re.compile("<?(-+|=+)>")
        if met_prod:
            if type(met_prod) != list: met_prod = [met_prod]
            for i, mid in enumerate(met_prod):
                r = Reaction("EX_CORDA_" + str(i))
                r.notes["mock"] = mid
                r.upper_bound = UPPER
                self.model.add_reaction(r)
                if type(mid) == str:
                    if arrow_re.search(mid):
                        r.build_reaction_from_string(mid)
                    else:
                        r.add_metabolites({mid: -1})
                elif type(mid) == dict:
                    r.add_metabolites(mid)
                else:
                    raise TypeError("metabolite test not string or dictionary")
                confidence[r.id] = 3

        convert_to_irreversible(self.model)

        # Map confidences from forward to backward reactions
        self.conf = {}
        for r in self.model.reactions:
            r.objective_coefficient = 0
            r.upper_bound = UPPER
            if r.id in confidence: self.conf[r.id] = confidence[r.id]
            elif "reflection" in r.notes:
                rev = self.model.reactions.get_by_id(r.notes["reflection"])
                self.conf[r.id] = confidence[rev.id]
            else:
                raise ValueError("{} missing from confidence!".format(r.id))

        self.__conf_old = self.conf.copy()
        self.built = False
        self.tflux = 1
        self.impossible = []
        self.n = n
        self.noise = 0.5
        self.support = support
        self.pf = penalty_factor
        self.solver = solver_dict[get_solver_name() if solver is None else solver]
        self.sargs = solver_kwargs

    def __perturb(self, lp, pen):
        noise = np.random.uniform(high=self.noise, size=len(pen))

        for i, p in enumerate(pen):
            if p > 0.0:
                self.solver.change_variable_objective(lp, i, p + noise[i])

    def __quiet_solve(self, lp, os):
        old = sys.stdout
        f = open(devnull, 'w')
        sys.stdout = f
        try:
            sol = self.solver.solve_problem(lp, objective_sense=os,
                **self.sargs)
        finally:
            sys.stdout = old
            f.close()
        return sol

    def __zero_objective(self, lp, m):
        for i in range(len(m.reactions)):
            self.solver.change_variable_objective(lp, i, 0.0)

    def associated(self, targets, conf=None, penalize_medium=True):
        """Gets the associated reactions for the target reactions"""

        if conf is None: conf = self.conf
        m = self.model

        penalties = []
        for r in m.reactions:
            if penalize_medium and conf[r.id] in [1,2]:
                penalties.append(1)
            elif conf[r.id] == -1:
                penalties.append(self.pf)
            else: penalties.append(0)

        lp = self.solver.create_problem(m)

        needed = {}
        for rid in targets:
            ti = m.reactions.index(rid)
            needed[rid] = np.array([], dtype=str)
            self.__zero_objective(lp, m)
            self.solver.change_variable_bounds(lp, ti, self.tflux, UPPER)
            for _ in range(self.n):
                self.__perturb(lp, penalties)
                sol = self.__quiet_solve(lp, "minimize")
                if sol != "optimal":
                    self.impossible.append(rid)
                    self.conf[rid] = -1
                    continue
                sol = self.solver.format_solution(lp, m)
                need = [r for r in sol.x_dict if sol.x_dict[r] > TOL \
                    and conf[r] in [-1, 1, 2] and r != rid]
                need = np.hstack([needed[rid], need])
                needed[rid] = np.unique(need)
            self.solver.change_variable_bounds(lp, ti, 0.0, UPPER)
        return needed

    def build(self):
        """Constructs a tissue-specific model"""
        if self.built: raise ValueError("Model has been constructed already!")

        # First iteration - find reactions required for high confidence
        include = [r.id for r in self.model.reactions if self.conf[r.id] == 3]
        need = self.associated(include)
        add = np.unique([x for v in need.values() for x in v])
        for a in add: self.conf[a] = 3

        # Second iteration - add the best no confidence and independent medium
        # confidence
        include = [r.id for r in self.model.reactions \
            if self.conf[r.id] in [1,2]]
        need = self.associated(include, penalize_medium=False)
        add = [x for v in need.values() for x in v if self.conf[x] == -1]
        count = Counter(add)
        add = [k for k in count if count[k] >= self.support]
        for a in add: self.conf[a] = 3

        not_included = [rid for rid in self.conf if self.conf[rid] == -1]
        for rid in not_included:
            self.model.reactions.get_by_id(rid).upper_bound = 0.0
        lp = self.solver.create_problem(self.model)
        for i, r in enumerate(self.model.reactions):
            if self.conf[r.id] == 1 or self.conf[r.id] == 2:
                self.solver.change_variable_objective(lp, i, 1.0)
                sol = self.__quiet_solve(lp, "maximize")
                if sol == "optimal":
                    sol = self.solver.format_solution(lp, self.model)
                    if sol.f > self.tflux: self.conf[r.id] = 3
            self.solver.change_variable_objective(lp, i, 0.0)

        # Third iteration block all non-included N+M add free reactions
        for rid, co in self.conf.items():
            if co == 1 or co == 2:
                self.model.reactions.get_by_id(rid).upper_bound = 0.0
            elif co == 0:
                self.conf[rid] = -1
        need = self.associated([k for k in self.conf if self.conf[k] == 3],
            penalize_medium=False)
        add = np.unique([x for v in need.values() for x in v])
        for a in add: self.conf[a] = 3

        self.impossible = np.unique(self.impossible).tolist()
        self.built = True

    def __str__(self):
        old_counts = Counter(self.__conf_old.values())
        if not self.built:
            out = "build status: not built\n" + \
                "#reactions (including mock): {}\n".\
                format(len(self.__conf_old)) + \
                "Reaction confidence:\n" + \
                " - unclear: {}\n".format(old_counts[0]) + \
                " - exclude: {}\n".format(old_counts[-1]) + \
                " - low and medium: {}\n".format(old_counts[1] + \
                old_counts[2]) + \
                " - high: {}\n".format(old_counts[3])
        else:
            rids = self.conf.keys()
            old = np.array([self.__conf_old[k] for k in rids])
            new = np.array([self.conf[k] for k in rids])
            med_inc = np.sum(((old == 1) | (old == 2)) & (new == 3))
            noc_inc = np.sum((old == -1) & (new == 3))
            free_inc = np.sum((old == 0) & (new == 3))
            high_inc = np.sum((old == 3) & (new == 3))
            out = "build status: reconstruction complete\n" + \
                "Inc. reactions: {}/{}\n".format(np.sum(new == 3), len(old)) +\
                " - unclear: {}/{}\n".format(free_inc, old_counts[0]) + \
                " - exclude: {}/{}\n".format(noc_inc, old_counts[-1]) + \
                " - low and medium: {}/{}\n".format(med_inc, old_counts[1] + \
                old_counts[2]) + \
                " - high: {}/{}\n".format(high_inc, old_counts[3])
        return out


    def cobra_model(self, name, reversible=True, bound=1000):
        new_mod = Model(name)
        for rid in self.conf:
            r = self.model.reactions.get_by_id(rid)
            if self.conf[rid] == 3 and "mock" not in r.notes:
                if r not in new_mod.reactions:
                    r.upper_bound = bound
                    new_mod.add_reaction(r)
                if "reflection" in r.notes:
                    rev = self.model.reactions.get_by_id(r.notes["reflection"])
                    if rev not in new_mod.reactions:
                        rev.upper_bound = bound
                        new_mod.add_reaction(rev)

        if reversible: revert_to_reversible(new_mod)
        still_valid = True
        for r in self.objective:
            still_valid &= (self.conf[r.id] == 3)
        if still_valid: new_mod.change_objective(self.objective)
        return new_mod
