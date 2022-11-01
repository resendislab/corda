#  corda.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

"""This is a Python implementation based on the paper of Schultz et. al.

Reconstruction of Tissue-Specific Metabolic Networks Using CORDA
<http://journals.plos.org/ploscompbiol/article/
 authors?id=10.1371%2Fjournal.pcbi.1004808>
"""

from cobra import Reaction
import numpy as np
from collections import Counter
import re
from optlang.symbolics import Zero
import warnings

UPPER = 1e6    # default upper bound
CI = 1.01      # cost increase for redundancy detection


class CORDA(object):
    """The reconstruction worker.

    CORDA will perform the reconstruction from a base model and a confidence
    mapping onto the reactions. Note that the reconstruction worker is *not*
    supposed to be recycled and a new instance should be created for several
    reconstruction runs. The object supports the `print` method in order to
    give information about the reconstruction status and included reactions.

    Args:
        model (cobra model): A cobra model used as the "universe" for the
            reconstruction.
        confidence (dict): A mapping of all reaction ids in `model` to an
            integer denoting the reaction confidence. Allowed confidence values
            are -1 (absent/do not include), 0 (unknown), 1 (low confidence),
            2 (medium confidence) and 3 (high confidence). See
            `reaction_confidence` for a way to construct this dictionary.
        met_prod (Optional[list]): Additional metabolic targets that have to be
            achieved by the model. Can be a single object or list of objects.
            List elements can be given in various forms:
            (1) A string naming a metabolite in the model will ensure that the
            metabolite can be produced.
            For instance: met_prod = ["pi_c", "atp_c"]
            (2) A dictionary of metabolite -> int will define an irreversible
            reaction that must be able to carry flux.
            For instance: met_prod = {"adp_c": -1, "pi_c": -1, "atp_c": 1}
            (3) A string representation of a reversible or irreversible
            reaction that must be able to carry flux.
            For instance: met_prod = "adp_c + pi_c -> atp_c"
        n (Optional[int]): The maximum amount of redundant pathways that can be
            detected for a given high confidence reaction. For example for
            n=3 (the default) CORDA will include *at most* 3 different pathways
            to activate each of the high confidence reactions.
            **Note that setting a high n can make the build much slower!**
        penalty_factor (Optional[float]): How much more to penalize -1
            confidence in comparison to low confidence. The default is to
            penalize 100 times more.
        support (Optional[int]): The reconstruction will include an absent (-1)
            reaction if it allows inclusion of at least `support` additional
            medium confidence reactions.
        solver (Optional[str]): The LP solver to use.

    Attributes:
        conf (dict): The updated confidence. All reactions with confidence = 3
            are included.
        impossible (list): A list of reaction IDs that were chosen as high
            confidence at some point but can not carry sufficient flux.
        redundancies (dict): A dictionary {rid: n} which defines how many
            redundant pathways there are in the reconstruction to activate
            the (irreversible) reaction `rid`.
    """

    def __init__(self, model, confidence, met_prod=None, n=3,
                 penalty_factor=100, support=5):
        """Initialize parameters and model."""
        self.model = model.copy()
        self.objective = model.problem.Objective(
            model.objective.expression,
            direction=model.objective.direction)
        self.bounds = {}
        for r in self.model.reactions:
            self.bounds[r.id] = r.bounds
        self.mocks = []
        self.tol = self.model.solver.configuration.tolerances.feasibility

        # Add metabolic targets as mock reactions
        arrow_re = re.compile("<?(-+|=+)>")
        if met_prod:
            if type(met_prod) != list:
                met_prod = [met_prod]
            for i, mid in enumerate(met_prod):
                r = Reaction("EX_CORDA_" + str(i))
                r.notes["mock"] = mid
                r.upper_bound = UPPER
                self.model.add_reactions([r])
                self.mocks.append(r.id)
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

        # Map confidences from forward to backward reactions
        self.model.objective = Zero
        self.model.objective.direction = "min"
        self.conf = {}
        self.redundancies = {}
        for r in self.model.reactions:
            if r.lower_bound < -self.tol:
                r.lower_bound = -UPPER
            if r.upper_bound > self.tol:
                r.upper_bound = UPPER
            if r.id in confidence:
                if confidence[r.id] not in [-1, 0, 1, 2, 3]:
                    raise ValueError("Not a valid confidence value!")
                else:
                    self.conf[r.id] = confidence[r.id]
                    self.conf[r.reverse_id] = confidence[r.id]
                    self.redundancies[r.id] = 0
                    self.redundancies[r.reverse_id] = 0
            else:
                raise ValueError("{} missing from confidences!".format(r.id))

        self.__conf_old = self.conf.copy()
        self.built = False
        self.tflux = 1
        self.impossible = []
        self.n = n
        self.support = support
        self.pf = penalty_factor

    def __corda_objective(self, pen):
        self.model.objective.set_linear_coefficients(pen)

    def __zero_objective(self):
        self.model.objective.set_linear_coefficients(
            {v: 0 for v in self.model.variables})

    def __reduce_conf(self, conf):
        rxns = self.model.reactions
        red_conf = {r.id: max(conf[r.id], conf[r.reverse_id]) for r in rxns}

        return red_conf

    def associated(self, targets, conf=None, penalize_medium=True,
                   redundancies=True):
        """Get the associated reactions for a list of target reactions.

        `associated` calculates the smallest subset of reactions that
        are necessary so that the targets can carry flux.

        Args:
            targets (list): A list of reaction IDs used as targets.
            conf (Optional[dict]): The confidences to use for penalty
                calculation.
            penalize_medium (Optional[bool]): Whether to penalize medium
                confidence reactions.
            redundancies (Optional[bool]): Whether to get redundancies.

        Returns:
            dict: A dictionary of str->np.array mapping the reactions IDs
                of the targets to a numpy array of the associated reactions.
        """
        if conf is None:
            conf = self.conf
        m = self.model

        penalties = {}
        for r in m.reactions:
            if penalize_medium and conf[r.id] in [1, 2]:
                pen = 1
            elif conf[r.id] == -1:
                pen = self.pf
            else:
                continue
            penalties[r.forward_variable] = pen
            penalties[r.reverse_variable] = pen

        max_iter = self.n if redundancies else 1

        needed = np.array([], dtype=str)
        for vid in targets:
            va = m.variables[vid]
            self.__zero_objective()
            old_bounds = (va.lb, va.ub)
            if va.ub < self.tol:
                self.impossible.append(vid)
                self.conf[vid] = -1
                continue
            else:
                va.lb = max(self.tflux, va.lb)
                va.ub = UPPER
            has_new = True
            pen = penalties.copy()
            iteration = 0
            needed_for_v = np.array([], dtype=str)
            while has_new and iteration < max_iter:
                self.__corda_objective(pen)
                sol = self.model.solver.optimize()
                iteration += 1
                if sol != "optimal":
                    self.impossible.append(vid)
                    self.conf[vid] = -1
                    break
                sol = self.model.solver.primal_values
                need = np.array([v for v in sol if sol[v] > self.tol
                                 and conf[v] in [-1, 1, 2] and v != vid])
                new = np.in1d(need, needed_for_v, assume_unique=True,
                              invert=True)
                has_new = new.any()
                if redundancies:
                    self.redundancies[vid] += has_new
                for vi in need[new]:
                    v = m.variables[vi]
                    if v in pen:
                        pen[v] *= CI
                needed_for_v = np.unique(np.hstack([needed_for_v, need]))
            needed = np.hstack([needed, needed_for_v])
            va.lb, va.ub = old_bounds
        self.__zero_objective()

        return needed

    def build(self):
        """Construct a tissue-specific model.

        This function will initiate the build process and is the only
        computation-heavy part of CORDA.

        Args:
            None.

        Returns:
            Nothing.
        """
        if self.built:
            raise ValueError("Model has been constructed already!")

        # First iteration - find reactions required for high confidence
        include = [i for i, c in self.conf.items() if c == 3]
        need = self.associated(include)
        for a in np.unique(need):
            self.conf[a] = 3

        # Second iteration - add the best no confidence and independent medium
        # confidence
        include = [i for i, c in self.conf.items() if c in [1, 2]]
        need = self.associated(include, penalize_medium=False)
        add = [x for x in need if self.conf[x] == -1]
        count = Counter(add)
        add = [k for k in count if count[k] >= self.support]
        for a in add:
            self.conf[a] = 3

        not_included = [vid for vid, c in self.conf.items() if c == -1]
        for vid in not_included:
            v = self.model.variables[vid]
            v.ub = max(0.0, v.lb)
        self.__zero_objective()
        for v in self.model.variables:
            if self.conf[v.name] == 1 or self.conf[v.name] == 2:
                self.model.objective.set_linear_coefficients({v: 1})
                sol = self.model.solver.optimize()
                if (sol == "optimal" and
                        self.model.objective.value > self.tflux):
                    self.conf[v.name] = 3
            self.model.objective.set_linear_coefficients({v: 0})

        # Third iteration block all non-included N+M add free reactions
        for vid, co in self.conf.items():
            if co == 1 or co == 2:
                self.model.variables[vid].ub = 0.0
            elif co == 0:
                self.conf[vid] = -1
        need = self.associated([k for k in self.conf if self.conf[k] == 3],
                               penalize_medium=False, redundancies=False)
        for a in np.unique(need):
            self.conf[a] = 3

        self.impossible = np.unique(self.impossible)
        self.built = True
        self.redundancies = {vid: v for vid, v in self.redundancies.items()
                             if self.conf[vid] == 3
                             and vid not in self.impossible}

    def __str__(self):
        """Obtain basic performance infos about the reconstruction.

        Generates an acceptably nice output describing which reactions
        from each confidence level were included in the reconstruction.

        Args:
            reversible (Optional[boolean]): Whether the statistics should be
            given for the model allwoing for reversible reactions or for the
            model where each reaction is split into its forward and backward
            reactions.

        Returns:
            A formatted output string.
        """
        conf = self.__reduce_conf(self.conf)
        conf_old = self.__reduce_conf(self.__conf_old)
        old_counts = Counter([conf_old[k] for k in conf_old])

        if not self.built:
            out = "build status: not built\n" + \
                "#reactions (including mock): {}\n".format(len(conf_old)) + \
                "Reaction confidence:\n" + \
                " - unclear: {}\n".format(old_counts[0]) + \
                " - exclude: {}\n".format(old_counts[-1]) + \
                " - low and medium: {}\n".format(
                    old_counts[1] + old_counts[2]) + \
                " - high: {}\n".format(old_counts[3])
        else:
            old = np.array([conf_old[k] for k in conf_old])
            new = np.array([conf[k] for k in conf_old])
            med_inc = np.sum(((old == 1) | (old == 2)) & (new == 3))
            noc_inc = np.sum((old == -1) & (new == 3))
            free_inc = np.sum((old == 0) & (new == 3))
            high_inc = np.sum((old == 3) & (new == 3))
            out = "build status: reconstruction complete\n" + \
                "Inc. reactions: {}/{}\n".format(np.sum(new == 3), len(old)) +\
                " - unclear: {}/{}\n".format(free_inc, old_counts[0]) + \
                " - exclude: {}/{}\n".format(noc_inc, old_counts[-1]) + \
                " - low and medium: {}/{}\n".format(
                    med_inc, old_counts[1] + old_counts[2]) + \
                " - high: {}/{}\n".format(high_inc, old_counts[3])
        return out

    @property
    def included(self):
        """Show which reactions are included/excluded from the reconstruction.

        Returns:
            A dictionary {rid: True/False} indicating whether the reaction with
            ID `rid` has been included in the reconstruction.
        """
        conf = self.__reduce_conf(self.conf)
        return {rid: v == 3 for rid, v in conf.items()}

    def cobra_model(self, name=None):
        """Construct a cobra model for the reconstruction.

        Returns:
            A cobra model containing the reconstruction. The original objective
            will be conserved if it is still included in the model.
        """
        mod = self.model.copy()
        if name:
            mod.name = name

        to_remove = []
        for rxn in mod.reactions:
            co = max(self.conf[rxn.id], self.conf[rxn.reverse_id])
            if co == 3 and rxn.id not in self.mocks:
                rxn.bounds = self.bounds[rxn.id]
            else:
                to_remove.append(rxn)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mod.remove_reactions(to_remove, remove_orphans=True)
        still_valid = all(v.name in mod.variables for v in
                          self.objective.variables)
        if still_valid:
            mod.objective = self.objective
        else:
            mod.objective = Zero
        return mod
