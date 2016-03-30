#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import unittest
from corda import *
from cobra import Model, Reaction, Metabolite
from cobra.manipulation import convert_to_irreversible, revert_to_reversible
from cobra.io import read_sbml_model

class TestConf(unittest.TestCase):

    def test_confidence(self):
        vals = {"g1": -1, "g2": 1, "g3": 2, "g4": 3}
        cases = [("g1 and g2 or g3", 2), ("g1 and (g2 or g3)", -1),
            ("g1 or g2 or g4 or g5", 3), ("g3 and g6", 0), ("", 0)]

        for rule, res in cases:
            conf = reaction_confidence(rule, vals)
            self.assertEqual(conf, res)

    def test_eval_safe(self):
        cases = ["print()", "A + B", "A ^ B"]
        for ca in cases:
            with self.assertRaises(TypeError): reaction_confidence(ca, {})

    def test_none(self):
        self.assertEqual(reaction_confidence("  ", {}), 0)

class TestRevert(unittest.TestCase):

    def setUp(self):
        self.model = Model("test model")
        A = Metabolite("A")
        r = Reaction("r")
        r.add_metabolites({A: -1})
        r.lower_bound = -1000
        r.upper_bound = 1000
        self.model.add_reaction(r)
        convert_to_irreversible(self.model)
        self.model.remove_reactions(["r"])

    def test_remove_breaks(self):
        self.assertRaises(KeyError, revert_to_reversible, self.model)

    def test_safe_revert(self):
        safe_revert_reversible(self.model)
        self.assertEqual(len(self.model.reactions), 1)
        self.assertEqual(self.model.reactions[0].id, "r_reverse")
        self.assertFalse(self.model.reactions[0].reversibility)

class TestCORDAsimple(unittest.TestCase):
    def setUp(self):
        A = Metabolite("A")
        B = Metabolite("B")
        C = Metabolite("C")
        r1 = Reaction("r1")
        r1.add_metabolites({A: -1, C: 1})
        r2 = Reaction("r2")
        r2.add_metabolites({B: -1, C: 1})
        r3 = Reaction("EX_A")
        r3.add_metabolites({A: 1})
        r4 = Reaction("EX_B")
        r4.add_metabolites({B: 1})
        r5 = Reaction("EX_C")
        r5.add_metabolites({C: -1})
        model = Model("test model")
        model.add_reactions([r1, r2, r3, r4, r5])
        conf = {"r1": 1, "r2": -1, "EX_A": 1, "EX_B": 1, "EX_C": 1}

        self.model = model
        self.conf = conf
        self.opt = CORDA(self.model, self.conf, met_prod=["C"])

    def test_mock_add(self):
        r = self.opt.model.reactions.get_by_id("EX_CORDA_0")
        self.assertTrue("mock" in r.notes)
        opt = CORDA(self.model, self.conf, met_prod={"C": -1})
        r = opt.model.reactions.get_by_id("EX_CORDA_0")
        self.assertTrue("mock" in r.notes)
        opt = CORDA(self.model, self.conf, met_prod="C ->")
        r = opt.model.reactions.get_by_id("EX_CORDA_0")
        self.assertTrue("mock" in r.notes)
        with self.assertRaises(TypeError):
            CORDA(self.model, self.conf, met_prod=[["C"]])

    def test_conf_check(self):
        conf = self.conf.copy()
        del conf["EX_A"]
        self.assertRaises(ValueError, CORDA, self.model, conf)

    def test_performance_metrics(self):
        self.assertTrue("not built" in str(self.opt))

    def test_impossible_req(self):
        model = self.model.copy()
        D = Metabolite("D")
        model.add_metabolites([D])
        opt = CORDA(model, self.conf, met_prod=["D"])
        need = opt.associated(["EX_CORDA_0"])
        self.assertTrue(len(need["EX_CORDA_0"]) == 0)
        self.assertTrue("EX_CORDA_0" in opt.impossible)

    def test_association_works(self):
        need = self.opt.associated(["EX_CORDA_0"])
        self.assertEqual(list(need["EX_CORDA_0"]), ["EX_A", "r1"])

    def test_redundancy_works(self):
        conf = self.opt.conf.copy()
        conf["r2"] = 2
        need = self.opt.associated(["EX_CORDA_0"], conf)
        self.assertEqual(len(need["EX_CORDA_0"]), 4)
        self.opt.n = 1
        need = self.opt.associated(["EX_CORDA_0"], conf)
        self.assertEqual(len(need["EX_CORDA_0"]), 2)

class TestCORDAlarge(unittest.TestCase):
    def setUp(self):
        model = read_sbml_model("data/cemet.xml")
        conf = {}
        for i, r in enumerate(model.reactions):
            if i % 2 == 0: conf[r.id] = -1
            else: conf[r.id] = 2
        conf["r60"] = 3
        self.model = model
        self.conf = conf

    def test_init_works(self):
        opt = CORDA(self.model, self.conf)
        self.assertTrue(len(opt.conf) > 0)

    def test_association_work(self):
        opt = CORDA(self.model, self.conf)
        need = opt.associated(["r60"])
        self.assertTrue(len(need["r60"]) > 0)

    def test_conf_vals(self):
        conf = self.conf.copy()
        for r in self.model.reactions: conf[r.id] = 1
        conf["r60"] = 3
        conf["r42"] = 0
        conf["r12"] = 1
        opt = CORDA(self.model, conf)
        opt.build()
        include = [c for c in opt.conf if opt.conf[c] == 3]
        self.assertTrue(len(include) > 3)

    def test_performance_metrics(self):
        opt = CORDA(self.model, self.conf)
        opt.build()
        self.assertTrue("reconstruction complete" in str(opt))

    def test_build_works(self):
        opt = CORDA(self.model, self.conf)
        opt.build()
        include = [c for c in opt.conf if opt.conf[c] == 3]
        self.assertTrue(len(include) > 3)
        rec = opt.cobra_model("reconstruction")
        rec.optimize()
        self.assertTrue(rec.solution.f > 1)

if __name__ == '__main__':
    unittest.main()
