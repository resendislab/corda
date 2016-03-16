#  tests.py
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

import unittest
from . import *
from cobra import Model, Reaction, Metabolite
from cobra.manipulation import convert_to_irreversible, revert_to_reversible

class TestConf(unittest.TestCase):
    
    def test_confidence(self):
        vals = {"g1": -1, "g2": 1, "g3": 2, "g4": 3}
        cases = [("g1 and g2 or g3", 2), ("g1 and (g2 or g3)", -1),
            ("g1 or g2 or g4 or g5", 3), ("g3 and g6", 0)]
        
        for rule, res in cases:
            conf = reaction_confidence(rule, vals) 
            self.assertEqual(conf, res)

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
        
class TestCORDA(unittest.TestCase):
    def setUp(self):
        A = Metabolite("A")
        B = Metabolite("B")
        r1 = Reaction("EX_A")
        r1.add_metabolites({A: -1})
        r1.lower_bound = -1
        r1.upper_bound= 1
        r2 = Reaction("EX_B")
        r2.add_metabolites({B: -1})
        r2.lower_bound = -1
        r2.upper_bound= 1
        r3 = Reaction("conversion")
        r3.add_metabolites({A: -1, B: 1})
        model = Model("test model")
        model.add_reactions([r1, r2, r3])
        conf = {"EX_A": 1, "EX_B": 3, "conversion": 2}
        
        self.model = model
        self.conf = conf
    
    
if __name__ == '__main__':
    unittest.main()
