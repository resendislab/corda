#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import pytest
from corda import reaction_confidence, test_model
from cobra import Reaction
from cobra.core.gene import GPR


class TestConf:
    @pytest.mark.parametrize(["case", "sol"], [
        ("g1 and g2 or g3", 2), ("g1 and (g2 or g3)", -1),
        ("g1 or g2 or g4 or g5", 3), ("g3 and g6", 0), ("", 0)
        ])
    def test_confidence(self, case, sol):
        vals = {"g1": -1, "g2": 1, "g3": 2, "g4": 3}
        r = Reaction("test")
        r.gpr = GPR.from_string(case)
        conf = reaction_confidence(r, vals)
        assert conf == sol

    @pytest.mark.parametrize("case", ["print() * 2", "A + B", "A ^ B"])
    def test_eval_safe(self, case):
        r = Reaction("test")
        with pytest.raises(TypeError):
            r.gpr = GPR.from_string(case)
            reaction_confidence(r, {})

    def test_none(self):
        r = Reaction("test")
        assert reaction_confidence(r, {}) == 0
        r.gpr = GPR.from_string("")
        assert reaction_confidence(r, {}) == 0


class TestMisc:
    def test_cemet(self):
        model = test_model()
        assert len(model.reactions) == 60
        assert len(model.metabolites) == 43


if __name__ == '__main__':
    pytest.main()
