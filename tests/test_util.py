#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import pytest
from corda import reaction_confidence, test_model
from cobra import Model, Reaction, Metabolite
from cobra.manipulation import convert_to_irreversible, revert_to_reversible


class TestConf:
    @pytest.mark.parametrize("case", [
        ("g1 and g2 or g3", 2), ("g1 and (g2 or g3)", -1),
        ("g1 or g2 or g4 or g5", 3), ("g3 and g6", 0), ("", 0)
        ])
    def test_confidence(self, case):
        vals = {"g1": -1, "g2": 1, "g3": 2, "g4": 3}
        conf = reaction_confidence(case[0], vals)
        assert conf == case[1]

    @pytest.mark.parametrize("case", ["print()", "A + B", "A ^ B"])
    def test_eval_safe(self, case):
        with pytest.raises(TypeError):
            reaction_confidence(case, {})

    def test_none(self):
        assert reaction_confidence("  ", {}) == 0


class TestMisc:

    def test_cemet(self):
        model = test_model()
        assert len(model.reactions) == 60
        assert len(model.metabolites) == 43


if __name__ == '__main__':
    pytest.main()
