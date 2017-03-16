#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import pytest
from corda import CORDA, reaction_confidence, test_model
from cobra import Model, Reaction, Metabolite
from cobra.manipulation import convert_to_irreversible, revert_to_reversible


@pytest.fixture
def model():
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
    mod = Model("test model")
    mod.add_reactions([r1, r2, r3, r4, r5])
    conf = {"r1": 1, "r2": -1, "EX_A": 1, "EX_B": 1, "EX_C": 1}

    return (mod, conf)


@pytest.fixture
def large():
    mod = test_model()
    conf = {}
    for i, r in enumerate(mod.reactions):
        if i % 2 == 0:
            conf[r.id] = -1
        else:
            conf[r.id] = 2
    conf["r60"] = 3

    return (mod, conf)


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

    def test_remove_breaks(self):
        model = Model("test model")
        A = Metabolite("A")
        r = Reaction("r")
        r.add_metabolites({A: -1})
        r.lower_bound = -1000
        r.upper_bound = 1000
        model.add_reaction(r)
        convert_to_irreversible(model)
        model.remove_reactions(["r"])
        with pytest.raises(KeyError):
            revert_to_reversible(model)

    def test_cemet(self):
        model = test_model()
        assert len(model.reactions) == 60
        assert len(model.metabolites) == 43


class TestCORDAsimple:

    def test_mock_add(self, model):
        mod, conf = model
        opt = CORDA(mod, conf, met_prod={"C": -1})
        r = opt.model.reactions.get_by_id("EX_CORDA_0")
        assert "mock" in r.notes
        opt = CORDA(mod, conf, met_prod="C ->")
        r = opt.model.reactions.get_by_id("EX_CORDA_0")
        assert "mock" in r.notes
        with pytest.raises(TypeError):
            CORDA(mod, conf, met_prod=[["C"]])

    def test_conf_check(self, model):
        mod, conf = model
        co = conf.copy()
        del co["EX_A"]
        with pytest.raises(ValueError):
            CORDA(mod, co)

    def test_valid_conf(self, model):
        mod, conf = model
        co = conf.copy()
        co["EX_A"] = 4
        with pytest.raises(ValueError):
            CORDA(mod, co)

    def test_performance_metrics(self, model):
        opt = CORDA(model[0], model[1])
        assert "not built" in str(opt)

    def test_impossible_req(self, model):
        mod, conf = model
        D = Metabolite("D")
        mod.add_metabolites([D])
        opt = CORDA(mod, conf, met_prod=["D"])
        need = opt.associated(["EX_CORDA_0"])
        assert len(need["EX_CORDA_0"]) == 0
        assert "EX_CORDA_0" in opt.impossible

    def test_association_works(self, model):
        mod, conf = model
        opt = CORDA(mod, conf, met_prod="C ->")
        need = opt.associated(["EX_CORDA_0"])
        solutions = (["EX_A", "r1"], ["EX_B", "r2"])
        assert list(need["EX_CORDA_0"]) in solutions

    def test_redundancy_works(self, model):
        mod, conf = model
        conf["r2"] = 2
        opt = CORDA(mod, conf, met_prod="C ->")
        need = opt.associated(["EX_CORDA_0"], conf)
        assert len(need["EX_CORDA_0"]) >= 2
        opt.n = 1
        need = opt.associated(["EX_CORDA_0"], conf)
        assert len(need["EX_CORDA_0"]) == 2


class TestCORDAlarge:

    def test_init_works(self, large):
        opt = CORDA(*large)
        assert len(opt.conf) > 0

    def test_association_work(self, large):
        opt = CORDA(*large)
        need = opt.associated(["r60"])
        assert len(need["r60"]) > 0

    def test_conf_vals(self, large):
        mod, conf = large
        for r in mod.reactions:
            conf[r.id] = 1
        conf["r60"] = 3
        conf["r42"] = 0
        conf["r12"] = 1
        opt = CORDA(mod, conf)
        opt.build()
        include = [c for c in opt.conf if opt.conf[c] == 3]
        assert len(include) > 3

    def test_performance_metrics(self, large):
        opt = CORDA(*large)
        opt.build()
        assert "reconstruction complete" in str(opt)
        assert "/60" in opt.info()

    def test_build_works(self, large):
        opt = CORDA(*large)
        opt.build()
        include = [c for c in opt.conf if opt.conf[c] == 3]
        assert len(include) > 3
        rec = opt.cobra_model("reconstruction")
        sol = rec.optimize()
        assert sol.f > 1

if __name__ == '__main__':
    pytest.main()
