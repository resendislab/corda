#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import pytest
from corda import CORDA
from cobra import Model, Reaction, Metabolite


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
        assert len(need["EX_CORDA_0"]) == 4
        assert opt.redundancies["EX_CORDA_0"] == 1
        opt = CORDA(mod, conf, met_prod="C ->", n=1)
        need = opt.associated(["EX_CORDA_0"], conf)
        assert len(need["EX_CORDA_0"]) == 2
        assert opt.redundancies["EX_CORDA_0"] == 0


if __name__ == '__main__':
    pytest.main()
