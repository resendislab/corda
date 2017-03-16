#  tests.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import pytest
from corda import CORDA, test_model
from cobra.util.solver import solvers

if "scipy" in solvers:
    del solvers["scipy"]


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


class TestCORDAlarge:

    def test_init_works(self, large):
        opt = CORDA(*large)
        assert len(opt.conf) > 0

    def test_association_work(self, large):
        opt = CORDA(*large)
        need = opt.associated(["r60"])
        assert len(need) > 0

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
        assert "/60" in str(opt)

    def test_build_works(self, large):
        opt = CORDA(*large)
        opt.build()
        include = [c for c in opt.conf if opt.conf[c] == 3]
        assert len(include) > 3
        rec = opt.cobra_model("reconstruction")
        sol = rec.optimize()
        assert sol.f > 1

    @pytest.mark.parametrize("solver", solvers)
    def test_benchmark_large(self, large, benchmark, solver):
        mod, conf = large
        mod.solver = solver

        def _():
            opt = CORDA(mod, conf)
            opt.build()
        benchmark(_)


if __name__ == '__main__':
    pytest.main()
