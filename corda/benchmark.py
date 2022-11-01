#  benchmark.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

"""Benchmark a model reconstruction by filling in only the growth."""

from cobra.util import linear_reaction_coefficients
from cobra import Model
from time import time_ns
from .corda import CORDA
from rich import print as rprint


def benchmark(model: Model, solver: str = "glpk", **corda_args) -> CORDA:
    """Benchmark a CORDA run by building the minimal model enabling growth."""
    rprint(f"Running [green]setup[/green] for model `{model.id}`.")
    conf = {}
    for r in model.reactions:
        conf[r.id] = -1
    for r in linear_reaction_coefficients(model):
        conf[r.id] = 3
    model.solver = solver

    rprint("Running CORDA setup...", end=" ")
    start = time_ns()
    opt = CORDA(model, conf, **corda_args)
    elapsed = time_ns() - start
    rprint(f"[green]:heavy_check_mark: [{1e-9 * elapsed:.3g} s]")

    rprint("Running CORDA build...", end=" ")
    start = time_ns()
    opt.build()
    elapsed = time_ns() - start
    rprint(f"[green]:heavy_check_mark: [{1e-9 * elapsed:.3g} s]")

    rprint("Running validation on reduced model...", end=" ")
    start = time_ns()
    reduced = opt.cobra_model()
    gr = reduced.slim_optimize()
    if gr is None or gr < 1e-6:
        rprint("[red]Failed model validation![/red] Please report this at https://github.com/resendislab/corda/issues.")
    elapsed = time_ns() - start
    rprint(f"[green]:heavy_check_mark: [{1e-9 * elapsed:.3g} s]")

    return opt
