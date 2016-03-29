#  util.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import re
from cobra.manipulation import revert_to_reversible
from cobra.core.Gene import parse_gpr
from ast import Name, And, Or, BoolOp, Expression

def format_gid(gid):
    return re.sub(r"\.\d*", "", gid)

def safe_eval_gpr(expr, conf_genes):
    """safely evaluate confidence of a gene-reaction rule"""
    if isinstance(expr, Expression):
        return safe_eval_gpr(expr.body, conf_genes)
    elif isinstance(expr, Name):
        fgid = format_gid(expr.id)
        if fgid not in conf_genes: return 0
        return conf_genes[fgid]
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            return max(safe_eval_gpr(i, conf_genes) for i in expr.values)
        elif isinstance(op, And):
            return min(safe_eval_gpr(i, conf_genes) for i in expr.values)
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    elif expr is None:
        return 0
    else:
        raise TypeError("unsupported operation  " + repr(expr))

def reaction_confidence(rule, conf_genes):
    """Calculates the confidence for the reaction based on a gene-reaction
        rule"""

    ast_rule, _ = parse_gpr(rule)
    return safe_eval_gpr(ast_rule, conf_genes)


def safe_revert_reversible(model):
    maybe_rev = (r for r in model.reactions if "reflection" in r.notes)

    for r in maybe_rev:
        try:
            model.reactions.get_by_id(r.notes["reflection"])
        except KeyError:
            r.notes.pop("reflection")

    revert_to_reversible(model)
