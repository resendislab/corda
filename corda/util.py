#  util.py
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

import re
from cobra.core.Gene import parse_gpr
from ast import Name, And, Or, BoolOp, Expression


def format_gid(gid):
    """Internal function to strip transcript dot-notation from IDs."""
    return re.sub(r"\.\d*", "", gid)


def safe_eval_gpr(expr, conf_genes):
    """Internal function to evaluate a gene-protein rule in an
    injection-safe manner (hopefully).
    """
    if isinstance(expr, Expression):
        return safe_eval_gpr(expr.body, conf_genes)
    elif isinstance(expr, Name):
        fgid = format_gid(expr.id)
        if fgid not in conf_genes:
            return 0
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
    rule.

    Args:
        rule (str): A gene-reaction rule. For instance "A and B".
        conf_genes (dict): A str->int map denoting the mapping of gene IDs
            to expression confidence values. Allowed confidence values are -1
            (absent/do not include), 0 (unknown), 1 (low confidence),
            2 (medium confidence) and 3 (high confidence).
    """
    ast_rule, _ = parse_gpr(rule)
    return safe_eval_gpr(ast_rule, conf_genes)


def test_model():
    """Gets a small test model.

    Returns a small test model for the central carbon metabolism.

    Args:
        None:

    Returns:
        cobra model: A model of the central carbon metabolism.
    """
    from os.path import split, join
    from cobra.io import read_sbml_model
    this_dir, _ = split(__file__)
    data_path = join(this_dir, "data", "cemet.xml")
    return read_sbml_model(data_path)
