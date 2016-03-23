#  util.py
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

import re
from cobra.manipulation import revert_to_reversible

class ci(int): 
    def __and__(self, other):
        return min(self, other)
    
    def __or__(self, other):
        return max(self, other)

def format_gid(gid):
    return re.sub(r"\.\d*", "", gid)

def reaction_confidence(rule, conf_genes):
    """Calculates the confidence for the reaction based on a gene-reaction
        rule"""
    
    if rule == '': return 0
    
    rule = rule.replace("and", "&").replace("or", "|")
    matches = re.finditer(r"\w*\d+\.*\d*", rule)
    shift = 0
    for m in matches:
        fgid = format_gid(m.group())
        rep = "ci(" + str(conf_genes[fgid]) + ")" if fgid in conf_genes \
            else "ci(0)"
        # replace won't work here because IDs may include other IDs 
        # for instance 1234.1 includes 34.1 and 234.1 etc.
        s, e = m.start() + shift, m.end() + shift
        rule = rule[:s] + rep + rule[e:]
        shift += len(rep) - len(m.group())
    
    return eval(rule)
    

def safe_revert_reversible(model):
    maybe_rev = (r for r in model.reactions if "reflection" in r.notes)
    
    for r in maybe_rev:
        try:
            model.reactions.get_by_id(r.notes["reflection"])
        except KeyError:
            r.notes.pop("reflection")
    
    revert_to_reversible(model)
