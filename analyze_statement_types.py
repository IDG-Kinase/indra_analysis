import pickle
from collections import Counter
from indra.statements import *


with open('all_kinase_statements_v6_patched.pkl', 'rb') as fh:
    stmts_dict = pickle.load(fh)


def get_position(stmt, name):
    if isinstance(stmt, Complex):
        return 0
    elif len(stmt.agent_list()) == 1:
        return 0
    elif isinstance(stmt, Conversion):
        return 1
    else:
        agent_names = [a.name if a is not None else None
                       for a in stmt.agent_list()]
        assert len(agent_names) == 2
        if name == agent_names[0]:
            return 1
        else:
            return -1


def get_downstream_phospho_targets(stmts, name):
    subs = set()
    for stmt in stmts:
        if not isinstance(stmt, Phosphorylation):
            continue
        if not stmt.enz or stmt.enz.name != name:
            continue
        subs.add(name)
    return sorted(list(subs))


stmt_type_counts = {}
upstream_vs_downstream = {}
downstream_phospho_targets = {}
for kinase, stmts in stmts_dict.items():
    stmt_type_counts[kinase] = Counter([s.__class__.__name__ for s in stmts])
    upstream_vs_downstream[kinase] = Counter([get_position(s, kinase)
                                              for s in stmts])
    downstream_phospho_targets[kinase] = \
        len(get_downstream_phospho_targets(stmts, kinase))



