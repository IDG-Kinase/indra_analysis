import pickle
from gilda.api import *
import requests
from collections import defaultdict
from indra.preassembler.grounding_mapper.adeft import _get_text_for_grounding
from indra.preassembler.grounding_mapper.gilda import *
from gilda.grounder import load_gilda_models

gilda_models = load_gilda_models()


if __name__ == '__main__':
    with open('all_kinase_statements_v6.pkl', 'rb') as fh:
        stmts_dict = pickle.load(fh)

    gilda_strs = {'AIM1',
                  'FLG',
                  'Ndr1',
                  'OSR1',
                  'Osr-1',
                  'Osr1',
                  'PDK1',
                  'SK2',
                  'Stk1',
                  'Warts',
                  'osr1',
                  'warts'}

    gilda_str_groundings = {
        s: gilda_models[s].pos_labels for s in gilda_strs
    }


    def get_stmt_hgnc_ids(stmt):
        agents = [a for a in stmt.agent_list() if a is not None]
        hgnc_ids = {'HGNC:%s' % a.db_refs.get('HGNC')
                    for a in agents if 'HGNC' in a.db_refs}
        return hgnc_ids

    # Stage 1: identify groundings that are affected by the given entity strings
    txt_to_sentence = defaultdict(list)
    for kinase, stmts in stmts_dict.items():
        for stmt in stmts:
            hgnc_ids = get_stmt_hgnc_ids(stmt)
            for ev in stmt.evidence:
                for gs in gilda_strs:
                    if ev.text and gs in ev.text and \
                            set(hgnc_ids) & set(gilda_str_groundings[gs]):
                        txt_to_sentence[gs].append((stmt, ev))

    # Stage 2: Identify groundings that are incorrect per a model and create a
    # list of statement/evidence pairs to remove
    to_remove = []
    for gs, cases in txt_to_sentence.items():
        for stmt, ev in cases:
            gr_text = _get_text_for_grounding(stmt, gs)
            terms = ground(gs, ev.text)
            top_term = terms[0].term
            hgnc_ids = get_stmt_hgnc_ids(stmt)
            if '%s:%s' % (top_term.db, top_term.id) not in hgnc_ids:
                to_remove.append((stmt, ev))
    # Save the removables for reference
    with open('to_remove.pkl', 'wb') as fh:
        pickle.dump(to_remove, fh)

    with open('to_remove.pkl', 'rb') as fh:
        to_remove = pickle.load(fh)

    # Create a lookup index to find statement/evidence pairs to remove
    to_remove_index = {(s.matches_key(), e.matches_key())
                       for s, e in to_remove}

    # Stage 3: Create a new statement and evidence structure with the removables
    # removed
    new_stmts_dict = {}
    for kinase, stmts in stmts_dict.items():
        new_stmts = []
        for stmt in stmts:
            new_ev = []
            for ev in stmt.evidence:
                if (stmt.matches_key(), ev.matches_key()) in to_remove_index:
                    print('Removing %s' % str((stmt, ev)))
                    continue
                new_ev.append(ev)
            if new_ev:
                new_stmts.append(stmt)
            else:
                print('No evidence left for %s' % str(stmt))
        new_stmts_dict[kinase] = new_stmts

    # Save the results
    with open('all_kinase_statements_v6_patched.pkl', 'wb') as fh:
        pickle.dump(new_stmts_dict, fh)
