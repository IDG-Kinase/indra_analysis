import time
import pickle
from indra.databases import hgnc_client
from indra.literature.pubmed_client import get_ids_for_gene, get_ids

prefix = 'dark_kinase_statements_v4'


def load_stmts(prefix):
    with open(f'{prefix}.pkl', 'rb') as fh:
        return pickle.load(fh)


def get_pmids_entrez(kinase):
    pmids = get_ids_for_gene(kinase)
    time.sleep(1)
    return pmids


def get_pmids_text(kinase):
    pmids = get_ids(kinase)
    time.sleep(1)
    return pmids


def get_raw_strings(all_stmts):
    raw_strings = {}
    for kinase, stmts in all_stmts.items():
        raw_strings[kinase] = []
        hgnc_id = hgnc_client.get_hgnc_id(kinase)
        if not hgnc_id:
            continue
        for stmt in stmts:
            for agent in stmt.agent_list():
                if agent is not None and agent.db_refs.get('HGNC') == hgnc_id:
                    text = agent.db_refs.get('TEXT')
                    if text:
                        raw_strings[kinase].append(text)
    return raw_strings


if __name__ == '__main__':
    all_stmts = load_stmts(prefix)
    '''
    stmt_counts = []
    pmid_counts = []
    for kinase, stmts in all_stmts.items():
        pmids = get_pmids_text(kinase)
        pmid_counts.append(len(pmids))
        stmt_counts.append(len(stmts))
        print('%s: %d PMIDs, %d stmts' % (kinase, len(pmids), len(stmts)))
    '''
    raw_strings = get_raw_strings(all_stmts)
