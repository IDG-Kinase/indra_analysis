import json
import tqdm
import pandas
import pickle
from collections import defaultdict
from indra.databases import hgnc_client
from indra.ontology.bio import bio_ontology
from indra.sources import indra_db_rest
from indra.tools import assemble_corpus as ac
from get_kinase_interactions import filter_out_medscan


kinase_pkl = 'all_kinase_statements_v6.pkl'


def get_fplx_stmts(fplx_id):
    ip = indra_db_rest.get_statements(agents=['%s@FPLX' % fplx_id],
                                      ev_limit=10000)
    stmts = filter_out_medscan(ip.statements)
    stmts = ac.filter_human_only(stmts)
    return stmts


if __name__ == '__main__':
    with open(kinase_pkl, 'rb') as fh:
        kinase_stmts = pickle.load(fh)
    fplx_by_kinase = defaultdict(set)
    kinase_by_fplx = defaultdict(set)
    kinase_counts = {}
    for kinase, stmts in kinase_stmts.items():
        hgnc_id = hgnc_client.get_hgnc_id(kinase)
        parents = bio_ontology.get_parents('HGNC', hgnc_id)
        fplx_by_kinase[kinase] |= {fplx_id for _, fplx_id in parents}
        for _, fplx_id in parents:
            kinase_by_fplx[fplx_id].add(kinase)
        kinase_counts[kinase] = len(stmts)
    kinase_by_fplx = dict(kinase_by_fplx)
    fplx_by_kinase = dict(fplx_by_kinase)

    fplx_stmts = {}
    fplx_counts = {}
    for fplx_id in tqdm.tqdm(kinase_by_fplx):
        fplx_stmts[fplx_id] = get_fplx_stmts(fplx_id)
        fplx_counts[fplx_id] = len(fplx_stmts[fplx_id])
    with open('fplx_stmts.pkl', 'wb') as fh:
        pickle.dump(fplx_stmts, fh)
    with open('fplx_counts.json', 'w') as fh:
        json.dump(fplx_counts, fh, indent=1)

    fplx_stats = []
    for fplx_id, kinases in kinase_by_fplx.items():
        child_counts = [kinase_counts[kinase] for kinase in kinases]
        parent_count = fplx_counts[fplx_id]
        fplx_stats.append(
            {'fplx_id': fplx_id,
             'children': ', '.join(sorted(kinases)),
             'family_count': parent_count,
             'total_child': sum(child_counts),
             'mean_child': sum(child_counts)/len(child_counts),
             'total_ratio': parent_count / sum(child_counts),
             'mean_ratio': parent_count / (sum(child_counts) /
                len(child_counts))}
        )
    df = pandas.DataFrame.from_records(fplx_stats)
    df.sort_values(by='total_ratio', inplace=True, ascending=False)
    df.to_csv('family_stats.tsv', sep='\t', index=None)
