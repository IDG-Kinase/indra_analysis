import os
import json
import pandas
from indra.statements import stmts_to_json
from indra.db.client import get_primary_db, get_statements_by_gene_role_type


def get_kianse_statements(kinases):
    all_statements = {}
    db = get_primary_db()
    for kinase in kinases:
        statements = get_statements_by_gene_role_type(kinase, db=db,
                                                      with_support=False)
        all_statements[kinase] = statements
    return all_statements


def export_tsv(statements):
    for kinase, stmts in sorted(stmts.items(), key=lambda x: x[0]):
        agents = stmt.agent_list()


def export_json(statements):
    stmts_json = {}
    for kinase, stmts in sorted(all_statements.items(), key=lambda x: x[0]):
        sj = stmts_to_json(stmts)
        stmts_json[kinase] = sj
    return stmts_json


if __name__ == '__main__':
    fname = 'Table_005_IDG_dark_kinome.csv'

    if not os.path.exists('dark_kinase_statements.pkl'):
        df = pandas.read_table(fname, sep=',')
        kinases = list(df['gene_symbol'])
        stmts = get_kianse_statements(kinases)
        with open('dark_kinase_statements.pkl', 'wb') as fh:
            pickle.dump(stmts, fh)

    stmts_json = export_json(stmts)
    with open('dark_kinase_statements.json', 'w') as fh:
        json.dump(stmts_json, fh, indent=1)
