import os
import json
import pickle
import pandas
from indra.statements import stmts_to_json
from indra.assemblers import TsvAssembler
from indra.db.client import get_primary_db, get_statements_by_gene_role_type


def get_kinase_statements(kinases):
    """Get all statements from the database for a list of gene symbols."""
    all_statements = {}
    db = get_primary_db()
    for kinase in kinases:
        statements = get_statements_by_gene_role_type(kinase, db=db,
                                                      with_support=False)
        all_statements[kinase] = statements
    return all_statements


def export_tsv(statements, fname):
    """Export statements into TSV."""
    def get_stmt_list(statements):
        stmt_list = []
        for kinase, stmts in sorted(statements.items(), key=lambda x: x[0]):
            stmt_list += stmts
        return stmt_list
    ta = TsvAssembler(get_stmt_list(statements))
    ta.make_model(fname)


def export_json(statements, fname):
    """Export statements into JSON."""
    stmts_json = {}
    for kinase, stmts in sorted(statements.items(), key=lambda x: x[0]):
        sj = stmts_to_json(stmts)
        stmts_json[kinase] = sj
    with open(fname, 'w') as fh:
        json.dump(stmts_json, fh, indent=1)


if __name__ == '__main__':
    prefix = 'dark_kinase_statements'
    fname = 'Table_005_IDG_dark_kinome.csv'

    # If we have a pickle just reuse that
    if not os.path.exists(f'{prefix}.pkl'):
        df = pandas.read_table(fname, sep=',')
        kinases = list(df['gene_symbol'])
        # Get all statements for kinases
        stmts = get_kinase_statements(kinases)
        with open(f'{prefix}.pkl' % prefix, 'wb') as fh:
            pickle.dump(stmts, fh)
    else:
        with open(f'{prefix}.pkl', 'rb') as fh:
            stmts = pickle.load(fh)

    # Export into JSON and TSV
    export_json(stmts, f'{prefix}.json')
    export_tsv(stmts, f'{prefix}.tsv')
