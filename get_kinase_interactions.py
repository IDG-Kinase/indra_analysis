import os
import json
import numpy
import pickle
import pandas
from indra.statements import stmts_to_json
from indra.assemblers.tsv import TsvAssembler
from indra.tools import assemble_corpus as ac
from indra_db.client import get_primary_db, get_statements_by_gene_role_type


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


def print_statistics(statements):
    counts = sorted([(k, len(s)) for k, s in statements.items()],
                    key=lambda x: x[1], reverse=True)
    raw_counts = [c[1] for c in counts]
    missing = [c[0] for c in counts if c[1] == 0]
    print(f'No statements for kinases: {", ".join(missing)}')
    print(f'{counts[0][1]} statements for the top kinase {counts[0][0]}')
    print(f'{numpy.mean(raw_counts)} statements on average per kinase')


def make_all_kinase_statements(fname, prefix, col_name):
    # If we have a pickle just reuse that
    if not os.path.exists(f'{prefix}.pkl'):
        df = pandas.read_table(fname, sep=',')
        kinases = list(df[col_name])
        # Get all statements for kinases
        stmts = get_kinase_statements(kinases)
        with open(f'{prefix}_before_assembly.pkl', 'wb') as fh:
            pickle.dump(stmts, fh)
        stmts = assemble_statements(stmts)
        with open(f'{prefix}.pkl', 'wb') as fh:
            pickle.dump(stmts, fh)
    else:
        with open(f'{prefix}.pkl', 'rb') as fh:
            stmts = pickle.load(fh)
    # Export into JSON and TSV
    export_json(stmts, f'{prefix}.json')
    export_tsv(stmts, f'{prefix}.tsv')

    print_statistics(stmts)
    return stmts


def assemble_statements(stmts):
    for kinase, kinase_stmts in stmts.items():
        stmts[kinase] = ac.filter_human_only(kinase_stmts)
    return stmts


if __name__ == '__main__':
    # Get all dark kinase Statements
    fname = 'Table_005_IDG_dark_kinome.csv'
    prefix = 'dark_kinase_statements_v2'
    col_name = 'gene_symbol'
    make_all_kinase_statements(fname, prefix, col_name)

    # Get all kinase Statements
    #fname = 'Table_001_all_kinases.csv'
    #prefix = 'all_kinase_statements'
    #col_name = 'gene_symbol'
    #make_all_kinase_statements(fname, prefix, col_name)
