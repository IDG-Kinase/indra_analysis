import os
import json
import numpy
import boto3
import pickle
import pandas
import logging
from indra.statements import stmts_to_json
from indra.assemblers.tsv import TsvAssembler
from indra.assemblers.html import HtmlAssembler
from indra.tools import assemble_corpus as ac
from indra.belief import BeliefEngine
from indra.databases import hgnc_client
from indra.sources.indra_db_rest import get_statements
from indra_db.client.statements import get_statements_by_gene_role_type


logger = logging.getLogger('get_kinase_interactions')


def get_statements_for_kinase_db_direct(kinase):
    logger.info('Getting statements for %s' % kinase)
    hgnc_id = hgnc_client.get_current_hgnc_id(kinase)
    if hgnc_id is None:
        logger.warning('Could not get HGNC ID for %s' % kinase)
        return None
    stmts = get_statements_by_gene_role_type(agent_id=hgnc_id,
                                             agent_ns='HGNC',
                                             with_support=False)
    stmts = filter_out_medscan(stmts)
    stmts = sorted(stmts, key=lambda x: len(x.evidence), reverse=True)
    return stmts


def get_statements_for_kinase_db_api(kinase):
    logger.info('Getting statements for %s' % kinase)
    hgnc_id = hgnc_client.get_current_hgnc_id(kinase)
    if hgnc_id is None:
        logger.warning('Could not get HGNC ID for %s' % kinase)
        return None
    ip = get_statements(agents=['%s@HGNC' % hgnc_id],
                        ev_limit=10000)
    stmts = filter_out_medscan(ip.statements)
    stmts = sorted(stmts, key=lambda x: len(x.evidence), reverse=True)
    return stmts


def get_kinase_statements(kinases, db_mode='api'):
    """Get all statements from the database for a list of gene symbols."""
    all_statements = {}
    for kinase in kinases:
        if db_mode == 'direct':
            stmts = get_statements_for_kinase_db_direct(kinase)
        elif db_mode == 'api':
            stmts = get_statements_for_kinase_db_api(kinase)
        else:
            raise ValueError('Invalid DB mode: %s' % db_mode)
        if stmts is None:
            continue
        all_statements[kinase] = stmts
    return all_statements


def filter_out_medscan(stmts):
    logger.info('Starting medscan filter with %d statements' % len(stmts))
    new_stmts = []
    for stmt in stmts:
        new_evidence = []
        for ev in stmt.evidence:
            if ev.source_api == 'medscan':
                continue
            new_evidence.append(ev)
        if new_evidence:
            new_stmts.append(stmt)
    logger.info('Finished medscan filter with %d statements' % len(new_stmts))
    return new_stmts


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


def make_all_kinase_statements(kinases, prefix, db_mode, exports=None):
    before_assembly_pickle = f'{prefix}_before_assembly.pkl'
    after_assembly_pickle = f'{prefix}.pkl'
    if os.path.exists(after_assembly_pickle):
        with open(after_assembly_pickle, 'rb') as fh:
            stmts = pickle.load(fh)
    else:
        if os.path.exists(before_assembly_pickle):
            with open(before_assembly_pickle, 'rb') as fh:
                stmts = pickle.load(fh)
        else:
            # Get all statements for kinases
            stmts = get_kinase_statements(kinases, db_mode=db_mode)
            with open(before_assembly_pickle, 'wb') as fh:
                pickle.dump(stmts, fh)
        stmts = assemble_statements(stmts)
        with open(after_assembly_pickle, 'wb') as fh:
            pickle.dump(stmts, fh)

    # Export into JSON and TSV
    exports = exports if exports else set()
    if 'json' in exports:
        export_json(stmts, f'{prefix}.json')
    if 'tsv' in exports:
        export_tsv(stmts, f'{prefix}.tsv')
    if 's3' in exports:
        dump_html_to_s3(stmts)

    print_statistics(stmts)
    return stmts


def get_kinase_list(fname, col_name):
    df = pandas.read_table(fname, sep=',')
    kinases = sorted(set(df[col_name]))
    return kinases


def assemble_statements(stmts):
    """Run assembly steps on statements."""
    be = BeliefEngine()
    for kinase, kinase_stmts in stmts.items():
        stmts[kinase] = ac.filter_human_only(kinase_stmts)
        be.set_prior_probs(stmts[kinase])
    return stmts


def dump_html_to_s3(stmts):
    s3 = boto3.client('s3')
    bucket = 'dark-kinases'
    for kinase, sts in stmts.items():
        fname = f'{kinase}.html'
        sts_sorted = sorted(sts, key=lambda x: len(x.evidence), reverse=True)
        ha = HtmlAssembler(sts_sorted)
        html_str = ha.make_model()
        url = 'https://s3.amazonaws.com/%s/%s' % (bucket, fname)
        print('Dumping to %s' % url)
        s3.put_object(Key=fname, Body=html_str.encode('utf-8'), Bucket=bucket,
                      ContentType='text/html')


if __name__ == '__main__':
    # Get all dark kinase Statements
    #fname = 'Table_005_IDG_dark_kinome.csv'
    #prefix = 'dark_kinase_statements_v4'
    #col_name = 'gene_symbol'
    #make_all_kinase_statements(fname, prefix, col_name, ev_limit=10000)

    # Get all kinase Statements
    # There are 722 kinases but only 709 unique here
    fname = 'allsources_HMS_it3_cleaned_manual.csv'
    prefix = 'all_kinase_statements_v7'
    col_name = 'HGNC_name'
    db_mode = 'api'
    kinases = get_kinase_list(fname, col_name)
    stmts = make_all_kinase_statements(kinases, prefix, db_mode=db_mode,
                                       exports={'s3'})