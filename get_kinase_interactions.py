import os
import json
import tqdm
import types
import numpy
import boto3
import pickle
import pandas
import logging
from collections import defaultdict
from indra.statements import stmts_to_json, Inhibition, Activation, \
    IncreaseAmount, DecreaseAmount, Phosphorylation, Dephosphorylation
import indra.statements.agent
from indra.assemblers.tsv import TsvAssembler
from indra.assemblers.html import HtmlAssembler
from indra.tools import assemble_corpus as ac
from indra.databases import hgnc_client
from indra.sources.indra_db_rest import get_statements
from indra_db.client.principal.curation import get_curations


logger = logging.getLogger('get_kinase_interactions')


iupac_to_chebi_dict = pickle.load(open('data/iupac_to_chebi_dict.p', 'rb'))
chebi_to_selleck_dict = pickle.load(open('data/chebi_to_selleck_dict', 'rb'))


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


def filter_out_medscan(stmts):
    logger.info('Starting medscan filter with %d statements' % len(stmts))
    new_stmts = []
    for stmt in stmts:
        new_evidence = []
        for ev in stmt.evidence:
            if ev.source_api == 'medscan':
                continue
            new_evidence.append(ev)
        stmt.evidence = new_evidence
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


def rename_chemical(agent):
    from indra.ontology.bio import bio_ontology
    if agent.db_refs['CHEBI'] in chebi_to_selleck_dict:
        selleckname = chebi_to_selleck_dict[agent.db_refs['CHEBI']]
        if selleckname:
            agent.name = selleckname
            return
    for db_ns, db_id in sorted(agent.db_refs.items()):
        if db_ns in {'TEXT', 'TEXT_NORM', 'CHEBI'}:
            continue
        name = bio_ontology.get_name(db_ns, db_id)
        if name:
            agent.name = name
            return


def remove_contradictions(stmts):
    logger.info('Filtering contradictions on %d statements' % len(stmts))
    stmts_by_agents = defaultdict(list)
    for stmt in stmts:
        if not len(stmt.real_agent_list()) == 2:
            continue
        stmts_by_agents[tuple(agent.name
                              for agent in stmt.real_agent_list())].append(stmt)
    to_remove = set()
    contradiction_pairs = [
        (Activation, Inhibition),
        (IncreaseAmount, DecreaseAmount),
        (Phosphorylation, Dephosphorylation),
    ]
    for agent_names, stmts_for_agents in stmts_by_agents.items():
        stmt_types = defaultdict(int)
        for stmt in stmts_for_agents:
            stmt_types[type(stmt)] += len(stmt.evidence)
        for stmt_type1, stmt_type2 in contradiction_pairs:
            if {stmt_type1, stmt_type2} <= set(stmt_types):
                if stmt_types[stmt_type1] < stmt_types[stmt_type2]:
                    to_remove |= {stmt.get_hash() for stmt in stmts_for_agents
                                  if isinstance(stmt, stmt_type1)}
                elif stmt_types[stmt_type2] < stmt_types[stmt_type1]:
                    to_remove |= {stmt.get_hash() for stmt in stmts_for_agents
                                  if isinstance(stmt, stmt_type2)}
    stmts = [s for s in stmts if s.get_hash() not in to_remove]
    logger.info('Finishing with %d statements' % len(stmts))
    return stmts


def assemble_statements(stmts, curs):
    """Run assembly steps on statements."""
    # Remove unary statements and ones with many agents
    stmts = [stmt for stmt in stmts
             if (1 < len(stmt.real_agent_list()) < 4)]
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_by_curation(stmts, curations=curs)
    stmts = unify_lspci(stmts)
    stmts = remove_contradictions(stmts)
    # Rename chemicals
    for stmt in stmts:
        for agent in stmt.real_agent_list():
            if agent.db_refs.get('CHEBI') and len(agent.name) > 25:
                rename_chemical(agent)
    # Remove long names
    stmts = [stmt for stmt in stmts if
             all(len(a.name) < 20 for a in stmt.real_agent_list())]
    # Remove microRNAs
    stmts = [stmt for stmt in stmts
             if not any('miR' in a.name for a in stmt.real_agent_list())]
    return stmts


def unify_lspci(stmts):
    logger.info('Unifying by LSPCI with %d statements' % len(stmts))
    orig_ns_order = indra.statements.agent.default_ns_order[:]
    indra.statements.agent.default_ns_order = ['LSPCI'] + \
        indra.statements.agent.default_ns_order
    from indra.ontology.bio import bio_ontology
    agents_by_lspci = defaultdict(list)
    from indra.statements.agent import default_ns_order
    ns_order = default_ns_order + ['CHEMBL', 'DRUGBANK', 'HMS-LINCS', 'CAS']
    for stmt in stmts:
        for agent in stmt.real_agent_list():
            if 'LSPCI' in agent.db_refs:
                agents_by_lspci[agent.db_refs['LSPCI']].append(agent)
            else:
                agent_gr = agent.get_grounding(ns_order=ns_order)
                if agent_gr[0] is None:
                    continue
                else:
                    parents = bio_ontology.get_parents(*agent_gr)
                    lspci_parents = [p[1] for p in parents if p[0] == 'LSPCI']
                    if len(lspci_parents) != 1:
                        continue
                    lspci_parent = lspci_parents[0]
                    agents_by_lspci[lspci_parent].append(agent)

    for lspci, agents in agents_by_lspci.items():
        lspci_name = bio_ontology.get_name('LSPCI', lspci)
        standard_name = lspci_name if lspci_name else agents[0].name
        for agent in agents:
            agent.db_refs['LSPCI'] = lspci
            agent.name = standard_name

    from indra.preassembler import Preassembler
    pa = Preassembler(bio_ontology, stmts)
    unique_stmts = pa.combine_duplicates()
    indra.statements.agent.default_ns_order = orig_ns_order
    logger.info('Finished unification with %d statements' % len(unique_stmts))
    return unique_stmts


def dump_html_to_s3(kinase, stmts):
    s3 = boto3.client('s3')
    bucket = 'dark-kinases'
    fname = f'{kinase}.html'
    sts_sorted = sorted(stmts, key=lambda x: len(x.evidence), reverse=True)
    ha = HtmlAssembler(sts_sorted, db_rest_url='https://db.indra.bio')
    html_str = ha.make_model(no_redundancy=True)
    url = 'https://s3.amazonaws.com/%s/%s' % (bucket, fname)
    print('Dumping to %s' % url)
    s3.put_object(Key=fname, Body=html_str.encode('utf-8'), Bucket=bucket,
                  ContentType='text/html')


def load_raw_stmts(kinase):
    fname = os.path.join('data', f'{kinase}.pkl')
    if os.path.exists(fname):
        logger.info('Loading cached %s' % fname)
        with open(fname, 'rb') as fh:
            stmts = pickle.load(fh)
    else:
        stmts = get_statements_for_kinase_db_api(kinase)
        with open(fname, 'wb') as fh:
            pickle.dump(stmts, fh)
    return stmts


if __name__ == '__main__':
    # Get all dark kinase Statements
    #fname = 'Table_005_IDG_dark_kinome.csv'
    #prefix = 'dark_kinase_statements_v4'
    #col_name = 'gene_symbol'
    #make_all_kinase_statements(fname, prefix, col_name, ev_limit=10000)

    # Get all kinase Statements
    # There are 722 kinases but only 709 unique here
    fname = 'data/allsources_HMS_it3_cleaned_manual.csv'
    kinases = get_kinase_list(fname, 'HGNC_name')
    curs = get_curations()
    all_stmts = {}
    for kinase in tqdm.tqdm(kinases):
        stmts = load_raw_stmts(kinase)
        stmts = assemble_statements(stmts, curs)
        dump_html_to_s3(kinase, stmts)
        all_stmts[kinase] = stmts
    with open('data/all_stmts.pkl', 'wb') as fh:
        pickle.dump(all_stmts, fh)