import os
import json
import pickle
from collections import Counter
import matplotlib.pyplot as plt
from indra.statements import *


if os.path.exists('sources.json'):
    with open('sources.json', 'r') as fh:
        sources = json.load(fh)
else:
    with open('all_kinase_statements_v6_patched.pkl', 'rb') as fh:
        stmts_dict = pickle.load(fh)


    sources = []
    for kinase, stmts in stmts_dict.items():
        for stmt in stmts:
            if isinstance(stmt, Modification) and stmt.enz is None:
                continue
            for ev in stmt.evidence:
                source = ev.source_api
                if ev.source_api == 'biopax':
                    source = 'biopax:%s' % ev.annotations['source_sub_id']
                sources.append(source)

    with open('sources.json', 'w') as fh:
        json.dump(sources, fh, indent=1)

cnt = Counter(sources)
source_order, source_counts = zip(*cnt.most_common())

readers = ['reach', 'sparser', 'medscan', 'rlimsp', 'trips', 'isi']

labels = {
    'reach': 'REACH',
    'sparser': 'Sparser',
    'medscan': 'Medscan',
    'rlimsp': 'RLIMS-P',
    'trips': 'TRIPS',
    'bel': 'BEL',
    'hprd': 'HPRD',
    'biopax:phosphositeplus': 'PhosphoSitePlus',
    'biopax:kegg': 'KEGG',
    'biopax:recon x': 'ReconX',
    'biopax:panther': 'Panther',
    'biopax:ctd': 'CTD',
    'biopax:humancyc': 'HumanCyc',
    'biopax:inoh': 'INOH',
    'signor': 'Signor',
    'biopax:pid': 'NCI-PID',
    'biopax:mirtarbase': 'MiRTarBase',
    'biopax:reactome': 'Reactome',
    'biopax:netpath': 'NetPath',
    'isi': 'ISI',
    'tas': 'TAS',
    'trrust': 'TRRUST',
    'biogrid': 'BioGRID'
    }

labeled_sources = [labels[s] for s in source_order]
source_colors = ['red' if s in readers else 'blue' for s in source_order]
plt.figure()
plt.bar(labeled_sources, source_counts, color=source_colors)
plt.xticks(range(len(labeled_sources)), rotation='vertical')
plt.yscale('log')
plt.ylabel('Number of INDRA Statement evidences')
plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.24)
plt.show()
