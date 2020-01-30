import sys
import pickle
from indra.databases import ndex_client
from indra.assemblers.cx import CxAssembler
from indra.assemblers.cx.hub_layout import add_semantic_hub_layout


if __name__ == '__main__':
    kinase = sys.argv[1]
    with open('dark_kinase_statements_assembled.pkl', 'rb') as fh:
        stmts = pickle.load(fh)
    cxa = CxAssembler(stmts[kinase])
    cxa.make_model()
    add_semantic_hub_layout(cxa.cx, kinase)
    model_id = cxa.upload_model()
    ndex_client.set_style(model_id)
