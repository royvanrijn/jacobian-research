"""Opt-in fitted model policy from the sealed finite-cancellation experiment.

The factor-free V3 default is unchanged. This policy takes only the current
equation, known subgroup and proposed centre, and constructs at most3 models.
The fixed12-case control has a measured advantage but misses its strict10%
aggregate CPU gate; this is an experimental mapping option.
"""
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path

from finite_cancellation_features import prepare
from half_lattice_pointed_sieve import linear_combination
from finite_cancellation_corpus import ROOT, canonical, digest

CAS=Path(__file__).resolve().parent
_base=SourceFileLoader('finite_policy_base_mapper',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
pari=_base.pari
POLICY=ROOT/'elliptic-curves/data/finite_cancellation_policy_v1.json'

def select(features,fit):
    base=features[0];scores=[]
    for f in features:
        z=[(f[k]-base[k]-mu)/sd for k,mu,sd in zip(fit['features'],fit['means'],fit['scales'])]
        scores.append(fit['coefficients'][0]+sum(a*b for a,b in zip(fit['coefficients'][1:],z)))
    return min(range(len(scores)),key=lambda i:scores[i]),scores

def mapping(model,points,centre,progress=lambda stage:None):
    policy=json.loads(POLICY.read_text())
    for name,h in policy['feature_sources'].items():
        if digest((CAS/name).read_bytes())!=h:raise ArithmeticError('frozen predictor feature source changed')
    progress('known_anchor')
    Q=linear_combination(model,points,centre['representative'])
    if Q is None:raise ArithmeticError('finite nonzero centre required')
    progress('factor_free_models_and_residue_features')
    prepared=prepare(list(map(str,model)),[list(map(str,Q))],0,_base)
    features=[m['features'] for m in prepared['models']]
    i,scores=select(features,policy['fit']);selected=prepared['models'][i]
    out=dict(selected['mapping']);out['centre']=centre
    # Retain arithmetic inputs directly; never hash timings into the sole
    # reconstruction witness (the earlier CPU diagnostic is preserved).
    arithmetic={k:v for k,v in prepared.items() if k not in ['base_map_cpu_seconds','preparation_cpu_seconds']}
    out['finite_cancellation_receipt']={'policy_sha256':digest(POLICY.read_bytes()),
        'arithmetic_preparation_sha256':digest(canonical(arithmetic)),'selected':selected['name'],
        'candidate_features':features,'scores':scores,'models_built':len(features),
        'preparation_cpu_seconds':prepared['preparation_cpu_seconds'],
        'status':'EXPERIMENTAL_STRICT_CPU_GATE_MISSED'}
    progress('exact_selected_map')
    return out
