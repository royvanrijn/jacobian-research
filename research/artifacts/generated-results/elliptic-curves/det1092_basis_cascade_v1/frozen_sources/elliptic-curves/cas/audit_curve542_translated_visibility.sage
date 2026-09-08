#!/usr/bin/env sage-python
"""Oracle-only translations of a certified missing direction in43 frozen charts."""
import argparse
from collections import Counter
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import QQ,ZZ,RealField,matrix,vector,pari,EllipticCurve
from fpylll import GSO,IntegerMatrix,Enumeration
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from search_observability import point_visibility
from research_runtime.store import checkpoint
geometry=SourceFileLoader('curve542_translation_geometry',str(CAS/'prospective_half_lattice.sage')).load_module()
AUDIT=ROOT/'artifacts/generated-results/elliptic-curves/curve542_initial_visibility_v1.json'
INPUT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1/a1-fibration-04/candidate-00/result.json'

def run(output):
    if output.exists():raise FileExistsError('preserve translated audit')
    audit=cert.read(AUDIT);source=cert.read(INPUT)
    if cert.hashed(INPUT)!=audit['input_sha256']:raise ArithmeticError('retained initial source changed')
    model=tuple(map(cert.F,audit['curve']));raw=audit['rank26_witness']['points'];points=[tuple(map(cert.F,p)) for p in raw]
    proof=audit['rank26_witness']['rank_certificate'];cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    E=EllipticCurve(QQ,[QQ(str(q)) for q in model]);basis=[E([QQ(str(x)),QQ(str(y))]) for x,y in points[:25]];P=E([QQ(str(q)) for q in points[25]])
    gram,asymmetry=geometry.canonical_height_gram(model,points);RF=RealField(256);full=matrix(RF,[[RF(str(x)) for x in row] for row in gram]);G=full[:25,:25]
    projection=G.inverse()*full[:25,25:26];scaled=matrix(ZZ,[[ZZ((x*2**20).round()) for x in row] for row in G.rows()])
    U=matrix(ZZ,pari(scaled).qflllgram()).transpose()
    if abs(U.det())!=1:raise ArithmeticError('LLL transport not unimodular')
    reduced=U*scaled*U.transpose();gso=GSO.Mat(IntegerMatrix.from_matrix([list(map(int,r)) for r in reduced.rows()]),gram=True,float_type='dd',update=True)
    mu=matrix(RF,[[RF(gso.get_mu(i,j)) if i>j else int(i==j) for j in range(25)] for i in range(25)]);inverse=U.inverse()
    result={'schema':'elliptic-curves.curve542-translated-visibility.v1','retrospective_only':True,
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),AUDIT,INPUT,CAS/'prospective_half_lattice.sage',CAS/'search_observability.py')},
        'curve':source['curve'],'basis':[list(map(str,p)) for p in points[:25]],'oracle_point':list(map(str,points[25])),
        'full_gram':[[str(x) for x in row] for row in gram],'maximum_asymmetry':str(asymmetry),'lll_matrix':[list(map(int,r)) for r in U.rows()],
        'proposals':[],'status':'RUNNING','claim_boundary':'One floating-CVP translation per sign and retained chart; exact group words and chart lifts. No CVP optimality, complete coset visibility, new curve, or rank upper bound. Oracle outputs never enter prospective selection or search.'}
    checkpoint(output,result);counts=Counter()
    for i,row in enumerate(source['charts']):
        current=row['search']['input']['subgroup']
        if [(p['x'],p['y']) for p in current]!=[tuple(p) for p in raw[:len(current)]]:raise ArithmeticError('chart subgroup not a prefix of the discovered25 basis')
        coefficients=row['search']['input']['centre']['coefficients'];centre=vector(QQ,coefficients+[0]*(25-len(coefficients)))
        for sign in (1,-1):
            target=(centre/2-sign*projection.column(0))*inverse;gs_target=tuple(map(float,target*mu));trial=vector(ZZ,[round(float(x)) for x in target]);delta=trial-target
            radius=float(delta*reduced*delta)+1;solutions=Enumeration(gso).enumerate(0,25,radius,0,target=gs_target)
            if not solutions:raise ArithmeticError('translation CVP returned no proposal')
            _,coordinates=solutions[0];small=vector(ZZ,[round(x) for x in coordinates])
            if any(abs(float(x)-int(y))>1e-7 for x,y in zip(coordinates,small)):raise ArithmeticError('nonintegral CVP proposal')
            word=small*U;translated=sign*P+sum((int(k)*p for k,p in zip(word,basis) if k),E(0))
            if translated.is_zero():raise ArithmeticError('certified escaping direction became zero')
            point=[str(translated[0]),str(translated[1])];v=point_visibility(row['search'],point);counts[v['status']]+=1
            result['proposals'].append({'chart_index':i,'sign':sign,'word':list(map(int,word)),'point':point,'visibility':v})
        checkpoint(output,result)
        if (i+1)%10==0:print('TRANSLATED542',i+1,dict(counts),flush=True)
    best=min(result['proposals'],key=lambda r:r['visibility'].get('minimum_affine_height') or 10**1000)
    result.update(status='COMPLETE_DECLARED_ORACLE_AUDIT',status_counts=dict(counts),best=best);checkpoint(output,result)
    print('BEST TRANSLATED542',best,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
