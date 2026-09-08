#!/usr/bin/env sage-python
"""98 oracle-only translations on the completed native rank29 control."""
import sys
from pathlib import Path
from collections import Counter
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,RealField,matrix,vector,pari,EllipticCurve
from fpylll import GSO,IntegerMatrix,Enumeration
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from search_observability import point_visibility
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native11952-rank28-coset-visibility-v1';OUT=ART/'native11952_rank28_coset_visibility_v1.json'
AUDIT=D/'input.json';INPUT=ROOT/'artifacts/local/elliptic-curves/native11952-pari49-control-v1/candidate-00/result.json'
geometry=SourceFileLoader('native29_translation_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
def sources():
    paths=[Path(__file__).resolve(),AUDIT,INPUT,CAS/'prospective_half_lattice_v2.sage',CAS/'search_observability.py',CAS/'memory_rank_certificate.py',CAS/'replay_native11952_rank28_coset_visibility.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}
def main():
    protocol=cert.read(D/'protocol.json')
    if protocol['sources']!=sources() or protocol['oracle_published_indices']!=[27]:raise ArithmeticError('frozen retrospective protocol differs')
    if OUT.exists():raise FileExistsError('preserve oracle audit')
    audit,source=cert.read(AUDIT),cert.read(INPUT);raw=audit['basis']+audit['oracle_points'];proof=audit['rank_certificate'];witness={'points':raw,'rank_certificate':proof};model=tuple(map(cert.F,audit['curve']));points=[tuple(map(cert.F,p)) for p in raw]
    for name,h in audit['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('rank28 coset input changed')
    if source['curve']!=audit['curve'] or len(audit['basis'])!=28 or len(audit['oracle_points'])!=1 or audit['oracle_published_indices']!=[27]:raise ArithmeticError('blind28 plus missing1 gate differs')
    checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if raw[:17]!=source['generic_points'] or len(source['charts'])!=49:raise ArithmeticError('generic17 centre transport differs')
    E=EllipticCurve(QQ,[QQ(str(q)) for q in model]);basis=[E([QQ(str(x)),QQ(str(y))]) for x,y in points[:28]];oracles=[E([QQ(str(x)),QQ(str(y))]) for x,y in points[28:]]
    gram,asymmetry=geometry.canonical_height_gram(model,points);RF=RealField(384);full=matrix(RF,[[RF(str(x)) for x in row] for row in gram]);G=full[:28,:28];projections=G.inverse()*full[:28,28:29];scaled=matrix(ZZ,[[ZZ((x*2**20).round()) for x in row] for row in G.rows()]);U=matrix(ZZ,pari(scaled).qflllgram()).transpose()
    if abs(U.det())!=1:raise ArithmeticError('nonunimodular LLL')
    reduced=U*scaled*U.transpose();gso=GSO.Mat(IntegerMatrix.from_matrix([list(map(int,r)) for r in reduced.rows()]),gram=True,float_type='dd',update=True);mu=matrix(RF,[[RF(gso.get_mu(i,j)) if i>j else int(i==j) for j in range(28)] for i in range(28)]);inverse=U.inverse()
    result={'schema':'elliptic-curves.native11952-rank28-coset-visibility.v1','status':'RUNNING','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'curve':source['curve'],'basis':raw[:28],'oracle_points':raw[28:],'rank29_witness':witness,'full_gram':[[str(x) for x in row] for row in gram],'maximum_asymmetry':str(asymmetry),'lll_matrix':[list(map(int,r)) for r in U.rows()],'proposals':[],'claim_boundary':'Separate retrospective oracle audit against the original completed H100000 chart records, with the independently recovered H125000 rank28 basis. One floating CVP proposal per missing published direction, sign and frozen chart, translated only by the blindly recovered28 subgroup. Exact group words and chart identities are replayed separately. No CVP optimality, coset completeness, new point search, new curve or rank upper bound; no oracle output enters prospective centre selection.'};checkpoint(OUT,result);counts=Counter()
    for i,row in enumerate(source['charts']):
        coefficients=row['centre']['representative'];centre=vector(QQ,coefficients+[0]*(28-len(coefficients)))
        if len(coefficients)!=17:raise ArithmeticError('fixed generic centre changed')
        for k,P in enumerate(oracles):
            for sign in (1,-1):
                target=(centre/2-sign*projections.column(k))*inverse;gs_target=tuple(map(float,target*mu));trial=vector(ZZ,[round(float(x)) for x in target]);delta=trial-target;radius=float(delta*reduced*delta)+1;solutions=Enumeration(gso).enumerate(0,28,radius,0,target=gs_target)
                if not solutions:raise ArithmeticError('no CVP proposal')
                _,coordinates=solutions[0];small=vector(ZZ,[round(x) for x in coordinates])
                if any(abs(float(x)-int(y))>1e-7 for x,y in zip(coordinates,small)):raise ArithmeticError('nonintegral proposal')
                word=small*U;translated=sign*P+sum((int(c)*p for c,p in zip(word,basis) if c),E(0))
                if translated.is_zero():raise ArithmeticError('independent oracle became zero')
                point=[str(translated[0]),str(translated[1])];v=point_visibility(row['search'],point)
                if v['status'] in ('VISIBLE_NOT_RECORDED','UNSEARCHED_INTERVAL'):raise ArithmeticError('completed control discrepancy')
                result['proposals'].append({'chart_index':i,'oracle_index':k,'published_index':protocol['oracle_published_indices'][k],'sign':sign,'word':list(map(int,word)),'point':point,'visibility':v});counts[v['status']]+=1
        checkpoint(OUT,result)
        if (i+1)%5==0:print('TRANSLATED NATIVE29',i+1,dict(counts),flush=True)
    best=[min((r for r in result['proposals'] if r['oracle_index']==k),key=lambda r:r['visibility'].get('minimum_affine_height') or 10**1000) for k in range(1)]
    result.update(status='COMPLETE_DECLARED_ORACLE_AUDIT',status_counts=dict(counts),best_per_direction=best);checkpoint(OUT,result)
    print('BEST TRANSLATED NATIVE29',[(r['published_index'],r['chart_index'],r['visibility']['minimum_affine_height']) for r in best],flush=True)
if __name__=='__main__':main()
