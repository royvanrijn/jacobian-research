#!/usr/bin/env sage-python
"""Test the opposite marked quartic point as an explicit third auxiliary direction."""
import argparse,json,sys
from dataclasses import asdict
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_retained_cloud_modl as odd
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/native-rank3-carrier-marked-point-v1';OUT=ART/'native_rank3_carrier_marked_point_v1.json'
INPUT=ART/'rank_jump_minimal_native_block_carrier_v1.json';PROOF=ART/'rank_jump_native_triple_carrier_verification_v1.json';HELPER=ROOT/'elkies-k3/scripts/prepare_r17_norm12_11952_alternate_v4_rank_one_bases.sage'
def sources():
    return {**odd.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,PROOF,HELPER)}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve marked-point protocol')
    v=cert.read(PROOF)
    if v['status']!='PASS' or v['auxiliary_Jacobian_exact_rank']!=3:raise ArithmeticError('rank3 carrier proof required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native-carrier-marked-point.v1','sources':sources(),'maximum_opposite_points':2,'moduli':[3,5],'prime_bound':997,'seconds_per_stage':90,'rss_bytes':1610612736,'gate':'The marked genus-one carrier has certified auxiliary rank3 but the recorded descent returns only two explicit points. The pointed quartic helper supplies the exact x-coordinate of the opposite marked rational point. Test its at most two rational lifts for a missing independent auxiliary direction.','scope':'Exact rational square roots and curve isomorphism only; no point search or additional descent. Combine the two recorded points with all opposite-point lifts and compute exact finite odd-prime quotient lower bounds. Odd torsion is excluded by a good-reduction witness; existing rational2-torsion is respected. A third independent point gives a full-rank rational subspace, not saturation or an integral Mordell-Weil basis. If the span remains two, retain that finite result without changing candidates.'})
def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('marked point sources changed')
    source=cert.read(INPUT);g=source['geometry'];desc=source['descent'];helper=SourceFileLoader('marked_carrier_helpers',str(HELPER)).load_module();R=PolynomialRing(QQ,'z');q=R(list(reversed(g['quartic_coefficients'])));v=q[0].sqrt();_,rebuilt,base,oppx,constants=helper.pointed_curve(list(q),QQ(0),v,'z')
    if rebuilt!=q:raise ArithmeticError('pointed quartic differs')
    J=EllipticCurve(QQ,g['minimal_Jacobian_model']);iso=base.isomorphism_to(J);opposite=base.lift_x(oppx,all=True)
    if not 1<=len(opposite)<=p['maximum_opposite_points']:raise ArithmeticError('finite marked-point roster differs')
    original=[J([QQ(x) for x in P]) for P in desc['points']];points=original+[iso(Q) for Q in opposite]
    a1,a2,a3,a4,a6=J.a_invariants();c2=a2+a1*a1/4;c4=a4+a1*a3/2;c6=a6+a3*a3/4
    model=tuple(cert.F(str(a)) for a in [0,0,0,c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27]);short=[]
    for P in points:
        X,Y=P.xy();short.append((cert.F(str(X+c2/3)),cert.F(str(Y+(a1*X+a3)/2))))
    if any(not cert.is_on_weierstrass_curve(model,P) for P in short):raise ArithmeticError('short auxiliary point transport differs')
    audits=[]
    for ell in p['moduli']:
        tp=odd.ml.find_no_rational_l_torsion_prime(model,modulus=ell);basis={};signatures=[];processed=[]
        for prime in odd.ml._primes_up_to(p['prime_bound']):
            if prime in (2,ell):continue
            try:sig=odd.finite.signature(model,short,prime,ell)
            except ValueError:continue
            before=len(basis)
            for row in sig.rows:odd.insert(basis,row,ell)
            if len(basis)>before:signatures.append(asdict(sig))
            processed.append(prime)
            if len(basis)>=3:break
        audits.append({'modulus':ell,'no_rational_ell_torsion_prime':tp,'finite_column_rank':len(basis),'independent_column_indices':sorted(basis),'signatures':signatures,'processed_primes':processed})
    lower=max(a['finite_column_rank'] for a in audits)
    if lower>3:raise ArithmeticError('finite bound conflicts with auxiliary upper bound')
    return {'schema':'elliptic-curves.native-carrier-marked-point-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'pointed_curve':list(map(str,base.a_invariants())),'opposite_x':str(oppx),'opposite_points':[list(map(str,P.xy())) for P in opposite],'Jacobian_model':list(map(str,J.a_invariants())),'Jacobian_points':[list(map(str,P.xy())) for P in points],'short_curve':list(map(str,model)),'short_points':[[str(v) for v in P] for P in short],'audits':audits,'rank_lower_bound':lower,'auxiliary_rank_upper_bound':3,'claim_boundary':p['scope']}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','build','check']);v=a.parse_args()
    if v.stage=='prepare':prepare()
    else:
        if v.stage=='build' and OUT.exists():raise FileExistsError('preserve marked-point audit')
        d=expected()
        if v.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('marked-point replay differs')
        else:checkpoint(OUT,d)
        print('MARKED AUXILIARY POINT SPAN',d['rank_lower_bound'],[(a['modulus'],a['finite_column_rank']) for a in d['audits']],flush=True)
