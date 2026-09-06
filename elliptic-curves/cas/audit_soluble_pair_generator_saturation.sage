#!/usr/bin/env sage-python
"""Finite certificates for2,3,5 saturation of the fixed rank2 carrier subgroup."""
import argparse,json,sys
from dataclasses import asdict
from pathlib import Path
from sage.all import QQ,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_compact_r17_ambiguous as finite
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/soluble-pair-generator-saturation-v1';OUT=ART/'soluble_pair_generator_saturation_v1.json'

def sources():
    paths=[Path(__file__).resolve(),ART/'rank_jump_global_pair_solubility_v1.json',ART/'rank_jump_global_carrier_verification_v1.json',ART/'soluble_pair_carrier_height_v1.json',CAS/'audit_compact_r17_ambiguous.py',CAS/'mod_l_reduction_independence.py',CAS/'mod2_reduction_independence.py',CAS/'certify_compact_r17_candidates.py',CAS/'elliptic_candidate_record.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed carrier saturation protocol')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.soluble-pair-generator-saturation.v1','sources':sources(),'case':'observed_positive','moduli':[2,3,5],'prime_bound':997,'wall_seconds':120,'rss_bytes':1073741824,'scope':'The twelve fixed carrier images produced no new model within the360-bit gate. Check whether the recorded free rank2 subgroup together with its rational2-torsion is2,3,5-saturated using exact finite quotients through997. Include the torsion column at2; at odd moduli require an exact no-rational-ell-torsion reduction witness. This is a finite group audit with no division-point search, new parameter, descent, height sweep or claim of full saturation. A failed finite injectivity gate remains UNKNOWN.'})

def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen saturation audit inputs changed')
    original=next(r for r in cert.read(ART/'rank_jump_global_pair_solubility_v1.json')['rows'] if r['id']==p['case']);raw=original['geometry']['minimal_Jacobian_model'];E=EllipticCurve(QQ,raw)
    roots=E.division_polynomial(2).roots(QQ,multiplicities=False)
    if len(roots)!=1 or original['descent']['rank_lower_bound']!=2 or original['descent']['rank_upper_bound']!=2:raise ArithmeticError('fixed rank2 and one rational2-torsion direction required')
    x=roots[0];T=E(x,-(E.a1()*x+E.a3())/2)
    if not T or 2*T:raise ArithmeticError('exact rational2-torsion point required')
    points=[tuple(map(cert.F,P)) for P in original['descent']['points']];model=tuple(map(cert.F,raw));inv=cert.weierstrass_invariants(model);short=(cert.F(0),cert.F(0),cert.F(0),-inv['c4']/48,-inv['c6']/864)
    def transport(P):
        a,b=P;return a+inv['b2']/12,b+(model[0]*a+model[2])/2
    free=[transport(P) for P in points];torsion=transport(tuple(cert.F(str(v)) for v in T.xy()))
    if any(not cert.is_on_weierstrass_curve(short,P) for P in free+[torsion]):raise ArithmeticError('short-model transport differs')
    audits=[]
    for ell in p['moduli']:
        columns=free+[torsion] if ell==2 else free;rows=[];signatures=[];processed=[];tp=None
        for prime in finite.ml._primes_up_to(p['prime_bound']):
            if prime in (2,ell):continue
            try:sig=finite.signature(short,columns,prime,ell)
            except ValueError:continue
            processed.append(prime)
            if ell!=2 and sig.group_order%ell and tp is None:tp=prime
            if len(finite.pivots(rows+list(sig.rows),ell))>len(finite.pivots(rows,ell)):
                rows.extend(sig.rows);signatures.append(asdict(sig))
            if len(finite.pivots(rows,ell))==len(columns) and (ell==2 or tp is not None):break
        full=len(finite.pivots(rows,ell))==len(columns)
        if tp is not None and not finite.ml.no_rational_l_torsion_reduction_certificate(short,tp,ell):raise ArithmeticError('odd torsion witness differs')
        audits.append({'modulus':ell,'columns':[[str(v) for v in P] for P in columns],'finite_column_rank':len(finite.pivots(rows,ell)),'required_column_rank':len(columns),'no_rational_ell_torsion_prime':tp,'signatures':signatures,'processed_good_primes':processed,'status':'PROVED_SATURATED' if full and (ell==2 or tp is not None) else 'UNKNOWN'})
    return {'schema':'elliptic-curves.soluble-pair-generator-saturation-result.v1','status':'PASS_FIXED_AUDIT','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'Jacobian_model':raw,'free_generators':original['descent']['points'],'rational_2_torsion_point':[str(v) for v in T.xy()],'short_model':list(map(str,short)),'audits':audits,'lemma':'Let H contain the two free generators and the full rational2-torsion. At2, independence of both free columns and the torsion column forces all coefficients in2Q=aP+bR+cT to be even; the remaining2-torsion lies in H. At an odd ell, multiplication by ell is invertible on the rational2-torsion; independent free columns force divisibility of a,b, and the no-ell-torsion witness then puts Q in H. These arguments prove only the stated prime saturations.','claim_boundary':p['scope']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','build','check']);v=a.parse_args()
    if v.stage=='prepare':prepare()
    else:
        if v.stage=='build' and OUT.exists():raise FileExistsError('preserve saturation result')
        d=expected()
        if v.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('finite saturation replay differs')
        else:checkpoint(OUT,d)
        print('CARRIER SUBGROUP SATURATION',[(a['modulus'],a['status'],a['finite_column_rank']) for a in d['audits']],flush=True)
