#!/usr/bin/env sage-python
"""Exact retained parent intake: full seed and good-reduction surface fingerprints."""
import sys
from pathlib import Path
from dataclasses import asdict
from sage.all import QQ,GF,PolynomialRing,EllipticCurve,legendre_symbol
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from memory_rank_certificate import checked_rank
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/mestre-parent-portfolio-intake-v2';OUT=ART/'mestre_parent_portfolio_intake_v1.json'

def model(data):
    R=PolynomialRing(QQ,'T');return R([QQ(v) for v in data['A_coefficients']]),R([QQ(v) for v in data['B_coefficients']])

def geometry(A,B):
    D=-16*(4*A**3+27*B**2);c4=-48*A
    assert A.degree()==8 and B.degree()==12 and D.gcd(c4).degree()==0
    factors=D.squarefree_decomposition();profile=sorted((int(f.degree()),int(m)) for f,m in factors)
    assert (D.degree()==20 and profile==[(2,2),(16,1)]) or (D.degree()==24 and profile==[(24,1)])
    rad=D//D.gcd(D.derivative())
    return D,rad,profile

def good_model(A,B,prime):
    D,rad,profile=geometry(A,B);F=GF(prime);R=PolynomialRing(F,'T')
    try:a=R([F(v) for v in A]);b=R([F(v) for v in B]);r=R([F(v) for v in rad])
    except (ValueError,ZeroDivisionError):return None
    disc=-16*(4*a**3+27*b**2)
    if a.degree()!=8 or b.degree()!=12 or disc.degree()!=D.degree() or r.degree()!=rad.degree() or r.gcd(r.derivative()).degree()!=0 or disc.gcd(a).degree()!=0:return None
    # Require the I4 infinity node to split; the X948 control has good infinity.
    if D.degree()==20:
        X=R.gen();f=X**3+a[8]*X+b[12];g=f.gcd(f.derivative())
        if g.degree()!=1:return None
        root=-g[0]/g[1]
        if legendre_symbol(int(3*root),prime)!=1:return None
    return a,b,disc,profile

def count_surface(A,B,prime):
    a,b,disc,profile=good_model(A,B,prime);F=GF(prime);total=0;repairs=[];smooth=0
    for t in list(F)+[None]:
        aa,bb=(a[8],b[12]) if t is None else (a(t),b(t))
        naive=prime+1+sum(int(legendre_symbol(int(z**3+aa*z+bb),prime)) for z in F)
        singular=4*aa**3+27*bb**2==0
        if not singular:
            assert int(EllipticCurve(F,[aa,bb]).cardinality())==naive;smooth+=1
        elif t is None:
            assert profile==[(2,2),(16,1)];naive+=3*prime;repairs.append({'base':'infinity','type':'split_I4','correction':3*prime})
        else:
            T=disc.parent().gen();multiplicity=0;tmp=disc
            while tmp(t)==0:tmp=tmp//(T-t);multiplicity+=1
            assert multiplicity in (1,2)
            if multiplicity==2:naive+=prime;repairs.append({'base':int(t),'type':'I2','correction':prime})
        total+=naive
    return {'prime':prime,'surface_point_count':total,'smooth_fibres_independently_counted':smooth,'resolution_corrections':repairs}

def main():
    if OUT.exists():raise FileExistsError('preserve parent certificate')
    paths=[Path(__file__).resolve(),D/'protocol.json',D/'ledger.json'];protocol=cert.read(D/'protocol.json');ledger=cert.read(D/'ledger.json')
    assert ledger['status']=='COMPLETE_DECLARED_ATTEMPTS' and all(r['supervision']['outcome']=='completed' and r['supervision']['returncode']==0 for r in ledger['rows'])
    rows=[];models=[]
    for metadata in protocol['rows']:
        path=D/metadata['id']/'intake.json';paths.append(path);d=cert.read(path);assert d['stage']=='complete';A,B=model(d);_,_,profile=geometry(A,B)
        curve=tuple(map(cert.F,d['curve']));halves=[tuple(map(cert.F,P)) for P in d['points']];origin=tuple(map(cert.F,d['covariant_images'][0]));points=[origin,*halves]
        E=EllipticCurve(QQ,[QQ(str(v)) for v in curve]);images=[E([QQ(v) for v in P]) for P in d['covariant_images']]
        for half,P in zip(halves,images[1:]):assert 2*E([QQ(str(v)) for v in half])==P-images[0]
        assert E.is_isomorphic(EllipticCurve(QQ,[0,0,0,A(1),B(1)]))
        cache=ReductionCache(MemoryFactStore());pivots={};sigs=[]
        for prime in _primes_up_to(997):
            if prime==2:continue
            try:s=signature(cache,curve,points,prime)
            except ValueError:continue
            before=len(pivots)
            for bits in s.rows:insert(pivots,bits)
            if len(pivots)>before:sigs.append(asdict(s))
        chosen=sorted(pivots);basis=[points[i] for i in chosen];proof=checked_rank(curve,basis,[s['prime'] for s in sigs],d['rank_certificate']['no_rational_2_torsion_prime'])
        rows.append({'id':metadata['id'],'outer_u':d['outer_u'],'fibre_T':'1','A_coefficients':d['A_coefficients'],'B_coefficients':d['B_coefficients'],'geometry_profile':profile,'curve':d['curve'],'points':[list(map(str,P)) for P in points],'independent_points':[list(map(str,P)) for P in basis],'independent_column_indices':chosen,'signatures':sigs,'rank_certificate':proof,'rank_lower_bound':len(basis),'relative_quartic_divisor_seed_lower_bound':d['rank_lower_bound']})
        models.append((A,B));print(metadata['id'],'full seed',len(basis),flush=True)
    atlas=ART/'compact_six_r17_atlas_v1.json';paths.append(atlas);r=next(r for r in cert.read(atlas)['families'] if r['family']=='11952');reference={'A_coefficients':r['A_coefficients_low_to_high'],'B_coefficients':r['B_coefficients_low_to_high']};models.append(model(reference))
    primes=[]
    for p in _primes_up_to(251):
        if p<101:continue
        if all(good_model(A,B,p) is not None for A,B in models):primes.append(p)
        if len(primes)==3:break
    if len(primes)!=3:raise ArithmeticError('three common good primes not found in fixed101..251 range')
    fingerprints=[[count_surface(A,B,p) for p in primes] for A,B in models]
    for r,f in zip(rows,fingerprints):r['surface_counts']=f
    separation=[]
    labels=[r['id'] for r in rows]+['X948-11952']
    for i in range(7):
        for j in range(i):
            differences=[p for p,a,b in zip(primes,fingerprints[i],fingerprints[j]) if a['surface_point_count']!=b['surface_point_count']]
            separation.append({'pair':[labels[j],labels[i]],'separating_good_primes':differences,'status':'PROVED_NOT_Q_ISOMORPHIC_K3' if differences else 'UNRESOLVED'})
    result={'schema':'elliptic-curves.mestre-parent-portfolio-intake.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'reference_parent':{'id':'X948-11952','surface_counts':fingerprints[-1]},'common_good_primes':primes,'pairwise_separation':separation,'point_search_boxes':0,'scope':'Six fixed outer-parameter K3 parents. Exact rational sections and14-point specialized seeds; lower bounds from exact finite quotients and torsion exclusion. Profile16 I1 +2 I2 +split I4 at infinity, checked overQQ and at common good primes. Counts include regular-resolution corrections and independently check every smooth fibre by Sage. Different good-reduction point counts prove distinct Q-isomorphism classes of smooth projective K3 surfaces, hence distinct Q-birational parents; no geometric-isomorphism or literature-novelty claim. No rank13 inference from the earlier mixed-specialization certificate, no score or elliptic point enumeration.'}
    checkpoint(OUT,result);print('PASS',primes,[(r['id'],r['rank_lower_bound']) for r in rows],sum(bool(r['separating_good_primes']) for r in separation),'of21 pairs separated',flush=True)

if __name__=='__main__':main()
