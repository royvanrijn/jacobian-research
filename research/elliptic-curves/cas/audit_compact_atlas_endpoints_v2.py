#!/usr/bin/env python3
"""Exact zero/infinity fibres omitted from the eleven compact nonzero parameter boxes."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.search_state import raw_state
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from memory_rank_certificate import checked_rank
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/compact-atlas-endpoints-v2'
ATLASES=[ART/'compact_six_r17_atlas_v1.json',ART/'compact_five_mw16_atlas_v1.json'];PRIOR=ART/'skew8_r17_results_v1.json';OUT=ART/'compact_atlas_endpoints_v2.json'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),*ATLASES,PRIOR,ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/research_runtime/search_state.py',ROOT/'elliptic-curves/cas/audit_compact_atlas_endpoints.py',ROOT/'artifacts/local/elliptic-curves/compact-atlas-endpoints-v1/protocol.json',ROOT/'artifacts/local/elliptic-curves/compact-atlas-endpoints-v1/run.supervisor.json']}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve endpoint audit protocol')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.compact-atlas-endpoints.v2','sources':sources(),'supersedes_failed_protocol':'artifacts/local/elliptic-curves/compact-atlas-endpoints-v1/protocol.json','correction':'The first audit failed before a rank certificate because some generic sections have nonconstant denominators. Preserve that failure. Evaluate rational-function Laurent leading terms exactly and transport pole sections to the group identity, checking compatible pole orders and leading terms.','endpoint_roster':['zero','infinity'],'maximum_fibres':22,'prime_bound':997,'wall_seconds':180,'rss_bytes':1073741824,'gate':'All eleven compact prospective rational-parameter boxes exclude zero and infinity. Audit these22 omitted endpoints directly from the family equations, without scoring, rank-labelled selection or parameter scanning. Equality against the current pinned catalogue and528 previous address-equations is checked after constructing the complete endpoint roster.','geometry':'At zero take constant terms. At infinity put t=1/v, X=v^4*x,Y=v^6*y; the integral fibre has short coefficients A_8,B_12. Transport rational sections by exact Laurent valuations and leading coefficients, including the identity for compatible poles. At infinity subtract degrees4 and6 from the coordinate growth orders. Check discriminants and every transported point exactly.','scope':'Only existing sections, bounded finite-reduction independence and equation comparisons. This is not a rational-point search, full rank, saturation, conductor or universal-novelty claim. Any endpoint point campaign requires a separate frozen protocol.'})
def coordinate(r,endpoint,weight):
    n=list(map(cert.F,r['numerator_coefficients_low_to_high']));d=list(map(cert.F,r['denominator_coefficients_low_to_high']))
    ni=[i for i,v in enumerate(n) if v];di=[i for i,v in enumerate(d) if v]
    if not di:raise ArithmeticError('zero section denominator')
    if not ni:return {'pole':False,'value':cert.F(0),'valuation':None,'leading':cert.F(0)}
    if endpoint=='zero':i,j=ni[0],di[0];valuation=i-j
    else:i,j=ni[-1],di[-1];valuation=weight+j-i
    leading=n[i]/d[j]
    return {'pole':valuation<0,'value':None if valuation<0 else leading if valuation==0 else cert.F(0),'valuation':valuation,'leading':leading}

def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('endpoint source binding changed')
    rows=[]
    for atlas in ATLASES:
        for f in cert.read(atlas)['families']:
            name=f.get('family',f.get('fibration_id'));A=list(map(cert.F,f['A_coefficients_low_to_high']));B=list(map(cert.F,f['B_coefficients_low_to_high']))
            for endpoint in p['endpoint_roster']:
                model=(cert.F(0),)*3+((A[0],B[0]) if endpoint=='zero' else (A[8] if len(A)>8 else cert.F(0),B[12] if len(B)>12 else cert.F(0)))
                row={'family':name,'endpoint':endpoint,'curve':list(map(str,model)),'source_atlas':str(atlas.relative_to(ROOT))}
                if not cert.weierstrass_invariants(model)['discriminant']:row['status']='SINGULAR_FIBRE';rows.append(row);continue
                points=[];transports=[]
                for i,section in enumerate(f['sections']):
                    x=coordinate(section['X'],endpoint,4);y=coordinate(section['Y'],endpoint,6)
                    if x['pole'] or y['pole']:
                        if not (x['pole'] and y['pole']) or 3*x['valuation']!=2*y['valuation'] or y['leading']**2!=x['leading']**3:raise ArithmeticError('incompatible section pole at smooth endpoint')
                        transports.append({'basis_index':i,'status':'POINT_AT_INFINITY','x_valuation':x['valuation'],'y_valuation':y['valuation'],'x_leading':str(x['leading']),'y_leading':str(y['leading'])});continue
                    xy=(x['value'],y['value'])
                    if not cert.is_on_weierstrass_curve(model,xy):raise ArithmeticError('endpoint section identity failed')
                    points.append(xy);transports.append({'basis_index':i,'status':'AFFINE_POINT','point':list(map(str,xy))})
                row['section_transports']=transports
                cache=QuotientOnlyReductionCache(MemoryFactStore());state=raw_state(model,points,cache=cache,prime_bound=p['prime_bound'])
                if state.no_two_torsion_prime is None:row.update(status='NO_RATIONAL_2_TORSION_WITNESS_UNKNOWN',generic_points=[list(map(str,P)) for P in points]);rows.append(row);continue
                proof=checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
                row.update(status='CERTIFIED_SPECIALIZED_SUBGROUP',generic_points=[list(map(str,P)) for P in points],points=[list(map(str,P)) for P in state.basis],rank_lower_bound=state.rank,rank_certificate=proof);rows.append(row)
    if len(rows)!=22:raise ArithmeticError('fixed22 endpoint roster differs')
    previous=cert.read(PRIOR);old=previous['previous_equations']+[{'address':PRIOR.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in previous['curves']];catalogue=previous['catalogue']
    for row in rows:
        if row['status']=='SINGULAR_FIBRE':continue
        row['catalogue_matches']=[r['id'] for r in catalogue['equations'] if cert.isomorphic(row['curve'],r['ainvs'])];row['previous_matches']=[r['address'] for r in old if cert.isomorphic(row['curve'],r['curve'])]
    pairs=[[j,i] for i,r in enumerate(rows) for j,q in enumerate(rows[:i]) if r['status']!='SINGULAR_FIBRE' and q['status']!='SINGULAR_FIBRE' and cert.isomorphic(r['curve'],q['curve'])]
    return {'schema':'elliptic-curves.compact-atlas-endpoints-result.v2','status':'COMPLETE_DECLARED_ENDPOINT_AUDIT','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'within_roster_isomorphic_pairs':pairs,'catalogue_snapshot_curves':catalogue['curve_count'],'catalogue_raw_sha256':catalogue['raw_sha256'],'previous_address_equations':len(old),'claim_boundary':p['scope']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('endpoint audit replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve endpoint audit')
            checkpoint(OUT,d)
        print('EXACT22 ENDPOINT AUDIT',[(r['family'],r['endpoint'],r['status'],r.get('rank_lower_bound'),r.get('catalogue_matches'),len(r.get('previous_matches',[]))) for r in d['rows']],flush=True)
