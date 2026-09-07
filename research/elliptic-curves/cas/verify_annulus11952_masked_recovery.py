#!/usr/bin/env python3
"""Independent rational relations certify four withheld generic-direction recoveries."""
from math import lcm
from pathlib import Path
import annulus11952_masked_controls_v2 as control
import certify_compact_r17_candidates as cert
from half_lattice_pointed_sieve import linear_combination_python
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
ROOT=control.ROOT;ART=control.ART;INPUT=ART/'annulus11952_masked_relations_v2.json';OUT=ART/'annulus11952_masked_recovery_replay_v1.json'
def main():
    p=control.protocol();data=cert.read(INPUT);ledger=cert.read(control.D/'controller/ledger.json');prepared=cert.read(control.D/'prepared.json')
    if ledger['status']!='PASS' or data['status']!='COMPLETE_BOUNDED_RELATION_AUDIT' or len(data['rows'])!=4 or prepared['status']!='PASS' or prepared['protocol_hash']!=digest(p):raise ArithmeticError('all fixed masked searches, geometry and oracle relation stages required')
    if any(cert.hashed(ROOT/n)!=h for n,h in data['sources'].items()):raise ArithmeticError('oracle/result bindings changed')
    paths={Path(__file__).resolve(),INPUT,control.D/'protocol.json',control.D/'controller/ledger.json',control.D/'prepared.json',ROOT/'elliptic-curves/cas/half_lattice_pointed_sieve.py',ROOT/'elliptic-curves/cas/alternate_quartic_covers.py',ROOT/'elliptic-curves/cas/memory_rank_certificate.py'};rows=[]
    for source,row in zip(p['rows'],data['rows']):
        folder=control.D/source['id'];raw=cert.read(folder/'result.json');oracle=cert.read(folder/'oracle.json');blind=cert.read(folder/'blind.json');paths.update(folder/n for n in ('result.json','oracle.json','blind.json','maps.json'))
        if source['id']!=row['id'] or prepared['oracle_hashes'][row['id']]!=cert.hashed(folder/'oracle.json') or prepared['maps_hashes'][row['id']]!=cert.hashed(folder/'maps.json') or row['completed_boxes']!=12 or len(raw['charts'])!=12:raise ArithmeticError('fixed48 completed boxes and sealed oracle/map hashes required')
        model=tuple(map(cert.F,raw['curve']));generic=tuple(tuple(map(cert.F,P)) for P in oracle['original_generic_points']);proof=oracle['original_independence']
        got=checked_rank(model,generic,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if digest(got)!=digest(proof) or len(generic)!=17 or blind['points']!=[dict(zip(('x','y'),map(str,P))) for P in generic[1:]]:raise ArithmeticError('independent generic17 proof or masked16 differs')
        points=[];seen=set()
        for chart in raw['charts']:
            if chart['search']['status']!='bounded_search_complete':raise ArithmeticError('complete masked boxes required')
            for item in chart['search']['finite_curve_points']:
                P=(cert.F(item['x']),cert.F(item['y']));key=(P[0],abs(P[1]))
                if key not in seen:seen.add(key);points.append(P)
        witness=row['recovery_witness']
        if row['recovery_status']!='RECOVERED_KNOWN_DIRECTION' or witness is None:raise ArithmeticError('exact known-direction witness required')
        i=witness['returned_point_index'];coeff=tuple(map(cert.F,witness['coefficients']))
        if not 0<=i<min(128,len(points)) or tuple(map(cert.F,witness['point']))!=points[i] or len(coeff)!=17 or not coeff[0]:raise ArithmeticError('retained raw point or nonzero withheld coefficient differs')
        den=lcm(*(c.denominator for c in coeff));lhs=linear_combination_python(model,[points[i]],[den]);rhs=linear_combination_python(model,generic,[int(den*c) for c in coeff])
        if lhs!=rhs:raise ArithmeticError('independent exact rational group equality differs')
        rows.append({'id':row['id'],'parameter':row['parameter'],'known_generic_rank':17,'retained_generic_rank':16,'completed_boxes':12,'point':list(map(str,points[i])),'coefficients':list(map(str,coeff)),'withheld_coefficient':str(coeff[0]),'status':'RECOVERED_KNOWN_DIRECTION'})
    result={'schema':'elliptic-curves.annulus11952-masked-recovery-replay.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in sorted(paths)},'rows':rows,'completed_boxes':48,'claim_boundary':'Four fixed new-annulus indices recover a previously withheld known generic direction. Independent finite-quotient generic17 proofs and rational group equalities with nonzero withheld coefficients certify this endpoint. All original candidate searches remain unchanged. This does not measure incidence or sensitivity for exceptional directions, add a new curve rank or prove rank upper bounds.'}
    if OUT.exists():
        if cert.read(OUT)!=result:raise ArithmeticError('saved independent masked replay differs')
    else:checkpoint(OUT,result)
    print('INDEPENDENT FOUR NEW-ANNULUS MASKED RECOVERIES AND48 BOXES PASS',flush=True)
if __name__=='__main__':main()
