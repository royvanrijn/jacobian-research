#!/usr/bin/env sage-python
"""Bounded finite-quotient word proposals, verified by the exact rational group law."""
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
from sage.all import EllipticCurve,QQ,ZZ,GF,matrix,lcm,prime_range
import annulus11952_masked_controls as control
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from mod_l_reduction_independence import mod_l_reduction_signature
from research_runtime.store import checkpoint
from search_observability import point_visibility
OUT=control.ART/'annulus11952_masked_relations_v1.json'
def inputs():
    p=control.protocol();ledger=cert.read(control.D/'ledger.json');prepared=cert.read(control.D/'prepared.json')
    if ledger['status']!='PASS' or len(ledger['rows'])!=4:raise ArithmeticError('all four searches and map replays must finish before oracle access')
    paths=[Path(__file__).resolve(),CAS/'mod_l_reduction_independence.py',CAS/'memory_rank_certificate.py',CAS/'search_observability.py',control.D/'protocol.json',control.D/'ledger.json',control.D/'prepared.json'];cases=[]
    for row,terminal in zip(p['rows'],ledger['rows']):
        folder=control.D/row['id'];raw=cert.read(folder/'result.json');oracle=cert.read(folder/'oracle.json');blind=cert.read(folder/'blind.json')
        if terminal['id']!=row['id'] or terminal['result_sha256']!=cert.hashed(folder/'result.json') or prepared['oracle_hashes'][row['id']]!=cert.hashed(folder/'oracle.json'):raise ArithmeticError('frozen result/oracle differs')
        for step in terminal['stages']:
            s=step['supervision']
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('search or replay not terminal successful')
        model=tuple(map(cert.F,raw['curve']));generic=oracle['original_generic_points'];proof=oracle['original_independence'];checked_rank(model,[tuple(map(cert.F,P)) for P in generic],[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if len(generic)!=17 or blind['points']!=[dict(zip(('x','y'),P)) for P in generic[1:]]:raise ArithmeticError('withheld subgroup differs')
        seen=set();points=[]
        for ch in raw['charts']:
            for P in ch['search']['finite_curve_points']:
                key=cert.F(P['x']),abs(cert.F(P['y']))
                if key not in seen:seen.add(key);points.append([P['x'],P['y']])
        paths.extend(folder/n for n in ('result.json','oracle.json','blind.json','maps.json'))
        cases.append((row,raw,model,generic,points))
    return p,{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},cases
def relation(E,basis,point,word):
    d=lcm(c.denominator() for c in word)
    return d*point==sum((ZZ(d*c)*P for c,P in zip(word,basis) if c),E(0))
def main(check=False):
    p,sources,cases=inputs()
    if check:
        out=cert.read(OUT)
        if out['sources']!=sources or len(out['rows'])!=4:raise ArithmeticError('masked audit binding differs')
    else:
        if OUT.exists():raise FileExistsError('preserve masked relation audit')
        out={'schema':'elliptic-curves.annulus11952-masked-relations.v1','status':'RUNNING','sources':sources,'moduli':[3,5,7,11],'prime_bound':997,'maximum_returned_points_per_curve':128,'rational_denominator_proposals':[1,2,3,4],'rows':[],'endpoint':p['endpoint'],'claim_boundary':'At most first128 sign-distinct returned points per curve enter finite-quotient word proposals; exact rational group equality and a nonzero withheld coefficient certify recovery of the known direction. All other relations remain UNKNOWN. A changed16-point reference subgroup and four fixed new-annulus curves do not calibrate exceptional-direction sensitivity or population rank incidence. No new curve rank, upper bound or saturation is claimed.'}
    for index,(row,raw,model,generic,points) in enumerate(cases):
        E=EllipticCurve(QQ,[QQ(str(a)) for a in model]);basis=[E([QQ(v) for v in P]) for P in generic];found=[E([QQ(v) for v in P]) for P in points[:128]];witness=None;attempts=[]
        if check:
            old=out['rows'][index];witness=old['recovery_witness'];attempts=old['proposal_attempts']
            if witness:
                i=witness['returned_point_index'];word=[QQ(c) for c in witness['coefficients']]
                if not 0<=i<len(found) or len(word)!=17 or witness['point']!=points[i] or witness['exact_group_relation'] is not True or not word[0] or not relation(E,basis,found[i],word):raise ArithmeticError('exact withheld-direction relation differs')
        else:
            allpoints=tuple(tuple(map(cert.F,P)) for P in generic+points[:128])
            for modulus in out['moduli']:
                rows=[]
                for prime in prime_range(3,998):
                    try:s=mod_l_reduction_signature(model,allpoints,int(prime),modulus)
                    except ValueError:continue
                    rows.extend(s.rows)
                    if rows and matrix(GF(modulus),rows)[:,:17].rank()==17:break
                if not rows or matrix(GF(modulus),rows)[:,:17].rank()!=17:
                    attempts.append({'modulus':modulus,'status':'INSUFFICIENT_FINITE_RANK'});continue
                m=matrix(GF(modulus),rows);b=m[:,:17];attempts.append({'modulus':modulus,'status':'FULL_FINITE_GENERIC_RANK','last_prime':int(prime)})
                for i,P in enumerate(found):
                    try:c=b.solve_right(m.column(17+i))
                    except ValueError:continue
                    for denominator in out['rational_denominator_proposals']:
                        if denominator%modulus==0:continue
                        nums=[int(denominator*v) for v in c];word=[QQ(v if v<=modulus//2 else v-modulus)/denominator for v in nums]
                        if word[0] and relation(E,basis,P,word):
                            witness={'returned_point_index':i,'point':points[i],'coefficients':list(map(str,word)),'exact_group_relation':True};break
                    if witness:break
                if witness:break
        visibility=[]
        x,y=map(cert.F,generic[0])
        for i,ch in enumerate(raw['charts']):
            r=ch['search'];record={**r,**({'completed_denominator':p['height']} if r['status']=='bounded_search_complete' else {})}
            for sign in (-1,1):visibility.append({'chart':i,'sign':sign,'observation':point_visibility(record,(x,sign*y))})
        result={'id':row['id'],'family':row['family'],'parameter':row['parameter'],'attempted_charts':len(raw['charts']),'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in raw['charts']),'returned_points_up_to_sign':len(points),'points_considered_for_relations':min(128,len(points)),'recovery_status':'RECOVERED_KNOWN_DIRECTION' if witness else 'UNKNOWN','recovery_witness':witness,'proposal_attempts':attempts,'direct_oracle_visibility':visibility}
        if check:
            if result!=old:raise ArithmeticError('masked audit row differs')
        else:out['rows'].append(result);checkpoint(OUT,out)
        print('MASKED NEW11952 RELATION',row['id'],result['recovery_status'],len(points),'points',flush=True)
    if not check:out['status']='COMPLETE_BOUNDED_RELATION_AUDIT';checkpoint(OUT,out)
    elif out['status']!='COMPLETE_BOUNDED_RELATION_AUDIT':raise ArithmeticError('audit incomplete')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
