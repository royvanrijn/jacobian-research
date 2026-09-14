"""Portable replay of three sealed next-direction target-location diagnostics.

No point search, minimisation, factorisation, selector repair or timing rerun.
Reuse certified primes, exact rational/Sturm bound replay and finite-rank proofs.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from verify_pointed_height_bounds import verify as verify_bound, require
from verify_pointed_minimal_neighbours import verify as verify_minimal
from verify_pointed_height_portfolio import run as verify_portfolio
from memory_rank_certificate import checked_rank
from research_runtime.store import digest

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elliptic-curves/height_model_next_direction_v1'


def read(p):return json.loads(Path(p).read_text())


def add(P,Q,A):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    slope=(3*x*x+A)/(2*y) if P==Q else (v-y)/(u-x)
    a=slope*slope-x-u
    return a,slope*(x-a)-y


def verify(folder):
    completion=read(folder/'completion.json')
    require(len(completion['arms'])==9,'incomplete arm ledger')
    require(sum(bool(r['success']) for r in completion['arms'])==3,'unexpected number of control gains')
    require(all(r['arm']=='v3_factor_free' for r in completion['arms'] if r['success']), 'unexpected successful arm')
    require(all(r['status']=='COLD_SUPPORT_UNRESOLVED' for r in completion['arms'] if not r['success']), 'changed failure boundary')
    outcomes=[]
    for rank in (27,29,30):
        case=folder/f'curve302-M{rank}';report=read(case/'report.json');gain=read(case/'rank-input.json')
        model=tuple(map(F,gain['curve']));points=[tuple(map(F,p)) for p in gain['points']]
        require(len(points)==rank+1 and list(map(str,points[-1]))==report['winning_point'],'wrong winning point')
        certificate=gain['rank_certificate']
        actual=checked_rank(model,points,[s['prime'] for s in gain['signatures']],certificate['no_rational_2_torsion_prime'])
        require(digest(actual)==digest(certificate),'rank certificate differs')
        x,y=points[-1];a,b=map(F,report['anchor']);A,B=model[3:]
        require(model[:3]==(0,0,0) and b*b==a**3+A*a+B,'wrong pointed anchor')
        R=add(add((x,y),(x,y),A),(a,-b),A)
        require(R is not None and list(map(str,R))==report['R_2P_minus_Q'],'2P-Q transport changed')
        Hx=max(abs(R[0].numerator),R[0].denominator)
        require(str(Hx)==report['x_height'],'wrong elliptic x-height')
        verify_minimal(case,ROOT)
        names=['factor_free','preconditioned_full','neighbours/neighbour-0','neighbours/neighbour-1']
        require(len(report['models'])==4,'incomplete model diagnostics')
        for name,row in zip(names,report['models']):
            packet=read(case/(name+'-bounds.json'));verify_bound(packet,ROOT)
            require(packet['input']==read(case/(name+'-input.json')),'model input differs')
            require(packet['input']['curve']==gain['curve'] and packet['anchor']==report['anchor'],'chart equation/anchor differs')
            n,d=map(int,row['coordinate']);require(math.gcd(n,d)==1 and d>=0,'nonprimitive parameter')
            u,v,w,z=map(F,packet['input']['mapping']['matrix'])
            require((u*n+v*d)*(x-a)==(w*n+z*d)*(y+b),'point slope transport differs')
            H=max(abs(n),abs(d));require(str(H)==row['parameter_height'],'parameter height differs')
            ev=lambda key:sum(F(c)*n**i*d**(4-i) for i,c in enumerate(packet[key]))
            N,D=ev('numerator'),ev('denominator')
            require(D and N/D==R[0],'quartic x-map differs')
            require(N.denominator==D.denominator==1,'nonintegral homogeneous map')
            g=math.gcd(N.numerator,D.numerator);C=int(packet['finite_gcd_divisor'])
            S=max(abs(N),abs(D))/H**4;L=F(packet['real']['lower']);bound=F(C)/L
            require(C%g==0 and max(abs(N),abs(D))/g==Hx,'cancellation identity differs')
            ratio=F(H**4,Hx)/bound
            require(0<ratio<=1 and ratio==L*g/(C*S),'height inequality/residual identity fails')
            require(row['D']==str(bound) and row['actual_over_sufficient_fourth_power']==str(ratio), 'diagnostic bound differs')
            require(row['finite_slack_factor']==str(F(C,g)) and row['real_slack_factor']==str(S/L),'local slack differs')
            require(row['inside_H125000']==(H<=125000) and row['guaranteed_at_H125000']==(F(Hx)<=125000**4/bound),'coverage flag differs')
        predicted=min(report['models'],key=lambda r:F(r['D']))
        observed=min(report['models'],key=lambda r:int(r['parameter_height']))
        require(report['best_bound_model']==predicted['model'] and report['smallest_actual_height_model']==observed['model'],'model ranking differs')
        if report['joint_selection'] is None:
            require(len(set(report['neighbour_primes']))==1 and 'same-prime' in report['joint_unknown_reason'],'unjustified joint UNKNOWN')
        else:
            from height_model_benchmark import selection_from_portfolio
            result=verify_portfolio(case,ROOT)
            choice=selection_from_portfolio(result,[read(case/(name+'-bounds.json')) for name in names[1:]])
            require(choice==report['joint_selection'],'joint selector differs')
        outcomes.append({'case':report['case'],'rank_lower_bound':rank+1,
            'parameter_heights':[r['parameter_height'] for r in report['models']],
            'best_bound_matches_smallest_height':predicted['parameter_height']==observed['parameter_height']})
    require(sum(r['best_bound_matches_smallest_height'] for r in outcomes)==1,'changed target-order result')
    return {'status':'PASS_PORTABLE_TARGET_HEIGHT_AND_RANK_REPLAY','rows':outcomes,
        'boundary':'Three literal winning points on one curve. Global sufficient bounds need not rank individual '
                   'target accessibility. Six cold-preparation timeouts made no point calls; no comparative CPU '
                   'claim follows. Rank replay shares one finite implementation; independent Sage receipts retained.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,default=DEFAULT)
    a=p.parse_args();print(json.dumps(verify(a.folder.resolve()),indent=2))
