#!/usr/bin/env sage-python
"""Independent equation/source-prefix check complementing the descent replay.

No producer/classifier import, seed packet, search artifact or later point.
The302 member remains calibrated; this does not certify prospective selection.
"""
import argparse,hashlib,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ART/'det1092_split_descent_v1'
OUT=D/'incidence-replay.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(b,s):
    if not b:raise ArithmeticError(s)
def point(E,p):return E(0) if p is None else E(list(map(QQ,p)))

def verify():
    protocol=read(D/'protocol.json');independent=read(D/'independent-replay.json')
    require(independent['status']=='PASS_INDEPENDENT_HALVING_CYCLE_REPLAY','descent replay')
    for path,h in independent['inputs'].items():require(sha(ROOT/path)==h,'descent binding')
    for path,h in protocol['inputs'].items():require(sha(ROOT/path)==h,'input binding')
    original=read(ART/'curve302_recovered_mw17_parent_v1.json')
    reduced=read(ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json')
    chart=read(ART/'det1092_reduced_parameter_chart_v1/generic-proof.json')
    require(chart['export_sha256']==sha(ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'),'reduced chart')
    conic=read(ART/'det1092_orbit8044_rank18_base_change_v2.json')
    carrier=read(ART/'det1092_norm8_seed_cover_v2/equation-only-cover-01.json')
    evaluation=read(ART/'det1092_blind_mw16_v1/evaluation.json')
    R=PolynomialRing(QQ,'t')
    def val(d,t):return R(d['numerator'])(t)/R(d['denominator'])(t)
    count=0;families=[]
    for family in ['302',*protocol['parameters']]:
        folder=D/family;frame=read(folder/'generic-frame.json')
        E=EllipticCurve(QQ,list(map(QQ,frame['curve'])))
        if family=='302':
            tau=QQ(0);old=EllipticCurve(QQ,[val(v,tau) for v in original['a_invariants']])
            require(E==EllipticCurve(QQ,[-old.c4()/48,-old.c6()/864]),'302 short equation')
            def forward(x,y):return E([x+old.b2()/12,y+(old.a1()*x+old.a3())/2])
            expected=[forward(val(x,0),val(y,0)) for x,y in original['basis_weierstrass_coordinates']]
        else:
            s=QQ(protocol['parameters'][family]);b=s.denominator()
            require(E.a_invariants()==tuple(val(v,s)*b**k for v,k in zip(reduced['a_invariants'],[2,4,6,8,12])),'reduced equation')
            expected=[E([val(x,s)*b**4,val(y,s)*b**6]) for x,y in reduced['basis_weierstrass_coordinates']]
            a,b0,c,d=map(QQ,chart['parameter_matrix']);tau=(a*s+b0)/(c*s+d)
            old=EllipticCurve(QQ,[val(v,tau) for v in original['a_invariants']])
            h=c*s.numerator()+d*s.denominator();scale=QQ(chart['weierstrass_u'])
        require(expected==[point(E,p) for p in frame['basis']],'literal generic17 specialization')
        branches=[point(E,read(folder/f'branch-{i}.json')['original']) for i in range(2)]
        for P in branches:
            if family=='302':
                m=carrier['maps'];cx0,cx1,cy0,cy1=[val(m[k],0) for k in ['x0','x1','y0','y1']]
                x=P[0]-old.b2()/12;y=P[1]-(old.a1()*x+old.a3())/2
                require(cx1!=0,'carrier coordinate patch');W=(x-cx0)/cx1
                require(W*W==R(carrier['quartic_coefficients'])(0) and y==cy0+cy1*W,'carrier incidence')
            else:
                x=P[0]*scale**2/h**4-old.b2()/12
                y=P[1]*scale**3/h**6-(old.a1()*x+old.a3())/2
                c0,c1,c2=[val(v,tau) for v in conic['lift']['residual_coefficients']]
                f0,f1,f2=[R(v)(tau) for v in conic['lift']['line_coefficients']]
                require(c2*x*x+c1*x+c0==0 and f0+f1*x+f2*y==0,'conic incidence')
            old([x,y]);count+=1
        trace=branches[0]+branches[1]
        expected_trace=expected[14]-expected[15] if family=='302' else sum((int(n)*P for n,P in zip(conic['lift']['trace_word'],expected)),E(0))
        require(trace==expected_trace,'inherited trace')
        for i,P in enumerate(branches):
            anti=point(E,read(folder/f'anti-trace-{i}.json')['original'])
            require(anti==2*P-trace,'anti-trace incidence');count+=1
        if family=='302':
            for arm in [r for r in evaluation['rows'] if r['success']]:
                raw=read(ART/'det1092_blind_mw16_v1/arms'/arm['arm_id']/'certificate.json')['candidates'][arm['first_success_candidate_id']]
                x,y=[val(v,0) for v in raw['coordinates']];P=forward(x,y)
                for suffix in ['core16','full17']:
                    saved=read(folder/(arm['arm_id']+'-'+suffix+'.json'))
                    require(point(E,saved['original'])==P,'old RR reconstruction incidence')
                    require(saved['omitted_basis_index']==(arm['omitted_basis_index_one_based']-1 if suffix=='core16' else None),'reference subgroup')
                    count+=1
        families.append(dict(family=family,generic_points=17,branch_points=2,anti_traces=2))
    require(count==38,'complete fixed roster')
    return dict(status='PASS_INDEPENDENT_SPLIT_INCIDENCE_AND_GENERIC_PREFIXES',decisions=count,
        families=families,checker_sha256=sha(Path(__file__)),
        descent_replay_sha256=sha(D/'independent-replay.json'),protocol_sha256=sha(D/'protocol.json'),
        boundary='Exact incidence and specialization check only. Descent replay independently certifies rational-span decisions. Calibrated302 member selection is not made prospective.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=verify();text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.write and not OUT.exists():
        with OUT.open('x') as f:f.write(text)
    else:require(OUT.read_text()==text,'immutable replay changed')
    print(result['status'],result['decisions'],flush=True)
