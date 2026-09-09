#!/usr/bin/env sage-python
"""Predict inherited conic splits from generic unit intersections only.

One old generic orbit; exactly578 signed words on at most two generators.
No control parameter, exceptional point, catalogue or search output input.
Freeze before construct. Every collision is checkpointed. Each action <=25s.
"""
import argparse,csv,hashlib,itertools,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,gcd
from sage.env import SAGE_VERSION
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_unit_intersection_obstructions_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
TABLE=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
CHART=ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def coeff(f):return list(map(str,f.list()))
def save(name,row):
    p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert read(p)==row,'immutable checkpoint changed'
    else:
        with p.open('x') as s:json.dump(row,s,indent=2,sort_keys=True);s.write('\n')
def freeze():
    save('protocol.json',dict(classification='generic construction of forced dependent splits; retrospective rule development',
        orbit=8044,selector='All signed words with support1 or2 and nonzero coefficients +/-1; retain exactly norm4 and pairing3 with w. No rational parameter in selection.',
        limits=dict(seconds_per_action=25,lattice_words=578,constructed_intersections_cap=16,
                    control_inputs=0,point_searches=0,parameter_scans=0,full_atlas_runs=0),
        sage_version=SAGE_VERSION,
        boundary='An exact finite subset of inherited split incidences, not all dependent specializations. '
                 'Support is basis-dependent; the intersection criterion is intrinsic. '
                 'Control comparison is a separate post-construction action.',
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,CHART,Path(__file__)]}))
def construct():
    protocol=read(OUT/'protocol.json')
    for p,h in protocol['inputs'].items():assert sha(ROOT/p)==h
    parent=read(PARENT);G=matrix(ZZ,parent['generic_height_gram'])
    with TABLE.open() as stream:
        row=next(v for v in csv.DictReader(stream,delimiter='\t') if int(v['orbit_mask'])==protocol['orbit'])
    assert row['category']=='rational'
    w=vector(ZZ,row['parent_MW17_w'].split());assert w*G*w==10
    selected=[];exposure=[]
    for size in [1,2]:
        for support in itertools.combinations(range(17),size):
            for signs in itertools.product([-1,1],repeat=size):
                v=vector(ZZ,17)
                for i,sgn in zip(support,signs):v[i]=sgn
                norm,pair=v*G*v,w*G*v
                exposure.append(dict(word=list(map(int,v)),norm=int(norm),pairing=int(pair),selected=bool(norm==4 and pair==3)))
                if norm==4 and pair==3:selected.append(v)
    assert len(exposure)==578 and 0<len(selected)<=protocol['limits']['constructed_intersections_cap']
    save('selection.json',dict(word=list(map(int,w)),selected=[list(map(int,v)) for v in selected],
        exposures=exposure,protocol_sha256=sha(OUT/'protocol.json')))
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    original=EllipticCurve(K,[dec(v) for v in parent['a_invariants']])
    A,B=R(-original.c4()/48),R(-original.c6()/864);E=EllipticCurve(K,[A,B])
    basis=[]
    for raw in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,raw);basis.append(E([x+original.b2()/12,y+(original.a1()*x+original.a3())/2]))
    trace=-sum((n*P for n,P in zip(w,basis)),E(0))
    h=R(trace[0].denominator()).sqrt().monic();N,V=R(trace[0]*h**2),R(trace[1]*h**3)
    assert h.degree()==3 and gcd(h,N)==1
    def exactdiv(a,b):
        v,r=a.quo_rem(b);assert not r;return v
    m=(-V*N.inverse_mod(h*h))%(h*h);g=exactdiv(m*N+V,h*h)
    b=exactdiv(m*m-N,h*h);k=exactdiv(m*b-2*g,h*h)
    q=exactdiv(4*m*k-3*b*b-4*A,h*h)
    assert q.degree()==2 and gcd(q,q.derivative())==1
    x0,x1,y0,y1=b/2,h/2,-h*k/2,-m/2
    one,U,V0=h.xgcd(m);assert one==1 and U*h+V0*m==1
    save('conic.json',dict(h=coeff(h),nx=coeff(N),ny=coeff(V),m=coeff(m),g=coeff(g),b=coeff(b),k=coeff(k),
        q=coeff(q),maps=[coeff(v) for v in [x0,x1,y0,y1]],bezout=[coeff(U),coeff(V0)],
        convention='Short original parent; W^2=q(t). Original parameter infinity uses X/t^4,Y/t^6,W/t.'))
    aa,bb,cc,dd=map(QQ,read(CHART)['parameter_matrix']);results=[]
    for index,v in enumerate(selected):
        P=sum((n*P for n,P in zip(v,basis)),E(0));X,Y=R(P[0]),R(P[1])
        assert X.degree()<=4 and Y.degree()<=6
        lift=U*(2*X-b)-V0*(2*Y+h*k)
        tests=[2*X-b-h*lift,2*Y+h*k+m*lift,lift*lift-q]
        divisor=gcd(tests).monic();assert divisor.degree() in [0,1]
        if divisor.degree()==1:
            tau=-divisor[0];omega=lift(tau)
            assert omega*omega==q(tau)
            assert X(tau)==x0(tau)+x1(tau)*omega and Y(tau)==y0(tau)+y1(tau)*omega
            reduced=None if cc*tau==aa else (bb-dd*tau)/(cc*tau-aa)
            branch=[X(tau),Y(tau)];other=[x0(tau)-x1(tau)*omega,y0(tau)-y1(tau)*omega]
            smooth=bool(4*A(tau)**3+27*B(tau)**2)
            if smooth:
                es=EllipticCurve(QQ,[A(tau),B(tau)])
                ev=[es([z(tau) for z in P.xy()]) for P in basis]
                assert es(branch)+es(other)==sum((n*P for n,P in zip(w,ev)),es(0))
        else:
            tau=None;omega=(X[4]-x0[4])/x1[3]
            assert omega*omega==q[2] and Y[6]==y0[6]+y1[5]*omega
            reduced=-dd/cc if cc else None
            branch=[X[4],Y[6]];other=[x0[4]-x1[3]*omega,y0[6]-y1[5]*omega]
            smooth=bool(4*A[8]**3+27*B[12]**2)
        result=dict(index=index,word=list(map(int,v)),companion_word=list(map(int,w-v)),
            norm=4,pairing=3,intersection_number=1,
            section=[coeff(X),coeff(Y)],ordinate_lift=coeff(lift),intersection_divisor=coeff(divisor),
            original_parameter=None if tau is None else str(tau),
            reduced_parameter=None if reduced is None else str(reduced),omega=str(omega),
            branch=list(map(str,branch)),companion=list(map(str,other)),smooth_fibre=smooth,
            conclusion='BOTH_BRANCHES_INHERITED' if smooth else 'INTERSECTION_OVER_SINGULAR_FIBRE')
        save('intersection-%02d.json'%index,result);results.append(result)
        print(index,'s =',result['reduced_parameter'],result['conclusion'],flush=True)
    summary=dict(status='PASS_GENERIC_UNIT_INTERSECTION_OBSTRUCTIONS',word=list(map(int,w)),
        count=len(results),distinct_parameters=len({r['original_parameter'] for r in results}),
        cases=[{k:r[k] for k in ['index','word','companion_word','original_parameter','reduced_parameter','conclusion']} for r in results],
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'selection.json',OUT/'conic.json',
            *[OUT/('intersection-%02d.json'%i) for i in range(len(results))]]},
        boundary='No exceptional coordinates or control parameters read. Exact inherited collisions, not a complete dependence locus or new seeds.')
    save('summary.json',summary)
    print(summary['status'],len(results),'intersections',summary['distinct_parameters'],'parameters',flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','construct']);args=parser.parse_args()
    signal.alarm(25);begun=time.monotonic()
    if args.action=='freeze':freeze()
    else:construct()
    print('seconds',round(time.monotonic()-begun,3),flush=True)
