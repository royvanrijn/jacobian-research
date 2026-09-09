#!/usr/bin/env sage-python
"""Two frozen good primes, full 17-section Gram, same 15-orbit panel.

Four-chart theorem prototype; no new orbit, target, prime, or point search.
The old prime157 is retained only as a failed global-geometry regression.
"""
import csv, hashlib, json, signal, time
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_uniform_modular_rr_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
TABLE=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
OLD=ART/'det1092_modular_bisection_direct_v1'
PRIMES=[149,151];CHARTS=['infinity',1,2,3]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def coeff(f):return list(map(int,f.list()))
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
OUT.mkdir(exist_ok=True)
old=read(OLD/'summary.json');masks=sorted(old['new_masks']+old['regression_masks'])
save('protocol.json',{'classification':'generic full-frame geometry and uniform four-chart RR theorem, fixed-panel application',
    'limits':{'seconds':25,'good_primes':PRIMES,'geometry_only_regression_primes':[157],
              'new_orbits':0,'old_orbits':15,'old_target':0,'old_target_trials':30,
              'generic_Gram_checks':2,'charts_per_trace_max':4,'rational_RR_solves':0,
              'point_searches':0,'class_or_unit_groups':0,'full_atlas_runs':0},
    'chart_order':CHARTS,'selection':'Use precisely the previous15 masks and first two old primes. No replacement or new target. Test the old157 only for the global24I1 hypothesis.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,OLD/'summary.json',Path(__file__)]}})
begin=time.monotonic();parent=read(PARENT);G=matrix(ZZ,parent['generic_height_gram'])
with TABLE.open() as stream:
    words={int(r['orbit_mask']):list(map(int,r['parent_MW17_w'].split()))
           for r in csv.DictReader(stream,delimiter='\t') if r['category']=='rational'}
assert len(words)==40917 and len(masks)==15
results=[];geometry=[]
for p in PRIMES+[157]:
    field=GF(p);R=PolynomialRing(field,'t');t=R.gen();K=R.fraction_field()
    def dec(r):return K(R([field(QQ(c)) for c in r['numerator']]))/R([field(QQ(c)) for c in r['denominator']])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']])
    A=R(-E.c4()/48);B=R(-E.c6()/864);D=-16*(4*A**3+27*B**2)
    geom={'prime':p,'A':coeff(A),'B':coeff(B),'discriminant_degree':int(D.degree()),
          'discriminant_derivative_gcd':coeff(D.gcd(D.derivative()))}
    good=D.degree()==24 and D.gcd(D.derivative()).degree()==0
    geom['rootless_24I1_gate']=good
    if p==157:
        assert not good
        geometry.append(geom);continue
    assert good and D.gcd(A).degree()==0
    short=EllipticCurve(K,[A,B])
    basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,row);basis.append(short([x+E.b2()/12,y+(E.a1()*x+E.a3())/2]))
    def height(P):
        if not P:return 0
        x=P[0];num=x.numerator().degree();den=x.denominator().degree()
        assert den%2==0
        excess=max(0,num-den-4);assert excess%2==0
        return 4+den+excess
    heights=[height(P) for P in basis]
    pair_heights=[[i,j,height(basis[i]+basis[j])] for i in range(17) for j in range(i)]
    gram=matrix(QQ,17,17)
    for i,h in enumerate(heights):gram[i,i]=h
    for i,j,h in pair_heights:gram[i,j]=gram[j,i]=QQ(h-heights[i]-heights[j])/2
    assert gram==G
    geom['basis_heights']=heights;geom['pair_sum_heights']=pair_heights
    geom['Gram']=list(map(list,[[int(v) for v in row] for row in gram.rows()]))
    geom['Gram_matches_generic']=True;geometry.append(geom)
    for mask in masks:
        w=vector(ZZ,words[mask]);assert w*G*w==10
        trace=-sum((n*P for n,P in zip(w,basis)),short(0));assert height(trace)==10
        X,Y=trace.xy();trials=[]
        for chart in CHARTS:
            if chart=='infinity':
                x,y=X,Y;aa,bb=A,B;target=field(0)
            else:
                sub=K(chart+1/t)
                x=K(t**4)*X(sub);y=K(t**6)*Y(sub)
                aa=R(t**8*A(sub));bb=R(t**12*B(sub));target=-field(1)/chart
            degree=int(x.denominator().degree())
            trials.append({'chart':chart,'denominator_degree':degree})
            if degree==6:break
        assert degree==6
        h=R(x.denominator()).sqrt().monic();nx=R(x*h*h);ny=R(y*h**3)
        assert h.degree()==3 and nx.degree()<=10 and ny.degree()<=15
        columns=[term*t**j for bound,term in [(9,h**3),(5,nx*h),(3,ny)] for j in range(bound+1)]
        mat=matrix(field,19,20,lambda i,j:columns[j][i]);pivots=list(mat.pivots())
        assert len(pivots)==19
        kernel=mat.right_kernel().basis()[0]
        f0=R(list(kernel[:10]));f1=R(list(kernel[10:16]));f2=R(list(kernel[16:]))
        g0,g1,g2=f0(target),f1(target),f2(target)
        assert -16*(4*aa(target)**3+27*bb(target)**2)
        if h(target):
            assert g2
            x0=nx(target)/h(target)**2
            a=-g2*g2;b=g1*g1+a*x0;c=2*g0*g1-g2*g2*aa(target)+b*x0
            assert -c*x0==g0*g0-g2*g2*bb(target)
            discriminant=b*b-4*a*c
            mode='finite-trace-quadratic';residual=[int(c),int(b),int(a)]
        else:
            assert not g2 and g1
            xc=-g0/g1
            discriminant=xc**3+aa(target)*xc+bb(target)
            mode='zero-trace-vertical-line';residual=[int(-discriminant),0,1]
        if not discriminant:status='UNKNOWN_REPEATED_REDUCTION'
        elif discriminant.is_square():status='LOCALLY_SPLIT_SIMPLE_UNKNOWN_OVER_Q'
        else:status='EXCLUDED_NONSQUARE'
        result={'mask':mask,'prime':p,'word':words[mask],'charts':trials,'selected_chart':chart,
                'target_coordinate':int(target),'pole':coeff(h),'nx':coeff(nx),'ny':coeff(ny),
                'kernel':list(map(int,kernel)),'pivot_columns':pivots,
                'pivot_minor_determinant':int(mat.matrix_from_columns(pivots).det()),
                'mode':mode,'residual_coefficients':residual,'square_test_value':int(discriminant),'status':status}
        save('orbit-%d-p%d.json'%(mask,p),result);results.append(result)
save('geometry.json',{'status':'PASS_TWO_ROOTLESS_PRIMES_AND_FULL_GENERIC_GRAMS','primes':geometry})
summary={'status':'PASS_UNIFORM_FOUR_CHART_RR_PROTOTYPE','classification':'written uniform theorem plus fixed-panel verification',
    'new_masks':old['new_masks'],'regression_masks':old['regression_masks'],'primes':PRIMES,
    'new_excluded':[m for m in old['new_masks'] if any(r['mask']==m and r['status']=='EXCLUDED_NONSQUARE' for r in results)],
    'new_unresolved':[m for m in old['new_masks'] if not any(r['mask']==m and r['status']=='EXCLUDED_NONSQUARE' for r in results)],
    'chart_counts':{str(c):sum(r['selected_chart']==c for r in results) for c in CHARTS},
    'status_counts':{s:sum(r['status']==s for r in results) for s in sorted({r['status'] for r in results})},
    'pole_gate_recovered':[{'mask':r['mask'],'prime':r['prime'],'chart':r['selected_chart'],'status':r['status']} for r in results if r['selected_chart']!='infinity'],
    'boundary':'No complete atlas execution or rational split/seed construction. Prime157 is not discarded from the old experiment; it fails this stronger global theorem hypothesis.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',OUT/'geometry.json',
        *[OUT/('orbit-%d-p%d.json'%(r['mask'],r['prime'])) for r in results]]}}
save('summary.json',summary)
if not (OUT/'runtime.json').exists():save('runtime.json',{'seconds':time.monotonic()-begin,'scope':'two full17-section Grams and30 old-orbit trials, not full-atlas timing'})
print(summary['status'],summary['pole_gate_recovered'],'unresolved',summary['new_unresolved'],'seconds',round(time.monotonic()-begin,3),flush=True)
