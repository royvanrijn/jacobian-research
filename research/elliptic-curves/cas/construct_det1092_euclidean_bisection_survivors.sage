#!/usr/bin/env sage-python
"""Three sealed survivors: Euclidean RR and factor-free exact incidence.

No point addition, nullspace, polynomial/integer factorization, new prime,
or point/parameter search. All traces and nine addresses were already saved.
"""
import hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_euclidean_bisection_survivors_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
PRE=ART/'det1092_modular_bisection_exclusion_v1'
SELECT=ART/'det1092_uniform_modular_rr_v2/summary.json'
CONTROLS=ART/'det1092_rational_bisection_index_v1/controls.json'
MASKS=[61,107,111]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    payload=json.dumps(data,indent=2,sort_keys=True)+'\n';path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as stream:stream.write(payload)
def coeff(f):return list(map(str,f.list()))
def square_certificate(q):
    q=QQ(q)
    if q<0:return {'square':False,'obstruction':'negative','value':str(q)}
    n,d=ZZ(q.numerator()),ZZ(q.denominator());a=n.isqrt();b=d.isqrt()
    square=a*a==n and b*b==d
    result={'square':bool(square),'value':str(q),'numerator':str(n),'denominator':str(d),
            'numerator_floor_sqrt':str(a),'denominator_floor_sqrt':str(b)}
    if square:result['root']=str(QQ(a)/b)
    return result
OUT.mkdir(exist_ok=True)
selection=read(SELECT);assert selection['new_unresolved']==MASKS
old_controls=read(CONTROLS)
roster=[{'label':r['label'],'parameter':r['parameter']} for r in old_controls['rows'][0]['cases']]
assert len(roster)==9 and len({r['parameter'] for r in roster})==9
assert sum(QQ(r['parameter'])==0 for r in roster)==1
save('protocol.json',{'classification':'generic equations for three previously unresolved orbits; retrospective unchanged-control evaluation',
    'masks':MASKS,'roster':roster,
    'limits':{'seconds':25,'orbits':3,'old_addresses':9,'incidence_tests':27,'new_parameters':0,
              'new_primes':0,'point_searches':0,'point_additions':0,'RR_matrix_solves':0,
              'polynomial_factorizations':0,'integer_factorizations':0,'full_atlas_runs':0},
    'rule':'Use exactly61,107,111, the remaining UNKNOWN masks of the fixed old panel. Construct each curve from its saved generic trace only. Then evaluate all nine unchanged addresses. No outcome-driven replacement.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,SELECT,CONTROLS,PRE/'summary.json',
         *[PRE/('trace-%d.json'%m) for m in MASKS],Path(__file__)]}})
begun=time.monotonic();parent=read(PARENT)
R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
def dec(r):return K(R(r['numerator']))/R(r['denominator'])
E=EllipticCurve(K,[dec(r) for r in parent['a_invariants']])
A=R(-E.c4()/48);B=R(-E.c6()/864)
G=matrix(ZZ,parent['generic_height_gram']);results=[]
for mask in MASKS:
    tr=read(PRE/('trace-%d.json'%mask));w=vector(ZZ,tr['word']);assert w*G*w==10
    h,nx,ny=[R(tr[k]) for k in ['pole','nx','ny']]
    assert h.is_monic() and h.degree()==3 and nx.degree()<=10 and ny.degree()<=15
    assert nx.gcd(h).degree()==0 and ny*ny==nx**3+A*nx*h**4+B*h**6
    f1=(-ny*nx.inverse_mod(h*h))%(h*h)
    f0,remainder=(-(f1*nx+ny)).quo_rem(h*h);assert not remainder
    f2=h
    branch_numerator=f1**4-6*nx*f1*f1-8*ny*f1-3*nx*nx-4*A*h**4
    q,remainder=branch_numerator.quo_rem(h**6);assert not remainder
    assert 1<=q.degree()<=2 and q.gcd(q.derivative()).degree()==0
    b,remainder=(f1*f1-nx).quo_rem(h*h);assert not remainder
    y0,remainder=(-(2*f0+f1*b)).quo_rem(2*h);assert not remainder
    x0=b/2;x1=h/2;y1=-f1/2
    assert y0*y0+y1*y1*q==x0**3+3*x0*x1*x1*q+A*x0+B
    assert 2*y0*y1==3*x0*x0*x1+x1**3*q+A*x1
    assert f0+f1*x0+h*y0==0 and f1*x1+h*y1==0
    curve={'status':'PASS_EUCLIDEAN_NORM10_CONIC','mask':mask,'word':tr['word'],
        'line':[coeff(f) for f in [f0,f1,h]],'branch_quadratic':coeff(q),
        'elliptic_maps':[coeff(f) for f in [x0,x1,y0,y1]],
        'map_convention':'X=x0+x1*W,Y=y0+y1*W,W^2=q(t); short parent coordinates',
        'degree_bounds':{'h':int(h.degree()),'nx':int(nx.degree()),'ny':int(ny.degree()),
                         'f0':int(f0.degree()),'f1':int(f1.degree()),'q':int(q.degree())},
        'trace_sha256':sha(PRE/('trace-%d.json'%mask))}
    save('orbit-%d.json'%mask,curve)
    exposures=[]
    for addr in roster:
        tt=QQ(addr['parameter']);assert -16*(4*A(tt)**3+27*B(tt)**2)
        cert=square_certificate(q(tt));row={**addr,'square_test':cert,
             'outcome':'RATIONAL_SPLIT_NEEDS_INDEPENDENCE' if cert['square'] else 'EXACT_NONSPLIT'}
        if cert['square']:
            rr=QQ(cert['root']);points=[]
            for sign in [1,-1]:
                xx=x0(tt)+sign*x1(tt)*rr;yy=y0(tt)+sign*y1(tt)*rr
                assert yy*yy==xx**3+A(tt)*xx+B(tt)
                points.append([str(xx),str(yy)])
            row['short_points']=points
        exposures.append(row)
    record={'mask':mask,'curve_sha256':sha(OUT/('orbit-%d.json'%mask)),'cases':exposures}
    save('incidence-%d.json'%mask,record);results.append(record)
    print('orbit',mask,'q degree',q.degree(),'rational splits',[(r['label'],r['parameter']) for r in exposures if r['square_test']['square']],flush=True)
summary={'status':'PASS_THREE_EXACT_EUCLIDEAN_CONICS_AND_27_INCIDENCES',
    'classification':'generic construction followed by retrospective control evaluation; independent replay required',
    'masks':MASKS,'trials':27,
    'splits':[{'mask':r['mask'],'label':v['label'],'parameter':v['parameter']} for r in results for v in r['cases'] if v['square_test']['square']],
    'nonsplit_count':sum(not v['square_test']['square'] for r in results for v in r['cases']),
    'boundary':'Only the three previously unresolved orbits, not a full atlas. A rational split is not automatically independent. No new prime or exceptional point input.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',*[OUT/('orbit-%d.json'%m) for m in MASKS],*[OUT/('incidence-%d.json'%m) for m in MASKS]]}}
save('summary.json',summary)
if not (OUT/'runtime.json').exists():save('runtime.json',{'seconds':time.monotonic()-begun,'scope':'three polynomial inverses and27 exact square decisions'})
print(summary['status'],'nonsplit',summary['nonsplit_count'],'seconds',round(time.monotonic()-begun,3),flush=True)
