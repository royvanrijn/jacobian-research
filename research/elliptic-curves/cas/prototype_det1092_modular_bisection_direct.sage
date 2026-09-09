#!/usr/bin/env sage-python
"""Same frozen15x3 panel; compute traces only over finite function fields.

A maximal denominator degree certifies integral lifting from the Gauss DVR.
Reject degree drops. No rational point addition or RR solve.25s cap.
"""
import csv,hashlib,json,signal,time
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_modular_bisection_direct_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
TABLE=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
CENSUS=ART/'curve302_parent_degree2_multisection_lattice_v1.json'
PRIMES=[149,151,157];REGRESSION=[8044,47755,103186]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def coeff(f):return list(map(int,f.list()))
OUT.mkdir(exist_ok=True)
save('protocol.json',{'classification':'same-panel finite-field-only method with a degree-preserving lifting gate',
    'selection':'Unchanged first12 saved norm10 masks plus8044,47755,103186 regressions. No outcome-driven replacement. Target original parameter0 only.',
    'prime_pool':PRIMES,'limits':{'seconds':25,'new_orbits':12,'regression_orbits':3,'prime_trials':45,
        'rational_point_additions':0,'rational_RR_solves':0,'point_searches':0,'complete_atlas_runs':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,CENSUS,Path(__file__)]}})
begun=time.monotonic();parent=read(PARENT);census=read(CENSUS)
assert census['orbits_tsv_sha256']==sha(TABLE)
stored={}
with TABLE.open() as stream:
    for row in csv.DictReader(stream,delimiter='\t'):
        if row['category']=='rational':stored[int(row['orbit_mask'])]=list(map(int,row['parent_MW17_w'].split()))
assert len(stored)==40917
selected=sorted(stored)[:12];masks=sorted(set(selected+REGRESSION));assert len(masks)==15
G=matrix(ZZ,parent['generic_height_gram'])
for mask in masks:
    w=vector(ZZ,stored[mask]);assert w*G*w==10
outcomes={mask:[] for mask in masks}
for p in PRIMES:
    field=GF(p);R=PolynomialRing(field,'t');t=R.gen();F=R.fraction_field()
    def dec(record):
        nn=[field(QQ(v)) for v in record['numerator']]
        dd=[field(QQ(v)) for v in record['denominator']]
        assert R(dd)
        return F(R(nn))/R(dd)
    E=EllipticCurve(F,[dec(v) for v in parent['a_invariants']])
    A=R(-E.c4()/48);B=R(-E.c6()/864)
    assert A.degree()<=8 and B.degree()<=12 and -16*(4*A**3+27*B**2)
    short=EllipticCurve(F,[A,B]);basis=[]
    for point in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,point);basis.append(short([x+E.b2()/12,y+(E.a1()*x+E.a3())/2]))
    for mask in masks:
        w=stored[mask];row={'prime':p}
        trace=-sum((n*P for n,P in zip(w,basis)),short(0))
        if not trace or trace[0].denominator().degree()!=6:
            row['status']='UNKNOWN_MAXIMAL_POLE_GATE';outcomes[mask].append(row);continue
        h=R(trace[0].denominator()).sqrt().monic()
        nx=R(trace[0]*h*h);ny=R(trace[1]*h**3)
        assert h.degree()==3 and nx.degree()<=10 and ny.degree()<=15
        row['pole']=coeff(h);row['nx']=coeff(nx);row['ny']=coeff(ny)
        if not h(0) or not -16*(4*A(0)**3+27*B(0)**2):
            row['status']='UNKNOWN_BAD_TARGET_OR_TRACE';outcomes[mask].append(row);continue
        columns=[term*t**i for bound,term in [(9,h**3),(5,nx*h),(3,ny)] for i in range(bound+1)]
        mat=matrix(field,19,20,lambda i,j:columns[j][i]);pivots=list(mat.pivots())
        if len(pivots)!=19:
            row['status']='UNKNOWN_RR_RANK_DROP';outcomes[mask].append(row);continue
        kernel=mat.right_kernel().basis()[0]
        row['kernel']=list(map(int,kernel));row['pivot_columns']=pivots
        row['pivot_minor_determinant']=int(mat.matrix_from_columns(pivots).det())
        f0=R(list(kernel[:10]));f1=R(list(kernel[10:16]));f2=R(list(kernel[16:]))
        if not f2(0):
            row['status']='UNKNOWN_AFFINE_CHART';outcomes[mask].append(row);continue
        x0=nx(0)/h(0)**2;y0=ny(0)/h(0)**3
        assert f0(0)+f1(0)*x0+f2(0)*y0==0
        a=-f2(0)**2;b=f1(0)**2+a*x0;c=2*f0(0)*f1(0)-f2(0)**2*A(0)+b*x0
        assert -c*x0==f0(0)**2-f2(0)**2*B(0)
        disc=b*b-4*a*c
        row['residual_coefficients']=[int(c),int(b),int(a)];row['discriminant']=int(disc)
        row['status']='EXCLUDED_NONSQUARE' if disc and not disc.is_square() else 'UNKNOWN_SPLIT_OR_REPEATED_REDUCTION'
        outcomes[mask].append(row)
for mask in masks:
    save('orbit-%d.json'%mask,{'mask':mask,'word':stored[mask],'trials':outcomes[mask],
        'excluded':any(row['status']=='EXCLUDED_NONSQUARE' for row in outcomes[mask])})
summary={'status':'PASS_FINITE_FIELD_ONLY_BISECTION_PROTOTYPE',
    'classification':'exact modular exclusion method; independent replay and written lifting proof required',
    'new_masks':selected,'regression_masks':REGRESSION,'primes':PRIMES,
    'new_excluded':[i for i in selected if any(r['status']=='EXCLUDED_NONSQUARE' for r in outcomes[i])],
    'new_unknown':[i for i in selected if not any(r['status']=='EXCLUDED_NONSQUARE' for r in outcomes[i])],
    'regression_excluded':[i for i in REGRESSION if any(r['status']=='EXCLUDED_NONSQUARE' for r in outcomes[i])],
    'status_counts':{s:sum(r['status']==s for v in outcomes.values() for r in v) for s in sorted({r['status'] for v in outcomes.values() for r in v})},
    'boundary':'Same45 trial panel only. Split reductions and all failed lifting/chart gates are UNKNOWN. No rational trace addition, rational RR solve, point search or full atlas run.',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,CENSUS,Path(__file__),*[OUT/('orbit-%d.json'%i) for i in masks]]}}
save('summary.json',summary)
seconds=time.monotonic()-begun
if not (OUT/'runtime.json').exists():save('runtime.json',{'seconds':seconds,'orbits':15,'trials':45,'classification':'small-sample timing, not a full-atlas guarantee'})
print(summary['status'],'excluded',summary['new_excluded'],'unknown',summary['new_unknown'],
      'counts',summary['status_counts'],'seconds',round(seconds,3),flush=True)
