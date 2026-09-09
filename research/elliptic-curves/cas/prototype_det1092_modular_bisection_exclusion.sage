#!/usr/bin/env sage-python
"""Frozen15-orbit modular RR prototype at three primes, no rational RR solve.

First12 saved norm10 masks plus three completed-map regressions. Target302
is used only for incidence evaluation. No point search or full atlas run.
"""
import csv,hashlib,json,signal,time,traceback
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_modular_bisection_exclusion_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
TABLE=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
CENSUS=ART/'curve302_parent_degree2_multisection_lattice_v1.json'
REG=ART/'det1092_rational_bisection_index_v1'
PRIMES=[149,151,157];REGRESSION=[8044,47755,103186]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rec(f):return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
def main():
    begun=time.monotonic();parent=read(PARENT);census=read(CENSUS)
    assert census['orbits_tsv_sha256']==sha(TABLE)
    stored={}
    with TABLE.open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            if row['category']=='rational':stored[int(row['orbit_mask'])]=list(map(int,row['parent_MW17_w'].split()))
    assert len(stored)==40917
    selected=sorted(stored)[:12];all_masks=sorted(set(selected+REGRESSION))
    assert len(all_masks)==15
    save('selection.json',{'rule':'First12 rational rows by saved mask, plus8044,47755,103186 as old exact-map regressions.',
        'new_masks':selected,'regression_masks':REGRESSION,
        'words':{str(i):stored[i] for i in all_masks},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [TABLE,CENSUS]}})
    R=PolynomialRing(QQ,'t');t=R.gen();F=R.fraction_field()
    def dec(v):return F(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(F,[dec(v) for v in parent['a_invariants']])
    A=R(-E.c4()/48);B=R(-E.c6()/864);short=EllipticCurve(F,[A,B])
    basis=[]
    for P in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,P);basis.append(short([x+E.b2()/12,y+(E.a1()*x+E.a3())/2]))
    G=matrix(ZZ,parent['generic_height_gram']);outcomes=[];trace_seconds=0;rr_seconds=0
    for mask in all_masks:
        tick=time.monotonic();word=vector(ZZ,stored[mask]);assert word*G*word==10
        trace=-sum((n*P for n,P in zip(word,basis)),short(0));assert trace
        h=R(trace[0].denominator()).sqrt()
        if h.leading_coefficient()<0:h=-h
        nx=R(trace[0]*h*h);ny=R(trace[1]*h**3)
        assert h.degree()<=3 and nx.degree()<=10 and ny.degree()<=15 and h(0)
        columns=[term*t**i for bound,term in [(9,h**3),(5,nx*h),(3,ny)] for i in range(bound+1)]
        assert all(col.degree()<=18 for col in columns)
        interpolation=matrix(QQ,19,20,lambda i,j:columns[j][i])
        trace_seconds+=time.monotonic()-tick
        trace_record={'mask':mask,'word':list(map(int,word)),
            'short_trace':[rec(v) for v in trace.xy()],
            'pole':list(map(str,h.list())),'nx':list(map(str,nx.list())),'ny':list(map(str,ny.list()))}
        save('trace-%d.json'%mask,trace_record)
        trials=[]
        for p in PRIMES:
            tick=time.monotonic();row={'prime':p}
            values=[*A.list(),*B.list(),*h.list(),*nx.list(),*ny.list(),*interpolation.list()]
            if any(c.denominator()%p==0 for c in values):
                row['status']='UNKNOWN_COEFFICIENT_DENOMINATOR';trials.append(row);continue
            field=GF(p);Rp=PolynomialRing(field,'t');hp=Rp(h)
            if not hp(0) or not field(-16*(4*A(0)**3+27*B(0)**2)):
                row['status']='UNKNOWN_BAD_TARGET_OR_TRACE';trials.append(row);continue
            mat=matrix(field,interpolation);pivots=list(mat.pivots())
            if len(pivots)!=19:
                row['status']='UNKNOWN_RR_RANK_DROP';trials.append(row);continue
            kernel=mat.right_kernel().basis()[0]
            f0=Rp(list(kernel[:10]));f1=Rp(list(kernel[10:16]));f2=Rp(list(kernel[16:]))
            row['kernel']=list(map(int,kernel));row['pivot_columns']=pivots
            row['pivot_minor_determinant']=int(mat.matrix_from_columns(pivots).det())
            if not f2(0):
                row['status']='UNKNOWN_AFFINE_CHART';trials.append(row);continue
            x0=Rp(nx)(0)/hp(0)**2;y0=Rp(ny)(0)/hp(0)**3
            assert f0(0)+f1(0)*x0+f2(0)*y0==0
            a=-f2(0)**2;b=f1(0)**2+a*x0
            c=2*f0(0)*f1(0)-f2(0)**2*field(A(0))+b*x0
            assert -c*x0==f0(0)**2-f2(0)**2*field(B(0))
            disc=b*b-4*a*c
            row['residual_coefficients']=[int(c),int(b),int(a)];row['discriminant']=int(disc)
            row['status']='EXCLUDED_NONSQUARE' if disc and not disc.is_square() else 'UNKNOWN_SPLIT_OR_REPEATED_REDUCTION'
            if mask in REGRESSION and disc:
                old=read(REG/('orbit-%d.json'%mask));q0=QQ(old['q'][0]);valuation=q0.valuation(p)
                assert valuation%2==0
                unit=q0/QQ(p)**valuation
                expected=field(unit).is_square()
                assert expected==disc.is_square()
                row['old_exact_squareclass_regression']='PASS'
            rr_seconds+=time.monotonic()-tick;trials.append(row)
        report={'mask':mask,'word':stored[mask],
                'excluded':any(row['status']=='EXCLUDED_NONSQUARE' for row in trials),
                'trials':trials,'trace_sha256':sha(OUT/('trace-%d.json'%mask))}
        save('orbit-%d.json'%mask,report);outcomes.append(report)
    summary={'status':'PASS_BOUNDED_MODULAR_RR_PROTOTYPE',
        'classification':'exact modular non-incidence certificates, not a full-atlas claim',
        'selected_new_masks':selected,'regression_masks':REGRESSION,'primes':PRIMES,
        'new_excluded':[r['mask'] for r in outcomes if r['mask'] in selected and r['excluded']],
        'new_unknown':[r['mask'] for r in outcomes if r['mask'] in selected and not r['excluded']],
        'regression_excluded':[r['mask'] for r in outcomes if r['mask'] in REGRESSION and r['excluded']],
        'trial_count':sum(len(r['trials']) for r in outcomes),
        'status_counts':{s:sum(row['status']==s for r in outcomes for row in r['trials']) for s in sorted({row['status'] for r in outcomes for row in r['trials']})},
        'boundary':'Only first12 new masks and three old regressions tested. A split reduction is UNKNOWN over Q. No rational RR kernel, bisection equation, point search or complete40917-orbit run was computed.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,CENSUS,OUT/'selection.json',Path(__file__),*[OUT/('orbit-%d.json'%i) for i in all_masks]]}}
    save('summary.json',summary)
    runtime={'total_seconds':time.monotonic()-begun,'trace_and_matrix_seconds':trace_seconds,
             'modular_RR_seconds':rr_seconds,'orbits':len(all_masks),'trials':45,
             'classification':'single small timing observation, not a full-run runtime guarantee'}
    if not (OUT/'runtime.json').exists():save('runtime.json',runtime)
    print(summary['status'],'new exclusions',summary['new_excluded'],'unknown',summary['new_unknown'],
          'counts',summary['status_counts'],'seconds',round(runtime['total_seconds'],3),flush=True)
if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    save('protocol.json',{'classification':'bounded generic-selected modular bisection certificate prototype',
        'selection':'First12 saved rational masks numerically, without target data; add exactly three completed-map regressions8044,47755,103186. Original parameter0 is used only for target incidence.',
        'prime_pool':PRIMES,'limits':{'seconds':25,'new_orbits':12,'regression_orbits':3,'prime_trials':45,
            'rational_RR_kernel_solves':0,'new_point_searches':0,'new_CVPs':0,'complete_atlas_runs':0},
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,TABLE,CENSUS,Path(__file__),*[REG/('orbit-%d.json'%i) for i in REGRESSION]]}})
    def expired(signum,frame):raise TimeoutError('FROZEN_25_SECOND_LIMIT')
    signal.signal(signal.SIGALRM,expired);signal.alarm(24)
    try:main()
    except BaseException as exc:save('failure.json',{'error':str(exc),'traceback':traceback.format_exc()});raise
    finally:signal.alarm(0)
