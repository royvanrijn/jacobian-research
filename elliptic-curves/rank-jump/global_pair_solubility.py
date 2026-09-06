#!/usr/bin/env python3
"""Complete local checks and bounded 2-descent on two frozen auxiliary carriers."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import traceback
import retrospective as r
from cover_experiment import sqrtq
from local_solubility_blocks import qp_square, value

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'GLOBAL_PAIR_SOLUBILITY_PROTOCOL.json'
OLD_INPUT=r.OUT/'rank_jump_local_solubility_block_inputs_v1.json'
OLD_RESULT=r.OUT/'rank_jump_local_solubility_blocks_v1.json'
INPUT=r.OUT/'rank_jump_global_pair_solubility_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_global_pair_solubility_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-global-pair-solubility-v1'


def capture():
    old=r.read(OLD_INPUT);results=r.read(OLD_RESULT);covers=old['covers']
    def height(i):return max(abs(x) for x in covers[i]['form'])
    def key(pair):
        i,j=pair['indices'];return max(height(i),height(j)),min(height(i),height(j)),covers[i]['label'],covers[j]['label']
    pair=min((p for p in results['pairs'] if not p['observed_together']),key=key)
    assert pair['indices']==[8,10]
    common=10;group=next(g for g in old['observed_groups'] if g['mask']>>common&1)
    partner=min((i for i in range(14) if i!=common and group['mask']>>i&1),key=lambda i:(height(i),covers[i]['label']))
    assert partner==4
    t=F(group['published_parameter']);f=covers[common]['form']
    u=sqrtq(F(value(f,t.numerator,t.denominator),t.denominator**2));assert u is not None
    cases=[]
    for name,j in [('cross_group',8),('observed_positive',partner)]:
        mask=(1<<common)|(1<<j)
        witnesses=[]
        for p in results['places']:
            w=next(w for w in p['maximal_proved_masks'] if w['mask']&mask==mask)
            witnesses.append({'prime':p['prime'],'base_point':w['rational_base_point']})
        cases.append({'id':name,'labels':[covers[common]['label'],covers[j]['label']],
                      'forms':[f,covers[j]['form']],'single_conic_anchor':[str(t),str(u)],
                      'old_local_witnesses':witnesses})
    r.write_new(INPUT,{'schema':'rank-jump.global-pair-solubility-inputs.v1','cases':cases,
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,OLD_INPUT,OLD_RESULT)},
                       'boundary':'Single-conic anchor is a retained cover lift, not a new specialization or exceptional-point input. Cross-group carrier has no observed common split in the 32-fibre cohort.'})


def worker(case):
    from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,Jacobian,pari
    from sage.version import version
    inp=r.read(INPUT);c=inp['cases'][case];name=c['id'];protocol=r.read(PROTOCOL)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    ring=PolynomialRing(QQ,'v');v=ring.gen();f,g=[ring(q) for q in c['forms']]
    t,u=map(QQ,c['single_conic_anchor']);assert u*u==f(t) and u
    a,b=f[2],f[1]
    D=v*v-a;N=t*v*v-2*u*v+a*t+b
    U=-u*v*v+f.derivative()(t)*v-u*a
    assert U*U==f[0]*D*D+f[1]*N*D+f[2]*N*N
    quartic=g[0]*D*D+g[1]*N*D+g[2]*N*N
    den=ZZ(quartic.denominator());integral=ring(den*den*quartic)
    assert integral.degree()==4 and integral.discriminant()
    e,d,c2,b2,a2=[integral[k] for k in range(5)]
    I=12*a2*e-3*b2*d+c2*c2
    J=72*a2*c2*e+9*b2*c2*d-27*a2*d*d-27*b2*b2*e-2*c2**3
    jac=EllipticCurve(QQ,[-27*I,-27*J])
    mult=PolynomialRing(QQ,names=('v','w'));vv,ww=mult.gens()
    other=Jacobian(ww*ww-sum(integral[k]*vv**k for k in range(5)))
    assert jac.is_isomorphic(other)
    E=jac.minimal_model();assert jac.is_isomorphic(E)
    bad=abs(ZZ(2*f.discriminant()*g.discriminant()*f.resultant(g)))
    factorization=list(bad.factor(proof=True))
    assert all(p.is_prime(proof=True) for p,e in factorization)
    geometry={'id':name,'forms':c['forms'],'single_conic_anchor':c['single_conic_anchor'],
              'parameter_numerator':list(map(str,N.list())),'parameter_denominator':list(map(str,D.list())),
              'conic_root_numerator':list(map(str,U.list())),'quartic_coefficients':list(map(str,integral.list())),
              'quartic_square_scaling':str(den),'invariants_I_J':[str(I),str(J)],
              'Jacobian_model':list(map(str,jac.a_invariants())),
              'Sage_Jacobian_model':list(map(str,other.a_invariants())),
              'minimal_Jacobian_model':list(map(str,E.a_invariants())),
              'bad_reduction_integer':str(bad),'bad_support_factorization':[[str(p),int(e)] for p,e in factorization],
              'factor_primality_proved':True,'genus':1,'torsor_period_divides':2,
              'software':{'sage':version,'pari':str(pari('version()'))}}
    r.write_new(WORK/(name+'_geometry.json'),geometry)
    print(name,'geometry complete; bad primes',[str(p) for p,e in factorization],flush=True)
    local=[]
    def works(T,Z,p):return all(qp_square(value(q,T,Z),p) for q in c['forms'])
    for p,e in factorization:
        p=int(p)
        candidates=[tuple(w['base_point']) for w in c['old_local_witnesses'] if w['prime']==p]
        candidates += [(int(t.numerator()),int(t.denominator()))]
        candidates += [(k,1) for k in [0]+[x for i in range(1,129) for x in (-i,i)]]
        candidates += [(1,p),(1,-p),(1,0)]
        witness=next(((T,Z) for T,Z in candidates if works(T,Z,p)),None)
        local.append({'prime':str(p),'status':'PROVED' if witness else 'UNKNOWN',
                      'base_point':list(witness) if witness else None,
                      'exact_form_values':[str(value(q,*witness)) for q in c['forms']] if witness else None})
    real=next(((k,1) for k in [0]+[x for i in range(1,129) for x in (-i,i)] if all(value(q,k,1)>0 for q in c['forms'])),None)
    lc={'id':name,'finite_bad_places':local,'real_witness':list(real) if real else None,
        'everywhere_locally_soluble':bool(real and all(x['status']=='PROVED' for x in local)),
        'other_places':'Smooth projective genus-one reduction outside the factored support: Weil bound supplies an Fp point and Hensel lifts it.'}
    r.write_new(WORK/(name+'_local.json'),lc);print(name,'local complete',lc['everywhere_locally_soluble'],flush=True)
    pari.allocatemem(protocol['limits']['pari_stack_bytes'],silent=True)
    ans=pari.ellrank(pari.ellinit(E.a_invariants()),0)
    nroots=len(E.division_polynomial(2).roots(QQ,multiplicities=False))
    assert nroots in (0,1,3);T={0:0,1:1,3:2}[nroots]
    lo,hi,ct=map(int,ans[:3]);pts=[[str(x) for x in P] for P in ans[3]]
    assert all(E([QQ(x) for x in P]) for P in pts)
    sha=ct if lo==hi else 'UNKNOWN'
    descent={'id':name,'raw_ellrank':str(ans),'rank_lower_bound':lo,'rank_upper_bound':hi,
             'CT_Sha2_mod_2Sha4_dimension':ct,'rational_2_torsion_dimension':T,
             'full_2_Selmer_dimension':hi+T+ct,'Sha_2_dimension':sha,'points':pts,
             'global_carrier_solubility_proved':bool(lc['everywhere_locally_soluble'] and sha==0),
             'reason':'An everywhere locally soluble period-dividing-two torsor is trivial when Sha(J)[2]=0.'}
    r.write_new(WORK/(name+'_descent.json'),descent)
    print(name,'descent',lo,hi,'Sha2',sha,'global',descent['global_carrier_solubility_proved'],flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);results=[]
    for k,c in enumerate(r.read(INPUT)['cases']):
        name=c['id'];log=WORK/(name+'.log')
        if not log.exists():
            with log.open('x') as f:
                try:
                    p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--case',str(k)],cwd=r.ROOT,stdout=f,stderr=f,timeout=60)
                    status='COMPLETE' if p.returncode==0 else 'FAILED'
                except subprocess.TimeoutExpired:status='TIMEOUT'
            r.write_new(WORK/(name+'_execution.json'),{'status':status})
        row={'id':name,'execution':r.read(WORK/(name+'_execution.json')),
             'log':log.read_text()}
        for stage in ('geometry','local','descent'):
            path=WORK/(name+'_'+stage+'.json')
            row[stage]=r.read(path) if path.exists() else {'status':'UNKNOWN'}
        results.append(row);print(name,row['execution'],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.global-pair-solubility.v1',
                        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py',HERE/'local_solubility_blocks.py')},
                        'rows':results,'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        try:worker(a.case)
        except Exception:traceback.print_exc();raise
