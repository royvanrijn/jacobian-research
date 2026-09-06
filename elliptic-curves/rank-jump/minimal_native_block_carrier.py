#!/usr/bin/env python3
"""A single retrospective genus-one carrier and its product-square quotient."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'MINIMAL_NATIVE_BLOCK_CARRIER_PROTOCOL.json'
FORMS=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
TRIPLE=r.OUT/'rank_jump_native_triple_intersection_v1.json'
INPUT=r.OUT/'rank_jump_minimal_native_block_carrier_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_minimal_native_block_carrier_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-minimal-native-block-carrier-v1'


def capture():
    case=next(c for c in r.read(FORMS)['cases'] if c['id']=='08234-003')
    covers=[next(c for c in case['covers'] if c['label']==label) for label in ('orbit-0b2d0','orbit-19e45')]
    r.write_new(INPUT,{'schema':'rank-jump.minimal-native-block-carrier-inputs.v1','covers':covers,
        'anchor':{'t':'-288/65','u':'44253/5','v':'3924473/65'},
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,FORMS,TRIPLE)}})


def worker():
    from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,Jacobian,pari
    from sage.env import SAGE_VERSION
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'z');z=R.gen();f,g=[R(c['form']) for c in inp['covers']]
    t,u,v=[QQ(inp['anchor'][k]) for k in ('t','u','v')]
    assert f.is_squarefree() and g.is_squarefree() and f.gcd(g)==1 and u*u==f(t) and v*v==g(t)
    # A line of slope z through the specified rational point of the first conic.
    D=z*z-f[2];N=t*z*z-2*u*z+f[2]*t+f[1];U=-u*z*z+f.derivative()(t)*z-u*f[2]
    assert U*U==sum(f[i]*N**i*D**(2-i) for i in range(3))
    quartic=sum(g[i]*N**i*D**(2-i) for i in range(3))
    den=ZZ(quartic.denominator());q=R(den**2*quartic);assert q.degree()==4 and q.is_squarefree()
    e,d,c,b,a=q.list();I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    raw=EllipticCurve(QQ,[-27*I,-27*J]);E=raw.minimal_model()
    RR=PolynomialRing(QQ,names=('z','w'));zz,ww=RR.gens()
    independent=Jacobian(ww*ww-sum(q[i]*zz**i for i in range(5)));assert independent.is_isomorphic(E)
    geometry={'forms':[x['form'] for x in inp['covers']],'anchor':inp['anchor'],
        'carrier_fibre_degree':4,'carrier_genus':1,'product_quotient_fibre_degree':2,'product_quotient_genus':1,
        'conic_parameter_numerator':list(map(str,N.list())),'conic_parameter_denominator':list(map(str,D.list())),
        'conic_root_numerator':list(map(str,U.list())),'quartic_coefficients':list(map(str,q.list())),
        'quartic_square_scale':str(den),'Jacobian_invariants_I_J':[str(I),str(J)],
        'Jacobian_model':list(map(str,raw.a_invariants())),'minimal_Jacobian_model':list(map(str,E.a_invariants())),
        'software':{'sage':SAGE_VERSION,'pari':str(pari.version())}}
    r.write_new(WORK/'geometry.json',geometry);print('Carrier geometry complete',flush=True)
    pari.allocatemem(r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    ans=pari.ellrank(pari.ellinit(E.a_invariants()),0)
    lo,hi,ct=map(int,ans[:3]);pts=[[str(x) for x in P] for P in ans[3]]
    assert all(E([QQ(x) for x in P]) for P in pts)
    nroots=len(E.division_polynomial(2).roots(QQ,multiplicities=False));torsion={0:0,1:1,3:2}[nroots]
    assert torsion>=1
    descent={'rank_lower_bound':lo,'rank_upper_bound':hi,'raw_ellrank':str(ans),
        'CT_Sha2_mod_2Sha4_dimension':ct,'rational_2_torsion_dimension':torsion,
        'full_2_Selmer_dimension':hi+torsion+ct,'points':pts,
        'boundary':'Auxiliary Jacobian only; original specialized rank and point independence are separate.'}
    r.write_new(WORK/'descent.json',descent);print('Auxiliary rank',lo,hi,flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
        r.write_new(execution,status)
    result={'schema':'rank-jump.minimal-native-block-carrier.v1','execution':r.read(execution),'log':log.read_text(),
        'geometry':r.read(WORK/'geometry.json') if (WORK/'geometry.json').exists() else {'status':'UNKNOWN'},
        'descent':r.read(WORK/'descent.json') if (WORK/'descent.json').exists() else {'status':'UNKNOWN'},
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}
    r.write_new(OUTPUT,result);print(result['execution']);print(result['log'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:worker()
