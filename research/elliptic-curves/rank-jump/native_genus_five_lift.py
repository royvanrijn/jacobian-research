#!/usr/bin/env python3
"""The fixed native genus-five cover, its branch algebra and elliptic factors."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_GENUS_FIVE_LIFT_PROTOCOL.json'
FORMS=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
PAIR=r.OUT/'rank_jump_minimal_native_block_carrier_inputs_v1.json'
INPUT=r.OUT/'rank_jump_native_genus_five_lift_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_native_genus_five_lift_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-genus-five-lift-v1'


def capture():
    case=next(c for c in r.read(FORMS)['cases'] if c['id']=='08234-003')
    covers=[next(c for c in case['covers'] if c['label']==label) for label in ('orbit-01333','orbit-0b2d0','orbit-19e45')]
    r.write_new(INPUT,{'schema':'rank-jump.native-genus-five-lift-inputs.v1','covers':covers,
        'retained_lift':{'t':'-288/65','roots':['848451138/65','44253/5','3924473/65']},
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,FORMS,PAIR)}})


def geometry():
    from sage.all import QQ,PolynomialRing,matrix,vector
    inp=r.read(INPUT);R=PolynomialRing(QQ,'t');t=R.gen();qs=[R(c['form']) for c in inp['covers']]
    assert all(q.degree()==2 and q.is_squarefree() for q in qs)
    assert all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
    t0=QQ(inp['retained_lift']['t']);roots=list(map(QQ,inp['retained_lift']['roots']))
    assert all(u*u==q(t0) and u for u,q in zip(roots,qs,strict=True))
    h,f,g=qs
    # Multiplication in Q[t,u,v]/(h(t),u^2-f(t),v^2-g(t)).
    # Basis t^i*u^j*v^k with each exponent 0 or 1.
    basis=[(i,j,k) for k in range(2) for j in range(2) for i in range(2)]
    def mul(a,b):
        out=[QQ(0)]*8
        for ia,ca in enumerate(a):
            for ib,cb in enumerate(b):
                if not ca or not cb:continue
                i=basis[ia][0]+basis[ib][0];j=basis[ia][1]+basis[ib][1];k=basis[ia][2]+basis[ib][2]
                pol=(t**i*f**(j//2)*g**(k//2))%h
                for e in range(2):out[basis.index((e,j%2,k%2))]+=ca*cb*pol[e]
        return vector(QQ,out)
    one=vector(QQ,[1]+[0]*7);element=vector(QQ,[0,1,1,0,1,0,0,0])
    cols=[mul(element,vector(QQ,[int(i==j) for i in range(8)])) for j in range(8)]
    M=matrix(QQ,cols).transpose();cp=M.charpoly();powers=[one]
    for _ in range(1,8):powers.append(mul(powers[-1],element))
    assert matrix(QQ,powers).rank()==8 and cp.is_squarefree()
    return {'pair_carrier_degree':4,'pair_carrier_genus':1,'triple_carrier_degree':8,'triple_carrier_genus':5,
        'relative_cover_degree':2,'geometric_ramification_points':8,'Prym_dimension':4,
        'quotient_characters':[{'labels':[inp['covers'][i]['label'] for i in indices],
            'polynomial':list(map(str,pol.list())),'genus':int((pol.degree()-2)//2)}
            for indices,pol in [([0,1],h*f),([0,2],h*g),([1,2],f*g),([0,1,2],h*f*g)]],
        'ramification_basis':[list(x) for x in basis],
        'ramification_primitive_element':'t+u+v; u^2=f(t), v^2=g(t), h(t)=0',
        'ramification_multiplication_matrix':[[str(x) for x in row] for row in M],
        'ramification_primitive_polynomial':list(map(str,cp.list())),
        'ramification_primitive_power_matrix':[[str(x) for x in row] for row in powers],
        'retained_lift':inp['retained_lift']}


def worker(i):
    from sage.all import QQ,PolynomialRing,EllipticCurve,Jacobian,pari
    from sage.env import SAGE_VERSION
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    geo=geometry();case=geo['quotient_characters'][i];R=PolynomialRing(QQ,'t');q=R(case['polynomial'])
    e,d,c,b,a=q.list();I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    E=EllipticCurve(QQ,[-27*I,-27*J]).minimal_model()
    RR=PolynomialRing(QQ,names=('t','y'));t,y=RR.gens()
    check=Jacobian(y*y-sum(q[i]*t**i for i in range(5)));assert check.is_isomorphic(E)
    row={'index':i,**case,'Jacobian_model':list(map(str,E.a_invariants())),
         'invariants_I_J':[str(I),str(J)],'software':{'sage':SAGE_VERSION,'pari':str(pari.version())}}
    r.write_new(WORK/f'{i}_geometry.json',row);print('Geometry',i,'complete',flush=True)
    pari.allocatemem(r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    ans=pari.ellrank(pari.ellinit(E.a_invariants()),0);lo,hi,ct=map(int,ans[:3]);pts=[[str(x) for x in P] for P in ans[3]]
    assert all(E([QQ(x) for x in P]) for P in pts)
    torsion={0:0,1:1,3:2}[len(E.division_polynomial(2).roots(QQ,multiplicities=False))]
    row.update({'rank_lower_bound':lo,'rank_upper_bound':hi,'raw_ellrank':str(ans),
                'CT_Sha2_mod_2Sha4_dimension':ct,'rational_2_torsion_dimension':torsion,
                'full_2_Selmer_dimension':hi+torsion+ct,'points':pts})
    r.write_new(WORK/f'{i}_descent.json',row);print(i,'rank',lo,hi,'CT',ct,flush=True)


def run():
    WORK.mkdir(parents=True,exist_ok=True);geopath=WORK/'cover_geometry.json'
    if not geopath.exists():r.write_new(geopath,geometry())
    rows=[]
    for i in range(3):
        log=WORK/f'{i}.log';execution=WORK/f'{i}_execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i)],stdout=out,stderr=out,timeout=60)
                    status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
            r.write_new(execution,status)
        row={'index':i,'execution':r.read(execution),'log':log.read_text()}
        for key in ('geometry','descent'):
            path=WORK/f'{i}_{key}.json';row[key]=r.read(path) if path.exists() else {'status':'UNKNOWN'}
        rows.append(row);print(row['execution'],row['log'],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.native-genus-five-lift.v1','geometry':r.read(geopath),'elliptic_factors':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:worker(a.case)
