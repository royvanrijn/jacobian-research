#!/usr/bin/env python3
"""Checkpointed exact first two Frobenius moments for one fixed native twist."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_TWIST_FROBENIUS_MOMENTS_PROTOCOL.json'
CPP=HERE/'native_twist_frobenius_count.cpp'
SOURCE=r.OUT/'rank_jump_solubility_first_inputs_v1.json'
COVERS=r.OUT/'rank_jump_local_solubility_block_inputs_v1.json'
INPUT=r.OUT/'rank_jump_native_twist_frobenius_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-twist-frobenius-v1'


def capture():
    old=r.read(SOURCE);cover=next(c for c in r.read(COVERS)['covers'] if c['label']=='orbit-1795d')
    r.write_new(INPUT,{'schema':'rank-jump.native-twist-frobenius-inputs.v1',
                       'A':old['A'],'B':old['B'],'q':cover['form'],
                       'removed_square_root':cover['removed_rational_square_root'],
                       'lift':old['split_lift_maps']['orbit-1795d'],
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,SOURCE,COVERS)}})


def geometry():
    from sage.all import QQ,GF,PolynomialRing
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(inp['A']);B=R(inp['B']);q=R(inp['q'])
    delta=-16*(4*A**3+27*B**2)
    assert delta.degree()==24 and delta.is_squarefree() and q.gcd(delta)==1
    x0,x1,y0,y1=[R(inp['lift'][k+'_coefficients']) for k in ('x0','x1','y0','y1')]
    s=QQ(inp['removed_square_root']);x1*=s;y1*=s
    assert y0*y0+q*y1*y1==x0**3+3*q*x0*x1*x1+A*x0+B
    assert 2*y0*y1==3*x0*x0*x1+q*x1**3+A*x1
    # P-sigma(P); this expression is invariant in x and anti-invariant in y.
    xa=K(y0*y0)/(x1*x1*q)-2*x0
    ya=K(y0)/(x1*q)*(x0-xa)-y1
    X=R(q*xa);Y=R(q*q*ya)
    assert Y*Y==X**3+A*q*q*X+B*q**3
    assert X.degree()<=6 and Y.degree()<=9
    Fp=PolynomialRing(GF(131),'t');ap,bp,qp=map(Fp,(A,B,q));dp=Fp(delta)
    assert dp.degree()==24 and dp.is_squarefree() and qp.degree()==2 and qp.discriminant()!=0 and qp.gcd(dp)==1
    assert 4*ap[8]**3+27*bp[12]**2!=0
    return {'prime':131,'twist_label':'orbit-1795d','known_twist_section':{'X':list(map(str,X.list())),'Y':list(map(str,Y.list()))},
            'known_twist_section_height':6,'height_source':'Theorem F2: height12 anti-invariant point on the quadratic pullback; divide by base-change degree2.',
            'good_reduction':{'original_discriminant_degree':int(dp.degree()),'original_discriminant_squarefree':True,
                              'branch_discriminant':int(qp.discriminant()),'branch_discriminant_gcd_degree':int(qp.gcd(dp).degree()),'infinity_smooth':True},
            'modular_coefficients':[[int(v) for v in f.list()] for f in (ap,bp,qp)]}


def run():
    WORK.mkdir(parents=True,exist_ok=True);geo_path=WORK/'geometry.json'
    if not geo_path.exists():r.write_new(geo_path,geometry())
    geo=r.read(geo_path);exe=WORK/'count';command=['c++','-O3','-std=c++17',str(CPP),'-o',str(exe)]
    log=WORK/'count.txt'
    if not log.exists():
        compile_result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
        r.write_new(WORK/'compile.json',{'command':command,'returncode':compile_result.returncode,'output':compile_result.stdout,'source_sha256':r.digest(CPP.read_bytes())})
        assert compile_result.returncode==0
        data=' '.join(str(x) for row in geo['modular_coefficients'] for x in row)+'\n'
        with log.open('x') as out:
            try:
                proc=subprocess.run([str(exe)],input=data,text=True,stdout=out,stderr=subprocess.PIPE,timeout=60)
                execution={'returncode':proc.returncode,'stderr':proc.stderr,'status':'COMPLETE' if proc.returncode==0 else 'FAILED'}
            except subprocess.TimeoutExpired:execution={'status':'TIMEOUT'}
        r.write_new(WORK/'execution.json',execution)
    assert r.read(WORK/'execution.json')['status']=='COMPLETE'
    lines=log.read_text().splitlines();assert lines[0]=='D 2'
    ledger=[list(map(int,line.split())) for line in lines[1:]]
    for n in (1,2):
        rr=[x for x in ledger if x[0]==n]
        assert sum(x[3] for x in rr)==131**n+1
        expected={(a,b) for b in range(1 if n==1 else 66) for a in range(131)}|{(-1,0)}
        assert {(x[1],x[2]) for x in rr}==expected and len(rr)==len(expected)
    rows=[]
    for name,column,N,lower in [('original',4,20,17),('native_twist',5,24,1)]:
        traces=[-sum(x[3]*x[column] for x in ledger if x[0]==n) for n in (1,2)]
        s1=F(traces[0],131);s2=F(traces[1],131**2);m2=(N+s2)/2
        assert -N<=s1<=N and m2>=s1*s1/N and m2<=N
        c=(s1-m2)/(N-s1)
        assert c!=1
        bound=(m2-2*c*s1+N*c*c)/(1-c)**2
        fejer=(N+F(4,3)*s1+F(2,3)*s2)/3
        upper=bound.numerator//bound.denominator
        assert lower<=upper and bound<=fejer
        rows.append({'id':name,'cohomology_dimension':N,'Frobenius_traces':traces,
                     'normalized_traces':list(map(str,(s1,s2))),'cosine_first_moment':str(s1),'cosine_second_moment':str(m2),
                     'quadratic_center':str(c),'exact_eigenvalue_one_bound':str(bound),'Fejer_degree_two_bound':str(fejer),
                     'arithmetic_generic_rank_lower_bound':lower,'arithmetic_generic_rank_upper_bound':upper})
    r.write_new(OUTPUT,{'schema':'rank-jump.native-twist-frobenius.v1','status':'PASS','geometry':geo,'rows':rows,'fibre_trace_ledger':ledger,
                        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,CPP,Path(__file__),HERE/'retrospective.py')},
                        'execution':r.read(WORK/'execution.json'),'compile':r.read(WORK/'compile.json'),
                        'boundary':r.read(PROTOCOL)['boundary']})
    for row in rows:print(row['id'],row['Frobenius_traces'],'rank upper',row['arithmetic_generic_rank_upper_bound'],'exact',row['exact_eigenvalue_one_bound'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run']);mode=p.parse_args().mode
    if mode=='capture':capture()
    else:run()
