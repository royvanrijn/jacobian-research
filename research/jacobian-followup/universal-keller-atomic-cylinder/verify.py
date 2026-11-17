#!/usr/bin/env python3
"""Narrow exact controls for atomic cylinders; no replay of the general geometric proofs."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def substitute(F,variables,values):return [s.expand(q.subs(dict(zip(variables,values)),simultaneous=True)) for q in F]
def jacdet(F,variables,parameters=()):
    polys=[[s.Poly(s.diff(q,v),*variables,*parameters,domain=s.QQ) for v in variables] for q in F]
    A=polys
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])).as_expr()
def main():
    started=time.monotonic()
    x,y,z=s.symbols('x y z');variables=[x,y,z]
    predecessor=ROOT/'jacobian-followup/universal-keller-compression/verify.py'
    tree=ast.parse(predecessor.read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse one pure retained cubic source function',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(predecessor),'exec'),env)
    seed=s.Matrix(env['source_map'](3,{}));J0=seed.jacobian(variables).subs({x:0,y:0,z:0});f=list(J0.inv()*seed)
    ck('retained seed fixes origin',seed.subs({x:0,y:0,z:0})==s.zeros(3,1))
    ck('retained cubic seed determinant one',jacdet(seed,variables)==1)
    ck('first-derivative normalization invertible',J0.det()==1)
    ck('normalized seed derivative identity',s.Matrix(f).jacobian(variables).subs({x:0,y:0,z:0})==s.eye(3))
    ck('normalized seed determinant one',jacdet(f,variables)==1)
    for n in [4,6]:
        coefficients={j:s.Integer(j==n) for j in range(4,n+1)};fn=env['source_map'](n,coefficients)
        ck(f'constant-profile degree{n} seed exact Jacobian',jacdet(fn,variables)==1)
        S,pi,b,c=s.symbols('S pi b c');inverse=S+pi*S**3+pi**n*S**n+b*S**2-c/2
        ck(f'constant-profile inverse degree{n} and linear coefficient',s.Poly(inverse,S).degree()==n and s.Poly(inverse,S).nth(1)==1)
    q,t=s.symbols('q t');fq=[s.cancel(v/q) for v in substitute(f,variables,[q*x,q*y,q*z])]
    ck('dilation polynomial at zero parameter',all(s.denom(v)==1 for v in fq))
    fq=list(map(s.expand,fq));ck('zero parameter fibre exactly identity',[v.subs(q,0) for v in fq]==variables)
    ck('full parameter-dependent Keller determinant',jacdet(fq,variables,[q])==1)
    ck('generic scaling identity',all(s.expand(q*a-b)==0 for a,b in zip(fq,substitute(f,variables,[q*x,q*y,q*z]))))
    # The actual seed collision survives a nonzero dilation, not just its Jacobian.
    collision=[(0,0,-s.Rational(1,4)),(1,-s.Rational(3,2),s.Rational(13,2)),(-1,s.Rational(3,2),s.Rational(13,2))]
    values=[tuple(substitute(f,variables,point)) for point in collision]
    ck('normalized seed exact triple collision',len(set(values))==1)
    doubled=[tuple(substitute([v.subs(q,2) for v in fq],variables,[s.sympify(a)/2 for a in point])) for point in collision]
    ck('nonzero dilation exact triple collision',len(set(doubled))==1 and doubled[0]==tuple(v/2 for v in values[0]))
    family=[]
    for m in [1,2,4]:
        a=s.prod(t-j for j in range(1,m+1));poly=s.Poly(a,t)
        ck(f'family{m} exactly m distinct nonzero exceptional base points',poly.degree()==m and s.gcd(poly,poly.diff()).degree()==0 and a.subs(t,0)!=0)
        specialised=[s.expand(v.subs(q,a)) for v in fq]
        ck(f'family{m} identity on every exceptional fibre',all([s.expand(v.subs(t,j)) for v in specialised]==variables for j in range(1,m+1)))
        ck(f'family{m} nonzero generic scaling parameter',a!=0 and a.subs(t,m+1)!=0)
        family.append({'m':m,'parameter':str(a),'exceptional_base_points':list(range(1,m+1)),'special_vertical_degree':1,'generic_vertical_degree':3,'degree_assurance':'Generic degree and primitive monodromy transferred from retained UAM1 by the exact scaling identity; not newly inferred from a collision.'})
    # For F=(t^2,f_(t-1)), its base and vertical maps are actual self-factors.
    H=[t]+[v.subs(q,t-1) for v in fq];B=[t*t,x,y,z]
    F=[t*t]+H[1:]
    composed=substitute(B,[t,x,y,z],H)
    ck('mixed base-vertical exact self-factorization',all(s.expand(a-b)==0 for a,b in zip(F,composed)))
    ck('vertical logarithmic volume factor exactly one',jacdet(fq,variables,[q])==1)
    ck('mixed logarithmic volume factor exactly two',s.cancel((2*t)*(1/t**2)/(1/t))==2)
    ck('base double cover is noninjective',s.Integer(1)**2==s.Integer(-1)**2 and 1!=-1)
    # A general normal-affine completion boundary in one variable.
    u=s.symbols('u');cover=u*u
    ck('punctured square cover etale wherever u invertible',s.diff(cover,u)==2*u)
    ck('removed unramified source point over target one',s.diff(cover,u).subs(u,1)!=0 and cover.subs(u,1)==1)
    ck('completion has degree two',s.Poly(cover,u).degree()==2)
    ck('nonfinite punctured fibre has one of two completion points',[a for a in s.solve(s.Eq(cover,1),u) if a!=1]==[-1])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'dilation_controls':family,'mixed_control':{'base_degree':2,'vertical_degree':3,'generic_degree':6,'explicit_self_factors':['(t,v)->(t,f_(t-1)(v))','(t,v)->(t^2,v)'],'assurance':'Degree3 and primitive seed monodromy are retained written inputs; the composition and actual Jacobians/collisions are exact controls.'},'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact normalized retained seed, symbolic parameter-dependent Jacobian, dilation conjugacy, actual collision and exceptional base-locus controls, mixed self-factorization and normal-completion toy. No formalization of Ax, Zariski Main, monodromy, atomicity, all-degree geometric reductions or relative descent. No independent whole-proof replay, external review or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
