#!/usr/bin/env python3
"""Selected exact controls for written response and rank-one feedback theorems."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
checks=[]
def ck(label,ok):
    assert ok,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    started=time.monotonic(); x,y,z=s.symbols('x y z'); xyz=[x,y,z]
    source=ROOT/'jacobian/universal-keller-compression/verify.py'
    nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map']
    ck('reuse only pure retained compiler',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),env);compiler=env['source_map']
    def sub(F,vs,vals):return [f.subs(dict(zip(vs,vals)),simultaneous=True) for f in F]
    K=s.Matrix(compiler(3,{})); DK=K.jacobian(xyz); ck('actual cubic Jacobian one',zero(DK.det()-1))
    pi,w,eps=s.symbols('pi w eps'); u,v,H0,H1,H2,H3=s.symbols('u v H0 H1 H2 H3')
    cases=[({4:u,5:v},{4:H0+H1*pi,5:H2+H3*pi*pi}),({4:u*pi**3,5:v*pi+u*u},{4:v*pi**3+u,5:1+u*pi*pi}),({4:pi*u,5:s.Integer(1)},{4:pi,5:s.Integer(-1)})]
    for index,(g,h) in enumerate(cases,1):
        coeff={j:g[j].subs(pi,K[0]) for j in (4,5)}; F=s.Matrix(compiler(5,coeff)); A=F.jacobian(xyz).subs(y,0)
        perturbed=s.Matrix(compiler(5,{j:coeff[j]+eps*h[j].subs(pi,K[0]) for j in (4,5)})); Vh=perturbed.diff(eps).subs({eps:0,y:0})
        response=s.expand((A.adjugate()*Vh)[0]); B=sum(h[j].subs(pi,z)*(x*z)**j for j in (4,5)); C0=sum(j*(j-2)*g[j].subs(pi,z)*(x*z)**j for j in (4,5)); expect=sum((j+1)*h[j].subs(pi,z)*(x*z)**j for j in (4,5))-3*x*(x*z)*B-C0*B/x
        ck(f'arbitrary direction {index} actual sliced Jacobian one',zero(A.det()-1))
        ck(f'arbitrary direction {index} exact response identity',zero(response-expect))
        ck(f'arbitrary direction {index} nonzero independent response',response!=0 and not s.Poly(s.expand(response.subs(z,w/x)*x**20),x).is_zero)
    # All polynomial degree directions are distinct, including equal w-powers.
    pairs=[(j,e) for j in range(4,13) for e in range(8)]
    ck('monomial encoding keeps all sampled coefficient directions distinct',len({(-e,j+e) for j,e in pairs})==len(pairs))
    S=x/(1+x*y); Eprime=1+2*K[1]*S+3*K[0]*S*S
    ck('actual inverse x reconstruction exact',zero(S/Eprime-x))
    ck('actual inverse y reconstruction exact',zero((1-Eprime)/S-y))
    a,b,c=s.symbols('a b c'); r1,r2,tau,T=s.symbols('r1 r2 tau T'); m1,m2,m3=s.symbols('m1 m2 m3')
    Q=T+b*T*T+a*T**3-c/2
    mu=s.Matrix([m1,m2,m3]); Delta=DK.adjugate()*mu
    ck('generic output direction lifted exactly',all(zero(q) for q in DK*Delta-mu))
    R=[m3*K[0]-m1*K[2],m3*K[1]-m2*K[2]]
    ck('generic chart two output integrals exact',all(zero(q) for q in s.Matrix(R).jacobian(xyz)*Delta))
    expr=Q.subs({a:r1/m3+tau*m1,b:r2/m3+tau*m2,c:tau*m3})
    ck('generic chart inverse equation exactly linear in time',zero(expr-(T+r2*T*T/m3+r1*T**3/m3)-tau*(m1*T**3+m2*T*T-m3/2)))
    directions=[(1,0,0),(0,1,0),(0,0,1),(2,-1,3),(s.Rational(1,3),-1,s.Rational(2,3))]
    for direction in directions:
        mu0=s.Matrix(direction); delta=DK.adjugate()*mu0
        pivot=next(i for i,m in enumerate(direction) if m!=0); forms=[]
        for j in range(3):
            if j!=pivot:forms.append(direction[pivot]*K[j]-direction[j]*K[pivot])
        ck('direction '+str(direction)+' polynomial output kernel generators',all(zero(q) for q in s.Matrix(forms).jacobian(xyz)*delta))
        ck('direction '+str(direction)+' nonzero inverse line coefficient',any(direction) and s.Poly(direction[0]*T**3+direction[1]*T*T-direction[2]/2,T).as_expr()!=0)
    # Actual N=2 source-coordinate controls, no assumed common scalar presentation.
    P=[s.symbols('p0:3'),s.symbols('q0:3')]; U=s.symbols('U0:3'); coords=list(P[0])+list(P[1])+list(U)
    KU=s.Matrix(sub(K,xyz,U)); mu0=s.Matrix([s.Rational(1,3),-1,s.Rational(2,3)]); R1=KU[1]+3*KU[0]; R2=KU[2]-2*KU[0]
    C=sum((3*p[0]*p[0]-p[0])/2 for p in P); beta=1+C+C*C; Vout=KU-beta*mu0
    recovered=[Vout[1]+3*Vout[0],Vout[2]-2*Vout[0]]
    ck('both independent cubic output controls recovered exactly',zero(recovered[0]-R1) and zero(recovered[1]-R2))
    ck('two recovered controls algebraically independent by nonzero Jacobian minor',s.Matrix([R1,R2]).jacobian(U).rank()==2)
    def make_map(allowed):
        out=[]
        for block in P:
            sig=sub(K,xyz,block)[0]
            g={4:R1+sig*R2,5:R2+sig*R1*R1} if allowed else {4:U[0],5:U[1]}
            out+=sub(compiler(5,g),xyz,block)
        return s.Matrix(out+list(Vout))
    allowed=make_map(True); bad=make_map(False)
    origin=dict(zip(coords,[0]*9)); point=dict(zip(coords,[1,0,1,0,0,0,0,0,0])); other=dict(zip(coords,[0,1,0,1,0,1,1,0,0]))
    bad_det=bad.jacobian(coords).subs(point).det()
    ck('actual multiple-source-control feedback not Keller',bad_det!=1)
    ck('actual multiple-source-control feedback determinant equals one at zero source',bad.jacobian(coords).subs(origin).det()==1)
    allowed_J=allowed.jacobian(coords)
    for index,q in enumerate([origin,point,other]):ck('actual recovered multiple-control feedback Jacobian control '+str(index),allowed_J.subs(q).det()==1)
    # Verify recovery of the effective profiles before compiler substitution; this is the exact QH identity.
    sig=s.symbols('sig'); gp=[R1+sig*R2,R2+sig*R1*R1]; gh=[recovered[0]+sig*recovered[1],recovered[1]+sig*recovered[0]**2]
    ck('both nonlinear signal-dependent coefficient profiles recover exactly',all(zero(f-g) for f,g in zip(gp,gh)))
    # Another selector lies outside b(C): it is Keller despite non-invariant coefficients, but factors forward.
    sig1=sub(K,xyz,P[0])[0]; signal_map=s.Matrix(list(bad[:6])+list(KU-sig1*mu0))
    ck('different signal-only selector allows actual non-invariant controls',signal_map.jacobian(coords).subs(point).det()==1)
    ck('different signal-only selector second actual Jacobian control',signal_map.jacobian(coords).subs(other).det()==1)
    ck('source controls not annihilated by cubic translation',list((DK.adjugate()*mu0).subs({x:0,y:0,z:0})[:2])==[s.Rational(1,3),s.Integer(2)])
    Hforward=list(bad[:6])+list(U); Qforward=list(P[0])+list(P[1])+list(KU-P[0][0]*mu0)
    ck('different selector exact whole-map forward factorization',all(zero(f-g) for f,g in zip(sub(Qforward,coords,Hforward),signal_map)))
    # A different auxiliary automorphism obeys the same output criterion.
    GA=s.Matrix([U[0]+U[1]**2,U[1]+U[2]**3,U[2]]); DA=GA.jacobian(U); ma=s.Matrix([1,1,0]); da=DA.adjugate()*ma; RA=s.Matrix([GA[0]-GA[1],GA[2]])
    ck('different auxiliary actual Jacobian one',DA.det()==1)
    ck('different auxiliary output invariants exact',all(zero(f) for f in RA.jacobian(U)*da))
    ck('different auxiliary time coordinate has unit derivative',zero(s.Matrix([GA[1]]).jacobian(U).dot(da)-1))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'arbitrary_direction_cases':len(cases),'constant_output_directions':[list(map(str,m)) for m in directions],'actual_bad_multiple_control_jacobian':str(bad_det),'actual_good_multiple_control_jacobian_points':3,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Selected exact arbitrary coefficient responses, inverse reconstruction, generic-direction cubic identities, actual nine-variable failed and allowed feedback maps, effective-profile recovery and another auxiliary map, and a different-selector forward-factor control. Not universal proof replay, formal verification of injectivity/geometric criteria, independent review or publication-priority evidence.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],'bad_jacobian':str(bad_det)}))
if __name__=='__main__':main()
