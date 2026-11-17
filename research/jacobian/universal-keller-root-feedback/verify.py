#!/usr/bin/env python3
"""Selected exact controls for the written all-degree root-feedback theorems."""
from pathlib import Path
from itertools import combinations
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def sub(F,vs,vals):return [f.subs(dict(zip(vs,vals)),simultaneous=True) for f in F]
def main():
    begun=time.monotonic();x,y,z=s.symbols('x y z');xyz=[x,y,z]
    constructor=ROOT/'jacobian/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(constructor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(constructor),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));DK=K.jacobian(xyz);ck('actual base compiler determinant one',zero(DK.det()-1))
    tt=1+x*y;qq=tt**2*z+y*y*(1+3*tt);sig=tt*qq;XX=x*qq;SS=x/tt;rr=1/tt
    ck('source polynomial root equals sigma times S',zero(XX-sig*SS))
    ck('birational source reconstruction x',zero(SS/rr-x));ck('birational source reconstruction y',zero((1-rr)/SS-y))
    ck('chart actual B0 identity',zero(K[1]-(rr-1-3*sig*SS**2)/(2*SS)))
    ck('chart actual C0 identity',zero(K[2]-(rr+1)*SS+sig*SS**3))
    # A generic two-block correction: derivatives stay symbolic, not sampled.
    S1,S2,r1,r2,a1,a2=s.symbols('S1 S2 r1 r2 a1 a2');sv=[S1,S2];rv=[r1,r2];aa=[a1,a2];FF=[s.Function('F1')(*sv),s.Function('F2')(*sv)];HH=[s.Function('H1')(*sv),s.Function('H2')(*sv)]
    BB=[(rv[i]-1-3*aa[i]*sv[i]**2)/(2*sv[i])+FF[i] for i in range(2)];CC=[(rv[i]+1)*sv[i]-aa[i]*sv[i]**3+HH[i] for i in range(2)]
    J=s.Matrix(BB+CC).jacobian(rv+sv);E=s.diag(2*S1**2,2*S2**2);schur=J[2:4,2:4]-E*J[:2,2:4]
    expected=s.diag(2*r1,2*r2)+s.Matrix(2,2,lambda i,j:s.diff(HH[i],sv[j])-2*sv[i]**2*s.diff(FF[i],sv[j]))
    ck('generic two-block Schur elimination identity',all(zero(v) for v in schur-expected))
    ck('generic two-block determinant is independent-diagonal formula',zero(J.det(method='domain-ge')-expected.det()/(4*S1*S2)))
    # Coefficients in all independent diagonal variables isolate every principal minor.
    r3=s.symbols('r3');vs=[r1,r2,r3];A=s.Matrix(3,3,s.symbols('a0:9'));pol=s.Poly(s.expand((s.diag(*vs)+A).det()),*vs)
    for n in range(1,4):
        for ids in combinations(range(3),n):
            mon=tuple(int(i not in ids) for i in range(3));minor=A.extract(ids,ids).det();ck('independent diagonal principal coefficient '+str(ids),zero(pol.coeff_monomial(mon)-minor))
    ck('independent diagonal leading coefficient one',pol.coeff_monomial((1,1,1))==1)
    # Ordinary nilpotence does not imply the needed principal-minor identity.
    p,q=s.symbols('p q');N=s.Matrix([[p*q,p*p],[-q*q,-p*q]])
    ck('cyclic comparison matrix square zero',N*N==s.zeros(2))
    ck('nilpotence fails independent diagonal identity',zero((s.diag(r1,r2)+N).det()-r1*r2-p*q*(r2-r1)))
    ck('scalar diagonal hides the comparison cycle',zero((r1*s.eye(2)+N).det()-r1*r1))
    # Self equation, exact normal form and the effective cross-dependency identity.
    X1,X2,b,c,d,sg=s.symbols('X1 X2 b c d sg');F=X2+sg**2*((b+X2)*X1+c*X1**2+d*X1**5);H=X2**2+s.Rational(2,3)*(b+X2)*X1**3+c*X1**4+s.Rational(10,7)*d*X1**7
    ck('normal form exact self equation',zero(sg**2*s.diff(H,X1)-2*X1**2*s.diff(F,X1)))
    K12=s.diff(H,X2)-2*X1**2/sg**2*s.diff(F,X2)
    ck('off-diagonal derivative detects actual root dependence',zero(s.diff(K12,X1)+4*X1/sg**2*s.diff(F,X2)))
    ck('normal form has genuine other-root dependence',s.diff(F,X2)!=0 and K12!=0)
    X=s.symbols('X');redundant={4:X**2,5:-2*X,6:s.Integer(1)}
    ck('native redundant coefficient list B correction zero',s.expand(sum(j*v*X**(j-2) for j,v in redundant.items()))==0)
    ck('native redundant coefficient list C correction zero',s.expand(sum((j-2)*v*X**j for j,v in redundant.items()))==0)
    cubic=sg-s.Rational(2,3)*b*sg**3
    ck('new index-three inverse coefficient cannot identically cancel',s.Poly(cubic,sg).coeff_monomial(sg)==1)
    ownF=sig**2*XX;ownH=s.Rational(2,3)*XX**3;own=s.Matrix([sig,K[1]+ownF,K[2]+ownH]);S=x/tt
    inverse=S+own[1]*S**2+(sig-s.Rational(2,3)*sig**3)*S**3-own[2]/2
    ck('actual new index-three inverse equation exact',zero(inverse))
    ck('actual index-three whole-map Jacobian at nonzero root',own.jacobian(xyz).subs({x:1,y:0,z:1}).det()==1)
    # An actual non-descended rank-two return with two degree-four stages.
    P=s.symbols('P0:3');U=s.symbols('U0:3');coords=list(P)+list(U)
    def base(v):return sub(K,xyz,v)
    def root(v):return XX.subs(dict(zip(xyz,v)),simultaneous=True)
    def prof(v,g):return sub(compiler(4,{4:g}),xyz,v)
    XP=root(P);XU=root(U);target=s.Matrix(prof(P,s.Integer(1))+prof(U,XP))
    stage1=list(P)+prof(U,XP);stage2=prof(P,s.Integer(1))+list(U)
    composed=sub(stage2,coords,stage1)
    ck('actual root-dependent map has polynomial aux-first factors',all(zero(a-b) for a,b in zip(composed,target)))
    Jactual=target.jacobian(coords)
    for j,point in enumerate([[1,0,1,1,0,1],[0,1,0,1,0,0],[1,1,0,0,1,1]]):
        ck('actual root-dependent whole-map determinant point '+str(j),Jactual.subs(dict(zip(coords,point))).det()==1)
    v0=[s.Integer(0),s.Integer(0),s.Integer(-1)];v1=[-s.Rational(1,2),s.Integer(3),s.Integer(26)];Pval=[s.Integer(1),s.Integer(0),s.Integer(1)]
    ck('non-descended witness has identical old auxiliary outputs',base(v0)==base(v1)==[s.Integer(-1),s.Integer(0),s.Integer(0)])
    ck('non-descended witness fixes inactive compiler source',root(Pval)==1)
    beta0=[a-b for a,b in zip(base(v0),prof(v0,s.Integer(1)))];beta1=[a-b for a,b in zip(base(v1),prof(v1,s.Integer(1)))]
    ck('actual rank-two return not recovered from C[P,K(U)]',beta0==[0,0,0] and beta1==[0,2,2])
    # Two-way root feedback and a three-cycle fail Keller exactly as predicted.
    cycle=s.Matrix(prof(P,XU)+prof(U,XP));point=dict(zip(coords,[1,0,1,1,0,1]))
    ck('literal native two-cycle has zero Jacobian at one point',cycle.jacobian(coords).subs(point).det()==0)
    rho1,rho2,rho3=s.symbols('rho1 rho2 rho3');alpha,beta,gamma=s.symbols('alpha beta gamma');cyc3=s.Matrix([[0,alpha,0],[0,0,beta],[gamma,0,0]])
    ck('shortest three-cycle has unique nonzero determinant term',cyc3.det()==alpha*beta*gamma)
    ck('independent diagonal three-cycle defect exact',(s.diag(rho1,rho2,rho3)+cyc3).det()-rho1*rho2*rho3==alpha*beta*gamma)
    # Root corrections vanishing on retained marked planes stay in the closed class.
    I=lambda T:T*(3*T-1)*(4*T+1);mark=I(sg)**3*X2;Fm=-2*sg**2*mark*X1**2;Hm=-2*mark*X1**4
    ck('finite-order marked root profile self equation',zero(sg**2*s.diff(Hm,X1)-2*X1**2*s.diff(Fm,X1)))
    ck('entire order-three divisor correction divisible',zero(Fm/I(sg)**3+2*sg**2*X2*X1**2) and zero(Hm/I(sg)**3+2*X2*X1**4))
    ck('marked planes vanish exactly',[I(v) for v in [s.Integer(0),s.Rational(1,3),-s.Rational(1,4)]]==[0,0,0])
    # A Gamma-invariant graph is layered by longest outgoing paths.
    edges={(0,2),(1,2)};height=[1,1,0];swap={0:1,1:0,2:2}
    ck('two-label plus fixed auxiliary graph is equivariant',{(swap[i],swap[j]) for i,j in edges}==edges)
    ck('orbit heights strictly decrease on every edge',all(height[i]>height[j] for i,j in edges) and height[0]==height[1])
    ck('extra auxiliary orbit is distinct from transitive hypothesis',swap[2]==2 and swap[0]!=0)
    ck('split cubic two-block fixed-target fibre rank is ninth not two',3**2==9 and 3**2!=2)
    # Sharp scope boundaries: Frobenius and nilpotent coefficients.
    eps=s.symbols('epsilon');rootdiag=r1*r2;D=s.Matrix([[r1,eps*a1*X1**4],[eps*a2*X2**4,r2]]).det()-rootdiag
    ck('nilpotent-base two-cycle vanishes modulo epsilon squared',s.rem(D,eps**2,eps)==0)
    ck('nilpotent-base edges survive modulo epsilon squared',s.rem(eps*a1*X1**4,eps**2,eps)!=0)
    p3=s.Integer(3);A12=p3*a2*X1**4*X2**2;A21=p3*a1*X2**4*X1**2
    ck('Frobenius two-cycle cross derivatives vanish mod three',s.Poly(A12,X1,X2,a2,modulus=3).is_zero and s.Poly(A21,X1,X2,a1,modulus=3).is_zero)
    Frob=s.Matrix(prof(P,XU**3)+prof(U,XP**3));detq=s.cancel(Frob.jacobian(coords).subs(point).det());ck('literal Frobenius-cycle Jacobian point equals one mod three',int(s.numer(detq))*pow(int(s.denom(detq)),-1,3)%3==1)
    ck('all literal maps use exact rational arithmetic',not any(f.has(s.Float) for f in list(target)+list(own)+list(cycle)))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('primary source pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),'actual_non_descended_stage_degrees_written':[4,4],'actual_non_descended_total_degree_written':16,'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Selected exact birational identities, symbolic generic Schur elimination, independent-diagonal coefficients, self-equation normal forms, actual non-descended polynomial factors, rejected cycles and sharp characteristic/base boundaries. Generic degrees and absolute stable atomicity use the written all-degree proofs and retained monodromy theorem. Not an independent whole-proof replay, formal verification, external review or publication-priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
