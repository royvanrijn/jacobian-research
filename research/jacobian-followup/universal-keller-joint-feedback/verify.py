#!/usr/bin/env python3
"""Selected exact joint-feedback/Bertini-hypothesis controls; not universal theorem replay."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def sub(F,vs,vals):return [f.subs(dict(zip(vs,vals)),simultaneous=True) for f in F]
def deriv(f,vs,V):return s.Matrix([f]).jacobian(vs).dot(V)
def main():
    started=time.monotonic();x,y,z=s.symbols('x y z');xyz=[x,y,z]
    source=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));DK=K.jacobian(xyz)
    ck('actual auxiliary cubic Jacobian one',zero(DK.det()-1))
    # Hypotheses of the imported field lemma, not sampled geometric irreducibility.
    S,a,p,B,C,X1,X2=s.symbols('S a p B C X1 X2')
    cases=[{4:a},{4:a*a+1,5:a**3+a},{4:(a*a+1)/(a-1),5:a**3}]
    for index,g in enumerate(cases):
        f=S+p*S**3+sum(gj*p**j*S**j for j,gj in g.items());u=S*S;J=s.Matrix([u,f]).jacobian([S,a]).det()
        ck('Bertini active independence '+str(index),not zero(J) and zero(J-2*S*s.diff(f,a)))
        four=s.Matrix([X1,X2,X1*u+X2*f,u]).jacobian([X1,X2,S,a]).det()
        ck('Bertini four-element independence '+str(index),not zero(four))
        E=f+B*S*S-C/2;poly=s.Poly(E.subs({a:2,p:1,B:0,C:2}),S,domain=s.QQ)
        ck('exact specialization certifies generic chart separability '+str(index),s.gcd(poly,poly.diff()).degree()==0)
        ck('generic inverse zero root excluded '+str(index),zero(E.subs(S,0)+C/2))
    fconst=S+p*S**3+7*p**4*S**4;ck('inactive profile fails required surface independence',s.Matrix([S*S,fconst]).jacobian([S,a]).det()==0)
    # Full joint horizontal identity in leaf coordinates, independent of the source-only argument.
    r=s.symbols('r');A=a*p**4*S**4+(a*a+a)*p**5*S**5;H=s.diff(A,a);Jc=S*S*s.diff(A,S,2)-S*s.diff(A,S)
    Bout=(r-1-3*p*S*S-s.diff(A,S))/(2*S);Cout=S*(r+1)-p*S**3-(S*s.diff(A,S)-2*A)
    W=s.Matrix([H/r,H/S-s.diff(H,S)+H*(-1/S+3*p*S+Jc/S**2)/r]);D=s.Matrix([-W[0],-W[1],1]);vs=[S,r,a]
    ck('joint horizontal derivation fixes second compiler output',zero(deriv(Bout,vs,D)))
    ck('joint horizontal derivation fixes third compiler output',zero(deriv(Cout,vs,D)))
    ck('joint horizontal derivation has polynomial time slice',deriv(a,vs,D)==1)
    h=Bout+Cout
    for c in [s.Integer(1),s.Integer(2)]:
        beta=(1-c)*a+h;det=s.Matrix([Bout,Cout,a-beta]).jacobian(vs).det()
        ck('joint leaf constant Jacobian '+str(c),zero((-S/r)*det-c))
        ck('joint leaf affine-slice criterion '+str(c),zero(deriv(beta,vs,D)-(1-c)))
    ck('auxiliary slice alone makes feedback singular',zero(s.Matrix([Bout,Cout,0]).jacobian(vs).det()))
    # Actual nine-variable mixed active/inactive cubic feedback and its full polynomial factors.
    P=[s.symbols('P0:3'),s.symbols('Q0:3')];U=s.symbols('U0:3');coords=list(P[0])+list(P[1])+list(U);KU=s.Matrix(sub(K,xyz,U));mu=s.Matrix([s.Rational(1,3),-1,s.Rational(2,3)]);delta=KU.jacobian(U).adjugate()*mu
    R1=KU[1]+3*KU[0];R2=KU[2]-2*KU[0];tt=-KU[1];sigma=[sub(K,xyz,q)[0] for q in P]
    ck('actual cubic time has derivative one',zero(deriv(tt,U,delta)-1))
    ck('actual cubic polynomial constants fixed',zero(deriv(R1,U,delta)) and zero(deriv(R2,U,delta)))
    F1=sub(compiler(4,{4:U[0]+sigma[1]*U[1]}),xyz,P[0]);F2=sub(compiler(5,{4:s.Integer(0),5:R1+sigma[0]*R2}),xyz,P[1]);hjoint=F1[1]+R1*F1[2]+P[1][0]**2+P[1][1]*R2**2
    point=dict(zip(coords,[1,0,1,0,0,0,0,0,0]));other=dict(zip(coords,[0,1,0,1,0,1,1,0,0]))
    ck('actual allowed selector genuinely uses joint variables',hjoint.has(*U) and hjoint.has(*P[0]))
    ck('actual allowed joint selector not an auxiliary first integral',deriv(hjoint,U,delta).subs(point)==-s.Rational(2,3))
    first=F1+list(P[1])+list(U);hrec=P[0][1]+R1*P[0][2]+P[1][0]**2+P[1][1]*R2**2;r1=U[1]+3*U[0];r2=U[2]-2*U[0];third=list(P[0])+sub(compiler(5,{4:s.Integer(0),5:r1+P[0][0]*r2}),xyz,P[1])+list(U)
    for c in [s.Integer(1),s.Integer(2)]:
        beta=(1-c)*tt+hjoint;T=s.Matrix(F1+F2+list(KU-beta*mu));second=list(P[0])+list(P[1])+list(KU-((1-c)*tt+hrec)*mu)
        composed=sub(sub(third,coords,second),coords,first)
        ck('actual full joint three-stage composition c'+str(c),all(zero(f-g) for f,g in zip(composed,T)))
        JT=T.jacobian(coords)
        ck('actual joint Jacobian first point c'+str(c),JT.subs(point).det()==c)
        ck('actual joint Jacobian second point c'+str(c),JT.subs(other).det()==c)
        ck('actual joint auxiliary stage Jacobian c'+str(c),s.Matrix(second).jacobian(coords).subs(point).det()==c)
    # Non-output-recoverable but polynomial auxiliary constants: two cubic blocks.
    PP=s.symbols('X0:3');BB=s.symbols('B0:3');VV=s.symbols('V0:3');vv=list(PP)+list(BB)+list(VV);KB=s.Matrix(sub(K,xyz,BB));KV=s.Matrix(sub(K,xyz,VV));rv1=KV[1]+3*KV[0];rv2=KV[2]-2*KV[0]
    g=BB[0]+rv1;hh=PP[0]+BB[1]*rv1+rv2*rv2;FP=sub(compiler(4,{4:g}),xyz,PP);Tprod=s.Matrix(FP+list(KB)+list(KV-hh*mu))
    stage1=list(PP)+list(BB)+list(KV-hh*mu);linearR=VV[1]+3*VV[0];stage2=sub(compiler(4,{4:BB[0]+linearR}),xyz,PP)+list(BB)+list(VV);stage3=list(PP)+list(KB)+list(VV)
    composed=sub(sub(stage3,vv,stage2),vv,stage1)
    ck('actual product-auxiliary three-stage factorization',all(zero(f-g) for f,g in zip(composed,Tprod)))
    origin=dict(zip(vv,[0]*9));sample=dict(zip(vv,[1,0,1,0,1,0,1,0,0]));JP=Tprod.jacobian(vv)
    ck('actual product-auxiliary Jacobian origin',JP.subs(origin).det()==1)
    ck('actual product-auxiliary Jacobian second point',JP.subs(sample).det()==1)
    for index,stage in enumerate([stage1,stage2,stage3]):ck('actual product-auxiliary stage Jacobian '+str(index),s.Matrix(stage).jacobian(vv).subs(origin).det()==1)
    eta=s.Matrix(list(BB)+[rv1,rv2]);ck('polynomial constant generators independent',eta.jacobian(list(BB)+list(VV)).subs(origin).rank()==5)
    q0=[s.Integer(0),s.Integer(0),-s.Rational(1,4)];q1=[s.Integer(1),-s.Rational(3,2),s.Rational(13,2)]
    ck('transverse constant not recoverable from original outputs',sub(K,xyz,q0)==sub(K,xyz,q1) and q0[0]!=q1[0])
    # Locally nilpotent slice example with an explicit source automorphism.
    v=s.symbols('v');cc=list(PP)+list(BB)+[v];hh0=PP[0]+BB[1]**2;c=s.Integer(2);theta=list(PP)+list(BB)+[c*v-hh0];theta_inv=list(PP)+list(BB)+[(v+hh0)/c]
    ck('LND slice source change actual polynomial automorphism',all(zero(f-g) for f,g in zip(sub(theta,cc,theta_inv),cc)))
    Fp=sub(compiler(4,{4:BB[0]}),xyz,PP);middle=Fp+list(BB)+[v];outer=list(PP)+list(KB)+[v];TL=s.Matrix(Fp+list(KB)+[c*v-hh0]);comp=sub(sub(outer,cc,middle),cc,theta)
    ck('LND slice exact unstabilized compiler-auxiliary factors',all(zero(f-g) for f,g in zip(comp,TL)))
    ck('LND slice example constant Jacobian two',TL.jacobian(cc).subs(dict(zip(cc,[0]*7))).det()==2)
    # Positive endpoint: one compiler block and an invertible auxiliary map, giving a literal LR equivalence.
    u0,u1=s.symbols('u0 u1');cc=list(PP)+[u0,u1];Fpos=sub(compiler(4,{4:u0}),xyz,PP);Hpos=Fpos+[u0,u1];Qpos=list(PP)+[u0-PP[1]-u1*u1,u1];Qinv=list(PP)+[u0+PP[1]+u1*u1,u1];Tpos=s.Matrix(Fpos+[u0-Fpos[1]-u1*u1,u1])
    ck('positive atomic endpoint auxiliary shear is a polynomial automorphism',all(zero(f-g) for f,g in zip(sub(Qpos,cc,Qinv),cc)))
    ck('positive atomic endpoint exact known-profile LR equivalence',all(zero(f-g) for f,g in zip(sub(Qpos,cc,Hpos),Tpos)))
    ck('positive atomic endpoint actual Jacobian one',Tpos.jacobian(cc).subs(dict(zip(cc,[0]*5))).det()==1)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'Bertini_hypothesis_controls':len(cases),'joint_constant_Jacobians':['1','2'],'mixed_cubic_stage_degrees':[4,3,5],'product_auxiliary_stage_degrees':[3,4,3],'polynomial_constant_kernel_example':'k[B_1,B_2,B_3,R_1(V),R_2(V)] strictly larger than the original polynomial output quotient, witnessed by two points. Full equality uses the written cubic-kernel argument.','python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact Bertini independence/chart controls, horizontal leaf identities, two actual nine-variable joint-Jacobian classes and literal polynomial factors, a polynomial-kernel example beyond output recovery, an actual LND source automorphism and a literal positive single-block/invertible-auxiliary LR equivalence. Not universal replay of Bertini, relative constants, polynomial saturation, stable-polynomial factorization or canonical fixed-point arguments. No independent whole-proof review, formalization or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
