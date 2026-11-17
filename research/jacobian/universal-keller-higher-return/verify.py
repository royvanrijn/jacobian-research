#!/usr/bin/env python3
"""Selected exact higher-return, marked-deformation and non-descent controls."""
from pathlib import Path
from itertools import product
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def sub(F,vs,vals):return [f.subs(dict(zip(vs,vals)),simultaneous=True) for f in F]
def main():
    started=time.monotonic();x,y,z=s.symbols('x y z');xyz=[x,y,z]
    source=ROOT/'jacobian/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));DK=K.jacobian(xyz)
    ck('actual cubic Jacobian one',zero(DK.det()-1));M=s.Matrix([[0,0],[1,0],[0,1]]);Delta=DK.adjugate()*M
    ck('actual rank-two translation frame contraction',all(zero(v) for v in DK*Delta-M))
    ck('actual rank-two quotient fixed',all(zero(v) for v in s.Matrix([K[0]]).jacobian(xyz)*Delta))
    ck('actual rank-two time derivative matrix identity',all(zero(v) for v in s.Matrix([K[1],K[2]]).jacobian(xyz)*Delta-s.eye(2)))
    bracket=Delta[:,1].jacobian(xyz)*Delta[:,0]-Delta[:,0].jacobian(xyz)*Delta[:,1]
    ck('actual rank-two auxiliary frame commutes',all(zero(v) for v in bracket))
    S=x/(1+x*y);ck('rank-two quotient inverse cubic relation exact',zero(S+K[1]*S*S+K[0]*S**3-K[2]/2))
    # Ordinary time nilpotence has invariant coefficients; both off-diagonal entries are nonzero.
    p,q,T1,T2=s.symbols('p q T1 T2');N=s.Matrix([[p*q,p*p],[-q*q,-p*q]]);tv=s.Matrix([T1,T2]);F=(s.eye(2)-N)*tv
    ck('cyclic polynomial coefficient matrix square zero',N*N==s.zeros(2));ck('cyclic matrix has both directions',N[0,1]!=0 and N[1,0]!=0)
    ck('cyclic time map actual Jacobian one',s.factor(F.jacobian([T1,T2]).det())==1)
    ck('cyclic time map actual polynomial inverse',all(zero(v) for v in (s.eye(2)+N)*F-tv))
    # Literal six-variable descended rank-two feedback factors through K.
    P=s.symbols('P0:3');U=s.symbols('U0:3');coords=list(P)+list(U);KU=s.Matrix(sub(K,xyz,U));R=KU[0];pp=P[0]+1;qq=R+1;NN=s.Matrix([[pp*qq,pp*pp],[-qq*qq,-pp*qq]]);times=s.Matrix([KU[1],KU[2]]);newtimes=(s.eye(2)-NN)*times
    FP=sub(compiler(4,{4:R}),xyz,P);target=s.Matrix(FP+[R,newtimes[0],newtimes[1]])
    stage1=list(P)+list(KU);p0=P[0]+1;q0=U[0]+1;N0=s.Matrix([[p0*q0,p0*p0],[-q0*q0,-p0*q0]]);time0=(s.eye(2)-N0)*s.Matrix([U[1],U[2]]);stage2=list(P)+[U[0],time0[0],time0[1]];stage3=sub(compiler(4,{4:U[0]}),xyz,P)+list(U)
    composed=sub(sub(stage3,coords,stage2),coords,stage1)
    ck('actual descended rank-two whole-map composition',all(zero(a-b) for a,b in zip(composed,target)))
    point=dict(zip(coords,[1,0,1,0,0,0]));other=dict(zip(coords,[0,1,0,1,0,1]));J=target.jacobian(coords)
    ck('actual descended rank-two Jacobian first point',J.subs(point).det()==1);ck('actual descended rank-two Jacobian second point',J.subs(other).det()==1)
    ck('actual descended time stage Jacobian one',s.Matrix(stage2).jacobian(coords).subs(point).det()==1)
    # Published Wang-Yan family, explicit inverse for the checked P(w)=w^2 member.
    sparam=s.symbols('s');w=y+x*x;Pw=w*w;H=s.Matrix([Pw,z-2*x*Pw,Pw*Pw]);WF=s.Matrix(xyz)+sparam*H;JH=H.jacobian(xyz)
    ck('published time displacement Jacobian trace zero',zero(s.trace(JH)))
    ck('published time displacement Jacobian nilpotent',all(zero(v) for v in JH*JH*JH))
    recovered=s.expand(WF[1]+WF[0]**2-sparam*WF[2]);ck('published time invariant exact',recovered==w)
    inverse_after=[WF[0]-sparam*recovered**2,recovered-(WF[0]-sparam*recovered**2)**2,WF[2]-sparam*recovered**4]
    ck('published time family actual polynomial inverse',all(zero(a-b) for a,b in zip(inverse_after,xyz)))
    ck('published time family actual Jacobian one',zero(WF.jacobian(xyz).det()-1))
    # Universal divisor-preserving coefficient construction.
    pi,a,beta=s.symbols('pi a beta');I=lambda v:v*(3*v-1)*(4*v+1);sig=K[0];GI=s.Matrix(compiler(4,{4:I(sig)}));const4=s.Matrix(compiler(4,{4:s.Integer(1)}))
    ck('profile correction is exactly divisible by retained divisor',all(zero(v) for v in GI-K-I(sig)*(const4-K)))
    ck('all marked first values annihilate divisor',[I(v) for v in [s.Integer(0),s.Rational(1,3),-s.Rational(1,4)]]==[0,0,0])
    Q=s.div(s.Poly(I(pi)-I(a),pi),s.Poly(pi-a,pi))[0].as_expr();ck('universal fibre ideal divided-difference identity',zero(I(pi)-I(a)-(pi-a)*Q))
    ck('idempotent selector kills auxiliary divisor universally',zero(I(beta/3)-beta*(beta-1)*(4*beta+3)/9))
    for i,pt in enumerate([[0,0,0],[1,0,1],[0,1,0]]):ck('actual marked-preserving degree-four Jacobian point'+str(i),GI.jacobian(xyz).subs(dict(zip(xyz,pt))).det()==1)
    p0=[0,0,-s.Rational(1,4)];p1=[1,-s.Rational(3,2),s.Rational(13,2)];p2=[-1,s.Rational(3,2),s.Rational(13,2)];rho=[-s.Rational(1,4),0,0]
    ck('all three old split marked points retained',all(sub(GI,xyz,pt)==rho for pt in [p0,p1,p2]))
    S0=s.symbols('S');Emu=S0-S0*S0+S0**3/3-s.Rational(1,3)
    ck('old empty-fibre inverse polynomial unchanged at every degree',I(s.Rational(1,3))==0 and s.expand(Emu-(S0-1)**3/3)==0)
    ck('empty-fibre derivative universally nilpotent',zero(s.diff(Emu,S0)**2-3*(S0-1)*Emu))
    # Rank-two non-descent witnessed on a complete old-output fibre away from the retained divisor.
    v0=[0,0,-1];v1=[-s.Rational(1,2),3,26]
    ck('non-descent witness has identical old cubic outputs',sub(K,xyz,v0)==sub(K,xyz,v1)==[-1,0,0])
    ck('non-descent witness has distinct new marked-preserving outputs',sub(GI,xyz,v0)==[-1,0,0] and sub(GI,xyz,v1)==[-1,24,24])
    ck('native deformation has exactly two return components',zero(GI[0]-K[0]) and not zero(GI[1]-K[1]) and not zero(GI[2]-K[2]))
    ck('generic degree4 is not total coordinate degree4',max(s.Poly(f,*xyz).total_degree() for f in const4)==26)
    # Exact hyperbola-complement leaf model and its volume form.
    aa=s.symbols('aa');t=1+x*y;zz=(aa-t*y*y*(1+3*t))/t**3;leaf=sub(K,xyz,[x,y,zz]);ck('leaf localization inverse satisfies retained signal',zero(leaf[0]-aa))
    qleaf=t*t*zz+y*y*(1+3*t);ck('leaf retained t is an actual unit',zero(t*qleaf/aa-1))
    B=-y/2-s.Rational(3,2)*aa*x/t;C=x*(t+1)/t**2-aa*x**3/t**3
    ck('leaf output formulas exact',zero(leaf[1]-B) and zero(leaf[2]-C))
    ck('leaf volume form exact inverse t cubed',zero(s.Matrix([B,C]).jacobian([x,y]).det()-1/t**3))
    # Full-rank cancellation gives a known atomic profile, not a scalar obstruction extension.
    FPactive=sub(compiler(4,{4:U[0]}),xyz,P);betavec=KU-s.Matrix(U);cancelled=s.Matrix(FPactive+list(KU-betavec))
    ck('full-rank return erases original auxiliary map literally',list(cancelled)==FPactive+list(U))
    ck('full-rank known-profile Jacobian one',cancelled.jacobian(coords).subs(point).det()==1)
    Z0=sub(FPactive,coords,[0,0,0]+p0);Z1=sub(FPactive,coords,[0,0,0]+p1)
    ck('full-rank selector fails processed-output descent on exact fibre',Z0==Z1 and sub(K,xyz,p0)==sub(K,xyz,p1) and sub(list(betavec),U,p0)!=sub(list(betavec),U,p1))
    # Native profile plus nonzero trace-filter direction requires three constant directions.
    bb,coef,ps=s.symbols('bb coef ps');mu=s.Matrix([s.Rational(1,3),-1,s.Rational(2,3)]);returnvec=bb*mu-s.Matrix([0,-2*coef*ps**4*S0**2,-2*coef*ps**4*S0**4]);cols=[returnvec.applyfunc(lambda f:s.expand(f).coeff(S0,j)) for j in [0,2,4]]
    ck('native combined return has full constant span',s.factor(s.det(s.Matrix.hstack(*cols)))==s.Rational(4,3)*bb*coef**2*ps**8)
    # Actual N=2 canonical deformation with a source-dependent invariant parameter.
    PP=[s.symbols('A0:3'),s.symbols('B0:3')];UU=s.symbols('V0:3');allvars=list(PP[0])+list(PP[1])+list(UU);KU2=s.Matrix(sub(K,xyz,UU));sigU=KU2[0];ell=lambda v:s.Rational(1,2)*(3*v*v-v);Ctrace=sum(ell(q[0]) for q in PP);lam=lambda T:-s.Rational(1,6)*T*(T-2)*(T-3)*(T-4);betatrace=1-lam(Ctrace)
    Gsource=sub(compiler(4,{4:I(sigU)*(Ctrace+1)}),xyz,UU);firstout=[v for q in PP for v in sub(K,xyz,q)];receiver=s.Matrix(firstout+list(s.Matrix(Gsource)-betatrace*mu));Haux=list(PP[0])+list(PP[1])+list(s.Matrix(Gsource)-betatrace*mu);outer=firstout+list(UU)
    ck('actual source-parameter canonical reverse factors exact',all(zero(f-g) for f,g in zip(sub(outer,allvars,Haux),receiver)))
    ck('actual source-parameter canonical Jacobian one',receiver.jacobian(allvars).subs(dict(zip(allvars,[0]*9))).det()==1)
    configs=list(product([p0,p1,p2],repeat=2));selected=[];all_idempotent=True;all_coeff_zero=True
    for i,qs in enumerate(configs):
        cc=sum(ell(q[0]) for q in qs);bv=s.cancel(1-lam(cc));all_idempotent &= bv in [0,1];all_coeff_zero &= I(bv/3)*(cc+1)==0
        if bv==0:selected.append(i)
    ck('canonical map uses exact rational arithmetic',not any(f.has(s.Float) for f in receiver))
    ck('all split N2 selector values remain idempotent',all_idempotent)
    ck('all split N2 auxiliary coefficient corrections vanish on fibre',all_coeff_zero)
    ck('canonical deformation keeps exactly natural two labels',len(selected)==2 and [configs[i] for i in selected]==[(p0,p1),(p1,p0)])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'descended_rank2_stage_degrees':[3,1,4],'non_descended_rank2_auxiliary_degrees':[3,4],'full_rank_cancelled_total_degree':4,'native_combined_return_rank':3,'canonical_N2_degree_written':36,'generic_degree_notice':'These are generic covering degrees from the written inverse/partial-map arguments and retained monodromy theorem, not total coordinate degrees or inferences from Jacobian samples.','python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact auxiliary frame, descended nilpotent time factors, a published nonlinear time member, universal divisor/idempotent identities, actual marked deformation maps, exact non-descent witnesses, hyperbola-leaf volume and literal canonical factors. Not universal replay of multi-derivation Bertini, stable-polynomial factors, all-degree atomicity or arbitrary-base canonical descent. No independent whole-proof verification, formalization, external review or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
