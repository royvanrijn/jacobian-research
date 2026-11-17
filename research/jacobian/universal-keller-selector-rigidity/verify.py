#!/usr/bin/env python3
"""Selected exact controls for common-selector rigidity, not a universal proof replay."""
from pathlib import Path
import ast,hashlib,json,platform,time,itertools
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def sub(F,vs,vals):return [f.subs(dict(zip(vs,vals)),simultaneous=True) for f in F]
def deriv(F,vs,V):return s.Matrix([F]).jacobian(vs).dot(V)
def bracket(V,Y,vs):return Y.jacobian(vs)*V-V.jacobian(vs)*Y

def main():
    started=time.monotonic();x,y,z=s.symbols('x y z');xyz=[x,y,z]
    source=ROOT/'jacobian/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),env);compiler=env['source_map']
    K=s.Matrix(compiler(3,{}));ck('source cubic actual Jacobian one',zero(K.jacobian(xyz).det()-1))
    S0=x/(1+x*y);r0=1/(1+x*y);sig=K[0]
    ck('leaf x reconstruction exact',zero(S0/r0-x));ck('leaf y reconstruction exact',zero((1-r0)/S0-y));ck('leaf z reconstruction exact',zero(sig*r0**3-((1-r0)/S0)**2*(r0*r0+3*r0)-z))
    S,r,p=s.symbols('S r p');gs=s.symbols('g4:7');hs=s.symbols('h4:7')
    A=sum(gs[j-4]*p**j*S**j for j in range(4,7));H=sum(hs[j-4]*p**j*S**j for j in range(4,7));Jc=S*S*s.diff(A,S,2)-S*s.diff(A,S)
    B=(r-1-3*p*S*S-s.diff(A,S))/(2*S);C=S*(r+1)-p*S**3-(S*s.diff(A,S)-2*A);J=s.Matrix([B,C]).jacobian([S,r]);V=s.Matrix([-s.diff(H,S)/(2*S),-(S*s.diff(H,S)-2*H)])
    ck('generic three-profile leaf Jacobian exact',zero(J.det()+r/S))
    W=s.Matrix([H/r,H/S-s.diff(H,S)+H*(-1/S+3*p*S+Jc/S**2)/r]);actual=J.adjugate()*V/(-r/S)
    ck('generic unrelated direction leaf response exact',all(zero(a-b) for a,b in zip(actual,W)))
    flip=s.Matrix([W[0].subs(r,-r),-W[1].subs(r,-r)]);Y=s.Matrix([H/r,H/S-s.diff(H,S)])
    ck('odd response under coordinate sign exactly isolated',all(zero(v) for v in (W-flip)/2-Y))
    ck('even response has only reciprocal vertical component',zero(((W+flip)/2)[0]))
    cfun=s.Function('c')(S);hfun=s.Function('h')(S);F=(S*S*s.diff(cfun,S,2)-S*s.diff(cfun,S))*hfun/S**2;ell=hfun/S-s.diff(hfun,S)
    YY=s.Matrix([hfun/r,ell]);VV=s.Matrix([0,F/r]);BB=bracket(VV,YY,[S,r]);expected=s.Matrix([-hfun*F/r**3,(ell*F-hfun*s.diff(F,S))/r**2])
    ck('formal reciprocal-field commutator identity',all(zero(v) for v in BB-expected))
    ck('formal two-vector determinant identity',zero(s.det(s.Matrix.hstack(YY,BB))-hfun*(2*ell*F-hfun*s.diff(F,S))/r**3))
    poly_controls=[]
    for c,h in [(S**4,S**4),(S**4+2*S**7,3*S**5+S**8),(S**5-S**8,S**4+S**9)]:
        jc=S*S*s.diff(c,S,2)-S*s.diff(c,S);f=jc*h/S**2;l=h/S-s.diff(h,S);P0=s.cancel(h**3*jc/S**4);det=h*(2*l*f-h*s.diff(f,S))/r**3
        ck('polynomial determinant log derivative '+str(len(poly_controls)),zero(det+h*h*f*s.diff(P0,S)/(r**3*P0)))
        ck('polynomial determinant nonzero '+str(len(poly_controls)),not zero(det))
        order=min(m[0][0] for m in s.Poly(P0,S).terms());ck('polynomial positive-order bound '+str(len(poly_controls)),order>=12)
        poly_controls.append({'c':str(c),'h':str(h),'P0_order':order})
    for n in range(1,5):
        SS=s.symbols('S0:'+str(n));rr=s.symbols('r0:'+str(n));ff=[SS[i]**(i+4)+SS[i]**(i+5) for i in range(n)];coords=list(SS)+list(rr);DV=s.Matrix([0]*n+[ff[i]/rr[i] for i in range(n)]);tt=rr[0]**2/(2*ff[0]);ds=[rr[i]**2/(2*ff[i])-tt for i in range(1,n)]
        ck('reciprocal field rank'+str(n)+' time derivative one',zero(deriv(tt,coords,DV)-1))
        ck('reciprocal field rank'+str(n)+' exact rational constants',all(zero(deriv(d,coords,DV)) for d in ds))
        ck('reciprocal field rank'+str(n)+' independent sign invariance',all(zero(d.subs(rr[j],-rr[j])-d) for d in ds for j in range(n)))
        # Exact finite valuation controls; arbitrary rank follows from the written unique-factor argument.
        tt0=s.symbols('time');branches=list(range(n));rad=[tt0+j for j in branches]
        def odd_branch(mask):
            poly=s.Poly(s.prod(rad[i] for i in range(n) if mask>>i&1),tt0)
            orders=[]
            for j in branches:
                factor=s.Poly(tt0+j,tt0);quotient=poly;order=0
                while quotient.eval(-j)==0:
                    quotient=quotient.exquo(factor);order+=1
                orders.append(order)
            return any(order%2 for order in orders)
        ck('rank'+str(n)+' actual polynomial square classes have odd branch',all(odd_branch(mask) for mask in range(1,2**n)))
    # The leading profile can be Delta-constant; the proof must still peel its block.
    u,v,t=s.symbols('u v t');SS=s.symbols('s0:3');rr=s.symbols('r0:3');pp=s.symbols('p0:3');AA=[u*pp[0]**4*SS[0]**4+v**20*pp[0]**5*SS[0]**5,u**3*pp[1]**4*SS[1]**4+v*pp[1]**6*SS[1]**6,v*pp[2]**4*SS[2]**4];HH=[s.diff(a,u) for a in AA];WW=[]
    for Si,ri,pi,ai,hi in zip(SS,rr,pp,AA,HH):
        ji=Si*Si*s.diff(ai,Si,2)-Si*s.diff(ai,Si);WW.extend([hi/ri,hi/Si-s.diff(hi,Si)+hi*(-1/Si+3*pi*Si+ji/Si**2)/ri])
    def weighted(poly):return s.Poly(s.expand(poly.subs({u:t**10,v:t},simultaneous=True)),t)
    ck('arbitrary-control active set exact',[i for i,h in enumerate(HH) if h!=0]==[0,1])
    top=max(weighted(f).degree() for f in WW if f!=0);topvec=[weighted(f).coeff_monomial(t**top) for f in WW]
    ck('first leading term isolates highest feedback product',top==50 and all(zero(f) for i,f in enumerate(topvec) if i!=3) and not zero(topvec[3]))
    nexttop=max(weighted(f).degree() for f in WW[:2]);ck('second leading product can start from invariant top profile',nexttop==20 and s.diff(AA[0].coeff(v,20),u)==0)
    ck('inactive block response exactly zero',WW[4:]==[0,0])
    # Actual mixed active/inactive feedback with signal cycles and an arbitrary allowed selector.
    P=[s.symbols('P0:3'),s.symbols('Q0:3')];U=s.symbols('U0:3');coords=list(P[0])+list(P[1])+list(U);KU=s.Matrix(sub(K,xyz,U));mu=s.Matrix([s.Rational(1,3),-1,s.Rational(2,3)]);R1=KU[1]+3*KU[0];R2=KU[2]-2*KU[0];sigma=[sub(K,xyz,q)[0] for q in P]
    F1=sub(compiler(4,{4:U[0]+sigma[1]*U[1]}),xyz,P[0]);F2=sub(compiler(5,{4:s.Integer(0),5:R1+sigma[0]*R2}),xyz,P[1]);beta=sigma[0]+P[1][0]**2+P[1][1]*P[1][2];T=s.Matrix(F1+F2+list(KU-beta*mu))
    first=F1+list(P[1])+list(U);betarec=P[0][0]+P[1][0]**2+P[1][1]*P[1][2];second=list(P[0])+list(P[1])+list(KU-betarec*mu);r1=U[1]+3*U[0];r2=U[2]-2*U[0];third=list(P[0])+sub(compiler(5,{4:s.Integer(0),5:r1+P[0][0]*r2}),xyz,P[1])+list(U)
    composed=sub(sub(third,coords,second),coords,first)
    ck('actual whole-source mixed three-stage composition exact',all(zero(f-g) for f,g in zip(composed,T)))
    ck('actual mixed profiles have a signal dependency cycle',F1[1].has(*P[1]) and F2[1].has(*P[0]))
    point=dict(zip(coords,[1,0,1,0,0,0,0,0,0]));other=dict(zip(coords,[0,1,0,1,0,1,1,0,0]));JT=T.jacobian(coords)
    ck('actual mixed arbitrary selector Jacobian first point',JT.subs(point).det()==1)
    ck('actual mixed arbitrary selector Jacobian second point',JT.subs(other).det()==1)
    for index,stage in enumerate([first,second,third]):ck('actual polynomial stage'+str(index)+' Jacobian control',s.Matrix(stage).jacobian(coords).subs(point).det()==1)
    bad=s.Matrix(F1+F2+list(KU-P[0][0]*mu));bad_det=bad.jacobian(coords).subs(point).det();ck('active source-coordinate selector fails exactly5over3',bad_det==s.Rational(5,3))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'polynomial_commutator_controls':poly_controls,'mixed_active_blocks':[0],'mixed_inactive_blocks':[1],'actual_forbidden_selector_jacobian':str(bad_det),'mixed_stage_degrees_written':[4,3,5],'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact leaf response/reconstruction, generic commutator identity, selected nonzero determinants and reciprocal invariants, leading-term boundary controls and a literal mixed nine-variable three-stage factorization. No universal replay of the multiquadratic constant-field lemma, all-profile common kernels, flat regular descent, monodromy or the atomic-chain assertion. No independent whole-proof verification, formalization, external review or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],'bad_jacobian':str(bad_det)}))
if __name__=='__main__':main()
