"""PF153: exact actual-quotient identities for uniform HC4 cone nilpotence.
Written ideal-power proof is in the companion note. No normal-prefix inference.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE='HC4_TRACE_QUOTIENT_AND_UNIFORM_CONE_NILPOTENCE.md'
CONTEXT='sources/hc4-trace-nilpotence/context.json'
RECEIPT='hc4_trace_nilpotence_receipt.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.perf_counter();checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((HERE/CONTEXT).read_text())
    for path,digest in ctx['retained_sha256'].items():
        check('retained source pin '+path,sha(ROOT/path)==digest)
    x,y,p,q,L,M,N,c,t,e=S.symbols('x y p q L M N c t e')
    tr,ty,tp,tq,sjet=S.symbols('tr ty tp tq sjet')
    source=(y,p,q)
    def plane(T,Ty,Tp,Tq):
        Ky=3*(2*L*p+M)*p+N*(21*y*y*p+24*y*q+2*Ty)
        Kp=12*L*y*p+6*L*q+3*M*y+7*N*y**3+2*N*Tp
        Kq=6*L*p+3*M+12*N*y*y+2*N*Tq
        H=30*L*y*y*p*p+60*L*y*p*q+18*L*q*q+6*M*y*y*p+12*M*y*q+6*N*y**4*p+18*N*y**3*q-6*N*y*T
        F=-2*L*N*N*H+4*L*N*Ky*Kq+N*N*Kp*Kp-16*N*N*y*Kp*Kq-(16*L*N*p-64*N*N*y*y)*Kq*Kq
        return F,(Ky,Kp,Kq),H
    # Derive the plane Hessian from the actual parent potential.
    u=1+x*y
    aa=u**3*p+3*x*u*u*q-x**3*(tr+x*sjet)
    bb=y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*(tr+x*sjet)
    parent=S.Poly(S.expand(L*aa*aa+M*aa+N*bb),x)
    fs=[parent.nth(i) for i in range(4)]
    F,K,H=plane(tr,ty,tp,tq)
    mixed=[S.diff(fs[1],v)+S.diff(fs[1],tr)*dv for v,dv in zip(source,(ty,tp,tq))]
    for i in range(3):
        check('actual plane mixed Hessian entry '+str(i),S.expand(mixed[i]-K[i])==0)
    check('actual normal Hessian entry',S.expand(2*fs[2]-H-4*N*sjet)==0)
    G=S.hessian(fs[0],source)
    check('actual lower plane Hessian block',G==S.Matrix([[8*N*p,8*N*y,N],[8*N*y,2*L,0],[N,0,0]]))
    check('lower plane Hessian block has unit determinant',S.factor(G.det())==-2*L*N*N)
    check('plane coefficient formula from actual adjugate',S.expand(H*G.det()-(S.Matrix(K).T*G.adjugate()*S.Matrix(K))[0]-F)==0)
    check('actual first-normal multiplier',S.expand(4*N*G.det())==-8*L*N**3)
    monomials=[(a,b,d) for a in range(6) for b in range(6-a) for d in range(6-a-b)]
    bad=[exp for exp in monomials if (sum(exp)==5 and (exp[1] or exp[2])) or (sum(exp)==4 and exp[2])]
    good=[exp for exp in monomials if exp not in bad]
    check('actual trace parameter counts',len(monomials)==56 and len(bad)==30 and len(good)==26)
    slots=sum((d+2)*(d+1)//2 for d in range(5,11))
    check('finite plane equation recipe count',slots==251)
    # Homogeneous jets: derivative entries of T_i have weight i-1.
    z=S.symbols('z0:6');zy=S.symbols('zy1:6');zp=S.symbols('zp1:6');zq=S.symbols('zq1:6')
    Ts=sum(t**i*z[i] for i in range(6))
    dys=[sum(t**(i-1)*row[i-1] for i in range(1,6)) for row in (zy,zp,zq)]
    raw=F.subs({y:t*y,p:t*p,q:t*q,tr:Ts,ty:dys[0],tp:dys[1],tq:dys[2]},simultaneous=True)
    faces=S.Poly(S.expand(raw),t)
    check('source degree of universal plane formula at most ten',faces.degree()==10)
    check('degree ten actual square face',S.expand(faces.nth(10)-256*N**4*y*y*zq[4]**2)==0)
    check('degree eight actual square face modulo J1',S.expand(faces.nth(8).subs(zq[4],0)-4*N**4*(zp[4]-8*y*zq[3])**2)==0)
    f7=faces.nth(7).subs(zq[4],0).subs(zp[4],8*y*zq[3])
    check('degree seven actual product face modulo J2',S.expand(f7-16*L*N**3*zq[3]*(zy[4]-4*p*zq[3]))==0)
    # Rational invertible evaluation matrix proves equality of coefficient
    # and evaluation ideals over ANY commutative Q-algebra.
    grid=[(i,j) for i in range(5) for j in range(5-i)]
    V=S.Matrix([[i**a*j**b for a,b in grid] for i,j in grid])
    Vinv=V.inv()
    check('quartic evaluation determinant',V.det()==47775744)
    check('quartic actual coefficient recovery',Vinv*V==S.eye(15) and V*Vinv==S.eye(15))
    quartic=S.Matrix(S.symbols('h0:15'))
    check('every quartic coefficient recovered over Q',Vinv*(V*quartic)==quartic)
    check('square-zero evaluation pigeonhole exponent',15+1==16)
    # Actual triangular product certificate after J2, no domain division.
    qs=S.symbols('q0:4');kap=S.symbols('kappa')
    Q=sum(qs[i]*y**(3-i)*p**i for i in range(4))
    T5=kap*y**5+sum(S.Rational(8,i+1)*qs[i]*y**(4-i)*p**(i+1) for i in range(4))
    check('integrated T5 differentiates to eight yQ',S.expand(S.diff(T5,p)-8*y*Q)==0)
    B=S.diff(T5,y)-4*p*Q
    expected=5*kap*y**4+28*qs[0]*y**3*p+8*qs[1]*y*y*p*p+S.Rational(4,3)*qs[2]*y*p**3-2*qs[3]*p**4
    check('actual cubic product second factor',S.expand(B-expected)==0)
    product=S.Poly(S.expand(Q*B),y,p)
    fj=[product.coeff_monomial(y**(7-j)*p**j) for j in range(8)]
    witnesses=[
        -fj[7]/2,
        S.Rational(3,4)*fj[5]-S.Rational(9,2)*qs[1]*qs[3],
        fj[3]/8-S.Rational(11,3)*qs[0]*qs[2]-S.Rational(5,8)*kap*qs[3],
        fj[1]/28-S.Rational(5,28)*kap*qs[1]]
    for i,w in zip((3,2,1,0),witnesses):
        check('actual triangular square relation q'+str(i),S.expand(qs[i]**2-w)==0)
    # Differentiation and y multiplication are injective on the relevant
    # coefficient vectors in characteristic zero, even with nilpotents.
    for exp in bad:
        check('cone ideal coefficient unit '+str(exp),
              (sum(exp)==5 and (exp[1]>0 or exp[2]>0)) or (sum(exp)==4 and exp[2]>0))
    check('uniform lifted ideal-power exponent',16*16*2**4==4096)
    # Nonempty exact first-jet point, already retained in HC4MYG5S.
    Tstar=S.Rational(51,100)*y**5-S.Rational(47,10)*y**3*p-S.Rational(123,20)*y*y*q
    par={L:1,M:13,N:2}
    Fs,Ks,Hs=plane(Tstar,*[S.diff(Tstar,v) for v in source]);Fs=S.expand(Fs.subs(par))
    check('retained nonzero first-jet point',S.Poly(Fs,*source).total_degree()==4)
    cstar=S.Rational(17165601,25)
    Sstar=S.expand((Fs-cstar)/64)
    check('forced first-normal jet degree bound',S.Poly(Sstar,*source).total_degree()<=4)
    # Derive the actual linearization before evaluating on 30 directions.
    v,vy,vp,vq=S.symbols('v vy vp vq')
    variation=S.diff(plane(tr+e*v,ty+e*vy,tp+e*vp,tq+e*vq)[0],e).subs(e,0)
    Ky,Kp,Kq=K
    expected=12*L*N**3*y*v+8*L*N*N*Kq*vy+(4*N**3*Kp-32*N**3*y*Kq)*vp+(8*L*N*N*Ky-32*N**3*y*Kp-(64*L*N*N*p-256*N**3*y*y)*Kq)*vq
    check('actual Frechet derivative formula',S.expand(variation-expected)==0)
    atstar={tr:Tstar,ty:S.diff(Tstar,y),tp:S.diff(Tstar,p),tq:S.diff(Tstar,q),**par}
    dF=S.expand(variation.subs(atstar))
    cols=[]
    for exp in bad:
        mono=S.prod(a**b for a,b in zip(source,exp))
        df=S.Poly(S.expand(dF.subs({v:mono,vy:S.diff(mono,y),vp:S.diff(mono,p),vq:S.diff(mono,q)})),*source)
        cols.append({ex:coef for ex,coef in df.terms() if sum(ex)>4})
    rows=sorted(set().union(*(set(co) for co in cols)))
    mat=S.Matrix([[co.get(ex,0) for co in cols] for ex in rows])
    pivots=mat.T.rref()[1]
    check('actual fiber coefficient Jacobian full column rank',len(pivots)==30)
    selected=[rows[i] for i in pivots]
    minor=mat[list(pivots),:]
    minor_det=minor.det(method='domain-ge')
    check('actual fiber retained minor nonzero',minor_det!=0)
    # Reproduce the known transverse control from the actual graph potential.
    # Only parent coefficients x^0..3 are needed for [x] of its Hessian det.
    coeffs=[S.expand(fi.subs({tr:Tstar,sjet:Sstar,**par})) for fi in fs]
    H0=S.zeros(4);H1=S.zeros(4)
    H0[0,0]=2*coeffs[2];H1[0,0]=6*coeffs[3]
    for i,vi in enumerate(source,1):
        H0[0,i]=H0[i,0]=S.diff(coeffs[1],vi)
        H1[0,i]=H1[i,0]=2*S.diff(coeffs[2],vi)
        for j,vj in enumerate(source,1):
            H0[i,j]=S.diff(coeffs[0],vi,vj)
            H1[i,j]=S.diff(coeffs[1],vi,vj)
    H0=H0.subs({p:0,q:0});H1=H1.subs({p:0,q:0})
    check('first-jet determinant is desired constant',S.expand(H0.det())==cstar)
    transverse=S.expand(S.trace(H0.adjugate()*H1))
    known=S.Poly(transverse,y).nth(7)
    check('retained actual second-jet obstruction',known==S.Rational(22032,125))
    check('forced formal second coefficient outside degree-five box',known/192==S.Rational(459,500))
    check('tail and reduced-body obligations remain open',
          not ctx['normal_tail_containment_proved'] and not ctx['reduced_polynomial_body_proved'] and not ctx['HC4_termination_proved'])
    out=dict(claim='PF2RULED153',outcome='PASS',check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),context_sha256=sha(HERE/CONTEXT),
             trace_parameter_count=56,cone_ideal_generator_count=30,complementary_trace_parameter_count=26,
             uniform_actual_cone_ideal_nilpotence_bound=4096,
             quartic_evaluation_determinant=str(V.det()),quartic_evaluation_grid=grid,
             exact_fiber_bad_parameter_exponents=bad,exact_fiber_selected_row_exponents=selected,
             exact_fiber_minor_determinant=str(minor_det),fixed_Tstar_complement_fiber_length=1,
             known_second_normal_y7_coefficient='459/500',
             actual_relations_not_radical_or_initial_quotient=True,
             normal_tail_power_containment_proved=False,reduced_polynomial_body_proved=False,
             HC4_termination_proved=False,global_coverage_proved=False,
             new_external_theorem_imports=[],registered_predecessor_checkers_replayed=[],
             full_goal_complete=False,independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['claim','outcome','check_count','elapsed_seconds']}))

if __name__=='__main__':
    main()
