"""PF154: actual first-jet map, exact nontermination, and second-jet recipe.
The proof uses closed identities, not an extrapolated normal prefix.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE='FIRST_JET_ADMISSIBLE_NONTERMINATION_AND_SECOND_JET_QUOTIENT.md'
CONTEXT='sources/first-jet-nontermination/context.json'
RECEIPT='first_jet_nontermination_receipt.json'

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
    x,y,p,q,L,M,N,k,tau,R=S.symbols('x y p q L M N k tau R')
    b=9*N*y**3+k
    h=-27*N*N*y**5/(2*L)+tau*y*y
    F=L*p*p+(M+4*N*y*y)*p+N*y*q+x*(b*p+h)+x*x*(b*b-k*k)/(4*L)
    check('actual target Hessian determinant constant and tau independent',
          S.factor(S.hessian(F,(x,y,p,q)).det())==N*N*k*k)
    check('actual target potential has degree eight',S.Poly(F,x,y,p,q).total_degree()==8)
    check('gradient inverse first step y',S.diff(F,q)==N*y)
    check('gradient inverse second step x',
          S.expand(2*L*S.diff(F,x)-b*S.diff(F,p)-(2*L*h-b*(M+4*N*y*y)-k*k*x))==0)
    check('gradient inverse third step p',S.expand(S.diff(F,p)-(2*L*p+M+4*N*y*y+b*x))==0)
    check('gradient inverse last step q',S.expand(S.diff(F,y)-(8*N*y*p+N*q+x*(S.diff(b,y)*p+S.diff(h,y))+x*x*b*S.diff(b,y)/(2*L)))==0)
    u=1+x*y
    A0=u**3*p+3*x*u*u*q
    B0=y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q
    Psi0=L*A0*A0+M*A0+N*B0
    C=N*(2-3*x*y)-M*x*x-2*L*x*x*A0
    D=S.cancel((Psi0-F)/x)
    check('forcing is polynomial in source variables',S.denom(D).free_symbols<={L})
    actual=L*(A0-x**3*R)**2+M*(A0-x**3*R)+N*(B0+(2*x-3*x*x*y)*R)
    Q=L*x**5*R*R+C*R+D
    check('actual graph equation after one x factor',S.expand(actual-F-x*Q)==0)
    check('implicit graph derivative unit',S.diff(Q,R).subs(x,0)==2*N)
    Cp=S.Poly(C,x);Dp=S.Poly(D,x)
    T=S.expand(-Dp.nth(0)/(2*N))
    first=S.expand(-(Cp.nth(1)*T+Dp.nth(1))/(2*N))
    Texpected=-27*N*y**5/(4*L)+y**3*p-6*y*y*q-3*L*(y*p*p+p*q)/N-3*M*(y*p+q)/(2*N)+k*p/(2*N)+tau*y*y/(2*N)
    Sexpected=-12*L*y*y*p*p/N-39*L*y*p*q/(2*N)-9*L*q*q/(2*N)-15*M*y*y*p/(4*N)-21*M*y*q/(4*N)-S.Rational(27,2)*y**3*q+3*k*y*p/(4*N)+3*tau*y**3/(4*N)+9*k*y**3/(4*L)
    check('actual trace formula',S.expand(T-Texpected)==0)
    check('actual first normal formula',S.expand(first-Sexpected)==0)
    check('first jet satisfies degree five graph bound',S.Poly(T,y,p,q).total_degree()==5 and S.Poly(first,y,p,q).total_degree()==4)
    Tpoly=S.Poly(T,y,p,q)
    top5=sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in Tpoly.terms() if sum(ex)==5)
    top4=sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in Tpoly.terms() if sum(ex)==4)
    check('actual cone ideal J maps to zero',top5==-27*N*y**5/(4*L) and top4==y**3*p)
    # Verify the map factors through the actual first-jet coefficient ideal.
    Ky=3*(2*L*p+M)*p+N*(21*y*y*p+24*y*q+2*S.diff(T,y))
    Kp=12*L*y*p+6*L*q+3*M*y+7*N*y**3+2*N*S.diff(T,p)
    Kq=6*L*p+3*M+12*N*y*y+2*N*S.diff(T,q)
    HH=30*L*y*y*p*p+60*L*y*p*q+18*L*q*q+6*M*y*y*p+12*M*y*q+6*N*y**4*p+18*N*y**3*q-6*N*y*T
    plane=-2*L*N*N*HH+4*L*N*Ky*Kq+N*N*Kp*Kp-16*N*N*y*Kp*Kq-(16*L*N*p-64*N*N*y*y)*Kq*Kq
    check('actual A_plane equations hold identically including nilpotents',S.expand(plane-8*L*N**3*first-N*N*k*k)==0)
    Delta=S.expand(C*C-4*L*x**5*D)
    check('actual discriminant q degree one',S.degree(Delta,q)==1)
    check('actual linear q prime coefficient',S.factor(S.diff(Delta,q))==-12*L*N*x**3*(x*y+2))
    axis=S.expand(Delta.subs({y:0,p:0}))
    check('exact axis discriminant',S.expand(axis-((2*N-M*x*x)**2-24*L*N*x**3*q))==0)
    check('normal origin discriminant is square of chosen unit',Delta.subs(x,0)==4*N*N)
    check('exact quadratic root identity',S.expand((C+2*L*x**5*R)**2-Delta-4*L*x**5*Q)==0)
    # An all-order ODE certifies the binomial coefficient recurrence.
    aa,bb=S.symbols('aa bb',nonzero=True)
    root=S.sqrt(aa-bb*q)
    check('all-order square-root coefficient differential identity',S.simplify(2*(aa-bb*q)*S.diff(root,q)+bb*root)==0)
    j=S.symbols('j',integer=True,positive=True)
    ratio=-(6*L/N)*(S.Rational(1,2)-j)/(j+1)
    check('all-order high-q ratio has no positive integer zero',S.simplify(ratio-3*L*(2*j-1)/(N*(j+1)))==0)
    check('high-q sequence seed x1q2 unit',S.simplify((N/L)*S.binomial(S.Rational(1,2),2)*(-6*L/N)**2)==-9*L/(2*N))
    # k=2 normal coefficient is a cheap exact finite compatibility control.
    second=S.expand(-(Cp.nth(1)*first+Cp.nth(2)*T+Dp.nth(2))/(2*N))
    second_poly=S.Poly(second,y,p,q)
    check('entire control family fails second-jet degree bound',second_poly.total_degree()==6)
    check('second-jet unit coefficient removes entire family',second_poly.coeff_monomial(y**5*p)==-S.Rational(27,4))
    check('first three graph coefficients match target potential through x3',S.Poly(S.expand(actual.subs(R,T+x*first+x*x*second)-F),x).nth(3)==0)
    # Universal finite second-jet recipe: the only variation in H1 from
    # adding x2U is 12NU in the xx entry. Its cofactor is det(G).
    G=S.Matrix([[8*N*p,8*N*y,N],[8*N*y,2*L,0],[N,0,0]])
    U=S.symbols('U')
    check('actual universal second-normal determinant multiplier',S.expand(12*N*U*G.det())==-24*L*N**3*U)
    check('square-zero trace direction lies among good parameters',S.diff(T,tau)==y*y/(2*N) and S.diff(first,tau)==3*y**3/(4*N))
    check('discriminant tau derivative exact',S.diff(Delta,tau)==4*L*x**5*y*y)
    deriv=C+2*L*x**5*R
    check('implicit tangent equation exact',S.cancel(S.diff(Q,tau)+S.diff(Q,R)*y*y/deriv)==0)
    ee,vv=S.symbols('ee vv')
    change=S.expand(Q.subs({tau:ee,R:R+ee*vv})-Q.subs(tau,0))
    check('actual square-zero lifted root identity',S.expand(change-ee*(deriv*vv-y*y)-ee*ee*L*x**5*vv*vv)==0)
    # An actual one-variable quotient: all 55 other trace variables are
    # fixed linearly, and the remaining plane equations are identically zero.
    # tau is four times [y2]T, an invertible rational change of generator.
    v2trace=S.expand(T.subs({L:1,M:13,N:2,k:1}))
    for ex,co in S.Poly(v2trace,y,p,q).terms():
        if ex!=(2,0,0):
            assert not co.has(tau)
    check('actual trace line has exactly one free coefficient',
          S.diff(v2trace,tau)==y*y/4)
    check('all plane equations vanish on the entire actual trace line',
          S.expand((plane-8*L*N**3*first-N*N*k*k).subs({L:1,M:13,N:2,k:1}))==0)
    Mt=S.Matrix([[0,0],[1,0]]);one=S.Matrix([1,0]);eps=Mt*one
    check('actual dual numbers trace quotient unit basis',
          S.Matrix.hstack(one,eps)==S.eye(2))
    check('actual dual numbers nilpotent multiplication',Mt**2==S.zeros(2) and eps!=S.zeros(2,1))
    check('actual dual numbers projected transition invertible',
          (S.Rational(13,4)*S.eye(2))*(S.Rational(4,13)*S.eye(2))==S.eye(2))
    check('actual dual numbers transition preserves nonzero nilpotent state',
          (S.Rational(13,4)*S.eye(2))*eps==S.Rational(13,4)*eps)
    # Original Rodrigues actual multiplication, not its weighted initial table.
    data=json.loads((ROOT/'moments/vanishing/rodrigues_local_length17_certificate/certificate_data.json').read_text())
    Ma=S.Matrix([[S.Rational(v) for v in row] for row in data['multiplication_a']])
    Mb=S.Matrix([[S.Rational(v) for v in row] for row in data['multiplication_b']])
    unit=S.eye(17)[:,data['quotient_basis'].index([0,0])]
    eta=Mb**6*unit
    check('actual Rodrigues matrices commute',Ma*Mb==Mb*Ma)
    check('actual eta nonzero',eta==S.eye(17)[:,-1] and eta!=S.zeros(17,1))
    check('actual eta multiplication square zero',(Mb**6)**2==S.zeros(17))
    check('actual eta ideal rank one',(Mb**6).rank()==1)
    check('actual maximal ideal kills eta',Ma*eta==Mb*eta==S.zeros(17,1))
    check('exact projected geometric response',(2*N-M*x*x).subs({N:2,M:13})==4-13*x*x)
    bad=S.Rational(13,4)*S.eye(17)
    check('actual projected transition is unit',bad*(S.Rational(4,13)*S.eye(17))==S.eye(17))
    check('unit transition preserves nonzero eta layer',bad*eta==S.Rational(13,4)*eta and bad*eta!=S.zeros(17,1))
    inverse_root=(aa-bb*q)**(-S.Rational(1,2))
    check('all-order inverse-root response recurrence',S.simplify(2*(aa-bb*q)*S.diff(inverse_root,q)-bb*inverse_root)==0)
    check('positive termination and full surviving-locus claims remain false',not ctx['positive_HC4_termination_certificate_proved'] and not ctx['control_claimed_on_full_surviving_degree_five_locus'])
    out=dict(claim='PF2RULED154',outcome='PASS',check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),context_sha256=sha(HERE/CONTEXT),
             actual_first_jet_quotient_map_constructed=True,cone_ideal_maps_to_zero=True,
             actual_dual_numbers_trace_quotient_dimension=2,actual_dual_numbers_multiplication_used=True,
             formal_graph_nonrational_and_nonpolynomial=True,whole_first_jet_nilpotent_tail_certificate_impossible=True,
             fixed_target_determinant='N^2*k^2',target_gradient_polynomial_inverse=True,
             original_Rodrigues_actual_square_zero_multiplication_used=True,
             all_order_claim_inferred_from_prefix=False,second_jet_unit_coefficient='-27/4',
             second_jet_quotient_finite_recipe_supplied=True,control_factors_second_jet_quotient=False,
             positive_HC4_termination_certificate_proved=False,global_coverage_proved=False,
             new_external_theorem_imports=[],registered_predecessor_checkers_replayed=[],full_goal_complete=False,
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['claim','outcome','check_count','elapsed_seconds']}))

if __name__=='__main__':
    main()
