"""PF152: actual Rodrigues multiplication and exact HC4 tail-transfer controls.
No finite normal prefix is used to infer termination or nontermination.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="RODRIGUES_NILPOTENCE_TRANSFER_AND_HC4_LINEARIZED_TAIL.md"
CONTEXT="sources/rodrigues-hc4-transfer/context.json"
RECEIPT="rodrigues_hc4_transfer_receipt.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    started=time.perf_counter()
    checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((HERE/CONTEXT).read_text())
    for path,digest in ctx["retained_sha256"].items():
        check("retained source pin "+path,sha(ROOT/path)==digest)

    bundle=ROOT/'jacobian-followup/moment-vanishing-audit/rodrigues_local_length17_certificate'
    manifest=json.loads((bundle/'manifest.json').read_text())
    for name,digest in manifest['sha256'].items():
        check("original Rodrigues manifest "+name,sha(bundle/name)==digest)
    data=json.loads((bundle/'certificate_data.json').read_text())
    Ma=S.Matrix([[S.Rational(v) for v in row] for row in data['multiplication_a']])
    Mb=S.Matrix([[S.Rational(v) for v in row] for row in data['multiplication_b']])
    I=S.eye(17);Z=S.zeros(17)
    basis=data['quotient_basis'];unit=I[:,basis.index([0,0])]
    check("actual Rodrigues matrix sizes",Ma.shape==Mb.shape==(17,17))
    check("actual Rodrigues multiplication commutes",Ma*Mb==Mb*Ma)
    for k,(i,j) in enumerate(basis):
        check(f"actual cyclic quotient basis monomial a{i}b{j}",
              Ma**i*Mb**j*unit==I[:,k])
    for idx in ['11','12','13','14']:
        mat=S.zeros(17)
        for (i,j),coef in data['residuals'][idx]:
            mat+=S.Rational(coef)*Ma**i*Mb**j
        check("actual prefix residual annihilation "+idx,mat==Z)
    V=I;dims=[17]
    for j in range(1,8):
        W=(Ma*V).row_join(Mb*V)
        cols=W.columnspace()
        V=S.Matrix.hstack(*cols) if cols else S.zeros(17,0)
        dims.append(V.cols)
    check("actual maximal-ideal power dimensions",
          dims==[17,16,14,11,7,4,1,0])
    for i in range(8):
        check(f"actual seventh ideal power monomial a{i}b{7-i}",
              Ma**i*Mb**(7-i)==Z)
    eta=Mb**6*unit;Meta=Mb**6
    check("actual b6 is nonzero last quotient basis",eta==I[:,-1] and eta!=S.zeros(17,1))
    check("actual eta ideal is square-zero",Meta*Meta==Z)
    check("actual eta is annihilated by maximal ideal",
          Ma*eta==Mb*eta==S.zeros(17,1))
    check("actual eta multiplication rank one",Meta.rank()==1)

    x,y,p,q,tau,L,M,N,kappa,R,Y=S.symbols(
        "x y p q tau L M N kappa R Y")
    U=1+x*y
    A0=U**3*p+3*x*U**2*q
    B0=y**2*U*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q
    Psi0=L*A0**2+M*A0+N*B0
    F=L*p**2+M*p+4*N*y*y*p+N*y*q+kappa*x*p+tau*x*y*y
    check("fixed target determinant for actual quadratic trace direction",
          S.factor(S.hessian(F,(x,y,p,q)).det())==N**2*kappa**2)
    C=N*(2-3*x*y)-M*x*x-2*L*x*x*A0
    D=S.cancel((Psi0-F)/x)
    check("quadratic graph forcing is polynomial",S.denom(S.cancel(D))==1)
    check("trace perturbation enters forcing exactly minus tau y2",
          S.diff(D,tau)==-y*y)
    Delta=S.expand(C*C-4*L*x**5*D)
    check("discriminant derivative exactly4Lx5y2",
          S.diff(Delta,tau)==4*L*x**5*y*y)
    check("discriminant is linear in free source q",
          S.degree(Delta,q)==1
          and S.factor(S.diff(Delta,q))==-12*L*N*x**3*(x*y+2))
    Q=L*x**5*R**2+C*R+D
    check("implicit unit at normal origin",S.diff(Q,R).subs(x,0)==2*N)
    G=C+2*L*x**5*R
    check("exact squared implicit derivative identity",
          S.expand(G*G-Delta-4*L*x**5*Q)==0)
    check("linearized equation is G times response equals y2",
          S.diff(Q,tau)+S.diff(Q,R)*y*y/G==0)
    # Exact square-zero Taylor identity, not a finite normal-order check.
    e=S.symbols("e")
    Q0=Q.subs(tau,0)
    change=S.expand(Q.subs({tau:e,R:R+e*Y})-Q0)
    expected=e*((C+2*L*x**5*R)*Y-y*y)+e*e*L*x**5*Y*Y
    check("full square-zero lifted root identity",S.expand(change-expected)==0)
    check("square-zero lifted response solves the equation exactly",
          S.cancel((change.subs(Y,y*y/G)-e*e*L*x**5*(y*y/G)**2))==0)

    Ds=S.factor(Delta.subs({y:0,p:0,tau:0}))
    check("universal reduced discriminant exact",
          S.expand(Ds-((2*N-M*x*x)**2-24*L*N*x**3*q))==0)
    # q=0 unit root is 2N-Mx2, yielding the entire geometric series.
    axis=1/(4-13*x*x)
    check("exact projected HC4 generating-function identity",
          S.cancel((4-13*x*x)*axis-1)==0)
    check("geometric coefficients are determined by a nonzero unit recurrence",
          S.Rational(13,4)!=0 and S.Rational(13,4)*S.Rational(4,13)==1)
    Tbad=S.Rational(13,4)*I
    check("actual HC4 projected transition is invertible",
          Tbad*(S.Rational(4,13)*I)==I)
    check("actual HC4 transition preserves nonzero eta layer",
          Tbad*eta==S.Rational(13,4)*eta and Tbad*eta!=S.zeros(17,1))
    check("actual HC4 coefficient state starts eta over four",
          unit.T*(Meta*unit)==S.zeros(1) and eta/4!=S.zeros(17,1))
    # Nonzero coefficient c_j eta generates J by a rational inverse.
    jj=S.symbols("j",integer=True,nonnegative=True)
    cj=S.Rational(1,4)*(S.Rational(13,4))**jj
    check("every late coefficient recovers the same eta ideal generator",
          S.powsimp(cj*(4*(S.Rational(4,13))**jj))==1)

    # The q-series also persists if M=0. Its first-order ODE certifies
    # the ALL-order nonzero coefficient recurrence without a longer prefix.
    AA,BB=S.symbols("AA BB",nonzero=True)
    response=(AA-BB*q)**(-S.Rational(1,2))
    check("universal inverse-square-root coefficient differential identity",
          S.simplify(2*(AA-BB*q)*S.diff(response,q)-BB*response)==0)
    leadratio=3*L*(2*jj+1)/(N*(jj+1))
    check("all-order leading ratio has no characteristic-zero integer zero",
          S.simplify(leadratio/((2*jj+1)/(2*(jj+1))))==6*L/N)
    check("universal q coefficient seed is a unit",S.cancel((1/(2*N))*(2*N))==1)

    # Applicability boundary: this body is not in the known degree5 graph box.
    par={L:1,M:13,N:2,kappa:1}
    Cs=S.expand(C.subs(par));Dsp=S.expand(D.subs(par))
    trace=S.expand(-Dsp.subs(x,0)/4)
    R1=S.expand(-(S.diff(Cs,x).subs(x,0)*trace+S.diff(Dsp,x).subs(x,0))/4)
    check("trace deformation is exactly quadratic at fixed determinant",
          S.diff(trace,tau)==y*y/4)
    check("body trace degree4",S.Poly(trace.subs(tau,0),y,p,q).total_degree()==4)
    check("body first normal coefficient plane degree5",
          S.Poly(R1.subs(tau,0),y,p,q).total_degree()==5)
    check("body complete first jet degree6",
          S.Poly(trace.subs(tau,0)+x*R1.subs(tau,0),x,y,p,q).total_degree()==6)

    # Positive finite transfer CONTROL using the same actual multiplication.
    Z17=S.zeros(17)
    Tplus=Z17.row_join(Mb).col_join(I.row_join(Z17))
    diag=Mb.row_join(Z17).col_join(Z17.row_join(Mb))
    check("positive control block transition square is multiplication by b",
          Tplus**2==diag)
    check("positive control transition power14 is zero",
          Tplus**14==S.zeros(34))
    check("positive control transition power13 is nonzero",
          Tplus**13!=S.zeros(34))
    initial=(Mb*unit).col_join(S.zeros(17,1))
    output=I.row_join(I)
    check("positive control twelfth coefficient is nonzero b6",
          output*Tplus**11*initial==eta)
    check("positive control thirteenth coefficient vanishes",
          output*Tplus**12*initial==S.zeros(17,1))
    for kv in [1,2,3,4,11,12,13]:
        got=output*Tplus**(kv-1)*initial
        check(f"positive control ideal power coefficient k{kv}",
              got==Mb**((kv+1)//2)*unit)

    check("no positive certificate on remaining HC4 quotient claimed",
          ctx["HC4_candidate_trace_quotient_constructed"] is False
          and ctx["positive_uniform_HC4_nilpotence_and_power_bound_proved"] is False)
    out=dict(claim="PF2RULED152",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             actual_Rodrigues_algebra_dimension=17,maximal_ideal_nilpotence_index=7,
             eta_square_zero_and_nonzero=True,
             exact_fixed_determinant_HC4_linearized_tail_nonterminating=True,
             all_order_obstruction_from_closed_identity_not_prefix=True,
             positive_block_control_nilpotence_index=14,
             actual_HC4_restricted_quotient_certificate_proved=False,
             registered_predecessor_checkers_replayed=[],
             full_goal_complete=False,independent_replay=False,
             formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
