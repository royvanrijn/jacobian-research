"""Exact controls for the all-degree Newton and critical-ramification gates.

Formal lifting, arbitrary-degree pole orders and the polynomial derivative
lemma are written proofs. These controls do not replace those arguments.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,expression):
        values = list(expression) if isinstance(expression,S.MatrixBase) else [expression]
        check(name,all(S.cancel(value.doit())==0 for value in values))

    x,p,z,t,U,Z,V,Q = S.symbols('x p z t U Z V Q')
    L,N,c = S.symbols('L N c',nonzero=True)
    M = S.symbols('M')
    # The empty-scheme example has a repeated whole curve, rather than
    # a spurious second component at its repeated leading roots.
    R = 3*p*p+S.Rational(3,2)*x*p**4
    F = 3+2*x*R
    zero('entirely repeated critical curve exact factorization',F-3*(1+x*p*p)**2)
    zero('relative derivative vanishes along every repeated factor',S.diff(R,p)-6*p*(1+x*p*p))
    zero('no-simple-root boundary of the empty-scheme control',R.subs(x,0)-3*p*p)
    zero('boundary p derivative cannot vanish in the empty-scheme control with M=1',
         (2*L*p+1+N*S.diff(R.subs(x,0),p)*Q).subs(p,0)-1)
    older = 3*p+S.Rational(3,2)*x*p*p
    zero('earlier repeated curve retains its simple boundary root',3+2*x*older-3*(1+x*p)**2)
    check('empty-scheme rule does not incorrectly include the earlier simple-boundary control',
          S.diff(older.subs(x,0),p)==3)
    # Rank-one nonemptiness and exact reconstruction, including a curved
    # Keller coordinate C whose zero fibre really is nonempty.
    A = U+Z*Z
    C = Z-A**3
    D0 = U**3+U*Z
    J = S.Matrix([A,C]).jacobian((U,Z)).det()
    zero('rank-one control residual pair has unit Jacobian',J-1)
    ss = (-S.diff(D0,U)*S.diff(C,Z)+S.diff(D0,Z)*S.diff(C,U))/J
    vv = (-S.diff(A,U)*S.diff(D0,Z)+S.diff(A,Z)*S.diff(D0,U))/J
    h = D0+x*A+V*C
    zero('rank-one critical equations reconstruct the entire nonempty plane zero curve',
         S.Matrix([S.diff(h,U),S.diff(h,Z)]).subs({x:ss,V:vv},simultaneous=True))
    zero('rank-one reconstructed component includes a literal plane point',C.subs({U:0,Z:0}))

    # Universal coefficients of the critical-Hessian leading term.
    universal_D = -(M-13*N/(2*x*x))**2/(4*L)
    zero('universal order-six coefficient of D second derivative',
         S.limit(x**6*S.diff(universal_D,x,2),x,0)+845*N*N/(4*L))
    Pprime,zeta,lam = S.symbols('Pprime zeta lam',nonzero=True)
    zero('simple Newton branch leading critical-Hessian coefficient is 1690',
         -32*L*N*N*(Pprime/2)**2*(-845*N*N/(4*L))-1690*N**4*Pprime**2)
    zero('large-slope leading coefficient is nonzero symbolically',
         -32*L*N*N*(Pprime/2)**2*(-N*zeta*(lam+2)*(lam+3)/4)
         -2*L*N**3*zeta*Pprime**2*(lam+2)*(lam+3))
    # Formal lifting on the slope-one, degree-two separated control.
    # Higher coefficient terms are retained, not only its Newton binomial.
    a = x+2*x*x
    P = p*p+3*p+5
    b = 7+x+x**3
    quadratic = a*P+b
    # F after p=z/x has simple initial face 3+2z^2.
    zero('separated quadratic initial face is the asserted binomial',
         S.limit((3+2*x*quadratic).subs(p,z/x),x,0)-(3+2*z*z))
    disc = (3*a)**2-4*a*(5*a+b+3/(2*x))
    on_curve = S.rem(S.diff(quadratic,p)**2-disc,quadratic+3/(2*x),p)
    zero('quadratic derivative squares equal the same exact rational discriminant',on_curve)
    # Leading root difference c/x contributes an unavoidable order-five
    # difference even after any common universal critical-value term.
    zero('quadratic branch-difference second derivative coefficient',
         S.diff((c/x)/x**2,x,2)-12*c/x**5)
    # Elementary all-degree derivative lemma: the identities used in its
    # divisibility proof, with exact polynomial controls for distinct degrees.
    PP = S.Function('PP')(z)
    zero('derivative lemma final identity forces PP third derivative to vanish',
         S.diff(S.diff(PP,z)**2-2*PP*S.diff(PP,z,2),z)+2*PP*S.diff(PP,z,3))
    face2 = 3+2*z*z
    zero('quadratic face can match the necessary squared derivative values',
         S.rem(S.diff(face2,z)**2,face2,z)+24)
    for degree in (3,5,8):
        face = 3+2*z**degree
        remainder = S.rem(S.diff(face,z)**2,face,z)
        check('squarefree binomial degree '+str(degree)+' has nonconstant squared-derivative remainder',
              S.degree(remainder,z)>0)
    # All monic p degrees >=2 are excluded in the proof. A nontrivial
    # degree-seven control also retains lower coefficients and base degree.
    monic = p**7+(x*x+1)*p**3+x**12*p+9
    zero('monic degree-seven formal face after x=t^7,p=z/t',
         (3+2*x*monic).subs({x:t**7,p:z/t},simultaneous=True).expand().subs(t,0)-(3+2*z**7))

    # Finite smooth ramification: R=p^e+1 gives x=-3/(2(1+p^e)).
    # Every control checks the actual critical Hessian, not a proxy Hessian.
    for degree in (2,3,5):
        xc = -S.Rational(3,2)/(1+t**degree)
        pc = t
        Rp = degree*t**(degree-1)
        Dc = -N*pc/(4*xc**2)-(M-13*N/(2*xc**2))**2/(4*L)
        delta_D = S.diff(Dc,t)/S.diff(xc,t)
        delta2_D = S.diff(delta_D,t)/S.diff(xc,t)
        Hcrit = -32*L*N*N*xc**6*Rp**2*delta2_D
        zero('finite ramification index '+str(degree)+' has the exact nonzero simple-pole coefficient',
             S.limit(t*Hcrit,t,0)-18*L*N**3*(1-degree))
    mixed = 3*p*p+S.Rational(3,2)*x*p**4+S.Rational(3,2)*x*p*p
    mixed_F = 3+2*x*mixed
    point = {x:-4,p:S.I/2}
    zero('repeated-face splitting control has a finite smooth ramification point',
         S.Matrix([mixed_F,S.diff(mixed_F,p)]).subs(point))
    zero('finite control has nonzero F_x=3',S.diff(mixed_F,x).subs(point)-3)
    zero('finite control has ramification index two via F_pp=-96',S.diff(mixed_F,p,2).subs(point)+96)
    zero('finite control still has the entirely repeated initial face',
         (mixed_F.subs({x:t*t,p:z/t},simultaneous=True)).expand().subs(t,0)-3*(1+z*z)**2)
    vertical = (x-1)*(p*p+1)
    zero('vertical-component caveat is real: F_p vanishes identically on that component',
         S.diff(vertical,p).subs(x,1))
    check('vertical germ is explicitly not a finite ramification control',S.diff(vertical,x).subs({x:1,p:0})!=0)

    dependencies = {}
    for receipt_name,proof_name in [('receipt-relative-critical.json','RELATIVE_CRITICAL_LOCUS.md'),
                                    ('receipt-critical-components.json','CRITICAL_COMPONENTS.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in [proof_name]:
            digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            check('retained prior written proof byte pin '+name,digest==retained['source_sha256'][name])
            dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input for empty-scheme theorem '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'NEWTON_CRITICAL_GATE.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['two projection ranks for empty relative critical schemes',
            'formal implicit lifting for any simple initial Newton root',
            'all rational Newton slopes and exact common critical-Hessian limits',
            'squarefree polynomial derivative lemma in every degree',
            'all separated coefficient degrees including the quadratic exception',
            'every finite nonvertical smooth ramification index'],
        every_separated_nonconstant_p_graph_forbids_unit_Hessian=True,
        separated_and_Newton_proofs_use_Tang=False,
        empty_scheme_Euler_proof_uses_Tang=True,
        every_monic_in_p_degree_at_least_two_graph_forbids_unit_Hessian=True,
        repeated_initial_face_alone_asserted_closed=False,
        vertical_components_mislabeled_finite_ramification=False,
        every_nonlinear_R_graph_closed=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
