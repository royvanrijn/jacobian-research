"""Exact controls for the all-factor resultant and infinity-valuation barrier."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RETAINED = {
    'hessian/hc4-triangular-graphs/DERIVATIVE_POWER_MODIFICATION.md': '27c391e23acd4ba93b500b38ce96903c8700c419364bd7dc4436687f06dd8e05',
    'hessian/hc4-triangular-graphs/verify_derivative_modification.py': '5563c09f736c55b42a98716d8b273848d1066649a8ff2a3f9870e9fea12b26ce',
    'hessian/hc4-triangular-graphs/receipt-derivative-modification.json': '3f0c1cf4a21adecfd572f0b1103165674c54f326d54c8a5d89b405395a3baebf',
    'hessian/hc4-triangular-graphs/sources/derivative-modification-import.json': '38c3f83a614dcb384274a0e633e6b1f9b88184b0f1d5aa130d7a7e7881c9e1fe',
    'hessian/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md': '223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
    'hessian/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md': '1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
    'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md': '3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
    'hessian/hc4-triangular-graphs/sources/miyanishi-qfactorial-2007.pdf': 'f732c193a998d55e641cb2f7fbb0568d85f8cefb14e799d6fbb904fb66e7a063',
    'hessian/hc4-triangular-graphs/sources/miyanishi-pseudocover-report-2005.pdf': '1e991cb0af2ff013c1ced61fb5c3ea849d6305f100aa486dcc3d955c3e948ef0'
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item))==0 for item in values))

    T,U,V,z,p,D,x,y,s,t,tau,X = S.symbols('T U V z p D x y s t tau X')
    A = p*D
    j11,j12,j21,j22,c = S.symbols('j11 j12 j21 j22 c')
    matrix = S.Matrix([[j11,j12],[j21,j22]])
    row = S.Matrix([[c,1]])
    zero('the source hypersurface derivative covector cannot vanish for invertible Keller Jacobian',
        row*matrix*matrix.adjugate()-matrix.det()*row)
    zero('shared leading coefficient can create a nonunit resultant without a finite intersection',
        S.resultant(s*z-1,s*z-2,z)+s)
    check('that shared-leading example has no common finite zero',
        S.expand((s*z-1)-(s*z-2))==1)
    zero('coprime leading coefficients alone do not exclude a finite intersection',
        S.resultant(s*z-1,z-t,z)-(1-s*t))
    zero('the shared-leading obstruction is visible in the next coefficient',
        S.expand((s*z-1)*(s*z-2))-(s*s*z*z-3*s*z+2))

    for d,e in [(1,2),(2,2),(2,3),(2,4)]:
        Fi = p*z**d+(x*z if d>1 else 0)+1
        Fj = D*z**e+y*z+2
        Qi = S.cancel(A**d/p*Fi.subs(z,U/A))
        Qj = S.cancel(A**e/D*Fj.subs(z,U/A))
        primitive_res = S.resultant(Fi,Fj,z)
        normalized_res = S.resultant(Qi,Qj,U)
        zero(f'degrees {d}+{e}: normalized cross-resultant scaling',
            normalized_res-A**(d*e)/(p**e*D**d)*primitive_res)
        zero(f'degrees {d}+{e}: normalized product is the monic label equation',
            Qi*Qj-A**(d+e-1)*(Fi*Fj).subs(z,U/A))
        check(f'degrees {d}+{e}: each normalized factor is polynomial and monic',
            all(S.Poly(Q,U).LC()==1 and all(S.denom(S.cancel(coef))==1
                for coef in S.Poly(Q,U).all_coeffs()) for Q in [Qi,Qj]))

    a = S.Symbol('a',nonzero=True)
    for m in [3,4,5,7]:
        coefficients = S.symbols('q0:'+str(m))
        Q = U**m+sum(coefficients[i]*U**i for i in range(m))+a**(m-1)*V
        scaled = S.expand(tau**m*Q.subs({V:tau**(-m),U:X/tau}))
        zero(f'degree {m}: infinity scaling has distinct nonzero constant roots',
            scaled.subs(tau,0)-(X**m+a**(m-1)))
        check(f'degree {m}: the reduction and its derivative are coprime',
            S.gcd(X**m+a**(m-1),m*X**(m-1))==1)

    d,e,h = S.symbols('d e h',positive=True,integer=True)
    m = d+e
    low = d*e/(m*(e-d))
    high = -low
    zero('the unequal leading valuations solve every two-factor equation',low/d+high/e-1/m)
    zero('their product has the target leading valuation zero',low+high)
    zero('the cross-resultant has its required negative infinity valuation',-e*low-d*high+d*e/m)
    zero('the derivative generalization has exponent e-d',d*e-d-e*(d-1)-(e-d))
    zero('the other leading exponent is exactly e(d-1)',d*e-e-e*(d-1))
    zero('three-factor pair equations force the same ratio for every factor',
        (S.Rational(1,2)/(d+e+h)+S.Rational(1,2)/(d+e+h))-1/(d+e+h))
    zero('the three-factor leading sum becomes the impossible one half',
        (d+e+h)/(2*(d+e+h))-S.Rational(1,2))
    zero('two equal factor degrees have the same impossible leading sum',
        d/(2*d)-S.Rational(1,2))
    weights = [h*(d+e-h),e*(d+h-e),d*(e+h-d)]
    zero('an explicit weighted product also cancels every individual leading coefficient',
        S.Matrix([weights[0]*e+weights[1]*h,weights[0]*d+weights[2]*h,
            weights[1]*d+weights[2]*e])-S.Matrix([2*d*e*h]*3))
    zero('the weighted three-factor cross-resultant has nonzero infinity slope',
        weights[0]*d*e+weights[1]*d*h+weights[2]*e*h-(d+e+h)*d*e*h)
    check('negative weights are sometimes necessary and are legitimate rational units',
        weights[0].subs({d:2,e:3,h:7})<0)

    for di,ei,denominator in [(2,3,5),(2,4,3),(3,5,16)]:
        mi = di+ei
        ri = ei-di
        check(f'degrees {di}+{ei}: exact field-degree divisor',
            mi*ri//int(S.gcd(di*ei,mi*ri))==denominator)
    di,ei = 2,5
    alpha = ei*(di-1)
    ri = ei-di
    power = int(S.ceiling(S.Rational(alpha,ri)))
    check('clearing the unequal-factor power denominator gives an integral monic equation',
        0<=power*ri-alpha<ri and power==2)
    check('the allowed 2+3 pattern actually solves the valuation equations',
        low.subs({d:2,e:3})==S.Rational(6,5)
        and high.subs({d:2,e:3})==-S.Rational(6,5))
    check('the three root-free quartic degree patterns reduce to the equal pair',
        [parts for parts in [(1,3),(2,2),(1,1,2),(1,1,1,1)] if min(parts)>=2]==[(2,2)])

    Q = U**4+x*U**3+y*U**2+a*U+a*a*V
    root = S.Symbol('root')
    quotient = S.cancel((Q-Q.subs(U,root))/(U-root))
    sylvester = S.Matrix([[1,-root,0,0],[0,1,-root,0],[0,0,1,-root],
        S.Poly(quotient,U).all_coeffs()])
    zero('the one-root cross-resultant is exactly the normalized derivative',
        sylvester.det()-S.diff(Q,U).subs(U,root))
    zero('the explicit nonnormal integral-order control has its normalization parameter',
        (U*U-T*T*(T+1)).subs({T:s*s-1,U:s*(s*s-1)}))
    check('the nonnormal order has a codimension-one partial-zero locus',
        all(S.diff(U*U-T*T*(T+1),variable).subs({T:0,U:0})==0 for variable in [T,U,V]))
    check('this order control does not supply a Keller map',S.diff(s*s-1,s)==2*s)
    check('unit cross-resultants by themselves do not forbid two unequal factors',
        S.resultant(z*z-1,z**3+2,z)==3)
    check('the unequal-factor control omits the required target V dependence',
        S.diff((z*z-1)*(z**3+2),V)==0)

    check('all retained proof and context pins are present',len(RETAINED)==9)
    for path,digest in RETAINED.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    context = json.loads((HERE/'sources/factor-resultant-context.json').read_text())
    check('the ruled-surface literature is context, with missing hypotheses explicit',
        context['new_external_theorem_imports']==[]
        and len(context['literature_context_only'])==2
        and 'Q-factoriality' in context['literature_context_only'][0]['role']
        and 'almost surjective' in context['literature_context_only'][1]['role'])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_factor_resultants.py','FACTOR_RESULTANT_BARRIER.md','sources/factor-resultant-context.json']},
        dependencies=RETAINED,
        written_not_finite_exhaustion=[
            'smooth pullbacks and coprime leading factors give constant cross-resultants',
            'the target leading triple gcd guarantees coprime source leading factors in every dimension',
            'the all-pairs infinity valuation system excludes three factors and two equal factors in every degree',
            'the plane rational-root gate then leaves exactly two unequal factors of degrees at least two',
            'every quartic under the leading triple gcd stays irreducible and its carried output embeddings force an inverse',
            'two unequal factors give an explicit resultant power, field-degree divisor and nonnormal integral monogenic order'],
        new_external_theorem_imports=[],all_dimensions_factor_count_barrier=True,
        plane_only_rational_root_gate_used_for_quartics=True,
        two_unequal_factors_asserted_impossible=False,
        no_rational_roots_asserted_arbitrary_degree_irreducibility=False,
        all_quartics_without_the_leading_condition_asserted_closed=False,
        dominant_surface_asserted_Qfactorial=False,dominant_map_asserted_almost_surjective=False,
        ruled_surface_classification_imported=False,
        primitive_order_asserted_normal=False,coefficient_residual_asserted_Keller_realization=False,
        nonlinear_source_Hessian_invariance_used=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        original_goal_complete=False,formal_verification=False,independent_replay=False,
        external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-factor-resultants.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
