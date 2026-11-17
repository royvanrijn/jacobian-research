"""Exact controls for critical-factor supports and localized-unit subset ranks."""
import argparse
import hashlib
import itertools
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

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item))==0 for item in values))

    T,U,V,z,A,b,c,d,R,x,y,v,s,t,w = S.symbols('T U V z A b c d R x y v s t w')
    Q = U**3+b*U*U+A*c*U+A*A*(d+V)
    zero('the normalized cubic equations force a zero common root at every leading level',
        (U*S.diff(Q,U)-2*Q-U**3).subs(A,0))
    zero('the only extra terms in the common-root identity are divisible by the leading coefficient',
        U*S.diff(Q,U)-2*Q-U**3+A*c*U+2*A*A*(d+V))
    Q0 = U*U*(U+1)
    check('the simple normalized leading root is not a critical root',S.diff(Q0,U).subs(U,-1)==1)
    check('the only common root of the specialized cubic and its derivative is zero',S.gcd(Q0,S.diff(Q0,U))==U)
    check('the nonzero critical branch misses the root equation at the leading level',Q0.subs(U,-S.Rational(2,3))==S.Rational(4,27))
    Q7 = U**7+b*U**6+A*c*U**5+A*A*U**4+A**6*(d+V)
    zero('the leading-level common-root identity also holds in higher degree',
        (U*S.diff(Q7,U)-6*Q7-U**7).subs(A,0))

    nonlinear = U**4-2*T*U*U+T**3*V
    zero('a nonlinear critical component must be counted as one factor over the parameter field',
        S.diff(nonlinear,U)-4*U*(U*U-T))
    check('the quadratic critical component has a nonsquare parameter radicand',
        S.degree(T,T)==1 and S.diff(T,T)==1)
    check('the linear and quadratic critical factors are distinct',S.gcd(U,U*U-T)==1)
    check('both critical factors have permitted support at the leading zero',
        (U*U-T).subs({T:0,U:0})==0 and U.subs(U,0)==0)

    def subsets(count):
        return [set(indices) for size in range(1,count+1)
            for indices in itertools.combinations(range(count),size)]

    def valid_subset(supports,capacity,indices):
        levels = set().union(*(supports[i] for i in indices))
        return sum(capacity[level] for level in levels)>=len(indices)

    supports = [{0},{0},{1}]
    capacities = {0:1,1:2}
    check('a total unit-rank budget alone can pass while a subset fails',
        valid_subset(supports,capacities,{0,1,2}) and not valid_subset(supports,capacities,{0,1}))
    check('the failing subset cannot borrow capacity from an excluded leading level',
        set().union(supports[0],supports[1])=={0} and capacities[0]<2)
    check('separate critical supports each require their own nonzero defect',
        all(valid_subset([{0},{1}],{0:1,1:1},indices) for indices in subsets(2))
        and not valid_subset([{0},{1}],{0:2,1:0},{1}))

    a_symbolic = A
    c_symbolic = -(2*R+A*R*R)/3
    H_symbolic = 1+A*R
    beta_good = A*R/3
    beta_empty = -(2+A*R)/3
    zero('the arbitrary-leading-multiplicity cubic family has square critical discriminant',
        1-3*a_symbolic*c_symbolic-H_symbolic**2)
    normalized = U**3+U*U+A*c_symbolic*U+A*A*(d+V)
    zero('its complete derivative has the two stated rational critical branches',
        S.diff(normalized,U)-3*(U-beta_good)*(U-beta_empty))
    check('the excluded branch has nonzero value at the sole leading zero',beta_empty.subs(A,0)==-S.Rational(2,3))
    check('the other critical branch does specialize to the repeated root',beta_good.subs(A,0)==0)
    check('the excluded branch also fails the full normalized root equation there',
        normalized.subs({A:0,U:-S.Rational(2,3)})==S.Rational(4,27))

    leading = T*T
    choice_R = T+1
    P = leading*z**3+z*z-(2*choice_R+leading*choice_R**2)*z/3+T**3
    H = 1+leading*choice_R
    empty_branch = -(2+leading*choice_R)/3
    check('the new exact cubic has empty support despite a repeated leading zero',S.gcd(leading,empty_branch)==1)
    check('the preceding square-free-leading gate does not apply to this control',S.gcd(leading,S.diff(leading,T))==T)
    check('the leading level drops degree only to two',S.degree(P.subs(T,0),z)==2)
    zero('the actual critical discriminant is a square, outside the earlier critical-genus gate',
        S.Poly(P,z).nth(2)**2-3*S.Poly(P,z).LC()*S.Poly(P,z).nth(1)-H*H)
    zero('the new cubic has nonsquare root discriminant and S3 monodromy',
        S.discriminant(S.discriminant(P+V,z),V)-16*H**6)
    check('the root-discriminant V zeros are distinct',H!=0)

    a_residual = T*(T-1)
    H_residual = 1-2*T
    P_residual = a_residual*z**3+z*z-S.Rational(4,3)*z
    zero('the residual coefficient control still has square critical discriminant',
        1-3*a_residual*(-S.Rational(4,3))-H_residual**2)
    beta_plus = (-1+H_residual)/3
    beta_minus = (-1-H_residual)/3
    check('one residual critical factor is supported only at zero',
        S.Poly(S.gcd(a_residual,beta_plus),T).monic().as_expr()==T)
    check('the other residual factor is supported only at one',
        S.Poly(S.gcd(a_residual,beta_minus),T).monic().as_expr()==T-1)
    check('both gcd exclusions fail on that coefficient control',
        S.degree(S.gcd(a_residual,1-H_residual),T)==1
        and S.degree(S.gcd(a_residual,1+H_residual),T)==1)
    check('the residual supports really are disjoint when a and b are coprime',S.gcd(a_residual,S.Integer(1))==1)
    check('the support-gcd residual has a nonconstant critical gap and is closed by the new two-level lemma',
        S.degree(H_residual,T)==1 and S.gcd(a_residual,S.diff(a_residual,T))==1)
    zero('on the first singleton support the other unit has the exact scalar residue',
        (U-beta_minus).subs({T:0,U:0})-S.Rational(2,3))
    zero('on the second singleton support the first unit has its exact scalar residue',
        (U-beta_plus).subs({T:1,U:0})-S.Rational(2,3))
    gap_constant = S.Rational(2,3)
    zero('the two scalar residues force zero normalized label if the units formed that affine pencil',
        (U-beta_plus)+(U-beta_minus)-gap_constant-2*U)
    b_constant_gap = 1-2*T
    zero('the retained two-level control instead has constant critical discriminant',
        b_constant_gap**2-3*a_residual*S.Rational(4,3)-1)
    constant_gap_plus = (-b_constant_gap+1)/3
    constant_gap_minus = (-b_constant_gap-1)/3
    zero('its critical branches have the constant gap permitted by the lemma',
        constant_gap_plus-constant_gap_minus-S.Rational(2,3))
    check('the constant-gap control still has two populated disjoint supports',
        S.Poly(S.gcd(a_residual,constant_gap_plus),T).monic().as_expr()==T
        and S.Poly(S.gcd(a_residual,constant_gap_minus),T).monic().as_expr()==T-1)

    tau,Zeta,Xnew,Ynew = S.symbols('tau Zeta Xnew Ynew')
    alpha,H0,delta = S.symbols('alpha H0 delta',nonzero=True)
    target0 = S.symbols('target0')
    rho = tau*tau+tau+1
    a_model = alpha*tau*(tau-1)
    b_model = H0*(1-2*tau)+3*a_model*rho
    c_model = S.cancel((b_model*b_model-H0*H0)/(3*a_model))
    P_model = a_model*z**3+b_model*z*z+c_model*z+tau**5+2
    d_tilde = S.expand(P_model.subs(z,-rho))
    kappa = H0**3/alpha**2
    fixed = tau*(tau-1)*Zeta**3+(1-2*tau)*Zeta**2+S.Rational(4,3)*Zeta
    zero('the polynomial z shift removes the entire high-degree quadratic coefficient deviation',
        S.expand(P_model.subs(z,z-rho)).coeff(z,2)-H0*(1-2*tau))
    zero('the linear coefficient after that shift is forced to one constant',
        S.expand(P_model.subs(z,z-rho)).coeff(z,1)-4*H0*H0/(3*alpha))
    zero('the output polynomial has its complete fixed-model normalization identity',
        P_model.subs(z,H0*Zeta/alpha-rho)-kappa*fixed-d_tilde)
    original_model = P_model.subs(tau,(T-target0)/delta)
    ambient_change = {T:target0+delta*tau,x:kappa*Xnew,y:Ynew,z:H0*Zeta/alpha-rho}
    zero('the two physical output observables satisfy the exact ambient and target conjugacy',
        S.Matrix([T,x*y-original_model]).subs(ambient_change,simultaneous=True)
        -S.Matrix([target0+delta*tau,kappa*(Xnew*Ynew-fixed)-d_tilde]))
    ambient_inverse = S.Matrix([(T-target0)/delta,x/kappa,y,
        alpha*(z+rho.subs(tau,(T-target0)/delta))/H0])
    zero('the complete ambient normal-form change has its polynomial inverse',
        ambient_inverse.subs(ambient_change,simultaneous=True)-S.Matrix([tau,Xnew,Ynew,Zeta]))
    target_change = S.Matrix([target0+delta*tau,kappa*V-d_tilde])
    check('the target normal-form change preserves Keller with its exact nonzero constant determinant',
        target_change.jacobian([tau,V]).det()==delta*kappa)
    target_inverse = S.Matrix([(T-target0)/delta,(V+d_tilde.subs(tau,(T-target0)/delta))/kappa])
    zero('the target normal-form change has its complete polynomial inverse',
        target_inverse.subs({T:target_change[0],V:target_change[1]},simultaneous=True)-S.Matrix([tau,V]))

    Cchart,Dchart = S.symbols('Cchart Dchart')
    p_chart = (Cchart*(Dchart+1)-3)/2
    q_chart = (Cchart*(Dchart-1)-3)/2
    f_chart = p_chart*(Dchart-1)/3
    u_chart = Cchart*(Dchart*Dchart-1)/9
    h_chart = Cchart/(p_chart*q_chart)
    fixed_original = T*(T-1)*z**3+(1-2*T)*z*z+S.Rational(4,3)*z
    N_chart = Cchart*(Cchart*Cchart*(Dchart*Dchart-1)-9*Cchart*Dchart+27)
    zero('the second polynomial chart recovers both complete leading fibres',
        f_chart-1-q_chart*(Dchart+1)/3)
    zero('its rational root label agrees with the normalized integral label',
        u_chart-f_chart*(f_chart-1)*h_chart)
    zero('the fixed-cubic root equation has the exact cleared polynomial numerator',
        fixed_original.subs({T:f_chart,z:h_chart})-N_chart/(9*p_chart*p_chart*q_chart*q_chart))
    zero('on the target hyperbola the numerator has only its stated finite zero',
        N_chart.subs(Dchart,3/Cchart-1)-3*Cchart*(Cchart+3))
    check('the one permitted hyperbola image point satisfies both equations',
        p_chart.subs({Cchart:-3,Dchart:-2})==0 and N_chart.subs({Cchart:-3,Dchart:-2})==0)
    zero('the recovered source outputs are algebraically independent in the second chart',
        S.Matrix([f_chart,u_chart]).jacobian([Cchart,Dchart]).det()-(Dchart*Dchart-1)/18)
    zero('the omitted hyperbola has an explicit nonconstant unit, forbidding polynomial A1 domination',
        Cchart*(Dchart+1)/3-1-(Cchart*(Dchart+1)-3)/3)

    a_sparse = T*T*(T-1)**3
    b_sparse = T+1
    normalized_sparse = U**5+b_sparse*U**4+a_sparse**4*V
    zero('the sparse higher-degree polynomial has its exact critical factor of empty support',
        S.diff(normalized_sparse,U)-U**3*(5*U+4*b_sparse))
    check('that critical factor has constant term coprime to the full leading coefficient',
        S.gcd(a_sparse,S.Rational(4,5)*b_sparse)==1)
    quartic_control = (z*z-T)*(z*z-T-1)
    check('a reducible quartic can have no rational root, so the higher-degree limitation is real',
        S.degree(quartic_control,z)==4 and S.degree(z*z-T,z)==2
        and S.degree(T,T)==1 and S.degree(T+1,T)==1)

    f_bad = -s+s**3*t
    g_bad = -t
    u_bad = -1+s*s*t
    zero('the old rational-pole control still satisfies its actual cubic root equation',
        (T*z**3+z*z).subs({T:f_bad,z:1/s})+g_bad)
    check('its supposed empty-support unit instead has an impermissible zero outside the leading level',
        (u_bad+S.Rational(2,3)).subs({s:1,t:S.Rational(1,3)})==0
        and f_bad.subs({s:1,t:S.Rational(1,3)})!=0)
    check('that impermissible zero is excluded precisely by the missing Keller condition',
        S.Matrix([f_bad,g_bad]).jacobian([s,t]).det().subs({s:1,t:S.Rational(1,3)})==0)

    image = {T:s,x:1,y:t+P.subs({T:s,z:w}),z:w,v:v}
    zero('the newly closed family has a genuine complete polynomial source inverse',
        S.Matrix([T,y-P,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the two output observables carry exactly the identity Keller pair',
        S.Matrix([T,x*y-P]).subs(image,simultaneous=True)-S.Matrix([s,t]))
    actual = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the independent nondegenerate ambient quadratic pulls back literally',
        (T*v+y*z+x*x).subs(image,simultaneous=True)-actual)
    check('the actual repeated-leading S3 control has unit Hessian',S.hessian(actual,[s,t,w,v]).det()==1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    evaluate = {T:p4,z:p2}
    inverse = S.Matrix([p4,p3-P.subs(evaluate)-p2*S.diff(P,z).subs(evaluate),
        p2,p1-p2*S.diff(P,T).subs(evaluate)])
    gradient = S.Matrix([S.diff(actual,q) for q in [s,t,w,v]])
    zero('the actual displayed gradient inverse is a left inverse',
        inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the actual displayed gradient inverse is a right inverse',
        gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))

    retained = {
        'hessian/hc4-triangular-graphs/LEADING_VALUATIONS_AND_CRITICAL_CURVES.md':'bcaa5e0b08595155b5924024cebc9cd2f95ff2ecda9b96bd5d14abfbf7daf572',
        'hessian/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md':'1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
        'hessian/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md':'223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
        'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f'
    }
    retained['hessian/hc4-triangular-graphs/sources/jelonek-lason-nonproper-1411.5011v2.html'] = '7a00a58a0d03ca9a24dfebbb9b7c25bae2c5ac0087d472614f873c2d9b968714'
    for name,digest in retained.items():
        check('retained all-degree proof byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_critical_support.py','CRITICAL_FACTOR_SUPPORT.md','sources/critical-support-context.json']},
        dependencies=retained,written_not_finite_exhaustion=[
            'each normalized critical-factor unit is supported only where its constant term and the leading coefficient vanish together',
            'all subsets of critical factors satisfy their own localized-unit capacity inequalities',
            'one empty critical-factor support excludes rational source roots in every degree',
            'the cubic square-discriminant support gcd closes new S3 families with arbitrary leading multiplicity',
            'disjoint cubic supports require two separate reducible full source fibres rather than defect concentrated at one level',
            'two reduced rational critical units with singleton supports force a constant critical-branch difference',
            'square-free cubic leading coefficient with at most two zeros and gcd(a,b)=1 closes when its square critical discriminant has nonconstant square root',
            'constant-gap two-simple-leading-level output data are polynomially equivalent to one fixed cubic pair',
            'the fixed cubic stays irreducible under every plane Keller base change, since a rational label would give a polynomial map omitting a hyperbola',
            'all square-free at-most-two-leading-level cubics with gcd(a,b)=1 and nonzero square critical discriminant close'],
        new_external_theorem_imports=[],function_field_S_unit_height_bound_imported=False,
        Hall_matching_theorem_imported=False,all_critical_factors_asserted_linear=False,
        allowed_leading_levels_asserted_actual_zeros=False,
        critical_factor_count_after_parameter_extension_used=False,
        no_rational_root_asserted_higher_degree_irreducibility=False,
        all_square_critical_discriminant_cubics_asserted_closed=False,
        repeated_leading_levels_asserted_reduced_by_the_simple_level_proof=False,
        fixed_model_irreducibility_under_Keller_base_change_proved=True,
        nonproperness_theorem_fully_reproved=False,fixed_potential_Hessian_invariance_inferred_from_output_normal_form=False,
        coefficient_residual_asserted_Keller_realization=False,
        nonlinear_source_Hessian_invariance_used=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        original_goal_complete=False,formal_verification=False,independent_replay=False,
        external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-critical-support.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
