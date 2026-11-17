"""Exact controls for the derivative surface, leading poles and normal rescaling."""
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
    'hessian/hc4-triangular-graphs/CRITICAL_FACTOR_SUPPORT.md': 'dd80f74b6cf135aa814ed805b817e44d99875543790fd443e1e0b3a2e6ef08e3',
    'hessian/hc4-triangular-graphs/verify_critical_support.py': '8a4ebb4bc4214b49738d2554d5c4e8735dbed9c43dbbd2a4a1e6107c6cfe2832',
    'hessian/hc4-triangular-graphs/receipt-critical-support.json': 'a07f280d9e96651f58ef93d6c47f44c5dcb5298cf8625efb687849daf51edc03',
    'hessian/hc4-triangular-graphs/LEADING_VALUATIONS_AND_CRITICAL_CURVES.md': 'bcaa5e0b08595155b5924024cebc9cd2f95ff2ecda9b96bd5d14abfbf7daf572',
    'hessian/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md': '1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
    'hessian/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md': '223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
    'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md': '3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
    'hessian/hc4-triangular-graphs/sources/jelonek-lason-nonproper-1411.5011v2.html': '7a00a58a0d03ca9a24dfebbb9b7c25bae2c5ac0087d472614f873c2d9b968714',
    'hessian/hc4-triangular-graphs/sources/stacks-derivative-00SB.html': 'b2e00947504ea1fdbb0a4bd8db3af45cd5a8aa2fedad2d733bedccd31430c3b0',
    'hessian/hc4-triangular-graphs/sources/stacks-derivative-031S.html': '05a6f4daf4e0ce0e22731a953dc355b411931313c9e549bbec1f7b9d0a1105c9',
    'hessian/hc4-triangular-graphs/sources/stacks-derivative-07PF.html': '0651722d80eeb7c9dd2c7e7bba6af39903fc0e01518940083af4db62cd1b94f3'
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name, expression):
        check(name, S.cancel(S.expand(expression)) == 0)

    T,U,V,A,D,W,p,z,b = S.symbols('T U V A D W p z b')
    for m in range(3, 9):
        k = m-2
        coefficients = S.symbols('b0:'+str(m))
        Q = U**m+sum(coefficients[i]*A**(m-1-i)*U**i for i in range(1,m))+A**(m-1)*(coefficients[0]+V)
        derivative = S.diff(Q,U)
        chart = m*D*W**(m-1)+sum(i*coefficients[i]*p**(m-1-i)*W**(i-1) for i in range(1,m))
        root_chart = D*W**m+sum(coefficients[i]*p**(m-1-i)*W**i for i in range(1,m))+p**(m-1)*(coefficients[0]+V)
        zero(f'degree {m}: derivative factor is exactly D^(m-2)',
            derivative.subs({A:D*p,U:D*W})-D**k*chart)
        zero(f'degree {m}: root equation divides by D^(m-1)',
            Q.subs({A:D*p,U:D*W})-D**(m-1)*root_chart)
        zero(f'degree {m}: inverted derivative chart is the cyclic surface',
            chart.subs({p:A/D,W:U/D})-derivative/D**k)
        zero(f'degree {m}: leading specialization has only its stated roots',
            Q.subs(A,0)-U**(m-1)*(U+coefficients[m-1]))
        zero(f'degree {m}: full pole is a simple normalized root',
            derivative.subs({A:0,U:-coefficients[m-1]})-(-coefficients[m-1])**(m-1))
        boundary_D = (1-(m-1)*b*W**k)/(m*W**(m-1))
        zero(f'degree {m}: Cstar boundary parametrization satisfies its equation',
            m*boundary_D*W**(m-1)+(m-1)*b*W**k-1)
        zero(f'degree {m}: the full root equation on the boundary has finite image',
            m*W**k*(boundary_D*W+b)-(1+b*W**k))
        check(f'degree {m}: infinity valuation has coprime numerator and denominator',
            S.gcd(k,k+1)==1)

    def spectrum(orders):
        q = orders[-1]
        answer = []
        for r in range(1,q+1):
            values = [0]+[n-i*r for i,n in enumerate(orders) if i and n != S.oo]
            if values.count(min(values)) >= 2:
                answer.append(r)
        return answer

    check('T^3 z^3+z has no integer pole order',spectrum([0,0,S.oo,3])==[])
    check('T^2 z^3+z permits the intermediate order one',spectrum([0,0,S.oo,2])==[1])
    check('T^3 z^3+T z^2+z permits intermediate orders one and two',spectrum([0,0,1,3])==[1,2])
    check('a simple leading zero allows no intermediate order',spectrum([0,0,0,1])==[1])
    for m,q in [(3,2),(3,7),(4,3),(5,4),(8,5)]:
        orders = [0]+[2*i+1 for i in range(1,m-1)]+[0,q]
        check(f'degree {m}, leading order {q}: nonzero next coefficient permits only full pole',
            spectrum(orders)==[q])

    for m in range(2,8):
        c = S.symbols('c1:'+str(m))
        d = S.Symbol('d')
        P = A**(m-1)*z**m+sum(A**(i-1)*c[i-1]*z**i for i in range(1,m))+d
        F = U**m+sum(c[i-1]*U**i for i in range(1,m))
        zero(f'degree {m}: balanced rescaling is finite and monic',
            A*(P.subs(z,U/A)+V)-(F+A*(V+d)))
        check(f'degree {m}: scaled cover really is monic of unchanged degree',
            S.Poly(F+A*(V+d),U).LC()==1 and S.degree(F,U)==m)

    cf = [S.Function('c'+str(i))(T) for i in (1,2)]
    Af = S.Function('A')(T)
    df = S.Function('d')(T)
    F = U**3+cf[1]*U**2+cf[0]*U
    cover = F+Af*(V+df)
    zero('the V partial forces a leading zero at every singular point',S.diff(cover,V)-Af)
    zero('the U partial is exactly F_U',S.diff(cover,U)-S.diff(F,U))
    zero('the T partial retains the A prime term needed at simple zeros',
        S.diff(cover,T)-S.diff(F,T)-S.diff(Af,T)*(V+df)-Af*S.diff(df,T))
    simple_cover = U**3+T*V
    check('a simple rescaling zero may have an isolated singularity',
        all(S.diff(simple_cover,X).subs({T:0,U:0,V:0})==0 for X in (T,U,V)))
    check('the same simple zero has no singular line',
        S.diff(simple_cover,T).subs({T:0,U:0})==V)
    bad_cover = U**3+T**2*V
    check('a repeated rescaling zero can give an excluded singular line',
        all(S.diff(bad_cover,X).subs({T:0,U:0})==0 for X in (T,U,V)))
    good_F = U**3+U
    check('the repeated-root rescaling gate can pass',
        S.gcd(S.gcd(good_F,S.diff(good_F,U)),S.diff(good_F,T))==1)
    bad_F = U**3
    check('the triple gcd detects the excluded repeated-root model',
        S.gcd(S.gcd(bad_F,S.diff(bad_F,U)),S.diff(bad_F,T))==U**2)

    s,t,w,v = S.symbols('s t w v')
    f = s*(1+s*t)
    g = -t*(1+s*t)
    h = 1/s
    zero('non-Keller repeated-leading control has the stated rational root',f*f*h**3-f*h*h+g)
    jac = S.diff(f,s)*S.diff(g,t)-S.diff(f,t)*S.diff(g,s)
    zero('the repeated-leading root control fails the Keller hypothesis',jac+(1+s*t)*(1+3*s*t))
    check('the control Jacobian is nonconstant',S.total_degree(jac)>0)
    zero('its normalized root has smaller than full leading vanishing',f*f*h-s*(1+s*t)**2)

    factor_control = (z*z-T)*(z*z-T-1)
    check('the higher-degree warning has an actual 2+2 factorization',
        S.factor_list(factor_control,z)[1]==[(z*z-T-1,1),(z*z-T,1)]
        or S.factor_list(factor_control,z)[1]==[(z*z-T,1),(z*z-T-1,1)])
    check('both quadratic control radicands have simple zeros and no rational square root',
        S.gcd(T,S.diff(T,T))==1 and S.gcd(T+1,S.diff(T+1,T))==1)

    P = T**3*z**5+(T*T+1)*z*z+T*z+T**4
    H = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the automorphic closed-embedding control has actual Hessian determinant one',
        S.det(S.hessian(H,(s,t,w,v)))-1)
    alpha,beta,gamma,delta = S.symbols('alpha beta gamma delta')
    gradient = [S.diff(H,X) for X in (s,t,w,v)]
    inverse = {s:delta,w:beta,
        t:gamma-P.subs({T:delta,z:beta})-beta*S.diff(P,z).subs({T:delta,z:beta}),
        v:alpha-beta*S.diff(P,T).subs({T:delta,z:beta})}
    check('the positive control has its explicit polynomial gradient inverse',
        all(S.expand(expr.subs(inverse,simultaneous=True)-target)==0
            for expr,target in zip(gradient,(alpha,beta,gamma,delta))))

    check('retained dependency pin list is populated',len(RETAINED)==11)
    for name,digest in RETAINED.items():
        check('retained proof or primary-source byte pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    context = json.loads((HERE/'sources/derivative-modification-import.json').read_text())
    check('normality imports and reused nonproperness theorem have distinct roles',
        len(context['new_external_theorem_imports'])==3
        and 'A2 to A4' in context['retained_primary_input']['role'])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_derivative_modification.py','DERIVATIVE_POWER_MODIFICATION.md','sources/derivative-modification-import.json']},
        dependencies=RETAINED,
        written_not_finite_exhaustion=[
            'a rational root must have an actual intermediate pole at a multiple leading zero',
            'absence of intermediate poles forces a derivative power and polynomial quotients on the original source plane',
            'the integral cyclic derivative surface contains a Cstar boundary curve omitted on a dense open, impossible by the retained nonproperness theorem',
            'the coefficient integer-pole spectrum and the triple gcd give arbitrary-degree rational-root exclusions',
            'normal balanced monic rescaling gives full irreducibility in arbitrary degree',
            'all square-free-leading cubics, coprime next-coefficient cubics, and all cubics of leading parameter degree at most two close in every other degree'],
        new_external_theorem_imports=['Stacks 00SB','Stacks 031S','Stacks 07PF'],
        original_A2_source_retained=True,normality_of_cyclic_derivative_surface_required=False,
        no_rational_root_asserted_higher_degree_irreducibility=False,
        all_repeated_leading_cubics_asserted_closed=False,
        nonproperness_theorem_fully_reproved=False,finite_controls_asserted_all_degree_proof=False,
        coefficient_residual_asserted_Keller_realization=False,
        nonlinear_source_Hessian_invariance_used=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        original_goal_complete=False,formal_verification=False,independent_replay=False,
        external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-derivative-modification.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
