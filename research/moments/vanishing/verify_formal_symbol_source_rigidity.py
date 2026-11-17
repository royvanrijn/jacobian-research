#!/usr/bin/env python3
"""Exact bounded transcription controls; universal rigidity is a written proof."""
from pathlib import Path
import hashlib
import json
import sympy as S

HERE = Path(__file__).resolve().parent
x, y, z, X, Y, Z, h, a, b, c, t = S.symbols('x y z X Y Z h a b c t')
variables = (x, y, z)
symbols = (X, Y, Z)


def action(symbol, polynomial):
    symbol = S.Poly(S.expand(symbol), *symbols)
    polynomial = S.Poly(S.expand(polynomial), *variables)
    assert not symbol.as_expr().has(S.Float) and not polynomial.as_expr().has(S.Float)
    assert polynomial.total_degree() <= 16, 'bounded control only'
    result = 0
    for powers, coefficient in symbol.terms():
        term = polynomial.as_expr()
        for variable, power in zip(variables, powers):
            term = S.diff(term, variable, power)
        result += coefficient * term
    return S.expand(result)


def iterate(symbol, polynomial, m):
    assert isinstance(m, int) and 0 <= m <= 4
    for _ in range(m):
        polynomial = action(symbol, polynomial)
    return polynomial


def gauge_coefficient(symbol, phase, m):
    # Independent finite Taylor contraction of lambda(s+partial_x)^m P^m.
    power = S.Poly(S.expand(phase**m), *variables)
    value = 0
    for powers, coefficient in power.terms():
        monomial = coefficient * S.prod(v**e for v, e in zip(variables, powers))
        # Differentiate each polynomial monomial by every supported beta.
        import itertools
        for beta in itertools.product(*(range(e + 1) for e in powers)):
            derivative = S.expand(symbol**m)
            output = monomial
            for v, q, k in zip(variables, symbols, beta):
                output = S.diff(output, v, k)
                derivative = S.diff(derivative, q, k) / S.factorial(k)
            value += derivative * output
    return S.expand(value)


def direct_gauge_coefficient(symbol, phase, m):
    # Exponential conjugation computed by binomially splitting each derivative.
    import itertools
    value = 0
    for powers, coefficient in S.Poly(S.expand(symbol**m), *symbols).terms():
        for beta in itertools.product(*(range(e + 1) for e in powers)):
            output = phase**m
            factor = coefficient
            for v, q, e, k in zip(variables, symbols, powers, beta):
                output = S.diff(output, v, k)
                factor *= S.binomial(e, k) * q**(e-k)
            value += factor * output
    return S.expand(value)


def main():
    models = [
        ('lacunary_graph_jet', X-Y**2-Y**6, z**3+x*y+1),
        ('algebraic_jet', X+Y*(1+Z/2-Z**2/8), x*z+y**2+1),
        ('nonzero_constant_phase', X+Y**2, S.Integer(3)),
        ('nonzero_constant_symbol', 2+X, x+y**2),
    ]
    controls = []
    for name, symbol, phase in models:
        degree = S.Poly(phase, *variables).total_degree()
        top = S.Add(*(coefficient*S.prod(v**e for v,e in zip(variables,powers))
                      for powers,coefficient in S.Poly(phase,*variables).terms()
                      if sum(powers)==degree))
        for m in (0, 1, 2):
            gauge = gauge_coefficient(symbol, phase, m)
            assert S.expand(gauge-direct_gauge_coefficient(symbol,phase,m)) == 0
            assert S.Poly(gauge,*variables).total_degree() <= degree*m
            scaled = S.expand(h**(degree*m)*gauge.subs({x:a/h,y:b/h,z:c/h}))
            assert all(exponent[0]>=0 for exponent,_ in S.Poly(scaled,h).terms())
            expected = S.expand(symbol**m*top.subs({x:a,y:b,z:c})**m)
            assert S.expand(scaled.subs(h,0)-expected)==0
            controls.append(dict(model=name,m=m,phase_degree=degree,
                                 gauge_equal=True,nonnegative_h=True,
                                 top_output_equal=True))
    eigen = []
    lam = S.symbols('lambda')
    A = lambda f: S.expand(t*S.diff(f,t,2)+S.diff(f,t))
    for m in range(1,8):
        assert A(t**m/S.factorial(m)**2)==t**(m-1)/S.factorial(m-1)**2
        eigen.append(m)
    H = sum((c*t*lam)**m/S.factorial(m)**2 for m in range(8))
    assert S.Poly(S.expand(A(H)-c*lam*H),t).terms()==[((7,),-c**8*lam**8/S.factorial(7)**2)]
    # Simulated relative specialization: a primitive operator survives h=0.
    O_coefficients = [h**3*(t+h),h**3*(1+h*t)]
    primitive = [S.cancel(v/h**3) for v in O_coefficients]
    assert [v.subs(h,0) for v in primitive]==[t,1]
    transverse=[]
    symbol=X-Y**2-Y**6
    for degree in (3,4):
        for marker in (S.Integer(1),x,y,x*y):
            qdegree=S.Poly(marker,*variables).total_degree()
            m=qdegree+1
            assert iterate(symbol,S.expand(marker*z**(degree*m)),m)==0
            transverse.append(dict(phase_degree=degree,marker=str(marker),m=m,zero=True))
    sharp=[]
    for q in (1,2,3):
        actual=iterate(symbol,x**q*z**(3*q),q)
        assert actual==S.factorial(q)*z**(3*q)
        sharp.append(q)
    gaps=[]
    for k in (3,4,5,6):
        left=int(S.factorial(k));right=int(S.factorial(k+1))
        gaps.append(dict(left=left,right=right,zeros=right-left-1))
        assert right-left-1>=k
    guards=0
    for args in [(X,S.Float(1)),(X,x**17)]:
        try:action(*args)
        except AssertionError:guards+=1
        else:raise AssertionError('guard did not refuse input')
    try:iterate(X,x,5)
    except AssertionError:guards+=1
    else:raise AssertionError('guard did not refuse index')
    dependencies=['verify_algebraic_formal_symbols.py','verify_holonomic_source_cutoffs.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
        sympy_version=S.__version__,gauge_degeneration_controls=controls,
        bessel_monomial_controls=eigen,relative_specialization_primitive_control=True,
        transverse_zero_controls=transverse,sharp_marker_indices=sharp,
        illustrative_factorial_gaps=gaps,guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        universal_rigidity_proved_by_finite_controls=False,
        lacunary_transcendence_proved_by_finite_controls=False,
        general_eventuality_to_linearity_decided=False,
        safe_nonholonomic_example_is_failure_of_linearity=False,
        all_ring_rigidity_claimed=False,
        general_algebraicity_recognizer_implemented=False,
        full_exceptional_MZ_decided=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    target=HERE/'formal_symbol_source_rigidity_verification.json'
    assert not target.exists(),'immutable receipt already exists'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,gauge_controls=len(controls),transverse_controls=len(transverse),
                          bessel_controls=len(eigen),sharp=sharp,guards=guards)))


if __name__=='__main__':main()
