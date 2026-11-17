"""Exact trace, principal-part and Fourier boundary controls; no universal oracle."""
from pathlib import Path
import hashlib,json
import sympy as s


def principal_part(expr,z,max_order=12):
    expr=s.sympify(expr)
    if expr.has(s.Float):raise ValueError('Exact coefficients required')
    if not expr.is_rational_function(z):raise ValueError('This calculator takes rational functions only')
    ser=s.series(expr,z,0,1).removeO().expand()
    terms=[]
    for term in s.Add.make_args(ser):
        power=term.as_powers_dict().get(z,s.Integer(0))
        if power<0:
            if power < -max_order:raise ValueError('Control pole cap exceeded')
            terms.append(term)
    return s.expand(sum(terms))


def run():
    z,t,u,v=s.symbols('z t u v')
    F=1+2*t+3*t*t+5*t**4
    def borel(p):
        p=s.Poly(p,t)
        return sum(c*s.factorial(a)/z**(a+1) for (a,),c in p.terms())
    B=borel(F)
    assert s.expand(borel(t*F)+s.diff(B,z))==0
    assert s.expand(borel(s.diff(F,t))-z*B+F.subs(t,0))==0
    # Trace in F(sqrt(1-1/z))/F: represented in basis1,u, u^2=1-1/z.
    a,b=1/z-1,1
    trace=(a+b*u+a-b*u)/2
    assert s.simplify(trace-a)==0 and principal_part(trace,z)==1/z
    deriv_u=s.diff(1-1/z,z)*u/(2*(1-1/z))
    derivative=s.diff(a,z)+s.diff(b,z)*u+b*deriv_u
    assert s.simplify((derivative+derivative.subs(u,-u))/2-s.diff(trace,z))==0
    # Cubic separable trace: u^3=z. Tr(u)=Tr(u^2)=0, Tr(1)=3.
    mult_u=s.Matrix([[0,0,z],[1,0,0],[0,1,0]])
    cubic=(s.eye(3)/z+2*mult_u+3*mult_u**2).trace()/3
    assert cubic==1/z
    # Two independent square roots: averaging four conjugates equals iterated trace.
    expr=1/z+(z+1)*u+(z-1)*v+u*v
    average=s.expand(sum(expr.subs({u:i*u,v:j*v},simultaneous=True) for i in [1,-1] for j in [1,-1])/4)
    assert average==1/z
    pp=[]
    for g in [1/z,z**-3+2/z+3+z*z,1/(z-1)+1/z**2,(z+1)/(z*z*(z-2)),z**3]:
        p=principal_part(g,z)
        assert s.simplify(principal_part(s.diff(g,z),z)-s.diff(p,z))==0
        difference=s.simplify(principal_part(z*g,z)-z*p)
        assert not difference.has(z)  # Zero or a polynomial constant, killed inW.
        pp.append(dict(rational=str(g),principal_part=str(p),multiplication_boundary=str(difference)))
    # Integer logarithms obstruct the projection: d(log z)=1/z.
    assert s.diff(s.log(z),z)==1/z
    integer_log_defect=principal_part(s.diff(s.log(z),z),z)
    assert integer_log_defect==1/z
    # Noninteger classes stay noninteger, including their logarithmic derivatives.
    noninteger=[]
    for a in [s.Rational(1,2),s.Rational(-3,2),s.Rational(2,3)]:
        assert not a.is_integer and not (a-1).is_integer and not (a+1).is_integer
        noninteger.append(dict(exponent=str(a),derivative_class=str(a-1),multiplication_class=str(a+1)))
    cauchy=s.log(z)-s.log(z-1)
    assert s.simplify((1-z)*s.diff(cauchy,z)-1/z)==0
    # Ordinary vs exponential normalization: n! gives a rational EGF but not
    # a convergent ordinary series. Only exact factorial coefficients are checked.
    factorials=[str(s.factorial(n)) for n in range(7)]
    # Algebraicity is stable under the formal divided-difference reconstruction.
    H=u**3+2*u
    quotient,rem=s.div(H.subs(u,z)-H,z-u,u);assert rem==0
    assert s.cancel(quotient/(H.subs(u,z)-H)-1/(z-u))==0
    rejected=0
    for bad in [s.log(z),s.sqrt(z),s.Float('0.5')/z,z**-13]:
        try:principal_part(bad,z)
        except ValueError:rejected+=1
        else:raise AssertionError('Expected exact rational guard refusal')
    return dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                borel_boundary_constant=str(F.subs(t,0)),quadratic_trace=str(trace),cubic_trace=str(cubic),tower_trace=str(average),
                principal_part_controls=pp,integer_log_projection_defect=str(integer_log_defect),noninteger_log_classes=noninteger,
                nonsplit_cauchy_control=str(cauchy),factorial_normalization_coefficients=factorials,guard_refusals=rejected,
                supporting_checker_hashes={},exact_principal_part_calculator_implemented=True,
                universal_theorem_proved_by_finite_controls=False,all_bad_markers_claimed_transcendental=False,
                exponential_and_ordinary_series_identified=False,generic_fibre_monodromy_substitution_used=False,
                full_exceptional_MZ_decided=False,general_algebraicity_or_monodromy_solver_implemented=False,
                independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)


if __name__=='__main__':
    receipt=run();output=Path(__file__).with_name('algebraic_resolvent_splitting_verification.json')
    if output.exists():assert json.loads(output.read_text())==receipt,'Preserve immutable receipt'
    else:output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=True,checker_sha256=receipt['checker_sha256'],trace_models=3,principal_parts=len(receipt['principal_part_controls']),noninteger_classes=len(receipt['noninteger_log_classes']),guard_refusals=receipt['guard_refusals'])))
