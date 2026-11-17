"""Exact internal replay for lower_exponent_initial_and_boundary.tex.

The universal statements are proved in that note. The bounded replay checks
parameter identities, polynomial-space ranks, exact annulus and Rouche inputs,
and moment gcds. It does not verify the entire finite support rectangle or
extrapolate the selected inputs.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import argparse
import json
import platform
import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = "sources/wolfs-2024-multiple-orthogonal-2401.08312v1.pdf"
SOURCE_HASH = "f067406cc088fbad56308b00d8add51117401379504a1109b30442b62cf82182"
TERMINAL_FIVE_CERTIFICATE = "lower_exponent_initial_and_boundary_terminal_five_certificates.json"
s, z, t = sp.symbols("s z t")

def rising(x, n):
    out = 1
    for j in range(n):
        out *= x + j
    return out

def U(n, d, shape_minus_one, zz):
    return sp.expand(sum((-1)**k * comb(n, k)
                         * rising(shape_minus_one + 1, d*k) * zz**(n-k)
                         for k in range(n + 1)))

def equation(v, n, d, zz):
    return sp.expand(zz*(v - v.subs(s, s+1))
        - rising(s+2, d-1)*((s+d*n+1)*v.subs(s, s+d)
                           - (s+d+1)*v.subs(s, s+d+1)))


def gap_bound(n):
    """Integer sufficient gap bound, for a specified terminal n>=2."""
    assert isinstance(n, int) and n >= 2
    exponent = (4*n*(n-1)-1).bit_length()
    return 2*(n-1)*exponent

def cluster_parameters(n, d):
    delta = min(Fraction(1, 4), Fraction(d, 4*n))
    cutoff = (n*d*d+d)*((2+delta)/delta)**n
    return delta, cutoff

def support_rectangle(n):
    d_bound = gap_bound(n)
    a_bound = (n*d_bound*d_bound+d_bound)*(4*n+1)**n
    return a_bound, d_bound

def ordinary_normalized_moment(j, a, d):
    # M_j/(aj)!; do not construct the potentially enormous factor (aj)!.
    return sp.Poly(sum(comb(j,k)*rising(a*j+1,d*k)*t**k
                       for k in range(j+1)), t)

def check_pairwise_moments(n, a, d, require_squarefree=False):
    polynomials = [ordinary_normalized_moment(j,a,d) for j in range(1,n+1)]
    pair_checks = 0
    squarefree_checks = 0
    for j, f in enumerate(polynomials):
        if require_squarefree:
            assert sp.gcd(f,f.diff()).degree() == 0
            squarefree_checks += 1
        for g in polynomials[:j]:
            assert sp.gcd(f,g).degree() == 0
            pair_checks += 1
    return pair_checks, squarefree_checks

def support_replay():
    gap_cases = 0
    cluster_cases = 0
    sparse_cases = 0
    pair_checks = 0
    squarefree_checks = 0
    cluster_rows = []
    # Direct exact comparisons of the annuli at the universal gap threshold.
    for n, a in product(range(2, 6), [1, 3, 20]):
        d = gap_bound(n)
        for j in range(2,n+1):
            fj = rising(a*j+1,d*j)
            fp = rising(a*(j-1)+1,d*(j-1))
            # This is precisely G_j/G_(j-1)>4j(j-1), cleared of roots.
            assert fj**(j-1) > (4*j*(j-1))**(j*(j-1))*fp**j
        pairs, _ = check_pairwise_moments(n,a,d)
        pair_checks += pairs
        gap_cases += 1

    # Actual Rouché error bounds at the stated large-a cutoff, all rational.
    for n, d in product(range(2, 7), [1, 2, 3, 5]):
        delta, cutoff = cluster_parameters(n,d)
        a = (cutoff.numerator+cutoff.denominator-1)//cutoff.denominator
        for j in range(1,n+1):
            error = sum(Fraction(comb(j,k))
                        *(Fraction(rising(a*j+1,d*k),(a*j)**(d*k))-1)
                        *(1+delta)**k for k in range(j+1))
            assert error < delta**j
        for j in range(2,n+1):
            ratio = Fraction(j,j-1)**d
            assert ratio > (1+delta)/(1-delta)
        if d >= 2 and d < gap_bound(n):
            a_bound, d_bound = support_rectangle(n)
            assert d < d_bound and cutoff <= a_bound
        pairs, _ = check_pairwise_moments(n,a,d)
        pair_checks += pairs
        cluster_rows.append({"terminal":n,"gap":d,"delta":str(delta),
                             "lower_exponent_cutoff":str(cutoff),
                             "tested_lower_exponent":a})
        cluster_cases += 1

    # Use a stronger integer condition to certify q<=1/(8n) exactly:
    # m=floor(d^2/(4(a+d))), e>2, and 2^m>=8n.
    for n,a,d in [(5,1,32),(6,2,40),(7,5,48),(5,10,80)]:
        m = d*d//(4*(a+d))
        assert 2**m >= 8*n
        q_bound = Fraction(1,2**m)
        for j in range(1,n+1):
            mass_bound = sum(comb(j,k)*q_bound**min(k,j-k)
                             for k in range(1,j))
            assert mass_bound < Fraction(2,7)
            fj = rising(a*j+1,d*j)
            for k in range(1,j):
                # b_(j,k)/binom(j,k), raised to degree j, compared directly
                # against the q-bound. All expressions here are integers.
                assert (rising(a*j+1,d*k)**j
                        *2**(m*j*min(k,j-k)) <= fj**k)
        pairs,squares = check_pairwise_moments(n,a,d,require_squarefree=True)
        pair_checks += pairs
        squarefree_checks += squares
        sparse_cases += 1

    # The cluster error criterion must reject an unsupported small-a input.
    n,d,a = 4,2,1
    delta,_ = cluster_parameters(n,d)
    error = sum(Fraction(comb(n,k))
                *(Fraction(rising(a*n+1,d*k),(a*n)**(d*k))-1)
                *(1+delta)**k for k in range(n+1))
    assert error >= delta**n
    return {"large_gap_annulus_support_cases":gap_cases,
            "large_lower_exponent_Rouche_support_cases":cluster_cases,
            "sparse_endpoint_mass_support_cases":sparse_cases,
            "direct_pairwise_moment_gcd_checks":pair_checks,
            "direct_individual_squarefreeness_gcd_checks":squarefree_checks,
            "large_lower_exponent_exact_inputs":cluster_rows}


def reciprocal_moment(n, gamma, d):
    return sp.Poly(sum(comb(n,k)*rising(gamma,d*k)*z**(n-k)
                       for k in range(n+1)), z)

def derivative_array_replay():
    identity_cases = 0
    for a,n,d,beta in product(range(1,4),range(1,5),range(1,4),
                              [sp.Rational(1),sp.Rational(1,2)]):
        m = a*n
        seed = reciprocal_moment(m,beta+m,d)
        derived = seed.diff((z,(a-1)*n))
        scale = Fraction(1,1)
        for k in range(n+1,m+1):
            scale /= k
        # n!/m!, computed without relying on the coefficient generator.
        target = reciprocal_moment(n,beta+a*n,d)
        assert sp.Poly(scale*derived.as_expr()-target.as_expr(),z).is_zero
        identity_cases += 1

    beta = sp.symbols("beta")
    seed = reciprocal_moment(6,beta+6,2)
    target = reciprocal_moment(3,beta+6,2)
    assert sp.Poly(sp.expand(seed.diff((z,3)).as_expr()/120
                              -target.as_expr()),z).is_zero

    parameter_permutations = 0
    for d,m in product(range(1,7),range(0,16)):
        q,r = divmod(m,d)
        # Remove the common beta/d from each upper parameter.
        left = sorted(Fraction(j,d)+q+(1 if j<r else 0) for j in range(d))
        right = sorted(Fraction(m+j,d) for j in range(d))
        assert left == right
        parameter_permutations += 1

    u = z+1
    family = [sp.Poly(1,z),sp.Poly(z+3,z),sp.Poly(u*u+1,z),
              sp.Poly(u**3+u/10+2,z),
              sp.Poly((u+sp.Rational(1,10))*(u**3+1),z)]
    family += [sp.Poly(sp.prod(z+10*m*m+j for j in range(1,m+1)),z)
               for m in range(5,11)]
    derivative_simple_checks = 0
    stable_minor_checks = 0
    for f in family:
        assert all(c>0 for c in f.all_coeffs())
        derivative = f
        while derivative.degree() >= 1:
            assert sp.gcd(derivative,derivative.diff()).degree() == 0
            # Independent exact leading ordinary Hurwitz minors of each
            # represented derivative, without numerical root approximations.
            degree = derivative.degree()
            coeff = derivative.all_coeffs()
            at = lambda k: coeff[k] if 0<=k<=degree else 0
            hurwitz = sp.Matrix([[at(2*(j+1)-(i+1))
                                 for j in range(degree)]
                                for i in range(degree)])
            assert all(hurwitz[:j,:j].det()>0 for j in range(1,degree+1))
            stable_minor_checks += degree
            derivative_simple_checks += 1
            derivative = derivative.diff()
    for f,g in zip(family,family[1:]):
        assert sp.gcd(f,g).degree() == 0
    diagonal_one = family[2].diff().as_expr()/2
    diagonal_two = family[4].diff((z,2)).as_expr()/12
    assert sp.expand(diagonal_one-u) == 0
    assert sp.expand(diagonal_two-u*u-u/20) == 0
    common = sp.gcd(sp.Poly(diagonal_one,z),sp.Poly(diagonal_two,z))
    assert common.degree() == 1 and common.eval(-1) == 0
    assert sp.gcd(sp.Poly(diagonal_two,z),sp.Poly(diagonal_two,z).diff()).degree() == 0
    return {"literal_derivative_embedding_cases":identity_cases,
            "symbolic_shape_derivative_embedding_cases":1,
            "fixed_base_parameter_permutations":parameter_permutations,
            "generic_counterexample_max_base_index":10,
            "generic_counterexample_derivative_simplicity_checks":derivative_simple_checks,
            "generic_counterexample_positive_Hurwitz_minor_checks":stable_minor_checks,
            "generic_counterexample_consecutive_base_gcd_checks":len(family)-1,
            "generic_counterexample_diagonal_gcd_degree":common.degree()}


def cumulant_separation_replay():
    r0,r1,r2 = sp.symbols("r0 r1 r2")
    assert sp.expand(r0*r1-r0*r0-r0*(r1-r0)) == 0
    third = r0*r1*r2-3*r0*r0*r1+2*r0**3
    positive_form = r0*(2*(r1-r0)**2+r1*(r2-2*r1+r0))
    assert sp.expand(third-positive_form) == 0

    gamma = sp.symbols("gamma")
    coefficient_profiles = 0
    for d in range(1,13):
        base = sp.Poly(sp.expand(rising(gamma,d)),gamma)
        first = sp.expand(base.as_expr().subs(gamma,gamma+d)-base.as_expr())
        second = sp.expand(first.subs(gamma,gamma+d)-first)
        assert all(c>=0 for c in sp.Poly(first,gamma).all_coeffs())
        assert all(c>=0 for c in sp.Poly(second,gamma).all_coeffs())
        variance = sp.Poly(sp.expand(base.as_expr()*first),gamma)
        kappa = sp.Poly(sp.expand(base.as_expr()
                        *(2*first**2+(base.as_expr()+first)*second)),gamma)
        for polynomial in [base,variance,kappa]:
            assert polynomial.TC() == 0 and not polynomial.is_zero
            assert all(c>=0 for c in polynomial.all_coeffs())
            assert polynomial.diff().TC()>0
        coefficient_profiles += 1

    rho,v,h,V,K,w,eta,lam = sp.symbols("rho v h V K w eta lam",real=True)
    kappa = rho**3+3*v*rho
    assert sp.expand((w+rho)*(w*w-rho*w+rho*rho+3*v)
                      -(w**3+3*v*w+kappa)) == 0
    A = h+rho/2
    quartic = (A+sp.I*eta)**4+6*(v+V)*(A+sp.I*eta)**2 \
              +4*(kappa+K)*(A+sp.I*eta)+lam
    imag = sp.expand(sp.im(sp.expand(quartic))/(4*eta)).subs(
                          eta**2,3*v+3*rho**2/4)
    certificate = K+3*v*rho+3*rho**3/4+h**3+3*rho*h*h/2 \
                    +3*V*h+3*V*rho/2
    assert sp.expand(imag-certificate) == 0
    assert all(c>0 for c in sp.Poly(certificate,rho,v,h,V,K).coeffs())
    assert sp.expand(certificate.subs({h:0,V:0,K:0})
                      -(3*v*rho+3*rho**3/4)) == 0
    cubic = (h+sp.I*eta)**3+3*(v+V)*(h+sp.I*eta)+(kappa+K)
    imag = sp.expand(sp.im(sp.expand(cubic))/eta).subs(eta**2,v)
    assert sp.expand(imag-(3*h*h+2*v+3*V)) == 0

    ordered_gcds = 0
    for d,x,increment in product(range(1,7),
                                 [sp.Rational(1,7),sp.Rational(1,2),1,3,20],
                                 [0,sp.Rational(1,3),1,4,17]):
        y = x+increment
        f = sum(comb(3,k)*rising(x,d*k)*t**k for k in range(4))
        g = sum(comb(4,k)*rising(y,d*k)*t**k for k in range(5))
        assert sp.Poly(sp.gcd(f,g),t).degree() == 0
        ordered_gcds += 1
    return {"abstract_central_moment_identities":2,
            "gamma_power_nonnegative_coefficient_profiles":coefficient_profiles,
            "abstract_cubic_factor_identity":1,
            "abstract_positive_imaginary_part_identities":2,
            "ordered_shape_three_four_gcd_cases":ordered_gcds}


def finite_gamma_convolution_replay():
    alpha1,theta1 = sp.Rational(125,9),sp.Rational(1,2)
    alpha2,theta2 = sp.Rational(1,3249),sp.Rational(19,2)
    assert min(alpha1,alpha2,theta1,theta2)>0
    # Derive raw moments literally from the two independent Gamma variables.
    moments = [sum(comb(n,k)*rising(alpha1,k)*theta1**k
                    *rising(alpha2,n-k)*theta2**(n-k) for k in range(n+1))
               for n in range(6)]
    cumulants = [0]
    for n in range(1,6):
        cumulants.append(sp.cancel(moments[n]-sum(
            comb(n-1,k-1)*cumulants[k]*moments[n-k] for k in range(1,n))))
    for n in range(1,6):
        factorial = 1
        for k in range(1,n):
            factorial *= k
        assert cumulants[n] == factorial*(alpha1*theta1**n+alpha2*theta2**n)
    assert cumulants[1:] == [sp.Rational(132,19),sp.Rational(7,2),
                            4,sp.Rational(81,4),582]
    w = sp.symbols("w")
    polynomials = [sp.Poly(sum(comb(n,k)*moments[k]*z**(n-k)
                               for k in range(n+1)),z) for n in [4,5]]
    centered = [sp.Poly(sp.expand(p.as_expr().subs(z,w-moments[1])),w)
                for p in polynomials]
    quadratic = w*w-w+19
    assert sp.expand(centered[0].as_expr()-quadratic*(w*w+w+3)) == 0
    assert sp.expand(centered[1].as_expr()-(w+2)*quadratic**2) == 0
    common = sp.gcd(*centered)
    assert common == sp.Poly(quadratic,w) and sp.discriminant(quadratic,w) == -75
    assert sp.gcd(centered[1],centered[1].diff()).degree() == 2
    assert sp.Poly(polynomials[1].diff().as_expr()-5*polynomials[0].as_expr(),z).is_zero
    # Equal-scale control collapses the distribution to one Gamma variable.
    single = [sp.Poly(sum(comb(n,k)*rising(alpha1+alpha2,k)*theta1**k*z**(n-k)
                          for k in range(n+1)),z) for n in [4,5]]
    assert sp.gcd(*single).degree() == 0
    return {"positive_shape_scale_parameters":[[str(alpha1),str(theta1)],
                                               [str(alpha2),str(theta2)]],
            "literal_independent_raw_moment_orders":list(range(0,6)),
            "moment_to_cumulant_inversion_orders":list(range(1,6)),
            "common_quadratic_centered":"w^2-w+19",
            "common_quadratic_discriminant":-75,
            "four_five_common_factor_degree":common.degree(),
            "fifth_moment_repeated_factor_degree":2,
            "single_equal_scale_Gamma_control_coprime":True,
            "HCM_density_asserted":False,
            "Gamma_power_or_factorial_counterexample":False}


def generate_terminal_five_certificates():
    a = sp.symbols("a")
    profiles = []
    for d in range(2,17):
        f = sp.Poly(sum(comb(4,k)*sp.rf(4*a+1,d*k)*t**k
                        for k in range(5)),t,domain=sp.ZZ[a])
        g = sp.Poly(sum(comb(5,k)*sp.rf(5*a+1,d*k)*t**k
                        for k in range(6)),t,domain=sp.ZZ[a])
        polynomial = sp.Poly(f.resultant(g),a)
        content,primitive = polynomial.primitive()
        assert content>0 and primitive.degree()==20*d
        coefficients = [int(primitive.nth(k)) for k in range(20*d+1)]
        assert all(c>0 for c in coefficients)
        profiles.append({"gap":d,"content":str(content),
                         "degree_bound":20*d,
                         "primitive_coefficients_ascending":list(map(str,coefficients))})
        print(json.dumps({"generated_terminal_five_gap":d,
                          "degree":20*d,"positive_coefficients":len(coefficients)}),flush=True)
    data = {"scope":"ordinary factorial moments; terminal five; fixed gaps2..16; every lower exponent",
            "polynomials":"F_j,d(a,t)=sum_k binom(j,k)(ja+1)_(dk)t^k",
            "resultant":"Res_t(F_4,d,F_5,d)",
            "coefficient_domain":"Z[a]",
            "generation":"SymPy exact resultant over ZZ[a]; generator is not trusted by verifier",
            "sympy":sp.__version__,
            "verification":"literal integer Sylvester Bareiss at20d+1 nodes; degree bound20d",
            "profiles":profiles}
    (HERE/TERMINAL_FIVE_CERTIFICATE).write_text(json.dumps(data,indent=2)+"\n",
                                               encoding="utf-8")

def bareiss_integer_determinant(matrix):
    work = [row[:] for row in matrix]
    n = len(work)
    previous,sign = 1,1
    for k in range(n-1):
        if work[k][k] == 0:
            pivot_row = next((i for i in range(k+1,n) if work[i][k]),None)
            if pivot_row is None:
                return 0
            work[k],work[pivot_row] = work[pivot_row],work[k]
            sign = -sign
        pivot = work[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = pivot*work[i][j]-work[i][k]*work[k][j]
                assert numerator % previous == 0
                work[i][j] = numerator//previous
            work[i][k] = 0
        previous = pivot
    return sign*work[-1][-1]

def literal_terminal_five_resultant(lower_exponent,d):
    rows = []
    for n in [4,5]:
        ratio,ascending = 1,[1]
        for k in range(1,n+1):
            for ell in range(d*(k-1)+1,d*k+1):
                ratio *= n*lower_exponent+ell
            ascending.append(comb(n,k)*ratio)
        rows.append(ascending[::-1])
    f,g = rows
    matrix = [[0]*i+f+[0]*(4-i) for i in range(5)]
    matrix += [[0]*i+g+[0]*(3-i) for i in range(4)]
    assert len(matrix)==9 and all(len(row)==9 for row in matrix)
    return bareiss_integer_determinant(matrix)

def verify_terminal_five_profile(profile,nodes=None):
    d = profile["gap"]
    assert isinstance(d,int) and 2<=d<=16
    content = int(profile["content"])
    coeff = list(map(int,profile["primitive_coefficients_ascending"]))
    assert profile["degree_bound"]==20*d and len(coeff)==20*d+1
    assert content>0 and all(c>0 for c in coeff)
    if nodes is None:
        nodes = range(20*d+1)
    checked = 0
    for a in nodes:
        value = 0
        for c in reversed(coeff):
            value = value*a+c
        assert literal_terminal_five_resultant(a,d) == content*value
        checked += 1
    return checked

def terminal_five_replay():
    data = json.loads((HERE/TERMINAL_FIVE_CERTIFICATE).read_text(encoding="utf-8"))
    profiles = data["profiles"]
    assert [row["gap"] for row in profiles] == list(range(2,17))
    count = sum(verify_terminal_five_profile(profile) for profile in profiles)
    assert count==2715
    b5,b4 = rising(6,85),rising(5,68)
    assert 63**4<2*53**4 and 31**5<2*27**5
    assert b5**4*40**20 > 1431**20*b4**5
    # Keep positivity intact while corrupting an actual coefficient:
    # the independent literal determinant must still reject the identity.
    bad = dict(profiles[0])
    bad["primitive_coefficients_ascending"] = list(bad["primitive_coefficients_ascending"])
    bad["primitive_coefficients_ascending"][0] = str(
        int(bad["primitive_coefficients_ascending"][0])+1)
    try:
        verify_terminal_five_profile(bad,nodes=[0])
    except AssertionError:
        pass
    else:
        raise AssertionError("corrupted positive resultant certificate accepted")
    return {"fixed_gaps":list(range(2,17)),
            "whole_resultant_polynomial_degree_bound":"20*d",
            "literal_integer_Bareiss_determinant_evaluations":count,
            "positive_primitive_coefficients":count,
            "all_lower_exponents_covered":"a>=0 real at fixed gaps2..16",
            "large_gap_cutoff":{"terminal":5,"gap_at_least":17,"lower_exponent_at_least":1,
                               "rational_annuli":["53/10","27/4"],
                               "factorial_integer_comparison":"pass"},
            "corrupted_positive_polynomial_certificate_rejected":True,
            "ordinary_binomial_terminal_five_proved":True,
            "all_positive_Gamma_shapes_terminal_five_proved":False}

def replay():
    assert sha256((HERE/SOURCE).read_bytes()).hexdigest() == SOURCE_HASH
    cases = 0
    for n, d in product(range(2, 7), range(1, 5)):
        v = U(n-1, d, s, z)
        m = d*(n-1)
        assert equation(v, n, d, z) == 0
        assert sp.Poly(v, s).degree() == m
        gap = sp.Poly(v - (-1)**(n-1)*rising(s+1, m), s)
        assert gap.is_zero or gap.degree() <= m-d
        # Independent coefficient-matrix nullspace check at nonzero and zero z.
        # The proof's degree argument, not these samples, gives all complex z.
        for zz in [0, 1, -2, sp.I]:
            columns = [sp.Poly(equation(s**j, n, d, zz), s)
                       for j in range(m+1)]
            max_degree = max(p.degree() for p in columns if not p.is_zero)
            mat = sp.Matrix([[p.nth(i) for p in columns]
                             for i in range(max_degree+1)])
            assert mat.rank() == m
            residual = mat * sp.Matrix([sp.Poly(v.subs(z, zz), s).nth(j)
                                        for j in range(m+1)])
            assert all(sp.expand(entry) == 0 for entry in residual)
        cases += 1

    inequality_cases = 0
    for d, x, delta in product(range(1, 9),
                               [Fraction(1, 7), Fraction(1, 2), 1, 3, 20],
                               [0, Fraction(1, 3), 1, 4, 17]):
        y = x + delta
        assert rising(y, 3*d) > rising(x, 2*d)*(3*rising(y, d)-2*rising(x, d))
        inequality_cases += 1

    moment_cases = 0
    quartic_moment_cases = 0
    for a, d, beta in product(range(0, 6), range(1, 7),
                              [sp.Rational(1, 3), sp.Rational(1, 2), 1, 3]):
        f2 = sum(comb(2, k)*rising(beta+2*a, d*k)*t**k for k in range(3))
        f3 = sum(comb(3, k)*rising(beta+3*a, d*k)*t**k for k in range(4))
        assert sp.Poly(sp.gcd(f2, f3), t).degree() == 0
        # The finite obstruction equations are direct consequences of a
        # hypothetical consecutive common zero; replay the underlying identity.
        n = 3
        ll = beta + a*(n-1)-1
        mm = ll + a
        assert sp.expand(U(n,d,mm,z) - z*U(n-1,d,mm,z)
                         + rising(mm+1,d)*U(n-1,d,mm+d,z)) == 0
        f4 = sum(comb(4,k)*rising(beta+4*a,d*k)*t**k for k in range(5))
        assert sp.Poly(sp.gcd(f3,f4),t).degree() == 0
        quartic_moment_cases += 1
        moment_cases += 1

    # Corrupt a selected polynomial without changing its leading coefficient.
    # The recurrence must reject it rather than accepting a degree-only control.
    assert equation(U(3, 3, s, z)+1, 4, 3, z) != 0
    support_counts = support_replay()
    derivative_counts = derivative_array_replay()
    cumulant_counts = cumulant_separation_replay()
    convolution_counts = finite_gamma_convolution_replay()
    fifth_counts = terminal_five_replay()
    report = {
        "kind": "internal bounded exact replay; universal proof is the TeX source",
        "claims": {
            "two_three_separation": "all integer d>=1 and real y>=x>0; elementary proof",
            "three_four_separation": "all integer d>=1 and real y>=x>0; monotone cumulant proof",
            "full_ordinary_objective_equivalence": "consecutive gcd in finite rectangles for every N>=6; infinitely many N remain",
            "lower_exponent_initial_block": "all a>=0,d>=1 integers and beta>0, terminal<=4",
            "parameter_polynomial_selection": "all N>=2,d>=1,z complex; elementary degree proof",
            "finite_support_reduction": "each terminal N>=2: a<A_N and 2<=d<D_N; ordinary factorials",
            "large_gap_pairwise_coprimality": "all a>=1, d>=D_N, orders 1..N",
            "large_lower_exponent_pairwise_coprimality": "all d>=1, a>=A(N,d), orders 1..N",
            "sparse_endpoint_mass_criterion": "d^2>=4(a+d)log(8N): simplicity and pairwise coprimality through N",
            "diagonal_derivative_embedding": "all a,d,N>=1 integers and beta>0; fixed lower-one seed S_m",
            "generic_preservation_shortcut": "counterexample despite simplicity and Hurwitz stability of every derivative and consecutive base coprimality",
            "generic_Gamma_convolution_simplicity": "explicit finite positive Gamma convolution counterexample at fifth order",
            "ordinary_binomial_terminal_five": "all integer a>=0,d>=1; monotone center cutoff and exact polynomial certificates",
            "all_positive_Gamma_shapes_terminal_five_closed": False,
            "all_index_binomial_sfc": False,
        },
        "coverage": {
            "symbolic_recurrence_cases": cases,
            "polynomial_kernel_rank_specializations_per_case": 4,
            "rational_inequality_cases": inequality_cases,
            "exact_moment_gcd_cases": moment_cases,
            "exact_three_four_moment_gcd_cases": quartic_moment_cases,
        },
        "support_bound_replay": support_counts,
        "diagonal_derivative_replay": derivative_counts,
        "cumulant_separation_replay": cumulant_counts,
        "finite_Gamma_convolution_replay": convolution_counts,
        "ordinary_terminal_five_replay": fifth_counts,
        "representation_source": {"file": SOURCE, "sha256": SOURCE_HASH,
                                  "theorem": "Wolfs arXiv2401.08312v1 Theorem4.6 p=0,q=d",
                                  "scope": "representation attribution only; coefficient identity rederived; no zero theorem imported",
                                  "cached_primary_source_hash_checked": True},
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "source_sha256": {
            name: sha256((HERE/name).read_bytes()).hexdigest()
            for name in ["lower_exponent_initial_and_boundary.tex",
                         "lower_exponent_initial_and_boundary.py", SOURCE, TERMINAL_FIVE_CERTIFICATE]
        },
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
        "passed": True,
    }
    return report

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate-terminal-five",action="store_true",
                        help="regenerate exact resultant inputs before the full replay")
    args = parser.parse_args()
    if args.generate_terminal_five:
        generate_terminal_five_certificates()
    report = replay()
    report["command"] = ("python research/factorial/image-rigidity/"
                         "lower_exponent_initial_and_boundary.py"
                         +(" --generate-terminal-five" if args.generate_terminal_five else ""))
    report["certificate_regeneration_command"] = (
        "python research/factorial/image-rigidity/"
        "lower_exponent_initial_and_boundary.py --generate-terminal-five")
    (HERE/"lower_exponent_initial_and_boundary_verification.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
