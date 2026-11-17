"""Exact internal replay for lower_exponent_initial_and_boundary.tex.

The universal statements are proved in that note. The bounded replay checks
parameter identities, polynomial-space ranks, exact annulus and Rouche inputs,
and moment gcds. It does not verify the entire finite support rectangle or
extrapolate the selected inputs.
"""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb, factorial
from pathlib import Path
import argparse
import json
import platform
import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = "sources/wolfs-2024-multiple-orthogonal-2401.08312v1.pdf"
SOURCE_HASH = "f067406cc088fbad56308b00d8add51117401379504a1109b30442b62cf82182"
HIGHER_ADJOINT_SOURCE = "HIGHER_GAP_ADJOINT_RANK.tex"
HIGHER_ADJOINT_HASH = "ec19415938c9b0b92e7d365fb4a4c66b701fdf36d73883b68f7cfb2971d7ec28"
ARITHMETIC_SOURCE = "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex"
ARITHMETIC_HASH = "7eaedb0f2f96fef505b7045d71c2f2a8ffb972d9a807272e5ac9defc62d22351"
CUBIC_SOURCE = "CUBIC_ADJOINT_STABILITY_SFC.tex"
CUBIC_SOURCE_HASH = "ea2510e511545b9dd7a49a0796d75c45575687e60770c3359d38c5bc99cb4125"
FRACTION_SOURCE = "sources/petreolle-sokal-zhu-1807.03271v2.pdf"
FRACTION_SOURCE_HASH = "37b51eb4c7b47951460c0287e88b295be366acf76db55b1bd594dc8120554c1c"
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



def gap_two_shape_offset_replay():
    """Exact symbolic controls for the all-order offsets2,3 proof.

    Universal root separation is the analytic argument in the TeX source;
    the finite Routh checks below are controls, not its justification.
    """
    m,MM,gamma = sp.symbols("m M gamma")
    aode = (2-2*MM)/z
    bode = MM*(MM-1)/z**2-1
    code = 2*m/z
    def derivative_vector(v):
        r,r1,r2 = v
        return [sp.diff(r,z)-code*r2,
                r+sp.diff(r1,z)-bode*r2,
                r1+sp.diff(r2,z)-aode*r2]
    sequence = [[sp.Integer(1),sp.Integer(0),sp.Integer(0)]]
    for j in range(5):
        v = sequence[-1]
        dv = derivative_vector(v)
        sequence.append([sp.cancel((MM+j)*v[i]-z*dv[i]) for i in range(3)])
    AA = {2:(MM+2)*(MM+3),3:(MM+2)*(MM+3)*(MM+4)}
    qq = {
        2:-MM*(MM+1)*AA[2]+(MM*(MM+1)-2*m*(2*MM+3))*z**2,
        3:-MM*(MM+1)*AA[3]+(MM+2)*(MM*(MM+1)-6*m*(MM+2))*z**2}
    bb = {2:2*z*(MM*AA[2]+(m+1)*z**2),
          3:2*z*(MM*AA[3]+(m+1)*(3*MM+4)*z**2)}
    cc = {2:-z**2*AA[2],3:-z**2*(AA[3]+2*(m+1)*z**2)}
    for h in [2,3]:
        reduced = [sp.cancel(z**2*sequence[h][i]-sequence[h+2][i])
                   for i in range(3)]
        assert all(sp.cancel(x-y)==0 for x,y in zip(reduced,[qq[h],bb[h],cc[h]]))

    b0,b1 = sp.symbols("B Bprime")
    for h in [2,3]:
        bp = (z+MM+1)*b0-z*b1
        bm = (-z+MM+1)*b0+z*b1
        if h==2:
            hp = AA[h]*bp+2*(m+1)*z**2*b0
            hm = AA[h]*bm+2*(m+1)*z**2*b0
        else:
            hp = (AA[h]+2*(m+1)*z**2)*bp+2*(m+1)*(MM+4)*z**2*b0
            hm = (AA[h]+2*(m+1)*z**2)*bm+2*(m+1)*(MM+4)*z**2*b0
        plus = cc[h]*(b1+(-1+(1-MM)/z)*b0)+(bb[h]-cc[h]*aode)*b0
        minus = cc[h]*(-b1+(1+(1-MM)/z)*b0)+(bb[h]-cc[h]*aode)*b0
        assert sp.cancel(plus-z*hp)==0
        assert sp.cancel(minus-z*hm)==0

    def transformed_moment(order,shape):
        return sp.expand(sum((-1)**k*comb(order,k)*rising(shape,2*k)
                             *z**(2*(order-k)) for k in range(order+1)))
    def auxiliaries(order,shape):
        values = [sp.Integer(1)]
        for j in range(order):
            values.append(sp.expand((z+2*j+shape+1)*values[-1]
                                     -z*sp.diff(values[-1],z)))
        return values

    literal_operators = literal_odes = literal_concomitants = 0
    for order in range(1,7):
        r = transformed_moment(order,gamma)
        parameter = gamma+2*order
        aa,bbb,ccc = aode.subs(MM,parameter),bode.subs(MM,parameter),code.subs(m,order)
        assert sp.cancel(sp.diff(r,z,3)+aa*sp.diff(r,z,2)
                         +bbb*sp.diff(r,z)+ccc*r)==0
        literal_odes += 1
        for h in [2,3]:
            substitutions = {m:order,MM:parameter}
            expected = rising(gamma,h)*transformed_moment(order+1,gamma+h)
            actual = (qq[h].subs(substitutions)*r
                      +bb[h].subs(substitutions)*sp.diff(r,z)
                      +cc[h].subs(substitutions)*sp.diff(r,z,2))
            assert sp.expand(expected-actual)==0
            literal_operators += 1
        if order<=4:
            base = auxiliaries(order,gamma)[-1]
            for sign in [1,-1]:
                b = base.subs(z,sign*z)
                ell = -sign+(1-parameter)/z
                f1 = sp.diff(b,z)+ell*b
                f2 = sp.diff(b,z,2)+2*ell*sp.diff(b,z)+(sp.diff(ell,z)+ell**2)*b
                concomitant = (b*sp.diff(r,z,2)+(aa*b-f1)*sp.diff(r,z)
                               +(f2-sp.diff(aa,z)*b-aa*f1+bbb*b)*r)
                assert sp.cancel(concomitant)==0
                literal_concomitants += 1
    # A changed derivative coefficient must fail the literal operator identity.
    r = transformed_moment(2,gamma)
    parameter = gamma+4
    substitutions = {m:2,MM:parameter}
    expected = rising(gamma,3)*transformed_moment(3,gamma+3)
    actual = (qq[3].subs(substitutions)*r
              +(bb[3].subs(substitutions)+z)*sp.diff(r,z)
              +cc[3].subs(substitutions)*sp.diff(r,z,2))
    assert sp.expand(expected-actual)!=0

    hurwitz_profiles = hurwitz_minors = moment_gcds = 0
    for h,order in product([2,3],range(1,7)):
        for shape in [2*order-2+sp.Rational(1,7),h*order+1,h*order+5]:
            parameter = shape+2*order
            values = auxiliaries(order+1,shape)
            base,next_base = values[-2:]
            if h==2:
                hh = ((parameter+2)*(parameter+3)*next_base
                      +2*(order+1)*z**2*base)
            else:
                aa = (parameter+2)*(parameter+3)*(parameter+4)
                hh = ((aa+2*(order+1)*z**2)*next_base
                      +2*(order+1)*(parameter+4)*z**2*base)
            poly = sp.Poly(hh,z)
            coefficients = poly.all_coeffs()
            assert all(c>0 for c in coefficients)
            at = lambda k:coefficients[k] if 0<=k<=poly.degree() else 0
            matrix = sp.Matrix([[at(2*(j+1)-(i+1)) for j in range(poly.degree())]
                                for i in range(poly.degree())])
            for size in range(1,poly.degree()+1):
                assert matrix[:size,:size].det(method="domain-ge")>0
                hurwitz_minors += 1
            hurwitz_profiles += 1
            lower = sp.Poly(transformed_moment(order,shape),z)
            upper = sp.Poly(transformed_moment(order+1,shape+h),z)
            assert sp.gcd(lower,upper).degree()==0
            moment_gcds += 1
    return {"all_order_theorem":"T_(m,gamma,2) andT_(m+1,gamma+h,2) coprime for h2,3 andgamma>2m-2",
            "all_positive_Gamma_families":["X^2(A+B X^2)","X^3(A+B X^2)"],
            "symbolic_Ore_reduction_component_identities":6,
            "symbolic_transformed_adjoint_identities":4,
            "literal_symbolic_shape_operator_identities":literal_operators,
            "literal_symbolic_shape_cubic_ODE_identities":literal_odes,
            "literal_symbolic_shape_zero_concomitants":literal_concomitants,
            "bounded_Hurwitz_profiles":hurwitz_profiles,
            "positive_exact_Hurwitz_minors":hurwitz_minors,
            "direct_transformed_moment_gcd_controls":moment_gcds,
            "corrupted_derivative_coefficient_rejected":True,
            "all_gap_two_lower_exponents_proved":False}


def mixed_quintic_replay():
    """Focused exact controls for the mixed rank theorem and 8/9 corollary."""
    m,MM,gamma,omega,uu = sp.symbols("m M gamma omega u")
    coeff = [3*m/z,-MM*(MM-1)*(MM+1)/z**3-1,
             3*MM*(MM-1)/z**2,-3*(MM-1)/z,sp.Integer(1)]
    def derivative(v):
        out = [sp.diff(v[i],z)-coeff[i]*v[3] for i in range(4)]
        for i in range(1,4):
            out[i] += v[i-1]
        return [sp.cancel(x) for x in out]
    seq = [[sp.Integer(1),sp.Integer(0),sp.Integer(0),sp.Integer(0)]]
    for j in range(5):
        v = seq[-1];dv = derivative(v)
        seq.append([sp.cancel((MM+j)*v[i]-z*dv[i]) for i in range(4)])
    AA = (MM+3)*(MM+4)
    expected = [rising(MM,5)+(MM*(MM+1)-6*m*(MM+2))*z**3,
                3*z*(-MM*(MM+1)*AA+(m+1)*z**3),
                3*MM*z**2*AA,-z**3*AA]
    assert all(sp.cancel(z**3*seq[2][i]+seq[5][i]-expected[i])==0
               for i in range(4))
    _,bb,cc,dd = expected
    B,B1,B2 = sp.symbols("B B1 B2")
    ell = -omega+(1-MM)/z
    f1 = omega*B1+ell*B
    f2 = omega**2*B2+2*ell*omega*B1+(sp.diff(ell,z)+ell**2)*B
    U = dd*(coeff[3]*B-f1)-cc*B
    V = dd*(f2-coeff[3]*f1+(coeff[2]-sp.diff(coeff[3],z))*B)-bb*B
    xi = omega*z
    C = (xi+MM+2)*B-xi*B1
    D = (xi+MM+5)*C-xi*(B+(xi+MM+1)*B1-xi*B2)
    H = D+3*(m+1)/AA*xi**3*B
    assert sp.cancel(U+z*z*AA*C)==0
    residual = sp.together(V+(3*MM+5)/z*U+z*AA*H).as_numer_denom()[0]
    assert sp.rem(residual,omega**3-1,omega)==0

    def transformed(order,shape):
        return sp.expand(sum(comb(order,k)*rising(shape,3*k)*z**(3*(order-k))
                             for k in range(order+1)))
    literal = 0
    for order in range(1,5):
        r = transformed(order,gamma);parameter = gamma+3*order
        sub = {m:order,MM:parameter}
        assert sp.cancel(sum(coeff[i].subs(sub)*sp.diff(r,z,i)
                             for i in range(5)))==0
        assert sp.expand(rising(gamma,2)*transformed(order+1,gamma+2)
                         -sum(expected[i].subs(sub)*sp.diff(r,z,i)
                              for i in range(4)))==0
        literal += 1
    def aux(order,shape):
        values = [sp.Integer(1)]
        for j in range(order):
            values.append(sp.expand((z+3*j+shape+2)*values[-1]
                                     -z*sp.diff(values[-1],z)))
        return values
    def cyclic(f):
        p = sp.Poly(f,z)
        return [sum(p.nth(j+3*k)*uu**k for k in range(p.degree()//3+1))
                for j in range(3)]
    def minors(left,right):
        return [sp.expand(left[i]*right[j]-left[j]*right[i])
                for i,j in [(0,1),(0,2),(1,2)]]
    values = aux(3,gamma)
    ap = (gamma+6)*(gamma+7)
    left = cyclic(values[2])
    right = cyclic(ap*values[3]+6*z**3*values[1])
    quotient = sp.Poly(3*(gamma+2)**4+36*(gamma+2)**3
                       +71*(gamma+2)**2+90*(gamma+2)+40,gamma)
    assert sp.gcd(quotient,sp.Poly(ap*rising(gamma,3),gamma)).degree()==0
    for minor in minors(left,right):
        assert sp.rem(minor.subs(uu,-rising(gamma,3)),quotient.as_expr(),gamma)==0
    bad = cyclic(ap*values[3]+7*z**3*values[1])
    assert any(sp.rem(x.subs(uu,-rising(gamma,3)),quotient.as_expr(),gamma)!=0
               for x in minors(left,bad))
    # Complete finite residues of the reduced A: A^3=A and positive powers
    # have period2, while gamma=2m+1 depends only on m mod3.
    for order in range(1,7):
        c0 = (2*order+3)%3
        matrix = sp.Matrix([[c0,0,0],[1,(c0-1)%3,0],[0,1,(c0-2)%3]])
        assert all(int(v)%3==0 for v in matrix**3-matrix)
        vectors = [matrix**(order+1)*sp.Matrix([1,0,0]),
                   matrix**(order+2)*sp.Matrix([1,0,0])]
        assert any(int(v)%3 for v in minors(list(vectors[0]),list(vectors[1])))
    for terminal in range(9):
        assert ((5*terminal-1)%9==0)==(terminal==2)
    tested = [3,4,8,10,12,13]
    for terminal in tested:
        assert terminal%9!=2
        lower = ordinary_normalized_moment(terminal-1,2,3)
        upper = ordinary_normalized_moment(terminal,2,3)
        assert sp.gcd(lower,upper).degree()==0
    return {"generic_Ore_component_identities":4,
            "generic_stripped_adjoint_identities":2,
            "literal_symbolic_Gamma_fourth_order_ODEs":literal,
            "literal_symbolic_Gamma_upper_operator_identities":literal,
            "algebraic_common_root_rank_control_minors":3,
            "corrupted_mixed_coefficient_rejected":True,
            "complete_mod3_reduced_operator_residue_controls":6,
            "mod9_support_residue_controls":9,
            "bounded_ordinary_terminal_gcd_controls":tested,
            "proved_ordinary_terminal_scope":"N>=2 andN not2mod9",
            "degree_at_most_five_equivalent_remaining_scope":"support{2,5}, N>=11,N2mod9; consecutive gcd only",
            "remaining_progression_closed":False}

def p3_index_interpolation_replay():
    """Exact Gauss-tail, congruence and Hensel controls for the interpolation theorem."""
    T,w=sp.symbols("T w")
    def vp(a):
        if a==0:return 10**9
        q=Fraction(a);n,d=q.numerator,q.denominator;v=0
        while n%3==0:v+=1;n//=3
        while d%3==0:v-=1;d//=3
        return v
    def vf(k):
        v=0
        while k:k//=3;v+=k
        return v
    def gauss(p):
        return min(vp(c) for c in sp.Poly(p,T).all_coeffs() if c)
    def rise(a,k):
        out=sp.Poly(1,T,domain=sp.QQ)
        for j in range(k):out*=sp.Poly(a+j,T,domain=sp.QQ)
        return out
    def modpoly(p,modulus):
        pp=sp.Poly(p,T,w,domain=sp.QQ)
        return sp.Poly.from_dict({idx:(int(c.p)%modulus)*pow(int(c.q),-1,modulus)%modulus
                                 for idx,c in pp.terms()},(T,w),domain=sp.ZZ).as_expr()
    def equalmod(left,right,modulus):
        return all(int(c)%modulus==0 for c in sp.Poly(left-right,T,w).coeffs())
    m=1+27*T;h=1+18*T
    f=[];g=[]
    for k in range(13):
        falling=sp.Poly(1,T,domain=sp.QQ)
        falling_up=sp.Poly(1,T,domain=sp.QQ)
        unit=sp.Poly(1,T,domain=sp.QQ);unit_up=unit
        for j in range(k):
            falling*=sp.Poly(m-j,T,domain=sp.QQ)
            falling_up*=sp.Poly(m+1-j,T,domain=sp.QQ)
            unit*=sp.Poly((3*(h+j)+1)*(3*(h+j)+2),T,domain=sp.QQ)
            unit_up*=sp.Poly((3*(h+j)+2)*(3*(h+j)+4),T,domain=sp.QQ)
        ff=falling*rise(h,k)*unit/sp.factorial(k)
        gg=falling_up*rise(h+1,k)*unit_up/sp.factorial(k)
        f.append(sp.Poly(ff,T).as_expr());g.append(sp.Poly(gg,T).as_expr())
        assert gauss(f[-1])>=0 and gauss(g[-1])>=0
    def aval(k):
        return sum((k-2)//3**r+1 for r in range(1,4))+sum(k//3**r for r in range(1,3))-vf(k)
    def bval(k):
        return sum((k-3)//3**r+1 for r in range(1,4))+sum((k+1)//3**r for r in range(1,3))-vf(k)
    for k in range(2,13):assert gauss(f[k])==aval(k)>=3
    for k in range(3,13):assert gauss(g[k])==bval(k)>=3
    assert Fraction(23*13-133,54)>3 and Fraction(23*13-135,54)>3
    F=modpoly(sum(f[k]*w**k for k in range(13)),27)
    G=modpoly(sum(g[k]*w**k for k in range(13)),27)
    assert equalmod(F,1+(20+9*T)*w,27)
    assert equalmod(G,1+(5+18*T)*w+(6+9*T)*w*w,27)
    W=4+9*T
    assert equalmod(F.subs(w,W),0,27)
    assert equalmod(G.subs(w,W),9*(1-T),27)
    assert not equalmod((G+9*w).subs(w,W),9*(1-T),27)

    def literal_coeff(n,gamma,K):
        modulus=3**K;out=[];raw=1
        for k in range(min(n+1,3*K)):
            if k:
                for j in range(3):raw*=gamma+3*(k-1)+j
            assert raw%3**k==0
            out.append(comb(n,k)*(raw//3**k)%modulus)
        return out
    def ev(c,x,modulus):
        out=0
        for a in reversed(c):out=(out*x+a)%modulus
        return out
    def residual(index,K):
        coeff=literal_coeff(index,2*index+1,K);root=1
        for e in range(2,K+1):
            modulus=3**e
            der=sum(k*coeff[k]*pow(root,k-1,modulus) for k in range(1,len(coeff)))%modulus
            root=(root-ev(coeff,root,modulus)*pow(der,-1,modulus))%modulus
        assert ev(coeff,root,3**K)==0
        return ev(literal_coeff(index+1,2*index+3,K),root,3**K)
    rootT=0;digits=[]
    for level in range(1,21):
        place=3**(level-1)
        candidates=[rootT+j*place for j in range(3)
                    if residual(1+27*(rootT+j*place),level+2)==0]
        assert len(candidates)==1
        rootT=candidates[0];digits.append({"precision":level,"T_residue":rootT})
    assert rootT%9==7
    x=sp.Symbol("x")
    def literal_moment(n,gamma):
        raw=1;coeff=[]
        for k in range(n+1):
            if k:
                for j in range(3):raw*=gamma+3*(k-1)+j
            coeff.append(comb(n,k)*raw)
        return sp.Poly.from_list(coeff[::-1],x,domain=sp.ZZ)
    endpoint=[literal_moment(28,57),literal_moment(29,59)]
    assert sp.gcd(*endpoint).degree()==0
    report={"scope":"actual quintic moment interpolation on m=1+27T, T,w in Z3; not all root-radius branches",
            "universal_tail_proof":{"lower_gauss_bound":"(23*k-133)/54",
                                   "upper_gauss_bound":"(23*k-135)/54",
                                   "cutoff":13,"valuation_at_least":3},
            "exact_low_Gauss_profiles":{"lower_k2_to12":[aval(k) for k in range(2,13)],
                                        "upper_k3_to12":[bval(k) for k in range(3,13)]},
            "coefficientwise_mod27_congruences_verified":True,
            "two_stage_contraction_proved_in_companion_source":True,
            "Hensel_digit_controls":digits,
            "positive_integer_endpoint_T1_orders":[28,29],
            "endpoint_gcd_degree":0,
            "T_star_is_ordinary_integer_established":False,
            "SFC_counterexample":False,"passed":True}
    return report

def formal_cumulant_fraction_replay():
    """Compare raw-moment cumulants, typed trees, formal ratios and Dyck paths.

    These are bounded exact symbolic controls of the universal formal proof;
    sampled minors do not establish the imported all-minors theorem.
    """
    gamma = sp.Symbol("gamma")
    zero,one = sp.Poly(0,gamma,domain=sp.ZZ),sp.Poly(1,gamma,domain=sp.ZZ)
    max_n = 6
    counts = {"symbolic_cumulant_forest_coefficients":0,
              "formal_ratio_coefficients":0,
              "independent_partial_Dyck_path_coefficients":0,
              "sampled_coefficientwise_nonnegative_Hankel_minors":0}
    def gamma_rising(shift,length):
        value = one
        for j in range(length):
            value *= sp.Poly(gamma+shift+j,gamma,domain=sp.ZZ)
        return value
    def convolution(left,right,degree):
        out = [zero]*(degree+1)
        for i in range(min(len(left),degree+1)):
            for j in range(min(len(right),degree+1-i)):
                out[i+j] += left[i]*right[j]
        return out
    for d in range(1,7):
        m = d-1
        @lru_cache(None)
        def weight(j):
            return d*gamma_rising(j+1,m)
        @lru_cache(None)
        def tree(j,n):
            if n==0:
                return one
            forest = [one]+[zero]*(n-1)
            for i in range(d):
                forest = convolution(forest,[tree(j+i,k) for k in range(n)],n-1)
            return weight(j)*forest[n-1]
        forest = [one]+[zero]*max_n
        for i in range(d):
            forest = convolution(forest,[tree(i,k) for k in range(max_n+1)],max_n)
        # Independently invert the exponential moment-cumulant relation.
        moment = [gamma_rising(0,d*n) for n in range(max_n+2)]
        cumulant = [zero]
        for r in range(1,max_n+2):
            cumulant.append(moment[r]-sum(
                (comb(r-1,k-1)*cumulant[k]*moment[r-k] for k in range(1,r)),zero))
        for n,u in enumerate(forest):
            quotient,remainder = sp.div(cumulant[n+1],factorial(n)*moment[1])
            assert remainder.is_zero and quotient==u
            assert u.degree()==m*n and all(c>=0 for c in u.all_coeffs())
            assert u.nth(0)>0
            leading = Fraction(d**n*comb(d*(n+1),n),n+1)
            assert leading.denominator==1 and u.LC()==leading.numerator
            assert cumulant[n+1].degree()==d*(n+1)-n
            counts["symbolic_cumulant_forest_coefficients"] += 1

        # Literal formal cross-products avoid using analytic moment functions.
        # Multiplication by n! leaves integer polynomial identities.
        for n in range(max_n+1):
            shifted_by_one = gamma_rising(1,d*n)
            shifted_by_d = gamma_rising(d,d*n)
            rhs_f = sum((factorial(n)//factorial(k)*moment[k]*tree(0,n-k)
                         for k in range(n+1)),zero)
            rhs_u = sum((factorial(n)//factorial(k)*moment[k]*forest[n-k]
                         for k in range(n+1)),zero)
            assert shifted_by_one==rhs_f and shifted_by_d==rhs_u
            counts["formal_ratio_coefficients"] += 2
        # Partial m-Dyck paths: rises +1, falls -m from height h have
        # weight d*(gamma+h-m+1)_m. They end at height m after d*n+m steps.
        if d>=2 and d<=5:
            for n in range(5):
                states = {0:one}
                for step in range(d*n+m):
                    following = {}
                    for height,value in states.items():
                        following[height+1] = following.get(height+1,zero)+value
                        if height>=m:
                            down = height-m
                            following[down] = (following.get(down,zero)
                                               +value*weight(height-m))
                    states = following
                assert states.get(m,zero)==forest[n]
                counts["independent_partial_Dyck_path_coefficients"] += 1

        for rows,cols in [([0],[1]),([0,1],[0,1]),([0,1],[1,2]),
                          ([0,2],[0,2]),([0,1,2],[0,1,2]),
                          ([0,1,2],[1,2,3])]:
            matrix = sp.Matrix([[forest[i+j].as_expr() for j in cols] for i in rows])
            minor = sp.Poly(matrix.det(method="domain-ge"),gamma)
            assert minor.is_zero or all(c>=0 for c in minor.all_coeffs())
            counts["sampled_coefficientwise_nonnegative_Hankel_minors"] += 1
        if d==3:
            bad = list(forest)
            bad[2] += one
            assert factorial(2)*moment[1]*bad[2] != cumulant[3]
            # A shifted fall weight must not pass the literal one-fall control.
            assert weight(1)+weight(2)+weight(3) != forest[1]

    counts.update({"gaps":list(range(1,7)),"maximum_forest_order":max_n,
                   "corrupted_series_and_shifted_weight_controls_rejected":True,
                   "universal_Hankel_scope":"imported PSZ Theorem9.12, ell=m=d-1; fixed gamma",
                   "two_shape_consecutive_resultant_transfer_proved":False})
    return counts

def replay():
    assert sha256((HERE/SOURCE).read_bytes()).hexdigest() == SOURCE_HASH
    assert sha256((HERE/FRACTION_SOURCE).read_bytes()).hexdigest() == FRACTION_SOURCE_HASH
    assert sha256((HERE/CUBIC_SOURCE).read_bytes()).hexdigest() == CUBIC_SOURCE_HASH
    assert sha256((HERE/HIGHER_ADJOINT_SOURCE).read_bytes()).hexdigest() == HIGHER_ADJOINT_HASH
    assert sha256((HERE/ARITHMETIC_SOURCE).read_bytes()).hexdigest() == ARITHMETIC_HASH
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
    fraction_counts = formal_cumulant_fraction_replay()
    gap_two_counts = gap_two_shape_offset_replay()
    mixed_counts = mixed_quintic_replay()
    interpolation_counts = p3_index_interpolation_replay()
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
            "gap_two_shape_offsets": "all m>=1,h2,3,gamma>2m-2; transformed adjoint Hurwitz proof",
            "gap_two_lower_exponents_two_three": "all consecutive orders and all positive Gamma shapes",
            "all_gap_two_lower_exponents_closed": False,
            "mixed_quintic_rank_equivalence": "all m>=1,gamma>0; exactC,H cyclic rank iff common moment",
            "ordinary_gap_three_lower_two": "all terminalN>=2 withN not2mod9",
            "ordinary_degree_at_most_five_equivalence": "remaining single support{2,5} atN>=11,N2mod9; infinitely many consecutive tests remain",
            "ordinary_all_quintic_binomials_closed":False,
            "p3_physical_index_interpolation_noncollision": "false: unique common zero in the stated integral Tate polydisk",
            "p3_interpolation_ordinary_integer_index_established":False,
            "formal_cumulant_branched_fraction": "all integer d>=1 and all orders; typed-tree formula over Q[gamma][[x]]",
            "normalized_cumulant_Hankel_nonnegativity": "all finite minors coefficientwise nonnegative in gamma; PSZ Theorem9.12 at ell=m=d-1",
            "branched_fraction_two_shape_resultant_transfer": False,
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
        "formal_cumulant_fraction_replay": fraction_counts,
        "gap_two_shape_offset_replay":gap_two_counts,
        "mixed_quintic_replay":mixed_counts,
        "p3_index_interpolation_replay":interpolation_counts,
        "mixed_quintic_sources":[
            {"file":HIGHER_ADJOINT_SOURCE,"sha256":HIGHER_ADJOINT_HASH,
             "scope":"gap3 order4 equation, zero concomitants and full rank3; no simplicity assumption"},
            {"file":ARITHMETIC_SOURCE,"sha256":ARITHMETIC_HASH,
             "scope":"prime-power lower-one theorem for degree<=5 classification gaps2,3,4"}
        ],
        "gap_two_adjoint_source":{"file":CUBIC_SOURCE,"sha256":CUBIC_SOURCE_HASH,
            "scope":"Theorem1 simplicity, half-plane bound4 and zero-concomitant identity10; internal elementary derivations",
            "cached_source_hash_checked":True},
        "Hankel_positivity_source": {"file":FRACTION_SOURCE,"sha256":FRACTION_SOURCE_HASH,
            "version":"arXiv1807.03271v2, May17 2020",
            "theorem":"Theorem9.12, m=d-1, ell=m; Eqs2.26/2.27 and Section2.3 product convention",
            "scope":"normalized Gamma-power cumulants at fixed shape; no resultant positivity or coprimality imported",
            "cached_primary_source_hash_checked":True},
        "representation_source": {"file": SOURCE, "sha256": SOURCE_HASH,
                                  "theorem": "Wolfs arXiv2401.08312v1 Theorem4.6 p=0,q=d",
                                  "scope": "representation attribution only; coefficient identity rederived; no zero theorem imported",
                                  "cached_primary_source_hash_checked": True},
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "source_sha256": {
            name: sha256((HERE/name).read_bytes()).hexdigest()
            for name in ["lower_exponent_initial_and_boundary.tex",
                         "lower_exponent_initial_and_boundary.py", SOURCE, TERMINAL_FIVE_CERTIFICATE, FRACTION_SOURCE, CUBIC_SOURCE, HIGHER_ADJOINT_SOURCE, ARITHMETIC_SOURCE]
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
