"""Exact bounded replay for BINOMIAL_IRREDUCIBILITY_REDUCTION.tex.

The all-parameter results are proved in the note. This script checks the
algebraic identities, literal moment coefficients, valuation inequalities,
and generic Eisenstein conditions on the declared finite domains. It does
not infer all-order irreducibility from a factorization search.
"""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb, factorial, gcd
from pathlib import Path
import json
import platform
import sympy as sp

HERE = Path(__file__).resolve().parent
TEX = "BINOMIAL_IRREDUCIBILITY_REDUCTION.tex"
PY = Path(__file__).name
t, x, n, b = sp.symbols("t x n b")

def rising(z, k):
    value = 1
    for j in range(k):
        value *= z + j
    return value

def polynomial(j, a, d, beta):
    return sp.Poly(sum(comb(j,k)*rising(a*j+beta,d*k)*t**k
                       for k in range(j+1)), t, domain=sp.QQ)

def valuation(z, p):
    assert isinstance(z,int) and z != 0
    count = 0
    while z % p == 0:
        count += 1
        z //= p
    return count

def identities():
    shift_cases = convex_cases = numerator_cases = 0
    for a,d in product(range(5),range(1,6)):
        r = a+d
        aa,bb = rising(x,r),rising(x,a)*rising(x,d)
        fd = rising(x+d,r)
        assert sp.expand(rising(x,d)*fd-aa*rising(x+r,d)) == 0
        assert sp.expand(aa*rising(x,d)-bb*rising(x+a,d)) == 0
        shift_cases += 1
        s = sp.symbols("s")
        second = sp.Poly(sp.diff(rising(x+s,r),s,2),x,s)
        assert all(v >= 0 for v in second.coeffs())
        convex_cases += 1
        cc = (n+1)*rising(x+a,d)-n*rising(x,d)
        numerator = rising(x+d*n,r)-rising(x,a)*cc
        assert sp.expand(numerator-(rising(x+d*n,r)-(n+1)*aa+n*bb)) == 0
        numerator_cases += 1
    clear = (5*sp.prod(5*n+j for j in range(1,5))
             -8*(2*n+1)*(n+1)*(2*n+3)*(4*n+5))
    assert sp.expand(clear-3*n*(3*n+2)*(9*n+7)*(37*n+17)) == 0
    return {"rising_shift_identity_cases":shift_cases,
            "convex_second_derivative_cases":convex_cases,
            "nondivisibility_numerator_identity_cases":numerator_cases,
            "quintic_symbolic_factor_identity":True}

def literal_coefficients():
    count = gcd_count = 0
    shapes = [sp.Rational(1,4),sp.Rational(1),sp.Rational(7,3)]
    for a,d,j,beta in product(range(4),range(1,5),range(1,7),shapes):
        f,g = polynomial(j,a,d,beta),polynomial(j+1,a,d,beta)
        shape = a*j+beta
        c = g.nth(1)-f.nth(1)
        ell = g.LC()/f.LC()
        assert c == (j+1)*rising(shape+a,d)-j*rising(shape,d)
        assert ell == rising(shape+d*j,a+d)/rising(shape,a)
        assert ell-c > 0
        count += 1
        if j<=5 and beta in shapes[:2]:
            assert sp.gcd(f,g).degree() == 0
            gcd_count += 1
    f3 = polynomial(3,0,1,sp.Rational(1,4))
    assert sp.expand(f3.as_expr()-(3*t+4)*(15*t*t+16)/64) == 0
    root = sp.Rational(-4,3)
    assert f3.eval(root) == 0
    assert polynomial(2,0,1,sp.Rational(1,4)).eval(root) == sp.Rational(8,9)
    assert polynomial(4,0,1,sp.Rational(1,4)).eval(root) == sp.Rational(32,9)
    assert polynomial(2,2,3,sp.Rational(1)).eval(sp.Rational(-1,60)) == 36
    return {"literal_positive_rational_shape_cases":count,
            "literal_consecutive_gcd_controls":gcd_count,
            "positive_algebraic_shape_reducibility_control":True,
            "quintic_first_pair_control":True}

def eisenstein():
    cases = coefficient_cases = 0
    for a,d,j in product(range(4),range(1,5),range(1,6)):
        coeff = [sp.expand(comb(j,k)*rising(a*j+b,d*k))
                 for k in range(j+1)]
        assert coeff[0] == 1
        for c in coeff[1:]:
            assert c.subs(b,-a*j) == 0
            coefficient_cases += 1
        assert sp.diff(coeff[-1],b).subs(b,-a*j) == factorial(d*j-1)
        cases += 1
    return {"generic_eisenstein_cases":cases,
            "nonleading_coefficient_divisibility_controls":coefficient_cases}

@lru_cache(None)
def fact(z):
    return factorial(z)

def prime_degree_cases():
    cases = []
    for p in [2,3,5,7]:
        for r in [1,2]:
            for d in range(1,p):
                for cofactor in range(1,p):
                    for a in range(p):
                        if (a*cofactor)%p+d*cofactor<p:
                            cases.append((a,d,p,r,cofactor))
    cases += [(100,3,7,1,1),(100,3,7,2,1),(37,1,3,2,1),
              (2,3,7,3,1),(0,1,2,3,1),(1,1,3,3,1),
              (1,2,5,3,1),(2,3,29,1,5),(2,3,73,1,2)]
    cases = list(dict.fromkeys(cases))
    coefficients = strict = chord_equalities = 0
    prime_power_orders = 0
    for a,d,p,r,cofactor in cases:
        pp = p**r
        j = cofactor*pp
        assert (a*cofactor)%p+d*cofactor<p
        assert d*cofactor < p
        ss = d*(pp-1)//(p-1)
        assert gcd(ss,pp) == 1
        assert pp*pp>d*j
        if cofactor==1:
            prime_power_orders += 1
        endpoint = Fraction(cofactor*ss,1)
        for k in range(j+1):
            literal_rf = fact(a*j+d*k)//fact(a*j)
            literal = comb(j,k)*literal_rf
            vrf, vc = valuation(literal_rf,p), valuation(literal,p)
            expected = sum((d*k)//p**h for h in range(1,r+1))
            assert vrf == expected == valuation(fact(d*k),p)
            chord = Fraction(ss*k,pp)
            assert vc >= chord
            if k in [0,j] or k%pp==0:
                assert vc == chord
                chord_equalities += 1
            else:
                assert vc > chord
                assert valuation(comb(j,k),p) == r-valuation(k,p)
                strict += 1
            coefficients += 1
        assert valuation(fact((a+d)*j)//fact(a*j),p) == endpoint
    # Dropping the no-carry condition can invalidate the exact endpoint formula.
    actual = valuation(fact(25)//fact(10),5)
    proposed = 3*(5-1)//(5-1)
    assert actual == 4 and proposed == 3 and actual != proposed
    # Literal application with neither index itself a prime power.
    assert 145 == 5*29 and 146 == 2*73
    assert (2*5)%29+3*5 == 25 < 29
    assert (2*2)%73+3*2 == 10 < 73
    assert 29*73 > 145
    assert sp.gcd(polynomial(145,2,3,sp.Rational(1)),
                  polynomial(146,2,3,sp.Rational(1))).degree() == 0
    return {"no_carry_factor_degree_cases":len(cases),
            "prime_power_order_cases":prime_power_orders,
            "literal_coefficient_valuation_controls":coefficients,
            "strict_interior_chord_controls":strict,
            "on_chord_coefficient_controls":chord_equalities,
            "largest_checked_order":max(cofactor*p**r for a,d,p,r,cofactor in cases),
            "large_lower_exponent_cases":[{"a":a,"d":d,"p":p,"r":r,"cofactor":cofactor}
                for a,d,p,r,cofactor in cases if a>=p],
            "no_carry_hypothesis_failure_control":{"a":2,"d":3,"p":5,"order":5,
                "actual_endpoint_valuation":actual,"unjustified_formula":proposed},
            "adjacent_nonprimepower_orders_control":[145,146]}

def replay():
    controls = {}
    controls.update(identities())
    controls.update(literal_coefficients())
    controls.update(eisenstein())
    controls.update(prime_degree_cases())
    return {"kind":"internal bounded exact replay; universal proofs are the TeX source",
            "scope":{
                "consecutive_nondivisibility":"all integers a>=0,d>=1,n>=1 and real beta>0",
                "generic_shape_irreducibility":"over Q(b), b an indeterminate; all a,d,n in stated integer domains",
                "transcendental_shape_coprimality":"all complex beta transcendental over Q and all orders",
                "ordinary_factor_degrees":"n=b*p^r with (a*b mod p)+d*b<p; all factor degrees divisible by p^r",
                "ordinary_primepower_irreducibility":"b=1 in the no-carry theorem",
                "ordinary_adjacent_parts":"both consecutive orders admit the stated no-carry representations",
                "all_ordinary_orders_irreducible":False,
                "all_ordinary_binomial_sfc":False,
                "bounded_controls_are_not_an_all_index_factorization_proof":True},
            "controls":controls,
            "python":platform.python_version(),"sympy":sp.__version__,
            "source_sha256":{name:sha256((HERE/name).read_bytes()).hexdigest()
                             for name in [TEX,PY]},
            "independent_replay":False,"formal_verification":False,
            "external_review":False,"passed":True}

if __name__=="__main__":
    report=replay()
    report["command"]="python research/factorial/image-rigidity/"+PY
    dest=HERE/"binomial_irreducibility_reduction_verification.json"
    dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
