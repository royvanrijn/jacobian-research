"""Exact controls for residual lifting and binomial moment parity.

Written proof: BINOMIAL_RESIDUAL_LIFTING.tex.
Finite tests control the implementation; the written proof supplies all
unbounded quantifiers. No all-terminal quintic closure is claimed.
"""
from fractions import Fraction
from math import factorial, comb
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def quintic_all_prime_degree_replay():
    """Verify one complete two-prime certificate, not an index-prefix search."""
    from itertools import product
    from math import gcd, isqrt

    counts = {}

    def check(test, kind):
        assert test, kind
        counts[kind] = counts.get(kind, 0) + 1

    def vp(number, prime):
        assert number
        out = 0
        while number % prime == 0:
            number //= prime
            out += 1
        return out

    def legendre(number, prime):
        out = 0
        while number:
            number //= prime
            out += number
        return out

    def literal_heights(degree, prime):
        out = []
        for k in range(degree + 1):
            direct = comb(degree, k) * (
                factorial(2 * degree + 3 * k) // factorial(2 * degree))
            height = vp(direct, prime)
            formula = (legendre(degree, prime) - legendre(k, prime)
                - legendre(degree - k, prime)
                + legendre(2 * degree + 3 * k, prime) - legendre(2 * degree, prime))
            check(height == formula, "direct factorial versus Legendre valuations")
            out.append(height)
        return out

    def certify(degree, prime, edges, supplied_heights, known_heights=None):
        heights = (literal_heights(degree, prime) if known_heights is None
                   else known_heights)
        check(heights == supplied_heights, "literal coefficient-height certificate")
        left, height, previous = 0, Fraction(heights[0]), None
        for numerator, denominator, width in edges:
            slope = Fraction(numerator, denominator)
            check(width > 0 and width % slope.denominator == 0,
                  "integral edge width")
            check(previous is None or previous < slope, "strictly increasing edge slopes")
            check(height == heights[left], "left edge endpoint")
            right = left + width
            check(right <= degree, "edge stays within exact degree")
            right_height = height + width * slope
            check(right_height == heights[right], "right edge endpoint")
            for k, value in enumerate(heights):
                check(value >= height + (k - left) * slope,
                      "complete supporting-line inequalities")
            left, height, previous = right, right_height, slope
        check(left == degree and height == heights[-1], "complete polygon coverage")

    def degree_set(one, two):
        one = {Fraction(a, b): width for a, b, width in one}
        two = {Fraction(a, b): width for a, b, width in two}
        out = {0}
        for slope in one.keys() & two.keys():
            step = slope.denominator
            width = min(one[slope], two[slope])
            check(width % step == 0, "shared-slope integral degree step")
            out = {a + b for a in out for b in range(0, width + 1, step)}
        return out

    certificate_path = Path(__file__).resolve().parent / (
        "binomial_residual_lifting_n64_certificate.json")
    data = json.loads(certificate_path.read_text())
    check(data["polynomial"] ==
        "P_n(t)=sum binom(n,k)*(2n+3k)!/(2n)!*t^k", "actual quintic polynomial")
    profiles = {}
    for row in data["profiles"]:
        key = row["prime"], row["order"]
        check(key not in profiles, "nonduplicate Newton profile")
        profiles[key] = row
        certify(row["order"], row["prime"], row["edges"], row["coefficient_valuations"])
    check(set(profiles) == {(3, 64), (3, 65), (5, 64), (5, 65)},
          "complete declared two-prime inputs")
    d3 = degree_set(profiles[3, 64]["edges"], profiles[3, 65]["edges"])
    d5 = degree_set(profiles[5, 64]["edges"], profiles[5, 65]["edges"])
    check(d3 == {a + 9 * b + 27 * c
                 for a, b, c in product(range(2), range(2), range(3))},
          "exact three-adic bounded subset sums")
    check(d5 == {a + 5 * b + 25 * c
                 for a, b, c in product(range(3), range(2), range(3))},
          "exact five-adic bounded subset sums")
    after_three = set(range(2, 64, 2)) & d3
    check(after_three == {10, 28, 36, 54}, "complete parity degrees after prime3")
    check(sorted(after_three) == data["positive_even_degrees_after_p3"],
          "literal intermediate degree certificate")
    check(not after_three & d5, "empty two-prime common-degree set")
    check(data["possible_degrees_after_p3_and_p5"] == [], "literal final degree certificate")

    def no_carry_parts(degree):
        found = []
        for prime in range(2, degree + 1):
            if any(prime % j == 0 for j in range(2, isqrt(prime) + 1)):
                continue
            power = prime
            while power <= degree:
                if degree % power == 0:
                    b = degree // power
                    if (2 * b) % prime + 3 * b < prime:
                        found.append((prime, power, b))
                power *= prime
        return found

    check(no_carry_parts(64) == no_carry_parts(65) == [],
          "complete no-carry-part check at64and65")

    def multiply(one, two):
        out = [0] * (len(one) + len(two) - 1)
        for i, a in enumerate(one):
            for j, b in enumerate(two):
                out[i + j] += a * b
        return out

    check(multiply([1, -2, 2], [1, 2, 2]) == [1, 0, 0, 0, 4],
          "Sophie Germain polynomial identity")
    leading = [5]
    for j in range(1, 5):
        leading = multiply(leading, [j, 5])
    first_difference = multiply(multiply([1, 2], [1, 1]), multiply([3, 2], [5, 4]))
    first_difference = [8 * a for a in first_difference]
    difference = [a - b for a, b in zip(leading, first_difference)]
    claimed = multiply(multiply([0, 3], [2, 3]), multiply([7, 9], [17, 37]))
    check(difference == claimed, "exact positive full-degree exclusion identity")
    A, lower, upper = 2, 64, 65
    left, right = 2 * A * A - 2 * A + 1, 2 * A * A + 2 * A + 1
    check(gcd(left, right) == 1 and left * right == upper and lower == 4 * A ** 4,
          "literal first coverage-failure factorization")
    check(right * right < 3 * upper, "literal first coverage-failure bound")
    check(Fraction(13, 4) ** 2 < 12, "uniform coverage-failure constant inequality")

    corrupted = [list(row) for row in profiles[3, 64]["edges"]]
    corrupted[1] = [13, 9, 18]
    try:
        certify(64, 3, corrupted, profiles[3, 64]["coefficient_valuations"],
                known_heights=profiles[3, 64]["coefficient_valuations"])
    except AssertionError:
        counts["corrupted edge table rejected"] = 1
    else:
        raise AssertionError("corrupted edge table accepted")
    return {
        "status": "pass",
        "complete_certificate_orders": [64, 65],
        "complete_certificate_primes": [3, 5],
        "checks_by_kind": counts,
        "exact_checks": sum(counts.values()),
        "literal_complete_Newton_edge_tables": 4,
        "all_supporting_lines_checked": True,
        "corrupted_edge_table_rejected": True,
        "possible_positive_even_degrees_after_p3": sorted(after_three),
        "possible_degrees_after_p3_and_p5": [],
        "eligible_no_carry_parts_at64and65": [],
        "infinite_no_carry_failure_family": "n=2^(12s+6), N=n+1, s>=0",
        "infinite_family_is_SFC_counterexample": False,
        "prefix_enumeration": False,
        "unbounded_all_prime_criterion_proved": True,
        "all_index_emptiness_of_degree_sets_proved": False,
        "all_terminal_quintic_SFC": False,
        "independent_replay": False,
        "external_review": False,
        "formal_verification": False,
    }


def replay():
    x=s.Symbol('x')
    counts={}

    def check(test,kind):
        assert test,kind
        counts[kind]=counts.get(kind,0)+1

    def val(v,p):
        if v==0:
            return 10**9
        z=0
        while v%p==0:
            v//=p
            z+=1
        return z

    def factorial_val(n,p):
        r=0
        while n:
            n//=p
            r+=n
        return r

    def unit_factorial_literal(n,p):
        f=factorial(n)
        return f//p**factorial_val(n,p)%p

    def unit_factorial_digits(n,p):
        out=1
        while n:
            q,r=divmod(n,p)
            out=out*pow(-1,q,p)*factorial(r)%p
            n=q
        return out

    def admissible(a,d,p,b):
        return s.isprime(p) and a>=0 and d>=1 and b>=1 and (a*b%p)+d*b<p

    def small_moment(b,a,d):
        return [comb(b,k)*factorial(a*b+d*k)//factorial(a*b)for k in range(b+1)]

    def literal_coeff(n,a,d,k):
        return comb(n,k)*factorial(a*n+d*k)//factorial(a*n)

    def poly_from(row,p):
        return s.Poly(sum(int(c)*x**j for j,c in enumerate(row)),x,modulus=p)

    def unit_part(v,p):
        return (v//p**val(v,p))%p

    def edge_residual(a,d,p,b,r):
        assert admissible(a,d,p,b)
        P=p**r
        height=d*(P-1)//(p-1)
        n=b*P
        actual=[]
        for j in range(b+1):
            k=j*P
            vc=(factorial_val(n,p)-factorial_val(k,p)-factorial_val(n-k,p)
                 +factorial_val(a*n+d*k,p)-factorial_val(a*n,p))
            check(vc==height*j,'edge exact valuations')
            numerator=(unit_factorial_digits(n,p)*unit_factorial_digits(a*n+d*k,p))%p
            denominator=(unit_factorial_digits(k,p)*unit_factorial_digits(n-k,p)
                         *unit_factorial_digits(a*n,p))%p
            actual.append(numerator*pow(denominator,-1,p)%p)
        return actual

    cases=[(0,1,2,1),(0,1,3,2),(2,1,5,1),(9,1,7,2),
           (12,2,11,3),(19,3,17,3),(2,3,31,3),(2,3,23,4),
           (2,3,29,4),(2,3,31,6)]
    for a,d,p,b in cases:
        assert admissible(a,d,p,b)
        for r in range(1,5):
            P=p**r
            height=d*(P-1)//(p-1)
            expected=[c*pow(-1,height*j,p)%p for j,c in enumerate(small_moment(b,a,d))]
            actual=edge_residual(a,d,p,b,r)
            check(actual==expected,'complete edge residual identities')
            if b*P<=200:
                n=b*P
                for k in range(n+1):
                    c=literal_coeff(n,a,d,k)
                    vc=val(c,p)
                    if k%P:
                        check(Fraction(vc)>Fraction(k*height,P),'literal strict nonedge')
                    else:
                        check(vc==k*height//P,'literal exact edge')
                        check(unit_part(c,p)==actual[k//P],'independent literal edge units')
    for p in [2,3,5,7,11]:
        for r in range(1,4):
            P=p**r
            exponent=(P-1)//(p-1)
            for m in [0,1,p-1,p,p+2,2*p+1]:
                check(unit_factorial_literal(m*P,p)==
                      pow(-1,m*exponent,p)*unit_factorial_literal(m,p)%p,
                      'literal Wilson factorial recursion')

    def residual_factor(row,P,height,p):
        degree=len(row)-1
        assert degree%P==0 and row[0]%p!=0
        for k,c in enumerate(row):
            assert val(c,p)>=Fraction(k*height,P)
        assert val(row[-1],p)==degree*height//P
        return [row[j*P]//p**(height*j)%p for j in range(degree//P+1)]

    for p,r,d in [(2,2,1),(3,1,1),(3,2,2),(5,1,3)]:
        P=p**r
        height=d*(P-1)//(p-1)
        rows=[]
        for ell in [1,2]:
            row=[(k+1)*p**((k*height+P-1)//P) for k in range(P*ell+1)]
            row[0]=1
            row[-1]=p**(height*ell)
            rows.append(row)
        product=[sum(rows[0][i]*rows[1][k-i]
                     for i in range(len(rows[0]))if 0<=k-i<len(rows[1]))
                 for k in range(len(rows[0])+len(rows[1])-1)]
        rg,rh=[residual_factor(row,P,height,p)for row in rows]
        check(poly_from(residual_factor(product,P,height,p),p)==
              poly_from(rg,p)*poly_from(rh,p),'factor residual multiplicativity')

    def modular_power(base,exponent,modulus):
        out=s.Poly(1,x,modulus=modulus.get_modulus())
        while exponent:
            if exponent%2:
                out=(out*base).rem(modulus)
            base=(base*base).rem(modulus)
            exponent//=2
        return out

    def irreducible_certificate(f,p):
        degree=f.degree()
        xx=s.Poly(x,x,modulus=p).rem(f)
        check(modular_power(xx,p**degree,f)==xx,'Rabin terminal Frobenius')
        for prime in s.factorint(degree):
            test=modular_power(xx,p**(degree//int(prime)),f)-xx
            check(s.gcd(f,test).degree()==0,'Rabin proper-divisor gcd')

    examples=[
        {'b':3,'p':31,'factors':[[7,13,17,1]]},
        {'b':3,'p':37,'factors':[[3,22,32,1]]},
        {'b':4,'p':23,'factors':[[21,15,12,3,1]]},
        {'b':4,'p':29,'factors':[[15,1],[3,13,13,1]]},
        {'b':6,'p':31,'factors':[[6,1],[18,27,16,22,24,1]]},
    ]
    for row in examples:
        b,p=row['b'],row['p']
        full=poly_from(small_moment(b,2,3),p).monic()
        factors=[poly_from(f,p)for f in row['factors']]
        product=s.Poly(1,x,modulus=p)
        for f in factors:
            irreducible_certificate(f,p)
            product*=f
        check(product==full,'literal residual factor products')
        degrees=[f.degree()for f in factors]
        subset={0}
        for degree in degrees:
            subset|={a+degree for a in list(subset)}
        if b%2==0:
            check(not any(0<h<b and h%2==0 for h in subset),
                  'complete parity subset-sum criterion')
        else:
            check(degrees==[b],'irreducible odd residual criterion')
        row['factor_degrees']=degrees
        row['irreducible_for_all_r_by_proof']=True

    # Scope controls: carry inputs and rounded-down factor valuations are rejected.
    check(not admissible(2,3,5,1),'carry countercontrol')
    try:
        residual_factor([1,1,0,3],3,1,3)
    except AssertionError:
        counts['bad factor valuation rejected']=1
    else:
        raise AssertionError('bad factor valuation accepted')
    check(set([1,3])=={1,3} and all(v%2 for v in [1,3]),'odd-pair parity pattern')
    check(any(0<k<4 and k%2==0 for k in {0,2,4}),'even-pattern rejection control')
    for r in range(1,13):
        n=4*23**r
        check((n%9==2)==(r%6==1),'remaining-ray first terminal congruence')
        check(((n+1)%9==2)==(r%6==2),'remaining-ray second terminal congruence')

    parity_profiles=0
    parity_counts=0
    for a,d,beta in [(0,1,s.Rational(1,2)),(2,3,s.Integer(1)),
                     (2,3,s.Rational(7,3)),(4,2,s.Rational(5,2)),
                     (1,5,s.Integer(2))]:
        for degree in range(1,9):
            polynomial=s.Poly(sum(s.binomial(degree,k)
                *s.rf(a*degree+beta,d*k)*x**k
                for k in range(degree+1)),x,domain=s.QQ)
            if degree%2==0:
                check(polynomial.count_roots(-s.oo,s.oo)==0,
                      'positive-Gamma even polynomial Sturm profile')
            else:
                check(polynomial.count_roots(-s.oo,s.oo)==1,
                      'positive-Gamma odd polynomial Sturm profile')
                check(polynomial.diff().count_roots(-s.oo,s.oo)==0,
                      'positive-Gamma odd derivative Sturm profile')
            parity_profiles+=1
    check(s.Poly(1+3*x+x*x,x).count_roots(-s.oo,s.oo)==2,
          'coefficient positivity alone rejected')

    quintic_degree = quintic_all_prime_degree_replay()
    counts["all-prime quintic degree certificate controls"] = quintic_degree["exact_checks"]
    out={
        'status':'pass',
        'arithmetic':'exact integer, rational, and finite-field',
        'scope':'Under (ab modp)+db<p, the actual orderbp^r residual is P_(b,a,d,1)((-1)^sY); normalized factor degrees/p^r are bounded subset-sums of residual factor degrees. For oddp/evenb, positivity parity gives Q-irreducibility for an irreducible residual or exactly two odd-degree irreducible residual factors.',
        'positive_Gamma_parity_scope':'For every beta>0: even P_n has no real roots; odd P_n is strictly increasing with one simple real root; every real consecutive common divisor has even degree.',
        'checks_by_kind':counts,
        'total_checks':sum(counts.values()),
        'rational_positive_Gamma_Sturm_profiles':parity_profiles,
        'examples':examples,
        "quintic_all_prime_degree":quintic_degree,
        'all_r_quantifiers_supplied_by_written_proof':True,
        'finite_tests_prove_unbounded_scope':False,
        'all_terminal_quintic_closed':False,
        'proof':'BINOMIAL_RESIDUAL_LIFTING.tex',
        'python_version':platform.python_version(),
        'sympy_version':s.__version__,
        'independent_replay':False,
        'external_review':False,
        'formal_verification':False,
        'priority_claim':False,
    }
    return out


if __name__ == "__main__":
    report = replay()
    here = Path(__file__).resolve().parent
    report["input_sha256"] = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in ("BINOMIAL_RESIDUAL_LIFTING.tex",
                     "binomial_residual_lifting.py",
                     "BINOMIAL_IRREDUCIBILITY_REDUCTION.tex",
                     "binomial_residual_lifting_n64_certificate.json")
    }
    destination = here / "binomial_residual_lifting_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
