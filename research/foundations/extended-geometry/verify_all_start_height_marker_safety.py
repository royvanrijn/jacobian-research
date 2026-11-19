#!/usr/bin/env python3
"""Exact controls for conditional all-start-height marker safety; stdlib only.

The unrestricted theorem is in ALL_START_HEIGHT_MARKER_SAFETY.md.
Finite moment controls do not replace its all-order purity hypothesis.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import json
import sys

if not __debug__:
    raise SystemExit("Run without -O: assertion controls must remain enabled.")

HERE = Path(__file__).resolve().parent

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p

def add(p, q):
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p): r[i] += c
    for i, c in enumerate(q): r[i] += c
    return trim(r)

def neg(p): return [-c for c in p]
def sub(p, q): return add(p, neg(q))
def scale(p, c): return trim([c*t for t in p])

def mul(p, q):
    if not p or not q: return []
    r = [0] * (len(p)+len(q)-1)
    for i, c in enumerate(p):
        for j, d in enumerate(q): r[i+j] += c*d
    return trim(r)

def power(p, n):
    r = [1]
    while n:
        if n & 1: r = mul(r, p)
        p = mul(p, p)
        n >>= 1
    return r

def product(items):
    r = [1]
    for p in items: r = mul(r, p)
    return r

def permutation_sign(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1,len(p))) % 2 else 1

def sylvester_identity(inputs):
    a2 = [mul([3,2],[5,2]), scale(mul([1,2],[5,2]),2), mul([1,2],[3,2])]
    linear = [[2*j+1,4] for j in range(5)]
    a4 = [scale(product(linear[:j]+linear[j+1:]),comb(4,j)) for j in range(5)]
    matrix = []
    for shift in range(4):
        matrix.append([[]]*shift + a2 + [[]]*(3-shift))
    for shift in range(2):
        matrix.append([[]]*shift + a4 + [[]]*(1-shift))
    det = []
    surviving = 0
    for p in permutations(range(6)):
        term = product(matrix[i][p[i]] for i in range(6))
        if term: surviving += 1
        det = add(det, scale(term, permutation_sign(p)))
    expected = scale(inputs["sylvester_positive_coefficients_ascending"],4096)
    assert det == expected, (det, expected)
    assert all(c > 0 for c in inputs["sylvester_positive_coefficients_ascending"])
    wrong = expected[:]
    wrong[0] += 1
    assert det != wrong
    return {"polynomial_degree":len(det)-1,"nonzero_permutation_terms":surviving,
            "determinant_coefficients_ascending":det,
            "positive_for_every_integer_a_ge_0":True,
            "corrupted_constant_rejected":True}

# Generic sparse multivariate polynomial arithmetic, independent of rational formulas.
def sp_clean(p): return {k:v for k,v in p.items() if v}
def sp_const(c,n): return {} if not c else {(0,)*n:c}
def sp_var(i,n):
    e = [0]*n; e[i] = 1
    return {tuple(e):1}
def sp_add(p,q):
    r = dict(p)
    for k,v in q.items(): r[k] = r.get(k,0)+v
    return sp_clean(r)
def sp_scale(p,c): return sp_clean({k:v*c for k,v in p.items()})
def sp_sub(p,q): return sp_add(p,sp_scale(q,-1))
def sp_mul(p,q):
    r = {}
    for a,c in p.items():
        for b,d in q.items():
            e = tuple(x+y for x,y in zip(a,b))
            r[e] = r.get(e,0)+c*d
    return sp_clean(r)
def sp_pow(p,n):
    r = sp_const(1,len(next(iter(p)))) if p else {}
    if n == 0: return r
    while n:
        if n&1: r = sp_mul(r,p)
        p = sp_mul(p,p); n >>= 1
    return r

def two_marker_identities():
    x,y,A,B,beta = [sp_var(i,5) for i in range(5)]
    two = sp_const(2,5)
    P = sp_add(sp_mul(x,A),sp_mul(y,B))
    R = sp_add(sp_mul(sp_pow(x,2),A),sp_mul(sp_pow(y,2),B))
    xy = sp_mul(x,y)
    # Without normalizing the weights, the recurrence has mass A+B.
    residual = sp_add(sp_sub(R,sp_mul(sp_add(x,y),P)),sp_mul(xy,sp_add(A,B)))
    assert not residual
    # Solving the normalized recurrence for y is exact after clearing denominators.
    Qformal = sp_sub(sp_mul(sp_add(x,y),P),sp_scale(xy,2))
    assert not sp_sub(sp_mul(y,sp_sub(P,sp_scale(x,2))),
                      sp_sub(Qformal,sp_mul(P,x)))
    # b0=-beta*t0 then Q=-(beta*P+2)*t0; verify algebra separately.
    t,p,b = [sp_var(i,3) for i in range(3)]
    q0 = sp_sub(sp_mul(sp_scale(sp_mul(b,t),-1),p),sp_scale(t,2))
    assert not sp_add(q0,sp_mul(sp_add(sp_mul(b,p),sp_const(2,3)),t))
    # j=-x/(1+beta*x): numerator of x+j+beta*x*j is zero.
    assert not sp_sub(sp_mul(x,sp_add(sp_const(1,5),sp_mul(beta,x))),
                      sp_add(x,sp_mul(beta,sp_pow(x,2))))
    # Negative control: changing total mass 2 to 1 changes the recurrence.
    assert sp_sub(Qformal,sp_sub(sp_mul(sp_add(x,y),P),xy))
    return {"weighted_recurrence_polynomial_identity":True,
            "rational_return_cleared_denominator_identity":True,
            "quadratic_field_reduction_identity":True,
            "mass_one_corruption_rejected":True}

def binom_integer(n,k):
    if k < 0: return 0
    ans = 1
    for i in range(k): ans = ans*(n-i)//(i+1)
    return ans

def kappa(a,m):
    ans = Q(2**m*factorial(m),1)
    for q in range(m+1): ans /= a*m+2*q+1
    return ans

def raw_unit(a,r,lam,m,j=0,k=0,ell=0):
    """Direct CT_x interval integration; no endpoint substitution."""
    target = m-j+k
    if target < 0: return Q(0)
    ans = Q(0)
    for i in range(m+1):
        exponent = m*(a*r+1)+2*r*i
        coeff = binom_integer(exponent,target)*lam**target
        for h in range(k+1):
            e = a*m+2*i+2*h+ell
            if e % 2 == 0:
                ans += Q(comb(m,i)*(-1)**i*comb(k,h)*(-1)**h,2**k*(e+1))*coeff
    return ans

def endpoint_unit(a,r,lam,m,j):
    if (a*m)%2 or m<j: return Q(0)
    return kappa(a,m)*binom_integer(r-1-j,m-j)*(-lam)**(m-j)

def raw_polynomial(C,D,a,m,j=0,k=0,ell=0):
    target = m-j+k
    if target < 0: return Q(0)
    ans = Q(0)
    for i in range(m+1):
        p = mul(power(C,m-i),power(D,i))
        coeff = p[target] if target<len(p) else 0
        for h in range(k+1):
            e = a*m+2*i+2*h+ell
            if e%2 == 0:
                ans += Q(comb(m,i)*comb(k,h)*(-1)**h,2**k*(e+1))*coeff
    return ans

def moment_controls(inputs):
    raw_count = sharp_count = negative_count = 0
    starts = inputs["starting_heights"]
    lambdas = [Q(s) for s in inputs["lambda_values"]]
    for a in starts:
        for lam in lambdas:
            for r in [1,2]:
                for m in range(1,inputs["maximum_moment"]+1):
                    for j in range(4):
                        actual = raw_unit(a,r,lam,m,j)
                        expected = endpoint_unit(a,r,lam,m,j)
                        assert actual == expected, (a,r,str(lam),m,j,actual,expected)
                        raw_count += 1
                        if r == 2 and a%2:
                            if j <= 1: assert actual == 0
                            if j == 2 and m%2 == 0:
                                assert actual == kappa(a,m)*lam**(m-2) != 0
                                sharp_count += 1
            # r=2 is deliberately impure for every even starting height.
            if a%2 == 0:
                assert raw_unit(a,2,lam,1) != 0
                negative_count += 1
    safe_count = safe_nonzero = 0
    for a in [0,1,2,3,5,6]:
        for j in range(3):
            for k in range(4):
                for ell in range(3):
                    for m in range(1,7):
                        c = raw_polynomial([1],[-1],a,m,j,k,ell)
                        if m > j-k: assert c == 0
                        elif c: safe_nonzero += 1
                        # C=x^2(1+x), D=x^2(2-x): minimum phase weight 1.
                        z = raw_polynomial([0,0,1,1],[0,0,2,-1],a,m,j,k,ell)
                        if m > k-j: assert z == 0
                        elif z: safe_nonzero += 1
                        safe_count += 2
    assert safe_nonzero > 0
    return {"raw_interval_vs_endpoint_identities":raw_count,
            "sharp_x_squared_nonzero_controls":sharp_count,
            "even_start_impurity_negative_controls":negative_count,
            "safe_weight_interval_controls":safe_count,
            "nonzero_before_safe_threshold_controls":safe_nonzero,
            "finite_controls_do_not_test_all_order_purity":True}

def scalar_contraction(p,j=0):
    value = 0
    for e,c in p.items():
        a,b,c1,d = e[0]+j,e[1],e[2],e[3]+j
        if a == c1 and b == d:
            value += c*factorial(a)*factorial(b)
    return value

def image_contraction(p):
    """The original map xi^alpha z^beta -> partial_z^alpha(z^beta)."""
    out = {}
    for (a,b,c,d),v in p.items():
        if c >= a and d >= b:
            e = (c-a,d-b)
            coefficient = v*factorial(c)//factorial(c-a)*factorial(d)//factorial(d-b)
            out[e] = out.get(e,0)+coefficient
    return sp_clean(out)

def eval_sparse(p,values):
    return sum(c*__import__("functools").reduce(lambda a,b:a*b,
                (v**e for v,e in zip(values,exponent)),1) for exponent,c in p.items())

def original_controls(inputs):
    x1,x2,z1,z2 = [sp_var(i,4) for i in range(4)]
    R = sp_add(sp_mul(x1,z1),sp_mul(x2,z2))
    Z = sp_mul(x1,z2)
    W = sp_scale(sp_mul(x2,z1),2)
    T = sp_sub(sp_mul(x1,z1),sp_mul(x2,z2))
    assert not sp_sub(sp_add(sp_pow(T,2),sp_scale(sp_mul(Z,W),2)),sp_pow(R,2))
    # A real negative control: E_2(z1)=z1 has nonzero negative weight.
    assert image_contraction(z1) == {(1,0):1}
    output_interval_checks = early_nonzero = 0
    marker_exponents = [(0,0,1,0),(1,0,0,0),(1,0,0,1),
                        (1,0,0,2),(0,0,2,1),(0,1,3,0)]
    for a in range(4):
        negative_phase = sp_scale(sp_mul(sp_pow(T,a),W),2)
        positive_phase = sp_mul(sp_pow(T,a),Z)
        for marker_e in marker_exponents:
            marker = {marker_e:1}
            A,B = sum(marker_e[:2]),sum(marker_e[2:])
            marker_weight = marker_e[0]-marker_e[2]
            for m in range(1,6):
                for direction,phase in [(-1,negative_phase),(1,positive_phase)]:
                    output = image_contraction(sp_mul(marker,sp_pow(phase,m)))
                    if A > B: assert not output
                    else:
                        assert all(sum(e)==B-A and -(B-A)<=-e[0]<=0 for e in output)
                        input_weight = marker_weight+direction*m
                        if input_weight < -(B-A) or input_weight > 0:
                            assert not output
                        elif output: early_nonzero += 1
                    output_interval_checks += 1
    assert early_nonzero > 0
    records = []
    for a in inputs["original_starting_heights"]:
        lam = 1
        U = sp_add(R,sp_scale(Z,lam))
        quotient = {}
        for j in range(1,5):
            quotient = sp_add(quotient,sp_scale(sp_mul(sp_pow(R,4-j),sp_pow(Z,j-1)),comb(4,j)*lam**j))
        assert not sp_sub(sp_mul(Z,quotient),sp_sub(sp_pow(U,4),sp_pow(R,4)))
        B = sp_sub(sp_scale(sp_mul(W,sp_pow(R,4)),2),sp_mul(sp_pow(T,2),quotient))
        F = sp_mul(sp_mul(sp_pow(T,a),sp_pow(U,2*a+1)),B)
        d = 3*a+6
        assert all(e[0]+e[1] == d and e[2]+e[3] == d for e in F)
        assert eval_sparse(F,[1,1,1,-1]) == -(2**(a+2))*lam**(2*a+5)
        assert scalar_contraction(F) == 0
        assert scalar_contraction(F,1) == 0
        F2 = sp_mul(F,F)
        assert scalar_contraction(F2) == 0
        assert scalar_contraction(F2,1) == 0
        actual = scalar_contraction(F2,2)
        expected = factorial(2*d+3)*kappa(a,2)
        assert actual == expected != 0
        assert actual != factorial(2*d+2)*kappa(a,2)  # wrong radial factorial
        records.append({"a":a,"balanced_degree":d,"ordinary_total_degree":2*d,
                        "F_monomials":len(F),"F_squared_monomials":len(F2),
                        "Z_squared_F_squared_contraction":str(actual),
                        "normalized_value":str(Q(actual,factorial(2*d+3)))})
    # Exact U_h-order residue, using T^2=R^2-2ZW and U_h=R+lambda Z.
    # After multiplying by Z the residue is R^6, hence at R=-lambda Z is lambda^6 Z^6.
    r,z,w = [sp_var(i,3) for i in range(3)]
    r6 = sp_pow(r,6)
    residue_times_z = sp_add(sp_scale(sp_mul(sp_pow(r,4),sp_mul(z,w)),2),
                    sp_mul(sp_pow(r,4),sp_sub(sp_pow(r,2),sp_scale(sp_mul(z,w),2))))
    assert residue_times_z == r6
    return {"fixtures":records,"quadratic_invariant_identity":True,
            "literal_divisibility_by_Z":True,"exact_U_residue_identity":True,
            "wrong_radial_factorial_rejected":True,
            "arbitrary_marker_output_interval_controls":output_interval_checks,
            "nonzero_inside_output_interval_controls":early_nonzero,
            "E2_z1_equals_z1_negative_weight_control":True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root",type=Path,default=Path.cwd())
    parser.add_argument("--output",type=Path,default=HERE/"all_start_height_marker_safety_verification.json")
    args = parser.parse_args()
    inputs_path = HERE/"all_start_height_marker_safety_inputs.json"
    ledger_path = HERE/"all_start_height_marker_safety_source_ledger.json"
    inputs = json.loads(inputs_path.read_text(encoding="utf-8-sig"))
    ledger = json.loads(ledger_path.read_text(encoding="utf-8-sig"))
    source_checks = []
    for entry in ledger["sources"]:
        path = (args.repo_root/entry["path"]) if entry["location"] == "repository" else (HERE/entry["path"])
        got = digest(path)
        assert got == entry["sha256"], (str(path),got,entry["sha256"])
        source_checks.append({"id":entry["id"],"sha256":got,"role":entry["role"]})
    checks = {
        "two_marker_rational_return":two_marker_identities(),
        "all_height_sylvester":sylvester_identity(inputs),
        "mixed_moment_controls":moment_controls(inputs),
        "original_contraction_controls":original_controls(inputs),
    }
    packet_files = ["ALL_START_HEIGHT_MARKER_SAFETY.md",Path(__file__).name,
                    inputs_path.name,ledger_path.name]
    receipt = {"status":"PASS","theorem":"AllStartHeightMarkerSafety",
        "scope":"Conditional on all-order purity; every integer a>=0 and every characteristic-zero field. Two original monomial markers Z,Z^2 suffice; even a needs Z only.",
        "assurance":"Exact polynomial identities and bounded rational controls supplement the written unrestricted proof; not formal verification or external human review.",
        "runtime":{"python":sys.version.split()[0],"standard_library_only":True},
        "source_checks":source_checks,"checks":checks,
        "packet_sha256":{name:digest(HERE/name) for name in packet_files}}
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","raw_moments":checks["mixed_moment_controls"]["raw_interval_vs_endpoint_identities"],
                     "original_lifts":len(checks["original_contraction_controls"]["fixtures"]),
                     "output":str(args.output)},sort_keys=True))

if __name__ == "__main__":
    main()

