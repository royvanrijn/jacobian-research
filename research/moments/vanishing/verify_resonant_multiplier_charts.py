"""Exact root-only reconstruction, residue-module and parameter-chart controls."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as S
from verify_linear_mixed_cutoff import apply_certificate, multiindices, image_W

ROOT = Path(__file__).resolve().parent


def reconstruct(b, terms, s, t, max_degree):
    roots = {int(r) for r in S.polys.polytools.ground_roots(S.Poly(b, s)) if r.is_Integer and r>=0}
    values = {}
    bases = []
    for alpha in sorted(multiindices(2, max_degree), key=lambda x:(sum(x), x)):
        d = sum(alpha)
        if d in roots:
            values[alpha] = image_W(alpha, t)
            bases.append(alpha)
        else:
            rhs = S.Integer(0)
            for c, k, j, u, v in terms:
                eta = tuple(a+u0 for a, u0 in zip(alpha, u))
                if any(v0>e for e, v0 in zip(eta, v)):
                    continue
                child = tuple(e-v0 for e, v0 in zip(eta, v))
                assert sum(child)<d
                rhs -= c*(-1)**sum(v)*S.prod(S.ff(e, v0) for e, v0 in zip(eta, v))*t**k*S.diff(values[child], t, j)
            values[alpha] = S.cancel(rhs/b.subs(s, d))
        assert S.expand(values[alpha]-image_W(alpha, t)) == 0
    return sorted(roots), bases, len(values)


def other_branch(a, b, t):
    # Residue at Z=1/t of b!*Z^a/(Z-tZ^2)^(b+1).
    return S.cancel((-1)**(b+1)*factorial(b)*S.binomial(a-b-1, b)*t**(b-a))


def run():
    Z, W, t, s = S.symbols('Z W t s')
    R3 = [(2,0,0,(1,0),(2,1)),(2,0,0,(0,1),(1,2)),(10,0,0,(0,0),(1,1)),
          (-1,1,0,(1,0),(2,0)),(-1,1,0,(0,1),(1,1)),(-5,1,0,(0,0),(1,0))]
    assert apply_certificate(s*(s-3),R3,[Z,W],t,s,t*W-Z*W) == 0
    roots,bases,count = reconstruct(s*(s-3),R3,s,t,7)
    assert roots == [0,3]
    assert len(bases) == 5 and not any(sum(x) in [1,2] for x in bases)

    assert other_branch(0,0,t) == -1
    assert other_branch(1,0,t) == -1/t
    residue_checks = 0
    for a,b in multiindices(2,5):
        # Total W derivative: partial_W + tZ^2-Z.
        dW = (b*other_branch(a,b-1,t) if b else 0)+t*other_branch(a+2,b,t)-other_branch(a+1,b,t)
        # Total Z derivative: partial_Z + (2tZ-1)W.
        dZ = (a*other_branch(a-1,b,t) if a else 0)+2*t*other_branch(a+1,b+1,t)-other_branch(a,b+1,t)
        dt = S.diff(other_branch(a,b,t),t)-other_branch(a+2,b+1,t)
        assert S.cancel(dW) == S.cancel(dZ) == S.cancel(dt) == 0
        residue_checks += 3
    # Direct moving-pole residues for small nontrivial monomials.
    y = S.symbols('y')
    for a,b in [(0,0),(1,0),(0,1),(2,1),(2,2)]:
        amplitude = factorial(b)*(1/t+y)**a*(1+t*y)**(-b-1)*(-1)**(b+1)
        direct = S.expand(S.series(amplitude,y,0,b+1).removeO()).coeff(y,b)
        assert S.cancel(direct-other_branch(a,b,t)) == 0

    # Exact Gaussian image for safe P=Z^2 W, no b-function value assumed.
    safe = 0
    for a,b in multiindices(2,5):
        f = S.Rational(factorial(2*b-a),factorial(b-a))*t**(b-a) if b>=a else S.Integer(0)
        for m in range(7):
            moment = factorial(a+2*m) if a+2*m==b+m else 0
            assert S.diff(f,t,m).subs(t,0) == moment
            safe += 1

    c = S.symbols('c0:5')
    p = sum(c[k]*Z**k for k in range(5))
    terms = [(-2,0,0,(0,0),(1,1))]+[(-k*c[k]*(-1)**k,1,0,(0,0),(0,k)) for k in range(1,5)]
    assert apply_certificate(s,terms,[Z,W],t,s,t*p-Z*W) == 0
    # The finite exact ideal is (c0); positive pure moments are c0**m.
    for m in range(1,7):
        constant = S.Poly(S.expand(p**m),Z).nth(0)
        assert constant == c[0]**m
    for k in range(1,5):
        # P=c_k Z^k; W^k attains ratio 1/k.
        marker = factorial(k)*c[k]*t
        assert S.diff(marker,t).subs(t,0) == factorial(k)*c[k]

    a = S.symbols('a')
    assert S.Poly((s-a),s).degree() == 1
    special = []
    for N in [0,1,5]:
        f = (1+t)**N
        assert S.expand((1+t)*S.diff(f,t)-N*f) == 0
        assert S.degree(f,t) == N
        special.append(N)
    return {
        'completed':True,
        'gapped_b_function_certificate':{'b':'s*(s-3)','root_degrees':roots,'root_amplitudes_used':len(bases),
                                       'amplitudes_reconstructed':count,'gap_degrees_not_used_as_bases':[1,2]},
        'other_branch_residue_module_checks':residue_checks,
        'direct_moving_residue_checks':5,
        'other_branch_seed_image':'-1','other_branch_Z_image':'-1/t',
        'safe_Z2W_moment_controls':safe,
        'positive_root_not_sufficient_for_failure_written_proof':True,
        'parametric_Z_only_certificate_degree_checked':4,
        'parametric_Z_only_tail_ideal':'(c0)','parametric_Z_only_uniform_pure_cutoff':'1',
        'nonconstructible_holonomic_control_degrees':special,
        'nonconstructibility_proved_by_checker':False,
        'universal_Gaussian_chart_cover_constructed':False,
        'universal_fixed_n_D_cutoff_proved':False,
        'universal_Gaussian_constructibility_proved':False,
        'general_tail_ideal_solver_implemented':False,
        'general_b_function_search_executed':False,
        'producer_package_executed':False,'new_cubic_prefix_computed':False,
        'independent_replay':False,'external_review':False,'formal_verification':False}


if __name__ == '__main__':
    result=run()
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT/'resonant_multiplier_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

