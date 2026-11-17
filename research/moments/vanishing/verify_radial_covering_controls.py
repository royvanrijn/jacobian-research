#!/usr/bin/env python3
"""Small exact controls; the dimension and covering proofs are written."""
from pathlib import Path
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent

def moment(poly,z,zb,rho):
    poly=S.sympify(poly)
    assert not poly.has(S.Float)
    assert S.Poly(poly,*(z+zb)).total_degree()<=16
    total=0
    for powers,c in S.Poly(S.expand(poly),*(z+zb)).terms():
        alpha=powers[:len(z)];beta=powers[len(z):]
        if alpha==beta:total+=c*rho(sum(alpha))*math.prod(math.factorial(a) for a in alpha)
    return S.simplify(total)

def main():
    a,b,aa,bb=S.symbols('a b aa bb');z=(a,b);zb=(aa,bb)
    rho=lambda m:S.Rational(1,math.factorial(m+1))
    u=(a+b)/S.sqrt(2);v=S.I*(a-b)/S.sqrt(2)
    ub=(aa+bb)/S.sqrt(2);vb=-S.I*(aa-bb)/S.sqrt(2)
    transforms=[]
    for powers in [(1,0,1,0),(0,1,0,1),(2,0,2,0),(1,1,1,1),(2,1,2,1),(2,0,1,1),(1,0,0,1)]:
        p,q,r,s=powers
        actual=moment(u**p*v**q*ub**r*vb**s,z,zb,rho)
        expected=moment(a**p*b**q*aa**r*bb**s,z,zb,rho)
        assert actual==expected
        transforms.append(dict(powers=list(powers),value=str(actual)))
    norms=[]
    for m in (1,2,3,4):
        actual=moment((a*aa+b*bb)**m,z,zb,rho)
        assert actual==1
        norms.append(dict(q=2,m=m,value=str(actual)))
    c,cc=S.symbols('c cc');z3=(a,b,c);zb3=(aa,bb,cc)
    rho3=lambda m:S.Rational(2,math.factorial(m+2))
    for m in (1,2,3):
        actual=moment((a*aa+b*bb+c*cc)**m,z3,zb3,rho3)
        assert actual==1
        norms.append(dict(q=3,m=m,value=str(actual)))
    # Liu Proposition7.1,D=6: rho_m=2*m!/(2m+2)!. Its radial
    # norm is nonconstant and therefore needs more than three source dimensions.
    projective=lambda m:S.Rational(2*math.factorial(m),math.factorial(2*m+2))
    r2=moment(a*aa+b*bb,z,zb,projective)
    r4=moment((a*aa+b*bb)**2,z,zb,projective)
    assert r2==S.Rational(1,6) and r4==S.Rational(1,30)
    variance=S.simplify(r4-r2**2);assert variance==S.Rational(1,180)
    phase=(1+bb)*(b-a*aa);circuits=[]
    for m in (1,2,3,4):
        assert moment(phase**m,z,zb,rho)==0
        for k in (1,2):
            actual=moment(bb**k*phase**m,z,zb,rho)
            expected=0 if k>m else S.Rational(math.comb(m-1,k-1),m+1)
            assert actual==expected
            circuits.append(dict(m=m,k=k,mixed=str(actual)))
    # For x->x^3 on uniform[-1,1], target mass[-eps,eps]=eps^(1/3),
    # exceeding a bounded target-density O(eps) by an unbounded factor.
    local=[]
    for scale in (2,4,8):
        epsilon=S.Rational(1,scale**3);mass=S.Rational(1,scale)
        assert mass/epsilon==scale**2
        local.append(dict(epsilon=str(epsilon),source_mass=str(mass),ratio=str(mass/epsilon)))
    guards=0
    for poly in [S.Float(1),a**17]:
        try:moment(poly,z,zb,rho)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['verify_seed_density_and_trace.py','verify_gaussian_radial_transfer.py']},
        exact_unitary_moment_controls=transforms,sphere_norm_controls=norms,
        projective_norm_second_moment=str(r2),projective_norm_fourth_moment=str(r4),projective_norm_variance=str(variance),
        source_circuit_controls=circuits,critical_map_scaling_controls=local,guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        covering_theorem_proved_by_samples=False,all_degree_exclusion_proved_by_samples=False,
        all_phase_invariant_circuits_excluded=False,actual_conjugacy_hypothesis_removed=False,
        source_sparse_minimality_certificate_replayed=False,full_exceptional_MZ_decided=False,
        gaussian_radial_bridge_invalidated=False,new_frontier_endpoint_claimed=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (HERE/'radial_covering_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['completed','checker_sha256','projective_norm_variance','guard_refusals']}))

if __name__=='__main__':main()
