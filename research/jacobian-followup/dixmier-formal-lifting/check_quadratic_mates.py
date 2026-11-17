#!/usr/bin/env python3
"""Finite exact controls for QSO2; the all-order and all-curve proof is written."""
from pathlib import Path
from math import comb
import hashlib
import importlib.util
import json
import sympy as s

HERE = Path(__file__).resolve().parent
LONG = HERE.parent / "dixmier-audit" / "verify_long_pbw.py"

def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

def main():
    jets = s.symbols("u0:16")
    def dx(f, n=1):
        for _ in range(n):
            f = s.expand(sum(s.diff(f, jets[i])*jets[i+1] for i in range(len(jets)-1)))
        return f

    # R_j are built by the quadratic diagonal-resolvent identity, with
    # no differential integration and no assumed resonance factorization.
    R = [s.Rational(1, 2)]
    for m in range(1, 6):
        known = 4*sum(R[i]*R[m-i] for i in range(1, m))
        known -= 4*jets[0]*sum(R[i]*R[m-1-i] for i in range(m))
        known -= 2*sum(R[i]*dx(R[m-1-i], 2) for i in range(m))
        known += sum(dx(R[i])*dx(R[m-1-i]) for i in range(m))
        R.append(s.expand(-known/4))
    for j in range(5):
        assert s.expand(dx(R[j+1])-dx(R[j],3)/4-jets[0]*dx(R[j])-jets[1]*R[j]/2) == 0
    for j in range(1, 6):
        assert s.diff(R[j],jets[2*j-2]) == s.Rational(1,4**j)
        pure = R[j].subs({x:0 for x in jets[1:]})
        assert pure == s.binomial(2*j,j)*jets[0]**j/(2*4**j)
        for monomial in s.Poly(R[j],*jets).monoms():
            assert sum((i+2)*e for i,e in enumerate(monomial)) == 2*j
        # Formal adjoint of the universal Frechet derivative.
        T = {i:s.diff(R[j],jets[i]) for i in range(2*j-1)}
        for power in range(2*j-1):
            adj = sum((-1)**i*comb(i,power)*dx(T[i],i-power)
                      for i in range(power,2*j-1))
            assert s.expand(adj-T[power]) == 0

    a, rho = s.symbols("a rho")
    c = [s.Rational(1,2)]
    d = [s.Integer(0)]
    for j in range(6):
        c.append(s.factor(s.Rational(2*j+1,2*j+2)*(a+j*(j+1))*c[j]))
        numerator = (rho-2*j-1)*((rho-2*j)*(rho-2*j-2)/4+a)*d[j]
        numerator += (rho-4*j-2)*c[j]/2
        d.append(s.cancel(numerator/(rho-2*j-2)))
    for j in range(1,6):
        evaluation = {jet:a*(-1)**i*s.factorial(i+1) for i,jet in enumerate(jets)}
        assert s.expand(R[j].subs(evaluation)-c[j]) == 0
        actual_d = sum(s.diff(R[j],jets[i]).subs(evaluation) *
                       s.prod(rho-2-q for q in range(i))
                       for i in range(2*j-1))
        assert s.cancel(actual_d-d[j]) == 0
    bases = []
    for k in range(1,6):
        expected = s.prod((rho+2*m+1)*(rho-(2*k+4+2*m)) for m in range(k))/4**(k+1)
        assert s.cancel(d[k+1].subs(a,-k*(k+1))-expected) == 0
        bases.append({"k":k,"roots":list(range(-(2*k-1),0,2))+list(range(2*k+4,4*k+3,2))})
    for g in range(1,6):
        for k in range(1,g+1):
            polynomial = s.Poly(s.cancel(d[g+1].subs(a,-k*(k+1))),rho)
            _, factors = s.factor_list(polynomial.as_expr(),rho)
            assert sum(s.degree(f,rho)*e for f,e in factors) == 2*g
            for factor, multiplicity in factors:
                assert s.degree(factor,rho) == 1
                root = s.solve(factor,rho)[0]
                assert root.is_integer

    # A separate generic Ore implementation checks B1,B3,B5 and sharp Airy mates.
    ore = module(HERE/"check_curve_and_module.py","qso2_ore")
    z = s.symbols("z")
    u = s.Function("u")(z)
    def ordinary(f):
        return f.subs({jet:s.diff(u,z,i) for i,jet in enumerate(jets)})
    L = {2:s.Integer(1),0:u}
    powers = [{0:s.Integer(1)}]
    for j in range(1,4):
        powers.append(ore.mul(powers[-1],L,z))
    for j in range(3):
        pieces = []
        for ell in range(j+1):
            r = ordinary(R[ell])
            pieces.append(ore.mul({1:2*r,0:-s.diff(r,z)},powers[j-ell],z))
        B = ore.add(*pieces)
        assert max(B) == 2*j+1 and B[2*j+1] == 1
        assert ore.equal(ore.comm(L,B,z),{0:-4*s.diff(ordinary(R[j+1]),z)})
    airy = {2:s.Integer(1),0:z}
    even_orders = []
    for m in [1,2,4]:
        power = {0:s.Integer(1)}
        for _ in range(m):
            power = ore.mul(power,airy,z)
        Q = ore.add({1:s.Integer(-1)},power)
        assert ore.equal(ore.comm(airy,Q,z),{0:1})
        assert max(Q) == 2*m
        even_orders.append(2*m)

    # Match the retained full Long images, rather than a different fiber seed.
    long = module(LONG,"qso2_long")
    x,q,p,zz,B,d0,v,r,eta = s.symbols("x q p z B d0 v r eta")
    def expanded(poly):
        return sum(s.Rational(c.numerator,c.denominator) *
                   x**e[0]*q**e[1]*p**e[2]*zz**e[3] for e,c in poly.items())
    images = long.long_images()
    aa = 1-3*x*q
    inverse = {p:2*aa*d0+3*q*q*B,zz:-3*x*x*d0+(1+3*x*q)*B/2}
    localized = {x:1/v,q:(2*v-r*v*v)/3,B:-3*v*eta}
    fiber_images = []
    for index in [3,1]:
        f = s.expand(expanded(images[index]).subs(inverse,simultaneous=True))
        assert s.diff(f,d0) == 0
        fiber_images.append(s.expand(f.subs(localized,simultaneous=True)))
    cshift = 3*v-s.Rational(5,3)*r*v*v+r*r*v**3/3
    aval = eta+cshift
    expectedS = 2*v+4*aval-3*r*aval**2
    expectedT = (r*aval**3-aval**2-v*aval)/2
    assert s.expand(fiber_images[0]-expectedS) == 0
    assert s.expand(fiber_images[1]-expectedT) == 0
    assert [s.degree(f,eta) for f in fiber_images] == [2,3]
    assert s.Poly(fiber_images[0],eta).LC() == -3*r
    assert s.Poly(fiber_images[1],eta).LC() == r/2

    files = [HERE/"QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md",Path(__file__),LONG,
             HERE/"check_curve_and_module.py"]
    result = {
        "claim":"QSO2","status":"PASS",
        "GD_polynomials_built":list(range(6)),
        "Lenard_identities_checked":list(range(5)),
        "universal_Frechet_selfadjointness_samples":list(range(1,6)),
        "indicial_recurrence_matches_differential_polynomials":list(range(1,6)),
        "stationary_base_resonance_samples":bases,
        "all_branches_integer_roots_sample_g":list(range(1,6)),
        "Ore_Lax_operator_orders_checked":[1,3,5],
        "sharp_even_Airy_orders":even_orders,
        "retained_Long_source_order":["R","T","D","S"],
        "actual_Long_scalar_orders":{"S":2,"T":3},
        "actual_Long_scalar_leading_coefficients":{"S":"-3*r","T":"r/2"},
        "r_convention":"r=R, rather than QGR1 and QRP1's r=R/2",
        "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "coverage_boundary":"Finite exact symbolic controls only. The all-hierarchy resonance proof, no-algebraic-solution theorem, all-mate classification, valued order-drop argument and Long corollary use the written proof. No finite sample establishes those universal statements.",
        "independent_replay":False,"formal_verification":False,"external_review":False
    }
    (HERE/"quadratic_mate_controls.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: GD/Lenard and Frechet controls through R5; resonance samples through g5; B1,B3,B5; Airy mates; retained Long orders (2,3)")
    print("Written universal theorem and boundedness assumptions remain separate from these finite controls.")

if __name__ == "__main__":
    main()
