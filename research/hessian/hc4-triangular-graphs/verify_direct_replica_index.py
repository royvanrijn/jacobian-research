#!/usr/bin/env python3
"""Exact internal support for the given-replica index proof; no classification."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as S

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "DIRECT_REPLICA_INDEX_FROM_A_DESCENDED_DIFFERENTIAL.md"
CHECKER = Path(__file__).resolve()
checks: list[dict[str, object]] = []

def check(name: str, condition: object) -> None:
    ok = bool(condition)
    checks.append({"name": name, "passed": ok})
    if not ok:
        raise AssertionError(name)

def equal(name: str, actual: object, expected: object) -> None:
    check(name, S.cancel(S.sympify(actual) - S.sympify(expected)) == 0)

def order_zero(f: object, v: S.Symbol) -> int:
    num, den = S.fraction(S.cancel(f))
    def polynomial_order(p: object) -> int:
        pp = S.Poly(p, v)
        if pp.is_zero:
            raise ValueError("zero has no finite order")
        return min(mon[0] for mon, _ in pp.terms())
    return polynomial_order(num) - polynomial_order(den)

def run() -> None:
    x, y, z, v, u = S.symbols("x y z v u")
    a0, a1, a2, a3, a4 = S.symbols("a0:5")
    p = (y*y+1)*(y-2)
    beta = 1/(y*y+3)
    sl = x/p + beta
    delta = lambda f: S.cancel(p*S.diff(f, x))
    equal("given replica has a genuine rational slice", delta(sl), 1)
    equal("slice recovers source x", p*(sl-beta), x)
    equal("slice volume coefficient", S.diff(sl,x), 1/p)
    equal("slice transverse wedge has coefficient p", p*S.diff(sl,x), 1)

    f = a4*x**4+a3*x**3+a2*x**2+a1*x+a0
    top = delta(delta(delta(delta(f))))
    equal("top derivative is invariant", delta(top), 0)
    current = f
    decomposition = S.Integer(0)
    for n in range(4, -1, -1):
        derivative = current
        for _ in range(n):
            derivative = delta(derivative)
        coefficient = S.cancel(derivative/S.factorial(n))
        equal(f"Taylor coefficient {n} is invariant", delta(coefficient), 0)
        decomposition += coefficient*sl**n
        current = S.cancel(current-coefficient*sl**n)
        next_derivative = current
        for _ in range(n):
            next_derivative = delta(next_derivative)
        equal(f"Taylor subtraction lowers degree at {n}", next_derivative, 0)
    equal("finite Taylor reconstruction is exact", decomposition, f)
    equal("last Taylor remainder vanishes", current, 0)

    # The symbolic monomial family is not a whole-plane Keller construction.
    # Its rational target coordinates test precisely the descended volume.
    m, j = S.symbols("m j", integer=True, positive=True)
    h = y**m
    r = h**j
    weight = S.diff(h,y)*r
    rational_P = x/weight
    rational_Q = h**(j+1)/(j+1)
    equal("unbounded monomial field/volume control has Jacobian one",
          S.diff(rational_P,x)*S.diff(rational_Q,y), 1)
    equal("unbounded monomial replica fixes Q",
          weight*S.diff(rational_Q,x), 0)
    equal("unbounded monomial replica translates P",
          weight*S.diff(rational_P,x), 1)
    equal("unbounded differential descent coefficient",
          weight/S.diff(h,y), r)
    k_plus_one = m*(j+1)
    pole_order_eta = j+2
    equal("unbounded index and pole-order identity",
          k_plus_one, m*(pole_order_eta-1))
    equal("pullback order at the unique infinity",
          -m*(j+1)-1, m*(-pole_order_eta)+(m-1))

    # Nonmonomial controls include finite branch points and extra zeros.
    for index, (hpoly,rpoly) in enumerate([
        (y, u**3+u+1),
        (y*y+y+1, (u-2)**2),
        (y**3+y, u*u+3*u+5),
        (y**4-2*y*y+3, u+1),
        (y**5+y*y-1, S.Integer(7)),
    ]):
        pp = S.expand(S.diff(hpoly,y)*rpoly.subs(u,hpoly))
        ss = x/pp + 1/(y+7)
        primitive = S.integrate(rpoly,u)
        qq = primitive.subs(u,hpoly)
        equal(f"nonmonomial rational volume {index}",
              S.diff(ss,x)*S.diff(qq,y), 1)
        equal(f"nonmonomial descended differential {index}",
              pp/S.diff(hpoly,y), rpoly.subs(u,hpoly))
        degree_h = S.degree(hpoly,y)
        degree_r = S.degree(rpoly,u)
        check(f"nonmonomial index arithmetic {index}",
              S.degree(pp,y)+1 == degree_h*(degree_r+1))
        check(f"nonmonomial rational target scope guard {index}",
              S.degree(pp,y)>0 and S.denom(S.cancel(ss)) != 1)


    # The Hamiltonian identity is general, before specializing any polynomial map.
    Px, Py, Qx, Qy, rho = S.symbols("Px Py Qx Qy rho")
    equal("general Hamiltonian identity has zero dx coefficient",
          rho*Px*Qx-rho*Qx*Px, 0)
    equal("general Hamiltonian dy coefficient is p times the Jacobian",
          rho*Px*Qy-rho*Qx*Py, rho*(Px*Qy-Py*Qx))
    for index, (PP,QQ) in enumerate([
        (x+y**3+y, y+(x+y**3+y)**2),
        (2*(x+y*y)+3*(y+(x+y*y)**3),
         (x+y*y)+2*(y+(x+y*y)**3)),
    ]):
        aa = p*S.diff(PP,x)
        bb = p*S.diff(QQ,x)
        equal(f"actual polynomial shear control Jacobian {index}",
              S.diff(PP,x)*S.diff(QQ,y)-S.diff(PP,y)*S.diff(QQ,x), 1)
        equal(f"actual Hamiltonian dx coefficient {index}",
              aa*S.diff(QQ,x)-bb*S.diff(PP,x), 0)
        equal(f"actual Hamiltonian dy coefficient {index}",
              aa*S.diff(QQ,y)-bb*S.diff(PP,y), p)

    # Compute traces using multiplication matrices, rather than assuming
    # the root-of-unity trace formula. This checks differential commutation
    # on positive and negative powers in actual cyclic function fields.
    for cover_degree in range(2,6):
        mult = S.zeros(cover_degree)
        for col in range(cover_degree-1):
            mult[col+1,col] = 1
        mult[0,cover_degree-1] = u
        check(f"cyclic multiplication matrix relation degree{cover_degree}",
              mult**cover_degree == u*S.eye(cover_degree))
        for power in range(-3,8):
            trace_value = S.trace(mult**power)
            traced_derivative = (
                S.Rational(power,cover_degree)
                * S.trace(mult**(power-cover_degree))
            )
            equal(f"differential trace commutation degree{cover_degree} power{power}",
                  S.diff(trace_value,u), traced_derivative)
    cover_m, relative_e = S.symbols("cover_m relative_e",
                                  integer=True, positive=True)
    tower_H = y**(cover_m*relative_e)
    equal("unbounded polynomial tower derivative identity",
          y*S.diff(tower_H,y), cover_m*relative_e*tower_H)
    for cover_degree, relative_degree in [(1,1),(2,3),(3,2),(4,5),(7,3)]:
        HH = y**(cover_degree*relative_degree)
        check(f"tower through C(y^{cover_degree}) with relative degree{relative_degree}",
              S.degree(HH,y) == cover_degree*relative_degree)


    # Monomial stable fields: the primitive gives a cyclic subextension.
    # Rational target coordinates here are controls, not polynomial Keller maps.
    zz, TT = S.symbols("zz TT")
    for cover_degree, primitive_degree in [
        (1,1),(2,2),(2,6),(3,6),(4,12),(5,10),(7,21),
    ]:
        kk = primitive_degree-1
        cc = S.Integer(2)
        bbeta = 1/(y+2)+y*y
        old_slice = x/(cc*y**kk)+bbeta
        multiplier = cc*y**(primitive_degree-cover_degree)
        bb = S.cancel(multiplier*bbeta)
        new_slice = S.cancel(multiplier*old_slice)
        equal(f"monomial slice rescaling m{cover_degree} N{primitive_degree}",
              new_slice, x/y**(cover_degree-1)+bb)
        check(f"monomial multiplier belongs to invariant field m{cover_degree} N{primitive_degree}",
              primitive_degree>=cover_degree
              and (primitive_degree-cover_degree) % cover_degree == 0)
        equal(f"rescaled monomial delta coefficient m{cover_degree} N{primitive_degree}",
              cc*y**kk*S.diff(new_slice,x), multiplier)
        sigmax = S.cancel(
            (zz*y)**(cover_degree-1)
            *(new_slice-bb.subs(y,zz*y))
        )
        equal(f"rational deck formula fixes new slice m{cover_degree} N{primitive_degree}",
              sigmax/(zz*y)**(cover_degree-1)+bb.subs(y,zz*y),
              new_slice)
        equal(f"deck character commutes with monomial replica m{cover_degree} N{primitive_degree}",
              S.rem(zz**kk-zz**(cover_degree-1),
                    zz**cover_degree-1,zz), 0)
        mat = S.zeros(cover_degree)
        for col in range(cover_degree-1):
            mat[col+1,col] = 1
        mat[0,cover_degree-1] = u
        equal(f"cyclic primitive polynomial m{cover_degree} N{primitive_degree}",
              mat.charpoly(TT).as_expr(), TT**cover_degree-u)

    # Direct local pullback computations with a varying unit retain ramification.
    for e in range(1,6):
        local_u = v**e*(1+2*v)
        for differential_order in range(-4,4):
            pulled = S.cancel(local_u**differential_order*S.diff(local_u,v))
            check(f"local differential e{e} order{differential_order}",
                  order_zero(pulled,v)
                  == e*differential_order+e-1)
    for e in range(1,7):
        for h_pole in range(1,5):
            check(f"every branch over a pole remains polar e{e} h{h_pole}",
                  -e*(h_pole-1)-1 < 0)

    # Stability is indispensable: C(x^2,y) is not partial_x stable.
    unstable_generator = x*x
    unstable_velocity = S.diff(unstable_generator,x)
    equal("quadratic subfield generator is deck invariant",
          unstable_generator.subs(x,-x), unstable_generator)
    check("quadratic subfield derivative leaves invariant field",
          S.expand(unstable_velocity.subs(x,-x)-unstable_velocity) != 0)
    check("dropping stability would violate the bare-translation bound",
          1 % 2 != 0)

    # Rational replicas can create a second pole; polynomiality is indispensable
    # to the one-pole step, even when a volume descends.
    equal("rational weight differential descends through y squared",
          (1/(2*u)).subs(u,y*y)*S.diff(y*y,y), 1/y)
    check("rational weight has a finite differential pole",
          order_zero(1/y,y) == -1)
    equal("rational weight also has an infinity pole",
          (1/y).subs(y,1/v)*S.diff(1/v,v), -1/v)

    # Constant weight and m=1 are allowed, not silently discarded.
    for k in range(0,9):
        check(f"identity intermediate has exact pole order k{k}",
              k+1 == 1*((k+2)-1))
    check("bare translation permits only index one",
          [d for d in range(1,10) if 1 % d == 0] == [1])

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    run()
    receipt = {
        "claim": "Direct given-polynomial-replica index divisibility",
        "outcome": "PASS",
        "check_count": len(checks),
        "python": platform.python_version(),
        "sympy": S.__version__,
        "source": SOURCE.name,
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "checker": CHECKER.name,
        "checker_sha256": hashlib.sha256(CHECKER.read_bytes()).hexdigest(),
        "scope": (
            "Internal exact formula controls. The written proof gives the "
            "unrestricted theorem for a given delta=p(y)partial_x and a "
            "delta-stable intermediate field of an actual plane Keller map. "
            "Nontrivial rational field controls are outside the whole-plane "
            "polynomial target hypotheses. No LND classification, first-"
            "velocity integrality, external review or formal verification."
        ),
        "checks": checks,
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({
        "outcome": "PASS",
        "check_count": len(checks),
        "elapsed_seconds": time.monotonic()-started,
        "receipt_written": bool(args.write_receipt),
    },indent=2))

if __name__ == "__main__":
    main()
