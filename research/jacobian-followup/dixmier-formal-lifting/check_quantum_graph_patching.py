#!/usr/bin/env python3
"""QPG1 exact central patching, inner-jet and quantum shear controls."""
from pathlib import Path
import hashlib
import itertools
import json
import math
import sympy as s
import check_projective_frame as pbw

HERE = Path(__file__).resolve().parent
x, p, t = pbw.x, pbw.p, pbw.t
h = s.Symbol("h")


def mul(f, g):
    return pbw.mul(f, g, h)


def mm(A, B):
    return s.Matrix(A.rows, B.cols,
                    lambda i, j: s.expand(sum(mul(A[i, k], B[k, j])
                                             for k in range(A.cols))))


def zero(A, label):
    assert all(s.cancel(c) == 0 for c in A), (label, A)


def scalar_zero(c, label):
    assert s.cancel(c) == 0, (label, c)


def U(c):
    return s.Matrix([[1, c], [0, 1]])


def V(c):
    return s.Matrix([[1, 0], [c, 1]])


def power(f, n):
    out = s.Integer(1)
    for _ in range(n):
        out = mul(out, f)
    return out


def t_valuation(c):
    c = s.cancel(c)
    if c == 0:
        return s.oo
    num, den = s.fraction(c)
    vn = min(term[0][0] for term in s.Poly(num, t).terms())
    vd = min(term[0][0] for term in s.Poly(den, t).terms())
    return vn-vd


def main():
    # Universal block identity tested over a genuinely noncommuting ring.
    X, P = s.symbols("X P", commutative=False)
    G = s.Matrix([[1+X*P, X], [P, 1]])
    Gi = s.Matrix([[1, -X], [-P, 1+P*X]])
    zero((G*Gi-s.eye(2)).applyfunc(s.expand), "free-ring inverse G Gi")
    zero((Gi*G-s.eye(2)).applyfunc(s.expand), "free-ring inverse Gi G")
    Z = s.zeros(2)
    def bU(B):
        return s.BlockMatrix([[s.eye(2), B], [Z, s.eye(2)]]).as_explicit()
    def bV(B):
        return s.BlockMatrix([[s.eye(2), Z], [B, s.eye(2)]]).as_explicit()
    J4 = s.BlockMatrix([[Z, -s.eye(2)], [s.eye(2), Z]]).as_explicit()
    free_error = bU(G)*bV(-Gi)*bU(G)*J4-s.diag(G, Gi)
    zero(free_error.applyfunc(s.expand), "noncommutative stabilized factorization")

    # An admitted scalar overlap unit with unbounded positive spatial tail.
    # b=1, N=5b+1=6. This is one exact constructive control; the general
    # precision bound follows from the written telescoping inequalities.
    q = 1-t*x
    g, gi = 1/(t*q), t*q
    scalar_zero(g*gi-1, "exact overlap inverse")
    b, N = 1, 6
    gn = sum(t**(j-1)*x**j for j in range(N+1))
    hn = gi
    assert t_valuation(g-gn) == N
    J = s.Matrix([[0, -1], [1, 0]])
    Hm = (U(gn)*V(-hn)*U(gn)*J).applyfunc(s.expand)
    Hmi = (J.inv()*U(-gn)*V(hn)*U(-gn)).applyfunc(s.expand)
    zero(Hm*Hmi-s.eye(2), "finite elementary inverse")
    zero(Hmi*Hm-s.eye(2), "finite elementary reverse inverse")
    Hplus = (s.diag(g, gi)*Hmi).applyfunc(s.cancel)
    minimum = min(t_valuation(c) for c in Hplus-s.eye(2))
    assert minimum >= N-5*b == 1
    # Idempotent agrees on the generic and formal patches.
    e = s.diag(1, 0)
    E = (Hm*e*Hmi).applyfunc(s.expand)
    Eformal = (Hplus.inv()*e*Hplus).applyfunc(s.cancel)
    zero(E-Eformal, "idempotent intersection identity")
    zero(E*E-E, "finite patched idempotent")
    for c in E:
        s.Poly(c, t)  # no negative parameter powers
    assert E.subs(t, 0) == e
    col = (t*Hm[:, 0]).applyfunc(s.expand)
    row = (Hmi[0, :]/t).applyfunc(s.expand)
    for c in list(col)+list(row):
        s.Poly(c, t)
    zero(row*col-s.eye(1), "integral rank-one split frame")
    zero(col*row-E, "integral rank-one image")

    # Full admitted right corner, with independent Weyl parameter h.
    rx = mm(mm(Hm, s.diag(x, 0)), Hmi)
    rp = mm(mm(Hm, s.diag(p, 0)), Hmi)
    for image in (rx, rp):
        for c in image:
            s.Poly(c, t)
        zero(mm(mm(E, image), E)-image, "finite corner admission")
    zero(mm(rp, rx)-mm(rx, rp)-h*E, "full corner Weyl relation")
    zero(rx.subs({t: 0, h: 0})-x*e, "special x graph")
    zero(rp.subs({t: 0, h: 0})-p*e, "special p graph")

    # The comparison has an unbounded representative but the generic action
    # and resulting rank-one algebraic module are free.
    formal_p = mul(mul(g, p), gi)
    scalar_zero(formal_p-(p-h*t/q), "unbounded inner representative")
    scalar_zero(mul(formal_p, x)-mul(x, formal_p)-h,
                "unbounded representative full relation")
    scalar_zero((formal_p-p)/(h*t)+1/q, "inner quantum jet threshold")
    # A different exact shear changes the first quantum jet; the written
    # unit classification proves it cannot be a B- or C-unit conjugation.
    shifted = p+h*x
    scalar_zero(mul(shifted, x)-mul(x, shifted)-h,
                "first quantum jet outer shear")
    assert t_valuation((shifted-p)/h) == 0

    # Polarization gives finite sums of powers over characteristic zero.
    linear = [x, x, p]
    polarized = 0
    for bits in itertools.product((0, 1), repeat=3):
        ell = sum(c*z for c, z in zip(bits, linear))
        polarized += (-1)**(3-sum(bits))*ell**3
    scalar_zero(polarized/math.factorial(3)-x*x*p,
                "repeated-form Hamiltonian polarization")
    # Exact mixed-position/momentum shear. Powers are Weyl powers, not
    # ordinary PBW symbols of a mixed linear form.
    ell = x+p
    correction = 6*t**2*power(ell, 3)
    sx, sp = x+correction, p-correction
    scalar_zero(sx+sp-ell, "mixed linear form remains fixed")
    scalar_zero(mul(sp, sx)-mul(sx, sp)-h, "exact mixed quantum shear")
    # The fixed linear form makes subtraction the exact inverse.
    scalar_zero(sx-correction-x, "mixed shear inverse x")
    scalar_zero(sp+correction-p, "mixed shear inverse p")

    inputs = [
        "CENTRAL_PATCHING_AND_QUANTUM_GRAPH_JETS.md",
        "check_quantum_graph_patching.py",
        "check_projective_frame.py",
        "INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "PROOF.md",
        "PROJECTIVE_FRAME_BOUNDARY.md",
        "FORMAL_GRAPH_BOUNDARY_MODEL.md",
    ]
    receipt = {
        "claim": "QPG1", "status": "PASS",
        "noncommutative_patching": {
            "free_ring_stabilized_identity": "PASS",
            "precision_control": {"pole_bound": b, "prefix_order": N,
                                  "minimum_formal_unit_contact": int(minimum)},
            "finite_idempotent": "PASS",
            "idempotent_spatial_degree": max(int(s.Poly(c, x).degree()) for c in E),
            "idempotent_has_no_t_poles": True,
            "two_sided_integral_rank_one_split": "PASS",
            "full_corner_Weyl_relation": "PASS",
            "special_idempotent": "e11",
            "overlap_unit": "t^-1*(1-t*x)^-1",
            "unbounded_positive_spatial_tail": True},
        "formal_graph_jets": {
            "inner_control": "p -> p-h*t/(1-t*x)",
            "inner_difference_divisible_by": "h*t",
            "noninner_first_quantum_control": "p -> p+h*x",
            "exact_mixed_linear_quantum_shear": "PASS",
            "repeated_form_polarization": "PASS",
            "all_orders_finite_jet_theorem": "written full-frame H1 proof",
            "all_ramification_orders": "written proof, no e census"},
        "sources": [
            {"url": "https://doi.org/10.3390/math14050861",
             "used": "Karoubi square patching context; the needed case is proved directly"},
            {"url": "https://stacks.math.columbia.edu/tag/0BNI",
             "used": "Commutative Beauville-Laszlo comparison only"}],
        "proof_boundary": (
            "Exact controls verify universal elementary identities and one full "
            "patched corner from an unbounded overlap frame. The general "
            "noncommutative patching, precision estimate, generic-to-integral "
            "rank-one freeness, inner finite-jet classification and finite shear "
            "reduction use written proofs. Exact overlap admission and C-freeness "
            "are required, not inferred from finite jets or division-ring frames. "
            "No finite generic Long action or Dixmier counterexample is constructed."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs}}
    (HERE/"quantum_graph_patching_controls.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QPG1", "status": "PASS",
                      "patch_prefix_order": N, "finite_corner_relation": "PASS",
                      "inner_jet_threshold": "t*a(t)", "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
