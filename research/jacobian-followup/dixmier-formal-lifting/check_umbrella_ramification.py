#!/usr/bin/env python3
"""QUR1 exact identities for the all-ramification umbrella obstruction."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
q, z, u, v, x, w, y, a = s.symbols("q z u v x w y a")


def zero(M, label):
    assert all(s.cancel(c) == 0 for c in M), (label, M)


def comm(A, B):
    return A*B-B*A


def curvature(B1, B2):
    return a*(B2.diff(q)-B1.diff(z))+comm(B1, B2)


def main():
    # Ambient order u,v,x,y. Right target multiplication reverses its bracket.
    J = s.Matrix([[0, 27, 0, 0], [-27, 0, 0, 0],
                  [0, 0, 0, -1], [0, 0, 1, 0]])
    T = s.Matrix([
        [s.Rational(1, 4), 0, s.Rational(3, 4), 0],
        [-s.Rational(2, 15), 0, s.Rational(3, 5), 0],
        [0, -s.Rational(4, 45), 0, s.Rational(8, 15)],
        [0, s.Rational(1, 9), 0, 1]])
    standard = s.Matrix([[0, 0, -1, 0], [0, 0, 0, -1],
                         [1, 0, 0, 0], [0, 1, 0, 0]])
    zero(T*J*T.T-standard, "exact ambient quantum canonical change")
    assert T.det() == s.Rational(1, 27)
    normalization = s.Matrix([4*w**2-3*x, 8*w**3-9*x*w, x, x*w])
    transformed = (T*normalization).subs(x, z+8*w**2/15)
    expected = s.Matrix([w**2, z, 4*z*w/3, 8*w**3/9])
    zero(transformed-expected, "open umbrella polarization")
    scaling = s.diag(s.Rational(1, 2), s.Rational(8, 3),
                     2, s.Rational(3, 8))
    zero(scaling*standard*scaling.T-standard, "standard umbrella scaling")
    zero(scaling*expected-s.Matrix([w**2/2, 8*z/3,
                                   8*z*w/3, w**3/3]),
         "Bogaevski-Ishikawa standard parametrization")

    # Both old and new finite-module descriptions are identified exactly.
    I3 = s.eye(3)
    W3 = s.Matrix([[0, 0, -v/4], [1, 0, 3*u/4], [0, 1, 0]])
    X = (4*W3**2-u*I3)/3
    Y = (2*u*W3-v*I3)/3
    zero(W3**3-3*u*W3/4+v*I3/4, "rank-three monic relation")
    zero(Y-X*W3, "identified special physical module")
    zero((u*I3+3*X)/4-W3**2, "new base q on old basis")
    zero(3*X/5-2*u*I3/15-(X-8*W3**2/15), "new base z")
    W = s.Matrix([[0, q], [1, 0]])
    I2 = s.eye(2)
    zero(W**2-q*I2, "rank-two basis 1,w")
    uu, xx = s.Rational(12, 5)*q-3*z, z+s.Rational(8, 15)*q
    VV = (s.Rational(16, 5)*q-9*z)*W
    YY = xx*W
    zero((uu*I2+3*xx*I2)/4-q*I2, "old variables on new basis")
    zero(-4*VV/45+8*YY/15-4*z*W/3, "new p momentum")
    zero(VV/9+YY-8*q*W/9, "new r momentum")

    # Generic traceless pair: these equations generate the two needed
    # regular-sequence syzygies. The divisibility argument is in the proof.
    aa, bb, ff, cc, dd, gg = s.symbols("aa bb ff cc dd gg")
    A = s.Matrix([[aa, bb], [ff, -aa]])
    B = s.Matrix([[cc, dd], [gg, -cc]])
    AB = comm(A, B)
    assert s.expand(AB[1, 0]-2*(ff*cc-gg*aa)) == 0
    assert s.expand(AB[0, 0]-(bb*gg-dd*ff)) == 0
    # Independent a denotes the whole quantum parameter, not a chosen jet.
    cr, cs = s.symbols("cr cs")
    C = s.Matrix([[cr, cs], [1, -cr]])
    uv = s.symbols("uv0:4")
    vv = s.symbols("vv0:4")
    U, V = s.Matrix(2, 2, uv), s.Matrix(2, 2, vv)
    full = comm(ff*C+a*U, gg*C+a*V)
    claimed = a*(ff*comm(C, V)-gg*comm(C, U))+a*a*comm(U, V)
    zero(full-claimed, "universal moving-zero factorization identity")
    zero(full.subs({ff: 0, gg: 0})-a*a*comm(U, V),
         "commutator divisible by a squared at the moving zero")

    alpha, beta = s.Rational(4, 3), s.Rational(8, 9)
    B1, B2 = alpha*z*W, beta*q*W
    F = B2.diff(q)-B1.diff(z)
    zero(F-s.Matrix([[0, 4*q/9], [-s.Rational(4, 9), 0]]),
         "nonzero special derivative matrix")
    assert F[1, 0] == -s.Rational(4, 9)
    f0, g0 = B1[1, 0], B2[1, 0]
    jacobian = s.Matrix([f0, g0]).jacobian([q, z])
    assert jacobian.det() == -s.Rational(32, 27)
    assert s.diff(g0, q)-s.diff(f0, z) == -s.Rational(4, 9)

    # Exact full localized positive controls, valid for symbolic a.
    K = s.diag(s.Rational(1, 6), -s.Rational(1, 6))
    zero(curvature(B1, B2+a*K/z), "z-local full flat connection")
    zero(curvature(B1+a*s.diag(0, 1/(2*q)), B2),
         "q-local full flat connection")
    N = s.Matrix([[0, 0], [1, 0]])
    zero(curvature(z*N, q*N), "closed scalar one-form positive control")

    # A genuine commuting deformation for which the old projector cannot
    # persist. This is not a test of any quantum ramification order.
    sh = s.Symbol("sh")
    Xt, Yt = X.subs(v, v-6*sh), Y.subs(v, v-6*sh)-2*sh*I3
    zero(comm(Xt, Yt), "moving-cusp commuting deformation")
    Xc = Xt.subs({u: 0, v: 0})
    Yc = Yt.subs({u: 0, v: 0})
    P = s.diag(0, 1, 0)
    zero(Yc, "deformed old-cusp Y")
    pv = s.symbols("pv0:9")
    P1 = s.Matrix(3, 3, pv)
    linear = comm(P, Xc.diff(sh))+comm(P1, Xc.subs(sh, 0))
    assert linear[0, 1] == -2 and linear[1, 2] == 2

    # The correspondence ring misses w and contains every other monomial:
    # w-even powers use w^2; odd powers >=3 use w^3, and x*w is present.
    # This exhaustive monomial argument is written, not sampled here.
    jac = normalization.jacobian([x, w])
    assert jac.subs({x: 0, w: 0}).rank() == 1
    assert jac.subs({x: 1, w: 0}).rank() == 2

    inputs = [
        "REPOLARIZATION_AND_UMBRELLA_RAMIFICATION.md",
        "check_umbrella_ramification.py",
        "FIBRE_RESIDUES_AND_MATRIX_CUSP_OBSTRUCTION.md",
        "FINITE_MATRIX_ACTIONS_AND_PRIMITIVE_ADMISSION.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "PROOF.md",
    ]
    out = {
        "claim": "QUR1", "status": "PASS",
        "canonical_change": {"matrix": s.sstr(T), "determinant": "1/27",
                             "all_six_brackets": "PASS",
                             "right_target_sign": -27},
        "normalization": {
            "image_ring": "k[x,w^2,xw,w^3]", "missing_monomial": "w",
            "conductor": "(x,w^2)", "old_left_rank": 3,
            "new_position_rank": 2, "standard_umbrella_scaling": "PASS",
            "singular_differential_rank": 1},
        "rank_two_obstruction": {
            "traceless_syzygy_equations": "PASS",
            "scalar_factor_jacobian": "-32/27",
            "curvature_lower_left": "-4/9",
            "universal_factorization": "PASS",
            "moving_zero_commutator_divisibility": "a^2",
            "derivative_valuation": "e",
            "commutator_minimum_valuation": "2e",
            "all_e_proof": "regular sequence modulo a and formal moving zero",
            "persistent_block_required": False,
            "spatial_degree_cutoff": None, "quantum_order_census": False},
        "positive_controls": {
            "q_local_full_flat_connection": "PASS",
            "z_local_full_flat_connection": "PASS",
            "zero_derivative_difference": "PASS"},
        "nonpersistent_projector": {
            "exact_commuting_deformation": "PASS",
            "unavoidable_linear_entries": [-2, 2]},
        "source": {
            "url": "https://msp.org/pjm/2002/203-1/pjm-v203-n1-p06-p.pdf",
            "use": "standard open Whitney umbrella parametrization only"},
        "proof_boundary": (
            "All quantum parameter orders for the single prescribed monic "
            "normalization module use the written regular-sequence, moving-zero "
            "and completed-repolarization proof. Exact controls verify the "
            "supporting universal identities. The higher-multiplicity case "
            "and arbitrary added summands are not covered by this rank-two "
            "argument. Localized connections allow poles and do not repair "
            "regular admission. The scalar forward Long graph is different "
            "and remains formally unobstructed. No actual finite Long carrier "
            "or Dixmier counterexample is claimed."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs}}
    (HERE/"umbrella_ramification_controls.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps({"claim": "QUR1", "status": "PASS",
                      "quantum_orders": "all e>=1",
                      "curvature_entry": "-4/9", "new_position_rank": 2,
                      "localized_full_connections": 2, "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
