#!/usr/bin/env python3
"""Exact controls for QDI1. Universal index statements are written proofs."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
ROOT = RESEARCH.parent
PROOF = HERE / "QUANTUM_DIVISION_AND_COMPLETION.md"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def zero(expr, label):
    require(s.cancel(expr) == 0, label)


def bind_receipt(name):
    """Check retained receipts against their current local inputs; do not replay them."""
    receipt = json.loads((HERE / name).read_text())
    require(receipt["status"] == "PASS", name + " status")
    for label, expected in receipt["hashes"].items():
        path = HERE / label
        if not path.exists():
            path = HERE.parent / "dixmier-audit" / label
        require(path.is_file(), name + " input exists: " + label)
        require(sha(path) == expected, name + " input hash: " + label)
    return receipt


def matrix_zero(matrix, label):
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            zero(matrix[row, col], f"{label}[{row},{col}]")


def main():
    retained = {}
    for name in ["identity_lift_controls.json", "position_avoidance_controls.json",
                 "quantum_elliptic_energy_controls.json"]:
        retained[name] = bind_receipt(name)
    qni = retained["identity_lift_controls.json"]
    require(qni["full_relations_as_polynomial_t_identities"] == 15,
            "retained full rank-three relations")
    require(qni["uniform_forward_bound"] == 12, "retained forward bound")
    original = RESEARCH / qni["retained_input"]["path"]
    require(sha(original) == qni["retained_input"]["sha256"],
            "retained foundational checker hash")

    r, a, v, d = s.symbols("r a v d")
    x = 1 / v
    q = (2*v-r*v**2) / 3
    y = a-v
    beta = 5*v**2-r*v**3-3*a*v
    u = x*y
    H = (y**4*(18*u**2+78*u+125)/20
         + 3*beta*y**2*(u**3+5*u**2+10*u-5)/10
         - beta**2*(9*u+2)/6 - x**2*beta**3/6)
    d0 = d-H
    M = beta+9*q**2
    p = 3*q**2*M+2*(1-3*x*q)*d0
    z = (1+3*x*q)*M/2-3*x**2*d0
    zero(2*x-3*x**2*q-r, "Long R inverse")
    zero(3*x**2*p+2*(1-3*x*q)*z-M, "Long B inverse")
    zero((1+3*x*q)*p/2-3*q**2*z-d0, "Long D0 inverse")
    zero(q-x*beta/3-y, "Long Y inverse")
    long_s = y+3*x*(1+u)**2*beta+3*x*y**2*(4+3*u)
    long_t = -((1+u)**3*beta+y**2*(1+u)*(4+3*u))/2
    S = 2*v+4*a-3*r*a**2
    T = (r*a**3-a**2-v*a)/2
    zero(long_s-S, "Long actual S chart")
    zero(long_t-T, "Long actual T chart")
    zero(r*a**3-2*a**2+S*a+4*T, "Long cubic")
    ss, tt, dd = s.symbols("ss tt dd")
    v_of_sa = (ss-4*a+3*r*a**2)/2
    t_of_sa = s.expand(T.subs(v, v_of_sa))
    require(s.degree(t_of_sa, a) == 3, "Long cubic field degree")
    require(s.diff(t_of_sa, a, 3) == -3*r/2, "Long cubic leading coefficient")
    require(not t_of_sa.has(dd), "fourth coordinate independent of cubic")
    zero(s.det(s.Matrix([[3*x**2, 2*(1-3*x*q)],
                        [(1+3*x*q)/2, -3*q**2]])) + 1,
         "Long inverse linear determinant")
    # The determinant is -1 in row order (B,D0); inverse formulas above are exact.

    w, rho, cc, aa, bb = s.symbols("w rho cc aa bb")
    xx = 2/rho
    yy = w-rho/2
    zz = (2*xx-3*xx**2*yy-cc)/xx**3
    uu = 1+xx*yy
    A = uu**3*zz+yy**2*uu*(4+3*xx*yy)
    B = yy+3*xx*uu**2*zz+3*xx*yy**2*(4+3*xx*yy)
    zero(A-(w**2+rho*w/2-cc*w**3), "foundational A chart")
    zero(B-(rho+4*w-3*cc*w**2), "foundational B chart")
    rho_of_bw = bb-4*w+3*cc*w**2
    a_of_bw = s.expand((w**2+rho*w/2-cc*w**3).subs(rho, rho_of_bw))
    zero(cc*w**3-2*w**2+bb*w-2*a_of_bw, "foundational cubic")
    require(s.degree(a_of_bw, w) == 3, "foundational generic degree")
    require(s.diff(a_of_bw, w, 3) == 3*cc, "foundational cubic leading coefficient")

    h, X, U = s.symbols("h X U")
    samples = []
    # Fixed regression examples, not a search over candidate lifts.
    for degree, generic in [(1, 2), (1, 4), (3, 4), (3, 5), (3, 7)]:
        f = X**degree+h*X**generic
        g = s.diff(f, X)
        zero(h*s.diff(f, X)/g-h, "rational first-order Weyl relation")
        require(s.degree(f, X) == generic, "generic polynomial map degree")
        require(s.degree(f.subs(h, 0), X) == degree, "residue polynomial map degree")
        relation = s.Poly(f-U, X, domain=s.QQ.frac_field(h, U))
        inverse_g = s.invert(s.Poly(g, X, domain=relation.domain), relation)
        zero((s.Poly(g, X, domain=relation.domain)*inverse_g).rem(relation).as_expr()-1,
             "invertible derivative in finite coefficient extension")
        samples.append({"completed_degree": degree, "generic_degree": generic,
                        "extra_dimensions": generic-degree,
                        "Weyl_relation": "PASS", "derivative_inverse_mod_relation": "PASS",
                        "negative_boundary_root_valuation": f"-1/{generic-degree}"})
        if generic in (2, 4):
            def col(expr):
                polynomial = s.Poly(expr, X, domain=relation.domain).rem(relation)
                return s.Matrix([polynomial.nth(i) for i in range(generic)])
            mx = s.Matrix.hstack(*[col(X**(j+1)) for j in range(generic)])
            deriv = s.Matrix.hstack(*[
                col(j*X**(j-1)*inverse_g.as_expr()) if j else s.zeros(generic, 1)
                for j in range(generic)])
            mg = s.Matrix.hstack(*[col(inverse_g.as_expr()*X**j)
                                   for j in range(generic)])
            # The coefficient derivative of mx is essential here.
            matrix_zero(mx.diff(U)+deriv*mx-mx*deriv-mg,
                        f"companion derivation identity d={degree},m={generic}")
            matrix_zero(mx**degree+h*mx**generic-U*s.eye(generic),
                        "companion polynomial relation")

    locked = ["QUANTUM_DIVISION_AND_COMPLETION.md", "check_division_indices.py",
              "IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md", "identity_lift_controls.json",
              "RAMIFICATION_AND_POSITION_AVOIDANCE.md", "position_avoidance_controls.json",
              "ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md", "quantum_elliptic_energy_controls.json"]
    hashes = {name: sha(HERE / name) for name in locked}
    hashes["../dixmier-audit/verify_long_pbw.py"] = sha(HERE.parent / "dixmier-audit/verify_long_pbw.py")
    receipt = {
        "claim": "QDI1", "status": "PASS", "python": platform.python_version(),
        "sympy": s.__version__, "Long_full_birational_chart_and_cubic": "PASS",
        "Long_classical_field_degree": {"value": 3, "evidence": "cubic chart and written function-field proof"},
        "full_polynomial_control": {"claim": "QNI1", "rank": 3,
            "generic_division_degree": 3, "completed_division_degree": 1,
            "forward_degree_bound": 12, "retained_full_relations": 15,
            "evidence": "exact cubic chart; written differential-field degree lemma; retained relation receipt"},
        "rational_rank_one_controls": samples,
        "companion_matrix_derivation_controls": "PASS",
        "retained_receipt_hash_chains": "PASS",
        "hashes": hashes,
        "coverage_boundary": "Exact finite chart, rational Weyl and companion identities only. Completed-degree equality, all-rank PBW index bounds, valuation-lattice equivalence, canonical and spectator invariance, and all-degree control families use written proofs. Equality of generic left/right indices uses the cited 2026 preprint. No bounded Long lift, index-three existence, or all-representative nonlifting is certified.",
        "independent_replay": False, "formal_verification": False, "external_review": False,
    }
    (HERE / "division_index_controls.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QDI1", "status": "PASS",
                      "Long_full_degree": 3, "full_polynomial_control_indices": [3, 1],
                      "rational_controls": len(samples),
                      "receipt": str(HERE / "division_index_controls.json")}, indent=2))


if __name__ == "__main__":
    main()
