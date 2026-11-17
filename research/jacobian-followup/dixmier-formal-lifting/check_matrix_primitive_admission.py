#!/usr/bin/env python3
"""QRA1 exact matrix controls. Universal assertions use the written proofs."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
x, p, h = s.symbols("x p h")


def zero(M, label):
    assert all(s.expand(e) == 0 for e in M), (label, M)


def companion(coeffs):
    m = len(coeffs)
    Z = s.zeros(m)
    for j in range(m-1):
        Z[j+1, j] = 1
    for j, coefficient in enumerate(coeffs):
        Z[j, m-1] = -coefficient
    Z[0, m-1] += x
    C = s.zeros(m)
    for j in range(1, m):
        C[j-1, j] = j
    return Z, C


def monic_controls():
    out = []
    for m in (2, 3, 5):
        coeffs = s.symbols("c0:"+str(m))
        Z, C = companion(coeffs)
        derivative = m*Z**(m-1)
        for j in range(1, m):
            derivative += j*coeffs[j]*Z**(j-1)
        derivative = derivative.applyfunc(s.expand)
        zero(derivative*Z-Z*derivative, "f'(Z) commutes with Z")
        zero(derivative*Z.diff(x)+C*Z-Z*C-s.eye(m), "chain rule identity")
        T = derivative*p+h*C
        # Both operators have momentum order at most one.
        comm = T*Z-Z*T+h*T.diff(p)*Z.diff(x)
        zero(comm-h*s.eye(m), "exact companion Weyl relation")
        degZ = max(s.Poly(e, x, p).total_degree() for e in Z if e != 0)
        degT = max(s.Poly(e, x, p).total_degree() for e in T if e != 0)
        assert degZ <= 1 and degT <= 2
        powerZ, powerC = companion([0]*m)
        zero(powerZ**m-x*s.eye(m), "power companion exact relation")
        inverse = powerZ**(m-1)/x
        zero(inverse*powerZ-s.eye(m), "rational inverse")
        assert s.cancel((inverse*powerC)[m-1, 1]-1/x) == 0
        leadingW = m*powerZ**(m-2)*p
        powerT = m*powerZ**(m-1)*p+h*powerC
        zero(powerT-powerZ*leadingW-h*powerC, "primitive correction")
        out.append({"m": m, "all_monic_coefficients_symbolic": True,
                    "full_Weyl_pair": "PASS", "spatial_degrees": [int(degZ), int(degT)],
                    "power_companion_primitive_defect": "h*C",
                    "inverse_primitive_correction_entry": "1/x"})
    return out


def elliptic_controls():
    u, v, b, pu, pv = s.symbols("u v b pu pv")
    relation = v*v-u**3+27*b
    def red(z):
        numerator, denominator = s.fraction(s.cancel(z))
        return s.cancel(s.rem(numerator, relation, v)/denominator)
    def rm(M):
        return M.applyfunc(red)
    def rz(M, label):
        zero(rm(M), label)
    I = s.eye(3)
    W = s.Matrix([[0, 0, -v/4], [1, 0, 3*u/4], [0, 1, 0]])
    Cd = s.Matrix([[0, u/2, v/2], [0, 0, -u/2], [0, -1, 0]])
    Ce = s.diag(0, 1/(6*b), 1/(3*b))
    X = (4*W**2-u*I)/3
    Y = (2*u*W-v*I)/3
    delta = lambda f: v*s.diff(f, u)+3*u**2*s.diff(f, v)/2
    energy = lambda f: u*s.diff(f, u)/(3*b)+v*s.diff(f, v)/(2*b)+s.diff(f, b)
    dm = lambda M: M.applyfunc(delta)
    em = lambda M: M.applyfunc(energy)
    comm = lambda A, B: A*B-B*A
    assert red(delta(relation)) == 0 and red(energy(relation)) == 0
    for f in (u, v, b):
        assert red(energy(delta(f))-delta(energy(f))-delta(f)/(6*b)) == 0
    rz(dm(W)+comm(Cd, W)-(u*I/2-W**2), "W delta")
    rz(em(W)+comm(Ce, W)-W/(6*b), "W energy")
    rz(em(Cd)-dm(Ce)+comm(Ce, Cd)-Cd/(6*b), "mixed curvature")
    rz(dm(X)+comm(Cd, X)+Y, "X delta")
    rz(dm(Y)+comm(Cd, Y)+3*X**2/2, "Y delta")
    rz(em(X)+comm(Ce, X)-X/(3*b), "X energy")
    rz(em(Y)+comm(Ce, Y)-Y/(2*b), "Y energy")
    rz(Y*Y-X**3-b*I, "source curve")
    zero(Y-X*W, "exact admitted primitive")

    # D_x/h=A*delta+B*E+C and D_y/h=D*delta+F*E+G.
    A = -Y/b
    B = -3*X**2
    C = rm(A*Cd+B*Ce)
    D = 2*X/(3*b)
    F = 2*Y
    G = rm(D*Cd+F*Ce)
    rz(comm(X, Y), "positions commute")
    rz(A*dm(X)+B*em(X)+comm(C, X)-I, "physical dx x")
    rz(A*dm(Y)+B*em(Y)+comm(C, Y), "physical dx y")
    rz(D*dm(X)+F*em(X)+comm(G, X), "physical dy x")
    rz(D*dm(Y)+F*em(Y)+comm(G, Y)-I, "physical dy y")
    # Coefficients of delta^2,E^2,delta*E,delta,E,1 in [D_x,D_y]/h^2.
    rz(comm(A, D), "second delta")
    rz(comm(B, F), "second energy")
    rz(comm(A, F)+comm(B, D), "mixed second")
    rz(A*dm(D)+B*em(D)-D*dm(A)-F*em(A)+comm(A, G)+comm(C, D)
       +(B*D-F*A)/(6*b), "first delta")
    rz(A*dm(F)+B*em(F)-D*dm(B)-F*em(B)+comm(B, G)+comm(C, F),
       "first energy")
    rz(A*dm(G)+B*em(G)-D*dm(C)-F*em(C)+comm(C, G), "constant")
    assert energy(b) == 1
    # The mate and exact primitive force the same localized x inverse.
    forced_inverse = rm(D*dm(W)+F*em(W)+comm(G, W))
    rz(X*forced_inverse-I, "forced right inverse")
    rz(forced_inverse*X-I, "forced left inverse")
    assert red(X.det()+4*b) == 0

    # A different full, global action uses monic source coordinates x,w.
    C0 = s.Matrix([[0, 1, 0], [0, 0, 2], [0, 0, 0]])
    Pxw = -3*pu*I-9*W*pv
    Pwx = 8*W*pu+(24*W**2-9*X)*pv+h*C0
    def star(A0, B0):
        # Each left factor here is first order in the target momenta.
        return (A0*B0+h*(A0.diff(pu)*B0.diff(u)+A0.diff(pv)*B0.diff(v))).applyfunc(s.expand)
    ops = [X, W, Pxw, Pwx]
    for i in range(4):
        for j in range(i+1, 4):
            kappa = -1 if (i, j) in ((0, 2), (1, 3)) else 0
            zero(star(ops[i], ops[j])-star(ops[j], ops[i])-h*kappa*I,
                 "global monic Weyl relation")
    q = v*pu+3*u**2*pv/2
    logP = u*pu/3+v*pv/2
    logM = logP*I+h*s.diag(0, s.Rational(1, 6), s.Rational(1, 3))
    zero((X*Pxw/3+W*Pwx/6)-logM, "global logarithmic energy")
    zero(-Y*Pxw+(W**2-3*X/2)*Pwx-(q*I+h*Cd), "global curve derivation")
    rz(b*em(W)+comm(b*Ce, W)-W/6, "global logarithmic W")
    rz(b*em(Cd)-dm(b*Ce)+comm(b*Ce, Cd)-Cd/6, "global logarithmic mixed")
    assert s.expand(b*energy(b)-b) == 0
    # b=(u^3-v^2)/27 has a zero inside the target affine plane.
    bpoly = (u**3-v**2)/27
    assert bpoly.subs({u: 1, v: 1}) == 0
    assert bpoly.subs({u: 0, v: 1}) != 0
    assert X.det().subs({u: 1, v: 1}) == 0

    return {
        "free_basis": ["1", "w", "w^2"], "matrix_size": 3,
        "Euler_connection": "diag(0,1/(6*b),1/(3*b))",
        "mixed_curvature": "PASS", "physical_surface_full_relations": 6,
        "noncentral_energy_relation": "[P_B,b]=h",
        "global_monic_full_relations": 6,
        "global_monic_positions": ["x", "w"],
        "global_physical_y_relation": "[D_w,y]=h*x",
        "global_logarithmic_energy": "[L_B,b]=h*b",
        "global_logarithmic_derivation_and_mixed_relation": "PASS",
        "physical_surface_requires_b_inverse": True,
        "exact_primitive": "Y=X*W", "forced_inverse_identity": "PASS",
        "special_x_determinant": "-4*b", "interior_zero_energy_point": {"u": 1, "v": 1},
        "special_left_rank": 3, "identified_rank_one_Long_graph": False,
    }


def main():
    monic = monic_controls()
    elliptic = elliptic_controls()
    inputs = [
        "FINITE_MATRIX_ACTIONS_AND_PRIMITIVE_ADMISSION.md",
        "check_matrix_primitive_admission.py",
        "ELLIPTIC_ISOGENY_MATRIX_MODULE.md",
        "ELLIPTIC_QUANTUM_ENERGY_MOTION.md",
        "INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "CENTRALIZER_INDEX_AND_RATIONAL_POLARIZATION.md",
        "QUANTUM_DIVISION_AND_COMPLETION.md",
    ]
    out = {
        "claim": "QRA1", "status": "PASS",
        "ramified_companion_actions": monic,
        "elliptic_surface_extension": elliptic,
        "unit_specialization_boundaries": {
            "inverse_pole_removal_requires_injective_special_map": True,
            "t_times_identity": "generic inverse has a pole; special map is not injective",
            "one_minus_t_x": "special unit does not imply a polynomial generic inverse",
            "classical_factorization_alone_excluded": False,
            "arbitrary_broken_quantum_primitive_excluded": False},
        "research_context": [
            {"source": "https://stacks.math.columbia.edu/tag/00UP",
             "used": "Finite-order differential operators extend under formally "
                     "etale base change; compares the smooth isogeny case. "
                     "The ramified matrix restriction formula is proved directly."}],
        "coverage_boundary": (
            "Finite-free restriction of scalars, the uniform monic degree "
            "bound, division-regular irreducibility, inverse-pole removal and "
            "exact-primitive obstruction are written universal proofs. The "
            "companion controls use symbolic coefficients at listed degrees, "
            "not a proof by matrix-size census. The explicit surface model "
            "satisfies six physical Weyl relations after inverting b. A "
            "different polynomial monic-coordinate model satisfies six global "
            "relations in x,w and has an Euler energy relation. Neither has "
            "the identified rank-one Long graph. The exact factorization "
            "hypothesis cannot be weakened to its classical reduction, as "
            "the positive power-companion corrections demonstrate. Rational "
            "compression is excluded over the original target division ring, "
            "not after all coefficient field extensions. No bounded full "
            "Long lift or DC2/JC2 resolution is claimed."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs},
    }
    (HERE/"matrix_primitive_admission_controls.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps({"claim": "QRA1", "status": "PASS",
                      "monic_degrees": [2, 3, 5], "physical_relations": 6,
                      "global_monic_relations": 6, "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
