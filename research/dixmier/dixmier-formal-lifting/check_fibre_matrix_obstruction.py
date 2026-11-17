#!/usr/bin/env python3
"""QFR1 exact cusp residue and exhaustive homogeneous correction controls."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
u, v = s.symbols("u v")


def unit(i, j):
    E = s.zeros(3)
    E[i, j] = 1
    return E


I = s.eye(3)
W = s.Matrix([[0, 0, -v/4], [1, 0, 3*u/4], [0, 1, 0]])
X = s.expand((4*W**2-u*I)/3)
Y = s.expand((2*u*W-v*I)/3)
D = s.expand(I-27*(Y.diff(v)*X.diff(u)-X.diff(v)*Y.diff(u)))


def L(A, B):
    return (B*X-X*B+Y*A-A*Y).applyfunc(s.expand)


def zero(M, label):
    assert all(s.cancel(z) == 0 for z in M), (label, M)


def residue(M):
    return s.expand(M[1, 1]).subs({u: 0, v: 0})


def homogeneous_basis(degree):
    """Every monomial matrix of this operator weight; no spatial cutoff."""
    out = []
    for row in range(3):
        for col in range(3):
            need = degree+col-row
            for b in range(max(0, need//3+1)):
                remainder = need-3*b
                if remainder >= 0 and remainder % 2 == 0:
                    M = s.zeros(3)
                    M[row, col] = u**(remainder//2)*v**b
                    out.append((M, row, col, remainder//2, b))
    return out


def expected_kernel(weight):
    if weight in (-5, -4):
        return []
    if weight == -3:
        return [(2*(unit(0, 1)+unit(1, 2)), u*unit(0, 2)),
                (s.zeros(3), I)]
    if weight == -2:
        return [(I, s.zeros(3)),
                (2*unit(0, 0)+2*unit(1, 1)+3*u*unit(0, 2)/2,
                 u*unit(0, 1)),
                (-2*unit(0, 0), u*unit(1, 2)),
                (3*unit(0, 0)/2,
                 -v*unit(0, 2)/4+unit(1, 0)+unit(2, 1))]
    raise AssertionError(weight)


def homogeneous_controls():
    records = []
    kernels = {}
    for weight in (-5, -4, -3, -2):
        aa = homogeneous_basis(2+weight)
        bb = homogeneous_basis(3+weight)
        columns = [L(z[0], s.zeros(3)) for z in aa]
        columns += [L(s.zeros(3), z[0]) for z in bb]
        coefficient_rows = []
        for i in range(3):
            for j in range(3):
                monomials = set()
                for M in columns:
                    monomials.update(s.Poly(M[i, j], u, v).monoms())
                for powers in sorted(monomials):
                    coefficient_rows.append((i, j, powers))
        system = s.Matrix([
            [s.Poly(M[i, j], u, v).coeff_monomial(u**powers[0]*v**powers[1])
             for M in columns] for i, j, powers in coefficient_rows])
        rank = system.rank()
        expected = expected_kernel(weight)
        assert len(expected) == len(columns)-rank
        encoded = []
        for A, B in expected:
            zero(L(A, B), "complete kernel basis")
            coefficients = [
                s.Poly(M[i, j], u, v).coeff_monomial(u**a*v**b)
                for M, basis in [(A, aa), (B, bb)]
                for _, i, j, a, b in basis]
            vector = s.Matrix(coefficients)
            assert system*vector == s.zeros(system.rows, 1)
            encoded.append(vector)
        if encoded:
            assert s.Matrix.hstack(*encoded).rank() == len(expected)
        pivot_columns = list(system.rref()[1])
        pivot_rows = list(system[:, pivot_columns].T.rref()[1])
        minor = system.extract(pivot_rows, pivot_columns)
        determinant = minor.det()
        assert determinant != 0 and minor.rows == rank
        records.append({
            "deformation_weight": weight, "unknowns": len(columns),
            "rank": rank, "kernel_dimension": len(expected),
            "original_system_rows": system.rows,
            "minor_rows": pivot_rows, "minor_columns": pivot_columns,
            "minor_determinant": str(determinant),
            "system_sha256": hashlib.sha256(str(system.tolist()).encode()).hexdigest(),
            "basis_complete": True})
        kernels[weight] = expected

    pairings = []
    for i, (A, B) in enumerate(kernels[-3]):
        for j, (C, E) in enumerate(kernels[-2]):
            value = residue(B*C-C*B+E*A-A*E)
            assert value == 0
            pairings.append({"weights": [-3, -2], "basis_indices": [i, j],
                             "period": 0})

    # Finite-copy trace check with generic 2-by-2, noncommuting multiplicity
    # matrices. The proof uses cyclic trace and applies to every finite size.
    copied_A = s.zeros(6)
    copied_B = s.zeros(6)
    for i, (A, B) in enumerate(kernels[-3]+kernels[-2]):
        variables = s.symbols("k"+str(i)+"_0:4")
        K = s.Matrix(2, 2, variables)
        copied_A += s.kronecker_product(A, K)
        copied_B += s.kronecker_product(B, K)
    comm = copied_B*copied_A-copied_A*copied_B
    middle_trace = sum(comm[i, i] for i in (2, 3))
    assert s.expand(middle_trace).subs({u: 0, v: 0}) == 0
    return records, pairings


def main():
    zero(Y*X-X*Y, "classical commuting physical matrices")
    zero(X.subs({u: 0, v: 0})-4*unit(2, 0)/3, "cusp X")
    zero(Y.subs({u: 0, v: 0}), "cusp Y")
    zero(D.subs({u: 0, v: 0})-s.diag(-8, 1, 7), "full cusp defect")
    assert residue(D) == 1
    for i in range(3):
        for j in range(3):
            assert residue(L(unit(i, j), s.zeros(3))) == 0
            assert residue(L(s.zeros(3), unit(i, j))) == 0
    # These 18 S-generator checks prove residue annihilation of the whole
    # S-linear correction map, not just corrections in a degree box.
    Au = 12*unit(0, 1)+s.Rational(21, 2)*unit(1, 2)
    Av = -s.Rational(9, 2)*unit(0, 0)
    Bv = -s.Rational(21, 4)*v*unit(0, 2)-3*unit(1, 0)
    zero(L(Au, s.zeros(3))-u*D, "u annihilator identity")
    zero(L(Av, Bv)-v*D, "v annihilator identity")
    zero(L(Au/u, s.zeros(3))-D, "localized first correction")
    b = (u**3-v**2)/27
    Ab, Bb = (u**2*Au-v*Av)/(27*b), -v*Bv/(27*b)
    zero(L(Ab, Bb)-D, "smooth-energy localized first correction")

    # Changing normal order to symmetric order changes the defect by L,
    # so the same first-order obstruction class and cusp period persist.
    symD = I-s.Rational(27, 2)*(
        Y.diff(v)*X.diff(u)-Y.diff(u)*X.diff(v)
        -X.diff(v)*Y.diff(u)+X.diff(u)*Y.diff(v))
    gaugeA = -s.Rational(27, 2)*X.diff(u).diff(v)
    gaugeB = -s.Rational(27, 2)*Y.diff(u).diff(v)
    zero(symD-D-L(gaugeA, gaugeB), "ordering coboundary")
    assert residue(symD) == 1

    weight_records, pairings = homogeneous_controls()
    assert [z["rank"] for z in weight_records] == [1, 3, 4, 5]
    assert [z["kernel_dimension"] for z in weight_records] == [0, 0, 2, 4]

    # Persistent-block compressed-commutator identity, with independent q
    # denoting t^e. It is an identity in q, not a census of e.
    q = s.Symbol("q")
    a = s.symbols("a0:9")
    c = s.symbols("c0:9")
    # Every off-middle-block entry is divisible by q; the middle block
    # is scalar, and its commutator cancels exactly.
    A = s.Matrix(3, 3, a)
    B = s.Matrix(3, 3, c)
    for i in range(3):
        for j in range(3):
            if (i == 1) != (j == 1):
                A[i, j] *= q
                B[i, j] *= q
    middle = s.expand((B*A-A*B)[1, 1])
    assert middle.coeff(q, 0) == middle.coeff(q, 1) == 0

    sx, sw = s.symbols("sx sw")
    normalization = s.Matrix([4*sw**2-3*sx, 8*sw**3-9*sx*sw, sx, sx*sw])
    jac = normalization.jacobian([sx, sw])
    assert jac.subs({sx: 0, sw: 0}).rank() == 1
    assert jac.subs({sx: 1, sw: 0}).rank() == 2

    inputs = [
        "FIBRE_RESIDUES_AND_MATRIX_CUSP_OBSTRUCTION.md",
        "check_fibre_matrix_obstruction.py",
        "FINITE_MATRIX_ACTIONS_AND_PRIMITIVE_ADMISSION.md",
        "ELLIPTIC_ISOGENY_MATRIX_MODULE.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "PROOF.md",
    ]
    out = {
        "claim": "QFR1", "status": "PASS",
        "cusp_residue": {
            "classical_pair_commutes": True,
            "projector": "E11", "target_cusp": ["u=0", "v=0"],
            "defect_diagonal": [-8, 1, 7], "period": 1,
            "full_correction_module_generators_checked": 18,
            "generated_obstruction_annihilator": "(u,v)",
            "annihilator_identities": ["u*D=L(Au,0)", "v*D=L(Av,Bv)"],
            "localized_first_order_solutions": ["u^-1", "b^-1"],
            "normal_to_symmetric_ordering_coboundary": "PASS"},
        "quadratic_lower_corrections": {
            "exhaustive_weight_systems": weight_records,
            "all_relevant_pairings": pairings,
            "generic_two_copy_trace": "PASS",
            "no_polynomial_degree_cutoff": True,
            "quantum_parameter_orders_excluded_unconditionally": [1, 2]},
        "all_parameter_persistent_block": {
            "hypothesis": "P(t) idempotent and commutes with evaluated pair modulo t^e",
            "generic_q_compressed_commutator_identity": "PASS",
            "automatic_for_all_ramification": False,
            "higher_orders_without_persistence_excluded": False},
        "normalization_point": {"origin_differential_rank": 1, "generic_rank": 2},
        "research_context": [
            {"source": "https://arxiv.org/html/1403.3493v2",
             "used": "Smooth Lagrangian line-bundle quantization does not "
                     "automatically apply to this nonimmersive normalization."},
            {"source": "https://arxiv.org/html/2604.16636v1",
             "used": "First-order Hochschild obstruction methodology is "
                     "comparison context, not an identification of classes."},
            {"source": "https://arxiv.org/html/2607.12910v1",
             "used": "Bundle and individual module-morphism lifting is "
                     "distinct from a simultaneous multiplicative Weyl action."}],
        "coverage_boundary": (
            "The first-order residue annihilates the entire S-linear "
            "correction module, and the annihilator is proved exactly "
            "by two polynomial identities and the period. The quadratic "
            "result uses every possible relevant homogeneous component; "
            "its finite matrices are exhaustive by weights. Ordinary "
            "nonlifting, two-fold quantum ramification, finite copies, "
            "identity spectators, actual canonical transport and "
            "regular spatial formal corrections use the written proofs. "
            "For arbitrary quantum ramification a persistent reducing "
            "block is an additional hypothesis; no automatic block "
            "lifting or all-ramified exclusion is asserted. The "
            "ordinary rank-one forward Long graph has zero compressed "
            "Poisson defect and is not excluded. No actual Long carrier "
            "or DC2/JC2 solution is claimed."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs},
    }
    (HERE/"fibre_matrix_obstruction_controls.json").write_text(
        json.dumps(out, indent=2)+"\n")
    print(json.dumps({"claim": "QFR1", "status": "PASS",
                      "residue_period": 1, "annihilator": "(u,v)",
                      "weight_kernel_dimensions": [0, 0, 2, 4],
                      "quadratic_pairings": len(pairings), "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
