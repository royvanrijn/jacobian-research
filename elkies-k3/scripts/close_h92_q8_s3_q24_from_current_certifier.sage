#!/usr/bin/env sage -python
"""
Close the H92 q8 -> q24 lattice/frame gap using the CURRENT locally passing
equation-level q8 NS certifier as the sole source of F8eq and S3std.

No q6-section constants are copied here.

Steps:
  1. execute the current certify_h92_q8_equation_ns_divisor.sage unchanged;
  2. choose the actual equation-level zero II*_E8_1;
  3. compute AJ(S3) in the resulting deterministic D13 frame;
  4. compare its pole profile with the direct modular 52/48, Z=24 result;
  5. form the primitive isotropic divisor D = O + P - k F8;
  6. classify its orthogonal child and require D12 / MW5;
  7. optionally bind the exact characteristic-zero Hensel section if present.

This deliberately does NOT use the stale physical-q8 or pinned-D13 q24 profile
to derive the answer.
"""

import json
import sys
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix,
    lcm, matrix, pari, vector
)


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


ROOT = locate_repo()
LOCAL = ROOT / "artifacts/local/elkies-k3"
CERT = ROOT / "elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
OUTPUT = LOCAL / "q8-s3-q24-equation-frame-closure.json"

if not CERT.exists():
    raise SystemExit(f"missing current q8 certifier: {CERT}")

# Execute the current local certifier exactly as-is, but isolate its argparse.
saved_argv = list(sys.argv)
scope = {"__name__": "__embedded_q8_equation_certifier__"}
try:
    sys.argv = [str(CERT)]
    exec(compile(CERT.read_text(), str(CERT), "exec"), scope)
finally:
    sys.argv = saved_argv

required = (
    "ns", "F8eq", "S3std", "selected", "isotropic_mate", "roots_and_data"
)
missing = [name for name in required if name not in scope]
if missing:
    raise SystemExit(
        "current q8 certifier passed but does not expose required variables: "
        + ",".join(missing)
    )

ns = scope["ns"]
F8eq = vector(ZZ, scope["F8eq"])
S3std = vector(ZZ, scope["S3std"])
selected = scope["selected"]
isotropic_mate = scope["isotropic_mate"]
roots_and_data = scope["roots_and_data"]

assert F8eq * ns * F8eq == 0
assert S3std * ns * S3std == -2
assert F8eq * ns * S3std == 52
assert roots_and_data(
    -(matrix(
        ZZ,
        [
            list(F8eq * ns),
            list(isotropic_mate(ns, F8eq) * ns),
        ],
    ).right_kernel_matrix()
      * ns
      * matrix(
          ZZ,
          [
              list(F8eq * ns),
              list(isotropic_mate(ns, F8eq) * ns),
          ],
      ).right_kernel_matrix().transpose())
)[2][0] == 13

print(
    "Q8S3CLOSE_SOURCE|current_local_q8_certifier=PASS|"
    "S3_degree=52|root=D13|status=PASS",
    flush=True,
)


def child_frame_with_zero(ns, fibre, zero):
    assert fibre * ns * fibre == 0
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    assert mate * ns * mate == 0
    assert mate * ns * fibre == 1
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in complement.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    U2 = matrix(ZZ, ((0, 1), (1, 0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(
            ZZ,
            [(i + 1) ** 2 + shift * (i + 1) + 1 for i in range(gram.nrows())],
        )
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [root for root in roots if regular * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root
        for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == rank
    return simple, simple * gram * simple.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4), invariants

    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i, i]) for i in range(13)) == (1,) * 13
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height0 = tail - coupling.transpose() * cartan.inverse() * coupling

    scale = ZZ(1)
    for value in height0.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale * height0).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ, 13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose() * root.inverse() * coupling
    assert height.det() == 237
    return basis, adapted, height


def class_order(dual):
    order = ZZ(1)
    for value in dual:
        order = lcm(order, ZZ(QQ(value).denominator()))
    return order


def d13_correction(dual, root):
    order = class_order(dual)
    expected = {ZZ(1): QQ(0), ZZ(2): QQ(1), ZZ(4): QQ(13) / 4}
    assert order in expected
    correction = expected[order]

    # Check the quadratic discriminant-form class, modulo 2Z.
    raw = QQ(dual * root * dual)
    mod2 = lambda x: QQ(x) - 2 * (QQ(x) / 2).floor()
    assert mod2(raw) == mod2(correction)
    return order, correction


def minimal_section_for_mw(adapted, z, cap=32768):
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height_gram = tail - coupling.transpose() * root.inverse() * coupling

    z = vector(ZZ, z)
    height = QQ(z * height_gram * z)

    base = vector(ZZ, [0] * 13 + list(z))
    pairings = vector(QQ, base * adapted[:, :13])
    dual = pairings * root.inverse()
    order, correction = d13_correction(dual, root)

    target_norm = height + correction
    assert target_norm in ZZ
    assert target_norm >= 4
    assert target_norm % 2 == 0

    lattice = IntegralLattice(root)
    iterator = lattice.enumerate_close_vectors(-dual)

    chosen = None
    for unused in range(cap):
        shift = vector(ZZ, next(iterator))
        candidate = base + vector(ZZ, list(shift) + [0] * 4)
        norm = ZZ(candidate * adapted * candidate)
        if norm == target_norm:
            chosen = candidate
            break
        if norm > target_norm:
            break

    assert chosen is not None, (
        "D13 CVP failed", tuple(z), height, correction, target_norm
    )

    pole = QQ(height + correction - 4) / 2
    assert pole in ZZ and pole >= 0
    pole = ZZ(pole)

    a = ZZ((target_norm - 2) / 2)
    section = vector(ZZ, [a, 1] + list(chosen))

    U2 = matrix(ZZ, ((0, 1), (1, 0)))
    child_ns = block_diagonal_matrix(U2, -adapted)
    F = vector(ZZ, [1, 0] + [0] * 17)
    O = vector(ZZ, [-1, 1] + [0] * 17)

    assert section * child_ns * section == -2
    assert section * child_ns * F == 1
    assert section * child_ns * O == pole

    return {
        "mw": z,
        "height": height,
        "class_order": order,
        "correction": correction,
        "pole": pole,
        "lift": chosen,
        "section": section,
    }


# ---------------------------------------------------------------------------
# 1. Equation-level q8 zero and D13 frame.
# ---------------------------------------------------------------------------

# Stored equation component degrees have E8 component 1 as multiplicity-one,
# and the local branch anchor identifies this actual curve as II*_E8_1.
O8 = vector(
    ZZ, selected["E8"]["simple_root_vectors_in_source_h3_ns"][0]
)
assert O8 * ns * O8 == -2
assert O8 * ns * F8eq == 1

child8, Bzero = child_frame_with_zero(ns, F8eq, O8)
assert roots_and_data(child8)[2] == (13, 312, 4)

A13, adapted, H = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt


def coords(curve):
    result = vector(QQ, curve) * Badapt.inverse()
    assert all(value in ZZ for value in result)
    return vector(ZZ, result)


# ---------------------------------------------------------------------------
# 2. Exact AJ of literal equation-level S3.
# ---------------------------------------------------------------------------

cS3 = coords(S3std)
degree = ZZ(cS3[1])
assert degree == 52

# For a degree-d divisor C, [C-dO] changes only U-plane coordinates.
# The positive-frame tail is therefore precisely the Pic^0/MW quotient.
zS3 = vector(ZZ, cS3[-4:])
profile = minimal_section_for_mw(adapted, zS3)
P = profile["section"] * Badapt

assert P * ns * P == -2
assert P * ns * F8eq == 1
assert P * ns * O8 == profile["pole"]

print(
    "Q8S3CLOSE_AJ|"
    f"mw={','.join(map(str,zS3))}|"
    f"height={profile['height']}|"
    f"class_order={profile['class_order']}|"
    f"correction={profile['correction']}|"
    f"PdotO={profile['pole']}|"
    "status=PASS_EXACT_EQUATION_FRAME_AJ",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Independent direct-Weierstrass profile comparison.
# ---------------------------------------------------------------------------

DIRECT = LOCAL / "q8-s3-direct-x-mod-100003.json"
direct_match = None
if DIRECT.exists():
    dm = json.loads(DIRECT.read_text())
    assert dm["status"] == "PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY"

    direct_z = ZZ(
        dm["weierstrass_structure"]["denominator_root_degree"]
    )
    xdeg = (
        ZZ(dm["x"]["numerator_degree"]),
        ZZ(dm["x"]["denominator_degree"]),
    )
    ydeg = (
        ZZ(dm["weierstrass_structure"]["y_abs_numerator_degree"]),
        ZZ(dm["weierstrass_structure"]["y_denominator_degree"]),
    )

    expected_x = (2 * profile["pole"] + 4, 2 * profile["pole"])
    expected_y = (3 * profile["pole"] + 6, 3 * profile["pole"])

    direct_match = (
        direct_z == profile["pole"]
        and xdeg == expected_x
        and ydeg == expected_y
    )

    print(
        "Q8S3CLOSE_DIRECT|"
        f"lattice_PdotO={profile['pole']}|direct_Z={direct_z}|"
        f"x={xdeg[0]}/{xdeg[1]}|y={ydeg[0]}/{ydeg[1]}|"
        f"expected_x={expected_x[0]}/{expected_x[1]}|"
        f"expected_y={expected_y[0]}/{expected_y[1]}|"
        f"match={int(direct_match)}|status=PASS",
        flush=True,
    )
    assert direct_match


# ---------------------------------------------------------------------------
# 4. The q24 isotropic divisor forced by this exact section.
# ---------------------------------------------------------------------------

# D = O + P - kF has
# D^2 = -4 + 2(P.O) - 4k,
# so k=(P.O-2)/2.
assert (profile["pole"] - 2) % 2 == 0
k = ZZ((profile["pole"] - 2) // 2)

D24 = O8 + P - k * F8eq
assert D24 * ns * D24 == 0
assert D24 * ns * F8eq == 2
assert gcd(tuple(D24)) == 1

mate24 = isotropic_mate(ns, D24)
orth24 = matrix(
    ZZ, [list(D24 * ns), list(mate24 * ns)]
).right_kernel_matrix()
child24 = -(orth24 * ns * orth24.transpose())
root24 = roots_and_data(child24)[2]

is_d12 = root24 == (12, 264, 4)
print(
    "Q8S3CLOSE_Q24|"
    f"twist={k}|square={D24*ns*D24}|"
    f"old_q8_degree={D24*ns*F8eq}|primitive={int(gcd(tuple(D24))==1)}|"
    f"root_data={root24[0]},{root24[1]},{root24[2]}|"
    f"D12={int(is_d12)}|MW_rank_if_rho19={17-root24[0]}|"
    f"status={'PASS' if is_d12 else 'NOT_D12'}",
    flush=True,
)
assert is_d12


# ---------------------------------------------------------------------------
# 5. Diagnostics against historical q24 and optional exact QQ(U) binding.
# ---------------------------------------------------------------------------

historical_same = None
historical = None
OLDQ24 = LOCAL / "q8-q24-effective-zero-choices.json"

if OLDQ24.exists():
    old = json.loads(OLDQ24.read_text())
    if old.get("status") == "PASS_EXACT_Q24_EFFECTIVE_ZERO_CHOICES":
        Dh = vector(ZZ, old["transport"]["q24_divisor_source_h3_ns"])
        delta = D24 - Dh
        historical_same = bool(D24 == Dh)
        historical = {
            "same_absolute_class": historical_same,
            "historical_degree_on_equation_F8": int(Dh * ns * F8eq),
            "difference_square": int(delta * ns * delta),
            "difference_vector": list(map(int, delta)),
        }
        print(
            "Q8S3CLOSE_HISTORICAL|"
            f"same_absolute_class={int(historical_same)}|"
            f"historical_degree_on_equation_F8={Dh*ns*F8eq}|"
            f"difference_square={delta*ns*delta}|"
            "status=PASS_DIAGNOSTIC",
            flush=True,
        )

exact_bound = False
EXACT = LOCAL / "q8-s3-direct-section-qq.json"

if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    assert exact["status"] == "PASS_EXACT_Q8_S3_DIRECT_SECTION"
    assert exact["verification"]["exact_weierstrass_identity"] is True

    assert ZZ(exact["profile"]["P_dot_O_from_denominator"]) == profile["pole"]
    assert exact["profile"]["x_degrees"] == [
        int(2 * profile["pole"] + 4),
        int(2 * profile["pole"]),
    ]
    assert exact["profile"]["y_degrees"] == [
        int(3 * profile["pole"] + 6),
        int(3 * profile["pole"]),
    ]

    exact_bound = True
    print(
        "Q8S3CLOSE_CHAR0|"
        f"PdotO={profile['pole']}|identity=1|"
        "status=PASS_EXACT_CHAR0_BINDING",
        flush=True,
    )


payload = {
    "schema": "elkies-k3.h92-q8-s3-q24-equation-frame-closure.v1",
    "status": "PASS_EXACT_AJ_Q24_LATTICE_CLOSURE",
    "source": {
        "current_local_equation_certifier": str(CERT.relative_to(ROOT)),
        "equation_q8_fibre_source_h3_ns": list(map(int, F8eq)),
        "equation_zero": "II*_E8_1",
        "equation_zero_source_h3_ns": list(map(int, O8)),
    },
    "AJ_S3": {
        "degree": int(degree),
        "mw_coordinates_in_equation_D13_basis": list(map(int, zS3)),
        "height": str(profile["height"]),
        "D13_class_order": int(profile["class_order"]),
        "D13_local_correction": str(profile["correction"]),
        "P_dot_O": int(profile["pole"]),
        "effective_section_source_h3_ns": list(map(int, P)),
        "direct_modular_profile_match": direct_match,
    },
    "q24": {
        "divisor_source_h3_ns": list(map(int, D24)),
        "decomposition": f"D=O+P-{k}F",
        "twist": int(k),
        "square": int(D24 * ns * D24),
        "old_q8_degree": int(D24 * ns * F8eq),
        "primitive": bool(gcd(tuple(D24)) == 1),
        "child_root_data": list(map(int, root24)),
        "child_root_lattice": "D12",
        "MW_rank_if_rho19": int(17 - root24[0]),
    },
    "historical_q24_diagnostic": historical,
    "exact_characteristic_zero_section_bound": exact_bound,
    "boundary": (
        "This closes the equation-frame lattice identification AJ_q8(S3) -> "
        "the primitive D12-producing q24 divisor. The remaining equation-level "
        "gap is to construct the explicit D12 genus-one pencil/Jacobian over QQ."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8S3CLOSE_RESULT|"
    f"height={profile['height']}|correction={profile['correction']}|"
    f"PdotO={profile['pole']}|twist={k}|"
    f"child=D12|MW=5|char0_bound={int(exact_bound)}|"
    "status=PASS_EXACT_AJ_Q24_LATTICE_CLOSURE",
    flush=True,
)
