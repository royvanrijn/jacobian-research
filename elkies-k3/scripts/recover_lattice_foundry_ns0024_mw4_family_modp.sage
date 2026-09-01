#!/usr/bin/env sage -python
"""Recover the optimally marked NS0024 MW4 source component modulo p.

The surface chart has split I7, I5 and I4 fibres at 0, 1 and infinity. Its
fibre equations are imposed directly through discriminant jets, keeping the
system cubic rather than expanding symbolic square-root branches. The four
sections are the exact minimum-pole, maximum-node-incidence basis certified by
``certify_lattice_foundry_ns0024_mw4_basis.sage``.

This script is finite-field only. It exports a low-degree system for msolve;
no characteristic-zero reconstruction is attempted here.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial


ROOT = Path.cwd().resolve()
DEFAULT_SOURCE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
DEFAULT_BASIS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"


def nonzero_coefficients(polynomial):
    return [entry for entry in polynomial if entry != 0]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--basis", type=Path, default=DEFAULT_BASIS)
parser.add_argument("--slice-a1", type=int)
parser.add_argument("--r1", type=int, help="fix the split I5 node chart")
parser.add_argument("--ri", type=int, help="fix the split I4 node chart")
parser.add_argument("--mw3-seed", type=Path)
parser.add_argument(
    "--mw3-only",
    action="store_true",
    help="recover the resolved three-polynomial-section locus before adjoining P4",
)
parser.add_argument(
    "--pair-gates",
    action="store_true",
    help=(
        "mark the unique smooth intersection of P4 with each polynomial "
        "basis section; this removes the wrong P4 incidence components"
    ),
)
parser.add_argument("--groebner", action="store_true")
parser.add_argument("--export-msolve", type=Path)
args = parser.parse_args()
if (args.r1 is None) != (args.ri is None):
    raise SystemExit("--r1 and --ri must be supplied together")
if args.mw3_only and args.pair_gates:
    raise SystemExit("--pair-gates requires P4")

source = json.loads(args.source.resolve().read_text())
basis = json.loads(args.basis.resolve().read_text())
assert source["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_MW4_SOURCE_AND_NIEMEIER_CERTIFICATE"
assert source["source"]["root_type"] == "A3+A4+A6"
assert basis["status"] == "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS"
assert [item["components_I7_I5_I4"] for item in basis["basis"]] == [
    [1, 0, 0], [2, 1, 3], [2, 1, 1], [1, 1, 1]
]

prime = ZZ(args.prime)
assert prime.is_prime() and prime not in (2, 3, 5, 7)
field = GF(prime)

surface_names = [f"a{index}" for index in range(1, 9)]
if args.r1 is None:
    surface_names += ["r1", "ri", "s1", "si"]
p1_names = [f"p1x{index}" for index in range(1, 5)] + [f"p1y{index}" for index in range(1, 7)]
p2_names = ["p2x2"] + [f"p2y{index}" for index in range(1, 5)]
p3_names = ["p3x2"] + [f"p3y{index}" for index in range(1, 5)]
p4_names = [] if args.mw3_only else (
    ["c", "c0_inverse", "c1_inverse"]
    + [f"p4x{index}" for index in range(1, 5)]
    + [f"p4y{index}" for index in range(1, 8)]
)
pair_gate_names = (
    sum(([f"z{index}", f"z{index}_inverse"] for index in range(1, 4)), [])
    if args.pair_gates else []
)
names = surface_names + p1_names + p2_names + p3_names + p4_names + pair_gate_names
coefficient_ring = PolynomialRing(field, names, order="degrevlex")
v = coefficient_ring.gens_dict()
polynomial_ring = PolynomialRing(coefficient_ring, "t")
t = polynomial_ring.gen()

a = [coefficient_ring(-3)] + [v[f"a{index}"] for index in range(1, 9)]
A = polynomial_ring(a)
if args.r1 is None:
    r1, ri, s1, si = (v[name] for name in ("r1", "ri", "s1", "si"))
else:
    r1, ri = field(args.r1), field(args.ri)
    if not r1 or not ri:
        raise SystemExit("fixed node values must be nonzero")
    s1, si = 1 / r1, 1 / ri
equations = []

# P1 has profile (1,0,0): it meets only the I7 node.
X1 = 1 + sum(v[f"p1x{index}"] * t**index for index in range(1, 5))
Y1 = sum(v[f"p1y{index}"] * t**index for index in range(1, 7))


def all_node_polynomial(prefix):
    """Resolved depth-two I7 section through all three labelled nodes.

    At t=0 the formal double-root branch satisfies
    ``center^2=-A/3`` and starts ``1-a1*t/6``.  Components +/-2 are the
    strict-transform chart ``X=center mod t^2``.  Substituting this before
    elimination removes the singular tangent direction left by mere node
    incidence.
    """
    x1, x2 = -a[1] / 6, v[f"{prefix}x2"]
    x3 = r1 - 1 - ri - x1 - x2
    X = 1 + x1 * t + x2 * t**2 + x3 * t**3 + ri * t**4
    ys = [v[f"{prefix}y{index}"] for index in range(1, 5)]
    y5 = -sum(ys)
    Y = sum(ys[index - 1] * t**index for index in range(1, 5)) + y5 * t**5
    return X, Y


X2, Y2 = all_node_polynomial("p2")
X3, Y3 = all_node_polynomial("p3")

if not args.mw3_only:
    # P4 has P.O=1 and profile (1,1,1). With h=t-c, node values determine
    # X(0), X(1), X(infinity) and Y(0), Y(1), Y(infinity).
    c = v["c"]
    h = t - c
    p4x = [v[f"p4x{index}"] for index in range(1, 5)]
    x5 = r1 * (1 - c)**2 - c**2 - ri - sum(p4x)
    X4 = c**2 + sum(p4x[index - 1] * t**index for index in range(1, 5)) + x5 * t**5 + ri * t**6
    p4y = [v[f"p4y{index}"] for index in range(1, 8)]
    y8 = -sum(p4y)
    Y4 = sum(p4y[index - 1] * t**index for index in range(1, 8)) + y8 * t**8
    equations += [c * v["c0_inverse"] - 1, (c - 1) * v["c1_inverse"] - 1]

# Eliminate all B coefficients through the first exact section identity.
B = Y1**2 - X1**3 - A * X1
assert B.degree() <= 12 and B[0] == 2
b = [B[index] for index in range(13)]
delta = 4 * A**3 + 27 * B**2
equations += [delta[index] for index in range(7)]
equations += [
    sum(delta[index] * binomial(index, jet) for index in range(jet, 25))
    for jet in range(5)
]
equations += [delta[index] for index in range(21, 25)]
equations += [
    sum(a) + 3 * r1**2,
    sum(b) - 2 * r1**3,
    a[8] + 3 * ri**2,
    b[12] - 2 * ri**3,
]
if args.r1 is None:
    equations += [r1 * s1 - 1, ri * si - 1]

equations += nonzero_coefficients(Y2**2 - X2**3 - A * X2 - B)
equations += nonzero_coefficients(Y3**2 - X3**3 - A * X3 - B)
if not args.mw3_only:
    equations += nonzero_coefficients(Y4**2 - X4**3 - A * X4 * h**4 - B * h**6)

# The lattice marking has P_i.P4=1 for i=1,2,3.  The raw section identities
# also contain collision branches at the reducible nodes; requiring a marked
# equality away from 0,1,c is the inexpensive pair gate used successfully in
# the earlier Q80 modular recovery.  Candidate families are still subjected
# to the exact chord-cancellation and group-law intersection verifier.
if args.pair_gates:
    for index, (Xi, Yi) in enumerate(((X1, Y1), (X2, Y2), (X3, Y3)), 1):
        z = v[f"z{index}"]
        equations += [
            (Xi * h**2 - X4)(t=z),
            (Yi * h**3 - Y4)(t=z),
            z * (z - 1) * (z - c) * v[f"z{index}_inverse"] - 1,
        ]
equations = [equation for equation in equations if equation != 0]

if args.mw3_seed is not None:
    if args.r1 is not None:
        raise SystemExit("--mw3-seed cannot be combined with a fixed node chart")
    fields = {}
    for item in args.mw3_seed.resolve().read_text().strip().split("|")[1:]:
        key, value = item.split("=", 1)
        fields[key] = value
    assert int(fields["p"]) == prime
    parsed = {
        key: [field(int(value)) for value in fields[key].split(",")]
        for key in ("A", "B", "P1X", "P1Y", "P2X", "P2Y", "P3X", "P3Y")
    }
    fixes = {v[f"a{index}"]: parsed["A"][index] for index in range(1, 9)}
    fixes.update({r1: field(int(fields["r1"])), ri: field(int(fields["ri"]))})
    fixes.update({s1: 1 / fixes[r1], si: 1 / fixes[ri]})
    for prefix in ("p1", "p2", "p3"):
        upper = prefix.upper()
        x_indices = range(1, 5) if prefix == "p1" else (2,)
        y_count = 6 if prefix == "p1" else 4
        for index in x_indices:
            fixes[v[f"{prefix}x{index}"]] = parsed[f"{upper}X"][index]
        if prefix != "p1":
            assert parsed[f"{upper}X"][1] == -parsed["A"][1] / 6
        for index in range(1, y_count + 1):
            fixes[v[f"{prefix}y{index}"]] = parsed[f"{upper}Y"][index]
    equations += [generator - value for generator, value in fixes.items()]
    assert all((B[index] - parsed["B"][index]).subs(fixes) == 0 for index in range(13))

if args.slice_a1 is not None:
    equations.append(a[1] - field(args.slice_a1))

ideal = coefficient_ring.ideal(equations)
print(
    f"NS0024MW4SYSTEM|p={prime}|variables={coefficient_ring.ngens()}"
    f"|equations={len(equations)}|slice_a1={args.slice_a1}"
    f"|pair_gates={int(args.pair_gates)}|mw3_only={int(args.mw3_only)}"
    f"|r1={r1}|ri={ri}",
    flush=True,
)

if args.export_msolve is not None:
    output = args.export_msolve.resolve()
    output.write_text(
        ",".join(names) + "\n" + str(prime) + "\n"
        + ",\n".join(str(equation) for equation in equations) + "\n"
    )
    print(f"NS0024MW4MSOLVE|output={output}", flush=True)

if args.groebner:
    gb = ideal.groebner_basis(algorithm="libsingular:slimgb")
    print(
        f"NS0024MW4GROEBNER|p={prime}|basis={len(gb)}"
        f"|dimension={ideal.dimension()}|unit={int(coefficient_ring.one() in ideal)}",
        flush=True,
    )
