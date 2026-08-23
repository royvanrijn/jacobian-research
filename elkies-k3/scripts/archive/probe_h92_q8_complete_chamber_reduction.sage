#!/usr/bin/env sage -python
"""
Complete the H92 q8 source Weyl/chamber reduction, alternating horizontal
positive-degree (-2)-curves with all old E7/E8 fibre components.

This corrects the earlier q8 "source-nef" normalization, which only reduced
against the 15 finite fibre roots.

Expected H92 behavior:
  * original D has old-fibre degree 18;
  * one irreducible degree-2 bisection has D.C=-1;
  * after reflecting in it, four E7 fibre components become fixed in sequence
        E7_2, E7_4, E7_3, E7_1;
  * the final movable isotropic class has old-fibre degree 16;
  * its generic horizontal support is 8*O + 8*(-P1);
  * no further negative horizontal (-2)-curve is found through degree 24.

Run:
  sage -python ~/Downloads/probe_h92_q8_complete_chamber_reduction.sage

Optional:
  --repo /path/to/jacobian-research
  --max-degree 24
  --max-horizontal-rounds 12
"""

import argparse
import json
from math import gcd
from functools import reduce
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
        home / "Documents" / "jacobian-research",
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
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "elkies-k3" / "data" / "fibrations").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--max-degree", type=int, default=24)
parser.add_argument("--max-horizontal-rounds", type=int, default=12)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
FRAME = ROOT / "elkies-k3" / "data" / "fibrations" / "kumar_e7e8_mw2_frame_3.txt"
AMBIENT = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"
OUTPUT = GEN / "zz-h92-q8-complete-chamber-reduction.json"

ambient = json.loads(AMBIENT.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"

G = load_gram(FRAME)
assert G.nrows() == 17 and G.is_positive_definite()
NS = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -G)

F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_P1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in twice_minuscule]
    + [0] * 8
    + [1, 0],
)
assert O * NS * O == -2
assert minus_P1 * NS * minus_P1 == -2
assert O * NS * F == minus_P1 * NS * F == 1

highest_e7 = (2, 2, 3, 4, 3, 2, 1)
highest_e8 = (2, 3, 4, 6, 5, 4, 3, 2)
zero19 = vector(ZZ, [0] * 19)
affine_e7 = F - sum(
    (c * simple[i] for i, c in enumerate(highest_e7)), zero19
)
affine_e8 = F - sum(
    (c * simple[7+i] for i, c in enumerate(highest_e8)), zero19
)

fiber_components = simple + (affine_e7, affine_e8)
fiber_names = (
    tuple(f"E7_{i}" for i in range(1,8))
    + tuple(f"E8_{i}" for i in range(1,9))
    + ("E7_affine", "E8_affine")
)
assert len(fiber_components) == 17
for C in fiber_components:
    assert C * NS * C == -2
    assert C * NS * F == 0


def pairings_with_fibres(D):
    return tuple(int(D * NS * C) for C in fiber_components)


def reflect(D, C):
    p = ZZ(D * NS * C)
    assert p < 0
    out = D + p * C
    assert out * NS * out == D * NS * D
    return vector(ZZ, out), int(p)


def reduce_fibre_chamber(D):
    """
    Deterministically reflect against negative finite/affine fibre components.
    Positive old-fibre degree makes the affine Weyl alcove reduction terminate.
    """
    D = vector(ZZ, D)
    steps = []
    for _ in range(10000):
        ps = pairings_with_fibres(D)
        negative = [(i, p) for i, p in enumerate(ps) if p < 0]
        if not negative:
            return D, steps
        i, p = min(negative, key=lambda x: (x[0], x[1]))
        before_df = int(D * NS * F)
        D, coeff = reflect(D, fiber_components[i])
        assert int(D * NS * F) == before_df
        steps.append({
            "component": fiber_names[i],
            "pairing": int(p),
            "reflection_coefficient": coeff,
        })
    raise RuntimeError("fibre chamber reduction did not terminate")


def enumerate_negative_horizontal(D, max_degree):
    """
    Exact bounded search for effective positive-degree (-2)-curves C with D.C<0.

    Let r=D.F and write D=(d0,r,d), C=(c0,n,v).
    Since D^2=0 and C^2=-2,

        2*r*n*(D.C) = Q(r*v-n*d) - 2*r^2.

    Thus D.C<0 iff Q(r*v-n*d) < 2*r^2.
    """
    r = ZZ(D * NS * F)
    assert r > 0
    d = vector(ZZ, D[2:])
    found = {}

    for n in range(1, int(max_degree) + 1):
        nd = n * d
        Gnd = G * nd.column()

        A = matrix(ZZ, 18, 18)
        A[:17, :17] = r*r * G
        cross = -r * Gnd
        for i in range(17):
            A[i,17] = cross[i,0]
            A[17,i] = cross[i,0]
        A[17,17] = ZZ(nd * G * nd) + 1
        assert A.is_positive_definite()

        # Negative pairing requires Qshift <= 2r^2-1, so augmented norm <=2r^2.
        raw = pari(A).qfminim(2*r*r)
        reps = matrix(ZZ, raw[2]).transpose() if int(raw[0]) else matrix(ZZ, 0, 18)

        for row in reps.rows():
            k = ZZ(row[17])
            if abs(k) != 1:
                continue
            rr = vector(ZZ, row)
            if k == -1:
                rr = -rr
            if rr[17] != 1:
                continue

            v = vector(ZZ, rr[:17])
            qv = ZZ(v * G * v)
            if (qv - 2) % (2*n):
                continue
            c0 = (qv - 2) // (2*n)
            C = vector(ZZ, [c0, n] + list(v))
            assert C * NS * C == -2
            assert C * NS * F == n

            p = ZZ(D * NS * C)
            if p >= 0:
                continue

            # Put the curve itself in the old-fibre chamber.  This turns
            # vertical variants into one irreducible horizontal representative.
            Cnef, steps = reduce_fibre_chamber(C)
            pnef = ZZ(D * NS * Cnef)
            if pnef >= 0:
                # This candidate was negative only because vertical fibre
                # components were included in the representative.
                continue

            key = tuple(Cnef)
            found[key] = {
                "degree": int(Cnef * NS * F),
                "D_pairing": int(pnef),
                "class": list(map(int, Cnef)),
                "fiber_pairings": list(pairings_with_fibres(Cnef)),
                "curve_fibre_reduction": steps,
            }

    return sorted(
        found.values(),
        key=lambda rec: (
            rec["D_pairing"],
            rec["degree"],
            tuple(rec["class"]),
        ),
    )


D0 = vector(ZZ, ambient["source_q8_lattice_class"])
assert D0 * NS * D0 == 0
assert D0 * NS * F == 18

D = vector(ZZ, D0)
history = []
fixed_horizontal = []
fixed_vertical = []

# First ensure the supplied class is in the fibre alcove.
D, initial_vertical = reduce_fibre_chamber(D)
fixed_vertical.extend(initial_vertical)
assert D == D0, "existing ambient unexpectedly was not fibre-nef"

for round_index in range(int(args.max_horizontal_rounds)):
    horizontal = enumerate_negative_horizontal(D, args.max_degree)
    print(
        "Q8FULLCHAMBER_ROUND|"
        f"round={round_index}|DF={D*NS*F}|D2={D*NS*D}|"
        f"horizontal_negative={len(horizontal)}|"
        f"fiber_pairings={','.join(map(str,pairings_with_fibres(D)))}",
        flush=True,
    )

    if not horizontal:
        history.append({
            "round": round_index,
            "action": "stop",
            "class": list(map(int, D)),
        })
        break

    chosen = horizontal[0]
    C = vector(ZZ, chosen["class"])
    before = vector(ZZ, D)
    D, coeff = reflect(D, C)
    fixed_horizontal.append({
        "round": round_index,
        "curve": chosen,
        "reflection_coefficient": coeff,
    })

    print(
        "Q8FULLCHAMBER_HORIZONTAL|"
        f"round={round_index}|degree={C*NS*F}|pairing={chosen['D_pairing']}|"
        f"DF_before={before*NS*F}|DF_after={D*NS*F}|"
        f"class={','.join(map(str,C))}",
        flush=True,
    )

    # Horizontal reflection may expose vertical fixed components.
    D, vertical_steps = reduce_fibre_chamber(D)
    fixed_vertical.extend(vertical_steps)
    for step in vertical_steps:
        print(
            "Q8FULLCHAMBER_VERTICAL|"
            f"round={round_index}|component={step['component']}|"
            f"pairing={step['pairing']}",
            flush=True,
        )

    history.append({
        "round": round_index,
        "action": "horizontal_then_vertical",
        "class_before": list(map(int, before)),
        "horizontal": chosen,
        "vertical_steps": vertical_steps,
        "class_after": list(map(int, D)),
        "old_fibre_degree_after": int(D * NS * F),
    })
else:
    raise RuntimeError("max horizontal rounds reached")

remaining = enumerate_negative_horizontal(D, args.max_degree)
assert not remaining

# Final generic-fibre support.  The discovered fixed horizontal component has
# old-fibre degree two and the same U/MW coordinates as O+(-P1), so after one
# horizontal reflection the expected restriction is 8O+8(-P1).
horizontal_degree = int(D * NS * F)
generic_support = None
vertical_difference = None
if horizontal_degree == 16:
    generic = 8*O + 8*minus_P1
    vertical_difference = D - generic
    assert vertical_difference[1] == 0
    assert tuple(vertical_difference[-2:]) == (0,0)
    generic_support = "8*O + 8*(-P1)"

pairing_gcd = reduce(gcd, [abs(int(x)) for x in NS * D if x] or [0])

print(
    "Q8FULLCHAMBER_FINAL|"
    f"DF={D*NS*F}|D2={D*NS*D}|primitive_gcd={pairing_gcd}|"
    f"horizontal_reflections={len(fixed_horizontal)}|"
    f"vertical_reflections={len(fixed_vertical)}|"
    f"fiber_pairings={','.join(map(str,pairings_with_fibres(D)))}",
    flush=True,
)

if vertical_difference is not None:
    print(
        "Q8FULLCHAMBER_SUPPORT|"
        f"generic={generic_support}|"
        f"vertical_difference={','.join(map(str,vertical_difference))}",
        flush=True,
    )

status = "PASS_BOUNDED_COMPLETE_Q8_CHAMBER_REDUCTION"
OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-complete-chamber-reduction.v1",
    "status": status,
    "horizontal_search_max_old_fibre_degree": int(args.max_degree),
    "initial_class": list(map(int, D0)),
    "initial_old_fibre_degree": int(D0 * NS * F),
    "fixed_horizontal": fixed_horizontal,
    "fixed_vertical": fixed_vertical,
    "history": history,
    "final_class": list(map(int, D)),
    "final_square": int(D * NS * D),
    "final_old_fibre_degree": int(D * NS * F),
    "final_fiber_pairings": list(pairings_with_fibres(D)),
    "final_primitive_pairing_gcd": int(pairing_gcd),
    "remaining_negative_horizontal_roots": remaining,
    "generic_fibre_support": generic_support,
    "vertical_difference_from_generic_support": (
        list(map(int, vertical_difference))
        if vertical_difference is not None else None
    ),
    "boundary": (
        "Horizontal nefness is certified only for negative (-2)-curves through "
        "the declared old-fibre degree bound. The fibre-component chamber "
        "reduction itself is exact. This artifact supersedes using the old "
        "degree-18 class directly as the movable q8 source divisor."
    ),
}, indent=2, sort_keys=True) + "\n")

print(
    "Q8FULLCHAMBER_SUMMARY|"
    f"initial_DF={D0*NS*F}|final_DF={D*NS*F}|"
    f"horizontal={len(fixed_horizontal)}|vertical={len(fixed_vertical)}|"
    f"remaining_negative={len(remaining)}|status={status}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
