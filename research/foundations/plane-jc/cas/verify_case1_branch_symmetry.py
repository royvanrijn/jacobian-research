#!/usr/bin/env python3
"""Verify the terminal exact data for the two Case-1 branches.

The raw degree-35 coefficient serialization is not directly equivariant under
the visible sign change.  Decode through the pinned archive's exact quintic
field implementation first, then apply the involution fixing h and negating
(u1,u2), together with the recorded row units.

Before using that quotient as a coefficient field, compare the archive's actual
modulus with the displayed quintic and reduce it modulo 67.  Exact modular
irreducibility then implies irreducibility over Q by Gauss's lemma.  The pinned
regression prime 71 is checked as well.  The archived image of q=u^7 is also
checked against the intrinsic degree-five eliminant.  Since the field degree is
the prime 5 and that image is non-rational, this identifies the quotient field
exactly.  The archived degree-35 eliminant is then checked to be h(u^7), making
its coefficient algebra the monic rank-seven extension L[u]/(u^7-q).  This
finite free extension is faithfully flat, so terminal unit identities replayed
over the old presentation descend to the intrinsic quintic system.

The new adjacent-minor proof, replayed by the following CI step, forces every
hard-ideal solution onto h=N=0.  This checker therefore also replays the pinned
serialized Nullstellensatz certificates excluding that special fibre in both
sign branches.  It intentionally does not replay the superseded 89 MB general
h-membership certificate.
"""
from __future__ import annotations

import pickle
import subprocess
import sys
from pathlib import Path

import sympy as sp

REPO = Path(__file__).resolve().parents[2]
EXACT_REPLAY = REPO / (
    "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)
sys.path.insert(0, str(EXACT_REPLAY))

from degree5_core import (  # noqa: E402
    L,
    MOD,
    decode_l,
    decode_poly,
    sign_substitution,
)

BRANCH1 = EXACT_REPLAY / "hne0_polred.pkl"
BRANCH2 = EXACT_REPLAY / "hne0_branch2_polred.pkl"
FIRSTBLOCK_OUTPUT = EXACT_REPLAY / "firstblock_Q_exact.out"
ROW_SCALES = (1, 1, -1, -1, -1, -1)
SERIALIZED_CHECKER = EXACT_REPLAY / "verify_serialized_certificates.py"
SERIALIZED_MARKERS = (
    "CASE2_SERIALIZED_EXACT_PASS",
    "s=c H0_SERIALIZED_EXACT_PASS",
    "s=-c H0_SERIALIZED_EXACT_PASS",
    "ALL_SERIALIZED_EXACT_CERTIFICATES_PASS",
)
QUOTIENT_ELIMINANT = (
    -1888043347611739526396142670327809715470336,
    586529490054134032292876680565455306752,
    591414847960503971284831143987840,
    265472843532245531128968765,
    62410476400737833472,
    9374377445732,
)
PHI = decode_l(
    {
        0: (-9725570295901, 12623962),
        1: (-1170753213563, 971074),
        2: (-387111042229, 12623962),
        3: (1578225240619, 12623962),
        4: (-469713794365, 6311981),
    }
)

w = sp.Symbol("w")
expected_minpoly = sp.Poly(
    w**5 - w**4 + 3 * w**3 + 3 * w**2 + 26,
    w,
    domain=sp.QQ,
)
archive_minpoly = sp.Poly(
    sum(
        sp.Rational(int(MOD[index].p), int(MOD[index].q)) * w**index
        for index in range(len(MOD))
    ),
    w,
    domain=sp.QQ,
)
assert archive_minpoly == expected_minpoly, (
    "the archive's Case-1 quintic modulus differs from the audited polynomial"
)
for prime in (67, 71):
    reduction = sp.Poly(archive_minpoly.as_expr(), w, modulus=prime)
    assert reduction.is_irreducible, (
        f"the Case-1 quintic is reducible at audited prime {prime}"
    )
print("CASE1_QUINTIC_FIELD_IRREDUCIBLE_PASS")
print("CASE1_PINNED_PRIME_FIELD_PASS")

quotient_residual = L(0)
for coefficient in reversed(QUOTIENT_ELIMINANT):
    quotient_residual = quotient_residual * PHI + coefficient
assert not quotient_residual, "the intrinsic quotient eliminant does not vanish at q=phi"
assert PHI.p.degree() > 0, "the quotient generator unexpectedly lies in Q"
print("CASE1_QUINTIC_DESCENT_PASS")

a7 = sp.Symbol("a7")
firstblock_line = next(
    line
    for line in FIRSTBLOCK_OUTPUT.read_text(encoding="utf-8").splitlines()
    if line.startswith("L[1]=")
)
archive_degree35 = sp.Poly(
    sp.sympify(
        firstblock_line.split("=", 1)[1].replace("^", "**"),
        locals={"a7": a7},
    ),
    a7,
    domain=sp.ZZ,
)
expected_degree35 = sp.Poly(
    sum(
        coefficient * a7 ** (7 * exponent)
        for exponent, coefficient in enumerate(QUOTIENT_ELIMINANT)
    ),
    a7,
    domain=sp.ZZ,
)
assert archive_degree35 == expected_degree35, (
    "the archived degree-35 eliminant is not the rank-seven pullback h(a7^7)"
)
print("CASE1_RANK_SEVEN_PULLBACK_PASS")

branch1 = [decode_poly(item) for item in pickle.loads(BRANCH1.read_bytes())]
branch2 = [decode_poly(item) for item in pickle.loads(BRANCH2.read_bytes())]
assert len(branch1) == len(branch2) == len(ROW_SCALES)

for index, (left, right, scale) in enumerate(
    zip(branch1, branch2, ROW_SCALES), start=1
):
    transported = sign_substitution(left, (1, -1, -1)) * scale
    assert transported.terms == right.terms, (
        f"degree-five Case-1 branch symmetry fails on residual {index}"
    )

print("SYSTEM_SYMMETRY_PASS")
print("CASE1_HARD_RESIDUAL_SYMMETRY_PASS")

completed = subprocess.run(
    [sys.executable, SERIALIZED_CHECKER.name],
    cwd=EXACT_REPLAY,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=False,
)
print(completed.stdout, end="")
assert completed.returncode == 0, (
    f"serialized exact certificate replay exited with {completed.returncode}"
)
lines = set(completed.stdout.splitlines())
missing = [marker for marker in SERIALIZED_MARKERS if marker not in lines]
assert not missing, f"missing serialized certificate markers: {missing}"

print("CASE1_SPECIAL_FIBRE_CERTIFICATES_PASS")
print("CASE1_SPECIAL_FIBRE_FAITHFUL_DESCENT_PASS")
print("CASE1_BRANCH_SYMMETRY_PASS")
