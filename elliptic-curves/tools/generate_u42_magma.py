#!/usr/bin/env python3
"""Emit exact Magma V2.29 code for the Nagao u=42 rank probes.

The stored points live on the rational short Jacobian, not on the integral
minimal model.  The emitted program therefore verifies them on that curve and
uses the forward isomorphism returned by ``MinimalModel`` before asking Magma
for rank or Selmer information.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[2]
SOURCE = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_u42_height_10000000.json"
)
SOURCE_SHA256 = "4fea0207fd637988bcc1147143657cbec5c2404cb81b4c4a487e2dde20cc43b8"

A_NUM = -74879150695093957092648257365083
A_DEN = 92236816
B_NUM = 121839825430716337244033564674334552301153773691
B_DEN = 442921190432

EXPECTED_MINIMAL = (
    1,
    1,
    1,
    -713297808023681194679243421,
    7164472716086503658294244200686835329779,
)


def magma_rational(value: str | int) -> str:
    q = Fraction(value)
    return f"(Q!({q.numerator})/Q!({q.denominator}))"


def load_points() -> list[tuple[str, str]]:
    raw = SOURCE.read_bytes()
    actual_sha256 = hashlib.sha256(raw).hexdigest()
    if actual_sha256 != SOURCE_SHA256:
        raise SystemExit(
            f"refusing changed input: expected {SOURCE_SHA256}, got {actual_sha256}"
        )
    data = json.loads(raw)
    basis = data["small_prime_saturation"]["saturated_basis"]
    points = [(point["jacobian_x"], point["jacobian_y"]) for point in basis]
    if len(points) != 17:
        raise SystemExit(f"expected 17 points, got {len(points)}")
    return points


def setup_code(points: list[tuple[str, str]]) -> str:
    point_rows = ",\n".join(
        f"    Eshort![ {magma_rational(x)}, {magma_rational(y)}, Q!1 ]"
        for x, y in points
    )
    expected = ", ".join(str(value) for value in EXPECTED_MINIMAL)
    return f"""// Nagao u=42, T=3631/14.  Exact input; no GRH class-group setting.
// Source SHA256: {SOURCE_SHA256}
Q := Rationals();
A := {magma_rational(Fraction(A_NUM, A_DEN))};
B := {magma_rational(Fraction(B_NUM, B_DEN))};
Eshort := EllipticCurve([ Q | 0, 0, 0, A, B ]);

Pshort := [
{point_rows}
];

// Point construction and these identities are exact over Q.
assert #Pshort eq 17;
assert &and[ P[3] eq 1 : P in Pshort ];
assert &and[ P[2]^2 eq P[1]^3 + A*P[1] + B : P in Pshort ];
print "EXACT_SHORT_POINTS_VERIFIED", #Pshort;

// The JSON points are not coordinates on the known integral minimal model.
// MinimalModel returns Eshort -> Emin as its second return value.
Emin, short_to_min, min_to_short := MinimalModel(Eshort);
assert IsMinimalModel(Emin);
Pmin := [ short_to_min(P) : P in Pshort ];
assert #Pmin eq 17;

Eknown := EllipticCurve([ Q | {expected} ]);
assert IsMinimalModel(Eknown);
assert IsIsomorphic(Emin, Eknown);
print "MAGMA_MINIMAL_A_INVARIANTS", aInvariants(Emin);
print "KNOWN_MINIMAL_MODEL_ISOMORPHIC", true;

independent := IsLinearlyIndependent(Pmin);
print "POINTS_INDEPENDENT_MOD_TORSION", independent;
assert independent;
print "RIGOROUS_RANK_LOWER_BOUND", #Pmin;
"""


def probe_code(mode: str) -> str:
    if mode == "verify":
        return 'print "PROBE", "verification only";\n'
    if mode == "rankbounds":
        return """procedure RunRankBounds()
lower, upper := RankBounds(Emin : Effort := 1);
assert upper ge #Pmin;
print "MAGMA_RANK_BOUNDS", lower, upper;
print "COMBINED_RANK_BOUNDS", #Pmin, upper;
end procedure;
RunRankBounds();
"""
    if mode == "twoselmer":
        return """procedure RunTwoSelmer()
// Bound=-1 is the documented unconditional class-group computation.
SetVerbose("TwoDescent", 2);
T2, T2map := TwoTorsionSubgroup(Emin);
S2, AtoS := TwoSelmerGroup(Emin : Bound := -1);
selmer_dimension := Ngens(S2);
torsion_dimension := Ngens(T2);
upper := selmer_dimension - torsion_dimension;
assert upper ge #Pmin;
print "TWO_TORSION_DIMENSION", torsion_dimension;
print "TWO_SELMER_INVARIANTS", Invariants(S2);
print "TWO_SELMER_DIMENSION", selmer_dimension;
print "RIGOROUS_RANK_BOUNDS_FROM_2_SELMER", #Pmin, upper;
end procedure;
RunTwoSelmer();
"""
    if mode == "twodescent":
        return """procedure RunTwoDescent()
// Factoring out the 17 known points avoids constructing 2^17 covers.
SetVerbose("TwoDescent", 1);
known := { P : P in Pmin };
covers := TwoDescent(Emin :
    RemoveTorsion := true,
    RemoveGens := known,
    WithMaps := false);
nclasses := #covers + 1; // TwoDescent omits the trivial covering.
residual_dimension := 0;
while nclasses gt 1 do
    assert nclasses mod 2 eq 0;
    nclasses div:= 2;
    residual_dimension +:= 1;
end while;
print "NONTRIVIAL_RESIDUAL_TWO_COVERS", #covers;
print "RESIDUAL_TWO_SELMER_DIMENSION", residual_dimension;
print "RIGOROUS_RANK_BOUNDS_FROM_QUOTIENT", #Pmin,
    #Pmin + residual_dimension;
end procedure;
RunTwoDescent();
"""
    raise AssertionError(mode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("verify", "twoselmer", "twodescent", "rankbounds"),
        default="twoselmer",
    )
    args = parser.parse_args()
    print(setup_code(load_points()), end="")
    print(probe_code(args.mode), end="")


if __name__ == "__main__":
    main()
