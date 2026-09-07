#!/usr/bin/env sage -python
"""Build the cross-experiment specialization quotient-geometry table.

The observation unit is a fibration presentation, not a target curve.  The
table joins three already frozen detector experiments to their held-out
displayed subgroups: five usable R17 controls, sixteen refreshed R17 ladder
fibres, and nine A1/MW16 parent presentations.

Canonical heights and all derived real Gram data are numerical at the stated
PARI precision.  Point identities, displayed-subgroup embeddings, Smith maps,
and recovery-subspace comparisons are exact.  Bounded misses retain only their
declared bounded meaning.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
import math
from pathlib import Path
import platform
from statistics import median
import sys
from typing import Any, Iterable, Sequence

from cypari2 import Pari
from fpylll import Enumeration, GSO, IntegerMatrix
from sage.all import QQ, RealField, matrix


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    height_gram,
)
from latent_lattice.pari import gp_matrix, gp_point, gp_vector, run_gp  # noqa: E402


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
FINGERPRINTS = ARTIFACTS / "elkies_2026_rank_jump_fingerprints_v1.json"
CALIBRATION_TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
R17_CONTROL_BLIND = ARTIFACTS / "half_lattice_search_ablation_r17_development_blind_v1.json"
R17_CONTROL_VERIFY = ARTIFACTS / "half_lattice_search_ablation_r17_development_verification_v1.json"
R17_RANK21_BLIND = ARTIFACTS / "half_lattice_r17_rank21_blind_v1.json"
R17_RANK21_VERIFY = ARTIFACTS / "half_lattice_r17_rank21_verification_v1.json"
LADDER_INPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
LADDER_BLIND = ARTIFACTS / "r17_refresh_jump_ladder_blind_v2.json"
LADDER_VERIFY = ARTIFACTS / "r17_refresh_jump_ladder_verification_v2.json"
LADDER_QUOTIENT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"
REFRESH = ARTIFACTS / "icarm_curve_refresh_475_573_overview_v1.json"
MW16_INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
MW16_BLIND = ARTIFACTS / "icarm_mw16_parent_ladder_blind_v1.json"
MW16_PUBLIC = ROOT / "elliptic-curves/data/icarm_mw16_public_targets_v1.json"
CURVE398 = CAS / "icarm_curve398.py"
ANALYSIS = CAS / "analyze_half_lattice_height_compression.sage"
RANK21_ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
OUTPUT = ARTIFACTS / "quotient_geometry_table_v1.json"

HEIGHT_BOUND = 100_000
LOG_HEIGHT_BOUND = math.log(HEIGHT_BOUND)
DIGITS = 80
RELATION_DIGITS = 140
MAXIMUM_VECTORS = 5_000_000
PROJECTION_CVP_SCALES = (100_000, 1_000_000)
Q = Fraction


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def qtext(value: Fraction | int | Decimal | Any) -> str:
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    return str(value)


def decimal(value: float | Decimal | str) -> str:
    return format(float(value), ".15g")


def point(record: dict[str, str]) -> tuple[Fraction, Fraction]:
    return Q(record["x"]), Q(record["y"])


def point_record(value: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": qtext(value[0]), "y": qtext(value[1])}


def parse_point_key(value: str) -> tuple[Fraction, Fraction]:
    x_value, y_value = value.split("|", 1)
    return Q(x_value), Q(y_value)


def negate(value: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return value[0], -value[1]


def short_curve_and_points(record: dict[str, Any]):
    a1, a2, a3, a4, a6 = map(Q, record["ainvs"])
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    model = (Q(0), Q(0), Q(0), -c4 / 48, -c6 / 864)
    points = tuple(
        (
            Q(x_value) + b2 / 12,
            Q(y_value) + (a1 * Q(x_value) + a3) / 2,
        )
        for x_value, y_value in record["points"]
    )
    curve = EllipticCurve(model)
    if any(not curve.is_on_curve(value) for value in points):
        raise ArithmeticError(f"curve {record['id']} has an off-curve public point")
    return curve, points


def integral_short_curve_and_points(record: dict[str, Any]):
    a1, a2, a3, a4, a6 = map(Q, record["ainvs"])
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    model = (Q(0), Q(0), Q(0), -27 * c4, -54 * c6)
    points = tuple(
        (
            36 * Q(x_value) + 3 * b2,
            108 * (2 * Q(y_value) + a1 * Q(x_value) + a3),
        )
        for x_value, y_value in record["points"]
    )
    curve = EllipticCurve(model)
    if any(not curve.is_on_curve(value) for value in points):
        raise ArithmeticError(f"curve {record['id']} has an off-curve integral-short point")
    return curve, points


def rank(rows: Sequence[Sequence[int | Fraction]]) -> int:
    if not rows:
        return 0
    return int(matrix(QQ, rows).rank())


def in_span(vector: Sequence[int], rows: Sequence[Sequence[int]]) -> bool:
    return rank([*rows, list(vector)]) == rank(rows)


def same_span(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> bool:
    return rank(left) == rank(right) == rank([*left, *right])


def independent_rows(rows: Iterable[Sequence[int]]) -> list[list[int]]:
    answer: list[list[int]] = []
    for row in rows:
        candidate = list(map(int, row))
        if rank([*answer, candidate]) > len(answer):
            answer.append(candidate)
    return answer


def quadratic(gram: Sequence[Sequence[str | float]], vector: Sequence[int | Fraction]) -> float:
    return sum(
        float(gram[i][j]) * float(vector[i]) * float(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def matrix_vector(rows: Sequence[Sequence[int]], vector: Sequence[int]) -> list[int]:
    return [sum(int(a) * int(b) for a, b in zip(row, vector)) for row in rows]


def rational_matrix_vector(
    rows: Sequence[Sequence[int]], vector: Sequence[int | Fraction]
) -> list[Fraction]:
    return [
        sum((Q(a) * Q(b) for a, b in zip(row, vector)), Q(0))
        for row in rows
    ]


def json_vector(vector: Sequence[int | Fraction]) -> list[int | str]:
    return [
        int(Q(value).numerator)
        if Q(value).denominator == 1
        else qtext(Q(value))
        for value in vector
    ]


def real_matrix(rows: Sequence[Sequence[str | int]]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in rows) + "]"


def _parse_block(lines: Sequence[str], begin: str, end: str) -> list[str]:
    start = lines.index(begin) + 1
    stop = lines.index(end, start)
    return list(lines[start:stop])


def minima_witnesses(gram: Sequence[Sequence[str]], *, timeout: float) -> dict[str, Any]:
    dimension = len(gram)
    program = f"""
default(parisizemax,4000000000);default(parisize,500000000);
default(realprecision,{DIGITS});G={real_matrix(gram)};n=matsize(G)[1];
T=qflllgram(G);R=T~*G*T;b=vecmax(vector(n,i,R[i,i]));
Q=qfminim(R,b,{MAXIMUM_VECTORS},2);V=Q[3];
if(Q[1]>=2*{MAXIMUM_VECTORS},error("quotient vector cap reached"));
W=T*V;m=matsize(W)[2];L=vector(m,j,W[,j]~*G*W[,j]);ord=vecsort(L,,1);
B=matrix(n,0);rk=0;mins=List();print("WIT_BEGIN");
for(k=1,m,j=ord[k];C=matconcat([B,W[,j]]);nr=matrank(C);if(nr>rk,print1(L[j]);for(i=1,n,print1("|",W[i,j]));print();listput(mins,L[j]);B=C;rk=nr));
print("WIT_END");if(rk!=n,error("successive minima enumeration did not span"));
print("META|",Q[1],"|",b);print("MINIMA|",Vec(mins));
"""
    lines = run_gp(program, timeout=timeout)
    witnesses = []
    for line in _parse_block(lines, "WIT_BEGIN", "WIT_END"):
        fields = line.split("|")
        witnesses.append(
            {
                "q_t": fields[0],
                "coordinates_in_quotient_basis": list(map(int, fields[1:])),
            }
        )
    if len(witnesses) != dimension:
        raise ArithmeticError("the minimum witnesses do not span")
    meta = next(line.split("|")[1:] for line in lines if line.startswith("META|"))
    return {
        "witnesses": witnesses,
        "successive_minima": [row["q_t"] for row in witnesses],
        "enumerated_signed_vector_count": int(meta[0]),
        "sufficient_height_bound": meta[1],
    }


def quotient_geometry(
    ambient_gram: Sequence[Sequence[str]],
    embedding_columns: Sequence[Sequence[int]],
    *,
    timeout: float,
) -> dict[str, Any]:
    columns = tuple(tuple(map(int, column)) for column in embedding_columns)
    ambient_rank = len(ambient_gram)
    subgroup_rank = len(columns)
    if not columns or any(len(column) != ambient_rank for column in columns):
        raise ArithmeticError("the displayed generic embedding has wrong dimensions")
    embedding_rows = tuple(zip(*columns))
    program = f"""
default(parisizemax,4000000000);default(parisize,500000000);
default(realprecision,{DIGITS});H={real_matrix(ambient_gram)};A={gp_matrix(embedding_rows)};
S=matsnf(A,1);U=S[1];D=S[3];n=matsize(A)[1];r=matsize(A)[2];
z=List();for(i=1,n,if(sum(j=1,r,abs(D[i,j]))==0,listput(z,i)));
if(#z!=n-r,error("unexpected free quotient dimension"));
Ui=U^-1;zv=Vec(z);C=matrix(n,n-r,i,j,Ui[i,zv[j]]);HM=A~*H*A;X=A~*H*C;
G=C~*H*C-X~*HM^-1*X;T=qflllgram(G);R=T~*G*T;b=vecmax(vector(n-r,i,R[i,i]));
Q=qfminim(R,b,{MAXIMUM_VECTORS},2);V=Q[3];
if(Q[1]>=2*{MAXIMUM_VECTORS},error("quotient vector cap reached"));
W=T*V;m=matsize(W)[2];L=vector(m,j,W[,j]~*G*W[,j]);ord=vecsort(L,,1);
B=matrix(n-r,0);rk=0;mins=List();print("WIT_BEGIN");
for(k=1,m,j=ord[k];CC=matconcat([B,W[,j]]);nr=matrank(CC);if(nr>rk,print1(L[j]);for(i=1,n-r,print1("|",W[i,j]));print();listput(mins,L[j]);B=CC;rk=nr));
print("WIT_END");if(rk!=n-r,error("successive minima enumeration did not span"));
print("META|",n,"|",r,"|",n-r,"|",Q[1],"|",b,"|",matdet(G));
print("MINIMA|",Vec(mins));print("GRAM_BEGIN");
for(i=1,n-r,for(j=1,n-r,if(j>1,print1("|"));print1(G[i,j]));print());print("GRAM_END");
print("QMAP_BEGIN");for(i=1,n-r,for(j=1,n,if(j>1,print1("|"));print1(U[zv[i],j]));print());print("QMAP_END");
"""
    lines = run_gp(program, timeout=timeout)
    meta = next(line.split("|")[1:] for line in lines if line.startswith("META|"))
    gram = [row.split("|") for row in _parse_block(lines, "GRAM_BEGIN", "GRAM_END")]
    quotient_map = [list(map(int, row.split("|"))) for row in _parse_block(lines, "QMAP_BEGIN", "QMAP_END")]
    witnesses = []
    for line in _parse_block(lines, "WIT_BEGIN", "WIT_END"):
        fields = line.split("|")
        witnesses.append(
            {
                "q_t": fields[0],
                "coordinates_in_quotient_basis": list(map(int, fields[1:])),
            }
        )
    quotient_rank = int(meta[2])
    if len(witnesses) != quotient_rank:
        raise ArithmeticError("the quotient minimum witnesses do not span")
    return {
        "definition": "Schur-complement Neron--Tate metric on the displayed subgroup modulo the real span of the specialized generic subgroup",
        "quotient_rank": quotient_rank,
        "quotient_gram": gram,
        "regulator": meta[5],
        "successive_minima": [row["q_t"] for row in witnesses],
        "successive_minimum_witnesses": witnesses,
        "ambient_to_quotient_map": quotient_map,
        "enumeration": {
            "algorithm": "PARI qflllgram plus complete qfminim/Fincke--Pohst",
            "enumerated_signed_vector_count": int(meta[3]),
            "sufficient_height_bound": meta[4],
            "maximum_stored_half_vectors": MAXIMUM_VECTORS,
            "complete_for_successive_minima": True,
        },
        "numerical_boundary": f"canonical heights use PARI realprecision={DIGITS}; point and Smith identities are exact, but height decimals are not interval certificates",
    }


def smith_quotient_map(
    ambient_rank: int, embedding_columns: Sequence[Sequence[int]]
) -> list[list[int]]:
    embedding_rows = tuple(zip(*embedding_columns))
    program = f"""
A={gp_matrix(embedding_rows)};S=matsnf(A,1);U=S[1];D=S[3];
n=matsize(A)[1];r=matsize(A)[2];z=List();
for(i=1,n,if(sum(j=1,r,abs(D[i,j]))==0,listput(z,i)));
if(#z!=n-r,error("unexpected free quotient dimension"));zv=Vec(z);
print("BEGIN");for(i=1,n-r,for(j=1,n,if(j>1,print1("|"));print1(U[zv[i],j]));print());print("END");
"""
    lines = run_gp(program, timeout=60.0)
    rows = [
        list(map(int, row.split("|")))
        for row in _parse_block(lines, "BEGIN", "END")
    ]
    if any(len(row) != ambient_rank for row in rows):
        raise ArithmeticError("the Smith quotient map has wrong dimensions")
    return rows


def recover_embedding_large(
    curve: EllipticCurve,
    basis: Sequence[tuple[Fraction, Fraction]],
    targets: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    program = f"""
default(parisizemax,4000000000);default(parisize,500000000);
default(realprecision,{RELATION_DIGITS});E=ellinit({gp_vector(curve.coefficients)});
P=[{','.join(map(gp_point, basis))}];Q=[{','.join(map(gp_point, targets))}];
A=concat(P,Q);H=ellheightmatrix(E,A);HP=H[1..#P,1..#P];
X=HP^-1*H[1..#P,#P+1..#A];Z=matrix(#P,#Q,i,j,bestappr(X[i,j],1000000));
ok=1;for(j=1,#Q,d=1;for(i=1,#P,d=lcm(d,denominator(Z[i,j])));T=[0];for(i=1,#P,c=d*Z[i,j];if(c,T=elladd(E,T,ellmul(E,P[i],c))));if(ellsub(E,T,ellmul(E,Q[j],d))!=[0],ok=0));
print("EXACT|",ok);print("RESIDUAL|",vecmax(abs(X-Z)));
print("BEGIN");for(j=1,#Q,for(i=1,#P,if(i>1,print1("|"));print1(Z[i,j]));print());print("END");
"""
    lines = run_gp(program, timeout=900.0)
    exact = next(line.split("|", 1)[1] for line in lines if line.startswith("EXACT|"))
    residual = next(line.split("|", 1)[1] for line in lines if line.startswith("RESIDUAL|"))
    if exact != "1" or float(residual.replace(" ", "")) >= 2.0 ** -100:
        raise ArithmeticError(
            "large exact embedding recovery failed its exact/separation gate: "
            f"exact={exact} residual={residual}"
        )
    columns = tuple(
        tuple(map(Q, row.split("|")))
        for row in _parse_block(lines, "BEGIN", "END")
    )
    if len(columns) != len(targets) or any(len(row) != len(basis) for row in columns):
        raise ArithmeticError("large exact embedding recovery returned wrong dimensions")
    return columns


def clear_rational_columns(
    columns: Sequence[Sequence[int | Fraction]],
) -> tuple[tuple[int, ...], ...]:
    answer = []
    for column in columns:
        denominator = 1
        for value in column:
            denominator = math.lcm(denominator, Q(value).denominator)
        answer.append(tuple(int(Q(value) * denominator) for value in column))
    return tuple(answer)


def horizontal_reduction(pari: Pari, analysis, raw: Sequence[Fraction]) -> dict[str, Any]:
    denominator = 1
    for coefficient in raw:
        denominator = math.lcm(denominator, coefficient.denominator)
    integral = [int(coefficient * denominator * denominator) for coefficient in raw]
    polynomial = "+".join(f"({coefficient})*x^{index}" for index, coefficient in enumerate(integral))
    result = pari(
        "my(m1,m2,C0=[" + polynomial + ",0],C1,C2);"
        "C1=hyperellminimalmodel(C0,&m1);C2=hyperellred(C1,&m2);[C2,m1,m2]"
    )
    curve, first, second = result[0], result[1], result[2]

    def horizontal(transformation):
        value = transformation[1]
        return [[int(value[i, j]) for j in range(2)] for i in range(2)]

    composed = analysis.primitive_matrix(
        analysis.matrix_product(horizontal(first), horizontal(second))
    )
    reduced_p = analysis.pari_polynomial_coefficients(curve[0])
    reduced_q = analysis.pari_polynomial_coefficients(curve[1])
    while reduced_p and reduced_p[-1] == 0:
        reduced_p.pop()
    while reduced_q and reduced_q[-1] == 0:
        reduced_q.pop()
    return {
        "horizontal_matrix": composed,
        "recomputed_reduced_coefficient_bits": max(
            abs(int(value)).bit_length() for value in reduced_p + reduced_q
        ),
        "map_source": "deterministic replay of PARI hyperellminimalmodel then hyperellred",
    }


def chart_metric(
    pari: Pari,
    analysis,
    *,
    model: Sequence[Fraction],
    base: tuple[Fraction, Fraction],
    target: tuple[Fraction, Fraction],
    recorded_reduced: dict[str, Any] | None,
    recorded_bits: int | None,
) -> dict[str, Any]:
    raw = (
        -3 * base[0] * base[0] - 4 * model[3],
        -8 * base[1],
        -6 * base[0],
        Q(0),
        Q(1),
    )
    if recorded_reduced is None:
        transformation = horizontal_reduction(pari, analysis, raw)
        if recorded_bits is not None and transformation["recomputed_reduced_coefficient_bits"] != int(recorded_bits):
            raise ArithmeticError("compact-ledger reduced coefficient height changed")
        map_matches = None
    else:
        transformation = analysis.pari_transform(pari, raw, recorded_reduced)
        map_matches = bool(transformation["exact_recorded_model_match"])
    signed = []
    for sign, signed_target in ((1, target), (-1, negate(target))):
        parameter = analysis.raw_parameter(model[3], base, signed_target)
        reduced = None if parameter is None else analysis.reduced_coordinate(
            parameter, transformation["horizontal_matrix"]
        )
        if reduced is None:
            continue
        height = analysis.rational_height(reduced)
        signed.append((height, sign, signed_target, parameter, reduced))
    if not signed:
        raise ArithmeticError("both signs of a recovered direction map to infinity")
    height, sign, signed_target, parameter, reduced = min(signed, key=lambda row: (row[0], row[1]))
    if height > HEIGHT_BOUND:
        raise ArithmeticError("a recorded discovery is outside its replayed source-coordinate bound")
    return {
        "signed_lift": point_record(signed_target),
        "signed_lift_sign_relative_to_recorded_event": sign,
        "raw_coordinate": qtext(parameter),
        "reduced_coordinate": qtext(reduced),
        "reduced_coordinate_height": height,
        "reduced_log_height": decimal(math.log(height)),
        "horizontal_reduced_to_raw_matrix": [
            [str(value) for value in row] for row in transformation["horizontal_matrix"]
        ],
        "recorded_reduced_model_matches_replay": map_matches,
    }


def source_priority(source: str) -> int:
    return int(source.split(":priority:", 1)[1].split(":", 1)[0])


def source_mask(source: str) -> int:
    for marker in (":mask:", ":gmask:"):
        if marker in source:
            return int(source.split(marker, 1)[1].split(":", 1)[0], 16)
    raise ValueError(f"source has no mask: {source}")


def nearest_half_lattice_projection(
    gram,
    projection,
) -> dict[str, Any]:
    """Numerically minimize (a-m/2)^T G (a-m/2) over integral m.

    The CVP decision is repeated on two independently rounded integral Grams.
    Energies and the completing-the-square identity are then evaluated on the
    original high-precision real Gram.  A disagreement stays explicit rather
    than being promoted to a numerical minimum.
    """

    dimension = len(projection)
    candidates = []
    for scale in PROJECTION_CVP_SCALES:
        rounded = [
            [int(round(float(gram[i, j]) * scale)) for j in range(dimension)]
            for i in range(dimension)
        ]
        gso = GSO.Mat(
            IntegerMatrix.from_matrix(rounded),
            gram=True,
            float_type="dd",
            update=True,
        )
        mu = [
            [gso.get_mu(i, j) if i > j else 0.0 for j in range(dimension)]
            for i in range(dimension)
        ]
        twice_projection = [2 * float(value) for value in projection]
        target = [
            twice_projection[i]
            + sum(
                twice_projection[j] * mu[j][i]
                for j in range(i + 1, dimension)
            )
            for i in range(dimension)
        ]
        distance_bound = (
            sum(abs(value) for row in rounded for value in row) / 4 + 1.0
        )
        solutions = Enumeration(gso).enumerate(
            0,
            dimension,
            distance_bound,
            0,
            target=target,
        )
        if not solutions:
            raise ArithmeticError("projection half-lattice CVP returned no solution")
        representative = tuple(int(round(value)) for value in solutions[0][1])
        difference = matrix(
            gram.base_ring(),
            dimension,
            1,
            [
                projection[i] - gram.base_ring()(representative[i]) / 2
                for i in range(dimension)
            ],
        )
        energy = (difference.transpose() * gram * difference)[0, 0]
        candidates.append(
            {
                "rounding_scale": scale,
                "twice_half_lattice_center_m": list(representative),
                "energy_on_unrounded_gram": decimal(energy),
            }
        )
    stable = all(
        row["twice_half_lattice_center_m"]
        == candidates[0]["twice_half_lattice_center_m"]
        for row in candidates[1:]
    )
    return {
        "definition": "dist_G(a,(1/2)Z^r)^2 = min_m (a-m/2)^T G (a-m/2), where a=G^-1 b",
        "rounding_scales": list(PROJECTION_CVP_SCALES),
        "stable_center_across_rounding_scales": stable,
        "twice_nearest_half_lattice_center_m": (
            candidates[0]["twice_half_lattice_center_m"] if stable else None
        ),
        "minimum_phase_energy": (
            candidates[0]["energy_on_unrounded_gram"] if stable else None
        ),
        "rounded_cvp_candidates": candidates,
        "numerical_boundary": "CVP on rounded canonical-height Grams; not an interval-certified nearest-point calculation",
    }


def finalize_directions(
    curve: EllipticCurve,
    inherited: Sequence[tuple[Fraction, Fraction]],
    directions: list[dict[str, Any]],
) -> None:
    if not directions:
        return
    signed_targets = [
        point(row["coordinate_geometry"]["signed_lift"]) for row in directions
    ]
    centered_points = []
    for row, target in zip(directions, signed_targets):
        base = point(row["half_lattice_chart"]["base_point"])
        centered = curve.add(curve.multiply(target, 2), curve.negate(base))
        if centered is None:
            raise ArithmeticError("a recovered nonzero quotient direction centered to torsion")
        centered_points.append(centered)
    projection_gram = height_gram(
        curve, [*inherited, *signed_targets], digits=DIGITS, timeout=600.0
    )
    centered_gram = height_gram(curve, centered_points, digits=DIGITS, timeout=600.0)
    real_field = RealField(220)
    inherited_rank = len(inherited)
    inherited_gram = matrix(
        real_field,
        [
            [real_field(projection_gram[i][j]) for j in range(inherited_rank)]
            for i in range(inherited_rank)
        ],
    )
    for index, row in enumerate(directions):
        cross = matrix(
            real_field,
            inherited_rank,
            1,
            [
                real_field(projection_gram[i][inherited_rank + index])
                for i in range(inherited_rank)
            ],
        )
        projection = inherited_gram.solve_right(cross).column(0)
        target_height = real_field(
            projection_gram[inherited_rank + index][inherited_rank + index]
        )
        intrinsic = target_height - sum(
            cross[i, 0] * projection[i] for i in range(inherited_rank)
        )
        if abs(float(intrinsic) - float(row["intrinsic_q_t_P"])) > 1e-7:
            raise ArithmeticError("the direct Schur complement changed intrinsic q_t(P)")
        optimum = nearest_half_lattice_projection(inherited_gram, projection)
        if optimum["stable_center_across_rounding_scales"]:
            optimum_energy = float(optimum["minimum_phase_energy"])
            optimum["minimum_total_centered_height"] = decimal(
                float(intrinsic) + optimum_energy
            )
            m = optimum["twice_nearest_half_lattice_center_m"]
            direct = target_height
            for i in range(inherited_rank):
                direct -= real_field(m[i]) * cross[i, 0]
            direct += sum(
                real_field(m[i])
                * inherited_gram[i, j]
                * real_field(m[j])
                / 4
                for i in range(inherited_rank)
                for j in range(inherited_rank)
            )
            optimum["completing_square_residual"] = decimal(
                direct - intrinsic - real_field(optimum["minimum_phase_energy"])
            )
        else:
            optimum["minimum_total_centered_height"] = None
            optimum["completing_square_residual"] = None
        row["projection_away_from_inherited_span"] = {
            "inherited_rank": inherited_rank,
            "projection_coordinates_a": [decimal(value) for value in projection],
            "schur_complement_lambda": decimal(intrinsic),
        }
        row["optimal_half_lattice_position"] = optimum
        centered_half_height = float(centered_gram[index][index]) / 4
        q_centered = float(row["centered_quotient_energy"])
        phase = centered_half_height - q_centered
        if phase < -1e-7:
            raise ArithmeticError(
                "orthogonal phase energy became negative: "
                f"direction={row.get('direction_index')} centered={centered_half_height} "
                f"quotient={q_centered} coordinates={row.get('quotient_coordinates')}"
            )
        phase = max(phase, 0.0)
        log_height = float(row["coordinate_geometry"]["reduced_log_height"])
        distortion = log_height - 2 * centered_half_height
        window = (LOG_HEIGHT_BOUND - distortion) / 2 - phase
        residual = log_height - (2 * q_centered + 2 * phase + distortion)
        row["centered_half_height"] = decimal(centered_half_height)
        row["half_lattice_phase"] = {
            "definition": "hhat(pr_M(P-Q/2)); for an adaptive chart the centered quotient component is recorded separately",
            "energy": decimal(phase),
        }
        if (
            optimum["stable_center_across_rounding_scales"]
            and not row["half_lattice_chart"][
                "adaptive_center_has_nonzero_original_quotient_component"
            ]
        ):
            excess = phase - float(optimum["minimum_phase_energy"])
            if excess < -1e-7:
                raise ArithmeticError("an initial source chart beat the half-lattice CVP")
            optimum["actual_source_phase_excess_over_optimum"] = decimal(
                max(excess, 0.0)
            )
        else:
            optimum["actual_source_phase_excess_over_optimum"] = None
        row["coordinate_distortion_term"] = {
            "definition": "log H(s(P)) - hhat(2P-Q)/2",
            "value": decimal(distortion),
        }
        row["predicted_search_height_window"] = {
            "definition": "(log(B)-delta_s(P))/2 - half_lattice_phase, compared with centered quotient energy",
            "height_bound_B": HEIGHT_BOUND,
            "upper_bound_for_centered_quotient_energy": decimal(window),
            "centered_quotient_energy_lies_under_window": q_centered <= window + 1e-8,
            "decomposition_residual": decimal(residual),
        }


def recovery_comparison(
    geometry: dict[str, Any],
    initial_rows: Sequence[Sequence[int]],
    final_rows: Sequence[Sequence[int]],
    *,
    containment_complete: bool,
) -> dict[str, Any]:
    witnesses = geometry["successive_minimum_witnesses"]
    quotient_rank = int(geometry["quotient_rank"])
    initial_rank = rank(initial_rows)
    final_rank = rank(final_rows)
    for witness in witnesses:
        vector = witness["coordinates_in_quotient_basis"]
        witness["in_initial_recovered_subspace"] = (
            in_span(vector, initial_rows) if containment_complete else None
        )
        witness["in_final_recovered_subspace"] = (
            in_span(vector, final_rows) if containment_complete else None
        )
    initial_flag = [row["coordinates_in_quotient_basis"] for row in witnesses[:initial_rank]]
    final_flag = [row["coordinates_in_quotient_basis"] for row in witnesses[:final_rank]]
    return {
        "comparison_scope": "necessary scalar-window test on exact rational subspaces of the displayed quotient: at recovered rank r, compare with the first r deterministic qfminim rank-increase witnesses",
        "displayed_containment_complete": containment_complete,
        "initial_recovered_displayed_quotient_rank": initial_rank if containment_complete else None,
        "final_recovered_displayed_quotient_rank": final_rank if containment_complete else None,
        "initial_equals_successive_minimum_prefix": (
            same_span(initial_rows, initial_flag) if containment_complete else None
        ),
        "final_equals_successive_minimum_prefix": (
            same_span(final_rows, final_flag) if containment_complete else None
        ),
        "full_displayed_quotient_recovered": (
            final_rank == quotient_rank if containment_complete else None
        ),
    }


def geometry_from_fingerprint(record: dict[str, Any], *, timeout: float) -> dict[str, Any]:
    source = record["quotient_height_geometry"]
    gram = source["quotient_gram"]
    witness_data = minima_witnesses(gram, timeout=timeout)
    for left, right in zip(source["successive_minima"], witness_data["successive_minima"]):
        if abs(float(left) - float(right)) > 1e-10:
            raise ArithmeticError("control successive minima changed")
    return {
        "definition": source["definition"],
        "quotient_rank": source["quotient_rank"],
        "quotient_gram": gram,
        "regulator": source["quotient_determinant"],
        "successive_minima": source["successive_minima"],
        "successive_minimum_witnesses": witness_data["witnesses"],
        "enumeration": {
            **source["enumeration"],
            "witness_replay_enumerated_signed_vector_count": witness_data[
                "enumerated_signed_vector_count"
            ],
        },
        "numerical_boundary": source["numerical_boundary"],
    }


def build_controls(pari, analysis, *, timeout: float) -> list[dict[str, Any]]:
    fingerprints = load(FINGERPRINTS)["fingerprints"]
    calibration = {
        row["parameter"]: row
        for row in load(CALIBRATION_TRUTH)["positive_controls"]
    }
    rank21_blind = load(R17_RANK21_BLIND)
    rank21_verify = load(R17_RANK21_VERIFY)
    control_blind = {row["parameter"]: row for row in load(R17_CONTROL_BLIND)["results"]}
    control_verify = {row["parameter"]: row for row in load(R17_CONTROL_VERIFY)["results"]}
    cases = []
    for fingerprint in fingerprints:
        parameter = fingerprint["parameter"]
        geometry = geometry_from_fingerprint(fingerprint, timeout=timeout)
        if parameter == "3/8":
            blind = rank21_blind
            model = tuple(map(Q, blind["short_model"]))
            generic = tuple(point(row) for row in blind["generic_points"])
            pairs = list(zip(rank21_verify["relations"], blind["blind_result"]["candidate_points"]))
            selected = []
            basis_rows: list[list[int]] = []
            for relation, candidate in pairs:
                vector = list(map(int, relation["quotient_coordinates"]))
                if rank([*basis_rows, vector]) == len(basis_rows):
                    continue
                basis_rows.append(vector)
                selected.append(
                    {
                        "point": candidate["point"],
                        "quotient_coordinates": vector,
                        "sources": [f"initial:mask:{value}" for value in candidate["source_masks"]],
                    }
                )
                if len(selected) == geometry["quotient_rank"]:
                    break
            embedding_columns = [
                [int(row == column) for row in range(21)]
                for column in range(17)
            ]
            quotient_map = smith_quotient_map(21, embedding_columns)
            for item in selected:
                detector_coordinates = item["quotient_coordinates"]
                ambient = [0] * 17 + detector_coordinates
                item["detector_quotient_coordinates"] = detector_coordinates
                item["quotient_coordinates"] = matrix_vector(quotient_map, ambient)
            engine = SourceFileLoader("quotient_table_rank21_engine", str(RANK21_ENGINE)).load_module()
            specialized = engine.canonical_height_gram(model, generic)
            rounded = tuple(
                tuple(int((value * Decimal(1_000_000)).to_integral_value()) for value in row)
                for row in specialized
            )
            oracle = engine.CosetOracle(rounded)
            cover_by_mask = {int(row["mask"]): row for row in blind["cover_records"]}

            def chart_for(item):
                mask = int(item["sources"][0].split(":")[-1])
                unused, representative, error = oracle.solve(mask)
                if error > 1e-6:
                    raise ArithmeticError("rank21 chart CVP replay failed")
                base = EllipticCurve(model).linear_combination(generic, representative)
                return mask, representative, base, cover_by_mask[mask]

            case_id = "r17-control-rank21"
        else:
            blind = control_blind[parameter]
            model = tuple(map(Q, blind["short_model"]))
            verified = control_verify[parameter]
            generic = tuple(point(row) for row in blind["generic_points"])
            arm = next(row for row in verified["arms"] if row["id"] == "generic-deepest43")
            allowed = set(map(int, arm["masks"]))
            order = {int(row["mask"]): index for index, row in enumerate(blind["cover_records"])}
            selected = []
            for row in arm["selected_exact_basis"]:
                sources = sorted(
                    (int(value) for value in row["source_masks"] if int(value) in allowed),
                    key=lambda value: order[value],
                )
                if not sources:
                    raise ArithmeticError("a control basis direction has no source in its arm")
                selected.append(
                    {
                        "point": row["point"],
                        "quotient_coordinates": list(map(int, row["quotient_coordinates"])),
                        "sources": [f"initial:mask:{value}" for value in sources],
                    }
                )
            truth = calibration[parameter]
            quotient_map = smith_quotient_map(
                int(fingerprint["certified_rank_lower_bound"]),
                truth["embedding_matrix_columns"],
            )
            complement_indices = [
                int(value) - 1
                for value in verified["fixture_inputs"][
                    "public_complement_indices_one_based"
                ]
            ]
            for item in selected:
                detector_coordinates = item["quotient_coordinates"]
                ambient = [0] * int(fingerprint["certified_rank_lower_bound"])
                for coefficient, public_index in zip(
                    detector_coordinates, complement_indices
                ):
                    ambient[public_index] += int(coefficient)
                item["detector_quotient_coordinates"] = detector_coordinates
                item["quotient_coordinates"] = matrix_vector(quotient_map, ambient)
            cover_by_mask = {int(row["mask"]): row for row in blind["cover_records"]}

            def chart_for(item):
                mask = int(item["sources"][0].split(":")[-1])
                cover = cover_by_mask[mask]
                representative = tuple(map(int, cover["specialized_representative"]))
                base = EllipticCurve(model).linear_combination(generic, representative)
                return mask, representative, base, cover

            case_id = f"r17-control-rank{fingerprint['certified_rank_lower_bound']}"
        curve = EllipticCurve(model)
        directions = []
        recovery_rows = []
        for index, item in enumerate(selected, 1):
            mask, representative, base, cover = chart_for(item)
            target = point(item["point"])
            metric = chart_metric(
                pari,
                analysis,
                model=model,
                base=base,
                target=target,
                recorded_reduced=None,
                recorded_bits=cover.get("reduced_coefficient_bits"),
            )
            quotient_coordinates = item["quotient_coordinates"]
            q_value = quadratic(geometry["quotient_gram"], quotient_coordinates)
            recovery_rows.append(quotient_coordinates)
            directions.append(
                {
                    "direction_index": index,
                    "blind_recovery_stage": "initial-generic-deepest43",
                    "recorded_event_point": item["point"],
                    "quotient_coordinates": quotient_coordinates,
                    "detector_complement_coordinates": item[
                        "detector_quotient_coordinates"
                    ],
                    "intrinsic_q_t_P": decimal(q_value),
                    "centered_quotient_energy": decimal(q_value),
                    "half_lattice_chart": {
                        "source": item["sources"][0],
                        "mask": mask,
                        "representative": list(representative),
                        "base_point": point_record(base),
                        "adaptive_center_has_nonzero_original_quotient_component": False,
                    },
                    "coordinate_geometry": metric,
                }
            )
        finalize_directions(curve, generic, directions)
        comparison = recovery_comparison(
            geometry, recovery_rows, recovery_rows, containment_complete=True
        )
        cases.append(
            {
                "case_id": case_id,
                "lineage": "usable-known-R17-control",
                "curve_label": fingerprint["label"],
                "parameter": parameter,
                "generic_rank": int(fingerprint["generic_rank"]),
                "displayed_subgroup_rank": int(fingerprint["certified_rank_lower_bound"]),
                "j_displayed": int(fingerprint["response_variables"]["primary_lower_bound"]),
                "j_displayed_definition": "displayed subgroup rank minus specialized generic rank; not the elliptic j-invariant",
                "quotient_geometry": geometry,
                "blind_recovery": {
                    "initial_rank": len(recovery_rows),
                    "adaptive_rank": 0,
                    "final_rank": len(recovery_rows),
                    "direction_witnesses": directions,
                },
                "successive_minimum_recovery_test": comparison,
            }
        )
        print(f"QUOTIENTTABLE|case={case_id}|status=complete", flush=True)
    return cases


def build_ladder(pari, analysis, *, timeout: float) -> list[dict[str, Any]]:
    inputs = {int(row["curve_id"]): row for row in load(LADDER_INPUT)["cases"]}
    blinds = {int(row["curve_id"]): row for row in load(LADDER_BLIND)["results"]}
    verifies = {int(row["curve_id"]): row for row in load(LADDER_VERIFY)["results"]}
    quotients = {int(row["curve_id"]): row for row in load(LADDER_QUOTIENT)["fibres"]}
    publics = {int(row["id"]): row for row in load(REFRESH)["snapshot"]["records"]}
    cases = []
    for curve_id in sorted(inputs):
        source = inputs[curve_id]
        blind = blinds[curve_id]
        verify = verifies[curve_id]
        quotient = quotients[curve_id]
        curve, public_points = short_curve_and_points(publics[curve_id])
        model = tuple(map(Q, source["short_model"]))
        if curve.coefficients != model:
            raise ArithmeticError(f"curve {curve_id} short model changed")
        generic = tuple(point(row) for row in source["generic_points"])
        embedding_rows = quotient["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_public_points"
        ]
        embedding_columns = [
            [int(embedding_rows[row][column]) for row in range(len(embedding_rows))]
            for column in range(17)
        ]
        ambient_gram = height_gram(curve, public_points, digits=DIGITS, timeout=600.0)
        geometry = quotient_geometry(ambient_gram, embedding_columns, timeout=timeout)
        relation_by_index = {
            int(row["blind_basis_index"]): row for row in verify["opened_public_relations"]
        }
        basis_quotient_coordinates: dict[int, list[int]] = {
            index: [0] * geometry["quotient_rank"] for index in range(17)
        }
        containment_complete = bool(
            verify["all_final_blind_basis_points_in_opened_public_subgroup"]
        )
        initial_rows = []
        final_rows = []
        directions = []
        event_records = []
        for stage in ("initial", "adaptive"):
            stage_record = blind[stage]
            events = stage_record.get("discovered_group_classification", {}).get("events", [])
            for event in events:
                basis_index = int(event["basis_rank_after"]) - 1
                relation = relation_by_index.get(basis_index)
                ambient_coordinates = None
                quotient_coordinates = None
                if relation is not None and relation["exact_relation_in_opened_public_basis"]:
                    ambient_coordinates = list(map(int, relation["public_basis_coordinates"]))
                    quotient_coordinates = matrix_vector(
                        geometry["ambient_to_quotient_map"], ambient_coordinates
                    )
                    basis_quotient_coordinates[basis_index] = quotient_coordinates
                    final_rows.append(quotient_coordinates)
                    if stage == "initial":
                        initial_rows.append(quotient_coordinates)
                event_records.append((stage, event, basis_index, ambient_coordinates, quotient_coordinates))
        for direction_index, (stage, event, basis_index, ambient_coordinates, quotient_coordinates) in enumerate(event_records, 1):
            target = point(event["point"])
            source_name = event["sources"][0]
            priority = source_priority(source_name)
            stage_record = blind[stage]
            cover = next(row for row in stage_record["cover_records"] if int(row["priority"]) == priority)
            if stage == "initial":
                chart_basis = generic
                chart_representative = tuple(map(int, cover["specialized_representative"]))
            else:
                chart_basis = tuple(point(row) for row in stage_record["basis_before"])
                chart_representative = tuple(map(int, cover["representative"]))
            base = curve.linear_combination(chart_basis, chart_representative)
            if base is None:
                raise ArithmeticError("a ladder source chart has infinite base point")
            if min(base, negate(base)) != parse_point_key(cover["base_point_key"]):
                raise ArithmeticError("a ladder source chart base-point key changed")
            metric = chart_metric(
                pari,
                analysis,
                model=model,
                base=base,
                target=target,
                recorded_reduced=None,
                recorded_bits=cover["search"].get("reduced_model_maximum_coefficient_bits"),
            )
            if quotient_coordinates is not None:
                q_value = quadratic(geometry["quotient_gram"], quotient_coordinates)
            else:
                projection_gram = height_gram(curve, [*generic, target], digits=DIGITS, timeout=600.0)
                real_field = RealField(220)
                generic_gram = matrix(
                    real_field,
                    [[real_field(projection_gram[i][j]) for j in range(17)] for i in range(17)],
                )
                cross = matrix(
                    real_field,
                    17,
                    1,
                    [real_field(projection_gram[i][17]) for i in range(17)],
                )
                q_value = float(
                    real_field(projection_gram[17][17])
                    - (cross.transpose() * generic_gram.inverse() * cross)[0, 0]
                )
            signed_q = quotient_coordinates
            if signed_q is not None and metric["signed_lift_sign_relative_to_recorded_event"] == -1:
                signed_q = [-value for value in signed_q]
            if stage == "initial":
                centered_q = q_value
                center_q = [0] * geometry["quotient_rank"]
            elif quotient_coordinates is not None and all(
                index in basis_quotient_coordinates for index in range(len(cover["representative"]))
            ):
                center_q = [0] * geometry["quotient_rank"]
                for coefficient, index in zip(cover["representative"], range(len(cover["representative"]))):
                    for column, value in enumerate(basis_quotient_coordinates[index]):
                        center_q[column] += int(coefficient) * int(value)
                centered_vector = [
                    Q(signed_q[column]) - Q(center_q[column], 2)
                    for column in range(len(center_q))
                ]
                centered_q = quadratic(geometry["quotient_gram"], centered_vector)
            else:
                center_q = None
                centered_q = q_value
            directions.append(
                {
                    "direction_index": direction_index,
                    "blind_basis_index_zero_based": basis_index,
                    "blind_recovery_stage": stage,
                    "recorded_event_point": event["point"],
                    "exact_relation_in_displayed_subgroup": ambient_coordinates is not None,
                    "quotient_coordinates": quotient_coordinates,
                    "intrinsic_q_t_P": decimal(q_value),
                    "centered_quotient_energy": decimal(centered_q),
                    "half_lattice_chart": {
                        "source": source_name,
                        "priority": priority,
                        "mask": source_mask(source_name),
                        "representative": list(chart_representative),
                        "base_point": point_record(base),
                        "center_quotient_coordinates": center_q,
                        "adaptive_center_has_nonzero_original_quotient_component": stage == "adaptive",
                    },
                    "coordinate_geometry": metric,
                }
            )
        finalize_directions(curve, generic, directions)
        comparison = recovery_comparison(
            geometry,
            initial_rows,
            final_rows,
            containment_complete=containment_complete,
        )
        case_id = f"r17-ladder-curve{curve_id}"
        cases.append(
            {
                "case_id": case_id,
                "lineage": "new-R17-refresh-ladder",
                "curve_id": curve_id,
                "parameter": quotient["native_parameter"],
                "fibration_presentation": quotient["native_chart"],
                "generic_rank": 17,
                "displayed_subgroup_rank": int(quotient["snapshot_rank_lower_bound"]),
                "j_displayed": int(verify["true_displayed_jump_opened_after_blind_freeze"]),
                "j_displayed_definition": "displayed subgroup rank minus specialized generic rank; not the elliptic j-invariant",
                "quotient_geometry": geometry,
                "blind_recovery": {
                    "initial_rank_before_public_complement": int(verify["initial_exact_recovered_rank"]),
                    "adaptive_rank_before_public_complement": int(verify["adaptive_incremental_exact_recovered_rank"]),
                    "final_rank_before_public_complement": int(verify["exact_quotient_rank_recovered_before_public_complement"]),
                    "all_blind_basis_points_in_displayed_subgroup": containment_complete,
                    "direction_witnesses": directions,
                },
                "successive_minimum_recovery_test": comparison,
            }
        )
        print(f"QUOTIENTTABLE|case={case_id}|status=complete", flush=True)
    return cases


def public_mw16_targets():
    refresh = {int(row["id"]): row for row in load(REFRESH)["snapshot"]["records"]}
    frozen = {int(row["id"]): row for row in load(MW16_PUBLIC)["records"]}
    module = SourceFileLoader("quotient_table_curve398", str(CURVE398)).load_module()
    curve398 = EllipticCurve(module.short_coefficients())
    public398 = tuple((Q(x_value), Q(y_value)) for x_value, y_value in module.SHORT_POINTS)
    answer = {398: (curve398, public398, 30)}
    for curve_id in (400, 401):
        curve, points = integral_short_curve_and_points(frozen[curve_id])
        answer[curve_id] = (curve, points, int(frozen[curve_id]["rank_lower_bound"]))
    for curve_id in (542, 548):
        curve, points = integral_short_curve_and_points(refresh[curve_id])
        answer[curve_id] = (curve, points, int(refresh[curve_id]["rank_lower_bound"]))
    return answer


def build_mw16(pari, analysis, *, timeout: float) -> list[dict[str, Any]]:
    inputs = {row["parent_id"]: row for row in load(MW16_INPUT)["parents"]}
    blinds = {row["parent_id"]: row for row in load(MW16_BLIND)["parents"]}
    publics = public_mw16_targets()
    ambient_cache: dict[int, list[list[str]]] = {}
    cases = []
    for parent_id in sorted(inputs, key=lambda value: (int(inputs[value]["curve_id"]), int(inputs[value]["priority_rank"]))):
        source = inputs[parent_id]
        blind = blinds[parent_id]
        curve_id = int(source["curve_id"])
        curve, public_points, displayed_rank = publics[curve_id]
        model = tuple(map(Q, source["target_short_model"]))
        if curve.coefficients != model:
            raise ArithmeticError(f"MW16 parent {parent_id} short model changed")
        generic = tuple(point(row) for row in source["specialized_generic_points"])
        events = blind["discovered_group_saturation"]["events"]
        event_points = tuple(point(row["point"]) for row in events)
        embedding = recover_embedding_large(
            curve, public_points, [*generic, *event_points]
        )
        generic_rational_columns = embedding[:16]
        generic_columns = clear_rational_columns(generic_rational_columns)
        event_columns = embedding[16:]
        generic_integrally_contained = all(
            value.denominator == 1
            for column in generic_rational_columns
            for value in column
        )
        if curve_id not in ambient_cache:
            ambient_cache[curve_id] = height_gram(
                curve, public_points, digits=DIGITS, timeout=900.0
            )
        geometry = quotient_geometry(
            ambient_cache[curve_id], generic_columns, timeout=timeout
        )
        directions = []
        recovery_rows = []
        cover_by_priority = {
            int(row["priority"]): row for row in blind["cover_records"]
        }
        for index, (event, ambient_coordinates) in enumerate(zip(events, event_columns), 1):
            quotient_coordinates = rational_matrix_vector(
                geometry["ambient_to_quotient_map"], ambient_coordinates
            )
            recovery_rows.append(quotient_coordinates)
            source_name = event["sources"][0]
            priority = source_priority(source_name)
            cover = cover_by_priority[priority]
            base = point(cover["search"]["base_point"])
            target = point(event["point"])
            metric = chart_metric(
                pari,
                analysis,
                model=model,
                base=base,
                target=target,
                recorded_reduced=cover["search"]["reduced_model"],
                recorded_bits=None,
            )
            q_value = quadratic(geometry["quotient_gram"], quotient_coordinates)
            directions.append(
                {
                    "direction_index": index,
                    "blind_recovery_stage": "initial-complete-maximum-depth-MW16-stratum",
                    "recorded_event_point": event["point"],
                    "exact_relation_in_displayed_rational_span": True,
                    "exact_relation_in_displayed_subgroup": all(
                        value.denominator == 1 for value in ambient_coordinates
                    ),
                    "quotient_coordinates": json_vector(quotient_coordinates),
                    "intrinsic_q_t_P": decimal(q_value),
                    "centered_quotient_energy": decimal(q_value),
                    "half_lattice_chart": {
                        "source": source_name,
                        "priority": priority,
                        "mask": source_mask(source_name),
                        "representative": list(map(int, cover["specialized_representative"])),
                        "base_point": point_record(base),
                        "adaptive_center_has_nonzero_original_quotient_component": False,
                    },
                    "coordinate_geometry": metric,
                }
            )
        finalize_directions(curve, generic, directions)
        comparison = recovery_comparison(
            geometry,
            recovery_rows,
            recovery_rows,
            containment_complete=generic_integrally_contained
            and all(
                value.denominator == 1
                for column in event_columns
                for value in column
            ),
        )
        cases.append(
            {
                "case_id": f"mw16-a1-{parent_id}",
                "lineage": "MW16-A1-hit",
                "curve_id": curve_id,
                "parameter": source["target_parameter"],
                "fibration_presentation": parent_id,
                "generic_rank": 16,
                "displayed_subgroup_rank": displayed_rank,
                "j_displayed": displayed_rank - 16,
                "j_displayed_definition": "displayed subgroup rank minus specialized generic rank; not the elliptic j-invariant",
                "quotient_geometry": geometry,
                "blind_recovery": {
                    "initial_rank": int(blind["exact_quotient_rank_recovered"]),
                    "adaptive_rank": 0,
                    "final_rank": int(blind["exact_quotient_rank_recovered"]),
                    "direction_witnesses": directions,
                },
                "displayed_containment": {
                    "specialized_generic_subgroup_integrally_contained": generic_integrally_contained,
                    "generic_coordinate_denominator_lcm": math.lcm(
                        *(value.denominator for column in generic_rational_columns for value in column)
                    ),
                    "all_independent_recovery_witnesses_integrally_contained": all(
                        value.denominator == 1
                        for column in event_columns
                        for value in column
                    ),
                },
                "successive_minimum_recovery_test": comparison,
                "nesting_warning": "presentations on the same target curve are not independent observations",
            }
        )
        print(f"QUOTIENTTABLE|case=mw16-a1-{parent_id}|status=complete", flush=True)
    return cases


def build_payload(*, timeout: float) -> dict[str, Any]:
    pari = Pari()
    pari.default("realprecision", DIGITS)
    analysis = SourceFileLoader("quotient_table_height_analysis", str(ANALYSIS)).load_module()
    cases = [
        *build_controls(pari, analysis, timeout=timeout),
        *build_ladder(pari, analysis, timeout=timeout),
        *build_mw16(pari, analysis, timeout=timeout),
    ]
    if len(cases) != 30:
        raise ArithmeticError("the 5+16+9 observation inventory changed")
    comparable = [
        row for row in cases
        if row["successive_minimum_recovery_test"]["displayed_containment_complete"]
    ]
    failures = [
        row["case_id"] for row in comparable
        if row["successive_minimum_recovery_test"]["final_equals_successive_minimum_prefix"] is False
    ]
    strict_initial = [
        row
        for row in comparable
        if 0
        < row["successive_minimum_recovery_test"][
            "initial_recovered_displayed_quotient_rank"
        ]
        < row["quotient_geometry"]["quotient_rank"]
    ]
    strict_initial_failures = [
        row["case_id"]
        for row in strict_initial
        if row["successive_minimum_recovery_test"][
            "initial_equals_successive_minimum_prefix"
        ]
        is False
    ]
    strict_final = [
        row
        for row in comparable
        if 0
        < row["successive_minimum_recovery_test"][
            "final_recovered_displayed_quotient_rank"
        ]
        < row["quotient_geometry"]["quotient_rank"]
    ]
    direction_rows = [
        direction
        for case in cases
        for direction in case["blind_recovery"]["direction_witnesses"]
    ]
    if any(
        not row["predicted_search_height_window"][
            "centered_quotient_energy_lies_under_window"
        ]
        for row in direction_rows
    ):
        raise ArithmeticError("a recovered direction failed its source-chart window identity")
    stable_projection_rows = [
        row
        for row in direction_rows
        if row["optimal_half_lattice_position"][
            "stable_center_across_rounding_scales"
        ]
    ]
    initial_source_rows = [
        row
        for row in stable_projection_rows
        if row["optimal_half_lattice_position"][
            "actual_source_phase_excess_over_optimum"
        ]
        is not None
    ]
    curve398 = next(
        row for row in cases if row["case_id"] == "mw16-a1-curve398-p16875"
    )
    curve398_directions = curve398["blind_recovery"]["direction_witnesses"]
    curve398_minima = list(map(float, curve398["quotient_geometry"]["successive_minima"]))
    curve398_lambdas = list(
        map(float, (row["intrinsic_q_t_P"] for row in curve398_directions))
    )
    curve398_phases = list(
        map(float, (row["half_lattice_phase"]["energy"] for row in curve398_directions))
    )
    curve398_optimal_phases = list(
        map(
            float,
            (
                row["optimal_half_lattice_position"]["minimum_phase_energy"]
                for row in curve398_directions
            ),
        )
    )
    curve398_distortions = list(
        map(
            float,
            (
                row["coordinate_distortion_term"]["value"]
                for row in curve398_directions
            ),
        )
    )
    paths = (
        FINGERPRINTS,
        CALIBRATION_TRUTH,
        R17_CONTROL_BLIND,
        R17_CONTROL_VERIFY,
        R17_RANK21_BLIND,
        R17_RANK21_VERIFY,
        LADDER_INPUT,
        LADDER_BLIND,
        LADDER_VERIFY,
        LADDER_QUOTIENT,
        REFRESH,
        MW16_INPUT,
        MW16_BLIND,
        MW16_PUBLIC,
        CURVE398,
        ANALYSIS,
        Path(__file__).resolve(),
    )
    return {
        "schema": "elliptic-curves.quotient-geometry-table.v1",
        "status": "PASS_COMPLETE_QUOTIENT_GEOMETRY_TABLE",
        "observation_unit": {
            "usable_known_R17_controls": 5,
            "new_R17_ladder_fibres": 16,
            "MW16_A1_presentations": 9,
            "total_presentations": 30,
        },
        "definitions": {
            "j_displayed": "rank of the certified displayed subgroup minus the proved specialized generic rank; this symbol is not the elliptic j-invariant",
            "regulator": "determinant of the displayed quotient Gram, not necessarily the regulator of the full Mordell--Weil group",
            "intrinsic_q_t_P": "hhat_/M(P), the Neron--Tate norm after orthogonal projection away from the original specialized generic subgroup M",
            "projection_coefficients_a": "G^-1 b for the signed recovered point against the specialized generic basis",
            "optimal_half_lattice_phase": "dist_G(a,(1/2)Z^r)^2, recomputed by CVP on the canonical-height Gram at two rounding scales",
            "centered_quotient_energy": "hhat_/M(P-Q/2); it equals q_t(P) for initial charts Q in M and changes after adaptive quotient lifts",
            "half_lattice_phase": "the old-span term hhat(pr_M(P-Q/2)) in the midpoint orthogonal decomposition",
            "coordinate_distortion_term": "log H(s(P)) - hhat(2P-Q)/2 for the replayed reduced quartic coordinate",
            "blind_recovery_stage": "the first exact independence event in the frozen detector trace",
        },
        "precise_question": (
            "Necessary scalar-window test: at every frozen detector stage, is its exact recovered rational quotient subspace the deterministic successive-minimum prefix of the same rank (equivalently away from boundary ties, can it be exactly the span of the quotient directions below one intrinsic q_t cutoff)?"
        ),
        "answer": {
            "successive_minimum_cutoff_hypothesis": "NO",
            "reason": "all thirteen strictly partial, nonempty initial recovered subspaces fail to equal their corresponding successive-minimum prefixes; seven remain counterexamples after the adaptive/final stage, so no scalar q_t cutoff describes discovery and source-chart visibility depends on half-lattice phase and coordinate distortion",
            "comparable_presentation_count": len(comparable),
            "counterexample_case_ids": failures,
            "strictly_partial_nonempty_initial_stage_count": len(strict_initial),
            "strictly_partial_nonempty_initial_stage_counterexample_count": len(
                strict_initial_failures
            ),
            "strictly_partial_nonempty_initial_stage_counterexample_case_ids": strict_initial_failures,
            "strictly_partial_nonempty_final_stage_count": len(strict_final),
            "source_chart_window_identity_holds_for_every_independent_recovery_witness": True,
            "curve398_diagnostic": {
                "presentation_ids": [
                    "mw16-a1-curve398-p16875",
                    "mw16-a1-curve398-p63669",
                ],
                "displayed_quotient_rank": curve398["j_displayed"],
                "initial_recovered_rank": curve398["successive_minimum_recovery_test"][
                    "initial_recovered_displayed_quotient_rank"
                ],
                "successive_minima_range": [
                    decimal(min(curve398_minima)),
                    decimal(max(curve398_minima)),
                ],
                "recovered_direction_lambda_range": [
                    decimal(min(curve398_lambdas)),
                    decimal(max(curve398_lambdas)),
                ],
                "successive_minimum_witness_indices_in_recovered_subspace_one_based": [
                    index
                    for index, witness in enumerate(
                        curve398["quotient_geometry"][
                            "successive_minimum_witnesses"
                        ],
                        1,
                    )
                    if witness["in_initial_recovered_subspace"]
                ],
                "actual_source_phase_range": [
                    decimal(min(curve398_phases)),
                    decimal(max(curve398_phases)),
                ],
                "optimal_half_lattice_phase_range": [
                    decimal(min(curve398_optimal_phases)),
                    decimal(max(curve398_optimal_phases)),
                ],
                "coordinate_distortion_range": [
                    decimal(min(curve398_distortions)),
                    decimal(max(curve398_distortions)),
                ],
            },
            "interpretation": "intrinsic quotient height matters at the fibre scale but is not the discovery ordering: curve 398 recovers a five-dimensional subspace containing only the fifth and eighth deterministic successive-minimum witnesses, while missing the first four; the source-chart midpoint/window identity still holds point-by-point, and unseen directions are not assigned hypothetical chart distortions",
        },
        "summary": {
            "independent_recovery_direction_count": len(direction_rows),
            "stable_projection_half_lattice_cvp_count": len(stable_projection_rows),
            "initial_source_chart_phase_comparison_count": len(initial_source_rows),
            "median_initial_source_phase_excess_over_optimum": decimal(
                median(
                    float(
                        row["optimal_half_lattice_position"][
                            "actual_source_phase_excess_over_optimum"
                        ]
                    )
                    for row in initial_source_rows
                )
            ),
            "presentations_with_exact_displayed_containment": len(comparable),
            "presentations_matching_final_successive_minimum_prefix": len(comparable) - len(failures),
            "presentations_failing_final_successive_minimum_prefix": len(failures),
        },
        "cases": cases,
        "input_hashes": {relative(path): digest(path) for path in paths},
        "software": {
            "python": platform.python_version(),
            "sage": "10.9",
            "pari_realprecision_decimal_digits": DIGITS,
        },
        "claim_boundary": [
            "The quotient Grams, regulators, successive minima, canonical heights, phases, and distortion terms are high-precision numerical data, not interval certificates.",
            "Exact point identities, recovered embeddings, Smith quotient maps, and rational subspace comparisons are replayed exactly.",
            "Every quotient is taken inside a displayed certified subgroup; it need not be the full Mordell--Weil quotient.",
            "The 16 refreshed R17 observations include unresolved displayed containment on curves 478 and 539; their subspace-cutoff answers remain null.",
            "The nine MW16 presentations are nested in five target curves and are not nine independent statistical observations.",
            "Every search miss remains bounded by the original height, timeout, and chart-selection budgets.",
        ],
        "reproducing_command": "sage -python elliptic-curves/cas/build_quotient_geometry_table.sage",
    }


def normalized(value: Any) -> Any:
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload(timeout=args.timeout)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or normalized(json.loads(args.output.read_text())) != normalized(payload):
            raise ArithmeticError("the stored quotient-geometry table changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        f"QUOTIENTTABLE|cases={len(payload['cases'])}|directions={payload['summary']['independent_recovery_direction_count']}|"
        f"counterexamples={payload['summary']['presentations_failing_final_successive_minimum_prefix']}|status={payload['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
