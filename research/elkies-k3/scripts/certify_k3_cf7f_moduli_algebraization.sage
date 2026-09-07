#!/usr/bin/env sage-python
"""Certify the exact Q-chart and conditional Shimura genus for det-714.

The marked ``I5+I7+I7`` MW1 equations already define an affine scheme over
``QQ``.  This certificate records that algebraization, checks its smooth
one-dimensional ``GF(7)`` branch, and independently evaluates Ogg's genus and
fixed-point formulas for the candidate Shimura datum ``(D,N)=(51,7)``.

The comparison between the marked coefficient component and the full
Atkin--Lehner quotient is deliberately left conditional: the integral even
Clifford order has not yet been proved locally conjugate to the Eichler order,
and no birational map from the coefficient chart to a genus-two equation is
known.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    BinaryQF_reduced_representatives,
    QQ,
    ZZ,
    kronecker,
    matrix,
    pari,
    prod,
)


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
ADAPTER = GEN / (
    "elkies-k3-k3-cf7f6c91a3a40d32-"
    "source-search-target-partner2-lattice-only-v1.json"
)
FORMAL = GEN / (
    "elkies-k3-k3-cf7f6c91a3a40d32-"
    "a4-2a6-mw1-formal-smoothness-v1.json"
)
HENSEL = GEN / (
    "elkies-k3-k3-cf7f6c91a3a40d32-"
    "a4-2a6-mw1-marked-gf7-hensel-v1.json"
)
DEFAULT_OUTPUT = GEN / (
    "elkies-k3-k3-cf7f6c91a3a40d32-"
    "a4-2a6-mw1-moduli-algebraization-v1.json"
)


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def factors(value):
    return list(ZZ(value).factor())


def phi(value):
    return ZZ(value) * prod(QQ(prime - 1) / prime for prime, _ in factors(value))


def psi(value):
    return prod((prime + 1) * prime ** (exponent - 1) for prime, exponent in factors(value))


def base_genus(discriminant, level):
    e2 = (
        0
        if level % 4 == 0
        else prod(1 - kronecker(-4, prime) for prime, _ in factors(discriminant))
        * prod(1 + kronecker(-4, prime) for prime, _ in factors(level))
    )
    e3 = (
        0
        if level % 9 == 0
        else prod(1 - kronecker(-3, prime) for prime, _ in factors(discriminant))
        * prod(1 + kronecker(-3, prime) for prime, _ in factors(level))
    )
    genus = 1 + phi(discriminant) * psi(level) / 12 - QQ(e2) / 4 - QQ(e3) / 3
    if genus.denominator() != 1:
        raise ArithmeticError("Shimura genus formula is not integral")
    return int(genus), int(e2), int(e3)


def squarefree_part(value):
    return prod(prime for prime, exponent in factors(value) if exponent % 2)


def order_class_number(discriminant):
    return len(BinaryQF_reduced_representatives(int(discriminant), primitive_only=True))


def local_embedding_count(discriminant, level, prime, field_discriminant, conductor):
    symbol = 1 if conductor % prime == 0 else kronecker(field_discriminant, prime)
    if discriminant % prime == 0:
        return 1 - symbol
    exponent = ZZ(level).valuation(prime)
    conductor_exponent = ZZ(conductor).valuation(prime)
    if exponent == 1:
        return 1 + symbol
    if exponent >= 2 + 2 * conductor_exponent:
        return 2 * (prime ** conductor_exponent + prime ** max(conductor_exponent - 1, 0)) if kronecker(field_discriminant, prime) == 1 else 0
    raise NotImplementedError("only the squarefree (51,7) case is certified here")


def cm_fixed_count(discriminant, level, involution, field_discriminant, conductor):
    complement = ZZ(discriminant * level // involution)
    local = prod(
        local_embedding_count(discriminant, level, prime, field_discriminant, conductor)
        for prime, _ in factors(complement)
    )
    return order_class_number(field_discriminant * conductor**2) * local


def involution_fixed_count(discriminant, level, involution):
    if involution == 2:
        return cm_fixed_count(discriminant, level, involution, -4, 1) + cm_fixed_count(
            discriminant, level, involution, -8, 1
        )
    squarefree = ZZ(squarefree_part(involution))
    conductor = ZZ(involution // squarefree).isqrt()
    if involution % 4 == 3:
        return cm_fixed_count(
            discriminant, level, involution, -squarefree, conductor
        ) + cm_fixed_count(discriminant, level, involution, -squarefree, 2 * conductor)
    return cm_fixed_count(
        discriminant, level, involution, -4 * squarefree, conductor
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

adapter = json.loads(ADAPTER.read_text())
formal = json.loads(FORMAL.read_text())
hensel = json.loads(HENSEL.read_text())
if adapter["surface_id"] != "K3-cf7f6c91a3a40d32" or adapter["determinant"] != 714:
    raise ValueError("wrong determinant-714 adapter")
if formal["status"] != "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY":
    raise ValueError("missing formal source branch")
if hensel["jacobian_certificate"]["rank_mod_prime"] != 39:
    raise ValueError("marked seed no longer has codimension 39")

t_gate = adapter["t_arithmetic_pre_solver_gate"]
t_gram = matrix(ZZ, t_gate["literal_transcendental_gram"])
if t_gram.det() != -714 or tuple(map(int, pari(t_gram).qfsign())) != (2, 1):
    raise ArithmeticError("unexpected transcendental lattice")
clifford = t_gate["clifford"]
if clifford["quaternion_discriminant"] != 51:
    raise ArithmeticError("unexpected quaternion discriminant")
if clifford["integral_even_clifford_order"]["local_level_index"] != 7:
    raise ArithmeticError("unexpected local level index")

D, N = 51, 7
top_genus, elliptic_2, elliptic_3 = base_genus(D, N)
full_group = [1, 3, 7, 17, 21, 51, 119, 357]
fixed_counts = {
    str(involution): int(involution_fixed_count(D, N, involution))
    for involution in full_group[1:]
}
fixed_total = sum(fixed_counts.values())
quotient_genus = QQ(2 * top_genus - 2 - fixed_total) / (2 * len(full_group)) + 1
if top_genus != 21 or fixed_total != 24 or quotient_genus != 2:
    raise ArithmeticError("(51,7) full Atkin--Lehner genus regression failed")

payload = {
    "schema": "elkies-k3.k3-cf7f-a4-2a6-mw1-moduli-algebraization.v1",
    "status": "PASS_EXACT_Q_AFFINE_SOURCE_CHART_AND_CONDITIONAL_GENUS_TWO_SHIMURA_QUOTIENT",
    "surface_id": "K3-cf7f6c91a3a40d32",
    "determinant": 714,
    "exact_affine_Q_chart": {
        "ambient": "A^40_Q with coordinates a0..a8,b0..b12,c0,n0..n6,m0..m9",
        "normalization": "a0=-27",
        "weierstrass": "y^2=x^3+A(t)x+B(t)",
        "marked_section": "x=N(t)/(t+c0)^2, y=M(t)/(t+c0)^3",
        "fibre_orders": {"t=0": 5, "t=1": 7, "t=infinity": 7},
        "component_depths": [1, 2, 1],
        "displayed_equations": 47,
        "independent_equations_on_certified_localization": 39,
        "formally_smooth_relative_dimension_at_GF7_seed": 1,
        "definition": (
            "Use the exact coefficient equations constructed by "
            "certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage; the formal "
            "certificate proves that its selected localized component is a curve."
        ),
    },
    "candidate_shimura_arithmetic": {
        "transcendental_gram": [list(map(int, row)) for row in t_gram.rows()],
        "quaternion_discriminant": D,
        "eichler_level_candidate": N,
        "base_curve": f"X_0^{D}({N})",
        "base_genus": top_genus,
        "elliptic_orbits_order_2": elliptic_2,
        "elliptic_orbits_order_3": elliptic_3,
        "full_atkin_lehner_group": full_group,
        "fixed_points_by_nonidentity_involution": fixed_counts,
        "fixed_point_total": fixed_total,
        "full_quotient_genus": int(quotient_genus),
        "formula": "g(Y)=1+(2g(X)-2-sum_{w!=1}#Fix(w))/(2|W|)",
    },
    "literature": {
        "reference": "Padurariu--Saia, Shimura curve Atkin--Lehner quotients of genus at most two, arXiv:2509.25368",
        "independent_table_entry": "(D,N,W,g)=(51,7,<w3,w7,w17>,2)",
        "code_snapshot": "https://github.com/fsaia/GenusAtMost2/commit/6cc368fe37aa67187783118f18d149b2b1fd6230",
    },
    "proof_boundary": {
        "proved": (
            "The normalized marked source is an explicit affine Q-scheme and the "
            "selected GF(7) point lies on a one-dimensional formally smooth component. "
            "Ogg's formulas independently give genus 21 for X_0^51(7) and genus two "
            "for its full Atkin--Lehner quotient."
        ),
        "not_proved": (
            "The even Clifford order is not yet certified to be Eichler, and no "
            "birational map identifies the marked coefficient component with the "
            "genus-two quotient. No hyperelliptic equation, Q-rational moduli point, "
            "or rational parameterization is claimed."
        ),
    },
    "inputs": {
        relative(ADAPTER): digest(ADAPTER),
        relative(FORMAL): digest(FORMAL),
        relative(HENSEL): digest(HENSEL),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_cf7f_moduli_algebraization.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = arguments.output.resolve()
if arguments.check:
    if not output.exists() or output.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output}")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
print(
    "K3CF7FMODULI|base_genus=21|full_AL_fixed=24|quotient_genus=2|"
    "marked_chart_dimension=1|status=PASS"
)
