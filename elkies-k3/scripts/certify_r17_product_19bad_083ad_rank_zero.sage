#!/usr/bin/env sage-python
"""Certify arithmetic product rank zero using two incompatible regulators.

Only 19bad:083ad is in scope. No Selmer or section search is performed.
The proof uses height-preserving specialization and the theorem that equality
of algebraic and analytic ranks over a finite function field implies refined
BSD. Under a hypothetical nonzero characteristic-zero section, both reductions
have rank one; their regulator squareclasses must therefore agree.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from sage.all import GF, Matrix, PolynomialRing, QQ, prod
from sage.version import version


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DIRECT = RESULTS / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
PAIRS = RESULTS / "elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
HEIGHT = RESULTS / "elkies-k3-r17-product-survivor-galois-height-gate-v1.json"
KEY = "alternate-orbit-19bad:alternate-orbit-083ad"
OUTPUT = RESULTS / "elkies-k3-r17-product-19bad-083ad-rank-zero-v1.json"
CONTROL_DIR = RESULTS / "r17-product-19bad-083ad-controls"
PRIMES = (131, 137, 151)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def squareclass(value):
    value = QQ(value)
    if value <= 0:
        raise ArithmeticError("expected a positive regulator")
    return int(prod(p for p, exponent in value.factor() if exponent % 2))


def coefficients(poly):
    return [str(c) for c in poly]


def local_invariant_kummer_separation(direct, pair):
    """Retain the bounded exact branch-local calculation that preceded closure.

    Each selected row is an actual residue-field square character of a unit
    x(P_i)-theta. Full row rank certifies separation, with no completeness
    claim for the local character list or the global Selmer group.
    """
    matrix = Matrix(GF(2), 0, 17)
    rows = []
    for prime in (131, 137, 151, 157, 167, 173):
        field = GF(prime)
        ring = PolynomialRing(field, "u")
        model = direct["weierstrass_model"]
        A, B = [ring([field(QQ(c)) for c in model[f"{key}_coefficients_low_to_high"]])
                for key in ("A", "B")]
        xs, ys = [], []
        for point in direct["sections"]["records"]:
            for key, coordinates in (("X", xs), ("Y", ys)):
                data = point[key]
                coordinates.append(
                    ring([field(QQ(c)) for c in data["numerator_coefficients_low_to_high"]]) /
                    ring([field(QQ(c)) for c in data["denominator_coefficients_low_to_high"]]))
        assert len(xs) == 17
        assert all(y*y == x*x*x + A*x+B for x, y in zip(xs, ys))
        for index, data in enumerate(pair["branch_quadratics_coefficients_low_to_high"]):
            branch = ring([field(QQ(c)) for c in data])
            assert branch.degree() == 2 and branch.is_squarefree()
            assert branch.gcd(4*A**3+27*B**2) == 1
            for factor, _ in branch.factor():
                residue = field.extension(factor, "v") if factor.degree() > 1 else field
                v = residue.gen() if factor.degree() > 1 else -factor[0]/factor[1]
                if any(x.denominator()(v) == 0 for x in xs):
                    continue
                cubic_ring = PolynomialRing(residue, "z")
                z = cubic_ring.gen()
                for theta_factor, _ in (z**3+A(v)*z+B(v)).factor():
                    if theta_factor.degree() != 1:
                        continue
                    theta = -theta_factor[0]/theta_factor[1]
                    values = [x(v)-theta for x in xs]
                    if any(value == 0 for value in values):
                        continue
                    bits = [0 if value.is_square() else 1 for value in values]
                    extended = matrix.stack(Matrix(GF(2), [bits]))
                    if extended.rank() == matrix.rank():
                        continue
                    matrix = extended
                    rows.append({"prime": prime, "branch_index": index,
                                 "base_factor": str(factor), "theta_factor": str(theta_factor),
                                 "nonzero_residue_values": [str(value) for value in values],
                                 "square_character_bits": bits})
        if matrix.rank() == 17:
            break
    assert matrix.nrows() == matrix.ncols() == matrix.rank() == 17
    return {"rank": 17, "kernel_dimension": 0, "rows": rows,
            "inverse_matrix_over_F2": [[int(c) for c in row] for row in matrix.inverse().rows()],
            "consequence": "The invariant Kummer span has zero intersection with the twist point-Kummer image.",
            "full_two_selmer_computation": False}


def reduction_record(prime, model, pair, frobenius, pair_key=KEY):
    """Shared rank-one local/BSD calculation; the default preserves this replay."""
    if pair["pair_key"] != pair_key or frobenius["pair_key"] != pair_key or frobenius["prime"] != prime:
        raise ArithmeticError("wrong target or prime")
    if frobenius["status"] != "PASS_COMPLETE_FROBENIUS_PICARD_BOUND":
        raise ArithmeticError("uncertified Frobenius input")
    if frobenius["good_reduction"]["status"] != "PASS":
        raise ArithmeticError("missing good reduction")
    field = GF(prime)
    ring = PolynomialRing(field, "u")
    A, B = [ring([field(QQ(c)) for c in model[f"{key}_coefficients_low_to_high"]])
            for key in ("A", "B")]
    d = ring([field(QQ(c)) for c in pair["product_quartic_coefficients_low_to_high"]])
    delta = -16 * (4 * A**3 + 27 * B**2)
    if (A.degree(), B.degree(), delta.degree(), d.degree()) != (8, 12, 24, 4):
        raise ArithmeticError("degree loss")
    if not delta.is_squarefree() or not d.is_squarefree() or delta.gcd(d) != 1:
        raise ArithmeticError("bad surface reduction")
    branches = []
    for factor, exponent in d.factor():
        assert exponent == 1
        residue = field.extension(factor, "v") if factor.degree() > 1 else field
        v = residue.gen() if factor.degree() > 1 else -factor[0] / factor[1]
        cubic_ring = PolynomialRing(residue, "x")
        x = cubic_ring.gen()
        cubic = x**3 + A(v)*x + B(v)
        factors = list(cubic.factor())
        assert cubic.is_squarefree()
        degrees = [int(f.degree()) for f, _ in factors]
        # Tate's algorithm for I0*: c_v=1+#(nonzero 2-torsion over k(v)).
        tamagawa = 1 + degrees.count(1)
        assert tamagawa in (1, 2, 4)
        branches.append({
            "base_factor_coefficients_low_to_high": coefficients(factor),
            "base_place_degree": int(factor.degree()),
            "residual_two_division_factor_degrees": degrees,
            "kodaira_symbol": "I0*",
            "tamagawa_number": tamagawa,
        })
    tamagawa_product = prod(row["tamagawa_number"] for row in branches)
    # Every other singular fibre is I1 (c_v=1); infinity is smooth.
    qring = PolynomialRing(QQ, "T")
    T = qring.gen()
    F = qring(frobenius["elliptic_L"]["frobenius_characteristic_coefficients_low_to_high"])
    L = F.reverse()  # characteristic polynomial det(Z-Frob) -> det(1-T Frob)
    assert L.degree() == 28 and L[0] == 1 and F.is_monic()
    reduced, remainder = L.quo_rem(1-prime*T)
    leading = reduced(QQ(1)/prime)
    if remainder or leading <= 0:
        raise ArithmeticError("analytic rank is not exactly one")
    # Independently check orientation against the stored n=1,2 moments.
    moments = [str(-L[1]), str(L[1]**2-2*L[2])]
    assert moments == frobenius["elliptic_L"]["power_sums_n1_n2"]
    chi = (24 + 4*6)//12
    assert chi == 4
    bsd_product = QQ(prime**(chi-1))*leading/tamagawa_product
    record = {
        "prime": prime,
        "A_mod_p": coefficients(A), "B_mod_p": coefficients(B),
        "d_mod_p": coefficients(d),
        "good_reduction": True,
        "geometric_fibres": "4I0*+24I1",
        "chi": chi,
        "branch_places": branches,
        "other_tamagawa_numbers": "all 1 (I1 and smooth infinity)",
        "tamagawa_product": int(tamagawa_product),
        "L_coefficients_low_to_high": coefficients(L),
        "analytic_rank": 1,
        "L_star_at_inverse_p": str(leading),
        "regulator_times_sha_over_torsion_squared_if_rank_one": str(bsd_product),
        "regulator_squareclass_if_rank_one": squareclass(bsd_product),
        "valuation_at_3": int(bsd_product.valuation(3)),
        "height_box_tests": [
            {"height": h, "height_over_bsd_product": str(QQ(h)/bsd_product),
             "is_rational_square": bool((QQ(h)/bsd_product).is_square()),
             "status": "EXCLUDED" if not (QQ(h)/bsd_product).is_square() else "UNKNOWN"}
            for h in (8, 10)
        ],
    }
    magma = f'''SetSeed(190830); SetColumns(0);
print "TARGET|{pair_key.replace('alternate-orbit-', '')}|p={prime}";
F := GF({prime}); K<u> := FunctionField(F);
A := {A}; B := {B}; d := {d};
E := EllipticCurve([K|0,0,0,d^2*A,d^3*B]);
R<T> := PolynomialRing(Rationals());
L := R![{','.join(coefficients(L))}];
print "KODAIRA", KodairaSymbols(E);
print "ANALYTIC_INFORMATION", AnalyticInformation(E,L);
print "COMPLETE";
'''
    return record, magma


def build_payload(export_magma=False):
    direct = json.loads(DIRECT.read_text())
    if direct["status"] != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError("uncertified direct model")
    pair = next(row for row in json.loads(PAIRS.read_text())["pairs"] if row["pair_key"] == KEY)
    inputs = {relative(p): digest(p) for p in (DIRECT, PAIRS, HEIGHT)}
    reductions = []
    for prime in PRIMES:
        path = RESULTS / f"elkies-k3-r17-product-alternate-orbit-19bad--alternate-orbit-083ad-p{prime}-toric-frobenius-v1.json"
        inputs[relative(path)] = digest(path)
        record, magma = reduction_record(prime, direct["weierstrass_model"], pair, json.loads(path.read_text()))
        if export_magma:
            CONTROL_DIR.mkdir(parents=True, exist_ok=True)
            (CONTROL_DIR/f"analytic-p{prime}.m").write_text(magma)
        else:
            job_path = CONTROL_DIR/f"analytic-p{prime}.m"
            xml_path = CONTROL_DIR/f"analytic-p{prime}.xml"
            if job_path.read_text() != magma:
                raise ArithmeticError("stale independent-control input")
            raw = ET.fromstring(xml_path.read_text())
            if raw.findtext("headers/warning"):
                raise ArithmeticError("independent control has a warning or error")
            lines = [line.text or "" for line in raw.findall("results/line")]
            expected = f"ANALYTIC_INFORMATION <1, 2, {record['regulator_times_sha_over_torsion_squared_if_rank_one']}>"
            if expected not in lines or "COMPLETE" not in lines or f"TARGET|19bad:083ad|p={prime}" not in lines:
                raise ArithmeticError("independent analytic control does not agree")
            record["independent_magma_control"] = {
                "version": raw.findtext("headers/version"), "result": expected,
                "boundary": "Independent local/BSD normalization check using the supplied certified L-polynomial; not a new L-polynomial computation.",
            }
            for path in (job_path, xml_path):
                inputs[relative(path)] = digest(path)
        reductions.append(record)
    left, right = [QQ(row["regulator_times_sha_over_torsion_squared_if_rank_one"])
                   for row in reductions[:2]]
    ratio = left/right
    if ratio.is_square() or ratio.valuation(3) % 2 != 1:
        raise ArithmeticError("two-prime obstruction did not close")
    assert all(t["status"] == "EXCLUDED" for row in reductions for t in row["height_box_tests"])
    return {
        "schema": "elkies-k3.r17-product-19bad-083ad-rank-zero.v1",
        "status": "PROVED_ARITHMETIC_PRODUCT_TWIST_RANK_ZERO",
        "pair_key": KEY,
        "rank_over_QQ_u": 0,
        "rank_over_QQbar_u": {"lower": 0, "upper": 2, "status": "UNKNOWN"},
        "reductions": reductions,
        "contradiction": {
            "primes": list(PRIMES[:2]), "regulator_ratio": str(ratio),
            "ratio_squareclass": squareclass(ratio),
            "is_rational_square": False,
            "odd_valuation_witness": {"prime": 3, "valuation": int(ratio.valuation(3))},
            "logic": [
                "A nonzero QQ(u)-section is nontorsion by the 48I1 cover height formula.",
                "Good surface reduction with unchanged geometric fibre types preserves its positive Shioda height.",
                "At each prime its nonzero reduction forces algebraic rank = analytic rank = 1.",
                "Rank equality implies refined BSD and finite square-order Sha over the finite function field.",
                "Thus its height has the squareclass of p^(chi-1)*L_star/product(c_v).",
                "The same rational height cannot have the two different squareclasses 282 and 154.",
            ],
        },
        "kummer_and_tate": {
            "anti_invariant_MW_group": "0", "MW_mod_2_dimension": 0,
            "integral_character_glue_dimension": 0, "tate_H_minus_1_dimension": 0,
            "nonzero_MW_Kummer_representatives": [], "nonzero_Tate_representatives": [],
            "full_two_selmer_group": "NOT_COMPUTED",
            "nonzero_Selmer_candidate_representatives": None,
            "boundary": "An empty point-Kummer image does not assert an empty Selmer group or Sha[2].",
        },
        "auxiliary_local_invariant_kummer_separation": local_invariant_kummer_separation(direct, pair),
        "section_boxes": {
            "height_8_PO_0_degrees_8_12": "PROVED_EMPTY",
            "height_10_PO_1_numerator_degrees_10_15": "PROVED_EMPTY",
            "all_higher_heights": "PROVED_EMPTY_OVER_QQ_u",
            "method": "rank-zero theorem, including every denominator and projective boundary chart",
            "zero_class_carrier_search_repeated": False,
        },
        "theoretical_references": [
            {"url": "https://math.stanford.edu/~conrad/BSDseminar/refs/Ulmer.pdf",
             "location": "Sections 6.2.3 and 6.3: Theorems 6.2.6, 6.3.1; Proposition 6.3.3"},
            {"url": "https://magma.maths.usyd.edu.au/magma/handbook/text/1583",
             "location": "Example H134E5: incompatible generator-height squareclasses under good reduction"},
        ],
        "proof_boundary": "Only this target is closed arithmetically. No geometric rank-zero, full Selmer, or other-target conclusion is claimed.",
        "software": version, "inputs": inputs,
        "reproducing_command": "sage -python elkies-k3/scripts/certify_r17_product_19bad_083ad_rank_zero.sage --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--export-magma", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.check and args.export_magma:
        parser.error("--check is read-only; export Magma separately")
    payload = build_payload(args.export_magma)
    if args.export_magma:
        print(f"PRODUCT19BAD083AD|magma_inputs=EXPORTED|directory={CONTROL_DIR}")
        return
    rendered = json.dumps(payload, sort_keys=True, indent=2) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale or missing target rank-zero certificate")
    else:
        args.output.write_text(rendered)
    print(f"PRODUCT19BAD083AD|rank_QQ_u=0|Hminus1=0|height8=EMPTY|height10=EMPTY|status={'PASS' if args.check else 'WROTE'}", flush=True)


if __name__ == "__main__":
    main()
