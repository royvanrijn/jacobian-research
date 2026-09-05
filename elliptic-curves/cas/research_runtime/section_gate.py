"""Theorem/regulator admission for direct invocations of section solvers."""
from fractions import Fraction
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

from .regulator import Surface


def surface_from_export(export, root):
    """Recover exact characteristic-zero scope from current or retained exports.

    Reduced coefficients alone cannot identify a rational surface. Older
    exports retain their exact input files; verify those before resolving labels.
    """
    documents = []
    root = Path(root).resolve()
    for name, expected in export.get("inputs", {}).items():
        path = (root/name).resolve()
        if not path.is_relative_to(root) or sha256(path.read_bytes()).hexdigest() != expected:
            raise ArithmeticError("section export input changed")
        if path.suffix == ".json":
            documents.append(json.loads(path.read_text()))
    if "exact_surface" in export:
        return Surface(**export["exact_surface"])
    models = []
    for document in documents:
        for row in (document, document.get("weierstrass_model", {}), document.get("representative", {})):
            if all(key in row for key in ("A_coefficients_low_to_high", "B_coefficients_low_to_high")):
                models.append(row)
    if len(models) != 1:
        return None
    model = models[0]
    candidate = export.get("candidate", {"kind": "ladder", "key": export.get("label")})
    kind, key = candidate["kind"], candidate["key"]

    def integral(values):
        values = tuple(Fraction(str(c)) for c in values)
        denominator = lcm(*(c.denominator for c in values))
        return tuple(int(c*denominator**2) for c in values)

    def branch(label):
        for document in documents:
            for row in document.get("bisections", []):
                if row.get("label") == label or str(row.get("lattice_orbit_mask")) == label:
                    values = (row["branch"]["numerator_coefficients"] if "branch" in row else
                              row["residual_chord"]["q_coefficients"])
                    return integral(values)
            for row in document.get("construction", {}).get("records", []):
                if row.get("label") == label:
                    return integral(row["branch_polynomial_q_coefficients_low_to_high"])
            for fibre in document.get("fibres", []):
                for row in fibre.get("records", []):
                    if row.get("label") == label and "branch_quadratic_coefficients_low_to_high" in row:
                        # The legacy ladder used the literal rational model.
                        return tuple(row["branch_quadratic_coefficients_low_to_high"])
        return None
    if kind in ("product", "direct_product"):
        labels = str(key).split(":")
        if len(labels) != 2:
            return None
        left, right = map(branch, labels)
        if left is None or right is None:
            return None
        d = [sum(Fraction(str(a))*Fraction(str(b)) for i, a in enumerate(left)
                 for j, b in enumerate(right) if i+j == k) for k in range(len(left)+len(right)-1)]
    else:
        d = branch(str(key))
    return None if d is None else Surface(model["A_coefficients_low_to_high"], model["B_coefficients_low_to_high"], d)


def guard_export(export, root, *, reduction_only=False, limits=None):
    """Run before equation solving, compilation or coefficient enumeration."""
    from production_search_gates import function_field_gate_record
    surface = surface_from_export(export, root)
    if surface is None:
        return {"status": "UNKNOWN", "reason": "exact rational surface unavailable in this legacy export",
                "mathematical_exclusion": False, "search_eligible": True}
    gate = function_field_gate_record(surface=surface, target_rank=1, search_limits=limits or {"gate_checks": 1})
    purpose = "finite_field_proof" if reduction_only or export.get("purpose") == "finite_field_proof" else "rational_section_search"
    if purpose == "rational_section_search" and not gate["search_budget_gate"]["bounded_search_authorized"]:
        raise SystemExit("EXCLUDED_BY_THEOREM before section solving: " +
                         ", ".join(gate["theorem_pruning"]["theorems"]))
    return {"purpose": purpose, "surface_key": surface.key, "gate": gate}
