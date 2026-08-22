#!/usr/bin/env sage -python
"""Derive actual E8 coefficient floors for the q=8 generic RR basis.

The q=8 generic ambient is built from the P1 chord

    m=(y-y(P1))/(x-x(P1)).

On the actual II* chart ``u=1/t, X=u^4*x, Y=u^6*y`` the already certified
q=6 calculation gives ``m=u^-2*Q``, with ``Q`` a unit at the E8 singular
point.  The q=8 E8 module is the actual resolved complete module

    u^9*(u^2, X, Y).

For a generic-basis monomial ``x^a*m^b``, this supplies the sharp minimal
base exponent at E8:

* if ``a=0``, ``u^(11+2b)*m^b = u^9*(u^2*Q^b)`` lies in the module;
* if ``a>=1``, ``u^(9+4a+2b)*x^a*m^b = u^9*X^a*Q^b`` lies in the module.

This is an actual resolved E8 condition, not a complete q=8 global ambient:
the E7 gluing, finite smooth conditions, and compatible global base bounds
remain to be derived.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
MODULE = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-complete-module.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-ambient-weights.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--module", type=Path, default=MODULE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
module = json.loads(args.module.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert module["status"] == "PASS_EXACT_Q8_E8_COMPLETE_MODULE"
assert ambient["chord"] == "m=(y-y(P1))/(x-x(P1))"
assert ambient["dimension"] == 18
assert module["module"]["q8_E8_module"] == "u^9*(u^2,X,Y)"
assert module["complete_ideal"]["generators"] == ["Y", "X", "u^2"]
assert module["complete_ideal"]["quotient_basis"] == ["1", "u"]

weights = []
for entry in ambient["basis"]:
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    if a == 0:
        exponent = 11+2*b
        local_representative = "u^2*Q^{}".format(b)
        ideal_generator = "u^2"
    else:
        exponent = 9+4*a+2*b
        local_representative = "X^{}*Q^{}".format(a, b)
        ideal_generator = "X"
    # In u^e*x^a*m^b = u^(e-4a-2b) X^a Q^b, the displayed exponent
    # gives u^11 Q^b=u^9*(u^2 Q^b) for a=0 and u^9 X^a Q^b for a>0.
    # The parenthesized factors lie in (u^2,X,Y).
    residual_u_exponent = exponent-4*a-2*b
    assert residual_u_exponent == (11 if a == 0 else 9)
    weights.append({
        "basis": {"kind": entry["kind"], "x_power": a, "m_power": b},
        "minimal_u_power": exponent,
        "trivialized_representative": local_representative,
        "complete_ideal_generator": ideal_generator,
    })

assert [row["minimal_u_power"] for row in weights[:10]] == list(range(11, 30, 2))
assert [row["minimal_u_power"] for row in weights[10:]] == list(range(13, 28, 2))

payload = {
    "schema": "elkies-k3.h92-q8-e8-ambient-weights.v1",
    "status": "PASS_EXACT_Q8_E8_AMBIENT_WEIGHTS",
    "inputs": {
        "generic_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "actual_e8_module": {"path": str(args.module.relative_to(ROOT)), "sha256": digest(args.module)},
    },
    "actual_chart": {
        "coordinates": "u=1/t, X=u^4*x, Y=u^6*y",
        "marked_chord": "m=u^-2*Q, Q=(Y-Y(P1))/(X-X(P1))",
        "Q_at_E8_singular_point": "unit",
        "target_module": "u^9*(u^2,X,Y)",
    },
    "basis_weight_floors": weights,
    "sharpness_certificate": {
        "for_m_powers": (
            "At one lower u-power the factor after u^9 is u*Q^b. "
            "Because Q is a unit and R/(u^2,X,Y) has nonzero class u, it "
            "does not lie in (u^2,X,Y)."
        ),
        "for_x_m_powers": (
            "At one lower u-power the expression is u^8*X*Q^b. Since Q "
            "is a unit and X is nonzero modulo u on the E8 surface germ, "
            "it is not in u^9*R and hence not in u^9*(u^2,X,Y)."
        ),
    },
    "conclusion": (
        "For each listed q=8 generic-basis generator x^a*m^b, multiplying "
        "by u to the displayed minimal power gives an element of the actual "
        "resolved E8 module. These are exact necessary E8 coefficient floors."
    ),
    "boundary": (
        "This gives only the E8 end of a q=8 global coefficient ambient. It "
        "does not establish compatible global functions, resolved E7 gluing, "
        "finite smooth conditions, a complete condition matrix, or h0=2."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E8WEIGHTS|basis=18|m_floors=11..29|xm_floors=13..27|"
    "status=PASS_EXACT_Q8_E8_AMBIENT_WEIGHTS",
    flush=True,
)
