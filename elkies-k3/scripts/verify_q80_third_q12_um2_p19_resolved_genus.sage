#!/usr/bin/env sage -python
"""Verify genus one for two fibres of the resolved p=19 third-q12 pencil."""

import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.env import SAGE_SHARE
from sage.interfaces.singular import singular


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / (
    "artifacts/generated-results/"
    "q80-third-q12-um2-p19-resolved-pencil.json"
)
OUTPUT = ROOT / (
    "artifacts/generated-results/"
    "q80-third-q12-um2-p19-resolved-genus.json"
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(INPUT.read_text())
if payload["schema"] != "elkies-k3.q80-third-q12-resolved-pencil-modp2.v1":
    raise ValueError("unexpected resolved-pencil schema")
if payload["status"] != (
    "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC"
):
    raise ValueError("resolved pencil is not certified")

base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "r_modulus")
r_modulus = modulus_ring.gen()
finite = GF(19**2, "r", modulus=r_modulus**2 + 12 * r_modulus + 3)
plane_ring = PolynomialRing(finite, names=("W", "x"))
W, x = plane_ring.gens()


def field_element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * finite.gen()


terms = payload["moving_equation"]["terms_T_W_x_coefficient_1_r"]


def specialize(parameter):
    answer = plane_ring.zero()
    for t_degree, w_degree, x_degree, coordinates in terms:
        answer += (
            finite(parameter) ** int(t_degree)
            * field_element(coordinates)
            * W**int(w_degree)
            * x**int(x_degree)
        )
    return answer


# Singular's algebraic finite-field parameter has the fixed name ``a``.
# The upstream hnoether ``extdevelop`` procedure also uses ``a`` as a local
# ideal name, which breaks on a nonsplit branch.  Apply the same fail-closed
# compatibility patch as the certified CM24 Brill--Noether worker.
singular.set_ring(plane_ring._singular_())
singular.lib("brnoeth.lib")
hnoether_source = Path(SAGE_SHARE) / "singular/LIB/hnoether.lib"
hnoether_text = hnoether_source.read_text()
procedure_start = hnoether_text.index("proc extdevelop (list l, int Exaktheit)")
procedure_end = hnoether_text.index("\nexample\n", procedure_start)
patched_extdevelop = hnoether_text[procedure_start:procedure_end]
renames = {
    "ideal a=hole(lastrow);": "ideal q80row=hole(lastrow);",
    "else { ideal a=lastrow; }": "else { ideal q80row=lastrow; }",
    "a[Q]=delt;": "q80row[Q]=delt;",
    "a[Q+1]=x;": "q80row[Q+1]=x;",
    "lastrow=zurueck(a);": "lastrow=zurueck(q80row);",
    "else { lastrow=a; }": "else { lastrow=q80row; }",
}
for old, new in renames.items():
    if patched_extdevelop.count(old) != 1:
        raise RuntimeError("unexpected hnoether extdevelop source")
    patched_extdevelop = patched_extdevelop.replace(old, new)
singular.eval("kill extdevelop;")
singular.eval(patched_extdevelop)

records = []
for record_index, parameter in enumerate((1, 7), 1):
    equation = specialize(parameter)
    factorization = tuple(equation.factor())
    irreducible = len(factorization) == 1 and int(factorization[0][1]) == 1
    if not irreducible:
        raise ArithmeticError(f"moving fibre T={parameter} is reducible")
    singular_equation = equation._singular_()
    curve_name = f"Q80CURVE{record_index}"
    singular.eval(
        "printlevel=-1; "
        f"list {curve_name}=Adj_div({singular_equation.name()});"
    )
    print(
        "ADJDEBUG",
        singular.eval(f"size({curve_name})"),
        flush=True,
    )
    genus = int(singular.eval(f"{curve_name}[2][2]"))
    if genus != 1:
        raise ArithmeticError(
            f"moving fibre T={parameter} has normalization genus {genus}"
        )
    records.append(
        {
            "new_base": parameter,
            "degree_W": int(equation.degree(W)),
            "degree_x": int(equation.degree(x)),
            "irreducible": True,
            "normalization_genus": genus,
        }
    )

output = {
    "schema": "elkies-k3.q80-third-q12-resolved-genus-modp2.v1",
    "status": "PASS_EXACT_THIRD_Q12_GENUS_ONE_MOD19_QUADRATIC",
    "specialization": {
        "u": "-2",
        "prime": 19,
        "extension_modulus": "r^2+12*r+3",
    },
    "fibres": records,
    "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": (
        "Two independent members of the exact resolved moving pencil are "
        "irreducible with normalization genus one. This validates the local "
        "quotient as a genus-one pencil at the pinned modular specialization; "
        "it does not yet supply its minimal Jacobian or child fibre marking."
    ),
    "reproduce": (
        "sage -python "
        "elkies-k3/scripts/verify_q80_third_q12_um2_p19_resolved_genus.sage"
    ),
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12GENUS|prime=19|extension=2|T=1,7|"
    "irreducible=1,1|genus=1,1|"
    "status=PASS_EXACT_THIRD_Q12_GENUS_ONE_MOD19_QUADRATIC",
    flush=True,
)
