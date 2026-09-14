#!/usr/bin/env sage-python
"""Remove the artificial monic index using a verified binary-sextic order.

Only the already frozen pair is admitted. The explicit order is used for LLL
generator reduction, never asserted maximal. Later nfinit/nfcertify must still
certify maximality. The old models and all previous outcomes are retained.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, matrix, pari, vector
from sage.version import version as sage_version

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
PREVIOUS = ART / "det1092_rr_global_pair_v2"
OUT = ART / "det1092_rr_global_pair_v3"
R = PolynomialRing(QQ, "x")
x = R.gen()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retain(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def coeff(f):
    return list(map(str, R(f).list()))


def coords(f):
    return vector(QQ, [f[j] for j in range(6)])


def freeze():
    inputs = [PREVIOUS / "protocol.json"] + [PREVIOUS / ("case-%02d-model.json" % i) for i in (8,9)]
    retain(OUT / "protocol.json", {
        "schema": "det1092.rr.binary-order-refinement.v1",
        "cases": [8, 9], "inputs": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        "source_sha256": sha(Path(__file__)),
        "purpose": "Remove the artificial index from making a nonmonic sextic monic before reattempting maximal-order certification on the same fields.",
        "limits": {"wall_seconds_per_case": 300, "rss_bytes": 2147483648,
                   "pari_stack_bytes": 536870912, "trial_prime_bound": 1000,
                   "polredbest_calls_per_case": 1, "new_curves": 0, "class_group_runs": 0},
        "boundary": "Ring closure and trace discriminants certify an integral order, not maximality. Only a successful separate maximal-order test permits class/unit work."})


def run(i):
    protocol = json.loads((OUT / "protocol.json").read_text())
    assert protocol["source_sha256"] == sha(Path(__file__))
    for p, h in protocol["inputs"].items():
        assert sha(ROOT / p) == h
    data = json.loads((PREVIOUS / ("case-%02d-model.json" % i)).read_text())
    f = R(data["primitive_sextic"])
    a = f[6]
    G = R(a**5 * f(x/a))
    assert G.is_monic() and all(v.denominator() == 1 for v in G)
    theta = x / a
    B = [R(1)] + [R(sum(f[6-j] * theta**(k-j) for j in range(k))) for k in range(1,6)]
    base = matrix(QQ, [coords(b) for b in B]).transpose()
    inv = base.inverse()
    products = []
    for j in range(6):
        row = []
        for k in range(6):
            c = inv * coords(B[j] * B[k] % G)
            assert all(v.denominator() == 1 for v in c)
            row.append(list(map(str, c)))
        products.append(row)
    matrices = [matrix(ZZ, products[j]).transpose() for j in range(6)]
    trace = matrix(ZZ, 6, 6, lambda j,k: (matrices[j]*matrices[k]).trace())
    assert trace.det() == f.discriminant() == G.discriminant()*base.det()**2
    start = time.monotonic()

    def emit(stage, record):
        record.update({"case_index": i, "stage": stage, "elapsed_seconds": time.monotonic()-start,
                       "protocol_sha256": sha(OUT/"protocol.json"),
                       "Selmer_dimension": None, "software": {"sage": sage_version, "pari": str(pari.version())}})
        retain(OUT / ("case-%02d-%s.json" % (i,stage)), record)
        print(json.dumps({k:v for k,v in record.items() if k in ["stage","status","elapsed_seconds","coefficient_bits","order_discriminant_bits","unresolved_cofactor_bits"]}),flush=True)

    emit("binary-order", {"status": "PASS_INTEGRAL_ORDER", "monic_polynomial": coeff(G),
        "basis": [coeff(b) for b in B], "multiplication_table": products,
        "trace_pairing": [[str(v) for v in row] for row in trace.rows()],
        "order_discriminant": str(trace.det()), "order_discriminant_bits": abs(ZZ(trace.det())).nbits()})
    pari.allocatemem(32*1024**2,512*1024**2,silent=True)
    print("START polredbest_on_verified_order",flush=True)
    reduced, alpha = pari.polredbest(pari([pari(G),pari([pari(b) for b in B])]),1)
    reduced = R(reduced)
    alpha = R(alpha.lift())
    assert G(alpha) % reduced == 0
    assert matrix(QQ, [coords(alpha**j % reduced) for j in range(6)]).det()
    emit("reduced-generator", {"status": "PASS_EXACT_FIELD_ISOMORPHISM", "defining_polynomial": coeff(reduced),
        "monic_root_in_reduced_field": coeff(alpha), "coefficient_bits": max(abs(ZZ(v)).nbits() for v in reduced),
        "polynomial_discriminant_bits": abs(ZZ(reduced.discriminant())).nbits()})
    print("START nfinit_partial_order",flush=True)
    nf = pari.nfinit(pari([pari(reduced),pari(1000)]))
    emit("order", {"status": "PARTIAL_ORDER_REQUIRES_CERTIFICATION", "defining_polynomial": coeff(R(nf.nf_get_pol())),
        "order_basis": [coeff(R(b)) for b in nf.nf_get_zk()], "order_discriminant": str(nf.disc()),
        "order_discriminant_bits": abs(ZZ(nf.disc())).nbits()})
    print("START nfcertify",flush=True)
    pending = list(pari.nfcertify(nf))
    emit("certification", {"status": "PASS_CERTIFIED_MAXIMAL_ORDER" if not pending else "UNRESOLVED_MAXIMAL_ORDER",
        "maximal_order_certified": not pending, "unresolved_cofactors": list(map(str,pending)),
        "unresolved_cofactor_bits": [ZZ(n).nbits() for n in pending], "class_group_complete": False,
        "units_complete": False, "supported_squareclasses_complete": False})


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--freeze",action="store_true")
    ap.add_argument("--case",type=int,choices=(8,9))
    args = ap.parse_args()
    if args.freeze:
        freeze()
    elif args.case is not None:
        run(args.case)
    else:
        ap.error("choose --freeze or --case")
