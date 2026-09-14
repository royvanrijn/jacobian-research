#!/usr/bin/env sage-python
"""Two frozen RR fields: bounded prerequisites for an unconditional descent.

This is an order/model probe, not a Selmer implementation. It never calls BNF
on an uncertified order, and a timeout never supplies an upper bound.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, gcd, matrix, pari, prime_range
from sage.version import version as sage_version

sys.set_int_max_str_digits(0)
CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
INPUTS = ART / "det1092_rr_full_selmer_inputs_v1"
OUT = ART / "det1092_rr_global_pair_v2"
LOCAL = ROOT / "artifacts/local/det1092-rr-global-pair-v2"
CASES = (8, 9)
BOUND = 1000
R = PolynomialRing(QQ, "x")
x = R.gen()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retain(path, value):
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    path.write_text(payload)


def coefficients(f):
    return list(map(str, R(f).list()))


def bits(f):
    return max((abs(ZZ(a)).nbits() for a in R(f)), default=0)


def freeze():
    rows = []
    for i in CASES:
        source = INPUTS / ("case-%02d.json" % i)
        prepared = json.loads(source.read_text())
        payload = INPUTS / prepared["worker_input"]
        rows.append({"case_index": i, "arm": prepared["arm"],
                     "prepared_case": str(source.relative_to(ROOT)),
                     "prepared_case_sha256": sha(source),
                     "worker_input": str(payload.relative_to(ROOT)),
                     "worker_input_sha256": sha(payload)})
    retain(OUT / "protocol.json", {
        "schema": "det1092.rr.global-pair-prerequisites.v1",
        "question": "Compare complete global Sel_2(J), then quotient by the inherited rank17 and known rank18 rational images.",
        "selection": "Exactly the historical first-unlock calibration and the source-only section0 control on the same Curve302 fibre; prior cases09 and08. Selection is retrospective.",
        "cases": rows,
        "limits": {"cases": 2, "wall_seconds_per_case": 300,
                   "rss_bytes": 2147483648, "pari_stack_bytes": 536870912,
                   "trial_prime_bound": BOUND, "point_searches": 0,
                   "remote_arithmetic": 0, "class_group_runs_in_this_probe": 0},
        "stages": ["exact_equation_replay", "model_minimization_at_trial_primes",
                   "partial_order_and_generator_reduction", "maximal_order_certification"],
        "completion_gate": "Only nfcertify=[] certifies the maximal order. Full Selmer additionally requires complete supported squareclasses, certified class/unit contributions, all necessary local maps, and the true/fake comparison. This probe does not implement those later steps.",
        "source_sha256": sha(Path(__file__)),
    })


def run_case(i):
    protocol = json.loads((OUT / "protocol.json").read_text())
    assert protocol["source_sha256"] == sha(Path(__file__))
    row = next(r for r in protocol["cases"] if r["case_index"] == i)
    payload_path = ROOT / row["worker_input"]
    assert sha(payload_path) == row["worker_input_sha256"]
    payload = json.loads(payload_path.read_text())
    assert set(payload) == {"schema", "P", "Q"}
    P, Q = R(payload["P"]), R(payload["Q"])
    assert all(a.denominator() == 1 for a in list(P) + list(Q))
    f = 4 * P + Q * Q
    assert f.degree() == 6 and f.gcd(f.derivative()) == 1
    clock = time.monotonic()

    def emit(stage, value):
        value.update({"case_index": i, "stage": stage,
                      "worker_input_sha256": row["worker_input_sha256"],
                      "protocol_sha256": sha(OUT / "protocol.json"),
                      "elapsed_seconds": time.monotonic() - clock,
                      "software": {"sage": sage_version, "pari": str(pari.version())},
                      "Selmer_dimension": None})
        retain(OUT / ("case-%02d-%s.json" % (i, stage)), value)
        print(json.dumps({k: value[k] for k in value if k in {
            "stage", "elapsed_seconds", "order_discriminant_bits",
            "coefficient_bits", "status", "unresolved_cofactor_bits"}}), flush=True)

    pari.allocatemem(32 * 1024**2, 512 * 1024**2, silent=True)
    D = ZZ(f.discriminant())
    primes = [p for p in prime_range(2, BOUND + 1) if D % p == 0]
    emit("input", {"status": "PASS_EQUATION_INPUT", "coefficient_bits": bits(f),
                   "polynomial_discriminant": str(D),
                   "model_minimization_primes": list(map(int, primes))})

    print("START model_minimization_at_trial_primes", flush=True)
    reduce_model = pari("(C,ps)->{my(m);my(W=hyperellminimalmodel(C,&m,ps));[W,m]}")
    W, trans = reduce_model(pari([pari(P), pari(Q)]), pari(primes))
    p, q = R(W[0]), R(W[1])
    e, H = QQ(trans[0]), R(trans[2])
    a, b, c, d = [QQ(trans[1][j, k]) for j, k in ((0, 0), (0, 1), (1, 0), (1, 1))]
    assert e and a * d - b * c
    PH = sum(P[j] * (a*x+b)**j * (c*x+d)**(6-j) for j in range(7))
    QH = sum(Q[j] * (a*x+b)**j * (c*x+d)**(3-j) for j in range(4))
    assert QH + 2 * H == e * q
    assert PH - H*H - QH*H == e*e * p
    f = 4*p + q*q
    content = gcd(ZZ(v) for v in f)
    primitive = f / content
    if primitive.leading_coefficient() < 0:
        primitive = -primitive
    assert primitive.degree() == 6
    emit("model", {"status": "PASS_EXACT_MODEL_MAP", "P": coefficients(p),
                   "Q": coefficients(q), "primitive_sextic": coefficients(primitive),
                   "coefficient_bits": bits(f),
                   "map_to_worker_input": {"e": str(e), "mobius": list(map(str, (a,b,c,d))), "H": coefficients(H)},
                   "boundary": "Minimal only at the supplied trial primes. Primitive sextic defines the descent field, not necessarily the original curve twist."})

    print("START partial_order_and_generator_reduction", flush=True)
    # With the explicit bound, this does not request full factorization.
    init = pari.nfinit(pari([pari(primitive), pari(BOUND)]), 1)
    if len(init) == 2:
        nf, old_root = init[0], init[1]
    else:
        nf, old_root = init, pari(x).Mod(pari(primitive))
    new_poly = R(nf.nf_get_pol())
    old_root_poly = R(old_root.lift())
    assert primitive(old_root_poly) % new_poly == 0
    root_powers = matrix(QQ, [list((old_root_poly**j % new_poly).list()) + [0]*(6-len((old_root_poly**j % new_poly).list())) for j in range(6)])
    assert root_powers.det() != 0
    order_disc = ZZ(nf.disc())
    order_basis = [R(b) for b in nf.nf_get_zk()]
    emit("order", {"status": "PARTIAL_ORDER_REQUIRES_CERTIFICATION",
                   "defining_polynomial": coefficients(new_poly),
                   "old_root_in_new_field": coefficients(old_root_poly),
                   "order_basis": [coefficients(b) for b in order_basis],
                   "order_discriminant": str(order_disc),
                   "order_discriminant_bits": abs(order_disc).nbits(),
                   "coefficient_bits": bits(new_poly),
                   "signature": list(map(int, nf.nf_get_sign())),
                   "boundary": "The displayed discriminant is an order discriminant until nfcertify passes."})
    LOCAL.mkdir(parents=True, exist_ok=True)
    binary = LOCAL / ("case-%02d.nf" % i)
    if binary.exists():
        raise FileExistsError(binary)
    pari.writebin(str(binary), nf)
    print("START maximal_order_certification", flush=True)
    unresolved = list(pari.nfcertify(nf))
    emit("certification", {"status": "PASS_CERTIFIED_MAXIMAL_ORDER" if not unresolved else "UNRESOLVED_MAXIMAL_ORDER",
                           "unresolved_cofactors": list(map(str, unresolved)),
                           "unresolved_cofactor_bits": [ZZ(n).nbits() for n in unresolved],
                           "maximal_order_certified": not unresolved,
                           "class_group_complete": False, "units_complete": False,
                           "supported_squareclasses_complete": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--case", type=int, choices=CASES)
    args = parser.parse_args()
    if args.freeze:
        freeze()
    elif args.case is not None:
        run_case(args.case)
    else:
        parser.error("choose --freeze or --case")
