#!/usr/bin/env sage-python
"""Replay the two RR model/field/order certificates, not a Selmer group.

No nfinit, nfcertify, factorization, reduction or class-group routine is used.
Trace pairings are reconstructed with Newton sums. Backend incompleteness is
a retained execution outcome, not a theorem that the calculation is impossible.
"""
import hashlib
import json
import signal
import sys
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, vector

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
R = PolynomialRing(QQ, "z")
z = R.gen()
CHECKED = {}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    CHECKED[str(p.relative_to(ROOT))] = sha(p)
    return json.loads(p.read_text())


def coords(f):
    return vector(QQ, [f[j] for j in range(6)])


def check_isomorphism(old, new, root):
    assert new.degree() == 6 and new.is_monic()
    assert old(root) % new == 0
    assert matrix(QQ, [coords(root**j % new) for j in range(6)]).det()


def check_order(f, basis, claimed_disc):
    assert f.degree() == 6 and f.is_monic()
    M = matrix(QQ, [coords(b) for b in basis]).transpose()
    inverse = M.inverse()
    one = inverse * coords(R(1))
    assert all(a.denominator() == 1 for a in one)
    # Newton sums reconstruct trace without calling the constructor's
    # multiplication-matrix trace computation or any number-field backend.
    sums = [QQ(6)]
    for k in range(1,6):
        sums.append(-sum(f[6-j]*sums[k-j] for j in range(1,k))-k*f[6-k])
    traces = []
    for bi in basis:
        row = []
        for bj in basis:
            c = coords(bi*bj % f)
            assert all(a.denominator() == 1 for a in inverse*c)
            row.append(sum(c[k]*sums[k] for k in range(6)))
        traces.append(row)
    D = matrix(QQ, traces).det()
    assert D == QQ(claimed_disc) != 0
    assert D == f.discriminant()*M.det()**2
    return abs(ZZ(D)).nbits()


def verify():
    v2 = ART / "det1092_rr_global_pair_v2"
    v3 = ART / "det1092_rr_global_pair_v3"
    protocol = read(v2 / "protocol.json")
    refined_protocol = read(v3 / "protocol.json")
    assert [r["case_index"] for r in protocol["cases"]] == [8,9]
    assert refined_protocol["cases"] == [8,9]
    for path,h in refined_protocol["inputs"].items():
        assert sha(ROOT/path) == h
    for name,p in (("probe_det1092_rr_global_pair.sage",protocol),
                   ("refine_det1092_rr_pair_order.sage",refined_protocol)):
        source = ROOT / "elliptic-curves/cas" / name
        assert sha(source) == p["source_sha256"]
        CHECKED[str(source.relative_to(ROOT))] = sha(source)
    results = []
    for row in protocol["cases"]:
        i = row["case_index"]
        source = ROOT/row["prepared_case"]
        assert sha(source) == row["prepared_case_sha256"]
        prepared = read(source)
        worker = ROOT/row["worker_input"]
        assert sha(worker) == row["worker_input_sha256"]
        eq = read(worker)
        assert set(eq) == {"schema","P","Q"}
        assert (eq["P"],eq["Q"]) == (prepared["reduced_P"],prepared["reduced_Q"])
        P,Q = R(eq["P"]),R(eq["Q"])
        f = 4*P+Q*Q
        # Reuse the two specified irreducibility places, not a prime search.
        prime = {8:79,9:83}[i]
        reduction = PolynomialRing(GF(prime),"z")(f)
        assert reduction.degree() == 6 and reduction.is_irreducible()
        inp = read(v2/("case-%02d-input.json" % i))
        assert f.discriminant() == ZZ(inp["polynomial_discriminant"])
        model = read(v2/("case-%02d-model.json" % i))
        p,q = R(model["P"]),R(model["Q"])
        trans = model["map_to_worker_input"]
        e = QQ(trans["e"])
        a,b,c,d = map(QQ,trans["mobius"])
        H = R(trans["H"])
        assert e and a*d-b*c
        pull = sum(f[j]*(a*z+b)**j*(c*z+d)**(6-j) for j in range(7))
        assert pull == e*e*(4*p+q*q)
        pullQ = sum(Q[j]*(a*z+b)**j*(c*z+d)**(3-j) for j in range(4))
        assert 2*H+pullQ == e*q
        primitive = R(model["primitive_sextic"])
        assert (4*p+q*q)/primitive in QQ
        degrees = []
        order = read(v2/("case-%02d-order.json" % i))
        field = R(order["defining_polynomial"])
        check_isomorphism(primitive,field,R(order["old_root_in_new_field"]))
        degrees.append(check_order(field,[R(b) for b in order["order_basis"]],order["order_discriminant"]))
        binary = read(v3/("case-%02d-binary-order.json" % i))
        G = R(binary["monic_polynomial"])
        assert G == primitive[6]**5*primitive(z/primitive[6])
        degrees.append(check_order(G,[R(b) for b in binary["basis"]],binary["order_discriminant"]))
        assert ZZ(binary["order_discriminant"]) == primitive.discriminant()
        reduced = read(v3/("case-%02d-reduced-generator.json" % i))
        better = R(reduced["defining_polynomial"])
        check_isomorphism(G,better,R(reduced["monic_root_in_reduced_field"]))
        refined = read(v3/("case-%02d-order.json" % i))
        assert R(refined["defining_polynomial"]) == better
        degrees.append(check_order(better,[R(b) for b in refined["order_basis"]],refined["order_discriminant"]))
        unresolved = []
        for directory in (v2,v3):
            cert = read(directory/("case-%02d-certification.json" % i))
            assert cert["status"] == "UNRESOLVED_MAXIMAL_ORDER"
            assert cert["maximal_order_certified"] is False
            assert cert["Selmer_dimension"] is None
            assert cert["unresolved_cofactors"]
            sizes = [ZZ(n).nbits() for n in cert["unresolved_cofactors"]]
            assert sizes == cert["unresolved_cofactor_bits"]
            unresolved.append(sizes)
        results.append({"case_index":i,"arm":row["arm"],
                        "irreducibility_prime":prime,"verified_order_discriminant_bits":degrees,
                        "retained_nfcertify_unresolved_bits":unresolved,
                        "Selmer_dimension":None,"residual_mod_inherited_dimension":None,
                        "residual_mod_known_rational_dimension":None})
    attempts = read(v3/"execution.json")
    assert len(attempts["runs"]) == 6
    for record in attempts["runs"]:
        p = ROOT/record["retained_log"]
        assert sha(p) == record["supervisor"]["log_sha256"]
        CHECKED[str(p.relative_to(ROOT))] = sha(p)
        assert record["supervisor"]["outcome"] == ("backend_failure" if record["version"]==1 else "completed")
    return {"status":"PASS_EXACT_MODELS_FIELDS_AND_INTEGRAL_ORDERS_ONLY",
            "global_comparison_status":"NOT_COMPUTED", "cases":results,
            "boundary":"No maximal-order, class-group, unit-completeness, full supported-squareclass or global Selmer certificate. Recorded backend failures do not prove infeasibility or equality of Selmer dimensions.",
            "inputs":CHECKED,"checker_sha256":sha(Path(__file__))}


if __name__ == "__main__":
    signal.alarm(60)
    result = verify()
    output = ART/"det1092_rr_global_pair_v3/replay.json"
    payload = json.dumps(result,indent=2,sort_keys=True)+"\n"
    if output.exists():
        assert output.read_text() == payload
    else:
        output.write_text(payload)
    print(result["status"])
