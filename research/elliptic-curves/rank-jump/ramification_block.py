#!/usr/bin/env python3
"""Global reciprocity forces a single coherent new-prime ramification block."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import affine_selmer as af

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"RAMIFICATION_BLOCK_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_ramification_block_v1.json"


def build(check=False):
    inp=r.read(lc.INPUT);collision=r.read(lc.OUTPUT)
    raw=r.read(af.INPUT);affine=r.read(af.OUTPUT)
    rows=[]
    for u in r.read(PROTOCOL)["parameters"]:
        old=next(x for x in inp["rows"] if int(x["parameter_u"])==u)
        local=next(x for x in collision["local_rows"] if x["u"]==u)
        chars=next(x for x in raw["cases"] if x["u"]==u)
        old_constraints=[v for place in local["local_places"] if not place["new"]
                         for v in place["constraint_basis"]]
        old_constraints.extend(r.transpose(old["real_local_condition"]["known_span_quotient_rows"]))
        test_space=lc.orthogonal(old_constraints,20)
        checks=[]
        for event in local["root_character_checks"]:
            assert event["eligible"] and event["matches"]
            p=event["prime"]
            rec=next(x for x in chars["local"] if x["place"]==p)
            assert p>2
            # At odd primes the signature order is (valuation, unit) per ideal.
            beta=rec["class_signature_rows"][:20]
            eta=rec["class_signature_rows"][20]
            k=len(rec["prime_decomposition"])
            assert len(eta)==2*k
            assert all(all(row[2*j]==0 for j in range(k)) for row in beta)
            pattern=eta[::2]
            assert any(pattern)
            ram_bits=[]
            for gen in rec["point_signature_rows"]:
                valuations=gen[::2]
                assert valuations==[0]*k or valuations==pattern
                bit=int(valuations==pattern)
                ram_bits.append(bit)
                pairing=r.pack([sum(b[2*j+1]*gen[2*j] for j in range(k))%2 for b in beta])
                assert pairing==(event["predicted_constraint"] if bit else 0)
            assert any(ram_bits)
            checks.append({"prime":p,"root_character":event["predicted_constraint"],
                           "local_ramification_pattern":pattern,
                           "point_basis_ramification_bits":ram_bits,
                           "local_pairing_identity_verified":True})
        matrix=[r.pack([(b&c["root_character"]).bit_count()%2 for c in checks])
                for b in test_space]
        kernel=lc.orthogonal(matrix,len(checks))
        assert kernel==[(1<<len(checks))-1]
        row=next(x for x in affine["cases"] if x["u"]==u)
        realized=row["affine_solution"]["consistent"]
        rows.append({"u":u,"new_prime_order":[x["prime"] for x in checks],
                     "test_space_basis":test_space,"test_space_dimension":len(test_space),
                     "local_checks":checks,"reciprocity_matrix_packed_rows":matrix,
                     "reciprocity_matrix_rank":r.rank(matrix),"allowed_ramification_basis":kernel,
                     "full_selmer_new_ramification_dimension_upper_bound":1,
                     "full_selmer_new_ramification_dimension_exact":1 if realized else None,
                     "realization_source":str(af.OUTPUT.relative_to(r.ROOT)) if realized else None,
                     "exact_dimension_status":"PROVED_ONE" if realized else "UNKNOWN_ZERO_OR_ONE"})
    bindings={str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
              for p in (PROTOCOL,Path(__file__),lc.INPUT,lc.OUTPUT,af.INPUT,af.OUTPUT)}
    out={"schema":"rank-jump.ramification-block.v1","bindings":bindings,"rows":rows,
         "scope":"Full Selmer classes satisfy the ramification constraint by global reciprocity. No full Selmer dimension or curve-rank upper bound is inferred."}
    if check:
        assert r.read(OUTPUT)==out
        print("PASS full-Selmer ramification block certificate")
    else:r.write_new(OUTPUT,out)
    for row in rows:
        print(row["u"],len(row["new_prime_order"]),row["test_space_dimension"],
              row["reciprocity_matrix_rank"],row["exact_dimension_status"])


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("build","check"))
    build(parser.parse_args().mode=="check")
