#!/usr/bin/env sage-python
"""Build one target-blind determinant-1092 low-shell wave geometry."""
import argparse
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import ZZ, matrix, pari

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
sys.path.insert(0, str(CAS))
import det1092_low_shell_cascade as control
import certify_compact_r17_candidates as cert

mapper = SourceFileLoader("low_shell_mapper", str(CAS / "prepare_fresh_r17_pari_batch.sage")).load_module()
mapper.mapping = SourceFileLoader("low_shell_factor_free_mapping", str(CAS / "factor_free_pari_mapping.sage")).load_module().mapping
mapper.geometry = SourceFileLoader("low_shell_geometry", str(CAS / "prospective_half_lattice_v3.sage")).load_module()


def main(index, wave):
    control.configure(index, wave)
    p, rank, out = control.protocol(), control.ROW["initial_rank"], control.D / "maps.json"
    if out.exists(): raise FileExistsError("preserve low-shell geometry")
    seed = cert.read(control.SEED)
    model = tuple(map(cert.F, seed["curve"])); points = tuple(tuple(map(cert.F, point)) for point in seed["points"])
    gram, asymmetry = mapper.geometry.canonical_height_gram(model, points)
    rounded = mapper.geometry.rounded_gram(gram, 1_000_000)
    g = matrix(ZZ, rounded); change = matrix(ZZ, pari(g).qflllgram()).transpose(); inverse = change.inverse()
    if abs(change.det()) != 1: raise ArithmeticError("non-unimodular low-shell metric transport")
    reduced = change * g * change.transpose(); oracle = mapper.geometry.CosetOracle(reduced.rows()); sample = []
    data = {"status":"RUNNING_SAMPLE", "protocol_hash":control.digest(p), "metric_gram":[[str(value) for value in row] for row in gram], "maximum_gram_asymmetry":str(asymmetry), "rounded_gram":[list(map(int,row)) for row in g.rows()], "change_of_basis":[list(map(int,row)) for row in change.rows()], "reduced_gram":[list(map(int,row)) for row in reduced.rows()], "sample":sample, "rows":[]}
    control.checkpoint(out, data)
    for mask in control.masks(p):
        residue = matrix(ZZ, 1, rank, [(mask >> bit) & 1 for bit in range(rank)])
        target = [int(value) % 2 for value in (residue * inverse).row(0)]
        norm, reduced_word, error = oracle.solve(target)
        word = list(map(int, (matrix(ZZ, 1, rank, reduced_word) * change).row(0)))
        if any((word[bit] - ((mask >> bit) & 1)) % 2 for bit in range(rank)) or sum(word[i] * g[i,j] * word[j] for i in range(rank) for j in range(rank)) != norm:
            raise ArithmeticError("exact low-shell parity transport differs")
        sample.append({"parity":mask,"representative":word,"metric_norm":norm,"cvp_error":error,"reduced_representative":list(reduced_word)})
    data["centres"] = sorted(sample, key=lambda row:(row["metric_norm"], row["parity"]))[:49]
    data["status"] = "RUNNING_MAPS"; control.checkpoint(out, data)
    mapper.pari.allocatemem(256000000, silent=True)
    for centre in data["centres"]:
        data["rows"].append(mapper.mapping(model, points, centre)); control.checkpoint(out, data)
    data["status"] = "COMPLETE_DECLARED_MAPS"; control.checkpoint(out, data)
    print("FROZEN LOW-SHELL 49 MAPS", control.ROW["id"], flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--index",type=int,required=True); parser.add_argument("--wave",type=int,required=True); args=parser.parse_args(); main(args.index,args.wave)
