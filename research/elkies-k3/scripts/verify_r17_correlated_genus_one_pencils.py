#!/usr/bin/env python3
"""Independent Fraction/integer replay of the 25 complete-pencil exclusion.

No Sage, producer imports, lattice census, point search or parameter cutoff.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import itertools
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1"


def require(value, message):
    if not value:
        raise ValueError(message)


def trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return a


def poly(a):
    return trim([Q(x) for x in a])


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, s):
    return trim([x*s for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    if not a or not b:
        return []
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def power(a, n):
    out = [Q(1)]
    for _ in range(n):
        out = mul(out, a)
    return out


def remainder(a, b):
    a = list(a)
    require(b, "zero polynomial divisor")
    while len(a) >= len(b):
        n, c = len(a)-len(b), a[-1]/b[-1]
        for i, x in enumerate(b):
            a[i+n] -= c*x
        a = trim(a)
    return a


def gcd_poly(a, b):
    while b:
        a,b = b,remainder(a,b)
    return scale(a,1/a[-1]) if a else []


def matmul(a, b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def transpose(a):
    return list(map(list, zip(*a)))


def rank_mod(a, p):
    a = [list(row) for row in a]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank,len(a)) if a[i][col] % p), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col] % p, -1, p)
        a[rank] = [x*inv % p for x in a[rank]]
        for i in range(rank+1,len(a)):
            c = a[i][col]
            a[i] = [(x-c*y) % p for x,y in zip(a[i],a[rank])]
        rank += 1
    return rank


def reduce_matrix(a, p):
    try:
        out = [[x.numerator*pow(x.denominator % p,-1,p) % p for x in row] for row in a]
    except ValueError:
        return None
    return out if rank_mod(out,p) == 5 else None


def projective_image(a, p):
    result = set()
    for u in range(p+1):
        v = [pow(u,k,p) for k in range(5)] if u < p else [0,0,0,0,1]
        w = [sum(x*y for x,y in zip(row,v)) % p for row in a]
        pivot = next((x for x in w if x), None)
        require(pivot is not None, "basepoint in projective branch map")
        inv = pow(pivot,-1,p)
        result.add(tuple(x*inv % p for x in w))
    require(len(result) == p+1, "branch map is not an embedding")
    return result


def source_projection(packet):
    for relative, expected in packet["source_hashes"].items():
        require(sha256((ROOT/relative).read_bytes()).hexdigest() == expected, "source hash mismatch")
    data = ROOT / "elkies-k3/data/fibrations"
    model = json.loads((data/"elkies_2026_published_r17_model.json").read_text())
    require(packet["A"] == model["A_coefficients_low_to_high"] and
            packet["B"] == model["B_coefficients_low_to_high"], "model projection mismatch")
    sections = json.loads((data/"elkies_2026_published_r17_sections.json").read_text())["sections"]
    reconstructed = []
    for i, row in enumerate(sections):
        require(row["basis_index"] == i, "source basis order")
        x = poly(row["x_coefficients_low_to_high"])
        if i == 0:
            y = poly(row["y_coefficients_low_to_high"])
        else:
            chord = row["chord"]
            rx, ry = reconstructed[chord["reference_basis_index"]]
            y = add(ry,mul(poly(chord["slope_coefficients_low_to_high"]),sub(x,rx)))
        reconstructed.append((x,y))
    require(reconstructed == [(poly(r["x"]),poly(r["y"])) for r in packet["basis"]], "generic basis projection mismatch")
    identification = json.loads((ROOT/"artifacts/generated-results/elkies-2026-published-r17-target.json").read_text())["pinned_identification"]
    require(identification["gram_identity_orientation"] == "M^T*Gpub*M=Gpinned", "Gram orientation")
    M = [[Q(x) for x in row] for row in identification["basis_change_matrix"]]
    G = [[Q(x) for x in row] for row in packet["gram"]]
    pinned = [[Q(x) for x in line.split()] for line in (ROOT/"elkies-k3/data/lattice/rank17_gram.txt").read_text().splitlines()
              if line.strip() and not line.lstrip().startswith("#")]
    require(matmul(matmul(transpose(M),G),M) == pinned, "generic Gram projection mismatch")


def verify(packet, pencils, result, check_sources=True):
    require(packet["schema"] == "r17-correlated-genus-one-input-v1", "input schema")
    require(result["schema"] == "r17-correlated-genus-one-result-v1", "result schema")
    if check_sources:
        source_projection(packet)
    G = [[Q(x) for x in row] for row in packet["gram"]]
    expected = []
    for i,j in itertools.combinations(range(17),2):
        for sign in (1,-1):
            if G[i][i]+G[j][j]+2*sign*G[i][j] == 8:
                word = [0]*17
                word[i],word[j] = 1,sign
                expected.append(word)
                break
    expected.sort()
    require(len(expected) == packet["selection_pool_count"] == 25, "selection pool mismatch")
    require(packet["selected_words"] == expected[:64], "frozen word selection mismatch")
    require(len(pencils) == result["pencils"] == 25, "pencil coverage mismatch")
    A,B = poly(packet["A"]),poly(packet["B"])
    basis = [(poly(r["x"]),poly(r["y"])) for r in packet["basis"]]
    require(len(basis) == 17, "basis size")
    for x,y in basis:
        require(power(y,2) == add(add(power(x,3),mul(A,x)),B), "generic section identity")
    matrices = []
    for index,(row,word) in enumerate(zip(pencils,expected)):
        require(row["index"] == index and row["word"] == word, "trace attachment mismatch")
        i,j = [k for k,n in enumerate(word) if n]
        x1,y1 = basis[i]
        x2,y2 = basis[j]
        y2 = scale(y2,word[j])
        den,num = sub(x2,x1),sub(y2,y1)
        require(den, "addition requires distinct x coordinates")
        xn = sub(power(num,2),mul(add(x1,x2),power(den,2)))
        yn = sub(mul(num,sub(mul(x1,power(den,2)),xn)),mul(y1,power(den,3)))
        h,Nx,Ny,M = [poly(row[k]) for k in ("h","Nx","Ny","M0")]
        require(len(h) == 3 and h[-1] == 1 and len(M) <= 4, "trace degree or normalization")
        require(mul(Nx,power(den,2)) == mul(xn,power(h,2)), "trace x group-law identity")
        require(mul(Ny,power(den,3)) == mul(yn,power(h,3)), "trace y group-law identity")
        require(not remainder(add(mul(M,Nx),Ny),power(h,2)), "RR line congruence")
        mat = [[Q(x) for x in r] for r in row["branch_matrix"]]
        require(len(mat) == 5 and all(len(r) == 5 for r in mat), "branch matrix shape")
        rhs = [sub(sub(sub(sub(power(M,4),scale(mul(power(M,2),Nx),6)),scale(mul(M,Ny),8)),scale(power(Nx,2),3)),scale(mul(A,power(h,4)),4)),
               mul(sub(sub(scale(power(M,3),4),scale(mul(M,Nx),12)),scale(Ny,8)),power(h,2)),
               scale(mul(sub(power(M,2),Nx),power(h,4)),6),
               scale(mul(M,power(h,6)),4),power(h,8)]
        for column,numerator in zip(zip(*mat),rhs):
            require(mul(trim(column),power(h,6)) == numerator, "branch coefficient identity")
        q0 = [r[0] for r in mat]
        derivative = [k*q0[k] for k in range(1,5)]
        require(q0[4] and len(gcd_poly(q0,derivative)) == 1,
                "lambda zero does not certify a smooth genus-one member")
        matrices.append(mat)
    remaining = set(itertools.combinations(range(25),2))
    require(result["pair_count"] == len(remaining) == 300, "pair count mismatch")
    require([r["prime"] for r in result["modular_witnesses"]] == [101,103,107], "prime roster")
    stats = []
    for row in result["modular_witnesses"]:
        p = row["prime"]
        images = {}
        for i,mat in enumerate(matrices):
            reduced = reduce_matrix(mat,p)
            if reduced is not None:
                images[i] = projective_image(reduced,p)
        require(sorted(images) == row["usable"], "good-reduction roster mismatch")
        rejected = sorted((i,j) for i,j in remaining if i in images and j in images and not images[i].intersection(images[j]))
        require([list(pair) for pair in rejected] == row["rejected_pairs"], "projective exclusion witness mismatch")
        remaining.difference_update(rejected)
        require(len(remaining) == row["remaining_pairs"], "remaining pair coverage mismatch")
        stats.append({"prime":p,"rejected_pairs":len(rejected)})
    require(not remaining and result["exact_pairs"] == [], "this certificate does not close its rational pair roster")
    return {"status":"PASS_COMPLETE_25_PENCIL_PAIR_EXCLUSION", "pencils":25,"pairs":300,
            "exact_generic_sections":17,"exact_trace_and_branch_identities":25,"prime_witnesses":stats,
            "smooth_genus_one_members_at_lambda_zero":25,
            "scope":"No two distinct pencils in this frozen bank share a rational projective binary quartic; no same-quadratic-cover two-gain construction is supplied."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=DEFAULT)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    packet,pencils,result = [json.loads((args.input_dir/name).read_text()) for name in ("input.json","pencils.json","result.json")]
    require(sha256((args.input_dir/"input.json").read_bytes()).hexdigest() == result["input_sha256"], "input binding mismatch")
    replay = verify(packet,pencils,result)
    replay["elapsed_seconds"] = round(time.monotonic()-started,3)
    replay["inputs"] = {name:sha256((args.input_dir/name).read_bytes()).hexdigest() for name in ("input.json","pencils.json","result.json")}
    if args.output:
        with args.output.open("x") as stream:
            json.dump(replay,stream,indent=2,sort_keys=True)
            stream.write("\n")
    print(json.dumps(replay,sort_keys=True))


if __name__ == "__main__":
    main()
