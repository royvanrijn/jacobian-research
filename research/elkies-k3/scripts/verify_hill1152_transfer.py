#!/usr/bin/env python3
"""Independent rational-polynomial/modular replay of the Hill1152 pilot.

Uses retained standard-library polynomial primitives, not Sage or the producer.
Replays all trace labels, branch matrices, conic points and finite local screens.
"""
import argparse
import ast
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import itertools
import json
from math import gcd, isqrt, lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


V = module("hill_pencil_polynomials", "verify_r17_correlated_genus_one_pencils.py")
M = module("hill_mestre_polynomials", "verify_r17_mestre_shared_twist.py")
require = V.require


def parse_poly(text):
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return [Q(node.value)] if node.value else []
        if isinstance(node, ast.Name) and node.id == "U":
            return [Q(0), Q(1)]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return V.scale(visit(node.operand), -1)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return V.add(a, b)
            if isinstance(node.op, ast.Sub):
                return V.sub(a, b)
            if isinstance(node.op, ast.Mult):
                return V.mul(a, b)
            if isinstance(node.op, ast.Div):
                require(len(b) == 1 and b[0], "nonconstant coefficient divisor")
                return V.scale(a, 1/b[0])
            if isinstance(node.op, ast.Pow):
                require(len(b) == 1 and b[0].denominator == 1 and 0 <= b[0] <= 20, "exponent")
                return V.power(a, int(b[0]))
        raise ValueError("unsupported polynomial source syntax")
    return visit(ast.parse(text.replace("^", "**"), mode="eval").body)


def rf(a, b=None):
    return V.trim(a), V.trim([Q(1)] if b is None else b)


def fadd(a, b):
    return rf(V.add(V.mul(a[0], b[1]), V.mul(b[0], a[1])), V.mul(a[1], b[1]))


def fneg(a):
    return rf(V.scale(a[0], -1), a[1])


def fsub(a, b):
    return fadd(a, fneg(b))


def fmul(a, b):
    return rf(V.mul(a[0], b[0]), V.mul(a[1], b[1]))


def fdiv(a, b):
    require(b[0], "zero rational-function divisor")
    return rf(V.mul(a[0], b[1]), V.mul(a[1], b[0]))


def fsquare(a):
    return fmul(a, a)


def fequal(a, b):
    return V.mul(a[0], b[1]) == V.mul(b[0], a[1])


def point_sum(P, Q0):
    x1, y1 = P
    x2, y2 = Q0
    slope = fdiv(fsub(y2, y1), fsub(x2, x1))
    x = fsub(fsub(fsquare(slope), x1), x2)
    return x, fsub(fmul(slope, fsub(x1, x)), y1)


def mobius(f, c, weight=0):
    # u^weight f(c+1/u), with f represented by rational polynomials.
    def hom(p):
        n = len(p)-1
        out = []
        for i, a in enumerate(p):
            out = V.add(out, V.scale(V.mul([Q(0)]*(n-i)+[Q(1)],
                                          V.power([Q(1), c], i)), a))
        return out
    num, den = f
    exponent = weight-len(num)+len(den)
    n, d = hom(num), hom(den)
    if exponent >= 0:
        n = [Q(0)]*exponent+n
    else:
        d = [Q(0)]*(-exponent)+d
    return rf(n, d)


def source_basis(data, c):
    basis = []
    for row in data["neighbor_sections"]:
        x = rf(parse_poly(row["x_numerator"]), parse_poly(row["x_denominator"]))
        y = rf(parse_poly(row["y_numerator"]), parse_poly(row["y_denominator"]))
        basis.append((mobius(x, c, 4), mobius(y, c, 6)))
    return basis


def trace_from_word(word, basis):
    selected = [(basis[i][0], basis[i][1] if n == 1 else fneg(basis[i][1]))
                for i, n in enumerate(word) if n]
    require(len(selected) in (1, 2), "unsupported word support")
    return selected[0] if len(selected) == 1 else point_sum(*selected)


def check_frame(row, word, basis, A, B):
    require(row["word"] == word, "word label")
    h, Nx, Ny, m = [V.poly(row[k]) for k in ("h", "Nx", "Ny", "M0")]
    h2, h4, h6 = [V.power(h, n) for n in (2, 4, 6)]
    tx, ty = trace_from_word(word, basis)
    require(fequal(tx, rf(Nx, h2)) and fequal(ty, rf(Ny, V.power(h, 3))), "trace group-law attachment")
    require(V.power(Ny, 2) == V.add(V.add(V.power(Nx, 3), V.mul(A, V.mul(Nx, h4))), V.mul(B, h6)),
            "trace Weierstrass identity")
    require(not V.remainder(V.add(V.mul(m, Nx), Ny), h2), "regular chord congruence")
    numerators = [
        V.sub(V.sub(V.sub(V.sub(V.power(m, 4), V.scale(V.mul(V.power(m, 2), Nx), 6)),
                         V.scale(V.mul(m, Ny), 8)), V.scale(V.power(Nx, 2), 3)), V.scale(V.mul(A, h4), 4)),
        V.mul(V.sub(V.sub(V.scale(V.power(m, 3), 4), V.scale(V.mul(m, Nx), 12)), V.scale(Ny, 8)), h2),
        V.scale(V.mul(V.sub(V.power(m, 2), Nx), h4), 6),
        V.scale(V.mul(m, h6), 4), V.power(h, 8)]
    return h, numerators


def check_mestre(data, gate):
    A = V.poly(data["jacobian_model"]["a4_low_to_high"])
    B = V.poly(data["jacobian_model"]["a6_low_to_high"])
    D = V.scale(V.mul(V.mul(A, B), V.add(V.scale(V.power(B, 2), 9261), V.scale(V.power(A, 3), 400))), -5)
    delta = V.add(V.scale(V.power(A, 3), 4), V.scale(V.power(B, 2), 27))
    require(gate["status"] == "PASS" and len(D)-1 == gate["cover_degree"] == 44, "Mestre scope")
    p = gate["smooth_disjoint_prime"]
    d, e = M.reduction(D, p), M.reduction(delta, p)
    for f in [d, e]:
        derivative = V.trim([i*x % p for i, x in enumerate(f)][1:])
        require(len(M.ff_gcd(f, derivative, p)) == 1, "squarefree reduction")
    require(len(M.ff_gcd(d, e, p)) == 1, "branch meets singular fibre")
    for name, f in [("A", A), ("B", B)]:
        row = gate["fields"][name]
        ell = row["irreducible_prime"]
        require(M.irreducible(M.reduction(f, ell), ell), "root field irreducibility")
        p, r = row["local_place"]["prime"], row["local_place"]["root"]
        reduced = M.reduction(f, p)
        derivative = V.trim([i*x % p for i, x in enumerate(reduced)][1:])
        require(M.evaluate(reduced, r, p) == 0 and M.evaluate(derivative, r, p), "simple local root")
        require(p % 4 == 3 and (name != "A" or p % 3 == 2), "cyclotomic exclusion")
    require(gate["fixed_u2_genus"] == 21 and gate["all_rational_auxiliary_genus_lower_bound"] == 9, "genus conclusion")
    return True


def local_values(q, prime, depth):
    den = lcm(*(x.denominator for x in q))
    coeff = [int(x*den*den) for x in q]+[0]*(3-len(q))
    def valuation(x):
        out = 0
        while x and x % prime == 0:
            out += 1
            x //= prime
        return out
    common = min(valuation(x) for x in coeff if x)
    coeff = [x//prime**(2*(common//2)) for x in coeff]
    modulus = prime**depth
    squares = {x*x % modulus for x in range(modulus)}
    addresses = [(x, 1) for x in range(modulus)]+[(1, prime*x) for x in range(modulus//prime)]
    return {i for i, (x, z) in enumerate(addresses)
            if (coeff[2]*x*x+coeff[1]*x*z+coeff[0]*z*z) % modulus in squares}


def verify(directory):
    packet = json.loads((directory/"input.json").read_text())
    result = json.loads((directory/"result.json").read_text())
    pencils = json.loads((directory/"pencils.json").read_text())
    source = ROOT/packet["source"]
    require(sha256(source.read_bytes()).hexdigest() == packet["source_sha256"], "parent source binding")
    require(result["input_sha256"] == sha256((directory/"input.json").read_bytes()).hexdigest(), "input binding")
    data = json.loads(source.read_text())
    c = Q(result["common_chart"]["c"])
    A, B = [V.poly(result["common_chart"][k]) for k in ("A", "B")]
    for key, f in [("A", A), ("B", B)]:
        degree = 8 if key == "A" else 12
        old = V.poly(data["jacobian_model"]["a4_low_to_high" if key == "A" else "a6_low_to_high"])
        require(fequal(mobius(rf(old), c, degree), rf(f)), "common chart")
    basis = source_basis(data, c)
    require(len(pencils) == result["pencil_count"] == len(packet["words8"]) == 35, "pencil coverage")
    matrices = []
    for row, word in zip(pencils, packet["words8"]):
        h, numerators = check_frame(row, word, basis, A, B)
        mat = [[Q(v) for v in r] for r in row["branch_matrix"]]
        require(len(mat) == 5 and all(len(r) == 5 for r in mat), "matrix shape")
        for col, n in zip(zip(*mat), numerators):
            require(V.mul(V.trim(col), V.power(h, 6)) == n, "branch matrix identity")
        matrices.append(mat)
    remaining = set(itertools.combinations(range(35), 2))
    require(result["pencil_pairs"] == len(remaining) == 595, "pair coverage")
    for row in result["projective_screens"]:
        p = row["prime"]
        require(p >= 2 and all(p % i for i in range(2, isqrt(p)+1)), "prime witness")
        images = {}
        for i, mat in enumerate(matrices):
            red = V.reduce_matrix(mat, p)
            if red is not None:
                images[i] = V.projective_image(red, p)
        require(sorted(images) == row["usable"], "usable pencil roster")
        rejected = sorted((i, j) for i, j in remaining if i in images and j in images and images[i].isdisjoint(images[j]))
        require([list(v) for v in rejected] == row["rejected_pairs"], "modular exclusion")
        remaining.difference_update(rejected)
        require(len(remaining) == row["remaining_pairs"], "modular coverage")
    require(not remaining and not result["unresolved_projective_pairs"], "finite pencil conclusion")
    covers = result["quadratic_covers"]
    require(len(covers) == 15 and len(packet["words10"]) == 13, "quadratic bank coverage")
    valid, rational = [], 0
    for i, (row, word) in enumerate(zip(covers, packet["words10"])):
        h, numerators = check_frame(row, word, basis, A, B)
        q = V.poly(row["q"])
        require(V.mul(q, V.power(h, 6)) == numerators[0], "quadratic chord identity")
        dq = [j*x for j, x in enumerate(q)][1:]
        smooth = len(q) in (2, 3) and len(V.gcd_poly(q, dq)) == 1
        if row.get("status") == "DEGENERATE_RESIDUAL":
            require(not smooth, "wrongly discarded smooth cover")
            continue
        require(smooth, "singular cover")
        x0, x1, y0, y1 = [V.poly(row[k]) for k in ("x0", "x1", "y0", "y1")]
        require(V.add(V.power(y0, 2), V.mul(V.power(y1, 2), q)) ==
                V.add(V.add(V.add(V.power(x0, 3), V.scale(V.mul(x0, V.mul(V.power(x1, 2), q)), 3)), V.mul(A, x0)), B),
                "quadratic section constant identity")
        require(V.scale(V.mul(y0, y1), 2) ==
                V.add(V.add(V.scale(V.mul(V.power(x0, 2), x1), 3), V.mul(V.power(x1, 3), q)), V.mul(A, x1)),
                "quadratic section radical identity")
        if row["rational_conic"]:
            t, z, s = map(Q, row["conic_point"])
            require((t or z) and s*s == q[2]*t*t+q[1]*t*z+q[0]*z*z, "conic rational point")
            rational += 1
        valid.append(i)
    for i, index in [(13, 4339), (14, 36773)]:
        old = json.loads((source.parent/f"cover-{index}.json").read_text())["quadratic_basechange"]["branch_coefficients"]
        require(fequal(mobius(rf(V.poly(old)), c, 2), rf(V.poly(covers[i]["q"]))), "published cover comparison")
        valid.append(i)
    remaining = set(itertools.combinations(valid, 2))
    for i, j in list(remaining):
        if len(V.gcd_poly(V.poly(covers[i]["q"]), V.poly(covers[j]["q"]))) != 1:
            remaining.remove((i, j))
    for row in result["local_sieves"]:
        p, depth = row["prime"], row["depth"]
        images = [local_values(V.poly(r["q"]), p, depth) for r in covers]
        rejected = sorted((i, j) for i, j in remaining if images[i].isdisjoint(images[j]))
        require([list(v) for v in rejected] == row["rejected_pairs"], "local squareclass sieve")
        remaining.difference_update(rejected)
        require(len(remaining) == row["remaining_pairs"], "local pair coverage")
    require([list(v) for v in sorted(remaining)] == result["locally_surviving_disjoint_pairs"], "local survivor roster")
    check_mestre(data, json.loads((directory/"mestre.json").read_text()))
    return {"schema": "hill1152-independent-pilot-replay-v1", "status": "PASS",
            "pencils": 35, "excluded_pencil_pairs": 595, "rational_quadratic_covers": rational,
            "degenerate_residuals": 13-(len(valid)-2), "locally_surviving_pairs": len(remaining),
            "mestre_genus_lower_bound": 9,
            "boundary": "Independent exact trace/chord/matrix/conic/local-witness arithmetic. No universal cover exclusion, complete Mordell-Weil basis, saturation, rational point on surviving paired bases, or new record is claimed."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = verify(args.input_dir)
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps(result, sort_keys=True))
