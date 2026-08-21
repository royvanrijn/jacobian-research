from sage.all import *
from pathlib import Path
import argparse

ap = argparse.ArgumentParser(description="Verify one exact U-neighbor transition between positive frame Grams.")
ap.add_argument("--parent", type=Path, required=True)
ap.add_argument("--child", type=Path, required=True)
ap.add_argument("--q", type=int, required=True)
ap.add_argument("--a", type=int, required=True)
ap.add_argument("--b", type=int, required=True)
ap.add_argument("--v", required=True, help="comma-separated vector in the parent-frame basis")
args = ap.parse_args()


def read_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()]
                       for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


M = read_matrix(args.parent)
expected = read_matrix(args.child)
v = vector(ZZ, [ZZ(x) for x in args.v.split(",")])
assert M.is_square() and M.is_symmetric() and M.is_positive_definite()
frame_rank = M.nrows()
assert expected.nrows() == frame_rank and expected.is_symmetric()
assert expected.is_positive_definite() and expected.det() == M.det()
assert len(v) == frame_rank and v * M * v == 2 * args.q
assert args.a * args.b == args.q

U = matrix(ZZ, [[0, 1], [1, 0]])
NS = block_diagonal_matrix(U, -M)
f = vector(ZZ, [args.a, args.b] + list(v))
assert f * NS * f == 0
assert gcd([abs(ZZ(x)) for x in NS * f]) == 1

p = list(NS * f)
cur = ZZ(0)
g = [ZZ(0)] * (frame_rank + 2)
for i, pi in enumerate(p):
    if pi == 0:
        continue
    gg, s, t = xgcd(cur, ZZ(pi))
    g = [s * x for x in g]
    g[i] += t
    cur = gg
assert abs(cur) == 1
if cur == -1:
    g = [-x for x in g]
g = vector(ZZ, g)
assert f * NS * g == 1
gsq = ZZ(g * NS * g)
assert gsq % 2 == 0
g0 = g - (gsq // 2) * f
assert g0 * NS * g0 == 0 and f * NS * g0 == 1

K = matrix(ZZ, [list(f * NS), list(g0 * NS)]).right_kernel_matrix()
actual = -(K * NS * K.transpose())
assert actual == expected
assert actual.is_positive_definite()

coeffs = []
for i in range(frame_rank):
    for j in range(i, frame_rank):
        coeffs.append(actual[i, i] // 2 if i == j else actual[i, j])
Q = QuadraticForm(ZZ, frame_rank, coeffs)
shells = Q.short_vector_list_up_to_length(2, True)
half = [vector(ZZ, r) for r in shells[1]]
signed = half + [-r for r in half]
root_module = matrix(ZZ, [list(r) for r in signed]).row_module()
root_basis = root_module.basis_matrix()
root_gram = root_basis * actual * root_basis.transpose()

print(
    "FIBNEIGHBOR|status=PASS"
    f"|q={args.q}|ab={args.a},{args.b}|root_rank={root_basis.rank()}"
    f"|roots={len(signed)}|rootdet={abs(root_gram.det())}|MW={frame_rank-root_basis.rank()}"
)
