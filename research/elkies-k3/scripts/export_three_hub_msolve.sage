from sage.all import *
from pathlib import Path
from collections import defaultdict, deque
import argparse, random

ap = argparse.ArgumentParser(description="Export the reduced Elkies rank-17 three-hub reconstruction system to msolve.")
ap.add_argument("--p", type=int, required=True, help="good prime, e.g. 101,103,107 (avoid 2,3,79)")
ap.add_argument("--slices", type=int, default=5, help="number of random affine linear slices; start with 4,5,6")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", required=True)
ap.add_argument("--slice-width", type=int, default=6, help="variables per random affine slice")
args = ap.parse_args()

p0 = args.p
if not is_prime(p0) or p0 in (2,3,79):
    raise SystemExit("--p must be a good prime (prime and not 2,3,79)")
F = GF(p0)
random.seed(args.seed)
set_random_seed(args.seed)

# Same connected -2 graph used by build_three_hub_system.sage.
edges = [
    (5,2),(5,3),(5,4),(6,5),(8,5),
    (9,3),(9,7),(9,8),(10,5),(10,9),
    (11,5),(11,9),(13,9),(14,3),(14,6),
    (14,7),(14,12),(14,13),(15,5),(15,9),
    (16,14),
]
vertices = sorted(set(sum(([a,b] for a,b in edges), [])))
assert vertices == list(range(2,17))

adj = defaultdict(list)
for a,b in edges:
    adj[a].append(b); adj[b].append(a)
root = 5
parent = {root: None}
tree = []
q = deque([root])
while q:
    a = q.popleft()
    for b in adj[a]:
        if b in parent: continue
        parent[b] = a; tree.append((a,b)); q.append(b)
assert len(tree) == 14
Tset = {tuple(sorted(e)) for e in tree}
cycle_edges = [e for e in edges if tuple(sorted(e)) not in Tset]

# Eliminate A(t),B(t) entirely. Unknowns are x_i(t), one root y_5(t), and
# the degree-2 slopes attached to the 21 -2 edges.
names = []
for i in vertices:
    names += [f"x{i}_{k}" for k in range(5)]
names += [f"y{root}_{k}" for k in range(7)]
for a,b in edges:
    lo,hi = sorted((a,b))
    names += [f"m{lo}_{hi}_{k}" for k in range(3)]

R = PolynomialRing(F, names, order='degrevlex')
g = R.gens_dict()
Rt = PolynomialRing(R, 't'); t = Rt.gen()

def poly(prefix,d):
    return sum(g[f"{prefix}_{k}"]*t**k for k in range(d+1))

def slope(a,b):
    lo,hi = sorted((a,b))
    return poly(f"m{lo}_{hi}",2)

def coeffs(P,d):
    return [R(P[k]) for k in range(d+1)]

x = {i: poly(f"x{i}",4) for i in vertices}
y = {root: poly(f"y{root}",6)}
children = defaultdict(list)
for a,b in tree:
    if parent[b] == a: children[a].append(b)
    elif parent[a] == b: children[b].append(a)
    else: raise RuntimeError("tree orientation mismatch")
qq = deque([root])
while qq:
    a = qq.popleft()
    for b in children[a]:
        y[b] = y[a] + slope(a,b)*(x[b]-x[a])
        qq.append(b)
assert set(y) == set(vertices)

# Pick one tree edge as the definition of A(t):
#   A = m(P+Q)*(yP+yQ) - (xP^2+xPxQ+xQ^2).
anchor = tree[0]
aa,bb = anchor
ma = slope(aa,bb)
A = ma*(y[aa]+y[bb]) - (x[aa]**2 + x[aa]*x[bb] + x[bb]**2)
assert A.degree() <= 8

E = []
# Every other tree edge must recover the same A(t).
for a,b in tree[1:]:
    m = slope(a,b)
    Ae = m*(y[a]+y[b]) - (x[a]**2 + x[a]*x[b] + x[b]**2)
    E += coeffs(Ae-A,8)

# Cycle edges impose both line compatibility and the same curve identity.
for a,b in cycle_edges:
    m = slope(a,b)
    E += coeffs(y[b]-y[a]-m*(x[b]-x[a]),6)
    Ae = m*(y[a]+y[b]) - (x[a]**2 + x[a]*x[b] + x[b]**2)
    E += coeffs(Ae-A,8)

# Remove literal zero equations after simplification.
E = [e for e in E if e != 0]
base_eq_count = len(E)

# Generic affine slices. Expected geometric dimension is plausibly 5:
# 3 from PGL2(base), 1 Weierstrass scaling, 1 K3 moduli parameter.
# Testing slices=4/5/6 across several primes is itself a useful dimension probe.
vars_ = list(R.gens())
for s in range(args.slices):
    chosen = random.sample(vars_, min(args.slice_width, len(vars_)))
    form = R(F.random_element())
    for v in chosen:
        c = F.random_element()
        while c == 0: c = F.random_element()
        form += R(c)*v
    E.append(form)

# msolve requires expanded polynomials with each monomial combined once.
def msolve_str(f):
    f = R(f)
    terms = []
    # Sage's str() is already canonical/combined; replace ** by ^ for msolve.
    return str(f).replace('**','^')

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w') as fh:
    fh.write(','.join(names) + '\n')
    fh.write(str(p0) + '\n')
    for i,e in enumerate(E):
        fh.write(msolve_str(e))
        fh.write(',\n' if i+1 < len(E) else '\n')

print(f"THREEHUBMSOLVE|p={p0}|seed={args.seed}|slices={args.slices}|vars={len(names)}|base_eqs={base_eq_count}|eqs={len(E)}|anchor={anchor}|out={out}")
print("THREEHUBMSOLVE|next=msolve -t <threads> -v 2 -g 1 -f %s -o %s.dim" % (out, out))
