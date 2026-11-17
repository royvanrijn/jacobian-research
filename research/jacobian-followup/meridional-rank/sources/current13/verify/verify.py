#!/usr/bin/env python3
"""Check the bridge-index certificates in data/certificates.json.

Uses only the Python standard library. For every knot it checks:

  (A) the PD code is well formed (KnotInfo convention: [a,b,c,d] lists the edges at a
      crossing counterclockwise from the incoming under-edge a; edges are numbered 1..2n
      along the orientation);
  (B) a lower bound, meridional rank >= 4: the 4x4 matrices over F_q attached to the edges
      are constant along arcs, satisfy the Wirtinger relation at every crossing, each one
      minus the identity has rank 1, and the hyperplanes they fix meet only in 0;
  (C) an upper bound, bridge index <= 4: starting from the four listed seed arcs, the
      coloring moves of Blair-Kjuchukova-Velazquez-Villanueva color every arc, so the
      Wirtinger number of the diagram is at most 4.

Optionally (--knotinfo) it also checks that each PD code is the one in the database_knotinfo
package and that KnotInfo lists the bridge index of the knot as [3,4].

Usage: python3 verify/verify.py [--knotinfo] [path/to/certificates.json]
"""
import json, os, sys

def rank(rows, q):
    A = [[x % q for x in r] for r in rows]
    r = 0
    for c in range(len(A[0])):
        p = next((i for i in range(r, len(A)) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        iv = pow(A[r][c], -1, q)
        A[r] = [x * iv % q for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % q for x, y in zip(A[i], A[r])]
        r += 1
    return r

def mul(X, Y, q):
    return [[sum(X[i][k] * Y[k][j] for k in range(4)) % q for j in range(4)] for i in range(4)]

def inverse(X, q):
    A = [[x % q for x in r] + [int(i == j) for j in range(4)] for i, r in enumerate(X)]
    for c in range(4):
        p = next(i for i in range(c, 4) if A[i][c])
        A[c], A[p] = A[p], A[c]
        iv = pow(A[c][c], -1, q)
        A[c] = [x * iv % q for x in A[c]]
        for i in range(4):
            if i != c and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % q for x, y in zip(A[i], A[c])]
    return [r[4:] for r in A]

def is_prime(q):
    return q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))

def arcs_of(pd):
    """Arcs as lists of edges, by walking the knot: an arc ends at the edge that runs under."""
    n2 = 2 * len(pd)
    ends_under = {c[0] for c in pd}
    e = pd[0][2]
    arcs, cur = [], []
    for _ in range(n2):
        cur.append(e)
        if e in ends_under:
            arcs.append(cur)
            cur = []
        e = e % n2 + 1
    if cur:
        arcs[0] = cur + arcs[0]
    return arcs

def check_pd(pd):
    n2 = 2 * len(pd)
    seen = sorted(e for c in pd for e in c)
    if seen != sorted(list(range(1, n2 + 1)) * 2):
        return 'each edge 1..2n must occur exactly twice'
    for a, b, c, d in pd:
        if c != a % n2 + 1:
            return f'under strand not consecutive at {[a, b, c, d]}'
        if not (d == b % n2 + 1 or b == d % n2 + 1):
            return f'over strand not consecutive at {[a, b, c, d]}'
    return None

def lower_bound(pd, q, mats):
    n2 = 2 * len(pd)
    I4 = [[int(i == j) for j in range(4)] for i in range(4)]
    M = {e: [[x % q for x in r] for r in mats[e - 1]] for e in range(1, n2 + 1)}
    arcs = arcs_of(pd)
    for arc in arcs:
        if any(M[e] != M[arc[0]] for e in arc):
            return 'matrix not constant along an arc'
    for a, b, c, d in pd:
        O = M[b]
        Oi = inverse(O, q)
        # Wirtinger relation; the convention is fixed once for all crossings
        lhs = mul(mul(Oi, M[a], q), O, q) if d == b % n2 + 1 else mul(mul(O, M[a], q), Oi, q)
        if lhs != M[c]:
            return f'Wirtinger relation fails at {[a, b, c, d]}'
    stack = []
    for arc in arcs:
        D = [[(M[arc[0]][i][j] - I4[i][j]) % q for j in range(4)] for i in range(4)]
        if rank(D, q) != 1:
            return 'an arc matrix is not a pseudo-reflection'
        stack += D
    if rank(stack, q) != 4:
        return 'the common fixed space is not zero'
    return None

def upper_bound(pd, seeds):
    arcs = arcs_of(pd)
    arc_of = {e: i for i, arc in enumerate(arcs) for e in arc}
    colored = {arc_of[e] for e in seeds}
    if len(colored) != 4:
        return 'seeds must lie on four distinct arcs'
    changed = True
    while changed:
        changed = False
        for a, b, c, d in pd:
            o, u, v = arc_of[b], arc_of[a], arc_of[c]
            if o in colored and (u in colored) != (v in colored):
                colored |= {u, v}
                changed = True
    if len(colored) != len(arcs):
        return 'seed coloring does not reach every arc'
    return None

def knotinfo_rows():
    import csv, ast
    import database_knotinfo
    f = os.path.join(os.path.dirname(database_knotinfo.__file__), 'csv_data', 'knotinfo_data_complete.csv')
    csv.field_size_limit(10 ** 9)
    r = csv.reader(open(f, newline=''), delimiter='|')
    hdr = next(r)
    I = {h: i for i, h in enumerate(hdr)}
    out = {}
    for x in r:
        if x[I['name']] == 'name' or not x[I['pd_notation']].startswith('[['):
            continue
        out[x[I['name']]] = (ast.literal_eval(x[I['pd_notation']]), x[I['bridge_index']])
    return out

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    path = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'certificates.json')
    certs = json.load(open(path))
    ki = knotinfo_rows() if '--knotinfo' in sys.argv else None
    failed = 0
    for cert in certs:
        name, pd, q = cert['name'], cert['pd'], cert['field']
        err = (None if is_prime(q) else 'field size must be prime') or check_pd(pd) \
            or lower_bound(pd, q, cert['matrices']) or upper_bound(pd, cert['wirtinger_seeds'])
        if ki is not None and err is None:
            if name not in ki:
                err = 'not in KnotInfo'
            elif ki[name][0] != pd:
                err = 'PD code differs from KnotInfo'
            elif ki[name][1] != '[3,4]':
                err = f'KnotInfo bridge index is {ki[name][1]}, not [3,4]'
        if err:
            failed += 1
            print(f'{name}: FAILED: {err}')
    print(f'{len(certs) - failed} of {len(certs)} certificates verified: bridge index = meridional rank = 4'
          + (' (PD codes and [3,4] entries checked against KnotInfo)' if ki is not None else ''))
    sys.exit(1 if failed else 0)

if __name__ == '__main__':
    main()
