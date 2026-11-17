"""Exact radial intersection and cyclic-chain regressions.

Universal geometry and all-length proofs are in PROOF.md. Finite samples
below support them; they do not enumerate possible Keller boundary graphs.
"""
from pathlib import Path
import hashlib
import json
import math
import platform
import sys
import time
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CAS = ROOT/'foundations/plane-jc/cas'
sys.path.insert(0, str(CAS))
from intrinsic_a2_boundary import symmetric_inertia

checks = []


def check(name, condition):
    assert condition, name
    checks.append(name)


def zero(name, expr):
    values = list(expr) if isinstance(expr, S.MatrixBase) else [expr]
    check(name, all(S.cancel(x) == 0 for x in values))


def chain(weights):
    size = len(weights)
    A = S.diag(*weights)
    for i in range(size-1):
        A[i, i+1] = A[i+1, i] = -1
    return A


def prefixes(weights):
    result = [1]
    previous = 0
    for b in weights:
        result.append(b*result[-1]-previous)
        previous = result[-2]
    return result


def integral(v):
    return all(x.is_Integer for x in v)


def main():
    start = time.monotonic()
    representatives = [[2], [3], [2, 2], [2, 2, 2], [2, 3, 2],
                       [3, 2, 4, 2], [2]*9, [2, 5], [4, 3, 2, 2, 3]]
    chain_receipts = []
    for weights in representatives:
        A = chain(weights)
        size = len(weights)
        K = prefixes(weights)
        delta = K[-1]
        tag = 'chain '+','.join(map(str, weights))
        check(tag+' determinant recurrence', delta == A.det())
        check(tag+' strict positive prefix increase', all(a < b for a, b in zip(K, K[1:])))
        check(tag+' negative exceptional inertia', symmetric_inertia(-A) == (0, size, 0))
        inverse = A.inv()
        check(tag+' inverse is strictly positive', all(x > 0 for x in inverse))
        for i in range(size):
            for j in range(size):
                low, high = min(i, j), max(i, j)
                tail = prefixes(weights[high+1:])[-1]
                zero(tag+f' Green inverse {i},{j}', inverse[i, j]-S.Rational(K[low]*tail, delta))
        minimum_degrees = []
        for j in range(size):
            unit = S.zeros(size, 1)
            unit[j] = 1
            h = delta//math.gcd(delta, K[j])
            check(tag+f' no sole degree-one horizontal at {j}', h >= 2 and not integral(inverse*unit))
            f = h*unit
            e = inverse*f
            check(tag+f' exact minimum sole-horizontal degree at {j}',
                  integral(e) and all(x > 0 for x in e)
                  and not integral(inverse*((h-1)*unit)))
            zero(tag+f' resolved intersection equation at {j}', A*e-f)
            zero(tag+f' local length quadratic identity at {j}', (e.T*f)[0]-(e.T*A*e)[0])
            minimum_degrees.append(h)
        samples = [S.ones(size, 1), S.Matrix([j % 3 for j in range(size)])]
        ends = S.zeros(size, 1)
        ends[0] += 1
        ends[-1] += 1
        samples.append(ends)
        for j, f in enumerate(samples):
            congruent = sum(K[i]*int(f[i]) for i in range(size)) % delta == 0
            check(tag+f' cyclic congruence versus direct inversion {j}', congruent == bool(integral(inverse*f)))
        chain_receipts.append(dict(weights=weights, determinant=delta,
                                   minimum_sole_horizontal_degrees=minimum_degrees))
    # Single-curve degree-one case, including self-intersection -1.
    for b in range(1, 8):
        e = S.Rational(1, b)
        check(f'single curve b={b} degree-one integrality', bool(e.is_Integer) == (b == 1))
    # A_1 model contains an e=1 radial divisor even though its cover ramifies.
    zero('A1 e=1 f=2 local block', S.Matrix([[-2]])*S.Matrix([1])+S.Matrix([2]))
    node_receipts = []
    a, b, t = S.symbols('a b t')
    for n in range(2, 13):
        tag = f'node n={n}'
        size = n-1
        Q = -chain([2]*size)
        e = S.Matrix([min(j, n-j) for j in range(1, n)])
        unresolved = []
        for i in range(n):
            u = a**(i+1)*b**i
            v = a**(n-i-1)*b**(n-i)
            w = a*b
            zero(tag+f' chart {i} surface equation', u*v-w**n)
            zero(tag+f' chart {i} volume order',
                 S.Matrix([u, v]).jacobian((a, b)).det()-n*a**(n-1)*b**(n-1))
            ua, ub, va, vb = i+1, i, n-i-1, n-i
            comparable = (ua <= va and ub <= vb) or (va <= ua and vb <= ub)
            if not comparable:
                unresolved.append(i)
            if i < n-1:
                aa, bb = 1/b, a*b*b
                next_u = aa**(i+2)*bb**(i+1)
                next_v = aa**(n-i-2)*bb**(n-i-1)
                zero(tag+f' chart {i} transition preserves u,v,w',
                     S.Matrix([next_u-u, next_v-v, aa*bb-w]))
        if n % 2 == 0:
            m = n//2
            check(tag+' radial ideal already principal on every chart', unresolved == [])
            f = S.zeros(size, 1)
            f[m-1] = 2
            zero(tag+' even radial labels', Q*e+f)
            # Middle chart ratio shows degree two, not just lattice support.
            i = m-1
            u = a**(i+1)*b**i
            v = a**(n-i-1)*b**(n-i)
            zero(tag+' horizontal residue map degree two', v/u-b*b)
            k = S.ones(size, 1)*(n-1)
        else:
            m = n//2
            check(tag+' unique middle unresolved radial basepoint', unresolved == [m])
            i = m
            u = a**(i+1)*b**i
            v = a**(n-i-1)*b**(n-i)
            zero(tag+' middle residual ideal generators a,b',
                 S.Matrix([u/(a**m*b**m)-a, v/(a**m*b**m)-b]))
            U = u.subs(b, a*t)
            V = v.subs(b, a*t)
            zero(tag+' new blowup chart has principal radial ideal', S.Matrix([U/a**n-t**m, V/U-t]))
            zero(tag+' new horizontal volume order',
                 S.Matrix([U, V]).jacobian((a, t)).det()-n*a**(2*n-1)*t**(n-1))
            expanded = S.zeros(size+1)
            expanded[:size, :size] = Q
            left, right = m-1, m
            expanded[left, left] -= 1
            expanded[right, right] -= 1
            expanded[left, right] = expanded[right, left] = 0
            expanded[left, size] = expanded[size, left] = 1
            expanded[right, size] = expanded[size, right] = 1
            expanded[size, size] = -1
            Q = expanded
            e = S.Matrix(list(e)+[n])
            f = S.zeros(size+1, 1)
            f[size] = 1
            k = S.Matrix([n-1]*size+[2*n-1])
            zero(tag+' complete odd radial labels after required blowup', Q*e+f)
            check(tag+' new horizontal is outside the all-b>=2 class', Q[size, size] == -1)
        check(tag+' resolved negative inertia', symmetric_inertia(Q) == (0, Q.rows, 0))
        check(tag+' pointwise local fibre length', (e.T*f)[0] == n)
        zero(tag+' local length intersection square', (e.T*f)[0]+(e.T*Q*e)[0])
        for j in range(Q.rows):
            if f[j] > 0:
                check(tag+f' horizontal canonical order {j}', k[j] == 2*e[j]-1)
            else:
                check(tag+f' vertical canonical order {j}', k[j] >= 2*e[j])
        node_receipts.append(dict(n=n, radial_basepoint_count_before_extra_blowup=len(unresolved),
                                 Q=[[int(x) for x in Q.row(j)] for j in range(Q.rows)],
                                 e=list(map(int, e)), f=list(map(int, f)), k=list(map(int, k)),
                                 local_length=int((e.T*f)[0]),
                                 local_control_is_not_a_Keller_map=True))
    manifest = json.loads((HERE/'sources/manifest.json').read_text())
    for item in manifest['sources']:
        path = HERE/'sources'/item['file']
        check('primary pin '+item['file'], hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'])
    result = dict(as_of='2026-11-13', outcome='PASS', check_count=len(checks), checks=checks,
                  runtime_seconds=round(time.monotonic()-start, 6),
                  checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  inherited_inertia_source_sha256=hashlib.sha256((CAS/'intrinsic_a2_boundary.py').read_bytes()).hexdigest(),
                  python=platform.python_version(), sympy=S.__version__,
                  chain_representatives=chain_receipts, node_representatives=node_receipts,
                  universal_claims='Written all-length congruence and both parity proofs, not finite enumeration',
                  global_geometry='Written exceptional negativity, completion splitting, radial degree and canonical-form arguments',
                  excluded_stratum='One radial horizontal of degree1 on a resolved chain with every weight at most -2',
                  unrestricted_pointed_Keller_gate='OPEN', formal_verification=False, external_review=False)
    (HERE/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
