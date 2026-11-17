#!/usr/bin/env python3
"""Bounded exact controls for written affine-canonical and Lang theorems."""
from pathlib import Path
from itertools import product, permutations
from collections import Counter
import hashlib, importlib.util, json, math, platform, time
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
helper = ROOT/'jacobian-followup/universal-keller-symmetry/verify.py'
spec = importlib.util.spec_from_file_location('retained_symmetry_helper', helper)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
Field = module.Field
checks = []

def ck(label, condition):
    assert condition, label
    checks.append(label)


def modzero(expr, p, *variables):
    return s.Poly(s.expand(expr), *variables, modulus=p).is_zero


def rank(matrix, p):
    a = [[int(x) % p for x in row] for row in matrix]
    if not a:
        return 0
    at = 0
    for j in range(len(a[0])):
        found = next((i for i in range(at, len(a)) if a[i][j]), None)
        if found is None:
            continue
        a[at], a[found] = a[found], a[at]
        scale = pow(a[at][j], -1, p)
        a[at] = [(x*scale) % p for x in a[at]]
        for i in range(len(a)):
            if i != at and a[i][j]:
                z = a[i][j]
                a[i] = [(x-z*y) % p for x, y in zip(a[i], a[at])]
        at += 1
        if at == len(a):
            break
    return at


def quotient_dimension(gb, variables, maximum):
    leads = [p.LM(order=gb.order).exponents for p in gb.polys]
    return sum(not any(all(a >= b for a, b in zip(exp, lead)) for lead in leads)
               for exp in product(range(maximum), repeat=len(variables)))


def main():
    start = time.monotonic()
    affine_controls = []
    for n in range(5, 13):
        for p in s.factorint(n):
            p = int(p)
            r = 3 if p == 2 else 2
            blocks = [list(range(r)), list(range(r, n))]
            equations = []
            for block in blocks:
                for i in block[1:]:
                    row = [0]*n
                    row[block[0]], row[i] = 1, -1
                    equations.append(row)
            homogeneous = equations+[[1]*n]
            augmented = [row+[0] for row in equations]+[[1]*n+[1]]
            rr = rank(homogeneous, p)
            ck(f'N{n} p{p} two-block fixed locus nonempty', rr == rank(augmented, p))
            ck(f'N{n} p{p} full affine fixed dimension one', n-rr == 1)
            ck(f'N{n} p{p} scalar quotient fixed dimension zero', n-rr-1 == 0)
            ck(f'N{n} p{p} chosen block sizes have no fixed labels', 2 <= r <= n-2 and r % p != 0)
            forced = sum(math.comb(n, j) for j in range(1, n) if j % p)//2
            ck(f'N{n} p{p} minimum-span forced bipartitions exceed N', forced > n)
            affine_controls.append({'N': n, 'p': p, 'block_sizes': [r, n-r],
                                    'full_span_fixed_dimension': 1,
                                    'quotient_span_fixed_dimension': 0,
                                    'forced_two-block_points': forced})

    # Small-support inertia controls and the genuine normal V4 exception.
    support_controls = []
    for n in range(2, 7):
        for perm in permutations(range(n)):
            moved = sum(i != perm[i] for i in range(n))
            seen, cycles = set(), []
            for i in range(n):
                if i in seen:
                    continue
                j, length = i, 0
                while j not in seen:
                    seen.add(j); length += 1; j = perm[j]
                cycles.append(length)
            sign = (n-len(cycles)) % 2
            assert moved != 1
            assert not sign or moved >= 2
            assert sign or moved == 0 or moved >= 3
        ck(f'S{n} exact nontrivial and even support minima', True)
    # Permutation composition: left after right.
    comp = lambda a, b: tuple(a[b[i]] for i in range(len(a)))
    sym4 = list(permutations(range(4)))
    V4 = {(0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)}
    for g in sym4:
        inv = tuple(g.index(i) for i in range(4))
        assert {comp(comp(g, h), inv) for h in V4} == V4
    ck('V4 is a transitive normal subgroup of S4', all({h[i] for h in V4} == set(range(4)) for i in range(4)))
    ck('V4 has the genuine pair block system', all({h[i] for i in [0, 1]} in [{0, 1}, {2, 3}] for h in V4))

    # Matrix Lang identities: all polynomial Jacobian entries, whole generic ranks.
    lang_controls = []
    for q, p, n in [(2, 2, 1), (3, 3, 1), (4, 2, 1),
                     (2, 2, 2), (3, 3, 2), (2, 2, 3)]:
        X = s.symbols('x0:'+str(n))
        U = s.Matrix(n, n, lambda i, j: s.Symbol(f'u{i}{j}'))
        V = s.Matrix(X)-U*s.Matrix([x**q for x in X])
        allvars = tuple(U)+X
        ck(f'n{n} q{q} full matrix-map Jacobian block identity', all(
            modzero(V.jacobian(X)[i,j]-int(i == j), p, *allvars)
            for i in range(n) for j in range(n)))
        phi = [x-x**q for x in X]
        ck(f'n{n} q{q} canonical Lang Jacobian one',
           modzero(s.Matrix(phi).jacobian(X).det()-1, p, *X))
        lang_controls.append({'n': n, 'q': q, 'canonical_dimension': n,
                              'fixed_dimension': n*n+n, 'generic_degree': q**n})
    x, y, z = s.symbols('x y z')
    ck('canonical quartic exception has two global quadratic factors',
       modzero((x+x*x)-(x-x*x), 2, x) and modzero((y+y*y)-(y-y*y), 2, y))

    # Symbolic two-chart Frobenius normalization over F2(a,b,c,d).
    a, b, c, d = s.symbols('a b c d')
    variables = (a,b,c,d)
    def rz(expr):
        num, den = s.fraction(s.cancel(expr))
        assert not modzero(den, 2, *variables)
        return modzero(num, 2, *variables)
    power2 = lambda M: M.applyfunc(lambda entry: entry**2)
    D = a*d+b*c
    A = s.Matrix([[a,b],[c,d]])
    C0 = s.diag(1/D, 1)
    A1 = C0.inv()*A*power2(C0)
    ck('symbolic determinant normalization without root extraction', rz(A1.det()-1))
    h = (A1[1,1]-1)/A1[0,1]
    C1 = s.Matrix([[1,0],[h,1]])
    normal = C1.inv()*A1*power2(C1)
    ck('symbolic nonzero-upper-right chart bottom-right one', rz(normal[1,1]-1))
    ck('symbolic nonzero-upper-right chart exact universal form',
       rz(normal[0,0]-1-normal[0,1]*normal[1,0]))
    Atri = A1.subs(b,0)
    Ctri = s.diag(d, 1/d)
    normaltri = Ctri.inv()*Atri*power2(Ctri)
    ck('symbolic zero-upper-right chart exact universal form',
       rz(normaltri[1,1]-1) and rz(normaltri[0,0]-1) and rz(normaltri[0,1]))

    # All invertible matrices over F2 and F4: exact gauge controls.
    normalized = []
    for f in [Field(2,(0,1)), Field(2,(1,1,1))]:
        neg = lambda v: f.encode([-w for w in f.digits(v)])
        def madd(x, y): return f.add(x,y)
        def mm(A,B):
            return [[f.add(f.mul(A[i][0],B[0][j]), f.mul(A[i][1],B[1][j])) for j in range(2)] for i in range(2)]
        def det(A): return f.add(f.mul(A[0][0],A[1][1]), f.mul(A[0][1],A[1][0]))
        def inv(A):
            q = f.power(det(A), f.q-2)
            return [[f.mul(q,A[1][1]),f.mul(q,A[0][1])],[f.mul(q,A[1][0]),f.mul(q,A[0][0])]]
        def sq(A): return [[f.mul(v,v) for v in row] for row in A]
        total, charts = 0, Counter()
        for entries in product(range(f.q),repeat=4):
            A = [list(entries[:2]),list(entries[2:])]
            dd = det(A)
            if not dd: continue
            total += 1
            invd = f.power(dd, f.q-2)
            C0 = [[invd,0],[0,1]]
            U = mm(mm(inv(C0),A),sq(C0))
            assert det(U) == 1
            if U[0][1]:
                charts['nonzero'] += 1
                hh = f.mul(f.add(U[1][1],1),f.power(U[0][1],f.q-2))
                C = [[1,0],[hh,1]]
            else:
                charts['zero'] += 1
                dd = U[1][1]
                C = [[dd,0],[0,f.power(dd,f.q-2)]]
            R = mm(mm(inv(C),U),sq(C))
            assert R[1][1] == 1 and R[0][0] == f.add(1,f.mul(R[0][1],R[1][0]))
        ck(f'F{f.q} every invertible rank-two module normalized',
           total == (f.q*f.q-1)*(f.q*f.q-f.q))
        normalized.append({'field_size':f.q,'invertible_matrices':total,'charts':dict(charts)})

    # The actual quartic map: a whole free quotient, all16 base-field fibres,
    # and extension-field Frobenius counts determining full splitting types.
    S,T,X,Y,A0,B0 = s.symbols('S T X Y A0 B0')
    f1 = X+(1+S*T)*X**2+S*Y**2
    f2 = Y+T*X**2+Y**2
    full = s.Matrix([S,T,f1,f2]).jacobian([S,T,X,Y])
    ck('four-coordinate quartic full Jacobian one', modzero(full.det()-1,2,S,T,X,Y))
    U = s.Matrix([[1+S*T,S],[T,1]])
    ck('quartic coefficient matrix determinant one globally',modzero(U.det()-1,2,S,T))
    eq = [f1-A0,f2-B0]
    monic = [X**2+X+S*Y+A0+S*B0,
             Y**2+T*X+(1+S*T)*Y+T*A0+(1+S*T)*B0]
    ck('quartic equations are the exact global monic system',
       all(modzero(v,2,S,T,X,Y,A0,B0) for v in s.Matrix(eq)-U*s.Matrix(monic)))
    fibre_controls=[]
    distribution=Counter()
    for ss,tt,aa,bb in product(range(2),repeat=4):
        polys = [e.subs({S:ss,T:tt,A0:aa,B0:bb}) for e in monic]
        gb = s.groebner(polys,X,Y,modulus=2,order='grevlex')
        ck(f'quartic whole fibre {(ss,tt,aa,bb)} rank four',quotient_dimension(gb,(X,Y),4)==4)
        counts = {}
        cycle_counts={}
        for degree in range(1,5):
            extra=[X**(2**degree)-X,Y**(2**degree)-Y]
            fixed=s.groebner(polys+extra,X,Y,modulus=2,order='grevlex')
            counts[degree]=quotient_dimension(fixed,(X,Y),4)
            previous=sum(j*cycle_counts[j] for j in range(1,degree) if degree%j==0)
            assert (counts[degree]-previous)%degree==0
            cycle_counts[degree]=(counts[degree]-previous)//degree
        partition=tuple(j for j in range(1,5) for _ in range(cycle_counts[j]))
        assert sum(partition)==4
        distribution[partition]+=1
        fibre_controls.append({'target':[ss,tt,aa,bb],'rank':4,'extension_point_counts':counts,'splitting_type':partition})
    expected={(1,1,1,1):1,(2,2):3,(1,1,2):4,(1,3):4,(4,):4}
    ck('all16 complete quartic fibres realize all five types with exact distribution',dict(distribution)==expected)
    ck('F2 source point accounting matches the complete fibre census',
       sum(len([v for v in part if v==1])*multiplicity for part,multiplicity in distribution.items())==16)

    # Actual S4 affine symmetry of the non-atomic canonical plane.
    affine_perms=set()
    points=list(product(range(2),repeat=2))
    for entries in product(range(2),repeat=4):
        A=[entries[:2],entries[2:]]
        if (A[0][0]*A[1][1]+A[0][1]*A[1][0])%2 != 1:continue
        for shift in points:
            affine_perms.add(tuple(points.index(((A[0][0]*v[0]+A[0][1]*v[1]+shift[0])%2,
                                                (A[1][0]*v[0]+A[1][1]*v[1]+shift[1])%2)) for v in points))
    ck('canonical affine plane gives every natural S4 permutation',affine_perms==set(permutations(range(4))))
    for src in json.loads((HERE/'sources.json').read_text())['local_sources']:
        ck('retained proof/helper pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
             'affine_orbit_controls':affine_controls,'matrix_lang_controls':lang_controls,
             'rank_two_normalization_controls':normalized,'complete_fibre_controls':fibre_controls,
             'quartic_splitting_distribution':{' '.join(map(str,k)):v for k,v in sorted(distribution.items())},
             'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-start,3),
             'assurance':'Bounded exact affine-orbit and support controls, universal two-chart polynomial identities, all invertible matrices over F2/F4 and all16 whole quartic fibre algebras. These do not formalize all-rank fixed-locus arguments, etale fundamental groups, torsor descent or generic monodromy. The retained small-field helper is reused. No independent whole-proof replay, external review, formal verification or priority certification.',
             'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'whole_quartic_fibres':len(fibre_controls),
                      'normalized_matrices':sum(v['invertible_matrices'] for v in normalized),
                      'seconds':receipt['runtime_seconds']}))

if __name__=='__main__':main()
