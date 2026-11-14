# Verificateur : reductions preliminaires du cone binaire (HC4 quintique).
# Lancer : python3 verify_hc4_prelim.py   (aucun argument)
import time, itertools
T0 = time.time()
from sage.all__sagemath_singular import *
import flint
import sympy as sp

def say(*a):
    print(*a, flush=True)

NP = [0, 0]
def check(label, ok):
    NP[0 if ok else 1] += 1
    say("  [OK]   " if ok else "  [ECHEC]", label, " (%.0f s)" % (time.time() - T0))

def mon(n, d):
    return [e for e in itertools.product(range(d + 1), repeat=n) if sum(e) == d]

def det4(M):
    tot = None
    for p in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    sg = -sg
        t = M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]]
        tot = (t if sg > 0 else -t) if tot is None else (tot + t if sg > 0 else tot - t)
    return tot

say("=" * 72)
say("Verificateur : reductions preliminaires du cone binaire")
say("=" * 72)

# --- P1 : le cas B2 est exclu ---
say("\n[P1] Cas B2 : f4 contient c (x1 x3 + x2 s)**2")
x1, x2, x3, s, c = sp.symbols('x1 x2 x3 s c')
p = sp.symbols('p0:6')
phi = sum(p[i] * x1 ** i * x2 ** (5 - i) for i in range(6))
top = phi + c * (x1 * x3 + x2 * s) ** 2
Vv = [x1, x2, x3, s]
H = sp.Matrix(4, 4, lambda i, j: sp.diff(top, Vv[i], Vv[j]))
q = x1 * x3 + x2 * s
check("P1 det Hess(phi + c q**2) = -16 c**3 q**2 (10 phi - 3 c q**2)",
      sp.expand(H.det() + 16 * c ** 3 * q ** 2 * (10 * phi - 3 * c * q ** 2)) == 0)

# --- ansatz B1 : f = phi + A4(X,x3) + s B3(X) + termes de degre 2 et 3 + e s**3 ---
M5 = mon(2, 5); A4m = mon(3, 4); B3m = mon(2, 3)
low = [ex for d in (2, 3) for ex in mon(4, d) if ex != (0, 0, 0, 3)]
names = ['x1', 'x2', 'x3', 's'] + ['p%d' % i for i in range(6)] + ['a%d' % i for i in range(15)] + \
        ['b%d' % i for i in range(4)] + ['l%d' % i for i in range(len(low))] + ['e']
ctx = flint.fmpq_mpoly_ctx.get(names, 'degrevlex')
g = ctx.gens()
X1, X2, X3, S = g[:4]
FV = dict(zip(names, g))
Pc = g[4:10]; Ac = g[10:25]; Bc = g[25:29]; Lc = g[29:-1]; Ec = g[-1]
Z = ctx.from_dict({})
phiF = sum((cc * X1 ** k[0] * X2 ** k[1] for cc, k in zip(Pc, M5)), Z)
A4F = sum((cc * X1 ** k[0] * X2 ** k[1] * X3 ** k[2] for cc, k in zip(Ac, A4m)), Z)
B3F = sum((cc * X1 ** k[0] * X2 ** k[1] for cc, k in zip(Bc, B3m)), Z)
def lowF(skip):
    return sum((cc * X1 ** k[0] * X2 ** k[1] * X3 ** k[2] * S ** k[3] for cc, k in zip(Lc, low) if str(cc) not in skip), Z)
def hessF(f):
    return [[f.derivative(i).derivative(j) for j in range(4)] for i in range(4)]
def wcomp(D, w):
    return ctx.from_dict({ex: co for ex, co in D.to_dict().items() if ex[0] + ex[1] + ex[2] + 2 * ex[3] == w})
L = dict((m, Lc[i]) for i, m in enumerate(low))
l10 = L[(0, 0, 1, 2)]; l11 = L[(0, 0, 2, 1)]

say("\n[P2] Composante de poids 10 (poids 1,1,1,2)")
f = phiF + A4F + S * B3F + lowF([]) + Ec * S ** 3
D = det4(hessF(f))
W10 = wcomp(D, 10)
K = [[(phiF + S * B3F).derivative(i).derivative(j) for j in range(2)] for i in range(2)]
detK = K[0][0] * K[1][1] - K[0][1] * K[1][0]
br = 6 * Ec * S * A4F.derivative(2).derivative(2) + 12 * Ec * S ** 2 * l11 - 4 * S ** 2 * l10 ** 2
check("P2 W10 = det Hess_X(phi + s B3) * (6 e s d33 A4 + 12 e s**2 l11 - 4 s**2 l10**2)", W10 == detK * br)
say("     (det Hess_X(phi + s B3) en s = 0 vaut det Hess phi, non nul : le crochet est nul,")
say("      d'ou e d33 A4 = 0 et 3 e l11 = l10**2 ; si e = 0 alors l10 = 0)")

say("\n[P3] e = 0, L = x1 : A4 lineaire en x3 (poids 9, phi non degeneree)")
f0 = phiF + A4F + S * B3F + lowF(['l10', 'l13', 'l0', 'l19']) + S ** 2 * X1
W9 = wcomp(det4(hessF(f0)), 9)
vn = names[4:] + ['T', 'U']
RR = PolynomialRing(QQ, vn, order='degrevlex')
RG = RR.gens(); RV = dict(zip(vn, RG))
coeffs = {}
for ex, co in W9.to_dict().items():
    coeffs.setdefault(ex[:4], {})[ex[4:]] = QQ(int(co.p)) / QQ(int(co.q))
eqs9 = [sum(cv * prod(RG[i] ** int(k) for i, k in enumerate(e) if k) for e, cv in dct.items()) for dct in coeffs.values()]
Sx = PolynomialRing(RR, 'Y1,Y2')
Y1, Y2 = Sx.gens()
phiS = sum(RV['p%d' % i] * Y1 ** k[0] * Y2 ** k[1] for i, k in enumerate(M5))
dHS = phiS.derivative(Y1, 2) * phiS.derivative(Y2, 2) - phiS.derivative(Y1).derivative(Y2) ** 2
gks = [RR(cv) for cv in dHS.coefficients()]
for tname in ['a0', 'a1', 'a2', 'a5', 'a6', 'a9']:
    ok = all(list(RR.ideal(eqs9 + [RV['T'] * gk - 1, RV['U'] * RV[tname] - 1]).groebner_basis()) == [1] for gk in gks)
    check("P3 %s = 0 (coefficient de A4 de degre >= 2 en x3)" % tname, ok)

say("\n[P4] e = 0, L = 0 : f = P + s R avec A4 complet")
NM = []
for pref, k in [('p', 6), ('c', 4), ('d', 5), ('q', 10), ('r', 6), ('b', 4), ('m', 6), ('n', 3), ('h', 6)]:
    NM += [pref + str(i) for i in range(k)]
ctx2 = flint.fmpq_mpoly_ctx.get(['x1', 'x2', 'x3', 's'] + NM, 'lex')
G2 = ctx2.gens(); V2 = dict(zip(['x1', 'x2', 'x3', 's'] + NM, G2))
y1, y2, y3, ys = G2[:4]
bm = lambda e: y1 ** e[0] * y2 ** e[1]
tm = lambda e: y1 ** e[0] * y2 ** e[1] * y3 ** e[2]
HM = [(2, 0, 2), (1, 1, 2), (0, 2, 2), (1, 0, 3), (0, 1, 3), (0, 0, 4)]
PP = sum(V2['p%d' % i] * bm(e) for i, e in enumerate(mon(2, 5))) + y3 * sum(V2['c%d' % i] * bm(e) for i, e in enumerate(mon(2, 3)))
PP += sum(V2['d%d' % i] * bm(e) for i, e in enumerate(mon(2, 4))) + sum(V2['q%d' % i] * tm(e) for i, e in enumerate(mon(3, 3)))
PP += sum(V2['r%d' % i] * tm(e) for i, e in enumerate(mon(3, 2))) + sum(V2['h%d' % i] * tm(e) for i, e in enumerate(HM))
RRr = sum(V2['b%d' % i] * bm(e) for i, e in enumerate(mon(2, 3))) + sum(V2['m%d' % i] * tm(e) for i, e in enumerate(mon(3, 2)))
RRr += sum(V2['n%d' % i] * tm(e) for i, e in enumerate(mon(3, 1)))
D2 = det4([[(PP + ys * RRr).derivative(i).derivative(j) for j in range(4)] for i in range(4)])
R2 = PolynomialRing(QQ, NM + ['TT'], order='degrevlex')
G2r = R2.gens(); V2r = dict(zip(NM + ['TT'], G2r))
cmap = {}
for ex, co in D2.to_dict().items():
    cmap.setdefault(ex[:4], {})[ex[4:]] = QQ(int(co.p)) / QQ(int(co.q))
POL = {k: sum(cv * prod(G2r[i] ** int(kk) for i, kk in enumerate(e) if kk) for e, cv in d.items()) for k, d in cmap.items()}
def sysl(sub):
    Ssub = {V2r[k]: R2(QQ(v)) for k, v in sub.items()}
    eqs = []
    cst = R2(0)
    for k, pl in POL.items():
        e = pl.subs(Ssub)
        if k == (0, 0, 0, 0):
            cst = e
        elif e != 0:
            eqs.append(e)
    return eqs, cst
def incoh(sub):
    eqs, cst = sysl(sub)
    return list(R2.ideal(eqs + [V2r['TT'] * cst - 1]).groebner_basis()) == [1]
B21 = {'b0': 0, 'b1': 0, 'b2': 1, 'b3': 0}
C3 = {'b0': 0, 'b1': 0, 'b2': 0, 'b3': 1}
def dd(*ds):
    out = {}
    for d in ds:
        out.update(d)
    return out
def implied_s2(sub, target):
    Ssub = {V2r[k]: R2(QQ(v)) for k, v in sub.items()}
    eqs = [pl.subs(Ssub) for k, pl in POL.items() if k[3] == 2]
    eqs = [e for e in eqs if e != 0]
    U2 = PolynomialRing(QQ, NM + ['TT', 'UU'], order='degrevlex')
    t = U2(R2(target).subs(Ssub))
    if t == 0:
        return True
    return list(U2.ideal([U2(e) for e in eqs] + [U2.gens()[-1] * t - 1]).groebner_basis()) == [1]
check("P4 B3 = x1**3 : m1**2 = 4 m0 m2 (couche s**2, ne depend que de R)",
      implied_s2(C3, V2r['m1'] ** 2 - 4 * V2r['m0'] * V2r['m2']))
check("P4 B3 = x1**3, m0 = 0 : m1 = 0", implied_s2(dd(C3, {'m0': 0}), V2r['m1']))
check("P4 B3 = x1**3, m0 = 1 (m1 = m3 = m5 = n0 = 0 normalises) : m2 = 0",
      implied_s2(dd(C3, {'m0': 1, 'm1': 0, 'm3': 0, 'm5': 0, 'n0': 0}), V2r['m2']))
check("P4 B3 = x1**2 x2, R = ... + x3 : incoherent", incoh(dd(B21, {'m0': 0, 'm1': 0, 'm3': 0, 'n0': 1})))
LX = {'m1': 0, 'm2': 0, 'm3': 0, 'm5': 0, 'n0': 0}
check("P4 B3 = x1**3, Q = x3**2 + x1 x2 : incoherent", incoh(dd(C3, LX, {'m0': 1, 'm4': 1})))
check("P4 B3 = x1**3, Q = x3**2 : incoherent", incoh(dd(C3, LX, {'m0': 1, 'm4': 0})))
LY = {'m0': 0, 'm1': 0, 'm4': 0, 'm5': 0, 'n0': 0, 'n2': 0}
check("P4 B3 = x1**3, Q = x1 x3 + x2**2 : incoherent", incoh(dd(C3, LY, {'m3': 1, 'm2': 1, 'n1': 0})))
check("P4 B3 = x1**3, Q = x1 x3, R lin = x2 : incoherent", incoh(dd(C3, LY, {'m3': 1, 'm2': 0, 'n1': 1})))
check("P4 B3 = x1**3, Q = x1 x3 : incoherent", incoh(dd(C3, LY, {'m3': 1, 'm2': 0, 'n1': 0})))
eqs, cst = sysl(dd(C3, {'m0': 0, 'm1': 0, 'm3': 0, 'n0': 1, 'n1': 0, 'n2': 0}))
Pq = R2.ideal(eqs + [V2r['TT'] * cst - 1]).minimal_associated_primes()
check("P4 B3 = x1**3, R = x1**3 + Q(X) + x3 : chaque composante a p0 = ... = p4 = 0 (phi degeneree)",
      all(Q.reduce(V2r['p%d' % i]) == 0 for Q in Pq for i in range(5)))

say("\n[P5] Cas e != 0, R = x1**3 + x3 avec P = a(x1,x3) + x2 b(x1)")
aa = sp.Function('a')(x1, x3); bb = sp.Function('b')(x1)
f5 = aa + x2 * bb + s * (x1 ** 3 + x3) + s ** 3
H5 = sp.Matrix(4, 4, lambda i, j: sp.diff(f5, Vv[i], Vv[j]))
check("P5 det Hess f = -(6 s a_33 - 1) b'**2",
      sp.expand(H5.det() + (6 * s * sp.diff(aa, x3, 2) - 1) * sp.diff(bb, x1) ** 2) == 0)
say("     (donc a_33 = 0 et b' constant : le coefficient p4 de x1**4 x2 est nul, phi = p5 x1**5)")

say("\n" + "=" * 72)
say("BILAN : %d verifications passees, %d echecs  (%.0f s)" % (NP[0], NP[1], time.time() - T0))
say("=" * 72)
