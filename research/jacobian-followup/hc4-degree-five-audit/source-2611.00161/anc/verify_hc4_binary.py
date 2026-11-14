# Verificateur : HC4 quintique, cone binaire (sommet phi(x1,x2) non degenere).
# Rejoue toutes les verifications algebriques du cas binaire.
# Lancer : python3 verify_hc4_binary.py   (aucun argument)
import time, itertools, sys
T0 = time.time()
from sage.all__sagemath_singular import *
import flint

def say(*a):
    print(*a, flush=True)

NPASS = [0]
NFAIL = [0]
def check(label, ok):
    if ok:
        NPASS[0] += 1
        say("  [OK]   ", label, " (%.0f s)" % (time.time() - T0))
    else:
        NFAIL[0] += 1
        say("  [ECHEC]", label)

# ---------------------------------------------------------------------------
# 1. Construction des composantes du determinant pour trois ansatz
# ---------------------------------------------------------------------------
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

NAMES = []
for pref, k in [('p', 6), ('c', 4), ('d', 5), ('q', 10), ('r', 6), ('b', 4), ('m', 6), ('n', 3)]:
    NAMES += [pref + str(i) for i in range(k)]

def build(kind):
    names = ['x1', 'x2', 'x3', 's'] + NAMES
    ctx = flint.fmpq_mpoly_ctx.get(names, 'lex')
    G = ctx.gens()
    V = dict(zip(names, G))
    x1, x2, x3, s = G[:4]
    bm = lambda e: x1 ** e[0] * x2 ** e[1]
    tm = lambda e: x1 ** e[0] * x2 ** e[1] * x3 ** e[2]
    P = sum(V['p%d' % i] * bm(e) for i, e in enumerate(mon(2, 5)))
    P += x3 * sum(V['c%d' % i] * bm(e) for i, e in enumerate(mon(2, 3)))
    P += sum(V['d%d' % i] * bm(e) for i, e in enumerate(mon(2, 4)))
    P += sum(V['q%d' % i] * tm(e) for i, e in enumerate(mon(3, 3)))
    P += sum(V['r%d' % i] * tm(e) for i, e in enumerate(mon(3, 2)))
    R = sum(V['b%d' % i] * bm(e) for i, e in enumerate(mon(2, 3)))
    R += sum(V['m%d' % i] * tm(e) for i, e in enumerate(mon(3, 2)))
    R += sum(V['n%d' % i] * tm(e) for i, e in enumerate(mon(3, 1)))
    if kind == 'cube':
        F = P + s * R + s ** 3
    elif kind == 'e0':
        F = P + s * R + s ** 2 * x1
    else:
        F = P + s * R
    H = [[F.derivative(i).derivative(j) for j in range(4)] for i in range(4)]
    D = det4(H)
    comps = {}
    for ex, c in zip(D.monoms(), D.coeffs()):
        key = (ex[3], ex[0] + ex[1], ex[2], ex[0], ex[1])   # (puissance de s, degre en X, degre en x3, monome)
        comps.setdefault(key, {})[tuple(ex[4:])] = QQ(int(c.p)) / QQ(int(c.q))
    return comps

RING = PolynomialRing(QQ, NAMES + ['T', 'TT', 'U'], order='degrevlex')
GV = dict(zip(NAMES + ['T', 'TT', 'U'], RING.gens()))
RG = RING.gens()

def topoly(d):
    return sum(c * prod(RG[i] ** int(k) for i, k in enumerate(ex) if k) for ex, c in d.items())

COMPS = {}
POLYS = {}
def setup(kind):
    if kind not in COMPS:
        COMPS[kind] = build(kind)
        POLYS[kind] = {k: topoly(v) for k, v in COMPS[kind].items()}

def system(kind, sub, extra=(), lay=None):
    S = {GV[k]: QQ(v) for k, v in sub.items()}
    eqs = []
    cst = None
    for k, p in POLYS[kind].items():
        e = p.subs(S)
        if k[0] == 0 and k[1] == 0 and k[2] == 0:
            cst = e
        elif e != 0 and (lay is None or k[0] in lay or (k[0], k[1]) in lay):
            eqs.append(e)
    eqs += [RING(x).subs(S) for x in extra]
    return eqs, cst, S

def is_one(I):
    return list(I.groebner_basis()) == [1]

def incoherent(kind, sub, extra=()):
    eqs, cst, S = system(kind, sub, extra)
    return is_one(RING.ideal(eqs + [GV['TT'] * cst - 1]))

PRCACHE = {}
def primes_of(kind, sub, extra, lay):
    key = (kind, tuple(sorted(sub.items())), tuple(str(x) for x in extra), lay)
    if key not in PRCACHE:
        eqs, cst, S = system(kind, sub, extra, lay)
        extra_c = [GV['TT'] * cst - 1] if lay is None else []
        PRCACHE[key] = (RING.ideal(eqs + extra_c).minimal_associated_primes(), S)
    return PRCACHE[key]

def implies(kind, sub, target, extra=(), lay=None):
    # target dans le radical (via les premiers associes minimaux)
    P, S = primes_of(kind, sub, extra, lay)
    t = RING(target).subs(S)
    return all(Q.reduce(t) == 0 for Q in P)

def implies_nd(kind, sub, target, extra=(), lay=None):
    # target = 0 des que phi est non degeneree : pour chaque coefficient g de det Hess phi,
    # 1 est dans (equations) + (T g - 1) + (U target - 1)
    eqs, cst, S = system(kind, sub, extra, lay)
    t = RING(target).subs(S)
    if t == 0:
        return True
    extra_c = [GV['TT'] * cst - 1] if lay is None else []
    for g in phi_coeffs(S):
        if g == 0:
            continue
        if not is_one(RING.ideal(eqs + extra_c + [GV['T'] * g - 1, GV['U'] * t - 1])):
            return False
    return True

def degenerate_lay(kind, sub, lay, extra=()):
    eqs, cst, S = system(kind, sub, extra, lay)
    return all(implies(kind, sub, g, extra, lay) for g in phi_coeffs(S) if g != 0)

def phi_coeffs(S):
    X1, X2 = PolynomialRing(RING, 'X1,X2').gens()
    phi = sum(GV['p%d' % i] * X1 ** e[0] * X2 ** e[1] for i, e in enumerate(mon(2, 5)))
    dH = phi.derivative(X1, 2) * phi.derivative(X2, 2) - phi.derivative(X1).derivative(X2) ** 2
    return [RING(c).subs(S) for c in dH.coefficients()]

def degenerate(kind, sub, extra=()):
    # phi non degeneree impossible : chaque coefficient de det Hess phi est dans le radical
    eqs, cst, S = system(kind, sub, extra)
    for g in phi_coeffs(S):
        if g == 0:
            continue
        if not is_one(RING.ideal(eqs + [GV['TT'] * cst - 1, GV['T'] * g - 1])):
            return False
    return True

say("=" * 72)
say("Verificateur HC4 quintique, cone binaire")
say("=" * 72)

# ---------------------------------------------------------------------------
# 2. Reductions symboliques
# ---------------------------------------------------------------------------
say("\n[A] Reductions generales")
import sympy as sp
x1, x2, x3, s, t = sp.symbols('x1 x2 x3 s t')
Pf, Rf, Lf = [sp.Function(n)(x1, x2, x3) for n in 'PRL']
f = Pf + s * Rf + s ** 2 * Lf + s ** 3
g = sp.expand(f.subs(s, t - Lf / 3))
tgt = t ** 3 + t * (Rf - Lf ** 2 / 3) + (Pf - Lf * Rf / 3 + 2 * Lf ** 3 / 27)
check("A1 translation s -> s - L/3 : f = P' + s R' + s**3", sp.simplify(g - sp.expand(tgt)) == 0)

# identite de Schur utilisee pour R = x1**3 + x3
h = sp.symbols('h11 h12 h13 h22 h23 h33')
H = sp.Matrix([[h[0], h[1], h[2]], [h[1], h[3], h[4]], [h[2], h[4], h[5]]])
adj11 = h[3] * h[5] - h[4] ** 2
adj13 = h[1] * h[4] - h[2] * h[3]
ident = h[3] * H.det() + adj13 ** 2 - adj11 * (h[0] * h[3] - h[1] ** 2)
check("A2 identite h22 det H + adj13**2 = adj11 (h11 h22 - h12**2)", sp.expand(ident) == 0)

# formule du determinant quand R ne depend pas de w (lemme direction morte)
P2 = sp.Function('P')(x1, x2, x3)
R2 = sp.Function('R')(x1, x2)
F2 = P2 + s * R2
Vv = [x1, x2, x3, s]
HH = sp.Matrix(4, 4, lambda i, j: sp.diff(F2, Vv[i], Vv[j]))
Q1 = sp.diff(R2, x2) ** 2 * sp.diff(R2, x1, 2) - 2 * sp.diff(R2, x1) * sp.diff(R2, x2) * sp.diff(R2, x1, x2) + sp.diff(R2, x1) ** 2 * sp.diff(R2, x2, 2)
Q0 = sp.diff(R2, x2) ** 2 * sp.diff(P2, x1, 2) - 2 * sp.diff(R2, x1) * sp.diff(R2, x2) * sp.diff(P2, x1, x2) + sp.diff(R2, x1) ** 2 * sp.diff(P2, x2, 2)
form = -sp.diff(P2, x3, 2) * (Q0 + s * Q1) + (sp.diff(R2, x1) * sp.diff(P2, x2, x3) - sp.diff(R2, x2) * sp.diff(P2, x1, x3)) ** 2
check("A3 det Hess(P + s R(x1,x2)) = -P_33 (Q0 + s Q1) + (R1 P_23 - R2 P_13)**2", sp.expand(HH.det() - form) == 0)

# lemme de Legendre : instance explicite avec f_zz = [[0,1],[1,2 X1]] (det -1) et b, c quelconques
X1s, X2s, z1, z2, y1, y2 = sp.symbols('X1 X2 z1 z2 y1 y2')
b1 = 3*X1s**3*X2s - X2s**2 + 2*X1s
b2 = X1s**2*X2s + 5*X2s**3 - X1s
c0 = X1s**5 - 2*X1s**2*X2s**3 + 7*X2s**4 + X1s*X2s
fq = z1*z2 + X1s*z2**2 + b1*z1 + b2*z2 + c0
Hq = sp.Matrix(4, 4, lambda i, j: sp.diff(fq, [X1s, X2s, z1, z2][i], [X1s, X2s, z1, z2][j]))
A = sp.Matrix([[0, 1], [1, 2*X1s]])
zs = A.inv() * sp.Matrix([y1 - b1, y2 - b2])
Phi = sp.expand((fq - y1*z1 - y2*z2).subs({z1: zs[0], z2: zs[1]}))
HPhi = sp.Matrix(2, 2, lambda i, j: sp.diff(Phi, [X1s, X2s][i], [X1s, X2s][j]))
lhs = sp.expand(Hq.det().subs({z1: zs[0], z2: zs[1]}))
check("A4 Legendre partiel, controle sur une instance (identite de Schur generale) : det Hess f = det f_zz * det Hess Phi_y sur z = z(X)",
      sp.expand(lhs - A.det()*HPhi.det()) == 0)

# ---------------------------------------------------------------------------
# 3. Branche e != 0 : f = P + s R + s**3
# ---------------------------------------------------------------------------
say("\n[B] Branche e != 0 : f = P + s R + s**3")
setup('cube')
GEN = {'b0': 1, 'b1': 0, 'b2': 0, 'b3': 1}
for v in ['m0', 'm1', 'm3', 'q0', 'q1', 'q4', 'r0']:
    check("B1 B3 = x1**3 + x2**3 : %s = 0" % v, implies('cube', GEN, GV[v], lay=(4, 3)))
G1 = dict(GEN, m0=0, m1=0, m3=0, q0=0, q1=0, q4=0, r0=0, m2=0, m5=0)
for v in ['c0', 'c1', 'c2', 'c3', 'q2', 'q5', 'q7', 'n0', 'r1', 'r3']:
    check("B1 (normalise m2 = m5 = 0) : %s = 0" % v, implies('cube', G1, GV[v], lay=(4, 3, 2)))
check("B1 alors f ne depend plus de x3 : incoherent",
      incoherent('cube', dict(G1, c0=0, c1=0, c2=0, c3=0, q2=0, q5=0, q7=0, n0=0, r1=0, r3=0)))

B21 = {'b0': 0, 'b1': 0, 'b2': 1, 'b3': 0}
for v in ['m0', 'm1']:
    check("B2 B3 = x1**2 x2 : %s = 0" % v, implies('cube', B21, GV[v], lay=(4,)))
check("B2 B3 = x1**2 x2 : m2 m3 = 0", implies('cube', dict(B21, m0=0, m1=0), GV['m2'] * GV['m3'], lay=(4,)))
check("B2a terme x1 x3 (m3 = 1, m2 = m4 = m5 = 0) : incoherent",
      incoherent('cube', dict(B21, m0=0, m1=0, m3=1, m2=0, m4=0, m5=0)))
B2b = dict(B21, m0=0, m1=0, m3=0, m4=0, m5=0)
for m2 in [1, 0]:
    for v in ['q0', 'q1', 'q4', 'r0']:
        check("B2b m2 = %d : %s = 0 (P lineaire en x3)" % (m2, v), implies('cube', dict(B2b, m2=m2), GV[v], lay=(4, 3)))
Bb1 = dict(B2b, m2=1, q0=0, q1=0, q4=0, r0=0)
for v in ['c0', 'c1', 'c2', 'c3', 'q2', 'q5', 'r3']:
    check("B2b m2 = 1 : %s = 0" % v, implies('cube', Bb1, GV[v], lay=(4, 3, 2)))
check("B2b m2 = 1 : r1 = 2 q7", implies('cube', Bb1, GV['r1'] - 2 * GV['q7'], lay=(4, 3, 2)))
check("B2b m2 = 1 : n0**2 + 12 q7**2 = 0", implies('cube', Bb1, GV['n0'] ** 2 + 12 * GV['q7'] ** 2, lay=(4, 3, 2)))
check("B2b m2 = 1 : q7 = 0 incoherent", incoherent('cube', dict(Bb1, q7=0)))
check("B2b m2 = 0 : incoherent (branche c1 != 0)",
      incoherent('cube', dict(B2b, m2=0, q0=0, q1=0, q4=0, r0=0, c0=0, c2=0, r1=0, q2=0, q5=0, n0=0, c1=1, c3=0, q7=0, r3=0)))
for v in ['c0', 'c2', 'r1', 'q2', 'q5', 'n0']:
    check("B2b m2 = 0 : %s = 0" % v, implies('cube', dict(B2b, m2=0, q0=0, q1=0, q4=0, r0=0), GV[v], lay=(4, 3, 2)))
check("B2b m2 = 0 : c1 c3 = 0 et c1 q7 = 0 et c1 r3 = 0",
      all(implies('cube', dict(B2b, m2=0, q0=0, q1=0, q4=0, r0=0, c0=0, c2=0, r1=0, q2=0, q5=0, n0=0), e, lay=(4, 3, 2))
          for e in [GV['c1'] * GV['c3'], GV['c1'] * GV['q7'], GV['c1'] * GV['r3']]))
check("B2b m2 = 0 : branche c1 = 0 incoherente",
      incoherent('cube', dict(B2b, m2=0, q0=0, q1=0, q4=0, r0=0, c0=0, c2=0, r1=0, q2=0, q5=0, n0=0, c1=0)))

# famille survivante B2b, m2 = 1 : relations, determinant, inverse explicite
K = QuadraticField(-3, 'w')
w = K.gen()
RK = PolynomialRing(K, NAMES + ['T', 'TT', 'U'], order='degrevlex')
GK = dict(zip(NAMES + ['T', 'TT', 'U'], RK.gens()))
def kimplies(sub, target):
    S = {GK[k]: K(v) for k, v in sub.items()}
    eqs = []
    cst = None
    for k, p in POLYS['cube'].items():
        e = RK(p).subs(S)
        if k[0] == 0 and k[1] == 0 and k[2] == 0:
            cst = e
        elif e != 0:
            eqs.append(e)
    tt = RK(target).subs(S)
    if tt == 0:
        return True
    return list(RK.ideal(eqs + [GK['TT'] * cst - 1, GK['T'] * tt - 1]).groebner_basis()) == [1]
for sg in [1, -1]:
    base = dict(B2b, m2=1, q0=0, q1=0, q4=0, r0=0, c0=0, c1=0, c2=0, c3=0, q2=0, q5=0, r3=0, q7=1, r1=2)
    base['n0'] = sg * 2 * w
    rel = [3 * GK['r2'] + sg * w * GK['n1'], 6 * GK['q8'] + sg * w * GK['n1'], GK['q6'], 9 * GK['q3'] + sg * 2 * w,
           GK['d3'], 3 * GK['d2'] + sg * w, 12 * GK['p4'] + sg * w, GK['p0'], GK['p1'], GK['p2'], GK['p3'], GK['d0'], GK['d1']]
    check("B2b famille (n0 = %s) : relations lineaires" % (sg * 2 * w), all(kimplies(base, r) for r in rel))
A7 = PolynomialRing(K, 'n1,n2,r4,p5,d4,q9,r5')
S7 = PolynomialRing(A7, 'x1,z,x3,s')
n1, n2, r4, p5, d4, q9, r5 = [S7(g) for g in A7.gens()]
X1, Z, X3, Ss = S7.gens()
X2 = Z - w * Ss
Pfam = (-w / 12) * X1 ** 4 * X2 + p5 * X1 ** 5 + (-w / 3) * X1 ** 2 * X2 ** 2 + d4 * X1 ** 4 + (-2 * w / 9) * X2 ** 3 + X1 ** 2 * X3 \
    + (-w * n1 / 6) * X1 ** 2 * X2 + q9 * X1 ** 3 + 2 * X2 * X3 + (-w * n1 / 3) * X2 ** 2 + r4 * X1 * X2 + r5 * X1 ** 2
Rfam = X1 ** 2 * X2 + X2 ** 2 + 2 * w * X3 + n1 * X2 + n2 * X1
ffam = Pfam + Ss * Rfam + Ss ** 3
fs = ffam.polynomial(Ss)
mm = Z + X1 ** 2 / 2
check("B2b famille : en (x1, z = x2 + w s, x3, s), f = 2 m x3 + G(x1,z) + s (kappa x1 - m**2 - n1 m)",
      ffam.derivative(Ss, 2) == 0 and ffam.derivative(Ss).subs({Ss: 0}) == (n2 - w * r4) * X1 - mm ** 2 - n1 * mm
      and ffam.derivative(X3) == 2 * mm)
# inverse explicite : y3 = 2 m, y4 = kappa x1 - m**2 - n1 m, puis (x3, s) par un systeme lineaire de det 2 kappa
kap = n2 - w * r4
J2 = matrix(S7, 2, 2, [ffam.derivative(Z).derivative(X3), ffam.derivative(Z).derivative(Ss),
                       ffam.derivative(X1).derivative(X3), ffam.derivative(X1).derivative(Ss)])
check("B2b famille : (f_z, f_x1) affine en (x3, s), matrice de det 2 kappa", J2.det() == 2 * kap and
      all(e.derivative(X3).derivative(X3) == 0 and e.derivative(Ss).derivative(Ss) == 0 and e.derivative(X3).derivative(Ss) == 0
          for e in [ffam.derivative(Z), ffam.derivative(X1)]))
Hf = matrix(S7, 4, 4, lambda i, j: ffam.derivative([X1, Z, X3, Ss][i]).derivative([X1, Z, X3, Ss][j]))
say("  (n0 = -2w : conjugaison galoisienne de Q(sqrt(-3)))")
check("B2b famille : det Hess f = 4 (n2 - w r4)**2", Hf.det() == 4 * (n2 - w * r4) ** 2)

# B3 = x1**3
C3 = {'b0': 0, 'b1': 0, 'b2': 0, 'b3': 1}
DQ = (-GV['m2'] * GV['m3'] ** 2 + GV['m1'] * GV['m3'] * GV['m4'] - GV['m0'] * GV['m4'] ** 2
      - GV['m1'] ** 2 * GV['m5'] + 4 * GV['m0'] * GV['m2'] * GV['m5'])
check("B3 B3 = x1**3 : det H_Q = 0 (couche s**4)", implies('cube', C3, DQ, lay=(4,)))
check("B3 Q = x3**2 + x1 x2 : incoherent", incoherent('cube', dict(C3, m0=1, m1=0, m2=0, m3=0, m4=1, m5=0)))
check("B3 Q = x2**2 + x1 x3 : incoherent", incoherent('cube', dict(C3, m0=0, m1=0, m2=1, m3=1, m4=0, m5=0)))
check("B3 R = x1**3 + n2 x1 : incoherent", incoherent('cube', dict(C3, m0=0, m1=0, m2=0, m3=0, m4=0, m5=0, n0=0, n1=0)))
check("B3 B3 = x1**3 : m1**2 = 4 m0 m2", implies('cube', C3, GV['m1'] ** 2 - 4 * GV['m0'] * GV['m2'], lay=(4,)))
Z0 = {'m0': 0, 'm1': 0, 'm2': 0, 'm3': 0, 'm4': 0, 'm5': 0}
check("B3 Q = x3**2 : phi degeneree", degenerate('cube', dict(C3, **dict(Z0, m0=1))))
X22 = dict(C3, **dict(Z0, m2=1))
for v in ['q0', 'q1', 'q4', 'r0']:
    check("B3 Q = x2**2 : %s = 0" % v, implies('cube', X22, GV[v], lay=(4, 3)))
X22b = dict(X22, q0=0, q1=0, q4=0, r0=0)
for v in ['c0', 'c1', 'c2', 'c3', 'q2', 'q5', 'q7', 'r3']:
    check("B3 Q = x2**2 : %s = 0" % v, implies('cube', X22b, GV[v], lay=(4, 3, 2)))
X22c = dict(X22b, c0=0, c1=0, c2=0, c3=0, q2=0, q5=0, q7=0, r3=0)
check("B3 Q = x2**2 : n0**2 + 3 r1**2 = 0", implies('cube', X22c, GV['n0'] ** 2 + 3 * GV['r1'] ** 2, lay=(4, 3, 2)))
check("B3 Q = x2**2 : r1 = 0 incoherent", incoherent('cube', dict(X22c, r1=0, n0=0)))
check("B3 Q = x2**2 (r1 = 1) : phi degeneree", degenerate('cube', dict(X22c, r1=1), extra=[GV['n0'] ** 2 + 3]))
check("B3 Q = x1 x3 : incoherent", incoherent('cube', dict(C3, **dict(Z0, m3=1))))
X12 = dict(C3, **dict(Z0, m4=1))
for v in ['q0', 'q1', 'q4', 'r0']:
    check("B3 Q = x1 x2 : %s = 0" % v, implies('cube', X12, GV[v], lay=(4, 3)))
X12b = dict(X12, q0=0, q1=0, q4=0, r0=0)
for v in ['c0', 'c1', 'c2', 'q2', 'q5']:
    check("B3 Q = x1 x2 : %s = 0" % v, implies('cube', X12b, GV[v], lay=(4, 3, 2)))
X12c = dict(X12b, c0=0, c1=0, c2=0, q2=0, q5=0)
check("B3 Q = x1 x2 : c3 r1 = 0", implies('cube', X12c, GV['c3'] * GV['r1'], lay=(4, 3, 2)))
check("B3 Q = x1 x2 : n0**2 + 12 r1 r3 = 0", implies('cube', X12c, GV['n0'] ** 2 + 12 * GV['r1'] * GV['r3'], lay=(4, 3, 2)))
check("B3 Q = x1 x2, c3 = 0 : r1 (2 q7 - 3 r1) = 0",
      implies('cube', dict(X12c, c3=0), GV['r1'] * (2 * GV['q7'] - 3 * GV['r1']), lay=(4, 3, 2)))
check("B3 Q = x1 x2, r1 = 0 (donc n0 = 0) : incoherent", incoherent('cube', dict(X12c, r1=0, n0=0)))
check("B3 Q = x1 x2, c3 = 0, r1 = 2, q7 = 3 : incoherent",
      incoherent('cube', dict(X12c, c3=0, r1=2, q7=3), extra=[GV['n0'] ** 2 + 24 * GV['r3']]))
check("B3 R = x1**3 + x2 : composante (P lineaire en x2 et en x3, x3 seulement via r3 x1 x3)",
      all(implies('cube', dict(C3, **dict(Z0, n0=0, n1=1, n2=0)), GV[v]) for v in
          ['p0', 'p1', 'p2', 'p3', 'c0', 'c1', 'c2', 'c3', 'd0', 'd1', 'd2', 'q0', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q7', 'r0', 'r1', 'r2']))

# R = x1**3 + x3 : coefficients en s du determinant (P generique)
Pg = sp.Function('P')(x1, x2, x3)
Rg = x1 ** 3 + x3
Fg = Pg + s * Rg + s ** 3
Hg = sp.Matrix(4, 4, lambda i, j: sp.diff(Fg, Vv[i], Vv[j]))
HP = sp.Matrix(3, 3, lambda i, j: sp.diff(Pg, [x1, x2, x3][i], [x1, x2, x3][j]))
dg = sp.Poly(sp.expand(Hg.det()), s)
a11 = HP[1, 1] * HP[2, 2] - HP[1, 2] ** 2
check("B3 R = x1**3 + x3 : coefficient de s**4 et s**3 nuls", sp.expand(dg.coeff_monomial(s ** 4)) == 0 and sp.expand(dg.coeff_monomial(s ** 3)) == 0)
check("B3 R = x1**3 + x3 : coefficient de s**2 = 36 x1 adj11", sp.expand(dg.coeff_monomial(s ** 2) - 36 * x1 * a11) == 0)
check("B3 R = x1**3 + x3 : coefficient de s = 6 det H_P - 6 x1 P_22",
      sp.expand(dg.coeff_monomial(s) - 6 * HP.det() + 6 * x1 * HP[1, 1]) == 0)

# B3 = 0
Zb = {'b0': 0, 'b1': 0, 'b2': 0, 'b3': 0, 'm0': 0, 'm1': 0, 'm2': 0, 'm3': 0, 'm4': 0, 'm5': 0, 'n0': 0, 'n1': 0, 'n2': 0}
def zb(**kw):
    d = dict(Zb)
    d.update(kw)
    return d
for nm, sb in [("x3**2 + x1**2", zb(m0=1, m5=1)), ("x3**2", zb(m0=1)), ("x1 x3", zb(m3=1)), ("x1 x2", zb(m4=1)),
               ("x1**2", zb(m5=1)), ("0", zb()), ("x3**2 + x1 x2", zb(m0=1, m4=1)), ("x1 x3 + x2**2", zb(m3=1, m2=1))]:
    check("B4 R = %s (sans partie lineaire) : incoherent" % nm, incoherent('cube', sb))
check("B4 R = x3**2 + x1**2 + x2 : P_22 = 0 (couches s**4, s**3)",
      all(implies('cube', zb(m0=1, m5=1, n1=1), GV[v], lay=(4, 3)) for v in
          ['p0', 'p1', 'p2', 'p3', 'c0', 'c1', 'd0', 'd1', 'd2', 'q2', 'q3', 'q6', 'r2']))
check("B4 R = x3**2 + x1**2 + x2 : puis p4 = 0 (phi degeneree)",
      implies('cube', zb(m0=1, m5=1, n1=1, p0=0, p1=0, p2=0, p3=0, c0=0, c1=0, d0=0, d1=0, d2=0, q2=0, q3=0, q6=0, r2=0), GV['p4'], lay=(4, 3, 2)))
check("B4 R = x3**2 + x1 : phi degeneree (couche s**2, degre 6 en X)", degenerate_lay('cube', zb(m0=1, n2=1), lay=((2, 6),)))
check("B4 R = x1 x3 + x2 : phi degeneree", degenerate('cube', zb(m3=1, n1=1)))
check("B4 R = x1 x2 + x3 : phi degeneree", degenerate('cube', zb(m4=1, n0=1)))
check("B4 R = x3 : phi degeneree (couche s**0, degre 6 en X)", degenerate_lay('cube', zb(n0=1), lay=((0, 6),)))
for nm, sb in [("x1**2 + x3", zb(m5=1, n0=1)), ("x1**2 + x2", zb(m5=1, n1=1))]:
    check("B4 R = %s : q0 = 0 si phi non degeneree (couche s**1)" % nm, implies_nd('cube', sb, GV['q0'], lay=(1,)))
    sb1 = dict(sb, q0=0)
    check("B4 R = %s : q1 = 0" % nm, implies('cube', sb1, GV['q1'], lay=(2,)))
    sb2 = dict(sb1, q1=0)
    check("B4 R = %s : c0 = 0" % nm, implies('cube', sb2, GV['c0'], lay=(2,)))
    sb3 = dict(sb2, c0=0)
    check("B4 R = %s, q4 = 1 : p0 = c1 = p1 = p2 = 0" % nm,
          all(implies('cube', dict(sb3, q4=1), GV[v], lay=(2,)) for v in ['p0', 'c1', 'p1', 'p2']))
    check("B4 R = %s, q4 = 1 : incoherent" % nm, incoherent('cube', dict(sb3, q4=1, p0=0, c1=0, p1=0, p2=0)))
    check("B4 R = %s, q4 = 0 : c1 = c2 = 0" % nm,
          all(implies('cube', dict(sb3, q4=0), GV[v], lay=(2,)) for v in ['c1', 'c2']))
check("B4 R = x1**2 + x3, q4 = 0 : phi degeneree", degenerate('cube', dict(zb(m5=1, n0=1), q0=0, q1=0, c0=0, q4=0, c1=0, c2=0)))
check("B4 R = x1**2 + x2, q4 = 0 : composante non degeneree triangulaire",
      all(implies_nd('cube', dict(zb(m5=1, n1=1), q0=0, q1=0, c0=0, q4=0, c1=0, c2=0), GV[v]) for v in
          ['r2', 'r1', 'r0', 'q7', 'q6', 'q5', 'q3', 'q2', 'd2', 'd1', 'd0', 'c3', 'p3', 'p2', 'p1', 'p0']))
Fx = Pg + s * x1 + s ** 3
Hx = sp.Matrix(4, 4, lambda i, j: sp.diff(Fx, Vv[i], Vv[j]))
check("B4 R = x1 : det Hess f = 6 s det H_P - adj11(H_P)", sp.expand(Hx.det() - 6 * s * HP.det() + a11) == 0)
say("\n[fin de la partie B : %.0f s]" % (time.time() - T0))

# ---------------------------------------------------------------------------
# 4. Branche e = 0, L = x1 : f = P + s R + s**2 x1
# ---------------------------------------------------------------------------
say("\n[C] Branche e = 0, L = x1 : f = P + s R + s**2 x1")
setup('e0')
check("C1 m0 != 0 : b0 = b1 = 0 et p0..p3 = 0",
      all(implies('e0', {'m0': 1}, GV[v], lay=(4, 3)) for v in ['b0', 'b1', 'p0', 'p1', 'p2', 'p3']))
E1 = {'m0': 1, 'm1': 0, 'm2': 0, 'm3': 0, 'b0': 0, 'b1': 0, 'p0': 0, 'p1': 0, 'p2': 0, 'p3': 0}
check("C1 m0 = 1 (normalise m1 = m3 = 0) : m2 = 0", implies('e0', {'m0': 1, 'm1': 0, 'm3': 0, 'b0': 0, 'b1': 0}, GV['m2'], lay=(4,)))
check("C1 m0 = 1 : b2 = 0", implies('e0', E1, GV['b2'], lay=(4, 3, 2)))
check("C1 m0 = 1 : p4 = 0 (phi degeneree)", implies('e0', dict(E1, b2=0), GV['p4'], lay=((1, 7),)))
check("C2 m0 = 0 : m1 = 0", implies('e0', {'m0': 0}, GV['m1'], lay=(4,)))
E2 = {'m0': 0, 'm1': 0, 'm3': 0}
check("C2 q0 = 1 (normalise q1 = q4 = r0 = 0) : phi = 0",
      all(implies('e0', dict(E2, q0=1, q1=0, q4=0, r0=0), GV[v]) for v in ['p0', 'p1', 'p2', 'p3', 'p4', 'p5']))
E3 = dict(E2, q0=0)
check("C3 q0 = 0 : q1 = 0", implies('e0', E3, GV['q1'], lay=((2, 0),)))
E3 = dict(E3, q1=0)
check("C3 q4 = 1 : c0 = c1 = p0 = p1 = p2 = b0 = b1 = 0",
      all(implies('e0', dict(E3, q4=1), GV[v], lay=((2, 4), (2, 2), (3, 2))) for v in ['c0', 'c1', 'p0', 'p1', 'p2', 'b0', 'b1']))
E4 = dict(E3, q4=1, c0=0, c1=0, p0=0, p1=0, p2=0, b0=0, b1=0)
check("C3 q4 = 1 : m2 = q2 = d0 = d1 = b2 = 0",
      all(implies('e0', E4, GV[v], lay=((3, 1), (2, 1), (2, 3), (2, 4), (1, 6), (1, 4))) for v in ['m2', 'q2', 'd0', 'd1', 'b2']))
check("C3 q4 = 1 : incoherent", incoherent('e0', dict(E4, m2=0, q2=0, d0=0, d1=0, b2=0)))
E5 = dict(E3, q4=0)
check("C4 q4 = 0 : c0 = c1 = c2 = 0", all(implies('e0', E5, GV[v], lay=((2, 4),)) for v in ['c0', 'c1', 'c2']))
E5 = dict(E5, c0=0, c1=0, c2=0)
check("C4 r0 = 1 : b0 = b1 = m2 = p0 = p1 = q2 = p2 = p3 = 0",
      all(implies('e0', dict(E5, r0=1), GV[v], lay=((3, 1), (3, 0), (2, 3), (2, 0))) for v in ['b0', 'b1', 'm2', 'p0', 'p1', 'q2', 'p2', 'p3']))
E6 = dict(E5, r0=1, b0=0, b1=0, m2=0, p0=0, p1=0, q2=0, p2=0, p3=0)
check("C4 r0 = 1 : b2 = d0 = d1 = q3 = 0", all(implies('e0', E6, GV[v], lay=(2, (1, 5))) for v in ['b2', 'd0', 'd1', 'q3']))
check("C4 r0 = 1 : p4 = 1 incoherent (phi degeneree)", incoherent('e0', dict(E6, b2=0, d0=0, d1=0, q3=0, p4=1)))
E7 = dict(E5, r0=0)
check("C5 r0 = 0 : c3 b0 = c3 b1 = 0", implies('e0', E7, GV['c3'] * GV['b0'], lay=((1, 6),)) and implies('e0', E7, GV['c3'] * GV['b1'], lay=((1, 6),)))
check("C5 r0 = 0, c3 = 1, n0 = 0 : incoherent", incoherent('e0', dict(E7, c3=1, b0=0, b1=0, n0=0)))
check("C5 r0 = 0, n0 = 0, c3 = 0 : composante non degeneree (lemme K(i), x3 morte, R = r3 x1)",
      all(implies_nd('e0', dict(E7, c3=0, n0=0), GV[v]) for v in
          ['m2', 'b1', 'b0', 'r1', 'q7', 'q5', 'q3', 'q2', 'd1', 'd0', 'p2', 'p1', 'p0']))
say("  (r0 = 0, n0 != 0 : f quadratique en (x3, s), det f_zz = -n0**2 : lemme de Legendre, A4)")

# ---------------------------------------------------------------------------
# 5. Branche e = 0, L = 0 : f = P + s R
# ---------------------------------------------------------------------------
say("\n[D] Branche e = 0, L = 0 : f = P + s R")
setup('lin')
check("D1 B3 = x1**3 + x2**3 : R ne depend pas de x3 (m0 = m1 = m3 = n0 = 0)",
      all(implies('lin', GEN, GV[v], lay=(2,)) for v in ['m0', 'm1', 'm3', 'n0']))
check("D2 B3 = x1**2 x2 : m0 = m1 = m3 = 0", all(implies('lin', B21, GV[v], lay=(2,)) for v in ['m0', 'm1', 'm3']))
say("  (cas D2/D3 : voir verify_hc4_prelim.py, partie P4, avec A4 complet)")
DELTA = (-GV['m2'] * GV['m3'] ** 2 + GV['m1'] * GV['m3'] * GV['m4'] - GV['m0'] * GV['m4'] ** 2
         - GV['m1'] ** 2 * GV['m5'] + 4 * GV['m0'] * GV['m2'] * GV['m5'])
BORD = (-GV['m4'] ** 2 * GV['n0'] ** 2 + 4 * GV['m2'] * GV['m5'] * GV['n0'] ** 2 + 2 * GV['m3'] * GV['m4'] * GV['n0'] * GV['n1']
        - 4 * GV['m1'] * GV['m5'] * GV['n0'] * GV['n1'] - GV['m3'] ** 2 * GV['n1'] ** 2 + 4 * GV['m0'] * GV['m5'] * GV['n1'] ** 2
        - 4 * GV['m2'] * GV['m3'] * GV['n0'] * GV['n2'] + 2 * GV['m1'] * GV['m4'] * GV['n0'] * GV['n2'] + 2 * GV['m1'] * GV['m3'] * GV['n1'] * GV['n2']
        - 4 * GV['m0'] * GV['m4'] * GV['n1'] * GV['n2'] - GV['m1'] ** 2 * GV['n2'] ** 2 + 4 * GV['m0'] * GV['m2'] * GV['n2'] ** 2)
ZB = {'b0': 0, 'b1': 0, 'b2': 0, 'b3': 0}
check("D4 B3 = 0 : v * det(Q) = 0 pour chaque coefficient v de R (couche s**2)",
      all(implies('lin', ZB, GV[v] * DELTA, lay=(2,)) for v in ['m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'n0', 'n1', 'n2']))
check("D4 B3 = 0 : bordant n**T adj(Q) n = 0 (couche s**2)", implies('lin', ZB, BORD, lay=(2,)))

say("\n" + "=" * 72)
say("BILAN : %d verifications passees, %d echecs  (%.0f s)" % (NPASS[0], NFAIL[0], time.time() - T0))
say("=" * 72)
