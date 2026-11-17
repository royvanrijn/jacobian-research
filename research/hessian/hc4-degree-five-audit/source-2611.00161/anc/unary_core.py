# Noyau commun des verificateurs du cone unaire (HC4 quintique) : f = x1**5 + f4 + f3 + f2.
import time, itertools
from sage.all__sagemath_singular import *
import flint

T0 = time.time()
def say(*a):
    print(*a, flush=True)

NP = [0, 0]
def check(label, ok):
    NP[0 if ok else 1] += 1
    say("  [OK]   " if ok else "  [ECHEC]", label, " (%.0f s)" % (time.time() - T0))

def bilan():
    say("\n" + "=" * 72)
    say("BILAN : %d verifications passees, %d echecs  (%.0f s)" % (NP[0], NP[1], time.time() - T0))
    say("=" * 72)

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

M4 = [e for e in mon(4, 4) if e[0] < 4]   # x1**4 retire par translation de x1
M3 = mon(4, 3)
M2 = mon(4, 2)
VN = ['a%d' % i for i in range(len(M4))] + ['b%d' % i for i in range(len(M3))] + ['c%d' % i for i in range(len(M2))]
EXP = dict([('a%d' % i, e) for i, e in enumerate(M4)] + [('b%d' % i, e) for i, e in enumerate(M3)] + [('c%d' % i, e) for i, e in enumerate(M2)])

_ctx = flint.fmpq_mpoly_ctx.get(['x1', 'x2', 'x3', 'x4'] + VN, 'lex')
_G = _ctx.gens()
_X = _G[:4]
_V = dict(zip(['x1', 'x2', 'x3', 'x4'] + VN, _G))
def _m(e):
    return _X[0] ** e[0] * _X[1] ** e[1] * _X[2] ** e[2] * _X[3] ** e[3]
_F = _X[0] ** 5 + sum(_V[v] * _m(EXP[v]) for v in VN)
_D = det4([[_F.derivative(i).derivative(j) for j in range(4)] for i in range(4)])
RING = PolynomialRing(QQ, VN + ['T', 'U', 'W'], order='degrevlex')
RG = RING.gens()
RV = dict(zip(VN + ['T', 'U', 'W'], RG))
_cm = {}
for ex, co in zip(_D.monoms(), _D.coeffs()):
    _cm.setdefault(tuple(ex[:4]), {})[tuple(ex[4:])] = QQ(int(co.p)) / QQ(int(co.q))
LAYERS = {}
for k, d in _cm.items():
    pol = sum(cv * prod(RG[i] ** int(kk) for i, kk in enumerate(e) if kk) for e, cv in d.items())
    LAYERS.setdefault(sum(k), []).append(pol)
say("determinant du cone unaire : %d composantes (%.0f s)" % (len(_cm), time.time() - T0))

def is_one(L):
    return list(RING.ideal(L).groebner_basis()) == [1]

def _eqs(sub, tot):
    S = {RV[k]: RING(v) for k, v in sub.items()}
    out = []
    for t in tot:
        for p in LAYERS.get(t, []):
            e = p.subs(S)
            if e != 0:
                out.append(e)
    return out

def peel(sub, nz=(), tot=(9, 8, 7), kmax=6):
    # algorithme deterministe : equations lineaires resolues, zeros forces par Rabinowitsch
    # sur les equations a au plus kmax variables ; nz : quantites supposees non nulles.
    sub = dict(sub)
    while True:
        eqs = _eqs(sub, tot)
        if any(e.is_constant() for e in eqs):
            return True, sub
        lin = [e for e in eqs if e.degree() == 1]
        if lin:
            e = lin[0]
            v = e.variables()[0]
            c = e.monomial_coefficient(v)
            sub[str(v)] = -(e - c * v) / c
            sub = {k: (RING(w).subs({v: sub[str(v)]}) if k != str(v) else w) for k, w in sub.items()}
            continue
        S = {RV[k]: RING(v) for k, v in sub.items()}
        nzp = [RV['T'] * prod(RING(z).subs(S) for z in nz) - 1] if nz else []
        small = [e for e in eqs if len(e.variables()) <= kmax]
        if not small:
            return False, sub
        I = small + nzp
        if is_one(I):
            return True, sub
        cand = sorted(set(str(x) for e in small for x in e.variables()))
        found = [v for v in cand if is_one(I + [RV['U'] * RV[v] - 1])]
        if not found:
            return False, sub
        for v in found:
            sub[v] = 0

def implied(sub, target, nz=(), tot=(9,)):
    # target est dans le radical des couches tot (avec les quantites nz non nulles)
    eqs = _eqs(sub, tot)
    S = {RV[k]: RING(v) for k, v in sub.items()}
    t = RING(target).subs(S)
    if t == 0:
        return True
    nzp = [RV['T'] * prod(RING(z).subs(S) for z in nz) - 1] if nz else []
    return is_one(eqs + nzp + [RV['U'] * t - 1])

# monomes qui doivent s'annuler pour que x4 soit une direction "morte"
# (f_44 constante et coefficient de x4 de degre <= 2 : lemmes A et L)
BAD_X4 = [v for v in VN if (EXP[v][3] >= 2 and not (v[0] == 'c' and EXP[v] == (0, 0, 0, 2)))
          or (EXP[v][3] == 1 and sum(EXP[v]) - 1 >= 3)]
# f quadratique en z = (x3, x4) avec f_zz = A0 + x1 A1
BAD_Z = [v for v in VN if EXP[v][2] + EXP[v][3] >= 3
         or (EXP[v][2] + EXP[v][3] == 2 and not (sum(EXP[v]) == 2 or (sum(EXP[v]) == 3 and EXP[v][0] == 1)))]

def leaf_dead_x4(label, sub, nz=(), tot=(9, 8, 7), kmax=6):
    inc, fin = peel(sub, nz, tot, kmax)
    rest = [v for v in BAD_X4 if fin.get(v, None) != 0]
    check(label + " : x4 direction morte (lemme A ou L)", (not inc) and not rest)
    if rest and not inc:
        say("        restants :", rest)

def leaf_incoherent(label, sub, nz=(), tot=(9, 8, 7, 6), kmax=7):
    inc, fin = peel(sub, nz, tot, kmax)
    check(label + " : incoherent", inc)

def leaf_zquad(label, sub, rank2, nz=(), tot=(9, 8, 7, 6, 5), kmax=8):
    inc, fin = peel(sub, nz, tot, kmax)
    rest = [v for v in BAD_Z if fin.get(v, None) != 0]
    check(label + " : f quadratique en (x3,x4), f_zz = A0 + x1 A1", (not inc) and not rest)
    if rest and not inc:
        say("        restants :", rest)
    inc2, _ = peel(dict(fin, **rank2), nz, tot, kmax)
    check(label + " : A1 de rang 2 incoherent", inc2)

ZERO_A4 = {'a%d' % i: 0 for i in range(len(M4)) if EXP['a%d' % i][0] == 0}

def implied_c(sub, target, nz=()):
    # target dans le radical de toutes les couches, avec det Hess f = constante NON NULLE
    S = {RV[k]: RING(v) for k, v in sub.items()}
    eqs = _eqs(sub, [t for t in LAYERS if t > 0])
    cst = sum(p.subs(S) for p in LAYERS[0])
    t = RING(target).subs(S)
    if t == 0:
        return True
    nzp = [RV['T'] * prod(RING(z).subs(S) for z in nz) - 1] if nz else []
    return is_one(eqs + nzp + [RV['W'] * cst - 1, RV['U'] * t - 1])

_HS = PolynomialRing(RING, 'u2,u3,u4')
def hess_coeffs(vars_and_exps):
    # coefficients de det Hess d'une forme en (x2,x3,x4) dont les coefficients sont des inconnues
    u = _HS.gens()
    F = sum(RV[v] * u[0] ** e[1] * u[1] ** e[2] * u[2] ** e[3] for v, e in vars_and_exps)
    H = matrix(_HS, 3, 3, lambda i, j: F.derivative(u[i]).derivative(u[j]))
    return [RING(c) for c in H.det().coefficients()]

def forms(prefix_filter):
    return [(v, EXP[v]) for v in VN if prefix_filter(v, EXP[v])]
