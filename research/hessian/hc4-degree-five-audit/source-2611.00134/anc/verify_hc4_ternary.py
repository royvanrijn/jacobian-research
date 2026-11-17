"""
verify_hc4_ternary.py  --  exact verification of the computational steps in the proof that
HC_4 holds for every quintic potential whose top form is a genuinely ternary cone.

Author: David Niedbala Giraudin (ORCID 0009-0009-1526-1178); computations assisted by Claude (Anthropic).
Requirements: python-flint, sympy, and Singular through Sage (passagemath-singular or a SageMath install).
Run without arguments. Every check prints PASS or FAIL; the run ends with a summary.
No exponentiation is written with a caret anywhere in this file.
"""
import itertools, time, sys
import flint, sympy as sp
from itertools import permutations
try:
    from sage.all__sagemath_singular import PolynomialRing, QQ, prod
except ImportError:
    from sage.all import PolynomialRing, QQ, prod

RESULTS = []
def report(name, ok):
    RESULTS.append((name, ok)); print(("PASS " if ok else "FAIL ") + name, flush=True)

# ---------- helpers ----------
def det4(M):
    tot = None
    for p in permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]: sg = -sg
        t = M[0][p[0]] * M[1][p[1]] * M[2][p[2]] * M[3][p[3]]
        tot = (t if sg > 0 else -t) if tot is None else (tot + t if sg > 0 else tot - t)
    return tot
def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
def adj3(M):
    c = lambda i, j: M[(i+1)%3][(j+1)%3]*M[(i+2)%3][(j+2)%3] - M[(i+1)%3][(j+2)%3]*M[(i+2)%3][(j+1)%3]
    return [[c(j, i) for j in range(3)] for i in range(3)]
def hess(f, idx): return [[f.derivative(i).derivative(j) for j in idx] for i in idx]
def monos(nv, d): return [e for e in itertools.product(range(d + 1), repeat=nv) if sum(e) == d]
def form(ctx, xs, coeffs, d):
    Z = ctx.from_dict({})
    for c, e in zip(coeffs, monos(len(xs), d)):
        t = c
        for x, k in zip(xs, e): t = t * x**int(k)
        Z = Z + t
    return Z
def coeff_eqs(D, nx, R):
    """split D (flint) by monomials in the first nx variables -> list of Sage polynomials in the rest"""
    G = R.gens(); bucket = {}
    for ex, co in D.to_dict().items(): bucket.setdefault(ex[:nx], []).append((ex[nx:], co))
    mk = lambda lst: sum(QQ(str(c)) * prod(G[i]**int(k) for i, k in enumerate(ex) if k) for ex, c in lst)
    return {k: mk(v) for k, v in bucket.items()}

t0 = time.time()
# ============ S1. ordinary-degree layers 11 and 10 (ternary cone top) ============
nm = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + ['g%d'%i for i in range(35)] + ['k%d'%i for i in range(20)] + ['q%d'%i for i in range(10)]
ctx = flint.fmpq_mpoly_ctx.get(nm, 'degrevlex'); g = ctx.gens(); x1, x2, x3, s = g[:4]
Pc = g[4:25]; F4 = g[25:60]; F3 = g[60:80]; F2 = g[80:90]
phi = form(ctx, [x1,x2,x3], Pc, 5)
f4 = form(ctx, [x1,x2,x3,s], F4, 4); f3 = form(ctx, [x1,x2,x3,s], F3, 3); f2 = form(ctx, [x1,x2,x3,s], F2, 2)
A = hess(phi, [0,1,2]); dA = det3(A)
def layer(D, d):
    return ctx.from_dict({ex: co for ex, co in D.to_dict().items() if sum(ex[:4]) == d})
D = det4(hess(phi + f4 + f3 + f2, [0,1,2,3]))
report("S1a [det Hess f]_12 = 0 for a cone top phi(x1,x2,x3)", layer(D, 12) == 0)
report("S1b [det Hess f]_11 = detHess3(phi) * d44 f4", layer(D, 11) == dA * f4.derivative(3).derivative(3))
# impose d44 f4 = 0 : f4 = g4(x') + s h(x')
nm2 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + ['g%d'%i for i in range(15)] + ['h%d'%i for i in range(10)] + ['k%d'%i for i in range(20)] + ['q%d'%i for i in range(10)]
ctx2 = flint.fmpq_mpoly_ctx.get(nm2, 'degrevlex'); g2 = ctx2.gens(); y1, y2, y3, u = g2[:4]
phi2 = form(ctx2, [y1,y2,y3], g2[4:25], 5); g4 = form(ctx2, [y1,y2,y3], g2[25:40], 4); h3 = form(ctx2, [y1,y2,y3], g2[40:50], 3)
f3b = form(ctx2, [y1,y2,y3,u], g2[50:70], 3); f2b = form(ctx2, [y1,y2,y3,u], g2[70:80], 2)
D2 = det4(hess(phi2 + g4 + u*h3 + f3b + f2b, [0,1,2,3]))
A2 = hess(phi2, [0,1,2]); dA2 = det3(A2); ad = adj3(A2); gr = [h3.derivative(i) for i in range(3)]
quad = ctx2.from_dict({})
for i in range(3):
    for j in range(3): quad = quad + gr[i]*ad[i][j]*gr[j]
lay10 = ctx2.from_dict({ex: co for ex, co in D2.to_dict().items() if sum(ex[:4]) == 10})
report("S1c [det]_10 = detA * d44 f3 - grad(h)**T adj(A) grad(h)", lay10 == dA2*f3b.derivative(3).derivative(3) - quad)
print("   (S1 %.0fs)" % (time.time()-t0), flush=True)

# ============ S2. weighted top (weights 1,1,1,2) ============
t1 = time.time()
nm3 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + ['h%d'%i for i in range(10)] + ['m%d'%i for i in range(3)]
low = [e for e in itertools.product(range(5), repeat=4) if e[0]+e[1]+e[2]+2*e[3] <= 4 and sum(e) >= 2]
nm3 += ['c%d'%i for i in range(len(low))]
ctx3 = flint.fmpq_mpoly_ctx.get(nm3, 'degrevlex'); g3 = ctx3.gens(); z1, z2, z3, w = g3[:4]
phi3 = form(ctx3, [z1,z2,z3], g3[4:25], 5); hh = form(ctx3, [z1,z2,z3], g3[25:35], 3); mm = form(ctx3, [z1,z2,z3], g3[35:38], 1)
lw = ctx3.from_dict({})
for c, e in zip(g3[38:], low): lw = lw + c*z1**e[0]*z2**e[1]*z3**e[2]*w**e[3]
psi = phi3 + w*hh + w**2*mm
D3 = det4(hess(psi + lw, [0,1,2,3]))
top = ctx3.from_dict({ex: co for ex, co in D3.to_dict().items() if ex[0]+ex[1]+ex[2]+2*ex[3] == 10})
report("S2  weight-10 part of det Hess f = det Hess psi, psi = phi + s h + s**2 m (general symbolic lower terms)", top == det4(hess(psi, [0,1,2,3])))
print("   (S2 %.0fs)" % (time.time()-t1), flush=True)

# ============ S3. branch m = 0 ============
t1 = time.time()
nm4 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + ['h%d'%i for i in range(10)]
ctx4 = flint.fmpq_mpoly_ctx.get(nm4, 'degrevlex'); g4_ = ctx4.gens(); a1, a2, a3, v = g4_[:4]
phi4 = form(ctx4, [a1,a2,a3], g4_[4:25], 5); h4 = form(ctx4, [a1,a2,a3], g4_[25:35], 3)
D4 = det4(hess(phi4 + v*h4, [0,1,2,3]))
sc = lambda D, k: D.to_dict() and D.__class__  # placeholder
def scoef(D, ctx_, k, sidx=3):
    return ctx_.from_dict({ex: co for ex, co in D.to_dict().items() if ex[sidx] == k})
report("S3a m=0: s**3, s**4 coefficients vanish", scoef(D4, ctx4, 3) == 0 and scoef(D4, ctx4, 4) == 0)
report("S3b m=0: s**2 coefficient = -(3/2) h detHess3(h)", 2*scoef(D4, ctx4, 2) == -3*h4*det3(hess(h4, [0,1,2]))*v**2)
# normal forms of binary cubics : x1**3, x1**2 x2, x1**3 + x2**3
def branch_phi(hexpr_fn, label):
    nm5 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)]
    c5 = flint.fmpq_mpoly_ctx.get(nm5, 'degrevlex'); gg = c5.gens(); b1, b2, b3, sv = gg[:4]
    ph = form(c5, [b1,b2,b3], gg[4:25], 5)
    Dd = det4(hess(ph + sv*hexpr_fn(b1,b2,b3), [0,1,2,3]))
    R = PolynomialRing(QQ, nm5[4:] + ['T'], order='degrevlex')
    eqs = list(coeff_eqs(Dd, 4, R).values())
    dAc = list(coeff_eqs(det3(hess(ph, [0,1,2])), 4, R).values())
    T = R.gens()[-1]
    ok = all(list(R.ideal(eqs + [1 - T*dd]).groebner_basis()) == [1] for dd in dAc)
    report("S3c m=0, h=%s : detHess3(phi) vanishes on every solution (branch excluded)" % label, ok)
branch_phi(lambda b1,b2,b3: b1**2*b2, "x1**2 x2")
branch_phi(lambda b1,b2,b3: b1**3 + b2**3, "x1**3+x2**3")
# h = x1**3 : det Hess psi = -9 x1**4 * detHess_(x2,x3)(phi)
nm6 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)]
c6 = flint.fmpq_mpoly_ctx.get(nm6, 'degrevlex'); gg6 = c6.gens(); d1, d2, d3, sv6 = gg6[:4]
ph6 = form(c6, [d1,d2,d3], gg6[4:25], 5)
H23 = ph6.derivative(1).derivative(1)*ph6.derivative(2).derivative(2) - ph6.derivative(1).derivative(2)**2
report("S3d m=0, h=x1**3 : det Hess psi = -9 x1**4 detHess_(x2,x3)(phi)", det4(hess(ph6 + sv6*d1**3, [0,1,2,3])) == -9*d1**4*H23)
print("   (S3 %.0fs)" % (time.time()-t1), flush=True)

# ============ S4. branch m = x1 ============
t1 = time.time()
nm7 = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + ['q%d'%i for i in range(10)]
c7 = flint.fmpq_mpoly_ctx.get(nm7, 'degrevlex'); gg7 = c7.gens(); e1, e2, e3, sv7 = gg7[:4]
ph7 = form(c7, [e1,e2,e3], gg7[4:25], 5); h7 = form(c7, [e1,e2,e3], gg7[25:35], 3)
D7 = det4(hess(ph7 + sv7*h7 + sv7**2*e1, [0,1,2,3]))
H23h = h7.derivative(1).derivative(1)*h7.derivative(2).derivative(2) - h7.derivative(1).derivative(2)**2
report("S4a m=x1 : s**4 coefficient = -4 detHess_(x2,x3)(h)", scoef(D7, c7, 4) == -4*H23h*sv7**4)
def m1_branch(label, hfn, extra_names, extra_eqs, mode):
    nmb = ['x1','x2','x3','s'] + ['p%d'%i for i in range(21)] + extra_names
    cb = flint.fmpq_mpoly_ctx.get(nmb, 'degrevlex'); gb = cb.gens(); f1, f2_, f3_, sb = gb[:4]
    E = dict(zip(extra_names, gb[25:]))
    ph = form(cb, [f1,f2_,f3_], gb[4:25], 5)
    Dd = det4(hess(ph + sb*hfn(f1,f2_,f3_,E) + sb**2*f1, [0,1,2,3]))
    R = PolynomialRing(QQ, nmb[4:] + ['T'], order='degrevlex'); G = R.gens(); T = G[-1]
    eqs = list(coeff_eqs(Dd, 4, R).values()) + [R(x) for x in extra_eqs]
    dAc = list(coeff_eqs(det3(hess(ph, [0,1,2])), 4, R).values())
    if mode == 'exclude':
        ok = all(list(R.ideal(eqs + [1 - T*dd]).groebner_basis()) == [1] for dd in dAc)
        report("S4 %s : detHess3(phi) vanishes on every solution (branch excluded)" % label, ok)
    else:
        R2 = PolynomialRing(QQ, nmb[4:], order='degrevlex')
        Pr = R2.ideal([R2(str(x)) for x in eqs]).minimal_associated_primes()
        dA2 = [R2(str(x)) for x in dAc]
        surv = [Q for Q in Pr if not all(d in Q for d in dA2)]
        expected = R2.ideal([R2('4*p17-1')] + [R2('p%d' % i) for i in range(17)])
        report("S4 %s : the only component with detHess3(phi) != 0 is phi = x1**3 x2**2/4 + p18 x1**4 x3 + p19 x1**4 x2 + p20 x1**5" % label,
               len(surv) == 1 and surv[0] == expected)
hn = lambda a,b,c,E: E['a']*b**3 + E['b']*a*b**2 + a**2*(E['c']*b + E['cp']*c) + E['e']*a**3
m1_branch("m=x1, a=1 (I)", hn, ['a','b','c','cp','e'], ['a-1'], 'exclude')
m1_branch("m=x1, a=0, b=1, c=0 (II)", hn, ['a','b','c','cp','e'], ['a','b-1','c'], 'exclude')
m1_branch("m=x1, h=x1**2 x2 (III.1)", lambda a,b,c,E: a**2*b, [], [], 'classify')
m1_branch("m=x1, h=e x1**3 (III.2)", lambda a,b,c,E: E['e']*a**3, ['e'], [], 'exclude')
print("   (S4 %.0fs)" % (time.time()-t1), flush=True)

# ============ S5. top (F): no completion with nonzero constant Hessian ============
t1 = time.time()
lowF = [e for e in itertools.product(range(5), repeat=4) if e[0]+e[1]+e[2]+2*e[3] <= 4 and sum(e) >= 2]
nmF = ['x1','x2','x3','s'] + ['c%d' % i for i in range(len(lowF))]
cF = flint.fmpq_mpoly_ctx.get(nmF, 'degrevlex'); gF = cF.gens(); w1, w2, w3, ws = gF[:4]
fF = w1*(ws + w1*w2)**2 + w1**4*w3
for c, e in zip(gF[4:], lowF): fF = fF + c*w1**e[0]*w2**e[1]*w3**e[2]*ws**e[3]
DF = det4(hess(fF, [0,1,2,3]))
RF = PolynomialRing(QQ, nmF[4:] + ['T'], order='degrevlex'); TF = RF.gens()[-1]
eqF = coeff_eqs(DF, 4, RF); constF = eqF.pop((0,0,0,0))
report("S5  top (F): the constant term of det Hess f vanishes on every solution", list(RF.ideal(list(eqF.values()) + [1 - TF*constF]).groebner_basis()) == [1])
print("   (S5 %.0fs)" % (time.time()-t1), flush=True)

# ============ S6. top (D), a = 0 : the four links, and Lemma A ============
t1 = time.time()
X1, X2, X3, U, a_ = sp.symbols('x1 x2 x3 u a'); XX = [X1, X2, X3, U]
Pf = sp.Function('P')(X1, X2, X3); Qf = sp.Function('Q')(X1, X2, X3)
Hs = lambda F: sp.Matrix(4, 4, lambda i, j: sp.diff(F, XX[i], XX[j]))
z = sp.Symbol('z'); Rr = Pf - Qf**2/(4*a_); Sz = Rr + z*Qf
HSz = sp.Matrix(3, 3, lambda i, j: sp.diff(Sz, XX[i], XX[j]))
report("S6a Lemma A (any Q): det Hess(P+uQ+a u**2) = 2a det Hess_x(R + zQ) at z = u + Q/(2a)",
       sp.simplify(Hs(Pf + U*Qf + a_*U**2).det() - 2*a_*HSz.det().subs(z, U + Qf/(2*a_))) == 0)
r = sp.symbols('r0:9')
Rc = X1**3 + r[0]*X1**2 + r[1]*X2**2 + r[2]*X3**2 + r[3]*X1*X2 + r[4]*X1*X3 + r[5]*X2*X3 + r[6]*X1 + r[7]*X2 + r[8]*X3
c2 = sp.Poly(sp.expand(Hs(Pf + U*Rc).det()).coeff(U, 2), X1, X2, X3).coeffs()
S9 = PolynomialRing(QQ, ['r%d' % i for i in range(9)] + ['T'], order='degrevlex'); T9 = S9.gens()[-1]
I9 = [S9(str(sp.expand(c))) for c in c2]
Mx = sp.Matrix([[2*r[1], r[5], r[3], r[7]], [r[5], 2*r[2], r[4], r[8]]])
mins = [sp.expand(Mx[:, [i, j]].det()) for i in range(4) for j in range(i+1, 4)]
report("S6b u**2 coefficient = 0  =>  (dR/dx2, dR/dx3) proportional (R depends on x1 and one form in x2,x3)",
       all(list(S9.ideal(I9 + [1 - T9*S9(str(m))]).groebner_basis()) == [1] for m in mins if m != 0))
Rt = sp.Function('Rt')(X1, X2)
br = sp.diff(Rt,X1)**2*sp.diff(Rt,X2,2) - 2*sp.diff(Rt,X1)*sp.diff(Rt,X2)*sp.diff(Rt,X1,X2) + sp.diff(Rt,X2)**2*sp.diff(Rt,X1,2)
report("S6c u-coefficient = -P_33 * bracket(R) for R = R(x1,x2)", sp.simplify(sp.expand(Hs(Pf + U*Rt).det()).coeff(U, 1) + sp.diff(Pf, X3, 2)*br) == 0)
w_, r0_, aa, bb, cc, dd = sp.symbols('w r0 aa bb cc dd')
Rq = X1**3 + r0_*X1**2 + aa*w_**2 + bb*X1*w_ + cc*X1 + dd*w_
brq = sp.Poly(sp.expand(sp.diff(Rq,X1)**2*sp.diff(Rq,w_,2) - 2*sp.diff(Rq,X1)*sp.diff(Rq,w_)*sp.diff(Rq,X1,w_) + sp.diff(Rq,w_)**2*sp.diff(Rq,X1,2)), X1, w_).coeffs()
T_ = sp.Symbol('T_')
report("S6d bracket = 0 for R = x1**3 + (deg<=2 in x1,w)  =>  R = g(x1)",
       all(sp.groebner(brq + [1 - T_*vv], aa, bb, dd, r0_, cc, T_).exprs == [1] for vv in (aa, bb, dd)))
gfun = sp.Function('g')(X1)
report("S6e R = g(x1): det Hess = -g'**2 detHess_(x2,x3)(P)",
       sp.simplify(Hs(Pf + U*gfun).det() + sp.diff(gfun, X1)**2*(sp.diff(Pf,X2,2)*sp.diff(Pf,X3,2) - sp.diff(Pf,X2,X3)**2)) == 0)
P0 = sp.Function('P0')(X1, X2); P1 = sp.Function('P1')(X1, X2)
J = sp.diff(P1,X1)*sp.diff(Rt,X2) - sp.diff(P1,X2)*sp.diff(Rt,X1)
report("S6f P = P0(x1,x2) + x3 P1(x1,x2), R = R(x1,x2): det Hess = J(P1,R)**2", sp.simplify(Hs(P0 + X3*P1 + U*Rt).det() - J**2) == 0)
print("   (S6 %.0fs)" % (time.time()-t1), flush=True)

# ============ S7. Ni's Proposition 6.1, branch a = 0 identities ============
t1 = time.time()
al = sp.Function('al')(X1, X2); be = sp.Function('be')(X1, X2)
report("S7a rank 2: det Hess(al + x3 be + u x1 x2) = (x1 be_1 - x2 be_2)**2",
       sp.simplify(Hs(al + X3*be + U*X1*X2).det() - (X1*sp.diff(be,X1) - X2*sp.diff(be,X2))**2) == 0)
mu = sp.Symbol('mu'); al3 = sp.Function('al')(X1, X3); be3 = sp.Function('be')(X1, X3)
report("S7b rank 1: det Hess(al + x2 be + u(x1 + mu x3**2/2)) = (be_3 - mu x3 be_1)**2",
       sp.simplify(Hs(al3 + X2*be3 + U*(X1 + mu*X3**2/2)).det() - (sp.diff(be3,X3) - mu*X3*sp.diff(be3,X1))**2) == 0)
report("S7c rank 0: det Hess(P + u x3) = -detHess_(x1,x2)(P)",
       sp.simplify(Hs(Pf + U*X3).det() + (sp.diff(Pf,X1,2)*sp.diff(Pf,X2,2) - sp.diff(Pf,X1,X2)**2)) == 0)
print("   (S7 %.0fs)" % (time.time()-t1), flush=True)


# ============ S8. complements requested by the independent review ============
t1 = time.time()
# S8a : Rabinowitsch confirmation of III.1 (no primary decomposition): every family-violating
# coordinate g satisfies g*d in rad(I) for every coefficient d of det A.
nmb = ['x1','x2','x3','s'] + ['p%d' % i for i in range(21)]
cb = flint.fmpq_mpoly_ctx.get(nmb, 'degrevlex'); gb = cb.gens(); f1, f2_, f3_, sb = gb[:4]
ph = form(cb, [f1,f2_,f3_], gb[4:25], 5)
Dd = det4(hess(ph + sb*f1**2*f2_ + sb**2*f1, [0,1,2,3]))
R8 = PolynomialRing(QQ, nmb[4:] + ['T'], order='degrevlex'); G8 = R8.gens(); T8 = G8[-1]
eqs8 = list(coeff_eqs(Dd, 4, R8).values())
dA8 = list(coeff_eqs(det3(hess(ph, [0,1,2])), 4, R8).values())
viol = [G8[i] for i in range(17)] + [4*G8[17] - 1]
ok = all(list(R8.ideal(eqs8 + [1 - T8*gv*dd]).groebner_basis()) == [1] for gv in viol for dd in dA8)
report("S8a III.1 by Rabinowitsch: off the family phi = x1**3 x2**2/4 + p18 x1**4 x3 + p19 x1**4 x2 + p20 x1**5, det A = 0", ok)
# S8b : det Hess3 of the family, and of the (D) top
p18, p19, p20, al_ = sp.symbols('p18 p19 p20 alpha')
phF = X1**3*X2**2/4 + p18*X1**4*X3 + p19*X1**4*X2 + p20*X1**5
H3 = lambda F: sp.Matrix(3, 3, lambda i, j: sp.diff(F, [X1,X2,X3][i], [X1,X2,X3][j]))
report("S8b det Hess3 of the III.1 family = -8 p18**2 x1**9", sp.expand(H3(phF).det() + 8*p18**2*X1**9) == 0)
Qb = sp.Function('Qb')(X1, X2)
report("S8c det Hess3(Q(x1,x2) + alpha x1**4 x3) = -16 alpha**2 x1**6 Q_22",
       sp.simplify(H3(Qb + al_*X1**4*X3).det() + 16*al_**2*X1**6*sp.diff(Qb, X2, 2)) == 0)
# S8d : normalisation of (F)
psiF = X1*U**2 + U*X1**2*X2 + phF
sub = {X3: (X3 - p19*X2 - p20*X1)/p18}
t_ = sp.expand(psiF.subs(sub, simultaneous=True).subs(X2, 2*X2))
report("S8d (F) normalisation: x3 <- (x3 - p19 x2 - p20 x1)/p18, x2 -> 2 x2 gives x1 (u + x1 x2)**2 + x1**4 x3",
       sp.expand(t_ - (X1*(U + X1*X2)**2 + X1**4*X3)) == 0)
# S8e : step (iv): for v = (0, v2, v3), d_v**2 of Q(x1,x2) + x1**4 x3 equals v2**2 Q_22
v2, v3 = sp.symbols('v2 v3')
top = Qb + X1**4*X3
dv = lambda F: v2*sp.diff(F, X2) + v3*sp.diff(F, X3)
report("S8e d_v**2 (Q(x1,x2) + x1**4 x3) = v2**2 Q_22 for v = (0,v2,v3)", sp.simplify(dv(dv(top)) - v2**2*sp.diff(Qb, X2, 2)) == 0)
print("   (S8 %.0fs)" % (time.time()-t1), flush=True)

nf = sum(1 for _, ok in RESULTS if not ok)
print("=" * 60)
print("SUMMARY: %d checks, %d failed, total %.0fs" % (len(RESULTS), nf, time.time() - t0))
print("ALL CHECKS PASSED" if nf == 0 else "SOME CHECKS FAILED")
