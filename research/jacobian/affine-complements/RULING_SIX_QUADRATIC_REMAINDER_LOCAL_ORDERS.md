# The quadratic remainder vanishes to half the cubic leading order

PF2RULED174. 17 November 2026. Actual complex polynomial plane Keller maps.

For every off-M h-root whose outer face contains four escaping
sheets, the coefficient d of z in PF171's quadratic remainder
satisfies ord_c d>=ord_c(h)/2. In particular, at every actual
A pole divisor of slope L,

    ord_c d>=L.

If d is nonzero, it is divisible by the product of these pole
divisors to their slopes, so deg d is at least their total slope.
The linear-remainder branch d=0 remains allowed by this theorem.

The argument keeps the entire cubic square and the full quartic
factor; no lower-ruling function is declared a target coordinate.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
JC2 remains open.

## 1. Whole scaled polynomial and actual input

Use the regular normal form of
[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md):

    w^2=N+M*z, u=t^2, dP wedge dQ=du wedge dw/M,
    P=H3^2+R2,
    H3=(h*w+K)*z+B0*w+e0,
    R2=d*z+e*w+f,                                       (1)

where all displayed coefficients belong to C[u], and H3,R2
remain outside the actual target field. The full original
normalizations S_E and S, and every one of their primes, remain
unchanged. Coefficient and companion degrees are unrestricted.

Let c be off M, v=u-c, and suppose the outer face consumes all
four escaping sheets. PF168 gives rational ell>0 with

    j=ord_c h=2ell,
    in(P*)=(aY^3+bY^2+c1Y)^2, a*c1!=0, Y=v^ell*w.       (2)

This includes one quartic root and two distinct double roots.
The square-root sign agrees with h/M.

Set S=v^ell*H3*(v,v^-ell*Y). Its leading coefficients are
regular. Its constant term is

    v^ell*e0-v^ell*K*N/M,

which has order at least ell. Therefore

    S=Y*A(v,Y)+O(v^ell),                                 (3)

with A a regular quadratic of unit leading coefficient a2(v),
and A(0,0)=c1!=0.

Assume d!=0 and put aord=ord_c d. We prove aord>=ell.
The initial square in (2) already forces d(c)=0: the initial
of S is the chosen cubic square root, and a unit d would add
a nonzero quadratic initial. Thus aord>0.

Suppose, for contradiction,

    0<aord<ell.                                         (4)

The WHOLE scaled equation is

    F(v,Y,p)=v^(2ell)*(P*(v,v^-ell*Y)-p)
      =S^2+(d/M)*Y^2+v^ell*e*Y
                          +v^(2ell)*(f-d*N/M-p).         (5)

Every coefficient is regular after clearing fractional powers
with v=s^alpha. The only parameter term is -p*v^(2ell).

Modulo v^ell, equations (3),(5) give
Y^2*(A^2+d/M). Normalize the unit leading coefficient a2^2
and Hensel-factor the finite pair from the four nonzero roots:

    F/a2^2=F2*G,
    F2 monicdegree2, G monicdegree4,
    F2(0,Y)=Y^2,
    G=(A/a2)^2+d/(M*a2^2) modulo v^ell.                  (6)

The factors are coprime at v=0 because A(0,0)!=0. The factorization
is over C[p][[s]], with constant-unit resultants, and retains
both possible double clusters.

## 2. Depressed quartic coefficients detect the small perturbation

Translate the monic quartic to remove its cubic term. Its EXACT
form is

    G=X^4+2B*X^2+q*X+B^2+Delta,                          (7)

where B is regular, ord q>=ell, and

    ord Delta=aord<ell.                                 (8)

Indeed a square of a monic quadratic becomes (X^2+B)^2,
and adding a constant adds exactly that constant to Delta.
Equation (6) preserves these identities modulo v^ell.
In particular the leading coefficient of Delta is nonzero and
independent of p. All the low-order coefficients used below
are independent of p, because p first enters the whole equation
at order 2ell and constant-unit factorization preserves that order.

We show that at least one of the four roots has

    mu=ord_v G'(X)<ell.                                 (9)

Put bval=ord B, allowing infinity. There are three complete cases.

If 2bval<aord, take a p-independent centre
X0=+sqrt(-B) or -sqrt(-B), to the finite precision required
for the following leading terms. Its order is bval/2.
The perturbation Delta splits that quadratic centre at displacement
order (aord-bval)/2. The term qX is higher because ord q>=ell>aord.
The reduced equation has a nonzero constant and a nonzero
quadratic coefficient, with simple roots. At either root,

    mu=(aord+bval)/2<3aord/4<ell.                         (10)

If 2bval>aord, including B=0, put X=v^(aord/4)*Z.
The full initial equation is Z^4+Delta_lead=0, with four
nonzero simple roots. The B and q terms are strictly higher.
At those roots

    mu=3aord/4<ell.                                      (11)

If 2bval=aord, the same substitution gives

    (Z^2+B_lead)^2+Delta_lead=0.                         (12)

Since Delta_lead!=0, this quartic has at least two nonzero
simple roots. The only possible multiple root is Z=0 when
B_lead^2+Delta_lead=0; the other pair remains simple.
At a simple nonzero root equation (11) holds again.

All parameter-independent leading roots, all Puiseux denominators
and the possible multiple initial root in (12) are retained.
No genericity hypothesis on B or Delta is imposed.

## 3. The small derivative forces an affine response and a graph

Choose the root in (9). Its original Y has nonzero constant term,
so the finite quadratic factor F2 is a unit. On that root,

    ord_v P*_w=mu-ell<0,
    ord_v Y_p=2ell-mu>ell.                              (13)

The coframe is holomorphic, with exact s order
alpha*(ell-mu)+alpha-1 after v=s^alpha.

The whole polynomial in (5), not just its quartic factor, is
F0(v,Y)-p*v^(2ell), with F0 independent of p and regular.
Let r=ord_v Y_p>ell and truncate Y just before its first
coefficient lambda outside C. The extended partial_p derivation
has constants C, so r is that exact first-coefficient index.

Every Taylor term of degree at least two at this p-independent
truncation has order at least 2r>2ell. The derivative changes
only in order at least r, strictly above mu, so the linear term
has exact order mu+r=2ell and nonzero constant leading slope.
The constant term cannot occur earlier. The initial equation is

    A0*lambda+B0const=p, A0 in C*.                       (14)

Thus the response is affine. Constant-unit Hensel recursion gives
all later coefficients in C[p] after clearing Puiseux powers.
Equivalently the PF159 full Taylor-leading lemma makes the
finite companion residue a polynomial in lambda, hence in p.

Since Y has nonzero limit, w and the original source x are polar;
the original t lifts over c are units. The original full finite
normalization therefore contains an omitted prime with WHOLE
closed graph image Q=q0(P). The retained Chau component theorem,
in its exact PF165 application, excludes that component.

This contradiction proves

    ord_c d>=ell=ord_c(h)/2.                             (15)

The d=0 case satisfies the valuation inequality trivially.

## 4. All actual pole divisors and global coefficient consequence

At every surviving A pole divisor ell=L is an integer,
ord h=2L and (15) gives ord d>=L. Index c_i by the DISTINCT
supporting base values u=c_i
whose quartic cluster carries at least one actual pole, with
slope L_i. Several normalized pole primes in the same cluster
contribute only ONCE to this product. Put

    A_pole(u)=product_i (u-c_i)^(L_i).

Then

    A_pole divides d; if d!=0, deg d>=sum_i L_i.        (16)

Other h-roots are not silently included in this divisor.
This theorem supplies a necessary zero order, not global
h-divisibility of d. It does not exclude the linear remainder
d=0,e!=0 which remains after PF172.

The checker verifies the full Hensel factor precision, the exact
depressed-quartic constant invariant, all three Newton cases,
the multiple-initial-root boundary, and exact response/coframe
orders. The checker passes 65 exact controls. The uniform proof is (5)--(15), not finite
coefficient testing.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_quadratic_remainder_orders.py

No new external classification is imported. All original primes
and source/target distinctions remain retained. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
Derivative integrality, automatic trace perfection, normal-tail termination
and HC4/JC2 remain open. Independent
whole-proof review, formal verification and external human
review remain unclaimed.
