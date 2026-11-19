# A rational residue hierarchy constrains every odd half-power

The subsequent [constant companion resonance theorem](RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md)
excludes the entire linear-remainder branch, including the finite-test
controls below. The [balanced quadratic theorem](RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md)
requires both d and e to be nonzero and their weights to agree at every
actual A support. The residue hierarchy proved here remains necessary;
its individual exactness conditions are not sufficient for a companion.

PF2RULED175. 17 November 2026. Actual complex polynomial plane Keller maps.

Every positive odd half-power of P has a formal residue whose quotient
by M has a rational primitive in u. This supplies an unrestricted
hierarchy of necessary equations. The first condition reproduces the
quadratic-coefficient equation of PF173; the next gives a new condition
for a linear remainder:

    R2=e(u)w+f(u) => e(u)^2*du/h(u) is rationally exact.

A separate actual-pole argument forces the quadratic remainder to
have a pole, or a nonconstant residue, at every A pole cluster.
It gives ord d<2L or ord e<L at each such supporting base value.
In the linear branch, a unit-K h-root forces ord e>=ord h;
the central Keller condition then forces h(0)!=0.

The equations and coefficient controls below are necessary conditions
on the A class. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes that class in the retained marked quadratic
model. JC2 remains open.

## 1. Actual inputs and a formal coefficient residue

Use the original fields, full normalizations and regular normal form
of PF171:

    E=C(u,w), F=C(P,Q), L=C(x,t), u=t^2,
    w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    P=H3^2+R2,
    H3 in O(V), deg_w H3=3, leading_w H3=h/M,
    R2=d(u)z+e(u)w+f(u), d,e,f in C[u].                  (1)

The coefficient and companion degrees are unrestricted. Neither H3
nor R2 is assumed to belong to F or the full finite target
normalization. The trace/perfection obstruction is not assumed.

Let K=C(u). In K((w^-1)) choose the unique square root

    Z=P^(1/2), leading_w Z=(h/M)*w^3.

Its polynomial part is H3. Define

    res_w A=[w^-1]A,                                     (2)

the coefficient of w^-1 in a formal Laurent series. This is a
coefficient residue; the geometric residue at infinity has the
opposite sign. All identities below use (2).

For every positive odd integer n put

    m_n=res_w P^(n/2),
    J_n=res_w (Q*P^(n/2)*P_w).                           (3)

Both are rational functions of u. Although the Laurent expansions
can have infinitely many negative terms, they have finite upper
degree, and every required coefficient is obtained from finitely
many algebraic operations in K.

## 2. The entire actual Jacobian gives an exact hierarchy

Formal differentiation commutes with res_w. Moreover the residue
of a w derivative is zero. Therefore

    dJ_n/du
      =res_w[Q_u*P^(n/2)*P_w
          +(n/2)*Q*P^(n/2-1)*P_u*P_w
          +Q*P^(n/2)*P_uw].

Subtract the zero residue of d_w(Q*P^(n/2)*P_u).
The middle and last terms cancel, giving the exact identity

    dJ_n/du
      =res_w[P^(n/2)*(Q_u*P_w-Q_w*P_u)]
      =-m_n/M.                                         (4)

Only the ACTUAL Jacobian P_u*Q_w-P_w*Q_u=1/M was used.

Consequently, for EVERY positive odd n,

    (m_n/M)*du has a rational primitive.                 (5)

Equivalently all its finite residues vanish. This equivalence follows
from partial fractions in C(u): polynomial terms and poles of order
at least two integrate rationally, whereas a nonzero simple-pole
coefficient does not. The residue at infinity then vanishes as well.

The primitive in (4) is explicit as a coefficient of the actual Q.
It is not asserted to be polynomial, a target function, or integral
over B. These coefficient residues are not traces in a substitute
normalization or in a normal-tail quotient.

For direct computation, the binomial expansion P=H3^2+R2 gives
the finite formula

    m_n=sum_(j=(n+1)/2)^floor((3n+1)/4)
             binom(n/2,j)*res_w(H3^(n-2j)*R2^j).          (6)

Terms with n-2j>=0 are polynomials and have zero residue.
For the other terms, w degree is at most 3n-4j; once this
is below -1 there is again no contribution. This proves the
complete finite range in (6) for every n, rather than a prefix
of the hierarchy. Even positive n give polynomial P^(n/2)
and zero m_n.

No converse to (5) is asserted. Passing finitely or infinitely
many residue equations does not itself supply a polynomial mate.

## 3. The first and third characters

The first identity is

    m_1=d/(2h).                                         (7)

Indeed the first negative term in
sqrt(P)=H3+R2/(2H3)-... has w^-1 coefficient d/(2h);
all later binomial terms have lower w degree.
Thus (5) for n=1 independently recovers the rational-character
condition established in
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md):

    d*du/(M*h) is rationally exact.                       (8)

For n=3, only the R2^2/H3 term contributes:

    m_3=(3/8)*res_w(R2^2/H3).                            (9)

For an explicit general quadratic formula write

    H3=a3*w^3+a2*w^2+a1*w+a0,
    R2=r2*w^2+r1*w+r0.

The coefficient in (9) is

    Xi=(r1-r2*a2/a3)^2/a3
                   +2r2*r0/a3-r2^2*a1/a3^2.             (10)

This follows by expanding 1/H3 through the two terms needed
after its leading w^-3. Therefore Xi*du/M must be rationally exact.

If d=0, so R2=e*w+f, equation (10) reduces to Xi=e^2/a3.
Since a3=h/M, this gives the new complete linear-remainder condition

    e^2*du/h is rationally exact.                        (11)

It is nonvacuous: PF172 excludes e=0 in this branch, but does
not make (11) automatic.

At every finite root c of h, let j=ord_c h and b=ord_c e.
If j=2b+1, (11) has a simple pole with nonzero residue.
Hence this WHOLE valuation stratum is excluded:

    d=0,e!=0 => ord_c h != 2ord_c e+1.                   (12)

In particular, every simple zero of h must be a zero of e
in an actual linear-remainder survivor. Condition (12) does
not replace the full residue equation at higher poles.

The whole-plane Keller condition adds a central restriction in the
linear branch. On t=0, choose x0 with H3(x0,0)=0, using its
nonconstant affine restriction from PF171. At that point
w=0, du=0 and dw=Sigma(0)*dt. Thus

    dP=e(0)*Sigma(0)*dt,

so e(0)!=0 is necessary. With (12) this already excludes
ord_0 h=1. A second whole-plane argument excludes EVERY
positive central h order in this branch.

More generally let c be any h-root where K, the z coefficient
in H3=(h*w+K)*z+B0*w+e0, is a unit. Put j=ord_c h>=1
and m=ord_c M>=0. There is an exact p-independent Laurent
zero wo of H3 with

    wo=-K/h+O((u-c)^m), ord wo=-j,
    ord H3_w(wo)=-j-m.

To prove it, put U=hw+K and multiply H3=0 by h^2*M:

    U*((U-K)^2-N*h^2)+B0*h*M*(U-K)+e0*h^2*M=0.

Its U derivative at U=u-c=0 is K(c)^2, a unit. Constant-unit
recursion gives U of order at least j+m and the displayed wo.

If a=ord_c e<j, the linear remainder at wo has exact negative
order a-j; f is regular and cannot cancel it. Substitute

    w=wo+(u-c)^((a+j+2m)/2)*Y.

The displacement is smaller than wo. The entire P*-p equation,
after scaling and clearing a possible denominator two, has
a nonzero quadratic plus a nonzero constant as initial.
Both roots are simple and constant. The exact orders are

    ord_v w_p=(3j+2m-a)/2,
    omega=(u-c)^((3j-a)/2)*(unit)*du.

Here omega=dQ on P=p. Its normalized order after clearing
v=s^alpha is alpha*(3j-a)/2+alpha-1>=0.
Thus the companion has polynomial p residue. The original
source escapes (also at c=0 after u=t^2), and the full closed
graph is forbidden by the retained Chau component theorem,
as in PF165. Consequently

    K(c)!=0, h(c)=0, d=0 => ord_c e>=ord_c h.             (12a)

The central K is a unit by PF162/PF171. If h(0)=0, (12a)
would force e(0)=0, contradicting the central Keller condition.
The ENTIRE linear branch therefore satisfies

    d=0 => e(0)!=0 and h(0)!=0.                          (12b)

This does not make h divide e at common zeros of h and K.

## 4. A supporting pole cannot have constant remainder residue

Fix a DISTINCT supporting base value c of an actual A pole cluster,
of slope L. Write v=u-c. All branches in its quartic cluster have

    w=v^-L*(c0+higher terms), c0 in C*, M(c)!=0.           (13)

At an actual pole, the first coefficient lambda of w outside C,
after a finite C-coefficient truncation, has AFFINE parameter
response H0(lambda)=p by the retained PF159 polynomial-mate lemma.

Suppose both

    ord_c d>=2L, ord_c e>=L,                             (14)

using infinite order for a zero polynomial. Then R2 is finite
on this branch and its residue belongs to C: every nonzero
order-zero contribution in d*z or e*w uses only the constant
leading coefficient c0 in (13); f contributes f(c).
No later parameter coefficient contributes at order zero.

Since P=p and R2 is finite, H3 is finite as well.
The full Taylor-leading lemma for a polynomial source function
makes its residue J(lambda), J in C[Y]. It cannot be constant,
because

    p=J(lambda)^2+constant.                              (15)

Equation (15) has degree at least two in lambda, contradicting
the actual affine response H0(lambda)=p. Equivalently lambda
is affine in p, so no nonconstant polynomial in p can have
square equal to p minus a constant.

Thus every supporting A base value satisfies

    ord_c d<2L OR ord_c e<L.                             (16)

This is a statement about the actual pole branches, not every
escaping branch in every specialization. All original primes
remain retained.

Combining (16) with PF174's ord_c d>=L and the simple-pole
exclusion ord_c d!=2L-1 from (8) gives, in particular,

    ord_c e>=L => d!=0 and L<=ord_c d<=2L-2.              (17)

In the linear branch d=0, equation (16) forces ord_c e<L.
No entire linear-remainder closure follows from this bound.

## 5. Independent controls at the current smallest leading degree

Use the positive leading data retained in PF173:

    M=u^2, h=u*(u-1)^10, G=M^2*h=u^5*(u-1)^10,
    D=deg(M*h^2)=24, one candidate A slope L=5.

The displayed polynomial rho in PF173 solves
rho'-(G'/(3G))*rho=1. Thus these leading coefficients pass
the leading primitive equation. Set d=0, so (8) is vacuous.

For e_bad=(u-1)^5, the next condition FAILS:

    e_bad^2/h=1/u.

Its residue at zero is one. These coefficients are excluded
for an actual companion, despite passing the first rational
character and the leading coefficient ODE.

For e_good=u*(u-1)^5,

    e_good^2/h=u

has primitive u^2/2. This is a positive control for (11).
It fails the separate actual-pole condition (16), since its
order at u=1 is L. It is not an actual Keller realization.

There are also positive controls for BOTH (11) and (16).
Put v=u-1 and

    e_*=u*(u-1)^4*(3-u).

Then ord_(u=1)e_*=4<L=5, and

    e_*^2/h=v-3+4/v^2
            =d/du(v^2/2-3v-4/v).                        (18)

It also vanishes at the central simple h-root. These data pass the stated leading, first/third-character and
A-pole-support tests. Both positive D=24 e controls vanish at
the central origin and fail the separate central Keller condition
e(0)!=0. In fact the entire linear branch for these specific
h=u*(u-1)^10 leading data is excluded by the central condition
and (11), whatever e is.

There remain coefficient controls passing ALL conditions just
listed. Take

    M=u^2, h=(u-1)^14, e=(u-1)^6, d=0, D=30, L=7.

Here e(0)=1, ord_(u=1)e=6<L, and e^2/h=1/(u-1)^2,
with primitive -1/(u-1). Here h(0)=1, so the stronger
central condition (12b) also passes. For G=M^2*h=u^4*(u-1)^14,
the exact polynomial

    rho=3u(u-1)*(243u^4-891u^3+1188u^2-660u+110)/110

solves rho'-(G'/(3G))*rho=1. Thus the leading ODE,
first/third rational characters, A support and central tests
all pass.

A paired D=30 BAD control passes the other listed guards but
fails only the displayed third character: replace e by

    e_bad30=(u-1)^6*(2u-1).

Its central value is -1, its A order is still six, but

    e_bad30^2/h=1/v^2+4/v+4

has residue four at v=0. Thus the new residue condition remains
a proper obstruction even after the stronger central and
A-support conditions.

These are coefficient controls, not actual companions or proofs
of global submersion, full end coverage or completion of the
hierarchy.

## 6. Verification and the remaining goal

The checker verifies the complete formal residue identity on
actual polynomial Jacobian controls, the finite binomial formula
for several odd n, both generic coefficients (7),(10), all
displayed bad/good rational functions, the full unit-K cubic
splittings and the local constant-residue obstruction. The checker
passes 110 exact controls. The actual polynomial Jacobian controls
outside the cubic class verify only the general residue identity.
The unbounded theorem is (4)--(6) and
the actual polynomial-mate argument (13)--(16).

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_residue_hierarchy.py

No new external classification is imported. Source functions
remain separate from target functions; no all-prime integrality,
codifferent multiplication or polynomial termination is inferred.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
Automatic trace perfection, normal-tail termination and HC4/JC2 remain
open. Independent whole-proof
review, formal verification and external human review remain unclaimed.
