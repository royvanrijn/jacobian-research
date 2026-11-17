# Specialization of affine planes and Frobenius shear obstructions

Claims HC4MYKK0/1/2. [MATH_STATUS.json](../../MATH_STATUS.json) is the
authority. The previous [divided-equation obstruction](RAMIFICATION_INDEPENDENT_LIFTING_OBSTRUCTIONS_FOR_HC4.md)
left arbitrary nonlinear source changes outside its scope. Here an entire
family of genuine nonlinear changes preserves the characteristic-two
Hessian and collision, and one member passes the displayed divided
cubic--quartic equation. Nevertheless none can lift to characteristic
zero. The new obstruction is global polynomial affineness in translation
directions, specialized through a saturated lattice. It uses the retained
HCTANG2 structural alternative, with that reproduction's existing assurance
limits. It does not prove HC4 or exclude every nonlinear recharting.

## 1. Joint affineness specializes: HC4MYKK0

For a polynomial h over a field, call an r-dimensional constant vector
space W jointly affine if

    h(z+t1*w1+...+tr*wr)

has TOTAL degree at most one in the translation variables t1,...,tr,
for a basis of W, as a polynomial identity also in z. This is affineness
on EVERY translate, not at one point. It is stronger than separate
affineness in a chosen basis. In positive characteristic one must use
these translation identities, or equivalently all relevant Hasse
coefficients; vanishing ordinary second derivatives is insufficient.

Let O be a DVR with fraction field K and residue field k. If h in O[z]
has a jointly affine r-plane W over K, then its reduction has one over k.
Indeed

    L=W intersect O^n

is a saturated rank-r submodule of O^n. The quotient is torsion-free,
hence free over the DVR; L is a direct summand and a basis reduces to r
independent vectors. Every translation coefficient of total t-degree at
least two vanishes in K[z] and therefore in O[z]. Reducing the identity
gives the required plane. A basis of W with poles is harmless: use this
saturated lattice, not an arbitrary rescaled basis that may lose rank.

If W exists only over the algebraic closure of K, it is defined over a
finite extension K'. Extend the valuation and apply the same argument to
the valuation DVR O' (or its completion). Its residue is an algebraic
extension of k. This proves the geometric specialization statement.
Alternatively the jointly affine r-planes form a closed incidence locus
in the Grassmannian: the finitely many translation coefficients vanish
on its tautological bundle. The lattice argument proves the precise
valuative fact needed here, without importing a separate properness
theorem. [Stacks 089R](https://stacks.math.columbia.edu/tag/089R) supplies
the standard Grassmannian/local-free-quotient convention.

Now take n=4, characteristic K=0, and det Hess(h) in O^*. The retained
[HCTANG2 alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md) says that
over an algebraically closed characteristic-zero field either grad(h) is
an automorphism or h has a constant jointly affine two-plane. Its stated
complex-field form transfers to any characteristic-zero coefficient field:
take a finitely generated field of coefficients, embed it in C, and use
the finite polynomial coefficient equations for a plane on a Grassmann
chart. A complex solution gives a solution over an algebraic closure of
the coefficient field, and hence a plane over a finite extension of K.
If the gradient is an automorphism over C, its unique formal inverse
descends to the coefficient field, so this branch transfers as well.

The automorphism branch also has good reduction. Translate the constant
gradient output to zero. Its linear Jacobian has determinant a unit in
O, so the unique formal inverse has coefficients in O, constructed
recursively using only that inverse linear matrix. Any polynomial inverse
over K or an algebraic extension must agree with it. It is consequently
a polynomial over O and reduces to an inverse of grad(h mod pi).

We obtain the all-degree reduction alternative:

    Every good reduction of a characteristic-zero four-variable
    constant-unit-Hessian polynomial has either an invertible gradient,
    or a geometric constant jointly affine two-plane.                 (1)

In particular, a residue potential with noninvertible gradient and NO
such plane cannot lift, over any ramification and with arbitrary higher
potential terms in the lift. This holds in any residue characteristic.
It excludes characteristic-zero lifts, not a stated uniform finite
precision. Different degrees at successive finite precisions do not
produce a polynomial characteristic-zero lift.

## 2. A nonlinear characteristic-two family: HC4MYKK1

Over any field k of characteristic two retain

    P=x+x^2*y+x^4+x^6*y^2,
    Q=y+x^5+x^6*y+x^7*y^2+x^8*y^3,
    H0=uP+vQ.

Let a be ANY polynomial in k[x,y], and b ANY polynomial in k[x,y,u]
of positive degree in u. There are no degree or coefficient bounds.
Set

    phi=(x,y,u+a(x,y)^2,v+b(x,y,u)^2),
    H_ab=H0 composed phi
        =(u+a^2)P+(v+b^2)Q.                                      (2)

This is a polynomial source automorphism. Its inverse is

    (x,y,U,V) |-> (x,y,U+a^2,V+b(x,y,U+a^2)^2).

All derivatives of the added squares vanish in characteristic two:
J(phi)=I and every component Hessian is zero. The exact chain rules give

    grad(H_ab)=grad(H0) composed phi,
    Hess(H_ab)=Hess(H0) composed phi,
    det Hess(H_ab)=1.                                              (3)

The noninjective gradient of H0 remains noninjective under this source
automorphism: take the preimages of its three retained collision points.
Thus these are genuine characteristic-two unit-Hessian examples, not
polynomials which have lost the required Hessian condition under an
arbitrary nonlinear change.

Nevertheless their complete set of constant affine LINE directions over
any extension field is exactly k*partial_v, including zero. Here a line
direction w means deg_t H_ab(z+tw)<=1 as a polynomial identity.

To prove this, the coefficient of the independent variable v in that
identity is Q(x+t*w_x,y+t*w_y). It must have degree at most one in t.
The highest homogeneous part of Q is the single monomial x^8*y^3.
If w_x or w_y is nonzero, its translation has a nonzero t coefficient of
degree at least three: its maximal t power is the sum of 8 and/or 3 over
the nonzero direction entries. This coefficient has the corresponding
remaining x/y monomial and cannot cancel with terms of lower total
degree of Q. Equivalently take the homogeneous total-degree-eleven
component in x,y,t. Therefore w_x=w_y=0, in every characteristic-two
field extension; no perfection or infinitude is required.

Now only u and v translate. Put d=deg_u b>=1 and let b_d(x,y) be its
nonzero leading coefficient. If w_u!=0, the coefficient of t^(2d) in
H_ab(z+tw) is

    b_d(x,y)^2 * w_u^(2d) * Q(x,y),                               (4)

which is nonzero. All other contributions are at most linear in u,v.
Thus w_u=0. Conversely every v direction is affine by (2). This proves
the exact line-direction statement. In particular H_ab has no jointly
affine two-plane, even geometrically. Notice that its ordinary second
directional derivative is zero for EVERY direction in characteristic
two; that weaker test would give the wrong answer.

Combining (1)--(4), NONE of the potentials (2) has a characteristic-zero
constant-unit-Hessian polynomial lift over any mixed (0,2) complete DVR.
The putative lift may have arbitrary degree, arbitrary noncotangent
corrections, and a cotangent plane defined only after a ramified field
extension with coefficients of negative valuation. Saturation still
forces a geometric residue plane, which (4) excludes. Affine source
rechartings, nonzero scalar multiplication and affine potential additions
are excluded as well. This does NOT classify all nonlinear rechartings
of H0: b of u-degree zero is outside this particular line-direction proof.

## 3. A shear passes the old divided equation: HC4MYKK2

Take a=x and b=u. Then

    Hdagger=(u+x^2)P+(v+u^2)Q.

Its quadratic part is q=xu+yv, its cubic part is

    C=x^3+y*u^2,

and its quartic part is x^2*y*u. In the KJ1 notation ell_x=ell_y=
ell_u=ell_v=0, A+B=1 and T=3 over Z. Hence

    A+B-T=-2=0 modulo2,
    [xy]det Hess(Hdagger)=-8=0 modulo8.                           (5)

This is a genuine escape from that ONE divided cubic--quartic equation.
It is not a certificate satisfying all integer-contracted equations,
nor proof of a full lift modulo eight. The original zero-cubic argument
cannot simply be reused after this nonlinear change. KK1 rules out every
characteristic-zero lift for a different, global reason.

There are even naive lower-precision lifts for the WHOLE family. Lift
a,b and P,Q coefficientwise to a mixed (0,2) DVR and use (2). The
integral triangular map phi still has determinant Jacobian one, and its
component Hessians are all divisible by two. The integral Hessian chain
rule writes

    Hess(H0 composed phi)=J(phi)^T*(Hess H0 composed phi)*J(phi)+2B,

with B symmetric. The first summand A is symmetric with even diagonal.
For every symmetric four-by-four A with even diagonal and every symmetric B,

    det(A+2B)=det(A) modulo4.                                    (6)

Indeed the linear correction is 2 tr(adj(A)B). Every principal
three-by-three minor of A is even, since an odd-dimensional alternating
matrix has zero determinant in characteristic two. Thus adj(A) has even
diagonal, and its off-diagonal contributions to the trace pair up;
the trace is even. All terms of correction order at least two already
contain four. This is an integer identity, not division in the residue
ring. As det Hess(H0)=Jac(P,Q)^2=1 modulo4 and det J(phi)=1, (6) gives

    det Hess(H_ab)=1 modulo4                                    (7)

for these naive representatives, for arbitrary a,b. They therefore lift
through precision 2e when e=val_pi(2). No sharp maximal finite precision
for Hdagger or for the general sheared family is asserted. The all-degree
characteristic-zero obstruction is the absence of a jointly affine plane,
not a new tower of increasingly large lifting calculations.

## 4. Scope and next useful gate

KK0 is a necessary good-reduction condition following from the retained
Tang alternative, not an independent proof of that classification.
KK1 closes all nonlinear square shears (2) with positive u-degree in b.
KK2 shows why passing a low divided equation or the whole mod4 Hessian
condition does not rescue that family. The old KJ proof/checker/receipts
are preserved. No arbitrary nonlinear orbit exclusion, full integer-ideal
sufficiency, characteristic-zero counterexample, native source alignment,
or full HC4/JC2 decision is claimed.

A different arithmetic candidate must now pass BOTH tests: its divided
integer Hessian equations and, if its gradient is noninvertible, the
existence of a geometric constant jointly affine two-plane. Merely
removing that plane by a Frobenius source shear makes the lifting route
impossible. A supplied residue plane is only necessary: the original
H0 has one and is still excluded by KJ1. The general JC2 inversion and
the native critical non-CI recognition questions remain open.

The primary Grassmannian conventions were checked at Stacks089R. The
recent [Adamus--Hajto paper](https://arxiv.org/html/2610.36736v1),
introduction through its Nambu/Legendre criteria, was re-screened from its retained source for an
alternative global inversion input. Its abstract and stated criteria do
not settle HC4/JC2; the doubling criterion introduces no independent
conditions. None of its theorems is imported here. All new specialization,
translation and parity arguments above are proved directly. The only
classification input is the already-retained HCTANG2 theorem, with its
existing assurance limits.

The [checker](verify_affine_plane_specialization.py),
[receipt](receipt-affine-plane-specialization.json) and
[context](sources/affine-plane-specialization-context.json) provide exact
symbolic transcription and scope controls. They do not certify Tang's
geometric proof or the all-degree lattice argument by finite computation.
Independent replay, formal verification and external review are not claimed.
