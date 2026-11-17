# Four simple conductor roots and two-node action rigidity

Let K=(P,Q) be an actual complex polynomial plane Keller map in a
source frame where P is monic of positive degree in x. Put
A=C[P,Q,y], R=C[x,y], and let kappa generate the exact conductor
(A:R) in R as in
[the generic-coordinate conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).

It is IMPOSSIBLE that kappa(x,c) has degree FOUR and four DISTINCT
roots for every c in C. The subsequent
[global degree-four theorem](GLOBAL_DEGREE_FOUR_CONDUCTOR_EXCLUSION.md)
excludes the whole degree-four case, including repeated roots and degree drops.
The subsequent [constant-degree-four theorem](CONSTANT_DEGREE_FOUR_CONDUCTORS_AND_PERIOD_ODE.md)
allows repeated roots and excludes every conductor whose fibre degree
stays four. This elementary simple-root argument alone does not exclude
a globally squarefree conductor with a branched projection.

The proof combines an elementary radical-conductor lemma with the
genus-two row of Neumann's unconditional low-genus classification.
Its action argument uses two fixed node pairs; it does not require
a uniform target automorphism for the whole family.

## 1. The specialized conductor is exact

The monogenic bases in the linked conductor theorem remain bases on
EVERY source line. Indeed P_c is monic of the original degree and
Q_c is primitive over C(P_c), by the polynomial-immersion argument.
Thus A/(y-c)=C[P_c,Q_c] inside C[x], with normalization C[x].

The same one-variable monogenic trace-dual calculation gives

    (C[P_c,Q_c]:C[x])=(F_V(P_c,Q_c,c)/P_c')*C[x]
                    =kappa(x,c)*C[x],                     (1)

up to a nonzero constant sign. The final equality specializes the
actual Keller identity F_V=-P_x*F_T. Consequently kappa_c is the
FIBRE conductor, without an assumed general commutation of
normalization with base change.

The index-square theorem also gives

    dim_C(C[x]/C[P_c,Q_c])=deg_x(kappa_c)/2.                (2)

Under the stated hypothesis this dimension is TWO on every fibre.

## 2. A squarefree conductor of length four gives exactly two nodes

Here is the elementary ring statement needed in each fibre. Let
A_c subset C[t] be a finite birational order, with squarefree
conductor I=f*C[t] of degree four, and suppose
dim_C(C[t]/A_c)=2. Then

    C[t]/I=C^4, dim_C(A_c/I)=2.

Every unital C-subalgebra of C^4 consists of functions constant on
the blocks of a partition of its four coordinates. To check this,
declare two coordinates equivalent when every subalgebra element
has the same value there; polynomial interpolation in finitely many
separating elements gives each block idempotent.

There are exactly two blocks because the subalgebra has dimension
two. Neither block is a singleton. A singleton coordinate idempotent
would generate a nonzero ideal of C^4 contained in A_c/I, contradicting
that I is the FULL conductor. Therefore the blocks have sizes TWO
and TWO.

It follows that A_c is EXACTLY the algebra obtained by identifying
two disjoint pairs of points of A^1. Each resulting singularity is
an ordinary node: its completed local ring is two power-series
branches with their constant terms identified. There is no hidden
jet condition, third branch, tangency or further singularity.

More generally, the same partition argument shows that a squarefree
conductor of degree 2*delta in a codimension-delta order has only
two-element blocks. Only the degree-four case is used below.

## 3. The four roots and their pairing become fixed

Write kappa=a_4(y)*x^4+... . Every fibre has degree four, so a_4
has no complex zero and is a nonzero constant. Divide by it to
make kappa monic. Its discriminant has no zero either, so it too is
a nonzero constant.

The four-root algebra is finite etale over A^1_y. Equivalently,
its roots form an unramified finite analytic covering of C.
Since C is simply connected, the roots are four entire functions
h_i(y). The elementary monic-root bound gives polynomial growth
in y, so Cauchy's estimates make every h_i a POLYNOMIAL.
Each h_i-h_j is a polynomial without a zero, hence a nonzero
constant. One source shear u=x-h_1(y) therefore makes the four
roots fixed distinct complex numbers a_1,...,a_4.

Their pairing is fixed as well. For each i, the three sets

    Z_ij={c:P(a_i,c)=P(a_j,c), Q(a_i,c)=Q(a_j,c)}, j!=i,

are closed subsets of A^1 and cover A^1. Irreducibility of A^1
makes one of them the entire line. Section2 says that a conductor
point has exactly one mate on every fibre, so a second partner
cannot appear even at a special parameter. Thus the SAME two
pairs are identified for every y.

In particular the corresponding endpoint equalities are polynomial
identities in y. This argument does not infer flatness from
unramifiedness or posit a locally constant pairing.

One may also recover the entire order. Let T be the fixed two-pair
pinching algebra over C[y]. Then A subset T and both contain
kappa*R. The quotient T/A is finite over C[y]. Since R/T is free
over C[y], specialization identifies T_c with the corresponding
fixed pinching algebra inside C[u]. Section2 gives A_c=T_c.
Every fibre of T/A is therefore zero, and Nakayama gives A=T.
The action contradiction only requires the endpoint identities.

## 4. Neumann's genus-two row forces a common midpoint

The exact primary inputs are retained in
[the two-node source audit](sources/neumann-two-node-genus-two-source-audit.json).
The cached sources and their precise roles are:

- Neumann, Invent. Math.98(1989), Section8, page482: the low-genus
  table classifies defining polynomials up to POLYNOMIAL ambient
  automorphisms and affine changes of the polynomial value.
  The only realizable genus-two row is the hyperelliptic quintic.
  The other genus-two knot row is explicitly unrealizable.
- Neumann-Rudolph, Math. Ann.282(1988), corrected Lemma7.1,
  pages349 and351: a reduced fibre with one place at infinity
  has a good defining polynomial.

The latter hypothesis holds because the normalization here is A^1
and has just one point at infinity. The twice-pinched affine line
has Euler characteristic -1. Smoothing its two ordinary nodes
decreases this to -3, while regularity at infinity retains one
puncture. Thus a nearby smooth fibre has genus TWO, so the table
applies. The singular fibre itself is not assumed regular.

After a polynomial target change the two-node fibre consequently has
equation Y^2=G(X), where G has degree five. Two ordinary nodes and
rational normalization force

    G(X)=(X-alpha)^2*(X-beta)^2*(X-gamma),

up to a nonzero leading constant, with alpha,beta,gamma distinct.
Its normalization is

    X=tau^2+gamma,
    Y=tau*(tau^2-R_1)*(tau^2-R_2),                         (3)

where R_1=alpha-gamma and R_2=beta-gamma are distinct and nonzero;
a nonzero target scaling absorbs the leading constant.

The isomorphism of affine-line normalizations has
tau=a*u+b, a!=0. Both node pairs in the tau coordinate have midpoint
zero. Their preimages in u therefore have the SAME midpoint -b/a.
Since the two node pairs are fixed, translate that common midpoint
once and write their endpoints as

    {-r,r}, {-s,s}, r*s*(r^2-s^2)!=0.                      (4)

Absorb the scale a from tau=a*u into the target coordinates of(3).
Every generator pair of this fixed order is then carried by some
polynomial target automorphism to the SAME canonical pair

    p_0=u^2,
    q_0=u*(u^2-r^2)*(u^2-s^2).                            (5)

This conclusion allows a different target automorphism on each fibre.
No classification of arbitrary higher-genus nodal curves is used.

For an independent check of the exact algebra in(5), split any
polynomial as E(u^2)+u*O(u^2). Equal values at both pairs in(4)
mean O is divisible by (u^2-r^2)*(u^2-s^2). Thus the fixed order is

    C[u^2] direct-sum
    u*(u^2-r^2)*(u^2-s^2)*C[u^2]
    =C[p_0,q_0].

## 5. A fixed nonzero period contradicts the Keller variation

For any fixed endpoints identified by a generator pair, a polynomial
target automorphism of constant Jacobian lambda multiplies its
action integral by lambda. This follows because

    F*dG-lambda*X*dY

is polynomial-exact; its primitive has equal endpoint values.
The proof is independent of the chosen pair of endpoints.

For the canonical pair(5), direct integration gives

    A_r=integral_(-r)^r p_0*dq_0
       =8*r^5*(3*r^2-7*s^2)/105,

    A_s=integral_(-s)^s p_0*dq_0
       =8*s^5*(3*s^2-7*r^2)/105.                          (6)

These numbers cannot BOTH vanish. Simultaneous vanishing would give
3*r^2=7*s^2 and3*s^2=7*r^2, forcing r=s=0. Choose one nonzero
canonical period, with radius v equal to r or s.

For every source fibre, the corresponding action

    I_v(y)=integral_(-v)^v P(u,y)*Q_u(u,y)*du

is therefore nonzero: the fibrewise target change in Section4
multiplies it by a nonzero constant. Yet the fixed endpoint
identities and constant source Jacobian J give

    I_v'(y)=[P*Q_y]_-v^v
             -integral_(-v)^v(P_u*Q_y-P_y*Q_u)*du
            =-2*v*J.                                     (7)

The source shear and translation preserve J!=0, and v!=0.
Thus I_v is an affine polynomial with a complex zero, contradicting
its nonvanishing on every fibre. This proves the theorem.

An INDIVIDUAL canonical period can vanish: r^2=7,s^2=3 gives A_r=0
although the algebra has exactly two ordinary nodes and A_s!=0.
The proof selects a fixed nonzero period; it does not claim
nonvanishing at every node of a two-node curve.

## 6. Finite unramified orders can have degree drops and repeated roots

The finite-order architecture alone does not remove the remaining
cases. Here is one exact control, with the actual source coordinate y:

    R0=C[x,y], s=x^2, D0=C[s,y]=C[P0,y],
    g(s,y)=(s-1)*(y*s-1),
    P0=x^2+y, Q0=x*g(x^2,y), Adeg=C[P0,Q0,y].             (8)

Its raw order and normalization are

    Adeg=D0 direct-sum g*x*D0,
    R0=D0 direct-sum x*D0, x^2=s.

The normalization is finite, birational and monogenic. Its conductor
is EXACTLY g*R0. Indeed for a+b*x in R0, membership of both this
element and its product with x in Adeg forces b and a respectively
to be divisible by g. The converse is immediate. The inclusion
matrix is diag(1,g), so its index is g(s,y), and the conductor norm
is the exact square g(s,y)^2.

This normalization is also UNRAMIFIED everywhere: the relative
derivatives of P0,Q0 in x have no common zero. If P0_x=2x is zero,
then x=0 and Q0_x=g(0,y)=1. The source pairs are exactly x and -x.
Their P0 values agree, and their Q0 values agree precisely when
g(x^2,y)=0; x=0 is never a conductor point. Thus the same-line
closed-collision support is exact, on every fibre.

The conductor in the source is

    g(x^2,y)=(x^2-1)*(y*x^2-1).

It has degree FOUR generically, degree TWO at y=0, and two double
roots at y=1. Every fibre is nonempty and finite; its degree is even.
For y different from 0 and1 the four roots are distinct. Its exact
x discriminant is 16*y*(y-1)^4. The index degree similarly drops
from two to one at y=0. This control therefore has both kinds of
degeneration left outside the theorem above.

The plane target map (P0,Q0) is NOT Keller. Its source Jacobian is

    J0=2x^6-(5y+2)x^4+3(y+1)x^2-1.                       (9)

Its monic relation, with independent target and pencil variables U,V,T,
is

    F0(U,V,T)=V^2-(U-T)*[(U-T-1)*(T*(U-T)-1)]^2.

The exact identities in the source are

    (F0_V/P0_x)(P0,Q0,y)=g(x^2,y),
    F0_T(P0,Q0,y)=-J0*g(x^2,y).                          (10)

Thus the trace-dual ratio still gives the conductor, but replacing it
by F0_T would introduce the nonconstant Jacobian factor. This isolates
the missing actual Keller condition; it does not refute a future
argument using that condition globally.

Flatness alone cannot repair the gap. In fact the ACTUAL setup already
makes R/A flat over C[y]. The specialized inclusion A/(y-c) into
R/(y-c) is injective for every c by Section1, so R/A has no
(y-c)-torsion. Over C it is therefore torsion-free over the PID C[y],
hence flat. Likewise R/(kappa) is torsion-free over C[y], because
kappa has no factor y-c. A torsion-free module over a PID is a
filtered union of finite free submodules, which verifies this flatness
without a finite-generation assumption.

Neither quotient has been shown FINITE over C[y]. In the control (8),
R0/Adeg is (D0/g*D0)*x and is flat over C[y], yet its fibre dimension
drops from two to one. The moving conductor component y*x^2=1
escapes to x infinity as y tends to zero; its image is the nonclosed
set y!=0. Its projection, and hence the full conductor projection,
is not finite. Constant fibre length would require an additional
properness or finiteness argument, not just flatness, finite
normalization over D0, or reduced normalization fibres.

## 7. Verification and remaining cases

The exact checker pins the two primary PDFs and the new source audit.
It checks the canonical two-node equation, conductor, index, period
formulas and non-simultaneous vanishing; explicit target changes with
unrestricted added shear degrees; root shears and the discriminant;
the complete four-point partition alternatives, and the exact finite
unramified degeneration control (8)--(10), including its conductor,
index square, repeated fibre and degree drop.

The classification and analytic-cover arguments are the written
proofs and attributed primary results, not consequences of finite
samples. The separate global period-ODE theorem excludes degree-four conductors
also when fibre degree drops. Higher conductor degrees,
automatic trace integrality, HC4 and JC2 remain unresolved.
Independent replay, formal verification and external review are
not claimed.
