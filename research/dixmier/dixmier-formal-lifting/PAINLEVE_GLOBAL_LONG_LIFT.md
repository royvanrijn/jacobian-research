# An explicit global Painleve-I formal lift of Long

14 November 2026. **QPL1.** Long's actual four classical images admit
an explicit all-orders polynomial-coefficient quantum lift with
fixed momentum orders
\[
                      (R,S_h,T_h,D_h):\quad(0,2,3,5).          \tag{1}
\]
All six Weyl relations hold. A single scalar Painleve-I recurrence
determines the construction. The fourth image is globally
polynomial at every order, including across R=0.

This is not a bounded Weyl endomorphism: for every n>=1 the h^(2n)
coefficient of S_h is nonzero and has total spatial degree exactly
8n-1. The recurrence even gives factorial-squared coefficient
growth. The construction therefore isolates position-degree growth
while keeping the momentum orders fixed. It does not prove that
every full formal representative has unbounded degrees.

The string-equation framework is classical; see
[Takasaki, Section 2](https://arxiv.org/html/nlin/0610073v2).
The result here is the explicit global four-image realization of
the retained Long map, its finite admission certificate and its
all-orders degree formula. No novelty priority is asserted.

## Global operators and conventions

Use the full Weyl algebra on positions X,Q and momenta P,Z,
with [P,X]=[Z,Q]=h and the other generator commutators zero.
Normal PBW coefficients are on the left. Define
\[
 \begin{split}
 R&=X(2-3XQ),\quad A_0=2-3XQ=R/X,\\
 B&=3X^2P+2(1-3XQ)Z,\\
 D_0&=(1+3XQ)P/2-3Q^2Z .
 \end{split}                                                   \tag{2}
\]
R/X is the displayed polynomial A_0; no X inverse is needed
to define it. The exact commutators are
\[
 [B,X]=3hX^2,\quad [B,R]=0,\quad [D_0,R]=h,\quad
 [D_0,X]=h(3-R/X)/2,\quad
 [D_0,B]=h(3/X-3R/(2X^2))B.                                  \tag{3}
\]
The rational expressions in the last line simplify to their
original polynomial expressions using R=X(2-3XQ).

In the polynomial R centralizer put
\[
 w=A_0^3-5A_0^2+9A_0-2-RXB-\tfrac32 RhX^2,\qquad
                   x_e=8A_0,\quad \epsilon=192R^3h.           \tag{4}
\]
These are global polynomials. Work temporarily where R and X
are invertible, so x_e is invertible as well. Equation (3) gives
[w,x_e]=epsilon/x_e. Thus
\[
                t=x_e^2/2,\qquad [w,t]=\epsilon.              \tag{5}
\]
Equivalently y_e=(x_e w+w x_e)/2 has [y_e,x_e]=epsilon and
w=(x_e^{-1}y_e+y_e x_e^{-1})/2. The classical pair (x_e,y_e)
is precisely QEI1's normalized input to the three-isogeny.

## The single scalar recurrence

Let
\[
 W(t;\epsilon)=\sum_{n\ge0}a_n\epsilon^{2n}x_e^{1-5n},
 \qquad a_0=-3,
\]
where
\[
 a_n={1\over6}\sum_{j=1}^{n-1}a_j a_{n-j}
       +{2\over9}\bigl(25(n-1)^2-1\bigr)a_{n-1}\quad(n\ge1).
                                                               \tag{6}
\]
The first coefficients are a_1=2/3 and a_2=98/27.
With derivative d/dt=x_e^{-1}d/dx_e, this is exactly the
coefficient equation for
\[
                 \tfrac32 W^2+2\epsilon^2 W''=27t.           \tag{7}
\]
At each positive order the coefficient of a_n is -9,
so the normalized homogeneous formal branch is unique.

The standard quadratic and cubic string operators are
\[
        U=4w^2+W,\qquad
        V=8w^3+\tfrac32(Ww+wW).                               \tag{8}
\]
Direct commutation, or the third Moyal derivative, gives
[V,U]=epsilon*(3WW'+2epsilon^2 W''')=27epsilon by (7).
No operator-order truncation is involved: W is a multiplication
series and (8) has orders two and three in w at every order.

Every a_n for n>=1 is positive. For n>=2,
\[
 a_n\ge {16\over3}(n-1)^2a_{n-1}
 \ge {2\over3}\left({16\over3}\right)^{n-1}((n-1)!)^2.          \tag{9}
\]
The first inequality uses the positive convolution and
25(n-1)^2-1>=24(n-1)^2. Hence the coefficients never terminate.
Over complex constants the one-variable series in
epsilon^2/x_e^5 has zero convergence radius. This is a formal
construction; convergence is not needed for its algebraic
relations.

## Every potential term is a global polynomial

The scaling in (4) clears every apparent denominator:
\[
 F_n:=a_n\epsilon^{2n}x_e^{1-5n}
       =8a_n(9/8)^n h^{2n}R^{n+1}X^{5n-1}\quad(n\ge1).
                                                               \tag{10}
\]
In particular F_1=6h^2R^2X^4 and
F_2=(147/4)h^4R^3X^9. Thus
W=-3x_e+sum_(n>=1)F_n is a global polynomial-coefficient
series. Its positive h orders never require new position
divisors.

Set
\[
       S_h={16-U\over12R},\qquad
       T_h={V+6U-32\over432R^2}.                              \tag{11}
\]
The finite base terms W=-3x_e give global polynomial
operators S_0,T_0, including their ordering corrections.
Their finite exact expansions are retained in the admission
certificate below.

For every positive n the added contributions are
\[
 (S_h)_n=-{F_n\over12R},\qquad
 (T_h)_n={F_n(w+2)+(w+2)F_n\over288R^2}.                      \tag{12}
\]
Both are polynomial: F_n contains R^(n+1), R is central
in this algebra, and n>=1. The S contribution is positional
and divisible by R^n; the T contribution has momentum order
one and is divisible by R^(n-1).

Equation (8) implies [S_h,T_h]=h. Both commute with R.
Their classical reductions are the actual retained S,T,
as follows from the normalized isogeny or the exact full
coefficient comparison in the checker.

## A fixed rational connection supplies the fourth image

Define K=D_0+C_K, where, in the local PBW order of (3),
\[
 \begin{split}
 C_K={}&\left(-{1\over2X^2}+{R\over6X^3}
                              -{1\over15RX}\right)B\\
       &-{3\over4X^4}+{R\over5X^5}+{3\over5RX^3}
                         +{1\over5R^2X^2}+{h\over10R}.
 \end{split}                                                   \tag{13}
\]
A direct use of (3) gives
\[
 [K,R]=h,\quad [K,x_e]={6h\over5R}x_e,\quad
 [K,w]={3h\over5R}w,\quad
                         [K,\epsilon]={3h\over R}\epsilon.   \tag{14}
\]
The weights (x_e,w,epsilon)=(2,1,5) make every term of U
homogeneous of weight two and every term of V of weight
three. Consequently
[K,U]=6hU/(5R), [K,V]=9hV/(5R), and
\[
 [K,S_h]/h=S_h/(5R)-8/(5R^2),\qquad
 [K,T_h]/h=-T_h/(5R)+S_h/(10R^2).
\]
Now put
\[
 \widehat A_h={2R(S_hT_h+T_hS_h)-32T_h-S_h^2\over20R^2},
                         \qquad D_h=K+\widehat A_h.           \tag{15}
\]
Using [S_h,T_h]=h, direct commutation of A_h with S_h,T_h
cancels the two displayed K commutators. Thus
[D_h,R]=h and [D_h,S_h]=[D_h,T_h]=0. Along with (11),
these are all six full rank-two Weyl relations, in every
formal order. The source order is (R,T_h,D_h,S_h).

This uses the exact elliptic normalizer of
[QEI1](ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md), but the
commutator cancellation above is also a direct proof.

## Global admission reduces to a finite exceptional seed

K and A_h were defined rationally. Here the polynomial
admission argument must be supplied; generic-fiber existence
alone would not suffice.

Write S=S_base+sum S_n, T=T_base+sum T_n as in (12).
The base is fixed, polynomial and finite. S_n is divisible
by R^n, and T_n by R^(n-1). In (15), every quadratic term
in positive indices i,j is therefore divisible, after the
R^2 denominator, by at least R^(i+j-2); this is nonnegative.

The terms linear in index n have, in order, factors at
least R^(n-1), R^(n-2), R^(n-3), or R^(n-2), depending
on whether they come from S_n*T_base, S_base*T_n,
the standalone T_n, or S_base*S_n. Hence all linear
terms with n>=3 are global polynomial operators.

Only the base and n=1,2 can have negative R powers or
uncancelled X poles. These are finite operator identities.
The [retained admission seed](painleve_long_seed.json) expands
the exact operators with W=-3x_e+F_1+F_2, including every
ordering term. In PBW order X^i Q^j B^k h^l the numbers
of nonzero terms are 21,48,160 for S,T,D-D_0 respectively,
and every i,j,k,l is nonnegative.

The [Fraction checker](check_painleve_long.py) verifies this
by complete exact coefficient arithmetic, not by sampled
positions or a support search. Its local multiplication is
the all-integer formula
\[
 \begin{split}
 &(R^aX^bB^ch^d)(R^{a'}X^{b'}B^{c'}h^{d'})\\
 &\quad=\sum_{j=0}^c {c\choose j}3^j(b')^{\overline j}
       R^{a+a'}X^{b+b'+j}B^{c+c'-j}h^{d+d'+j},
 \end{split}
\]
where (b')^(overline j)=b'(b'+1)...(b'+j-1), with empty
product one. This is valid for all integer position exponents.

The retained seed also compares the entire h=0 four-image
tuple with the actual Long source, including its corrected
D, not the uncompletable D_0. Thus the admitted D_h has
the required classical reduction.

Combining this finite exceptional certificate with the
all-n factors proves global polynomial admission in every
order. There is no generic-period-only inference and no
unchecked higher-order admission condition.

## Fixed momentum orders, sharp spatial growth

S has B order exactly two, T exactly three and D exactly
five; their classical leading orders are nonzero. The
positive-index S_n have order zero and T_n order one.
The quadratic expression (15) cannot raise D above five.
Since B is first-order in the original P,Z, these are also
the original ambient momentum orders in (1).

The base S has no normal-PBW h powers above two. At h^2
the exact coefficient, including its ordering term, is
\[
                         [h^2]S_h=-{17\over4}RX^4.           \tag{16}
\]
For every n>=2,
\[
 [h^{2n}]S_h=-{2\over3}a_n(9/8)^n R^nX^{5n-1},
 \qquad \deg([h^{2n}]S_h)=8n-1.                              \tag{17}
\]
Indeed R has highest total degree three and its leading
monomial cannot cancel in that nonzero product. Equation
(16) has the same degree formula at n=1. The coefficient
growth in (9) is therefore visible in a single explicit
global image, not just in an auxiliary gauge.

The cubic energy V^2-U^3 also has scalar differential order exactly
two throughout this formal construction. For arbitrary W, its top
w coefficient is -3W^2-4epsilon^2 W'', which equals -54t by (7);
all w powers above two cancel identically. Its coefficient series
still have unbounded spatial support. Thus finite momentum order,
a formal quadratic energy and full global fourth-image admission
do not supply the algebraic generic coefficient field required by
[QSO2](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md) and the bounded
energy restriction of QEI1. The construction lies in the completed
coefficient field, outside that algebraic category.

This coefficient-field distinction can be proved, rather than
inferred just from unbounded support. Let A(z)=sum a_n z^n.
The lower bound (9) makes A nonalgebraic over C(z), since an
algebraic power series has positive convergence radius.
Nonalgebraicity remains after any characteristic-zero constant
field extension: a polynomial relation of fixed degrees would
give a homogeneous linear system with rational entries for
its finitely many coefficients, and a nonzero solution after
extension would already give one over the original field.

The positional tail of S is, up to a nonzero rational factor
and finitely many initial terms, A((9/8)h^2 R X^5). If it were
algebraic over K(R)((h))(X), clear the polynomial denominators.
The resulting identity can also be read X-adically: after
one common h-power is removed, each finite expression is a
bivariate formal series with the same coefficients in either
order. The change z=(9/8)h^2 R X^5 is finite algebraic over
the scalar field K(R)((h)). It would make A(z) algebraic over
that field's rational z-line, contradicting the constant
extension argument. Thus this is genuinely outside the
generic algebraic coefficient category, even though its
individual h coefficients are polynomials.

The tuple has no common spatial degree bound. It cannot
be made bounded by regular finite source and target Weyl
automorphisms whose inverses also have finite spatial
support: otherwise composing the bounded tuple with those
finite inverses would make this original tuple bounded.

Arbitrary unbounded formal Hamiltonian gauges are a larger
class. [QFL1](PROOF.md) compares full formal lifts using
that class, and [QNI1](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md)
explains why boundedness cannot be transferred through it.
A full formal representative reached through an unbounded
gauge might still be bounded; (17) does not decide that.

## Energy motion and what remains

In the normalized (x_e,y_e) rational Weyl chart, the first
corrections of this construction are
U_2=-1/(3x_e^4), V_2=-4y_e/x_e^5. Its cubic energy starts
\[
 V^2-U^3=-27(y_e^2-x_e^3)-9\epsilon^2/x_e^2
                                               +O(\epsilon^4).
\]
Thus it escapes
[QEM1](ELLIPTIC_QUANTUM_ENERGY_MOTION.md) by actual
transverse energy motion; the energy commutator starts
at the required third order. This is not a fixed-energy
quantization and not an algebraic-coefficient quadratic
mate of the kind excluded by QSO2.

The new constructive information is an explicit global,
fully admitted, bounded-momentum all-orders carrier.
Its precise failure is position-degree growth. The next
question is whether a structurally different finite
higher-order carrier can avoid that growth, or whether
a comparison theorem under additional verified hypotheses
can force every bounded candidate into this Painleve
class. No such comparison is supplied here, and none
follows just from formal gauge equivalence.

## Verification boundary

Exact controls verify the forced two exceptional Painleve
coefficients, the complete finite polynomial admission seed,
the exact Euler identities, every full relation through
h^6, the generic string identity and the complete classical
Long tuple. All-orders relations, admission beyond the
seed, positivity, divergence and degree growth use the
written proofs above. Checking h^6 alone is not claimed
to establish them.

No common spatial bound, finite Weyl endomorphism from Long,
all-representative nonlifting, DC2 solution, independent
theorem replay, formal verification, external review or
novelty priority is claimed.
