# No compatible quartic coordinate has pole order two

14 November 2026. PF2RULED5 excludes a target coordinate of pole order
two on the compatible F1 completion, without bounding the pole order
of the other coordinate. Pole order one was already excluded by
triangular reduction and the intersection bound B^2=3.

Thus every coordinate in every polynomial target presentation has pole
order at least three. This removes (2,5) and every later (2,n) pair;
the subsequent [complete (3,5) proof](POLE_THREE_FIVE_EXCLUSION.md)
is now complemented by [PF2RULED7](POLE_SEVEN_EXCLUSION.md), excluding
every height through seven. The first open pairs are (6,8),(7,8).
It does not exclude arbitrary low-pole elements of the target algebra
that are not members of a polynomial coordinate system.

## Setup and the square-root leading symbol

Use S=F1, B=C0+2f, E=C0+3f, with E smooth and B.E=4 at one point q.
Let N=S\B, and suppose pi=(P,Q):N -> A2 is finite of generic degree four,
with

\[
 \operatorname{div}(dP\wedge dQ)=E-3B.
 \tag{1}
\]

Write V_m=H0(S,O(mB)) and L=O_B(B)=O_P1(3).
Assume delta(P)=2. If delta(Q) is even, the
[divisible-pole rule](POLE_REDUCTION.md) subtracts a polynomial in P
to lower it. Repeating gives delta(Q)=1 or an odd n=2r+1>=3.
The first case reduces both coordinates to V1 and contradicts degree four.

In the odd case, the leading symbols are, after nonzero target scalings,

\[
 \overline P=\tau^2,\qquad\overline Q=\tau^{2r+1},
 \qquad\tau\in H0(L)\setminus\{0\}.
\]

Lift tau to T in V1 and write P=T^2+A, A in V1. Let a be the leading
symbol of A, allowing a=0. Every zero of tau has multiplicity at most
three. The proof below keeps all of them, including the marked point and
the point at infinity.

## A polar centralizer lemma excludes proportional symbols

If a=mu*tau, then A=mu*T+a0 and P is a quadratic polynomial in T.
More generally, a component f(T) of polynomial degree d>=2 in a
function T of pole one cannot have an independent companion Q whose
differential has pole at most three.

Indeed, if delta(Q)=j>=1 with symbol sigma, the leading differential
coefficient must vanish since d+j+1>3. It gives

\[
 (\,\sigma/\tau^j\,)'=0,\qquad \sigma=c\tau^j.
\]

Subtract c*T^j from Q. This preserves dP wedge dQ, because P is a
polynomial in T, and strictly lowers the pole of Q. These subtractions
are an algebraic argument, not an assertion of target equivalence.
They terminate at a constant Q, contradicting the nonzero differential.

We may therefore assume a is not proportional to tau.

## A uniform principal-part calculation

For odd n=2r+1 define the global polynomial

\[
 F_{2r+1}(T,A)=
 \sum_{j=0}^{r}\binom{r+\tfrac12}{j}T^{2r+1-2j}A^j.
\]

Put c_r=binomial(r+1/2,r), which is nonzero in characteristic zero.
The binomial recurrence gives the exact identity

\[
 dP\wedge dF_{2r+1}
 =-c_r A^r\,dT\wedge dA.
 \tag{2}
\]

Its pole is at most r+3. Starting with Q-F_(2r+1), eliminate every
remaining pole j>r as follows. The leading differential equation forces
the symbol c*tau^j. For even j subtract c*P^(j/2), which does not
change dP wedge dQ. For odd j subtract c*F_j from the remainder and
record that term in the representation of Q. Its differential has pole
at most (j-1)/2+3<=r+2. These steps strictly lower the remainder's pole.

We obtain, up to the permitted polynomial target shears,

\[
 Q=F_{2r+1}+\sum_{\substack{r<j<2r+1\\j\ {\rm odd}}}b_j F_j+D,
 \qquad D\in V_r.
\]

Only the first term and D can contribute to pole r+3 in the differential.
Writing d in H0(L^r) for the symbol of D, that coefficient vanishes:

\[
 -c_r a^r(\tau'a-\tau a')
 +2r\tau\tau'd-2\tau^2d'=0.
\]

Equivalently,

\[
 -2\tau^{r+2}\left(
 \frac d{\tau^r}
 -\frac{c_r}{2(r+1)}\frac{a^{r+1}}{\tau^{r+1}}
 \right)'=0.
\]

Hence

\[
 d=\frac{c_r}{2(r+1)}\frac{a^{r+1}}{\tau}
   +\lambda\tau^r.
 \tag{3}
\]

All symbols are regular sections, so (3) forces
tau to divide a^(r+1). In particular **a vanishes at every zero of tau**.
The finite checker tests binomial instances and the general differential
identity; the recurrence and the strictly lowering construction above
prove the statement for every odd n.

## All nonproportional zero configurations

Since tau and a are cubic sections, every zero of tau is a zero of a,
and a is not proportional to tau, their local multiplicities at zeros of
tau have only these possibilities:

| Zeros of tau | At the multiple zero | Other zero |
| --- | --- | --- |
| double plus simple | k=2, h=1 | k=1, h>=1 |
| triple | k=3, h=1 | none |
| triple | k=3, h=2 | none |

Here k=ord(tau), h=ord(a). If tau were squarefree, vanishing of a at
all three zeros would make it proportional. In the double-plus-simple
case, h>=2 at the double zero would likewise force proportionality.
These observations include a=0 in the already excluded branch.

## A generic first-coordinate fibre and its branches at infinity

Take a generic value c of P. The normalization of the completed fibre
has meromorphic Q of total degree four, summed over its components if
necessary. This follows by base change from the finite degree-four map.
We do not assume geometric irreducibility of the generic fibre.

All points at infinity of this fibre lie over zeros of tau, since
the pencil numerator restricts to tau^2 on B. At such a zero choose
local coordinates z=0 for B and u along B. Write

\[
 T=\tau(u)/z+\text{regular},\qquad
 A=a(u)/z+\text{regular}.
\]

The fibre numerator has the form

\[
 z^2(P-c)=\tau(u)^2+
 z\big(a(u)+2\tau(u)t_0(u)\big)+z^2 C(u,z,c).
 \tag{4}
\]

If 0<h<k, Newton balance gives exactly two smooth local branches:

\[
 \begin{array}{c|c|c}
 &e=\operatorname{ord}_u z&\operatorname{ord}_u P_z\\ \hline
 \mathrm I&2k-h&3h-4k\\
 \mathrm{II}&h&-h .
 \end{array}
 \tag{5}
\]

For branch I, set z=u^(2k-h)*v in (4). After division by u^(2k),
the equation at u=0 is tau0^2+a0*v=0 with nonzero derivative a0,
so the implicit-function theorem gives one branch. For branch II,
the generic nonzero coefficient C(0,0,c) gives the other solution
z=u^h*v, with nonzero leading v. Their intersection multiplicity is h;
their multiplicities sum to the full quadratic Weierstrass degree in z,
so there are no further branches.

At the simple zero k=1 with h>=1, a generic c gives two distinct
branches with e=1 and ord(P_z)=-1. These conclusions are stable under
all regular terms in (4); generic c avoids the quadratic discriminant
and zero leading coefficients.

## The volume divisor determines Q's poles

Locally write

\[
 dP\wedge dQ=U(z,u)\,\ell(z,u)z^{-3}dz\wedge du,
\]

where U is a unit and ell is an equation of E. On a fibre P=c,

\[
 \frac{dQ}{du}=\frac{U\ell z^{-3}}{P_z}.
 \tag{6}
\]

On branch I of (5), its order is ord(ell)-2k. On branch II,
it is ord(ell)-2h. At a simple zero it is ord(ell)-2.

If the boundary point is not q, ell is a unit. The total pole orders
of Q are therefore:

| Configuration | Poles at the multiple zero | Poles at other zeros | Total |
| --- | --- | --- | --- |
| double plus simple | 3,1 | 1,1 | 6 |
| triple, h=1 | 5,1 | none | 6 |
| triple, h=2 | 5,3 | none | 8 |

A pole of order p has derivative of order -p-1 in characteristic zero.
Every total contradicts degree four when q is not a zero of tau.

Suppose instead the zero is q. Contact four means
ell=z-phi(u), up to a unit, with ord(phi)=4. On a branch with e<4,
ord(ell)=e; for e>4 it is four; for e=4 it is at least four.

In the double-plus-simple configuration, putting q at the double zero
gives derivative order -1 on branch II (e=h=1). Putting q at the
simple zero gives the same order on each simple branch. A derivative
of a meromorphic function cannot have a nonzero u^-1 coefficient.
Both positions are impossible.

For a triple zero with h=1, branch II again gives derivative order -1,
so q cannot be that point.

For a triple zero with h=2, branch II has e=2 and Q has pole order one.
Branch I has e=4 and ord(ell)>=4. Equation (6) gives derivative order
at least -2, so its Q pole has order at most one. This includes possible
leading cancellation in ell: it can only lower the pole, or give the
impossible derivative order -1. Since these are all the points at infinity,
the total degree is at most two, again contradicting four.

This exhausts the position of q and every zero configuration.

## Consequence for the programme

There is no compatible quartic map with a target coordinate of pole two,
at any height of the other coordinate. Pole one is excluded by the
earlier triangular reduction. Therefore

\[
 \boxed{\delta(P)\geq3,\qquad\delta(Q)\geq3
 \text{ in every polynomial target presentation}.}
\]

Together with the complete pole-height-four exclusion, the first open
height-five pairs were (3,5) and (4,5). The subsequent
[PF2RULED6 proof](POLE_THREE_FIVE_EXCLUSION.md) closes (3,5);
[PF2RULED7](POLE_SEVEN_EXCLUSION.md) then excludes every height through
seven, leaving (6,8),(7,8) as the first open reduced pairs.
The unrestricted higher-pole
problem and actual ruled compatibility remain open.

The proof uses the actual quartic fibre degree and the full divisor
E-3B; an arbitrary abstract A2 complement does not supply either.
The [checker](verify_low_coordinate.py) verifies the uniform identities,
local Newton models, residue obstruction and degree tables. The
all-height principal-part construction, coverage, generic-fibre argument
and geometric interpretation are written proofs, without formal
verification, independent whole-proof implementation or external review.


The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
