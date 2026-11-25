# Quintic sphere phases with an unstable leading binary form

Let \(K\) have characteristic zero and work in
\[
A_K=K[x,y,t]/(t^2+2xy-1),\qquad
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
f\bigl(x,(1-t^2)/(2x),t\bigr)\,dt.
\]
The integral is defined by its rational values on monomials. A pure phase
has \(\Lambda(f^m)=0\) for every \(m>0\). A safe phase has eventual mixed
vanishing \(\Lambda(gf^m)=0\) for each fixed \(g\in A_K\).

**Theorem.** Let \(f\) have least sphere degree five. If its leading binary
decic has a geometric root of multiplicity at least six, purity implies
safety. Here the leading binary form is the degree-five symbol evaluated at
\((x,y,t)=(u^2,-v^2/2,uv)\) on the homogeneous cone.

The proof uses the [rank-two cyclotomic ordinary-prime criterion](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md) for the
four curve orders \(7,9,11,15\) isolated below. The remaining cases use
ordinary reduction in smaller dimension, exact binomial coefficients, and
two polynomial shears. The leading-form hypothesis allows arbitrary lower
terms of \(f\).

## The supporting lines

The linear trace-preserving \(\mathrm{SL}_2\) action on the cone and whole
quadric is given in [Spherical leading binary forms](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#2-the-full-sphere-filtration-and-trace-preserving-changes).
Move the root to \(v=0\) after finite scalar extension. A degree-five normal
monomial of angular exponent \(i\) has \(v\)-exponent \(5-i\) in the binary
form. Thus every degree-five term has negative angular exponent. The normal
monomials are \(y^b t^r\) for \(b\ge1\) and \(x^a t^r\) for \(a\ge0\).
On substituting for \(y\), a monomial \(y^b t^r\) can contribute to the upper
boundary only through its highest term \((-1/2)^b x^{-b}t^{2b+r}\).
The allowed upper-boundary nodes \((i,j)\) of \(x^it^j\) therefore satisfy
\[
 -5\le i\le4,\qquad \max(0,-2i)\le j\le (5\text{ if }i<0\text{ else }4)-i.
\tag{1}
\]

If no positive term occurs, every power moment is the interval moment of
the angular core. The [interval polynomial theorem](https://arxiv.org/abs/1910.12148v1),
Theorem 2.1, forces that core to be zero under purity. The remaining support
has strictly negative angular exponents, which gives safety against every
fixed marker.

Otherwise the support has both angular signs. Let \((0,H)\) be the highest
point of the Newton polygon at angular exponent zero. The supporting line
can be written \(j+\rho i=H\). Since the degree-five nodes all have negative
angle, \(H<5\). Since at least one degree-five node is occupied, \(\rho>1\).
If the supporting face has angle zero only at an occupied core node and all
other occupied nodes have a single strict sign, its top angular-constant
coefficient in the \(m\)th power is the core coefficient to power \(m\).
This gives the usual prime denominator obstruction. If a face crossing both
signs has exactly two occupied nodes, with angles \(-b,a\) and heights
\(j_-,j_+\), use \(m=(a+b)(p-1)/(a j_-+b j_+)\), taking \(p=1\) modulo
\(a j_-+b j_+\). The top coefficient is a binomial times a nonzero monomial
in the two phase coefficients. The indices are less than \(p\) because
\(a j_-+b j_+\ge 2ab\ge a+b\). In both arguments the only angular-constant
height whose interval denominator is divisible by \(p\) is \(p-1\), so
the top coefficient gives a nonzero moment at all sufficiently good such
primes. These arguments apply after a number-field specialization that
preserves the selected support and its nonzero coefficients.

The following table gives every line with both signs, at least three allowed
nodes, and at least one degree-five node on or below it. The actual face must
occupy at least three of the displayed nodes and both signs. The proof of
exhaustiveness is to form the line for every pair of opposite-sign nodes in
(1), apply the two strict inequalities and the degree-five-node test, and
count all equality nodes in (1). The [exact checker](verify_unstable_quintic_faces.py)
performs these integer and rational tests.

| \(\rho\) | \(H\) | Allowed angular exponents | Test |
| ---: | ---: | :--- | :--- |
| \(5/3\) | \(5/3\) | \(-5,-2,1\) | Order five, or a shear when the quadratic repeats |
| \(2\) | \(2\) | \(-3,-2,-1,0,1\) | Quadratic height |
| \(3/2\) | \(5/2\) | \(-5,-3,-1,1\) | Order five, or a shear when the cubic is a cube |
| \(3/2\) | \(3\) | \(-4,-2,0,2\) | Order three after angular compression |
| \(2\) | \(3\) | \(-2,-1,0,1\) | Order three |
| \(3\) | \(3\) | \(-1,0,1\) | Order three |
| \(4/3\) | \(10/3\) | \(-5,-2,1\) | Order ten |
| \(3/2\) | \(7/2\) | \(-3,-1,1\) | Order seven |
| \(4/3\) | \(11/3\) | \(-4,-1,2\) | Order eleven |
| \(5/4\) | \(15/4\) | \(-5,-1,3\) | Order fifteen |
| \(4/3\) | \(4\) | \(-3,0,3\) | Order four after angular compression |
| \(3/2\) | \(4\) | \(-2,0,2\) | Order four after angular compression |
| \(2\) | \(4\) | \(-1,0,1,2\) | Order four |
| \(3/2\) | \(9/2\) | \(-1,1,3\) | Order nine |

## The cyclic-cover test

Let the angular Laurent polynomial on a face have the form
\(z^{-b}P(z^q)\), with occupied extreme angles \(-b,a\). For a noninteger
height write \(H=n/q\) in reduced form, where the table gives the same
\(q\) as the denominator of \(\rho\). The raw homogeneous degree-\(2n\)
completion of \(z^{n-b}P(z)^q\) has endpoint exponents \(n-b,n-a\) and at
each root of \(P\) has multiplicity \(q\) times its multiplicity in \(P\).
For \(p=nk+1\) use \(m=qk\). The face coefficient is
\[
 [z^{bk}]P(z)^{qk}=[z^{p-1}]\bigl(z^{n-b}P(z)^q\bigr)^k.
\tag{2}
\]
At integer height \(H=n\) first divide all angular exponents by their
common positive step when the table indicates compression. This leaves
constant terms unchanged. Now use the same formula with \(q=1\) and the
compressed \(a,b\).

If every nonzero raw branch exponent is less than \(n\), the normalized
cyclic cover \(w^n=z^{n-b}P(z)^q\) has character-one sheaf \(O(-2)\) and
Frobenius coefficient (2) at absolute degree-one good primes \(p=1\) modulo
\(n\), as in [the normalized-cover computation](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#4-the-exact-normalized-cover-frobenius-line).
Every cover used below is connected, as checked by the greatest common
divisor of \(n\) and the nonzero exponents. Ordinary reduction makes (2)
nonzero. Multiplication by \(p\) and reduction of the original moment gives
exactly (2), since the top angular-constant height is \(p-1\). Work over
a number field containing the cover data, its roots and the required roots
of unity. For coefficients in an arbitrary characteristic-zero field,
a single number-field specialization preserves all the polynomial moment
identities and the finitely many nonvanishing conditions for the chosen
face and branch partition.

### Integer heights and existing smaller covers

* At \(\rho=H=2\) all normal terms are of weight at most two for weights
  \((x,y,t)=(2,0,1)\). After exchanging \(x,y\), the literal Laurent height
  is at most two, because substitution for the new \(y\) increases height
  by twice its exponent. The [quadratic-height classification](QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md)
  gives safety for every pure phase of least degree five in this class.
* At height three the tabulated angular supports, after compression when
  indicated, are subsets of \(\{-2,-1,0,1\}\). For two extremes \(-1,1\),
  use the tied order-three quadratic-cover test in the
  [stripped Newton reduction](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#9-a-threefold-input-excludes-the-lower-degree-even-perturbation-sector).
  If the extremes are \(-2,1\), the exponents for a distinct-root cubic are
  \((1,1,1,1,2)\). The genus is three and its Jacobian has
  \(\mathbb Q(\zeta_3)\) action. Fité's Theorem 1 supplies ordinary primes
  after adjoining \(\zeta_3\), including absolute degree-one primes.
  The double-and-simple partition has exponents \((1,2,1,2)\) and genus two.
  For a triple root the coefficient is a unit without using the covering
  formula, since it is an endpoint unit times \(\binom{3k}{2k}\) with
  \(3k=p-1\).
* At height four the angular possibilities normalize to \(-1,0,1\) or a
  subset of \(\{-2,-1,0,1\}\) by compression and reversal. The full-cover
  order-four tests for these two possibilities, including repeated roots,
  are proved in the [stripped Newton reduction](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#9-a-threefold-input-excludes-the-lower-degree-even-perturbation-sector).
* The order-ten case has exponents \((5,3,3,9)\) or \((5,6,9)\), up to
  permutation. These are precisely the order-ten covers proved ordinary at
  suitable primes in the [stripped Newton reduction](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#6-excluding-the-sigma2-and-sigma43-faces).

### Order five and the shears

For height \(5/3\), distinct quadratic roots give exponents \((3,3,4)\),
genus two. For height \(5/2\), the negative extreme is \(-5\) or \(-3\),
since at least three nodes must be occupied. The positive extreme is \(1\).
Distinct roots give exponents \((4,2,2,2)\) in either case, and one double
root merges two exponents. These give ordinary primes by the four-root
and three-root cases in [Spherical leading binary forms](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md).
For clarity, this uses their cyclic-cover ordinary-reduction proofs,
not a claim that the original phase is governed by its global leading decic
in the present weighted frame.

The two omitted cases have a raw root exceeding order five. Define
\[
 \theta_{c,r}(y)=y,\quad \theta_{c,r}(t)=t+cy^{r+1},\quad
 \theta_{c,r}(x)=x-cy^r t-\tfrac12c^2y^{2r+1}.
\tag{3}
\]
These maps are trace-preserving with inverses \(\theta_{-c,r}\). Indeed their
derivation in the literal Laurent \(y\)-chart is \(cy^{r+1}\partial_t\).
Its trace is zero by the integration boundary formula, since at \(t=\pm1\)
the original \(x\) vanishes and every regular polynomial evaluates to a
polynomial in \(y\). The trace is invariant under its finite exponential.
This also gives invariance of safety, by transforming the marker.

At height \(5/3\), the face in the original variables is
\(Bx+Ay^2t+Cy^5\). A repeated root of its face quadratic is equivalent to
\(A^2+2BC=0\), so the face equals the top part of \(B\theta_{c,2}(x)\) for
\(c=-A/B\). Apply its inverse. For the integer weights \((x,y,t)=(5,1,3)\),
every lower term stays of weight less than five under the shear and normal
reduction. The resulting phase is
\[
 Bx+P_4(y)+(a+b y)t,\qquad \deg P_4\le4.
\]
The upper line from the positive node of \(Bx\) detects \(y^4\), if present,
as a unique opposite extreme. The next possible tied largest slope is
\(3/2\), shared by \(y^3\) and \(yt\). One such term alone is again a
two-term obstruction. Both terms together give the order-three quadratic
with exponents \((2,2,2)\) when the roots are distinct. If that quadratic
repeats, the coefficients of \(Bx+b y t+d y^3\) satisfy \(b^2+2Bd=0\),
and the face is \(B\theta_{-b/B,1}(x)\). The inverse shear leaves the lower
terms strictly below its face for weights \((3,1,2)\). In every exceptional
pure case the result has least degree at most two.

At height \(5/2\) the remaining raw-root collision is a triple root of the
face cubic with extreme angles \(-5,1\). On normal reduction
\[
 \theta_{c,1}(xt)=xt-\tfrac32cyt^2-\tfrac32c^2y^3t
                       -\tfrac12c^3y^5+\tfrac12cy.
\]
Its face polynomial \(P(u)\) for weights \((3,1,2)\) is
\((u+c/4)^3\). Any triple-root face is therefore the leading part of
\(B\theta_{c,1}(xt)\) for a suitable \(c\), where \(B\) is its coefficient
of \(xt\). Apply the inverse to obtain
\[
 x(Bt+a)+P_4(y)+tP_2(y)+d t^2,\qquad \deg P_i\le i.
\]
The slopes from the node of \(Bxt\) to \(y^4,y^2t,y^3\) are respectively
\(7/5,4/3,5/4\), each strictly greater than that of every other node after
the earlier ones have been removed. The two-term test removes each in
that order. The result has least degree at most two.

Every pure phase of degree at most two is safe by the
[complete quadratic sphere criterion](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md#the-full-quadratic-sphere-family-a-sharp-six-moment-criterion).
Its six-moment identities and its algebraic normal forms descend through a
finitely generated coefficient field embedded in \(\mathbb C\), since every
fixed mixed equality is rational in the coefficients. Thus both shear cases
give safety over the original coefficient field.

### Orders seven, nine, eleven and fifteen

The remaining face polynomials are quadratic with both endpoints occupied.
Their raw exponents and primitive signatures are as follows. Signature
coordinates are in the increasing order of indices prime to the order.

| Order | Distinct roots | Repeated root | Primitive signature for distinct roots |
| ---: | :--- | :--- | :--- |
| \(7\) | \((4,2,2,6)\) | \((4,4,6)\) | \((1,1,2,0,1,1)\) |
| \(9\) | \((8,2,2,6)\) | \((8,4,6)\) | \((1,1,2,0,1,1)\) |
| \(11\) | \((7,3,3,9)\) | \((7,6,9)\) | \((1,1,2,0,0,2,2,0,1,1)\) |
| \(15\) | \((10,4,4,12)\) | \((10,8,12)\) | \((1,1,0,2,0,2,1,1)\) |

The normalized eigensheaf formula gives every entry as
\(-1+\sum_j\{i a_j/n\}\). The checker verifies the signatures and the two
balanced pairs in each row. All primitive Jacobian factors for repeated-root
covers are potential CM. The nonprimitive factors are of dimensions at
most two in the distinct-root covers, with none of larger dimension in
the repeated-root cases, by the same character calculation. They have
ordinary reduction after finite extension on a density-one set.

In each distinct-root case the primitive factor has rank two over
\(\mathbb Q(\zeta_n)\), since the two multiplicities in every conjugate pair
sum to two. There are exactly two balanced conjugate pairs. The rank-two
criterion applies directly, including extra endomorphisms and nonsimple
primitive factors. Intersect with the density-one ordinary set for lower
order factors and the absolute degree-one primes after adjoining
\(\zeta_n\). This supplies the primes for (2) in all four orders.

The table exhausts mixed quintics in the unstable leading branch. Its
surviving pure phases have either one strict angular sign, quadratic
height, or a safe quadratic normal form after a shear.

The arithmetic inputs in the smaller-cover cases are Sawin,
[*Ordinary Primes for Abelian Surfaces*](https://arxiv.org/abs/1506.04784v2)
and Fité, [*Ordinary primes for some varieties with extra endomorphisms*](https://arxiv.org/abs/2005.10185v3),
Theorem 1 and Remark 12. Fité's theorem applies to the primitive
threefold in the genus-three order-three case through its imaginary
quadratic action, independently of full geometric endomorphisms.

## Exhaustion of the coefficient space

The [period obstruction](STABLE_LEADING_SPHERE_PHASE_PERIOD_OBSTRUCTION.md)
excludes leading decics with all root multiplicities at most four, and the
[exact-half-root argument](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#4-proof-of-the-elementary-obstructions)
excludes any leading root of multiplicity exactly five. Every remaining
quintic has a root of multiplicity at least six and falls under the theorem
above. Thus purity implies full-marker safety on the entire degree-five
sphere coefficient space, by reductions that allow arbitrary lower terms.
