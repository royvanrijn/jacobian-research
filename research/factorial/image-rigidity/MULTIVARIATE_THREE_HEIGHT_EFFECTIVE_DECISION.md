# Effective finite reduction for three and four independent unit-gap heights

Fix a comparable ordinary core: every unequal exponent row in
it is increasing from the first monomial to the second. The core
may be empty. Append three unit-gap row types
\[
       (A_1,A_1+1),\quad (A_2,A_2+1),\quad (A_3,A_3+1)
\]
with positive multiplicities \(r_1,r_2,r_3\) and independently
unbounded nonnegative integer heights \(A_i\).
Put \(D=\sum_i r_i\) and \(N=n+1\).

**Theorem.** For a fixed pair of adjacent orders \(n,N\ge2\),
suppose the exact rational polynomials
\[
                Q_n^C(n^D z),\qquad Q_N^C(N^D z)           \tag{1}
\]
of the core are coprime. Then there is an effectively computable
finite list of integer height triples containing every collision
for the appended profile. The hypothesis is automatic for the
empty core. Thus all-degree coprimality is effectively decidable
for the entire cooriented three-height unit-gap family with
empty core and arbitrary positive multiplicities. The decision
uses the [absolute eventual cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and the existing [two-height specialization](MULTIVARIATE_TWO_HEIGHT_FACTOR_SPECIALIZATION.md)
on finitely many bounded slices.

The same conclusion holds for four independently variable
increasing unit-gap rows with multiplicity one and empty core.
A first-order calculation handles the twisted comparator at
the one retained row. The four heights may agree.

The argument bounds at least one height by inspecting the
rescaled two-polynomial limit. This supplies an explicit slice
bound without computing any Puiseux expansions.

## A direct bound for one of the heights

As usual, strip the factorial at the smaller endpoint of every
appended row, so the coefficient at \(z^k\) becomes
\[
 c_{j,k}\prod_i(j A_i+1)_k^{r_i},
                    \qquad c_{j,k}=[z^k]Q_j^C(z).
\]
If every \(A_i\ge H>DN\), rescale by
\(z=(\prod_i A_i^{-r_i})x\). The normalized coefficient is
\(c_{j,k}j^{Dk}\) times the product of
\((1+h/(j A_i))^{r_i}\) over \(i\) and \(1\le h\le k\).
Its relative displacement from the coefficient of (1)
is at most
\[
       (1+1/H)^{DN}-1\ \le\ \frac{DN/H}{1-DN/H}.            \tag{2}
\]
Indeed every \(h\le k\le j\), and
\(\binom{DN}{b}\le(DN)^b\).

For the empty core, (1) is \((1+n^D x)^n,(1+N^D x)^N\).
Set
\[
 \delta=n^{-D}-N^{-D},\qquad \rho=\delta/3,\qquad
 \varepsilon=\min_{j=n,N}\left(\frac{j^D\rho}{2+j^D\rho}\right)^j .
\]
Take any integer \(H>2DN/\varepsilon\).
On \(\lvert x+j^{-D}\rvert=\rho\) the empty-core comparator has
modulus \((j^D\rho)^j\), whereas the perturbation is strictly
below \(\varepsilon(2+j^D\rho)^j\) by (2) and the binomial
sum. Rouché's theorem puts every root in its comparator
disk, and the two disks are disjoint. Thus at a collision
at least one height lies in \(0,\ldots,H-1\).

For a general fixed core satisfying (1), there is a rational
determinant bound instead. Let \(M\ge1\) bound the absolute
values of the coefficients of (1), set \(L=n+N\), and
let \(R\ne0\) be their resultant. Put
\[
 \varepsilon=\min\left(\tfrac12,
             \frac{|R|}{2L!\,L\,2^{L-1}M^L}\right).
\]
Use the same \(H>2DN/\varepsilon\). The original and rescaled
Sylvester matrices use fixed nominal degrees. Every entry
is at most \(2M\) in modulus, and each entry changes by
less than \(M\varepsilon\). Telescoping the difference
of each \(L\)-fold determinant product shows that the
resultant changes by less than
\(L!\,L\,2^{L-1}M^L\varepsilon\le |R|/2\).
It stays nonzero. This also bounds at least one collision
height in \(0,\ldots,H-1\).

## The finite two-height slices

For each of the three height indices and each integer in the
bounded interval just obtained, adjoin that row type with its
multiplicity to the fixed core. The enlarged core remains
comparable. The two-height specialization theorem applies to
the other two independent heights and both of their positive
multiplicities. It applies to the saturated degree-at-least-two
common-factor locus with both moment degrees full. Every
ordinary collision belongs to it. In fact one of the two
moments is an even moment and is strictly positive on the
real line. A common ordinary root is nonreal and has its
distinct conjugate as another common root. This is the
point at which the stronger rank condition is required.

Each slice gives an effective finite list of height pairs
by rational polynomial elimination. Their union is a finite
list for the triple, and actual common roots can be checked
using the factorial polynomials at the listed triples.
This proves the fixed-terminal assertion.

For all orders, first use the absolute eventual cutoff
\(J\). Terminal two is disjoint by even real positivity.
There are only finitely many adjacent pairs below \(J\).
For unrestricted counts, the
[signed degree-imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
makes \(D+\Delta_C\) bounded at a possible ordinary collision,
where \(\Delta_C\) is the fixed signed gap of the core.
For the empty core this leaves finitely many triples of
positive multiplicities and the hypothesis (1) is automatic
at every terminal. For a fixed comparable core, the same
decision applies if (1) passes its exact finite checks at
all relevant pairs and counts.

More generally, for any fixed number of comparable unit-gap
height types, apply (2) and the determinant check recursively.
Every nonzero comparator fixes one height in a bounded slice.
When two heights remain, apply the two-height theorem. The
finite recursion either gives a complete decision or leaves a
finite list of bounded cores and affected terminals whose
explicit twisted comparators in (1) have a common root.
This reduction gives no implication from such a comparator
root to a collision in an actual integer lift.

## A first-order test at simple comparator contacts

The determinant test can also cover a vanishing comparator
resultant. Suppose its common roots are all simple in both
comparators. Write \(B_j(z)=Q_j^C(j^D z)\), and let
\(\vartheta=z\,d/dz\). For each common root \(\zeta\), define
\[
       V_j(\zeta)=-
       \frac{((\vartheta^2+\vartheta) B_j)(\zeta)}
                  {2j B_j'(\zeta)},\qquad j=n,N .
\]
If \(V_n(\zeta)\ne V_N(\zeta)\) at every common root, at least
one height is again bounded by an effective constant.

Indeed, set \(\epsilon_i=1/A_i\) and \(H_1=\sum_i r_i\epsilon_i\).
The rescaled parent is a polynomial in the \(\epsilon_i\) and \(z\).
Its constant and degree-one parts in \(\boldsymbol\epsilon\) are
\[
        F_j(\boldsymbol\epsilon,z)
        =B_j(z)+\frac{H_1}{2j}(\vartheta^2+\vartheta)B_j(z)
                     + \text{terms of degree at least two}.
\]
This is the linear part of the exact product of
\((1+h\epsilon_i/j)^{r_i}\), since \(\sum_{h=1}^k h=k(k+1)/2\).
At a common simple root the analytic root in a direction
\(\boldsymbol\epsilon=t\mathbf d\) thus has first derivative
\((\sum_i r_i d_i)V_j(\zeta)\).

Let \(R(\boldsymbol\epsilon)\) be the resultant of the
rescaled parents, and let \(g\) count their distinct common
roots at \(\boldsymbol\epsilon=0\). The first nonzero
homogeneous part of \(R\) is exactly \(c H_1^g\), with
\(c\in\mathbb Q\setminus\{0\}\). To see this, evaluate along
a generic line \(\boldsymbol\epsilon=t\mathbf d\) using the
product of the upper polynomial evaluated at the lower
roots, together with the leading coefficient factor. Each
common simple root contributes exactly one power of \(t\),
with nonzero coefficient equal to its velocity difference
times \(B_N'(\zeta)\sum_i r_i d_i\), up to sign. Each other
root contributes a nonzero limit. This remains true when
a noncommon root is repeated: all roots in its small disk
give the same nonzero limit to their evaluations in the
other parent. Comparing homogeneous parts on generic
\(\mathbf d\) proves the claim, and \(c\) is rational by
specializing a single coordinate of \(\boldsymbol\epsilon\).
It is positive: along \(\epsilon_i=1/M\) for all \(i\) and
positive integers \(M\) tending to infinity, the arrays
come from ordinary moments and the even parent is real
positive. Their resultants are nonnegative by the
conjugate-pair product formula, with no sign on reversing
the roles of consecutive parents.

Let \(C_*\) be the sum of the absolute coefficients of
the terms of \(R\) of total degree greater than \(g\).
It is rational and computable from the exact resultant.
For nonnegative \(\epsilon_i\) with \(0<H_1\le1\), those
terms have absolute sum at most \(C_* H_1^{g+1}\), since
each \(\epsilon_i\le H_1\). Thus if all heights are at least
an integer \(H>\max(D,2D C_*/c)\), the resultant is nonzero.
This explicit bound requires only common-root simplicity
and the stated velocity inequalities. It can be used in
the same recursive procedure as (2), wherever this test
holds at a node with a vanishing comparator resultant.

### Four independently variable ordinary rows

For one retained increasing unit-gap row at height \(A\),
its stripped moments are
\(P_j(z)=\sum_k\binom jk(jA+1)_k z^k\).
The coefficient ratio gives the exact differential identity
\[
             \vartheta P_j
                  =z(j-\vartheta)(jA+1+\vartheta)P_j .
\]
They have nonzero constant terms. They are simple at every
root, by uniqueness for this second-order ordinary differential
equation away from \(z=0\).
For \(B_j(z)=P_j(j^D z)\), the same identity at a root
\(\zeta\) gives
\[
  \frac{((\vartheta^2+\vartheta)B_j)(\zeta)}{B_j'(\zeta)}
                  =(j-jA)\zeta-j^{-D}.
\]
Thus \(V_j(\zeta)=1/(2j^{D+1})-(1-A)\zeta/2\).
At a shared root the difference \(V_n-V_N\) is the
nonzero constant \((n^{-D-1}-N^{-D-1})/2\).

With four unit-gap rows of multiplicity one and empty core,
first apply the disk bound to all four. At least one row is
then fixed in a bounded slice. For the other three, either
the comparator is coprime or the preceding first-order test
applies with the single retained row. This bounds a second
height. The two-height specialization applies to the two
remaining heights with their comparable fixed core.
Doing this at the finitely many terminals below \(J\)
gives the four-height decision.

## One retained unit row with tails of either orientation

**Proposition.** Fix one increasing unit-gap row of height
\(A\ge0\), with multiplicity one, and two consecutive orders \(n,N\).
Append any fixed positive multiplicities \(r_i\) of unit-gap rows
in prescribed orientations \(\epsilon_i\in\{1,-1\}\).
If the signed tail gap \(D=\sum_i\epsilon_i r_i\) differs from \(-1\),
a physical collision requires at least one tail height to lie
below an effective bound. The number of height types is arbitrary.

Set \(\epsilon_i'=1/A_i\) when all tail heights are large.
For each decreasing row strip its scalar \( (j A_i)^{j r_i}\),
and rescale \(z\) by \(\prod_i A_i^{-\epsilon_i r_i}\).
The limiting arrays are \(B_j(z)=P_j(j^D z)\), with \(P_j\)
the one-row polynomials above. Write
\(H_+=\sum_{\epsilon_i=1} r_i\epsilon_i'\) and
\(H_-=\sum_{\epsilon_i=-1} r_i\epsilon_i'\).
The degree-one expansion of the rescaled parent now has the form
\[
 \frac{H_+}{2j}(\vartheta^2+\vartheta)B_j
 +\frac{H_-}{2j}(j-\vartheta)(j-\vartheta+1)B_j .
\]
This follows from the sums \(1+\cdots+k\) and
\(1+\cdots+(j-k)\) on the two orientations.
At a shared root \(\zeta\), the first velocity is the
one already computed. The one-row differential equation
gives the second velocity
\[
 V_j^-(\zeta)=\frac1{2j^{D+1}}+\frac{1+A+2/j}{2}\zeta .
\]
Therefore the difference of the two root velocities is
\[
     \frac{n^{-D-1}-N^{-D-1}}2(H_++H_-)
                  +(1/n-1/N)\zeta H_- .                       \tag{3}
\]
The root \(\zeta\) is nonreal since the even one-row moment
is positive at every real argument, including after dilation.
With \(H=H_++H_->0\) and \(\lambda=H_-/H\in[0,1]\),
the expression (3) divided by \(H\) is nonzero on the entire
closed interval: for \(\lambda>0\) the imaginary part is nonzero,
and at \(\lambda=0\) the real term is nonzero by \(D\ne-1\).

As in the first-order resultant proof, if there are \(g\)
common roots, the homogeneous part of degree \(g\) of the
resultant is a nonzero scalar times the product of the \(g\)
forms (3), one per shared simple root. The other homogeneous
parts start at degree \(g+1\). Divide the leading part by
\(H^g\) and set \(\lambda=H_-/H\). This gives a polynomial
in \(\lambda\) with algebraic coefficients and a strictly
positive minimum of its absolute value on \([0,1]\).
A rational positive lower bound is effective by univariate
root isolation and interval subdivision. The higher terms
are bounded by their absolute coefficient sum times \(H^{g+1}\)
for \(H\le1\), since each \(\epsilon_i'\le H\).
As \(H\le(\sum_i r_i)/\min_i A_i\), the resulting bound is
uniform in the relative tail heights. If the limiting
comparators are coprime the ordinary determinant bound already
suffices. In a mixed recursive slice, subsequent two-height
tests still depend on the orientations of the remaining rows.

## The corresponding Runge condition

Although the slice argument already gives effectivity, these
three-height curves also satisfy a Runge criterion when (1)
is coprime. The saturated complex degree-at-least-two locus
has dimension at most one. Otherwise, project a higher-dimensional
component with a nonconstant height to that coordinate.
Generic fibers in the full-degree chart are positive-dimensional.
In particular, such fibers occur at nonnegative integer
heights, contrary to the two-height specialization theorem.
Any curve component on which a height is constant cannot
contain an ordinary integer point in the full-degree chart:
the constant would then be a nonnegative integer, and the
entire curve would lie in a forbidden two-height slice.

For a geometrically irreducible curve component defined over
\(\mathbb Q\), pass to its smooth projective normalization.
Each height is nonconstant, so each has a pole. No point can
be a pole of all three heights. At such a point the rescaling
preceding (2), along a sequence in the full-degree chart,
would give the coprime limits (1). The saturated rank bound
is closed in the coefficient arrays of those full degrees,
which would force degree at least two for their common factor.
Consequently at least two distinct proper sets of labeled
height coordinates occur as pole sets. These sets are Galois
invariants, so the points at infinity have at least two
Galois orbits.

An integer linear combination of the three height functions
can be chosen with a pole at every point at infinity.
At each such point, merely avoid cancellation of the
largest-order pole, a proper linear condition. Its values
at ordinary integer triples are integers. Runge's theorem
with \(K=\mathbb Q\) and only the archimedean place then
gives effective finiteness on the normalized component.
Rational points on intersections of non-rational conjugate
components and singular affine points are a finite algebraic
exception. [Levin, Theorem 1.2](https://www.numdam.org/item/10.5802/jtnb.634.pdf)

## An exact test of the ordered heights

For the smallest nontrivial algebraic example take the empty
core, three rows with multiplicities one, and orders \(2,3\).
Write \(C(t)=\prod_i(t+A_i)\) and
\[
 u=8C(1/2),\quad v=8C(1),\quad a=27C(1/3),\quad b=27C(2/3).
\]
The stripped lower array is
\(1+2uz+uvz^2\) and the upper is
\(1+3az+3abz^2+(27/8)abv z^3\).
On the full-degree chart a common factor of degree at least
two requires divisibility of the upper by the lower. Its
linear cofactor has constant one and linear coefficient
\(27ab/(8u)\). Comparing coefficients gives
\[
 2u+\frac{27ab}{8u}=3a,\qquad uv=-\frac{15}{4}ab.
\]
In particular this algebraic locus has no ordinary nonnegative
points. Nevertheless its infinity branches illustrate the
difference between the ordered height chart and symmetric
coordinates.

Put \(q=a/u\) and \(D_0=513q^2-1050q-128\).
Cubic interpolation at \(1/3,1/2,2/3,1\), using the two
divisibility identities, gives \(u=-108q/D_0\) and
\[
 6D_0 C(t)=6D_0t^3+A_1t^2+3B_1t+2C_1,
\]
where
\[
 A_1=-5049q^2+11982q+1408,\quad
 B_1=909q^2-2382q-256,\quad
 C_1=-243q^2+642q+64.
\]
The quadratic \(D_0\) is irreducible over \(\mathbb Q\),
with discriminant \(36\cdot13\cdot2917\).
Both \(A_1\) and \(9B_1^2-8A_1C_1\) are coprime to \(D_0\)
already modulo eleven. Hence above each root of \(D_0\)
one height has a simple pole and the other two have
distinct finite limits. This follows by applying the
implicit function theorem to the reciprocal variable
at the vanishing leading coefficient, and directly to
the two simple solutions of the limiting quadratic.
There are six ordered branches at each root of \(D_0\),
one for each ordering of those three roots. There are
no other height poles, since the coefficients of the
monic cubic \(C(t)\) have only those two poles and
are finite at \(q=\infty\).

For clarity these branches belong to a geometrically
irreducible ordered curve. At \(q=1\) the specialized cubic
\(6D_0C(t)\) is, modulo eleven,
\(3t^3+3t^2+5t+2\), which has no root in \(\mathbb F_{11}\).
This gives irreducibility over \(\mathbb Q(q)\) by monic
specialization where \(D_0\ne0\).
The discriminant numerator of the displayed cubic is
of degree eight and squarefree modulo thirteen, so its
discriminant is nonsquare even over \(\overline{\mathbb Q}(q)\);
the denominator for making the cubic monic has even
discriminant exponent. The arithmetic group is therefore
\(S_3\), and its normal geometric subgroup contains an odd
permutation, hence is also \(S_3\). The ordered root cover
is geometrically irreducible. Its branches distinguish
which labeled height has a pole. The symmetric coordinates
alone have just the single pair of conjugate poles of
\(D_0^{-1}\), and therefore one Galois orbit at infinity.
