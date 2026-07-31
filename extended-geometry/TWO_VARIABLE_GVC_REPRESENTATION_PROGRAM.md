# Two-variable GVC and rank-efficient Laplacian polarization

## 1. Status

The nonhomogeneous termination program has been reduced to one exact
promotion problem.  The
[finite-trace digit-separation theorem](BINARY_GVC_UNIFORM_FACE_TERMINATION.md)
splits every scale-compatible finite-character trace, but the Hall--jet
filtration has not yet been proved to promote each affine,
prime-dependent carry shell to such a trace.  Consequently unrestricted
\(\operatorname{GVC}(2)\) remains open.  The earlier balanced
homogeneous target is proved independently.  The
[balanced cubic theorem](TWO_VARIABLE_CUBIC_GVC_THEOREM.md) proves that
the first four pure moments force every cubic Segre point into the
one-sided nullcone, with the explicit mixed cutoff \(m>\deg Q\).
More generally, the
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) proves the GVC
conclusion in every balanced homogeneous degree: every binary operator
symbol splits into linear factors, and translated complete polarization
converts its pure values to constant terms of powers of one Laurent
polynomial.  Retaining the translation variable removes the earlier
degree restriction: homogeneous binary operators satisfy GVC for
arbitrary \(P\).

The later
[separable escape obstruction](SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
closes a further nonhomogeneous class.  If \(r\) is the lowest positive
homogeneous order of a binary constant-coefficient operator \(\Lambda\),
then GVC holds for every \(P\) of degree at most \(r\).  It also proves
GVC for arbitrary \(P\) when \(\Lambda\) is a homogeneous split factor
times a differential operator with nonzero constant term.  It further proves
that separated multiplicative auxiliary specialization, rank-one
dilation, and nonlinear Segre substitution cannot convert the new
two-pair witness: the output remains rank one, whereas the witness has
rank five.  Coefficient extraction can sum rank-one channels only by
losing multiplicativity, so it supplies no automatic all-order transfer.
Moreover, fixed linear translation followed by one diagonal coefficient
represents exactly a product of powers of linear symbols.  The first
remaining operator architectures must therefore use coupled coefficients,
nonlinear translation, or an irreducible nonhomogeneous symbol.

There are now two continuing targets: prove scale-compatible carry
promotion—most sharply, prove translation-curvature rigidity for every
rank-one Cartesian isoperiodic cluster—or find a genuine binary
counterexample, and replace the third-order Dvorsky operator by a quadratic
constant-coefficient operator using as few auxiliary variables as
possible.  Polynomial degree seven is a regression for the first target.

The binary partial results do not change the known counterexample upper
bounds: they remain five variables for unrestricted constant-coefficient
GVC and forty variables for ordinary-Laplacian GVC.  They do not yet
establish a counterexample-free range through two variables.
The separate
[three-variable tagged-lift analysis](THREE_VARIABLE_GVC_TAGGED_LIFT.md)
tests whether the two-pair Image-Mathieu witness can lower the former
bound.  It derives an exact power-compatible channel formula but excludes
the minimal Long-tag architecture at pure moment five over
\(\mathbb Q\).  Its three-channel continuation also closes every
rank-three Dvorsky parallelogram on the first positive
order/degree-\((2,3,4)\) grading and proves the smallest persistent
endpoint repairs terminal.  Its complete odd-quartic chart has
moment-six radical \((A,S,RU)\); both components have all-order mixed
cutoffs, including elimination of the degree-two operator endpoint at
every field-valued pure-zero point.  It neither proves nor disproves
GVC(3).

The first target has an additional associated-graded gate.  For arbitrary
nonhomogeneous \(\Lambda\), every fixed number of leading homogeneous
layers of \(\Lambda^m(QP^m)\) vanishes eventually under the pure premise.
Consequently a counterexample cannot live on a fixed leading face: its
mixed defect must move to unbounded depth as \(m\) grows.
If the operator has only linear and quadratic pieces and its linear part
is nonzero, the first two moments close arbitrary \(P\).  In normalized
coordinates, the unique highest \(y\)-degree term of the second moment is
the square contribution from \(C^2P_{yy}^2\), forcing \(P=ay+b\).
The same product-defect argument closes every separated drift
\(\partial_x+h(\partial_y)\), even for formal \(h\), since derivative
series act locally finitely on polynomials.  Formal Weierstrass division
factors every binary symbol with nonzero linear part as
\(U(\xi,\eta)(\xi+q(\eta))\), where \(U(\partial)\) is a locally finite
differential automorphism.  This closes the entire lowest-order-one
frontier for arbitrary \(P\).  Exact cubic and quartic jet computations
give independent finite regressions of the factorization argument.
Before the uniform termination theorem, this reduced any possible
GVC(2) counterexample to lowest positive operator order at least two.
The first three pure equations also close
the entire \(r=2,\deg P\le3\) cell: the double-line quadratic orbit has
one apparent second-moment cancellation, killed exactly at moment three,
and every surviving branch in both quadratic orbits has a strict weighted
degree separator.  Thus the first \(r=2\) target has \(\deg P\ge4\).
The quartic analysis below closes that endpoint as well.
In that quartic target, the distinct-root quadratic orbit and the
double-line \(xy^3\) top form close through moment three and weighted
cutoffs.  Only the double-line \(y^4\) top form remains; eliminating one
coefficient from its second and third moments leaves finitely many ratios
on an explicit sextic.  Its fourth-moment octavic is coprime to that
sextic.  Moreover, weights
\(\operatorname{wt}(x)=2,\operatorname{wt}(y)=1\) show that the full pure
moments equal those of this three-term weighted face: every higher
operator jet has excess weight and every lower polynomial term has
deficient weight.  Hence the entire \(r=2,\deg P\le4\) cell is closed,
and the first quadratic-leading target has degree at least five.
For \(r=3,\deg P=4\), the triple-root cubic-symbol orbit closes through
moment four and three weighted-face reductions.  The double-root orbit
closes through moment three.  For the squarefree orbit, the first two
leading moments leave three fourth-power tips in one root-permutation
orbit; its \(x^4\) tip has the same weight-four correction face as the
double-root case and dies at moment three.  Hence every binary
constant-coefficient operator is GVC-safe for \(\deg P\le4\).  A
counterexample must have degree at least five; initially its lowest
positive operator order could only be \(2,3,\) or \(4\).  For \(r=2\),
the first two
leading equations leave only
\((\partial_x\partial_y,x^5)\),
\((\partial_x^2,xy^4)\), and
\((\partial_x^2,y^5)\).  The full correction calculation closes all
three: the distinct-root tip has a linear-factor/high-order cutoff, the
\(xy^4\) tip dies on two nested weighted faces, and the \(y^5\) tip
reduces to six ratios, four killed at moment three and two already
one-sided.  Hence the first quadratic-leading target has degree at least
six, and a degree-five counterexample must have \(r=3\) or \(4\).
For \(r=3\), the first two leading equations leave eight top-form normal
forms: four triple-root, three double-root, and one squarefree.  The
[binary degree-five frontier theorem](BINARY_DEGREE_FIVE_GVC_FRONTIER.md)
closes all eight correction systems.  It also closes the \(r=4\)
squarefree-quartic row uniformly in its cross-ratio.  The
[quadruple-root quartic theorem](BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md)
now closes the partition \((4)\), including arbitrary lower symbol terms
and arbitrary higher operator jets.  The
[triple-plus-simple theorem](BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md)
and [double-root theorem](BINARY_QUARTIC_DOUBLE_ROOT_GVC.md)
close the other three quartic root partitions.  Thus every binary
constant-coefficient operator satisfies GVC through polynomial degree
five; the next nonhomogeneous degree frontier begins at six.  The
[complete quintic-leading sextic theorem](BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md)
now closes the entire \((r,\deg P)=(5,6)\) row.  A Hall-matching argument
classifies the leading locus for all seven quintic root partitions, and
local exact calculations close root multiplicities one through five with
arbitrary higher jets.  The longest defect chain, on the quintuple-root
cell, reaches operator order ten and dies at pure moment six.  The
quartic-leading row is closed by
[GVC2D6R4](BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md), and the
[complete cubic-leading sextic theorem](BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md)
closes \(r=3\) for all three cubic root partitions.  Finally, the
[complete quadratic-leading sextic theorem](BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md)
closes both quadratic root partitions.  Its double-line pure-sixth-power
endpoint branches through the exact slope-two radical
\((r,q,A,Bz)\), and both axes terminate on common-threshold
coordinate-deficit faces.  Thus every binary constant-coefficient operator
satisfies GVC through polynomial degree six, and the next genuinely
nonhomogeneous polynomial-degree frontier begins at seven.
The
[uniform face-termination theorem](BINARY_GVC_UNIFORM_FACE_TERMINATION.md)
now removes two previously degree-specific parts of this analysis.  For
arbitrary \(r<d\), every leading Hall component is localized at a single
root of multiplicity \(\mu\) and has
\[
 (X^\mu Y^{r-\mu},\,y^{d-\mu+1}C_{\mu-1}).
\]
Moreover, every unequal-weight common-threshold face is automatically
terminal: a prime-dilated least endpoint proves that its operator and
polynomial Newton segments are disjoint, hence one coordinate has a
linear derivative deficit.  The shifted version excludes one
unequal-weight exposed endpoint on every fixed rational normalized output
ray.  Exact terminal radicals are therefore no longer part of the
all-degree problem.

What remains is narrower but genuine.  At depth \(\lfloor\rho m\rfloor\),
coefficient extraction can couple several selection distributions with
positive limiting frequency.  The split-symbol separator also shows that
\(o(m)\) excursions from one ordinary-homogeneous face cannot repair its
linear gap.  Hence a surviving coefficient must use at least two distinct
faces with positive limiting density.  It is not the power of one fixed
rank-one face.  Liu--Sun homogeneous binary factorial rigidity removes
every one-channel slope-\(-1\) Pareto tie in all degrees, but it does not
apply after a genuine Hall selection constant term: already for
\(G=S U_1+S^{-1}U_2\),
\[
 \operatorname{CT}_S G=0,\qquad
 \operatorname{CT}_S G^2=2U_1U_2.
\]
After the beta integral, the exact residual is a beta--torus moment
\[
 \int_0^1\operatorname{CT}_S
       G(S,t,1-t)^m\,dt,
\]
not the moment of \((\operatorname{CT}_S G)^m\).  A complete termination
proof cannot assert that the unrestricted beta--torus kernel is
Mathieu--Zhao: Long's circuit
\[
 (1-S^{-1})((1-t)+tS)
\]
has every pure moment zero and a nonzero \(S^{-1}\)-mixed moment for every
power.  In the GVC realization it is only the linear Hall annihilator
\[
 (b\partial_x+a\partial_y)(dx+cy)=bd+ac=0;
\]
the bad Laurent character is not a polynomial multiplier.  Thus the
correct restricted theorem must show that the convolution has an exposed
one-radial lower face, descends to smaller Hall/jet support, or acquires a
split-symbol separator whose gap controls the polynomial-multiplier cone.
The minimal two-by-two Bernstein circuit is proved to have exactly this
last form in the uniform termination note.  That note now also closes the
first genuine higher-jet circuit.  For primitive endpoints
\((r,0),(0,s)\), every nondegenerate four-channel parallelogram has a
third-moment obstruction
\[
 T_r-T_s+\frac92(C_s-C_r)(C_r+C_s-2)\ne0
 \qquad(r\ne s).
\]
Together with the centered three-level theorem and a moment-two sparse
return argument, this closes every two-operator-endpoint/four-polynomial-
channel circuit.  A survivor therefore needs at least five polynomial
channels or at least three operator endpoints.  Finite-prefix inheritance
is false: an exact five-channel pair has its first three scalar moments
zero while every four-channel minor fails, and only its fourth moment
provides the Bernstein pivot.  Exact saturation through moment four now
closes every support with at most five polynomial channels for
\((r,s)=(1,2),(1,3),(1,4),(2,3)\); only 181 of 403,960 support ideals
require the fourth pivot.  One mixed fourth-pivot graph is now uniform:
the half-bridge support at every even \((1,n)\) has an explicit nonzero
central-multinomial fourth moment.  Canonical balance rows compress all
181 fourth-pivot supports to 14 return matrices, already visible at
\((1,2)\), and primitive-return combinatorics proves that list exhaustive
in the reduced early-entry regime.  The coefficient equations pair these
into only eight factorial-ratio obstructions, with no exact zero through
unequal endpoint orders 100.  The transverse return equation proves that
all order-four rows are generated by the primitive quadratic and cubic
rows, so there is no special primitive order-four branch.  Both
double-quadratic determinants are uniformly nonzero, and an exact
coefficient-positive successor certificate closes
\(\mathcal H_{0,0}\) for every unequal endpoint pair.  Coupled endpoint-
ratio cones also close \(\mathcal H_{0,3}\) and
\(\mathcal D_{0,2}\) uniformly.  Monotone ordered-tail cones and fixed
endpoint rays close \(\mathcal H_{0,1},\mathcal H_{0,2}\), and
\(\mathcal H_{1,1}\) as well.  On the final wedge \(r>s\ge4\), the
increasing factorial product \(L_nM_n/C_n^3\) gives a nonlinear cone
that closes \(\mathcal H_{1,0}\) and \(\mathcal D_{0,1}\).  Thus the
eight-obstruction arithmetic residue is complete.  Moreover, a channel
entering after a core witness cannot change that witness; cofinite
prime/Bessel witnesses give an explicit entry-order cutoff.  The missing
descent must now pass from whole-convolution vanishing to this closed
two-endpoint list.  The normalized return cone sharpens that target:
every extreme return circuit uses at most four total operator/polynomial
variables; two incomparable circuits over two operator endpoints use at
most six polynomial channels, with equality only for two disjoint
centered triples attached to opposite endpoints.  Prime dilation
separates unequal radial-degree profiles.  The maximal opposite
three-by-three packet is terminal once exposed: its normalized series
factors, and the coefficients in degrees \(d\) and \(2d\) contradict
all-order vanishing by strict central-binomial supermultiplicativity.
The remaining theorem must expose a pure-zero proper packet inside the
equal-profile convolution or produce support loss or a split-symbol
separator, also when three or more operator endpoints are active.  The
same-prime filtration now has the exact score “high-digit jet excess
plus binary carry penalty”; the latter is nonnegative and vanishes only
on Frobenius returns.  On an ordinary-homogeneous radial face the first
correction has only centered triples and two-by-two
beta/parallelograms as its inclusion-minimal low-digit support patterns.
The centered ghost diagonal has universal factor \(X-1\), and the beta
diagonal has universal factors \(X(X+1)\), recovering the known terminal
relations.  Additional roots occur at individual primes, so fixed-prime
triangularity is false.  The beta diagonal also has the persistent
algebraic factor \(X^2+X+1\) at every prime at least five, so single-row
cross-prime avoidance is false.  After adjoining the ordinary rows,
however, the isolated atom blocks are terminal: \(1+X\) reduces the
beta block to the Hall value \(X=-1\), while the centered Bessel rows
\(U,U^2+2V\) have determinant \(2\) and force support loss.  The
outstanding step is the global compatibility theorem showing that
these pivots occur with the same normalized high-digit quotient, in
triangular order, inside the full ghost shell and its next unit
correction.  A circuit-only quotient reduction is false after coordinate
projection: the first primitive is
\(R_3B_1B_2=R_0B_3^2\), with support five.  Restoring \(R_2\) or the
reversed level \(B_0\) gives a two-step circuit path; if neither is
active, radial digit separation exposes the two-state block, so this
first projected obstruction is terminal.  Rational-normal-scroll
Gröbner bounds make the remaining quotient list finite for every fixed
support.  Every whole exposed color-count/radial-profile fiber is also
terminal by a one-Laurent-variable constant-term reduction and the
Duistermaat--van der Kallen theorem.  In particular, the first completed
scrolls not controlled by Graver-equals-universal-Gröbner visibility,
\(S(6)\) and \(S(5,4)\), create no new terminality obstruction after
their complete profiles are exposed.  The normalized prime-power rows
of every exposed profile form an integral \(p\)-typical Witt ghost
sequence in all heights: both the Laurent constant-term sequence and
the signed radial factorial unit satisfy Gauss congruences.  Thus the
next-prime-power correction does not cause migration.  Ghost injectivity
alone does not split a vanishing sum of profile Witt vectors, but
repeated equal base-\(p\) digits do: the rows
\(m(1+p+\cdots+p^{k-1})\) recover the first \(q\) power sums of any
\(q\) finite trace components, and Newton identities split them
componentwise.
After an entire oriented radial vector is
exposed, summing over every operator and polynomial color count gives
\[
 [X^{N\rho}]\lambda(X)^{Nd}[Y^{N\rho}]P(Y)^{Nd}
 =
 \operatorname {CT}_{X,Y}
 \bigl(X^{-\rho}Y^{-\rho}\lambda(X)^dP(Y)^d\bigr)^N.
\]
Hence Duistermaat--van der Kallen closes the complete radial union
without any color-count idempotents.  Achievable color counts need not
be saturated—the minimal parity family has persistent odd holes—so the
inheritance problem is precisely exposure of the complete
oriented union rather than a carry-selected proper subset, together with
the coordinate-reversed tie.
For the first reversal width, the tie is also terminal: exact saturation
of the four endpoint charts for
\(G(z)=\sum_{-2}^{2}c_kz^k\) closes the opposite coefficient sums by
rows \(2,4,4,8\).  In general, reversal and a scale-compatible carry
class \(\pi(p,q)=Nc\) both become identity-coefficient problems in
\(\mathbb C[\mathbb Z^r\times F]\), where the finite abelian factor
records reversal or carry characters.
The exact reduction uses the regular representation of the finite
character group:
\[
 \operatorname {CT}_{\mathbb Z^r\times F}(u^N)
 =
 |F|^{-1}\operatorname {CT}_{\mathbb Z^r}
 \operatorname {Tr}(\operatorname {Reg}_F(u)^N).
\]
Equivalently, the full moment series is the constant logarithmic
derivative of
\(\det(I-t\operatorname {Reg}_F(u))\).  For primes congruent to one
modulo \(\exp(F)\), Frobenius fixes every torsion character.  The
finite-trace digit theorem is stronger than either proposed attack: it
splits arbitrary Laurent character components, not only the rank-one
Cartesian ones.  Componentwise Duistermaat--van der Kallen separation
then closes every packet already promoted to a fixed scale-compatible
trace.

What remains is the promotion itself.  A fixed affine residue has
Fourier weights outside the \(N\)-th power.  Mixed repeated digits show
that such a weighted trace cancels exactly in isoperiodic character
clusters; the dilation pair
\(z+z^{-1}\), \(z^2+z^{-2}\) shows that this mechanism can retain a
mixed row.  The radial factorial also blocks a direct shortcut:
for \(C=y^2+4xy+2x^2\), the \(p=11\), \(N=12\) factorial moment has
valuation three rather than the naive digit-factorized valuation two.
Thus the remaining theorem must show that every proportional-depth
Hall shell becomes scale-compatible, is Hall-terminal, or loses
support.  Generic translation makes this more precise.  The Taylor
channels satisfy
\(\partial_{z_i}p_\beta=(\beta_i+1)p_{\beta+e_i}\).
If the ratio of two character twists is constant along each Taylor
direction, it is a scalar--torus character, so its moving rows differ by
an \(N\)-th-power phase and are already scale-compatible.  Otherwise a
minimal collinear triple or unit square has nontrivial multiplicative
curvature.  The immediate target is to expose that curvature block by
differentiating the all-order period identity while preserving the
common high-digit quotient; the isolated augmented Bessel and beta
blocks are already terminal.  Exact \(C_2\) searches find no
non-symmetry collision through degree 12 in one Taylor direction or on
the complete \((2,2)\) and \((3,2)\) binomial Taylor rectangles, but
these are bounded experiments.
Equivalently, generic translation identifies the balanced radial polytope
with
\[
 \operatorname{Newt}(\lambda)\cap
 \operatorname{conv}\{\beta:\partial^\beta P\ne0\}.
\]
The multiradial prime theorem closes an empty intersection and every
intersection with a componentwise least point.  A remaining pair must
therefore mix at least two incomparable Pareto-minimal endpoints.  The
degree-at-most-six calculations verify that alternative only in their
finite support ranges.  A toric blow-up does not close the gap by itself:
although it makes incomparable exponents comparable, it changes \(a!b!\)
by an exponent-dependent factorial ratio.  The required Hall reduction
must therefore be factorial-compatible.
Although a squarefree tangent cone formally factors into smooth branches,
this does not transport powers: a multiplicative conjugate of a
derivation is a derivation, so a constant-coefficient umbral conjugate
must have linear symbol.  Exploiting formal branches therefore requires a
new multilinear identity, not an algebra-coordinate straightening.

The later
[two-pair SIC counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
lies in bidegree \((4,4)\), but its \(5\times5\) coefficient matrix has
determinant \(48\).  It is therefore a full-rank nonseparable point of
\(\operatorname{End}(\operatorname{Sym}^4)\), not a rank-one Segre point
\(A\otimes P\).  It falsifies unrestricted SIC(2) without changing either
GVC target above.

## 2. Balanced homogeneous GVC is a Segre slice of SIC(2)

Let \(k\) be a characteristic-zero field and let \(U\) be two-dimensional.
A homogeneous constant-coefficient operator of order \(d\) has a symbol

\[
 A\in\operatorname{Sym}^d U,
 \qquad
 \Lambda=A(\partial),
 \]

while a homogeneous polynomial of degree \(d\) is an element

\[
 P\in\operatorname{Sym}^d(U^*).
 \]

Introduce dual variables \(\zeta\) and coordinate variables \(z\), and put

\[
 f(\zeta,z)=A(\zeta)P(z).
 \tag{2.1}
\]

Then \(f\) is a rank-one point of

\[
 V_d
 =\operatorname{Sym}^d U\otimes\operatorname{Sym}^d(U^*)
 =\operatorname{End}(\operatorname{Sym}^d U).
 \tag{2.2}
\]

If \(\mathcal E_2\) is the two-pair contraction, the definitions give the
exact identities

\[
 \boxed{
 \mathcal E_2(f^m)=\Lambda^m(P^m),\qquad
 \mathcal E_2(Q(z)f^m)=\Lambda^m(QP^m).
 }
 \tag{2.3}
\]

Thus balanced homogeneous two-variable GVC is not merely analogous to the
SIC(2) problem.  It is the restriction of the SIC moments to the affine
Segre cone

\[
 \Sigma_d
 =\{A\otimes P\}\subset V_d.
 \tag{2.4}
\]

The multiplier class in GVC is smaller: \(Q\) depends only on the
coordinate variables, whereas SIC allows an arbitrary multiplier in both
sets of variables.  Consequently a moment--nullcone theorem on all of
\(V_d\) proves the balanced homogeneous GVC conclusion, but the restricted
GVC problem on \(\Sigma_d\) may be strictly easier.

After self-dualizing the binary representation, Clebsch--Gordan gives

\[
 A\otimes P
 \longmapsto
 (F_0,F_2,\ldots,F_{2d})
 \in\bigoplus_{r=0}^d\operatorname{Sym}^{2r}U.
 \tag{2.5}
\]

The components \(F_{2r}\) are the joint transvectants of \(A\) and \(P\),
and the quantities in (2.3) are their invariant contractions.  The
one-sided nullcone criterion from the
[two-pair moment--nullcone program](TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
becomes

\[
 F_0=0,\qquad
 F_{2r}=L^{r+1}G_{r-1}\quad(1\leq r\leq d)
 \tag{2.6}
\]

for one common linear form \(L\).  This supplies the same two-stage proof
architecture:

1. force the lowest nonzero transvectant into its binary-form nullcone; and
2. synchronize the unique high-multiplicity root through every higher
   transvectant.

The additional rank-one equations defining \(\Sigma_d\) should be imposed
before elimination.  Eliminating in the full \((d+1)^2\)-dimensional SIC
space discards the main advantage of the GVC specialization.

## 3. The cubic Segre slice is closed

The complete bidegree-\((2,2)\) theorem proves SIC(2) on all of \(V_2\), so
it proves the balanced quadratic GVC conclusion a fortiori.  It also shows
that any bidegree-preserving linear pair compression of the compact
three-pair witness to two pairs cannot work: linear pair substitution
preserves bidegree \((2,2)\), and the resulting form lies in the proven-safe
stratum.  This statement does not cover a linear substitution that mixes
dual and coordinate variables and thereby destroys the bigrading.

The next case was the cubic Segre slice.  It now satisfies the exact
equality

\[
 \boxed{
 \Sigma_3\cap
 V(\mu_1,\mu_2,\mu_3,\mu_4)
 =
 \Sigma_3\cap N_3.
 }
 \tag{3.1}
\]

Indeed, normalize the nonzero cubic symbol to its triple-root,
double-root, or squarefree orbit.  The first orbit closes after moment one,
the second after moment three, and the squarefree orbit after moment four.
In the squarefree orbit the surviving polynomial is one of three pure
cubes, each annihilated by one linear factor of the operator.  Every branch
then gives

\[
 \Lambda^m(QP^m)=0\qquad(m>\deg Q).
 \tag{3.2}
\]

The full proof and exact moment identities are in
[`TWO_VARIABLE_CUBIC_GVC_THEOREM.md`](TWO_VARIABLE_CUBIC_GVC_THEOREM.md).
Since \(\Sigma_3\) has dimension seven inside the sixteen-dimensional
\(V_3\), this also proves that the extra semistable component of the first
thirteen full \((3,3)\) moments does not meet the Segre cone.  The first
open balanced homogeneous separable case is now degree four.

The degree-four frontier has since narrowed further.  The
[low-root theorems](TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md) prove every
balanced degree when the operator symbol has at most two distinct roots,
using the one-variable Laurent constant-term theorem.  In degree four,
moments through order five also cut out the one-sided locus on the
\((2,1,1)\) three-root orbit, with five exact fourth-power certificates.
On a nonempty Zariski-open set of squarefree quartic cross-ratios, moments
through order six cut out the four annihilator lines; this follows from an
exact reduced four-point fiber at cross-ratio \(2\) and proper-family upper
semicontinuity.
The later [split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) closes the
GVC conclusion on the one-parameter squarefree quartic orbit and in every
higher balanced degree.  What remains on that quartic orbit is the stronger
finite-moment nullcone question, not eventual coordinate-multiplier
vanishing.

## 4. What the SIC(2) counterexample changes

The bidegree-\((2,2)\) theorem closes more than the natural four-parameter
compression family: it closes every two-pair form in that stratum.
The successful witness uses balanced bidegree \((4,4)\) and is genuinely
nonseparable.  Thus the earlier list of possible escape mechanisms is now
resolved as follows:

1. bidegree \((2,2)\) remains completely safe;
2. full bidegree \((3,3)\) remains undecided;
3. bidegree \((4,4)\) already fails at a full-rank nonseparable point; and
4. unequal/nonhomogeneous degree and nonlinear polarization are no longer
   needed to falsify SIC(2), though they remain possible sources of smaller
   support or degree.

A single unequal bihomogeneous component is already one-sided by its
nonzero central weight.  The viable unequal-degree option is therefore a
nonhomogeneous mixture whose different central weights can interact in
powers; unequal bidegree by itself is not an open counterexample stratum.

For ordinary homogeneous GVC, (2.1) is necessarily separable.  The new
counterexample therefore belongs only to the larger SIC problem, while the
Segre restriction remains the structural advantage available to GVC.

Moreover, the rank obstruction is stable under every separated
multiplicative conversion.  Evaluation after auxiliary polarization,
separated restriction of a rank-one dilation, and nonlinear substitutions
on the dual and coordinate factors all preserve the product form.  They
therefore cannot produce the rank-five matrix (1.6).  A sum needs at least
five rank-one channels.  Auxiliary coefficient extraction is not
multiplicative, so a construction using it must prove new identities for
every power rather than inherit the witness identities by specialization.

For a nonhomogeneous operator
\(\Lambda=\Lambda_r+\cdots+\Lambda_d\), the lowest-order theorem further
shows that \(\deg P\leq r\) is safe.  In the pure premise, all higher
operator pieces are killed by degree and one obtains
\(\Lambda_r^m(P_r^m)=0\).  In the mixed expression only a bounded number
of higher-order pieces can occur; the split-symbol Newton gap for
\(\Lambda_r\) absorbs these bounded defects.  Thus nonhomogeneity alone is
not an escape mechanism.

## 5. Rank-efficient quadraticization of the Dvorsky operator

The Dvorsky--Long witness uses

\[
 \Lambda
 =\partial_t(\partial_a\partial_d-\partial_b\partial_c),
 \qquad
 P=(t+c)(ad+bt),
 \tag{5.1}
\]

so its operator symbol is the factored cubic

\[
 \sigma_\Lambda
 =\xi_t(\xi_a\xi_d-\xi_b\xi_c).
 \tag{5.2}
\]

This does not directly give ordinary-Laplacian GVC, whose operator symbol
is a nondegenerate quadratic form.  The promising conversion problem is to
find a polynomial lift to \(N\) variables with a nondegenerate quadratic
operator \(\widetilde\Delta\), a polynomial \(\widetilde P\), and a fixed
multiplier \(\widetilde Q\), together with nonzero constants \(c_m,d_m\),
such that the all-order identities

\[
\begin{aligned}
 \widetilde\Delta^m(\widetilde P^m)
   &=c_m\,\iota\!\left(\Lambda^m(P^m)\right),\\
 \widetilde\Delta^m(\widetilde Q\,\widetilde P^m)
   &=d_m\,\iota\!\left(\Lambda^m((-c)P^m)\right)
\end{aligned}
\tag{5.3}
\]

hold for every \(m\), for an explicit injective polynomial substitution
\(\iota\).  After a linear change over an algebraic closure, a
nondegenerate quadratic constant-coefficient operator is an ordinary
Laplacian.

The optimization objective is the auxiliary rank, hence \(N\), rather than
the size of a generic polarization.  The factorization (5.2) should be
retained as a \(1+2\) tensor network: the determinant factor has rank two
and should be quadraticized as one block, not expanded and polarized
monomial by monomial.

The useful machinery from the \(HC_4\) program is methodological:

1. choose a nondegenerate quadratic pivot;
2. eliminate auxiliary blocks by a Schur complement;
3. stratify the residual coupling by rank and constant-kernel dimension;
4. use source/dual bidegrees to isolate faces that cannot cancel; and
5. synchronize any surviving moving kernel before increasing the ansatz.

A constant Schur complement of quadratic symbols is still quadratic, so it
cannot by itself reproduce the cubic operator.  Any successful
construction must act jointly on the operator and polynomial sides, or use
a nonlinear polarization.  This is the first structural gate.

The first calculation should therefore solve the universal coefficient
problem for one auxiliary quadratic block while preserving both lines of
(5.3).  That calculation is now complete for every polynomial or formal
hyperplane lift of the canonical one-pair completion.  With one new
variable \(s\), consider the
nondegenerate completion
\[
 \widetilde\Delta
 =\partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s.
\tag{5.4}
\]
For every polynomial or formal \(F\) restricting to \(P\) on \(s=0\), the
first pure equation \(\widetilde\Delta F=0\) forces
\[
 \left.
 \widetilde\Delta^2(F^2)
 \right|_{a=b=c=d=s=0}
 =12t^2-8\rho t
\tag{5.5}
\]
for one transverse first-jet coefficient \(\rho\).  Thus every such
six-variable hyperplane lift fails at moment two.  The proof and exact
certificate are in
[the one-pair obstruction](DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md).
Degree mixing cannot change the coefficient \(12\).

The next ansatz must therefore add a second auxiliary block, replace the
hyperplane restriction by a nonlinear polarization, or use a different
quadratic completion with genuinely new cross terms.
Only after all-order identities are derived should the block rank be
minimized.  A finite moment match, or preservation of the pure identity
without the mixed defect, does not lower the ordinary-Laplacian endpoint.

## 6. Laplacian promotion rule

The cubic Segre theorem and the one-pair obstruction have exact checkers.
No successful quadratic ordinary-Laplacian lift is currently attached to
the second branch.  Such a result may change the Laplacian status ledger
only after it supplies:

1. an explicit characteristic-zero construction;
2. an all-order proof of the pure vanishing;
3. an all-order or infinite-tail proof of the fixed-multiplier defect;
4. a proof that the quadratic operator is nondegenerate and hence linearly
   equivalent to the ordinary Laplacian; and
5. the exact number of variables and auxiliary-block rank.
