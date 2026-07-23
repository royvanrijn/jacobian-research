# GMC(2) research program

## Scope

The Gaussian Moments Conjecture is settled negatively in every real Gaussian
dimension \(n\geq 3\): Long's explicit three-variable polynomial

\[
 P=(1+Z)\left(W-\frac{2+Z}{2}T^2\right)
\]

has the five-term expansion

\[
 P=W+WZ-T^2-\frac32ZT^2-\frac12Z^2T^2,
\]

all pure Gaussian moments vanish, and \(Q=Z\) has nonzero mixed moments.
Adding unused independent Gaussian coordinates gives the same conclusion in
every higher dimension.  The exact finite regressions and the all-order
coefficient identity used in Long's proof are transcribed and checked in
finite ranges by
[`verify_long_gaussian_moments.py`](../scripts/verify_long_gaussian_moments.py).

The lower-face prime theorem now proves \(\mathrm{GMC}(2)\).  Together with
the known one-dimensional result and Long's counterexamples, the
dimension-by-dimension classification is therefore:
\[
\operatorname{GMC}(n)\quad\Longleftrightarrow\quad n\leq2.
\]
Accordingly, there is no remaining Gaussian-dimension problem.

The next finite problem is optimization in the first failing dimension:
classify minimal three-real-variable witnesses.  The
[minimal-failures program](THREE_REAL_GMC_MINIMAL_FAILURES.md) proves that
Long's five-term quartic is the unique nontrivial template in the natural
rank-one ansatz \(P=W(1+aZ)+(b_0+b_1Z+b_2Z^2)T^2\).  The first three pure
moments cut out the family scheme-theoretically, and a forced formal square
identity proves all-order vanishing.  Global degree and support minimality
outside that ansatz remain open.

## Resolution of the former frontier

For arbitrary
\[
P=\sum_kT^kB_k(U),
\]
form the radial-order Newton points
\((k,\operatorname{ord}_U B_k)\).  If the rotational support straddles
zero, the lower convex envelope over weight zero has an exposed contact
face whose Laurent polynomial also has zero in its Newton interval.  The
Duistermaat--van der Kallen theorem supplies a nonzero constant term \(c\)
in some power of that face.  At the corresponding prime-dilated moment,
Frobenius and factorial divisibility isolate \(c^p\), contradicting
pure-moment vanishing.  Thus every pure-moment-zero polynomial has strictly
one-sided rotational support, which gives GMC(2).  See the
[lower-face theorem](TWO_REAL_GMC_LOWER_FACE_THEOREM.md).

The route to the theorem produced several independently useful intermediate
results.  The exact three-level family
\(P=WA(U)+C(U)+ZB(U)\) with rotational support \(\{-1,0,1\}\) is now closed
in every degree.  After centering, its Bessel--factorial functional equation
has a two-case prime-endpoint proof: the moments of orders \(p\) and \(2p\)
isolate the two endpoint invariants according as
\(\operatorname{ord}_U D\geq2\operatorname{ord}_U C\) or the reverse
inequality, where \(D=UAB\).  A direct reduction of the finitely generated
coefficient ring replaces the preliminary algebraic-coefficient
specialization.  The earlier 6, 10, and 15 chart calculations in degrees
four, five, and six remain finite-cutoff regressions.  See the
[standalone theorem](PRIME_ENDPOINT_BESSEL_FACTORIAL_RIGIDITY.md) and
[support-graph analysis](TWO_REAL_GMC_SUPPORT_GRAPH_EXPLORATION.md).

The same prime method also closes the first genuine circuit
\(\{-2,-1,1,2\}\).  Frobenius gives
\(\operatorname{CT}(P^{kp})\equiv\operatorname{CT}(P^k)^p\pmod p\).
The toric relation leaves only unique minima, adjacent two-way ties, and a
four-way tie; invariant moments through degrees 12 separate all of them.
There is no surviving valuation cone.

It also closes every unit star
\(\{0,1,-d_1,\ldots,-d_q\}\) and its reflection.  The zero-weight relation
monoid is free on the leaf invariants, and minimum normalized \(U\)-order
plus prime dilation isolates one of them.  This includes the smallest star
\(\{-2,-1,0,1\}\) in all degrees.

The unit-star proof also identifies why literal graph leaf removal is not
the final formulation: arbitrary weights have mixed Hilbert-basis relations
such as \(5-2-3=0\).  The lower-face theorem packages all such relations at
once, so no arbitrary-star or multi-cycle case remains.

The differential systems remain structural interpretations rather than an
active proof route.
For \(G=((1-tC)^2-4t^2D)^{-1/2}\), integration by parts reduces every radial
moment \(\mathcal L(U^kG)\) to finitely many initial radial moments and turns
the two algebraic differential identities for \(G\) into a finite
meromorphic Pfaffian system.  In the three-level case the prime theorem says
that its distinguished factorial initial-value solution cannot remain in
the hyperplane \(\mathcal L(G)=1\).  The support-graph analysis gives the
exact recurrence, a degree-\((2,3)\) regression, and the associated
Laplace--Bessel and Laguerre interpretations for possible reuse.

## Retained high-dimensional material

High-dimensional results remain in the active reference corpus when their
purpose is logical transport rather than witness discovery:

- the fixed-dimensional
  \(\mathrm{GMC}(2r)\Rightarrow\mathrm{SIC}(r)\Rightarrow\) Keller
  invertibility implication;
- the BCW consequence chains to failures of GMC(42) and GMC(158);
- the four-real weighted-seed bridge and its transport of seed geometry into
  moment coordinates; and
- the associated Image/Vanishing consequences.

These are examples of how structure moves between conjectures and
categories.  They are not active attempts to optimize a Gaussian witness
dimension, because the five-term witness already settles every \(n\geq3\).

## Archive

The former three-real weighted-family searches are preserved under
[`archive/high-dimensional-gmc`](../archive/high-dimensional-gmc/README.md).
Their exact finite exclusions remain reproducible historical data, but their
surviving ansatzes are not an open research queue.
