# Complete-product residue images at arbitrary prime-power depths

For a common target \(C\), the two positive complementary products
fill their natural joint residue range at every sufficiently large
prime, once there are at least three slots. This includes primes
that do not divide \(C\). The joint image over an arbitrary
composite modulus retains a fixed positive fraction of the
natural range. When the modulus is small relative to specified
intervals, every occupied pair is realized by actual distinct
positive primitive rows in those intervals.

## Joint image and interval statement

Fix \(k\ge3\), a positive target \(C\ge2\), and a modulus \(Q\ge1\).
Consider
\[
   U=\prod_{i=1}^k a_i,\qquad V=\prod_{i=1}^k(C-a_i)\pmod Q,
\]
with each \(a_i\) a unit at the primes common to \(Q\) and \(C\).
Put \(D=\gcd(Q,C)\) and
\[
                      Q_T=\prod_{\substack{p^e\parallel Q\\p\mid C}}p^e .
\]
The natural range consists of pairs with \(U,V\) units at \(Q_T\)
and \(V\equiv(-1)^k U\bmod D\). The image lies in this range
by multiplying the congruences \(C-a_i\equiv-a_i\bmod D\).
The natural range has cardinality
\[
                       \frac{\varphi(Q_T)}{Q_T}\frac{Q^2}{D}.
\]

**CR1.** At each sufficiently large prime \(p\), the image
modulo \(p^e\) equals its natural range for every depth \(e\)
and every \(C\). The prime threshold can be taken absolute.
For arbitrary \(Q\) the image has cardinality at least
\[
                    c_k\frac{\varphi(Q_T)}{Q_T}\frac{Q^2}{D}
                                                               \tag{CR1}
\]
for a constant \(c_k>0\) depending only on \(k\).
At \(p\nmid C\) the local natural range is all residue pairs,
including zero or nonunit coordinates.

**CR2.** Fix \(\eta>0\) and prescribe an interval \(I_i\subseteq[1,C-1]\)
of length at least \(Q C^\eta\) for each slot. For all sufficiently
large \(C\), every pair in the image in (CR1) can be realized
with \(A_i\in I_i\), \(B_i=C-A_i\), \((A_i,B_i)=1\), and the
\(A_i\) pairwise distinct. The statement is uniform in \(Q\)
and the intervals.

## Primes away from the target

For \(k=3\) and \(p\nmid C\), first take unit values of \(U,V\)
modulo \(p\). Scale the inputs by \(C\), so the prescribed
outputs become \(u=U/C^3,\ v=V/C^3\). The equations are
\(xyz=u\) and \((1-x)(1-y)(1-z)=v\).
Eliminating \(z\) gives the plane curve
\[
                  (1-x)(1-y)(xy-u)-vxy=0.
\]
It has a rational smooth point at \((0,1)\), where the
\(y\)-derivative equals \(u\). The unique geometric component
through this point has multiplicity one and is defined over
\(\mathbb F_p\), because Frobenius must fix that component.
The component is not contained in a boundary line. Its
rational smooth point is affine, and neither \(x=0\) nor
\(y=1\) divides the equation. The same substitution checks
\(x=1\) and \(y=0\).

Take the smooth projective normalization of that component.
Its genus is bounded absolutely by the plane degree bound,
and the Hasse–Weil estimate gives \(p+1+O(\sqrt p)\) points.
[Weil, Part II, §IV, Corollary 3 to Theorem 13](https://denisevellachemla.eu/Weil-courbes-varietes.pdf)
Bézout bounds the number on the boundary, at intersections
with other components, or over singular points, by an
absolute constant. The last assertion uses that, for primes
larger than the degree, a partial derivative of the reduced
component is nonzero and has lower degree. The number of
branches at any one point is also bounded by the degree.
For all large enough \(p\) there is consequently an affine
point with \(x,y\) different from zero and one.
Then \(z=u/(xy)\) is nonzero and the equation forces \(z\ne1\).
Discard the at most three solutions with \(x=y=z\), since
these have \(x^3=u\).

In the original unscaled inputs, when \(U,V\) are units the
Jacobian minor in slots \(i,j\) is
\[
        \frac{UV C(a_i-a_j)}{a_i a_j(C-a_i)(C-a_j)} .
\]
Thus away from the target the rank is two if and only if
two inputs differ. Keep the third input fixed and lift the
two chosen coordinates by Hensel's lemma to any prescribed
pair of lifts \(U,V\bmod p^e\).

If \(U=0\) and \(V\ne0\bmod p\), set \(a_1=0\) and
choose \(a_2,a_3\) neither zero nor \(C\), with
\((C-a_2)(C-a_3)=V/C\). The first factor can be any
nonzero element except \(C\) and the one choice making the
second factor equal to \(C\), so a choice exists for
\(p\ge5\). The minor in slots \(1,2\) is a unit, since
\(\partial U/\partial a_1\) and \(\partial V/\partial a_2\)
are units while \(\partial U/\partial a_2=0\).
The case \(V=0\), \(U\ne0\) follows by replacing every
\(a_i\) by \(C-a_i\). If both coordinates vanish, take
\(a_1=0,\ a_2=C,\ a_3\notin\{0,C\}\). The same
minor is diagonal with unit entries. These witnesses
therefore also lift at all depths.

## Primes on the target

For \(p\mid C\), let
\[
               W=\frac{V-(-1)^k U}{C},
\]
an integral polynomial in the \(a_i\) with \(C\) fixed.
On unit inputs it satisfies
\[
             W\equiv (-1)^{k-1}U\sum_i a_i^{-1}\pmod p .
\]
For \(k=3\), the data \(U,W\bmod p\) thus prescribe the
nonzero product \(P\) and the sum \(S\) of the three
nonzero inverses \(y_i=a_i^{-1}\). Eliminate the third
inverse using the sum. The projective plane cubic for the
first two is
\[
                       XY(SZ-X-Y)-P Z^3=0.
\]
At \([1:0:0]\) its \(Y\)-derivative equals \(-1\).
As above there is a geometric component defined over
\(\mathbb F_p\) through that rational smooth point.
None of \(X=0,Y=0,Z=0\) is a component, by substitution
and \(P\ne0\). The smooth normalization and Hasse–Weil
bound supply an affine point with \(X,Y\ne0\), outside the
finitely many points lying over singularities or intersections.
The third inverse \(S-X-Y\) is nonzero by the product equation.
The all-equal solutions again form a set of at most three
points and can be discarded.

For product and sum, the minor in two inverse slots is
\(P(1/y_i-1/y_j)\), a unit for unequal entries.
The changes \(a_i\mapsto a_i^{-1}\) and
\((U,W)\mapsto(U^{-1},(-1)^{k-1}W/U)\) have unit
Jacobian determinants on the unit region. The rank for
\((U,W)\) is therefore two whenever the chosen inputs
differ. For \(e>v_p(C)\), prescribe \(U\bmod p^e\)
and choose \(W\bmod p^e\) extending its prescribed value
modulo \(p^{e-v_p(C)}\). Hensel lifting gives the desired
\(U,V\bmod p^e\). If \(e\le v_p(C)\), every unit \(U\)
occurs by taking all but one input equal to \(1\).

For \(k>3\), fix the additional inputs with both \(a_i\)
and \(C-a_i\) units at any of these sufficiently large
primes. The product coordinates are independently
rescaled by units, and their scaling factors have relative
sign \((-1)^{k-3}\) modulo the common target divisor.
The \(k=3\) conclusions then give the full natural range.

## Small primes and integer lifts

For primes on \(C\), the two-input local argument of (CT4)
gives a fixed positive fraction at each depth. Fixing further
inputs equal to \(1\) independently rescales the two coordinates
by units. For an odd prime \(p\ge5\) away from \(C\), one can
take two distinct roots modulo \(p\) neither zero nor \(C\).
The \((p-2)(p-3)/2\) unordered choices give distinct pairs,
with a unit minor by the calculation above. Every lift of each
pair occurs. Further inputs can be fixed with both roles units.

For \(p=3\) away from \(C\), choose two integer roots
\(a,b\equiv-C\bmod3\) with \(b-a=3\). Prescribe a pair
\((U,S)\equiv(ab,a+b)\bmod27\), with \(S\) the sum of the
two moving roots. The discriminant \(S^2-4U\) is \(9\)
times a unit congruent to \(1\bmod3\), so it is a square at
all required depths. This realizes at least \(3^{-6}\) of
the pairs since \(V=U+C^2-CS\) is an invertible affine
change from \((U,S)\) at this prime. Below depth three
project the chosen classes. For \(p=2\) away from \(C\),
the witness \(a=0,b=C\bmod2\) has a unit minor, so every
even-even pair lifts for two inputs. Fixing further inputs
equal to \(2\) multiplies \(V\) by an odd unit and reduces
the size of the \(U\)-image by at most \(2^{k-2}\).
These are fixed positive fractions for the finitely many
small primes.
[Local target-prime counts for two factors (CT4)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#local-range-of-positive-complete-products)

The Chinese remainder theorem combines the local input
choices and the local output images. All sufficiently large
prime factors contribute their full natural range, so the
product of the losses is bounded below by a constant
depending only on \(k\). This proves (CR1).

For integer lifts, fix an allowed input residue \(a_i\bmod Q\).
It is already coprime to all primes of \(C\) that divide \(Q\).
Inclusion–exclusion over the other primes of \(C\) gives
\[
  \#\{A_i\in I_i:A_i\equiv a_i\bmod Q,\ (A_i,C)=1\}
  =\frac{|I_i|}{Q}\prod_{\substack{p\mid C\\p\nmid Q}}\left(1-\frac1p\right)
        +O(2^{\omega(C)}).
\]
Here the main term follows by solving the additional coprimality
tests on the progression, since \(Q\) is invertible modulo
each squarefree combination of these other primes. Under the
interval hypothesis the main term is \(C^{\eta-o(1)}\),
while the error is \(C^{o(1)}\). Choose \(A_i\) sequentially
and avoid the at most \(k-1\) previously selected values.
Then \(B_i=C-A_i>0\) and \((A_i,B_i)=1\). This proves
(CR2) and allows all summands to be tall when intervals
are chosen in a fixed interior range of \([1,C-1]\).

## Relevance to the small-budget image bound

These results impose the additive equation, primitivity and
prescribed real intervals. A power saving for the abc image
must also account for the radical caps on both free letters.
The complete-product range above can have a constant fraction
of its ambient size at large composite moduli even with all
three numerical heights comparable.

The combinatorial large sieve of Croot, Mao, Pohoata, Sheffer
and Yip assumes pointwise control of nonzero difference
multiplicities and small prime-reduction images. Its saving
for a fixed reduction density is subpower in the ambient
height. A height-first bound for the occupied abc image
would need a different source of loss that depends on the
free-letter arithmetic. [Combinatorial large sieve, Theorem 2.2](https://arxiv.org/html/2606.17487v1)

Yang's simultaneous-progression result concerns the
distribution of smooth numbers averaged over two modulus
ranges, with a fixed residue in one of them. The
relative logarithmic error and the marginal smooth-number
hypotheses do not provide a power saving for two simultaneous
shifted free-letter radical caps. [Yang, Theorem 1.1](https://arxiv.org/html/2611.01522v1)

Weighted-box determinant estimates over a fixed number field
count points on a prescribed algebraic surface or curve.
They can complete equations once the other power factors
are fixed. The needed price for varying the diffuse factors
and their occupied residues is outside that completion.
[Viglino, Theorems 1.1–1.2](https://arxiv.org/html/2610.35030v1)
