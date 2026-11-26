# Fractional-linear recovery of projective reciprocal pairs

The split genus-two recovery extends from affine equations in \(F\)
to rational fractional-linear equations. A potentially large
exception must preserve the conjugate branch pair. Such an
exception can lift only at an exact cube determinant. Shared
prime support gives a separate bound even on this exceptional
family of transformations.

## Recovery and the branch-pair exceptions

Let \(F(z)=z-z^{-1}+(1-z)^{-1}\). Use the two Möbius coordinates
\[
 U(z)=\frac{z-\eta}{z-\bar\eta},\quad
 T(t)=\frac{t-\alpha}{t-\bar\alpha},\quad
 \eta=\frac{1+\sqrt{-3}}2,\quad
 \alpha=\frac{1+3\sqrt{-3}}2 .
\]
Over \(E=\mathbb Q(\sqrt{-3})\) one has \(T(F(z))=U(z)^3\).
Let \(M\in\operatorname{PGL}_2(\mathbb Q)\) with a primitive
integer matrix of height at most \(X^K\), for fixed \(K\).
If \(M\) does not preserve \(\{\alpha,\bar\alpha\}\), then
\[
            F(x)=M(F(y)),\qquad x,y\in(0,1/2]\cap\mathbb Q
 \tag{FM1}
\]
has \(O_{K,\varepsilon}(X^\varepsilon)\) ordered solutions.
There is no height cap on \(x,y\) in this bound.

Indeed, set \(t=F(x)\) and \(v=U(x)U(y)\). Then
\[
                    v^3=T(t)\,T(M^{-1}(t)).
\]
Its four branch locations are \(\alpha,\bar\alpha,M(\alpha),
M(\bar\alpha)\). They are nonreal and finite. Any overlap
of the two conjugate pairs would give equality of the
pairs. The four branch orders are \(1,-1,1,-1\),
so the smooth cubic cover has genus two. Its branch
positions and the constant at infinity have polynomial
height, as do their denominators and the norms of the
nonzero branch differences after clearing denominators.
Use the split \(1-\zeta_3\) descent in (GR1), including
the primes of the constant at infinity in the finite
excluded set. Uniform Mordell–Lang bounds the \(E\)-points
by \(\exp(O(1+\omega(D)))\) with a polynomial-height
integer \(D\). Each point has at most nine lifts
to ordered arguments of \(F\), giving (FM1).
[Split genus-two descent and recovery](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#affine-pair-recovery-over-a-fixed-quadratic-field)

If \(M\) preserves the branch pair, then for some
\(\gamma\in E^\times\) of norm one and \(e=1\) or \(-1\),
\[
                       T(M(t))=\gamma T(t)^e.
\]
A rational solution requires \(\gamma\) to be a cube in \(E\).
Conversely, if \(\gamma=\delta^3\), all possible solutions
are on the three rational fractional-linear graphs
\[
                     U(x)=\zeta_3^j\delta\,U(y)^e,\qquad j=0,1,2.
\]
Each graph is defined over \(\mathbb Q\). In fact \(\delta\)
has norm one, and complex conjugation in the \(U\)-coordinate
acts by inversion on elements coming from rational arguments.
For an identity on the rational projective line this description
also allows \(U^{-1}(1)=\infty\). Restrict to the given interval
when counting solutions.

## Exact cube determinant

A branch-preserving primitive integer matrix with a solution
must have determinant of absolute value an exact cube. After
composing with \(t\mapsto 1-t\) if it swaps the branch points,
it fixes them individually. The polynomial of \(\alpha\) is
\(t^2-t+7\), so the matrix now has the form
\[
       \begin{pmatrix}c+d&-7c\\c&d\end{pmatrix},
       \qquad (c,d)=1,\qquad
       \det M=d^2+cd+7c^2=N_{E/\mathbb Q}(d+c\alpha).
 \tag{FM2}
\]
Writing \(z=d+c\alpha\), the multiplier in the \(T\)-coordinate
is \(\bar z/z\). Away from \(3\), no prime ideal divides both
\(z\) and \(\bar z\), because their difference is
\(3c\sqrt{-3}\) and \((c,d)=1\). A cube multiplier therefore
forces every exponent of \(z\) away from the prime over \(3\)
to be divisible by three.

For the prime over \(3\), put \(\omega=\zeta_3\), so
\(\alpha=2+3\omega\) and \(z=d+2c+3c\omega\).
If \(3\mid N(z)\), then \(c\) is a \(3\)-adic unit and
\(z/3=(d+2c)/3+c\omega\) is integral. Its order at the
ramified prime is zero or one, since its reduction modulo
\(3\mathbb Z[\omega]\) is nonzero. Thus \(v_3(N(z))\) is
zero, two or three. If it is two, Euclidean factorization
in \(\mathbb Z[\omega]\) gives \(z/3=u\,w^3\) for a unit
\(u\) and a \(3\)-adic unit \(w\). The cube condition on
\(\bar z/z\) requires \(\bar u/u=1\), since a nontrivial
third root of unity is not a cube in \(E\). Hence \(u=\pm1\).
But \(w^3\) is rational modulo \(3\mathbb Z[\omega]\),
by cubing its two coefficients in characteristic three.
This contradicts the nonzero \(\omega\)-coefficient
\(c\bmod3\) of \(z/3\). The only possible exponents in
the determinant are therefore zero or three at \(3\).
Composition with \(t\mapsto 1-t\) preserves the absolute
determinant and primitivity of the matrix.

There are \(O(H^{2/3})\) branch-preserving primitive
matrices of height at most \(H\) for which the cube
condition holds. Indeed \(z=u(1-\omega)^j w^3\) with
\(j=0\) or \(3\) and one of six units \(u\), so
\(|w|=O(H^{1/3})\) in the Eisenstein lattice. The swap
case contributes at most another factor of two.
In particular determinants of absolute value one yield
only the identity and \(t\mapsto1-t\) among the possible
exceptions, by (FM2). Neither admits distinct arguments
in the sorted interval of (FM1). Thus (FM1) has a uniform
subpower bound for distinct sorted arguments whenever the
absolute primitive determinant is not a cube, and also for
absolute determinant one.

## Möbius records with shared radical support

The new pair bound can be used when a family of arithmetic
pairs comes with a dictionary of fractional-linear maps.
Let \(n=X^L\) and use sorted primitive positive rows
\(a+b=c\), \(X/2<c\le X\), with
\(\operatorname{rad}(abc)\le n^B\). For fixed
\(0<\rho\le B\), suppose we count distinct row pairs
whose two values \(F(a/c)\) and \(F(a'/c')\) satisfy an
equation (FM1) for at least one map in a dictionary
of size at most \(n^\eta\) and of uniformly polynomial
height. Require the rows to have common radical at least
\(n^\rho\):
\[
                   \operatorname{rad}\gcd(abc,a'b'c')
                                      \ge n^\rho .
\]
Then for arbitrary fixed positive \(\xi,\varepsilon\) and
sufficiently small fixed budget their number is
\[
                    \ll X^{\{\eta+(B-\rho)/3+\xi\}L+\varepsilon}.
 \tag{FM3}
\]
The count is per pair, even if there are several maps.
On a single prescribed complete-letter fibre with total
moving radical cap \(n^B\), the same proof gives
\(\eta+(B-\rho)/4\), with the shared radical taken
outside the prescribed prime support of the fixed
letter and any initial pool on the two free letters.

For a non-branch-preserving map use (FM1). For a
branch-preserving map with bottom-left matrix entry
\(c_M=0\), the two possibilities are identity or
reflection, ruled out for distinct sorted arguments.
For any remaining map, all common primes lie in the
fixed support of \(c_M\det M\). In fact the reduced
denominator of \(F(a/c)\) is \(abc\). At a common prime
outside the determinant both reduced projective values
specialize to infinity, so the reduced map fixes
infinity and forces \(p\mid c_M\). The initial row
therefore has outside radical cap \(n^{B-\rho}\)
after adjoining this fixed pool. Count it by the
effective radical bound globally and the outside
complete-letter bound in a fibre. The second row is
then prescribed by injectivity of \(F\) and primitivity.
Summing over maps proves (FM3). In particular a
dictionary of size exponent \(\eta<\rho/3\) globally
or \(\eta<\rho/4\) on the fibre gives a strict saving
from the one-row bound. The dictionary assumption
is an extra input. The proof does not construct
such a dictionary from unlabelled support overlaps.
[Effective radical counting with fixed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Complete-letter bound](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Source comparison

The simplest cubic field-isomorphism literature gives
Tschirnhausen transformations for cyclic cubic
specializations, including formulas for coincidence
of their splitting fields. Here the splitting fields
of the input specializations are trivial, and the
varying datum is a prescribed fractional-linear
relation between their target values. The split
genus-two descent gives the uniform count on the
non-stabilizer part. The exact cube filter on integral
matrix determinants and the shared support estimate
address the additional arithmetic constraints in
the pair problem.
[Hoshi and Miyake, cubic generic polynomials and the field intersection problem](https://arxiv.org/abs/0710.0287)
· [Hoshi, simplest cubic fields and cubic Thue equations](https://arxiv.org/abs/0810.3374)
## Near-quadratic reciprocal expansion on the arithmetic image

Let \(A=\{F(x):x\in\mathcal A\}\), where the rational arguments
\(\mathcal A\subset(0,1/2]\) have height at most \(X\). For every
nonzero rational \(\lambda\) of height at most \(X^K\), fixed \(K\),
\[
                         |A+\lambda A^{-1}|
                          \gg_{K,\varepsilon}|A|^2 X^{-\varepsilon}.
 \tag{FM4}
\]
Also, if \(R\) is a rational Möbius map of polynomial coefficient
height, defined and nonzero on \(A\), then
\[
                             |A R(A)|
                          \gg_\varepsilon |A|^2 X^{-\varepsilon},
 \tag{FM5}
\]
with the constant allowed to depend on the coefficient-height
exponent. There is no radical constraint in either bound.

Indeed \(F\) has no rational zero by the rational-root test on its
cubic numerator. Fix a possible output \(t\) of the sum, necessarily
of polynomial height. The relation is \(F(x)=M_t(F(y))\)
with \(M_t(v)=t-\lambda/v\). Formula (FM2), or its composition
with reflection, gives the only branch-preserving possibilities
\((t,\lambda)=(1,7),(0,-7)\). In both cases the primitive
determinant has absolute value \(7\), ruling out solutions by the
cube condition. Now (FM1) bounds every output multiplicity by
\(O_{K,\varepsilon}(X^\varepsilon)\), proving (FM4).

For (FM5) a nonzero product output \(t\) uses \(M_t(v)=t/R(v)\).
There can be at most one branch-preserving \(M_t\), because
\(M_tM_{t'}^{-1}\) is multiplication by \(t/t'\), and the identity
is the only rational scaling preserving \(\{\alpha,\bar\alpha\}\).
That output has at most \(|A|\) ordered representations by
injectivity of \(R\). All other outputs have subpower multiplicity
by (FM1). This gives (FM5), with bounded \(|A|\) covered by
nonemptiness when \(A\ne\varnothing\).
[Nonvanishing of the scalar projective coordinate](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#fixed-product-recovery-for-the-same-projective-coordinate)

Yao Zhi's Theorem 1.1 gives the lower bound \(|A|^{10/9}\) in
characteristic zero for arbitrary finite sets in the reciprocal
sumset problem, and for product sets with a Möbius map whose four
matrix entries are nonzero and with no poles or zeros on \(A\).
It also gives positive-characteristic estimates. The above
near-quadratic bounds specialize to the image of the fixed
arithmetic function \(F\) at rational arguments of bounded
height. They use a coefficient-uniform bound on exact output
multiplicity from cyclic-cover descent.
[Yao Zhi, Theorem 1.1](https://arxiv.org/html/2611.26661v1#S1)
