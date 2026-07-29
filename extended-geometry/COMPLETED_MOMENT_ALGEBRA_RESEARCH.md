# Completed moment coordinates and the all-power Casimir ladder

## 1. Status and answer

Let
\[
 V_d=\operatorname{End}(\operatorname{Sym}^d),\qquad
 R_d=\mathbb Q[V_d]^{\operatorname{SL}_2},\qquad
 {\cal A}_d=\mathbb Q[\mu_1,\mu_2,\ldots].
\]
This note tests whether a small set of quadratic Casimirs and one
apolar-odd invariant repairs the pure moment coordinates in degrees
\(d=3,4,5,6\), and studies the known moment-zero witness families in
arbitrary degree.

The quadratic-completion and power-witness statements below are
all-degree theorems.  Finiteness and global reconstruction remain open
and are supported only by the explicitly labelled bounded calculations.

1. Completing the degree-two piece of \(R_d\) requires exactly \(d-1\)
   added quadratic invariants.  Thus the exact numbers are two, three,
   four, and five for \(d=3,4,5,6\).
2. Once the moment field is known to equal the apolar-fixed field, one
   nonzero odd invariant always generates the remaining quadratic
   extension.  This part needs no further computation.
3. The moment prefixes have exact good-prime Jacobian ranks
   \(13,22,33,46\) in degrees \(3,4,5,6\), certifying full
   transcendence degree there.
4. Low-weight searches find no rational formula for the missing
   quadratics or for the square of the first odd invariant, even after
   adjoining the natural quadratic repairs.  These are bounded
   nonexistence statements for a linear-denominator ansatz.
5. The single invariant \(q_2\) removes the primitive all-moment-zero
   witness \(F_4\) and all its radial multiples.  It does not remove
   \(F_4^3\), so a \(q_2\)-only proposal cannot be uniform in \(d\).
6. On the complete diagonal slices in all three degrees, the full moment
   field is exactly the apolar-reversal fixed field and has generic
   degree two.  This is an exact finite-fiber theorem, not merely a
   Jacobian or finite-field test.
7. The same fixed-field statement holds on certified one-phase quotient
   slices in every nonzero phase for \(d=3,5\).  On these slices the first
   \(d+3\) moments are already finite at the origin.  Cross-direction
   slices in quintic phases one and two have nonzero \(c_{234}\), so they
   also certify the remaining odd orientation step on the invariant
   quotient.
8. For every \(m\geq1\), the first nonzero quadratic Casimir on
   \(F_4^m\) is \(q_{2\lceil m/2\rceil}\).  Hence in ambient degree \(d\)
   the prefix through \(q_{2\lfloor(d+4)/8\rfloor}\) separates every
   known power-radial witness.

Consequently
\[
 {\cal B}_d={\cal A}_d[q_2,c_d]
 \tag{1.1}
\]
is a viable low-degree candidate but not an all-degree one.  A
witness-complete sparse candidate is
\[
 \boxed{
 {\cal B}^{\mathrm{pow}}_d=
 {\cal A}_d[
 q_2,q_4,\ldots,q_{2\lfloor(d+4)/8\rfloor},c_d].}
 \tag{1.2}
\]
Here \(c_3\) is a degree-four apolar-odd invariant and
\(c_d=c_{234}\) has degree three for \(d\geq4\).  The larger algebra
\[
 {\cal A}_d[q_2,q_4,\ldots,q_{2d-2},c_d]
 \tag{1.3}
\]
is the canonical *quadratically completed* candidate.  Neither algebra
is yet proved finite or birational in general; witness separation is not
a nullcone or reconstruction theorem.

## 2. The exact quadratic completion

Write the multiplicity-free Casimir decomposition as
\[
 A=A_0+A_2+\cdots+A_{2d}
 \tag{2.1}
\]
and put
\[
 q_{2r}=\operatorname{tr}(A_{2r}^2),\qquad 0\leq r\leq d.
 \tag{2.2}
\]
Schur's lemma gives
\[
 (R_d)_2=\langle q_0,q_2,\ldots,q_{2d}\rangle,
 \qquad \dim(R_d)_2=d+1.
 \tag{2.3}
\]
The two moment-generated directions are
\[
 q_0=\frac{\mu_1^2}{d+1}
 \tag{2.4}
\]
and
\[
 \boxed{
 \mu_2=\sum_{r=0}^d
 \binom{2d+1}{d-r}q_{2r}.}
 \tag{2.5}
\]
For the three degrees in question, the coefficient rows in (2.5) are
\[
\begin{array}{c|l}
d&\text{coefficients of }(q_0,q_2,\ldots,q_{2d})\\ \hline
3&(35,21,7,1)\\
4&(126,84,36,9,1)\\
5&(462,330,165,55,11,1).
\end{array}
\tag{2.6}
\]
Hence
\[
 \dim (R_d)_2/({\cal A}_d)_2=d-1.
\tag{2.7}
\]
Adding \(q_2,q_4,\ldots,q_{2d-2}\) is one minimal linear completion;
\(q_{2d}\) is then recovered from (2.5).

This distinguishes two meanings of “smallest.”  To reproduce the full
degree-two vector space, the answer is rigidly \(d-1\).  To recover only
the fraction field or to cut out the nullcone, fewer quadratics may
suffice because higher moments can impose nonlinear relations.  In
particular, the proposed \(d=4\) pair \(q_2,q_6\) leaves one polynomial
quadratic direction missing, although it may still be enough
rationally.

## 3. Orientation is exactly one quadratic step

Let \(K_d=\operatorname{Frac}R_d\).  The apolar adjoint \(\tau\) is a
nontrivial involution of \(K_d\), so
\[
 [K_d:K_d^\tau]=2.
\tag{3.1}
\]
If \(0\ne c\in K_d\) satisfies \(\tau(c)=-c\), then
\[
 c^2\in K_d^\tau,\qquad c\notin K_d^\tau,
\]
and therefore
\[
 \boxed{K_d=K_d^\tau(c),\qquad
 \min_{K_d^\tau}(c)=T^2-c^2.}
\tag{3.2}
\]
Thus the sharp target
\[
 \operatorname{Frac}{\cal A}_4=K_4^\tau
\tag{3.3}
\]
automatically implies that one \(c_{234}\) generates the remaining
quadratic extension.  The hard part is only (3.3), not the second step.

Convenient first odd invariants are
\[
\begin{array}{c|c|c}
d&\deg c_d&c_d\\ \hline
3&4&\operatorname{tr}(A_2A_4A_6^2)\\
4&3&\operatorname{tr}(A_4A_6A_8)\\
5&3&\operatorname{tr}(A_4A_6A_8).
\end{array}
\tag{3.4}
\]
The \(d=3\) odd space in degree four is three-dimensional; the displayed
trace is merely one convenient nonzero choice.

## 4. Generic reconstruction tests

The quotient dimensions are
\[
 \dim R_3=13,\qquad\dim R_4=22,\qquad\dim R_5=33.
\tag{4.1}
\]
At deterministic integral points, reduction modulo \(1000003\) gives
\[
\begin{array}{c|c|c}
d&\text{moment prefix}&\text{Jacobian rank}\\ \hline
3&\mu_1,\ldots,\mu_{13}&13\\
4&\mu_1,\ldots,\mu_{22}&22\\
5&\mu_1,\ldots,\mu_{33}&33.
\end{array}
\tag{4.2}
\]
A nonzero minor modulo one prime is an exact characteristic-zero
algebraic-independence certificate.  It proves generic finiteness of the
moment field extension, but not its discrete degree.

The bounded relation search uses unrestricted weighted monomials in the
chosen base invariants.  For a quadratic \(q\), it tests
\[
 Q_D(\mu,\mathbf q)+qP_{D-2}(\mu,\mathbf q)=0,
\tag{4.3}
\]
and for an odd invariant of degree \(e\), it tests
\[
 Q_D(\mu,\mathbf q)+c_d^2P_{D-2e}(\mu,\mathbf q)=0.
\tag{4.4}
\]
For \(d=3,4,5\), every tested intersection is zero through total weight
ten.  The quadratic targets were tested over
\({\cal A}_d\) and \({\cal A}_d[q_2]\).  The odd-square target was tested
over all four bases
\[
 {\cal A}_d,\quad {\cal A}_d[q_2],\quad
 {\cal A}_d[q_2,q_6],\quad
 {\cal A}_d[q_2,\ldots,q_{2d-2}]
\tag{4.5}
\]
Thus none of (4.3)--(4.4) exists over \(\mathbb Q\) with the corresponding
tested support.  The separate
degree-four search extends the moments-only \(c_{234}^2\) test through
weight sixteen.

This is evidence about the *complexity* of a rational reconstruction,
not evidence against fixed-field equality.  A formula may require a
higher-weight denominator, and if the moment field is smaller than the
fixed field then \(c_d^2\) can have a higher minimal polynomial.

On every diagonal slice for \(d=3,4,5\), the exact result is stronger:
the full restricted moment field equals the reversal-fixed field and has
degree two.  The first \(d+1\) moments give a parameter quotient of length
\((d+1)!\), and the first \(d+2\) moments through the selected integral
point cut out exactly that point and its reversal:
\[
\begin{array}{c|c|c}
d&\text{parameter length}&\text{next-moment fiber length}\\ \hline
3&24&2\\
4&120&2\\
5&720&2.
\end{array}
\tag{4.6}
\]
For \(d=5\), a homogeneous quotient of length \(720\) modulo \(32003\)
proves that the characteristic-zero projective zero fiber is empty:
otherwise its proper image over \(\operatorname{Spec}\mathbb Z\) would
meet every special fiber.  The six characteristic-zero forms are
therefore a regular sequence and have quotient length \(1\cdots6=720\).
After an invertible linear change sending the point and reversal to
\(s=\pm1\), the exact rational first-seven-moment ideal is
\[
 (y_0,y_1,y_2,y_3,y_4,s^2-1).
\tag{4.7}
\]
The same certificate works in degrees three and four.  The degree-four
development is also described in
[`DEGREE_FOUR_MOMENT_FIELD.md`](DEGREE_FOUR_MOMENT_FIELD.md).

There is now an exact off-diagonal extension for \(d=3,5\).  Fix a
nonzero phase \(h\), matching apolar eigendirections \(B,C\) of phases
\(+h,-h\), and consider
\[
 C=\operatorname{diag}(a_0,\ldots,a_d)+bB+cC.
\tag{4.8}
\]
After the residual diagonal torus quotient the coordinates are
\[
 (a_0,\ldots,a_d,z=bc),
\tag{4.9}
\]
and apolar adjunction reverses the \(a_i\) and fixes \(z\).  The \(+1\)
eigenspaces are used except in the extreme odd phase \(h=d\), where the
unique directions are necessarily apolar odd; both \(b,c\) then change
sign, so \(z\) is still fixed.

For one exact direction pair in every phase, the first \(d+2\) moments
have a nonzero Jacobian determinant modulo \(32003\).  Their special
moment-origin fiber is nevertheless one-dimensional.  Adding
\(\mu_{d+3}\) makes it zero-dimensional:
\[
\begin{array}{c|c|c|c|c}
d&\text{phases}&
\dim Z_{\mathbb F_{32003}}(\mu_1,\ldots,\mu_{d+2})&
\operatorname{length}Z_{\mathbb F_{32003}}
 (\mu_1,\ldots,\mu_{d+3})&
\text{fixed-target fiber length}\\ \hline
3&1,2,3&1&54&2\\
5&1,2,3,4,5&1&1934&2.
\end{array}
\tag{4.10}
\]
Because these equations are weighted homogeneous, the finite special
origin fiber makes the characteristic-zero weighted-projective zero
fiber empty by properness.  Thus the restricted coordinate ring is
integral over the first \(d+3\) moments.

For the target through \((2,3,5,7,221)\) in degree three and
\((2,3,5,7,11,13,221)\) in degree five, an affine midpoint-direction
change gives the special-fiber ideal
\[
 (y_0,\ldots,y_{d-1},w,s^2-1).
\tag{4.11}
\]
Weighted homogenization has no points at infinity by the origin
certificate.  Properness and Nakayama bound the characteristic-zero
fiber length by two, while the two explicit rational reversal points
give equality.  Hence on every displayed raw slice the full moment field
is exactly the apolar-fixed field.

The quotient-level interpretation needs one further check.  The chosen
cubic slices and the aligned quintic slices have zero first apolar-odd
invariant, so their two raw points can lie in one
\(\operatorname{SL}_2\)-orbit.  In quintic phases one and two, cross
direction pairs instead give
\[
 c_{234}=-\frac{273686400}{7}
\tag{4.12}
\]
at the chosen lift and the opposite value at its reversal.  Those are
genuinely distinct invariant-quotient points.  On these two
off-diagonal slice families, moments recover the fixed field and one
odd invariant generates the remaining quadratic extension.  All these
claims remain slice theorems and do not determine
\(\operatorname{Frac}({\cal A}_d)\subset\operatorname{Frac}(R_d)\)
globally.

## 5. Finiteness and nullcone detection are one test

Let \(B\subset R_d\) be a finitely generated homogeneous subalgebra.  The
following are equivalent:

1. \(R_d\) is finite over \(B\);
2. \(\sqrt{B_+R_d}=(R_d)_+\);
3. the common zero of the chosen invariants on
   \(\operatorname{Spec}R_d\) is only the vertex;
4. their common zero on \(V_d\) is exactly the
   \(\operatorname{SL}_2\)-nullcone, which is the pair-linear one-sided
   nullcone.

Thus “finiteness at the moment origin” and “nullcone detection” do not
require separate global eliminations.  One finite augmented prefix with
the correct zero fiber proves both.

Hilbert-series multiplication supplies a necessary test for possible
parameter degrees.  It does not prove the zero-fiber equality.  The
bounded calculations give the following useful candidates.

\[
\begin{array}{c|l|l}
d&\text{unaugmented consecutive prefix}&
 \text{one-}q_2\text{ candidate}\\ \hline
3&
1,\ldots,13\text{ fails at degree }63&
\mu_1,\ldots,\mu_{12},q_2\\
4&
1,\ldots,22\text{ passes the Hilbert test}&
\mu_1,\ldots,\mu_{21},q_2\\
5&
1,\ldots,33\text{ first fails at degree }483&
\mu_1,\mu_2,\mu_3,\ldots,\mu_{30},\mu_{32},\mu_{33},q_2.
\end{array}
\tag{5.1}
\]
Each displayed one-\(q_2\) system has full modular Jacobian rank in the
reproducible calculation.  “Passes” means only that the computed
Hilbert numerator is nonnegative and has no tail through the recorded
cutoff.  It is not a system-of-parameters theorem.

The \(d=5\) gaps in (5.1) are forced by the Hilbert test: among corrections
of the natural one-\(q_2\) degree sequence with total degree increase at
most fifteen, the least compatible correction has increase two and
replaces the final moment pair \(31,32\) by \(32,33\).

## 6. The known moment-zero witnesses

Let \(F_4\) be the exact bidegree-\((4,4)\) all-moment-zero witness and
\(F_5=RF_4\) its propagated degree-five form.  Exact Casimir projection
gives
\[
\begin{array}{c|l|c}
d&(q_0,q_2,\ldots,q_{2d})(F_d)&c_d(F_d)\\ \hline
4&(0,-864,2016,0,0)&0\\
5&(0,-24192,48384,0,0,0)&0.
\end{array}
\tag{6.1}
\]
Therefore \(q_2\) removes both known semistable points from the augmented
moment origin.  The odd invariant cannot do so: at these witnesses the
apolar adjoint is an \(\operatorname{SL}_2\)-translate, so every odd
invariant vanishes.

Equation (6.1) explains why an even repair is logically necessary for
finiteness in degrees four and five.  It does not prove that \(q_2\)
removes every semistable component.  For \(d=3\), whether the all-moment
zero fiber already equals the nullcone remains open, so there is not yet
an analogous all-order point to test.

## 7. Degree-by-degree frontier

### Degree three

The moments have full transcendence degree, but the first thirteen do not
form a parameter system.  Adding one independent quadratic produces a
Hilbert-compatible, full-rank candidate
\[
 (\mu_1,\ldots,\mu_{12},q_2).
\tag{7.1}
\]
The full zero fiber of (7.1), and hence finiteness, is open.  Generic
birationality after adjoining one degree-four odd invariant is also open.
On the four-dimensional diagonal slice, however, moments already recover
the reversal-fixed field and the odd orientation completes it
birationally.  The five-dimensional one-phase quotient slices now give
the same raw fixed-field theorem in every phase, but the tested first odd
invariant vanishes there, so they do not yet give a quotient-moving
orientation test.

### Degree four

The all-moment algebra is not integral in \(R_4\) because of \(F_4\).
The single \(q_2\) removes that witness, and
\[
 (\mu_1,\ldots,\mu_{21},q_2)
\tag{7.2}
\]
is a full-rank Hilbert-compatible candidate.  The suggested
\((q_2,q_6,c_{234})\) augmentation is therefore plausible but not
minimal by any proved field criterion; polynomial degree-two completion
would require one further quadratic.  The full fixed-field equality is
open, while its diagonal-slice analogue is proved with parameter length
\(120\) and an exact two-point fiber.

### Degree five

The new rank-\(33\) certificate proves full moment transcendence degree.
The consecutive first thirty-three degrees first fail the necessary
Hilbert numerator test very late, at degree \(483\), which explains why a cutoff at
degree \(100\) misses the obstruction.  The corrected one-\(q_2\)
candidate in (5.1) has full rank and passes the same necessary test.
The propagated witness is removed by \(q_2\).  Its zero fiber and generic
degree remain open.  On the six-dimensional diagonal slice, the first six
moments are finite parameters and the first seven have exactly the two
reversal-related points in the certified fiber, proving fixed-field
equality on that slice.  The seven-dimensional one-phase quotient slices
also have finite moment algebra and degree-two fixed field in every
phase.  Cross slices in phases one and two have nonzero \(c_{234}\), so
they certify the odd quadratic orientation extension on nonempty
off-diagonal quotient families.

## 8. Next decisive computations

The most informative next steps are:

1. **Finiteness first.**  Compute the zero fibers of (7.1), (7.2), and
   the \(d=5\) system in (5.1), branchwise under the lowest nonzero
   Casimir component.  A unit certificate off the nullcone would prove
   both finiteness and nullcone detection.  The first quartic
   nonzero-\(F_2\) normal-jet calculation is now recorded in
   [`DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md`](DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md):
   four linear pivots leave quadratic and cubic jet dimensions six and
   four, so quartic/quintic synchronization and the \(F_2=0\) boundary
   remain necessary.
2. **One complete generic fiber.**  On a rational quotient chart, find a
   finite fiber of the augmented map and prove that it consists of one
   orbit.  For moments without the odd coordinate, prove that it consists
   of exactly the two \(\tau\)-paired orbits.
3. **Rational even reconstruction.**  Express a transcendence basis of
   \(R_d^\tau\), rather than every low-degree generator, in the moment
   field.  For \(d=4\), the twenty-two even parameters in
   [`DEGREE_FOUR_MOMENT_FIELD.md`](DEGREE_FOUR_MOMENT_FIELD.md) are the
   natural target.
4. **Do not enlarge the quadratic set prematurely.**  The full
   degree-two completion is canonical, but current evidence makes
   \({\cal A}_d[q_2]\) the sharper finiteness experiment.  Add \(q_6\) or
   the remaining Casimirs only if an explicit semistable residual survives.

The bounded calculations are replayed by
[`research_completed_moment_algebra.py`](../scripts/research_completed_moment_algebra.py).
With the repository virtual environment, prime \(1000003\), maximum
weight ten, and three extra evaluation points, its
[`completed_moment_algebra_bounded_tests.json`](../artifacts/generated-results/completed_moment_algebra_bounded_tests.json)
output has SHA-256
`590fe262178bc4e8f11f3b633be9649ae82050afab962768e6e7946a1b15aa7c`.
It records every tested support and the Hilbert and Jacobian certificates.
The exact three-degree diagonal theorem is replayed by
[`verify_completed_moment_diagonal_fields.py`](../scripts/verify_completed_moment_diagonal_fields.py).
Its
[`completed_moment_diagonal_fields.json`](../artifacts/generated-results/completed_moment_diagonal_fields.json)
artifact has SHA-256
`17fdcbd88f261c7aff209e86207d1a3d9170fd5099db592e66cb004da063ff10`.
The exact cubic and quintic single-phase extension is replayed by
[`verify_completed_moment_single_phase_fields.py`](../scripts/verify_completed_moment_single_phase_fields.py).
Its
[`completed_moment_single_phase_fields.json`](../artifacts/generated-results/completed_moment_single_phase_fields.json)
artifact has SHA-256
`3626114820cf156c16327827ec9d8659c216964bd74401c2ea1ef8cc647b8f1e`.
The stronger degree-four-only tests are in
[`research_degree_four_moment_field.py`](../scripts/research_degree_four_moment_field.py),
[`verify_degree_four_tau_even_parameters.py`](../scripts/verify_degree_four_tau_even_parameters.py),
and
[`verify_degree_four_diagonal_moment_field.py`](../scripts/verify_degree_four_diagonal_moment_field.py).

## 9. Automatic missing-invariant scan through degree six

The low-degree comparison can be made uniform without first guessing a
trace word.  Refine the weight-zero-minus-weight-two character calculation
by the apolar parity
\[
 \tau|_{\operatorname{Sym}^{2r}}=(-1)^r.
\tag{9.1}
\]
The connected contraction-preserving group is
\(\operatorname{PGL}_2\).  After taking its invariant ring, (9.1) is the
canonical residual \(C_2\)-character used here.  No classification of
larger nonlinear automorphisms of the moment structure is needed or
claimed.
This gives the exact dimensions of \(R_{d,n}^+\) and \(R_{d,n}^-\).
Because the displayed moment prefixes have full Jacobian rank, the moment
monomials through polynomial degree six are independent and have dimension
equal to the partition number \(p(n)\).  Subtracting their span from the
even character gives the missing even directions; the entire odd character
is missing because every moment is even.

For \(d=3,4,5,6\), the first result is:
\[
\begin{array}{c|c|c|c|c|c}
d&\dim R_{d,2}&\dim({\cal A}_d)_2&
\dim R_{d,2}/({\cal A}_d)_2&
\text{first odd degree}&\text{odd dimension there}\\ \hline
3&4&2&2&4&3\\
4&5&2&3&3&1\\
5&6&2&4&3&2\\
6&7&2&5&3&5.
\end{array}
\tag{9.2}
\]
Thus the first missing invariant is always quadratic, always
apolar-even, and has multiplicity \(d-1\).  An explicit basis of the
quadratic invariant space is
\[
 q_0,q_2,\ldots,q_{2d},
\tag{9.3}
\]
while the moment subspace is spanned by \(\mu _1^2,\mu _2\).  The
canonical quotient basis may be represented by
\[
 \boxed{q_2,q_4,\ldots,q_{2d-2}.}
\tag{9.4}
\]
This proves that \(q_2\) is the first member of a uniform
moment-completion family, not an isolated quartic phenomenon.

The odd cubic calculation also becomes uniform.  For three distinct
components, a cubic contraction \(c_{abc}\) exists precisely when
\(a+b\geq c\), and it is odd precisely when \(a+b+c\) is odd.  The scan
finds
\[
\begin{array}{c|l}
d&\text{odd cubic component triples}\\ \hline
4&(2,3,4)\\
5&(2,3,4),(2,4,5)\\
6&(2,3,4),(2,4,5),(2,5,6),(3,4,6),(4,5,6).
\end{array}
\tag{9.5}
\]
These lists have sizes \(1,2,5\), agreeing independently with the refined
Hilbert character.  Degree three has no odd cubic; its first odd
three-dimensional space occurs in degree four.

The exact modular moment-Jacobian ranks are
\[
\begin{array}{c|c|c}
d&\dim R_d&\operatorname{rank}d(\mu_1,\ldots,\mu_{\dim R_d})\\ \hline
3&13&13\\
4&22&22\\
5&33&33\\
6&46&46.
\end{array}
\tag{9.6}
\]
The new \(d=6\) minor modulo \(1000003\) is therefore a
characteristic-zero algebraic-independence certificate, not a numerical
rank estimate.

For the propagated all-order witnesses \(F_d=R^{d-4}F_4\), exact Casimir
projection gives
\[
\begin{array}{c|l|c}
d&(q_0,q_2,\ldots,q_{2d})(F_d)&c_{234}(F_d)\\ \hline
4&(0,-864,2016,0,0)&0\\
5&(0,-24192,48384,0,0,0)&0\\
6&(0,-967680,1741824,0,0,0,0)&0.
\end{array}
\tag{9.7}
\]
Hence one invariant, \(q_2\), removes every explicit all-order
semistable witness presently recorded in these degrees.  This is only a
witness-coverage statement.  The semistable moment-zero components have
not been classified, so (9.7) does not prove that \(q_2\) separates every
component.  In degree three there is no recorded all-order semistable
moment-zero point to evaluate.

The resulting candidate hierarchy is:
\[
\begin{array}{c|l|l}
d&\text{one-}q_2\text{ parameter candidate}&
\text{orientation-completed algebra}\\ \hline
3&(\mu_1,\ldots,\mu_{12},q_2)&{\cal A}_3[q_2,c_3]\\
4&(\mu_1,\ldots,\mu_{21},q_2)&{\cal A}_4[q_2,c_{234}]\\
5&(\mu_1,\ldots,\mu_{30},\mu_{32},\mu_{33},q_2)&
{\cal A}_5[q_2,c_{234}]\\
6&(\mu_1,\ldots,\mu_{45},q_2)&{\cal A}_6[q_2,c_{234}].
\end{array}
\tag{9.8}
\]
Every system in the middle column has full Jacobian rank modulo
\(1000003\) and passes the recorded bounded Hilbert-numerator necessary
test.  Neither condition proves that its zero fiber is the nullcone.  If
polynomial equality in degree two is required rather than sparse
finiteness, replace the last column by the canonical completion
\[
 {\cal A}_d[q_2,q_4,\ldots,q_{2d-2},c_d].
\tag{9.9}
\]
If the moment field is eventually proved equal to the \(\tau\)-fixed
field, any one nonzero odd invariant then supplies the remaining quadratic
orientation extension; the multiplicities in (9.2) and (9.5) do not
require adjoining every odd direction.

The scan is replayed by

```bash
.venv/bin/python scripts/research_completed_moment_algebra.py \
  --degrees 3 4 5 6 --invariant-cutoff 6 \
  --skip-relation-tests --power-witness-cutoff 12 \
  --output artifacts/generated-results/automatic_missing_invariants_d3_d6.json
```
The generated artifact has SHA-256
`a526b52feaeeba5cfca4d4903232e931443b0e0e6e92709c54e84e486573fb7c`.
It stores the even and odd invariant dimensions in every polynomial
degree through six, the moment-monomial quotient dimensions, the Hilbert
tests, the parameter Jacobians, and the exact witness evaluations.
It does not serialize a full coefficient basis for every
\(R_{d,n}\) through degree six.  It constructs the exact character
dimensions and explicit bases for the first missing quadratic spaces and
the first odd cubic spaces.  A complete basis-producing version should
next solve the sparse raising-operator kernel in each multidegree, then
reduce the resulting transvectants against moment monomials.

This division of labor matches the standard computational invariant-theory
pipeline: multigraded Poincaré series for direct sums of binary forms
control the target dimensions
([Bedratyuk](https://arxiv.org/abs/1108.1555)), while candidate parameter
degrees must satisfy independent invariant-ring restrictions
([Brouwer--Draisma--Popoviciu](https://doi.org/10.1007/s00031-015-9335-8)).
Neither source supplies the contraction-specific moment quotient or the
witness evaluations in (9.7).

## 10. Radial recurrence and the Casimir ladder

The \(d=4,5,6\) witness values are instances of an exact all-degree
recurrence.  Let
\[
 T_d:V_d\longrightarrow V_{d+1},\qquad T_d(f)=Rf.
\tag{10.1}
\]
It is \(\operatorname{SL}_2\)-equivariant and preserves the component
label \(r\).  In operator coordinates \(A=C^TD_d\), direct coefficient
conversion gives
\[
 (T_dA)_{ba}=(d+1-a)A_{ba}+aA_{b-1,a-1},
\tag{10.2}
\]
where an out-of-range entry is zero.

Let \(E_d,F_d\) be the raising and lowering matrices on
\(\operatorname{Sym}^d\).  The highest and lowest vectors in the
\(\operatorname{Sym}^{2r}\) component can be taken as
\(E_d^r,F_d^r\).  Equation (10.2) gives
\[
 T_d(E_d^r)=(d-r+1)E_{d+1}^r,\qquad
 T_d(F_d^r)=(d-r+1)F_{d+1}^r.
\tag{10.3}
\]
The trace normalization is
\[
 \operatorname{tr}(E_d^rF_d^r)
 =(r!)^2\binom{d+r+1}{2r+1}.
\tag{10.4}
\]
Schur's lemma and the ratio of (10.4) in consecutive degrees therefore
prove
\[
 \boxed{
 q_{2r}^{(d+1)}(Rf)
 =(d-r+1)(d+r+2)q_{2r}^{(d)}(f).}
\tag{10.5}
\]

Applying (10.5) to the quartic seed yields, for every \(d\geq4\),
\[
\boxed{
\begin{aligned}
q_2(F_d)&=-\frac{(d-1)!(d+2)!}{5},\\
q_4(F_d)&=\frac{(d-2)!(d+3)!}{5},\\
q_{2r}(F_d)&=0\qquad(r=0,3,4,\ldots,d),
\end{aligned}}
\qquad F_d=R^{d-4}F_4.
\tag{10.6}
\]
In particular
\[
 q_4(F_d)=-\frac{d+3}{d-1}q_2(F_d).
\tag{10.7}
\]
The binomial coefficient ratio between the \(q_2\) and \(q_4\) terms of
\(\mu_2\) is the same \((d+3)/(d-1)\), explaining the exact moment
cancellation.  Thus \(q_2\) separates the entire infinite radial chain,
not only the three degrees in (9.7).

The power witnesses show why the rest of the Casimir ladder should not be
discarded.  Exact sparse Casimir projection gives
\[
\begin{array}{c|c}
m&\text{first nonzero tested quadratic on }F_4^m\in V_{4m}\\ \hline
1,2&q_2\\
3,4&q_4\\
5,6&q_6\\
7,8&q_8\\
9,10&q_{10}\\
11,12&q_{12}.
\end{array}
\tag{10.8}
\]
All earlier quadratics in each row vanish exactly, and the displayed one
is nonzero over \(\mathbb Q\).  Radial propagation by (10.5) preserves
this pattern for \(R^kF_4^m\).  In fact the pattern is all-degree:
\[
 \boxed{\text{first nonzero quadratic on }F_4^m
 =q_{2\lceil m/2\rceil}.}
\tag{10.9}
\]
Here is a direct proof.  On the chart \(R=1\), put
\[
 x=Z,\qquad u=T^2,\qquad W=\frac{1-u}{2x}.
\]
The quartic seed reduces to
\[
 \boxed{2F_4=x^{-1}(1+x)\bigl(1-u(1+x)^2\bigr).}
\tag{10.10}
\]
Consequently
\[
 2^mF_4^m=\sum_{a=-m}^{2m}x^aH_{m,a}(u),
\qquad
 H_{m,a}(u)=
 \sum_{j=0}^m(-1)^j\binom mj
 \binom{m+2j}{m+a}u^j.
\tag{10.11}
\]

The phase \(-a\) line in the \(\operatorname{Sym}^{2r}\) component,
for \(a\geq0\) and \(r-a\) even, is represented on this chart by
\[
 x^a C_{r-a}^{a+1/2}(t),\qquad u=t^2,
\tag{10.12}
\]
up to a nonzero scalar.  The opposite phase is represented by
\[
 x^{-a}(1-t^2)^a C_{r-a}^{a+1/2}(t).
\tag{10.13}
\]
These are the standard associated-Legendre, or equivalently Gegenbauer,
weight vectors on the adjoint quadric.

To project \(H_{m,a}\) onto (10.12), it suffices by Gegenbauer
orthogonality to evaluate
\[
 M_{m,a,\ell}=
 \int_{-1}^1
 t^{2\ell}(1-t^2)^aH_{m,a}(t^2)\,dt
\tag{10.14}
\]
for \(0\leq\ell\leq(r-a)/2\).  Substitution of (10.11) gives
\[
 M_{m,a,\ell}
 =
 \sum_{j=0}^m(-1)^j\binom mj\binom{m+2j}{m+a}
 B\left(j+\ell+\frac12,a+1\right).
\tag{10.15}
\]
Let \(s=\lceil m/2\rceil\).  If \(\ell<s-a\), then
\[
\binom{m+2j}{m+a}
B\left(j+\ell+\frac12,a+1\right)
\tag{10.16}
\]
is a polynomial in \(j\) of degree \(m-1\).  Indeed, after writing the
binomial as
\[
 \frac{(2j-a+1)_{m+a}}{(m+a)!},
\]
all \(a+1\) denominator factors from the beta function cancel because
\[
 2\ell+2a+1\leq m.
\]
The alternating sum (10.15) is therefore the \(m\)-th finite difference
of a polynomial of degree \(m-1\), so
\[
 \boxed{M_{m,a,\ell}=0\qquad(\ell<s-a).}
\tag{10.17}
\]

At the boundary \(a=s,\ell=0\), all beta denominator factors except
\(2j+2s+1\) cancel.  Polynomial division writes (10.16) as a polynomial
of degree at most \(m-1\) plus a nonzero constant divided by
\(2j+2s+1\).  The polynomial part again disappears, while
\[
 \sum_{j=0}^m\frac{(-1)^j\binom mj}{2j+2s+1}
 =\frac12B\left(s+\frac12,m+1\right)\ne0.
\tag{10.18}
\]
Thus phase \(-s\) is nonzero in component \(s\).  Phase \(+s\) is also
nonzero: its coefficient is proportional to
\[
 \int_{-1}^1H_{m,-s}(t^2)\,dt.
\tag{10.19}
\]
The same division argument leaves a nonzero remainder at the single
uncancelled factor \(2j+1\).

For \(r<s\), equation (10.17) kills every nonpositive phase of
\(P_r(F_4^m)\), so that component is strictly one-sided and
\(q_{2r}(F_4^m)=0\).  In component \(s\), it kills every nonpositive
phase except \(-s\); equations (10.18)--(10.19) show that the pairing of
the two extreme phases is nonzero.  Hence
\[
\boxed{
q_{2r}(F_4^m)=0\ (r<\lceil m/2\rceil),\qquad
q_{2\lceil m/2\rceil}(F_4^m)\ne0,}
\tag{10.20}
\]
proving (10.9) for every \(m\geq1\).

The computation through \(m=12\) records the stronger full phase-support
pattern
\[
\operatorname{phases}P_r(F_4^m)=
\begin{cases}
\{1,3,\ldots,r\},&r<s,\ r\text{ odd},\\
\{2,4,\ldots,r\},&r<s,\ r\text{ even},\\
\{-s\}\cup\{1,3,\ldots,s\},&r=s,\ s\text{ odd},\\
\{-s\}\cup\{2,4,\ldots,s\},&r=s,\ s\text{ even}.
\end{cases}
\tag{10.21}
\]
Only the one-sided vanishing and the two extreme coefficients needed for
(10.20) are proved in all degrees.  Equality in the stronger support
formula (10.21) remains bounded evidence.

There is also a closed count for the odd cubic space.  Adding a new
largest component \(c=d\) contributes
\(\binom{\lfloor d/2\rfloor}{2}\) odd triangle triples.  Indeed, write
\(u=d-b\) and \(a=u+v\).  The triangle and parity conditions become
\(v=2w+1\) and
\[
 u\geq1,\qquad w\geq0,\qquad
 u+w\leq\lfloor d/2\rfloor-1,
\]
which has exactly that many solutions.  Summation gives
\[
\boxed{
\dim (R_{d,3})^-=
\begin{cases}
d(d-1)(d-2)/24,&d\text{ even},\\
(d+1)(d-1)(d-3)/24,&d\text{ odd}.
\end{cases}}
\tag{10.22}
\]
This recovers \(1,2,5\) for \(d=4,5,6\) and shows cubic orientation
directions grow cubically.  One nonzero odd invariant can still generate
the quadratic field extension over the full fixed field, but it cannot
span the odd polynomial space.

Together, (10.5)--(10.22) sharpen the augmentation strategy:

1. use \(q_2\) first, since it is the minimum-degree separator and kills
   the complete radial chain;
2. on a surviving branch, test \(q_4,q_6,\ldots\) in order rather than
   adding all of them immediately;
3. add \(c_{234}\) only for apolar orientation—it cannot replace the
   even Casimir tests on the power witnesses.

For a fixed degree \(d\), the known power-radial witnesses are
\[
 R^{d-4m}F_4^m,\qquad1\leq m\leq\lfloor d/4\rfloor.
\tag{10.23}
\]
Consequently, the all-power theorem (10.20) shows in every degree that
the sparse prefix
\[
\boxed{
q_2,q_4,\ldots,q_{2s_d},\qquad
s_d=\left\lfloor\frac{d+4}{8}\right\rfloor}
\tag{10.24}
\]
separates every witness in (10.23).  This is much smaller than the full
\(d-1\)-element
quadratic completion, but unlike \(q_2\) alone it responds to the complete
known power family.  Every apolar-odd invariant still vanishes on these
witnesses because apolar adjunction is a
\(\operatorname{PGL}_2\)-translate on \(F_4\), hence on all its powers and
radial multiples.

The updated automatic scan verifies (10.4)--(10.6), regresses the
beta-sum vanishing and boundary nonvanishing in (10.17)--(10.19) through
\(m=32\), verifies the cubic count, and records the exact stronger
phase-support table (10.21) through \(m=12\).  The finite beta regression
is a check on the identities used in the proof; the finite-difference
argument above is what proves (10.20) for arbitrary \(m\).
