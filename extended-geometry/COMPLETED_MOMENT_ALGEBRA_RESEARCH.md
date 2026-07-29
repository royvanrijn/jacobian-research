# Completed moment coordinates in degrees three through five

## 1. Status and answer

Let
\[
 V_d=\operatorname{End}(\operatorname{Sym}^d),\qquad
 R_d=\mathbb Q[V_d]^{\operatorname{SL}_2},\qquad
 {\cal A}_d=\mathbb Q[\mu_1,\mu_2,\ldots].
\]
This note tests whether a small set of quadratic Casimirs and one
apolar-odd invariant repairs the pure moment coordinates for \(d=3,4,5\).

The present answer is a sharp reduction plus bounded evidence, not a
global reconstruction theorem.

1. Completing the degree-two piece of \(R_d\) requires exactly \(d-1\)
   added quadratic invariants.  Thus the exact numbers are two, three,
   and four for \(d=3,4,5\).
2. Once the moment field is known to equal the apolar-fixed field, one
   nonzero odd invariant always generates the remaining quadratic
   extension.  This part needs no further computation.
3. The moment prefixes have exact good-prime Jacobian ranks
   \(13,22,33\) in degrees \(3,4,5\).  The \(d=5\) rank is a new bounded
   exact certificate that the moments have full transcendence degree
   there.
4. Low-weight searches find no rational formula for the missing
   quadratics or for the square of the first odd invariant, even after
   adjoining the natural quadratic repairs.  These are bounded
   nonexistence statements for a linear-denominator ansatz.
5. The single invariant \(q_2\) removes the known all-moment-zero
   semistable witnesses in both \(d=4\) and \(d=5\).  No computation here
   proves that it removes every semistable point.

Consequently the smallest serious uniform candidate is
\[
 \boxed{{\cal B}_d={\cal A}_d[q_2,c_d],}
 \tag{1.1}
\]
where \(c_3\) is a degree-four apolar-odd invariant and
\(c_d=c_{234}\) has degree three for \(d=4,5\).  The larger algebra
\[
 {\cal A}_d[q_2,q_4,\ldots,q_{2d-2},c_d]
 \tag{1.2}
\]
is the canonical *quadratically completed* candidate.  Neither algebra
is yet proved finite or birational for \(d=3,4,5\).

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

On the diagonal \(d=4\) slice, the exact result is stronger: the full
restricted moment field equals the reversal-fixed field and has degree
two.  See
[`DEGREE_FOUR_MOMENT_FIELD.md`](DEGREE_FOUR_MOMENT_FIELD.md).

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
open, while its diagonal-slice analogue is proved.

### Degree five

The new rank-\(33\) certificate proves full moment transcendence degree.
The consecutive first thirty-three degrees first fail the necessary
Hilbert numerator test very late, at degree \(483\), which explains why a cutoff at
degree \(100\) misses the obstruction.  The corrected one-\(q_2\)
candidate in (5.1) has full rank and passes the same necessary test.
The propagated witness is removed by \(q_2\).  Its zero fiber and generic
degree remain open.

## 8. Next decisive computations

The most informative next steps are:

1. **Finiteness first.**  Compute the zero fibers of (7.1), (7.2), and
   the \(d=5\) system in (5.1), branchwise under the lowest nonzero
   Casimir component.  A unit certificate off the nullcone would prove
   both finiteness and nullcone detection.
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
The stronger degree-four-only tests are in
[`research_degree_four_moment_field.py`](../scripts/research_degree_four_moment_field.py),
[`verify_degree_four_tau_even_parameters.py`](../scripts/verify_degree_four_tau_even_parameters.py),
and
[`verify_degree_four_diagonal_moment_field.py`](../scripts/verify_degree_four_diagonal_moment_field.py).
