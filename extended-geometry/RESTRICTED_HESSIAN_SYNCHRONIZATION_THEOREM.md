# Restricted Hessian synchronization: the generic theorem and exact frontier

Work over an algebraically closed field \(k\) of characteristic zero.  For
a proper decomposition word \(\mathbf d\) of \(N\), let \(H_{\mathbf d}\)
be its canonical Hessian residual ideal and let
\(\lambda_{\mathbf d}\) be the reconstructed coefficient of \(x\).  For a
set \(D\) of requested words put

\[
 H_D=\sum_{\mathbf d\in D}H_{\mathbf d}.
\]

There are two different synchronization statements:

\[
\begin{array}{rcl}
\text{generic synchronization}
&\Longleftrightarrow&
\lambda_{\mathbf d}-\lambda_{\mathbf e}\in\sqrt{H_D},\\[2mm]
\text{scheme-theoretic synchronization}
&\Longleftrightarrow&
\lambda_{\mathbf d}-\lambda_{\mathbf e}\in H_D.
\end{array}                                                   \tag{0.1}
\]

The square-free relation-graph hypothesis proves the first statement
completely.  It does not, without an additional infinitesimal theorem,
prove the second.  This note gives the noncircular proof, then records the
remaining exact frontier.

## 1. The missing-line lemma

The key observation is independent of Ritt normal forms.

> **Lemma 1.1 (Abhyankar--Moh missing-line obstruction).**  Let \(K\) be
> a characteristic-zero field and let
> \[
>  P=A\circ B,\qquad Q=C\circ E
> \]
> be monic polynomials of the same degree.  Suppose
> \[
>  P-Q=\delta x+\epsilon,\qquad \delta\ne0.
> \]
> Then either \(\deg B\mid\deg E\) or \(\deg E\mid\deg B\).

**Proof.**  The displayed identity gives

\[
 x=\delta^{-1}\bigl(A(B)-C(E)-\epsilon\bigr)\in K[B,E].
\]

Thus \(K[B,E]=K[x]\).  The algebraic Abhyankar--Moh epimorphism
theorem says that two nonconstant polynomials generating \(K[x]\) have
comparable degrees.  Hence \(\deg B\mid\deg E\) or
\(\deg E\mid\deg B\).  \(\square\)

The point is that this argument starts with two potentially different
linear lifts.  It therefore avoids the circular use of an ordinary
collision theorem, which starts with one already synchronized polynomial.

## 2. Generic square-free synchronization

Assume that the requested incidences are the normalized decomposition words
encoded by the relation graph: every displayed factor cut is part of the
incidence data, rather than a refinement inferred from equality of two
coarse lifts.  (Inferring such a refinement before synchronization would be
circular.)  The strongly connected components form a chain.  A singleton
component has one fixed position.  Inside a nontrivial square-free
component, any two requested orders are joined by adjacent swaps of distinct
vertices.

Consider one such swap.  Cancel the common prefix and retain the common
suffix \(S\).  The two right-component degrees at the swap are

\[
 uS\quad\text{and}\quad vS,
\]

where \(u\ne v\) are vertices in the square-free block.  Gcd refinement
makes \(u\) and \(v\) coprime after common decorations have been removed.
Consequently neither \(uS\) nor \(vS\) divides the other.  Lemma 1.1 says
that two decomposable lifts differing by a nonzero linear term cannot
realize this swap.  Their linear coefficients are therefore equal.

Following a swap path and then the chain of singleton components proves:

> **Theorem 2.1 (restricted generic synchronization).**  Suppose every
> requested incidence is a normalized decomposition word in the relation
> graph, and every strongly connected component is either an acyclic
> singleton or one tame power/Dickson block whose degree multiset is
> square-free.  Then, for every requested pair,
> \[
>  \lambda_{\mathbf d}-\lambda_{\mathbf e}\in\sqrt{H_D}.       \tag{2.1}
> \]
> Equivalently, the Hessian intersection is generically synchronized.

**Proof.**  Let \(\mathfrak p\) be any prime containing \(H_D\) and pass to
the fraction field of \(R_N/\mathfrak p\), then to its algebraic closure.
The canonical reconstructions give, for each requested word, a monic
decomposable polynomial

\[
 F+\lambda_{\mathbf d}x
\]

with the same coefficients in degrees at least two.  The adjacent-swap
argument above makes all of their linear coefficients equal.  Hence every
difference belongs to every prime over \(H_D\), which is exactly (2.1).
\(\square\)

This proof also supplies component completeness for the reduced question:
there cannot be an extra unsynchronized minimal component off the
power/Dickson charts.  Ziegler's multi-collision theorem is needed only
after synchronization, to classify the resulting field-valued collision.

The theorem is slightly stronger than the proposed “generically” clause:
the Abhyankar--Moh obstruction applies at every field-valued point, not
only at the generic points of the standard charts.

For an intersection specified only by coarse outer cuts, Theorem 2.1 applies
after the corresponding refined factor incidences have been established.
The ordinary gcd-refinement theorem cannot establish them here, because it
assumes that the two decompositions belong to one polynomial.  The direct
degree-thirty arguments below do not make this assumption.

## 3. Why radical membership is not ideal membership

Theorem 2.1 proves

\[
 \lambda_{\mathbf d}-\lambda_{\mathbf e}\in\sqrt{H_D}.
\]

It proves membership in \(H_D\) if \(H_D\) is radical.  Without
radicality, an additional infinitesimal argument is necessary.  A
power/Dickson parametrization kills every nilpotent on its reduced support,
so testing those parametrizations cannot by itself detect an embedded
linear defect.

The tempting formal argument has a real missing premise.  On a synchronized
field-valued Ritt chart, write a common right decoration as \(T(x)\).  The
free \(K[T]\)-module decomposition

\[
 K[x]=K[T]\oplus xK[T]\oplus\cdots\oplus x^{r-1}K[T]
\]
shows that the first-order missing-line equation lies in the \(x\)-summand.
For an adjacent tame collision it has the form

\[
 A'(B)u-C'(E)v=\delta.                                       \tag{3.1}
\]

In both the power and Dickson normal forms, the two derivative factors in
(3.1) have a common monic factor of positive degree, so \(\delta=0\).
This proves transverse first-order rigidity on every smooth normal-form
chart.

What is not supplied by the field-valued Ritt theorem is formal
completeness: every Artin-ring deformation of the synchronized
intersection would have to come from those charts.  At a chart overlap or
an embedded primary component, higher-order deformation terms need not be
controlled by the first-order calculation.  Therefore (3.1) is not a
proof that

\[
 \lambda_{\mathbf d}-\lambda_{\mathbf e}\in H_D              \tag{3.2}
\]
in complete generality.

An all-degree proof of (3.2) needs one of the following genuinely new
inputs:

1. a scheme-theoretic Ritt normal-form theorem over local Artin
   characteristic-zero algebras;
2. radicality of the restricted Hessian intersections;
3. a direct universal remainder identity showing that the linear defect
   belongs to the Hessian residual ideal.

## 4. A universal quadratic scheme lemma

One substantial part of (3.2) does have a direct ring proof.

> **Lemma 4.1 (quadratic edge over a ring).**  Let \(S\) be a
> \(\mathbb Q\)-algebra.  Suppose \(P=A\circ B\) and \(Q=C\circ E\) are
> monic of degree \(2m\), with outer degrees \(2\) and \(m\),
> respectively.  If \(P-Q\) has degree at most one, then its linear
> coefficient is zero.

**Proof.**  Complete the square in the outer quadratic and translate the
source so that the inner quadratic in the second decomposition is even.
For a monic degree-\(m\) polynomial \(G(y)\), the equation becomes

\[
 G(y)^2-\text{an even polynomial}=\delta y+\epsilon.
\]

Writing \(G=G_{\mathrm{ev}}+G_{\mathrm{odd}}\) and subtracting the equation
at \(-y\) from the equation at \(y\) gives

\[
 4G_{\mathrm{ev}}G_{\mathrm{odd}}=2\delta y.                 \tag{4.1}
\]

If \(m\) is even, \(G_{\mathrm{ev}}\) is monic of degree \(m\); if \(m\)
is odd, \(G_{\mathrm{odd}}\) is monic of degree \(m\).  Multiplication by
a monic polynomial is degree-faithful over an arbitrary ring.  Since
\(m\ge2\), (4.1) forces the other parity part to vanish and then
\(\delta=0\).  \(\square\)

Applied to degree thirty, Lemma 4.1 proves the exact membership for the
primitive pair

\[
 \lambda_2-\lambda_{15}\in H_2+H_{15}.                       \tag{4.2}
\]

Unlike a reduced normal-form calculation, this proof is valid over every
\(\mathbb Q\)-algebra and therefore controls all primary thickenings.

## 5. The cubic remainder certificate

The other primitive degree-thirty pair admits a compact exact remainder
calculation.  Center the inner cubic and the outer cubic:

\[
 D(y)=y^3+py+q,\qquad T(G)=G^3+uG.
\]

Write the monic degree-ten polynomial \(G\) uniquely in the free
\(\mathbb Q[D]\)-module basis \(1,y,y^2\):

\[
 G=g_0(D)+y\,g_1(D)+y^2g_2(D).
\]

The Hessian residual equations say that the entire \(y^2\)-remainder of
\(G^3+uG\) vanishes and that all positive-\(D\) coefficients of its
\(y\)-remainder vanish.  The possible synchronization defect is the
remaining constant coefficient of that \(y\)-remainder.

[`verify_cubic_remainder_synchronization.py`](../scripts/verify_cubic_remainder_synchronization.py)
constructs these equations in thirteen variables and performs the exact
characteristic-zero reduction.  Singular returns a basis of size \(184\)
and reduces the defect to zero.  Hence

\[
 \lambda_3-\lambda_{10}\in H_3+H_{10}.                        \tag{5.1}
\]

This is a universal ideal identity, not a sample-point or radical test.

## 6. Updated degree-thirty frontier

The original five open pairs split as follows.

| cuts | status | reason |
|---|---|---|
| \(\{2,15\}\) | exact | universal quadratic parity lemma |
| \(\{3,10\}\) | exact | cubic free-module remainder certificate |
| \(\{2,3\}\) | generically exact; scheme open | degree-six core transported through a quintic |
| \(\{2,5\}\) | generically exact; scheme open | degree-ten core transported through a cubic |
| \(\{3,5\}\) | generically exact; scheme open | degree-fifteen core transported through a quadratic |

Thus twelve of the fifteen degree-thirty pair memberships are now proved
scheme-theoretically.  The remaining three all have a nontrivial common
right decoration.  They would follow from the precise transport statement:

> **Required transport lemma.**  If two tame decomposable lifts differ by
> a polynomial of degree below the common right degree, then the common
> right factor furnished by Ritt theory on the reduced fiber lifts strongly
> enough to force the difference to be a polynomial in that factor.

Over fields this follows indirectly from Theorem 2.1 and the ordinary
greatest-common-right-component theorem.  Over nonreduced rings it is the
unproved scheme-theoretic step.

## 7. Reproducible status

The rigorous result is now:

> The square-free relation-graph hypothesis proves generic Hessian
> synchronization in all degrees.  The stronger ideal-membership theorem is
> proved when the relevant Hessian ideal is radical, for every primitive
> quadratic edge, and for the degree-thirty cubic remainder edge.  It is not
> yet proved for arbitrary primary thickenings transported through
> singleton blocks.

The main external inputs are the algebraic
[Abhyankar--Moh theorem](https://encyclopediaofmath.org/wiki/Abhyankar%E2%80%93Moh_theorem),
Ziegler's
[Tame Decompositions and Collisions](https://arxiv.org/abs/1402.5945),
Müller--Zieve's
[On Ritt's Polynomial Decomposition Theorems](https://arxiv.org/abs/0807.3578),
and Bodin's ring-valued approximate-root construction in
[Decomposition of polynomials and approximate roots](https://arxiv.org/abs/0910.1676).
