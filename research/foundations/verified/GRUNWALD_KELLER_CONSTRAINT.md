# Grunwald problems with a Keller constraint

The Keller condition adds no new arithmetic obstruction when the polynomial
map may depend on the answer.  The genuine difficulty is exactly the
classical Grunwald--inverse-Galois problem in a chosen permutation
representation.  Once the required number field has been constructed, the
[local-to-global Keller theorem](LOCAL_GLOBAL_KELLER_FIBERS.md) compiles its
stem polynomial into a connected complete determinant-one fiber without
changing the field, its completions, or its normal closure.

This note makes that reduction precise, separates two commonly conflated
forms of the local data, records the cases that are already unconditional,
and gives an exact \(A_5\) example at \(2,3,5\).

## 1. The correct Grunwald datum

Let \(G\) act transitively and faithfully on a set \(\Omega\) of size
\(N\), and let \(H\) be a point stabilizer.  Thus
\(\Omega=G/H\) and the core of \(H\) is trivial.

A classical local datum at a place \(v\) is a continuous homomorphism

\[
 \phi_v:\Gamma_{\mathbb Q_v}\longrightarrow G
\]

up to conjugacy.  Its image is a decomposition subgroup \(D_v\leq G\).
Equivalently it is a \(G\)-Galois étale algebra, not necessarily a field.
The induced degree-\(N\) stem algebra is the finite étale
\(\mathbb Q_v\)-algebra attached to the \(D_v\)-set \(G/H\).

This distinction is essential for nonsolvable \(G\).  Every finite Galois
group over a \(p\)-adic field is solvable: wild inertia is a \(p\)-group,
tame inertia is cyclic, and the residue Galois group is cyclic.  Hence there
is no \(A_5\)-Galois **field** over \(\mathbb Q_p\).  An \(A_5\)-Grunwald
datum at \(p\) instead has solvable image \(D_p\leq A_5\), and the associated
quintic stem algebra may be a field or a product.

If the desired connected Keller fiber itself is required to be a
\(G\)-Galois field, use the regular action \(H=1\), of degree \(|G|\).
For the usual phrase “an \(A_5\)-field,” this note uses the standard quintic
meaning: a degree-five field whose Galois closure has group \(A_5\).

## 2. Exact reduction to the classical problem

> **Keller--Grunwald transfer principle.**  Let \(G\leq S_N\) be transitive
> and faithful, with \(N\geq3\).  Suppose a degree-\(N\) field \(K/\mathbb Q\)
> realizes prescribed stem algebras at finitely many places and the Galois
> closure of \(K\) has group exactly \(G\) in the given action.  Then there
> is an explicit determinant-one polynomial map
>
> \[
>  F:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q}
> \]
>
> of geometric degree \(N\) with a rational target \(y\) such that
>
> \[
>  F^{-1}(y)\simeq\operatorname{Spec}K
> \]
>
> scheme-theoretically.  The fiber is connected and complete, has the
> prescribed completions, and its Galois closure still has group \(G\).

Choose a monic primitive-element polynomial \(P\) for \(K\).  The
quadratic-gauge compiler in
[`jcsearch.keller_fiber`](../jcsearch/keller_fiber.py) produces \(F,y\) and
the displayed scheme isomorphism.  Because the quotient algebra is exactly
\(\mathbb Q[T]/(P)\), every completion and the normal closure are unchanged.
That proves the principle.

Conversely, a connected finite étale complete fiber is already the spectrum
of a number field.  Asking that its normal closure have group \(G\) and that
its completions match the local data is precisely the original arithmetic
problem.  Thus, with the map allowed to vary, the Keller-constrained and
classical existence problems are equivalent for the chosen representation.

The transfer also preserves specified inertia and solutions of embedding
problems, provided those structures have first been certified in the normal
closure.  The current compiler preserves the field; it does not by itself
construct or label its normal closure.

## 3. Cases that are already unconditional

### 3.1 Symmetric groups

For every \(N\geq3\), arbitrary compatible rank-\(N\) finite étale
\(\mathbb Q_p\)-algebras can be realized with global Galois closure exactly
\(S_N\).

Start with the coefficientwise construction in
[the local-to-global theorem](LOCAL_GLOBAL_KELLER_FIBERS.md).  At three new
large good primes impose the squarefree factor types

\[
 (N),\qquad (N-1,1),\qquad (2,1,\ldots,1).
\]

The first makes the polynomial irreducible.  The second gives an
\((N-1)\)-cycle; its fixed-point stabilizer is transitive on the other
\(N-1\) letters, so the global group is doubly transitive.  A doubly
transitive group containing a transposition contains every transposition,
and is \(S_N\).  For \(N=3\), transitivity and a transposition already give
\(S_3\).  Compilation therefore proves the full \(S_N\) Keller--Grunwald
theorem constructively.

### 3.2 Groups with a generic extension

Saltman's generic-extension theorem gives a positive Grunwald theorem for
every group admitting a generic Galois extension over the base field.  A
generic polynomial provides local parameter points; étale/Krasner stability
turns them into open local neighborhoods, weak approximation combines the
parameters, and Hilbert irreducibility with weak approximation retains the
full global group.

Consequently every such group has the Keller--Grunwald property in every
faithful transitive representation of degree at least three.  This includes
all \(S_N\) and \(A_5\).  Gene Ward Smith gives an explicit generic
\(A_5\) polynomial over \(\mathbb Q\), including specialization to every
characteristic-zero \(A_5\)-extension.

Primary references:

- D. J. Saltman,
  [*Generic Galois extensions and problems in field theory*](https://pmc.ncbi.nlm.nih.gov/articles/PMC348470/),
  *Advances in Mathematics* 43 (1982), especially Theorem 5.9.
- G. W. Smith,
  [*A Generic Polynomial for the Alternating Group \(A_5\)*](https://arxiv.org/abs/1210.4991),
  especially Theorem 2.
- J. König and D. Neftin,
  [*The local dimension of a finite group over a number field*](https://arxiv.org/abs/2007.05383),
  for the modern specialization and Hilbert--Grunwald formulation.

## 4. An exact \(A_5\) field prescribed at \(2,3,5\)

The following three stem algebras are prescribed:

\[
\begin{aligned}
 A_2&=\mathbb Q_2(\sqrt[3]{2})\times U_{2,2},\\
 A_3&=U_{3,5},\\
 A_5&=\mathbb Q_5[U]/
 \left(U^5+15U^4+50U^3-200U^2-1375U-1915\right),
\end{aligned}
\]

where \(U_{p,d}\) denotes the unramified degree-\(d\) extension of
\(\mathbb Q_p\).  The last displayed polynomial is Eisenstein at \(5\).

Put

\[
 H(T)=T^5-40T^3-110T^2-40T+32
\]

and \(K=\mathbb Q[T]/(H)\).  Then

\[
 \operatorname{Disc}(H)=580000^2.
\]

Modulo \(3\), \(H\) is irreducible.  Modulo \(23\), its factor degrees are
\((3,1,1)\).  Hence the transitive Galois group lies in \(A_5\), contains a
5-cycle and a 3-cycle, and is exactly \(A_5\).

At \(2\), the coefficient valuation points of \(H\) have lower Newton
polygon

\[
 (0,5)\longrightarrow(2,1)\longrightarrow(5,0).
\]

The slopes are \(-2\) and \(-1/3\).  Their residual polynomials are
\(Y^2+Y+1\) and \(Y+1\), respectively.  The first gives the unramified
quadratic and the second gives the unique totally tamely ramified cubic over
\(\mathbb Q_2\), namely \(\mathbb Q_2(\sqrt[3]{2})\).  Thus
\(K\otimes\mathbb Q_2\simeq A_2\).

At \(3\), irreducibility modulo \(3\) and the 3-adic unit discriminant give
\(K\otimes\mathbb Q_3\simeq U_{3,5}=A_3\).

At \(5\), translating \(T=U+3\) gives exactly

\[
 H(U+3)=U^5+15U^4+50U^3-200U^2-1375U-1915,
\]

so Eisenstein's criterion gives
\(K\otimes\mathbb Q_5\simeq A_5\).

### 4.1 Compilation into one full Keller fiber

The smaller inverse presentation

\[
 P(S)=S^5-5S^3+4S+\frac45
\]

defines the same field.  If \(\beta\) is the class of \(T\) modulo \(H\),
then

\[
 \alpha=
 \frac{-\beta^4+2\beta^3+36\beta^2+38\beta-36}{20}
\]

satisfies \(P(\alpha)=0\).  Both \(H\) and \(P\) are irreducible modulo
\(3\), so \(\mathbb Q(\alpha)=\mathbb Q(\beta)=K\).

Set

\[
 t=1+xy,\qquad
 q=t^2z-\frac45y^2(1+3t)
\]

and

\[
\begin{aligned}
 \Pi&=tq,\\
 B&=y-\frac{15}{4}xq+\frac54t^2x^3q^5,\\
 C&=x(5-3t)+\frac54x^3z-\frac34(xq)^5.
\end{aligned}
\]

The map

\[
 \boxed{\widetilde F=(\Pi,-B/2,C)}
\]

has determinant one, geometric degree five, and coordinate degrees
\((7,32,30)\).  At the rational target

\[
 \boxed{y=\left(1,0,-\frac25\right)}
\]

the quadratic-gauge reconstruction gives the scheme-theoretic identity

\[
 \boxed{
 \widetilde F^{-1}(y)
 \simeq\operatorname{Spec}\mathbb Q[S]/(P)
 \simeq\operatorname{Spec}K.
 }
\]

This is a connected complete Keller fiber, its global Galois closure is
exactly \(A_5\), and its completions at \(2,3,5\) are the three algebras
prescribed above.

The exact checker uses no Galois-group oracle and no numerical \(p\)-adic
factorization:

```bash
.venv/bin/python scripts/verify_a5_grunwald_keller_fiber.py
```

## 5. What remains open

The remaining frontier is arithmetic, not Keller-theoretic.

| case | present status | next verifiable theorem |
|---|---|---|
| \(S_N\) | unconditional and constructive | package the three auxiliary Frobenius witnesses in the local-global API |
| \(A_5\) | unconditional via a generic polynomial; explicit \(2,3,5\) example above | synthesize from user-supplied local \(A_5\)-torsors and emit a labelled decomposition certificate |
| groups with a generic polynomial | unconditional by Saltman, then compile | implement parameter-space weak approximation plus a group certificate |
| solvable \(G\) | odd-order and several broader classes are positive; even-order Wang/Brauer--Manin obstructions prevent a blanket statement | implement one obstruction-aware family, beginning with supersolvable groups |
| \(A_N,\ N\geq6\) | regular realizations give broad unramified specialization results; this note does not assert the full ramified Grunwald theorem | choose a specific regular family and prove a precise allowed-prime/inertia theorem |
| embedding problems with inertia | compilation is automatic after the global embedding problem is solved | certificate schema for \(I_v\triangleleft D_v\leq G\), ramification filtration, and the stem orbits |

Recent boundary results should be read with their prime restrictions:

- F. Motte,
  [*Hilbert irreducibility, the Malle conjecture and the Grunwald problem*](https://www.numdam.org/item/10.5802/aif.3567.pdf),
  proves quantitative unramified specialization results for regular groups.
- J. König,
  [*The Grunwald problem and specialization of families of regular Galois extensions*](https://arxiv.org/abs/1710.05548),
  treats families and prescribed decomposition/inertia behavior.
- E. Boughattas and D. Neftin,
  [*The Grunwald problem and homogeneous spaces with nonsolvable stabilisers*](https://aif.centre-mersenne.org/articles/10.5802/aif.3784/),
  includes \(A_5\)-kernel families away from primes dividing the group order.
- J. L. Demeio,
  [*Solvable Descent and the Grunwald Problem for Solvable Groups*](https://arxiv.org/abs/2604.18099),
  gives the 2026 solvable-group result up to the necessary Brauer--Manin
  obstruction.

For specified inertia, a local pair must first be admissible: in the tame
case \(I_v\) is cyclic and Frobenius acts by the residue-cardinality power;
in the wild case the wild inertia is a normal \(p\)-subgroup with the usual
ramification filtration.  A degree-\(N\) factorization alone records only
the orbit sizes of \(D_v\) on \(G/H\).  The next implementation therefore
needs normal-closure data, not a larger Keller compiler.
