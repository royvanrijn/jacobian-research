# Cellular Postnikov transitivity for ideal flags

> **Status.** The transitivity and multi-flag statements in Sections 1--3
> are general algebraic theorems.  They apply to arbitrary connective
> derived rings, and their ideal formulas require only ordinary surjective
> ring maps.  The finite-module tower and cellular associated graded are
> implemented exactly over \(\mathbb Q\).  What remains conjectural is the
> Hessian--Ritt comparison identifying the actual completed intersection
> with this cellular totalization uniformly in all degrees.

This note extracts the degree-thirty and degree-forty-two calculations in
[the cellular cotangent prototype](HESSIAN_RITT_CELLULAR_COTANGENT_PROTOTYPE.md)
into a reusable transitivity theorem.  The central distinction is between

* higher cotangent homology entering the first Postnikov sequence;
* non-split extensions inside that sequence; and
* ordinary Tor introduced by a non-flat finite base change.

These are different phenomena and should not be inferred from the same
finite-jet dimension count.

## 1. Transitivity in a connective derived category

Let

\[
 A\longrightarrow C\longrightarrow B
\]

be composable maps of connective commutative derived rings.  Cotangent
transitivity is the cofiber sequence

\[
 X:=B\otimes_C^{\mathbf L}L_{C/A}
 \longrightarrow
 Y:=L_{B/A}
 \longrightarrow
 Z:=L_{B/C}
 \overset{\partial}{\longrightarrow}X[1].                    \tag{1.1}
\]

Assume \(H_0(X)=0\).  The relevant part of the homology sequence is

\[
 H_2(Y)\longrightarrow H_2(Z)\longrightarrow H_1(X)
 \longrightarrow H_1(Y)\longrightarrow H_1(Z)\longrightarrow0.
                                                                    \tag{1.2}
\]

Define the **first Postnikov overlap**

\[
 \mathcal O_{A,C,B}
 =
 \operatorname{im}\bigl(H_2(Z)\longrightarrow H_1(X)\bigr).
                                                                    \tag{1.3}
\]

Exactness immediately gives the canonical short exact sequence

\[
 0\longrightarrow H_1(X)/\mathcal O_{A,C,B}
 \longrightarrow H_1(Y)
 \longrightarrow H_1(Z)\longrightarrow0.                    \tag{1.4}
\]

Thus three assertions have distinct strengths:

1. \(\mathcal O_{A,C,B}=0\) says that higher cotangent homology has zero
   image in the first transitivity sequence.
2. Splitting (1.4) is an additional extension problem in the heart of the
   standard \(t\)-structure.
3. Splitting the full triangle (1.1) is stronger still and forces the
   connecting morphism \(\partial\) to vanish.

In particular, a non-split map on \(H_1\) proves
\(\partial\ne0\), while overlap vanishing alone does not make the triangle
split.

## 2. The ideal formula

Let \(R\) be an ordinary commutative ring with ideals

\[
 I\subset J\subset K
\]

and put

\[
 A=R/I,\qquad C=R/J,\qquad B=R/K.
\]

All three maps are surjective, so their relative cotangent complexes have
zero \(H_0\).  Their first homology modules are the conormal modules

\[
\begin{aligned}
 H_1(L_{B/A})&=K/(I+K^2)=:N_I,\\
 H_1(L_{B/C})&=K/(J+K^2)=:N_J.
\end{aligned}                                                \tag{2.1}
\]

The base-change spectral sequence for the left term of (1.1) has only one
contribution to total degree one, because \(H_0(L_{C/A})=0\).  Therefore

\[
\begin{aligned}
 H_1(B\otimes_C^{\mathbf L}L_{C/A})
 &=
 B\otimes_C J/(I+J^2)\\
 &=
 J/(I+KJ)=:S_J,                                             \tag{2.2}
\end{aligned}
\]

where \(J^2\subset KJ\).  The map \(S_J\to N_I\) is induced by
\(J\subset K\).  Its kernel is exactly

\[
 \boxed{
 \mathcal O(I,J,K)
 =
 \frac{J\cap(I+K^2)}{I+KJ}.}                                \tag{2.3}
\]

Indeed, a class represented by \(j\in J\) vanishes in \(N_I\) precisely
when \(j\in I+K^2\).  Combining (1.4) and (2.3) proves:

> **Flag--cotangent theorem.** For every ideal flag
> \(I\subset J\subset K\), there is a canonical exact sequence
> \[
> 0\to S_J/\mathcal O(I,J,K)\to N_I\to N_J\to0.
> \]
> It is short exact with sector term \(S_J\) if and only if
> \[
> J\cap(I+K^2)=I+KJ.                                        \tag{2.4}
> \]

No completion, smoothness, characteristic-zero, or Noetherian hypothesis
is needed for this theorem.

## 3. Arbitrary ideal flags

Now take a finite flag

\[
 I_0\subset I_1\subset\cdots\subset I_r=K.                  \tag{3.1}
\]

For \(1\le i\le r\), define

\[
\begin{aligned}
 N_i&=K/(I_i+K^2),\\
 S_i&=I_i/(I_{i-1}+KI_i),\\
 \mathcal O_i
 &=
 \frac{I_i\cap(I_{i-1}+K^2)}
      {I_{i-1}+KI_i}.
\end{aligned}                                                \tag{3.2}
\]

Applying the flag--cotangent theorem at every stage gives

\[
 0\longrightarrow S_i/\mathcal O_i
 \longrightarrow N_{i-1}
 \longrightarrow N_i\longrightarrow0.                      \tag{3.3}
\]

Since \(N_r=0\), these sequences filter
\(N_0=H_1(L_{(R/K)/(R/I_0)})\).  If every overlap vanishes, its successive
layers are exactly the sector modules \(S_i\):

\[
 0\to S_i\to N_{i-1}\to N_i\to0.                            \tag{3.4}
\]

Each stage has an extension class

\[
 e_i\in\operatorname{Ext}^1_{R/K}(N_i,S_i).                 \tag{3.5}
\]

For a flag of length greater than two, the adjacent classes are not the
whole story.  Compatibility of their iterated cones is measured by higher
Toda or Massey-type operations.  In a cellular Ritt diagram, commuting and
braid two-cells are the first locations where these compatibility
conditions must be imposed.

## 4. Completion, Nakayama, and finite base change

Suppose \(R\) is local, \(B=R/K\), and the overlap in (2.3) is a finite
\(B\)-module.  For any ideal \(\mathfrak a\) in the Jacobson radical,

\[
 \mathcal O/\mathfrak a\mathcal O=0
 \quad\Longrightarrow\quad
 \mathcal O=0                                                \tag{4.1}
\]

by Nakayama's lemma.  This makes one finite quotient a valid overlap
certificate, provided the quotient module is computed before truncation.

If (3.4) is exact, tensoring with \(Q=B/\mathfrak a\) gives

\[
 \operatorname{Tor}_1^B(Q,N_i)\longrightarrow
 S_i\otimes_BQ\longrightarrow
 N_{i-1}\otimes_BQ\longrightarrow
 N_i\otimes_BQ\longrightarrow0.                             \tag{4.2}
\]

Consequently a kernel created only after finite base change is an image of
ordinary Tor.  It is not evidence that \(\mathcal O_i\ne0\).  Conversely,
any \(B\)-linear section of \(N_{i-1}\to N_i\) remains a section after
tensoring.  A single non-split finite quotient therefore proves formal
non-splitting.

## 5. Cellular totalization

Let \(\Gamma\) be a finite relation complex and suppose every cell carries
a tower (3.3), functorially under incidence maps.  The associated graded
places each sector layer \(S_i\) on the appropriate relative path block.
For a three-edge half-braid relative to its endpoints, that block has
cellular cohomology

\[
 H^0=0,\qquad H^1=S_i,\qquad H^2=0.                          \tag{5.1}
\]

The executable model accepts an arbitrary finite tower of equivariant
finite modules.  It verifies:

* commutation of all coordinate actions;
* equivariance and surjectivity of every adjacent map;
* exact restricted actions on every kernel layer;
* existence or nonexistence of compatible sections; and
* the associated-graded cellular totalization of all layers.

If the actual cotangent diagram satisfies derived cellular descent, its
skeletal filtration has the usual totalization spectral sequence

\[
 E_2^{p,q}
 =
 H_{\mathrm{cell}}^p(\Gamma;\mathcal H_q)
 \Longrightarrow
 H_{q-p}(L_{\mathrm{total}}).                               \tag{5.2}
\]

Equation (5.2) is the correct all-cell theorem target.  The algebraic
transitivity statements above do not by themselves prove the comparison
between \(L_{\mathrm{total}}\) and the completed Hessian intersection.

## 6. Degree-thirty and degree-forty-two regressions

In degree thirty, \(I_\partial=K\).  The tower has dimensions

\[
 1\longrightarrow0,
\]

one sector layer, and no spectator.  With the two-dimensional Dickson base,
its filled-braid totalization has

\[
 (\dim H^0,\dim H^1,\dim H^2)=(2,1,0).                      \tag{6.1}
\]

For the degree-forty-two conormal modules modulo
\((\tau,z)^2\), the actual nine-coordinate action matrices give the tower

\[
 6\longrightarrow2\longrightarrow0.                         \tag{6.2}
\]

Its exact kernel dimensions and splitting profile are

\[
 (4,2),\qquad(\mathrm{non\mbox{-}split},\mathrm{split}).     \tag{6.3}
\]

The first four-dimensional kernel is the effective sector after the
non-flat base change.  The completed sector source has dimension six
modulo the same base square; HRCELL5 proves that the two-dimensional
difference is the Tor image in (4.2), while the completed overlap (2.3)
vanishes.

The implementation also verifies a synthetic three-layer tower

\[
 3\longrightarrow2\longrightarrow1\longrightarrow0
\]

to ensure that the code is not specialized to one- or two-layer flags.

## 7. Scope of the generalization

The following are now degree-independent theorems:

1. the derived first-Postnikov sequence (1.4);
2. the ideal overlap formula (2.3);
3. the arbitrary-flag filtration (3.3);
4. Nakayama overlap detection and the Tor distinction (4.1)--(4.2).

The remaining Hessian--Ritt work is geometric:

1. construct the ideal flags functorially for every labelled move and cell;
2. prove overlap vanishing uniformly for the universal move types;
3. prove compatibility of the extension classes on commuting and braid
   cells;
4. establish the derived cellular-descent comparison behind (5.2); and
5. compute or eliminate the higher algebraic \(H^2\) rows.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_cellular_postnikov_transitivity.py
```

The command writes
`artifacts/generated-results/cellular_postnikov_transitivity.json`.
It consumes the exact degree-forty-two conormal action matrices from
`artifacts/generated-results/degree42_ritt_conormal_transitivity.json`.
