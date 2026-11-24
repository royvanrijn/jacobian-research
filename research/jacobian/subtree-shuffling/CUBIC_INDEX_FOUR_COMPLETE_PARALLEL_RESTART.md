# Complete cubic parallel restarts have a sharp inverse-degree ceiling

Over a field of characteristic zero, add **every** admissible cubic helper to
the nine-coordinate pre-\(S\) syzygy block, in parallel. A further admissible
cubic may use these helpers nonlinearly and mix them with arbitrary quadratic
coefficients. Nevertheless, it cannot use the new \(S\) direction, and its
inverse has the same degree ceiling as the existing degree-207 construction.
This is a complete finite-dimensional classification of this restart layer.
It does not bound arbitrary cubic maps of Jacobian nilpotency index four or
an arbitrary number of successive restart layers.

## Statement and exact domain

Let \(K\) be any characteristic-zero field, and let \(H_0\) be a homogeneous
cubic map on a finite-dimensional coordinate space. Choose distinct coordinate
functions \(u,w\), and assume
\[
 J_0=JH_0,\qquad J_0^3=0,\qquad (H_0)_u=0,\qquad
 r=dwJ_0^2\text{ and }du\text{ are linearly independent over }K(x).
 \tag{1}
\]
Adjoin the nine independent coordinates
\(p_1,p_2,p_3,A,B,E,D,C,R\), with precisely the nine components before
\(P_S\) in [the degree-207 construction](CUBIC_INDEX_FOUR_SYZYGY58.md).
Write
\[
 A_0=K[u,w,p_1,p_2,p_3,A,B,E,D,C,R].
\]
All helpers and subsequent components below use only these displayed old
coordinates and their specified new helper coordinates; they cannot use
other coordinates of the base.

Let \(\delta,\epsilon\) be the two polynomial derivations of
[the pre-\(S\) row criterion](COMMUTING_SYZYGY_SLICES_AND_SCALAR_OBSTRUCTION.md),
fixing \(u,w\), and put
\[
 V=\{Y\in (A_0)_3:\delta Y=\epsilon Y=0\}.
 \tag{2}
\]
The subscript denotes ordinary homogeneous degree. Take an arbitrary finite
family \(Y_1,\ldots,Y_m\in V\), including repetitions or zero polynomials,
and adjoin independent coordinates \(y_i\) with components \(H_{y_i}=Y_i\).
These additions are parallel: no \(Y_i\) uses another \(y_j\).
Then adjoin any finite number of further coordinates, also in parallel,
with ordinary homogeneous cubic components
\[
 Q\in K[u,w,p_1,p_2,p_3,A,B,E,D,C,R,y_1,\ldots,y_m]_3.
 \tag{3}
\]
Only components \(Q\) preserving fourth Jacobian power zero are admitted.

For the next statement, use the auxiliary weights
\[
 (1,29,31,59,87,89,117,89,119,147,177)
 \tag{4}
\]
on the eleven coordinates of \(A_0\), in the displayed order. These weights
serve to split a finite linear system; they impose no hypothesis on the
weights of \(H_0\) or of the chosen helpers. The weight-207 subspace of
\(V\) is exactly \(K P_S\), where \(P_S\) is the displayed last cubic in
the degree-207 construction. Define \(\lambda(Y)\in K\) by
\[
 Y_{[207]}=\lambda(Y)P_S.
\]

**Theorem CompleteParallelRestart.** Under (1)--(3), every admissible \(Q\)
satisfies
\[
 \left(\sum_{i=1}^m\lambda(Y_i)\partial_{y_i}\right)Q=0.
 \tag{5}
\]
In particular, when a basis of \(V\) is used as the helper family, with its
last member \(P_S\), every admissible cubic is independent of that last
helper coordinate \(S\).

Assume additionally that \(F_0=\operatorname{id}-H_0\) is a polynomial
automorphism, with inverse \(G_0\), and put
\[
 D=\deg G_0,\qquad B=\deg (G_0)_w\ge1.
\]
Every map obtained by the two parallel layers above is a polynomial
automorphism and has inverse degree at most
\[
 \boxed{\max\{D,7B+4\}}.
 \tag{6}
\]
Its Jacobian has nilpotency index exactly four. For the retained
48-coordinate cube-zero base, \(D=B=29\); the sharp ceiling is therefore
207. It is attained by the existing component \(P_S\).

**Corollary UniversalParallelRestartCeiling.** For every base satisfying
(1), with no additional invertibility hypothesis, the entire specified
two-parallel-layer architecture has inverse degree at most 207. This
bound is independent of the base dimension, the number of helpers, and
the number of final parallel coordinates, and it is sharp.

This corollary imports **CubicCubeZeroSharpInverseDegree**, Section 1 of
[the retained sharp cube-zero inverse theorem](CUBIC_CUBE_ZERO_SHUFFLE_EXTENSION.md).
Its hypotheses are exactly a homogeneous cubic map over a characteristic-zero
field with \(J_0^3=0\); it gives polynomial invertibility of
\(F_0\) and \(D\le29\). Since \(B\le D\), (6) gives
\(\max\{D,7B+4\}\le207\). The literal sharp instance proves equality
can occur. The imported cube-zero theorem is a premise only of this
corollary; the direct bound (6), finite restart classification, and
no-\(S\)-reuse statement have their separate proofs below.

Equation (5) allows differences of duplicate helpers. For instance, two
coordinates with the same component \(P_S\) allow
\(Q=u^2(S-S')\). Their difference is fixed by the triangular map, and this
component has inverse value of degree three. A literal claim
\(Q_S=0\) for families with duplicate \(S\) helpers would be false.

## Three covectors give the complete condition

Put \(s=dwJ_0\) and
\[
 \omega=w^2dp_1-2uwdp_2+u^2dp_3.
\]
The third Jacobian rows on the nine old added coordinates are
\((\delta z)r+(\epsilon z)du\). They have zero entries in the three
\(p_i\) columns. The row \(\omega\) has entries only in those columns and
satisfies \(\omega J=0\), since
\[
 w^2d(u^2w)-2uwd(uw^2)+u^2d(w^3)=0.
\]
Consequently \(r,du,\omega\) are independent over the enlarged rational
function field: the first two are independent by (1), and the last has
nonzero fibre entries.

For clarity, the formulas producing every helper row can be checked without
an assumed continuation of the old derivations. Index the nine components
as \(P_z\) in their triangular order. Define
\[
 \eta z=\sum_{j\text{ moving}}(P_z)_j(P_j)_u,
 \qquad
 (F_z)_k=\sum_{j\text{ moving}}(P_z)_j(P_j)_k,
 \quad k\text{ moving}.
 \tag{7}
\]
The literal nine cubic formulas give
\(\delta z=(P_z)_w\), \(\delta P_z=0\), and
\[
 e_zJ^2=(\delta z)s+(\eta z)du+F_z.
 \tag{8}
\]
All \(F_z\) have entries only at \(p_1,p_2,p_3\). Multiplying by \(J\)
gives \(F_zJ=(\epsilon z)du\); its \(dw\) coefficient vanishes.
Thus (8) proves the old third-row formula over every base satisfying (1).
These are polynomial identities in independent coordinates, with no
localization at \(u\) or \(w\).

Define the old derivation
\[
 \gamma=9uw\partial_D-36w^2\partial_C
                    -9(6wp_2-4up_3)\partial_R,
 \qquad \eta u=\eta w=\gamma u=\gamma w=0.
 \tag{9}
\]
For \(Y\in V\), summing (8) against its gradient yields the complete helper
third row
\[
 dYJ^2=Y_wr+\eta(Y)du+\gamma(Y)\omega.
 \tag{10}
\]
Indeed the coefficient of \(s\) is \(\delta Y=0\), and the remaining
fibre row is \(\gamma(Y)\omega\). The latter identity is also derived in
[the third-covector proof](THIRD_COVECTOR_CUBIC_IMAGE_AND_RESTART_OBSTRUCTIONS.md).
The present checker reconstructs (7)--(9) directly from the literal cubic
components and checks all three fibre entries on a complete basis of \(V\),
which suffices for every cubic \(Y\) used here.

On \(A_0[y_1,\ldots,y_m]\), introduce **new** derivations
\(\mathcal D,\mathcal E,\mathcal C\) by
\[
\begin{array}{c|ccc}
 &\mathcal D&\mathcal E&\mathcal C\\\hline
 z\in A_0&\delta z&\epsilon z&0\\
 y_i&(Y_i)_w&\eta(Y_i)&\gamma(Y_i).
\end{array}
 \tag{11}
\]
In particular, \(\mathcal C\) is zero on \(A_0\); it is not the old
\(\gamma\). No commutation assertion is needed for (11).

The fourth row of a further component \(Q\) is exactly
\[
 (\mathcal DQ)r+(\mathcal EQ)du+(\mathcal CQ)\omega.
 \tag{12}
\]
Independence proves the necessary and sufficient condition
\[
 \mathcal DQ=\mathcal EQ=\mathcal CQ=0.
 \tag{13}
\]
Every first-layer helper is valid by (2), and every second-layer component
is valid precisely under (13). The row at \(p_1\) in the cube of the
Jacobian is \(u^2r\ne0\), so the final index is exactly four.

## Complete finite classification

Choose a weight-homogeneous basis \(Y_1,\ldots,Y_{21}\) of \(V\), with its
unique weight-207 member normalized to \(P_S\). The checker constructs this
basis by exact column elimination on all \(\binom{13}{3}=286\) cubic
monomials in the eleven old variables. The matrix of \((\delta,\epsilon)\)
has rank 265. Its full 21-element kernel, including the normalization of
\(P_S\), is retained in the replay receipt.

Let
\[
 K_d=K[A_0,y_1,\ldots,y_{21}]_d
       \cap\ker\mathcal D\cap\ker\mathcal E\cap\ker\mathcal C.
\]
Here \(K[A_0,y_1,\ldots,y_{21}]\) denotes the polynomial ring on the eleven
old coordinates and the 21 independent helper coordinates. There are 32
variables in total. The following is a complete exact computation, not a
sample or a bound on a subset of monomials.

**Lemma CompleteRestartKernelCertificate.** The dimensions and inverse
support bounds are:

| Ordinary degree \(d\) | All monomials | Rank | \(\dim K_d\) | Maximum exponent of \(W\) after inverse substitution | Maximum ordinary degree after inverse substitution |
|---|---:|---:|---:|---:|---:|
| 0 | 1 | 0 | 1 | 0 | 0 |
| 1 | 32 | 20 | 12 | 1 | 3 |
| 2 | 528 | 439 | 89 | 5 | 8 |
| 3 | 5,984 | 5,462 | 522 | 7 | 11 |

Every vector in each of these four kernels is independent of \(S=y_{21}\).
The inverse substitution in the table keeps \(U,W\) independent, sends the
nine old auxiliary coordinates to their triangular inverse values, and
sends every helper coordinate to its target coordinate plus its component
evaluated on those old inverse values. All auxiliary target coordinates
remain independent; none is set to zero.

For a reproducible, slightly more detailed certificate, the maximum total
degree outside \(W\), indexed by the exponent of \(W\), is
\[
\begin{array}{c|rrrrrrrr}
 d\backslash\deg_W&0&1&2&3&4&5&6&7\\\hline
1&3&0&&&&&&\\
2&6&5&4&3&4&3&&\\
3&9&8&7&6&7&6&5&4.
\end{array}
 \tag{14}
\]
A blank means no such monomial occurs.

The verification enumerates all degree-\(d\) monomials with repetition and
performs rational Gaussian elimination with the original-column
combinations retained. Each dependent column supplies a kernel vector with
coefficient one at its own new column and zero at later columns. These
vectors are linearly independent and account for every dependent column;
rank plus kernel dimension equals the complete domain size. Each resulting
vector is independently substituted into all three derivations.

For efficiency, the auxiliary weights (4) split the matrices into finite
blocks. The checker verifies that every coordinate image in (11) has
weight equal to its coordinate weight minus 29, 1, or 89, respectively.
Thus different input-weight blocks have disjoint row supports within each
tagged derivation and no relation is lost by the split. The weights do not
provide the inverse-degree conclusion: exact substitution and cancellation
in the full 32-variable target polynomial ring establish (14). The checker
also verifies both inverse compositions on all 32 coordinates.

All ranks, kernel vectors, support restrictions, and identities use
\(\mathbb Q\). Characteristic-zero scalar extension preserves ranks and
kernels, so the same conclusions hold over every field \(K\) in the theorem.

## Arbitrarily many helpers and the degree bound

A finite helper family need not be a basis. First add missing basis
components until its span is all of \(V\); this does not change the rows or
inverse coordinates already present. An invertible linear change among the
helper coordinates then gives the 21 basis components and additional
coordinates with zero component. Call the latter coordinates \(f_\alpha\).
They are fixed under the map and under all three derivations (11).

Expand a total homogeneous cubic \(Q\) as a polynomial in the \(f_\alpha\).
Each coefficient of a frozen-variable monomial of degree \(3-d\) lies in
\(K_d\), by (13). The finite classification therefore covers every finite
helper family, with no bound on its length. All four kernels ignore the
normalized basis coordinate \(S\), which transforms back into precisely
the direction (5).

The same expansion proves that the inverse value of any such \(Q\) has
\(W\)-exponent at most seven and ordinary total degree at most eleven,
including the frozen target variables. For the degree-\(d\) coefficient,
the maximum ordinary degree after multiplication by the frozen monomial is,
respectively, 3, 5, 9, or 11. First-layer helper inverse coordinates obey
the same envelope directly, as does the old nine-coordinate block.

Now substitute \(U=(G_0)_u=u\) and \(W=(G_0)_w\). Every other target
variable in this partial inverse has degree one. A monomial of total degree
\(t\le11\) and \(W\)-exponent \(b\le7\) consequently has degree at most
\[
 (t-b)+bB=t+b(B-1)\le11+7(B-1)=7B+4.
\]
The original base coordinates have inverse degree at most \(D\), proving
(6). Polynomial invertibility follows by successive substitution over the
assumed polynomial inverse \(G_0\); no formal-to-polynomial passage is left
unproved.

For sharpness, the retained literal base has \(D=B=29\) and
\(G_w(se_0)=c s^{29}\) with nonzero rational \(c\). With the auxiliary
targets zero, the inverse value of \(P_S\) is
\[
 -2U^4W^7=-2c^7s^{207}.
\]
The checker revalidates the literal base's positive weights, triangularity,
zero Jacobian cube, independence in (1), and this exact inverse coefficient.
Thus the bound is attained inside the stated family. Allowing the complete
first parallel helper space and one further arbitrary cubic parallel layer
cannot improve the existing lower bound 207.

## What this closes, and what remains open

The two previously excluded simple restarts selected one helper and one
quadratic coefficient at a time. The present theorem permits the complete
21-dimensional helper space simultaneously, arbitrary mixed coefficients,
products involving two or three helpers, and every old cubic correction.
Its conclusion comes from a complete finite polynomial space.

It does not permit helpers involving other coordinates of the base, change
the fixed nine pre-\(S\) components, or iterate a new layer after the final
\(Q\) coordinates. Those are different mathematical domains. The theorem
also does not assert whole-piece vanishing in an index-four Engel algebra;
that approach is excluded by the separate classical cyclic example.

The classical inverse interpretation is in Umirbaev, *Polarization algebras
and their relations*, MPIM 2014-60, Lemma 7, Lemma 8, and Corollary 2,
printed pages 7--8 (retained PDF pages 11--12). Those actual statements
identify nilpotent Jacobians with the Engel identity and polynomial inverse
termination with eventual vanishing of the weighted word sums. They do not
supply this finite restart classification or a dimension-uniform bound at
index four. The four-dimensional example on printed page 12 further warns
that nonzero individual words need not force nonzero inverse sums. The
primary source is context; the new matrix and substitution computations
are direct. No worldwide novelty, independent mathematical review, formal
verification, or external acceptance is asserted.

## Reproduction and evidence

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_cubic_index_four_complete_parallel_restart.py --repo-root .

The checker uses only the Python standard library and exact rational
arithmetic. It rebuilds the complete helper space and the degree-one,
degree-two, and degree-three restart kernels, verifies every kernel vector
and its full inverse substitution, and checks the sharp literal base. It
rejects an invalid helper, a direct forbidden reuse of \(S\), and the false
ceiling 206. A duplicate-helper positive control checks the precise scope
of (5).

The [source ledger](cubic_index_four_complete_parallel_restart_source_ledger.json)
and [replay receipt](cubic_index_four_complete_parallel_restart_verification.json)
accompany this proof. The deterministic monomial enumeration, pinned
literal cubic input, exact elimination, and stated support checks completely
specify the finite computation; no sampled tests or missing generated
certificate are premises.
