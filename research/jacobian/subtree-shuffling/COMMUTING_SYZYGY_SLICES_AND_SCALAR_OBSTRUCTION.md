# Commuting syzygy derivations and the next scalar-extension obstruction

Work over a characteristic-zero field \(K\), and write \(L=K(u,w)\).
Use the cubic extension and coordinate names
\(p_1,p_2,p_3,A,B,E,D,C,R,S\) of
[Cubic index four, inverse degree 207](CUBIC_INDEX_FOUR_SYZYGY58.md).
All statements below concern that specified coordinate algebra. They do not
classify arbitrary index-four cubic maps.

## The two pre-\(S\) derivations

Let \(J_0^3=0\) for the base, with zero component at \(u\), and put
\(r=dw\,J_0^2\). On the polynomial algebra
\[
 \mathcal A=L[p_1,p_2,p_3,A,B,E,D,C,R]
\]
define derivations \(\delta,\epsilon\), fixing \(L\), by the table
\[
\begin{array}{c|c|c}
z&\delta z&\epsilon z\\ \hline
p_1&u^2&0\\
p_2&2uw&0\\
p_3&3w^2&0\\
A&3up_2-18wp_1&0\\
B&6wp_2-2up_3&0\\
E&3up_2&0\\
D&-9p_1p_2-uA&15u^2w^4\\
C&12uB+8wA+24p_1p_3-21p_2^2&-24uw^5\\
R&6p_2A+24p_1B&
  f=6w^3(-21uwp_2+10u^2p_3+12w^2p_1).
\end{array}
\]
The third Jacobian row at each pre-\(S\) new coordinate \(z\) is
\((\delta z)r+(\epsilon z)du\). The actual literal base has \(r,du\)
independent over its rational function field, as certified in the degree-207
proof and checker.

**Lemma PreSScalarSyzygyCriterion.** For that actual base, appending a scalar
polynomial
\(Q\in K[u,w,p_1,p_2,p_3,A,B,E,D,C,R]\) preserves fourth Jacobian power zero
if and only if \(\delta Q=\epsilon Q=0\).
This is the exact block-matrix identity
\(dQ\,J^3=(\delta Q)r+(\epsilon Q)du\), and independence proves the converse.
The assertion is restricted to polynomials in these displayed coordinates,
not polynomials involving the other base variables.

Both derivations are locally nilpotent: their coordinate images are
triangular, and each coordinate is killed by some iterate. They commute.
Indeed \(\epsilon\) kills all variables occurring in \(\delta z\);
\(\delta(\epsilon z)=0\) is immediate for \(D,C\), and for \(R\) it is the
identity
\[
 \delta f=6w^3(-42u^2w^2+30u^2w^2+12u^2w^2)=0 .
\]
Checking the commutator on the generators proves it vanishes on all of
\(\mathcal A\). Neither derivation is asserted to have a uniformly bounded
nilpotency exponent on all polynomials.

## Complete localized normal forms

Set \(t=p_1/u^2\), so \(\delta t=1\). The finite projection
\[
 \pi(z)=\sum_{k\ge0}\frac{(-t)^k}{k!}\delta^kz
\]
on any generator is well defined by local nilpotence. It satisfies
\(\delta\pi(z)=0\); this follows by differentiating the sum and cancelling
successive terms. Define
\[
\begin{aligned}
 I_2&=p_2-2wp_1/u,& I_3&=p_3-3w^2p_1/u^2,\\
 I_A&=A-3p_1p_2/u+12wp_1^2/u^2,\\
 I_B&=B-6wp_1p_2/u^2+2p_1p_3/u+3w^2p_1^2/u^3,\\
 I_E&=E-3p_1p_2/u+3wp_1^2/u^2 .
\end{aligned}
\]
These are the respective \(\pi\)-images. Write
\(J_D=\pi(D),J_C=\pi(C),J_R=\pi(R)\). In particular
\[
 J_D=D+p_1A/u+3p_1^2p_2/u^2+4wp_1^3/u^3,\qquad
 \epsilon J_D=15u^2w^4 .
\]
Since the derivations commute and \(\epsilon t=0\), the projection commutes
with \(\epsilon\). Thus
\(\epsilon J_C=-24uw^5\) and \(\epsilon J_R=f\).
Here \(f\) is already \(\delta\)-invariant and is expressed in the new
coordinates as
\(6w^3(-21uwI_2+10u^2I_3)\).

**Theorem JointLocalizedSyzygyNormalForm.** The common kernel is the polynomial
algebra over \(L\) on the following seven algebraically independent generators:
\[
 \ker\delta\cap\ker\epsilon
 =
 L\left[
 I_2,I_3,I_A,I_B,I_E,\
 J_C+\frac{8w}{5u}J_D,\
 J_R-\frac{f}{15u^2w^4}J_D
 \right].
\]

To prove completeness, first the triangular coordinate change
\[
 (p_1,p_2,p_3,A,B,E,D,C,R)
 \longleftrightarrow
 (t,I_2,I_3,I_A,I_B,I_E,J_D,J_C,J_R)
\]
is invertible over \(L\): each normalized coordinate is its original
coordinate plus a polynomial in preceding coordinates, with nonzero
constant diagonal coefficient. In these coordinates \(\delta=\partial_t\).
Now set \(v=J_D/(15u^2w^4)\). The remaining triangular change replaces
\(J_D,J_C,J_R\) by \(v\) and the two normalized joint generators displayed
above. It makes \(\epsilon=\partial_v\), while all seven joint generators
and \(t\) are fixed by \(\epsilon\). Over characteristic zero the kernels
of the two partial derivatives are exactly the polynomials independent of
\(t,v\). This proves both completeness and algebraic independence.

Consequently the ordinary polynomial scalar extensions in the lemma are
exactly the intersection of this localized polynomial algebra with
\[
 K[u,w,p_1,p_2,p_3,A,B,E,D,C,R].
\]
Polynomiality at \(u=0,w=0\), and ordinary cubic degree, remain necessary.
A rational invariant with a parameter denominator is not automatically a
valid polynomial cubic extension.

For the earlier five moving coordinates alone one obtains the simpler
complete kernel
\[
 \ker\delta\big|_{L[p_1,p_2,p_3,A,B]}=L[I_2,I_3,I_A,I_B].
\]
On the zero-auxiliary inverse input, with \(w=W\), the four invariants have
values
\[
 (I_2,I_3,I_A,I_B)=(-uW^2,-2W^3,4u^2W^3,0).
\]
Thus a zero-value invariant is present; omitting every such compensator from
a subsequent search is not justified.

## A genuinely new covector at \(S\)

The preceding two-derivation description is a pre-\(S\) statement.
For the degree-207 extension, direct exact multiplication gives
\[
\begin{aligned}
 (J^3)_{S,p_1}&=-108w^4(up_2-2wp_1),\\
 (J^3)_{S,p_2}&=216uw^3(up_2-2wp_1),\\
 (J^3)_{S,p_3}&=-108u^2w^2(up_2-2wp_1).
\end{aligned}
\]
Every earlier third-power row has zero entries in these fibre columns.
In particular this row is not in the span of \(r,du\). Its fibre portion is
\[
 -108w^2(up_2-2wp_1)\,
 \omega,\qquad
 \omega=w^2\,dp_1-2uw\,dp_2+u^2\,dp_3 .
\]
The complete identity \(\omega J=0\) follows directly from
\(P_1=u^2w,P_2=uw^2,P_3=w^3\), and is rechecked against the full 58-variable
matrix. No assertion that these displayed entries give the entire row is
needed.

**Theorem PostSScalarExtensionObstruction.** Any single appended polynomial
component \(Q\) whose enlarged Jacobian still has fourth power zero must
satisfy \(Q_S=0\). This holds for the displayed ten-coordinate extension
over any base independent of the new variables; it does not use the
particular base's \(r,du\) independence.

The \(p_1\)-column of the new last fourth-power row is
\[
 Q_S[-108w^4(up_2-2wp_1)].
\]
Its second factor is a nonzero polynomial in independent coordinates, and
the coefficient ring is a domain. Hence \(Q_S=0\), and characteristic zero
makes \(Q\) independent of \(S\).
The three column formulas and their uniqueness depend only on the displayed
ten added cubics: an arbitrary base block has zero columns at the new
variables. This proves the stated base-independent scope.

**Corollary PostSCompleteScalarCriterion.** For the actual base and
\(Q\in K[u,w,p_1,p_2,p_3,A,B,E,D,C,R,S]\), a scalar extension preserves
fourth Jacobian power zero if and only if \(Q_S=0\) and \(Q\), regarded
as a polynomial in the pre-\(S\) algebra, belongs to the ordinary polynomial
intersection of the seven-generator common kernel above.
Necessity first follows from the post-\(S\) column obstruction and then from
the pre-\(S\) criterion. Conversely, if these conditions hold, \(dQ\)
annihilates every pre-\(S\) third-power row, which is unchanged by adjoining
\(S\), and it has zero coefficient at the new row. Thus the last fourth-power
row vanishes. This argument does not define an extended two-derivation table
on \(S\).

The condition \(Q_S=0\) by itself is insufficient:
\(Q=p_1^3\) has \(\delta Q=3u^2p_1^2\ne0\).
The actual new fourth-power row has coefficient
\(1179049050096099039/3404635167332\) in column one at
\(x_0^3x_1^2x_{30}p_1^2\), verified directly.

The obstruction permits parallel additions whose new third-power rows can
cancel. It supplies neither an unbounded family nor a finite bound for all
cubic index-four maps. In particular one must not keep using the same two
derivations after adjoining \(S\).

## Classical comparison and exact evidence

The finite slice transformations above are proved directly. Their first-stage
affine-linear derivation, homogenized by a fixed coordinate over \(L\), has
Jordan block sizes \(3,2,1\): adjoining a fixed coordinate \(h\) replaces
the three constant images by \(u^2h,2uwh,3w^2h\).
Its image is spanned by \(h,3up_2-18wp_1,6wp_2-2up_3\), of dimension three;
its square has one-dimensional image since
\(\delta^2A=-12u^2wh\), \(\delta^2B=6uw^2h\), and its cube is zero.
It fits the general Weitzenbock setting discussed
by Drensky and Makar-Limanov, arXiv0804.2933v1, Introduction, pages1--2, and by
Bedratyuk, *Serdica Math. J.*36 (2010), pages99--102, Sections1--2.
Drensky--Makar-Limanov's explicit determinant generating theorem is restricted
to Jordan blocks of size two; that special theorem does not describe this
mixed-block action. Bedratyuk's classical binary-covariant correspondence
concerns linear locally nilpotent actions, not an unproved infinite
prolongation of the new cubic coordinate blocks.

The attached standard-library checker expands the Laurent polynomial
derivations, commutator, slices and normal forms exactly. It also reconstructs
the entire literal Jacobian from the pinned canonical degree-207 input and
checks all three new-column formulas, column uniqueness and \(\omega J=0\).
The retained polynomial intersection and infinite-iteration obligations are
mathematical scope limits, not computations the checker has solved.


## Reproduction

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_commuting_syzygy_slices.py \
      --output /tmp/commuting_syzygy_slice_review.json

The [inputs](commuting_syzygy_slice_inputs.json),
[source ledger](commuting_syzygy_slice_source_ledger.json), and
[exact receipt](commuting_syzygy_slice_verification.json) accompany this proof.
The ledger also pins the existing literal degree-207 map and its proof.
All 522 pre-S row/column identities are compared with the literal matrix;
the negative converse uses the actual nonzero fourth-power row of p1 cubed.
The two primary PDFs are retained under sources/ as context, not imported
premises of the slice proof. No external human or formal review is claimed.

Updated 20 November 2026.
