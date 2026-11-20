# A cubic index-four map with inverse degree 207

**Theorem CubicIndexFourSyzygy58.** Over every characteristic-zero field \(K\), there is a homogeneous cubic map
\(H:K^{58}\to K^{58}\) for which
\[
 (JH)^4=0,\qquad (JH)^3\ne0,\qquad
 \deg(\operatorname{id}-H)^{-1}=207.
\]
The complete literal cubics and an exact rational checker are attached.
This finite example implies that any dimension-independent inverse-degree
bound for cubic maps with fourth Jacobian power zero is at least 207.
It does not prove an unbounded family or settle the fixed-index subtree
shuffle conjecture.

## The map

Take the literal 48-coordinate cubic map \(H_0\) reproduced in the attached
input, and put \(u=x_0,\ w=x_{47}\). Adjoin ten coordinates
\(p_1,p_2,p_3,A,B,E,D,C,R,S\), with cubic components
\[
\begin{aligned}
 P_1&=u^2w,&P_2&=uw^2,&P_3&=w^3,\\
 P_A&=u^2p_3+3uwp_2-9w^2p_1,\\
 P_B&=3w^2p_2-2uwp_3,\\
 P_E&=3uwp_2-2u^2p_3,\\
 P_D&=3up_2^2-9wp_1p_2-uwA,\\
 P_C&=12uwB+4w^2A+24wp_1p_3-21wp_2^2,\\
 P_R&=A(6wp_2-4up_3)+B(-12up_2+24wp_1)
                         +9p_2^3-12p_1p_2p_3,\\
 P_S&=-2uwR+3(up_2-2wp_1)C+(8up_3-12wp_2)D\\
    &\quad-8uAB-2wA^2-8wE^2-12p_1p_2B-8p_1p_3A
                                      +9p_2^2A+12p_2^2E .
\end{aligned}
\]
Every component uses only earlier coordinates. The complete map has 1,248
nonzero cubic monomials. All coefficients of the added components are integers.

The exact input contains the 48 original coordinates as literal cubics,
rather than importing a quotient or a construction algorithm. The checker
reconstructs the full polynomial Jacobian \(J\), multiplies it over
\(\mathbb Q[x_0,\ldots,x_{57}]\), and verifies \(J^4=0\). It also checks
\[
 [x_0^2x_2^3x_{10}](J^3)_{48,0}=36,
\]
so the nilpotency index is exactly four. These are polynomial identities and
a nonzero rational coefficient; they remain valid over every field of
characteristic zero.

## Inverse degree

The old positive coordinate weights have unique minimum one at \(x_0\),
maximum twenty-nine at \(x_{47}\), and the added weights are
\[
 (31,59,87,89,117,89,119,147,177,207).
\]
The input and checker validate every coordinate weight and every monomial:
its output weight equals the sum of its three input weights, and all input
indices are smaller than the output index. Hence
\(F=\operatorname{id}-H\) is a triangular polynomial automorphism. Recursive
substitution produces a polynomial inverse \(G\) preserving these weights.
All weights are positive and at most 207, so \(\deg G\le207\).

The exact inverse line of the 48-coordinate base is
\[
 G_{47}(s e_0)=c\,s^{29},\qquad
 c=\frac{8803690708443294969}{95329784685296}\ne0.
\]
This is rechecked by literal substitution in the present checker.
Put \(U=s\), \(W=c s^{29}\). With all ten new target coordinates zero,
their inverse values are respectively
\[
 U^2W,\ UW^2,\ W^3,\ -5U^2W^3,\ UW^4,\
 U^2W^3,\ -U^3W^4,\ -5U^2W^5,\ -U^3W^6,\ -2U^4W^7.
\]
They follow in order by substituting in the displayed cubics.
The last value has nonzero degree-207 coefficient \(-2c^7\).
Thus \(\deg G=207\).

**Lemma Syzygy58LineCriterion.** For this weight profile and any
weight-preserving triangular automorphism, inverse degree 207 is equivalent to a
nonzero coefficient of \(x_0^{207}\) in the final inverse coordinate,
equivalently a nonzero degree-207 restriction to \(s e_0\).
The reason is that \(x_0\) is the unique weight-one input and the final
coordinate is the unique weight-207 output; every other monomial or output
has ordinary degree below 207. This is a genuine instance-level equivalence.
The consequence that every dimension-independent upper bound is at least
207 is one-way.

## The parallel cancellation

The following identities explain the choice of the extra cubics. They are
not a substitute for the complete literal multiplication above.
For a cube-zero base write \(J_0=JH_0\),
\(r=dw\,J_0^2\) and \(s=dw\,J_0\); here \(du\,J_0=0\).
Before adjoining \(E,D,C,R,S\), the third-power rows on the new coordinates
\(p_1,p_2,p_3,A,B\) are
\[
 r\,(u^2,\ 2uw,\ 3w^2,\ 3up_2-18wp_1,\ 6wp_2-2up_3).
\]
**Lemma CubicIndexFourSyzygyExtension.** For any homogeneous cubic base
with \(J_0^3=0\) and zero component at \(u\), adjoining only
\(p_1,p_2,p_3,A,B,R\) with the displayed cubics preserves fourth-power
Jacobian vanishing. Thus the six-coordinate construction is a universal
method lemma, rather than a second principal finite-degree example.

The component \(P_R\) is killed by the polynomial derivation having precisely
these five coefficients. All base third-power rows are zero, and the last
fourth-power row is that derivation applied to \(P_R\), multiplied by
\(r\). This proves the lemma. Its own third-power row is
\[
 (6p_2A+24p_1B)r+f\,du,\qquad
 f=6w^3(-21uwp_2+10u^2p_3+12w^2p_1).
\]
The new third-power rows of \(E,D,C\) have \(du\)-coefficients
\(0,\ 15u^2w^4,\ -24uw^5\), respectively.
The three partial derivatives of \(P_S\) at \(D,C,R\) therefore give the
cancellation
\[
 (8up_3-12wp_2)(15u^2w^4)
 +(3up_2-6wp_1)(-24uw^5)-2uwf=0.
\]
The remaining \(r\)-coefficient cancels as well. Thus several parallel rows
permit reuse of \(R\). For the literal 48-variable base, a single scalar
extension of the preceding six-coordinate construction must instead be
independent of \(R\). Indeed \(r\) and \(du\) are independent over the
rational function field: the column-one entry of \(r\) has coefficient
\(393016350032033013/3404635167332\) at \(x_0x_1^2x_{30}\).
An added scalar \(Q\) has fourth-power row
\(\alpha r+Q_Rf\,du\); independence and \(f\ne0\) force \(Q_R=0\).
This last obstruction is for the specified base, not for every cube-zero
base in the universal method lemma.

There is no proved repeatable amplification rule here. A finite collection
of polynomial syzygies or growing finite inverse degrees does not establish
an infinite family.

## Evidence and context

The standard-library checker uses only rational arithmetic. It validates the
full data schema, strict triangularity and positive weights; reconstructs all
Jacobian products; checks the exact nonzero cube coefficient and inverse line;
and rejects malformed data and a literal coefficient perturbation making the
fourth power nonzero. The proof, literal input and checker are hash-pinned by
the source ledger. No CAS, tree quotient, or old checker is a premise.

Polarization, the Engel/Jacobian correspondence, and inverse sums of
nonassociative words are classical. See Umirbaev, *Polarization algebras and
their relations*, MPIM 14-60, Lemmas 7--8 and Corollary 2, printed pages 7--8.
Its example on printed pages 9--10 also explains why non-nilpotence of
the full polarization algebra need not give arbitrarily high inverse sums.
The present finite proof imports no theorem from that paper and makes no
novelty claim about the general method.

Run `python3 research/jacobian/subtree-shuffling/verify_cubic_index_four_syzygy58.py --output /tmp/cubic_index_four_syzygy58_review.json`.

The [source ledger](cubic_index_four_syzygy58_source_ledger.json),
[exact replay receipt](cubic_index_four_syzygy58_verification.json), and
[literal cubic input](cubic_index_four_syzygy58.json) accompany this proof.
The [primary context source](sources/umirbaev-polarization-algebras-2014.pdf)
is retained for the exact locators above; it is not a premise of the finite
polynomial identity checks.

Updated 20 November 2026.
