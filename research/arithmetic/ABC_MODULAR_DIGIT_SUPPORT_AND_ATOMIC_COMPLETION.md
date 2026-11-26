# Modular digit supports and recovery of the atomic factors

A difference of adjacent binary digits of a valuation can be
read from three projected orders. Its numerator and denominator
are squarefree and have product at most the original radical
cap of the two letters. Thus many such differences can be fixed
at a modulus on the radical scale. At the balanced endpoint
one can use still smaller masks to recover the nonatomic
supports and then use a full divisor to release the remaining
atomic factors.

## Binary differences at the radical scale

Let \(U,V,T\) be a fixed ordering of positive primitive
original letters at most \(X^{B_0}\). Set \(n=X^L\) with \(L>0\)
fixed, and suppose that
\(\operatorname{rad}(U)\operatorname{rad}(T)\le n^\varrho\).
For a prescribed pool \(P\) of polynomial prime product and an
integer \(j\ge1\), put
\[
 X_j(W)=\prod_{p\notin P} p^{\lfloor v_p(W)/j\rfloor},
 \qquad R_j=X_j(U)/X_j(T).
\]
Choose a modulus \(D\) coprime to the outside parts of \(U,T\).
For \(j=0,\ldots,J\) consider the residues of
\[
                  M_j=R_{2^j}R_{2^{j+1}}^{-3}R_{2^{j+2}}^2
                       \pmod D .                           \tag{MD1}
\]
The record can be computed from the individual residues even when
all of them have been multiplied by a common unit, allowed to
depend on the row.

Write \(b_j(e)\in\{0,1\}\) for the \(j\)-th binary digit of \(e\).
The identity
\[
 \lfloor e/2^j\rfloor-3\lfloor e/2^{j+1}\rfloor
           +2\lfloor e/2^{j+2}\rfloor=b_j(e)-b_{j+1}(e)
\]
gives every prime valuation of \(M_j\), with a minus sign
for a prime of \(T\). The numerator and denominator are
coprime and squarefree. Their product is at most \(n^\varrho\).

**MD2.** If \(D>4n^\varrho\) and \(2^{J+1}>B_0\log_2 X\),
each record (MD1) gives at most \(X^\varepsilon\) original rows,
with any \(\varepsilon>0\), provided \(J=O(\log\log X)\).
For \(D>n^{2\varrho}\) and \(P\) empty, it gives at most one.

Indeed, consider one fraction of numerator \(A\) and denominator
\(B\) with \(AB\le n^\varrho\). Split \(A,B\) into dyadic boxes.
Within one box, the two possible cross products in the equality
of residues are at most \(4n^\varrho\), so the determinant must
vanish. At \(D>n^{2\varrho}\) the same reconstruction applies
without a box partition. Over \(O(\log\log X)\) fractions the
box choices number at most
\(\exp(O((\log\log X)^2))\).

For a prime outside \(P\), \(b_{J+1}(e)=0\) at every
occurring valuation. Hence the signed digits are recovered by
\[
                   \sigma_p b_j(e)=\sum_{i=j}^J v_p(M_i),
          \qquad \sigma_p=1 \text{ on }U,\quad -1\text{ on }T .
\]
This fixes the role and every original prime valuation outside
\(P\). Restore the \(P\)-parts at subpower cost on fixed support.
The original signed sum determines the third letter. The
argument does not use a modulus dividing an original letter.
For example \(P=\{2\}\) permits a prescribed power of two on
the radical scale for every row.
[Uniform fixed-support enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For an order list fixed at each \(L\), one can instead choose
\(J\) with \(2^{J+1}\ge L^{-2}\). If a positive valuation has no
visible outside bit difference, then its lowest \(J+2\) digits
are all equal, and it is at least \(2^{J+2}-1\). Thus the
unseen radical of the two letters costs at most
\(X^{2B_0/(2^{J+2}-1)}\), and support restoration gives
\(n^{O_{B_0}(L)}X^\varepsilon\) completions per record. This
simpler binary bank can use orders larger than \(1/L\).
The alternative logarithmic bank with non-dyadic probes
keeps every \(jL\) vanishing.
[Small-order bank (MR2)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md#a-logarithmic-list-of-twisted-residues-at-radical-precision)

## Completing all balanced profiles after small-modulus support recovery

Here take a height-first sequence of balanced boxes with limiting
degree measures \(\nu_U,\nu_T,\nu_V\), normalized total radical
of each letter \(1/3\), and pairwise sums of nonatomic masses at
least \(1/3\). Choose \(V\) with positive nonatomic mass. At
least two letters have positive nonatomic mass, so such a
choice is available. Use the exact-degree refinements and
fixed-support preparation in the mixed-residue results. The
finitely many prescribed exceptional coefficient and packed-root
factors may be frozen at vanishing normalized cost.
[Balanced degree measures](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

Fix any positive \(\eta,\delta\). For all sufficiently late
budgets in the sequence, there are row-dependent records of
two divisors \(D,D'\mid V\), with the following properties.
The numerical exponent of \(D'\) is bounded below by a positive
constant that can depend on \(\eta\) and the \(V\) profile.
It is an exact block of full prime powers with depths comparable
to \(1/L\). The divisor \(D\mid D'\) satisfies
\[
                       n^{\delta}<D\le n^{\delta+\eta}.
\]
Their joint dictionary has at most \(n^\eta X^\varepsilon\)
values after shrinking the raw price chosen for the block.
For \(D'\) choose an atomless compact carrier in the positive
degree distribution of \(V\), and freeze a block of mass
strictly between zero and \(\eta/2\). As in (AP1), selecting
root slots there has arbitrarily small overshoot, since the
maximum slot mass near that carrier tends to zero. The
combined numerical weight stays positive. Inside the prescribed
\(D'\), take \(D\) by multiplying its prime powers until the
threshold is crossed, using a partial final prime power if
needed. Each selected prime has normalized radical mass \(o(1)\),
so the partial choice increases \(\log_n D\) by at most \(o(1)\).
Alternatively enumerate the divisors of every \(D'\) at
divisor cost. These choices are made on original powers.
[Affordable full divisor selection (AP1)](ABC_AFFORDABLE_DIFFUSE_PROJECTED_PLANES.md#a-short-multiplier-from-the-third-radical)

**MD3.** For every \(\eta,\delta>0\), along such a balanced
sequence with \(L\) sufficiently late, a suitable list of
contrasts of projected ratios modulo \(D\), together with the
two divisor records, has fibre size
\[
                               O(n^\eta X^\varepsilon).
                                                               \tag{MD3}
\]
Use the bounded list of parity contrasts in (MR4), with each
admitted mask having total normalized radical at most
\(\delta/2\). The number of orders is bounded independently
of \(L\) along the limiting profile at fixed margins. It
costs only a fixed choice at each \(L\) to specify the list.
The row count
uses original rows and their exact powers. The residues
can be twisted by a common row-dependent unit.

For the proof choose finite lists of positive atoms of
\(\nu_U,\nu_T\) leaving omitted positive atomic mass less
than \(\eta/16\) in total. Retain the zero atom separately.
For each of the two nonatomic parts choose a compact
positive carrier disjoint from the atoms, leaving total
nonatomic mass less than \(\eta/16\) outside the carriers.
The mixed-support lemma applies with these carriers even
in the presence of atom masses. Its cheap masks have height
product at most \(n^{\delta/2+o(1)}\); the fixed residues
determine their rational values up to inner subpower box
choices. Their support union contains all but arbitrarily
small normalized mass of the two carriers. The prime pool
supplied by the recovered supports has polynomial product,
since every mask support is contained in the support of \(UT\).
Freeze the full parts of both \(U,T\) on that pool at inner
subpower cost. The pool may depend on the fixed transcript.

For the remaining primes of \(U,T\), enclose the finite atom
lists in shrinking neighborhoods and take a shrinking
interval at zero. Choose the widths to shrink sufficiently
slowly for weak convergence along the box sequence. Outside
these neighborhoods the radical mass, after removal of the
recovered mask supports, is at most \(\eta/4\) for late
budgets, by the finite atom and carrier choices and the support
lemma. Freeze those remaining outside parts using their radical
cap, at cost at most \(n^{\eta/4}X^\varepsilon\). All frozen
factors are units modulo \(D'\). For each retained positive
atom \(a\), combine its remaining squarefree roots into an
integer root of height \(X^{O(L)}\), raised to a representative
order with \(L j_a\to a\). The signed exponent corrections
have total numerical height \(X^{o(1)}\) because the window
widths shrink. The factors in the zero window also have
numerical height \(X^{o(1)}\). Combine these into a rational
wide correction in the ratio \(U/T\). The finitely many
narrow roots have disjoint supports. These are precisely
the amplified (GW) hypotheses with modulus \(D'\), including
the allowed overlap with the correction. The original signed
sum supplies an exact unit residue for \(U/T\bmod D'\).
The amplified count is \(n^{o(1)}X^\varepsilon\).
If there are no positive atom groups, rational reconstruction
applies directly to the short correction. An exact ratio
determines the coprime positive values \(U,T\) and then \(V\).
Combining the frozen costs proves (MD3), taking fixed inner
allowances small and \(L\) sufficiently late. The number of
atom windows is fixed after choosing \(\eta\), before the
small-budget passage.
[Diffuse support recovery with bounded probes (MR4)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md#a-bounded-number-of-small-modulus-probes)
· [Bounded root count and signed amplification (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

Thus a full-exponent balanced family needs a contrast-record
dictionary of exponent \(1/3-o(1)\), even at a divisor
of numerical size \(n^{\delta+o(1)}\), after paying arbitrarily
small fixed costs for the two divisors. For a fixed-power
exclusion it suffices to give a genuine joint dictionary saving
on these small-modulus contrasts in one orientation with a
nonatomic modulus letter. The support and completion arguments
apply to all balanced profiles in that orientation. The
radical caps alone have yet to supply the dictionary saving.
