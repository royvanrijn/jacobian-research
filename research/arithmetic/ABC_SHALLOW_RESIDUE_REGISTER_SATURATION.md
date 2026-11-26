# Shallow residue saturation in balanced abc charts

Exact valuation layers give a square-root count for *lists* of
unit residues already at radical-scale moduli. For a nearly
full-size complete-letter dictionary, the resulting collision
bound forces every sufficiently shallow projection to have
nearly maximal entropy in units of \(\log n\). At zero complete
multiplicity, three actual locks can consequently be truncated
to three residue alphabets of size \(n^{1/6+o(1)}\). Almost
\(n^{1/3}\) genuine rows persist with each pair of shallow
labels determining the row.

## Counting an arbitrary list of shallow residues

Fix \(L>0\), write \(n=X^L\), and consider distinct positive
complete values of the form
\[
                     U=c\prod_{i=1}^k x_i^{e_i},
        \qquad 1\le x_i\le n^{a_i},\qquad
              \rho=\sum_i a_i,\quad a_*=\max_i a_i .
\]
The positive fixed coefficient \(c\) can depend on \(X\).
The number of positions, the positive integer exponents and
the nonnegative caps are fixed in the height limit at each
budget. Positions of cap zero can be omitted. Assume the
values are units at an integer modulus \(2\le M\le X^B\)
and put \(\sigma=\log_n M\). For a chosen set \(I\) of
positive-cap positions let \(t=\sum_{i\in I}a_i\). The choice
\(I=\varnothing\) is allowed with \(t=0\).

**SRL1.** If \(I\ne\varnothing\) and \(t+\max_{i\in I}a_i<\sigma\),
then for any \(K\) prescribed unit residue classes the number of
complete values in those classes is at most
\[
                  O_\varepsilon\!\left(
                       K^{1/2} n^{\rho-t/2} X^\varepsilon
                                  \right).                  \tag{SRL1}
\]
The same estimate holds with \(t=0\). The count is uniform
in the coefficient, the modulus and the list. Its constant
may depend on the fixed budget, the box data, \(B\) and
\(\varepsilon\). For any fixed \(\kappa>0\), a greedy
sublist can give
\[
                  t\ \ge\
        \max\{0,\min(\rho,\sigma)-2a_*-\kappa\}.
\]

Fix all complementary roots at a cost at most \(n^{\rho-t}\).
Their product and \(c\) multiply the allowed classes by a
unit. For a selected position write
\(S_i(\chi)=\sum_{x_i}\chi(x_i^{e_i})\), where the sum can
use any allowed unit subalphabet. Give the character group
of \((\mathbb Z/M\mathbb Z)^\times\) uniform probability.
Take \(h_i=\lceil t/a_i\rceil\), so
\[
          \sum_{i\in I}h_i^{-1}\le1,\qquad
                       h_i a_i\le t+a_i<\sigma .
\]
The moment \(\mathbb E_\chi |S_i(\chi)|^{2h_i}\) counts pairs
of collapsed products whose \(e_i\)-th powers agree modulo
\(M\). The quotient of the products is one of the \(e_i\)-th
roots of unity modulo \(M\). There are \(X^{o_X(1)}\) such
roots uniformly for \(M\le X^B\), since a prime-power factor
contributes at most \(2e_i\), and \(\omega(M)=O_B(\log X/
\log\log X)\). A collapsed product is smaller than \(M\),
so for one partner and one torsion residue there is at most
one other product. Each product has only height-subpower
many ordered factorizations at fixed \(h_i\). The moment
norm is therefore at most \(n^{a_i/2}X^{o_X(1)}\).
Cauchy–Schwarz on the character sum of the residue list
gives \(K^{1/2}\) by orthogonality. Hölder on the selected
moments, with a factor \(1\) for any spare reciprocal
exponent, gives the remaining factor \(n^{t/2}X^{o_X(1)}\).
For the greedy statement use a sublist of total cap at most
\(\min(\rho,\sigma-a_*-\kappa)\) when that quantity is
positive. The loss at the last position is at most \(a_*\).

This is the character comparison behind the one-class
bound (RO3) with the residue-list norm restored. The
condition on the actual moment heights is
\(t+\max_{i\in I}a_i<\sigma\), which can be weaker than
\(2t<\sigma\) for a fine exact-layer box. A separate
elementary bound covers a slot with \(a_i\ge\sigma\).
After fixing the other slots, each prescribed power residue
permits at most \(X^{o_X(1)}(1+n^{a_i}/M)\) root choices.
In that case the list count is
\(O_\varepsilon(K n^{\rho-\sigma}X^\varepsilon)\).
[One-class shallow moments (RO3)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#a-residue-bound-at-several-depths)
[Positive numerical modulus and residue lists (RB2)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)

## Collision and conditional information bounds

Use entropy normalized by \(\log n\). For a random variable
\(Z\) with law \((p_z)_z\), put
\(H_{2,n}(Z)=-\log_n\sum_z p_z^2\).
Suppose the law of \(U\) satisfies the pointwise bound
\(\Pr(U=u)\le n^{-\alpha}\). This includes the uniform law on
at least \(n^\alpha\) distinct complete values.
In a uniform row law with \(N\) rows and maximum degree
\(D_U\) on that role, one can take \(\alpha=\log_n(N/D_U)\).
Write \(\Theta_M=U\bmod M\).

**SRL2.** With the \(t\) of (SRL1), the \(j\)-th largest
residue probability is at most
\[
                 O_\varepsilon\!\left(
                   j^{-1/2} n^{\rho-\alpha-t/2}X^\varepsilon
                                \right).
\]
In particular,
\[
 \begin{split}
   \Pr(\Theta_M=\Theta'_M)
      &\ll_\varepsilon n^{2(\rho-\alpha)-t}
                                    X^\varepsilon,\\
   H_{2,n}(\Theta_M)&\ge
                  t-2(\rho-\alpha)-\eta
 \end{split}                                                   \tag{SRL2}
\]
for any fixed \(\eta>0\) at sufficiently large heights at
each fixed budget. The prime on \(\Theta'_M\) denotes an
independent draw. The Shannon entropy obeys the same
lower bound. This follows by applying (SRL1) to the \(j\)
most probable classes and then summing \(j^{-1}\) up to
\(M\). The harmonic factor is height-subpower because
\(M\le X^B\). For example, if the collision probability is
at most \(n^{-\omega+\eta}\), residues of mass greater than
\(n^{-\omega+\xi}\) carry total probability at most
\(n^{-\xi+\eta}\). This last assertion follows by Markov's
inequality on the probability of the sampled residue.

An entropy version allows an arbitrary side record \(R\),
with no pointwise probability condition on the conditional
law. In each state let the fixed coefficient, residual
exact-layer caps and modulus be specified by the record.
Write \(\rho_R,t_R\) for the resulting parameters of (SRL1).
Assume the box data range over a fixed finite mesh at the
given budget, with uniform polynomial height bounds. Then
\[
    H(U\mid R)\le
        \mathbb E\bigl(\rho_R-t_R/2\bigr)
                     +\tfrac12 H(\Theta_M\mid R)+o(1) .        \tag{SRL3}
\]
The error includes an arbitrarily small normalized counting
allowance. In each state order the possible residues by
conditional probability and bin their ranks in dyadic
intervals. The bin index costs \(O(\log\log X)/\log n\) in
entropy. If a sampled residue has rank \(j\), its
conditional probability is at most \(1/j\). The expected
logarithm of the dyadic list size is thus bounded by the
conditional residue entropy with an \(O(1)/\log n\) error.
Apply (SRL1) to each bin and average. An actual full
supported part fixed on \(U\) can be stripped as a unit
coefficient before this argument. The collision conclusion
in a conditional fibre additionally needs a bound on its
maximum complete-value probability.
[Conditional entropy pricing on full divisors (LFR3)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#conditional-costs-of-exact-valuation-layers)

## Several actual probes at zero multiplicity

In a regular balanced cohort of actual rows, take a common
exact-layer box within the height-first conventions of
(RB1)–(RB5). Write \(d=1/3\), and let \(u\) be the common
complete-letter degree exponent. The row exponent is
\(d-o(1)\), the complete maxima of degrees are \(n^{u+o(1)}\),
and \(\rho_T=d+o(1)\) for each role. Suppose the maximum
exact-layer cap on the roles being used tends to zero.
With the usual tight dyadic cap refinement, a wholly
nonatomic limiting measure has this property after fixing a
tail of arbitrarily small normalized radical cost. Indeed a positive-size bounded scaled-order
layer along a subsequence would give a limiting atom.
Scaled orders tending to infinity carry arbitrarily little
root mass by the uniform numerical first-moment bound.
The argument can also use explicit same-order subdivisions
into small subroots when those subdivisions have the
stated total cap and small maximum cap.

Take the uniform law on the actual rows. At \(u=0\), for
any unit modulus \(M\le X^B\) with \(\sigma=\log_n M\),
(SRL2) gives
\[
                  H_{2,n}(U\bmod M)
                         =\min(d,\sigma)+o(1) .               \tag{SRL4}
\]
The upper bound uses that the variables are functions of the
row and have at most \(M\) possible values. Before passing to
zero multiplicity, a chart lower count \(n^{d-\gamma-o(1)}\)
and degree cap \(n^{u+o(1)}\) give the lower bound
\[
            \min(d,\sigma)-2u-2\gamma-O(a_*+\kappa)-o(1).
\]
Here the last allowance includes the root caps tending to
\(d\) and the arbitrarily small counting errors. Choose
the margins \(\kappa,\gamma\) and the chart allowances first,
then take sufficiently small fixed \(L\) and take the
height limit before decreasing those allowances. For a fixed finite list of coprime unit moduli, apply
(SRL4) to their products. Each sublist with total normalized
logarithm \(\sigma'\) has joint Shannon and collision entropy
\(\min(d,\sigma')+o(1)\). The same statement for a noncoprime
compatible list uses the least common multiple. In particular
coprime sublists with total below \(d\) have additive Shannon
entropies to this precision. For nested probes the conditional
Shannon entropy is the difference of the two capped depths.

For a prime-power terminal modulus one can make the
upper estimates pointwise on a retained set of terminal
residues. Suppose \(p\nmid U\) throughout the chart and
\(M_N=p^N\le X^B\). Put \(\sigma_j=\log_n p^j\),
\(q_j=\min(d,\sigma_j)\), and take \(\xi>0\) larger than
the common collision error, including the complete-degree
and chart losses. Apply (SRL2) at every \(p^j\), using
\(j=0\) for the trivial register. Outside a power-small
proportion of chart rows, their class probabilities at
every depth \(j\) are between
\(n^{-q_j-\xi}\) and \(n^{-q_j+\xi}\).
The upper tail follows from (SRL2). For the lower tail,
the number of available classes is at most the smaller of
\(p^j\) and the chart row count \(n^{d+o(1)}\).
There are \(O_B(\log X)\) depths, so the union of the
exceptional events is still power-small.

Use one representative per surviving residue at the
terminal depth and call the resulting set \(S\subset
(\mathbb Z/p^N\mathbb Z)^\times\). Then
\(|S|\ge (1-o(1))n^{q_N-\xi}\).
Each terminal class represented in a ball modulo \(p^j\)
has original probability at least \(n^{-q_N-\xi}\).
A nonempty such ball has original probability at most
\(n^{-q_j+\xi}\), so it contains at most
\(n^{q_N-q_j+2\xi}\) representatives. Thus
\[
 \frac{\#\{s\in S:s\equiv z\bmod p^j\}}{|S|}
                 \ \ll\ n^{-q_j+3\xi}
                 \ \le\ n^{3\xi}p^{-j q_N/\sigma_N}
                   \qquad (0\le j\le N).
\]
The last inequality uses \(q_j\ge \sigma_j q_N/\sigma_N\).
For a fixed positive terminal exponent this gives
nonconcentration at all \(p\)-adic scales with arbitrarily
small power allowances in the diffuse zero-multiplicity
limit. This is the upper ball bound in Shen's Definition 1.5.
[Prime-power nonconcentration, Definition 1.5](https://arxiv.org/pdf/2607.12251v1#page=4)

Suppose a cheap lock \(D_W\) on a marked role \(W\) is fixed
in the chart, with all its exponents at least \(\kappa_W/L\)
for fixed \(\kappa_W>0\). Suppose also its normalized raw
mass \(\mu_W=\log_n\operatorname{rad}D_W\) is bounded below
by a fixed positive constant and above by an arbitrarily
small fixed allowance. Such a lock is obtained from a
compact positive band of nonatomic mass. For a fixed probe
size \(\sigma_0>0\), choose
\[
                   M_W=\operatorname{rad}(D_W)^{j_W},
                 \qquad j_W=\lceil\sigma_0/\mu_W\rceil .
\]
The positive lower bound for \(\mu_W\) is fixed after
the lock tolerance is selected. It bounds \(j_W\) independently
of height and of later budgets, so for sufficiently small
fixed budgets \(M_W\mid D_W\). Its logarithmic size in base
\(n\) is in \([\sigma_0,\sigma_0+\mu_W]\). Its absolute register is the
common residue, up to sign, of the two opposite complete
values. Low fixed depth thus already accounts for the
full normalized information of a deep register when
\(\sigma_0\ge d\).
[Balanced dictionaries and original registers (RB3–RB7)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## A three-register partial Latin square

**SRL5.** Suppose all three original roles supply cheap
locks on one chart and have the preceding layer granularity.
Choose \(\sigma_0=d/2\) for each truncated modulus, and write
\(J_A,J_B,J_C\) for the three actual shallow registers.
The three moduli are pairwise coprime by primitivity on
every original row. Each pair of registers is determined
by the third complete role, and, after adjusting the
fixed signs, is exactly its reduction at the product of
the two selected moduli. At \(u=0\), in the indicated
limits,
\[
         H_{2,n}(J_T)=d/2+o(1),\qquad
         H_{2,n}(J_T,J_S)=d+o(1)\quad (T\ne S).
\]
The same identities hold for Shannon entropy. In particular
every pair has vanishing normalized mutual information
and determines the actual row up to vanishing normalized
conditional entropy. The three pairwise independence
assertions are at the exponent scale. Their joint row
entropy is only \(d+o(1)\).

For a cardinal version, retain actual rows whose pair
probabilities in the original uniform chart law are at most
\(n^{-d+2u+O(\delta)}\), using (SRL2). Here the raw lock
allowance, layer granularity, chart loss and counting
margin are chosen small in terms of a fixed \(\delta>0\).
Markov's inequality shows that the deletion is a power-small
proportion when the threshold includes a further fixed
\(\delta\) margin. Each pair of shallow labels now occurs
on at most \(n^{2u+O(\delta)}\) rows. Select a subfamily in
which no such pair repeats by greedily keeping one row and
removing rows sharing any of its three pairs. It has
\(n^{d-2u-O(\delta)}\) rows and all three shallow alphabets
have size at most \(n^{d/2+O(\delta)}\). At zero
multiplicity this is a partial Latin square with the full
pair-count exponent. Its rows are actual primitive sums
with their three original radical caps and all original
eligibility restrictions. Its vertices are the shallow
residue labels.

## Sharpness and the lifting requirement

The entropy exponents in (SRL4) are attained by actual
one-letter dictionaries. Use the prime-paletted model
(RO1) with all slots independently chosen. It gives
\(n^{d-o(1)}\) distinct tall integers with radical cap \(n^d\)
and wholly diffuse limiting valuation measure
\(d\,\mathrm{Unif}[2,4]\). The caps in one exact degree
tend to zero. Fix a prime \(p_0\) outside the growing
palettes and take \(M=p_0^{\lfloor \sigma\log_{p_0}n\rfloor}\).
Then (SRL4) matches the upper bound from the dictionary and
the number of residues. This sharpness concerns exponents.
For example, one can put every palette in the progression
\(1\bmod p_0^j\) for any \(j\) fixed at the budget. The prime
number theorem in that fixed progression preserves the
dictionary exponents but every resulting complete value
is \(1\bmod p_0^j\). The residue law need not be close to uniform
in total variation.
[Actual prime-paletted dictionaries (RO1)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#an-arithmetic-dictionary-from-finite-field-evaluations)

There is also a direct integer lifting example for the
register algebra. Fix positive \(w_T\) with \(\sum_Tw_T<1\)
and a large fixed integer \(j_0\), and set \(\mu=d/(2j_0)\).
Choose three distinct primes \(p_T\asymp n^\mu\), for example
in three fixed disjoint relative intervals. At each small
fixed \(L\), set \(e_T=\lfloor w_T/(L\mu)\rfloor\),
\(D_T=p_T^{e_T}\) and \(M_T=p_T^{j_0}\). Thus the locks have
arbitrarily small fixed positive raw masses, positive scaled
orders, and numerical weights \(w_T+O(L\mu)+o_X(1)\)
whose sum stays below one. Also \(M_T\mid D_T\) and
\(\log_n M_T=d/2+o_X(1)\).
A cyclic Latin square of order \(Q=\min_T\varphi(M_T)\)
maps injectively into the three shallow unit groups. Its
\(Q^2=n^{d+o_X(1)}\) rows can be lifted to integers with
all three \(D_T\)'s fixed, and with each full register
injective on the rows. For this, lift every shallow
symbol on its repeated rows to distinct full unit
residues. Each lift has \(D_T/M_T\) choices, more than
\(Q\) at sufficiently small fixed \(L\).

At each \(D_T\), prescribe the signed coordinates
\((A,B,-C)\) to be zero on the marked role and opposite
units on the other two. These classes determine \(A,B\)
modulo \(D_0=D_AD_BD_C\). Choose \(A\) in its class in
\([X/5,X/4]\) and \(B\) in its class in \([X/3,2X/5]\).
Every prime common to \(A,B\) lies outside \(D_0\).
For the part \(A'\) of \(A\) outside \(D_0\),
inclusion-exclusion on the \(B\)-progression gives about
\(X\varphi(A')/(D_0 A')\) coprime choices, with error
\(O(2^{\omega(A')})\). The main term exceeds the error
at large heights because \(D_0\le X^{\sum w_T+o_X(1)}\) and
\(\varphi(A')/A'\ge 2^{-\omega(A')}=X^{-o(1)}\).
Take such a \(B\) and set \(C=A+B\). All triples
are positive and primitive, with \(X/2<C<X\).
Their complete coordinates are distinct on each role
because each has two opposite shallow residues and the
Latin square uses each pair once. All three specified
deep locks remain. This construction imposes no small
radical cap on the complete coordinates.

Thus shallow register algebra and even exact full-lock
integer lifting can coexist with an entire Latin square.
The new necessary structure keeps the three radical
caps on a putative exceptional family. An improvement
still requires a restriction on lifting the dense shallow
register graph to these three capped complete values.
For example, Shen's restricted-graph expansion requires
one common ambient pair of sets and a multiscale
nonconcentrated unit-coefficient set. Shallow marginal
entropy supplies neither a dense restricted graph on
that common pair for many coefficients nor a sparse
output at the same terminal scale.
[Shen, Theorems 1.7–1.8](https://arxiv.org/pdf/2607.12251v1#page=4)
