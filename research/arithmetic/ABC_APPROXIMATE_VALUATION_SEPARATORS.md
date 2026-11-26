# Fractional packings of approximate valuation separators

A prime projection can contribute to a counting bound even when it
fails to determine the complete letter or primitive ratio. The loss is
its conditional entropy. In particular a projection that leaves at
most logarithmic rank unseen has vanishing normalized loss. The
three radical constraints still enter the count through the
fixed-letter and two-letter inequalities. The hypotheses concern
primitive quotients, so arbitrary common multipliers play no part.

## Entropy with incomplete prime projections

Let \(Z\) be supported on positive reduced rationals of height at most
\(X^{H_0}\), with \(H_0\) fixed. For a prime set \(S\) let \(Y_S\) be
the vector \((v_p(Z))_{p\in S}\), and write
\(\mathcal R_S(Z)=\prod_{p\in S,\,v_p(Z)\ne0}p\).
The same entropy argument as in (VS4), now on the possible values of
\(Y_S\), gives
\[
             H(Y_S)\le \mathbb E\log\mathcal R_S(Z)
                                   +o_{H_0}(\log X).
\]
Indeed, on each projected radical there are only height-subpower
many signed valuation records. The inverse-radical weights on all
records have total at most \(X^{o_{H_0}(1)}\), since the projected
radical is at most \(X^{2H_0}\) and the harmonic sum over possible
radicals is \(O_{H_0}(\log X)\). The error is uniform in \(S\) and
in the law. [Separating valuation entropy and fixed-support enumeration](ABC_VALUATION_SPAN_COSETS.md#fractional-packing-of-separating-prime-sets)

For a finite collection of such sets take weights \(w_S\ge0\)
with \(\sum_{S\ni p}w_S\le1\) for each prime \(p\), and put
\(P=\sum_S w_S>0\). Repeated use of the entropy chain rule gives
\[
 H(Z)\le \frac1P\mathbb E\log\mathcal R(Z)
             +\frac1P\sum_S w_SH(Z\mid Y_S)+o_{H_0}(\log X).
                                                               \tag{AV1}
\]
The projection errors are averaged with weights \(w_S/P\).
The number of sets and the number of possible valuation records
may grow with height. The inequality also holds conditional on
a recorded datum, if the projected sets and their packing can be
chosen separately in each fibre and scaled to the same total
weight \(P\). Apply the uniform statement in each fibre and
average.

## Insertion into the three radical constraints

Use the uniform law on a family of positive primitive triples in
\(X/2<c\le X\) with \(\operatorname{rad}(abc)\le X^\lambda\).
Fix the small budget \(\lambda\) before the height limit and use
\(\lambda\log X\) as the unit of entropy. Let \(d\) be the normalized
logarithm of the family size and define the mean radical masses
\[
 t_W=\frac{\mathbb E\log\operatorname{rad}(W)}{\lambda\log X},
 \qquad W=A,B,C .
\]
They satisfy \(t_A+t_B+t_C\le1\). Suppose a packing of fixed
weight \(P\ge2\) has average loss at most \(\omega\) in these units,
in the sense that
\[
         \frac1P\sum_S w_SH(Z\mid Y_S)\le\omega\lambda\log X .
\]
For \(Z=T\), a complete letter, and for \(Z=U/V\), a primitive ratio
of two distinct letters, respectively, the count satisfies
\[
\begin{split}
 d &\le \min\left\{\frac{t_T}{P}+\frac{1-t_T}{4}+\omega,\
                              \frac{1-t_T}{2}\right\}
                                  +\eta+o_X(1),\\
 d &\le \min\left\{\frac{t_U+t_V}{P}+\omega,\
                              \frac12-\frac{t_U+t_V}{4}\right\}
                                  +\eta+o_X(1).
\end{split}                                                     \tag{AV2}
\]
Here any fixed \(\eta>0\) is allowed at sufficiently small fixed
budgets, and the errors in the height limit are uniform. In
particular on families with individual normalized radical caps
\(1/3+o(1)\) in the iterated balance limit, (AV1) and the fixed
fibre count give respectively
\[
              d\le \frac16+\frac{1}{3P}+\omega+o(1),
              \qquad d\le\frac{2}{3P}+\omega+o(1).
\]
A fixed positive gap below \(1/3\) in either bound gives a
uniform saving on the associated sector.

To prove (AV2), bin the three normalized radical masses on a
fixed mesh much smaller than \(\eta\), and denote the resulting
record by \(J\). Its entropy cost is \(O_\eta(1)\).
Writing \(\mathcal F\) for the random row, the chain rule gives
\[
 H(\mathcal F\mid T)\le H(J)+H(\mathcal F\mid T,J).
\]
In a fixed bin and at each fixed \(T\), (FL1) bounds the fibre
entropy by half the lesser of the two free radical caps, up to
the stated allowances. This is at most one quarter of their
sum. Averaging over \(T,J\) gives the bound
\((1-t_T)/4\) for the normalized fibre entropy, up to losses.
Applying (GD3) in each bin with the two free letters gives the
second bound \((1-t_T)/2\), again averaging and paying \(H(J)\).
For a ratio \(U/V\), primitivity recovers the two complete
letters and hence the row. It also gives
\(\mathcal R(U/V)=\operatorname{rad}(UV)\).
The other bound follows by averaging the two (GD3) inequalities
for the pairs \((U,T)\) and \((V,T)\). The projection loss in
(AV1) is used unconditionally. There is no requirement that its
bound hold separately in each radical bin.
[Fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

For \(0\le\omega<1/4\), maximizing (AV2) over the permitted mean
masses gives the following unconditional coefficients on the
packing sector, before the allowance \(\eta\):
\[
 \gamma_{\rm letter}(P,\omega)=
 \begin{cases}
    (2+2P\omega)/(P+4),&2\le P\le4,\\
    1/4+\omega,&P\ge4,
 \end{cases}
 \qquad
 \gamma_{\rm ratio}(P,\omega)=
     \min\{1/P+\omega,\ (2+P\omega)/(P+4)\}.
\]
The estimates specialize to (VS8) at zero loss. If each letter
admits its own packing of weight at least \(P\), with respective
losses \(\omega_W\), adding the three pair entropy inequalities
improves the coefficient to \(2/(3P)+(2/3)\sum_W\omega_W\).
If each pair ratio admits such a packing with losses
\(\omega_{UV}\), adding the three ratio entropy bounds gives
\(2/(3P)+(1/3)\sum_{\{U,V\}}\omega_{UV}\). Packings for different
letters or ratios need no disjointness across those choices.
[Exact packing bound (VS8)](ABC_VALUATION_SPAN_COSETS.md#counting-with-several-separating-sets)

## Missing valuation rank

Let \(V\) be the rational span of differences of valuation vectors
of \(Z\). For a projection onto \(S\) suppose its kernel on \(V\)
has dimension \(O((\log X)^\kappa)\) for fixed \(\kappa>0\). Set
\(\beta_\kappa=(1-1/\kappa)_+\). For every fixed
\(\sigma\in(\beta_\kappa,1)\),
\[
 H(Z\mid Y_S)\le \sigma\,\mathbb E
          \log\mathcal R_{\bar S}(Z)+o_{H_0,\sigma,\kappa}(\log X).
                                                               \tag{AV3}
\]
The bound on the kernel rank may vary with the selector subject to
a uniform constant. Given \(Y_S\), remove the observed prime
powers. The outside components still have bounded height. Their
valuation differences lie in the stated kernel, so some basis
of at most \(O((\log X)^\kappa)\) outside prime coordinates
determines the component in that fibre. Give each component
weight \(Q^{-\sigma}\), where \(Q\) is its radical on the chosen
basis. The total weight is subpower in height. In fact the
signed-valuation multiplicity for each \(Q\) is subpower, and
the Euler product for the possible \(Q\)'s is
\[
 \prod_{p\ {\rm in\ the\ basis}} (1+p^{-\sigma})
   \le \exp\!\left(O_\sigma\big((1+\dim\ker)^{1-\sigma}\big)\right)
   =X^{o(1)} .
\]
Relative entropy bounds the conditional entropy by
\(\sigma\mathbb E\log Q+o(\log X)\). Use
\(Q\le\mathcal R_{\bar S}(Z)\), and average the uniform bound.
[Coordinate bases and the coset count (VS1)](ABC_VALUATION_SPAN_COSETS.md#counting-low-rank-cosets-by-separating-valuations)

Combining (AV3) and the projected entropy inequality gives
\(H(Z)\le \sigma\,\mathbb E\log\mathcal R(Z)
 +(1-\sigma)\mathbb E\log\mathcal R_S(Z)+o(\log X)\) for one set.
If all selectors in a fractional packing of weight \(P\) satisfy
the kernel bound, average this inequality over the packing.
The result is the zero-loss entropy inequality with effective
weight arbitrarily close to
\[
                 P_{\rm eff}
                   =\frac{P}{1+(P-1)\beta_\kappa}.
\]
Whenever \(P_{\rm eff}>2\), it therefore gives the strict
saving of (VS8) at that weight up to arbitrary allowances.
For \(\kappa\le1\), the missing rank is free at the exponent
level. For instance, three disjoint selectors each leaving
\(O((\log X)^\kappa)\) dimensions still give a strict saving
for \(\kappa<4/3\). No requirement is imposed on the number
or sizes of the primes on which those residual dimensions
are supported.

## Obstruction to extracting cheap projections

The [polynomial evaluation model for target multiplicities](ABC_TARGET_MULTIPLICITY_POINTWISE_OBSTRUCTION.md)
also tests the extraction of valuation projections. Recall that
its rows are indexed by \(x\in G\), \(h\in H\), with labels
\(x,x+h,2x+h\), dimensions \(a=s-b,b\), and each of the \(s\)
slot symbols encoded by a prime from a disjoint list of size
\(q\). The primes have size \(\asymp_s q\log q\). There are
\(q^s\) rows and \(q^a\) distinct complete labels on each part.
The rows obey the joint radical and numerical-height bounds of
that construction. Here and below \(q\) tends to infinity with
\(s\) fixed before \(s\) tends to infinity.

Fix \(0\le\alpha<1\) and take an arbitrary \(S\) with
\(\sum_{p\in S}p^{-\alpha}=o(\log X)\). It contains at most
\(q^{\alpha+o_q(1)}\) of the model's primes. Let
\(k=\lfloor s/4\rfloor\). For a complete letter or a pair ratio,
consider rows whose projection on \(S\) has nonzero valuations
in at least \(k\) distinct slots. Choose one marked prime in
each selected slot. Given \(h\), the symbols of \(x\) at those
slots are independent and uniform since \(k\le a\), and the
other two letters' symbols at any one slot are affine in \(x\)
with unit coefficient. The probability of hitting the \(k\)
prescribed primes is \(q^{-k}\). Summing over marks and slots
counts at most \(q^{s-(1-\alpha)k+o_q(1)}\) such rows.

On the other rows, the number of projection records is at most
\(q^{\alpha k+o_q(1)}\) for a letter and
\(q^{2\alpha k+o_q(1)}\) for a ratio, as fewer than \(k\)
slots can be marked and each ratio slot uses at most two primes.
If the projection separates the values on a subfamily, these
records give at most \(q^{b+\alpha k+o_q(1)}\) rows for a letter
and at most \(q^{2\alpha k+o_q(1)}\) for a ratio. For a letter
the factor \(q^b\) is its complete fibre size. For a ratio
the primes on its two letters are disjoint and the ratio
determines their two polynomial labels, hence the row. For
\(s\ge8\), combining the two parts bounds any such separated
subfamily by
\[
                     q^{s-(1-\alpha)s/8+o_q(1)} .
\]
This bound includes subfamilies and prime sets chosen after
observing the model. A cover by prescribed such projections
with uniform partition-sum bound multiplies it by the number
of pieces. Every prime occurs with probability \(1/q\), the
correct occurrence scale up to logarithmic factors. In any
fixed target fibre with \(b>0\), a free-letter slot has the
same occurrence law, and this persists on remaining slots
after fixing fewer than \(b\) other slots in that letter.
[Prime occurrence entropy (EB7–EB8)](ABC_EFFECTIVE_RADICAL_BUDGET.md#occurrence-frequencies-and-anatomic-entropy-at-saturation)

The model admits exactly two disjoint separating sets for
each ratio, up to the vanishing fractional-packing margin
in (VS5). Choose all primes of \(a\) slots on its first
letter and \(b\) slots on its second. The first \(a\) symbols
give its entire polynomial, and the other \(b\) evaluations
recover \(h\). The complementary collection has \(a\) slots
on the other letter and \(b\) on the first, and also recovers
the row. The upper bound \(2+o_q(1)\) follows from ratio
entropy \(s\log q\) and total ratio radical
\((2s+o_q(1))\log q\). For one letter the analogous optimal
fractional packing is \(s/a+o_q(1)\), by averaging the
separating sets of all \(a\)-slot collections and using (VS5).

Valuation bases give a different count in this example.
For \(a\ge2\) the one-letter slot symbols are pairwise
independent. For \(b\ge1\), the \(2s\) slot symbols of any
pair ratio are pairwise independent as well. Across different
slots this follows by conditioning on \(h\); at one slot,
different letter symbols are invertible linear combinations
of the independent uniform pair \(x(t_j),h(t_j)\). A
linear relation on the one-hot valuation coordinates would
express a sum of functions of these symbols as a constant.
After centering each function, pairwise independence makes
the variance of the sum equal the sum of the variances.
Every function must therefore be constant. The valuation
difference span is the direct sum of the block sum-zero
spaces. Its rank is \(s(q-1)\) for one letter and
\(2s(q-1)\) for a ratio. Even allowing logarithmic-rank
deficiency, a linear coordinate selector must contain all
but a vanishing fraction of the relevant coordinate primes.
Thus a fractional packing of such selectors has weight
at most \(1+o_q(1)\).

The model uses genuine prime-powered integer labels and satisfies
the pointwise pool bounds, while its addition relation holds on
the polynomial labels. The exact equation on integer labels
remains an independent obligation in an extraction argument.
