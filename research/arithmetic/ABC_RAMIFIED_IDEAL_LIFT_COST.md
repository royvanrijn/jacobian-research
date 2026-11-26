# Discriminant and counting costs of ramified powerful lifts

A ramified extension raises valuations of the original letters while
preserving coprimality of their principal ideals. The root discriminant
records the decrease in the normalized ideal radical. An explicit field
attains the minimum discriminant cost up to a factor depending only on
the prescribed orders. In a near-full family, even fields of unbounded
degree that achieve a fixed radical reduction have to vary over a
power-sized dictionary. A separate fibre bound quantifies this field
cost for minimally repaired targets.

## The local cost and a field attaining it

Let \(a+b=c\) be a positive primitive triple, and put \(R=\operatorname{rad}(abc)\).
For a number field \(F\) of degree \(d_F\), write
\[
 \operatorname{rd}(F)=|\operatorname{Disc}(F)|^{1/d_F},\qquad
 R_F=\left(N_{F/\mathbb Q}\prod_{\mathfrak p\mid (abc)}\mathfrak p\right)^{1/d_F}.
\]
Here \(e_{\mathfrak p},f_{\mathfrak p}\) are the ramification index and
residue degree over a rational prime \(p\). Put
\[
 s_F(p)=\frac1{d_F}\sum_{\mathfrak p\mid p}f_{\mathfrak p},\qquad
 \delta_F(p)=1-s_F(p).
\]
Thus \(\log R_F=\sum_{p\mid R}s_F(p)\log p\). The local different
exponent is at least \(e_{\mathfrak p}-1\), with equality in the tame
case. The identity \(\sum_{\mathfrak p\mid p}e_{\mathfrak p}f_{\mathfrak p}=d_F\) gives
\[
 \sum_p\delta_F(p)\log p\le \log\operatorname{rd}(F),\qquad
 \operatorname{rd}(F)R_F\ge R.
 \tag{RL1}
\]
The second inequality is an equality if all primes on \(R\) are
tamely ramified or unramified and the field is unramified away from
\(R\). These are the local contributions to the normalized truncated
counting function and discriminant. [Stevenhagen, Theorem 4.13](https://websites.math.leidenuniv.nl/algebra/localfields.pdf) · [Vojta, §1](https://swc-math.github.io/aws/1998/98Vojta.pdf)

Prescribe positive integer orders \(M_a,M_b,M_c\). At each \(p\mid R\)
let \(T\in\{a,b,c\}\) be its unique supported letter and put
\[
 m_p=\left\lceil \frac{M_T}{v_p(T)}\right\rceil,\qquad
 C_{\mathbf M}=\prod_{p\mid R}p^{1-1/m_p}.
\]
The principal ideal \((T)\) is \(M_T\)-full if every nonzero prime-ideal
valuation is at least \(M_T\). Any \(F\) making the three ideals full
at the prescribed orders has \(e_{\mathfrak p}\ge m_p\) at every
place above \(p\mid R\). It follows that \(s_F(p)\le1/m_p\) there and
\[
 \operatorname{rd}(F)\ge C_{\mathbf M},\qquad
 C_{\mathbf M}\ge R\prod_T T^{-1/M_T}.
 \tag{RL2}
\]
The second bound uses \(1/m_p\le v_p(T)/M_T\) and also allows orders
equal to one. The ideals remain pairwise coprime over every field.

Let \(d=\operatorname{lcm}_{p\mid R}m_p\). For \(d=1\) take \(K=\mathbb Q\).
Otherwise take \(K=\mathbb Q(\alpha)\) with
\[
 \alpha^d=A,\qquad A=\prod_{m_p>1}p^{d/m_p},\qquad \alpha>0.
\]
For every place of \(K\) above \(p\mid A\), the valuation of
\(\alpha\) is \(e_{\mathfrak p}/m_p\), so \(m_p\) divides
\(e_{\mathfrak p}\). Summing \(e_{\mathfrak p}f_{\mathfrak p}\) over that
prime shows that \(m_p\mid[K:\mathbb Q]\). Taking the least common
multiple gives \([K:\mathbb Q]=d\), since \(\alpha\) solves a degree-\(d\)
polynomial. This also proves fullness at all primes, including primes
dividing \(d\).

If \(p\nmid d\) and \(p\mid A\), pass to a finite unramified local extension
containing the \(d\)-th roots of unity and a \(d\)-th root of the unit in
\(A=p^{d/m_p}u\). These roots exist after a finite residue-field extension
and Hensel lifting. All the roots of \(Z^d-A\) then lie in the extension
obtained by adjoining \(p^{1/m_p}\). That last extension has ramification
index \(m_p\) by Eisenstein, and unramified base change preserves the
ramification index. Hence \(e_{\mathfrak p}=m_p\) at every place above \(p\).
For \(p\nmid dA\) the same unit argument gives unramifiedness. At primes
dividing \(d\), the polynomial discriminant \(\pm d^d A^{d-1}\) bounds
the field discriminant exponents in terms of \(d\), since \(v_p(A)\le d\).
Consequently
\[
 \operatorname{rd}(K)\asymp_d C_{\mathbf M},\qquad
 R_K\asymp_d R/C_{\mathbf M}.
 \tag{RL3}
\]
At fixed orders \(d\) is bounded by
\(\operatorname{lcm}(1,\ldots,\max_T M_T)\). In particular the constants
in (RL3) are fixed before the height tends to infinity. At equal orders
\(M\ge K_0/\lambda\), (RL2) costs at least
\(R c^{-3\lambda/K_0}\) in root discriminant. For lifting only the target,
take the other orders to be one. The corresponding bound is
\(\operatorname{rad}(c)c^{-\lambda/K_0}\).

## Counting with the full fractional ramification weights

Fix a function \(g(X)=o(\log X\log\log X)\) with \(g(X)\ge0\).
At height \(X\) the field \(F\) may vary with \(X\) and its degree may
be unbounded. Assume \(\log\operatorname{rd}(F)\le g(X)\). For every
fixed \(\eta>0\) and sufficiently small fixed \(\lambda>0\), the number
of positive primitive \(a+b=c\) in \(X/2<c\le X\) with
\(R<c^\lambda\) and \(R_F\le X^{b_0}\) is
\[
 \ll_{g,\lambda,\eta,\varepsilon} X^{b_0/3+\eta\lambda+\varepsilon},
 \qquad 0\le b_0\le\lambda.
 \tag{RL4}
\]
Constants are uniform over the fields and \(b_0\). In particular a
discount \(R_F\le R X^{-\chi\lambda}\), for fixed \(0<\chi<1\), costs
\(\chi\lambda/3\) in the counting exponent up to \(\eta\lambda+\varepsilon\).

Use the prime price \(w_X(p)=s_F(p)\log_X p\). Convexity on
\(0\le\delta_F(p)\le1\) gives
\[
 X^{-w_X(p)}=p^{-1+\delta_F(p)}
             \le p^{-1}+\delta_F(p).
\]
Split the sum of \(\delta_F(p)\) at \(p=\sqrt{\log X}\) and use (RL1)
on the larger primes. This proves
\[
 \sum_{p\le X}X^{-w_X(p)}
 \ll \log\log X+\sqrt{\log X}
                 +\frac{2g(X)}{\log\log X}=o(\log X).
\]
The effective radical price of the row is \(\log_X R_F\). Apply the
uniform effective-radical bound (EB4). This proof charges the actual
fractional residue-degree weights, including primes where only some
places ramify. [Effective radical count (EB3)–(EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

For instance, if all three ideals are \(M\)-full in \(F\), then
\(R_F\le(abc)^{1/M}\le X^{3/M}\). With \(M\ge K_0/\lambda\) and
\(K_0>3\), each field accounts for at most
\(X^{(1/K_0+\eta)\lambda+\varepsilon}\) rows. Covering a family of
exponent \((1/3-\rho)\lambda\) by such fields requires at least
\(\gg X^{(1/3-1/K_0-\rho-\eta)\lambda-\varepsilon}\) distinct fields
under the same root-discriminant bound. The explicit field in (RL3)
has \(\log\operatorname{rd}(K)\le\lambda\log X+O_d(1)\) on this shell,
so such a dictionary can always be formed with variable fields.

## One lifted target in the minimal-repair count

Fix \(K_0>0\), \(0<\xi<1/3\), and \(M=\lceil K_0/\lambda\rceil\).
Use the central sector
\[
 c^{(1/3-\xi)\lambda}\le\operatorname{rad}(c)
   \le c^{(1/3+\xi)\lambda},\qquad
 \operatorname{rad}(ab)\le c^{(2/3+\xi)\lambda}.
\]
Let \(F\) be prescribed with \(\log\operatorname{rd}(F)\le g(X)\),
where \(g\) is as in (RL4). Among the original primitive triples in
\(X/2<c\le X\), those for which \((c)\) is \(M\)-full in \(F\) number
\[
 \ll_{g,\xi,\lambda,\eta,\varepsilon}
 X^{\,1/M+(1/6+\xi/4+\eta)\lambda+\varepsilon}.
 \tag{RL5}
\]
Here \(\eta>0\) is arbitrary and \(\lambda\) is sufficiently small
and fixed. Define \(\mathcal P_F\) to be the rational primes all of
whose places in \(F\) have \(e_{\mathfrak p}\ge2\). Its cardinality
is \(o(\log X)\) uniformly, since \(\delta_F(p)\ge1/2\) there and the
sum bound in the proof of (RL4) applies with the same cutoff, also for
primes above \(X\). Every prime of \(c\) outside this pool already has
\(v_p(c)\ge M\). Thus there are at most \(X^{1/M}\) possible outside
radicals. On each such support, united with \(\mathcal P_F\), the
number of complete values \(c\le X\) is \(X^{o(1)}\) uniformly. In
fact the union has \(o(\log X)\) primes and the nonnegative exponent
vectors have coordinate sum at most \(\log_2 X\). At fixed \(c\), one
of the two summands has radical at most \(X^{(1/3+\xi/2)\lambda}\).
The fixed-letter bound (FL1) gives the remaining factor in (RL5).
[Fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

For the minimally repaired height \(H_M(c)=\prod_{p\mid c}p^{\max(M,v_p(c))}\le Y\),
the bound becomes
\[
 \ll_{g,K_0,\xi,\lambda,\eta,\varepsilon}
 Y^{\,\frac{1/6+\xi/4+\eta+1/K_0}{1/3-\xi}\,\frac1M+\varepsilon}
 \tag{RL6}
\]
for each prescribed field with \(\log\operatorname{rd}(F)\le g(Y)\).
Indeed \(c\le Y^{1/[M(1/3-\xi)\lambda]}\) by the lower target radical
cap. Write \(\lambda'=(1+2\xi)\lambda\). For
\(c\le Y^{1/(4M\lambda')}\), elementary enumeration of the whole
radical and the exponent vectors gives at most
\(Y^{1/(4M)+\varepsilon}\) rows. For larger heights
\(\log Y\ll_{M,\lambda,\xi}\log c\), so the field bound at scale \(Y\)
satisfies the uniform \(o(\log X\log\log X)\) bound on each
original height shell with endpoint \(X\). Sum (RL5) over these shells
and use \(1/M\le\lambda/K_0\). The exponent ratio in (RL6) tends to
\(1/2\) as \(K_0\to\infty\) and \(\eta,\xi\to0\). A dictionary covering a
critical-size minimally repaired family must therefore have nearly
\(Y^{1/(2M)}\) fields in that parameter range. Each input is counted
by its primitive quotient and determines its unique least repair.
[Least-target-repair equivalence](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)

## Families with a common ramification bound

A prescribed family \(\mathcal F_X\) of fields can be handled as one
record if its pointwise envelope
\[
 \delta_*(p)=\sup_{F\in\mathcal F_X}\delta_F(p)
 \qquad\text{satisfies}\qquad
 \sum_p\delta_*(p)\log p\le g(X).
 \tag{RL7}
\]
For (RL4), use the single price \((1-\delta_*(p))\log_X p\).
Its effective budget bounds from below the normalized ideal radical
budget in each member field. The same convexity estimate proves
(EB3), so (RL4) counts rows that qualify in at least one member of
the family with no extra factor. For (RL5), use the common pool
\(\{p:\delta_*(p)\ge1/2\}\). Outside this pool, no member field has
ramification index at least two at every place over the prime.
Thus the same target radical enumeration proves (RL5), and then
(RL6) with the corresponding envelope bound at scale \(Y\).

For instance, let a prescribed algebraic extension, possibly
infinite, have \(\log\operatorname{rd}(F)\le g(X)\) at every finite
subfield \(F\). Those subfields satisfy (RL7). For finite fields
\(F\subseteq E\), summing the relative residue degrees over a place
of \(F\) gives at most \([E:F]\), so \(s_E(p)\le s_F(p)\) at each
rational prime \(p\). To test (RL7) on any finite prime set, choose
a subfield approaching each supremum and take their compositum
inside the extension. Its weights dominate the chosen values,
and (RL1) gives the bound. Taking the supremum over finite prime
sets proves (RL7).

## A common field on an actual shared-value box

Use an ordered box of size \(m\) from (LB1). Its \(3m^2\) labels are
positive differences between the three groups of potentials. Fix a
spanning tree \(\mathcal T\) on these \(3m\) indexed vertices. Prescribe
a height cap \(X^{h_e}\) and an integer order \(M_e\ge1\) for each
tree-edge label \(L_e\). Suppose there exists a field \(F\), chosen
separately for each box, such that
\(\operatorname{rd}(F)\le n^u\) and the ideal \((L_e)\) is \(M_e\)-full
for every \(e\in\mathcal T\), where \(n=X^L\). Then the number of
such boxes is
\[
 \ll_{m,\varepsilon}
 n^{\,u+\sum_{e\in\mathcal T}h_e/(LM_e)}X^\varepsilon.
 \tag{RL8}
\]
The bound allows every degree and counts all admissible fields
together. The tree can be chosen separately for each box at bounded
combinatorial cost.

The finite-place proof of (RL1) applies to any positive integer,
in particular the product \(D_{\mathcal T}=\prod_{e\in\mathcal T}L_e\).
The normalized ideal radical is submultiplicative, and each tree
label satisfies \(R_F(L_e)\le L_e^{1/M_e}\), where \(R_F(D)\) denotes
the normalized norm of the radical ideal of \((D)\). Hence
\[
 \operatorname{rad}(D_{\mathcal T})
    \le\operatorname{rd}(F)R_F(D_{\mathcal T})
    \le n^{\,u+\sum_{e\in\mathcal T}h_e/(LM_e)}.
\]
Enumerate this single squarefree support at the displayed cost.
Its prime product is also at most \(X^{3m-1}\) since the original
labels are at most \(X\). At fixed support each tree label has
\(X^{o(1)}\) possibilities by uniform supported-value enumeration.
The exact tree lengths and the signs given by the three roles fix
all potentials up to translation. They determine all box entries
and their residue labels. This proves (RL8).
[Actual shared-value boxes (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem)
· [Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For a surviving simultaneous packet tuple write \(d=s_1+s_2+s_3\)
and \(t=1/6-d>0\). The forced ordered-box count in (LB4) has exponent
\(md+(3m-1)t\) up to the stated fixed tolerance losses. Suppose all
original \(3m^2\) label ideals must be \(M\)-full over a common field,
with \(M\ge K_0/L\) for a fixed \(K_0>0\). For every fixed margin
\(\zeta>0\), the boxes admitting such a field with
\[
 \log_n\operatorname{rd}(F)
       \le md+(3m-1)t-(3m-1)/K_0-\zeta
\]
are power fewer than the forced boxes when the right side is
nonnegative. Use \(h_e=1\) on any spanning tree, and take the packet
tolerances and the inner allowance \(\varepsilon/L\) small in terms
of \(m,\zeta\). The orders and all packet parameters are fixed at
the budget before taking the height limit.

For the whole box, there is an explicit field attaining the
pointwise root-discriminant minimum up to the fixed-order factor
from (RL3). At each rational prime dividing a label, set \(m_p\) to
the maximum of \(\lceil M_e/v_p(L_e)\rceil\) over the labels it
divides, now allowing prescribed orders on every box label.
The local necessity \(e_{\mathfrak p}\ge m_p\) and the construction
of (RL2)–(RL3) apply verbatim. This extends the local conclusion
to repeated and noncoprime labels in distinct rows.
[Forced count (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

## Scope

Browning–Verzobio count primitive powerful triples over \(\mathbb Q\),
and their lopsided-box theorem counts primitive integer points on a
weighted generalized Fermat surface. A ramified lift supplies the
ideal fullness with no shared ideal factor, but its field becomes
part of the arithmetic data. Equations (RL4) and (RL6) show the cost
of any argument that covers a near-full family by a small list of
such fields and gains from normalized ideal radicals separately on
each field. An unrestricted improvement would need, for example, an
aggregate estimate over the varying fields using all the original
radical constraints. The counting estimates here are consequences
of (EB4) and (FL1), and do not bound that aggregate.
[Browning–Verzobio, Theorems 1.1, 1.2 and 4.1](https://arxiv.org/html/2608.24512v1)
