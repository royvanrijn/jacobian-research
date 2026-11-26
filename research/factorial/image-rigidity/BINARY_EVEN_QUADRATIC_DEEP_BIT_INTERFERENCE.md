# Exact-bit dyadic separation for even binary quadratics

Normalize the constant coefficient to one and write
\[
 f=ax^2+bxy+cy^2+1,\qquad M_j=\Lambda_2(f^j),\qquad
 A=v_2(4a),\quad C=v_2(4c),\quad Q=v_2(4b).
\]
Work over a valued characteristic-zero field with \(v_2(2)=1\).

**Theorem.** Suppose \(a\ne c\), \(d=2^e\ge8\), and \(A,C\in(1/d,2/d)\). If \(M_n,\dots,M_{n+3}\) vanish for some \(n\ge1\), then \(A=C=r\), \(v_2(n+3)=e\), and
\[
                    \frac{2}{d}\le h:=2-dr<\frac14.             \tag{1}
\]
Let \(\zeta\) denote the residue of \(c/a\). It is a nontrivial \((d-1)\)-st root of unity. If \(h\) differs from \(2^{-j}\) for every \(0\le j<e\), define \(t\) as the number of integers \(0\le j<e\) with \(2^j h<1\). The proposed collision then requires
\[
                           \zeta^{\,2^{t-1}-1}=1.                \tag{2}
\]
It also requires one of the two strengths in (8). Consequently \(d=8\) has no four-window in this sector. When \(e\) is prime, every nonresonant \(h\) is excluded.

## The exact-bit branch

The [deeper dyadic boundary reduction](BINARY_EVEN_QUADRATIC_INTERMEDIATE_DIGIT_SECTORS.md#a-uniform-extension-to-deeper-dyadic-intervals) shows \(A=C=r\), \(v_2(N)=e\) for \(N=n+3\), and the asserted properties of \(\zeta\). For the boundary value \(g_N(y)=\Lambda_1((1+y x_1^2)^N)\), the coefficient of \(Z^i\) in \(g_N(Z/4)\) has valuation \(s(N-i)-s(N)\), with \(s\) the binary digit sum. At \(v_2(Z)=r\), every nonconstant term has positive valuation and the unique least nonconstant term is \(i=d\), of valuation \(dr-1\). Indeed \(N/d\) is odd. For \(i=dv+w\) with \(0\le w<d\), the block lower bound is \(v(dr-1)\) if \(w=0\) and
\[
                     v(dr-1)+wr+s(d-w)-1
\]
otherwise. The latter exceeds \(dr-1\) since \(s(d-w)\ge 2(d-w)/d>(d-w)r\).

The constant part in \(g_N(a)+g_N(c)\) is two. The two terms with \(i=d\) have valuation \(dr-1<1\); their residue ratio is \(\zeta^d=\zeta\ne1\). Hence \(v_2(g_N(a)+g_N(c))=dr-1\). Both \(g_N(a)\) and \(g_N(c)\) reduce to one, so \(2(a g_N(a)+c g_N(c))\) has valuation \(r-1\). The [weighted boundary equation](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation) now gives
\[
                       b=\frac{2(a g_N(a)+c g_N(c))}
                                    {g_N(a)+g_N(c)},\qquad Q=r+h. \tag{3}
\]

For an original moment term, let \(i,\ell,m\) be its counts of \(a,b,c\), and \(k=i+\ell+m\). Define
\[
 \eta_j(k)=kr+s(j-k)-s(j),\qquad \phi_h(\ell)=h\ell-s(\ell).
\]
The exact term valuation is
\[
 \eta_j(k)+\phi_h(\ell)
  +\operatorname{car}(i,\lfloor\ell/2\rfloor)
  +\operatorname{car}(m,\lfloor\ell/2\rfloor),                    \tag{4}
\]
where \(\operatorname{car}\) counts binary carries. This is the [multinomial digit formula](BINARY_EVEN_QUADRATIC_INTERMEDIATE_DIGIT_SECTORS.md#coefficient-valuations-and-a-useful-integer-ratio). The unique minimum of \(\eta_j\) is at \(k=u=j\bmod d\), by the same deeper dyadic calculation. For fixed \(k,\ell\), the residue of the sum of its carry-free terms, after division by \(a^{k-\ell}b^\ell\) and the common power of two in their integer coefficients, is
\[
                       \sum_{\text{carry-free }m}\zeta^m.         \tag{5}
\]
Indeed the integer coefficients after this division are all odd. This normalization applies to a fixed pair \((k,\ell)\), even if some other pair has the same minimum valuation.

## Elimination for \(h\ge1/4\)

Each bit \(2^q\) contributes \(2^q h-1\) to \(\phi_h\) when present and zero when absent. This identifies its global minima.

- If \(h>1/2\), the sole minimizing \(\ell\) is one. Use \(j=N-3\), so \(k=u=d-3\) at the minimum of \(\eta_j\). Every split \(i+m=k-1\) is carry-free, and (5) equals \((1+\zeta^{d-3})/(1+\zeta)\ne0\). Its numerator could vanish only if \(\zeta^2=1\).
- If \(h=1/2\), the minimizing counts are one and three. Use \(j=N-2\), so \(u=d-2\) is even. At \((k,\ell)=(u,3)\) a carry-free split would require both \(i,m\) even, though their sum is odd. At \((u,1)\) the residue equals \((1+\zeta^{d-2})/(1+\zeta)\ne0\).
- If \(1/4<h<1/2\), the sole minimizing count is three. Use \(j=N-1\) and \(u=d-1\). Both axis counts must be even. Their residue sum is \((1+\zeta^{d-2})/(1+\zeta^2)\ne0\).
- If \(h=1/4\), the minimizing counts are three and seven. At \(j=N-3\) and \(u=d-3\), count seven cannot occur carry-free since both axis counts would be multiples of four and sum to \(d-10\) (if seven fits at all). At count three the carry-free residue is \((1+\zeta^{d-4})/(1+\zeta^2)\), which is nonzero unless \(\zeta^3=1\).

Consider the remaining endpoint \(h=1/4,\ \zeta^3=1\), and take \(j=N-2\), with \(u=d-2\). For each odd \(k\), \(j-k\) is odd, so
\[
                        \eta_j(k)-\eta_j(k+1)=1-r>1/2.
\]
Since \(\min\phi_h=-5/4\), every odd \(k\) lies strictly above \(\eta_j(u)-3/4\). For even \(k\), counts three and seven each cost at least one carry. Of the remaining counts, one and five alone attain \(\phi_h=-3/4\), and all other counts are higher. Thus the only least pairs are \((u,1)\), \((u,5)\), each with carry-free splits. For \((u,5)\) the carry-free \(m\) are \(4v,4v+1\) for \(0\le v\le(d-8)/4\), since both axis counts must avoid the bit of value two and sum to \(d-7\). The residue is
\[
                       (1+\zeta)\frac{1+\zeta^{d-4}}{1+\zeta^4}=0.
\]
For \((u,1)\) it is \((1+\zeta^{d-2})/(1+\zeta)\ne0\). This proves a nonzero moment at the endpoint.

## Elimination using the final odd moment

Use \(j=N-1\), so \(u=d-1\), and write \(k=dv+w\), \(\ell=dv_\ell+l\), with \(0\le w,l<d\), \(v_\ell\le v\). There is no borrow in the low block when subtracting \(k\) from \(j\). Its digit loss is \(-s(w)\), and the high block loses at most \(s(v)\) by digit subadditivity. With \(\mu=\min_{0\le l<d}\phi_h(l)\), the term valuation minus \(\eta_j(u)+\mu\) is therefore at least
\[
 v(dr-1)+v_\ell(dh-1)
 +s(d-1-w)-(d-1-w)r
 +\phi_h(l)-\mu.                                                \tag{6}
\]
The first two terms are nonnegative and strictly positive if \(v>0\). For \(dh\ge1\) this uses \(dr>1\). For \(dh<1\), use \(v_\ell\le v\) and \(d(r+h)-2=(d-1)h>0\). The next two terms together are strictly positive unless \(w=d-1\), since each low bit has size at most \(d/2\). Consequently only \(k=d-1\) and the low minimizers of \(\phi_h\) can attain \(\eta_j(u)+\mu\); they do so exactly when their carry counts vanish.

If \(h<2/d\), every low bit is profitable. The only minimizing count is \(d-1\), which gives the pure cross term and proves \(M_{N-1}\ne0\).

For the nonresonant claim, in view of the already excluded ranges one has \(t\ge2\). The sole minimizing count is \(\ell=2^t-1\). Put \(P=2^{t-1}\). Both axis counts are multiples of \(P\), since \(\lfloor\ell/2\rfloor=P-1\). They sum to \(d-2^t\). Their residue sum is
\[
                   \frac{1+\zeta^{d-P}}{1+\zeta^P}
                    =\frac{1+\zeta^{1-P}}{1+\zeta^P}.           \tag{7}
\]
The denominator is nonzero by the odd order of \(\zeta\). This yields (2). If \(e\) is prime, \(1\le t-1<e\) implies
\[
              \gcd(2^e-1,2^{t-1}-1)=2^{\gcd(e,t-1)}-1=1,
\]
so (2) is impossible.

Equation (6) retains two possible minimizing cross counts at reciprocal powers of two within (1). For composite \(e\) and nonresonant \(h\), (7) can vanish in a proper residue subfield. These are necessary conditions for the remaining exact-bit cases in the common annulus.

## The next valuation layer in a proper residue subfield

Suppose \(h\) is nonresonant, (1) holds, and the leading group (7) vanishes. Thus \(3\le t\le e-1\) and \(\zeta^P=\zeta\) for \(P=2^{t-1}\). In fact only the two strengths
\[
                     h=\frac{2}{2P+1}
                     \quad\text{or}\quad h=\frac{2}{4P-1}       \tag{8}
\]
can survive in this case.

To see this, put \(Z=4a,\ X=c/a\), and let \(m_0\) be the odd order of \(\zeta\). In \(g_{N-1}(Z/4)\) the unique least term has index \(u=d-1\), and the margin of every other term is at least \(h/2\). Indeed, for the terms of index below \(d\), omitting the bit of value \(d/2\) raises \(\eta_{N-1}\) by \(1-dr/2=h/2\), and omitting any other bit costs more. For indices at least \(d\), the high block alone costs at least \(dr-1>h/2\). In \(g_N(Z/4)\) every nonconstant term has valuation at least \(dr-1>h/2\). Therefore, after division by its least monomial, the boundary quotient is one plus a correction of valuation at least \(h/2\), on the full annulus of axis valuation \(r\). Comparing the quotients at the two axes shows \(v_2(X^{1-d}-1)\ge h/2\). Since \(m_0\mid d-1\) and \((d-1)/m_0\) is odd, factoring a difference of powers at \(X^{m_0}\) gives
\[
                              v_2(X^{m_0}-1)\ge h/2 .
\]

At the leading pair \((u,2^t-1)\), the sum of the carry-free coefficients after the normalization in (5) agrees with \((X^{d-P}-1)/(X^P-1)\) up to valuation at least one. The denominator is a unit, and \(2m_0\mid d-P\) because \(m_0\mid d-P\) is odd and \(d-P\) is even. Squaring \(X^{m_0}=1+E\) with \(v_2(E)\ge h/2\) raises the valuation of its difference from one strictly above \(h/2\). Raising to any further positive integer power keeps that lower bound. The whole contribution at this pair therefore lies strictly above the level \(\eta_{N-1}(u)+\mu+h/2\); terms with carries add at least one.

It remains to find the first surviving pair below or at that level. In (6), high blocks cost at least \(dr-1>h/2\) since \(dh\ge2\). Among low blocks the unique smallest positive axis cost is \(h/2\), at \(k=d/2-1\). The only cross modifications whose cost can be at most \(h/2\) are dropping \(P\), of cost \(1-Ph\), or adding \(2P\), of cost \(2Ph-1\). Every other single bit flip costs at least \(1/2\), and the two displayed costs cannot both be at most \(h/2\) since their combination twice the first plus the second is one and \(3h/2<1\). All bit costs are positive, so no simultaneous modification involving a positive axis cost can reach this level.

At \(k=d/2-1\) and the original cross count \(\ell=2^t-1\), a carry-free split again uses axis counts that are multiples of \(P\), now summing to \(d/2-2P\). Its residue numerator is \(1+\zeta^{d/2-P}\), which is nonzero because \(\zeta^{d/2}\) is the square root of \(\zeta\) and \(\zeta^P=\zeta\). At \(k=u\) either modified cross count has the form (7) with \(P\) replaced by \(P/2\) or \(2P\). Their residue numerators are nonzero since \(\zeta^{P/2}\) and \(\zeta^{2P}\) both differ from \(\zeta\). It follows that the first surviving residue is nonzero unless an axis cost and a cross cost coincide at \(h/2\). Solving these two equalities gives (8).

## Exact checks

The [finite replay](verify_binary_even_quadratic_deep_bit_interference.py) checks 1,000 coefficient valuations, 45 exact-bit axis minima, 99 least-pair fixtures, 132 residue identities and order checks, and 40 next-layer checks. The proof supplies the uniform bounds for larger degrees.
