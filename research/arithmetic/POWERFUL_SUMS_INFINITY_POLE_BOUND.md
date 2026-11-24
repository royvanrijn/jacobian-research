# Poles at infinity in the powerful-sums surface estimate

Browning–Verzobio's determinant-cover proof counts each Kummer component by its degree over a generic plane coordinate. The degree of the proper power in the cover is restricted by the poles at infinity of the two coordinate monomials. For nearby distinct exponents this restricts the exceptional-slice loss to order \(1/r\) in the box-side exponent. [Primary proof, Lemma 2.17 and Remark 3.1](https://arxiv.org/html/2608.24512v1#S3) · [Previous divisor-height estimate](POWERFUL_SUMS_SLICE_POWER_TEST.md)

## The pole obstruction

Let \(C\) be the smooth projective normalization of a geometrically integral plane curve over \(\overline{\mathbb Q}\), of degree \(d\), with nonconstant coordinates \(x,y\). Let \(p,q,r\ge2\), and let \(a_1,a_2,a_3\) be nonzero. Suppose that
\[
 A=-(a_1x^p+a_2y^q)/a_3=h^k\ne0,\qquad k\mid r,\quad k\ge2
\]
in the function field. Put \(P=\max(p,q)\). Then
\[
 \frac{k}{dr}\le
 \Lambda(p,q;r):=
 \max\left\{\frac{\gcd(p,q)}P,\frac{\gcd(p,r)}r,
                              \frac{\gcd(q,r)}r\right\}.                 \tag{1}
\]
This holds regardless of whether the supports of the two coordinate divisors are comparable.

At a point \(Q\) over the line at infinity put \(e_x=\max(0,-\operatorname{ord}_Q(x))\) and define \(e_y\) similarly. At least one is positive, and \(e_x,e_y\le d\). Indeed the pullback of the line at infinity has degree \(d\), and the order of its local defining section is \(\max(e_x,e_y)\), since one of the two projective coordinate numerators is a local unit. If \(p e_x=q e_y\) at any point at infinity, both are positive, so \((e_x,e_y)\) is a positive integer multiple of \((q,p)/\gcd(p,q)\). In particular \(d\ge P/\gcd(p,q)\), which proves the first alternative in (1). Otherwise pick a point at infinity. The larger pole order, say \(p e_x\), is the pole order of \(A\). It is divisible by \(k\), so \(k/\gcd(k,p)\) divides \(e_x\). This gives \(k/d\le\gcd(k,p)\le\gcd(r,p)\). The case of larger \(q e_y\) yields the third term of (1).

## Equal plane exponents

Suppose \(p=q\), with \(u=x/y\) nonconstant. The map \(u\) has degree at most \(d\), by the pencil of lines through the origin. Now \(A\) is a nonzero constant times \(y^p(u^p+a_2/a_1)\). The \(p\) roots of the second factor are distinct nonzero constants. If \(k\mid p\), at every point over any of these roots the local multiplicity in the \(u\)-fibre is divisible by \(k\), since the order of \(A\) and the contribution from \(y^p\) are divisible by \(k\). Hence \(k\mid\deg u\), and \(k/(dr)\le1/r\).

Even if \(k\nmid p\), at most \(2d\) of these root values can meet the support of the divisor of \(y\), since each of its zero and pole divisors has degree at most \(d\). If \(p>2d\), there is a fibre disjoint from that support and the same degree divisibility follows. If \(p\le2d\), then \(k/(dr)\le1/d\le2/p\). In particular, when the two plane exponents agree the uniform alternative for (1) is
\[
 \Lambda(p,p;r)=
 \begin{cases}
    1/r,&r\mid p,\\[2pt]
    \max(1/r,2/p),&r\nmid p.
 \end{cases}                                                         \tag{1a}
\]
A constant \(x/y\) puts the slice on a line and is already treated by the small-degree estimate.

## Rectangular and nearly equal exponent consequences

Keep the signed primitive point count, the boxes \(X,Y,Z\ge2\), \(M=\max(X,Y)\) and \(W=\exp(\sqrt{\log X\log Y/r})\) of the [previous surface proof](POWERFUL_SUMS_SLICE_POWER_TEST.md#coefficient-uniform-rectangular-bound). The count is
\[
 \ll_{\varepsilon,p,q,r}(XYZ)^\varepsilon
 \left(W^2+W M^{1/r}
   +W M^{\min\{1/3,\,\sqrt{r/((r-1)P)},\,\Lambda(p,q;r)\}}\right).         \tag{2}
\]
To see this, keep the determinant cover, resultant elimination, and special-slice estimates in the previous proof. For an intermediate slice with incomparable supports, \(A\neq0\). Take the maximal \(k\mid r\) for which \(A=h^k\), allowing \(k=1\), so each of the \(k\) components has function field degree \(r/k\) over \(\overline{\mathbb Q}(C)\). For \(z\) nonconstant on that component, the generic projection in Lemma 2.17 and the lopsided curve count give \(O_{\varepsilon,p,q,r}(B^\varepsilon M^{k/(dr)})\), where \(B=\max(X,Y,Z)\). A constant \(z\) lies on a component with \(k=r\), and the plane-curve bound gives the same exponent \(1/d=k/(dr)\) directly. Apply the new bounds (1) and (1a) and the earlier divisor-height bound to \(k\ge2\), and count the \(k=1\) or otherwise pretreated slices by the original estimates. The uniform plane-curve estimate caps the intermediate-slice exponent at \(1/3\). Summing over the same covering gives (2), uniformly in the nonzero integer coefficients. [Primary proof of the cover and projection](https://arxiv.org/html/2608.24512v1#S3)

If \(p,q,r\) are pairwise coprime and \(r\le P\), then \(\Lambda=1/r\). Thus only \(W^2+W M^{1/r}\) is needed. In particular the bound is \(O_{\varepsilon,p,q,r}((XYZ)^\varepsilon W^2)\) when \(\log\min(X,Y)\ge\log M/r\). For three distinct exponents in \([s,s+A]\), where \(A\ge1\) is fixed, every pairwise gcd is at most \(A\) because it divides the nonzero difference. Therefore \(\Lambda\le A/s\) for any choice of distinguished exponent. The same \(W^2\) conclusion holds when \(\log\min(X,Y)\ge A^2r\log M/s^2\). More generally the \(W M^{A/s}\) term suffices for the reducible slices even in strongly lopsided boxes. If two exponents agree, choose them as the plane exponents and use (1a) to get exponent \(O(1/s)\). When all three agree, any distinguished variable gives \(W^2+W M^{1/r}\).

## Retaining the actual side lengths in fixed-offset boxes

Consider the [competitive fixed-offset boxes](../../unpublished/powerful-number-representations/main.tex) for orders \((s+u,s+v,s)\), with fixed \(u\ge v\ge0\). In the canonical factorization of a \(b\)-full number, the coordinates have exponents \(b,b+1,\ldots,2b-1\). A prime with exponent \(qb+j\), \(0\le j<b\), goes to the first coordinate with exponent \(q\) if \(j=0\). Otherwise it goes to the first coordinate with exponent \(q-1\) and to the \(b+j\) coordinate once. Thus the factorization is unique. Let \(S_i\) be each side's total logarithmic coordinate mass in base \(X\). The size constraints give \(S_i\le1/s\). Fixing each pair of complete canonical tuples determines the third, so the box counting exponent \(\Delta\) is at most the least of the three pair sums. In a competitive box with \(\Delta\ge2/s-O_{u,v}(s^{-2})\) this yields \(S_i=1/s+O_{u,v}(s^{-2})\). The size constraint on each side further says \(S_i(s+m_i)\le1\), where \(m_i\) is the offset-plus-index mean after normalizing by \(S_i\). In particular all \(m_i\) are bounded. Let \(\alpha\le\beta\le\gamma\) be the three maximum normalized masses on their respective sides and \(T=\alpha+\beta+\gamma\). Their maxima have bounded indices \(O_{u,v}(1)\). Indeed each normalized vector has bounded first moment, so at least half of its mass lies in a bounded initial interval. Its maximum is bounded below by a positive \(u,v\)-dependent constant. The index carrying that mass is bounded by its first moment divided by the mass.

Fix the remaining coordinates. Their contribution is at most \(X^{\sum_i S_i(1-\alpha_i)}\), where \(\alpha_i\) is the free maximum on that side. The free coordinates satisfy a primitive diagonal equation with fixed nonzero integer coefficients because each free coordinate divides its own primitive summand. Use the dimension indexed by \(\gamma\) as the distinguished one in (2). The logarithms of the other two side lengths to the base of the target bound are \(\alpha/s+O_{u,v}(s^{-2})\), \(\beta/s+O_{u,v}(s^{-2})\). The divisor-height version of (2) then gives an exponent loss at most
\[
 (\sqrt{\alpha\beta}+\beta)s^{-3/2}+O_{u,v}(s^{-2}).
\]
The \(W^2\) term fits this estimate since \(\sqrt{\alpha\beta}\le\beta\). The maximum of \((\sqrt{\alpha\beta}+\beta)/T\) on \(\alpha\le\beta\le\gamma\) is \((1+\sqrt3)/4\). To verify it, set \(\alpha=t^2\beta\), \(0\le t\le1\), and use \(T\ge(t^2+2)\beta\). The function \((1+t)/(2+t^2)\) has its maximum at \(t=\sqrt3-1\), with equality also requiring \(\beta=\gamma\). Combining the surface count with the fixed coordinates gives
\[
 s\Delta\le 3-T+\frac{1+\sqrt3}{4\sqrt{s}}T+O_{u,v}(s^{-1}),
 \qquad
 T\le1+\frac{1+\sqrt3}{4\sqrt{s}}+O_{u,v}(s^{-1}).                        \tag{3}
\]
Here and below an arbitrarily small dyadic exponent loss is absorbed by first taking it small at fixed \(s,u,v\).

When the three selected exponents are distinct, the bounded-index observation puts them in \([s,s+A]\) for a fixed \(A=A(u,v)\). This time the \(W M^{A/s}\) term costs at most \(\sqrt{\alpha\beta}\,s^{-3/2}+O_{u,v}(s^{-2})\). The \(W^2\) term costs \(2\sqrt{\alpha\beta}\,s^{-3/2}+O_{u,v}(s^{-2})\) and hence governs this estimate. Since \(\alpha\le\beta\le\gamma\), one has \(\alpha+\beta\le2T/3\) and \(2\sqrt{\alpha\beta}\le2T/3\). Consequently,
\[
 T\le1+\frac{2}{3\sqrt{s}}+O_{u,v}(s^{-1})
 \qquad\text{when the three chosen powers are distinct}.              \tag{4}
\]
The refinement (3) applies to the surface step for all fixed offsets. Bound (4) also holds if the two smaller sides have the same selected exponent, by (1a). This includes all boxes with three equal selected exponents. The \(s^{-1/2}\) mass slack in this route still comes from the covering term \(W^2\). The short-congruence loss is a separate estimate in \(\texttt{GAP-ARITH-POWERFUL-FIXED-OFFSET-SHARP-ERROR}\).


## Explicit surface-only error for offsets \((u,0)\)

For every fixed \(u\ge0\) and sufficiently large \(s\), the primitive positive triples of orders \((s+u,s,s)\) satisfy
\[
 N_{s+u,s,s}(X)\ll_{u,s,\varepsilon}
 X^{\,1/(s+u)+1/s-2/s^2+
         \frac{3(1+\sqrt3)}4\,s^{-5/2}+O_u(s^{-3})+\varepsilon}.          \tag{5}
\]
Indeed denote the three cost means, including their offsets, by \(m_A,m_B,m_C\), and let \(D=\max(m_A+m_B,m_A+m_C,m_B+m_C)\). For a probability vector with maximum \(z\), truncating the indices at \(j\ge h\) gives mean at least \(h-\binom{h+1}{2}z\). Use \(h=2\) on \(A\) and \(h=3\) on \(B,C\). This gives
\[
 D\ge m_A+(m_B+m_C)/2\ge u+5-3T .
\]
Combine (3) with \(\Delta\le2/s-D/s^2+O_u(s^{-3})\), obtained by bounding each \(S_i\le1/(s+m_i)\) in the pair-count inequalities. Boxes below the competitive threshold \(2/s-(u+2)/s^2\) already meet (5). Within the boxes in (4), or with equal selected exponents on the smaller two sides, the coefficient of \(s^{-5/2}\) improves to \(2\). This estimates the finite-order loss of the surface route, using the same elementary truncations as the fixed-offset manuscript.
