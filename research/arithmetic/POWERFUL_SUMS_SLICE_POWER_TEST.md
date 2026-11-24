# Divisor degrees in the powerful-sums surface estimate

The coefficient-uniform rectangular estimate in Theorem 1.2 of Browning–Verzobio, comes from a determinant cover followed by a separate analysis of its plane slices. A sharper estimate for geometrically reducible slices propagates through that proof to the full surface bound. The change concerns the stated proof of arXiv:2608.24512v1. [Retained primary source](sources/browning-verzobio-powerful-sums-2608.24512v1.pdf) · [Lemma 2.11 and Remark 3.1](https://arxiv.org/html/2608.24512v1#S3)

## Exact power test on a plane slice

Let \(C\) be the smooth projective normalization of an absolutely irreducible plane curve \(c(x,y)=0\) of degree \(d\) over \(\overline{\mathbb Q}\), with nonconstant coordinate functions \(x,y\). Let \(p,q\ge2\), \(a_1a_2a_3\ne0\), and put
\[
 A=-\frac{a_1x^p+a_2y^q}{a_3},\qquad U=\frac{a_1x^p}{a_2y^q},\qquad H=\deg(U:C\longrightarrow\mathbb P^1).
\]
Assume \(U\) is nonconstant. If \(A=h^k\) in \(\overline{\mathbb Q}(C)\) for \(k\ge2\), then
\[
 (1-1/k)H\ \le\ 2g(C)-2+
 \bigl|\operatorname{supp}(\operatorname{div}x)\cup
       \operatorname{supp}(\operatorname{div}y)\bigr|\ \le\ d^2.
 \tag{1}
\]
Indeed, write \(S_0\) for the union of those supports. All points above \(0\) and \(\infty\) for \(U\) are in \(S_0\). At a point outside \(S_0\) above \(-1\), the identity \(1+U=-a_3h^k/(a_2y^q)\) shows that the multiplicity in the \(-1\) fibre is divisible by \(k\). The sum of multiplicities over that fibre is \(H\), so altogether at most \(|S_0|+H/k\) points lie above the three values. Riemann–Hurwitz gives total ramification \(2g(C)-2+2H\), with at least \(3H-|S_0|-H/k\) accounted for by the three fibres. Finally the zeros of \(x,y\) contribute at most \(2d\) points, their poles lie above the line at infinity and contribute at most \(d\) points, and \(g(C)\le(d-1)(d-2)/2\).

The individual bound \(|\operatorname{ord}_P(U)|\le H\) follows because its zero and pole divisors **each** have degree \(H\). The proof step in Lemma 2.11 uses \(2H\), the degree of their union counted with multiplicities. If the two supports are incomparable, there is a point contributing solely an order of \(x\), and another solely an order of \(y\). Hence \(H\ge\max(p,q)\) and
\[
 \boxed{\ \max(p,q)\le \frac{k}{k-1}d^2\ } . \tag{2}
\]
For individual slices the first inequality of (1), using the actual function degree and support size, can exclude additional powers. In the incomparable case \(A\ne0\) automatically. Over \(\overline{\mathbb Q}(C)\), Capelli's criterion makes the degree test applicable to reducibility of \(T^r-A\) by testing primes dividing \(r\). The exceptional Capelli case for \(4\mid r\) is already a square over this field. [Primary proof and Capelli statement](https://arxiv.org/html/2608.24512v1#S2.SS3)

## Coefficient-uniform rectangular bound

Let \(p,q,r\ge2\) be fixed integers, \(a_1,a_2,a_3\) nonzero integers, and \(X,Y,Z\ge2\). Put \(B=\max(X,Y,Z)\), \(M=\max(X,Y)\), \(P=\max(p,q)\) and \(W=\exp(\sqrt{\log X\log Y/r})\). The number \(S\) of primitive signed integer triples in this rectangle satisfying \(a_1x^p+a_2y^q+a_3z^r=0\) satisfies
\[
 \boxed{\ S\ll_{\varepsilon,p,q,r} (XYZ)^\varepsilon
 \left(W^2+W M^{\min\{1/3,\sqrt{r/((r-1)P)}\}}+W M^{1/r}\right).\ } \tag{3}
\]
The constant is independent of \(a_1,a_2,a_3\), and all three arrangements of the distinguished variable are available. Equation (3) replaces the term \(M^{2/\sqrt{\max(p,q,36)}}\) in Theorem 1.2 and the term \(M^{\min\{1/3,\sqrt{2r/((r-1)P)}\}}\) in Remark 3.1.

The determinant cover in Section 3 supplies \(O_{\varepsilon,p,q,r}(B^\varepsilon W)\) auxiliary equations and \(O_{\varepsilon,p,q,r}(B^\varepsilon W^2)\) exceptional points. Eliminating \(z\) by a resultant gives nonzero plane equations with total factor degree \(O_{\varepsilon,p,q,r}(\log B)\) per auxiliary equation. Keep the original treatment of degrees \(1,2\), degrees at least \(\max(p,q,r)\), and the special slices of Lemmas 2.8–2.10. The original bounds there depend on the degrees and box, uniformly in the coefficients. [Determinant cover and resultant](https://arxiv.org/html/2608.24512v1#S3) · [Side-condition theorem, Theorem 1.1](https://arxiv.org/abs/2408.11453)

For an intervening absolutely irreducible slice, take the largest \(k\mid r\) such that \(A=h^k\), allowing \(k=1\). Over \(\overline{\mathbb Q}(C)\) the \(k\) factors \(T^{r/k}-\zeta h\), \(\zeta^k=1\), are irreducible. Otherwise Capelli would make \(\zeta h\) an \(\ell\)-th power for a prime \(\ell\mid r/k\), which would make \(A\) a \((k\ell)\)-th power since constants have all roots. The primitive-element projection of Lemma 2.17 uses \(w=x+\lambda y\), with \(\lambda\) chosen from a bounded set and \([\overline{\mathbb Q}(C):\overline{\mathbb Q}(w)]=d\), so each projected irreducible component has degree \(dr/k\) in \(z\). Its degree and height are bounded by the same resultant after the coefficient-height reduction in Lemma 2.8. A component defined only over an extension contributes \(O_{p,q,r}(1)\) rational points by intersection with a distinct conjugate. Heath-Brown's rectangular curve estimate on the other components gives \(O_{\varepsilon,p,q,r}(B^\varepsilon M^{k/(dr)})\) points. Multiplicity above a projected \((w,z)\) is bounded by \(d\). [Lemmas 2.8 and 2.17 and Remark 3.1](https://arxiv.org/html/2608.24512v1) · [Heath-Brown, Theorem 15](https://arxiv.org/abs/math/0405392)

For \(k=1\) this fits \(B^\varepsilon M^{1/r}\). For \(k\ge2\), (2) gives
\[
 \left(\frac{k}{dr}\right)^2\le \frac{k^3}{(k-1)r^2P}\le\frac{r}{(r-1)P}.
\]
The second inequality uses monotonicity of \(k^3/(k-1)\) for \(k\ge2\). The uniform plane-curve estimate gives also \(B^\varepsilon M^{1/d}\) by counting plane points first and then the bounded number of \(z\) values. In the intervening range \(d\ge3\) this caps the exponent at \(1/3\). Summing the slices absorbs their logarithmic count into \(B^\varepsilon\), proving (3).

## Transfer to the fixed-offset surface step

In the [fixed-offset manuscript, Section 6.1](../../unpublished/powerful-number-representations/main.tex), fix all but a maximum-scale coordinate on each side of \(a+b=c\). The remaining exponents are all at least the base exponent \(s\), the coordinate bounds are \(O_s(L^{1/s})\), and the three coordinates are primitive. The coefficient bound used there, \(L^{3s^{-3/2}+\varepsilon}\), can be replaced by
\[
 L^{s^{-3/2}(1+\sqrt{s/(s-1)})+\varepsilon}.
 \tag{4}
\]
In (3) the distinguished exponent is at least \(s\), so \(W\ll_s L^{s^{-3/2}}\), \(\sqrt{r/((r-1)P)}\le 1/\sqrt{s-1}\), and \(1/r\le1/s\). Each of the three terms gives (4). If the three normalized maximum masses are \(\alpha,\beta,\chi\), the exponent inequality in the manuscript accordingly improves to
\[
 \Delta\le S_A+S_B+S_C-S_A\alpha-S_B\beta-S_C\chi
             +s^{-3/2}(1+\sqrt{s/(s-1)})
\]
with the usual arbitrarily small dyadic exponent loss. The leading surface loss changes from \(3s^{-3/2}\) to \((2+O(s^{-1}))s^{-3/2}\). After normalizing, the maximum-mass slack is still \(O_{u,v}(s^{-1/2})\). Thus this proof step strengthens \(\texttt{ARITH-POWERFUL-FIXED-OFFSET-SURFACE-MOMENTS}\). For the former \(\texttt{GAP-ARITH-POWERFUL-FIXED-OFFSET-DOUBLE}\), the leading grouped-constraint obstruction persists on this surface route. The manuscript reaches the leading saving by its separate monomial argument. The determinant cover applies to fixed-coefficient surfaces. Summing it over variable coefficient tuples would require a separate saving, as in the pointwise square-full kernel-pair problem.
