# Moments of the odd radial perturbation

On \(2xy+t^2=1\), write \(z=x/t^2\) and set
\[
 Q(z)=\frac{(1+z^2)^4}{z},\qquad
 S(t)=t^{16}(t^2-1),\qquad g_b=S(t)Q(z)+bt^{17}.
\]
The terms of negative angle together give \(-2yt^{18}\),
so the phase is regular and has least sphere degree 19.

For \(d\) an even nonnegative integer, put
\[
 J(d,h)=\int_0^1 t^d(1-t^2)^h\,dt
       =\frac{2^h h!}{\prod_{v=0}^h(d+1+2v)}.
\]
The equality follows by integration by parts on \(h\). The
sphere functional takes the angular constant term and integrates
the even height terms over \(0\le t\le1\).

**Proposition (two unmarked orders).**
For every \(b\) over any characteristic-zero field, at least one of
\(\mathcal A(g_b^2)\) and \(\mathcal A(g_b^4)\) is nonzero.

Odd powers of \(Q\) have zero angular constant term. The constant
terms in \(Q^2,Q^4\) are \(\binom81=8\) and \(\binom{16}{2}=120\).
Writing \(M_j=\mathcal A(g_b^j)\), binomial expansion therefore
gives
\[
 M_2=\frac{b^2}{35}+\frac{64}{42735},\qquad
 M_4=\frac{b^4}{69}
       +\frac{128b^2}{109411}+\frac{3072}{103831039}.
\]
Taking the remainder of \(M_4\) modulo \(M_2\) gives the identity
\[
 M_4-\left(\frac{35}{69}b^2+\frac{5754560}{400772493}\right)M_2
       =\frac{3753620480}{464386710041397}.
\]
The right side is a nonzero rational number. The degree-19 Hopf
lift therefore also has a nonzero contraction moment at order 2
or 4 by the [sphere normalization](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md).

## The coherent marker and the correction columns

Allow
\[
 H=t^{12}(a_0+a_1z)(1+z^2)^2+t^{11}C(z),\qquad \deg C\le5.
\]

**Theorem (finite coherent detection for every radial coefficient).**
Over any characteristic-zero field, the simultaneous vanishing of
\(\mathcal A(Hg_b^k)\) at \(k=0,1,3,5,7,9\) forces \(a_0=a_1=0\).
If \(b\ne0\), adding rows \(2,4,6\) forces \(C=0\) as well. If
\(b=0\), every correction marker \(t^{11}C(z)\) is invisible in all
rows by reflection in \(t\).

Row zero is \(a_0/13\). Suppose \(a_0=0\) and \(a_1\ne0\).
If \(b=0\), row 1 is \(-2a_1/899\). For \(b\ne0\) put
\(X=b^2\) and \(D_j=b[z^{2j}]C/a_1\) for \(0\le j\le2\).
Only these three correction coefficients occur in the odd rows.
Indeed, if \(h\) factors are taken from \(S(t)Q(z)\), the
coherent constant term comes from \(h\) odd and forces a
binomial index \((h-1)/2\). The correction term has even height
only when \(h\) is even, and its angle is zero only when the
marker angle is even as well.

For \(k\) odd, divide the row by \(a_1\). It is
\(U_k+\sum_{j=0}^2V_{k,j}D_j\), where
\[
 \begin{aligned}
 U_k&=-\!\!\sum_{\substack{1\le h\le k\\ h\text{ odd}}}
   \binom{k}{h}\binom{4h+2}{(h-1)/2}
               X^{(k-h)/2}J(17k-h+12,h),\\
 V_{k,j}&=\!\!\sum_{\substack{2j\le h<k\\ h\text{ even}}}
   \binom{k}{h}\binom{4h}{h/2-j}
               X^{(k-h-1)/2}J(17k-h+11,h).
 \end{aligned}
\]
Rows 1, 3 and 5 solve for \(D_0,D_1,D_2\) in order.
Their values are
\[
 D_0=2/31,\quad
 D_1=427X/93-204/2077,\quad
 D_2=-707X^2/6+551306X/6231-360474/213931 .
\]
The two resulting residual rows at 7 and 9 have degrees
three and four. After primitive integral normalization with
negative leading coefficients, their reductions at 31 are
\[
     8X^3+27X+15,\qquad 25X^4+5X^2+20X .
\]
Their degrees are preserved and the successive Euclidean
remainders, dividing the quartic by the cubic first, are
\(2X^2+8X\) and \(15\). Thus the rational resultant is
nonzero. The residual rows have no common zero and \(a_1=0\).

With the coherent component zero, the rows \(1,\ldots,6\)
are triangular on the correction angles \(0,\ldots,5\).
The first entry at angle \(j\) requires \(j\) factors of
\(S(t)Q(z)\) and one radial odd factor. Its coefficient is
\((j+1)b(-1)^jJ(16j+28,j)\), which is nonzero if \(b\ne0\).
This recovers every correction coefficient.

The exact [checker](verify_cm_newton_odd_radial_perturbation.py)
replays the two unmarked moments and the polynomial eliminations.

```sh
sage -python research/foundations/extended-geometry/verify_cm_newton_odd_radial_perturbation.py
```

The [sphere period proof](CM_NEWTON_SPHERE_PERIOD_CYCLES.md)
gives a geometric version of coherent detection near the
separable phase by tracking actual monodromy cycles.

## Polynomial reconstruction across the radial line

**Corollary (Smith form and a polynomial left inverse).**
Let \(M(b)\) have rows \(0,\ldots,7,9\) and columns indexed by
\(a_0,a_1,c_0,\ldots,c_5\), with \(c_j=[z^j]C\).
Over \(\mathbb Q[b]\) its nonzero Smith entries are \(1,1\) and
six copies of \(b\). With \(u=b^2\),
\[
              M(b)=N(u)\operatorname{diag}(1,1,b,b,b,b,b,b),
\]
where the normalized matrix \(N(u)\) has a left inverse over
\(\mathbb Q[u]\). It reconstructs \((a_0,a_1,bc_0,\ldots,bc_5)\)
at every specialization, including zero.

In any row \(k\) let \(h\) count factors from \(S Q\). For the
coherent columns the even-height condition is \(h\equiv k\bmod2\).
For a correction column it is \(h\equiv k-1\bmod2\), so \(k-h\)
is positive and odd. This gives the asserted factorization.

Write \(N=(W\mid V)\) with coherent block \(W\), and let \(E\)
consist of rows \(1,\ldots,6\). The lower triangular block \(V_E\)
has constant nonzero diagonal by the preceding pivot formula
with \(b\) removed. Thus \(T=V_E^{-1}\) is polynomial. For
an arbitrary vector \(y\) of row values, put
\[
 \begin{aligned}
 a_0&=13y_0,&
 q&=y_E-W_{E,0}a_0,\\
 r_k&=y_k-W_{k,0}a_0-V_kTq \quad (k=7,9),&
 s_k&=W_{k,1}-V_kTW_{E,1}.
 \end{aligned}
\]
On \(y=N(u)(a_0,a_1,bc_0,\ldots,bc_5)^t\) these satisfy
\(r_k=s_k a_1\). The rows \(s_7,s_9\) are the two coprime
residual polynomials above, up to scalar normalization.
Take a Bézout identity \(\alpha s_7+\beta s_9=1\). Then
\(a_1=\alpha r_7+\beta r_9\) and the weighted correction
vector is \(T(q-W_{E,1}a_1)\). This is the polynomial left
inverse checked by full matrix multiplication.

Taking also \(-s_9r_7+s_7r_9\) as the ninth output augments
these operations to an invertible row transformation.
On the last two residual rows this uses determinant one.
The other operations have determinant \(13/\det V_E\),
a rational unit. Substituting \(u=b^2\) thus gives the
claimed Smith form of \(M(b)\). In particular the six
elementary factors have order exactly one at zero.
For the all-row marked-response operator the kernel is zero
for \(b\ne0\) and the six-dimensional correction space for
\(b=0\), since parity then kills the corrections in every row.
