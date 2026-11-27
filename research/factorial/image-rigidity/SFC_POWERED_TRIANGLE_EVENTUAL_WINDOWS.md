# Uniform high-order nonvanishing for powered factorial triangles

For \(F_D=a(xy)^D+b(yz)^D+c(zx)^D\) with integer \(D\ge2\), the actual support-sized factorial moment window is nonzero at every start. The proof assigns each term of a moment to a coefficient of maximal multiplicity. The three resulting contributions have a common holomorphic dependence on the selected coefficient. Divided differences then give an invertible matrix even when two or three coefficients approach one another.

This is the window condition of [Edo and van den Essen, Definition 2.7 and Conjecture 2.8](https://arxiv.org/pdf/1304.3956v2). The [separating-row theorem](SFC_SEPARATING_ROW_AMPLIFICATION.md#3-an-intrinsic-zero-pattern-theorem) already handles every start at \(D\ge36\). Its comparison controls the sum of the mixed terms. Here the three mixed-term contributions are kept attached to their respective pure powers across the whole window.

## Exact decomposition of the moments

Use the normalization
\[
          M_{m,D}(a,b,c)=\frac{\Lambda(F_D^m)}{((Dm)!)^2}.
\]
Each pure coefficient is one. If the exponent of \(a\) in a term is \(m-l\) and those of \(b,c\) are \(j,l-j\), the coefficient is
\[
 R_{m,D}(l,j)=
      \binom m l\binom l j
      \frac{(Dl)!\,(D(m-j))!\,(D(m-l+j))!}{((Dm)!)^2}.
 \tag{1}
\]
Divide each term equally among the coefficient positions with largest multiplicity, ignoring positions with smaller multiplicity. Write \(w_m(l,j)\in[0,1]\) for the share given to \(a\). It is zero unless \(m-l\ge\max(j,l-j)\), and is symmetric in \(j,l-j\). In particular \(l\le2m/3\) on that share.

Let \(s=a+b+c\), \(h=ab+bc+ca\). For a nonzero complex variable \(t\), let \(u,v\) be the unordered roots of
\[
                   z^2-(s-t)z+h-st+t^2=0 .
\]
Define
\[
 U_{m,D}(t;s,h)=1+
     \sum_{1\le l\le 2m/3}\sum_{j=0}^l
       w_m(l,j) R_{m,D}(l,j)\frac{u^j v^{l-j}}{t^l}.
 \tag{2}
\]
The numerator of each level is symmetric in \(u,v\) and is a polynomial in \(s-t,h-st+t^2\), so \(U\) is holomorphic off \(t=0\). At \(t=a\), the roots are \(b,c\). The exact decomposition is the sum of \(x^m U_{m,D}(x;s,h)\) over the three coefficient positions \(x=a,b,c\), with the expression defined by continuity at zero. Indeed each nonzero share has \(m-l\ge m/3>0\), so a zero coefficient contributes zero there.

## Bounds on the factorial layers

The log-convexity of the factorial gives
\[
  (D(m-j))!\,(D(m-l+j))!\le(Dm)!\,(D(m-l))!.
\]
Moreover \(\binom{Dm}{Dl}\ge\binom{2m}{2l}\binom m l^{D-2}\) by choosing in blocks of sizes \(2m,m,\ldots,m\). Consequently the absolute sum of level \(l\) in (2), when \(\max(|u/t|,|v/t|)\le\rho\), is at most
\[
        (2\rho)^l\,\frac{\binom m l}{\binom{2m}{2l}}\,
                                  \binom m l^{-(D-2)}.              \tag{3}
\]
For \(m/2<l\le2m/3\), all assigned terms have \(2l-m\le j\le m-l\). On this interval the factorial pair is largest at either endpoint, again by log-convexity. At \(j=m-l\), put \(k=l-j=2l-m\). The extra factor relative to the displayed factorial bound is at most
\[
       \frac{(D(m-j))!\,(D(m-k))!}{(Dm)!\,(D(m-l))!}
       =\prod_{v=1}^{Dj}\frac{D(m-l)+v}{D(m-j)+v}
       \le(1-k/m)^{Dj}.
 \tag{4}
\]
The last ratio in the product is its maximum. Thus (3) has this extra factor in the indicated range.

Write \(H(x)=-x\log x-(1-x)\log(1-x)\). For \(0<l<m\) the beta integral gives
\[
 \frac{\binom m l}{\binom{2m}{2l}}
   =(2m+1)\binom m l\int_0^1 z^{2l}(1-z)^{2m-2l}\,dz
   \le \frac{2m+1}{m+1}\exp(-mH(l/m)).
 \tag{5}
\]
Indeed \(z^l(1-z)^{m-l}\) is at most \(\exp(-mH(l/m))\), and its integral is \(1/((m+1)\binom m l)\), by integration by parts.

For \(D=2\), take \(\rho=3/2\). Denote the bound (3) by \(A_l\), ignoring (4) for small \(l\). The ratio \(A_{l+1}/A_l=3(2l+1)/(2m-2l-1)\) is at most \(1/2\) for \(l\le m/8\) and \(m\ge128\). The sum over \(4\le l\le m/8\) is at most \(2A_4\le17010/(2m-7)^4\). For \(x=l/m\) in \([1/8,1/2]\), (5) gives an exponential saving with exponent \(H(x)-x\log3\). On \([1/2,2/3]\), (4) adds \(-2(1-x)\log(2(1-x))\) to this exponent. The two functions are concave on their respective intervals. Their values at \(1/8,1/2,2/3\) are at least \(1/7\). These endpoint inequalities follow from
\[
       \frac{8^8}{3\cdot7^7}>e^{8/7},\qquad
       \frac43>e^{2/7},\qquad
       \frac{27}{16}>e^{3/7}.
\]
For example \(e^x<(1-x)^{-1}\) for \(0<x<1\) bounds the first and third exponentials by \((7/6)^8\) and \((14/11)^2\), respectively. For the second use \(e^x\le1+x+x^2/(2(1-x/3))\). The contribution of \(l\ge4\) to \(U_{m,D}\) therefore has modulus bounded by
\[
              E_2(m)=\frac{17010}{(2m-7)^4}+2m e^{-m/7}
                     \qquad(m\ge128).
 \tag{6}
\]

For \(D\ge3\), take \(\rho=5/3\). The bound (3) is at most \((10/3)^l/\binom m l^2\). For \(l\le m/5\) the ratios of these bounds are at most \(1/3\) when \(m\ge48\). The sum for \(4\le l\le m/5\) is at most \((3/2)(10/3)^4/\binom m4^2\). On \(l/m\in[1/5,2/3]\), use (5) and \(1/\binom m l\le(m+1)e^{-mH(l/m)}\). The latter follows because the term at \(l\) is a mode of the binomial probability distribution at parameter \(l/m\). The exponent \(2H(x)-x\log(10/3)\) is concave and is at least \(9/20\) at the two endpoints. At \(1/5\) this follows from \((5^5/4^4)^2/(10/3)>e^{9/4}\), using \(e^{9/4}<(4/3)^9\). At \(2/3\) it follows from \(6561/1600>e^{27/20}\), using the preceding bound for \(e^x\) at \(x=9/20\) and cubing. Hence the contribution of \(l\ge4\) has modulus bounded by
\[
           E_3(m)=\frac{3(10/3)^4}{2\binom m4^2}
                              +2m(m+1)e^{-9m/20}
                    \qquad(m\ge48).
 \tag{7}
\]

## Invertibility through coincident coefficients

**Theorem (Uniform terminal bound).** Let \(F_D=a(xy)^D+b(yz)^D+c(zx)^D\) be nonzero over a characteristic-zero field, and let \(r\) be its actual number of monomials. At least one moment in orders \(n,\ldots,n+r-1\) is nonzero for all \(n\ge128\) when \(D=2\) and for all \(n\ge48\) when \(D\ge3\).

Work over \(\mathbb C\) and normalize \(\max(|a|,|b|,|c|)=1\). Call a coefficient large if its modulus is at least \(\theta=9/10\) for \(D=2\) or \(\theta=4/5\) for \(D\ge3\). For each large value \(t_0\), use the complex disk \(|t-t_0|\le1/50\). On \(|z|=13/10\), the quadratic defining \(u,v\) at \(t_0\) has modulus at least \(9/100\). Moving by \(\delta=t-t_0\) changes its value by \(\delta(z-s+2t_0)+\delta^2\), of modulus at most \((1/50)(13/10+3+1/50)<9/100\). By Rouché's theorem, \(|u|,|v|<13/10\) on the disk. Also \(|t|\ge\theta-1/50\); these bounds give the respective ratio limits \(3/2\) and \(5/3\) used above.

Let \(t_1,\ldots,t_d\) be the distinct large values, with \(|t_d|=1\), and set \(g_j(t)=t^j U_{n+j,D}(t;s,h)\) for \(0\le j<d\). Write \(g(t)\) for their vector. Replace the value columns at the nodes by the successive Newton divided differences \(g[t_1],g[t_1,t_2],\ldots,g[t_1,\ldots,t_d]\). For the pure terms \(1,t,t^2\) the resulting matrix is the top-left block of
\[
 H=\begin{pmatrix}1&0&0\t_1&1&0\t_1^2&t_1+t_2&1\end{pmatrix}.
\]
The maximum absolute row sum of its inverse is at most four.

To bound the perturbation, first consider \(1\le l\le3\), where the largest multiplicity \(m-l\) is strict. Write \(B_l(m)=\binom m l/\binom{2m}{2l}\). The coefficients of \(u^l+v^l\) in level \(l\) before division by \(t^l\) are at most \(B_l(m)\) for \(D=2\) and \(B_l(m)/\binom m l\) for \(D\ge3\). For \(l=2,3\), the accompanying coefficients of \(uv\) and \(uv(u+v)\) are respectively at most \(2\) and \(3\) times those same bounds, by (1). For \(l=1\) the numerator is just the coefficient of \(u+v\) times \(u+v\).

The Laurent polynomials for these levels, before multiplying by \(t^j\), are
\[
\begin{array}{ll}
 l=1:& (u+v)/t=s/t-1,\
 l=2:& (u^2+v^2)/t^2=(s^2-2h)t^{-2}-1,\qquad
               uv/t^2=h t^{-2}-s t^{-1}+1,\
 l=3:& (u^3+v^3)/t^3=(s^3-3hs)t^{-3}+3h t^{-2}-3s t^{-1}+2,\qquad
         uv(u+v)/t^3=hs t^{-3}-(h+s^2)t^{-2}+2s t^{-1}-1 .
\end{array}
\]
For a negative power \(t^{-v}\), \(v>0\), the absolute divided difference on \(\ell+1\) nodes is at most \(\binom{v-1+\ell}{\ell}\theta^{-v-\ell}\). This follows by expanding the divided difference of \(1/(t-z)\) as a series in \(z\) and taking its \(z^{v-1}\) coefficient. For \(t^0,t^1,t^2\) the row-sum bounds are at most four, by the explicit monomial divided differences. Therefore the Laurent monomials with exponents \(-3,\ldots,2\) have row-sum bounds \(17\) at \(\theta=9/10\) and \(28\) at \(\theta=4/5\). For exponents \(-1,\ldots,2\) these bounds can be replaced by \(4,5\), respectively. With \(|s|,|h|\le3\), the resulting row sums of the low-level perturbation are bounded by
\[
\begin{split}
 D=2:\quad&16 B_1(n)+510 B_2(n)+2686 B_3(n),\
 D\ge3:\quad&20 B_1(n)/n+840 B_2(n)/\binom n2+
                                      4424 B_3(n)/\binom n3 .
\end{split}                                                       \tag{8}
\]
For instance the sum of the bounds on the scalar coefficients from the \(l=3\) Laurent polynomials is \(74+3\cdot28=158\) times the pure-layer bound. In each row \(m=n+j\), and the coefficient bounds decrease with \(m\).

For the remaining levels use (6) and (7), multiplying by \((51/50)^2\) for \(t^j\). For any column, integrate the error divided by \(\prod_i(z-t_i)\) on the boundary of the union of equal-radius disks around the nodes in that column. Choose a radius \(R\) between \(19/1000\) and \(1/50\) avoiding tangencies. Each boundary point has distance at least \(R\) from every included node and the total arc length is at most \(2\pi R\) times the number of nodes. The contour formula for divided differences thus gives a row-sum factor at most \(1+2/R+3/R^2<8500\). The full perturbation is bounded by the applicable row of (8) plus \(8500(51/50)^2 E_i(n)\), \(i=2,3\). Multiplying by four gives less than \(3/4\). The bounds decrease in \(n\) on both asserted ranges, and the lower endpoints follow by the rational checks in the [exact checker](verify_sfc_powered_triangle_finite_windows.py). In those checks, one can bound \(e^{-128/7}<20^{-6}\) from \(e^3>20\), and \(e^{-9\cdot48/20}<11^{-9}\) from \(e^{12/5}>11\). The resulting divided-difference matrix is invertible, with inverse row-sum norm at most \(16\).

A small coefficient position contributes at most \(2\theta^m\) at moment \(m\). Indeed bound its \(m-l\) factors by \(\theta^{m-l}\) and the others by \(1\), and apply the same layer estimates at \(\rho=1/\theta\), which lies below the respective ratio limit. The sum of the pure contribution, the first three bounds in (3), and (6) or (7) is less than two on the asserted ranges. Thus the vector of small contributions has maximum norm at most \(4\theta^n\). If all \(r\) moments vanish, the coefficients of the sum of the large contributions in the divided-difference basis have magnitude at most \(64\theta^n<1/n\); the last inequality is checked at the lower endpoints and \(n\theta^n\) decreases thereafter.

Each large value carries weight \(w_i=e_i t_i^n\), where \(e_i\) is its positive multiplicity among the coefficients. Newton interpolation gives coefficient \(\sum_i w_i\) in the first divided-difference column and \(w_d\prod_{i<d}(t_d-t_i)\) in the last. If all active positions are large, small contributions vanish exactly, and the latter nonzero coefficient contradicts the zero window. Otherwise there are at most two large positions. With one distinct large value its total weight has modulus at least one, a contradiction. With two distinct large values each has multiplicity one, and the last divided-difference coefficient gives \(|t_2-t_1|<1/n\). Then \(|t_2^n-t_1^n|\le n|t_2-t_1|<1\) and \(|t_2|=1\), so the first coefficient \(|t_1^n+t_2^n|\) exceeds one, also a contradiction. The matrix uses at most \(r\) rows. The field extension follows by embedding the finitely generated field of the coefficients into \(\mathbb C\).

## The first two windows and the finite residual

**Proposition (Initial window).** For every integer \(D\ge1\), the support-sized factorial window at \(n=1\) is nonzero for every nonzero powered triangle.

The first moment is \(s\), and the second is \(s^2+(u-2)h\) with \(u=2/\binom{2D}{D}<2\). With two active terms, zero first and second moments would give \(s=h=0\), impossible. With three active terms, the same two zeros force \(s=h=0\) and \(r_0=abc\ne0\). The normalized third moment is then \(r_0(3-3v+w)\), where \(v=3/\binom{3D}{D}\) is the coefficient of each \(a^2b\) type and \(w=6((2D)!)^3/((3D)!)^2\) is the coefficient of \(abc\). Since \(v\le1\) and \(w>0\), the third moment is nonzero. A single active term has a nonzero first moment.

**Proposition (Second window).** For every integer \(D\ge2\), the support-sized factorial window at \(n=2\) is nonzero for every nonzero powered triangle.

Retain \(u,v,w\) from the initial window. When \(P_2=M_{2,D}=0\), a zero value of \(s\) gives \(h=0\). For two active terms this is impossible. For three terms with \(s=h=0\) the third moment is nonzero by the preceding calculation. Thus assume \(s=1\) after scaling, and \(h=1/(2-u)\). The third moment is
\[
                       1+(v-3)h+(3-3v+w)r_0 .
\]
With two active terms \(r_0=0\), it cannot vanish since \(3-v>2-u\). With three terms, its zero determines \(r_0=((3-v)h-1)/(3-3v+w)>0\). Put
\[
 \gamma=\frac4{\binom{4D}D},\qquad
 \eta=\frac6{\binom{4D}{2D}},\qquad
 \tau=\frac{12(2D)!\,((3D)!)^2}{((4D)!)^2}.
\]
These are respectively the coefficients of \(a^3b,a^2b^2,a^2bc\) in the normalized fourth moment. Newton's identities give that moment as
\[
               1+(\gamma-4)h+(2-2\gamma+\eta)h^2+
                                     (4-\gamma-2\eta+\tau)r_0 .
 \tag{7}
\]
At \(D=2\) the substitutions give \(h=3/5,\ r_0=17/64\), and (7) equals \(264307/784000>0\).

For \(D\ge3\), the bounds \(1/2\le h\le3/5\), \(v\le1/28\), \(w\le1/50\), \(\gamma\le1/55\), and \(\eta\le1/154\) follow by monotonicity of the binomial coefficients and of \(w=6((2D)!)^3/((3D)!)^2\). For \(w\), in the ratio at \(D+1,D\) each factor \((2D+1)(2D+2)/(3D+i)^2\), \(i=1,2,3\), is at most one. Thus \(4-\gamma-2\eta+\tau\ge39/10\) and \(r_0\ge((83/28)h-1)/(151/50)\). In (7) the contribution \(\gamma(h-2h^2)\) is at least \(-3/1375\), and \(\eta h^2\ge0\). These estimates give the lower bound
\[
               1-4h+2h^2-\frac3{1375}
                         +\frac{195}{151}\left(\frac{83}{28}h-1\right).
\]
It is increasing on \([1/2,3/5]\) and positive at \(1/2\). This proves the second window. A single active term again has a nonzero moment.


## Exact closure of the finite range

**Corollary (Finite reduction for powered triangles).** For integer \(D\ge2\), the all-order claim is equivalent to the support-sized window test for \(D=2,\ 3\le n\le127\) and \(3\le D\le35,\ 3\le n\le47\).

The analytic theorem and the first two windows leave \(3\le n\le127\) for \(D=2\) and \(3\le n\le47\) for \(3\le D\le35\). For \(D\ge36\) the [separating-row theorem](SFC_SEPARATING_ROW_AMPLIFICATION.md#3-an-intrinsic-zero-pattern-theorem) applies: the exponent columns have the same positive-entry multiset and two rows distinguish their binary columns. It provides all windows at those exponents, including the smaller supports.

**Theorem (All-order nonvanishing for powered triangles).** Over every characteristic-zero field, every nonzero \(F_D=a(xy)^D+b(yz)^D+c(zx)^D\) with integer \(D\ge2\) satisfies the support-sized Strong Factorial window at every start \(n\ge1\).

The finite range is certified in \(\mathbb F_{100003}\). For each order \(m\) in a tested window, the polynomial with coefficients
\[
         w_{m,i}w_{m,j}w_{m,k},\qquad
         w_{m,i}=\frac{(D(m-i))!}{i!},\quad i+j+k=m
 \tag{9}
\]
is the actual moment polynomial up to a nonzero scalar modulo the prime. All factorial arguments in (9) are less than the prime. Use these homogeneous forms in the following tests.

- For two active monomials, the chart \((a,b,c)=(1,b,0)\) suffices by symmetry. The two univariate polynomials of consecutive orders have gcd one.
- For the three-moment window, use the symmetric coordinates \(s=a+b+c,\ h=ab+bc+ca,\ r_0=abc\). Where \(s\ne0\), normalize to \(s=1\). The three polynomials in \((h,r_0)\) generate the unit ideal.
- Where \(s=0\), some coefficient can be scaled to \(a=1\) by symmetry. With \(c=-1-b\), one has \(h=-1-b-b^2\) and \(r_0=h+1\). The three polynomials after substituting \(s=0,\ r_0=h+1\) have gcd one in \(h\).

The symmetric formulas in these charts can be generated directly from (9). Group terms by the sorted triple \(i\ge j\ge k\). If all are equal, the orbit contributes \(r_0^k\). Otherwise remove \(r_0^k\) and put \(I=i-k,\ J=j-k\). For \(J=0\) the orbit sum is \(p_I=a^I+b^I+c^I\). For \(J>0\) it is \(p_Ip_J-p_{I+J}\), divided by two when \(I=J\). The \(p_I\) follow the cubic recurrence with coefficients \(s,h,r_0\), starting with \(p_0=3,\ p_1=s,\ p_2=s^2-2h\). The checked univariate gcds and the Gröbner bases \([1]\) are exact over the finite field. Together the tests cover every field extension of \(\mathbb F_{100003}\), including projective points with some coordinates zero.

A complex counterexample yields one over \(\overline{\mathbb Q}\) by the weak Nullstellensatz in a nonzero-coefficient chart. At a place over \(100003\), choose integral projective coordinates with at least one unit. Reduction preserves the homogeneous moment equations and the nonzero projective point. For a two-term counterexample use the two homogeneous equations restricted to its two coefficient positions, so its reduction is covered by the two-term test. A one-term polynomial has every moment nonzero by the factorial formula. The finite-field certificates therefore establish the desired characteristic-zero exclusions.

The [exact checker](verify_sfc_powered_triangle_finite_windows.py) verifies all \(125+33\cdot45=1610\) starts in the finite range with SageMath 10.9. It also evaluates the rational constant inequalities at the analytic cutoffs. To replay both parts from the repository root:

~~~sh
sage -python research/factorial/image-rigidity/verify_sfc_powered_triangle_finite_windows.py
~~~

For \(D=1\), the [quadratic triangle recurrence](SFC_TRIANGLE_OVERLAP_MOMENTS.md#remaining-sector) supplies separate windows and reductions.
