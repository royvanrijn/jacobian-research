# Relative quadratic Hermite bounds for simplex transfers

Let \(T=(T_i)\) be the central regular \(k\)-simplex fan in
dimension \(d\ge k-1\), extended cylindrically. Let
\(f:\mathbb R^d\to\Delta_k\) be a measurable fractional assignment
with masses \(p_i\le B/k\), for a fixed \(B\ge1\).
Put \(e_i=f_i-1_{T_i}\), let \(P_m\) denote Gaussian Hermite projection
of degree \(m\), and write \(B_f,B_T\) for the matrices of cell
first moments. The linear-score regret is
\[
 \mathcal G=\|B_T\|_F^2-\langle B_T,B_f\rangle_F
    =\int\left(\max_i b_i^T\cdot x-\sum_i f_i(x)b_i^T\cdot x\right)
                     d\gamma_d(x),
\]
where \(b_i^T\) is the \(i\)-th row of \(B_T\).
For constants \(C_0,C_1\) depending only on \(B\), for all sufficiently
large \(k\) and every \(m\ge0\),
\[
 \sum_i\|P_m e_i\|_2^2
        \le C_0(C_1\log k)^m\mathcal G.
 \tag{RT1}
\]
The constants are independent of \(d\). In particular, for the
signed remainder \(R_r=\sum_{m\ge3}(-r)^mP_m\),
\[
 \left|\sum_i\langle e_i,R_r e_i\rangle\right|
     \le C r^3(\log k)^3\mathcal G
       \qquad\left(0<r\le\frac1{k-1}\right).
 \tag{RT2}
\]
The estimate concerns the quadratic error in the remainder expansion.
Replacing the remaining linear term against the simplex by an
estimate relative to its exact score cost is a separate step.

## Coordinates for ordered transfers

Work in the active label coordinates \(\mathbf1^\perp\), adding an
independent standard Gaussian in the common-coordinate direction so
that the simplex is the argmax of \(k\) independent standard
coordinates. Adding this coordinate leaves all Hermite degrees,
probabilities and regrets of the original functions unchanged.
Fix \(i\ne j\), let \(t=(x_i-x_j)/\sqrt2\), and write \(x=(t,y)\)
with \(y\) in the orthogonal complement of the unit normal.
Gaussian probability factors as \(\phi(t)dt\,d\gamma(y)\).
The ordered transfer and its regret are
\[
 h_{ij}=1_{T_i}f_j,\qquad
 g_{ij}=\ell\int t h_{ij}\,d\gamma,\qquad
 \ell=\|b_i^T-b_j^T\|=\sqrt2 k\, I_2(k).
\]
Here \(t\ge0\) on the transfer support and
\(I_j(k)=\int\phi(y)^j\Phi(y)^{k-j}\,dy\) as in
[(VH1)–(VH2)](NEGATIVE_VOTING_HESSIAN.md#the-exact-two-row-variance).
In particular \(z=\Phi^{-1}(1-1/k)\asymp\sqrt{\log k}\),
\(I_2(k)\asymp z/k^2\), and \(\ell\asymp z/k\).
Symmetry and the Gaussian flux formula give
\((b_i^T-b_j^T)\cdot x=\ell t\).
Thus \(\mathcal G=\sum_{i\ne j}g_{ij}\) and
\(e_i=\sum_{j\ne i}(h_{ji}-h_{ij})\) almost everywhere.

The part of \(h_{ij}\) with \(0\le t\le1/z\) projects onto a subset of
\[
 \Omega_{ij}=\{Z_l\le Y+1/(\sqrt2 z)\text{ for every }l\ne i,j\}
\]
in the tie plane. Here \(Y\) is the average of the two tied
coordinates in \(y\), distributed as \(N(0,1/2)\), and the \(Z_l\)
are independent standard Gaussians, independent of \(Y\).
Indeed \(x_i=Y+t/\sqrt2\) and \(x\in T_i\).
There is an absolute bound
\[
 \gamma(\Omega_{ij})
   =2\sqrt\pi\int \phi(Y)^2
         \Phi\big(Y+1/(\sqrt2 z)\big)^{k-2}\,dY
   \le C I_2(k) .
 \tag{RT3}
\]
For the integral put \(Y'=Y+1/(\sqrt2 z)\).
On \(|Y'|\le2z\), the ratio
\(\phi(Y'-1/(\sqrt2 z))^2/\phi(Y')^2\) is bounded by an
absolute constant. For \(Y'>2z\), the integral of the shifted
Gaussian square is \(o(I_2(k))\) since \(z^2\sim2\log k\).
For \(Y'<-2z\), use \(\Phi(Y')^{k-2}\le2^{-(k-2)}\).
The input density of \(Y\) is \(2\sqrt\pi\phi(Y)^2\), giving
the first equality.

## Relative estimates on each side of the transfer strip

Set \(q=\log k\) for large enough \(k\) that \(q\ge4\).
Let \(P\) be any element of Gaussian chaos of degree \(m\) with
\(\|P\|_2=1\). In the normal coordinate expand
\[
 P(t,y)=\sum_{j=0}^m H_j(t)P_{m-j}(y),
\]
where \(H_j\) are the orthonormal probabilist's Hermites and
the tangential functions \(P_{m-j}\) lie in chaos of degree \(m-j\).
Their squared \(L^2\) norms sum to one. The generating function
\(e^{tu-u^2/2}=\sum_{j\ge0}H_j(t)u^j/\sqrt{j!}\),
Cauchy's bound on the complex circle \(|u|=\sqrt j\),
and \(j!\le j^j\) give
\[
 \sup_{|t|\le1}|H_j(t)|
 \le \sqrt{j!}\,j^{-j/2} e^{\sqrt j+j/2}
 \le e^{3j/2}\quad (j\ge1).
\]
For \(j=0\) use \(H_0=1\).
By Gaussian hypercontractivity on the tangential space,
\(\|P_{m-j}\|_q\le(q-1)^{(m-j)/2}\|P_{m-j}\|_2\).
This is the chaos form of
[Nelson's inequality as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573).
The triangle inequality, then Cauchy–Schwarz in \(j\),
give
\[
 \left\|\sup_{|t|\le1/z}|P(t,\cdot)|\right\|_q
            \le (C\sqrt q)^m ,
 \tag{RT4}
\]
including \(m=0\), for which the left side is one.
Indeed the additional factor \(\sqrt{m+1}\) for \(m\ge1\)
can be absorbed into \(C^m\).

In the near part of an ordered transfer set
\[
 M(y)=\int_0^{1/z} h_{ij}(t,y)\phi(t)\,dt .
\]
The density in this integral is bounded by \(\phi(0)\).
For any density of mass \(M\) on the positive halfline with this
bound, its first moment is at least \(M^2/(2\phi(0))\).
Filling the interval closest to zero at density \(\phi(0)\)
minimizes the first moment, since the signed difference from that
density has zero integral and nonnegative pairing with the
distance from its right endpoint. Thus
\(\int M(y)^2\,d\gamma(y)\le 2\phi(0)g_{ij}/\ell\).
Apply Cauchy–Schwarz, then Hölder on \(\Omega_{ij}\) with (RT3)–(RT4).
Since \(I_2(k)^{-2/q}=O(1)\), the resulting bound is
\[
 \begin{split}
 \left|\int_{\{t\le1/z\}} h_{ij}P\,d\gamma\right|
 &\le \left(\frac{2\phi(0)g_{ij}}{\ell}\right)^{1/2}
       \left(\int_{\Omega_{ij}}
              \sup_{|t|\le1/z}|P(t,y)|^2\,d\gamma(y)\right)^{1/2}\\
 &\le C\left(\frac{(Cq)^m}{k}\,g_{ij}\right)^{1/2}.
 \end{split}
 \tag{RT5}
\]
The same bound applies to incoming ordered transfers by taking
the inward normal of their source label. For any fixed label,
summing (RT5) over its incoming or outgoing edges and using
Cauchy–Schwarz in at most \(k-1\) terms cancels the apparent
\(1/k\) in the bound on its squared interaction.

For the far part \(t\ge1/z\), group the incoming transfers of
a fixed label \(i\) into a single density. Its integral is at most
\(p_i\), since the reference source cells of the transfers are
disjoint. Their sum is bounded by one pointwise.
Likewise, the sum of the outgoing densities of \(i\) integrates
to at most \(1/k\) and is bounded by one.
Each individual far transfer has pointwise regret cost at least
\(\ell/z\), so the sum of its integrals is at most \(z/\ell\)
times the corresponding sum of \(g_{ij}\).
If \(H\) is such an incoming or outgoing far density, Hölder and
hypercontractivity give
\[
 \int H P^2\,d\gamma
   \le \left(\int H\,d\gamma\right)^{1-2/q}\|P\|_q^2
   \le \frac{C_B}{k}(Cq)^m .
\]
Cauchy–Schwarz against \(H\), together with \(z/\ell\le C k\),
therefore bounds the squared pairing of this grouped far density
with \(P\) by \(C_B(Cq)^m\) times its grouped regret.

Combine the two near and two far bounds for \(e_i\).
The squared pairing \(|\int e_iP\,d\gamma|^2\) is bounded by
\(C_B(Cq)^m\sum_{j\ne i}(g_{ij}+g_{ji})\).
Taking the supremum over the unit chaos and then summing in \(i\)
gives (RT1), since each \(g_{ij}\) is counted twice.
Finally \(C_1 r\log k\to0\) uniformly on the voting interval,
so summing the series for \(m\ge3\) gives (RT2).

For a possible improver in the polar alignment of
[(V18)–(V20)](NEGATIVE_VOTING_INTERVAL.md#dimension-uniform-localization-on-the-voting-scale),
\(\mathcal G=O(k\log(k)\Delta)\). Thus (RT2) makes the quadratic
part of the remainder \(o(r\Delta)\) uniformly throughout the
voting interval. This conclusion uses the actual deviations
\(e_i\) and does not require a minimum size for the label changes.
