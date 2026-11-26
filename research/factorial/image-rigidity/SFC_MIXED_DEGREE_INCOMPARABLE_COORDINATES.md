# Mixed-degree saddle windows with incomparable coordinate supports

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer exponent
columns, \(d_i=\sum_\nu b_{\nu i}\), and \(S_i=\{\nu:b_{\nu i}>0\}\).
Suppose that whenever \(d_i\ne d_j\), neither \(S_i\) nor \(S_j\)
contains the other.

**Theorem.** For all sufficiently large integers \(D\), the ordinary
factorial moments of \(\sum_i c_iX^{D b_i}\) satisfy the actual
support-size window at every positive starting order, for every
coefficient vector over a characteristic-zero field. The threshold
depends on the columns. The same conclusion holds uniformly for
actual integer exponent columns \(a_i=Db_i\) at a real scale \(D\)
when positive normalized entries range in a fixed compact positive
interval, centers have a fixed positive separation, and distinct
normalized total degrees have a fixed positive gap.

Incomparable coordinate sets make the saddle phase at each peak
equal to the phase of its own homogeneous block. A foreign block
vanishes to high order in at least one inactive coordinate there.
The resulting block-specific order dependence can be incorporated
into a finite-difference filter by using an interpolation matrix
whose groups occupy separated annuli.

## Two local estimates

Use the peak quantities \(P_i,q_i(t)\) of
[the separated-support theorem](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md)
and a peak separation constant \(\gamma>0\) in its estimate
\(|P_j/P_i|\le e^{-4\gamma}\). We may reduce \(\gamma\) by any fixed
factor. Let \(g>0\) be a lower bound on distinct normalized degree
gaps. Discard zero coefficients and call their remaining number
\(s'\). By the
[exponential-range mixed-degree theorem](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order),
it suffices to consider windows beginning at \(m>e^{\eta D}\) for
some fixed \(\eta>0\). If there is a single degree block the
[homogeneous mixed-active theorem](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md)
already applies.

For \(0\le\ell<s'\) write \(t=m+\ell\), \(Q=\max_i|q_i(m)|\),
and \(d_*=\max_i d_i\). A peak is eligible if
\(|q_i(m)|\ge e^{-\gamma D}Q\). Over these shifts every ratio
\(q_i(t)/q_i(m)=(t/m)^{Dd_i}\) has modulus \(1+O(D/m)\).
Consequently all noneligible boxes and the exterior of the
peak boxes give \(O(e^{-cmD})\) after normalization by \(Q^t\)
in the scaled integral.

For an eligible \(i\), let \(\Phi_i^0\) be the phase in its
own degree group, and \((\theta_i,0)\) its saddle with inactive
coordinates zero. It is independent of \(t\). In the polydisc
there, every foreign quotient power \((P_j/P_i)^D\) vanishes to
order \(a_{\nu j}=Db_{\nu j}\) in an inactive coordinate \(\nu\)
of \(i\). On a smaller polydisc, Cauchy's formula gives the
bound \(C_L e^{-4\gamma D}|y_\nu|^L\) for any fixed \(L\) and
sufficiently large \(D\). This follows by dividing the
holomorphic quotient power by \(y_\nu^{a_{\nu j}}\), bounding
on the circle of the larger radius, and then retaining only
the first \(L\) powers of the radius ratio. The perturbation
\(D(\Phi_i(y,t)-\Phi_i^0(y))\) is therefore bounded by
\(C_L e^{-3\gamma D}(\sum_{\nu\notin S_i}|y_\nu|)^L\) on all
local contours. Its real part is absorbed by the linear
inactive-coordinate loss on the translated contour, since
the inactive coordinates stay real and nonnegative. The
absorption also holds along the straight segment of phase
interpolations. The contour sides retain their fixed
exponential loss.

Use the common contour for the two integrals and write the
difference of the two exponentials as a segment integral.
This gives a local error bounded by
\[
 C_L e^{-3\gamma D} t(tD)^{-L}
     |\hat q_i|^t (tD)^{-\rho_i}
     \exp(tD\operatorname{Re}\Phi_i^0(\theta_i,0)),
 \qquad \rho_i=r-\tfrac12|S_i|,
\]
where \(\hat q_i=q_i(t)/Q\). The factor \((tD)^{-L}\) follows
by integrating the \(L\)-th inactive-coordinate moment against
its endpoint exponential. Side contours satisfy the same
bound after increasing \(C_L\). No derivative in the order
parameter is used for the error. The homogeneous term has the
usual analytic sector amplitude
\[
 A_i(t)=(tD)^{-\rho_i}B_i(tD),\qquad
 B_i(\lambda)=K_i+O(|\lambda|^{-1/2})
\]
from the
[homogeneous sector estimate](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#a-uniform-analytic-saddle-amplitude).
In particular
\[
 |\Delta^h A_i(m)|\le C_h(mD)^{-\rho_i}m^{-h},\qquad
 \frac{m^h\Delta^h A_i(m)}{(mD)^{-\rho_i}}
        =(-1)^h(\rho_i)^{\overline h}K_i+o(1),                 \tag{1}
\]
and \(K_i\) tends to its positive pure amplitude. These bounds
are uniform on the compact shape sets in the theorem.

After removing the nonzero scalars \((tD)^rQ^t
 e^{D d_*\ell}\) from the original moment, its normalized
value \(M_{m+\ell}\) thus has the expansion
\[
 \begin{split}
 M_{m+\ell}
   &=\sum_{i\ {\rm eligible}} u_i^m v_i^\ell
                  \sigma_{G_i}(\ell) A_i(m+\ell) +R_\ell+\sum_i E_{i\ell},\\
 u_i&=\frac{q_i(m)}{Q}e^{D\Phi_i^0(\theta_i,0)},\qquad
 v_i=u_i e^{D(d_i-d_*)},\\
 \sigma_G(\ell)&=\exp\!\left(Dd_G\big[(m+\ell)
                      \log(1+\ell/m)-\ell\big]\right).
                                                               \tag{2}
 \end{split}
\]
Here \(G_i\) denotes the degree group of \(i\).
The bounds are \(R_\ell=O(e^{-cmD})\) and, for any fixed \(L'\),
\[
 |E_{i\ell}|\le C_{L'}e^{-3\gamma D}(mD)^{-L'}
                |u_i|^m |v_i|^\ell(mD)^{-\rho_i},
\]
by choosing \(L\) sufficiently large. The factors
\(\sigma_G(\ell)=1+O(D/m)\) are uniformly bounded above and below.
All \(v_i\) in group \(G\) lie in the annulus
\[
         \tfrac12 e^{-\gamma D} R_G \le |v_i| \le 2 R_G,
           \qquad R_G=e^{D(d_G-d_*)}.
                                                               \tag{3}
\]
There is a raw maximizing eligible label with
\(|u_i|\ge e^{-C e^{-3\gamma D}}\).

## Interpolation on separated annuli

The following observation applies to any multiset of at most \(s'\)
roots \(w_j\) from the annuli (3), each assigned its degree group.
Multiplicity is treated by requiring derivatives to vanish.
Let \(n\) be the size of this multiset, and choose a monic
coefficient vector \(p_0,\ldots,p_n\) by requiring
\[
                    P_G(z)=\sum_{\ell=0}^n
                                  p_\ell\sigma_G(\ell)z^\ell
\]
to have the prescribed roots in each group \(G\). Take \(\gamma\)
sufficiently small in terms of \(s,g\). This system has a unique
solution with uniformly bounded coefficients. If \(n_G\) is the
number prescribed in \(G\) and \(L_G=\sum_{K<G} n_K\), with groups
ordered by degree, then on a fixed enlargement of its annulus
\[
 P_G(z)= c_G z^{L_G}
              \prod_{K>G}\prod_{j\in K}(-w_j)
              \prod_{j\in G}(z-w_j)(1+O(e^{-\kappa D})),
                                                               \tag{4}
\]
where \(c_G\) is a positive real number bounded above and below.
The holomorphic relative error is uniform on the annulus.

Here is the determinant proof. For initially distinct roots
within each group, the coefficient matrix has rows
\((\sigma_G(\ell) w_j^\ell)_{\ell=0}^{n-1}\). Divide its
determinant by the product of Vandermondes within groups and
expand over partitions of column indices into subsets of sizes
\(n_G\). Within any block the alternant in \(w_j\) is divisible
by its Vandermonde, since it vanishes when two variables agree.
The quotient is a fixed polynomial of degree
\(\sum_{\ell\ {\rm assigned}}\ell-n_G(n_G-1)/2\). In the
partition assigning consecutive exponents to each group in
increasing order this quotient is \(\prod_{j\in G}w_j^{L_G}\).
After setting \(w_j=R_G x_j\), every other partition has radial
scale smaller by at least \(e^{-Dg}\): swapping any inverted
pair of exponents strictly increases the radial exponent by at
least that much. The quotient polynomials in the \(x_j\) are
bounded on \(|x_j|\le 3\), and the magnitude of the consecutive
term at worst loses \(C_s e^{-s^2\gamma D}\). The \(\sigma_G\)
factors are bounded above and below. Taking \(2s^2\gamma<g\)
leaves an exponential relative error.

The augmented determinant with a row
\((\sigma_G(\ell)z^\ell)_{\ell=0}^n\), divided by the original
one, is \(P_G(z)\). Apply the same estimate after assigning
the new row to group \(G\). The additional Vandermonde factors
are exactly \(\prod_{j\in G}(z-w_j)\), and the consecutive
radial monomials give \(z^{L_G}\prod_{K>G}\prod_{j\in K}w_j\).
The ordering signs yield the displayed minus signs. Ratios of
products of \(\sigma\) give \(c_G\). This proves (4) and its
holomorphic error after canceling the exact root factors.
Cramer's rule gives the coefficient bounds by the same
expansion with exponent sets \(\{0,\ldots,n\}\setminus\{\ell\}\).
The leading, consecutively assigned sets now have at most one
omission, so their quotients contain the corresponding
denominator monomials, with only bounded extra powers and
fixed polynomial factors on bounded \(w\). Every other
assignment has the same \(e^{-Dg}\) radial loss, with an
additional nonincreasing radial factor because every assigned
exponent in the ordered list is at least its denominator
counterpart and \(R_G\le1\). Finally all divided determinants
extend continuously to repeated roots by within-block
polynomial division. Their nonzero leading estimates persist.
Taking limits of the original divisibilities gives the
multiplicity conditions.

## Relative clusters and the window contradiction

Suppose a sequence \(D\to\infty\) gives vanishing windows in
(2). After subsequences fix eligibility and group labels and
positive pure-amplitude limits. Define
\[
 d_{ij}=|v_i-v_j|+|v_i|/m,\quad V_i=\prod_{j\ne i}d_{ij},\quad
 W_i=|u_i|^m(mD)^{-\rho_i}V_i.
\]
Choose a maximizing label \(p\), fixed on a subsequence.
Let \(I\) be the indices where \(m|v_i/v_p-1|\) stays bounded,
and take all other such distances to infinity by subsequences.
The cluster \(I\) belongs to the single degree group of \(p\)
by (3). Keep the limits \(w_i=\lim m(v_i/v_p-1)\) for \(i\in I\).
Since \(u_i/u_p=v_i/v_p\) there, maximizing \(W_p\) forces
\(\rho_i\ge\rho_p\) in \(I\), and terms with larger \(\rho_i\)
have vanishing relative amplitude. This is the relative-cluster
comparison of the homogeneous proof.

For each \(0\le k<|I|\), prescribe outside frequencies \(v_j\)
as roots, in their respective groups, and prescribe \(k\)
further copies of \(v_p\) in its group. The resulting
interpolation coefficient vector has degree at most \(s'-1\),
so its scalar product with the vanishing window is zero.
On the \(i\)-th mode this scalar product is exactly the
ordinary Newton difference formula for the polynomial
\(P_{G_i}\), because \(\sigma_{G_i}(\ell)\) is already in its
coefficients. In detail its value is
\[
                  u_i^m\sum_{a=0}^{n}
          \frac{v_i^a P_{G_i}^{(a)}(v_i)}{a!}\Delta^a A_i(m).
                                                               \tag{5}
\]
Let \(C_{p,k}\) be the value at \(v_p\) of
\(P_{G_p}(z)/(z-v_p)^k\). By (4) its modulus is comparable to
the product of all outside distances from \(p\). Divide the
tested window by
\(u_p^m(mD)^{-\rho_p} C_{p,k}(v_p/m)^k\). Its modulus is
comparable to \(W_p (|v_p|/m)^{k-|I|+1}\), as in the
homogeneous proof.

For outside \(i\), any nonzero differentiated term in (5)
must remove the factor \(z-v_i\) of its group. Derivatives of
(4) from the product of ordinary within-group factors give
their corresponding padded distances. The monomial part and
the holomorphic relative error have Cauchy bounds on a
circle of radius \(|v_i|/4\), so each derivative there costs
at most \(C/|v_i|\) and comes with the damping
\(|v_i|/m\) from (1). The external factors in (4) are
comparable to the products of cross-group distances by (3).
Consequently the outside contribution is bounded by
\[
           C W_i\, (|v_i|/m)\,
                 d_{ip}^k\big/\prod_{j\in I} d_{ij}.
\]
Upon normalization it tends to zero, since \(W_i\le W_p\),
\(d_{ij}\sim d_{ip}\) for \(j\in I\), and
\((|v_i|/m)/d_{ip}\), \((|v_p|/m)/d_{ip}\) both tend to
zero for outside \(i\).

Inside \(I\), the other factors of \(P_{G_p}\) have normalized
values tending to one and scaled derivatives of positive order
tending to zero. This follows from (4) and Cauchy's bounds for
its holomorphic factor and from the outside-distance condition
for the same-group roots. Thus after normalization the limit of
(5), using (1), is
\[
      \sum_{\substack{i\in I\\\rho_i=\rho_p}}\!
             K_i^*e^{w_i}\sum_{a=0}^k
             \binom ka\,w_i^{k-a}(-1)^a(\rho_p)^{\overline a},
                                                               \tag{6}
\]
where \(K_i^*>0\).

The errors of (2) also vanish after normalization.
Indeed all eligible \(|v_i|\ge e^{-C D}\), hence
\(V_i\ge e^{-C_sD}m^{-(s-1)}\). The bounds on the filter
coefficients and \(W_i\le W_p\) control the local errors
as \(O(e^{C_sD}m^{s-1}(mD)^{-L'}W_p)\), which is \(o(W_p)\)
by \(m>e^{\eta D}\) and a sufficiently large fixed \(L'\).
For \(R_\ell\), a raw maximizing index gives
\(W_p\ge c(mD)^{-r}e^{-C_sD}m^{-(s-1)}
 e^{-C m e^{-3\gamma D}}\), which dominates the exponential
remainder. The extra factor
\((|v_p|/m)^{k-|I|+1}\) in the normalizer is bounded below.

All expressions (6) must therefore vanish. Their inner
polynomials are monic of consecutive degrees \(0,\ldots,|I|-1\).
After grouping equal \(w_i\), the ordinary Vandermonde matrix
forces each weight \(e^{w_i}\sum K_j^*\) to vanish. Each is
nonzero, giving the contradiction.

The excluded very large crossover for general mixed supports
can involve foreign monomials nonzero at a saddle of another
degree group. Those terms introduce order-dependent saddle
values. The filter argument here uses the cross-group
coordinate hypothesis to make each saddle value independent
of the starting-order shift.

