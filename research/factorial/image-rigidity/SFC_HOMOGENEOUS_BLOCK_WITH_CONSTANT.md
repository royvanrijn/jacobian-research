# An arbitrary constant adjoining a dilated homogeneous support

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer exponent
columns with one positive total degree \(d\). For all sufficiently
large integers \(D\), depending only on the columns, the ordinary
factorial moments of
\[
                    c_0+\sum_{i=1}^s c_i X^{Db_i}
\]
have a nonzero moment in every positive window of the actual
support size, uniformly over all coefficients in characteristic
zero. If \(c_0=0\), this is the
[homogeneous theorem](SFC_HOMOGENEOUS_MIXED_ACTIVE_SADDLE_WINDOWS.md).
The same existence of a uniform threshold holds on compact
normalized homogeneous shape families with positive lower bounds
on nonzero entries and on pairwise center distances.

The constant contributes at every saddle of the homogeneous
block. The active critical value nevertheless changes by a
universal function of the original saddle value, because the
perturbation changes only the radial critical equation. The
exact mode filter can incorporate this function even for
arbitrarily close effective frequencies. In particular the
proof covers the large crossover where the constant has
coefficient size depending without restriction on the starting
order.

## Radial change of the full saddle

Work over \(\mathbb C\) and remove zero coefficients. The
constant-only case is immediate. The
[exponential-range theorem](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order)
reduces the problem to a possible sequence \(D\to\infty\) of
vanishing windows starting at \(m>e^{\eta D}\), with fixed
active support size \(s'\le s+1\) and fixed zero patterns.
Set \(t=m+\ell\), where
\(0\le\ell<s'\), and otherwise use the scaled peak variables
\(P_i,q_i(t)\), the normalization \(Q=\max_i|q_i(m)|\), and
the eligibility threshold \(e^{-\gamma D}Q\) of the
[incomparable-coordinate theorem](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates).
The constant peak is at the origin and \(q_0=c_0\). The
peak separation exponent \(\gamma>0\) may be chosen
arbitrarily small in terms of the fixed shape and \(\eta,d,s\).
After removing the common factor \((tD)^r Q^t e^{Dd\ell}\)
the exterior and noneligible contributions are \(O(e^{-cmD})\).
In fact all raw ratios \(q_i(t)/q_i(m)\) vary by \(1+O(D/m)\)
across the shifts, so the same real peak bounds give this
estimate, including when all eligible indices lie in one block.

In the normalized integral before the scalar \(e^{Dd\ell}\),
write the full local function as
\[
 e^{-D\sum_\nu y_\nu}
       \left((t/m)^{Dd} U(y)+c_0/Q\right),\qquad
 U(y)=\sum_{i\ge1}\frac{c_i(mD)^{Dd}}{Q}y^{Db_i}.
 \tag{1}
\]
At an eligible positive-degree label \(i\), let \(\theta_i^0\)
be the saddle of the \(U\) contribution, with inactive
coordinates zero. The logarithmic derivatives in every
active coordinate satisfy
\(\partial_\nu U(\theta_i^0)=D U(\theta_i^0)\).
The block saddle is nonzero by the local relative peak
bound. Euler's identity for \(U\) implies
\(\sum_\nu\theta_{i,\nu}^0=d\). Put
\[
       w_i=e^{-Dd}U(\theta_i^0),\qquad
       \kappa=(c_0/Q)e^{-Dd}.
\]
Then \(w_i=(q_i(m)/Q)\exp(O(e^{-3\gamma D}))\)
and \(|\kappa|\le e^{-Dd}\).

At the full saddle (1), the active point is
\((d+s_i)/d\) times \(\theta_i^0\).
Indeed its tangential logarithmic derivatives vanish by
homogeneity, and the active radial equation is
\[
                  s_i=-\frac{d B}{(1+s_i/d)^{Dd}+B},
          \qquad B=(\kappa/w_i)(t/m)^{-Dd}.
 \tag{2}
\]
Here and below the letter \(s_i\) is a displacement.
The solution \(s(B)=O(B)\) is unique and holomorphic for
\(|B|\le a/D\) with \(a\) small enough, by contraction on
\(|s|\le 2d|B|\); the derivative of the right side is
\(O(D|B|)\). The resulting point is the unique local
saddle from the active contraction on the peak box. Set
\[
          g(B)=\log((1+s(B)/d)^{Dd}+B)-D s(B).
\]
Uniformly on smaller such discs \(g(B)=O(B)\),
\(g'(B)=O(1)\). For the first bound separate
\(Dd\log(1+s/d)-Ds=O(D|s|^2)\) from the remaining
logarithm. The second follows by Cauchy's formula.
The full saddle value inside (1) is exactly
\[
                          (t/m)^{Dd}w_i \exp(g(B)).
 \tag{3}
\]

Define \(B_i=\kappa/w_i\) and remove the linear dependence
on the offset in \(t g(B_i(t/m)^{-Dd})\). The resulting
base and effective frequency are
\[
                 U_i=w_i e^{g(B_i)},\qquad
                 v_i=w_i e^{g(B_i)-Dd B_i g'(B_i)} .
 \tag{4}
\]
On \(e^{-2\gamma D}<|w|<8\), the map
\(w\mapsto w\exp(g(\kappa/w)-Dd(\kappa/w)g'(\kappa/w))\)
differs from the identity by \(O(D|\kappa|)\) and its
derivative differs by \(O(D|\kappa|/|w|)\).
It has a holomorphic inverse \(w=w(v)\) on a fixed radial
enlargement of the eligible \(v_i\) annulus. In fact
\(w=v-(v(w)-w)\) is a contraction in discs of radius
\(C D|\kappa|\) centered at \(v\). The branches agree by
uniqueness, and stay in the original annulus. All these
annuli can also be taken to include the two circles used
for Cauchy integration below, since every eligible
\(|w_i|\ge \frac12 e^{-\gamma D}\) for large \(D\).

Writing \(T=\ell/m\) and \(B=\kappa/w(v)\), define
\[
 H(v,T)=(1+T)g(B(1+T)^{-Dd})
       -g(B)-T[g(B)-Dd B g'(B)] .
 \tag{5}
\]
This is holomorphic in \(v,T\), with zero constant and
linear coefficients in \(T\). On \(|T|\le a_1/D\)
it is \(O(e^{-cD})\) on an enlarged annulus in \(v\),
for some \(c>0\), since \(|B|\le e^{-(d-2\gamma)D}\)
and \((1+T)^{-Dd}\) is bounded there. There are uniform
Cauchy bounds for any fixed number of logarithmic \(v\)
derivatives on smaller annuli, using discs of radius
proportional to \(|v|\). Put
\[
 \begin{split}
 F_\ell(v)&=v^\ell \sigma_0(\ell)\exp(mH(v,\ell/m)),\\
 \sigma_0(\ell)&=\exp\!\left(Dd[(m+\ell)\log(1+\ell/m)-\ell]\right).
 \end{split}                                                     \tag{6}
\]

The exact eligible local contribution from \(i\ge1\)
to the normalized moments now takes the form
\[
                         U_i^m F_\ell(v_i) A_i(m+\ell).
 \tag{7}
\]
Here \(A_i(t)=(tD)^{-\rho_i} \mathcal B_i(t)\),
\(\rho_i=r-\# S_i/2\), with
\[
 \mathcal B_i(t)=K_i(t)+O((mD)^{-1/2})
                \quad\text{on } |t-m|<a_2 m/D.
\]
The prefactor \(K_i(t)\) differs from the positive pure
amplitude by \(O(e^{-cD})\). These are the complex
translated-contour bounds of the
[homogeneous sector estimate](SFC_HOMOGENEOUS_DILATED_SADDLE_WINDOWS.md#a-uniform-analytic-saddle-amplitude),
applied to the slightly varying saddles (2). As in the
[singleton phase estimate](SFC_MIXED_DEGREE_SINGLETON_FOREIGN_BLOCKS.md#the-shifted-expansion),
the losses on the translated box and its side contours are
uniform when \(|t-m|<a_2 m/D\). Spatial logarithms are
holomorphic on the relative boxes and \(t\) dependence can be
taken on the undeformed integral. Cauchy estimates for
any fixed \(h\ge1\) make the derivatives of the peeled
amplitude, after multiplication by \(m^h\), of size
\(O(D^h((mD)^{-1/2}+e^{-cD}))=o(1)\). Thus
\[
 \frac{m^h\Delta^h A_i(m)}{(mD)^{-\rho_i}}
       =(-1)^h(\rho_i)^{\overline h}K_i^*+o(1),\qquad K_i^*>0 ,
 \tag{8}
\]
for fixed \(h\ge0\), taking a subsequential limit if the
columns vary. Here \(\Delta\) denotes forward difference
and the fixed-shift argument integrates the \(h\)-th
derivative over a unit cube.

For an eligible constant peak, (7) is replaced by
\(U_0^m v_0^\ell A_0(m+\ell)\) with \(U_0=c_0/Q\),
\(v_0=U_0 e^{-Dd}\), and pure endpoint amplitude
\(A_0(t)=(tD)^{-r}(1+O(e^{-ctD}))\).
Every nonconstant monomial vanishes to high order in at
least one inactive coordinate at the constant peak.
The local comparison of the
[incomparable-coordinate theorem](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates)
therefore bounds the remaining term by
\(C_L(mD)^{-L}|U_0|^m|v_0|^\ell(mD)^{-r}\)
for any fixed \(L\). It also gives the usual sector
bounds for the pure endpoint amplitude. All sequences in
(7) and the constant mode are normalized with the same
nonzero scalar at each shift.


## Interpolation of analytic modes on two annuli

The effective frequencies from (4) lie in
\(\frac14 e^{-\gamma D}\le |v_i|\le 2\).
If eligible, the constant node lies in
\(e^{-(d+\gamma)D}\le |v_0|\le e^{-Dd}\).
Write \(a_D=\frac1{16}e^{-\gamma D}\).
Choose \(\gamma\) sufficiently small after \(\eta\) is
fixed. At the finitely many shifts needed for any filter,
(5) and (6) give
\[
 \sup_{a_D\le |v|\le 4}
     |F_\ell(v)-v^\ell|
        \le C D/m + C D^2 e^{-cD}/m .                         \tag{9}
\]
The inverse defining \(H\) exists on this annulus by the
same contraction. Indeed \(H\) has a double zero in \(T\)
and its Cauchy bound on a disc of radius proportional to
\(1/D\) gives the second term in (9).

Consider a list of at most \(s'-1\) prescribed roots,
with multiplicities, drawn from the eligible positive
nodes and at most one copy of the constant node. Let
\(n\) be the total number, and choose \(p_n=1\), seeking
the other coefficients so that \(\sum_{\ell=0}^n
p_\ell F_\ell(v)\) has the prescribed positive roots, and
\(\sum_{\ell=0}^n p_\ell v^\ell\) vanishes at the
constant node if prescribed. Denote these two test functions
by \(P_H(v)\) and \(P_0(v)\).
Compare the determinant of the interpolation system to
that of the ordinary polynomial root system. Divide by
the positive-node Vandermonde, with confluent interpretation.
The divided matrix has entries that are divided differences
of the \(F_\ell\); the constant row is unchanged.
The divided difference at \(j+1\) positive nodes is
\[
 \frac{1}{2\pi i}\int_{\partial\{a_D<|\zeta|<4\}}
     \frac{F_\ell(\zeta)}
           {\prod_{h=0}^{j}(\zeta-w_h)}\,d\zeta,
\]
where \(w_h\) are those nodes. At distinct nodes, replacing
rows successively by Newton divided differences divides the
determinant by the usual Vandermonde. The contour formula
extends each row holomorphically as nodes repeat. At repeated
nodes its residues are the corresponding confluent data,
so the divided systems describe the prescribed derivative
conditions as well. The distance from the nodes to either
circle is at least \(a_D\), so the error in
each entry after replacing \(F_\ell\) by \(\zeta^\ell\) is
bounded by the right side of (9) times \(C a_D^{-s'}\).
The polynomial divided differences stay bounded because
the roots are bounded. After division, the ordinary
determinant equals one up to sign if the constant was
absent, and otherwise equals, up to sign,
\(\prod_{j\in H}(w_j-v_0)\) of modulus at least
\(c a_D^{s'}\). By (9) and \(m>e^{\eta D}\), for small
enough \(\gamma\) the perturbation relative to this
determinant tends exponentially to zero. Thus the system
is invertible. The same estimates on the Cramer minors
give bounded coefficients, since the ordinary monic
polynomial of these roots has bounded coefficients.

The augmented determinant with an evaluation row gives
\[
 \begin{split}
  P_H(v)&=(1+o(1))\prod_{j\ {\rm prescribed}}(v-w_j),
      && \tfrac18 e^{-\gamma D}\le |v|\le 3,\\
  P_0(v)&=(1+o(1))\prod_{j\ {\rm prescribed}}(v-w_j),
      && |v-v_0|\le |v_0|/4\quad\text{if eligible}.
 \end{split}                                                     \tag{10}
\]
The relative errors are holomorphic and uniform in these
sets, including roots after factoring. For the first line,
include the evaluation point among the positive nodes when
forming divided differences, which extracts its differences
with those nodes; its difference from the constant node is
separated by \(c a_D\). For the second line, if the constant
was prescribed, subtract its row and divide by \(v-v_0\)
before estimating the augmented determinant. These quotients
of ordinary powers are bounded on the constant disc, and
the remaining ordinary products with positive nodes are
separated by \(c a_D\). This proves the uniform relative
bounds and gives Cauchy bounds for those holomorphic factors
on discs of radius a fixed multiple of the node modulus.
The case \(n=0\) uses empty products.


## From shift powers to analytic derivatives

Testing the amplitude involves powers of the shift \(\ell\).
These can be transferred to scaled derivatives of the analytic
test function. Write \(\mathcal X=(v\partial_v)/m\) and
\[
     L_h(v,T)=\mathcal X^hF_\ell(v)/F_\ell(v),\qquad T=\ell/m.
\]
These extend holomorphically in \(T\) through the recurrence
\[
 L_0=1,\qquad L_{h+1}=(T+v\partial_v H)L_h+\mathcal X L_h .    \tag{11}
\]
For any fixed \(J\), the Taylor coefficient matrix of
\(L_0,\ldots,L_J\) through degree \(J\) in \(T\) equals
the identity plus \(o(1)\) uniformly on a fixed radial
enlargement of the eligible nodes. Indeed
\(H=O(e^{-cD})\) on \(|T|<a_1/D\) and its fixed
logarithmic \(v\) derivatives have the same bound on
nested enlarged annuli. The recurrence makes
\(L_h-T^h=O(e^{-c'D})\) on a smaller product domain for
each fixed \(h\). Its Taylor coefficients and the
Taylor remainder at \(T=\ell/m\) are therefore bounded
respectively by \(O(D^J e^{-c'D})\) and
\(O(D^{J+1}e^{-c'D}m^{-J-1})\). Inverting the finite
Taylor matrix yields
\[
  (\ell/m)^k F_\ell(v)
    =\sum_{h=0}^J(\mathbf1_{h=k}+o(1))\mathcal X^hF_\ell(v)
       +O(m^{-J-1}|F_\ell(v)|),\qquad 0\le k\le J.             \tag{12}
\]
The coefficients of the inversion are independent of the
shift \(\ell\).

Apply \(p_\ell\) to the contribution (7). For
\(\ell\le n\), use the exact Newton expansion
\(A_i(m+\ell)=\sum_{h=0}^n\binom{\ell}h\Delta^h A_i(m)\).
After division by \((mD)^{-\rho_i}\), the coefficients of
this polynomial in \(T=\ell/m\) tend by (8) to those of
\(\sum_{h=0}^n K_i^*(-1)^h(\rho_i)^{\overline h}T^h/h!\).
By (12), the filtered contribution is
\[
 U_i^m(mD)^{-\rho_i}
       \left[\sum_{h=0}^J(\lambda_{ih}+o(1))
                             \mathcal X^hP_H(v_i)
                   +O(m^{-J-1})\right],                       \tag{13}
\]
where \(\lambda_{ih}=K_i^*(-1)^h(\rho_i)^{\overline h}/h!\)
for \(h\le n\) and is zero for \(h>n\).
The bounds on \(p_\ell\) and \(F_\ell(v_i)\) justify the
error. The transfer uses finite Taylor inversion and
keeps the exact root conditions in (10).


## Contradiction in a vanishing window

Over eligible nodes put
\[
 d_{ij}=|v_i-v_j|+|v_i|/m,\qquad
 V_i=\prod_{j\ne i}d_{ij},\qquad
 W_i=|U_i|^m(mD)^{-\rho_i}V_i ,
\]
including \(\rho_0=r\) if the constant is eligible.
Choose a maximum at \(p\) and take a subsequence on which
\(p\) and the relative cluster \(I\) are fixed:
\(i\in I\) if \(m(v_i/v_p-1)\to x_i\) is finite, and the
modulus tends to infinity outside \(I\). The cluster is
within the positive block or consists of the constant
alone. Prescribe all outside eligible nodes and \(k\)
copies of \(v_p\), for \(0\le k<|I|\), in (10). The degree
is at most \(s'-1\) and the tested vanishing window has
enough terms. Write \(\Pi(v)=P_H(v)/(v-v_p)^k\) if the
cluster is positive, and \(\Pi=P_0\) for a constant
cluster, where \(k=0\). Its test normalization is
\[
              Z_k=U_p^m(mD)^{-\rho_p}\Pi(v_p)(v_p/m)^k .
\]
By (10), its modulus is comparable along the subsequence
to \(W_p(|v_p|/m)^{k-|I|+1}\).

The error of (13) divided by \(Z_k\) tends to zero by
\(W_i\le W_p\) and
\(V_i\ge e^{-C_s D}m^{-(s'-1)}\), choosing \(J>s'+C_s/\eta\).
The local constant error is also negligible with \(L\)
large enough. For the exterior error, a raw maximizing
eligible label has \(|U_i|\ge e^{-C e^{-cD}}\); this follows
from its saddle bound or from \(|U_0|=1\). Hence
\(W_p\ge c(mD)^{-r} m^{-(s'-1)}e^{-C_s D-Cm e^{-cD}}\),
which dominates \(O(e^{-cmD})\). The extra factor
\((|v_p|/m)^{k-|I|+1}\) in \(Z_k\) is bounded below.

For an outside positive node, (10) has the exact factor
\(v-v_i\). In \(\mathcal X^h P_H(v_i)\) some derivative
must remove it. Expanding
\[
 \mathcal X^h=m^{-h}\sum_{a=0}^h
                 S(h,a)v^a\partial_v^a
\]
with the integer coefficients \(S(h,a)\) of iterated
Euler differentiation, charge a scale \(|v_i|/m\) to
each differentiated factor and use the padded distances
for the others. Derivatives of the holomorphic relative
factor from (10) cost \(O(|v_i|^{-1})\) per derivative by
Cauchy's formula. This bounds (13), before its negligible
error, by
\[
              C W_i (|v_i|/m)d_{ip}^{\,k}
                                     \big/\prod_{j\in I}d_{ij}.
\]
If \(i\) is the outside constant, the pure Newton
difference bound and the same factor calculation for
\(P_0\) give this estimate. After normalization all these
terms tend to zero, since \(d_{ij}\sim d_{ip}\) for
\(j\in I\) and the ratios
\((|v_i|/m)/d_{ip}\) and \((|v_p|/m)/d_{ip}\) vanish.

If the cluster is positive, formula (4) gives
\(U(v)/v=\exp(Dd B g'(B))\), with logarithmic derivative
in \(v\) of size \(O(D e^{-cD})\). Thus
\((U_i/U_p)^m\to e^{x_i}\) for \(i\in I\).
The padded products \(V_i/V_p\) stay bounded above and
below. By maximality, \(\rho_i\ge\rho_p\) in \(I\),
and the terms with strict inequality vanish after
normalization. Within \(I\), (10) gives
\(\Pi(v_i)/\Pi(v_p)\to1\); every positive derivative
of \(\Pi\), scaled by the corresponding power of
\(v_p/m\), tends to zero relative to \(\Pi(v_p)\).
This follows from the outside separations and the
holomorphic factor bounds. Consequently the limit of
\(\mathcal X^h P_H(v_i)/(\Pi(v_p)(v_p/m)^k)\)
is \(k^{\underline h}x_i^{k-h}\) for \(h\le k\) and
zero for \(h>k\). Here \(k^{\underline h}\) is the
falling factorial. Equation (13) gives the limit
equations
\[
  \sum_{\substack{i\in I\\\rho_i=\rho_p}}
      K_i^* e^{x_i}\sum_{h=0}^k
       \binom kh x_i^{k-h}(-1)^h(\rho_p)^{\overline h}
     =0,\qquad 0\le k<|I|.
\]
The inner polynomials are monic of consecutive degrees.
Grouping equal \(x_i\) yields the Vandermonde contradiction,
since each grouped weight is \(e^{x_i}\) times a positive
sum. If \(p\) is the constant alone, \(k=0\), and its
positive prefactor survives whereas its positive differences
vanish in the same normalization. That also contradicts
vanishing of the window.

This proves the claimed existence of a uniform threshold.
The support-size window convention is
[Edo and van den Essen, Conjecture 2.8](https://arxiv.org/abs/1304.3956).
A coefficient field of characteristic zero is handled by
embedding the finitely generated field of coefficients
into \(\mathbb C\).




