# Deepest rational order centers for factorial windows

Let \(F=cG\) be a nonmonomial polynomial over a finite extension of
\(\mathbb Q_p\), with \(v(p)=1\), and write
\[
 G=X^a+\sum_{\nu\ne a}t_\nu X^\nu,\qquad
 \delta=\min_{\nu\ne a}\left(v(t_\nu)+\frac{|\nu|-|a|}{p-1}\right)>0.
\]
Put \(U(n)=\Lambda(F^n)/(c^n\prod_i(a_i n)!)\). The dimension, degree
and finite support are arbitrary. The width \(w\) is any fixed positive
integer. In particular one can take \(w\) equal to the support size.

For every start \(N\), an effective search among rational numbers of
bounded denominator and numerator \(O(1+\log_p(N+w))\) finds a depth
\(L\). At every requested precision \(P\ge0\), all \(w\) normalized
moments are given modulo valuation \(P\) by counts of size
\(O(1+P+L)\). When \(L\) is large the search also gives a unique common
center \(\rho\), of height \(O(1+L)\), for the entire window. For every
fixed \(J\ge0\), with \(z=N-\rho\) and \(v(z)=L\), one has
\[
 v\left(z^{m(\rho+j)}U(N+j)-
           \sum_{b=0}^J B_b(\rho+j)z^b\right)
 \ \ge\ (J+1)L-C_J(1+\log_p(L+2)),\qquad 0\le j<w .                 \tag{1}
\]
The \(B_b\) are the regularized Taylor coefficients of the same
count terms as in the
[logarithmic tower expansion](SFC_LOGARITHMIC_TOWER_RESIDUES.md).
Here the centers can move and there is no size condition relating
\(N\) to \(p^L\). The depth must be chosen by the rule below.
An exact zero window at large \(L\) imposes simultaneous
divisibility on the jets at this one center, including on their
minors. Bounded \(L\) admits finite congruence certificates.
These conclusions use only the isolated effective slope.

The arithmetic input is the same factorial counting behind
[Granville's Section 2](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node3-an.shtml).
The fixed-center tower proof discards count tails using an error
depending on the size of the order. Here the tail is instead
controlled by the largest distance valuation to a rational root
whose height is small relative to that valuation.

## 1. A count bound uniform in the order

For counts \(k_\nu\ge0\), put \(K=\sum_\nu k_\nu\) and
\(\Delta_i=\sum_\nu k_\nu(\nu_i-a_i)\). Write
\[
 R_k(T)=\frac{(T)_{\underline K}}{\prod_\nu k_\nu!}
             \prod_\nu t_\nu^{k_\nu}
             \prod_i\phi_{\Delta_i}(a_iT),                          \tag{2}
\]
where \(\phi_b(T)\) is the rational factorial ratio \((T+b)!/T!\),
interpreted as a product or reciprocal product. Then
\(U(n)=\sum_{K\le n}R_k(n)\). Put \(A=\max(1,\max_i a_i)\).
A denominator root in row \(i\) is \(b/a_i\) with
\(0\le b<-\Delta_i\le a_iK\). Rows with \(a_i=0\) never lose
factors. For the window define the nested finite grids
\[
 \mathcal P_K=\{b/a_i-j\in\mathbb Z_p:
       a_i>0,\ 0\le b<a_iK,\ 0\le j<w\},\qquad K\ge1 .             \tag{3}
\]
Each \(\rho\in\bigcup_K\mathcal P_K\) has denominator at most
\(A\). Let \(k(\rho)\) be the least \(K\) with \(\rho\in\mathcal P_K\).
The grid has \(O(K)\) points, and any two distinct points of
\(\mathcal P_K\) differ by a rational with
\(p\)-adic valuation at most \(\log_p(2A^2(K+w+1))\).
Indeed after multiplying the difference by the two denominators,
its nonzero integer numerator has absolute value at most
\(2A^2(K+w+1)\). Roots outside \(\mathbb Z_p\) cannot be \(p\)-adically
close to an integral order.

For any rational \(\alpha\in\mathbb Z_p\), define \(m(\alpha)=0\)
when \(\alpha<0\). Otherwise let
\(I_\alpha=\{i:a_i>0,\ a_i\alpha\in\mathbb Z\}\). Let \(h_\alpha\)
be the maximum number of coordinates in \(I_\alpha\) strictly
decreased from \(a\) by a convex combination of the other support
columns. Put
\[
 m(\alpha)=\max(0,h_\alpha-\boldsymbol 1_{\alpha\in\mathbb Z})
 \quad(\alpha\ge0),\qquad
 M=\max_{\alpha\in\mathbb Q\cap\mathbb Z_p} m(\alpha).
\]
The maximum for \(M\) is a finite support calculation: \(I_\alpha\)
depends only on a reduced denominator dividing some \(a_i\), and
nonintegral centers with empty \(I_\alpha\) contribute zero.
At \(\alpha\) each row can contribute at most one linear denominator
factor with root \(\alpha\). If \(\alpha\) is a nonnegative integer
and there is a pole row, then \(\alpha<K\) and the falling factorial
cancels one such factor. Thus \(R_k\) has pole order at most
\(m(\alpha)\) at \(\alpha\).

For \(1\le K\le n=N+j\), the following estimate holds uniformly:
\[
 v(R_k(n))\ge \delta K-C(1+\log_p(K+w+1))
                 -M\max_{\substack{\rho\in\mathcal P_K\\ \rho<N}}
                                 v(N-\rho).                         \tag{4}
\]
The maximum over an empty set is zero. For the proof, take
\(\lambda_K=1+\lceil\log_p(2A^2(K+w+1))\rceil\) and
\(b_*=\max_{a_i>0}v_p(a_i)\), with empty maximum zero.
Count divisibility up to level \(\lambda_K+b_*\) in each interval
of linear factors in (2). An interval of length \(q\), with
consecutive integer entries, has truncated valuation
\(q/(p-1)+O(1+\log_p(K+w+1))\). In fact at each level the number
of multiples differs from \(q/p^\ell\) by at most one, and all
the interval lengths here are \(O(K)\). This argument also works
at \(p\)-integral order arguments where the interval is a translate
of consecutive integers in \(\mathbb Z_p\).
Legendre's formula for the count factorials gives the same
precision. The linear parts in \(K\) from the falling factorial and
the count factorials cancel, leaving at least \(\delta K\),
with the indicated logarithmic error.

Any denominator factor above the truncation level has root in
\(\mathbb Z_p\) and distance valuation greater than \(\lambda_K\).
It belongs to a single root across the window, by (3). If \(q\)
denominator factors have that root in a term, the excess cost is
at most \(q v(N-\rho)\). If the order center \(\rho+j\) is an
integer, the falling factorial supplies the same valuation above
the truncation level, up to the truncation level itself. This loss
is then at most \((q-1)v(N-\rho)+O(1+\log_p(K+w+1))\).
The resulting coefficient of \(v(N-\rho)\) is at most \(M\).
Every denominator root of an actual term satisfies
\(\rho+j<K\le N+j\), hence \(\rho<N\). This proves (4).
More precisely one can use the distance to the single root
above the truncation level, with zero loss if none occurs.
The proof also gives (4) for \(R_k(x+j)\) at \(x\in\mathbb Z_p\)
outside the grids, for every \(K\), if the maximum is taken over
all of \(\mathcal P_K\). The same cancellation in the falling
factorial applies at integer roots irrespective of \(x\).

## 2. The deepest center and the finite residue sum

Fix \(\theta>0\) with \(M\theta<\delta\), allowing any positive
\(\theta\) when \(M=0\). Given \(N\), call a center eligible if
\[
 \rho<N,\qquad \rho\in\bigcup_K\mathcal P_K,\qquad
                  v(N-\rho)>\theta k(\rho).                          \tag{5}
\]
Let \(L\) be the greatest distance valuation among eligible centers,
with \(L=0\) when none exists. Since every such \(\rho\) lies in
\([-(w-1),N)\), multiplying \(N-\rho\) by its denominator gives
a nonzero integer of size at most \(A(N+w)\). Thus
\(L\le \log_p(A(N+w))\) and eligibility can be decided using (3)
at \(K\le\lceil\log_p(A(N+w))/\theta\rceil\).
A center attaining \(L\) satisfies \(k(\rho)<L/\theta\), so its
absolute value is \(O(1+L)\). If two centers attain it, their
difference has valuation at least \(L\), whereas by (3) its
valuation is \(O(1+\log_p(L+2))\). The maximizing center is
therefore unique for all sufficiently large \(L\).

For every \(K\), the maximum in (4) is at most
\(\max(L,\theta K)\). It follows that
\[
 v(R_k(N+j))\ge (\delta-M\theta)K
                    -M L-C(1+\log_p(K+w+1)).                         \tag{6}
\]
For \(K>C'(1+L+P)\), the right side is at least \(P\), with an
effective constant \(C'\) independent of \(N,j,P,L\).
This proves the finite sum assertion. In (6) the actual count
restriction remains \(K\le N+j\).

In particular, fix a depth bound \(D\) and a precision \(P>0\).
For \(L\le D\), use a fixed cutoff \(K_0\) for that assertion.
Every denominator root in the cutoff has distance valuation
at most \(\max(D,\theta K_0)\) for \(N>K_0\); roots outside
\(\mathbb Z_p\) have constant factor valuation. On the complement
of the corresponding open balls around \(\mathcal P_{K_0}\) in
\(\mathbb Z_p\), the finite sum in (2), for each \(j\), is
continuous and its residues at valuation \(P\) are decided by
finitely many congruence classes. One can obtain the modulus
directly by bounding the denominator valuations in the finite
products and testing their numerators at that precision.
All starts \(N>K_0\) of depth at most \(D\) belong to this
complement. If no congruence class there has all its finite
sums of valuation at least \(P\), there is no exact zero window
among those starts. Starts \(N\le K_0\) are finite and can be
tested directly with their actual count bounds.

## 3. Uniform Taylor expansion at the selected center

Let \(L\) be large and \(\alpha=\rho+j\), \(0\le j<w\).
For each \(b\ge0\) define
\[
 B_b(\alpha)=\sum_{k}[Z^b]\bigl(Z^{m(\alpha)}R_k(\alpha+Z)\bigr).
                                                                         \tag{7}
\]
For fixed \(b\) the valuation of each coefficient in the sum
is at least
\[
         \delta K-C_b(1+\log_p(K+|\alpha|+2)).                          \tag{8}
\]
To check uniformity in \(\alpha\), write \(\alpha=u/q\) in lowest
terms with \(q\le A\) and \(p\nmid q\). After multiplying by \(q\),
the offsets of the factors in (2) at \(Z=0\) form intervals in
arithmetic progressions with step \(q\), whose entries have
absolute values \(O(K+|\alpha|+1)\). Removing zero changes the
count at each precision by at most one. The interval count
therefore controls the valuation of the first nonzero coefficient
with error \(O(1+\log_p(K+|\alpha|+2))\); the stripped zero
factors contribute only fixed integer slopes. A positive row
with a zero increases the power of \(Z\).
Removing the first coefficient and its power of \(Z\), all
remaining factors are \(1+dZ/e\) or inverses of this expression,
where \(d\) is integral and \(e\ne0\) is one of the offsets.
The coefficient of degree at most \(b\) loses at most \(b\)
times the largest offset valuation. This proves (8), convergence
of (7), and
\[
                   v(B_b(\alpha))\ge-C'_b(1+\log_p(L+2)).
                                                                         \tag{9}
\]
To deduce (9), use \(|\alpha|=O(1+L)\) and bound the infimum
over \(K\) after separating \(\log_p(K+|\alpha|+2)\) into a
sum of \(O(1+\log_p(L+2))\) and \(O(\log_p(K+2))\).

Take a cutoff \(K_*= \lceil C''_J(L+1)\rceil\). With \(C''_J\)
sufficiently large, (6) discards every actual term above the
cutoff to valuation \((J+1)L\). The same cutoff works for
degrees at most \(J\) of the regularized formal tail by (8).
For \(K\le K_*\) and \(L\) large, every nonzero offset in the
product at \(\alpha\) has valuation \(O(1+\log_p(L+2))<L\).
The remainder of the convergent expansion at \(Z=z\) after
degree \(J\) has valuation at least
\[
     (J+1)L+\delta K-C'''_J(1+\log_p(K+|\alpha|+2)).
\]
This follows from the geometric factors used in (8): each
additional degree beyond \(J+1\) adds at least \(L\) minus
the largest offset valuation. If the initial power of \(Z\) already exceeds \(J\),
use the same bound directly. Finally, since \(q(N-\rho)\)
is a nonzero integer of valuation \(L\), its absolute value
is at least \(p^L\). Thus \(N\ge p^L/A-|\rho|\), and all
counts up to \(K_*\) are actual counts at every offset for
large \(L\). Combining the estimates proves (1).

At nonnegative integer centers (7) agrees with the
logarithmic functional (7) of the fixed-center tower proof.
At positive noninteger centers it uses the regularization of
Section 4 there. Negative centers have no denominator pole
and use (7) directly.

## 4. Conditions on an exact zero window

Let \(0\le J<w\), and let \(j_0,\ldots,j_J\) be distinct offsets
in the window. If all the moments vanish at a start of large
depth, (1) implies
\[
  v\left(\det\big(B_b(\rho+j_i)\big)_{0\le i,b\le J}\right)
                  \ge (J+1)L-C^{*}_J(1+\log_p(L+2)).                 \tag{10}
\]
Indeed replace the first column by its sum with the other
columns weighted by \(z^b\), which leaves the determinant
unchanged. The new column has the lower bound from (1).
All entries in the other columns satisfy (9), so expansion
along the new column proves (10).
For \(J=1\) the necessary condition uses the cross-difference
of the leading pole coefficient and the regular coefficient
at two consecutive centers. The estimate uses those
coefficients jointly, without requiring either leading
coefficient to be nonzero in isolation. If a leading column
vanishes identically, (1) can instead be divided by its
known vanishing power of \(z\) before applying the determinant
argument.

A sufficient local route to a full width-\(w\) window theorem
for an isolated-slope phase is therefore to find fixed
\(J<w\) and a uniform upper bound on some minor in (10)
strictly below its lower bound at every eligible large-depth
center, then certify the bounded-depth range by Section 2.
The individual-pole argument alone supplies no upper bound
on those minors. The window problem here calls for a joint
estimate with the center height growing linearly in the
distance valuation, rather than just nonflatness at each
fixed center.

## 5. Fixed centers and a size threshold for possible zeros

Fix \(\sigma\in\mathbb Q\cap\mathbb Z_p\), and let
\(b=v(N-\sigma)\) tend to infinity along positive integer
starts. A sufficiently small constant \(c>0\), depending on the
fixed data and \(\sigma\), makes
\[
                    \log_p(N+w)\le c p^b                              \tag{11}
\]
sufficient for the fixed-center expansion through all the offsets,
with error valuation \((J+1)b-O_{\sigma,J}(1)\).
To see this, every eligible center has \(k(\rho)=O(1+\log_p(N+w))\).
If \(\rho\ne\sigma\) and \(v(N-\rho)\ge b\), the difference
\(\rho-\sigma\) has valuation at least \(b\). Its numerator after
multiplying by the denominators has size \(O(1+\log_p(N+w))\)
and is nonzero. Condition (11) with small enough \(c\) makes
that size strictly less than \(p^b\), which is impossible.
Hence \(L\le b\), whether or not \(\sigma\) itself occurs in
the grids. Use (6) to truncate at \(O_{\sigma,J}(b+1)\) and
then repeat Section 3 at the fixed \(\sigma\). The Taylor
remainder has valuation at least \((J+1)b+\delta K-
C_{\sigma,J}(1+\log_p(K+2))\), whose last two terms have a
fixed lower bound. For each fixed \(\sigma\) the inequality
\(N\ge p^b/\operatorname{den}(\sigma)-|\sigma|\) ensures that
all low counts are actual counts for large \(b\).
Here (7) also applies when \(\operatorname{den}(\sigma)>A\):
one uses that fixed denominator in the offset estimate (8).
In particular (11) includes \(\log_p(N+w)=o(p^b)\).

If some germ at \(\sigma+j\) has a nonzero coefficient, the
first nonzero coefficient dominates at large \(b\) under (11).
Thus exact zero windows approaching that fixed center must satisfy
\(\log_p(N+w)>c p^b\) for all sufficiently large \(b\).
Equivalently their distance valuation from \(\sigma\) is bounded
by \(\log_p\log_p(N+w)+O(1)\). The corresponding fixed-center
nonflatness obligation alone gives no estimate at the moving
centers in Section 4.

## 6. Joint residue frequencies at arbitrary starts

The same estimate gives a limiting distribution of the normalized
windows at every fixed \(p\)-adic precision. For \(T\ge1\) put
\[
 \mathcal E_T=\bigcup_{\rho\in\cup_K\mathcal P_K}
     \{x\in\mathbb Z_p:v(x-\rho)>\theta\max(T,k(\rho))\}.
\]
There are \(O(1)\) new grid centers at each \(K\) and \(O(T)\)
centers up to \(T\). The normalized Haar measure of
\(\mathcal E_T\) is therefore \(O((T+1)p^{-\theta T})\).
Outside it, all \(R_k(x+j)\) for \(K\ge T\) have valuation at
least \((\delta-M\theta)K-C(1+\log_p(K+w+1))\), provided \(x\)
is outside the countable grid. Thus
\[
                   \Phi(x+j)=\sum_k R_k(x+j),\qquad 0\le j<w ,
\]
is defined almost everywhere and has uniform count tails outside
\(\mathcal E_T\) at the stated precision.

For the actual integer windows with \(1\le N\le X\), the same
count-tail conclusion holds outside the balls in \(\mathcal E_T\)
with centers \(\rho<N\). Only centers with
\(k(\rho)\le\log_p(A(X+w))/\theta\) can affect these counts,
since the numerator of the nonzero difference is at most
\(A(X+w)\). Each of these balls is one residue class of some
integer prime-power modulus. Counting integers in those classes
gives the upper bound
\[
 O\bigl(X(T+1)p^{-\theta T}\bigr)+O(1+\log_p(X+w))
\]
on the number of exceptional actual starts. Allowing the centers
\(\rho=N\) would spoil this argument, because such centers do
not occur in the actual denominator factors.

Fix \(P\ge0\). For every set \(S\) of residue vectors in the
quotient by \(\{u:v(u)\ge P\}\) in each component, the natural
frequency of actual windows with residue vector in \(S\)
exists. It equals the Haar measure of the corresponding
\(\Phi\)-window event. To prove this, first use the finite sum
through \(K<T\). Its poles in \(\mathbb Z_p\) are contained in
the finite grid \(\mathcal P_T\), and one can delete balls of
arbitrarily small total measure around those finitely many
points. On the compact complement, reduction of the finite sum
at precision \(P\) is constant on a finite collection of
congruence classes, so its integer frequency equals its Haar
measure there. Now take \(T\) large enough that the tail bound is
at least \(P\). The discrepancy between the finite-sum
frequency and the actual frequency is bounded in upper density
by the estimate for \(\mathcal E_T\); a fixed finite initial
segment due to \(N<T\) is immaterial. The Haar discrepancy
between the finite and limiting sums has the same bound.
Sending \(T\) to infinity proves the frequency assertion.

Consequently the upper density of exact zero windows is at
most the Haar measure of the common zero set of the \(\Phi\)
window, by taking \(S=\{0\}\) and letting \(P\) grow. For a
simple unconditional bound at a strongly isolated slope take
\(\theta=\delta/(2\max(1,M))\). If \(\delta\) is large enough
that \(\delta K/2>C(1+\log_p(K+w+1))\) for every \(K\ge1\),
then outside the integer exceptional balls at \(T=1\) every
nonanchor term has strictly positive valuation and every
normalized moment in the window is a unit. For \(\theta\ge1\)
the exceptional upper density is \(O(p^{-\theta})\) with
constant depending only on the support, the width and the
prime. The remaining exact-window question concerns starts inside
those exceptional congruence classes.

## 7. A sharper selection using the support and the slope costs

The threshold can use the actual slope cost of a pole. Write
\(\delta_\nu=v(t_\nu)+(|\nu|-|a|)/(p-1)\) and
\(D(k)=\sum_\nu k_\nu\delta_\nu\).
For a root \(\alpha\ge0\) let \(q(k,\alpha)\) be the pole order
of \(R_k\) after canceling the falling factorial, taking zero
when it is regular. At a nonnegative root the numerator factors
from increasing rows have no zero, so this order is exactly
the number of pole rows minus the falling-factorial cancellation,
with a minimum of zero. Define
\[
   \mu(\rho)=\min_{\substack{0\le j<w,\ k\\q(k,\rho+j)>0}}
                            \frac{D(k)}{q(k,\rho+j)},                 \tag{12}
\]
with an empty minimum \(+\infty\). Each finite minimum is
attained since \(D(k)\ge\delta K\) and \(q(k,\rho+j)\le M\).
In particular \(\mu(\rho)\ge \delta k(\rho)/M\) when \(M>0\).

For \(0<\epsilon<1\) replace (5) by
\[
       \rho<N,\qquad v(N-\rho)>(1-\epsilon)\mu(\rho).
\]
If \(M=0\), all the terms are regular on the \(p\)-integral
orders and the weighted rule selects no centers. Otherwise the
search still examines \(O(1+\log_p(N+w))\) grid centers. In
deciding eligibility only count vectors with
\(D(k)\le M\log_p(A(N+w))/(1-\epsilon)\) are needed. The
maximizing center still has height \(O(1+L)\) and is unique
at large depth by the same rational separation.

Apply the more precise single-root form of (4). A term with a
deep root \(\rho\) and multiplicity \(q=q(k,\rho+j)>0\) loses at
most \(qv(N-\rho)\) beyond the logarithmic error. If the root
fails the weighted threshold, that loss is at most
\((1-\epsilon)q\mu(\rho)\le (1-\epsilon)D(k)\). If it passes,
the loss is at most \(ML\). Thus the weighted version of (6)
has the lower bound
\[
              \epsilon D(k)-ML-C(1+\log_p(K+w+1)).
\]
All conclusions of Sections 2–4 now follow with this bound.
For the bounded-depth finite sums first cancel the removable
roots of each head term. Each remaining \(p\)-integral root
in the head has finite \(\mu(\rho)\) and its valuation distance
is bounded by \(\max(D,(1-\epsilon)\mu(\rho))\) when \(L\le D\).
This provides the finite congruence modulus exactly as before.

For example, when \(w\ge2\) a sufficient large-depth result is
an upper bound of the form
\[
 \min_{0\le i<j<w}
 v\!\left(B_0(\rho+i)B_1(\rho+j)-B_0(\rho+j)B_1(\rho+i)\right)
 \ \le\ (2-\eta)\mu(\rho)+C_*(1+\log_p(\mu(\rho)+2))                  \tag{13}
\]
at all centers with finite cost, with fixed \(0<\eta<2\).
Choose \(\epsilon<\eta/2\). At an eligible maximizing center
one has \(\mu(\rho)<L/(1-\epsilon)\), so (13) contradicts
(10) for all sufficiently large \(L\). The same calculation
works for minors of size \(J+1\), with \(2\) replaced by
\(J+1\). This reduces the high-depth proof to a quantitative
upper bound referenced to the first possible pole cost.
In two variables at integer centers the genuine pole order is
one, so (12) there is the least slope cost of counts descending
past the same integer center in both rows of some offset.
A uniform bound such as (13), or a replacement using higher
layers when its minors vanish, remains to be proved.
