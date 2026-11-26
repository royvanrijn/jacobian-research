# Common slope removal in retained factorial rows

In an ordinary two-monomial profile, strip the factorial of each
row's smaller exponent from each moment. Equal rows then contribute
one. All comparisons below are in the physical coefficient ratio.
The retained rows may have arbitrary independent heights, gaps,
orientations, and multiplicities.

## Arbitrary retained rows at one prime

Fix consecutive orders \(n,N=n+1\). Write \(s_{j,1},\ldots,s_{j,j}\)
for the ordered unit slopes of the lower Newton polygon of the
selected rows at prime \(p\). The selected coefficient arrays include
the binomial coefficients. Each slope is listed according to width.
For a retained row \(i\), write \(\ell_i=\min(a_i,b_i)\),
\(M_i=\max(a_i,b_i)\), \(d_i=M_i-\ell_i\), and
\(\sigma_i=\operatorname{sgn}(b_i-a_i)\). Define
\[
 L_i=\operatorname{lcm}\{j\ell_i+1,\ldots,jM_i:\ j\in\{n,N\}\},\quad
 R_{i,p}=v_p(L_i),\quad R_p=\sum_{i\ {\rm retained}}R_{i,p}.
\]
An empty least common multiple is one.

**Lemma.** The full polygon slopes differ in order from
\(s_{j,k}+\beta_p\) by at most \(R_p\), where
\[
 \beta_p=\sum_{i\ {\rm retained}}\sigma_i d_i
                                 \sum_{e=1}^{R_{i,p}}p^{-e}
\]
is common to both orders. Join a lower-order and an upper-order
selected unit slope when their absolute difference is at most
\(2R_p\). The size of a maximum matching bounds the degree of a
common rational factor. In particular a bound of one gives
physical coprimality.

**Proof.** At index \(k\), row \(i\) samples its first
\(d_i K_{j,i}(k)\) factors, where \(K_{j,i}(k)=k\) in the
increasing orientation and \(K_{j,i}(k)=j-k\) in the decreasing
orientation. The multiplicity of \(p\) in this product is
\[
 V_{j,i,p}(k)=\sum_{e=1}^{R_{i,p}}
       \left(\left\lfloor\frac{j\ell_i+d_iK_{j,i}(k)}{p^e}\right\rfloor
                      -\left\lfloor\frac{j\ell_i}{p^e}\right\rfloor\right).
\]
Each summand counts the multiples of \(p^e\) in the sampled
interval. Higher powers contribute zero by definition of \(L_i\).
Subtracting the slope \(\sigma_i d_i\sum_{e=1}^{R_{i,p}}p^{-e}\)
times \(k\) leaves, up to a constant depending on \(j,i,p\),
a negative sum of \(R_{i,p}\) fractional parts. Its oscillation,
the maximum minus the minimum over the indices, is at most \(R_{i,p}\).
After subtracting \(\beta_p k\) and a separate intercept at each
order, the retained ordinates therefore have size at most \(R_p/2\).

The lower convex minorant at an integer \(k\) is the minimum of
the ordinate there and the interpolants over all chords bracketing
\(k\). A perturbation of ordinates by at most \(R_p/2\) changes
the minorant's ordinates by at most \(R_p/2\) and its unit slopes
by at most \(R_p\). The
[Newton polygon product rule, Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)
puts the slopes of a common factor in both full slope lists,
counting widths. The common factor thus supplies the claimed
matching in the selected lists. For distinct monomials the even
ordinary moment is strictly positive on the real axis by its
independent exponential integral. A common rational factor must
have even degree. Rational coprimality extends over every
characteristic-zero field by a polynomial Bezout identity.

If the actual residual oscillations of the retained ordinates
after subtracting \(\beta_p k\) are \(R'_n,R'_N\), the matching
tolerance can be taken to be \(R'_n+R'_N\). This includes
cancellation between different retained rows.

For a finite prime set \(\mathcal P\), let \(L_{i,\mathcal P}\)
be the part of \(L_i\) supported on \(\mathcal P\), with
multiplicities. For \(N<J\), the retained charge satisfies
\[
 \mathcal L_{\mathcal P}:=\sum_i\log L_{i,\mathcal P}
  =\sum_{p\in\mathcal P}R_p\log p
  \ \le\sum_{i\ {\rm retained}}\min\{2Jd_i,|\mathcal P|\}\log(JM_i).
                                                               \tag{1}
\]
The first bound for a row counts the factors in the two ramps.
The second charges each prime only its maximum valuation times
its logarithm, which is at most \(\log(JM_i)\). When
\(JM_i<\min\mathcal P\), the exact row charge is zero.

## Growing-gap bands with arbitrary long retained ramps

Use the hypotheses and constants \(A,B,C,D_0,H,Q,J\) of the
[retained-scale band theorem](MULTIVARIATE_ARBITRARY_RETAINED_SCALE_DECOUPLING.md).
There are \(r\) selected rows with smaller exponents
\(c_iM+\ell_i\), \(1\le c_i\le C\), \(|\ell_i|\le BD\),
and \(D\le d_i\le AD\). Assume \(D\ge D_0\). At a terminal
\(4\le N<J\), take \(\mathcal P\) to be the primes greater than
\(HD\) dividing the selected band's product of marks. A physical
common root at this pair forces
\[
                    \mathcal L_{\mathcal P}\ge
                         \frac{r D\log M}{256 C J Q}.           \tag{2}
\]

Indeed, write \(m\) for the even marked order, \(m'\) for its
neighbor, \(t\ge r/Q\) for the active-box threshold, and
\(r_0\le t/(4J)\) for the count omitted within the selected group.
At \(p\in\mathcal P\), let \(E=v_p(x)\) for the unique mark
divisible by \(p\). By the cited construction, order \(m\) has
at least \(t\) active copies of the mark outside an interval
\(I\) of width at most one, and no active copy at a chosen
\(k_*\in I\). Order \(m'\) has no active copy at all. In either
order at most \(r_0\) omitted copies occur at each index, and
these copies and the active copies account for all selected
prime valuations.

Remove the common slope \(\beta_p\) from the full ordinates.
The binomial coefficients are units because \(p>HD\). The
order-\(m'\) slopes after removal have absolute value at most
\(r_0E+R_p\), a bound on the oscillation of all its ordinates.
After subtracting any of these slopes times the index, each
order-\(m\) ordinate outside \(I\) exceeds the one at \(k_*\)
by at least
\[
             (t-(m+1)r_0)E-(m+1)R_p
                     \ \ge\ \frac{3t}{4}E-JR_p .
\]
If this is positive, shared slope width is at most one. Thus
a physical contact forces \(E\le4JR_p/(3t)\). Each
\(p\in\mathcal P\) divides only one mark, so the large-prime
part \(P_{>HD}\) of their product satisfies
\(\log P_{>HD}\le4J\mathcal L_{\mathcal P}/(3t)\).
Combine this with the cited lower bound
\(\log P_{>HD}\ge D\log M/(160C)\) to obtain (2).

In particular, suppose \(D\to\infty\), \(M\le D^K\), and the
retained heights satisfy \(M_i\le M^L\), where \(K,L\) are fixed.
Then \(o(r)\) retained rows suffice for all-degree coprimality
without any additional retained gap bound. In fact there are at
most \(D\) marks of size \(O_{A,B,C,J}(M)\) and all the charged
primes exceed \(HD\), so \(|\mathcal P|=O(D)\). Equation (1)
gives \(o(rD\log M)\), contradicting (2). The first two pairs
and all terminals \(N\ge J\) follow from the cited band theorem's
low order and absolute angular inputs.

## Certified affine templates with arbitrary high retained rows

Fix a bounded-gap affine template and a compact admissible
fraction polytope passing the
[count-simplex test](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md)
for every \(4\le N<J\). Let \(r\) be its row count and \(M\) the
affine parameter. For \(M\) above a computable constant, a
physical common root forces
\[
       \mathcal L_{\mathcal P}\ge c\,r\log M,
                  \qquad \mathcal P=\{p:p\mid x(M)\},           \tag{3}
\]
where \(x\) is a certifying mark of the selected fraction vector
and \(c>0\) depends only on the template data and \(J\).

Let \(L_0\) be divisible by all the positive affine multipliers
and all nonzero cross-offsets in that test. Let \(D\) bound the
template gaps, and let \(\delta_{n,N}>0\) be the uniform
two-match distance. For a mark depth \(E=v_p(x(M))>v_p(L_0)\),
the selected ordinates before binomial coefficients differ
from \(rE\) times the nominal ones by at most
\(rJDv_p(L_0)\). Including the binomial coefficients adds at
most \(v_p(J!)\). Remove \(\beta_p k\) from the retained
ordinates and center the remaining interval at each order,
which adds error bounded by \(R_p\). The lower minorant bound
above implies that after division by \(rE\) and common slope
translation, every cross-order slope difference changes by at
most
\[
           \frac4E\left(JD v_p(L_0)+\frac{R_p+v_p(J!)}r\right).
\]
If this is below \(\delta_{n,N}/2\), two distinct slope matches
are impossible. Thus at a common root each prime of the mark
satisfies
\[
         v_p(x(M))\le C_1
           \left(JD v_p(L_0)+\frac{R_p+v_p(J!)}r\right),\qquad
         C_1=\max_{4\le N<J}\max(1,8/\delta_{N-1,N}).
\]
This inequality is automatic when \(E\le v_p(L_0)\).
If \(J=4\) the low and high order theorems suffice. Otherwise
summing against \(\log p\) gives
\[
  \log x(M)\le C_1\left(JD\log L_0+
                      \frac{\mathcal L_{\mathcal P}+\log(J!)}r\right).
\]
For all sufficiently large \(M\), every mark is at least \(M\)
and \(\log M\ge2C_1(JD\log L_0+\log(J!))\). This proves
(3) with \(c=1/(2C_1)\).

A fixed affine mark is \(O(M)\). The product of \(s\) distinct
primes is at least \((s+1)!\), so \(|\mathcal P|\) is
\(O(\log M/\log\log M)\). For a fixed \(L\), if retained heights
are at most \(M^L\), then
\[
                r_{\rm retained}
                    =o\!\left(r\,\frac{\log\log M}{\log M}\right)
\]
suffices for all-degree coprimality by (1) and (3). These rows
may have gaps greater than \(r\).

For example, append to any fixed ordinary core
\(\lfloor M/16\rfloor\) copies each of \((M-1,M)\) and \((M,M-1)\),
and two opposite rows of gap \(M^2\) with their smaller exponents
chosen independently in \([M^4,2M^4]\). The selected types pass
the [common-height test](MULTIVARIATE_ARBITRARY_RETAINED_SCALE_DECOUPLING.md#bounded-gaps-from-any-certified-affine-template).
The retained rows carry curvature bounded away from zero.
Each would be charged order \(M^2\log M\) by the full product
estimate. Equation (3) proves all-degree coprimality uniformly
over these independent high exponents for all sufficiently
large \(M\).
