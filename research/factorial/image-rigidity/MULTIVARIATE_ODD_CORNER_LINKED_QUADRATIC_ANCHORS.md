# Linked quadratic specializations at odd-terminal corners

Let the consecutive orders be \(m=2q\) and \(n=m+1\), where \(m\ge6\). Consider an original ordinary two-monomial profile in the simultaneous corner of the [envelope classification](MULTIVARIATE_UNIVERSAL_ENVELOPE_MARKS.md#an-exact-mark-window-in-every-profile). The rows may have any heights, gaps, and multiplicities permitted by that classification. After a possible reversal of both columns, its conditions can be written as follows. Let \(B=1+\max_i a_i\) for the smaller exponents \(a_i\), and put \(t_i=B-a_i\). Then \(t_i\ge1\), and
\[
 \begin{array}{ll}
 d_i\le t_i & \text{for increasing rows},\\
 2d_i\le m t_i & \text{for decreasing rows}.
 \end{array}                                                   \tag{1}
\]
There is an increasing row \((t_i,d_i)=(1,1)\) and a decreasing row \((t_i,d_i)=(1,q)\).

**Theorem.** For every such profile, at least one of two explicitly selected affine moduli \(aB+b\) must divide a nonzero integer determined just by the deficits \(t_i\), the gaps, orientations, and multiplicities if the physical moments have a common root. The two chosen moduli are coprime. More generally, a list of \(m-1\) rational specializations gives simultaneous necessary divisor conditions on the original profile, including exact prime-power conditions from two scalar remainders. In particular there are only finitely many exceptional translations with these row data fixed. Their possible values lie in an explicit divisor list.

This treats every simultaneous corner above terminal five without any analytic estimates. The resulting divisors can grow with the moving row data, so there is still no uniform exclusion over arbitrary profiles.

## The marked polynomial always drops to a quadratic

Translate all smaller exponents to \(X-t_i\), retaining the gaps and orientations. Reverse the common coefficient variable, so index \(k\) counts the decreasing rather than increasing prefix lengths. Put
\[
 C_{j,k}(X)=\binom jk
 \prod_{i:\,+}\prod_{r=1}^{(j-k)d_i}(j(X-t_i)+r)
 \prod_{i:\,-}\prod_{r=1}^{kd_i}(j(X-t_i)+r),\qquad
 F_j(X,z)=\sum_{k=0}^j C_{j,k}(X)z^k.                            \tag{2}
\]
These have integer polynomial coefficients in \(X\). At \(X=B\) they are precisely the ordinary reversed parents after stripping the smaller-end factorials. At every integer \(X\ge\max_i t_i\) the profile is likewise ordinary and both endpoints in \(z\) are nonzero.

Choose a marked order \(s\in\{m,n\}\), and \(1\le h\le3q-s\), so there are \(q\) choices in the even order and \(q-1\) in the odd order. Write \(j'\) for the other order, \(X_0=-h/s\), \(g=\gcd(s,h)\), and \(v(X)=(sX+h)/g\). At \(X_0\) a literal factor of the marked order is zero only at \(r=s t_i+h\) in row \(i\). Increasing rows cannot reach that factor by (1). At indices \(k=0,1,2\) the decreasing rows cannot reach it either, as \(2d_i\le m t_i\le s t_i\). For \(k\ge3\) the anchored decreasing row of gap \(q\) does reach it because \(s+h\le3q\). Consequently
\[
             P=F_s(X_0,z)\text{ is quadratic, with all three coefficients nonzero}.
\]
The other specialization \(Q=F_{j'}(X_0,z)\) has no zero coefficient: any zero factor would require \(s\mid j'h\), impossible since the orders are coprime and \(0<h<s\). Thus \(Q\) has degree \(j'\).

## Divisors from the linked specializations

Form the square Sylvester matrix of \(F_s(X,z)\) and \(F_{j'}(X,z)\), using their original degrees \(s,j'\), and put \(D_{s,h}\) equal to the gcd of the reduced numerators of all its minors of size \(s+j'-1\) evaluated at \(X_0=-h/s\). Zero entries are allowed, with gcd zero if all are zero. A common root at any ordinary integer \(X\) forces
\[
                             v(X)\mid D_{s,h}.                 \tag{3}
\]
Indeed the even ordinary moment has no real roots, by its independent exponential integral. A rational common factor therefore has degree at least two. The Sylvester map on multipliers of degrees less than \(j'\) and \(s\) has kernel of at least that dimension, so all the indicated minors vanish at \(X\). For any one of these integer polynomials \(f\), evaluation gives \(f(X)\equiv f(X_0)\pmod{v(X)}\) in \(\mathbb Z[1/(s/g)]\). This follows from \(v(X)=(s/g)X+h/g\) and \(\gcd(v(X),s/g)=1\). In particular \(v(X)\) divides the numerator of \(f(X_0)\) whenever \(f(X)=0\). This proves (3), with full prime powers.

The same padded Sylvester map at \(X_0\) has kernel of dimension exactly \(\deg\gcd(P,Q)\). In fact, if \(H=\gcd(P,Q)\), every relation is given by multipliers \((Q/H)U,-(P/H)U\). The first degree bound is precisely \(\deg U<\deg H\), and the second allows all such \(U\) since \(\deg P=2\le s\). It follows that \(D_{s,h}\ne0\) exactly when \(P\) does not divide \(Q\) over \(\mathbb Q\). No general nonvanishing assumption on the full specialized resultant is required.

The \(m-1\) moduli from the fan have pairwise disjoint prime support above \(m(m+1)\). In fact for two distinct choices \((s,h)\), \((t,k)\), the common primes divide \(t h-s k\), which is nonzero. This follows by subtraction of the raw linear forms. For \(s=t\) it is immediate, and for \(s\ne t\) equality would force \(s\mid h\). Thus one can impose (3) simultaneously even on primes that are not isolated in the ordinary factorial ramps.

There is an exact scalar version at primes where the leading coefficients of \(P,Q\) are units. Fix a prime \(\ell\mid v(X)\) with that property, and let \(r_0+r_1z\) be the remainder of \(Q\) upon division by \(P\). Then
\[
           v_\ell(D_{s,h})=\min\{v_\ell(r_0),v_\ell(r_1)\},
           \qquad
           v_\ell(v(X))\le\min\{v_\ell(r_0),v_\ell(r_1)\}
                    \quad\text{at a common root}.              \tag{4}
\]
Use \(v_\ell(0)=+\infty\). To see the equality, work in the local ring at \(\ell\). In the padded map, the columns \(z^2Q,\ldots,z^{s-1}Q\) are triangular on the highest \(s-2\) output coefficients with diagonal the leading coefficient of \(Q\). The other columns lie in degrees at most \(j'+1\). Eliminating these highest outputs gives an identity block and the ordinary Sylvester map of \(P,Q\) of degrees \(2,j'\). Divide its output by \(P\), using the unit leading coefficient. The columns \(P,zP,\ldots,z^{j'-1}P\) become an identity block. The two remaining remainder columns are multiplication by \(Q\) on the module with basis \(1,z\) modulo \(P\). Together the two identity blocks have size \(s+j'-2\). The minors of size \(s+j'-1\) generate the ideal of entries of the remaining two-by-two block. Its first column consists of \(r_0,r_1\), and the second is obtained from these by multiplication by \(z\) modulo \(P\). The ideal is precisely \((r_0,r_1)\), proving (4). This gives a two-remainder check without forming full minors at the specified primes.

## One of two nonzero divisors at every corner

Set \(\Delta=\sum_{i:\,+}d_i-\sum_{i:\,-}d_i\). Take any prime \(p\mid s/g\), and put \(e=v_p(s/g)\). Every linear factor of \(Q\) has valuation \(-e\) at \(X_0=-h/s\): its numerator when written over \(s/g\) is \(-j'h/g\) modulo \(p\) and is a unit. Thus the coefficient ordinates are, up to their binomial valuations, a common affine function of slope \(e\Delta\) in the reversed index. The binomial valuations are nonnegative and the endpoint contributions are zero. The Newton polygon of \(Q\) has only slope \(e\Delta\).

Let
\[
 W_s(h)=\frac{C_{s,2}(-h/s)}{C_{s,0}(-h/s)}
      =\binom s2\,
        \frac{\displaystyle\prod_{i:\,-}
                 \ \prod_{r=1}^{2d_i}(s t_i+h-r)}
             {\displaystyle\prod_{i:\,+}
                 \ \prod_{r=0}^{2d_i-1}(s(t_i-d_i)+h+r)}.       \tag{5}
\]
All factors here are positive by (1). If \(P\) divides \(Q\), both slopes of \(P\), counted with width, must equal \(e\Delta\). In particular
\[
                             v_p(W_s(h))=2e\Delta.             \tag{6}
\]
We use this necessary condition at the two marked anchors selected below. The Newton polygon product rule over \(\mathbb Q_p\) is unchanged when denominators of the polynomials are cleared. [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

For the even order \(s=m\), the change from \(h\) to \(h+1\) in \(v_p W_m\) is
\[
 \sum_{i:\,-}
     \bigl(v_p(m t_i+h)-v_p(m t_i+h-2d_i)\bigr)
 +\sum_{i:\,+}
     \bigl(v_p(m(t_i-d_i)+h)-v_p(m(t_i-d_i)+h+2d_i)\bigr).
                                                               \tag{7}
\]
It is obtained by cancelling the overlapping interval factors in (5).

- If an odd prime \(p\) divides \(m\), choose \(h=p-2,p-1\). Both belong to \([1,q]\), because \(p\mid q\). Both are \(p\)-units, so they give the same \(e=v_p(m)\) and the same target in (6). At the increment \(p-2\) every unshifted term in (7) is a unit, so all contributions are nonpositive. The increasing unit row contributes a strictly negative term, because the shifted argument equals \(p\) at that row.
- Otherwise \(m\) is a power of two. Since \(m\ge6\), we have \(4\mid m\) and \(q\ge4\). Take \(p=2\) and \(h=1,3\), which again give the same target in (6). The increment \(1\) contributes zero in (7), as all arguments are odd. At the increment \(2\), every unshifted argument has valuation one, and every shifted argument is even. The contributions are nonpositive and the increasing unit row gives a strict drop, as its shifted argument is four.

In both cases \(v_p W_m\) strictly decreases between the selected anchors. Hence (6) cannot hold at both. At least one selected anchor has \(D_{m,h}\ne0\). The two raw moduli are consecutive in the odd-prime case. In the power-of-two case they differ by two and both are odd. Thus the selected primitive moduli are coprime. For a fixed deficit and gap profile, choose a nonzero \(D_{m,h}\). Every exceptional integer \(X\) is obtained by solving \((mX+h)/g=d\) for a positive divisor \(d\) of \(D_{m,h}\), with \(X\ge\max_i t_i\). This proves the finite translation claim.

The argument uses the exact increasing unit row forced by the corner classification. It does not bound its multiplicity or restrict any of the other rows. At terminal five the known [ordinary angular separation](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md) supplies coprimality directly. At even terminals, the universal envelope mark supplies the alternative arithmetic route. Neither route yet bounds the possible moving-prime overlap or the nonzero divisor in (3) uniformly over arbitrary tails.

## Every fixed relative profile has a finite divisor list

For an arbitrary finite ordinary template, fix its offsets, gaps, orientations, and multiplicities. Translate every entry by the same nonnegative integer \(T\). Given the [absolute eventual cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md), only finitely many terminal orders need certificates. The [first two pairs](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md) and terminal five are already separated.

At a remaining pair without a simultaneous corner, choose a permissible envelope mark \(x\) at \(T=0\). At the specialization \(T=-x/s\), the marked polynomial is supported exactly at the interval \(I(x)\) of width at most one, and the other polynomial has all coefficients nonzero. These claims follow from the same literal-factor comparisons, using the original rather than reversed index if desired. After factoring out a power of \(z\), the marked polynomial has degree at most one. The other polynomial has nonzero constant, so the padded Sylvester kernel has dimension at most one. The gcd of the corresponding rank-minor numerators is nonzero. A physical common root at an ordinary translation forces \((sT+x)/\gcd(s,x)\) to divide it, by the evaluation argument of (3).

At every simultaneous corner use the nonzero linked-anchor divisor, with \(X=T+B\) where \(B\) comes from the template. The union of the finite divisor lists over the tested terminals contains every translation at which the ordinary profile fails consecutive coprimality. Each listed translation can then be checked by exact ordinary polynomial gcds. This finite arithmetic reduction is for fixed relative data, including the multiplicities. It places no bound on the range of those data across all original supports.

## Exact checks

The [finite checker](verify_multivariate_odd_corner_linked_quadratic_anchors.py) verifies 1,330 quadratic specializations, 70 strict linked comparisons, and 140 unmarked valuation arrays. It checks 60 nonzero anchor remainders directly, eight padded rank computations, 40 translated minor congruences, and six local remainder ideals, including three depth-two or deeper cases. The exact comparisons above prove the statements for arbitrary corner profiles.
