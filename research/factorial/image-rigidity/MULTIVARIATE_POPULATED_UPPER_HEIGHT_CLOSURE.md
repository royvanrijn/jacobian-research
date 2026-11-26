# Deterministic coprimality from a populated upper height interval

Let \(J\ge4\) be an integer at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md). In an original ordinary two-monomial profile, remove equal rows and let \(H\) be the largest remaining exponent. A row with smaller exponent \(a\) and gap \(d>0\) straddles the integer \(b\) if \(a\le b<a+d\). A height \(b\) is called favorable if, for either orientation, the number of rows of that orientation starting at smaller exponent \(b\) is strictly larger than the number of rows of the opposite orientation straddling \(b\). The favorable orientation may vary with the height. The counts include all row multiplicities.

**Theorem.** In a family with \(H\to\infty\), suppose there is an integer interval or a fixed-stride arithmetic progression of \(T\) candidate heights in
\[
              [(1-1/(2J))H,H-1],
              \qquad\frac{T\log H}{H}\longrightarrow\infty .
\]
In the progression case, write its heights as \(a_0+d r\) for \(0\le r<T\), with fixed positive step \(d\), and require \(\gcd(Na_0+1,d)=1\) for every \(2\le N<J\). Suppose all but \(o(T/\log H)\) of the candidate heights are favorable. Then every sufficiently large profile in the family satisfies physical SFC(2) in every degree over every characteristic-zero field. The other rows may have arbitrary lower heights, gaps, orientations, and multiplicities. No assumptions are imposed on their curvature.

## An upper-order prime excludes every other-order factor

For each \(2\le N<J\), select a favorable candidate height \(a\) for which \(p=Na+1\) is prime. Such a choice exists for all sufficiently large \(H\). In the progression case, the primes are in the reduced class \(Na_0+1\bmod Nd\), since this residue is coprime to \(N\) and to \(d\). The fixed-modulus prime number theorem gives \(\gg_{J,d}T/\log H\) primes uniformly over the indicated range. More precisely, subtract the estimates at its endpoints. The main term has size comparable to \(T/\log H\), and the endpoint errors are \(O_{J,d}(H/(\log H)^2)=o(T/\log H)\). Modulus two reduces to excluding the prime two from the ordinary estimate. Discarding the unfavorable candidates still leaves a prime. [Bennett, Martin, O'Bryant and Rechnitzer, Theorem 1.3 and (1.15)](https://personal.math.ubc.ca/~gerg/papers/downloads/EBPAP.pdf)

The upper height bound is what permits arbitrary other rows. Indeed,
\[
            p>(N-1)H,\qquad 2p>NH,
\]
because \(a\ge(1-1/(2J))H\) and \(N<J\). Thus every factor in the order-\((N-1)\) full stripped ramps is a \(p\)-unit. In the order-\(N\) ramps the only possible nonunit is the literal factor \(p\), with valuation one. A row's full ramp contains it exactly when that row straddles \(a\).

Reverse the columns if necessary so that the favorable orientation is increasing. Every increasing row starting at \(a\) contributes the literal \(Na+1\) at each index \(k\ge1\) and contributes none at index zero. The other increasing rows cannot decrease the valuation relative to zero. The decreasing rows can lose at most one copy each, and only if they straddle \(a\). Therefore every positive-index ordinate of the upper coefficient array exceeds its ordinate at zero. The binomial coefficients are units since \(p>J\). The smallest slope of the upper lower Newton polygon is the minimum of its secants from index zero, so every slope is strictly positive. The lower-order polygon is horizontal. The Newton polygon product rule excludes a rational common factor. A rational Bezout identity gives the same conclusion over any characteristic-zero field. The absolute angular theorem handles \(N\ge J\). [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

The fixed interval is used only to guarantee the existence of the primes. For a single explicit profile the test is exact as soon as such a favorable prime is found at each of its remaining terminals.

## Sparse deterministic regularization of arbitrary ordinary shapes

This gives two explicit extensions of any family of ordinary profiles with largest exponent \(M\to\infty\). Each extension has the same signed degree difference as the original, has maximum exponent \(M(1+o(1))\), adds curvature tending to zero, and satisfies coprimality in every degree for sufficiently large \(M\). Original rows and coefficients are retained.

Write \(T=\lceil M\log\log M/\log M\rceil\), for sufficiently large \(M\).

**Bounded number of added rows.** For each \(2\le N<J\), take the least \(a_N\in[M+T,M+2T-1]\) such that \(Na_N+1\) is prime. Such an \(a_N\) exists by the same prime estimate. Add an increasing unit-gap row at \(a_N\) and a decreasing unit-gap row at \(M\). The new maximum \(H\le M+2T\), and \(a_N\) lies in the requisite upper height range for all sufficiently large \(M\). No decreasing row straddles \(a_N\), since every original row ends at or below \(M\) and the added decreasing rows end at \(M+1\). Thus the witness works for each pair. There are at most \(2(J-2)\) new rows and their total curvature is at most \(2(J-2)/(M+1)\).

**Predetermined rows with no search for witness primes.** Let \(A\) be the least even integer greater than \(M\). At each height \(a=A+2r\), \(0\le r<T\), add one increasing unit-gap row at smaller exponent \(a\) and one decreasing unit-gap row at smaller exponent \(a-1\). Every such \(a\) is favorable. The decreasing row immediately below it ends at \(a\), the next one starts at \(a+1\), and the original rows all end below \(A\). The progression has step two and \(Na+1\equiv1\bmod 2N\). Its heights fit the upper range for all sufficiently large \(M\), so the theorem applies. The added curvature is \(O(T/M)\).

Both extensions also preserve the limiting constant-normalized coefficient arrays in the physical common variable. For a fixed order \(j\), divide the original and extended moment polynomials by their respective constant coefficients. One added pair, increasing at smaller exponent \(a\) and decreasing at smaller exponent \(b\), multiplies coefficient \(k\) of the constant-normalized array by
\[
                    \prod_{r=1}^k\frac{ja+r}{jb+j-k+r}.         \tag{1}
\]
For the bounded extension \(b=M\), and the ratio in each factor lies between \(1\) and \(1+2T/M\) for large \(M\). Multiplying over the bounded number of pairs gives \(1+o(1)\) uniformly in \(k\) at each fixed \(j\). In the predetermined extension \(b=a-1\), each numerator in (1) exceeds its denominator by \(k\). Thus the logarithm of the full coefficient multiplier is between zero and \(O_j(T/M)\). The multipliers therefore tend to one for every fixed finite collection of orders without rescaling the common variable.

In particular, when a fixed core and diffuse tails approach any positive Gaussian strength, including an exceptional strength, these extensions preserve the limiting curvature, signed degree, and any convergent finite physical coefficient arrays. They give coprime exact integer profiles approaching the same limits. For unrestricted original tails without the favorable upper-height condition, the existence of such extensions leaves coprimality open.

## Exact checks

The [finite checker](verify_multivariate_populated_upper_height_closure.py) tests three moving mixed profiles with broad straddling rows and favorable heights in both directions, and both forms of regularization for three base profiles. It verifies 90 private-prime certificates and 1,260 coefficient valuations, as well as 2,250 bounded-extension and 450 predetermined-extension coefficient multipliers. These checks use the original factorial valuations and independent row products.
