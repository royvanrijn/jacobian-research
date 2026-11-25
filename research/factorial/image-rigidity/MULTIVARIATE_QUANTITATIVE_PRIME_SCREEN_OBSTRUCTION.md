# A quantitative obstruction to shallow simultaneous prime screens

Fix a terminal \(N\ge3\) and any nonempty open interval
\(I\subset(0,\infty)\). There are constants \(C_{N,I}>0\) and
\(Y_N\ge2\) such that for every integer \(Y\ge Y_N\) and every
integer \(L\ge1\), one can choose a fully shared ordinary profile
with distinct columns, zero signed total gap, and curvature in \(I\)
for which
\[
 \frac{M_j(z)}{M_j(0)}\equiv(1+z)^j\pmod {p^L}
 \qquad(j=N-1,N,\ \text{all primes }p\le Y).
 \tag{1}
\]
The row count \(r\) and maximum exponent \(H\) satisfy
\[
         \max(r,H)\le C_{N,I}\,2^{2Y\max(1,L-1)}.
 \tag{2}
\]
All rows can have the same positive absolute gap, depending only on
\(N\), with no reversed unequal row pairs.

This strengthens the [simultaneous fixed-prime obstruction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
by bounding the size of profiles as more primes are included. At
each of these primes, the normalized pair has the horizontal
Newton segment for both orders, with a shared degree-two
residual divisor. At any fixed precision \(L\), the obstruction
therefore persists for primes through a positive constant times
\(\log\max(r,H)\) on this family.

## Fixed valuation cancellation and enlargement to all small primes

Use the common gap \(D\) chosen from the fractions \(k/j\) in the
fixed-prime proof. It depends only on \(N\). Set \(Y_N=ND\).
Apply its valuation-cancellation construction only to the fixed
primes \(p\le ND\). It yields a fixed nonempty multiset of mirrored
rows, with smaller bases \(A_a\), and all factorial coefficient
ratios to the constant term have valuation zero at those primes.
Keep this multiset, of size \(r_0\), before the unit-neutralizing
repetitions in that proof.

Choose a fixed period \(B_0\) divisible by a sufficiently high
power of each prime \(p\le ND\). More precisely, the power has
exponent greater than every \(v_p(jA_a+h)\) for every original
row \(a\), \(j=N-1,N\), and \(1\le h\le jD\). Replacing any
\(A_a\) by a nonnegative integer congruent to it modulo \(B_0\)
preserves all those valuations. In each row, the factorial
coefficient ratio is a product or quotient of sampled factors
\(jA_a+h\) for \(1\le h\le jD\), so its valuation is unchanged.

For \(Y\ge ND\) put
\[
        B=B_0\prod_{ND<p\le Y}p .
\]
Independently in every row, choose a new smaller base congruent
to its old value modulo \(B_0\), and congruent to zero modulo every
prime in the product. The Chinese remainder theorem gives such a
base in \([1,B]\). Put increasing bases in \([1,B]\) and decreasing
bases in \([B+1,2B]\) by adding \(B\) to the latter. This keeps
all factorial coefficient ratios at valuation zero for the old
primes. At any new prime \(p>ND\), every sampled factor
\(jA_a+h\) reduces to \(h\), with \(1\le h\le ND<p\), so each
row factor is a unit already. These statements remain true after
adding a common nonnegative multiple of \(B\) to all bases.

## Unit precision and curvature within the size bound

Put
\[
 Q=\operatorname{lcm}\left(\operatorname{lcm}(1,\ldots,Y),
                         \prod_{p\le Y}p^{L-1}\right).
\]
For every \(p\le Y\), \(Q\) is divisible by \(p^{L-1}(p-1)\):
the coprime factors \(p^{L-1}\) and \(p-1\) each divide \(Q\).
Also \(B\mid B_0 Q\), because each new prime divides
\(\operatorname{lcm}(1,\ldots,Y)\). Repeat the shifted multiset
\(R=K B_0 Q\) times, with a fixed integer \(K\) to be chosen.
Every unit factorial ratio raised to the \(R\)-th power is one
modulo \(p^L\), by the order of the unit group
\((\mathbb Z/p^L\mathbb Z)^\times\). This also works for a
rational unit, by reducing its numerator and denominator after
canceling their \(p\)-powers. Including the binomial coefficients
then gives (1).

To choose the common shift, fix \(s\in I\) and \(\varepsilon>0\)
such that \([s/(1+\varepsilon),s/(1-\varepsilon)]\subset I\)
and \(\varepsilon<1\). Write \(u=r_0 R D^2/s\) and set
\(t=\lfloor u/B\rfloor\). Add \(tB\) to each base. The new
row maxima \(tB+A_a+D\) differ from \(u\) by at most \(2B+D\).
Choose \(K\) depending only on \(N,I\) so large that
\(2B+D\le\varepsilon u\); this is possible since \(R\ge KB\).
The curvature is the sum of \(r_0R\) terms \(D^2/(tB+A_a+D)\),
so it lies in the required interval. The shift preserves the
disjoint base intervals of the two orientations. Repetition
preserves zero signed total gap and the common absolute gap.
It gives \(r=r_0R\) and \(H\le (1+\varepsilon)r_0RD^2/s\).

It remains to bound \(Q\). For integer \(Y\ge2\), put
\(\Lambda(Y)=\operatorname{lcm}(1,\ldots,Y)\). The exponent of
each prime in \(Q\) is at most \(\max(1,L-1)\) times its exponent
in \(\Lambda(Y)\). The elementary bound
\(\Lambda(Y)\le2^{2Y-2}\) follows by induction. In fact
\(\Lambda(Y)/\Lambda(\lceil Y/2\rceil)\) divides
\(\binom{Y}{\lfloor Y/2\rfloor}\), since any newly included
prime power in \((\lceil Y/2\rceil,Y]\) contributes a carry in
the factorial valuation of this binomial coefficient. The ratio
contains at most one factor of each prime, because two successive
powers cannot both be in that interval. The binomial coefficient
is at most \(2^Y\), and at most \(2^{Y-1}\) when \(Y\) is odd,
using its equal middle partner. The induction gives the bound
of \(2^{2Y-2}\) in both parities, starting at \(Y=1\).
This bounds \(R\) and proves (2).

The common residual divisor in (1) is \((1+z)^2\), since
\(N-1\ge2\). The two normalized endpoint units place the entire
length of each Newton polygon at slope zero. The conclusion
concerns the valuation and bounded-precision residual tests.
It leaves open methods using higher precision depending on the
profile, or primes beyond the indicated range.
