# All-degree coprimality for rationally separated upper bands

Fix integers \(A\ge1\), \(B\ge0\), \(L\ge1\). Let \(J\ge4\)
be any integer at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Suppose for every consecutive pair \(n,N=n+1\) with
\(4\le N<J\), at least one of \(n,N\) does not divide \(L\).
This condition holds whenever \(L\) is a prime power, and for \(L=1\).

In a fixed ordinary profile, omit equal rows and select a nonempty
group of unequal rows. Write \(M_i=\max(a_i,b_i)\) and
\(d_i=|b_i-a_i|\). In the selected group choose multipliers
\(c_i\mid L\) and an integer \(M>0\) with
\[
             c_iM_i=M-e_i,\quad 0\le e_i\le BD,\qquad
             D\le c_i d_i\le AD .
\]
Here \(D\) is the minimum scaled gap in the group. These multipliers
are used for the arithmetic proof. The factorial arguments and
exponents in the profile remain the original ones.

Set \(K=J(A+B+L)\), \(H=JK\), and
\[
       \mathcal B=\sum_{\rm retained}d_i\log(JM_i),\qquad
       \delta=\frac1{64J^3L}.
\]
For this sum one can omit all unselected rows with \(JM_i\le HD\).
One can also omit any unselected row for which a multiplier \(c_i\mid L\)
places both exponents \(c_i a_i,c_i b_i\) in \([M-D/6,M]\).
The number and total gaps of these omitted rows are unrestricted.

**Theorem.** For \(D\ge D_0\), with \(D_0\) computable from \(A,B,L,J\),
the condition
\[
                            \mathcal B\le \delta D\log M
\]
implies ordinary SFC(2) in every degree over every characteristic-zero
field. There is no bound on the number of selected rows.
In particular the selected upper exponents may lie in several
rationally related bands separated by a distance of order \(M\).
The value of the absolute cutoff \(J\) is unevaluated.

## An aligned progression of marks

The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
allows total curvature below \(12\) throughout. For any selected
row,
\[
       \frac{d_i^2}{M_i}
           =\frac{(c_i d_i)^2}{c_i(c_iM_i)}
           \ge \frac{D^2}{LM},
\]
so \(M>D^2/(12L)\).

Fix \(4\le N<J\). Mark an order \(m\in\{n,N\}\) with \(m\nmid L\),
and write \(m'\) for the other order. For \(0\le k\le m\)
form the integer envelopes of the
[discrete mark construction](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-discrete-interval-of-eligible-marks):
\[
 f_+(k)=\min_{\rm selected,\ increasing}\bigl(me_i+(m-k)c_i d_i\bigr),
 \quad
 f_-(k)=\min_{\rm selected,\ decreasing}\bigl(me_i+k c_i d_i\bigr),
 \quad
 F=\max_k\min(f_+(k),f_-(k)).
\]
A missing orientation has minimum \(+\infty\).
The middle index gives \(F\ge\lfloor m/2\rfloor D\), and comparison
with any selected row gives \(F\le\min_i m(e_i+c_i d_i)\).
Thus \(F\le m(A+B)D\).

Take integers \(u\) in
\[
           [\max(F-D,mD/6),F)
\]
with the simultaneous congruences
\[
                     u\equiv mM\pmod L,\qquad
                     u\equiv g\pmod m,\qquad g=\gcd(m,L).
\]
The interval has length at least \(D/2\).
For \(m\ge4\) the lower bound \(F-D\ge mD/6\) holds, and for
\(m=3\) one uses \(F\ge D\).
The congruences are compatible modulo \(g\), and \(0<g<m\).
They define a progression of step \(\operatorname{lcm}(m,L)\le JL\).
In particular \(m\nmid u\). Put
\[
                            x_u=mM-u,\qquad y_u=x_u/L .
\]
Every \(y_u\) is integral, the \(y_u\)'s form a progression
of step size \(m/g\), and each is coprime to that step, since
\((L/g)y_u\equiv-1\pmod{m/g}\).

After stripping each row's smaller-end factorial,
a selected row has full order-\(m\) ramp \(mM_i-h\),
\(0\le h<md_i\). Its potential marked factor is \(x_u/c_i\),
which is integral. The associated scaled offset
\(u-me_i\) is a multiple of \(c_i\). It has not passed the full
ramp because \(u<F\le m(e_i+c_i d_i)\). The factor occurs
at coefficient index \(k\) precisely if that row's threshold
in \(f_+(k)\) or \(f_-(k)\) is at most \(u\).
Hence the indices with no mark form the nonempty interval
\(\{k:\min(f_+(k),f_-(k))>u\}\). Its width is at most one.
Indeed three consecutive such indices would give, at the middle
one, values greater than \(u+D\ge F\) in both finite branches.
Use the rightmost index for \(f_+\) and the leftmost one for \(f_-\).

## Large primes of the original factors

At order \(m\), another selected factor \(mM_i-h\) satisfies
\[
                    c_i(mM_i-h)-x_u=u-me_i-c_i h .
\]
At the other order,
\[
        mc_i(m'M_i-h)-m'x_u=m'u-mm'e_i-mc_i h .
\]
The latter offset is nonzero modulo \(m\), because
\(\gcd(m,m')=1\) and \(m\nmid u\).
All nonzero magnitudes here are at most \(HD\).

A prime \(p>HD\) dividing \(y_u\) divides \(x_u\) with the same
valuation. It is coprime to \(m,L\). All selected factors at order
\(m'\) and all nonmarked selected factors at order \(m\) are
\(p\)-units, by the two offset identities. Every marked factor
\(x_u/c_i\) has valuation \(E=v_p(y_u)\).
The retained rows omitted because \(JM_i\le HD\) also have unit
factors. For a retained row in the scaled top strip
\([M-D/6,M]\), all scaled factors at order \(m\) are strictly
greater than \(m(M-D/6)\ge x_u\). Thus this row contributes
no mark. Its offsets at the two orders are still bounded in
magnitude by \(HD\) and nonzero by the same argument.
All its factorial factors at both orders are \(p\)-units.

Define \(B_j\) as \(j!\) times the product of full order-\(j\)
factorial ramps of the rows counted in \(\mathcal B\).
The [deep-block comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
with valuation bounds \(v_p(B_m),v_p(B_{m'})\) now shows
\[
                E\le v_p(B_m B_{m'}^m)
\]
at a possible common root. In fact a larger \(E\) would
confine all shared slopes to the unmarked interval, since
every index outside it has at least one factor of valuation
\(E\), and any index is within distance \(m\) of the interval.
Width at most one contradicts even-degree parity.
The span of the \(x_u\)'s is at most \(D\), so a prime above
\(HD\) divides at most one mark. Therefore the product of
the parts of all \(y_u\)'s supported on such primes divides
the corresponding part of \(B_m B_{m'}^m\).

## Completion by the retained prime budget

Write \(k\) for the number of marks. For \(D\ge4JL\),
\[
                         D/(4JL)\le k\le KD .
\]
For each prime \(p\le HD\) dividing a \(y_u\), remove one term
of maximal \(p\)-valuation. If \(p\mid m/g\), then all
\(y_u\)'s are units. Otherwise the indices divisible by \(p^\nu\)
form a single residue class modulo \(p^\nu\). The remaining
product has \(p\)-valuation at most
\(\sum_{\nu\ge1}\lfloor(k-1)/p^\nu\rfloor=v_p((k-1)!)\).
The entire remaining product therefore divides
\((k-1)!\,(B_m B_{m'}^m)_{>HD}\), where the subscript retains only
prime factors above \(HD\).

Use \(D_0\) large enough that \(\log D\ge288JLH\).
Then the [elementary prime-count bound](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root)
\(\pi(HD)\le9HD/\log(HD)\) gives at most \(k/8\) removed
marks. The order factorials are \(p\)-units for \(p>HD\),
so only the retained ramps enter the large-prime part.
Each counted row contributes at most
\((m+mm')d_i\le J^2d_i\) factors of size at most \(JM_i\).
Thus the logarithm of the remaining product is at most
\[
        k\log(KD)+J^2\mathcal B
             \le \left(\frac23+\frac1{16}\right) k\log M .
\]
Here \(D\ge144K^3 L^2\) gives \(KD\le M^{2/3}\), and
\(J^2\delta D\le k/16\). On the other hand \(D\ge24L(A+B)\)
gives \(x_u> M\), and \(D\ge 12L(JL)^8\) gives
\(\log(M/L)\ge(15/16)\log M\).
The remaining \(y_u\)'s therefore give a logarithm greater than
\((7/8)(15/16)k\log M\). This contradicts the upper bound.
For example one can take
\[
 D_0=\max\{4JL,24L(A+B),144K^3L^2,12L(JL)^8,2^{1024JLH}\}.
\]

Terminal two follows from even positivity, terminal three from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients),
and \(N\ge J\) from the absolute cutoff. Rational Bézout identities
transfer the result to every characteristic-zero field.

For a balanced example let \(D\) be even and \(M=2D^2\).
Take one increasing row \((M-D,M)\) and two decreasing rows
\((M/2,M/2-D/2)\). The multipliers \(1,2,2\) give exactly
the same scaled upper endpoint and scaled gap in all three rows.
Their physical signed gap is zero and their total curvature
is one. One may retain any fixed balanced core, for example
\((1,2),(2,1)\), preserving zero signed gap and keeping
moderate curvature in this example. Its contribution
\(\mathcal B\) is bounded as \(D\) grows. The theorem separates
these original profiles in every degree for all sufficiently large
\(D\), even though the two upper bands are separated by \(M/2\).
