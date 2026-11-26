# All-degree coprimality from the large-prime part of retained factorials

Consider an original ordinary two-monomial profile with distinct
exponent columns, after omitting equal rows. Select a nonempty
group of rows, with gaps \(d_i\in[D,AD]\), where \(D\) is the
least selected gap and \(A\ge1\) is an integer. Let \(M\) be the
largest exponent among the selected rows. Their larger exponents
are \(M-e_i\), where
\[
                              0\le e_i\le B D
\]
for an integer \(B\ge0\). Either orientation may occur, and the
retained rows need not have exponents at most \(M\).
Take \(J\ge4\) at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and put \(H=J^2(A+B)\). For a retained row let \(M_i=\max(a_i,b_i)\)
and \(d_i=|b_i-a_i|\). Define the retained log height above the
cutoff by
\[
                 \mathcal L=\sum_{\substack{i\ {\rm retained}\\
                                             J M_i>H D}}
                                     d_i\log(JM_i).
\]

**Theorem.** Suppose
\[
                   D\ge D_0(J,A,B),\qquad
                   \mathcal L\le\frac{D\log M}{64J^3},
 \qquad
 D_0=\max\{1152A^3,\ 64J^3,\ 2^{512JH}\}.
 \tag{1}
\]
Then the original binomial satisfies physical SFC(2) in every
degree over every characteristic-zero field. The criterion has
no row-count bound. Retained rows with \(J M_i\le HD\) have
no contribution to it. In particular, if one selected row has
gap \(D\) and every other unequal row has upper exponent at
most \(JD\), the original profile is regular in every degree
for all sufficiently large \(D\), regardless of their number
or their total signed or unsigned gap. One may use \(A=1,B=0\)
for this case.

If the selected largest exponent \(M\) is also the global
largest, a sufficient alternative to the second condition of
(1) is that the sum of retained gaps with \(J M_i>HD\) be at
most \(D/(128J^3)\). The argument uses the exact prime parts
of the original coefficients. The numerical value of \(J\)
is unevaluated. The proof also works with the separate conditions
\[
             D\ge\max(1152A^3,64J^3),\qquad \log D\ge144JH
\]
in place of \(D\ge D_0\). Thus the parameters \(A,B\) need
not be fixed as \(D\) grows.

## A discrete interval of eligible marks

Take \(m=N\) in an adjacent pair with \(4\le N<J\), and
write \(m'=m-1\). Strip the smaller-end factorial of each row
at each order. At the marked order, row \(i\) in the group has
possible factors \(m(M-e_i)-h\), \(0\le h< m d_i\).
We will mark the integer \(x_u=mM-u\) wherever it occurs.

For integer \(k\in[0,m]\), set
\[
 f_+(k)=\min_{i\ {\rm increasing}}(me_i+(m-k)d_i),\qquad
 f_-(k)=\min_{i\ {\rm decreasing}}(me_i+k d_i),
 \qquad f(k)=\min(f_+(k),f_-(k)).
\]
A missing orientation has minimum \(+\infty\).
Let \(F\) be the largest value of \(f\) on this integer grid.
Each finite branch \(f_+\) decreases by at least \(D\) per
step, and each finite branch \(f_-\) increases by at least \(D\)
per step. Moreover
\[
             \lfloor m/2\rfloor D \le F
                  \le \min_i m(e_i+d_i)\le m A D.
\]
The lower bound uses a middle index and \(e_i\ge0\). The
first upper bound holds at every index by comparing \(f(k)\)
with any one of the rows. The last uses a selected row with
\(e_i=0\).

For each
\[
                           F-D\le u<F,\qquad u\equiv1\pmod m,
 \tag{2}
\]
we have \(0\le u<\min_i m(e_i+d_i)\). A selected row samples
\(x_u\) at index \(k\) precisely when its threshold in the
definition of \(f_\pm(k)\) is at most \(u\). In fact this
threshold condition is exactly the condition that the local
offset \(u-me_i\) lies in the relevant rising block. The
bound on \(u\) ensures that it has not passed the end of any
row's full range. Rows whose local offsets are negative
contribute no mark.

Consequently the indices where no mark occurs form the
nonempty interval \(I=\{k:f(k)>u\}\cap\mathbb Z\), of
width at most one. To see the width bound, if three consecutive
indices belonged to \(I\), then at their middle index each
finite branch would be greater than \(u+D\ge F\).
Use the right-hand index for the decreasing branch and the
left-hand index for the increasing branch. This contradicts
the definition of \(F\). The interval property follows because
the superlevel set of the decreasing branch is a prefix
and that of the increasing branch is a suffix.

Let \(k_0\) be the number of integers in (2). They give an
arithmetic progression of marks of step \(-m\) and span less
than \(D\), with \(\gcd(x_u,m)=1\). They satisfy
\[
                        D/m-1\le k_0\le D/m+1.
 \tag{3}
\]
This construction also applies to a group with only one sign.

## A product bound at a possible common root

Put
\[
 B_j=j!\prod_{i\ {\rm retained}}\prod_{h=1}^{j d_i}
                          (j\ell_i+h),\qquad
                    \ell_i=\min(a_i,b_i)
                       \quad (j=m,m').
\]
For a positive integer \(b\), write \(b_{>HD}\) for the part
supported on primes larger than \(HD\), with multiplicities.

Suppose a prime \(p>HD\) divides a mark \(x_u\).
Then it divides neither \(m\) nor \(m'\). At order \(m\)
the difference of any other group factor from \(x_u\)
is \(u-me_i-h\), of absolute value at most \(HD\).
So only copies of the marked integer itself contribute
\(p\)-valuation. At the other order,
\[
 m(m'(M-e_i)-h)=m'(mM-u)+(m'u-mm'e_i-mh).
\]
The constant term on the right is nonzero because \(m\nmid u\),
and its absolute value is at most \(HD\), by the bounds
\(u<mAD\), \(e_i\le BD\), and \(0\le h<m'd_i\).
In particular all the group factors at order \(m'\) are units.

Let \(E=v_p(x_u)\), \(U=v_p(B_m)\), and \(L=v_p(B_{m'})\).
The cleaned marked ordinates, with all copies of \(x_u\)
removed, lie in \([0,U]\), and the other polygon's slopes lie
in \([-L,L]\). Outside the interval \(I\), at least one marked
factor occurs. If \(E>U+mL\), comparison with the nearest
index of \(I\), after subtracting a slope in \([-L,L]\),
excludes every outside index from a minimizing face.
Thus the marked slopes usable by a common factor have
total width at most one. The
[deep factorial-block comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
then excludes any common factor, because its degree would
have to be even by positivity of the even moment on the real
axis. Accordingly a possible common root forces
\[
                   E\le v_p(B_m B_{m'}^{m}).
\]
A prime \(p>HD\) divides at most one mark, since the span
of the progression is less than \(D\). Therefore the large-prime
part of the product of all marks divides
\((B_m B_{m'}^m)_{>HD}\).

For each prime \(p\le HD\) dividing at least one mark, remove
one mark of highest \(p\)-valuation in the original progression.
The marks are units at primes dividing \(m\). At every other
prime, for each \(p^r\), the indices of divisible terms lie
in one class modulo \(p^r\). After the removal, at most
\(\lfloor(k_0-1)/p^r\rfloor\) such indices remain.
Thus the small-prime part of the remaining product divides
\((k_0-1)!\), by Legendre's formula. The two arguments together
give a necessary condition for a common root when \(k_0>\pi(HD)\):
\[
       \bigl(\min_u x_u\bigr)^{\,k_0-\pi(HD)}
                       \ \le\ (k_0-1)!\,(B_m B_{m'}^m)_{>HD}.
 \tag{4}
\]
The small-prime count can be replaced by the number of those
primes actually dividing the progression. This yields a
certificate in terms of the exact retained prime product,
even outside the sufficient log-height bound in (1).

## Uniform estimates for the certificate

The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
handles \(S\ge12\). Otherwise \(M>D^2/12\) by the
contribution of a single selected row. The choice of \(D_0\)
gives \(M>2AD\), so all marks exceed \(M\).
Equations (2) and (3) give \(k_0\ge D/(2J)\) and \(k_0\le2AD\),
hence \(\log k_0\le(2/3)\log M\). For the latter bound it
suffices that \(D\ge(2A)^3\,12^2\).

The elementary estimate \(\pi(x)\le9x/\log x\) for \(x>1\)
follows by summing the dyadic binomial estimates
\(\sum_{2^{i-1}<p\le2^i}\log p\le\log\binom{2^i}{2^{i-1}}\)
and splitting the prime count at \(\sqrt{x}\).
In particular \(\log D\ge144JH\) implies
\(\pi(HD)\le D/(16J)\le k_0/8\).
This log bound holds at the stated \(D_0\).

A retained row with \(JM_i\le HD\) and either order factorial
has no factors contributing to \((B_m B_{m'}^m)_{>HD}\).
Each remaining row contributes at most \(m^2 d_i\) factors
to \(B_m B_{m'}^m\), each at most \(JM_i\).
Thus
\[
 \log (B_m B_{m'}^m)_{>HD}
     \le J^2\mathcal L
     \le \frac{D\log M}{64J}
     \le \frac{k_0\log M}{32}.
\]
Taking logarithms in (4) now gives a lower bound at least
\((7/8)k_0\log M\) and an upper bound at most
\((2/3+1/32)k_0\log M\), a contradiction.
If \(M\) is globally maximal, then \(M>J\) and
\(\log(JM_i)\le2\log M\), proving the alternative sufficient
condition.

## Retained rows near the upper endpoint

In addition to the retained rows below \(HD/J\), one may omit
from \(\mathcal L\) all retained rows whose two exponents lie
in \([M-D/5,M]\). Their numbers, orientations and total gaps
are unrestricted. For every \(m\ge4\) the eligible marks have
\[
                     u\ge F-D\ge(\lfloor m/2\rfloor-1)D
                                      \ge mD/5 .
\]
Every rising factor of such a row at order \(m\) is strictly
greater than \(m(M-D/5)\ge x_u\), since its smaller exponent
is at least \(M-D/5\). The difference from \(x_u\) is positive
and at most \(u\le mAD\). At order \(m'\), multiplying a
factor by \(m\) and subtracting \(m'x_u\) gives an integer
congruent to \(m'u\pmod m\), which is nonzero and of magnitude
at most \(mm'AD\). The magnitude bound uses that the factor
lies between \(m'(M-D/5)\) and \(m'M\). Both orders therefore
give units at every \(p>HD\) dividing any mark. These rows
can be removed from the retained product in the large-prime
part of (4). The small-prime proof already uses only the
progression of marks. More generally, at any marked pair
the same omission applies to a row whose endpoints lie in
\([M-(A+B)D,M]\) and whose full sampled interval of stripped
factors at order \(m\) misses the progression in (2).
The bounds on both offset magnitudes remain \(HD\).

The first two pairs follow from
[even positivity and quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
The absolute cutoff handles all \(N\ge J\). The argument applies
to the same variable in both original moments, and rational Bézout
identities give the result over every characteristic-zero field.
