# Unpaired factorial envelopes and prime witnesses at several height scales

Use the [effective absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
and take \(J\ge16\) at least that cutoff. The
[terminal-fourteen theorem](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md)
leaves only \(15\le N<J\), with \(n=N-1\). The row counts below are
independent nonnegative integers and the core is arbitrary and fixed.

**Theorem.** Fix finitely many row types with smaller exponents
\(\ell_i(T)=c_iT^{e_i}+a_i\), where \(c_i,e_i>0\) and \(a_i\ge0\) are
integers, and fixed gaps \(d_i>0\) and orientations. Assume that
different exponents \(e_i\) have distinct \(2\)-adic valuations,
or that the stripped ramp polynomials from different degree blocks
are pairwise coprime over \(\mathbb Q[T]\) at the residual orders.
At each pair \(n,N\) the following tests give coprimality for all
\(T\) above a computable threshold independent of the row counts.

- For some present nonlinear block, at one order \(j\in\{n,N\}\)
  choose a prime dividing \(j\), and retain the present types of
  maximal slope valuation at this prime in that block. If \(c_*\)
  is the gcd of their slopes, require \(d_i/c_i\ge2/c_*\) for
  these types.
- For some present block, make the same test at each of the two
  orders, choosing a prime divisor of each order separately.
  This permits an affine block and arbitrary integral slopes.

Each gap hypothesis may be relaxed to the envelope test below.
Other blocks may have arbitrary gaps and counts. For common
slopes in a selected block, the condition is simply \(d_i\ge2\).
In particular no reflection or signed balance is required.

There is also a template-free version. Fix \(C\ge1\). Give each
appended unequal row an arbitrary integer scale \(1\le c_i\le C\).
Suppose both its exponents divided by \(c_i\) lie in \([H,H+W]\),
with \(H,W\) positive integers. Write
\[
                    \Theta(H)=
           \frac{\log_2 H\,\log_3 H}{\log_4 H}
\]
with iterated logarithms truncated below by one. There is a
computable \(\kappa_{J,C}>0\) such that the original moments are
coprime at every adjacent pair for all sufficiently large \(H\)
if \(W\le \kappa_{J,C}\Theta(H)\) and the envelope test passes
at both orders of each residual pair. In particular it passes if
all appended gaps are at least \(2C\), or all have a common scale
and gaps at least two. The thresholds are uniform
in the orientations, gaps, number of row types and multiplicities.
At \(C=1\) the same conclusion holds for all-unit gaps with arbitrary
counts and orientations.

## The normalized envelope test

At order \(j\), after stripping the scalar \((j\ell_i)!\), the
ramp of type \(i\) contributes the factors \(j\ell_i+h\) with
\(1\le h\le d_i k_i\), where \(k_i=k\) for an increasing
row and \(k_i=j-k\) for a decreasing row. In one degree
block of a fixed template, write \(A_i=a_i\). In a
template-free scaled band, write \(A_i=\ell_i\).
Choose a prime \(\pi\mid j\). Keep only the present types
with largest \(v_\pi(c_i)\) in the block, or in the entire
band in the template-free case. On this tier put
\[
    f(k)=\max_i (jA_i+d_i k_i)/c_i,\qquad F=\min_{0\le k\le j} f(k).
\]
Let \(c_*=\gcd_i c_i\) over this tier. Set \(x=u/c_*\), where
\(u\) is the least integer strictly greater than \(c_*F\) with
\(\pi\nmid u\). Then \(0<x-F\le2/c_*\). Every \(c_i x\) in the
tier is an integer and every \(c_i/c_*\) is a \(\pi\)-unit,
since all slopes have the same \(\pi\)-valuation. Consequently
\(h_i=c_i x-jA_i\) is positive and \(\pi\nmid h_i\) for every
tier row. Positivity follows since \(x\) exceeds the normalized
lower endpoints, which are bounded by \(F\).

A sufficient envelope test is that at each interior \(k\), at least
one attaining row has \(d_i/c_i\ge2/c_*\). More generally it suffices
that its next argument in the increasing direction reach \(x\).
This ensures that at least one grid index reaches \(x\). Every
such index samples the common normalized mark, by integrality of
\(c_i x\) on this tier. Thus \(I=\{k:f(k)<x\}\) is the nonempty
integer interval of unmarked indices. If it contained three
consecutive indices, the attaining row at their middle index
would reach \(x\) at one neighbor by the hypothesis. This is
impossible, so the width is at most one.

At the marked order \(j\), a proportional factor from any type
in this block has normalized offset \(x\). Proportionality to
a mark from type \(i\) means \(c_i(jA_l+h')=c_l(jA_i+h)\).
The right side has valuation \(v_\pi(c_l)\), so \(l\) is
in the maximal tier. A companion-order proportional factor
would instead satisfy
\[
                    j c_i(qA_l+h')=q c_l(jA_i+h).
\]
The left valuation is at least \(v_\pi(j)+v_\pi(c_i)\),
while the right valuation is at most \(v_\pi(c_i)\). This
uses \(\gcd(j,q)=1\) and the maximal tier. The companion
has no proportional factor at this mark.

A prime isolating this mark therefore yields a flat companion
Newton polygon and a horizontal width at most one in the marked
polygon. The valuations at \(I\) vanish and those outside
are positive. The Newton polygon product rule makes every
shared slope zero with cumulative width at most one. A common
factor over \(\mathbb Q\) must have degree at least two by the
strict real positivity of the even factorial moment. The argument
of the original moments is unchanged when their stripped scalar
factors are removed. The use of Newton polygons is as in the
[central-mark theorem](MULTIVARIATE_CENTRAL_AFFINE_MULTISCALE_MARKS.md#exact-central-coefficient-minima)
and [Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016).

## Good primes in a fixed template

Choose a finite prime set \(S\) containing the prime divisors
of order factorials, stripped core factors, leading ramp
coefficients and nonzero resultants of pairs of literal ramp
polynomials. These are from a fixed finite list. The different
degree blocks are coprime under the asserted parity condition:
a common zero of \(T^e+\alpha\) and \(T^f+\beta\) with positive
\(\alpha,\beta\) would give opposite signs on raising the
equations to powers \(f/\gcd(e,f)\) and \(e/\gcd(e,f)\).
Within a block, the only possible common zero of distinct
literals is between proportional factors. Hence a prime outside
\(S\) dividing the selected mark divides only proportional
literal factors.

For a nonlinear mark the polynomial \(j c_i T^{e_i}+j a_i+h\)
has at least two distinct roots. The effective bound of
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/pdf/1708.08290)
gives an outside prime for every sufficiently large integer \(T\):
an \(S\)-smooth value would satisfy
\(|G(T)|\le K|G(T)|^{1-\delta}\) with fixed \(\delta>0\).
For a two-order witness in one block, the two selected marks are
nonproportional by the companion-order test. As linear forms in
\(T^e\), their values obey a fixed equation \(A u+B v=1\) with
\(A,B\in\mathbb Q^\times\). If both were \(S\)-smooth this would
be an \(S\)-unit equation of bounded effective height by
[Bérczes, Evertse and Győry, Theorem 2.1](https://math.leidenuniv.nl/reports/files/2008-16.pdf).
Both mark values tend to infinity, so one has an outside prime
for every sufficiently large \(T\). There are finitely many
possible selected marks and pairs, so the thresholds can be
taken independently of the counts.

## Independent rows in a scaled band

Use the envelope test at orders \(j=n,N\) with
\(A_i=\ell_i\) and potentially different choices of \(\pi\).
Write \(x_j=j\ell_i+h\) for the actual integer mark selected
at each order and \(c_j^*\) for its selected type scale. Let
\(L_C=\operatorname{lcm}(1,\ldots,C)\) and put
\[
             B=\max(3,J,JM_C,C,J^2 C^2 W),
\]
where \(M_C\) is the largest core exponent, zero for the empty
core. A prime \(p>B\) dividing \(x_j\) isolates the mark.
In fact, for a factor \(G=q\ell_l+h'\) at either order \(q\),
\[
         D=j c_j^* G-q c_l x_j
\]
has absolute value at most \(jq C^2 W\): the normalized
quantities \(G/(q c_l)\) and \(x_j/(j c_j^*)\) are in
\([H,H+W]\). If \(p\mid G\) then \(p\mid D\), forcing this
determinant to vanish. At the marked order, equality means a
literal with the same normalized mark. At the other order
equality is excluded by the maximal-tier argument. Core factors
and binomial coefficients are units by the choice of \(B\).

A collision would make both \(x_n,x_N\) \(B\)-smooth. Then
\[
                 u_j=\frac{nN}{j}\frac{L_C}{c_j^*}x_j
                        \qquad (j=n,N)
\]
would be distinct \(B\)-smooth integers in
\([nN L_C H,nN L_C(H+W)]\). Distinctness follows from the
companion-order test. Their difference is at most \(J^2 L_C W\).
To exclude this, apply
[Stewart, Theorem 1 and inequality (6)](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/on_sequences-2.pdf)
to the increasing sequence whose largest prime factor is at most
\[
     y(x)=\max(3,J,JM_C,C,J^2 C^2\kappa_{J,C}\Theta(x)).
\]
Use a monotone adjustment on a fixed initial interval. For large
\(H\), \(y(u_j)\ge B\). Stewart gives a lower bound for any
consecutive pair \(u<v\) of that sequence of the form
\(u/(\log u)^{\delta(c y(v))}\), with computable \(c>0\) and
\(\delta(t)=\exp(t\log_2 t/\log t)\). Choose
\(c J^2 C^2\kappa_{J,C}<1/4\). The ratio
\(\log\delta(cy(v))/\log_2 v\) tends to
\(c J^2 C^2\kappa_{J,C}\), so beyond a computable bound it is
less than \(1/2\). The computability follows by elementary bounds
on each successive logarithm, as in the
[growing-height-band argument](MULTIVARIATE_CENTRAL_AFFINE_MULTISCALE_MARKS.md#varying-row-types-in-a-growing-height-band).
Throughout the interval between the proposed values, the
resulting lower bound on a consecutive smooth gap exceeds
\(H^{1/2}\), whereas \(J^2 L_C W<H^{1/2}\) eventually.
This is impossible. A rational Bézout identity then gives
coprimality over every characteristic-zero field.

### Unscaled unit-gap bands

Suppose \(C=1\) and every appended gap is one. At each smaller
height \(L\), let \(U,V\) be the two orientation counts. If
some \(U\ne V\), use that height and choose \(x_j=jL+1\) at
both orders. At the marked order the only proportional factors
come from that height with offset one, so the ordinates in
units of the mark valuation are \(V,U+V,\ldots,U+V,U\).
The lower polygon is the endpoint chord, which has nonzero
slope. The companion is flat because \(j\nmid x_j\).
Thus either outside prime proves coprimality, and the two
marks cannot both be smooth by the same interval argument.
If \(U=V\) at each height, take \(x_j\) one above the largest
central factorial argument at both orders. At a maximizing
height both orientations occur and the offset of this mark
is \(\lceil j/2\rceil+1\), which is nonzero modulo \(j\).
All central indices avoid the mark and each exterior index
samples it. This again gives horizontal width at most one.
The same interval argument finishes the proof.

The earlier [widening-band theorem](MULTIVARIATE_NARROW_HEIGHT_BAND_NATIVE_MARKS.md)
uses an unscaled envelope and a divisor threshold with a fixed
gap cap. The two-order construction uses smooth gaps to control
an increasing gap cap in the stated regimes. The
[reflected central-mark theorem](MULTIVARIATE_CENTRAL_AFFINE_MULTISCALE_MARKS.md)
also treats blocks containing mixed gaps when a reflected
central maximizer is present. The envelope condition here can
hold for mixed gaps, with a different finite test.
