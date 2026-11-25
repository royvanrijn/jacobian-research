# Two high roots at an arbitrary fixed shift

Let \(X\) bound two positive coprime integers \(U,V\) satisfying a
signed equation
\[
                        \epsilon_U U+\epsilon_V V=t,\qquad t\ne0.
\]
The target is fixed in each count and can vary with \(X\). Its
absolute value is at most \(2X\). There is no radical hypothesis
on \(t\).

## Dictionary bounds

Write each term as a fixed positive integer coefficient times a
product of positive integral roots in distinct indexed valuation
positions. The positions are pairwise coprime in each term, and
the two complete terms are coprime. The coefficients may share
primes with roots on their own side. Fix the position lists and
their exponents before the height limit at each fixed budget
\(\lambda\). Work in dyadic boxes for the roots. The actual raw
root height exponents are \(\alpha_i\), with
\(\sum_i\alpha_i\le C_0\lambda\) for a fixed \(C_0>0\). In position
\(i\) the dictionary has at most \(X^{W_i+o_X(1)}\) values, with
\(W_i\le\alpha_i\). Write \(W_U,W_V\) for the two sums of
these costs. Let \(D\) be the count exponent of distinct pairs.

1. If \(\max(U,V)\) is comparable to \(X\) and the fixed coefficients
   have numerical height exponent tending to zero in the budget
   limit, then
   \[
                             D\le (W_U+W_V)/3+o(\lambda).           \tag{FT1}
   \]
2. For any fixed \(a>0\), let \(B\) be the summed dictionary cost of
   the positions on both terms with degree below \(a/\lambda\).
   Allow arbitrary fixed coefficients and any heights at most \(X\).
   Then
   \[
                        D\le \max\{\max(W_U,W_V)/2,\ B\}+o(\lambda).
                                                                    \tag{FT2}
   \]
   In particular the \(B\) term vanishes when all free degrees are
   at least \(a/\lambda\).

Here the error means an arbitrary fixed positive multiple of
\(\lambda\) at sufficiently small budgets, plus an arbitrary
inner \(X\)-exponent allowance. Thresholds depend on \(C_0\),
the requested margin and the fixed height and cutoff margins in
the respective statement. The number of positions and their
degrees can depend on \(\lambda\). Fix them before the height
limit. The dictionary bounds concern individual original
coordinates.

### One fixed residue forces concentration

Work in the sequential limits, dividing all dictionary costs and
the count exponent by \(\lambda\). Denote the limits by \(w_U,w_V,d\).
Fix a block of positions of positive limiting numerical weight
on \(U\) having degrees at least \(\gamma/\lambda\) for some fixed
\(\gamma>0\). Its full product \(M\) divides \(U\); the other term
lies in the fixed unit residue \(\epsilon_V t\bmod M\). If the
block costs \(f\) in normalized units, direct multiplicative
character orthogonality and the singleton moments give
\[
                              d\le f+w_V/2 .                       \tag{FT3}
\]
The sum of the raw root caps on \(V\) is \(O(\lambda)\).
For a nonzero-cap coordinate use the \(2h_i\) moment with
\(h_i=\lceil(\sum_j\alpha_j)/\alpha_i\rceil\) on that side. Fix
its first \(h_i\) values, collapse the rest to a single product of
height \(X^{O(\lambda)}\), and use (FC2). The product has divisor-many
ordered decompositions, so the \(2h_i\) norm costs \(W_i/2\)
plus \(o(\lambda)/(2h_i)\). Hölder bounds the direct character
average by the product of these norms. Zero raw caps are frozen
at inner subpower cost. The value of \(V\) determines \(U\) at
the fixed target. The result does not involve the original
target support.
[Singleton congruence moments](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#diagonal-moments-for-a-joint-dictionary)

If \(d>\min(w_U,w_V)/2\), name the terms so that \(d>w_V/2\).
Every positive-weight high block on \(U\) has cost at least
\(\theta=d-w_V/2>0\). The following concentration argument uses
only the variable root height. Set aside the boundedly many
positions of cost at least \(\theta/4\). For any fixed
\(\gamma>0\), pack the remaining positions of degree at least
\(\gamma/\lambda\) into greedy blocks of cap \(\theta/2\) in
normalized units. Every closed block except possibly the last
costs more than \(\theta/4\), so their number is bounded. Each
has vanishing numerical weight by (FT3). All the remaining
degrees below \(\gamma/\lambda\) contribute at most \(C_0\gamma\).
Letting \(\gamma\) tend to zero, the weight outside the
set-aside positions vanishes. Pass to a subsequence fixing
their number and weight limits. Drop those of vanishing weight
into the same remainder. Thus \(U\), apart from its fixed
coefficient, consists of boundedly many narrow roots of height
\(X^{O(\lambda)}\) and a cofactor of numerical height exponent
\(o(1)\). This includes the empty set of narrow roots.

Now fix any positive-weight high block on \(V\). Apply (GW) to
the narrow roots and wider cofactor of \(U\) modulo that block,
with the fixed coefficient entering the unit target. The
count after fixing the block costs \(o(\lambda)\). For the
amplification, the exact unfixed product determines \(U\) and
then \(V\). Ordered decompositions of a product of a fixed
number of such integers are divisor-bounded. The same bound
covers representations of the single product by its roots and
wide factor. Fix the amplification order first and then send
\(\lambda\) to zero. An empty narrow list can be handled as
a single wide variable or by adding a root identically one.
Consequently every such block on \(V\) costs at least \(d\).
The concentration argument with this threshold applies to \(V\).
Interchanging the two sides and using (GW) again gives:
\[
 \text{every original coordinate of positive limiting numerical
 weight on either side has normalized cost at least }d.           \tag{FT4}
\]
Such a coordinate has degree at least a fixed positive multiple
of \(1/\lambda\), since its raw root height exponent is
\(O(\lambda)\). The multiplicative count is uniform in the
fixed coefficients, which are units at the opposite block.
[Wide-cofactor release](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

### Releasing two original roots

When both terms have positive variable numerical height, select
one positive-weight original coordinate on each term, with
degrees \(p,q\). Fix all the other positions. The resulting
normalized count cost of this record is at most \(w_U+w_V-2d\)
by (FT4). Any remaining pair solves
\[
                             A x^p+B y^q=t
\]
for fixed nonzero signed integer coefficients. This affine
plane curve is absolutely irreducible: over \(\mathbb C(x)\)
the polynomial in \(y\) is Eisenstein at any simple zero of
\(t-Ax^p\), and Gauss's lemma applies. The inputs \(x,y\) have
height at most \(X^{C_0\lambda+o_X(1)}\). With
\(\max(p,q)\gg1/\lambda\), the coefficient-uniform integral
plane-curve bound in the square gives count
\(O(X^{O(\lambda^2)+\varepsilon})\). At each fixed budget the
degrees are fixed for this application. Hence
\(d\le w_U+w_V-2d\).
[Coefficient-uniform plane-curve bound](https://arxiv.org/html/2608.24512v1#S2.SS1)

To finish (FT1), the nontrivial case is \(d>(w_U+w_V)/3\), which
implies \(d>\min(w_U,w_V)/2\). If one free term has vanishing
variable height, it also has vanishing complete height under the
coefficient hypothesis. The other term is tall by the shell
condition and has ratio to \(|t|\) within \(O(X^{-\tau})\) of
one for any fixed \(0<\tau<1\). Its narrow roots and wide
cofactor can be released using the real strip estimate (RUW1),
with the fixed coefficient absorbed into the centre. An exact
value determines the pair, with divisor-bounded multiplicity
from root representations. The chart parameter can be taken
arbitrarily small after fixing it before the budget limit.
This gives \(d=0\). When both variable heights are positive the
plane-curve count already contradicts the assumed failure.
[Unequal-height real strip](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

For (FT2), suppose \(d>\max\{\max(w_U,w_V)/2,b\}\), where \(b\)
is the limiting normalized low-degree cost. The concentration
argument applies. A side of vanishing variable height has raw
cost \(o(\lambda)\) in the positions of degree at least
\(a/\lambda\), and the other positions cost at most
\((b+o(1))\lambda\). Listing that complete value then contradicts
\(d>b\). Thus each side has a positive-weight original
coordinate. By (FT4) each has exactly one, because two such
positions would cost at least \(2d>\max(w_U,w_V)\). Every other
high-degree position combined has raw cost \(o(\lambda)\),
since its contribution to numerical height is at least
\(a/\lambda\) times its raw cap. Freeze all these complements
at normalized cost at most \(b+o(1)\) and release the two
roots by the same high-degree curve bound. This gives
\(d\le b\), proving (FT2). All numerical limits in this
argument are taken after the height limit at the corresponding
fixed budget.

## Sums of two highly full integers at every target

Let \(r_k(N)\) count ordered positive representations \(N=m_1+m_2\)
with both summands \(k\)-full, including representations with a
common divisor. For every \(\eta>0\), all sufficiently large fixed
integers \(k\), and every \(\varepsilon>0\),
\[
               r_k(N)\ll_{k,\eta,\varepsilon}
                       N^{(1/2+\eta)/k+\varepsilon}
                         \qquad (N\ge1).                            \tag{FT5}
\]

First consider coprime pairs at fixed coefficients, with the
variable integer part of each summand \(k\)-full. At
\(\lambda=1/k\), use the anatomic decomposition only on those
\(k\)-full parts. Choose its loss smaller than both
\(1/k\) and an arbitrarily small fixed multiple of \(\eta/k\).
In the proof of BBLT Lemma 2.2 all large exponents are packed
into a degree at least \(k\), by taking the truncation sufficiently
large so that its packing index is at least \(k\).
The exceptional coefficients are enumerated within the chosen
loss. The other variable degrees are the original degrees,
which are at least \(k\). Their raw cap sum on each side is at
most \(1/k+o_X(1)\) directly from the numerical height bound.
The packing construction has no overlap among distinct
positive-index positions. The original fixed coefficients,
even when of positive numerical weight, remain fixed and
are allowed in (FT2). Apply that bound with \(a=1\), absorbing
the dyadic and coefficient losses in \(\eta/k\) and \(\varepsilon\).
This gives (FT5) for these coprime coefficient-weighted sums,
uniformly in the coefficients and target.
[BBLT, Lemma 2.2 and its packing construction](https://arxiv.org/html/2410.12234v2#S2)

For an arbitrary representation, put \(h=\gcd(m_1,m_2)\mid N\)
and divide by \(h\). The coprime quotients have target \(N/h\).
Fix the parts of both quotients supported on the primes of \(h\).
These have \(O_\varepsilon(N^\varepsilon)\) possibilities,
with any smaller prescribed exponent available, by the uniform
fixed-support count. Such fixed parts are the coefficients of
the previous case. The primes outside \(h\) retain their
original exponents, so the remaining variable parts are
\(k\)-full. Sum the coprime bound over these choices and over
the divisor-many \(h\). This proves (FT5), with a constant
independent of \(N\).
[Arbitrary common divisor in the powerful sum problem](POWERFUL_SUM_GRAPHS.md)

For example, for any fixed \(\eta<1/2\) this improves the
leading exponent in the signed canonical-coordinate graph
bound \(1/(k+1)\) for sufficiently large \(k\). It gives an
arithmetic transfer from the small-budget modulus moments
to the powerful-sum problem with no condition on the
prime factorization of the integer being represented.
The all-target square-full problem corresponds to \(k=2\).
