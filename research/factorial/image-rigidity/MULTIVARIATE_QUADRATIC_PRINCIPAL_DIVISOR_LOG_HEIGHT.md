# A prime-support logarithmic height bound for a quadratic factor

A common irreducible quadratic of two ordinary factorial moments
has a root \(\alpha\) in an imaginary quadratic field \(K\).
Its direction \(\gamma=\alpha/\overline\alpha\) is a norm-one
element of \(K\). The valuations of \(\gamma\) determine the
reduced denominator \(D\) of \(\gamma+\gamma^{-1}\).
The number and size of the split primes in \(D\), and their
ideal classes in \(K\), give an effective quantitative
restriction near an algebraic diffuse contact. The estimate
does not bound the row counts, heights, gaps, or retained core
size. It applies to moving fields and moving split primes.

## The bound in terms of the split primes

Let \(\gamma_*\) be a fixed algebraic number of modulus one.
Suppose \(\gamma\ne\gamma_*\). For each \(p\mid D\) choose
one of the primes \(\mathfrak p_p\) of \(K\) above \(p\), and set
\[
        r=\#\{p:p\mid D\},\qquad
        \kappa=\left|\left\langle[\mathfrak p_p]^2:
                     p\mid D\right\rangle\right|.
\]
The choice of the prime above \(p\) does not change \(\kappa\).
With an effective constant depending only on \(r,\gamma_*\),
\[
 \log^+\frac{1}{|\gamma-\gamma_*|}
 \ \ll_{r,\gamma_*}\
       \kappa \prod_{p\mid D}\log(3p)\,
                         \log\!\bigl(2+\kappa\log(3D)\bigr).
                                                               \tag{1}
\]
If \(D=1\), the direction belongs to the fixed finite set
of roots of unity in imaginary quadratic fields. In
particular the nonattained distances are bounded below
in this case as well. Since \(\kappa\le h_K\), one may
use the class number in place of the indicated subgroup size.

Write \(s_p,t_p\) for the two root valuations at \(p\)
in the original common variable. If they differ, \(p\)
splits in \(K\), and \(v_p(D)=|s_p-t_p|\) by the
[reciprocal trace identity](MULTIVARIATE_PRIMITIVE_QUADRATIC_NORM_GATE.md#primitive-scale-and-the-norm-equation).
At equal valuations the quotient is a unit. Thus every
prime in (1) is an actual prime of spread, and its positive
exponent can be computed from a given Newton allocation.
Combining (1) with a local allocation test means using the
same allowed allocation for the support and exponent
estimates. The denominator may be computed in any common
rational argument variable.

At a simple transverse diffuse reference, the actual
factor direction satisfies \(|\gamma-\gamma_*|=O(\eta)\)
for the maximum added relative gap \(\eta\), in the
notation of the
[directional height comparison](MULTIVARIATE_TRANSVERSE_QUADRATIC_DIRECTIONAL_HEIGHT.md#the-exact-first-order-direction).
Any uniform bound \(|\gamma-\gamma_*|=O(\eta^a)\), \(a>0\),
can be used instead. It makes the right side of (1)
at least a fixed positive multiple of \(\log(1/\eta)\)
up to a bounded error. For a fixed finite set of possible
fields and a fixed finite set of possible spread primes,
(1) then gives the stronger necessary growth
\[
                         \log(3D)\gg \eta^{-b}
                         \quad\text{for some }b>0 .            \tag{2}
\]
This excludes a nonattained contact if all allowed spreads
keep \(D\) at most a fixed power of \(1/\eta\).
For bounded \(r\) and moving primes, the right side of (1)
is a different test. It can still be too small when one
deep prime accounts for most of \(\log D\). By contrast,
a spread supported on a prime \(p\) growing as a fixed
power of \(1/\eta\) can already pay a logarithmic cost
through \(\log p\). If \(r,\kappa\) are bounded and all allowed
denominators satisfy \(D\le \eta^{-C}\) for a fixed \(C\), at least
one spread prime must satisfy
\[
 \log(3p)\gg
 \left(\frac{\log(1/\eta)}{\log\log(1/\eta)}\right)^{1/r}
\]
near a nonattained reference. Indeed the logarithm on the right
of (1) is \(O(\log\log(1/\eta))\) under these bounds, and the
product has at most \(r\) factors. This excludes support made
only of polylogarithmically growing primes in this regime.
The estimate gives no bound of the
required kind on the actual prime allocations for all
unrestricted profiles. Exact equality \(\gamma=\gamma_*\)
must also be considered separately. It is impossible
if \(\gamma_*\) belongs to none of the allowed fields.

## Simultaneous principalization with one class-group factor

Use the signed exponents \(e_p\) of
\(\mathfrak p_p/\overline{\mathfrak p_p}\) in the divisor
of \(\gamma\). Thus \(|e_p|=v_p(D)\).
The exponent vector \(e\) lies in the full lattice
\[
       L=\left\{v\in\mathbb Z^r:
                \prod_{p\mid D}[\mathfrak p_p]^{2v_p}=1
                    \text{ in }\operatorname{Cl}(K)\right\}.
\]
Its index is \(\kappa\), since complex conjugation
inverts the split prime class. Set \(w_p=\log(3p)\) and
\(\|v\|_w=\max_p w_p|v_p|\). After applying the diagonal
map with entries \(w_p\), the lattice has determinant
\(\kappa\prod_p w_p\). Minkowski's second theorem applied
to the unit cube gives independent vectors \(v_1,\ldots,v_r\)
with
\[
                     \prod_i\|v_i\|_w\le\kappa\prod_p w_p .     \tag{3}
\]
For each \(v_i\) choose a generator \(\lambda_i\) of the
corresponding principal fractional ideal. Its norm is one,
and its absolute logarithmic Weil height is
\(\frac12\sum_p |v_{i,p}|\log p\). Hence
\[
           \max\{1,h(\lambda_i),|\log\lambda_i|\}
                          \ll_r \|v_i\|_w
\]
for the principal complex logarithm. This uses
\(|\lambda_i|=1\) and \(\|v_i\|_w\ge\log 6\).

Let \(V\) have columns \(v_i\), and put \(d=|\det V|\).
Determinant expansion and (3) show \(d\le r!\kappa\).
Cramer's rule gives integers \(b_i\) with \(Vb=d e\),
and replacing column \(i\) by \(e\) in the weighted
determinant gives
\[
          |b_i|\le
            \frac{r!\|e\|_w\prod_{j\ne i}\|v_j\|_w}
                 {\prod_p w_p}
          \le\frac{r!\kappa\|e\|_w}{\|v_i\|_w}.
\]
In particular \( |b_i|\ll_r\kappa\log(3D)\), since
\(\log(3p)\le(\log 6/\log2)\log p\).
The divisor of \(\gamma^d/\prod_i\lambda_i^{b_i}\)
is zero. It is a unit of \(K\), and its twelfth power
is one, since imaginary quadratic unit orders divide
twelve. With \(A=12d\) this gives
\[
                          \gamma^A=\prod_i\lambda_i^{12b_i}.
                                                               \tag{4}
\]
Using the lattice rather than treating each prime ideal
independently is what gives the single factor \(\kappa\)
in (3).

If \(|\gamma/\gamma_*-1|>1/(4A)\), the claimed logarithmic
bound follows from \(A\ll_r\kappa\). Otherwise, put
\(x=\gamma/\gamma_*\). A nontrivial \(A\)-th root of unity
is at distance at least \(4/A\) from one for \(A\ge2\),
so \(x^A\ne1\). Moreover
\[
            |x^A-1|\le A|x-1|(1+|x-1|)^{A-1}
                                  <\tfrac12 .
\]
The principal logarithm of \(x^A\) is nonzero, is at
most \(2|x^A-1|\) in modulus, and by (4) is a linear form
\[
        \sum_i 12b_i\log\lambda_i-A\log\gamma_*
                                        -2k\log(-1)
\]
for some integer \(k\).
The imaginary parts of the principal logs yield
\(|2k|\le1+12\sum_i|b_i|+A\).
The field of the logarithms has degree at most
\(2[\mathbb Q(\gamma_*):\mathbb Q]\), independently of \(K\).
Apply the complex form of
[Matveev's Corollary 2.3](https://www.mathnet.ru/links/11696330cc866d273c92c771a63da235/im314_eng.pdf),
which does not require logarithmic independence. Its product
of height parameters is \(O_{r,\gamma_*}(\kappa\prod_p w_p)\)
by (3), and its logarithm of the coefficient bound is
\(O_r(\log(2+\kappa\log(3D)))\).
Dividing the resulting lower bound by \(2e^{1/4}A\)
and absorbing \(\log A\) proves (1). The degree and
height hypotheses here depend on the reference target,
not on the size of the factorial parents.
[Minkowski's second theorem, Theorem 3.3](https://people.maths.ox.ac.uk/greenbj/papers/addcomb2009-3.pdf)

## The scaled root and norm use the same argument

One can instead take \(z=q\alpha\) for any positive rational
scale \(q\) in the physical comparison. Let \(S\) be the set
of prime ideals of \(K\) supporting the fractional principal
ideal \((z)\), put \(s=|S|\), and write
\[
 \begin{split}
    E &=\sum_{\mathfrak p\in S}
        |v_{\mathfrak p}(z)|\log(3N\mathfrak p),\\
    \kappa_S&=\left|\left\langle[\mathfrak p]:
                                 \mathfrak p\in S\right\rangle\right|.
 \end{split}
\]
Here \(v_{\mathfrak p}\) is the integer ideal valuation.
For fixed nonzero algebraic \(z_*\) and \(z\ne z_*\),
\[
 \log^+\frac1{|z-z_*|}
 \ \ll_{s,z_*}\
  \kappa_S\prod_{\mathfrak p\in S}\log(3N\mathfrak p)\,
                       \log\!\bigl(2+\kappa_S(1+E)\bigr).
                                                               \tag{5}
\]
Indeed use the principal exponent lattice on the prime ideals
themselves. Its index is \(\kappa_S\). For a generator
with divisor vector \(v\), the two Archimedean conjugates
have the same modulus, whose logarithm is
\(\frac12\sum_{\mathfrak p}v_{\mathfrak p}\log N\mathfrak p\).
Its height and absolute logarithm are therefore
\(O_s(\|v\|_w)\) with weights \(w_{\mathfrak p}=\log(3N\mathfrak p)\).
The rest of the proof, including the comparison of powers
for \(x=z/z_*\) near one, is unchanged.
If \(S\) is empty, \(z\) is again in the fixed finite
unit set.

For \(q\) built from exact factorial ratios as in the
[common rational chart](MULTIVARIATE_SIMULTANEOUS_CONTACT_HEIGHT.md),
the prime support and \(E\) for a proposed factor follow
from the shifted root valuations and its quadratic field.
At a nonsplit prime, multiply the extended \(p\)-adic root
valuation by the local ramification index to obtain the
integer ideal valuation. At a split prime use the two
shifted valuations directly. Thus (5) can use absolute
norm restrictions as well as spread. It can exclude exact
angular attainment when the scaled target \(z_*\) is
outside the possible quadratic fields. When choosing \(q\)
only to make \(q\alpha\) primitive in \(\mathcal O_K\), its
norm is \(DR\) from the
[primitive norm equation](MULTIVARIATE_PRIMITIVE_QUADRATIC_NORM_GATE.md#primitive-scale-and-the-norm-equation).
The corresponding support then consists of one prime ideal
over each \(p\mid D\) and one over each ramified \(p\mid R\).
A comparison to a fixed scaled target must use a scale
with a known real limit. Primitivity alone does not supply
that limit.
