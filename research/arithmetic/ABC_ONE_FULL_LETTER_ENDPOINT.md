# One powerful letter in a small-budget \(abc\) count

A positive integer is \(p\)-full when each positive prime valuation
is at least \(p\). Fix the position of one letter in an ordered
primitive positive triple \(a+b=c\le X\).

## Theorem. A saving above the single-letter threshold

For every \(\tau>0\) there are \(\delta_\tau,\lambda_\tau>0\)
such that for all fixed \(0<\lambda\le\lambda_\tau\),
integers \(p\) with \(p\lambda\ge3+\tau\), and
\(\varepsilon>0\), the number of such triples with the prescribed
letter \(p\)-full and \(\operatorname{rad}(abc)<c^\lambda\) is
\[
O_{\lambda,p,\tau,\varepsilon}
  \bigl(X^{(1/3-\delta_\tau)\lambda+\varepsilon}\bigr).
\]
The constants may be chosen for all three possible positions of the
prescribed letter.

For a contradiction, choose failing fixed-parameter counts with
\(\lambda_j\to0\) and count exponents at least
\((1/3-o(1))\lambda_j\). Take the height limit for each such
count before letting \(j\to\infty\). Dyadic splitting and the
anatomic factorization give boxes of root mass
\(L_j\le\lambda_j+o(\lambda_j)\) with count exponent
\(D_j\ge (1/3-o(1))\lambda_j\). The one-third box bound
forces \(L_j/\lambda_j\to1\) and \(D_j/L_j\to1/3\).
The truncation loss is \(o(\lambda_j)\). At each fixed \(j\)
choose the anatomic compression index at least \(p_j\), so each
power coordinate on the prescribed letter has degree at least
\(p_j\) whenever its mass is positive.
[The anatomic reduction and order of limits](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#theorem-a-strict-saving-from-endpoint-diffusion) · [The one-third box bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

The normalized root mass \(l_{T,j}\) of this letter tends to \(1/3\)
by endpoint balance. Its numerical height is at most \(1+o(1)\)
and is at least \(p_jL_j l_{T,j}\). Therefore
\(\limsup p_jL_j\le3\), contradicting
\(p_jL_j\ge3+\tau-o(1)\). [Endpoint balance](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)

## Corollary. The degree measures at the single-letter threshold

At any limiting \(abc\) box with normalized count exponent \(d=1/3\),
every letter's weak degree measure \(\nu_T\) satisfies
\[
 \nu_T([0,s))\ \ge\ \frac13-\frac1s\qquad(s>3).
\]
If one letter has \(\nu_T([0,3))=0\), that measure is
\(\frac13\delta_3\), the other two measures are wholly nonatomic
of mass \(1/3\), and the total nonatomic mass is \(2/3\).
In particular, the latter profile is necessary for a saturating
sequence of single-\(p_j\)-full counts with \(p_j\lambda_j\to3\).

The balance and first-moment bounds give total mass \(1/3\) and
\(\int r\,d\nu_T(r)\le1\). Thus \(\nu_T([s,\infty))\le1/s\).
A measure supported on \([3,\infty)\) can satisfy this first-moment
bound only by placing all its mass at three. The endpoint diffusion
theorem says that a letter with zero nonatomic mass forces both
other measures to be wholly nonatomic. For the saturating
single-full sequence, the argument above gives \(L_j/\lambda_j\to1\),
so its \(p_jL_j\to3\) gives the required support.
[Endpoint measures and the zero-nonatomic-mass case](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
