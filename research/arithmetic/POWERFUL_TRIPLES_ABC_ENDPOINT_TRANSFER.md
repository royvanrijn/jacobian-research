# High-powerful triples from the small-budget \(abc\) bounds

Let \(T_{p,q,r}(X)\) count ordered primitive positive triples
\(a+b=c\le X\) for which \(a,b,c\) are respectively \(p\)-full,
\(q\)-full and \(r\)-full. A positive integer is \(s\)-full if each
nonzero prime valuation is at least \(s\). Write
\[
\lambda_0=\frac1p+\frac1q+\frac1r.
\]

## Corollary. Transfer at the limiting one-third coefficient

For each fixed \(\eta>0\), every sufficiently small fixed \(\lambda_0\)
and each \(\varepsilon>0\),
\[
T_{p,q,r}(X)\ll_{p,q,r,\eta,\varepsilon}
   X^{(1/3+\eta)\lambda_0+\varepsilon}.
\]
Indeed, \(\operatorname{rad}(a)\le a^{1/p}\), with the analogous
inequalities for the other terms. Since \(a,b<c\),
\[
\operatorname{rad}(abc)
\le a^{1/p}b^{1/q}c^{1/r}<c^{\lambda_0}.
\]
Apply the cumulative \(abc\) count at \(\lambda_0\). At fixed offsets
\(p=r+u,\ q=r+v\), the exponent before the arbitrary \(\varepsilon\)
is at most \((1+o(1))/r\) as the base order tends to infinity.
[Powerfulness and the radical budget](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#radical-and-powerfulness-scopes) · [The limiting \(abc\) count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

## Theorem. A strict saving from endpoint diffusion

There are absolute constants \(\delta,\lambda_*>0\) such that for all
fixed \(p,q,r\ge2\) with \(\lambda_0\le\lambda_*\) and each
\(\varepsilon>0\),
\[
T_{p,q,r}(X)\ll_{p,q,r,\varepsilon}
   X^{(1/3-\delta)\lambda_0+\varepsilon}.
\]
In particular, for \(p=r+u,\ q=r+v\) and fixed
\(u\ge v\ge0\), the exponent is at most \((1-3\delta)/r+\varepsilon\)
for all sufficiently large \(r\). The constants \(\delta,\lambda_*\)
come from a limiting argument.

Suppose the assertion fails for arbitrarily small positive \(\delta\)
and \(\lambda_*\). This gives a sequence of orders
\(s_{A,j},s_{B,j},s_{C,j}\), with \(\lambda_{0,j}=\sum_T1/s_{T,j}\to0\),
for which there are arbitrarily large heights attaining a count
exponent at least \((1/3-o(1))\lambda_{0,j}\). Take the height limit at
each fixed triple of orders, and take the limit \(j\to\infty\) after
forming its box parameters.

Use the anatomic factorization
\(n=n_0\prod_{e=1}^{d_0}n_e^e\) with truncation loss
\(\varepsilon_{0,j}=o(\lambda_{0,j})\). In the construction of
Lemma 2.2 of Bernert, Browning, Lichtman and Teräväinen, valuations
through \(d_0\) stay at their own index and valuations above \(d_0\)
go into a power at the compression index \(k\), with a remainder
in \(n_0\). At each fixed \(j\) choose \(\varepsilon_{0,j}\) small
enough that \(k\ge\max_T s_{T,j}\). In particular, for a letter of
order \(s_{T,j}\), each power coordinate below \(s_{T,j}\) is \(1\).
The exceptional coefficients are fixed at a total count loss
\(o(\lambda_{0,j})\) in the \(X\)-exponent.

Dyadically split the heights, radicals and power coordinates and retain
a box supporting the claimed exponent up to \(o(\lambda_{0,j})\).
Write \(L_j\) for its total logarithmic root mass, \(L_{T,e,j}\) for
a coordinate's mass and \(D_j\) for its count exponent, with heights
measured relative to the shell of \(c\). The radical comparison in
Lemma 2.2 and \(\operatorname{rad}(abc)<c^{\lambda_{0,j}}\) give
\(L_j\le\lambda_{0,j}+o(\lambda_{0,j})\). The box has
\(D_j\ge (1/3-o(1))\lambda_{0,j}\), and the one-third box bound gives
\(\limsup D_j/L_j\le1/3\). These inequalities force
\(L_j/\lambda_{0,j}\to1\) and \(D_j/L_j\to1/3\). Choices of coefficients
and boxes have already been paid for at fixed \(j\) before taking
these limits. [Anatomic reduction, Proposition 2.1 and Lemma 2.2](sources/bernert-browning-lichtman-teravainen-2410.12234v2.pdf) · [The limiting \(abc\) box bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

By the endpoint diffusion result, each letter's normalized root mass
\(l_{T,j}=\sum_e L_{T,e,j}/L_j\) tends to \(1/3\). The weak limits of
\[
 \nu_{T,j}=\sum_e\frac{L_{T,e,j}}{L_j}\delta_{eL_j}
\]
must have total nonatomic mass at least \(1/2\) when summed over the
three letters. Here \(\delta_x\) is unit mass at \(x\). Since all
nonzero coordinates of \(T\) have \(e\ge s_{T,j}\), its total numerical
height bounds \(s_{T,j}L_jl_{T,j}\) by \(1+o(1)\). Thus
\(\limsup s_{T,j}\lambda_{0,j}\le3\) for each letter. The three numbers
\(1/(s_{T,j}\lambda_{0,j})\) sum to one and each has limit inferior
at least \(1/3\). Each tends to \(1/3\), hence \(s_{T,j}L_j\to3\)
for all three letters.

Every weak limit is therefore supported on \([3,\infty)\). It has
mass \(1/3\) and first moment at most one, by the height bound and
lower semicontinuity. The moment is already at least one on that
support. Equality forces \(\nu_T=\frac13\delta_3\) on every letter.
The total nonatomic mass is zero, contrary to the endpoint
condition. This contradiction gives absolute \(\delta,\lambda_*>0\).
[The endpoint diffusion theorem and its measures](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## Corollary. A restricted fixed-target bound

Fix \(p,q,r\) as in the theorem. For every positive \(r\)-full \(N\),
the number of ordered coprime positive pairs \(a+b=N\) with \(a\)
\(p\)-full and \(b\) \(q\)-full is
\(O_{p,q,r,\varepsilon}(N^{(1/3-\delta)\lambda_0+\varepsilon})\).
Each pair is counted by \(T_{p,q,r}(N)\). The proof uses the
fullness of \(N\) and the coprimality of the summands. The graph
bounds for general prescribed targets use different hypotheses.
[Fixed-target graph bounds](POWERFUL_SUM_GRAPHS.md)

## Corollary. The earlier fixed-offset count threshold

For fixed \(u\ge v\ge0\), the strict saving theorem eventually implies
\[
T_{r+u,r+v,r}(X)\ll_{u,v,r,\varepsilon}
 X^{1/(r+u)+1/(r+v)-2/r^2+\varepsilon}.
\]
Indeed, the exponent preceding \(\varepsilon\) on the right is
\(2/r+O_{u,v}(r^{-2})\), while the strict saving gives at most
\((1-3\delta)/r\). This also gives the count with an allowed
\(O_{u,v}(r^{-5/2})\) error in the exponent. In the aligned
common-neighbor argument, improving the \(O_{u,v}(r^{-11/5})\)
error for unequal offsets is a separate estimate in that
method. [The aligned graph estimate](POWERFUL_MONOMIAL_TAYLOR_REFINEMENT.md)

## Corollary. Low scaled degrees in at least two letters

At any limiting \(abc\) box with normalized count exponent \(d=1/3\),
at least two of the three weak degree measures give positive root mass
to \([0,3)\). If exactly two do, those two measures are wholly
nonatomic and the remaining measure is \(\frac13\delta_3\).

Indeed, every letter measure has mass \(1/3\) and first moment at most
one. A measure with zero mass below three is then
\(\frac13\delta_3\) by the same height argument as in the strict
saving theorem. The endpoint diffusion theorem implies that when
one letter has zero nonatomic mass, the other two letters have
nonatomic mass \(1/3\) each. Neither can equal \(\frac13\delta_3\),
so each must have positive mass below three.
[The endpoint diffusion theorem](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

## Corollary. A small-budget saving with two full letters

There are absolute constants \(\tau,\delta,\lambda_*>0\)
with the following property. For fixed \(0<\lambda\le\lambda_*\),
impose \(p\)-fullness and \(q\)-fullness on any prescribed two
letters, with
\[
 \min(p,q)\lambda\ge3-\tau.
\]
The number of ordered positive primitive \(abc\) triples
with these two fullness conditions, \(c\le X\), and
\(\operatorname{rad}(abc)<c^\lambda\) is
\[
 O_{\lambda,p,q,\varepsilon}
 \bigl(X^{(1/3-\delta)\lambda+\varepsilon}\bigr).
\]
The third letter can have arbitrary valuation degrees.

If the bound failed for every choice of the three positive
constants, take a sequence with \(\lambda_j\to0\),
\(\min(p_j,q_j)\lambda_j\ge3-o(1)\) and count
exponent at least \((1/3-o(1))\lambda_j\). Fix the
positions of the two full letters on a subsequence.
Form the anatomic boxes with truncation loss
\(o(\lambda_j)\), as in the three-letter proof.
The radical budget gives \(L_j\le(1+o(1))\lambda_j\),
and the limiting one-third box bound gives
\(\limsup D_j/L_j\le1/3\). Thus
\(L_j/\lambda_j\to1\) and \(D_j/L_j\to1/3\).
The two full letters have no roots at indices below
their respective orders, taking the compression index
above both orders. Since \(p_jL_j,q_jL_j\ge3-o(1)\),
their two weak degree measures are supported on
\([3,\infty)\). The preceding low-degree corollary
requires positive mass below three on at least two
letters, giving the contradiction.

The same limiting argument applies to boxes in which
each of two letters has normalized root mass at most
\(\tau\) below scaled degree \(3-\tau\). For otherwise
a diagonal sequence with \(\tau\to0\) and \(d\to1/3\)
would again give two measures supported on
\([3,\infty)\). In particular the gap for boxes does
not require exact absence of the low degree coordinates.
