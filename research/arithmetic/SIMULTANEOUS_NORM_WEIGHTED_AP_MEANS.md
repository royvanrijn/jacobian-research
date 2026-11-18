# Weighted simultaneous norms in arithmetic progressions

The simultaneous Gaussian--Eisenstein weight has a uniform one-variable
mean in reduced arithmetic progressions. Its twisted Dirichlet series is a
product of four Dirichlet \(L\)-functions with integer exponents, divided by
an Euler product analytic to the right of \(1/2\). An elementary hyperbola
partition therefore gives a power-saving error for nonexceptional twists.
This is a distribution estimate for one weight, not a pointwise binary
correlation theorem.

## Statement and normalization

Let \(h(n)=r_G(n)r_E(n)/24\), where \(r_G,r_E\) count all ordered integer
representations by \(x^2+y^2\) and \(u^2-uv+v^2\). Put

\[
 h_6(n)=h(n)\mathbf1_{(n,6)=1},\qquad
 D_p=\sum_{e\ge0}\frac{h(p^e)}{p^e},\qquad
 E(q)=\prod_{\substack{p\mid q\\p\ge5}}D_p^{-1}.
\]

Write
\(\rho_6=\log(2+\sqrt3)/(3\sqrt3)>0\), the residue of the Dirichlet
series of \(h_6\) at one. All implied constants below concern this fixed
weight; in particular they are independent of the modulus and residue.

**Weighted AP theorem.** Let \(q\ge1\), \((a,q)=1\), and \(d=(q,12)\).
If \(a\not\equiv1\pmod d\), the AP sum is identically zero. Otherwise,
uniformly for \(x\ge1\),

\[
 \sum_{\substack{n\le x\\n\equiv a\ (q)}}h_6(n)
 =\rho_6\frac{\varphi(d)}{\varphi(q)}E(q)x
   +O\!\left(qx^{3/4}(\log(2x))^2\right).
 \tag{1}
\]

In particular, if \((q,12)=1\), every fixed \(A,B>0\) gives

\[
 \sum_{\substack{n\le x\\n\equiv a\ (q)}}h_6(n)
 -\frac1{\varphi(q)}\sum_{\substack{n\le x\\(n,q)=1}}h_6(n)
 \ll_{A,B}\frac{x}{(\log x)^A}
 \quad\left(q\le(\log x)^B\right).
 \tag{2}
\]

When the modulus shares factors with 12, the compatible-class factor in
(1) must be retained. A principal-character-only formula would be false.

## The twisted Dirichlet series

Let \(\mathcal X_{12}\) be the four characters **modulo 12**, all extended
by zero on nonunits, including its principal character. For a character
\(\chi\pmod q\), Euler-factor comparison gives

\[
 F_\chi(s):=\sum_{n\ge1}\frac{h_6(n)\chi(n)}{n^s}
 =\frac{\displaystyle\prod_{\psi\in\mathcal X_{12}}L(s,\chi\psi)}
        {L(2s,\chi^2\chi_{12})},\qquad\Re s>1.
 \tag{3}
\]

Here \(\chi_{12}\) is the quadratic character of discriminant 12.
Each product \(\chi\psi\) is the actual, possibly imprimitive, character
modulo \(\operatorname{lcm}(q,12)\). At 2 and 3 all Euler factors in
(3) are one. At every other prime, substitute \(z=\chi(p)p^{-s}\) into

\[
 \sum_{e\ge0}h(p^e)z^e=
 \begin{cases}
 (1+z)/(1-z)^3,&p\equiv1\ (12),\\
 (1+z^2)/(1-z^2)^2,&p\equiv5,7\ (12),\\
 (1-z^2)^{-1},&p\equiv11\ (12).
 \end{cases}
 \tag{4}
\]

This proves (3), including primes dividing the conductor of \(\chi\).
If primitive representatives \(1,\chi_{-4},\chi_{-3},\chi_{12}\) are used
instead of characters modulo 12, the numerator requires the extra factors
\((1-\chi(2)^22^{-2s})(1-\chi(3)^23^{-2s})\).

The reciprocal denominator in (3) is analytic for \(\Re s>1/2\), bounded
uniformly on \(\Re s\ge1/2+\varepsilon\) by \(\zeta(1+2\varepsilon)\).
Its square-divisor coefficients are

\[
 \gamma_\chi(m)=\mu(m)\chi(m)^2\chi_{12}(m),\qquad
 |\gamma_\chi(m)|\le1.
 \tag{5}
\]

All four numerator exponents are the positive integer one. Their zeros
create no fractional-power branch issue.

## Nonexceptional twists: an elementary power saving

Call \(\chi\) exceptional if its underlying primitive character is one of
\(1,\chi_{-3},\chi_{-4},\chi_{12}\). Otherwise every \(\chi\psi\) in
(3) is nonprincipal. Each has period
\(P=\operatorname{lcm}(q,12)\le12q\), absolute value at most one, and
zero complete-period sum. If \(a_\chi(n)\) denotes the coefficients of
the numerator product, then

\[
 \sum_{n\le x}a_\chi(n)\ll qx^{3/4}(\log(2x))^2.
 \tag{6}
\]

To verify this directly, set \(y=x^{1/4}\) in the fourfold product sum.
The cube with every coordinate at most \(y\) has product at most \(x\).
Its character sum factors; bound one factor by \(P\) and the other three
by \(y\), giving \(Py^3\). Outside the cube, partition by the first
coordinate larger than \(y\). Fixing the other three coordinates, with
product \(m\), leaves a character sum over \(y<n\le x/m\), of absolute
value at most \(2P\). Such a sum requires \(m<x/y=x^{3/4}\). The number
of ordered triples with product at most \(t\) is
\(O(t(\log(2t))^2)\). Summing the four pieces proves (6).
Convolution with (5) now gives

\[
 \sum_{n\le x}h_6(n)\chi(n)
 =\sum_{m\le\sqrt x}\gamma_\chi(m)
                   \sum_{n\le x/m^2}a_\chi(n)
 \ll qx^{3/4}(\log(2x))^2,
 \tag{7}
\]

since \(\sum m^{-3/2}\) converges. This argument uses character periodicity,
not a zero-free region.

## Exceptional twists and the AP main term

The [case study](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#an-unconditional-global-count-averaged-over-targets)
proves

\[
 H_6(x):=\sum_{n\le x}h_6(n)=\rho_6x+O(x^{3/4}).
 \tag{8}
\]

The imported input is [Debaene, arXiv:1611.10103v2, Corollary 2](https://arxiv.org/html/1611.10103v2):
for a fixed degree-four field, ideal counting has main term
\(\operatorname{Res}_{s=1}\zeta_K(s)x\) and error \(O_K(x^{3/4})\).
The case study applies it only to \(K=\mathbb Q(\zeta_{12})\), transfers
through the absolute square-divisor convolution, and removes the 2 and 3
Euler factors. The latter multiply the residue by
\((1-1/4)(1-1/9)=2/3\), giving the stated \(\rho_6\).

Every supported value of \(h_6\) is \(1\pmod{12}\), so an exceptional
character satisfies
\(h_6(n)\chi(n)=h_6(n)\mathbf1_{(n,q)=1}\). Finite Euler inversion in
(8) gives, uniformly in \(q,x\),

\[
 \sum_{n\le x}h_6(n)\chi(n)=\rho_6E(q)x+O(qx^{3/4}).
 \tag{9}
\]

Here is uniform control of the inversion. Its local inverse factors from
(4) are \((1-z)^3/(1+z)\), \((1-z^2)^2/(1+z^2)\), and \(1-z^2\);
their coefficients have absolute value at most eight. If \(\delta_q(m)\)
are the coefficients of the finite inverse product, then

\[
 \sum_{m\ge1}\frac{|\delta_q(m)|}{m^{3/4}}
 \le\prod_{p\mid q,\ p\ge5}\left(1+\frac8{p^{3/4}-1}\right)
 \le\prod_{p\mid q,\ p\ge5}p\le q.
 \tag{10}
\]

Indeed \(p^{3/4}>3\) for \(p\ge5\), so each parenthesis is less than
five. Convolving (8) gives the asserted error; extending its main coefficient
sum beyond \(m=x\) also costs at most \(qx^{3/4}\).

The exceptional characters modulo \(q\) are exactly those induced from
characters modulo \(d=(q,12)\). Their number is \(\varphi(d)\), and their
sum at \(a\) is \(\varphi(d)\) when \(a\equiv1\pmod d\), and zero
otherwise. Their twisted sums (9) are identical. Character orthogonality,
(7) and (9) prove (1). The incompatible classes are identically empty
because the support of \(h_6\) is \(1\pmod{12}\). When \((q,12)=1\)
only the principal character is exceptional, which proves (2).

## Short intervals and the pointwise limitation

Subtract (1) at \(X+H\) and \(X\), for \(0<H\le X\). In every
compatible reduced class,

\[
 \sum_{\substack{X<n\le X+H\\n\equiv a\ (q)}}h_6(n)
 =\rho_6\frac{\varphi((q,12))}{\varphi(q)}E(q)H
   +O\!\left(qX^{3/4}(\log(2X))^2\right).
 \tag{11}
\]

Its leading coefficient is at least
\(c/[q(\log\log(3q))^4]\). Indeed \(h(p^e)\le d_4(p^e)\), where
\(d_4\) is the four-divisor function, so \(D_p\le(1-1/p)^{-4}\).
The [prime harmonic estimate of Loughran--Matthiesen, Lemma 2.4(2)](https://arxiv.org/html/1904.12845v2#S2),
splitting prime divisors at \(\log(3q)\), gives
\(E(q)\gg(\log\log(3q))^{-4}\). Thus for fixed
\(0<\theta<1\), \(\kappa\ge0\) with \(\theta>3/4+2\kappa\),
(11) is a uniform positive asymptotic when \(H=X^\theta\) and
\(q\le X^\kappa\). In particular the full-length weighted AP level
proved here includes every fixed exponent below \(1/8\).

These estimates are uniform in a growing reduced residue \(a\). They do
not control the second weight in
\(\sum_nh_6(n)h_6(N-n)\): that reflected weight is not an AP indicator.
The [natural pointwise lower bound and fixed-prime correlations](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md)
therefore remain open on every target. The separate
[full weighted correlation theorem](SIMULTANEOUS_NORM_POWER_CORRELATIONS.md)
now proves them outside a power exceptional set, using a further Fourier
argument and exact local-model comparison. Roughness support and the parameter uniformity of any
separate bilinear-distribution theorem also require their own proofs.
The existing unweighted variance alone supplies its own count scale, rather
than the distribution of coordinate weights across every target's pairs.

Debaene's actual Corollary 2 and the residue normalization were checked.
The prime harmonic estimate used for the uniform lower coefficient was also
checked in the cited primary lemma. The remaining arguments use elementary
character periodicity, convolution and counting. No novelty determination, external mathematical review or formal
verification is claimed; no finite coverage search is used.
