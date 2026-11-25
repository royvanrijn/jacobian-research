# Extended near-flat separation at terminal six for shared support

Let all unequal exponent rows of a distinct ordinary two-monomial profile
be shared, so both exponents in each such row are positive. Let \(u\) be
the order-five normalized outer geometric mean. If \(u\ge9/10\), the
order-five and order-six polynomials have disjoint root arguments after
independent positive rescalings, so their specialized rational Bézout
identity descends to every characteristic-zero field. The dimension is
arbitrary.

Use the scaled Cayley transform, notation, and exact polynomial expansions
(2) and (3) of the
[original near-flat proof](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md).
The [extension checker](verify_multivariate_terminal_six_near_flat_extension.py)
uses exact rational interval arithmetic. Replay from the repository root
with

```sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_six_near_flat_extension.py
```

## Linked coefficient bounds

Take \(0\le\lambda\le19/60\), \(p=\lambda^2\), and \(q=r p\), and use
the rational caps \(p<101/1000\), \(q<107/1000\).
For actual ordinary fully shared profiles, \(r\) lies between \(517/500\)
and \(53/50\). The lower bound uses the same inequality as the original
near-flat proof,
\[
 r \ge 1+\frac{u(1+u)}{48(1+(1-u^2)/24)^{3/2}},
\]
with \(u\ge1-p_{\max}\) and \(1+(1-u^2)/24\le1+p_{\max}/12\).
For the upper bound \(x_5\le((1-p_{\max})^{-2}-1)/4<1/16\), hence
\(d/a\le 5x_5/(1-x_5)\le1/3\) for the maximal reciprocal in its row.
Thus the order ratio is at most \(20/19<53/50\), and the Bernoulli
comparison in the initial proof applies. Every actual \(u\ge9/10\)
satisfies \(\lambda\le19/60\).

As in the base note, take \(x_5^*=(u^{-2}-1)/4\) and
\(x_6^*=(w^{-2}-1)/5\). The functions
\((2-t)/(4(1-t)^2)\) and \((2-t)/(5(1-t)^2)\) increase.
Use the endpoint caps \(x_5^*\le.588p,\ x_6^*\le.475q\).
Then the absolute skews are bounded by
\[
 |\alpha|\le.346p^2,\quad |\chi|\le.174p^2,\quad
 |\alpha-2\chi|\le.072p^4,\qquad
 |\beta|\le.377q^2,\quad|\eta|\le.303q^2,\quad |4\beta-5\eta|\le.306q^4.
\]
This is the same derivation by local log pair ratios and linked bounds
as in the original proof. The mean defects obey \(v-u^{3/2}\le.23p^3\),
\(y-w^{8/5}\le.282q^3\). Their proof again uses
\(G_{1/2}(s)\le s^3/6\), \(G_{3/5}(s)\le s^3/5\),
\(-\log(1-t)\le t/(1-t)\), and \(e^t\le(1-t)^{-1}\).
The central log defect \(J=-\log z+\tfrac32\log y-\tfrac35\log w\) satisfies
\(0\le J\le (6/5)(x_6^*)^5\le .03q^5\).

Substitute these enlarged bounds into the stable coefficient
expressions of the base note. In its notation (4) use
\[
\begin{array}{c|c@{\quad}c|c}
 A_5&[0,.06]& B_5&[0,.016]\\
 R_5&[3/8,3/8+.304p]&H&[0,.282]\\
 A_6&[0,.072]&B_6&[0,.046]\\
 R_y&[.064,.064+.03q]&R_z&[.048,.048+.02q]\\
 D&[-.098,.098]&&
\end{array}
\]
The pure cubic remainder for \(u^{3/2}\) is at most \(.074p^3\).
The fourth-order pure bounds for \(w^{8/5},w^{9/5}\) are from
their fourth derivatives on \(w\ge.893\). The central difference
uses exactly the same linked expansion as the original near-flat proof,
now with \(w^{-7/5}\le1.18\), \(w^{-3/5}\le1.08\), and
\((1-w^{1/5})/q\le.22\). Consequently
\[
 \big|\delta_z-\tfrac32\delta_y\big|
 \le \big(\tfrac32(.22)(.282)+\tfrac38(1.18)(.282)^2(.107)^2
                      +(1.08)(.03)(.107)\big)q^4 < .098q^4.
\]

The absolute opposite-parity coefficient bounds after scaling are
\[
\begin{array}{c|ccc}
 j=5,\text{degree}&4&2&0\\
 &.435\lambda^3&.436\lambda&.058\lambda\\
 j=6,\text{degree}&5&3&1\\
 &.57r^2\lambda^3&.57r^2\lambda&.179r^3\lambda .
\end{array}
\]
The first two degrees for each order follow directly from (2),(3)
using \(|c|,|d|\le .347p^2,.175p^2\) and
\(|c|,|d|\le .378q^2,.304q^2\), respectively.
For the linked degrees,
\[
 |c-2d|\le p^3(.175+.072p+.053p^3),\qquad
 |4c-5d|\le q^3(.912+.306q+.36q^3),
\]
by the skew link, the cubic sinh remainder and the bounds
\(u-v\le p/2,\ w-y\le3q/5\).
The extension checker uses these substitutions to bound each scaled
even coefficient from these Taylor substitutions. The leading term
gives \(|X|<2\), so all arguments still lie in the left half-plane.

## Moving rays and the flat limit

The four scaled ray parameters are
\(b_1=527/500+(7/20)\lambda^2\), \(b_2=43/50\),
\(b_3=53/100\), and \(b_4=39/100\), in increasing order of physical
argument \(\pi-2\arctan(\lambda b_i)\). For \(X=x+iY\), a root on one
of these rays must satisfy
\[
 (Y-b_i)(1+\lambda^2 b_iY)=-b_i\lambda^2x^2,\qquad 0\le Y\le b_i .
\]
This also specifies the limiting ray \(Y=b_i(0)\) at \(\lambda=0\).
The checker bounds \(Y\) by applying this equation twice. It suffices
to check \(0\le x\le2\), since reflection \(X\mapsto-\overline X\)
only reverses the signs of the coefficients of parity opposite to
the order and those coefficients have symmetric enclosures.

Each rational rectangle in \((\lambda,r,x)\) is accepted when some
projection \(\Im(P_j(X)\overline X^{\,h})\), \(0\le h\le j\), is
bounded strictly away from zero. Otherwise the rectangle is bisected.
The exact replay covered the four rays in \(687,81,361,147\)
rectangles at order five and \(71557,6513,44709,3613\) at order six.
It also verifies the leading-term dominance on \(|X|=2\), keeping
all physical roots in the open left half-plane, and the strict
alternation of the limiting Hermite roots about the four ray positions
for every permitted \(r\).

For an actual profile, multiply every Gauss reciprocal at both
orders by the same parameter and decrease it from one to zero.
Then \(u\) increases, the factorwise order ratios are unchanged,
and the bounds for \(r\) still hold. The scaled polynomials converge to
\[
 X^5+\tfrac54 X^3+\tfrac{15}{64}X,\qquad
 X^6+\tfrac32 rX^4+\tfrac9{20}r^2X^2+\tfrac3{200}r^3,
\]
where the limiting \(r\) is the ratio of the first derivatives
of \(1-w\) and \(1-u\) along the path. In the limit the upper roots
in the five successive ray sectors belong to orders six, five, six,
five, and six. The sole real quintic root lies in the last
sector after that sector is joined to its reflection through the
negative real axis. Continuity and the ray exclusions preserve the
root counts in the four outer upper and lower sectors and in that
joined sector. The actual sextic has positive values on the real
axis by the [even-order positivity argument](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md).
In the joined sector the actual quintic's sole root is real by
conjugacy, while the two sextic roots form a nonreal pair. This
proves disjoint arguments. The [characteristic-zero descent](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md)
then applies.
