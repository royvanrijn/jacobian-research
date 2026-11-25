# Terminal-six separation for two private linear rows

Let \(a,b\) be ordinary exponent columns with a private
linear variable on each side and no other private rows.
Delete rows with equal entries. Write
\[
 S_0=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\max(a_i,b_i)}
\]
for the shared curvature. Each of the two private rows
has curvature one.

**Theorem.** The ordinary factorial moments of orders
five and six are coprime for every such profile, in
arbitrary dimension and over every characteristic-zero
field. Their independently normalized polynomials have
disjoint root arguments.

For \(S_0\ge19/5\), the weighted curvature is
\(S+P=S_0+3\ge34/5\), so the
[weighted-curvature theorem](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md)
applies. We prove the complementary range by retaining
the same value of \(S_0\) in both coefficient arrays.
The two private factorials \(k!(j-k)!\) cancel the binomial
factor. Under the positive endpoint normalization the
remaining coefficient is
\[
 c_{j,k}=
 \prod_{\rm shared}
 \frac{((j-k)a_i+kb_i)!}
      {((j a_i)!)^{1-k/j}((j b_i)!)^{k/j}},\qquad
 c_{j,0}=c_{j,j}=1 .
\]

## A two-sided loss estimate in both orders

For one increasing shared row write \(a_i=l\),
\(b_i=l+d\), \(r=l+d\), \(s=d^2/r\), where \(l,d\ge1\).
The logarithms of its raw coefficients have second
difference at index \(k\) equal to
\[
 \sum_{h=1}^d
 \log\left(1+\frac{d}{j l+(k-1)d+h}\right).
\]
For \(1\le k<j\), this quantity lies between \(s/j\)
and \(d^2/(jl)\). The lower bound uses
\(\log(1+x)\ge x/(1+x)\) and
\(jl+kd+h\le jr\). The upper bound uses
\(\log(1+x)\le x\). Reflection gives the same bounds
for decreasing rows. Subtracting the endpoint chord
and comparing with the endpoint-zero quadratic gives
a lower bound \(s k(j-k)/(2j)\) on the normalized log
loss at index \(k\).

The opposite bound \(13s k(j-k)/(12j)\) holds at orders
five and six for every row with \(s\le19/5\).
If \(l\ge51\) and \(d/l\ge1/3\), then
\[
 s=l\,\frac{(d/l)^2}{1+d/l}\ge l/12>19/5 ,
\]
so \(r/l<4/3\) and the upper second-difference bound
suffices. For \(l\le50\) the condition \(s\le19/5\)
forces \(d\le15\), since \(d^2/(50+d)\) increases
and its value at \(d=16\) exceeds \(19/5\).
The [exact checker](verify_multivariate_terminal_six_two_private_linear.py)
tests the remaining rows and indices by the inequalities
\[
 \frac{((jl)!)^{j-k}((j(l+d))!)^k}{((jl+kd)!)^j}
 \ \le\ \sum_{h=0}^{64}\frac{u^h}{h!},
 \qquad
 u=\frac{13}{6}\,s\,\frac{k(j-k)}2 .
 \tag{1}
\]
The right side is a lower bound on \(e^u\), and the
left side is the reciprocal of the \(j\)-th power
of the normalized coefficient. This is a finite
integer and rational calculation on the specified
row range. Summing the two log-loss bounds yields
\[
 e^{-(13/6)S_0 k(j-k)/(2j)}
       \le c_{j,k}\le e^{-S_0 k(j-k)/(2j)} \qquad(j=5,6).
 \tag{2}
\]

## A linked curvature band and its symmetric anchor

Cover \([0,19/5]\) by the 152 intervals
\([b/40,(b+1)/40]\). In a fixed band, use (2) with the
right endpoint for the lower coefficient bound and the
left endpoint for the upper. The anchor in that box
has coefficients \(c_{j,k}^*=\exp(-s_* k(j-k)/(2j))\),
where \(s_*\) is the band midpoint. Both orders use the
same midpoint and the same four decreasing ray cosines
\[
 \frac{14}{25}+\frac38\,\frac{s_*}{s_*+3/2},\qquad
 \frac15,\qquad-\frac{27}{100},\qquad-\frac{13}{20}.          \tag{3}
\]
The coefficients are at most one. The constant term
dominates on \(|z|=1/5\) and the leading term dominates
on \(|z|=5\), by the geometric sums.

The checker applies the
[exact rational ray projections](MULTIVARIATE_TERMINAL_SIX_SMALL_SHARED.md#exact-angle-certificates)
to every band and ray. It encloses the negative
exponentials with reciprocal Taylor bounds, using
\[
 \sum_{h=0}^{64}\frac{u^h}{h!}
 \ \le\ e^u\ \le\
 \sum_{h=0}^{64}\frac{u^h}{h!}
       +\frac{u^{65}}{65!(1-u/66)}
 \qquad (0\le u<66).
\]
All arguments in the boxes meet the displayed range.
Thus the ray exclusions hold for the entire box,
including every linear interpolation from its
anchor to an actual pair.

The anchor root counts follow from real sign tests.
For \(u=z+z^{-1}\), the quintic and sextic anchors
factor respectively as
\[
 \begin{aligned}
 (1+z)z^2\big[u^2+(c_{5,1}^*-1)u+c_{5,2}^*-c_{5,1}^*-1\big],\\
 z^3\big[u^3+c_{6,1}^*u^2+(c_{6,2}^*-3)u
                                  +c_{6,3}^*-2c_{6,1}^*\big].
 \end{aligned}
\]
At twice the four cosines (3), the first bracket has
signs \(+,-,-,+\) and the second has signs
\(-,-,+,+\), checked with the same rational exponential
bounds. The sextic bracket is positive at \(u=2\) and
negative at \(u=-2\). For the latter sign use
\(2-2c_{6,1}^*+2c_{6,2}^*-c_{6,3}^*>0\), which follows
from \(1\ge c_{6,1}^*\ge c_{6,2}^*\ge c_{6,3}^*>0\).
Each indicated interval contains exactly one root of
the corresponding real bracket by its degree.
Since these roots are in \((-2,2)\), they yield the
upper and lower unit-circle roots of the anchor.

It follows that the five successive upper sectors
have, in order, one sextic, one quintic, one sextic,
one quintic and one sextic anchor root. The negative
quintic root is \(-1\). Under the interpolation, the
sectors about the real axes are counted with their
conjugates. In particular the last sector centered
on the negative real axis has only one quintic root,
so it stays real. The actual sextic is positive on
the real axis by its factorial integral, so it cannot
share that root's argument. The other sectors have
disjoint root arguments by their counts. This
proves the theorem, and rational Bézout transfers
coprimality to every characteristic-zero field.

The exact replay checks 2192 row-index inequalities
for order five and 2740 for order six. The ray checks
certify 6552 and 8088 radial intervals, respectively.
