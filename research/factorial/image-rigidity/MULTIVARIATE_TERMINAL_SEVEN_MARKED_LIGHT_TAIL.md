# Terminal-seven separation with one or two marked (1,2) rows

Let all unequal rows of two distinct ordinary exponent columns
be shared. Suppose exactly one or two are \((1,2)\) or \((2,1)\),
and on every other unequal row \(d/(6\ell+1)\le1/8\),
where \(\ell=\min(a_i,b_i)\) and \(d=|b_i-a_i|\).
If
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}<\frac{17}{4},
\]
the order-six and order-seven factorial moment polynomials have
disjoint root arguments up to independent positive rescalings.
Rational Bézout gives coprimality over every characteristic-zero field.

The [exact replay](verify_multivariate_terminal_seven_marked_light_tail.py)
certifies the three orientation patterns, up to exchanging the
columns. The patterns are one increasing row, two increasing
rows, and one row of each orientation.

## Residual outer-mean bounds

Write \(u_R\) for the order-six outer geometric mean of the
unmarked rows and \(S_R=S-m/2\), where \(m\) is the number of
marked rows. For a residual row put \(r=d/\ell\). The reciprocal
cap gives \(r\le6/7\) since \(d/(6\ell+d)\le1/8\).
Its contribution to \(-\log u_R\), divided by its curvature,
is at most
\[
 F(r)=\frac{1+r}{2r}\log(1+5r/6)
                \le \frac{13}{12}\log(12/7).
\]
The first inequality bounds each of the \(d\) outer-mean factor
logarithms using its reciprocal at most \(d/(6\ell)\).
For the second, differentiating gives
\[
 2r^2F'(r)=\frac{(5/6)r(1+r)}{1+(5/6)r}-\log(1+(5/6)r)\ge0
\]
by \(\log(1+x)\le x\) and \(1+r\ge1+(5/6)r\). Hence the exact
bounds
\[
 u_R>(7/12)^{65/16}>.110 \quad (m=1),\qquad
 u_R>(7/12)^{169/48}>.149 \quad (m=2)
\]
follow from the respective curvature assumptions. The replay
verifies the rational-power endpoint comparisons.

Use the following stronger bounds on the residual reciprocals
on the indicated ranges.

| Range of \(u_R\) | Residual reciprocal cap |
| --- | ---: |
| Up to \(.620\) | \(1/8\) |
| \([.620,.662]\) | \(2/19\) |
| \([.662,.900]\) | \(2/25\) |
| \([.900,1]\) | \(1/20\) |

For \(d=1\), a residual row has \(\ell\ge2\), hence its reciprocal
is at most \(1/13\). For \(d=2\), use the exact row mean squared
\((6\ell+1)(6\ell+2)/[(6\ell+11)(6\ell+12)]\).
It increases in \(\ell\); its values at \(\ell=2,3\) are too small
for the second and third bands, respectively. Thus \(\ell\ge3,4\)
there. For \(d\ge3\) the row outer-mean bound
\[
 u_R\le\left(\frac{6\ell+d}{6\ell+6d}\right)^{d/2}
\]
gives \(x\le(1-u_R^{2/3})/(6u_R^{2/3}-1)\). The replay verifies
this bound against both smaller caps at their lower endpoints.
For the final band, use
\(x\le(u_R^{-2}-1)/5\), the bound from the residual outer product.

## Projections and pure marked anchors

The order-\(j\) multiplier at index \(k\) of a single increasing
marked row is
\[
 M_{j,k}=\frac{(j+k)!}{(j!)^{\,1-k/j}((2j)!)^{\,k/j}} .
\]
For a decreasing marked row use \(M_{j,j-k}\).
The checker encloses these \(j\)-th roots of rational numbers
by integer comparisons. At each residual outer-mean interval,
multiply the Gaussian reference coefficients from the
[paired terminal-seven proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md)
by the marked-row multipliers for the given orientation.
Use exactly the same shared residual log contributions and
exponential remainder estimate as in that proof. The marked
multipliers are fixed, so they multiply both the nominal
coefficients and every shared linear projection.

The linked Bernstein projections prove strict ray exclusions
throughout each residual band, with rational affine ray cosines
and outward dyadic arithmetic. All normalized interior
multipliers are at most one on the actual and reference
polynomials and on the coefficient segments between them.
Thus \((11/10)^7<2\) gives the annulus \(1/10<|t|<10\)
used by the checker.

At any residual \(u_R,\rho\), deform to the marked array times
the symmetric residual with those same means. The segment
lies in the verified linked projection enclosures. Vary
\(\rho\) to \(36/35\), then move \(u_R\) to one along the
symmetric residual. At the endpoint the residual multiplier
array is identically one, so the anchors are the pure marked
arrays. The checker verifies their sector counts directly by
disjoint rational disks of radius \(1/500\). At each center
\(c\) it bounds the order-\(j\) polynomial by
\[
 |P'_j(c)|r>|P_j(c)|+\sum_{k=2}^j|P_j^{(k)}(c)|r^k/k!,
 \qquad r=1/500 .
\]
Rouché's theorem gives one root in each disk. The signed
distance from the disk center to each ray puts the whole
disk in the prescribed sector. The successive upper roots
have orders \(7,6,7,6,7,6\); the remaining septic root lies
on the negative real line.

The exact ray exclusions preserve those counts at the actual
profile. After joining the last sector with its conjugate,
its unique septic root is real by conjugacy. Its two sextic
roots are nonreal by positivity of the actual sixth-power
factorial moment on the real line. This proves disjoint
arguments.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_marked_light_tail.py
~~~
