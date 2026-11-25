# Terminal-six separation with a shared (1,2) row for outer mean at least 13/20

Suppose all unequal rows of two distinct ordinary exponent columns
are shared and one of them is \((1,2)\) or \((2,1)\). Let \(u\) be
the order-five normalized outer geometric mean. If \(u\ge13/20\),
the order-five and order-six factorial polynomials have disjoint
root arguments under independent positive rescalings. This gives
coprimality over every characteristic-zero field.

Only \(13/20\le u\le18/25\) needs a new check, since the preceding
[intermediate](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md)
and [near-flat](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT_EXTENSION.md)
proofs cover \(u\ge18/25\). Orient the marked row from smaller to
larger exponent by exchanging the columns if necessary. It is
unique, because two such rows would give \(u\le3/5\).

## Integer bounds on the residual rows

The marked row contributes \(\sqrt{3/5}\) to \(u\). Write
\(u_R=u/\sqrt{3/5}\) and \(p_R=1-u_R\). If \(u\ge13/20\),
then \(u_R>.839\) by squaring, so \(p_R<.161\).
The marked part of the existing certificate covers \(p_R\le.071\).
For the extension, use the following caps on the order-five
reciprocals of the residual rows.

| Range of \(p_R\) | Reciprocal cap |
| --- | ---: |
| \([.071,.083]\) | \(1/26\) |
| \([.083,.105]\) | \(1/21\) |
| \([.105,.143]\) | \(1/16\) |
| \([.143,.161]\) | \(1/11\) |

For a unit-gap row with lower exponent \(\ell\), the row mean is
\(\sqrt{(5\ell+1)/(5\ell+5)}\). Using \(u_R\) as a lower bound
on every row mean forces \(\ell\ge5,4,3,2\), respectively.
For a gap \(d\ge2\), use the row bound
\[
 u_R\le\left(\frac{5\ell+d}{5\ell+5d}\right)^{d/2}
       \le\frac{5\ell+d}{5\ell+5d}
\]
to obtain \(x_5=d/(5\ell+h)\le(1-u_R)/(5u_R-1)\). Its value
at the lower endpoint of \(u_R\) is below each listed cap.
The checker verifies every endpoint inequality exactly.

## Ray continuation

Use the marked-row linear contributions and the outer-log ratio
interval (1) from the [intermediate proof](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md).
The [extension checker](verify_multivariate_terminal_six_marked_row_extension.py)
calls the same Cayley coefficient expansion, interval logarithm,
ray equation, and rotated projections. It uses no numerical
root computations. Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_six_marked_row_extension.py
~~~

The four ray cosines interpolate linearly between the displayed
rational knots. At \(p_R=.071\) they coincide with the earlier
certificate's rays. Each cosine stays at most \(-1/4\).
The marked row's interior multipliers are at most \(4/5\);
all residual multipliers on the actual and nominal arrays are
at most one. The earlier dominance bound thus keeps all roots in
\(1/8<|t|<8\), including along their straight coefficient segments.
With the unscaled Cayley transform \(t=(X-1)/(X+1)\), this gives
\[
 |\Re X|\le
 \frac{1-(1/8)^2}{1+(1/2)(1/8)+(1/8)^2}
 =\frac{21}{23}
\]
on every checked ray. For radii above one, the same bound follows
by replacing the radius by its reciprocal. The checker uses this
larger real interval for each ray and both orders. The nine
deficit intervals are verified by outward dyadic interval
arithmetic and strict signed projections.

At fixed \(p_R,\rho\), deform the actual arrays to the marked
array times the nominal symmetric residual. The certified joint
linear intervals contain this entire straight segment. Then vary
\(\rho\) to \(25/24\), which belongs to every outer-log ratio
interval in the replay, and decrease \(p_R\) to .071 along that
nominal family. At .071 the earlier marked certificate continues
to its pure marked-row anchor. The counts in the sectors
and their conjugates are consequently the same as in that
certificate. The last sector joined to its conjugate contains
one real quintic root and two sextic roots. The actual sextic
is strictly positive on the real line by the factorial moment
representation. The sextic roots in that sector form a nonreal
pair, proving the required separation.
