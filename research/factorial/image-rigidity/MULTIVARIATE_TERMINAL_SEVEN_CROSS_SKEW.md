# Linked sextic and septic skew for terminal-seven separation at outer mean .48 to .62

Let the unequal rows of distinct ordinary exponent columns
\(a,b\) be shared. Put \(\ell_i=\min(a_i,b_i)\),
\(d_i=|b_i-a_i|\) on those rows and let \(u\) be the
normalized order-six outer geometric mean. If
\[
 .48\le u\le .62,\qquad \max_i d_i/(6\ell_i+1)\le1/8 ,
\]
the order-six and order-seven factorial moment polynomials
have disjoint root arguments under independent positive
rescalings. Their rational Bézout identity proves coprimality
over every characteristic-zero field.

## Transfer across one matched factor

Write \(\eta\) for the sextic inner skew and \(\gamma\) for
the septic innermost skew, in the notation of the
[paired terminal-seven proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md).
For a single pair of matched factors \(x,x'\), their magnitudes are
\[
\begin{aligned}
 E(x)&=\frac{2\log(1+2x)+2\log(1+3x)
      -\log(1+x)-\log(1+4x)-\log(1+5x)}6,\\
 G(x')&=\frac{6\log(1+3x')-\log(1+x')-\log(1+2x')
      -\log(1+4x')-\log(1+5x')-\log(1+6x')} {14}.
\end{aligned}
\]
Both signs are negative in increasing orientation and positive
in decreasing orientation. This follows from the single-factor
log-skew formulas and convexity of the magnitudes in the local
strength, with value and derivative zero at zero. In particular
\(G\) is increasing in the reciprocal on the range used here.
The
[exact skew-convexity check](verify_multivariate_terminal_seven_private_linear.py)
supplies this fact.

The cap \(x\le1/8\) and the matching relation
\(x'=6x/(7-xh/d)\) give \(6x/7\le x'\le48x/55\).
The [checker](verify_multivariate_terminal_seven_cross_skew.py)
verifies
\[
       \frac{107}{200}E(x)\le G(6x/7),\qquad
       G(48x/55)\le\frac{29}{50}E(x).
\]
For each inequality take the derivative of the difference.
Multiply by the positive product of its linear denominators
and divide out one factor \(x\), as the derivative vanishes
at zero. The resulting polynomial is strictly positive by
its Bernstein coefficients on \([0,1/8]\). Both original
differences vanish at zero, so this proves the inequalities.
Summing with the orientation signs gives the joint bound
\[
             \big|\gamma-\tfrac{14}{25}\eta\big|
           \le\frac1{40}\sum_{\rm factors} E(x).
\]
This is a bound on a signed difference, so cancellations
of opposite orientations do not invalidate it.

## Skew-conditioned Bernstein boxes

Use the [cap-extension model](MULTIVARIATE_TERMINAL_SEVEN_LIGHT_CAP_EXTENSION.md)
with \(t_6=1/8\). Write \(U=-\log u\) and let \(K_6\) be the
upper coefficient satisfying \(\sum E(x)\le K_6 U\)
from its single-factor outer and linked inner-skew bounds.
After exchanging columns if necessary, take \(\eta\ge0\).
The checker subdivides \(u\) and the scalar \(\eta/(K_6 U)\)
in \([0,1]\). In each cell the order-six primary skew uses
that scalar multiplied by \(K_6U\). The order-seven primary
skew uses
\((14/25)\eta+[-K_6U/40,K_6U/40]\), intersected with its
separate skew bound. The other shared log contributions and
the coefficientwise exponential remainder are as in the paired
model. The five ray cosines interpolate affinely in \(u\);
they are independent of the skew subdivision.

All actual and symmetric interior multipliers are at most one.
Deform the actual to the symmetric profile by scaling its
skew and defect contributions in the **logarithms** of the
coefficients from one to zero while keeping \(u,\rho\) fixed.
Their homogeneous linear bounds, including the transfer,
remain valid in every intermediate profile. The multipliers
on this path are geometric interpolations of values at most
one, so the annulus dominance bound \((11/10)^7<2\) applies.
Then vary \(\rho\) to \(36/35\) and move \(u\) to .48 in
the symmetric Gaussian arrays. The ray cosines there match
the preceding cap-extension certificate exactly.

The outward dyadic Bernstein replay proves strict weighted
ray signs in every accepted box. It preserves the preceding
anchor counts. The final sector joined to its conjugate
contains a real septic root and a nonreal sextic pair, by
conjugacy and strict positivity of the actual sextic
factorial moment. All root arguments are disjoint.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_cross_skew.py
~~~
