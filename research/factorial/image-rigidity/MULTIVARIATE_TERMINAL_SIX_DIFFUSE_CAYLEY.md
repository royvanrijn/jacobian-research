# Shared terminal six with outer mean between .60 and .72 and reciprocal cap 1/8

Let \(a,b\) be distinct ordinary exponent columns whose unequal
rows are all shared. Let \(u\) be their order-five normalized
outer geometric mean, and put \(x_*=\max_i d_i/(5\ell_i+1)\),
where \(d_i=|b_i-a_i|\) and \(\ell_i=\min(a_i,b_i)\)
on the unequal rows. If
\[
                       \frac35\le u\le\frac{18}{25},
 \qquad                 x_*\le\frac18,
\]
the order-five and order-six factorial polynomials have disjoint
root arguments, hence are coprime over every characteristic-zero
field.

Apply the [intermediate Cayley proof](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md)
with \(p=1-u\in[.28,.40]\), \(\lambda=\sqrt p\),
and \(t=1/8\) as the order-five Gauss reciprocal cap. Use the
same order-ratio interval (1), central defect and linked skew
bounds, and polynomial expansions. The
[replay](verify_multivariate_terminal_six_diffuse_cayley.py)
changes the following estimates for the larger deficit.

- The order-six mean satisfies \(w>57/100\) in the checked boxes.
  Use \(w^{-3/5}<71/50\) and \((3/8)w^{-7/5}<17/20\)
  for the lower and upper errors in the central difference.
  The replay checks both rational power comparisons.
- The larger skew arguments have absolute value at most \(3/50\).
  On this interval the series give
  \(\cosh s\le1.002\), \(|\sinh s|\le1.001|s|\),
  \(\cosh s-1\le .501s^2\), and
  \(|\sinh s-s|\le|s|^3/5\).
  The replay verifies the tail bounds rationally. Thus the
  linear linked-error factor \(1.001\) in each bound for
  \(\epsilon_5,\epsilon_6\) becomes \(1.002\).
  The cubic factors \(1.01,1\) still apply.

The actual outer multipliers are at most \(.72/(1-.06)<19/20\)
by these skew bounds and \(w\le u\); log convexity gives the
same cap on all interior multipliers. The nominal arrays and
their segments to the actual ones obey the same bound. Thus
the earlier radial dominance inequality still places every
root in \(1/8<|t_{\rm phys}|<8\).

The four scaled ray parameters interpolate at the rational
knots in the replay. For a scaled parameter \(b\), the original
cosine is \(c=(p b^2-1)/(1+p b^2)\). The replay computes its
negative upper bound \(C\) on each initial knot interval.
On the verified annulus,
\[
 |\lambda\Re X|\le\frac{63}{65-16C}
\]
by the Cayley real-coordinate formula and radial monotonicity,
using the maximal value at radius \(1/8\) for \(R\le1\)
and reciprocity for larger \(R\).
The same exact projection test certifies each ray over the
enlarged real-coordinate interval. The seven knot intervals
are checked by outward dyadic arithmetic, and the only
acceptance criterion is a strictly signed projection.

For the root counts, deform each actual coefficient vector
to the nominal symmetric vector at fixed \(p,\rho\)
as in the intermediate proof. Next vary \(\rho\) to \(25/24\)
and decrease \(p\) to .28 along the nominal family. Its
ray parameters at .28 equal those of the intermediate proof
exactly. That proof connects this same nominal vector to its
verified Gaussian anchor, so its sector counts apply here.
The quintic has one real root in the final joined sector,
whereas the sextic pair is nonreal by positivity of the
actual sextic factorial moment on the real line. This proves
disjoint arguments.

Run from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_six_diffuse_cayley.py
~~~
