# Terminal-eleven Cayley separation from outer mean .88 to .94

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the order-ten outer mean after positive endpoint
normalization and removal of the binomial factors. If \(22/25\le
u\le47/50\) and every order-ten Gauss reciprocal is at most \(1/70\),
their tenth and eleventh factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
Rational Bézout gives coprimality over characteristic-zero fields.

Use the [coupled scales, moment bounds and Cayley enclosure
of the middle-band proof](MULTIVARIATE_TERMINAL_ELEVEN_MIDDLE_CAYLEY.md).
The reciprocal caps and moment inequalities apply here with the
outer ratio and deficit bounds
\[
  r=\frac{1-w}{1-u}\in[10094/10000,4043/4000],
      \qquad q=1-w\le 1213/10000.
\]
The checker obtains these bounds from the lower matched
outer-ratio inequality and the upper logarithmic ratio bound
in the middle-band proof. It proves leading coefficient dominance
in the Cayley disk of radius \(12/5\), which maps into the
physical left half-plane since \((12/5)^2q<1\).

The first eleventh-order ray on \(p=1-u\in[.060,.065]\) uses
the narrower interval \(r\in[50487/50000,20209/20000]\) and
\(q\le66/1000\). The lower ratio again follows from the matched
inequality. For the upper ratio, if \(L\) is the enclosure's
upper logarithmic ratio bound, then
\[
 r\le\frac{1-u^L}{1-u}\le\frac{1-(47/50)^L}{1-47/50}.
\]
The first expression increases in \(u\), by convexity of
\(u^L\). The checker encloses the last rational power with
integer power comparisons. It uses a common affine increment
in the square-root deficit to enclose the pure Gaussian
directions on this ray, as in the middle-band proof. The
moment and remainder bounds use \(q_0=66/1000\) for this
subband and \(q_0=1213/10000\) elsewhere.

The checker proves imaginary-axis anchor signs at \(p=.060\)
and verifies every piecewise affine ray by strict centered
interval projection. To use the narrow rectangle along the
continuation, start at the actual \(p,r\), connect to the
Gaussian references at the same \(p,r\), and then move
\(r\) to \(101/100\) at fixed \(p\) in the reference family.
The selected ratio belongs to both admissible intervals.
At fixed \(r=101/100\), move \(p\) down to the anchor.
The ten upper sectors retain the alternating root counts
from the anchor, with its additional negative real odd-order
root. The Cayley height and conjugation argument of the
middle-band proof applies throughout. Positivity of the even
factorial moment on the real axis gives separation in the
joined sector through the negative real axis.

The [exact checker](verify_multivariate_terminal_eleven_upper_cayley.py)
certifies 283232 projection cells in 216 order, ray, and subband
checks. The exact counts are frozen in the checker.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_upper_cayley.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_upper_cayley.py 11 0 60 65
~~~
