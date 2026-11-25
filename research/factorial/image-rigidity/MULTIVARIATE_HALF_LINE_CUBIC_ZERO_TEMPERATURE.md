# Zero separation at zero temperature for the half-line cubic

Let
\[
                F(H)=\int_0^\infty e^{-y^3/3+Hy}\,dy
                      =\pi\operatorname{Hi}(H).
\]

**Theorem.** All zeros of \(F\) are simple, and any two
distinct zeros have distinct imaginary heights. The
order derivative of the
[half-line cubic family](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md)
is nonzero at every zero at \(T=0\).

The proof uses exact exterior estimates and an enclosed
winding computation on three small rectangles. The
published [Scorer zero study](https://ir.cwi.nl/pub/10698/10698D.pdf)
gives asymptotics and numerical tables, which are not
used as the zero counts below.

## Exterior derivative and exclusion

The decreasing-density sine-transform argument in the
[half-line cubic proof](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#all-zeros-at-large-positive-temperature)
puts all upper zeros in the open right quadrant. Let
\(\omega=e^{2\pi i/3}\) and suppose \(H\) lies in that
quadrant, with \(R=|H|\ge8\). The exact decomposition is
\[
 F(H)=\omega F(\omega H)+2\pi i\omega^2\operatorname{Ai}(\omega^2 H).
\]
For \(q=\omega H\), one has \(\Re q\le-R/2\). Comparing
the integrand for \(F(q)\) with \(e^{qy}\), and using
\(1-e^{-y^3/3}\le y^3/3\), gives
\[
  |F(q)+1/q|\le 32 R^{-4},\qquad |F'(q)|\le4 R^{-2}.
\]
Consequently \(F(q)\ne0\) and
\(|\omega F'(q)/F(q)|\le5/R\).

The argument of \(u=\omega^2 H\) lies in
\([-2\pi/3,-\pi/6]\). The leading Ai and Ai-prime terms
have relative errors at most
\((25/24)R^{-3/2}\) and \((49/48)R^{-3/2}\), respectively.
These are the explicit primary-sector bounds of the
[Airy envelope estimate](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4).
In particular Ai is nonzero there and
\[
 \left|\frac{\omega^2\operatorname{Ai}'(u)}
                    {\operatorname{Ai}(u)}-\sqrt H\right|
       \le 3/R,
]
since \(-\omega^2\sqrt u=\sqrt H\) on this quadrant.
It follows that the logarithmic derivative of the
saddle-to-endpoint ratio has real part at least
\(\sqrt{R/2}-8/R>0\).
At a zero the ratio of these two terms is \(-1\),
so this derivative proves simplicity. Two exterior
zeros at the same imaginary height would lie on a
horizontal segment that stays exterior, on which
the ratio modulus strictly increases. Hence they
cannot share a height.

There are no zeros with \(\Re H\ge16\) and
\(0\le\Im H\le8\). Indeed \(\arg H<\pi/6\) there, so
the exponential of the saddle term has modulus
at least one and its leading modulus is at least
\(R^{-1/4}/2>2/R\). The endpoint term is at most
\(2/R\). This also rules out a same-height partner
to a root below height eight outside the rectangle
counted next.

## Bounded winding certificate

The verifier
[verify_half_line_cubic_zero_temperature.py](verify_half_line_cubic_zero_temperature.py)
gives the following counterclockwise winding counts.

| Rectangle in \((\Re H,\Im H)\) | Winding | Certified boundary segments |
| --- | ---: | ---: |
| \([0,16]\times[0,8]\) | 2 | 930 |
| \([1.31,1.33]\times[3.92,3.94]\) | 1 | 37 |
| \([2.71,2.73]\times[6.25,6.27]\) | 1 | 42 |

For reproducibility, run from the repository root with
the python-flint package available.

```text
.venv/bin/python research/factorial/image-rigidity/verify_half_line_cubic_zero_temperature.py
```

Here is the enclosure argument implemented by the verifier.
The Taylor coefficients \(a_k\) follow from
\(F''=HF+1\), with
\(a_0=3^{-2/3}\Gamma(1/3)\), \(a_1=3^{-1/3}\Gamma(2/3)\),
and \(a_2=1/2\). Thereafter
\(a_{k+3}=a_k/[(k+2)(k+3)]\). At \(|H|\le18\),
terms of the value and derivative tails after term
index \(319\) decrease with ratio less than \(1/2\) on
each of the three residue classes. The verifier
therefore uses twice the sum of the first three
omitted terms to bound each tail, with all numbers
enclosed by Arb arithmetic.

For a rectilinear segment of length \(d\) starting at
\(z_0\), the differential system \(F''=HF+1\) bounds
the variation of \(F\) along the segment by
\[
           d\,(|F(z_0)|+|F'(z_0)|+1)e^{18d}.
\]
This follows by integrating the 1-norm growth of
\((F,F',1)\) on the disk of radius \(18\).
The verifier subdivides until this bound plus the
two endpoint midpoint errors is strictly less
than a lower bound for \(|F(z_0)|\).
Thus the true image of each subsegment and its
straight rational midpoint edge lie in the same
disk avoiding zero. A piecewise straight homotopy
reduces the winding to the rational polygon.
Ray crossings of that polygon are counted by exact
rational arithmetic. These checks certify the
boundary avoidance and give the argument principle
counts with multiplicity in the table.

The two inner rectangles exhaust the outer count.
Each contains one simple root, and their imaginary
intervals are disjoint. Any other upper root has
height greater than eight by the exterior exclusion,
so the derivative estimate handles the rest.
Conjugation gives the lower zeros. At \(T=0\)
the order derivative reduces to \(2H\), which is
nonzero at every zero by \(F(0)>0\). \(\square\)

By local finiteness in the
[temperature separation theorem](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md),
there is therefore an open interval about zero
where all three properties continue to hold.
This yields the conditional-model exclusion in the
[bounded moving window](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md#bounded-moving-fields)
whenever its limiting temperature lies in this
interval and both scaled fields are bounded.
