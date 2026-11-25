# A terminal-eleven Cayley bridge from outer mean .94 to .95

For fully shared distinct ordinary exponent columns, let \(u\) denote
the normalized order-ten reflected outer mean. On
\(47/50\le u\le19/20\), the tenth and eleventh factorial moments
have disjoint root arguments under independent positive rescalings.
There is no reciprocal cap.

Put \(p=1-u\), \(q=1-w\), and \(r=q/p\). The
[linked deficit bounds](MULTIVARIATE_TERMINAL_ELEVEN_STABLE_UPPER.md)
with \(u\ge47/50\) give
\[
 100919/100000\le r\le101145/100000,\qquad
 \sqrt p\le2449490/10^7,\qquad q\le61/1000.
\]
The exact checker verifies the rational comparisons. They hold along
the same common reciprocal contraction to the profile with \(p=.05\).

Use the weight-seven Cayley models at orders ten and eleven.
For either deficit \(h\) and \(q_0=61/1000\), the elementary series
bound
\[
 -\log(1-h)\le h+\frac{h^2}{2(1-h)}
      \le h\left(1+\frac{q_0}{2(1-q_0)}\right)
\]
gives the common constant used in the log moment and raw remainder
estimates. Its proof is \(1/k\le1/2\) for all terms of order
\(k\ge2\). All the grouped directions and Taylor cancellations
are the same exact identities as for the stable upper range.
Leading-coefficient dominance gives the radius-two disk and
\(4q_0<1\) keeps its physical image in the left half-plane.

Nine affine ray levels are certified on each of the two deficit
subbands \([.05,.055]\), \([.055,.06]\). At \(p=.05\) they agree
exactly with the certified last levels of the stable upper proof.
Thus the sector counts there apply to the actual reciprocal
contraction, and the ray and root-disk exclusions preserve them
up to the desired profile. The order-eleven root in the sector
joined through the negative ray is real by conjugacy. The
order-ten pair is nonreal by strict positivity of the real even
factorial moment. This also separates the root arguments in that
sector.

The [exact checker](verify_multivariate_terminal_eleven_upper_seam.py)
replays 170,835 projection cells across 36 order, ray and subband
checks and checks the exact match at the common endpoint.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_upper_seam.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_upper_seam.py 11 8
~~~
