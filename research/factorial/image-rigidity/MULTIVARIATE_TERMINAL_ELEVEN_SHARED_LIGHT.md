# Terminal-eleven light shared profiles through outer mean .76

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their endpoint-normalized order-ten reflected outer
mean, after removing the binomial factors. If \(u\) lies in a band
below and all reciprocals \(d/(10\ell+1)\) are at most its listed cap,
then the tenth and eleventh factorial moment polynomials have disjoint
root arguments under independent positive rescalings. Thus they are
coprime over every characteristic-zero field.

| Outer mean \(u\) | Reciprocal cap |
| --- | ---: |
| \([0,.42]\) | \(1/8\) |
| \([.42,.44]\) | \(93/1000\) |
| \([.44,.57]\) | \(9/100\) |
| \([.57,.64]\) | \(3/50\) |
| \([.64,.70]\) | \(3/40\) |
| \([.70,.73]\) | \(11/200\) |
| \([.73,.76]\) | \(9/200\) |

## Exact outer-mean and coefficient bounds

For a shared row with smaller exponent \(\ell\) and gap \(d\), its
squared order-ten mean is
\[
 u_{\ell,d}^{,2}
   =\prod_{h=1}^d\frac{10\ell+h}{10\ell+9d+h}.
\]
If its reciprocal exceeds a floor \(t>0\), the exact mean satisfies
\(u_{\ell,d}^2<(1+9t/(1+t))^{-d}\), by
\(10\ell+d<d(1+1/t)\). Also each row mean is at least the product
mean. At a slab starting at \(u_->0\) this gives a finite search over
all rows with reciprocal above \(t=1/40\). The checker intersects
the declared cap with the resulting ordinary cap. For example the
extra knot \(u=.475\) drops the maximal reciprocal \(3/41\) once the
associated squared mean is smaller than \(.475^2\).

The matched upper-order reciprocal is \(x'=10x/(11-ax)\), where
\(a=h/d\in(0,1]\). The upper cap associated to a lower cap \(c\) is
\(10c/(11-c)\). For \(U_j=-\log u_j\) the ratio \(\rho=U_{11}/U_{10}\)
is bounded by
\[
 \frac{\log(1+(100/11)c)}{\log(1+9c)}
      \le\rho\le
 1+\frac{\log(1/(1-c/11))}{\log(1+9c)} .
\]
The inequalities follow from
\(1+(100/11)x\le 1+10x'\le (1+9x)/(1-x/11)\) and the logarithm
quotient comparisons from the
[terminal-ten light proof](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#ordinary-caps-and-exact-coefficient-enclosures).
The [all-order outer-ratio lemma](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
gives \(\rho\ge100/99\) in the slabs that use the tighter lower bound.
The reference ratio \(100/99\) lies in every interval used.

The [matched tanh moment lemma](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
gives the signed second-moment transport with center
\((100/99)^2(1+\eta^2)\) and error coefficient \(2(100/99)^2\eta\),
where \(\eta=c/(22+98c)\). A common column reversal lets the
order-ten signed moment parameter range from zero to one in units
of its maximum magnitude. The skew and mean-defect directions
are the exact rational directions of the
[terminal-ten light model](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison)
evaluated at orders ten and eleven. The early matched bands keep
the quadratic exponential enclosure and the later ones keep the
cubic enclosure. The independent reflected mean and skew bounds
also apply at both orders, since the checker verifies the cap
condition \(T_*^2<1/2\) needed in their monotonicity calculation.

## Ray exclusions and root counts

The checker uses independent boxes on \(0\le u\le.04\). They
include the Gaussian references and the endpoint arrays
\(1+z^{10},1+z^{11}\) at \(u=0\). In the other slabs it checks
the grouped coefficient models along nine rays with piecewise
affine rational cosines. The comparison path at fixed \(u,\rho\)
scales the signed log skews and mean defects jointly to zero.
The Jensen lower bounds remain valid because scaling the
signed moment by \(s\in[0,1]\) multiplies their lower sides by
powers at least two, while the actual unsigned mean-defect
parameters scale by \(s\). The Gaussian comparison then moves
\(\rho\) to \(100/99\) and \(u\) along the slab to its anchor.
The independent boxes used for root moduli contain this path.

Leading and constant coefficient domination gives the root
annulus for the first bands. Where this fails, the checker
instead verifies strict adjacent-ratio inequalities for the
independent boxes throughout the band and uses the
Eneström–Kakeya annulus. The exact Bernstein sign test for the
rays uses common displacements of the outer mean and angle as
in the terminal-ten light proof. The native path uses outward
dyadic interval arithmetic. At order eleven a common positive
scale keeps the Bernstein transform entries within the exact
integer bounds. Some bands also enclose small projected groups
in the independent base before retaining the complete model as
a fallback. Both enclosures keep the same mean-ray coupling.

At each lower band endpoint the checker isolates the two
comparison polynomials' roots in disjoint rational disks or
uses the exact binomial cosine signs. In increasing upper
argument the first ten sectors have roots of orders
\(11,10,11,10,11,10,11,10,11,10\). The order-eleven root in
the joined sector through the negative ray stays real under
the comparison. The order-ten factorial moment is strictly
positive on the real axis by its Gamma integral, and thus its
pair in that sector is nonreal. The ray exclusions and the
annulus therefore give disjoint arguments at the actual arrays.

The [exact checker](verify_multivariate_terminal_eleven_shared_light.py)
replays 276,777 exact radial cells across the fixed rational ray and
root-disk tables. It verifies the expected native cell counts.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_shared_light.py --native
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_shared_light.py --native 0 .47 .70
~~~
