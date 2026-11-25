# Terminal-thirteen fixed-ray separation for shared loss and private rows

Let \(a,b\) be distinct ordinary exponent columns. Write
\(S_0=\sum_i d_i^2/(\ell_i+d_i)\) over shared unequal rows,
where \(\ell_i\) is the smaller exponent and \(d_i\) is the gap.
Let \(p,q\) be the two private degrees and put \(S=S_0+p+q\).

**Theorem.** The order-twelve and order-thirteen factorial moment
polynomials have disjoint root arguments under independent positive
rescalings if \(S\ge63/8\), \(p+q\ge4\), or one of the private
allocations below satisfies its shared-curvature threshold.
In particular the two moment polynomials are coprime over every
characteristic-zero field in these cases.

| Increasing private gaps | Decreasing private gaps | Required shared curvature \(S_0\) |
| --- | --- | ---: |
| none | none | \(63/8\) |
| \(1\) | none | \(6\) |
| \(1,1,1,1\) | none | \(0\) |
| \(1,1\) | \(1,1\) | \(0\) |
| \(1,1,1\) | \(1\) | \(0\) |
| \(3\) | none | \(3/8\) |
| \(1,1,1\) | none | \(9/8\) |
| \(1,2\) | none | \(3/4\) |
| \(1,1\) | none | \(27/8\) |
| \(2\) | none | \(23/8\) |
| \(1\) | \(1\) | \(13/4\) |
| \(1,1\) | \(1\) | \(7/8\) |
| \(2\) | \(1\) | \(3/8\) |

Column reversal preserves each assertion. Each box also applies
after adding private factors or increasing their gaps.

## Coefficient boxes

Use the [terminal-eleven coefficient boxes](MULTIVARIATE_TERMINAL_ELEVEN_HIGH_CURVATURE.md#coefficient-boxes)
at orders \(j=12,13\). Positive rescaling makes both endpoints
one. A private gap \(d\) in the increasing orientation contributes
\[
                f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}},
\]
and its decreasing version uses \(j-k\). The private factor
decreases with its gap and satisfies \(f_{j,k}(d)\le f_{j,k}(1)^d\).
Hence the zero-curvature cases of the table cover a private degree
at least four on one side, private degrees at least two on both
sides, or a private degree at least three with positive opposite
private degree. These possibilities exhaust \(p+q\ge4\).

The shared normalized factors at index \(k\) have product at most
\(\exp(-S_0 k(j-k)/(2j))\), by their logarithmic second-difference
bounds. This gives interval caps with the listed private gaps.
The same second-difference lower bound with smaller exponent
zero gives \(d/j\) for a private gap \(d\). Thus the total-curvature
cap \(\binom jk\exp(-S k(j-k)/(2j))\) is valid for all profiles and
the first row of the table also certifies \(S\ge63/8\).
If \(E\) is the curvature exponent in the \(j\)-th power of a cap,
the checker lower-bounds \(e^E\) by the fourth power of the finite
positive Taylor sum for \(e^{E/4}\) through degree 64. Integer
power comparisons then round the caps outward. The exact annulus checks bound all normalized
roots between radii \(1/8\) and \(8\), uniformly over each box.

## Exact rays and parity at the negative ray

The [checker](verify_multivariate_terminal_thirteen_high_curvature.py)
uses eleven ordered rational ray cosines for each case.
The reference polynomials \(1+z^{12}\) and \(1+z^{13}\)
belong to their coefficient boxes. The
[exact Bernstein projection](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
excludes a root on each ray in the closed annulus, including during
straight interpolation from either reference to its covered array.

At a sector boundary, the reference sign polynomial in the ray
cosine has value \(\cos(6\theta)\) at order twelve and
\(\cos(13\theta/2)/\cos(\theta/2)\) at order thirteen, using the
limiting value for the latter at \(\theta=\pi\).
Both start positive at \(\theta=0\). At boundary index
\(i=0,\ldots,12\) the checker verifies that the sign is
\((-1)^{\lfloor(i+\mathbf{1}_{j=13})/2\rfloor}\) at order \(j\).
The twelve alternating upper sectors contain the six upper
nonreal roots of each reference, one per sector in successive
orders thirteen and twelve. The order-thirteen reference has
one further root at \(-1\).

Join the first and last sectors to their reflections over the
real axis. All full-sector counts remain unchanged along the
coefficient interpolation, since the rays and annulus boundaries
stay root-free. The last joined sector contains exactly one
order-thirteen root, so conjugation fixes it on the negative ray.
The order-twelve actual moment is strictly positive on the real
axis by the Gamma integral of an even power of a nonzero binomial.
Its two roots in that sector are therefore nonreal. No common
argument is possible in the last joined sector or in any of the
other sectors. The coefficient polynomials are integral, so a
rational Bézout identity transfers their coprimality to every
characteristic-zero field.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_high_curvature.py
~~~

The fixed table checks 13 cases at both orders, 26 annulus tests,
338 anchor signs, and 286 ray exclusions with 1,612 radial cells.
