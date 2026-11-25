# Terminal-ten separation for one private linear row and for retained shared rows

Let \(a,b\) be distinct ordinary exponent columns. Delete equal rows.
On each shared unequal row, let \(\ell\) be the smaller exponent,
\(d\) its gap, and \(s=d^2/(\ell+d)\). Let \(S_0=\sum s\) over these
rows. The order-nine reciprocal of a Gauss factor is
\(x_h=d/(9\ell+h)\) with \(1\le h\le d\).

**Theorem.** The order-nine and order-ten moment polynomials have
disjoint root arguments under independent positive rescalings if
there is just one private row and its nonzero exponent is one, or
if all unequal rows are shared and some \(d/(9\ell+1)>1/8\).
Their rational Bézout identity gives coprimality over every
characteristic-zero field.

Together with the [larger private-degree theorem](MULTIVARIATE_TERMINAL_TEN_PRIVATE_DEGREES.md),
this covers every terminal-ten profile with a private row.
For fully shared profiles it also gives an unconditional separation
bound at \(S_0\ge3/2\) in combination with the light proof. To see this,
a Gauss factor gives
\(\log(1+8x_h)\ge 8d/(9\ell+h+8d)\ge 8d/(9(\ell+d))\).
Hence the total outer mean satisfies
\(u\le e^{-4S_0/9}\le e^{-2/3}<13/25\).
The strict numerical comparison follows from the exponential Taylor
sum through degree three. If every reciprocal is at most \(1/8\),
the [terminal-ten light theorem](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md)
covers this outer range. Otherwise the theorem above applies.

## Exhaustive retained lists

The [fixed-ray certificate](MULTIVARIATE_TERMINAL_TEN_PRIVATE_DEGREES.md)
covers the single private row when \(S_0\ge57/10\), and all shared
profiles when \(S_0\ge77/10\). Use these cutoffs as budgets for the
remaining two families. In each family retain all shared rows
with \(d/(9\ell+1)>c\), where \(c=1/8\).
This implies \(\ell<d/(9c)-1/9\), hence
\(d<S_0(1+1/(9c))\). The [exact checker](verify_multivariate_terminal_ten_retained_profiles.py)
enumerates those integer bases and gaps and all their multiplicities
subject to the budget. In the fully shared family the list is nonempty,
and simultaneous reversal of all orientations leaves the question
invariant. The private row fixes the increasing orientation on one
side in the other family, so both shared orientations are kept there,
and the list may be empty. There are 988 fully shared and 280
single-private lists.

For each retained list the core in both orders is the exact product
of its normalized factorial factors, with the private factor also
included in the private family. Write \(B\) for the remaining
shared-curvature budget and \(v\) for the order-nine reflected outer
mean of the light residual tail. All residual reciprocals at this
order are at most \(c\).

## Tail strength and the two linked polynomials

Put \(U=-\log v\). For a residual row with \(x=d/(9\ell+1)\),
\[
 \frac{U_{\rm row}}s
 \le \frac{d}{2s}\log(1+8x)
 \le \frac12\left(1+\frac1{9x}\right)\log(1+8x)
 \le \frac{17}{18}\log 2 .
\]
The last function increases in \(x\). To verify this, put \(z=8x\).
The derivative has the sign of
\(z(1+(9/8)z)-(1+z)\log(1+z)\), whose derivative is
\((9/4)z-\log(1+z)>0\) for \(z>0\).
Consequently \(2^{-17B/18}\le v\le1\). A subinterval with
left endpoint \(v_-\) also enforces \(x_h\le(v_-^{-2}-1)/8\)
for all residual factors, by the product formula for \(v^{-2}\).

For the normalized residual multipliers \(F_{j,k}\), removing
also the binomial coefficient, the
[adjacent endpoint bounds](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
give
\[
 F_{10,k}\le F_{9,k}\,v^{k^2/126},\qquad
 F_{10,k}\le F_{9,k-1}\,v^{(10-k)^2/126},
                 \qquad 1\le k\le9.
\]
Indeed fully shared rows have maximal exponents at least two, so their
attenuation in strength is at least \(\log(190/189)\ge1/190\).
The tail curvature is at least \(U/((17/18)\log2)\), and the checker
verifies \(190(17/18)\log 2<126\) by a rational log enclosure.
The checker intersects the independent coefficients with these upper
bounds, multiplying by the exact order-ten core at the corresponding
index. The comparison Gaussians with any admissible outer log ratio
\(\rho\) satisfy these inequalities as well. The ratio's lower bound
is at least \(\log(161/80)/\log2>1+1/126\), by evaluating the
decreasing lower ratio bound at cap \(1/8\). For the left endpoint,
the Gaussian exponent difference at ratio \(81/80\) is \(k^2/80\),
and its decrease at a smaller ratio is
\((81/80-\rho)k(10-k)/9\). The bound on \(\rho\) gives
the needed \(k^2/126\). Reflection gives the other endpoint.

Each checker cell uses the same core at both orders and one residual
outer interval. It applies the
[exact terminal-ten light coefficient model](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md),
retaining linked mean and skew parameters and using the quadratic
version of the grouped model. The signed tail parameter may run over
\([-1,1]\) because the core fixes the orientation. The core factors
multiply each grouped coefficient at its own degree with outward
bounds. For the independent coefficient boxes the checker uses the
outer-loss intervals and attenuation just proved. A ray can use
either exact enclosure at either order. The signed parameter
is subdivided when a grouped replay first fails.

## Common homotopies and sector verification

At the actual residual mean, scale the signed log skew and mean
defects by \(t\) from one to zero, leaving the core and outer log
ratio fixed. The tanh moment and tail bounds remain valid, including
the Jensen lower bounds since \(t^{2l}\le t\) for \(0\le t\le1\).
The logarithmic attenuation inequalities and the independent
logarithmic coefficient intervals are convex, so they also remain
valid in this step. At zero defects, move the ratio to \(81/80\),
then move the residual mean to the cell midpoint. The same bounds
hold throughout, including the attenuation of the Gaussians.
The annulus \(1/8<|z|<8\) follows from domination using the normalized
private and retained core times the binomial coefficients. The
residual normalized multipliers along this path are at most one.

At each midpoint, the checker isolates the comparison roots by
outward complex balls. It chooses rational ray levels strictly between
the certified upper root arguments, and requires the two orders to
alternate throughout the upper half-plane at this anchor. It applies
strict interval projection signs to exclude the selected rays along
the entire cell. The negative real nonic root remains real in the
sector joined across the negative real axis. The actual tenth moment
has positive value at every real argument, by its even Gamma integral,
which also separates its two roots from the nonic root in that sector.

The exact replay covers 395 outer cells and 24,871 radial cells for
one private row, and 1,334 outer cells and 76,274 radial cells for
fully shared retained rows. A shard is a deterministic slice of the
full finite enumeration. The optional native path uses the same
checked Bernstein arithmetic as the light-ray certificate.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_retained_profiles.py --native private
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_retained_profiles.py --native shared
~~~
