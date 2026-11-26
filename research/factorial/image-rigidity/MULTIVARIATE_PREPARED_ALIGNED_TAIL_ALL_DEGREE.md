# Arbitrarily weak core extensions with all-degree aligned tail stability

Fix a finite ordinary two-monomial core \(C\), possibly empty,
and take an integer \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For every \(2\le n<N=n+1<J\), choose one prime \(p_n\mid n\).
Take a common modulus \(M\) divisible by every \(p_n^{E_n}\)
with \(p_n^{E_n}>N\). The modulus can be chosen solely from
\(J\). Choose also a fixed prime gap \(D>J^2\).

**Theorem.** One can append balanced pairs of ordinary shared
rows of gap \(D\) to \(C\), with arbitrarily small total
curvature, to form a core \(C'\) having the following
uniform property. There is an effective integer \(S_0\)
such that every tail of at least \(S_0\) unit-gap rows
of each orientation, all with smaller exponent
\(a\equiv-1\pmod M\), gives SFC(2) with \(C'\) in every
degree. The tail heights may vary independently and be
unbounded. There is no constraint on the signed tail
gap. The number of rows added to prepare the core can
be bounded in terms of \(J\) alone.

The additional core rows have signed gap zero. Their
normalized real coefficient contributions and signed center
corrections can be made arbitrarily small. Thus the theorem
applies near any finite core without requiring that the
original core have a uniform real Gaussian separation.

## Independent affine variations of the central slopes

It suffices to force the fixed central slope sets from the
[aligned unit-gap theorem](MULTIVARIATE_ALIGNED_UNIT_GAP_CENTRAL_NEWTON.md)
to be disjoint at each designated prime. By the
[absolute signed-gap bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md),
the only tail gaps \(h\) to check satisfy
\[
                     |h+\Delta_C|<146J-73 ,
\]
where \(\Delta_C\) is the original total signed gap.
This is a finite list and is unaffected by balanced
preparation. To deal with negative \(h\), include the same
central comparison for the reversed core and gap \(-h\).

Fix one of the finitely many assigned primes \(p\).
Use the union of the order sets \(n,N\) assigned to it.
For each pair choose the marks of the
[fixed-gap interpolation basis](MULTIVARIATE_UNIT_GAP_LOCAL_COEFFICIENT_INTERPOLATION.md#a-fixed-larger-gap-for-arbitrary-primes).
For clarity, at the order \(j_0\) coprime to \(p\), take a
noninteger \(u=h_0/j_0\) in the interior of each block
\(((k-1)D/j_0,kD/j_0)\), with integer \(h_0\). At the
other order take an integer \(u\) strictly inside each
block. All these \(u\) are \(p\)-integral. They avoid
block boundaries in every order \(j<J\), not just
their assigned pair. For an integer mark this is
because \(D\) is prime and \(0<u<D\). For a noninteger
mark a boundary would give \(j h_0=j_0 k D\), which
would force \(D\mid h_0\), contrary to its location
strictly between consecutive multiples of \(D\).

Use the union of these marks, dropping duplicates.
For a mark \(u\) let \(v_u\) be the vector of indicators
of the blocks containing the integer offset \(ju\),
over all the orders assigned to \(p\). Use one
increasing row with smaller exponent \(p\)-adically
close to \(-u\) and one decreasing row with smaller
exponent close to \(-(D-u)\).
Take the two depths above a fixed bound depending only
on the unmarked offsets, and allow their difference
\(x_u\) to vary over arbitrary integers. All unmarked
factor valuations remain fixed. The paired marked
factors contribute exactly \(x_u v_u\) to the vector
of adjacent coefficient ratio valuations, since each
reflected mark is in the same ratio coordinate and the
two factors \(j\) cancel. Thus the total vector, including
the original core and fixed unmarked contributions, has
the affine form
\[
                         y=y_0+\sum_u x_u v_u .              \tag{1}
\]
The projection of the linear part onto the ratio coordinates
of any assigned pair is surjective as a rational linear
map. Its corresponding marks include the single coordinate
directions at \(j_0\) and the sums of one coordinate from
each order at the integer marks. Subtracting the single
directions gives the full basis at the other order.

At any particular \(h\) under consideration, each
prospective central slope at order \(j\) is the secant
between two of its internal indices. Its expression
in the core valuations is a nonzero linear form in
the adjacent valuations of that order plus the fixed
tail correction for \(h\). Equality of a prospective
lower-order slope and a prospective upper-order slope
therefore defines a proper affine hyperplane in the
variables \(x_u\), by surjectivity onto the two blocks
of coordinates. This holds also after reversal, which
negates and reverses the adjacent valuations.
Include all such pairs of secants, for every specified
\(h\) and every pair assigned to \(p\). Avoiding these
hyperplanes suffices regardless of which secants are
actual lower edges and regardless of the unit residues.

If there are \(B\) hyperplanes, one can choose integer
values \(x_u\in\{0,\ldots,B\}\) avoiding all of them.
Indeed a proper affine hyperplane contains at most
\((B+1)^{m-1}\) points of the \(m\)-dimensional box,
by solving for one coordinate having nonzero coefficient.
The union contains fewer points than the full box.
This gives an effective finite avoidance computation.

## Realizing all primes and the all-degree conclusion

Use the largest of the required mark counts to set
a common list of row labels at all primes, including
one row of each orientation per label. At a prime
where a label is unused, choose both bases zero
modulo a sufficiently large prime power and include
their fixed valuations in \(y_0\). The prescribed
depth differences in (1) are realized at every
prime by choosing two high individual depths with
the required difference and arbitrary unit parts.
Apply the Chinese remainder theorem separately to the
base of each row. The resulting ordinary smaller
exponents can all exceed any prescribed lower bound.
The valuation avoidance is exact, while the total
curvature of these finitely many gap-\(D\) rows tends
to zero. Their endpoint-normalized coefficients and
relative center corrections also tend to those of
the empty addition, by the shared-row curvature and
center bounds. Their signed gap is identically zero.

For the prepared core \(C'\), every aligned unit-gap
tail with a specified allowed \(h\) and sufficiently large
smaller orientation count \(S\) has no shared central
slope at each tested pair. At its corresponding prime,
the exterior slopes are already separated by the
aligned unit-gap theorem. There can then be no common
root at that pair. Compute \(S_0\) by taking the maximum
of the central theorem's explicit thresholds for the
finitely many pairs and allowed \(h\)'s, including reversal.
The first pair, terminals at least \(J\), and signed
gaps outside the selected menu are covered by the
ordinary positivity and absolute theorems. This proves
the stated all-degree regularity.

The construction excludes exact ordinary tails of this
residue class even near real Gaussian contacts. Its
constraint is the common congruence class of the small
endpoints and the unit gaps of the appended tail. Arbitrary
ordinary tails are outside this family.

