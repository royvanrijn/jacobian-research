# Local coefficient interpolation by balanced ordinary rows

Fix an ordinary two-monomial background, with adjacent moment
coefficient ratios in the original physical variable. The background
can contain arbitrary rows, including private rows. At every prime
\(p>nN\), an ordinary unit-gap addition with zero signed gap can
prescribe all the adjacent coefficient ratios at orders \(n,N=n+1\)
to arbitrary \(p\)-adic unit precision and with exactly prescribed
valuations. The original background stays in this statement. The
addition needs only \(n+N-1\) rows of each orientation, and all
their smaller exponents can exceed any given bound.
Allowing a common larger fixed gap \(D\) extends adjacent-pair
interpolation to every prime, with \(n+N\) rows of each orientation.

The prescriptions can be simultaneous at any finite set of primes.
The pair can differ with the prime. One list of labelled rows works,
with the maximum of the required counts in each orientation. This
generalizes the
[four-row exterior matching](MULTIVARIATE_FOUR_ROW_PADIC_EXTERIOR_OBSTRUCTION.md),
which concerns the two extreme roots of a mirrored profile.

## The coefficient identity and the finite local graph

For a finite set \(\mathcal J\) of orders, denote the background's
adjacent ratios by \(B_{j,k}=c_{j,k}/c_{j,k-1}\),
\(j\in\mathcal J,\ 1\le k\le j\). All are positive rationals.
Adding increasing unit-gap rows of smaller bases \(x_i\) and
decreasing ones of bases \(y_i\), in equal counts, changes that
ratio by the exact multiplier
\[
 \frac{\prod_i(x_i+r)}{\prod_i(y_i+s)},\qquad
                     r=k/j,\quad s=(j-k+1)/j.                 \tag{1}
\]
The powers of \(j\) from the individual ramp factors cancel
because the counts are equal. In particular (1) does not require
the two orientations to share bases.

Make a bipartite graph on two copies of
\(\mathcal R=\{k/j:j\in\mathcal J,\ 1\le k\le j\}\).
The coordinate \((j,k)\) is an edge from left \(r\) to right \(s\).
A system of nonzero target multipliers \(Z_{j,k}\) over a field
has the form \(L_r/R_s\) if and only if the alternating product
of its edge values is one on every cycle. Necessity is
telescoping. For sufficiency choose a value at one vertex
in each component and propagate along a spanning tree.
The remaining edges obey the targets by the cycle condition.

For adjacent \(n,N=n+1\) with \(n\ge2\), the nonterminal fractions
of the two orders are distinct by coprimality. The only repeated
vertex in either copy is \(1\), which occurs twice.
The two nontrivial components are the two stars at \(1\) in the
respective copies. No edge joins their centers, since \(k=j\)
and \(k=1\) cannot both hold for \(j\ge2\). All other edges are
isolated. In particular there are no cycle conditions and
\(|\mathcal R|=n+N-1\). For example, the larger order set
containing \(3,4,6,12\) has a cycle giving
\(Z_{12,3}Z_{3,1}=Z_{6,2}Z_{4,1}\).

For the more general order set assume that \(p\) makes the
denominators in \(\mathcal R\) units and makes the distinct
elements of \(\mathcal R\) distinct modulo \(p\). The bound
\(p>(\max\mathcal J)^2\) suffices, while \(p>nN\) suffices for
an adjacent pair. Target ratios \(\widehat B_{j,k}\in
\mathbb Q_p^\times\) are admissible when the quotients
\(\widehat B_{j,k}/B_{j,k}\) satisfy the cycle condition.
Choose corresponding \(L_r,R_s\) in \(\mathbb Q_p^\times\).
In each graph component multiply all its \(L_r,R_s\) by
the same power of \(p\) so they all have valuation at least
the desired precision \(L\ge1\). This does not change the
target quotients.

Use one increasing row per left vertex \(r\), prescribing its base by
\[
       x_r\equiv-r+p^{e_r}u_r\pmod{p^{e_r+L}},
                   \qquad e_r=v_p(L_r)\ge L.
\]
Here \(u_r\) is a unit to be chosen. All congruences for rational
fractions are interpreted using their unit denominators.
For any \(r'\in\mathcal R\), the factors \(x_r+r'\) with \(r\ne r'\)
are units congruent to \(r'-r\) modulo \(p^L\). The single
remaining factor has valuation \(e_{r'}\) and unit
\(u_{r'}\) modulo \(p^L\). Choose
\[
          u_{r'}\equiv
             \frac{p^{-e_{r'}}L_{r'}}
                        {\prod_{r\ne r'}(r'-r)}\pmod{p^L}.
\]
Then \(\prod_r(x_r+r')/L_{r'}\equiv1\pmod{p^L}\).
The decreasing bases \(y_s\) are chosen independently by the
same formula for the targets \(R_s\). Equation (1) now gives
\[
       \frac{(c_{j,k}/c_{j,k-1})_{\mathrm{augmented}}}
            {\widehat B_{j,k}}\equiv1\pmod{p^L}
\]
at every edge. This congruence includes the exact equality
of valuations. The same construction gives exact interpolation
by bases in \(\mathbb Z_p\): after fixing the initial residues
of the units, a change by \(p^q\) in a unit changes its own
normalized product at order \(p^q\) by a unit times that change
in the digit. It changes all other normalized products only
at order at least \(p^{q+1}\). The diagonal equations lift
one digit at a time.

For a finite list of primes and order sets, reserve as many
row labels per orientation as the largest \(|\mathcal R|\).
At a prime requiring fewer labels, give all unused bases
the class \(0\bmod p^L\). Include their known unit contributions
in the background. These contributions have the form (1),
so they preserve the cycle condition.
Apply the Chinese remainder theorem separately to every label.
Each smaller exponent can then be taken arbitrarily large
in its prescribed residue class, independently of all the others.

## A fixed larger gap for arbitrary primes

Fix adjacent orders \(n,N=n+1\) with \(n\ge2\), any prime \(p\),
and a prime gap \(D>N^2\). Use one list of \(n+N\) marks \(u\in
(0,D)\cap\mathbb Z_p\), with rational values. One row in each
orientation will be associated with every mark. For a mark \(u\),
let \(v_u\) be the vector with one in coordinate \((j,k)\) if
\(ju\in\mathbb Z\) and \((k-1)D<ju<kD\), and zero otherwise.
Choose all marks avoiding the two orders' block boundaries.

To obtain the marks, choose an order \(j_0\) prime to \(p\)
and let \(j_1\) be the other order. For each block \(k\) of
\(j_0\), choose an integer \(h\) strictly inside
\(((k-1)D,kD)\) and not divisible by \(j_0\). Put \(u=h/j_0\).
It is integral at \(p\), and \(j_1 u\notin\mathbb Z\) since
\(\gcd(j_0,j_1)=1\). This gives the single coordinate direction
at every block of order \(j_0\). For each block of \(j_1\),
choose an integer \(u\) strictly inside
\(((k-1)D/j_1,kD/j_1)\). Such an integer exists because
\(D>j_1\), and \(0<u<D\). The prime \(D>j_0\) ensures that
\(j_0u\) is also strictly inside a block. The direction is
then the sum of the desired coordinate of \(j_1\) and one
of the available order-\(j_0\) coordinates. The resulting
\(n+N\) vectors are an integer basis, by subtracting the
already available directions.

Take the increasing row at a smaller base
\(x=-u+p^{H_+}z_+\) and the decreasing row at
\(y=-(D-u)+p^{H_-}z_-\) in \(\mathbb Z_p\), with unit
\(z_\pm\). The reflected factor in the decreasing row
belongs to the same adjacent ratio as the marked factor
in the increasing row: an offset \(ju\) strictly in block
\(k\) reflects to \(jD-ju\) strictly in block \(j+1-k\).
The increasing factor contributes \(j p^{H_+}z_+\)
and the decreasing factor divides by \(j p^{H_-}z_-\).
Their net marked ratio is \(p^{H_+-H_-}z_+/z_-\) in
each coordinate of \(v_u\).

Every unmarked factor is a nonzero rational offset
\(a=h-ju\) or \(a=h-j(D-u)\) plus \(j p^{H_\pm}z_\pm\).
For all \(H_\pm\ge L+\max v_p(a)\) over these offsets,
it has the valuation of \(a\) and agrees with \(a\) to
relative precision \(p^L\). Divide the desired ratio array
by the background and the fixed unmarked contributions.
The integer basis \(v_u\) gives an exact multiplicative
representation of the resulting vector using a value in
\(\mathbb Q_p^\times\) for each mark. First assign
\(H_+-H_-\) to match its valuation while keeping both
depths sufficiently large. Then choose \(z_+/z_-\) with the
prescribed unit residue. These conditions give ordinary
positive bases by reducing them modulo a sufficiently high
power of \(p\) and lifting those classes to integers.

The same prime gap \(D\) can serve any finite list of
pairs by taking \(D\) larger than the square of the largest
terminal. To use one list of labels at all prescribed
primes, pad pairs with smaller label requirements by bases
congruent to zero at sufficiently high precision. Their
sampled unsigned factors agree with the fixed offsets
\(h\) to relative precision \(p^L\), so include their
contributions in the background. The Chinese remainder
theorem combines the resulting prescriptions. All labels
are balanced, including the padding, and their bases may
again increase independently by multiples of the moduli.

## Forced large valuation spreads at matched exterior roots

At a prescribed prime for a pair \(n,N=n+1\), fix any positive
integers \(E,L\). Use the adjacent ratios of the following target
polynomials:
\[
 \widehat F_j(z)=
     (z+p^E)(1+p^E z)\sum_{k=0}^{j-2} z^k .
                                                               \tag{2}
\]
The endpoints have valuation \(E\) and all interior coefficients
are units at every prime. The two exterior roots are
\(-p^E,-p^{-E}\), and the two zero-valuation residue parts after
removal of their monomials are proportional to the consecutive
geometric sums \(S_{n-2}\) and \(S_{n-1}\), where
\(S_b(z)=\sum_{k=0}^b z^k\). They are coprime in any residue
characteristic by \(S_{n-1}-z S_{n-2}=1\).

The factorial pair has the same ordinates up to overall translations
and the same reduced unit ratios modulo \(p^L\).
Indeed normalize each transformed polynomial by a scalar
so its constant coefficient equals that of (2), and multiply
the successive ratio congruences. This compares every coefficient
with its counterpart in (2), at relative precision \(p^L\).
Consequently at each exterior valuation the factorial polynomials
have a unique root. After rescaling \(z=p^E w\) or taking the
reciprocal variable and making the same rescaling, the two
normalized roots agree with \(-1\) modulo \(p^L\).
The normalized polynomial has a simple residue root there,
so this follows by the elementary single-digit lifting for
a simple root. The two horizontal residual polynomials
remain coprime. Thus at all precisions up to \(L\) the
[physical truncated-face tests](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)
leave precisely the two exterior branches, each of width one.

Any rational factor surviving these tests at this prime must
be quadratic if it is a nonconstant factorial common factor.
Its two valuations must be \(E,-E\), since the even moment is
strictly positive on the real line and its rational common factors
have even degree. The
[trace denominator identity](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md#the-total-depth-and-real-factor-bound)
then contributes \(2E\) to \(v_p(D_H)\).
Its reciprocal trace has valuation \(-2E\). This prime supplies
no integral trace-cluster depth, although the exterior normalized
roots at the two orders match to arbitrary specified precision.
The conclusion concerns permitted local allocations and does not
assert that a rational common factor exists.

## Arbitrary real diffuse limits with these local prescriptions

Fix a finite ordinary core, a signed integer gap \(h\) for
the full added profile, and a desired \(\tau>0,\ |\lambda|\le\tau\).
Before the local interpolation, add the same number \(t\) of
unit rows in each orientation. For \(T_\pm=(\tau\pm\lambda)/2>0\),
choose the smaller exponents in that orientation to be
\[
    2\bigl(\lfloor A_\pm t\rfloor+i\bigr)+b_\pm,\qquad
    1\le i\le t,\quad
    A_\pm=(e^{2T_\pm}-1)^{-1},\quad b_+=0,\ b_-=1 .
\]
Their respective curvature sums tend to
\(\tfrac12\log((A_\pm+1)/A_\pm)=T_\pm\) by Riemann sums.
If \(T_\pm=0\), take bases of order \(t^2\) instead.
The parity distinction keeps the bases of the two orientations
disjoint. Add \(|h|\) correction rows of the required sign
at still larger distinct bases. This gives signed gap \(h\).
Such rows contribute \(o(1)\) to the curvature and signed bias.

Treat this current profile including the core as the background.
Its valuations at the prescribed primes are unrestricted:
the interpolation adjusts their exact values. Choose the new
local row bases so large that their curvature contribution
tends to zero and their maximum relative gap is no larger
than that of the earlier added rows. This keeps
\(\eta\asymp1/t\) and makes the total added curvature and bias
tend to \(\tau,\lambda\), respectively. The number of distinct
added types diverges. The primes, \(E\)'s and precision depths
\(L\) may vary with \(t\), even after the earlier added rows
have been fixed. In particular \(E\log p/\log(1/\eta)\) can
grow arbitrarily fast. Any finite list of adjacent pairs,
including all pairs below a fixed absolute terminal cutoff,
can be prescribed at distinct primes simultaneously.

The same examples can have disjoint consecutive root arguments
in every degree. Fix an [absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
\(J\), and reserve a second distinct prime \(q>nN\) for every
pair \(2\le n<N=n+1<J\). At its two orders \(j\), prescribe the
ratios of \(1+q\sum_{k=1}^{j-1}z^k+z^j\).
These extra conditions coexist with all the original primes.
They give reductions \(1+z^j\) after constant normalization.
Because \(q>N\), the reductions are squarefree. Quotients
of distinct roots reduce to nontrivial \(j\)-th roots of unity,
so the coprime degrees give disjoint quotient sets. A common
nonreal root argument would give a common quotient by taking
each root divided by its complex conjugate. A common real
argument is excluded by the even moment's strict positivity.
That positivity also covers the first pair, and the cutoff
covers the remaining degrees. Thus the permitted local
quadratic allocation can occur at any designated primes even
on profiles that are regular in all degrees.

This construction rules out a universal \(O(\log(1/\eta))\)
upper bound for the locally permitted valuation spreads based
only on the diffuse parameters and bounded-precision checks at
the prescribed primes. It applies even when the fixed core and
limiting real strengths and bias are specified in advance.
A separate prime, a sufficiently deep profile-dependent check,
or an arithmetic relation across primes can still exclude a
candidate global factor.

## One rational quadratic compatible with a prescribed finite screen

Suppose at the chosen real diffuse limit both Gaussian
references have a common nonreal root in the same physical
variable. Given any prescribed finite screen of primes and
precisions, the same construction admits a single rational
quadratic whose conjugate pair is arbitrarily close to the
reference pair at the real place, and whose two exterior roots
agree with the factorial roots at each screened prime to
the prescribed precision. The quadratic can depend on the
screen and the profile. This also applies when the screen
contains every prime up to a bound chosen as a function of
the final maximum added relative gap.

For the proof, first choose the real diffuse background of the
previous section. Fix a positive rational common argument
scale tending, after normalization, to the lower endpoint
scale used in the
[physical reference](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
For example, the inverse of the lower tail's first coefficient
ratio, with its fixed Gaussian factor absorbed into the limiting
variable, gives such a rational scale by the exact comparison
of the two arrays. Now choose a monic rational quadratic \(H\)
so that in this common scale it approximates the monic
quadratic of the reference conjugate pair over the reals.
At each screened prime \(p\), independently require \(H\) to
approximate \((z+p^E)(z+p^{-E})\) in the original variable,
as closely as needed so its two normalized roots agree with
the target units. The real and finite-prime conditions are
simultaneously possible for both rational coefficients.
Indeed put both \(p\)-adic approximation conditions into
integer residue classes by clearing a fixed denominator.
The solutions can be divided by an arbitrarily large positive
integer that is one modulo all the resulting prime powers.
The resulting real mesh of the common residue class tends
to zero. Strict negativity of the real discriminant persists,
so this is one irreducible rational quadratic.

At each screened prime use the adjacent ratios of \(H(z)\)
times the consecutive geometric cofactors from (2).
These targets have nonzero coefficients when the local
approximation of \(H\) is sufficiently close, and they leave
precisely the two required exterior branches. Interpolate
these ratios by the ordinary balanced rows. The real
coefficient comparison can be kept at arbitrarily high
precision in the fixed rational scale. For this, after the
simultaneous residue classes for all the new bases have been
chosen, take every base at bounded distance from a common
large real height \(B\), letting \(B\) tend to infinity.
Each added factorial block at order \(j\) then consists of
\(D\) factors asymptotic to \(jB\), where \(D\) is the chosen
common gap, possibly one. Equal counts of the two orientations
make the quotient in each adjacent coefficient tend to one
over the reals. Thus the ordinary real pair stays near the
same reference in the same rational scale. The interpolation
rows can all have relative gaps below that of the preliminary
diffuse rows, so the screen may be chosen after this maximum
relative gap is already fixed.

The extra primes giving disjoint consecutive arguments can
still be chosen outside the prescribed finite screen. This
strengthens the limitation to tests that ask for simultaneous
compatibility of their finite data with one global rational
quadratic of unrestricted height. An exact common divisor,
a bound on its height, and tests that use all primes at once
are distinct requirements.
