# A paired Newton-binomial certificate for consecutive angles

For two adjacent ordinary factorial moment polynomials there
is an angular certificate at a single common prime. Neither
endpoint difference has to be coprime to its degree. This is
useful even for horizontal polygons, including mirrored profiles.

## The paired criterion

At a prime \(p\), call a lower Newton edge binomial if only
its two endpoints attain its supporting line. Its width is
the difference of the two indices on that line.
All coefficients here are the original rational coefficients,
with no need to compute real endpoint normalizations.

**Theorem.** Fix \(n\ge2\) and \(N=n+1\). If the upper
Newton polygon consists of one binomial edge and all lower
edges are binomial with widths coprime to \(N\), then the
two ordinary moments have disjoint root arguments. In
particular, if both Newton polygons consist of one binomial
edge, the certificate applies for every prime \(p\), including
primes dividing one of the two degrees.

Embed the splitting field of both moments, taken in
\(\mathbb C\), into \(\overline{\mathbb Q}_p\). Normalize the
valuation by \(v(p)=1\), and use its algebraically closed
residue field. For an edge of slope \(s\), rescale \(z\) by
one algebraic element of valuation \(-s\) and divide all
coefficients by a scalar of the minimum valuation. The
reduction at this edge is a monomial times
\(u+v x^w\), for two nonzero residue coefficients and edge
width \(w\). A root with valuation \(-s\) rescales to a
unit. Its residue is a nonzero root of this binomial.
Thus for two roots with that valuation, their quotient
reduces to an element of multiplicative order dividing
\(w\). This statement concerns the quotient itself and so
is independent of the choice of the rescaling. If \(p\)
does not divide \(w\), distinct unit roots of the rescaled
polynomial have distinct residues. Indeed two unit roots
with the same residue would give a repeated root of the
reduction by factoring out both of their integral monic
linear factors.

All upper roots have the same valuation by the Newton
polygon root-count rule. Suppose there is a common quotient
of distinct roots, one pair from each degree. The common
quotient is a unit, so the two lower roots have the same
valuation and correspond to one of the lower edges. The
order of the common residue divides both \(N\) and that
edge's width \(w\). It is therefore one. Since \(N,w\) are
coprime, at least one is prime to \(p\). Its binomial
reduction has simple nonzero roots, which contradicts the
coincident residues of the two distinct roots in that degree.

If two nonreal roots have the same complex argument,
their quotients by their respective complex conjugates
are equal. The preceding argument excludes this.
Both polynomials have positive coefficients, and the even
factorial moment is strictly positive at every real
argument by the positive integral and distinctness of
the monomials. This excludes both real rays. Argument
separation also holds under independent positive
rescalings, and coprimality transfers to every
characteristic-zero field by a rational Bézout identity.

The certificate uses both exact arrays. The
[upper-only native-prime criterion](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md)
can use an upper primitive slope at a prime where the lower
polygon is unrestricted. The paired test dispenses with
primitivity by comparing the two reductions.

## A lower residue test with a nonbinomial edge

If \(p\nmid N\), the lower binomial assumption can be replaced
by an exact residue test. For each lower edge compute its
rescaled residual polynomial, dropping any initial monomial,
and let \(q(x)\) be the squarefree part. Require that
\[
                 \operatorname{Res}_x(q(x),y-x^N)
\]
be squarefree in the residue field. Equivalently, the map
\(\xi\mapsto\xi^N\) is injective on the distinct nonzero roots
of \(q\). If an upper quotient equalled a lower one, the
lower roots would belong to one edge and the quotient of
their residual roots would be an \(N\)-th root of unity.
The upper binomial is separable because \(p\nmid N\), so
that root of unity would be nontrivial. The two distinct
lower residues would therefore have equal \(N\)-th powers,
contrary to the test. No lower residual separability is
required before taking the squarefree part.

For example, use the single ordinary mirrored pair
\((9,10),(10,9)\), \(n=4\), \(N=5\), \(p=23\).
At the lower order the unscaled residual polynomial is
\(11x^4-3x^3+8x^2-3x+11\) over \(\mathbb F_{23}\).
The upper residual is a nonzero scalar times \(1+x^5\),
since both first blocks hit \(5\cdot9+1=46\).
The monic residue image polynomial is
\[
                       y^4-4y^3-3y^2-4y+1.
\]
Its gcd with its derivative is one modulo \(23\).
The image polynomial follows by a resultant with
\(y-x^5\), so the certificate applies although the
lower edge has all its interior points on the chord.

## Marking first blocks at both orders

Write \(\ell\) for a smaller row exponent and \(d>0\) for
the absolute gap. At order \(m\) remove the scalar
\((m\ell)!\) in that row. Its remaining factor at index \(k\)
uses the consecutive factors starting at \(m\ell+1\) and
ending at \(m\ell+d k\) or \(m\ell+d(m-k)\).
Take \(p>m\) and assume no factor in any row's range
\(m\ell+d+1,\ldots,m\ell+md\) is divisible by \(p\).
Let \(E,F\) be the total first-block valuations in increasing
and decreasing rows respectively. Then the normalized
coefficient valuations are
\[
        F \quad(k=0),\qquad
        E+F \quad(0<k<m),\qquad
        E \quad(k=m).
\]
The binomial coefficient is a unit. Every interior
ordinate is strictly above the endpoint chord if \(E+F>0\).
Thus a prime with only first-block hits and at least one
hit in each order gives the paired certificate.

For a horizontal example take the two mirrored pairs
\((7,8),(8,7)\) and \((23,24),(24,23)\) at orders four and
five with \(p=29\). The first pair has its only hits at
order four in \(4\cdot7+1\), and the second pair has its
only hits at order five in \(5\cdot23+1\). Both endpoint
differences in the valuations are zero.

## All-degree applications for movable rows

Let \(J\ge4\) be any integer at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Its numerical value is unevaluated. Fix any finite collection
of ordinary rows, whether or not its two moment sequences
are known to be coprime. Independently movable rows below
have smaller exponents that may be any positive integers
and fixed positive gaps and orientations.

**Two rows of arbitrary prescribed gaps.** With two
independently movable rows of prescribed gaps \(d_1,d_2\ge1\)
and either orientation for each, there are two simultaneous
arithmetic progressions of the smaller exponents for
which every profile has disjoint consecutive arguments in
every degree. Both added contributions to curvature can
be arbitrarily small. One can also move two mirrored pairs
at once, using the same smaller exponent in each pair.

For each \(2\le n<N=n+1<J\) take a distinct prime \(p\)
larger than all fixed factorial arguments in these orders
and larger than \(nN\max(d_1,d_2)+N\).
Require the smaller exponents to satisfy
\[
                n\ell_1+1\equiv0\pmod p,\qquad
                N\ell_2+1\equiv0\pmod p.
\]
At their designated orders the two rows each hit only
\(h=1\). A hit in the other order would require
\(nh\equiv N\pmod p\) for the first row or
\(Nh\equiv n\pmod p\) for the second. The bound on \(p\)
makes these equalities over the integers, and each is
impossible. Thus the paired certificate applies. Its
calculation is unchanged if one or both moving rows are
replaced by mirrored pairs. Chinese remaindering gives
the two simultaneous progressions over all the indicated
orders. No arithmetic input is needed at \(N=2\), by
strict positivity of the quadratic on the real line,
or at \(N\ge J\), by the absolute angular theorem.
The added curvature tends to zero as the two bases
tend independently to infinity in their progressions.
In particular opposite unit gaps preserve the original
signed total degree difference exactly.

**One row or mirrored pair of larger gap.** For a single
movable row, or a single movable mirrored pair, of fixed
gap \(d\ge J-1\), the same conclusion holds at a density-one
set of smaller exponents, including within every fixed
arithmetic progression of positive smaller exponents.

Choose a bound \(P_0\) larger than all the fixed arguments
below \(J\) and than \(J(d+1)\). It suffices that
\(\ell+1\) have a prime divisor \(p>P_0\). For every order
\(2\le m<J\) there is exactly one hit, at
\(m\ell+h\) with \(h=m\le d\), in the moving row or
in each row of the moving pair. Every fixed row is a
unit after stripping its scalar. Thus this same prime
certifies every adjacent pair below \(J\).
The exceptions have \(\ell+1\) supported on primes at most
\(P_0\), giving \(\mathrm{O}((1+\log H)^{\pi(P_0)})\)
exceptions up to \(H\) by counting prime exponents. This
is \(o(H)\) also relative to any fixed progression.
The conclusion includes mirrored families when the
fixed part is mirrored, and approaches arbitrarily small
added curvature without any angular approximation.
