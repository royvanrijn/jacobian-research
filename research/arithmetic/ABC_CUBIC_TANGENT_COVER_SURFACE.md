# The rational surface behind the cubic tangent coverings

The covering of the common tangent Selmer class has
a rational total surface. An explicit change of
coordinates turns the genus-one family into a pencil
on a split singular cubic surface. Its rational
parameters have three coprime denominator factors,
with cancellation governed by one fixed pure cubic
field.

## A simple pencil for the covering

Write \(t=p/h\),
\(D(t)=-3t^4+6t^2+1\),
\(u(t)=2t^3-3t^2-1\), and
\(n(t)=t(t-1)\). In the generic
[CS4 covering](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#the-explicit-genus-one-curve-on-the-third-twist),
with norm coefficient \(\nu=1/2\), the equation is
\[
 U^3+3u U^2 V+3D UV^2+uD V^3=4W^3 .
\]
Use the changes
\[
        X=U+uV,\qquad B=nV,\qquad A=X/2-tB .
\]
They are invertible for \(t\ne0,1\).
The norm identity \(u^2-D=4n^3\) gives the
depressed form \(X^3-12nXB^2+8uB^3=4W^3\).
Substituting \(X=2(A+tB)\) now gives
\[
                 W^3=2(A^3-B^3)+6t\,AB(A+B).            \tag{1}
\]
For the specialized original variables at
\((p,h)\), first replace \(V\) by \(h^2 V\)
on the normalized side. This absorbs the
degree-four field parameter and degree-two
coefficient ratio in CS4.

**Theorem TC1.** The tangent covering has a
rational total surface over \(\mathbb Q\).
Its generic genus-one fibration is birational
to the plane pencil \(Z=tB\) on the cubic
surface
\[
                   W^3-2A^3+2B^3-6ZA(A+B)=0
                   \quad\text{in }\mathbb P^3.            \tag{2}
\]
The surface has one split rational \(A_2\) point,
and the pencil base line meets it in the
closed point with residue field
\(\mathbb Q(\sqrt[3]2)\).

Equation (1) is the section \(Z=tB\) of (2).
Its base line \(Z=B=0\) meets (2) in
\(W^3=2A^3\), three distinct points in a
single cubic orbit. The line is disjoint
from the singular point \([Z:A:B:W]=[1:0:0:0]\).
To check uniqueness of that singularity,
the partial derivatives imply \(W=0\) and
\(A(A+B)=0\). For \(A=0\) they force
\(B=0\). If \(A=-B\ne0\), the
\(B\)-derivative forces \(Z=A\),
and the \(A\)-derivative is \(-12A^2\ne0\).
At the singular point the tangent cone
is the split form \(-6A(A+B)\).
The nondegenerate quadratic part in
\(A,B\), together with the independent
\(W^3\) term, gives type \(A_2\).
Projection from this point is given by
\([A:B:W]\) and solves for \(Z\)
rationally. This proves rationality.
Blowing up the three base marks and
resolving the double point gives a
smooth projective model of the
genus-one pencil. The generic fibre
still has period and index three
by CS4.

The example from [CS5](ABC_CUBIC_TANGENT_COMMON_SELMER_CLASS.md#an-explicit-failure-of-the-hasse-principle)
therefore has the simpler equation
\[
                     W^3=2(A^3-B^3)+17AB(A+B).           \tag{3}
\]
The rational change above transports
its local points, its point over the
pure cubic field and its failure
of rational solubility.

**Corollary TC2.** On a smooth projective
total model the Brauer group is the
image of \(\operatorname{Br}(\mathbb Q)\).
For any finite collection of completions,
parameters with a rational point on their
fibre are dense among simultaneous local
parameter choices at smooth soluble fibres.

The Brauer assertion is birational
invariance for smooth proper surfaces.
For density, a local point on a
smooth fibre has neighborhoods in
the local surface mapping into
arbitrarily small neighborhoods of
its parameter. Each meets the
birational plane chart by smoothness.
Weak approximation in that chart
gives a rational point of the surface
meeting all the neighborhoods.
This statement has no prime or
squarefree density hypothesis on
the rational parameters.

## Exact denominator restrictions at a prime gap

**Theorem TC3.** For any \(p\) and \(h=6k\)
in the squarefree tangent prime family,
a rational point on the covering
corresponds to a primitive integral
triple \((A,B,W)\) with
\[
       k W^3=2k(A^3-B^3)+p\,AB(A+B),\qquad
       \gcd(A,B)=1.                                    \tag{4}
\]
Conversely any primitive triple
satisfying (4) gives a rational
point of the covering. The three
integers \(A,B,A+B\) are nonzero
and pairwise coprime. Exactly one
of them is divisible by \(k\).
For every odd prime \(q\ne k\)
dividing one of them, \(2\) is
a cube modulo \(q\). If \(2\)
is a noncube modulo \(k\), the
valuation of the unique divisible
factor at \(k\) is exactly one.

For \(t=p/(6k)\) the invertible
change to (1) gives (4).
If a prime divides both \(A,B\),
the right side of (4) is
divisible by its cube, and its
valuation on \(kW^3\) forces it
to divide \(W\), also when the
prime is \(k\).
This contradicts primitivity.
None of the three factors can
be zero: on \(A=0,B=0,A+B=0\)
the equation for the quotient
of the remaining coordinates is
a cube of \(-2,2,4\), respectively,
which is impossible over \(\mathbb Q\).
Reduction modulo \(k\) and
\(\gcd(p,k)=1\) give its unique
factor. For \(q\ne2,k\) the
reductions on the three lines
again give cubes of \(-2,2,4\).
If \(k^2\) divides its factor,
divide (4) by \(k\) and reduce
modulo \(k\) to get the same
local test. Thus when \(2\)
is a noncube at \(k\) its
valuation must be one.

At \(2\) there is an additional
exact restriction. Since \(p,k\)
are odd, \(W\) is even. If
\(A\) or \(B\) is even its
valuation is one, because the
first term on the right of (4)
has valuation one and must cancel
the second. If both are odd,
\(A+B\) has valuation two.
Indeed \(v_2(A^3-B^3)=v_2(A-B)\)
and of \(v_2(A-B),v_2(A+B)\) one
is one and the other is at least
two. Cancellation of the smaller
valuation on the right of (4)
requires \(v_2(A-B)=1\) and
\(v_2(A+B)=2\).

These restrictions can be expressed
in the fixed pure cubic field
\(L=\mathbb Q(\sqrt[3]2)\).
Every prime admitted for
\(A,B,A+B\) has a prime above
it of residue degree one in
\(L\), except possibly \(k\).
This is the cube condition for
unramified primes and holds at
\(2,3\) by ramification.
The field has discriminant
\(-108\) and class number one:
the Minkowski bound is less
than three, and the prime over
two is generated by \(\sqrt[3]2\).
Thus all three factors are
signed norms of integral
elements of \(L\), after
dividing the unique \(k\)-factor
by \(k\) when \(2\) is a
noncube modulo \(k\).

This is an exact constraint on a
rational solution, including a
pairwise coprime additive triple
of pure cubic norms. An unbounded
fibre count still needs control
of their joint divisibility or of
smallest rational representative
heights as \(p,k\) vary.

## A power saving for small rational representatives

**Theorem TC4.** Fix \(0<\eta<1/6\).
Among marked prime-k pairs \(p,k\asymp Z\),
with \(p>6k\), the number having a rational
point of height at most \(Z^{7/6-\eta}\)
in the primitive coordinates \((A,B,W)\)
of (4) is
\[
               O_{\eta,\epsilon}
                    (Z^{\,2-(12/7)\eta+\epsilon}).
\]
For either sign and each sufficiently
large dyadic window \([X,2X)\), there
are \(\gg X^{1/2}/(\log X)^4\) squarefree
fields supplied with an everywhere
locally soluble marked tangent fibre
with no rational point up to a fixed
constant times \(X^{7/24-\eta/4}\)
in these coordinates.

Take the three linear forms
\(L_1=A,L_2=B,L_3=A+B\) and write
\(G_i=AB(A+B)/L_i\).
By TC3 there is exactly one \(i\)
with \(k\mid L_i(A,B)\).
Set \(Z_i=pL_i(A,B)/k\in\mathbb Z\).
Then
\[
                   W^3-2A^3+2B^3=Z_i G_i(A,B).          \tag{5}
\]
For each \(i\) this is a single
fixed geometrically integral
cubic surface in \(\mathbb P^3\).
Indeed it is linear in \(Z_i\)
and its coefficient \(G_i\)
is coprime over \(\overline{\mathbb Q}\)
to \(W^3-2A^3+2B^3\).
The surface has only one
singular point, of type \(A_2\),
by the same derivative test as
in (2), applied to each choice
of \(G_i\). In particular it
has only finitely many geometric
lines. The tuple in (5) has
projective height \(O(Z^{7/6-\eta})\),
since \(p/k\) is bounded and
\(|L_i(A,B)|\) is bounded by
twice the point height.
It recovers \(p,k\) from the
reduced quotient \(Z_i/L_i(A,B)=p/k\).
There is no summation over the
integer multiple \(L_i(A,B)/k\).

A rational line with vanishing
neither base coordinate is
vertical or maps with degree one
to the pencil \(Z_i/L_i\).
At finite smooth parameters
its fibre is the same plane
cubic (1), so a vertical
line cannot meet an eligible
smooth fibre in this chart.
A horizontal rational line
would give a section of the
genus-one pencil, contrary to
its generic index three.
The only rational points on
nonrational lines are the
finitely many possible
intersections with their
conjugates. These contribute
only boundedly many parameters.

Now apply the
[Salberger bound for cubic surfaces outside their lines](https://research.chalmers.se/en/publication/236437),
\(O_\epsilon(H^{12/7+\epsilon})\) at
projective height \(H\), to the
three fixed surfaces. This gives
the asserted power saving.
The sieve of CT4 gives
\(\gg Z^2/(\log Z)^4\) prime-k
pairs in either signed squarefree
box, and the bounded essential
representation orbits distinguish
\(\gg Z^2/(\log Z)^4\) squarefree
fields after removing the
power-small exceptional set.
With \(Z\asymp X^{1/4}\),
the height threshold becomes
the displayed field exponent.
The absence of a rational
point below a given height
does not exclude a point
of higher height.
