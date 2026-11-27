# Virtual norms force power values at their divisor fibers

A fixed rational function on a quadratic curve can yield
virtual elements at many specializations only if its
geometric divisor is divisible by the same power.
[CK4](MARKED_CUBIC_BASE_CHANGE_RATIONAL_THREE_TORSION.md#virtuality-on-a-positive-proportion-of-parameters)
detects the obstruction at depth one. Tracking the full
valuation identifies exact power conditions on binary
forms. This gives a power-saving bound, applies inside
genuine prime progressions, and sometimes gives finiteness.

## The full local condition

**Proposition VB1.** Let \(C\) be the smooth projective
model of \(Y^2=f(t)\), with \(f\in\mathbb Q[t]\) squarefree
and nonconstant. Fix \(a\in\mathbb Q(C)^\times\) and an
integer \(m>1\). For each closed point \(b\) of the
parameter line below a zero or pole of \(a\), let
\(\Phi_b(A,B)\in\mathbb Z[A,B]\) be its primitive
irreducible homogeneous form, with \(\Phi_\infty=B\).
Put
\[
 d_b=\operatorname{lcm}_{Q\mapsto b}
                  \frac{m}{\gcd(m,\operatorname{ord}_Q(a))},
 \tag{1}
\]
where the points \(Q\) are geometric. There is a fixed
finite set of rational primes \(S\) with the following
property. For coprime integers \(A,B\), \(B>0\), set
\(s=A/B\) and \(K_s=\mathbb Q(\sqrt{f(s)})\), and
assume this is quadratic and the evaluations of \(a\)
are defined and nonzero. If \(a_s\) has valuations
divisible by \(m\) away from \(S\), then, for every
\(b\) with \(d_b>1\),
\[
                 \Phi_b(A,B)=e_b z_b^{d_b}
                 \quad (z_b\in\mathbb Z\setminus\{0\}),
 \tag{2}
\]
with \(e_b\) in a finite set of signed integers supported
on \(S\). When \(m\) is odd, (2) for all such \(b\) is
also sufficient for virtuality away from \(S\).

Enlarge \(S\) to separate the closures of the
geometric support points of the divisor and the branch
fibers in integral models, and to exclude bad leading
coefficients of the finite \(\Phi_b\), their polynomial
discriminants and the prime two. In this construction
one can choose equations over the residue fields of
the points and remove the primes below all excluded
places. Away from \(S\), if \(p\mid \Phi_b(A,B)\) for
finite \(b\), then \(B\) is a unit. The point \(s\)
reduces to one simple root of \(\Phi_b(t,1)\), which
lifts to a root \(b_0\in\mathbb Q_p\). All other factors
are units and
\[
  \nu=v_p(\Phi_b(A,B))=v_p(s-b_0).
\]
For \(b=\infty\) take \(w=1/t\) and
\(\nu=v_p(B)=v_p(w(s))\). Over a given \(b_0\), the
geometric points have the same list of orders as in
(1), because they are obtained by embedding the
residue field of \(b\).

At an unramified point above \(b_0\), the local
coordinate \(t-b_0\), or \(w\) at infinity, has order
one. The two points in that fiber are either split or
defined over the unramified quadratic extension at
\(p\). The specialized quadratic algebra has the same
type since the normalized square root equation has
unit right side with the same reduction. On the
corresponding branch the normalized valuation of
\(a_s\) is \(n\nu\), where \(n\) is the divisor order.
This applies to each branch, including when the
unramified quadratic place is unique. In the split
case the conjugate evaluation interchanges the places.

At a ramified point above \(b_0\), the base uniformizer
is a square of a curve uniformizer times a unit.
The local specialization is quadratic ramified for
\(\nu\) odd, and split or unramified quadratic for
\(\nu\) even. This follows from the simple zero of the
local right side, using the reciprocal chart at an odd
degree infinity. Thus the normalized valuation of
\(a_s\) is \(n\nu\) for \(\nu\) odd and \(n\nu/2\)
for \(\nu\) even. In all cases divisibility by \(m\)
implies \(m/\gcd(m,n)\mid\nu\), and for odd \(m\) this
is equivalent. Away from the support closures the
evaluation is a unit. In the odd case, at fibers with
\(d_b=1\) the divisor orders are multiples of \(m\)
and the valuations are again divisible by \(m\).
This proves both local assertions. Taking the sign
and the exponents at \(p\in S\) modulo \(d_b\) gives
the finite set for \(e_b\).

In particular, if a rational exceptional fiber is
the zero of one of the primitive homogeneous linear
forms in a genuine prime constellation, an output on
that form is a prime to exponent one. At such a prime
outside \(S\) it violates (2). For example on a
symmetric line of [SL1](ABC_SYMMETRIC_PRIME_LINES_CYCLIC_OUTER_FIELDS.md)
the three directions are \(t=0,1/r,-1/r\).
No choice of local sheet can remove this obstruction
when an evaluation must work at every place of its
quadratic field.

## A power saving in rational parameter height

**Corollary VB2.** If one geometric multiplicity of
\(\operatorname{div}(a)\) is nonzero modulo \(m\), the
number of quadratic parameters \(s=A/B\) of height
\(\max(|A|,B)\le H\) where \(a_s\) has valuations divisible
by \(m\) outside any fixed finite set is
\[
                            O_{f,a,m,S}(H^{3/2}\log H).
 \tag{3}
\]
If such a divisor point lies above a rational parameter
fiber, the same count is \(O(H^{1+1/d_b})\).

Apply VB1 with the exceptional set enlarged to the given
finite set. Fix one \(b\) with \(d=d_b>1\), and a possible
\(e_b\). The polynomial \(Z^d-\Phi_b(A,B)/e_b\) is
irreducible over \(\mathbb Q(A,B)\), by Eisenstein at
\(\Phi_b\). It is separable in characteristic zero and
the generic Galois group is transitive on its roots.
A specialized rational root gives a proper
specialization group. The zero locus of \(\Phi_b\)
contributes at most \(O(H)\) pairs and the other
parameters are covered by the integral two-variable
Cohen–Serre bound, as recorded in
[Zywina, Theorem 1.2](sources/zywina-quantitative-hit-2010.pdf).
No regularity of the splitting field is required.
Summing over the fixed finite set of \(e_b\) proves
(3). For a linear \(\Phi_b\) it suffices to count its
\(O(H^{1/d_b})\) possible power values and \(O(H)\) integer
pairs on each value line.

## When the constraints give only finitely many parameters

**Corollary VB3.** For \(m=3\), let
\(\Phi_1,\ldots,\Phi_j\) be distinct forms in VB1 with
\(d_b=3\). Suppose there exist \(e_i\in\{1,2\}\)
such that \(\sum_i e_i\deg\Phi_i\equiv0\bmod3\),
and the total number \(r=\sum_i\deg\Phi_i\) of geometric
fibers is at least four. Then \(a_s\) is virtually cubic
away from a fixed finite set at only finitely many
rational \(s\). In particular four rational exceptional
fibers suffice. Three rational exceptional fibers
instead constrain the parameters to finitely many
fixed genus-one covers of the line, giving an upper
bound of the form \(O((\log H)^B)\) for fixed data.

In (2), multiply the equations with exponents \(e_i\).
Because the product has degree divisible by three,
dehomogenization and division by a power of \(B\)
give
\[
             v^3=\delta\prod_i\Phi_i(s,1)^{e_i}
\]
for a finite set of \(\delta\in\mathbb Q^\times\).
The degree-three cover of \(\mathbb P^1\) is geometrically
integral and totally ramified at the \(r\) stated
geometric points, including infinity if one form is
\(B\). Riemann–Hurwitz gives genus \(r-2\).
For \(r\ge4\), [Faltings, Satz 7](sources/faltings-1983-endlichkeitssaetze.pdf)
gives finiteness on each fixed smooth projective model. For four rational forms one
takes exponents \(1,1,2,2\). For three rational forms
one takes \(1,1,1\). Each resulting genus-one model
with a rational point is an elliptic curve. From the
equation, the height of \(v\) is bounded by a fixed
power of \(H(s)\). An elliptic abscissa is a fixed
rational function of \(s,v\), so its logarithmic
height is \(O(\log H)\) away from finitely many
points. The canonical height differs from a fixed
multiple of this height by a bounded amount. The
Mordell–Weil lattice count in a ball of radius
\(O(\sqrt{\log H})\) gives the polylogarithmic bound.

For an example on the original marked octic, take
\[
                  a=\frac{R-110t^2}{R+110t^2}.
\]
Its norm is one, and
\[
 F(t)-(110t^2)^2
      =(t^2-4)(t^2-9)(49t^4+1387t^2+1764).
\]
At \(t=\pm2,\pm3\) the two points have \(R=\pm110t^2\).
These roots of the displayed norm identity are simple.
The four fibers have divisor orders \(+1,-1\) on the
two sheets. Thus this additional rational norm-one
function can be virtual at only finitely many quadratic
parameters, even when any specified finite set of
primes is ignored.

## Applications to the existing field counts

**Corollary VB4.** The following numbers of known
marked directions exhaust all fixed functions with a
divisor divisible by three and cube norm on their
respective curves.

| family | marked directions |
|---|---:|
| The cubic base change \(t=3u^3/7\) of the marked octic, MO3 | 4 |
| The genuine symmetric prime line of each fixed \(r\equiv1\bmod6\), SL1–SL2 | 2 |
| The genuine cubic tangent-square progression, CTQ1–CTQ3 | 1 |

For each family and each finite list of fixed rational
functions on its curve, only \(O(H^{3/2}\log H)\)
parameters of height at most \(H\) can specialize any
of those functions to an additional independent
virtual cubic class. In particular the bound is
smaller than the established distinct-field lower
counts at their construction heights.

The first row is CK3. In the even sextic of SL1
the degree-two maps to the two elliptic quotients
are the two signs for \(t\mapsto-t\). As in CK2,
diagonal push-pull is two and off-diagonal
push-pull is zero, giving a prime-to-three isogeny
and thus an isomorphism on \(3\)-torsion. Each
elliptic quotient has a marked rational order-three
point and exactly one rational direction, by the
Weil pairing. The generic independence in SL1
therefore gives a basis for the rational
\(3\)-torsion of the sextic Jacobian. The norm
functions have only triple poles and zero orders
divisible by three, by the displayed norm identities
and the two distinct branches at \(t=0,1/r,-1/r\).
In CTQ1 the degree-four map between genus-one curves
gives an isogeny on Jacobians of degree four. The
target has a single rational \(3\)-torsion direction
by the same argument. The nonzero marked direction
of CTQ1 exhausts it.

On each curve, the divided-divisor map on the
cube-norm function classes is injective by the
constant norm argument of CK3. An arbitrary fixed
rational function with divisor divisible by three
has norm a constant times a rational-function cube.
At any virtual evaluation the rational norm must
be a cube, so the constant must be a cube. Such a
function is already in the marked span. All other
fixed functions are covered by VB2. Discarding the
finite specialization exceptions in the cube
identities has no effect on (3).

For MO4, the pairs in the proof of
[Stewart–Top, Theorem 2](sources/stewart-top-ranks-of-twists.pdf)
have height \(O(X^{1/24})\) for the degree-24 binary form
with nonzero discriminant. Passing to a primitive
pair only lowers height and preserves the parameter
and its square class. Thus the number of possible
exceptional fields through these height-bounded
witnesses is \(O(X^{1/16}\log X)\), versus the lower
bound \(\gg X^{1/12}/(\log X)^2\). In SL2 the prime
tuples used for the field count have height
\(O_r(X^{1/6})\), and in CTQ3 they have height
\(O(X^{1/4})\). The respective exceptional witness
counts are \(O_r(X^{1/4}\log X)\) and
\(O(X^{3/8}\log X)\), below the field lower exponents
\(1/3\) and \(1/2\) with logarithmic losses.
A field accounted for by one of these exceptional
witnesses contributes at most one to the comparison.
Hence the full lower counts survive deletion of the
exceptional witnesses for any fixed finite list.
These statements are at the prescribed construction
heights. They leave open extra classes from individually
chosen field elements or witnesses of other heights.
