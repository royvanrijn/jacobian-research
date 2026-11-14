# The compact Lie group Mathieu property, including disconnected groups

14 November 2026. This is a written deduction from the locally retained
[GVC2FTS finite-trace theorem](../../archive/non-elliptic/extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md)
and the [audited connected classification](COMPACT_GROUPS.md). The prime-digit
proof of GVC2FTS, including its number-field specialization, split-prime step
and Newton endpoint, was read again for this application. Its existing exact
regression passed. This is an internal proof, without external review or
formalization; no novelty or priority claim is made.

**Theorem.** Let G be a compact Lie group, H=G^0 its identity component,
and R(G) its complex representative-function algebra. The kernel of
normalized Haar integration on R(G) is a Mathieu--Zhao subspace if and only
if H is abelian, equivalently a torus (including the trivial torus).
In particular every finite group and every compact Lie group with toral
identity component has the property. Noncommutativity of the finite component
group, and nontrivial conjugation of H, do not change this conclusion.

## Restriction and zero extension

H is open, normal and of finite index N. For coset representatives g_i,
restriction by h -> g_i h gives an injective algebra map

    R(G) -> product_(i=1)^N R(H),
    I_G(f)=(1/N)*sum_i I_H(f_i).

It is also surjective. Here is the finite-type step, which cannot be
replaced by merely extending a continuous function. If
f(h)=lambda(rho(h)v) is a matrix coefficient on H, consider the induced
representation consisting of functions F:G->V with
F(hg)=rho(h)F(g). It has dimension N dim V and G acts by right
translation. The vector F_v supported on H, with F_v(h)=rho(h)v, is
continuous because H is clopen. The matrix coefficient obtained by evaluating
lambda(F_v(g)) at the identity is exactly f on H and zero off H.
Thus zero extension belongs to R(G). Its left translates and finite sums
realize arbitrary tuples on the cosets. Normalized Haar measure gives every
component mass 1/N by translation invariance. This proves the algebra and
functional identification, with no splitting of G over its component group.

## Positive direction

If H is a torus of dimension r, identify R(H) with
C[z_1^(+/-1),...,z_r^(+/-1)] and I_H with constant term. Let f have all
positive pure Haar moments zero. Its restrictions f_i then satisfy

    sum_i CT(f_i^m)=0 for every m>=1.

GVC2FTS forces CT(f_i^m)=0 for every i,m. Duistermaat--van der Kallen
and strict separation give CT(q_i f_i^m)=0 for all large m for each fixed
component q_i. There are finitely many components, so the maximum of their
thresholds works for q on G. Averaging proves the Mathieu--Zhao property.
For r=0 the same proof reduces to Newton identities for N scalar values:
all values are zero. This includes nonabelian finite groups.

The equal positive component weights matter. Weighted complex cancellation
from GVC2WTD is not an obstruction: Haar supplies precisely the equal trace.
Also use separate support separators on the different components. For
example f_+=z+z^2 and f_-=z^(-1)+2z^(-3) have zero pure moments componentwise,
but the union of their supports contains zero in its convex hull. Requiring
one global separator would wrongly discard this valid case. The new checker
covers this opposite-component example through m=11 with multiplier tuple
(z^(-7),z^8); all mixed terms vanish for m>8.

## Negative direction

If H is nonabelian, the connected theorem supplies p,q in R(H) with all
positive pure moments of p zero and infinitely many mixed moments nonzero
(in fact the source gives every positive index). Zero extend both functions
to G. The preceding induction argument proves they remain representative.
For m>=1,

    I_G(p_tilde^m)=(1/N)I_H(p^m)=0,
    I_G(q_tilde*p_tilde^m)=(1/N)I_H(q*p^m).

The second sequence retains its nonzero values, proving failure on G.
This completes both directions.

## Relation to the new literature

Liu's [Coherent-State Circuit Cancellation, v2.2](https://zenodo.org/records/21689760)
record describes a different connected construction and reduces the virtually
toral remainder to the identity-coefficient functional on C[Z^r x C_N].
Only its public record was inspected here; its PDF proof and producer code
were not audited. The finite-trace theorem already in this repository closes
that reduced remainder. Our proof above supplies the representative-function
identification directly and depends on Long for the connected negative case.
Minimal representation or polynomial degrees, non-Lie compact groups and
other coefficient algebras are outside the theorem.
