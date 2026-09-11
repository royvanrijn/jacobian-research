# The prescribed rank16 bridge fails integral gluing

Research-priority update: the [accessibility closure and independent11952
rank-bound target](ACCESSIBILITY_CLOSURE_AND_11952_RANK_BOUND_2026-09-11.md)
now govern the rank-jump investigation. The carrier and gluing results below
remain valid; this note is not the current experiment handoff.

**The retained determinant4100 common-core witness does not extend to the
marked K3 lattice.** Its embeddings in class6/embedding25 and
class1/embedding33 have different order-two gluing classes. The core has
only the automorphisms \(+1,-1\), so changing its identification cannot
repair this mismatch.

This excludes simultaneous geometric realization of **that particular
marked overlap**, even before effectivity, nefness or rationality. It
does not exclude class1, class6, other overlaps, or their existing rational
fibrations. The currently realized degree-two class1 fibration has a
different literal common core, of determinant13104.

For that existing realization, all fourteen seed-universality directions
transport exactly to distinct B-base parameters. Each transported point
adds an eighteenth independent direction to the specialized generic17.
None is in the rational span of generic B-sections. The fourteen B-fibres
give genus-one quadratic carriers over A; all91 joint normalizations have
genus5. More generally, no two distinct smooth members of this B-fibre
family have the same branch divisor.

Each tested carrier supplies a primitive height8 twist section. The full
twist ranks remain unknown. A separate automorphism-return theorem below
limits one possible propagation mechanism, without explaining the finite
13-of14 recovery experiment.

## The determinant4100 core obstruction

Let \(L_6,L_1\) be the positive definite rank17 frames in the retained
Niemeier packet, embeddings25 and33. Both have determinant1092. Their
common primitive rank16 lattice \(K\) has determinant4100 and the exact
common basis specified by the
[class1 marking](../../artifacts/generated-results/elliptic-curves/x1092_class1_realization_marking_v1.json).

In each frame, let \(v_i\) generate the integral orthogonal complement of
\(K\). Exact coordinates give

\[
v_i^2=1119300,\qquad
\operatorname{div}_{L_i}(v_i)=546,\qquad
[L_i:K\oplus\mathbf Z v_i]=2050.
\]

Although the abstract orthogonal complements agree, the integral
overlattices differ. Project \(L_i\) orthogonally to \(K^\vee\), and write

\[
H_i=\pi_K(L_i)/K\subset K^\vee/K.
\]

Both \(H_i\) are cyclic of order2050. Their unique elements of order two
are \(a_i/2\), in the common ordered basis of \(K\), where

~~~text
a_6 = (0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,1)
a_1 = (1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,1).
~~~

Thus \(H_6\ne H_1\). This is already an obstruction to extending the
literal common-core identification.

It also survives every possible isometry of \(K\). There are exactly
1032 norm4 vectors, or516 antipodal lines, spanning rank16. Color
refinement on their absolute inner products gives248 classes and then516
singleton classes. Every isometry therefore fixes every minimal-vector
line. The graph joining lines with nonzero inner product is connected,
forcing all their signs to agree. Consequently

\[
O(K)=\{+I,-I\}.
\]

The independent checker proves completeness of the norm4 enumeration by
an exact rational LDL sphere recursion. It does not use the producer's
short-vector enumeration to establish completeness, or assume a reported
automorphism-group order.

Now consider the stable lattices
\(N_i=U\oplus(-L_i)\). Adding \(U\), which is orthogonal to \(K\), does
not change \(H_i\). An integral isometry \(N_6\to N_1\) matching the two
embedded copies of \(K\) would restrict to an element of \(O(K)\) and
carry \(H_6\) to \(H_1\). Both possible restrictions fix each subgroup.
The displayed mismatch is a contradiction.

Hence the two retained core embeddings cannot be the same integral
rank16 sublattice in simultaneous marked frame realizations. This
statement concerns the specified primitive inclusions, not merely the
abstract frame types. The earlier
[class1 realization](X1092_CLASS1_RATIONAL_MW17_REALIZATION_2026-09-10.md)
explicitly did not assert that its low-degree transport extended this
common-core identification; its realization theorem remains valid.

Indeed, for the actual realization \(D_B=(3,2,w)\), with the old zero
retained, the common frame intersection is \(w^\perp\subset L_A\).
Here \(w^2=12\) and \(\operatorname{div}(w)=1\), so its determinant is

\[
1092\cdot12=13104.
\]

## Exact transport of the same fourteen directions

Use the original public rank31 subgroup \(D\), the specialized full
generic MW17 subgroup \(M\), and the fourteen literal public words from
the seed-universality input. Adjoining those words to the31-by17 generic
embedding gives determinant **−1**. Thus these are an integral basis of
\(D/M\), not a substituted mod2 complement or newly chosen easier points.

Their transport uses the realized class1 trace/Riemann–Roch pencil,
pointed quartic, pointed cubic, and compact polynomial parent. Every
source point remains on the original A-fibre \(t=0\); each receives its
own B-parameter \(s_i\). All fourteen \(s_i\) are distinct.

At each parameter, the17 saturated generic B-sections and the transported
point have independent finite mod2 images of rank18, together with a
modular certificate excluding rational2-torsion. Infinite descent proves
their independence over \(\mathbf Z\), hence over \(\mathbf Q\). Therefore

\[
P_i^B\notin
\operatorname{span}_{\mathbf Q}
\{S_1(s_i),\ldots,S_{17}(s_i)\},\qquad i=1,\ldots,14.
\]

This excludes arbitrary multiples/translates of generic rational
B-sections for these points, not merely a bounded coefficient search.
Each B-fibre has rank at least18; no exact rank is claimed. Since
fibrewise B-translations preserve the base, they also cannot put these
fourteen distinct-parameter points in one orbit.

The independent replay uses Python rational arithmetic for the original
public words, inverse maps for the transport, and exhaustively enumerated
finite elliptic groups instead of the producer's cubic Legendre characters.
The complete class1 realization checker also passes separately.

This is a result for the existing class1 realization with core13104.
It cannot be reinterpreted as a calculation on the incompatible4100
marked bridge. A point on A302 is a surface point, not a generic A-section;
the map between fibrations does not preserve the group law on fixed fibres.

## The branch map excludes shared smooth B-fibre covers uniformly

Fixing \(u=u_i\) in the class1 pencil gives an actual A-multisection

\[
C_i:\ v^2=q_i(t),\qquad \deg q_i=4.
\]

The carrier packet includes the exact rational functions

\[
x=x_0+x_1v,\qquad y=y_0+y_1v
\]

on the short A-model \(y^2=x^3+A(t)x+B(t)\), and verifies both coefficient
identities after reducing modulo \(v^2-q_i(t)\). The original point gives
the rational point \((0,v_i)\), with \(v_i\ne0\).

A single certificate modulo79 proves that all14 quartics are squarefree
and all91 pairwise resultants are nonzero. Thus their branch divisors
are disjoint, with no branch at infinity. Each \(C_i\) has genus1, and
Riemann–Hurwitz for the connected biquadratic cover gives

\[
g\!\left(\widetilde{C_i\times_{\mathbf P^1_A}C_j}\right)
=1+1+3=5\quad(i\ne j).
\]

These genus5 curves do have the known rational point above \(t=0\).
Their genus is a complexity statement, not an insolubility claim.

There is also a uniform result for the whole smooth B-fibre family.
Clear the common denominator of the five quartic coefficients \(c_i(u)\)
and remove their common polynomial content. The resulting polynomials
\(p_0(u),\ldots,p_4(u)\) all have degree4. Their five-by-five coefficient
matrix \(M\) is invertible, certified modulo17 and replayed by an exact
rational inverse. In homogeneous coordinates \(u=V/U\), the branch map is

\[
[U:V]\longmapsto
M[U^4:U^3V:U^2V^2:UV^3:V^4].
\]

This is a projectively transformed degree-four Veronese embedding.
Different parameters give different projective binary quartics.
For smooth members those quartics are their branch divisors, so no two
distinct members define the same quadratic extension of the A-base,
even up to a constant twist. Their joint normalization has genus at
least2. The91 tested pairs have the stronger genus5 conclusion.

The uniform statement excludes this smooth B-fibre route to a shared
genus0/1 cover. It does not cover singular-fibre normalizations, other
multisections, or additional sections on one fixed twist.

## Explicit primitive height8 sections on all fourteen twists

The deck involution on each A-carrier satisfies
\(\sigma(P)=T-P\), where \(T=P_w\) is the generic trace section.
This follows from the residual chord through \(-T\). Thus the subgroup
generated by \(M,P,\sigma(P)\) adds exactly one rank direction.

An explicit anti-invariant section \(R=P-\sigma(P)=2P-T\) is obtained by

\[
r=\frac{y_0^2}{x_1^2q}-2x_0,\qquad
z=\frac{y_0(x_0-r)}{x_1q}-y_1.
\]

It gives the rational point \(X=qr,\ Y=q^2z\) on the standard twist

\[
E_A^{(q)}:\quad Y^2=X^3+q^2A X+q^3B.
\]

For every tested \(q_i\), \(X,Y\) simplify to polynomials of degrees8,12.
A certificate modulo149 proves that the original discriminant has24
simple finite roots, that each \(q_i\) is disjoint from that discriminant,
and that \(\gcd(q_i,X_i)=1\). The twist therefore has24 fibres \(I_1\),
four fibres \(I_0^*\), and a smooth fibre at infinity. Its holomorphic
Euler characteristic is4.

The displayed section has no intersection with the zero section, and
meets the identity component at all four \(I_0^*\) fibres. The Shioda
height formula gives height8. Every nonzero section has height at least4:
the four \(D_4\) fibre corrections contribute at most1 each. This also
excludes torsion and shows that the displayed height8 section is primitive,
since a nontrivial division would have height at most2.

The geometric trivial lattice has rank18 and \(h^{1,1}=40\), giving

\[
1\le \operatorname{rank}E_A^{(q_i)}(\overline{\mathbf Q}(t))\le22.
\]

In particular the known rank-one subgroup is saturated, but the **full
twist rank is UNKNOWN**. A second independent section would be a real
shared-carrier result; translations or multiples of this one are not.

The invariant/anti-invariant rank decomposition and section-translation
automorphisms are standard; see
[Schütt–Shioda, §§7.6 and7.12](https://arxiv.org/html/0907.0298v3).
The polynomial sections and local fibre checks here are exact computations
on the frozen inputs.

## A limit on automorphism-based propagation

Let \(S\to\mathbf P^1\) be the A-fibration and \(C=E_{302}\) its smooth
fibre at0. Its \(j\)-invariant is neither0 nor1728. Let \(G\) be the
specialization of the full generic Mordell–Weil group.

**Theorem.** If a \(\mathbf Q\)-automorphism \(g\) of the K3 surface satisfies \(g(C)=C\),
then on \(C(\mathbf Q)/G\) it acts by \(P\mapsto P\) or \(P\mapsto-P\).

**Proof.** The complete linear system \(|C|\) is the elliptic pencil.
Preserving \(C\) preserves its divisor class and this pencil, inducing
a base automorphism fixing0. The curve \(g(O)\) is a rational section,
because it has fibre intersection1. Subtract that section from \(g\).
On \(C\), the result is an automorphism fixing the zero. Since the
\(j\)-invariant is not special, it is \(+1\) or \(-1\). Restoring the
section changes the result by an element of \(G\). ∎

Consequently, any fixed word in the A/B section-translation automorphisms
that preserves the entire302 fibre cannot amplify one seed into independent
quotient directions. If a fixed word does not preserve that fibre, its
returns lie in the finite intersection \(C\cap g^{-1}(C)\).
A finite collection of such words therefore cannot give uniform
amplification on a Zariski-dense set of seeds.

This does not rule out the observed finite13-of14 recovery, seed-dependent
correspondences, noninvertible maps, or a genuinely new rationality theorem.
It identifies an implication that a propagation proof cannot assume.

## Replay and next decision

The [manifest](../../artifacts/generated-results/elliptic-curves/curve302_class1_bridge_manifest_v1.json)
pins inputs, producers, independent checkers and this note. The complete
portable replay is

~~~sh
python3 research/elliptic-curves/cas/verify_curve302_class1_bridge_bundle.py
~~~

The first point replay repeated large rational point identities at every
finite prime and reached its180-second cap after nine rows. Its refined
checker verifies each exact point once, then evaluates independently
enumerated finite groups. Fourteen checkpointed cases all passed under
60-second caps. The carrier, twist, Veronese and integral-core replays
also passed. Coercion/serialization development failures and the capped
attempt remain locally preserved; they did not change the point roster
or mathematical inputs. Sage10.9/PARI2.17.3 were used.

The previous [descent and matched-fibre results](CURVE302_RELATIVE_DESCENT_AND_MATCHED_FIBRES_2026-09-10.md)
remain current arithmetic inputs. Three support gaps, exact control
half-ideal images, and full class/Selmer upper bounds remain open.

**Next priority:** seek a second independent section on a fixed certified
twist, or a different low-genus carrier construction with compatible
integral geometry. The prescribed4100 overlap and repeated smooth
B-fibre covers are now excluded routes. No rank-jump predictor or
self-propagating rank-jump theorem has yet been proved.
