# Geometry of projective linkage in rank two

Let \(A\) be an abelian variety over a number field with a
Rosati-compatible CM field action by \(E\) of relative rank two.
Assume that \(A\) is not potentially CM and that at least one
conjugate Hodge pair is balanced. Write \(b\) for the number of
balanced pairs. After finite extension define the projective trace
field \(P\subset E^+\) as in the
[projective trace field theorem](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md).
Set \(c=[E^+:P]\). The projective classes on conjugate pairs
are the restriction fibres to \(P\), and every class has a
constant balanced or definite type.

Projective linkage is strictly weaker than isomorphism of connected
monodromy summands, including for absolutely simple varieties with
full geometric endomorphism field \(E\). Nevertheless, the
projective count of failed pairs can be recovered from a smaller
rank-two system after enlarging the coefficient field and making a
CM character twist. The reduction uses absolute Hodge theory. For
one balanced projective class the resulting ordinarity is already
the simple-signature case of Cantoral Farfán, Li, Mantovan, Pries and
Tang, or the dimension-at-most-two case when all pairs are balanced.

## The exact linear-linkage test

Fix an algebraic closure of \(\mathbb Q\), its embedding in
\(\mathbb C\), and a normal closure \(N\) of \(E\) in it.
Write \(d_\sigma\in\{0,1,2\}\) for the multiplicity on the tangent
space at an embedding \(\sigma:E\to N\). Choose the same convention
on slopes as in the split coefficient argument. The *oriented
signature profile* of \(\sigma\) is the function
\(\gamma\mapsto d_{\gamma\sigma}\) on \(\operatorname{Gal}(N/\mathbb Q)\).
This keeps the distinction between multiplicities zero and two.

**Proposition 1.** After extending the base to define the geometric
endomorphisms and obtain connected monodromy, two projectively
linked rank-two summands \(V_\sigma,V_\tau\) are isomorphic precisely
when their oriented signature profiles agree.

**Proof.** Use good absolute degree-one primes of a further
extension containing \(N\), where the rational residue prime \(p\)
splits in \(N\). The relative determinant is \(\delta_v\in E\).
The ratio \(\sigma(\delta_v)/\tau(\delta_v)\) has complex absolute
value one at every embedding of \(N\). It is a unit away from \(p\)
because every Frobenius eigenvalue and its polarized dual multiply
to \(p\) and are integral. At the places over \(p\) the ratio's
valuations are the differences \(d_{\gamma\sigma}-d_{\gamma\tau}\),
with \(\gamma\) varying over \(\operatorname{Gal}(N/\mathbb Q)\).
The determinant valuation of each component is its prescribed
dimension. Equal profiles therefore put every ratio in the fixed
finite group of roots of unity of \(N\), by Kronecker's theorem.
Chebotarev and connectedness make
\(\det(V_\sigma)/\det(V_\tau)\) the constant function one on the
algebraic monodromy group. Indeed, a regular function taking values in
a finite set on its Frobenius elements takes those values on the whole
group, and the connected reductive group is geometrically irreducible.

Projective equality and equal determinants give equality of squared
traces. The trace sum has value four at the identity. The
coordinate ring is a domain, so the two trace functions agree.
Semisimplicity gives an isomorphism. Conversely an isomorphism
implies equal determinants, and their valuations at the split
places give all the profile equalities. \(\square\)

In the non-potentially-CM setting each \(V_\sigma\) is irreducible,
as in the isotypic argument of the linked criterion. Let \(E_{\rm H}\)
be the intermediate field of \(E/\mathbb Q\) corresponding to the
partition by equal oriented profiles. This is a Galois-invariant
partition of a transitive Galois set, so its classes are the
restriction fibres to a unique intermediate field. Intersecting
those fibres with the projective fibres shows that the linear
classes are the fibres of restriction to \(P E_{\rm H}\) inside \(E\).
Faltings's endomorphism comparison identifies this compositum with
the center of \(\operatorname{End}^0(A_{\overline{\mathbb Q}})\).
More explicitly, the centralizer of \(E\) in this endomorphism
algebra is \(E\), since each split \(V_\sigma\) has scalar
centralizer. Thus the center sits in \(E\). Its embeddings index
the isotypic blocks, so its restriction fibres are the linear
classes. If all pairs are balanced, \(E_{\rm H}=\mathbb Q\) and
projective linkage is exactly linear linkage.
[Faltings, Finiteness theorems for abelian varieties over number
fields, Theorems 3–4](https://doi.org/10.1007/978-1-4613-8655-1_2).

The use of \(\operatorname{tr}^2/\det\) is the rank-two form of
the adjoint trace construction. Pink's local trace characterization
and Ribet's squared-trace argument for modular forms are antecedents.
They address fields and images of representations. Proposition 1
adds the oriented split Hodge test for which projective links are
linear in this setting.
[Pink, Compact subgroups of linear algebraic groups, Proposition 0.6](https://people.math.ethz.ch/~pinkri/ftp/LastVersion.pdf),
[Ribet, On \(\ell\)-adic representations attached to modular forms II, proof of Theorem 3.1](https://math.berkeley.edu/~ribet/Articles/rankin.pdf).

## Realization of coefficient twists

We use two elementary operations on polarizable weight-one Hodge
structures of abelian varieties with CM field action.

First, for a CM extension \(E'/E\), scalar extension
\(V\otimes_E E'\) of \(V=H^1(A,\mathbb Q)\) is \(H^1(A^{[E':E]},\mathbb Q)\)
as a rational Hodge structure. The natural \(E'\)-action is
adjoint to CM conjugation for a polarization. One can use the
trace pairing \(\operatorname{Tr}_{E'/\mathbb Q}(u\bar v)\) on \(E'\)
viewed as a weight-zero space and restrict the tensor polarization
to the \(E\)-balanced tensor summand. This summand is cut out
by the diagonal idempotent of \(E\otimes_{\mathbb Q}E\),
fixed by simultaneous conjugation. The relative rank is
unchanged and the old summand is repeated at each extending
embedding.

Second, assign integers \(k_\sigma\in\{-1,0,1\}\) to the embeddings
of \(E'\), with \(k_{\bar\sigma}=-k_\sigma\). There is a rank-one
\(E'\)-linear weight-zero Hodge structure \(T\) with type
\((k_\sigma,-k_\sigma)\) at \(\sigma\). Choose two CM types of
\(E'\) whose indicators differ by \(k_\sigma\). Then
\(T=H^1(Y,\mathbb Q)\otimes_{E'}H^1(Z,\mathbb Q)^\vee\) for
the corresponding CM abelian varieties \(Y,Z\). The choices
of \(k_\sigma\) used below are zero on balanced pairs. On
definite pairs they either vanish or exchange which conjugate
has the two \((1,0)\) entries. Then \(V\otimes_{E'}T\) is polarizable
and still of weight one with types \((1,0),(0,1)\). The tensor
polarization again restricts to the diagonal summand, and
continues to make the \(E'\)-action Rosati compatible.
[Deligne and Milne, Tannakian categories, Theorem 6.20](https://www.jmilne.org/math/xnotes/tc.pdf),
[Deligne, Hodge cycles on abelian varieties, §4 and Main Theorem 2.11](https://jmilne.org/math/Documents/Deligne82.pdf).

For use at arithmetic places, these weight-one twists are realized
over number fields with the prescribed scalar tensor action on
\(\ell\)-adic cohomology. Here is a field-of-definition argument.
Suppose \(A,Y,Z\) are defined over \(\overline{\mathbb Q}\).
CM theory gives such models for \(Y,Z\). The Riemann construction
first produces an abelian variety \(A'\) over \(\mathbb C\).
The Hodge isomorphism between \(H^1(A')\) and the tensor summand is
a rational Hodge tensor on the product of the abelian varieties.
By Deligne's theorem it is absolute Hodge. If
\(s\in\operatorname{Aut}(\mathbb C/\overline{\mathbb Q})\), its
transport identifies \(H^1((A')^s)\) with the same input tensor
summand, so \((A')^s\) is isogenous to \(A'\).
The isogeny class is countable, by the countable choices of
finite subgroups. Endow \(A'\) with a polarization and level
structure and use the fine moduli scheme over
\(\overline{\mathbb Q}\) of that type. The orbit of its moduli
point is countable, since each member of the isogeny class has
countably many polarizations and finitely many level structures
of fixed level. A point on an affine chart of a finite-type scheme
over \(\overline{\mathbb Q}\) with a transcendental coordinate has
uncountable orbit under these field automorphisms. Thus \(A'\)
has a model over \(\overline{\mathbb Q}\), hence over a number field.
The same absolute Hodge tensor, in \(\ell\)-adic cohomology,
is fixed after a finite extension of that number field. For example,
its Galois orbit is contained in the countable set of rational
Betti Hodge tensors under comparison. A countable compact
homogeneous quotient of a profinite group is finite by the
Baire theorem. This gives the claimed \(\ell\)-adic twist and
requires no assertion that the Hodge tensor is algebraic.

At good absolute degree-one primes of a sufficiently large extension
where \(E'\) splits, the CM rank-one components have slope zero or
one equal to their Hodge multiplicity. The twist at any balanced
embedding is a ratio of two such components of equal Hodge type,
so it has slope zero there. It preserves failure or ordinarity on
every balanced component. It only changes the orientation of definite
pairs, which are ordinary. It preserves the projective trace functions
on connected algebraic monodromy because the twist is scalar.

## Reduction to a smaller rank-two system

**Proposition 2.** Suppose that there is at least one definite pair.
There is an absolutely simple abelian variety \(C\) over a number
field with full geometric endomorphism field a CM quadratic
extension \(K/P\). Its relative rank over \(K\) is two, its
Hodge signature has \(b/c\) balanced pairs, and at the common
good absolute degree-one primes of a finite extension the
number of nonordinary balanced pairs of \(A\) is \(c\) times
that number for \(C\).

**Proof.** Temporarily choose any CM quadratic \(K/P\), put
\(E'=E K\), and extend scalars from \(E\) to \(E'\).
The projective partition is still exactly restriction to \(P\).
On balanced fibres leave the twist trivial. On each definite
\(P\)-fibre choose one of the two \(K\)-orientations and exchange
precisely those definite pairs necessary to make the oriented
multiplicity constant on each \(K\)-embedding fibre. The
preceding construction realizes the twisted variety \(A'\).
Its linear linkage classes are exactly the \(K\)-fibres by
Proposition 1. Within each \(K\)-fibre the profiles agree.
Two opposite \(K\)-fibres over one embedding of \(P\) differ
after translating by some \(\gamma\) that takes this embedding
to a definite class. Different \(P\)-fibres are already
projectively distinct. Thus the center of the full geometric
endomorphism algebra of \(A'\) is \(K\). It has one geometric
isotypic factor, so write \(A'\sim C^n\).

We can choose \(K\) so that
\(\operatorname{End}^0(C)=K\). To see
this, examine the adjoint Hodge group. Faltings identifies
the linear intertwiners over each \(K\)-fibre with algebraic
endomorphisms, which also intertwine the Hodge group. Polarization
identifies the projective actions across conjugates. Thus the
linked \(P\)-fibres agree projectively for the Hodge group as well.
Conversely Hodge linkage implies \(\ell\)-adic linkage, since
connected \(\ell\)-adic monodromy is contained in the Mumford–Tate
group by Deligne's theorem. Every derived projection of the Hodge
group is \(\mathrm{SL}_2\), since it is semisimple and contains
the corresponding nontrivial \(\ell\)-adic derived projection.
As in the linked criterion, the adjoint simple factors are
therefore indexed by the embeddings of \(P\) as a Galois set.
Descent from the factor at the structural embedding of \(P\)
gives a \(P\)-form of \(\mathrm{PGL}_2\), so the adjoint Hodge
group is \(\operatorname{Res}_{P/\mathbb Q}
\mathrm{PGL}_1(D_0)\) for a quaternion algebra \(D_0/P\), with
the matrix algebra allowed. Here the quaternion algebra is
the degree-two central simple algebra obtained by descending
the matrix algebra under the same inner descent datum.
The adjoint Hodge group, with this indexing, is unchanged
by coefficient extension and the weight-zero CM twist. For
the latter claim, take the Hodge group of the direct sum of
the original Hodge structure and the CM inputs. Its projection
to the original Hodge group is surjective, its derived group
acts trivially on the CM inputs, and on the weight-one twist
its derived action is the original action with possible
repeated coefficient embeddings. Hence the adjoint derived
group and its geometric factor labels are preserved.

Now choose \(K/P\) totally imaginary quadratic and nonsplit
at each finite place where \(D_0\) ramifies. Weak approximation on
a nonsquare or negative defining element and the quaternion
local–global splitting criterion show that \(K\) splits \(D_0\).
Make the preceding scalar extension and flattening with this
\(K\). Write \(D=\operatorname{End}^0(C)\), with center \(K\) and
Schur index \(t\). As a \(K\)-space the canonical \(K\)-component
of \(H^1(C)\) restricts to \(t\) copies of the standard
representation of the simple factor at the structural
embedding of \(P\), after scalar extension to an algebraic
closure. That factor in the simply connected cover becomes
\(\mathrm{SL}_2\) over \(K\), because \(D_0\) splits there.
The representation over \(K\) is therefore
\(K^2\otimes_K K^t\), by extending scalars in the intertwiner
space from the split standard representation. Its commutant
is \(M_t(K)\). The algebraic endomorphism algebra \(D\) embeds
in this commutant and both have \(K\)-dimension \(t^2\).
As \(D\) is division, \(t=1\). This gives full endomorphism
field \(K\), and relative rank two, for \(C\).

The \(K\)-signature of \(C\) has one balanced pair per balanced
projective class. At split good degree-one primes the balanced
coefficient twist has slope zero and every \(E'\)-embedding over
a \(K\)-embedding has the rank-two slopes of the corresponding
component of \(C\). Scalar extension just repeated each old
component, so one bad pair of \(C\) accounts for exactly \(c\)
bad pairs of \(A\). \(\square\)
[Deligne, Hodge cycles on abelian varieties, §§2 and 6](https://jmilne.org/math/Documents/Deligne82.pdf),
[Voight, Quaternion Algebras, Proposition 14.6.7](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_14).

If \(b=c\), \(C\) has simple signature. Cantoral Farfán et al.,
Theorem 4.18, applies to \(C\) with full endomorphism field \(K\).
If \(b=2c\), the earlier two-balanced-pair theorem applies to
\(C\). For general \(b\) the same rank-two bound on \(C\)
recovers the projective bound for \(A\).
[Cantoral Farfán et al., Assumptions 2.1 and 2.5 and Theorem 4.18](https://arxiv.org/html/2508.11174v1#S4.SS5),
[Rank-two ordinary-reduction theorem](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md#ordinary-primes-for-rank-two-multiplication).

If all pairs are balanced, the linear test makes this a direct
isogeny reduction. The endomorphism center is then \(P\),
and \(A\sim C^n\) with \(nt=2c\) for the Schur index \(t\) of
\(\operatorname{End}^0(C)\). The dimension of \(H^1(C)\) over
\(P\) is \(2t\), hence \(t\mid2\). If \(t=1\), \(C^2\) admits
the Rosati-compatible rank-two action of \(P(i)\) by the
standard rotation on the two copies of \(C\) and the product
polarization. If \(t=2\), take a nonzero Rosati-skew element
of the quaternion algebra \(\operatorname{End}^0(C)\)
over \(P\). It has reduced trace zero, its square is
central, and positivity makes that square totally negative.
It generates a Rosati-compatible CM quadratic extension
of \(P\) acting on \(C\) in relative rank two.
In both cases the resulting rank-two system has
\([P:\mathbb Q]=b/c\) balanced pairs. The two conjugate
representations over each embedding of \(P\) occur in the same
linear block and their Hodge multiplicities are equal and
complementary, hence one on each side. The split Frobenius slopes
in each such pair are the slopes of the corresponding \(P\)-block,
so the same multiplicity \(c\) transfers the count of failures.
For \(b=c\), \(C\) has dimension \(t\le2\) and the known
potential ordinary-reduction theorem applies. For \(b=2c\),
the two-balanced-pair theorem applies to the smaller rank-two
system.
[Sawin, Ordinary primes for abelian surfaces, Corollary 2.2](https://www.numdam.org/item/10.1016/j.crma.2016.01.025.pdf).

## A projective-only example with full CM endomorphism field

Take \(K=\mathbb Q(\zeta_5)\) and a non-CM abelian fourfold \(B\)
in the family of Jacobians of \(y^5=x(x-1)(x-u)\).
Such a fourfold exists over a number field. The family has positive
dimension, with generic Hodge signature one mixed pair and one
definite pair; Noot's specialization theorem in its connected
Shimura variety supplies a number-field realization of the
non-CM Mumford–Tate group. Its full endomorphism field is \(K\), as
explained in the example in §5.2 of Cantoral Farfán et al.
Let \(M=\mathbb Q(\zeta_7+\zeta_7^{-1})\), a cyclic real cubic,
and put \(E=K M\). The two fields are linearly disjoint, so
\(\operatorname{Gal}(E/\mathbb Q)\cong C_4\times C_3\).

Extend scalars to \(E\), giving \(B^3\) up to isogeny.
The two \(K^+\)-classes both have three conjugate pairs. Exactly
one class is balanced, and the other class is entirely definite.
Make a weight-zero \(E\)-CM twist as above that exchanges the
definite types on exactly one pair and leaves all other types
unchanged. The resulting variety \(X\) has dimension twelve
and three balanced pairs. Its projective trace field is \(K^+\).
Indeed, it is the same after the twist, with three copies of
each original pair of \(B\), and the distinct original mixed
and definite pair could not be projectively linked by parity.

The oriented signature function on \(C_4\times C_3\) has
trivial stabilizer under translation. An odd translation in
\(C_4\) would exchange balanced and definite classes. The
nontrivial even translation in \(C_4\) interchanges zero and two
on the definite class, where their multiplicities differ.
With the \(C_4\) coordinate fixed, an element of \(C_3\)
must preserve the unique flipped pair, so is trivial.
Consequently distinct \(E\)-embeddings have different
oriented profiles and no linked summands are linearly
isomorphic. Proposition 1 and Faltings's comparison give
\(\operatorname{End}^0(X_{\overline{\mathbb Q}})=E\), so \(X\)
is absolutely simple. This gives projective-only linkage even
with full geometric endomorphism field \(E\). Its ordinary
reduction can also be proved by transporting that of \(B\)
through the coefficient extension and CM twist.
[Cantoral Farfán et al., §5.2](https://arxiv.org/html/2508.11174v1#S5.SS2),
[Noot, Compos. Math. 97 (1995), 161–171, Theorem 1.7 and proof](https://www.numdam.org/item/CM_1995__97_1-2_161_0.pdf#page=6).

## Comparison with positive-density criteria

Fité's Galois coefficient-field theorem gives the weak trace-norm
threshold in relative rank two, and his real-cubic threefold
theorem imposes a noncontainment condition between the cubic
field and a field of definition of its action. Wang and Zhang
prove positive density for specified fields of definition and
degrees, and in Theorem 1.3 impose an auxiliary splitting
condition and exclude Albert type IV. These form separate routes
to positive density. The projective trace field theorem gives
density one over the action field for the linked non-potential-CM
cases, including nonsplit primes. For the potential conclusion,
its geometric reduction, with
absolute Hodge CM twists allowed, accounts for its one-class
ordinary case using the earlier simple-signature theorem, and
for its two-class ordinary case using the two-balanced-pair
theorem. Neither theorem supplies ordinary reduction in the
unlinked three-class case. The sources examined here do not
establish priority for the projective reduction or the exact
oriented-profile criterion.
[Fité, Theorems 2–3](https://arxiv.org/html/2005.10185v3#S1),
[Wang–Zhang, Theorems 1.1–1.3](https://arxiv.org/html/2503.21111v2#S1).
