# A stable hypersurface obstruction over the Keller base

14 November 2026. The target-dependent matrix hypothesis in
[HC4MYHI3](RECONSTRUCTION_EMBEDDING_OBSTRUCTION.md#4-a-target-dependent-affine-dual-frame-cannot-repair-a-nonautomorphic-base-pair-hc4myhi3)
can be removed entirely. Neither source-dependent matrices nor nonlinear
dual dependence evade the obstruction **when the first two physical
coordinates are still the plane Keller outputs**.

Let K=(K1,K2):A2_z->A2_t be a complex plane Keller map, r>=0,
and let H=(H1,...,H_(r+1)) be arbitrary polynomials in z,w, with r
independent w variables. Then

\[
 j(z,w)=(K(z),H(z,w)):A^{2+r}\longrightarrow A^{3+r}
 \quad\hbox{closed embedding}\quad\Longrightarrow\quad K\text{ is invertible}. \tag{1}
\]

There is no restriction on degrees, affine dependence, the shape of an
individual fibre embedding or the availability of an ambient inverse.
The inverse required by the hypothesis is precisely that j is a closed
embedding, equivalently that its coordinate functions generate the
whole source polynomial algebra. A polynomial substitution alone is
insufficient. This note supplies a full proof; no novelty priority is
asserted.

## 1. Disjoint conjugate components form one pencil (HC4MYHJ0)

We first isolate two algebraic steps, including the denominator issue.

**Disjoint equal-degree components.** Let Omega be algebraically
closed and let e1,...,ed in Omega[D1,...,D_(r+1)] be distinct
irreducible polynomials. Suppose their zero sets are pairwise disjoint,
each quotient Omega[D]/(e_i) has only constant units, and all e_i have
the same total degree. Normalize their common leading monomial to
coefficient one, when they are conjugate over a field F.

For i!=j, disjointness says that e_j is a unit modulo (e_i).
The constant-unit hypothesis gives

\[
                         e_j-c_{ij}=q_{ij}e_i,
                         \qquad c_{ij}\in\Omega^\times.       \tag{2}
\]

The degree equality makes q_ij a scalar. Common leading coefficient
one makes that scalar one. Consequently all nonconstant coefficients
of the e_i agree. If the components form one conjugate orbit over F,
these coefficients are fixed by every F-automorphism of an algebraic
closure, so they belong to F in characteristic zero. Thus

\[
                            e_i=h_0(D)-a_i,
                 \qquad h_0\in F[D],\quad a_i\in\Omega.      \tag{3}
\]

The a_i are distinct. Galois invariance gives a monic polynomial
R0(T)=product_i(T-a_i) in F[T]. Any defining polynomial over F
having exactly these reduced components is a nonzero F-scalar times
R0(h0). No claim is made that h0 is a coordinate polynomial: an
embedding rectification theorem is not used.

**Primitive composition has no hidden denominator.** Let A be a
UFD, F its fraction field, and h in A[D] have constant coefficient
zero and gcd of its nonconstant coefficients equal to one. Then

\[
                              F[h]\cap A[D]=A[h].             \tag{4}
\]

For any irreducible pi in A, localize at the prime (pi). At least
one nonconstant coefficient of h remains nonzero modulo pi, so its
reduction is transcendental over Frac(A/(pi)). Suppose coefficients
of a polynomial Q(T) in F[T] have negative minimum pi-valuation.
Multiply Q by the compensating power of pi so its coefficients lie
in A_(pi), with at least one nonzero modulo pi. The resulting
nonzero polynomial cannot vanish after substituting the nonconstant
reduction of h. Thus Q(h) has the same negative valuation. It cannot
belong to A[D]. Since this holds for every irreducible pi, all
coefficients of Q belong to A. This proves (4), in all degrees.

Starting with h0 in (3), subtract its constant coefficient, clear
denominators and divide the gcd of the remaining coefficients. This
gives h satisfying (4), with F[h]=F[h0]. In particular, if E in
A[D] has the components (3), then

\[
                         E=R(t,h(D)),\qquad R\in A[T],
                         \quad\deg_T R=d.                    \tag{5}
\]

The factorization is over the fraction field first; (4) is exactly
what upgrades it to a polynomial identity over the full target.

## 2. Arbitrary stable hypersurface embeddings force an inverse (HC4MYHJ1)

Put A=C[t1,t2], B0=C[z1,z2], F=C(t1,t2), and L=C(z1,z2),
using t=K(z) for the A-algebra structure. The Keller condition makes
K1,K2 algebraically independent and L/F finite separable. Moreover

\[
                         B0\otimes_A F=L.                    \tag{6}
\]

Indeed this localization is a finitely generated F-algebra, algebraic
over F, hence a finite-dimensional domain over F and therefore a
field with fraction field L.

Assume j in (1) is a closed embedding. Its image has prime
codimension-one ideal (E) in the UFD A[D1,...,D_(r+1)], and

\[
                       A[D]/(E)\simeq B0[w].                 \tag{7}
\]

Its image is smooth, being affine space. Its generic target fibre is
Spec L[w]. After extending F to an algebraic closure Omega this is
the disjoint union of d=[L:F] copies of A^r_Omega. Hence E factors
over Omega into d distinct irreducible polynomials e_i, with
pairwise comaximal ideals. Each quotient is Omega[w], whose units
are Omega nonzero. Galois acts transitively on these factors, so
they have the same total degree and the same leading monomial.

Apply Section 1. Normalize the factors, obtain their common
nonconstant part over F, and make it primitive over A. Then (5)
holds globally. If R factored nontrivially in A[T], substituting the
nonconstant polynomial h would factor E nontrivially in A[D].
Thus R is irreducible. No total-degree bound is involved.

By construction h evaluated on any geometric source fibre is the
constant a_i labeling its component. Therefore

\[
                  h(K(z),H(z,w))=a(z),\qquad a\in B0.        \tag{8}
\]

This is a polynomial identity globally: the left side is polynomial
and its positive-w coefficients vanish generically. Equivalently
a(z)=h(K(z),H(z,0)). From E=R(t,h) we obtain

\[
 E_{D_j}(j)=R_T(K,a)\,h_{D_j}(j)\quad\hbox{for all }j.      \tag{9}
\]

The scalar R_T(K,a) has no source zero. If it vanished at some
(z,w), all E_D derivatives would vanish there. Differentiating
E(K(z),H(z,w))=0 in z would then give

\[
                                DK(z)^T E_t(j)=0.
\]

Since DK is invertible, E_t would also vanish, contradicting the
smoothness of the image hypersurface. Thus R_T(K,a) is a unit in
the whole source polynomial ring, so it is a constant c in C nonzero.

The irreducible relation R is the minimal relation of a over F:
R(K,a)=0, and irreducibility and Gauss's lemma rule out a smaller
relation. The polynomial R_T-c also vanishes at (K,a), but has
T degree strictly smaller than R. Consequently

\[
                        R_T=c,\qquad R=cT+b(t),\qquad d=1.   \tag{10}
\]

Here characteristic zero is essential. Thus L=F: the plane Keller
map is birational. The retained birational Keller case in
[HC4MYOR2](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length)
then gives its polynomial inverse. This proves (1).

The same pencil and primitive-composition argument forces generic
degree one for a Keller base of any dimension. The registered HC4
conclusion uses the retained plane birational endpoint only. For
r=0, the argument recovers the classical Formanek single-generator
criterion reproduced in HC4MYHI0, rather than asserting that
classical result is new.

## 3. Exact stable generator minimum and the graph consequence

For a hypothetical nonautomorphic plane Keller map K and every r>=0,
the polynomial algebra C[z,w] needs **exactly r+2 additional algebra
generators over C[K]**. The functions z1,z2,w1,...,wr supply r+2.
If r+1 or fewer sufficed, pad by zero if necessary; their pullbacks
would give (1), a contradiction. This is a relative algebra-generator
minimum over the stipulated Keller target, not an absolute embedding
dimension. It is stable under adding arbitrary polynomial variables.

In the HC4 case r=2, no nonautomorphic plane Keller pair admits any
closed embedding

\[
                         (z,w)\longmapsto(K(z),H(z,w))
                                   \quad\hbox{into }A^5,    \tag{11}
\]

even when H is fully nonlinear in all four source variables. This
closes the source-dependent affine matrices left outside HC4MYHI3,
and every nonlinear extension over the same physical Keller base.
Changing only the dual dependence cannot repair that encoding.

Under Tang's hypothetical nonautomorphic HC4 form
f=D(z)+w^T K(z), the two w-directional gradient components are K.
Therefore no actual closed embedding of the source in A5 can have
these two components as restrictions of an ambient coordinate pair
(a pair completing to an actual ambient polynomial automorphism).
Compose with that ambient automorphism and apply (11). This needs
only the cotangent pair to be ambient coordinates, rather than all
four gradient components as in HC4MYHI1.

These conclusions also apply after any actual source coordinate
change or ambient automorphism that preserves the stipulated output
base interpretation. A general graph in the original Meng--Yang
coordinates need not have its physical base equal to K(z), nor must
the intrinsic cotangent gradient pair be an ambient coordinate pair.
The existing unit-ideal and singleton-zero constructions use different
physical-base data; neither is excluded by (11).

## Hypotheses, sources and controls

The proof uses closed embedding, polynomial source units, separability,
the constant-unit property of each geometric generic component, and
the primitive-composition step (4). It does not use an inverse for
the Keller map before the birational conclusion, an embedding
rectification theorem for A2 in A3, or Hessian invariance under a
nonlinear coordinate change.

The basic Galois normal-closure and fixed-field facts are
[Stacks, Lemma 9.21.5 and Theorem 9.21.7](https://stacks.math.columbia.edu/tag/09DU).
Basic separable fibre and base-change facts are
[Stacks, Lemmas 29.37.4 and 29.37.7](https://stacks.math.columbia.edu/tag/02GH).
The [source record](sources/stable-hypersurface-import.json) retains
the inspected primary bytes. The pencil and denominator lemmas are
proved here, and no new classification theorem is imported.

The [exact checker](verify_stable_hypersurface.py) and
[receipt](receipt-stable-hypersurface.json) keep finite algebraic
controls separate from this written arbitrary-degree argument.
Negative controls show why a ramified base, a source with nonconstant
units, nonconstant component units or an unnormalized denominator
cannot be silently admitted. Positive controls include a genuinely
nonlinear fibre embedding and exact source recovery.

The preceding turn closed exact A/B encodings and target-dependent
affine matrices. The present result removes **every condition on H**
over a fixed Keller output base. What remains is different physical-base
data or the lack of an ambient coordinate extension of the intrinsic
cotangent pair, together with the universal SL2 row equation and the
distinct unit-ideal and singleton-zero plane endpoints. General graph
embeddings and HC4/JC2 remain unresolved. No formal verification,
independent replay, external review or novelty priority is claimed.
