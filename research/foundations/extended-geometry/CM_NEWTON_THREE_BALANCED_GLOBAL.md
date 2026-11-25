# Three balanced pairs in relative rank two

Let \(A\) be an abelian variety over a number field with a geometric CM
field action by \(E\), compatible with a polarization in the sense that
Rosati acts on \(E\) by complex conjugation. Write \([E:\mathbb Q]=2e\) and
assume \(\dim A=2e\). Label a conjugate pair of embeddings of \(E\) balanced
if its Hodge multiplicities are \((1,1)\). This note assumes there are
exactly three balanced pairs. It gives density-one ordinarity after finite
extension when two different conjugate pairs have linked projective
monodromy, and gives a separate positive-density sufficient condition at
primes with a suitable nonsplit coefficient field. The remaining split
case includes real-multiplication threefolds over arbitrary number fields.

## Inputs and notation

The [rank-two trace argument](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md#ordinary-primes-for-rank-two-multiplication)
gives the relative trace \(a_v\in O_E\) of Frobenius at a good prime of
absolute degree one, and
\[
 |N_{E/\mathbb Q}(a_v)|\le 2^{2e}p^e,
 \qquad v_p(N_{E/\mathbb Q}(a_v))\ge e-3+2r                              \tag{1}
\]
at split primes, where \(p\) is the residue characteristic and \(r\) is the
number of nonordinary balanced pairs. The first inequality is the Weil
bound. Each definite pair contributes at least one factor of \(p\), and
each nonordinary balanced pair contributes at least two. On a density-one
set after finite extension, \(r\le1\). For this conclusion it suffices to
discard trace-zero primes by algebraic Chebotarev and compare the two
bounds for \(r\ge2\).

We also use the relative determinant \(d_v\in O_E\). The trace-pairing
construction applied to Frobenius and its square gives both relative
traces in \(E\), and the rank-two identity
\(2d_v=a_v^2-\operatorname{Tr}_E(\pi_v^2)\) gives compatibility of the
determinant between crystalline and étale realizations. Polarization gives
\(d_v\bar d_v=p^2\) and \(\bar a_v= p a_v/d_v\).

## The projective trace criterion

**Proposition 1.** Pass to a finite extension \(L\) defining the action,
polarization and geometric endomorphisms, containing a normal closure of
\(E\), and with connected \(\ell\)-adic algebraic monodromy group \(G\).
For each embedding \(\sigma:E\hookrightarrow\overline{\mathbb Q}_\ell\),
let \(V_\sigma\) be its rank-two summand of \(H^1_{\mathrm{et}}\). Define
\[
 B_\sigma(g)=\frac{\operatorname{tr}(g\mid V_\sigma)^2}
                   {\det(g\mid V_\sigma)}.                              \tag{2}
\]
If \(B_\sigma=B_\tau\) on \(G\) for embeddings in different conjugate pairs,
then \(A_L\) has density-one ordinary reduction.

**Proof.** Let \(\sigma\sim\tau\) mean equality of the functions (2), and
identify their indices with embeddings into \(\overline{\mathbb Q}\) through
one embedding into \(\overline{\mathbb Q}_\ell\). The relation is equivalent
to \(\sigma(a_v^2/d_v)=\tau(a_v^2/d_v)\) at every good absolute degree-one
prime \(v\) of \(L\), because these Frobenius elements are dense for
regular class functions on \(G\) by Chebotarev. Polarization gives
\(\bar a_v^2/\bar d_v=a_v^2/d_v\). Thus the projective trace field
\[
 E_0=\mathbb Q\bigl(a_v^2/d_v : v\text{ good and of absolute degree one}\bigr)
       \subseteq E^+
\]
describes the relation. Its classes of conjugate pairs are the fibres of
restriction of the embeddings of \(E^+\) to \(E_0\), so they all have
cardinal \(c=[E^+:E_0]\).

The same parity and valuation argument applies with any positive number
\(b\) of balanced pairs. It gives \(c\mid b\). The number of nonordinary
balanced classes is at most \(\lfloor(b/c-1)/2\rfloor\), since each
failed class accounts for \(c\) failed pairs. Thus any linkage between
distinct pairs forces density-one ordinarity after finite extension
when \(b\) is prime. This numerical consequence does not classify when
the projective linkage occurs.

Take a large absolute degree-one good prime of \(L\) with \(a_v\ne0\).
Such primes form a density-one set, since the trace functions are nonzero
at the identity of connected \(G\). Identify a fixed algebraic closure
with an embedded subfield of \(\overline{\mathbb Q}_p\) for the computation.
The determinants have valuations zero and two in a definite pair, and
one in a balanced pair. The traces have integral valuations because the
coefficient field splits completely. Thus the valuation of (2) is odd at
every balanced embedding and even at every definite embedding. Each
equivalence class is entirely balanced or entirely definite. Since there
are exactly three balanced pairs, \(c\mid3\). If any link exists between
different pairs, \(c=3\), and all balanced pairs are in one class. In a
balanced component the trace is a unit precisely when the Newton slopes
are \((0,1)\). Equality of (2) makes trace divisibility simultaneous in
the three balanced pairs. Any failure would therefore have \(r=3\),
which contradicts (1) at all sufficiently large primes with \(a_v\ne0\).
\(\square\)

In the non-CM cases below every \(V_\sigma\) is an irreducible \(G\)-module.
Then two different conjugate pairs are linked exactly when their derived
actions use the same simple factor of \(G^{\mathrm{der}}\) up to central
isogeny. Each projection of the derived group is \(\mathrm{SL}_2\) because
it is nontrivial and semisimple on a two-dimensional space. Its map is
nontrivial on exactly one simple ideal of the derived Lie algebra.
When two projections use the same ideal, their restrictions differ only
by an inner change of basis of the standard \(\mathrm{SL}_2\) representation.
On the connected center they act by scalars, which cancel in (2). When
they use different ideals, (2) varies on a torus of one factor and is
constant on that torus for the other. Consequently the unlinked case has
the full product of \(e\) independent \(\mathrm{SL}_2\) factors on the
conjugate pairs, up to central isogeny. This uses Faltings's semisimplicity
and endomorphism comparison and the elementary structure of a subgroup
that projects onto factors of type \(A_1\), with no Mumford–Tate conjecture.

## Reduction through geometric endomorphisms

**Proposition 2.** Outside the density-one cases in Proposition 1 and the
potential-CM case, there are only the following geometric endomorphism
possibilities.

1. \(A\) is absolutely simple with \(\operatorname{End}^0(A_{\bar L})=E\)
   and \(e\ge4\).
2. \(A\) is geometrically isogenous to the square of an absolutely simple
   threefold with geometric endomorphism field a totally real cubic
   field.
3. \(A\) is an absolutely simple sixfold with geometric endomorphism
   algebra a totally indefinite quaternion division algebra over a
   totally real cubic field.

In all three possibilities that remain, the derived projections onto
different conjugate pairs are independent as above. These are the
remaining cases for the split-prime density-one question.

**Proof.** A projection onto any geometric isotypic factor retains a
faithful field action. The component's homology dimension is a positive
multiple of \([E:\mathbb Q]=2e\), and the total homology has dimension
\(4e\). If there are two distinct isotypic factors, both are of rank one
over \(E\), hence potential CM. The torus argument in the
[rank-two proof](RANK_TWO_CM_ORDINARITY_AND_QUINTIC_FACE_EXCLUSION.md#ordinary-primes-for-rank-two-multiplication)
gives density-one ordinary reduction after finite extension for a
potential-CM product.

Otherwise \(A\sim C^m\) geometrically, with \(C\) absolutely simple. Put
\(D=\operatorname{End}^0(C)\), let \(F\) be its center and let \(t\) be its
Schur index over \(F\). Faltings identifies the connected monodromy
commutant with \(M_m(D)\) after scalar extension. Over
\(\overline{\mathbb Q}_\ell\) its blocks are indexed by embeddings of \(F\), and
each has size \(mt\). The corresponding monodromy-isotypic summand is a
sum of \(mt\) copies of an irreducible representation, of dimension
\(4e/([F:\mathbb Q]mt)\). The dimension of the entire summand is uniform
because rational homology is an \(F\)-vector space. Each \(V_\sigma\) is a
sum of irreducibles of this common dimension, so that dimension is one or
two. If it is one, the monodromy is a torus and the same rank-one splitting
argument gives potential CM. Otherwise each \(V_\sigma\) is irreducible,
and \([F:\mathbb Q]mt=2e\).

The centralizer of \(E\) in \(M_m(D)\) has rational dimension \(2e\),
since after scalar extension it has one scalar block for every
\(V_\sigma\). Hence it equals \(E\) and contains the center \(F\). In
particular \([E:F]=mt\). Embeddings over the same embedding of \(F\) give
isomorphic monodromy representations, so are linked by (2). In the
unlinked case a fibre over \(F\) therefore has size at most two, and when
it has size two the embeddings are conjugates. If \([E:F]=1\), then
\(m=t=1\) and the full endomorphism algebra is \(E\). If \([E:F]=2\),
then \(F=E^+\) is totally real. The isomorphic representations for the
two embeddings over a place of \(F\) give equal determinant valuations
at split good primes. Their sum is two, so every pair is balanced.
Thus \(e=3\). Now \(mt=2\) gives either \((m,t)=(2,1)\), the
real-cubic square in case 2, or \((m,t)=(1,2)\), the quaternion case.

In the quaternion case the Rosati involution is symplectic if the algebra
is ramified at any real place. Indeed on the real Hamilton quaternions the
positive involution of the first kind is the canonical involution. An
orthogonal involution would be canonical conjugation composed with inner
conjugation by a pure quaternion. It would fix a nonzero pure quaternion
anticommuting with that quaternion, whose square has negative reduced
trace, contrary to positivity. Involution type is preserved by scalar
extension. Over an algebraic
closure of \(\mathbb Q_\ell\) each \(F\)-component has the form
\(U\otimes W\) with \(\dim U=\dim W=2\) and the quaternion matrix algebra
acting on \(U\). Compatibility with Rosati makes the polarization on this
component the tensor of the alternating form on \(U\) and a symmetric
nondegenerate form on \(W\). The connected derived group acts on \(W\) through
\(\mathrm{SO}_2\), so acts trivially. This would give toric monodromy,
contrary to the assumed two-dimensional irreducibles. The quaternion
algebra must split at every real place.

Finally, full endomorphism field \(E\) is impossible when all its pairs
are balanced in relative rank two. On the rational weight-one Hodge
structure, the \(E\)-determinant has type \((1,1)\) at every embedding.
The Hodge group acts trivially on these determinant lines. Its
representations on the two conjugate summands are dual by polarization,
and become isomorphic because a rank-two determinant-one representation
is self-dual. This forces off-diagonal operators in the complex Hodge
commutant. The rational Hodge commutant is the rational algebraic
endomorphism algebra for an abelian variety. Indeed, a rational Hodge
endomorphism becomes integral after multiplication by an integer, then
induces an analytic homomorphism of complex tori, which is algebraic.
Formation of the commutant commutes with scalar extension, so it cannot
have been \(E\). This eliminates \(e=3\) from case 1. \(\square\)

## A positive-density nonsplit criterion

**Proposition 3.** Suppose \(A\) is not potential CM. Let \(K\) be a
number field defining \(A\), the \(E\)-action and the polarization. Let
\(N\) be a normal closure of \(E\) in \(\overline{\mathbb Q}\). If an element
of \(\operatorname{Gal}(KN/K)\) cycles the three balanced conjugate pairs,
then ordinary primes of \(K\) contain a set of positive natural density.

**Proof.** The balanced set is Galois-stable over \(K\), since the
multiplicities of the action on the tangent space are defined over \(K\).
By finite Chebotarev, the chosen conjugacy class gives a positive-density
set of unramified primes of \(K\), of which all but a density-zero subset
have absolute degree one. Let \(p\) be the residue characteristic of one
of them and choose an embedding into \(\overline{\mathbb Q}_p\) over the
degree-one place. For a coefficient place \(\lambda\mid p\) of \(E\), the
crystalline \(E_\lambda\)-component has rank two and is linear over
\(E_\lambda\) because the residue field of the prime of \(K\) is
\(\mathbb F_p\). The Hodge multiplicity is constant over its embeddings
into \(\overline{\mathbb Q}_p\). Its determinant valuation is that common
multiplicity: this follows by dividing the dimension of the
\(p\)-divisible component, equivalently the \(\mathbb Q_p\)-determinant
valuation, by \([E_\lambda:\mathbb Q_p]\). The Newton slopes are also
constant over these embeddings. The arguments for the lower bound in (1)
thus apply to each embedding, counting the residue degree of \(\lambda\)
in the norm if the prime is nonsplit. In particular a nonordinary
balanced pair has both traces of positive valuation, which is an integer
because the coefficient place is unramified. Definite pairs are ordinary.

The set of nonordinary balanced pairs is invariant under the Frobenius
permutation of coefficient embeddings. If one such pair fails, all
three fail, since the permutation cycles them. Bound (1) then forces
\(a_v=0\) at all sufficiently large nonordinary primes of the chosen
Chebotarev class.

The trace-zero primes have density zero even over this \(K\), where
connectedness has not been assumed. In fact the normalized trace-norm
class function \(B=\chi^{-e}\prod_\sigma\operatorname{tr}_{V_\sigma}\)
is nonconstant on every component of the monodromy over \(K\).
The preceding isotypic argument shows that each rank-two summand in the
non-CM case restricts to an irreducible representation of \(G^0\), so each
derived projection is \(\mathrm{SL}_2\). Choose a cocharacter \(\lambda\)
of a maximal torus of \((G^0)^{\rm der}\) whose two weights on every
component are \(n_\sigma,-n_\sigma\) with \(n_\sigma\ne0\).
For a component \(gG^0\), choose \(h\in G^0\) so all prospective leading
coefficients of \(\operatorname{tr}((gh)_\sigma\lambda_\sigma(z))\) in
\(z\) are nonzero. Each coefficient is a diagonal entry of \((gh)_\sigma\)
in the weight basis, and its nonvanishing defines a nonempty open set:
the projected \(\mathrm{SL}_2\) spans the full matrix algebra and
\(g_\sigma\) is invertible. The finite intersection of these opens is
nonempty in the connected group \(G^0\). Every \(g\) preserves each
summand because the action is defined over \(K\). On \(gh\lambda(z)\) the
normalized product has nonzero leading coefficient of degree
\(\sum_\sigma|n_\sigma|>0\), because \(\chi(\lambda(z))=1\).
Hence all fixed level sets of \(B\), including zero, have density-zero
Frobenius by algebraic Chebotarev. This proves the proposition.
\(\square\)

This condition gives ordinary primes in the selected class over \(K\).
Adjoining \(N\) removes the chosen nonsplit class from the absolute
degree-one primes. The use of nontrivial
coefficient-field residue degrees is closely related to Fité's
[Theorem 3 and Section 6](https://arxiv.org/pdf/2005.10185v3) for cubic real
multiplication, and Wang–Zhang's
[Theorems 1.3 and 1.7–1.8](https://arxiv.org/html/2503.21111v2#S1).
Wang–Zhang's abelian-variety theorem requires an absolutely simple variety
of the first kind and the specified inert or equal-degree two-prime
splitting condition. Proposition 3 instead uses the orbit of the
balanced pairs under a Rosati-compatible CM action.

The same calculation gives quantitative bounds over the original field
of definition. Let \(H\subseteq S_3\) be the image of
\(\operatorname{Gal}(KN/K)\) on the balanced pairs. At almost every
absolute degree-one prime the number of failing pairs is at most one
by (1). Moreover at such a nonordinary prime the trace-norm valuation
is exactly \(e-1\). A larger valuation makes \(p^{-e}N_{E/\mathbb Q}(a_v)\)
one of finitely many bounded integers, whose level sets of \(B\) have
density zero. Thus the two bad traces have valuation one, all other
balanced traces are units, and the slope-zero traces of definite pairs
are units. The unique failing pair is
a fixed point of the Frobenius permutation in \(H\). Every three-cycle
class consists of ordinary primes up to a density-zero subset. Thus the
lower natural density of ordinary primes is at least \(2/3\) when
\(H=C_3\), and at least \(1/3\) when \(H=S_3\). Transitivity of \(H\) is
equivalent to the absence of a balanced embedding of \(E^+\) with image
in \(K\): an intransitive group on three letters has an individually
fixed letter, and an unordered conjugate pair is exactly an embedding
of \(E^+\). In particular, if \(H\) is generated by a transposition, any
failure on its nontrivial Frobenius class is confined to the same one
balanced pair. The identity class leaves all three individual choices
available.

For a real-cubic RM threefold, the published Fité theorem applies when
the multiplication is defined over a field containing no copy of the
cubic coefficient field. Section 6 uses an absolute degree-one prime
inert in the cubic extension after an auxiliary quadratic base change.
Wang–Zhang's condition (P1) in Theorem 1.3 gives the same inert-prime
route there. Their condition (P2), two distinct places of equal inertia
degree, cannot occur at an unramified cubic coefficient prime. The
formal estimate under (P1) uses three places' worth of (p)-valuation
in the coefficient norm. At a split prime one failing place supplies
one, while the Weil bound for the real-cubic trace norm has power
\(p^{3/2}\). The published statements and proofs are in
[Fité, DOI 10.5565/PUBLMAT6812402, Theorem 3 and §6](https://doi.org/10.5565/PUBLMAT6812402)
and [Wang–Zhang, DOI 10.1112/mtk.70081, Theorems 1.3, 1.7–1.8 and §6.1](https://doi.org/10.1112/mtk.70081).

## Exact obstructions to extending the split-prime proof

The global question includes arbitrary real-multiplication threefolds.
Let \(B\) be a polarized threefold with action by a totally real cubic
field \(F\), and assume Rosati fixes \(F\). In the unresolved simple
geometric endomorphism-field case this assumption holds automatically.
Take \(A=B^2\) with its product polarization. The diagonal \(F\)-action
and \(\left(\begin{smallmatrix}0&-1\\1&0\end{smallmatrix}\right)\) give
a Rosati-compatible action by \(E=F(i)\). Each \(F\)-embedding occurs
once on the Lie algebra of \(B\), because its homology has rank two
over \(F\) and real embeddings have equal conjugate Hodge multiplicities.
Thus \(E\) has three balanced pairs on \(A\), whose reduction is ordinary
exactly when \(B\)'s is ordinary. An unconditional density-one result for
all three-balanced-pair actions would therefore settle density-one
ordinarity after extension for these threefolds without Fité's field
hypothesis. The simple-signature condition in Cantoral Farfán–Li–Mantovan–
Pries–Tang's [Theorem 4.18](https://arxiv.org/html/2508.11174v1#S4.SS5)
allows only one mixed pair. Noot's [Theorem 2.2](https://www.numdam.org/item/CM_1995__97_1-2_161_0.pdf)
uses an odd tensor power of the standard two-dimensional factors. A
direct sum of separate rank-two components does not meet that hypothesis.

The missing power of \(p\) in (1) for a single failing pair can occur for
actual Frobenius over finite fields at arbitrary large split primes.
Take \(F=\mathbb Q(\zeta_7+\zeta_7^{-1})\), with discriminant \(49\), and
take a large prime \(p\equiv1\pmod{28}\). For one place \(\mathfrak q\mid p\)
of \(F\), the real embedding lattice of \(\mathfrak q\) has covolume
\(7p\). Minkowski's convex body theorem on the box of radius
\(3p^{1/3}\) supplies a nonzero \(a\in\mathfrak q\) with all conjugates
of absolute value at most \(3p^{1/3}<2\sqrt p\). Its norm is a nonzero
multiple of \(p\), of magnitude at most \(27p\). Thus for \(p>27\) its
ideal is divisible by exactly one place over \(p\), with exponent one.
Also \(a\) generates \(F\), since the bound excludes a nonzero rational
integer divisible by \(p\) and \(F\) has prime degree.

Let \(\pi\) be a root of \(X^2-aX+p\). It is a \(p\)-Weil integer of
degree six, since it generates \(F\) through \(\pi+p/\pi=a\) and the
quadratic discriminant is totally negative. At the two unit-trace places
of \(F\) the slopes are \((0,1)\). At the remaining place the discriminant
has valuation one, giving a ramified quadratic extension and slopes
\((1/2,1/2)\). The local invariants of the division algebra in the
Honda–Tate classification are zero, since at the ramified quadratic
place the local degree multiplies \(1/2\) by two and the other slopes
are integers. There are no real places. The resulting simple isogeny
class over \(\mathbb F_p\) is that of a threefold \(B_p\) with commutative
endomorphism algebra \(\mathbb Q(\pi)\), containing \(F\) fixed by
Rosati. This uses [Tate's exposé, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0.pdf).
For \(A_p=B_p^2\) the action of \(E=F(i)\) as above has one nonordinary
conjugate pair and
\[
 p^{-3}N_{E/\mathbb Q}(a)=k^2/p,\qquad 1\le |k|\le27.
\]
For a concrete instance with \(p=13\), use \(a=\beta+3\), where
\(\beta^3+\beta^2-2\beta-1=0\). The three roots modulo \(13\) are
\(-6,-5,-3\), and \(N_{F/\mathbb Q}(a)=13\). This yields the exact
normalized value \(1/13\). The rational Frobenius polynomial of \(B_{13}\),
obtained by taking the norm of \(X^2-aX+13\), is
\[
 X^6-8X^5+58X^4-221X^3+754X^2-1352X+2197.
\]
It reduces to \(X^4(X^2-8X+6)\) modulo \(13\), also giving \(13\)-rank two.
These varying finite-field isogeny classes
show that the local divisibility and archimedean bound alone cannot
put the normalized trace norm in a fixed finite set in the
one-failing-pair case. They assert no failure of a density statement for
a fixed abelian variety over a number field.

The same examples also allow the maximal multiplicative rank for three
independent Frobenius pairs. The three conjugate discriminants
\(\sigma(a)^2-4p\) have valuation one at different places of the cyclic
cubic field \(F\) above \(p\), and zero at the other two places. Their
square classes are independent, so the splitting field over \(F\) has
an automorphism swapping either root pair while fixing the other two.
Each ratio of the two roots in a pair is nontorsion, since at another
place its valuations are zero and one. Dividing a multiplicative
relation by its image under one root swap therefore forces the two
exponents of that pair to be equal. Applying this to every pair shows
that all multiplicative relations with \(p\) are generated by the three
determinant identities. Thus maximal rank of these Frobenius tori
alone does not rule out the exact one-pair slope pattern.

## Simultaneous archimedean and auxiliary Frobenius data

**Proposition 4.** Fix \(F=\mathbb Q(\zeta_7+\zeta_7^{-1})\) and its three
real embeddings. Take a sufficiently large prime \(p\equiv1\pmod{28}\), a
prime \(\mathfrak q\mid p\) of \(F\), and any vector \(u\in[-2,2]^3\).
For finitely many auxiliary primes \(\ell\ne p\), fix positive integers
\(N_\ell\) and matrices
\[
 g_\ell\in\operatorname{GL}_2(O_F/\ell^{N_\ell}O_F),\qquad
 \det(g_\ell)=p\bmod\ell^{N_\ell}O_F.
\]
Put \(M=\prod_\ell\ell^{2N_\ell}\). Provided \(M p^{-1/6}\) is sufficiently
small, there is a simple threefold \(B_p/\mathbb F_p\) with real
multiplication by \(F\), nonordinary at exactly \(\mathfrak q\), whose
relative trace \(a\) and auxiliary Frobenius matrices satisfy
\[
 \max_\sigma\big|\sigma(a)/\sqrt p-u_\sigma\big|\ll_F M p^{-1/6},
 \qquad
 \pi\mid B_p[\ell^{N_\ell}]\sim g_\ell .                              \tag{3}
\]
The second assertion is conjugacy using the free rank-two
\(O_F/\ell^{N_\ell}O_F\)-module structure. On \(B_p^2\) the action by
\(E=F(i)\) has three balanced pairs, one of which is nonordinary.

**Proof.** For each \(\lambda\mid\ell\) choose a lift \(g^{(0)}_\lambda\)
over \(O_{F_\lambda}\) of the indicated residue matrix. Its upper-right
entry \(b\) can be chosen nonzero with
\(v_\lambda(b)\le v_\lambda(\ell^{N_\ell})\), by taking \(b=\ell^{N_\ell}\)
when the entry is zero in the quotient. Multiply the second row by
the unit \(p/\det(g^{(0)}_\lambda)\), which is one modulo \(\ell^{N_\ell}\).
The resulting matrix still lifts the residue matrix and has determinant
exactly \(p\). Write its trace as \(a^{(0)}_\lambda\). Prescribe
\(a\equiv a^{(0)}_\lambda\pmod{\ell^{2N_\ell}O_{F_\lambda}}\) at each such
place and require \(a\in\mathfrak q\). The Chinese remainder theorem
gives a nonempty coset modulo \(M\mathfrak q\).

The real embedding lattice of \(\mathfrak q\) has a basis with coordinates
bounded by \(O_F(p^{1/3})\). Indeed the cubic field of discriminant \(49\)
has class number one by Minkowski's bound, and Dirichlet's unit theorem
allows the logarithms of the absolute values of a generator's three
embeddings to be made uniformly close to \((\log p)/3\). Multiplying a
fixed integral basis by that generator gives the stated lattice basis.
The chosen coset has covering radius at most \(C_F M p^{1/3}\) in the
sup norm, by translating a fundamental parallelepiped. Approximate
\(\sqrt p\,u\) by a lattice point in the coset, first scaling the target
toward zero by a factor \(1-O_F(Mp^{-1/6})\) to ensure all resulting
conjugates lie strictly inside \((-2\sqrt p,2\sqrt p)\). If the chosen
point is zero, move by one lattice basis vector and take the initial
constant in the scaling large enough for this move as well. This proves
the first bound in (3).

Now \(0<|N_{F/\mathbb Q}(a)|<8p^{3/2}<p^2\). Since \(a\in\mathfrak q\),
it has valuation exactly one there and is a unit at the other two
primes over \(p\). The archimedean bound excludes a nonzero rational
integer divisible by \(p\), so \(a\) generates the cubic field. The
Honda–Tate argument above gives the claimed simple isogeny class, with
relative Frobenius polynomial \(X^2-aX+p\). Choose an initial member
with maximal rational endomorphism order by
[Waterhouse, Theorem 6.1](https://www.numdam.org/item/ASENS_1969_4_2_4_521_0.pdf).

It remains to realize the matrices. At \(\lambda\mid\ell\) write the
chosen determinant-\(p\) lift as \(\left(\begin{smallmatrix}x&b\\c&y\end{smallmatrix}\right)\).
Replace \(y\) by \(y+a-a^{(0)}_\lambda\) and \(c\) by
\(c+x(a-a^{(0)}_\lambda)/b\). The congruence at \(\lambda\) makes both
changes zero modulo \(\ell^{N_\ell}\). The resulting matrix has exact
trace \(a\) and determinant \(p\). Its characteristic polynomial is
separable and identifies it over \(F_\lambda\) with multiplication by
\(\pi\), so it specifies a Frobenius-stable free lattice in the rational
Tate module. After scaling, the direct sum of these lattices at
\(\lambda\mid\ell\) contains the original Tate lattice. Its quotient by
that lattice is a finite Frobenius-stable subgroup of \(B_p[\ell^\infty]\).
Quotienting by these subgroups for the finitely many auxiliary primes
realizes the specified lattices simultaneously. Their \(O_F\)-stability
retains the integral \(F\)-action. This changes only the member of the
isogeny class, so the Newton pattern and the archimedean trace survive.
\(\square\)

In particular, on any sequence of such split primes tending to infinity,
one may prescribe arbitrary limiting joint distributions of normalized
traces together with compatible auxiliary matrix conjugacy data at every
fixed finite level. To see this, specify target sequences for the two
data sets and impose a finite collection of levels expanding slowly
enough that \(M=o(p^{1/6})\). For any fixed collection, (3) eventually
matches its conjugacy data exactly and makes the error against every
continuous test function on the compact trace cube tend to zero. All
these finite-field examples still have one failing pair. The construction
concerns varying isogeny classes. Its implication is that weak
archimedean equidistribution and fixed-level auxiliary Frobenius
distributions, even jointly, leave room for the one-failure behavior at
split primes. A proof for a fixed characteristic-zero variety must use
additional compatibility across residue primes.
