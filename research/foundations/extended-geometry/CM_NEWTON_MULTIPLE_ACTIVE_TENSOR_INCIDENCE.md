# Ordinarity from support incidences of adjoint factors

Let \(A/K\) be a polarized abelian variety with a
Galois-stable CM coefficient field \(E\) on which Rosati
acts by conjugation. Assume that the geometric simple
adjoint factors of its Betti Mumford–Tate group \(G\) have
type \(A_1\) and that each \(E\)-embedding block of \(H^1\)
is irreducible under \(G^{\mathrm{der}}\). In this note a
Hodge-active label is a geometric adjoint factor on which
the Hodge cocharacter acts nontrivially. Assume that each
\(\mathbb Q\)-simple adjoint factor contains exactly one
Hodge-active label. The full geometric endomorphism
algebra can be larger than \(E\).

Let \(R=\prod_i R_i\) be the centroid of
\(\operatorname{Lie}(G^{\mathrm{der}})\), with the factors
indexed by the rational simple summands. Write \(S\) for
the active subset of \(\coprod_i\operatorname{Emb}(R_i)\),
and \(C_\lambda\) for the set of supporting factors of the
block labelled by an embedding \(\lambda\) of \(E\).
All embeddings are taken in a fixed
\(\overline{\mathbb Q}\subset\mathbb C\).
For \(g\) in the absolute Galois group of \(K\), let
\(\theta_g\) denote postcomposition by \(g\) on these
embedding sets and put
\[
                    S_g=\bigcup_{j\in\mathbb Z}\theta_g^j S .
\]

**Theorem.** The relative rank is \(2^k\) for some \(k\).
Every block is a tensor of \(k\) standard representations,
and every mixed pair is balanced. Outside a set of
natural density zero, at good absolute degree-one
primes \(v\) with arithmetic Frobenius \(g\),
reduction is \(\mu\)-ordinary relative to \(E\) exactly when
\[
                   \left|C_\lambda\cap S_g\right|\le 1
                    \quad\hbox{for every }\lambda .
\]
It is ordinary when this condition holds and \(\theta_g\)
preserves the Hodge multiplicities on the \(E\)-labels.
The two conditions give exact natural densities by finite
Chebotarev. In particular, ordinary density is positive
over \(K\) and is one after defining the tensors.
The full Newton vector in each block is
determined by the positive numbers \(1/|\operatorname{Orb}_g(\rho)|\)
for the labels \(\rho\) in \(C_\lambda\cap S_g\).

The criterion allows several mixed pairs with intersecting
supports. It specializes to the single-active calculation
and uses the rationality of Frobenius traces on each rational
adjoint summand separately.
[Single-active crystalline calculation](CM_NEWTON_TENSOR_ADJOINT_CRYSTALLINE.md),
[Kisin–Zhou, Theorem 1.2, arXiv v1](https://arxiv.org/pdf/2505.02165v1).

## Support and crystalline calculation

Each block is a tensor product of symmetric powers of
the standard modules of its supporting \(A_1\) factors.
The two possible weights on \(H^1\) imply that its support
contains at most one active label, with standard exponent
one. Indeed the number of distinct weights is at least
one plus the sum of the symmetric-power exponents over
active factors of the block. This follows by ordering
the weights in each factor and adding them successively
in the tensor. Every occurrence of a given active label
therefore has exponent one. Galois transitivity within
each rational simple factor, and equivariance of support
and exponent, give exponent one at every label.
Transitivity on the \(E\)-labels gives a common support
size \(k\). Hence a block is a standard tensor and its
Hodge vector is definite or balanced. It is balanced
exactly when its support meets \(S\).

The centroid and the projector onto the derived Lie
algebra are absolute Hodge tensors, as in the single-active
proof. Deligne–Blasius comparison lets us take the
crystalline realizations of these tensors at a good
degree-one prime. Over a sufficiently large finite
unramified \(p\)-adic field, write \(\Phi\) for the
semilinear crystalline Frobenius. It commutes with the
centroid and transports the labels by the arithmetic
postcomposition \(\theta_g\). A postcomposition orbit
cannot cross rational components, because each label is
one on its own rational primitive idempotent and zero
on the others. Thus each crystalline orbit has at most
one active label. On an orbit of length \(a\), the
adjoint slopes in each rank-three factor are
\((-d,0,d)\). Weak admissibility of the orbit summand
gives \(d=0\) on inactive orbits and \(ad\le1\) on
active ones, using Hodge weights \((-1,0,1)\) at the one
active label and zero weights elsewhere. If \(d>0\), the
first return has a unique lowest slope line, so \(ad\)
is integral. The only nonflat possibility on an active
orbit is \(d=1/a\).
[Blasius, Theorem (0.3)](https://archive.mpim-bonn.mpg.de/id/eprint/3859/1/preprint_1992_59.pdf),
[Crystalline calculation with endpoints and normalization](CM_NEWTON_TENSOR_ADJOINT_CRYSTALLINE.md#crystalline-slopes-of-the-adjoint-factors).

The next section excludes, up to density zero, complete
flatness in each rational adjoint summand. As there is at
most one active crystalline orbit in each such summand,
we can then write \(d(\rho)\) for the average of
\(\mathbf 1_S\) along the orbit of \(\rho\). On a coefficient
block, its Newton slopes are the Hodge mean of the
block's orbit plus all the signed sums
\[
                 \frac12\sum_{\rho\in C_\lambda}
                              \epsilon_\rho d(\rho),
                 \qquad \epsilon_\rho\in\{-1,1\}.
\]
The formula follows by taking a linear Frobenius power
that preserves the factors and lifting the inner
automorphisms of those factors to their two-dimensional
standard modules. The common shift agrees with the
averaged Hodge mean by endpoint equality. Supports are
equivariant under \(\theta_g\), so
\(\sum_{\rho\in C_\lambda}d(\rho)\) is the average of
\(|\theta_g^j(C_\lambda)\cap S|\), precisely the fraction
of mixed labels along the coefficient orbit. Its averaged
Hodge vector consists of the same extreme slopes, with
equal multiplicity whenever the extremes are distinct.
The signed sums attain only those extremes exactly when
at most one gap \(d(\rho)\) in the block is positive.
This proves the incidence test. The averaged Hodge vector
is integral exactly when its \(0\)-\(1\) Hodge vectors are
constant along the orbit. Since Newton lies above Hodge
with equal endpoints, any fractional Hodge reference
rules out ordinarity.

## Density over the original field

Choose a finite Galois extension \(L/K\) defining the
tensors and giving the \(G\)-valued representations of
Kisin–Zhou. Vasiu's theorem gives derived Mumford–Tate
comparison in type \(A_1\). For an auxiliary \(\ell\), let
\(\mathcal G_\ell\) be the full monodromy group over \(K\).
It normalizes the derived Lie algebra and permutes its
centroid labels. Fix \(M\) divisible by the exponent of
\(\operatorname{Gal}(L/K)\) and by the exponents of the
monodromy component group and the finite permutation
group on the centroid labels.

For each rational simple summand \(\mathfrak d_i\),
consider the function
\[
        u\longmapsto
        \operatorname{Tr}\bigl(\operatorname{Ad}(u^M)|_{\mathfrak d_i}\bigr)
        \quad\hbox{on }\mathcal G_\ell .
\]
At a good unramified prime \(v\), for any \(w\) over \(v\)
in \(L\), the \(M\)-th power of Frobenius at \(v\) is a
power of Frobenius at \(w\). Kisin–Zhou's \(G\)-conjugacy
compatibility at \(w\) makes the trace rational and equal
to its crystalline value. This applies even if \(w\)
has residue degree greater than one. The adjoint
eigenvalues are ratios of eigenvalues on \(H^1\), so they
are units away from \(p\) and have absolute value one
under all complex embeddings. If the crystalline
summand \(\mathfrak d_i\) is flat, the trace is also
\(p\)-integral. The trace is then an integer bounded by
\(\dim\mathfrak d_i\) in absolute value.
[Kisin–Zhou, Theorem 1.2 and its definition of rationality, arXiv v1](https://arxiv.org/pdf/2505.02165v1),
[Vasiu, Theorem 1.3.4(a) and Example 1.3.5(iv)](https://arxiv.org/pdf/math/0212066v7).

Each of the trace functions is nonconstant on every
component. Multiply a representative \(u_0\) by a
derived group element to make its action on the
derived group preserve a pinning. Its \(M\)-th power
then acts trivially on the derived group, since its
label permutation is killed and type \(A_1\) has no
diagram automorphism. The image of the sum of all
simple coroots is a cocharacter \(\nu\) fixed by this
pinned action. On \(u_0\nu(z)\) the trace on
\(\mathfrak d_i\) is
\([R_i:\mathbb Q](1+z^{2M}+z^{-2M})\). The union of
the finite-level loci over \(i\) is stable under
conjugation and contains no component. Algebraic
Chebotarev gives the required density-zero conclusion
over \(K\) itself.
[Pink, §7](https://people.math.ethz.ch/~pink/ftp/AMT-v3.pdf),
[Rajan, Theorem 3](https://doi.org/10.1155/S1073792898000142).

The support test only uses a finite Galois action:
take the normal closures over \(K\) of the images
of the coefficient and centroid fields. For conjugacy
invariance, let \(\beta(h)\) and \(\gamma(h)\) be the
actions of \(h\in G_K\) on the absolute Hodge centroid
and on \(E\). The simultaneous permutation of labels
given by postcomposition by \(h\) and precomposition
by \(\beta(h)^{-1}\) or \(\gamma(h)^{-1}\), respectively,
preserves support incidence and the Hodge signatures.
It preserves \(S\), since that subset is specified by
the positive part of the adjoint Hodge filtration.
Precomposition commutes with postcomposition, so this
simultaneous permutation transports the orbit test for
\(g\) to the test for \(hgh^{-1}\). Finite Chebotarev
therefore computes both densities by the proportions
of classes passing their respective tests.
The identity always passes. After defining the tensors,
\(\theta_g\) preserves \(S\) and the coefficient Hodge
multiplicities, since \(\beta(g)\) and \(\gamma(g)\)
are then trivial. The two tests always pass on that base.

## A fixed example with two mixed pairs and split coefficients

Let \(B=\operatorname{Jac}(C_{7,2})/\mathbb Q\) be the
verified Mumford-type example with adjoint centroid
\(\mathbb Q(\zeta_7)^+\). Give \(A=B^4\) the action of
\(E=\mathbb Q(\zeta_5)\) by its regular representation on
the multiplicity space. Equip that space with the
positive definite form \((x,y)\mapsto
\operatorname{Tr}_{E/\mathbb Q}(x\bar y)\).
After scaling to an integral form in a rational basis,
it gives a polarization on \(B^4\) whose Rosati involution
on \(E\) is conjugation. Each relative block is a copy
of the rank-eight tensor of \(B\), so there are two
balanced \((4,4)\) pairs. The entire coefficient action
is defined over \(\mathbb Q\).

The coefficient and adjoint fields are Galois and
disjoint by their degrees four and three. Among primes
split completely in \(E\), the two nonidentity Frobenius
classes of the cubic adjoint field still have relative
density \(2/3\). They give three positive gaps on every
block outside a set of density zero, so fail the
incidence test. Thus even at split coefficient primes,
ordinary and \(\mu\)-ordinary reduction here have relative
density \(1/3\).
[Exact verification of the fixed Jacobian and its reduction law](CM_NEWTON_TENSOR_ADJOINT_CRYSTALLINE.md#ordinary-density-and-an-exchanging-cubic-example).

## Two Shimura surfaces

The next constructions give fixed abelian varieties over
number fields whose mixed pairs have distinct supports.
The first realizes two active labels in a single rational
simple adjoint factor. The second realizes three mixed
pairs with intersecting supports under the theorem.

**Two active labels in one rational component.** Let
\(R_0=\mathbb Q(\zeta_{13})^+\), cyclic of degree six, and
let \(F\) be its cubic subfield. Number its real embeddings
cyclically modulo six. Take a quaternion algebra \(D/R_0\)
split at the real labels \(0,1\), ramified at the other
four real labels, and split at all finite places. Its
corestriction \(C=\operatorname{Cor}_{R_0/F}(D)\)
is a central simple algebra of degree four over \(F\).
At real labels of \(F\) its tensor supports are the
fibers \(\{i,i+3\}\). The Brauer class is split by
\(E=F(i)\): at finite places the quaternion factors
are split, and above real places \(E\) is complex.
Use a simple \(C\otimes_F E\)-module \(W=E^4\).
The represented group
\(\operatorname{Res}_{R_0/\mathbb Q}\mathrm{SL}_1(D)\)
acts through the corestriction on \(W\).
The two split labels occur in different supports.
The resulting abelian varieties have dimension \(12\),
with two \((2,2)\) pairs and one definite pair.
In this construction the single rational simple summand
has two active labels. Potential ordinary density one
follows from the bound for two active labels.
[Adjoint trace bound for two active labels](CM_NEWTON_ADJOINT_PROJECTIVE_TRACE_BOUNDS.md).

**Intersecting supports in two rational components.**
Let \(F=\mathbb Q(\zeta_{15})^+\), cyclic of degree four
with generator \(\iota\). Use two independent copies
\(D_1,D_2\) of a quaternion algebra over \(F\) split at
real label \(0\), ramified at the other three real
places and at one prime above \(59\), and split
elsewhere. These ramification sets are permitted by
quaternion reciprocity. The rational prime \(59\) splits
completely in \(F\), since it is \(-1\) modulo \(15\),
and \(-1\) is nonsquare modulo \(59\).
Thus \(E=F(i)\) splits both quaternion algebras and their
twists by \(\iota\). Take the central simple algebra
\[
           C=D_1\otimes_F {}^\iota D_2
                   \otimes_F {}^{\iota^2}D_2
\]
using twists indexed so that the label \(j\) has second
and third factors \(j+1,j+2\). Use its split simple
module \(W=E^8\) after extending to \(E\). The two
groups \(\operatorname{Res}_{F/\mathbb Q}\mathrm{SL}_1(D_i)\)
act by the tensor of their indicated factors. Write
\(U_j,V_j\) for the labels of these two independent
groups. The four conjugate coefficient supports are
\[
       \{U_j,V_{j+1},V_{j+2}\},\qquad j\in\mathbb Z/4 .
\]
The active labels are \(U_0,V_0\). Each support contains
at most one of them. Three pairs have signature \((4,4)\)
and the fourth is definite. The abelian dimension is
\(32\). Different supports can share a \(V\)-factor.

**Hodge and arithmetic realization.** For either
representation, take the image of the displayed
simply connected groups and the commuting torus
\[
    T=\{t\in\operatorname{Res}_{E/\mathbb Q}\mathbb G_m
                        \mid t\bar t\in\mathbb G_m\}
\]
acting by multiplication. On \(\mathbb C^\times\),
at every split real quaternion factor take the standard
rotation by \(z/|z|\), and take the identity at compact
factors. On the torus coordinates take \(|z|\) at each
mixed pair and \(z\) at each definite pair. This satisfies
the common norm condition on \(T\). Every mixed block
has one split factor, so its action is the real standard
matrix for \(z\) on that factor times the identity on
the remaining factors. The resulting map to
\(\operatorname{GL}(W_\mathbb R)\) is algebraic over
\(\mathbb R\) and gives the Hodge weights of \(H^1\).

On a mixed real block use the standard alternating
form on the split factor tensored with the positive
real inner product on the compact tensor of the other
factors. On a definite block use the imaginary part of
a positive Hermitian form invariant under the compact
factors. Choose the signs so that the complex structure
is polarized. Both prescriptions give the adjoint
\(e\mapsto\bar e\) on \(E\), and the similitude on the
central torus is its common norm. The space of forms
with the required invariance and adjoint is rational,
so rational approximation to this real form in it
gives a rational polarization. The adjoint action of
the Hodge homomorphism gives the Cartan involution at
each split factor and the identity at compact factors.
Its weight is the rational scalar homomorphism. Every
rational adjoint factor has a split active label.
This gives Shimura data of Hodge type, with two active
real factors in each construction. The same
corestriction and subset constructions for one
noncompact label appear in Moonen's classification.
[Moonen, Examples 3.1–3.2 and Construction 5.6](https://doi.org/10.4153/S0008414X24001159).

At sufficiently small level their canonical models
carry the abelian schemes from the Hodge embeddings.
The connected geometric monodromy is full derived.
Indeed arithmetic uniformization gives a lattice
whose projection to the noncompact real factors is
Zariski dense by Borel density. The rational connected
Zariski closure projects onto every geometric factor
by transitivity within each rational simple group.
A subdirect product could identify two simple factors
only by a diagonal. No split real factor can be
identified with a compact one, and projection to
the noncompact factors rules out diagonals between
the split factors. The rational closure therefore
contains every split factor separately, and then
every factor by Galois transitivity.
[Borel density](https://doi.org/10.2307/1970150),
[Borel–Harish-Chandra, arithmetic lattices](https://www.mathi.uni-heidelberg.de/~wienhard/retreat10/references/borel_harishchandra_2.pdf).

Hilbert irreducibility supplies number-field fibers
retaining this connected derived monodromy at an
auxiliary prime \(\ell\). One may use an open normal
pro-\(\ell\) subgroup \(U\) of the arithmetic image and
the finite quotient by its open Frattini subgroup.
After an étale chart over affine space, specialize
the corresponding connected finite cover to a
connected fiber, allowing the residue field to be
finite over the field of the chosen component.
Surjectivity onto the finite quotient implies that
the specialized image contains \(U\). The retained
derived group and the Hodge cocharacter give geometric
endomorphism algebra exactly \(E\) at these fibers.
Distinct coefficient supports distinguish blocks
over different real labels of \(F\). Within a
conjugate pair the two central characters differ
at the definite label where the Hodge action is the
two conjugate scalar characters. Galois transitivity
then distinguishes the orientations at every pair.
[Noot, §5 for the specialization argument](https://irma.math.unistra.fr/~noot/publications/jpaa.pdf).
