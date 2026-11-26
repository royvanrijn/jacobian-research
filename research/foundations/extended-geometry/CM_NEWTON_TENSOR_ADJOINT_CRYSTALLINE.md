# Newton slopes for one mixed coefficient pair of type A1

Let \(A/K\) be a polarized abelian variety over a number field.
Let \(E\subset\operatorname{End}^0(A_{\overline K})\) be a
Galois-stable CM field on which Rosati acts by complex
conjugation, and put \(F=E^+\). Use a fixed complex embedding
of \(K\) for the Hodge signatures. Assume exactly one
conjugate coefficient pair is mixed and every other pair
is definite. Assume that the geometric simple factors of
the adjoint Betti Mumford–Tate group are of type \(A_1\),
and that each \(E\)-embedding block of \(H^1\) is irreducible
under the derived group over \(\overline{\mathbb Q}\). The coefficient
field can be part of a larger endomorphism algebra.

**Theorem.** Write \(R\) for the centroid of the derived
Betti Lie algebra. It is a field containing \(F\). For
\(k=[R:F]\), the relative rank is \(2^k\), and the mixed
pair has signature \((2^{k-1},2^{k-1})\). There is exactly
one active adjoint Hodge label \(\tau:R\hookrightarrow
\overline{\mathbb Q}\). The field of definition of the
adjoint labels is the finite Galois extension \(K\tau(R)\).
The connected projective spectrum field of the coefficient
blocks is \(F\).

At a good absolute degree-one prime \(v\) outside a finite
set, write \(\beta\) for its Frobenius action on \(R\) and put
\[
 f=\operatorname{ord}(\beta),\qquad
 s=\operatorname{ord}(\beta|_F),\qquad r=f/s.
\]
Then \(r\mid k\). There are two possible relative Newton
patterns. In the first, each coefficient block has all
slopes equal. This case occurs only on a set of natural
density zero. In the second, for each coefficient block
over the restricted orbit of the active label, the slopes are
\[
        m+\frac{r-2j}{2f},
        \qquad 0\le j\le r,
        \qquad \text{with multiplicity }2^{k-r}\binom rj .
\]
Here \(m\) is the block's averaged Hodge mean. Every other
block has all slopes equal to its averaged Hodge mean.
In particular, \(\mu\)-ordinary reduction has natural
density equal to the proportion of \(\beta\) in
\(\operatorname{Gal}(K\tau(R)/K)\) with \(f=s\).

The binomial Newton polygons for a single Mumford tensor
over a base containing the field of its adjoint labels
are classical. Noot treats the cubic case and Sheng–Zuo
give the general tensor formula. Pink's density-one result
after extension for geometric endomorphism ring \(\mathbb Z\)
includes the classical \(A_1\) case and uses finite ranges
of Frobenius traces, including adjoint traces. The calculations
here use the rank-three crystalline adjoint factors and apply
to several CM coefficient blocks over the original base.
[Noot, Proposition 2.2](https://irma.math.unistra.fr/~noot/publications/jpaa.pdf),
[Sheng–Zuo, Theorem 1.2](https://arxiv.org/pdf/1106.3505),
[Pink, §7](https://people.math.ethz.ch/~pink/ftp/AMT-v3.pdf).

## Betti support and label fields

An irreducible block is the tensor product of the
\(\operatorname{Sym}^{a_i}\) of the standard modules of its
supporting \(A_1\) factors, with \(a_i>0\). One can take a
power of the Hodge cocharacter to lift it to the product
of the connected center and the simply connected factors
when counting weights. An active factor on such a block
contributes \(a_i+1\) distinct weights. A sum of two finite
ordered weight sets has at least the sum of their sizes
minus one distinct weights: add the first element of the
second set along the first, then add its remaining elements
to the last element of the first. Iterating, a block with
active factors \(i\) has at least \(1+\sum_i a_i\) weights.

Every active factor is visible in some block by faithfulness.
It can occur only in the mixed pair, where there are exactly
two weights. The conjugate blocks have the same support by
polarization. It follows that there is one active factor,
and its exponent in the tensor is one. A rational simple
adjoint quotient that omits this factor kills the Hodge
cocharacter and its complex conjugate, hence the whole
Hodge homomorphism. Mumford–Tate minimality excludes such
a quotient. Thus the centroid is a field \(R\), with the
geometric simple factors labelled by \(\operatorname{Emb}(R)\).

The active label is incident to precisely the mixed
conjugate coefficient pair. Betti support incidence and
the exponents \(a_i\) are Galois-equivariant over
\(\mathbb Q\). Transitivity on the labels of \(R\) therefore
gives exactly one incident conjugate pair and exponent one
at every label. The incidence map to
\(\operatorname{Emb}(F)\) is a map of transitive Galois sets.
Its fibers are the supports of the blocks. The stabilizer
inclusion identifies \(F\) with a subfield of \(R\), and
the fibers have size \(k=[R:F]\). This proves the relative
rank, and a single active standard factor gives the two
Hodge weights with equal multiplicity. This is an incidence
of the kind appearing in the subset construction for
Shimura curves.
[Moonen, Example 3.2 and Remark 3.3](https://doi.org/10.4153/S0008414X24001159).

Let \(\mathfrak d\subset\operatorname{End}(H^1)\) be the
derived Lie algebra. The trace pairing on endomorphisms is
nondegenerate on \(\mathfrak d\): its simple factors are
orthogonal and each has nonzero trace index. The orthogonal
projector onto \(\mathfrak d\) and the centroid acting on
its image are Hodge tensors. In type \(A_1\), Vasiu's
Mumford–Tate theorem gives Betti–étale comparison of the
derived groups at every auxiliary prime. Galois consequently
fixes the projector and acts on the centroid by field
automorphisms \(\beta(g)\). On the de Rham realization the
nonzero positive filtration piece of \(\mathfrak d\) is
supported at the unique label \(\tau\). Since the filtration
is defined over \(K\), the label is fixed under the
semilinear action, so
\[
                g\circ\tau=\tau\circ\beta(g).
\]
Thus the fixed field of the kernel of \(\beta\) is
\(K\tau(R)\). The support incidence also gives
\(\beta(g)|_F=\gamma(g)|_F\) for the coefficient action
\(\gamma\) on \(E\).
[Deligne, Proposition 2.9 and Theorem 2.11](https://www.jmilne.org/math/Documents/Deligne82.pdf),
[Vasiu, Theorem 1.3.4(a) and Example 1.3.5(iv)](https://arxiv.org/pdf/math/0212066v7),
[Label field argument in rank four](CM_NEWTON_TENSOR_ADJOINT_LABEL_FIELD.md#field-of-the-adjoint-labels).

The projective action on a block is the standard tensor of
the factors supporting that block. Thus the two conjugate
blocks above a given \(F\)-label have the same projective
spectrum. Different \(F\)-labels have disjoint supports, so
a torus in one supporting factor distinguishes their spectra.
The fibers of the [connected projective spectrum field](CM_NEWTON_PROJECTIVE_SPECTRUM_FIELDS.md)
are exactly the restriction fibers to \(F\).

## Crystalline slopes of the adjoint factors

Take a good absolute degree-one prime with \(K_v=\mathbb Q_p\).
Exclude the finite set ramifying in normal closures of the
coefficient and label fields and the fields of definition
of the tensors. Choose a finite unramified \(U/\mathbb Q_p\)
containing their local images and all embeddings of \(E\)
and \(R\). Write \(\sigma\) for its arithmetic Frobenius.
By crystalline comparison at good reduction, the module
\(D_{\mathrm{cris}}(H^1)\) identifies with the crystalline
cohomology of the reduction. Over \(U\), crystalline
Frobenius \(\Phi\) is \(\sigma\)-semilinear. The projector
onto the derived Lie algebra and the centroid tensors
compare with their filtered de Rham realizations by the
compatibility with absolute Hodge classes on abelian
varieties. In particular \(\Phi\) commutes with the \(R\)-action
as a semilinear map and takes the label \(\rho\) to
\(\sigma\circ\rho\).
[Tsuji, Theorem 2.2.1 in the smooth proper case](https://emis.de/ft/9986),
[Blasius, Theorem (0.3)](https://archive.mpim-bonn.mpg.de/id/eprint/3859/1/preprint_1992_59.pdf).

The Hodge weights on the active adjoint factor are
\((-1,0,1)\), since the active factor acts on its standard
module with a Hodge weight difference of one. On every
other factor the weights are \((0,0,0)\).
The direct sums over label orbits are filtered
\(\Phi\)-stable summands and are weakly admissible.
The defining inequality for each summand, together with
equality of the total endpoints, gives equality of the
endpoints for each summand.
[Colmez–Fontaine, Corollary following Théorème 4.3](https://www.imo.universite-paris-saclay.fr/~fontaine/CF.pdf).

On an orbit of length \(a\), the slopes of the first
return \(\Phi^a\) on one rank-three factor are
\((-ad,0,ad)\) for \(d\ge0\). Indeed a power of the
return that is linear acts by an automorphism of
\(\mathfrak{sl}_2\). Its eigenvalues have the form
\(u,1,u^{-1}\). The orbit summand for \(\Phi\) thus
has slopes \((-d,0,d)\), each repeated \(a\) times.
If \(d>0\), the first return has a rank-one lowest
slope summand over \(U\). To see this, take the
least-valuation eigenspace of a linear power. Its
eigenvalue is the unique root of that valuation, so it
lies in \(U\), and the line is stable under the
semilinear first return. The return on the line has
integral slope because \(U/\mathbb Q_p\) is unramified.
Thus \(ad\) is an integer. The lowest slope subspace
for \(\Phi\) in the full orbit has dimension \(a\).
Its Hodge number is at least zero if the orbit is
inactive, and at least \(-1\) on the active orbit.
Weak admissibility gives \(ad=0\) off the active orbit
and \(ad\le 1\) on it.

By \(g\tau=\tau\beta(g)\), the active crystalline orbit is
\(\{\tau\circ\beta^j\}\), where \(\beta\) is the arithmetic
Frobenius action at \(v\). Its length is \(f\) because
\(\tau\) is injective. The restricted orbit in
\(\operatorname{Emb}(F)\) has length \(s\). Each of its
fibers receives \(r=f/s\) active orbit labels. Moreover
\(\beta^s\) acts freely with order \(r\) on each fiber of
\(\operatorname{Emb}(R)\to\operatorname{Emb}(F)\)
by precomposition, so \(r\mid k\).

To recover the \(H^1\) slopes, use the \(N\)-th power of
\(\Phi\), with \(N=[U:\mathbb Q_p]\), on each coefficient
block. It preserves its individual supporting factors.
It acts, up to a scalar, as the tensor product of lifts
of the \(N\)-th adjoint action on those factors. This follows
from the inner automorphisms of \(\mathrm{PGL}_2\) and
Schur's lemma. On a factor with slopes \((-d,0,d)\)
the two lifted slopes have difference \(Nd\).
In the nonflat case \(d=1/f\) on all \(f\) labels of the
active orbit and zero elsewhere. Dividing by \(N\) and
counting the \(r\) active two-dimensional factors yields
the binomial formula. In both cases the scalar mean
agrees with the averaged Hodge mean by endpoint equality
on each coefficient orbit.

Every oriented coefficient orbit above the restricted
active orbit has mixed fraction \(1/s\), since every
restricted label occurs equally often and there is just
one mixed pair. Its averaged Hodge vector has slopes
\(m-1/(2s)\) and \(m+1/(2s)\), each of multiplicity
\(2^{k-1}\). The nonflat binomial vector equals this
exactly when \(r=1\). If \(r>1\), it has \(r+1\)
distinct slopes. In the flat case the active orbit
still has a two-slope reference, so this case also
fails \(\mu\)-ordinarity.

## The first pattern has density zero

Fix an auxiliary \(\ell\). Choose one positive integer \(M\)
annihilating the two finite Galois actions on \(R\) and
\(E\). At a good prime away from \(\ell\), let \(\pi_v\)
be the Frobenius endomorphism of the special fiber.
The power \(\pi_v^M\) commutes with its geometric action
of \(E\). Its relative characteristic polynomial is
defined over \(E\) and is the same in the étale and
crystalline realizations. One way to see this is to
recover each relative power trace from the rational
traces of \(e\pi_v^{Mj}\), for \(e\) in a rational basis
of \(E\), by the trace pairing on \(E\). For an actual
endomorphism \(u\) the characteristic polynomial on
\(H^1\) is determined by \(\det(n-u)=\deg(n-u)\) for
large integers \(n\), using the exterior algebra
description of the cohomology of an abelian variety.
Scaling clears denominators for a rational endomorphism.
Newton identities recover the relative polynomial from
these relative traces.

On the auxiliary étale realization, the adjoint Lie
algebra belongs to the centralizer of \(E\) in
\(\operatorname{End}(H^1)\). The eigenvalues of
\(\operatorname{Ad}(\pi_v^M)\) on it are therefore ratios
of roots within individual relative characteristic
polynomials. If all relative Newton blocks are flat,
each such ratio has valuation zero under every
embedding into \(\overline{\mathbb Q}_p\), since
conjugating a relative polynomial gives another relative
polynomial. The ratios are units at other finite places
because the roots are Weil numbers and are units away
from \(p\). All complex conjugates of a ratio have absolute
value one. By Kronecker, the ratios are roots of unity.
Their degrees are bounded by the degree of a splitting
field of the absolute Frobenius polynomial, hence by
\((2\dim A)!\). Thus their sum, the adjoint trace of
\(\pi_v^M\), lies in a fixed finite set.

The regular function \(g\mapsto\operatorname{Tr}_{\mathfrak d}
(\operatorname{Ad}(g^M))\) on the auxiliary algebraic
monodromy is nonconstant on every component.
Over an algebraic closure, multiply a representative of
the component by a derived element so that its adjoint
action preserves a pinning. Denote this representative by
\(g_0\). The automorphism \(\operatorname{Ad}(g_0^M)\)
fixes every factor and the pinning. Each factor has type
\(A_1\), so this automorphism is trivial. Let \(\nu\) be
the sum of the simple coroots on the simply connected
cover, followed by the map to the derived group. It is
fixed by the pinned action. At \(g_0\nu(z)\) the trace is
\([R:\mathbb Q](1+z^{2M}+z^{-2M})\).
Algebraic Chebotarev now excludes the finite set of trace
levels outside a set of natural density zero. Finite
Chebotarev for \(K\tau(R)/K\) gives the stated density
since the order condition is invariant under conjugation.
[Pink, §7 on the finite adjoint-trace method](https://people.math.ethz.ch/~pink/ftp/AMT-v3.pdf),
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

## Ordinary density and an exchanging cubic example

If \(F=\mathbb Q\), every averaged Hodge vector is balanced
with integral slopes, so ordinary and \(\mu\)-ordinary
reduction coincide. If \(F\ne\mathbb Q\), let \(L_E/K\)
define the \(E\)-action. Ordinary reduction has natural
density
\[
                 \frac{1}{[L_EK\tau(R):K]}.
\]
Indeed on \(F\)-labels the Hodge-preserving permutation
\(\theta\gamma^{-1}\) fixes the sole mixed label, where
\(\theta\) denotes \(\Phi\)'s label permutation and \(\gamma\)
is its compatible coefficient action. If the action on
\(F\) is nontrivial, the corresponding orbit contains a
definite label. The resulting fractional reference slopes
rule out ordinarity by Newton above Hodge. If the action
on \(F\) is trivial but \(\gamma\) is complex conjugation,
the definite labels have fractional reference slopes
\(1/2\) in their orbits. For \(\gamma=1\) all reference
slopes are integral, and \(\beta|_F=1\). Then, up to the
set of density zero above, ordinary reduction is equivalent
to \(\gamma=\beta=1\). Applying the same formula after
any finite extension of \(K\) shows that \(L_EK\tau(R)\)
is the least extension over which ordinary reduction
has density one. For \(F=\mathbb Q\), the least such
extension is \(K\tau(R)\), since \(s=1\) in that case.
[Rank-four ordinary-density calculation](CM_NEWTON_TENSOR_ADJOINT_LABEL_FIELD.md#exact-density-in-the-quadratic-corestriction-case).

In the first cubic-support case \(k=3\), every \(\beta^s\)
on a support fiber is the identity or a three-cycle.
The latter components give the four centered slopes
\(-1/(2s),-1/(6s),1/(6s),1/(2s)\) with multiplicities
\(1,3,3,1\) on each mixed orbit, away from the first
pattern of density zero.

**Example over \(\mathbb Q\).** Let \(B=\operatorname{Jac}(C_{7,2})\)
be the \(t=2\) member of the hyperelliptic family of
Bouchet, Hanselman, Pieper, and Schiavone. Their
Theorem 4.15 gives Mumford type for its non-CM members,
with adjoint centroid \(R=\mathbb Q(\zeta_7)^+\).
[Bouchet–Hanselman–Pieper–Schiavone, equation (1.1) and Theorem 4.15](https://arxiv.org/html/2510.00093v2).

The [exact PARI/GP replay](../scripts/verify_mumford_tensor_specialization.gp)
verifies that \(B\) has good ordinary reduction at \(13\)
and \(29\). Both Frobenius polynomials are irreducible of
degree eight with splitting field degree \(24\). Each
nontrivial root-of-unity ratio of two distinct roots
would have order \(n\ge2\) with \(\varphi(n)\mid24\).
The bound \(\varphi(n)^2\ge n/2\) follows by multiplying
the prime-power factors of \(n/\varphi(n)^2\).
Thus \(n\le1152\), and the replay checks that the
characteristic polynomials of the \(n\)-th powers remain
irreducible for all such \(n\). Hence no distinct roots
coalesce under any positive power. Each reduction is
absolutely simple and its geometric endomorphism algebra
is the degree-eight Frobenius field: over every extension
its Frobenius still generates that field, and Tate's
centralizer theorem gives equality. The two field
discriminants are \(139546236594961\) and
\(234444145241761\). Specialization of endomorphisms
is injective at good reduction. It also precludes a proper
geometric isogeny factor of \(B\), whose rational idempotent
would give one on the absolutely simple reduction. Hence
geometric CM for \(B\) would identify one degree-eight endomorphism field with
both of these distinct fields. This proves that \(B\) is
non-CM.
[Tate, main theorem in §1](https://pazuki.perso.math.cnrs.fr/index_fichiers/Tate66.pdf).

Take \(A=B^2\) with its product polarization and
\(E=\mathbb Q(i)\) acting by \(i(x,y)=(-y,x)\). Rosati
acts on \(E\) by complex conjugation. Both relative
blocks of \(H^1(A)\) are copies of the irreducible
rank-eight tensor \(H^1(B)\). The unique pair is of
signature \((4,4)\) and \(K\tau(R)=\mathbb Q(\zeta_7)^+\)
is cyclic of degree three. Thus \(A\) and \(B\) have
ordinary density \(1/3\). At primes in the two nonidentity
classes, outside a set of density zero, the \(H^1(B)\)
slopes are \(0,(1/3)^3,(2/3)^3,1\), where the exponents
denote multiplicities. The example realizes a cubic
support exchange in relative rank eight.
