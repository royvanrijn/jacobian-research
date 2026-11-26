# The adjoint factor field in the rank-four tensor case

Let \(A/K\) be a polarized abelian variety over a number field
with a Rosati-compatible geometric CM coefficient field \(E\)
preserved by \(G_K\). Write \(F=E^+\). Assume relative rank
four, one conjugate Hodge pair of type \((2,2)\), and every
other pair definite. Write \(M\) for its Betti Mumford–Tate group
for a fixed embedding \(K\hookrightarrow\mathbb C\), and define
\(\mathfrak d_B=\operatorname{Lie}(M^{\mathrm{der}})\).
Let \(R\) be its centroid.
Assume that at an auxiliary prime \(\ell\), the derived
connected algebraic monodromy group is \(M^{\mathrm{der}}_{\mathbb Q_\ell}\)
under the Betti–étale comparison. Assume also that every
relative block is irreducible under \(M^{\mathrm{der}}\)
over \(\overline{\mathbb Q}\).

The centroid is the algebra of endomorphisms of a Lie algebra
commuting with all its adjoint operators. Since \(\mathfrak d_B\)
is semisimple, \(R\otimes\overline{\mathbb Q}\) is the product
of scalars on its simple factors. Thus the embeddings of
\(R\) label these factors.

**Theorem.** The centroid \(R\) is a field. There is a
distinguished embedding \(\tau:R\hookrightarrow\overline{\mathbb Q}\) at which the
adjoint Hodge cocharacter is nontrivial. It is the only
such label. The field of definition of the simple-factor
labels in the auxiliary monodromy is exactly
\(K\tau(R)\), a Galois extension of \(K\).

In the quadratic corestriction case, \(F\subset R\) has
\([R:F]=2\) and every relative block is the tensor of the
two standard representations at the labels of \(R\) above
the block's restriction to \(F\). For an element of the finite Galois action on
the factor labels write \(\beta\in\operatorname{Aut}_{\mathbb Q}(R)\).
Then \(\mu\)-ordinary reduction has natural density
equal to the proportion of \(\beta\) satisfying
\[
              \operatorname{ord}(\beta)
                 =\operatorname{ord}(\beta|_F).
\]
Here \(\mu\)-ordinarity is relative to the given \(E\)-action.
In particular \(\mu\)-ordinary reduction has density one
whenever \(K\) contains the adjoint reflex field \(\tau(R)\).
If the Galois action on \(E\) is contained in
\(\{1,c\}\), with \(c\) complex conjugation, this formula
is \(1/[K\tau(R):K]\).

When all the geometric adjoint factors are of type \(A_1\),
the [quadratic support incidence](#support-incidence-in-type-a1)
follows from the Hodge and block assumptions, and the
[derived comparison](#derived-comparison-for-rank-one-factors)
holds at every \(\ell\).

The field criterion gives a test for tensor components
in the rank-four setting of the
[linked Pfaffian theorem](CM_NEWTON_LINKED_RANK_FOUR_PFAFFIAN.md).
It uses the field of the simple factors as well as the
field of definition of the geometric endomorphisms.
The use of a distinguished Hodge eigenspace to recognize
an arithmetic label field has a precedent for the
transcendental endomorphism fields of \(K3\) surfaces.
[Elsenhans–Jahnel, Theorem 4.5 and §§6–7](https://arxiv.org/abs/2304.07585).

## Field of the adjoint labels

The trace pairing \((X,Y)\mapsto\operatorname{tr}(XY)\) on
\(\operatorname{End}(H^1_B(A))\) is nondegenerate on
\(\mathfrak d_B\). Indeed distinct simple factors are
orthogonal, and its restriction to each factor is a
nonzero multiple of the Killing form, by faithfulness
and the trace index of a nontrivial representation.
Let \(P\) be the rational orthogonal projector onto
\(\mathfrak d_B\). It is a Hodge tensor because \(M\)
normalizes \(\mathfrak d_B\), preserves the pairing
and thus preserves the orthogonal direct sum. All
elements of the centroid are Hodge tensors on
\(\mathfrak d_B\) because the connected group \(M\)
acts trivially on the finite algebra automorphism
group of \(R\). Composing with \(P\) regards them as
tensors on \(H^1_B(A)\).

Deligne's theorem gives absolute Hodge realizations
of these tensors. On the \(\ell\)-adic realization,
\(G_K\) preserves the derived Lie algebra because it
normalizes the connected monodromy group. It also
preserves the trace pairing and therefore fixes
\(P\). The realization on absolute Hodge tensors
is faithful, so \(P\) is \(G_K\)-invariant and
its de Rham image is defined over \(K\). In particular
the de Rham realization of \(\mathfrak d_B\) has
a filtration over \(K\). The centroid tensors are
stable under Galois because they are exactly the
endomorphisms of this Lie algebra commuting with
its adjoint operators. Write \(\beta(g)\) for Galois
conjugation on their rational algebra \(R\).
On the \(\ell\)-adic realization this gives exactly
the permutation of the simple-factor labels.
[Deligne, Proposition 2.9 and Theorem 2.11](https://www.jmilne.org/math/Documents/Deligne82.pdf).

The Hodge cocharacter acts by scalars on the
definite coefficient blocks. A nonzero projection
on a simple factor is visible projectively on any
block that sees it. Indeed the projective homomorphism
on the supporting adjoint factors is faithful.
Thus only factors supporting the balanced conjugate
pair can be Hodge-active. There is at least one
active factor, since the mixed representation has
two distinct weights. A block of dimension four
has at most two supporting factors. If it has two,
its representation is a tensor product of two
standard two-dimensional representations. A
cocharacter with nonzero projection on both would
give at least three distinct weights on the tensor,
which is impossible for \(H^1(A)\).
The balanced conjugate block sees the same factors
by polarization. This proves uniqueness of the
active factor, giving \(\tau\).
[Hui, §3.1 for the four-dimensional possibilities](https://arxiv.org/pdf/2208.04004).

A rational simple adjoint quotient on which the Hodge
cocharacter is trivial also kills the complex conjugate
cocharacter. The Hodge homomorphism would then factor
through the kernel of that quotient, a proper rational
subgroup of \(M\). Mumford–Tate minimality rules this out.
Consequently every rational adjoint simple factor contains
an active label. There is just one such label, so the
geometric simple factors form one Galois orbit, and the
centroid \(R\) is a field.

The first positive Hodge filtration piece of the
de Rham Lie algebra is nonzero and contained in
the single \(\tau\)-eigenspace for \(R\). For
\(g\in G_K\), application of the semilinear de Rham
action to a coefficient eigenvector sends its label
\(\tau\) to \(g\circ\tau\circ\beta(g)^{-1}\).
The filtration is defined over \(K\), so uniqueness
of the label gives \(g\circ\tau=\tau\circ\beta(g)\).
As \(\tau\) is injective, the kernel of \(\beta\)
is exactly the subgroup fixing \(\tau(R)\)
pointwise. Its fixed field is \(K\tau(R)\).

## Support incidence in type A1

Assume all the geometric adjoint factors are of type \(A_1\).
The mixed block cannot be a symmetric cube of a standard
two-dimensional representation. Such a representation gives
four distinct weights for every noncentral Hodge cocharacter,
whereas \(H^1(A)\) has two weights. The mixed conjugate pair
therefore has tensor-standard support. The active factor
belongs to this support and to no definite block, since its
noncentral cocharacter would give distinct weights on any
block that it supports. Polarization identifies the supports
of conjugate blocks, so the active factor is incident to
exactly one conjugate coefficient pair.

The incidence relation is Galois-equivariant over \(\mathbb Q\),
because both \(M^{\mathrm{der}}\) and the coefficient action on
the Betti representation are rational. Galois transitivity
on the labels of \(R\) now gives exactly one incident
conjugate pair at every factor. If a block were a symmetric
cube, one could transport its sole support factor to the
active label and obtain a symmetric cube on a mixed block.
Thus every block is tensor-standard. Sending a factor to
its incident pair gives a Galois-equivariant map
\(\operatorname{Emb}(R)\to\operatorname{Emb}(F)\) with
two-element fibers. For a chosen label and its image, the
inclusion of their stabilizers identifies \(F\) with a
subfield of \(R\) of degree two, and the map is restriction.
For the description of support subsets and their fields
in Shimura curves, see
[Moonen, Example 3.2 and Remark 3.3](https://doi.org/10.4153/S0008414X24001159).

## Exact density in the quadratic corestriction case

Let \(\gamma\) denote the coefficient action on \(E\) of
the same Galois element as \(\beta\). The support incidence
identifies the projection of the permutation of the factors
above a given \(F\)-label with the permutation of \(F\)-labels
induced by \(\gamma\). These actions are by precomposition
on embeddings, so \(\beta\) preserves \(F\) and
\(\beta|_F=\gamma|_F\). Put \(s=\operatorname{ord}(\beta|_F)\).
Every \(F\)-label has orbit length \(s\), since an embedding
is injective. On each two-element fiber of the label map,
\(\beta^s\) acts either as the identity on every fiber
or exchanges the two factors on every fiber. These occur
when \(\operatorname{ord}(\beta)=s\) and \(2s\), respectively.

Blocks over one \(F\)-label have equivalent projective actions
on the connected monodromy, given by the standard tensor
representation of the same two factors. Blocks over distinct
labels have different supports, so their projective spectra
are distinct as functions on the connected group. The
connected projective spectrum field is therefore \(F\).
For example, to distinguish two different supports,
specialize to a torus on a factor in just one support.
The kernel on \(F\) of the action generated by \(\gamma\)
and \(c\), or by the full Galois coefficient image and \(c\),
is exactly \(\{1,c\}\).

On a component of the first kind, each pair of supporting
factors lies in separate orbits. The
[linked Pfaffian theorem](CM_NEWTON_LINKED_RANK_FOUR_PFAFFIAN.md)
gives \(\mu\)-ordinary reduction there away from a set of
natural density zero. On a component of the second kind,
each pair of supporting factors lies in one orbit. The
[tensor obstruction](CM_NEWTON_LINKED_RANK_FOUR_PFAFFIAN.md#a-tensor-obstruction-on-the-other-components)
forces equal middle slopes on every block. At least one
coefficient orbit contains a balanced label. Its averaged
Hodge vector has a strictly positive middle jump, equal to
the fraction of balanced labels in the orbit. Equal middle
slopes therefore exclude \(\mu\)-ordinarity at every good
degree-one prime in this second case. The degree-one primes
have density one. The order condition is conjugacy-invariant,
and finite Chebotarev now gives the asserted density in the
Galois group of \(K\tau(R)/K\).

**Ordinary density when \(F\ne\mathbb Q\).** Let \(L_E/K\)
be the field of definition of the \(E\)-action. Then
ordinary reduction has natural density
\[
                  \frac{1}{[L_EK\tau(R):K]}.
\]
Write \(\theta\) for Frobenius on the embedding labels.
The Hodge-preserving permutation \(\theta\gamma^{-1}\)
fixes the sole balanced \(F\)-label. If \(\gamma|_F\ne 1\),
this label moves under \(\theta\), and its orbit contains
a definite label. Each \(E\)-orbit above this \(F\)-orbit
visits a balanced and a definite label. Averaging their
Hodge vectors gives at least one fractional reference slope.
The ordinary slope vector at the same endpoint would
have a strictly smaller prefix, contrary to Newton above
Hodge. Hence ordinary reduction requires \(\gamma|_F=1\),
so \(\gamma=1\) or \(c\).

When \(\gamma=1\), the Hodge multiplicities are
constant on each coefficient orbit, so all reference slopes
are zero or one. When \(\gamma=c\), each orbit that contains
a definite label has only definite labels, with the two
orientations equally frequent. Its reference slopes are
all \(1/2\). Here one uses that \(\theta\gamma^{-1}\)
preserves Hodge multiplicities and \(\gamma=c\) commutes
with every coefficient label permutation. There is a
definite label because \(F\ne\mathbb Q\).
The fractional slopes rule out ordinarity by the same
prefix comparison. For \(\gamma=1\), ordinary and
\(\mu\)-ordinary reduction coincide and \(\beta|_F=1\).
The preceding density theorem gives ordinary reduction
precisely when \(\gamma=\beta=1\), up to a set of density
zero. These are the primes split completely in
\(L_EK\tau(R)\), so Chebotarev gives the stated density.

If \(F=\mathbb Q\), every reference vector is
\((0,0,1,1)\). Ordinarity then coincides with
\(\mu\)-ordinarity.

For example, if \(K\) has a specified embedding into
\(\mathbb R\), \(R\) is totally real, \(F\ne\mathbb Q\),
and \([K\tau(R):K]=2\), assume the coefficient action
is contained in \(\{1,c\}\). This gives conditional
densities \(1/2\) for \(\mu\)-ordinary and \(1/4\) for
ordinary reduction.
Complex conjugation at the real place acts nontrivially on
\(E\): a definite pair has complex structure a scalar
from \(E\otimes\mathbb R\) on its real block, whereas
complex conjugation reverses the complex structure.
Consequently \(L_E/K\) is a quadratic field that is
nonreal at the chosen place, while \(K\tau(R)\) is real
there. These quadratic extensions are distinct.

## Scope of the field condition

When the derived factors have type \(A_1\), the field
\(\tau(R)\) is the adjoint reflex field of the Hodge
cocharacter. For \(\mathrm{PGL}_2\), conjugacy of a
cocharacter of fixed nonzero weight is determined by
that weight, and the Weyl group exchanges its signs.
The Galois stabilizer of the adjoint conjugacy class
is thus the stabilizer of its single active label.
Every number-field point of a model over this adjoint
reflex field satisfies the no-exchange condition,
provided it has the stipulated derived monodromy.

Moonen's classification of Shimura curves includes
standard tensor actions with this incidence. In his
type IV construction, use an orbit of two disjoint
two-element subsets of the four labels of a totally
real biquadratic field \(R\), with each subset the
fiber over an embedding of a fixed quadratic subfield
\(F\). Its subset field is \(F\). The quaternion algebra
over \(R\) is split at the single active real place,
and the CM extension \(E/F\) can be chosen to split
the corresponding Brauer class by choosing it nonsplit
at every finite ramified place. The partial CM type
consists of one extension of the definite \(F\)-label.
It is primitive: any two distinct \(E\)-labels over
one \(F\)-label are conjugate, and an absolute Galois
element sending one of them to the selected label
sends the other to its distinct conjugate.
The relative rank is
\(2^2=4\), with one balanced and one definite pair.
[Moonen, Examples 3.1–3.2 and Construction 5.6](https://doi.org/10.4153/S0008414X24001159).

**A finite-place obstruction to descent.** Suppose the
Betti adjoint group is
\(\operatorname{Res}_{R/\mathbb Q}\mathrm{PGL}_1(B)\)
for a quaternion algebra \(B/R\). At any auxiliary
\(\ell\) where the derived comparison holds, the
image of \(\beta\) preserves the set of ramified places
of \(B\) above \(\ell\). Indeed the
\(\mathbb Q_\ell\)-Lie algebra of the derived group
decomposes over the factors of its centroid as
\[
 \bigoplus_{w\mid\ell}
       \operatorname{Lie}(\mathrm{PGL}_1(B\otimes_R R_w)).
\]
The monodromy normalizer supplies a semilinear isomorphism
between the summands at \(w\) and \(\beta(w)\). A
split and a division quaternion algebra give different
forms of the rank-one Lie algebra, as is seen from
isotropy of the Killing form. Thus \(\beta\) preserves the
ramification pattern at the chosen \(\ell\).

For a concrete case with \(R/\mathbb Q\) biquadratic, take
a rational \(\ell\) fully split in \(R\), and choose \(B\)
ramified at the three nondistinguished real places and
at exactly one prime above \(\ell\). This ramification set
is permitted by quaternion Brauer reciprocity. Its
\(\ell\)-part has trivial stabilizer under
\(\operatorname{Aut}_{\mathbb Q}(R)\). One can take
\(R=\mathbb Q(\sqrt2,\sqrt3)\), \(F=\mathbb Q(\sqrt2)\),
\(\ell=23\), and \(E=F(i)\). Both \(2\) and \(3\)
are square residues modulo \(23\), while \(-1\)
is not. Thus \(E\) splits the Brauer class of
\(\operatorname{Cor}_{R/F}(B)\) at both the real and
finite places. The resulting
Moonen construction has the required rank-four tensor
incidence. Vasiu's rank-one case supplies the derived
comparison at \(\ell=23\), so any abelian variety with
these Betti data has \(R\subset K\). The \(\mu\)-ordinary
density is one even without assuming in advance that
the field of definition contains the adjoint reflex field.

The conditional densities for \(K\) below the adjoint
reflex field require a fixed abelian variety with this
descent. A specialization over a field containing the reflex field
does not realize the exchanging component. The existing
explicit genus-four examples over \(\mathbb Q\) with
Mumford's three \(A_1\) factors have a tensor of dimension
eight and a different coefficient representation.
[Bouchet–Hanselman–Pieper–Schiavone, Theorems 1.1–1.2](https://arxiv.org/html/2510.00093v2).

## Derived comparison for rank-one factors

Vasiu's Mumford–Tate theorem includes abelian varieties
whose adjoint simple factors are of type \(A_1\), by the
non-special \(A_n\) case with \(n+1=2\). Thus the derived
comparison used above in the quadratic corestriction case
is unconditional at every \(\ell\), including when the
chosen coefficient field is part of a larger geometric
endomorphism algebra.
[Vasiu, Theorem 1.3.4(a) and Example 1.3.5(iv)](https://arxiv.org/pdf/math/0212066v7).

There is also a direct argument for the required derived
equality when the Hodge cocharacter acts nontrivially on
exactly one geometric adjoint factor and all the factors
have type \(A_1\). This argument uses no coefficient-field
assumption. The rational-quotient argument above gives
that the centroid \(R\) is a field.

Fix a rational prime \(\ell\) and use Betti–étale comparison
to view the connected \(\ell\)-adic monodromy \(G_\ell^\circ\)
inside \(M_{\mathbb Q_\ell}\). The inclusion follows from
absolute Hodge tensors, and \(G_\ell^\circ\) is reductive by
Faltings's semisimplicity theorem. After a finite extension
that defines all geometric endomorphisms and makes the
monodromy connected, Faltings's endomorphism theorem and
the identification of rational Hodge endomorphisms of
\(H^1\) with geometric endomorphisms give
\[
 \operatorname{End}_{G_\ell^\circ}(H^1_\ell)
 =
 \operatorname{End}_{M_{\mathbb Q_\ell}}(H^1_\ell).
\]
It follows that the center of \(G_\ell^\circ\) centralizes
\(M_{\mathbb Q_\ell}\), and thus lies in its center.
Projection of \(\operatorname{Lie}(G_\ell^\circ)\) to
\(\mathfrak d_B\otimes\mathbb Q_\ell\) in the central-plus-derived
decomposition is therefore precisely the Lie algebra of
\((G_\ell^\circ)^{\mathrm{der}}\).
[Faltings, Satz 3–4](https://doi.org/10.1007/BF01388432).

Let \(L\subset\overline{\mathbb Q}\subset\mathbb C\) be a
finite extension of the base that contains \(\tau(R)\) and
defines the projector \(P\) and a basis of \(R\) as absolute
Hodge tensors. Such an extension exists before any derived
comparison. For each embedding
\(i:L\hookrightarrow\overline{\mathbb Q}_\ell\), extend it
to \(\overline{\mathbb Q}\) and compare the realizations at
the induced decomposition group. The filtered
\(B_{\mathrm{dR}}\) comparison preserves \(P\) and the
centroid tensors by the de Rham property of absolute Hodge
classes. In particular it identifies the Hodge–Tate weights
on the simple factor labelled \(\rho:R\to\mathbb C_\ell\)
with the degrees of the transported Hodge filtration on
that factor. Those degrees are all zero unless
\(\rho=i\circ\tau\). For \(\rho=i\circ\tau\) there is a
nonzero degree. Equivalently, the local Hodge–Tate
cocharacter has adjoint action just on \(i\circ\tau\).
Here the embeddings \(i\) vary while the Betti tensors
and their étale realizations on the geometric fiber are
fixed.
[Blasius, Theorem (0.3)](https://archive.mpim-bonn.mpg.de/id/eprint/3859/1/preprint_1992_59.pdf),
[Noot, Proposition 2.2, using Wintenberger's Proposition 7](https://irma.math.unistra.fr/~noot/publications/jpaa.pdf).

The local Sen operator on \(H^1_\ell\) belongs to
\(\operatorname{Lie}(G_\ell^\circ)_{\mathbb C_\ell}\).
Its action on the derived Lie algebra is the adjoint
action of its derived projection, so the projection is
nonzero just on the factor \(I\) labelled \(i\circ\tau\).
Write \(\mathfrak h=\operatorname{Lie}((G_\ell^\circ)^{\mathrm{der}})_{\mathbb C_\ell}\).
The image of \(\mathfrak h\) in \(I\simeq\mathfrak{sl}_2\)
is nonzero and semisimple, hence is all of \(I\). The
intersection \(\mathfrak h\cap I\) is a nonzero ideal in
that image: bracketing an element of the intersection
with a lift from \(\mathfrak h\) agrees with bracketing it
with the corresponding element of \(I\). Simplicity
gives \(I\subset\mathfrak h\). Every embedding of \(R\)
into \(\mathbb C_\ell\) arises as \(i\circ\tau\), so
\(\mathfrak h\) contains all factors. This is the derived
comparison at the chosen \(\ell\).
[Sen, Theorem 1](https://annals.math.princeton.edu/1973/97-1/p06).
