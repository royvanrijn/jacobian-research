# Laurent separation and binary Gaussian vanishing

## Laurent constant terms

Let \(F\) be a nonzero Laurent polynomial over a characteristic-zero field.
The following conditions are equivalent.

- \(0\notin\operatorname{Newt}(F)\).
- \(\operatorname{CT}(F^m)=0\) for every positive integer \(m\).
- \(\operatorname{CT}(F^m)=0\) for every sufficiently large integer \(m\).
- For every fixed Laurent polynomial marker \(Q\), \(\operatorname{CT}(QF^m)=0\)
  for every sufficiently large integer \(m\).

Duistermaat and van der Kallen prove the implication from vanishing at every
positive power to Newton separation over \(\mathbb C\). Their arbitrary-rank
argument first reduces to the minimal face containing zero and uses the
compactification and residue analysis in Sections 3–5. The same implication
holds over every characteristic-zero field by embedding the finitely
generated coefficient subfield into \(\mathbb C\). The embedding preserves
each nonzero coefficient and each constant-term identity. See their
[primary paper](https://webspace.science.uu.nl/~kalle101/powers.pdf) and the
[field reduction in the binary GVC paper](../../../results/generalized-vanishing-two-variables/main.tex).

If the constant terms vanish for every \(m\ge m_0\), choose \(q\ge\max(1,m_0)\)
and apply the positive-power theorem to \(F^q\). All of its positive powers
have zero constant term. The identity
\(\operatorname{Newt}(F^q)=q\operatorname{Newt}(F)\) follows because each vertex
is exposed by a form for which its term is the unique maximizer, so its
coefficient survives taking powers. Thus \(0\notin\operatorname{Newt}(F)\).
The same argument needs only vanishing on all positive multiples of a single
fixed \(q\).

For the converse and the marker bound, choose a real linear form \(\ell\)
with \(\delta=\min_{\operatorname{supp}(F)}\ell>0\). For each nonzero \(Q\) put
\(b=-\min_{\operatorname{supp}(Q)}\ell\). Every exponent in \(QF^m\) then has
\(\ell\)-value at least \(m\delta-b\), so its constant term vanishes when
\(m\delta>b\). This gives every positive pure power by taking \(Q=1\).
Eventual vanishing for all markers implies eventual pure vanishing by the
same choice.

## One radial variable and an arbitrary angular torus

In the [pinned manuscript of Christopher D. Long](https://github.com/octonion/mathematics/blob/2f1ba8ae5d3d282141aca2285326a2d31ce20cd2/gmc/gmc2_stronger_arbitrary_torus.tex),
the functional on \(\mathbb C[u,z_1^{\pm1},\ldots,z_r^{\pm1}]\) sends
\(u^jz^\mu\) to \(j!\) for \(\mu=0\), and to zero otherwise. In this setting
the same four conditions are equivalent with \(\operatorname{CT}\) replaced
by this weighted functional and the Newton polytope replaced by the convex
hull of the angular exponents alone. The same marker bound from separation
holds for the angular weights. This statement uses exactly one radial
variable.

Long's stronger nonvanishing assertion gives, when the origin belongs to
the angular convex hull, one integer \(k\ge1\) for which the weighted value
of \(F^{kp}\) is nonzero for every sufficiently large prime \(p\). The proof
uses an exposed face minimizing radial degree at angular weight zero. The
angular projection of that face is injective on its support, so there is
no coefficient cancellation before applying Duistermaat–van der Kallen.
A fixed power of the face polynomial yields a nonzero balanced coefficient
at the minimal radial degree \(d\). Dividing the prime-dilated factorial
expression by \((pd)!\), then reducing in characteristic \(p\), preserves
that coefficient to the \(p\)-th power and kills the higher factorial
ratios.

The finite-type coefficient domain is localized at the distinguished
coefficient. By [Chevalley's theorem](https://stacks.math.columbia.edu/tag/00FE)
the image of its spectrum in \(\operatorname{Spec}\mathbb Z\) is
constructible and contains the generic point. It therefore contains all
but finitely many primes. The fibre over each such prime is a nonzero finite-type \(\mathbb F_p\)-algebra.
Taking a maximal ideal in that fibre gives a finite residue field by
Zariski's lemma and preserves the inverted coefficient. This proves the mathematical good-reduction
statement used here and in the written binary argument. The corresponding
Lean interface still requires a formal proof.

## Two standard Gaussian variables

For \(Z=(X+iY)/\sqrt2\), \(W=(X-iY)/\sqrt2\) and a complex polynomial \(P\),
the substitution \(Z\mapsto z,\ W\mapsto uz^{-1}\) is injective and
intertwines Gaussian expectation with Long's functional at \(r=1\).
Indeed \(Z^aW^b\) maps to \(u^b z^{a-b}\), with expected value
\(\delta_{ab}a!\). The four equivalent binary conditions are one strict
sign of every rotational weight \(a-b\), vanishing of all positive pure
moments, eventual vanishing of the pure moments, and eventual vanishing
for each polynomial marker. If the support has weight zero in its
convex hull, Long's prime-dilation nonvanishing also transfers. The zero
phase satisfies the vanishing conditions directly.

There is a second route from eventual pure vanishing to one strict sign,
using only the existing positive-power binary classification. Choose
\(q\) beyond the onset and apply that classification to \(P^q\).
The lowest and highest rotational weight parts of \(P\), raised to the
\(q\)-th power, are exactly the corresponding extreme weight parts of
\(P^q\) and are nonzero. One strict sign for \(P^q\) therefore implies the
same strict sign for \(P\). This route applies to the
[written lower-face theorem](../../foundations/extended-geometry/TWO_REAL_GMC_LOWER_FACE_THEOREM.md)
or to [Wilson's Theorem 1](https://arxiv.org/html/2607.23887v1), with the
[repository's binary Gaussian paper](../../../results/gaussian-moments-two-variables/main.tex)
giving a separately written filtered Frobenius proof.

If \(\sigma w(P)>0\) on the support for \(\sigma=1\) or \(-1\), set
\(\delta=\min\sigma w(P)\) and \(b=-\min\sigma w(Q)\) for a nonzero marker
\(Q\). Then \(\mathbb E(QP^m)=0\) whenever \(m\delta>b\). In particular
\(\delta\ge1\) and \(b\le\deg Q\) recover \(m>\deg Q\).
