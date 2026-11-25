# Simultaneous lifts of split one-failed-pair reductions

Fix \(F=\mathbb Q(\zeta_7+\zeta_7^{-1})\) and use auxiliary level
\(3\) in the [finite-field trace construction](CM_NEWTON_SPLIT_POLYNOMIAL_OBSTRUCTION.md).
It supplies absolutely simple \(\mathcal O_F\)-abelian threefolds
\(B_p/\mathbb F_p\) with rational \(3\)-torsion. Their \(p\)-divisible
groups have slopes \((0,1)\) at two coefficient places and
\((1/2,1/2)\) at the remaining one. The rational primes in the construction
split in \(F(i)\) and satisfy \(p\equiv1\bmod 9\).

**Proposition.** Given a finite set of distinct primes with chosen models
\(B_p\) of this kind, there is a number field \(K\) in which each of those
rational primes splits completely, and an abelian threefold \(B/K\)
with an \(\mathcal O_F\)-action, rational \(3\)-torsion and
\(\operatorname{End}^0(B_{\overline K})=F\). At every place of \(K\)
over a specified \(p\), it has good reduction isomorphic to the given
\(B_p/\mathbb F_p\) with its \(\mathcal O_F\)-action. In particular its
three rank-two projective monodromy factors are independent. One can
require good reduction away from a fixed finite set of rational primes
and require \(K\) to be linearly disjoint over \(K_0=\mathbb Q(\zeta_3)\)
from any prescribed finite extension of \(K_0\). The number field and
the threefold may depend on the finite set of prescriptions.

## A common fine moduli space

Each of the three coefficient summands of \(B_p[p^\infty]\) has height
two. Indeed \(a=\pi_p+p/\pi_p\) generates \(F\), and on the rational
Dieudonné module each of the three eigenvalues of \(a\) occurs on its
corresponding quadratic Frobenius pair. The dimension of each summand
is one from its height and Newton slopes, so the Lie algebra is free of
rank one over \(\mathcal O_F\otimes\mathbb F_p\). This is the Rapoport
condition. In particular the polarization module \(\mathcal P(B_p)\)
gives a polarization in the sense of Deligne–Pappas. It is a constant
projective \(\mathcal O_F\)-module of rank one with the positive cone
defined by polarizations. [Yu, Proposition 2.8](https://www.numdam.org/article/AIF_2003__53_7_2105_0.pdf)
· [Rapoport, Proposition 1.17 and Variant 1.18](https://www.numdam.org/item/CM_1978__36_3_255_0.pdf).

The class group of \(F\) is trivial by the Minkowski bound \(14/9\),
using discriminant \(49\). The three conjugates of
\(\zeta_7+\zeta_7^{-1}\) are units with signatures obtained by permuting
one plus and two minuses. Together with \(-1\) these give all unit
signatures. The narrow class group is therefore trivial, so we can fix
a positive isomorphism of \(\mathcal P(B_p)\) with the inverse different
\(\mathfrak d_F^{-1}\) at every prescribed prime. The chosen module has
the usual totally positive cone. [Rapoport, §1.15 and Theorem 1.20](https://www.numdam.org/item/CM_1978__36_3_255_0.pdf).

Take the fine Hilbert moduli scheme with this fixed polarization module
and determinant-fixed full level \(3\) over
\(\mathcal O_{K_0}[1/21]\), where \(K_0=\mathbb Q(\zeta_3)\).
It is smooth, with geometrically connected generic fibre and hence
geometrically irreducible generic fibre. Smoothness and representability
follow from the fine level and the invertibility of the field
discriminant. The choice of \(\zeta_3\) fixes the determinant component
in the full-level space. [Rapoport, Theorems 1.20, 1.22 and 1.28 and Lemma 1.23](https://www.numdam.org/item/CM_1978__36_3_255_0.pdf)
· [Deligne–Pappas, Theorem 2.2 and §§2.4, 2.10, 2.12–2.13](https://www.numdam.org/item/CM_1994__90_1_59_0.pdf).

The prime \(3\) is inert in \(F\), so \(B_p[3]\) is free of rank two
as an \(\mathcal O_F/3\)-module, by its rank-six \(\mathbb F_3\)-dimension.
Its Weil pairing with the chosen polarization module gives the
perfect alternating determinant pairing of the level moduli problem.
Because the torsion is rational, one can choose a rational basis and
rescale one vector by a unit of \(\mathcal O_F/3\) to match the prescribed
determinant. Every \(p\) under consideration splits in \(K_0\).
For each place of \(K_0\) over \(p\), use the reduction of its chosen
\(\zeta_3\) to fix this determinant, so the model gives a rational point
on the special fibre of the same scheme at each such place.
[Rapoport, §§1.21–1.22](https://www.numdam.org/item/CM_1978__36_3_255_0.pdf).

## Prescribing the reductions over one number field

Smoothness gives a lift of each special-fibre point to the corresponding
complete discrete valuation ring \(\mathbb Z_p\). The subset of local
points with that specified integral reduction is nonempty and open in
the \(p\)-adic topology. Apply Moret–Bailly's theorem to the generic
fibre over \(K_0\) with these finitely many local conditions. It gives an
extension \(K/K_0\) split completely at all the specified local places
and a point whose images in their completions lie in the specified open
sets. The moduli scheme is fine, so the point defines the asserted
abelian variety with the specified reductions. As the conditions were
imposed at both places of \(K_0\) above each \(p\), these rational primes
split completely in \(K\). [Moret–Bailly, Theorem 1.3 and Remarks 1.5 and 1.7](https://www.numdam.org/item/10.24033/asens.1582.pdf).

It remains to guarantee the exact geometric endomorphism field. Add
two auxiliary primes of the same kind and their corresponding models
to the finite list before applying the approximation. Write
\(E_j=F(\pi_{p_j})\) for their quadratic CM fields over \(F\). The
one-failed-pair ramification calculation shows absolute simplicity and
shows that the powers of the six Frobenius roots stay distinct. The
simple Honda–Tate isogeny class over each finite residue extension
then still has characteristic polynomial of degree six and exponent
one. Its rational endomorphism algebra is the same field \(E_j\), so
this is also the geometric endomorphism algebra. The second auxiliary
prime can be chosen unramified in \(E_1\). The field \(E_2\) ramifies
above that prime by the valuation-one discriminant at its failed pair.
Thus \(E_1\) and \(E_2\) are nonisomorphic over \(F\).
[Frobenius ramification and absolute simplicity](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group)
· [Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/).

Specialization of geometric endomorphisms at each good reduction is
injective and respects the given action of \(F\). Thus the geometric
endomorphism algebra of \(B\) is a field containing \(F\) that embeds
over \(F\) in both quadratic extensions \(E_1,E_2\). It equals \(F\).
Lombardo's image theorem at sufficiently large split auxiliary primes
now gives the full common-determinant group, with derived algebraic
group \(\mathrm{SL}_2^3\). This verifies the monodromy assertion.
[Yu, Lemma 2.11](https://www.numdam.org/article/AIF_2003__53_7_2105_0.pdf)
· [Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).

## Controlling the base field and good reduction

For the integral form of Moret–Bailly's theorem, take the moduli scheme
over the ring of integers of \(K_0\) with a fixed finite set of rational
primes and all prescribed local places inverted. The fixed set can be
chosen once to make the smooth model surjective, since its image is
open and contains the generic point. Each component of the smooth
model meets the generic fibre, so the model is irreducible. The local
opens at the prescribed
places specify good reduction there. The integral point supplied by the
theorem gives good reduction at all remaining places outside the fixed
set. [Moret–Bailly, Theorem 1.3 and Remark 1.5](https://www.numdam.org/item/10.24033/asens.1582.pdf).

For linear disjointness take a finite Galois extension \(L/K_0\)
containing the prescribed extension. Choose one auxiliary place in
each conjugacy class of \(\operatorname{Gal}(L/K_0)\), outside the
places used so far. At sufficiently large places the fine moduli
scheme has geometrically integral reduction by spreading out, and
Lang–Weil supplies rational points. Add their smooth local reduction
tubes to the approximation. Complete splitting of these places in
\(K\) means each selected Frobenius acts trivially on the embeddings
of \(K\cap L\) into \(L\). Since every conjugacy class was selected,
this action is trivial, so \(K\cap L=K_0\). As \(L/K_0\) is Galois,
this gives linear disjointness. The added tubes also prescribe good
reduction at their places. [Lang and Weil, Theorem 1](https://doi.org/10.2307/2372655).

The same construction gives an all-balanced rank-two CM action on
\(B^2\) by adjoining \(i\) as the rotation of its two real-multiplication
factors. The global construction permits any finite set of these exact
one-failed-pair reductions on one characteristic-zero variety with
independent projective factors. It gives no density conclusion for the
split primes of a fixed threefold.
