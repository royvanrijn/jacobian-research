# The colored Cox-packet frontier for the corrected \(F_{20}\) cover

## Status and outcome

This note constructs the first genuinely colored **combinatorial** Cox
architecture for the corrected Lecacheux \(F_{20}\) boundary.  It proves a
general packet-Hilbert-basis theorem and applies it to all sixty-three exact
colors.

Six primitive packets break the six proportional-row obstructions, but they
do not contain the whole derivative divisor.  Closing under the complete
positive derivative support gives sixteen Galois-stable packets.  Their
indicator columns contain the derivative target integrally and in the
nonnegative semigroup, with one unique normalized model.  Those primitive
columns admit a compact three-divisor allocation

\[
 \operatorname{div}(P_X)=3D_d+D_q+D_r,              \tag{0.1}
\]

whose unique nonnegative exponent vector is \((3,1,1)\).  All three columns
break the required row classes collectively.  The discrete \(q\)-conductor
lattice also passes.

This is not yet a global Cox construction.  The current fan certifies the
base rays and the finite color atlas, not a global regular SNC resolution of
the source.  The global normalization and its total conductor parity are now
exact, but conductor valuations do not determine the residue cocycle of the
individual Cox sections.  The pipeline therefore stops before entrywise
inverse-adjugate polynomiality and affine-space recognition.

The exact certificate is
[`verify_f20_colored_cox_packets.py`](../scripts/verify_f20_colored_cox_packets.py),
with generated output
[`f20_colored_cox_packets.json`](../artifacts/generated-results/f20_colored_cox_packets.json).

## 1. Packet Cartier divisors

Let \(\widetilde X\) be a regular resolved model and let
\(C_1,\ldots,C_m\) be distinct boundary prime divisors.  Suppose a finite
arithmetic or conductor symmetry partitions a selected set of colors into
disjoint packets

\[
 \mathcal O_1,\ldots,\mathcal O_p.
\]

Put

\[
 D_j=\sum_{C_i\in\mathcal O_j}C_i.                 \tag{1.1}
\]

> **Regular packet-Cartier lemma.**  Every \(D_j\) is an effective Cartier
> divisor on \(\widetilde X\).  If the packet is stable under the descent
> symmetry, its divisor class is invariant.  Invariance of the divisor does
> not by itself descend a chosen section: that additionally requires a
> line-bundle and conductor-residue cocycle.

The Cartier assertion is the standard fact that a pure codimension-one
closed subscheme on a regular scheme is Cartier; see
[Stacks Project, Lemma 31.16.9](https://stacks.math.columbia.edu/tag/0BXH).
Equivalently, on a smooth toric chart every invariant Weil divisor is
Cartier, with its integral support function providing the local equation;
see Cox--Little--Schenck, *Toric Varieties*, Section 4.2.

This lemma is conditional in the present \(F_{20}\) application.  The exact
residual calculations prove the listed generic colors and exceptional
profiles, but the repository does not yet contain a global regular SNC
source model realizing all their intersections.  Smoothness of the six-ray
**base** fan is not a substitute: several source colors lie over one base
ray.

## 2. The packet Hilbert-basis theorem

Let \(e_j\in\mathbf N^m\) be the indicator vector of
\(\mathcal O_j\).  Because the packets are nonempty and disjoint, the matrix

\[
 E=(e_1\ \cdots\ e_p)                               \tag{2.1}
\]

has rank \(p\), and its columns are the Hilbert basis of the monoid of
effective divisors which are supported on the packets and constant on every
packet.

> **Packet Hilbert-basis theorem.**  A target
> \(\tau\in\mathbf N^m\) supported on the packet union lies in
> \(E\mathbf N^p\) if and only if it is constant on every packet.  When it
> does, the expression is unique:
> \[
>  \tau=\sum_{j=1}^p \tau|_{\mathcal O_j}\,e_j.      \tag{2.2}
> \]
> The same expression proves integral-span membership.

This is elementary but useful: it distinguishes “we found columns breaking
the old witnesses” from “the full target actually lies in their semigroup.”

There is also an immediate affine-class screen.  Append \(r\) independent
base-character columns and \(z\) zero-order selectors.  If the resulting
\(m\)-row valuation map has rank \(p+r\), contains a unimodular maximal
minor, and is the complete boundary map of a factorial core, then

\[
 \operatorname{rank}\mathcal O(U)^*/k^*=z,
 \qquad
 \operatorname{rank}\operatorname{Cl}(U)=m-p-r.     \tag{2.3}
\]

Thus packet compression can solve a derivative divisor and still be far
from affine space.  Formula (2.3) is conditional on the same complete-core
hypotheses as the boundary-package class-group screen.

## 3. Valuations do not certify conductor gluing

The conductor gate has genuinely nonlinear residue data.  For the nodal
ring

\[
 A=k[x,y]/(xy),
 \qquad
 \widetilde A=k[x]\oplus k[y],                       \tag{3.1}
\]

the two normalization sections \((1,1)\) and \((1,\lambda)\),
\(\lambda\ne1\), have identical order zero on both branches.  The first
matches at the conductor point, while the second does not.  No valuation
matrix can distinguish them.

> **Valuation--residue separation theorem.**  Equality of branch orders is
> necessary for conductor descent of a section, but is not sufficient.
> Even a saturated or unimodular puncture-unit lattice only removes the
> integral character obstruction.  A conductor-fiber identification, or an
> exact residue cocycle, remains necessary.

This is why the compiler must not promote the packet semigroup survivor to a
geometric Cox survivor.

## 4. Six breakers are not enough

For the six witnesses in the
[colored divisor-span theorem](COLORED_DIVISOR_SPAN_OBSTRUCTION.md), take the
primitive packet columns

\[
\begin{gathered}
 D_{d,4},\quad D_{q,\mathrm{cross}},\quad D_{r,2},
 \quad D_{\mathrm{triple}\ E_1,4},\quad
 D_{\mathrm{triple}\ E_2,\mathrm{cluster}},\quad
 D_{qr,E_1,A}.                                      \tag{4.1}
\end{gathered}
\]

Each column has a nonzero mismatch on its corresponding proportional-row
relation.  The enlarged matrix has rank nine, however, while adjoining the
derivative target raises the rank to ten.  These six columns repair the six
displayed certificates but leave positive derivative colors at the
\(q\)-node, the \(r\)-cusp, and the remaining \(q\)-\(r\) exceptional
packets uncovered.

This is an exact warning against treating a list of broken witnesses as a
sufficiency theorem.

## 5. The sixteen-packet closure

The complete positive support of \(\nu(P_X)\) splits into the following
sixteen declared rational packets.  The last column is the unique exponent
in (2.2).

| packet | geometric colors | exponent |
|---|---:|---:|
| \(d\)-ramification | 1 | 3 |
| \(q\)-crossing | 2 | 1 |
| \(q\)-node slopes | 4 | 1 |
| \(r\)-ramification | 2 | 1 |
| \(r\)-cusp \(E_1\) | 1 | 4 |
| \(r\)-cusp \(E_2\) | 1 | 8 |
| \(r\)-cusp \(E_3\), unramified | 1 | 2 |
| \(r\)-cusp \(E_3\), ramified | 2 | 4 |
| \(r\)-cusp \(E_4\) | 5 | 4 |
| triple \(E_1\), index four | 2 | 7 |
| triple \(E_2\), cluster | 8 | 3 |
| \(qr\ E_1\), \(A\)-ramified | 3 | 2 |
| \(qr\ E_1\), \(A\)-unramified | 3 | 1 |
| \(qr\ E_1\), \(B\)-ramified | 3 | 1 |
| \(qr\ E_2\), \(A\)-sheets | 9 | 2 |
| \(qr\ E_2\), \(B\)-sheets | 6 | 1 |

Together with the three base-factor columns and the zero finite-order
\(q\)-selector, this gives a \(63\times20\) matrix of rank nineteen.  The
derivative-augmented rank remains nineteen.  A displayed \(19\times19\)
minor has determinant one, and (2.2) gives target class order one.

The nonnegative model is unique after normalizing the selector exponent:
the target-zero rows `d_unramified`, `q_residual_1`, and `r_unramified`
force the three base exponents to zero; every packet pivot then forces the
exponent in the table.

### 5.1 Three different-factor Cartier columns

The primitive packet basis exposes a much smaller candidate.  Define three
effective packet divisors by

\[
\begin{aligned}
D_d={}&D_{d,4}+D_{\mathrm{triple}\ E_1,4}
              +D_{\mathrm{triple}\ E_2,\mathrm{cluster}},\\
D_q={}&D_{q,\mathrm{cross}}+D_{q,\mathrm{node}}
       +2D_{\mathrm{triple}\ E_1,4}
       +D_{qr,E_1,A,\mathrm{ram}}
       +D_{qr,E_1,A,\mathrm{unr}}+D_{qr,E_2,A},\\
D_r={}&D_{r,2}+4D_{rE_1}+8D_{rE_2}+2D_{rE_3,\mathrm{unr}}
       +4D_{rE_3,\mathrm{ram}}+4D_{rE_4}\\
    &+2D_{\mathrm{triple}\ E_1,4}
       +D_{qr,E_1,A,\mathrm{ram}}+D_{qr,E_1,B,\mathrm{ram}}
       +D_{qr,E_2,A}+D_{qr,E_2,B}.
                                                               \tag{5.1}
\end{aligned}
\]

The allocation follows the exact local different factors: the \(q\)-node
is in \(D_q\), the ramphoid cusp is in \(D_r\), the first triple exceptional
has contribution \(3+2+2=7\), the second triple cluster has contribution
three from \(D_d\), and the \(q\)-\(r\) residual discriminants split as
\(q^2r\) on the \(A\)-cluster and \(r\) on the \(B\)-cluster.

Coefficientwise on all sixty-three rows,

\[
 \boxed{\nu(P_X)=3\nu(D_d)+\nu(D_q)+\nu(D_r).}       \tag{5.2}
\]

The three-column matrix has rank three and the augmented rank remains three.
After adjoining the base columns and the zero-order selector the ranks are
both six.  Exhaustion of \(0\leq a,b,c\leq8\), with the generic rows already
forcing the values, gives the unique model

\[
 (a,b,c)=(3,1,1).                                   \tag{5.3}
\]

On a global regular colored resolution the divisors in (5.1) are Cartier by
Section 1, so they are the first compact Cartier-compatible Cox columns.
What remains conditional is the existence of that global resolution and of
compatible global Cox sections, not the sixty-three-row identity.

The first global algebra test is now complete in
[`F20_GLOBAL_MULTI_REES_COX_ALGEBRA.md`](F20_GLOBAL_MULTI_REES_COX_ALGEBRA.md).
The natural incidence ideals on the nonnormal root hypersurface do **not**
realize these columns: at a triple-(E_1) ramified color their (d)-ideal
has order two, while (D_d) requires order one and the original local value
semigroup has a gap at one.  Thus the Cox construction must be performed
after normalization; this sharpens, rather than removes, the conditional
status above.

The next
[`F20_NORMALIZED_COX_CONDUCTOR_FRONTIER.md`](F20_NORMALIZED_COX_CONDUCTOR_FRONTIER.md)
computes that normalization.  Its two nontrivial global module generators
still have triple-(E_1) order two, so the missing value-one column is an
exceptional variable on a controlled transform, not a normalization
generator.

## 6. Conductor and principal-selector gates

The packet orders are constant on every declared Galois orbit, and the two
\(q\)-crossing orders agree.  The exact conductor pullback matrix still has
Smith diagonal \((1,1,1)\), and adjoining \(w-1\) gives determinant \(-1\).
Hence there is no integral unit-lattice obstruction.

The total conductor equalizer is now explicit.  On the connected rational
cover its slope-residual discriminant is a square, its distinguished square
root is anti-invariant, and the selector completion remains unimodular.  The
[exceptional Cox atlas](F20_EXCEPTIONAL_COX_ATLAS.md) gives one exact
parity-compatible residue factorization for the compact three-column
compression on the punctured conductor.  There is not yet a residue cocycle
for all sixteen individual Cox sections, nor certified overlap gluing of the
three compact frames through the controlled charts, so the global
sectionwise conductor gate remains `uncertified` by Section 3.

For the natural pre-normalization incidence control, the missing residue is
now exact: the degree-((3,1,1)) product does not contain (P_X), its cyclic
residue has length fifty-seven, its base projection has length thirty-three,
and its reduced support is exactly the eight known finite collision centers.
Primary decomposition splits it into the node, cusp, conjugate-triple,
transverse, and tangency packets.  Extending the natural product to the
global normalization still does not contain (P_X).  These are not the
residue cocycles of the still-unconstructed exceptional Cox sections.

The obvious polynomial sheet separators do not avoid Cox data.  Exact norms
give

\[
 \operatorname{Nm}(4X-1)=-4d(t-2)^2,                \tag{6.1}
\]

and the natural linear \(q\)-collision selector has norm \(q^2H(s,t)\) with
nonconstant \(H\).  Both create interior divisors away from \(dqr=0\); they
are not units on the discriminant complement.  This does not exclude more
complicated principal units, but it proves that the packet columns are not
already realized by the two natural linear choices.

## 7. Downstream status

At the divisor level the Cox monomial in (2.2) has the same sixty-three
orders as \(P_X\), so divisorial inverse-adjugate cancellation passes
formally.  Local derivative-denominator polynomiality now passes on all
thirteen controlled chart types, but entrywise polynomiality is not a
divisor-only or chartwise-only statement.  It is not reached without the
global Cox algebra and its overlap cocycle.

The three-column compression fails the independent conditional affine-space
screen (2.3) even more strongly.  With the three base characters and the
zero-order selector its rank is six, so the conditional ranks are

\[
 \operatorname{rank}\mathcal O(U)^*/k^*=1,
 \qquad
 \operatorname{rank}\operatorname{Cl}(U)=63-6=57.  \tag{7.1}
\]

Expanding to the complete sixteen-packet basis raises the rank to nineteen
but still gives

\[
 \operatorname{rank}\mathcal O(U)^*/k^*=1,
 \qquad
 \operatorname{rank}\operatorname{Cl}(U)=63-19=44. \tag{7.2}
\]

This conclusion applies only if the twenty columns are the complete unit
basis of a certified factorial core.  It says that the naive packet slice
cannot be the desired affine space under those standard hypotheses.

Adjoining one independent coordinate for every geometric color replaces the
packet matrix by \(I_{63}\) and removes the unit/class obstruction at the
lattice level.  It does not solve Galois and conductor descent, Cox relations
and irrelevant loci, the dimension-preserving quotient, or affine-space
recognition.  That full-color construction is therefore the next honest
frontier, not a realized Keller map.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_f20_colored_cox_packets.py \
  --output artifacts/generated-results/f20_colored_cox_packets.json
```

The checker verifies the six broken relations, the rank-nine/rank-ten
failure of the six-column attempt, the complete sixteen-packet partition,
the compact identity \(3D_d+D_q+D_r=\operatorname{div}(P_X)\) with unique
model \((3,1,1)\), the primitive unique nonnegative model, the unimodular
rank-nineteen minor, the \(q\)-conductor lattice, and the two
principal-selector norm failures.
