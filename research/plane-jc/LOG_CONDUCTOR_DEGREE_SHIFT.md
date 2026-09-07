# Log-conductor degree shift and determinant insufficiency

> **Status.**  This note proves the natural cohomological placement of a
> conductor mismatch, proves that the normalized logarithmic determinant has
> zero conductor mismatch on every resolved boundary, and gives two genuine
> local polynomial Jacobian matrices with the same determinant divisor and
> generic branch profiles but different nodal cokernels.  It corrects the
> proposed degree-zero comparison in
> [`UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md`](UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md).
> The result does **not** prove `JC(2)`: it shows that a terminal type-I
> determinant is insufficient and identifies the missing invariant as the
> nodal `Fitt_1`/localized-second-Chern profile of the full logarithmic
> differential.

The exact local models and the `(75,125)` terminal normalization are checked
by
[`verify_log_conductor_degree_shift.py`](../scripts/verify_log_conductor_degree_shift.py).

## 1. Two bases and two cohomological degrees

The previous programme combined objects living on two different spaces.
The conductor matching module in the wild-boundary atlas is a coherent
module over a **coefficient base** `S`.  The logarithmic cotangent cokernel

\[
 \mathcal T_f^{\log}=\operatorname{coker}\left(
 f^*\Omega_Y^1(\log D_Y)\longrightarrow
 \Omega_X^1(\log D_X)
 \right)                                           \tag{1.1}
\]

is a sheaf on a **resolved source surface** `X`.  A comparison between them
requires a resolved family over `S` and a derived pushforward; it cannot be
an unqualified inclusion of one module into the other.

There is a second distinction.  A section of a pure curve module which is
supported at finitely many points belongs to `H_Z^0`.  A mismatch between
regular sections on the normalized branches is a principal part and belongs
naturally to `H_Z^1`.  The normalization sequence makes this degree shift
exact.

## 2. The conductor degree-shift theorem

Let `C` be a reduced Noetherian curve with finite normalization

\[
 \nu:\widetilde C\longrightarrow C,
\]

and let `N` be a coherent torsion-free `O_C`-module.  Put

\[
 \widetilde N=
 \nu_*\bigl(\nu^*N/\text{torsion}\bigr),
 \qquad
 \mathcal C_N=\widetilde N/N.                    \tag{2.1}
\]

The first map in (2.1) is injective because its kernel has finite support and
`N` is torsion-free.  Let `Z` be any finite closed set containing the support
of `C_N`.

### Theorem 2.1 -- normalization mismatches live in degree one

There is a functorial short exact sequence

\[
 \boxed{
 0\longrightarrow\mathcal C_N
 \mathop{\longrightarrow}^{\delta}
 \mathcal H_Z^1(N)
 \longrightarrow
 \mathcal H_Z^1(\widetilde N)
 \longrightarrow0.}                              \tag{2.2}
\]

Moreover,

\[
 \mathcal H_Z^0(N)=\mathcal H_Z^0(\widetilde N)=0,
 \qquad
 \operatorname{Hom}_{\mathcal O_C}(\mathcal C_N,N)=0. \tag{2.3}
\]

Consequently a nonzero conductor mismatch has no nonzero
`O_C`-linear lift to `H_Z^0(N)`.  If `s_tilde` is a regular normalized-branch
section, its class

\[
 [s_{\rm tilde}]\in\mathcal C_N                 \tag{2.4}
\]

vanishes exactly when it descends to a regular section of `N`.

#### Proof

Torsion-freeness excludes finite-support submodules of both `N` and
`N_tilde`, proving the two degree-zero vanishings.  The quotient `C_N` is
supported on `Z`, so

\[
 \mathcal H_Z^0(\mathcal C_N)=\mathcal C_N,
 \qquad
 \mathcal H_Z^i(\mathcal C_N)=0\quad(i>0).
\]

Apply local cohomology to

\[
 0\longrightarrow N\longrightarrow\widetilde N
 \longrightarrow\mathcal C_N\longrightarrow0.
\]

The resulting long exact sequence gives (2.2).  The image of a map from the
finite-support module `C_N` to `N` would be a finite-support submodule of
`N`, hence is zero.  Finally (2.4) is zero precisely when its representative
belongs to the image of `N`.  \(\square\)

On an affine neighbourhood, this is the sheaf version of the standard
localization exact sequence

\[
0\to H_Z^0(N)\to\Gamma(C,N)\to
\Gamma(C-Z,N|_{C-Z})\to H_Z^1(N)\to0.
\tag{2.5}
\]

See the Stacks Project,
[*Local Cohomology*, Lemma 8.2](https://stacks.math.columbia.edu/tag/0BK0).

## 3. The split-node control

Let

\[
 A=k[u,v]/(uv),\qquad
 \widetilde A=k[u]\oplus k[v].                  \tag{3.1}
\]

The normalization embeds `A` as the pairs with equal constant terms, so

\[
 0\longrightarrow A\longrightarrow\widetilde A
 \mathop{\longrightarrow}^{(f,g)\mapsto f(0)-g(0)}
 k\longrightarrow0.                              \tag{3.2}
\]

The class of `(1,0)` is nonzero in the conductor quotient and is killed by
the node ideal `(u,v)`.  It does not define a section in `H^0_(u,v)(A)`,
which is zero.  Instead, (2.2) sends it injectively to

\[
 H^1_{(u,v)}(A)
 =\frac{A_u\oplus A_v}{A}.                       \tag{3.3}
\]

Thus even the simplest nodal mismatch is a degree-one principal part.
`S1` gives uniqueness of descent across the puncture, not existence of a
descent.

## 4. Normalized logarithmic determinant descent

Let

\[
 f:(X,D_X)\longrightarrow(Y,D_Y)                 \tag{4.1}
\]

be a morphism of smooth integral SNC surface pairs whose logarithmic
differential is generically invertible.  Its determinant is a nonzero
section

\[
 j_f\in H^0(X,\mathcal L_f),\qquad
 \mathcal L_f=
 \det\Omega_X^1(\log D_X)\otimes
 f^*\det\Omega_Y^1(\log D_Y)^{-1}.               \tag{4.2}
\]

Write `Delta_f=div(j_f)`.

### Theorem 4.1 -- complete-different conductor closure

The section `j_f` gives a canonical isomorphism

\[
 \mathcal L_f(-\Delta_f)\simeq\mathcal O_X.       \tag{4.3}
\]

For every reduced boundary subcurve `B` and its normalization, the
normalized determinant `j_f/s_Delta` has zero conductor mismatch.  At an SNC
node with local equations `uv=0`, write

\[
 j_f=u^a v^b h(u,v)\,e,\qquad h(0,0)\ne0,        \tag{4.4}
\]

in a local frame `e`.  After the complete divisor has been removed, the two
branch restrictions have the same conductor value `h(0,0)`.

This remains true after every boundary blowup, provided the new logarithmic
discrepancy is included in the transformed different.  Blowing up a node is
log crepant.  Blowing up a smooth boundary point contributes the exceptional
divisor once to `K_X+D_X`, and the logarithmic Jacobian acquires exactly that
exceptional factor.

#### Proof

A nonzero section of a line bundle identifies that line bundle with the
Cartier divisor of the section, which proves (4.3).  Restriction of the
resulting global unit to the two normalized branches agrees in their common
fiber, proving (4.4).  For a boundary blowup `pi`, determinants in the
composition of logarithmic differentials multiply.  Hence

\[
 \mathcal L_{f\circ\pi}
 =\pi^*\mathcal L_f\otimes\mathcal L_\pi,
 \qquad
 j_{f\circ\pi}=\pi^*j_f\,j_\pi,                  \tag{4.5}
\]

and the divisor of the right side is the complete transformed different.
Applying the first assertion again proves the claim.  \(\square\)

The theorem closes the scalar determinant comparison: once the divisor is
complete, its conductor residue is identically zero.  It cannot, however,
identify the cokernel of the full rank-two logarithmic differential.

## 5. Determinant and generic Smith data are insufficient

Work over a characteristic-zero field and let `R=k[u,v]_(u,v)`.  Consider
the two polynomial maps

\[
 F_{\rm glue}=\left(\frac{u^2v}{2},v\right),
 \qquad
 F_{\rm split}=\left(\frac{u^2}{2},\frac{v^2}{2}\right). \tag{5.1}
\]

Their Jacobian matrices are

\[
 A_{\rm glue}=
 \begin{pmatrix}uv&u^2/2\\0&1\end{pmatrix},
 \qquad
 A_{\rm split}=
 \begin{pmatrix}u&0\\0&v\end{pmatrix}.          \tag{5.2}
\]

Both have determinant `uv`, and after the full determinant divisor is
removed both normalized determinants are `1`.  Along the generic point of
either branch their cokernels have the same rank-one Smith profile.
Nevertheless a regular row operation gives

\[
 \operatorname{coker}A_{\rm glue}\simeq R/(uv),
 \qquad
 \operatorname{coker}A_{\rm split}
 \simeq R/(u)\oplus R/(v).                       \tag{5.3}
\]

The second module is the normalization of the first, and their difference
is the one-dimensional node quotient (3.2).  Equivalently,

\[
 \operatorname{Fitt}_1(\operatorname{coker}A_{\rm glue})=R,
 \qquad
 \operatorname{Fitt}_1(\operatorname{coker}A_{\rm split})=(u,v). \tag{5.4}
\]

### Theorem 5.1 -- determinant-insufficiency theorem

The determinant divisor, its normalized unit, and the generic Smith profile
on every irreducible branch do not determine the nodal conductor gluing of
the cotangent cokernel, even for matrices which are genuine polynomial
Jacobians.  A universal comparison theorem must retain at least the nodal
`Fitt_1` profile, or equivalently the corresponding localized codimension-two
Chern/Smith defect.

The proof is the pair (5.1)--(5.4).  Notice that this is stronger than an
arbitrary-matrix countermodel: both presentations satisfy the integrability
relations of an actual differential.

## 6. Exact effect on the terminal type-I programme

For the selected `(75,125)` terminal block, put

\[
 s=X^{17}y^5,
\]

\[
 P=X^4y(1+s),\qquad
 Q=-X\left(1+3s+\frac95s^2\right).               \tag{6.1}
\]

The exact calculation already recorded in
[`F2_TERMINAL_RESIDUE_COVER.md`](F2_TERMINAL_RESIDUE_COVER.md) is

\[
 [P,Q]_{X,y}=X^4.                                \tag{6.2}
\]

Thus the type-I nonvanishing supplies the determinant divisor `4(X=0)`.
After that complete divisor is removed,

\[
 X^{-4}[P,Q]_{X,y}=1.                            \tag{6.3}
\]

By Theorem 4.1 its scalar conductor mismatch is zero at every terminal
attachment.  Therefore terminal type-I nonvanishing does **not** imply the
nonzero finite-support residue postulated in the previous programme.  The
degree-six `A_6` residue cover describes the map along the terminal divisor;
it does not by itself give a nodal extension defect of the rank-two
logarithmic cokernel.

## 7. Corrected complete-chain closure criterion

For a map-decorated realized chain `tau`, the canonical local object is the
full logarithmic matrix, not its determinant.  Let `R_tau` be the selected
reduced curve in the one-dimensional support and first form its torsion-free
restriction

\[
 \mathcal N_\tau=
 \left(\mathcal T_\tau^{\log}\otimes\mathcal O_{R_\tau}\right)
 /\text{torsion}.
\]

Then form

\[
 \mathcal C_\tau=
 \nu_*\left(\nu^*\mathcal N_\tau/\text{torsion}\right)
 /\mathcal N_\tau.                               \tag{7.1}
\]

Any nilpotent thickening of the determinant support is additional matrix
data and must be retained separately; (7.1) records the reduced branch
gluing only.

Theorem 2.1 places `C_tau` canonically in `H_Z^1`, while Theorem 5.1 shows
why it cannot be recovered from the determinant ledger.  A valid universal
contradiction theorem must prove all of the following.

1. **Matrix realization.**  Each Newton transformation is related to an
   actual boundary blowup chart, and the complete `2 x 2` logarithmic
   differential is transported up to invertible row and column operations
   and the explicit discrepancy factor.
2. **Nodal profile.**  Compute `Fitt_1` and the normalization defect (7.1) at
   every attachment, not only `Fitt_0=(j_f)` at generic components.
3. **Residue identification.**  Show that a distinguished Laurent class is
   the normalization mismatch of an actual branchwise section of (7.1).
4. **Nonzero nodal class.**  Prove that this class is nonzero in `C_tau`.
   Type-I bracket nonvanishing alone does not do this.
5. **Descent contradiction.**  Independently prove that Keller geometry
   makes the same branchwise section descend, equivalently that its class in
   `C_tau` is zero.

Items 3--5 give an immediate contradiction by Theorem 2.1.  This criterion
is degree-independent, but it is not automatic: the split matrix (5.2)
shows that nontrivial nodal behavior is compatible with the same determinant
data.

The next viable invariant is therefore precisely
[`Attack E`](FRONTIER_CLOSING_ATTACKS.md#attack-e--logarithmic-second-chern-defect):
compile the nodal Smith profiles and compare their localized second Chern
length with the global Chern identity.  Degree-specific carrier systems may
serve as regressions for that compiler, but support saturation of
`T_f^log` cannot replace it.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_log_conductor_degree_shift.py
```

The checker verifies (5.1)--(5.4), the finite node normalization quotient,
the smooth-point and node blowup factors, every labelled tree incidence
matrix through six vertices, and (6.2)--(6.3).  The all-curve degree-shift
and determinant-descent statements are the written proofs above, not bounded
computer searches.

## References used by the proof audit

- The Stacks Project,
  [local-cohomology localization](https://stacks.math.columbia.edu/tag/0BK0)
  and [cohomological-dimension bounds](https://stacks.math.columbia.edu/tag/0DXC).
- J. A. Guccione, J. J. Guccione, R. Horruitiner, and C. Valqui,
  [*Some algorithms related to the Jacobian Conjecture*](https://arxiv.org/abs/1708.07936).
