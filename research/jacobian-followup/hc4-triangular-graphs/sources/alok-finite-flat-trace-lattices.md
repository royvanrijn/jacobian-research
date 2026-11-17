# Trace lattices on the smooth degree-six normalization

## Claim boundary

Let

\[
R=\mathbb C[P,Q],\qquad Y=\operatorname{Spec}R\cong\mathbb A^2,
\]

and let

\[
\rho:W\longrightarrow Y
\]

be one of the smooth finite-flat rank-six normalizations left by the
degree-six reduction.  Put

\[
A=\Gamma(W,\mathcal O_W)=\Gamma(Y,\rho_*\mathcal O_W).
\]

Thus `A` is a finite free `R`-algebra of rank six.  In the one-boundary
cases there is a source boundary \(D\cong\mathbb A^1\), with normalization map

\[
\nu:D\longrightarrow B\subset Y,
\]

of degree one.  In the saturated two-boundary case there are disjoint source
curves \(D_2,D_3\cong\mathbb A^1\) mapping with degree one to distinct target
curves `B2,B3`.
Write `b`, or `b2,b3`, for reduced irreducible equations of those target
curves.

This note extracts the exact trace-algebra consequence of the already proved
canonical-divisor formulas.  It does **not** eliminate any of the three
passports.  Its output is a family of free lattices and polynomial matrix
factorizations on which a later compactification or multiplication-algebra
argument can act.

All determinant equations below are understood up to multiplication by an
element of `C*`.  This qualification is necessary: choosing bases of free
`R`-modules changes the determinants by units of `R`, and the only such units
are the nonzero constants.

## 1. Relative duality, the different, and its sign

Because `W` and `Y` are smooth over `C`, the graph of `rho` is a regular
immersion into \(W\times Y\).  Hence `rho` is a local-complete-intersection
morphism.  It is finite flat, so it is syntomic and Gorenstein.  Finite
duality gives

\[
\omega_{W/Y}
\cong
\widetilde{\operatorname{Hom}_R(A,R)}
\cong
\omega_W\otimes\rho^*\omega_Y^{-1}.
\tag{1.1}
\]

Fixing the target form \(dP\wedge dQ\), the previously established
Riemann--Hurwitz calculation identifies this line with

\[
\omega_{W/Y}\cong\mathcal O_W(R_\rho),
\tag{1.2}
\]

where `R_rho` is the ramification divisor.  The three exact possibilities are

\[
R_\rho=D,\qquad R_\rho=2D,
\qquad R_\rho=D_2+2D_3.
\tag{1.3}
\]

There is an important sign convention here.  The relative dualizing line is
the **inverse different**.  Under (1.2), the trace element is the canonical
section

\[
\tau:\mathcal O_W\longrightarrow\mathcal O_W(R_\rho).
\tag{1.4}
\]

Its zero divisor is `R_rho`.  Consequently the different ideal is

\[
\mathfrak D_{W/Y}=\mathcal O_W(-R_\rho),
\tag{1.5}
\]

not `O(R_rho)`.  The different ideal defines the effective divisor
\(R_\rho\), whose norm is the discriminant.  Equivalently, in characteristic
zero,

\[
\operatorname{div}_Y(\det T)=\rho_*R_\rho,
\tag{1.6}
\]

where

\[
T:A\longrightarrow A^*=\operatorname{Hom}_R(A,R),
\qquad
a\longmapsto\bigl(a'\mapsto\operatorname{Tr}_{A/R}(aa')\bigr)
\tag{1.7}
\]

is the trace Gram map.  Thus its determinant is respectively

\[
b,\qquad b^2,\qquad b_2b_3^2.
\tag{1.8}
\]

This recovers the expected discriminant multiplicities, but the lattice
factorization below remembers more than (1.8).

## 2. The half-different lattice

In all three cases write

\[
R_\rho=2E+O,
\tag{2.1}
\]

where `O` is the coefficient-odd part of the ramification divisor.  Concretely,

\[
\begin{array}{c|c|c}
\text{case}&E&O\\ \hline
\text{one-boundary }S_6&0&D\\
\text{one-boundary }A_6&D&0\\
\text{saturated }S_6&D_3&D_2.
\end{array}
\tag{2.2}
\]

Set

\[
L=\mathcal O_W(E),\qquad M=\rho_*L.
\tag{2.3}
\]

The line `L` is projective of rank one over `A`.  Since `A` is finite free over
`R`, `M` is projective of rank six over `R`; Quillen--Suslin makes it free.
The same argument applies to every divisorial lattice

\[
M(Z)=\rho_*\mathcal O_W(Z)
\tag{2.4}
\]

appearing below.

Finite duality gives the exact dual lattice

\[
M^*
\cong
\rho_*\mathcal Hom_W(L,\omega_{W/Y})
\cong
\rho_*\mathcal O_W(E+O).
\tag{2.5}
\]

There are now two natural maps between free rank-six `R`-modules.

First, the canonical divisor section gives

\[
\Phi:A=\rho_*\mathcal O_W\longrightarrow
M=\rho_*\mathcal O_W(E).
\tag{2.6}
\]

Second, multiplication of two sections of `L`, followed by the inclusion
\(\mathcal O_W(2E)\to\mathcal O_W(2E+O)=\omega_{W/Y}\) and the Grothendieck
trace, gives a symmetric pairing on `M`.  Its Gram map is

\[
H:M\longrightarrow M^*.
\tag{2.7}
\]

The ordinary trace map (1.7) factors exactly as

\[
\boxed{T=\Phi^*H\Phi.}
\tag{2.8}
\]

Indeed, the pairing represented by the right-hand side sends \(a,a'\in A\) to
the Grothendieck trace of their product, which is the ordinary algebra trace.
This also proves the determinant identity

\[
\det T=(\det H)(\det\Phi)^2.
\tag{2.9}
\]

Equation (2.9) is the algebraic separation between the coefficient-odd part
of the ramification and its doubled part.  It is compatible with the sign
character of monodromy; it is not by itself a splitting theorem for `A`.

## 3. Cokernels and determinant equations

For an effective Cartier divisor `C` on `W`, the divisor sequence is

\[
0\longrightarrow\mathcal O_W(Z)
\longrightarrow\mathcal O_W(Z+C)
\longrightarrow\mathcal O_C(Z+C)
\longrightarrow0.
\tag{3.1}
\]

The map `rho` is finite, so pushing (3.1) forward remains exact.  If
`C=D_i` maps birationally to `B_i`, the last term is the pushforward by the
normalization map of a line bundle on `D_i`.  Its generic length over `B_i` is
one.  Therefore the determinant of the corresponding inclusion of free
rank-six modules has divisor exactly `B_i`.

Under the inverse-different realization of finite duality, `H` is precisely
the pushforward of

\[
\mathcal O_W(E)\longrightarrow\mathcal O_W(E+O),
\tag{3.2}
\]

and \(\Phi^*\) is the pushforward of the remaining inclusion into
\(\mathcal O_W(2E+O)\).  Thus the divisor sequences apply to the actual maps
in (2.8), not only to abstractly isomorphic source and target modules.

Applying this observation to (2.6) and (2.7) gives the complete table.

| passport | `R_rho=2E+O` | `coker(Phi)` | `det Phi` | `coker(H)` | `det H` | `det T` |
|---|---|---|---|---|---|---|
| one-boundary `S6` | `D` | `0` | `1` | `nu_* O_D(D)` | `b` | `b` |
| one-boundary `A6` | `2D` | `nu_* O_D(D)` | `b` | `0` | `1` | `b^2` |
| saturated `S6` | `D2+2D3` | `(nu3)_* O_D3(D3)` | `b3` | `(nu2)_* O_D2(D2+D3)` | `b2` | `b2*b3^2` |

In the last row `D2` and `D3` are disjoint.  Hence

\[
\mathcal O_{D_2}(D_2+D_3)\cong\mathcal O_{D_2}(D_2),
\tag{3.3}
\]

although retaining `D2+D3` in the table records which global lattice produced
the quotient.

The entries `1` in the table mean a nonzero constant after bases are chosen.
For example, in the `A6` row `H` is a unimodular symmetric form on the free
module `M`; the theorem does not assert that a chosen polynomial Gram matrix
for `H` is literally the identity.

There is also a useful one-boundary filtration.  Put

\[
M_i=\rho_*\mathcal O_W(iD).
\tag{3.4}
\]

If `R_rho=(e-1)D`, finite duality and (3.1) give

\[
M_i^*\cong M_{e-1-i},
\tag{3.5}
\]

and

\[
0\longrightarrow M_{i-1}\longrightarrow M_i
\longrightarrow\nu_*\mathcal O_D(iD)\longrightarrow0.
\tag{3.6}
\]

Each inclusion in (3.6) has determinant `b`, up to `C*`.  For `e=3`, this is
the self-dual chain

\[
A=M_0\subset M_1=M_1^*\subset M_2=A^*.
\tag{3.7}
\]

The first inclusion has matrix `Phi`; the second is its dual, and (2.8) is
their composition through the unimodular middle form.

## 4. Symmetry and the boundary theta line

In the coefficient-odd cases, `H` is not merely a square presentation; it is
symmetric.  Applying `Hom_R(-,R)` to its two-term free resolution identifies
its cokernel with its codimension-one dual.  On the boundary this is the
adjunction identity.

For the one-boundary `S6` case, set

\[
N=\nu_*\mathcal O_D(D).
\]

Then

\[
N\cong\operatorname{Ext}^1_R(N,R),
\tag{4.1}
\]

while finite duality for `nu` identifies the right-hand side with

\[
\nu_*\bigl(\mathcal O_D(-D)\otimes\omega_D\bigr).
\]

Under the geometric identifications that produced `H`, this self-duality is
the one supplied by

\[
\mathcal O_D(2D)\cong\omega_D,
\tag{4.2}
\]

which follows directly from
\(\omega_D=(\omega_W\otimes\mathcal O_W(D))|_D\) and `K_W=D`.  The saturated index-two
component has the same identity because `D3` is disjoint from `D2`.

On the affine line every line bundle is abstractly trivial, so (4.2) is not a
parity contradiction.  Its value is that any projective compactification must
extend a specific symmetric normalization-module presentation.  Degrees and
puncture corrections on that compactification are information that the
affine Picard calculation deliberately forgets.

## 5. Exact matrix coranks

The cokernels in Section 3 determine pointwise coranks without choosing any
matrix entries.  Suppose

\[
0\longrightarrow R^6\mathop{\longrightarrow}^{\Psi}R^6
\longrightarrow\nu_*\mathcal L\longrightarrow0
\tag{5.1}
\]

is one of the divisor-step presentations, where `L` is a line bundle on a
smooth normalization curve.  At a target point `y`, tensoring with `k(y)`
gives

\[
\operatorname{corank}\Psi(y)
=\sum_{t\in\nu^{-1}(y)}
\dim_{\mathbb C}
\frac{\mathcal O_{D,t}}{\mathfrak m_y\mathcal O_{D,t}}.
\tag{5.2}
\]

If a local branch is parametrized by `(p(z),q(z))`, its summand in (5.2) is

\[
\min\{\operatorname{ord}_z p,\operatorname{ord}_z q\},
\tag{5.3}
\]

the multiplicity of that branch.  Consequently:

| matrix | target event | exact corank contribution |
|---|---|---|
| one-boundary `S6` trace matrix `H=T` | generic point of `B` | `1` |
| one-boundary `S6` trace matrix `H=T` | isolated `T(2,3)`, `T(3,4)`, or `T(4,5)` jump branch | `2`, `3`, or `4` |
| one-boundary `S6` trace matrix `H=T` | collision of smooth normalization branches | number of colliding preimages |
| one-boundary `A6` square-root matrix `Phi` | generic point of `B` | `1` |
| one-boundary `A6` square-root matrix `Phi` | the forced `T(2,5)` cusp | `2` |
| one-boundary `A6` square-root matrix `Phi` | the separate `3+3` collision | `2` |
| saturated `S6` matrix `H` | generic point, a two-branch self-collision, or a `2+2+2` three-branch collision of `B2` | `1`, `2`, or `3` |
| saturated `S6` matrix `Phi` | generic point or a two-branch self-collision of `B3` | `1` or `2` |

If a jump and another normalization branch share a target value, their branch
multiplicities add in (5.2).  Likewise, at an intersection of `B2` and `B3`,
the ranks of `H` and `Phi` must be read separately; there is no general rule
that the corank of their product is the sum of their coranks.

The trace matrix itself has a complementary fiber interpretation.  Over a
complex point `y`, the radical of the trace form of the zero-dimensional
fiber algebra is its nilradical.  Write \(D_{\mathrm{bdry}}\) for `D` in a
one-boundary case and for \(D_2\cup D_3\) in the saturated case.  Therefore

\[
\operatorname{corank}T(y)
=6-\#\rho^{-1}(y)
=\sum_{t\in D_{\mathrm{bdry}}\cap\rho^{-1}(y)}(\mu_t-1),
\tag{5.4}
\]

with the evident sum over both boundary curves in the saturated case.  Thus
the `A6` trace matrix has generic corank two, corank four at its degree-five
cusp, and corank four at a `3+3` collision, even though its square-root matrix
`Phi` has coranks one, two, and two.  This distinction is essential.

Since `6` is invertible, the unit and trace split

\[
A=R\cdot1\oplus\ker(\operatorname{Tr}_{A/R}).
\tag{5.5}
\]

The trace-zero summand is free of rank five, and the trace matrix is the
orthogonal block sum of the unit `(6)` and a symmetric rank-five matrix.
Hence the same discriminant and the same trace cokernel admit a rank-five
symmetric presentation.  This sharpening still gives no contradiction: all
the forced trace coranks in (5.4) are at most four.

## 6. Non-Frobenius and non-global-standard-syntomic consequences

The relative dualizing classes are nonzero:

\[
[\omega_{W/Y}]=[D],\quad 2[D],\quad
[D_2]+2[D_3]
\tag{6.1}
\]

in, respectively,

\[
\operatorname{Pic}(W)=\mathbb Z[D]
\quad\text{or}\quad
\operatorname{Pic}(W)=\mathbb Z[D_2]\oplus\mathbb Z[D_3].
\tag{6.2}
\]

This has four immediate algebraic consequences.

1. `A/R` is Gorenstein but not Frobenius: `Hom_R(A,R)` is an invertible
   `A`-module, but it is not a free rank-one `A`-module.

2. `A` is not monogenic as an `R`-algebra.  A monogenic finite flat algebra
   has a presentation `R[z]/(f)` with `f` monic, hence is a global
   hypersurface and has trivial relative dualizing module.

3. More generally, `A` admits no single global square complete-intersection
   presentation

   \[
   A\cong
   R[z_1,\ldots,z_n]/(f_1,\ldots,f_n)
   \tag{6.3}
   \]

   with `(f1,...,fn)` a regular sequence.  In such a presentation both the
   ambient relative differential module and the conormal module have the
   displayed global bases.  The determinant of the relative cotangent
   complex, and hence the relative dualizing line under the Tate map, is
   trivial.  This contradicts (6.1).

4. The natural open overring

   \[
   \Gamma(W\setminus D_{\mathrm{bdry}},\mathcal O)=\mathbb C[x,y]
   \tag{6.4}
   \]

   is not `A[f^-1]` for one element `f`, nor a localization at any finite set
   of elements.  A finite set can be replaced by its product.  If its
   principal open were exactly the affine-plane interior, the zero divisor of
   that product would be `nD` in a one-boundary case, or
   `n2*D2+n3*D3` with both coefficients positive in the saturated case.  This
   would give a nonzero relation in the free Picard group (6.2).  In
   particular, the target branch equation is not a hidden boundary equation
   in `A`: its pullback also has unramified components inside the affine
   plane.

The third statement must not be shortened to “`rho` is not syntomic.”  It is
syntomic and has square complete-intersection presentations locally on `W`.
What fails is one global standard presentation over the whole affine source
and target.  The failure is witnessed by the nontrivial relative-dualizing
class; triviality of that class alone would not conversely construct a global
presentation.

This also explains why a primitive element for the field extension
\(\mathbb C(P,Q)\subset\mathbb C(x,y)\) does not solve the integral problem:
the integral normalization algebra cannot be generated globally by that one
field element.

## 7. Euler characteristic is an identity, not an obstruction

The smooth pair has

\[
\chi_c(W)=
\begin{cases}
\chi_c(\mathbb A^2)+\chi_c(\mathbb A^1)=2,
&\text{one boundary},\\
\chi_c(\mathbb A^2)+2\chi_c(\mathbb A^1)=3,
&\text{two boundaries}.
\end{cases}
\tag{7.1}
\]

On the other hand, finite pushforward computes `chi_c(W)` by integrating the
number of geometric points in the fibers.  If `t` is a boundary point of
local degree `mu_t`, all affine points are reduced because the original
Keller map is etale.  Length six gives

\[
6-\#\rho^{-1}(y)
=\sum_{t\in\rho^{-1}(y)\cap D_{\mathrm{bdry}}}(\mu_t-1).
\tag{7.2}
\]

For one-boundary `S6`, integrating (7.2) gives

\[
(2-1)\chi_c(D)+\sum_t(\mu_t-2)=1+3=4.
\tag{7.3}
\]

For one-boundary `A6` it gives

\[
(3-1)\chi_c(D)+\sum_t(\mu_t-3)=2+2=4.
\tag{7.4}
\]

Thus both yield `chi_c(W)=6-4=2`.  In the saturated case every local degree
is constant along its boundary and the defect is

\[
(2-1)\chi_c(D_2)+(3-1)\chi_c(D_3)=1+2=3,
\tag{7.5}
\]

giving `chi_c(W)=6-3=3`.  Normalization collisions merely group several
summands of (7.2) over one target value; they do not change the integral.

Accordingly, ordinary Euler characteristic cannot remove any of these
passports.  A useful next step must retain the multiplication law, the
filtered degrees supplied by a compactification, or more refined monodromy
than the scalar fiber deficit.

## 8. What the theorem does and does not provide

The strongest concrete new object is the factorization

\[
T=\Phi^*H\Phi
\]

by free polynomial matrices with prescribed determinants, normalization
cokernels, symmetry, and exact rank-drop strata.  It gives a finite target
for symbolic or compactification-based work:

- in the `A6` case, extend the unimodular middle form and the `det Phi=b`
  normalization presentation across target infinity;
- in the one-boundary `S6` case, extend the symmetric theta presentation of
  `nu_*O_D(D)`;
- in the saturated case, couple the symmetric `B2` matrix to the two dual
  `B3` steps at target self-collisions and cross-intersections.

None of the following tempting conclusions is valid.

- `Hom_R(A,R)` is free over `R`, but that does not make it free as an
  `A`-line.
- A square discriminant in the `A6` case records the trivial sign character;
  it does not split the rank-six algebra.
- The norm of every line bundle lands in `Pic(A2)=0`; this does not make the
  source line bundle torsion.
- \(D\cong\mathbb A^1\) makes its line bundles abstractly trivial, but does not erase the
  normalization-module action or its extension across infinity.
- The rank bounds in Section 5 are all compatible with the surviving local
  germs.

Explicit symmetric matrices now show that this boundary is real rather than
provisional.  The [`trace-hostile-matrices.md`](trace-hostile-matrices.md)
note constructs exact `A6` and `S6` normalization-module presentations with
all the determinant, corank, theta, and trace-unit properties above.  They
omit precisely the commutative associative multiplication whose regular
trace would have to equal `T`.

The trace-lattice theorem is therefore a sharpened reduction, not a proof of
the two-dimensional Jacobian conjecture.

## Sources

- The Stacks Project, [Discriminant of a finite locally free
  morphism](https://stacks.math.columbia.edu/tag/0BVH).
- The Stacks Project, [the different, its norm, and the relative dualizing
  line](https://stacks.math.columbia.edu/tag/0BWA).
- The Stacks Project, [quasi-finite syntomic morphisms and invertible relative
  dualizing modules](https://stacks.math.columbia.edu/tag/0DWJ).
- The Stacks Project, [the Tate map from the determinant of the cotangent
  complex to the relative dualizing module](https://stacks.math.columbia.edu/tag/0FKB).
