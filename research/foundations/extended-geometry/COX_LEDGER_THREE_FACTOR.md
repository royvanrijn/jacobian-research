# The first three-factor Cox ledger

This note implements the first tractable case of the Cox-ledger proposal.
The result is a constant-Jacobian morphism between smooth affine boundary
complements, not a polynomial self-map of affine space.

The construction deliberately retains the nonconstant boundary unit.  It
also identifies an integral obstruction which is invisible in the rank
count: for a general degree triple, one free unit direction need not admit a
one-row unimodular completion because the boundary-class map can have finite
cokernel.

Work over a characteristic-zero field `k`.

## 1. What is constructed

Let

\[
 L_i=u_iX+v_iY\qquad (i=1,2,3)
\]

and put

\[
 r_{ij}=u_iv_j-v_iu_j,\qquad
 P=L_1L_2L_3=p_0X^3+p_1X^2Y+p_2XY^2+p_3Y^3.
\]

Choose the irreducible `(1,1,1)` boundary

\[
 m(P)=p_0+p_3=u_1u_2u_3+v_1v_2v_3.
\]

Define the normalized ordered-factor space

\[
 Y=\left\{
 r_{13}=r_{23}=m=1,\quad r_{12}\ne0
 \right\}\subset\mathbb A^6.                         \tag{1}
\]

The target is

\[
 T=\left\{
 p_0+p_3=1,\quad \operatorname{Disc}(P)\ne0
 \right\}\simeq
 \operatorname{Spec}
 k[p_0,p_1,p_2,\operatorname{Disc}(P)^{-1}].          \tag{2}
\]

Multiplication gives `mu:Y -> T`.  Introduce one primitive coordinate `z`
and set

\[
 \boxed{
 \widehat\mu:Y\times\mathbb A^1_z\longrightarrow
 T\times\mathbb A^1_Z,\qquad
 (L_1,L_2,L_3,z)\longmapsto
 \left(P,\frac z{r_{12}}\right).
 }                                                    \tag{3}
\]

The quotient `z/r_(12)` is regular because `r_(12)` is retained as a unit,
not set equal to a constant.

### Theorem 1.1

With the residue volume form described in Section 4, (3) has Jacobian
determinant `1`, up to the harmless simultaneous choice of source and target
orientation.  It is a finite etale morphism of degree six.

Thus (3) is a **torsor-Keller morphism**: a constant-Jacobian morphism
between smooth affine fourfolds with nonconstant units.  It is not called a
Keller map of `A^4`.

The extra suspension dimension is one, exactly the rank of
`O(Y)^*/k^*`.

## 2. The full boundary lattice

On

\[
 X=(\mathbb P^1)^3
\]

remove the four prime divisors

\[
 R_{12},R_{13},R_{23},E=(m=0).
\]

The polynomial `m` is irreducible: viewed as a linear form in `(u_1,v_1)`,
its two coefficients `u_2u_3` and `v_2v_3` are coprime, so Gauss's lemma
applies.  In the standard basis of `Pic(X)=Z^3`, the full boundary-class
matrix is

\[
 B=
 \begin{pmatrix}
 1&1&0&1\\
 1&0&1&1\\
 0&1&1&1
 \end{pmatrix}.                                      \tag{4}
\]

It has

\[
 \ker B=\mathbb Z(1,1,1,-2),\qquad
 \operatorname{coker}B=0.                            \tag{5}
\]

Consequently the projective complement

\[
 U=X\setminus(R_{12}\cup R_{13}\cup R_{23}\cup E)
\]

has

\[
 \mathcal O(U)^*/k^*
 =\langle\varepsilon\rangle\simeq\mathbb Z,\qquad
 \varepsilon=\frac{r_{12}r_{13}r_{23}}{m^2},
 \qquad
 \operatorname{Pic}(U)=0.                            \tag{6}
\]

The three columns selected by `(r_(13),r_(23),m)` have determinant `-1`.
Therefore every projective factor triple in `U` has a unique representative
with

\[
 r_{13}=r_{23}=m=1.
\]

Explicitly, if raw representatives have values `r_(13),r_(23),m`, rescale
them by

\[
 \lambda_1=\frac{r_{23}}m,\qquad
 \lambda_2=\frac{r_{13}}m,\qquad
 \lambda_3=\frac m{r_{13}r_{23}}.                    \tag{7}
\]

This proves

\[
 Y\simeq U,\qquad \varepsilon|_Y=r_{12}.              \tag{8}
\]

In particular, passing to the normalized Cox slice does not erase the
extra unit.

There is a terminology point here.  Since `Pic(U)=0`, the
Neron--Severi universal torsor of `U` is trivial.  The useful object in this
calculation is the restricted Cox characteristic space together with the
unimodular slice (7).  Calling (1) a nontrivial universal torsor would be
incorrect.

The full four-column ledger becomes unimodular after adjoining the primitive
row which reads the `R_(12)` direction:

\[
 \widetilde B=
 \begin{pmatrix}
 1&1&0&1\\
 1&0&1&1\\
 0&1&1&1\\
 1&0&0&0
 \end{pmatrix},
 \qquad \det\widetilde B=1.                           \tag{9}
\]

The last row is realized geometrically by `z -> z/r_(12)` in (3).

## 3. The ordered cubic cover

The multiplication morphism

\[
 \mu:Y\longrightarrow T,\qquad
 (L_1,L_2,L_3)\longmapsto L_1L_2L_3                 \tag{10}
\]

is the ordered-factor cover of the squarefree binary-cubic locus.
Indeed,

\[
 \operatorname{Disc}(P)
 =(r_{12}r_{13}r_{23})^2=r_{12}^2
 \quad\hbox{on }Y.                                   \tag{11}
\]

Every geometric squarefree cubic has six orderings of its three projective
linear factors.  For a fixed ordering, the three scaling conditions

\[
 \lambda_1\lambda_2\lambda_3=1,\qquad
 \lambda_1\lambda_3r_{13}=1,\qquad
 \lambda_2\lambda_3r_{23}=1
\]

have a unique solution because their exponent matrix

\[
 \begin{pmatrix}
 1&1&1\\
 1&0&1\\
 0&1&1
 \end{pmatrix}
\]

has determinant `-1`.  Hence (10) is finite of degree six.

The factor slice before imposing `m=1` has a particularly transparent
description.  The equations `r_(13)=r_(23)=1` make `(L_1,L_3)` an
`SL_2` frame and force

\[
 L_2=L_1+tL_3,\qquad t=r_{12}.
\]

After inverting `r_(12)`, that slice is `SL_2 times G_m`.  Equation `m=1`
then cuts out (1).

## 4. Exact determinant ledger

Consider the ambient square map

\[
 \Theta:
 \mathbb A^6\longrightarrow\mathbb A^6,\qquad
 (L_1,L_2,L_3)\longmapsto
 (p_0,p_1,p_2,m,r_{13},r_{23}).                      \tag{12}
\]

Direct differentiation gives

\[
 \boxed{
 \det D\Theta=-r_{12}r_{13}^2r_{23}^2.
 }                                                    \tag{13}
\]

Take the Poincare residue form `omega_Y` obtained by dividing the standard
six-form by

\[
 d(m-1)\wedge d(r_{13}-1)\wedge d(r_{23}-1)
\]

with the compatible orientation.  Restricting (13) to (1) gives

\[
 \mu^*(dp_0\wedge dp_1\wedge dp_2)
 =-r_{12}\omega_Y.                                   \tag{14}
\]

Now

\[
 d\left(\frac z{r_{12}}\right)
 =\frac{dz}{r_{12}}-\frac{z\,dr_{12}}{r_{12}^2}.
\]

The second summand wedges to zero with the top form in (14), and the first
one cancels the retained unit.  Thus

\[
 \widehat\mu^*
 (dp_0\wedge dp_1\wedge dp_2\wedge dZ)
 =-\omega_Y\wedge dz.                                \tag{15}
\]

Reversing one orientation changes `-1` to `1`.  This proves the
constant-Jacobian assertion.

The new coordinate does not change the fibers: once an ordered
factorization is chosen, `z=r_(12)Z` is forced.  Therefore (3) remains
finite etale of degree six.

Equivalently, the deck group `S_3` acts on `r_(12)=epsilon` by the sign
character after renormalization.  Letting it act on `z` by the same
character makes `Z=z/r_(12)` invariant.  The primitive suspension
trivializes exactly this alternating boundary character.

## 5. The integral correction for arbitrary degrees

For three factors of degrees `(a,b,c)`, with irreducible `E`, the full
boundary matrix is

\[
 B_{a,b,c}=
 \begin{pmatrix}
 b&c&0&1\\
 a&0&c&1\\
 0&a&b&1
 \end{pmatrix}.                                      \tag{16}
\]

The rank of its kernel is always one, but its finite cokernel has order

\[
 \boxed{
 g(a,b,c)=\gcd\left(
 2abc,\,
 c(a+b-c),\,
 b(a+c-b),\,
 a(b+c-a)
 \right).
 }                                                    \tag{17}
\]

The primitive unit relation is

\[
 \frac1g
 \left(
 c(a+b-c),\,
 b(a+c-b),\,
 a(b+c-a),\,
 -2abc
 \right).                                            \tag{18}
\]

Formula (17) is the gcd of the maximal minors of (16).  If one appends any
single integral row to (16), the determinant of the resulting square matrix
is an integral combination of those minors.  Therefore

\[
 \text{a one-row unimodular completion exists}
 \quad\Longleftrightarrow\quad g(a,b,c)=1.            \tag{19}
\]

For `(a,b,c)=(1,1,1)`, `g=1` and (9) is a completion.  In contrast,

\[
 g(1,2,3)=4.
\]

Thus “one free unit direction” is only a rational rank statement.  A
general three-factor construction also has to remove the finite Smith
obstruction, for example by a finite root cover or a stacky Cox chart.

## 6. What this does and does not settle

This construction proves the first sharp rank-one ledger:

1. all four prime boundary columns are retained;
2. the unique boundary unit survives on the source;
3. one primitive coordinate cancels its Jacobian character;
4. the augmented integral ledger is unimodular; and
5. the resulting morphism is finite etale with constant determinant.

It does **not** yet produce a polynomial Keller self-map of affine
four-space.  The source has

\[
 \mathcal O(Y\times\mathbb A^1)^*/k^*
 \simeq\mathbb Z,
\]

so it cannot be affine four-space.  The reciprocal `1/r_(12)` is regular on
the boundary complement but does not extend across `R_(12)=0`.

It also does not yet produce several dicritical divisors for a Keller map of
affine space.  In the natural compactification, the three ordered collision
divisors remain distinct upstairs but all lie over the cubic discriminant.
Extending (3) across them requires a new affine modification; the present
primitive coordinate has a pole there.

Finally, the determinant identity alone does not prove the proposed lower
bound in unit rank for `s>3`.  If a Jacobian is one unit `u`, a single
triangular coordinate `z/u` cancels it even when `u` is a product of many
independent unit generators.  A rank lower bound needs the stronger
reconstruction requirement that the added coordinates separately detect
all independent boundary characters.  For `s=3` the unit rank is one, so
this distinction does not affect the construction above.

The finite-field `C=0` excess of the weighted family is absent because the
chart is `m=1` and has no distinguished `C=0` plane.  Nevertheless (3) is
geometrically a degree-six ordered-factor cover, so it is not itself a
permutation construction.  Arithmetic twisting or a nonsymmetric
factorization core is still needed for that objective.

The construction for every number of ordered linear factors, including the
distinction between separated and compressed unit cancellation, is proved in
[the all-arity Cox-ledger theorem](COX_LEDGER_LINEAR_FACTORS.md).
Allowing the retained collision unit to vanish while adjoining its oriented
square root to the target gives the stronger
[polynomial oriented cubic chart](ORIENTED_CUBIC_COX_CHART.md), with constant
residue Jacobian and two distinct dicritical divisors.

The
[affine-source triple-root chart](AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md)
goes further on the source side: its toric normalization fills to `A^3`,
and orienting the cleared determinant gives four dicritical components with
distinct target images.  Its target remains a five-branch Cox hypersurface
rather than affine three-space.

The optimized endpoint is the
[Danielewski multi-dicritical family](DANIELEWSKI_MULTI_DICRITICAL_FAMILY.md).
It replaces the singular five-branch target by a smooth target `cw=P(x)`,
retains one dicritical affine plane per nonzero root of `P`, and has
bijective reductions when those roots form nonsplit Frobenius orbits.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_three_factor_cox_ledger.py
```

The checker verifies the full and augmented boundary matrices, the general
degree torsion formula on an integer grid, normalization covariance, the
binary-cubic discriminant, both ambient Jacobian identities, and the
constant determinant after primitive suspension.  It also enumerates the
normalized cover over `F_5`, where every split squarefree target has six
preimages.
