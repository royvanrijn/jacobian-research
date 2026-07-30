# Rank-four collision frames and the cross-ratio obstruction

The rank-four collision tower separates two issues that coincide in rank
three.

1. Ordered pairs do not label a quartic fiber completely.  The ordered
   three-root cover is the full `S_4` frame torsor.
2. Even after all four roots are labeled, two primitive presentations need
   not be related by a projective root-coordinate change.  Their exact defect
   is one cross-ratio equation.

For a primitive change

\[
 u=q_0+q_1r+q_2r^2+q_3r^3,                            \tag{0.1}
\]

write `e_1,e_2` for the first two elementary symmetric functions of the four
old roots.  The projective defect is

\[
\boxed{
 \Psi(q;r)
 =q_2^2-q_1q_3+q_2q_3e_1+q_3^2e_2.}                 \tag{0.2}
\]

On the squarefree framed overlap, the unique projective transformation
matching the first three roots matches the fourth if and only if `Psi=0`.
Thus the collision frame closes the finite `S_4` labeling problem, but full
projective presentation descent holds only on the hypersurface (0.2).

This is an exact obstruction to the canonical `PGL_2` root transport.  It is
not a proof that the rank-four Keller incidence itself cannot descend:
a genuinely nonlinear ambient Keller equivalence need not arise from a
projective transformation of the root line.

Work over a characteristic-zero field, or over a reduced base on which the
displayed discriminants and normalization coefficients are units.

## 1. From the collision sheet to the full quartic frame

Let `E -> S` be finite etale of rank four.  The ordered off-diagonal
collision cover

\[
 \operatorname{Off}_2(E/S)
 =(E\times_SE)\setminus\Delta_E                      \tag{1.1}
\]

has rank twelve.  A geometric ordered pair leaves two roots unlabeled, and
its pointwise stabilizer in `S_4` is `S_2`.  Therefore `Off_2` is not the
normal closure or frame torsor.

The tensor collision algebra and diagonal/off-diagonal interface in (1.1)
are credited to Chloe van der Vlugt's
*Collision Ideals and Off-Diagonal Sheets* as recorded in the
[external audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md).  The higher
configuration and cross-ratio deductions in this note are results of this
repository.

Pass instead to

\[
 \operatorname{Conf}_3(E/S)
 =
 E^3\setminus\bigcup_{i<j}\Delta_{ij}.                \tag{1.2}
\]

Three disjoint sections of a rank-four finite-etale cover have a unique
rank-one open-and-closed complement.  Appending that fourth section gives

\[
\boxed{
 \operatorname{Conf}_3(E/S)
 \simeq
 \operatorname{Isom}_S(\{1,2,3,4\}_S,E).}             \tag{1.3}
\]

Thus `Conf_3` has rank `4!=24` and is canonically the `S_4` frame torsor.
This is the precise quartic replacement for the cubic identity
`Off_2=Conf_3`.

## 2. Three-point interpolation and the fourth-root residual

On the frame torsor, let `(r_1,r_2,r_3,r_4)` and
`(u_1,u_2,u_3,u_4)` be the roots of two primitive presentations, paired by
the underlying algebra isomorphism.  Interpolate the first three pairs by
the signed-minor construction from the
[rank-three audit](RANK_THREE_COLLISION_DESCENT.md):

\[
 M=
 \begin{pmatrix}
 r_1&1&-u_1r_1&-u_1\\
 r_2&1&-u_2r_2&-u_2\\
 r_3&1&-u_3r_3&-u_3
 \end{pmatrix},
\qquad
 g=\begin{pmatrix}m_0&m_1\\m_2&m_3\end{pmatrix},     \tag{2.1}
\]

where `m_j=(-1)^j det(M_hat_j)`.  Then

\[
 g(r_i)=u_i\quad(i=1,2,3),\qquad
 \det g=V(r_1,r_2,r_3)V(u_1,u_2,u_3).                 \tag{2.2}
\]

Now impose (0.1).  Pairwise differences factor as

\[
 u_i-u_j=(r_i-r_j)D_{ij},                             \tag{2.3}
\]

where

\[
 D_{ij}
 =q_1+q_2(r_i+r_j)
   +q_3(r_i^2+r_ir_j+r_j^2).                         \tag{2.4}
\]

Let

\[
 V_4(r)=\prod_{1\le i<j\le4}(r_i-r_j).
\]

Exact maximal-minor expansion gives the fourth-point residual

\[
\boxed{
 m_0r_4+m_1-u_4(m_2r_4+m_3)
 =-V_4(r)\Psi(q;r).}                                  \tag{2.5}
\]

Since `V_4(r)` is a unit on the frame torsor, `g(r_4)=u_4` exactly when
`Psi=0`.

## 3. Cross-ratio interpretation

Use the labeled cross-ratio

\[
 \chi(r)
 =\frac{(r_1-r_3)(r_2-r_4)}
        {(r_1-r_4)(r_2-r_3)}.                         \tag{3.1}
\]

Equation (2.3) gives

\[
 \frac{\chi(u)}{\chi(r)}
 =\frac{D_{13}D_{24}}{D_{14}D_{23}}.                 \tag{3.2}
\]

The numerator of the difference factors as

\[
\boxed{
 D_{13}D_{24}-D_{14}D_{23}
 =-(r_1-r_2)(r_3-r_4)\Psi(q;r).}                     \tag{3.3}
\]

All root differences and all `D_ij` are units on the framed primitive
overlap.  Consequently

\[
\boxed{
 \chi(u)=\chi(r)
 \Longleftrightarrow
 \Psi(q;r)=0
 \Longleftrightarrow
 \text{one }g\in PGL_2\text{ sends every }r_i\text{ to }u_i.}  \tag{3.4}
\]

When (3.4) holds, interpolation using any three labels produces the same
projective transformation.  Simultaneous `S_4` relabeling therefore leaves
it unchanged, and it descends from the frame torsor with the usual
composition cocycle.  When `Psi!=0`, different choices of the omitted fourth
label give incompatible projective interpolants.  Framing alone cannot
remove that defect.

## 4. Primitive and projective boundaries are different

The determinant of the basis change from
`(1,r,r^2,r^3)` to `(1,u,u^2,u^3)` is the Vandermonde ratio

\[
\boxed{
 \Theta_4(q;r)
 =\prod_{i<j}D_{ij},\qquad
 V_4(u)=V_4(r)\Theta_4(q;r).}                          \tag{4.1}
\]

Thus `u` remains primitive exactly when `Theta_4` is a unit.  If
`f_j=e_j(u_1,u_2,u_3,u_4)`, the normalized quartic presentation open also
requires `f_1f_3` to be a unit: `f_3` normalizes the linear coefficient and
`f_1` keeps the normalized cubic coefficient nonzero.

The projective condition `Psi=0` is independent of these presentation-open
conditions.  In particular, for a genuinely quadratic change `q_3=0`,

\[
\boxed{\Psi=q_2^2.}                                   \tag{4.2}
\]

Over a field, no quadratic change with `q_2!=0` is projective on four
distinct points, even when it remains a perfectly good primitive generator.
The scheme equation retains the double structure `q_2^2=0`; only its reduced
support is the affine locus `q_2=0`.

When `q_3` is nonzero, non-affine projective changes do exist.  Their
coefficient locus is

\[
 q_1
 =\frac{q_2^2+q_2q_3e_1+q_3^2e_2}{q_3}.             \tag{4.3}
\]

These are cubic polynomial representatives, modulo the quartic relation, of
Möbius transformations on the four selected roots.

## 5. Equation on the universal quartic Keller chart

The rank-four inverse polynomial in the
[universal relative Keller map](UNIVERSAL_RELATIVE_KELLER_MAP.md) is

\[
 E_{u,\pi,b,c}(R)
 =u\pi^4R^4+\pi R^3+bR^2+R-\frac c2.                 \tag{5.1}
\]

Writing it as `a_4 prod(R-r_i)`, with `a_4=u pi^4`, gives

\[
 e_1=-\frac1{u\pi^3},\qquad
 e_2=\frac b{u\pi^4},\qquad
 e_3=-\frac1{u\pi^4},\qquad
 e_4=-\frac c{2u\pi^4}.                              \tag{5.2}
\]

Multiplying (0.2) by the invertible leading coefficient produces the
polynomial equation

\[
\boxed{
 \widehat\Psi
 =
 u\pi^4(q_2^2-q_1q_3)
 -\pi q_2q_3
 +bq_3^2.}                                           \tag{5.3}
\]

On the universal open `u*pi!=0`,

\[
 \widehat\Psi=0\Longleftrightarrow\Psi=0.             \tag{5.4}
\]

Thus the projectively transportable part of the rank-four presentation
overlap is one explicit hypersurface in the actual Keller parameter--target
coordinates.  The constant coordinate `q_0` and target constant `c` do not
enter: cross-ratio is translation-invariant and depends only on the first two
symmetric functions.

## 6. Exact witness cards

Take the ordered roots `(1,2,3,4)`.  Their normalized universal target is

\[
 (u,\pi,b,c)
 =\left(-\frac{25}{2},\frac15,-\frac7{10},\frac{24}{25}\right).
                                                               \tag{6.1}
\]

Two primitive changes separate the loci:

\[
\begin{array}{c|c|c|c}
q(r)&(u_1,u_2,u_3,u_4)&\Theta_4&\Psi\\ \hline
r+r^2&(2,6,12,20)&\ne0&1\\
35r+r^3&(36,78,132,204)&\ne0&0.
\end{array}                                             \tag{6.2}
\]

The first is a literal primitive-presentation overlap on which the
projective transport fails.  The second is a non-affine point of the
projective hypersurface and the interpolating `PGL_2` map matches all four
roots.

## 7. Consequence for Keller descent

The collision input is now complete in rank four:

- `Off_2` gives ordered pairs but retains an `S_2` stabilizer;
- `Conf_3` is the full `S_4` frame torsor;
- on that torsor, `Psi` is the exact obstruction to descending the canonical
  projective root transport.

The
[all-rank projective-descent theorem](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md)
identifies this as the `N=4` determinant of the universal matrix with columns
`1,r,u,r*u`.  In general its minors cut out a smooth codimension-`N-3`
projective locus, so `Psi` is the first member of a uniform family of
projective-moduli obstructions.

This changes the rank-four research question.  It is not enough to imitate
the rank-three target-localized `PGL_2` lift on the whole presentation
groupoid.  Such a lift can only cover \(\widehat\Psi=0\).  Away from that
hypersurface, any Keller-incidence descent must realize a genuinely
nonprojective Tschirnhaus change or use an enhancement that forgets the root
embedding while retaining enough ambient data.

No global nonexistence conclusion follows.  A polynomial source--target
equivalence of the rank-four Keller family could act on the finite fiber by a
nonprojective cubic interpolation.  The theorem above supplies a minimal
exact falsification card and a natural two-stratum program:

1. on \(\widehat\Psi=0\), attempt the target-localized projective lift and audit
   its boundary denominators;
2. on \(\widehat\Psi\ne0\), test the witness `q(r)=r+r^2` against stable boundary,
   conormal, and marked-normalization invariants of the rank-four Keller
   family;
3. only promote a global descent claim if those nonprojective transitions
   acquire an exact polynomial cocycle.

## 8. Exact regression

Run

```bash
.venv/bin/python scripts/verify_rank_four_collision_cross_ratio.py
```

The checker verifies:

- the ranks and stabilizers of `Conf_2` and `Conf_3`;
- the signed-minor three-point interpolation determinant;
- the exact fourth-root residual `-V_4*Psi`;
- the equivalent cross-ratio numerator factorization;
- the primitive determinant `Theta_4`;
- the normalized universal-quartic coefficient formulas and cleared equation
  (5.3); and
- both witness cards in (6.2).

Every calculation is an exact symbolic identity.  The checker does not test
or assume that all Keller-incidence equivalences arise from `PGL_2`.
