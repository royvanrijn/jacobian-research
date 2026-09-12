# Flat three-scalar theories with every finite number of vacua

This note combines the all-degree complete-fiber construction with Bin
Zhu's pullback construction for scalar field theories.  The algebraic input
is the [all-degree rational-fiber theorem](ALL_DEGREE_RATIONAL_FIBERS.md) and
the [exact real-sheet spectrum](REAL_FIBER_SPECTRUM.md).  The physical
framework and the global operator test are due to Zhu:

> Bin Zhu, *Non-injective field redefinitions and quantum inequivalence in
> scalar theories*, [arXiv:2607.18166](https://arxiv.org/abs/2607.18166).

The result stays in three scalar fields for every inverse degree.  The
integer `N` counts vacua and generic inverse sheets; it is not the number of
fields.

## The all-degree scalar theorem

Fix an integer `N>=3`, a spacetime dimension `d`, and positive constants
`f_phi,m`.  There are an explicit polynomial map

\[
 F_N:\mathbb R^3\longrightarrow\mathbb R^3,\qquad
 \det DF_N=1,
\]

and an explicit rational target `Q_N^*` such that the Euclidean three-scalar
action

\[
 S_N[\phi]=\frac{f_\phi^2}{2}\int d^dX\,\sqrt h\,
 \left[
 h^{\mu\nu}g^{(N)}_{ab}(\phi)
       \partial_\mu\phi^a\partial_\nu\phi^b
 +m^2|F_N(\phi)-Q_N^*|^2
 \right],
 \tag{1}
\]

where

\[
 g^{(N)}=(DF_N)^TDF_N,                                  \tag{2}
\]

has the following properties.

1. The field-space metric is positive, flat, and has unit volume density:
   \[
   \det g^{(N)}=(\det DF_N)^2=1.
   \]
2. The potential has exactly `N` isolated, nondegenerate vacua.  Every
   vacuum has rational field coordinates.
3. On a neighborhood of each vacuum, the variables
   \[
   q^I=f_\phi(F_N^I(\phi)-Q_N^{*I})
   \]
   are three free scalar modes of mass `m`.  In fact the exact field
   equations in every local inverse sheet are
   \[
   (-\nabla^2+m^2)(F_N^I(\phi)-Q_N^{*I})=0.
   \]
4. The metric is geodesically incomplete.
5. On the complete regular target locus, the real inverse-sheet counts
   include every value
   \[
   N,N-2,N-4,\ldots,N\bmod 2.
   \tag{3}
   \]
   Each count occurs on a nonempty open chamber, and every chamber in the
   displayed chain has a rational target witness.

Consequently, in Zhu's Schrödinger-space test, the commuting position tuple
`M_{F_N}` has nonconstant joint spectral multiplicity.  It is not globally
unitarily equivalent to the ordinary multiplicity-one coordinate tuple, and
the formal cotangent-lift momenta do not exponentiate with these positions
to a regular representation of the finite-dimensional Weyl relations.

The last sentence is the precise operator-theoretic conclusion proved by
Zhu's argument.  A nonperturbative field-theory completion still requires
boundary conditions at the incomplete ends; the unit-Jacobian identity alone
does not select those boundary conditions.

## Explicit algebraic data

Write `N=2k+1` in odd degree and `N=2k` in even degree, and put

\[
\begin{aligned}
 {\cal S}_{2k+1}
 &=\{2\}\cup\{j,1-j:3\le j\le k+1\},\\
 {\cal S}_{2k}
 &=\{3,4\}\cup\{j,1-j:5\le j\le k+2\},
\end{aligned}
\]

with empty ranges omitted.  Define

\[
 G_N(W)=W(W+1)\prod_{\rho\in{\cal S}_N}(W-\rho),
 \qquad
 \lambda_N=G_N'(0)-G_N'(1),
 \tag{4}
\]

\[
 H_N(W)=\frac{G_N(W)-G_N'(0)W}{\lambda_N}.              \tag{5}
\]

The closed scale formulas are

\[
 \lambda_{2k+1}=(-1)^k k(k!)^2,\qquad
 \lambda_{2k}=\frac{(-1)^k}{144}(7k+2)((k+1)!)^2,
 \tag{6}
\]

so the construction is effective in every degree.  The polynomial `H_N` is
an admissible weighted seed.  Let `F_N=G_{H_N}` be its explicit weighted
suspension from the
[tangent-map core theorem](TANGENT_MAP_CORE.md#weighted-suspension-theorem).
Here the construction constants are `b_0=c=1`, so
`\det DF_N=1`.  Its inverse pencil is

\[
 E_{A,B,C}(W)=H_N(W)-BCW+AC^2.                          \tag{7}
\]

Take

\[
 Q_N^*=(A_N,B_N,C_N)
 =\left(0,-\frac{G_N'(0)}{\lambda_N},1\right).
 \tag{8}
\]

Then

\[
 E_{Q_N^*}(W)=\frac{G_N(W)}{\lambda_N},                 \tag{9}
\]

whose `N` roots

\[
 R_N=\{0,-1\}\cup{\cal S}_N
\]

are distinct integers.  For each `r in R_N`, set

\[
 \gamma_r=-\frac{G_N'(r)}{\lambda_N},\qquad
 a_0=-\frac{1+H_N''(1)}{2+H_N''(1)}
\]

and

\[
 v_r=\left(
 \frac1{\gamma_r},\
 r-\gamma_r,\
 \frac{\gamma_r-1-a_0(r/\gamma_r-1)}
      {(1/\gamma_r)^2}
 \right).
 \tag{10}
\]

All entries are rational, and reconstruction proves the equality of sets

\[
 F_N^{-1}(Q_N^*)=\{v_r:r\in R_N\}.                      \tag{11}
\]

Thus (11) is a complete fiber, not merely a list of `N` colliding points.

## Proof of the scalar statements

Because `DF_N` is invertible everywhere, (2) is the pullback of the
Euclidean metric by a local diffeomorphism.  It is therefore positive and
flat, and its determinant is one.

The potential is

\[
 V_N(\phi)=\frac{f_\phi^2m^2}{2}|F_N(\phi)-Q_N^*|^2.
\]

Its critical-point equation is

\[
 DF_N(\phi)^T(F_N(\phi)-Q_N^*)=0.
\]

Invertibility of `DF_N` makes this equivalent to
`F_N(\phi)=Q_N^*`.  Equation (11) therefore gives exactly `N` critical
points.  At a vacuum `v_r`,

\[
 \operatorname{Hess}_{v_r}V_N
 =f_\phi^2m^2 g^{(N)}(v_r),
\]

which is positive definite.  The vacua are isolated and nondegenerate, and
the kinetic and mass matrices show that every local normal mode has mass
`m`.  Multiplying the Euler--Lagrange equations by the pointwise inverse
Jacobian gives the exact free equations stated above.

If the pullback metric on `R^3` were complete, the local isometry

\[
 F_N:(\mathbb R^3,g^{(N)})\longrightarrow
     (\mathbb R^3,\delta)
\]

would be a covering map.  The Euclidean target is simply connected, so the
cover would be one-to-one.  The `N`-point fiber contradicts this for every
`N>=3`; hence the metric is incomplete.

Finally, the real-root theorem for the pencil (7), together with the
reconstruction isomorphism on `C!=0`, gives the chamber spectrum (3).  In
particular there are positive-measure open regions with `N` and `N-2`
preimages.  The area formula gives the direct-integral decomposition

\[
 L^2(\mathbb R^3_\phi)
 \simeq
 \int_{\mathbb R^3_Q}^{\oplus}
 \mathbb C^{\,N_{F_N}(Q)}\,d^3Q,
 \qquad
 N_{F_N}(Q)=\#F_N^{-1}(Q).
 \tag{12}
\]

The joint spectral multiplicity is therefore nonconstant.  Zhu's
Stone--von Neumann argument now supplies the stated global unitary and Weyl
obstructions.

## Programmable rational fibers

There are two different meanings of “prescribed vacua.”

The theorem above prescribes the number of vacua and gives their rational
coordinates explicitly.  A stronger arithmetic variant prescribes a
split squarefree root configuration.  Given `N` distinct rational numbers,
form their root polynomial and reroot at a root `rho` for which the cubic
coefficient is nonzero.  Such a marking always exists: otherwise the
degree-`N-3` third derivative would vanish at all `N` roots.  The
denominator-free quadratic-gauge construction in
[Root-engineered quadratic gauges](../cancellation/ROOT_ENGINEERED_QUADRATIC_GAUGE.md#4-the-denominator-free-all-degree-map)
then produces a polynomial map whose complete distinguished fiber is
parameterized by precisely those roots.  A rational linear normalization of
one target coordinate changes its constant nonzero Jacobian to one, after
which the same scalar construction applies.

This is control of the fiber polynomial and its reconstructed rational
vacua.  It should not be restated as the currently unproved claim that an
arbitrary preassigned configuration of points in `R^3` can be realized as
the source fiber.

## Attribution boundary

- Zhu supplies the scalar pullback construction, the interpretation of the
  multiple preimages as vacua, the metric-incompleteness argument, and the
  spectral-multiplicity/Weyl obstruction.
- Alexis Gallagher supplies the weighted-lift viewpoint and public
  all-degree seed construction.
- The complete all-degree rational-fiber theorem, rational reconstruction,
  exact real-sheet spectrum, and programmable arithmetic refinement supply
  the `N`-vacuum and arithmetic statements used here.

Because the physics theorem and its global interpretation are Zhu's, any
paper presenting the combined all-degree physics result should either cite
his argument theorem-by-theorem with this division explicit or be prepared
jointly with him.

