# Projective-gradient Segre machinery in all dimensions

## Status and outcome

This note extracts the projective-degree/Segre calculation from the
four-dimensional `HC4` atlas and makes it a repository-wide invariant.  The
central object is always the compactification of the **actual affine map**:

\[
 \Gamma_F=[X_0^m:F_1^h:\cdots:F_n^h]\colon
 \mathbf P^n\dashrightarrow\mathbf P^n,              \tag{0.1}
\]

where \(F=(F_1,\ldots,F_n)\colon\mathbf A^n\to\mathbf A^n\) has polynomial
degree \(m\), and every \(F_i^h\) is homogenized to the common degree \(m\).
For a potential \(\Psi\), take \(F=\nabla\Psi\).  This is not, in general,
the full polar map of the homogenization of \(\Psi\).

The implementation is
[`jcsearch/projective_gradient_segre.py`](jcsearch/projective_gradient_segre.py).
It provides:

- the all-dimensional bijection
  \((g_0,\ldots,g_n)\leftrightarrow(\sigma_1,\ldots,\sigma_n)\);
- constructors for \(\Gamma_F\), \(\Gamma_\Psi\), and the separately named
  full polar map;
- exact homogeneous-gradient integrability checks and Euler
  reconstruction; and
- the all-dimensional smooth-essential normal-slice model, including its
  filtered missing-generator and unit-penultimate laws;
- the singular-stratum join and DVR-torsion profile controlling the first
  componentwise Segre multiplicity; and
- validated records that distinguish a complete multidegree computation
  from a top-degree-only control.

The generated family registry is
[`artifacts/generated-results/projective_gradient_segre_registry.json`](artifacts/generated-results/projective_gradient_segre_registry.json).
It attaches the invariant, at the strongest currently proved level, to the
explicit triangular, cotangent, Meng--Yang/Schur, and restricted-minima HN
families discussed below.  A missing individual Segre degree is recorded as
missing, not inferred from a collision count, Hessian rank, or normalization
packet.

## 1. The all-dimensional transform

Let \(B\) be the base scheme of (0.1), and use the convention

\[
 i_*s(B,\mathbf P^n)=\sum_{k=1}^n\sigma_kH^k,\qquad
 i_*s(B,\mathbf P^n)
 =\pi_*\!\left(\frac E{1+E}\right),                 \tag{1.1}
\]

on the blow-up of the base ideal.  The graph linear system is \(mH-E\).
Consequently its projective degrees are

\[
\boxed{
 g_i=m^i-\sum_{k=1}^i
 \binom{i}{k}m^{i-k}\sigma_k,\qquad 0\leq i\leq n,} \tag{1.2}
\]

with \(g_0=1\).  This triangular relation has the recursive inverse

\[
\boxed{
 \sigma_i=m^i-g_i-
 \sum_{k=1}^{i-1}\binom{i}{k}m^{i-k}\sigma_k.}       \tag{1.3}
\]

Equivalently,

\[
 \sigma_i=(-1)^{i+1}
 \sum_{j=0}^i(-1)^j\binom{i}{j}m^{i-j}g_j.          \tag{1.4}
\]

The checker tests (1.2)--(1.3) as inverse integer transforms for every
\(1\leq n\leq12\) and \(1\leq m\leq7\).  These bounds are regression
bounds only; the functions have no dimension cutoff.

If \(F\) is dominant, then

\[
 g_n=\deg_{\rm aff}(F),\qquad
 \sum_{k=1}^n\binom nk m^{n-k}\sigma_k
 =m^n-\deg_{\rm aff}(F).                            \tag{1.5}
\]

Equation (1.5) is only one weighted aggregate.  It does not recover the
individual \(\sigma_k\).

## 2. Universal restrictions and gradient-specific restrictions

The projective degrees of a dominant rational self-map satisfy

\[
 g_0=1,\qquad 0<g_i\leq m^i,\qquad
 g_i^2\geq g_{i-1}g_{i+1}.                          \tag{2.1}
\]

The last inequalities are the log-concavity of graph multidegrees.  They
are general rational-map restrictions; they do not use integrability.
The intersection-theoretic background is recorded, for example, in
Huh's work on multidegrees and in the mixed-multiplicity formulation of
projective degrees:

- June Huh,
  [*Milnor numbers of projective hypersurfaces and the chromatic polynomial of graphs*](https://arxiv.org/abs/1008.4749);
- Yairon Cid-Ruiz,
  [*Mixed multiplicities and projective degrees of rational maps*](https://arxiv.org/abs/2001.00547).

For an actual gradient \(F=\nabla\Psi\), write

\[
 \Psi=h_{m+1}+h_m+\cdots+h_0
\]

with \(h_j\) homogeneous.  On \(X_0=0\), the base ideal restricts to

\[
 I_{\infty,\mathrm{red}}\ \text{supported on}\
 V(\partial_1h_{m+1},\ldots,\partial_nh_{m+1})
 \subset\mathbf P^{n-1}.                            \tag{2.2}
\]

Thus the reduced infinity support is the singular support of one
homogeneous hypersurface.  The leading generators also obey

\[
 \partial_jF_i^{(m)}=\partial_iF_j^{(m)}             \tag{2.3}
\]

for every \(i,j\), and Euler reconstruction gives

\[
 h_{m+1}=\frac1{m+1}\sum_i x_iF_i^{(m)}.             \tag{2.4}
\]

Equations (2.2)--(2.4) are the first genuine integrability sieve.  A tuple
of same-degree generators with nonzero curl is forbidden as the leading
infinity tuple of a gradient map, even if its abstract Segre numbers obey
(1.2) and (2.1).

This is not yet a forbidden **numerical Segre class**.  The transform
\((g_i)\leftrightarrow(\sigma_i)\) forgets the generator presentation and
its curl syzygies, and many nonisomorphic ideals have the same Segre
degrees.  The current search therefore produces a gradient-realization
condition on base ideals, not an additional inequality in the
\(\sigma_i\) alone.  Any claimed numerical exclusion must eliminate every
integrable ideal realizing that vector.

In characteristic zero, a nonzero \(h_{m+1}\) has at least one nonzero
partial derivative.  Any fixed divisor of (0.1) would have to divide
\(X_0^m\), while that nonzero leading partial makes one homogenized
gradient coordinate indivisible by \(X_0\).  Thus \(\Gamma_\Psi\) has no
fixed divisor and

\[
 g_1=m,\qquad \sigma_1=0.                            \tag{2.5}
\]

The first nonzero Segre component is the effective fundamental cycle of
the top-dimensional normal cone.  Later \(\sigma_k\) are signed and need
not be nonnegative.

For a constant-Hessian potential, the highest determinant layer adds

\[
 \det\operatorname{Hess}(h_{m+1})=0.                 \tag{2.6}
\]

This is a restriction on the leading hypersurface and its Jacobian ideal,
not a formula for the full Segre class.  Low-dimensional
Gordan--Noether/cone theorems and the `HC4` quintic Hessian-discriminant
analysis refine (2.6) in their stated dimensions.  They must not be
relabelled as projective polar-degree theorems.

The full scheme \(B\), not merely its support (2.2), is defined by

\[
 (X_0^m,(\partial_1\Psi)^h,\ldots,
             (\partial_n\Psi)^h).                   \tag{2.7}
\]

Lower homogeneous layers of \(\Psi\) control the \(X_0\)-adic
thickening in (2.7), hence can change every later \(\sigma_k\) without
changing the reduced support.

### 2.1 Smooth-essential normal slices

There is nevertheless an all-dimensional normal-cone calculation when the
top potential is a smooth cone.  Suppose, after a linear change, that

\[
 h_{m+1}=h_{m+1}(u_1,\ldots,u_r),\qquad 1\le r<n,
                                                               \tag{2.8}
\]

and that \((h_{m+1}=0)\subset\mathbf P^{r-1}\) is smooth.  The constant
kernel has dimension \(n-r\), and the reduced infinity base is its vertex

\[
 V=\mathbf P^{\,n-r-1}\subset\mathbf P^{n-1},\qquad
 \operatorname{codim}_{\mathbf P^n}V=r+1.           \tag{2.9}
\]

At the generic point of \(V\), let \(K=k(V)\), put
\(\epsilon=X_0\), and use \(u_1,\ldots,u_r\) as transverse coordinates.
The active homogenized gradient components are

\[
 G_i=\partial_i h_{m+1}
   +\epsilon\partial_i h_m+\cdots
   +\epsilon^{m-1}\partial_i h_2.                  \tag{2.10}
\]

Their special fiber is the Artinian Jacobian complete intersection

\[
 B=K[u_1,\ldots,u_r]/
   (\partial_1h_{m+1},\ldots,\partial_rh_{m+1}).    \tag{2.11}
\]

Consequently

\[
\begin{aligned}
 \operatorname{Hilb}_B(z)&=(1+z+\cdots+z^{m-1})^r,\\
 \dim_KB&=m^r,\\
 \operatorname{socdeg}B&=r(m-1).                  \tag{2.12}
\end{aligned}
\]

In the complete local ring \(K[[\epsilon,u]]\), the \(G_i\) form a regular
sequence.  Their quotient is one-dimensional Cohen--Macaulay, and
\(\epsilon\) is a parameter and hence a nonzerodivisor.  Finite special
fiber and completeness make this quotient finite free of rank \(m^r\) over
\(K[[\epsilon]]\).  After adjoining the compactifying generator
\(\epsilon^m\), the active transverse algebra \(A\) therefore satisfies

\[
 \operatorname{gr}_\epsilon A
 \simeq B\otimes_K K[\epsilon]/(\epsilon^m),\qquad
 \dim_KA=m^{r+1}.                                  \tag{2.13}
\]

The \(n-r\) missing kernel-gradient components are all divisible by
\(\epsilon\).  Suppose one of them has initial form

\[
 \operatorname{in}_\epsilon(G_t)=\epsilon^q s,
 \qquad 1\le q<m,\quad 0\ne s\in B.                \tag{2.14}
\]

The associated graded ideal of \(AG_t\) contains
\((\epsilon^q s)\).  Hence

\[
 \dim_K(AG_t)\ge(m-q)\dim_K(Bs).                   \tag{2.15}
\]

Since \(V\) has degree one, the first Segre component is the generic
transverse length of the full base along \(V\).  Adding the other missing
components can only shorten that quotient, so

\[
 \boxed{\sigma_{r+1}
 \le m^{r+1}-(m-q)\dim_K(Bs).}                     \tag{2.16}
\]

There are two useful exact corollaries.

First, if the penultimate potential layer restricts nontrivially to the
kernel vertex,

\[
 h_m|_V\ne0,                                       \tag{2.17}
\]

then some kernel derivative is nonzero at the generic point.  Its
homogenized component is \(\epsilon\) times a unit.  It kills
\(\epsilon\), after which every missing component vanishes, so

\[
 \boxed{\sigma_{r+1}=m^r.}                         \tag{2.18}
\]

Second, if the kernel vertex is a point, equivalently \(n-r=1\), the base
is zero-dimensional and \(r+1=n\).  Writing
\(\delta=\deg_{\rm aff}\nabla\Psi\), the top projective-degree formula and
the principal-ideal exact sequence give

\[
 \delta=m^n-\sigma_n=\dim_K(AG_t)
 \ge(m-q)\dim_K(Bs).                               \tag{2.19}
\]

If \(s\) is not in the socle of \(B\), then \(\dim(Bs)\ge2\), giving

\[
 \delta\ge2(m-q).                                  \tag{2.20}
\]

> **Theorem `PGS2` -- All-dimensional smooth-essential normal slice.**
> Under (2.8), equations (2.9)--(2.20) hold for the actual affine-gradient
> compactification.  They determine the active transverse Hilbert series
> and length in every dimension, give a filtered missing-generator bound
> for the leading Segre multiplicity, give the exact unit-penultimate law
> (2.18), and give an affine-degree lower bound only when the kernel vertex
> is zero-dimensional.  They do not determine later Segre components or
> apply unchanged to a singular essential top.

The reusable implementation is
`SmoothEssentialGradientNormalSlice` in
[`jcsearch/projective_gradient_segre.py`](jcsearch/projective_gradient_segre.py).
The exact ledger
[`scripts/verify_projective_gradient_normal_slices.py`](scripts/verify_projective_gradient_normal_slices.py)
checks every \(2\le n\le10\), \(2\le m\le7\), and \(1\le r<n\);
these are regression bounds only.  The independent Macaulay2 checker
[`scripts/verify_projective_gradient_normal_slices.m2`](scripts/verify_projective_gradient_normal_slices.m2)
calibrates the complete-intersection and filtered-length formulas.

### 2.2 Singular-essential strata and boundary torsion

For a singular essential top, the support still has an all-dimensional
description, but the smooth length \(m^{r+1}\) must be replaced by a
one-parameter module profile.  Let

\[
 C\subset\operatorname{Sing}(h_{m+1}=0)\subset\mathbf P^{r-1}
                                                               \tag{2.21}
\]

be an irreducible component of dimension \(s\ge0\) and degree \(d\).
If \(V=\mathbf P^{n-r-1}\) is the kernel vertex, the corresponding
infinity component is

\[
 Z_C=\operatorname{Join}(V,C),\qquad
 \dim Z_C=n-r+s,\quad
 \operatorname{codim}_{\mathbf P^n}Z_C=r-s,\quad
 \deg Z_C=d.                                        \tag{2.22}
\]

This follows directly from the top gradient equations: they do not involve
the kernel coordinates, and Euler identifies their nonvertex projective
zeros with the singular locus of \(h_{m+1}=0\).

At the generic point of \(Z_C\), let \(K=k(Z_C)\), take transverse active
coordinates \(z_1,\ldots,z_{r-s-1}\), and let

\[
 B_C=
 K[[z_1,\ldots,z_{r-s-1}]]/
 (\partial_1h_{m+1},\ldots,\partial_rh_{m+1}),
 \qquad \mu_C=\dim_KB_C.                            \tag{2.23}
\]

The \(r\) displayed generators need not be a regular sequence: this is
precisely the excess that distinguishes a singular top.  Keeping the lower
active gradient layers and writing \(R=K[[\epsilon]]\), define

\[
 M_C=K[[\epsilon,z_1,\ldots,z_{r-s-1}]]/
       (G_1,\ldots,G_r).                            \tag{2.24}
\]

Its special fiber is \(B_C\), so completeness and finite special fiber
make \(M_C\) a finite \(R\)-module.  The structure theorem over the DVR
\(R\) gives

\[
 M_C\simeq R^{\rho_C}\oplus
 \bigoplus_{j=1}^{\mu_C-\rho_C}R/(\epsilon^{a_{C,j}}),
 \qquad a_{C,j}\ge1.                               \tag{2.25}
\]

Thus the active algebra after adjoining \(\epsilon^m\) has exact length

\[
 L_C(m)=m\rho_C+
 \sum_{j=1}^{\mu_C-\rho_C}\min(m,a_{C,j}).          \tag{2.26}
\]

After quotienting by the missing kernel-gradient components, let
\(\lambda_C\) be the final generic transverse length.  Every missing
component is divisible by \(\epsilon\), so the special fiber remains
\(B_C\).  Hence

\[
 \mu_C\le\lambda_C\le L_C(m)\le m\mu_C.             \tag{2.27}
\]

The contribution of \(Z_C\) to the codimension-\((r-s)\) Segre component
is

\[
 d\lambda_C.                                       \tag{2.28}
\]

Contributions from components of the same codimension must be summed.
If the active profile has \(\rho_C=0\) and every \(a_{C,j}=1\), then
\(L_C(m)=\mu_C\) and (2.27) is already an equality.  The same equality
holds if some missing kernel-gradient component is \(\epsilon\) times a
unit at the generic point of \(Z_C\), because it kills \(\epsilon\).
Thus either order-one mechanism gives the exact law

\[
 \boxed{d\lambda_C=d\mu_C.}                         \tag{2.29}
\]

If instead the active module is flat, so \(\rho_C=\mu_C\), the filtered
argument from `PGS2` survives.  For
\(\operatorname{in}_\epsilon(G_t)=\epsilon^q s\) it gives

\[
 d\lambda_C\le
 d\bigl(m\mu_C-(m-q)\dim_K(B_Cs)\bigr).             \tag{2.30}
\]

Without flatness, the top class \(s\) alone does not determine how many
\(\epsilon\)-layers survive; the torsion orders in (2.25) are essential.

> **Theorem `PGS3` -- All-dimensional singular-stratum torsion law.**
> Equations (2.21)--(2.30) hold for every nonempty singular component of
> an essential top potential.  The support codimension is \(r-s\), while
> its leading Segre contribution is governed by the transverse Jacobian
> multiplicity, component degree, and the DVR profile
> \((\rho_C;a_{C,j})\).  Singularity dimension or Hessian rank alone does
> not determine that contribution.

The necessity of the torsion profile already appears for the binary
quintic

\[
 h_5=x^3y^2
\]

at the repeated root \([0:1]\), where \(B_C=\mathbf Q[x]/(x^2)\) and
\(\mu_C=2\).  With \(m=4\), the three lower quartics

\[
 h_4=0,\qquad h_4=xy^3,\qquad h_4=y^4
\]

give respectively the profiles

\[
 R^2,\qquad R/(\epsilon^2)\oplus R/(\epsilon),\qquad
 R/(\epsilon)\oplus R/(\epsilon),
\]

and truncated active lengths \(8,3,2\).  The top singularity and
\(\mu_C\) are identical.

The reusable implementation is
`SingularEssentialGradientNormalSlice` in
[`jcsearch/projective_gradient_segre.py`](jcsearch/projective_gradient_segre.py).
The dimension-free regression and binary-quintic certificate are
[`scripts/verify_projective_gradient_singular_slices.py`](scripts/verify_projective_gradient_singular_slices.py);
[`scripts/verify_projective_gradient_singular_slices.m2`](scripts/verify_projective_gradient_singular_slices.m2)
independently checks the three lengths over \(\mathbf Q\).

## 3. What is stable, and what is not

Affine source changes and affine target changes extend to projective
linear transformations preserving the affine chart.  They preserve the
entire projective-degree list of \(\Gamma_F\), and therefore the entire
Segre-degree list.

The top degree is also preserved by several repository constructions:

1. adding identity coordinates or quadratic spectator variables;
2. the cotangent lift
   \[
   \Psi_F(x,y)=y^TF(x)+H(x)
   \]
   when the remaining gradient equations solve uniquely over a generic
   \(F\)-fiber; and
3. the Meng--Yang Schur descent, by the exact Keller-linear solve in
   [`HC4_PROJECTIVE_POLAR_GEOMETRY.md`](HC4_PROJECTIVE_POLAR_GEOMETRY.md).

The interior projective degrees are not stable under these operations.
Two exact cotangent calibrations make this visible.  For

\[
 F_r(x,y)=(x+y^r,y),\qquad
 \Psi_r(x,y,t,u)=t(x+y^r)+uy,
\]

Macaulay2 gives

\[
\begin{array}{c|c|c|c}
r&(g_i)(\Gamma_{F_r})&(g_i)(\Gamma_{\Psi_r})
 &(\sigma_i)(\Gamma_{\Psi_r})\\ \hline
2&(1,2,1)&(1,2,3,2,1)&(0,1,0,-9)\\
3&(1,3,1)&(1,3,5,3,1)&(0,4,-12,8).
\end{array}                                         \tag{3.1}
\]

Thus cotangent lift preserves \(g_{\rm top}=1\) in these examples but
creates new interior degrees.  There is no rule that pads the plane
degree list or recovers the lift from the plane Segre list alone.

Quadratic stabilization of the four-variable triangular potential gives
a second calibration:

\[
\begin{array}{c|c|c}
m&(g_i)(\Gamma_\Psi)&
(g_i)(\Gamma_{\Psi+z^2/2})\\ \hline
2&(1,2,2,2,1)&(1,2,2,2,2,1)\\
3&(1,3,3,3,1)&(1,3,3,3,3,1).
\end{array}                                         \tag{3.2}
\]

At the level of infinity support this is expected: after stabilization,
the new derivative homogenizes to \(X_0^{m-1}z\), so it vanishes on
\(X_0=0\).  The old singular support acquires a projective join/cone and
a new vertex contribution.  Segre classes under joins are structured,
but require the actual extended ideal; see Paolo Aluffi,
[*Tensored Segre classes*](https://arxiv.org/abs/1605.09393) and
[*The Segre zeta function of an ideal*](https://arxiv.org/abs/1606.03098).

Therefore the repository-wide stable invariant currently justified for
these lifts is the top affine degree together with the construction data
\((n,m)\).  The complete vector \((g_i)\), and hence \((\sigma_i)\), must
be recomputed from the lifted base ideal.

## 4. Cotangent boundary incidence

Suppose \(F=(F_1,\ldots,F_q)\) has degree \(m\), and the leading part of
the cotangent potential is

\[
 p_{m+1}(x,y)=y^TF^{(m)}(x).
\]

Then its leading gradient tuple is

\[
 \left((DF^{(m)}(x))^Ty,\ F^{(m)}(x)\right).         \tag{4.1}
\]

Consequently the reduced infinity support lies on the incidence

\[
 F^{(m)}(x)=0,\qquad (DF^{(m)}(x))^Ty=0.             \tag{4.2}
\]

This is the projective conormal/kernel incidence naturally missing from a
plane finite-normalization packet.  Formula (4.2) explains both sides of
the plane-to-`HC4` gap:

- a plane packet can determine the generic fiber count transferred to
  the cotangent gradient;
- it does not determine the rank strata of \(DF^{(m)}\), the projective
  \(y\)-fibers in (4.2), or the \(X_0\)-adic thickening from lower
  layers.

Therefore boundary normalization data does not determine an infinity
Segre class.  To compute one needs the actual representative \(F\), its
common coordinate degree, the full ideal (2.7), and a Rees/Segre
calculation.

## 5. Schur descent and Hessian discriminants

For the Meng--Yang controls, the top-degree transport is exact:

\[
\begin{array}{c|c|c|c}
\text{family}&n&m&g_n\\ \hline
\text{doubled potential}&6&7&3\\
\text{Schur-descended potential}&5&13&3.
\end{array}
\]

Hence their weighted corrections are

\[
 7^6-3=117646,\qquad 13^5-3=371290.                 \tag{5.1}
\]

Schur descent changes both \(n\) and \(m\), and its nonlinear graph
substitution is not a projective linear equivalence.  Equation (5.1)
does not give a transform between the two interior projective-degree
lists.

Likewise, a Hessian discriminant belongs first to the leading potential
or Schur denominator.  It affects the possible rank strata in (2.6) and
the geometry supporting (2.2), but it is not itself:

- the discriminant of the full polar map;
- the base ideal (2.7);
- a projective degree; or
- a Segre class.

The `HC4` squarefree rank-three obstruction is therefore recorded as a
restriction on possible leading base supports.  Translating it into
individual forbidden Segre vectors requires a local normal-cone
calculation on the surviving nonsquarefree strata.

## 6. Family attachment ledger

The machine-readable registry uses three evidence levels.

| family | \(n\) | \(m\) | attached invariant | status |
|---|---:|---:|---|---|
| triangular constant-Hessian, \(r=2\) | 4 | 2 | full \(g,\sigma\); separate full-polar comparison | exact Macaulay2 |
| triangular constant-Hessian, \(r=3\) | 4 | 3 | full \(g,\sigma\); separate full-polar comparison | exact Macaulay2 |
| plane triangular cotangent, \(r=2,3\) | 4 | 2,3 | full \(g,\sigma\) and source-plane \(g,\sigma\) | exact Macaulay2 |
| one-variable quadratic stabilization, \(r=2,3\) | 5 | 2,3 | full \(g,\sigma\) before and after | exact Macaulay2 |
| plane quartic packet target | 4 | 2 or 3 | \(g_4=4\), aggregate only | conditional; no representative |
| Meng--Yang doubled control | 6 | 7 | \(g_6=3\), aggregate only | proved top degree |
| Meng--Yang Schur control | 5 | 13 | \(g_5=3\), aggregate only | proved top degree |
| homogeneous cotangent HN witness | 38 | 3 | construction metadata | explicit; \(g\) uncomputed |
| rank-reduced cotangent HN witness | 44 | 3 | construction metadata | explicit; \(g\) uncomputed |
| nonhomogeneous cotangent HN witness | 40 | 3 | construction metadata | explicit; \(g\) uncomputed |

For the HN entries, noninjectivity and nonzero Jacobian imply a dominant
non-birational gradient Keller map, but the exact generic degree has not
been computed in their canonical certificates.  The registry therefore
does not guess it.

Universal gradient-map coefficient schemes acquire an additional
projective layer as follows:

1. impose the usual coefficient equations for symmetric Jacobian and
   constant determinant;
2. extract the common top tuple \(F^{(m)}\);
3. impose or verify the curl equations (2.3);
4. form the universal homogeneous ideal
   \((X_0^m,F_1^h,\ldots,F_n^h)\); and
5. stratify by projective degrees, equivalently by the Segre transform
   (1.2).

This is a constructible stratification problem for the universal graph,
not a claim that one Segre vector is constant on the whole coefficient
scheme.

## 7. The quintic `HC4` top-gradient consumer

The first coefficient-level consumer is now implemented in
[`scripts/analyze_hc4_quintic_infinity_rees.py`](scripts/analyze_hc4_quintic_infinity_rees.py).
It constructs the universal \(56\)-coefficient four-variable \(h_5\),
checks its curl, Euler, Hessian--Euler, Koszul, and midpoint-parity
relations, removes generic Hessian rank four using the constant-determinant
top face, and builds the essential rank-one/two/three strata.

On their smooth essential loci the top gradient ideals are complete
intersections of linear type.  Their Rees ideals contain only Koszul
relations and inactive-target linear equations; this is independently
certified by
[`scripts/verify_hc4_quintic_infinity_rees_strata.m2`](scripts/verify_hc4_quintic_infinity_rees_strata.m2).
The reduced base support nevertheless gives an exact atlas filter:

\[
\begin{array}{c|c|c}
\text{top type}&\operatorname{codim}B_\infty&
\text{affine-degree }2,3\text{ rows}\\ \hline
\text{rank 1}&2&260,249\\
\text{smooth essential rank 2}&3&58,57\\
\text{smooth essential rank 3}&4&1,1.
\end{array}
\]

Singular binary or ternary tops feed the lower-codimension columns
according to the dimension of their singular locus.  The generated
artifact
[`hc4_quintic_infinity_rees_strata.json`](artifacts/generated-results/hc4_quintic_infinity_rees_strata.json)
keeps the pure-top complete-intersection Segre vectors separate from the
unknown lower-layer normal-cone multiplicities.

The first lower-layer consequence is theorem `HC4PPG7`.  On the unique
vertex of a smooth essential rank-three top quintic, the three active
gradient components form a flat rank-\(64\) family over the compactifying
parameter; truncation at order four has length \(256\).  The missing
component begins with \(\epsilon s_3\), and the degree-nine socle of the
ternary \((4,4,4)\) Jacobian complete intersection forces its ideal to have
length at least six when \(s_3\ne0\).  Thus affine degrees two and three
cannot realize the codimension-four Segre signatures.  This is a lower
normal-cone restriction that the numerical transform alone cannot see.
It is exactly the \((n,m,r)=(4,4,3)\) specialization of `PGS2`.

The next consumer, theorem `HC4PPG8`, treats the codimension-three
packets.  A nonzero restriction of \(h_4\) to the rank-two kernel line
synchronizes a constant Hessian-kernel direction.  The squarefree binary
Hessian-discriminant branch then reaches the common-direction obstruction;
on the nonsquarefree remainder the generic transverse Segre multiplicity is
forced to \(\sigma_3=16\).  For rank-three tops with isolated ordinary
singularities, the Schur cubic must vanish at every singular point.  The
shared invariant therefore carries both numerical and coefficient-incidence
data without identifying the two packet types.  The exact
\(\sigma_3=16\) normal-slice law is the \((n,m,r)=(4,4,2)\)
specialization of `PGS2`.

The first direct singular-top consumer is theorem `HC4PPG9`.  For an
essential binary quintic root of multiplicity \(e\), the transverse
Jacobian length is \(e-1\).  On the open stratum where a redundant active
gradient has order-one unit coefficient, `PGS3` makes that root contribute
exactly \(e-1\) to \(\sigma_2\).  Thus a root partition with \(q\)
distinct roots forces \(\sigma_2=5-q\); the generic double-root packet
shrinks from \(260,249\) numerical rows to \(51,50\).  The complementary
higher-torsion locus remains open.

## 8. Reproduction and next searches

For cleanup and provenance review only, without importing SymPy, invoking
Macaulay2, or rebuilding a ledger, run:

```bash
python3 scripts/audit_projective_gradient_segre_artifacts.py
```

The audit pins the registry, smooth-slice, and singular-slice artifacts at
whole-file SHA-256
`1678eac19cc8e59a123ec84836f8f2a89f3b697a29c241e3b26e6987180fd00f`,
`5853c8fa609879663b31f680591a5e612ab944b1637902de5dcd115c9400837b`,
and
`c6971874b5359e4aed11a8918328804f9ffdd6e67811f49c9ff79b2a8c5d7b72`.
It also checks the imported helper and the four Macaulay2 calibration files.
The three Python verifiers now compare their generated serialization with the
committed bytes by default; rewriting requires an explicit `--write`.

Run the dimension-free algebra, constructors, integrability test, family
registry, and independent Macaulay2 calibrations with:

```bash
.venv/bin/python scripts/verify_projective_gradient_segre_machinery.py
M2 --script scripts/verify_projective_gradient_segre_families.m2
.venv/bin/python scripts/verify_projective_gradient_normal_slices.py
M2 --script scripts/verify_projective_gradient_normal_slices.m2
.venv/bin/python scripts/verify_projective_gradient_singular_slices.py
M2 --script scripts/verify_projective_gradient_singular_slices.m2
.venv/bin/python scripts/analyze_hc4_quintic_infinity_rees.py
M2 --script scripts/verify_hc4_quintic_infinity_rees_strata.m2
.venv/bin/python scripts/verify_hc4_rank3_vertex_colength.py
M2 --script scripts/verify_hc4_rank3_vertex_colength.m2
.venv/bin/python scripts/verify_hc4_codim3_gradient_strata.py
M2 --script scripts/verify_hc4_codim3_gradient_strata.m2
.venv/bin/python scripts/verify_hc4_binary_root_partition_segre.py
M2 --script scripts/verify_hc4_binary_root_partition_segre.m2
```

Use `--write` on one of the first three Python commands only after an
intentional theorem, regression-range, or schema change.

The existing `HC4` atlas now imports the same transform:

```bash
.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py
```

The next searches have sharply separated inputs and outputs.

1. **Exceptional codimension-three packets.**  Treat \(h_4|_K=0\) in
   essential rank two, the nonsquarefree binary-Hessian
   \(\sigma_3=16\) row, and the lower normal cones at isolated ternary
   singularities subject to the Schur-cubic incidence.
2. **Exceptional codimension-two torsion.**  On the repeated-binary-root
   packet, impose the constant-Hessian determinant faces on the failure of
   the `HC4PPG9` active-unit condition and compute the remaining DVR
   torsion orders.
3. **Cotangent incidence strata.**  Compute Segre classes of (4.2) for
   explicit plane representatives, stratified by the leading Jacobian
   rank, before adding lower-layer thickenings.
4. **Stable join law.**  Use the actual stabilized ideal and Segre-zeta
   machinery to determine when identity/quadratic padding has a closed
   transform on \(g\)-vectors.
5. **Schur pairs.**  Compute both graph multidegree lists for an explicit
   source/descended pair; compare more than the already proved common top
   degree.
6. **HN minima.**  Start with the 38-variable homogeneous witness, whose
   construction is smaller than the 40/44-variable internal records, and
   compute \(g_n\) before attempting all intermediate degrees.
7. **Boundary versus normal cone.**  For a fixed finite-normalization
   packet, vary polynomial representatives and measure which
   \(X_0\)-adic base thickenings, hence which Segre vectors, actually
   occur.

These are computations and realization problems.  No new exclusion of a
Keller map follows merely from the existence of the registry.
