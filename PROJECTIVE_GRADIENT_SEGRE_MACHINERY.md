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

## 7. Reproduction and next searches

Run the dimension-free algebra, constructors, integrability test, family
registry, and independent Macaulay2 calibrations with:

```bash
.venv/bin/python scripts/verify_projective_gradient_segre_machinery.py
M2 --script scripts/verify_projective_gradient_segre_families.m2
```

The existing `HC4` atlas now imports the same transform:

```bash
.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py
```

The next searches have sharply separated inputs and outputs.

1. **Integrable Segre sieve.**  Enumerate log-concave degree lists, but
   realize only base ideals whose leading generators pass (2.3)--(2.4).
2. **Cotangent incidence strata.**  Compute Segre classes of (4.2) for
   explicit plane representatives, stratified by the leading Jacobian
   rank, before adding lower-layer thickenings.
3. **Stable join law.**  Use the actual stabilized ideal and Segre-zeta
   machinery to determine when identity/quadratic padding has a closed
   transform on \(g\)-vectors.
4. **Schur pairs.**  Compute both graph multidegree lists for an explicit
   source/descended pair; compare more than the already proved common top
   degree.
5. **HN minima.**  Start with the 38-variable homogeneous witness, whose
   construction is smaller than the 40/44-variable internal records, and
   compute \(g_n\) before attempting all intermediate degrees.
6. **Boundary versus normal cone.**  For a fixed finite-normalization
   packet, vary polynomial representatives and measure which
   \(X_0\)-adic base thickenings, hence which Segre vectors, actually
   occur.

These are computations and realization problems.  No new exclusion of a
Keller map follows merely from the existence of the registry.
