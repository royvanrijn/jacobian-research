# Fitting-denominator extraction for `HC4`

## Status

Registry entry: `HC4QSE5` (**partial**).

This note records an **exact specialized and finite-field module
experiment**, not a completed characteristic-zero Fitting computation.
It corrects the proposed module presentation, identifies a new rational
nilpotence-jump point, and explains why cube-certificate torsion cannot by
itself equal the reduced exceptional Schur locus.

The even-block integral annihilator and function-field lift computations
described below each reached their declared 900-second Singular timeout.
The relation extraction needed before the full zeroth Fitting ideal was
therefore not reached.  No partial output from those runs is used as a
certificate.  The reproducible finite-field scans and exact rational
specializations are recorded in
[`hc4_fitting_denominator_extraction.json`](artifacts/generated-results/hc4_fitting_denominator_extraction.json).
The companion bounded fourth-power ledger is
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json).
For maintenance without any algebra or scan, run
`python3 scripts/audit_hc4_fitting_denominator_artifacts.py`; it only verifies
the committed bytes, generating-source hash, and the fail-closed scope below.

## 1. The canonical module

Put

\[
A=\mathbb Q[\mu,\nu],\qquad
S=A[s_0,\ldots,s_{14}],
\]

and let \(f_1,\ldots,f_{114}\in S_2\) be the primitive quadratic Schur
equations obtained from the two-parameter calculation.  The 114 quadrics
do not canonically define a map \(A^{114}\to A^{15}\): their coefficient
vectors lie in

\[
S_2\simeq A^{120}.
\]

The coefficient cubes occur one degree later.  The canonical
multiplication presentation is

\[
\Phi:A^{114}\otimes_A S_1\simeq A^{1710}
   \longrightarrow S_3\simeq A^{680},
\qquad e_j\otimes s_i\longmapsto s_i f_j.       \tag{1.1}
\]

Let

\[
D:A^{15}\longrightarrow A^{680},
\qquad e_i\longmapsto s_i^3.
\]

The intrinsic cube-certificate torsion is

\[
T=\frac{\operatorname{im}\Phi+\operatorname{im}D}
        {\operatorname{im}\Phi}
  =\left\langle[s_0^3],\ldots,[s_{14}^3]\right\rangle_A
     \subset (S/(f_1,\ldots,f_{114}))_3.         \tag{1.2}
\]

If

\[
K=\ker\left(A^{15}\xrightarrow{D}\operatorname{coker}\Phi\right),
\]

then \(K\to A^{15}\to T\to0\) is the desired 15-generator
presentation.  Its number of relation columns is not a priori 114.  In
particular,

\[
\operatorname{Fitt}_0(T)=I_{15}(K),\qquad
\operatorname{Ann}_A(T)=(\operatorname{im}\Phi:\operatorname{im}D).
                                                        \tag{1.3}
\]

Any direct \(A^{114}\to A^{15}\) matrix requires noncanonical choices of
linear multipliers and is already a choice of cube certificates.

## 2. Primitive rows and chart boundary

The 114 cleared equations have parameter contents

| content | multiplicity |
|---|---:|
| \(1\) | 2 |
| \(2\) | 10 |
| \(4\) | 4 |
| \(8\) | 4 |
| \(16\) | 1 |
| \(\nu\) | 1 |
| \(2\nu\) | 13 |
| \(4\nu\) | 14 |
| \(8\nu\) | 19 |
| \(16\nu\) | 2 |
| \(4\nu^2\) | 16 |
| \(8\nu^2\) | 26 |
| \(8\nu^3\) | 2 |

The calculation divides out these contents before constructing (1.1).
This is harmless on the original pivot chart \(D(\nu)\), but the resulting
extension across \(\nu=0\) is chart-dependent.  Indeed its character-block
ranks cease to be permutation-symmetric on that boundary.  Consequently
the module constructed from these 114 equations must be interpreted over
\(A[\nu^{-1}]\).  Fermat lies outside this chart and needs a separate
boundary presentation.

If one nevertheless uses the primitive extension over all of
\(\operatorname{Spec}A\), several cube classes are supported on the whole
line \(\nu=0\).  This is consistent with the separate one-parameter
calculation, where fourth powers, rather than cubes, give the uniform
nilpotence certificate.  It is not evidence for a reduced exceptional
curve.

## 3. Character blocks and determinantal ranks

The sign-character splitting reduces (1.1) to one even block and three
nontrivial blocks.  A representative even block has size

\[
A^{441}\longrightarrow A^{191}
\]

and generic rank 190.  A representative nontrivial block has size

\[
A^{423}\longrightarrow A^{163}
\]

and generic rank 163.  Thus the maximal determinantal ideals to study are
\(I_{190}(\Phi_{\mathrm{even}})\) and
\(I_{163}(\Phi_\chi)\), together with the corresponding augmented
matrices \([\Phi_\chi\mid D_\chi]\).

Exact specialization gives the following ranks:

| parameter | even rank | nontrivial rank |
|---|---:|---:|
| generic, e.g. \((1,1)\) | 190 | 163 |
| radial \((1/5,1/10)\) | 130 | 120 |
| mixed point \((-5/3,-1/6)\) | 189 | 163 |

At the radial point every tested cube-orbit representative leaves the
image.  At the mixed point the \(x^2y^2\)-coefficient cube leaves the
image, while the \(x^4\)-coefficient cube remains in it.  The mixed point
is therefore visible in an augmented determinantal comparison even though
most cube coordinates do not detect it.

## 4. Finite-field support reconstruction

Complete scans of the four coefficient-monomial orbit representatives

\[
x^4,\qquad x^3y,\qquad x^2y^2,\qquad x^2yz
\]

were run over \(\mathbb F_p^2\).  On \(\nu\ne0\), the \(x^4,x^3y,x^2yz\)
cubes fail only at the radial point.  The \(x^2y^2\) cube fails at the
radial point and one additional point:

| \(p\) | radial point | additional point |
|---:|---:|---:|
| 11 | \((9,10)\) | \((2,9)\) |
| 13 | \((8,4)\) | \((7,2)\) |
| 17 | \((7,12)\) | \((4,14)\) |
| 19 | \((4,2)\) | \((11,3)\) |

The second column is the reduction of \((1/5,1/10)\).  Chinese
reconstruction of the third column gives unambiguously

\[
(\mu,\nu)=\left(-\frac53,-\frac16\right).        \tag{4.1}
\]

This yields the candidate parameter primes

\[
\mathfrak p_{\mathrm{rad}}=(5\mu-1,\,10\nu-1),
\qquad
\mathfrak p_{\mathrm{mix}}=(3\mu+5,\,6\nu+1).    \tag{4.2}
\]

Thus the candidate support over \(A[\nu^{-1}]\) is

\[
V(\mathfrak p_{\mathrm{rad}})
\cup V(\mathfrak p_{\mathrm{mix}}),              \tag{4.3}
\]

not Fermat plus radial.  Equation (4.3) is a stable modular
reconstruction, not yet an exact computation of
\(\sqrt{\operatorname{Fitt}_0(T)}\) or
\(\operatorname{Ass}_A(T)\).

## 5. Exact mixed fiber

At (4.1), exact rational row reduction gives

\[
\operatorname{rank}\Phi_{\mathrm{even}}=189
\]

and a nonzero remainder for the \(x^2y^2\)-coefficient cube.  The
specialized 114-equation ideal has a 157-element reduced basis and quotient
vector-space dimension 60.  Precisely the three \(x_i^2x_j^2\)
coefficient cubes survive, while the fourth powers of all fifteen
coefficients reduce to zero.

Therefore the reduced fiber is still the coefficient origin.  The point
(4.1) is a jump from cube nilpotence to fourth-power nilpotence, not a
nonzero reduced Schur quartic and not a mixed-character curve.

For geometric comparison, write
\(R=x^2+y^2+z^2\), \(P_2=x^2y^2+x^2z^2+y^2z^2\), and
\(P_3=x^2y^2z^2\).  At (4.1),

\[
h=\frac{R^3-8RP_2-32P_3}{30},
\]

and exact factorization gives the irreducible Hessian discriminant

\[
\det\operatorname{Hess}(h)=\frac1{27}\left(
3R^6+40R^4P_2-1664R^3P_3+64R^2P_2^2
+5120RP_2P_3-512P_2^3-8192P_3^2\right).       \tag{5.1}
\]

The sextic and its discriminant retain the manifest projective signed
permutation subgroup \((\mu_2^2)\rtimes S_3\), of order \(24\).  A full
automorphism-group calculation is intentionally not promoted here:
(4.1) is not a reduced exceptional Schur component.  For the same reason
there is no lower-face prolongation to run at this point.

This disproves the proposed identification

\[
\text{cube-certificate torsion support}
=\text{reduced exceptional Schur locus}.
\]

The torsion is still useful: it detects the precise nonflat nilpotent
thickening that made transformation certificates acquire denominators.

The follow-up full degree-four scans in
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json)
test all fifteen coefficient fourth powers modulo \(7,11,13\).  They
certify that every nonradial \(\mathbb F_p\)-rational parameter point on
\(D(\nu)\) has empty reduced projective Schur fiber.  In particular the
point (4.1) is reduced-empty in all three characteristics, while the
radial reduction is the unique surviving rational parameter point.

## 6. Remaining exact calculation

The integral character-block computation should now target (4.2), rather
than search for an unspecified curve:

1. compute \(K=\operatorname{modulo}(D,\Phi)\) over
   \(A[\nu^{-1}]\);
2. compute \(I_{15}(K)=\operatorname{Fitt}_0(T)\);
3. saturate its contraction by \(\nu\);
4. test equality of its radical with
   \(\mathfrak p_{\mathrm{rad}}\cap\mathfrak p_{\mathrm{mix}}\);
5. compute the primary thickness at both points and the associated primes
   of \(T\);
6. construct a separate \(\nu=0\) boundary module if a global parameter
   object including Fermat is desired.

The even character-block annihilator and function-field lift each exceeded
900 seconds.  The relation module needed for the exact zeroth Fitting ideal
was not obtained.  This leaves that ideal, its primary multiplicities, and
the associated-prime equality in step 5 open.
