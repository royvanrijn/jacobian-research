# Axial B, zero cubic contraction and a linear eigenvector

17 November2026. Claim `SIC2B33AXIALCOVEIG20261117`.
The **entire cyclic-vector rank-at-most-one stratum** below is excluded
by a finite purely algebraic certificate. The complementary chart's
unconditional first-eight ordinary unit remains uncomputed. Its
all-order pure points are now excluded by the stronger
[top-two-boundary-jet proof](CUBIC_TOP_TWO_JET_GATE_AND_FULL_AXIAL_CONTRACTION.md).

Let \(B\ne0\) be semisimple, trace-free and self-adjoint, with a
repeated eigenvalue. Keep all linear and cubic coefficients and impose
only the declared conditions

\[
 \mathcal C=\Delta^2(H_2H_3)=0,\qquad v\wedge Bv=0,
 \quad H_1=\langle v,x\rangle.
\]

**Theorem.** These phases have no common zero of their first eight
moments over a characteristic-zero field. In particular none is pure.
No minimal prefix is asserted. The new exclusions are certified by
directly multiplied rational polynomial identities, without a
standard-basis oracle. There is no msolve dependency.

A consequence is a closed pure-safe container with arbitrary linear
coefficients: require that \(B\) lie in the closed harmonic-square
quadratic cone, and impose the two equations above. Every pure point is
one-sided, so the full original-multiplier cutoff is
\(m>\deg_{Z,Y}Q\). The container is not entirely pure and is not
asserted to contain every pure phase. The
[machine-readable equations](axial_contraction_zero_eigenvector_safe_container.json)
drop the old \(H_1=0\) restrictions and add the actual eigenvector and
contraction equations.

Run [the checker](verify_axial_contraction_zero_eigenvectors.py). Its
[receipt](axial_contraction_zero_eigenvector_verification.json) records
the exact evidence. Imported cone and boundary classifications retain
their existing assurance scopes; no independent specialist, external or
formal verification is claimed.

## 1. Complete field normalization

The full axial quadratic stratum normalizes to a nonzero multiple of
\(H_2=t^2-zw\). In this frame the operator has eigenvalues
\((-1/2,1,-1/2)\). Direct differentiation gives

\[
 \Delta^2((t^2-zw)H_3)=24\gamma_{-1}w+72\gamma_0t+24\gamma_1z.
\]

Thus, after the contraction condition, the **entire** remaining phase is

\[
 F=az+bt+cw+(t^2-zw)+dw^3+ew^2t+fz^2t+gz^3.       \tag{1}
\]

All seven displayed coefficients stay free before the further field
conditions are used. The [exact generator](axial_contraction_zero.py)
integrates each charge-balanced spatial product. It retains the complete
[moment input through twelve](axial_contraction_zero_prefix12.json),
with independent direct spatial checks through six.

The dual vector is \(v=(c,b,a)\). Its wedge with \(Bv\) is
\((-3ab/2,0,3bc/2)\). Over a field, the eigenvector condition is therefore
exactly the union of \(a=c=0\) and \(b=0\), including \(v=0\).
Equivalently \(\operatorname{span}(v,Bv,B^2v)\) has dimension at most
one. No unavailable independent rotation of coefficient channels is used.

The first equations for (1) are

\[
\begin{aligned}
 \mu_2&=(2ac+b^2)/3+1/5+4ef/105+4dg/35,\\
 \mu_3&=2(b^2-ac)/5+2/35-4dg/35.
\end{aligned}                                                    \tag{2}
\]

## 2. The simple eigenaxis: a=c=0

Equations (2) solve the **products**, without dividing by individual
cubic coefficients:

\[
 ef=-(27+77b^2)/4,\qquad dg=(1+7b^2)/2.
\]

For the next two moments every surviving cubic monomial is a product
of \(dg\) and \(ef\). Substitution gives the primitive constraints

\[
 -3227b^4+1314b^2+297=0,\qquad
 -173b^4+306b^2+27=0.
\]

They are coprime. The [exact Bezout row](axial_contraction_zero_central_unit.json)
is directly multiplied, so this entire eigenaxis case is excluded by
moment five, including \(b=0\) and every vanishing cubic coefficient.

## 3. The repeated eigenspace: b=0

If \(a=c=0\), the previous case applies. Otherwise an actual oriented
swap puts \(c\ne0\).

If \(ac=0\), this gives \(a=0\). Use amplitude and the stabilizing
torus to put both the quadratic coefficient and \(c\) equal to one.
Equations (2) give \(dg=1/2\), \(ef=-27/4\), so \(d,e\) are units.
Put \(g=1/(2d)\), \(f=-27/(4e)\); clearing only these declared units
gives primitive residuals

\[
\begin{aligned}
 R_4&=3159d^2+198de^2+286e^2,\\
 R_5&=-729d^2+27de^2-143e^2,\\
 R_6&=-166212d^3+22464d^2e^2+1020357d^2-12920de^2+5168e^4.
\end{aligned}
\]

The [six exact multiplier entries](axial_contraction_zero_transverse_isotropic_witness.log)
prove both \(e^2,d^2\in(R_4,R_5,R_6)\) by direct multiplication.
This contradicts their being units. Hence moment six excludes the entire
isotropic repeated-eigenspace linear case.

If \(ac\ne0\), use amplitude and the same torus to put \(a=c=1\).
This leaves the quadratic coefficient \(\lambda\ne0\) variable; it
must not also be silently normalized. Set \(u=\lambda^2\). Odd moments
are divisible by \(\lambda\), because the remaining terms are odd under
the antipodal map. Divide them by this declared unit. The first equations
now give

\[
 u=7+2dg,\qquad 217+54dg+4ef=0.
\]

The remaining moments4..8 are exact polynomials in the four unrestricted
cubic coefficients. The
[input](axial_contraction_zero_linear_normalized_charts.json) and
[six rational multiplier rows](axial_contraction_zero_transverse_nonisotropic_unit_witness.log)
give

\[
 1=A_0(217+54dg+4ef)+\sum_{n=4}^{8}A_nR_n.
\]

The multipliers have term counts \((138,51,51,29,29,3)\). The checker
parses and multiplies them exactly with FLINT rational polynomials.
No cubic coefficient is inverted in this case. This excludes the
entire nonisotropic repeated-eigenspace case by moment eight and
completes the eigenvector stratum.

## 4. The closed safe-cone consequence

The earlier [harmonic-square cone](AXIAL_QUADRATIC_ZERO_LINEAR_EXCLUSION.md)
has a finite lifted presentation
\(\operatorname{rank}(B+\eta I)\le1\), with \(\operatorname{tr}B=0\),
and seven cubic equations in the quadratic coefficients. Its nonzero
non-nilpotent points are exactly semisimple axial; at \(\eta=0\) one
has rank-at-most-one nilpotent \(B\).

Intersect that closed cone with \(v\wedge Bv=0\) and \(\mathcal C=0\).
The semisimple axial part has no pure points by the theorem. At rank-one
nilpotent \(B\), the rank-one part of the
[nilpotent contraction-zero classification](NILPOTENT_CONTRACTION_ZERO_AND_THREE_GLOBAL_TARGETS.md)
gives a one-sided pure phase, with arbitrary linear harmonic retained.
At \(B=0\), the all-order odd theorem applies. Hence every pure point
of this entire closed container is one-sided and full-Image safe with
cutoff coefficient one. The existing odd sharp examples lie at \(B=0\),
so the coefficient-one bound is sharp on the container.

The [exported container](axial_contraction_zero_eigenvector_safe_container.json)
uses seven quadratic-cone equations, three cubic eigenvector equations
and three quadratic contraction equations in the full fifteen-variable
harmonic coefficient space. It imposes no zero-linear condition and no
extra restriction on cubic coefficients. No ordinary coefficient-ring
safety is asserted from this field-valued proof.

## 5. Complete input for the remaining mixed-eigenspace case

For an unconditional finite-prefix axial contraction-zero exclusion,
the remaining linear case is \(b\ne0\) with a nonzero transverse part, equivalently cyclic
rank two. Choose the swap so \(c\ne0\), then use amplitude and the
torus to set \(b=c=1\), retaining the nonzero quadratic amplitude
\(\lambda\). Equations \(\mu_2=0\), \(\mu_3/\lambda=0\) give

\[
 a=(112-54dg-4ef)/217,\qquad
 \lambda^2=(8dg-4ef-105)/31.
\]

This is a complete four-variable input in \((d,e,f,g)\). It retains
every value of \(a\), including isotropic linear vectors and degeneracy
of the restricted form on the rank-two cyclic space. The
[normalized chart data](axial_contraction_zero_linear_normalized_charts.json)
and exact source generator are preserved. A unit in this entire
four-variable input would exclude the remaining chart, without needing
to invert the displayed squared-amplitude polynomial.

Finite-field first-eight calculations return a unit. Rational std and
slimgb attempts stopped at45/35 seconds. The earlier five-variable
normalization, its channel block order, the one-parameter generic
calculation and modular reconstruction also stopped without rational
certificates. All source/output/run receipts are preserved. None is a
characteristic-zero exclusion. The verified eigenvector theorem does
not depend on these scouts.

For the homogeneous mixed input, the modular first unit power is
\(h^{20}\). This is **not** a rational witness-degree bound. A
Hilbert-guided QQ run reached degree eighteen before its150-second cap.
A separate degree-at-most-seventeen calculation completed and retained
190 exact rational polynomials in
`axial_contraction_zero_mixed_linear_homogeneous_degree17.log`.
This is an explicit reusable **partial basis**, not a complete Groebner
basis or a rational unit proof. Subsequent source files define a new
ideal from these polynomials, so a stale `isSB` flag is not reused.
Compressed source/output checkpoints retain the larger coefficient data.

The degree-eighteen resume from the saved basis stopped at its
ten-minute cap without a completed checkpoint. Its compressed source,
progress output and run receipt remain preserved. Macaulay21.22's
experimental `LinearAlgebra` algorithm failed immediately on the affine
QQ input with an array-order error and on the homogeneous QQ input with
a segmentation fault. These are implementation failures and provide no
mathematical exclusion. Do not restart these unchanged calculations.

Acceptance of any resumed full homogeneous certificate requires exact
standard-basis verification, exact reduction of the original input,
the proposed pure homogenizer power, and matching leading ideals with
the retained modular original ideal. Degreewise integer matrix ranks
then prove ideal equality over \(\mathbb Q\). Until those checks pass,
the mixed chart's ordinary first-eight unit remains uncomputed.
The later top-two-jet proof excludes its all-order pure points by a
different method; the failed large computations remain failures.

The three global radical memberships and the full(3,3) Image theorem
also remain OPEN. The overarching goal stays ACTIVE.
