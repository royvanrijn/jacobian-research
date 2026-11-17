# The complete nilpotent quadratic boundary with zero linear harmonic

17 November 2026. **SIC2B33NILPQUAD20261117.** The field-valued
two-pair(3,3) Image conclusion holds on the complete closed container

\[
 H_1=0,\qquad B\text{ nilpotent},\qquad H_3\text{ unrestricted}.
\]

Here “holds on the container” means every all-order pure point is safe
for every original multiplier. It does not mean every point is pure.
The cutoff is \(m>\deg_{Z,Y}Q\).

The new rank-two statement is more explicit. If \(B\ne0\) has rank two
and is nilpotent, choose an oriented linear isometry normalizing
\(H_2=zt\). Then

\[
\begin{split}
F\text{ pure}
&\iff \mu_1=\cdots=\mu_8=0\\
&\iff F=zt+a z(t^2-zw/2)+b z^2t+c z^3.
\end{split}
\]

Every surviving term has positive charge, so these phases are in the
SO3 nullcone. Eight is the least consecutive **set-theoretic** prefix
on the nonzero rank-two nilpotent stratum. Ordinary moment-ideal
stabilization is not asserted.

The largest exclusion below retains an **exact rational Singular
dependency**. Two other exclusions have directly checked polynomial unit
witnesses. Direct witness exports for the largest chart timed out in
three bounded attempts; that stronger audit is incomplete. There is no
msolve dependency or independent specialist/formal verification.

## Complete normalization

For a nonzero rank-two nilpotent self-adjoint operator, \(B^3=0\) and
\(B^2\ne0\). Choose \(v\) with
\(c=\langle v,B^2v\rangle\ne0\); this is possible since the nonzero
self-adjoint operator \(B^2\) defines a nonzero quadratic form. In the
Jordan chain \((B^2v,Bv,v)\), the Gram matrix is

\[
 \begin{pmatrix}0&0&c\\0&c&d\\c&d&e\end{pmatrix}.
\]

Replacing \(v\) by \(v+\alpha Bv+\beta B^2v\), with

\[
 \alpha=-\frac d{2c},\qquad
 \beta=-\frac{e+2\alpha d+\alpha^2c}{2c},
\]

reduces it to \(c\begin{psmallmatrix}0&0&1\\0&1&0\\1&0&0\end{psmallmatrix}\).
Scaling the chain and applying a diagonal isometry normalizes the
quadratic to \(zt\). Orientation can be corrected by \(-\mathrm{Id}\),
which fixes the quadratic. The only divided quantity is the chosen
nonzero \(c\). Work over an algebraic closure and descend the resulting
polynomial vanishing to the original characteristic-zero field.

Retain all seven cubic coefficients:

\[
\begin{aligned}
H_3={}&g_{-3}w^3+g_{-2}w^2t+g_{-1}w(t^2-zw/2)\\
 &+g_0(t^3-3ztw)+g_1z(t^2-zw/2)+g_2z^2t+g_3z^3.
\end{aligned}
\]

The torus \((z,t,w)\mapsto(\rho z,t,\rho^{-1}w)\), combined with phase
amplitude \(\rho^{-1}\), fixes \(zt\) and sends
\(g_j\mapsto\rho^{j-1}g_j\). If a nonpositive-charge coefficient is
nonzero, choose its lowest index. The four exhaustive cases are
\(-3,-2,-1,0\). Its nonzero coefficient can be normalized to one,
since \(j-1\ne0\). There is no extra stabilizer or omitted cubic
coefficient.

The primitive first two nontrivial moments are

\[
 M_2=15g_0^2+5g_1g_{-1}+4g_2g_{-2}+12g_3g_{-3},\qquad
 M_3=g_0g_{-1}+g_1g_{-2}+2g_2g_{-3}.
\]

If the lowest index is zero, \(M_2=15g_0^2\ne0\), excluding that case.

## Three explicit elimination charts

The remaining charts use only constant-pivot eliminations of \(M_2,M_3\).

| Lowest cubic index | Normalization and elimination | Free coordinates | Excluding prefix |
| --- | --- | --- | --- |
| \(-3\) | \(g_{-3}=1\), \(g_2=-(g_0g_{-1}+g_1g_{-2})/2\), \(g_3=-(15g_0^2+5g_1g_{-1}-2(g_0g_{-1}+g_1g_{-2})g_{-2})/12\) | \(g_0,g_1,g_{-1},g_{-2}\) | 8 |
| \(-2\) | \(g_{-3}=0,g_{-2}=1\), \(g_1=-g_0g_{-1}\), \(g_2=-(15g_0^2-5g_0g_{-1}^2)/4\) | \(g_0,g_3,g_{-1}\) | 7 |
| \(-1\) | \(g_{-3}=g_{-2}=0,g_{-1}=1\), \(g_0=g_1=0\) | \(g_2,g_3\) | 8 |

The [machine-readable certificate](boundary_nilpotent_quadratic_certificate.json)
contains every primitive residual polynomial, the exact coordinate
substitutions and the available polynomial unit multipliers. The checker
regenerates them from the complete phase, rather than trusting stored
elimination input.

On the lowest-\(-1\) chart the fourth and eighth residuals alone give
the small identity

\[
 \frac{6g_3+171}{6991}(39-2g_3)
 +\frac2{6991}(6g_3^2+54g_3+161)=1.
\]

The lowest-\(-2\) unit identity uses four multipliers with term counts
\((26,31,10,13)\); their full rational coefficients are retained and
checked by polynomial multiplication.

For the four-variable lowest-\(-3\) chart, an exact characteristic-zero
standard-basis computation on the residuals of moments4 through8 gives
the basis \((1)\). Its complete regenerated input and output are
[retained](boundary_nilpotent_largest_chart_verified_input.sing), with
[output](boundary_nilpotent_largest_chart_verified_output.log) and hashes
in the [verification receipt](boundary_nilpotent_quadratic_verification.json).
The receipt identifies Singular4.3.2 and the 30-second replay bound.
This last exclusion currently depends on that exact CAS calculation;
the note does not claim an independently expanded unit witness for it.

## Minimal prefix, marked safety and the closed container

The exact phase

\[
 F=zt+w(t^2-zw/2)+\frac{39}{2}z^3
\]

has moments1 through7 zero, but

\[
 h(F^8)=\frac{111856}{9561123}\ne0.
\]

This is independently checked by spatial integration. Thus eight is
necessary. The exhaustive exclusions above make it sufficient: only
positive cubic charges remain, and these are pure at every order by
torus balance.

Each factor of a surviving phase has charge at least one. A spherical
marker of degree at most \(b\) has charge at least \(-b\), so every
marked moment is zero for \(m>b\). The retained Hopf coefficient-balancing
identity handles every coefficient of the original polynomial Image
output, giving the cutoff in the statement. For example,
\(h((zt+z^3)^2w^2)=2/105\ne0\), so the strict uniform cutoff is attained.

For a trace-free three-by-three operator, nilpotence is equivalent to
\(\operatorname{tr}(B^2)=\det B=0\). Nonzero nilpotent ranks are exactly
one and two. Combine this result with the
[complete rank-one classification](RANK_ONE_QUADRATIC_BOUNDARY.md)
and the [odd theorem](ODD_HOPF_CUBIC_IMAGE_TRANSFER.md) at \(B=0\).
This proves all-original-multiplier safety at every pure field point of
the full closed container \(H_1=0,B\text{ nilpotent}\).

For nonzero nilpotent \(B\), the maximum of the two proved consecutive
prefixes is10; rank-one false positives show that10 is minimal across
those nonzero ranks. No prefix assertion for \(B=0\) is added here.

## Replay and limits

Run `python verify_boundary_nilpotent_quadratic.py` with SymPy and Singular.
The [moment engine](boundary_nilpotent_quadratic.py) uses exact torus
balance and a finite binomial integration formula. Its moments are
cross-checked by direct integration through6; the first-seven false
positive is checked through8. The checker verifies the full Jordan-chain
Gram adapter, torus normalization, constant-pivot eliminations, both
available unit identities, and the exact largest-chart CAS calculation.

The failed direct-lift and tracked-lift inputs, logs and run receipts
are preserved alongside the successful inputs. Their timeouts neither
invalidate the retained standard-basis computation nor constitute a
complete independent unit-witness audit.

The other zero-linear quadratic strata, nonzero-linear \(D=0\) strata
and unrestricted cyclic quotient remain open. This theorem does not
complete the [overarching goal](IMAGE_33_OVERARCHING_GOAL.md).
