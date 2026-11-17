# Zero linear harmonic: complete semisimple axial exclusion

17 November 2026. **SIC2B33AXIALQUAD20261117.** No all-order pure phase
with \(H_1=0\) can have a nonzero semisimple trace-free quadratic
operator with a repeated eigenvalue. The harmonic cubic is unrestricted.
The first ten moments already exclude this entire stratum. Minimality
of this prefix and ordinary ideal stabilization are not claimed.

The proof retains an exact rational Singular dependency on its largest
chart. A separate three-variable chart has a directly checked polynomial
unit identity; the smaller charts have constant exclusions. No
finite-field calculation is used to establish characteristic-zero
exhaustion. A 45-second tracked-witness export on the largest chart timed
out and remains preserved; its independent polynomial-multiplier audit
is incomplete.

## Complete normalization and charts

The nonzero repeated-eigenvalue spectrum is \((\lambda,\lambda,-2\lambda)\)
with \(\lambda\ne0\). Semisimplicity and self-adjointness make the
eigenspaces orthogonal and nondegenerate. Over an algebraic closure an
oriented isometry and phase scaling put the quadratic in the form

\[
 H_2=t^2-zw,\qquad q=2zw+t^2.
\]

Retain the full cubic, with the same seven harmonic basis vectors as
the [rank-one boundary proof](RANK_ONE_QUADRATIC_BOUNDARY.md):

\[
 H_3=\sum_{j=-3}^{3}g_j C_j,
\]

where \(C_j\) has charge \(j\). The torus preserves \(H_2\). The
oriented reversal \((z,t,w)\mapsto(w,-t,z)\) also preserves \(H_2\)
and exchanges charges \(j,-j\), with a nonzero sign. Thus choose the
largest absolute nonzero cubic charge \(J\), reverse if necessary, and
normalize \(g_J=1\) by a torus root. This uses no division by a
coefficient not already declared nonzero.

There are exactly four cases: \(J=3,2,1,0\). In the \(J=3\) chart all
six remaining cubic coefficients remain free before the moment equations.
In the lower charts, both coefficients of larger absolute charge vanish
by the definition of \(J\). The zero cubic is included in the central
case. No independent channel symmetry or hidden coefficient slice is used.

The primitive moment polynomials are

\[
\begin{aligned}
M_2&=15g_0^2+5g_1g_{-1}+4g_2g_{-2}+12g_3g_{-3}+21,\\
M_3&=4g_0^2+g_1g_{-1}-4g_3g_{-3}+2.
\end{aligned}
\]

They differ from the raw moments only by nonzero rational factors.

If \(J=0\), \((4M_2-15M_3)/54=1\), excluding the entire central
family already at moment three.

## Constant-pivot elimination

For the other cases put \(y=g_0\). The following substitutions, together
with the displayed quadratic relation, are exactly equivalent to
\(M_2=M_3=0\) on their charts:

| Chart | Substitutions | Retained relation |
| --- | --- | --- |
| \(J=3\) | \(g_3=1,\quad g_{-3}=-\frac12-\frac5{108}g_1g_{-1}-\frac4{27}g_2g_{-2}\) | \(27y^2+8g_1g_{-1}+4g_2g_{-2}+27=0\) |
| \(J=2\) | \(g_{\pm3}=0,\quad g_2=1,\quad g_{-2}=-(5g_1g_{-1}+54)/16\) | \(4y^2+g_1g_{-1}+2=0\) |
| \(J=1\) | \(g_{\pm3}=g_{\pm2}=0,\quad g_1=1,\quad g_{-1}=-54/5\) | \(5y^2-11=0\) |

Every pivot is a nonzero rational constant. Dividing the retained relation
by its constant leading coefficient makes it monic in \(y^2\). Reduce
all later moments by that relation, retaining the relation in the ideal.

In the \(J=1\) case the primitive reduced fourth moment is exactly one,
so no point survives.

The \(J=2\) certificate gives an explicit rational unit combination of
the reduced fourth through seventh moments. Its coefficient polynomials
depend only on \(g_1g_{-1}\), have at most three nonzero terms each,
and are checked by direct polynomial multiplication. The retained
quadratic relation and moments eight through ten have zero multipliers
in this identity.

On \(J=3\) the remaining variables are
\((y,g_1,g_2,g_{-1},g_{-2})\). The exact rational standard basis of the
retained quadratic relation and reduced moments four through ten is
\((1)\). The complete regenerated
[input](boundary_axial_largest_chart_verified_input.sing),
[output](boundary_axial_largest_chart_verified_output.log), and their
hashes are retained in the
[verification receipt](boundary_axial_quadratic_verification.json).
The successful replay uses Singular4.3.2 and takes about32 seconds for
this chart; its replay limit is60 seconds. This is the explicit exact
CAS dependency of the exclusion.

Together these four cases exhaust the complete family.

## A closed pure-safe cone

The result has a useful closed formulation. Let \(Y_{\rm sq}\) be the
quadratic cone obtained by harmonic projection of squares of linear forms:

\[
 H_2(x)=\langle n,x\rangle^2-\frac{q(n)}3q(x).
\]

Equivalently, introduce \(\lambda\) and impose
\(\operatorname{rank}(B+\lambda I)\le1\). Trace zero gives

\[
 \lambda^2=\frac{\operatorname{tr}(B^2)}6,
\]

so the parameter is integral and the projection to the quadratic
coefficient space is closed. Exact elimination gives **seven cubic
equations** for this cone; all are exported in
[the safe-container input](boundary_axial_safe_container.json).
The small rational elimination is replayed separately in SymPy.

At \(\lambda\ne0\), rank one and trace \(3\lambda\) give
\((B+\lambda I)^2=3\lambda(B+\lambda I)\). Thus \(B\) is semisimple
with eigenvalues \((2\lambda,-\lambda,-\lambda)\), and the theorem
excludes every pure point when \(H_1=0\).

At \(\lambda=0\), \(B\) has rank at most one and is nilpotent. Its pure
points are already covered by the
[rank-one theorem](RANK_ONE_QUADRATIC_BOUNDARY.md), with the
[odd theorem](ODD_HOPF_CUBIC_IMAGE_TRANSFER.md) at \(B=0\). Hence every
pure field point in the **closed** container
\(H_1=0,H_2\in Y_{\rm sq}\), with unrestricted \(H_3\), is safe for
every original multiplier, with cutoff \(m>\deg_{Z,Y}Q\).

This container need not consist entirely of pure points. The semisimple
exclusion does not establish a coefficient-ring theorem or a fixed
nilpotence exponent for \(\lambda\).

## Replay and remaining cases

Run `python verify_boundary_axial_quadratic.py`. The
[exact moment engine](boundary_axial_quadratic.py) is independently
cross-checked by spatial integration through moment six. The checker
verifies every torus/reversal adapter, constant-pivot elimination,
quadratic remainder and exported unit witness, replays the largest exact
rational unit calculation, and checks all seven closed-cone equations.
The [machine-readable certificate](boundary_axial_quadratic_certificate.json)
retains every residual polynomial and the small unit identity.

The distinct-eigenvalue rank-two attempt, normalized to
\(H_2=z^2+w^2\), remains **OPEN in characteristic zero**. Its exact
input through12 is retained, together with rational std/slimgb/sba
timeouts and an integer-strategy run stopped at120 seconds. Finite-field
units and a degree16 homogeneous modular membership do not prove the
rational exclusion; all tracking timeouts are retained. These results
must not be promoted to a theorem.

Other distinct-eigenvalue quadratic shapes, the nonsemisimple repeated
eigenvalue case, and all remaining nonzero-linear strata are also open.
The [overarching(3,3) goal](IMAGE_33_OVERARCHING_GOAL.md) stays active.
There is no msolve dependency, independent specialist review, external
review or formal verification for this result.
