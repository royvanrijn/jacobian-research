# Two degree-three MW14 pencils and their inverse tests against 302

The [completed ten-pencil dictionary](CURVE302_SHORTWORD_TRIANGLES_2026-09-07.md)
is the active construction-recovery summary. This note and its artifacts
are retained as exact proofs and calibration inputs for the larger result.

Authority: `EC-K3-CURVE302-TRIANGLE-MW14`. This is one explicit rational
pencil on the determinant-948 K3, with arithmetic generic rank exactly 14.
It does **not** specialize to 302. The second pencil below has the same
generic rank, is inequivalent even over Qbar, and also excludes302.
No alternative parent has been recovered.

The input is the full rank-17 rootless basis of the
[11952 fibration](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md),
authority `EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION`.
Use its one-based section numbering. Put

\[
D=O+P_4+P_{11},\qquad Z=P_7.
\]

Each of the three components has square −2, and each pair intersects once.
Thus `D²=0`, `D·Z=1`, and `D` is nef: its intersection with each component
is zero and with every other irreducible curve is nonnegative. It is
primitive because of `Z`. The primitive nef genus-one pencil theorem on a
K3 therefore gives a Jacobian elliptic fibration over Q with section `Z`.
Its intersection with the old fibre is three.

## Complete rank calculation

In the old integral basis `(O,F,v1,...,v17)`, the Gram matrix is
`diag([[-2,1],[1,0]],-G)`. The section corresponding to a MW word `w` is
`(1,h(w)/2,w)`. The certificate computes the integral orthogonal complement
of `D,Z`, a positive frame of rank 17 and determinant 948.

PARI enumeration and a separate enumeration using exact rational LDL
bounds both find precisely eight signed roots. Their simple-root lattice
is `A2+A1`, of rank three and determinant six. Its Smith factors in the
frame are all one, so the root sublattice is primitive. All of the source
Néron–Severi group is rational. Consequently the new fibration has trivial
torsion and

\[
\operatorname{rank} E(\mathbb Q(s))=19-2-3=14,
\qquad \det\operatorname{MW}=948/6=158.
\]

The certificate gives an integral basis of the quotient frame/root lattice,
its representatives in the old NS basis, and its projected height Gram.
These are a complete **abstract** MW basis and lattice. Fourteen explicit
Weierstrass section coordinates have not been compiled.

## Explicit rational pencil

Write the old equation as `y²=x³+A(u)x+B(u)`. Let `P=P4`, `Q=P11`, and let
`rP,rQ,lambda` be the distinct simple parameters where `P`, `Q`, and `P−Q`
meet the old zero. For `R=P,Q`, define

\[
f_R=\frac{y+y_R}{x-x_R},\quad
a_R=(u-r_R)\frac{-y_R}{x_R},\quad
k_R=\frac{a_R(r_R)}{u-r_R}+a_R'(r_R).
\]

Then a second generator of `H⁰(D)`, besides 1, is

\[
s=c_1f_P+c_2f_Q+c_0,
\]

where

\[
c_1=\frac{\lambda-r_P}{(u-r_P)(u-\lambda)},\quad
c_2=-\frac{\lambda-r_Q}{(u-r_Q)(u-\lambda)},\quad
c_0=\frac{k_P(\lambda)-k_Q(\lambda)}{u-\lambda}-c_1k_P-c_2k_Q.
\]

At `rR`, the difference `fR−kR` vanishes along the old fibre: `xR,yR` have
orders −2,−3 and the displayed subtraction removes the pole and constant.
At `lambda`, the two sections agree and the residues cancel. At infinity,
`fR=O(u²)` and `c1,c2=O(u⁻²)`. Thus there are no vertical poles. The only
poles are simple poles on `O,P,Q`, each present. This verifies the pencil
formula over Q. Its restriction to `Z` is an exact Möbius function of `u`,
giving the rational zero section on the new base.

Eliminating `y` gives the pinned cubic in `x` over `Q(u,s)`. The complete
rational coefficients of the pencil and the trigonal equation are in the
[certificate](../artifacts/generated-results/elkies-k3-curve302-triangle-mw14-qq-v1.json).
This implicit genus-one equation is not being reported as a short
Weierstrass equation.

## Exact inverse obstruction

The frozen modular budget is three primes `1009,1013,1021`, at most 64
sample parameters and 300 seconds per prime. At each prime the source has
squarefree degree-24 discriminant. The characteristic-zero pencil reduces
exactly to the modular pencil, with its degree-one zero section retained.

For each prime, 52 distinct sample fibres normalize to genus one. Exact
Riemann–Roch computations give their Weierstrass models and j-invariants.
The K3 bound gives degree at most 24 for the j-map. Rational interpolation
from 49 samples determines it uniquely: two such maps agreeing at 49
non-poles have cross-product of degree at most 48. Three further samples
check the result. In each case the reduced map has degree exactly 24.
Therefore there is no degree loss from the characteristic-zero map, and
the resulting homogeneous inverse test includes the parameter at infinity
and rational parameters with denominators divisible by the test prime.

Modulo 1009, `j(s)=j302` has the single residue `s=37`; this is only a
necessary condition. Modulo **1013** and independently modulo **1021**,
the degree-24 homogeneous comparison has no projective root. Hence no
rational parameter on this one fibration has j-invariant equal to 302.
There is no specialization-height bound in this exclusion.

This candidate lies outside the earlier degree-two `mA1` frame census:
its frame contains `A2`, and its displayed old degree is three. This does
not establish a new surface isomorphism class, novelty in the literature,
or any exclusion of all rank-14 parents. Since the inverse obstruction
already closes this candidate, compiling its Weierstrass coordinates is
not required for the 302 search.

## A second, inequivalent triangle

Authority: `EC-K3-CURVE302-TRIANGLE-6-8-MW14`. In the same source numbering,
take

\[
D'=O+P_6+P_8,\qquad Z'=P_7.
\]

The three components again meet pairwise once, `D'^2=0`, and `D'.Z'=1`.
Thus the same effective-nef argument gives a Q-Jacobian fibration. Its
exact frame has eight signed roots of type `A2+A1`. Independent PARI and
rational LDL enumerations agree; the latter visits322 nodes. The root
lattice is primitive, so the arithmetic generic rank is exactly14,
torsion is zero and the full abstract MW determinant is158.

The preceding rational pencil formula applies with `P=P6`, `Q=P8`.
All three intersection parameters are distinct rational numbers; every
vertical-pole cancellation and the degree-one restriction to `P7` are
checked over Q. The
[new certificate](../artifacts/generated-results/elkies-k3-curve302-triangle-6-8-mw14-qq-v1.json)
contains the rational pencil, implicit trigonal equation, full abstract
MW Gram and quotient representatives. It does not supply fourteen
Weierstrass section coordinates.

Good reduction at1013 preserves this pencil and its zero section. The
[modular calculation](../artifacts/generated-results/elkies-k3-curve302-triangle-6-8-mw14-mod1013-v1.json)
uses52 genus-one fibre normalizations,49 interpolation samples and three
additional checks. The resulting j-map has degree24 and its homogeneous
comparison with302 has no projective root. Therefore **no rational
parameter on this second pencil specializes to302**. This entire modular
calculation replays in about34 seconds, followed by the characteristic-zero
checker. No parameter-height restriction is used.

## The two presentations are not duplicates

Authority: `EC-K3-CURVE302-TRIANGLE-BRANCH-SEPARATION`. For a degree24
j-map `n(s)/d(s)`, form the degree24 binary form `n_h-J*d_h`. Its binary
discriminant is a polynomial in the target coordinate J, of degree at
most46. Under a PGL2 change of the source coordinate, this polynomial
changes only by a nonzero constant. Its projective coefficient vector
therefore records an invariant of the j-map.

Both pinned reductions at1013 have degree24, and both discriminant
polynomials are nonzero. Consequently equality of their characteristic-zero
branch discriminants up to scalar would force the same equality after
reduction. This avoids assuming that a hypothetical source-coordinate
change has good reduction.

The [branch certificate](../artifacts/generated-results/elkies-k3-curve302-triangle-branch-separation-v1.json)
computes each polynomial from47 evaluations and checks three additional
values. The normalized polynomials have degree43 and factor as

\[
J^{16}(J-1728)^{12}H(J),\qquad \deg H=15,
\]

with H monic. The first pencil has `H(0)=794 mod1013`; the second has
`H(0)=97 mod1013`. They are unequal, so the two j-maps cannot be related
by PGL2 even over Qbar. In particular, the two elliptic fibrations are
inequivalent. They remain fibrations on the **same** K3 surface.

As a positive control, the checker applies `s -> (2s+3)/(s+2)` to the
first map and recovers exactly the same normalized branch polynomial.
It explicitly tests the parameter where the transformed polynomial loses
its leading term, using the homogeneous binary-discriminant convention.
The branch computation replays. It shares the two modular j-map inputs;
it is not an independent recomputation of their genus-one normalizations.

## Replay

```
sage -python elkies-k3/scripts/probe_curve302_triangle_mw14.sage --prime 1009 --check
sage -python elkies-k3/scripts/probe_curve302_triangle_mw14.sage --prime 1013 --check
sage -python elkies-k3/scripts/probe_curve302_triangle_mw14.sage --prime 1021 --check
sage -python elkies-k3/scripts/certify_curve302_triangle_mw14.sage --check
sage -python elkies-k3/scripts/probe_curve302_triangle_6_8_mw14.sage --prime 1013 --check
sage -python elkies-k3/scripts/certify_curve302_triangle_6_8_mw14.sage --check
sage -python elkies-k3/scripts/certify_curve302_triangle_branch_separation.sage --check
```

The characteristic-zero checker takes about half a second locally. Modular
replay takes about 35–45 seconds per prime. The finite-field conversion and
interpolation share implementation; only the root census has a second
algorithm. Logs and sample checkpoints are under
`artifacts/local/elkies-k3/curve302-a2-mw15/`; that historical directory name
does not assert MW15. Failed interface runs are not exclusion evidence.
The second pencil's logs are
`artifacts/local/elkies-k3/curve302-triangle-6-8-mod1013.log` and
`curve302-triangle-6-8-mod1013-replay.log` in the same parent directory;
its sample checkpoint is in `curve302-triangle-6-8-mw14/`.
These computations used SageMath10.9. The earlier QQ-only checkpoint and
the branch certificate before adding the leading-term-loss control are
preserved under the local evidence directory.
