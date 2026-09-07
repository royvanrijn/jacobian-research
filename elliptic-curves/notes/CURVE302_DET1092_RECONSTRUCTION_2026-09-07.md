# Curve 302: a uniquely recovered determinant-1092 form

Authority: `EC-CURVE302-DET1092-RECONSTRUCTION`.

**Lattice-stage certificate, now realized by the [explicit full MW17 parent](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md).**
This note and its unchanged checker retain the earlier lattice-only result.
Their no-parent boundaries describe that stage; the linked theorem supplies
the equation, full generic basis and exact302 specialization.
Removing one uniquely determined word from the existing 401 integral-point
words recovers a positive even rank-17 Gram matrix of determinant 1092 and
minimum four. Its discriminant passes a primitive K3-lattice embedding test
and the local tests for Elkies's lattice `L_546`. This gives a specific
arithmetic moduli curve to investigate. MW14–16 alternatives remain in scope;
rank17 was not imposed as the construction endpoint.

## Exact recovery from the point words

The inputs are the existing `target_core` in
[`curve302_integral_shell_inputs_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_integral_shell_inputs_v1.json).
They were selected using specialized numerical heights and exact integrality,
not this candidate Gram. The 31-by-17 integer embedding into the displayed
independent group `D` is primitive.

Write the unknown quadratic form as a symmetric matrix `G`, with one extra
unknown constant `c`. Each word `v` gives the linear equation

\[
v^{\mathsf T}Gv=c.
\]

The 401-by-154 evaluation matrix has rank154. Its left kernel has exactly one
identically zero coordinate, at zero-based row356. Thus row356 is the **only**
row whose deletion can lower the rank. After deleting it, the rank is153 and
the right kernel has dimension one, with nonzero constant coordinate.
Normalizing `c=4` gives the displayed integral Gram in the certificate. This
argument identifies the outlier without the exploratory integer optimizer.
It proves uniqueness within the single-outlier constant-norm model; it does
not prove that integrality means generic height four.

The recovered form satisfies:

- all17 diagonal entries are4;
- 400 input words have norm4, while row356 has norm6;
- determinant1092 and Smith factors `1` sixteen times, followed by1092;
- positive definiteness and no norm2 vectors, checked by exact rational LDL
  enumeration (538 nodes);
- exactly2436 signed norm4 vectors (20,688 enumeration nodes).

The discriminant form has a cyclic generator with norm
`4997/1092 = 629/1092 mod 2Z`. There is no proper even integral overlattice:
its index could only be2, but the unique order-two discriminant element has
odd norm `4997*273`. This is a statement about the abstract form, **not** a
proof of a saturated generic MW group or saturated specialization.

An exploratory mixed-height integer optimization first found the form. Its
control recovered the withheld determinant948 Gram even after143 norm6 words
were added to572 norm4 words. Those local scripts and outputs are retained
under `artifacts/local/elliptic-curves/curve302-mixed-height-*`. The exact
certificate below does not rely on a numerical optimization result or an
optimality claim.

## Additional point diagnostic

Freeze the recovered form and evaluate all1218 norm4 rays in `D` on302.
Of these,400 were used in the fit and818 were not. Exactly four of the818
additional rays give integral points; their words and exact coordinates are
in the certificate. The total is404 integral rays.

This is a limited diagnostic, not an independent section-lifting certificate.
In the known native determinant948 control, the full1313-ray norm4 shell
has572 integral rays, all already in the fit; the other741 are nonintegral.
Consequently, integrality of every predicted minimum vector is not a valid
requirement in these specialized integral-model coordinates. The target
candidate was frozen before these additional words were evaluated.

## Primitive K3-lattice embedding

Set `N=U + (-G)`, of signature `(1,18)`. A compatible ternary lattice is

\[
T=\begin{pmatrix}-2&1&0\\1&2&2\\0&2&220\end{pmatrix},
\qquad \det T=-1092,\qquad \operatorname{sig}T=(2,1).
\]

Its cyclic discriminant generator has norm `5/1092`. The identity

\[
629\cdot25^2\equiv5\pmod{2184}
\]

gives the anti-isometry between the discriminant forms of `N` and `T`.
The checker constructs the graph glue explicitly, verifies an even
unimodular rank22 Gram of signature `(3,19)`, and verifies that both `N` and
`T` embed primitively. This removes the abstract integral-embedding obstacle.
It does not construct a K3 surface over `Q` or identify a fibre through302.

The even Clifford algebra of `T` has Hilbert parameters `(5,2184/5)` and
quaternion discriminant546. These numbers are distinct from the determinant948
source; no equality of discriminants or subgroup matching with that source
is being assumed.

## A specific arithmetic moduli lead

For each odd prime `p | 546`, the positive frame's dual contains a vector
of norm `c/p` satisfying

\[
\left(\frac c p\right)=-\left(\frac{1092/p}{p}\right).
\]

The exact checker verifies the three tests at3,7,13. Together with signature
and determinant, these are Elkies's characterization of `L_546` (expressed
using the positive essential lattice). His moduli description associates
`L_546`-marked K3 surfaces with `X(546)/<w546>`.
See [Elkies, section2, pp6–7](https://arxiv.org/pdf/0802.1301).

The quotient is the elliptic curve546c2, with a published model

\[
C:\quad y^2+xy+y=x^3-137x+380.
\]

See [Rotger, thesis Table6.4](https://web.mat.upc.edu/victor.rotger/docs/Tesi.pdf).
The identification also agrees with the local Cremona database. It is a
literature input, not an independently reconstructed period map.

The point `P=(-9,-26)` has infinite order. Good reductions at5,11,17 have
orders8,16,20, whose gcd is4; the checker verifies `4P != O`. Thus `C(Q)`
is infinite. Rational CM points on this fixed Shimura quotient are finite,
so its rational points cannot all be CM. **The CM status of the displayed
point P itself remains unknown here.** No universal marked K3 equation or
explicit non-CM marked specialization has been recovered.

The useful next construction problem is now specific: recover the K3 family
and its rational marking over this elliptic moduli curve, transport the
recovered `U` and section lattice, and solve for a fibre Q-isomorphic to302.
An equation for the moduli curve is not an equation for the desired elliptic
surface. A rational point on it is not yet the required specialization `t0`.
This candidate does not establish the discoverers' provenance.

## Replay

```
OPENBLAS_NUM_THREADS=1 sage -python elliptic-curves/cas/verify_curve302_det1092_recovery.sage
```

The [certificate](../../artifacts/generated-results/elliptic-curves/curve302_det1092_recovery_v1.json)
pins the inputs and checker. Replay verifies all401 original point words,
linear-algebra uniqueness, the two complete finite shells, the four new
integral points, primitive embeddings, discriminant forms and gluing, local
`L_546` conditions, and infinite order on the cited moduli curve. The shell
routine is independent of the exploratory PARI enumeration. The literature
moduli identification is not independently replayed. The initial point
selection remains finite and numerically bounded.
