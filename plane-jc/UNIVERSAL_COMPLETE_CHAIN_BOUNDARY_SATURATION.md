# Universal complete-chain boundary saturation

> **Status.**  This note proves a degree-independent purity theorem for the
> logarithmic cotangent cokernel of every resolved plane Keller map and gives
> exact nodal-tree countermodels to the stronger claim that an arbitrary
> conductor/gauge matching cokernel must be `S1` or have positive collision
> height.  It does **not** prove `JC(2)`.  The subsequent
> [`log-conductor degree-shift theorem`](LOG_CONDUCTOR_DEGREE_SHIFT.md)
> shows that the degree-zero comparison originally proposed here is not the
> natural conductor map: a branch mismatch embeds in `H_Z^1`, while the
> terminal type-I determinant becomes a unit after its full divisor is
> removed.  The correct next invariant is the nodal `Fitt_1`/localized-`c_2`
> profile of the full logarithmic differential.

The purpose of this note is to test the proposed replacement of
degree-by-degree coefficient elimination by a theorem of the form

\[
H^0_Z(M)=0
\tag{0.1}
\]

for every admissible complete chain.  The conclusion is mixed but useful:

1. there is a canonical, chain-independent module for which (0.1) is an
   unconditional theorem;
2. the matching cokernel currently proposed in the boundary atlas is not
   automatically that module;
3. even a constant nodal tree with simultaneous normalization can produce a
   matching cokernel that fails either relative height or `S1`; and
4. the universal programme cannot identify a conductor mismatch with a
   finite-support section of the logarithmic cokernel; it must retain the
   full nodal matrix and its degree-one normalization defect.

The elementary countermodels are checked by
[`verify_universal_boundary_saturation.py`](../scripts/verify_universal_boundary_saturation.py).

## 1. What the complete-chain theorem supplies

The Guccione--Guccione--Horruitiner--Valqui complete-chain theorem in
[*Some algorithms related to the Jacobian Conjecture*](https://arxiv.org/abs/1708.07936)
starts with a minimal standard pair and produces a finite sequence of Laurent
corners, directions, root choices, and transformations.  Conversely, their
algorithm enumerates combinatorial chains satisfying the displayed
geometric and arithmetic conditions.  It does not construct a resolved
proper morphism of surfaces, a simultaneous boundary normalization, or a
coherent conductor matching map.  This is explicit in the construction:
the transformations take place in the Laurent rings `L^(l)` and are driven
by roots of the successive edge polynomials.

Thus a theorem quantified over complete chains has to distinguish:

\[
\text{admissible combinatorial chain}
\quad\text{from}\quad
\text{chain realized by a map-decorated boundary package}.
\tag{1.1}
\]

Only the second object has a boundary sheaf on which local cohomology is
defined.  The distinction is also recorded in
[`NEWTON_BOUNDARY_DICTIONARY.md`](NEWTON_BOUNDARY_DICTIONARY.md).

For families, normalization does not silently commute with specialization.
One must work on strata carrying a simultaneous normalization and a flat
conductor quotient.  Constancy of the delta invariant is the standard
equinormalization signal for families of reduced curve singularities; see
Greuel--Pfister, [*The Delta Invariant and Simultaneous Normalization for
Families of Isolated Non-Normal Singularities*](https://arxiv.org/abs/2107.07012).
The conductor/contact-loss theorem in
[`CONDUCTOR_JET_TRUNCATION.md`](CONDUCTOR_JET_TRUNCATION.md) already uses
exactly this stratified formulation.

## 2. The canonical logarithmic module

Let

\[
F=(P,Q):U=\mathbb A^2\longrightarrow\mathbb A^2
\]

be a dominant Keller map over an algebraically closed characteristic-zero
field.  Choose smooth projective completions and resolve the rational map so
as to obtain a morphism of smooth surfaces

\[
f:X\longrightarrow Y
\tag{2.1}
\]

with simple-normal-crossing boundaries

\[
D_X=X\setminus U,\qquad D_Y=Y\setminus\mathbb A^2,
\qquad f^{-1}(D_Y)\subseteq D_X.
\tag{2.2}
\]

The pullback of logarithmic differentials gives a morphism of rank-two
vector bundles

\[
\theta_f:
f^*\Omega_Y^1(\log D_Y)
\longrightarrow
\Omega_X^1(\log D_X).
\tag{2.3}
\]

Define the logarithmic cotangent cokernel

\[
\boxed{\mathcal T_f^{\log}=\operatorname{coker}(\theta_f).}
\tag{2.4}
\]

Its zeroth Fitting section is the logarithmic Jacobian

\[
j_f=\det(\theta_f)
\in
H^0\!\left(
X,
\det\Omega_X^1(\log D_X)
\otimes
f^*\det\Omega_Y^1(\log D_Y)^{-1}
\right).
\tag{2.5}
\]

On `U`, this is the nonzero constant ordinary Jacobian, multiplied only by
the boundary units introduced by the logarithmic bases.  In particular,
`j_f` is not the zero section.

### Theorem 2.1 -- universal logarithmic `S1` theorem

For every resolved plane Keller map (2.1)--(2.2):

1. `theta_f` is injective;
2. `T_f^log` has the two-term locally free resolution
   \[
   0\longrightarrow f^*\Omega_Y^1(\log D_Y)
   \mathop{\longrightarrow}^{\theta_f}
   \Omega_X^1(\log D_X)
   \longrightarrow\mathcal T_f^{\log}
   \longrightarrow0;
   \tag{2.6}
   \]
3. \(\operatorname{Fitt}_0(\mathcal T_f^{\log})=(j_f)\);
4. unless it is zero, `T_f^log` is a pure one-dimensional
   Cohen--Macaulay module and hence satisfies `S1`; and
5. for every zero-dimensional closed subset \(Z\subset X\),
   \[
   \boxed{H_Z^0(\mathcal T_f^{\log})=0.}
   \tag{2.7}
   \]

#### Proof

The assertion is local on `X`.  Let `R` be a regular local ring of the
smooth surface `X`, and represent (2.3) by a `2 x 2` matrix `A`.  Its
determinant `d` is nonzero because (2.5) is nonzero on the dense open `U`.
If `Av=0`, multiplication by the adjugate gives `dv=0`.  The regular local
ring `R` is a domain, so `v=0`.  This proves injectivity and (2.6).

The presentation gives

\[
\operatorname{Fitt}_0(\operatorname{coker}A)=(\det A)=(d),
\tag{2.8}
\]

so the support is the effective Cartier divisor `V(d)`.  At a closed point
of this divisor the cokernel has projective dimension one.  Auslander--
Buchsbaum over the two-dimensional regular local ring gives depth one,
which equals the local dimension of its support.  At a generic point of the
divisor the localized support has dimension zero, so the Cohen--Macaulay
condition is automatic.  Therefore the nonzero cokernel is Cohen--Macaulay
of pure dimension one and is `S1`.

A zero-dimensional set contains no irreducible component of this pure
one-dimensional support.  The `S1` support-saturation theorem now gives
(2.7).  If the cokernel is zero, (2.7) is immediate.  \(\square\)

The proof uses neither coordinate degrees nor the length of a complete
chain.  It applies on every resolved boundary model separately; no single
Noetherian parameter space containing all chains is required.

### Corollary 2.2 -- no isolated logarithmic determinant defect

A section of `T_f^log` which vanishes away from finitely many boundary
collision points vanishes identically.  Equivalently, the logarithmic
Jacobian has only its Cartier divisorial defect; after that divisor is
accounted for, it has no residual codimension-two defect.

This is the universal boundary theorem that follows from surface geometry.
It is stronger and more canonical than asserting `S1` separately for a
different coefficient matrix attached to every degree pair.

## 3. Why an arbitrary conductor matching cokernel is different

The boundary-atlas proposal starts instead with a family of reduced boundary
curves `B -> S`, its normalization `nu:Btilde -> B`, and

\[
\mathcal Q_B=\nu_*\mathcal O_{\widetilde B}/\mathcal O_B.
\tag{3.1}
\]

After a determinant-line twist, source, target, and gauge jets give

\[
\Phi:\mathcal J\longrightarrow\mathcal Q_B\otimes\mathcal L,
\qquad M=\operatorname{coker}\Phi.
\tag{3.2}
\]

Even if `Q_B` and `J` are vector bundles on `S`, the cokernel of their
matrix can be an arbitrary finitely presented module.  The following two
examples show that the desired depth and height properties do not follow
from a nodal tree or simultaneous normalization.

### Proposition 3.1 -- `S1` does not give collision height

Let `A=k[a,b]` and let `B/A` be the constant family of a split node.  The
normalization sequence is

\[
0\longrightarrow\mathcal O_B
\longrightarrow\nu_*\mathcal O_{\widetilde B}
\longrightarrow A\longrightarrow0,
\tag{3.3}
\]

so the conductor quotient is the trivial line `Q=A`.  Take the one-jet
matching map

\[
\Phi:A\mathop{\longrightarrow}^{a}A.
\tag{3.4}
\]

Then

\[
M=A/(a).
\tag{3.5}
\]

This module is Cohen--Macaulay and `S1`, but for the collision ideal
`I=(a)` one has

\[
\operatorname{ht}_M(I)=0,
\qquad
H_I^0(M)=M\ne0.
\tag{3.6}
\]

Thus `S1` alone cannot prove the proposed vanishing.  Geometrically, the
whole matching defect is carried by the collision divisor.

### Proposition 3.2 -- a nodal tree does not force `S1`

Take a constant connected curve whose normalization is a chain of three
affine lines.  Glue the first and second lines at one marked point and the
second and third at another.  Explicitly, its coordinate ring is the
subring of

\[
A[t_0]\oplus A[t_1]\oplus A[t_2]
\]

defined by

\[
f_0(0)=f_1(0),\qquad f_1(1)=f_2(0).
\tag{3.7}
\]

The normalization quotient is `Q=A^2`, one copy for each node.  Define

\[
\Phi:A^3\longrightarrow A^2,
\qquad
\Phi=
\begin{pmatrix}
a&0&0\\
0&a&b
\end{pmatrix}.
\tag{3.8}
\]

Then

\[
M\simeq A/(a)\oplus A/(a,b).
\tag{3.9}
\]

The support has minimal prime `(a)`, while the second summand makes `(a,b)`
an embedded associated prime.  Hence `M` is not `S1`.  With
`I=(a,b)`, the class of the second conductor basis vector is nonzero and is
killed by `I`, so

\[
0\ne A/(a,b)\subseteq H_I^0(M).
\tag{3.10}
\]

The dual graph in this example is a tree, the normalization is simultaneous,
and the conductor quotient and all jet modules are free over the coefficient
base.  Therefore a proposed "tree-conductor perfection theorem" is false
without an additional restriction on the actual matching matrix.

## 4. A sufficient theorem for the original matching module

The original cokernel strategy remains valid under a checkable stronger
hypothesis.

### Proposition 4.1 -- perfect matching criterion

Let `A` be a regular Noetherian ring and let `M` be a finite `A`-module.
Assume that `M` is perfect of pure codimension `c`; equivalently, at every
prime of its support its projective dimension and support codimension are
both `c`.  Then `M` is Cohen--Macaulay and hence `S1`.  Consequently, for
an ideal `I` satisfying `ht_M(I)>=1`,

\[
H_I^0(M)=0.
\tag{4.1}
\]

#### Proof

Localize at a prime in the support.  Auslander--Buchsbaum gives

\[
\operatorname{depth}M
=\dim A-\operatorname{pd}M
=\dim A-c
=\dim\operatorname{Supp}M.
\tag{4.2}
\]

Thus `M` is Cohen--Macaulay at every localization.  It is therefore `S1`,
and (4.1) follows from positive relative height and support saturation.
\(\square\)

For a concrete presentation this criterion can be proved by a
Buchsbaum--Eisenbud, Hilbert--Burch, or Buchsbaum--Rim grade certificate.
Theorem 2.1 is its rank-two, codimension-one surface instance.  Proposition
3.2 shows that purity/expected grade is an actual hypothesis, not formal
language for the existence of a finite presentation.

## 5. Correction: conductor comparison is degree one

Let `N` be a torsion-free coherent module on a reduced boundary curve and
let

\[
 \mathcal C_N=
 \nu_*(\nu^*N/\text{torsion})/N
\tag{5.1}
\]

be its normalization defect.  The
[`log-conductor degree-shift theorem`](LOG_CONDUCTOR_DEGREE_SHIFT.md)
proves the functorial exact sequence

\[
0\longrightarrow\mathcal C_N
\longrightarrow\mathcal H_Z^1(N)
\longrightarrow
\mathcal H_Z^1\!\left(\nu_*(\nu^*N/\text{torsion})\right)
\longrightarrow0,
\tag{5.2}
\]

while `H_Z^0(N)=0`.  In particular every `O`-linear map from the
finite-support module `C_N` to `N` is zero.  A nonzero conductor mismatch
therefore has no faithful degree-zero lift of the form previously proposed.

The normalized logarithmic determinant does not repair this.  If
`Delta_f=div(j_f)` is the **complete** logarithmic different, then

\[
 \mathcal L_f(-\Delta_f)\simeq\mathcal O_X,
 \qquad j_f/s_{\Delta_f}=1
\tag{5.3}
\]

canonically, so its scalar conductor mismatch is zero on every resolved
boundary.  Two integrable local Jacobian matrices can nevertheless have the
same determinant, normalized unit, and generic branch Smith profiles while
their nodal cokernels are respectively the glued node ring and its split
normalization.  Their `Fitt_1` ideals distinguish them.

Thus a corrected complete-chain contradiction requires a class

\[
 \rho_\tau\in\mathcal C_\tau
\tag{5.4}
\]

for the normalization defect of the **full** logarithmic matrix, together
with two independent theorems: the terminal data make (5.4) nonzero, and
global Keller descent makes the same class zero.  Type-I determinant
nonvanishing proves neither assertion.

If a coefficient-base matching problem is genuinely derived rather than a
single cokernel,
the comparison must be made for the full perfect complex.  In that case the
single condition `H_Z^0(H^2)=0` must be supplemented by

\[
\operatorname{grade}(I_Z,H^q)\ge3-q,
\qquad q=0,1,2,
\tag{5.5}
\]

unless the lower cohomology modules vanish.  This is the derived boundary
criterion proved in
[`BOUNDARY_OBSTRUCTION_THEORY.md`](../extended-geometry/BOUNDARY_OBSTRUCTION_THEORY.md).

## 6. Audit of the current degree programmes

### `(72,108)`

The two coefficient systems are already empty.  The recovered lower bands
make an alternate-chart residue constructible, but no full logarithmic node
matrices, nodal `Fitt_1` profiles, normalization-defect modules, or classes
as in (5.4) have been compiled.  The row therefore supplies a regression
fixture, not evidence for a universal nodal comparison.

### `(75,125)`

The carrier and terminal calculations determine exact residue covers,
attachment points, and lower-Laurent cokernel coordinates.  They do not form
the conductor map (3.2).  In particular, the `294+53` and `7+6` blocks are
equations in a specialized coefficient circuit; they are not yet a coherent
module over a chain parameter base.  The double-carrier locus is already
zero-dimensional after the Wronskian specialization, so positive relative
height cannot be inferred by working only over that specialized base.

This is precisely the failure exhibited abstractly by Proposition 3.1.  A
universal coefficient-base module must be constructed before the carrier
specialization.  On the resolved surface, the exact terminal identity
`[P,Q]_(X,y)=X^4` becomes the normalized unit `1` after the full different is
removed, so it supplies no scalar conductor mismatch.  Any new obstruction
must come from the nodal `Fitt_1`/localized-`c_2` profile of the full matrix.

## 7. Revised research order

The degree-independent route is now:

1. **Chain-to-log realization.**  Turn a chain arising from an actual
   minimal standard pair into a resolved morphism of SNC surface pairs,
   retaining the three-section linear series and every contracted and
   dicritical component.
2. **Divisorial completeness.**  Identify the full logarithmic Jacobian
   divisor.  A residual height-one different is not a class supported at
   collision points.
3. **Full node matrices.**  Transport the complete `2 x 2` logarithmic
   differential through every boundary blowup, including the discrepancy
   factor at smooth boundary centers.
4. **Nodal Smith/Chern profile.**  Compute `Fitt_1` and the normalization
   defect (5.1) at every attachment and compare their total localized second
   Chern length with the global Chern identity.
5. **Residue and descent.**  Identify a terminal Laurent class in that
   normalization defect, prove it is nonzero, and separately prove that
   Keller geometry makes it descend.

Theorem 2.1 remains the correct purity statement for the surface cokernel,
but it does not eliminate the degree-one mismatch in step 5.  Another carrier
classification is useful only if it supplies the full node matrices needed
by steps 3--4.

<!-- status-consumer: LCDS1 5b4d92acd50d6c41 -->

## 8. Reproduction

Run the exact algebraic countermodel checks with

```bash
.venv/bin/python scripts/verify_universal_boundary_saturation.py
```

The checker verifies the presentation ideals, explicit local-cohomology
witnesses, maximal-minor ideal, and the pure-divisor control model.  The
general logarithmic purity theorem is the proof in Section 2, not a bounded
computer verification.

## References used by the proof audit

- J. A. Guccione, J. J. Guccione, R. Horruitiner, and C. Valqui,
  [*Some algorithms related to the Jacobian Conjecture*](https://arxiv.org/abs/1708.07936).
- G.-M. Greuel and G. Pfister,
  [*The Delta Invariant and Simultaneous Normalization for Families of
  Isolated Non-Normal Singularities*](https://arxiv.org/abs/2107.07012).
- The Stacks Project,
  [Auslander--Buchsbaum](https://stacks.math.columbia.edu/tag/090U),
  [Cohen--Macaulay modules](https://stacks.math.columbia.edu/tag/00N2), and
  [perfect complexes](https://stacks.math.columbia.edu/tag/0656).
