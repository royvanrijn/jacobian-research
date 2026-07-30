# Open problems for extending the map constructions

This is a working list of theorem targets suggested by the repository's
current map constructions.  It is not a status ledger: verified claims and
their audit levels remain recorded in
[`MATH_STATUS.json`](MATH_STATUS.json) and [`STATUS.md`](STATUS.md).

The problems below are useful filters for proposed extensions.  A new family
is most valuable when it supplies one of the missing modules, obstruction
classes, support reductions, or saturation certificates described here,
rather than only adding another bounded coefficient calculation.

## 1. Resolved GMC(2) module and retained design lesson

For circular coordinates \(Z,W\), put \(U=ZW\) and consider

\[
 P=W A(U)+C(U)+ZB(U),\qquad D(U)=U A(U)B(U),
\]

with all three rotational levels nonzero.  Constant-term extraction reduces
the vanishing of every pure Gaussian moment to

\[
 \mathcal L\!\left(e^{tC(U)}
 I_0\!\left(2t\sqrt{D(U)}\right)\right)=1,
 \qquad \mathcal L(U^j)=j!.
\]

The former primary theorem target was:

> **Three-level rigidity.**  If the displayed formal identity holds for
> polynomials \(C,D\), then \(C=D=0\).

This target is now proved in every degree by the prime-endpoint theorem, and
the lower-face prime theorem proves the stronger arbitrary-support statement
and hence all of GMC(2).  The meromorphic Pfaffian system

\[
 G(t,U)=\bigl((1-tC)^2-4t^2D\bigr)^{-1/2}.
\]

is therefore retained as a structural interpretation, not an open proof
route.  The reusable design lesson is to expose a finite weight decomposition,
take the lower radial-order Newton face over weight zero, and combine a
nonzero face constant term with prime dilation and factorial divisibility.
This replaces support-tree leaf removal and automatically includes mixed
Hilbert-basis relations.

See
[`extended-geometry/GMC2_RESEARCH_PROGRAM.md`](extended-geometry/GMC2_RESEARCH_PROGRAM.md)
and
[`extended-geometry/TWO_REAL_GMC_SUPPORT_GRAPH_EXPLORATION.md`](extended-geometry/TWO_REAL_GMC_SUPPORT_GRAPH_EXPLORATION.md).

## 2. Minimal cubic Keller classification: cotangent saturation

For a finite-normalization diagram \(A\to B\), let
\(\Omega_{B/A}\) have a finite free presentation

\[
 F_1\longrightarrow F_0\longrightarrow\Omega_{B/A}\longrightarrow0,
\]

write \(N=\operatorname{im}(F_1\to F_0)\), and let \(I\) be the collision
ideal used in the cubic frontend.  After the canonical-bidual defect
\(C/T\) vanishes, the remaining closed-point obstruction is

\[
 H_Z^0(\Omega_{B/A})
 \simeq (N:I^\infty)/N.
\]

The direct theorem target is therefore

\[
 \boxed{N:I^\infty=N.}
\]

This is weaker than proving that the entire normalization is flat.  By the
support-saturation principle, it is enough to prove any one of:

- no associated prime of \(\Omega_{B/A}\) contains \(I\);
- \(\operatorname{grade}(I,\Omega_{B/A})\geq1\);
- \(I\) contains an \(\Omega_{B/A}\)-regular element; or
- the displayed presentation is already \(I\)-saturated.

A useful extension of the cubic maps should preserve a primitive conormal
generator in codimension one while making one of these depth or
associated-prime conditions transparent.  This may close intrinsic
curvilinearity without constructing a complete Cartier boundary atlas.

See
[`cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md`](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md),
[`cancellation/CUBIC_CLOSURE_ATTACKS.md`](cancellation/CUBIC_CLOSURE_ATTACKS.md),
and
[`verified/SUPPORT_SATURATION_PRINCIPLE.md`](verified/SUPPORT_SATURATION_PRINCIPLE.md).

## 3. Smallest false Dixmier rank and the `DC_2` obstruction program

The inverse-Jacobian construction already gives a non-surjective endomorphism
of `A_3`, hence of `A_n` for every `n>=3`.  Connecting the counterexample to
Dixmier is therefore complete.  The dimensional problem is to determine the
smallest false `DC_n`, with `DC_2` as the immediate target; `DC_1` remains a
separate question.

The rank-two constructions produce explicit classical symplectic symbols
and symbol-specific failures of filtered Weyl quantization.  At correction
order \(k\), the recursive equation has the form

\[
 d_G(\Delta_k)=-\mathcal O_k.
\]

The natural invariant target is

\[
 [\mathcal O_k]\in\operatorname{coker}(d_G),
\]

or equivalently a cohomology class in the filtered deformation complex or
mapping cone controlling lifts of the fixed classical map.

The unrestricted Hochschild or Poisson cohomology of the ambient affine
symplectic/Weyl algebra is not the right receptacle: the obstruction survives
because corrections are restricted by filtration, parity, degree, pole
order, and prescribed principal symbol.  The immediate theorem target is:

> Convert the computed residue functionals into dual cocycles of the
> restricted deformation complex, and prove that their pairing with
> \(\mathcal O_k\) is invariant under every allowed lower-order gauge change.

This target is now achieved at two certified fibers of the parity-restricted
complex.  The 16-term all-pole \(\hbar^5\) functional is a dual \(2\)-cocycle
and is constant on a \(1075\)-dimensional affine superset of all allowed
finite-Laurent lower lifts at \((\kappa,\tau)=(0,1)\).  On the explicit cubic
fifth-order branch, \([X^{18}]\) is a one-term \(\hbar^7\) dual cocycle
constant on the complete 20-dimensional fifth-order correction torsor.
The coherent kernel \(\mathcal P_n\) which contains these strong cocycles,
and the distinction between fiberwise, vertical, and parameter-uniform
classes, are defined in
[`extended-geometry/RESTRICTED_QUANTUM_DEFORMATION_COCYCLES.md`](extended-geometry/RESTRICTED_QUANTUM_DEFORMATION_COCYCLES.md).
What remains on this symbol-specific branch is to globalize these fiberwise
sections over its parameter and Fitting strata and to treat order seven over
the full reduced fifth-order lift component.

The first globalization chart is now explicit.  On the standard
16-monomial support, a primitive \(15\times16\) polynomial presentation has
two maximal minors with coprime nonboundary factors of degrees 34 and 35.
Thus its rank-drop locus has no divisorial component.  All sixteen maximal
minors give a zero-dimensional length-218 scheme over each of three good
finite fields.  The three saturated bases also have the same 21-generator
leading ideal, with saturation exponent 12; its staircase contains exactly
218 standard monomials.  Those 21 rational basis elements have now been
rebuilt over \(\mathbf Q\).  Exact adjacent-pair reductions prove that the
rebuilt candidate is a Gröbner basis, and an exact unit computation proves
that it is saturated by \(a(a+1)H\).  Equality with the saturated
maximal-minor ideal is not yet proved: the next finite task is to reconstruct
fraction-free quotient identities for input containment and for
\((a(a+1)H)^{12}G\subseteq I_{15}(M_\Sigma)\).  Canonical ordered division
has now been audited at 613 support-stable good primes: its 11,701 quotient
coefficients have prohibitive height, with only 30 surviving an independent
18,116-bit balanced reconstruction test.  For input containment, the next
bounded construction is therefore a fixed syzygy-normalized
Macaulay/RREF lift with smaller coefficient height, followed by integer
identity checking.  Comparing the alternate supports comes after that
certificate.

The primary branch is broader.  The Ore localization already gives exact
Darboux coordinates and reduces the problem to a rank-one fiber Weyl algebra
over the central parameter, followed by a Hamiltonian connection.  Therefore:

1. use the simple marked-root incidence presentation, or its normalized
   factorization presentation, as the classical object over a separate
   seed/symbol parameter base;
2. quantize the rank-one fiber relatively over that base;
3. solve the Hamiltonian connection before imposing global polynomiality;
4. compute valuations on the root-at-infinity chart and test pole
   cancellation modulo homogeneous corrections and Hamiltonian gauge;
5. use support saturation to eliminate obstruction classes supported only
   on the boundary; and
6. use the normalization conductor to glue the localized solutions.

This is organized by the
[`gauge -> corrections -> defects` complex](extended-geometry/UNIFIED_DEFORMATION_COMPLEX.md).
The relative obstruction module and its Fitting/Kuranishi loci retain the
horizontal classical-symbol directions which a fixed-symbol matrix omits.
Candidate components may be discovered across several good primes, but they
must be reconstructed and verified over \(\mathbb Q\).

This replacement has now been carried out for the complete normalized
degree-five two-parameter family.  The relative complex, coherent
strong-cocycle module, Fitting and Kuranishi loci, and root-at-infinity
valuation filtration are assembled in
[`QUANTUM_RESIDUE_OBSTRUCTION.md`, Section 12](extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md#12-relative-family-package).
Four exact period sections give a parameter-uniform order-five obstruction
off a zero-dimensional scheme.  Its only interior reduced support is a
vertical rational point, already obstructed, and a cubic closed point with a
genuine \(27\)-dimensional fifth-order lift scheme.  A constant six-column
order-seven pivot now globalizes the \(X^{18}\) cocycle across that affine
space; its value is a unit in the cubic residue field, so the complete
component is obstructed.  The next relative run should therefore use a
different classical-symbol family.

The direct two-PBW-correction search remains a useful independent control.
It should use cross-prime component reconstruction and exact Ore--Gröbner
generation certificates, not a larger coefficient alphabet.  Its
one-monomial predecessor is already closed within its declared bounds.

Even complete elimination of a displayed symbol does not settle `DC_2`;
settling `DC_2` requires either a genuine non-surjective `A_2` endomorphism or
a proof that every `A_2` endomorphism is an automorphism.

Accordingly, higher-support continuation on the current degree-five symbol
is secondary.  Useful map extensions should branch across classical symbols,
keep the incidence family explicit, and make the correction complex
functorial in its parameters.

See
[`extended-geometry/RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md`](extended-geometry/RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md)
and
[`extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md`](extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md).

## 4. Plane JC(2): boundary local-cohomology obstruction

The plane boundary compiler now supplies valuation, conductor, dicritical,
and intersection data.  The remaining Case-1 obstruction should be packaged
as a section of a coherent boundary module that vanishes away from the
finitely many surviving boundary strata.

The construction target is:

1. define a boundary residue/conductor matching map;
2. let \(M\) be its cokernel, or the corresponding reflexive quotient;
3. realize the surviving Case-1 residue as a class
   \(\rho\in H_Z^0(M)\); and
4. prove
   \[
   H_Z^0(M)=0
   \]
   by depth, normality, intersection theory, or presentation saturation.

This would place the plane-JC residue problem in the same
support-saturation architecture as the degree-forty-two synchronization
defect and the cubic cotangent defect.

There is a prerequisite: the archived Case-1 certificate omits lower Newton
bands that contribute to the alternate residue.  Those bands must first be
recovered, or a truncation lemma must show that the desired residue class is
independent of them.  Local cohomology does not remove the need to define
the class and its ambient module.

Useful map extensions should produce a functorial compactified boundary
complex whose conductor/residue cokernel is coherent and whose associated
primes can be read from the compiled boundary graph.

See
[`plane-jc/FRONTIER_CLOSING_ATTACKS.md`](plane-jc/FRONTIER_CLOSING_ATTACKS.md),
[`plane-jc/LOG_BOUNDARY_COMPILER.md`](plane-jc/LOG_BOUNDARY_COMPILER.md),
and
[`plane-jc/PLANE_BOUNDARY_EXCLUSION.md`](plane-jc/PLANE_BOUNDARY_EXCLUSION.md).

## Suggested priority

1. Prove the cubic presentation saturation \(N:I^\infty=N\).
2. Build the restricted filtered deformation complex for the rank-two
   quantization obstruction.
3. Recover the missing Case-1 bands and construct the plane boundary module.

The first is the smallest exact algebraic certificate already exposed by
the repository; the second has a strong invariant reformulation; and the
third is the broadest conceptual unification but still has a prerequisite
data gap.
