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
and intersection data.  The numerical `(72,108)` Case-1 systems have since
been closed by the quotient-first determinantal argument in
[`JC2_72_108_DETERMINANTAL_CLOSURE.md`](plane-jc/JC2_72_108_DETERMINANTAL_CLOSURE.md).
The remaining conceptual target is to package an alternate-chart Case-1
residue, for this or a future boundary stratum, as a section of a coherent
boundary module that vanishes away from its finite support.

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

The characteristic-divisible wild-boundary atlas now supplies an explicit
candidate.  If
`nu:Btilde->B` normalizes the compiled boundary, use the conductor quotient
`Q_B=nu_*O_Btilde/O_B`, twist it by the determinant line with the compiled
different removed, and take the cokernel `M` of the source-, target-, and
gauge-jet matching map.  The resulting supported residue belongs to
`H_Z^0(M)=0:_M I_Z^infinity`; see
[`PLANE_WILD_BOUNDARY_ATLAS.md`, Section 9](extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md#9-characteristic-zero-jc_2-boundary-module).

The same atlas now supplies a necessary preliminary gate.  For its canonical
rows in characteristics `3,5,7`, deleting the fierce boundary leaves a whole
tame ramification divisor with different lengths `1,3,5`.  A node-supported
local-cohomology class cannot cancel that divisorial defect.  Any future
boundary compiler must therefore verify vanishing of the relative different
away from the proposed omitted boundary before forming the finite-support
module below.

The balanced replacement `PQT -> P^(N-1)QT` now passes this preliminary gate:
its normalized complement is étale by the integral-slope Newton polygon and
purity.  Its natural birational affine-plane chart has Jacobian
`-u^(2N-4)`, however, so only `N=2` is Keller in that chart.  The divisor-class
localization now closes the alternative-chart loophole as well:
`Cl(C-E)=Z/(N-1)`, namely `Z/2,Z/4,Z/6` in the rows `N=3,5,7` and
`Z/5` in both mixed-residue `N=6` controls.  More precisely, the named
boundary class `[L1]` has exact order `N-1`.
Thus these complements are not affine planes despite having the same
base-field point counts as `A^2`.  A future presentation must cancel both the
residual different exposed by the original gluing and the boundary-class
torsion exposed by balanced gluing.  The monomial-gluing dichotomy proves that
no intermediate coefficient `P^a*Q*T`, `1<=a<=N-1`, escapes: the different
vanishes only at the balanced endpoint.  Any next candidate must therefore
change more than the `P`-valuation of a single `Q*T` gluing monomial.  The
extension to all `a>=0` closes the omitted powers as well: `a=0` has a free
unit and `Z/(N+1)` class torsion, while `a>=N` has a wild index-`N` branch.
Moreover an arbitrary `C(P)` with a factor away from `P=0` creates another
height-one different component.  Hence, under the same target-support
hypothesis, the one-variable coefficient search is now exhausted rather than
merely bounded.  With a general target coefficient `C(P,Q)`, every factor
away from `P*Q=0` is likewise excluded, reducing the search exactly to
`C=cP^aQ^b`.  This quadrant is now closed.  For balanced `a=N-1`, the
prime-to-characteristic part of `b+1` changes the compactly supported Euler
characteristic, pure `p`-powers with `N>2` preserve `[L1]` by push--pull, and
the residual `p=N=2` tower has
`Cl(U_(2,c))=Cl(D(xu))=(Z/2)^2` at every Frobenius height; the source-fill
valuation matrix is unimodular.  For the other `a`-ranges, the nonconstant
unit or the generic
`P=0` tame/wild ramification row is unchanged.  Thus the original
characteristic-two cubic is the unique affine-plane Keller row under this
one-omitted-fierce-boundary and target-support hypothesis.  A next candidate
must change the boundary support or the gluing architecture itself.
For balanced multi-retained gluing, the exact identity
`A*H_T-A'*H=P^(N-1)*Q*(A-T*A')` now forces
`A=a0+T*B(T^p)` under the same target-support condition.  Hence only cover
degrees congruent to one modulo `p` survive this preliminary gate.  The
normalization has the exact root-count formula
`#(C_A-E_A)(F_q)=q^2+(n_q(A)-1)q` and, after the retained roots split,
`chi_c(C_A-E_A)=deg(A)`.  Every nonlinear retained polynomial is therefore
geometrically non-affine-plane.  The linear row is already excluded by
`Cl=Z/(N-1)` for `N>2`, so the characteristic-two cubic is the unique row in
the balanced squarefree-retained architecture.  The former support-only
queue through degree `15`,
`(p,d)=(3,7),(3,10),(3,13),(5,11),(7,15)`, is now closed rather than awaiting
source certificates.
The degree-seven `F_3` row is now completely scanned over its six admissible
retained polynomials.  All normalizations are smooth and pass the different
gate; four fail the `F_3` point count, leaving only
`A=T^4+T+1,T^4+2T+1` at that field.  Both continue to match `A^2` over
`F_9` but have `810` open points over `F_27` instead of `729`; the extra
`3*27` comes from the cubic retained roots.  Hence the complete first odd row
has no geometric survivor.  The next odd-characteristic work must change the
boundary presentation, use a nonsquarefree/colliding retained boundary with
its own conductor ledger, or leave balanced single-polynomial gluing.

This calculation is now a general search theorem rather than a family-only
trick.  For any normal candidate open `U` with a certified dense class-
trivial core `W` (for example an affine UFD), free based unit lattice
`M=Gamma(W,O_W)^*/k^*`, and complement
whose codimension-one primes are `D_i`, the full valuation map

\[
 M\longrightarrow\bigoplus_i\mathbb ZD_i
\]

has kernel `Gamma(U,O)^*/k^*` and cokernel `Cl(U)`.  It also computes the
exact order of a named reflexive class from an augmented matrix.  For the
special core `G_m^2`, an affine-plane candidate must have exactly two fill
primes and a unimodular valuation matrix.  The hypotheses, proof, and
executable gate are in the
[boundary-lattice prefilter](plane-jc/BOUNDARY_LATTICE_PREFILTER.md#dual-torus-core-localization).
This applies to other finite-cover, conductor-fill, Cox, and affine-
modification searches whenever such a class-trivial core can be proved; it is not
valid for a selected boundary subset or an unproved unit-lattice presentation.

There is also a reusable multiple-fibre shortcut.  If prime fibres have
multiplicities `m_i`, their reduced sum is principal, and the geometric
generic-fibre unit group modulo base constants has rank one, their vertical
classes form `(direct_sum Z/m_i)/<diagonal>`.  Its order is
`prod(m_i)/lcm(m_i)`, so it is nonzero exactly when two multiplicities share
a prime.  This rejects affine-space reconstruction without a full core
class-group computation; if the generic class group and all other vertical
prime classes vanish, it computes the full class group.  Multiplicities alone
are not enough.

The filter also extends to cores with nonzero finitely presented class group.
If `Cl(W)=coker(R)` and rational witnesses for the lifted core relations have
boundary-valuation matrix `A`, then the complete reconstruction class group is
the cokernel of `[[V,A],[0,R]]`.  This is strictly stronger than computing the
boundary quotient and `Cl(W)` separately: the same `V=R=(2)` gives either
`Z/2 + Z/2` or `Z/4` according as `A=0` or `A=1`.  Thus searches may use
singular or nonfactorial cores when their class presentation and relation
lifts are exact; they need not discard them merely because `Cl(W)` is nonzero.

<!-- status-consumer: BL1 e86cdcd66993bccc -->

The abstract truncation prerequisite is now solved by
[`CONDUCTOR_JET_TRUNCATION.md`](plane-jc/CONDUCTOR_JET_TRUNCATION.md).  On a
normalization branch with conductor exponent `c` and certified expression
contact loss `lambda`, the conductor matching map and residue class depend
only on inputs modulo `t^(c+lambda)`; the bound is sharp in general.  Its
dependency-sensitive refinement tests `n_j >= c+lambda_(alpha,j)` only for
actual input-to-output paths.  A normal-valuation adapter compiles `n_j` from
a certified valuation frontier of omitted Newton exponents.  The archived
Case-1 conductor application remains uncertified because its matching map,
expression paths, and band-to-normal valuation data have not been compiled.
The independent exact continuation in
[`CASE1_FULL_BAND_CONTINUATION.md`](plane-jc/CASE1_FULL_BAND_CONTINUATION.md)
has nevertheless recovered all eleven lower bands directly.  The remaining
Case-1 work is alternate-chart transport, construction of the matching map
and residue, and the local-cohomology calculation; no input-band deficit
remains.  Local cohomology still does not manufacture the class or its
ambient module.

Useful map extensions should produce a functorial compactified boundary
complex whose conductor/residue cokernel is coherent and whose associated
primes can be read from the compiled boundary graph.  The compiler output
must now also separate three layers: the residual different divisor, the
normal-core valuation/Smith module (including lifted class relations when
`Cl(W)` is nonzero), and only then the finite-support local-cohomology residue.

See
[`plane-jc/FRONTIER_CLOSING_ATTACKS.md`](plane-jc/FRONTIER_CLOSING_ATTACKS.md),
[`plane-jc/LOG_BOUNDARY_COMPILER.md`](plane-jc/LOG_BOUNDARY_COMPILER.md),
[`plane-jc/PLANE_BOUNDARY_EXCLUSION.md`](plane-jc/PLANE_BOUNDARY_EXCLUSION.md),
and
[`extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md`](extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md).

<!-- status-consumer: PWB1 4ce9a0bf6d277321 -->
<!-- status-consumer: PWB2 2346d64d0f1eaa07 -->
<!-- status-consumer: PWB3 a6bcf405759ddd5d -->
<!-- status-consumer: PWB4 ebddf245e65b62a7 -->
<!-- status-consumer: PWB5 142b02344181fed3 -->
<!-- status-consumer: PWB6 35636805d73e0bec -->
<!-- status-consumer: CJT1 afb70f90ff10f3d7 -->

## Suggested priority

1. Prove the cubic presentation saturation \(N:I^\infty=N\).
2. Build the restricted filtered deformation complex for the rank-two
   quantization obstruction.
3. Transport the completed Case-1 bands to the alternate chart, compile the
   matching map, and construct the plane boundary module and residue class.

The first is the smallest exact algebraic certificate already exposed by
the repository; the second has a strong invariant reformulation; and the
third is the broadest conceptual unification; its former Newton-band data gap
is closed, while the map-and-module construction remains open.

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->
