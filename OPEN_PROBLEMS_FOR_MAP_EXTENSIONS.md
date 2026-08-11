# Open problems for extending the map constructions

This is a working list of theorem targets suggested by the repository's
current map constructions.  It is not a status ledger: verified claims and
their audit levels remain recorded in
[`MATH_STATUS.json`](MATH_STATUS.json) and [`STATUS.md`](STATUS.md).

The problems below are useful filters for proposed extensions.  A new family
is most valuable when it supplies one of the missing modules, obstruction
classes, support reductions, or saturation certificates described here,
rather than only adding another bounded coefficient calculation.
The cubic, plane-JC, and restricted-Weyl targets share the checked
[`G0`--`G4` support-saturation ledger](verified/SUPPORT_SATURATION_PATHS.json);
their stage labels below are machine-validated and remain open until every
upstream gate for the same module and ideal has passed.

<!-- status-consumer: SST1 12c5cb15e8b6de26 -->

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

The extracted
[`S_1` boundary theorem](verified/SUPPORT_SATURATION_PRINCIPLE.md#cubic-keller-normalization)
packages the intended structural route. Once \(C/T=0\), the collision set
has relative height two on the pure ramification support, so it is enough to
prove that \(\Omega_{B/A}\) is \(S_1\), for example by a certified
Cohen--Macaulay/perfect presentation. Normality of \(B\) and generic
conormal generation do not imply this module condition.

Operationally, Proposition 1.17 supplies `C0`; the active targets are `C1`
(the support-hull/relative-height certificate) and `C2` (the `S_1` or
associated-prime certificate).  `C3` is then a formal consequence of
`SST1`, not an additional search.

The first quartic nongauge layer is now completely stratified for all six
singular squarefree symbols.  On deterministic quotient complements of
dimensions $2,4,4,6,6,8$, strict Rees certificates for the cotangent
presentation and the cokernel of $B\to\Omega^{\oplus3}$ prove flat base
change of the cotangent module, its annihilator, and the intrinsic support
module.  Consequently every geometric parameter fiber passes `C2` and
fails `C1` by the same square-zero multiplicity-six `Ext^2` support-hull
module.  The intrinsic Kähler different also has
`dim_k J/nJ=6` on every geometric fiber, so it is not locally principal
anywhere on these quartic families.  This does not settle the global gates:
a Keller-compatible
normalization must exclude these models, or higher formal terms must change
the support defect.  Repeating quartic cotangent-saturation or specialization
searches inside the same complements is therefore no longer an active route.

Proposition 1.15a of the cubic frontend supplies the precise
different/conductor bridge under a restrictive but natural hypothesis.  If
the Kähler different is Cartier, $J=dB$, then

\[
 T^{[2]}/T\simeq (0:_{H^2_{\mathfrak n}(B)}d).
\]

Thus `C1` vanishes exactly when the local threefold normalization is
Cohen--Macaulay; normality/`S_2` and a codimension-one Cartier different do
not suffice.  The geometric programme must now prove this CM condition from
minimal-boundary geometry, or lift the certified six-generator different
through every compatible higher correction.  For the nodal symbol this
persistence is now exact to all formal orders.  Recursive gauge elimination
gives `h_nod+f(y,z)*eta`; a universal weight-one coefficient and monic graph
specialization prove that `J/nJ` remains six-dimensional for every such
formal tail.  The same multi-coefficient graph argument now closes all five
remaining singular-squarefree symbols.  Thus the entire higher-order
filtered Nakayama queue is closed.  Formal rigidity adds the smooth row, so
every compatible formal tail with any squarefree cubic symbol has a
six-generated non-Cartier intrinsic different.  A Cartier-different theorem
from boundary-minimal geometry would therefore force the leading symbol
into the double-line, triple-line, or zero rows, where generic étaleness is
the next gate.  Otherwise the active problem remains global normality and
Keller-open compatibility for the formal models.

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

<!-- status-consumer: KDCD3 c6e56b39bbb498d8 -->

<!-- status-consumer: KDSQ6 cd423f625f1f3cd2 -->

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

The [restricted-Weyl specialization of the support-saturation theorem](verified/SUPPORT_SATURATION_PRINCIPLE.md#restricted-weyl-deformation)
separates the inputs: localized solvability and conductor-compatible
primitives first put the class in boundary local cohomology; \(S_1\) plus
positive relative boundary height on \(E_k\), or on its Rees module, then
kills it. A Fitting component contained in the boundary fails the height
hypothesis and remains an exceptional component to analyze.

This is the `W0`--`W5` path.  The currently certified all-pole obstruction
stops at `W1` because it remains nonzero on the localized easy chart.
Support saturation applies only to a future locally soluble incidence family
that also passes Rees strictness `W2`, relative height `W3`, and `S_1` `W4`.

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

## 4. Plane JC(2): log-conductor and localized-Chern obstruction

The degree-independent depth theorem is now proved in
[`UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md`](plane-jc/UNIVERSAL_COMPLETE_CHAIN_BOUNDARY_SATURATION.md).
For any resolved morphism of smooth SNC surface pairs attached to a plane
Keller map, the logarithmic differential is an injective square map of
rank-two bundles.  Its cokernel is perfect of codimension one, hence pure
Cohen--Macaulay and `S1`.  Therefore

\[
H_Z^0(\mathcal T_f^{\log})=0
\]

for every zero-dimensional boundary collision locus, uniformly in the
coordinate degrees and the length of the complete chain.

The same note gives exact constant nodal-tree countermodels to the stronger
claim that the conductor/gauge matching cokernel is automatically `S1` or
has positive collision height.  Simultaneous normalization, a tree dual
graph, and free conductor and jet modules do not control the associated
primes of an arbitrary matching matrix.

The
[`log-conductor degree-shift theorem`](plane-jc/LOG_CONDUCTOR_DEGREE_SHIFT.md)
corrects the proposed comparison.  A normalization mismatch of a
torsion-free curve module embeds canonically in `H_Z^1`, and admits no
nonzero `O`-linear lift into `H_Z^0`.  Subtracting the complete logarithmic
different makes the scalar normalized determinant a global unit with zero
conductor mismatch.  Nevertheless two genuine local polynomial Jacobian
matrices can have the same determinant and generic branch Smith profiles
while their nodal cokernels are the glued node ring and its split
normalization.  Their `Fitt_1` ideals distinguish them.

The primary conceptual target is therefore the **full nodal logarithmic
matrix theorem**:

1. transport the complete `2 x 2` logarithmic differential through every
   Newton/blowup chart, including boundary discrepancies;
2. compute its nodal `Fitt_1`, normalization defect, and localized second
   Chern contribution;
3. identify a terminal Laurent class in that defect and prove it is nonzero;
4. independently prove that Keller geometry makes the same class descend.

The last two assertions would contradict each other.  Terminal type-I
determinant nonvanishing alone does not supply item 3: in the `(75,125)`
terminal block, `[P,Q]_(X,y)=X^4` normalizes to `1`, so its scalar conductor
mismatch is zero.  Existing degree calculations are regression fixtures for
the full-matrix/Chern compiler.

The first full local profiles are now known.  In the `(75,125)` terminal
packet, all three interior target-node attachments have transverse order
`-2`, require two source blowups, and become log-etale at both nodes of the
resulting attachment chain.  The endpoint `s=0` is log-etale as well.  At the
source endpoint `s=infinity` over a smooth target-boundary point, Keller
boundary support forces the unit-`Fitt_1` cokernel to be exactly `R/(w^3)`
on a smooth reduced support.  Hence none of the five terminal marked slots
supplies item 3.

The carrier-local continuation is also complete for both certified F2 rows.
The two squarefree simple-root spectators require one blowup each; the
double-row fivefold point requires four carrier-centered blowups and one
fan-alignment ray; and every principal arm requires six further common-fan
components.  All resulting exponent determinants are `1`, `3`, or `5`, so
these nodes are log-etale in characteristic zero.  Together with the
terminal audit this strengthens the source-boundary lower bounds to `27/48`
components, but still produces no normalization defect.

The upstream carrier-extraction audit changes the problem qualitatively.  Its
carrier-zero ladder is unimodular, but at the node between the first
exceptional and the strict line at infinity the exact model is
`R/(W^3 U^18)`.  Its canonical branchwise splitting has quotient
`R/(W^3,U^18)` of length `54`.  This is a forced nonzero degree-one matching
class, while `H_Z^0` remains zero.  The main target is no longer to find a
local class: it is to prove a global Keller/localized-`c_2` identity that
cannot accommodate this class after the affine purity row and other possible
cancellation centers are compiled.  The outgoing terminal-tail theorem now
maps the remaining source rays unimodularly to the existing target fan and
proves that they contribute no defect or correction.

The affine-purity frontier is now typed exactly.  Purity forces a new
boundary divisor with `e>1` over an affine nonproperness curve, because none
of the certified components has that generic profile.  The source-component
floors therefore rise from `27/48` to `28/49`.  The possible degrees remain
`6..9375` and `12..9375`, every ramified row obeys
`2<=e<=d-1` and `f<=floor((d-1)/e)`, and the target curve has a polynomial
parametrization of degree at most `124`.  A coarse row `(e,f)=(2,1)` with
affine degree `d-2` survives at every such `d`, so generic purity cannot
select the row or raise the degree floor.  The live input is the actual
nonproperness curve and the complete factorization of its pullback.

The target-curve theorem makes that input finite and collision-theoretic.
Every F2 nonproperness component has normalization degrees `(3k,5k)` for
`1<=k<=24`, projective degree `5k<=120`, unique infinity point `[0:1:0]`,
and top implicit form `P^(5k)`.  It must be singular, so its two normalization
divided differences generate a proper collision/critical ideal.  The next
computation is therefore a 24-chart collision stratification followed by
implicitization and pullback factorization, not an unrestricted
degree-124 curve search.

The first chart is already explicit.  For `k=1`, target and parameter
normalization give `p=t^3+a*t` and
`q=t^5+b*t^4+c*t^2+d*t`; every collision is indexed by a root of
`u^4+b*u^3+a*u^2+(2*a*b-c)u-(a^2+d)`.  A nonempty open subchart has four
ordinary affine nodes, and with the `(2,5)` infinity cusp its conductor
ledger exhausts the genus-six quintic budget.  The unresolved step is now
the source pullback over these marked quartic roots and the degenerate
quartic strata, followed by the remaining `k=2..24` charts.

The generic `k=1` target equation is now explicit: one twelve-support
quintic `F(P,Q)`.  Its normalization pullback satisfies
`(F_P,F_Q)=C(t)*(q',-p')`, with one degree-eight polynomial `C` cutting out
the eight preimages of the four nodes.  The remaining operation is no longer
target implicitization but factorization of this fixed equation after the
F2 Laurent substitution.

<!-- status-consumer: PF2K1I1 a7582c1e36140840 -->

The fixed-coordinate pullback interface is explicit too.  Every inverse
target normalization is encoded by five scalars, gives a denominator-free
quintic, and computes the carrier residue and all eight relevant infinity
jets.  Because a Keller map is étale on the affine plane, the source
pullback is automatically reduced and its affine singular/conductor scheme
is exactly the base change of the four target nodes.  Each affine node has
normalization-defect length one and conductor-divisor degree two; only four
finite fiber counts remain.  At geometric
degree `d`, each count is at most `d-1`, so the total affine nodal conductor
is at most `4(d-1)` (`20/44` at the two degree floors).  The unresolved
factorization is therefore confined to the source boundary, where the
purity divisor and its logarithmic point corrections live.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The conductor calculation now includes every degenerate `k=1` stratum.
Coprimality of degrees `3,5` forces birationality, the fixed `(2,5)` infinity
cusp leaves affine delta invariant `4`, and the same monic degree-eight
resultant is the exact conductor divisor even when nodes coalesce into
cusps, tacnodes, or higher multiple fibers.  Under a degree-`d` Keller
pullback, the affine normalization defect is at most `4(d-1)` and the affine
conductor-divisor degree at most `8(d-1)`.  Degenerations can change only the
distribution and finite fiber counts, not the total target conductor.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The conductor length is not itself a logarithmic point correction.  For an
fs tame Kummer toroidal packet over a resolved target node, both the split
model and the collided cyclic model `z^e=x*y` (including its full `(-2)`
resolution chain) are log-étale.  Their logarithmic cokernel and localized
`ch_2` vanish.  The remaining boundary calculation must therefore test the
actual completed source map against the toroidal chart criterion; only a
failure of that test can survive.  Even then, the general SNC matrix gives a
second gate: full-rank exponent data force zero cokernel, while rank one can
have singular determinant support only if two explicit logarithmic-unit
first jets vanish simultaneously.

<!-- status-consumer: PF2K1TN1 521fb57f7e6abc1f -->

This separation is now complete on the affine plane.  Any embedded
resolution of any reduced target curve pulls back through a Keller map to a
strict étale log resolution, so nodes, cusps, tacnodes, and higher multiple
points all have zero relative logarithmic cokernel.  Their ordinary
conductor lengths remain finite-fiber and boundary-escape data, but neither
those lengths nor carrier-parameter corank can be inserted into the source
`Fitt_1/ch_2` module.  Only compactification-boundary attachments remain;
under the minimal transverse unibranch SNC hypotheses their exact point
length is locally `q_p*m_C` and totals `m_C*f` over a complete
residue-degree-`f` fiber.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

The puncture is now located as well.  Every chart meets the extracted target
divisor `(5,2)` with contact `k` and residue
`lambda=p_lead^5/(-q_lead)^3`.  A formal direct comparison with the terminal
cover would force `lambda=125/729,e=3`, but that source neighborhood is
already a resolved morphism: every further exceptional divisor over it maps
to one point and cannot dominate the affine target curve.  Thus the terminal
passport does not determine the affine index, and the former conditional
`29/50` count is unavailable.  The unconditional floors stay `28/49`.

The blowup-conservation theorem identifies the correct birational quantity.
Raw matching length changes under further boundary blowups, but

\[
\mathcal Q(D)=\frac12\sum m_i^2C_i^2+\sum_{i\text{--}j}m_im_j
=\frac12D^2
\]

is conserved.  For the extraction-root cycle `D=3E+18L` it equals `27`.
The kernel-line theorem identifies the remaining cyclic twist exactly.  If
`K` is the kernel of the restricted logarithmic differential, then the
cokernel line is `K tensor O_D(D)` and
`ch_2(coker)=deg_D(K)+D^2/2`.  On a packet contracted to a target point, the
kernel line is tautological for a Gauss map `gamma:D->P^1`, so
`deg_D(K)=-deg(gamma^*O(1))`.  Thus the F2 root contribution is
`27-e_root<=27` for a nonnegative Gauss degree.  The tangential-coordinate
theorem now computes it: `f^*z=W^3U^18*unit` makes the fixed target covector
`dz` generate the kernel modulo the full thickened Cartier ideal, so
`e_root=0` and the cyclic root contribution is exactly `27`.  Noncyclic
nodal terms, the outgoing and affine rows, and the global cancellation
identity are now the central gaps.

<!-- status-consumer: UCBS1 824720a8f727bdf8 -->

<!-- status-consumer: LCDS1 5b4d92acd50d6c41 -->

<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->

<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->

<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->

<!-- status-consumer: LCBBC1 b3eb4679f781c55f -->

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

<!-- status-consumer: LKGD1 8a357250b5005186 -->

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

An alternative route may retain the original matching cokernel `M`.  The
[plane specialization of the general theorem](verified/SUPPORT_SATURATION_PRINCIPLE.md#plane-jc-conductor-residue)
then shows exactly what must follow construction: \(M\) must be \(S_1\) and
\(Z\) must have positive relative height on its support.  `UCBS1` proves that
neither property follows merely from a nodal tree, simultaneous
normalization, or free conductor and jet modules.  They would have to be new
Keller-specific theorems about the actual matching map.  If the entire
nonzero matching cokernel is already supported on \(Z\), support saturation
cannot annihilate it; one must prove componentwise surjectivity or enlarge
the family.

For this alternative route, the order remains `P0` (construct the
truncation-independent matching module), `P1` (construct and localize the
residue), `P2` (positive relative height), `P3`
(`S_1`/associated-prime exclusion), then conditional conclusion `P4`.  No
finite-support wording may skip `P1`'s divisorial-different check.  This
coefficient-base route is distinct from the surface normalization defect:
purity of the latter does not prove `P2`--`P3` for the former.

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
3. Compile the full nodal logarithmic matrices and compare their localized
   second Chern length with the global Chern identity, using the completed
   `(72,108)` and `(75,125)` packages as exact regression fixtures.

The first is the smallest exact algebraic certificate already exposed by
the repository; the second has a strong invariant reformulation; and the
third is the broadest conceptual unification.  The determinant and generic
Smith data are now proved insufficient; the nodal `Fitt_1` profile is the
first invariant capable of seeing the required extension defect.

<!-- status-consumer: C1FBC1 0f14ef01fff25097 -->
