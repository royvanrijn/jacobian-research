# Support-saturation and boundary extension

Let \(S\) be a commutative Noetherian ring, let
\(\mathfrak a\subset S\) be an ideal,
and let \(M\) be a finite \(S\)-module.  Write
\[
 H^0_{\mathfrak a}(M)
 =\{m\in M:\mathfrak a^n m=0\text{ for some }n\}.
\]
This module is the universal obstruction to extending a vanishing statement
from \(\operatorname{Spec}S\setminus V(\mathfrak a)\) across
\(V(\mathfrak a)\).

> **Status and scope.**  The module theorem, its presentation form, the
> linear conductor corollary, and the Rees corollary below are proved under
> their displayed hypotheses.  The exact associated-prime theorem and the
> presentation-saturation equivalence are kernel-checked in
> [`formal/support-saturation`](../formal/support-saturation/README.md).  The
> geometric \(S_1\)-to-associated-prime bridge and the conductor argument are
> ordinary commutative-algebra proofs.  The derived three-row criterion is a
> separate proved result imported from the boundary-obstruction note.  All
> three programme specializations are conditional: this note identifies their
> missing hypotheses but does not assert that those hypotheses are known.

<!-- status-consumer: SST1F 838e558b5fcb9d81 -->

## The theorem

> **Theorem 1 (support-saturation principle).**  The following conditions
> are equivalent:
>
> 1. \(H^0_{\mathfrak a}(M)=0\);
> 2. \(0:_M\mathfrak a^\infty=0\);
> 3. no associated prime of \(M\) contains \(\mathfrak a\);
> 4. \(\operatorname{grade}(\mathfrak a,M)\ge1\);
> 5. \(\mathfrak a\) contains an \(M\)-regular element.
>
> If \(F\) is finite free and \(M=F/N\), these are also equivalent to
> \[
>  N:_F\mathfrak a^\infty=N.                              \tag{1}
> \]

**Proof.**  Noetherianity makes the ascending chain
\[
 0:_M\mathfrak a\subseteq0:_M\mathfrak a^2\subseteq\cdots
\]
stationary, and its union is \(H^0_{\mathfrak a}(M)\).  This proves the
equivalence of the first two conditions.  The associated primes of this
submodule are exactly
\[
 \operatorname{Ass}(M)\cap V(\mathfrak a).
\]
Thus it vanishes exactly when no associated prime of \(M\) contains
\(\mathfrak a\).  Prime avoidance for the finite set
\(\operatorname{Ass}(M)\) then identifies this with the existence of an
\(M\)-regular element in \(\mathfrak a\), equivalently positive grade.
Finally, an element \(f+N\in F/N\) is killed by a power of
\(\mathfrak a\) exactly when
\(\mathfrak a^n f\subseteq N\) for some \(n\), which gives
\[
 H^0_{\mathfrak a}(F/N)
 \simeq (N:_F\mathfrak a^\infty)/N
\]
and proves (1).  \(\square\)

## Defects and completion

> **Corollary 2 (extension across a support).**  If
> \(d\in M\) vanishes after localization at every prime outside
> \(V(\mathfrak a)\), then
> \[
>  d\in H^0_{\mathfrak a}(M).
> \]
> Consequently any condition in Theorem 1 forces \(d=0\).

**Proof.**  The cyclic module \(Sd\) has support in
\(V(\mathfrak a)\).  Hence
\(\mathfrak a\subseteq\sqrt{\operatorname{Ann}(d)}\), and finite
generation of \(\mathfrak a\) gives
\(\mathfrak a^n d=0\) for some \(n\).  \(\square\)

Let \(\widehat S\) be a flat adic completion of \(S\).  Since the
annihilator chain stabilizes, flat base change gives
\[
 H^0_{\mathfrak a}(M)\otimes_S\widehat S
 \simeq
 H^0_{\mathfrak a\widehat S}(M\otimes_S\widehat S).         \tag{2}
\]
In particular, the precompletion saturation equality (1) is sufficient
for the completed module to have no \(\mathfrak a\)-torsion.

## Structural shortcuts

Theorem 1 is deliberately weaker than flatness.  Useful sufficient
conditions are:

- \(M\) is torsion-free over an integral base and
  \(\mathfrak a\) contains a nonzero base element;
- \(M\) has no embedded associated primes and every minimal component
  avoids \(V(\mathfrak a)\);
- \(M\) is a finite maximal Cohen--Macaulay module over a regular base:
  Auslander--Buchsbaum makes it locally free, hence torsion-free.

The last statement requires module-finiteness and full depth over the
regular base.  Cohen--Macaulayness of an algebra without these hypotheses
does not by itself imply support saturation.

## Two applications

1. In the degree-forty-two Hessian residual problem, the synchronization
   defect vanishes away from an explicitly determined closed support.
   The remaining global target is the saturation of the residual ideal by
   the ideal defining that support.
2. In the cubic-normalization frontend, the closed-point cotangent torsion
   is already presented as
   \[
   (N:_F\mathfrak a^\infty)/N.
   \]

They are the algebra and module versions of the same theorem: a defect
known to vanish off a closed set extends across it exactly when the
ambient module has no torsion supported there.

## The geometric depth theorem

Here \(M\) satisfies Serre's condition \(S_1\) when

\[
 \operatorname{depth}_{S_{\mathfrak p}}M_{\mathfrak p}
 \geq
 \min\!\left(1,\dim\operatorname{Supp}_{S_{\mathfrak p}}
                         M_{\mathfrak p}\right)
 \quad(\mathfrak p\in\operatorname{Supp}M).
\]

For a finite \(S\)-module \(M\), define the relative height of
\(\mathfrak a\) on \(M\) by

\[
 \operatorname{ht}_M(\mathfrak a)
 =\inf_{\substack{\mathfrak p\in\operatorname{Supp}M\\
                   \mathfrak a\subseteq\mathfrak p}}
   \dim M_{\mathfrak p},
\tag{3}
\]

with value \(+\infty\) if the indexing set is empty.  The inequality
\(\operatorname{ht}_M(\mathfrak a)\geq1\) says exactly that no irreducible
component of \(\operatorname{Supp}M\) is contained in
\(V(\mathfrak a)\).  This is the precise meaning of a sufficiently deep
support ideal in the module theorem.  It is a relative-height condition,
not a statement about a high power of the ideal.

> **Theorem 3 (the \(S_1\) boundary support-saturation theorem).**  Let
> \(S\) be Noetherian, let \(M\) be a finite \(S\)-module satisfying
> Serre's condition \(S_1\), and let
> \(\mathfrak a\subset S\) satisfy
> \(\operatorname{ht}_M(\mathfrak a)\geq1\).  Then
> \[
>  H^0_{\mathfrak a}(M)=0.
> \tag{4}
> \]
> Consequently every section of \(M\) which vanishes off
> \(V(\mathfrak a)\) vanishes globally.  If \(M=F/N\) with \(F\) finite
> free, then
> \[
>  N:_F\mathfrak a^\infty=N.
> \tag{5}
> \]

**Proof.**  The \(S_1\) condition says that \(M\) has no embedded
associated primes:

\[
 \operatorname{Ass}(M)=\operatorname{Min}(\operatorname{Supp}M).
\]

Relative height at least one says that no prime on the right contains
\(\mathfrak a\).  Theorem 1 now gives (4) and (5), while Corollary 2 gives
the assertion about sections.  \(\square\)

The theorem has two useful certificate forms.

1. If \(M\) is Cohen--Macaulay and every component of its support has
   positive relative height along \(\mathfrak a\), then (4) holds.  Over a
   regular ring, a finite free resolution proved perfect of the expected
   codimension at every support prime is one standard way to certify this
   Cohen--Macaulay hypothesis.  An lci or matrix-factorization description
   is useful only when it actually gives this depth statement for \(M\).
2. For a finite presentation of \(M\),
   \(\operatorname{Supp}M=V(\operatorname{Fitt}_0(M))\).  Hence the
   relative-height condition is checked by proving that no minimal prime of
   \(\operatorname{Fitt}_0(M)\) contains \(\mathfrak a\).  The \(S_1\)
   condition is still separate; generic local freeness and the Fitting
   support alone do not exclude embedded primes.

The separation is necessary.  Let \(R\) be a regular domain and let
\(I\subset R\) be any proper ideal, of arbitrarily large height.  Then

\[
 M=R\oplus R/I
\tag{6}
\]

is generically free, but
\(H_I^0(M)\supseteq R/I\ne0\).  This example survives with the trivial
normalization and trivial conductor.  It disproves any version in which
finite normalization, generic local freeness, or ambient codimension
replaces \(S_1\).  Replacing \(I\) by a high power also changes nothing,
since \(H^0_{I^n}(M)=H^0_I(M)\).

## Finite extensions and conductor-assisted extension

The normalization and conductor do a different job from depth: they show
that locally constructed corrections agree away from a residual collision
locus.  The \(S_1\) theorem then prevents that residual locus from carrying
a defect.

Let \(A\subset B\) be a finite injective ring extension and let

\[
 \mathfrak c=\operatorname{Ann}_A(B/A)
\]

be the conductor.  It is also an ideal of \(B\), and there is an exact
conductor sequence

\[
 0\longrightarrow A\longrightarrow
 B\oplus A/\mathfrak c\longrightarrow
 B/\mathfrak c\longrightarrow0.
\tag{7}
\]

The maps are \(a\mapsto(a,\bar a)\) and
\((b,\bar a)\mapsto\bar b-\bar a\).  Thus (7) simply states
\(A=B\times_{B/\mathfrak c}A/\mathfrak c\).  A finite normalization is the
application of interest, but normality is not needed for this gluing lemma.

> **Theorem 4 (finite-extension boundary saturation).**  Let
> \(d:C^1\to C^2\) be a map of finite projective \(A\)-modules and put
> \(E=\operatorname{coker}d\).  Let \(I\subset A\), let
> \(o\in C^2\), and write \(e=[o]\in E\).  Suppose that for every
> \(\mathfrak p\notin V(I)\) there are correction primitives
> \[
>  x_{\nu,\mathfrak p}\in
>       C^1_{\mathfrak p}\otimes_A B,
>  \qquad
>  x_{c,\mathfrak p}\in
>       C^1_{\mathfrak p}\otimes_A A/\mathfrak c
> \]
> satisfying \(d(x_{\nu,\mathfrak p})=o\),
> \(d(x_{c,\mathfrak p})=o\), and having the same image over
> \(B/\mathfrak c\).  Here both equations mean equality after the indicated
> scalar extension.  If
>
> 1. \(E\) satisfies \(S_1\); and
> 2. \(\operatorname{ht}_E(I)\geq1\),
>
> then \(e=0\).  In fact \(H_I^0(E)=0\), so the conclusion holds for every
> defect which has such conductor-compatible local corrections.

**Proof.**  Tensor (7) with \(C^1_{\mathfrak p}\).  Projectivity preserves
exactness, so the compatible pair
\((x_{\nu,\mathfrak p},x_{c,\mathfrak p})\) glues to an element
\(x_{\mathfrak p}\in C^1_{\mathfrak p}\).  The element
\(d(x_{\mathfrak p})-o\) maps to zero in both terms on the right of (7)
after tensoring with \(C^2_{\mathfrak p}\).  Exactness and projectivity of
\(C^2\) give \(d(x_{\mathfrak p})=o\).  Thus
\(e_{\mathfrak p}=0\) outside \(V(I)\), and Corollary 2 puts \(e\) in
\(H_I^0(E)\).  Theorem 3 makes this module zero.  \(\square\)

The equality of the two primitives over the conductor quotient is the
linear form of the required gauge/correction compatibility.  For a perfect
complex it should be replaced by compatible null-homotopies.  The
[perfect-complex conductor theorem](../extended-geometry/BOUNDARY_OBSTRUCTION_THEORY.md#4-theorem-b-perfect-conductor-descent)
gives the corresponding homotopy-Cartesian statement.  For a three-term
complex, genuinely derived vanishing also has lower cohomology rows: the
same note proves that it is enough to have

\[
 \operatorname{grade}(I,H^q(C^\bullet))\geq3-q
 \qquad(q=0,1,2).
\tag{8}
\]

A geometric sufficient package for (8) is

\[
 H^q(C^\bullet)\text{ satisfies }S_{3-q},
 \qquad
 \operatorname{ht}_{H^q(C^\bullet)}(I)\geq3-q.
\tag{9}
\]

Thus the module theorem is the \(q=2\) edge of the derived theorem, not a
replacement for the lower gauge and correction compatibilities.

## Rees modules and filtered defects

Let \(F_\bullet A\) and \(F_\bullet E\) be good filtrations and write

\[
 \mathcal R(A)=\bigoplus_nF_nA\,t^n,
 \qquad
 \mathcal R(E)=\bigoplus_nF_nE\,t^n.
\]

> **Corollary 5 (Rees support saturation).**  Suppose
> \(\mathcal R(E)\) is finite over the Noetherian ring
> \(\mathcal R(A)\), satisfies \(S_1\), and a homogeneous boundary ideal
> \(J\subset\mathcal R(A)\) has
> \(\operatorname{ht}_{\mathcal R(E)}(J)\geq1\).  Then
> \[
>  H_J^0(\mathcal R(E))=0.
> \tag{10}
> \]
> In particular, a homogeneous filtered defect which vanishes off
> \(V(J)\) is zero before specialization at \(t=1\).

This is Theorem 3 applied to the Rees module.  The hypotheses deliberately
refer to the actual Rees presentation.  Associated-graded exactness alone
does not prove Rees strictness, and \(S_1\) of one special fiber does not
automatically imply \(S_1\) of the Rees module.

## Operational gate contract

Every distinguished-defect application below must pass the same five gates.
An application which proves the entire module $H_I^0(E)$ vanishes may omit
`G1`, because there is no separately chosen defect whose localized vanishing
must first be established; the cubic path is of this type.  Its remaining
gates and every declared exception are recorded in the machine-readable
programme ledger
[`SUPPORT_SATURATION_PATHS.json`](SUPPORT_SATURATION_PATHS.json).

| Gate | Required certificate | Failure meaning |
|---|---|---|
| `G0` finite obstruction module | A finite module \(E\) over a Noetherian coefficient ring, its boundary ideal \(I\), and a finite presentation when a colon computation is used | There is no coherent global receptacle to saturate |
| `G1` localized vanishing | The defect is zero off \(V(I)\); normalization-chart corrections must agree over the conductor quotient | The class is a genuine open-locus obstruction, not boundary torsion |
| `G2` positive relative height | No minimal component of \(\operatorname{Supp}E\) is contained in \(V(I)\) | A whole obstruction component lies on the boundary and cannot be removed by saturation |
| `G3` no embedded boundary support | \(E\) is \(S_1\), or directly no associated prime of \(E\) contains \(I\) | An embedded \(I\)-supported class survives |
| `G4` support saturation | Apply Theorem 1/3 to obtain \(H_I^0(E)=0\), equivalently \(N:I^\infty=N\) | This gate is conditional until `G0`--`G3` refer to the same module and ideal |

Conductor descent is one way to discharge `G1`; it does not discharge
`G2` or `G3`.  For filtered problems, Rees strictness is an additional
adapter before `G2`--`G4`: those gates must be checked on the actual Rees
obstruction module.  Finite-tail effectivity remains separate even after
`G4`, because vanishing produces a global correction but does not bound its
filtration degree.

## Three conditional specializations

The following are the exact programme interfaces.  They are corollaries
under displayed hypotheses, not claims that those hypotheses have already
been proved in the three applications.

### Cubic Keller normalization

Use the notation of the
[cubic normalization frontend](../cancellation/CUBIC_NORMALIZATION_FRONTEND.md):

\[
 Q=\Omega_{B/A},\qquad
 T=B/\operatorname{Ann}_B(Q),\qquad
 C=T^{[2]},\qquad L=C/T,
\]

and let \(Z=V(I)\) be the finite collision locus.  Assume \(L=0\), so the
pure two-dimensional ramification support is \(S_2\).  If \(Q\), which has
rank one and full support over \(T\), satisfies \(S_1\), then

\[
 \operatorname{ht}_Q(I)\geq2,
 \qquad
 H_I^0(Q)=0.
\tag{11}
\]

For a presentation \(F_1\to F_0\to Q\to0\), with image \(N\), this is

\[
 N:_{F_0}I^\infty=N.
\tag{12}
\]

Proposition 1.16 of the frontend then identifies the primitive-conormal
cokernel with this local cohomology module, so cotangent cyclicity and
point-flatness follow.  The new theorem isolates the missing input exactly:
after \(L=0\), prove \(S_1\) for \(\Omega_{B/A}\), or certify it through a
Cohen--Macaulay/perfect presentation.  Normality of \(B\) and generic
codimension-one conormal generation do not by themselves prove this module
depth statement.

The exact quartic-model frontier now sharpens that warning.  On deterministic
complements of dimensions $2,4,4,6,6,8$ to the formal-gauge image for the
six singular squarefree cubic symbols, strict weighted-Rees presentations
commute the cotangent module, its annihilator, and the intrinsic support
module with every geometric parameter specialization.  Every fiber has
cotangent saturation but retains the same square-zero multiplicity-six
`Ext^2_A(T,A)` support-hull obstruction.  Thus all quartic fibers pass `C2`
and fail `C1`.  The global theorem still requires Keller geometry to exclude
the support defect and does not promote either programme gate beyond the
quartic layer.

The different/conductor interface is now exact under a Cartier hypothesis.
If locally `Ann_B(Q)=dB`, Proposition 1.15a of the frontend gives

\[
 L\simeq(0:_{H^2_{\mathfrak n}(B)}d),
\]

so `C1` vanishes precisely when the normal threefold $B$ is
Cohen--Macaulay at the collision.  This records the additional depth input
which codimension-one conductor control does not supply.  The same quartic
certificate computes `J/nJ` as a locally free rank-six parameter module, so
these six model families are nowhere Cartier-different.  Proposition 1.15c
isolates the all-orders mechanism: flatness of `Q` and the
annihilator-cokernel packet commutes `J` with base change, while local
freeness of `J/nJ` preserves the six-generator count by Nakayama.  See
[Exact computation 1.8f](../cancellation/CUBIC_NORMALIZATION_FRONTEND.md#exact-computation-18f----singular-squarefree-nongauge-double-saturation).

For the nodal symbol this persistence now holds to all formal orders.  The
cyclic quotient `ker(C)/im(G_nod)=Q[y,z](-3)` gives the recursive normal
form `h_nod+f(y,z)*eta`.  For the universal coefficient `u`, the relative
Nakayama module is `Q[u]^6`; weight-one strict Rees packets identify the
associated graded cotangent and annihilator-control modules with their
central polynomial extensions.  Every graph equation `u-f` has monic
initial form, so filtered regularity commutes the intrinsic annihilator
with arbitrary polynomial and formal graph specialization.  This is the
all-orders nodal instance of the relative persistence theorem.

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

The multi-coefficient form of the same argument closes all six singular
squarefree symbols.  Their exact gauge-cokernel generator counts are
`1,2,2,3,3,4`.  In every universal normal-coefficient family the relative
Nakayama module is free of rank six, while weight-one Rees packets for the
cotangent and annihilator-control cokernel are strict with central initial
modules.  Successive monic graph equations therefore commute the intrinsic
annihilator with every formal tail.  This proves all-orders non-Cartier
persistence without requiring the relative `Ext^2` presentation itself to
be constant.  The smooth row is the zero-coefficient case by formal
rigidity, and its central different is also six-generated.  Hence the
all-orders non-Cartier statement covers every squarefree cubic symbol.

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

<!-- status-consumer: KDCD3 c6e56b39bbb498d8 -->

### Plane-JC conductor residue

Let \(M=\operatorname{coker}\Phi\) and
\(\rho\in M\) be the matching module and residue proposed in
[the plane boundary atlas](../extended-geometry/PLANE_WILD_BOUNDARY_ATLAS.md#9-characteristic-zero-jc_2-boundary-module).
Once the conductor/contact-loss ledger proves that the finite presentation
is independent of the omitted tail and \(\rho\) vanishes off
\(Z=V(I_Z)\), Theorem 4 gives

\[
 M\text{ is }S_1,\qquad
 \operatorname{ht}_M(I_Z)\geq1
 \quad\Longrightarrow\quad
 \rho=0.
\tag{13}
\]

A restrictive sufficient certificate is that \(M\) be Cohen--Macaulay of
positive relative dimension along \(Z\), with no minimal Fitting component
contained in \(Z\).  This also states a necessary warning: if \(M\) itself
is a nonzero finite module supported entirely on \(Z\), then

\[
 H_Z^0(M)=M,
\]

and no support-saturation argument can kill \(\rho\).  One must instead
prove that the matching map is surjective on that component or enlarge the
family so that \(Z\) has positive relative height.

### Restricted Weyl deformation

At filtered order \(m\), let

\[
 E_m=\operatorname{coker}(d_{1,m})
\]

be the coherent obstruction module on a classical-symbol or incidence base,
and let \(\mathfrak b\) be its boundary ideal.  Suppose the localized
correction equation is solvable, the normalization-chart primitives satisfy
the conductor compatibility of Theorem 4, and the chosen finite filtered
presentation is Rees-strict.  Then either of the following packages kills
the global order-\(m\) obstruction:

\[
 E_m\text{ is }S_1,\qquad
 \operatorname{ht}_{E_m}(\mathfrak b)\geq1,
\tag{14}
\]

or the corresponding \(S_1\) and height conditions (10) on the Rees
obstruction module.  Minimal primes are visible on the appropriate Fitting
support.  An exceptional Fitting component contained in the boundary fails
the height hypothesis and must be analyzed rather than removed by
saturation.

For a derived correction torsor, (14) controls only the \(H^2\) edge; the
lower \(H^1\) and \(H^0\) modules require (8) or explicit compatible
primitives.  The already certified all-pole degree-five Weyl obstructions do
not satisfy the localized-solvability hypothesis, so this corollary does not
alter their status.  It gives the precise depth/Fitting target for future
incidence families that are locally solvable.

## What the theorem separates

The reusable chain is therefore

\[
\begin{array}{c}
\text{normalization/open solution}\\
\text{plus conductor-compatible corrections}
\end{array}
\Longrightarrow
 e\in H_I^0(E)
\xRightarrow[\operatorname{ht}_E(I)\ge1]{E\text{ is }S_1}
 e=0.
\tag{15}
\]

Finite normalization and the conductor control where a defect can live.
Fitting ideals and generic local freeness identify its components.  The
\(S_1\)/Cohen--Macaulay hypothesis excludes embedded boundary support, and
positive relative height excludes a boundary component.  Rees strictness
ensures that these statements concern the actual filtered deformation
module.  None of these four jobs implies the others.

## Verification and formalization ledger

| Layer | Mathematical status | Verification | What remains |
|---|---|---|---|
| Associated-prime support saturation | Proved | Lean: [`Core.lean`](../formal/support-saturation/SupportSaturation/Core.lean) | Nothing at this abstraction level |
| Quotient/presentation saturation | Proved | Lean: [`Presentation.lean`](../formal/support-saturation/SupportSaturation/Presentation.lean) | Nothing at this abstraction level |
| \(S_1\) geometric wrapper | Proved | Written reduction using \(S_1\iff\) no embedded associated primes | Add a reusable module-depth/Serre API to Mathlib |
| Linear conductor gluing | Proved | Written proof from the explicit fibre-product sequence (7) | Formalize the conductor square and base-changed cokernel class in Lean |
| Three-row derived criterion | Proved separately | [`BOUNDARY_OBSTRUCTION_THEORY.md`](../extended-geometry/BOUNDARY_OBSTRUCTION_THEORY.md) and its checker | Formalize local cohomology and the hypercohomology spectral sequence |
| Rees corollary | Proved conditionally on the actual Rees hypotheses | Direct application of Theorem 3 | Each programme must certify finiteness, \(S_1\), height, and strictness for its Rees module |
| Cubic, plane, Weyl interfaces | Conditional corollaries | Hypotheses displayed in (11)--(14); machine-readable routes in [`SUPPORT_SATURATION_PATHS.json`](SUPPORT_SATURATION_PATHS.json), checked by [`verify_support_saturation_paths.py`](../scripts/verify_support_saturation_paths.py) | Prove the listed depth, height, local-solvability, conductor-compatibility, and strictness inputs in each programme |

Replay the cross-program routing and claim-boundary check with

```bash
make verify-support-saturation-paths
```

The next formal milestones should therefore be kept separate:

1. formalize a finite-module \(S_1\) predicate and prove its equivalence with
   absence of embedded associated primes;
2. formalize the conductor fibre product for a finite injective extension and
   Theorem 4 for a two-term complex;
3. formalize the local-cohomology spectral sequence only after the two-term
   statement is stable;
4. build application-specific certificates for \(\Omega_{B/A}\), the plane
   matching cokernel, and the restricted-Weyl Rees obstruction modules.

The final item is where new mathematics remains.  The abstract theorem cannot
replace those certificates: it turns each programme's boundary problem into a
finite list of depth and support obligations.

## Standard references

- The equivalence between \(S_1\) and absence of embedded associated primes is
  [Stacks Project, Lemma 10.157.2](https://stacks.math.columbia.edu/tag/031O).
- Associated-prime conventions and finiteness are collected in
  [Stacks Project, Section 10.63](https://stacks.math.columbia.edu/tag/00L9).
- Flat base change for local cohomology is
  [Stacks Project, Lemma 51.5.7](https://stacks.math.columbia.edu/tag/0EF5).
- The Fitting-support identity
  \(V(\operatorname{Fitt}_0(M))=\operatorname{Supp}M\) is
  [Stacks Project, Lemma 15.8.4](https://stacks.math.columbia.edu/tag/07ZA).
- The Lean representation of \(H_I^0(M)\) is Mathlib's
  [`Ideal.primaryComponent`](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Algebra/Module/Torsion/PrimaryComponent.html).
