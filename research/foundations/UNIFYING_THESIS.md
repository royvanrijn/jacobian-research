# The unifying thesis: marked roots, boundary cancellation, and algebraization

This note gives the conceptual spine shared by the weighted, cancellation,
boundary, arithmetic, stable-equivalence, and quantization parts of the
repository.  It is a synthesis of proved results and an explicit statement
of the remaining classification problem.  It is not a new assertion that
every noninvertible Keller map has already been classified.

The thesis is:

> **Boundary-cancelled marked-root algebraization.**  Begin with the finite
> normalization of a polynomial incidence cover whose generic inverse is
> obtained by selecting a root of an equation `Psi(T;y)=0`.  Choose a
> weighted, birational, or reciprocal source chart in which that root is a
> coordinate.  Along the selected ramified boundary, the marked-root
> derivative `partial_T Psi` supplies the ramification zero and the chart
> and reconstruction functions supply the complementary pole.  After the
> target chart is included, their valuations balance in the determinant
> ledger, so the algebraized polynomial map has constant nonzero Jacobian.
> Global polynomial inversion nevertheless fails because affine source
> space is only the regular-reconstruction open in the normalized root
> cover.  Across its omitted boundary, the marked root may remain finite
> while another reconstruction coordinate has negative valuation.

Thus “opposite valuation” means opposite **net** valuation in the full
source--core--target ledger; it need not be a literal two-factor identity in
every presentation.  In the reciprocal cancellation chart it is the direct
product `D^r D^{-r}=1`.  In the weighted chart the source and target vertical
weights distribute the same cancellation among three determinants.

Once the controlled suspension square and its selected critical
normalization have been intrinsically identified, the two verified branches
have different logarithmic types:

| branch | selected critical normalization | logarithmic type | status |
|---|---|---|---|
| weighted tangent | `A^1` | one puncture in its `P^1` completion | verified |
| cancellation/quadratic | `G_m` | two punctures in its `P^1` completion | verified |
| possible new mechanisms | a rational curve with at least three punctures, several selected boundary components, a nonrational normalization, or nontrivial conductor gluing | higher or nonclassical logarithmic boundary | open |

This is a taxonomy of the selected ramified locus after suspension
extraction.  It is not the claim that a presentation-dependent plane slice
is itself a stable invariant; the canonical three-dimensional boundary
strata and their decorations are the objects used for stable comparison.

The package directly explains the constant Jacobian and nonproperness, and
it organizes the additional family-specific results on complete collision
fibres, monodromy, stable inequivalence, and moduli.  It also explains why
formal or local inverses can exist while a global polynomial inverse fails,
and why conductor and nonreduced collision data can distinguish cases that
have the same reduced discriminant.  Section 6 records exactly which of
these consequences are direct and which require extra theorems.

The recurring dichotomy is therefore

\[
\boxed{
\text{formal or local solvability}
\quad\text{versus}\quad
\text{global polynomial algebraization}.
}
\tag{T}
\]

The first side is controlled by a Jacobian identity, an etale inverse, a
formal automorphism, or a completed correction complex.  The second asks
whether those local solutions assemble into one finite polynomial object
on the required affine space.

## 1. The finite marked-root package

Let `k` be a field of characteristic zero, let `Y` be a normal affine
`k`-variety, and let

\[
F:U=\mathbb A_k^d\longrightarrow Y
\]

be dominant and quasi-finite.  Put

\[
\bar X_F=\operatorname{Norm}_Y k(U),\qquad
\pi_F:\bar X_F\longrightarrow Y.
\tag{1.1}
\]

Zariski's Main Theorem gives a unique open immersion

\[
j_F:U\hookrightarrow\bar X_F,\qquad
F=\pi_F\circ j_F.
\tag{1.2}
\]

The **normalization boundary** is

\[
\partial_F=(\bar X_F\setminus j_F(U))_{\mathrm{red}}.
\tag{1.3}
\]

This finite cover, its distinguished affine open, and its boundary are
intrinsic to `F`.  A primitive inverse polynomial is useful for computation,
but it is not itself invariant under changing the primitive element.

### Marked-root presentation

A marked-root presentation adds a separable irreducible equation

\[
\Psi(T;y)=0\qquad\text{over }k(Y)
\tag{1.4}
\]

which generates `k(U)/k(Y)`, together with rational reconstruction functions

\[
x_1=X_1(T,y),\ldots,x_d=X_d(T,y)\in k(\bar X_F).
\tag{1.5}
\]

The **regular-reconstruction open** is the largest open on which the
functions in (1.5), the target functions, and both inverse identities are
regular.  In the constructions of this repository it is exactly `j_F(U)`,
and the functions `x_i` identify it with affine source space.

There are consequently three distinct algebraic objects:

1. the finite normalization `\bar X_F`, which retains all generic root
   sheets;
2. the regular-reconstruction open `U`, on which a marked root gives finite
   polynomial source coordinates;
3. the boundary `\partial_F`, where at least one reconstruction coordinate
   has a pole.

Calling `F` a “finite cover” would conflate the first two objects.  The map
`\pi_F` is finite; the nonproper Keller map `F` is its restriction to the
distinguished open.

## 2. Controlled-boundary suspensions

The exact reusable statement is the
[boundary-cancelled incidence lemma](cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md#1-boundary-cancelled-incidence-lemma).
It packages seven obligations: finite root cover, source-chart valuation,
controlled divisor, polynomiality, determinant ledger, boundary
normalization, and collision fibre.

Its determinant square is

```text
 U = A^3  -------- F -------->  Y = A^3
    |                              |
    | alpha                        | beta
    v                              v
 Z -------------- Phi ----------> T
```

If `J_Phi=uD^r`, its one common determinant condition is

\[
\boxed{
\operatorname{div}(\det D\alpha)
+r\operatorname{div}(D\circ\alpha)
=F^*\operatorname{div}(\det D\beta).
}
\tag{2.1}
\]

Pullback through a rational chart is interpreted in the common function
field, so negative coefficients record source-chart poles. The lemma proves
the constant determinant only after the polynomiality obligation has also
been discharged, and it proves completeness of a collision only after every
specialized root has been shown to lie in the reconstruction open.

The weighted, cancellation, and quadratic-gauge families are now derived as
three explicit dictionaries under that lemma. Their determinant arguments
are identical at this level; their genuine differences are respectively a
distributed polynomial ledger, a reciprocal `D^{-r}` chart with a finite
jet gate, and a reciprocal `D^{-1}` marked-line chart with a
coefficient-weight gate.

## 3. Boundary--reconstruction and nonproperness theorem

The next result makes “escape to infinity” intrinsic.  It combines the
affine-versus-boundary selector in
[stable normalization functoriality](verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
with the finite Zariski--Main factorization.

### Theorem 3.1

Let `F:A^d->Y` be as in Section 1, and assume that the normalization in
(1.1) is finite.  Write `x_1,...,x_d` for the affine source coordinates.

1. For every prime divisor `E` of `\bar X_F`,
   \[
   E\subset\partial_F
   \quad\Longleftrightarrow\quad
   \min_i v_E(x_i)<0.
   \tag{3.1}
   \]
   Thus a divisorial branch is omitted from affine source space exactly when
   at least one source reconstruction coordinate has a pole there.  Moreover,
   every point of `\partial_F` lies on such a divisor; in particular, the
   irreducible components of the boundary have pure codimension one.
2. Define `S_F` to be the set of points `y in Y` for which no Zariski
   neighborhood `V` of `y` makes
   `F^{-1}(V)->V` proper.  Then
   \[
   \boxed{S_F=\pi_F(\partial_F).}
   \tag{3.2}
   \]
3. If `k=C`, every generic complex point of a component of
   `\pi_F(\partial_F)` is the limit of targets `F(p_n)` with the source
   points `p_n` escaping every compact subset of `A^d(C)`.

#### Proof

A divisorial valuation `v` of `k(x_1,...,x_d)` has a center on `A^d` exactly
when it is nonnegative on `k[x_1,...,x_d]`.  This is equivalent to
`v(x_i)>=0` for every coordinate generator: sums and products of
nonnegative-valuation elements still have nonnegative valuation.

If the generic point of `E` lies in `j_F(U)`, every polynomial source
coordinate is regular there, so all `v_E(x_i)` are nonnegative.  Conversely,
if all are nonnegative, the valuation ring defines a center on `A^d`.  Its
image under the open immersion `j_F` is a center on `\bar X_F`.  Since
`\bar X_F` is separated, a valuation has at most one center; that center must
be the generic point of `E`.  Hence the generic point lies in `j_F(U)`.
This proves (3.1).

Now take any point `e` of the boundary.  If every `x_i` were regular at `e`,
then on some neighborhood `W` they would define a morphism
`g:W->A^d`.  The two maps \(j_F\circ g\) and
\(\operatorname{id}_W\) to `\bar X_F` agree on the dense open
\(W\cap j_F(U)\).  Separatedness makes them equal on `W`, which
would put `e` in `j_F(U)`, a contradiction.  Hence some `x_i` is not in the
normal local ring `O_(\bar X_F,e)`.  A normal noetherian domain is the
intersection of its height-one valuation rings, so that `x_i` has negative
order on a prime divisor through `e`.  By (3.1), that divisor lies in the
boundary.  This proves the purity assertion.

The finite image `\pi_F(\partial_F)` is closed.  Over its complement,
`\pi_F^{-1}(Y\setminus\pi_F(\partial_F))` lies entirely in `j_F(U)`, so `F`
is the restriction of a finite map with no missing points and is therefore
finite, hence proper.

Now let `y` lie in `\pi_F(\partial_F)` and choose
`e in \partial_F` over `y`.  If `F` were proper over a neighborhood `V` of
`y`, then the `V`-morphism

\[
j_F:F^{-1}(V)\longrightarrow\pi_F^{-1}(V)
\]

would be proper because its source is proper over `V` and its target is
separated over `V`.  Its image would be closed.  But it is also a dense open
of the integral scheme `\pi_F^{-1}(V)` and omits `e`, a contradiction.
This proves (3.2).

Over `C`, choose a curve in `\bar X_F` through a general point of `E` which
is not contained in the boundary.  Its puncture lies in `U`.  By (3.1), some
`x_i` has negative order along the curve, so its absolute value is unbounded,
while finiteness of `\pi_F` makes the target converge to `\pi_F(E)`.  This
proves the last assertion.  QED

### Corollary 3.2

If `F:A^d->A^d` is Keller, then it has no affine critical locus.  If it is
nonproper, its failure of invertibility is nevertheless visible in the
nonempty finite-normalization boundary.  In particular,

\[
\text{etale on affine source}
\quad\not\Longrightarrow\quad
\text{proper or globally polynomially invertible}.
\tag{3.3}
\]

The obstruction is a pole of global reconstruction, not a zero of `det DF`
at a finite source point.

## 4. The three verified realizations

The framework above is a theorem for the weighted, cancellation, and
root-engineered quadratic-gauge families. Exhaustiveness beyond the stated
families is not asserted.

### Theorem 4.1 -- realization theorem

Every admissible weighted marked-root map, every polynomial cancellation
map, and every root-engineered quadratic-gauge map is a
boundary-cancelled incidence in the sense of Lemma 1.2.

The proof is the three-family dictionary in the
[common lemma](cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md#3-the-three-established-families-are-instances).
The determinant part is now one argument. The remaining proof obligations
are family-specific:

| family | polynomiality input | normalized critical curve | collision input |
|---|---|---|---|
| weighted | admissible weighted support | `A^1` | a squarefree specialized tangent pencil with regular reconstruction |
| cancellation | `L_(m,r)(h)=0` | `G_m` | the squarefree integral pencil at `(1,0,0)` |
| quadratic gauge | coefficient-weight identity for `G_P` | `G_m` | the prescribed squarefree seed `G` at `(1,0,0)` |

The exact weighted incidence and Fitting calculation is in the
[tangent-map core](verified/TANGENT_MAP_CORE.md). The finite jet gate and
reconstruction are in the
[cancellation construction](cancellation/CONSTRUCTION.md). The marked-line
polynomiality identity, dual normalization, and programmable fibre are in the
[quadratic-gauge theorem](cancellation/ROOT_ENGINEERED_QUADRATIC_GAUGE.md).
QED

## 5. The foundational cubic in all three families

Put `u=1+xy` and

\[
F_0(x,y,z)=
\left(
\begin{aligned}
&u^3z+y^2u(4+3xy),\\
&y+3xu^2z+3xy^2(4+3xy),\\
&2x-3x^2y-x^3z
\end{aligned}
\right).
\tag{5.1}
\]

This is the smallest example in which the polynomial weighted type and both
reciprocal constructions give the same polynomial Keller map up to
polynomial left--right equivalence.

### Weighted chart

Take

\[
H(W)=W^2(1-W).
\tag{5.2}
\]

In the invariant-plane coordinates of the foundational weighted map,

\[
\frac Q4=H'(W)+\gamma,\qquad
\frac P4=W(H'(W)+\gamma)-H(W).
\tag{5.3}
\]

Thus (5.3) is the tangent incidence (4.1).  On the generic discriminant
boundary,

\[
\partial_W\Psi=-\gamma,\qquad x=\frac C\gamma,
\]

so

\[
v_E(\gamma)=1,\qquad v_E(x)=-1.
\tag{5.4}
\]

The bottom core ramifies, the full threefold map remains Keller, and the
selected affine source coordinate escapes.

### Cancellation chart

For `(m,r)=(1,1)`, the unique cancellation parameter and normalized jet are

\[
q=3,\qquad h(A)=3+9A.
\tag{5.5}
\]

Set

\[
A=1+xy,\qquad
B=A^2z+y^2(3+9A),\qquad
P=AB,\qquad Q=y+xB,
\tag{5.6}
\]

\[
R_C=C\int_0^{x/A}\{1-t(Q-Pt)\}\,dt.
\tag{5.7}
\]

The exact identity

\[
F_0
=\operatorname{diag}(1/3,1,2/C)
\circ(P,Q,R_C)\circ(x,y,3z)
\tag{5.8}
\]

uses only polynomial source and target automorphisms.  On the generic
critical boundary,

\[
D=A^{-1},\qquad x=sD^{-1},
\]

and therefore

\[
v_E(D)=1,\qquad v_E(x)=-1.
\tag{5.9}
\]

### Quadratic gauge

The root-engineered seed `G(S)=S^3+S` makes the normalized all-degree
quadratic-gauge formula exactly `F_0`. Its marked-line core has derivative
`D`, its reciprocal source chart contributes `D^{-1}`, and its
coefficient-weight polynomiality sum is empty. Thus it is a third incidence
presentation of the same cubic, but it belongs to the same reciprocal
suspension type as the cancellation chart.

In the weighted and cancellation charts the unique ramified boundary prime over the cubic
discriminant has ramification index two and carries a simple source pole.
The left--right equivalence (5.8), together with functoriality of the finite
normalization package, identifies these as the same intrinsic boundary
valuation.  What changes is the ledger used to cancel it: polynomial
vertical weights in (5.3), a reciprocal rational chart in (5.7).

The complete branch-collapse calculation is
[Theorem 4.1 of the minimal-boundary note](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md#4-cubic-collapse-inside-the-two-branches).

This cubic overlap is exhaustive for the cancellation and root-engineered
quadratic-gauge families.  In every degree at least four, the
ramified-stratum Fitting supports have affine ranks one and two,
respectively; independently, the target-boundary contact has nilpotency
index `mr(m+1)` and `2`.  Thus even the higher `m=1` stationary-point
ladder does not meet the root-engineered quadratic-gauge family stably.  See
the
[quadratic-gauge/cancellation stable-intersection theorem](verified/QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md).

## 6. What follows, and what only echoes the thesis

The following table records the logical status of the repository's main
themes.  “Direct” means that the feature is part of the marked-root,
determinant, or normalization package.  “Derived” means that an additional
family-specific theorem is required.  “Analogy” means that the same
local/global dichotomy reappears, not that the result is a formal corollary
of Sections 1--4.

| Feature | Status | Precise relation to the thesis |
|---|---|---|
| inverse pencils and reconstruction denominators | direct | The pencil presents the generic finite extension; negative valuations of the reconstruction functions cut out the omitted affine boundary. |
| discriminant normalization and Hessian Fitting divisor | direct for the tangent core | Repeated marked roots give the ramification divisor.  Its restriction normalizes the discriminant, and its differential module gives the full Fitting divisor. |
| nonproperness without finite ramification | direct | Theorem 3.1 identifies nonproperness with the finite image of `partial_F`, while the determinant ledger makes `F` etale on affine source. |
| boundary valuations, differents, and conductors | direct after the relevant strata are intrinsically selected | They refine the finite normalization and its scheme-theoretic gluing; stable functoriality is supplied by `S1`. |
| finite root covers and complete fibers | derived | The normalization supplies the sheets.  Explicit seed and target choices are additionally needed to keep every sheet in the reconstruction open and to prescribe its field of definition. |
| nonreduced collision intersections | derived | Family-specific elimination or completed-local calculations show that cancelled orders survive as contact multiplicity and nilpotence. |
| stable distinction and stable moduli | derived | One must prove that the boundary or marked-Fitting decorations are intrinsic and faithful before applying stable normalization functoriality. |
| Ritt-component intersections | derived, not formal | Decomposition strata of the seed polynomial induce special loci in the inverse family.  Their nonreduced overlaps and conductor equalizers require separate Ritt and intersection theorems. |
| arithmetic variation of fiber splitting | derived | Arithmetic monodromy acts on the same finite root cover, but symmetric monodromy, Chebotarev, and local engineering are additional inputs. |
| formal source-orbit triviality versus stable polynomial equivalence | structural analogy | Every Artin deformation is formally source-trivial; one global polynomial automorphism family may still fail to exist or have unbounded degree. |
| formal symplectic solvability versus polynomial quantization | structural analogy | Formal-etale correction complexes may be solvable while finite polynomial or Laurent corrections remain obstructed by poles, filtration, or residue classes. |

The distinction in the last two rows is essential.  Marked-root geometry
motivates the algebraization question, but it does not by itself prove a
quantum obstruction.  The exact sources are
[formal orbit triviality](extended-geometry/FORMAL_ORBIT_TRIVIALITY.md),
[smallest-false-rank `DC_2` program](extended-geometry/FIXED_RANK_DIXMIER_REDUCTION.md),
and the
[quantum residue obstruction](extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md).

## 7. The classification conjecture

The proved realization theorem runs from explicit formulas to intrinsic
geometry:

\[
\text{weighted or cancellation formula}
\Longrightarrow
\text{controlled ledger and marked-root boundary}.
\tag{7.1}
\]

The open direction is to reconstruct the formula from the intrinsic
Zariski--Main package.  The sharpened classification strategy is therefore

\[
\boxed{
\text{operational intrinsic boundary invariants}
\Longrightarrow
\text{verified gateway and controlled suspension square}
\Longrightarrow
\text{formula and polynomial algebraization}.
}
\tag{7.2}
\]

In particular, one should classify the selected normalization, its
punctures, boundary valuations, multiple components, and conductor gluing
before classifying coordinate formulas.

### Minimal-boundary gateway and suspension conjecture

Let

\[
F:\mathbb A^3\longrightarrow\mathbb A^3
\]

be a nonproper Keller map whose canonical finite-normalization package passes
the following proposed intrinsic gateway tests:

1. a canonical, left--right-functorial rule selects one geometrically
   integral rational critical-boundary normalization;
2. saturated affine links with no undeclared zero or pole support;
3. an intrinsically oriented positive link is boundary-monotone;
4. at most two places at infinity on the critical normalization; and
5. a divisor-minimal ledger has no unrecorded graph valuation, and deleting
   a recorded boundary destroys the saturated link or ledger.

Then `F` is polynomially left--right equivalent to either a weighted tangent
suspension or a cancellation suspension.

This is the existing
[minimal-boundary gateway and classification conjecture](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md),
restated here to show its role in the thesis.  The five items are proposed
gateway conditions, not a currently available decision procedure.  In
particular, the repository still has to define or extract canonically the
selector, orientation, graph audit, and deletion test without first
recognizing one of the desired suspension formulas.  The statement is
therefore a classification conjecture relative to a conjectural intrinsic
gateway, not a classification theorem with independently checkable
hypotheses and not a conjecture about every nonproper Keller map.

The required implication chain is:

\[
\boxed{
\text{operational invariants of the canonical finite normalization}
\Longrightarrow
\text{verified gateway conditions}
\Longrightarrow
\text{coordinate-preserving controlled suspension square}
\Longrightarrow
\text{suspension formula}.
}
\tag{7.3}
\]

The first arrow includes canonical marking selection and an exhaustive
graph-divisor audit.  The second includes noncontraction and primitive
conormal extraction.  The last is conditional on positive-chart
straightening in the positive branch.  None of these steps may use
coordinates obtained from an already recognized weighted or cancellation
presentation.

Once the square is extracted, the number of places on the normalized
critical curve separates the expected branches:

\[
\mathbb A^1
\rightsquigarrow
\text{weighted tangent core},\qquad
\mathbb G_m
\rightsquigarrow
\text{reciprocal cancellation core}.
\tag{7.4}
\]

The remaining algebraization problem is then branch-specific:

- on the positive branch, straighten the polynomial vertical chart to the
  weighted form;
- on the reciprocal branch, extract the primitive two-place marking and
  prove the finite cancellation condition;
- at closed collision points, exclude the normalization/Fitting defects not
  visible in the height-one ledger.

In geometric degree three, the two extracted branches already collapse to
the same foundational map by Section 5.  The active task is therefore not to
compare two cubic answers.  After making the overarching intrinsic gateway
operational, it is enough to complete either cubic frontend: suspension
extraction, or flat binary-cubic normalization followed by coefficient-gauge
straightening.

## 8. A practical claim discipline

Future results can be placed in this framework by answering five questions.

1. **Finite object:** What is the intrinsic finite normalization or
   incidence cover?
2. **Affine chart:** Which valuations distinguish the polynomial source from
   its omitted boundary?
3. **Local cancellation:** What divisor appears in the core Jacobian, and
   where is its reciprocal factor recorded in the determinant ledger?
4. **Algebraization:** What finite divisibility, saturation, degree, or
   conductor condition makes the rational construction polynomial?
5. **Logical level:** Is the conclusion a direct consequence of the package,
   a family-specific derived theorem, a structural analogy, or an open
   classification assertion?

This checklist keeps the unifying thesis falsifiable.  A new family would
support it only after exhibiting both its determinant ledger and its
regular-reconstruction boundary.  A construction with formal cancellation
but no polynomial algebraization would illustrate the dichotomy without
producing a Keller map.

## 9. Canonical sources

The theorem authorities used in this synthesis are:

- [Marked-root Keller maps](MARKED_ROOT_KELLER_MAPS.md) for the common inverse
  framework;
- [Tangent-map core](verified/TANGENT_MAP_CORE.md) for the weighted plane
  incidence and determinant square;
- [Stable normalization functoriality](verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
  for the intrinsic Zariski--Main package and valuation selector;
- [Cancellation construction](cancellation/CONSTRUCTION.md) for the
  reciprocal chart, polynomiality criterion, and reconstruction;
- [Boundary geometry](cancellation/BOUNDARY_GEOMETRY.md) for complete
  divisorial exhaustion and nonreduced contact;
- [Controlled-boundary suspensions](cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md)
  for the determinant ledger;
- [Minimal-boundary classification](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md)
  for the scoped conjecture and cubic branch collapse; and
- [Exact image and nonproperness](verified/IMAGE_AND_NONPROPERNESS.md) for the
  foundational analytic escape model.

This note organizes those results; it does not replace their proofs or their
status entries in `MATH_STATUS.json`.
