# The unifying thesis: marked roots, boundary cancellation, and algebraization

This note gives the conceptual spine shared by the weighted, cancellation,
boundary, arithmetic, stable-equivalence, and quantization parts of the
repository.  It is a synthesis of proved results and an explicit statement
of the remaining classification problem.  It is not a new assertion that
every noninvertible Keller map has already been classified.

The thesis is:

> A noninvertible Keller map in the present families is built from a finite
> marked-root incidence whose plane core has an apparent ramification
> divisor.  A vertical map or rational source chart contributes the
> complementary Jacobian factor, so the suspended polynomial map is etale.
> Global invertibility still fails because affine source space is only the
> regular-reconstruction open in the finite normalized root cover.  At its
> omitted boundary, a reconstruction coordinate has negative valuation and
> source points escape to infinity.

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

The determinant mechanism is most naturally a commutative square

```text
 U = A^3  -------- F -------->  Y = A^3
    |                              |
    | alpha                        | beta
    v                              v
 Z -------------- Phi ----------> T
```

in which `alpha` and `beta` are dominant rational maps of smooth
threefolds, `Phi` is a plane map with one parameter coordinate, and

\[
\beta\circ F=\Phi\circ\alpha.
\tag{2.1}
\]

Suppose the only nonunit factor of the core Jacobian is a divisor `D=0` and

\[
\det D\Phi=uD^r,\qquad u\in k^\times.
\tag{2.2}
\]

Taking determinants in (2.1) gives

\[
(\det D\beta\circ F)\det DF
=u(D\circ\alpha)^r\det D\alpha.
\tag{2.3}
\]

The corresponding equality of principal Weil divisors is

\[
\boxed{
\operatorname{div}(\det D\alpha)
+r\operatorname{div}(D\circ\alpha)
=F^*\operatorname{div}(\det D\beta).
}
\tag{2.4}
\]

This is the **determinant ledger**.  Pullback through a rational chart is
interpreted in the common function field, so negative coefficients are
allowed.  When (2.4) holds, `det DF` is a global unit; on affine space it is
a nonzero constant.

A **controlled-boundary marked-root suspension** consists of:

1. a marked-root package (1.1)--(1.5);
2. a square (2.1) with core factorization (2.2);
3. the balanced determinant ledger (2.4);
4. a finite list of boundary valuations supporting all poles of the rational
   charts and reconstruction functions; and
5. a polynomial-algebraization proof showing that the displayed target
   coordinates pulled back to `U` are polynomials.

The ledger proves local etaleness of the final map.  The fifth condition is
logically separate.  In the cancellation family it is a finite jet
congruence; without it, the same rational determinant identity does not
produce a polynomial map of affine spaces.

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

## 4. The two verified realizations

The framework above is a theorem for the weighted and cancellation families.
Exhaustiveness beyond the stated families is not asserted.

### Theorem 4.1 -- realization theorem

Every admissible weighted marked-root map and every polynomial cancellation
map is a controlled-boundary marked-root suspension.

#### Weighted suspension

For a seed `H`, the plane tangent core is

\[
\Phi_H(W,\gamma)
=\bigl(H'(W)+c\gamma,\,
W(H'(W)+c\gamma)-H(W)\bigr)
\tag{4.1}
\]

with

\[
\det D\Phi_H=-c^2\gamma.
\tag{4.2}
\]

For

\[
\gamma=1+a_0xy+b_0x^2z,\qquad
W=(1+xy)\gamma,\qquad C=x\gamma,
\]

the source and target vertical determinants are

\[
\det D\alpha=b_0x^3\gamma^2,\qquad
\det D\beta=-cC^3.
\tag{4.3}
\]

Since `C=x gamma`, the three powers of `gamma` in the core/source product
match the three powers in the target ledger, and

\[
\det DF=b_0c.
\tag{4.4}
\]

The inverse incidence is

\[
\Psi_H(W;A,B,C)=H(W)-BCW+cAC^2=0,
\tag{4.5}
\]

and

\[
\partial_W\Psi_H=H'(W)-BC=-c\gamma,\qquad
x=-\frac{cC}{\partial_W\Psi_H}.
\tag{4.6}
\]

Thus the repeated-root equation, the core Jacobian zero, and the generic
reconstruction pole are one divisor viewed in three coordinate systems.
The exact normalization, including the additional `C=0` charts, is proved in
the [weighted marked-root theorem](verified/WEIGHTED_SEED_THEOREM.md).

#### Cancellation suspension

Put

\[
D(s,P,Q)=1-s(Q-Ps)^m
\]

and

\[
R(s,P,Q)=C\int_0^sD(t,P,Q)^r\,dt.
\tag{4.7}
\]

The plane core over `P` has

\[
\det\frac{\partial(P,Q,R)}{\partial(s,P,Q)}=CD^r.
\tag{4.8}
\]

On the rational source chart,

\[
D=A^{-1},\qquad
\det\frac{\partial(s,P,Q)}{\partial(x,y,z)}
=-A^r=-D^{-r}.
\tag{4.9}
\]

The factors cancel directly and give

\[
\det DF=-C.
\tag{4.10}
\]

The inverse equation and reconstruction are

\[
\Psi_{m,r}(T)
=C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R=0,
\tag{4.11}
\]

\[
y=Q-TP,\qquad A=D^{-1},\qquad x=TD^{-1},
\tag{4.12}
\]

with the displayed formula for `z` in the
[cancellation construction](cancellation/CONSTRUCTION.md).  Here

\[
\partial_T\Psi_{m,r}=CD^r.
\tag{4.13}
\]

The determinant identity holds for every parameter polynomial `h`, but the
last target coordinate is polynomial on `A^3` exactly when the finite
cancellation operator satisfies

\[
\mathcal L_{m,r}(h)=0.
\tag{4.14}
\]

Equation (4.14) is therefore a literal polynomial-algebraization condition:
it cancels the possible pole left by the rational chart, not a remaining
Jacobian zero.

The formulas (4.1)--(4.14) prove the theorem.  Their canonical sources are
the [tangent-map core](verified/TANGENT_MAP_CORE.md), the
[controlled-boundary ledger](cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md),
and the [cancellation construction](cancellation/CONSTRUCTION.md).  QED

## 5. The foundational cubic in both charts

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

This is the smallest example in which the two suspension types give the same
polynomial Keller map up to polynomial left--right equivalence.

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

In both charts the unique ramified boundary prime over the cubic
discriminant has ramification index two and carries a simple source pole.
The left--right equivalence (5.8), together with functoriality of the finite
normalization package, identifies these as the same intrinsic boundary
valuation.  What changes is the ledger used to cancel it: polynomial
vertical weights in (5.3), a reciprocal rational chart in (5.7).

The complete branch-collapse calculation is
[Theorem 4.1 of the minimal-boundary note](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md#4-cubic-collapse-inside-the-two-branches).

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
[fixed-rank Dixmier reduction](extended-geometry/FIXED_RANK_DIXMIER_REDUCTION.md),
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
Zariski--Main package.

### Minimal-boundary suspension conjecture

Let

\[
F:\mathbb A^3\longrightarrow\mathbb A^3
\]

be a nonproper Keller map whose intrinsic finite-normalization package has:

1. one geometrically integral rational critical-boundary normalization;
2. saturated affine links with no undeclared zero or pole support;
3. boundary-monotone positive links;
4. at most two places at infinity on the critical normalization; and
5. a divisor-minimal ledger with no unrecorded graph valuation.

Then `F` is polynomially left--right equivalent to either a weighted tangent
suspension or a cancellation suspension.

This is the existing
[minimal-boundary classification conjecture](cancellation/MINIMAL_BOUNDARY_CLASSIFICATION.md),
restated here to show its role in the thesis.  It is deliberately not a
conjecture about every arbitrary Keller map without hypotheses.

The first extraction theorem to prove is:

\[
\boxed{
\text{intrinsic minimal boundary package}
\Longrightarrow
\text{coordinate-preserving controlled suspension square}.
}
\tag{7.2}
\]

Once the square is extracted, the number of places on the normalized
critical curve separates the expected branches:

\[
\mathbb A^1
\rightsquigarrow
\text{weighted tangent core},\qquad
\mathbb G_m
\rightsquigarrow
\text{reciprocal cancellation core}.
\tag{7.3}
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
compare two cubic answers, but to prove either intrinsic cubic gateway:
suspension extraction, or flat binary-cubic normalization followed by
coefficient-gauge straightening.

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
