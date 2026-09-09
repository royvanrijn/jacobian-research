# The first302 seed survives ideal-parity normalization

**Latest refinement.** The Artin calculation below proves seven split
`Z/2` ideal-class factors already represented by inherited virtual units.
The seed supplies a seventh unramified character, but no detected increase
in ideal-image rank. Its unit-versus-new-ideal-class ambiguity reduces
exactly to two explicit integral cubic norm equations. Both are locally
soluble everywhere; neither integral equation has been solved.

The subsequent [common-quadric reduction](DET1092_SEED_NORM_QUADRIC_AND_INTEGRAL_GATE_2026-09-09.md)
gives explicit rational parametrizations while retaining the two integral
lattices and all chart exceptions. The rational norm surfaces are isomorphic
regardless of principality; this simplification does not solve either equation.

**Bank-only obstruction.** A subsequent [exact relation audit](DET1092_SAVED_UNIT_RELATION_OBSTRUCTION_2026-09-09.md)
proves that arbitrary multiplicative words in all567 saved reduction
multipliers, even with rational rescaling, yield no units beyond `+/-1`.
This closes an overlooked-relation possibility in that bank, not either
norm equation or the unit-versus-ideal-class ambiguity.

## Result and boundary

**Verified application.** The historical first302 seed has the same
bad-prime ideal-parity footprint as a product of eight inherited classes.
Ideal-parity rank is9 before and after adjoining it, but full cubic
squareclass rank rises from17 to18. Valuation parity alone loses this seed.

**New exact construction.** Multiplying by that inherited product and an
explicit square gives an integral cubic element `beta` with

\[
 [\beta]\notin H,\qquad (\beta)=J^2.
\]

Here `H` is the generic Kummer image, not its ideal-class image. Exact
coefficients, ideal bases and square multipliers are certified below.
This is a *virtual unit*: all finite prime-ideal valuations are even.
It need not be an actual unit, a locally trivial squareclass, or an
unramified quadratic extension at primes above2.

**New lower-bound deduction.** The associated totally real cubic field has

\[
 \boxed{\dim_{\mathbf F_2}\operatorname{Cl}(K)[2]\ge7.}
\]

No class group or unit group was computed. This does not prove that this
individual `J` is nonprincipal or outside the inherited ideal-class image.
The construction is seed-derived, not a prospective seed selector or a
null-panel discriminator. That boundary remains open.

## Exact inputs

**Verified application.** Use the cubic and certified maximal order from
the immutable [constructor arithmetic](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_strict_constructor_arithmetic_v1.json).
In the original302 model, `X=4x`, `Y=8y+4x+4` and `Y^2=f(X)`. Put
`K=Q(theta)`, where

```text
f(z) = z^3 + 5*z^2 + b*z + c
b = -20555644225817083652508762176252556301867322902653588194515587713000
c = 35863572573072725712438997689722579989297270955836506038914807628765229483647299767251594888128395600
```

The cubic is irreducible modulo31 and has positive discriminant, so `K`
is totally real. Its defining-order index is
`9923258571576383661788160`; calculations use the certified maximal order,
not the nonmaximal power order.

**Retrospective input.** The historical first seed in the
[marked transport certificate](../../artifacts/generated-results/elliptic-curves/det1092_marked_kummer_transport_v1.json)
has integral representative

\[
 \alpha_*=
 3192654717569013009686535431767390554676-1234321\theta.
\]

Its norm is the square of
`97650149700983343619818748772115777389811711955964311816`.
This is not a later cascade point or the different autonomous first point.
Number the17 inherited integral representatives `alpha_0,...,alpha_16`
in the exact order of the constructor arithmetic artifact.

## What the valuation map sees

**Verified application.** Record valuations modulo2 at the40 prime ideals
above the20 already frozen bad rational primes. The odd valuations of
`alpha_*` occur precisely at

```text
(23,1), (41,0), (41,1), (73,1), (131,1), (167,1).
```

Prime indices are zero-based in the certified decomposition; their exact
HNF ideals are retained. The generic control `alpha_0` also has odd
valuations, at `(7,1),(11,1),(19,1),(23,1),(37,1)`.

**Verified linear algebra.**

| Space | Squareclass dimension | Ideal-parity rank | Even-valuation kernel |
|---|---:|---:|---:|
| Generic `H` |17|9|8|
| `H + <alpha_*>` |18|9|9|

A canonical linear solve, without support optimization, gives

\[
 W=\{0,1,2,4,5,7,9,10\},\qquad
 v(\alpha_*)=\sum_{i\in W}v(\alpha_i).
\]

Thus `alpha=alpha_* product_(i in W) alpha_i` has even valuations
everywhere. This is global, not a finite-prime sample: for every integral
representative `a` with `N(a)=y^2`, form `I=(a,y)`. Direct integral-lattice
calculations verify that `I^2/(a)` is integral and its norm is supported
only at the frozen bad primes. Outside that set `(a)=I^2`. Exact
bad-prime valuations finish the argument, without factoring the norm.

All49 frozen good-prime character blocks remain usable; no replacement
prime is selected. Independent finite-field square tests prove ranks17
and18. The cancellation word therefore cannot make the seed generic.

## The explicit square ideal

**Retrospective construction, exact verification.** Correct the product
half-ideal at the known bad primes to obtain `J_0^2=(alpha)`. On
`J_0^{-1}` use the positive definite form `Tr(u*v)` and take just the
first vector `gamma` of one exact three-dimensional LLL reduction. Set
`beta=alpha*gamma^2` and `J=(gamma)J_0`. The result is

```text
beta =
  5415611000709889833034620814517180079747088650548440003575813129631494971/224486984292346776
 -1143366181533958899624861845000515957026633457/62020366072352397886176 * theta
 +1092016605499/310101830361761989430880 * theta^2

N(J) = 4041308965738402063503944863521386373242705642959299925097038321747198596
```

Despite its power-basis denominators, `beta` is integral in the maximal
order. The [certificate](../../artifacts/generated-results/elliptic-curves/det1092_seed_half_ideal_v1/virtual-unit.json)
contains `gamma`, both exact ideal bases, and the identities
`beta=alpha*gamma^2`, `J_0^2=(alpha)` and `J^2=(beta)`.

**Explicit limitation.** This reduction did not produce a unit generator.
A reduced ideal different from the unit ideal does not prove
nonprincipality. This is not a shortest-representative or low-cost chart
certificate. The earlier two-class weighted-trace diagnostic is retained,
including nonintegral outputs and the `alpha^2` fallback when the
`Tr(alpha*u*v)` form was indefinite.

## Class-group consequence

**Established literature.** For virtual units
`V={a in K*/K*2 : all finite valuations even}`, the standard sequence is

\[
 0\longrightarrow\mathcal O_K^*/\mathcal O_K^{*2}
 \longrightarrow V\longrightarrow\operatorname{Cl}(K)[2]
 \longrightarrow0,
 \quad [a]\longmapsto[J],\quad(a)=J^2.
\]

See [Schaefer--Stoll, section7, the sequence for D(S,p)](https://www.mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf),
with `S` empty and `p=2`. No novelty is claimed for this sequence.

**New application; proof.** All18 displayed classes have square rational
norm. Their9-dimensional valuation kernel lies in `V` with square norm.
A totally real cubic has unit rank2 and torsion `{+1,-1}`. Since
`N(-1)=-1`, norm-positive unit squareclasses have dimension exactly2.
The kernel's image in `Cl(K)[2]` consequently has dimension at least
`9-2=7`.

In general, `r` independent norm-square classes in a totally real cubic,
whose valuation-parity map has rank `v`, prove
`dim Cl(K)[2] >= r-v-2`. Here the17 inherited classes alone give the
weaker guaranteed bound6. These are lower bounds, not a claim that the
actual class-group dimension changes on adjoining a point.

## Meaning and remaining obstruction

**New deduction.** An ordinary class-group bridge becomes meaningful
after cancelling the inherited ideal-parity footprint. Its possible unit
kernel cannot be discarded: the individual `[J]` might still be inherited
or trivial, even though `[beta]` is not generic.

**Verified distinction.** Even valuations do not imply the older strict
local conditions vanish. The first seed is already certified outside
`H + K_strict`; generic multiplication preserves that property. This
virtual unit is therefore still non-strict modulo `H`. Moreover it is a
Kummer class of a rational point (seed plus inherited points), so its
Sha image is zero. It is not a newly inferred Sha class.

**Open obstruction.** Determine whether `[J]` is outside the ideal-class
image of the eight-dimensional inherited virtual-unit space, or whether
its nongeneric information can be carried by a unit after an inherited
correction. Then produce the relevant class from equations without
`alpha_*`. Neither step is achieved here. No candidate-free invariant
has been shown to distinguish302 from the eight null fibres, and no
amplification prediction follows from the class-group lower bound.

## Checkpoints and independent replay

**Verified computation.** Each process completed in under one second,
with a25-second hard cap. The independent checker imports no constructor
and makes no PARI number-field or ideal calls. It uses rational cubic
multiplication, integer Hermite forms and finite-field square tests to
check720 exact prime-ideal valuations, global defect support, square-ideal
identities, squareclass transport and ranks. The earlier certified
maximal order and prime decomposition are explicit dependencies.

The first replay compared different HNF conventions directly and failed;
canonicalizing both ideal bases corrected the checker. No mathematical
identity or construction output changed.

- [All protocols and arithmetic](../../artifacts/generated-results/elliptic-curves/det1092_seed_half_ideal_v1/)
- [Independent replay result](../../artifacts/generated-results/elliptic-curves/det1092_seed_half_ideal_v1/independent-replay.json)
- [Independent checker](../cas/verify_det1092_seed_virtual_unit.sage)

```bash
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_seed_virtual_unit.sage
```

**Scope.** No later302/rank21 point, V3 artifact, new parameter, point
search, unit/class-group calculation, unrestricted factorization,
production change or detached process was used.

## Artin refinement: characters are not ideal classes

**Audit of completed work.** The older
[constructor checkpoint](../rank-jump/CURVE302_CONSTRUCTOR_TRANSFER_CHECKPOINT.md)
already tested567 fixed archimedean reductions and obtained only units
`+/-1`. That experiment was not rerun or enlarged. The old
[generic-character certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_only_class_anchors_v1.json)
provides six ordinary unramified quadratic characters on302. They come
from generic sections only. A new, separately frozen check of the eight
canonical generic valuation-kernel basis ideals used one trace-LLL vector
per ideal; it produced no unit relation or equal reduced-ideal pair.
These failures remain bounded, not nonprincipality proofs.

**Established literature.** A quadratic extension unramified at all finite
places and split at all real places gives a character of the ordinary
ideal class group by Artin reciprocity. See
[Milne, Class Field Theory, Theorem0.4 and the Hilbert class field](https://www.jmilne.org/math/CourseNotes/CFT.pdf).
Such a character kills every principal ideal, but need not detect every
nonprincipal ideal. Quadratic characters in particular kill `2Cl(K)`.

**New exact construction, retrospective.** Append the seed to the17
generic classes and impose even finite valuations, ordinary dyadic
unramifiedness, and positivity at the three real embeddings. Constraint
rank stays11, so the unramified-character dimension grows from6 to7.
The additional character is represented by

\[
 a_7=\alpha_*\prod_{i\in\{0,4,5,6,7,9,10\}}\alpha_i.
\]

Its complete cubic coefficients and all local tests are in
[characters.json](../../artifacts/generated-results/elliptic-curves/det1092_seed_artin_v2/characters.json).
At each dyadic completion, either `a_7` or `5a_7` is a square, proving
ordinary unramifiedness there; all other finite valuations are even and
all three real signs are positive. The six older characters are retained
in their original order. No claim that the seventh was selected without
the seed is made.

**Verified Artin evaluation.** Let `G_0,...,G_7` be the eight inherited
half-ideal classes from the canonical generic valuation-kernel basis,
and `J_*` the seed-normalized ideal of the preceding section. The matrix
below has these nine ideal classes as rows and the seven unramified
characters as columns, with `1` denoting Artin sign `-1`:

```text
G0   0 0 1 1 1 0 1
G1   1 1 1 0 1 1 1
G2   1 1 0 0 1 0 1
G3   0 0 0 0 1 1 0
G4   0 1 1 0 1 1 1
G5   1 0 1 0 0 1 1
G6   0 1 1 1 1 0 0
G7   1 1 1 1 1 1 1
J*   0 1 0 0 1 0 1
```

There is no norm factorization. For each ideal, remove only its exact
factors above the20 already frozen bad primes. Every remaining ideal has
an odd cyclic quotient `O_K/I=Z/N`, on which all18 marked factors are
units. Their Artin product is therefore an ordinary rational Jacobi
symbol. Local symbols at the removed primes supply the remaining factors.
All nine cases pass the fixed quotient/unit gate, without retries.

**New deduction: split factors, not just a lower bound.** The first seven
rows form an invertible seven-by-seven matrix over `F_2`. Each `G_i` has
order dividing2 because its square is principal. Hence `G_0,...,G_6`
generate a subgroup `C` isomorphic to `(Z/2)^7`, and the seven characters
give a retraction onto `C` after inverting this matrix. Consequently

\[
 \operatorname{Cl}(K)\simeq(\mathbf Z/2\mathbf Z)^7\oplus A
\]

for an unspecified finite abelian group `A`. This is stronger than merely
`2-rank >=7`; these seven split factors are represented by inherited
ideals. The seventh *character* is seed-derived, but its detection of
those ideals does not make the ideals exceptional.

The eighth inherited row and the seed row both lie in the span of the
first seven. Thus the measured ideal-image rank remains7. Vanishing in
this finite character quotient does not imply a unit representation.

**Additional existence consequence, not a constructor.** The generic
virtual-unit space has dimension8 and its Artin image has rank7. Its
intersection with norm-one units modulo squares has dimension at most1.
The full norm-one unit squareclass space has dimension2. Therefore at
least one norm-one unit squareclass lies outside `H`. No such unit is
explicitly constructed here, and this does not prove it is a point-Kummer
class, Selmer class, or soluble elliptic covering.

## Exactly two remaining principality questions

**New exact reduction.** All Artin-compatible inherited corrections of
the seed are the following two words in the eight `G_i`:

\[
 w_0=G_1+G_2+G_4,\qquad
 w_1=G_0+G_1+G_2+G_3+G_5+G_7.
\]

Their difference is the unique inherited Artin-kernel word
`G_0+G_3+G_4+G_5+G_7`. This exhausts the possibilities by exact linear
algebra, not by enumerating256 combinations. Form and reduce once the
two relative ideals `R_j=J_* / product(G_i in w_j)`. Exact multipliers
are retained, and the resulting integral representatives satisfy
`R_j^2=(b_j)` and `N(b_j)=N(R_j)^2`.

Their norms are

```text
N(R0) = 14359835166045906374689557773242739873751525483578816266508557674617313784
N(R1) = 26293910318995873645829948921293546250151485289367390394100625408239040903
```

For the saved integral basis `e_0,e_1,e_2` of each `R_j`, define

\[
 F_j(u,v,w)=\frac{N_{K/\mathbf Q}(ue_0+ve_1+we_2)}{N(R_j)}.
\]

The ten coefficients are all integers; they and the bases are fully
expanded in [norm-form-00.json](../../artifacts/generated-results/elliptic-curves/det1092_seed_artin_v2/norm-form-00.json)
and [norm-form-01.json](../../artifacts/generated-results/elliptic-curves/det1092_seed_artin_v2/norm-form-01.json).

\[
 \boxed{[J_*]\in\langle[G_0],\ldots,[G_7]\rangle
 \ \Longleftrightarrow\
 F_0(\mathbf Z^3)\ni1\ \text{or}\ F_1(\mathbf Z^3)\ni1.}
\]

Indeed, an element of `R_j` with norm `N(R_j)` generates it, by equality
of ideal indices. Conversely, a principal ideal has a generator with
absolute norm `N(R_j)`; change its sign if necessary, since the field
degree is odd. This is an exact two-equation criterion, not a solution
of either equation. Both single reductions returned nonunit ideals.

## Why rational or local points cannot resolve this fork

**New deduction, with explicit witnesses.** Every integral ideal `R` in a
cubic field with `R^2=(b)` and `N(b)=N(R)^2` has the rational norm witness

\[
 a=N(R)/b,\qquad N(a)=N(R).
\]

Its coordinates in the basis of `R` give an explicit rational point on
`F=1`. Both witnesses are saved and exactly substituted into the forms.
The condition of interest is integer coordinates, not rational coordinates.

**Proof of everywhere local integral solubility.** Fix a rational prime
`p`. If `K tensor Q_p` is a cubic field with ramification and residue
degrees `e,f`, and `v_P(R)=m`, then
`v_P(N(R)/b)=ef*m-2m=m`, because `ef=3`. The explicit rational witness
therefore belongs to `R tensor Z_p` and generates it locally.

Otherwise the cubic algebra has a `Q_p` factor and a complementary
degree-two etale algebra. Choose a generator of the ideal in each
complementary factor. Set the `Q_p` coordinate equal to `N(R)` divided
by their product norm. Its valuation is exactly the remaining ideal
valuation, by the ideal-norm identity. This supplies a local integral
solution of `F=1`. The argument includes the fully split case. The real
equation already has the displayed rational point.

Thus both norm equations have points over every `Z_p` and over `R`,
regardless of whether their ideals are principal. No finite local
solubility obstruction can settle these two principality questions.
This is an obstruction to that diagnostic, not a proof of global
insolubility or a Selmer/Sha calculation for the elliptic curve.

## Refined checkpoint and boundary

**Verified computation.** The independent
[Artin/norm-form replay](../cas/verify_det1092_seed_artin.sage)
uses local square tests in place of the constructor's Hilbert symbols,
its own integer Jacobi algorithm, and rational cubic multiplication with
integer Hermite forms in place of PARI ideal multiplication. It checks
the full matrix, the invertible minor, both relative class identities,
both expanded norm forms, and their rational points. Runtime is under
one second, under a25-second cap; Sage10.9 and PARI2.17.3 were used.

The first Artin script stopped on a Sage matrix-versus-vector coercion
before evaluating any Artin column. Its source and failure record are
preserved under `det1092_seed_artin_v1`; v2 changes that interface call,
not the mathematical rule or limits. All subsequent fixed cases pass.

**Open boundary.** The seed's new unramified character is explicit, but
its selection still uses the first seed. The two integral norm equations
are unresolved. No prospective nongeneric rational point, control-panel
discriminator, additional strict class, or amplification mechanism has
been inferred. No running search, class-group/unit-group computation,
parameter scan or point search was started.
