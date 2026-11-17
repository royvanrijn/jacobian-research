# Logarithmic target counts on the elliptic exceptional slice

The fixed Keller map has an exceptional noncentered slice on which the
matching-resolvent target count grows logarithmically. This count concerns
targets with matching quadratic and cubic resolvents; the additional local
conditions for Hasse failure remain separate.

Retain the slice `(pi,ell)=(128/25,9/10)` and its target map `y(v)` from
[the resolvent-fibration note](FIXED_MAP_RESOLVENT_FIBRATION.md), equations
(1.2)--(1.3). Its normalization is

\[
 E:\quad T^2=X^3-10758382X+6600759744,
 \qquad v=\frac{3(X+3984)}{6400}.                    \tag{1}
\]

The admissible open excludes `D Delta Res(q,h)=0`, with
`D=ell^2-4v`; these are finitely many parameters. The two signs of `T`
give the same parameter and target. The target height is always the
actual projective height `H(y)=H([1:y_1:y_2:y_3])`.

The [arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md) proves that
`E(Q)` has rank two and torsion of order two, that all admissible
quadratic factors are irreducible, and that only finitely many admissible
cubic factors are reducible. It also proves that the possible quadratic
fields are exactly

\[
 \mathcal K=\{\mathbb Q(\sqrt{-3}),\mathbb Q(\sqrt{-6}),
    \mathbb Q(\sqrt{-40637595}),
    \mathbb Q(\sqrt{-81275190})\}.                  \tag{2}
\]

Those rank and finiteness statements are proof inputs here; they are not
inferred from a numerical rank calculation or a sampled target count.

## Theorem EEC1: a leading asymptotic for matching-resolvent targets

Normalize the canonical height by

\[
 q_x(P)=\lim_{n\to\infty}4^{-n}
                  \log H(X([2^n]P)),               \tag{3}
\]

with height one at the identity. Thus `q_x` differs from `log H(X)` by
a bounded amount; it is twice the alternative normalization that uses
one half of the logarithmic `X`-height. If `P_1,P_2` are any integral
basis of the free group `E(Q)/E(Q)_tors`, let `M` be the matrix of its
quadratic form:

\[
 M_{ii}=q_x(P_i),\qquad
 M_{12}=M_{21}=
 \frac{q_x(P_1+P_2)-q_x(P_1)-q_x(P_2)}2,
 \qquad R_x=\det M>0.                              \tag{4}
\]

The determinant `R_x` is independent of the choice of integral basis.
No assertion that the previously displayed points form such a basis is
needed. Write `pi_circ` for the circle constant, to distinguish it from
the fixed inverse parameter `pi=128/25`. Let `M_E(B)` count distinct admissible targets on this slice with
matching squareclasses and height at most `B`. Then

\[
 \boxed{M_E(B)=\frac{\pi_{\rm circ}}{3\sqrt{R_x}}\log B
                       +O_E(\sqrt{\log B})}.        \tag{5}
\]

The same asymptotic holds if one requires irreducibility of both factors
and hence a full rank-five quadratic-times-cubic field algebra.
Write `C_E=pi_circ/(3 sqrt(R_x))` for its leading coefficient.

**Proof.** The exact height comparison in Proposition 3.3 of the
fibration note gives

\[
 \log H(y(v))=3\log H(v)+O_E(1).
\]

The affine linear change in (1) gives `log H(v)=log H(X)+O_E(1)`.
The canonical-height theorem therefore yields

\[
 \log H(y(P))=3q_x(P)+O_E(1).                       \tag{6}
\]

The exact imported interface is J. S. Milne, *Elliptic Curves* (2006),
Chapter IV, Proposition 4.5 and Theorem 4.7, printed pp.121--122, for
the bounded difference from logarithmic `X`-height, and Theorem 6.1,
printed p.131, for positive definiteness of the height pairing:
[author-hosted book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The normalization in (3) is the one used there. The associated pairing
is twice the matrix in (4), so no hidden factor of two is absorbed into
`R_x`.

Write a rational point as a torsion point plus `n_1P_1+n_2P_2`.
Torsion does not change `q_x`, so

\[
 \#\{P\in E(\mathbb Q):q_x(P)\le R\}
       =\frac{2\pi_{\rm circ} R}{\sqrt{R_x}}+O_E(\sqrt R+1).    \tag{7}
\]

Indeed the ellipse `n^t M n<=R` has area `pi_circ R/sqrt(R_x)`.
Covering its boundary by the fixed unit fundamental squares gives
`O_E(sqrt R+1)` error, and there are two torsion translates. By (6),
the height condition lies between ellipses of radii squared
`(log B)/3+O_E(1)`, which have the same main term and stated error.

Except at finitely many parameters, the target map has precisely the
two preimages `P,-P`. The map `E -> P^1_v` has degree two. The rational
map from the parameter line to its target cubic curve is birational:
`b(v)` and `c(v)` have degrees two and three, so the degree of
`Q(v)/Q(b,c)` divides both. The exceptional identifications of distinct
parameters are finite, as for any finite birational normalization of a
projective integral curve. Removing the finite inadmissible parameters
and dividing (7) by two proves (5).

All admissible quadratics are irreducible, and the arithmetic note's
genus-four root cover and Faltings argument leave only finitely many
reducible cubics. Requiring the rank-five field product therefore
changes the count by `O_E(1)`. \(\square\)

## Theorem EEC2: equal leading counts for the four quadratic resolvents

For each `K` in (2), let `M_{E,K}(B)` count the full rank-five matching
targets of height at most `B` with quadratic field `K`. Then

\[
 \boxed{M_{E,K}(B)=\frac{\pi_{\rm circ}}{12\sqrt{R_x}}\log B
                         +O_E(\sqrt{\log B})}.      \tag{8}
\]

In particular each of these four fields occurs in one quarter of the
matching targets asymptotically. This distribution concerns this one
fixed slice, not the entire target space.

**Proof.** The descent homomorphism in the arithmetic note has exactly
four squareclasses. Its fibres are four cosets of an index-four subgroup
of `E(Q)`, and the quadratic field is exactly `Q(sqrt(-3 delta(P)))`.
Every such coset has one quarter of the main term in (7), with the same
boundary-error order. To see this without choosing a saturated Mordell--Weil
basis, project the kernel to the free lattice and retain its finite
torsion intersection: the product of the lattice index and the torsion
quotient index is four. Translates of the resulting full lattice have
the same area and boundary estimates. Negation preserves each descent
class, so dividing by the two geometric preimages of a target preserves
the quarter proportion. Finite exclusions and cubic-root exceptions do
not change the main term. \(\square\)

For every fixed finite set of places `S` and every quadratic
squareclass `a`, [EFA1--EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md)
give the full finite-place asymptotic

\[
 N_{S,a}(B)=\mu_{S,a}C_E\log B+o_{E,S}(\log B).
\]

Here `mu_{S,a}` is normalized Haar measure on the actual joint closure
of `E(Q)`, including its squareclass coordinate. It retains the
correlations between places and lies between zero and `1/4`.
For the certified class `a_0=10*2709173`, with
`K_0=Q(sqrt(-81275190))`, it is positive for every fixed finite `S`.
At `S={2,3,17}`, [J2174](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md)
gives the exact value `mu_{S,a_0}=1/8`.

## Corollary EEC3: bounds and the remaining local condition

On the matching-resolvent slice, globally soluble targets form a finite
set. Its Hasse-target count `A_E(B)` satisfies, for all sufficiently large `B`,

\[
 1\le A_E(B)\le \frac{2081912199403}{9893175318390}C_E\log B
                         +O_E(\sqrt{\log B}).       \tag{9}
\]

More precisely the everywhere locally soluble count on this restricted
family is `A_E(B)+O_E(1)`. No leading asymptotic or infinitude of
`A_E(B)` follows from the matching-target asymptotic.

For the certified quadratic field `K0=Q(sqrt(-81275190))`, one also has

\[
 A_{K_0}(B)\le\frac{13361}{114520} C_E\log B
                         +O_E(\sqrt{\log B}).       \tag{10}
\]

The [complete six-prime prefix at two, three, five, eleven, seventeen and fifty-three](ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_53_SIEVE.md),
J533, gives (9)--(10), using the square-root error theorem ESQ3.
Fifty-three multiplies the previous certified-class coefficient by
`13361/14315` and the class-one coefficient by
`312824741/312825695`; the other two coefficients are unchanged.
J531 proves five-saturation of the known rank-two subgroup from
the unit character matrix at eleven. Together with the existing
two-/three-saturation, this identifies the full finite image at
fifty-three without asserting a displayed Mordell--Weil basis.
J532 retains its correlations with the old eleven- and seventeen-
adic tests explicitly.

The earlier [five-prime prefix at two, three, five, eleven and seventeen](ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md)
retains its exact coefficient `116302241/531619620` for that set.
Write `q_0=2709173`. Eleven multiplies its preceding class-one and class-`5q_0` coefficients by
`533/665`, giving `6724/8860327` and `533/55860`; the other two
coefficients remain unchanged. The earlier
[four-prime prefix at two, three, five and seventeen](ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md),
J235173, retains its exact coefficient `2300089/10392564` for that
specified set. Its class-one coefficient is `820/866047`, obtained
by multiplying the three-prime coefficient by `16/31`. The direct
finite-quotient certificate J2172 still gives the valid, weaker
coefficient `1687/6936` without the logarithm input.

For comparison, the [complete three-adic classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md)
gives coefficient `55/156` for solubility at three alone. The
[two-adic and joint classification](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md)
gives `4/7` at two alone and `289/1092` at two and three together.
The joint `2,3,17` coefficient is `5321/23946`. Each of these is an
exact leading coefficient relative to `C_E log B`, now with
`O_E(sqrt(log B))` error by ESQ3, for its specified finite local set. These proportions cannot
be multiplied as independent events and do not count all-places
Hasse failures.

[EFA3](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md) identifies the
all-primes limit of these finite-place measures and the exact additional
tail condition needed to obtain its coefficient in the Hasse count.
For the certified class the decreasing limit is between zero and
`13361/114520`; if it is zero, the class's Hasse count is `o(log B)`.
Its value and, in the positive case, the required tail control remain
open. The total-class version has the same distinction. Thus the full
finite-`S` asymptotics do not establish Hasse infinitude or an all-places
leading asymptotic.

[EQA2](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md)
proves that every fixed finite-place coefficient is rational.
ESQ3 gives error `O_(E,S)(sqrt(log B))` for each such count, while
[EPE2](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md) erases automatic
prime tests exactly. The elementary EQA3 bound retains its separate
proof without the logarithm input. A finite-local field
target set is either empty or has a positive rational coefficient,
by EQA2a. These statements retain their fixed finite-prime scope.

These counts use the actual joint rational closure throughout.

[ECO1](ELLIPTIC_EXCEPTION_FINITE_LOCAL_CLASS_ONE.md) proves that the
class `Q(sqrt(-3))` has a positive coefficient for every finite `S`,
using formal neighbourhoods rather than an all-places seed. Together
with the certified class, two of the four quadratic fields therefore
have this finite-place positivity property. The
[triangular prime sieve](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md)
now supplies an exact conditional hazard at every formal depth and
a weighted pathwise bound. The arithmetic weighted-return divergence
and all-primes tail remain unresolved. Separately,
[EOP1--EOP3](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md) prove that all
four classes pass every individual prime, and that every feasible
finite set admits all but finitely many one-prime extensions. These
facts do not settle simultaneous positivity for all finite sets in
the other two classes. The same result proves a uniform absolute
coefficient loss `O_E(1/log p)` when a large new good prime is added,
and `mu_({p},a)=1/4+O_E(1/log p)`. It does not give a uniform
target-counting error or a summable all-primes tail.

[The Kummer branch-exclusion theorem](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md)
identifies a positive-density set of eligible primes where the
numerator is prime to that prime for every rational point on the
slice. Only pole returns contribute to the local hazard there.
This restricts the actual prime support without changing the
finite-prefix upper coefficient in (9) or deciding the all-primes
measure.

[EFO2--EFO3](ELLIPTIC_EXCEPTION_FINITE_OBSTRUCTION_CERTIFICATES.md)
give exact finite congruence-cover and fatal-branch criteria for
simultaneous infeasibility, under the certified-quotient input
hypotheses. They do not exhibit an obstruction in the two unresolved
classes. [ECM1--ECM3](ELLIPTIC_EXCEPTION_COMMON_QUOTIENT_MOMENTS.md)
compute the relevant finite correlations and give a sufficient
conditional mean-variance criterion for zero all-primes measure.
Its asymptotic moment hypotheses remain unproved; the upper bounds
(9)--(10) therefore retain their present scope.

[ESR1--ESR2](ELLIPTIC_EXCEPTION_SMALL_REDUCTION_IMAGES.md) strengthen
the marginal loss estimate: for every fixed `A>0`, outside a set of
`O_(E,A)(X/(log X)^(2A+1))` primes up to `X`, the absolute coefficient
loss on adding `p` is at most `3(log p)^A/sqrt(p)`, uniformly over
the finite prefix and quadratic class. A product of denominators of
height-bounded difference points proves the exceptional-prime bound.
This does not bound relative losses for shrinking prefixes or the
conditional covariances in ECM3.

[The finite dependency core](ELLIPTIC_EXCEPTION_DEPENDENCY_CLOSED_CORE.md)
removes the unspecified initial cutoff from feasibility in the two
remaining classes. The core is exactly `2,3,5,1889`, with concrete
full-field witnesses. A fixed permitted prime pool of density
`13/16+delta_Kum` has the following property: every finite subset
retains a positive rational logarithmic target coefficient in all
four quadratic fields. Each remaining class separately has a pool
of density `7/8+delta_Kum`. There is no coefficient lower bound
uniform in the finite subset.

[The rank-three branch estimate](ELLIPTIC_EXCEPTION_UNIFORM_OBSTRUCTION_PRIMES.md)
improves numerator-return masses to `2(log p)^A/p^(3/5)` on a set
of prime density one, for every fixed `A>0`. It gives an effective
finite candidate set for any new prime killing a permitted prefix
of mass at least `eta`, with weighted size `O_E(eta^(-5/3))`.
The numerator marginals are square summable over all primes, with
an explicit logarithmically decaying tail. That tail concerns two
independent Haar vectors; it does not give the single-vector tail
or a rational-pair counting asymptotic.

[The pole variance obstruction](ELLIPTIC_EXCEPTION_POLE_VARIANCE_OBSTRUCTION.md)
shows why unconditional decorrelation cannot replace the actual
prefix estimates: every nonnegative unconditional pole sum has
variance at least one seventh of its squared mean. Affine-prefix
criteria retain their compatibility signs and remain conditional.

[The complete bounded depth certificate](ELLIPTIC_EXCEPTION_BOUNDED_FATAL_PRIME_CERTIFICATE.md)
proves that every finite local set supported on primes at most
`100000` has a positive coefficient in all four classes. The point
`Q` gives formal depth one at every one of the 1722 potentially
fatal eligible primes in that range; this is not an unbounded
depth theorem.

[Affine-coset pole means](ELLIPTIC_EXCEPTION_AFFINE_COSET_POLE_MEANS.md)
prove that every fixed feasible finite-place prefix has expected
unweighted pole returns at least a positive constant times
`log log Y`, where `Y` is the prime cutoff. This does not bound
the target-height counting error, the growing-prefix loss or the
conditional variance. Its construction uses the exact
[three-torsion disjointness](ELLIPTIC_EXCEPTION_THREE_TORSION_DISJOINTNESS.md).
The [class-group bridge](ELLIPTIC_EXCEPTION_WIEFERICH_CLASS_GROUP_BRIDGE.md)
turns deep formal images into specific equivariant unramified
division-field class-group quotients; no distribution estimate for
those quotients is asserted.

Under its stated good-prime and bad-local exclusions and the image
hypothesis `Gal(Q(E[p])/Q)=GL_2(F_p)`, the
[strict Selmer converse](ELLIPTIC_EXCEPTION_CLASS_GROUP_SELMER_CONVERSE.md)
identifies the level-one class-group maps with Selmer classes whose
localization at `p` is zero. A second `E[p]` quotient copy occurs
exactly when `e_p>1` or `Sha(E/Q)[p]!=0`. No finiteness or density
claim for `Sha` is assumed or proved. The
[adelic image theorem](ELLIPTIC_EXCEPTION_ADELIC_OPEN_IMAGE.md)
supplies those image hypotheses outside a fixed finite prime set.
The [multiplicity refinement](ELLIPTIC_EXCEPTION_CLASS_GROUP_MULTIPLICITIES.md)
then recovers `dim Sha[p]=2 floor((c_p-1)/2)` from the actual
class-group quotient multiplicity `c_p`. A third copy occurs exactly
when `Sha[p]!=0`; exactly two copies certify deep formal image and
zero `Sha[p]`, without assuming Sha finite.

[Finite-primary separation](ELLIPTIC_EXCEPTION_FINITE_PRIMARY_SEPARATION.md)
removes any prescribed finite set of primary correlations after
refining a congruence cell and selecting a positive-density prime
set. Its abstract countermodel shows that this separation and
divergent fixed-cell means still do not prove recurrence: large
common-quotient clustering remains uncontrolled. The
[fixed-character theorem](ELLIPTIC_EXCEPTION_FIXED_CHARACTER_DENSITIES.md)
now gives arithmetic bounds for each fixed Fourier mode and, after
forgetting finitely many exceptional primary coordinates, a weak
limit of the normalized reweighted pole measures equivalent to Haar.
It does not bound the complete finite-time Fourier energy or prove
uniform integrability. The abstract countermodel also has Haar weak
limit, so this averaged result alone still does not imply recurrence.

Beyond this slice, [the fixed-plane local count](FIXED_MAP_PLANE_LOCAL_COUNTS.md)
gives `O_P(B^2)` Hasse targets on every fixed rational affine plane.
The two section planes `u=0,c=0` are exactly those with order `B^3`
locally soluble targets. The constants depend on the plane, so this
does not supply a uniform sum over moving planes or an ambient
Hasse asymptotic.

The [ambient Haar calculation](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md)
gives measure zero both to the all-soluble set and to the rational
closure in the full ambient product. The intrinsic all-primes
coefficient remains undetermined. EFA4 expresses its positivity
exactly as summability of the actual conditional prime losses.

**Proof.** The quadratic is never split on the admissible family. A
rational source point therefore requires a rational root of the cubic,
and the arithmetic note proves that there are only finitely many such
targets. Partition the locally soluble targets into the globally soluble
ones and the Hasse failures. The
[explicit certified Hasse target](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md)
gives the lower bound once `B` exceeds its height. The upper bounds are
J533, using the five-prime prefix, its correlated finite image at
fifty-three and the independent new formal coordinate, together
with ESQ3.
The rational local images are proved from two-descent,
finite reductions at `41,53`, and formal-depth witnesses. Their
correlations are retained; no full Mordell--Weil basis is assumed.
\(\square\)

These are counts of target points, without quotienting targets by
isomorphism of their cubic fields or fibre algebras. The genus-four
argument makes the number of targets representing any one fixed cubic
field finite; it does not give a uniform bound adequate for a distinct-field
counting asymptotic. The regulator constant is intrinsic and positive;
no numerical regulator, full generator basis or saturation is claimed.

The unrestricted counting proofs are the height reduction and elementary
lattice argument above. Exact rank, torsion, descent classes and branch
calculations are checked in the arithmetic note's retained checker; its
finite controls do not certify the canonical-height or Faltings imports.
No formal verification or external review is claimed.
