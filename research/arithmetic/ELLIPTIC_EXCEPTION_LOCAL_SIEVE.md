# Finite-place local sieves on the existing elliptic Hasse slice

The existing elliptic slice already supplies infinitely many full
quadratic-times-cubic field fibres that pass any prescribed finite set
of local tests. Such a subfamily has a positive logarithmic target
count. The unresolved step is simultaneous control of the new bad
primes that vary with the cubic field, so this does not yet prove
infinitely many everywhere locally soluble fibres on the slice.

The full set passing any fixed finite collection is counted by
[EFA1--EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md), with a
Haar-measure leading constant. The finite-index subfamilies below
retain their sharper square-root logarithmic error.

## 1. Fixed inputs and an exact local seed

Retain the elliptic curve, parameter and target map

\[
 E:\quad T^2=X^3-10758382X+6600759744,\qquad
 v=\frac{3(X+3984)}{6400},\qquad y=y(v)
\]

at `(pi,ell)=(128/25,9/10)` from
[the resolvent fibration](FIXED_MAP_RESOLVENT_FIBRATION.md).
The target depends only on `v`, so `y(P)=y(-P)`.
The [elliptic arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md) proves:

- `E(Q)` has rank two and torsion subgroup of order two;
- every admissible quadratic factor is irreducible;
- only finitely many admissible points have a reducible cubic factor;
- each fixed cubic-field isomorphism class occurs at only finitely many
  admissible targets.

Fix either of the two certified rational points `P_0` above

\[
 v(P_0)=\frac{39976591281}{36871680400}.               \tag{1.1}
\]

Its full fibre is Hasse-failing, by
[the exact elliptic example](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md),
and its quadratic field is `K_0=Q(sqrt(-81275190))`.
The [height count](ELLIPTIC_EXCEPTION_COUNT.md), Theorem EEC1, gives
the unrestricted matching-target coefficient

\[
 C_E=\frac{\pi}{3\sqrt{R_x}}>0,                     \tag{1.2}
\]

where `R_x` is the determinant of the canonical quadratic-height matrix
on an integral basis of the free Mordell--Weil group, in its explicitly
stated normalization.

## 2. Every finite set of places admits a logarithmic field subfamily

**Theorem ELS1 (finite-place local sieve).** For every finite set `S`
of places of `Q`, there is a finite-index subgroup
`Gamma_S subset E(Q)` such that the admissible targets from

\[
 U_S=(P_0+\Gamma_S)\cup(-P_0+\Gamma_S)              \tag{2.1}
\]

are locally soluble at every place of `S`. The subgroup can be chosen
inside the kernel of the squareclass map in Theorem 2.1 of the
arithmetic note; all admissible targets in (2.1) then have quadratic
field `K_0`.

After removing finitely many parameters, these are full rank-five
`K_0 times L_3` field fibres of the same fixed map `Phi`, with no
rational source point. If `M_S(B)` counts their distinct target images
of actual projective height at most `B`, then

\[
 \boxed{M_S(B)=c_S\log B+O_S(\sqrt{\log B}),
             \qquad c_S>0.}                        \tag{2.2}
\]

In particular they include infinitely many distinct cubic-field and
complete fibre-algebra isomorphism classes. No quantitative
distinct-class asymptotic is asserted.

**Proof.** A real cubic always has a real root, so the real place
imposes no further condition. Fix a finite prime `p` in `S`. At `P_0`,
one of the two separable factors `q,h` has a root in `Q_p`. Its
derivative at that root is nonzero because the fibre is etale. The
coefficients are continuous functions of `v`. Consequently this root
persists for `v` in a sufficiently small `p`-adic neighbourhood of
`v(P_0)`.

This last fact follows directly from the strict Hensel criterion.
After a fixed scaling of the variable and coefficients to an integral
model, evaluate at the original root. The polynomial value tends to
zero as `v` tends to `v(P_0)`, while the derivative remains nonzero
with fixed valuation. On a small enough neighbourhood, the strict
inequality `v_p(H(a))>2v_p(H'(a))` holds and gives a nearby local
root. The exact classical input is Milne, *Algebraic Number Theory*,
version 3.00, Theorem 7.32, p. 114:
[primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf).

The preimage in `E(Q_p)` is thus an open neighbourhood of `P_0`.
Choose a sufficiently small open subgroup `V_p` of `E(Q_p)` whose
translate `P_0+V_p` lies in this neighbourhood. The elementary local
group facts needed here are as follows. The projective curve `E(Q_p)`
is compact. At its identity, a local analytic parameter identifies a
neighbourhood with a `p`-adic ball; addition and inversion have
convergent series with linear terms `t+u` and `-t`. On a sufficiently
small ball, the higher terms preserve the ball, so this ball is an
open subgroup. Shrinking the ball gives arbitrarily small such
subgroups. Compactness then implies that each open subgroup has
finite index, since its quotient is compact and discrete.

The pullback `E(Q) intersect V_p` consequently has finite index in
`E(Q)`. Intersect these subgroups over the finitely many primes in
`S`, and further intersect with the kernel of the finite squareclass
map. This gives a finite-index subgroup `Gamma_S`. Every point of
`P_0+Gamma_S` has the desired local roots and the same quadratic
squareclass as `P_0`. The negatives have exactly the same parameters
and targets, so (2.1) has both properties as well.

Remove the finite inadmissible set and the finite cubic-root exceptions
from the genus-four/Faltings argument in the arithmetic note. The
remaining fibres are exactly products of irreducible quadratic and
cubic fields with matching resolvents. Neither component is `Q`, so
there is no rational source point.

Here is the target count, including negation. Write
`I_S=[E(Q):Gamma_S]`, and let `r_S` be the number of distinct cosets
in (2.1), which is one or two. The canonical-height lattice argument
of EEC1 counts every finite-index coset with `1/I_S` of the full
point main term and the same `O_S(sqrt(log B))` boundary error.
This remains true with torsion: the index is the product of the free
lattice index and the torsion quotient index. The union (2.1) is
symmetric under negation. Away from finitely many exceptions, each
target has exactly the two preimages `P,-P`, because the parameter
curve is birational to its target curve. Dividing by these two
preimages gives

\[
 c_S=\frac{r_S}{I_S}C_E>0.                          \tag{2.3}
\]

Finite deleted sets do not affect (2.2). Finally, the arithmetic note's
Corollary 3.3 shows that any one fixed cubic-field class yields only
finitely many targets, by applying Faltings over that number field to
the genus-four incidence curve. An infinite subfamily here must
therefore contain infinitely many cubic-field and fibre-algebra
classes. \(\square\)

This proof uses the existing canonical-height, rank, torsion and
genus-four inputs. It does not require a new parameter family or a
claim that the displayed rational points form a saturated basis.

**Corollary ELS1a (finite-local target density).** Fix any finite set
`S` of rational primes, and let `N_{S,K_0}(B)` count admissible targets
on this elliptic slice of actual projective height at most `B`, with
quadratic field `K_0=Q(sqrt(-81275190))`, no rational source point,
and a source point over every `Q_p`, `p in S`. Then

\[
 \boxed{N_{S,K_0}(B)=\Theta_S(\log B).}             \tag{2.4}
\]

The subfamily may be chosen as the target image of
`P_0+[M]E(Q)` for one suitably divisible even integer `M` depending
on `S`. It retains rank two and the same quadratic field.

**Proof.** The finite quotient `E(Q)/Gamma_S` in ELS1 is annihilated
by some integer. Choose an even multiple `M` of that integer. Then
`[M]E(Q) subset Gamma_S`, so the coset lies in every required local
neighbourhood of `P_0`. Evenness also directly makes its squareclass
equal to that of `P_0`, since the two-isogeny squareclass homomorphism
has exponent two.

For the target count use the symmetric union of this coset and its
negative, which has the same target image as the single coset.
Since `E(Q)` has rank two and torsion of order two,
`[E(Q):[M]E(Q)]=2M^2`. Let `r_M` be one or two according as those
two cosets coincide or differ. The lattice and negation count of
ELS1 gives this subfamily, after finite exclusions, the asymptotic

\[
 \frac{r_M}{2M^2}C_E\log B+O_S(\sqrt{\log B}).      \tag{2.5}
\]

The coefficient is positive. Removing the finitely many inadmissible
or reducible-cubic exceptions leaves full field fibres with no
rational point. This proves the lower bound in (2.4); EEC1's full
matching-target count gives its upper bound. No lower coefficient
uniform in `S` or in growing `M` is asserted. \(\square\)

For the whole finite-place set, [EFA1--EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md)
strengthen the order statement (2.4). With `p_0=2709173`, for every
`a in {1,2,5p_0,10p_0}` the full field targets of quadratic class
`a` passing a fixed finite set `S` satisfy

\[
 N_{S,a}(B)=\mu_{S,a}C_E\log B+o_{E,S}(\log B).
\]

The constant is normalized Haar measure on the actual joint closure,
including the finite squareclass coordinate, and satisfies
`0<=mu_{S,a}<=1/4`. For the certified class `a_0=10p_0` it is
positive for every fixed finite `S`. At `S={2,3,17}` its exact value
is `1/8`, by [J2174](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md).

## 3. The exact remaining all-primes criterion on this slice

**Corollary ELS2 (a split-prime discriminant sieve).** For every full
field fibre on this elliptic slice, let `L_3` be its cubic field and
write its maximal-order field discriminant as

\[
 \operatorname{Disc}(L_3)=d_K f^2.
\]

Then the fibre is Hasse-failing if and only if

\[
 \boxed{f\text{ is squarefree and every prime dividing }f
                          \text{ is }1\pmod3.}     \tag{3.1}
\]

The value `f=1` is allowed. Polynomial discriminants cannot replace
the maximal-order field discriminant in this test.

**Proof.** The four possible quadratic fields have radicands
`-3,-6,-15p_0,-30p_0`, where `p_0=2709173` is prime and is not
three. Thus each field ramifies at three, with `v_3(d_K)=1`.
None splits over `Q_3`. Theorem RDF1 of
[the field-discriminant classification](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md)
therefore forbids `3` in `f`. Its remaining primes have exponent one
and must be `1 mod3`. Conversely exactly those conditions satisfy
RDF1. \(\square\)

## 4. The next arithmetic problem

The focused unresolved question is whether infinitely many points of
this **existing** elliptic slice pass (3.1), equivalently whether its
everywhere locally soluble, globally insoluble target count `A_E(B)`
tends to infinity. The certified point (1.1), EEC1 and the joint local sieve give

\[
 1\le A_E(B)\le\frac{1687}{6936}C_E\log B+O_E(\sqrt{\log B})
                \qquad\text{for sufficiently large }B.         \tag{4.1}
\]

The [joint local sieve](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md)
gives the exact local-only coefficient `5321/23946` at `2,3,17`.
The [sharp error theorem](ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md),
ESQ3--ESQ4, upgrades this to the Hasse upper bound
`A_E(B)<=(5321/23946) C_E log B+O_E(sqrt(log B))`.
The [five-adic refinement](ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md)
improves the coefficient to `2300089/10392564` at the same error
order. It multiplies the class-one coefficient by `16/31` and
leaves the other three coefficients unchanged.
The certified quadratic field has upper coefficient `1/8` at the
same error order. The [explicit local criterion](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md)
replaces the varying tame tests by `v_p(128x-25865)=0 mod3` at every
`p=2 mod3`, `p>=11`, together with the three local tests at `2,3,5`.
Fixing a cubeclass of that rational function produces a genus-four
curve and only finitely many rational points, so a construction with
finitely many such cubeclasses cannot establish Hasse infinitude.
The tests at `2,3,5` are now completely explicit in the
[three-adic classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md)
and the [two/five-adic classification](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md).

[EQA2--EQA3](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md)
prove rationality of every finite-place coefficient, with an elementary
quantitative error involving the nonautomatic-prime count `r_a`.
ESQ3 sharpens the error to `O_(E,S)(sqrt(log B))` using the
Hirata-Kohno--Takada logarithm bound. The
[eligibility sieve](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md) requires
the class `a` locally square for a denominator obstruction; a
numerator obstruction additionally requires `51730` locally square.
For each nontrivial class the eligible prime sets have respective
densities `1/4` and `1/8`. These are prime eligibility densities.

Theorem ELS1 supplies finite-index subfamilies, while EFA1--EFA2 count
the whole set passing any fixed finite collection. The additional
all-primes step is isolated in
[EFA3](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md). For the certified
quadratic class, the positive finite-place measures decrease to a limit
`mu_infinity` between zero and `1/8`. This limit bounds the normalized
all-places count from above. If it vanishes, that fixed-class count is
`o(log B)`, which still leaves its infinitude open. If it is positive,
the leading asymptotic requires the exact tail condition in EFA3:
targets passing the small-prime tests but failing at a larger prime
must have negligible logarithmic count as the cutoff grows. The limit's
value and this tail control remain undetermined. The same distinction
holds after summing over all four classes. No interchange of the
prime-cutoff and target-height limits is asserted. EFA4 identifies
positivity of the limiting measure with summability of actual
conditional prime losses. The
[ambient Haar theorem](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md)
shows why independent ambient local fractions cannot determine this
intrinsic measure: the rational closure itself is ambient-null.

The all-depth conditional kernel and weighted pathwise theorem in
[ETP4--ETP5](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md) avoid a
separate depth-one distribution hypothesis. An almost-sure divergent
sum of their predictable hazards would give a zero limiting measure;
this arithmetic divergence has not been established. The
[one-prime positivity theorem](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md)
shows that every individual local test is feasible in every class.
Its cofinite-extension result leaves possible simultaneous finite-set
obstructions in the two remaining classes confined to an explicit
finite exceptional set relative to any supplied feasible coset.

All polynomial, field, rank and local-seed identities are already
covered by the retained checkers linked above. The new finite-place
sieve is the written openness and finite-index lattice argument. No
additional finite experiment is promoted to an all-primes conclusion,
and no formal verification or external review is claimed.
