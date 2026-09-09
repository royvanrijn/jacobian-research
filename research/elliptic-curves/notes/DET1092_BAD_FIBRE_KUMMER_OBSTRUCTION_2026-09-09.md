# Singular-fibre ramification cannot carry the primitive302 seed class

## Theorem and scope

**New deduction, with verified equation-only applications.** On the fixed
determinant1092 parent write

\[
 E/\mathbf Q(t):Y^2=f_t(X),\quad
 L=\mathbf Q(t)[\theta]/f_t(\theta),\quad \Delta(t)=\operatorname{disc}_X(f_t).
\]

Let `Gamma` be the smooth projective cubic spectral curve. Let `A_bad`
consist of norm-square classes in `L*/L*2` with even valuation over **every
smooth parameter fibre**, including infinity. Ramification over singular
fibres is allowed. With `H=M17`, put

\[
 \boxed{\beta=-\Delta f'_t(\theta).}
\]

Then

\[
 \boxed{A_{\rm bad}=\delta(H)\ \sqcup\ \beta\delta(H),
             \qquad N_{L/\mathbf Q(t)}(\beta)=\Delta^4.}
\tag{1}
\]

This is an all-degree classification of this specified construction family,
not a generic Selmer-dimension computation. On302 **and all eight frozen
controls**, every class in the second coset fails a local Kummer condition.
Consequently a lift of the certified primitive302 seed class must ramify
over at least one **smooth** parameter fibre. Singular-fibre ramification
alone is insufficient, however complicated its rational-function presentation.

This does not exclude independent points whose mod-2 class is inherited,
nor any specialized rank jump. The productive anti-traces remain explicit
examples of that distinction. Nor does it exclude covers obtained after a
base change or lifts ramified at smooth parameter places.

## 1. Audit of completed inputs

**Verified prior applications.** The
[unramified Kummer theorem](DET1092_UNRAMIFIED_KUMMER_OBSTRUCTION_2026-09-09.md)
proves that every everywhere-unramified norm-square class equals `delta(H)`.
It also proves that two norm-square classes with identical parity divisors
differ by a generic Kummer class.

The earlier
[polynomial-lift theorem](DET1092_SEED_KUMMER_COVER_2026-09-08.md#polynomial-lift-obstruction-an-entire-family-not-one-copied-abscissa)
already proves that `Delta` is irreducible of degree24. Its modular factor
patterns at149 and151 have common subset sums only0 and24. It also constructs
the double root `r` and simple root `s` over `Q[t]/Delta`, with

\[
 f_t(X)=(X-r)^2(X-s)\pmod\Delta,\qquad r\ne s.
\tag{2}
\]

These results are reused, not rediscovered by rational factorization.
The new independent checker replays just their needed polynomial and
finite-field certificates, not the degree24 square-root boundary calculation.

## 2. Why there is precisely one additional coset

**New deduction.** Since the discriminant has one rational orbit and is
squarefree, `Gamma` has exactly two closed points above that orbit: a
ramified point `R` of index2 and an unramified point `S` of index1. Both
residue fields equal `Q[t]/Delta` by (2).

For `alpha in A_bad`, the norm-square condition imposes

\[
 v_R(\alpha)+v_S(\alpha)\equiv0\pmod2.
\]

There are therefore only two possible residue patterns: `(0,0)` and
`(1,1)`. This uses residue-field degrees, which are both1 over the base
discriminant place; ramification indices must not be inserted in this norm
valuation formula. The zero pattern is exactly `delta(H)` by the completed
unramified theorem.

For a monic cubic, `N(f'(theta))=-Delta`. Also `N(-Delta)=-Delta^3`, so
the displayed `beta` has norm `Delta^4`. At `R`, the derivative has valuation1
and `Delta` has valuation2. At `S`, the derivative is a unit and `Delta`
has valuation1. Thus

\[
(v_R(\beta),v_S(\beta))=(3,1).
\]

There are no other finite zeros or poles. At infinity use the minimal
coordinate `X'=t^-4 X`. The discriminant has pole24, and the derivative
has pole8 with nonzero leading value on every branch, since the infinity
cubic is smooth. Therefore `beta` has pole32 at each geometric infinity
point, an even order. It realizes the unique nonzero allowed pattern.
Taking its quotient with any other such class proves (1).

No irreducibility claim about the infinity cubic is needed. The optional
checks modulo149 and151 both find a root and do **not** certify irreducibility.
Infinity is smooth and is included in the unramifiedness hypothesis.

## 3. An elementary odd-nodal-prime obstruction

**New deduction from local cubic algebra.** Let `p>3`, and suppose a
specialized monic cubic has integral coefficients, unit `c4`, and odd
positive discriminant valuation `n`. Then every local rational elliptic
Kummer image has even valuation in both factors of the cubic algebra,
whereas `beta` has odd valuation in both.

Here is a direct proof, without assuming a Tate parametrization. The
reduction has a double root and a distinct simple root. Hensel lifting
the simple root gives an integral factorization

\[
 f(X)=(X-s)((X-r)^2-d),\quad r-s\text{ a unit},\quad v_p(d)=n\text{ odd}.
\]

The cubic algebra is `Qp` times the ramified quadratic field
`F=Qp(sqrt(d))`. Use its normalized valuation, so `v_F(p)=2`.
For a local elliptic point with nonintegral abscissa `x`, the monic equation
makes `v_p(x)` even. Its Kummer valuations in both factors are even.
For integral `x`, consider the residue classes:

- Away from `r,s`, both Kummer factors are units.
- Near `s`, the quadratic factor is a unit; `y^2=f(x)` forces
  `v_p(x-s)` even and the quadratic Kummer factor is a unit.
- Near `r`, put `m=v_p(x-r)`, allowing `m=infinity`. Since `n` is odd,
  `v_p((x-r)^2-d)=min(2m,n)`. A square ordinate forces `2m<n`.
  The quadratic Kummer valuation is then `v_F(x-r-sqrt(d))=2m`, even,
  and the simple factor is a unit.

The identity has trivial class. At the rational2-torsion point `(s,0)`,
the extended Kummer value uses `f'(s)` in the simple factor; it is a unit,
as is the quadratic factor. Thus these exceptional point cases also pass.

On the simple factor, `v_p(beta)=n`. On the quadratic factor,
`f'(theta)=2(theta-r)(theta-s)` has valuation `n`, giving

\[
\boxed{v_F(\beta)=3n,\qquad v_p(\beta|_{X=s})=n.}
\tag{3}
\]

Both are odd. Hence `beta` is not a local Kummer image, and multiplying it
by **any** rational Kummer image cannot repair it. In particular this holds
for every generic section combination, not only a sampled finite subgroup.
Its2-covering is not locally soluble and its class is not in Selmer; it is
not a Sha class. The standard relation between covering solubility and
local Kummer images is described in
[Cremona--Fisher--O'Neil--Simon--Stoll, sections2--3](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf).

## 4. A factor-free equation-only test

**New deduction/application.** For a smooth rational specialization, let
`B` be the product of2,3, all cubic-coefficient denominators, the denominator
of `Delta`, and numerator and denominator of its nonzero `c4`. Starting
from `abs(numerator(Delta))`, repeatedly remove its gcd with `B` until
the remainder `R` is coprime to `B`.

If `R` is not a square, some prime outside this excluded support has odd
discriminant valuation. It meets all hypotheses of section3. This proves
the obstruction **without identifying or factoring that prime**.
The remainder, its floor square root, all gcd removals and a Bezout
coprimality witness are sufficient for independent checking.

**Verified application.** This uniform test passes on all nine unchanged
parameters. No exceptional point coordinates, later V3 point, catalogue
rank or search result is an input. The source roster predates this class.

| Frozen parameter label | Real obstruction | Odd-nodal-prime obstruction |
|---|---|---|
| scale-0131232 | Not supplied | Proved |
| scale-0257585 | Proved | Proved |
| scale-0487239 | Not supplied | Proved |
| scale-0177036 | Proved | Proved |
| scale-0043332 | Proved | Proved |
| scale-0590501 | Proved | Proved |
| scale-0290097 | Proved | Proved |
| scale-0748009 | Proved | Proved |
|302, at `t=0`| Proved | Proved |

The real test is separate: if `Delta(t)>0`, order the three real roots
`e1<e2<e3`. Every real elliptic Kummer image has positive first component,
including the extended value at `(e1,0)`. But
`-Delta*f'(e1)<0`, so `beta` and its generic multiples are excluded.
This real argument alone does not cover the two negative-discriminant
controls. The exact Sturm calculation gives18 real roots of the parent
discriminant, so we explicitly do **not** claim a constant real sign on
the whole parameter line.

The class was defined from the generic equation before this comparison.
The real test was evaluated first; after its two uncovered controls, the
odd-nodal-prime lemma was derived and applied uniformly to all nine.
This adaptive proof history is retained; the two tests are not presented
as a prospectively frozen joint predictor. No production policy changed.

## 5. Consequence for the seed and the other controls

**New obstruction.** At any of these nine parameters, a class in `A_bad`
can specialize to a rational Kummer image only through the inherited coset
in (1). Thus it cannot equal the primitive302 seed class, which the existing
certificate places outside that image.

The earlier degree26 opening for polynomial linear-norm lifts is narrowed
in a precise way. If `k(t)` is polynomial of degree greater than4, then
`f_t(k)(k-theta)` has infinity pole `4 deg(k)`, hence even. If it has no
ramification over finite smooth parameter fibres, it belongs to `A_bad`.
It therefore cannot specialize to the primitive302 seed class, at **any**
such polynomial degree. This is not a proof that all these polynomial
abscissas are generic, and does not exclude an independent doubled point
with inherited Kummer image.

The nine blinded M16 recoveries concern the smaller reference subgroup;
their full-M17 classes are inherited. Productive conic branches can remain
new, while their anti-traces have inherited mod-2 classes. The dependent
conic cycles in that same inherited image. Nothing in (1)--(3) changes
the completed halving-or-cycle decisions or predicts amplification.

The old calibrated302 continuation `f_t(k0)(k0-theta)` ramifies over12
smooth parameter values, so it lies **outside** this newly excluded class
of constructions. The theorem explains why smooth-fibre ramification there
is indispensable; it does not turn the copied abscissa into a generic selector.

## 6. Exact checkpoint and open endpoint

**Verified computation.** The
[packet](../../artifacts/generated-results/elliptic-curves/det1092_bad_fibre_kummer_v1/)
retains the construction, all nine local certificates, the two failed
Sage API preflights and the independent replay. Each process has a25-second
cap. The control constructor takes about half a second. No point search,
new parameter, integer factorization or group census was run.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_bad_fibre_kummer.sage
```

The proof is hybrid: explicit polynomial/valuation calculations and the
previous written cohomological theorem, not a formal proof assistant result.

**Open constructive endpoint.** A prospective selector still has to produce
a class with permitted ramification over smooth parameter fibres and prove
that a rational specialization is a primitive seed. Nonzero ramification,
norm-square status and distinct covering equations do not by themselves
prove local solubility, a rational point, or independence. The broad
singular-fibre-only route is now closed on the entire fixed panel; the
seed-incidence problem itself remains open.
