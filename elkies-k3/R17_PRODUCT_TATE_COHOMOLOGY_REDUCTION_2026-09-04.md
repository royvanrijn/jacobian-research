# The 2-primary product-character quotient: exact reduction and input gate

<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-ZERO-TATE-CLASS-EXCLUSION 9e1c09d47fcf0bde -->
<!-- status-consumer: EC-K3-R17-PRODUCT-SURVIVOR-GALOIS-HEIGHT-GATE 32b87c35a2573768 -->

## Status

This note closes the formal reduction of the product-character question to
finite integral linear algebra and a two-descent/Kummer calculation.  The
complete two-prime Frobenius classification now closes twelve targets more
strongly: each product twist has geometric Mordell--Weil rank zero, so
`A-=0`, `Gamma_d=0`, and `H_d=0`.  For the five persistent Tate-factor
survivors, the full Mordell--Weil lattice after quadratic base change, or an
equivalent complete two-Selmer calculation, is not present in the current
certificates.  The exact 49-class norm-twelve trace-parity
certificate and the no-hit inversion of all `49 * 17 = 833` residual cases do
exclude the zero Tate class for a height-eight section under the stated direct
polynomial and local-component height hypotheses.  Their five quotients
`H_d`, the existence of any height-eight section, and all possible nonzero
classes remain `UNKNOWN`.

The targets are exactly the seventeen rank-one rows selected by
[`elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json)
from
[`elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json).
Their pair keys, in the stored order, are

```text
alternate-orbit-1463f:alternate-orbit-19bad
alternate-orbit-19bad:alternate-orbit-083ad
alternate-orbit-11ae6:alternate-orbit-0f82c
alternate-orbit-0f82c:alternate-orbit-025be
alternate-orbit-025be:alternate-orbit-13dbe
alternate-orbit-1ad20:alternate-orbit-1b24d
alternate-orbit-1a465:alternate-orbit-19b4e
alternate-orbit-19ead:alternate-orbit-146dc
alternate-orbit-19b4e:alternate-orbit-17b71
alternate-orbit-11ee2:alternate-orbit-0c36e
alternate-orbit-0c36e:alternate-orbit-02bf1
alternate-orbit-0fda0:alternate-orbit-1a6c8
alternate-orbit-1059f:alternate-orbit-1db8d
alternate-orbit-0fda0:alternate-orbit-1037d
alternate-orbit-0c10b:alternate-orbit-17a1a
alternate-orbit-1ede3:alternate-orbit-1c364
alternate-orbit-13dbe:alternate-orbit-1019b
```

The product quartics themselves are stored in the shortlist; duplicating
their large coefficients here would create a second authority surface.

## 1. The exact Tate-cohomology sequence

Fix one product squareclass `d`, and put

```text
K=QQ(u),       L=K(sqrt(d)),       G=<sigma>,
A=E(L),
A+ = ker(1-sigma)=E(K),
A- = ker(1+sigma).
```

The standard twist isomorphism identifies `A-` with `E^(d)(K)`.  All four
roots of `d` avoid the `24 I1` discriminant in the certified targets.  Hence
the base-changed elliptic surface has `48 I1` fibres, arithmetic genus four,
and no reducible fibres.  A nonzero torsion section would have canonical
height

```text
2*chi + 2(P.O) = 8 + 2(P.O) > 0,
```

which is impossible.  Thus `A` is torsion-free, and the Mordell--Weil group
with its height pairing is an honest free integral `G`-module.

Define the integral character-glue group

```text
Gamma_d = A/(A+ direct_sum A-).
```

Both eigensublattices are primitive and orthogonal.  The identity

```text
2P=(1+sigma)P+(1-sigma)P
```

shows that `Gamma_d` is an `F_2`-vector space.  There is a canonical exact
sequence

```text
0 -> Gamma_d --partial--> A-/2A-
   -> Hhat^(-1)(G,A) -> 0,                         (1)

partial([P]) = (1-sigma)P mod 2A-.                (2)
```

Indeed, (2) is unchanged after adding an invariant point and changes by
twice an anti-invariant point after adding an anti-invariant point.  If
`(1-sigma)P=2T` with `T in A-`, then `P-T` is invariant, so `[P]=0` in
`Gamma_d`; hence `partial` is injective.  Its image is exactly
`(1-sigma)A/2A-`.  Taking the cokernel proves

```text
H_d = Hhat^(-1)(G,E(L))
    = ker(1+sigma)/(1-sigma)E(L)
    = (A-/2A-)/partial(Gamma_d).                  (3)
```

In particular

```text
dim_F2(H_d) = rank(A-) - dim_F2(Gamma_d).          (4)
```

This is the precise relation between the obstruction and integral
character glue.  It uses the **full** base-change Mordell--Weil group.  The
known orthogonal sublattice `R17(2)` and a hypothetical anti-invariant line
do not determine `Gamma_d`; inserting only visible sublattices in (3) would
assume the saturation statement that is to be proved.

### Matrix presentation

Given an integral basis of `A` and the integral involution matrix `S`, let
`B_-` be a saturated column basis of `ker_Z(1+S)`.  There is a unique integral
matrix `D` with

```text
(1-S) = B_- D.
```

Then

```text
H_d = coker(D: Z^rank(A) -> Z^rank(A-)).           (5)
```

The Smith invariants in (5) are all `1` or `2`.  Equivalently,

```text
H_d = F_2^rank(A-) / column_span_F2(D mod 2).       (6)
```

Equations (5)--(6) are the smallest exact certificate format for a completed
integral-lattice computation.

[`compute_involution_tate_hminus1.sage`](scripts/compute_involution_tate_hminus1.sage)
implements (5)--(6) for a supplied full integral lattice, checks the Smith
invariants, and emits quotient functionals for candidate anti-invariant
vectors.  Its two synthetic controls are replayed by

```bash
sage -python \
  elkies-k3/scripts/compute_involution_tate_hminus1.sage --self-test
```

## 2. Kummer form of the same quotient

For `F=K,L`, the multiplication-by-two Kummer sequence

```text
0 -> E[2] -> E(Fbar) --[2]--> E(Fbar) -> 0
```

gives an injection

```text
kappa_F: E(F)/2E(F) -> H^1(F,E[2]).
```

Inside `H^1(L,E[2])`, put

```text
U+ = res_(L/K)(kappa_K(E(K)/2E(K))),
U- = kappa_L(A-/2A-).
```

The natural maps from both eigenspaces modulo two to `A/2A` are injective.
For example, if `T in A-` and `T=2P` in `A`, then
`2(1+sigma)P=(1+sigma)T=0`; torsion-freeness gives `(1+sigma)P=0`, so already
`P in A-`.  The invariant argument is identical with `1-sigma`.  For
`T in A-`, the following are equivalent:

```text
T=(1-sigma)P for some P in A;
there is tau in E(K) with T+tau in 2A;
kappa_L(T)=res_(L/K)(kappa_K(tau)).                 (7)
```

The forward implication uses `tau=(1+sigma)P`; the converse uses
`P=(T+tau)/2`.  Therefore

```text
partial(Gamma_d) = U- intersection U+,
H_d = U-/(U- intersection U+).                    (8)
```

Here (8) means the corresponding point-Kummer subspaces, not the whole
cohomology group.  Replacing a Mordell--Weil Kummer image by a Selmer group
produces an upper candidate space, not an equality unless the descent is
complete.

For the short twist model

```text
E^(d): Y^2 = X^3 + d^2 A X + d^3 B,
```

let `theta_d` be the image of `X` in the cubic etale algebra

```text
R_d=K[theta_d]/(theta_d^3+d^2*A*theta_d+d^3*B).
```

Away from the usual two-torsion exceptional cases, the explicit
two-descent representative is

```text
kappa_K(X,Y) = X-theta_d in R_d^*/R_d^(*)2,
Norm(X-theta_d)=Y^2.                               (9)
```

Under the twist isomorphism over `L`, division by `d` changes (9) by a square
because `d` is a square in `L`.  Thus (9) can be compared directly with the
seventeen known invariant Kummer generators `x(P_i)-theta` after restriction
to `L`.

## 3. Exactly which height-eight classes remain

Let

```text
C_d(8) = {[T] in H_d : T in A-, height_(E^(d))(T)=8}.
```

The zero-class carrier calculation has two exact layers.  The completed
norm-eight bisection inversion rejects all `63,917 * 17` minimal-trace target
comparisons.  The exact trace-parity certificate proves that the only other
possible trace minima are the 49 norm-twelve parities.  The deep trace
inversion then rejects all `49 * 17 = 833` residual comparisons.  Consequently
zero is not in `C_d(8)` for any of the seventeen targets, under the direct
polynomial and local-component height hypotheses used by those certificates.
These results are recorded in
[`elkies-k3-r17-norm12-11952-product-tate-parity-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json)
and
[`elkies-k3-r17-norm12-11952-product-deep-trace-inversion-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-product-deep-trace-inversion-v1.json).

Every height-eight representative is primitive.  Its image in `A` has
height sixteen; if it were `nQ` for `|n|>=2`, then `Q` would have height at
most four on the base-changed surface, contradicting the lower bound eight
from the no-reducible-fibre height formula.  Hence a surviving minimal
section gives a genuine nonzero primitive class in (3), not the multiple of
a shorter section.

Once a full anti-invariant Gram matrix and the glue map are known,
`C_d(8)` is finite and is enumerated exactly as follows:

1. enumerate every integral vector of twist height eight (equivalently cover
   height sixteen) in `A-`;
2. map its coordinate vector through the quotient (6);
3. discard the zero class and identify duplicates (sign is already invisible
   over `F_2`).

If in addition `A-` is the primitive rank-one lattice `<16>` generated by a
height-eight twist section, the exclusion of the zero class gives the
conditional corollary

```text
H_d=F_2, and the height-eight generator is its unique nonzero class.
```

In higher rank, (6) is still complete.  The geometric Shioda--Tate bound on
the `chi=4` twist surface gives `rank(A-)<=22` (four `I0*` fibres contribute
root rank sixteen), so even the unsliced ambient quotient has at most
`2^22` classes.  A Selmer calculation should make the actual list much
smaller.

## 4. Fail-closed computation for the seventeen targets

For each stored product quartic `d`, the exact route is:

1. Build the finite `S`-squareclass space in `R_d^*/R_d^(*)2`, with `S`
   containing the places over `2*d*Delta` and infinity, impose square norm,
   and impose the local Kummer image at every place.  This gives the complete
   two-Selmer group only after the global `S`-unit/class-group and all local
   image computations are certified.
2. Compute the span `U+` from the saturated seventeen-section basis of
   `E(K)`.  Restrict these classes to `L` and quotient the anti-invariant
   candidate space by its intersection with `U+`, as in (8).
3. Report separately:

   ```text
   MW Kummer image / invariant intersection     exact H_d;
   Selmer candidate / invariant intersection    upper candidate space only.
   ```

   A zero Selmer quotient proves `H_d=0`.  A nonzero Selmer quotient does not
   prove that a Mordell--Weil class exists.
4. Apply the height-eight local conditions to the nonzero representatives.
   Component data at the four `I0*` fibres and the degree bound at infinity
   should be applied before any polynomial solve.

This calculation is independent for the seventeen `d`; no rank-one theorem
for the genus-one **base Jacobian** computes the Mordell--Weil rank of the
product twist.

## 5. Class-sliced direct section equations

A disjoint section of the product twist lies in the direct box

```text
deg(X)<=8,       deg(Y)<=12,
Y^2=X^3+d^2*A*X+d^3*B.                            (10)
```

For a residual Kummer representative `alpha_c in R_d^*/R_d^(*)2`, membership
in class `c` is the exact parity constraint

```text
X-theta_d = alpha_c * Z^2 in R_d,                 (11)
```

for `Z=z0+z1*theta_d+z2*theta_d^2`.  Choose the finitely many denominator
profiles of `Z` allowed by the divisor of `alpha_c` and the height-eight
bounds, clear only the declared denominators, and equate the three
coefficients in (11).  Equations (10)--(11) are the direct polynomial-section
system sliced by one global mod-two class.  Merely checking valuation parities
at a subset of places is only a necessary filter unless those places are
proved to be a complete squareclass-coordinate basis.

At a good odd prime, reduce `d`, the cubic algebra, every `alpha_c`, and the
known invariant Kummer span.  Use the prime only if

```text
d*Delta and the cubic discriminant have good reduction,
the selected residual Kummer classes remain distinct modulo the known span,
all denominator profiles and the infinity chart remain represented.
```

Then solve only the class-sliced systems (10)--(11).  To turn an empty
finite-field result into a characteristic-zero nonexistence theorem, the
export must include a saturated projective compactification and all boundary
charts.  An empty affine coefficient chart alone has the same boundary noted
in the earlier `msolve` experiment and is not an upper bound.  Any modular
survivor must be lifted over `QQ` and checked by literal substitution,
Kummer-class equality, height, and non-torsion.

## 6. Independent finite-field upper-bound route

For a good reduction of the twist surface, an exact Frobenius polynomial on
`H^2` gives an unconditional upper bound on the geometric Picard number by
counting eigenvalues of the form `p` times a root of unity.  Shioda--Tate then
gives

```text
rank E^(d)(QQbar(u))
 <= rho(surface over Fpbar) - 2 - 4*rank(D4)
 =  rho(surface over Fpbar) - 18.                  (12)
```

Thus one prime with upper bound `rho<=18` proves geometric product-twist rank
zero and closes the target without computing (3).  Fibral `n=1` trace
averages do not determine this bound: the full, exactly reconstructed
Frobenius polynomial (or an equivalent certified cohomological upper bound)
is required.

The current finite-field audit
[`audit_r17_norm12_11952_product_twist_finite_field_bounds.sage`](scripts/audit_r17_norm12_11952_product_twist_finite_field_bounds.sage),
with certificate
[`elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json),
proves good reduction and computes the exact Frobenius power sums only for
`n=1,2` at `p=131,137`.  The nonconstant elliptic `L`-polynomial has degree
28, so even with its functional equation exact reconstruction requires power
sums through `n=14`.  The audit therefore gives no finite-field
Mordell--Weil upper bound.

Before using target-specific component Galois, the general geometric height
inequality on a `chi=4` surface with four `I0*` fibres gives

```text
8 = 8 + 2(P.O) - sum_v contr_v(P),
sum_v contr_v(P) <= 4,
```

hence only `P.O<=2`.  For the five surviving product twists, the subsequent
exact residue-field calculation is sharper: every two-division cubic is
irreducible at the quadratic branch place, so no rational section can meet a
nonidentity `I0*` component.  All corrections vanish and
`height(P)=8+2(P.O)`.  Thus a rational height-eight section has `P.O=0`; one
trivial smooth-fibre specialization would exclude the whole height-eight box,
while it would say nothing about height at least ten.

## Current hard input boundary

The repository currently supplies the full invariant lattice `E(K)=R17`, the
seventeen product quartics, the complete norm-eight inversion, the exact
49-class norm-twelve trace complement, and the no-hit deep inversion of all
`49 * 17 = 833` residual trace/target pairs.  Thus the zero Tate class is
excluded for a height-eight section under the stated direct polynomial and
local-component height hypotheses.  It does not supply any of the following
target-specific inputs:

```text
a full integral basis of E(L) with sigma;
a complete two-Selmer group for E^(d)/K;
a full H^2 Frobenius polynomial giving rho<=18.
```

The subsequent full toric calculations supply complete degree-28 polynomials
at `p=131,137` for these five targets, but their Tate degree two gives only
the geometric bound `rank<=2`.  The component-Galois gate makes the
height-eight problem exactly `P.O=0`; it does not compute the nonzero Kummer
classes in that box or address height at least ten.

Therefore the quotient `H_d`, its nonzero classes, and the existence of a
height-eight section remain `UNKNOWN`.  Equations (3), (8), and (12) are exact
closing criteria for those stronger questions, but none may yet be marked as
passed.

## Subsequent one-target closure

The finite-field boundary above has now been crossed for one of the seventeen
targets.  For shortlist rank `55`, character pair

```text
alternate-orbit-0fda0 : alternate-orbit-1037d,
```

the complete degree-28 Frobenius polynomial at `p=131` has no Tate root.
The independent verifier reconstructs the full degree-46 `H^2` polynomial,
the exact degree-18 `U+4D4` Frobenius factor, and a Weil degree-28 quotient.
Therefore `rho<=18`, and (12) proves geometric product-twist rank zero.  See
[`R17_ALTERNATE_Q80_PRODUCT_TWIST_RANK_ZERO_2026-09-04.md`](R17_ALTERNATE_Q80_PRODUCT_TWIST_RANK_ZERO_2026-09-04.md).

<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-RANK55-GEOMETRIC-RANK-ZERO 61d8e223f0ef7d0f -->

For this target the twist identification gives `A-=0`.  Since the pulled-back
`48I1` Mordell--Weil group `A` is torsion-free, (1) has zero middle term and
forces both `Gamma_d=0` and
`Hhat^(-1)(<sigma>,A)=0`.  Thus neither the zero nor a nonzero Tate class can
carry a product-character section.

Accordingly, the preceding `UNKNOWN` conclusion no longer applies to the
rank-55 product twist or its Tate-cohomology quotient.  The complete campaign
below reduces the remaining input boundary further, to five targets.
<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-RANK55-GEOMETRIC-RANK-ZERO 61d8e223f0ef7d0f -->

## Complete seventeen-target classification

The same hardened computation has now been run at `p=131` for all seventeen
targets and at `p=137` for every survivor.  Twelve targets have Tate degree
zero and hence satisfy `A-=Gamma_d=H_d=0`.  Exactly five retain Tate degree
two at both primes and have only the geometric upper bound `rank(A-)<=2`:

```text
alternate-orbit-19bad : alternate-orbit-083ad
alternate-orbit-11ae6 : alternate-orbit-0f82c
alternate-orbit-0f82c : alternate-orbit-025be
alternate-orbit-11ee2 : alternate-orbit-0c36e
alternate-orbit-0c10b : alternate-orbit-17a1a.
```

These five are the complete descent/class-sliced worklist.  A Tate factor at
both primes is not a positive-rank theorem.  The full table and replay are in
[`R17_ALTERNATE_Q80_ALL17_PRODUCT_TWIST_CLASSIFICATION_2026-09-04.md`](R17_ALTERNATE_Q80_ALL17_PRODUCT_TWIST_CLASSIFICATION_2026-09-04.md).

<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-ALL17-TORIC-CLASSIFICATION cfb2417a30fab18d -->

## Five-survivor component-Galois gate

At all ten quadratic branch places of the five survivors, the smooth fibre's
two-division cubic is irreducible over the residue field.  Since the three
nonidentity components of the twisted `I0*` fibre are indexed by the three
nonzero two-torsion points, none is residue-Galois fixed.  A `QQ(u)`-rational
section therefore has zero local correction at every branch fibre and

```text
height(P)=8+2(P.O).
```

Consequently every height-eight section is in the direct `P.O=0` box and,
by the zero-class exclusion above, would represent a nonzero Tate class.  The
constructor-aware order starts with `19bad:083ad`, whose primitive base
generator maps to a 124-bit `u`-value.  Its complete independently audited
`p=151` quotient again has Tate factor `(Z-1)(Z+1)`, so it persists as a
candidate but no section is proved.  At `p=137` all five have the same
normalized factor.  A rational section specializes into its unique
Frobenius-fixed direction, proving arithmetic product rank at most one for
each survivor.  Hence one nonzero rational section would prove that rank is
exactly one.  See
[`R17_PRODUCT_SURVIVOR_GALOIS_HEIGHT_GATE_2026-09-04.md`](R17_PRODUCT_SURVIVOR_GALOIS_HEIGHT_GATE_2026-09-04.md).

<!-- status-consumer: EC-K3-R17-PRODUCT-SURVIVOR-GALOIS-HEIGHT-GATE 32b87c35a2573768 -->
