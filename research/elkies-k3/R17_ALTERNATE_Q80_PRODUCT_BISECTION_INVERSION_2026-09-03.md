# Alternate-Q80 product-character bisection inversion (2026-09-03)

<!-- status-consumer: EC-K3-R17-NORM12-11952-PRODUCT-BISECTION-INVERSION 6cfef74eb08601a6 -->

## Exact outcome

The eight-variable product-twist `msolve` pilot has been replaced by a
complete finite inversion through the norm-eight genus-one bisection layer of
the direct alternate-Q80 fibration.

The alternate Mordell--Weil lattice has exactly `63,917`
section-nonnegative parity classes of minimum norm eight.  For each class the
regular chord pencil was compared with each of the seventeen product
quartics belonging to exact rank-one `V4` bases.  There are no matches.

The exact result is

```text
norm-eight translation classes       63,917
exact rank-one product targets            17
projective coefficient comparisons 1,086,589
modular survivors                         0
exact squareclass hits                     0
```

Prime `131` obstructs `63,915` classes.  Two traces have bad reduction for
this calculation at `131`; both are obstructed at `137`.  The search includes
the member at `lambda=infinity`, and a zero modular coefficient vector is
retained as a survivor rather than discarded.

This is an exact negative result for the complete
norm-eight/pole-order-zero **bisection image**.  The integral word is
essential: the result does not by itself prove that the product twists have
no minimal section, because the inverse of the bisection construction has a
possible 2-primary descent obstruction.

## 1. Precise dictionary

Put

```text
K=QQ(u),       L=K(sqrt(d)),       Gal(L/K)=<sigma>,
```

where `d` is a squarefree quartic, and identify the quadratic twist
`E^(d)_L` with `E_L`.  Write

```text
E(L)^- = {T in E(L) : sigma(T)=-T}.
```

### Proposition

Modulo translation by a generic section in `E(K)`, genus-one bisections with
normalization `L/K` are equivalent to pairs

```text
(T,[tau]),
T in E(L)^- / {+/-1},       [tau] in E(K)/2E(K),
T+tau in 2E(L).                                      (1)
```

The bisection is disjoint from the zero section and lies in the minimal
rootless layer exactly when representatives can be chosen with

```text
height_E^(d)(T)=8,       height_E(tau)=8.             (2)
```

In the marked alternate frame, the second condition in (2) is precisely the
minimum-norm-eight parity layer

```text
D_tau=(2,2,w),       (w,w)=8,       w modulo 2M.      (3)
```

Thus the unconditional finite dictionary is:

> Minimal product-character sections which are integral coboundaries are
> equivalent, modulo generic-section translation and sign, to members of the
> finite set of disjoint genus-one bisection pencils whose branch squareclass
> is the product character.

Removing “which are integral coboundaries” requires a separate proof that
the relevant 2-primary Tate-cohomology quotient vanishes.  Rational
surjectivity after tensoring with `QQ` is not an integral saturation theorem.

### Proof

Let `C` be a bisection and choose one of its two points `P in E(L)` over the
generic old fibre.  Its conjugate is `sigma(P)`.  Set

```text
tau=P+sigma(P) in E(K),       T=P-sigma(P) in E(L)^-.
```

Translation by `S in E(K)` changes `P` to `P+S`, hence changes `tau` to
`tau+2S` and leaves `T` unchanged.  Interchanging the two sheets changes
`T` to `-T`.  This gives (1) in the forward direction.

Conversely, (1) gives

```text
P=(T+tau)/2 in E(L),       sigma(P)=(-T+tau)/2.
```

The unordered pair descends to a degree-two curve over `K`; changing `tau`
by `2S` translates that curve by `S`.  If two choices give the same `T`,
their two halves differ by a `sigma`-fixed point, so their traces differ by
twice a point of `E(K)`.  Hence the translation class is exactly `[tau]`.

On the degree-two base change, the invariant and anti-invariant height
subspaces are orthogonal and both source heights double.  Under (2),

```text
height(P)=(16+16)/4=8.
```

The base-changed rootless surface has `chi=4`, so Shioda's formula gives
`height(P)=8+2(P.O)`.  Therefore `P.O=0`.  Conversely, for a disjoint smooth
genus-one bisection the two lifted sections have height eight and meet at the
four simple ramification points.  Their sum and difference have height
sixteen on the base change, hence height eight on their invariant and twist
models.  This proves (2).

Finally, for `NS=U+(-M)`, a degree-two isotropic class disjoint from `O` is
`D_w=(2,2,w)` with `(w,w)=8`.  Translation is `w -> w+2x`, and

```text
D_w.S_x = (w-2x,w-2x)/4 - 2.
```

It is nonnegative on every section exactly when its parity coset has no
representative of norm below eight.  This proves (3) and the finite
dictionary. QED.

## 2. Complete norm-eight layer

The alternate frame is enumerated in its deterministic LLL row basis and
transported through

```text
short basis
  -> historical alternate frame
  -> direct compiled alternate frame
  -> saturated equation section basis.
```

The exact shell and parity counts are:

| norm | signed vectors | parity cosets hit |
|---:|---:|---:|
| 2 | 0 | 0 |
| 4 | 2,626 | 1,313 |
| 6 | 53,290 | 26,645 |
| 8 | 460,360 | 63,917 |

After removing the lower-norm cosets, all `63,917` surviving classes have
minimum norm exactly eight.  The cheapest representative in the saturated
equation basis has coefficients of absolute value at most five; its declared
score uses group additions, support, maximum coefficient, `L1` norm, and the
lexicographic word.

The certificate is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json),
SHA-256
`3525a09e52398242a935ad8fdbd36c6911c75863de51411783b731ab9c581aa1`.
Its complete table has SHA-256
`85f19177d9da3eac695b432cb5440653f70ee72d20d19c71bc77fc14ba8cdb69`.

## 3. Chord inversion

For a norm-eight trace

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=2,
```

the regular chord representatives are

```text
M=M0+lambda*h^2,       M0*Nx+Ny == 0 mod h^2.
```

Write `q(lambda)=sum_(j=0)^4 q_j lambda^j`.  Expanding the universal chord
discriminant gives

```text
q_0 = (M0^4-6M0^2 Nx-8M0 Ny-3Nx^2-4A h^4)/h^6,
q_1 = (4M0^3 h^2-12M0 h^2 Nx-8h^2 Ny)/h^6,
q_2 = (6M0^2 h^4-6h^4 Nx)/h^6,
q_3 = 4M0,
q_4 = h^2.
```

Every `q_j` has degree at most four in `u`.  Equality with a target quartic,
up to scalar, is therefore projective equality of two five-coefficient
vectors in the single parameter `lambda`.  The final scalar must be a square
in `QQ` for squareclass equality.

The modular gate homogenizes in `[lambda:mu]` and tests every point of
`P1(F_p)`.  Therefore a rational solution would reduce to a tested point at
every prime where the trace calculation is regular.  A class is rejected
only when one good prime has no projective match.  If a class survives the
declared primes, the checker forms the four exact coefficient minors over
`QQ[lambda]`, takes their gcd, tests every rational root and infinity, and
then tests the rational scalar squareclass.

## 4. Synthetic positive control

Before reading the seventeen targets, the checker constructs the member

```text
orbit 0x0c00d,       priority rank 1,       lambda=0.
```

Its chord quartic is squarefree, irreducible as a cover branch, coprime to the
alternate `24 I1` discriminant, and its lifted bisection satisfies the two
coefficient identities exactly.  The inversion path recovers this declared
member at `lambda=0 mod 131`.

This controls the lattice transport, trace group law, regular chord,
old-base chart, projective quartic normalization, and parameter comparison
used by the negative search.  It is not a target selected from a failed
search.

## 5. Exact negative result and boundary

For the seventeen exact rank-one pairs, no target survives:

```text
p=131 first obstructions    63,915
p=137 first obstructions         2
unresolved classes               0
squareclass hits                 0
```

The main certificate is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-product-bisection-inversion-v1.json),
SHA-256
`49809d43c1347339dcf4fc6e3aa07bda0e5cd02a68813dd3a4a1d5fde5d836c3`.
The complete per-class ledger has SHA-256
`582c0decbf0c8a98bb8fe68a595f2b572ce55cbecb07ed73dbc8606a5256e365`.

Consequently none of the selected product characters is realized by a
disjoint member of any section-nonnegative norm-eight genus-one bisection
pencil.  Equivalently, none supplies a minimal product-character section in
the integral coboundary image with norm-eight trace.

The following stronger sentence is **not yet proved**:

> None of the selected product twists has a minimal product-character
> section.

To prove it from this calculation one must additionally show that every
height-eight product-twist section lies in
`(1-sigma)E(QQ(u)(sqrt(q_iq_j)))`.  The cokernel is killed by two but need not
vanish integrally.  A non-coboundary minimal section would be invisible to
every bisection inversion, even though twice that section is a coboundary.
Likewise, the present negative does not imply that every remaining section
has positive zero intersection or comes from a higher-arithmetic-genus
curve.

Thus the next exact gate is a 2-primary saturation/descent calculation, not a
return to the twelve-equation product-twist `msolve` systems.

## Replay

```bash
sage -python \
  elkies-k3/scripts/rank_r17_norm12_11952_alternate_norm8_pencils.sage

sage -python \
  elkies-k3/scripts/rank_r17_norm12_11952_alternate_norm8_pencils.sage \
  --check

sage -python \
  elkies-k3/scripts/search_r17_norm12_11952_product_bisection_inversion.sage

sage -python \
  elkies-k3/scripts/search_r17_norm12_11952_product_bisection_inversion.sage \
  --check
```

The historical `msolve` exports and timeout manifests are retained as
reproducibility evidence, but they are no longer the active product-character
search route.
