# H3 q=8 module-intersection ledger — 2026-08-22

## Purpose and status

> **Historical ledger (superseded 2026-08-23).** The q8 equation is now
> certified over `QQ` and gives the `D13/MW4` child.  The canonical primitive
> child section is `S=Pmap+Qmap`, of height `24` and collision degree `10`.
> The degree-46 child route discussed below used `2*S` (and later an omitted
> `Dx` factor) and is retained only to explain discarded experiments.  See
> [`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md) for the active
> equation frontier.

This note consolidates the source-side H3 `q=8` equation work performed after
commit `8512558cdc1aac3d2dfe313c8cd61b005232dc10`.

The goal is the second exact H3 neighbor

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4.
```

The first `q6` equation is already exact.  The remaining `q8` task is to
compute the global two-dimensional Riemann--Roch space of the source-nef q8
class without replacing resolved line-bundle lattices by termwise endpoint
bounds.

This ledger is deliberately conservative.  It records successful exact or
two-prime calculations, but it also records the false starts and later
normalization corrections.  In particular, the phrases **reduced E7 lattice**
and **resolved E8 lattice** below refer to the q8 divisor with the global
fibre twist kept separate; they must not be confused with the ninth power of
a q6 representative carrying its own `-F` factor.

## 1. Exact source data already in the repository

The source-nef q8 class has old generic-fibre degree 18.  Relative to
`9*O+9*(-P1)` its exact vertical difference in the pinned
`[U,E7,E8,P1,P2]` frame is

```text
(-11,0,2,3,4,6,5,5,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0).
```

Thus the generic fibre is exactly

```text
9*O + 9*(-P1),
```

with generic Riemann--Roch basis

```text
B = 1,m,...,m^9, x,x*m,...,x*m^7,
m = (y-y(P1))/(x-x(P1)).
```

The function field is quadratic over `QQ(t)(m)`, with `x` satisfying a monic
quadratic, so these 18 functions are the complete generic-fibre space.

Relevant exact artifacts/scripts:

```text
scripts/derive_h92_q8_generic_rr_ambient.sage
scripts/derive_h92_q8_smooth_collision_frame.sage
scripts/derive_h92_q8_smooth_line_bundle_lattice.sage
scripts/derive_h92_q8_e7_local_target.sage
scripts/derive_h92_q8_actual_e7_gluing.sage
scripts/derive_h92_q8_actual_e7_power_module.sage
scripts/derive_h92_q8_e8_local_target.sage
scripts/derive_h92_q8_e8_complete_module.sage
```

The four smooth `O.(-P1)` collision fibres have collision polynomial `h` of
degree four.  In their actual regular algebra use

```text
p = y(P1)/x(P1),
q = (m-p)/h,
X = h^2*x.
```

The exact local line-bundle lattice there is

```text
<1,q,...,q^9, X,X*q,...,X*q^7>.
```

This is a **local h-adic frame**.  It is not a global basis over the old base.
That distinction became decisive below.

## 2. Child-side component-nef shortcut is closed as a dead end

Before returning to the literal source q8 divisor, the q6-child
component-nef route was pushed through several finite/infinity guesses.

The exact q-regular finite residue calculation found a nonzero degree-three
residue, but the bounded infinity sweep gave no branch-degree-four pencil.
The last safe fractional global-intersection experiment had

```text
ambient = 95
rank = 47
kernel = 48
```

at both `43` and `59`.  The 48-dimensional kernel was a base-twist ladder,
not `H0(D)`: it consisted of the expected sequence of shifted polynomial
coefficients associated with the missing degree-46 collision/gluing divisor.

Consequences:

- stop treating child finite and infinity modules as independent guessed
  fractional ideals;
- local regularity in `q_regular` can hide a nontrivial global base transition;
- the previous assumption "P0 is smooth at IV*, therefore P0 lies on the
  identity component" is invalid; the IV* component group is nontrivial;
- do not resume this child-side route until the exact resolved translation and
  Cartier transition are available.

This negative result motivated the source-side lattice-intersection approach.

## 3. E7-pole slack sweep: first smooth kernel at exactly degree(h)

The existing smooth-principal-part probe was run with `extra_h=0` and
`extra_e7=R`, for both good primes `43` and `59`.

The two primes agreed exactly:

| R | columns | smooth rank | smooth kernel |
|---:|---:|---:|---:|
| 0 | 54 | 54 | 0 |
| 1 | 72 | 72 | 0 |
| 2 | 90 | 90 | 0 |
| 3 | 108 | 108 | 0 |
| 4 | 126 | 124 | 2 |
| 5 | 144 | 140 | 4 |
| 6 | 162 | 156 | 6 |
| 7 | 180 | 172 | 8 |
| 8 | 198 | 188 | 10 |
| 9 | 216 | 203 | 13 |
| 10 | 234 | 218 | 16 |
| 11 | 252 | 229 | 23 |
| 12 | 270 | 240 | 30 |

The first kernel occurs at `R=4=deg(h)`, which correctly diagnosed that the
old termwise marked-E7 bound was too restrictive: a full h-factor can be
hidden by cancellation.

However, this did **not** produce the q8 pencil.

### R=4

The two-dimensional smooth kernel is killed by two exact singleton generic
E7 conditions at both primes:

```text
smooth kernel = 2
singleton constraints = 2
survivor = 0.
```

The full 2348-row generic residue calculation is unnecessary once these
necessary rows kill the space.

### R=5

The four-dimensional smooth kernel is killed by four singleton rows, again at
both primes.  All four obstructions are on `E7_1` and form a staircase in the
highest `u` coefficients of the `1` and `m` families:

```text
E7_1 ord -38 : m^0, i=19, k=2
E7_1 ord -37 : m^1, i=19, k=2
E7_1 ord -36 : m^0, i=18, k=2
E7_1 ord -35 : m^1, i=18, k=2.
```

Conclusion: increasing `extra_e7_pole` manufactures smooth-collision
cancellations by pushing poles toward E7_1, and E7_1 kills the new directions
one-for-one.  Do not continue this as a blind slack search.

## 4. Why a direct global q-frame looked promising, and why it fails

The exact h-adic frame suggested parameterizing globally in

```text
1,q,...,q^9, X,Xq,...,Xq^7.
```

Using the sharp E8 floor and the generic E7_1 order leaves only eleven pure-q
candidates:

```text
u^23*q^7,
u^25*q^8, u^26*q^8, u^27*q^8, u^28*q^8,
u^27*q^9, u^28*q^9, u^29*q^9, u^30*q^9, u^31*q^9, u^32*q^9.
```

All six actual resolved E7 node tests were then applied.

- the finite-corner obstruction is vacuous: rank 0, kernel 11;
- the stronger global product-ideal test `(surface,t^5)` certifies every one
  of the 11 candidates on all six E7 node charts at both `43` and `59`.

That was initially surprising, because a primitive nef isotropic K3 divisor
cannot have eleven global sections.  The missing condition was not at the E7
nodes.

### The missing degree-10 finite transition

Write the reduced rational function

```text
p = y(P1)/x(P1).
```

Its denominator factors as

```text
den(p) = h * u^2 * d0(u),
deg(d0)=10,
gcd(d0,u*h)=1.
```

At both `43` and `59`, `d0` has good squarefree reduction.  Away from `h=0`,
use the safe chord coordinate `r=m/h`.  Then

```text
q = r - c/d0,
```

with `c` a unit modulo `d0`.

The exact finite-transition matrix on the 11 q-frame candidates has

```text
rows = 450
rank = 11
kernel = 0
```

at both primes.

There is also a triangular characteristic-zero explanation:

```text
r^8 forces d0 | A9;
deg(d0)=10 > deg(A9/u^27)<=5, hence A9=0;
r^7 then forces A8=0;
r^6 then forces A7=0.
```

Therefore all eleven candidates vanish.

**Main lesson:** `q=(m-p)/h` is an excellent local h-adic coordinate but is
not a global base-regular function.  Promoting a local regular frame to a
global basis is exactly the error the module-intersection compiler must avoid.
The degree-10 `d0` divisor is a regression test for future frame changes.

## 5. Diagnostic marked-E7 layer scan under the q6^9 helper normalization

A uniform marked-chart scan was performed in the 18-generator frame

```text
m^b/t^6,       b=0..9,
x*m^b/t^8,     b=0..7.
```

This frame is convenient because it is individually admissible in the marked
chart.  The exact seven-component generic residue calculations gave identical
results at `43` and `59`.

The first nonzero generic leading layer appears at round 10:

```text
round 9  : rank 18, kernel 0
round 10 : rank 16, kernel <m^8,m^9>
round 11 : kernel <m^6,m^7,m^8,m^9>
round 12 : kernel <m^4,...,m^9, x*m^6,x*m^7>
round 13 : kernel <m^2,...,m^9, x*m^4,...,x*m^7>
round 14 : kernel <m^0,...,m^9, x*m^2,...,x*m^7>
round 15+: kernel all 18 directions.
```

The two round-10 vectors `t^4*m^8,t^4*m^9` also pass the composed q6-cover
membership calculation over `QQ` in that helper normalization.

This scan is useful diagnostic data, but **its t-exponents are not the final
q8 E7 elementary divisors**.  The normalization inherited the q6 module's
`-F` factor and then took its ninth tensor power.  The actual q8 divisor has a
single global fibre coefficient `-11F`; exceptional-cycle comparison alone
does not authorize counting the q6 fibre factor nine times and then applying
`-11F` again.

## 6. True reverse saturation exposed the adapted coordinate x-m^2

Starting from a safe mixed lattice in the helper normalization, the correct
reverse operation was tested: divide the lattice by one `t` and ask which
linear combinations remain E7-admissible.

### Saturation step 1

At both primes the obstruction matrix has

```text
rank = 10
kernel = 8.
```

The eight vectors are exactly

```text
x*m^j - m^(j+2),  j=0,...,7.
```

Thus the natural adapted coordinate is

```text
z = x - m^2.
```

This is not a random modular coincidence: the coefficients are `-1,+1` at
both primes.  It is also geometrically natural on the affine E7 component,
where `y^2=x^3` and `m=y/x`, hence `x=m^2`.

### Saturation step 2

After changing to the z-family, only two further directions divide by one
more `t`:

```text
t^4*z*m^3,
t^3*z*m^5.
```

Again the kernels are identical at `43` and `59`.

### Saturation step 3

The next reverse-saturation obstruction has

```text
rows = 119
rank = 18
kernel = 0
```

at both primes.

This proves saturation of the **helper-normalized mixed lattice** with respect
to the generic E7 conditions.  The structural basis change `z=x-m^2` is
reusable.  The absolute t-exponents of that helper lattice are superseded by
the fibre-normalization audit described next.

## 7. Fibre-normalization correction: keep -11F global and separate

The source q8 divisor has the literal vertical difference

```text
-11F + E7 correction + E8 correction.
```

The q6 E7 marked module, however, is a representative for

```text
O+(-P1)-F.
```

Its ninth power therefore carries `-9F` before the exceptional integral twist
is applied.  The exact identities

```text
c8 = 9*c6 + (2,5,6,4,6,3,5)     [E7 exceptional cycles]
```

and the analogous E8 cycle relation concern the exceptional classes.  They do
not by themselves identify the common base-fibre representative.

Consequences:

1. strip the helper q6 fibre factor when reconstructing the **reduced** E7 and
   E8 local lattices;
2. carry the q8 `-11F` exactly once as a global base twist;
3. do not use determinant values obtained by combining helper `t^9/u^9`
   factors with the literal `-11F` divisor;
4. the earlier provisional E7 determinant target `-98` is withdrawn;
5. the earlier statement that the helper-normalized E7 lattice was the final
   q8 E7 lattice is withdrawn.

This correction is essential before any global lattice intersection.

## 8. Reduced E7 normalization: exact first saturation layer

The reduced E7 reconstruction was restarted from the generic horizontal frame

```text
B = 1,m,...,m^9,x,xm,...,xm^7,
```

with the q6 helper fibre factor stripped.

Test `t^-1 B` against:

- the resolved exceptional E7 conditions with the helper `-9F` contribution
  removed; and
- the affine E7 component.

At the affine component `t=0` one has exactly

```text
y^2=x^3,
m=y/x,
x=m^2.
```

The first reduced saturation test gives, at both primes,

```text
exceptional rows = 0
affine rows = 10
rank = 10
kernel = 8.
```

The kernel is exactly

```text
x*m^j - m^(j+2),  j=0,...,7,
```

i.e. the same adapted coordinate `z=x-m^2`, now in the correctly reduced
normalization.

This is the **current valid E7 frontier**.  The next operation is reduced E7
saturation step 2 in the mixed `(m^b,z*m^b)` frame.  Do not reuse the absolute
helper-normalized exponents from section 6 without rederiving them here.

## 9. E8 target and singular integral module: what is exact, and what is not

The source q8 E8 component degrees and resolved cycle are exact.  In actual
chart order `(B1,B2,B3,B4,N3,N40,N4B,N4inf)` the exceptional cycle is

```text
(-2,-4,-6,-10,-4,-7,-5,-8).
```

In the integral II* coordinates

```text
u=1/t,
X=u^4*x,
Y=u^6*y,
Q=u^2*m=(Y-Y(P1))/(X-X(P1)),
```

the exceptional twist has complete singular-ring ideal

```text
I=(u^2,X,Y),
R/I has basis 1,u.
```

Inside the integral horizontal frame

```text
C=1,Q,...,Q^9,X,XQ,...,XQ^7,
```

put `s=Y(P1)/X(P1)`.  Since

```text
Q-s = (X(P1)*Y-Y(P1)*X)/(X(P1)*(X-X(P1))) in (X,Y),
```

the full preimage of the **singular integral ideal** `I` has basis

```text
u^2,
Q^b-s^b      (b=1..9),
X*Q^b        (b=0..7).
```

The exact checker reported

```text
preimage(I) determinant contribution = 2
C/B coordinate contribution = 178
helper u^9 contribution = 162
reported helper-normalized total = 342.
```

The algebraic preimage statement is valid.  The interpretation
"PASS_EXACT_Q8_E8_SATURATED_LATTICE" as the **resolved** E8 pushforward
lattice is too strong and is withdrawn.  The resolved II* surface can contain
functions not regular in the singular coordinate ring (for example valuation
checks already show that apparent divisions such as `X/u` can be regular on
all exceptional components).  E8 therefore needs a resolved-valuation
saturation step analogous to E7.

Likewise, the `u^9` factor is the q6^9 helper normalization and must be kept
separate from the literal global q8 `-11F` twist.

## 10. What the failed determinant checks taught us

Several provisional determinant checks were useful precisely because they did
not add up correctly.  They exposed two independent issues:

- individual endpoint floors are not the determinant of a saturated local
  lattice;
- exceptional-cycle tensor identities do not automatically align the common
  fibre representative.

Do not use the withdrawn `-98` E7 target or the helper-normalized E8 value
`342` as global pushforward determinants.

A valid determinant checksum should be applied only after:

1. reduced resolved E7 saturation is complete;
2. reduced resolved E8 saturation is complete;
3. both lattices use the same generic frame `B` and exclude helper fibre
   powers; and
4. the single global `-11F` twist is applied once.

Then the K3 identity `D^2=0` and `chi(O(D))=2`, together with generic fibre
rank 18, gives

```text
deg(pi_* O(D)) = 2 - 18 = -16.
```

That is the correct final vector-bundle checksum, not a shortcut for deriving
any local lattice.

## 11. Reusable compiler lesson

The central reusable result of this session is methodological.

A Riemann--Roch ambient should be assembled as an intersection of **local
lattices with explicit transition matrices**, not as independent degree
bounds on a chosen generic basis.

For this q8 problem the useful frames are:

```text
generic/open base:  B=(m^b, x*m^b)
E7 reduced:         adapted mixed frame using z=x-m^2
E8 reduced:         resolved II* frame in (u,X,Y,Q), still to saturate
h=0 collisions:     (q,X), q=(m-p)/h, X=h^2*x
```

The degree-10 `d0` example is the regression test: a frame can be perfectly
regular at one divisor and still introduce a hidden pole at another divisor.
A future module-intersection compiler should refuse to promote a local frame
to a global basis unless all transition denominators are accounted for.

## 12. Current frontier / next exact gates

As of this checkpoint, the next work should be:

1. **Reduced E7 saturation step 2.**  Start from the first reduced lift
   `z=x-m^2`; divide the resulting mixed lattice by `t` and compute the next
   exact/two-prime obstruction.  Continue until the reduced lattice is
   saturated.  The helper-normalized exponent table is diagnostic only.
2. **Resolved E8 saturation.**  Use the actual E8 blow-up valuation atlas and
   the exact q8 target cycle.  Start from the singular integral preimage basis
   above, enlarge by divisions that remain regular on every resolved
   component, and stop at a full-rank reverse-saturation obstruction.
3. **Choose one common reduced generic frame.**  Express the reduced E7 and
   E8 lattices in `B=(m^b,xm^b)` over `QQ(t)` with no inherited q6 fibre
   powers.
4. **Intersect with the exact h-adic lattice.**  Use the local `(q,X)` frame
   only at `h=0`; keep its full transition matrix.  Do not use `q` globally.
   In the safe `B` frame, the degree-10 `d0` divisor is a transition warning,
   not an extra ad hoc condition.
5. **Apply the global `-11F` twist once.**  Then compute the global section
   space.  The expected exact result is dimension two.
6. If `h0=2`, form the q8 pencil, eliminate to a genus-one quartic/Jacobian,
   verify the expected `D13/MW4` child, and only then resume the downstream
   H3 chain and rootless bisection-collision program.

## 13. Claim boundary

### Exact / structurally certified

- source q8 class, vertical difference, generic degree-18 RR basis;
- smooth h-adic q/X line-bundle lattice;
- E7 and E8 resolved target cycles;
- E7 actual component/function-field residue machinery;
- E8 singular integral ideal `(u^2,X,Y)` and its quotient;
- degree-10 `d0` finite transition and the rank-11 rejection of the fake
  global q-frame at both `43` and `59`;
- reduced-E7 first saturation kernel
  `<x*m^j-m^(j+2): 0<=j<=7>` at both primes;
- adapted coordinate `z=x-m^2`;
- helper-normalized reverse-saturation steps and their prime-independent
  kernels, as diagnostics only.

### Retracted / superseded

- treating `extra_e7_pole` as a route to the global q8 pencil;
- promoting the local h-adic q-frame to a global basis;
- the child-side fractional finite/infinity shortcut;
- `P0 smooth at IV* => identity component`;
- the provisional determinant target `-98`;
- the claim that the helper-normalized E7 saturation exponents were the final
  q8 E7 lattice;
- the claim that the singular-integral E8 preimage with helper `u^9` factor
  was already the complete resolved E8 pushforward lattice;
- using helper-normalized determinant `342` as a global q8 invariant.

### Still open

- complete reduced resolved E7 lattice;
- complete reduced resolved E8 lattice;
- exact common global lattice intersection and `h0(D)=2`;
- q8 pencil/equation and D13 child over `QQ`;
- equation-level continuation to the rootless H3 model;
- any rootless bisection extension collision or generic rank-19 claim.
