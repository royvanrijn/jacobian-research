# Elkies rank-17 K3 reconstruction — research update (2026-08-20)

> **Superseded source marking (2026-08-21).** The H2/`H8 cap H237` route
> described below remains a valid downstream comparison, but it is not the
> recovered source polarization.  The exact source is the third Kumar frame
> `H3=[[21/2,3],[3,46]]`: its `H21 cap H92` component has been normalized over
> `QQ` to the published level-474 genus-two curve.  The marked H3 corridor and
> direct rootless endpoint gate are now complete; the current programme gate
> is residual 2-descent in the compact published `t` chart.  See
> [`README.md`](README.md) and
> [`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).

## 2026-08-21 upstream correction: backtrack through Kumar, not a guessed CM chart

The primary sources determine the order of construction more tightly than the
earlier CM-deformation programme used.  The non-CM point on `X(6,79)` first
gives a principally polarized QM abelian surface, whose Dolgachev--Kumar K3 has
a canonical `E7+E8` fibration of MW rank two and regulator `474`.  Neighbor
transformations then lead to the rootless MW-rank-17 fibration.

An exact finite classification of determinant-474 binary height lattices and
their `E7` discriminant glue leaves exactly three frames in the recovered
determinant-948 genus:

```text
[[5/2,1],[1,190]], [[4,0],[0,237/2]], [[21/2,3],[3,46]].
```

Each has root system `E7+E8` and cyclic discriminant group `Z/948`.  The middle
one is the unique frame with an extra height-lattice involution.  Decomposing
the two elliptic quotients of Elkies's genus-two curve identifies that
involution as `w2=w237`, its hyperelliptic involution, so the middle frame is
the Kumar anchor used by the quotient construction.  See
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).

Consequently, the `E8+A2^3` deformation remains a useful secondary CM chart,
but it is no longer the primary reconstruction assumption.  The decisive
missing datum is now the rational `t`-line to `H2` Kumar/Clebsch--Igusa map
and the section descent over `u`.

The backtrack now also fixes both CM boundary equations.  Primitive closure
of the middle Kumar frame gives `E8+E8+A2`, MW zero, at discriminant `-3`, and
`E8+E8`, MW lattice `diag(4,6)`, at discriminant `-24`.  Their standard Inose
models are

```text
Y^2 = X^3 + T^5*(T-1)^2,
Y^2 = X^3 - 51*T^4*X + T^5*(T^2-92*T+1),
```

respectively, with the second determined up to rational quadratic twist.  For
the H2 comparison, the remaining problem is local deformation/interpolation
from these marked boundary fibers and recovery of the twist and section
descent selected by `u`, not another broad equation search.

The H2 height decomposition identifies its ambient modular geometry even more
precisely.  The height-4 section with the `E7` frame has determinant 8, while
the height-`237/2` section alone has determinant 237.  This comparison locus
is therefore `H8 cap H237`, not the corrected H3 source curve.
Elkies--Kumar's published ancillary data supply the
complete two-parameter `H8` Kumar equation, its Clebsch--Igusa map, and its
oriented double cover.  The unrecovered datum has dropped to one equation: the
discriminant-237 divisor in that explicit plane, normalized against the known
`(t,u)` Shimura model.

The first useful downstream step from this H2 anchor is also pinned.
Constraining the high-section coordinate makes the norm-120 shell exact and
small: 56 sign-pairs in 441 enumeration nodes.  A `q=60` neighbor with
`(a,b)=(5,12)` gives `E8+E6/MW3`, with reduced height Gram
`[[4,0,0],[0,20/3,1],[0,1,12]]`.  Thus the corrected geometric route begins

```text
Kumar E7+E8/MW2 --q=60--> E8+E6/MW3,
```

not with the guessed split E6 chart.  The remaining gap at this arrow is the
birational Weierstrass transformation, not its integral lattice existence.

## Decisive correction: the all-IV CM deformation is obstructed

The abstract inherited root frame `E8+A2^3` remains exact, but its previous
promotion to `II*+3 IV+II` is impossible for the target.  That Weierstrass
family has `j=0`; its order-three CM automorphism forces the geometric
Mordell--Weil rank to be even.  Shioda--Tate requires the target rank to be
`19-2-14=3`.  See
[`E8_A2_KODAIRA_CORRECTION.md`](E8_A2_KODAIRA_CORRECTION.md).

An `A2` root lattice may be `I3` or `IV`.  At least one of the three target
`A2` fibers must be `I3`, and the correct deformation must allow nonzero
`A(t)`.  The old all-IV section systems are historical only.

A viable mixed lift has since been derived directly:

```text
D=t(t-1),
A=-3r^2D^2,
B=D^2((t-lambda)^3-2r^3D),
```

with fibers `II*+2 IV+I3+3 I1`, nonconstant `j`, and the correct CM endpoint.
See [`E8_A2_MIXED_FAMILY.md`](E8_A2_MIXED_FAMILY.md).  Its target
determinant-948 sublocus is not yet identified.

The exact NS glue supplies the missing section profiles in the reduced basis:
`(1,1,0)`, `(0,2,0)`, and `(0,0,0)` across the three `A2` factors.  Their
nonzero counts recover heights `8/3,10/3,4` and force all three generators to
be polynomial sections (`P_i.O=0`); their pairwise section intersections are
all `2`.

## Final correction: the exact E6 neighbor path is recovered

The actual frame chain is now pinned:

```text
rank17 --q=90--> MW7 --q=4--> MW4 --q=4--> E6/MW3.
```

An exact verifier reconstructs each committed child frame from its primitive
isotropic vector and checks a determinant-one composite Neron--Severi
transport.  See [`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).  The remaining
gap is geometric: the earlier split E6 Weierstrass chart was inferred from the
fiber type and was not obtained by applying these neighbor operations to the
original K3.  Further large scans of that chart are therefore stopped in
favor of backtracking the genuine chain.

## Late correction: exact frame glue and the missing height gate

The E6 frame itself remains exact, but the first modular searches did not
enforce the off-diagonal Mordell-Weil heights.  Direct recovery from the
17-dimensional frame selects the canonical component orbit and excludes the
other height-only orbit.  For the selected P1/P2 profiles, the missing target
condition is `(P1+P2).O=1`, equivalently `<P1,P2>=-5/6`.

Compiled exhaustive rational-chart scans over
`GF(5),GF(7),GF(11),GF(13),GF(17)` tested 406,655,040 cores.  All 69 sections
that passed the older P2 component tests fail the full height gate.  Thus the
next E6 step is an algebraic backtrack to a corrected P1/P2 system, not a
larger run of the old scan.  Exact scripts and counts are in
[`E6_MW3_ATTACK.md`](E6_MW3_ATTACK.md).

This note records the current state of the `elkies-k3` branch of the project after the recent reconstruction work. It is intentionally split into exact/computationally verified results, failed approaches, and the current frontier.

**Later-day consolidation.** The arithmetic identification and lattice results
below remain current.  Both equation-level shortcuts were later corrected:
the E6 split chart was not geometrically transported, and the all-IV
`E8+A2^3` chart has the parity obstruction above.  The consolidated frontier
is in [`RECONSTRUCTION_PROGRESS.md`](RECONSTRUCTION_PROGRESS.md).

## Executive summary

The project has moved from a broad “recover Elkies's missing rank-17 K3 model from its Mordell–Weil lattice” search into a much more constrained reconstruction problem.

The main progress is:

1. The recovered rank-17 Mordell–Weil lattice has been connected exactly to the quaternionic Shimura datum `(D,N)=(6,79)` via an explicit Eichler order and inverse-Clifford construction.
2. The distinguished CM endpoint has been identified arithmetically with CM discriminant `-3`; transporting the corresponding Gross vector back into the K3 transcendental lattice gives a rank-20 specialization with transcendental lattice `A2`, determinant `3`.
3. The discriminant-3 CM K3 is therefore the classical singular K3 `X3`, for which explicit elliptic fibrations are known. Utsumi No.1,

   `y^2 = x^3 + t^5 s^5 (t-s)^2`,

   has fibers `II* + II* + IV`, ADE type `E8^2 + A2`, and MW group `0`.
4. The generic discriminant-948 rank-19 Néron–Severi lattice embeds in this CM lattice as the orthogonal complement of a primitive class of square `-316`. In the inherited elliptic fibration, one `E8` breaks to `A2 + A2`, leaving generic fiber root lattice

   `E8 + A2 + A2 + A2`

   and generic MW rank `3`.
5. The exact rank-3 MW height lattice has now been recovered, including the missing 3-adic glue. In a reduced basis it is

   ```text
   (1/3) * [ 8 -1  0 ]
           [-1 10  0 ]
           [ 0  0 12 ]
   ```

   with determinant `316/9`.
6. A formerly proposed small family,

   `y^2 = x^3 + [t(t-1)(t-lambda)]^2 (t-mu)`

   does specialize to Utsumi No.1, but is now exactly rejected for the target:
   it has constant `j=0` and hence even geometric MW rank, whereas the target
   frame requires rank `3`.  The correct family must realize at least one
   `A2` as `I3` and retain a nonzero `A(t)`.

This is a substantial reduction from the original 17-section reconstruction problem.

---

## 1. Recovered rank-17 lattice

The recovered positive-definite rank-17 Mordell–Weil lattice has determinant

`948 = 2^2 * 3 * 79`

and the expected short-vector fingerprint: `1311` pairs of height-4 sections.

The lattice data is in:

- `elkies-k3/data/lattice/rank17_gram.txt`
- `elkies-k3/data/lattice/all_1311_in_short_basis.txt`
- `elkies-k3/data/lattice/short_vector_basis_gram.txt`

The relation structure among the 1311 pairs is also highly nontrivial and produced dense additive motifs. In particular a rank-3 motif generated by three minimal vectors `P,Q,R` contains the ten sections

`P,Q,R,-P,-Q,-R, +/- (P+Q), -(P+R), -(P+Q+R)`.

That motif was useful for reducing the original polynomial reconstruction systems, but not enough to make direct Groebner solving practical.

---

## 2. Direct polynomial reconstruction attempts that failed

Two direct approaches were tested and should currently be regarded as exhausted:

### Three-hub system

The relation-compressed section system had roughly

- 167 variables
- 251 polynomial equations

A finite-field `msolve` attack over `GF(101)` timed out after two hours even on the first sliced probes.

### Rank-3 motif system

Using only `P,Q,R` and the associativity diamond

`(P+Q)+R = (P+R)+Q`

reduced the problem to about

- 48 variables
- 67 equations

but this still timed out over `GF(101)`.

Conclusion: direct reconstruction of the rank-17 model from many polynomial sections is the wrong representation. More compute is unlikely to fix this efficiently.

---

## 3. Generic transcendental lattice and exact quaternionic model

For a K3 surface with Picard rank `19`, the transcendental lattice has rank `3`, signature `(2,1)`, and determinant `-948`.

A search over even ternary forms with:

- determinant `-948`
- signature `(2,1)`
- cyclic discriminant group of order `948`
- even Clifford algebra ramified at `2,3`

produced many coordinate representatives. Exact discriminant-form checks and later explicit isometries showed that these representatives lie in the same integral class relevant to the recovered K3 lattice.

The exact quaternion algebra is

`B = QuaternionAlgebra(6) = (-1,3)_Q`,

ramified at the finite primes `2` and `3`.

Because Sage 10.9 cannot directly construct an Eichler order of level `79` in a quaternion algebra ramified at two finite primes, the order was built manually as the intersection of two maximal orders. An element of norm `-79`

`alpha = -2 - 2*i + 5*j + 2*k`

produced an Eichler order of index exactly `79` in a maximal order:

```text
O = <
  1/2 + 1/2*i + 1/2*j + 37/2*k,
  i + 60*k,
  j + 56*k,
  79*k
>
```

This is an exact `(D,N)=(6,79)` Eichler order.

---

## 4. CM endpoint: discriminant -3

Eichler optimal-embedding counts were computed for imaginary quadratic orders. Up to `|Delta| <= 20000`, the only orders producing exactly four optimal-embedding classes were

- `Delta = -3`
- `Delta = -24`

Both occur in the exact Gross lattice of the `(6,79)` Eichler order.

The exact Gross lattice `O^T` was constructed and primitive representatives found:

```text
Delta=-3:
  beta = -3*i - j + k
  norm(beta)=3

Delta=-24:
  beta = -78*i + 38*j - 24*k
  norm(beta)=24
```

The Atkin–Lehner normalizer action provided a useful distinction: the norm-3 vector is fixed by the `w_3` representative, while the norm-24 vector is not.

Both embeddings have now been transported through the exact inverse-Clifford
correspondence into the recovered K3 transcendental lattice.

For `T=0`, using an explicit integral isometry from the inverse-Clifford lattice to the recovered ternary lattice, the CM functional produced the primitive K3 vector

```text
v = (81,95,-52)
```

with

```text
v^2 = -316
q_T(v) = -158
div(v) = 316
```

and orthogonal complement

```text
[26 7]
[ 7 2]
```

which reduces to the binary form

```text
(1,1,1)
```

with discriminant `-3` and determinant `3`.

The determinant identity checks exactly:

`det(v^perp) = 948 * 316 / 316^2 = 3`.

Therefore the relevant CM specialization is the singular K3 surface with transcendental lattice

```text
T_CM = [2 1]
       [1 2]
```

namely `A2`, determinant `3`.

This is now an arithmetic/computational identification, not merely a heuristic CM-label match.

For the norm-24 vector the same exact transport gives

```text
v = (70,86,-3)
v^2 = -158
q_T(v) = -79
div(v) = 79
v^perp = [4 0]
          [0 6]
```

The complement has determinant `24`, binary form `(2,0,3)`, and discriminant
`-24`; the determinant identity is

`det(v^perp) = 948 * 158 / 79^2 = 24`.

The coordinate geometry identifies the two CM orbit shapes.  In the chart
`s=1/t, v=u/t^3`, `w3` fixes `(s,v)=(0,+/-4)`, matching the `w3`-fixed
norm-3 class.  The four points `(t,u)=(+/-2,+/-32)` have no stabilizer in the
visible Klein four-group and match the norm-24 class.  Thus the correct Kumar
component has exact singular-K3 anchors

```text
t=infinity : Delta=-3,  T=[[2,1],[1,2]]
t=+/-2     : Delta=-24, T=[[4,0],[0,6]].
```

Relevant scripts/results:

- `construct_exact_gross_lattice.sage`
- `atkin_lehner_cm_orbits.sage`
- `map_clifford_to_k3_T.sage`
- `transport_cm_delta3_to_k3.sage`
- `artifacts/local/elkies-k3/cm-delta3-k3-vector.txt`
- `artifacts/local/elkies-k3/cm-delta24-k3-vector.txt`

---

## 5. Explicit discriminant-3 K3 endpoint

The discriminant-3 singular K3 is well known and has explicit Jacobian fibrations.

Two were independently checked.

### Utsumi No.1 — preferred deformation anchor

Homogeneous model:

```text
y^2 = x^3 + t^5 s^5 (t-s)^2
```

Affine `s=1`:

```text
y^2 = x^3 + t^5 (t-1)^2
```

Discriminant:

```text
Delta = -432 * t^10 * (t-1)^4
```

Fibers:

- `II*` at `t=0`
- `II*` at `t=infinity`
- `IV` at `t=1`

so the ADE root lattice is

`E8 + E8 + A2`

and MW group is trivial. The trivial lattice already has discriminant `3`.

### Utsumi No.2 — independent cross-check

```text
y^2 = x^3
      - 3 t^2(t^6-16t^3+16)x
      + 2 t^3(t^3-2)(t^6+32t^3-32)
```

Fibers:

- `I12*`
- `I3`
- `3 I1`

ADE type `D16 + A2`, MW torsion `Z/2`; Shioda's discriminant formula again gives `|disc(NS)|=3`.

These are verified by:

`elkies-k3/scripts/verify_disc3_k3_fibrations.sage`.

---

## 6. Generic NS as a codimension-one sublattice of the CM surface

For Utsumi No.1,

`NS_CM = U + E8(-1)^2 + A2(-1)`

with determinant `-3`.

An explicit primitive vector of square `-316` and divisibility `1` was chosen inside one `E8` summand:

```text
(9,-2,-10,-7,-4,4,3,-7)
```

The orthogonal complement inside `NS_CM` has:

- rank `19`
- determinant `948`
- cyclic discriminant form isometric to the recovered generic target `U + (-M_17)`.

This gives an explicit realization of the generic rank-19 lattice as a codimension-one lattice polarization inside the discriminant-3 CM K3.

Crucially, the chosen class breaks one of the two `E8` root systems. The surviving roots inside the broken `E8` form exactly

`A2 + A2`.

The other `E8` and the original `A2` remain untouched.

Therefore the inherited generic elliptic fibration has reducible fiber root lattice

`E8 + A2 + A2 + A2`,

root rank `14`.

By Shioda–Tate:

`rank(MW) = 19 - 2 - 14 = 3`.

This is a much more useful intermediate fibration than either the CM rank-0 fibration or the final rootless rank-17 fibration.

Verified by:

- `embed_generic_ns_in_disc3_cm.sage`
- `analyze_inherited_rank3_fibration.sage`

---

## 7. Exact rank-3 Mordell–Weil lattice

The first orthogonal-intersection computation produced the integral rank-3 essential lattice

```text
[-4  -4   20]
[-4 -34  -73]
[20 -73 -412]
```

or positive version

```text
[ 4  4 -20]
[ 4 34  73]
[-20 73 412]
```

with determinant `2844`.

At first this appeared inconsistent with Shioda's expected regulator

`948 / 27 = 316/9`.

The ratio is exactly

`2844 / (316/9) = 81 = 9^2`.

The resolution is that the orthogonal intersection is an index-9 sublattice of the actual MW lattice; the missing structure is glue from the full NS lattice.

The exact quotient was then computed directly:

```text
S / (Triv + C3) ~= (Z/3)^2
```

with index `9`.

Projecting all nine NS glue cosets into the essential space gives the true MW lattice. A basis in the old `C3` coordinates is

```text
[1/3  0   2/3]
[ 0  1/3  2/3]
[ 0   0    1 ]
```

and the exact height Gram matrix is

```text
[524/3 586/3 268]
[586/3 658/3 299]
[  268   299  412]
```

with determinant

`316/9`.

After multiplying by 3 and lattice reduction, this becomes

```text
[ 8 -1  0]
[-1 10  0]
[ 0  0 12]
```

so a convenient reduced MW height Gram is

```text
(1/3) * [ 8 -1  0]
        [-1 10  0]
        [ 0  0 12]
```

This is currently the most useful arithmetic fingerprint for reconstructing the intermediate explicit family.

Verified by:

`elkies-k3/scripts/recover_rank3_mw_via_ns_glue.sage`.

---

## 8. Rejected all-IV family candidate

The abstract inherited root configuration is

`E8 + A2^3`.

It was previously promoted to the short Weierstrass family

```text
y^2 = x^3 + [t(t-1)(t-lambda)]^2 (t-mu).
```

Generically this has:

- `II*` at infinity
- `IV` at `t=0`
- `IV` at `t=1`
- `IV` at `t=lambda`
- `II` at `t=mu`

At the CM endpoint

```text
lambda=0
mu=0
```

it specializes to

```text
y^2 = x^3 + t^5(t-1)^2,
```

exactly Utsumi No.1.  Nevertheless, this family cannot contain the intended
exact `E8+A2^3`, Picard-rank-19 fibration.  Its `j`-invariant is zero and the
order-three automorphism makes the geometric MW rank even, while the target
rank is three.  The error was treating every `A2` root factor as an `IV`
fiber; `I3` is the other possibility.  See
[`E8_A2_KODAIRA_CORRECTION.md`](E8_A2_KODAIRA_CORRECTION.md).

---

## 9. Historical all-IV section equations

The equations in this section are retained to document the rejected
experiment.  They must not be used to reconstruct `X(6,79)`.

A useful degree-reduction identity for polynomial sections is:

```text
x = q(t)^2 + r(t)
y = q(t)^3 + s(t)
```

with

- `deg q = 2`
- `deg r <= 1`
- `deg s <= 1`.

For

```text
y^2 = x^3 + [t(t-1)(t-lambda)]^2(t-mu)
```

this gives exactly eight coefficient equations in nine unknowns

```text
lambda, mu,
q0,q1,q2,
r0,r1,
s0,s1.
```

This is already dramatically smaller than the failed 48/167-variable reconstruction systems.

The coefficient system also has useful triangular structure:

- the highest coefficient equation (`k=7`) is linear in `s1`;
- after substituting `s1`, `k=6` is linear in `mu`;
- after eliminating those, `k=5` is only quadratic in `lambda`.

A triangular exporter/solver has therefore been added:

- `export_rank3_jump_triangular.sage`
- `run_rank3_jump_triangular.py`

This elimination is no longer an active experiment because its ambient family
is parity-obstructed.

---

## 10. Superseded 2026-08-20 frontier

The exact lattice backtrack has recovered the rank17-to-E6 neighbor chain, and
the parity audit has rejected the all-IV CM deformation.  These together
locate the missing step: derive the non-isotrivial Kodaira lift of the exact
`E8+A2^3` root frame.

The first lift is now explicit: the two-parameter
`II*+2 IV+I3+3 I1` family in
[`E8_A2_MIXED_FAMILY.md`](E8_A2_MIXED_FAMILY.md).

The then-current goal was:

1. Impose on the mixed two-parameter family the full reduced
   MW height/glue matrix, not merely one arbitrary section.
2. Identify the one-dimensional component with the exact determinant-948
   Neron--Severi lattice and `X(6,79)` period data.
3. Use the pinned neighbor chain to produce and cross-check the E6 and
   rootless rank-17 fibrations.
4. Retain the other non-isotrivial `I3/IV` distributions as fallbacks if the
   mixed family fails the exact lattice gates.
5. Reproduce E29 as a mechanism calibration, then search for specializations
   beyond the now-public rank-at-least-30 record.

This is now a concrete algebraic reconstruction programme rather than an unconstrained search.

---

## 11. What not to spend compute on right now

The following approaches produced enough negative evidence that they should remain paused:

- direct point searches on the fixed E29 curve;
- larger rational-point boxes on fake-2-descent covers;
- direct 167-variable rank-17 section systems;
- direct 48-variable rank-3 motif Groebner systems;
- broader ternary-lattice genus searches;
- any larger scan of the split E6 chart;
- any further solve of the all-IV `j=0` rank-jump system;
- larger Heegner vector bounds without a new arithmetic discriminator.

The current bottleneck is choosing and deriving the correct non-isotrivial
Kodaira lift, not raw compute.

---

## 12. Important exact checkpoints

For quick restart/review, the strongest exact checkpoints are:

```text
rank-17 MW lattice determinant        = 948
quaternion discriminant               = 6
Eichler level                         = 79
exact Eichler-order index             = 79
CM order used                         = Delta -3
CM K3 transcendental determinant      = 3
CM K3 transcendental lattice          = A2
extra CM class square                 = -316
generic NS determinant                = 948
broken E8 root system                 = A2 + A2
inherited generic ADE                 = E8 + A2^3
inherited generic MW rank             = 3
NS glue quotient                      = (Z/3)^2
C3 -> true MW index                   = 9
exact MW regulator                    = 316/9
reduced 3*height Gram                 = [[8,-1,0],[-1,10,0],[0,0,12]]
```

The main conceptual reduction is:

```text
1311 height-4 sections / rank-17 inverse problem
    -> exact (6,79) quaternion/Eichler arithmetic
    -> CM discriminant -3 endpoint
    -> explicit discriminant-3 K3
    -> rank-3 E8+A2^3 abstract intermediate frame
    -> reject the all-IV j=0 lift by MW-rank parity
    -> explicit II*+2IV+I3+3I1 two-parameter family
    -> determinant-948 one-dimensional X(6,79) locus (current target)
```
