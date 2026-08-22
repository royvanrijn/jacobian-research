# Q80 low-q alternate corridor — 2026-08-22

## Status

This note records a new exact low-q continuation from the generic determinant-948 Q80 frame. It is a **secondary/fallback route** and does **not** replace the certified six-step rootless path in [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md).

The canonical route remains the only Q80 path currently certified all the way to the rootless MW17 frame. The new route is interesting because its first two newly discovered moves have much simpler equation geometry than the old pinned q12 move: both reduce to old-fibre degree two and have explicit CM24 binary-quartic models.

The exact generic lattice prefix now is

```text
E6+D5+A3 / MW3
  q4 (2,2)
D9+A4 / MW4
  q4 (2,2)
D7+D5 / MW5
  q6 (2,3)   [new escape]
D7+D4 / MW6
  q4 (2,2)   [orbit 424]
A6+A4 / MW7
  q4 (2,2)   [orbit 1222]
A6+A3 / MW8
```

The first two q4 arrows are the existing Q80 path. The last three arrows are the new corridor. Their machine-readable vectors are pinned in [`data/fibrations/kumar_q80_lowq_alternate_prefix.tsv`](data/fibrations/kumar_q80_lowq_alternate_prefix.tsv), and [`scripts/verify_q80_lowq_alternate_prefix.sage`](scripts/verify_q80_lowq_alternate_prefix.sage) replays the three new arrows exactly.

## 1. Reusable compiler result

The Q80 third-q12 local modules exposed an avoidable restriction in [`scripts/elliptic_neighbor_compiler.sage`](scripts/elliptic_neighbor_compiler.sage): `quotient_condition`, `resolved_chart_quotient_condition`, and `compile_resolved_conditions` assumed coefficient field `QQ`, even though `finite_ambient_image_condition` was already field-generic.

A field-generic compatibility layer is now stored in [`scripts/elliptic_neighbor_compiler_field_generic.sage`](scripts/elliptic_neighbor_compiler_field_generic.sage). It preserves historical `QQ` behavior and lets resolved/module intersections run over compatible exact number fields. The local regression which motivated it used

```text
L = QQ(sqrt(-6), sqrt(-3))
ambient dimension = 9
D7 rows = 4
D7+E6 rows = 5
all rows = 7
rank = 7
nullity = 2
```

and the selected/conjugate quadratic-I2 choices produced the same two-dimensional kernel. This is a reusable gain from the Q80 fallback even if the alternate route is later abandoned.

## 2. Second q4 is universal at equation level

The second Q80 q4 move was checked directly on the entire four-parameter Q80 ambient chart, rather than inferred from the pinned lattice frame.

With the transformed finite additive place

```text
L = W - 9(d-1)
```

one has identically

```text
ord_L(A,B,Delta) = (2,3,7)  -> I1* -> D5,
ord_infinity(A,B,Delta) = (2,3,9) -> I3* -> D7.
```

One exact `QQ(u)` witness has residual discriminant degree eight and squarefree residual factor. Hence the generic second child is universally

```text
I3* + I1* + 8 I1
ADE = D7 + D5
MW = 5.
```

Unlike the first q4, no extra collision condition is needed for the second q4. The first q4 still needs the rank-19 `I4 -> I5` collision.

## 3. Why the obvious low-q MW shortcuts failed

A saturated rank-five MW-quotient scan on `D7+D5/MW5` showed many low-height classes with q2/q3/q4/q5 shortest lifts. The pinned third-q12 MW class itself has height eight and a shortest full-frame lift of norm eight, i.e. q4.

However, scoring the cheapest representatives showed that these shortest lifts merely loop back to the same root system. In particular, no cheaper representative reproduced the productive root-rank-11 child of the pinned q12 move.

The important lesson was that the productive representative can be a **non-shortest root correction inside the same MW coset**.

## 4. q6 escape in the pinned q12 MW coset

Fix the pinned third-q12 MW coset in the generic `D7+D5/MW5` frame. Enumerating affine root-coset representatives by increasing norm gives the first escape already at q6:

```text
q = 6
(a,b) = (2,3)
v = (-5,-3,6,6,-8,-4,2,4,-1,8,-16,-1,0,3,5,-2,-2)
child = D7 + D4
root rank = 11
MW = 6
roots = 108
root determinant = 16.
```

This vector has the **same generic MW projection** as the old pinned q12 divisor.

### Chamber geometry

The raw q6 divisor has old-fibre degree three and `D.O=-1`. Deterministic chamber reduction gives

```text
reduced = (3,2,-5,-5,6,0,-4,-6,4,2,3,10,-14,-1,0,5,5,-2,-2)
D.F = 2
D.O = 1
MW norm = 8.
```

The horizontal section is exactly the same section as for the pinned q12 move, but the reduced divisor is now

```text
D = S + O + F + R
```

with integral vertical correction. Therefore the equation gate is a **marked chord / binary quartic**, not a marked trisection.

This is the key structural improvement over the old third q12 move.

## 5. CM24 q6 specialization and explicit equation

At CM24 the generic q6 child specializes to

```text
A1 + A1 + D6 + D8
root rank = 16
MW = 2.
```

The reduced q6 divisor specializes as

```text
D = Q + O + 2F + R,
Q = P1 + 3P2,
```

and the integral vertical correction is supported **only on the E6 fibre** of the old CM24 second child:

```text
E6 coefficients = -(2,3,4,3,2,2),
all A1 and D7 coefficients = 0.
```

Let `s^2=-6`, let the old second-child base be `W`, and let

```text
u = W + 27/2.
```

The known polynomial section `Q=P1+3P2` has chord `z_Q`. The exact local E6 jet forces

```text
lambda = -(39/4) s.
```

Indeed the branch discriminant of

```text
z_Q = 162 s + lambda u + V u^2
```

has its first three `u` coefficients zero automatically, while the `u^3` coefficient is

```text
(729/2) (4*(lambda/s)+39)^3.
```

Hence the exact new base coordinate is

```text
V = (z_Q - 162 s + (39/4)s*u) / u^2.
```

The residual genus-one quartic is

```text
(V^4 + 16/9*V^2 - 128/243*s*V - 64/243)*u^4
+ (-39*s*V^3 - 228*V^2 + 200/3*s*V + 992/27)*u^3
+ (648*s*V^3 + 18387/4*V^2 - 1710*s*V - 1174)*u^2
+ (6561/8*s*V + 2187/2)*u
+ 531441/64.
```

A short Weierstrass Jacobian is

```text
A(V) =
  68024448 V^6 - 160849476 s V^5 - 916676676 V^4
+ 446944068 s V^3 + 713569986 V^2 - 98644392 s V
- 33250608,

B(V) =
 -88159684608 s V^9 - 1876148288064 V^8
+2875406171544 s V^7 + 15062360165208 V^6
-8283115953504 s V^5 - 17889909720408 V^4
+4221255729708 s V^3 + 3780760503600 V^2
-324124340832 s V - 72965135232.
```

Its fibres are exactly

```text
V = 13s/18 : I2
V = 5s/9   : I2
V = 2s/9   : I4*
quadratic factor : 2 I1
infinity : I2*
```

so the specialized root system is `D8+D6+2A1`, root rank 16, MW2, exactly as predicted by the specialized lattice.

The exact CM24 equation check is consolidated in [`scripts/verify_q80_lowq_cm24_equations.sage`](scripts/verify_q80_lowq_cm24_equations.sage).

## 6. Exact low-q landscape from D7+D4/MW6

A full vector-by-vector q-shell scan is unnecessarily large. Quotienting by the exact Weyl group of the parent root system gives the following complete q<=4 result:

| q | raw vectors up to sign | Weyl orbits | productive orbits |
|---:|---:|---:|---:|
| 2 | 2,952 | 37 | 0 |
| 3 | 45,996 | 156 | 0 |
| 4 | 344,363 | 453 | 12 |

Every q2 and q3 orbit loops to `D7+D4`. Every `(1,4)` q4 presentation also loops. The productive moves occur only for `(a,b)=(2,2)` and all reduce root rank `11 -> 10`.

The twelve productive q4 orbits have four child types:

```text
A6+A3+A1   (2 orbits)
A6+A4      (4 orbits)
A7+A3      (2 orbits)
D6+A3+A1   (4 orbits).
```

The summary is pinned in [`data/fibrations/kumar_q80_lowq_weyl_summary.tsv`](data/fibrations/kumar_q80_lowq_weyl_summary.tsv).

## 7. Preferred q4 from D7+D4: orbit 424

Equation-geometry scoring of all twelve productive q4 orbits showed that every productive divisor reduces to old-fibre degree two.

The preferred generic choice is orbit 424:

```text
q = 4
(a,b) = (2,2)
v = (32,48,-21,28,8,-52,-34,0,18,5,-23,43,9,-18,16,-6,-6)
child = A6 + A4
root rank = 10
MW = 7
roots = 62
root determinant = 35
D.F = 2
D.O = 0
MW height = 13/4
shortest section norm = 6
section P.O = 1
fiber twist = 1.
```

Its generic vertical correction is integral and supported on both old root components, with total coefficient L1 norm 12.

### CM24 specialization: the marked section becomes 2-torsion

The CM24 q6 source has fibres `D8+D6+2A1`. Orbit 424 specializes to

```text
A7 + A7
root rank = 14
MW = 4.
```

The generic marked section survives as a section but its specialized Shioda projection is **zero**. Thus it becomes torsion. A node-constrained search over `GF(73)` found a unique compatible `Y=0` polynomial section, and exact factorization over `QQ(sqrt(-6))` identifies it as the rational 2-torsion point

```text
T_y = 0,
T_x = 1944*s*V^3 + 12150*V^2 - 4401*s*V - 3036.
```

Its reduction under `s -> 33 mod 73` is

```text
T_x = 30 + 37V + 32V^2 + 58V^3,
```

exactly the unique modular candidate.

Hence the degree-two generic-fibre Riemann--Roch basis is simply

```text
1,
z_T = Y/(X-T_x).
```

The orbit-424 vertical correction uses the other finite I2 fibre, at `V=13s/18`. Therefore an exact new base is

```text
U = Y / ((X-T_x)*(V-13s/18)).
```

The resulting residual quartic is

```text
-204073344 V^4
+ (-11664*s*U^2 + 206907696*s) V^3
+ (U^4 - 72900*U^2 + 483965604) V^2
+ (-13/9*s*U^4 + 26406*s*U^2 - 84925584*s) V
- 169/54*U^4 + 18216*U^2 - 33662304.
```

Its Jacobian has

```text
A(U) =
 -27 U^8 - 4251528 U^6 - 81616583016 U^4
 -27113235502176 U^2 - 9882774340543152,

B(U) =
 54 U^12 + 12754584 U^10 + 746946702792 U^8
 +6181817694496128 U^6 + 5751774666196114464 U^4
 +1556181178759286886528 U^2
 +378152026438506713426304.
```

The discriminant factors as

```text
U^8 *
(U^8 + 314199/2 U^6 + 1905215985 U^4
 + 509070522546 U^2 + 183014339639688),
```

with the second factor squarefree. Infinity contributes another order eight. Therefore the fibres are exactly

```text
I8 at U=0,
I8 at infinity,
8 I1,
```

so the specialized child is `A7+A7/MW4`, exactly matching the lattice specialization.

## 8. Exact low-q landscape from A6+A4/MW7

The same exact Weyl-orbit strategy gives:

| q | raw vectors up to sign | Weyl orbits | productive orbits |
|---:|---:|---:|---:|
| 2 | 2,812 | 78 | 0 |
| 3 | 46,610 | 466 | 0 |
| 4 | 345,555 | 1,668 | 6 |

Again all q2 and q3 moves loop, all `(1,4)` q4 presentations loop, and the first productive moves are `(2,2)` q4 neighbours. Every productive move drops root rank `10 -> 9`, hence MW `7 -> 8`.

The six productive q4 orbits are

```text
A6+A3          : orbits 205, 718, 1222, 1632
A5+A2+2A1      : orbits 305, 720.
```

This repeated pattern is now a useful search heuristic: after the exceptional q6 escape, `(2,2)` q4 neighbours form a root-peeling corridor while q2/q3 and `(1,4)` q4 moves are loops in the two frames checked so far.

It is not yet a theorem and should not be extrapolated without exact scans on later frames.

## 9. Preferred rank-9 continuation: orbit 1222

Generic equation-geometry scoring gives two interesting candidates:

```text
orbit 1222:
  child A6+A3 / MW8
  MW height 122/35
  section P.O=1
  two vertical root components
  vertical L1 = 9

orbit 720:
  child A5+A2+2A1 / MW8
  MW height 208/35
  section P.O=2
  one vertical component, coefficient 1.
```

CM24 specialization resolves the choice. Both marked sections specialize to ordinary non-torsion sections with `P.O=1`, so orbit 720 loses its hoped-for horizontal simplification.

For orbit 1222:

```text
v = (10,53,-192,-114,29,-256,-170,-12,-14,74,-32,-14,-6,-26,-58,84,-28)
generic child = A6+A3 / MW8
special child = 2A6 + 3A1
special section P.O = 1
special section MW height = 25/8
special D.F = 2
special D.O = 0
special vertical support = one A7 fibre
coefficients = (-1,-1,-1,0,0,0,0).
```

For orbit 720 the specialized vertical support is also one A7 fibre but has five nonzero coefficients. Therefore orbit **1222 is the retained candidate**.

The CM24 root signature `2A6+3A1` is the same as the old pinned third-q12 CM24 child. This is an observation only. It is **not yet proved** that the two constructions reach the same CM24 elliptic fibration or differ merely by a base/MW automorphism.

## 10. What is exact, and what remains open

### Exact now

- the original q4,q4 prefix;
- universal equation-level `D7+D5` structure of the second q4;
- the q6 escape vector and generic child `D7+D4/MW6`;
- q6 chamber reduction to degree two with the same marked horizontal section as the old q12;
- the CM24 q6 binary-quartic coordinate and exact Weierstrass model;
- the complete q<=4 Weyl-orbit landscape from `D7+D4/MW6`;
- orbit-424 generic q4 lattice/chamber data;
- orbit-424 CM24 2-torsion specialization and exact q4 equation;
- the complete q<=4 Weyl-orbit landscape from `A6+A4/MW7`;
- orbit-1222 generic lattice/chamber data and its CM24 specialization;
- field-generic exact module-intersection machinery via the compatibility layer.

### Still open

1. Continue the new `A6+A3/MW8` branch to rootless and compare total complexity with the certified old suffix.
2. Test whether the repeated `(2,2)` q4 root-peeling behavior continues.
3. Prove or disprove that orbit-1222's CM24 `2A6+3A1` child is the same fibration as the old q12 CM24 `2A6+3A1` child.
4. Algebraize the q6 marked horizontal section on the generic characteristic-zero Q80 family; the CM24 equation alone does not supply the generic `QQ(u)` rational function.
5. Lift the later q4 models from CM24 to the generic family once the q6 generic marked section is available.
6. Prove full nefness for any retained later corridor divisors when the current chamber scoring has only used the declared old fibre/section machinery.
7. Do not replace the canonical rootless Q80 route until the new corridor has both a complete rootless lattice suffix and a credible characteristic-zero equation path.

## 11. Recommended live branch

For subsequent work, use

```text
D7+D5/MW5
  -- q6 (2,3), escape vector --> D7+D4/MW6
  -- q4 (2,2), orbit 424    --> A6+A4/MW7
  -- q4 (2,2), orbit 1222   --> A6+A3/MW8
```

as the preferred low-q exploratory branch.

The old q12 route remains the certified fallback to rootless.
