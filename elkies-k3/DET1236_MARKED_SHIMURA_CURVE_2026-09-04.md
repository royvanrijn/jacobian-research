# Determinant 1236: exact marked Shimura curve and remaining lift obstruction

<!-- status-consumer: EC-K3-DET1236-GENUS2-RATIONAL-POINTS 5a3c84eb9f7f0604 -->
<!-- status-consumer: EC-K3-DET1236-MARKED-SHIMURA-CURVE 185d31609e7702fc -->

Date: 2026-09-04.

Status: **UNRESOLVED_FOR_EXPLICIT_REASON**.

## Arithmetic certificate

The selected rootless-MW17 row is

```text
surface       K3-6d288cfad55e0d15
legacy NS     NS0035
determinant   1236
frame         K3-6d288cfad55e0d15-F001
```

Its exact projective stable discriminant-kernel period curve is

```text
C_1236 = X_0^6(103)/<w_618>,     genus 6.
```

This is no longer an unidentified marked curve. It has no cusps and has two
certified rational CM points of discriminant `-3`. Those points are
Picard-rank-20 specializations, so neither realizes the requested saturated
rank-19 Neron--Severi lattice. No rational non-CM point on `C_1236` is
currently certified.

The rational points on the displayed genus-two quotient are now complete:
there are exactly fourteen. The precise remaining problem is to construct
the degree-two map from `C_1236` to that quotient and decide which of its
twelve non-fixed rational points have rational non-CM lifts. A positive lift
would make the row `ARITHMETICALLY_REALIZABLE`; a proof that the two CM
points exhaust `C_1236(QQ)` would make it `ARITHMETICALLY_EXCLUDED`. Until
then no equation work is authorized.

## Literal lattice and discriminant form

The transcendental lattice is literally

```text
T = [ -2   0    1 ]
    [  0   4    0 ]
    [  1   0  154 ],

signature(T) = (2,1),       det(T) = -1236,       content(T)=1.
```

Its Smith form is `diag(1,1,1236)`, so `A_T=Z/1236Z`. For the exact dual
generator

```text
g = (1/309, -1/4, 2/309),
```

the discriminant values are

```text
b(g,g) = 317/1236  mod ZZ,
q(g)   = 317/2472  mod ZZ.
```

Direct enumeration gives

```text
O(A_T) = {1,205,413,617,619,823,1031,1235},
```

acting by multiplication on the cyclic generator.

The catalogue frame `F001` has rank 17, determinant 1236, minimum squared
norm four, and no roots. Its Gram hash is

```text
1fb6824ca519210cbfe32af3c0c4c15fb1d3867c5847cca32878e5a4929b0682.
```

That geometric information selects the row but is not used to infer a
rational marking.

## Literal even Clifford order

In the basis `1,e0e1,e0e2,e1e2`, the literal even Clifford order has reduced
trace pairing

```text
[ 2   0    1     0 ]
[ 0   4    0     2 ]
[ 1   0  155     0 ]
[ 0   2    0  -308 ].
```

Its determinant is `-618^2`. Its rational quaternion algebra has Hilbert
symbol `(2,309/4)`, ramifies exactly at `2` and `3`, and has discriminant
`D=6`. The order reduced discriminant is `618=6*103`; the local
squarefree-discriminant classification therefore makes it maximal at the
ramified primes and the Borel/Eichler order of level `103` at the sole extra
split prime. Thus the literal order, not merely a similar primitive order,
has exact Eichler pair `(D,N)=(6,103)`. There is no similarity gap.

## Full discriminant action and stable kernel

The checker constructs three integral Clifford normalizer elements:

| class | coordinates in `1,e0e1,e0e2,e1e2` | norm | action on `A_T` |
|---|---|---:|---:|
| `w_2` | `(-14,-3,2,1)` | 2 | 619 |
| `w_6` | `(-10,-15,8,-6)` | 6 | 1031 |
| `w_618` | `(-4,-3,8,6)` | 618 | 1235 = -1 |

Conjugation preserves the literal order and gives integral determinant-one
isometries of `T`. Products give the full action

```text
w_1     1          w_103   205
w_2     619        w_206   823
w_3     413        w_309   617
w_6     1031       w_618   1235.
```

This is all of `O(A_T)`. The norm-one group of the literal order is stable.
The squarefree Eichler normalizer theorem says that the displayed
Atkin--Lehner cosets are all remaining projective components.

Central inversion also acts by `-1` on `A_T` but trivially on the period
line. Hence the projective class of `(-1)w_618` is stable and acts on the
period domain as `w_618`. No other nontrivial Atkin--Lehner class has action
`+/-1`. The exact marked curve is therefore `C_1236`, not the coarse
`X_0^6(103)` and not the full Atkin--Lehner quotient.

## Signatures and quotient maps

Ogg's formulas give

```text
X_0^6(103): genus 17, cusps 0, e2=0, e3=4, mu=208.
```

The fixed-point counts are

```text
involution   w2  w3  w6  w103  w206  w309  w618
fixed points  0   4   4     0     0    12    12.
```

Riemann--Hurwitz gives the useful tower:

| quotient of `X_0^6(103)` | degree | genus |
|---|---:|---:|
| `<w_618>` (exact marked curve) | 2 | 6 |
| `<w_2,w_309>` | 4 | 2 |
| `<w_3,w_206>` | 4 | 3 |
| `<w_6,w_103>` | 4 | 3 |
| full Atkin--Lehner group | 8 | 1 |

The exact marked signature is

```text
genus 6, cusps 0, e2=12, e3=2, mu=104.
```

Its first useful arithmetic map is

```text
C_1236 -> X_0^6(103)/<w_2,w_309>,     degree 2,
```

with six geometric branch points.

The fixed-point calculation identifies the branch locus exactly at the CM
level: all twelve `w_309`-fixed points upstairs have order discriminant
`-1236`, and `w_618` pairs them into the six branch points on the genus-two
quotient. What remains unknown is their coordinate divisor on the displayed
model of `B`, equivalently the resulting squareclass in `QQ(B)`.

## The low-genus quotient

Padurariu and Saia determine the exact model

```text
B: y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9.
```

There are exactly fourteen rational points:

```text
(0,+/-3),              (+/-1,+/-48),
(+/-1/3,+/-8/3),       (+/-1/5,+/-312/125).
```

The quotient by the involution induced by `w_3` is

```text
E: Y^2 + X*Y = X^3 - 185*X + 1401,       Cremona 618f1,
```

under

```text
X = 54*x^2 + 4,       Y = 9*y - 27*x^2 - 2.
```

The curve `E` has rank one, trivial torsion, and generator `G=(10,-29)`.
The displayed points map to `+/-G`, `+/-3G`, `+/-4G`, and `+/-10G`.
This is useful structure, but it does not lift through the nontrivial
degree-two marking cover automatically.

### Complete rational-point certificate for `B`

The completeness statement is an exact quadratic-Chabauty plus
Mordell--Weil-sieve certificate, not a bounded point search. The two elliptic
quotients are `618f1` and `618e1`; both have certified rank one, trivial
torsion, and saturated generators. At the good ordinary primes `7` and `11`,
the pinned bielliptic quadratic-Chabauty implementation gives `621` compatible
height classes. Restoring all four images under `x -> -x` and `y -> -y`, then
combining coefficients modulo `7^4` and `11^4`, leaves `231760` mock
Mordell--Weil cosets across `405` height classes.

Exact finite-field reductions eliminate every mock coset. The primes that
strictly reduce the candidate set are

```text
31, 37, 193, 199, 227, 449, 503, 743,
1093, 1427, 1733, 2647, 3347, 3539, 5273, 6599.
```

The final eight cosets are eliminated at `6599`. The checker asserts at least
four digits of coefficient precision, restores the automorphism orbits using
the exact actions

```text
x -> -x: (n1,n2) -> (n1,-n2),
y -> -y: (n1,n2) -> (-6-n1,-n2),
```

At each Chabauty prime it also asserts that the recognized rational roots are
exactly the three known nonzero-`x` orbits. The exceptional `x=0` orbit
`(0,+/-3)` is checked directly in both elliptic quotients, and the nonsquare
leading coefficient excludes rational points at infinity. The sieve enumerates
all reductions of `B(F_q)`, including points at infinity.
Consequently

```text
B(QQ) = {(0,+/-3), (+/-1,+/-48),
         (+/-1/3,+/-8/3), (+/-1/5,+/-312/125)}.
```

The quadratic-Chabauty code is pinned to Balakrishnan's
[`QC_bielliptic`](https://github.com/jbalakrishnan/QC_bielliptic) commit
`84af22e9cd1244c3d44e3c083073b44b8d728159`. Two representation-only changes
adapt its formal point-at-infinity tests to Sage 10.9; both source hashes and
the patch boundary are recorded in the generated certificate. The underlying
method is Bianchi--Padurariu,
[*Rational points on rank 2 genus 2 bielliptic curves in the LMFDB*](https://doi.org/10.1090/conm/796/16003).

The model is from [Padurariu--Saia, *Shimura curve Atkin--Lehner quotients
of genus at most two*](https://arxiv.org/abs/2509.25368), with machine data
in [`GenusAtMost2`](https://github.com/fsaia/GenusAtMost2) at commit
`6cc368fe37aa67187783118f18d149b2b1fd6230`.

## Rational-point separation

There are no cusps. For the CM order of discriminant `-3`, the residue-field
inputs are

```text
h(-3)=1,       D(R)=2,       N*(R)=103,
m_R=3,         Q=206=D(R)N*(R),       m=618.
```

The auxiliary quaternion algebra `(-3,206)` ramifies exactly at `2,3`.
Gonzalez--Rotger's formula makes every image of these points on the `w_618`
quotient rational. The four order-three points upstairs form two Fricke
pairs, producing exactly two rational CM points on `C_1236`.

The other two `w_3`-fixed phenomena are also separable. The four points of
order discriminant `-24` fixed by `w_6` upstairs have

```text
h(-24)=2,       D(R)=1,       N*(R)=103,
m_R=6,          m/m_R=103.
```

Corollary 5.14 gives a quadratic residue field on `C_1236`, not `QQ`: the
four top points become one quadratic closed point, or two conjugate geometric
points, on the marked curve. Thus they supply no rational marked point.

On the genus-two model, the involution induced by `w_3` is `x -> -x` and
its rational fixed locus is exactly `(0,+/-3)`. Up to interchanging those two
signs, one is the image of the pair of rational discriminant `-3` points on
`C_1236`; the other is the image of the conjugate discriminant `-24` pair.
Consequently both fixed fibers are understood. The remaining displayed
points, beginning with `x=+/-1/3`, `+/-1`, and `+/-1/5`, are the first
non-fixed fibers whose lift squareclasses must be computed.

| point class on `C_1236` | result |
|---|---|
| rational cusps | exactly 0 |
| certified rational CM points | exactly 2 of discriminant `-3` |
| certified nonrational CM points | one quadratic closed point of discriminant `-24` |
| other rational CM points | not determined |
| rational non-CM points | not determined |
| certified desired rank-19 markings | 0 |
| certified higher-Picard specializations | 2 |

The two known points add an algebraic class and have geometric Picard rank
20. They are precisely the higher-Picard points that the arithmetic-first
gate must not promote.

The residue-field input is Corollary 5.14 of Gonzalez--Rotger,
[*Non-elliptic Shimura curves of genus one*](https://doi.org/10.2969/jmsj/1179759530).
The signature calculation uses Ogg,
[*Real points on Shimura curves*](https://doi.org/10.1007/978-1-4757-9284-3_12).

## Jacquet--Langlands cover precheck

The exact weight-two newspace at classical level `618` has dimension `17`,
matching the genus of `X_0^6(103)`. At the ramified primes `2,3`, geometric
Atkin--Lehner signs are the negatives of the classical signs; at the
Eichler-level prime `103` they agree. The geometric `w_618`-invariant part
therefore has the exact isogeny decomposition

```text
Jac(C_1236) ~ 618a1 x 618b1 x 618c1 x 618d1 x 618e1 x 618f1.
```

The deck involution of `C_1236 -> B` is geometric `w_2=-W_2`. Its invariant
and anti-invariant parts are consequently

```text
Jac(B)                         ~ 618e1 x 618f1,
Prym(C_1236/B)                 ~ 618a1 x 618b1 x 618c1 x 618d1.
```

All six elliptic factors have Mordell--Weil rank one. Hence

```text
rank Jac(C_1236)(QQ) = 6 = genus(C_1236).
```

Classical Chabauty does not pass its strict rank bound. On the other hand,
the six pairwise nonisogenous factors give `rho(Jac(C_1236)) >= 6`, so the
necessary quadratic-Chabauty dimension inequality passes:

```text
6 = rank J(QQ) < genus(C_1236) + rho(J) - 1 >= 11.
```

This does not decide `C_1236(QQ)`: quadratic Chabauty still needs an explicit
model of the degree-two cover. The complete rational-point input on `B` is now
available. The factor accounting does identify
the four-dimensional Prym carrying the missing marking character and rules
out a classical rank-less-than-genus shortcut.

The factor accounting is replayed from modular symbols in the checker. Its
interpretation uses the classical Jacquet--Langlands correspondence and the
ramified-place Atkin--Lehner sign normalization.

## Exact next task

The next computation is narrow, not a broad rank or rational-point search:

1. construct `C_1236 -> B` explicitly;
2. compute its squareclass in `QQ(B)^*/QQ(B)^{*2}`;
3. evaluate the twelve non-fixed rational fibers in the now-complete set
   `B(QQ)`;
4. distinguish rational non-CM lifts from CM points or points realizing a
   saturated overlattice.

One non-CM rational lift is a positive certificate. A complete proof that
only the two CM points lift is a negative certificate. Points on `B`,
including the complete `+/-G`, `+/-4G`, and `+/-10G` fibers, remain
quotient-level evidence until their fibers are evaluated.

## Replay and independence

The checker is
[`scripts/certify_det1236_marked_shimura_curve.sage`](scripts/certify_det1236_marked_shimura_curve.sage),
and the generated certificate is
[`../artifacts/generated-results/elkies-k3-det1236-marked-shimura-curve-v1.json`](../artifacts/generated-results/elkies-k3-det1236-marked-shimura-curve-v1.json).
The quotient rational-point checker and certificate are
[`scripts/certify_det1236_genus2_rational_points.sage`](scripts/certify_det1236_genus2_rational_points.sage)
and
[`../artifacts/generated-results/elkies-k3-det1236-genus2-rational-points-v1.json`](../artifacts/generated-results/elkies-k3-det1236-genus2-rational-points-v1.json).

```bash
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --fresh
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage --check
```

The replay reads only the lattice catalogue and its transcendental-arithmetic
ledger. It does not read curve 356, curve 385, or any frozen prospective
experiment artifact. It performs no K3 equation work and no broad rank
search.
