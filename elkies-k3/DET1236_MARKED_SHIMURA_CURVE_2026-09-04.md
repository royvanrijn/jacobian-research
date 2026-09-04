# Determinant 1236: exact marked Shimura curve and remaining lift obstruction

<!-- status-consumer: EC-K3-DET1236-GENUS2-RATIONAL-POINTS 5a3c84eb9f7f0604 -->
<!-- status-consumer: EC-K3-DET1236-MARKED-SHIMURA-CURVE e482668e1208f764 -->
<!-- status-consumer: EC-K3-DET1236-V4-LOCAL-CONSISTENCY deb09ed7326145f9 -->
<!-- status-consumer: EC-K3-DET1236-CANDIDATE-DOUBLE-COVER fecec75b4ff1e8e1 -->
<!-- status-consumer: EC-K3-DET1236-RATIONAL-CM-LOCUS bd6ab0e86ca70ab2 -->

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

This is no longer an unidentified marked curve. It has no cusps and has
exactly ten rational CM points: two of discriminant `-3` and four each of
discriminants `-43` and `-67`. They are Picard-rank-20 specializations, so
none realizes the requested saturated rank-19 Neron--Severi lattice. No
rational non-CM point on `C_1236` is currently certified.

The rational points on the displayed genus-two quotient are complete: there
are exactly fourteen.  An exact degree-two cover candidate has now been
constructed and all twelve non-fixed fibres evaluated.  Its displayed twist
would give eight rational lifts; its `-3` twist gives none.  A separate exact
CM residue-field calculation proves that the marked curve has exactly ten
rational CM points: two of discriminant `-3` and four each of discriminants
`-43` and `-67`.  Thus all eight conditional non-fixed lifts are CM, and the
no-nonfixed-lift twist cannot be the actual descent.  The precise remaining
problem is only to identify the candidate's cubic branch orbit with the
order-discriminant-`-1236` Shimura branch orbit in characteristic zero.  If
that succeeds, the row is `ARITHMETICALLY_EXCLUDED`, not realizable. Until
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
| certified rational CM points | exactly 10: 2 of discriminant `-3`, 4 of `-43`, and 4 of `-67` |
| certified nonrational CM points | one quadratic closed point of discriminant `-24` |
| other rational CM points | none |
| rational non-CM points | not determined until the candidate cover is identified |
| certified desired rank-19 markings | 0 |
| certified higher-Picard specializations | 10 |

Every CM point adds an algebraic class and has geometric Picard rank 20. None
of the ten certified rational CM points can pass the exact-rank-19 gate.

### Complete rational CM locus on the marked curve

Gonzalez--Rotger Corollary 5.14 applies directly to
`C_1236=X_0^6(103)/<w_618>`.  Apart from the exceptional two-involution case,
the residue degree is at least the order class number `h(R)`; in that case it
is at least `h(R)/2`.  Therefore a rational CM image forces `h(R)<=2`.
Checking the complete class-number-one and class-number-two order lists with
the exact local optimal-embedding factors leaves precisely

```text
order discriminant   top CM points   rational points on C_1236
-3                         4                       2
-43                        8                       4
-67                        8                       4
```

For `-3`, `(D(R),N*(R),m_R,m/m_R)=(2,103,3,206)`; for `-43` and
`-67` it is `(6,103,1,618)`.  In each case Corollary 5.14 gives residue
field `QQ`.  The `w_618` action is free on the `-43` and `-67` loci and pairs
the top-curve points.  No class-number-two order reaches the only residue
configuration that could be rational.  Hence the displayed ten points are
the complete rational CM locus, not merely known examples.

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

This abstract isogeny accounting does not by itself identify the explicit
quadratic characters in a published quotient model.  The corrected local
audit below is compatible with the `ab|cd` character assignment and supplies
a very rigid characteristic-zero candidate, but the final modular
identification of its CM branch orbit is still required.

## Double-cover reconstruction audit

The current Atkin--Lehner factor assignment suggests a `V_4` formulation of
the lift gate.  Put

```text
D_3 = C_1236/<w_3>.
```

Under that assignment, `D_3 -> E` is a degree-two genus-three cover, its
Prym factors are `618c1 x 618d1`, and `C_1236` is the normalization of
`B x_E D_3`.  Thus a rational point of `B` would lift to `C_1236` exactly
when its image on `E` lifts to `D_3`.  The expected rational branch point of
`D_3 -> E` is one of `+3G,-3G`; the other three geometric branch points are
the images of the order-discriminant `-1236` branch divisor.

This gives a finite local reconstruction.  On

```text
E: Y^2+X*Y=X^3-185*X+1401
```

a squareclass with four branch points has a representative

```text
b = A(X) + B(X)*Y in L(6*O),

Norm(b) = A(X)^2-X*A(X)*B(X)
          -(X^3-185*X+1401)*B(X)^2.
```

For the already certified quotient map

```text
B -> E,   X=54*x^2+4,   Y=9*y-27*x^2-2,
```

the desired squareclass on `B` is therefore literally

```text
b_B(x,y) = b(54*x^2+4, 9*y-27*x^2-2)
            in QQ(B)^*/QQ(B)^{*2}.
```

In particular `(x,y)` and `(-x,y)` have the same value.  The twelve
non-fixed rational evaluations consequently reduce, without loss, to six
evaluations at the images `+/-G`, `+/-4G`, and `+/-10G` on `E`.

After prescribing the rational branch point, the remaining even divisor is
represented by a double zero `Q`.  Since `E(QQ)=ZZ*G`, its reduction lies in
`<G>` at every good prime.  The exceptional case `Q=O` is represented in
`L(4*O)`.  If one of the other branch points specializes to `O`, the
corresponding odd-pole charts are `L(5*O)` and, when also `Q=O`, `L(3*O)`.
Enumerating all four projective charts over `F_p` and imposing the degree-one,
degree-two, and degree-three Frobenius power traces of `618c1 x 618d1` gives
the fail-closed local screen.

The known cover `B -> E` has squareclass

```text
h = (X-4)/54
```

and Prym `618e1`.  At the good prime `p=5`, every genus-three quadratic
extension of `F_5(E)` has a squareclass representative in `L(6*O)`.  The
exhaustive audit enumerates all

```text
#P(L(6*O))(F_5) = (5^6-1)/(5-1) = 3906
```

projective classes.  Exactly `1040` have branch degree four, and `304` have
branch degree four for both `b` and `h*b`.  For each class it computes over
`F_25` the two character sums belonging to `b` and `h*b`.  The degree-two
Frobenius targets forced by the asserted Pryms are

```text
618a1 x 618b1 : 15,
618c1 x 618d1 : 11.
```

There are `24` compatible classes.  At `p=7`, the corresponding corrected
counts are `3744`, `848`, and `56`.  The earlier zero counts were caused by
evaluating `chi(b(P))=0` at every zero.  That is wrong at the auxiliary
double zero `2Q`: its contribution is the quadratic character of the leading
local unit after the square uniformizer is removed.  Literal-square unit
tests now guard this convention.  The old upstream-inconsistency conclusion
is withdrawn.

The corrected degree-one-through-three `V_4` screen leaves four classes at
`p=5,7` and two classes at each of `p=11,13,17,19`, one for each choice of
rational branch point.  This is local compatibility, not yet a global cover.

## Exact characteristic-zero candidate

Put

```text
K = QQ[a]/(a^3-a^2+4*a+12),    disc(K)=-1236.
```

This is not merely a cubic field with a suggestive discriminant.  Exact
class-field arithmetic gives class group `C6 x C2` for
`QQ(sqrt(-1236))`; its Hilbert class field has relative degree `12` and
absolute degree `24`.  Its three embedded cubic subfields all have
discriminant `-1236` and are `QQ`-isomorphic to `K`.  Thus `K` is the unique
cubic residue-field isomorphism class available to the CM branch image.
What this does not yet identify is the particular point of `E(K)`.

On `E/K`, take

```text
G=(10,-29),
A=(-3*a^2-9*a-20, 15*a^2-57*a-155).
```

Exact addition in the splitting field gives `Tr(A)=-16G`.  The unique small
global class matching the local survivors at `p=11,13,17,19` is

```text
P = 9G+2A
  = (18*a^2-48*a+130, 270*a^2-672*a+2059),
Tr(P) = -5G.
```

With rational branch point `-3G` and auxiliary double zero `Q=4G`, the sum
`-3G+Tr(P)+2Q` is zero.  Exact linear algebra in `L(6O)` gives

```text
b0 = -X^3-328*X^2-2772*X+66512+(32*X+1600)*Y,

div(b0) = (-3G)+Orbit(P)+2*(4G)-6O,

Norm(b0) = (X-4)(X-58)^2
           (X^3-216X^2-6924X-62224).
```

The last cubic is the minimal polynomial of `x(P)` and has discriminant
`-1236*24192^2`.  The twist selected by the observed Jacquet--Langlands trace
signs is

```text
D_candidate: z^2 = 2*b0.
```

This is much stronger than a local survivor.  For each of
`p=5,11,13,17,19,23,29,31`, exact character sums over `F_(p^n)` for
`n=1,2,3` agree with `618c1 x 618d1`; multiplying by
`h=(X-4)/54` agrees with `618a1 x 618b1`.  The smooth canonical plane quartic
of the candidate is

```text
87616 U^4 - 43784 U^3 V + 7944 U^2 V^2 - 608 U V^3 + 16 V^4
+ 268 U^2 W^2 - 118 U V W^2 + 10 V^2 W^2 + W^4 = 0,

[U:V:W]=[X-58:Y-403:z].
```

These finitely many Euler-factor identities are exact evidence, but they are
not a characteristic-zero isomorphism or modular-function proof.

Pulling `2*b0` back to `B` gives the particularly small squareclass

```text
C_candidate:
  Z^2 = 3*(-81*x^6-534*x^4+8*x^2*y-177*x^2+8*y+24).
```

Indeed `2*b0(X(x),Y(x,y))=3888*q=3*36^2*q`.  Its exact normalization-fibre
norm is

```text
Norm_B(q) = 81*x^2*(x^2-1)^2
            *(81*x^6-306*x^4-239*x^2-48).
```

The last factor is squarefree of degree six.  Thus the candidate cover of the
genus-two curve has six geometric branch points and normalization genus six,
as required.  Its exact normalization-fibre evaluations are:

| points on `B` | image on `E` | squareclass | rational fibre? |
|---|---:|---:|---|
| `(+-1,48)` | `4G` | `1` (even zero, leading unit `36`) | yes |
| `(+-1,-48)` | `-4G` | `-2` | no |
| `(+-1/3,8/3)` | `-G` | `1` | yes |
| `(+-1/3,-8/3)` | `G` | `-11` | no |
| `(+-1/5,312/125)` | `-10G` | `37` | no |
| `(+-1/5,-312/125)` | `10G` | `-383` | no |

Thus the displayed candidate has four non-fixed rational fibres, hence eight
rational normalization points.  The `-3` constant twist has no rational
fibre at any of the twelve non-fixed points.  The two twists exchange the
fixed CM behaviour: for the displayed candidate `(0,3)` has value `144`,
whereas at `(0,-3)` the normalized leading unit has squareclass `-3`; the
`-3` twist swaps the rational and quadratic fixed fibres.  Because the
published model identifies those two CM image classes only up to `y`-sign,
the fixed-point controls alone do not choose between these outcomes.

The complete rational CM calculation does choose the outcome conditionally
on the branch divisor.  The actual marked curve already has four rational
CM points of discriminant `-43` and four of discriminant `-67`, besides its
two discriminant-`-3` points.  Therefore the no-nonfixed-lift `-3` twist
cannot be the marked descent.  For the displayed twist the total rational
point count is `2+8=10`, exactly the complete rational CM count.  Consequently
all eight conditional lifts are CM; none realizes the desired rank-19
marking.

The status remains `UNRESOLVED_FOR_EXPLICIT_REASON`.  One exact modular step
is missing: the cubic CM residue-field class is now matched, but one must
prove that the particular point `P=9G+2A` is the image of the
order-discriminant-`-1236` CM orbit on this specified genus-one quotient. The
rational CM count then forces its rational descent to have the displayed
constant twist. This can be
closed by an explicit characteristic-zero model/map for the genus-three
quotient or by a Borcherds--Schofer CM-value calculation on the genus-one
base. After that identification, the complete quotient-point and CM-locus
certificates give `ARITHMETICALLY_EXCLUDED` immediately.

## Exact next task

The next computation is narrow, not a broad rank or rational-point search:

1. identify the displayed canonical plane quartic with
   `X_0^6(103)/<w_3,w_618>` over `QQ`, or compute the corresponding CM
   modular value, thereby certifying the cubic orbit and constant twist;
2. invoke the certified complete rational CM locus to identify all eight
   non-fixed lifts as the discriminant-`-43` and `-67` points;
3. issue `ARITHMETICALLY_EXCLUDED` for determinant `1236`.

The squareclass and all twelve evaluations are exact for the candidate, and
the CM/non-CM separation is complete conditionally on its identification.
Only the marked Shimura branch-orbit identification remains open.

## Replay and independence

The checker is
[`scripts/certify_det1236_marked_shimura_curve.sage`](scripts/certify_det1236_marked_shimura_curve.sage),
and the generated certificate is
[`../artifacts/generated-results/elkies-k3-det1236-marked-shimura-curve-v1.json`](../artifacts/generated-results/elkies-k3-det1236-marked-shimura-curve-v1.json).
The quotient rational-point checker and certificate are
[`scripts/certify_det1236_genus2_rational_points.sage`](scripts/certify_det1236_genus2_rational_points.sage)
and
[`../artifacts/generated-results/elkies-k3-det1236-genus2-rational-points-v1.json`](../artifacts/generated-results/elkies-k3-det1236-genus2-rational-points-v1.json).
The fail-closed local squareclass screen is
[`scripts/explore_det1236_double_cover_local_gate.sage`](scripts/explore_det1236_double_cover_local_gate.sage).
The exhaustive `p=5` `V_4` consistency audit and its generated certificate
are
[`scripts/audit_det1236_v4_local_consistency.sage`](scripts/audit_det1236_v4_local_consistency.sage)
and
[`../artifacts/generated-results/elkies-k3-det1236-v4-local-consistency-v1.json`](../artifacts/generated-results/elkies-k3-det1236-v4-local-consistency-v1.json),
with the independent `p=7` output in
[`../artifacts/generated-results/elkies-k3-det1236-v4-local-consistency-p7-v1.json`](../artifacts/generated-results/elkies-k3-det1236-v4-local-consistency-p7-v1.json).
The exact candidate checker and fail-closed certificate are
[`scripts/certify_det1236_candidate_double_cover.sage`](scripts/certify_det1236_candidate_double_cover.sage)
and
[`../artifacts/generated-results/elkies-k3-det1236-candidate-double-cover-v1.json`](../artifacts/generated-results/elkies-k3-det1236-candidate-double-cover-v1.json).
The complete rational CM-locus checker and certificate are
[`scripts/certify_det1236_rational_cm_locus.sage`](scripts/certify_det1236_rational_cm_locus.sage)
and
[`../artifacts/generated-results/elkies-k3-det1236-rational-cm-locus-v1.json`](../artifacts/generated-results/elkies-k3-det1236-rational-cm-locus-v1.json).

```bash
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --fresh
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --check
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage --check
sage elkies-k3/scripts/explore_det1236_double_cover_local_gate.sage 5 7
sage elkies-k3/scripts/audit_det1236_v4_local_consistency.sage write
sage elkies-k3/scripts/audit_det1236_v4_local_consistency.sage 7 write
sage elkies-k3/scripts/search_det1236_cm_orbit_cover.sage 10
sage elkies-k3/scripts/certify_det1236_rational_cm_locus.sage check
sage elkies-k3/scripts/certify_det1236_candidate_double_cover.sage
sage elkies-k3/scripts/certify_det1236_candidate_double_cover.sage check
```

The replay reads only the lattice catalogue and its transcendental-arithmetic
ledger. It does not read curve 356, curve 385, or any frozen prospective
experiment artifact. It performs no K3 equation work and no broad rank
search.
