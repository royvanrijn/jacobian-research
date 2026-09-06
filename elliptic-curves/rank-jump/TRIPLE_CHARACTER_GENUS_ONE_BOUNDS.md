# A genus-one three-point construction has only two unresolved mixed directions

The frozen three-generic-point test gives the same bounds on the retained
MW16-04 high/zero-gain pair:
\[
\boxed{
4\le \operatorname{rank}E_u/\overline{\mathbb Q}(C)\le6,\qquad
3\le \operatorname{rank}E_u/\mathbb Q(C)\le5.
}
\]
Two of the four mixed character spaces have **exact geometric rank zero**
on each anchor. The other two each have rank zero or one. Both genus-one
bases have infinitely many rational points.

This supplies a viable small simultaneous construction and a certified
ceiling on its generic size. It does **not** explain the high fibre's +9
quotient: the three guaranteed directions specialize to already marked
generic points at the anchor, and no remaining mixed section or exceptional
specialization image has been constructed.

| Anchor, same original MW16-04 family | Original generic rank | Certified witness quotient | New construction's arithmetic generic rank |
|---|---:|---:|---|
| `-1647/91` | 16 | 9 | 3–5 |
| `-2177/2397` | 16 | observed 0, censored | 3–5 |

The last column is for a **different family through the anchor**. It neither
replaces the original generic rank 16 nor bounds either anchor's rank.
The observed-zero row remains a lower-bound observation.

## 1. Frozen construction and scope

The [protocol](TRIPLE_CHARACTER_PROTOCOL.json) fixes the same matched pair
as [the two-point test](MIXED_CHARACTER_PAIRED_BOUNDS.md), taking generic
indices 0, 1 and 2. The
[input projection](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_inputs_v1.json)
contains only the short equations and these three generic points.
No exceptional coordinates enter the construction or its point counts.
The high/low selection is retrospective, so this is not validation of a
prospective rank predictor.

Write the three points as \((a_i,y_i)\) on
\(E_0:y^2=f(x)=x^3+Ax+B\), and put
\[
g_i=1-a_i u,\quad
C:\ z_i^2=g_i,\quad
d_I=\prod_{i\in I}g_i.
\]
All abscissas are nonzero and distinct; `B` is nonsquare on both anchors.
Use the same fixed-cubic pencil
\[
E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]

Each `g_i` gives the explicit rank-one arithmetic twist from
[the linear-twist theorem](LINEAR_TWIST_SOLUBLE_BLOCKS.md).
The trivial character contributes geometric rank one and arithmetic
rank zero. Thus only masks `3,5,6,7`, respectively
`g0*g1, g0*g2, g1*g2, g0*g1*g2`, remain to be determined.
The guaranteed three arithmetic directions occupy distinct character
spaces and are independent over \(\mathbb Q(C)\).

## 2. A positive-rank genus-one base, before exceptional points

There are four branch points: the three private zeroes and the shared
point at infinity. Their inertia vectors are
\(100,010,001,111\). The cover is geometrically connected of degree eight,
and
\[
2g(C)-2=-16+4\cdot4=0.
\]
For `a=a0,b=a1,c=a2`, its smooth projective model is the intersection
\[
b z_0^2-a z_1^2=(b-a)h^2,\qquad
c z_0^2-a z_2^2=(c-a)h^2.
\]
The point \((1:1:1:1)\) is rational and smooth.

The norm map
\[
C\longrightarrow C_N:\quad
v^2=(1-au)(1-bu)(1-cu),\qquad
v=z_0z_1z_2/h^3
\]
is an unramified map of degree four. Both covers of the `u`-line have
ramification index two at the same four branch points. After choosing
rational origins, its geometric kernel is the even-sign `V4`, hence the
entire 2-torsion. Factoring through multiplication by two gives an
isomorphism \(C\simeq C_N\) over \(\mathbb Q\).

Set \(L=-abc\). The change \(X=Lu,\ Y=Lv\) identifies `C_N` with
\[
E_N:\quad
Y^2=X^3+(ab+ac+bc)X^2+abc(a+b+c)X+(abc)^2.
\]
The image of the rational origin is \(Q=(0,L)\). Its non-torsion proof
uses only the fixed generic-point coordinates:

| Anchor | Good-prime point orders | Hypothetical torsion order divides | Exact check |
|---|---|---:|---|
| high | `ord17(Q)=4, ord29(Q)=4` | 4 | `[4]Q != O` |
| low | `ord17(Q)=12, ord19(Q)=6` | 6 | `[6]Q != O` |

The bound follows from prime-to-residue-characteristic torsion injection
at both primes. The nonzero multiples are evaluated exactly and retained
in [the geometry certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_geometry_v1.json).
Thus both parameter curves have positive rational rank and infinitely many
rational parameters. This conclusion uses no rational point search.

The three square conditions therefore are simultaneously realizable on an
infinite rational parameter set. This is a genuine **solubility construction**
for the three known singleton sections. It is not a certificate for the
mixed characters.

## 3. All four remaining character surfaces are K3

Put \(D=1+Au^2+Bu^3\) and \(\delta=-4A^3-27B^2\).
For a conventional quadratic twist by `d=d_I`,
\[
\Delta=16\delta D^2d^6,\qquad
c_4=16d^2(A^2u^2-9Bu-3A).
\]
Every weight-two or weight-three mask has fibres
\[
3I_2+3I_0^*
\]
and trivial lattice rank 17. For weight two, two `I0*` fibres are finite
and the third is at infinity. For weight three, all three are finite and
infinity is smooth. The symbolic verifier checks these configurations,
including coprimality and infinity valuations. Shioda–Tate gives
\(\operatorname{rank}E_u^{d_I}=\rho_I-17\); the surface tools are standard
[Schütt–Shioda, sections 5–6](https://arxiv.org/pdf/0907.0298).

For the triple mask, write \(k=-a_0a_1a_2\). Its smooth infinity fibre is
\[
Y^2=X^3+2AkX^2+A^2k^2X-B^2k^3.
\]
Its character sum **must** be included in the surface count.
For example, at `p=17` on the high anchor, the finite-base sums for
degrees 1,2,3 are `28,-186,6355`; infinity contributes `3,25,-126`.
The correct complementary traces are `31,-161,6229`.
Treating infinity as a cusp would give the wrong Frobenius calculation.

As in the two-point test, if `r_q` counts roots of `f` over \(\mathbb F_q\),
the trivial lattice has trace \(q(5+4r_q)\). Removing it leaves a
five-dimensional complement. Three field counts determine its reciprocal
Frobenius polynomial; normalized cyclotomic factors give the reduction
Picard bound. An independent companion-matrix calculation checks all traces,
cyclotomic dimensions and Artin–Tate determinants.

## 4. Four exact zero-rank certificates

Only two new untwisted prime counts were needed: `p=11` on each anchor.
The four previous prime arrays are reused byte for byte. The protocol
fixes high-anchor primes `11,13,17` and low-anchor primes `11,17,19`,
with degrees at most three. A mask is omitted at a prime if its branch
points collide or its required coefficient units fail; no replacement
prime is added.

The following pairs of rank-18 reductions have incompatible
Néron–Severi discriminant squareclasses:

| Anchor | Mask / twist | First witness | Second witness | Characteristic-zero Picard rank | Mixed MW rank |
|---|---|---|---|---:|---:|
| high | `5: g0*g2` | `p=11: -5` | `p=17: -105` | 17 | 0 |
| high | `6: g1*g2` | `p=11: -5` | `p=13: -1` | 17 | 0 |
| low | `5: g0*g2` | `p=11: -5` | `p=17: -113` | 17 | 0 |
| low | `6: g1*g2` | `p=11: -5` | `p=17: -17` | 17 | 0 |

If the characteristic-zero Picard rank were 18, equal-rank specialization
would force equal discriminant squareclasses. The mismatch, together with
the rank-17 trivial lattice, proves rank 17. This is the
[van Luijk specialization argument](https://arxiv.org/pdf/math/0506416);
the discriminants are computed by Artin–Tate on the elliptic K3 reductions,
as explained in the preceding paired note.

Masks `3` and `7` each have a rank-18 good reduction, giving geometric
MW rank at most one. The frozen data do not close either exact rank.
For high-mask `7`, only `p=17` is eligible. For low-mask `7`, `p=17`
has reduction Picard rank 20 and `p=19` has rank 18, so these are not a
two-place rank-17 witness. All omitted and unresolved rows remain in
[the report](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_v1.json).

The four non-singleton intervals are therefore
\[
(\operatorname{rank}_3,\operatorname{rank}_5,
  \operatorname{rank}_6,\operatorname{rank}_7)
\in [0,1]\times\{0\}\times\{0\}\times[0,1]
\]
on **both** anchors. Adding the known characters proves the opening bounds.

## 5. Verification and replay

The [raw arrays](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_counts_v1.json)
bind all inputs and source versions. The
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_character_verification_v1.json)
checks 2,926 new base-parameter counts directly with NumPy, using no elliptic
point counter, and binds the 20,056 previously independently checked counts.
It recomputes every new character weight and every infinity contribution.
Sage/PARI supplies the original elliptic cardinalities; Sage exact linear
algebra independently verifies the reconstructed Frobenius polynomials.

The frozen bounds were not enlarged after seeing the results. No production
curve was generated, no new exceptional point was sought, and no search
policy or live output was changed.

From the repository root:

```sh
sage -python elliptic-curves/rank-jump/triple_character.py check
python3 elliptic-curves/rank-jump/verify_triple_character.py
sage -python elliptic-curves/rank-jump/verify_triple_character_geometry.py
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_triple_character.py
```

## 6. Mechanism ranking and what remains missing

1. **Solubility:** a positive-rank genus-one base can make three generic-point
   covers soluble together. This is now explicit on both controls. It does
   not distinguish their exceptional quotients.
2. **Incidence:** mixed-character blocks remain possible, but only two
   geometric directions are unresolved in this construction. Two other
   candidate directions are rigorously absent. A nonzero mixed section,
   its field of definition, and its specialization modulo the original
   generic subgroup are the remaining requirements.
3. **Weakened explanation:** low genus, an infinite rational base, and many
   available character labels do not by themselves explain a large jump.
   The first two properties hold on the low control too; half the mixed
   spaces vanish on both.
4. **Visibility:** this experiment measures no chart or point-search feature.
   Its generic-rank bounds cannot veto any individual production fibre.

Do not grow this construction by adding a fourth point merely to obtain
more labels: four independent linear characters with shared infinity have
genus five. The next useful work is to resolve an existing mixed section
and its specialization image, or to close the independent cup-product
solubility calculation on the retained class blocks. Either would address
a missing implication. An enlarged cover count alone would not.

Agent 1 can use these results as construction constraints: this particular
three-point base carries at most two mixed directions beyond its three
known singleton directions. It supplies no validated candidate score.
The unresolved target remains a condition, available before exceptional
points, that forces several new rational directions **modulo the original
generic subgroup**, not merely sections of a newly chosen family.
