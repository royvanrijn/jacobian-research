# Determinant1092: the literal marking and the rational CM locus

The [curve302 lattice lead](CURVE302_DET1092_RECONSTRUCTION_2026-09-07.md)
passes a new integral marking audit. Its full projective stable period curve
is **X(546)/<w546>** under the named ternary-spin/marked-K3 correspondence.
An independent algebra implementation verifies the finite arithmetic. There
is no additional marking cover in this identification.

A separate argument determines the rational CM locus abstractly: exactly
sixteen points, eight each of discriminants **-67 and -163**. Their coordinates
on the published elliptic model remain unknown. Neither assertion constructs
a K3 parent, a rational generic section basis, or a new elliptic fibre.

Authorities: `EC-DET1092-LITERAL-MARKING-20260907` and
`EC-DET1092-RATIONAL-CM-LOCUS-20260907` in `MATH_STATUS.json`.

## Integral arithmetic before equation work

Use the literal ternary Gram, with no similarity rescaling:

```
T = [[-2,1,0],[1,2,2],[0,2,220]].
```

In its even Clifford algebra use the integral basis
`1, e0*e1, e0*e2, e1*e2`, where `ei²=Tii/2` and
`ei*ej+ej*ei=Tij`. Exact multiplication closes this order. The reduced
trace pairing is

```
[[2,1,0,2],[1,3,2,0],[0,2,220,-110],[2,0,-110,-216]].
```

Its determinant is `-546²`. The algebra has Hilbert presentation
`(5,546/5)`, equivalent to the earlier `(5,2184/5)`, and discriminant546.
Equality of order and algebra reduced discriminants proves maximality;
the Eichler level is1.

The discriminant group is cyclic of order1092, with quadratic value
`5/2184` on the recorded generator. All sixteen orthogonal units are
enumerated exactly. A fixed search over `(b,c,d) in [-24,24]^3`, solving
the scalar coordinate by an integer square test, finds the following
normalizers after9085 triples:

| Norm / Atkin–Lehner label | Integral even-order coordinates | Action on A_T |
|---:|---|---:|
|2|(228,-24,-23,-5)|547|
|3|(227,-23,-19,11)|365|
|7|(173,-23,-24,-11)|937|
|13|(178,-21,14,-5)|337|

Their conjugations preserve both the integral order and T. The certificate
includes every matrix and verifies determinant one, the T isometry equation,
and the action on the dual generator. Products represent all sixteen
normalizer classes and exhaust O(A_T). Only labels1 and546 act as ±1.
Projectively, these give the full Fricke quotient.

The maximal-order normalizer theorem and the ternary-spin period
correspondence are external inputs; the finite matrices alone are not a
proof of those theorems. This identification also agrees directly with
[Elkies's L_N moduli description, section2](https://arxiv.org/pdf/0802.1301).
The [published elliptic model](https://web.mat.upc.edu/victor.rotger/docs/Tesi.pdf)
is `y²+xy+y=x³-137x+380`, Cremona546c2. Its infinite rational point set was
already certified in the reconstruction note.

## Complete CM discriminants, without coordinates

Let R be the CM order, h its class number, and D(R) the product of inert
primes among2,3,7,13. A rational point downstairs lifts to degree at most2.
González–Rotger's Theorem5.12 gives upstairs degree2h when D(R)>1 and h
otherwise. Thus the first case forces h=1. In the second case all four
primes ramify; quadratic genus theory makes the field class number divisible
by8, hence h>=8, excluding degree at most2.

Of the complete thirteen class-number-one orders, only discriminants-67
and-163 pass nonsplitting at all four primes. All four are inert for each
survivor. Proposition5.6 gives sixteen CM points upstairs per discriminant;
Corollary5.14 makes their full-Fricke images rational. The involution is free
there, giving eight images per discriminant, sixteen in total.
See [González–Rotger, appendix5](https://web.mat.upc.edu/victor.rotger/docs/ShimuraGenusOne.pdf).

Completeness uses the class-number-one theorem, including nonmaximal
orders; checking small discriminants alone would not suffice. The
[Sage reference list](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/cm.html)
records these thirteen orders. The checker computes their primitive reduced
forms and splitting symbols on the underlying fields, so an order conductor
does not hide a field-splitting obstruction.

## What construction still needs

The seventeen multiples `0P,...,16P` of `P=(-9,-26)` are distinct rational
points, so at least one is non-CM under any fixed moduli isomorphism. This
does **not** identify which one. Infinite order on the elliptic model does
not imply a non-CM K3 specialization.

The published isomorphism class does not by itself provide a coordinate
map to a K3 family or a positioned CM divisor. Translations of the elliptic
model preserve its isomorphism class while changing the coordinates of that
divisor. The next useful input is an explicit moduli map/CM coordinate
identification, followed by an individually certified non-CM rational point
and the rational marked K3 construction. The foundry's individual-point
gate remains open. No rootless-frame or coefficient search is authorized
by the finite-list existence statement alone.

## Replay and cost

The [normalizer certificate](../../artifacts/generated-results/elliptic-curves/det1092_marking_v1.json)
has an [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_marking_independent_v1.json)
that uses explicit Clifford word reduction instead of Sage's Clifford
implementation or repository algebra code. It also reconstructs the Hilbert
presentation. The [CM arithmetic certificate](../../artifacts/generated-results/elliptic-curves/det1092_cm_locus_v1.json)
relies on the named completeness and residue-field theorems; a separate
independent CM arithmetic implementation is not claimed.

```
sage -python elliptic-curves/cas/verify_det1092_marking.sage \
  --input artifacts/generated-results/elliptic-curves/det1092_marking_v1.json \
  --output NEW_REPLAY.json
python3 elliptic-curves/cas/report_det1092_marking.py --check
```

Three sequential stages used **2.031180749 supervised seconds** in total,
with fixed wall caps120/60/60seconds and2GiB each. Protocols, source hashes,
logs and terminal ledgers are retained under the three local paths listed
in the [bound report](../../artifacts/generated-results/elliptic-curves/det1092_marking_followup_v1.json).
There were no point-search boxes and no high-rank inventory changes.
