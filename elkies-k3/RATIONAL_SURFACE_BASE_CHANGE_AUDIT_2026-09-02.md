# Rational-surface base changes behind the repeated-fibre sources

## Outcome

The two promoted repeated-fibre patterns both have an exact hidden quadratic
base-change explanation in the model and field where their equations are
currently certified.

| source model | deck involution | quotient coordinate | rational elliptic quotient | section result |
|---|---|---|---|---|
| rational `3I6+6I1` Golay-chart specialization | `t -> 1/t` | `u=t+1/t` | `I6+I3+3I1`, MW rank 1, torsion `Z/3` | invariant rank 1 plus anti-invariant rank 1 |
| marked NS0031 model 157 over `GF(7)`, `I2+2I8+6I1` | `t -> t/(t-1)` | `u=t^2/(t-1)` | `I8+4I1`, MW rank 1 | one exact invariant trace and one exact anti-invariant direction |

In both cases the degree-24 `j`-map is the composition of a degree-12
rational-surface `j`-map with the displayed quadratic quotient.  More
strongly, the short-Weierstrass coefficients themselves descend with weights
four and six.  This rules out a coincidental equality of `j`-maps or an
unaccounted quadratic twist.

No cubic construction explains a promoted marked source.  The complete
normalized `3I6` fibre charts do contain models with `C3` or `S3` `j`-symmetry,
but none of the `GF(7)` models carrying one of the 24 correct marked MW2 pairs
has an order-three stabilizer.  The three exact rational points found in the
bounded Golay parameter scans all have only a `C2` stabilizer.  For the
weighted `I2+2I8` support, no order-three Möbius permutation is possible.

The canonical machine certificate is
[`../artifacts/generated-results/elkies-k3-repeated-fibre-rational-base-change-audit-v1.json`](../artifacts/generated-results/elkies-k3-repeated-fibre-rational-base-change-audit-v1.json).
Its SHA-256 is
`7bb12308cf379b091b7c55f767df253690cad659284b3b9c51a4dbf55b8b178c`.
Replay it with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_repeated_fibre_rational_base_changes.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_repeated_fibre_rational_base_changes.sage --check
```

## Why the fibre profiles suggest these quotients

For a base change ramified to order `e` at a multiplicative fibre, `I_n`
pulls back to `I_(en)`.  Thus the observed Euler profiles have only one
quadratic rational-surface explanation compatible with the normalized
reducible supports:

```text
I6 + I3 + 3I1
  -- degree 2, branch at I3 and one smooth fibre -->
3I6 + 6I1,

I8 + 4I1
  -- degree 2, branch at one I1 and one smooth fibre -->
2I8 + I2 + 6I1.
```

Both quotient profiles have Euler number 12 and short-Weierstrass degree
bounds `(4,6)`, so they are rational elliptic surfaces.  Their root ranks are
seven, hence Shioda--Tate gives MW rank one.

This numerical match is only a candidate generator.  The certificate also
requires a support-preserving Möbius involution, exact `j`-invariance, full
discriminant covariance, coefficient-level descent, the quotient fibre
factorization, and section identities.

## Rational `3I6` model

Let

```text
E_t: y^2=x^3+A(t)x+B(t)
```

be the exact rational `s6=10` model in
[`../artifacts/generated-results/elkies-k3-golay-det720-3a5-source-qq-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-3a5-source-qq-v1.json).
Its three reducible supports are `0,1,infinity`, so their unweighted Möbius
stabilizer is `S3`.  Testing all six elements gives only

```text
j(t)=j(1/t).
```

The two order-three maps fail.  The stronger coefficient identities are

```text
t^8*A(1/t)=A(t),
t^12*B(1/t)=B(t),
t^24*Delta(1/t)=Delta(t).
```

With `u=t+1/t`, there are exact polynomials `a,b` such that

```text
A(t)=t^4*a(u),
B(t)=t^6*b(u).
```

They are

```text
a(u) = -27*u^4 - 648/5*u^3 - 57672/25*u^2
       + 184032/125*u + 3075408/625,

b(u) = 54*u^6 + 1944/5*u^5 + 36936/5*u^4
       + 295488/25*u^3 + 12900384/125*u^2
       - 1624841856/3125*u + 7131556224/15625.
```

After `x=t^2*X`, `y=t^3*Y`, the quotient is

```text
Y^2=X^3+a(u)X+b(u),

4*a^3+27*b^2
 = (557256278016/125)
   * (u-2)^3
   * (u^3 + 18/5*u^2 + 2268/25*u - 21384/125).
```

The missing discriminant degree is six at infinity.  Hence the quotient
fibres are exactly `I6+I3+3I1`.  The quadratic map is branched at `t=1,-1`,
with quotient values `u=2,-2`; the `I3` at `u=2` becomes the third `I6`, while
the other branch fibre is smooth.

### Sections

Use `T` for the exact 3-torsion section and `H` for the rational half of the
displayed section `Q`, as certified in
[`../artifacts/generated-results/elkies-k3-golay-det720-3a5-saturation-rejection-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-3a5-saturation-rejection-v1.json).
For the induced deck action

```text
sigma(x(t),y(t))=(t^4*x(1/t),t^6*y(1/t)),
```

the exact group-law identities are

```text
sigma(T)=T,
sigma(H)=H,
P+sigma(P)=2T,
sigma(P-T)=-(P-T).
```

Thus `T` and `H` descend, while `P-T` is a section on the quadratic twist of
the rational quotient.  Explicitly,

```text
Hbar = (3*u^2 + 36/5*u + 2988/25, 41472/25),
Tbar = (3*u^2 + 36/5*u + 108/25, (1728/5)*(u-2)).
```

The quotient has root lattice `A5+A2`, MW rank one, and torsion `Z/3`.
The pullback height of `Hbar` is one, so `Hbar` has height `1/2` and generates
the free quotient MW group.  The source free rank two is therefore explained
exactly as

```text
rank 1 invariant pullback + rank 1 anti-invariant quadratic-twist section.
```

This makes the hidden base change useful structurally, but it does not rescue
the model for the determinant-720 target: its full saturated NS determinant
is still `-20`, as already proved by the torsion and halving certificate.

## Marked NS0031 model over `GF(7)`

The exact positive source precursor is model 157 in
[`../artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json),
with its two marked MW2 pairs stored in
[`../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json).

The weighted supports `0:I2,1:I8,infinity:I8` have one possible nontrivial
Möbius automorphism:

```text
sigma(t)=t/(t-1).
```

It passes all exact gates:

```text
j(t)=j(t/(t-1)),
(t-1)^8*A(t/(t-1))=A(t),
(t-1)^12*B(t/(t-1))=B(t),
(t-1)^24*Delta(t/(t-1))=Delta(t).
```

With `u=t^2/(t-1)`, one has

```text
A(t)=(t-1)^4*a(u),
B(t)=(t-1)^6*b(u),

a(u)=4*u^4+3*u^3+6*u^2+4,
b(u)=5*u^6+3*u^5+3*u^3+6*u^2+6*u+2
```

over `GF(7)`.  The quotient discriminant is

```text
2*u*(u+5)*(u^2+3*u+6).
```

Its degree-four finite part and order eight at infinity give `I8+4I1`.
The branch points are `t=0,2`, with quotient values `u=0,4`; the `I1` at
`u=0` pulls back to the displayed `I2`, while the other branch fibre is
smooth.

For one displayed marked pair `P,Q`, exact finite-field group law gives

```text
P+sigma(P) != O,
sigma(P+sigma(P))=P+sigma(P),
sigma(3P-2Q)=-(3P-2Q).
```

Thus the pair already exposes a nonzero invariant trace and a nonzero
anti-invariant direction.  This is the expected rational-surface/twist split
behind an MW2 quadratic base change.  It is not yet a full MW decomposition:
the complete `GF(7)(t)` MW group is not certified, and no characteristic-zero
or rational NS0031 source equation has been constructed.

## Complete normalized-chart census

The symmetry is selective rather than forced by the repeated fibre type.

| fibre chart | prime | all squarefree models | nontrivial `j` stabilizer | breakdown |
|---|---:|---:|---:|---|
| `3I6+6I1` | 5 | 73 | 31 | `24 C2`, `4 C3`, `3 S3` |
| `3I6+6I1` | 7 | 237 | 87 | `66 C2`, `16 C3`, `5 S3` |
| `I2+2I8+6I1` | 5 | 71 | 15 | `15 C2` |
| `I2+2I8+6I1` | 7 | 271 | 43 | `43 C2` |

All six `GF(7)` Golay fibre models carrying the 24 correct marked MW2 pairs
have `C2`, not `C3`, symmetry.  A seventh model carries both individual
section types and has full `S3` `j`-symmetry, but it carries no correctly
paired MW2 basis.  The NS0031 positive model is one of the 43 quadratic
models in its complete `GF(7)` chart.

Consequently repeated supports are a productive filter, not a proof of
descent.  The same audit should be run on any new source before a
base-change construction is used to predict sections.

## Theorem and status boundary

The section splitting is an instance of the `k=1` character decomposition in
Theorem F4 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md):
over a quadratic extension, the rational MW space is the direct sum of the
invariant group and the corresponding quadratic-twist group.  No new general
theorem is asserted here.

This note records exact computations over the stated fields.  It does not
change `MATH_STATUS.json`:

- the rational `3I6` model was already rejected from the determinant-720 NS
  class;
- the NS0031 result remains a finite-field equation precursor;
- absence of cubic symmetry is proved only inside the complete normalized
  weighted-support automorphism tests for the promoted models, not for every
  possible chart or unrelated source.

## Literature bridge

- R. Miranda and U. Persson,
  [*On extremal rational elliptic surfaces*](https://www.math.colostate.edu/~miranda/preprints/Miranda-Persson1986_Article_OnExtremalRationalEllipticSurf.pdf),
  especially the base-change table `I_M -> I_(NM)` and the relation between
  the degrees of the two `j`-maps.
- M. Schuett and T. Shioda,
  [*Elliptic surfaces*](https://arxiv.org/abs/0907.0298), for rational elliptic
  surfaces, Shioda--Tate, heights, and Mordell--Weil lattices.
- Y. Kimura,
  [*F-theory models on K3 surfaces with various Mordell--Weil ranks*](https://arxiv.org/abs/1802.05195),
  for explicit K3 constructions by quadratic base change of rational elliptic
  surfaces and the injection of the rational-surface MW group after pullback.
