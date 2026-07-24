# Prime-to-characteristic Keller counterexamples in every positive characteristic

> **Status.** This is a proposed broader-paper draft and theorem note. It does
> not modify the compiled manuscript
> [`Finite Étale Algebras as Keller Fibers`](main.tex). The polynomial
> identities are checked by
> [`../../scripts/verify_prime_to_characteristic_realization.py`](../../scripts/verify_prime_to_characteristic_realization.py).

## Proposed broader title

**Finite Étale Keller Fibers and Prime-to-Characteristic Counterexamples**

## Proposed abstract

A Keller map is a polynomial endomorphism of affine space with nonzero
constant Jacobian determinant. Over characteristic-zero fields, the companion
manuscript realizes every finite étale algebra of rank one or at least three as
a complete fiber of an explicit Keller map of affine three-space. This note
records the positive-characteristic existence theorem suggested by the same
marked-root viewpoint.

For every field `k` of characteristic `p>0`, we give an explicit noninjective
Keller map

```text
A^3_k -> A^3_k
```

whose generic function-field extension is separable and has degree prime to
`p`. Two formulas suffice. An integral normalized-factorization cubic has
Jacobian one, generic separable degree three, and a complete split
three-point fiber after base change to every field. It handles every
characteristic except `p=3`. In characteristic three, a ten-monomial
weight-redistributed quartic has Jacobian one, generic separable degree four,
and a visible prime-field collision. Consequently the usual Adjamagbo, or
separable, positive-characteristic Jacobian conjecture fails in dimension
three in every positive characteristic. Identity stabilization gives the same
conclusion in every dimension at least three.

The theorem is an existence result, not yet a positive-characteristic analogue
of the full finite-étale realization theorem: arbitrary prescribed finite
étale algebras in bad characteristic require additional polynomiality and
marking arguments.

## 1. Main theorem

Let `k` be a field of characteristic `p>0`. For a dominant polynomial map
`F:A^3_k -> A^3_k`, write

```text
gdeg(F) = [k(x,y,z):k(F_1,F_2,F_3)].
```

The formulation of the separable Jacobian conjecture used here asks that a
Keller map be invertible when `p` does not divide its generic degree. This is
the formulation used in the recent characteristic-two counterexample of
Huq-Kuruvilla and goes back to Adjamagbo.

### Theorem 1.1 — prime-to-characteristic realization

For every field `k` of characteristic `p>0`, there is an explicit polynomial
map

```text
F_p:A^3_k -> A^3_k
```

such that:

1. `det DF_p=1`;
2. `F_p` has two distinct `k`-rational points with the same image;
3. the generic extension is separable; and
4. its degree is

   ```text
   gdeg(F_p)=3,  if p!=3,
   gdeg(F_p)=4,  if p=3.
   ```

In particular `p` does not divide `gdeg(F_p)`. Hence the separable Jacobian
conjecture is false in dimension three over every field of positive
characteristic.

The proof uses the integral cubic in Section 2 when `p!=3` and the
characteristic-three quartic in Section 3 when `p=3`.

### Corollary 1.2 — all higher dimensions

For every `n>=3`, adjoining `n-3` identity coordinates gives a counterexample
of the same generic degree on `A^n_k`.

## 2. One integral cubic for every characteristic

The first construction is an integral chart of the normalized
linear-times-quadratic factorization model. It is related to the repository's
[normalized factorization model](../../verified/NORMALIZED_FACTORIZATION_MODEL.md),
but unlike the displayed characteristic-zero chart there, it has no
fractions.

Put

```text
beta  = xz-2y,
gamma = -xz+3y,
```

and define

```text
a = x,
b = 1+x beta,
c = 1+x gamma,
d = -(beta+gamma)-x beta gamma,
e = -z-2 beta^2-4 beta gamma-2x beta^2 gamma.
```

Now set

```text
Phi_3(x,y,z) = (P,Q,R) = (ac, ae+bd, be).                 (2.1)
```

All coefficients are integers, so this is one map over `Z` and may be reduced
or base-changed to any field.

### Proposition 2.1 — normalized factorization

The identities

```text
ad+bc=1,
a^2e-abd+b^2c=1                                      (2.2)
```

hold in `Z[x,y,z]`. Consequently

```text
P T^3+T^2+Q T+R = (aT+b)(cT^2+dT+e).                   (2.3)
```

The second equation in (2.2) is the resultant-one condition for the displayed
linear and quadratic factors.

There is also a global integral inverse from the normalized factorization
variety

```text
X = V(ad+bc-1, a^2e-abd+b^2c-1) subset A^5_Z          (2.4)
```

to the source chart. Namely `x=a` and

```text
y = 2ad^2+ace+6abd^2+3abce-4ae-bd,                    (2.5)

z = 4d^2+2ce+12bd^2+6bce-9e.                          (2.6)
```

Direct substitution in both directions proves

```text
X is isomorphic to A^3_Z.                               (2.7)
```

Thus (2.1) is normalized multiplication of a marked linear factor and a
quadratic factor, expressed on an integral affine-three-space chart.

For context, the full coefficient-resultant map

```text
Theta(a,b,c,d,e)
 = (ac, ad+bc, ae+bd, be, a^2e-abd+b^2c)
```

satisfies

```text
det DTheta = -(a^2e-abd+b^2c)^2.                       (2.8)
```

This is the characteristic-free coefficient-resultant étaleness mechanism.

### Proposition 2.2 — Keller identity and collision

Direct differentiation gives the integral identity

```text
det D Phi_3 = 1.                                        (2.9)
```

Moreover

```text
Phi_3(-1, 3,-8)
 = Phi_3( 0,-1,16)
 = Phi_3( 1,-2,-5)
 = (0,1,0).                                             (2.10)
```

The three source points remain distinct after base change to every field. In
fact (2.10) is the complete split fiber. Homogenizing the target cubic gives

```text
T^2S+TS^2 = TS(T+S),                                   (2.11)
```

whose three projective roots are distinct in every characteristic. A point of
the normalized factorization fiber is the choice of one of these three linear
factors; the resultant-one normalization fixes its scale. The three choices
are exactly the three points in (2.10).

### Proposition 2.3 — exact generic degree and separability

Let

```text
K = k(P,Q,R),       L = k(x,y,z).
```

On `x!=0`, put

```text
tau   = -b/x,
delta = 3P tau^2+2tau+Q.                                (2.12)
```

Equation (2.3) gives

```text
P tau^3+tau^2+Q tau+R=0.                                (2.13)
```

The reconstruction identities are

```text
x = delta^(-1),                                         (2.14)
y = P delta^2-2delta-tau,                               (2.15)
z = delta(2P delta^2-5delta-3tau).                      (2.16)
```

Indeed, as polynomial identities in the source function field,

```text
x delta=1.                                               (2.17)
```

Hence `L=K(tau)`.

The universal cubic

```text
D(T)=P T^3+T^2+Q T+R                                   (2.18)
```

is irreducible over `K`: in `k[P,Q,T][R]` it is monic of degree one in the
independent variable `R`, and its coefficients are primitive as a polynomial
in `T` because the coefficient of `T^2` is one. Gauss's lemma applies.
Therefore

```text
[L:K]=3.                                                 (2.19)
```

Finally

```text
D'(tau)=delta=1/x != 0,                                  (2.20)
```

so the extension is separable in every characteristic.

Thus `Phi_3` itself is a separable degree-three counterexample over every
field. For the prime-to-characteristic conjecture it is sufficient whenever
`p!=3`.

## 3. The missing characteristic-three quartic

Let `k` now have characteristic three. Define

```text
Psi_4(x,y,z) =
(
  -y-x^2y^2,
  xy-z+xz^2+x^5y^3-x^6y^3z+x^7y^3z^2,
  x+x^2z
).                                                       (3.1)
```

This formula has ten monomials across its three coordinates.

### Proposition 3.1 — visible collision

Direct substitution gives

```text
Psi_4( 1,0, 0) = (0,0,1),
Psi_4(-1,0,-1) = (0,0,1).                               (3.2)
```

The two points are distinct in characteristic three.

### Proposition 3.2 — weight-redistributed suspension

Put

```text
u     = 1+x^2y,
gamma = 1+xz,
W     = u gamma,
C     = x gamma.                                         (3.3)
```

Use the characteristic-three seed

```text
H(W)=2W^2+W^4,       H'(W)=W+W^3.                       (3.4)
```

Define

```text
S=H'(W)+gamma,
T=WS-H(W).                                               (3.5)
```

Then both quotients in

```text
(A,B,C)=(T/C^2, S/C, C)                                 (3.6)
```

are polynomials, and their expansions are exactly (3.1). Polynomiality may be
seen without expansion:

```text
T/C^2 = -uy,
S/C   = (u+1+u^3 gamma^2)/x.                             (3.7)
```

The numerator in the second quotient vanishes at `x=0` and expands to the six
terms displayed in (3.1).

The determinant is the product of three elementary factors:

```text
det d(W,gamma,C)/d(x,y,z) = x^3 gamma^2,
det d(T,S,C)/d(W,gamma,C) = gamma,
det d(A,B,C)/d(T,S,C)     = C^(-3).                     (3.8)
```

Since `C=x gamma`,

```text
det D Psi_4=1.                                           (3.9)
```

### Proposition 3.3 — exact generic degree and separability

For target coordinates `(A,B,C)`, the marked inverse coordinate satisfies

```text
E(W)=W^4+2W^2-BCW+AC^2=0.                               (3.10)
```

This quartic is irreducible over `k(A,B,C)`. Indeed, over `k(B,C)` it is a
polynomial of degree one in the independent parameter `A`, with nonzero
coefficient `C^2`; equivalently its quotient ring is a polynomial ring in
`W`. Gauss's lemma then gives irreducibility over the fraction field.

On the source incidence,

```text
E'(W)=W^3+W-BC=-gamma.                                   (3.11)
```

For a generic root, reconstruction is

```text
gamma = BC-(W+W^3),
x     = C/gamma,
u     = W/gamma,
y     = (u-1)/x^2,
z     = (gamma-1)/x.                                    (3.12)
```

Consequently

```text
k(x,y,z)=k(A,B,C)(W),
[k(x,y,z):k(A,B,C)]=4,                                  (3.13)
```

and (3.11) proves separability.

Since `3` does not divide `4`, (3.1) is a counterexample to the
prime-to-characteristic formulation in characteristic three.

## 4. Proof of the main theorem

If `p!=3`, use the base change of the integral cubic `Phi_3`. It has Jacobian
one, a prime-field collision, and generic separable degree three. Because
`p` does not divide three, it satisfies the prime-to-characteristic
hypothesis.

If `p=3`, use `Psi_4`. It has Jacobian one, a prime-field collision, and
generic separable degree four. Because three does not divide four, it also
satisfies the hypothesis.

These two formulas cover every field of positive characteristic and prove
Theorem 1.1.

## 5. Relation to the active finite-étale-fiber paper

The active manuscript proves a stronger *fiber-realization* theorem in
characteristic zero: every finite étale algebra of rank one or at least three
occurs as a complete Keller fiber. Its principal construction is the
root-engineered quadratic gauge. The theorem above has a different strength:
it is uniform in the characteristic but only asserts the existence of a
prime-to-characteristic noninvertible Keller map.

The common principle is marked-root algebraization:

1. represent an inverse fiber by one polynomial equation;
2. retain one root on the source;
3. arrange that the derivative of the inverse polynomial is the
   reconstruction coordinate;
4. cancel its zero in the incidence Jacobian against a reciprocal source or
   target chart; and
5. forget the root on the target.

The two cases use different realizations of this principle:

| characteristic | inverse degree | realization | relation to existing work |
|---|---:|---|---|
| `p!=3` | `3` | integral normalized linear-times-quadratic factorization | characteristic-free integral chart of the foundational factorization model |
| `p=3` | `4` | weight-redistributed tangent suspension | specialization of the complementary chart already verified by `verify_hasse_typical_seed_recovery.py` |

The characteristic-two reduction of the integral cubic and the map of
Huq-Kuruvilla arise from the same normalized cubic-factorization mechanism.
The characteristic-three quartic is not a cubic reparametrization: its generic
degree is four, so it is not polynomially left-right equivalent to a
degree-three map.

### What is not proved here

This draft should not yet replace the characteristic-zero scope of the active
paper by the claim that every finite étale algebra of rank at least three is
realizable in every characteristic. Remaining issues include:

- bad coefficients in the quadratic gauge when the characteristic divides
  the normalization constants;
- existence of a suitable rational translation over finite base fields;
- global scheme-theoretic reconstruction for arbitrary prescribed seeds in
  the complementary charts; and
- the correct rank classification in positive characteristic.

A safe broader manuscript would retain the characteristic-zero realization
theorem and add Theorem 1.1 as a separate positive-characteristic headline.

## 6. Literature and provenance boundary

Adjamagbo proposed the positive-characteristic formulation excluding the
classical inseparable degree-`p` pathology:

- K. Adjamagbo, *On separable algebras over a U.F.D. and the Jacobian
  conjecture in any characteristic*, in **Automorphisms of Affine Spaces**,
  1995, pp. 89--103,
  <https://doi.org/10.1007/978-94-015-8555-2_5>.

Huq-Kuruvilla gave the recent characteristic-two degree-three counterexample:

- I. Huq-Kuruvilla, *An Explicit Characteristic-2 Counterexample to the
  Separable Jacobian Conjecture*, arXiv:2607.20968,
  <https://arxiv.org/abs/2607.20968>.

That paper states that, to the author's knowledge, no other counterexample to
the separable formulation was then known. A dated search on 24 July 2026 for
`separable Jacobian conjecture`, `Adjamagbo counterexample`,
`characteristic 3`, `prime-to-p`, and `degree 4` found no earlier
characteristic-three or all-positive-characteristic theorem matching the two
formulas above. This is a qualified literature statement, not an absolute
priority claim.

The integral cubic is an integral chart derived here from the repository's
normalized factorization model; no claim of literature-wide priority is made
for that chart. The characteristic-three formula is a new small specialization
of the repository's already proved weight-redistributed suspension theorem,
not a new suspension mechanism.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_prime_to_characteristic_realization.py
```

The checker verifies:

- the integral cubic Jacobian over `Z`;
- both normalized factorization equations;
- the full coefficient-resultant determinant `-Res^2`;
- both compositions of the integral source-chart isomorphism;
- the cubic factorization and denominator-cleared reconstruction identities;
- the complete three-point split fiber;
- the characteristic-three Jacobian and ten-monomial count;
- the characteristic-three collision;
- all three suspension determinant factors;
- the inverse quartic and derivative identities; and
- the denominator-cleared quartic reconstruction formulas.

The two irreducibility proofs are structural degree-one-in-a-parameter
arguments and do not depend on a bounded computer search.

## 8. Recommended manuscript integration

A future revision of the active paper could use the following structure:

1. characteristic-zero finite étale realization and rank classification;
2. normalized marked-root mechanisms;
3. arithmetic and Hasse-principle consequences;
4. the prime-to-characteristic theorem above;
5. comparison with the characteristic-two cubic paper; and
6. open positive-characteristic finite-étale realization problems.

Until the arbitrary prescribed-fiber theorem has been extended across bad
characteristics, this Markdown draft should remain separate from the compiled
paper.
