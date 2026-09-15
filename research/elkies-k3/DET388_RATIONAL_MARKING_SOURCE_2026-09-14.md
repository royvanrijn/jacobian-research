# Determinant 388: a full rational rank-19 marking source

The row K3-f5168aef816b1b25 passes the full arithmetic marking gate.
There exists a projective K3 X/Q with geometric Picard rank 19, all its
geometric divisor classes represented over Q, and actual transcendental
lattice

\[
T=\begin{pmatrix}-2&-6&-1\\-6&2&1\\-1&1&10\end{pmatrix}.
\]

This is an **existence result with exact moduli and lattice data**, not an
explicit K3 equation or an MW17 fibration. The selected rational period is
x=2 on the full marked curve X_0^194(1)/<w_194>, identified with P1/Q.
The subsequent [global root obstruction](DET388_GLOBAL_ROOT_OBSTRUCTION_2026-09-14.md)
proves that this NS has no rootless frame and no MW17 fibration. The full
rational marking remains valid; the geometric MW17 branch is now closed.

## 1. Full stable group, with the literal arithmetic level

The even Clifford order has integral basis
(1,e0e1,e0e2,e1e2) and reduced trace pairing

\[
\begin{pmatrix}2&-6&-1&1\\-6&38&7&-1\\-1&7&11&29\\1&-1&29&-9\end{pmatrix}.
\]

Its reduced discriminant is 194. Its rational quaternion algebra is
(10,97/20), ramified precisely at 2 and 97. Hence this literal order is
maximal: an overorder would divide its reduced discriminant by its index,
below the algebra's discriminant. There is no hidden Eichler or principal
level condition.

The discriminant group of T is cyclic of order 388. In the Smith
coordinates used by the checker, its quadratic generator has value
8781/776 modulo Z. Its orthogonal automorphisms are multiplication by
1,193,195,387. Exact positive-norm normalizer elements are:

| Atkin–Lehner label | Coordinates in the even-order basis | Norm | Action on A_T |
|---|---|---:|---:|
| 2 | (-5,3,-2,-4) | 2 | 195 |
| 97 | (403,56,67,-15) | 97 | 193 |
| 194 | (-1,1,-6,2) | 194 | 387=-1 |

Their conjugations preserve both the order and T integrally. They exhaust
the maximal-order normalizer cosets, since the only ramified primes are
2 and 97. The norm-one curve is X_0^194(1).

It remains essential to check the units, including local determinant
units. The checker forms the quadratic polynomial numerator of
(Ad(q)-I)G^-1, using q-bar divided by nrd(q) for q^-1. Every coefficient
of this numerator is integral. Thus every local order unit acts trivially
on A_T at that prime. This is a polynomial identity, not an enumeration
of a few units. Reduced norms of maximal-order local units cover Z_p^*.

The Atkin–Lehner signs at 2 and 97 are independent before imposing the
marking. Projectivizing the full stable orthogonal group permits precisely
the two simultaneous signs (+,+) and (-,-). The latter is supplied by
the determinant-minus-one stable reflection -Ad(w_194): it is the
reflection in e0, whose square is -2. The mixed signs are excluded. There
are no further cosets. Consequently

\[
C_T=X_0^{194}(1)/\langle w_{194}\rangle.
\]

This is also the canonical curve over Q. The local units are the entire
maximal-order unit groups, with their full reduced-norm images; there is
no additional Cartan or congruence marking. Forgetting the determinant
orientation is the single global quadratic operation w_194, rather than
independent choices at the two primes. The arithmetic period-map inputs
are [Rizov, Section 3.9 and Theorem 3.16](https://arxiv.org/abs/math/0508018)
and the stable-group description in
[Dolgachev, Proposition 3.3](https://arxiv.org/abs/alg-geom/9502005).
This uses the full orthogonal group, not merely its norm-one cover.

## 2. A rational non-CM period, without extrapolating a CM table

[Guo–Yang](https://arxiv.org/pdf/1510.06193), Table 1 and the D=194
equation and CM normalization, give a genus-nine hyperelliptic model
y²=f(x) of X_0^194(1), with w_194 its hyperelliptic involution. Its
quotient coordinate x is defined over Q. The polynomial coefficients
in increasing order are

```text
-19,92,-286,592,-921,1016,-872,-460,1545,-1752,34,
1752,1545,460,-872,-1016,-921,-592,-286,-92,-19.
```

The [retrieved model](../artifacts/generated-results/elkies-k3-det388-marking-v1/published-model.m)
also gives the quotient as Y²=XT, with projection [x²:x:1]. Its bytes,
source URL and hash are retained. The canonical identification is a
published theorem input; the checker verifies the polynomial and every
subsequent arithmetic calculation.

At x=2 the quotient point is [4:2:1] and

\[
f(2)=-315538147=-7411\cdot42577.
\]

Both factors are prime. Thus a lift P has quadratic residue field
F=Q(sqrt(-7411*42577)), ramified at 7411.

Suppose P were CM by an order R in an imaginary quadratic field K.
[González–Rotger, Theorem 5.8](https://arxiv.org/pdf/math/0612732)
gives H_R=K Q(P), where H_R is the ring class field. Hence h(R)=[H_R:K]
is at most two, and F is a subfield of H_R. The complete class-number-one
and class-number-two classification bounds the absolute fundamental
discriminant of K by 427; see
[Watkins's account](https://magma.maths.usyd.edu.au/~watkins/papers/CLASSNO1.pdf).
The ring class number formula also rules out 7411 in the conductor:
such a prime would contribute at least (7411-1)/3>2, where three is the
largest possible unit-index denominator. Thus H_R is unramified at 7411,
contradicting its subfield F. P, and hence x=2, is non-CM.

This argument uses the complete small-class-number theorem and class-field
ramification. It does not infer non-CM status from absence in a bounded
table. The selected point avoids the Picard-rank-20 locus.

## 3. The actual primitive NS and descent to Q

Let F17 be the positive definite Gram matrix of the retained source frame
K3-f5168aef816b1b25-F001, and put S=U+(-F17). The certificate stores the
full matrix. It has signature (1,18) and determinant 388. The Smith
quadratic generators of S and T have values

\[
q_S=-1054533/776,\qquad q_T=8781/776.
\]

Multiplication by 75 supplies a gluing anti-isometry:
q_S+75²q_T is integral. The checker constructs the graph overlattice of
S+T, of index 388, and verifies its integral even Gram matrix is
unimodular of signature (3,19). It is therefore the K3 lattice. The
embedded S and T remain primitive and are actual orthogonal complements;
no unsaturated subgroup is being substituted for NS.

Choose an ample chamber in S and a very small cone within it. The
algebraic lattice-polarized moduli stack over Q is separated in
characteristic zero; see
[Bragg–Brakkee–Várilly-Alvarado, Sections 3 and 5](https://arxiv.org/html/2510.11477v1).
The Torelli period description, together with the canonical arithmetic
period map above, identifies its non-CM component with the corresponding
open part of C_T. The complex period-surjectivity exclusions are precisely
extra algebraic classes here, hence CM points. Thus x=2 lies in the image.

There is no residual automorphism obstruction to descending this object.
An automorphism fixing S fixes an ample class and therefore has finite
order. On the rank-three transcendental Hodge structure, a finite Hodge
isometry other than +I or -I would have a nonzero rational eigenspace
of type (1,1): its real eigenvalue is +1 or -1, and any eigenvalue on
H^(2,0) different from it leaves that rational eigenspace orthogonal to
the period plane. This contradicts NS=S. The alternative -I is excluded
by primitive gluing because it acts nontrivially on A_T of exponent 388,
whereas S is fixed. Hence the automorphism acts trivially on all H² and
is the identity by Torelli. The marked moduli stack has trivial inertia
at this point, so its rational coarse point is an actual object over Q.
Equivalently, the unique marked isomorphisms give an effective descent
datum. We obtain a projective K3 over Q with NS_bar=S and every class
Galois invariant.

## 4. Invariant classes are actual rational divisor classes

There is one more arithmetic distinction to discharge. The first two
frame coordinates give two classes r,s in S with

\[
r²=s²=-2,\qquad r\mathbin{\cdot}s=1.
\]

For a Galois-invariant geometric line bundle L, its obstruction in Br(Q)
has index dividing chi(L): its cohomology spaces carry the corresponding
twisted descent, and their dimensions are multiples of that index.
K3 Riemann–Roch gives chi(r)=chi(s)=1, so both line bundles descend.
Their intersection supplies a zero-cycle class of degree one on X.
Restriction and corestriction at its closed points therefore make
Br(Q) to Br(X) injective. The Picard descent exact sequence now shows
that **every** invariant geometric line bundle descends. Thus the full
saturated NS, rather than only a finite-index subgroup, is rational.

## Endpoint and next gate

This proves existence of a fully rational rank-19 marking with determinant
388, different from the 948 and 1092 comparison surfaces. It does not yet
give the requested second MW17 mechanism. The retained frame is labelled
2A1+A3 and MW rank 12; a single frame is not an exhaustive classification
of primitive U embeddings. Rootless-frame existence was left UNKNOWN by this arithmetic admission.
The subsequent global root obstruction now settles it negatively. No new
Weierstrass equation or list of17 sections is claimed.

```sh
sage -python research/elkies-k3/scripts/certify_det388_marking.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det388-marking-v1/certificate.json)
replays the integral order, full normalizer action, all-local-unit
polynomial identity, non-CM arithmetic, primitive gluing, and descent
intersection pair. Moduli descent, class-field theory, and the published
canonical model are explicit written theorem inputs. Independent
implementation, formal verification, external review, and newness in the
literature are unclaimed.

## Bounded rootless preflight

After arithmetic admission, a deterministic 128-neighbour probe at primes
3 and5 found no rootless frame. Every transport and nonempty root set
replays; the smallest root count was five sign pairs. This is not an
exhaustive genus or U-embedding classification and alone leaves MW17 UNKNOWN; the later global theorem closes that question.
The [packet](../artifacts/generated-results/elkies-k3-det388-rootless-probe-v1/checkpoint.json)
retains all Grams and transports. An initial enumeration-control failure
(`qfminim` with zero storage capacity) is explicitly withdrawn in its
[failure record](../artifacts/generated-results/elkies-k3-det388-rootless-probe-v1/failure.json);
the corrected run checks total root counts, stored representatives and the
known eight root pairs of the source. This false hit was never promoted to
a mathematical claim. No larger search is inferred from this miss.
