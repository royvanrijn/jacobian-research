# The residual fourth lift has no constant quadratic carrier

The fourth-root field over the degree-eleven component of the successful
relation is **not** obtained by adjoining one constant quadratic field.
Its degree-22 field has no quadratic subfield. Its normal closure has
the full signed permutation group

\[
\boxed{C_2\wr S_{11}=C_2^{11}\rtimes S_{11},}
\]

of order 81,749,606,400. The eleven conjugate fourth-value squareclasses
are independent over the splitting field of the degree-eleven polynomial.
This is a statement about fields of definition of the finite construction,
**not eleven independent Mordell–Weil directions**.

The [preceding fourth-lift calculation](FOURTH_DIRECTION_SPLITS_ONLY_ON_THE_RATIONAL_COMPONENT.md)
gave the complete finite component pattern 1+1+22. The two rational
components at t=−288/65 still supply the two signs of one additional
quotient direction. This result characterizes the remaining degree-22
component; it does not change that rational-point or rank accounting.

## The proposed shared-field mechanism and its exact test

Let F be the irreducible degree-eleven factor of the positive relation
polynomial. Put

\[
K=\mathbf Q[t]/F,\quad b=q_C(\bar t),\quad
L=K(\sqrt b),\quad N=N_{K/\mathbf Q}(b),\quad D=\operatorname{disc}(F).
\]

Here C is the retained `13109` cover. The previous certificate proves
that b is nonsquare, [K:Q]=11 and [L:Q]=22. We test whether the same
fourth root can be supplied by a constant quadratic extension M/Q:

\[
L=KM\quad\Longleftrightarrow\quad b/N\in K^{*2}. \tag{1}
\]

To prove (1), if L=K(sqrt d) for rational d, then b/d is square in K.
Taking norms shows N/d^11, hence N/d, is a rational square. Thus b/N
is square in K. Conversely, if b/N is square then L=K(sqrt N).
The odd degree of K ensures that a nonsquare rational d remains
nonsquare in K.

This also tests every possible quadratic subfield of L: if L contained
such an M, its compositum with K would already have degree 22 and equal L.

## A degree-one local place disproves (1)

The first admissible prime in the frozen bounded test is 59. The exact
reductions are

\[
F(4)=0,\quad F'(4)=51,\quad q_C(4)=7,\quad N=8
\qquad\text{in }\mathbf F_{59}.
\]

All coefficients are 59-integral and the discriminant and N are units.
The simple root gives a completion of K isomorphic to Q59. At this
place,

\[
\frac{b}{N}\equiv\frac78\equiv23\pmod{59},\qquad
23^{29}\equiv-1\pmod{59}.
\]

So b/N is nonsquare in K. More concretely, b itself has a local square
root because 19²≡7, while N≡8 is nonsquare. The fourth lift splits at
this place; the only possible constant-field candidate K(sqrt N) does
not. This certifies that L has **no quadratic subfield** without any
class-group computation.

## Why all eleven conjugate squareclasses are independent

Let E be the splitting field of F, with Galois group S11. The previous
S11 certificate is replayed using irreducibility modulo 73 and factor
degrees 1,2,3,5 modulo 79. The latter supplies a transposition by taking
the fifteenth power of Frobenius. Transitivity of prime degree gives
primitivity, and a primitive group containing a transposition is S11.

The independent norm/discriminant calculation gives

\[
N>0,\qquad D<0,\qquad N\notin\mathbf Q^{*2}.
\]

Therefore N and D are different rational squareclasses. S11 has only
one nontrivial quadratic quotient, its sign, so the only quadratic
subfield of E is Q(sqrt D). In particular N is not square in E.

Write b1,…,b11 for the conjugates of b and consider the S11-equivariant
map

\[
\mathbf F_2^{11}\longrightarrow E^*/E^{*2},
\qquad e_i\longmapsto[b_i].
\]

Let R be its kernel. A submodule containing a nonconstant vector also
contains a vector ei+ej: subtract the vector obtained by exchanging two
coordinates with different bits. The orbit of ei+ej spans the even-weight
augmentation subspace U, of dimension ten. Consequently the only
submodules are 0, the constant line, U, and the whole space. This uses
that eleven is odd. The verifier checks the transposition argument on
all 2,046 nonconstant binary vectors.

The constant vector maps to [N]≠0 in E*/E*². Thus R can only be 0 or U.
If R=U, every [bi] is the same and [b1]=[N], so b/N is square in E.
We now rule that out as well.

The stabilizer of K inside S11 is S10. Its unique subgroup of index two
is A10, so the only quadratic extension of K contained in E is K(sqrt D).
If b/N were square in E but not in K, it would therefore differ from D
by a square in K. Taking norms would give

\[
N_{K/\mathbf Q}(b/N)=N^{-10}
\equiv D^{11}\pmod{\mathbf Q^{*2}},
\]

which is impossible: the left side is a rational square, while D is
nonsquare (indeed negative). The local test already excluded a square
in K. Hence b/N is not square in E, R≠U, and R=0.

Kummer theory now gives

\[
[E(\sqrt{b_1},\ldots,\sqrt{b_{11}}):E]=2^{11}.
\]

The normal-closure group embeds in the signed permutation group,
surjects onto S11, and has its full sign kernel. Its order equals the
order of C2 wr S11, so it is that full group.

The norm field Q(sqrt N) does occur in the **normal closure**, through
the product of all eleven fourth roots. It does not occur in the
degree-22 field L. Confusing these two fields would create a false
shared-quadratic explanation.

## Consequence for the solubility mechanism

The successful finite construction separates a rational parameter
component from a degree-eleven companion. Its fourth root splits over
the rational component, while the companion has the largest possible
sign kernel over its S11 normal closure. There is no simultaneous
collapse of all companion fourth-root fields to one quadratic extension.

This excludes one specific **shared-descent-field explanation** of the
residual component. It does not prove probabilistic independence of the
rational splitting events. Nor does it exclude a different geometric
construction accounting for the original quotient directions.

The positive conclusion remains the previously certified conditional
chain: a rational realization of the fixed successful triple relation
forces t=−288/65, where all four native covers split and their quotient
span has rank three. The missing prospective mechanism is still what
forces a rational parameter component together with its split fourth
lift. The companion's large Galois group does not supply that event,
and its eleven squareclasses are not a rank feature for Agent 1.

This work concerns **solubility and field structure**. It neither alters
incidence accounting nor introduces a visibility score. The other four
retained quotient directions of the +7 fibre remain unexplained by this
construction.

## Replay

The [protocol](FOURTH_LIFT_SHARED_FIELD_PROTOCOL.json) permits only primes
at most 199 and the one fixed residue field. A witness was found at 59,
the first admissible prime; no number-field class group or online job
was needed.

```sh
sage -python elliptic-curves/rank-jump/fourth_lift_shared_field.py check
sage -python elliptic-curves/rank-jump/verify_fourth_lift_shared_field.py check
```

Certificates `rank_jump_fourth_lift_shared_field_v1.json` and
`rank_jump_fourth_lift_shared_field_verification_v1.json` are retained
under `artifacts/generated-results/elliptic-curves/`. The verifier
recomputes norms and discriminants by multiplication determinants and
checks the local arithmetic using integer modular operations. It also
replays the two base Galois patterns and the binary-module proof gates.
