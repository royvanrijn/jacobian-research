# Every degree-two change of fibration has arithmetic MW rank at most 13

Authority: `EC-MESTRE-ALL-DEGREE-TWO-MW13-BOUND-20260907`.

On each of the six determinant468 Mestre K3 surfaces, let `F` denote the
original fibre. **Every Jacobian elliptic fibration over Q whose fibre D
satisfies D.F=2 has arithmetic generic Mordell–Weil rank at most13.**
Consequently no such change of fibration supplies an MW14–16 parent of302.
This applies to all rational divisor classes of that degree, beyond the
nine compiled chord presentations and beyond any bounded section dictionary.

The proof uses the already certified
[full rational NS matrix](MESTRE_FULL_NS_GRAMS_2026-09-07.md), of rank18,
and exhibits three independent rational vertical root directions for every
eligible parity class. It does not classify higher-degree fibrations,
exclude the determinant468 surfaces as possible parents at higher degree,
or exclude other surfaces. It constructs no302 parameter or new equation.

## Reducing all degree-two classes to a finite roster

The supplied integral basis splits the full rational divisor lattice as

\[
 N=U\oplus L(-1),\qquad
 e=F,\quad h=O+F,\quad e^2=h^2=0,\quad e.h=1,
\]

where `L` is positive definite, even, of rank16 and determinant468.
All six parents have the same matrix in their supplied rational divisor
bases. Coordinates below use the reduced integral basis of `L` stored in
the certificate; its change of basis is checked to be unimodular.

Any integral isotropic class with `D.e=2` has the form

\[
 D=k e+2h+w,\qquad k=\frac{(w,w)_L}{4}.
\]

Its divisibility in `N` is

\[
 \gcd\bigl(2,k,(w,L)\bigr).
\]

A fibration with a rational section must have divisibility one. Thus
classes with nonintegral `k` or divisibility two may be excluded from this
Jacobian-fibration question.

For every `z in L`, the following integral isometry fixes `e`:

\[
 e\mapsto e,\quad
 h\mapsto h+z+\frac{(z,z)_L}{2}e,\quad
 x\mapsto x+(x,z)_L e\quad(x\in L).
\]

It sends the `L` coordinate `w` of D to `w+2z`. Therefore representatives
of all degree-two isotropic classes can be chosen with all sixteen
coordinates of `w` in `{0,1}`. These are lattice isometries; no claim that
they preserve the nef cone or arise from surface automorphisms is needed.
They transport integral orthogonal roots and their independence, which
suffices for the upper bound on any actual nef fibre class.

The verifier checks all **65,536** binary words. Exactly **32,509** have
integral `k` and divisibility one. This is a complete roster of parity
representatives, not a count of distinct elliptic fibrations. Of those,
32,508 are nonradical for the pairing modulo2. The remaining word has mask
25805 and squared norm12; it is included in the same witness verification.

## Explicit roots, without trusting an enumeration-completeness claim

For every admitted word the compact witness archive supplies three vectors
`n in L` of norm8. Each has either parity zero or parity `w`.

For parity zero, the verifier sets

\[
 b=0,\qquad x=n/2,\qquad a=(w,n)_L/4.
\]

For parity `w`, it sets

\[
 b=1,\qquad x=(n+w)/2,\qquad
 a=((w,n)_L+2k)/4.
\]

Every numerator divisibility is checked exactly. In both cases the resulting
integral rational divisor class

\[
 R=a e+b h+x
\]

satisfies `R^2=-2` and `R.D=0`, verified directly in the original NS matrix.
No assertion that R itself is an irreducible curve is necessary.

On `D^perp`, the map `a e+b h+x -> x-(b/2)w` has kernel `Q D`.
The three witness images are `n/2`. A supplied three-column minor has
determinant nonzero modulo101, independently proving their rational
linear independence. Hence the roots are independent modulo `Q D`.

For an actual nef fibre D, K3 Riemann–Roch makes either R or `-R` effective.
Its zero intersection with D then makes it vertical. All these classes
are rational divisors, so they give at least three rational vertical
directions beyond the fibre class. Taking rational Galois invariants in
the Shioda–Tate decomposition yields

\[
 \operatorname{rank}E/\mathbb Q(s)
 \leq 18-2-3=13.
\]

The standard divisor and rank results are reviewed in
[Schütt–Shioda, *Elliptic surfaces*, Corollary6.13 and Section12](https://arxiv.org/abs/0907.0298).
The parity reduction and explicit witnesses above specialize those tools
to the pinned rational NS lattice. They do not assume a Tate-conjecture
converse, full splitting of the new singular fibres, or an exact rank
for any new presentation.

## Reproduction

The producer
[`certify_mestre_degree_two_rank_gate.sage`](../cas/certify_mestre_degree_two_rank_gate.sage)
finds witnesses among old norm2 and norm8 vectors, under a120-second,
one-worker cap. It writes the
[certificate](../../artifacts/generated-results/elliptic-curves/mestre_degree_two_rank_gate_v1.json)
and the132-KB compressed witness archive linked there.
The producer never overwrites either output.

The standalone verifier imports no producer code and performs no
short-vector enumeration. It independently checks the integral U splitting,
all65,536 admissions, all97,527 root identities and all32,509 independence
minors. Completeness of the producer's PARI enumeration is unnecessary:
three exhibited witnesses per admitted class prove the result.

```sh
sage -python elliptic-curves/cas/verify_mestre_degree_two_rank_gate.sage
```

The fresh-directory
[replay record](../../artifacts/generated-results/elliptic-curves/mestre_degree_two_rank_gate_replay_v1.json)
binds the verifier, its three input files and its passing transcript.

The preceding exploratory gates inspected64 degree-two classes and64 each
at old degrees3,4,5, finding only rank6–11 frame configurations. Their
inputs and scripts are retained as finite exploratory records. The present
theorem replaces the sampled degree-two gate with a complete rank-at-most13
bound; it does not turn the higher-degree samples into exclusions.

For MW14–16 construction on this cohort, a different fibration must now
have old degree at least3. Degree one is impossible for two distinct
elliptic fibrations: it would give a degree-one map from a smooth old
genus-one fibre to the new base P1. The higher-degree equation and302
specialization problems remain open.
