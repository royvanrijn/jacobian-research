# Limits of finite-field detection of global invertibility

This note separates two superficially similar questions for weighted Keller
models over localizations of the integers:

1. can a noninvertible characteristic-zero map look bijective on rational
   points modulo many primes? and
2. can inequivalent characteristic-zero maps have indistinguishable
   reductions modulo many primes?

For the repository's weighted marked-root family, the first phenomenon is
impossible at every good prime.  The second can be engineered at any
prescribed finite set of good primes.

## 1. The three point-set tests are equivalent

For a self-map

\[
 F:\mathbb F_p^3\longrightarrow\mathbb F_p^3,
\]

the source and target are finite sets of the same size.  Consequently the
following conditions are equivalent:

* `F` is injective;
* `F` is surjective; and
* every rational target has an `F_p`-rational preimage.

They should not be treated as three independent search objectives.

## 2. Weighted marked-root maps never pass the test

For every admissible weighted seed, the third target coordinate is

\[
 C=x\gamma.
\]

The source divisor above `C=0` is the disjoint union of `x=0` and `gamma=0`.
Over `F_q` these contain, respectively, `q^2` and `q(q-1)` points.  Thus
exactly

\[
 2q^2-q
\]

rational source points map into the target plane `C=0`, which has only `q^2`
rational points.  The root-one chart maps onto that target plane, so the
restriction is surjective but noninjective.  Hence the full map is
noninjective and, because it is a self-map of a finite set, nonsurjective.

This is the same boundary identity recorded in
[`FINITE_FIELD_CHEBOTAREV.md`](FINITE_FIELD_CHEBOTAREV.md):

\[
 \sum_j B_j(q)=q^2,\qquad B_0(q)=0,\qquad
 \sum_j jB_j(q)=2q^2-q.
\]

Therefore an exact finite-field bijectivity test does not give a false
positive anywhere in this family.  The obstruction is elementary and does
not require Chebotarev.  Universal `S_N` monodromy gives a second, asymptotic
obstruction: for inverse degree `N>1`, the image density tends to

\[
 1-\frac{D_N}{N!}<1,
\]

where `D_N` is the derangement number.  Large good fields therefore expose a
positive density of omitted targets, not a nearly invisible exceptional set.

The same argument applies after replacing `F_p` by any finite extension.  It
also explains why being a permutation on infinitely many extensions is not
the relevant distinction here: these weighted maps are permutations on no
finite extension of a good residue field.

## 3. A finite set of primes cannot detect equivalence

The degree-five stable-moduli family has

\[
 H_\lambda(W)=
 \frac{W^2(W-1)(3W^2-(5\lambda+1)W+3\lambda)}{60}.
\]

Its stable Hessian invariant is

\[
 J(\lambda)=
 \frac{64(50\lambda^2-40\lambda+17)^3}
 {75(10\lambda^2-8\lambda+1)^2}.
\]

Stable polynomial left-right equivalence implies equality of `J`; see
[`DEGREE_FIVE_STABLE_MODULI.md`](DEGREE_FIVE_STABLE_MODULI.md).

Let `S` be any finite set of primes at which the model with parameter
`lambda` has certified good reduction, and put

\[
 M=\prod_{p\in S}p,\qquad \mu=\lambda+kM
\]

for a positive integer `k`.  Only finitely many choices of `k` put `mu` on the
explicit exceptional locus or make `J(mu)=J(lambda)`, so a valid `k` always
exists.  The two resulting characteristic-zero Keller maps are stably
inequivalent.  Nevertheless `mu=lambda mod p` for every `p in S`.  Because
the family coefficients are rational functions of the parameter with fixed
denominators and parameter-unit denominators, the two polynomial maps reduce
coefficientwise to the same map at each selected good prime.  In particular
they have identical:

* point functions on every `F_(p^r)^3` over those residue characteristics;
* full rational fiber-cardinality histograms;
* squarefree residue-degree/factorization distributions; and
* all statistics derived from any bounded collection of extension fields.

This construction defeats **any predetermined finite battery of reduction
tests**, not only the five-prime example below.  Its cost is arithmetic
height: the second parameter grows at least as the product of the tested
primes.  It does not produce an infinite or positive-density coincidence for
one fixed pair.

## 4. Reproducible example

The default experiment takes

\[
 S=\{11,13,17,19,23\},\qquad
 \lambda=2,\qquad
 \mu=2+11\cdot13\cdot17\cdot19\cdot23=1062349.
\]

It verifies the two good-reduction certificates, coefficientwise equality of
the reductions, inequality of the characteristic-zero Hessian invariants,
and the complete fiber-cardinality histogram at every selected prime:

```bash
.venv/bin/python scripts/explore_finite_field_detection.py
```

The conclusion is narrower than the proposed invertibility heuristic.  In
this repository, exact modular bijectivity is a reliable negative test because
the maps fail it at every good prime.  What fails spectacularly is the
stronger inference that agreement of many finite-field fingerprints is
evidence for characteristic-zero equivalence.  A finite test set can always
be hidden in the parameter congruences.

## 5. Fixed pairs and the remaining open search

The first direction suggested by the congruence construction is now solved.
The involution `mu=4/5-lambda` in
[`ARITHMETIC_INDISTINGUISHABILITY.md`](ARITHMETIC_INDISTINGUISHABILITY.md)
produces a fixed stably inequivalent pair whose exact full fiber data agree at
every common good prime and over every finite residue-field extension.  The
agreement comes from an affine isomorphism of the inverse pencils on `C!=0`
and a separate exact identification of the direct `C=0` fiber equations.

This leaves the genuinely different direction:

* leave the weighted marked-root family and seek an integral etale Keller
  family without the `C=0` collision excess, where permutation reductions at
  isolated primes are not ruled out structurally.

Universal `S_N` leading terms cannot distinguish same-degree seeds, so any
positive-density result must control the lower-order boundary and
discriminant contributions, not merely the Chebotarev main term.
