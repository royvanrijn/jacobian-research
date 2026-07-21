# Open problems

## Direct all-`k` transfer-block basis after the ribbon counterexample

Let `Z_k` be the completion of the monic square/cube factorization scheme
`U^2=V^3` along `(U,V)=(S^3,S^2)`.  The proposed identification with the
divided-power symmetric product of the split cusp ribbon is false: its norm
map kills `X^2` in the `k=2` collided fiber.  See
[C22_CONDUCTOR_RIBBON_AUDIT.md](C22_CONDUCTOR_RIBBON_AUDIT.md).

**Question.** Is `Z_k` nevertheless finite flat of rank `2^k` over the
extracted `S`-base for every `k` in characteristic zero?

A proof must construct a spanning set over the full completed base and rule
out discriminant-supported torsion without using the refuted norm
surjection.  A counterexample may instead exhibit rank jump, torsion, or an
additional Artin deformation for some `k`.  Exact presentations are available
through the bounded low-degree cases, with collided-fiber evidence for
`k=5,6`; they do not establish the all-`k` quantifier.

## Polynomial equivalence of `C24_(m,1)` and generic weighted seeds

Fix an algebraically closed characteristic-zero field and `m>1`.  Let
`F_m:A^3 -> A^3` be a polynomial C24 cancellation map with `r=1`, including
an arbitrary allowed tail `h_q(A)+A^2g(A)`.  Let `G_H` be a generic admissible
weighted-seed map whose inverse degree is `n=m+2`.

**Question.**  Do there exist polynomial automorphisms `alpha,beta` of
`A^3` such that

\[
 G_H=\beta\circ F_m\circ\alpha?
\]

The same question may be asked after adjoining the same number of identity
coordinates.

### What is known

The intrinsic data computed in
[RESOLVENT_RAMIFICATION_SIGNATURE.md](RESOLVENT_RAMIFICATION_SIGNATURE.md)
agree:

- generic degree `m+2` and monodromy `S_(m+2)`;
- generic critical partition `1^(m+1)`;
- one discriminant boundary prime with ramification index and sheet loss two;
- `m-1` geometric unramified boundary primes over a second target divisor,
  with total loss `m-1`.

For `m=1`, the explicit linear transformation in
[MASTER_CANCELLATION_CONSTRUCTION.md](MASTER_CANCELLATION_CONSTRUCTION.md)
identifies C24 with C01.  No transformation is known for `m>1`.  The
ramification argument which settles `r>=2` supplies no obstruction here.

### Exact unresolved condition

Let `bar(X)_F` and `bar(X)_G` be the normalizations of the target affine
spaces in the two source function fields, with their canonical affine-source
opens.  Left--right equivalence requires an isomorphism of finite normal
pairs, after a polynomial target change,

\[
 (\bar X_F,A^3)\simeq(\bar X_G,A^3),
\]

whose restriction to the affine opens is a polynomial source automorphism.
The signature proves equality only of the generic codimension-one numerical
data.  It does not compare:

1. conductor ideals and completed local rings where the discriminant meets
   the second nonproperness component;
2. intersection and specialization of the `m-1` unramified boundary primes;
3. the filtered valuation semigroups `v_E(k[x,y,z]\setminus0)`;
4. dependence of any of these data on the tail `g` or on the weighted seed
   coefficients.

An explicit polynomial pair `(alpha,beta)`, a mismatch in one of items 1--3,
or a proof that the normal-pair isomorphism cannot restrict polynomially would
settle the question.  Bounded symbolic agreement of discriminants or fiber
counts is not sufficient for an all-`m` conclusion.

## Nonmonomial solutions of the generalized cancellation equation

Use the notation of
[GENERALIZED_CANCELLATION_MECHANISM.md](GENERALIZED_CANCELLATION_MECHANISM.md).
For `e=1`, every polynomial leading solution is controlled by

\[
 I_n(q)=\int_0^1u\{1-q(1-u)\}^n\,du.
\]

A genuinely nonmonomial leading term requires a common root of `I_n` and
`I_l` for distinct positive integers `n,l`.  No such common root occurs for
`1<=n,l<=12`, but this is only a bounded exact computation.

**Questions.** Are the polynomials `I_n` pairwise coprime in characteristic
zero?  If not, does a common-root solution have an invertible linearization
and a full cancellation jet with coefficients in `k[y]` rather than merely
in `k[y,f(y)^(-1)]`?  Finally, is any resulting map inequivalent to a C24
monomial member?

For `e>=2`, even the polynomial solutions of the leading functional equation
are unclassified.  No existence claim for a new family is made.
