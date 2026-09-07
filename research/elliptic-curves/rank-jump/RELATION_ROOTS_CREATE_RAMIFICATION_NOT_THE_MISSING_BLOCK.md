# A relation root leaves the old span, but adds no Selmer class

Follow-up: [all generic chord and tangent roots on four controls](GENERIC_CHORD_BLOCKS_HAVE_PRIVATE_RAMIFICATION.md)
have independent private good-prime obstructions. The 1123-form test includes
the fresh 103b2 pair and historic +12/+14 fibres; every whole span adds zero
Selmer classes despite shared three-root local cancellations.

The explicit square root of the sole inherited relation in the retained
equation-only norm dictionary creates a **new ramified squareclass**. It
does not create an unramified class, even after multiplication by arbitrary
products of all 4134 dictionary elements and all 16 generic classes.

This closes one constructive possibility left open by
[exhaustion of the original dictionary](EARLY_RELATION_POOL_EXHAUSTION.md).
It is a single-fibre constructor audit on MW16-05 at t=3/17, not a new
high-versus-low comparison. The class-creation event required by the
[large-jump capacity theorem](CLASS_CREATION_REQUIRES_NEW_UNRAMIFIED_COVERS.md)
remains unconstructed. No exceptional point is an arithmetic input.

## Why this operation needed its own test

Let K be the cubic two-division field, and put pi(a)=N(a)a. Let G be the
16-dimensional generic Kummer subgroup and V the span of G together with
the 4134 retained norm projections. The previous exhaustive parity test
has rank 4133. Its sole coefficient relation is inherited:

\[
 \pi(\alpha)\prod_{i\in I}\gamma_i=w^2,
 \qquad
 I=\{2,5,8,12,13,15\}.
\]

Indices here are one-based. The retained primitive address of alpha is
(20941,464), the generic mask is 22674, and the exact coefficients of w
are in the earlier
[root certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_retained_norm_inherited_hit_v1.json).
The earlier result excludes new classes **within V**. It does not exclude
pi(w): taking an actual square root is not taking another F2-linear
combination of squareclasses. Also N(pi(w))=N(w)^4, so this is a legitimate
norm-square candidate in H1(Q,E[2]). The sign of w is immaterial because
pi(-w)=pi(w) in a cubic field.

This construction uses fixed representatives. Replacing the parent
product by its product with z^2 replaces pi(w) by pi(w)pi(z). Thus it is
not a canonical operation on squareclasses and does not exhaust every
possible choice of representatives. In particular, allowing arbitrary z
would smuggle an unrestricted new class into the construction.

## The exact obstruction is already at good primes

The [first protocol](RELATION_ROOT_CLASS_PROTOCOL.json) freezes the root
and at most the first 32 good rational primes in the retained parent
alpha's complete ideal support. There are five such primes. No new norm
is factored. Every displayed tuple lists all prime ideals above p, in the
retained decomposition order.

| p | Residue degrees | v_P(w) | v_P(pi(w)) |
|---:|---|---|---|
| 113 | (1,1,1) | (1,2,1) | (5,6,5) |
| 163 | (1,2) | (2,1) | (6,5) |
| 277 | (1,1,1) | (1,1,2) | (5,5,6) |
| 787 | (1,2) | (2,1) | (6,5) |
| 929 | (1,2) | (2,1) | (6,5) |

These primes are outside S and the cubic polynomial discriminant. Each
generic class is a unit at every displayed prime ideal. The seven odd
valuations therefore survive all 2^16 generic corrections. At good odd
places the elliptic Kummer image is unramified, so this excludes the
candidate from Selmer, before any rational-versus-Sha question.

The [completed certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_relation_root_class_v2.json)
contains the exact element and ideal data. Its
[independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_relation_root_class_verification_v1.json)
recomputes the root and norm identities using rational multiplication
matrices, then checks 24 valuations by membership in P^v and nonmembership
in P^(v+1). Ideal powers are independently multiplied integral lattices;
the replay does not use PARI's valuation or ideal-power routines. It also
checks the 192 generic unit conditions. The initial v1 UNKNOWN record is
preserved: arithmetic finished, but JSON serialization rejected a Sage
integer. The completion changes serialization only.

## The entire dictionary cannot repair the defect

The [second protocol](RELATION_ROOT_AFFINE_PROTOCOL.json) tests an affine
system against the whole retained dictionary, not a selection of further
elements. On its 10606 prime-ideal columns over 4569 rational primes,
let B be the outside-S parity matrix of the 4134 norm projections and b
the parity vector of pi(w). The
[result](../../artifacts/generated-results/elliptic-curves/rank_jump_relation_root_affine_v1.json)
is

\[
 \operatorname{rank}B=4133,
 \qquad \operatorname{rank}\begin{pmatrix}B\\b\end{pmatrix}=4134.
\]

Thus no dictionary correction cancels b, even on this restricted prime
set. Untested primes cannot repair an inconsistency on tested primes.
There are 22 odd coordinates in b on this set.

The certificate supplies a dual functional lambda supported on 1992
prime ideals over 1821 rational primes, with

\[
 \lambda(B_i)=0\quad(1\le i\le4134),\qquad
 \lambda(G)=0,\qquad \lambda(b)=1.
\]

The [independent dual replay](../../artifacts/generated-results/elliptic-curves/rank_jump_relation_root_affine_verification_v1.json)
checks every dictionary pairing directly from its sparse principal-ideal
factorization. It checks the root and generic valuations on the dual
support by 33864 independent lattice membership tests. It does not use
the worker's elimination or dual-construction algorithm. The retained
complete dictionary factorizations were already independently certified
by the earlier pool replay.

Consequently the entire coset pi(w)V is disjoint from Selmer. In fact the
previous kernel computation gives V intersect Sel2 = G: every outside-S
unramified combination is generic. Hence, with V'=span(V,pi(w)),

\[
 \boxed{V'\cap\operatorname{Sel}_2(E/\mathbb Q)=G,
 \qquad V'\cap U=G\cap U.}
\]

This is an exact exclusion of all combinations in one explicitly defined
space. It is not a bound on the full Selmer group or the class group.
The known additional strict directions of this fibre remain outside it.

## A general mod-four gate for this constructor

For any explicitly supplied relation Pi=w^2 in a cubic field, write
e_P=e(P/p). The elementary norm identity gives

\[
 v_P(\pi(w))
 =\frac{v_P(\Pi)+e_Pv_p(N\Pi)}{2}.
\]

Thus its root projection is unramified at a good odd prime P exactly when

\[
 v_P(\Pi)+e_Pv_p(N\Pi)\equiv0\pmod4.
\]

The ordinary mod-two relation does not imply this mod-four condition.
Here seven tested prime ideals have defect 2 modulo 4. Cancellation by
the retained dictionary also fails, by the dual certificate above.
To certify a strict new class through this construction would require
the mod-four condition at **all** good primes, local square certificates
at S including 2 and infinity, and independence modulo the inherited
pool. None of those conditions may be replaced by the existence of the
initial square-root relation.

This gate is labelled **incidence**. It is not a rational-solubility test,
a point-visibility score, or a demonstrated condition on the family
parameter t. It does not show how to obtain many successful relations
simultaneously. The desired parameter-dependent construction of the
[standard S3 unramified class block](STRICT_CLASS_CREATION_IS_A_STANDARD_S3_BLOCK.md)
is still missing.

The useful lesson for candidate selection is a rejection rule for this
specific constructor: certify good-prime parity before assigning any
significance to a new norm-square direction. There is no new rank
predictor to pass to Agent 1. No search settings or live outputs changed.

Replay the independent certificates:

```sh
timeout 30 sage -python elliptic-curves/rank-jump/verify_relation_root_class.py check
timeout 60 sage -python elliptic-curves/rank-jump/verify_relation_root_affine.py check
```
