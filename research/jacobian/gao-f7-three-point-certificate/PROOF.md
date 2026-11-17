# Gao F7: the cubic-cube coefficient obstruction

14 November 2026. Registered claim **GAOF7CUBE1**.

For exactly Gao's F7 over C, its degree-twelve eliminant never factors as
R=Q3^3*S3 for arbitrary monic cubics Q3,S3. Repeated roots and roots shared
between the two cubics are included. Thus every eliminant satisfies
sum_r floor(m_r/3)<=2, where m_r are its root multiplicities. In particular
there are at most two distinct triple-or-higher roots, and no root has
multiplicity greater than eight; sharpness is not asserted.

For a three-point fibre this excludes the nonsimple multiplicity lists
(9), (6,3), (3,3,3). This obstruction alone leaves five other lists. The
subsequent [GAOF7MULT1 proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md)
sharpens the individual multiplicity bound to five and also excludes(7,2);
four lists remained at that stage. The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 fivefold
targets and leaves only(4,3,2),(3,2,2,2) open. The separate [verifier](verify_rank_integer.py) completed
[761 exact checks](receipt-integer.json) in 724.77 seconds, including every
predecessor identity and the terminal full coefficient-space determinant.
This is a computational proof with a written transport argument, without
independent end-to-end review, formal verification or a novelty-priority claim.

## The map and the necessary coefficient equations

We work with exactly [Gao, arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2)
over C. On the c!=0 target chart use the monic degree-twelve eliminant

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t;Y)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

Precisely its simple roots lift to actual F7 source points. On c=0,
[GAOF7BND1](../gao-f7-fibres/PROOF.md) gives only one- and two-point fibres.
Thus any three-point fibre must lie in this chart.

Suppose, without any squarefreeness or disjointness assumptions, that

    R=Q3^3*S3,
    Q3=t^3+u*t^2+v*t+w,
    S3=t^3+z*t^2+h*t+k.

Write P=Q3^3*S3. Its coefficient of k in degree nine is exactly one.
The eliminant coefficients in degrees 11,10,9 therefore force

    Y1=-[t^11]P/4,
    Y2=([t^10]P-4Y1^2)/6,
    k=-12Y1*Y2-2-([t^9]P)|k=0.

After substituting this k, degrees eight and six force

    Y4=(6Y1+9Y2^2-[t^8]P)/2,
    Y3=([t^6]P-14Y1*Y2+6Y2*Y4-2)/2.

All divisions are by nonzero rational constants. There is no localization,
Delta denominator, genericity condition or discarded exceptional stratum.
The degree-12,11,10,9,8,6 residual coefficients vanish identically.

The six residual equations in degrees 7,3,5,4,1,2 are retained, as
primitive integer polynomials f1,...,f6 in Q[u,v,w,z,h], in
[cubic-six-input.json](cubic-six-input.json). The degree-zero equation
is not used. The checker reconstructs each of the six equations directly
from the displayed R and P and matches it up to a nonzero rational scalar.
Thus any proposed factorization must be a common zero of these six f_i.

## Rational predecessor identities

The [certificate](cubic-six-rank-dag.json.gz) contains 734 polynomials g_j.
Six initial nodes are nonzero rational multiples of the six f_i. The
other 728 nodes have literal polynomial identities

    g_j = sum_l c_jl * u^e1*v^e2*w^e3*z^e4*h^e5 * g_parent(j,l),

where c_jl is rational, all exponents are nonnegative integers and every
parent index is strictly smaller than j. Replaying those identities and applying
induction proves that every g_j belongs to I=(f1,...,f6) in Q[u,v,w,z,h].
The format avoids expanding all expressions back to the original generators.

The checker independently reconstructs the source equations and streams
these nodes in topological order. Its integer form verifies precisely the
same rational identities after clearing nonzero constant denominators:
if g_j=G_j/d_j with G_j integral and d_j>0, choose a positive common
multiple D of d_j and all denominators of c_jl/d_parent. Then check

    -(D/d_j)*G_j
      + sum_l [D*c_jl/d_parent] * monomial_l * G_parent = 0

as an integer polynomial. All bracketed coefficients are checked integers.
This is equivalent to the rational identity, with no division by a variable
or unverified specialization. The original direct rational checker is also
retained. No Groebner-basis computation is part of either replay.

## A conclusive determinant covering the entire coefficient space

The terminal witness specifies 402 monomial multiples h_k of already
proved g_j. Clear each h_k's rational coefficient denominators by a
positive integer d_k, obtaining an integer polynomial H_k=d_k*h_k.
Every H_k remains in I over Q.

Let E be the set of **all actual coefficient monomials** appearing in
these H_k, before reduction modulo any prime. The verifier computes E
from the replayed characteristic-zero polynomials; it does not trust a
finite-field support list. In this certificate |E|=402, and E contains
the constant monomial. Form the complete integer coefficient matrix

    M[e,k] = coefficient of monomial e in H_k,
    e in E, k=1,...,402.

It is square and its determinant satisfies

    det(M) = 401892858 mod 1000000007.

The modulus is verified prime by complete trial division. FLINT and a
separate plain-Python modular elimination both compute this residue.
Since it is nonzero, the integer det(M) is nonzero. Therefore M is
invertible over Q and its columns span the entire space with basis E.
In particular the constant polynomial1 is a rational linear combination
of the H_k. Hence 1 belongs to I.

Equivalently, at a common complex zero of the f_i, every H_k would vanish.
The vector of evaluations of E would then lie in the kernel of M^T.
Invertibility would make that vector zero, contradicting its constant
coordinate1. Thus the six f_i have no common complex zero.

**The complete coefficient support is essential.** A modular unit-ideal
report alone does not prove a characteristic-zero exclusion. Here a
nonzero determinant of an explicitly square integer matrix proves that
the rational columns span every actual coefficient row. If any extra
characteristic-zero monomial were present, even with every coefficient
divisible by the modulus, the verifier would include that row. The square
shape or nonzero-determinant check would then have to succeed on the
larger actual matrix; no row is silently dropped.

The direct reconstruction of the old trace's 735th explicit unit node
is therefore unnecessary. The 734 predecessor identities plus this
full coefficient-space witness prove the same ideal-membership conclusion.
The complete replay therefore excludes every R=Q3^3*S3 globally.

## Global root-multiplicity bound and three-point consequences

For an arbitrary eliminant, write its root multiplicities as m_r.
If sum_r floor(m_r/3)>=3, choose integers e_r with
0<=e_r<=floor(m_r/3) and sum_r e_r=3. Then Q3=product_r(t-r)^e_r
is a monic cubic with Q3^3 dividing R; the quotient is a monic cubic S3.
Therefore the obstruction gives the global
bound

    sum_r floor(m_r/3) <= 2.

In particular there are at most two distinct triple-or-higher roots, and
no root has multiplicity greater than eight. These are bounds, not claims
that the extremal multiplicity eight occurs.

A three-point fibre has precisely three simple roots; their monic
product S3 is squarefree and disjoint from all nonsimple factors. The
nonsimple multiplicities sum to nine. The exhaustive lists are

    (9), (7,2), (6,3), (5,4), (5,2,2), (4,3,2),
    (3,2,2,2), (3,3,3).

The [earlier reduction](../gao-f7-three-point/PROOF.md) groups them as
(t-r)^3*T3^2*S3 or Q3^3*S3. Allowing Q3 to have repeated roots is crucial:

- (9) gives Q3=(t-r)^3;
- (6,3) gives Q3=(t-r)^2*(t-s);
- (3,3,3) gives Q3=(t-r)*(t-s)*(t-q).

Thus the full cubic-cube obstruction excludes all three lists.
The remaining lists are exactly

    (7,2), (5,4), (5,2,2), (4,3,2), (3,2,2,2).

The subsequent [sharp maximum-five proof](../gao-f7-linear-cube/HIGH_ROOT_PROOF.md)
excludes(7,2). GAOF7FIVE1 subsequently excludes both fivefold lists,
leaving only(4,3,2),(3,2,2,2). All five listed here were
covered by the linear-cube/cubic-square family, with possible collisions
and repeated roots in T3. No distinctness condition
on T3 may be inserted to discard them. Their exclusion or an exact
example is still needed to decide whether the global cardinality three
occurs. F7 surjectivity, the ten previously occurring cardinalities and
the nonproper-locus results are separate established inputs. Nothing
here establishes JC2, injectivity or properness.

## Generation and the assurance boundary

The construction uses Groebner.jl0.10.6, with a pinned Julia project.
The public [change-matrix interface](https://sumiya11.github.io/Groebner.jl/interface/)
and [Demin's 2026 F4 tracing paper](https://arxiv.org/abs/2607.06372v3)
provide computational context. Our local exporter is a retained
adaptation of the GPL-v2 matrix-conversion routine; its original source
hash and license are preserved. It emits local relations to older basis
polynomials instead of expanding global generator expressions. That
backend is untrusted by the proof checker.

Direct global change-matrix attempts over Q and modulo32003 each reached
their90-second caps without exporting an identity. The local exporter
finished each modular trace in seconds. The needed local identities were
replayed independently in FLINT over the first prime. CRT reconstruction
using128 primes then recovered240 exact rational nodes. Its gain diminished,
so no further prime-generation campaign was pursued; one61-bit pilot is
retained but is not part of that CRT state.

For the remaining nodes, an earlier rational polynomial is already known.
The trace suggests monomial multiples of those parents and the intended
output support. Select independent rows and columns modulo a prime, solve
the selected square system exactly over Q, and then verify every required
coefficient constraint and the full rational local polynomial identity.
This last replay, not the modular selection, validates the node. Clearing
each parent polynomial's denominators columnwise is an invertible rational
change of unknowns. Shared systems can use several right-hand sides;
independent solves use at most eight numerical workers. These are computer
algebra workers, not separate research agents.

All completed local identities are journaled before proceeding. An early
checkpoint serialization failure is preserved and repaired with explicit
integer exponents, atomic replacement and subsequent append-only journals.
Serial, batched, Dixon and integral-column attempts are retained with their
inputs and checkpoint records; no unfinished attempt is treated as proof.
The modular outputs are compressed losslessly with original-byte hashes.

Generation is checkpointed and bounded: modular batches120seconds with an
8GB virtual-address cap and100MB per-file output cap; local batches use
explicit240- or300-second caps. Later batched and parallel runs have a4GB
address-space cap per process; the earliest serial runs had no separate
address-space cap. After a checkpoint-write MemoryError at703 completed
nodes, the streamed runner uses an8GiB parent cap and at most40GiB
aggregate worker address space, with at most8GiB per worker. The independent
checker runs neither F4 nor CRT nor the exact matrix solvers. It rederives
the necessary equations, replays all734 predecessor identities and verifies
the full coefficient-space determinant. The integer-cleared replay passed
all761 checks in724.77seconds within its1200-second time cap. Its receipt
records6 source nodes,728 derived local identities,265257 polynomial
terms and166552 relation terms. The original direct-rational replay reached
its600-second cap after last reporting700 nodes; no complete PASS is claimed
for that alternative. The compressed certificate is about2.08GB; replay streams
its nodes with the pinned iterative JSON parser rather than loading the
entire document into memory.

This passing exact certificate does not constitute formal verification,
an independent end-to-end review, external review, author endorsement or
novelty priority. MATH_STATUS.json remains the mathematical-status authority.

The [replay environment](replay-environment.json) and
[pinned dependencies](requirements-replay.txt) record Python3.12.3,
SymPy1.14.0, python-flint0.9.0 and ijson3.6.0. On Linux, from the repository
root with those dependencies installed, replay using:

```sh
timeout 1200s .venv/bin/python research/jacobian/gao-f7-three-point-certificate/verify_rank_integer.py
```

The receipt binds the exact checker, input and compressed certificate by
SHA256. [Generation metadata](terminal-rank-generation.json) records the
artifact size2077090980 bytes and its original pre-replay candidate status;
the complete PASS receipt supersedes that historical pending status.
