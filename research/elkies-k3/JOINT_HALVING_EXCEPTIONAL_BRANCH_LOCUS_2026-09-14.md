# Simultaneous branch conditions have a finite exceptional locus

**Written theorem with a finite-premise checker; no positive two-gain
construction.** On a fixed elliptic K3 surface with24 fibres of type I1,
two independent generic halving conditions can hold at only finitely many
closed base points of degree at most4 over a fixed number field.
The same finiteness holds for one primitive generic halving condition
together with nonzero residue-field2-torsion. These are arithmetic finite
sets, not necessarily empty or effectively enumerated sets.

For Q80, this makes the
[branch-specialization condition](Q80_BRANCH_TRACE_SPECIALIZATION_GATE_2026-09-14.md)
more precise. Within the retained good131 branch-reduction scope:

- there are only finitely many genus0 quadratic extensions with gain at
  least2, including their literal constant squareclasses;
- a genus1 cover with gain at least2 either belongs to the finite collection
  with at least two realized trace parities, or has a branch point in a
  finite exceptional half-and-torsion locus.

The [required construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open. Finiteness does not give the exceptional points, a covering
base with infinitely many rational points, or global trace-torsor sections.

## 1. Hypotheses and incidence curves

Let k be a number field and E/k(t) the generic fibre of an elliptic K3
surface whose geometric singular fibres are exactly24 I1 fibres. Assume
all of the full geometric Mordell--Weil group M=E(kbar(t)) is defined over
k(t). It is finitely generated and torsion-free. Fix T,U in M whose
classes in M/2M are linearly independent, and let V=E[2].

Over the smooth part of the base, consider the finite etale schemes

```
D_T       = {R : 2R=T},
J_(T,U)   = {(R,S) : 2R=T, 2S=U},
J_(T,2)   = {(R,e) : 2R=T, 0!=e in E[2]}.
```

Use their smooth projective normalizations. The last scheme needs only
T nonzero in M/2M. The results are

| Curve | Degree over P1 | Geometric genus | Points of degree <=4 over k |
|---|---:|---:|---|
| D_T | 4 | 9 | No finiteness claim here |
| J_(T,U) | 16 | 57 | Finite |
| J_(T,2) | 12 | 49 | Finite |

The hypotheses hold for the retained
[direct11952 Q80 parent](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
Its geometric/rational saturated basis, not an arbitrarily chosen subgroup,
is essential. Its irreducible degree24 discriminant also means that a
k=Q branch divisor of degree2 or4 cannot include a singular fibre.

One affine chart of D_T, for T=(a,b), is the trace quartic

```
H_T(m)=m^4-6*a*m^2-8*b*m-3*a^2-4*A=0.
```

Thus J_(T,U) has the fibre-product presentation H_T(m)=H_U(n)=0;
J_(T,2) has H_T(m)=0 and e^3+A*e+B=0. These charts need normalization
and completion. At poles or at specializations T=O, use the group-scheme
definition above rather than assuming the affine chart stays valid.

## 2. Full monodromy and connectedness

The single-division proof is the same24-I1 argument already written for
[determinant1092](../elliptic-curves/notes/DET1092_ALL_PRIME_DIVISION_AND_MULTISECTION_THEOREM_2026-09-09.md).
Here its relevant steps are included to specify the hypotheses.

Geometric mod2 monodromy is S3=GL2(F2). A fixed nonzero vector would give
a rational2-isogeny and a factorization of j through X0(2). The width-two
cusp would pull back to a pole of order at least2, contradicting the24
simple j-poles. The local monodromy at I1 is a transvection. A subgroup
of S3 containing a transposition and without a fixed nonzero vector is S3.

For T primitive modulo2M, the four-half torsor is connected: an orbit of
size1 would be a generic half; an orbit of size2 would give a generic
nonzero2-torsion point by taking the difference of its two halves.
Every nontrivial partition of4 has a part of size1 or2. Its affine
monodromy is therefore the full V semidirect S3=S4: its invariant
translation kernel is zero or V, and a group of order6 cannot act
transitively on4 points.

For the joint torsor, let G be its affine monodromy inside
(V direct_sum V) semidirect S3 and W its translation kernel. The invariant
subspaces of V direct_sum V are zero, the three planes V direct_sum0,
0 direct_sum V and the diagonal, and the whole space. This elementary
classification is independently enumerated in the finite checker.
If W were proper, a nonzero linear combination of its two coordinates
would vanish on W. The affine monodromy of the corresponding half torsor
for T, U or T+U would then have trivial translation kernel, contradicting
the preceding paragraph. Hence

```
G=(V direct_sum V) semidirect S3, of order96.
```

It acts transitively on16 joint halves. The stabilizer of a fixed pair is
S3. The stabilizer of its first half is V semidirect S3=S4; its S3 point
stabilizer is maximal. Consequently J_(T,U) -> D_T has degree4 and no
proper intermediate curve, even geometrically.

For J_(T,2), full affine S4 acts transitively on the12 pairs (half,
nonzero torsion vector). Forgetting the torsion vector gives a degree3
map to D_T, so this map also has no proper intermediate curve.

## 3. Local cycle counts and genus

At an I1 fibre each section has a local half over the strictly henselian
trait: it reduces into the smooth multiplicative group, where doubling
is etale and surjective over the algebraically closed residue field.
The same holds simultaneously for T and U. Thus local affine monodromy
can be translated to the diagonal linear transvection. There is no
ramification at smooth fibres.

On V it has2 fixed points and one2-cycle. On V direct_sum V it has4
fixed points and six2-cycles. On V times (V minus0) it has2 fixed points
and five2-cycles. Riemann--Hurwitz gives

```
g(D_T)     = 1-4 +24*1/2 = 9,
g(J_(T,U)) = 1-16+24*6/2 = 57,
g(J_(T,2)) = 1-12+24*5/2 = 49.
```

These are genera of connected curves; disconnected fibre products cannot
be substituted into this computation. The checker verifies transitivity,
the invariant subspaces, the maximal stabilizer and all inertia cycles.

## 4. No maps of degree at most4 to genus0 or1

Use the Castelnuovo--Severi inequality for two maps that generate the
source function field:

```
g(X) <= n*g(D_T)+d*g(Y)+(n-1)*(d-1),
```

where X -> D_T has degree n and X -> Y has degree d. The standard
function-field generation condition is crucial; see the proof and
statement in [Primitive algebraic points on curves, Theorem14](https://link.springer.com/article/10.1007/s40993-024-00543-4).

For X=J_(T,U), n=4. If d<=4 and g(Y)<=1, the right side is at most49,
strictly below57. If the maps do not generate the full field, primitivity
of X -> D_T forces kbar(Y) to be contained in kbar(D_T). Then d is a
multiple of4; the only possibility d=4 would identify D_T birationally
with Y, contradicting genus9.

For X=J_(T,2), n=3. The bound is at most37, strictly below49.
The same primitivity argument leaves only d=3 and a birational map
D_T -> Y, again impossible. These exclusions are geometric, so they
persist over every finite extension of k.

## 5. Finiteness of low-degree points

This step uses a published theorem, not a height search. By
[Kadets--Vogt, Theorems1.2(1) and1.3](https://arxiv.org/pdf/2208.01067),
a curve of genus greater than7 with infinitely many points of degree at
most4 must admit a map of degree at most4 to a genus0 or genus1 curve.
For clarity, this consequence follows by taking the smallest degree d
with infinitely many points: Theorem1.3 gives genus bounds1,2,4,7 for
d=1,2,3,4 unless the points descend along a nontrivial cover. For d<=3
that cover goes directly to a curve with infinitely many rational points.
For d=4 its only other possibility is a double cover of a curve with
infinitely many quadratic points; Theorem1.2(1) gives a second double
cover to genus0 or1. Section4 excludes every resulting map.

Thus both incidence curves have finitely many degree-at-most4 points
over each fixed number field. This does not assert finiteness when the
number field is allowed to vary without fixing it.

## 6. Consequences for common branch halving

For a smooth-branch quadratic cover C -> P1, any trace Q+sigma(Q)
specializes at each closed branch value b to twice a k(b)-rational point.
Choose fixed representatives of M/2M. Subtracting an inherited section
from Q changes its trace by twice that section, so every realized parity
can be represented by one of those fixed representatives.

If two independent trace parities are realized, every branch value lifts
to J_(T,U) over exactly its residue field. When g(C)<=1 its degree is at
most4. There are finitely many parity pairs, and each joint curve has only
finitely many such points. Therefore **all branch values of all these
covers lie in one finite set**, for the fixed parent and number field.
It follows that their branch divisors form a finite set.

There are also only finitely many literal quadratic extensions in this
case. Fix one branch divisor and a representative d(t). Other quadratic
covers with that divisor are w²=c*d(t), c in k*/k*². Over kbar their
pullback elliptic surfaces are isomorphic. The anti-invariant geometric
Mordell--Weil group has finite rank and is finitely generated, since the
elliptic curve remains nonisotrivial. Its generators are defined over a
finite constant extension. Twisting c twists this Galois representation
by the corresponding quadratic character. A nonzero invariant direction
can occur for only finitely many such characters. Two realized trace
parities require a nonzero anti-invariant group, so only finitely many c
can occur. This finiteness argument provides no list or effective bound.

Likewise, if one nonzero trace parity is realized and an anti-invariant
section has a nonzero2-torsion value at a branch b, that branch lifts to
J_(T,2). Such marked branch values form a finite set. This does **not**
put every other branch value in that finite set.

For Q80 use the retained [exact rank formula and good131 bounds](Q80_GOOD_BRANCH_CODE_AND_TRACE_NORMS_2026-09-14.md):

```
gain = dim(realized trace parities)+dim(actual branch image),
genus0 good131: branch dimension0,
genus1 good131: branch dimension at most1.
```

A genus0 gain of at least2 therefore belongs to the finite two-trace
collection. In genus1, gain at least2 either has at least two traces,
or has exactly one trace and a nonzero actual branch image. In the latter
case at least one branch lies in the finite half-and-torsion locus.
The good131 hypothesis is retained exactly; bad branch reduction does
not inherit these numerical branch-image bounds.

## 7. Evidence and next construction gate

The [finite checker](scripts/verify_joint_halving_incidence.py) and its
[certificate](../artifacts/generated-results/elkies-k3-joint-halving-incidence-v1/result.json)
verify the tiny group actions and integer genus/inequality premises in
under0.01 seconds, with10 CPU seconds and512MiB allowed. They do not
verify the published diophantine theorem, the entire written proof, or an
exceptional-point list. No formal or external review is claimed.

```
python3 research/elkies-k3/scripts/verify_joint_halving_incidence.py --record /tmp/joint-halving-incidence-replay.json
```

The next constructive step is to obtain an actual low-degree point on
one of these incidence curves from generic-only data, then test global
trace-torsor solubility on the resulting common divisor. Neither a
modular incidence nor this noneffective finiteness theorem supplies that
point. The known exact rank18 control supplies only one half condition;
it does not populate the two-condition exceptional locus.
