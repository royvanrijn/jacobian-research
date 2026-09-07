# The remaining genus-two factor has no elliptic quotient

For the fixed successful native triple `01333,0b2d0,19e45`, the
genus-two descent base

\[
H:\ y^2=h(t)f(t)g(t)
\]

has an **absolutely simple Jacobian**. In particular, H admits no
nonconstant map to an elliptic curve, even after a finite extension of
Q. This closes the possibility that the residual genus-two factor in
the [simultaneous-lift decomposition](GENUS_TWO_DESCENT_ORGANIZES_THE_NATIVE_TRIPLE_LIFT.md)
is secretly another product of elliptic curves.

The rational rank of J(H) is still **UNKNOWN**. A single bounded request
for an unconditional 2-Selmer computation received an explicit server
offline response before mathematical execution. No rank bound, Selmer
dimension or Chabauty conclusion is inferred from that attempt.

## Exact Frobenius certificate

The sextic has degree six and is squarefree modulo 131. Exact counts are

\[
\#H(\mathbf F_{131})=142,\qquad
\#H(\mathbf F_{131^2})=17330.
\]

The second count is computed twice: using Sage's finite-field arithmetic,
and using integer pairs in F131[w]/(w²−d), where d is the first positive
quadratic nonresidue. The second method evaluates the sextic by pair
arithmetic and tests its norm in F131. It includes both points at
infinity over the quadratic extension.

These counts give the Frobenius polynomial

\[
\boxed{P(X)=X^4+10X^3+134X^2+1310X+17161.}
\]

It is irreducible over Q. Irreducibility alone proves simplicity over
F131, but does not suffice for absolute simplicity. We therefore check
all possible roots of unity among ratios of distinct Frobenius roots.

Let α1,…,α4 be the roots of P. The monic resultant

\[
\operatorname{Res}_X(P(X),P(zX))
\]

has roots αj/αi. Dividing its monic normalization by (z−1)^4 gives a
degree-twelve polynomial R. The calculation verifies R(1)≠0 and

\[
\gcd(R,\Phi_n)=1
\quad\text{for every }n>1\text{ with }\varphi(n)\leq12.
\]

This is a finite exhaustive test, not a numerical angle calculation.
For each prime power p^e dividing n,
\(\varphi(p^e)^2/p^e=p^{e-2}(p-1)^2\); every factor is at least one
except the factor 1/2 for 2^1. Thus φ(n)²≥n/2, and φ(n)≤12 implies
n≤288. The 25 tested orders are

```text
2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,21,22,24,26,28,30,36,42.
```

No distinct αi,αj can therefore satisfy αi^m=αj^m for any m>0.
The four conjugates of α1^m remain distinct, so Q(α1^m)=Q(α1) and
the Frobenius polynomial over every finite extension remains irreducible.
A proper abelian subvariety over such an extension would give a proper
factor of that polynomial. Hence the reduction is absolutely simple.
This is the Frobenius-field approach to absolute simplicity used by
[Howe and Zhu](https://arxiv.org/abs/math/0002205); here the required
no-collision condition is certified by an explicit resultant.

Finally, a nontrivial decomposition of J(H) over Qbar would be defined
over a finite number field. At a place above this good prime its abelian
factors would specialize to a nontrivial decomposition of the reduction,
contradicting absolute simplicity. Thus J(H) is absolutely simple in
characteristic zero as well. A nonconstant map H→E would induce a
nonzero quotient J(H)→E, which is now excluded.

## What the bounded descent did and did not establish

The frozen program constructs H, verifies its retained rational point,
checks rational 2-torsion order four, and requests `TwoSelmerGroup(J)`.
It contains no `RankBounds`, `RationalPoints`, prospective parameter
enumeration or non-rigorous class-group setting. A completed full Selmer
group of order 2^s would give rank J(H)(Q)≤s−2. The relevant semantics
are in the [Magma handbook](https://docs.magma-maths.org/ArithmeticGeometry/HyperellipticCurves/jacobians-number-fields.html).

The repository's existing calculator transport was imported unchanged.
The entire server response was

```xml
<calculator><offline>The Magma calculator is temporarily disabled due to electrical work in the building.</offline></calculator>
```

The wrapper's successful HTTP execution is therefore distinct from
mathematical completion. The fail-closed parser returns UNKNOWN because
there is no Selmer result or completion marker. The response, exact
program, source hashes and local execution record are retained. No second
request was made and no remote computation is being treated as live.

The preflight also showed a potential future cost: the h branch field
is Q(sqrt(-128751105131280719677590)). A direct full descent can require
substantial class-group work. This is a cost observation, not evidence
that the unexecuted descent would time out or that the rank is large.

## Consequence for the solubility mechanism

The exact global condition remains

\[
i(H(\mathbf Q))\cap\psi(A(\mathbf Q)),
\]

with the specific degree-four isogeny ψ from the three quadratic factors.
The genus-two curve-point condition cannot be replaced by an elliptic
quotient of H. In the genus-five Jacobian decomposition, the three
elliptic factors have ranks 3,2,3, while this fourth factor is a simple
abelian surface with unknown rational rank. Information from the three
elliptic descents does not determine the missing surface component.

This is a **solubility-method exclusion**, not an original-fibre rank
prediction. It does not rule out covering methods: the genus-five
unramified cover already has elliptic quotients, and further covers may
be useful. Nor does it exclude the original varying elliptic-surface
point maps; those are not maps from H to a single fixed elliptic curve.

The highest-value unresolved calculation is still a certified rank bound
or a descent on the particular genus-two isogeny class. The frozen
program can be executed when a working Magma installation is available.
A direct (2,2)-isogeny descent over rational squareclasses is another
possible route; its local images and global Selmer quotient have not yet
been computed here. No additional genus-five elliptic-factor search is
justified by this result.

## Replay

The protocol is [GENUS_TWO_RANK_AND_SIMPLICITY_PROTOCOL.json](GENUS_TWO_RANK_AND_SIMPLICITY_PROTOCOL.json).
Certificates are `rank_jump_genus_two_absolute_simplicity_v1.json` and
`rank_jump_genus_two_rank_probe_v1.json`, with the latter's frozen input,
under `artifacts/generated-results/elliptic-curves/`.

```sh
sage -python elliptic-curves/rank-jump/genus_two_absolute_simplicity.py check
python3 elliptic-curves/rank-jump/genus_two_rank_probe.py check
```

The second replay is entirely offline and verifies the recorded UNKNOWN
outcome. Neither replay submits a calculator request.
