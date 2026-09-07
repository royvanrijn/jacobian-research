# Additive norm carriers and the ten-section node

The generic branch calculation finds one exact shared block: in the MW16
family, **ten generic sections meet the I2 node at t=-2**, and the 120
triples drawn from those sections all acquire the common norm factor
**(t+2)^4**. It is a fourth power, so it contributes no branch point to
the quadratic norm cover. After removing it, all branch polynomials in
all five tested families are squarefree and pairwise coprime.

This is a concrete shared-divisor identity. It is present in the family
before any smooth specialization is chosen, and therefore is not an event
distinguishing its high-gain and low-gain fibres. It reduces some necessary
norm-cover genera from 11 to 9; it produces no unramified class by itself.

The broader result is a necessary-condition theorem for this fixed additive
constructor. Its unramified specializations lie on finitely many quadratic
twists of explicit curves of genus **9–41**. Simultaneous norm conditions
have larger, exactly computed genera. The twist set is defined by fixed
resultants and coefficient constants; it has not been expanded or factored.
These are **necessary incidence carriers**, not sufficient rational
solubility carriers for the actual jump classes.

## Generic inputs and exact branch calculation

The [input projection](../../artifacts/generated-results/elliptic-curves/rank_jump_additive_branch_geometry_inputs_v1.json)
contains only A(t), B(t), and the marked generic section functions for
074d9, 103b2, 11952, a1-fibration-01, and the published R17 model. The
production data come from the immutable compact atlases. No specialized
parameter, outcome rank, exceptional point or candidate score enters
the branch workers.

For every unordered triple I of section x coordinates, use the previous
[additive constructor](ADDITIVE_COLLISION_BLOCKS_HAVE_FULL_PRIVATE_RAMIFICATION.md):

\[
 q_I(X)=-3X^2+2(x_i+x_j+x_k)X
       +x_i^2+x_j^2+x_k^2-2(x_ix_j+x_ix_k+x_jx_k).
\]

The coefficients may be rational functions of t. Clear their common
denominator and remove their polynomial gcd to obtain a primitive
Q_I(t,X)∈Q[t][X]. Rational scaling does not change the norm-projected
Kummer class. Compute

\[
 N_I(t)=\operatorname{Norm}_{\mathbb Q(t)[X]/(X^3+A(t)X+B(t))}
                  Q_I(t,X),
\]

and remove rational scalar content. No new elliptic point is sought.
The exact norm is a degree-at-most-two cubic-algebra determinant.

The successful certificates prove every claimed squarefreeness and gcd
by a degree-preserving reduction at one of 1009,1013,1019,1021,1031.
Pairwise gcd one at such a prime proves gcd one over Q, hence disjoint
geometric zero sets over Qbar. Degree preservation is checked for both
polynomials; reductions that lose degree are not accepted.

The initial input adapter accepted the published polynomial sections but
did not parse the compact atlas's numerator/denominator objects. That
v1 run records four UNKNOWN interface failures and the successful
published row. The [completion protocol](ADDITIVE_BRANCH_GEOMETRY_COMPLETION_PROTOCOL.json)
fixes the interface and handles actual section poles, retaining the same
triples, primes and 60-second per-family limit. It reuses the successful
published row without rerunning it. The original outputs are preserved.

## Section poles change the norm-cover degree

The compact bases are not all integral polynomial section bases. Generic
sections in 103b2, 11952 and MW16 have finite poles. Clearing these poles
in Q_I increases the degree of its norm. Assuming degree24 for every
compact-family triple would therefore give incorrect genera.

After the node correction described below, write H_I for the primitive
squarefree branch polynomial. All its degrees are even, so infinity is
unramified in the double cover z²=d H_I(t), for any nonzero rational d.
The smooth projective curve has genus (deg H_I−2)/2.

| Family | Norms | Branch degree: number of triples | Single-cover genus range |
|---|---:|---|---:|
| 074d9 | 680 | 24:680 | 11 |
| 103b2 | 680 | 24:220; 36:330; 48:120; 60:10 | 11–29 |
| 11952 | 680 | 24:10; 36:100; 48:245; 60:220; 72:95; 84:10 | 11–41 |
| MW16 a1-01 | 560 | 20:84; 24:371; 32:36; 36:69 | 9–17 |
| Published R17 | 680 | 24:680 | 11 |

All **3280** residual branch polynomials are squarefree and coprime to
the corresponding elliptic discriminant. All **1,079,960** pairs within
the same family are coprime. This statement does not compare polynomials
from different families.

An independent post-computation transport checks all thirteen fresh
frozen high/low fibres: the specialized family equation is isomorphic to
the masked short equation, and its full set of generic x coordinates
matches the frozen generic x coordinates after scaling and permutation.
The parameters are used only for this transport check, after constructing
all branch certificates. No rank label affects it. The three historical
rows have no parameter in that manifest and are not claimed as new
transport checks here.

## The common fourth power comes from the I2 node

For a1-fibration-01, gcd(Δ,Δ') has the single root t=-2. At that value,
c4 is nonzero and Δ has order two, giving a multiplicative I2 fibre.
Its cubic is

    f_-2(X) = (X-e)^2*(X+2e),       e != 0.

The exact node coordinate and all sixteen section values are in the
[node certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_additive_branch_evidence_v1.zip).
The zero-based generic indices passing through (e,0) are

    0, 1, 2, 5, 6, 9, 11, 12, 13, 15.

They give exactly C(10,3)=120 triples. These are precisely the triples
whose norm polynomial shares a factor with Δ or is not squarefree.
The 7140 unresolved pairwise gcds in the first rational-section run are
exactly C(120,2), the pairs inside this subset.

For each of those 120 norms the exact multiplicity at t=-2 is **four**;
for every other norm it is zero. Thus

\[
 N_I(t)=c_I(t+2)^4 H_I(t)\quad\text{on those 120 triples},
\tag{1}
\]

with c_I∈Q×, and N_I=c_I H_I elsewhere. The independent verifier checks
the multiplicity by vanishing of derivatives of orders zero through
three and nonvanishing of the fourth derivative. It checks the full
polynomial identity (1), not just the order of vanishing.

All 560 residual H_I are squarefree, coprime to Δ, and mutually coprime.
No shared odd branch divisor remains. The 120 factors split as 84 norms
of degree24 and 36 of degree36; removing the fourth power gives branch
degrees20 and32 respectively. This accounts for the genus drop exactly.

The calculation identifies the common singular-fibre source and the
divisor multiplicity. It does not infer a new ideal class from the
fourth power, nor does it equate this norm calculation with the
full ramification of the cubic quadratic extensions. The actual
specializations previously tested still have full private odd
ramification in the candidate span.

## Exact simultaneous norm-cover genus

Choose k different branch polynomials H_1,...,H_k and fixed rational
twists d_1,...,d_k. Their squareclasses over Qbar(t) are independent:
each H_i has simple zeros belonging to none of the other H_j. Hence
the normalization of

    z_i^2 = d_i*H_i(t),       i=1,...,k,

is a geometrically connected cover of P1 of degree 2^k. There are
sum_i deg H_i branch points, each with inertia order two; infinity
is unramified. Riemann–Hurwitz gives

\[
 \boxed{g=1+2^{k-2}\left(\sum_i\deg H_i-4\right).}
\tag{2}
\]

| Family | One constraint | Two constraints | Three constraints | Four constraints |
|---|---:|---:|---:|---:|
| 074d9 / published R17 | 11 | 45 | 137 | 369 |
| 103b2 | 11–29 | 45–117 | 137–353 | 369–945 |
| 11952 | 11–41 | 45–165 | 137–497 | 369–1329 |
| MW16 a1-01 | 9–17 | 37–69 | 113–209 | 305–561 |

The ranges are attained by selecting the indicated number of smallest
or largest branch degrees, not inferred from an average degree.

This is also a minimal-genus statement **for these specified simultaneous
norm constraints with the same parameter t**. Any smooth curve mapping
nonconstantly to all the double covers with compatible t-maps factors
through their normalized fibre product. Since its genus is at least two,
Riemann–Hurwitz makes the source genus at least (2). In particular there
is no rational or genus-one common parametrization of those constraints.
This does not claim minimality among all possible constructions of an
exceptional Mordell–Weil block.

## A finite family of necessary conditions on t

There is a useful uniform consequence of the disjoint branch calculation.
Fix one of the five families and its finite list of primitive generators
Q_I. There is a fixed finite set of rational primes T, independent of t0,
with this property:

\[
 \boxed{\prod_I [N(Q_I(t_0))Q_I(t_0)]^{c_I}
       \text{ unramified at odd good primes},\quad c_I=1
       \ \Longrightarrow\ [H_I(t_0)]\in\mathbb Q(T,2).}
\tag{3}
\]

Here Q(T,2) denotes rational squareclasses with even valuation outside T.
Exclude the finitely many singular fibres, zero norm values and undefined
section specializations. Adding a generic Selmer word to the product
does not change (3), because that word is unramified at the tested good
primes.

The set T can be defined from the following fixed nonzero constants:

- two, coefficient denominators and rational norm-normalization factors;
- leading coefficients of the branch polynomials;
- Res(H_I,Δ) and Res(H_I,H_J) for distinct triples;
- constants obtained by clearing denominators in Bezout identities for
  the three coefficients of each primitive Q_I.

The last constants exist because those three polynomials have gcd one
over Q[t]. Include the primes dividing all these constants. No assertion
that this set is small is made. The enormous product of resultants and
its prime factorization are **not computed**; their nonvanishing follows
from the exact coprimality certificates. The finite twist family is
defined by the recorded polynomials and these operations, not enumerated.

To prove (3), write t0=m/n in lowest terms and suppose p∉T has odd
valuation in H_I(t0). Since deg H_I is even and its leading coefficient
is a p-unit, p cannot divide n. The coefficient Bezout identity makes
Q_I(t0) primitive over Zp. The nonzero resultants ensure that p divides
neither Δ(t0) nor another H_J(t0). The node fourth power, where present,
does not alter parity; its zero is already a discriminant zero.

Thus p has odd valuation in N(Q_I(t0)), while the other norm values
are p-units. The primitive-degree-two lemma from the
[private-ramification proof](ADDITIVE_COLLISION_BLOCKS_HAVE_FULL_PRIVATE_RAMIFICATION.md)
gives a nonzero valuation-parity coordinate for this candidate and zero
for every other candidate. It cannot cancel in a word containing I.
This contradiction proves (3).

Consequently each participating I requires a rational point with the
same t0 on one of finitely many curves

\[
                 z_I^2=d_I H_I(t),\qquad d_I\in\mathbb Q(T,2).
\tag{4}
\]

Their genera and the simultaneous fibre-product genera are the ones
above. This is a genuine necessary arithmetic condition on t. Passing
it still leaves ideal valuation parity and the complete elliptic local
conditions unchecked, and says nothing about rational-versus-Sha
solubility of a resulting Selmer class.

Since every single curve in (4) has genus at least nine, Faltings's
finiteness theorem implies that only finitely many rational t can support
any nonempty unramified word from this **fixed finite constructor**.
See [Faltings (1983), Satz 7, p.365](https://pazuki.perso.math.cnrs.fr/index_fichiers/Faltings83.pdf).
This statement gives no effective list or height bound, and is not a
finiteness claim about rank jumps in the elliptic family. Other class
generators are not restricted to this list.

## What the result changes

The node gives an actual common-divisor structure involving 120 covers.
Its fourth-power multiplicity explains a modest genus reduction, while
the residual divisors remain disjoint. This separates shared algebraic
structure from new unramified incidence: the former is present, but
does not supply the latter on the successful fibres.

The strongest remaining route is a generator construction whose odd
branch divisors satisfy relations before specialization, leaving a
small simultaneous norm carrier. Such an identity would need to survive
primitive normalization and ideal-by-ideal local tests. The current
data supply neither that constructor nor a new strict class outside
the global root pool. The implication from a necessary norm-carrier
point to a multi-dimensional unramified class block remains missing;
global rational solubility is a further implication after that.

These calculations are **incidence** gates on a specified construction.
They produce no **solubility** certificate distinguishing rational points
from Sha and no **visibility** feature or candidate score. No active
search policy is changed. The next experiment should test a specific
shared odd-divisor identity, not enlarge the radical coefficient search
or rerun the same high/low correlation panel.

## Reproduction and verification

The [manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_additive_branch_evidence_v1.json)
summarizes the final branch degrees and genera. The
[archive](../../artifacts/generated-results/elliptic-curves/rank_jump_additive_branch_evidence_v1.zip)
preserves v1, v2 and the node certificate byte for byte, including every
norm polynomial and modular gcd witness. The initial four input failures
and the intermediate 120 unresolved MW16 norms are retained there.

The [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_additive_branch_geometry_verification_v1.json)
reconstructs primitive quadratics from the generic rational functions.
For each it bounds the norm degree and verifies the norm polynomial
up to rational scalar at enough exact integer t-values to prove the
identity. These are polynomial identity checks, not curve or point
searches. Its six-term integer determinant expansion is independent
of the worker's symbolic Sage matrix determinant.

All **110,620** interpolation evaluations, **3280** norm identities,
**1,079,960** branch-pair checks, **120** fourth-power identities and
**13** fresh generic-model transports passed. The verifier first
completed the arithmetic but encountered a Sage-integer JSON serialization
error; the incomplete JSON is preserved locally, and explicit integer
conversion was fixed before writing the valid immutable certificate.

```sh
sage -python elliptic-curves/rank-jump/verify_additive_branch_geometry.py check
```

The verifier reads the archive directly. To restore the original JSONs,
run `package_additive_branch_evidence.py unpack`. The construction scripts
and all three protocols are retained; workers have 60-second limits and
checkpoints. No integer factorization, class-group computation, new
rational-point search or new prospective parameter is used. Only new
rank-jump-specific files are committed.
