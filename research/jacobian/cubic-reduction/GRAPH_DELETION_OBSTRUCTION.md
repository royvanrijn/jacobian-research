# Arbitrary-degree elementary target obstructions for nine fixed graph deletions

Let \(K\) be the [integral36-term twelve-variable cubic](PROOF12.md), KCUB3.
For each \(j=4,\ldots,12\), its carrier coordinate has the form
\(K_j=x_j+g_j\), where \(g_j\) is independent of \(x_j\).
Put \(x_j=-g_j\), delete that source coordinate and delete output \(K_j\).
Write the eleven retained raw outputs as \(G^{(j)}_a\), keeping their old
indices.

**Claim KCUB4.** For each of these nine fixed slices, one retained raw
component \(G_i^{(j)}\) cannot be made degree at most three by

\[
G_i^{(j)}\longmapsto
G_i^{(j)}-P\bigl(G_a^{(j)}:a\ne i,j\bigr)
\]

for any polynomial \(P\) of any degree over any characteristic-zero field.

Thus increasing only the target correction degree cannot turn one of
these raw graph slices into a cubic map. The claim does not exclude
different source coordinates, different graph levels, prior transformations
that mix the obstructed component into other outputs, or general sequences
of target changes. It proves no minimum cubic dimension.

Both known colliding points lie on every carrier level \(K_j=0\), so
these are genuine candidate dimension reductions. The source graph
straightening has determinant one; the raw eleven-variable restriction
is Keller and inherits the collision. The unresolved question was whether
its higher-degree coordinates could be repaired into cubics.

## A finite functional can obstruct every target degree

For each slice, the certificate supplies a rational coefficient functional

\[
\Lambda(f)=\sum_{e\in E}\lambda_e[x^e]f
\]

with finite support, satisfying

\[
\Lambda(G_i^{(j)})=1,\qquad
\Lambda(f)=0\quad(\deg f\le3),\qquad
\Lambda(P(G_a^{(j)}:a\ne i,j))=0
\]

for every polynomial \(P\). These identities immediately prove the claim.

The last condition is finite to verify. Let \(T\) be the set of all
monomial divisors of the certificate support \(E\). Its complement
generates an upward-closed monomial ideal \(I_T\). Projection to the
finite algebra \(\mathbb Q[X]/I_T\) preserves every coefficient sampled
by \(\Lambda\).

Every retained raw output vanishes at the origin. If
\(D=\max_{e\in T}|e|\), a product of more than \(D\) such outputs is zero
in this quotient, because its source order exceeds \(D\). It is therefore
enough to check every commutative target monomial of total degree at
most \(D\), together with the constant. This is a bound on the
finite certificate check, not a restriction on the allowed degree of \(P\).

The checker performs exact sparse multiplication actions in the quotient.
It enumerates nondecreasing target words, verifies that every nonzero
word is annihilated, and accounts combinatorially for every pruned zero
branch. The covered number is

\[
\sum_{r=1}^{D}\binom{10+r-1}{r}=\binom{10+D}{D}-1.
\]

Words of higher degree vanish by positive source order. Rational identities
extend to every characteristic-zero coefficient field.

## Certified cases

| Deleted carrier \(j\) | Obstructed raw component \(i\) | Certificate order \(D\) | Nonzero functional weights |
|---:|---:|---:|---:|
|4|2|6|34|
|5|3|7|105|
|6|1|4|1|
|7|1|4|5|
|8|2|6|53|
|9|3|7|86|
|10|2|6|35|
|11|3|7|96|
|12|2|5|38|

For example, deleting \(K_6=x_6+x_1^2\) gives

\[
G_1=x_1-3x_1^2x_2-x_1^3x_3.
\]

The certificate is simply \(\Lambda=-[x_1^3x_3]\).
It evaluates to one on \(G_1\) and annihilates every polynomial in the
other ten raw outputs. The other eight certificates are retained
explicitly in the case files; none relies on a numerical rank tolerance.

## Production and independent verification

The initial [degree-three target screen](graph36.py) uses exact SymPy
polynomials and FLINT rational RREF, retaining dual certificates.
[Quartic own-variable jets](graph36_jet.py) already close two slices
against all target degrees.

The [adaptive producer](graph36_adaptive_wide.py) enlarges small
downward-closed monomial quotients until it separates the higher-degree
defect from the raw-output algebra. All nine cases terminate with exact
separators. The source support orders are at most seven.

The [independent checker](verify_graph36.py) uses Python standard-library
Fractions and a separate quotient multiplication-action implementation.
It does not use SymPy, FLINT, any rank assertion or the producer's
intermediate membership computations. Starting from the pinned KCUB3
map, it reconstructs each graph, checks the origin jets and the raw
generator set, verifies the nonzero defect, and checks the complete
target-word coverage for the final functional.

- [Independent receipt](graph36-validation.json)
- [Scoped integration checks](graph36-integration.json)
- [Primary complete manifest](adaptive36-wide/checkpoint.json)
- [Case4](adaptive36-wide/pivot-04.json), [case5](adaptive36-wide/pivot-05.json),
  [case6](adaptive36-wide/pivot-06.json), [case7](adaptive36-wide/pivot-07.json),
  [case8](adaptive36-wide/pivot-08.json), [case9](adaptive36-wide/pivot-09.json),
  [case10](adaptive36-wide/pivot-10.json), [case11](adaptive36-wide/pivot-11.json),
  [case12](adaptive36-wide/pivot-12.json)

The checker also replays the successful BCR2 square completion as a
positive control:

\[
(x_4+x_1^2x_2^2)-(x_8+x_1x_2)^2
=x_4-2x_1x_2x_8-x_8^2.
\]

A functional that sees the degree-four defect cannot annihilate this
target square. This guards the distinction between a genuine dual
obstruction and an incomplete target basis.

From the repository root:

~~~sh
.venv/bin/python research/jacobian/cubic-reduction/graph36.py
.venv/bin/python research/jacobian/cubic-reduction/graph36_jet.py
.venv/bin/python research/jacobian/cubic-reduction/graph36_adaptive_wide.py
.venv/bin/python research/jacobian/cubic-reduction/verify_graph36.py
.venv/bin/python research/jacobian/cubic-reduction/check_integration_graph36.py
~~~

The earlier finite quotients and capped adaptive run are retained in
all-degree36 and adaptive36. Their inconclusive statuses are not proof
premises. The final result uses only the explicit separating functionals,
the independently checked finite coverage, and the positive-order argument.

## Consequence for the next construction attempt

Any sequence of changes confined to the other ten raw outputs preserves
their generated polynomial algebra. Such changes cannot defeat the
certificate for \(G_i\). A successful continuation must leave that algebra:
for example, alter the source graph or mix \(G_i\) into another output
before repairing it.

A separate [raw target-relation search](stages36.py) through target degree
three finds no new cubic-pullback kernel direction involving the obstructed
component in these frames. Its [checkpoint](stages36/checkpoint.json) is
exploratory; it is not an independently certified restriction on general
target sequences or intermediate degrees.

The next attempt may therefore need temporary degree growth, another
source coordinate, a different carrier realization or another construction.
The current positive bounds remain dimension twelve with36 displayed
cubic nonlinear occurrences, and [10,48] in the separate normalized
three-variable degree-seven domain. The global simplicity objective is open.

