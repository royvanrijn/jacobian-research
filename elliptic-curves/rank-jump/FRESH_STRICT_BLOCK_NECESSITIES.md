# How much of a recorded jump must be strict?

The frozen sixteen-fibre panel now has complete equation-defined local
boundary bounds on its nine fully factored rows. Four fresh +10 fibres
require at least **9, 5, 6, 8 additional strict rational directions**;
historic curve 398's +14 requires at least **5**. These are necessities
obtained by joining existing rank lower bounds *after* the masked
arithmetic. They are not independent measurements of additional classes.

The distinction matters in the controls: three observed-zero fibres
already possess inherited strict rational subspaces of dimensions 9, 10
and 5. A large soluble strict subspace can belong to the generic subgroup.
The relevant incidence quantity is the **excess beyond that subgroup**,
not the total strict dimension or the degree of a generic governing field.

This extends the [masked governing/CT panel](FRESH_RANK27_GOVERNING_AND_CT_COMPARISON.md).
The full additional-class CT comparison remains UNKNOWN. Seven rows,
including the fresh +11 and both historic +12 controls, still lack complete
factorization in the frozen inputs and receive no complete boundary bound.
The [comparison certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_strict_boundary_coordinate_comparison_v1.json)
retains all sixteen rows and their missing values.

## Accounting before and after the rank labels

Let G be the marked generic Kummer subgroup of dimension m, U the strict
Selmer group, and k=dim(G∩U). Here strict means locally square above 2,
infinity and every bad rational prime, and unramified elsewhere in the
cubic field K. The retained [strict-class identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md)
gives

\[
 c_S=\dim U=\dim\operatorname{Cl}(\mathcal O_{K,S_K})/2.
\]

Let I be the localized full Selmer image and let ell be the dimension
of the product of the local point images at S. The derivative witness
below proves dim I≤h=ell−1 on each completed row. Generic sections give
dim loc(G)=g=m−k. Thus, with a=h−g,

\[
 0\to U/(G\cap U)\to\operatorname{Sel}_2(E)/G
 \to I/\operatorname{loc}(G)\to0,
\]
\[
 \boxed{\dim\operatorname{Sel}_2(E)/G=(c_S-k)+e,
 \qquad0\le e\le a.}
\tag{1}
\]

All of m,k,g,ell,h,a are computed before reading the rank labels.
The exact value of c_S, and generally e, is not computed. With a retained
rank lower bound R, no rational 2-torsion, and recorded gain J=R−m,
linear algebra then forces

\[
 \boxed{\dim (W\cap U)/(G\cap U)\ge\max(0,J-a),}
\tag{2}
\]

where W is the full rational Kummer image. The total strict rational
dimension, and hence c_S, is at least max(k,R−h). The cubic irreducibility
and generic mod-two independence are already certified in the masked panel.

For a family parameter t, this supplies a necessary incidence inequality

\[
 \operatorname{rank}E_t\ge m+J
 \quad\Longrightarrow\quad c_{S(t)}-k(t)\ge J-a(t).
\tag{3}
\]

It does not supply a sufficient condition on t for simultaneous rational
solubility. Even a computed c_S−k is a Selmer incidence quantity; the
rational subspace and Sha still have to be distinguished.

## Completed panel rows

The rank and gain columns remain lower-bound accounting. An observed zero
is censored and is not an exact-rank assertion. Each row has its own cubic
field and its own c_S.

| Case / parameter | Family | Rank ≥ / m | Gain | Inherited strict k | ell | Additional boundary cap a | Additional strict rational dimension forced by R |
|---|---|---:|---:|---:|---:|---:|---:|
| low, 2824/885 | 074d9 | 17 / 17 | 0 | 9 | 9 | 0 | ≥0 |
| new-71, 3726/881 | 103b2 | 27 / 17 | +10 | 0 | 19 | 1 | ≥9 |
| low, −1049/2296 | 103b2 | 17 / 17 | 0 | 0 | 20 | 2 | ≥0 |
| new-41, −2448/11 | 11952 | 27 / 17 | +10 | 0 | 23 | 5 | ≥5 |
| low, −1171/1683 | 11952 | 17 / 17 | 0 | 10 | 8 | 0 | ≥0 |
| new-188, 110314/102227 | 11952 | 27 / 17 | +10 | 0 | 22 | 4 | ≥6 |
| low, 130349/28916 | 11952 | 17 / 17 | 0 | 5 | 13 | 0 | ≥0 |
| new-48, 2828/2015 | 11952 | 27 / 17 | +10 | 0 | 20 | 2 | ≥8 |
| ICARM398 | recovered MW16 parent | 30 / 16 | +14 | 0 | 26 | 9 | ≥5 |

The upper cap a is not the measured boundary contribution e. In particular,
curve 398 may have more than five strict rational directions; the present
constraints simply do not force them. It would be incorrect to conclude
that its +14 has a smaller actual strict block than a +10 example.

There is a useful exact consequence on three controls. Since a=0,

\[
 \dim\operatorname{Sel}_2(E)/G=c_S-9,\quad c_S-10,\quad c_S-5
\]

for the 074d9 low, compact 11952 low, and larger 11952 low respectively.
Every additional Selmer direction there must be represented by a strict
class modulo G. Their class-group excess remains unknown. Independently
certifying c_S=9,10,5 respectively would close those controls at exact
rank 17 and Sha[2]=0. Those upper certificates have not been obtained.

## Paired lessons

1. **103b2, 3726/881 versus −1049/2296.** Both have k=0. The high fibre's
   additional boundary capacity is only one, while the control's is two.
   More available boundary dimension does not explain this successful +10.
   At least nine of its rational quotient directions must be strict.
2. **11952, −2448/11 versus −1171/1683.** The high fibre requires at least
   five additional strict directions and leaves at most five boundary
   directions. The control already has ten strict rational directions,
   all inherited, and zero capacity for an additional boundary contribution.
   Its large known soluble block is not evidence of a large jump.
3. **11952, 110314/102227 versus 130349/28916.** The larger-parameter
   comparison requires at least six additional strict directions on the
   high fibre; the low has five inherited ones. Both facts are compatible
   with unknown excess class dimensions on either curve. Total block counts
   cannot resolve the difference.
4. **11952, 2828/2015 versus −1171/1683.** A second high matched to the
   same compact control requires at least eight additional strict directions,
   despite having the same recorded gain as −2448/11. Even within one family,
   the equation-defined boundary allows different decompositions of +10.
   Reusing the control does not create an independent matched replication.

Matching and exposure limitations are unchanged from the original panel.
These are structural deductions on selected fibres, not an aggregate
prospective discrimination test.

## Exact derivative witnesses, including the complex cubic fields

For the monic cubic f, delta=disc(f), theta a root, set

\[
 \beta=-\delta f'(\theta),\qquad N\beta=\delta^4.
\]

The capture verifies even valuations at every polynomial-discriminant prime
omitted from S. At all other odd good primes beta is a unit. Thus global
reciprocity makes its local pairing vanish on I.

For three real roots, beta has signs (−,+,−), and the nontrivial real point
class has signs (+,−,−). Their nonzero Hilbert pairing certifies a nonzero
functional on the full local point product.

For one real root, that argument gives no witness. The coordinate method
uses a finite place where the generic local classes span the complete point
image L_p. This image is a maximal isotropic subspace for the perfect local
Tate pairing, expressed by the product of cubic Hilbert symbols on
norm-square classes. Therefore

\[
 \beta_p\notin L_p=L_p^\perp
 \quad\Longrightarrow\quad
 \langle\beta_p,-\rangle|_{L_p}\ne0.
\]

Local squareclass coordinates certify this nonmembership at **13** for the
074d9 control, at **2** for the compact 11952 control, and at **3** for the
larger-parameter 11952 high. This closes the missing real-place argument
without exceptional points or local point enumeration. The derivative
class itself is excluded from Selmer at the witness place; it is a
constraint, not a newly constructed soluble direction.

## Verification, failures and next priority

The [direct-Hilbert protocol](FRESH_STRICT_BOUNDARY_PROTOCOL.json) completed
three of nine eligible rows; six hit its 45-second cap. Seven further rows
were skipped for incomplete existing factorization. The
[direct certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_strict_boundary_v1.json)
preserves these failures. Direct symbols use
[PARI's documented local nfhilbert interface](https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html#nfhilbert).

The separate [coordinate protocol](FRESH_STRICT_BOUNDARY_COORDINATE_PROTOCOL.json)
uses no direct Hilbert calls and a 30-second cap per eligible row. All nine
complete. Its [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_strict_boundary_coordinates_v1.json)
binds equations, prior masked local data, scripts and limits. No failed
factorization is retried and no class-group calculation or point search runs.

The [verification](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_strict_boundary_verification_v1.json)
passes on all nine completed rows: independent rational norm identities,
PARI versus Sage bad-prime coverage, omitted-prime valuations, splitting
dimensions, matrix ranks, real witnesses and label-joining inequalities.
Finite-coordinate replay uses the same LocalSquareclasses backend; the
three completed direct-Hilbert rows agree with the coordinate bounds.
It does not pretend the seven unknown rows were verified complete.

Replay:

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_fresh_strict_boundary.py check
```

The priorities are now:

1. **Incidence:** independently compute or bound c_S−k, rather than count
   the inherited block again. Formula (3) supplies a precise necessary
   threshold once a class bound is available.
2. **Solubility:** construct independent additional classes and evaluate
   their CT obstruction; then identify the missing higher-descent or rational
   carrier criterion. Zero CT alone is not a rationality theorem.
3. **Coverage:** recover any already certified equation-only factor data for
   the seven unfinished rows before considering new factorization work.
4. **Weak explanations:** generic governing degree, inherited strict size,
   and local-boundary capacity alone fail to distinguish the observed pairs.
   The latter even orders the 103b2 pair in the opposite direction.

No visibility feature or new selection score is proposed for Agent1. The
new result locates a necessary part of each successful block; it still does
not explain what makes that block rationally soluble at the specialization.
