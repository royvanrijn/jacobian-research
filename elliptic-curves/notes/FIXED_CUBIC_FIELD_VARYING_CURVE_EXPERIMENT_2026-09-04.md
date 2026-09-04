# Fixed cubic field, varying elliptic curve

Status: **at `u=-1`, the inherited pairing has rank 16 and radical
dimension 2; rank at least one is certified outside the inherited span**.
Only the three nonzero radical classes remain inherited point-solving
candidates; their realization is open.

This is a mechanism experiment, not a replacement for the MW16/MW17
production lanes.  It fixes a cubic 2-division field and a certified space of
global Kummer classes, varies the elliptic curve, and asks which combinations
remain everywhere locally admissible.

## Fixed-field family

Let

\[
 E_0:y^2=f(x),\qquad f(x)=x^3+Ax+B,
 \qquad K=\mathbf Q(\theta),\quad f(\theta)=0.
\]

For

\[
 \alpha_u=\theta+u\theta^2
\]

the norm polynomial is

\[
\begin{aligned}
F_u(x)={}&x^3+2Au x^2+(A+3Bu+A^2u^2)x\\
          &+B+ABu^2-B^2u^3.
\end{aligned}
\]

Exact expansion gives

\[
 \operatorname{disc}(F_u)
 =\operatorname{disc}(f)(1+Au^2+Bu^3)^2.
\]

Write \(D_u=1+Au^2+Bu^3\).  The experiment excludes \(D_u=0\).  When
\(D_u\ne0\), the inverse change of cubic generator is explicit:

\[
 \theta=\frac{-2Bu^2+(1-Au^2)\alpha_u-u\alpha_u^2}{D_u}.
\]

Thus the roots are not merely roots in isomorphic cubic fields: the formula
fixes the labelled Galois-equivariant identification
\(E_0[2]\simeq E_u[2]\).

## Anchor and fixed global span

The pilot uses the pinned Fermigier--Mestre rank-at-least-20 specialization at
family parameter \(28917/20\).  Its certified short model is

\[
y^2=x^3+Ax+B
\]

with

```text
A = -5750886029903523759416717668139307
B = 167347710468055045100164888198438918505621536951206.
```

The repository already certifies that the twenty pinned points have
20-dimensional image in \(E_0(\mathbf Q)/2E_0(\mathbf Q)\).  Since \(f\) is
irreducible, \(E_0(\mathbf Q)[2]=0\), and their Kummer representatives

\[
 \beta_i=x(P_i)-\theta
\]

form a 20-dimensional subspace

\[
 W\subset H^1(\mathbf Q,E_0[2])
       \simeq H^1(\mathbf Q,E_u[2]).
\]

Every norm identity
\(N_{K/\mathbf Q}(\beta_i)=y(P_i)^2\) is replayed exactly.  This establishes
the global cohomology classes; it does not assert that they lie in the local
Kummer image for a new \(E_u\).

## Full-span local computation

The first bounded policy is all integers \(|u|\le2\).  For each \(u\), the
runner checks:

1. every finite prime dividing \(2\operatorname{disc}(f)D_u\), including all
   newly bad primes of \(E_u\);
2. the real place;
3. the restriction of all twenty \(\beta_i\) to the exact local squareclass
   space of \(K\otimes\mathbf Q_v\);
4. an explicit local-point Kummer basis for \(E_u(\mathbf Q_v)/2E_u(\mathbf
   Q_v)\).

At odd \(p\), the required local-image dimension is
\(\dim E_u(\mathbf Q_p)[2]\), and at \(p=2\) it is one larger.  The number of
local factors of the fixed cubic gives the torsion dimension.  Local points
are searched on both the raw and a \(p\)-minimal model, but a search result is
accepted only when exact local squareclass independence reaches this known
dimension.  Consequently the bounded point enumeration is a witness finder,
not a completeness assumption.  A failure to reach the dimension suppresses
the entire \(W_u\) claim.

Outside the declared support, both integral cubics have good reduction at odd
primes and the fixed classes are unramified, so the good-prime local condition
is automatic.  No class group or global Selmer basis is computed.

The local quotient rows are concatenated and one linear map

\[
 W\longrightarrow\bigoplus_v
 H^1(\mathbf Q_v,E_u[2])/
 \delta_v(E_u(\mathbf Q_v)/2E_u(\mathbf Q_v))
\]

is formed.  The reported space is its kernel.  In particular, the code does
not discard basis elements one at a time: combinations of individually
inadmissible \(\beta_i\) can and do remain eligible for the kernel.

## First result

All local Kummer images completed, and the five curves have pairwise distinct
rational \(j\)-invariants.

| \(u\) | checked finite places | newly bad primes | local-condition rank | \(\dim W_u\) |
|---:|---:|---:|---:|---:|
| \(-2\) | 17 | 5 | 7 | 13 |
| \(-1\) | 14 | 2 | 2 | 18 |
| \(0\) | 12 | 0 | 0 | 20 |
| \(1\) | 15 | 3 | 7 | 13 |
| \(2\) | 16 | 4 | 7 | 13 |

The \(u=0\) positive control recovers the whole 20-dimensional space.  The
sharpest new-curve observation is \(u=-1\): eighteen independent fixed global
classes remain everywhere locally admissible.  This is exactly the signal the
experiment was designed to isolate.

It is not a rank lower bound.  The exact sequence still allows surviving
classes to be represented by nontrivial elements of \(\Sha(E_u)[2]\).

## Point-realization covers

For a surviving span class, use a representative

\[
 \beta=\prod_{i\in I}\beta_i=b_0+b_1\theta+b_2\theta^2
\]

and write

\[
 \gamma=a+b\theta+c\theta^2.
\]

Reducing \(\beta\gamma^2\) modulo
\(\theta^3+A\theta+B\), the projective covering in variables
\((a:b:c:d)\) is the intersection of the two quadrics

\[
 [\theta](\beta\gamma^2)+d^2=0,
 \qquad
 [\theta^2](\beta\gamma^2)+u d^2=0.
\]

For \(d\ne0\), a rational point gives

\[
 x=\frac{[1](\beta\gamma^2)}{d^2},
 \qquad
 y=\frac{s_\beta N(\gamma)}{d^3},
 \qquad s_\beta^2=N(\beta).
\]

The artifact stores an exact basis of every \(W_u\), including the anchor
basis masks, the reduced power-basis coefficients of \(\beta\), and
\(s_\beta\). The point-realization experiment below uses these exact inputs.

## Cassels--Tate obstruction and the remaining candidates

The [exact pairing computation](FIXED_CUBIC_U_MINUS1_CASSELS_TATE_2026-09-05.md)
now gives rank 16 on the eighteen-dimensional surviving space. Every class
outside its two-dimensional restricted radical is provably nontrivial in
Sha, including all eighteen original basis elements. Point solving should
therefore use only the radical generators with anchor masks `1047173` and
`596921`, and their sum `450876`. All three have exact reduced-cover inputs;
their point-or-Sha status remains unknown. The complete matrix and arithmetic
witnesses are in the linked canonical note.
<!-- status-consumer: EC-FIXED-CUBIC-U-MINUS1-CASSELS-TATE df45391a84f0e3c9 -->

## A certified point outside the inherited space

At \(u=-1\), direct substitution gives

\[
 Q=(A+1,A-B+1),\qquad
 2Q=\left(\frac14,B+\frac A2+\frac18\right).
\]

The actual Kummer representatives are

\[
 \delta(Q)=\eta=\theta^2-\theta+A+1,\qquad
 \delta(2Q)=(\theta-1/2)^2,
 \qquad (\theta+1)\eta=1+A-B.
\]

For the pinned anchor,

```text
Q.x = -5750886029903523759416717668139306
Q.y = -167347710468055050851050918101962677922339205090512
```

This point is independent of the entire original twenty-dimensional global
Kummer span, not just the eighteen surviving directions. At 19 the fixed
cubic factors as

\[
 f(T)=(T+1)(T^2+18T+16)\pmod {19}.
\]

The two unramified prime ideals have residue degrees 1 and 2. The valuations
of \(\eta\) at them are respectively \((0,1)\); every one of the twenty
anchor representatives has even valuation at both. This gives an exact
linear functional separating \(\eta\) from the whole anchor span modulo
squares. In particular \(\delta(Q)\ne0\). Irreducibility of \(f\) excludes
rational 2-torsion, so all rational torsion has odd order and trivial Kummer
class. Therefore \(Q\) has infinite order and

\[
 \operatorname{rank}E_{-1}(\mathbf Q)\ge1.
\]

The locally admissible space now contains the direct sum
\(W_{-1}\oplus\langle\eta\rangle\), of dimension 19. This is a Selmer
subspace lower bound, not an exact Selmer dimension or a rank-19 claim.
No nonzero direction of \(W_{-1}\) has yet been realized.

## Reduction, search, and immediate certification

The runner is
[`run_fixed_field_point_realization.py`](../cas/run_fixed_field_point_realization.py).
Writing the three coefficients of \(\beta\gamma^2\) as \(q_0,q_1,q_2\), it
eliminates \(d\) to obtain the ternary conic \(q_2-uq_1=0\). PARI
`qfminimize`, `qfsolve`, and `qfparam` produce a parametrization
\(\gamma=M(s^2,st,t^2)^t\), retaining the rational transformation matrices.
The remaining equation is a binary quartic
\(d^2=-q_1(\gamma)\). Its integral model is minimized with
`hyperellminimalmodel` and reduced with `hyperellred`, as specified in the
[PARI reference](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html).
This minimizes an equivalent quartic presentation; it does not claim
minimality in the separate invariant theory of degree-four models.

Every conic identity and every quartic coordinate change is replayed over
\(\mathbf Q\). `hyperellratpoints` searches the reduced model with the
declared height argument. The rational points above the parameter infinity
are checked separately, and inverse changes use homogeneous coordinates so
that a pole in an affine coordinate never loses a point. Each hit must
satisfy the two original quadrics, the elliptic equation, and the actual
cubic-field identity

\[
 \frac{x-\alpha_{-1}}{\beta}=(\gamma/d)^2.
\]

The realized masks are immediately row-reduced in the certified anchor
Kummer basis. Their rank certifies independent Mordell--Weil directions;
the separating valuation adds the independent point \(Q\). Dependent
Kummer masks can conceal additional independent points, so this method
only certifies a lower bound.

There is also a second presentation for each selected class: replace
\(\beta\) by \(\beta\eta\). A point on that cover maps to \(P\), and
\(P-Q\) realizes the original class. The chord identity gives an explicit
square root of \((x(P-Q)-\alpha)/\beta\), and the runner verifies it rather
than simply trusting the class label.

Soluble classes form a subgroup of the Selmer group. What is not linear is
the success of a bounded point search in particular coordinate models.
Consequently a miss on every chosen basis representative does not justify
skipping their sums or translated presentations.

The `u=0` controls retain ten returned points, including parameter-infinity
points, across an anchor class and a two-anchor product. Their original
cover maps and Kummer identities replay exactly. Separate regressions check
the translation formula and reject altered point and coordinate-change
certificates.

## Completed bounded attack at `u=-1`

| Classes/presentations | Covers | PARI height argument | Returned points |
|---|---:|---:|---:|
| All basis, pair, and triple masks in \(W_{-1}\) | 987 | 100,000 | 0 |
| Basis and pair masks multiplied by \(\eta\) | 171 | 100,000 | 0 |
| Original basis, deeper pass | 18 | 1,000,000 | 0 |
| Translated basis, deeper pass | 18 | 1,000,000 | 0 |

All **1,194** searches completed, including minimization/reduction and the
parameter-infinity checks, with **zero timeouts and zero errors**. The
original masks comprise 18 basis vectors, 153 pairwise sums, and 816 triple
sums. There are 987 distinct nonzero classes of \(W_{-1}\) in this campaign;
the translated and deeper passes are additional searches for those same
realization questions.

The exact records are:

- [Original basis/pair/triple covers](../../artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_point_realization_v1.json).
- [Translated basis/pair covers](../../artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_translated_point_realization_v1.json).
- [Original basis at height one million](../../artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_basis_h1000000_v1.json).
- [Translated basis at height one million](../../artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_translated_basis_h1000000_v1.json).
- [Positive-control point maps](../../artifacts/generated-results/elliptic-curves/fixed_field_point_realization_positive_controls_v1.json).

Every stored model map and point certificate passes the exact replay, and
all seven point-realization regressions and nine original fixed-field
regressions pass. **The certified new rank lower bound is one, supplied by
\(Q\) outside \(W_{-1}\). The realized dimension inside \(W_{-1}\) is zero
in this experiment. All eighteen inherited directions retain classification
`UNKNOWN`; no nontrivial Sha class or rank upper bound is certified.**

## Reproduce the earlier bounded point searches

The searches use SageMath 10.9, PARI 2.17.3, PARI random seed 1, 100 decimal
digits for reduction, and a 30-second wall limit per cover. A completed row
is checkpointed with its exact inputs and output. Rerunning a command
replays and reuses its completed checkpoints; an interrupted or incomplete
row is retried. Timeouts and errors retain their last completed stage and
never count as point-search misses.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 sage -python elliptic-curves/cas/run_fixed_field_point_realization.py --max-weight 3 --height 100000 --timeout 30 --jobs 4
sage -python elliptic-curves/cas/run_fixed_field_point_realization.py --translate --max-weight 2 --height 100000 --timeout 30 --jobs 2 --output artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_translated_point_realization_v1.json
sage -python elliptic-curves/cas/run_fixed_field_point_realization.py --max-weight 1 --height 1000000 --timeout 30 --jobs 2 --output artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_basis_h1000000_v1.json
sage -python elliptic-curves/cas/run_fixed_field_point_realization.py --translate --max-weight 1 --height 1000000 --timeout 30 --jobs 2 --output artifacts/generated-results/elliptic-curves/fixed_field_u_minus1_translated_basis_h1000000_v1.json
```

Use a fresh `--checkpoint-dir` to repeat the actual searches. For cheap
algebraic verification use `--check --output PATH` on each artifact and
run the dedicated tests:

```sh
sage -python -m unittest elliptic-curves/tests/test_fixed_field_point_realization.py
```

The replayer checks the conic and quartic identities, elliptic minimal-model
transport, actual point/Kummer witnesses, the separating valuation for
\(Q\), and the declared class-mask coverage. It does not rerun the bounded
point enumeration and does not independently prove optimality of the
PARI-minimized models.

## Reproduction and claim boundary

```sh
sage -python elliptic-curves/cas/run_fixed_cubic_field_curve_family.sage --check
python3 -m unittest elliptic-curves/tests/test_fixed_cubic_field_curve_family.py
```

The canonical output is
[`fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json`](../../artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json).

The original local-pilot output certifies the family identities, fixed Kummer span, complete local
intersections, newly bad-prime coverage, and explicit cover inputs.  It does
not compute the full class group, the full 2-Selmer group of any \(E_u\), a
rational point on a new curve, a new rank lower bound, or an exact rank.
