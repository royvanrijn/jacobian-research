# Branch cancellation before constructing the two sections

The positive [correlated-cover objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains **OPEN**. There is a sharper algebraic gate, and it closes the same
Mestre identity on three further retained parents: alternate Q80, Curve302's
MW17 parent, and X1092 class1. On each, **every rational auxiliary function
still forces branch degree at least20 and genus at least9**. The fixed `u=2`
pair has height matrix `diag(24,24)` on a genus21 cover on each parent.

The new local criterion below says exactly when a simple coefficient divisor
*can* be removed by this identity. It is a necessary-and-sufficient local
existence criterion, not a sufficient condition for a globally low-genus cover.
The last section gives a different, explicit system with a common quartic
factor imposed from the start. That system has not been solved on these parents.

## 1. The exact global condition

For a nonzero rational function `R=N/H`, with coprime polynomials `N,H`, the
quadratic field is represented by `NH`, since `R=(NH)/H^2`. Write its actual
squareclass as

```
NH = d(t)*s(t)^2,          d squarefree, including its rational scalar.
```

If `n=degree(d)>0`, the branch degree on the original projective t-line is
`n+(n mod 2)`, including infinity. Thus the connected, nonconstant quadratic
cover has genus at most one **if and only if `1<=degree(d)<=4`**.
A constant nonsquare is a constant-field extension, outside this objective.
For two candidate discriminants, the equations to impose are

```
N1*H1 = d*s1^2,           N2*H2 = d*s2^2,
gcd(d,d')=1,              1<=degree(d)<=4.                 (1)
```

The same literal `d` must occur. Proportional branch polynomials with a
nonsquare scalar ratio do not give the same quadratic extension over Q(t).
Removing denominators, infinity, or multiplicity parity from (1) loses a
necessary condition. This is the usual double-cover calculation; see
[Rubin--Silverberg, Lemma2.5](https://www.maths.tcd.ie/EMIS/journals/EM/expmath/volumes/10/10.4/Rubin.pdf).

For the established Mestre identity, put `u=U/V` in lowest terms and

```
S=U^2+V^2,          K=U^4+U^2*V^2+V^4,
N_M=-A*B*S*(B^2*K^3+A^3*U^4*V^4*S^2).                    (2)
```

Then `D=N_M/V^14`. The denominator is already a square. The precise global
equation is consequently `N_M=d*s^2`, with the conditions in (1). One must
solve this equation, rather than merely make some resultant vanish.
Mestre's two-section identity is recorded in
[Rubin--Silverberg, Theorem3.7](https://www.maths.tcd.ie/EMIS/journals/EM/expmath/volumes/10/10.4/Rubin.pdf)
and the [previous moving-parent proof](R17_MESTRE_CORRELATED_SECTIONS_AND_GENUS_GATE_2026-09-13.md).

## 2. Sharp local cancellation criterion for Mestre

Let `A,B in Q[t]` be nonzero, squarefree and coprime. Let `pi` be a monic
irreducible factor of one coefficient, `k_pi=Q[t]/(pi)`, and use `ord_pi`
normalized by `ord_pi(pi)=1`. Exclude `D=0`, which defines no quadratic field.
There exists `u in Q(t)^*` making `ord_pi(D)` even exactly under these conditions:

| Divisor | Necessary and sufficient condition for some auxiliary function |
|---|---|
| `pi | A` | `k_pi` contains `i` or a primitive cube root of unity |
| `pi | B` | `k_pi` contains `i` |

When these conditions hold for several coefficient divisors, one polynomial
`u` can remove all those divisors simultaneously. No assertion is made about
the new branch points created away from them.

### Necessity and exact valuation cases

Write `v=u^2`, `s=v+1`, `k=v^2+v+1`, and

```
G=B^2*k^3+A^3*v^2*s^2,             D=-A*B*s*G.
```

If `n=ord_pi(u)` is nonzero, the two summands of `G` have unequal valuations:

| Divisor | `n<0` | `n>0` |
|---|---:|---:|
| `pi | A` | `ord_pi(D)=1+14*n` | `ord_pi(D)=1` |
| `pi | B` | `ord_pi(D)=3+14*n` | `ord_pi(D)=3` |

All entries are odd. Zeros or poles of `u` cannot remove these divisors.
For `n=0`, put `j=ord_pi(s)` and `l=ord_pi(k)`. They are nonnegative and
cannot both be positive, since `s=0` implies `k=1` in the residue field.

At `pi | A`:

* If `j=l=0`, then `ord_pi(D)=1`.
* If `j>0`, then `G` is a unit and `ord_pi(D)=1+j`.
* If `l>=2`, the second summand of `G` has valuation3, strictly below `3*l`,
  and `ord_pi(D)=4`.
* If `l=1`, both summands have valuation3. Writing `e=ord_pi(G)-3>=0`,
  one has `ord_pi(D)=4+e`. Cancellation requires even `e`.

Thus removal requires a residue with `u^2=-1` or `u^4+u^2+1=0`. The latter
condition is equivalent to the residue field containing a primitive cube
root of unity: one direction follows from `u^2`, and in the other direction
one can take `u` itself to be a primitive cube root.

At `pi | B`, if `j=0` the second summand of `G` is a unit, so
`ord_pi(D)=1`, regardless of `l`. For `j>=2`, the first summand has valuation2,
strictly below `2*j`, giving `ord_pi(D)=3+j`. For `j=1`, write
`e=ord_pi(G)-2>=0`; then `ord_pi(D)=4+e`. Therefore removal requires
`u^2=-1` in `k_pi`.

This also explains why a permitted cyclotomic residue is insufficient for
a *specified* `u`: its order of contact and, in the tied cases, the leading
coefficient cancellation still matter.

### Sufficiency and simultaneous construction

Each of the residue roots just used is simple in characteristic zero.
Hensel lifting in `Q[t]/(pi^r)` lets us prescribe an exact finite contact order:

* at `pi | A` with `i in k_pi`, choose `ord_pi(u^2+1)=1`, giving valuation2;
* at `pi | A` with a cube root, choose `ord_pi(u^4+u^2+1)=2`, giving valuation4;
* at `pi | B` with `i in k_pi`, choose `ord_pi(u^2+1)=3`, giving valuation6.

Choose the next coefficient to avoid a higher order. These are finite
polynomial congruences. The Chinese remainder theorem combines them at
distinct factors of `AB`; further polynomial perturbations preserving the
congruences avoid identically zero `u` or `D`. This proves sufficiency.

Consequently define the *forced coefficient divisor* to consist of the
factors of `A` whose residue fields contain neither cyclotomic element,
together with the factors of `B` whose residue fields do not contain `i`.
Its degree is a uniform lower bound on the branch degree for every `u`.
**Degree greater than four closes this identity for the present objective.**
Degree at most four merely passes this necessary gate; equation (2) still
controls all remaining branch points.

## 3. Application to the retained parents

The frozen input projects only generic equation coefficients from the retained
sources. Long Weierstrass models are put in short form using
`A=-c4/48`, `B=-c6/864`; this is a coordinate change over the same parameter
field. No exceptional point, specialization value, or specialized rank is used.

For each parent, the reductions in the following table prove that `A` is
irreducible of degree8 and `B` irreducible of degree12. Each displayed local
root is simple and the leading coefficient is a unit.

| Parent | Irreducibility prime for A | Local A place `(p,root)` | Irreducibility prime for B | Local B place `(p,root)` | Smoothness prime |
|---|---:|---:|---:|---:|---:|
| Published R17 control |29|`(59,51)`|59|`(107,83)`|131|
| Alternate Q80, direct11952 |107|`(71,36)`|23|`(127,77)`|131|
| Curve302 MW17 parent |47|`(239,167)`|83|`(79,43)`|149|
| X1092 class1 |71|`(131,100)`|19|`(127,123)`|151|

Irreducibility makes each coefficient define a single root field. A simple
root gives an embedding of that field into `Q_p`. All A-place primes satisfy
`p=11 mod12`, so `Q_p` contains neither `i` nor a primitive cube root. All
B-place primes satisfy `p=3 mod4`, so `Q_p` contains no `i`. These exclusions
follow already by reducing prime-to-p roots of unity to `F_p^*`.
The entire degree20 coefficient divisor is therefore forced on all four
parents, proving the genus lower bound9 for arbitrary rational `u`.

The last column independently verifies that

```
D_2=-5*A*B*(9261*B^2+400*A^3)
```

is squarefree of degree44, that `4*A^3+27*B^2` is squarefree of degree24,
and that they are coprime. There is no branch at infinity and the original
infinity fibre is smooth. For `w^2=D_2` the exact sections are

```
P=(-21*B/(5*A),  w/(25*A^2)),
Q=(-21*B/(20*A), w/(200*A^2)),
x(P+Q)=-(7/15)*B/A-(20/81)*A^2/B.                         (3)
```

The pullback has `chi=4` and only irreducible fibres. Both points meet zero
once over each of the eight A roots; their sum meets zero once over each
of the twenty A and B roots. Exact noncancellation and the degree bounds
exclude other zero intersections, including infinity. The Shioda heights are
therefore24,24,48, so the cross-pairing is zero. See the
[previous detailed pole proof](R17_MESTRE_CORRELATED_SECTIONS_AND_GENUS_GATE_2026-09-13.md#exact-independence-modulo-the-inherited-subgroup)
and [Schuett--Shioda, section11](https://arxiv.org/abs/0907.0298)
for the height convention. The arithmetic replay verifies these hypotheses
on every new equation, including both section equations and the addition law.

Both points are anti-invariant and hence orthogonal to the inherited group.
The established MW17 parent theorems give rank at least19 over these genus21
covers. No exact total rank or improved rational specialization is asserted.
The result covers these specific fibrations; other X1092 classes and different
identities are outside it. The old R17 certificate is retained unchanged.

## 4. A different system with the common quartic built in

There is a particularly concrete sufficient construction on a rootless
`24I1` K3 parent with smooth infinity. Seek polynomials

```
degree(d)=4,     gcd(d,d')=gcd(d,4*A^3+27*B^2)=1,
degree(x_i)<=4,  degree(r_i)<=4,   r_i!=0,   x_1!=x_2,
x_i^3+A*x_i+B=d*r_i^2       (i=1,2).                      (4)
```

This is a pair of coefficient identities with the same quartic already
present. Equivalently impose the first identity and the factor identity

```
(x_1-x_2)*(x_1^2+x_1*x_2+x_2^2+A)
    =d*(r_1-r_2)*(r_1+r_2).                              (5)
```

Equations (4) or (5) are an exact target for symbolic elimination. The
nonzero and coprimality conditions must survive elimination. This does not
claim that these overdetermined equations have a solution for a chosen parent.
Higher-pole sections and nonzero trace parities are outside this polynomial
chart; it is a sufficient family, not a complete parametrization of all gains.

If a solution exists, the points `P_i=(x_i,r_i*w)` on `w^2=d` are
anti-invariant, integral at finite places and have no zero intersection at
infinity by the displayed degree bounds. The pullback has `chi=4` and no
reducible fibres, so each has height8 and the torsion group is trivial.
If they were rationally dependent, equality of heights would give
`P_1=+/-P_2` modulo torsion. That would force `x_1=x_2`. Thus **(4) alone
proves two independent new directions modulo the inherited group**.

These are singular higher-arithmetic-genus bisections of precisely the kind
left open by the smooth-Q80 exclusion. Their images have class `2O+4F`,
self-intersection8 and arithmetic genus5. The connected quadratic
normalization has genus1, so the total singularity delta invariant is4.
No assertion that all four singularities are ordinary nodes is required.
This avoids requiring a second smooth genus-one bisection image.

The rational-base obligation remains separate: exhibit a rational point
on `w^2=d` and a certified nontorsion point on its pointed genus-one model.
Together with (4) and the nonconstant parent, the usual specialization theorem
would then supply infinitely many distinct rational parameters with rank at
least19. **No solution of (4), nor such a positive-rank quartic, is produced
in this note.**

For arithmetic-genus-two carriers, reuse the existing
[singular-member halving gate](../elliptic-curves/notes/DET1092_RR_NET_SINGULAR_MEMBER_GATE_2026-09-08.md)
and [uniform nonzero-parity genus-nine theorem](../elliptic-curves/notes/DET1092_TRACE_PARITY_DESCENT_2026-09-09.md).
Those results already rule out a rational or elliptic parametrization of the
stated dense-open node locus. They do not exclude isolated rational nodes,
exceptional charts, or the present genus5 image construction. They are prior
work, not new theorems or new rational-point computations here.

## 5. Evidence and replay boundary

The [input](../artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json)
fixes all four sources, their hashes, primes `5<=p<2000`, `u=2`, a40 CPU-second
limit and4GiB address space before the gate. The
[result](../artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/result.json)
and four parent checkpoints retain every witness. Missing witnesses are
`UNKNOWN`; an interrupted run has no final certificate.
The [independent replay](../artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/independent-replay.json)
uses Python fractions and integers, with no Sage or producer import. It
reconstructs short models by b-invariant formulas, verifies (3) coefficientwise,
and uses Rabin irreducibility, polynomial gcds and simple-root tests to replay
the finite gates. Parent ranks and saturated bases are inherited dependencies.

The nine small tests include altered source projections and local witnesses,
incomplete coverage, degree loss, the exact zero/pole valuation table, and
cyclotomic examples showing both cancellation and its failure at an incorrect
contact order. These examples check the sharp local criterion's boundary;
the general local proof, global double-cover formula, and conditional
construction (4) remain written mathematics, not formal verification.

```sh
python3 research/elkies-k3/scripts/verify_mestre_parent_branch_cancellation.py
python3 -m unittest discover -s research/tests -p 'test_mestre_parent_branch_cancellation.py' -q
```

Discovery is reproducible in a fresh output directory with
`gate_mestre_parent_branch_cancellation.sage freeze --output DIR`, then
`run --output DIR`, both under `sage -python` with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. Outputs refuse overwrites.
The [execution receipt](../artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/execution.json)
records the commands and component timings. No bisection enumeration,
auxiliary-function search, or rational-point campaign was performed.
