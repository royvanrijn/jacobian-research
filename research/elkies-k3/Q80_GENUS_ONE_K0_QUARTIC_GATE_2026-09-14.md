# Quartic norms confine Q80 genus-one k0 to nodal denominators

On the literal direct11952 alternate-Q80 MW17 parent, a polynomial genus-one
common-quartic pair with different nonzero2-torsion choices at all four branch
points cannot have both abscissas integral at131. Every such pair must have

```
D_0 mod131 proportional to (t^2+62*t+88)^2,
at least one abscissa of negative Gauss coefficient valuation.           (1)
```

This is a necessary boundary, not an exclusion of the five `k=0,j=0,...,4`
allocations. Both the mixed integral/negative case and the two-negative case
remain **UNKNOWN**. The branch field can be ramified over the unramified
quadratic extension of Q_131. The earlier k2 orientation contradiction
cannot be applied without its agreement and local-field hypotheses.
No infinite rank-at-least19 source is constructed.

## 1. Scope and inherited facts

Write the full fixed-parent chart as

```
x_i^3+A*x_i+B = eta*D_0*r_i^2,
deg x_i<=4, deg r_i<=4,
D_0 squarefree of binary degree4, gcd(D_0,Delta)=1,
gcd(D_0,x_1-x_2)=1.                                                     (2)
```

The characteristic-zero branches lie over smooth parent fibres. Scalar
classes, coefficient denominators, infinity and higher contacts are retained;
the covering genus-one curve need not have a rational point. Normalize D_0
to a primitive integral binary form and eta to valuation epsilon in `{0,1}`.

The [k1 denominator proof](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md#2-coefficient-integrality-for-k1)
already proves that any negative abscissa valuation in this chart requires
the reduction in (1). It uses the actual stable-map image cycle, not just
parity of the reduced coefficient identity. We reuse that written theorem.
It remains to exclude pairs with both abscissas integral.

Then both r_i are integral. If epsilon were one, the reduction of (2)
would give a nonzero rational2-torsion section of the reduced parent.
The reduced surface has24I1 fibres and every nonzero section has height
`4+2*(T.O)>0`, excluding such torsion. The same argument prevents r_i
from reducing identically to zero. Thus eta is a unit.

The complete [degree-two norm table](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md)
and [norm-four section roster](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
are retained inputs, with their original completeness and replay receipts.
All132 projective rational fibres at131 are smooth. The reduced discriminant
has simple irreducible factors of degrees2,3,19. Its quadratic factor is

```
q_0=t^2+62*t+88.
```

There are110 pairwise distinct nonzero rational root codes. Sixteen rational
fibres split fully. The complete8,515 nonrational quadratic base-orbit roster
has exactly one orbit with a zero root norm: at the nodal q_0 fibre the
three codes, preserving double-root multiplicity, are `(13412,13412,0)`.
The simple root is `S=73+6*t mod q_0`, and the node is `N=29+128*t mod q_0`.

## 2. Both selected-root norms must vanish

For each of the17 inherited sections T, use the regularized character
`X_T-e`: at an ordinate zero its value is `3*e^2+A`, and at an abscissa
pole its squareclass is1. The earlier reciprocity proof gives zero total
norm code on the degree-four selected-root divisor of each section.
As the root-character product is trivial, the unused-root divisor has
zero total code as well. Thus the three componentwise XOR sums are

```
C_D(e_1)=C_D(e_2)=C_D(e_3)=0.                                           (3)
```

This is stronger than one unused-root zero. At a smooth residual fibre,
integral polynomial abscissas choosing different actual roots must reduce
to different roots, by uniqueness of the simple-root Hensel lift. Each
choice is constant throughout one branch residue cluster. Norms of the
unramified root-character classes acquire the total local branch degree
as exponent, also for ramified branch fields. At nodal residues use the
good cubic-root-curve model and retain root multiplicities, as in the
degree-two proof. An inherited section never meets the node, so its
character there is a unit; a simple-root zero has the stated regularization.

If D_0,bar is neither squarefree nor a scalar multiple of a square, its
odd-multiplicity divisor has binary degree2. Even multiplicities disappear from (3). That odd
divisor is either two distinct rational residues or one quadratic residue.
In the rational case, injectivity of the110 nonzero codes precludes even
one selected-root norm cancellation between distinct sites. In the
quadratic case, the retained roster has no fibre with all three root codes
zero, including the nodal `(13412,13412,0)` case. Hence these reductions
are excluded. Every remaining reduction is either squarefree or a scalar
multiple of a square H^2 of a binary quadratic.

## 3. All squarefree degree-four norm patterns are empty

The five factor-degree partitions are

```
1+1+1+1,   1+1+2,   2+2,   1+3,   4.                                  (4)
```

For the first three, use the complete retained rational/quadratic tables.
At each full-splitting site, permute the three root labels. The signature
of a labelled pair of rational sites is the sorted triple of componentwise
XORs. All4,320 ordered rational-pair choices are included. Two disjoint
rational pairs never have the same signature; no rational-pair signature
equals a quadratic signature; and distinct quadratic orbits never have
the same signature. These assertions exclude the first three partitions
in (4). The nodal quadratic is included with its two node-root slots.
The independent checker quotients by the simultaneous S3 relabelling
before forming its joins.

For `1+3`, the cubic-orbit signature must equal the signature of one of the
sixteen fully split rational fibres. The new complete cubic run checks
all749,320 genuine cubic base orbits in `F_131[theta]/(theta^3+theta+3)`.
Exactly125,595 smooth cubics split fully and one cubic orbit is nodal.
None has a permitted rational signature. Independent generic-polynomial
replay uses repeated p-Frobenius substitution, different factor splitters
and determinant norms. It reproduces the complete counts and empty match
set. The earlier k1 run alone did not retain this signature condition;
the new run is separately frozen.

For `4`, use the exact field

```
F_131[i,j],  i^2=-1,  j^2=1+i,
encode((a+b*i)+(c+d*i)*j)=a+131*b+131^2*c+131^3*d.                       (5)
```

Here -1 is nonsquare in F_131. The norm of1+i to F_131 is2, also nonsquare,
so the second extension is quadratic. Equivalently j has irreducible
polynomial `j^4-2*j^2+2`. The genuine quartic elements are exactly those
with `(c,d)!=(0,0)`. Their p-Frobenius sends the j coefficient by

```
L(c,d)=(117*(c+d),117*(c-d)),      L^2=-1.
```

Choose the least encoded pair in each four-cycle under L. There are4,290
such pairs, and all17,161 choices of `(a,b)` for each pair give exactly one
representative of every genuine quartic base orbit. Consequently the
complete domain has

```
4290*17161 = (131^4-131^2)/4 = 73,620,690 orbits.                        (6)
```

The subfield elements are excluded: inherited characters become squares
there after an additional even-degree field extension, which would create
irrelevant positive controls. There is no nodal quartic orbit because the
bad-fibre degrees are2,3,19.

For a cubic `x^3+a*x+b`, put `delta=(b/2)^2+(a/3)^3`. Since this field
contains the third roots of unity, full splitting requires delta square.
For a square root s, choose a nonzero resolvent `U=-b/2+s` or `-b/2-s`.
Writing `q=131^4`, the order `q-1=3*98166640` has exactly one factor3.
The exponent `d=65444427` satisfies `3*d=2*(q-1)/3+1`. Thus `u=U^d` is
a cube root exactly when `u^3=U`; this test is necessary and sufficient
for U to be a cube. Put `v=-a/(3*u)` and use

```
u+v,    zeta*u+zeta^2*v,    zeta^2*u+zeta*v,
zeta=(65+19*i).
```

These are the three cubic roots when the test succeeds. Conversely the
Lagrange resolvents of three field-valued cubic roots give the needed cube
roots of `-b/2 +/- s`, proving completeness of this splitting test. Every
returned root is checked against the cubic. The producer uses tower square
roots and tower norms. The independent program uses monomial arithmetic
modulo `j^4-2*j^2+2`, Tonelli--Shanks square roots, the opposite resolvent,
Frobenius-product norms, and quotient evaluation with field inverses.

Both programs enumerate all73,620,690 orbits. Exactly12,267,623 smooth
fibres split fully, and **none** satisfies (3). Each fibre is rejected as
soon as one of its17 inherited characters fails; a full all-zero candidate
would have to pass all17. Counts by first failing character, root candidates
and checksums agree in every one of ten exact, disjoint checkpoint ranges.
The counts by first failing character, followed by the all-pass count, are

```
9201540,2299578,574929,143902,35609,9116,2196,563,141,
34,7,3,4,0,1,0,0,0.
```

Thus all squarefree reductions are excluded, without a characteristic-zero
coefficient-height or denominator cutoff.

## 4. Square reductions have no integral pair

Suppose `D_0,bar=c*H^2`. Absorb c into the unit scalar. The two reductions
give sections `(x_i,bar,H*r_i,bar)` on the same constant twist of the
reduced24I1 parent. Both miss the zero section, so both have height4.
The constant nonsquare twist has arithmetic MW rank at most one by the
retained Frobenius bound, and no torsion. Two height-four points in this
rank-one lattice have the same abscissa. At every smooth branch residue
Hensel uniqueness contradicts k0. If the support is only q_0, an ordinary
section cannot meet the node: differentiating its equation at the node
would give `0=f_t`, although f_t is a unit there. Both sections would
therefore choose the same simple-root lift. This also contradicts k0.

In the square scalar case, Hensel's lemma absorbs the actual unit scalar
over Q_131. The reductions belong to the full retained arithmetic norm-four
lattice: exactly2,626 signed points, or1,313 distinct abscissas. They have
polynomial coordinates of degrees at most4 and6. H must divide both
homogeneous ordinate forms. At each root of H their abscissas must be
different. A nodal root cannot occur in H for a pair, by the preceding
ordinary-section argument.

The complete unordered pair census examines861,328 pairs. Exactly69 have
a common homogeneous ordinate divisor of degree at least2. Removing every
support where the two abscissas agree leaves exactly30 quadratic-contact
pairs. Their H are all affine monic quadratics, including split quadratics;
no infinity or repeated-linear case survives. Only point732 in the full
ordinary norm-four roster has q_0 dividing its ordinate, and its abscissa
reduces to the simple root S.

The producer computes common ordinate gcds. The independent checker first
factors each ordinate in degrees1 and2, with multiplicities, using prime-field
polynomial arithmetic. It forms all binary quadratic divisors, including
infinity, and joins point incidences by H. The entire30-pair list agrees.
Ordinate signs do not change either the contact condition or lift consistency.

For each pair, normalize the actual quartic D to have leading coefficient1.
Its residue is H^2. Use24 coefficient variations: five for each x_i, five
for each r_i, and four for D after fixing its leading coefficient. In each
of the two coefficient identities the linearization is

```
(3*x_i^2+A)*delta x_i - 2*D*r_i*delta r_i - r_i^2*delta D.
```

There are26 coefficient equations. Evaluate the literal rational parent
modulo131^2 and retain the divided error vector. Every one of the30 matrices
has rank24; augmenting by that error vector gives rank25. No solution lifts
even to131^2. The independent checker rebuilds the equations with another
column order and exact integer Gaussian elimination, including the error
vector. A full-rank matrix alone would not give this exclusion.

This eliminates all pairs of integral abscissas. Combining with the inherited
denominator theorem proves (1). No claim is made that both abscissas must be
negative, that the remaining branch fields split as two quadratics over
Q_131, or that the remaining k0 allocations are empty.

## 5. Evidence and replay boundary

The [quartic input](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1/quartic-input.json)
was frozen after a one-pair benchmark, before either complete run. Its
[producer/replay comparison](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1/quartic-comparison.json)
binds all ten ranges from each implementation. Every range has a90 CPU
second and512MiB cap, with a95-second wall timeout. Full producer CPU was
45.18691 seconds and full independent replay CPU69.05172 seconds.

The [integral-gate input](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1/input.json)
then freezes these results with the composite and contact programs.
The [independent receipt](../artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1/independent-replay.json)
rechecks factor incidence, all first lifts, the complete cubic signature run
and the completed quartic replay's source-bound coverage. By default it
reuses the completed large quartic calculation; `--replay-quartic` requests
a fresh complete arithmetic replay. Drivers are capped at30 CPU seconds
and2GiB; each cubic worker has the same90-second/512MiB caps as a quartic
range. Eight targeted tests pass, including fully split positive controls,
irreducible-cubic negatives, zero regularization, multiplicities, infinity,
a missing contact and a changed coefficient modulo131^2.

An initial square-lift preview incorrectly asserted that all contact H were
irreducible. The assertion stopped that preview. Its exact source and failure
receipt are retained; the corrected preflight and frozen calculation include
split H. No failed frozen certificate was rewritten as a pass.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_genus_one_k0_gate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m unittest discover -s research/tests -p 'test_q80_genus_one_k0_gate.py' -q
```

The finite computations have independent arithmetic replay. The reciprocity,
norm specialization, constant-twist height and stable-map deductions are
written mathematics, not formal verification or external review. The older
genus-zero norm-eight singular-polynomial replay gap is neither used nor
upgraded. All five k0 allocations and the required actual-MW17 infinite
quadratic-base construction remain open.
