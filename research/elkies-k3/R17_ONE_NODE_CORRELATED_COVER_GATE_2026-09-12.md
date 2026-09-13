# One-node carriers and the common quadratic-cover gate

The requested [two-gain construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains **OPEN**. All **1,675 pairs** of the frozen 67 height-six regular chord
nets and 25 height-eight smooth genus-one pencils are excluded, for every
rational parameter and node position. This extends the earlier smooth-pencil
exclusion to a different carrier family. It is not a global obstruction to
two new sections on a quadratic cover.

## Construction and independence criterion

Work on a rootless elliptic K3 in short form
`E: y^2=x^3+A(t)x+B(t)`. For a height-six trace write

```
T=(X,Y)=(Nx/h^2,Ny/h^3),   deg h=1,
M0*Nx+Ny=0 mod h^2,       deg M0<2,
M=M0+(a+b*t)*h^2,
q=(M^4-6*M^2*Nx-8*M*Ny-3*Nx^2-4*A*h^4)/h^6.
```

The [regular chord calculation](RANK_MUTATION_AND_LIFT_THEOREMS.md#proposition-f11-the-height-eight-genus-one-bisection-pencil)
gives a polynomial sextic. Its arithmetic-genus-two double cover has genus-one
normalization when, in binary-form notation,

```
q=ell_r^2*d,       d a squarefree binary quartic.
```

Here `r` is rational: the removed degree-two divisor is twice one point,
and is defined over Q. This includes a double factor at infinity and the
triple-root case where the remaining quartic also vanishes at `r`.
At finite `r`, the lift on `s^2=d(t)` is

```
x0=(M^2-Nx)/(2*h^2),
y0=(M/h)*(x0-Nx/h^2)-Ny/h^3,
x=x0+h*(t-r)*s/2,
y=y0+M*(t-r)*s/2.
```

At a smooth finite fibre with `h(r)!=0`, the repeated-root condition is
constructive. The residual chord is tangent at a rational half `Q` of `T(r)`.
Its slope is a rational root of

```
F(m)=m^4-6*X(r)*m^2-8*Y(r)*m-3*X(r)^2-4*A(r).
```

Conversely each such root gives

```
x_Q=(m^2-X(r))/2,
y_Q=m*(x_Q-X(r))-Y(r),       2Q=T(r).
```

To see this, divide the Weierstrass cubic minus the squared chord by the
known root `x-X`. The residual quadratic has discriminant `F(m)`. Moreover
`F_m=8*y_Q`, which is nonzero: a point of order two could not double to the
finite trace `T(r)`. Differentiation therefore determines the unique first
jet of the tangent slope:

```
m'=(6*X'*m^2+8*Y'*m+6*X*X'+4*A')/(4*m^3-12*X*m-8*Y),
m_base=M0/h,
L=(m-m_base(r))/h(r),
b=(m'-m_base'(r)-L*h'(r))/h(r),
a=L-b*r.
```

All derivatives in these expressions are evaluated at `r`. This constructs
`q(r)=q'(r)=0` without selecting an exceptional specialized point. One must
still check squarefreeness of `d`, its constant squareclass, and its rational
points. A rational singular point need not give a rational point on the
normalization.

There is a useful conditional independence proof. Suppose a height-eight
smooth carrier and this height-six carrier define the same quadratic
extension. Require that the four branch points lie over smooth original
fibres and that both lifted sections are disjoint from zero; the regular
polynomial coordinate bounds provide the latter check. The pullback has
`chi=4`, no reducible fibres, and both lifts have height 8. If
`V=2P-T=P-sigma(P)`, trace and deck invariance give

```
height(V)=4*height(P)-2*height_parent(T).
```

Thus the two anti-invariant heights are **16 and 20**. They are orthogonal
to the inherited subgroup. A rational dependence between them would force
their height ratio `5/4` to be a rational square. It is not, so they would
supply two independent directions modulo the inherited subgroup. This is
a sufficient criterion, not an actual common-cover construction.

## The complete projective comparison

The source words are all 67 parity-distinct height-six words supported on
two published basis vectors, with the same deterministic sign preference
as the [frozen smooth bank](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).
Only generic equation, basis and Gram fields enter the packets.

Expand the source sextic in `L=a+b*t`:

```
q=q0+L*q1+L^2*q2+L^3*q3+L^4*q4,
q1=(4*M0^3-12*M0*Nx-8*Ny)/h^4,
q2=6*(M0^2-Nx)/h^2,     q3=4*M0,     q4=h^2.
```

Homogenization gives a degree-four parameter map from `P2` to binary sextics:

```
Q(a:b:c)=sum_j c^(4-j)*(a+b*t)^j*qj(t).
```

The target map uses the already verified invertible branch matrix `B`:

```
H(lambda,r)=ell_r^2*B*Veronese_4(lambda),
(lambda,r) in P1 x P1.
```

Equality of the quadratic extensions implies `[Q]=[H]` in `P6`. This forgets
the constant twist, so it is only a necessary test for accepting a match.
It is sufficient for an exclusion.

At a prime with integral coefficients and invertible target matrix, enumerate
**all** `p^2+p+1` source parameters and all `(p+1)^2` target parameters. A source
family is usable only when none of its rational reduction points is a base
point of `Q`. Every rational parameter tuple then has primitive `Z_p`
coordinates, and both evaluated coefficient vectors remain primitive.
A proportionality over Q reduces to a proportionality over `F_p`. Disjoint
projective images therefore exclude the whole rational pair.

This argument includes `c=0`, infinite `lambda`, and infinite `r`. At `c=0`
the source is the degenerate square `h^2*(a+b*t)^4`; at infinite `lambda` the
target quartic is also a square. These points must be included because a
smooth rational characteristic-zero candidate may reduce to one of them.
An affine miss or singular reduction alone is not an exclusion.

| Prime | New excluded pairs | Pairs remaining |
|---:|---:|---:|
| 101 | 948 | 727 |
| 103 | 598 | 129 |
| 107 | 110 | 19 |
| 109 | 12 | 7 |
| 113 | 6 | 1 |
| 127 | 1 | 0 |

The initial frozen prime list ended at 113. Its one surviving pair was
source `P6-P12` and target `P11+P12`. All its retained local matches had a
singular quartic, which did not decide the rational question. A separately
frozen continuation allowed exactly this pair at primes 127, 131 and 137.
At 127 its 16,257 source parameters contain no base point and no target-image
match. The continuation stopped there. The original partial certificate and
its diagnostic are preserved.

The independent Python checker reconstructs the exact source coefficients,
replays the generic trace and target-pencil identities, and recomputes every
projective image. It evaluates a binomial coefficient tensor by Horner's rule;
the Sage producer uses direct polynomial convolution. The final replay checks
3,820,642 source parameters and exact coverage of all 1,675 exclusions.

Artifacts: [input](../artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1/input.json),
[initial result](../artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1/result.json),
[continuation input](../artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1/continuation-input.json),
[last exclusion](../artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1/continuation-result.json),
[complete independent replay](../artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1/independent-replay-complete.json).

## Retained construction attempts

Before the projective comparison, an exact construction panel tested the
67 traces at all rational projective parameters of height at most 32,
their 67 trace poles, and infinity. After deduplication this is 1,361 fibres
and **91,187 cases**:

| Outcome | Cases |
|---|---:|
| No rational half, certified by good finite reduction | 87,911 |
| No rational root of the exact tangency quartic | 3,205 |
| Pole or infinity equations subsequently excluded | 70 |
| Rational half producing a split carrier | 1 |

At a trace pole, `q(r)=0` is a cubic in `a+b*r`. At infinity, the leading
sextic coefficient is a quartic in `b`. All 70 retained boundary equations
have an exact no-root witness. For every negative rational-root decision,
the portable checker clears denominators and content and checks a prime
with no **projective** root. It never trusts a factorization status flag.
Its finite tables are reconstructed through the tangency quartic, independently
of the producer's enumeration of the image of doubling on finite curves.

The sole half occurs for `T=P3-P16` at `r=-21/4`. The normalized polynomial is

```
d=442567226524897748498605056
  *(t^2+(430623/94393)*t+348480/94393)^2.
```

The scalar is a rational square. Exact generic group-law identities identify
the two branches as `P1+P3` and `-P1-P16`. Thus this is a split inherited
control and adds no direction.

An additional fixed experiment sought node positions from `T=2Q`, where
`Q` ranges over the 34 signed published basis sections. All 2,278 generic
word pairs have no finite rational intersection. An independent proof uses
1,139 numerator polynomials for `x(2Q)-x(T)` and projective no-root witnesses,
so this result has no parameter-height bound. It does not cover other
inherited words or noninherited rational halves.

Retained evidence: [finite panel and branch identities](../artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1/independent-replay-with-branch-identities.json),
[compressed panel checkpoints](../artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1/trace-records.zip),
[generic-word intersection replay](../artifacts/generated-results/elkies-k3-r17-inherited-halving-nodes-v1/independent-replay.json),
[compressed generic-word checkpoints](../artifacts/generated-results/elkies-k3-r17-inherited-halving-nodes-v1/trace-records.zip).

## Replay and remaining boundary

```
python3 research/elkies-k3/scripts/verify_r17_one_node_complete_pairs.py
python3 research/elkies-k3/scripts/verify_r17_inherited_halving_nodes.py
python3 -m unittest discover -s research/tests -p 'test_r17_one_node_*.py' -q
```

The first command also verifies the fixed fibre panel and the prior 25-pencil
certificate. The eleven small failure controls cover projective infinity,
vanishing source vectors, nonintegral reductions, denominator content,
composite witnesses, nonexact division, constant scaling, altered split-section
data, and a missing final-pair certificate. The two construction workers and main comparison
each had a 120 CPU-second, 4 GiB cap; the single-pair continuation had a
20-second cap. All use one BLAS thread. Execution receipts preserve the
actual component timings and the initial input-serialization failure.

To package a fresh constructor output for the portable replay, run
`python3 research/elkies-k3/scripts/pack_r17_trace_checkpoints.py --input-dir OUTPUT`
and then the corresponding verifier with `--input-dir OUTPUT --export`.
The finite-panel boundary worker takes `--input-dir OUTPUT` and runs before
that panel's portable replay. Use a fresh output directory for each run;
every producer, archive writer and certificate exporter refuses overwrites.

The complete comparison excludes these 67 nets paired with these 25 pencils.
It does not exclude all 63,925 norm-eight translation classes, pairs of two
other singular carriers, higher-pole sections, arbitrary twists, alternate
parents or other arithmetic sources. No bound on a specialized rank follows.

A next construction must supply a different carrier family or an explicit
two-section twist identity and establish its branch count before a parameter
campaign. General two-section twist constructions are available, but do not
by themselves meet the low-genus covering-base condition; see
[Rubin–Silverberg, Theorem 3.7 and Remark 2.12](https://www.maths.tcd.ie/EMIS/journals/EM/expmath/volumes/10/10.4/Rubin.pdf).
For their constant-parent setting, the twist rank is bounded by the genus
of its hyperelliptic cover. Applying an auxiliary-parameter identity to a
moving high-rank parent needs a separate branch and independence proof.
No further campaign is scheduled by this note.

There is still no constructed shared quadratic cover with two new independent
sections and infinitely many rational base points.
