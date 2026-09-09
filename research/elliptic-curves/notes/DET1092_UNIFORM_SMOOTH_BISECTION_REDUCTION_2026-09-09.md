# All norm-ten bisections have smooth reduction at149 and151

The later [Euclidean conic formula](DET1092_EUCLIDEAN_BISECTION_FORMULA_2026-09-09.md)
constructs the RR line without a matrix solve and resolves the three rational
UNKNOWN masks below: all are nonsplit at302 and the eight controls. The
thirty-trial modular record here remains an unchanged historical checkpoint.

## Theorem and scope

**New deduction from completed geometry; independently verified arithmetic
application.** Every one of the40,917 minimum norm-ten smooth rational
bisections has smooth geometrically irreducible reduction at both149 and151.
This is a uniform theorem, not an equation-level enumeration of those curves.
The extra geometric Mordell--Weil direction in either reduction cannot have
the height-six vector that a broken bisection would require.

At either prime, four fixed base charts always suffice for the maximal-pole
RR construction. Its interpolation matrix has rank19, for every one of the
40,917 classes. A vertical-line formula handles a trace specializing to the
zero point. Thus the earlier coordinate-failure cases can be removed by
generic geometry rather than by increasing a search budget.

**Verified bounded application.** Both full17-section height Grams and the
same fifteen old orbits were checked at149 and151. All30 trials pass. The
reciprocal chart repairs orbit61 at149; it is locally split there, still
UNKNOWN over `Q`. There are **no additional rational non-incidence or seed
claims** from this replay. Prime157 is retained as a failed global-geometry
regression, not retroactively removed from the old experiment.

No full-atlas run, new prime, rational parameter, exceptional point, point
search, class/unit group, production change or detached job was used. The
two old integral norm equations were audited but not recomputed or solved.

## 1. Verified reduction inputs

The [completed parent geometry](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md)
and [arithmetic Picard calculation](DET1092_SURFACE_BRAUER_TRIVIALITY_2026-09-09.md)
provide smooth proper K3 reductions with24 irreducible `I1` fibres at149
and151, and characteristic polynomials

\[
 R_p(T)=(T-p)^{19}(T+p)(T^2-a_pT+p^2),
 \qquad a_{149}=-248,\quad a_{151}=-244.                     \tag{1}
\]

**Established literature applied.** Tate's theorem for elliptic K3 surfaces
identifies the divisor ranks from these eigenvalues; see
[Schuett--Shioda, section14.6](https://arxiv.org/pdf/0907.0298).
The quadratic factor contributes no root-of-unity eigenvalue: otherwise
`a_p/p` would be a rational algebraic integer, whereas its denominator is
nontrivial. Hence the arithmetic Picard rank is19, the geometric rank is20,
and every geometric divisor class is defined over `Fp^2`.

The inherited arithmetic lattice `N=U+(-M17)` has determinant1092. An
equal-rank integral enlargement could only have index2. Its unique possible
order-two discriminant class has square `-45`, so it cannot lie in the even
Picard lattice of a K3. Thus the arithmetic Picard lattice, and therefore
the full `Fp(t)` Mordell--Weil lattice, are exactly the inherited ones.
There is no torsion because the height of every nonzero section on a24-I1
K3 is `4+2(P.O)>0`.

**Independent new check.** At each prime, compute the heights of all17
basis sections and all136 pair sums using only the rational-function pole
orders. Both reconstructed Grams equal the saved integral Gram exactly.
For any section other than `O`, in the original short chart,

\[
 h(P)=4+\deg\operatorname{den}(X(P))+
 \max(0,\deg\operatorname{num}(X(P))-
              \deg\operatorname{den}(X(P))-4).              \tag{2}
\]

The checker also verifies that the corresponding `Y` pole orders are in
ratio3:2. This pins specialization of every generic word, not just the
fifteen tested trace words. No new finite-field extension counts are made;
the completed exact counts underlying (1) remain explicit dependencies.

At157 the degree24 discriminant has a derivative gcd of degree3. This
prime does not satisfy the24-I1 hypothesis. Its earlier modular exclusions
remain valid under the old individual lifting gates, but the uniform
theorem here must not be applied to it.

## 2. What could happen to a norm-ten curve?

**New geometric lemma.** Work over an algebraically closed field of
characteristic greater than3 on a K3 with24 irreducible fibres. For a
section word of height10 put

\[
 C=2O+4F+\phi(w),\qquad C^2=-2,\quad C.F=2,\quad C.O=0.
\]

Riemann--Roch makes `C` effective; `-C` cannot be effective because its
fibre degree is negative. There is no vertical component. Indeed, subtract
`nF`, `n>=1`. An irreducible horizontal curve would have square
`-2-4n<-2`, contradicting adjunction. If the horizontal part is a sum of
two sections `S_u+S_v`, then

\[
 u+v=w,\qquad h(u)+h(v)=8-2n\le6.
\]

Two nonzero sections each have height at least4; if one is zero, the other
has height10. Both alternatives are impossible. A doubled section is also
impossible since `h(2u)=10` and the section lattice is even. Irreducibility
of the surface fibres accounts for every possible vertical component.

Consequently `C` is either a smooth irreducible genus-zero curve or a pair
of distinct sections. In the latter case necessarily

\[
 h(u)=h(v)=4,\quad\langle u,v\rangle=1,\quad
 S_u.S_v=1,\quad h(u-v)=6.                                 \tag{3}
\]

Both alternatives have a unique effective divisor in their class. The
irreducible `(-2)` curve is rigid; for the pair, intersection of `C` with
each component is `-1`, forcing both components to be fixed. Neither
alternative contains `O`, so `C.O=0` means genuine disjointness from `O`.

In particular the trace relation

\[
 C+P_{-w}=3O+9F                                             \tag{4}
\]

has a unique RR line in any valid minimal base chart. This already proves
rank19 of the19-by20 interpolation matrix after its trace pole degree is
made maximal. It does not require a full geometric MW basis.

## 3. Frobenius excludes the broken alternative for all40,917 classes

**New deduction.** Let `C` be the reduction of a minimum norm-ten class.
If both components in (3) were rational over `Fp`, they would belong to
the complete inherited lattice. But

\[
 u-v\equiv u+v=w\pmod{2M},\qquad h(u-v)=6,
\]

contradicting the certified minimum10 of that parity class. Thus any
splitting would exchange the components under Frobenius. The section
`u-v` would be Frobenius-anti-invariant and have height6.

The anti-invariant geometric MW space has dimension one by (1). Let `d`
be its positive rational height-form squareclass. Orthogonal decomposition
of the geometric Picard space gives

\[
 |\operatorname{disc}NS(X_{\overline{\mathbf F}_p})|
       \equiv1092d\pmod{\mathbf Q^{\times2}}.
\]

Over `Fp^2`, all20 divisor classes are defined. The reciprocal Frobenius
polynomial and Artin--Tate, with its square Brauer order, give

\[
 |\operatorname{disc}NS|\equiv4p^2-a_p^2
                    \pmod{\mathbf Q^{\times2}},\qquad
 d\equiv\frac{4p^2-a_p^2}{1092}.                            \tag{5}
\]

**Established literature applied.** The normalization and square-order
ambiguity are the [Artin--Tate formula](https://kskedlaya.org/weil-cohom/chapter-18.html),
applied over the squared field, not over `Fp`. There is no assertion that
the geometric discriminant or the squared-field Brauer order is determined
exactly by (5).

**Verified exact values.**

| Prime | `4p^2-a_p^2` | Anti-height squareclass representative | Height6 possible? |
| ---: | ---: | ---: | --- |
|149|27300|25, hence1|No: `6/25` is nonsquare|
|151|31668|29|No: `6/29` is nonsquare|

Every nonzero rational vector in a one-dimensional quadratic space has
the same norm squareclass. The required height6 vector therefore cannot
exist. This rules out the last splitting alternative, simultaneously for
all40,917 classes. The unique specialized effective divisor is a smooth
geometrically irreducible conic. Closure in the good smooth K3 model then
gives smooth reduction of each characteristic-zero conic at these primes.

This uses the completed minimum-parity census but does not construct or
evaluate all its curves. The new information is a universal obstruction to
their degeneration, not a new enumeration or a seed-incidence conclusion.

## 4. Four charts and a projective target formula

**New constructive deduction.** A height10 trace satisfies `P.O=3` over
either good finite function field. Among the four distinct base points

\[
 \infty,\quad1,\quad2,\quad3,
\]

at least one is not an intersection with `O`. Move that point to infinity.
In the displayed order, use the original chart first and then

\[
 t=a+1/u,\quad X'=u^4X,\quad Y'=u^6Y,
 \quad A'=u^8A(a+1/u),\quad B'=u^{12}B(a+1/u),
 \quad a=1,2,3.                                            \tag{6}
\]

All three intersections are now finite, so the reduced `X'` denominator
has degree6. The [maximal-pole lifting lemma](DET1092_MODULAR_BISECTION_INCIDENCE_GATE_2026-09-09.md)
applies. The bounds on the pole polynomial and numerators remain3,10,15;
the RR line degrees remain9,5,3. Its rank19 is guaranteed by section2.
The chart choice tests at most four pole degrees and performs only one RR
solve. The proof is geometric; the fixed chart ordering is not claimed
to be invariant under changes of base coordinate.

At the original target `t=0`, the transformed targets are `u=-1/a`, finite
and nonzero. More generally a chosen chart whose infinity is the target
requires a projective evaluation; the present four-chart target-zero
implementation avoids that case explicitly.

Write the specialized RR line as `l0+l1 X+l2 Y=0`.

- If the trace is finite, `l2!=0`: otherwise the vertical line would leave
  `O` in the residual divisor, contradicting `C.O=0`. Divide out the finite
  trace root from cubic elimination, as in the previous certificate.
- If the trace is `O`, then `l2=0`, `l1!=0`. The residual fibre is instead
  exactly
  \[
  X=-l_0/l_1,\qquad Y^2=X^3+A(t_0)X+B(t_0).                \tag{7}
  \]

The line cannot vanish identically on a whole fibre because (4) has no
vertical component. The line at infinity is also impossible: it would
leave a residual `2O`. Thus these cases cover every good target. Formula
(7) is checked symbolically; no old-panel trial needed that branch.

A nonzero nonsquare residual value certifies no local point. A nonzero
square gives two simple local branches by Hensel lifting, not rational
branches over `Q`. A repeated reduction remains UNKNOWN for nearby lifts.

## 5. Uniform local balance, not a global seed discriminator

**New corollary.** Over `Fp`, each conic is `P1`, and its separable degree-two
map to the original base has either zero or two rational branch values.
Write their number as `b`. Counting rational source points gives

\[
 \#\{\text{two-point fibres}\}
 =\#\{\text{zero-point fibres}\}
 =\frac{p+1-b}{2},\qquad b\in\{0,2\}.                    \tag{8}
\]

Indeed `2s+b=p+1` counts source points, while `s+b+n=p+1` counts base
points. Thus each cover has74 or75 nonsplit projective base residues at149,
and75 or76 at151. This is a statement about all base residues, including
any with singular original elliptic fibre; it is not a conditional count
restricted to good elliptic fibres. It uses no independence hypothesis.

On unramified residue classes, smooth proper reduction and Hensel lifting
make the split/nonsplit decision exact over `Q_p`. Near a branch residue,
higher precision may be needed. This explains why an ample rational conic
supply still has many exact local incidence misses. It does **not** predict
which fixed covers hit302 over `Q`, assume independent behaviour across
covers or primes, distinguish productive from inherited rational splits,
or explain single-seed amplification.

## 6. Bounded verification and next obstruction

**Verified application.** Retain exactly the fifteen old masks, target0,
and primes149,151. The full frame uses17 individual heights and136 pair
heights per prime. All30 RR trials have rank19. Orbit61 at149 changes from
`UNKNOWN_MAXIMAL_POLE_GATE` to a simple local split in chart `t=1+1/u`.
The three global UNKNOWN masks61,107,111 remain unchanged; so do all prior
nonsquare conclusions. No additional mask or prime was admitted.

Producer internal runtime0.273s; independent replay below one second. Each
process was capped at25 seconds. The independent implementation uses
manual reversed-order addition, its own two-coordinate height check,
saved-minor/kernel verification, polynomial division and complete residue
root tests. It also replays (1),(5), the index-two parity obstruction and
the vertical-line identity. It does not re-count the old extension fields
or formally verify the geometric arguments. Repeat replay is stable.

The original version finished the arithmetic but failed while serializing
a Sage integer into `geometry.json`. Its source,30 orbit checkpoints,
partial geometry and [failure record](../../artifacts/generated-results/elliptic-curves/det1092_uniform_modular_rr_v1/failure.json)
are preserved. Version2 changes only integer conversion and pre-encodes
JSON before file creation; the mathematical inputs and selection are unchanged.

The seed-construction goal remains open. A complete atlas incidence test
is now supported by a uniform geometric preflight, but it has not been
launched. Even a complete test at these two primes would leave rational
incidence unresolved for the survivors. No new authority for a larger
campaign is inferred from this small prototype.

### Evidence and replay

- [Frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_uniform_modular_rr_v2/protocol.json)
- [Geometry and full Grams](../../artifacts/generated-results/elliptic-curves/det1092_uniform_modular_rr_v2/geometry.json)
- [Thirty-trial summary](../../artifacts/generated-results/elliptic-curves/det1092_uniform_modular_rr_v2/summary.json)
- [Independent replay](../../artifacts/generated-results/elliptic-curves/det1092_uniform_modular_rr_v2/independent-replay.json)

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_uniform_modular_rr.sage
```

Checker SHA256:
`8eab6c38a4c5fa0992a459f8a9906668edb4a2f19442598890f50019761067de`.
Independent replay SHA256:
`82bd786a8b74d6d05837bb4ce3c7f3faaa59570b7a70b492864515b46c603153`.
The uniform result is a new deduction/application of the cited standard
geometry, not a claim of literature-wide novelty.
