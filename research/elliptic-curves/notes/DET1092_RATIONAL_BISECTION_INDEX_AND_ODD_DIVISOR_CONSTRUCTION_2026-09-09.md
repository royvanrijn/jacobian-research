# Determinant 1092: rational bisection index and constructive descent

## Result and scope

**New deduction; verified application.** All **40,917** geometrically rational
bisection translation orbits in the completed degree-two quotient have
representatives isomorphic to **`P1` over `Q`**, not merely over an extension.
The generic lattice forces an odd-degree rational divisor on each conic.
This strengthens the [completed geometric census](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md);
it does not repeat its CVP enumeration.

Two additional covers, plus the old orbit8044 regression, now have explicit
rational parametrizations and elliptic maps. Each separately gives rank
at least18 over `Q(u)`. All three are nonsplit at302 and each of the eight
unchanged control parameters, with exact local certificates. Rationality of
the auxiliary curve is therefore settled for this supply; **rational incidence
above the requested parameter is not**. Prospective302 member selection is open.

The written proof and algebraic applications below are not formally verified
or externally reviewed. An independent Sage checker replays the applications.

## 1. Uniform index-one argument

**Verified input.** The parent has24 irreducible `I1` fibres and full
rational/geometric MW lattice `M` of rank17, with even integral Gram matrix
`G`, determinant1092 and `NS(X)=U+(-M)`. The complete census gives40917 norm10
minimum parity representatives `w`; each divisor

\[
 C_w=2O+4F+\phi(w),\qquad C_w^2=-2,\quad C_w.F=2
\]

is a geometrically irreducible smooth genus-zero curve defined over `Q`.
For the rational basis section `S_i=O+(G_ii/2)F+phi(e_i)`,

\[
 C_w.S_i=G_{ii}-(Gw)_i.                                      \tag{1}
\]

**Verified application.** Exact reduction of the saved Gram matrix gives

\[
 \operatorname{rank}_{\mathbf F_2}(G)=16,\quad
 r=(1,1,0,0,0,0,0,1,0,1,0,1,1,1,1,0,1),\quad
 Gr\equiv0\pmod2,\quad r^tGr=180.
\]

**New deduction.** If `Gw` were even, `w mod2` would be zero or `r`. Norm
modulo4 is unchanged on adding twice an integral vector, since

\[
 (v+2x)^tG(v+2x)-v^tGv=4x^tGv+4x^tGx.
\]

Both radical classes have norm0 modulo4. A norm10 vector cannot belong to
the radical. Thus some `(Gw)_i` is odd, so (1) is an odd positive integer
`2m+1`. The rational line bundle

\[
 \left.\mathcal O_X(S_i-mF)\right|_{C_w}
\]

has degree1. Genus-zero Riemann--Roch gives an effective rational divisor of
degree1, hence a rational point and `C_w ≅ P1_Q`.

**Established literature.** The odd-degree-line-bundle characterization is
[Stacks, Proposition53.10.4](https://stacks.math.columbia.edu/tag/0C6U).
Here the lattice supplies its hypothesis; no displayed basis section is
assumed to intersect the conic in degree1.

The proof is basis invariant: it concerns the radical on `M/2M` and its norm
modulo4. Rational section translations preserve the conclusion. It covers
all40917 orbits without constructing all their equations. This property of
the common parent does **not distinguish302 from its controls**.

## 2. Frozen bounded construction

**Verified application.** Before constructing or testing new covers, select
exactly the three rows tied at minimum stored coefficient `l1=3` in the
immutable quotient table, ordered by orbit mask. In the saved one-based basis:

| Orbit | Stored word | First odd intersection used | Descent degrees |
| ---: | --- | --- | --- |
|8044|`-e2+e10-e14`|`C.S1=7`|`7 → 5 → 3 → 1`|
|47755|`e2+e12-e16`|`C.S3=3`|`3 → 1`|
|103186|`e13-e16-e17`|`C.S1=7`|`7 → 5 → 3 → 1`|

This is **stored-coordinate selection**, not a basis-invariant scheduler or
minimization over every representative. The uniform theorem is intrinsic;
the example selection is not. Later sections meet47755 in degree1, but the
frozen first-odd rule was retained.

For each new row the generic RR interface solves one19-by20 system in
`H0(3O+9F)` through `P_{-w}`. Its unique relation is

\[
 f_0(t)+f_1(t)x+f_2(t)y=0,\qquad
 (\deg f_0,\deg f_1,\deg f_2)\le(9,5,3).
\]

Remove the known trace factor from elliptic elimination. For the residual
`a x²+b x+c`, polynomial factorization gives

\[
 b^2-4ac=h(t)^2q(t),\qquad \deg q=2,\quad\operatorname{disc}(q)\ne0.
\]

The full rational scalar is retained in `q`; deleting a nonsquare scalar
would change the conic. The exact elliptic lift is

\[
 W^2=q(t),\quad x=x_0+x_1W,\quad y=y_0+y_1W,
\]
\[
 x_0=-b/(2a),\quad x_1=h/(2a),\quad
 y_0=-(f_0+f_1x_0)/f_2,\quad y_1=-f_1x_1/f_2.                \tag{2}
\]

All four coefficients in (2) are polynomials, with `gcd(x1,y1)=1`, in each
of the three examples. Old8044 equations were reused; only47755 and103186
required new RR systems. No exceptional point or later cascade artifact
is an input.

### Constructing a rational point by polynomial arithmetic

**New constructive deduction; verified on these three examples.** Polynomial
Bezout for `x1,y1`, followed by numerator gcds at the first odd-intersection
generic section, gives a monic `g(t)` of odd degree `d=C.S_i` and an ordinate
`v(t)` of degree below `d`, with `g | v²-q`. Replay verifies the section
incidence modulo `g` and that coordinate denominators are coprime to `g`.

For `d>1`, replace

\[
 g\longleftarrow\operatorname{monic}((v^2-q)/g),\qquad
 v\longleftarrow v\bmod g.                                  \tag{3}
\]

The new degree is positive, odd and at most `d-2`: divisibility and the
nonsquare quadratic exclude `deg(v)≤1` for `d≥3`, and the quotient degree
is `2 deg(v)-d`. Thus (3) terminates at a linear `g`. Its root `t0`, with
`W0=v(t0)`, is a rational conic point. These are polynomial calculations,
with no integer factorization, black-box conic solver, point search or class/unit group.
Incidence extraction was replayed on these three examples; no all40917
equation/implementation claim is made.

Writing `q=q2 t²+q1 t+q0`, a rational parametrization is

\[
 T(u)=t_0+\frac{q'(t_0)-2W_0u}{u^2-q_2},\qquad
 S(u)=W_0+u(T(u)-t_0).                                      \tag{4}
\]

The saved functions have `deg(T)=2` and satisfy `S²=q(T)` exactly. Equations
(2),(4) give the maps into the parent. The linked certificates contain every
coefficient, rational point and descent step; arrays use ascending coefficient
order. A sufficient finite excluded set consists of zeros of the reduced
denominator of `T` and numerator of the parent discriminant `Delta(T)`.
Outside it the polynomial elliptic maps are defined on smooth finite parent
fibres. This supplies rational points on the cover over `Q`, with no constant
number-field extension.

**Verified application of a completed theorem.** The
[all-prime multisection theorem](DET1092_ALL_PRIME_DIVISION_AND_MULTISECTION_THEOREM_2026-09-09.md)
makes each genuine geometrically irreducible bisection independent of the
inherited generic rational span. Each degree-two base change therefore gives
rank at least18 over `Q(u)`. This does **not** prove three jointly independent
extras or independence at every rational specialization. The old dependent
split remains a counterexample to the latter claim.

## 3. Exact negative exposure and interpretation

**Verified application.** After construction, evaluate `q(t)` only at the
unchanged302/eight-null original parent parameters. All27 values are nonsquare.
At302 the saved scalar convention gives:

| Orbit | Prime | Valuation of `q(0)` | Leading unit | Obstruction |
| ---: | ---: | ---: | ---: | --- |
|8044|29|2|21|Nonsquare residue|
|47755|5|59|1|Odd valuation|
|103186|7|23|3|Odd valuation|

All27 witnesses and exact values are in `controls.json`. Large valuations
are valid: the rational scalar is not square-normalized. The independent
checker uses finite residue squares, not the producer's square-test result.
Rational section translation preserves the base map and splitting. Thus the
misses exclude the **entire three translation orbits** on this panel, not
only their representatives. They do not exclude the other40914 orbits.

**New deduction / explanatory boundary.** Three gates must remain separate:

1. Rationality of the multisection itself: settled positively for all40917.
2. Rational incidence above the requested parameter: the square condition
   `q_w(t) ∈ Q²`, explicitly false for these three covers at302.
3. Independence of a resulting specialized point: a further requirement,
   handled for specified points by the completed
   [halving-or-cycle classifier](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md),
   not by conic splitting alone.

An ample rational geometric supply can therefore coexist with the control
misses. This does not explain why302 has its first seed or its amplification.
Whether a rational bisection orbit supplies a new302 seed is **UNKNOWN**;
the known genus-one carrier does not answer that question. No conjectural
implication from these results to a fibre rank upper bound is asserted.

## 4. Certificates and replay

Immutable packet:
[`det1092_rational_bisection_index_v1`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/).

- [`protocol.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/protocol.json)
  and [`index-and-selection.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/index-and-selection.json):
  input hashes, radical and frozen selection.
- [`orbit-8044.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/orbit-8044.json),
  [`orbit-47755.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/orbit-47755.json),
  [`orbit-103186.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/orbit-103186.json):
  exact equations and maps, with separate preconstruction protocols.
- [`controls.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/controls.json)
  and [`independent-replay.json`](../../artifacts/generated-results/elliptic-curves/det1092_rational_bisection_index_v1/independent-replay.json):
  all27 negative exposures and independent checks.

The independent checker imports no producer and does not repeat rational
kernel or factor discovery. It verifies hashes, Gram data, table selection,
RR uniqueness via kernel and a modular rank19 lower bound, trace elimination,
conic/elliptic identities, odd divisor incidence, every descent step, rational
maps and all local obstructions. Generic independence uses the cited theorem.

From the repository root:

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/audit_det1092_rational_bisection_index.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/construct_det1092_minimal_stored_bisections.sage 8044
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/construct_det1092_minimal_stored_bisections.sage 47755
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/construct_det1092_minimal_stored_bisections.sage 103186
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/audit_det1092_minimal_bisection_controls.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_rational_bisection_index.sage
```

Each internal computation finished below one second in Sage10.9; every
process is capped at25 seconds. A verifier preflight reached report creation
with an undefined selection-path constant. That failure is retained; the
constant was fixed and the successful immutable replay repeated. No
mathematical selection or equation changed. There was no new parameter/point
search, CVP census, broad descent, production mutation or detached job.
