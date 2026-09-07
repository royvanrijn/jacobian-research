# Curve 302: explicit alternative parent with full arithmetic MW17

Authority: `EC-CURVE302-RECOVERED-MW17-PARENT`.

**Geometric closure:** the [two-prime follow-up and search adapter](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md) now prove full geometric MW17 and Picard19, independently reconstruct the full149 Frobenius polynomial, and supply reusable trace tables. The original certificate below is retained unchanged.

**The alternative-parent construction is complete.** The supplied elliptic
surface over `Q(t)` has arithmetic generic Mordell–Weil rank exactly **17**,
trivial torsion, and an explicit full saturated basis of height determinant
**1092**. Its fibre at **t=0 is literally the displayed equation of curve302**.
The seventeen basis images generate the entire previously identified
primitive rank-17 core in the displayed independent rank-31 group `D`.

The family is genuinely different from the completed MW9 K3 and the
determinant948 production surface. Original discoverer provenance remains
`UNKNOWN`; this result constructs an alternative parent and does not identify
an unpublished generation procedure. No new rational-fibre rank is claimed.

## Equation, basis and immediate use

The model is

\[
 E_t:\quad y^2+xy+y=x^3+x^2+a_4(t)x+a_6(t),
 \qquad \deg a_4=8,\quad\deg a_6=12.
\]

Every coefficient of these two polynomials, all seventeen rational-function
point pairs, the exact height Gram and the 31-by-17 specialization matrix are
expanded in
[`curve302_recovered_mw17_parent_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json).
Each rational function is encoded by numerator and denominator coefficient
lists in **ascending degree**, with rational numbers stored as strings.
This is an explicit coefficient list, not an implicit equation-solving task.

From the repository root, in Sage:

```python
load('elliptic-curves/cas/load_curve302_recovered_parent.sage')
E, basis, t0 = load_curve302_recovered_parent()
# t0 == 0; len(basis) == 17
E0 = EllipticCurve(QQ, [a(t0) for a in E.a_invariants()])
images = [E0([c(t0) for c in P.xy()]) for P in basis]
```

The coefficients at zero are exactly

```text
a1 = a2 = a3 = 1
a4 = -1284727764113567728281797636015784768866707681415849262157224232063
a6 = 560368321454261339256859338901915312332769858684945406858043869199456710681989058863306170127006181
```

Thus no final twist or unresolved parameter identification is involved.
The large rational coefficients reflect the recovered quartic coordinates;
a smaller presentation is not required for the proof and remains optional.

## Construction from the recovered height form

The [preceding reconstruction](CURVE302_DET1092_RECONSTRUCTION_2026-09-07.md)
recovered a unique minimum-four determinant1092 form `G` from400 of401
integral-point words. That was a candidate generic height form. The present
construction and exact section heights now realize it.

Choose a degree-four polarization `D=O+P+Q`, with old fibre intersection
`D.F=3`. Fifteen predicted line sections have degree one for this
polarization. Their images on the plane model of302 are obtained from
the degree-three linear system `O+P+Q`. If the specialized plane coordinates
are `(Xi,Yi,1)`, write the candidate line in projective three-space as

\[
 (X_i z+u_i w,\;Y_i z+v_i w,\;z,\;w).
\]

Every required line incidence gives the linear equation

\[
 (X_i-X_j)(v_i-v_j)-(Y_i-Y_j)(u_i-u_j)=0.
\]

The rational kernel has dimension four, including the three trivial
concurrent-line directions. The nontrivial configuration has all27 required
incidences and all78 required nonincidences. Requiring the quartic to contain
the fifteen lines and have the prescribed node gives a unique quartic,
written

\[
 F=C(X,Y,Z)L(X,Y,Z)+WQ_3(X,Y,Z,W),
 \qquad L=X-Y+\frac{727317917032131068613045358848720}{6403254827996423}Z.
\]

The plane pencil `W=tL` leaves the explicit residual cubic

\[
 C(X,Y,Z)+tQ_3(X,Y,Z,tL)=0
\]

with rational origin `(1:1:0)`. A tangent-frame and Cremona transformation
give an explicit Weierstrass model. A constant change with
`(u,r,s,tau)=(-6,15,-3,-108)` normalizes its zero fibre to the literal302
model above.

The fifteen lines and two first conics yield seventeen displayed sections
but only **rank15**. Two further rational conics, selected outside this span
using `G`, lift to sections and raise the rank to17. Seventeen of these
nineteen sections have a unimodular word matrix in the recovered core.
This distinction is checked exactly; the number of displayed sections was
never used as an independence proof.

The frozen quartic and four conics are retained in
[`curve302_recovered_mw17_construction_inputs_v1.json`](../../artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_construction_inputs_v1.json).
The [construction replay](../cas/reconstruct_curve302_recovered_parent.sage)
rebuilds the residual cubic, its anchored Weierstrass transformation and all
seventeen basis coordinates, then compares them with the expanded endpoint.
The final theorem does not assume that an abstract lattice or moduli point
must lift: the resulting equations and points are checked directly.

## Picard capacity and exact generic rank

The short model has degrees `(8,12)` and a squarefree discriminant of degree24.
Its leading discriminant coefficient is nonzero. Therefore its minimal
elliptic surface is a K3 with precisely24 geometric `I1` fibres and smooth
infinity. There are no reducible-fibre height corrections.

The Frobenius gate is applied to this different surface. At `p=149`, use
the integral chart `t=v/(v+1)`, with the usual Weierstrass weights. The reduced
short coefficients, in ascending degree, are

```text
A = [33,100,70,78,93,99,16,80,147]
B = [6,39,148,4,69,116,63,121,85,128,23,22,91]
```

Squarefreeness of `B` and the discriminant, their nonzero endpoint values,
and the binomial face conditions verify nondegeneracy on every Newton face.
The reflexive simplex has primitive Hodge numbers `(1,18,1)`.
The pinned [ToricControlledReduction implementation](https://github.com/edgarcosta/ToricControlledReduction/tree/74cda9e8148cd8e9a3928fc15a558c9a70b67cc1)
computes the primitive degree20 characteristic polynomial

\[
 P_{149}(T)=(T-149)^{17}(T+149)(T^2+248T+22201).
\]

Adding the two ambient rational classes gives multiplicity19 for eigenvalue
149 in full `H2`. Good-reduction divisor specialization therefore bounds
the rational Néron–Severi rank by19. Shioda–Tate bounds arithmetic generic
MW rank by17. The seventeen explicit sections have a positive definite
height Gram, attaining the bound:

\[
 \operatorname{rank}E(\mathbb Q(t))=17,\qquad
 \operatorname{rank}\operatorname{NS}(X_{\overline{\mathbb Q}})^{G_{\mathbb Q}}=19.
\]

This uses a Frobenius **upper bound**, not a Tate-conjecture converse, BSD or
GRH. The `-149` eigenvalue contributes only to the geometric upper bound.
This first reduction alone gives geometric MW bounds17–18 and Picard bounds19–20.
The linked two-prime follow-up now proves exact geometric MW17 and Picard19.

Independent PARI elliptic-fibre counts reproduce:

| Field | Primitive trace | Surface point count |
|---|---:|---:|
| `F149` | 2136 | 24636 |
| `F149^2` | 416720 | 493345524 |

The full controlled-reduction calculation has **one implementation**. These
two independent moments originally corroborated it without reconstructing
all twenty coefficients. The linked follow-up now uses the19 verified divisor
directions and orthogonality to reconstruct the full polynomial from those
counts, independently of controlled reduction. The raw input, raw output, original backend log,
commit and hashes are retained beside the endpoint JSON. The default replay
checks their binding to the exact equation and recomputes both moments;
it does not silently rerun the external backend.

## Full basis and saturation

For a nonzero section on this24-I1 K3,

\[
 \langle P,P\rangle=4+2(P.O)
 =\max(\deg\operatorname{num}x(P),\;\deg\operatorname{den}x(P)+4).
\]

This formula, applied also to every `P_i-P_j`, verifies the entire exact
17-by-17 Gram. Its determinant is1092 and its discriminant group is cyclic
of order1092. The word matrix identifies it with the previously recovered
form, up to a unimodular basis change.

All nonzero sections have positive even integral height, excluding torsion.
Any larger group of the same rank would be an even integral overlattice.
Since `1092=2^2*3*7*13`, its first prime index could only be two. The unique
order-two discriminant class has a representative of norm **45**, which is
odd. Hence no such even overlattice exists. The seventeen sections form the
**full arithmetic generic MW basis**, and their span is also primitive in
the geometric MW group, regardless of the unresolved geometric rank.

The full rational NS lattice has absolute determinant1092. Its rank and
determinant distinguish this parent over `Q` from the completed MW9 K3
(rational rank11) and the production948 K3 (rank19, determinant948).

## What specialization explains

Let `D` denote the displayed independent rank-31 group on302, whose
independence is certified by `ECR31`. The endpoint contains an exact integer
31-by-17 column embedding. Every specialized point is checked against its
indicated combination of the public points. All Smith factors are one,
and the basis change in the previous core is unimodular. Consequently

\[
 \operatorname{sp}(M_{17})=\text{the complete primitive17-core in }D,
 \qquad D/\operatorname{sp}(M_{17})\cong\mathbb Z^{14}.
\]

This establishes seventeen generic directions and fourteen independent
directions in the displayed fibre quotient. It does not prove that the full
rational group of302 has rank exactly31 or that `D` has global index one.
Unlike the completed MW9 example, there is no finite specialization index
defect **inside this displayed group**.

## Reproduction and continuation boundary

With SageMath10.9, from the repository root:

```sh
sage -python elliptic-curves/cas/reconstruct_curve302_recovered_parent.sage
sage -python elliptic-curves/cas/verify_curve302_recovered_mw17_parent.sage
```

Both are bounded single-worker exact replays, with120-second caps. The
[compact proof output](../../artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_proof_v1.json)
pins all arithmetic verifier inputs. To recompute the full external
Frobenius polynomial, build the pinned backend commit and pass the retained
`curve302_recovered_mw17_frobenius_p149_input_txt.txt` to its
`build/examples/readfile.exe`, using a fresh output path. Compare the parsed
polynomial with the retained output; do not overwrite original evidence.

The requested surface/full-basis/302 endpoint is closed. The construction
now supplies a separate calibrated rank-search family, but no production
specialization sweep is part of this certificate. Original provenance,
equation simplification and future fibre search
are separate questions. The stronger different-NS foundry target also asks
for a positive-rank carrier with an independent pullback section; this result
does not by itself close that additional requirement.

The [Elkies Shimura/K3 framework](https://arxiv.org/abs/0802.1301) motivated
the determinant1092 arithmetic lead. The explicit equation and full MW proof
above are verified without selecting an individually identified non-CM
point on its moduli curve.
