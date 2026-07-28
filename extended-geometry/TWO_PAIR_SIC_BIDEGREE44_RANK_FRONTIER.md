# The coefficient-rank frontier in bidegree \((4,4)\)

## 1. Certified answer

Write
\[
 F=\sum_{i,j=0}^4c_{ij}
 \xi _1^i\xi _2^{4-i}z_1^jz_2^{4-j},
 \qquad C=(c_{ij}).
\]
The principal complexity invariant here is the ordinary matrix rank of
\(C\), not expanded support.

The present exact answer is the interval
\[
 \boxed{2\leq \operatorname{rank}C_{\min}\leq5.} \tag{1.1}
\]
The lower endpoint is a theorem: no rank-one bidegree-\((4,4)\) SIC
counterexample exists. The upper endpoint is attained by the known
sixteen-term witness, whose coefficient determinant is \(48\). Ranks two,
three, and four remain open. Thus (1.1) is a rigorous lower bound, not a
claim that the exact minimum is known.

No argument here uses local analysis or deformation around the existing
full-rank point.

## 2. Rank one is impossible

If \(\operatorname{rank}C=1\), then
\[
 F=A(\xi _1,\xi _2)P(z_1,z_2) \tag{2.1}
\]
for binary quartics \(A,P\). Put
\(\Lambda=A(\partial_{z_1},\partial_{z_2})\). Then
\[
 \mathcal E_2(F^m)=\Lambda^m(P^m). \tag{2.2}
\]
Over an algebraic closure every binary quartic \(A\) splits into linear
factors. The
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) says that vanishing of
(2.2) for every \(m\geq1\) implies
\[
 \Lambda^m(QP^m)=0\qquad(m\gg0) \tag{2.3}
\]
for every fixed coordinate polynomial \(Q\).

This gives the full SIC multiplier statement. Write an arbitrary
multiplier as \(g=\sum_\nu B_\nu(\xi)Q_\nu(z)\). Then
\[
 \mathcal E_2(gF^m)
 =\sum_\nu B_\nu(\partial_z)
       \left(\Lambda^m(Q_\nu P^m)\right). \tag{2.4}
\]
Every term is eventually zero by (2.3).

> **Theorem 2.1.** Every bidegree-\((4,4)\) SIC counterexample has
> coefficient-matrix rank at least two.

This is an all-order characteristic-zero theorem, not a bounded moment
exclusion, and it uses no assumption on the root partition of \(A\).

## 3. Determinantal rank charts

For \(1\leq r\leq4\), the rank-at-most-\(r\) locus is parametrized by
\[
 C=UV^{\mathsf T},\qquad
 U,V\in\operatorname{Mat}_{5\times r}. \tag{3.1}
\]
The generic \(\mathrm{GL}_r\) gauge
\[
 (U,V)\longmapsto(UG,VG^{-\mathsf T}) \tag{3.2}
\]
has dimension \(r^2\), so
\[
 \dim\{C:\operatorname{rank}C\leq r\}=r(10-r). \tag{3.3}
\]

| rank bound | factor parameters | gauge | dimension |
|---:|---:|---:|---:|
| \(1\) | \(10\) | \(1\) | \(9\) |
| \(2\) | \(20\) | \(4\) | \(16\) |
| \(3\) | \(30\) | \(9\) | \(21\) |
| \(4\) | \(40\) | \(16\) | \(24\) |

These are global determinantal charts and do not privilege the full-rank
counterexample.

## 4. Pure and mixed moment equations

Introduce
\[
 \Phi_C(x,y)=\sum_{i,j=0}^4c_{ij}x^iy^j. \tag{4.1}
\]
Direct contraction gives
\[
 \boxed{\mu_m(C)=
 \sum_{I=0}^{4m}(4m-I)!\,I!\,
 [x^Iy^I]\Phi_C(x,y)^m.} \tag{4.2}
\]
After (3.1), this is homogeneous of degree \(2m\) in \(U,V\) and invariant
under the gauge (3.2).

For
\[
 Q_{ab}=\xi _1^a\xi _2^{1-a}z_1^bz_2^{1-b},
 \qquad a,b\in\{0,1\},
\]
the bilinear mixed tests are
\[
 \boxed{\nu_{m,a,b}(C)=
 \sum_I (I+a)!\,(4m+1-I-a)!\,
 [x^Iy^{I+a-b}]\Phi_C(x,y)^m,} \tag{4.3}
\]
omitting coefficients outside \(0,\ldots,4m\).

The apolar pairing has diagonal matrix
\[
 D=\operatorname{diag}(24,6,4,6,24). \tag{4.4}
\]
Thus \(L_C=CD\) is, up to transpose convention, the endomorphism attached
to \(C\). On a fixed-flag one-sided chart \(C\) is strictly triangular, so
\[
 \operatorname{tr}(L_C^k)=0\qquad(1\leq k\leq5). \tag{4.5}
\]
A nonzero trace is a useful semistability screen. It is not itself an SIC
defect: one must still exhibit a single fixed \(Q\) whose mixed
contractions fail to vanish eventually.

## 5. Stratum ledger

| exact rank | result | status |
|---:|---|---|
| \(1\) | excluded by Theorem 2.1 | proved |
| \(2\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(3\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(4\) | chart (3.1), tests (4.2)--(4.3) | open |
| \(5\) | explicit witness with \(\det C=48\) | exact counterexample |

The strata are nested: a rank-two counterexample would settle every upper
stratum at once. Conversely, failure to solve finitely many equations on
(3.1) is not a rank lower bound.

An exploratory real rank-two least-squares screen, normalized away from
zero and constrained by nonzero \(\operatorname{tr}(L_C^2)\), produced
points with very small residuals for the first twenty-four pure equations.
They are not counterexamples: there is no exact reconstruction, all-order
identity, or certified nonvanishing mixed tail. This makes exact
reconstruction on the rank-two quotient the first computational target,
not a reason to exclude rank two.

The conditioning audit explains why that screen produced deceptively
small residuals. On a rank-two factor chart, the singular values of the
first-twelve-moment Jacobian at the numerical point range from about
\(3.1\) to \(2.2\cdot10^{-15}\). The seventh normalized residual is the
active one; later raw moments are small largely because of scale. The
numerical point is therefore not suitable for exact reconstruction as it
stands.

## 6. Exact rank-two invariant-quotient obstruction

There is nevertheless an exact reason for truncated semistable rank-two
leads. Let
\[
 X_2=\{C\in\operatorname{Mat}_5:\operatorname{rank}C\leq2\}.
\]
It has dimension \(16\), and its generic diagonal
\(\mathrm{SL}_2\)-quotient has dimension \(13\). On the global chart with
pivot rows of \(U\) equal to the identity, take
\[
 U=\begin{pmatrix}
 1&0\\0&1\\14&17\\18&4\\6&13
 \end{pmatrix},\qquad
 W=\begin{pmatrix}
 8&10&1&8&4\\19&1&4&6&17
 \end{pmatrix},\qquad C=UW.
 \tag{6.1}
\]
Exact evaluation modulo \(1000003\) gives
\[
 \operatorname{rank}
 d(\mu_1,\ldots,\mu_{13})|_C=13. \tag{6.2}
\]
Thus the first thirteen moments are algebraically independent on \(X_2\)
and attain the invariant-quotient dimension.

They cannot be a homogeneous system of parameters. The degree-\(n\)
coordinate ring of the rank-two determinantal variety has the Cauchy
decomposition
\[
 \mathbb Q[X_2]_n
 =\bigoplus_{\substack{\lambda\vdash n\\\ell(\lambda)\leq2}}
 S_\lambda(\operatorname{Sym}^4)^*
 \otimes S_\lambda(\operatorname{Sym}^4). \tag{6.3}
\]
For an \(\mathrm{SL}_2\)-module with even weights, invariant multiplicity
is weight-zero multiplicity minus weight-two multiplicity. Applying this
to (6.3) gives the invariant Hilbert coefficients
\[
 1,1,5,13,53,149,483,1274,3370,7994,18398,39472,81962,
 161896,\ldots. \tag{6.4}
\]
The invariant ring is Cohen--Macaulay: \(X_2\) is the
\(\mathrm{GL}_2\)-quotient of the polynomial factor space \((U,W)\), and
the diagonal \(\mathrm{SL}_2\) action commutes with that gauge, so the
ring is an invariant ring of the reductive product
\(\mathrm{GL}_2\times\mathrm{SL}_2\) acting on a polynomial ring.

If degrees \(1,\ldots,13\) were parameters, the Hilbert numerator would
have nonnegative coefficients. Exact expansion instead gives
\[
 \boxed{
 [t^{69}]H_{\mathbb Q[X_2]^{\mathrm{SL}_2}}(t)
 \prod_{m=1}^{13}(1-t^m)=-5266.} \tag{6.5}
\]
Consequently:

> **Proposition 6.1.** The common zero fiber of
> \(\mu_1,\ldots,\mu_{13}\) on \(X_2\) contains a semistable point.

This proposition does not prove that the extra point has exact rank two:
it could lie on the rank-one boundary, whose uniform finite
quartic moment--nullcone cutoff is not known. Nor does a zero of thirteen
moments satisfy the all-order SIC premise.

The boundary problem is itself tightly calibrated. On the rank-one Segre
cone, \(\mu_1,\ldots,\mu_6\) have exact Jacobian rank six, the dimension
of its generic invariant quotient. Through degree \(100\), their required
Hilbert numerator is
\[
\begin{aligned}
1+t^2+t^3+4t^4+2t^5+7t^6+5t^7+8t^8+5t^9
 +7t^{10}+2t^{11}\\
{}+4t^{12}+t^{13}+t^{14}+t^{16},
\end{aligned} \tag{6.6}
\]
with zero coefficients after degree \(16\) through the checked cutoff. It
is nonnegative and has coefficient sum \(50\). Thus the first six moments
are a Hilbert-compatible parameter candidate on rank one.
The remaining geometry is exactly the finite set of exceptional
squarefree-quartic cross-ratios left by the generic theorem, together with
a uniform finite cutoff on the at-most-two-root orbits. A direct generic
Gröbner basis over \(\mathbb Q(\lambda)\) did not finish within 180
seconds; a Fitting determinant or exceptional-fiber resultant is the
better next calculation.

The least degree repair is
\[
 \boxed{\mu_1,\ldots,\mu_{12},\mu_{14}.} \tag{6.7}
\]
At the same point (6.1), these corrected moments also have exact Jacobian
rank thirteen. Their Hilbert numerator is nonnegative through degree
\(100\), with last observed nonzero coefficient in degree \(82\). This is
a necessary-test result, not a proof that (6.6) is a system of parameters.

## 7. Next exact gates

The efficient order is:

1. start with rank two, whose determinantal quotient has expected affine
   dimension \(16-3=13\) when the diagonal
   \(\mathrm{SL}_2\)-stabilizer is finite;
2. normalize a non-null lowest Clebsch--Gordan component and remove the
   residual torus before elimination;
3. use the Hilbert-compatible corrected set
   \(\mu_1,\ldots,\mu_{12},\mu_{14}\), rather than the impossible
   consecutive set, and determine its projective zero fiber;
4. apply (4.3) with fixed low-bidegree multipliers on every survivor; and
5. only after a global rank-two exclusion, repeat for ranks three and
   four.

Promotion to a counterexample requires an exact coefficient field,
all-order pure vanishing, and one fixed mixed defect for infinitely many
\(m\). Promotion of the lower bound beyond two requires a global proof on
the entire determinantal variety.

## 8. Reproduction

Run

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_frontier.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_invariants.py
```

The dependency-free checker verifies the chart dimensions, exact
representatives of ranks one through four, the pure and all four bilinear
mixed formulas through order four, the nilpotent trace screen, and the
known rank-five determinant and moment formulas. It does not turn the
exploratory rank-two residuals into a theorem.

The second checker constructs the rank-two moment Jacobian directly,
proves both displayed rank-thirteen statements modulo a good prime, and
computes the invariant Hilbert coefficients from the Cauchy and
Jacobi--Trudi formulas. It verifies the exact obstruction (6.5) and the
corrected necessary test through degree \(100\).
