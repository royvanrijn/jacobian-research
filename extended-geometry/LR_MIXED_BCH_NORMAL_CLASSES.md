# Balanced mixed BCH classes in the filtered LR quotient

This note carries the rooted-tree compiler of
[rooted-tree LR normal classes](LR_ROOTED_TREE_NORMAL_CLASSES.md) into the
actual mixed BCH series of the degree-five map \(F_2\).  It proves that one
balanced multihomogeneous BCH coefficient has a nonzero saturated
associated-graded normal residue in every odd order.

The result closes the tree-cancellation issue for this sector, but not the
universal lower-jet problem.  The coefficient is present only when both
opposite-weight target amplitudes are active.  It vanishes on the valid
lower-jet locus where either amplitude is zero, so it is an all-order penalty
for using those directions rather than an unconditional LR obstruction.

The exact certificate is
[`lr_mixed_bch_classes.json`](../artifacts/generated-results/lr_mixed_bch_classes.json).

## 1. The mixed word

Use the map-composition pre-Lie product and bracket

\[
 P\mathbin{\triangleright}Q=(DP)Q,\qquad
 \{P,Q\}_{\rm map}=P\mathbin{\triangleright}Q-
                    Q\mathbin{\triangleright}P.
                                                               \tag{1.1}
\]

Let

\[
 D_B=\ell_{F_2}(\partial_B),\qquad
 D_C=\ell_{F_2}(\partial_C).
\]

Their weights are \(1,-1\).  Since the target fields are constant and
\(\ell_{F_2}\) is a Lie homomorphism,

\[
 \{D_B,D_C\}_{\rm map}=0.                            \tag{1.2}
\]

Normalize the leading weight-zero deformation direction by removing its
irrelevant scalar \(435/7\):

\[
 X=N(x,0,-3z),\qquad N=v^6S^4.                       \tag{1.3}
\]

In logarithmic invariant coordinates
\((\delta x/x,\delta u,\delta\gamma)\), this is

\[
 q(X)=
 \left(
 N,\ vN,\ -\left(\frac87v+S\right)N
 \right).                                            \tag{1.4}
\]

Define

\[
 \boxed{
 W_k=
 \bigl(\operatorname{ad}_{D_B}\operatorname{ad}_{D_C}\bigr)^kX
 =
 \{\{\cdots\{\{X,D_B\},D_C\}\cdots,D_B\},D_C\}.
 }                                                    \tag{1.5}
\]

The two expressions agree with the bracket convention (1.1): each right
bracket contributes a minus sign, and a \(B,C\) pair contributes two.

## 2. Why this is an actual BCH coefficient

Introduce a bookkeeping parameter \(\epsilon\) on \(X\).  The part of

\[
 \operatorname{BCH}_{\rm map}(-\epsilon X,-D),
 \qquad D=sD_B+tD_C,
\]

which is linear in \(\epsilon\) is

\[
 -\epsilon\,
 \frac{\operatorname{ad}_D}
      {1-\exp(-\operatorname{ad}_D)}X.               \tag{2.1}
\]

Consequently the term with \(2k\) copies of \(D\) has coefficient

\[
 -\frac{B_{2k}}{(2k)!}\operatorname{ad}_D^{2k}X.
                                                               \tag{2.2}
\]

Because (1.2) implies that
\(\operatorname{ad}_{D_B}\) and \(\operatorname{ad}_{D_C}\) commute, the
coefficient of \(s^kt^k\) in (2.2) is the single word

\[
 \boxed{
 -\frac{B_{2k}}{(2k)!}
 \binom{2k}{k}s^kt^k W_k.
 }                                                    \tag{2.3}
\]

Thus all rooted-tree placements in this multihomogeneous sector have already
been summed.  There is no remaining tree cancellation inside (2.3).
Moreover \(B_{2k}\ne0\) for every \(k\ge1\), so the scalar in front of
\(W_k\) never vanishes.

For \(k=1\), (2.3) gives \(-W_1/6\), agreeing with the order-three term
\(\{\{X,D_B\},D_C\}\) predicted by the filtered-coset BCH lemma.

## 3. Exact boundary-face transfer

Put

\[
 w=u\gamma.
\]

After the first mixed pair, the relevant boundary face of a weight-zero field
has the form

\[
 q=
 \left(
 u^{d+1}F(w),\
 u^{d+2}F(w),\
 u^dH(w)
 \right).                                            \tag{3.1}
\]

The exact affine-coordinate compiler, not a coordinate-invariant surrogate,
shows that the operator

\[
 T(V)=\{\{V,D_B\},D_C\}
\]

maps (3.1) to

\[
 \left(
 u^{d+5}F_{\rm new}(w),\
 u^{d+6}F_{\rm new}(w),\
 u^{d+4}H_{\rm new}(w)
 \right).                                            \tag{3.2}
\]

The complete two-by-two polynomial operator is stored in the certificate.
At the point \(w=0\), its triangular part is

\[
\begin{aligned}
 F_{\rm new}(0)
 &=-540(d+1)
   \left(17(d-1)F(0)+317H(0)\right),\\
 H_{\rm new}(0)
 &=-9180(d+1)(d+3)H(0).                              \tag{3.3}
\end{aligned}
\]

The second formula is independent of \(F\) and of every positive \(w\)-jet.
It therefore closes the third normal residue by itself.

For \(W_k\), \(k\ge1\), one has

\[
 d_k=4k+11,
\qquad
 [q_\gamma(W_k)]_{\gamma=0}
 =c_ku^{d_k}+\text{lower powers of }u,               \tag{3.4}
\]

with seed

\[
 c_1=\frac{14438891520}{2401}.                       \tag{3.5}
\]

Substituting \(d_k\) into (3.3) gives

\[
 \boxed{
 c_{k+1}=-73440(k+3)(2k+7)c_k\ne0.
 }                                                    \tag{3.6}
\]

This is an all-order proof.  The explicit calculations through \(k=3\), or
BCH order seven, are regression checks only.

### Proposition 3.1

For every \(k\ge1\), the coefficient linear in \(X\) and of bidegree
\((k,k)\) in \((D_B,D_C)\) in
\(\operatorname{BCH}_{\rm map}(-X,-D_B-D_C)\) has nonzero image in the third
saturated normal summand \(R/(\gamma)\).  Its associated-graded symbol has
\(u\)-degree \(4k+11\).

### Proof

Equation (3.6) makes the leading coefficient of \(W_k\) nonzero.  Equation
(2.3) multiplies it by a nonzero Bernoulli/binomial scalar.  The third row of
the logarithmic differential matrix is \((\gamma,0,1)\), so this leading
coefficient is exactly the leading term of the third saturated normal
residue.  QED

## 4. Descent audit

There are two distinct descent questions.

### 4.1 Linear target descent succeeds

The target image in the third saturated normal coordinate is the ideal
\((\gamma)\).  Proposition 3.1 is computed after quotienting by that ideal.
Hence no new linear target correction at order \(2k+1\) kills the displayed
class.  This is a genuine class in the same saturated normal module used by
the quadratic Rees/SAGBI computation.

### 4.2 Universal lower-jet descent fails for this sector alone

With target amplitudes \(s,t\), formula (2.3) contains the factor

\[
 s^kt^k.
\]

It vanishes on \(s=0\) or \(t=0\), including the valid lower-jet choice
\(D_1=0\).  Therefore Proposition 3.1 does not define a nowhere-vanishing
section on the full lower-jet solution scheme \(Z_{2k}\).  It cannot by itself
prove the universal class (6.13)--(6.14) of
[complexity-filtered contact](COMPLEXITY_FILTERED_CONTACT.md).

What it does prove is a sharp dichotomy component:

> Activating both minimal opposite-weight constant target directions incurs
> a nonzero saturated normal penalty in every odd BCH order.

This converts the old possible-cancellation concern into a lower-jet
stratification problem.

## 5. Use in the existing filtered-LR chain

The result belongs in the chain

\[
 \text{target-jet reduction}
 \longrightarrow
 \text{filtered-coset BCH}
 \longrightarrow
 \text{finite target modules}
 \longrightarrow
 \text{lower-jet stratification in `OP-CCDM`}.
\]

It should not become a separate programme.  The next useful steps are:

1. split the lower-jet scheme by whether the minimal opposite-weight
   amplitudes satisfy \(st=0\) or \(st\ne0\);
2. use Proposition 3.1 on the open stratum \(st\ne0\);
3. on \(st=0\), return to the weight-zero kernel and the
   \(v^{6m}S^{4m}\) torus recurrence;
4. prove that higher semi-invariant target generators either reduce to these
   strata or are killed by the finite annihilator cutoff from the
   LR Rees/SAGBI module; and
5. only then test whether the two branches give one universal degree-growth
   alternative.

The missing theorem is now precise: show that every low-resource lower jet
lies either in a branch controlled by the torus recurrence or in a branch
where a mixed coefficient such as Proposition 3.1 is forced to be nonzero.

## 6. Reproduction

Generate the exact boundary operator and regression data with

```bash
.venv/bin/python scripts/compile_lr_mixed_bch_classes.py --max-k 3
```

Replay the leading recurrence and the first Bernoulli/binomial coefficients
using only the Python standard library with

```bash
python3 scripts/audit_lr_mixed_bch_classes.py
```

The generated JSON stores the complete polynomial boundary operator, not
only its triangular specialization (3.3).
