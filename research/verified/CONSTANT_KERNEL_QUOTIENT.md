# Constant-kernel quotients and essential cubic input

## The theorem

Work over a characteristic-zero field `k`.  Let

\[
 F=I+H:\mathbb A^n\longrightarrow\mathbb A^n
\]

be a polynomial map and let

\[
 K\subseteq\bigcap_{x\in\mathbb A^n}\ker JH(x)
\]

be a constant `d`-dimensional subspace.  Choose a quotient and section

\[
 B:k^n\to k^{n-d},\qquad \ker B=K,
 \qquad C:k^{n-d}\to k^n,\qquad BC=I.
\]

Define

\[
 f(q)=B F(Cq)=q+B H(Cq).                            \tag{1}
\]

### Constant-kernel quotient theorem

There is a linear isomorphism `S=(C|K)` and a polynomial map
`g:A^(n-d)->A^d` such that

\[
 \boxed{S^{-1}FS(q,r)=(f(q),r+g(q)).}               \tag{2}
\]

Consequently:

1. `det DF=det Df`; hence `F` is Keller exactly when `f` is Keller.
2. If `H` is homogeneous of degree `e`, then the nonlinear part of `f` is
   homogeneous of degree `e`.
3. For every target `(a,b)` in the coordinates of (2), projection to `q`
   gives an isomorphism of fiber schemes
   \[
   (S^{-1}FS)^{-1}(a,b)\simeq f^{-1}(a),
   \qquad r=b-g(q).                                  \tag{3}
   \]
4. If two points in one fiber of `F` have distinct `B`-projections, then `f`
   is noninjective.  Any finite stored collision descends provided its
   projected points remain pairwise distinct.

When `K` is the full intersection of the Jacobian kernels, `f` is the
**essential-input quotient** of this type: its nonlinear part has no remaining
constant translation-invariance direction.

## Proof

For `k in K`, the polynomial in one variable `t -> H(x+tk)` has derivative

\[
 JH(x+tk)k=0.
\]

Characteristic zero makes it constant.  Thus `H` is constant on the affine
`K`-cosets.  Since `I-CB` has image `K`,

\[
 H=HCB.                                               \tag{4}
\]

Write a source point uniquely as `Cq+Kr`.  Applying `S^-1` to
`Cq+Kr+H(Cq)` gives (2), with its first block equal to (1).  Differentiating
(2) gives

\[
 D(S^{-1}FS)=
 \begin{pmatrix}Df&0\\Dg&I_d\end{pmatrix},          \tag{5}
\]

which proves the determinant assertion.  Homogeneity follows from (1).
Finally, the equations of the fiber over `(a,b)` are `f(q)=a` and
`r=b-g(q)`; eliminating `r` proves the scheme isomorphism (3), including
nilpotent structure and multiplicities.

## Relation to Gorni--Zampieri pairing

The matrices `B,C`, the identity `BC=I`, and the induced formula `f=BFC`
place this construction beside Gorni--Zampieri pairing.  Classical GZ pairing
passes between cubic-homogeneous and higher-dimensional cubic-linear maps and
transfers injectivity, surjectivity, polynomial invertibility, and the Keller
property.  See Gorni and Zampieri,
[*On cubic-linear polynomial mappings*](https://arxiv.org/abs/1204.4026).

The theorem above is narrower in one direction and stronger in another: it
does not require the ambient map to be cubic-linear, but assumes a literal
constant translation kernel and therefore obtains the explicit triangular
conjugacy (2) and the fiber-scheme statement (3).  It should be called
“GZ-type” rather than identified with every classical pairing theorem.

## The verified 24-to-22 quotient

For the repository's rank-compressed 24-variable map, the full constant
kernel has basis

\[
 -3e_6+e_{10},\qquad
 -\frac2{11}e_{12}-\frac3{11}e_{13}+e_{15}.
\]

The quotient has dimension 22, is cubic homogeneous and Keller, and retains
the three distinct rational collision points.  Its constant Jacobian kernel
is zero.  Two complementary checks certify the construction:

- [`verify_constant_kernel_bcw_22_route.py`](../scripts/verify_constant_kernel_bcw_22_route.py)
  constructs the quotient and generated artifact symbolically;
- [`audit_constant_kernel_bcw_22_independent.py`](../scripts/audit_constant_kernel_bcw_22_independent.py)
  independently parses the 24-dimensional source artifact and replays
  `BK=0`, `BC=I`, `H=HCB`, homogeneity, collision descent, and (5) using only
  the Python standard library.

The separate
[`audit_bcw_22_linear_quotients.py`](../scripts/audit_bcw_22_linear_quotients.py)
proves that this particular quotient has no further collision-preserving
linear quotient: its only proper common invariant row module is the
homogenizing covector, which is constant on the collision.

## Mandatory cubicization search protocol

Every completed degree-reduction trace must be scored by the following
pipeline:

1. extract `K=X+Q+C` and certify the identity linear part;
2. rank-factor the cubic output `C=B_c c`;
3. form the rank-compressed cubic-homogeneous map;
4. compute the full constant kernel of its homogeneous Jacobian;
5. quotient that kernel and recheck that the projected collision is separated;
6. compute common invariant-row-module diagnostics on the quotient and run an
   exact characteristic-zero module classification for every shortlisted
   improvement;
7. score the final essential quotient dimension.

The executable experiment
[`search_essential_bcw.py`](../scripts/search_essential_bcw.py) implements
steps 1--5 exactly for every completed beam trace and computes a fast
good-prime coordinate-cyclic row-module diagnostic at step 6.  It emits an
artifact only below the configured certified incumbent dimension.  Any emitted plan
still requires a frozen symbolic generator, an independent replay, and the
full rational invariant-module audit before a theorem statement changes.

The tangent-core pruning rules are orthogonal and should be applied before
this terminal pipeline.  The exact
[`tangent-core audit`](../extended-geometry/TANGENT_CORE_BCW_BYPASS.md) proves
that `r_t` cannot survive as a terminal coordinate, no four residual rows are
admissible, and only the triples `(r_gamma,r_W,r_C)` and
`(r_gamma,r_C,r_s)` retain full row rank at all collision lifts.

## Current frontier

No certified dimension is claimed minimal.  Improvement must presently
come from a better pre-homogenization stable trace, a polynomial circuit/DAG
search, a trace satisfying the two-parameter determinant identity, or a
nonlinear skew-product decomposition.  Dimension 20 would imply `not GMC(40)` by the same
fixed-dimensional bridge.

The bounded regression

```bash
.venv/bin/python scripts/search_essential_bcw.py --width 2 --max-steps 17
```

profiles three completed traces and finds no dimension below 22.  This tests
the pipeline only; it is not evidence of minimality.  The repository makes no
claim to an externally indexed “smallest cubic-homogeneous counterexample”
benchmark.

The faster rank-aware search also has a `kernel-aware` ordering.  It reranks
a bounded prebeam using the modular essential rank of the current
quadratic--cubic truncation.  Widths 24 and 64 produce 44 and 112 terminal
traces; every terminal trace has modular essential rank 21, which safely
excludes an essential dimension below 21 over `Q`.  Exact necessary
two-parameter tests on the 44 width-24 terminals find no pass (the best fail
four of eight samples).  These finite runs are pruning evidence only.

```bash
.venv/bin/python scripts/search_rank_aware_bcw.py --width 64 --max-steps 17 \
  --incumbent 21 --score-mode kernel-aware
.venv/bin/python scripts/search_rank_aware_bcw.py --width 24 --max-steps 17 \
  --incumbent 21 --score-mode kernel-aware --check-two-parameter
```

A larger bounded run,

```bash
.venv/bin/python scripts/search_essential_bcw.py --width 24 --max-steps 17 --incumbent-dimension 22
```

profiles 39 completed traces and emits
[`essential_bcw_candidate.json`](../artifacts/generated-results/essential_bcw_candidate.json)
with the profile

\[
 (s,\operatorname{rank}C,n_{\rm hom},\dim K,n_{\rm ess})
 =(14,6,24,3,21).
\]

The exact symbolic pre-audit
[`verify_essential_bcw_candidate.py`](../scripts/verify_essential_bcw_candidate.py)
replays all 17 stable cancellations, verifies the rank-compressed determinant
bridge, recomputes the three-dimensional constant kernel, checks that the 21D
quotient has zero constant kernel, and verifies three distinct projected
points with one common image.

The plan is now frozen by
[`verify_essential_bcw_21_route.py`](../scripts/verify_essential_bcw_21_route.py),
which writes
[`essential_bcw_21_counterexample.json`](../artifacts/generated-results/essential_bcw_21_counterexample.json).
The dependency-free
[`audit_essential_bcw_21_independent.py`](../scripts/audit_essential_bcw_21_independent.py)
reconstructs the entire route from the original three-dimensional map using
only sparse `Fraction` arithmetic.  It independently recovers rank six, the
24-dimensional homogenization, the three-dimensional constant kernel, the
21-dimensional quotient, its collision, and the triangular determinant
factorization.  The fixed-dimensional bridge therefore gives
`not GMC(42)`.  No minimality is claimed.

The exact
[`audit_bcw_21_linear_quotients.py`](../scripts/audit_bcw_21_linear_quotients.py)
then classifies the common invariant row modules.  The homogenizing covector
is the only proper one and is constant on the stored collision; modulo it the
coefficient algebra is the full `M_20(Q)`, while a common invariant
hyperplane is excluded by the zero common kernel of the 64 nilpotent
coefficient matrices.  Thus no further collision-preserving linear quotient
of this 21-dimensional map exists.

There is also no affine-vector-field extension of the constant-kernel
mechanism.  For an affine infinitesimal translation `V(x)=Ax+b`, cubic
homogeneity separates `JH(x)V(x)=0` into `JH(x)b=0` and
`JH(x)Ax=0`.  The first equation has only `b=0`.  The independent
[`affine symmetry audit`](../scripts/audit_bcw_21_affine_vector_symmetries.py)
expands the second into 2484 rational coefficient equations in the 441
entries of `A`; their matrix has full column rank modulo the good prime
1000003, hence full column rank over `Q`.  Therefore `A=0` as well.  This
rules out affine translation symmetries as seeds for a nonlinear
skew-product reduction, but it does not exclude general nonlinear
decompositions.
