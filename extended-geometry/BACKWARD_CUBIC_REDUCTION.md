# Backward cubic reduction

## Result and scope

The repository now treats the reverse of rank-compressed homogenization as a
first-class reduction.  If

\[
 F(x)=x+Q(x)+Bc(x)                                             \tag{1}
\]

has quadratic \(Q\), cubic vector \(c\), and constant matrix \(B\), its
rank-compressed cubic-homogeneous parent is

\[
 V(x,y,t)=
 \bigl(x+tQ(x)+t^2By,\;y-c(x),\;t\bigr).                       \tag{2}
\]

The backward route has two exact steps:

1. restrict (2) to the fixed level \(t=1\);
2. cancel all companion variables \(y\) by polynomial left--right
   equivalence.

On A. MacFarlane's external certificate this gives the fully checked chain

\[
 \boxed{G_{20}\ \longrightarrow\ M_{19}\ \sim\ F_{13}\times I_6}.      \tag{3}
\]

This is a calibration theorem and a search architecture.  It does **not**
construct a degree-three Keller counterexample in dimension twelve.

The implementation is
[`backward_cubic.py`](../jcsearch/backward_cubic.py), checked by
[`verify_backward_cubic_reduction.py`](../scripts/verify_backward_cubic_reduction.py).
The status-level checker is the complete
[`verify_backward_cubic_suite.py`](../scripts/verify_backward_cubic_suite.py),
which also runs the map-specific audits and current-field regressions.
The MacFarlane formulas and direct linear obstruction audit remain in
[`MACFARLANE_G20_DIMENSION_REDUCTION_AUDIT.md`](MACFARLANE_G20_DIMENSION_REDUCTION_AUDIT.md).

## 1. Collision-compatible fixed-covector restriction

Let

\[
 \Phi=I+H:\mathbb A^N\longrightarrow\mathbb A^N,\qquad
 \det D\Phi=1,
\]

and let \(\lambda\ne0\) be a constant covector satisfying

\[
 \lambda H=0.                                                   \tag{4}
\]

Then \(\lambda\Phi=\lambda\), so every affine hyperplane
\(\lambda x=a\) is invariant.  Choose linear coordinates whose last
coordinate is \(\lambda x\).  In those coordinates

\[
 \Phi(u,s)=(\phi_s(u),s),\qquad
 D\Phi=
 \begin{pmatrix}
 D\phi_s&*\\
 0&1
 \end{pmatrix}.                                                 \tag{5}
\]

Consequently

\[
 \det D\phi_s=1.                                                \tag{6}
\]

If two points \(p\ne q\) in one fiber of \(\Phi\) satisfy
\(\lambda p=\lambda q=a\), their restrictions give a collision for
\(\phi_a\).  More generally, if the fixed-covector space has basis
\(\lambda_1,\ldots,\lambda_r\), the collision-compatible part is the kernel
of

\[
 (\lambda_1,\ldots,\lambda_r)^T(p_i-p_0)
\]

over the selected collision points.  This is exact linear algebra.

The important bookkeeping correction is that only one pair must survive.
For a collision set \(p_0,\ldots,p_m\), every pair is tracked separately.
A projection or restriction is admissible when at least one pair remains
distinct and has one common image.  Requiring all three foundational points
to remain distinct is stronger than noninjectivity and can discard valid
reductions.

## 2. Companion cancellation theorem

Specializing (2) at \(t=1\) gives

\[
 M(x,y)=\bigl(x+Q(x)+By,\;y-c(x)\bigr).                         \tag{7}
\]

Define triangular polynomial automorphisms

\[
 S_c(x,y)=(x,y-c(x)),\qquad
 A_B(u,v)=(u+Bv,v).                                             \tag{8}
\]

Then direct substitution gives

\[
 \boxed{
 M=A_B\circ(F\times I_{\operatorname{rank}c})\circ S_c.
 }                                                               \tag{9}
\]

Thus:

1. \(M\) and \(F\times I\) are polynomially left--right equivalent;
2. \(\det DM=\det DF\);
3. a lifted collision \((p,y_p),(q,y_q)\) with
   \(y_p=c(p)\), \(y_q=c(q)\) descends to the collision \(F(p)=F(q)\);
4. all companion coordinates have a certified lifetime ending at \(t=1\).

The repository already used (9) implicitly inside the forward determinant
bridge.  The new point is to emit \(F\) as a terminal object rather than
treating it only as input to (2).

## 3. The whole homogenizing line

The cancellation is not confined to \(t=1\).  Put

\[
 E_t(x)=x+tQ(x)+t^2Bc(x).
\]

Over the full affine \(t\)-line, define the target shear

\[
 A_{t^2B}(u,v,t)=(u+t^2Bv,v,t).
\]

The same componentwise calculation gives the relative identity

\[
 \boxed{
 V=A_{t^2B}\circ(E_t\times I)\circ S_c.
 }                                                               \tag{10}
\]

Because \(Q\) and \(c\) are homogeneous of degrees two and three,

\[
 E_t(x)=t^{-1}F(tx) \qquad (t\ne0).                              \tag{11}
\]

Thus the family is isotrivial over \(\mathbb G_m\): every nonzero fiber is
stably left--right equivalent to a linearly scaled copy of \(F\).  At the
special fiber,

\[
 V|_{t=0}(x,y)=(x,y-c(x)),                                      \tag{12}
\]

which is a triangular polynomial automorphism with inverse
\((x,v)\mapsto(x,v+c(x))\).

This classifies every possible collision of the rank-compressed parent.
Since the last output is \(t\), colliding points have one common parameter
value.  Equation (12) excludes value zero.  Equations (10)--(11) then scale
the collision to \(t=1\).  Consequently:

> Searching separate nonzero homogenizing levels is redundant.  Every
> parent collision can be normalized to \(t=1\), and the \(t=0\) branch can
> be discarded by theorem rather than bounded testing.

The checker transports the published MacFarlane collision to a second exact
rational slice \(t=2\) and verifies its common image.  This is a regression
of the fiber classification, not a second counterexample.

There is a scheme-level consequence.  Over
\(\mathbb Q[t,t^{-1}]\), apply the source shear, target shear, and scaling to
both points in the fiber product \(V\times_VV\).  The collision equations
become the collision equations of \(F\), equality of the residual companion
coordinates, and one free nonzero parameter.  Thus, up to the displayed
stable coordinate changes,

\[
 \operatorname{Coll}(V)|_{\mathbb G_m}
 \simeq
 \operatorname{Coll}(F)\times\mathbb G_m\times\mathbb A^r.     \tag{13}
\]

The same statement holds for the off-diagonal collision scheme.  Rank
compression therefore transports collision geometry; it does not create a
new nonzero-parameter collision component.  A genuinely smaller construction
must change the upstream circuit or use a transformation that mixes the
homogenizing direction before this product structure becomes terminal.

For a base collision point \(p\), the canonical point on the \(t\)-fiber is

\[
 \left(t^{-1}p,\ c(t^{-1}p),\ t\right)
 =
 \left(t^{-1}p,\ t^{-3}c(p),\ t\right).                         \tag{14}
\]

Hence the collision escapes every affine compact set as \(t\to0\), explaining
why the triangular special fiber contains no collision.  The weighted
coordinates

\[
 X=tx,\qquad Y=t^3y
\]

have a finite boundary limit \((p,c(p))\).  Rank-compressed homogenization is
therefore a Rees-like degeneration whose weighted exceptional data remembers
the dehomogenized collision.  Collision-ideal calculations should saturate
by \(t\) before specializing, or use the weight ledger

\[
 \operatorname{wt}(t,x,y)=(1,-1,-3).
\]

Ordinary specialization at \(t=0\) erases the relevant component; weighted
boundary extraction recovers it.  For a degree-\(d\) companion block, the
same argument predicts companion pole order \(d\).

The relative cancellation identity itself is degree-agnostic.  The scaling
formula (11) uses precisely the quadratic/cubic weights, suggesting the
weighted analogue for higher-degree reduction systems.

## 4. Two objective functions, not one

For a terminal \(F=x+Q+C\) in dimension \(n\), let

\[
 r=\dim\operatorname{span}\{C_1,\ldots,C_n\}.
\]

There are two different minimization problems:

\[
\begin{aligned}
 \text{arbitrary degree-three:}\quad &(n,\ n+r+1,\ r),\\
 \text{cubic-homogeneous:}\quad &(n+r+1,\ n,\ r).
\end{aligned}                                                    \tag{15}
\]

The first coordinate in each row is the primary objective.  A state with
\((n,r)=(12,10)\) beats \((13,5)\) for the degree-three problem, even though
it is much worse for the homogeneous problem.  A beam ordered only by
\(n+r\) can therefore delete the desired dimension-twelve state.

The implementation exposes both
`BackwardTerminalProfile.direct_cubic_key` and
`BackwardTerminalProfile.homogeneous_key`.  Future searches must retain
separate Pareto archives for the two classes.

Nilpotency and the generic rank of the homogeneous correction remain useful
secondary data for the homogeneous archive.  They are not admissibility
conditions for the direct degree-three archive: MacFarlane's \(F_{13}\) is
Keller although its nonlinear Jacobian is not nilpotent.

## 5. Exact MacFarlane calibration

For the formulas reproduced in the companion audit,

\[
 F_{13}=x+R+B\gamma,\qquad \operatorname{rank}B=6,
\]

and

\[
 G_{20}(x,w,\tau)
 =
 \bigl(x+\tau R+\tau^2Bw,\;w-\gamma,\;\tau\bigr).
\]

The checker proves:

1. the sole fixed linear covector of \(G_{20}\) is \(\tau\);
2. both collision points lie on \(\tau=1\);
3. restriction gives
   \[
   M_{19}(x,w)=(x+R+Bw,\;w-\gamma);
   \]
4. (9) holds in all nineteen components with \(c=\gamma\);
5. the source shear sends both lifted collision points to companion value
   zero;
6. the remaining base collision is exactly the published collision of
   \(F_{13}\);
7. the two terminal keys are
   \[
   (13,20,6)\quad\hbox{and}\quad(20,13,6).
   \]

The generated certificate is
[`backward_cubic_reduction_calibration.json`](../artifacts/generated-results/backward_cubic_reduction_calibration.json).

This regression is mandatory for a future compiler.  A search representation
which cannot reconstruct (3) is not expressive enough to search for
dimension twelve.

## 6. Compiler state and transitions

A backward-search state must contain:

1. the placed polynomial circuit and exact current coordinate expressions;
2. source/target ownership for every exposed gate;
3. the lifted points for every still-surviving collision pair;
4. the full left kernel of the nonlinear component vector;
5. the full constant right kernel of its Jacobian;
6. every rank factorization \(C=Bc\);
7. a lifetime flag for companion variables introduced by \(c\);
8. both objective keys in (15).

The admissible transitions are:

1. determinant-preserving BCW source and target shears;
2. grouped product cancellations with exact ownership;
3. restriction to a collision-compatible fixed-covector level;
4. companion cancellation after fixing the homogenizing coordinate;
5. constant-kernel quotient, with pairwise collision separation rechecked;
6. exact linear cleanup and removal of literal identity factors.

Every transition must carry an exact factorization identity.  Modular rank
may order or prune candidates, but it cannot certify a terminal.

## 7. Application to the active restricted-minima field

The active circuit search in
[`search_restricted_bcw_circuits.py`](../scripts/search_restricted_bcw_circuits.py)
now consumes the backward interface at every terminal.  It:

1. computes the exact direct and raw-homogeneous keys in (15);
2. records all surviving collision pairs before and after the iterated
   constant-kernel quotient;
3. rejects a terminal only when every collision pair collapses;
4. emits separate best direct and raw-homogeneous terminal archives among
   all states reached by the bounded beam;
5. retains the existing rank/index/Hessian Pareto objective as a separate
   restricted-minima question.

The forward helper
[`rank_compressed_bcw_homogenization.py`](../scripts/rank_compressed_bcw_homogenization.py)
also invokes companion cancellation inside every parametric-factorization
check.  Consequently the established \(16\to24\) BCW regression now checks
both directions:

\[
 K_{16}\longrightarrow V_{24}
 \quad\hbox{and}\quad
 V_{24}|_{t=1}\sim K_{16}\times I_7.
\]

The exact current-archive audit examines the thirteen Pareto representatives
retained across ten active restricted-minima searches.  Their best direct
key is

\[
 (18,26,7),
\]

their best raw-homogeneous key is \((26,18,7)\), and their smallest
post-kernel homogeneous dimension is \(22\).  Two representatives—the best
direct source and the best final quotient—are reconstructed from their
stored circuit plans, and all three collision pairs survive exactly.

This is useful negative information about the present archive, but it is not
a dimension lower bound.  Those representatives were selected by
rank/index/Hessian objectives, so a direct-dimension winner could have been
discarded before serialization.  The result proves that a separate
direct-dimension archive is necessary; it does not prove that the old search
space contains no twelve-dimensional state.

There is also an immediate certified-frontier consequence.  With the pinned
external determinant certificate, MacFarlane \(G_{20}\) improves the ambient
cubic-homogeneous upper bound to

\[
 n_{\rm cub}\le20.
\]

Its standard homogeneous cotangent lift is a quartic HN counterexample in
forty variables, so

\[
 n_{\rm HN,4}\le40.
\]

These external-certificate endpoints improve the repository ledger from
\(21\) and \(42\).  The old witnesses remain the smallest internally
generated, dependency-free independent replays.

The generated record is
[`backward_cubic_current_applications.json`](../artifacts/generated-results/backward_cubic_current_applications.json).

## 8. Certification layers

The theorem is checked at four levels:

1. a symbolic non-coordinate fixed-covector example verifies both ambient
   and restricted Jacobian determinant one;
2. a generic quadratic--cubic companion fixture verifies
   \[
   \det D(x+Q+By,y-c)=\det D(x+Q+Bc)
   \]
   as a polynomial identity;
3. the pinned MacFarlane formulas verify the complete
   \(20\to19\sim13+6\) collision chain;
4. current BCW and restricted-minima terminal pipelines replay the reverse
   transition on repository witnesses.

The suite pins the SHA-256 of the reusable module, both external-map audits,
both backward verifiers, and the two integrated search/forward helpers before
executing them.  A changed component therefore invalidates the status checker
until the change and its new pin are reviewed together.

These are exact rational symbolic checks.  They are not a Lean
formalization or an external independent review.

## 9. Resolution of `BCR-OPEN1`

> Construct a determinant-preserving upstream schedule from the foundational
> three-variable collision to a square map
> \[
> K=I+Q+C:\mathbb A^n\to\mathbb A^n
> \]
> with \(n\le12\), degree at most three, and at least one exact surviving
> collision pair.

This target is now met by the
[twelve-variable coordinate-pair reduction](../verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md).
The
[direct MacFarlane audit](MACFARLANE_G20_DIMENSION_REDUCTION_AUDIT.md)
already closes constant input kernels, fixed linear output coordinates,
Keller affine-hyperplane restrictions through the stored collision, and
polynomial invariant coordinates through degree three for the displayed
\(F_{13}\).  The missing operation was not another pullback-fixed invariant.
It was a **source-target coordinate pair**:

\[
s=F_{13,13}=x_{13}+x_2^2,\qquad
y_4\longmapsto y_4-y_8^2.
\]

After the source change and target square completion, the transformed map is

\[
\bigl(K_1,K_2,K_3,K_4+s(2z_{12}-z_1^2),K_5,\ldots,K_{12},s\bigr).
\]

Restriction to \(s=0\) gives the exact twelve-variable map.  Its cubic-output
rank remains six, so rank-compressed homogenization gives dimension nineteen.
The successful compiler lesson is that the state must track candidate source
coordinates \(h\), target coordinates \(g\), and identities
\(g\circ F=h\), not only fixed invariants \(P\circ F=P\).  Degree cleanup may
also require a target completion before the graph coordinate is deleted.

## Reproduction

```bash
make verify-backward-cubic-reduction
make verify-macfarlane-f12
```

The first target verifies the generic backward operations and the complete
\(20\to19\sim13+6\) calibration.  The second verifies the
\(13\to12\to19\) theorem and its independent replay.
