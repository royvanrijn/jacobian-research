# Universal cubic quartic-kernel saturation frontier

## Status

The full 24-parameter question is open.  No exceptional parameter has been
found, but the calculations below do not prove that the exceptional locus is
empty.

This note separates three levels which must not be conflated:

1. exact polynomial-family calculations on specified parameter subspaces;
2. exact calculations at isolated dense parameter points;
3. the unresolved calculation over
   \(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\).

Only the first level gives a flatness statement for every point of a tested
subspace.

## 1. Universal family

Let

\[
 A=\mathbb Q[x,y,z],\qquad
 M=\operatorname{coker}\left(
 A\xrightarrow{(z,-y,x)^{\mathsf T}}A^3
 \right).
\]

For a squarefree ternary-cubic representative \(h\), write \(\Phi_h\) for
the homogeneous generalized triple-cover tensor on \(M\).  The order-four
compatibility equations have the fixed primitive integral basis
\(\psi_1,\ldots,\psi_{24}\) constructed by
[`verify_cubic_symbol_double_saturation.py`](../scripts/verify_cubic_symbol_double_saturation.py).
Put

\[
 S=\mathbb Q[u_1,\ldots,u_{24}],\qquad
 R=S[x,y,z],
\]

and form the exact universal tensor

\[
 \Phi_{\mathrm{univ}}
 =\Phi_h+\sum_{i=1}^{24}u_i\psi_i.                 \tag{1.1}
\]

The checker reconstructs from (1.1) the associative multiplication on
\(R\oplus(M\otimes_A R)\), the 31-column presentation \(N\) of
\(\Omega_{B/R}\), the ramification-support module

\[
 T=B/\operatorname{Ann}_B(\Omega_{B/R}),
\]

and

\[
 E=\operatorname{Ext}^2_R(T,R).                    \tag{1.2}
\]

The desired universal theorem is:

\[
\begin{aligned}
&(N:(x,y,z)^\infty)/N=0,\\
&\sqrt{\operatorname{Fitt}_0^R(E)}=(x,y,z),\\
&E\text{ is }S\text{-flat of rank }6.
\end{aligned}                                      \tag{1.3}
\]

For the finite \(S\)-module obtained from the verified
\((x,y,z)^2E=0\) truncation, flat rank six is equivalent to

\[
 \operatorname{Fitt}_6^S(E)=(1),\qquad
 \operatorname{Fitt}_5^S(E)=(0).                   \tag{1.4}
\]

Thus the parameter exceptional locus is scheme-theoretically contained in
the union of the nonunit locus of \(\operatorname{Fitt}_6^S(E)\), the
support of \(\operatorname{Fitt}_5^S(E)\), the parameter projection of the
cotangent-saturation quotient, and the locus where the radical in (1.3)
changes.  Computing those ideals over all 24 parameters is the unresolved
discriminant calculation.

## 2. Exact full-support planes and lines

The strongest current mixed-support calculation is KDFP6.  With
\(\psi_+=\sum_i\psi_i\) and
\(\psi_-=\sum_i(-1)^i\psi_i\), it treats

\[
 \Phi_h+u\psi_++v\psi_-
\]

over \(\mathbb Q[u,v,x,y,z]\) for every squarefree symbol.  The cotangent
presentation is saturated and the pruned rank-three Ext presentation is
pulled back from the origin with multiplicity six.  On
\(D(u^2-v^2)\), all 24 fixed basis coordinates are nonzero.  Thus basis
sparsity is not a necessary condition for retaining the defect, although
one parameter plane is not an open subset of \(\mathbb A^{24}\).

Four fixed coefficient vectors with all 24 entries nonzero are recorded in
[`research_universal_cubic_quartic_kernel_saturation.py`](../scripts/research_universal_cubic_quartic_kernel_saturation.py).
For each vector \(a=(a_i)\) and each of the seven squarefree cubic-symbol
orbits, the checker works over \(\mathbb Q[t,x,y,z]\) with

\[
 \Phi_h+t\sum_i a_i\psi_i.                          \tag{2.1}
\]

All 28 polynomial families satisfy, exactly:

\[
\begin{array}{c|c}
\text{invariant}&\text{value}\\ \hline
\text{relative cotangent saturation quotient}&0\\
\operatorname{mult}(E)&6\\
t\text{-torsion of }E&0\\
\sqrt{\operatorname{Fitt}_0(E)}&(x,y,z)\\
\text{difference from the central Ext presentation}&0.
\end{array}                                         \tag{2.2}
\]

Consequently (1.4) holds after restriction to each recorded line, and the
universal exceptional locus has empty intersection with all 28 lines.
This is stronger than checking a dense endpoint: it excludes every scalar
on each line.  A finite set of lines is not Zariski dense enough to settle
the 24-parameter question.

## 3. Exact higher-dimensional subspace

For the smooth cubic symbol, the same checker treats the full coordinate
subspace on the first ten basis tensors:

\[
 \Phi_h+\sum_{i=1}^{10}p_i\psi_i
 \quad\text{over}\quad
 \mathbb Q[p_1,\ldots,p_{10},x,y,z].                \tag{3.1}
\]

The relative cotangent presentation is saturated, the radical support of
\(E\) is exactly the parameter ten-space at \(x=y=z=0\), its multiplicity
is six, and its pruned rank-three presentation is pulled back from the
parameter origin.  Hence the parameter discriminant restricts to the empty
scheme on this \(\mathbb A^{10}\).

This extends the previously recorded smooth coordinate-three-space result
for one specified nested subspace.  It does not cover all coordinate
ten-spaces or an arbitrary ten-dimensional linear subspace.

## 4. Universal cotangent input reduction

There is now an exact calculation over the full ring
\(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\), but it is an input reduction rather
than the desired saturation theorem.  Comparing the universal 31-column
cotangent presentation \(N_{\mathrm{univ}}\) with its central value \(N_0\)
entry by entry gives only the parameter/collision bidegrees

\[
 (1,3),\qquad (1,5),\qquad (2,6).                 \tag{4.1}
\]

More precisely, 147 matrix entries change.  Their expanded terms comprise
1,019 terms of bidegree \((1,3)\), 933 of bidegree \((1,5)\), and 4,800 of
bidegree \((2,6)\).  In particular,

\[
 N_{\mathrm{univ}}\equiv N_0\pmod{(x,y,z)^3}.      \tag{4.2}
\]

The collision-degree bound follows structurally from the construction:
the quartic tensor changes the trace-free multiplication in degree three,
its cross term with the cubic multiplication enters the scalar part in
degree five, and its square enters in degree six.  Thus the same bound
applies to every cubic-symbol row, although the serialized matrix hash in
the artifact uses the smooth representative.

Before any standard-basis computation, six parameter-independent unit
pivots split off from the \(12\)-by-\(31\) presentation.  Exact elementary
row and column operations therefore replace it by a cokernel-equivalent
\(6\)-by-\(25\) presentation.  The six successive pivot values are
\(2,1,1,2,1,2\), hence are units over \(\mathbb Q\) on all of parameter
space.  The generated artifact records the pivot positions and a SHA-256
hash of the reduced universal matrix.

Equations (4.1)--(4.2) explain why the verified two-layer Ext module is
insensitive on many parameter slices: no quartic parameter occurs in the
two-jet of the cotangent presentation.  They do **not** prove saturation or
flatness.  Taking the annihilator of the cotangent module and then resolving
its support can convert higher-order relations into new low-order
syzygies.  Any exceptional locus must therefore enter through this
syzygy-lifting step, rather than through the raw two-jet.

## 5. Resolution-tail absorption and the Fitting reduction

There is a sharper finite reduction on the full-support plane of Section 2.
For every squarefree symbol, the minimal support resolution has tail

\[
 R^3\mathop{\longrightarrow}^{d_2}R^7
 \longrightarrow R^4\longrightarrow T\longrightarrow0.       \tag{5.1}
\]

After the exact minimal-resolution basis chosen by Singular, rows two
through seven of \(d_2\) are parameter-independent linear triples.  Let
\(L\subset R^3\) be the module they generate.  The checker verifies

\[
 \dim_{\mathbb Q}(R^3/L)=6,\qquad
 (x,y,z)^2(R^3/L)=0.                                  \tag{5.2}
\]

The remaining row lies in \((x,y,z)^2R^3\).  Its central part is quadratic,
its parameter-dependent part lies in \((x,y,z)^3R^3\), and that part is
linear in the two plane parameters.  Equation (5.2) therefore makes the
entire seventh row redundant:

\[
 \operatorname{coker}(d_2^{\mathsf T})
 \simeq R^3/L.                                        \tag{5.3}
\]

This proves the Fitting assertion on all seven planes without computing
large minors.  After truncation by \((x,y,z)^2\), the parameter module has
12 generators and six independent constant relations.  It is free of rank
six, so

\[
 \operatorname{Fitt}_6=(1),\qquad
 \operatorname{Fitt}_5=(0).                           \tag{5.4}
\]

More importantly, (5.2) identifies the exact missing universal lemma.  It
is enough to prove that the universal support resolution is
**Rees-strict at its last differential**: its linear strand must be the
same six-row strand, and every additional row must lie in
\((x,y,z)^2R^3\).  If so, all additional rows are automatically absorbed
by (5.2), (5.3) holds over
\(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\), and the universal Fitting
discriminant is empty.

The checker
[`verify_cubic_quartic_ext_tail_absorption.py`](../scripts/verify_cubic_quartic_ext_tail_absorption.py)
proves this tail statement on the seven full-support planes.  It does not
yet prove universal Rees strictness; a higher-order support syzygy could in
principle change the minimal tail away from all tested planes.

## 6. Universal elimination bottleneck

The complete smooth universal input has 27 variables and approximately
0.5 MB of exact Singular source.  Two direct routes were tested:

- a polynomial ring in all 27 variables;
- the rational function field
  \(\mathbb Q(u_1,\ldots,u_{24})[x,y,z]\).

Both bottleneck before saturation or Ext, while standardizing the initial
31-relation cotangent module.  The polynomial-ring route did not reach its
first invariant after ten minutes.  The rational-function-field route
exceeded 4 GB resident memory before the same point.  A block order with
\((x,y,z)\) first reduced memory but did not remove the elimination
bottleneck.  On the smooth nested coordinate family, ten parameters
complete.  The corresponding first-eleven-parameter calculation did not
complete within a fresh 300-second bound; an earlier first-twelve-parameter
calculation did not complete within 900 seconds and reached approximately
1.27 GB resident memory.  These timings are failed computational routes,
not evidence for or against flatness.

Unit pruning reduces the universal input to the \(6\)-by-\(25\) matrix of
Section 4 and approximately 0.2 MB of Singular source.  With
\((x,y,z)\) placed in the first elimination block, this reduced calculation
still did not return its first standard-basis diagnostic after six minutes
and reached approximately 1.75 GB resident memory.  This is another failed
direct route, not an exceptional parameter.

The useful homological reduction is the known two-layer condition
\((x,y,z)^2E=0\).  A universal proof should avoid the full module standard
basis and construct the finite parameter presentation directly from the
minimal resolution of \(T\).  Its six linear relation columns act on the
nine-dimensional space

\[
 (x,y,z)\,S^3/(x,y,z)^2S^3.                         \tag{6.1}
\]

The next certificate is therefore a parameter-only matrix:

1. construct the universal minimal second differential without localizing
   at a parameter pivot;
2. prove universally that the quadratic action vanishes;
3. compute (1.4) from the resulting finite \(S\)-presentation;
4. add the independent cotangent-saturation certificate.

This route computes the requested discriminant directly.  Pivot charts are
acceptable only if their exceptional ideals are retained and the charts
are proved to cover \(\mathbb A^{24}\).

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/research_universal_cubic_quartic_kernel_saturation.py
.venv/bin/python scripts/verify_cubic_quartic_ext_tail_absorption.py
```

The generated record is
[`universal_cubic_quartic_kernel_saturation_frontier.json`](../artifacts/generated-results/universal_cubic_quartic_kernel_saturation_frontier.json).

The calculation requires Singular 4.4.1.  It is an exact finite-subspace
computation.  The arbitrary 24-parameter combination, its global Fitting
ideals, and its parameter discriminant remain open.
