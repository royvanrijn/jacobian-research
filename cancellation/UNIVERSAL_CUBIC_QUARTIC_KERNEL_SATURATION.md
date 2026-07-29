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

## 5. Canonical-different complex and the Fitting reduction

Write the trace-free and scalar multiplication components as
\(\mu_{ij}\in M\) and \(s_{ij}\in R\).  Over the full 24-parameter ring,
form the seven-column canonical-different matrix

\[
 d_1=\left[
 (0,z,-y,x)^{\mathsf T},
 (s_{ij},2\mu_{ij})^{\mathsf T}_{0\leq i\leq j\leq2}
 \right].                                             \tag{5.1}
\]

The exact checker constructs a universal \(7\)-by-\(3\) matrix \(d_2\).
For \(r=(z,-y,x)\), associativity and the coefficient-module relation give

\[
 \sum_i r_i s_{ij}=0,\qquad
 2\sum_i r_i\mu_{ij}=q_jr.                            \tag{5.2}
\]

The top row of \(d_2\) is \((-q_0,-q_1,-q_2)\); its other six rows are the
fixed incidence coefficients expressing
\(\sum_i r_i(s_{ij},2\mu_{ij})\).  Thus \(d_1d_2=0\) identically.  Each
\(q_j\) has only parameter/collision bidegrees \((0,2)\) and \((1,3)\).
The fixed lower \(6\)-by-\(3\) block has maximal minors containing
\(x^3,y^3,z^3\), while \(d_1\) has a nonzero central \(4\)-by-\(4\) minor
for every squarefree symbol.  The Buchsbaum--Eisenbud grade criterion
therefore proves exactness of

\[
 0\longrightarrow R^3\mathop{\longrightarrow}^{d_2}R^7
 \mathop{\longrightarrow}^{d_1}R^4.                  \tag{5.3}
\]

Let \(T_\Delta=\operatorname{coker}(d_1)\).  Transposing the last
differential computes \(\operatorname{Ext}^2_R(T_\Delta,R)\).  The six
fixed linear rows alone generate a module \(L\subset R^3\) satisfying

\[
 \dim_{\mathbb Q}(R^3/L)=6,\qquad
 (x,y,z)^2(R^3/L)=0.                                  \tag{5.4}
\]

The varying top row lies in \((x,y,z)^2R^3\), so it is redundant.  Hence

\[
 \operatorname{Ext}^2_R(T_\Delta,R)
 \simeq (R^3/L)
 \simeq E_0\otimes_{\mathbb Q}
 \mathbb Q[u_1,\ldots,u_{24}].                        \tag{5.5}
\]

After truncation by \((x,y,z)^2\), the parameter module has 12 generators
and six independent constant relations.  It is free of rank six, proving
universally for the canonical-different support

\[
 \operatorname{Fitt}_6=(1),\qquad
 \operatorname{Fitt}_5=(0).                           \tag{5.6}
\]

The remaining issue is now an identification, not a Fitting-minor
calculation.  The ramification support in (1.2) uses
\[
 T=B/\operatorname{Ann}_B(\Omega_{B/R}),
\]
whereas (5.1) defines \(T_\Delta\) from the seven canonical different
generators.  On the full-support plane of Section 2, exact module reduction
proves that these seven generators span the complete annihilator for every
parameter and every squarefree symbol.  The resulting actual minimal
support resolution has tail

\[
 R^3\mathop{\longrightarrow}^{d_2}R^7
 \longrightarrow R^4\longrightarrow T\longrightarrow0.       \tag{5.7}
\]

After the exact minimal-resolution basis chosen by Singular, rows two
through seven of \(d_2\) are parameter-independent linear triples.  Let
\(L\subset R^3\) be the module they generate.  The checker verifies

\[
 \dim_{\mathbb Q}(R^3/L)=6,\qquad
 (x,y,z)^2(R^3/L)=0.                                  \tag{5.8}
\]

The remaining row lies in \((x,y,z)^2R^3\).  Its central part is quadratic,
its parameter-dependent part lies in \((x,y,z)^3R^3\), and that part is
linear in the two plane parameters.  Equation (5.2) therefore makes the
entire seventh row redundant:

\[
 \operatorname{coker}(d_2^{\mathsf T})
 \simeq R^3/L.                                        \tag{5.9}
\]

Thus (5.6) already closes the universal Fittings for \(T_\Delta\), and
(5.9) checks their identification with the requested Fittings on seven
full-support planes.  To close them globally for the actual \(T\), it is
enough to prove the universal annihilator--different equality

\[
 \operatorname{Ann}_B(\Omega_{B/R})
 =
 \left((0,z,-y,x),(s_{ij},2\mu_{ij})\right).          \tag{5.10}
\]

Equivalently, one must exclude additional annihilator generators supported
over a proper parameter locus.  No further determinant calculation is
needed once (5.10) is established.

There is now an exact conditional closure of this last step.  The universal
Deligne--Faddeev cubic algebra on a free trace-free rank-two module has

\[
 \operatorname{Fitt}_0^B(\Omega_{B/R})
 =
 \operatorname{Ann}_B(\Omega_{B/R}).                  \tag{5.11}
\]

The checker computes both ideals in the universal algebra over
\(\mathbb Q[a,b,c,d]\) and reduces them to the same three-generator ideal.
On \(D(x)\cup D(y)\cup D(z)\), the Koszul trace-free module is locally free,
so (5.11) identifies the canonical different in (5.1) with the actual
annihilator.

The exact complex (5.3) also shows that \(T_\Delta\) has projective
dimension at most two.  If a prime contains \((x,y,z)\), its height is at
least three; Auslander--Buchsbaum therefore gives depth at least one for
\((T_\Delta)_\mathfrak p\).  Hence

\[
 H^0_{(x,y,z)}(T_\Delta)=0.                           \tag{5.12}
\]

Assume now the requested universal relative cotangent saturation
\(H^0_{(x,y,z)}(\Omega_{B/R})=0\).  The canonical different annihilates
\(\Omega\) off the collision axis, so its action on \(\Omega\) is
\((x,y,z)\)-torsion and therefore vanishes globally.  Thus
\(\Delta\subseteq\operatorname{Ann}(\Omega)\).  Their quotient is supported
on the collision axis by (5.11) and injects into \(T_\Delta\); (5.12)
forces the quotient to vanish.  Consequently

\[
 H^0_{(x,y,z)}(\Omega)=0
 \quad\Longrightarrow\quad
 T=T_\Delta,\quad
 \operatorname{Fitt}_6(E)=(1),\quad
 \operatorname{Fitt}_5(E)=0.                         \tag{5.13}
\]

Thus there is no independent Ext-Fitting exceptional set: it is contained
in the cotangent-saturation failure locus.  Closing universal cotangent
saturation closes the requested Fittings at the same time.

The checker
[`verify_cubic_quartic_ext_tail_absorption.py`](../scripts/verify_cubic_quartic_ext_tail_absorption.py)
proves the actual-support statement on the seven full-support planes.
The checker
[`verify_universal_cubic_quartic_different_complex.py`](../scripts/verify_universal_cubic_quartic_different_complex.py)
proves (5.2)--(5.6) over all 24 parameters.  Neither checker proves the
universal equality (5.10) unconditionally.  The checker
[`verify_universal_cubic_kahler_annihilator.py`](../scripts/verify_universal_cubic_kahler_annihilator.py)
proves (5.11), which gives the conditional implication (5.13).

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

The canonical-different construction of Section 5 completes the
parameter-only matrix that this bottleneck originally suggested.  Its six
linear relation columns act on the nine-dimensional space

\[
 (x,y,z)\,S^3/(x,y,z)^2S^3.                         \tag{6.1}
\]

The Fittings of this matrix are now closed by (5.6).  The next certificate
is instead the equality (5.10).  A useful formulation is to compute the
finite quotient

\[
 \operatorname{Ann}_B(\Omega_{B/R})/\Delta,
 \qquad
 \Delta=((0,z,-y,x),(s_{ij},2\mu_{ij})),             \tag{6.2}
\]

and prove that it is zero.  It vanishes on every recorded line, plane, and
the smooth coordinate ten-space.  By (5.13), however, this quotient
vanishes automatically once the cotangent-saturation quotient does.
Therefore the remaining universal calculation is the single relative
cotangent saturation test; the Fittings no longer require a separate
elimination.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/research_universal_cubic_quartic_kernel_saturation.py
.venv/bin/python scripts/verify_cubic_quartic_ext_tail_absorption.py
.venv/bin/python scripts/verify_universal_cubic_quartic_different_complex.py
.venv/bin/python scripts/verify_universal_cubic_kahler_annihilator.py
```

The generated record is
[`universal_cubic_quartic_kernel_saturation_frontier.json`](../artifacts/generated-results/universal_cubic_quartic_kernel_saturation_frontier.json).

The calculations require Singular 4.4.1.  The canonical-different complex
and its Fitting ideals are computed over all 24 parameters.  The locally
free universal cubic calculation and the depth argument prove that
cotangent saturation implies equality with the actual annihilator and
hence the same Fittings.  The unconditional universal cotangent saturation
and its parameter discriminant remain open.
