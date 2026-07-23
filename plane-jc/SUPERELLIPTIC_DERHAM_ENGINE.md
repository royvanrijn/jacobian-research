# A superelliptic de Rham engine for plane-JC leading blocks

## Status

The exact reduction kernel described here is implemented in
[`cas/superelliptic_derham.py`](cas/superelliptic_derham.py), with regression
tests in [`cas/test_superelliptic_derham.py`](cas/test_superelliptic_derham.py).
It reproduces the structural facts of the audited `(72,108)` first block and
works for non-hyperelliptic superelliptic examples.

The compiler IR, explicit `(72,108)` tail-obstruction certificate, and
source-aware `(96,144)` record are in
[`NEWTON_DERHAM_COMPILER.md`](NEWTON_DERHAM_COMPILER.md) and
[`cas/newton_derham_compiler.py`](cas/newton_derham_compiler.py).  The engine
also computes character-wise Gauss--Manin connection matrices.

Published Laurent polygons can now be compiled mechanically to the leading
block by [`cas/laurent_band_frontend.py`](cas/laurent_band_frontend.py); see
[`LAURENT_BAND_FRONTEND.md`](LAURENT_BAND_FRONTEND.md).  The principal missing
interface has narrowed to deriving an exhaustive normalized polygon list from
an admissible corner chain.

This is an engine for a leading coefficient block once its Newton bands and
supports have been derived.  It does **not** derive Laurent polygons from an
admissible corner chain.  That is the principal missing interface for the
frontier pairs `(96,144)` and `(75,125)`.

## 1. The differential identity

Let

\[
 C_A:\quad y^a=A(t)
\]

with `A` squarefree over a characteristic-zero field.  The weighted
Wronskian equation

\[
 aAD'-bA'D=R                                             \tag{1}
\]

is equivalent to

\[
 d\!\left(\frac D{y^b}\right)
 =\frac{aAD'-bA'D}{a y^{a+b}}\,dt
 =\frac{R}{a y^{a+b}}\,dt.                              \tag{2}
\]

Thus the first block is an exactness problem with a prescribed class of
primitive, including support conditions on `D`.  For `(72,108)`,

\[
 (a,b,R)=(2,3,t^2),\qquad \deg A=8,
\]

so the smooth projective model has genus three.

## 2. Finite-pole Hermite reduction

For `m>a`, put `k=m-a`.  Given

\[
 \omega=\frac{P(t)}{y^m}\,dt,
\]

choose the unique polynomial `S` of degree less than `deg(A)` satisfying

\[
 S\equiv-\frac ak P(A')^{-1}\pmod A.                    \tag{3}
\]

Squarefreeness makes `(A')^{-1} mod A` well-defined.  Then

\[
 \frac P{y^m}\,dt
 =d\!\left(\frac S{y^k}\right)
  +\frac{P_1}{y^k}\,dt,                                 \tag{4}
\]

where

\[
 P_1=\frac{P-AS'+(k/a)A'S}{A}.                          \tag{5}
\]

The numerator in (5) is exactly divisible by `A`.  For a nontrivial
character, repetition lowers the denominator exponent to

\[
 r\equiv m\equiv b\pmod a,\qquad 1\le r<a.              \tag{6}
\]

There is a separate invariant branch when `b=0 mod a`.  Then `r=0` and

\[
 \frac{P(t)}{y^m}\,dt
 =\frac{P(t)}{A(t)^{m/a}}\,dt,
\]

so the differential descends to an ordinary rational differential on the
quotient line.  It must be treated by rational Hermite reduction and residue
tests there, rather than inserted into a nontrivial character space.  The
implementation returns this descended quotient differential explicitly.

For `(72,108)`, the two steps are `y^5 -> y^3 -> y`.  The first modular
inverse is the coordinate-free version of the archived interpolation step;
the second is the correction layer that exposes cohomology.

## 3. Reduction at infinity

At the final denominator `y^r`, polynomial degree is reduced with

\[
 d\!\left(Qy^{a-r}\right)
 =\left(AQ'+\frac{a-r}{a}A'Q\right)\frac{dt}{y^r}.       \tag{7}
\]

This leaves the unique affine remainder

\[
 \left(c_0+c_1t+\cdots+c_{n-2}t^{n-2}\right)\frac{dt}{y^r},
 \qquad n=\deg A.                                       \tag{8}
\]

Hence every nontrivial `mu_a` character has `n-1` affine obstruction
coordinates.  These spaces sum to

\[
 (a-1)(n-1)=\dim H^1_{\mathrm{dR}}(C_A\setminus\{\text{points at infinity}\}).
                                                                    \tag{9}
\]

The important correction to a naive generalization is that a fixed weighted
Wronskian sees the single character `r=b mod a`, not all of compact
`H^1_dR(C_A)`.  Hyperelliptic curves have only one nontrivial character, so
the distinction is invisible in `(72,108)`.

The displayed eigenspace language is literal after extending the ground
field to contain `mu_a`.  Over a characteristic-zero field that need not
contain those roots of unity, the construction is intrinsic: use the
corresponding eigensheaves after finite étale scalar extension, with Galois
descent.  None of the reduction formulas requires choosing a primitive
`a`-th root of unity.

## 4. Compact cohomology and residues at infinity

Let

\[
 \delta=\gcd(a,n).
\]

There are `delta` points over infinity and

\[
 2g=(a-1)(n-1)-(\delta-1).                              \tag{10}
\]

The character `r` contains one logarithmic residue direction precisely when

\[
 \frac a\delta\mid r.                                  \tag{11}
\]

Its compact dimension is therefore `n-2`; all other nontrivial characters
have compact dimension `n-1`.  Summing the character dimensions gives (10).

The implementation computes the residue functional exactly from the reverse
expansion

\[
 \left(\frac{A(t)}{\operatorname{lc}(A)t^n}\right)^{-r/a}
 \quad\text{at }t=\infty.                              \tag{12}
\]

For `(a,n,r)=(2,8,1)`, its pivot is `t^3dt/y`, and (12) gives exactly

\[
\begin{aligned}
0={}&c_3-\frac{a_7}{2}c_4
 +\left(\frac{3a_7^2}{8}-\frac{a_6}{2}\right)c_5\\
&+\left(-\frac{a_5}{2}+\frac{3a_7a_6}{4}
-\frac{5a_7^3}{16}\right)c_6.                         \tag{13}
\end{aligned}
\]

Deleting the pivot after imposing (13) gives the six compact obstruction
coordinates and the six basis representatives recorded in
[`WEIGHTED_WRONSKIAN_FIRST_BLOCK.md`](WEIGHTED_WRONSKIAN_FIRST_BLOCK.md).
This basis is canonical relative to the normalized Newton coordinate `t` and
the chosen monic model; it is not claimed to be intrinsic under arbitrary
changes of curve coordinates.

## 5. Supported primitives versus abstract exactness

Hermite reduction answers whether the differential class vanishes among
rational functions on the curve.  A Newton block asks the sharper question
whether it has a primitive

\[
 D/y^b
\]

with a prescribed finite support for `D`.  The engine therefore exposes two
layers:

1. `SuperellipticDeRham.reduce_weighted_wronskian(R,b)` computes the exact
   character-wise de Rham remainder.
2. `weighted_wronskian_compatibility(...)` solves only constant-pivot
   supported coefficients of `D` and returns the remaining compatibility
   equations without dividing by parameter strata.

For the normalized `(72,108)` support

\[
 A=t+a_2t^2+\cdots+a_7t^7+t^8,qquad
 D=d_2t^2+\cdots+d_{12}t^{12},
\]

the second layer solves all 11 `d_i` and returns six equations.  The first
layer identifies those six equations as compact de Rham coordinates after
clearing the modular-inverse denominator.  The old raw coefficient solve is
therefore retained only as a support certificate and a regression oracle.

There is now a denominator-free representative certificate for this
identification.  The supported solve gives

\[
 2AD_*'-3A'D_*-t^2=\sum_{i=0}^5 C_i t^{13+i}.
\]

The six tail differentials `t^(13+i)dt/y^5` form a compact de Rham basis:
the low-coefficient map on every possible primitive `D/y^3`, `deg D<=12`,
is triangular with nonzero determinant `prod_(k=0)^12(2k-3)`.  Hence the
compact coordinates are exactly `-C_i/2` in that basis.  The compiler
regression checks the identity, determinant, support bounds, and all six
coordinates.

## 6. Scaling quotient before elimination

For `(72,108)`, the endpoint normalization leaves

\[
 a_k\longmapsto\zeta^{k-1}a_k,\qquad \zeta^7=1.         \tag{14}
\]

The invariant chart

\[
 q=a_7^7,\qquad x_k=a_ka_7^{k-1}\quad(2\le k\le6)       \tag{15}

\]

turns the degree-35 eliminant into an irreducible `S_5` quintic in `q`; the
remaining degree seven is the Kummer equation `a_7^7=q`.  This validates the
order of operations:

\[
 \text{Newton bands}
 \longrightarrow \text{character de Rham remainder}
 \longrightarrow \text{residual-scaling invariants}
 \longrightarrow \text{reduced elimination}.           \tag{16}
\]

For a new chain the finite diagonal action must be derived from its own
endpoint normalization.  The order seven in (14) is not a universal feature.

## 7. Gauss--Manin connection

For a parameter `u` and a basis representative `p(t,u)dt/y^r`, differentiation
followed by the same Hermite reduction gives

\[
 \nabla_{\partial_u}\left[p\frac{dt}{y^r}\right]
 =\operatorname{Red}_r\left[
 \frac{A\partial_up-(r/a)p\partial_uA}{y^{a+r}}dt
 \right].
\]

`SuperellipticDeRham.gauss_manin_connection(...)` returns these matrices in
the deterministic compact basis.  A two-parameter elliptic regression checks
flatness and discriminant-supported poles.  A one-parameter degree-eight
plane-block slice additionally produces a dense `6x6` genus-three connection
matrix with poles only on its discriminant.  This promotes the fiberwise de
Rham reducer to a family-level engine; the plane-JC obstruction itself is a
generally non-horizontal section.

## 8. Frontier experiments

### `(96,144)`

The 2017 table contains the chain

\[
 (8,40)\longrightarrow(8,28)\longrightarrow(11/4,7)
\]

as a raw length-two complete-chain row.  The common-tail
divisibility calculation forces `q1=d0=4` and a residual cubic vertical
factor.  The two root partitions containing a simple root yield the forbidden
last lower corner `(8,4)`.  The only remaining row is

```text
R=kappa*x^2*y^7*(y-lambda)^3.
```

Translating that root gives the edge `(8,40)->(8,12)`.  Its complete-chain
length is at most three, and the permissive exact enumeration has open counts
`1,6,3,0` with no final corner.  Thus this repeated-tail row is excluded
before the following reuse data need to be recorded:

- the Laurent coordinate change and exact band supports at `(8,28)`;
- the leading pair `(a,b)` and right-hand monomial `R`;
- which coefficients have already been fixed by the added predecessor;
- the residual diagonal scaling weights.

The other five published `(96,144)` chains require separate band derivations
and are not eliminated by this no-escape result.

### `(75,125)`

The chain starts at `(5,20)` and ends at `(7/5,2)` with ratio `3:5`.  It is the
best architecture test because neither the old tail nor the old ratio is
available.  Derive its Laurent bands first, then feed the leading identity to
the engine.  The expected diagnostic is not necessarily genus three or six
conditions: it is the character dimension determined by the derived
`(a,b,deg A)` and formulae (10)--(11).

## 9. Reproduction

With the repository's pinned SymPy dependency installed, run

```bash
python3 plane-jc/cas/test_superelliptic_derham.py
```

Expected output:

```text
PASS: superelliptic Hermite/de Rham reduction
PASS: (72,108) gives 11 solved coefficients and 6 compact obstructions
PASS: character dimensions sum to compact H^1_deRham
PASS: trivial character descends to a rational quotient differential
PASS: elliptic Gauss--Manin matrices are flat with discriminant poles
PASS: genus-three plane slice has a dense 6x6 Gauss--Manin matrix
```

The compiler certificate is checked separately:

```bash
python3 plane-jc/cas/test_newton_derham_compiler.py
```

Expected output:

```text
PASS: (72,108) support equations are certified de Rham tail coordinates
PASS: the six tail classes have an exact triangular independence certificate
PASS: compiler IR refuses to invent missing frontier Laurent bands
```

The existing structural checker remains complementary:

```bash
python3 plane-jc/cas/weighted_wronskian_first_block.py
```

It verifies the archived `mu_7` grading, quotient quintic, irreducibility, and
Galois group.  The new engine does not rerun that Gröbner calculation.
