# All-degree rational fibers for weighted Keller maps

This note extracts a uniform consequence of the
[weighted-suspension theorem](TANGENT_MAP_CORE.md).  It replaces the bounded
integer-root atlas by one construction valid in every inverse degree.

## Theorem

For every integer `N>=3`, there is a polynomial Keller map

\[
 F_N:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q}
\]

with a target whose complete fiber consists of exactly `N` distinct rational
points.  Over the reals, that target has an open neighborhood on which `F_N`
has `N` disjoint real inverse sheets.

## Integer-root seeds

Put `N=2k+1` in odd degree and `N=2k` in even degree.  Define the following
sets of `N-2` integers:

\[
 \begin{aligned}
 S_{2k+1}&=\{2\}\mathbin\cup
   \{j,1-j:3\le j\le k+1\},\\
 S_{2k}&=\{3,4\}\mathbin\cup
   \{j,1-j:5\le j\le k+2\}.
 \end{aligned}
\]

Empty ranges are omitted.  Thus the formulas include `N=3` and `N=4`.
Let

\[
 G_N(W)=W(W+1)\prod_{\rho\in S_N}(W-\rho),\qquad
 \lambda_N=G_N'(0)-G_N'(1),
\]

and normalize

\[
 H_N(W)=\frac{G_N(W)-G_N'(0)W}{\lambda_N}.                 \tag{1}
\]

The roots `0,-1` and the elements of `S_N` are pairwise distinct integers,
so `G_N` has exactly `N` simple rational roots.

## Normalization and admissibility

The paired factors give `G_N(1)=G_N'(0)`.  If

\[
 D_N=\frac{G_N'(1)}{G_N(1)},
\]

then direct logarithmic differentiation gives

\[
 D_{2k+1}=\frac1{k+1},\qquad
 D_{2k}=\frac5{12}+\frac1{k+2}.
\]

The corresponding endpoint products are

\[
 G_{2k+1}'(0)=(-1)^k(k+1)(k!)^2,qquad
 G_{2k}'(0)=\frac{(-1)^k}{12}(k+2)((k+1)!)^2.
\]

Consequently

\[
 \boxed{\lambda_{2k+1}=(-1)^k k(k!)^2},\qquad
 \boxed{\lambda_{2k}=\frac{(-1)^k}{144}
             (7k+2)((k+1)!)^2}.                           \tag{2}
\]

In particular, the scale never vanishes.  Equations (1)--(2) imply

\[
 H_N(0)=H_N'(0)=H_N(1)=0,\qquad H_N'(1)=-1.
\]

It remains only to exclude the weighted-chart value `H_N''(1)=-2`.
Write `\mathcal H_m^{(2)}=\sum_{i=1}^m i^{-2}`.  In odd degree, logarithmic
differentiation a second time gives
`Q_{2k+1}=\mathcal H_k^{(2)}+\mathcal H_{k+1}^{(2)}` and hence the exact
expression

\[
 H_{2k+1}''(1)=-\frac{2(k+1)}k\mathcal H_k^{(2)}<-2.      \tag{3}
\]

For even degree put

\[
 Q_{2k}=\sum_{r:\,G_{2k}(r)=0}\frac1{(1-r)^2}
 =\frac{29}{18}+
   \sum_{j=5}^{k+2}\left(\frac1{(j-1)^2}+\frac1{j^2}\right).
\]

Since `5/12<D_{2k}<=2/3`, one has

\[
 Q_{2k}\ge\frac{29}{18}>
 1+\left(\frac7{12}\right)^2
 \ge 1+(1-D_{2k})^2.
\]

Using `G_N''(1)/G_N(1)=D_N^2-Q_N` and
`lambda_N=G_N(1)(1-D_N)` now yields

\[
 H_{2k}''(1)=\frac{D_{2k}^2-Q_{2k}}{1-D_{2k}}<-2.         \tag{4}
\]

Thus every `H_N` is an admissible weighted seed with `c=b_0=1`.  The
weighted-suspension theorem supplies a polynomial map `F_N` over `Q` with
`det DF_N=1` and inverse pencil

\[
 E_{A,B,C}(W)=H_N(W)-BCW+AC^2.                            \tag{5}
\]

## The complete rational fiber

Take

\[
 (A_N,B_N,C_N)=
 \left(0,-\frac{G_N'(0)}{\lambda_N},1\right).
\]

Then (1) and (5) give

\[
 E_{A_N,B_N,C_N}(W)=\frac{G_N(W)}{\lambda_N}.             \tag{6}
\]

All `N` roots are simple and rational.  Because `C_N=1`, the reconstruction
theorem applies to every root.  Explicitly,

\[
 \gamma=-E'(W),\qquad x=\frac1\gamma,qquad
 y=W-\gamma,qquad
 z=\frac{\gamma-1-a_0(W/\gamma-1)}{x^2},
\]

where `a_0=-(1+H_N''(1))/(2+H_N''(1))`.  These expressions are rational and
finite because `E'(W)\ne0`.  Distinct roots give distinct source points since
`W=(1+xy)\gamma` is recovered from the source.  Conversely, reconstruction
over `C!=0` identifies the entire fiber with the simple roots of (6), so
there are no additional source points.

Finally `det DF_N=1`, so the real inverse-function theorem supplies a real
inverse branch around each of the `N` points.  After shrinking their common
target neighborhood, the roots of (5) remain distinct and real.  The degree
`N` pencil and reconstruction therefore give exactly `N` disjoint real
sheets there.

## Verification

Run

```bash
.venv/bin/python scripts/verify_all_degree_rational_fibers.py
```

The checker audits the symbolic all-`k` product, logarithmic-derivative, and
inequality identities, then retains the former degrees `3,...,100` as a
finite exact regression rather than as the proof of the theorem.
