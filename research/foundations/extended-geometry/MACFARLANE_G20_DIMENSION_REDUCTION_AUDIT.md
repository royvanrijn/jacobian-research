# MacFarlane \(G_{20}\): exact dimension-reduction audit

## Status

This note does **not** construct a 19-variable cubic-homogeneous Keller
counterexample.  It identifies two exact construction gates that would do so
and closes the direct linear quotient and hyperplane-restriction routes for
MacFarlane's \(G_{20}\).

The external input is A. MacFarlane's
[`keller-counterexamples-13-20` at commit
`dad6090`](https://github.com/Amacfa/keller-counterexamples-13-20/tree/dad6090bf4f01b3cdad04048fbe16f3be52b485c)
certificate.  In its notation,

\[
 F_{13}(x)=x+R(x)+B\gamma(x)
\]

has degrees one through three and a rational collision, while

\[
 G_{20}(x,w,\tau)
 =
 \bigl(x+\tau R(x)+\tau^2Bw,\;w-\gamma(x),\;\tau\bigr)                 \tag{1}
\]

is a cubic-homogeneous Keller collision.  The checker below reconstructs
both collisions from the displayed rational formulas.  It does not replace
the external determinant certificate.

## 1. Why the count is \(20\)

The six components of \(\gamma\) are linearly independent, and \(B\) has
rank six.  Hence the cubic component vector \(B\gamma\) of \(F_{13}\) has
coefficient-row rank six.  Formula (1) is exactly the repository's
rank-compressed homogenization

\[
 n+\operatorname{rank}(C)+1=13+6+1=20.                              \tag{2}
\]

Consequently either of the following is sufficient for dimension 19:

1. construct a collision-preserving degree-three source in dimension 12
   whose cubic-output rank is at most six;
2. find a stable-equivalent 13-variable source whose cubic-output rank is at
   most five.

This is the useful positive reduction.  It couples the cubic-homogeneous
problem to the separate \(13\to12\) degree-three problem without conflating
the two classes.

## 2. Constant-kernel quotient is already terminal

Let \(H_{20}=G_{20}-I\).  Exact coefficient comparison gives

\[
 \bigcap_z\ker JH_{20}(z)=0.                                         \tag{3}
\]

Thus the constant-kernel quotient used in the repository's \(24\to22\) and
\(24\to21\) reductions removes nothing from \(G_{20}\).  The 20 components
of \(H_{20}\) span a 19-dimensional polynomial space; its only fixed linear
covector is \(\tau\).  For \(F_{13}-I\), the corresponding common input
kernel is also zero and the 13 nonlinear components have full output span.

These facts exclude collision-preserving linear semiconjugate quotients of
the standard constant-kernel type.  They do not exclude nonlinear
semiconjugacies.

## 3. No Keller hyperplane restriction of \(G_{20}\)

There is a slightly broader linear possibility than a constant-kernel
quotient: restrict \(G_{20}\) to an invariant hyperplane through its
collision.

Let

\[
 L=a\mathbin{\cdot}x+b\mathbin{\cdot}w+c\tau
\]

and suppose \(\{L=0\}\) is invariant.  Since the full determinant is one,
the tangent restriction can have constant Jacobian only if

\[
 L(H_{20})\ \text{is divisible by}\ L^2.                             \tag{4}
\]

Indeed, invariance first gives \(L(H_{20})=Lq_2\).  Along \(L=0\), the
normal Jacobian factor is \(1+q_2\), while the tangent determinant is its
reciprocal.  Constancy and evaluation at the origin force
\(q_2|_{L=0}=0\), which is (4).

Now

\[
 L(H_{20})
 =\tau\,aR+\tau^2aBw-b\gamma.                                       \tag{5}
\]

The right side has degree at most one in \(w\).  If \(b\ne0\), a nonzero
polynomial divisible by \(L^2\) would have \(w\)-degree at least two.
The correction components have output span 19, so the exceptional zero
polynomial would make \(L\) a multiple of the sole fixed covector \(\tau\);
that hyperplane misses the collision.  Hence \(b=0\).

The collision again excludes \(L\) proportional to \(\tau\), so
\(\gcd(L,\tau)=1\).  Equation (4) and (5) imply

\[
 aB=0,\qquad aR=d(a\mathbin{\cdot}x+c\tau)^2                         \tag{6}
\]

for a constant \(d\).  If \(d\ne0\), comparison of \(\tau\)-terms gives
\(c=0\).  The only square monomials occurring in a linear combination of
the components of \(R\) are \(x_1^2,x_2^2\), so (6) forces
\(a\in\langle e_1,e_2\rangle\).  But \(aB=0\) kills both coefficients.
If \(d=0\), exact row reduction gives

\[
 \ker_{\rm left}R
 =
 \left\langle e_6,\;-\frac13e_5+e_{13}\right\rangle,
\]

whose intersection with \(\ker_{\rm left}B\) is zero.  Thus no such
hyperplane exists.

## 4. The analogous \(F_{13}\) slice and low-degree invariants

The same normal-factor argument applies to an affine hyperplane
\(\ell=a\cdot x+d=0\) for \(F_{13}\).  If \(d\ne0\), divisibility of the
nonlinear polynomial \(a(F_{13}-I)\) by \(\ell^2\), together with its zero
constant and linear parts, forces that polynomial to vanish.  Full nonlinear
output span then gives \(a=0\).

For \(d=0\), the first collision point forces \(a_3=0\).  At quadratic
order, either \(aR\) is a nonzero multiple of \((a\cdot x)^2\), in which
case square support reduces to \(a_1,a_2\) and the second collision plus the
terms \(x_{11}x_{12}\) and \(x_8x_9\) force \(a=0\); or \(aR=0\).  In the
second case the two collision equations cut the displayed left kernel of
\(R\) to

\[
 a=-2e_5+13e_6+6e_{13}.
\]

Its cubic defect is

\[
 aB\gamma=13x_1^2x_2-2x_1x_2x_3,
\]

which cannot be divisible by \((a\cdot x)^2\).  Therefore this affine
restriction route also fails.

Finally, the equivariant pullback calculation for
\(P\mapsto P\circ F_{13}\) uses the exact torus weights

\[
 (-1,1,2,0,2,-1,3,0,1,2,1,-2,2).
\]

Full-column rank modulo \(1000003\) in every weight sector supplies nonzero
minors over \(\mathbb Q\), and shows

\[
 \{P:\deg P\le3,\ P\circ F_{13}=P\}=\mathbb Q.                       \tag{7}
\]

Equation (7) closes only polynomial invariant coordinates through degree
three.  Higher-degree invariants and general nonlinear quotients remain open.

## 5. Reproduction and next experiment

The generic restriction/cancellation theorem, pair-aware collision policy,
the classification of every homogenizing slice, and the two distinct
terminal objectives are now canonicalized in
[`BACKWARD_CUBIC_REDUCTION.md`](BACKWARD_CUBIC_REDUCTION.md).  The present
note remains the map-specific obstruction audit.

Run

```bash
.venv/bin/python scripts/audit_macfarlane_g20_dimension_reduction.py
.venv/bin/python scripts/audit_macfarlane_f13_low_degree_invariants.py
```

The scripts write
[`macfarlane_g20_dimension_reduction_audit.json`](../artifacts/generated-results/macfarlane_g20_dimension_reduction_audit.json)
and
[`macfarlane_f13_low_degree_invariants.json`](../artifacts/generated-results/macfarlane_f13_low_degree_invariants.json).
They verify the collision formulas, ranks, common kernels, hyperplane
linear algebra, and every pullback sector in (7).

The next search should not widen the old monomial beam or begin from the
expanded \(F_{13}\) alone.  It should replay the upstream homogeneous
restriction, dehomogenization, and companion-block cancellation as
first-class operations.  In particular it should:

1. use Thompson \(G_{24}\to\) MacFarlane \(G_{20}\to F_{13}\) as a
   regression case;
2. compute collision-compatible left kernels before dehomogenization and
   right kernels after every homogenization;
3. track the three collision pairs separately, since noninjectivity requires
   only one surviving pair;
4. give every \(B,c\) companion block an explicit lifetime ending at
   \(t=1\), where
   \[
   (x+Q+By,\;y-c)
   =
   A_B\circ((x+Q+Bc)\times I)\circ S_c
   \]
   cancels all companion variables by stable left--right equivalence; and
5. optimize the lexicographic target

   \[
   (n,n+\operatorname{rank}C+1,\operatorname{rank}C)
   \]

   for the direct degree-three archive, while retaining the separate
   homogeneous objective
   \[
   (n+\operatorname{rank}C+1,n,\operatorname{rank}C).
   \]

A direct dimension-twelve endpoint succeeds regardless of cubic-output
rank.  The secondary dimension-nineteen homogeneous gates are
\((n,\operatorname{rank}C)=(12,\le6)\) and \((13,\le5)\).  Modular candidates
do not change either endpoint; an exact collision, determinant bridge, and
independent sparse replay are required.

The generic theorem and its MacFarlane calibration are recorded as `BCR1`
in `MATH_STATUS.json`.  The map-specific obstruction results in this audit
remain supporting computations rather than a separate status entry.
