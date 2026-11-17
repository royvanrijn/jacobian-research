# The exceptional Cox atlas for the corrected \(F_{20}\) cover

## Status and outcome

The controlled-transform calculation now supplies the primitive value-one
variables that are absent from both the root-incidence ring and its affine
normalization.  Thirteen exact chart types cover all forty-eight positive
exceptional colors at the \(q\)-node, ramphoid \(r\)-cusp, conjugate triple
orbit, and cubic \(q\)-\(r\) tangency orbit.  The remaining five positive
colors are the already-certified generic \(d\)-, \(q\)-, and
\(r\)-ramification colors.

On every controlled chart, the compact divisors \(D_d,D_q,D_r\) have
principal local equations \(\tau^a,\tau^b,\tau^c\), with

\[
 3a+b+c=\nu_\tau(P_X).                              \tag{0.1}
\]

More strongly, if \(F\) is the exact strict transform, then

\[
 \frac{P_X}{(\tau^a)^3\tau^b\tau^c}=\frac{\partial F}{\partial Y}. \tag{0.2}
\]

Thus literal degree-\((3,1,1)\) derivative cancellation passes on every
exceptional chart, not merely at the level of divisors.

On the punctured normalized conductor cover there is also one explicit
parity-compatible distribution of the total derivative residue among the
three compact columns.  Its \((3,1,1)\) product is the previously certified
anti-invariant residue \(L(w)\).

This is a controlled **atlas**, not yet a global Cox ring.  The repository
does not yet contain all overlap maps gluing these principal local frames to
one regular SNC source and one global multisection algebra.  Global
\(H^0\)-membership of \(P_X\), full entrywise inverse-adjugate polynomiality,
and affine-space recognition therefore remain open.

The exact checker is
[`verify_f20_exceptional_cox_atlas.py`](../scripts/verify_f20_exceptional_cox_atlas.py),
with generated certificate
[`f20_exceptional_cox_atlas.json`](../artifacts/generated-results/f20_exceptional_cox_atlas.json).

## 1. General controlled-transform Cox theorem

Let \(A\) be a regular chart with an exceptional prime
\(E=(\tau=0)\).  For effective Cartier divisors

\[
 D_j=a_jE,
\]

their local multi-Rees algebra is

\[
 A[\tau^{a_1}T_1,\ldots,\tau^{a_m}T_m]
 \simeq A[Z_1,\ldots,Z_m].                          \tag{1.1}
\]

The full exceptional Cox atlas additionally contains the primitive variable

\[
 e=\tau T_E,                                       \tag{1.2}
\]

even when every global regular function pulled back from the affine
normalization has \(E\)-value at least two.  Formula (1.2) explains exactly
where the missing value-one column appears: it is a local frame of
\(\mathcal O(-E)\), not a new regular function on the normalization.

Suppose a hypersurface equation is pulled back by

\[
 X=X_0+\tau^rY
\]

and has exact factorization

\[
 P=\tau^NF(\tau,z,Y).                               \tag{1.3}
\]

Differentiating with respect to \(Y\) gives

\[
 \tau^rP_X=\tau^N F_Y,
 \qquad P_X=\tau^{N-r}F_Y.                         \tag{1.4}
\]

> **Controlled-transform Cox theorem.** If a proposed multidegree
> \(n=(n_j)\) has \(m=\sum_jn_ja_j\leq N-r\), then its local Cox monomial
> divides \(P_X\) literally, with quotient
> \(\tau^{N-r-m}F_Y\).  Equality \(m=N-r\) leaves no unused exceptional
> factor and gives the quotient \(F_Y\).  If \(m>N-r\), divisibility is
> impossible on that chart.

This theorem is formal once (1.3) is proved.  The substantive computation is
the exact divisibility of every pulled-back \(P\) and the compatibility of
the integers \(a_j\) across all colors.

## 2. The thirteen controlled chart types

The table records the exact strict-transform order \(N\), root-coordinate
weight \(r\), compact orders \((a,b,c)\), derivative order \(N-r\), and
number of geometric colors represented by the chart.

| chart | \(N\) | \(r\) | \((a,b,c)\) | \(N-r\) | colors |
|---|---:|---:|---:|---:|---:|
| \(q\)-node slopes | 2 | 1 | \((0,1,0)\) | 1 | 4 |
| \(r\)-cusp \(E_1\) | 5 | 1 | \((0,0,4)\) | 4 | 1 |
| \(r\)-cusp \(E_2\) | 10 | 2 | \((0,0,8)\) | 8 | 1 |
| \(r\)-cusp \(E_3\), unramified | 3 | 1 | \((0,0,2)\) | 2 | 1 |
| \(r\)-cusp \(E_3\), ramified | 5 | 1 | \((0,0,4)\) | 4 | 2 |
| \(r\)-cusp \(E_4\) | 5 | 1 | \((0,0,4)\) | 4 | 5 |
| triple \(E_1\), index four | 10 | 3 | \((1,2,2)\) | 7 | 2 |
| triple \(E_2\), quartic cluster | 4 | 1 | \((1,0,0)\) | 3 | 8 |
| \(qr\ E_1\), \(A\)-ramified | 3 | 1 | \((0,1,1)\) | 2 | 3 |
| \(qr\ E_1\), \(A\)-unramified | 2 | 1 | \((0,1,0)\) | 1 | 3 |
| \(qr\ E_1\), \(B\)-ramified | 2 | 1 | \((0,0,1)\) | 1 | 3 |
| \(qr\ E_2\), \(A\)-cluster | 3 | 1 | \((0,1,1)\) | 2 | 9 |
| \(qr\ E_2\), \(B\)-cluster | 2 | 1 | \((0,0,1)\) | 1 | 6 |
| **total** |  |  |  |  | **48** |

Every exceptional residual is generically separable.  The chart fields are
either rational or one of the three exact residue extensions

\[
 a^2-3a+1=0,\qquad b^2+i=0,\qquad
 8\alpha^3+16\alpha^2+2\alpha-7=0.                 \tag{2.1}
\]

The quadratic, quartic, and quintic residuals therefore encode whole
Galois-stable color packets without choosing numerical roots.

### 2.1 The missing triple value-one variable

At triple \(E_1\), the controlled substitution is

\[
\begin{aligned}
 s&=2i+\tau^4,\\
 t&=-\frac34+\frac i2+z\tau^4,\\
 X&=1+i+\tau^2(b+\tau Y),\qquad b^2+i=0.
\end{aligned}                                      \tag{2.2}
\]

Exact reduction gives

\[
 P=\tau^{10}F,\qquad P_X=\tau^7F_Y.                \tag{2.3}
\]

The three compact local coefficients are

\[
 Z_d=\tau T_d,qquad Z_q=\tau^2T_q,qquad
 Z_r=\tau^2T_r.                                    \tag{2.4}
\]

Consequently \(Z_d^3Z_qZ_r\) has coefficient \(\tau^7\).  The variable
\(e=\tau T_E\) is primitive and has value one, resolving the semigroup gap
found on the affine normalization.

At triple \(E_2\), only \(D_d\) occurs: \(P_X=\tau^3F_Y\) and
\(Z_d^3\) has coefficient \(\tau^3\).  At every cusp chart only \(D_r\)
occurs, with exponent exactly equal to \(\nu(P_X)\).  The five
\(q\)-\(r\) chart identities similarly reproduce the allocations
\(q^2r\), \(q\), and \(r\) from the residual discriminants.

## 3. A conductor-residue distribution

On the rational conductor cover put

\[
 A_-=w^2-2w+5,\qquad A_+=w^2+2w+5.
\]

Work on the punctured ring

\[
 C^\circ=\mathbf Q\left[w,
 \frac1{(w^2-1)(w^2+3)}\right].                    \tag{3.1}
\]

One exact choice of residue frames is

\[
 \rho_d=1,\qquad
 \rho_q=\frac{(3w^2+5)A_-A_+}{8(w^2-1)},\qquad
 \rho_r=\frac w4.                                  \tag{3.2}
\]

They satisfy

\[
 \rho_d(-w)=\rho_d(w),\qquad
 \rho_q(-w)=\rho_q(w),\qquad
 \rho_r(-w)=-\rho_r(w),                            \tag{3.3}
\]

and, exactly,

\[
 \rho_d^3\rho_q\rho_r
 =\frac{w(3w^2+5)A_-A_+}{32(w-1)(w+1)}
 =L(w).                                             \tag{3.4}
\]

Thus \(D_d,D_q,D_r\) carry characters \((+,+,-)\), whose product is the
anti-invariant derivative character.  The previous unit-plus-selector
completion has determinant \(-1\).

This distribution is not unique: invariant units can move between the
three frames, and even factors can be redistributed.  The choice (3.2) is
packet-compatible in a useful minimal sense.  The \(q\)-node and
\(q\)-\(r\)-tangency factors are coprime to all three numerators; the triple
factor \(A_-A_+\) lies in \(\rho_q\); the transverse factor \(w\) lies in
\(\rho_r\); and the \(r\)-cusp is off this conductor.  Extending these
frames through the controlled charts requires explicit overlap maps and is
the remaining global gluing problem.

## 4. What has and has not passed

The compact degree \((3,1,1)\) passes two stronger tests than the earlier
valuation matrix:

1. on each of the thirteen exceptional chart types, \(P_X\) is literally
   the compact Cox monomial times \(F_Y\);
2. on the punctured conductor cover, the three residue frames multiply
   literally to \(L(w)\) with the correct involution character.

These statements do not yet prove

\[
 P_X\in
 H^0\!\left(\widetilde X,
 \mathcal O(-3D_d-D_q-D_r)\right)                  \tag{4.1}
\]

for one global regular \(\widetilde X\), because that scheme and the overlap
transition functions have not been assembled.  Accordingly:

- the local derivative-denominator gate passes on all positive colors;
- global entrywise inverse-adjugate polynomiality is not reached;
- affine-space recognition is not reached.

The subsequent
[`F20_EXCEPTIONAL_COX_CORNERS.md`](F20_EXCEPTIONAL_COX_CORNERS.md) constructs
seven exact two-parameter corner types, proves the complementary \(q\)-node
transition, and deletes the false triple \(E_1\)-\(E_2\) edge by a
root-center separation certificate.  The remaining exact object is the Čech
transition matrix for the strict-boundary attachments, generic ramification
charts, and conductor-to-corner charts.  Its multidegree-\((3,1,1)\)
cocycle must be trivial before a global Cox algebra or downstream affine
claim can be made.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_f20_exceptional_cox_atlas.py \
  --output artifacts/generated-results/f20_exceptional_cox_atlas.json
```

The checker proves exact strict-transform divisibility and derivative
identities for all thirteen chart types, verifies generic separability of
every exceptional residual, counts the forty-eight positive exceptional
colors, and verifies (3.2)--(3.4).  It uses exact SymPy arithmetic over the
three residue fields in (2.1).
