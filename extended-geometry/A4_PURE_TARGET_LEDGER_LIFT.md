# The pure-target \(A_4\) ledger lift

## 1. Outcome

The residual \(WL\) ledger of the polynomial \(A_4\) cone can be used
polynomially.  Adjoin one source coordinate \(z\) and define

\[
\boxed{
\widehat\Phi(U,V,W,z)=
\left(
WN_1,\ WN_2,\ WH,\ \frac{WL}{4}z
\right).
}                                                     \tag{1.1}
\]

Let

\[
\begin{aligned}
\mathcal B(P,Q,R)
={}&P^3-3PQ^2+2Q^3-9PQR+9Q^2R\\
   &-27PR^2+27QR^2+27R^3.
\end{aligned}
\]

The target-ledger identity

\[
\mathcal B(WN_1,WN_2,WH)=W^3K^3L^2
\]

and \(\det D\Phi=4W^2K^3L\) give

\[
\boxed{
\det D\widehat\Phi
=W^3K^3L^2
=\mathcal B\bigl(\widehat\Phi_1,
                  \widehat\Phi_2,
                  \widehat\Phi_3\bigr).
}                                                     \tag{1.2}
\]

Thus the determinant defect is now a pure target pullback.  Equivalently,
for the target logarithmic volume form

\[
\Omega_{\mathcal B}
=\frac{dP\wedge dQ\wedge dR\wedge dS}
       {\mathcal B(P,Q,R)},
\]

one has

\[
\boxed{
\widehat\Phi^*\Omega_{\mathcal B}
=dU\wedge dV\wedge dW\wedge dz.
}                                                     \tag{1.3}
\]

This is a polynomial log-Keller realization on affine four-space.  It is
not yet an ordinary Keller map because \(\mathcal B\) is nonconstant.

## 2. The \(A_4\) cover is unchanged

Write

\[
(P,Q,R,S)=\widehat\Phi(U,V,W,z).
\]

The first three coordinates are the original cone map.  At the generic
point,

\[
 z=\frac{4S}{WL}.                                    \tag{2.1}
\]

Therefore adjoining \(z\) and \(S\) adds the same purely transcendental
coordinate to the source and target function fields:

\[
\mathbb Q(U,V,W,z)
=\mathbb Q(U,V,W)(S),
\]

and

\[
\mathbb Q(P,Q,R,S)
=\mathbb Q(P,Q,R)(S).
\]

Consequently

\[
\operatorname{gdeg}(\widehat\Phi)=4,
\]

and the geometric and arithmetic generic inverse monodromy remain the
natural \(A_4\)-action.  Over any target point with
\(\mathcal B(P,Q,R)\ne0\), the fourth coordinate gives exactly one \(z\)
above each point of the original four-sheet fiber.

## 3. Why the obvious target factorization still fails

The standard polynomial target chart with Jacobian \(\mathcal B\) is

\[
\beta(P,Q,R,T)=(P,Q,R,\mathcal B(P,Q,R)T).           \tag{3.1}
\]

If \(\widehat\Phi\) factored polynomially as

\[
\widehat\Phi=\beta\circ F,
\]

then the first three outputs of \(F\) would remain \(\Phi\), and its fourth
coordinate would have to be

\[
F_4
=\frac{(WL/4)z}{W^3K^3L^2}
=\frac{z}{4W^2K^3L}.                                \tag{3.2}
\]

This is not polynomial.  Hence the equality of Jacobian divisors (1.2)
does not by itself supply a polynomial factorization through the
multiplicative target modification.

More generally, any factorization which retains \(P,Q,R\) as the first
three target coordinates is block triangular.  Its determinant remains
divisible by \(\det D\Phi\).  A successful factorization must alter at
least two of the three cone outputs.

## 4. The three-puncture factor has a double-incidence encoding

The residual cubic admits a polynomial encoding adapted to its two finite
puncture characters.  Put

\[
r=\frac{U-V}{3}
\]

and retain \(V\).  Then

\[
\ell=\frac{L}{27}
=r(r-1)V+r^3-3r+1.                                  \tag{4.1}
\]

Define

\[
\begin{aligned}
u&=3-r^2-(r-1)V,\\
v&=rV+r^2+r-2.
\end{aligned}                                       \tag{4.2}
\]

Direct calculation gives the two incidence identities

\[
\boxed{
1-ru=\ell,\qquad
1-(r-1)v=-\ell.
}                                                     \tag{4.3}
\]

On \(L=0\),

\[
u=r^{-1},\qquad v=(r-1)^{-1},                       \tag{4.4}
\]

so the two independent units of

\[
\mathbb Q[r,r^{-1},(r-1)^{-1}]
\]

already have polynomial ambient representatives.  This explains the
power \(L^2\) in the target pullback \(\mathcal B(\Phi)\): the two
finite-puncture incidence defects are the same plane divisor with opposite
orientation.

It also identifies a restriction.  The separated double-incidence
primitive from the existing three-puncture core has determinant

\[
(1-ru)^a(1-(r-1)v)^b
=(-1)^b\ell^{a+b},
\qquad a,b\ge1.                                     \tag{4.5}
\]

Thus a separated primitive which records both puncture characters begins
with exponent two.  It naturally fits the \(L^2\) target ledger, but cannot
by itself realize the single relative-canonical copy of \(L\) in
\det D\Phi.  The source and target modifications must interact.

## 5. At least two cone outputs must be coupled

Every cone coordinate has the form

\[
F_i=W A_i(U,V).
\]

For any two constant linear combinations \(f,g\) of the three cone
coordinates,

\[
df\wedge dg
\]

is divisible by \(W\).  Therefore, in a polynomial stabilization of any
dimension, retaining two independent cone combinations as output
coordinates forces the full Jacobian to be divisible by \(W\).

Hence:

\[
\boxed{\text{An ordinary Keller completion must modify at least two of
the three cone outputs.}}                            \tag{5.1}
\]

This is stronger than saying that the added variables must feed back into
the cover block: one modified output is still insufficient.

## 6. Constant-direction masks give no new construction

Let \(x=(U,V,W)\), let \(z=(z_1,\ldots,z_m)\), and consider an
affine-linear mask with constant coefficient matrix:

\[
G(x,z)=H(x)+Cz,
\]

where \(G\) has \(3+m\) outputs and \(C\) has rank \(m\).  A constant target
linear change splits off the \(m\) mask directions.  The remaining three
outputs are three polynomial functions \(G_0(x)\), and

\[
\det DG=\text{constant}\cdot\det DG_0,\qquad
\operatorname{gdeg}(G)=\operatorname{gdeg}(G_0).     \tag{6.1}
\]

Thus a constant-direction mask is only a stabilization of a
three-dimensional polynomial map.  If it produced a degree-four
\(A_4\) Keller map, that map was already present in the three unmasked
coordinates.

In particular, the next search must use source-dependent mask coefficients.
Together with (5.1), the smallest genuinely new skeleton has two auxiliary
variables entering at least two cone outputs with coefficients depending
on \(U,V,W\).

## 7. The surviving factorization problem

The problem is now a polynomial factorization problem, rather than an
unorganized determinant search.

Construct polynomial maps

\[
F:\mathbb A^{3+m}\longrightarrow\mathbb A^{3+m},
\qquad
\beta:\mathbb A^{3+m}\longrightarrow\mathbb A^{3+m}
\]

such that:

1. \(\det DF\in\mathbb Q^*\);
2. \(\det D\beta=\mathcal B\) in a target-ledger chart;
3. a coupled analogue of
   \[
   \widehat\Phi=\beta\circ F
   \]
   holds after eliminating the modification variables;
4. at least two cone outputs are modified;
5. the mask coefficients are source-dependent;
6. elimination recovers the original extension
   \[
   \mathbb Q(U,V,W)/\mathbb Q(WN_1,WN_2,WH),
   \]
   rather than exposing \(U,V\), or \(W\);
7. the two puncture characters in (4.3) remain distinct.

The minimal live ansatz therefore starts in dimension five with two
source-dependent masks.  A one-variable or constant-mask search can be
discarded before coefficient elimination.

The [two-mask factorization frontier](A4_TWO_MASK_FACTORIZATION_FRONTIER.md)
constructs the puncture-adapted birational target blowdown with Jacobian
\(\mathcal B\), excludes every one-primitive triangular rechart, and checks
all \(120\) assignments of the existing five outputs to its base and mask
coordinates.  The first surviving coefficient system is its conjugation by
a source-dependent two-mask shear.

An ordinary solution would be a degree-four counterexample to the
Jacobian conjecture.  Identity (1.2) is nonetheless a genuine positive
advance: all source ramification has been converted into one explicit
target divisor while preserving exact \(A_4\) monodromy.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_a4_pure_target_ledger.py
```

The checker verifies the pure-target determinant identity, logarithmic
volume cancellation, preservation of the generic fourth coordinate, the
double-incidence formulas, and the \(W\)-divisibility of every pair of cone
coordinate differentials.
