# The degree-forty-two relative Ritt cotangent cone

Work over \(\mathbb Q\) on the normalized complete-decomposition chart

\[
 2\circ7\circ3.
\]

The two half-braids to \(3\circ7\circ2\) omit respectively the composite
cut \(6\) and the prime cut \(7\).  This note constructs their completed
path-to-boundary conormal flag along the Dickson component and separates the
sector defect from the common spectator layer.

The result corrects both naive possibilities suggested by the first tangent
calculation.  The raw degree-forty-two sector is not the degree-thirty Artin
defect unchanged, but it is also not an extra higher-dimensional branch.
It is a thicker divisor-supported formal layer.

## 1. Ideals and graph-normal coordinates

Let \(S\) be the nine-variable coefficient ring of the
\(2\circ7\circ3\) factor chart.  Write

\[
 I_6\subset S,\qquad I_7\subset S,\qquad I_\partial\subset S,\qquad
 K\subset S
\]

for the composite-omitting path, prime-omitting path, full braid boundary,
and reduced Dickson graph ideals.  The ideal convention is contravariant:
all path schemes contain the Dickson surface.

At the monomial point choose graph-normal coordinates

\[
 \widehat S=\mathbb Q[[n_1,\ldots,n_7,\tau,z]],
                                                               \tag{1.1}
\]

where \(\tau=t-1\), \(z\) is the Dickson parameter, and
\(K=(n_1,\ldots,n_7)\).  The coordinate map is polynomial.  The inner cubic
coefficients become

\[
 3((1+\tau)^2-z),\qquad 3(1+\tau),
\]

and every other factor coefficient is its Dickson value plus one normal
coordinate.  Thus (1.1) is an exact triangular coordinate change, not a
formal approximation.

All calculations below use the `7 normal | 2 base` block order and the
ordinary polynomial residuals, including the degree-one coefficient.

## 2. The exact ideal flag

Exact characteristic-zero reduction gives

\[
 \boxed{I_6\subsetneq I_7=I_\partial\subsetneq K.}            \tag{2.1}
\]

The equality \(I_7=I_\partial\) is the first clean separation statement:
the path omitting the prime cut already contains the entire common
spectator layer.  Adding the missing prime cut changes nothing.

All three distinct schemes in (2.1) have dimension two:

\[
 \dim S/I_6=\dim S/I_\partial=\dim S/K=2.                    \tag{2.2}
\]

Consequently the extra degree-forty-two structure is supported over the
Dickson component; it is not a three-dimensional branch passing through it.

## 3. The conormal flag

At \((n,\tau,z)=0\), the Jacobian ranks of
\(I_6,I_7=I_\partial,K\) are

\[
 5,\quad6,\quad7.                                            \tag{3.1}
\]

Equivalently, their tangent dimensions are \(4,3,2\).  Hence the completed
conormal flag has two one-dimensional successive quotients:

\[
 \frac{I_\partial+\mathfrak m^2}
      {I_6+\mathfrak m^2}
 \cong\mathbb Q,
 \qquad
 \frac{K+\mathfrak m^2}
      {I_\partial+\mathfrak m^2}
 \cong\mathbb Q.                                             \tag{3.2}
\]

The first is the **relative sector direction**.  The second is the
**spectator direction** shared by the prime-omitting path and the full
boundary.

Put

\[
 A_6=\widehat S/I_6,\qquad
 A_\partial=\widehat S/I_\partial,\qquad
 B=\widehat S/K.
\]

The quotient maps \(A_6\to A_\partial\to B\) give the canonical cotangent
transitivity triangle

\[
 B\otimes^{\mathbf L}_{A_\partial}L_{A_\partial/A_6}
 \longrightarrow
 L_{B/A_6}
 \longrightarrow
 L_{B/A_\partial}
 \overset{+1}{\longrightarrow}.                             \tag{3.3}
\]

Thus the separation is intrinsic even without a splitting: the left term
is the path-to-boundary sector cone and the right term is the
boundary-to-Dickson spectator cone.

## 4. Exact base-annihilation exponents

The two ideal quotients underlying the first conormal layers are

\[
 M_{\rm sec}=I_\partial/I_6,\qquad
 M_{\rm sp}=K/I_\partial.                                    \tag{4.1}
\]

Repeated exact reduction against the cached standard bases gives

\[
\begin{aligned}
 \min\{m\ge0:z^m I_\partial\subset I_6\}&=8,\\
 \min\{m\ge0:z^m K\subset I_\partial\}&=1.
\end{aligned}                                                \tag{4.2}
\]

Therefore

\[
 z^8M_{\rm sec}=0,\quad z^7M_{\rm sec}\ne0,
 \qquad
 zM_{\rm sp}=0,\quad M_{\rm sp}\ne0.                         \tag{4.3}
\]

Both layers are supported on the monomial divisor \(z=0\), but their
thicknesses are different.  The spectator is a first-order
divisor-supported layer.  The relative sector persists through seven
further \(z\)-adic orders.

This is the precise sense in which the spectator has been separated.  It is
not obtained by subtracting tangent dimensions: it is the quotient
\(K/I_\partial\) in the exact completed ideal flag, with annihilator exponent
one.  The sector is the preceding quotient \(I_\partial/I_6\), with exponent
eight.

## 5. Completed local jets

Let

\[
 \mathfrak m=(n_1,\ldots,n_7,\tau,z).
\]

Exact vector-space lengths of the first three completed jets are

\[
\begin{array}{c|ccc}
q&
\ell\bigl(\widehat S/(I_6+\mathfrak m^q)\bigr)&
\ell\bigl(\widehat S/(I_\partial+\mathfrak m^q)\bigr)&
\ell\bigl(\widehat S/(K+\mathfrak m^q)\bigr)\\
\hline
1&1&1&1\\
2&5&4&3\\
3&14&9&6.
\end{array}                                                  \tag{5.1}
\]

The prime-omitting path has the same column as \(I_\partial\), by (2.1).
The successive differences

\[
\begin{array}{c|ccc}
q&1&2&3\\
\hline
\text{sector }(I_\partial/I_6)&0&1&5\\
\text{spectator }(K/I_\partial)&0&1&3
\end{array}                                                  \tag{5.2}
\]

show that the two rank-one conormal directions acquire different higher
relations immediately.  They are not two copies of one universal
square-zero module.

## 6. Comparison with degree thirty

For the degree-thirty cut-\(6\) sector, the full boundary is already the
reduced Dickson surface and the defect annihilator is \((z^4)\).  In degree
forty-two, after removing the common spectator layer, the relative sector
has exact exponent \(8\).

Thus the following data survive changing the middle prime \(5\to7\):

* the composite-omitting orientation;
* one excess conormal direction;
* support on the monomial divisor;
* equality of the prime-omitting path with the full boundary.

The \(z\)-adic thickness does **not** survive:

\[
 4\quad\longrightarrow\quad8.                               \tag{6.1}
\]

Only two spectator degrees have been computed, so no general exponent
formula is claimed.  The values are compatible with \(2(p-3)\) for middle
prime \(p\), but this is presently a conjectural pattern, not a theorem.

## 7. Derived scope

Equations (2.1)--(5.2) construct the completed ideal flag, its conormal
quotients, their exact base support, and the cotangent transitivity triangle.
They do not prove that (3.3) splits, do not identify the full higher
homology of either relative cotangent complex, and do not compute the
\(L_\infty\) extension class coupling the two layers.

The next derived calculation is therefore finite and specific:

1. present \(I_\partial/I_6\) and \(K/I_\partial\) as
   \(\mathbb Q[[\tau,z]]\)-modules;
2. compute their first syzygies and the induced two-term relative cotangent
   complexes;
3. evaluate the extension class in the transitivity triangle (3.3);
4. compare that class with the degree-thirty cut-\(6\) sector after replacing
   \(z^4\) by \(z^8\).

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_degree42_ritt_relative_cotangent_cone.py
```

The checker constructs the ordinary polynomial path ideals, transports them
to (1.1), caches the four standard bases, verifies (2.1)--(3.2), finds the
minimal exponents in (4.2) by repeated exact reduction, and computes the
jets (5.1).
