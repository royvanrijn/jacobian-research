# F2 `j=1` derivation audit for `(75,125)`

> **Status: forced skeleton, not a normal-form certificate.**  The published
> complete-chain theorem determines two consecutive edges and one Puiseux
> translation.  It does not determine the rest of either Laurent polygon.
> Consequently this note does not claim to eliminate `(75,125)`, and the
> compiler must continue to reject it as a complete front end.

The executable record is
[`cas/f2_75_125_frontend.py`](cas/f2_75_125_frontend.py).  It emits the
machine-readable residual certificate
`plane-jc.f2-75-125-residual.v1`.

## 1. What the F2 table and chain theorem force

Family `F2` at `j=1` has

\[
 A_0=(5,20),\quad A'_0=(1,0),\quad
 A_1=(7/5,2),\quad (m,n)=(3,5).
\]

Thus `dir(A0-A0')=(5,-1)`.  The complete-chain theorem passes from
`L=K[x,y]` to `L^(5)` by choosing a nonzero root whose multiplicity is `2m`
in the `P` edge polynomial and applying

\[
 y\longmapsto y+\lambda x^{-1/5},\qquad \lambda\ne0.
\]

Put `X=x^(1/5)`.  The bracket changes by the chain rule:

\[
 [P,Q]_{X,y}=5X^4[P,Q]_{x,y},
\]

so scalar normalization makes the transformed bracket exactly `X^4`.
The two consecutive forced edges in the integral `(X,y)` lattice are

| edge | `P` endpoints | `Q` endpoints | weight |
| --- | --- | --- | --- |
| translated type II | `(75,60)` to `(21,6)` | `(125,100)` to `(35,10)` | `(1,-1)` |
| final type I | `(21,6)` to `(4,1)` | `(35,10)` to `(1,0)` | `(5,-17)` |

Every displayed endpoint is an actual vertex, hence its coefficient is
nonzero.  These facts are consequences of the regular-corner endpoints; they
do **not** say that the four displayed points are the complete polygons.

## 2. Complete normalization of the terminal edge

Write `s=X^17 y^5`.  Torus rescaling of `X,y,P,Q`, preserving the normalized
bracket, makes the terminal edge

\[
 P_{\rm I}=X^4y(1+s),\qquad
 Q_{\rm I}=-X\left(1+3s+\frac95s^2\right).
\]

Direct differentiation gives

\[
 [P_{\rm I},Q_{\rm I}]_{X,y}=X^4.
\]

Equivalently, before normalization, if the endpoint coefficients are
`a,b` for `P` and `c,d,e` for `Q`, the bracket equations are

\[
 -ac=1,\qquad 2ad-6bc=0,\qquad 5ae-3bd=0.
\]

All five coefficients are nonzero and the last three are uniquely fixed once
the two `P` vertices are normalized.  This block therefore has obstruction
rank zero: it is the mandatory terminal type-I bracket, not the new family
obstruction.

## 3. The forced common-power band

For

\[
 t=Xy,\qquad z=y^{-1},\qquad [t,z]_{X,y}=-1,
\]

the translated type-II edge is the common-power band

\[
 C_{\rm top}=t^7H(t)z^5,\quad \deg H=18,
\]

with nonzero constant and leading coefficients of `H`, and

\[
 P_{\rm top}=t^{21}H(t)^3z^{15},\qquad
 Q_{\rm top}=t^{35}H(t)^5z^{25}.
\]

Its bracket vanishes identically.  The normalized right-hand side is
`X^4=t^4z^4`; hence the missing descent spans 35 bracket layers, from the
formal top layer `15+25-1=39` to layer `4`.  No weighted-Wronskian or de Rham
class at the first nonzero layer can be computed without the intervening
band supports.

## 4. Why the older `(50,75)` calculation does not fill the gap

Section 5 of the 2014 polynomial-system paper treats the `j=0`, `(m,n)=(2,3)`
member.  It states two preliminary cases `gamma=2,3`, but explicitly says
that no proofs are supplied for that first reduction.  Its later coefficient
systems use `P=C^2`, `Q=C^3+...`; replacing those exponents by `3,5` is not a
consequence of the printed argument and would not establish exhaustive lower
supports.

## 5. Residual system and next obstruction

The JSON certificate contains three proof obligations:

1. prove support control after `y -> y+lambda*x^(-1/5)`;
2. determine, prove exhaustive, and normalize every `gamma` branch for
   `(m,n)=(3,5)` (the old values `2,3` may not simply be assumed);
3. enumerate every band through bracket layer 4 and only then compile the
   first nonzero weighted-Wronskian block.

This is the next genuinely new obstruction.  The terminal edge is already
solved and contributes no de Rham obstruction; the missing object is a
theorem-level exhaustive lower-boundary classification, not a larger
Gröbner-basis calculation.

## 6. Reproduction

```bash
python3 plane-jc/cas/test_f2_75_125_frontend.py
python3 plane-jc/cas/f2_75_125_frontend.py
```

The second command prints the complete machine-readable partial certificate.
