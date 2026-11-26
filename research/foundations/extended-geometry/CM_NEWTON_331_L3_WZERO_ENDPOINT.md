# The endpoint \(w=0\) of \(L_3\)

**Proposition (exclusion of the \(w=0\) endpoint).**
At an algebraic common zero of the seven normalized residuals \(S_7,\ldots,S_{13}\),
the weighted reduction cannot be the endpoint \(w=0\) of \(L_3\).

Use the chart \(v=1\) of the weighted compactification. Write
\((r,L,H,w,\tau)=(r,\ell-8,h-17w,w,T_{\mathrm{hom}}^2)\) for
the local coordinates, where \(T_{\mathrm{hom}}\) is the original
homogenizer. All five coordinates are in the maximal ideal at this
endpoint. The homogenizer is nonzero at an affine candidate. Let
\(d=\min(\nu(L),\nu(H),\nu(w),\nu(\tau))>0\).

## Two normalizations

For weights \((2,1,1,1,1)\) on \((r,L,H,w,\tau)\) there is an exact
integral residual combination with no terms of weight below eight.
Its weight-eight form is \(r^4\) modulo \(331\) and \(\tau\).
The checker forms its multipliers from the weight-eight truncated
coefficient matrix. It solves the \(238\)-pivot system to cancel
all lower coefficients exactly over \(\mathbb Q\), keeping the
weight-eight target modulo \(331\).

If \(0<\nu(r)<2d\), a monomial of weight at least eight with \(i\le4\)
powers of \(r\) has valuation at least \(8d-i(2d-\nu(r))\). This
is strictly greater than \(4\nu(r)\) except for \(r^4\) itself.
For \(i>4\) its valuation is again strictly greater. The combination
therefore forces \(\nu(r)\ge2d\). In these weights, use the lowest
forms of the scalar-echelon combinations of the actual rows,
with exact cancellation over \(\mathbb Q\) and all pivot scalars
units at \(331\). The forms have unit ideal on the \(\tau=1\)
residue chart. Consequently \(\nu(\tau)>d\), and the weight-eight
combination now forces \(\nu(r)>2d\). With \(r=\tau=0\) the lowest
forms imply \((L-H/8)^3=0\). Thus \(U=L-H/8\) satisfies
\(\nu(U)>d=\min(\nu(H),\nu(w))\).

In \((r,U,H,w,\tau)\), use weights \((3,2,1,1,2)\). Set
\(d'=\min(\nu(r)/3,\nu(U)/2,\nu(H),\nu(w),\nu(\tau)/2)\).
To exclude the charts \(r=1\) and \(\tau=1\) with this shorter
prefix, take monomial multiples of the rows through weight nine.
The truncated matrix is \(231\times342\), of rank \(62\) modulo
\(331\). Choose its pivot rows and columns modulo \(331\) and solve
the resulting square system over \(\mathbb Q\). All entries earlier
than each pivot vanish exactly, and the resulting lowest forms
generate the unit ideal on either chart modulo \(331\).
The scalar-echelon leading forms exclude \(U=1\) with
\(r=\tau=H=w=0\).
It follows that \(d'=d\), \(\nu(r)>3d\), \(\nu(\tau)>2d\),
and \(\nu(U)\ge2d\). The chart \(H=1,\ w=0\) is also excluded.
Hence \(d=\nu(w)\). Divide each scalar-echelon row by the power
of the actual \(w\) given by its lowest weight. Use the new
coordinates
\[
R=r/w^3,\quad X=U/w^2,\quad Y_0=H/w,\quad Z=\tau/w^2.
\]
Here \(R,Z\) are in the maximal ideal. The radical of the common
residue equations in these coordinates gives two directions
\[
    (Y_0,X)=(0,0),\qquad (-17,-289/512).
\]


## The L3-tangent direction

At \((Y_0,X)=(0,0)\), the first scaled row has unit \(R\)-pivot
\(7860239/288028122\). Solve through ordinary degree seven in
\((X,Y_0,Z)\) with \(w\) as a parameter. The remaining rows have
exact orders at least three. Their cubics at \(w=0\) force
\(X=0\) whenever \(Z=0\) in projectivized normal coordinates.
If \(Z\) is least, the chart \(Z=1\) of the initial cubics gives
a simple cubic lift, reducing to
\[
u^3+50u^2+56u+274 .
\]
The computed exact lift has coefficients integral when \(w\)
is in the maximal ideal and satisfies all initial equations.
Two pivots give the first correction, and the next residual
has obstruction norm \(324\bmod331\) after correction.

If \(Z\) is not least, set \(s=Y_0\), \(x=X/s\), \(z=Z/s\).
All \(x,z,s\) are small. Divide by \(s^3\). The resulting
rows start at degree two. Retain the first quadratic and
cancel the exact proportional quadratic portions in the other
rows. The new leading equations exclude \(s=0\) in both charts.
In the \(s=1\) chart the leading algebra has an exact simple
sextic lift, reducing to
\[
u^6+316u^5+286u^4+221u^3+74u^2+124u+255 .
\]
The two unit pivots and first correction now give obstruction
norm \(285\bmod331\).

In each test lift the selected simple residue root by Hensel's
lemma, and evaluate the exact leading solution at the actual
\(w\). The integral implicit-jet approximations at a unit
Jacobian differ from the true coordinates by valuation at
least the valuation of their error in the equations, by the
polynomial difference formula. Degree seven is sufficient
before the division by \(s^3\) and the final division by
\(s^2\) or \(s^3\), giving an error of order at least two
in the last normalized test. The unit norm then contradicts
vanishing.

## The L2-tangent direction

At \((Y_0,X)=(-17,-289/512)\), take the Hensel root \(a\) near
\(8\) of \(a^2-8a+17w=0\), and write \(w=a(8-a)/17\).
The centre along \(L_2\) is \(Y_0=-17,\ X=-289/(8a^2)\).
Use variables \(Y=Y_0+17=h/w\) and \(Z=\tau/w^2\), of weights
one and two, respectively. The first two scaled residuals
solve for \(R\) and \(X+289/(8a^2)\), using the unit Jacobian
whose reduction is \(255\) at \(a=8\).

The checker computes the implicit jets through weight six
over \(\mathbb Q(a)\), and verifies integrality at \(a\equiv8\).
The remaining rows vanish exactly below weight five. Their
weight-five forms are proportional over \(\mathbb Q(a)\).
One has reduced weight-five form
\(Y(37Z^2+165Y^2Z+142Y^4)\). Cancel weight five in the
selected next row using its exact multiple of this pivot.
The resulting weight-six form reduces to \(119Z^3\).

Write \(s=\min(\nu(Y),\nu(Z)/2)\). The implicit approximation
has error of valuation at least \(7s\), so these two exact
leading-weight tests apply after rescaling. If \(Y=0\) in
the scaled residue, the weight-six row excludes a nonzero
\(Z\). If \(Z=0\) there, the weight-five row excludes nonzero
\(Y\). If both are nonzero use \(Y=1\) and the resultant
\(325\bmod331\).

All valuations considered are at a place over \(331\). Unit
normalizations and ramified scalings can be performed after a
finite valued extension without changing any valuation
comparison. Both second-normalization directions are
excluded, eliminating the entire endpoint.

## Exact replay

    sage -python research/foundations/extended-geometry/verify_cm_newton_331_l3_wzero_endpoint.py
