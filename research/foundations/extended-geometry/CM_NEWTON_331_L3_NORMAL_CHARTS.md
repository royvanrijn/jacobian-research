# The nonendpoint stratum of the weighted boundary at 331

Let \(S_7,\ldots,S_{16}\) be the normalized conductor-18 residuals
in the coordinates \((r,\ell,h,v,w)\) of weights \((2,1,2,1,2)\).
The weighted reduction of a putative algebraic common zero lies on
the three curves of the [projective boundary calculation](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#the-projective-boundary-at-331).

**Proposition (the nonendpoint part of the third weighted component).**
At a place above \(331\), an algebraic common zero
cannot reduce to a point on
\(L_3=\{r=0,\ \ell=8v,\ h=17w\}\) where both \(v\) and \(w\) are
nonzero in the weighted residue coordinates.

Normalize by \(v=1\) in the weighted projective chart. Let \(t\)
be the homogenizer and write \(\tau=t^2\),
\(L=\ell-8\) and \(H=h-17w\). The elements \(L,H,r,\tau\)
are in the maximal ideal, \(\tau\ne0\), and \(w\) is a unit.
In the chart, \(S_7\) has \(r\)-derivative
\(-12829696w^2/194828780125125\) at \(L=H=r=\tau=0\),
a unit. Retain the exact solution jet for \(r\) through
ordinary total degree seven in \(L,H,\tau\).
Write the other rows on this graph. Their degrees zero
through two vanish exactly.

The jet comparisons used below follow from polynomial Taylor
expansion over the valuation ring. If two coordinate vectors have
the same reduction at a unit Jacobian, the minimum valuation of the
difference of their images equals that of the coordinate difference.
Apply this first to \(S_7\) with the degree-seven approximate solution.
If \(d=\min(\nu(L),\nu(H),\nu(\tau))\), the true \(r\) and its
approximation differ by valuation at least \(8d\). Hence the
residual chart expansions agree through degree seven for the
subsequent valuation tests.

## The chart where \(\tau\) has least valuation

Set \(L=\tau x,\ H=\tau y\). The residuals divided by \(\tau^3\)
have leading polynomials given by the cubic forms at \((x,y,1)\).
In the residue ring with \(w\) inverted their common zeros satisfy,
with \(K=wx\),
\[
 P_3(K)=K^3-136K^2-27K+75=0,\qquad
 y=98K^2+(8/w-101)K+62 .
\]
The checker obtains these inclusions by polynomial ideal
membership with an inverse variable for \(w\). It computes an
exact triangular lift of the cubic leading system over
\(\mathbb Q(w)\), verifies that all cubic rows vanish on this
lift, and checks that its coefficients are integral when \(w\)
is a unit. The reduction of its monic polynomial is \(P_3\),
which is square-free.

At an algebraic common zero, lift the required simple root of
\(P_3\) to the exact monic polynomial in the completed valued
field, adjoining the root if needed. This gives a leading pair
\((x_0,y_0)\) of the same reduction for which all cubic rows
vanish exactly. Rows 8 and 9 have Jacobian determinant whose
reduced norm in the cubic algebra is \(191w^6\).
In particular these are unit pivots at every leading pair.
Solving their order-\(\tau\) equations gives an approximate
pair with error of valuation at least \(2\nu(\tau)\) by the
Taylor comparison. The order-\(\tau\) coefficient of row 10
on the pair has reduced norm \(68\) in the cubic algebra,
so it is a unit and row 10 cannot vanish.

## The remaining normal chart

If \(\tau\) does not have least valuation, the projective initial
cubic equations are evaluated at \(\tau=0\). The checker verifies
over the residue ring with \(w\) inverted that their zeros have
\(H\ne0\) and \(L/H=1/8\) in this projectivized direction.
Thus put \(s=H,\ X=L/H-1/8,\ Z=\tau/H\).
The elements \(X,Z,s\) are in the maximal ideal.

Divide the residuals by \(s^3\). They start in degree two
in \(X,Z,s\); the exact quadratic forms are multiples of the
quadratic form of row 8. Cancel these quadratic forms in shifts
9 through 16 by their exact multiples of row 8. The coefficients
used in the cancellations are integral for unit \(w\).
The checker uses the quadratic form of row 8 and the cubic forms
of the cancelled rows. On the projective boundary \(s=0\), the
charts with \(X=1\) and \(Z=1\) both have unit ideal over the
residue ring with \(w\) inverted. Therefore \(s\) has least
valuation. Put \(x=X/s,\ z=Z/s\).

In this chart the initial equations imply the triangular system
\[
 \begin{aligned}
 P_6(z)&=z^6-15z^5-45z^4-110z^3+74z^2+124z-76=0,\\
 x&=-41z^5+74z^4+18z^3+140z^2-81z+32 .
 \end{aligned}
\]
The checker verifies the two inclusions by explicit ideal
membership over the residue ring with \(w\) inverted.
As in the first chart, it verifies the exact characteristic-zero
lift and every leading equation in the corresponding finite
algebra. Since \(P_6\) is square-free, a leading pair lifts to
a pair of the same residue satisfying these exact leading
equations. In the second chart, normalize row 8 by one
further factor \(s^2\) and the cancelled rows by \(s^3\).
The pivot determinant for row 8 and cancelled row 9 has
reduced norm \(245w^{12}\), a unit. The first-order correction
given by these pivots yields in cancelled row 10 an
order-\(s\) coefficient of reduced norm \(158w^6\).
This coefficient is a unit. Taylor comparison again makes
the higher terms have valuation at least \(2\nu(s)\),
contradicting vanishing.

This excludes the two normal charts over every weighted point
with \(vw\ne0\) on \(L_3\). In combination with the earlier
boundary calculations, only four reduced points remain
possible. Two are the points of \(L_2\) with \(v/\ell=2/17\)
or \(5/56\), with the additional transverse restrictions of the
off-component chart calculation. The other two are the endpoints
\(v=0\) and \(w=0\) of \(L_3\).

## Exact replay

    sage -python research/foundations/extended-geometry/verify_cm_newton_331_l3_normal_charts.py
