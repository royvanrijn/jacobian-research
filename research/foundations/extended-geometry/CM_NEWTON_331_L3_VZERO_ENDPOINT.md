# The \(v=0\) endpoint of the weighted component \(L_3\)

**Proposition (exclusion of the \(v=0\) endpoint).**
A normalized algebraic common zero of \(S_7,\ldots,S_{13}\)
cannot reduce above \(331\) to the \(v=0\) endpoint of \(L_3\).

The [normal-chart calculation](CM_NEWTON_331_L3_NORMAL_CHARTS.md)
excludes the nonendpoint part of \(L_3\). Use the weighted chart \(w=1\), adjoining a square root of the
normalizing unit if necessary. Put \(\tau=t^2\), \(L=\ell-8v\) and
\(H=h-17\) in these normalized coordinates. At the endpoint
\(r,L,H,\tau,v\) have strictly positive valuation. Also \(\tau\ne0\)
at a generic-fibre point.

The exact derivative of the homogenized row 7 in \(r\) at
\(r=L=H=\tau=0\) is \(-12829696/194828780125125\).
Thus the unique implicit formal solution \(r=R(L,H,\tau,v)\)
near this axis has \(331\)-integral coefficients. Here \(v\)
can be treated as a parameter, since \(R(0,0,0,v)=0\)
and the derivative on the axis is the displayed constant.
For a proposed solution, the implicit series converges at the
maximal-ideal coordinates, and equals its \(r\)-coordinate by
unit-Jacobian uniqueness. The verifier substitutes its degree-eight
jet in \(L,H,\tau\) into the remaining rows. They begin at degree
three. In every row the omitted terms have valuation at least
\(9\min\{\nu(L),\nu(H),\nu(\tau)\}\).
The substitution preserves the parity under simultaneous sign
change of \(L,v\). In each even-index row, every monomial therefore
has an odd combined exponent in \(L,v\).

## The first chart when \(\tau\) is least

Set \(x=L/\tau,\ y=H/\tau\) and divide the residuals by \(\tau^3\).
At \(v=0\) in residue, row 9 forces
\(y^3-69y^2-57y+56=0\), while row 8 is \(x(124y^2+100y-101)\).
The coefficient of \(x\) is a unit at every allowed \(y\), so
\(x\) is in the maximal ideal. More precisely row 8, evaluated on
the exact implicit solution, is \(xU+vV\) where \(U\) is a unit
at these residues. This follows by parity, with the leading
coefficient of \(x\) as displayed. If \(v=0\) exactly, the same
equation gives \(L=0\), hence \(\ell=v=0\) in the original
affine chart. The exact [zero-slope identity through row 13](CM_NEWTON_SEVENTEEN_ROW_LOCUS.md#exact-coherent-blindness-through-shift-12)
excludes this slice. Thus \(v\ne0\) in this chart and \(x=vK\) with
\(K\) integral. Parity makes the extra division by \(v\) in the
even-index rows valid for the exact series. It preserves the
omitted-term bound in \(\tau\) after the substitution \(L=\tau vK\).

The leading normalized equations have a triangular simple cubic
in \(K,y\) over the residue field. The verifier checks an exact lift
with integral coefficients for small \(v\), the unit norm
\(-191\) for the Jacobian of rows 8 and 9, and the norm \(68\)
for the order-\(\tau\) obstruction in row 10 after their first
correction. Taylor comparison over the valuation ring now
excludes vanishing of all three rows.

## The other first chart

If \(\tau\) is not least, row 9 at \(\tau=0\) forces \(H\) out
of the least-valuation direction. Put \(s=L,\ X=H/s,\ Z=\tau/s\).
Both \(X\) and \(Z\) lie in the maximal ideal. Divide the residuals
by \(s^3\). The leading degrees in \(X,Z,v,s\) are two for the
even-index rows and three for the odd-index rows. At \(s=0\)
the radical of the reduced leading equations forces \(X=8v\)
and \(Z=0\), as checked by ideal membership of powers.

### When \(s\) is least among \(X,Z,v,s\)

Put \(X=sx,\ Z=sz,\ v=s\eta\). The parameter \(\eta\) is integral.
Divide the rows by the additional powers \(s^2\) or \(s^3\)
corresponding to their leading degrees. Their geometric leading
solutions in the residue field are
\[
 P_6(z)=z^6+33z^5+47z^4-113z^3-154z^2+124z+101=0,
 \qquad
 x=139z^5+74z^4+159z^3+148z^2-14z-165+8\eta .
\]
The verifier checks the inclusions over the polynomial ring in
\(\eta\), an exact triangular lift over \(\mathbb Q(\eta)\) with
\(331\)-integral polynomial coefficients, and exact vanishing of
all leading rows on that lift. The sextic is square-free in the
residue field. Thus the chosen residue root lifts at the given
\(\eta\) to an exact leading solution. The Jacobian of rows 8 and
9 has reduced norm \(26\). Solve for the first two coefficients
of the correction in \(s\), and substitute them in row 10.
Its coefficient at order \(s\) is exactly zero, and its
order-\(s^2\) coefficient has reduced norm \(182\).
All omitted terms have valuation strictly greater than
\(2\nu(s)\), even after the divisions by \(s^3\) and \(s^2\)
or \(s^3\). The nonzero order-\(s^2\) term excludes this chart.

### When \(s\) is not least among \(X,Z,v,s\)

The projective test above forces \(v\) to be least and to have
the indicated residue direction. Write
\[
 X=v(8+A),\qquad Z=vB,\qquad s=vT.
\]
Here \(A,B,T\) lie in the maximal ideal. For each row divide
out the corresponding power \(v^2\) or \(v^3\). As polynomials
in \(A,B,T\) with coefficients in the \(v\)-ring, the quadratic
forms are exact multiples of the row-8 form. Cancel them in the
other rows using integral multiples of row 8. The initial
quadratic and cancelled cubic forms at \(T=0\) have unit ideal
on either projective chart \(A=1\) or \(B=1\) over the residue
field. Hence \(T\) must be least. Set \(A=Tx,\ B=Tz\), and
divide row 8 and the other rows by \(T^2\) and \(T^3\)
respectively.

The same residue sextic and graph as above, now with
\(\eta=0\), give an exact leading lift. The verifier shows
that every positive power of \(T\) in these normalized
equations has a coefficient divisible by \(v^2\).
The Jacobian of rows 8 and 9 again has reduced norm \(26\).
The displacement from the exact leading pair consequently
has valuation at least \(\nu(T)+2\nu(v)\), by the unit-Jacobian
Taylor comparison. Divide the first correction by \(Tv^2\).
After eliminating it with rows 8 and 9, the normalized
row-10 obstruction at that scale has reduced norm \(182\).
Terms with higher powers of \(T\) have strictly greater
valuation. For the omitted terms from the original implicit
jet, the worst denominator is \(s^3 v^3 T^3\);
their valuation is at least \(3\nu(T)+3\nu(v)\),
which is again strictly greater than \(\nu(T)+2\nu(v)\).
Thus row 10 cannot vanish.

The three charts exhaust this endpoint over every residue
extension at \(331\). Together with the previous exclusions,
a putative characteristic-zero blind point can now reduce
only to the \(w=0\) endpoint of \(L_3\) or to the two specified
off-component points.

## Replay

    sage -python research/foundations/extended-geometry/verify_cm_newton_331_l3_vzero_endpoint.py
