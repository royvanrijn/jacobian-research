# Regularity at the cusp divisor in an aligned sextic profile

**Theorem KellerSquarefreeWeightedCuspExclusion.**
Let \(k\) have characteristic zero, equipped with positive integer
source weights \((w_1,w_2,w_3)\). Let \(W,H\) be coprime weighted
homogeneous forms of weights \(2,b\) with \(b\ge4\) and
\(10+b-\sum w_i>9\). Assume over the algebraic closure that
\(W\) is squarefree, \((W,H)\) is onto the affine plane, and
\(k(W,H)\) is relatively algebraically closed in the source
function field. Then there is no Keller completion with highest
weighted tuple \((W^2,W^3,H)\).

**Theorem KellerWeightedCuspNonpureQuinticExclusion.**
Use weights \((1,2,2)\) on \(x,y,z\), and put
\(W=x^2+y,\ h=y^3z^2\).
There is no characteristic-zero three-variable Keller map
whose highest weighted tuple is \((W^2,W^3,h)\).

**Theorem SexticKellerWeightedSquareCoordinate.**
Let \(F\) be a characteristic-zero three-variable Keller map.
Suppose a component has highest form \((x^2+\ell(y,z))^2\) in
weights \((1,2,2)\), with \(\ell\ne0\) linear. That component
has a polynomial coordinate completion with inverse of ordinary
degree at most four. If \(F\) has ordinary degree at most
twenty-five, then it is invertible.

**Corollary SexticKeller654AlignedLinearPowerInversion.**
Let a characteristic-zero three-variable Keller map of
component degrees \((6,5,4)\) have highest forms
\((x^6,H,x^4)\) in a linear source and scalar frame.
Suppose \(H(0,y,z)\ne0\) and \(H/x^5\) is a closed projective
parameter over the algebraic closure. Then the map is
invertible. This handles the aligned linear-power alternative
of the irreducible-quintic reduction.

**Corollary SexticKeller654IrreducibleQuinticInversion.**
A characteristic-zero three-variable Keller map of ordinary
component degrees \((6,5,4)\) is invertible if its quintic highest
form is irreducible and coprime to its quartic highest form
over the algebraic closure.

The [irreducible-quintic reduction](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#an-irreducible-quintic-with-a-possibly-nonreduced-quartic)
leaves either the aligned linear powers of the preceding
corollary or the ordinary tuple \((W^3,H,W^2)\) with \(W\)
quadratic and \(H\) the quintic. In the latter case it proves
that \((W,H)\) is onto and that the pair field is relatively
closed. Also \(W\) is squarefree: otherwise over the algebraic
closure it would be proportional to \(x^2\), and the closed
projective parameter \(H^2/W^5\) from the reduction would be
the square of a nonconstant projective function. Apply
KellerSquarefreeWeightedCuspExclusion with ordinary weights
and \(b=5\), giving \(10+b-\sum w_i=12\).
The linear-power branch is invertible over the algebraic
closure, and the inverse descends by uniqueness of formal
inversion.

## Polynomial constants at the nonpure quintic boundary

Work over an algebraic closure \(k\).
The rational kernel of \(D=\det J(W,h,\cdot)\) is \(k(W,h)\).
To see this, fix generic values \(W=s_0,\ h=b\). After
adjoining \(\sqrt b\), on \(y\ne0\) one can put
\(y=r^2,\ z=\sqrt b/r^3\), with inverse
\(r=\sqrt b/(yz)\). The generic curve is the nonempty open
part \(r\ne0\) of the smooth conic \(x^2+r^2=s_0\).
This proves geometric integrality of the generic fibre
and relative algebraic closedness of its constant field.
The partials \(dW,dh\) are independent generically, so the
Jacobian criterion identifies the indicated rational kernel.
The pair \((W,h)\) is onto: choose \(x\) with
\(y=s_0-x^2\ne0\), and solve \(y^3z^2=b\).
Thus the polynomial kernel is \(k[W,h]\).
Indeed a reduced denominator in \(k[W,h]\) would have
a zero off the numerator on the parameter plane and
surjectivity would contradict the cross-multiplied identity.

For comparison, set \(R_0=xyz\).
The direct identities
\[
 \det J(W,h,R_0)=2Wh,\qquad R_0^2+h=W y^2z^2
\]
give the rational primitive \(R_0/(2Wh)\) for \(D\).
At the generic point of \(W=0\), the coefficient of its
simple pole cannot be removed by an invariant in \(k(W,h)\).
The residue of \(R_0\) there squares to \(-h\), which is
not a square in \(k(h)\). A correction in \(k(W,h)\)
with a simple pole has residue in \(k(h)\), and one with
a higher pole retains that higher pole.

The conductor of \(k[W^2,W^3,h]\subset k[W,h]\) is
\(W^2 k[W,h]\), by the powers of \(W\) that occur in the
cusp ring. The trace-dual comparison uses the finite free
extension \(B=k[W^2,h]\subset B\oplus W^3B\).
Its trace matrix in the basis \(1,W^3\) is
\(\operatorname{diag}(2,2W^6)\), hence its trace-dual
generates \(B[W^{-3}]=k[W,W^{-1},h]\).
This is the étale open in the
[trace-dual reconstruction theorem](../../hessian/hc4-triangular-graphs/CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md#1-full-reconstruction-theorem-hc4myku0).
The argument below works at the removed divisor and
uses regularity of the actual weighted rescaling there.

## Divisibility of each square correction

Work over the algebraic closure. Write \(F=(P,Q,T)\).
The rational kernel of \(D=\det J(W,H,\cdot)\) is \(k(W,H)\),
by the Jacobian criterion and relative closedness. The
polynomial kernel is \(k[W,H]\) by surjectivity and the
same reduced-denominator argument as above.

Use a reference
\[
                         (V^2,\varphi(V)+\lambda T,T)
\]
with \(V\) of highest weighted form \(W\) and \(\varphi\)
monic cubic. Start with \(V=W,\ \varphi(s)=s^3,\ \lambda=0\).
The reference determinant is zero. At the first weight loss
\(r\ge1\) from the reference, write \(E_P,E_Q\) for the
error terms of weights \(4-r,6-r\).
The rescaled determinant has order
\(N=10+b-\sum w_i>9\). At loss \(r\le6\), its
linear error equation is
\[
              D\big(W(2E_Q-3W E_P)\big)=0 .
\]
Here the first two errors are contracted with \(2W,3W^2\),
the derivatives of the leading reference functions of \(W\).
Since \(D(W)=0\), the polynomial \(2E_Q-3W E_P\)
belongs to \(k[W,H]\). Its weight \(6-r\le5\) and the
bound \(b\ge4\) allow a polynomial in \(W\) and possibly
a scalar multiple of \(H\).
Absorb half of this difference into \(\varphi(V)+\lambda T\).
The first errors then satisfy
\[
                            E_Q=\tfrac32 W E_P .
 \tag{1}
\]

Assume first that the \(P\)-error has positive weight, so
\(1\le r\le3\). Rescale with \(\epsilon\) using the
source and component weights. Fix \(T_\epsilon=S\), an
indeterminate, and use \(s=V_\epsilon\) as reference coordinate.
There is a regular formal chart at each geometric generic
component of \(W=0,\ H=S\). In fact \(W\) is reduced, and
\(H\) restricts to a nonconstant function on every
irreducible component of \(W=0\) by coprimality and homogeneity.
Thus \(dW,dH\) are independent at the generic points of the
geometric generic levels there, by separability in
characteristic zero. Choose a linear form to complete these
two differentials to coordinates at the selected generic
point. After fixing the third coordinate as well as
\(T_\epsilon=S\), the equations for the source variables
have invertible derivative in the residue field at that point.
Formal implicit solution gives regular functions of
\(s,\epsilon\) over the residue field.

On this chart write
\[
       P_\epsilon=s^2+U,\qquad
       Q_\epsilon=c(s,\epsilon)+V_* .
\]
Here \(c(s,\epsilon)\) is the second reference component
\(\varphi(V)+\lambda T\) rescaled on the fixed level.
It has \(c_0=s^3\), quadratic coefficient of order at
least two in \(\epsilon\), and linear coefficient of
order at least four. The latter orders are the degree
losses from weight six of \(s^2\) and \(s\).
Both errors start at order \(r\).
By (1) their leading coefficients satisfy
\((V_*)_r=(3/2)s U_r\).

The generic normalization gives a series
\(L(s,\epsilon)\) with scalar Laurent coefficients in
\(k(S)((s))\) satisfying on all the charts
\[
       c+V_*=\sum_{m\ge0}\frac{U^m}{m!}\delta^m L
                       \pmod{\epsilon^N},
             \qquad \delta=\frac{1}{2s}\partial_s.
 \tag{2}
\]
To obtain it, use the generic formal first integrals
\(P_\epsilon^{1/2}\), with leading \(W\), and \(T_\epsilon\).
At each unmatched \(\epsilon\)-order below \(N\), the
coefficient of \(Q_\epsilon\) after subtracting functions
of the two integrals is in the rational kernel \(k(W,H)\)
by the scaled determinant equation. Lift that coefficient
by evaluating the rational function at the formal integrals.
In each chart \(P_\epsilon^{1/2}=(s^2+U)^{1/2}\) with leading
\(s\), so the substitution gives (2). The same scalar
coefficients apply to every geometric component because
they were constructed in the generic field before
specializing at \(s=0\).

For \(i\le r\), \(L_i=c_i\), where subscripts denote
coefficients of \(\epsilon\). At order \(r\) this follows
from (1), since the change in evaluating \(c_0=s^3\)
is \(\delta(c_0)U_r=(3/2)sU_r\).
For every \(i<2r\), \(L_i\) is regular at \(s=0\):
any term involving \(U\) uses at most a single derivative
of \(L_j=c_j\) with \(j<r\), and these terms have no
linear coefficient in \(s\). The order \(2r\)
coefficient of (2) then forces at most a simple pole
of \(L_{2r}\) and gives its scalar coefficient
\[
                  \rho=[s^{-1}]L_{2r}=-\tfrac38\,u^2,
                      \qquad u=U_r(0).
\]
Here the only other possible pole at order \(2r\) is from
\(\tfrac12 U_r^2\delta^2(s^3)=3U_r^2/(8s)\), since
\(L_r=c_r\) also has no linear coefficient for \(r\le3\).

At order \(3r<N\), the terms that can have a pole of order
three use \(\delta L_{2r}\) and \(\delta^3(s^3)\).
In the terms with two Taylor derivatives the index of
\(L\) is at most \(r\), so they have at worst a simple pole
because the reference linear coefficients vanish through
\(r\). Terms with a single derivative and a lower index
also have at worst a simple pole. The cubic pole
coefficient therefore gives
\[
                 [s^{-3}]L_{3r}-\rho u/2-u^3/16
                     =[s^{-3}]L_{3r}+u^3/8=0 .
\]
Both \(u^2\) and \(u^3\) are thus scalars in \(k(S)\),
the same on all charts. Either \(u=0\) on every component
or \(u=u^3/u^2\) is the same element of \(k(S)\) on each.

It follows that \(E_P=f(H)\) on every component of
\(W=0\) for a single rational \(f\).
The restriction of \(H\) to each component is onto:
that component is a weighted cone, so scaling a point
where \(H\ne0\) by the weights reaches every nonzero value,
and the origin gives zero. Hence \(f\) has no denominator
by cross-multiplication on the component, and is a
polynomial. The weighted degree \(4-r\) is positive and
less than \(b\), so \(f=0\).
Reducedness now gives \(W\mid E_P\) in the polynomial ring.

Consequently the adjustment
\[
                          V\longmapsto V+\frac{E_P}{2W}
\]
is polynomial, regular across the cusp divisor. Its
change in the reference at the error weight is exactly
\(E_P\) and the right side of (1), and its added weight
is \(2-r<2\). Once no positive-weight \(P\)-error remains,
absorb a possible constant \(P\)-error and absorb the
remaining \(Q\)-terms by the same kernel calculation
with \(E_P=0\). The weighted loss strictly increases at
each step, so the map factors through the polynomial
pair \((V,T)\), contradicting the Keller condition.

For \(W=x^2+y,\ H=y^3z^2\), the hypotheses were checked
above and \(N=15\). In this particular case the cubic
pole comparison is unnecessary: at \(W=0,\ H=S\) one
can set \(y=-x^2,\ z=\gamma/x^3,\ \gamma^2=-S\), with
\(x\) free. A correction of weight \(j=4-r\in\{1,2,3\}\)
restricts to a linear combination of \(x^j\) and
\(\gamma x^{j-5}\), with the latter absent for \(j=1\).
Constancy of its square already forces both coefficients
to vanish. The general argument also explains why a
difference of signs between sheets cannot sustain a
nonconstant correction.

## Coordinate completion for any weighted square

Choose a linear frame with \(\ell=y\) and put
\(X=x,\ Y=x^2+y,\ Z=z\), a weighted polynomial
coordinate change. The component \(P\) in the weighted-square
theorem must have the form
\[
       P=Y^2+a(X)Y+b(X)+c(X)Z,\qquad
       \deg a,\deg c\le1,\quad \deg b\le3.
\]
Indeed all other terms involving \(Y,Z\) would have
weight at least four. Complete the square with
\(Y_1=Y+a(X)/2\) and \(b_0=b-a^2/4\).
The row \(dP\) cannot vanish anywhere over the
algebraic closure, by the Keller condition.

If \(c=c(X)\) is nonconstant, take its root \(X_0\).
At \(X_0,Y_1=0\) the remaining derivative
\(b_0'(X_0)+c'(X_0)Z\) can be made zero by solving
for \(Z\), giving a critical point. Thus \(c\) is
constant. If \(c\ne0\), \((P,X,Y_1)\) completes \(P\)
and its inverse has degree at most three, from
\[
          x=X,\qquad y=Y_1-a(X)/2-X^2,\qquad
          z=(P-Y_1^2-b_0(X))/c .
\]
If \(c=0\), then \(b_0'(X)\) has no zero, as any
zero at \(Y_1=0\) would be critical. Hence \(b_0\)
is linear with nonzero slope. Use
\((P,Y_1,Z)\) with inverse
\[
           x=b_0^{-1}(P-Y_1^2),\qquad
           y=Y_1-a(x)/2-x^2,\qquad z=Z .
\]
Its ordinary degree is at most four.
The calculation uses only the weight-four highest form and
the absence of critical points of \(P\). For a Keller map
of ordinary degree \(d\le25\), the resulting inverse completion
degree \(e\) satisfies \(de\le25\cdot4=100\), so the
[bounded coordinate completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
gives invertibility.
The aligned-linear-power corollary follows from the
[first-weight theorem](SEXTIC_KELLER_ALIGNED_LINEAR_POWERS.md),
which gives precisely such a weighted square for the
quartic component after its indicated target normalization.
