# Residue exclusions from the double-cover Cartier condition

Work over a characteristic-zero field with \(U=1+x\) and the regular
quartic-leading winding-one phase
\[
 f=\frac{U}{x}P(x,Ut),\qquad
 P(x,v)=v^3+d(x)v^2+b(x)v+c(x),\qquad
 \deg d\le3,\quad\deg b\le5,\quad\deg c\le7.
\]
For the normalized sphere functional \(\Lambda\), the
[double-cover Cartier obstruction](WINDING_ONE_CUBIC_QUARTIC_LEADING_CARTIER_TWO_FORM.md)
shows that an eventually vanishing positive moment tail forces
\[
 C_2\left(\frac{dx\wedge dv}{U Y}\right)=0,
 \qquad Y^2=xU P(x,v),                                           \tag{1}
\]
at all sufficiently large good primes of a number-field model. Here
\(C_2\) is Cartier on rational two-forms. Residues along two types of
divisor give additional characteristic-zero conditions.

**Proposition (Interior weighted residue).** Write
\(p=b-d^2/3\) and \(q=c-bd/3+2d^3/27\). Suppose that at a point
\(x_0\ne0,-1\), \(p\) vanishes to order at least four and \(q\) to
order at least six. Let \(p_4,q_6\) be the corresponding Taylor
coefficients. If \((p_4,q_6)\ne(0,0)\), infinitely many positive
sphere moments of \(f\) are nonzero.

**Proposition (Boundary residue).** Suppose there is a constant \(v_0\)
such that
\[
             P(x,v_0+Uw)=U^3 Q(U,w),                            \tag{2}
\]
where \(Q(U,w)\) has polynomial coefficients and is monic cubic
in \(w\). Put \(Q_0(w)=Q(0,w)\) and \(Q_1(w)=\partial_UQ(0,w)\).
Assume \(Q_0\) is squarefree. Eventual purity requires
\[
                 Q_1=3Q_0-w Q_0'+\alpha Q_0'                    \tag{3}
\]
for some scalar \(\alpha\). In particular, if
\(Q_1=\alpha Q_0'+\beta(3Q_0-w Q_0')\) with \(\beta\ne1\), then
infinitely many positive sphere moments are nonzero.

## Cartier and a Laurent residue

In the two charts below the completed local function field along the
specified divisor has the form \(L((s))\), where \(s\) is a transverse
parameter and \(L\) is the function field of the indicated residual
curve. At its generic point the residual square root is nonzero, so
the square-root equation lifts uniquely as a formal series after
choosing the residual root. This works in characteristic zero and
after the good reductions used below.

For a Laurent expansion \(ds\wedge \sum_j s^j \eta_j\) with rational
one-forms \(\eta_j\) over \(L\), the residue is \(\eta_{-1}\).
In characteristic \(p\), the formula for Cartier on Laurent terms
gives
\[
     \operatorname{Res}_s C_2(\omega)=C_1(\operatorname{Res}_s\omega).
                                                                    \tag{4}
\]
Indeed, the only input power of \(s\) that becomes \(s^{-1}\) is
\(s^{-1}\) itself. The formula uses \(s,w\) as a \(p\)-basis, with
\(w\) separating on the residual curve. It agrees with the same
formula in the rational function field inside the completion. This
argument also applies when the form has poles of order larger than
one. It uses the [Cartier formula in Katz, Theorem 7.2](https://web.math.princeton.edu/~nmk/old/nilpconn.pdf).

## The interior divisor

Put \(s=x-x_0\), \(v+d/3=s^2 w\), and \(Y=s^3 W\). Then
\[
 W^2=xU\left(w^3+\frac{p}{s^4}w+\frac{q}{s^6}\right),\qquad
 \omega=\frac{ds}{s}\wedge\frac{dw}{U W}.
\]
At \(s=0\) the residual differential is \(dw/(U_0 W_0)\), with
\(U_0=1+x_0\) and
\(W_0^2=x_0U_0(w^3+p_4w+q_6)\). For a squarefree residual cubic it
is a nonzero holomorphic differential on an elliptic curve.

Over a number field, pass to an extension defining the elliptic
curve's full \(3\)-torsion. Choose good absolute degree-one primes
of this extension, for example above rational primes split completely
in its normal closure. There are infinitely many by Chebotarev.
At such a prime \(p\ne3\) the elliptic trace \(a_p\) is \(2\) modulo
three because Frobenius is the identity on \(3\)-torsion. The Hasse
bound \(|a_p|\le2\sqrt p\) implies \(p\nmid a_p\) for \(p>4\).
On a cubic model \(y^2=g(w)\) the coefficient of \(w^{p-1}\) in
\(g(w)^{(p-1)/2}\) is \(a_p\) modulo \(p\), by summing the quadratic
character over \(w\in\mathbb F_p\). It is also the selected Cartier
coefficient of \(dw/y\), as seen by factoring out \((1/y)^p\).
Thus the holomorphic residue has nonzero Cartier at infinitely many
good primes. This contradicts (1) and (4).

If the residual cubic is singular, the assumption on \(p_4,q_6\)
excludes a triple root. Write it as \((w-r)^2(w+2r)\), with
\(r\ne0\). On its normalization, \(W_0/(w-r)\) is a square root
of \(x_0U_0(w+2r)\). The residual differential has two simple poles
over \(w=r\) with nonzero residues. These survive at good primes, so
the differential is not exact there. A one-form with zero Cartier
in a function field is exact by the same \(p\)-basis calculation as
for (1), or equivalently Cartier takes the \(p\)th root of its simple
residue. This again contradicts (1) and (4).

## The divisor at \(U=0\)

Under (2), put \(W=Y/U^2\). Then
\[
   W^2=(U-1)Q(U,w),\qquad
   \omega=\frac{dU\wedge dw}{U^2 W}.
\]
At \(U=0\), let \(W_0^2=-Q_0(w)\) and \(\eta_0=dw/W_0\). Taking
the linear term in \(1/W\) gives the residue
\[
        \eta_1=\frac12\left(1-\frac{Q_1}{Q_0}\right)\eta_0.      \tag{5}
\]
For example, if \(Q_1=\alpha Q_0'+\beta(3Q_0-wQ_0')\), then
\[
     \eta_1=\frac{1-\beta}{2}\eta_0
                 +\alpha\,d(1/W_0)-\beta\,d(w/W_0).
\]
If \(\beta\ne1\), the ordinary-prime argument just given shows that
Cartier is nonzero on (5) at infinitely many primes.

For general \(Q_1\), equations (1) and (4) instead force (5) to be
exact in the residual curve's function field at all sufficiently large
good primes. Apply [Bost, Theorem 2.9](https://numdam.org/item/10.1007/s10240-001-8191-3.pdf)
to the additive torsor \(dz=\eta_1\) over the smooth open curve
where \(\eta_1\) is regular. It has a flat connection and solvable
structure group \(\mathbb G_a\), and the exact reductions have
zero \(p\)-curvature. The connection therefore has a horizontal
section on a finite étale cover. Its field trace, divided by the
covering degree, gives a rational primitive of \(\eta_1\) on the
original curve.

The differential (5) changes sign under \(W_0\mapsto-W_0\). Its
primitive can be taken anti-invariant by subtracting its conjugate
and dividing by two, so it has the form \(A(w)/W_0\) with \(A\)
rational. The only possible poles of (5) are double poles at the
simple branch roots of \(Q_0\), and (5) is regular at infinity.
Thus its primitive has at most simple poles at the finite branch
points and no poles elsewhere. This forces \(A\) to be a polynomial
of degree at most one. Indeed a finite pole of \(A\) at a branch root
has order at least two upstairs, before the pole of \(1/W_0\) is
added. At infinity, \(w\) has pole order two and \(W_0\) order three.

Write \(A(w)=aw+b\). Comparing its derivative with (5) gives
\(2aQ_0-(aw+b)Q_0'=Q_0-Q_1\). Since \(Q_1\) has degree at most two,
the cubic coefficient forces \(a=-1\), which gives (3) with
\(\alpha=b\). Conversely (3) makes (5) the derivative of
\((\alpha-w)/W_0\), so (3) is exactly the condition obtained from
this residue.

For a degree-seven example in the excluded boundary subcase, take
\(v_0=0\) and \(Q(U,w)=w^3+2U^2w^2-w-2U^2\). The resulting
coefficients are \(d=2U^3\), \(b=-U^2\), and \(c=-2U^5\).
Regularity holds at \(x=0\), and \(D=U^3d\) has degree six,
giving exact sphere degree seven by the
[degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).
Here \(Q_0=w^3-w\) is squarefree and \(Q_1=0\), which violates (3).

All arguments may be made over a number-field specialization of the
finitely generated coefficient algebra. For the interior case retain
\(x_0U_0\ne0\) and one nonzero Taylor coefficient among \(p_4,q_6\).
For the boundary case retain squarefreeness and any proposed failure
of (3), which is detected by a nonzero minor between the coefficient
vectors \(Q_1-(3Q_0-wQ_0')\) and \(Q_0'\). The pure moment identities
are rational polynomial identities and persist under specialization.
The contradictions or necessary condition over a number field
therefore prove the statements over characteristic zero.
