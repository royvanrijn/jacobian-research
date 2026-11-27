# Unequal quartic parameters in the ternary quintic pencil

**Theorem QuinticKellerTernaryPencilUnequalParameters.** In the
[ternary quintic pencil equations](QUINTIC_KELLER_TERNARY_PENCIL_FIRST_EQUATIONS.md),
a Keller map with \(a=0\) necessarily has \(v_1=0\) and \(u_0=u_4\).
The proof uses exact coefficient operations over \(\mathbb Q\) and otherwise
works over \(\mathbb C\). It applies in characteristic zero by embedding
the finitely generated coefficient field.

## Degree-nine compatibility and the unequal case

Use the \(H_4\) parameters of the first-equations theorem with \(a=0\) and
\(u_1=u_3=v_3=0\). The linear operator for \(H_2\) in the degree-nine equation
is \(5z^3D(n\cdot H_2)\), with \(D,n\) as in that theorem. Its inhomogeneous
column is the sum of \(\det JH_4\) and the six determinants using one column
from each of \(JH_5,JH_4,JH_3\). The
[exact checker](verify_quintic_ternary_pencil_unequal_parameters.py)
constructs the full affine space of cubic solutions from the degree-ten
system, then takes the left nullspace of the degree-nine operator.
The resulting equations include, up to nonzero rational factors,
\[
 v_1^2(3v_0-v_4)=0,\qquad
 (u_0-u_4)u_0(3b+10u_0)=0,\qquad
 (u_0-u_4)(6bu_0-3bu_4+10u_0^2)=0 .
 \tag{1}
\]
For the last two equations it suffices to set \(v_1=u_5=v_2=v_5=0\) in
the checker.

The source shear \((x,y,z)\mapsto(x+\lambda z,y-2\lambda x-\lambda^2z,z)\)
preserves \(p,q\) and the leading vector. It preserves \(v_1\), adds
\(-2\lambda v_1\) to \(v_0\), and adds \(\lambda v_1\) to \(v_4\).
Apply the first equation of (1) before and after the shear with
\(\lambda=1\). The difference gives \(v_1=0\).
Source translations normalize \(u_5=v_2=v_5=0\). Indeed the \(z\)
translation adds a multiple of \(U=(5/2)z^2,\ V=q/2+yz\),
while the \(x,y\) translations add multiples of \(V=2xz,\ V=z^2\),
respectively, and leave \(U\) unchanged.

Suppose \(u_0-u_4\ne0\). The last two equations of (1) give \(b=0\)
and then \(u_0=0\). Now \(u_4=k\ne0\). The same shear changes
\(u_2\) by \(2\lambda(u_0-u_4)\), so we may arrange \(u_2=0\).
Renormalize \(u_5,v_2,v_5\) by translations as above. Neither
\(u_0-u_4\) nor \(u_2\) changes under these translations.

## A second weighted leading part

Give \(x,y,z\) the weights \(1,2,0\), respectively. In the normalized
unequal case all three components have weight at most four.
With \(s=q/y=x^2/y+z\), their weight-four parts are \(y^2\Phi(s,z)\), where
\[
 \begin{aligned}
 \Phi_1&=k^2z+A,\\
 \Phi_2&=(kz+f)s+khz+B,\\
 \Phi_3&=(z+2v)s^2+(2hz+g)s+h^2z+C .
 \end{aligned}
 \tag{2}
\]
Here \(h=v_4-v_0,\ v=v_0\), the constants \(A,B,C\) are the coefficients
of \(y^2\) in \(H_2\), and \(f,g\) are the coefficients of \(x^2y\)
in the second and third cubic components. The checker verifies this
weighted shape from the full cubic solution space. Terms of smaller
ordinary degree cannot affect (2) except for the three indicated
coefficients of \(H_2\).

The highest weighted determinant coefficient is zero, which in the
chart \((x,y,z)=(tX,t^2,z)\) is equivalent to
\(\det(\Phi_s,\Phi_z,\Phi)=0\), since \(s=X^2+z\) there and \(ds/dX=2X\)
is nonzero in the function field. The coefficients of \(s^2z,s^2\)
first give \(A=2kf-2k^2v\) and \(2k(f-2kv)^2=0\).
Thus \(f=2kv,\ A=2k^2v\). The coefficients of \(sz,z\) then give
\(g=2B/k,\ C=2Bh/k-2vh^2\). The constant coefficient is now
\(2k(B-2kvh)^2\), so \(B=2kvh\).
Consequently (2) is \((z+2v)N(k,s+h)\), with \(N(u,v)=(u^2,uv,v^2)\).
Translate \(z\) to \(z'=z+2v\). The weight-four leading vector is
\[
               \mathcal H=z' N(k y,\ Q),\qquad
               Q=x^2+y(z'+\eta),\quad \eta=h-2v .
 \tag{3}
\]

## First integrals of the new weighted leading vector

Suppose \(R\ne0\) is weighted homogeneous of positive weight \(j\le7\) and
\(\det J(\mathcal H_1,\mathcal H_2,R)=0\).
In the chart \((x,y,z')=(tX,t^2,z')\), use \(r=X^2+z'+\eta\) and
\(w=t^4z'\). The first two leading coordinates are \(k^2w,krw\).
The wedge-product argument of the
[conic-pencil exclusion](QUINTIC_KELLER_CONIC_PENCIL_FIRST_INTEGRALS.md)
applied to \((t^jR(X,1,z'))^4/w^j\) gives
\[
               R(X,1,z')^4/(z')^j\in\mathbb C(r).
\]
In fact \(w_t\ne0\), and \(\mathbb C(X,z')=\mathbb C(r)(X)\), so a binary
rational function algebraic over \(\mathbb C(r)\) lies in \(\mathbb C(r)\).
The polynomial \(z'=r-\eta-X^2\) has simple zeros over an algebraic
closure of \(\mathbb C(r)\). At either such zero the valuation equation
gives \(4\mid j\), hence \(j=4\). Applying the same argument directly
to \(t^4R(X,1,z')/w\) yields \(R=z'y^2 a(Q/y)\)
for \(a\in\mathbb C(T)\). The polynomials \(Q,y\) are coprime.
A finite pole of \(a\) would force \(Q-cy\mid z'\), and polynomial
degree greater than two would force \(y\mid z'\). Both are impossible,
so \(R\) lies in the span of (3).

Finally \(G=F_2^2-F_1F_3\) has weight at most seven and
\(\det J(F_1,F_2,G)=-\kappa F_1\), with \(\kappa=\det JF\ne0\).
A positive leading weight \(j\) of \(G\) gives a determinant coefficient
of weight \(j+5\) from its leading term and \(\mathcal H_1,\mathcal H_2\).
This coefficient is zero since the right side has weight four.
The first-integral result therefore allows subtraction of a linear
combination of \(F\) to leave \(G'\) of weight zero, meaning a polynomial
in \(z'\) alone. If \(G\) already has weight zero, take \(G'=G\).
In both cases the determinant with \(G'\) is \(-\kappa F_1\) plus
at most a constant. But for nonconstant \(G'\) it also equals
\[
                       (G')'(z')\,\{F_1,F_2\}_{x,y},
\]
which has weight five, since the leading bracket from (3) is
\(-4(z')^2 k^3xy^2\). If \(G'\) is constant this determinant vanishes.
Both possibilities contradict the weight-four right side. Thus
\(u_0-u_4\ne0\) is impossible, proving the theorem.
