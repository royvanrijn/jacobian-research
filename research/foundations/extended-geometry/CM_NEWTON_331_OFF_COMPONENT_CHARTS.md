# Off-component boundary charts at 331

The [full-width strip note](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md) reduces
the conductor-18 boundary outside \(L_3\) to two points of \(L_2\).
Here the ten residual equations are the existing rows \(S_7,\ldots,S_{16}\)
for \(\deg B\leq5\) and \(D=-z^{-1}+c+ez\).

**Proposition.** Let an algebraic common zero of these residual equations
have weighted projective reduction at a place above 331 on \(L_2\)
outside \(L_3\). Normalize \(\ell=1\) and put
\(d=17w+1-8v\). Write \(\varepsilon=v-\alpha\), with \(\alpha\) the indicated
rational lift of the residue of \(v\). Let \(J(v)\) be the \(r\)-derivative
of the weighted homogenization of \(S_7\) along \(L_2\), in the
coordinates \((r,h,d,\tau)\) with \(\tau\) the square of the homogenizer.
Then \(y=d/J(v)^2\) is nonzero, and the residue of
\(E=\varepsilon/y\) must satisfy the following table.

| \(\alpha\) | Required residue of \(E\) |
|---|---:|
| \(2/17\) | \(275\) |
| \(5/56\) | \(170\) |

These two residues and the component \(L_3\) still require analysis.

## The two valuation charts

Work in the completion at the chosen place. Here \(J(v)\) is a unit.
Shift 7 solves uniquely for \(r\) with \(r\) in the maximal ideal.
Write \(h=J(v)^2 x\), \(d=J(v)^2 y\), and
\(\tau=J(v)^2 z\). The remaining rows are integral formal series in
\(x,y,z,\varepsilon\). Their lowest-degree forms, including those of
combinations cancelling the quadratic terms exactly, have radical
\[
       (z,x-\kappa y),\qquad
       \kappa=\begin{cases}17/18&\alpha=2/17,\\7/9&\alpha=5/56.\end{cases}
\]
The checker obtains these forms from cubic jets. Its quadratic pivot
matrix has the same rank over \(\mathbb Q\) and \(\mathbb F_{331}\), and its
adjugate gives integral multipliers for the exact cancellations.

Scale the coordinates by one of least valuation. If \(\varepsilon\)
has strictly smaller valuation than all of \(x,y,z\), use the first
chart below. Otherwise the displayed radical forces \(y\) to have least
valuation, \(x/y\) to reduce to \(\kappa\), and \(z/y\) to reduce to zero.
This gives the second chart. Since the original zero is affine,
\(\tau\ne0\).

### The chart with \(\varepsilon\) of strictly least valuation

Set \(s=\varepsilon\), \(x=sX\), \(y=sa\), and \(z=sb\). Now
\(X,a,b\) lie in the maximal ideal. Divide every remaining row by
\(s^2\). Shift 8 has a unit \(X\)-derivative and solves for \(X\)
as an integral series. Denote the other residuals by \(R_j\).

The exact \(s\)-jets through degree 3 of every \(R_j\) belong to
\[
       (a^5s^3,\ a^3bs^2,\ ab^2s,\ b^3s).
\]
The checker verifies this on numerators whose denominators are powers
of a polynomial that is a unit at \(a=b=0\). Independently, the exact
jets for the weights \(\operatorname{wt}(a,b)=(1,2)\) vanish below
weight 5. Thus every term in \(s\)-degree at least 4 is also divisible
by one of the four displayed monomials and has strictly greater
valuation than that monomial. All four leading coefficients of \(R_9\)
are units. The coefficient columns of \(ab^2s,b^3s\), in \(R_9,R_{10}\),
are independent modulo 331.

For \(a=0\) the \(b^3s\) term has strictly least valuation. For \(a\ne0\),
put \(A=\nu(a)\), \(B=\nu(b)\), and \(S=\nu(s)\). The four monomials have
a unique minimum except at \(B=A\) and \(B=S+2A\). At \(B=A\) the
independence of the last two coefficient columns excludes cancellation.

For \(B=S+2A\), put \(b=sa^2U\), with \(U\) a unit. In this calculation,
the unscaled free coordinates are the new \(a,b\) multiplied by the
square of the shift 8 pivot derivative at the origin in \(X,a,b\).
This derivative is a unit in \(\mathbb Z_{(331)}[[s]]\), so the change
preserves the valuation comparison and redefines \(U\) by a unit
factor. Write \(p_j(s)=[a^5]R_j(a,0)\) in the rescaled coordinates.
The exact \(s^3\)-coefficients after substitution are proportional as
rational functions of \(a,U\), with ratio given by the leading
coefficients of \(p_j(s)\). Therefore \(R_9\) is divisible by \(s^3\)
and \(G=p_9R_{10}-p_{10}R_9\) is divisible by \(s^7\). The two exact
normal-weight jets show that \(a^5\mid R_9\) and \(a^6\mid G\) after
substitution, so the powers of \(s\) and \(a\) can be divided out
simultaneously in the integral series ring. The reductions of the two
quotients at \(s=a=0\) are polynomials of degrees 2 and 3 in \(U\).
Their resultants are \(163\) at \(2/17\) and \(104\) at \(5/56\)
modulo 331. This excludes the balanced valuation as well.

### The chart with \(y\) of least valuation

Set \(s=y\), \(x=s(\kappa+X)\), \(z=sZ\), and
\(\varepsilon=sE\). Here \(X,Z,s\) lie in the maximal ideal and
\(E\) is integral. Divide the remaining rows by \(s^2\). The
shift 8 \(X\)-derivative at \(X=Z=s=0\) reduces to a nonzero multiple
of \(E-275\) at \(2/17\) and \(E-170\) at \(5/56\). Outside the
indicated residue it solves uniquely for \(X\). When \(s=0\) all
the rows vanish on an exact common linear graph in \(X,Z\).
Uniqueness identifies that graph with the implicit solution at \(s=0\).
Thus every remaining row on the implicit branch can be divided by \(s\).

The divided series start at degree 2 in \(Z,s\). The checker uses
their degree-2 forms and the degree-3 forms of exact pairwise
combinations cancelling degree 2. It checks the projective charts
\(Z=1\) and \(s=1\), using an inverse variable for the pivot derivative.
Both ideals are the unit ideal over \(\mathbb F_{331}\), which excludes
common zeros over its algebraic closure. Scaling by a coordinate
of least valuation then gives the required contradiction.

## Exact replay

Run with SageMath 10.9 from the repository root.

```sh
sage -python research/foundations/extended-geometry/verify_cm_newton_331_off_component_charts.py
```

The checker takes the ten rows from the
[reconstruction certificate checker](verify_cm_newton_full_width_reconstruction_certificates.py)
and uses the same weighted homogenization as the
[projective boundary checker](verify_cm_newton_full_width_projective_boundary.py).
All implicit jets and cancellations are computed over \(\mathbb Q\)
before reducing the leading forms modulo 331.
