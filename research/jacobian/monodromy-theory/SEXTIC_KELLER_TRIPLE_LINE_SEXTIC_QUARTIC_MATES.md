# The sextic-quartic free-degree pair in the sextic triple-line shear

**Lemma KellerParallelScalarSexticQuarticLeadingPowerExclusion.**
Let \(r=x+x^2v\). In characteristic zero there are no polynomial
Keller mates \(P,Q\) of free degrees \((6,4)\) in \(y\)
whose leading coefficients are \(p=\ell^3,\ q=H(r)\ell^2\),
up to nonzero scalars. Here \(\ell(x)\) is any polynomial with
\(1\le\deg\ell\le2\), and \(H\) is a nonzero scalar or
a nonzero scalar times \(r-c_0\).

## Leading powers and centered equations

Embed the finitely generated coefficient field in \(\mathbb C\).
Write \(k=\mathbb C(r)\), \(D=-x^2\partial_x|_r\).

Put \(w=\sqrt\ell\), let \(L=k(x,w)\), and change free variable
to \(Z=w(y+B_3/(4q))\), where \(B_3=[y^3]Q_{\rm orig}\).
Write \(Q_{\rm orig}=H Q\), so \(Q\) is monic centered
quartic in \(Z\), and the other mate \(P\) is monic sextic.
Their bracket in these variables is
\[
 D(P)Q_Z-P_ZD(Q)=\frac{\kappa}{Hw},\quad
 P=Z^6+aZ^5+bZ^4+eZ^3+fZ^2+gZ+i,\quad
 Q=Z^4+cZ^2+dZ+h .
 \tag{1}
\]
The coefficient of \(Z^8\) makes \(a\) invariant.
Let \(b_0,e_0,f_0,g_0,i_0\) be the invariant integration constants
from the next five coefficients. With \(K=h-c^2/4\), the
remaining two positive coefficients in (1) make these
expressions invariant:
\[
\begin{aligned}
L_2={}&-5ac^3/32-5acK/4-5ad^2/8
             -3c^2e_0/8-cg_0-2df_0-3dK-3e_0K,\\
L_1={}&5ac^2d/32-5adK/4+3cd^2/4+3cd e_0/4
             -dg_0-2f_0K-3K^2/2 .
\end{aligned}
\tag{2}
\]
The constant coefficient is
\[
\begin{aligned}
\frac{\kappa}{Hw}={}&
 \big(-15ac^3/64-5acK/8+5ad^2/16-3c^2e_0/8
                     -cg_0/2+df_0/2+3dK/4\big)Dc\\
 &+\big(5acd/16+3d^2/4+3de_0/4\big)Dd
   +\big(-15ac^2/32-5aK/4-3ce_0/4-g_0\big)DK .
\end{aligned}
\tag{3}
\]
These identities are replayed by the
[exact checker](verify_sextic_triple_line_sextic_quartic.py).

## The nonsquare case

Suppose \(\ell\) is not a square. Since \(\ell\) has a simple
root over \(\mathbb C\), the quadratic function-field extension
has constant field \(k\). Its involution sends \(Z\) to \(-Z\).
Consequently odd coefficients of \(P,Q\) are anti-invariant and
even coefficients invariant under the involution. The top
equations in (1) make \(a,e_0,g_0\) anti-invariant, so those
three \(D\)-constants vanish. Set \(T=K+2f_0/3\).
Equation (2) now makes \(dT\) invariant under \(D\), and \(dT\)
is anti-invariant under the involution, so \(dT=0\).
If \(d=0\) then (3) vanishes. Thus \(T=0\), \(cd^2\in k\),
and \(D(d^3/4)=\kappa/(Hw)\).

At \(x=0\), if \(\ell(0)\ne0\) the differential
\(dx/(x^2 w)\) has a double pole at each point above zero.
A primitive would have a simple pole, incompatible with
it being a scalar times \(d^3\). Hence \(\ell\) has a simple
zero at \(x=0\), and can be written \(x(\alpha x+\beta)\)
with \(\beta\ne0\). Put \(t=w/x\), so
\(t^2=\alpha+\beta/x\) and \(L=k(t)\).
The differential becomes
\[
                      \frac{dx}{x^2w}
                         =-\frac{2}{\beta^2}(t^2-\alpha)\,dt .
\]
Its primitive is proportional to \(t^3/3-\alpha t\)
up to an additive invariant. Since it must be \(d^3\)
up to a nonzero scalar, \(d\) is linear in \(t\).
The vanishing quadratic coefficient forces zero constant
term in that linear form. The linear coefficient then
gives \(\alpha=0\). Rescale to \(w^2=x\) and write
\(d=\delta/w\). Equation (3) gives \(\delta^3\) equal
to a nonzero scalar times \(1/H\). In \(k\) this forces
\(H\) constant, since its only other possible form has
a simple zero. Thus \(\delta\) is a nonzero scalar.

With \(S=Z^2+c/2\), the centered expressions in this
case satisfy, by the top equations of (1),
\[
       Q=S^2+dZ-2f_0/3,\qquad
       P-b_0Q=S^3+\tfrac32 dZS+3d^2/8+i_0 .
\]
Since \(cd^2\in k\), write \(S=x(Y^2+\lambda(r))\) with
\(Y=y+B_3/(4q)\). Here \(f_0,b_0,i_0,\lambda\) are
polynomial in \(r\): the centered coefficients are in
\(\mathbb C[r,w,w^{-1}]\), whose intersection with the
constant field \(k\) is \(\mathbb C[r]\).
At the original \(y=0\), if \(Y\) has a pole of order
\(j\ge1\) at \(x=0\), then \(Q\) has its unique dominant
pole of order \(4j-2>j\) from \(S^2\).
If \(Y\) is regular there, \(P-b_0Q\) instead has
the pole of \(3\delta^2/(8x)\). Both contradict
polynomiality.

## The square case

If \(\ell\) is a square, it is proportional to
\((x-a_*)^2\). The formal residue identity for the
fourth root of \(q\) requires a rational primitive of
\(dx/(x^2(x-a_*))\), hence \(a_*=0\).
Rescale so \(w=x\). All coefficients in (1) are then
Laurent in \(x\) over \(k\).

We use the weighted leading bracket calculation of the
[coprime pole lemma](SEXTIC_KELLER_TRIPLE_LINE_COPRIME_POLE_LEMMA.md),
whose initial pole bound uses no coprimeness.
If \(Q\) has no coefficient poles at either end of the
Laurent line, the top bracket equations give none for
\(P\) either, impossible. At an end with any pole of
\(Q\), let \(\omega\) be the maximum of pole order divided
by coefficient drop. Then \(\omega\ge1/4\).
After dividing \(D\) by the signed factor \(x\), the
threshold \(9\omega\) exceeds the pole order of the right
side, which is at most two. Therefore the normalized
leading polynomials satisfy \(P_*^4=Q_*^6\).
Centering gives \(Q_*=(Z^2+b_*)^2\) and its nonzero
offset gives \(b_*\ne0\). At this end \(c\) has pole
order \(u=2\omega>0\). If \(v_*,z_*\) denote the pole
orders of \(d,K\), respectively, then
\[
                         v_*<3u/2,\qquad z_*<2u .
\]
An exponent of a zero coefficient is treated as below
every finite pole order.

If \(a\ne0\), invariance of \(L_2\) gives
\(v_*+z_*=3u\) and \(dK\sim-5ac^3/96\) in leading parts,
since \(cd^2\) does not occur there and the remaining
terms have order below \(3u\). In particular
\(v_*>u,\ z_*>3u/2\). In \(L_1\) only \(cd^2\) and
\(K^2\) can balance at the highest order, so
\(u+2v_*=2z_*\). Thus \(v_*=5u/4,\ z_*=7u/4\).
After division by the same signed factor \(x\),
the maximal order in (3) is \(4u\). It comes
from the two coefficient terms \(-15ac^3/64+3dK/4\)
of \(Dc\), with sum of leading parts
\(-35au c^4/128\ne0\) after differentiation.
This contradicts the pole bound of two.

If \(a=0\) and \(e_0\ne0\), \(L_2\) similarly gives
\(v_*+z_*=2u\), \(dK\sim-e_0c^2/8\).
Both orders are positive by the strict bounds above.
The leading balance in \(L_1\) gives \(u+2v_*=2z_*\),
so \(v_*=3u/4,\ z_*=5u/4\). Now (3), after the common
factor is removed, has dominant order \(3u\)
from \((-3c^2e_0/8+3dK/4)Dc\), with nonzero leading
part \(-15e_0u c^3/32\). This also contradicts the bound.

If \(a=e_0=0\), put \(T=K+2f_0/3\). The invariant
expressions and the bracket reduce, up to invariant
constants, to
\[
  -g_0c-3dT,\qquad 3cd^2/4-g_0d-3T^2/2,\qquad
 \frac{\kappa}{Hx}=(-g_0c/2+3dT/4)Dc+3d^2Dd/4-g_0DT .
\]
If \(g_0\ne0\), let \(z_*\) now denote the order of \(T\).
The first invariant gives \(v_*+z_*=u\) and
\(dT\sim-g_0c/3\). If \(v_*\le0\), \(T^2\)
dominates the second invariant, while if \(z_*\le0\),
\(cd^2\) dominates it. Hence both are positive and
the balance is \(u+2v_*=2z_*\), which gives
\(v_*=u/4,\ z_*=3u/4\).
The leading bracket now has order \(2u\) after
division by the common factor, with leading part
\(-3g_0u c^2/4\). The integer pole order \(v_*\ge1\)
makes \(2u>2\), again impossible.

Finally let \(a=e_0=g_0=0\).
Then \(dT=A_0\) and \(cd^2-2T^2=B_0\) are invariant,
and (3) integrates to
\[
              \frac34 A_0 c+\frac14 d^3
                         =\frac{\kappa}{2H x^2}+U_0,
              \qquad U_0\in k .
\]
If \(A_0=0\) the left side is a cube times a scalar
and cannot have a pole of order two.
Otherwise \(d\ne0\) and substitution of
\(c=B_0/d^2+2 A_0^2/d^4\) shows its poles occur
only where \(d\) has a pole or zero. Those pole orders
are respectively multiples of three and of four,
again impossible at \(x=0\). This finishes the
square case and the lemma.

## Triple-line corollary

**Theorem SexticTripleLineShearSexticQuarticExclusion.**
A three-variable Keller map of ordinary total degree at most six
in characteristic zero with a quartic component \(g(x)+x^2(z+xy)\)
cannot have free-variable degrees \((6,4)\) after the shear
to \(r=x+x^2v\).

Use the [triple-line support and leading equation](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md#leading-coefficients-and-equal-degrees).
Since \(p\) is univariate in \(x\), it contains no full fibre
factor. After removing any full fibre factor from \(q\),
unique factorization in \(p^4/q^6\in k\) gives scalar
multiples of \(\ell^3,\ell^2\) for the remaining coefficients,
with \(\ell\) univariate and nonconstant.
The weight bounds give \(1\le\deg\ell\le2\) and at most
one full fibre factor in \(q\). This is precisely the lemma.

Together with the
[cubic-mate exclusion](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md),
[quartic-quintic exclusion](SEXTIC_KELLER_TRIPLE_LINE_QUARTIC_QUINTIC_MATES.md),
and [sextic-quintic exclusion](SEXTIC_KELLER_TRIPLE_LINE_COPRIME_POLE_LEMMA.md#sextic-application),
this rules out the quartic triple-line scalar \(g(x)+x^2(z+xy)\)
in every three-variable Keller map of ordinary degree at most six.
