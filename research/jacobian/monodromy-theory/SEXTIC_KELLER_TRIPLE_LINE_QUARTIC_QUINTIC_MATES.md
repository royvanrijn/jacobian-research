# The quartic-quintic free-degree pair in the sextic triple-line shear

**Theorem SexticTripleLineShearQuarticQuinticExclusion.** A
three-variable Keller map of ordinary total degree at most six
in characteristic zero with a quartic component \(g(x)+x^2(z+xy)\)
cannot have free-variable degrees \((5,4)\) after the shear
to \(r=x+x^2v\).

Use the [triple-line support and cofactor bounds](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md#leading-coefficients-and-equal-degrees)
after embedding the coefficient field in \(\mathbb C\).
Put \(k=\mathbb C(r)\), \(D=-x^2\partial_x|_r\).
For leading free-variable coefficients \(p,q\), the leading
equation gives \(p^4/q^5\in k\). Since \(p\) is univariate
in \(x\) and nonconstant, it has no full fibre factor.
After stripping a possible full fibre factor from \(q\),
factorization gives scalar multiples of \(\ell^5,\ell^4\)
with \(\ell\) univariate linear by the weight bound of six.
A full fibre factor times \(\ell^4\) in \(q\) would exceed
that bound. The formal residue argument in
[Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1)
on the leading fourth root of \(q\) requires that
\(dx/(x^2\ell)\) have a rational primitive over \(k\).
A nonzero root of \(\ell\) would give a simple pole.
Thus after rescaling \(p=x^5,q=x^4\).

Set \(Z=xy+B_3/(4x^3)\) for \(B_3=[y^3]Q\).
The bracket over \(\mathbb C[r,x,x^{-1}]\), with
\(D\) holding \(Z\) fixed, is
\[
 D(P)Q_Z-P_ZD(Q)=\kappa/x,\quad
 P=Z^5+aZ^4+bZ^3+eZ^2+fZ+g,\quad Q=Z^4+cZ^2+dZ+h .
\]
Comparison at powers seven through three gives \(a\in k\)
and
\[
 b=5c/4+b_0,\quad e=ac+5d/4+e_0,\quad
 f=5c^2/32+3b_0c/4+ad+5h/4+f_0,\quad
 g=5cd/16+e_0c/2+3b_0d/4+ah+g_0
\]
with \(b_0,e_0,f_0,g_0\in k\). The coefficients of \(Z^2,Z\)
make the following two Laurent polynomials invariant:
\[
\begin{aligned}
 L_2 &=5c^3/32-5ch/4-5d^2/8
              +3b_0c^2/8-3b_0h-f_0c-2e_0d,\\
 L_1 &=15c^2d/32-5dh/4+e_0c^2/2-2e_0h+3b_0cd/4-f_0d .
\end{aligned}
\]
The constant bracket coefficient is
\[
 \frac{\kappa}{x}
 = (5d^2/16+e_0d/2)Dc
   +(3b_0d/4+5cd/16)Dd
   -(3b_0c/4+5c^2/32+f_0+5h/4)Dh .
 \tag{1}
\]
The [exact checker](verify_sextic_triple_line_quartic_quintic.py)
replays these identities and the dominant coefficients below.

Consider poles at zero or infinity over \(k\).
If \(d\) has a pole and \(c\) has none, invariance of
\(L_1\) makes \(h\) regular there since \(dh\) would
otherwise dominate. Then \(d^2\) dominates \(L_2\),
a contradiction. If \(c,d\) have pole orders \(u,v>0\)
at the same end, \(L_1\) forces \(h\) to have order \(2u\)
with leading part \(3c^2/8\). In \(L_2\) the leading
\(c^3,ch\) terms combine with coefficient \(-5/16\)
on \(c^3\), so \(2v=3u\) and \(d^2\) has leading part
\(-c^3/2\). To compare (1), divide the derivative terms
by their common signed factor \(x\) for that end.
The leading derivatives of \(c,d,h\) then become
\(u c,v d,2u h\), respectively.
At order \(4u\), the sum of dominant terms is
\(-55u c^4/64\) in leading parts. The right side
after the common factor has pole order at most two,
whereas \(4u>2\). This too is impossible.
Thus \(d\) is invariant.

If \(c\) has a pole of order \(u>0\), then
\(L_2\) forces \(h\) of order \(2u\) with leading
part \(c^2/8\). In (1) the dominant order after
the same common factor is \(4u\), with nonzero
leading part \(-5u c^4/64\). Thus \(c\) is invariant
as well.

Now (1) integrates to a quadratic polynomial of
\(h\) equal to \(\kappa/(2x^2)+H(r)\), with its
quadratic coefficient equal to \(-5/8\).
Completing the square and comparing Laurent ends
gives \(h=h_0(r)+\eta/x\), with
\(\eta^2=-4\kappa/5\ne0\).
As the coefficients of \(Q\) in \(Z\) lie in
\(\mathbb C[r,x,x^{-1}]\), the coefficients
\(c,d,h_0\) are polynomial in \(r\), and \(\eta\)
is a nonzero scalar. Setting the original \(y=0\)
gives \(Z=B_3/(4x^3)\). If \(Z\) is regular at
\(x=0\), the simple pole \(\eta/x\) persists.
If \(Z\) has a pole, its fourth power dominates
the other terms of \(Q\) there. Each possibility
contradicts polynomiality of \(Q\).
