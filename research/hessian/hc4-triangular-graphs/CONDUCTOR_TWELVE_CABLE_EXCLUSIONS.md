# The three cable exclusions in conductor degree twelve

**Theorem (Cable sequences in conductor degree twelve).** Let \(F=(P,Q)\) be an actual complex polynomial plane Keller map with
\(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\). Assume \(P\) is monic of positive
degree in \(x\) and the exact coordinate-order conductor has global
\(x\)-degree twelve. The reduced one-place characteristic sequence of
the generic coordinate order cannot be any of
\[
(6,4,9),\qquad (9,6,4),\qquad (10,4,5).
\]
These are the three sequences with a nontrivial degree-lowering
approximate root. The result concerns the actual polynomial generators
and retains every specialized source fibre. It makes no claim about the
remaining coprime generator sequences.

## 1. Common source and target identities

Write \(K=\mathbb C(y)\). The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives a birational generic order of defect six in \(K[x]\).
At every finite parameter in the hypothetical noninvertible case its
specialized conductor is a nonzero nonconstant polynomial. As in the
[cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md),
at least one generic matched pair has a finite distinct pair of limits
at each finite parameter, after extending to the endpoint branches.
Indeed a generic conductor root approaches a finite specialized root,
and its mate stays finite by monicity of the original \(P\). A pair
with equal limits would force \(P_x=Q_x=0\) there by their divided
differences. At each finite parameter at least one physical gap is
therefore a unit in the endpoint valuation.

Use the reduced generators of the exact singular generic order, obtained
by the polynomial target reductions described in the cubic–quintic proof.
For each sequence below the
[one-place approximate-root theorem of Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf)
gives a characteristic approximate root in these actual generators. The
standing hypotheses are as checked in the cubic–quintic proof: the
normalization is the source line and there is a single place at infinity;
the monic relation has the required bounds on its coefficient degrees.

Let the reduced, normalized generators be \(U,V\), of degrees \(m,n\)
in a longitudinal variable \(u=x-h(y)\) with \(h\in K\), obtained by a
polynomial target change over \(K\) of determinant \(\lambda\in K^*\).
Choose \(\rho^{m+n}=\lambda\), put \(u=\rho z\), and scale the two
generators by \(\rho^{-m},\rho^{-n}\). This gives a target determinant
one in an algebraic extension of \(K\), so
\[
 U_z V_y-U_y V_z=J\rho+\partial_z H(U,V)
\tag{1}
\]
for a polynomial \(H\) in the target coordinates over that extension.
The parameter derivative of the target change, composed with its inverse,
is a polynomial field of divergence zero and hence Hamiltonian.
The original map in \((z,y)\) has Jacobian \(J\rho\), giving (1).
A physical collision gap in \(x\) equals \(\rho\) times the gap in \(z\).

We will also use a boundedness consequence of monicity. Suppose two
polynomial functions of \(z\) in the normalized order, of coprime degrees,
have coefficients regular at a parameter place and nonzero leading
coefficients there. If either endpoint of a matched pair escaped to
infinity, so would the other, with the same valuation. The ratio of their
leading terms would be a root of unity for both degrees and hence one.
The divided difference of either polynomial would then have a unique
nonzero leading term, which is impossible. Thus the matched endpoints
are bounded in this situation.

## 2. Sequence \((6,4,9)\)

The approximate root is monic cubic in the degree-four generator.
Normalize \(U,V\) to be monic of degrees four and six, respectively.
Terms with exponents zero, one and two in \(U\) have distinct weights
modulo six for weights \((4,6)\). Consequently terms of weight greater
than twelve cannot cancel in the approximate root, and the weight-twelve
terms are \(U^3-V^2\). The possible product \(UV\) disappears upon
completing the square in \(V\); lower terms then have \(z\)-degree at
most eight. Thus \(\deg(V^2-U^3)=9\).

After depressing the source quartic, write \(U=z^4+az^2+bz+c\).
The polynomial part of its monic formal \(3/2\) power at infinity is
\[
 [U^{3/2}]_{\ge0}=\tfrac32 U L-\tfrac12 L^3+\tfrac38 b^2,
 \qquad L=z^2+a/2 .
\]
The difference of \(V\) from this polynomial part has degree three,
with nonzero leading coefficient \(p\). Translating \(U\) by a constant
in the source variable increases the \(z^2\) coefficient of its
polynomial \(3/2\) power by \(3/2\) times that constant. With this
translation and a translation of \(V\) we can arrange
\[
 V=[U^{3/2}]_{\ge0}+p z^3+r z .
\tag{2}
\]
These are allowable parameter-dependent target changes. Apply the scale
of Section 1, here \(\rho^{10}=\lambda\), and use the same notation for
the resulting scaled coefficients. Equation (2) remains valid.

The left side of (1) has degree at most six, so the filtered order
with semigroup \(\langle4,6,9\rangle\) gives
\(H\in\operatorname{span}\{1,U,V\}\) after evaluation. The degree-six
coefficient gives \(p'=0\) and the degree-four coefficient gives \(r'=3pa'/4\).
Thus \(p\in\mathbb C^*\) and \(r=R+3pa/4\), \(R\in\mathbb C\).
Since \(p\) has scaling weight three and \(\rho^{10}\in K\),
\(\rho^3\in K\) as well and hence \(\rho\in K\).

Put \(\Delta=a^2-4c\). After comparing the degree-five and degree-three
terms, the coefficients of \(V,U\) in \(H\) are respectively \(0\) and
\(-3pb'/4\). Degrees two, one and zero give
\[
\begin{split}
 A&=(b+p)\Delta-pa^2/2-4Ra/3\ \in\mathbb C,\\
 B&=\Delta^2-8ab(b+p)+32Rb/3\ \in\mathbb C,\\
 G&=b^3/4+3pb^2/8-3Aa/16-5pa^3/32
        +3pa\Delta/16-3Ra^2/8+R\Delta/4,\qquad G'=J\rho .
\end{split}
\tag{3}
\]
The
[exact checker](verify_conductor_twelve_cable_exclusions.py)
verifies these coefficient comparisons and their first integrals
without division by a coefficient that may vanish.

For a distinct matched pair let \(s\) be its sum and \(g\) its gap.
The divided differences of (2), for \(s\ne0\), give
\[
 \begin{split}
 0&=2s^6+4as^4-(2b+4p)s^3+3\Delta s^2-(2pa+8R)s+4pb,\\
 g^2&=-s^2-2a-2b/s.
 \end{split}
\tag{4}
\]
An identically zero pair sum would give \(b=0\) by the quartic
divided difference. If \(b\) is constant, (3) makes \(a,\Delta\) constant:
for \(b+p\ne0\), substituting \(\Delta=(A+pa^2/2+4Ra/3)/(b+p)\)
in the second invariant gives a nonzero polynomial equation for \(a\)
with constant coefficients and nonzero quartic coefficient. At \(b=-p\),
the first invariant makes \(a\) constant and the second makes \(\Delta\)
constant. Either way \(G'=0\) is impossible. Thus \(b\) is nonconstant
and every generic pair has nonzero sum.

A pole of \(a\) at a place where \(b\) is bounded is impossible:
the first equation in (3) would make the pole order of \(\Delta\) at
least twice that of \(a\), and the second would make it at most half
that of \(a\). If \(b\) has a pole and \(a\) is bounded, then (3) gives
\(\Delta=O(b^{-1})\), \(a=O(b^{-1})\), \(G\sim b^3/4\).
At a finite place \(\rho\) consequently has a pole. By the Newton polygon
of (4), three sums satisfy \(s^3\sim b\), with \(g^2\sim-3s^2\);
the other three have \(s^3\sim 2p\), with \(g^2\sim-2b/s\).
All physical gaps diverge.

If \(b,a\) both have poles, the invariants give
\[
 a\sim \eta b^{4/3},\quad \eta^3=32/p^2,\qquad
 \Delta\sim(p\eta^2/2)b^{5/3},\qquad G\sim-5pa^3/32 .
\]
There are three Newton polygon groups in (4), with two sums each.
In the largest group \(s\sim\sigma b^{2/3}\), \(\sigma^2=-2\eta\).
Eliminating \(s^2+2a\) in the gap formula gives leading gap squared
\[
 b^{1/3}\left(-3/\sigma+3p\eta^2/(4\sigma^2)\right).
\]
The coefficient cannot vanish since it would give \(\sigma=p\eta^2/4\)
and hence \(\sigma^2=2\eta\). In the middle group
\(s^2\sim-3\Delta/(4a)\), so \(g^2\sim-2a\).
In the smallest group \(s\sim \xi b^{-1/3}\) with
\(\tfrac32(\eta\xi)^2-2\eta\xi+4=0\).
The leading coefficient of \(g^2/b^{4/3}\) is \(-2\eta-2/\xi\);
its vanishing would give \(\eta\xi=-1\), which is not a root.
Thus each gap diverges, and so does the physical gap at a finite
pole since \(G'=J\rho\). The unit-gap property excludes all these finite
poles of the coefficients.

Now \(a,b,\Delta,G\) are polynomials in \(y\), with \(\rho=G'/J\).
For a putative finite zero of \(\rho\), the normalized sums remain bounded
by their monic equation in (4). The equation divided by \(s\) shows that
\(b/s\) is bounded as well, since \(p\ne0\). Then all normalized gaps
are bounded, and every physical gap would tend to zero. Hence
\(\rho\) is a nonzero constant and \(G\) is affine nonconstant.
But \(b\) is a nonconstant polynomial. The same asymptotics at its pole
at infinity give \(\deg_y G=3\deg_y b\) or \(4\deg_y b\).
This contradiction excludes \((6,4,9)\).

## 3. Sequence \((9,6,4)\)

In this case use \(U,V\) monic of degrees six and nine. The degree-four
approximate root, monic cubic in \(U\), gives after completing the
square and translating \(U\)
\[
 W=V^2-(U^3+\alpha U+\beta),\qquad \deg W=4.
\tag{5}
\]
Indeed weights \((6,9)\) are distinct modulo nine for the exponents zero,
one and two in \(U\), and lower terms after square completion have the
indicated form. Depress \(U\) in the source and use the scaling of
Section 1, here \(\rho^{15}=\lambda\). Let \(S\) denote the monic
formal square root of \(U^3+\alpha U+\beta\) at \(z=\infty\).
Then \(V-S=O(z^{-5})\). In (1), the contribution of \(S\) to the
polynomial part is \([\alpha' U_zU^{-1/2}/2]_{\ge0}\).
The term \(\beta'\) has strictly negative degree, and the contribution
of \(V-S\) is constant, six times the derivative of its \(z^{-5}\)
coefficient. Thus the left side has degree at most two. The semigroup
\(\langle4,6,9\rangle\) implies that \(H\) is constant as an element of
the order. The degree-two coefficient forces \(\alpha'=0\).

Write \(U=L^2+M\), where \(L=z^3+kz+\ell\) is monic depressed cubic
and \(M=r z^2+s z+t\) has degree at most two. The expansion of \(S\)
through degree minus five gives
\[
 S-V=\frac{3\,\operatorname{rem}(M^2,L)/8+\alpha/2}{L}
               -\frac{M^3}{16 L^3}+O(z^{-7}),
\tag{6}
\]
where \(V=[S]_{\ge0}\); \(\operatorname{rem}\) denotes the polynomial
remainder. The coefficients in degrees \(-1,-2,-3,-4\) of (6) vanish.
The degree-two and degree-one coefficients of the remainder give
\(s^2+2rt-r^2k=0\) and \(2st-r^2\ell-2rsk=0\).
The degree \(-4\) comparison gives \(r^2s=0\).
If \(r=0\) then \(s=0\) and the degree \(-3\) comparison makes the
negative tail in (6) have no \(z^{-5}\) term, contrary to (5).
Therefore \(s=\ell=0\) and \(t=rk/2\). The remaining coefficients give
\[
 \alpha=r^2(2r-3k^2)/16,\qquad
 [z^{-5}](S-V)=kr^3/32,\qquad k r\ne0 .
\tag{7}
\]
In particular \(J\rho=-3(kr^3)'/16\). The explicit model is
\[
 L=z^3+kz,\quad M=r(z^2+k/2),\quad
 U=L^2+M,\quad V=L^3+\tfrac32 LM+\tfrac38r^2z.
\tag{8}
\]
Up to a constant the degree-four function in (5) is
\[
 W=-\tfrac{r^3}{64}\big(2k^3+8k^2z^2+4kr+4kz^4-rz^2\big).
\tag{9}
\]

Let \(X=kr\); then \(3X^2=2r^3-16\alpha\).
If \(\alpha\ne0\), it has weight twelve under the scale, so
\(\rho^{12},\rho^{15}\in K\) and \(\rho^3\in K\).
The weight-six product \(X\) is therefore rational in \(y\).
It cannot have a finite pole: such a pole forces \(r\sim3k^2/2\)
and both \(k,r\) polar. Rescale \(z\) locally by \(z=\sqrt{k}\,w\)
and scale \(U,V,W\) by \(k^{-3},k^{-9/2},k^{-9}\) respectively.
The limiting functions in \(w\) are (8) and (9) at \(k=1,r=3/2\),
with a possible constant in \(W\) immaterial to matching.
The orders of the coefficients of the nonconstant terms are
nonnegative in this scaling. The two different degrees four
and nine of \(W,V\), with nonzero limiting leading coefficients,
make every rescaled matched gap bounded. None of those gaps
can tend to zero because the limiting \(U,V,W\) have no common
critical point. In fact the zeros of \(W_w\) have \(w=0\) or
\(w^2=-13/16\); \(V_w(0)=63/32\) and at the latter zeros
\(U_w/(2w)=315/256\). Thus every normalized gap in \(z\) has
the magnitude of \(\sqrt{k}\), and the pole of \(\rho\) from
\((kr^3)'\) makes every physical gap diverge. This violates
the finite unit-gap property.

Consequently \(X\) is polynomial, and it is nonconstant since
\(J\rho\ne0\). Differentiating the invariant gives
\[
 (kr^3)'=(Xr^2)'=X'(7X^2+16\alpha)/(2r).
\]
There is a finite value of \(y\) with \(7X^2+16\alpha=0\).
Here \(X,r,k\) each have a finite nonzero value by the invariant,
so the degree-four coefficient of (9) is nonzero.
The normalized gaps are bounded by the coprime degrees four and nine,
yet \(\rho=0\). This again violates the finite unit-gap property.

If \(\alpha=0\), (7) gives \(r=3k^2/2\). Put \(\tau^2=k\)
and \(w_0=\rho\tau\). All nonzero normalized matched gaps are nonzero
constant multiples of \(\tau\), by scaling the fixed homogeneous
shape (8). The finite unit-gap condition therefore forces the rational
function \(w_0^2=\rho^2k\in K\) to have neither a finite zero nor a pole.
So \(w_0\) is a nonzero constant. Since \(kr^3\) is a nonzero constant
multiple of \(\tau^{14}\), equation \(J\rho=-3(kr^3)'/16\) makes
\(\tau^{15}\) an affine nonconstant polynomial in \(y\).
The inverse determinant \(\lambda^{-1}=w_0^{-15}\tau^{15}\) has
a finite zero \(y_0\).

The original shift \(h(y)\) is polynomial: otherwise every matched
pair would escape at a pole of \(h\), since all actual shifted endpoints
are fixed finite numbers under \(w_0\). The unscaled reduced order is
fixed in \(u=x-h(y)\); denote it by \(A_*\). The actual generators
are in \(\mathbb C[y]\otimes A_*\), by coefficientwise projection to
\(\mathbb C[u]/A_*\). Their specialized order has colength at most six by
the exact specialized conductor and is contained in \(A_*\) of colength
six, so it equals \(A_*\), in particular at \(y_0\).
The physical matched actions \(\int P\,dQ\) at these fixed shifted
endpoints are polynomial functions of \(y\). Generically they equal the
fixed canonical actions scaled by \(\lambda^{-1}\), so all vanish at \(y_0\).

This is impossible for generators of the specialized order. To see
this, scale the fixed shape to \(k=1,r=3/2\) and take \(W=V^2-U^3\).
The degree-at-most-six filtered space has basis \(1,W,U\), all even,
so a reduced pair of degrees four and six cannot generate \(A_*\).
The characteristic sequences of
\(\langle4,6,9\rangle\) with reduced initial degrees are \((6,4,9)\)
and \((9,6,4)\). Indeed the conductor formula and the prefix estimate
of the cubic–quintic proof give \(2d_2\le 11\) if \(h\ge3\).
The gcds then require \(d_2=4,d_3=2\). Freeness gives
\(2r_2\in\langle r_0,r_1\rangle\) and \(r_2\equiv2\bmod 4\);
since \(r_1\ge8\), this forces \(r_2\ge6\). The resulting prefix
already contributes \(d_3\mu_3\ge8+6-2=12\), before the strictly
positive final contribution. Thus \(h=2\) and the three generators
must be the three minimal generators. The ordering and gcd condition
give the two sequences stated.
After target reduction of the actual specialized pair, the degrees
must therefore be six and nine. The lower degree generator is
\(U_*=a_*U+b_*W+c_*\) with \(a_*\ne0\); the higher is
\(V_*=d_*V+\text{an even function}\) with \(d_*\ne0\), since
\(1,W,U,W^2,V\) is the basis through degree nine. On each reflected
pair \((-z,z)\) with \(V(z)=0\), the even part contributes zero to
\(\int U_*\,dV_*\). Write \(T=z^2\) and
\[
 g(T)=V(z)/z=T^4+3T^3+21T^2/4+35T/8+63/32.
\]
This polynomial has discriminant \(35721/4096\) and nonzero constant
coefficient. Modulo \(g\), the reflected actions divided by \(z\)
for \(U\,dV\) and \(W\,dV\) are respectively
\[
 \frac{189}{160},\qquad
 -\frac{1161}{36608}T^3-\frac{8883}{56320}T^2
 -\frac{4473}{18304}T-\frac{200151}{585728}.
\]
They give linearly independent vectors on the four roots of \(g\).
Hence the specialized actions \(d_*(a_*\int U\,dV+b_*\int W\,dV)\)
cannot all vanish. Target reduction has nonzero determinant and scales
all matched actions by that determinant, giving the required
contradiction for the physical actions. This excludes \((9,6,4)\).

## 4. Sequence \((10,4,5)\)

Here \(U,V\) are monic of degrees four and ten, and the degree-five
approximate root is monic of degree five in \(U\). In the terms with
lower exponents in \(U\), weights \(4j\) for \(j=0,\ldots,4\) are
distinct modulo ten. As before the maximal terms are \(U^5-V^2\).
Completing the square in \(V\) by a target polynomial in \(U\) and
translating \(U\) gives an approximate root
\[
 W=V^2-(U^5+\alpha U^3+\beta U^2+\gamma U+\delta),\qquad \deg W=5 .
\]
Depress the source quartic and apply the determinant-one scaling with
\(\rho^{14}=\lambda\). The square root \(S\) of the displayed polynomial
in \(U\) satisfies \(V-S=O(z^{-5})\).
Here \(U_z\) has degree three, so this error contributes no polynomial
term to (1). The left side of (1) is precisely
\[
 [\,U_z(\alpha' U^{1/2}+\beta' U^{-1/2})/2\,]_{\ge0}.
\]
The other coefficient derivatives and the lower square-root terms have
strictly negative degree. In particular the degree is at most five.
The semigroup \(\langle4,5\rangle\) forces
\(H\in\operatorname{span}\{1,U,W\}\) after evaluation. The degree-five
term forces \(\alpha'=0\). For the depressed monic quartic
\([U_z U^{-1/2}]_{\ge0}=4z\); hence the left side becomes \(2\beta'z\).
The degree-four and degree-three terms of (1) force the coefficients
of \(W\) and \(U\) in \(H\) to vanish. The constant term is now
\(J\rho=0\), impossible. This excludes \((10,4,5)\) over the generic
parameter field itself.

The companion
[exact checker](verify_conductor_twelve_cable_exclusions.py)
checks the parameter coefficients in the quartic–sextic case, the
collision and critical-point identities, and the fixed reflected actions. The implications from the
primary approximate-root theorem and the finite-fibre valuation
arguments are the written parts of the proof.
