# The cable sequence \((6,4,11)\) in conductor degree fourteen

**Theorem (The \((6,4,11)\) conductor-fourteen exclusion).** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\). If its exact
coordinate-order conductor has global \(x\)-degree fourteen, the
generic order has no reduced characteristic sequence \((6,4,11)\).
The proof uses its actual degree-four and degree-six generators and
allows drops in the conductor degree on the specialized source lines.

## 1. Parameter equations for the actual generators

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\in\mathbb C^*\).
The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the finite birational order, its exact specialized conductor
and its collision support at every finite parameter. In the
hypothetical noninvertible case this conductor is a nonzero nonconstant
polynomial on each line.

As in the
[degree-twelve cable proof](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md),
reduce the actual target generators over \(K\) by elementary polynomial
automorphisms. The
[Assi–García-Sánchez approximate-root theorem](sources/assi-garcia-sanchez-1407.0490v1.pdf),
Propositions 2 and 4, gives the degree-eleven approximation for the
reduced sequence \((6,4,11)\). The standing hypotheses for this
one-place theorem are verified in the
[cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md):
the normalization is the source line and the reduced monic relation
satisfies the requisite coefficient bounds. Normalize the generators to
be monic of degrees four and six. Their leading cancellation is
\(V^2-U^3\) of degree eleven. Depress the quartic by translating the
source coordinate by \(h(y)\in K\). The polynomial part of
\(U^{3/2}\) at source infinity is
\[
 V_{\mathrm{base}}=\frac32 U D-\frac12 D^3+\frac38 b^2,\qquad
 U=z^4+az^2+bz+c,\quad D=z^2+a/2.
\]
The leading term of \(V-V_{\mathrm{base}}\) is \(p z^5\), with
\(p\ne0\). Subtract a multiple of \(U\) from \(V\) to remove the
degree-four correction. Translating \(U\) changes the
degree-two coefficient of \(V_{\mathrm{base}}\) by \(3/2\)
times that translation. This absorbs the degree-two correction,
and translating \(V\) removes the constant correction.

As in the degree-twelve proof, if the target determinant is
\(\lambda\in K^*\), choose \(\rho^{10}=\lambda\), scale the source
by \(x=h(y)+\rho z\), and scale the two target generators with
weights four and six. The target determinant is then one. In
the normalized coordinates we have
\[
 V=V_{\mathrm{base}}+p z^5+r z^3+t z,\qquad
 U_zV_y-U_yV_z=J\rho+\partial_z H(U,V).
\tag{1}
\]
Here \(H\) is polynomial in the target coordinates because the
parameter derivative of the target automorphism is a polynomial
vector field with divergence zero. The left side has degree at most
eight in \(z\), so \(H(U,V)\) lies in the filtered span of
\(1,U,V,U^2\) for semigroup \(\langle4,6,11\rangle\).
Comparing the degree-eight and degree-seven terms gives \(p'=0\)
and a zero coefficient of \(U^2\) in \(H\). The unscaled \(p\)
has weight one, so \(\rho\in K\) and all the normalized coefficients
are rational in \(y\).

The next coefficient comparisons give
\[
 \begin{split}
 r&=R+5pa/4,\\
 t&=T+5pa^2/32+3Ra/4+5pc/4-25p^2b/24,
 \end{split}
\]
with \(R,T\in\mathbb C\). The coefficients of \(V,U\) in
\(H\) are \(-5pb'/6\) and \(5p(ab)'/16-3Rb'/4\),
respectively. Write \(\Delta=a^2-4c\) and define
\[
\begin{split}
 A={}&3\Delta R/4+5\Delta ap/16-3Ra^2/8-Ta-5a^3p/32
       +15b^2p/8+b(3\Delta/4+5Rp/2+25ap^2/24),\\
 B={}&-3\Delta^2/32+3ab^2/4
       +b(-5\Delta p/16+3Ra/4-T+5a^2p/32),\\
 G={}&-5\Delta^2p/128+\Delta T/4+5\Delta a^2p/128
       -Ra^3/16-Ta^2/8-15a^4p/512\\
    &\quad +b^3/4+b^2(3R/8-5ap/16-125p^3/288)\\
    &\quad +b(-3\Delta a/16-25\Delta p^2/96+5Tp/6+25a^2p^2/192).
\end{split}
\]
After removing \(\partial_zH\) the coefficients of degrees two,
one and zero are \(A',B',G'+aA'/4\). Thus
\[
 A,B\in\mathbb C,\qquad G'=J\rho.
\tag{2}
\]
The
[exact checker](verify_conductor_fourteen_cable_exclusion.py)
verifies all these coefficient comparisons and identities.

## 2. Collision sums and the two local pole regimes

If \(b\) is identically zero then \(B=-3\Delta^2/32\) is constant,
so \(\Delta\) is constant. The first invariant would then make \(a\)
constant since its cubic coefficient in \(a\) is \(-5p/32\ne0\).
This contradicts (2). Thus \(b\not\equiv0\). Over the generic
parameter a distinct matched pair cannot have sum \(s=0\),
by the divided difference of the quartic. For such a pair with
gap \(g\), the quartic gives the product \((s^3+as+b)/(2s)\), and
the sextic gives
\[
\begin{split}
 0={}&2s^7+2ps^6+4as^5+(-4R+3ap-2b)s^4+(3\Delta+8bp)s^3\\
     &+\big(\tfrac52\Delta p-2Ra-8T-\tfrac34a^2p
                         +\tfrac{25}3bp^2\big)s^2
       +(4Rb+abp)s-2pb^2,\\
 g^2={}&-s^2-2a-2b/s .
\end{split}
\tag{3}
\]
Only necessity of these equations for each actual pair is used.
The physical gap in \(x\) is \(\rho g\).

At a pole of one or more of \(a,b,\Delta\), use a local uniformizer and
let \(k>0\) be the maximum of their pole orders divided by \(2,3,4\),
respectively. Denote the leading coefficients in these three weighted
orders by \(a_*,b_*,\Delta_*\), allowing zero if the coefficient has
smaller pole order. The top parts of the constant invariants
give \(b_*\Delta_*=0\) and \(\Delta_*^2=8a_*b_*^2\).
Thus \(\Delta_*=0\), and exactly one of \(a_*,b_*\) is nonzero.
All root counts below are counts with multiplicity over a local
extension. The symbol \(\sim\) denotes equality of nonzero leading
terms in that valuation.

**The \(a\)-leading regime.**
Here \(a\) has pole order \(2k\). The term \(-5pa^3/32\) in \(A\)
can only be balanced by \(3b\Delta/4\), so
\(b\Delta\sim(5p/24)a^3\). Since \(\Delta\) has order less than \(4k\),
\(b\) has order greater than \(2k\). In \(B\) this makes
\(ab^2\) dominate the \(a^2b\) term and all lower corrections other
than \(\Delta^2\). Hence \(\Delta^2\sim8ab^2\). These two
balances give pole orders \(5k/2,7k/2\) for \(b,\Delta\),
and \(G\sim-35pa^4/512\).

The valuations of the coefficients of the sum equation in (3),
in order from degree zero through seven, are
\[
 -5k,\ -9k/2,\ -4k,\ -7k/2,\ -5k/2,\ -2k,\ 0,\ 0 .
\]
Its lower Newton polygon has segments from degrees zero to three,
three to five, and five to seven. The first four coefficients
are on the first segment, the endpoints of the other two segments
give their faces, and all face endpoints are nonzero. Thus the roots
split into the following three groups.

- Two largest sums have \(v(s)=-k\) and \(s^2\sim-2a\).
  Eliminate \(s^2+2a\) by the sum equation. The leading terms of
  \(g^2\) are \(-3b/s+3\Delta/(2s^2)\), of pole order \(3k/2\).
  Cancellation would give \(\Delta\sim2bs\), contradicting
  \(\Delta^2\sim8ab^2\) and \(s^2\sim-2a\).
- Two middle sums have \(v(s)=-3k/4\), with
  \(s^2\sim-3\Delta/(4a)\). Here \(g^2\sim-2a\), since
  both \(s^2\) and \(b/s\) have smaller pole order than \(a\).
- The three smallest sums have \(v(s)=-k/2\).
  Their face equation is
  \(3\Delta s^3-(3/4)a^2p\,s^2+abp\,s-2pb^2=0\)
  in leading coefficients. The two terms \(-2a-2b/s\) in the
  gap equation have order \(2k\), and cancellation would require
  \(s\sim-b/a\). The face evaluated there has leading value
  \(-35pb^2/8\ne0\), using the balance \(b\Delta\sim5pa^3/24\).
  Repeated roots of the face would still have the same valuations
  and cannot change this noncancellation.

Thus every normalized gap diverges.

**The \(b\)-leading regime.**
Here \(b\) has pole order \(3k\). The \(b^2\) term in \(A\)
forces \(\Delta\sim-(5p/2)b\).
If \(a\) were polar, \(ab^2\) in \(B\) would dominate every other
term since \(a\) has order less than \(2k\). Hence \(a\) is bounded,
and comparison at weight \(b^2\) in \(B\) even gives
\(a\to-25p^2/96\). We have \(G\sim b^3/4\).

The lower hull of the sum equation has segments from degrees zero to
four and four to seven. The coefficients at degrees zero, four and
seven have valuations \(-6k,-3k,0\), respectively. All others lie
strictly above these segments. Three sums have \(s^3\sim b\) and
give \(g^2\sim-3s^2\). The other four have \(s^4\sim-pb\) and give
\(g^2\sim-2b/s\). Again all normalized gaps diverge.

## 3. The finite-parameter and infinity contradictions

The finite distinct-pair property from the degree-twelve proof
applies unchanged. At any finite parameter a generic conductor root
limits to a finite root of the specialized conductor. Monicity of the
original \(P\) keeps its mate finite. The two limits must be distinct,
since otherwise their divided differences would force \(P_x=Q_x=0\)
at the limit. Hence at least one physical gap is a unit at every
finite parameter.

In either pole regime of Section 2, \(G\) has a pole and therefore
\(\rho=G'/J\) has a pole at a finite parameter, where the derivative
raises the pole order by one. Since every normalized gap diverges,
every physical gap diverges, contradicting the unit gap. Thus
\(a,b,\Delta,G\) are polynomials in \(y\), and \(\rho\) is polynomial.

If \(\rho\) vanishes at a finite parameter, the sums \(s\) stay
bounded by their polynomial equation in (3) with constant leading
coefficient. Dividing that equation by \(s^2\) expresses a quadratic
equation for \(b/s\) with leading coefficient \(-2p\), linear
coefficient \(4R+ap\), and bounded constant coefficient. Thus \(b/s\)
stays bounded even when \(s\) tends to zero. By (3) all normalized
gaps are bounded and hence all physical gaps tend to zero. This
again contradicts the unit gap. It follows that \(\rho\) is a
nonzero constant and \(G\) is affine nonconstant.

If any of \(a,b,\Delta\) is nonconstant, the same two coefficient
pole regimes apply at infinity, giving \(G\) degree \(8k\) or \(9k\).
These are greater than one, since respectively \(2k\) or \(3k\) is
a positive integral pole order. If all three coefficients are
constant then \(G\) is constant. Both possibilities contradict
\(G'=J\rho\ne0\) with \(\rho\) constant. This closes the
\((6,4,11)\) case.
