# The quartic–decic sequence \((10,4,11)\) in conductor degree eighteen

**Theorem (The \((10,4,11)\) conductor-eighteen exclusion).** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\). If its exact
coordinate-order conductor has global \(x\)-degree eighteen, its
reduced generic characteristic sequence cannot be \((10,4,11)\).
The proof allows degree drops on individual source lines.

## 1. Laurent equations and the actual parameter invariants

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the exact finite birational generic and specialized orders.
Under the conductor-degree hypothesis the generic order has colength
nine. Its reduced generators of degrees four and ten have an
approximate root of degree eleven by
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf).
The one-place hypothesis and monic coefficient conditions for that
theorem are those verified in the
[cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).

Normalize the generators \(U,V\) to be monic. The approximate root
is monic of degree five in \(U\); the weights \(4i\), \(0\le i<5\),
are distinct modulo ten, so no terms above weight twenty can cancel.
The top terms are \(U^5-V^2\). Completing the square in \(V\)
and translating \(U\) by a constant in the parameter field gives
\[
 W=V^2-(U^5+\alpha U^3+\beta U^2+\gamma U+\delta),\qquad \deg W=11.
\tag{1}
\]
These are polynomial target changes, as in the
[preceding quartic–decic exclusion](CONDUCTOR_SIXTEEN_QUARTIC_DECIC_EXCLUSION.md).
Depress the quartic by a source translation \(h(y)\in K\).
For the determinant \(\lambda\in K^*\) of the target changes choose
\(\rho^{14}=\lambda\), write \(x=h(y)+\rho z\), and rescale the
generators by \(\rho^{-4},\rho^{-10}\), respectively. The resulting
target determinant is one. In the normalized frame write
\[
 U=z^4+az^2+bz+c,\qquad \Delta=a^2-4c .
\]

Let \(S\) be the monic formal square root of the polynomial in \(U\)
in (1), expanded at infinity. Since \(V+S\) has degree ten, the
Laurent difference is
\[
 V-S=pz+\sum_{i\ge1}e_i z^{-i},\qquad p\ne0 ,
\]
after a constant target translation of \(V\). Such a translation
may be made before scaling and does not affect the determinant.
In this expansion \(V=[S]_{\ge0}+pz\). The determinant-one parameter
identity is
\[
 U_z V_y-U_y V_z=J\rho+\partial_z H(U,V) .
\tag{2}
\]
As in the
[degree-twelve cable proof](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities),
the derivative of the target change composed with its inverse is
a polynomial target field of divergence zero, hence Hamiltonian.

The polynomial contribution from \(S\) to the left side of (2)
is
\([U_z(\alpha'U^{1/2}+\beta'U^{-1/2})/2]_{\ge0}\).
The term \(pz\) contributes \(p'zU_z-pU_y\).
The negative tail contributes
\(4e_1'z^2+4e_2'z+4e_3'+2ae_1'+a'e_1\).
The semigroup \(\langle4,10,11\rangle\) makes \(H(U,V)\) a linear
combination of \(1,U\), since the left side has degree at most five.
Its degree-five and degree-four coefficients force
\(\alpha'=p'=0\); the degree-three comparison then makes \(H\)
constant as an element of the order. Using
\([U_z U^{-1/2}]_{\ge0}=4z\), the remaining comparisons show
\[
 4e_1'-pa'=0,\qquad 4e_2'+2\beta'-pb'=0,\qquad
 (4e_3+ae_1+pa^2/8-pc)'=J\rho.
\]
Here and below primes denote derivatives in \(y\).

The binomial expansion
\(S=U^{5/2}+\alpha U^{1/2}/2+\beta U^{-1/2}/2+O(z^{-6})\)
gives
\[
 \begin{split}
 -e_1 &= b(15\Delta^2-40ab^2+64\alpha)/256,\\
 -e_2-\beta/2 &= 5(-\Delta^3+24\Delta ab^2-8b^4)/1024-\alpha\Delta/16,\\
 e_3 &= -ae_1/2-5\Delta b^3/128 .
 \end{split}
\]
Consequently, with \(\kappa=-e_1+pa/4\in\mathbb C\), there is a
constant \(C\) such that
\[
 \begin{split}
 a(40b^3-64p)&=b(15\Delta^2+64\alpha)-256\kappa,\\
 C&=-\Delta^3+24\Delta ab^2-8b^4-64\alpha\Delta/5+256pb/5,\\
 G&=a\kappa-3pa^2/8+p\Delta/4-5\Delta b^3/32,\qquad G'=J\rho .
 \end{split}
\tag{3}
\]
The
[exact checker](verify_conductor_eighteen_quartic_decic_exclusion.py)
verifies these coefficient identities and the collision equations
below. The unscaled coefficient of the correction \(pz\) is
\(\rho^9p\in K\), so \(\rho^9,\rho^{14}\in K\) and
\(\rho\in K\) by coprimality. Thus \(a,b,c,\Delta,G\) are rational
functions of \(y\).

## 2. The collision equations

If \(b\equiv0\), the first two equations of (3) force \(a,\Delta\)
constant and then \(G'=0\). We may therefore take \(b\ne0\) in \(K\).
Write \(s\) for the sum and \(g\) for the gap of a distinct matched
pair in the normalized \(z\)-coordinate. The quartic divided difference
makes \(s\ne0\) and gives
\[
 g^2=-s^2-2a-2b/s .
\tag{4}
\]
The decic divided difference, after substituting (4) and multiplying by
128, is the polynomial
\[
\begin{split}
 128p+s\big(&8s^8+32as^6-8bs^5+(20\Delta+32a^2)s^4-16abs^3\\
            & +(40\Delta a-8b^2)s^2-20\Delta b s
               +15\Delta^2-40ab^2+64\alpha\big).
\end{split}
\tag{5}
\]
It vanishes at each paired sum. In particular it has constant leading
coefficient and nonzero constant coefficient \(128p\).
Using the first parameter invariant and dividing by \(4s\) gives
the equivalent sum equation
\[
 \begin{split}
 F(s)={}&2s^4(s^2+2a)^2+5\Delta s^2(s^2+2a)-2bs^3(s^2+2a)\\
       & -2b^2s^2-5\Delta bs+(64\kappa-16pa)/b+32p/s=0 .
 \end{split}
\tag{6}
\]
Writing \(Y=-g^2=s^2+2a+2b/s\) in this equation gives
\[
 2s^4Y^2+(5\Delta s^2-10bs^3)Y
       +10b^2s^2-15\Delta bs+(64\kappa-16pa)/b+32p/s=0 .
\tag{7}
\]

At any finite parameter at least one generic matched pair has a finite
distinct pair of limits in the original \(x\)-coordinate. This is the
finite-unit-gap argument of the lower cable proofs. The specialized
conductor is nonzero and nonconstant, so a generic conductor branch
reaches a finite specialized root. Its mate stays finite by monicity
of \(P\), and equal limits would force \(P_x=Q_x=0\) there by divided
differences. Hence some physical gap \(\rho g\) is a unit at every finite
parameter, after extending to the endpoint branches.

## 3. The normalized coefficients have no finite poles

Consider first a finite pole of \(b\). The first parameter equation
shows that if \(\Delta\) were bounded then \(a=O(b^{-2})\), which
would make \(-8b^4\) dominate the second invariant. Therefore
\(\Delta\) has a pole, \(a\sim3\Delta^2/(8b^2)\), and
\(\Delta^3\sim b^4\). Here \(\sim\) denotes equality of nonzero leading
terms at the chosen place. The pole orders of \(b,\Delta,a\) are
\(3k,4k,2k\) for some \(k>0\). We have \(G\sim-5\Delta b^3/32\), so
\(\rho\) has pole order \(13k+1\). Under \(z=\sigma w\) with
\(\sigma=\Delta/b\), the quartic and decic scaled by
\(\sigma^{-4},\sigma^{-10}\) converge to
\[
 U_*=w^4+\tfrac38 w^2+w-\tfrac{55}{256},\qquad
 V_*=[U_*^{5/2}]_{\ge0}.
\]
The normalized correction \(pz\) vanishes in this limit. The two
limiting derivatives have resultant \(208000\ne0\), as checked in
the exact replay. A unit physical gap would give a scaled gap
in \(w\) tending to zero since \(\rho\sigma\) has a pole.
Such close scaled endpoints must be bounded by the monic quartic
divided difference: two unbounded endpoints at a vanishing gap
would have equal leading terms, giving a nonzero leading divided
difference. Bounded close endpoints would give a common critical point
of \(U_*,V_*\). Thus \(b\) has no finite pole.

If \(b\) tends to zero and \(\Delta\) has a pole, the first invariant
gives \(a=O(1+b\Delta^2)\). The term \(-\Delta^3\) in the second
invariant would then dominate, which is impossible. If both
\(b,\Delta\) are bounded and \(b\to0\), the first invariant makes \(a\)
bounded as well. It remains to examine poles with \(b\to b_0\ne0\).

**A pole of \(\Delta\) with \(b\to b_0\ne0\).**
The second invariant gives \(a\sim\Delta^2/(24b_0^2)\);
the first then forces \(b_0^3=-p/5\). Put \(n=-v(\Delta)>0\)
at the finite parameter. Thus \(v(a)=-2n\) and
\(G\sim-3pa^2/8\), so \(v(\rho)=-4n-1\).

Use the polynomial \(sF(s)\). In the list of coefficients indexed by
the powers of \(s\), the lower Newton polygon connects
\((0,0),(1,-2n),(5,-4n),(9,0)\). Its first face has the
leading terms \(32p-16pa s/b\), the second face has
degrees one, three and five, and the last has leading terms
\(2s^5(s^2+2a)^2\). It gives one sum of valuation \(2n\),
four of valuation \(n/2\) and four of valuation \(-n\).
The small sum satisfies \(s\sim2b/a\) and \(g^2\sim-3a\).
The middle sums give \(g^2\sim-2a\).
At the large sums, write \(s\sim q\Delta/b\), so \(q^2=-1/12\).
In (7) the coefficients of \(Y^2,Y,1\) have valuations
\(-4n,\ge-3n,-2n\), respectively. The last leading coefficient
relative to \(\Delta^2\) is
\[
 10q^2-15q-\frac{2p}{3b_0^3}=\frac52-15q\ne0 ,
\]
since \(q^2=-1/12\). Thus \(v(Y)=n\) in this group.
Even this shrinking normalized gap has a diverging physical gap because
\(v(\rho)=-4n-1\). All groups contradict the finite-unit-gap property.

**A pole of \(a\) with \(\Delta\) bounded and \(b\to b_0\ne0\).**
The second invariant makes \(\Delta a\) bounded, hence \(\Delta\to0\).
The first requires \(b_0^3=8p/5\). Put \(n=-v(a)>0\);
here \(G\sim-3pa^2/8\) and \(v(\rho)=-2n-1\).
The lower Newton polygon of \(sF(s)\) connects
\((0,0),(1,-n),(5,-2n),(9,0)\). Its first and last faces are again
\(32p-16pa s/b\) and \(2s^5(s^2+2a)^2\); its middle face has
degrees one and five. Thus one sum has valuation \(n\), four
have valuation \(n/4\), and four have valuation \(-n/2\).
The small and middle sums give \(g^2\sim-3a\) and \(-2a\),
respectively. The large sums satisfy \(s^2\sim-2a\). In (7)
the coefficients have valuations \(-2n,\ge-3n/2,-n\); the last
leading term is \((-20b_0^3-16p)a/b_0=-48pa/b_0\ne0\).
Thus \(v(Y)=n/2\) and the physical gap still diverges in this group.
This excludes the last finite coefficient pole.

## 4. The scale and the pole at infinity

It follows that \(a,b,\Delta,G\) are polynomials, so
\(\rho=G'/J\) is polynomial. At any finite parameter the paired sums
remain bounded and nonzero in the limit by (5), which has bounded
coefficients, nonzero fixed leading coefficient and nonzero fixed
constant term. Equation (4) then bounds all normalized gaps. If
\(\rho\) vanished there every physical gap would collapse, contrary
to the finite-unit-gap property. Thus \(\rho\) is a nonzero constant
and \(G\) is affine nonconstant.

A nonconstant polynomial \(b\) has a pole at infinity. The same
first pole balance above gives \(\deg G=13k>1\) there.
If \(b\) is constant and \(40b^3-64p\ne0\), then the first
equation of (3) expresses \(a\) as a quadratic polynomial in \(\Delta\)
with nonzero quadratic coefficient. If \(\Delta\) is nonconstant,
\(-3pa^2/8\) gives \(G\) degree \(4\deg\Delta>1\); otherwise
\(a,\Delta,G\) are constant. If \(b\) is constant and
\(40b^3-64p=0\), the first equation instead forces \(\Delta\)
constant. A nonconstant \(a\) would make \(G\) have degree
\(2\deg a>1\), and a constant \(a\) again gives constant \(G\).
All possibilities contradict \(G'=J\rho\ne0\) constant. This
excludes the \((10,4,11)\) embedding.
