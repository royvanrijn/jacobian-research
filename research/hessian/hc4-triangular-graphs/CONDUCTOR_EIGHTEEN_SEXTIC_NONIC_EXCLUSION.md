# Exclusion of the sextic–nonic conductor-eighteen sequence

**Theorem ConductorEighteenSexticNonicExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map in a frame with \(P\)
monic in the longitudinal coordinate \(x\). Its generic exact
directional order cannot have reduced sequence \((9,6,7)\).

We use the [exact conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
and the [finite-pair and parameter identities](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities).
At each finite parameter place there is a matching pair with
distinct finite limits in the original coordinate. Its physical
gap therefore has valuation zero. We may extend the local field
for endpoints and fractional scales. We write \(O(q)\) for a
quantity whose ratio to \(q\) is integral at the place, and
\(o(q)\) when the ratio has positive valuation.

## The two parameter equations

Set \(K_0=\mathbb C(y)\). By
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf),
the degree-seven approximate root is monic of degree three in
the degree-six generator \(W\). The exponents \(0,1,2\) of \(W\)
have distinct weights modulo nine, so terms above weight
eighteen cannot cancel. Completing the square in the monic
degree-nine generator \(V\) gives \(V^2=W^3+O(x^7)\)
with nonzero degree-seven coefficient in the difference.
Center the source so \(W\) has no degree-five term. In this
coordinate take \(\kappa\) to be the degree-six coefficient of
\(V^2-W^3\), giving \(V^2=W^3+\kappa W+T\) with \(T\)
of degree seven and zero degree-six coefficient.
If \(\lambda\) is the target determinant of the changes
before longitudinal rescaling, use \(x-h(y)=\rho z\) with
\(\rho^{15}=\lambda\), and scale the generators by \(\rho^{-6}\)
and \(\rho^{-9}\). These preserve the normal forms and give
target determinant one. Write
\[
 W=z^6+a z^4+b z^3+\cdots,\qquad
 T=\sigma z^7+t_5 z^5+t_4 z^4+\cdots,\qquad \sigma\ne0 .
\]
The parameter identity is
\(\operatorname{Jac}_{z,y}(W,V)=J\rho+\partial_z H(W,V)\),
where \(J\ne0\) is the original Jacobian and \(H\) is a
polynomial over the extended parameter field.

The monic Laurent expansion of \(V\) has the form
\[
 V=W^{3/2}+\frac{\kappa}{2}W^{-1/2}+
             \frac{T}{2W^{3/2}}+O(z^{-13}).
\]
The \(W^{3/2}\) part has zero Jacobian with \(W\).
The \(\kappa\) term contributes
\(\kappa'\partial_z[\sqrt W]_{\ge0}\) to its polynomial part.
For the \(T\) term the relevant coefficients in degrees
\(-2,-4,-5\) are \(\sigma/2,\ t_5/2-3\sigma a/4,\
t_4/2-3\sigma b/4\), respectively, with degree-\(-3\)
coefficient zero. Thus the polynomial Jacobian has degree
at most three. The semigroup \(\langle6,7,9\rangle\) permits
only constants in the order through degree four, so
\(\partial_z H=0\). Comparing degrees three and two gives
\(\sigma'=\kappa'=0\). The terms in degrees one and zero now
give
\[
 B=3t_5-\tfrac72\sigma a\in\mathbb C,
 \qquad G=3t_4-\tfrac72\sigma b,\qquad G'=J\rho .
 \tag{3}
\]
In particular \(\sigma\in\mathbb C^\times\), where the
constant field of the finite parameter extension is
\(\mathbb C\). The coefficient \(\sigma\) has longitudinal
weight eleven, so \(\rho^{11}\in K_0\) by comparison
with the unscaled leading coefficient of \(T\).
Since \(\rho^{15}\in K_0\) as well, \(\rho\in K_0\).
All normalized coefficients below are therefore rational in \(y\).

Write \(W=L^2+M\) with \(L=z^3+kz+l\) monic and
\(M=r z^2+s z+t\). The polynomial part of \(W^{3/2}\)
is exactly \(V\), from the expansion above. In these variables
\[
 V=L^3+\tfrac32 LM+\tfrac38(r^2z+2rs),\qquad
 W^{3/2}-V=\frac{3\,\operatorname{rem}(M^2,L)}{8L}
                    -\frac{M^3}{16L^3}+O(z^{-7}).
 \tag{4}
\]
For \(N_j=[z^{-j}]W^{3/2}\), the coefficients of \(V^2-W^3\)
and (3) give
\[
 N_1=0,\quad N_2=-\sigma/2,\quad N_3=-\kappa/2,\quad
 B=-6N_4+2k\sigma,\quad G=-6N_5+9k\kappa+2l\sigma .
\]
Equations \(N_1=0\) and \(N_2=-\sigma/2\) rule out
\(r\equiv0\) and give, with \(X=s/r\),
\[
 t=r(k-X^2)/2,\qquad l=\frac{4\sigma}{3r^2}-X(k+X^2).
\]
Set \(v=k+3X^2\). Substituting in the next Laurent
coefficients gives
\[
 3r^2 v^2-2r^3-32\sigma X+16\kappa=0,\qquad
 \frac98 r^3 X-\sigma k=B,
 \tag{5}
\]
and
\[
 G=6\kappa k-\frac3{16}kr^3+\sigma kX+
        \frac9{16}r^3X^2+\sigma X^3-\frac{4\sigma^2}{3r^2}.
 \tag{6}
\]
If \(r\) is a nonzero constant, eliminate \(k\) using the
second equation of (5). The first becomes a polynomial
equation for \(X\) over \(\mathbb C\), with nonzero quartic
coefficient \(27r^2\). Thus all coefficients are constant
and \(G'=0\), impossible. Hence \(r\) is nonconstant.

## Matching equations after rescaling

For a local scale \(z=\alpha Z\), write
\(\bar X=X/\alpha\), \(\bar k=k/\alpha^2\),
\(\bar l=l/\alpha^3\), \(\delta=r/\alpha^4\), and put
\[
 \bar L=Z^3+\bar k Z+\bar l,\qquad
 f=Z^2+\bar X Z+(\bar k-\bar X^2)/2.
\]
The scaled polynomials are
\[
 \bar W=\bar L^2+\delta f,\qquad
 \bar V=\bar L^3+\tfrac32\delta\bar L f
               +\tfrac38\delta^2(Z+2\bar X).
\]
Their difference of powers starts in degree two in \(\delta\):
\[
 \bar V^2-\bar W^3
  =\tfrac34\delta^2\bar L^2((Z+2\bar X)\bar L-f^2)
   +\delta^3(\tfrac98\bar Lf(Z+2\bar X)-f^3)
   +\tfrac9{64}\delta^4(Z+2\bar X)^2 .
 \tag{7}
\]
The polynomial in the first parentheses is affine, of slope
\(\bar l+\bar X(\bar k+\bar X^2)=4\sigma/(3r^2\alpha^3)\).

For a matched pair in \(Z\) let \(w\) be its midpoint and
\(g\) its gap. Define
\[
 Q=\bar L'(w)+g^2/4,\qquad
 S_L=2\bar L(w)+3wg^2/2 .
\]
Here primes on \(\bar L,f\) are source derivatives. These
are the divided difference and the sum, respectively, of
the two values of the cubic \(\bar L\).
The exact divided differences of \(\bar W,\bar V\) give
\[
 Q S_L+\delta f'(w)=0,\qquad
 Q f(w)+\delta/4+g^2\left(Q/4+Q^3/(6\delta)\right)=0 .
 \tag{8}
\]
The second follows by subtracting \(3S_L/4\) times the
first divided difference from the divided difference of
\(\bar V\), then dividing by \(3\delta/2\).
The physical gap is \(\rho\alpha g\).

## Finite zeros of \(r\)

Suppose \(r\) vanishes with order \(p>0\). If \(X\) is
bounded, then \(k\) is bounded by (5). If \(X\) has a
pole, \(v\sim 3X^2\) and (5) gives
\(r^2 X^3\to32\sigma/27\); again \(k\) is bounded.
Take \(\alpha=r^{-2/3}\). Then \(\delta=r^{11/3}\),
and the two residue profiles \((\bar L_0,f_0)\) are
\[
 \bigl(Z^3+4\sigma/3,\ Z^2\bigr),\qquad
 \bigl(Z^3+X_0^3/8,\ Z^2+X_0 Z-X_0^2/2\bigr)
\]
respectively, with \(X_0=\lim r^{2/3}X\ne0\) in the
second and \(X_0^3=32\sigma/27\).
In (7) the slope of the affine bracket is \(4\sigma/3\).
Consequently \((\bar V^2-\bar W^3)/\delta^2\) has
regular coefficients with limiting degree seven and
nonzero leading coefficient. Its original order has
equal values at the endpoints, as does the monic sextic
\(\bar W\). These coprime degrees bound all matched
endpoints in \(Z\). Indeed if a pair escaped, its two
leading terms would have ratio a common root of unity
for both degrees, hence one; a divided difference would
then have nonzero leading term.

The leading term of \(G\) in (6) is
\(-4\sigma^2/(3r^2)\) in the bounded case and
\(-4\sigma^2/(27r^2)\) in the unbounded case.
Thus \(\rho=G'/J\) has pole order at least \(2p\).
A physical gap of valuation zero would require
\(\operatorname{ord}g^2\ge16p/3> \operatorname{ord}\delta\).
The endpoints in \(Z\) would have the same limit \(w_0\),
which must be a root of \(\bar L_0\) or of \(\bar L_0'\)
by the first equation in (8). At any root of \(\bar L_0\),
\(Q\) is a unit and \(f_0(w_0)\ne0\). In the second
profile the latter assertion follows from the resultant
\(-27/64\) of \(Z^3+1/8\) and \(Z^2+Z-1/2\) after
scaling by \(X_0\). Then \(Q f(w)\) is a unit in (8)
while both other terms tend to zero, impossible.

The remaining possible limit is \(w_0=0\), where \(S_L\)
is a unit. The first equation in (8) gives \(Q=O(\delta)\).
The second, after substitution and division by \(\delta\),
forces the limit of \(1/4-f(w)f'(w)/S_L\) to be zero,
because its error is \(O(g^2)\). That limit is \(1/4\)
in the bounded profile and \(9/4\) in the unbounded
profile, again impossible. Thus no finite zero of \(r\)
can occur.

## Finite poles of \(r\)

Let \(r\) have pole order \(p>0\). The second equation
in (5) gives \(k=(9r^3 X/8-B)/\sigma\), so
\(v=3X^2+9r^3X/(8\sigma)-B/\sigma\). If the pole order
of \(X\) exceeded \(3p\), the \(X^2\) term would
dominate and \(3r^2v^2\) could not balance (5).
If it is below \(3p\), the \(r^3X\) term controls any
pole of \(v\). Equation (5) requires \(v^2\sim 2r/3\),
so \(v\) has pole order \(p/2\) and \(X\) vanishes to
order \(5p/2\). In that case \(k\sim v\),
\(r\sim 3k^2/2\), \(s,l\) are bounded, and
\(t\sim rk/2\). If the pole order of \(X\) equals \(3p\),
the two leading terms in \(v\) have to cancel and
\[
 X/r^3\ \longrightarrow -3/(8\sigma),
 \qquad v^2/r\ \longrightarrow -10/3 .
 \tag{9}
\]
These exhaust the possibilities. An identically zero \(X\)
would make \(k\) constant and force \(r\) constant by (5).

### The branch with \(X\) tending to zero

Scale by \(\alpha=\sqrt{k}\). The limit has
\(L_0=Z^3+Z,\ M_0=\tfrac32(Z^2+1/2)\), so
\[
 W_0=L_0^2+M_0,\qquad
 V_0=L_0^3+\tfrac32 L_0M_0+\tfrac{27}{32}Z .
\]
The polynomial \(T_0=V_0^2-W_0^3\) has degree four,
with derivative \(-27Z(16Z^2+13)/512\).
This is the immersed limiting profile of the
[degree-nine cable proof](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#3-sequence-964).
By (6), \(G\sim-3kr^3/16\) is polar, so \(\rho\alpha\)
is polar. A unit physical gap would force the scaled gap to
tend to zero. The endpoints would then be bounded: if they
were unbounded with the same leading term, the monic sextic
\(\bar W\) with bounded coefficients would have a nonzero
leading divided difference. The scaled difference of powers
converges coefficientwise to \(T_0\), so a common endpoint
limit would be a simultaneous critical point of \(W_0,V_0,T_0\).
But \(V_0'(0)=63/32\) and \(W_0'/(2Z)=315/256\) at
\(16Z^2+13=0\), so their derivatives have no common zero.

### The branch with \(X\) of pole order \(3p\)

Scale by \(\alpha=X\), and put
\[
 \varepsilon=v/X^2,\qquad
 \delta=r/X^4,\qquad \delta_L=4\sigma/(3r^2X^3).
\]
Then
\[
 \bar L=(Z-1)^2(Z+2)+\varepsilon(Z-1)+\delta_L,\qquad
 f=(Z-1)(Z+2)+\varepsilon/2 .
\]
The valuations of \(\varepsilon,\delta\) are \(11p/2,11p\),
and (9) gives
\[
 \delta/\varepsilon^2\longrightarrow -3/10,\qquad
 \delta_L/\delta\longrightarrow -1/2,\qquad
 \delta_L/\varepsilon^2\longrightarrow 3/20.
\]
The affine bracket in (7) is
\(\delta_L Z+2\delta_L-\varepsilon^2/4\). Accordingly
\((\bar V^2-\bar W^3)/\delta^3\) has regular coefficients
and has limiting degree seven, with nonzero leading
coefficient from \(\delta_L/\delta\). Together with the
monic sextic \(\bar W\) this bounds all scaled endpoints.

By (6) and (9), \(G\sim-5\sigma X^3\). The scaling
factor \(\rho\alpha\) has pole order at least \(12p\)
and strictly greater than \(12p\) at a finite place,
because \(G\) has pole order \(9p\) and differentiation
increases it. Thus a unit physical gap would satisfy
\(\operatorname{ord}(g^2)>24p\). Its midpoint limit must
be one of \(-2,1,-1\), the zeros of
\(\bar L_0\bar L_0'\).

- At \(w\to-2\), \(Q\to9\). The first identity in (8)
  gives \(\bar L(w)=O(\varepsilon^2)\), and hence
  \(w=-2+\varepsilon/3+O(\varepsilon^2)\).
  Then \(f(w)=-\varepsilon/2+O(\varepsilon^2)\).
  This term in the second equation of (8) cannot
  cancel: \(\delta\) has order \(11p\), and
  \(g^2Q^3/\delta\) has order greater than \(13p\).
- At \(w\to1\), write \(t_*=w-1\). If \(t_*\) has
  smaller valuation than \(\varepsilon\), then
  \(Q f(w)\sim18t_*^2\) dominates \(\delta/4\) and
  the \(g^2\) terms of the second equation in (8).
  Otherwise \(t_*=O(\varepsilon)\), so
  \(Q=O(\varepsilon)\) and \(S_L=O(\varepsilon^2)\).
  In the first equation \(\delta f'(w)\) has exact
  order \(11p\) because \(f'(w)\to3\), and cannot
  be canceled by \(Q S_L=O(\varepsilon^3)\).
- At \(w\to-1\), \(\bar L(w)\to4\) and \(f'(w)\to-1\).
  The first equation of (8) gives \(Q=-\delta f'(w)/S_L\),
  of exact order \(11p\), and
  \(\bar L'(w)\sim-3\varepsilon^2/80\). Expanding
  \(3(w^2-1)+\varepsilon\) gives
  \[
  w=-1+\varepsilon/6+(29/1440)\varepsilon^2
          +o(\varepsilon^2).
  \]
  With this midpoint,
  \[
  \bar L(w)/2-f(w)f'(w)
      = (3/32)\varepsilon^2+o(\varepsilon^2).
  \]
  Yet the second equation of (8), after substituting
  for \(Q\) and dividing by \(\delta\), requires
  \(1/4-f(w)f'(w)/S_L=O(g^2)\), since
  \(Q/\delta\) and \(Q^3/\delta^2\) are integral.
  Replacing \(S_L\) by \(2\bar L(w)\) changes this
  expression by \(O(g^2)\). Its actual leading
  term is \(3\varepsilon^2/256\), of order \(11p\),
  a contradiction.

This excludes a finite pole of \(r\). The rational function
\(r\) has neither a finite zero nor a finite pole and must
be constant, contrary to (3)–(6). This proves the exclusion.

## Exact controls

The [checker](verify_conductor_eighteen_sextic_nonic_exclusion.py)
verifies the Laurent coefficients in (4)–(6), the polynomial
difference in (7), the exact midpoint equations (8), and the
nonzero resultant and local midpoint computations. The
valuation case analysis and its use of finite physical gaps
are the written parts of the proof. Replay from the repository
root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_eighteen_sextic_nonic_exclusion.py
~~~
