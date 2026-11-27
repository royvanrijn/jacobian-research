# The quartic–decic sequence \((10,4,13)\) in conductor degree twenty

**Theorem ConductorTwentyQuarticDecicExclusion.** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\). Its exact
directional coordinate order cannot have reduced generic characteristic
sequence \((10,4,13)\), even when the specialized conductor degree drops.

## 1. The cubic correction and the parameter invariants

Write \(K=\mathbb C(y)\), and let
\(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The [coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
supplies the exact finite birational conductor on every source line.
In the hypothetical noninvertible case it is nonzero and nonconstant
on each finite line.

By the [one-place approximate-root theorem](sources/assi-garcia-sanchez-1407.0490v1.pdf),
Propositions 2 and 4, the monic reduced generators \(U,V\) of degrees
four and ten have an approximate root of source degree thirteen.
The monic coefficient hypotheses are verified for exact directional
orders as in the [cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).
As in the [preceding quartic–decic exclusion](CONDUCTOR_EIGHTEEN_QUARTIC_DECIC_EXCLUSION.md),
the weights of \(U^i\), \(0\le i<5\), are distinct modulo ten.
Thus completing the square in \(V\) and translating \(U\) expresses
the approximate root as \(V^2-R_0(U)\), where \(R_0\) is monic quintic
with no quartic term.

Depress \(U\) by a longitudinal translation \(h(y)\in K\).
For the determinant \(\lambda\in K^*\) of the target changes choose
\(\rho^{14}=\lambda\), put \(x=h(y)+\rho z\), and scale the
two generators by \(\rho^{-4},\rho^{-10}\). Continue with \(U,V\)
for the scaled monic generators. We have
\[
 U=z^4+az^2+bz+c,\qquad \Delta=a^2-4c.
\]
Since \(V-R_0(U)^{1/2}\) has highest exponent three in \(z\),
\[
 V=[U^{5/2}]_{\ge0}+p z^3+q z^2+r z
\tag{1}
\]
after a constant translation in the target, with \(p\ne0\).
Indeed the polynomial part of \(R_0(U)^{1/2}-U^{5/2}\) has
degree at most two, so is already included in the correction
in (1).

The normalized target determinant is one. Its parameter-dependent
target field is Hamiltonian as in the
[degree-twelve cable proof](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities),
giving
\[
 U_z V_y-U_y V_z=J\rho+\partial_z H(U,V).
\tag{2}
\]
The contribution of \(U^{5/2}\) cancels by the chain rule, so the
left side has degree at most six. In the semigroup
\(\langle4,10,13\rangle\) this forces the evaluated \(H\) to belong
to the span of \(1,U\). Degrees six, five and four in (2) give
\(p'=q'=0\) and \(r'=3pa'/4\). Hence
\[
 p\in\mathbb C^*,\qquad q=Q_0\in\mathbb C,\qquad
 r=R+3pa/4,\quad R\in\mathbb C.
\]
The degree-three coefficient determines the coefficient of \(U\)
in \(H\). For the remaining comparisons put
\[
 \begin{split}
 A&=-\bigl(15\Delta^2b-48\Delta p+128Q_0b+64Ra
                         +24a^2p-40ab^3\bigr)/64,\\
 B&=\bigl(5\Delta^3+128\Delta Q_0-120\Delta ab^2-256Rb
                         +192abp+40b^4\bigr)/256,\\
 G&=\bigl(8\Delta R+6\Delta ap-5\Delta b^3-12Ra^2
                          -5a^3p+12b^2p\bigr)/32-aA/4 .
 \end{split}
\tag{3}
\]
The [exact checker](verify_conductor_twenty_quartic_decic_exclusion.py)
verifies that after subtracting \(\partial_zH\) the coefficients
of degrees two, one and zero are \(A',B',G'+aA'/4\).
Consequently \(A,B\) are constants and
\[
 G'=J\rho .
\tag{4}
\]
The nonzero constant \(p\) has scale weight seven, so its unscaled
value \(\rho^7p\) lies in \(K\). Thus \([K(\rho):K]\) divides seven.
The coefficients of the normalized pair belong to \(K(\rho)\).

## 2. Matched-pair equations and the finite unit gap

If \(b\equiv0\), the constant \(B\) in (3) makes \(\Delta\) constant,
and \(A\) then makes \(a\) constant since \(p\ne0\).
This contradicts (4). For distinct matched endpoints in the \(z\)
coordinate let \(s\) be their sum and \(g\) their gap.
The quartic divided difference makes \(s\ne0\) over the generic
parameter and gives \(g^2=-s^2-2a-2b/s\).
Substituting the resulting product in the decic divided difference
and multiplying by \(128s\) gives the polynomial equation
\[
\begin{split}
 F(s)={}&8s^6(s^2+2a)^2+20\Delta s^4(s^2+2a)-8bs^5(s^2+2a)
                -8b^2s^4\\
       &+(64p-20\Delta b)s^3+(15\Delta^2+128Q_0-40ab^2)s^2
                +(128R+32ap)s-64bp=0 .
\end{split}
\tag{5}
\]
It is of degree ten with constant leading coefficient. Equivalently,
with \(Y=-g^2=s^2+2a+2b/s\), the equation is
\[
\begin{split}
 0={}&8s^6Y^2+(20\Delta s^4-40bs^5)Y
      +40b^2s^4+(64p-60\Delta b)s^3\\
    &+(15\Delta^2+128Q_0-40ab^2)s^2+(128R+32ap)s-64bp .
\end{split}
\tag{6}
\]
The checker verifies both identities directly.

At every finite parameter at least one generic pair has finite distinct
limits in the original \(x\)-coordinate. Indeed a generic conductor root
can be chosen to reach a finite root of the specialized conductor,
and monicity of the original \(P\) keeps its mate finite. Equal limits
would then force both \(P_x\) and \(Q_x\) to vanish by their divided
differences. This is the finite-unit-gap property used in the preceding
cable exclusions. The physical gap is \(\rho g\), so it must be a unit
for this pair.

## 3. The three possible coefficient-pole shapes

At any place of \(K(\rho)\) where a normalized coefficient has a pole,
choose the positive weight scale by the maximum of the pole orders
of \(a,b,\Delta\) divided by \(2,3,4\), respectively.
The highest-weight parts of the two constant invariants (3) are
\[
 b_*(15\Delta_*^2-40a_*b_*^2)=0,\qquad
 5\Delta_*^3-120\Delta_*a_*b_*^2+40b_*^4=0 .
\tag{7}
\]
The starred coefficients are leading coefficients in those three
weighted orders, with zero allowed for a smaller order.

If \(b_*\ne0\), equations (7) give
\(\Delta_*^3=b_*^4\) and \(a_*=3\Delta_*^2/(8b_*^2)\).
Thus all three are nonzero. The pole orders of \(a,b,\Delta\)
are \(2k,3k,4k\) with a positive integer \(k\), and
\[
 G\sim-5\Delta b^3/32 .
\tag{8}
\]
Here \(\sim\) denotes equality of nonzero leading terms at the given
place. This is the \(b\)-leading shape.

If \(b_*=0\), then \(\Delta_*=0\) by (7); only \(a\) has the
top weighted order. Let \(n=-v(a)>0\), \(m=-v(b)<3n/2\)
and \(d=-v(\Delta)<2n\). The invariants force \(m=n/3\).
Indeed, if \(m>n/3\), the term \(ab^3\) in \(A\) must balance
\(\Delta^2b\), yielding
\(\Delta^2\sim(8/3)ab^2\) and \(d=n/2+m\).
But then the dominant terms of \(B\) require
\(\Delta^2\sim24ab^2\). If \(m<n/3\), the dominant \(pa^2\)
in \(A\) must balance \(\Delta^2b\), giving
\(2d=2n-m\). Then \(\Delta^3\) strictly dominates \(B\).
Both cases are impossible. These comparisons also cover bounded
or vanishing \(b\).

Now \(m=n/3\) and \(d\le5n/6\), since a larger \(d\) would
make \(\Delta^2b\) dominate \(A\). When \(d<5n/6\) the leading
balance in \(A\) yields \(b^3\sim(3/5)pa\). In \(B\), the
term \(\Delta ab^2\) would dominate if \(d>-n/3\).
At the other extreme, if \(d<-n/3\), the combined leading
coefficient \(192abp+40b^4\sim216abp\) cannot cancel.
Hence \(d=-n/3\) and \(\Delta\sim9p/(5b)\).
When \(d=5n/6\), the leading terms of \(B\) give
\(\Delta^2\sim24ab^2\), and the leading terms of \(A\)
then give \(b^3\sim-3pa/40\).
The two \(a\)-leading shapes are consequently
\[
 \begin{array}{lll}
 \mathrm{I}:& b^3\sim(3/5)pa,&\Delta\sim9p/(5b),\\
 \mathrm{II}:& b^3\sim-3pa/40,&\Delta^2\sim24ab^2 .
 \end{array}
\tag{9}
\]
In both, (3) makes \(G\sim-5pa^3/32\), of pole order \(3n\ge9\),
since \(n=3m\) with \(m\) a positive integer.

## 4. None of the coefficient poles is finite

Suppose first that a finite coefficient pole is \(b\)-leading.
By (8), \(G\) has pole order \(13k\), so \(\rho=G'/J\) has a pole.
Rescale \(z=\sigma w\), where \(\sigma=\Delta/b\) has pole order \(k\).
The polynomials \(\sigma^{-4}U(\sigma w)\) and
\(\sigma^{-10}V(\sigma w)\) have bounded coefficients and limits
\[
 U_*=w^4+\tfrac38w^2+w-\tfrac{55}{256},\qquad
 V_*=[U_*^{5/2}]_{\ge0}.
\]
The terms correcting the polynomial part in (1) vanish in this scale.
The limiting derivatives have resultant \(208000\ne0\), as checked
by the exact replay. A unit physical gap would make the gap in
\(w\) tend to zero, since \(\rho\sigma\) has a pole. Such close
endpoints must be bounded: two unbounded endpoints at a vanishing gap
would have equal leading terms and their monic quartic divided
difference would have nonzero leading term. Bounded close endpoints
would give a common critical point of \(U_*,V_*\), a contradiction.

Now suppose the finite pole is \(a\)-leading, in one of (9).
The coefficient valuations of (5), divided by \(n\), are the
following, listing the two shapes separately.
\[
\begin{array}{c|rrrrrrrrrrr}
 \deg_s&0&1&2&3&4&5&6&7&8&9&10\\ \hline
 \mathrm{I}&-1/3&-1&-5/3&0&-2/3&-4/3&-2&-1/3&-1&+\infty&0\\
 \mathrm{II}&-1/3&-1&-5/3&-7/6&-11/6&-4/3&-2&-1/3&-1&+\infty&0
\end{array}
\]
In both rows the lower hull joins degrees zero, two, six, ten.
There are two small sums of order \(2n/3\), four of order \(n/12\),
and four large sums of order \(-n/2\). Multiplicities are counted
after local extension. For the smallest sums, write \(s\sim (b/a)t\).
The first face in either case, after division by \(pb\), gives
\(-24t^2+32t-64=0\). Its roots are nonzero and none is \(-1\),
so the leading terms of \(-2a-2b/s\) in the gap do not cancel.
The four intermediate sums have \(g^2\sim-2a\) directly.

For the large sums \(s^2\sim-2a\) by the last face of (5).
Use equation (6). The quadratic coefficient in \(Y\) has
valuation \(-3n\) and the linear coefficient valuation at least
\(-17n/6\). The constant coefficient has valuation exactly
\(-8n/3\). In shape I its leading value is
\(240a^2b^2\). In shape II the leading expression is
\(-480a^2b^2-60\Delta b s^3\). Cancellation would force
\(\Delta s\sim4ab\), by \(s^2\sim-2a\); squaring would contradict
\(\Delta^2s^2\sim-48a^2b^2\).
Thus the Newton polygon in \(Y\) is the straight segment
of slope \(-n/6\), so \(v(Y)=n/6\). Even these shrinking
normalized gaps diverge physically, since (4) and
\(G\sim-5pa^3/32\) make \(\rho\) have pole order greater than \(3n\)
at a finite place. This excludes the \(a\)-leading coefficient poles.

## 5. The scale cannot account for a pole at infinity

The normalized coefficients are now regular at every finite place.
If \(\rho\) had a finite zero, (5) would keep all normalized sums
bounded. Dividing that equation by \(s\) bounds \(b/s\): its
coefficient is the nonzero constant \(-64p\) and all other terms
are polynomials in \(s\) with bounded coefficients.
Thus the normalized gaps are bounded and every physical gap
would shrink to zero, contrary to the finite unit gap.

Since \(\rho^7\in K\) has no finite zero, it has no pole at infinity:
a rational function of \(y\) with no finite zero has constant
numerator in lowest terms. In particular \(\rho\) has no pole
over infinity in \(K(\rho)\). If some normalized coefficient is
nonconstant, it must have a pole on the smooth complete parameter
curve, and all such poles lie over \(y=\infty\).
At that place the preceding classification forces \(G\) to have
pole order at least thirteen in the \(b\)-leading case or at
least nine in the \(a\)-leading case. The ramification index \(e\)
of \(y\) at that place is at most seven. Differentiation of a pole
with respect to \(y\) decreases its pole order by \(e\);
hence \(G'\) still has a pole. This contradicts (4) and the absence
of a pole of \(\rho\) there. If all coefficients are constant,
then \(G'=0\), which also contradicts (4). The sequence is excluded.
