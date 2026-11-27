# Exclusion of the sextic–nonic conductor-twenty sequence

**Theorem ConductorTwentySexticNonicExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map in a frame with \(P\)
monic in the longitudinal coordinate \(x\). Its exact generic
directional order cannot have reduced sequence \((9,6,8)\).

We use the [finite-pair property](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
of the [exact conductor](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md):
at every finite parameter place some generic matching pair has
distinct finite original limits. Its physical gap has valuation
zero. Valuations below may be extended to the fields of the pair
and of a fractional scale. At such a place \(O(a)\) means the
ratio to \(a\) is integral and \(o(a)\) means the ratio has positive
valuation.

## The Laurent equations

Put \(\mathcal K=\mathbb C(y)\). For the reduced generators
\(W,V\) of degrees six and nine,
[Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf)
give the degree-eight approximate root, monic of degree three
in \(W\). The powers of \(W\) of exponents zero, one and two
have distinct weights modulo nine, so no terms of weight
above eighteen can cancel. Normalize both generators to
be monic, complete the square in \(V\), and center the
source so \(W\) has no degree-five term. The resulting
generators satisfy \(V^2-W^3\) of degree eight, and
\[
 W=z^6+a z^4+b z^3+c z^2+\cdots
\]
in the normalized source variable. For this normalization
let \(\lambda\) be the determinant of the target change and use
\(x-h(y)=\rho z\), \(\rho^{15}=\lambda\), with the generators
scaled by \(\rho^{-6},\rho^{-9}\). The target determinant
becomes one. The parameter identity is
\[
 \operatorname{Jac}_{z,y}(W,V)=J\rho+\partial_z H(W,V)
 \tag{1}
\]
for a polynomial \(H\) over the finite parameter extension,
where \(J\ne0\) is the original Jacobian.

The monic Laurent expansion of \(W^{3/2}\) has polynomial
part \(V\), because \(V^2-W^3\) has degree eight.
Write \(N_j=[z^{-j}]W^{3/2}\).
The chain rule cancels the entire \(W^{3/2}\) part of the
Jacobian, leaving the polynomial part of
\(-W_z\partial_y\sum_j N_jz^{-j}
  +W_y\partial_z\sum_j N_jz^{-j}\).
This has degree at most four, and the semigroup
\(\langle6,8,9\rangle\) makes \(H(W,V)\) constant in \(z\).
Comparing degrees four through zero in (1) gives
\[
\begin{gathered}
 \sigma=-2N_1\in\mathbb C^\times,\qquad \tau=-2N_2\in\mathbb C,\\
 K=6N_3-a\sigma/2\in\mathbb C,\qquad
 B=6N_4-a\tau-b\sigma/2\in\mathbb C,\\
 G=-6N_5-Ka/2-7\sigma a^2/24+\tau b+\sigma c/2,\qquad
 G'=J\rho .
 \tag{2}
\end{gathered}
\]
For example, after \(N_1,N_2\) become constant, the
coefficient of degree two is \(-6N_3'-a'N_1\).
The coefficient of degree zero before substituting \(K\) is
\(-6N_5'-4aN_3'-3a'N_3-2b'N_2-c'N_1\).
The constant field of the finite parameter extension is
\(\mathbb C\). The coefficient \(\sigma\) is the nonzero
coefficient of \(z^8\) in \(V^2-W^3\).

Write \(W=L^2+M\) with \(L=z^3+kz+l\) monic and
\(M=r z^2+s z+t\). As in the
[conductor-eighteen sextic–nonic proof](CONDUCTOR_EIGHTEEN_SEXTIC_NONIC_EXCLUSION.md),
\[
 V=L^3+\tfrac32 LM+\tfrac38(r^2z+2rs),\qquad
 W^{3/2}-V=\frac{3\,\operatorname{rem}(M^2,L)}{8L}
                    -\frac{M^3}{16L^3}+O(z^{-7}).
 \tag{3}
\]
If \(r\) vanished identically, the four constants
\(\sigma,\tau,K,B\) would successively fix \(s,t,k,l\)
by (2)–(3). In fact \(\sigma=-3s^2/4\) implies \(s\ne0\).
This would make \(G\) constant, so \(r\ne0\) in the
parameter field.

With \(X=s/r\), \(v=k+3X^2\), the first two tail
coefficients give
\[
 t=\frac r2(k-X^2)-\frac{2\sigma}{3r},\qquad
 l=\frac{4(\tau-\sigma X)}{3r^2}-X(k+X^2).
 \tag{4}
\]
The next two equations and the reduced primitive are
\[
\begin{aligned}
0={}&9r^4v^2+8\sigma r^2v+96\sigma r^2X^2-96\tau r^2X
              -6r^5+16\sigma^2-16Kr^2,\\
0={}&(48\sigma X-24\tau)v-96\sigma X^3+72\tau X^2+27r^3X
              +64\sigma(\sigma X-\tau)/r^2+24B,\\
G={}&3\sigma X^4-2\tau X^3+\tfrac98 r^3X^2-2\sigma X^2v+\tau Xv
       +\frac{4\sigma\tau X}{3r^2}
       -\tfrac3{16}r^3v-\tfrac14r\sigma+\tfrac13\sigma v^2
       -\frac{4\tau^2}{3r^2}.
\end{aligned}
\tag{5}
\]
If \(r\) is a nonzero constant and \(X\) nonconstant, choose
a pole of \(X\) on the complete parameter curve. The second
equation implies \(v\sim 2X^2\) there. The quartic term
\(9r^4v^2\) in the first equation then dominates all the
others, impossible. If both \(X,r\) are constant, \(v\)
also is constant by the first equation of (5), again
contradicting \(G'\ne0\). Thus \(r\) is nonconstant.

## Rescaling the matching equations

Use a local scale \(z=\alpha Z\), and put
\(\bar X=X/\alpha,\ \bar k=k/\alpha^2,\ \bar l=l/\alpha^3,\
\delta=r/\alpha^4\). With
\[
 \bar L=Z^3+\bar k Z+\bar l,\qquad
 f=Z^2+\bar X Z+(\bar k-\bar X^2)/2-\frac{2\sigma}{3r^2\alpha^2},
\]
the rescaled pair is
\[
 \bar W=\bar L^2+\delta f,\qquad
 \bar V=\bar L^3+\tfrac32\delta\bar L f
               +\tfrac38\delta^2(Z+2\bar X).
\]
Its difference of powers is
\[
 \bar V^2-\bar W^3
 =\tfrac34\delta^2\bar L^2((Z+2\bar X)\bar L-f^2)
 +\delta^3\bigl(\tfrac98\bar L f(Z+2\bar X)-f^3\bigr)
 +\tfrac9{64}\delta^4(Z+2\bar X)^2 .
 \tag{6}
\]
For a matching pair in \(Z\) write \(w\) for its midpoint,
\(g\) for its gap, \(Q=\bar L'(w)+g^2/4\), and
\(S_L=2\bar L(w)+3wg^2/2\). The exact divided differences
are
\[
 Q S_L+\delta f'(w)=0,\qquad
 Q f(w)+\delta/4+g^2\left(Q/4+Q^3/(6\delta)\right)=0 .
 \tag{7}
\]
Indeed \(Q\) is the divided difference of \(\bar L\) and
\(S_L\) its sum at the two endpoints. Subtract \(3S_L/4\)
times the divided difference of \(\bar W\) from that of
\(\bar V\) and divide by \(3\delta/2\). This gives (7)
for any monic quadratic \(f\), including its present
constant term. The physical gap is \(\rho\alpha g\).

If polynomials with regular coefficients and unit leading
coefficients of degrees with gcd one agree at a matching
pair, both endpoints are bounded. Otherwise their leading
terms have a ratio that is a common root of unity for all
these degrees, hence one. The divided difference of any
one of them would then have nonzero leading term. We use
this with degrees six, eight and nine.

## Places with \(r\) small

Suppose \(r\) vanishes to order \(p>0\). Write the pole
order of \(X\) as \(q\) when \(X\) has a pole. If \(q>p\),
the first equation of (5) forces \(v\) to have pole
order \(p+q\), by balancing its quadratic term with
\(96\sigma r^2X^2\). The term \(-96\sigma X^3\) in
the second equation then dominates, impossible.
For \(q<p\) or bounded \(X\), the first equation gives
a residue \(D=\lim r^2v\) with
\(9D^2+8\sigma D+16\sigma^2=0\).
A pole of \(X\) in this range would make the leading
coefficient of the second equation, in units of
\(X/r^2\), equal \(16\sigma(3D+4\sigma)\ne0\).
Thus \(X\) is bounded in this case.

It remains to consider \(q=p\), with \(u=rX\) and
\(D=r^2v\) integral. The exact equations become
\[
\begin{aligned}
0={}&9D^2+8\sigma D+96\sigma u^2-96\tau r u-6r^5
            +16\sigma^2-16Kr^2,\\
0={}&(48\sigma u-24\tau r)D-96\sigma u^3+72\tau r u^2
             +27r^5u+64\sigma(\sigma u-\tau r)+24Br^3.
\end{aligned}
\tag{8}
\]
If \(C=\lim u\ne0\) then \(D_0=\lim D=2C^2-4\sigma/3\).
Substituting into the first equation gives
\(C^2=-4\sigma/9\) or \(C^2=-4\sigma/3\).

### The removable profile

If \(C^2=-4\sigma/3\), \(D_0=3C^2\).
The two equations (8) have invertible Jacobian in \(u,D\)
at \(r=0,u=C,D=3C^2\), with determinant \(-1728C^6\).
Successive coefficient comparison, or formal implicit
solution, gives
\[
\begin{aligned}
 u&=C+\frac{2\tau}{3C^2}r+
          \frac{KC^2-3\tau^2}{3C^5}r^2+O(r^3),\\
 D&=3C^2+\frac{4\tau}{C}r+
          \frac{4(KC^2-3\tau^2)}{3C^4}r^2+O(r^3).
\end{aligned}
\]
Therefore \(r^2k=D-3u^2=O(r^2)\) and
\(r^3l=4(\tau r-\sigma u)/3-u(D-2u^2)=O(r^3)\).
Also \(s=u\) is regular and
\(t=rk/2+(-u^2/2-2\sigma/3)/r\) is regular.
Thus this is a regular chart for \(W,V,G\), even though
\(X=s/r\) has a pole.

### The two polar profiles

In the bounded-\(X\) case and the case \(C^2=-4\sigma/9\),
take \(\alpha=1/r\), so \(\delta=r^5\).
The residue pairs \((\bar L_0,f_0)\) are, respectively,
\[
 \bigl(Z^3+D Z,\ Z^2+D/2-2\sigma/3\bigr),\qquad
 \bigl(Z^3+2C^2Z,\ Z^2+CZ+2C^2\bigr),
\]
where \(D=\lim r^2v\) in the first. Equation (6) shows
that \((\bar V^2-\bar W^3)/\delta^2\) has regular
coefficients. Its actual degree is eight with leading
coefficient \(\sigma\), since the leading coefficient of
\(V^2-W^3\) is \(\sigma\). With \(\bar W,\bar V\),
the degrees six, eight and nine bound all scaled endpoints.

By (5), the leading coefficients of \(r^4G\) at the two
profiles are \(\sigma D^2/3\) and \(64\sigma^3/243\),
respectively. They are nonzero. At a finite place,
differentiation of a pole increases its order, so
\(\rho\alpha\) has pole order greater than \(5p\).
A unit physical gap would satisfy
\(\operatorname{ord}(g^2)>10p>\operatorname{ord}\delta\).
Its midpoint would tend to a root of \(\bar L_0\) or
\(\bar L_0'\) by (7). At any root of \(\bar L_0\),
\(Q\) is a unit and \(f_0\ne0\). For the bounded
profile, this uses \(D\ne0,\ \pm4\sigma/3\), by its
quadratic equation. In the second profile it follows
by substitution at \(0\) and at \(Z^2=-2C^2\).
The second equation of (7) then cannot hold.

At a root of \(\bar L_0'\), \(S_L\) is a unit.
Substitution from the first equation of (7) into the
second requires \(1/4-f(w)f'(w)/S_L=O(g^2)\).
Its limit is \(\sigma/D\ne0\) in the bounded profile.
For the second, put \(Z=Cz\) and use \(z^2=-2/3\);
the limit is \(1/4-11/8=-9/8\).
Thus neither polar profile occurs over a finite parameter.

## Poles of \(r\)

Let \(r\) have pole order \(p>0\), and let \(q\) be the
pole order of \(X\), allowing \(q\le0\).
The first equation of (5) bounds the pole order of \(v\)
by \(\max(p/2,q-p)\): above this bound, \(9r^4v^2\)
would strictly dominate every other term. If \(q>3p/2\),
the second equation would have unmatched leading term
\(-96\sigma X^3\). If \(q<3p/2\), the first equation
gives \(v^2\sim2r/3\), and the second then gives
\(X=O(r^{-5/2})\) by comparison with \(27r^3X\).
This includes \(X\equiv0\).
In this case \(k\sim v\), \(r\sim3k^2/2\), and (4)
gives \(l,s\) bounded and \(t\sim rk/2\).

If \(q=3p/2\), balancing \(-96\sigma X^3\) with
\(27r^3X\) in the second equation of (5), and substituting
in the first, yields
\[
 X^2/r^3\longrightarrow9/(32\sigma),\qquad
 v^2/r\longrightarrow-7/3 .
 \tag{9}
\]
These cases exhaust all poles of \(r\).

### The nonsquare limiting profile

In the case \(X=O(r^{-5/2})\), take \(\alpha=\sqrt k\).
The scaled limit is
\[
 L_0=Z^3+Z,\quad M_0=\tfrac32(Z^2+1/2),\quad
 W_0=L_0^2+M_0,\quad
 V_0=L_0^3+\tfrac32L_0M_0+\tfrac{27}{32}Z .
\]
The primitive satisfies \(G\sim-3r^3v/16\), so
\(\rho\alpha\) has a pole and a unit physical gap
requires \(g\to0\). Endpoints with this gap cannot
escape: the divided difference of the monic sextic
\(\bar W\) with bounded coefficients would have a
nonzero leading term if both endpoints escaped with
ratio tending to one. Therefore such endpoints would
coalesce at a common critical point of \(W_0,V_0\)
and \(V_0^2-W_0^3\). The
[degree-nine cable profile](CONDUCTOR_EIGHTEEN_SEXTIC_NONIC_EXCLUSION.md#the-branch-with-x-tending-to-zero)
has no such point. Its difference-of-powers derivative
is \(-27Z(16Z^2+13)/512\), while \(V_0'(0)=63/32\)
and \(W_0'/(2Z)=315/256\) at \(16Z^2+13=0\).
Thus this pole is impossible at a finite parameter.

### The square limiting profile

In case (9), take \(\alpha=X\) and set
\[
 \epsilon=v/X^2,\qquad \delta=r/X^4,\qquad
 \delta_L=4(\tau-\sigma X)/(3r^2X^3),\qquad
 \mu=-2\sigma/(3r^2X^2).
\]
Then the two polynomials used in (6)–(7) are
\[
 \bar L=(Z-1)^2(Z+2)+\epsilon(Z-1)+\delta_L,\qquad
 f=(Z-1)(Z+2)+\epsilon/2+\mu ,
\]
and
\[
 \operatorname{ord}\epsilon=5p/2,\quad
 \operatorname{ord}\delta=5p,\quad
 \delta/\epsilon^2\longrightarrow-3/7,\quad
 \delta_L/\epsilon^2\longrightarrow9/56,\quad
 \mu/\epsilon^2\longrightarrow9/112.
\]
The polynomial \((Z+2)\bar L-f^2\) has all its
coefficients \(O(\delta)\): for \(\mu=0\) it is
\(\delta_L Z+2\delta_L-\epsilon^2/4\), and the
additional terms contain \(\mu\). Hence (6) makes
\((\bar V^2-\bar W^3)/\delta^3\) regular, with actual
degree eight and leading coefficient
\(\sigma X^2/r^3\to9/32\). Coprime degrees collectively
bind the endpoints as before.

Here \(G\sim7\sigma X^4\) by (5). At a finite place
\(\rho\alpha\) has pole order greater than \(15p/2\).
A unit physical gap would require \(\operatorname{ord}g^2>15p\).
Its limiting midpoint is \(-2,1\) or \(-1\), a zero
of \(\bar L_0\bar L_0'\).

- At \(w\to-2\), \(Q\to9\). The first equality
  in (7) gives \(\bar L(w)=O(\epsilon^2)\), hence
  \(w=-2+\epsilon/3+O(\epsilon^2)\). Then
  \(f(w)=-\epsilon/2+O(\epsilon^2)\). Its term
  \(Q f(w)\) in the second equation of (7) has
  order \(5p/2\), while \(\delta/4\) has order \(5p\)
  and the \(g^2\) terms have order greater than \(10p\).
- At \(w\to1\), write \(t_*=w-1\). If it has
  smaller valuation than \(\epsilon\), \(Q f(w)\sim
  18 t_*^2\) is the sole leading term of the second
  equation in (7). Otherwise \(Q=O(\epsilon)\)
  and \(S_L=O(\epsilon^2)\). Then the first equation
  cannot hold because \(\delta f'(w)\) has exact
  order \(5p\) while \(Q S_L=O(\epsilon^3)\).
- At \(w\to-1\), \(\bar L(w)\to4\), \(f'(w)\to-1\).
  The first equation of (7) gives \(Q=-\delta f'(w)/S_L\)
  of exact order \(5p\), so
  \(\bar L'(w)\sim-3\epsilon^2/56\). Since
  \(\bar L'(w)=3(w^2-1)+\epsilon\), this gives
  \[
    w=-1+\epsilon/6+(23/1008)\epsilon^2+o(\epsilon^2).
  \]
  At that midpoint,
  \(\bar L(w)/2-f(w)f'(w)=(3/16)\epsilon^2+o(\epsilon^2)\).
  Substituting \(Q\) into the second equation of (7)
  and dividing by \(\delta\) requires
  \(1/4-f(w)f'(w)/S_L=O(g^2)\), since \(Q=O(\delta)\).
  Replacing \(S_L\) by \(2\bar L(w)\) changes this
  by \(O(g^2)\). Its actual leading term is
  \(3\epsilon^2/128\) of order \(5p\), contrary to
  \(\operatorname{ord}g^2>15p\).

This excludes the second type of finite pole as well.

## End of the parameter argument

At a place where \(r\) is a unit, \(X\) is integral.
Indeed a pole of \(X\) forces \(v\sim2X^2\) from the
second equation of (5), and the quartic term in the
first then dominates. With \(X\) integral the first
equation also bounds \(v\). Together with the removable
profile and the finite exclusions, this proves that all
normalized coefficients are regular at every finite
place. Then \(G\) in (2) is regular there by the Laurent
coefficient expressions, so \(G'=J\rho\) is regular at
unramified places. More directly, a pole of \(\rho\)
at a finite place would make the physical gaps diverge
unless a normalized pair coalesced. To rule out such a
pole even at a ramified place, use the following descent.

Before normalization every coefficient of weight \(j\)
is rational after multiplication by \(\rho^j\).
Since \(\sigma\ne0\) has weight ten and \(\rho^{15}\)
is rational, \(\rho^5\in\mathcal K\). If any of
\(\tau,K,B\) is nonzero, its weight \(11,12,13\),
respectively, gives \(\rho^j\in\mathcal K\) with
\(\gcd(5,j)=1\), so \(\rho\in\mathcal K\). All places
are then unramified. If \(\rho\notin\mathcal K\),
the extension has degree five and \(\tau=K=B=0\).
A ramified place has index five. Inertia acts on the
normalized coefficients \(k,l,r,s,t\) with weights
\(2,3,4,5,6\), respectively. Their regularity therefore
gives residues \(k=l=r=t=0\), with \(s\ne0\) since
\(\sigma=-3s^2/4\) there. These residues give
\[
 W_0=z^6+Cz,\qquad V_0=z^9+\tfrac32Cz^4,\qquad C\ne0.
\]
One has
\(W_0'=6z^5+C\) and \(V_0'=9z^8+6Cz^3\).
They have no common zero: at \(z^5=-C/6\),
\(V_0'= (9C/2)z^3\ne0\).
The scaled endpoints are bounded by the degrees six,
eight and nine, the last extra polynomial being
\(V^2-W^3\) with leading coefficient \(\sigma\).
A pole of \(\rho\) would therefore require a coalescing
pair at a common critical point, which is impossible.
This excludes poles also at ramified finite places.

Finally a finite zero of \(\rho\) would collapse all
physical gaps, since the same degrees keep the normalized
pairs bounded at every finite place. Thus \(\rho^{15}\)
is a rational function with neither finite zeros nor
finite poles, hence \(\rho\) is constant.
Now the coefficients are rational in \(y\) and regular
at every finite place, so they are polynomials.
The nonconstant polynomial \(r\) has a pole at infinity.
The two pole profiles above give \(G\) pole order either
\(7p/2\) or \(6p\), each greater than one. This is
incompatible with \(G'=J\rho\in\mathbb C^\times\).

## Exact controls

The [checker](verify_conductor_twenty_sextic_nonic_exclusion.py)
verifies the parameter equations, the removable profile
and the second-order midpoint coefficient at the square
pole. The matching equations are the quadratic identities
from the [sextic–nonic conductor-eighteen proof](CONDUCTOR_EIGHTEEN_SEXTIC_NONIC_EXCLUSION.md#matching-equations-after-rescaling).
Replay from the repository root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_twenty_sextic_nonic_exclusion.py
~~~
