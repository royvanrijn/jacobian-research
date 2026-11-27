# The degree-nine cable exclusion at conductor fourteen

**Theorem ConductorFourteenDegreeNineCableExclusion.** Let \(F=(P,Q)\)
be an actual complex polynomial plane Keller map with
\(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\). Suppose \(P\) is monic of
positive degree in \(x\) and the global \(x\)-degree of the exact
coordinate conductor is fourteen. The reduced one-place characteristic
sequence of the generic coordinate order cannot be \((9,6,5)\).

The degree-five approximate root supplies a nonzero parameter
invariant. Its value makes the longitudinal scale rational.
The remaining parameter has neither a finite pole nor a finite zero,
as shown by the collision gaps in the original source coordinate.

## Actual generators and two first integrals

Put \(K=\mathbb C(y)\). The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives normalization \(K[x]\) and identifies the exact conductor
on each source line. For every finite \(y_0\) at least one generic
matched pair has finite distinct limits in \(x\), after local extension
to its two endpoints. A conductor-root branch approaches a finite root
of the nonconstant specialized conductor. Monicity of the original
\(P\) keeps its mate finite, and equal limits would force \(P_x=Q_x=0\)
by the divided differences of the two original generators. Thus at
every finite parameter some physical collision gap is a unit.

Reduce the actual generators by triangular target changes over \(K\).
For sequence \((9,6,5)\) the monic generators \(U,V\) have degrees six
and nine. The one-place theorem gives an approximate root monic cubic
in \(U\), of longitudinal degree five. The relation has the coefficient
degree bounds from normalization by the source line, as in the
[conductor-twelve degree-nine case](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#3-sequence-964).
The source theorem is Assi–García-Sánchez, Propositions 2 and 4.
[Retained primary source](sources/assi-garcia-sanchez-1407.0490v1.pdf)

For weights \((6,9)\), the exponents zero, one and two in \(U\)
have distinct weights modulo nine. The approximate root can therefore
have no terms of weight above eighteen, and its highest terms are
\(U^3-V^2\). Completing the square in \(V\) and translating \(U\)
to remove the \(U^2\) term gives, up to a constant,
\[
               W=V^2-U^3-\alpha U,\qquad \deg_z W=5.        \tag{1}
\]
Depress \(U\) in \(u=x-h(y)\), \(h\in K\). If the target
determinant of these changes is \(\lambda\in K^\times\), take
\(\rho^{15}=\lambda\), set \(u=\rho z\), and scale \(U,V\)
by \(\rho^{-6},\rho^{-9}\). The resulting target determinant
is one. In this frame the parameter derivative of the target change
is a polynomial Hamiltonian field, and gives the identity
\[
                U_zV_y-U_yV_z=J\rho+\partial_z H(U,V)       \tag{2}
\]
for a target polynomial \(H\) over \(K(\rho)\).
This follows by differentiating the target determinant as in the
[parameter equation for conductor eight](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).

Write \(S=(U^3+\alpha U)^{1/2}\) for the monic formal series at
\(z=\infty\), so \(S-V=O(z^{-4})\). In (2) the contribution of \(S\)
has polynomial part \([\alpha' U_z U^{-1/2}/2]_{\ge0}\).
The contribution of \(V-S\) has degree at most one, since \(U\)
is monic depressed of degree six. Therefore the left side of (2)
has degree at most two. The pole semigroup \(\langle5,6,9\rangle\)
has no positive element at most three, so \(H(U,V)\) is constant
in \(z\). At degree two we obtain \(\alpha'=0\). At degree one
we obtain that \([z^{-4}](S-V)\) is constant in \(y\), and
at degree zero we obtain \(J\rho=-6([z^{-5}](S-V))'\).

Write \(U=L^2+M\) with \(L=z^3+kz+\ell\) and \(M=rz^2+sz+t\).
The formal expansion through \(z^{-5}\) gives \(V=[S]_{\ge0}\)
and
\[
 S-V=
 \frac{3\,\operatorname{rem}(M^2,L)/8+\alpha/2}{L}
                  -\frac{M^3}{16L^3}+O(z^{-7}).
\]
The coefficients in degrees \(-1,-2,-3\) vanish. These
give \(s^2+2rt-r^2k=0\), \(2st-r^2\ell-2rsk=0\), and
\(3(t^2-2rs\ell)/8+\alpha/2=r^3/16\), respectively.
If \(r=0\), the first equation gives \(s=0\) and the tail
has zero coefficient at \(z^{-4}\), contrary to the degree-five
condition. Set \(T=t+2s^2/r\) and use the three equations to write
\[
\begin{aligned}
 t&=T-2s^2/r,&
 k&=2T/r-3s^2/r^2,&
 \ell&=-2sT/r^2+2s^3/r^3,\\
 \alpha&=r^3/8-3T^2/4,&
 V&=L^3+\tfrac32 LM+\tfrac38(r^2z+2rs).
\end{aligned}                                             \tag{3}
\]
With \(A=r^2s\) and \(X=r^2T-3rs^2\), the next terms are
\[
 [z^{-4}](S-V)=-3A/16,\qquad [z^{-5}](S-V)=X/16 .
\]
Thus
\[
 A\in\mathbb C^\times,\qquad \alpha\in\mathbb C,\qquad
                    J\rho=-3X'/8 .                       \tag{4}
\]
The [exact checker](verify_conductor_fourteen_degree_nine_cable_exclusion.py)
replays these comparisons, including the residual bracket identity.

The scaling weight of \(r^2s\) is thirteen, since \(r,s\)
have weights four and five. The nonzero constant \(A\) gives
\(\rho^{13}\in K\), and \(\rho^{15}\in K\) by construction.
Thus \(\rho\in K\), and the coefficients in (3) are rational in \(y\).
In particular
\[
                   s=A/r^2,\qquad X=r^2T-3A^2/r^3 .       \tag{5}
\]

## A finite pole of \(r\)

At a pole of \(r\), equations (3) and (5) give
\(T^2\sim r^3/6\), \(k\sim 2T/r\), and \(r/k^2\to3/2\).
The terms \(s\) and \(\ell\) tend to zero. Scale \(z=\sqrt{k}\,w\)
and the three functions \(U,V,W\) by \(k^{-3}, k^{-9/2}, k^{-9}\),
respectively. Their limits are the functions in (3) and (1) at
\(k=1,\ell=s=0,r=3/2,t=3/4,\alpha=0\). The limit of \(W\)
has degree four, and that of \(V\) has degree nine.

For a bounded scaled pair, coalescence would require a common
critical point of the three limiting functions. The only zeros
of \(W_w\) are \(w=0\) and \(w^2=-13/16\). At the first
\(V_w(0)=63/32\). At the latter, \(U_w/(2w)=315/256\).
Thus any bounded scaled matching pair has a unit gap in \(w\).
If one endpoint escapes in \(w\), the scaled \(V\) remains
monic of degree nine with bounded lower coefficients, so the other
endpoint escapes with the same valuation. Scale both by their common
pole. The divided difference of the leading ninth powers is nonzero
on the nonzero diagonal. Consequently the two leading endpoints are
distinct and the gap in \(w\) has a pole. This argument allows the
degree-five coefficient of the scaled \(W\) to vanish in the limit.

In either case the gap in \(z\) has a pole because it equals
\(\sqrt{k}\) times the gap in \(w\). Also \(X\sim r^2T\) has a pole,
so (4) gives a pole of \(\rho\). All physical gaps diverge.
This contradicts the finite unit-gap condition.

## A finite zero of \(r\)

Normalize the local valuation by \(v(y-y_0)=1\) and put
\(q=v(r)\ge1\). Further local extensions may give fractional valuations.
By (5), \(X\sim-3A^2/r^3\), so \(v(\rho)=-3q-1\).
Let
\[
 c=s/r=A/r^3,\qquad \epsilon=r/c^4,\qquad \mu=T/(rc^2).
\]
Then \(v(\epsilon)=13q\). If \(\alpha\ne0\) we have \(v(\mu)=5q\),
and if \(\alpha=0\) we have \(v(\mu)=13q/2\), by the equation for
\(\alpha\) in (3). In \(w=z/c\) write
\[
 L_\mu=w^3+(2\mu-3)w+2(1-\mu),\qquad N_\mu=w^2+w+\mu-2 .
\]
The functions \(\widetilde U=c^{-6}U\), \(\widetilde V=c^{-9}V\),
\(\widetilde W=c^{-18}\epsilon^{-3}W\) are
\[
\begin{aligned}
 \widetilde U &=L_\mu^2+\epsilon N_\mu,\\
 \widetilde V &=L_\mu^3+\tfrac32\epsilon L_\mu N_\mu
                         +\tfrac38\epsilon^2(w+2),\\
 \widetilde W &=(N_\mu^3-3\mu^2N_\mu-L_\mu^2)/8
                         +\epsilon(9(w+2)^2/64-N_\mu/8).
\end{aligned}                                             \tag{6}
\]
For the last identity, \(N_\mu^2-(w+2)L_\mu=\mu^2\) and
\(\alpha/c^{12}=\epsilon^3/8-3\epsilon^2\mu^2/4\). In the limit
\[
 L_0=(w-1)^2(w+2),\qquad
 \widetilde V_0=L_0^3,\qquad
 \widetilde W_0=\tfrac38(w-1)^3(w+2)^2 .
\]
The last two degrees are nine and five, with their leading
coefficients units throughout this limit. Thus every matched \(w\)-pair
is bounded. Otherwise, scaling any escaping pair by its largest
coordinate would give a common nonzero zero of the divided differences
of the two monic powers. Coprimality excludes this off the diagonal,
and their derivatives exclude it on the nonzero diagonal.
A coalescing pair can meet
only at \(w=1\) or \(w=-2\), the common zeros of the two limiting
derivatives. The physical scale satisfies
\(v(\rho c)=-6q-1\). In particular a pair with distinct limiting
\(w\)-coordinates diverges in the original source.

### Collapse at \(w=1\)

Set \(a=w-1\) and choose \(e^3=\epsilon\), so \(v(e)=13q/3\)
and \(v(\mu)>v(e)\). Here
\[
 L_\mu=a(a^2+3a+2\mu),\qquad N_\mu=a^2+3a+\mu .
\]
If some matched endpoint had \(v(a)<v(e)\), take an endpoint
coordinate \(a_*\) of least valuation in its pair and set
\(a=a_*z_*\). After dividing \(\widetilde U,\widetilde W\) by
\(a_*^4,a_*^3\), their limits would be \(9z_*^4,27z_*^3/8\).
Indeed \(v(\epsilon/a_*^3)>0\) and \(v(\mu/a_*)>0\).
One endpoint has \(z_*=1\), so the two monomial divided
differences exclude a matching pair at that scale.

Hence both \(a/e\) are bounded. For \(a=e z_*\) the functions
\(\widetilde U/e^4,\widetilde W/e^3\) have limits
\[
              9z_*^4+3z_*,\qquad 27z_*^3/8+81/64 .
\]
They have no common critical point. Thus the normalized \(w\)-gap
has valuation \(13q/3\). Its physical valuation
\(13q/3-6q-1\) is negative.

### Collapse at \(w=-2\)

Use midpoint \(m\) near \(-2\) and half-gap \(g\). The divided
difference of \(\widetilde U\) is polynomial in \(m,g^2\),
and its derivative in \(m\) at \(m=-2,\epsilon=\mu=g=0\)
is \(\widetilde U_0''(-2)=162\). Solve for \(m\) as a formal
series in \(\mu,\epsilon,g^2\). At \(g=0\) this gives the critical
point \(\eta\) of \(\widetilde U\) near \(-2\), with
\[
 \eta=-2+\tfrac23\mu+O(\mu^2,\epsilon),\qquad
 N_\mu(\eta)=-\mu+O(\mu^2,\epsilon).
\]
Here the terms linear in \(\mu\) follow from the root of \(L_\mu\)
near \(-2\), since at \(\epsilon=0\) that root is the critical point.
At \(\eta\) the identity
\[
       \widetilde V'-\tfrac32 L_\mu\widetilde U'
                =\tfrac38\epsilon(4L_\mu' N_\mu+\epsilon)
\]
gives leading term \(-27\epsilon\mu/2\) for
\(\widetilde V'(\eta)\), using \(v(\mu)<v(\epsilon)\).
The divided difference of \(\widetilde V\), after solving for
the midpoint, has coefficient of \(g^2\) specializing to
\(\widetilde V_0'''(-2)/6=729\). The term from the midpoint
derivative vanishes since \(\widetilde V_0''(-2)=0\).
The residue of that coefficient is a unit. Consequently any
collapsing matching pair here has
\[
                    v(g)=(v(\epsilon)+v(\mu))/2 .
\]
The two possibilities are \(9q\) and \(39q/4\), both strictly
larger than \(6q+1\) for \(q\ge1\). The corresponding
physical gaps tend to zero. The explicit finite-limit identities
and these first jets are covered by the exact checker.

No pair at a finite zero of \(r\) has a unit physical gap. Therefore
\(r\in\mathbb C(y)^\times\) has neither a finite zero nor a finite
pole and is a nonzero constant. Then \(s\) is constant by (5)
and \(T\) is constant by (3). This makes \(X'=0\), contradicting
(4), and proves the theorem.
