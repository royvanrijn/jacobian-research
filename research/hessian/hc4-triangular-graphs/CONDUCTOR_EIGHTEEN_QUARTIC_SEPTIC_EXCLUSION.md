# The quartic–septic conductor-eighteen exclusion

**Theorem ConductorEighteenQuarticSepticExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic of positive
degree in the longitudinal coordinate \(x\). If its exact directional
conductor has global degree eighteen, its reduced generic
characteristic sequence cannot be \((7,4)\).

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the normalization of the generic order and its specialized
conductors. As in the
[cubic–septic proof](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch),
at each finite parameter at least one generic matching pair has
distinct finite original limiting coordinates. A finite conductor
root gives such a branch, monicity of \(P\) keeps its mate finite,
and a collapsed pair would force \(P_x=Q_x=0\) by the two divided
differences. In particular some physical gap must be a unit at
every finite parameter.

## The parameter identities

Reduce the actual generators over \(K\) by triangular target
changes. In type \((7,4)\) they have degrees four and seven.
Make both monic, depress the quartic by centering
\(u=x-h(y)\), and remove its constant. Subtract a multiple
of the quartic and a constant from the septic to remove
coefficients in degrees four and zero. Let \(\lambda\in K^\times\)
be the determinant of the resulting target change.
Adjoin \(\rho^{11}=\lambda\), set \(u=\rho z\), and rescale
the target generators by \(\rho^{-4}\) and \(\rho^{-7}\).
The normalized pair has the form
\[
 U=z^4+az^2+bz,\qquad V=z^7+cz^6+dz^5+ez^3+fz^2+gz .
 \tag{1}
\]
The target determinant is now one. The derivative of the
parameter-dependent target change, composed with its inverse,
is a polynomial target field of divergence zero. Writing it
as a Hamiltonian field gives
\[
           U_z V_y-U_y V_z=J\rho+\partial_z H(U,V).         \tag{2}
\]
The Jacobian on the left has degree at most nine. The only
pole degrees through ten in \(\langle4,7\rangle\) are
\(0,4,7,8\), so \(H\) is a linear combination of
\(1,U,V,U^2\) after evaluation.

At degrees nine and eight in (2) one obtains \(c'=0\) and
\(d'=7a'/4\). Write \(d=7a/4+D\) with \(D\in\mathbb C\).
At degrees seven, six and three the coefficients of \(U^2,V,U\)
in \(H\) are forced, respectively, to be
\[
 \begin{split}
 \theta_{U^2}&=-(6ca'+7b')/8,\\
 \theta_V&=(16e'-20Da'-21aa'-24cb')/28,\\
 \theta_U&=(6a^2ca'+7a^2b'+4af'-4fa'+2be'-6eb')/8 .
 \end{split}
\]
The degree-five and degree-four comparisons integrate to
\[
\begin{aligned}
 f={}&F-\tfrac{15}{14}Dca+\tfrac54 Db-\tfrac{27}{16}ca^2
       -\tfrac7{16}ab-\tfrac97 c^2b+\tfrac67 ce,\\
 g={}&G-\tfrac{25}{28}D^2a-\tfrac54 Da^2-\tfrac{15}{14}Dcb
       +\tfrac57 De-\tfrac{35}{64}a^3-\tfrac{15}{8}cab
       +\tfrac34 ae-\tfrac{35}{32}b^2 ,
\end{aligned}                                             \tag{3}
\]
with \(F,G\in\mathbb C\).

To specify the remaining integrals without expanding all their
terms, use the polynomial differential \(\mathrm d\) in
\((a,b,e)\) with \(c,D,F,G\) held constant. In the three
\(\theta\)'s replace each prime by this differential, including
\(f'=\mathrm d f\) from (3). Define \(\omega_2,\omega_1,\omega_0\)
by the exact identity of polynomial one-forms
\[
 U_z\,\mathrm dV-V_z\,\mathrm dU
 -\theta_{U^2}(U^2)_z-\theta_V V_z-\theta_U U_z
                =\omega_2 z^2+\omega_1 z+\omega_0 .
\]
For a polynomial one-form \(\omega\) set
\[
 {\cal I}(\omega)
   =\int_0^1\omega_{(ta,tb,te)}(a,b,e)\,\mathrm dt.
\]
This is an algebraic prescription for a polynomial, obtained by
integrating each power of \(t\). Define
\(I_8={\cal I}(\omega_2)\), \(I_9={\cal I}(\omega_1)\), and
\(G_{10}={\cal I}(\omega_0+I_8\,\mathrm da/4)\).
The [exact checker](verify_conductor_eighteen_quartic_septic_exclusion.py)
verifies
\[
 \mathrm dI_8=\omega_2,\qquad \mathrm dI_9=\omega_1,\qquad
 \mathrm dG_{10}=\omega_0+I_8\,\mathrm da/4 .
\]
Hence on the parameter curve \(R_8=I_8\) and \(R_9=I_9\) are
constant, and
\[
           \Gamma=G_{10}-R_8a/4,\qquad \Gamma'=J\rho .      \tag{4}
\]
The weights of \(a,b,e\) are \(2,3,4\); those of the constants
\(c,D,F,G,R_8,R_9\) are \(1,2,5,6,8,9\). The integrals \(I_8,I_9\)
and \(\Gamma\) have weights \(8,9,10\), respectively. These are
the longitudinal scaling weights, so inertia acts by inverse
powers of its multiplier on \(\rho\).

## The highest-weight shapes

When the constants of positive weight are set to zero, the
equations \(I_8=I_9=0\) become
\[
 1225a^4-2688a^2e-1176ab^2+1536e^2=0,\qquad
 b(49a^3-48ae-14b^2)=0 .                                  \tag{5}
\]
The weight-ten polynomial for \(\Gamma\) in this limit is
\[
 \Gamma_*=
 -(833a^5-2240a^3e+2352a^2b^2+1536ae^2-2688b^2e)/7168 .
\]
Every nonzero shape in (5) has \(a\ne0\), since \(a=0\)
would make \(b=e=0\). By a constant scaling take \(a=1\).
The limiting pair is then
\[
 U_*=z^4+z^2+bz,\qquad
 V_*=z^7+\tfrac74z^5+ez^3-\tfrac7{16}bz^2+
             \bigl(\tfrac34e-\tfrac{35}{64}-\tfrac{35}{32}b^2\bigr)z .
\]
The two possible parameter loci and the reductions of \(\Gamma_*\)
and the longitudinal derivative resultant are as follows.
In the second row the substitution \(b^2=(49-48e)/14\)
has been made before reduction modulo the displayed quadratic.

| Branch | Equation for \(e\) | \(\Gamma_*\) | \(\operatorname{Res}_z(U_*',V_*')\) |
|---|---|---|---|
| \(b=0\) | \(1536e^2-2688e+1225=0\) | \(-(8e-7)/128\) | \(-245(48e-49)/128\) |
| \(b\ne0\) | \(1536e^2+1344e-2891=0\) | \(13(160e-161)/512\) | \(15925(12576e-12607)/128\) |

The resultant of each reduced expression with its defining
quadratic is nonzero, as the checker verifies. Thus
\(\Gamma_*\ne0\), and the two limiting generators in every
shape are immersed together.

## Finite poles and nonhomogeneous constants

Suppose some \(a,b,e\) has a pole at a finite place. Take the
positive scale \(q\) to be the maximum of their pole orders
divided by \(2,3,4\), respectively. After local extension, rescale
\(z\) by a factor of pole order \(q\). Rescale the monic generators
by its inverse fourth and seventh powers. In the limit at least
one scaled \(a,b,e\) is nonzero, and the constants have zero
residue after their respective scaling. The scaled coefficient
equations are precisely (5). Also \(\Gamma_*\ne0\), so
\(\Gamma\) has pole order \(10q\). By (4) the factor \(\rho\)
has a pole.

Every matched pair is bounded in the rescaled coordinate because
the limiting degrees four and seven are coprime and the two
leading coefficients stay units. Indeed, an escaping pair
would give a common nonzero zero of the divided differences
of the two monic powers after rescaling by the largest endpoint.
Coprimality rules this out off the diagonal, and their derivatives
rule it out on the nonzero diagonal. The limiting immersion
then prevents bounded matching endpoints from coalescing.
Thus every normalized gap in \(z\) has pole order \(q\).
All physical gaps have poles once multiplied by \(\rho\),
contradicting the finite-pair condition.

Consequently \(a,b,e\) and \(\Gamma\) are regular at all finite
places. Monicity in the coprime degrees again bounds all normalized
matching pairs at such places. A finite zero of \(\rho\) would
shrink every physical gap to zero. If any of
\(c,D,F,G,R_8,R_9\) is nonzero, its weight is coprime to eleven.
Its unscaled value lies in \(K\), so its constant scaled value
gives the corresponding power of \(\rho\) in \(K\).
Together with \(\rho^{11}\in K\) this gives \(\rho\in K\).
Regularity of \(\Gamma\) and (4) then exclude finite poles of
\(\rho\). Hence \(\rho\) is constant and \(a,b,e\) are polynomials
in \(y\). Any pole of one of these polynomials at infinity
satisfies the same initial equations (5) and gives
\(\deg_y\Gamma=10q>1\), since the pole of \(a\) determines \(2q\).
If all three are constant, so is \(\Gamma\). Both conclusions
contradict \(\Gamma'=J\rho\in\mathbb C^\times\).

## Homogeneous constants and the specialized tangent determinant

Suppose the six constants are all zero. Equations (5) hold
exactly, and \(a\ne0\), for otherwise (4) would fail.
Choose \(\tau^2=a\). After scaling by \(\tau\), the coefficients
belong to one of the finite shapes in the table and hence are
fixed complex constants. In the normalized coordinate \(z/\tau\)
all matching pairs are fixed. They are finite because the field
index divides both four and seven. The physical gaps are their
nonzero constant gaps multiplied by \(w=\rho\tau\).
The finite-pair condition makes \(w\) a unit at every finite
place. Since \(w^2\) is the unscaled quartic coefficient of
degree two, \(w^2\in K\). It follows that \(w\in\mathbb C^\times\).
Now \(\Gamma=\gamma_0\tau^{10}\) with \(\gamma_0\ne0\) by the
table. Equation (4) gives
\[
                 (\tau^{11})'=11Jw/(10\gamma_0).
\]
Thus the determinant \(\lambda^{-1}=w^{-11}\tau^{11}\)
has a finite zero \(y_0\).

In \(u=x-h(y)\) the unscaled pair is now fixed, since the
source scale \(w\) is constant. All its matching coordinates
are fixed in \(u\), so a finite pole of \(h\) would take every
physical pair to infinity. Hence \(h\) is polynomial. Write
\(R_*\subset\mathbb C[u]\) for the fixed canonical order and
\(C_*\subset\mathbb A^2\) for the image curve of the two
fixed generators, with irreducible equation \(L=0\).
The original \(P,Q\) in \(u,y\) belong coefficientwise to this
fixed order. In fact in each bounded \(u\)-degree, membership
in \(R_*\) imposes constant linear equations, so generic
membership gives the coefficientwise assertion. Choosing
target-polynomial representatives for the finitely many
coefficients gives a lift \(\Psi_y\) of the original pair
in \(\mathbb C[y][U_*,V_*]\).

At every finite parameter the actual specialized order lies
in \(R_*\). Its colength is at most nine by exact conductor
specialization, whereas \(R_*\) has colength nine. Hence the
two orders agree at \(y_0\). Choose a point \(p_*\) of \(C_*\)
given by a fixed matched pair of distinct normalization points.
Such a pair exists by the finite-pair property. The image
point is singular, because the normalization has two distinct
points over it, so \(\mathrm dL(p_*)=0\).

Let \(\Phi_y\) denote the generic inverse target change from
the fixed pair to the original pair. The components of
\(\Phi_y-\Psi_y\) are multiples of \(L\) over \(K\). Their
target derivatives at \(p_*\) vanish, so
\[
                   \det D\Psi_y(p_*)=\lambda^{-1}.
\]
At \(y_0\) the map of image curves induced by \(\Psi_{y_0}\)
is an isomorphism, by equality of the orders. The two
singular tangent spaces are both the full two-dimensional
ambient tangent spaces. Hence \(D\Psi_{y_0}(p_*)\) is invertible,
contrary to the zero of \(\lambda^{-1}\). This is the
[singular tangent determinant argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md)
applied to the fixed quartic–septic order. The contradiction
excludes the homogeneous case and proves the theorem.
