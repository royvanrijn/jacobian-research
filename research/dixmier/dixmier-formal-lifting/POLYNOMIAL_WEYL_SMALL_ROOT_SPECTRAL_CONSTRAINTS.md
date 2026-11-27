# Spectral constraints at small roots of Weyl leading factors

Let \(P,Q\in\mathbb C[x][D_x]\) satisfy \([Q,P]=1\), of orders
\(ga,gb\) with \(\gcd(a,b)=1\). Write the leading coefficients as
\(p=\alpha h^a,\ q=\beta h^b\). A root \(x_0\) of \(h\) is small
if its multiplicity \(u\) satisfies \(0<u<g\). Set \(k=g-u\),
\(m=ak\), \(n=ag\) for \(P\). For \(Q\) use \(bk,bg\).
All roots are geometric.

## Exact subprincipal coefficients at small roots

**Theorem PolynomialWeylSmallRootIndicialMean.** At each small root
both spectral equations are regular singular for every finite spectral
value. Their lowest Euler polynomials have degrees \(ga,gb\) and
root means \((ak-1)/2,(bk-1)/2\). For the next differential
coefficients \(p_1,q_1\) this gives
\[
 \operatorname{ord}_{x_0}\left(p_1-\frac{ga}{2}p'\right)\ge au,\qquad
 \operatorname{ord}_{x_0}\left(q_1-\frac{gb}{2}q'\right)\ge bu. \tag{1}
\]
The fixed local monodromy spaces consist of the holomorphic germs
and have dimensions \(ak,bk\).

Translate \(x_0\) to zero. At weight \((-l,r)\) with coprime
\(0<l<r\), any maximizing point \((i,j)\) of \(P\) has
\(j-i\ge n-au=m>0\), by comparison with \((au,n)\) and \(j\le n\).
The weight is positive. [Guccione–Guccione–Valqui, Theorem
4.1(4)](sources/guccione-guccione-valqui-1111.6100v3.pdf) supplies
a polynomial bracket mate of weight \(r-l\). The lattice forces
the form \(xy\,G(x^r y^l)\). Write the face as
\(x^iy^jL(x^r y^l)\), with \((i,j)\) its lower endpoint and
\(L(0)\ne0\). If \(\deg G=b_0>0\), the top coefficient in the
bracket equation is proportional to
\[
         (i+r\deg L)-(j+l\deg L)+(li-rj)b_0.
\]
Both terms are strictly negative. Thus \(G\) is constant,
and the bracket equation then forces \(L\) constant because
\(r-l>0\). Every rational slope in \((0,1)\) has a unique
maximizer. At small positive slope it is \((au,n)\), so it
stays there to slope one. This gives \(j-i\le m\) on the support.
Likewise \(Q\) satisfies \(j-i\le bk\).

With \(E=xD_x\), the equation becomes
\[
        x^m(P-\lambda)=f(E)+\sum_{j>0}x^j f_j(E)-\lambda x^m ,
        \qquad \deg f=n,\quad (E)_m\mid f . \tag{2}
\]
The sum is finite and \(\deg f_j\le n\). The leading coefficient in
\(E\) is a unit at zero, so (2) is regular singular. The lowest
Euler components \(P_0=x^{-ak}f(E)\), \(Q_0=x^{-bk}g_0(E)\) commute,
and \(P_0^b,Q_0^a\) are proportional by the [rational period
argument](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#the-spectral-monodromy-in-arbitrary-orders).

For the moving calculation set \(\lambda=\tau^{-m},\ x=z\tau\),
and scale \(P,Q\) by \(\tau^m,\tau^{bk}\).
The spatial companion is analytic in \(\tau\) along fixed ordinary
paths in \(\mathbb C^\times_z\). The leading coefficient in
(2) has nonzero constant term at \(\tau=0\), so its inverse
expands over \(\mathbb C(z)[[\tau]]\). These are the hypotheses
used in the [moving spectral proof](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md#moving-spatial-point),
along with the limit \(f(zD_z)y=z^m y\) and the simultaneous
deck rotation preserving \(x,\lambda\). The same spatial
diagonalization at formal \(z=\infty\), \(\tau\)-adic recursion,
and finite polar descent give a power series unit spectral gauge
at \(z=1\) with rank-one factors of one common residue \(\mu\).
The finitely many descent and rotation paths avoid the other
roots for small enough \(\tau\).

The fixed-point Euler-jet matrix \(H(x,\lambda)\) for
\(\partial_\lambda\) is polynomial in \(\lambda\) and rational
in \(x\). The moving matrix has trace
\[
 -m\tau^{-m-1}\operatorname{tr}H(z\tau,\tau^{-m})
          +\tau^{-1}\operatorname{tr}A(z,\tau),
\]
where \(A\) is the spatial Euler matrix. The first term
contributes to the residue only from \(x^{m(j+1)}\lambda^j\)
in the local Laurent expansions, \(j\ge0\). Its contribution
is a polynomial in strictly positive powers of \(z\).
The second contributes the sum of the roots of \(f\).
The trace residue is constant in \(z\) by flatness, since
\(A\) is regular in \(\tau\). Hence the positive-power
polynomial vanishes, and \(\mu\) is the root mean of \(f\).

The [polynomial dual module argument](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#an-algebraic-model-for-the-entire-solution-connection)
works locally. With \(\lambda=P^\vee\), the polynomial space
\(\mathbb C[t]\) is free of rank \(m\), since \(P^\vee\) raises
degree by exactly \(m\) above a finite cutoff, and the scalar
commutator excludes polynomial eigenvectors by producing
unbounded generalized eigenvector dimension. Its dual
fibres are the Taylor solution germs. Their series converge
by (2), since the finitely many lag ratios are bounded above
the last nonnegative integer root of \(f\).
The free-module coordinates of \(t^j\) have \(\lambda\)-degree
at most \((j+C)/m\) by degree reduction. After multiplication
by a common power of \(\tau\), the scaled Taylor series have
polynomial \(\tau\)-coefficients and converge normally near
\((0,0)\) by the uniform finite-lag bound. Spatial transport
then gives meromorphic moving evaluation at \(z=1\).

Choose the last nonnegative integer root \(\nu_r\) of \(f\)
in each class \(r\bmod m\), \(0\le r<m\). Their Taylor
coefficients are coordinates on the limiting holomorphic
solutions. The free polynomial submodule generated by
\(t^{\nu_r}\) has quotient dimension
\(L=\sum_r(\nu_r-r)/m\), by the distinct leading degrees
\(\nu_r+jm\) of the \(P^\vee\)-iterates.
The determinant of the selected Taylor derivatives of
the evaluated dual basis is
\(\tau^{\sum_r\nu_r}\Delta(\tau^{-m})\) with \(\deg\Delta=L\).
For the [exterior-product comparison](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md#exact-indicial-means-and-the-mixed-weyl-coefficient),
normalize by the selected derivatives. Any pole of the
normalized columns would give a leading limiting solution
with zero selected derivatives. They therefore specialize
to independent solutions, also after transport to \(z=1\).
The exterior product of the original columns has valuation
\(\sum_r r=m(m-1)/2\). In the unit spectral gauge every
nonzero component of this exterior product has logarithmic
derivative of residue \(m\mu\), since the polynomial dual
connection has zero trace residue after pullback. Thus
\(\mu=(m-1)/2\). The two highest coefficients of \(f\)
give (1).

The limiting fixed monodromy solutions are meromorphic
by regular singularity. In each integer class the recurrence
forces all coefficients before the last root to vanish,
and that root is nonnegative. The limiting fixed rank is
therefore \(m\). The local monodromy is conjugate for every
finite \(\lambda\), by flatness of \(\partial_\lambda+Q\),
and the scaled monodromy converges to the limiting one.
Its fixed rank is at most \(m\) by specialization and at
least \(m\) by the Taylor module. This proves the fixed
germ assertion.

## Common spectral phases and the opposite Fourier coefficient

**Corollary PolynomialWeylSmallRootGlobalPhase.** Let
\(e=\gcd(g,\{u_i\})\), using the multiplicities of all roots of
nonconstant \(h\), and suppose there is at least one small root.
The number \(r_P\) of distinct higher phases is the same in every
small-root moving chart, and \(a\mid r_P\mid ae\).
Likewise \(b\mid r_Q\mid be\). Each small root imposes
\(r_P\mid ac_i,\ r_Q\mid bc_i\), where \(c_i\) is the
translation stabilizer order of its primitive indicial multiset
modulo \(g-u_i\).

At fixed ordinary \(x\), use \(Z^n=\lambda\).
Spatial diagonalization over the cover containing \(p^{1/n}\)
has leading eigenvalues \(\omega_i Zp^{-1/n}\). Flatness
makes the phase coefficients above the linear \(Z\)-term constant
in \(x\). The [leading-factor residue theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
gives a polynomial \(B\) satisfying \(B'-(p'/np)B=1\).
The linear phase coefficient is \(\omega_i B/p^{1/n}\) since
this has the right derivative and spatial deck character.
The difference of the two primitives is constant, and must
vanish because the spatial deck on \(p^{1/n}\) has order
\(g/e>1\). This deck fixes the higher coefficients in \(x\)
and rotates the labels by order \(g/e\). The higher phases
are therefore \(F(\omega_i Z)\) with powers in \(F\)
multiples of \(g/e\).

To compare with the moving chart at any small root, put
\(H=x^{k/g}Z\) in the translated coordinate. The
\(H^{-j}E^j\) spatial jets have a leading matrix \(H\) times
a companion with simple eigenvalues at \(x=0\). The eigenbasis
and diagonal recursion lie over \(\mathbb C[[x]][[H^{-1}]]\)
with invertible constant term. The \(Z\partial_Z\) matrix
has a uniform lower bound on the \(x\)-powers before gauge,
since its finite reduction in (2) divides by a unit leading
coefficient, and the gauge preserves this bound. Therefore
\(x=z\tau,\ H=z^{k/g}\) is defined in the double Laurent
expansion at formal spatial infinity, as in the [joint
phase proof](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md#the-same-phase-polynomial-controls-both-presentations).
The term \(ZB/p^{1/n}\) has only nonnegative \(\tau\)-powers
after substitution, since \(B\) has a simple zero at a small
root. Thus the moving irregular phases are exactly
\(F(\omega_i\tau^{-k/g})\), giving the same count at every
small root. The leading phase has degree \(ga+gb\) with
nonzero coefficient from the top term of \(Q\), so the deck
orbit on leading types has size \(a\). Spatial invariance
gives \(r_P\mid ae\). The sharper \(ac_i\) bound follows
by applying the [rational spatial-projector calculation](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md#indicial-classes-bound-algebraic-spatial-refinements)
to these two Euler limits, with projectors supplied by the
finite polar descent over \(\mathbb C(z)\) just as above.

**Corollary PolynomialWeylRegularRootSubprincipalSupport.** Define
\[
        S_P=\frac{p_1-(ga/2)p'}{ga\,p},\qquad
        S_Q=\frac{q_1-(gb/2)q'}{gb\,q}.
\]
In the presence of a small root the two rational functions are
equal and have poles only at roots of multiplicity greater than
\(g\). For the opposite Fourier presentation of a centered
pair with zero symmetric subprincipal symbols, they satisfy
\(S_P=S_Q=O(x^{-2})\) at infinity. If this Fourier leading
factor has only one high root and at least one small root,
the residue at the only possible pole is zero.

In symmetric ordering the next differential coefficient
of the commutator after the leading bracket is the sum
of the two Poisson terms involving the subprincipal symbols.
The cubic contraction has lower differential order.
Using \(p'/((ga)p)=q'/((gb)q)\) gives
\[
        (S_Q-S_P)'+\frac{p'}{ga\,p}(S_Q-S_P)=0.
\]
A nonzero solution is a constant times \(p^{-1/(ga)}\),
which is not rational at a small root. Both functions are
regular there by (1), and multiplicity \(g\) is excluded by
the leading-factor residue equation. For the Fourier claim,
the original highest-momentum coefficient is a monomial and
the next symmetric row vanishes. After swapping position
and momentum, the next symmetric differential coefficient
has degree at least two below the leading coefficient.
The one-pole assertion follows from the residue at infinity.


## Indicial constraints for each monodromy eigenvalue

**Theorem PolynomialWeylSmallRootSectorCongruence.** Let
\(\mathcal R\) be the common primitive indicial multiset modulo
\(k=g-u\), so the multisets of \(f,g_0\) modulo \(k\) are
\(a\mathcal R,b\mathcal R\). For each class \(\eta\bmod\mathbb Z\)
let \(\mathcal R_\eta\) be its submultiset and \(d_\eta\)
its size. Then
\[
                   \sum_{\rho\in\mathcal R_\eta}\rho
                    \equiv d_\eta(k-1)/2\pmod{k}.            \tag{3}
\]
In particular the finite spatial monodromy at a small root is
quasi-unipotent. For a nonidentity eigenvalue \(\zeta\) with
primitive multiplicity \(d_\eta\le u\) one has
\(\zeta^{d_\eta}=(-1)^{d_\eta(k-1)}\).
Every singleton noninteger class is \((k-1)/2\bmod k\).

The root-divisor equation from \([P_0,Q_0]=0\) makes
the \(f\)-multiset modulo \(ak\) invariant under translation
by \(bk\) and hence by \(k\). The \(g_0\)-multiset modulo
\(bk\) is invariant under \(k\) as well.
From \(P_0^b=\mathrm{const}\,Q_0^a\), the multisets modulo
\(k\) are \(a\mathcal R,b\mathcal R\) for one multiset
\(\mathcal R\) of size \(g\). In particular the \(f\)-classes
modulo \(m=ak\) are all \(a\) lifts of the classes in
\(\mathcal R\). The \(k\) integer classes occur in
\(\mathcal R\) by \((E)_{ak}\mid f\).

For a representative \(\eta\), let \(V\) be the generalized
eigenspace for \(\exp(2\pi i\eta)\) of the monodromy around
the small root. Its rank is \(D=a d_\eta\). Its germs have
Frobenius expansions
\[
       x^\eta\sum_\ell x^\ell U_\ell(\log x),
\]
with the logarithmic degrees bounded by \(n-1\).
In the recurrence from (2) the operator on \(U_\ell\) is
\(f(\eta+\ell+\partial_{\log x})\), triangular with diagonal
\(f(\eta+\ell)\). Only finitely many indices are exceptional.
From a cutoff before the first root to one above the last,
a finite polynomial matrix in \(\lambda\) determines such
germs. Above the cutoff the same recurrence determines
the tails, with bounded lag ratios on compact \(\lambda\)-sets:
the logarithmic derivatives of \(f_j\) only lower their degrees.
Thus the formal germs converge. Constant monodromy type
makes this kernel locally free over the whole spectral line,
so it is free over \(\mathbb C[\lambda]\). The operators
\(\partial_{\log x}\) and \(\partial_\lambda+Q\) preserve
the space and act in its finite coordinates by polynomial
operations. In particular the first is nilpotent and the
connection has polynomial trace. The same finite-coordinate
description shows that every constant-rank image of
\(\partial_{\log x}\) is a polynomial subconnection.

Use a free basis of \(V\). At \(x=z\tau,\lambda=\tau^{-m}\),
replace \(\log x\) by \(\log z\) and remove the common factor
\(\tau^\eta\). This normalized moving evaluation is meromorphic
in \(\tau\): the recurrence gives
\(\deg_\lambda U_\ell\le(\ell+C)/m\), and each logarithmic
coefficient has the uniform finite-lag bound near \(z=0\).
Along ordinary paths in \(z\) the evaluation extends by spatial
transport. Replacing \(\log x\) by \(\log z\) amounts to a
basis change by the exponential of \(-\log\tau\) times the
nilpotent logarithmic derivative. Its determinant is one.
Consequently the exterior product of the normalized frame
obeys the moving exterior equation with spectral trace
subtracted and with the extra term \(-D\eta/\tau\).
In the unit spectral gauge its valuation is
\(j=D(\mu-\eta)\), since each nonzero component obeys a
scalar logarithmic derivative equation with this residue.

The leading exterior coefficient spans the determinant of
the limiting generalized eigenspace. Indeed the scaled
monodromy is analytic in \(\tau\) with constant characteristic
polynomial, and its generalized eigenspace projectors are
polynomials in the monodromy, obtained by separating the
distinct eigenvalues with their powers bounded by \(n\).
They extend at constant rank through \(\tau=0\).
Let \(\xi=\exp(2\pi i/m)\). Under
\((z,\tau)\mapsto(\xi z,\xi^{-1}\tau)\) the normalized exterior
product acquires the factor \(\exp(2\pi i D\eta/m)\):
the other factor from the logarithm shift again has determinant
one. The leading exterior vector therefore transforms under
\(z\mapsto\xi z\) by \(\exp(2\pi i(j+D\eta)/m)\).

In the limiting equation \(f(zD_z)y=z^m y\) the rotation
is one turn of the regular singular equation in \(t=z^m\),
namely \(f(m tD_t)y=t y\). Its eigenvalues come from the
indicial roots of \(f\) modulo \(m\). On the selected
generalized eigenspace those roots are
\(\rho+ik\), for \(\rho\in\mathcal R_\eta\) and \(0\le i<a\).
The determinant of rotation there is consequently
\[
 \exp\!\left(\frac{2\pi i}{m}
     \left(a\sum_{\rho\in\mathcal R_\eta}\rho+
                       d_\eta k a(a-1)/2\right)\right).
\]
Comparing the two determinants and using \(\mu=(ak-1)/2\)
gives (3). Reduction modulo \(\mathbb Z\) gives the stated
root-of-unity restriction.

## Monodromy at simple and double leading roots

**Corollary PolynomialWeylSimpleRootMonodromy.** For \(u=1\)
one has \(k=g-1\), and \(\mathcal R\) consists of
\(0,\ldots,k-1\) and the class \((k-1)/2\bmod k\).
At every finite spectral value the \(P\)-monodromy is
unipotent for even \(g\), with \(a\) size-two blocks and
\(a(k-1)\) size-one blocks. For odd \(g\) it is semisimple,
with multiplicities \(ak\) and \(a\) at \(1\) and \(-1\).
The \(Q\)-formulas replace \(a\) by \(b\).
For a small double root, the two extra primitive classes
are both integral or both half-integral.

Summing (3) gives that the sum of all extra classes is
\(u(k-1)/2\bmod k\). This pins the one extra class for \(u=1\).
The multiplicity of each Jordan size at each eigenvalue
is divisible by \(a\) in the actual \(P\)-monodromy.
The [leading moving block
argument](POLYNOMIAL_WEYL_ALL_DEGREE_SPECTRAL_REDUCTION.md#the-spectral-monodromy-in-arbitrary-orders)
applies locally: the leading matrix of the \(s\)-connection,
\(s=\tau^{-1}\), is \(-ak[Q_0]_{P_0-1}\), semisimple by the
power identity. The deck cycles its \(a\) leading
eigenvalues. Local spatial monodromy is horizontal
and analytic in \(\tau\), preserves these blocks, and
has the same Jordan data on all of them. For even \(g\),
all eigenvalues are 1 and the nontrivial Jordan contribution
is \(ag-ak=a\), so there must be \(a\) size-two blocks.
For odd \(g\), the \(-1\) multiplicity is \(a\), so it
has only size-one blocks, and the fixed rank exhausts
the 1-primary multiplicity.

At a double root, the sum of the two extras is integral
modulo \(\mathbb Z\). If there is a singleton noninteger
class it must be half-integral by (3), and the other
root would have the same noninteger class. Any pair
in one noninteger class has half-integral class since
its sum is integral. This proves the double-root statement.
For comparison, the commuting Euler control with
\(g=6,k=4\) and primitive roots
\(0,1,2,3,1/3,8/3\) has the exact whole indicial mean,
but fails (3) in both singleton noninteger sectors.
The root choices \(0,1,2,3,1,2\) and
\(0,1,2,3,3/2,3/2\) pass this test.

## The four-simple-root order-twelve profile

**Corollary PolynomialWeylOrderTwelveSimpleRootSpectralProfile.**
In the reduced \((12,18)\) boundary profile
\(h=A\,t^{26}(t^4-\lambda_0)\) with \(A\lambda_0\ne0\),
every simple root gives extra indicial classes \(2,7\bmod10\)
for \(P\) and \(2,7,12\bmod15\) for \(Q\), in addition to
one copy of each integer class. The local monodromies have
two and three unipotent size-two blocks, respectively.
The moving phase counts are two and three. With
\(p=\alpha h^2,\ q=\beta h^3\), the subprincipal numerators
satisfy
\[
 (t^4-\lambda_0)^2\mid p_1-6p',\qquad
 (t^4-\lambda_0)^3\mid q_1-9q'.
\]
The equal normalized subprincipal functions can have finite
poles only at the multiplicity-26 root.

The [leading factor formula](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-degree-bounds-and-a-real-refinement)
gives the four simple roots and one high root. Here \(g=6\),
\(k=5\), \(a=2\), and \(b=3\). Apply the simple-root conclusion
with \((k-1)/2=2\), the phase bound with \(e=1\), and (1)
at all four roots. No existence of a Weyl pair with this
profile is implied.


## Fractional powers in the differential mate

**Corollary PolynomialWeylSmallRootFractionalPowerMate.** In
the presence of a small root, put \(e=\gcd(g,\{u_i\})\).
Let \(L=P^{1/(ag)}\) be a formal pseudodifferential root
with a chosen leading coefficient \(p^{1/(ag)}\), and let
\((\cdot)_+\) retain the nonnegative differential powers.
Then
\[
              Q=\sum_{j=0}^{be}\gamma_j
                    \bigl(P^{j/(ae)}\bigr)_+,\qquad
              \gamma_j\in\mathbb C,\quad \gamma_{be}\ne0, \tag{4}
\]
where \(P^{j/(ae)}=L^{jg/e}\).
In particular, a simple root gives \(e=1\). For the
order-\((12,18)\) four-simple-root profile, after
subtracting a scalar and a multiple of \(P\) from \(Q\),
one has
\[
            Q=\gamma_3(P^{3/2})_+ +\gamma_1(P^{1/2})_+,
            \qquad \gamma_3\ne0 .
\]

Use a scalar formal wave in one of the fixed ordinary
spatial lines for \(Z^n=\lambda\), \(n=ga\).
Its leading spatial logarithmic derivative is
\(\omega_i p^{-1/n}Z\). The pseudodifferential operator
\(L\) acts on this wave, since \(D_x^{-1}\) can be
expanded by division by that nonzero leading term.
It preserves the spatial solution equation and the
same formal line. Its eigenvalue there is \(\omega_i Z\):
the leading term has that value, and its \(n\)-th power
is \(Z^n\). Uniqueness of the spatial line follows from
the simple leading spatial eigenvalues.

The joint formal wave can be chosen horizontal for
\(\partial_\lambda+Q\). Its exponent is
\(F(\omega_i Z)+\omega_i Bp^{-1/n}Z\), up to the
logarithmic term. The logarithmic derivative in \(Z\)
of its remaining Laurent prefactor has \(Z\partial_Z\)
degree at most zero. Thus the nonnegative terms of
\(Qy/y\) are precisely those of
\(-Z^{-n}(Z\partial_Z)F(\omega_i Z)/n\).
All their powers of \(Z\) are divisible by \(g/e\).
Subtract the corresponding constant-coefficient sum
of powers of \(L\) from \(Q\). If the difference had
nonnegative differential order, its leading symbol
would act on this wave with a nonzero nonnegative
power of \(Z\). All such powers have been canceled.
Taking differential parts gives (4). In this formula
\(P^{1/(ae)}\) itself has rational pseudodifferential
coefficients: its leading coefficient is a constant
times the polynomial \(h^{1/e}\), and the root recursion
uses rational functions of \(x\). The fixed-wave comparison
is the local version of the diagonal calculation in
[Beals and Sattinger, Lemmas 8.2 and 9.3](https://arxiv.org/abs/solv-int/9801010).


## Exact indicial polynomials for ratio two to odd

**Corollary PolynomialWeylSimpleRootTwoOddPolynomials.**
At a simple common leading root for orders \(2g,bg\), with
\(b\ge3\) odd, let \(k=g-1\) and \(\theta=(k-1)/2\).
Besides the falling-factorial roots \(0,\ldots,2k-1\)
of \(f\), its only two extra roots are
\[
                  \theta-jk,\quad \theta+(j+1)k,\qquad
                         0\le j\le (b-1)/2.                 \tag{5}
\]
For the same \(j\), the extra roots of \(g_0\) are the
\(\theta+ik\) indexed by the nonzero coefficients of
\[
            (1+T+\cdots+T^{b-1})
                     \sum_{i=-j}^{j}(-1)^{i+j}T^i.        \tag{6}
\]
Each coefficient in (6) is zero or one. In the order-\((12,18)\)
boundary case the extra lists are \(2,7\) and \(2,7,12\),
or \(-3,12\) and \(-3,7,17\), at each simple root.

Let \(F_0,G_0\) be the root divisors of \(f,g_0\), and \(S\)
translate a root by \(k\). Define \(A_b=1+S+\cdots+S^{b-1}\).
Commutation gives \(A_b F_0=(1+S)G_0\), by canceling
\(1-S\) on finite divisors. Write \(A_b=1+(1+S)H\),
using \(H=S+S^3+\cdots+S^{b-2}\). Then \(R_0=G_0-HF_0\)
satisfies \(F_0=(1+S)R_0,\ G_0=A_b R_0\).
Subtract \(\sum_{i=0}^{k-1}[i]\) from \(R_0\) to get a
signed divisor \(R_1\). The divisors \((1+S)R_1\) and
\(A_bR_1\) are effective by the falling-factorial
divisibilities and have respective masses two and \(b\).
Their supports modulo \(k\) lie in the one extra class
\(\theta\). The same holds for \(R_1\) since \(1+S\) is
injective on finite signed divisors in each translation orbit.

On \(\theta+k\mathbb Z\), denote the shift by \(T\).
The two unit masses of \((1+T)R_1\) occur at
distinct positions of opposite parity, since its Laurent
polynomial vanishes at \(T=-1\). If those positions
are \(p<q\), then \(R_1\) is the alternating segment
starting at \(p\) of length \(q-p\), beginning and ending
with \(+1\). A length greater than \(b\) would give
coefficient \(-1\) at \(p+b\) in \(A_b R_1\).
For a length at most \(b\), every coefficient is the
sum of a prefix, a suffix, or the full alternating
segment, and is zero or one.
The exact mean of \(f\) gives first moment \(g(k-1)/2\)
for \(R_0\), so the first moment of \(R_1\) is \(\theta\).
The segment is therefore centered at zero in the
\(T\)-index and has odd length \(2j+1\le b\). This
proves (5) and (6). These are local necessary
conditions; the commuting Euler data do not realize
the scalar commutator.
