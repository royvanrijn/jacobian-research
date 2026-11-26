# A finite menu of diffuse extensions for arbitrary ordinary cores

Fix a positive rational interval \((a,b)\) and an integer \(J\ge4\)
at least the [effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
There is an explicit finite menu of rational strengths in \((a,b)\)
and there are computable positive tolerances \(\eta_*,\delta_*\),
depending only on the interval and \(J\), with the following property.

**Theorem.** For every finite ordinary two-monomial core, one can
select a strength \(\tau_*\) from the menu such that any appended
shared ordinary rows with curvature \(\tau\) and maximum relative
gap \(\eta\) satisfying
\[
                    |\tau-\tau_*|\le\delta_*,\qquad \eta\le\eta_*
\]
give disjoint consecutive root arguments in every degree.
Here an appended row \((u,v)\) has curvature
\((v-u)^2/\max(u,v)\) and relative gap \(|v-u|/\max(u,v)\).
The result allows an empty core. There is no bound on the core
dimension, private exponents, gap sum, or curvature. The strength
selection uses the core, while the menu and tolerances are uniform.
The interval can be arbitrarily close to zero.

Every menu strength can be realized exactly by ordinary shared
unit-gap rows with any prescribed integer signed total gap and
arbitrarily small relative gaps. If the prescribed gap is zero,
one may keep the raw ratio of leading and constant factorial
coefficients of the core fixed at every order by using mirrored
added rows. Thus the construction gives physical SFC(2) as well,
over every characteristic-zero field.

## A core-independent bound for the exceptional strengths

For \(2\le n<N=n+1<J\), normalize the two core polynomials to have
both endpoints one by independent positive scalars and changes of
variable. Denote these polynomials by \(P_j\), \(j=n,N\). Their
coefficients are at most \(\binom jk\) by the factorial chord bound.
Set \(q=e^{-\tau/(2nN)}\) and deform each coefficient \(k\) of
\(P_j\) by \(q^{(nN/j)k(j-k)}\). In degree \(j\), form the monic
polynomial of the ordered off-diagonal root quotients, adjoining
one factor \(u-1\) only when \(j\) is odd. Write \(R_n(q)\) for
the resultant of these two monic quotient polynomials.
[Common-clock quotient resultant](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)

The bound \(\deg R_n<D=J^6\) holds uniformly in the normalized
coefficients, and
\[
                    |R_n(0)|=j_{\rm ev}^{\,j_{\rm ev}}\ge1 ,
\]
where \(j_{\rm ev}\) is the even member of \(n,N\).
The degree bound follows from the Sylvester formulas for the
quotient polynomials and their resultant. At \(q=0\) the two
polynomials are \(1+z^n,1+z^N\). The nontrivial quotient sets
are disjoint by \(\gcd(n,N)=1\), and the dummy root \(1\)
from the odd order accounts for the displayed value.
[Degree bound and normalization](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-entire-one-sided-neighborhood-of-zero-added-curvature)

Put \(M=1+(J-3)D\) and use the menu
\[
                   \tau_i=a+\frac{i(b-a)}{M+1},
                   \qquad 1\le i\le M.
\]
For every low pair the consecutive menu clocks
\(q_{n,i}=e^{-\tau_i/(2nN)}\) are separated by at least
\(e^{-b/12}(b-a)/(2J^2(M+1))\), by the mean-value theorem.
Take a positive rational \(\epsilon\le1/2\) with \(4\epsilon\)
strictly below that bound.

To see the uniform effect of the menu, factor any \(R_n\) as
\(R_n(0)\prod_\ell(1-q/r_\ell)\) over its nonzero complex roots,
counted with multiplicity. If \(0<q\le1\) and
\(\lvert R_n(q)\rvert<\epsilon^D\), some factor has modulus
less than \(\epsilon\). For that factor
\[
           |r_\ell|<2,\qquad |q-r_\ell|<2\epsilon .
\]
Its disk contains at most one of the menu clocks for this
pair. Thus fewer than \(D\) menu strengths can fail the
bound \(|R_n(q_{n,i})|\ge\epsilon^D\) for this pair. There
are \(J-3\) low pairs, so at least one menu strength satisfies
the bound simultaneously for all of them. This uses only the
degrees and the fixed nonzero constant values, regardless of
the arithmetic of the normalized coefficients.

For a finite ordinary core the choice is effective.
Compute the quotient resultants and evaluate them with
arbitrarily sharp rational enclosures of the algebraic
coefficients and exponential clocks. Select any menu
point certified to have \(|R_n(q_{n,i})|>\epsilon^D/2\)
for all the low pairs. At least one point has
\(|R_n(q_{n,i})|\ge\epsilon^D\), so the search terminates.

The same disks give a uniform measure estimate on \([a,b]\).
For the pair \(n,N\) the exceptional set
\(|R_n(e^{-\tau/(2nN)})|<\epsilon^D\)
has length at most \(8nN D e^{b/(2nN)}\epsilon\).
Indeed each disk cuts out at most \(4\epsilon\) in the
real clock variable, and \(|d\tau/dq|\le 2nN e^{b/(2nN)}\)
on the interval. Over all low pairs the bound is at most
\[
                   8J^2(J-3)D e^{b/12}\epsilon .
\]
The estimates also apply to any fixed countably retained core,
taking the limits of the normalized coefficient arrays.
They control its low-order Gaussian references. The use of
the absolute cutoff below is for actual ordinary extensions.

## Uniform transfer to the appended rows

The endpoint-normalized references have roots of modulus
between \(2^{-j}\) and \(2^j\), so all root quotient moduli
are at most \(2^{2N}\). Let \(m_j=j(j-1)+\mathbf1_{j\ {\rm odd}}\).
For roots \(z,w\) of the two references the quotient resultant
product contains \(z/\bar z-w/\bar w\), with the odd dummy
quotient available at real roots. Even reference polynomials
are strictly positive on the real axis at positive strength
by lognormal averaging of the ordinary nonnegative even
moments. Factoring the difference of the squares of the unit
projections shows that a certified menu point has reference
chord gap at least
\[
            g=\min_{2\le n<N<J}
               \frac{\epsilon^D}
                    {4(2^{2N+1})^{m_nm_N-1}}>0 .               \tag{1}
\]

For appended shared rows with \(\eta\le1/2\), the log
coefficient difference between the ordinary extension and
the reference at strength \(\tau_*\) is at most
\[
                          j\tau\eta/4+j|\tau-\tau_*|/8.
\]
The first term is the
[matched diffuse Gamma chord bound](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#uniform-transfer-to-the-added-ordinary-rows).
The second comes from replacing the Gaussian strength
\(\tau\) by \(\tau_*\). Both bounds factor exactly from the
core contributions. If the resulting coefficient log
difference is at most \(e_0\le1\), the chord distance from
an ordinary extended root to a reference root direction is
at most \(A_J e_0^{1/(J-1)}\), where
\[
                  A_J=\max_{2\le j<J}2^{j+1}(1+2^j).
\]
Indeed at an ordinary root the monic reference polynomial
has modulus at most \(e_0(1+2^j)^j\). Factoring it gives a
reference root at distance \((1+2^j)e_0^{1/j}\), and the
unit projection costs at most \(2^{j+1}\) by the root
annuli.

Choose positive rational \(\delta_*,\eta_0\) with
\(\delta_*<(b-a)/(2(M+1))\), \(\eta_0\le1/2\), and
\[
 e_0=J(b\eta_0/4+\delta_*/8)\le1,\qquad
                       2A_J e_0^{1/(J-1)}<g.
\]
Use \(\eta_*=\eta_0\). Then the two ordinary angle sets still
separate at every low pair. Terminal two separates by
the even real positivity and the real linear root. The
absolute ordinary cutoff covers all \(N\ge J\), independently
of the core and the added rows. There is no limitation on
the physical signed degree drift because the entire transfer
compares arguments after positive rescalings.


## Exact ordinary placement of a selected strength

Let \(\tau_*>0\) be a selected rational strength, let
\(h\in\mathbb Z\) be the desired signed degree contribution
from the appended rows, and let \(\lambda\in[-\tau_*,\tau_*]\).
The last number is optional and prescribes a limiting signed
log ratio of the new exponent endpoints.

Choose a large positive integer \(m\) with
\(M_0=m\tau_*\) an integer. Start with \(M_0\) unit-gap rows
of smaller exponent \(m-1\). Orient \(M_+\) of them
increasing and \(M_-=M_0-M_+\) decreasing, where
\[
              M_+ = m(\tau_*+\lambda)/2+O(1),
                    \qquad 1\le M_+\le M_0-1 .
\]
Such a choice exists for all large \(m\), including the
endpoint values of \(\lambda\). A row of gap one whose
larger exponent is \(x\) contributes curvature \(1/x\).
Replacing it by two rows with larger exponent \(2x\)
in its own orientation preserves both curvature and
signed curvature, but increases the number of rows of that
orientation by one. Starting with both orientations,
repeat this splitting in the required orientation until
the signed degree contribution is exactly \(h\).
It takes finitely many steps for every \(m\).

The tail curvature remains exactly \(\tau_*\) and the
largest relative gap is at most \(1/m\). Its signed
curvature is \(2M_+/m-\tau_*\), which tends to \(\lambda\).
For a unit row the unsigned log endpoint ratio is
\(-\log(1-1/x)\), and
\[
       0\le-\log(1-1/x)-1/x
                    \le\frac{1/x}{2(m-1)}\qquad(x\ge m>1).
\]
Summing over the split rows, the difference of the signed log
endpoint ratio from the signed curvature is at most
\(\tau_*/[2(m-1)]\). Thus the signed log ratio tends to
\(\lambda\) with the exact signed integer gap held at \(h\).
For the compactified physical center ratio, this gives the
usual signed correction \(-\lambda/(2nN)\) in addition to
the integer contribution \(h\log(N/n)\) at each fixed pair.
[Center drift in the exact compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)

When \(h=0\) and no nonzero signed log ratio is requested,
choose \(m\) with \(m\tau_*\) even and simply use mirrored
pairs at \(m-1\), with no splitting. Each pair multiplies
the two raw endpoint coefficients by the same positive
factor in every order. This preserves their ratio exactly.

The selected extensions can have arbitrarily small positive
curvature by choosing the menu interval near zero. This
does not imply that a core is regular: the uniform gap (1)
depends on the interval and can vanish as the interval
shrinks to zero. The theorem selects a perturbation with
uniform stability, rather than bounding the unperturbed
factorial pair.

