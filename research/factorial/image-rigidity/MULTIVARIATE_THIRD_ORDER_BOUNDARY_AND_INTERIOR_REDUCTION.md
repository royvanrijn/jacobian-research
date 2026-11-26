# Third-order boundary tests and reduction to interior Gaussian contacts

Fix adjacent orders \(n,N=n+1\), where \(n\ge2\), and the integer
degree difference of the full profile. Consider a reference in the
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
with Gaussian strength \(T_*>0\) and signed bias
\(\lambda_*=\epsilon T_*\), \(\epsilon\in\{1,-1\}\). The retained core
may be finite or countable with summable curvature. Let \(z_j(T)\)
be the endpoint-normalized simple nonreal contacting roots, continued
for \(T\) near \(T_*\) before applying the positive physical scale.
Write
\[
 C_R=\Delta_0\log(N/n)+\sum_{\rm retained}E(w),\qquad
 H(T)=\log(z_N(T)/z_n(T))-C_R+\epsilon T/(2nN).
\]
Choose the logarithm real at the contact, so \(H(T_*)=0\). Here
\(E(w)\) is the exact signed adjacent-center correction of the
individual row. The physical reference on the boundary has scale
\(\exp(C_R-\epsilon T/(2nN))\).

Assume first that the contact is transverse in angle, meaning
\(p=H'(T_*)\) satisfies \(\Im p\ne0\). Put
\[
 \Pi(x)=\epsilon\left(\Re x-\frac{\Re p}{\Im p}\Im x\right).
\]
The two coefficients \(A_2,A_3\) and the quadratic response \(D\)
are defined below using just the two exact reference polynomials
and their simple root branches.

**Theorem.** The following statements hold at the specified pair
of roots.

1. For a finite ordinary core, if \(A_2<0\), all sufficiently
   diffuse ordinary tails near the specified contact are excluded.
   The same holds if
   \[
               A_2=0,\qquad A_3<0,\qquad A_3+T_*D<0 .
   \]
   The latter is a third-order boundary test with an unrestricted
   tail row count and unrestricted tail heights and gaps.
2. If the retained multiset is infinite and infinitely many rows
   have signed gap different from \(\epsilon\), there are transverse
   interior Gaussian contacts of finite ordinary cores tending
   to the original boundary contact, with the same total degree
   difference. The same conclusion holds if infinitely many
   retained rows have gap \(\epsilon\) and
   \((A_2,A_3)\ne(0,0)\).
3. For any retained core, finite or countable, if the first nonzero
   entry in \((A_2,A_3)\) is positive, finite ordinary cores with
   transverse interior contacts tend to the boundary contact.
   The finite cores can contain all prescribed finitely many of
   the original retained rows.

Here convergence includes strength, bias, normalized coefficients,
and the relative center. Any strict cap on the initial total
curvature holds eventually for the constructed references.
Part 1 concerns actual integer factorial tails. Parts 2 and 3
produce reference contacts with Gaussian tails.

## Exact row expansions and the projected coefficients

For an unequal shared row of signed gap \(d=b-a\) and midpoint
\(m=(a+b)/2\), put \(u=1/m\), \(w=d^2u\), and \(\ell=du\).
A row is called aligned in this proof if \(d=\epsilon\).
At \(s=k/j\) its endpoint-normalized log coefficient is
\(jL_j(s)\), where \(L_j(0)=L_j(1)=0\). Write
\[
 \phi(s)=\tfrac12s(1-s),\quad
 \psi(s)=s(1-s)(2s-1)/12,\quad
 \xi(s)=\big((s-\tfrac12)^4-\tfrac1{16}\big)/12 .
\]
For aligned rows the precise third-order compensation is
\[
 \begin{split}
 j(L_j+\!w\phi)(k/j)&=u^2P_j(k)+u^3R_j(k)+O_j(u^4),\\
 P_j(k)&=j\epsilon\psi(k/j)+\tfrac12\phi(k/j),\\
 R_j(k)&=j\xi(k/j)-\epsilon\psi(k/j)-\frac{\phi(k/j)}{6j},\\
 E(w)&=-\epsilon u/(2nN)+\epsilon c u^2
                    -\epsilon u^3/(24nN)+O_{n,N}(u^4),
 \quad c=(n^{-2}-N^{-2})/12 .
 \end{split}                                                    \tag{1}
\]

To check the chord formula set \(x=m+d(s-\frac12)\). Its exact
second derivative is \(d^2j\psi_1(1+jx)\), where \(\psi_1\)
denotes the trigamma function. Summing the Laplace transforms
of the terms in the [trigamma series](https://dlmf.nist.gov/5.15.E1)
gives
\[
 j\psi_1(1+jx)=\int_0^\infty e^{-xv}
                 \frac{v/j}{e^{v/j}-1}\,dv
 =\frac1x-\frac1{2jx^2}+\frac1{6j^2x^3}+O_j(x^{-5}).
\]
The error bound follows from the Taylor polynomial of
\(y/(e^y-1)\) through degree three, whose cubic coefficient
is zero. Its error is bounded by a constant times \(y^4\)
on \([0,1]\) by analyticity and on \([1,\infty)\) by
\(0<y/(e^y-1)\le1\).
Expanding at \(m\), the terms through order \(m^{-3}\) of the
second derivative for fixed \(d\) are
\[
 d^2u+d^2u^2[-d(s-\tfrac12)-1/(2j)]
       +d^2u^3[d^2(s-\tfrac12)^2+(d/j)(s-\tfrac12)+1/(6j^2)].
\]
The three zero-endpoint primitives have second derivatives
\(\phi''=-1,\ \psi''=-(s-\frac12),\ \xi''=(s-\frac12)^2\).
This gives the chord formula for \(d=\epsilon\).
For the center, use the exact integral
\[
 E(w)=\int_n^N\int_a^b
                 \bigl(t\psi_1(1+xt)-1/x\bigr)\,dt\,dx .
\]
The expansion of the integrand has terms
\(-1/(2x^2t)+1/(6x^3t^2)\) and remainder \(O_{n,N}(t^{-4})\).
The signed integrals of \(1/t\) and \(1/t^2\) across a fixed-gap
row are \(d u+d^3u^3/12+O_d(u^5)\) and
\(d u^2+O_d(u^4)\), respectively. This proves (1).

For a pair of multiplier arrays \(L=(L_j(k))\), write
\(\mathcal R_T(L)\) for the derivative at zero of the logarithm
of the root ratio after multiplying the two reference coefficient
arrays by \(\exp(v L_j(k))\).
Use \(K(T)\) for half the second derivative of the log root ratio
for the particular direction \(L=P\). All these are analytic in
\(T\) near \(T_*\), since the two roots are simple and nonzero.
Define
\[
 \begin{split}
 J_2(T)&=\mathcal R_T(P)-\epsilon c,\qquad
 J_3(T)=\mathcal R_T(R)+\epsilon/(24nN),\\
 A_i&=\Pi(J_i(T_*)),\qquad
 t_1=-\Im J_2(T_*)/\Im p,\\
 D&=\Pi\left(K(T_*)+t_1 J_2'(T_*)
                              +\tfrac12t_1^2 H''(T_*)\right).
 \end{split}                                                    \tag{2}
\]
The coefficient \(A_2\) is the leading boundary sign in the
[second-order test](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#transverse-contacts-on-the-signed-bias-boundary).

## Exact ordinary exclusion on the aligned cone

For the actual added tail define
\[
 W=\sum_i w_i,\qquad \lambda=\sum_i \ell_i,\qquad
 \delta=W-\epsilon\lambda,\qquad
 \eta=\max_i|d_i|/m_i,\qquad
 Z=\sum_{\rm aligned}u_i^2,\quad U=\sum_{\rm aligned}u_i^3 .
\]
Every ordinary row has \(w_i-\epsilon\ell_i=w_i(1-\epsilon/d_i)\).
The rows outside the aligned class therefore have total curvature
at most \(2\delta\). Their compensated log coefficients are
\(O_j(\eta\delta)\), and the difference of their signed center
correction from \(-\sum\ell_i/(2nN)\) is \(O_{n,N}(\eta\delta)\).
Indeed, \(x=m_i+d_i(s-\frac12)\) stays comparable to \(m_i\),
and \(|j d_i^2\psi_1(1+jx)-w_i|\le C_j\eta w_i\) for
\(\eta\le1/2\), by the integral above and
\(m_i^{-1}\le|d_i|/m_i\). Integrating with zero endpoints gives
the coefficient bound. The center estimate follows from
\(\log(b_i/a_i)=\ell_i+O(\eta w_i)\) and from the integral
estimate replacing \(t\psi_1(1+xt)-1/x\) by \(-1/(2x^2t)\).
In particular no higher signed-gap assumption enters this bound.

The aligned rows contribute their two arrays \(ZP+UR\) up to
an error \(O_j(\eta U)\), and the corresponding sum of the
center expansions in (1). The inequalities used below are
\[
        Z\le\eta W,\qquad U\le\eta Z,\qquad
                                    Z^2\le WU .
\]
The last is Cauchy--Schwarz on the aligned heights.
In a sufficiently small neighborhood any ordinary contact
continuing the specified branch pair must satisfy
\[
 0=H(W)-\frac{\epsilon\delta}{2nN}
         +ZJ_2(W)+UJ_3(W)+Z^2 K(W)+O_{C,n,N}(\eta(U+\delta)).
                                                               \tag{3}
\]
For this equation multiply each retained coefficient array by the
exact endpoint-normalized tail factor. Analytic root continuation
gives the first-order effects \(Z\mathcal R_W(P)+U\mathcal R_W(R)\)
and the second-order effect \(Z^2 K(W)\). All other terms belong
to the stated error, since \(ZU=O(\eta U)\) and
\(Z^3=O(\eta U)\) by \(Z^2\le WU\). The center scale of the full
ordinary pair is \(C_R+\sum_{\rm tail}E(w_i)\) at the fixed
total integer \(\Delta_0\). Inserting
\(\lambda=\epsilon(W-\delta)\) gives (3).

Write \(x=W-T_*\). The imaginary equation first gives
\(x=O(Z+\eta\delta)\), by shrinking the neighborhood so that
the Taylor remainder in \(H\) can be absorbed using \(\Im p\ne0\).
Substitution in the same equation gives
\(x=t_1 Z+O(U+\eta\delta)\) using \(Z^2\le WU\).
Apply \(\Pi\) to (3). It annihilates \(p\), so the linear
term \(px\) disappears. Taylor expansion and these two
bounds give
\[
 \frac{\delta}{2nN}
        =A_2 Z+A_3 U+D Z^2+O_{C,n,N}(\eta(U+\delta)).
                                                               \tag{4}
\]
For the quadratic term of \(H\), replacing \(x^2\) by \(t_1^2Z^2\)
costs \(O(\eta(U+\delta))\). For \(ZJ_2(W)\) the analogous
replacement in the linear Taylor term costs the same. Cubic
remainders are absorbed by \(x=O(Z+\eta\delta)\), bounded \(W\)
and \(Z^2\le WU\). This accounts for all terms in \(D\).

When \(A_2<0\), the remaining \(U,Z^2\) terms have size
\(O(\eta Z)\); equation (4) excludes \(\delta\ge0\) for small
\(\eta\). When \(A_2=0\) and both test signs in part 1 are
strictly negative, the affine expression
\(A_3+D Z^2/U\) is strictly negative uniformly for
\(0<U\) and \(0\le Z^2/U\le W\) with \(W\) close to \(T_*\).
It follows from (4) that \(\delta/(2nN)\le -c_0 U+
C_0\eta(U+\delta)\) for some \(c_0>0\). With \(\eta\) sufficiently
small this contradicts \(\delta\ge0\). If \(U=0\), all rows
are outside the aligned class, so \(\delta\ge W/2>0\), also
impossible in (4) when \(\eta\) is sufficiently small.

For a finite ordinary core all the arguments give effective
neighborhoods. The reference clock \(q_*=\exp(-T_* /(2nN))\)
and root coordinates are algebraic after accounting for the
positive algebraic endpoint scales. Differentiation in \(T\)
at a simple root can be carried out by algebraic operations
on the coefficient and root values, since each coefficient
is a fixed constant times a known exponential.
Consequently \(p,A_2,A_3,D\) are algebraic. The strict test
sign involving \(T_*\) is decidable by rational intervals:
if \(D\ne0\), equality is impossible by the transcendence of
\(T_*\) at positive algebraic \(q_*\).
Taylor and root-continuation bounds can be obtained over
isolating root disks and rational coefficient neighborhoods
as in the second-order test.

## Inward perturbations of retained cores

For a selected weak retained row, allow either adding one
copy or deleting one existing copy and denote the signed
change by \(\sigma=1,-1\), respectively. At the same time
replace the Gaussian strength \(T\) by \(T-\sigma w\).
This keeps the leading coefficient curvature fixed. The
perturbed boundary expression is
\[
 H_\sigma(T)=\log(z_{N,\sigma}/z_{n,\sigma})
       - C_R-\sigma E(w)+\epsilon(T-\sigma w)/(2nN).
\]
At an angular zero of \(\Im H_\sigma\), the unique bias
giving physical contact has distance from the signed boundary,
measured in the inward direction, equal to
\(2nN\epsilon\Re H_\sigma(T)\). The opposite boundary remains
a positive distance away near this initial positive contact.

For any sequence of weak retained rows of nonaligned gap,
summability gives \(w\to0\), \(m\to\infty\), and \(\eta\to0\).
The compensated coefficient correction is \(O_j(\eta w)=o(w)\)
and the center correction is \(-\ell/(2nN)+o(w)\) by the
same row estimates used in (3). Delete that row, choosing
\(\sigma=-1\). Then
\[
      H_{-}(T)=H(T)+\frac{\epsilon}{2nN}(w-\epsilon\ell)+o(w)
\]
uniformly near \(T_*\), with uniform convergence of the \(T\)
derivatives of the logarithmic root contributions.
Since the displayed leading change is real, the nearby
angular zero moves by \(o(w)\), by transversality.
There its inward slack is \(w-\epsilon\ell+o(w)\ge w/2+o(w)>0\).
For each selected weak row, truncate the rest of the countable
retained multiset with error \(o(w)\), retaining that row
before deletion. Absolute convergence of the log coefficients
and signed center corrections permits such truncations. The
resulting finite cores give interior contacts converging to the
original one.

For an aligned unit row, (1) gives instead
\[
        H_\sigma(T)=H(T)+
               \sigma u^2 J_2(T)+\sigma u^3 J_3(T)+O(u^4).
\]
Uniform expansion holds after one \(T\) derivative by analytic
root continuation. The imaginary equation has a unique nearby
zero with displacement \(O(u^2)\). Projecting by \(\Pi\) at
this zero, its inward slack is
\[
                 \sigma\,2nN(A_2u^2+A_3u^3)+O(u^4).
\]
Indeed the linear term in the displacement is annihilated
and every term involving its square or its product with \(u^2\)
is \(O(u^4)\). Thus the sign of the first nonzero projected
coefficient determines which of addition or deletion produces
strict inward slack. For addition, ordinary aligned unit
rows at arbitrarily large half-integral midpoint are available.
For deletion, take a sequence of retained aligned unit rows
with midpoints tending to infinity. Truncate the other retained
rows with errors \(o(u^3)\) while preserving each designated
row before it is deleted. The simple angular zeros remain
transverse under these perturbations.

Throughout the construction, the total integer \(\Delta_0\)
is kept in \(C_R\). For a finite retained multiset the integer
gap assigned to the Gaussian tail is \(\Delta_0\) minus the
sum of the retained signed gaps. The exact row center corrections
are thus changed only by \(\sigma E(w)\), as required. The
Gaussian strength stays positive and the strict curvature cap
persists: the compensation changes total curvature by a
quantity tending to zero, since \(d^2/\max(a,b)\) and \(w\)
differ by \(O(\eta w)\) for the selected shared row.


## Stationary boundary contacts with a nonzero radial derivative

Keep the boundary setup but now suppose the angular order is an
integer \(P\ge2\), so
\(\Im H(T_*+x)=c_\theta x^P+O(x^{P+1})\), \(c_\theta\ne0\).
Assume \(r=\epsilon\Re H'(T_*)\ne0\).
Each of the following hypotheses supplies finite ordinary cores
with simple transverse interior contacts tending to the original
boundary contact at the same total integer degree difference.

- For an allowable sign \(\sigma\), the polynomial
  \(c_\theta y^P+\sigma\Im J_2(T_*)\) has a nonzero simple real
  root \(y\) with \(r y>0\). Here addition makes \(\sigma=1\)
  allowable; deletion makes \(\sigma=-1\) allowable if there
  are infinitely many retained aligned unit rows.
- \(P=2\), \(\Im J_2(T_*)=0\), and for an allowable \(\sigma\)
  the numbers \(c_\theta\) and \(\sigma\Im J_3(T_*)\) have opposite
  signs.
- The imaginary part of
  \(\mathcal R_{T_*+x}((j\psi(k/j))_{j,k})\) has a first term
  \(b x^S\) with \(b\ne0\), and \(P\ge2S+3\).
  This case uses addition of one growing-gap ordinary row and
  permits either orientation.

For the first case, add or remove an aligned unit row of inverse
midpoint \(u\to0\) and use the compensated strength from above.
At \(T=T_*+u^{2/P}y\), the imaginary part of \(H_\sigma\)
divided by \(u^2\) converges to the stated polynomial, with
its first \(y\) derivative. Finite truncations can have error
\(o(u^2)\). A nearby angular zero is transverse and has inward
slack
\[
             2nN r u^{2/P}(y+o(1)).
\]
The real boundary residual from the row is \(O(u^2)=o(u^{2/P})\)
because alignment cancels the first-order variance and bias
contributions to that residual. The second case uses
\(T=T_*+u^{3/2}y\) and divides the angular equation by \(u^3\).
The term \(u^2(\Im J_2(T)-\Im J_2(T_*))=O(u^{7/2})\) is
negligible and truncation error is chosen \(o(u^3)\).
The limiting polynomial \(c_\theta y^2+\sigma\Im J_3(T_*)\)
has two opposite nonzero simple roots, so one gives inward
slack as above.

For the third case choose integers \(a,b_0>0\) with
\[
             \frac{2S}{P-2S}< b_0/a <P-S-2 .
\]
There is room since \(P-2S\ge3\) and \(P-S-2\ge S+1\).
For integers \(t\to\infty\) use the ordinary row of midpoint
\(m=t^{2b_0+a}\) and signed gap \(d=\pm2t^{b_0}\).
Its variance is \(w=4t^{-a}\). Put \(v=8t^{-b_0-2a}\) and
\(\nu=P/(P-S)<2\).
The [compensated chord expansion](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md#ordinary-row-perturbations-at-fixed-effective-strength)
gives \(j(L_j+w\phi)=\pm v j\psi+o(v^\nu)\) on the two
sampling grids. Indeed the terms after the first cubic term
are \(O(t^{-2b_0-2a})\), and
\(2b_0+2a>(b_0+2a)\nu\) by the lower bound on \(b_0/a\).
The quadratic root-response error \(O(v^2)\) is also \(o(v^\nu)\).
Choose finite truncation errors \(o(v^\nu)\). At
\(T=T_*+v^{1/(P-S)}y\) the normalized angular equation tends
with one derivative to \(c_\theta y^P\pm b y^S\).
At a nonzero root \(y\), its last factor is
\(c_\theta y^{P-S}\pm b\). If \(P-S\) is even, one orientation
gives two opposite simple real roots. If \(P-S\) is odd, each
orientation gives one nonzero real simple root and changing
orientation reverses its sign. In either case choose a simple
root with \(r y>0\).

For these growing-gap rows the leading real boundary residual
costs \(O(w)\). The upper bound on \(b_0/a\) gives
\(w=o(v^{1/(P-S)})\). The rest of the row and center changes in
\(H_\sigma\) is \(O(w)\) as well, since \(|\ell|\le w\) and
the compensated angular and radial root displacements are
\(O(v)\). Hence the inward slack at the perturbed angular zero
is \(2nN r\,v^{1/(P-S)}(y+o(1))>0\).
All cases keep the Gaussian strength positive, and the
inward slack establishes the admissible physical bias just
as for the transverse boundary reduction. Convergence of
parameters, normalized coefficients, and center follows
from the vanishing row curvature and the arbitrarily accurate
finite truncations.

