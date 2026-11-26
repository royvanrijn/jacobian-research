# Finite cores detect nondegenerate interior compactified contacts

Fix adjacent orders \(n,N=n+1\), \(n\ge2\), and an integer
total degree difference \(\Delta_0\). Use the
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data).
The retained core is a finite or countable multiset \(R\) of
ordinary unequal rows with summed curvature
\(\sum_R (b-a)^2/\max(a,b)<\infty\). Its endpoint-normalized
reference polynomials, with Gaussian strength \(T>0\), are
\[
 F_j(T,z)=\sum_{k=0}^j a_{j,k}\exp[-jT\phi(k/j)]z^k,\qquad
 \phi(s)=\tfrac12s(1-s),\qquad j=n,N .
\]
Here \(a_{j,k}\) includes the binomial coefficient and every
retained row is kept exact. The relative physical argument
scale is \(\exp(\Gamma-\lambda/(2nN))\), where
\(\Gamma=\Delta_0\log(N/n)+\sum_R E(w)\) uses the center
correction \(E(w)\) in the compactification.

Suppose \(F_n(T,z)\) and \(F_N(T,e^{\Gamma-\lambda/(2nN)}z)\)
have a physical common root at \((T_*,\lambda_*)\) with
\(|\lambda_*|<T_*\). Assume both contacting roots are simple.
Write \(z_j(T)\) for their branches in their respective
endpoint-normalized variables, and choose their continuous
arguments so that the difference
\[
                 \theta(T)=\arg z_N(T)-\arg z_n(T)
\]
vanishes at \(T_*\). For real polynomials \(H\) on \([0,1]\)
with zero endpoint values, define the phase response
\[
 \mathcal L_T(H)=
  \left.\frac{d}{du}
       \big(\arg z_{N,u}(T)-\arg z_{n,u}(T)\big)\right|_{u=0},
 \quad
 F_{j,u}(T,z)=\sum_{k=0}^j [z^k]F_j(T,z)\,e^{u j H(k/j)}z^k .
\]
Let
\[
 \Psi(s)=\frac{s(1-s)(2s-1)}{12},\qquad
 \Xi(s)=\frac{(s-1/2)^4-1/16}{12}.
\]
These are the signed cubic and the even quartic log-coefficient
directions, respectively.

**Theorem.** In each of the following cases there is a sequence
of finite ordinary cores and simple sign-changing interior
Gaussian contacts at the same pair of orders and the same total
integer \(\Delta_0\), tending to the specified compactified
contact in normalized coefficients, center and parameters.

1. The original \(\theta\) changes sign at \(T_*\).
2. \(\mathcal L_{T_*}(\Psi)\ne0\). The new contacts may be chosen
   transverse.
3. The retained multiset is infinite,
   \(\theta'(T_*)=0\), \(\mathcal L_{T_*}(\Psi)=0\), and the
   common angular velocity \((\arg z_N)'(T_*)=(\arg z_n)'(T_*)\)
   is nonzero.
4. The retained multiset has unbounded absolute gaps,
   both angular velocities and \(\mathcal L_{T_*}(\Psi)\) vanish,
   and \(\mathcal L_{T_*}(\Xi)\ne0\).

If the original total curvature is strictly below a specified
cap, the same cap holds eventually for these references.
All the finite cores use ordinary integer rows. The Gaussian
strength and center correction remain real diffuse parameters.

## Contact order and finite truncation

The series of retained log-coefficients and center corrections
converge absolutely by the estimates in the compactification.
Their finite truncations converge, with every \(T\) derivative of
the reference polynomials on a bounded neighborhood of \(T_*\).
For each simple root choose a small disk with no other root of
its polynomial at \(T_*\), and restrict \(T\) to a neighborhood
where the disk stays isolating. The roots and their logarithms
are analytic in \(T\) and in small log-coefficient
perturbations. This follows from the implicit function theorem
with \(z_j F_j'\ne0\). The chosen branches and their first
\(T\) derivatives therefore converge under finite truncation.

The contacting roots are nonreal. The even core polynomial
is nonnegative on the real axis as a limit of
endpoint-normalized ordinary even moment polynomials. Gaussian
attenuation at positive strength is an average over positive
lognormal rescalings, so the even reference is strictly positive
on the real axis. Its constant coefficient is one, and for each
nonzero real argument the rescaling misses its finite real
zero set almost surely.

The real analytic function \(\theta\) has a zero of finite order.
In fact, if it vanished identically near \(T_*\), the two
nonreal roots and their conjugates would give equal
off-diagonal root quotients throughout the interval. The
[off-diagonal quotient resultant](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts)
for the two polynomials is a polynomial in
\(q=e^{-T/(2nN)}\). This construction uses only fixed nonzero
endpoint coefficients, so it applies to the retained real arrays.
At \(q=0\), the quotient sets of the two endpoint binomials
are disjoint as \(n,N\) are coprime. The resultant is therefore
nonzero. This rules out a locally identically zero difference.

If \(\theta\) changes sign, choose a successively shrinking
sequence of intervals on which its two endpoint signs are
opposite. They persist under sufficiently accurate finite
truncation. The perturbed angular difference is real analytic,
and hence has a sign-changing zero on each interval. This
establishes the angular part of case 1.

## Ordinary row perturbations at fixed effective strength

For a shared unequal row let \(m=(a+b)/2\), \(d=b-a\) be signed,
\(w=d^2/m\), \(A=d^2/m^2\), and \(\eta=|d|/m\).
Write the normalized log-coefficient contribution at index \(k\)
as \(j L_j(k/j)\), so \(L_j(0)=L_j(1)=0\).
Using \(g(x)=\log\Gamma(1+x)\), the second derivative of the
exact endpoint chord is
\[
 L_j''(s)=j d^2 g''\big(j(m+d(s-1/2))\big).
\]
For fixed orders and \(\eta\le1/2\), writing \(\xi=s-1/2\) gives
\[
 L_j''(s)=w-Ad\,\xi-\frac{A}{2j}+Aw\,\xi^2
               +O_j\big(A(w\eta+\eta+m^{-1})\big).             \tag{1}
\]
To check the error, expand \(g''(x)=x^{-1}-\tfrac12x^{-2}
+O(x^{-3})\) at the indicated argument.
The first term contributes the first, second and fourth terms
on the right of (1) with a residual \(O(Aw\eta)\).
The second term contributes \(-A/(2j)\) with residual
\(O_j(A\eta)\), and the trigamma remainder contributes
\(O_j(A/m)\). The trigamma estimate follows directly from
\(\psi_1(1+x)=\int_0^\infty e^{-xu}u/(e^u-1)\,du\), using
\(u/(e^u-1)=1-u/2+O(u^2)\) for \(u\ge0\). The integral
identity results from integrating the nonnegative terms in
\(\sum_{h\ge1}(x+h)^{-2}\).

On adding the row, decrease the Gaussian strength from \(T\)
to \(T-w\), and on removing it increase the strength to \(T+w\).
The retained plus Gaussian log-coefficients then change by
\(\pm j R_j(k/j)\), respectively, where
\[
 R_j(s)= L_j(s)+w\phi(s)
       = Ad\,\Psi(s)+\frac{A}{2j}\phi(s)+Aw\,\Xi(s)
                  +O_j\big(A(w\eta+\eta+m^{-1})\big).           \tag{2}
\]
Indeed \(\Psi''=-\xi\), \(\phi''=-1\), \(\Xi''=\xi^2\).
All terms have zero endpoints. A bounded second derivative
error with zero endpoint values gives a bounded function error
on this interval, by quadratic comparison.

Each fixed small log-coefficient perturbation has a linear
angular response with a quadratic error, uniformly near \(T_*\)
and also after one \(T\) derivative. This follows from
analyticity of the two root logarithms on one common small
neighborhood. In particular the first-order phase response of
the clock term \(\frac{A}{2j}\phi\) in (2) is
\[
             \frac{A}{2}
               \left(\frac{(\arg z_n)'(T)}{n}
                          -\frac{(\arg z_N)'(T)}{N}\right),    \tag{3}
\]
because multiplication by \(\exp(-jT\phi)\) defines the reference
clock. The other two displayed responses are
\(Ad\,\mathcal L_T(\Psi)\) and \(Aw\,\mathcal L_T(\Xi)\).
The quadratic error is \(O((A|d|)^2)\) along any weak row sequence
with \(w\to0\).

## Signed unfolding

For case 2, add either \((t^3-t,t^3+t)\) or its reversal with
positive integer \(t\to\infty\). The variance compensation is
\(w=4/t\). Put \(\varepsilon=8/t^3\). In (2) the leading
direction is \(\pm\varepsilon\Psi\) and all other terms are
\(o(\varepsilon)\), uniformly in the log-coefficient arrays.
First truncate \(R\) accurately enough that its coefficient
errors are \(o(\varepsilon)\). Let \(p\ge1\) be the finite
order of the original angular zero and
\(\theta(T_*+s)=c s^p+O(s^{p+1})\), \(c\ne0\).
The resulting angular difference, evaluated at
\(T=T_*+\varepsilon^{1/p}y\) and divided by \(\varepsilon\),
converges on bounded \(y\) intervals, including its first
\(y\) derivative, to
\[
                         c y^p\ \pm\ \mathcal L_{T_*}(\Psi).
\]
For even \(p\) choose the orientation giving opposite signs
of the nonzero constant and \(c\). For odd \(p\) either
orientation works. The limit polynomial then has a nonzero
simple real root. The convergence and the intermediate value
theorem give a nearby zero with nonzero first derivative.
The physical strength of this finite reference is \(T-w\),
which tends to \(T_*\).

For case 3, pick a sequence of retained rows with \(w\to0\)
and \(m\to\infty\), available from summability and integrality.
One may remove one copy or add one copy of each row. Use
the positive or negative perturbation in (2) accordingly.
At \(T_*\) the cubic response vanishes, the quartic term is
\(o(A)\), and
\[
 (A|d|)^2/A= w^2\longrightarrow0.
\]
Thus (3) gives
\[
  \theta_{\rm pert}(T_*)=
                 \pm A\, \frac{v}{2nN}+o(A),
  \qquad v=(\arg z_n)'(T_*)=(\arg z_N)'(T_*)\ne0 .
\]
The normalized coefficient errors from finite truncation
can be made \(o(A)\) while retaining the selected row before
modifying its count. If the original angular zero has even
order, choose the new sign at \(T_*\) opposite to the sign
at sufficiently close fixed interval endpoints. Those signs
at the endpoints persist as the perturbations tend to zero.
Real analyticity gives a sign-changing zero inside. By taking
smaller intervals and accurate truncations this gives a
sequence tending to the original contact. Odd-order zeros
were already handled in case 1.

For case 4, take retained rows with \(|d|\to\infty\) and
\(w\to0\). The two leading responses vanish at \(T_*\).
The remainder in (2) is \(o(Aw)\), using
\(\eta/w=1/|d|\) and \(1/(mw)=1/d^2\).
The quadratic error divided by \(Aw\) is \(O(w)\).
The first nonzero phase change is therefore
\[
              \pm Aw\,\mathcal L_{T_*}(\Xi)+o(Aw).
\]
Choose the sign and the finite truncations as in case 3,
using accuracy \(o(Aw)\). This again gives sign-changing
zeros tending to \(T_*\).

## Recovering the physical scale

For each angular zero of the finite reference, set
\[
     \lambda_{\rm fin}
      =2nN\bigl(\Gamma_{\rm fin}
                    -\log|z_{N,\rm fin}/z_{n,\rm fin}|\bigr),
 \qquad
     \Gamma_{\rm fin}
        =\Delta_0\log(N/n)+\sum_{R_{\rm fin}} E(w).
\]
Then the two roots agree after the required positive
physical rescaling. A row added or removed above makes a
center contribution \(O_{n,N}(w)\), by the compactification
estimate. Hence the centers converge and
\(\lambda_{\rm fin}\to\lambda_*\). The strict interior
condition persists. When the retained core is finite,
the raw-core twist integer is \(\Delta_0-\sum_{R_{\rm fin}}d\).
The individual root simplicity persists throughout the
construction. The strength, retained curvature and physical
pair tend to the original ones, which also proves the
assertion about a strict curvature cap.
