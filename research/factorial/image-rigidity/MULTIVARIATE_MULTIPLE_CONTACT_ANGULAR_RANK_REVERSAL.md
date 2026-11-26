# Angular rank reversal at multiple interior Gaussian contacts

Fix adjacent orders \(n,N=n+1\), with \(n\ge2\), and an ordinary finite
or countably infinite retained core \(R\) of summable curvature.
Use the [fixed-terminal reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at a physical interior contact \((T_*,\lambda_*)\), where \(T_*>0\)
and \(|\lambda_*|<T_*\). The integer total degree difference is
\(\Delta_0\). Write its endpoint-normalized polynomials as
\[
 F_j(T,z)=\sum_{k=0}^{j}a_{j,k}
              e^{-T k(j-k)/(2j)}z^k,\qquad j=n,N .
\]
Let \(\xi_n,\xi_N\) be logarithms of the two contacting roots in
their own endpoint scales, with
\(\operatorname{Im}\xi_n=\operatorname{Im}\xi_N=\vartheta_*\).
These roots are nonreal. At each order define the local coefficients
\(g_j\ne0,A_j,B_j\) and multiplicity \(m_j\) by
\[
 e^{-j(\xi_j+x)/2}F_j(T_*,e^{\xi_j+x})
       =g_j x^{m_j}\bigl(1+A_jx+B_jx^2+O(x^3)\bigr).
\]
Set \(V_j=\operatorname{Im}(A_j)/j\).

**Theorem.** Suppose at least one multiplicity is at least two.
If \((V_n,V_N)\ne(0,0)\), the given contact is a limit of simple
transverse interior Gaussian contacts of finite ordinary cores.
There is no upper bound on the multiplicities in this test, and
the other root may be simple. The consecutive orders and
\(\Delta_0\) are preserved, the new strengths and biases converge
to \(T_*,\lambda_*\), and any strict total curvature bound persists.

There is a second test when \(V_n=V_N\). Write \(h_{j,r}\),
\(1\le r\le m_j\), in increasing order for the real roots of
\(e^{-\partial_X^2/(2j)}X^{m_j}\). Let \((u_{j,r})\) be the
increasing rearrangement of the real list
\[
       \left(
          \frac{\operatorname{Im}(B_j-A_j^2/2)}{j}\,h_{j,r}
       \right)_{r=1}^{m_j}.
\]
If for some \(r,s\) one has
\[
             (u_{n,r}-u_{N,s})(h_{n,r}-h_{N,s})<0,             \tag{1}
\]
the same conclusion holds.

## Ordinary perturbations with simple transverse angular contacts

For a fixed \(n,N\), there is a proper algebraic subset of the space
of pairs of constant-normalized real coefficient arrays such that
outside it, throughout the real Gaussian clock, nonreal roots
are simple and every same-angle contact of nonreal roots is transverse.
The assertion concerns arbitrary such arrays with nonzero endpoint
coefficients. It also holds for finite ordinary cores obtained
by adjoining \(n+N\) increasing rows of gap two of independently
chosen, arbitrarily large integer smaller exponents.

To prove the algebraic assertion, use the nonzero clock variable
\(q=e^{-T/(2nN)}\). Each coefficient multiplier is then an integer
power of \(q\). The coefficients are complexified, with constants
fixed at one. Two distinct nonzero double roots at one order
impose four independent affine conditions on those coefficients at
fixed \(q\) and fixed roots. Independence follows by varying the
coefficients of \(z,z^2,z^3,z^4\). After dividing by \(z\), the
homogeneous kernel is a polynomial of degree at most three with
two double roots and is zero. The case of an order below four
has no such pair. Allowing \(q\) and two roots to vary gives at most
three extra dimensions. Projection and Zariski closure therefore
give a proper exceptional set. A nonreal double root of a real
polynomial implies two distinct double roots by conjugation.

For transversality, let \(P,Q\) denote the two Gaussian polynomials
at clock \(q\), and impose
\[
 P(\zeta)=P(\omega)=0,\qquad Q(r\zeta)=Q(r\omega)=0,
            \qquad r\zeta\omega(\zeta-\omega)\ne0.
\]
The four evaluation conditions are independent for fixed \(q,r,\zeta,\omega\)
since the coefficients of \(z,z^2\) in each polynomial can be solved
from them. The full incidence is irreducible of dimension \(n+N\).
Indeed the parameters \(q,r,\zeta,\omega\) range through a complex
Zariski open set and the unsolved coefficients are free. Work on
its open set with simple displayed roots and nonzero leading
coefficients. At a simple root \(z\) of the order \(j\) polynomial,
the heat identity
\(\partial_TF_j=(D^2-jD)F_j/(2j)\), \(D=z\partial_z\), gives
\[
        \partial_T\log z
            =-\frac{1-j+zF_j''(z)/F_j'(z)}{2j}.
\]
Consequently, at a real same-angle contact with \(\omega=\bar\zeta\)
and \(r>0\), equal angular velocities require
\[
 \frac1n\left(\frac{\zeta P''(\zeta)}{P'(\zeta)}
             -\frac{\omega P''(\omega)}{P'(\omega)}\right)
 =
 \frac1N\left(\frac{r\zeta Q''(r\zeta)}{Q'(r\zeta)}
             -\frac{r\omega Q''(r\omega)}{Q'(r\omega)}\right).
                                                               \tag{2}
\]
This is a nonidentical rational equation on the incidence.
To see it, fix the lower polynomial, the clock and the two upper
roots, and vary one of the remaining upper roots \(v\), which
exists since \(N\ge3\). The right side contains the term
\[
          \frac2N\left(\frac{r\zeta}{r\zeta-v}
                         -\frac{r\omega}{r\omega-v}\right)
\]
and its other terms are fixed. Normalize the constant coefficient
after choosing roots, which leaves both derivative ratios unchanged.
The varying term has distinct poles with nonzero residues. Thus (2)
cuts at least one dimension, and again the projection and Zariski
closure are proper. The two algebraic exceptions just constructed
can be united.

The [paired coefficient map](MULTIVARIATE_EXCEPTIONAL_CODIMENSION_TWO.md#differential-of-the-paired-coefficient-map)
from \(n+N\) independent bases of rows of gap two is dominant.
For completeness, in logarithmic consecutive coefficient ratios
the column from a base \(a\) has entries
\(1/(a+(2k-1)/j)+1/(a+2k/j)\). Each of these functions has a
noninteger pole unique to its own order and step, since \(n,N\)
are coprime. The \(n+N\) functions are independent, so one can
choose \(n+N\) evaluation points with independent columns.
Multiplying coefficient ratios by any fixed retained core factors
does not change dominance. The inverse image of the exceptional
algebraic set is therefore proper. Arbitrarily large independent
integer bases avoid it by Zariski density of the integer grid.
As each chosen base tends to infinity, the added endpoint-normalized
log coefficients, center corrections and curvatures tend to zero.
Positive endpoint rescaling is independent of \(T\) and changes
no angles or their velocities. Thus regular finite ordinary
approximations can preserve any finite list of strict angular
inequalities in isolating root disks. The same construction
works after finite truncation of a countable ordinary core.

In particular, suppose some common ranked angular difference in
isolating logarithmic root disks takes strictly opposite signs
at two nearby strengths. The disk boundaries may be kept free
of roots throughout the interval by choosing the interval small
enough about the original contact. After the regular perturbation,
each cluster has continuously varying roots, and its angles
in increasing order are continuous. A zero of the ranked
difference is now a contact of simple nonreal roots, and its
two analytic branches have distinct angular derivatives by (2).
This observation also applies to strict endpoint signs from a
simple sign-changing contact.

## The two sides of a heat collision

In a logarithmic coordinate the clock operator, up to a nonzero
scalar, is \(e^{(T-T_*)\partial_x^2/(2j)}\). This follows directly
by comparing the exponential frequencies \(k-j/2\), and is the
heat identity used in the [unitary Hermite calculation, equation (1.4)](https://arxiv.org/pdf/2203.05533).
For \(\sigma=\pm1\) and \(T=T_*+\sigma\epsilon^2\), put
\(P_{\ell}=e^{a\partial_X^2}X^\ell\) with \(a=\sigma/(2j)\).
After writing the logarithmic displacement as \(\epsilon X\) and
dividing by \(g_j\epsilon^{m_j}\) and the scalar, the resulting
local polynomial has the expansion
\[
          P_{m_j}(X)+\epsilon A_j P_{m_j+1}(X)
                       +\epsilon^2 B_j P_{m_j+2}(X)
                       +O(\epsilon^3).                       \tag{3}
\]
The expansion and its \(X\) derivatives hold uniformly on compact
disks. It follows by Taylor expansion of the heat operator on
the finitely many frequencies, followed by expansion of the local
power series at \(\xi_j\). A truncation after the power
\(K\) of \(T-T_*\) has error \(O(\epsilon^{2K+2})\) before
division, so \(K\) can be chosen for every fixed multiplicity.
These polynomials are rescaled Hermite polynomials, with simple
real roots for \(a<0\) and simple purely imaginary roots for \(a>0\).
For the real-root assertion use orthogonality with a positive
Gaussian weight, which follows by integrating the Rodrigues
formula. Purely imaginary roots follow by substituting \(X=iY\).

At a root \(h\) of \(P_m\) the heat identities
\(P_{\ell}'=\ell P_{\ell-1}\),
\(P_{\ell+1}=(X+2a\partial_X)P_\ell\), and
\(2aP_m''+XP_m'-mP_m=0\) give
\[
 \frac{P_{m+1}(h)}{P_m'(h)}=2a,\quad
 P_{m+1}'(h)=0,\quad
 \frac{P_{m+2}(h)}{P_m'(h)}=2ah,\quad
 \frac{P_m''(h)}{P_m'(h)}=-\frac{h}{2a}.
\]
The analytic root of (3) near \(h\) is therefore
\(h-2aA_j\epsilon-2ah(B_j-A_j^2/2)\epsilon^2+O(\epsilon^3)\).
On the negative-clock side the angles are
\[
          \vartheta_*+\epsilon^2 V_j
                +\epsilon^3 h_{j,r}
                    \frac{\operatorname{Im}(B_j-A_j^2/2)}j
                +O(\epsilon^4).                              \tag{4}
\]
On the positive-clock side their increasing list is
\(\vartheta_*+\epsilon h_{j,r}+O(\epsilon^2)\).
If \(V_n\ne V_N\), the two clusters are strictly ordered
on the negative side. On the positive side each multiple
cluster has both negative and positive leading displacements.
A simple cluster has zero leading displacement. Thus some
ranked difference takes opposite strict signs. When the
velocities agree, (4) gives the same reversal under (1).
Regular finite approximation now yields a transverse contact
between the two strengths. Letting \(\epsilon\) decrease to zero
gives the claimed convergence, with the physical scale addressed
below.

## Equal nonzero velocities

Suppose \(V_n=V_N=V\ne0\). Use the
[reflected pair compensation](MULTIVARIATE_MULTIPLE_ROOT_INTERIOR_HEAT_UNFOLDING.md#exact-reflected-row-compensation)
with a fixed \(0<w<T_*-|\lambda_*|\). Take \(r_M=\lfloor w M/2\rfloor\)
reflected unit-gap pairs at half-integral midpoint \(M\) and put
\(u=M^{-1}\), \(w_M=2r_Mu\), and \(\epsilon^2=2r_Mu^2\).
At effective strength \(T=T_*+\epsilon^2 y\) the actual Gaussian
strength is \(T-w_M\). In the scaled \(X\) coordinate the leading
operators and polynomials are
\[
 \begin{split}
     S_j(y)&=a_j(y)\partial_X^2+b_j\partial_X^4,\qquad
     a_j(y)=\frac{y-1/(2j)}{2j},\qquad
     b_j=\frac1{12j^3w},\\
     Q_{j,\ell}(X,y)&=e^{S_j(y)}X^\ell .
 \end{split}
\]
The local normalized polynomial has expansion
\[
                    Q_{j,m_j}+\epsilon A_jQ_{j,m_j+1}
                                    +O(\epsilon^2).          \tag{5}
\]
Here the error and its local \(X\) derivatives are uniform for \(y\)
in a fixed compact interval. The improvement to a first correction
in (5) follows from reflection parity. In the paired row log
multiplier the order-\(u^h\) term is even in \(k-j/2\), of degree at
most \(2\lfloor(h+1)/2\rfloor\). A term of even degree \(d\)
from that order acts after dilation with coefficient
\(O(r_Mu^h\epsilon^{-d})=O(u^{h-1-d/2})\). For \(h\ge4\)
this is \(O(u)\). The quadratic contribution at \(h=3\)
also has that bound, and \(w_M-w=O(u)\) affects the limiting
quartic coefficient by \(O(u)\). Since \(u=O_w(\epsilon^2)\),
the first correction in (5) comes only from the first higher root
jet \(A_j\). Fixed sufficiently long Taylor truncations of the log
and its exponential on the finitely many original frequencies
justify the expansion as in (3).

For \(y\) sufficiently negative, \(Q_{j,m_j}\) has simple real roots.
Indeed as \(a_j\to-\infty\), dilation by \(\sqrt{-a_j}\)
gives convergence to \(e^{-\partial_X^2}X^{m_j}\).
Reality of the coefficients and simplicity of the limiting roots
give the stated real roots. At a root \(h\) of \(Q_{j,m_j}\),
commuting multiplication by \(X\) with \(e^{S_j}\) gives
\[
 \frac{Q_{j,m_j+1}(h,y)}{Q'_{j,m_j}(h,y)}
       = 2a_j+4b_j\frac{Q'''_{j,m_j}(h,y)}{Q'_{j,m_j}(h,y)}
       = 2a_j+O_{w,n,N}(|y|^{-1}).
\]
The derivative ratio estimate follows from the same dilation and the
nonzero derivative of the limiting Hermite polynomial at each root.
The angle from (5) consequently equals
\[
       \vartheta_*+\epsilon^2 V(1/(2j)-y)
                    +\epsilon^2 O_{w,R,n,N}(|y|^{-1})
                    +O_{w,R,n,N,y}(\epsilon^3).
\]
Choose a fixed sufficiently negative \(y\), so the difference of
the two displayed \(V\)-terms dominates both \(O(|y|^{-1})\)
errors. The two angular clusters are then strictly ordered.
For fixed sufficiently large positive \(y\), the analogous dilation
by \(i\sqrt{a_j}\) gives simple purely imaginary roots of
\(Q_{j,m_j}\). A multiple cluster again has both signs of leading
imaginary displacement, and a simple cluster has zero displacement.
Thus some ranked difference reverses. For a countable core,
truncate with coefficient errors \(o(\epsilon^{\max(m_n,m_N)+2})\)
and center errors tending to zero. This preserves (5) and the
strict signs. Add the arbitrarily weak ordinary perturbation of
the first part of the proof to ensure transversality. For each
fixed \(w\) the resulting strengths approach \(T_*-w\).
A sequence of positive \(w\) decreasing to zero, with the
ordinary midpoints chosen sufficiently large for each \(w\),
gives convergence to \(T_*\).

## Physical recovery

At every angular contact choose the bias
\[
       \lambda_{\rm fin}=2nN\left(
           \Gamma_{\rm fin}-\log|z_{N,\rm fin}/z_{n,\rm fin}|
                         \right),\qquad
       \Gamma_{\rm fin}=\Delta_0\log(N/n)+\sum_{R_{\rm fin}} E .
\]
The exact retained center contributions \(E\) use the
[compactified center formula](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data).
They cancel in each reflected pair. The additional
genericity rows have center corrections tending to zero.
Keep \(\Delta_0\) fixed by assigning the difference between
\(\Delta_0\) and the retained signed gap to the residual
Gaussian data. Since the selected logarithmic disks can
shrink to the original roots, \(\lambda_{\rm fin}\to\lambda_*\).
The stated restriction on \(w\) and the original strict
interior condition give \(|\lambda_{\rm fin}|<T_{\rm fin}\).
The transferred variance agrees with the added curvature up
to \(o(1)\), and all genericity rows have curvature tending to
zero. Thus the total curvature converges to the original
total curvature. This finishes the finite-core reduction.

The theorem is a reduction for compactified Gaussian contacts.
It leaves the exact ordinary exclusion near these contacts open.
When both velocities vanish, (1) is only a sufficient test.
