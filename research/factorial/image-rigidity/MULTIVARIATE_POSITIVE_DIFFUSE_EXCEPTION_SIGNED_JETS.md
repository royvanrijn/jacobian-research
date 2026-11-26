# Matched factorial tests at positive diffuse exceptional strengths

Fix an ordinary two-monomial core and one adjacent pair
\(n,N=n+1\), \(n\ge2\). The
[physical diffuse exception theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
gives a finite set of candidate contacts at positive Gaussian
strength, each with algebraic clock \(q_*\), algebraic relative
bias factor \(v_*\), and integer signed tail gap \(h\).
This note gives two exact local tests that can exclude a
neighborhood of such a contact for arbitrary ordinary diffuse
tails. Both tests concern simple reference roots. One uses the
boundary of the ordinary signed bias cone and the other uses
a stationary Gaussian angular gap.

## A matched expansion with a second-order ordinary remainder

For added unequal shared rows put \(m_i=(a_i+b_i)/2\),
\(d_i=b_i-a_i\), and \(\eta=\max_i|d_i|/m_i\). Define
\[
 \begin{gathered}
 W=\sum_i d_i^2/m_i,\quad \lambda=\sum_i d_i/m_i,\quad
 K=\sum_i d_i^3/m_i^2,\quad Z=\sum_i d_i^2/m_i^2,\\
 B_4=\sum_i d_i^4/m_i^3,\qquad L_2=\sum_i d_i/m_i^2,\qquad
 h=\sum_i d_i .
 \end{gathered}
\]
The midpoint quantities \(W,\lambda\) have the same limits as
the leading strength and signed bias using the larger exponent,
since the differences are at most \(\eta W/2\) in modulus.
The ordinary constraints include
\[
     |\lambda|\le W,\quad |K|\le\eta W,\quad Z\le\eta W,\quad
         B_4\le W Z,\quad K^2\le W B_4,\quad |L_2|\le Z .
                                                               \tag{1}
\]
The third bound uses \(1/m_i\le |d_i|/m_i\). The fourth compares
each summand of \(B_4\) with \(Z_i d_i^2/m_i\), and the fifth
is Cauchy–Schwarz with weights \(d_i^2/m_i\).

Use \(t=s-1/2\) and let \(L_j(s)\) denote the normalized tail
log coefficient per degree, with zero endpoint values.
Uniformly for \(W\) in a fixed compact positive interval,
the exact curvature is
\[
       L_j''(s)=W-tK+t^2B_4-Z/(2j)+O_j(\eta Z),
                           \qquad 0\le s\le1.                 \tag{2}
\]
Indeed, with \(x=m_i+t d_i\) the summand is
\(d_i^2j\psi'(1+jx)\). Expand
\(j\psi'(1+jx)=1/x-1/(2jx^2)+O_j(x^{-3})\)
using the trigamma integral. The geometric remainder beyond
the quadratic term in the expansion of \(1/x\) costs at most
a constant times \(\eta B_4\). The variation of \(1/x^2\)
costs at most a constant times \(\eta Z\). The remaining
trigamma term costs \(O_j(\sum_i d_i^2/m_i^3)=O_j(\eta Z)\).
For example \(x\ge3m_i/4\) when \(\eta\le1/2\), so all
constants follow from fixed rational bounds on these denominators.
The bound \(B_4\le W Z\) gives (2).

Let \(F_j,G_j,H_j\) be the index polynomials in the
[weak-tail response](MULTIVARIATE_WEAK_TAIL_SECOND_ORDER_CONTACT_GATE.md#exact-second-order-response-of-a-simple-core-root)
and put
\[
                I_j(k)=\frac{j}{12}\left((k/j-1/2)^4-1/16\right).
\]
Integrating (2) with zero endpoints gives the log multiplier
\(WF_j(k)+K H_j(k)+Z G_j(k)+B_4 I_j(k)+O_j(\eta Z)\).
The fourth-power term is the zero-endpoint primitive of \(t^2\).
This term can be the same size as \(Z\) at positive Gaussian
strength and is kept in the test.

The ratio of the two tail argument scales \(D_N/D_n\) satisfies
\[
 \log(D_N/D_n)=h\log(N/n)-\frac{\lambda}{2nN}
                            +c L_2+O_{n,N}(\eta Z),\qquad
        c=\frac{n^{-2}-N^{-2}}{12}.                            \tag{3}
\]
Here \(D_j\) normalizes the tail alone, so the raw core
argument is otherwise unchanged. To check the further term,
the exact integral for the order derivative of a row's log
scale is \(\int_{a_i}^{b_i}x\psi'(1+jx)\,dx\), with signed
endpoints. The trigamma expansion through its next term gives
\[
 x\psi'(1+jx)=j^{-1}-(2j^2x)^{-1}+(6j^3x^2)^{-1}+O_j(x^{-3}).
\]
The scalar bound for this expansion follows by expanding
\(u/(e^u-1)=1-u/2+u^2/12+O(u^3)\) in the Laplace integral.
An absolute constant in the \(O(u^3)\) works on all
\(u\ge0\): for \(0\le u\le1\) divide the Taylor expansion
of \((e^u-1)/u\), with denominator at least one, and for
\(u\ge1\) use \(0<u/(e^u-1)\le1\).
Integrating in \(j\) gives the coefficients of \(1/x\)
and \(1/x^2\) in (3). In the signed integrals in \(x\),
replace these by \(d_i/m_i\) and \(d_i/m_i^2\), respectively.
The first midpoint error is \(O(|d_i|^3/m_i^3)\), because
its linear term integrates to zero. The second is
\(O(d_i^2/m_i^3)\) by the mean-value bound.
All errors sum to \(O(\eta Z)\), proving (3).

## Root coefficients at an algebraic exceptional clock

Use the raw Gaussian core polynomials
\[
       Q_j^G(w;W)=\sum_{k=0}^j [w^k]Q_j(w)\,e^{WF_j(k)}w^k .
\]
Let \(\zeta_j(W)\) be the simple root branch through a chosen
reference contact at \(W_*>0\). Use brackets for operator
evaluation at \(\zeta_j(W)\) on \(Q_j^G(\cdot;W)\), let
\(\mathcal D=w\,d/dw\), and set \(q_1=[\mathcal D]\ne0\),
\(q_2=[\mathcal D^2]\). Write
\[
 \begin{split}
 \alpha_j&=-[F_j]/q_1,\quad
 \kappa_j=-[H_j]/q_1,\quad
 \gamma_j=-[G_j]/q_1,\quad \sigma_j=-[I_j]/q_1,\\
 \beta^H_j&=-\bigl(q_2\kappa_j^2/2+
                      \kappa_j[\mathcal D H_j]+[H_j^2]/2\bigr)/q_1 .
 \end{split}
\]
These are functions of \(W\); \(\alpha_j\) is the logarithmic
root derivative in Gaussian strength.
They and every derivative needed below are algebraic at
\(W=W_*\) when \(q_*=e^{-W_*/(2nN)}\) is algebraic.
Indeed the reference coefficients are integer powers of
\(q_*\) times raw core coefficients, and implicit derivatives
at a simple algebraic root involve only algebraic operations.

If \(w_j\) is the nearby exact root after normalizing the tail,
the root perturbation formula is
\[
 \log(w_j/\zeta_j(W))
       =K\kappa_j+Z\gamma_j+B_4\sigma_j+K^2\beta^H_j+O(\eta Z).
                                                               \tag{4}
\]
Expand the coefficient multiplier of (2) around its Gaussian
part and solve the local equation in the logarithmic root
coordinate, as in (2)–(3) of the weak-tail proof. The coefficient
of the unknown displacement is \(q_1\ne0\). The term
\(K^2\beta^H_j\) comes from the second response to the skew
operator. The remaining products of the parameters are
\(O(\eta Z)\) by (1), including
\(|K|^3\le\eta W^3 Z\). Uniform effective constants come from
isolating the simple reference roots on disjoint disks
avoiding zero, applying Rouché on their boundaries, and
bounding the analytic root map on smaller coefficient
polydisks by Cauchy's formula. The reference parameter can
be varied in a sufficiently small fixed interval about
\(W_*\) with all these bounds uniform.

Write \(\Delta f=f_N-f_n\) for differences of root responses,
and choose the local logarithm of \(\zeta_N(W)/\zeta_n(W)\)
which is real at the contact. With
\[
 g_{\mathbb C}(W)=\log(\zeta_N(W)/\zeta_n(W))
                    -h\log(N/n)+\lambda_* /(2nN),
\]
we have \(g_{\mathbb C}(W_*)=0\) and
\(g_{\mathbb C}'(W_*)=\Delta\alpha(W_*)\).
A physical contact near the chosen simple pair must satisfy
\[
 0=g_{\mathbb C}(W)+\frac{\lambda-\lambda_*}{2nN}
     +K\Delta\kappa+Z\Delta\gamma+B_4\Delta\sigma
                   +K^2\Delta\beta^H-c L_2+O(\eta Z).          \tag{5}
\]

## Transverse contacts on the signed bias boundary

Suppose \(\lambda_*=\epsilon W_*\) with \(\epsilon=\pm1\) and
\(\Im\Delta\alpha(W_*)\ne0\). Define the algebraic numbers
\[
           p=\Delta\alpha(W_*)+\epsilon/(2nN),\qquad
           J_*=\epsilon\Delta\kappa(W_*)+\Delta\gamma(W_*)-\epsilon c .
\]
A sufficient condition to exclude the contact for all
sufficiently diffuse tails in a neighborhood is
\[
 \epsilon\left(\Re J_*-\frac{\Re p}{\Im p}\Im J_*\right)<0 .
                                                               \tag{6}
\]
This condition uses the ordinary unit-gap saturation of the
bias boundary.

To prove it, put \(\delta=W-\epsilon\lambda\ge0\).
A row with signed gap \(d_i=\epsilon\) contributes zero to
\(\delta\). All other rows contribute at least half their
curvature, since \(1-\epsilon/d_i\ge1/2\).
The sum of their curvatures is thus at most \(2\delta\).
On the aligned unit rows \(K_i=L_{2,i}=\epsilon Z_i\) and
\(B_{4,i}\le\eta Z_i\). As \(1/m_i\le\eta\), the other rows
give
\[
 K=\epsilon Z+O(\eta\delta),\qquad
 L_2=\epsilon Z+O(\eta\delta),\qquad
                  B_4=O(\eta Z+\eta^2\delta).
\]
Putting \(x=W-W_*\) in (5) and using
\(\lambda-\lambda_*=\epsilon(x-\delta)\) gives
\[
              0=p x-\epsilon\delta/(2nN)+J_* Z
                              +O(x^2+|x|Z+\eta(Z+\delta)).
\]
Here the \(K^2\) term is part of the error, since \(Z\le\eta W\)
and \(\delta\le2W\). The imaginary equation, using
\(\Im p\ne0\), first gives \(|x|=O(Z+\eta\delta)\) by choosing
a small neighborhood and absorbing the quadratic term.
Eliminate \(x\) from the real equation to obtain
\[
 \delta=2nN\epsilon
       \left(\Re J_*-\frac{\Re p}{\Im p}\Im J_*\right) Z
                          +O(\eta(Z+\delta)).
\]
Since \(Z>0\) for any nonempty tail and \(\delta\ge0\), the
strict negative sign in (6) excludes a contact. The estimates
give effective neighborhood and diffuseness tolerances, by
choosing \(\eta\) so the uniform error coefficient is smaller
than the strict sign margin and than one.

## Stationary contacts with a common phase sign

Let \(g(W)=\Im g_{\mathbb C}(W)\).
This analytic function has a finite order of vanishing at the
chosen contact. Otherwise the branches would share a ray on
an interval, and the off-diagonal Gaussian quotient resultant
would vanish identically as a polynomial in \(q\), contrary
to its nonzero specialization at \(q=0\).
Suppose the order is even, \(2r\), with leading coefficient
\(a\ne0\) in \(W-W_*\). Assume \(A(W)=\Im\Delta\kappa(W)\)
vanishes through order \(r-1\), and let \(b\) be its coefficient
of order \(r\), allowing zero. Set
\[
        B=\Im\Delta\gamma(W_*),\qquad
        C=\Im\Delta\sigma(W_*),\qquad
        D=\Im\Delta\beta^H(W_*).
\]
A sufficient exclusion condition, for arbitrary signed biases,
is
\[
 \operatorname{sgn}(a)
       \left[B+C u+\left(D-\frac{b^2}{4a}\right)v\right]>0
       \quad\text{when}\quad 0\le u\le W_*,\quad 0\le v\le W_*u .
                                                               \tag{7}
\]
For reciprocal core arrays \(A\) vanishes identically at
unit-circle references by reciprocal conjugation of the
skew perturbation, so \(b=0\).

Indeed at a putative contact the imaginary part of (5), with
(1), gives
\(|W-W_*|^{2r}=O(\sqrt Z\,|W-W_*|^r+Z)\).
Thus \(y=(W-W_*)^r/\sqrt Z\) stays bounded.
Put \(k=K/\sqrt Z\), \(u=B_4/Z\), and \(v=k^2\); by (1),
\(0\le u\le W\) and \(0\le v\le Wu\).
Divide the phase equation by \(Z\) to obtain
\[
                      0=a y^2+b k y+B+C u+D k^2+o(1)
\]
uniformly along any contact sequence. The boundedness of
\(y,k,u\), the local Taylor errors, and \(O(\eta Z)\) in
(5) give this remainder. Completing the square in \(y\)
proves exclusion under (7). The sign margin on the compact
triangle and Taylor bounds give effective uniform tolerances.
To decide (7) it suffices to check its three vertices
\((0,0),(W_*,0),(W_*,W_*^2)\). Their values are polynomials
of degree at most two in \(W_*\) over the algebraic numbers.
Since \(W_*=-2nN\log q_*\) is positive and transcendental
by Hermite–Lindemann, a nonzero such polynomial has a
decidable strict sign by rational interval refinement.
Identically zero polynomials fail the strict test.

These are exclusion tests for individual positive reference
contacts of the exact core. If every contact in a fixed compact positive
strength and bias range satisfies one of the tests, the
existing compact gap away from the finite exceptional set
and these local neighborhoods give uniform separation for
all sufficiently diffuse ordinary tails in that range.
For transverse interior contacts the ordinary relation
\(|\lambda|\le W\) permits both directions of displacement;
neither boundary sign test supplies an exclusion there.


## Exact coefficient and elimination controls

The following identities check the endpoint primitives, the
boundary normal elimination for both orientations, and the
stationary square completion.
~~~python
import sympy as sp

s, j = sp.symbols("s j", nonzero=True)
t = s-sp.Rational(1, 2)
F = j*s*(s-1)/2
G = -F/(2*j)
H = -j*t*s*(s-1)/6
I = j*(t**4-sp.Rational(1, 16))/12
for primitive, curvature in ((F, 1), (G, -1/(2*j)), (H, -t), (I, t*t)):
    assert primitive.subs(s, 0) == primitive.subs(s, 1) == 0
    assert sp.simplify(sp.diff(primitive, s, 2)/j-curvature) == 0
pr, pi, jr, ji, Z, c = sp.symbols("pr pi jr ji Z c", nonzero=True)
x, delta = sp.symbols("x delta")
for orientation in (-1, 1):
    sol = sp.solve([pi*x+ji*Z, pr*x-orientation*delta/c+jr*Z],
                   [x, delta])
    assert sp.cancel(sol[delta]-orientation*c*(jr-pr*ji/pi)*Z) == 0
a, b, B, C, D, u, y, k = sp.symbols("a b B C D u y k", nonzero=True)
assert sp.expand(a*y*y+b*k*y+B+C*u+D*k*k
                 -a*(y+b*k/(2*a))**2-B-C*u-(D-b*b/(4*a))*k*k) == 0
~~~
