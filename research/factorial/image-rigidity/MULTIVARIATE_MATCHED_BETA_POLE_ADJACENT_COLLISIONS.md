# Physical consecutive collisions for matched Beta pole pairs

For each terminal \(N\ge4\), there is an integer gap \(d\ge1\)
and a real algebraic \(\gamma>1\) with the following property.
Take independent \(X,Y\) with law \(\operatorname{Beta}(\gamma,1)\)
and the two monomials
\[
                         U=XY^{1+d},\qquad V=X^{1+d}Y.
\]
The moment polynomials \(\mathbb E(U+zV)^{N-1}\) and
\(\mathbb E(U+zV)^N\) have a common nonreal unit-circle root
at the same physical coefficient ratio \(z\).
The exponents are fixed along the path in \(\gamma\).
The tilted log-ratio laws are Pólya-frequency of infinite order
and satisfy both folded likelihood-ratio comparisons at all
adjacent orders throughout the path.

For each \(N\ge4\) the same physical collision occurs with two
independent Gamma inputs, common positive shape, and positive
real exponents on the corresponding two reflected rows. The
Gamma construction uses arbitrarily small real exponent scales
and a proportional common shape. It does not give an ordinary
factorial profile with integer exponents and unit Gamma shapes.

## Exact coefficient coupling and the folded law

More generally take shared integer lower exponent \(\ell\ge1\)
and integer gap \(d\ge1\) in the reflected row pair. For Beta
inputs write \(Q_j\) for its moment polynomial divided by its
constant coefficient. Its two endpoints equal one, and
\[
 [z^k]Q_j(z)=\binom jk\,\frac{X_j}{X_j+k(j-k)},\qquad
 B_j=\frac{\gamma+j\ell}{d},\qquad X_j=B_j(B_j+j).
                                                               \tag{1}
\]
Indeed \(\mathbb E X^a=\gamma/(\gamma+a)\); at index \(k\) the
product of the two denominators is
\(d^2(B_j+k)(B_j+j-k)\).
In particular \(B_{j+1}-B_j=\ell/d\). This is the same affine
order schedule as a single Gauss parameter \((j\ell+h)/d\),
with fixed shift \(\gamma/d\). Formula (1) uses one reciprocal
linear factor per orientation. The ordinary factorial symbol
instead uses the Gamma lattice of rising factors.

Under the geometric mean tilt by \((UV)^{j/2}\), the logarithms of
the two independent Beta variables are negatives of independent
exponentials of rate \(\gamma+j(\ell+d/2)\). Consequently
\(S=\log(V/U)\) under that tilt has symmetric Laplace density
\[
              \tfrac12\lambda_j e^{-\lambda_j|s|},
      \qquad \lambda_j=B_j+j/2 .
\]
A Pólya-frequency density of infinite order is one whose translation
kernel has nonnegative ordered minors of every size. Each one-sided
exponential density has this property. After extracting positive
row and column factors from its minors, the kernel is a threshold
matrix \(1_{x\ge y}\), or its reflection. Its rows on increasing
grids have nondecreasing lengths of initial ones. A square matrix
of this form has determinant zero unless the row lengths are
\(1,\ldots,h\), in which case it is triangular with determinant
one. Reflection reverses both grids. The two exponential kernels
compose to the Laplace kernel, and the integral Cauchy–Binet
identity preserves the minor signs.

The law of \(|S|\) is exponential of rate \(\lambda_j\). Since
\(\lambda_j<\lambda_{j+1}<((j+1)/j)\lambda_j\), its folded laws
obey
\[
       (j/(j+1))\mu_j^+\preceq_{\rm LR}
                             \mu_{j+1}^+\preceq_{\rm LR}\mu_j^+.
\]
Here dilation acts on the positive variable, and the ratio of the
density of a smaller-rate exponential to that of a larger-rate
one is increasing. Thus both comparisons used in the
[folded-law comparison](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
hold along the whole path. This construction involves a single
shared pair, with no added private pair to remove.

## An exact terminal-four collision

Take \(\ell=d=1\), and abbreviate \(X=X_3\) and \(Y=X_4\).
The order-three polynomial has its factor \(1+z\) and a
reciprocal quadratic whose nontrivial unit-circle pair corresponds
to
\[
                          x_3=\frac{2(1-X)}{X+2}\in(-2,2)
\]
under \(x=z+z^{-1}\).
The order-four polynomial divided by \(z^2\) becomes
\[
                E(x)=x^2+\frac{4Y}{Y+3}x+\frac{6Y}{Y+4}-2.
\]
Substituting \(x_3\) gives
\[
 E(x_3)=\frac{12(2X^2-6XY-16X+3Y^2+9Y-4)}
                         {(X+2)^2(Y+3)(Y+4)}.
\]
For \(X=(\gamma+3)(\gamma+6)\) and \(Y=(\gamma+4)(\gamma+8)\)
the numerator inside the factor \(12\) is
\[
                  -(\gamma+5)(\gamma^3+13\gamma^2+32\gamma-52).
\]
It equals \(36\) at \(\gamma=1\) and \(-504\) at \(\gamma=2\).
The positive root of the cubic in \((1,2)\) gives the collision.
The map \(x_3=z+z^{-1}\) yields the common nonreal pair on the
unit circle.

For comparison, terminal two is separated for all these positive
laws by real-axis positivity of the even moment. At terminal
three, a common nonreal pair would make the monic order-two
polynomial equal to the quadratic obtained from order three by
dividing by \(1+z\). By (1), equality of middle coefficients requires
\(X_3=3X_2+1\). But \(B_3\le(3/2)B_2\), so
\(X_3\le(9/4)X_2\). Thus terminal four is the first possible
terminal for this Beta family.

## The start of the path at every higher terminal

Fix \(n=N-1\ge4\), use \(\ell=1\), and take
\[
                            d>4\cdot 2^N(N+1)^2nN .
\]
At \(\gamma=1\), (1) gives \(X_j\le2(j+1)^2/d\) for \(j=n,N\).
On the unit circle the real centered polynomial
\(e^{-ij\theta/2}Q_j(e^{i\theta})\) differs from
\(2\cos(j\theta/2)\) in absolute value by at most
\(2^j X_j/(j-1)\), since \(k(j-k)\ge j-1\) at every
interior index. Put \(\varepsilon=\pi/(4nN)\).
At angular displacement \(\pm\varepsilon\) from any binomial
root of degree \(j\), the absolute value of the latter centered
function is \(2\sin(j\varepsilon/2)\ge j/(2nN)\).
The displayed choice of \(d\) makes the perturbation strictly
smaller. Thus \(Q_j\) has one root in each short arc about a
root of \(1+z^j\), and the degree count makes these roots simple.
The degree-\(n\) and degree-\(N\) binomial arguments in the
positive cut \(0<\theta<2\pi\) alternate, with endpoints of
degree \(N\). Any two arguments from different degrees differ
by at least \(\pi/(nN)\). The certified arcs therefore preserve
the strict alternation.

## Finite loss of the circle at every degree \(j\ge5\)

Put \(z=(w-1)/(w+1)\), so circular roots away from \(z=1\) become
purely imaginary \(w\). The identity
\[
 e^{-S/2}\bigl((w+1)+(w-1)e^S\bigr)
                    =2\bigl(w\cosh(S/2)-\sinh(S/2)\bigr)
\]
shows that the Cayley polynomial, up to a positive constant,
has coefficients
\[
    \sum_{r=0}^{\lfloor j/2\rfloor}
       \binom j{2r} w^{j-2r}\,\widehat\mu_{2r},
 \qquad
 \widehat\mu_k=
 \frac{\mathbb E_{\rho_j}\!\left[
          \cosh(S/2)^j\tanh(S/2)^k\right]}
      {\mathbb E_{\rho_j}\cosh(S/2)^j}.
                                                               \tag{2}
\]
Here \(\rho_j\) is the tilted Laplace law. Odd moments vanish
by reflection. The expectations exist since \(\lambda_j>j/2\).

If all roots of (2) were imaginary, after removing its factor \(w\)
for odd \(j\) the polynomial in \(y=w^2\) would have all nonpositive
roots. Write \(m=\lfloor j/2\rfloor\). Comparing its first two
elementary symmetric functions by Cauchy–Schwarz gives
\[
 \left(\binom j2\widehat\mu_2\right)^2
                   \ge\frac{2m}{m-1}\binom j4\widehat\mu_4.
\]
Indeed for \(m\) nonnegative numbers \(r_i\), one has
\((\sum r_i)^2\ge (2m/(m-1))\sum_{i<l}r_i r_l\).
Thus the necessary moment ratio bound is
\[
 \frac{\widehat\mu_4}{\widehat\mu_2^2}\le \kappa_j,\qquad
 \kappa_j=
 \begin{cases}
      3(j-1)/(j-3),&j\ \text{even},\\
      3j/(j-2),&j\ \text{odd}.
 \end{cases}                                                    \tag{3}
\]
In degrees \(j\ge5\) one has \(\kappa_j\le5\).

There is an exact rational violation throughout
\(\lambda_j\ge j^3\). To see this set \(q=\lambda_j^{-2}\),
\(A=j/2\), and \(C_b=\mathbb E_{\rho_j}\cosh(S/2)^b\) for
\(b=j,j-2,j-4\). The Laplace integral gives
\[
       C_b=2^{-b}\sum_{h=0}^b\binom bh
                         \frac{1}{1-q(h-b/2)^2}.
\]
The coefficients at orders \(q,q^2\) are \(b/4\) and
\((3b^2-2b)/16\), the second and fourth moments of a sum of
\(b\) independent symmetric signs divided by two. The remainder
after order two lies between zero and \(2q^3 A^6\) since
\(qA^2\le1/2\). The geometric remainder after order one also gives
\[
 \begin{split}
 C_j-C_{j-2}&\le q/2+2A^4q^2,\\
 C_j-2C_{j-2}+C_{j-4}&\ge 3q^2/2-4A^6q^3>0 .
 \end{split}
\]
The identity
\(\widehat\mu_2=(C_j-C_{j-2})/C_j\) and its fourth-moment
counterpart give the violation. In fact, using \(C_j\ge1\),
\[
 \begin{split}
 C_j(C_j-2C_{j-2}+C_{j-4})
           -\kappa_j(C_j-C_{j-2})^2
 \ \ge\
 q^2\left(\frac{6-\kappa_j}{4}
             -4A^6q-2\kappa_j A^4q-4\kappa_j A^8q^2\right)>0 .
 \end{split}
\]
For \(q\le j^{-6}\), the three subtractions in parentheses
are at most \(1/16,\,1/40,\,1/8000\), respectively, while the
positive term is at least \(1/4\). This proves (3) fails
without any root computations.

Take, for example, \(\gamma_1=dN^3\) for the pair of orders
\(n,N\). Then \(\lambda_N\ge N^3\), and the upper polynomial has
an off-circle root at \(\gamma_1\). Starting at the strictly
alternating circular roots at \(\gamma=1\), consider the first
parameter where either polynomial ceases to have all its roots
simple and circular. Both properties persist locally as long as
all circular roots are simple, since reciprocity and conjugation
fix each nearby simple branch under \(z\mapsto1/\bar z\).
Roots depend continuously on the coefficients. At the first
failure some two consecutive roots of the same degree coalesce
on the circle. If there has already been a common root the claim
holds. Otherwise the intervening root of the other degree is
squeezed to the same point, giving a common root at the failure.

The angular cut at \(1\) cannot cause an exception, since all
coefficients are positive there. An even moment is strictly positive
at all real coefficient ratios. For odd degree the derivative on
the real axis is \(j\mathbb E[V(U+zV)^{j-1}]>0\), so its real
root at \(-1\) stays simple. The contact is therefore nonreal.
All coefficients in (1) are rational functions of \(\gamma\)
with positive denominators on the path. The initial interlacing
makes the resultant nonzero there, so any contact parameter is
a root of a nonzero polynomial over \(\mathbb Q\).

## Keeping the matched Gauss shift fixed

The same collisions can be obtained with fixed uniform inputs
\(\gamma=1\) by varying the positive real lower exponent \(\ell\).
In this version \(B_j=(j\ell+1)/d\) is precisely the first Gauss
parameter of the exponent row for every degree, with fixed
shift \(1/d\). At terminal four, take \(d=1\). The numerator
inside the factor \(12\) in the evaluation of \(E(x_3)\) is
\(-36\) at \(\ell=0\) and \(36\) at \(\ell=1\), giving an
interior positive lower exponent. At every higher terminal use
the same large integer gap as above. The initial alternation
at \(\ell=1\) persists until a common root, or until the first
loss of simplicity and circular location gives one. Such a
loss must occur by \(\ell=dN^2\), since then
\(\lambda_N\ge N^3\). Thus the obstruction also applies when the
paired affine coordinate is kept exactly on the first-shift path.
The coefficients still come from the reciprocal factors (1).

## Any fixed finite number of matched poles

Fix an integer \(M\ge1\). At every terminal \(N\ge2M+3\) the
same physical collision occurs for independent
\(\operatorname{Beta}(\gamma,M)\) inputs on fixed shared integer
rows \((1,1+d),(1+d,1)\), for a positive algebraic shape
\(\gamma>1\). It also occurs with fixed \(\gamma=1\) and
varying positive real lower exponent \(\ell\), using an
integer \(d\ge M\). In that version all \(M\) pole shifts are
exactly the first \(M\) matched shifts \(h/d\), \(1\le h\le M\),
at every order.

The moment \(\mathbb E X^a=(\gamma)_M/(\gamma+a)_M\)
follows from the Beta integral. Thus the normalized coefficient
at \(k\), after the binomial factor is removed, is a product of
\(M\) factors of the form in (1), with
\[
                B_{j,r}=\frac{\gamma+j\ell+r}{d},
                       \qquad 0\le r<M .
\]
All the factors are at most one. At \(\gamma=\ell=1\) the
first factor alone gives the same strict initial arc signs as
the single-pole proof, using the same large-gap condition.

At order \(j\) the tilt by \((UV)^{j/2}\) changes the first
Beta shape to \(\gamma+j(\ell+d/2)\). The log ratio \(S\)
then has bilateral moment transform
\[
                  \prod_{r=0}^{M-1}\frac{\lambda_{j,r}^2}
                                          {\lambda_{j,r}^2-u^2},
     \qquad
       \lambda_{j,r}=\frac{\gamma+j(\ell+d/2)+r}{d}.
\]
Indeed the transform of minus \(d\) times the logarithm of
one tilted Beta variable is the product of the corresponding
one-sided exponential transforms. Hence \(S\) is a sum of
independent symmetric Laplace variables of these rates.
In particular its translation kernel still has nonnegative
ordered minors of every size.

The folded likelihood-ratio comparisons also compose here. For
an even Laplace density \(f_\lambda\), the folded convolution
kernel on \(x,y\ge0\) is
\[
 f_\lambda(x-y)+f_\lambda(x+y)
    =\lambda e^{-\lambda(x+y)}
          H_\lambda(\min(x,y)),\qquad
     H_\lambda(u)=e^{\lambda u}\cosh(\lambda u).
\]
For \(x_1<x_2\), the ratio of the row at \(x_2\) to the row
at \(x_1\) is nondecreasing in \(y\), because \(H_\lambda\)
is increasing and the interval where exactly one minimum varies
is \(x_1<y<x_2\). Thus its order-two minors are nonnegative.
For positive folded densities \(p\preceq_{\rm LR}q\), expansion
of the transformed order-two comparison gives the integral over
\(y_1<y_2\) of the kernel minor times
\(p(y_1)q(y_2)-q(y_1)p(y_2)\), which is nonnegative. Folding
the full sum of even variables therefore preserves the order
as one factor is replaced at a time. Each factor has the rate
inequalities
\(\lambda_{j,r}<\lambda_{j+1,r}<(j+1)\lambda_{j,r}/j\).
Apply the folded convolution comparison to each inequality,
using that dilation commutes with convolution. This proves the
full folded sandwich for the \(M\)-pole sum.

At fixed \(j,d,\ell\), write
\[
      S=\sum_{r=0}^{M-1}\frac{L_r}{\lambda_{j,r}},
\]
where the \(L_r\) are independent standard symmetric Laplace
variables. As \(\gamma\to\infty\), the scaled variable
\((\gamma/d)S\) tends to \(\sum_r L_r\) in this coupling
and is bounded in absolute value by \(Z=\sum_r|L_r|\).
For \(\gamma\ge jd\), the positive tilt factor
\(\cosh(S/2)^j\) is bounded by \(e^{Z/2}\), while
\((2\gamma/d)|\tanh(S/2)|\le Z\).
The products \(Z^k e^{Z/2}\) are integrable for \(k=0,2,4\)
because the \(|L_r|\) are independent rate-one exponentials.
Dominated convergence therefore applies to the tilted moments
in (2), including their normalizing denominator.
A standard Laplace variable has second moment two and fourth
moment \(24\), so the sum has second moment \(2M\) and fourth
moment \(24M+12M(M-1)=12M(M+1)\). Consequently
\[
             \lim_{\gamma\to\infty}
                       \frac{\widehat\mu_4}{\widehat\mu_2^2}
                         =3+\frac3M .
\]
In degree \(j=N\) this strictly exceeds (3) when \(N\ge2M+3\):
for odd \(N\) the requirement is \(N-2>2M\), and for even
\(N\) it is \(N-3>2M\). The first-loss argument gives the
common physical root for a finite shape. The initial strict
alternation makes the rational resultant nonzero, so the
contact shape is algebraic.

For fixed \(\gamma=1\) and \(\ell\to\infty\), the same reasoning
works with \(j\ell/d\) in place of \(\gamma/d\). The rates then
retain the exact shift \(h/d\), \(h=1,\ldots,M\), in their
\(B_{j,h-1}\) coordinates; \(d\ge M\) ensures that these are
Gauss shifts of the exponent row. Thus a fixed finite truncation
of the matched poles, even in exact paired coordinates, cannot
give uniform consecutive separation.

## Passage to Gamma inputs with real exponents

For a fixed pair \(\ell,d\) above, use exponents
\(t(\ell,\ell+d)\) and their reflection with \(t>0\), and
independent unit-rate Gamma inputs of common shape \(t\gamma\).
The binomial-stripped normalized coefficient of order \(j\)
at index \(k\) is
\[
 \frac{\Gamma(t(\gamma+j\ell+dk))
             \Gamma(t(\gamma+j(\ell+d)-dk))}
      {\Gamma(t(\gamma+j\ell))\Gamma(t(\gamma+j(\ell+d)))}.
\]
It tends to the Beta coefficient in (1) as \(t\downarrow0\),
because \(\Gamma(x)\) has a simple pole of residue one at
zero. [DLMF Gamma properties](https://dlmf.nist.gov/5.2)

For \(N\ge5\), at the two chosen Beta endpoint parameters the
initial simple alternation and the final off-circle upper root
persist for sufficiently small \(t>0\). The polynomials remain
reciprocal, and their positive moment law gives the same real-axis
exclusions. The first-failure argument gives a physical Gamma
contact at an intermediate \(\gamma\) for each such fixed \(t\).
For \(N=4\) the strict endpoint sign reversal persists under the
same Gamma approximation. The nontrivial order-three pair remains
circular: its reciprocal quadratic factor has no real root,
since the odd moment has exactly one simple real root. Its
continuous \(x\)-coordinate therefore supplies the same
intermediate-value test against the order-four quotient.
For the \(M\)-pole version the same approximation applies
factor by factor. Use \(M\) independent mirrored Gamma pairs
with shapes \(t(\gamma+r)\), \(0\le r<M\), and the same scaled
exponent rows in each pair.

The Beta construction keeps integer exponent rows and their exact
joint order schedule, with a single pole in each orientation. The
Gamma construction uses the full Gamma function with proportional
small shapes and fractional exponents. In the ordinary factorial
class the shapes are one and the exponents are nonnegative
integers. The adjacent collisions above distinguish the larger
positive-measure classes from that arithmetic setting.
