# Exact angle boxes with private variables on both sides

Consider an ordinary two-monomial support with at least
one private variable in each monomial. Each of these
variables is absent from the opposite monomial. The
private exponents are arbitrary positive integers, and
the remaining rows may have either orientation.

**Theorem.** Every factorial moment polynomial has simple
roots. For degree \(m\), its \(k\)-th upper-half-plane root
has argument in
\[
        \left(\frac{(2k-1)\pi}{m},\
                 \frac{2k\pi}{m+1}\right),
                 \qquad 1\le k\le\lfloor m/2\rfloor,
                                                               \tag{1}
\]
unless the only unequal rows are one private unit row
on each side. Its other nonreal roots are the conjugates,
and there is exactly one negative real root in odd degree.
In the private-unit exception the roots are those of a
geometric sum, with arguments at the right endpoints of
(1).

Consequently, a common nonreal root at consecutive
orders \(n,n+1\) must have the same upper-half-plane
label in both degrees, after conjugation if necessary.
Its label \(k\) satisfies \(4k<n+2\), and its real part
is strictly positive. In particular all ordinary
profiles in this class have consecutive separation
in the closed left half-plane.

The angle boxes come from a positive integral
representation. They do not require the roots to
have a common modulus.

## The geometric resolvent without reflection symmetry

Reserve one private row of exponent \(d\ge1\) in each
monomial, allowing the two exponents to differ.
For each reserved row and every \(j\ge0\), the identity
\[
                 \frac{(dj)!}{j!}
                       =d^{dj}\prod_{h=1}^{d-1}(h/d)_j
\]
follows by grouping the factorial into residue classes.
The right side is the \(j\)-th moment of \(d^d\)
times a product of independent unit-scale Gamma
variables of shapes \(h/d\), with an empty product for
\(d=1\). Thus each reserved row supplies one \(j!\)
for the binomial cancellation and leaves a positive
independent random factor in its monomial.

Write \(U,V\) for the remaining monomial values,
including these residual Gamma products. The other
ordinary inputs are independent standard exponentials.
The two reserved factorials cancel the binomial
coefficient up to \(m!\). After normalization the
degree-\(m\) polynomial is
\[
                  P_m(z)=\mathbb E_m\sum_{j=0}^m(zT)^j,
               \qquad T=V/U,\qquad
         d\mathbb P_m=\frac{U^m}{\mathbb E U^m}\,d\mathbb P .
                                                               \tag{2}
\]
The ratio is constant exactly in the private-unit
exception stated above. Otherwise an unequal exponent
in an independent variable, or a nonconstant residual
Gamma product, makes its logarithm nonconstant.
Only positivity and a finite \(m\)-th moment of \(T\)
are needed in the proof below. Thus it also proves
(1) for every polynomial whose coefficients of indices
\(0,\ldots,m\) are the moments of a positive
nonconstant random variable.

For \(z=r e^{i\theta}\), \(r>0\), \(0<\theta<\pi\),
write \(s=rT\) and \(D=1-2s\cos\theta+s^2\). The exact
finite geometric sum gives
\[
 \begin{split}
         P_m(z)&=A(1-ae^{-i\theta})
                         +e^{im\theta}B(1-be^{i\theta}),\\
       A&=\mathbb E_m[D^{-1}],&
       B&=\mathbb E_m[s^{m+2}/D],\\
       a&=\frac{\mathbb E_m[s/D]}{\mathbb E_m[D^{-1}]},&
       b&=\frac{\mathbb E_m[s^{m+1}/D]}{\mathbb E_m[s^{m+2}/D]} .
 \end{split}                                                     \tag{3}
\]
All integrals are finite. At infinity \(D\sim s^2\), and
\(D\) is bounded below by a positive number for fixed
\(0<\theta<\pi\).

The successive moment ratios of the measure
\(D^{-1}d\mathbb P_m\) strictly increase whenever
\(T\) is nonconstant, by Cauchy–Schwarz. In particular
\(ab<1\). Put
\[
                       p_\theta(u)=\arg(1-u e^{-i\theta}),
                                         \qquad u>0 .
\]
This function strictly increases from zero to
\(\pi-\theta\), and
\(p_\theta(a)+p_\theta(1/a)=\pi-\theta\).
For the last identity, the product of the two complex
factors is a negative real scalar times \(e^{-i\theta}\).
The sum of the two arguments lies between zero and
\(2(\pi-\theta)\), which fixes the branch. Therefore
\[
                         0<p_\theta(a)+p_\theta(b)<\pi-\theta.
\]
At a root of (3) we have
\[
                m\theta=(2k-1)\pi+p_\theta(a)+p_\theta(b)
\]
for some integer \(k\). The bounds on the phase give
exactly the open boxes (1).

To count the roots, replace \(T\) in (2) by \(T^t\) for
\(0\le t\le1\), keeping the fixed measure \(\mathbb P_m\).
The moment coefficients are continuous, by domination
by \(1+T^m\), and the leading coefficient stays positive.
At \(t=0\) the geometric sum has simple roots, at the
right endpoints of the indicated boxes and their
reflections, plus \(-1\) in odd degree. For small
positive \(t\), the strict sector exclusion just proved
puts one root in each corresponding open box. For
\(t>0\) each allowed upper or lower box is separated
from the others and from the real axis. Roots therefore
keep their counts with multiplicity throughout the
homotopy. The geometric sum's remaining odd root
stays real by conjugation and cannot be positive by
coefficient positivity. This proves simplicity,
including that of the odd real root.

## Consecutive labels and the remaining comparison

In the geometric exception the nontrivial roots of unity
at consecutive orders have coprime denominators \(n+1,n+2\),
so those roots are distinct. For a nonconstant ratio, the lower \(k\)-th box ends before the upper
\((k+1)\)-st box starts since
\(2k\pi/(n+1)<(2k+1)\pi/(n+1)\).
The upper \(k\)-th box ends before the lower
\((k+1)\)-st box starts since
\(2k\pi/(n+2)<(2k+1)\pi/n\). A common argument
thus has the same label \(k\). Same-label overlap
requires
\[
                         \frac{(2k-1)\pi}{n}
                               <\frac{2k\pi}{n+2},
\]
which is equivalent to \(4k<n+2\).
As \(4k\le n+1\), the lower box then ends at or
before \(\pi/2\), and its angles are strictly
less than \(\pi/2\). The negative real root occurs
only in odd degree, so it cannot be common to two
consecutive moments.

There is a one-function form of the remaining
adjacent comparison. On the upper half-plane define
\[
       F_m(z)=-\frac{\mathbb E_m[(1-zT)^{-1}]}
                       {\mathbb E_m[T^{m+1}(1-zT)^{-1}]} .
\]
The denominator has strictly positive imaginary part
and is finite. In terms of (3), this is
\(z\, K_+(z)/K_-(z)\), where
\(K_+(z)=\mathbb E_m[(1-zT)^{-1}]\) and
\(K_-(z)=\mathbb E_m[T^m/(1-(zT)^{-1})]\).
Consequently \(\arg F_m(z)=\theta+p_\theta(a)+p_\theta(b)\)
lies strictly between \(\theta\) and \(\pi\) for
nonconstant \(T\). The root equation is
\[
                              z^{m+1}=-F_m(z).
\]
In particular a consecutive common root must solve
\(F_{n+1}(z)=zF_n(z)\) in the indicated right-half
angular boxes. In the mirrored case the
[matched folded resolvent comparison](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
excludes that equality at every root. Arbitrary unmatched
rows change both endpoint weights in \(F_m\), and their
remaining same-label comparison is open here.
