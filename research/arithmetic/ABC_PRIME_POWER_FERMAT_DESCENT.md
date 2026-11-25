# A support-uniform primitive count for weighted prime-power equations

Let \(\ell\) be a prime and \(d=\ell^j\), \(j\ge1\), with
\((\ell-1)(d-2)/2\ge2\). For arbitrary nonzero integers \(A,B,C\)
one has
\[
 \#\{(x,y,z)\in\mathbb Z^3:
             \gcd(x,y)=1,\ Ax^d+By^d=Cz^\ell\}
 \le \exp\!\big(O_{\ell,d}(1+\omega(|ABC|))\big).
 \tag{FP1}
\]
There is no height bound on the solutions or the coefficient
valuations. For polynomial-height coefficients this gives a
subpower bound for each fixed coefficient record, uniformly in
the coefficients. It applies when \(d\mid d_1,d_2\) and
\(\ell\mid d_3\) by replacing the roots in
\(Ax^{d_1}+By^{d_2}=Cz^{d_3}\) with their appropriate powers.
Signs may be carried by the coefficients.

The coefficient-uniform point estimates of Browning–Verzobio
apply to general signatures and give positive powers of box
sizes. The support bound above uses the special prime-power
splitting field and the primitivity of the two weighted
power roots. On this signature it replaces the height factor
by a function only of the coefficient support.
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1)

## Uniform rational-point bound from the splitting field

Let \(F\in\mathbb Z[T]\) be squarefree of degree \(d\) divisible
by \(\ell\), with leading coefficient \(a_d\). Let \(K\) contain
its roots and \(\zeta_\ell\), with \([K:\mathbb Q]\le D\).
For a nonzero integer \(n\), let \(\mathcal C\) be the smooth
projective model of \(n V^\ell=F(U)\), with Jacobian \(J\).
Take \(S\) to be the finite places of \(K\) above the primes
dividing \(\ell n a_d\operatorname{disc}(F)\).
When \(g=(\ell-1)(d-2)/2\ge2\),
\[
 \#\mathcal C(K)\le
 \exp O_{d,\ell,D}\!\left(
                1+|S|+\log|\mathrm{Cl}(K)[\ell]|\right).
 \tag{FP2}
\]
The root divisors yield rational torsion for the cyclic-cover
endomorphism. The rank cost therefore comes from the
\(\ell\)-part of the class group of the splitting field.

Indeed, the cover of the \(U\)-line ramifies just at the \(d\)
distinct roots, each with index \(\ell\), so Riemann–Hurwitz
gives the stated genus. Write \(\sigma\) for the endomorphism
of \(J\) from the deck transformation
\((U,V)\mapsto(U,\zeta_\ell V)\), and \(\phi=1-\sigma\).
The norm \(1+\sigma+\cdots+\sigma^{\ell-1}\) factors through
the Jacobian of the line and vanishes. The cyclotomic
identity in \(\mathbb Z[\zeta_\ell]\) gives
\([\ell]=u\phi^{\ell-1}\) for a unit endomorphism \(u\).
Thus \(\phi\) is an isogeny of degree
\(\ell^{2g/(\ell-1)}=\ell^{d-2}\).

Write \(P_i\) for the point over a root \(\alpha_i\).
The classes of degree-zero divisors on the \(P_i\) lie in
\(J[\phi]\), and are killed by \(\ell\) because
\(\operatorname{div}((U-\alpha_i)/(U-\alpha_1))
=\ell(P_i-P_1)\). If \(\sum a_iP_i\) with \(\sum a_i=0\)
is principal, the function defining it is a deck eigenfunction:
its transform divided by itself has trivial divisor. It has
the form \(V^h R(U)\) in the cyclic function field, so its
orders at the branch points all equal \(h\bmod\ell\).
Conversely, for such a coefficient vector the function
\[
                V^h\prod_i(U-\alpha_i)^{(a_i-h)/\ell}
\]
has these branch orders and order zero at all points over
infinity. The latter follows from the order \(-d/\ell\)
of \(V\) and from \(\sum a_i=0\).
The kernel on the space of branch distributions modulo
\(\ell\) is therefore the span of the all-ones vector.
Their classes have order \(\ell^{d-2}\), exhausting
\(J[\phi]\) over \(K\).
[Poonen–Schaefer, Proposition 6.2](https://math.mit.edu/~poonen/papers/descent.pdf)

The Kummer map injects \(J(K)/\phi J(K)\) into
\(H^1(K,J[\phi])\cong(K^\times/K^{\times\ell})^{d-2}\).
At a place outside \(S\), the affine equation has separate
simple branch points at unit leading coefficient on the
residue field. At infinity set \(U'=1/U\) and
\(V'=V(U')^{d/\ell}\). The resulting two monic
finite covers of the line charts glue to a smooth proper
curve over the local valuation ring. Smoothness follows
from the \(V\)-derivative away from the branch points and
the \(U\)-derivative at them. Its Jacobian has good
reduction, and the action and isogeny extend on the
abelian scheme. Since \(\phi\) divides \([\ell]\) and
\(\ell\) is invertible, \(\phi\) is finite étale there.
A \(K\)-point of the Jacobian extends to the valuation
ring. Its inverse image under \(\phi\) is finite étale,
so the local Kummer class is unramified.

In the power-class group with exponent \(\ell\),
the unramified classes have valuations divisible by
\(\ell\) outside \(S\). Their group has size at most
\(\ell^{O_D(1)+|S|}|\mathrm{Cl}(K)[\ell]|\): take the
ideal over the ring of \(S\)-integers, divide its valuations
by \(\ell\), then use the \(S\)-unit rank for the kernel
and the class group for its image. The class group of
\(S\)-integers is a quotient of \(\mathrm{Cl}(K)\) and its
\(\ell\)-torsion has at most the same size.
This bounds the logarithm of \(\#J(K)/\phi J(K)\).
The filtration by successive images under \(\phi\) yields
\[
 \ell^{\operatorname{rank}J(K)}
     \le \#J(K)/\ell J(K)
     \le \big(\#J(K)/\phi J(K)\big)^{\ell-1}.
\]
Embed \(\mathcal C\) into \(J\) using any \(P_i\in\mathcal C(K)\).
An Abel–Jacobi curve of genus at least two contains no
positive-dimensional abelian translate, and its degree in a
multiple of the canonical polarization depends only on the
genus. Uniform Mordell–Lang gives
\(O_g(1)^{1+\operatorname{rank}J(K)}\) points, proving (FP2).
[Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)

## Prime-power splitting fields

For \(F(U)=AU^d+B\) with \(d=\ell^j\), take
\[
                K=\mathbb Q(\zeta_d,(-B/A)^{1/d}).
\]
It is Galois of \(\ell\)-power degree over the fixed
\(K_0=\mathbb Q(\zeta_\ell)\), since both the cyclotomic
step and the Kummer step have \(\ell\)-power degree.
At all finite primes away from \(\ell AB\) the polynomial
has separable reduction with unit leading coefficient, so
its roots are acquired by an unramified extension. Factor
\(K/K_0\) into a bounded number of cyclic degree-\(\ell\)
steps. Each has \(O_d(1+\omega(|AB|))\) ramified primes.
The cyclic invariant-class bound for each step gives
\[
           |\mathrm{Cl}(K)[\ell]|
                      \le \exp O_{\ell,d}(1+\omega(|AB|)).
\]
In detail, for a step \(L/E\), the invariant ideal classes
are generated by extensions of classes of \(E\) and the
ramified ideals, up to an index at most the unit norm
index. This follows by applying Hilbert 90 to the multiplier
of any invariant class with trivial norm obstruction.
The \(\ell\)-rank of \(\mathrm{Cl}(L)[\ell]\) is at most
\(\ell\) times the invariant rank by nilpotence of
\(\sigma-1\). Thus the estimate is uniform in the
coefficient valuations and does not require irreducibility
of the binomial.
[Cyclic-step class group argument](ABC_PRIMITIVE_NORM_RAY_PACKING.md#prime-power-towers-remove-the-class-number-cost)

The discriminant of \(AU^d+B\) has prime factors only
among those dividing \(dAB\). With \(n=C\), (FP2) yields
the right side of (FP1). If \(y\ne0\), take
\(U=x/y,\ V=z/y^{d/\ell}\). Primitivity of the pair
determines \((x,y)\) from \(U\) up to common sign, and
for each pair the equation determines at most \(\ell\)
integers \(z\). If \(y=0\), then \(x=\pm1\), so there are
only boundedly many solutions. This proves (FP1).

In particular, equal exponents \(k\ge5\) other than \(6,12\)
have the support bound with a constant depending on \(k\).
Such \(k\) has a prime divisor at least five, or is divisible
by \(8\) or \(9\), and the same prime-power divisor can be
used on the two roots with its prime used on the third.
For polynomial-height coefficients and complete values of
height at most \(X\), the subpower fixed-coefficient bound
can be uniform in these \(k\). For \(k\ge 4/\varepsilon\),
enumerate the two power roots in \(O(X^{2/k})\) and
determine the third. All smaller admitted \(k\) form a
fixed finite list.

## A three-cap coefficient-defect inequality

In a primitive \(abc\) triple, choose any two letters as
\(a_1,a_2\) and write \(a_1=A x^d,\ a_2=B y^d\)
using their exact \(d\)-free coefficients. On the remaining
letter \(a_3=C z^\ell\), use its exact \(\ell\)-free
coefficient. After arranging the signed equation,
\(\gcd(x,y)=1\). If
\[
       \operatorname{rad}(AB)\le X^{f\lambda},\qquad
       \operatorname{rad}(C)\le X^{h\lambda},
\]
then the number of such triples up to height \(X\) is
\[
             O_{\ell,d,\lambda,\varepsilon}
                                (X^{(f+h)\lambda+\varepsilon}).
 \tag{FP3}
\]
For each of the two squarefree support labels there are
at most \(O_{\ell,d}(1)^{\omega(\operatorname{rad}(ABC))}\)
coefficient choices, since the residue exponents are
bounded and one colors each prime of \(AB\) by its letter.
Equation (FP1) contributes another factor of the same
form. The support caps and the elementary bound on
\(\omega(n)\) give (FP3).

Equivalently let \(s_i\) be the cap for the original
radical of \(a_i\) in units of \(\lambda\), and \(v_i(m)\)
the logarithmic radical mass at actual valuations divisible
by \(m\), in the same units. In a box fixing those masses,
(FP3) has coefficient
\[
                 \sum_i s_i-v_1(d)-v_2(d)-v_3(\ell).
 \tag{FP4}
\]
All three valuations enter. In the endpoint normalization
with total radical mass at most one, a fixed margin over
\(2/3\) in the sum of recovered divisibility masses gives
an absolute saving in that sector. Exact integer divisibility
is essential. For example spread each letter's mass over
the prime degrees in \([J,2J]\), putting equal mass \(1/3\)
in each complete letter and using \(\lambda=3/\bar j\)
with \(\bar j\) the mean prime degree. The limiting scaled
degree distributions are diffuse, with height weight one,
while any one choice of \(\ell,d\) recovers mass from
at most one degree on each letter. This is a permissible
profile datum and does not construct any triples.

On the balanced minimal-repair sector with
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\), a
repaired height \(H_M(c)\le Y\) implies
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\). Thus (FP3) gives
the upper exponent
\[
                   \frac{f+h}{M(1/3-\xi)}+\varepsilon.
\]
It is strictly below \(1/M\) when \(f+h<1/3-\xi\).
The reduction runs on the primitive quotient. No common
multiplier is included in (FP1) or in the support count.
[The minimal-target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
