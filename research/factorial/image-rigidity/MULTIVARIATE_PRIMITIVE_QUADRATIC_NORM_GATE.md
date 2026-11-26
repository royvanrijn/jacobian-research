# A primitive norm gate for common quadratic factors

Let \(F,G\) be nonzero rational polynomials. Strip powers of \(z\)
and replace each by its primitive integral squarefree part. Suppose
a common irreducible factor \(H(z)=z^2-uz+c\) has an upper-half-plane
root \(\alpha\) with \(0<\arg\alpha<\Theta<\pi/2\). Thus \(u,c>0\).
At each prime, the Newton polygons of \(F,G\) bound the possible
unordered pairs of valuations of \(\alpha,\bar\alpha\).
The [field discriminant bounds](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md)
give further constraints. The following restriction combines the
field, absolute valuations, and angle cap in an integral norm
equation with a class-group test. It applies to ordinary factorial
moment polynomials with any support, once a common quadratic in
such an angular sector is under consideration.

## Primitive scale and the norm equation

Write \(K=\mathbb Q(\alpha)\), with fundamental field discriminant
\(D_K<0\). For each rational prime \(p\) write \(s_p,t_p\) for the
valuations of \(\alpha,\bar\alpha\) in an algebraic closure of
\(\mathbb Q_p\), normalized by \(v_p(p)=1\). At a split prime both
are integers. At an inert prime they are equal integers. At a
ramified prime they are equal elements of \(\frac12\mathbb Z\).
These facts follow from the respective local ramification indices
and uniqueness of the extension of \(v_p\) in the nonsplit case.

Define the positive rational number
\[
 q=\prod_p p^{-\lfloor\min(s_p,t_p)\rfloor}.
\]
These powers are signed. Only finitely many are nonzero, since at
every prime where either primitive integer polynomial has unit
endpoints all its roots are units. The element
\(\beta=q\alpha\) is integral and no rational prime divides it in
\(\mathcal O_K\). Indeed its valuations are all nonnegative and at
each prime at least one of the two valuations is less than one.
The same condition forces every exponent of \(q\), so \(\beta\)
is the unique integral element with this property in the positive
rational scaling class of \(\alpha\). Replacing \(\alpha\) by its
conjugate replaces \(\beta\) by its conjugate.

Set
\[
 \begin{split}
 D&=\prod_{s_p\ne t_p}p^{e_p},\qquad e_p=|s_p-t_p|,\\
 R&=\prod_{\substack{p\ \mathrm{ramified\ in}\ K\\
                         s_p=t_p\in\frac12+\mathbb Z}}p .
 \end{split}
\]
Here \(D\) is supported on split primes, \(R\) is squarefree and
divides \(|D_K|\), and \(\gcd(D,R)=1\).
The valuations of \(\beta\) at a split prime dividing \(D\) are
\(0,e_p\), and at a prime dividing \(R\) they are both \(1/2\).
All other valuations of \(\beta\) vanish. Therefore
\(N_{K/\mathbb Q}(\beta)=DR\).

Write \(\beta=(x+y\sqrt{D_K})/2\), using the square root with
positive imaginary part. The basis
\(1,(D_K+\sqrt{D_K})/2\) of \(\mathcal O_K\) gives the exact
conditions
\[
 x,y\in\mathbb Z_{>0},\qquad x\equiv D_Ky\pmod 2,\qquad
 \gcd\!\left((x-D_Ky)/2,y\right)=1 .
\]
The basis follows by writing \(K=\mathbb Q(\sqrt{k})\) for negative
squarefree \(k\). For an integral element \(a+b\sqrt{k}\) the trace
gives \(2a\in\mathbb Z\), and \(k(2b)^2\in\mathbb Z\) follows from
its norm. Squarefreeness gives \(2b\in\mathbb Z\). The norm condition
modulo \(4\) gives the basis \(1,(1+\sqrt{k})/2\) if
\(k\equiv1\pmod4\) and \(1,\sqrt{k}\) otherwise. Their discriminants
are \(k\) and \(4k\), respectively, and changing the second basis
element by an integer gives the displayed basis.

In these coordinates the norm and the scale-invariant angular defect
\(\delta=4-u^2/c\) satisfy
\[
  4DR=x^2+|D_K|\,y^2,\qquad
  D\delta=(|D_K|/R)y^2,\qquad
  u=x/q,\qquad c=DR/q^2 .                                      \tag{1}
\]
At a prime with unequal root valuations the trace has valuation
\(\min(s_p,t_p)\), so \(v_p(\delta)=-e_p\). At an equal-valuation
prime the trace has valuation at least \(s_p=t_p\), so \(\delta\)
is \(p\)-integral. Thus \(D\) is the reduced denominator of
\(\delta\) and \(\gcd(y,D)=1\).

Since \(0<\delta=4\sin^2(\arg\alpha)<4\sin^2\Theta\), (1) gives
\[
    D>\frac{|D_K|}{4R\sin^2\Theta},\qquad
    0<y^2<\frac{4DR\sin^2\Theta}{|D_K|}.                         \tag{2}
\]
One may use any explicit upper bound for \(\sin^2\Theta\).
For example, \(\Theta=2\pi/(n+2)\) with \(n\ge3\) gives
\(D>(|D_K|/R)(n+2)^2/160\).
For a given field and slope package, test whether
\(4DR-|D_K|y^2\) is a positive square for the remaining \(y\).
The parity and coprimality checks enforce integrality and
primitivity. Conversely, an integral primitive solution of the
norm equation puts the entire exponent of each \(p\mid D\)
on one of its two prime ideals, since positive exponents on
both would make \(p\) divide \(\beta\). At \(p\mid R\) the
exponent in the norm is one. Therefore every retained solution
has the prescribed unordered slopes after rescaling by \(q^{-1}\).

## The exact discriminant index

Write \(q=f/g\) with positive coprime integers \(f,g\), and put
\(M=DR\). The polynomial
\[
                  f^2 z^2-fg xz+g^2 M
\]
has content \(h=\gcd(f,M)\). A prime away from \(f\) does not
divide the leading coefficient. At a prime \(p\mid f\) with
\(p\nmid M\) the constant coefficient is a unit. If \(p\mid f\) and \(p\mid R\),
the constant coefficient has valuation one and all coefficients
are divisible by \(p\). If \(p\mid f\) and \(p\mid D\), the two valuations
of \(\beta\) are unequal and their minimum is zero, so
\(v_p(x)=0\). The three coefficient valuations then have
minimum \(\min(v_p(f),v_p(M))\). This proves the content formula,
including signed rescalings where \(g>1\).

It follows that the primitive integral factor \(H_0\) and its
discriminant are
\[
 H_0(z)=\frac{f^2 z^2-fg xz+g^2M}{h},\qquad
       \operatorname{disc} H_0
                 =D_K \left(\frac{fg}{h}y\right)^2 .            \tag{3}
\]
Gauss's lemma and the discriminant product formula imply
\[
 |D_K|\left(\frac{fg}{h}y\right)^2
       \ \bigm|\ \gcd(|\operatorname{disc} F|,
                     |\operatorname{disc} G|).
\]
The chosen physical scale therefore constrains the candidate
\(y\) in addition to the norm and angle bounds. A surviving
candidate can be checked for its angle and by exact division
of \(F,G\) by \(H_0\).

## The ideal-class check

For each \(p\mid D\) let \(\mathfrak p_p\) be the prime above \(p\)
on which the nonzero valuation of \(\beta\) lies. For each
\(p\mid R\) let \(\mathfrak r_p\) be the unique prime above \(p\).
Factorization of ideals gives the necessary condition
\[
 (\beta)=\prod_{p\mid D}\mathfrak p_p^{e_p}
          \prod_{p\mid R}\mathfrak r_p,\qquad
 \prod_{p\mid D}[\mathfrak p_p]^{e_p}
          \prod_{p\mid R}[\mathfrak r_p]=1
                    \quad\hbox{in }\mathrm{Cl}(K).              \tag{4}
\]
One may test (4) for either choice of prime above each split
prime. The other choice inverts its class. Since
\(\mathfrak r_p^2=(p)\), squaring (4) eliminates the ramified
classes. In particular, if \(D=p^{e_p}\) and the split prime
class has order \(h_p\), then \(h_p\mid 2e_p\). If \(h_p\) is odd,
then \(h_p\mid e_p\).

Conversely, when the ideal on the right of (4) is principal,
any generator has the specified prime valuations and norm
\(DR\). Its associates differ by units. In an imaginary quadratic
field the units have norm one and there are finitely many, as
the unit circle meets the integral lattice in a finite set.
The norm test enumerates exactly the generators with positive
real and imaginary parts as the split choices vary, and the
first angle box and exact division impose the remaining tests.

For an exact illustration, take \(D_K=-23\) and sole split spread
\(e_2=7\). The integer \(\omega=(1+\sqrt{-23})/2\) satisfies
\(\omega^2-\omega+6=0\), so \(2\) splits with \(\omega\) reducing
to \(0\) and \(1\) at its two primes. Neither prime is principal,
since \(x^2+23y^2=8\) has no integral solution. The element
\(1+\omega=(3+\sqrt{-23})/2\) has norm \(8\) and reduces to zero
only at the prime with \(\omega=1\). Its principal ideal is
therefore that prime cubed, and each split class has order
three. This package is impossible for either value of \(R\),
since \(3\nmid 14\). The size condition by itself permits
\(D=128,R=1\) at \(\Theta=1/4\), since
\(512\sin^2(1/4)>512(95/384)^2>23\).
The example specifies a field and slopes without requiring a
factorial profile that realizes them.
