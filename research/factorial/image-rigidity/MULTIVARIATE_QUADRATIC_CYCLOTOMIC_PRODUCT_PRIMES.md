# Fresh primes and paired cyclotomic products for a common quadratic

An irreducible common quadratic with noncyclotomic root quotient
forces common residue pairs at primes of arbitrarily large size.
The order of the quotient modulo each such prime can be prescribed
as any prime \(q>30\). Two integral products constructed from the
original parent polynomials make this a finite test for each \(q\).
Their exact divisibility links the angle numerator, reduced trace
denominator, and primitive endpoint product. It does not assume
a bound on the number, sizes or orientations of any exponent rows.

## Integral pair products

For an integer \(q\ge2\) let
\[
       \Psi_q(x,y)=\sum_{k=0}^{q-1}x^{q-1-k}y^k .
\]
For an integral polynomial \(F\) of degree \(m\ge1\) with
leading coefficient \(a\) and roots \(\alpha_i\), and an
integral polynomial \(G\) of degree \(n\ge1\) with leading
coefficient \(b\) and roots \(\beta_j\), use
\[
 \begin{aligned}
  I_q(F)&=a^{(q-1)(m-1)}
                          \prod_{i<j}\Psi_q(\alpha_i,\alpha_j),\\
  J_q(F,G)&=a^{(q-1)n}b^{(q-1)m}
                          \prod_{i,j}\Psi_q(\alpha_i,\beta_j).
 \end{aligned}                                                    \tag{1}
\]
Both quantities are integers. In the first line the product
is a symmetric polynomial in the \(\alpha_i\), with integer
coefficients and degree at most \((q-1)(m-1)\) in each variable.
A symmetric polynomial with per-variable degree at most \(B\)
is an integer combination of products of at most \(B\)
elementary symmetric functions. To see this, successively
subtract a multiple of \(\prod_j e_j^{\lambda_j-\lambda_{j+1}}\)
to remove a leading monomial with decreasing exponents
\(\lambda_1,\ldots,\lambda_m\) in lexicographic order.
This product has the same leading monomial with coefficient
one and uses exactly \(\lambda_1\le B\) factors, so the
degree bound persists. The elementary symmetric functions
of the \(\alpha_i\) have denominators dividing \(a\).
Applying the same argument separately in the two sets of
roots proves integrality of the second line of (1).

These products can be evaluated without solving for roots.
Define \(F^{[q]}(z)=a^q\prod_i(z-\alpha_i^q)\), an integral
polynomial, and similarly for \(G\). The same symmetric
degree argument proves integrality. The standard
discriminant and resultant formulas give
\[
 \operatorname{disc}(F^{[q]})=\operatorname{disc}(F)I_q(F)^2,
 \qquad
 \operatorname{Res}(F^{[q]},G^{[q]})
             =\operatorname{Res}(F,G)J_q(F,G).
\]
These are identities including repeated roots and zero
resultants. Alternatively, \(J_q(F,G)\) is
\(\operatorname{Res}_z(F(z),\operatorname{Res}_w(G(w),\Psi_q(z,w)))\).
This last expression does not divide by a resultant. The
usual resultant convention with the polynomial having roots
listed first gives the displayed identity.

## Exact coupling of the three integer channels

Let \(H=A z^2+Bz+C\) be a primitive irreducible common factor
of the ordinary parents, with \(A,C>0\). Its roots
\(\alpha,\overline\alpha\) are nonreal by positivity of the
even parent. For now assume \(B\ne0\), and set
\[
         g=\gcd(B^2,AC),\qquad D=AC/g,\qquad U=B^2/g,
                      \qquad V=4D-U>0 .
\]
Then \(T=U/D-2\) is the reciprocal trace
\(\gamma+\gamma^{-1}\) of \(\gamma=\alpha/\overline\alpha\).
For any prime \(q>3\) set \(d=(q-1)/2\) and use the monic
integer polynomial \(C_q(t)\) characterized by
\[
         C_q(Y+Y^{-1})=Y^{-d}\Phi_q(Y).
\]
Its integrality follows by expanding the right side as a
sum of \(Y^\ell+Y^{-\ell}\) and one. Define the nonzero
integer \(L_q=D^d C_q(T)\). It is nonzero for a rational
trace \(T\) when \(q>3\), since the degree of the real
cyclotomic trace is \(d\). One has exactly
\[
 \begin{aligned}
  g^d L_q
     &\mid \gcd\bigl(I_q(F),I_q(G)\bigr),\\
  q^2 (AC)^{q-1} g^{q-1} L_q^2
     &\mid J_q(F,G).                                            \tag{2}
 \end{aligned}
\]
Zero in an invariant is allowed in this divisibility statement.

Indeed \(\Psi_q(\alpha,\overline\alpha)=(\alpha\overline\alpha)^d
 C_q(T)\), so \(I_q(H)=(AC)^d C_q(T)=g^d L_q\).
By Gauss's lemma the two parents have integral quotients by
\(H\). Formula (1) is multiplicative in each argument of
\(J_q\), and \(I_q(HF_1)=I_q(H)I_q(F_1)J_q(H,F_1)\)
when \(F_1\) has positive degree. For a constant quotient
\(F_1=c\), the self product is \(c^{q-1}I_q(H)\) directly.
The cross product with a constant is determined by its
leading coefficient power in (1). Finally the two diagonal
terms in \(J_q(H,H)\) are \(q\alpha^{q-1}\) and
\(q\overline\alpha^{\,q-1}\), giving
\(J_q(H,H)=q^2(AC)^{q-1} I_q(H)^2\).
This proves (2).

One can use all the prime valuations of (2) simultaneously.
At a prime \(p\), with \(v_p(0)=+\infty\), the candidate
must obey
\[
 \begin{split}
 v_p(g^d L_q)\le \min\biggl\{&
 v_p(I_q(F)),\,v_p(I_q(G)),\\
 &\frac{v_p(J_q(F,G))-2v_p(q)-(q-1)v_p(AC)}2\biggr\}.           \tag{3}
 \end{split}
\]
In particular the same primitive endpoint valuations
determine both the available norm power and the angular
divisor in the mixed product. These are the endpoint
valuations of the proposed quadratic in the original
variable, so they can use the
[joint quadratic allocation](MULTIVARIATE_QUADRATIC_JOINT_PARITY_DISCRIMINANT_BUDGET.md).
At any distinct primes \(q,r>3\), the integers \(L_q,L_r\)
are coprime. A rational prime \(p\mid D\) divides neither,
since \(C_q\) and \(C_r\) are monic and the rational trace
has reduced denominator \(D\). Otherwise solve
\(Y+Y^{-1}=T\) in the algebraic closure of \(\mathbb F_p\).
At \(p\ne q,r\), a shared zero would have simultaneously
orders \(q\) and \(r\), which is impossible. At \(p=q\)
a zero of \(\Phi_q\) must be \(Y=1\), where
\(\Phi_r(1)=r\ne0\); the case \(p=r\) is symmetric.

## A mandatory fresh prime of prescribed quotient order

Suppose \(\gamma\) is not a root of unity.
Take \(\lambda=A\alpha/\sqrt{g}\) and
\(\lambda'=A\overline\alpha/\sqrt{g}\) with the same
positive square root of \(g\). Their product is \(D\)
and their squared sum is \(U\). Both are algebraic
integers, since they satisfy the monic polynomial
\(X^4+(2D-U)X^2+D^2\). The integers \(U,D\) are nonzero
and coprime. Thus \((\lambda,\lambda')\) is a Lehmer
pair, and \(L_q=(\lambda^q-\lambda'^q)/(\lambda-\lambda')\)
for odd prime \(q\). The quantity \((\lambda^2-\lambda'^2)^2\) is \(-UV\).
[Bilu, Hanrot and Voutier, Theorem 1.4](https://doi.org/10.1515/crll.2001.080)
therefore gives a primitive divisor \(p\mid L_q\) for
every prime \(q>30\), with \(p\nmid UV\).
Also \(p\nmid D\) by the denominator argument above.
The units \(\lambda,\lambda'\) remain distinct at a place
above \(p\) and their quotient reduces to an element
of order exactly \(q\). Indeed \(p\nmid UV\) makes their
square-difference a unit and \(p\nmid D\) makes each a
unit. The \(q\)-th powers agree, and the ratio is not one.
In particular \(p\ne q\).

The discriminant squareclass determines whether this
ratio reduces in \(\mathbb F_p\) or in its norm-one
quadratic extension: its trace \(T\) reduces rationally,
and its quadratic discriminant has class \(-UV\). Hence
\[
           q\mid p-\chi_K(p),\qquad p\ge q-1,                 \tag{4}
\]
where \(\chi_K(p)\) is \(1\) or \(-1\) according as \(p\)
splits or is inert in \(K\). Such a prime divides both
\(I_q(F),I_q(G)\) and its square divides \(J_q(F,G)\)
by (2). No prime can serve two different odd primes \(q\),
because their \(L_q\)'s are coprime.

For an original ordinary profile, take \(M\) to bound all
factorial arguments and both orders. Every prime of any
individual coefficient is at most \(M\), even after
passing to primitive integer parents. Choosing an odd
prime \(q>\max(30,M+1)\) gives a required prime outside
this entire coefficient support. Both parent reductions
then have full degree and unit endpoints. The two roots
of the common quadratic reduce to distinct nonzero roots
of both parents, and their quotient has order \(q\).
In particular, disappearance of the required large primes
from the common divisor of the self products excludes
every noncyclotomic common quadratic. The angular
discriminant squareclass, when known, can restrict the
required congruence in (4) further.

The test may be used at any degree. A sufficient way to
make the self products nonzero before computing them is to
take \(q-1>m(m-1),n(n-1)\). An exact zero would give a
primitive \(q\)-th root of unity as some root ratio and
its trace of degree \((q-1)/2\) as a root of the rational
pairwise-trace polynomial of degree \(m(m-1)/2\) or
\(n(n-1)/2\). If the quotient is cyclotomic, only the orders
\(2,3,4,6\) can occur for a nonreal conjugate pair in an
imaginary quadratic field. Those possibilities can be
tested separately with the corresponding products.
A uniform argument that at least one required fresh
prime fails to occur for every remaining ordinary pair
is still needed.

## The sign forces an inert-prime parity

For a prime \(q>30\) put \(\delta_q=\mathbf 1_{q\mid V}\)
and \(L_q^\circ=L_q/q^{\delta_q}\).
Every prime \(p\mid L_q\) other than \(q\) is primitive
in the Lehmer sequence. To see this, \(p\nmid D\) by the
monic trace polynomial argument. Reducing
\(Y+Y^{-1}=T\) modulo \(p\) for \(Y\) the root quotient,
one obtains \(\Phi_q(Y)=0\). As \(p\ne q\), the order
of \(Y\) is exactly \(q\). In particular \(Y^2\ne1\),
so \(p\nmid UV\), and no earlier Lehmer number vanishes
modulo \(p\).

In characteristic \(q\) one has
\(C_q(t)\equiv(t-2)^d\) because
\(\Phi_q(Y)=(Y-1)^{q-1}\). Also \(C_q(2)=q\).
Thus
\[
             L_q\equiv(-V)^d\pmod q .
\]
If \(q\mid V\), then \(D\) is a unit at \(q\) and
\(T-2=-V/D\) is divisible by \(q\).
The Taylor expansion
\(C_q(2+X)=X^d+qR(X)\) with \(R(0)=1\)
shows \(L_q/q\equiv D^d\pmod q\) since \(d\ge2\).
In particular \(v_q(L_q)=\delta_q\).

Combining (4) over the prime factors of \(L_q^\circ\)
now yields the exact sign restriction
\[
 (-1)^{\sum_{p\ \mathrm{inert},\ p\mid L_q^\circ}v_p(L_q^\circ)}
    =\operatorname{sgn}(L_q)
       \begin{cases}
            \left(\frac{-V}{q}\right),&q\nmid V,\\[2mm]
            \left(\frac{D}{q}\right),&q\mid V .
       \end{cases}                                              \tag{5}
\]
Indeed each prime in that product is \(1\) or \(-1\)
modulo \(q\) according to its splitting symbol, and
Euler's criterion evaluates the preceding residues.
The real sign in (5) is the sign of \(\sin(q\theta)\)
for \(\theta=\arg\alpha\in(0,\pi)\), by
\(L_q=D^d\sin(q\theta)/\sin\theta\).
Thus a negative right side forces an inert prime
\(p\equiv-1\bmod q\) in the common divisor of the
two self products. A strict angle interval can already
fix the real sign. This refinement relates it to the
discriminant complement modulo \(q\), or to the reduced
norm denominator in the index-prime case.

## An ordinary pair permitting a quadratic at the first residue

Take the mirrored rows \((13,12),(12,13)\) at orders
\(3,4\). The primitive parent arrays in increasing order
are
\[
                     (13,37,37,13),\qquad
                     (221,833,1225,833,221).
\]
Both have unit endpoints at \(p=3\); their common residual
polynomial there is \(1+z^2\). This first test permits a
quadratic. The gcd of the self products is one at both
\(q=31\) and \(q=59\), while the products for \(q=2,3,4,6\)
are nonzero. The prime bound \(q-1=58\) at \(q=59\)
already exceeds the largest factorial argument \(52\).
Thus no noncyclotomic or cyclotomic common quadratic is
possible. Since the lower degree is three, even positivity
then excludes a common factor of any positive degree.
This is an exact comparison of tests on one original
profile. It does not require a finite enumeration of
terminal orders.
