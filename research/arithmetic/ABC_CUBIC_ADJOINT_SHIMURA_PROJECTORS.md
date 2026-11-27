# The cubic adjoint comparison with Shimura projector denominators

Let \(a+b=c\) be primitive and positive, put \(R=\operatorname{rad}(abc)\),
\(A=\log R\), \(h=\log c\), and use \(E:y^2=x(x-a)(x+b)\). Write
\(N=N_E\) for the conductor, so \(N/R\) is a bounded power of \(2\),
and put \(S=\{2,3,7\}\). Fix either cubic character \(\chi\) modulo
\(7\), set \(K=\mathbb Q(\zeta_3)\), and use the Tate twist \(s=1\)
throughout the local calculations. Let \(X_E\), \(C\), \(\mathcal D\),
and \(Z_E=\mathrm N_{K/\mathbb Q}(X_E)/(C\mathcal D)\) have their
earlier definitions. The value \(X_E\) and hence \(Z_E\) are used
for triples other than \((1,1,2)\).
[Critical value and Néron period](ABC_FIXED_CUBIC_ADJOINT_VALUE.md).
[Off-diagonal factors](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).
[Diagonal factors](ABC_CUBIC_ADJOINT_DIAGONAL_FACTORS.md).

## The rational modular degree and the adjoint quotient

For a modular parametrization \(\phi:X_0(N)\to E\) over \(\mathbb Q\),
write \(\phi^*\omega_E=m_\phi\,2\pi i f_E(z)\,dz\), where \(f_E\)
is the normalized newform and \(\omega_E\) a Néron differential.
Set \(B_E=\deg(\phi)/m_\phi^2\). Pullback of the area form gives
\[
 B_E=\frac{4\pi^2\|f_E\|_{\Gamma_0(N)}^2}{P(E)}
     =\frac{N L^{\rm nv}(1,\operatorname{Ad} f_E)}{2\pi P(E)}.
\tag{SP1}
\]
Here \(L^{\rm nv}\) is the adjoint Euler product with the factor
set to one at primes \(p^2\mid N\); it retains the correct factors
elsewhere. In this weight-two, trivial-character setting this
can also be seen directly by unfolding the Rankin integral.
The residue of the Eisenstein series at the infinity cusp is
\(3/(\pi[\mathrm{PSL}_2(\mathbb Z):\overline{\Gamma_0(N)}])\).
The residue of \(\sum a_n^2 n^{-(s+1)}\) at \(s=1\) is
\(L^{\rm nv}(1,\operatorname{Ad} f_E)\zeta(2)^{-1}
 \prod_{p\mid N}(1+p^{-1})^{-1}\).
Indeed the local series at a multiplicative prime equals
the adjoint Euler factor, and at an additive prime it is one.
Compared with \(\zeta(s)L^{\rm nv}(s,\operatorname{Ad}f_E)/
\zeta(2s)\), both cases remove the factor \(1+p^{-s}\).
The index \(N\prod_{p\mid N}(1+p^{-1})\) in the Eisenstein
residue cancels that product, yielding (SP1).
[Watkins, opening modular-degree formula](https://arxiv.org/abs/math/0408126).
[Diamond–Flach–Guo, Lemma 8.7](https://arxiv.org/abs/2512.02348).

In particular \(B_E\) is independent of the parametrization.
Let \(E_0\) be the optimal classical quotient, and
\(\delta_{1,N}\) its modular degree. Choose a rational
isogeny \(E_0\to E\) of degree \(d\le163\), and write \(b\)
for its differential multiplier on Néron differentials.
The multiplier of the dual isogeny gives \(b\mid d\).
If \(c_0\) is the Manin constant of the optimal quotient,
then
\[
 B_E=\frac{d}{b^2c_0^2}\,\delta_{1,N}.
\tag{SP2}
\]
The conductor of the Frey family is squarefree away from \(2\),
so Pasten's bound on \(c_0\) makes the numerator and denominator
of the ratio in (SP2) uniformly bounded.
[Pasten, §3 and Corollary 10.2](https://arxiv.org/abs/1705.09251).

Watkins's lower bound for the primitive adjoint at \(1\) is
\(\gg1/\log(2N)\). The adjoint conductor is \(\ll N^2\)
here. The usual convexity interpolation between
\(-1/\log(3+N)\) and \(1+1/\log(3+N)\), with the Euler
bound on the right and the functional equation on the left,
gives an upper bound by a fixed power of \(\log(3+N)\).
The sole possible modified Euler factor is at \(2\), and its
value and reciprocal at \(1\) are bounded by the rank and the
unit-modulus Euler parameter bound. The bounded-conductor
cases can be absorbed into the constants. Therefore
\[
 \log\delta_{1,N}=h+A+O(\log(h+2)).
\]
This is a coefficient-two form of the classical Frey
modular-degree equivalence, using the exact Frey period
estimate. Comparing the twisted value with (SP1) also gives
\[
 \Theta_E:=\frac{N X_E}{B_E}\in K^\times,\qquad
 \big|\log|\Theta_E|\big|\ll\log(A+2).
\tag{SP3}
\]
Indeed \(\Theta_E=2G(\bar\chi)^2
L(1,\operatorname{Ad} f_E\otimes\chi)/
L^{\rm nv}(1,\operatorname{Ad}f_E)\), and the twisted edge estimate
is already uniform. This archimedean estimate concerns the signed
norm and allows cancellation between its finite prime valuations.
[Watkins, Lemma 3.4](https://arxiv.org/abs/math/0408126).
[Twisted edge estimate and period size](ABC_FIXED_CUBIC_ADJOINT_VALUE.md#a-two-sided-bound-at-the-edge).

## Exact rational local packages away from \(2,3,7\)

For \(p\mid R\) outside \(S\), put \(n_p=2v_p(abc)\).
Let \(C_p=\prod_{\lambda\nmid p}c_{p,\lambda}\) be the integer
off-diagonal package at the geometric prime \(p\), and let
\(\mathcal D_p\) be its diagonal factor product. In the notation
of the Tate calculation, \(u_{p,\ell}\in\mathbb Z_\ell\)
is the upper entry of the arithmetic Frobenius lift on the
Tate module before taking the adjoint.

When \(\chi(p)=1\), define
\[
 g_p=\prod_{\ell\ne p}
 \ell^{\min\{v_\ell(2n_p),v_\ell(p-1),v_\ell(2u_{p,\ell})\}} .
\]
This positive integer divides \(\gcd(2n_p,p-1)\) and is
independent of the compatible Tate coordinates. In this case
\[
 C_p=(n_{p,p'}g_p)^2,\qquad
 p^2 C_p\mathcal D_p=(n_pg_p)^2,
\tag{SP4}
\]
where \(n_{p,p'}\) is the part of \(n_p\) prime to \(p\).
Indeed \(n_p\) is even, so the relations
\(2n_pe=0\) and \(n_pd=n_p^2e\) from inertia give the
direct sum \((O/(2n_p))e\oplus(O/n_p)d\) at every coefficient
prime \(\ell\ne p\), including \(2\).
For Tate twist \(1\), Frobenius minus one maps this sum to
its first factor by the row \((p-1,-2u_{p,\ell})\).
The kernel has order
\(|O/n_p|\cdot|O/(2n_p,p-1,2u_{p,\ell})|\).
Multiplication over the coefficient places above \(\ell\)
gives the square of the corresponding \(\mathbb Z_\ell\)-index.
Equation (DT2) contributes \(p^{2(v_p(n_p)-1)}\), giving (SP4).

When \(\chi(p)\ne1\), the corresponding formula is
\[
 C_p=b_p,\qquad p^2 C_p\mathcal D_p=b_p,\qquad
 b_p=3^{\mathbf 1_{3\mid n_p}}\gcd(n_p,p^2+p+1).
\tag{SP5}
\]
For odd \(\ell\ne3,p\), the second diagonal entry of the
Frobenius-minus-one matrix is \(\chi(p)-1\), a unit. The
first is \(p\chi(p)-1\). If \(\ell\mid p^2+p+1\), the two
coefficient places split, and exactly one contributes
\(\ell^{\min(v_\ell(n_p),v_\ell(p^2+p+1))}\).
Otherwise neither contributes. At \(\ell=2\) both diagonal
entries are units on the two cyclic factors, so the kernel
is trivial. At \(\ell=3\), in the quadratic ramified
coefficient ring \(v_\lambda(\chi(p)-1)=1\), and
\(v_\lambda(p\chi(p)-1)\) is one if \(p\equiv1\bmod3\)
and zero otherwise. If \(3\mid n_p\), then \(v_\lambda(n_p)\ge2\)
and the full determinant valuation contributes. If \(3\nmid n_p\)
the inertia torsion is trivial. The \(3\)-valuation of \(p^2+p+1\)
is exactly one for \(p\equiv1\bmod3\) and zero otherwise.
This gives the additional factor of \(3\) in (SP5).
The diagonal product at \(p\) is \(p^{-2}\) in this case.
For example, \((1,124,125)\) has \(n_5=6\) and \(b_5=3\).
The triple \((1,31^3-1,31^3)\) has \(n_{31}=6\) and
\(b_{31}=9\). Both contributions come from the coefficient
place above \(3\).
[Local Tate inertia matrix, (CT2)](ABC_CUBIC_ADJOINT_TAMAGAWA_FACTORS.md).
[Diagonal multiplicative factors, (DT2)](ABC_CUBIC_ADJOINT_DIAGONAL_FACTORS.md).

Write \(N_S\) for the part of \(N\) supported on \(S\) and
\(C_S=\prod_{p\in S}C_p\) for the off-diagonal product there.
All good primes outside \(S\) contribute one. Their removal gives
the explicit positive integer
\[
 \mathcal B_E:=N^2 C\mathcal D
   =N_S^2 C_S
   \prod_{\substack{p\mid R,\ p\notin S\\ \chi(p)=1}}(n_pg_p)^2
   \prod_{\substack{p\mid R,\ p\notin S\\ \chi(p)\ne1}}b_p.
\tag{SP6}
\]
In particular the Tate exponents at the primes that split in
the fixed cubic field are present with their full square, even
when \(p\mid n_p\).

## Uniform Shimura transfer

Let \(N=DM\) be admissible with \(D\) odd and \(M\) divisible
by at least two odd primes. Thus \(D\) is a squarefree product
of evenly many primes of multiplicative reduction. Write
\(\delta_{D,M}\) for the modular degree of the optimal
elliptic quotient of \(J^D_0(M)\) in the isogeny class of \(E\),
and \(J_D=\prod_{p\mid D}n_p\).
Pasten's comparison theorem gives
\[
 \delta_{1,N}=\gamma_D\,\delta_{D,M}J_D,
 \qquad h_{\mathbb Q}(\gamma_D)=O(\omega(D)).
\]
Here \(h_{\mathbb Q}\) is rational logarithmic height.
The numerator and denominator of \(\gamma_D\) are supported
on primes at most \(163\); the hypotheses on \(M\) are used
for the denominator. This is a statement for the original
Frey curve, whose discriminant exponent at odd \(p\mid D\)
is \(n_p\). In the displayed formula use \(\gamma_1=1\).
[Pasten, Theorem 6.1(b.2) and Corollary 6.2](https://arxiv.org/abs/1705.09251).

Combine this equality with (SP2), put
\(\Gamma_D=d\gamma_D/(b^2c_0^2)\), and take the quadratic norm
in (SP3). The exact bridge is
\[
 Z_E=\delta_{D,M}^{\,2}\,
             \Gamma_D^{\,2}\,\mathrm N_{K/\mathbb Q}(\Theta_E)
             \frac{J_D^{\,2}}{\mathcal B_E}.
\tag{SP7}
\]
All signed valuations in this identity are over \(\mathbb Q\).
The height of \(\Gamma_D\) is \(O(1+\omega(D))\), with support
in a fixed set of primes independent of the triple and factorization.

If \(D\) is supported on primes outside \(S\) with \(\chi(p)=1\),
then \(J_D^2\mid\mathcal B_E\) by (SP6). One can also include
primes outside \(S\) with \(\chi(p)\ne1\) and \(v_p(abc)=1\).
Their only contribution to \(J_D\) is a power of \(2\).
For these choices (SP7) gives, outside the fixed support of
\(\Gamma_D\) and \(2\),
\[
 v_\ell(Z_E)-v_\ell(\mathrm N_{K/\mathbb Q}\Theta_E)
       \le 2v_\ell(\delta_{D,M}).
\tag{SP8}
\]
At the archimedean place these choices give the uniform
upper comparison
\[
 \tfrac12\log Z_E\le
 \log\delta_{D,M}+O\big(1+\omega(D)+\log(A+2)\big).
\]
The splitting condition is exactly what allows the local
Tate exponent in a Shimura transfer to be divided out of
the cubic package. For other choices (SP7) retains the
factor \(J_D^2/\mathcal B_E\).

Without a splitting restriction, the degree comparison,
\(\log\delta_{1,N}=h+A+O(\log(h+2))\), and
\(\sum_{p\mid R}\log n_p=O(h/\log(h+2))\) give
\[
 \log\delta_{D,M}
       = h+A+O\!\left(\frac{h}{\log(h+2)}\right)
\tag{SP9}
\]
uniformly for all the admissible factorizations just described.
The bound on the sum follows from the elementary divisor
bound used in the off-diagonal product estimate, and
\(\omega(D)\) is also bounded by that divisor logarithm.
By absorbing the \(o(h)\) error with a slightly smaller
\(\varepsilon\), the sharp uniform bound for any one such
choice of factorization per triple is
\(\log\delta_{D,M}\le(2+\varepsilon)A+O_\varepsilon(1)\).
It is equivalent to abc on the family with at least two
odd primes in \(R\).

## A small-level projector target

The restriction to two odd primes loses no cases of abc.
With only one odd prime \(p\) in the support, pairwise
coprimality of \(a,b,c\) reduces the nontrivial case to
\(p^v=2^u+1\) or \(2^u=p^v+1\). An odd \(v>1\) gives an
odd factor greater than one in \(p^v-1\) or \(p^v+1\).
An even \(v\) in the first equation makes \(p^{v/2}-1\)
and \(p^{v/2}+1\) consecutive even powers of two differing
by two, hence \(p^v=9\). In the second equation even \(v\)
forces \(2^u\equiv2\bmod8\), which is impossible there.
For \(v=1\) we have \(c\le p+1\). The case with no odd
prime has just \((1,1,2)\).

For a triple with \(k\ge2\) distinct odd primes in \(R\),
keep the two smallest odd primes in \(M\) when \(k\) is
even and the three smallest when \(k\) is odd, along with
the whole \(2\)-part of \(N\). Put every other prime in \(D\).
This is admissible and
\[
 \log M\le \frac{3}{k}\log R+O(1).
\]
The inequality follows by comparing the average logarithm
of the kept primes with that of all \(k\) primes. Thus a
projector-denominator bound \(\delta_{D,M}\ll_\varepsilon
R^{2+\varepsilon}\) for this prescription alone would
prove abc. The converse follows from (SP9). The degree
\(\delta_{D,M}\) is exactly the least positive integer
clearing the denominator of the rational Hecke projector
onto the Frey quotient in the endomorphism ring of the
Shimura Jacobian. In particular the target allows \(M\)
supported on at most three odd primes even when the radical
has unbounded support.
[Pasten, Lemma 5.6 on the projector denominator](https://arxiv.org/abs/1705.09251).

Pasten's effective height bound for admissible factorizations
is \(h_{\mathrm F}(E)\ll M\varphi(D)\log\ell\) for a prime
\(\ell\nmid N\). For every factorization the expression is at least a constant
times \(N/\log\log(3+N)\), since \(M\varphi(D)
=N\prod_{p\mid D}(1-p^{-1})\ge\varphi(N)\).
The optimized bound is thus still larger than logarithmic. Thus the missing gain concerns the
Frey projector's integral denominator or a sharper treatment of
its congruences. The relative identities above remove the
specified Tate factors before making that comparison.
[Pasten, Theorem 1.1](https://arxiv.org/abs/2609.17390).

