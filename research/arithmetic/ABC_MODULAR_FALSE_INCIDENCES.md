# False modular incidences with all three radical caps and bounded heights

Fix an integer \(d\ge3\), a positive integer \(k\), and
positive scales \(\beta_1,\ldots,\beta_k\) with
\[
                      1\le\beta=\sum_i\beta_i<3/2.
\]
Take distinct prime tests \(q_i\asymp N^{\beta_i}\), independently
in fixed dyadic intervals apart from the distinctness restriction,
and let \(Q=\prod_i q_i\).
There are dictionaries of actual integers of heights
\(\asymp_d N^d\) with all radicals at most \(2N\), such that for a
positive proportion of the test tuples there are
\(\gg_{d,k,\boldsymbol\beta} N^{3-\beta}\) distinct three-radical
records of primitive triples satisfying all the congruences
\(A+B\equiv C\bmod q_i\), with the **same representatives**.
There is no true additive equality anywhere in the dictionaries.
This applies at arbitrarily small fixed budgets \(\lambda=3/d\).
It is a lower-cap obstruction for total test modulus up to
\(N^{3/2-\eta}\), not for the separate large-product regime
\(\prod_i q_i>N^{2+\eta}\).

## Dictionaries and the homogeneous cone Fourier estimate

Put \(s=1/(100d)\). Take integer roots in three positive
intervals \(I_1,I_2,I_3\) of lengths \(\asymp_d N\),
with
\[
 I_1=I_2=[N,(1+s)N],\qquad
 I_3=[(1+2s)N,(1+3s)N].
\]
Use
\[
                   A=x^d,\qquad B=y^d,\qquad C=z^d .
\]
Thus \(C>\max(A,B)\) but \(A+B>C\) for every candidate,
since \((1+3s)^d<2\). When \(x,y,z\) are squarefree and
pairwise coprime, the complete labels are pairwise coprime
and their exact radicals are \(x,y,z\), all \(\asymp_d N\).
In fact
\(\operatorname{rad}(ABC)=xyz<z^3=C^{3/d}\).
For any fixed \(M\ge d\), minimal target repair gives
\(H_M(C)=z^M\), so these candidates also have the correct
target repair height scale.

For a prime \(q\nmid d\), let \(\mathcal F_q\) be the set
of points of \(x^d+y^d=z^d\) in \(\mathbf F_q^3\), with all
three coordinates nonzero. The number of projective points with all coordinates
nonzero is \(q+O_d(\sqrt q)\). Normalize \(z=1\) and
expand the number of \(d\)-th roots of each unit \(t\)
as \(\sum_{\chi^d=1}\chi(t)\). The resulting count is
\(\sum_{\chi^d=\psi^d=1}J_q(\chi,\psi)\). Its principal
term is \(q-2\), and the other \(O_d(1)\) terms have
absolute value at most \(\sqrt q\).
Consequently the allowed affine cone has size
\(q^2+O_d(q^{3/2})\).
[Jacobi sum conventions and bounds](ABC_MODULAR_JACOBI_INCIDENCE.md#jacobi-sums-and-the-removal-of-bounded-order-modes)

For any nonzero frequency vector \(t\in\mathbf F_q^3\),
the complete Fourier coefficient of the allowed cone
has magnitude \(O_d(q)\). Indeed it is a union of
nonzero scalar multiples of projective points, and the
scalar exponential sum is \(q-1\) when the projective
point lies on the line \(t\cdot(x,y,z)=0\), and \(-1\)
otherwise. The Fermat curve is smooth for \(q\nmid d\), hence
geometrically irreducible by Bézout. The line therefore
contains at most \(d\) of its projective points. This gives the precise expression
\(q\,\#(\mathcal C_q^\circ\cap\{t\cdot x=0\})-\#\mathcal C_q^\circ\).
Removing the coordinate planes thus preserves the estimate.

Let \(\mathcal F_Q\) be the simultaneous allowed
cone condition at the \(k\) primes. By CRT its zero
Fourier mode is
\[
                 Q^2\left(1+O_{d,k}\left(\sum_i q_i^{-1/2}\right)\right).
\]
At any nonzero frequency \(t\bmod Q\), its Fourier
coefficient is at most
\(O_{d,k}(Q\,\gcd(t_1,t_2,t_3,Q))\), since each zero local
frequency contributes at most \(O_d(q_i^2)\) and each
nonzero one at most \(O_d(q_i)\).

It follows that for any three arithmetic progressions
\(J_i\) of sizes at most \(O(N)\), whose steps are units
modulo \(Q\), the number of triples in their product
satisfying \(\mathcal F_Q\) is
\[
 \frac{|J_1||J_2||J_3|}{Q}
 \left(1+O_{d,k}\left(\sum_i q_i^{-1/2}\right)\right)
       + O_{d,k}(Q\log^3(2Q)).                              \tag{FI1}
\]
Here \(Q\asymp N^\beta\), so \(Q\gg N\).
To verify the error, for each proper divisor \(h\mid Q\),
the sum of absolute one-coordinate interval transforms
over nonzero frequencies divisible by \(h\) is
\(O((Q/h)\log(2Q))\), by the geometric progression bound.
The zero coordinate frequency has size at most \(O(N)\).
Consequently the sum of the three-coordinate products
over \(t\ne0\) with \(h\mid t_1,t_2,t_3\) is at most
\[
 O\!\left((Q/h)\log(2Q)
                 (N+(Q/h)\log(2Q))^2\right).
\]
Counting a frequency with \(\gcd(t_1,t_2,t_3,Q)=h\)
in this larger set, multiplying by
\(O(Qh)/Q^3\), and summing over the at most \(2^k\)
proper divisors proves the error in (FI1).
Intervals may be viewed as multisets of residues if
\(Q\) is only comparable to \(N\); the geometric bound
still holds.

## Preserving squarefreeness and primitivity

Take a large fixed cutoff \(P_0\). Impose that each
root is squarefree modulo \(p^2\), and that at most
one root is divisible by \(p\), for each \(p\le P_0\).
The allowed CRT density of triples modulo the product
of these \(p^2\)'s is bounded below uniformly in
\(P_0\): there is an allowed class for each small prime,
and the excluded local density is \(O(p^{-2})\).
Apply (FI1) separately to all fixed residue classes.
The total count before the larger-prime exclusions is
at least \(c_d N^3/Q\) for all sufficiently large \(N\),
with \(c_d>0\) independent of \(P_0\).
The error here is harmless because \(2\beta<3\).

Choose a fixed small \(\rho>0\) with
\(\rho<\min(1/2,\min_i\beta_i)\) and \(\rho+2\beta<3\), and put
\(U=N^\rho\). For each prime \(P_0<p\le U\) and each
square-divisibility condition \(p^2\mid x\) or common
condition \(p\mid x,y\), apply (FI1) again, now simply
on the unrestricted interval boxes with that condition.
The steps \(p,p^2\) are units at all test primes.
Its main term is \(O(N^3/(Qp^2))\); the accumulated
error is \(O_{d,k}(UQ\log^3 Q)=o(N^3/Q)\).
Summing over these primes and over the six types of
conditions gives a loss
\(O_{d,k}(N^3/(QP_0))+o(N^3/Q)\).

For the remaining tail \(p>U\), count the bad triples
first as ordinary integer triples without any modular
condition. A square divisor \(p^2\mid x\) removes at
most \(O(N/U+\sqrt N)\) possibilities for that root.
A common prime of two roots removes at most
\(O(N^2/U)\) pairs, by summing \((N/p)^2\) for
\(U<p\ll N\). Thus at most
\[
                       O(N^3/U+N^{5/2})
\]
ordinary tuples are bad in this tail.
For every such tuple the integer
\(D=x^d+y^d-z^d\) is nonzero and \(O_d(N^d)\).
It has at most \(O_{d,\min\beta_i}(1)\) prime factors
in the ranges allowed for the tests, so it contributes
to only boundedly many test tuples. There are
\(\asymp_{\boldsymbol\beta,k}N^\beta/(\log N)^k\)
distinct prime test tuples. The average larger-prime
tail loss is therefore
\[
 O_{d,k,\boldsymbol\beta}
       \left(\frac{(N^3/U+N^{5/2})(\log N)^k}{N^\beta}\right)
                  =o(N^{3-\beta}).
\]
With \(P_0\) sufficiently large first, Markov's inequality
shows that a positive proportion of test tuples retain
\(\gg N^{3-\beta}\) distinct primitive modular
records. Each such record has \(x,y,z\) as its exact
radical labels and has no integer equality, as asserted.
The argument bounds averaged tail loss and does not
require a prime-tuple sieve beyond the elementary
prime number theorem.



## Comparison with the generic false-positive upper bound

For this family, for every fixed tuple of polynomial prime-test
scales, regardless of the size of the product, the average number
of false incidences is at most
\[
 O_{d,k,\boldsymbol\beta}\!
       \left(N^{3-\beta}(\log N)^k\right).
\]
All the integer differences are nonzero, and each has boundedly
many eligible prime factors, so the same divisor argument as in the
sieve tail now applies to all \(O(N^3)\) ordinary triples. Thus in
the range of the lower theorem the generic false-positive exponent
\(3-\beta\) is actually attained. In the full occurrence reduction
one must instead work at product exponent greater than two to make
this noise a power below \(N\). In that regime this example is
consistent with a power saving: it isolates the false-incidence
barrier for testing with smaller aggregate modulus.

