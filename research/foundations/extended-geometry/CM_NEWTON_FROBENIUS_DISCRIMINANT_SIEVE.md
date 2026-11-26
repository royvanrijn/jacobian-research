# The discriminant cofactor at split real-cubic primes

Let \(B/K\) be an abelian threefold with geometric rational endomorphism
algebra the totally real cubic field \(F\), and suppose its endomorphisms
are defined over \(K\). Its real multiplication is compatible with the
Rosati involution of a polarization over \(K\). Let \(L\) be the normal
closure of \(F\).
We consider good primes \(v\) of \(K\) of absolute residue degree one,
with residue characteristic \(p\) split completely in \(F\). This is
the identity class for \(KL/K\) after discarding higher-degree places.
At such a prime write the relative Frobenius polynomial, its
discriminant, and the relative quadratic field as
\[
 X^2-a_vX+p,\qquad D_v=a_v^2-4p,\qquad E_v=F(\sqrt{D_v}).
\]
Write \(\mathfrak d_v=\mathfrak d_{E_v/F}\) and \(d_v=N_{F/\mathbb Q}\mathfrak d_v\).
The Weil bound and the odd valuation of \(4p\) at the unramified
coefficient primes imply that \(D_v\) is totally negative. In particular
\(E_v/F\) is quadratic. The order with basis \(1,\pi_v\), for a root
\(\pi_v\) of the relative polynomial, has relative discriminant
\((D_v)\). Consequently
\[
  -N_{F/\mathbb Q}(D_v)=d_v h_v^2\leq (4p)^3
\]
for an integer \(h_v\geq1\), the norm of the relative index ideal.
By the [local ramification calculation](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group),
exactly one failed pair at a split prime gives \(v_p(d_v)=1\). At
these primes let \(c_v\) be the positive squarefree part of \(d_v/p\).
The rational quadratic field \(\mathbb Q(\sqrt{N_{F/\mathbb Q}(D_v)})\)
is then \(\mathbb Q(\sqrt{-p c_v})\), since \(p\nmid c_v\). Its absolute
discriminant is either \(p c_v\) or \(4p c_v\).

**Theorem.** Assume GRH for the Dedekind zeta functions of the composita
\(L K(B[m])\) and their subfields, for integers \(m\geq1\). Given
\(0<\delta<1/38\), the one-failed-pair identity-class primes with
\(c_v\leq p^\delta\) have a power-saving counting bound.
More precisely, for every \(\theta\) with \(\delta<\theta<1/38\),
their count up to \(X\) is \(O_{B,K,\theta}(X^{1-\theta+\delta})\).
Under the same GRH assumption, the number of prime ideals of \(F\)
ramified in \(E_v/F\), normalized by subtracting \(\log\log p\) and
dividing by \(\sqrt{\log\log p}\), converges to the standard normal
distribution on the identity class.

The square-sieve argument adapts the use of torsion representations in
Frobenius-field counting. Burungale, Hida and Lai's theorem fixes the
isomorphism type of the entire Frobenius field and proves a uniform
upper bound for that field. The test here concerns the squareclass of
the norm discriminant after removing \(p\), so it can be summed over a
growing range of possible cofactors. The prime-factor count adapts Liu's
moment argument for the sizes of reductions of elliptic curves.
[Burungale–Hida–Lai, Theorem 1.3 and Remark 1.5](https://arxiv.org/pdf/2402.07935)
· [Liu, Theorems 1 and 3](https://doi.org/10.1016/j.jnt.2005.10.014).

## Effective auxiliary distributions on the identity class

For all sufficiently large \(\ell\) and integers \(e\geq1\) the image
on \(B[\ell^e]\) is the group of tuples over
\(\mathcal O_{F,\mathfrak l}/\ell^e\) for \(\mathfrak l\mid\ell\), with
determinants equal to one common scalar in
\((\mathbb Z/\ell^e\mathbb Z)^\times\). Use Lombardo's theorem on each
factor over the finite independence extension in Serre's theorem.
One can then restrict further to the compositum of that extension
with \(KL\). Indeed the image of a finite-index subgroup in a compact
product of the Tate module images is closed of finite index, hence
open and contains the unrestricted product outside a finite set of
factors. It follows that the kernel of the projection to
\(\operatorname{Gal}(KL/K)\) still surjects onto any finite product of
the large auxiliary images. Thus imposing the identity class gives
precisely the product uniform distribution on those groups, with no
change in the individual marginal distributions. The same argument
allows the projective quotients used below.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703)
· [Serre, Theorem 1 and §3.1](https://www.college-de-france.fr/media/jean-pierre-serre/UPL8570957518086251784_independance.pdf).

For a product modulus \(m=\prod\ell^e\) of sufficiently large primes,
the full image has size at most \(m^{10}\). At each \(\ell\), the scalar
determinant accounts for one power of \(\ell^e\), and the three
coefficient degrees counted with residue degrees account for nine
powers. At split \(\ell\), the image in
\(\operatorname{PGL}_2(\mathbb F_\ell)^3\) has size at most \(\ell^9\).
For the projective product at squarefree split modulus \(m\), use the
bound \(m^9\) instead. Set \(b=10\) for full images and \(b=9\) for
these projective ones. The compositum with \(KL\) has absolute degree
\(n\ll_{K,L}m^b\) and is unramified outside a fixed finite set and the
primes dividing \(m\). Its logarithmic absolute discriminant is
\(O_{B,K}(n(\log(2m))^2)\). In fact at a ramified rational prime
\(\ell\) the exponent of the absolute discriminant is at most
\(n(1+\log_\ell n)\), by the local different bound
\(e-1+e v_\ell(e)\) for ramification index \(e\). The latter follows
by taking an Eisenstein uniformizer over the maximal unramified
subextension and bounding the valuation of its minimal-polynomial
derivative by that of its leading term. Distinct derivative terms
have distinct valuations modulo \(e\).

Effective Chebotarev for a bounded class function on this finite
group gives, after summing over conjugacy classes, an unnormalized
error at most
\[
 O_{B,K}\bigl(X^{1/2}m^b(\log(2Xm))^2
                  +m^{2b}(\log(2m))^2\bigr).
\]
The first term uses the class-size weights in the GRH error. For the
second term the number of conjugacy classes is at most the group
order. Removing higher-degree places costs \(O_K(X^{1/2})\). This
bound applies up to \(X\) or on \((X/2,X]\). After imposing the identity
class and normalizing by its count, which is of order \(X/\log X\),
it estimates the difference between the empirical mean and the product
group mean. All tests below have bounded absolute value independent of
the modulus. [Lagarias–Odlyzko, Theorem 1.1](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf).

## A square sieve for the cofactor

For split \(\ell\) write \(\chi_\ell\) for the quadratic character of
\(\mathbb F_\ell\), extended by zero. At a tuple of matrices with
traces \(t_i\) and common determinant \(d\), set
\[
 H_\ell =\chi_\ell\!\left(-d^{-1}\prod_{i=1}^3(t_i^2-4d)\right).
\]
This is a class function of the projective tuple. Multiplying a matrix
by a scalar changes its discriminant and determinant by squares,
and the determinant squareclass is the same on all three factors.
On a Frobenius for \(v\) of residue characteristic different from
\(\ell\), this is \(\chi_\ell(-N(D_v)/p)\).

For fixed determinant \(d\), the number of matrices in
\(\operatorname{GL}_2(\mathbb F_\ell)\) of trace \(t\) is
\(\ell(\ell+\chi_\ell(t^2-4d))\). For two distinct eigenvalues
this follows by dividing by the split or nonsplit torus centralizer.
For a repeated eigenvalue there are one scalar matrix and a Jordan
class of size \(\ell^2-1\). Since
\(\sum_t \chi_\ell(t^2-4d)=-1\), the average of the discriminant
character in one factor is \(-r_d/(\ell^2-1)\), with \(r_d\) the
number of solutions of \(t^2=4d\). To see the character sum, count
the \(\ell-1\) pairs satisfying \((t-u)(t+u)=4d\). Given the
common determinant the three matrices are independent, so
\(|\mathbb E H_\ell|\ll \ell^{-6}\). This is also its mean on the
projective image.

Work in the residue interval \(X/2<p\leq X\). Choose
\(Y=X^\theta\), and let \(\mathcal P\) be the rational primes in
\([Y,2Y]\) split completely in \(F\). Its size is of order
\(Y/\log Y\). For each positive integer \(c\leq X^\delta\), form
\[
   S_v(c)=\sum_{\ell\in\mathcal P}\chi_\ell(c)H_\ell(v).
\]
If \(-N(D_v)/p\) has squareclass \(c\), then every term equals one
except possibly those with \(\ell\mid cN(D_v)\). There are
\(O_\theta(1)\) of those by the Weil bound. Thus
\(S_v(c)\geq |\mathcal P|-O_\theta(1)\) on the given squareclass.

The diagonal terms in the mean of \(S_v(c)^2\) contribute at most
\(|\mathcal P|\). For distinct \(\ell,\ell'\), the product mean
is \(O(Y^{-12})\). The normalized Chebotarev error for these two
projective factors is
\[
 O_{B,K,\theta}\bigl(X^{-1/2+18\theta}(\log X)^3
                          + X^{-1+36\theta}(\log X)^3\bigr).
\]
When \(\theta<1/38\), the sum over \(O(|\mathcal P|^2)\) pairs
is \(o(|\mathcal P|)\). The second moment is therefore
\(O_{B,K,\theta}(|\mathcal P|)\), uniformly in \(c\). Each such
squareclass has proportion \(O_{B,K,\theta}(|\mathcal P|^{-1})\),
or count \(O_{B,K,\theta}(X^{1-\theta})\) in the interval. The
one-failed-pair primes with \(c_v\leq p^\delta\) belong to these
classes, because \(-N(D_v)/p=(d_v/p)h_v^2\) has squareclass \(c_v\).
Summing over \(c\) and the dyadic residue intervals proves the counting
bound.

Consequently, proving \(c_v\leq p^\delta\) on all but a density-zero
subset of the one-failed-pair identity class, for any fixed
\(\delta<1/38\), would give the split density-one result under GRH.
Here density zero refers to the full set of identity-class primes.
The available pointwise archimedean bound \(d_v\leq(4p)^3\) gives
\(c_v\leq64p^2\).

For any fixed bound on \(c_v\), the density-zero conclusion is
unconditional. For a fixed \(c\), every auxiliary \(\ell\nmid c\)
requires either \(H_\ell=0\) or \(H_\ell=\chi_\ell(c)\). This event
has probability \(1/2+O(1/\ell)\), since the character mean is
\(O(\ell^{-6})\) and the zero set has probability \(O(1/\ell)\).
Apply the product distribution to as many fixed split primes as
needed, then use a finite union over the bounded cofactors.

## Number of ramified prime ideals

At an odd unramified auxiliary prime \(\ell\), a coefficient place
\(\mathfrak l\mid\ell\) ramifies in \(E_v/F\) exactly when
\(v_{\mathfrak l}(D_v)\) is odd. Use the proxy event
\(v_{\mathfrak l}(D_v)=1\), which is visible modulo \(\ell^2\).
Write \(q=\ell^{f(\mathfrak l/\ell)}\). Given a fixed scalar
determinant modulo \(\ell^2\), the residual determinant has either
zero or two repeated-root traces. The nonscalar matrices at each such
trace have proportion \(1/q\) in their fixed determinant group.
For a nonscalar matrix with a repeated root, varying the lift of its
trace with fixed determinant is surjective over \(\mathbb F_q\).
Indeed the two linear forms on matrix increments given by the trace
and differential of the determinant are independent unless the base
matrix is scalar. The discriminant has nonzero trace derivative at
these repeated-root traces because the determinant is a unit and
\(\ell\ne2\). Thus the proxy has conditional probability
\(2(1-1/q)/q\) when the residual determinant is square in
\(\mathbb F_q^\times\), and zero otherwise. Scalar matrices have
discriminant divisible by \(\ell^2\) for every lift, but they have
proportion \(O(q^{-3})\). Including nonscalar lifts, the probability
that \(v_{\mathfrak l}(D_v)\geq2\) is \(O(q^{-2})\).

Use one indicator \(U_\ell\) for at least one coefficient place
satisfying the proxy event. Conditional on the common determinant,
events at distinct coefficient places are independent. The chance of
two proxy events at the same rational prime is \(O(\ell^{-2})\).
For residue degree one the common determinant is square with
probability \(1/2\), while higher residue degrees contribute
\(O(\ell^{-2})\) in total. Thus if \(n_1(\ell)\) counts the degree-one
prime ideals of \(F\) above \(\ell\),
\[
  \mu_\ell=\mathbb E U_\ell=n_1(\ell)/\ell+O(\ell^{-2}).
\]
The prime ideal theorem and partial summation give
\(\sum_{\ell\leq Y}n_1(\ell)/\ell=\log\log Y+O(1)\).
Excluding higher-residue-degree ideals in this sum changes it by a
bounded amount.

Take \(Y=X^{1/(\log\log X)^{1/4}}\) and consider the empirical
distribution up to \(X\), leaving out the fixed exceptional auxiliary
primes and discarding the negligible residue primes with
\(p\leq X^{1/2}\). Any fixed moment of the centered sum
\(\sum_{\ell\leq Y}(U_\ell-\mu_\ell)\) differs by \(o(1)\) from
that for independent Bernoulli variables with these means. To see
this, expand the moment of order \(r\) and apply the full-image
Chebotarev estimate to each summand at modulus
\(m=\prod \ell^2\leq Y^{2r}=X^{o(1)}\). Excluding prime norms
\(\leq X^{1/2}\) preserves the error \(X^{-1/2+o(1)}\); there are
at most \(Y^r=X^{o(1)}\) summands. The independent sum has variance
\(\log\log Y+O(1)\), and its fixed-order cumulants are
\(O_r(\log\log Y)\). Its centered and normalized moments therefore
converge to the Gaussian moments, which determine the limit law.

Replacing the proxies by the actual count of ramified coefficient
prime ideals over \(\ell\leq Y\) introduces an error of bounded
expected absolute value. The only discrepancy comes from the events
\(v_{\mathfrak l}(D_v)\geq2\) and the multiple proxy events, with
summable local probabilities and vanishing total Chebotarev error.
The fixed excluded primes contribute at most a constant. The number
of ramified prime ideals over rational primes larger than \(Y\) is
at most \(3\log(4X)/\log Y=O((\log\log X)^{1/4})\), since their
product of norms divides \(d_v\). These errors vanish in probability
after normalization by \(\sqrt{\log\log X}\).
Both \(\log\log Y-\log\log X\) and \(\log\log p-\log\log X\)
are negligible on that scale, proving the stated limit law.

The count of ramified primes changes by just one if the characteristic
prime is removed in the one-failed-pair case. The normal law therefore
does not decide the density of this diagonal ramification event.
