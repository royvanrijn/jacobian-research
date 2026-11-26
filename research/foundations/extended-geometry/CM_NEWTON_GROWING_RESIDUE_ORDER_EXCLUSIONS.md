# Growing residue-order exclusions in the real-cubic case

Fix a threefold \(A/K\) with exact real-cubic multiplication by
\(F\), defined over \(K\). Use good absolute degree-one places
of residue characteristic \(p\), excluding the finite set of
ramified primes and primes dividing the endomorphism-order index
or a polarization degree. Write \(a_v\in O_F\) for its relative
Frobenius trace,
\[
         s_1=\operatorname{Tr}_{F/\mathbb Q}(a_v),\qquad
         b_v=N_{F/\mathbb Q}(a_v)/p.
\]
The Weil bound gives \(|s_1|\le6\sqrt p\). In the sole
degree-one failure cases below, \(b_v\) is a nonzero
integer with \(|b_v|\le8\sqrt p\). Write
\(h_v=f(K_v(A[p])/K_v)\) for the residue degree of the full
local \(p\)-division field. The split case of its local formula is
proved in the [residue-prime torsion calculation](CM_NEWTON_RESIDUE_PRIME_TORSION_RAMIFICATION.md#small-residue-degrees-in-the-full-torsion-field).

## Split primes with relations between the ordinary eigenvalues

**Proposition (Split residue-order interpolation).** At split places
with exactly one failed factor, let \(u_1,u_2\in\mathbb F_p^\times\)
be the two ordinary residue eigenvalues. On \(X/2<p\le X\) the
count with \(h_v\le H\) and \(\operatorname{ord}(u_2/u_1)\le R\),
excluding \(s_1\in\{-2,-1,0,1,2\}\), is
\[
                              O_K(\sqrt X\,H^2 R).
\]
Consequently \(h_v\le p^\theta\) and
\(\operatorname{ord}(u_2/u_1)\le p^\rho\) cut out relative density
zero in the full split degree-one class when \(2\theta+\rho<1/2\).
There is also a bound \(O_K(\sqrt X\,H^2 |J|)\) with the same
five exceptions, if the ratio condition is replaced by
\(u_2=u_1^j\) for some \(j\) in a prescribed finite set \(J\)
of integers. The bound is uniform in the choices of exponents.

Put \(d=h_v\) and \(m=\operatorname{ord}(u_2/u_1)\). Then
\(m\mid d\mid p-1\). Since the third trace vanishes modulo \(p\),
\(s_1\equiv u_1+u_2\bmod p\). The roots of unity of orders
dividing \(d\) lift their residues by Hensel's lemma, so
\(p\mid Q_{d,m}(s_1)\) for
\[
 Q_{d,m}(Z)=
       \prod_{\substack{\zeta^d=1\\\operatorname{ord}(\xi)=m}}
                       \bigl(Z-\zeta(1+\xi)\bigr)\in\mathbb Z[Z].
\]
Its roots are bounded by two, and the degree is \(d\varphi(m)\).
Its coefficients are algebraic integers fixed by Galois.
A nonzero value at \(|s_1|\le6\sqrt X\) has at most
\(O(d\varphi(m))\) prime divisors in the interval, by bounding
their product by \(|Q_{d,m}(s_1)|\le(6\sqrt X+2)^{d\varphi(m)}\).
An exact zero forces the integer \(s_1\) to have absolute value
at most two. Sum over the \(O(\sqrt X)\) integers \(s_1\) and use
\[
 \sum_{d\le H}\sum_{\substack{m\mid d\\m\le R}}d\varphi(m)
   =\sum_{m\le\min(R,H)}
               m\varphi(m)\sum_{j\le H/m}j
   \ll H^2R .
\]
There are at most \([K:\mathbb Q]\) places above a rational prime.
For the exponent condition use instead
\(\prod_{\zeta^d=1}(Z-\zeta-\zeta^j)\), of degree \(d\) for
each \(j\in J\). Summing over \(d\le H\) gives the asserted bound.

The five exact trace loci have density zero. At an auxiliary
split prime the fixed-determinant trace distribution in each
factor has point probabilities \(O(1/\ell)\). By conditioning on
the other two traces and the determinant, \(s_1=j\) has the same
upper bound on the full residual image. The class-wise effective
argument in the [trace cofactor sieve](CM_NEWTON_TRACE_COFACTOR_SIEVE.md#counting-small-trace-cofactors)
only needs this full-image fraction after summing over classes,
so also gives the bound
\[
 O_{A,K}\!\left(\frac{X(\log\log X)^{1/6}}{(\log X)^{7/6}}\right)
\]
for the five loci. The Chebotarev count of the full split class is
of order \(X/\log X\), and higher absolute residue degrees are
negligible. Thus the dyadic bound gives relative density zero
also for moving thresholds with \(H^2R=o(\sqrt X/\log X)\).
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703)
· [Thorner–Zaman, Theorem 1.4](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-p.pdf#page=6).

Taking \(R=H\) recovers the previous \(\theta<1/6\) bound. The
stronger \(\theta<1/4\) range for bounded \(R\) requires the
extra ratio hypothesis. A polynomial containing all nonzero complex
sums of two \(d\)-th roots of unity must have degree of order
\(d^2\). For a nonzero sum \(z\), the equations
\(|w|=|z-w|=1\) allow at most two values of \(w\), giving just one
unordered pair. Requiring the two roots to generate the full
cyclic group leaves \(d^2\prod_{\ell\mid d}(1-\ell^{-2})\) ordered
pairs before excluding antipodal pairs. Hence primitivity alone
retains degree of order \(d^2\) in this polynomial argument.

## Coefficient splitting type \(1+2\)

**Proposition (Nonsplit residue-order exclusion).** Suppose \(p\)
has coefficient splitting type \(1+2\) and the degree-one factor
is the sole failed coefficient factor. Let \(u\in\mathbb F_{p^2}^\times\)
be the ordinary residue eigenvalue on the étale quotient at the
degree-two coefficient place. Then
\[
                             h_v=\operatorname{lcm}(2,\operatorname{ord}(u)).
\]
On \(X/2<p\le X\) the count with \(h_v\le H\), outside the same
five exact trace loci, is \(O_K(\sqrt X\,H^2\log(2H))\).
When the \(1+2\) class occurs, \(h_v\le p^\theta\) gives relative
density zero in that class for any fixed \(\theta<1/4\).
All other nonordinary \(1+2\) places of sufficiently large
characteristic lie in the zero trace locus.

If the degree-two coefficient place fails, the trace norm is
divisible by \(p^2\), which exceeds the Weil bound \(8p^{3/2}\)
for \(p>64\). Hence \(a_v=0\) in this case. Otherwise the ordinary
height-four factor has étale \(p\)-torsion of rank one over
\(\mathbb F_{p^2}\), with arithmetic Frobenius \(u\).
This gives an unramified quotient of order \(\operatorname{ord}(u)\).
The nonordinary height-two factor has the nonsplit Cartan as
its inertia image with normalizer as decomposition image,
giving a residue quotient of order two. Both orders divide \(h_v\).
[Serre, §1.9, Proposition 9 and §1.11, Proposition 12](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf).

For the opposite divisibility, take \(U/\mathbb Q_p\) unramified
of degree \(\operatorname{lcm}(2,\operatorname{ord}(u))\), prime to
\(p\). Over \(U\) the ordinary étale quotient is the constant
group \((\mathbb Z/p)^2\), and the connected group is \(\mu_p^2\)
by Cartier duality and the polarization. Lifts of a basis of
the quotient give torsors under \(\mu_p^2\) with classes in
\(O_U^\times/O_U^{\times p}\). Over \(U(\mu_p)\) they generate
an elementary abelian Kummer extension. Restriction to this
tame base is injective on the unit classes, and each nonzero
combination is represented by \(q=1+p\beta\) with \(\beta\in O_U^\times\).
The Eisenstein and different calculation in the
[residue-prime ramification theorem](CM_NEWTON_RESIDUE_PRIME_TORSION_RAMIFICATION.md#the-full-torsion-field)
shows that each nontrivial cyclic subextension is totally
ramified. The whole ordinary extension is therefore totally
ramified of \(p\)-power degree, since any unramified residue
quotient of an elementary abelian group supplies an unramified
cyclic subextension. The nonordinary field is totally tame after
the even unramified base change, of degree \(p+1\) over \(U(\mu_p)\).
The compositum has ramification index divisible by the ordinary
and tame indices. They are coprime, so the compositum is totally
ramified. This proves the formula for \(h_v\).

For the count put \(d=\operatorname{ord}(u)\le H\) and
\(e=p\bmod d\). Then \(e^2=1\bmod d\). The two ordinary residues
are \(u,u^p\), so \(s_1\equiv u+u^p\bmod p\).
Lifting \(u\) to a primitive \(d\)-th root gives
\[
 p\mid
   \prod_{\operatorname{ord}(\zeta)=d}(s_1-\zeta-\zeta^e).
\]
The corresponding monic integer polynomial in \(s_1\) has
degree \(\varphi(d)\) and roots bounded by two. For each \(s_1,d,e\)
the number of prime divisors in the interval is \(O(\varphi(d))\)
unless the value is zero, which forces \(|s_1|\le2\).
There are at most \(2^{\omega(d)+1}\) possible \(e\) for each
\(d\), by the prime-power solutions of \(e^2=1\). Since
\(2^{\omega(d)}=\sum_{j\mid d}\mu^2(j)\), summing gives
\[
 \sum_{d\le H}\varphi(d)\,2^{\omega(d)}
           \le H\sum_{d\le H}2^{\omega(d)}
           \ll H^2\log(2H).
\]
This proves the count and gives density zero for
\(H^2\log(2H)=o(\sqrt X/\log X)\) on dyadic intervals.
The \(1+2\) condition is a union of Frobenius conjugacy classes
in the coefficient normal closure compositum over \(K\), of
positive density when nonempty. The auxiliary-prime argument for
the five exact trace loci applies to every splitting type.

The factor \(\log(2H)\) can be replaced by
\(\log(2\min(R,H))\) under the extra condition
\(\operatorname{ord}(u^p/u)\le R\). Indeed
\(m=\operatorname{ord}(u^p/u)=d/\gcd(d,e-1)\) divides \(d\).
Knowing \(m\) determines \(e\) at every odd prime-power divisor
of \(d\), since it determines whether the local solution to
\(e^2=1\) is \(+1\) or \(-1\). There are at most four choices
at the power of two. The degree sum is therefore at most
\(4\sum_{d\le H}d\sum_{m\mid d,\,m\le R}1
\ll H^2\log(2\min(R,H))\).

## Orders of integral trace cofactors

**Proposition (Small trace-cofactor order).** At split places with
exactly one failed factor, the count with
\(\operatorname{ord}_{\mathbb F_p^\times}(b_v)\le H\) on
\(X/2<p\le X\), outside \(b_v=\pm1\), is \(O_K(\sqrt X\,H^2)\).
In particular the order bound \(p^\theta\) with
\(\theta<1/4\) gives relative density zero in the split
degree-one class.

For sufficiently large \(p\), \(0<|b_v|\le8\sqrt p<p\), so its
residue is a unit. If the exact order is \(d\le H\), then
\(p\mid\Phi_d(b_v)\), where \(\Phi_d\) is the cyclotomic
polynomial, because \(p\nmid d\). For each of the \(O(\sqrt X)\)
possible integers \(b_v\) and each \(d\), the nonzero value
has at most \(O(\varphi(d))\) prime divisors in the interval
by \(|\Phi_d(b_v)|\le (8\sqrt X+1)^{\varphi(d)}\).
The only exact zeros at an integer are \(b_v=\pm1\).
The sum of \(\varphi(d)\) for \(d\le H\) is of order \(H^2\),
and these two exceptional values have density zero by the
[trace cofactor bound](CM_NEWTON_TRACE_COFACTOR_SIEVE.md#counting-small-trace-cofactors).
Using primitive cyclotomic polynomials saves the imprimitive
factors, and the bound gives density zero whenever
\(H^2=o(\sqrt X/\log X)\).

At a split one-failure place let \(\kappa\) be the residue of
\(a_v/p\) in the failed \(\mathbb Q_p\) coefficient factor.
In the other factors the residues are \(u_1,u_2\), so
\(b_v\bmod p=\kappa u_1u_2\).
The torsion residue degree involves \(u_1,u_2\). The cofactor
order also involves the next \(p\)-adic digit of
the failed trace. These exclusions leave open split one-failure
places with both multiplicative quantities large.
