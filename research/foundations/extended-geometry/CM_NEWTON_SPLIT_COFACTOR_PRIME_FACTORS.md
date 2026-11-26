# Prime factors of the relative Frobenius cofactors

Let \(B/K\) be a polarized abelian threefold whose geometric rational
endomorphism algebra is a totally real cubic field \(F\), with all
endomorphisms defined over \(K\). At a good absolute degree-one prime
\(v\) of characteristic \(p\), write \(X^2-a_vX+p\) for its relative
Frobenius polynomial. Put
\[
 U_v=|N_{F/\mathbb Q}(a_v)|,\qquad
 V_v=-N_{F/\mathbb Q}(a_v^2-4p).
\]
The Weil bounds give \(U_v\le 8p^{3/2}\) and \(0<V_v\le64p^3\).
Equality in a conjugate Weil bound for \(a_v\) would imply
\(a_v^2=4p\) in the odd-degree field \(F\), which is impossible.
For a positive integer \(u\) write \(u_{\mathrm{sf}}\) for the
positive squarefree representative of its rational squareclass.
Write \(\omega_y(u)\) for the number of distinct prime divisors
\(\ell\le y\), and \(\omega(u)\) for the full count.

We count in any nonempty union of Frobenius classes in a fixed
finite Galois extension of \(K\), retaining the degree-one primes.
In particular this includes the condition that \(p\) split
completely in \(F\). The events \(U_v=0\) have density zero in
every such set.

**Theorem 1.** Unconditionally, take
\[
 y=\exp\bigl((\log\log X)/\log\log\log X\bigr),\qquad
 L_3=\log\log\log X .
\]
On primes of norm at most \(X\) in the chosen set with \(U_v\ne0\),
the vector
\[
 \frac{1}{\sqrt{L_3}}
 \bigl(\omega_y((U_v)_{\mathrm{sf}})-L_3,+       \omega_y((V_v)_{\mathrm{sf}})-L_3\bigr)
\]
converges in distribution to two independent standard normal
variables. The same holds for \(U_v,V_v\) in place of their
squarefree parts. The number of common prime divisors at most
\(y\) is bounded in probability.

**Theorem 2.** Assume GRH for the Dedekind zeta functions of the
absolute normal closures of the torsion fields of \(B\) composed
with the fixed conditioning field. Then the same joint normal law
holds for \(\omega\) and \(L_2=\log\log X\) in place of
\(\omega_y\) and \(L_3\). The same holds for \(U_v,V_v\) in place
of their squarefree parts. The number of common prime divisors is
\(o(\sqrt{L_2})\) in probability.

The normal law is a relative-trace counterpart of the
Erdős–Kac theorem for the rational first Frobenius trace in
Cojocaru, Davis, Silverberg and Stange, Theorem 3 and §5.
Their theorem assumes open image in the full symplectic
similitude group and uses GRH. Here the image is the
common-determinant real-multiplication group. Finite-level
counts compare the norm of the relative trace with the
squarefree norm of its quadratic discriminant, including at
auxiliary primes that do not split in the coefficient field.
[Cojocaru–Davis–Silverberg–Stange](https://www.alinacarmencojocaru.com/WIN2-IMRN2017-published.pdf).

## Local probabilities

Let \(C_0\) be the chosen background condition, defined using a
fixed finite Galois extension \(K_0/K\). Omit a fixed finite set
of auxiliary primes, including \(2\), the ramification in \(F\),
and the integral endomorphism order index. For every remaining
\(\ell\), the full \(\ell\)-adic image is
\[
 G_\ell=\left\{(g_\lambda)_{\lambda\mid\ell}\in
      \prod_{\lambda\mid\ell}\mathrm{GL}_2(\mathcal O_{F_\lambda}):
      \det g_\lambda=d\in\mathbb Z_\ell^\times
      \text{ is independent of }\lambda\right\}.
\]
These groups occur independently, also independently of the
projection to \(\operatorname{Gal}(K_0/K)\). Choose a finite
extension of \(K_0\) over which the Tate systems are independent
by Serre. Lombardo's image theorem applies over this field at
every sufficiently large \(\ell\), since all geometric
endomorphisms were defined over \(K\). Over \(K\), after fixing
any background Galois element, varying by the image of the
absolute Galois group of this extension gives all choices in
the product of the large local factors.
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703)
· [Serre, Theorem 1 and §3.1](https://www.college-de-france.fr/media/jean-pierre-serre/UPL8570957518086251784_independance.pdf).

Write \(r_\ell=\#\{\lambda\mid\ell:f_\lambda=1\}\). A uniform
element of \(G_\ell\) reduced modulo \(\ell\) has uniform common
determinant \(d\in\mathbb F_\ell^\times\). Given \(d\), its
components are independent and uniform in their determinant
fibers. Over the residue field of a component, of order \(q\),
the fixed-determinant mass of trace \(t\) is
\[
        \frac{q+\chi(t^2-4d)}{q^2-1},
\]
where \(\chi\) is the quadratic character extended by zero.
For the count, take diagonal entries \(x,t-x\) and solve
\(bc=x(t-x)-d\) for the off-diagonal entries. A nonzero
right-hand side gives \(q-1\) choices and zero gives \(2q-1\).
There are \(1+\chi(t^2-4d)\) zeros as \(x\) varies.

Modulo \(\ell\), the divisibility indicators \(I_\ell\) of \(U_v\)
and \(J_\ell\) of \(V_v\) test whether at least one component
has trace zero or discriminant zero, respectively. At a
degree-one component the trace-zero probability is
\(1/\ell+O(\ell^{-2})\) uniformly in \(d\). The probability
of discriminant zero is \(2\ell/(\ell^2-1)\) if \(d\) is a
square and zero otherwise. Its average is
\(\ell/(\ell^2-1)\). Higher-degree components contribute
\(O(\ell^{-2})\) to either event. Two prescribed components
are independent given \(d\) and both hit probabilities are
uniformly \(O(\ell^{-1})\). The two tests cannot hold in a
single component because \(d\ne0\). Consequently
\[
 \mathbb E I_\ell=r_\ell/\ell+O(\ell^{-2}),\qquad
 \mathbb E J_\ell=r_\ell/\ell+O(\ell^{-2}),\qquad
 \mathbb E(I_\ell J_\ell)=O(\ell^{-2}).                 \tag{1}
\]

Both events \(\ell^2\mid U_v\) and \(\ell^2\mid V_v\) have mass
\(O(\ell^{-2})\) at level \(\ell^2\). A hit of residue degree
at least two or simultaneous hits at two places have this
bound by the same calculation. In a single degree-one
component, fix the determinant lift and fix \(\bar g\)
modulo \(\ell\). For a nonscalar \(\bar g\), the trace lifts
are uniform in that determinant fiber because the two
linear forms on first-order corrections
\(X\mapsto\operatorname{Tr}X\) and
\(X\mapsto\operatorname{Tr}(\operatorname{adj}(\bar g)X)\)
are independent. If the residual trace is zero, one trace
lift is zero modulo \(\ell^2\). If the residual discriminant
is zero, the residual trace is nonzero, and one trace lift
has discriminant zero modulo \(\ell^2\) by the derivative
\(2\operatorname{Tr}\bar g\). There are at most two scalar
matrices for each determinant, of total mass \(O(\ell^{-3})\).

## Equidistribution with growing torsion levels

For a product \(m\) of auxiliary prime powers, let
\(M_m/\mathbb Q\) be the absolute normal closure of
\(K_0K(B[m])\). Its degree is \(O_{B,K_0}(m^c)\) for a
fixed \(c\), from its inclusion in the compositum of the
finitely many conjugate torsion fields over the fixed normal
closure of \(K_0\). Each factor has degree bounded by
\(|\mathrm{GL}_6(\mathbb Z/m\mathbb Z)|\).
Ramification is confined to a fixed set and primes dividing
\(m\). Thus
\(\log |\operatorname{disc}M_m|\ll_{B,K_0}[M_m:\mathbb Q]\log(2m)\)
by the local different bound
\(v_\ell(\operatorname{disc}M_m)\le n(1+v_\ell(n))\),
where \(n=[M_m:\mathbb Q]\). For the bound, pass locally
to the totally ramified extension of degree \(e\) over
the maximal unramified subextension. In the derivative
of an Eisenstein polynomial of a uniformizer, the leading
term has uniformizer valuation \(e-1+e v_\ell(e)\).
The derivative's nonzero terms have distinct valuations
modulo \(e\), so the leading term bounds the different
exponent above. Taking the norm yields the stated bound.

One can apply effective Chebotarev over \(\mathbb Q\) and
keep the exact multiplicity for degree-one places of \(K\).
In \(G=\operatorname{Gal}(M_m/\mathbb Q)\), let
\(H=\operatorname{Gal}(M_m/K)\). For an unramified rational
prime with Frobenius \(\sigma\), its degree-one places in
\(K\) correspond to cosets \(gH\) fixed by \(\sigma\).
The Frobenius over \(K\) at such a place is represented by
\(g^{-1}\sigma g\in H\). The sum of any bounded
conjugacy-invariant test over these places is thus a class
function on \(G\). Its mean on \(G\) is the original mean on
\(H\) by interchanging the sum over \(\sigma\) with the sum
over cosets. In particular the mean of a product of
divisibility indicators and the indicator of \(C_0\) is
the background mass times the product of the local means
at the distinct auxiliary primes.

Unconditionally, if \(m\le(\log X)^{o(1)}\), the degree,
the log discriminant and the root discriminant of \(M_m\)
are all \((\log X)^{o(1)}\). Theorems 1.3–1.4 of
Lagarias–Odlyzko apply with their degree/discriminant
threshold and their exceptional-zero bound for an absolute
normal extension. Both the exceptional term and the
error for each conjugacy class are
\(O(X\exp(-c'(\log X)^{1/3}))\) for large \(X\),
with the same bound after decreasing \(c'\) and
summing over classes. Under GRH, their Theorem 1.1 gives
error \(X^{1/2+o(1)}\) in the same counting problem for
\(m\le X^{o(1)}\). Removing primes of norm at most
\(X^{1/2}\) at most adds \(O_K(X^{1/2})\) to each count.
These statements are uniform for the indicator tests just
described, since each class-function multiplicity is at most
\([K:\mathbb Q]\).
[Lagarias–Odlyzko, Theorems 1.1, 1.3 and 1.4](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf).

## The joint limit

Put \(\mu=\log\log y\). By Chebotarev for the fixed normal
closure of \(F\) and partial summation,
\(\sum_{\ell\le y}r_\ell/\ell=\mu+O_F(1)\).
The mean number of fixed embeddings in the transitive
permutation action is one. The prime-counting error for
this fixed extension, including any fixed exceptional zero,
integrates to a bounded reciprocal-prime error. Hence in
the independent model of (1), both summed means and both
summed variances are \(\mu+O(1)\), and the summed covariance
is \(O(1)\).

Fix real \(s,t\) and take the \(k\)-th moment of the sum of
the centered variables
\(s(I_\ell-\mathbb EI_\ell)+t(J_\ell-\mathbb EJ_\ell)\),
normalized by \(\sqrt{\mu}\). In the independent model, any
prime that occurs only once in the expansion gives zero.
A prime occurring at least twice has moment
\(O_{k,s,t}(1/\ell)\). Terms with fewer than \(k/2\) distinct
primes therefore contribute \(o(1)\). When \(k\) is even,
the terms with exactly \(k/2\) distinct primes are the
pairings. Their normalized sum tends to
\((k-1)!!(s^2+t^2)^{k/2}\), since their local variance
sums to \((s^2+t^2)\mu+O(1)\), and excluding repeated
pair indices changes the result by \(o(1)\).
The resulting moments determine the centered normal law
with variance \(s^2+t^2\). This argument for all \(s,t\)
gives the joint normal law.

For the unconditional cutoff in Theorem 1, each term uses
\(m\le y^{2k}=(\log X)^{o(1)}\). Dividing the Chebotarev
error by the background count and summing over the
\(O_k(y^k)\) terms in any fixed moment gives \(o(1)\).
Thus the same moments hold on the prime set. The shift
from \(\mu\) to \(L_3\) is \(o(\sqrt{L_3})\).
Omitting a fixed set of auxiliary primes changes the
counts by a bounded amount. For \(U_v=0\) all the trace
divisibility indicators equal one. The moment comparison
includes these primes and the variance bound forces them
to have density zero.

To pass to squarefree parts, the differences in the
truncated counts are bounded by the numbers of
\(\ell^2\)-divisibility events. Each such count has bounded
expectation by the local \(O(\ell^{-2})\) bounds and uniform
Chebotarev at modulus \(\ell^2\). The same applies to the
count of common prime divisors by (1). Markov's inequality
gives tightness of these counts, which proves Theorem 1.

Under the GRH of Theorem 2, use
\(y=X^{1/(\log\log X)^{1/4}}\). At this cutoff
\(y^{2k}=X^{o(1)}\) for each fixed \(k\), so the same
moment comparison and bounded expectations hold. The
shift from \(\mu\) to \(L_2\) is \(o(\sqrt{L_2})\).
The Weil bounds imply that the number of prime divisors
larger than \(y\) of either nonzero integer is
\(O((\log\log X)^{1/4})\), uniformly for \(p\le X\),
because their product divides this integer. This also
bounds the tail for the squarefree parts and for the
common divisor count. It proves Theorem 2.

## Independent cofactor squareclasses

**Corollary.** Fix a finite subgroup \(H_0\subset
\mathbb Q^\times/\mathbb Q^{\times2}\). Outside a relative
density-zero subset of the split degree-one primes, at
every prime with exactly one failed factor and \(a_v\ne0\),
the squareclasses of the trace cofactor \(b_v=N(a_v)/p\)
and discriminant cofactor \(d_v=V_v/p\) are linearly
independent modulo \(H_0\). The fraction of exceptions up to
\(X\), measured in the ambient split prime set, is
\(O_{B,K,H_0}(\log\log\log X/\log\log X)\).
For the relation \([b_vd_v]\in H_0\) alone the bound is
the square of this fraction bound. Under the GRH of
Theorem 2, with the split condition as background, these
bounds hold with \(\log\log X/\log X\) in place of
\(\log\log\log X/\log\log X\).

At one failed split factor, exactly one coefficient place
divides \(a_v\), and \(a_v^2-4p\) has valuation one there
and zero at the other coefficient places. Thus \(b_v,d_v\)
are integers. The invariants at issue differ from those
for \(U_v,V_v\) only in sign and at the single prime \(p\).
Omit the finite union of supports of \(H_0\) from the
auxiliary primes. For the individual classes a witness
to nonmembership in \(H_0\) is \(v_\ell(U_v)=1\) or
\(v_\ell(V_v)=1\), respectively. The corresponding model
probability is \(r_\ell/\ell+O(\ell^{-2})\) by (1) and
the level-\(\ell^2\) calculation. For the product class,
use a witness with one valuation equal to one and the
other zero. The two orders of this condition are disjoint
and have total probability \(2r_\ell/\ell+O(\ell^{-2})\).
These witnesses are unchanged by dividing by \(p\)
on the retained range \(p>X^{1/2}\), \(\ell\le y<X^{1/2}\).

For each of the three witness tests let \(\rho_\ell\)
be its model probability, and \(j=1\) or \(2\) as appropriate.
Here use the new cutoff \(y=\exp(\kappa\log\log X/
\log\log\log X)\) unconditionally, or
\(y=\exp(\kappa\log X/\log\log X)\) under GRH, with
\(\kappa>0\) fixed and sufficiently small. Put \(\mu=\log\log y\)
and take the even integer \(k=2\lceil20\mu\rceil\).
The even Bonferroni sum of order \(k\) bounds the
fraction with no witness at the indicated primes.
In the independent model it differs from the full product
\(\prod_{\ell\le y}(1-\rho_\ell)\) by at most the next
elementary symmetric sum, bounded above by
\(\bigl(\sum_{\ell\le y}\rho_\ell\bigr)^{k+1}/(k+1)!\).
Since \(\sum\rho_\ell=j\mu+O(1)\), the full product is
\(O((\log y)^{-j})\) and the remainder is
\(O((\log y)^{-3})\) by \(h!\ge(h/e)^h\).

There are at most \((y+1)^k\) terms in this sum and
each uses a torsion modulus \(m\le y^{2k}\). In the
unconditional case these upper bounds are
\((\log X)^{(40+o(1))\kappa}\) and
\((\log X)^{(80+o(1))\kappa}\), respectively.
Choose \(\kappa\) so that the latter modulus range gives
degree, log discriminant and root discriminant at most
\((\log X)^{1/10}\), using the fixed polynomial bounds in
the equidistribution argument. Then the threshold in
Lagarias–Odlyzko Theorem 1.3 is met and Theorem 1.4
still absorbs the exceptional zero with an exponential
saving in \((\log X)^{1/3}\). After summing classes and
Bonferroni terms the error is \(o((\log y)^{-2})\) relative
to the background count. Under GRH the two displayed
bounds for terms and moduli hold with \(X\) in place of
\(\log X\). Choose \(\kappa\) so that the term count and
each field degree are at most \(X^{1/20}\).
Theorem 1.1 gives at most \(O(X^{1/2}n\log X)\) error
for the weighted sum over classes of a field of degree
\(n\), since the local class multiplicity is bounded
and \(\log |\operatorname{disc}M_m|\ll n\log(2m)\).
Summing all terms still saves a power of \(X\).
Consequently the failure fractions are \(O((\log y)^{-j})\),
giving the stated bounds. The density of the entire
one-failure locus remains undetermined.
