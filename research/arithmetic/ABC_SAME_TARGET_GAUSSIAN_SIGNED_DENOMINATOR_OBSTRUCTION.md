# Exact-power targets with large signed denominators

An exact powerful target can support several primitive square-full
representations whose free prime supports intersect only at fixed
primes. At the natural radical scale of the target, almost all
members of a two-parameter family retain a large denominator in
every signed word of positive scalar powers. The construction
tests the role of the two free-letter caps in the scalar moment
argument.

## A power-scale lower bound

Fix integers \(r,m\ge2\) and \(K>0\). Let \(\mathcal P\) be any
prime pool with product at most \(X^K\), prescribed for the whole
family at height \(X\). For each \(\varepsilon>0\) and all
sufficiently large \(X\) there are
\[
                  \gg_{r,m,K,\varepsilon} X^{1/r-\varepsilon}
\]
distinct exact \(r\)-th powers \(X/2<N\le X\) with \(m\)
distinct unordered primitive representations by positive squares.
Write these representations in some order as \(A_i+B_i=N\),
with \((A_i,B_i)=1\), and set
\[
           f_i=F(A_i/N),\qquad
           F(z)=z-z^{-1}+(1-z)^{-1}.
\]
For each such target the representations can be chosen to satisfy
the following simultaneously. For every nonempty
\(I\subseteq\{1,\ldots,m\}\), all positive integers \(q_i\), and
all nonzero rational \(w_i\) whose reduced numerators and
denominators are supported on \(\mathcal P\), one has
\[
   \operatorname{rad}_{(\mathcal P\cup\operatorname{supp}N)^c}
    \left(\operatorname{den}\left(\sum_{i\in I} w_i f_i^{q_i}\right)\right)
                      \ge X^{1/r-\varepsilon}.                 \tag{GT1}
\]
Here \(\operatorname{den}\) means the reduced positive denominator.
All bounds are uniform in the prescribed pool. The parameters
\(r,m,K,\varepsilon\) are fixed before \(X\) tends to infinity.

## Representations with controlled support intersections

Choose distinct rational split primes \(p_1,\ldots,p_j\) in
\(\mathbb Z[i]\) with \(2^j\ge m\). Put \(b_0=\prod_{h=1}^j p_h\)
and choose one Gaussian prime \(\pi_h\) above each \(p_h\).
Choose \(m\) distinct conjugation patterns. For each pattern form
\(\beta_i\) by choosing one of \(\pi_h^r,\bar\pi_h^r\) at each
\(h\) and taking their product. All \(\beta_i\) have norm \(b_0^r\).
By unique factorization their pairwise quotients are not Gaussian
units.

Choose a fixed \(1<\kappa<2^{1/(2r)}\) and take
\[
      T=\frac{X^{1/(2r)}}{\sqrt{2b_0}\kappa}.
\]
Use parents \(x+iy\) with \(T<x,y\le \kappa T\),
\[
           x\equiv0\pmod{2b_0},\qquad
           y\equiv1\pmod{2b_0},\qquad (x,y)=1 .
\]
There are \(\gg_{r,m}T^2\) parents. In fact Möbius inversion
on common divisors coprime to \(2b_0\) gives a positive main
coefficient proportional to
\(\prod_{p\nmid 2b_0}(1-p^{-2})\) and error \(O_{r,m}(T\log T)\).

Write \(\alpha=(x+iy)^r\), \(z_i=\alpha\beta_i\), and let
\(A_i=(\Re z_i)^2\), \(B_i=(\Im z_i)^2\). Their target is
\[
                           N=[b_0(x^2+y^2)]^r
\]
in the required shell. The parent and its conjugate are relatively
prime in \(\mathbb Z[i]\), since a common Gaussian prime would
divide both \(2x\) and \(2iy\), and the parent has odd norm.
Its norm is coprime to \(b_0\). Also \(\beta_i\) is coprime to
its conjugate since it uses only one choice at each split prime.
These facts show that \(z_i\) is coprime to its conjugate.
In particular its two coordinates are coprime and nonzero.
Equality of two unordered square pairs would make \(z_i\) a
unit times \(z_{i'}\) or \(\bar z_{i'}\). The first conflicts
with the distinct patterns. The second conflicts with the
nonempty prime support of \(\alpha\), since \(\alpha\) is
coprime to its conjugate and to the norm of every \(\beta_i\).

To control support overlap, write \(\alpha=s+it\) and
\(\beta_i=u_i+iv_i\). The two coordinates of \(z_i\) are
the linear forms with coefficient vectors \((u_i,-v_i)\)
and \((v_i,u_i)\) in \((s,t)\). Vectors from distinct rows
all have the same Euclidean length. If two were proportional,
they would differ only by sign, making their \(\beta\)'s
Gaussian unit multiples. Thus each of the four determinants
between distinct rows is nonzero. Let \(D\ge1\) be the
product of their absolute values over all row pairs.
If a prime divides any coordinate of \(z_i\) and any
coordinate of \(z_{i'}\), it divides the corresponding
determinant, because \((s,t)=1\). Hence all overlaps
of free supports between rows are contained in the fixed
prime support of \(D\).

## Counting free radical support and signed denominators

Enlarge \(\mathcal P\) by the primes dividing \(D\), to form
\(\mathcal P^+\). Fix a row index \(i\) and \(0<\eta<1\).
The number of parents for which
\[
       \operatorname{rad}_{(\mathcal P^+)^c}(A_iB_i)
                                       \le T^{2-\eta}
\]
is at most \(T^{2-\eta+o(1)}\). Enumerate the outside
radical at cost at most \(T^{2-\eta}\), then enumerate
the possible complete products \(A_iB_i\) on the resulting
polynomial-height prime pool, at subpower cost. The latter
enumeration is uniform in \(\mathcal P\) for fixed \(K\).
Each product determines its pair \(A_i,B_i\) up to
divisor-many choices. The pair determines \(z_i\) up to
sign choices for its two coordinates, then \(\alpha\) by
division by \(\beta_i\), and the parent up to at most
\(r\) Gaussian roots. This proves the bound.
[Uniform enumeration on polynomial-height prime pools](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Discard these parents for each of the fixed \(m\) indices.
The remaining parents still number \(\gg T^2\). Each value
of \(N\) has at most \(4\tau(x^2+y^2)=T^{o(1)}\) parents,
by Gaussian factorization of the norm. Thus there are
\(T^{2-o(1)}\) distinct targets represented by the retained
parents. For any chosen slot \(i\) of a retained parent its
free radical outside \(\mathcal P^+\) is greater than
\(T^{2-\eta}\). All these primes are absent from both \(N\)
and the free letters of every other slot.

The formula
\[
           f_i=\frac{A_i^2B_i+(A_i-B_i)N^2}{A_iB_iN}
\]
is in lowest terms. Its numerator is a unit at each prime
dividing its denominator, by primitivity and \(A_i+B_i=N\).
For any active slot \(i\), a prime of \(A_iB_i\) outside
\(\mathcal P^+\) therefore gives \(w_i f_i^{q_i}\) a pole
at that prime, whereas all other active summands are integral
there. Its valuation is strictly smaller than the valuations
of the other summands, so it survives in the reduced denominator
of the sum. It is outside the original pool and target support.
Taking \(\eta>0\) small in terms of \(r,\varepsilon\) proves
(GT1) and the count of targets. The same parents work for all
choices of the indicated powers and weights.

The target is already \(r\)-full, so its least repair to order
\(r\) is \(1\). The square-full summands are primitive even
before repair. The construction imposes the \(1/r\) radical
budget on the target alone. An extraction of low-denominator
words for the small-budget exceptional set must use the
additional original radical bounds on the free letters.
[Scalar moment bounds with two free radical caps](ABC_SCALAR_RECIPROCAL_POWER_AMPLIFICATION.md#statement)
· [One-target least repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Comparison with shifted-product estimates

Fang, He and Zhu bound multiplicative energy averaged over
translates of one common set using rich-line estimates for its
Cartesian square and diagonal pins. Their result concerns
the complete translates of that set. A same-target exceptional
fibre selects its own subset of each translate, and the radical
caps do not give a small product-set estimate for these inputs.
In fact a set of integers \(S\subseteq[1,X]\) satisfies
\(|S^k|\ge |S|^k X^{-o(1)}\) for fixed \(k\) by the ordered
divisor bound. This gives the same fact for the fractions
\(a/N\) of a fixed target. Thus a lower bound on the exact
product set does not by itself limit that target's multiplicity.
[Fang–He–Zhu, Theorems 2–4 and §4](https://arxiv.org/html/2611.30174v1)

Bloom's quantitative sparse-prime results give lower bounds
on sum or product expansion. In particular their asymptotic
forms permit \(o(\log |S|/\log\log |S|)\) distinct primes per
element. The radical cap by itself gives only a bound of order
\(\log |S|/\log\log |S|\) when \(|S|\) has a fixed positive
height exponent. The exact product-set issue above already
applies even on families where the stronger sparse-prime
hypothesis is available.
[Bloom, Theorems 7–8](https://arxiv.org/html/2512.04931v1)

Zhang, Zeng and Wang obtain mixed-sumset expansion from a
low finite-valuation rank coset cover. For a subgroup of
\(\mathbb Q^\times\) with valuation rank \(d\), select \(d\)
prime coordinates that are injective on its valuation space.
Each coset then contains at most \(2(O(\log X))^d\)
rational labels of polynomial height in \(X\), by restricting
the possible integer values of those valuations. Thus for a
dictionary of size \(X^s\), with fixed \(s>0\), the range
\(d=o(\log|S|/\log\log|S|)\) requires \(|S|^{1-o(1)}\)
cosets. The actual complete labels give no fixed-power
compression of the sort required by that mixed-sumset input.
[Zhang–Zeng–Wang, Theorem 1.4](https://arxiv.org/html/2610.36930v1)
