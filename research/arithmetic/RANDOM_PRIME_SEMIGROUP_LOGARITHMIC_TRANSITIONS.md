# Logarithmic modifications of random prime semigroups

Always admit \(2\), select the other primes independently, and let \(C\)
be the generated multiplicative semigroup, including 1. Write \(r_k(n)\)
for unordered representations by \(k\) members of \(C\), with repetitions
allowed. The power-law threshold and its critical gap bound are in the
[preceding result](RANDOM_PRIME_SEMIGROUP_BASIS_THRESHOLD.md).

**Theorem.** Let \(0<\gamma<1\), \(\lambda>0\), and
\(p_q=\min(1,\lambda q^{-\gamma})\) for odd \(q\). Put
\(\alpha=1-\gamma\). Almost surely,
\[
 N_C(x):=|C\cap[1,x]|\ \sim\
 A_{\omega,\lambda} x^\alpha(\log x)^{\lambda-1},
 \qquad A_{\omega,\lambda}=\frac{G_\omega(\alpha)}{\alpha\Gamma(\lambda)}>0 .
 \tag{1}
\]
On \(\Re s>\alpha\) the semigroup Dirichlet series is
\(\zeta(s+\gamma)^\lambda G_\omega(s)\), where \(G_\omega\) is
holomorphic and nonvanishing on \(\Re s>\alpha/2\) almost surely.
The conclusion holds after any fixed finite changes to the prime
selection probabilities.

**Theorem.** Fix \(k\ge2\) and \(\gamma=1-1/k\) in the multiplier
model of (1). If \(0<\lambda<1\), \(kC\) almost surely has density zero.
At \(\lambda=1\), almost surely the lower density of the complement is at least
\(\exp(-\kappa_k A_{\omega,1}^k)\), where \(\kappa_k=\Gamma(1+1/k)^k/k!\).
If \(\lambda>1\), \(kC\) has density one almost surely. More precisely,
for all large \(n\),
\[
 \Pr(n\notin kC)\ll
 \exp\!\left\{-c\,\frac{(\log n)^{k(\lambda-1)}}
                         {(\log\log n)^{k\lambda}}\right\}.
 \tag{2}
\]
Thus \(kC\) is almost surely cofinite if \(\lambda>1+1/k\).
If \(1<\lambda<1+1/k\), almost surely the gap fraction in dyadic
blocks \((x,2x]\) is at least
\(\exp(-(\kappa_k A_{\omega,\lambda}^k+o(1))(\log x)^{k(\lambda-1)})\).
In particular every sufficiently large dyadic block contains gaps.
At \(\lambda=1+1/k\), the probability of cofiniteness is strictly between zero and one.

**Theorem.** Fix \(k\ge2\), put \(\gamma=1-1/k\), and instead use
\(p_q=\min(1,q^{-\gamma}(\log q)^\beta)\) for odd \(q\).
For \(\beta<0\), \(kC\) has density zero almost surely. At \(\beta=0\),
its complement has positive lower density almost surely.
For every \(\beta>0\), \(kC\) is almost surely cofinite. For every fixed
\(L>0\) it satisfies \(\Pr(n\notin kC)\ll_{k,\beta,L} \exp(-c_L(\log n)^L)\).

**Theorem.** In the multiplier model of (1), for each fixed \(k\ge2\)
and \(0<\gamma<1-1/k\), almost surely for all sufficiently large \(n\),
\[
 r_k(n)\gg
 n^{k-1-k\gamma}
 \frac{(\log n)^{k(\lambda-1)}}{(\log\log n)^{k\lambda}}.
 \tag{3}
\]
For every \(0<\gamma<1\) and every fixed modulus \(m\ge1\), the following
averaged representation law holds almost surely, for every residue \(r\):
\[
 \frac{\displaystyle\sum_{\substack{n\le x\\ n\equiv r\bmod m}} r_k(n)}
      {x^{k\alpha}(\log x)^{k(\lambda-1)}}
 \ \longrightarrow\
 \frac{\Gamma(1+\alpha)^k}{k!\,\Gamma(1+k\alpha)}
       A_{\omega,\lambda}^k\,\nu_m^{*k}(r).
 \tag{4}
\]
Here \(\nu_m\) is a random probability measure modulo \(m\). For \(q^f\Vert m\),
draw independent capped valuations \(J_q\). For a selected \(q\) their
probabilities are \((1-q^{-\alpha})q^{-j\alpha}\) for \(0\le j<f\) and
\(q^{-f\alpha}\) for \(j=f\). For an excluded \(q\), set \(J_q=0\).
Given these valuations, use the uniform law on the corresponding
residue classes, and take the additive convolution modulo \(m\).

## Counting and residue classes

Let \(X_q\) indicate selection, with \(X_2=1\). Away from finitely many
primes, the centered first-order terms of \(\log G_\omega(s)\) are
\((X_q-\lambda q^{-\gamma})q^{-s}\). Kolmogorov's maximal inequality
on dyadic intervals gives
\(\sum_{q\le u}(X_q-p_q)=O_{\omega,\epsilon}
(u^{\alpha/2+\epsilon})\), since the summed variances are
\(O(u^\alpha)\). Partial summation gives normal convergence on
\(\Re s>\alpha/2\). The higher prime-power terms are absolutely
locally summable almost surely there, by taking expectations of their
absolute majorants and using \(\gamma+2\Re s>1\).
The finitely many changed local terms, including the one at 2, are
holomorphic logarithms on \(\Re s>0\). Exponentiating gives the
nonvanishing \(G_\omega\), real and positive at \(\alpha\), and the
factorization on \(\Re s>\alpha\), where the prime sum converges
absolutely. Take the branch of \(\zeta(s+\gamma)^\lambda\) defined
there by the Euler product. Classical nonvanishing on \(\Re s=1\)
gives continuation locally across \(\Re s=\alpha\) away from \(\alpha\).
Near \(\alpha\), \(((s-\alpha)\zeta(s+\gamma))^\lambda\) has
analytic branch equal to 1 at \(\alpha\).
Apply Delange's theorem to the ordinary Laplace transform
\(F_\omega(s)/s\) of \(N_C(e^t)\). Its leading singularity has
coefficient \(G_\omega(\alpha)/\alpha\) and exponent \(\lambda\).
[Delange, Theorem III](sources/delange-ikehara-1954.pdf)
[Soprounov, Claims 1, 2 and 4](../foundations/extended-geometry/sources/soprounov_pnt_arithmetic_progressions.pdf)

For fixed \(m\), restrict to \(b\in C\) with \((b,m)=1\). Inserting a
Dirichlet character \(\chi\bmod m\) in this series gives a factor
\(L(s+\gamma,\chi)^\lambda\) times a holomorphic factor on
\(\Re s>\alpha/2\). The same variance and partial summation argument
applies with the bounded coefficients \(\chi(q)\), and character
orthogonality gives nonnegative series for the reduced classes.
Only the principal term has the singularity at \(\alpha\).
The same nonvanishing and Delange's theorem give, in each reduced
class, the asymptotic (1) times
\(\prod_{q\mid m,\,X_q=1}(1-q^{-\alpha})/\varphi(m)\).
This is obtained by first removing the selected Euler factors at
primes dividing \(m\) from the full principal series.

For a fixed tuple of allowed exponents \(e_q\ge0\) at \(q\mid m\),
put \(d=\prod_{q\mid m}q^{e_q}\) and write \(a=d b\) with \((b,m)=1\).
The mass of this tuple among \(a\in C\) tends to
\(\prod_{q\mid m,\,X_q=1}(1-q^{-\alpha})d^{-\alpha}\), by applying
the reduced-class result at \(x/d\). Multiplication by \(d\) takes
the uniform units to the uniform classes with the specified capped
valuations. The allowed exponent weights sum to 1, so truncating
them and then increasing the truncation proves the claimed residue
law. The assertions hold simultaneously at fixed moduli because
there are countably many characters and moduli.

For each residue, the counting measures of \(a/x\) in \([0,1]\),
normalized by \(x^\alpha(\log x)^{\lambda-1}\), thus converge
weakly to \(A_{\omega,\lambda}\nu_m(r)\alpha t^{\alpha-1}\,dt\).
Use the \(k\)-fold product measures on \(\sum_i t_i\le1\), whose
boundary has measure zero, and sum over residue vectors adding
to \(r\). The simplex integral is
\(\Gamma(1+\alpha)^k/\Gamma(1+k\alpha)\) by the beta integral.
Repeated coordinates contribute at most
\(\binom{k}{2}N_C(x)^{k-1}=o(x^{k\alpha}(\log x)^{k(\lambda-1)})\).
Dividing the ordered count by \(k!\) proves (4).

## Conditioning on smaller selected primes

Here \(0<\gamma\le1-1/k\); put \(u=k-1-k\gamma\). Choose a small constant
\(\theta>0\) such that the fixed-dimension lower sieve at level
\(n^{1/2}\) is positive for \(z=n^\theta\).
Choose \(0<\eta<\theta\) with \(k\eta<1/4\),
\((2k+k\gamma)\eta<1\) and \(k\eta\gamma<\theta/k\).
Choose a small fixed \(\tau<\eta\), specified below, and set
\(U=n^\tau\) and \(Y=(\log n)^B\), where \(B\alpha>
k|\lambda-1|+4\). Expose primes at most \(z\).
Allow only selected primes in \((Y,U]\) with \(q\nmid n\) to enter
pairwise coprime squarefree multipliers \(d_1,\ldots,d_k\le n^\eta\).
The permitted prime pool and the multiplier tuples are thus determined
by the exposure. Write \(P=\prod_i d_i\).

With conditional data satisfying an event of probability
\(1-O(\exp(-c(\log n)^2))\), the multiplier weight satisfies
\[
 \sum_{(d_1,\ldots,d_k)} P^{-\alpha}
       \gg_{k,\lambda,\gamma,B}
           \frac{(\log n)^{k\lambda}}{(\log\log n)^{k\lambda}}.
 \tag{5}
\]
In fact without the size bound its value is \(\prod_q(1+k q^{-\alpha})\)
over the allowed pool. The mean of \(\sum_q X_q q^{-\alpha}\) over the
deterministic eligible primes is
\(\lambda(\log\log U-\log\log Y)+O(1)\) by Mertens. Omitting the prime
divisors of \(n\) costs \(O(\log n/(Y\log Y))=o(1)\).
The mean of \(\sum_q X_q q^{-\alpha}\log q\) is \(O_\lambda(\log U)\).
The latter follows from \(\sum_{q\le U}(\log q)/q\ll\log U\), or by
partitioning successively at square roots in Mertens' product bound.
The respective variances are \(O(Y^{-\alpha})\) and
\(O(Y^{-\alpha}(\log Y)^2)\), by summing over integers the upper bounds
\(p_q q^{-2\alpha}\) and \(p_q q^{-2\alpha}(\log q)^2\).
The maximum summand sizes are \(O(Y^{-\alpha})\) and
\(O(Y^{-\alpha}\log Y)\). For a centered sum with variance \(V\) and
maximum summand bound \(b\), deviation by 1 in either direction has
probability \(O(\exp(-c/(V+b)))\). This follows by bounding each
exponential moment by \(\exp(O(t^2\mathbb E Z^2))\) when \(|tZ|\le1\)
and choosing \(t=c/(V+b)\), with a sufficiently small absolute \(c\).
Our choice of \(B\) makes both exceptional probabilities
\(O(\exp(-c(\log n)^{k|\lambda-1|+2}))\).
The logarithm of the pool product differs from
\(k\sum_q X_q q^{-\alpha}\) by \(O(Y^{-\alpha}\log\log n)\)
on this event. Under the normalized product weight, each prime is
independently assigned to one of the \(k\) multipliers, with individual
probability \(q^{-\alpha}/(1+kq^{-\alpha})\), or unassigned.
The expected \(\sum_i\log d_i\) is at most
\(k\sum_q X_q q^{-\alpha}\log q\). Take \(\tau\) small enough
that this is at most \((\eta/2)\log n\) on the same event, for
all large \(n\). Markov's inequality retains at least half the
product weight with all \(d_i\le n^\eta\), giving (5).
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

For each admissible multiplier vector, we count decompositions
\(a_i=d_i b_i\) of \(n\). As in the binary and higher-order proofs,
the first \(j=k-1\) summands lie in
\([n/k,n/k+n/(4k^2)]\) and the last is determined, so all have
size comparable to \(n\). Impose fixed congruences modulo
\(H=4\prod_{3\le q\le2k}q\) on the first \(j\) summands so that
all summands are units at odd primes up to \(2k\), and all are odd
except the last, which is \(2\bmod4\) when \(n\not\equiv k\bmod2\).
Such residues exist as any element modulo an odd prime is a sum of
two units. The remaining parity can be set by adjusting one odd
residue modulo 4. All multipliers are coprime to \(H\).

For completeness, write \(b_i=r_i+H v_i\) for \(i\le j\), with
\(r_i\) fixed by these congruences. The requirement \(d_k\mid
n-\sum_{i\le j}d_i b_i\) gives \(v_j=\ell(\mathbf v')+d_k w\),
where \(\mathbf v'=(v_1,\ldots,v_{j-1})\), and \(\ell\) is an integral
affine lift of the solution modulo \(d_k\). This uses \((H d_j,d_k)=1\).
Then \(b_k=(n-\sum_{i\le j}d_i b_i)/d_k\) is integral and affine
in the free variables \((\mathbf v',w)\). For any modulus
\(r\) coprime to \(H\), the free residue vectors correspond bijectively to
solutions of \(\sum_i d_i b_i=n\bmod r\). Given \(\mathbf v'\),
changing \(w\) changes \((b_j,b_k)\) by \((H d_k,-H d_j)\).
The congruence \(d_j b_j+d_k b_k=n-\sum_{i<j}d_i b_i\bmod r\)
has exactly \(r\) solutions, since \((d_j,d_k)=1\),
and the indicated vector has order \(r\).
For each residue vector of the free variables the number of box points
is
\[
                 N/r^j+O_k((n/r)^{j-1}),\qquad
                 N\asymp_k n^j/P,\quad 1\le r\le n.
\]
Indeed the first \(j-1\) free variables range over intervals of lengths
\(\asymp n/d_i\). The last interval, for \(w\), has length
\(\asymp n/(d_jd_k)\) independent of \(\mathbf v'\), with endpoints
depending on \(\mathbf v'\). Each of the \(j\) progression counts
has an \(O(1)\) error, giving the displayed bound with \(N\) the product
of interval lengths.

Sieve \(\prod b_i\) over \(2k<q<z\). Since the \(d_i\) are pairwise
coprime, at least two coefficients in \(\sum d_i b_i=n\bmod q\)
are units when \(k\ge3\); each coordinate then ranges uniformly.
For \(k=2\) and \(q\mid d_i\), the other coordinate is forced to a
nonzero value, since \(q\nmid n\), and \(b_i\) ranges freely.
In all cases the root proportion is at most \(k/q\). By the CRT
and the box estimate, the remainder at a squarefree \(r\le n^{1/2}\)
is \(O_k(n^{j-1}k^{\omega(r)})\).
Mertens gives the dimension-\(k\) product bound and
\(V(z)\gg_k(\log z)^{-k}\). The lower-bound fundamental lemma
leaves \(\gg_k N/(\log n)^k\) tuples, uniformly over multipliers.
The summed remainder is \(O_k(n^{j-1+1/2}(\log n)^{k-1})\),
which is negligible because \(P\le n^{k\eta}\).
[Cao–Huang–Zhang, Theorem 3.2](sources/cao-huang-zhang-2509.23320v1.pdf)
[Rosser–Schoenfeld, Theorems 7–8](sources/rosser_schoenfeld_1962.pdf)

Keep just the tuples with squarefree odd product of the \(b_i\).
All offending primes are \(q\ge z\), hence do not divide \(H P\).
For \(q^2\mid b_i\) the modulus \(q^2\) has root proportion \(1/q^2\),
and total error \(O_k(n^{j-1})\) by the box estimate. Summing to
\(\sqrt n\) gives \(O_k(N/z+n^{j-1/2})\).
For \(k\ge3\), two zero coordinates modulo \(q\) have proportion
\(1/q^2\) and error \(O_k(n^{j-1}/q)\). Summing to \(n\) gives
\(O_k(N/z+n^{j-1}\log n)\).
For \(k=2\), a shared \(q\) divides \(n\), so only \(O_\theta(1)\)
can occur, with count \(O(N/q+1)\) per prime. All these errors are
\(o(N/(\log n)^k)\).

## Conditional representation counts and overlap

For each remaining tuple \(a\), write \(R_a\) for the odd product of
the \(b_i\). It is squarefree, has at most \(J=\lceil k/\theta\rceil\)
prime factors, and \(R_a\asymp_k n^k/P\).
The conditional success probability is
\(w_a=\lambda^{\omega(R_a)}R_a^{-\gamma}
       \asymp n^{-k\gamma}P^\gamma\).
Each ordered tuple has a unique factorization into the specified
\(d_i\) and \(b_i\), by the separation \(U<z\) of the odd prime
ranges and the fact that the \(d_i\) are odd.
Using the tuple count and (5), the conditional mean of the number
of successful tuples is at least a constant times
\[
 n^u L_n,\qquad
 L_n=\frac{(\log n)^{k(\lambda-1)}}{(\log\log n)^{k\lambda}}.
 \tag{6}
\]

Discard any tuple for which \(\sum_i c_i b_i=0\) for some
nonzero integer vector \((c_i)\) with \(\max_i|c_i|\le n^\eta\).
For fixed multipliers and \((c_i)\), the equation
\(\sum_i d_i b_i=n\) and the homogeneous relation either are
inconsistent (if the vectors are proportional), or have at most
\(O_k(n^{k-2})\) solutions by fixing all but two coordinates.
Over all choices this costs \(O_k(n^{k-2+2k\eta})\) tuples and
conditional mean \(O(n^{u-1+(2k+k\gamma)\eta})=o(n^uL_n)\).
For a fixed \(R_a\) there are at most \(k^J\) assignments of its
primes to the \(b_i\), with the parity choices fixing the factor of 2.
For each assignment, at most one remaining multiplier vector solves
the sum equation, since two would give a discarded homogeneous relation.
Thus a full common support contributes \(O(\mu)\) to the ordered
dependency sum, where \(\mu\) is the surviving conditional mean.

For a squarefree divisor \(d\) of some \(R_a\), write \(D(d)\)
for the number of surviving tuples for which \(d\mid R_a\).
Give each prime of \(d\) to one of the \(k\) coordinates \(a_i\),
producing \(e_i\mid a_i\) with \(\prod e_i=d\). There are at most
\(k^J\) assignments. Fix all coordinates except the two with
smallest \(e_i,e_\ell\), and use the sum and the CRT in those two.
They have \(O(n/(e_ie_\ell)+1)\) choices, since \((e_i,e_\ell)=1\).
Their product is at most \(d^{2/k}\). Including the other
coordinates proves
\[
                D(d)\ll_{k,J} \frac{n^{k-1}}{d}
                                     \left(1+\frac{d^{2/k}}n\right).
\]
Whenever the gcd \(d>1\) of two supports is proper in at least
one of them, \(z\le d\le n^k/z\), since a missing prime is at least
\(z\). On this whole interval, using \(\gamma\le1-1/k\) gives
\[
 d^\gamma D(d)
 \ll n^{k-1}\big(d^{-1/k}+d^{1/k}/n\big)
 \ll n^{k-1}z^{-1/k}.
\]
For a fixed overlapping pair the conditional joint success
probability is \(w_aw_{a'}d^\gamma\lambda^{-\omega(d)}\).
The last \(\lambda\)-factor is bounded in terms of \(\lambda,J\).
For the proper overlaps, sum over all divisors in \([z,n^k/z]\)
of each \(R_a\), of which there are at most \(2^J\), and bound the
other success weight by \(O(n^{-k\gamma+k\eta\gamma})\).
Together with the full-support estimate this gives
\[
        \Delta\ll \mu\big(1+n^{u+k\eta\gamma}z^{-1/k}\big).
\]
Thus Janson's inequality bounds the probability that the candidate
success count is below \(\mu/2\) by
\[
 \exp\!\left[-c L_n\min(n^u,n^{\theta/k-k\eta\gamma})\right].
 \tag{7}
\]
The number of unordered representations is at least the number of
successful ordered tuples divided by \(k!\).
[Janson–Warnke, introductory equation (2)](sources/janson-warnke-1406.1248v1.pdf)

For \(u>0\), (7) plus the exceptional probability in (5) is summable
over \(n\), proving (3). For \(u=0\), \(\lambda>1\), the exceptional
probability is smaller than the bound (7), giving (2). It is
summable over \(n\) if \(k(\lambda-1)>1\).
If \(k(\lambda-1)>0\), it always sums over dyadic blocks after
normalization by the block length. Indeed the expected gap fraction
in \((2^{j-1},2^j]\) is \(O(\exp(-c j^{k(\lambda-1)}/
(\log j)^{k\lambda}))\). Markov and Borel–Cantelli show almost-sure
density one. All these representations use only 2 and odd primes
above \(Y\).

For the logarithmic enhancement with \(\beta>0\), for each fixed
\(\lambda\) its prime probabilities dominate the multiplier model
at primes \(q>Y\) for all sufficiently large \(n\).
Couple the selections by the same independent uniforms at the primes.
Given any fixed \(L\), choose \(k(\lambda-1)>L\) in (2). Since
the construction uses only primes above \(Y\) and 2, the logarithmic
enhancement inherits the stated tail. For \(\beta<0\),
\(\sum_q p_q q^{-\alpha}<\infty\), by partial summation of the
prime reciprocal sum. Therefore \(\sum_{a\in C}a^{-\alpha}<\infty\)
almost surely by the expected Euler product, and
\(N_C(x)=o(x^\alpha)\) by bounding the contribution of its tail.
It follows that \(|kC\cap[1,x]|=o(x)\). The case \(\beta=0\)
is the unmodified multiplier model \(\lambda=1\), treated below.
In the multiplier model, (1) similarly gives density zero for
\(\lambda<1\).

## Gaps below the logarithmic transition

Fix critical \(\gamma=1-1/k\) and \(1\le\lambda<1+1/k\). Put
\(a=k(\lambda-1)\in[0,1)\) and \(\kappa_k=\Gamma(1+1/k)^k/k!\).
For each integer \(l\ge1\) the counting moment bound is
\[
        \mathbb E[N_C(x)^l]\ll_l x^{l\alpha}(1+\log x)^{l(\lambda-1)}.
 \tag{8}
\]
It also holds, uniformly in \(q\ge3\), with any one prime \(q\)
forced to be selected. The local series for the expected simultaneous
membership indicators in \(l\) coordinates is
\(1+p_q(\prod_i(1-t_i)^{-1}-1)\).
Multiply by \(\prod_i(1-q^{-\gamma}t_i)^\lambda\).
Its linear coefficients cancel away from the finitely many capped
or deterministic terms. The weighted absolute coefficient norm of
the nonconstant part at \(t_i=q^{-\alpha}\) is
\(O(q^{-\gamma-2\alpha})\), so the product of norms is finite.
Forcing one prime changes one factor to a norm bounded uniformly
over the primes. The factors removed in each coordinate have
nonnegative coefficients \(\tau_\lambda(c)c^{-\gamma}\), defined
by \(\zeta(s+\gamma)^\lambda\). By Delange their partial sums are
\(O(x^\alpha(1+\log x)^{\lambda-1})\). Multivariate convolution
gives (8), since \(\lambda\ge1\) and the residual coefficient
function has finite weighted absolute norm.
[Delange, Theorem III](sources/delange-ikehara-1954.pdf)

Take \(a<b<1\) and \(Y=\exp((\log x)^b)\).
Let \(\mathcal F_Y\) describe the choices up to \(Y\), and \(C_Y\)
the semigroup using only those selected primes and 2. Let \(I_x\)
be either the integers in \([1,x]\) or those in \((x,2x]\), for
integer \(x\), and put
\[
 Z_x=x^{-1}|I_x\setminus kC|,\qquad
 V_x=\frac{\sum_{n\in I_x}r_k(n)}{x(\log x)^a}.
\]
Along each fixed geometric sequence \(x_j=\lfloor v^j\rfloor\) with
\(v>1\), almost surely \(N_{C_Y}(2x)=o(x^\alpha)\).
Indeed put \(s=\alpha-1/\log Y\) and use the expected Euler product
over primes through \(Y\). Its odd-prime factors differ from 1 by
\(p_q/(q^s-1)\ll q^{-1+1/\log Y}\ll1/q\), and its dyadic factor
is bounded. Its value is \(O((\log Y)^C)\), so
\(\mathbb E N_{C_Y}(2x)\ll(2x)^s(\log Y)^C\). Markov bounds the
probability that \(N_{C_Y}(2x)\) exceeds
\(x^\alpha\exp(-(\log x)^{1-b}/2)\) by
\(O(\exp(-c(\log x)^{1-b}))\). This is summable along the sequence.
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

For an unexposed prime \(q\le 2x\), fix the other choices and let
\(C^+\) be the semigroup with \(q\) included. Any new representation
after switching on that prime has an element of \(C^+\) divisible
by \(q\), whose count through \(2x\) equals \(N_{C^+}(2x/q)\).
Thus the change in \(Z_x\) is bounded by
\(k N_{C^+}(2x/q)N_{C^+}(2x)^{k-1}/x\), and the change in
\(V_x\) has the same bound divided by \((\log x)^a\).
By Hölder and the forced-prime version of (8) with \(l=2k\), the
mean squared changes are \(O(q^{-2\alpha}(\log x)^{2a})\) and
\(O(q^{-2\alpha})\), respectively. Reveal the independent bits
above \(Y\) one at a time. Orthogonality of the differences of
successive conditional expectations, followed by Jensen on the
future bits, bounds each conditional variance by the sum of the
conditional mean squared changes times \(p_q(1-p_q)\). Averaging
over the exposure gives
\[
 \mathbb E\operatorname{Var}(Z_x\mid\mathcal F_Y)
 +\mathbb E\operatorname{Var}(V_x\mid\mathcal F_Y)
 \ll (\log x)^{2a}\sum_{q>Y} p_q q^{-2\alpha}
 \ll (\log x)^{2a}\exp(-\alpha(\log x)^b),
 \tag{9}
\]
since \(\gamma+2\alpha=1+\alpha\). Chebyshev and Borel–Cantelli
therefore give almost surely, along any such geometric sequence,
\[
 |Z_x-\mathbb E[Z_x\mid \mathcal F_Y]|
 +|V_x-\mathbb E[V_x\mid \mathcal F_Y]|
      =O_\omega\!\left(\exp(-(\alpha/4)(\log x)^b)\right).
 \tag{10}
\]

For either interval \(I_x\), (4) with \(m=1\) gives
\(V_x\to\kappa_k A_{\omega,\lambda}^k\), so (10) gives the same
limit for its conditional mean. In that exposure discard the
targets already represented from \(C_Y\), at most \(N_{C_Y}(2x)^k=o(x)\)
along the sequence. Every other representation has conditional
success probability at most \(\lambda Y^{-\gamma}\), or zero if
it uses a known excluded prime. Harris and
\(-\log(1-w)\le w/(1-\lambda Y^{-\gamma})\), followed by Jensen over
the retained targets, give
\[
 \mathbb E[Z_x\mid \mathcal F_Y]\ge
 (1-o(1))\exp\!\left\{-
 \frac{(\log x)^a\,\mathbb E[V_x\mid\mathcal F_Y]}
      {(1-o(1))(1-\lambda Y^{-\gamma})}\right\}.
\]
Because \(b>a\), the error for \(Z_x\) in (10) is negligible
compared with this lower bound. Taking \(v=2\) and \(I_x=(x,2x]\)
gives the asserted dyadic bound and hence infinitely many gaps
for \(1<\lambda<1+1/k\).
At \(\lambda=1\), take \(I_x=[1,x]\) and a countable sequence
of ratios \(v>1\) decreasing to 1. Monotonicity of the unnormalized
gap count between successive grid points loses a factor at most
\(v\) in the limiting lower bound. Letting \(v\downarrow1\)
gives the asserted lower density.
[Janson–Warnke, introductory equation (3)](sources/janson-warnke-1406.1248v1.pdf)

## Positive probability at the exact multiplier

Now take \(\gamma=1-1/k\) and \(\lambda=1+1/k\).
Use the same fixed \(\theta,\eta,\tau\) of the upper-gap proof,
but take a fixed lower cutoff \(Y_0>2k\) large enough that
\(p_q=\lambda q^{-\gamma}\) for \(q>Y_0\). Call a target \(n\)
suitable if
\[
             \sum_{\substack{Y_0<q\le n^\tau\\q\mid n}}
                         X_q q^{-\alpha}\le1 .
\]
This condition is determined by the exposure up to \(z=n^\theta\).
The full pool product without the exclusion \(q\mid n\) satisfies
\[
      \prod_{Y_0<q\le U}(1+k q^{-\alpha})^{X_q}
         \sim D_\omega(\log U)^{k\lambda}
\]
almost surely, with \(0<D_\omega<\infty\).
Indeed \(\sum_{q>Y_0}(X_q-p_q)\log(1+kq^{-\alpha})\) converges
almost surely by the finite sum of variances. Subtracting
\(k\lambda/q\) from its deterministic log-factor means leaves
an absolutely convergent sum, and
\(\sum_{q\le U}1/q-\log\log U\) has a finite limit.
The variable \(D_\omega\) has unbounded upper support. Selecting
every prime from \(Y_0\) through a fixed large \(R\) has positive
probability, and on this event the corresponding finite sum of
centered log factors grows without bound as \(R\) increases.
The tail centered factors are independent of it and have
variance tending to zero as \(R\) increases.
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

Excluding primes dividing a suitable target costs a factor at
most \(e^k\) in this product, by \(\log(1+kq^{-\alpha})\le kq^{-\alpha}\).
The size truncation \(d_i\le n^\eta\) still retains at least half
the product weight eventually almost surely. For this use
\(\sum_{Y_0<q\le U}X_q q^{-\alpha}\log q=O_\lambda(\log U)+O_\omega(1)\),
since the centered series with the extra \(\log q\) also has summable
variances. Take \(\tau\) small as before. Thus for any prescribed
large constant \(K\) there is a positive-probability event on which
the permitted multiplier weight through \(U=n^\tau\) exceeds
\(K(\log n)^{k\lambda}\) for every sufficiently large suitable \(n\).

Given an exposure where that weight bound holds for a suitable
\(n\), the tuple count and overlap argument give
\(\mu\gg K\log n\) and \(\Delta\ll\mu\).
The constants here are independent of \(K\). Choose \(K\) so the
conditional failure bound is \(O(n^{-2})\). Borel–Cantelli shows
that on the positive-probability event above, every sufficiently
large suitable \(n\) is represented.
There is therefore a fixed integer \(L>Y_0\) such that the event
that all suitable \(n>L\) are represented has positive probability.
This event is increasing in the prime selections. Adding primes
can only increase the set of represented targets and shrink the
set of suitable targets. Hence it still has positive probability
after requiring all primes through \(L\) to be selected.

On that event, start with \(n>L\) and divide out one selected bank
prime whenever the current target exceeds \(L\) and is unsuitable.
Such a prime exists by definition. The procedure either stops at a
represented suitable target or at an integer at most \(L\).
In the first case, multiply the representation by the removed
factors. In the second case, the final integer and all removed
primes belong to \(C\). At least one removed prime \(q>Y_0\) can
itself be represented. If \(q>L\), it is suitable as
\(q^\tau<q\). If \(q\le L\), use \(k-1\) copies of 1 and \(q-k+1\),
which is a positive integer at most \(L\). Multiplying that
representation by \(n/q\in C\) covers \(n\) in the second case.
This proves cofiniteness on an event of positive probability.

## Small constants and gaps at the exact multiplier

**Theorem.** Fix \(k\ge2\), \(\alpha=1/k\), \(\gamma=1-\alpha\), and
\(\lambda=1+1/k\) in the constant-multiplier model. The constant
\(A=A_{\omega,\lambda}\) in (1) has positive probability of lying
below any positive number and above any finite number.
For each \(\varepsilon>0\) there is a fixed finite prime pattern
\(\mathcal E\) of positive probability such that, for all sufficiently
large integer \(x\),
\[
 \mathbb E\big[ |(x,2x]\setminus kC|\mid\mathcal E\big]
 \gg_\varepsilon x^{1-\varepsilon}.
 \tag{11}
\]
Consequently, for the extended largest gap
\(M=\sup((\mathbb N\setminus kC)\cup\{0\})\),
\[
              \Pr(M>x)\gg_\varepsilon x^{-\varepsilon}
                        \qquad(\varepsilon>0).
 \tag{12}
\]
Every positive moment of \(M\) is infinite.

At \(s=\alpha\), the factor \(G_\omega\) in (1) is
\[
            G_\omega(\alpha)
     = \prod_q (1-1/q)^\lambda(1-q^{-\alpha})^{-X_q},
\]
including the mandatory choice \(X_2=1\). For \(y\) beyond the
capped primes, write \(R_y\) for the logarithm of the product of
the factors with \(q>y\). Expanding the logarithms and using
\(p_q=\lambda q^{-\gamma}\) gives
\[
 |\mathbb E R_y|+\operatorname{Var}(R_y)\ll y^{-\alpha},
 \qquad
 \mathbb E e^{tR_y}=\exp(O_t(y^{-\alpha}))\quad(t\in\mathbb R\text{ fixed}).
 \tag{13}
\]
For the exponential moment, the individual factor is
\[
 (1-1/q)^{t\lambda}
       \big(1+p_q[(1-q^{-\alpha})^{-t}-1]\big)
         =1+O_t(q^{-1-\alpha}).
\]
The partial products of the expectations at \(2t\) are bounded,
so uniform integrability justifies passing to the infinite product.
In particular the tail factor tends to 1 in probability and in
each fixed positive moment.

If every optional prime through \(y\) is excluded, the first
part of the product for \(G_\omega(\alpha)\) is
\(\asymp(\log y)^{-\lambda}\) by Mertens, with a constant from the
finitely many mandatory primes. If every prime through \(y\)
is included, the logarithm of this initial factor tends to
infinity, since the positive terms are \(\gg q^{-\alpha}\) and
\(\sum_q q^{-\alpha}=\infty\), whereas the negative term at each
large prime is \(O(1/q)\). Both finite patterns have positive
probability and are independent of \(R_y\). This proves the support
assertion and, for the all-optional-excluded pattern \(\mathcal E_y\),
\[
                   \mathbb E[A^k\mid\mathcal E_y]
                     \ll (\log y)^{-k\lambda}.
 \tag{14}
\]
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

The proof of (8) applies for every \(\lambda\ge1\), including this
endpoint, since it uses only the convergence of the multivariate
Euler factors and the nonnegative univariate coefficients. Higher
choices of \(l\) give uniform integrability of every fixed power
of \(N_C(x)/(x^\alpha(\log x)^{\lambda-1})\).
Thus (1) holds in all finite \(L^r\). Subtracting the two cumulative counts in (4) with \(m=1\), and using
the same higher moment bounds, gives
\[
          \frac{\sum_{x<n\le2x}r_k(n)}{x\log x}
             \longrightarrow \kappa_k A^k
             \quad\text{almost surely and in each finite }L^r.
 \tag{15}
\]
Here the block counts use the unordered representation convention
of (4). The passage to \(L^r\) follows from the bound of the
numerator by \(N_C(2x)^k\). These limits also hold conditional on any
fixed prime pattern of positive probability, as the conditional
expectations of nonnegative errors are at most their unconditional
values divided by the probability of the pattern.

Condition on \(\mathcal E_y\), with \(y\) fixed. The finitely many
mandatory primes can supply representations for only
\(O_y((1+\log x)^{O_y(1)})\) targets through \(2x\).
Every representation of a remaining target either is impossible
under \(\mathcal E_y\) or has conditional success probability
at most \(\lambda y^{-\gamma}\). Harris bounds its target's failure
probability below by the exponential of minus the conditional
mean representation count divided by \(1-\lambda y^{-\gamma}\).
Jensen over the \(x-o(x)\) retained targets and (15) show that
\[
 \mathbb E\big[ |(x,2x]\setminus kC|\mid\mathcal E_y\big]\ \ge\
 (1-o(1))x
 \exp\!\left[-\left(
 \frac{\kappa_k\mathbb E[A^k\mid\mathcal E_y]}{1-\lambda y^{-\gamma}}
          +o(1)\right)\log x\right].
\]
Use (14) and choose \(y\) so the displayed coefficient is less
than \(\varepsilon\). This proves (11).
The count in (11) is at most \(x\,\mathbf 1_{\{M>x\}}\) for
integer \(x\), giving (12) after multiplying by
\(\Pr(\mathcal E_y)>0\). For each \(r>0\), choose
\(\varepsilon<r\) in \( \mathbb E[M^r]\ge x^r\Pr(M>x)\).
[Janson–Warnke, introductory equation (3)](sources/janson-warnke-1406.1248v1.pdf)

## Conditioning at subpower and polynomial cutoffs

At equality, use the notation \(Z_x,V_x\) of (9) with \(a=1\) and
\(I_x=(x,2x]\). For \(Y=\exp((\log x)^b)\), \(0<b<1\),
the proofs of the Rankin bound and (9)–(10) still apply.
Thus at dyadic \(x\), the conditional mean of \(V_x\) tends
almost surely to \(\kappa_k A^k\), and only \(o(x)\) targets
are covered using exposed primes alone. The conditional
Harris–Jensen lower bound has the form
\[
                \mathbb E[Z_x\mid\mathcal F_Y]
                          \ge x^{-\kappa_k A^k-o(1)} .
 \tag{16}
\]
For every fixed \(A>0\), the additive concentration error in (10)
is larger in scale than the lower bound (16). A polynomial cutoff
\(Y=x^\theta\), \(0<\theta<1\), gives the stronger conditional
variance bound \(O((\log x)^2x^{-\alpha\theta})\) after averaging
over the exposure. The following count determines the scale of
representations already fixed at that cutoff.

**Theorem.** At the same endpoint, fix \(0<\theta<1\). Let
\(C_{x^\theta}\) be generated by the selected primes through
\(x^\theta\) and by \(2\), and let \(r_{k,x^\theta}\) count its
unordered \(k\)-term representations. There is a continuous
positive function \(W_\lambda\) on \(0<\theta\le1\) such that
almost surely
\[
 \begin{aligned}
 |C_{x^\theta}\cap[1,x]|
   &\sim A W_\lambda(\theta)x^\alpha(\log x)^{\lambda-1},\\
 \sum_{x<n\le2x}r_{k,x^\theta}(n)
   &\sim \kappa_k A^k W_\lambda(\theta)^k x\log x .
 \end{aligned}                                                    \tag{17}
\]
In particular the polynomial exposure already determines
representation tuples at order \(x\log x\).

By inclusion–exclusion over the selected primes exceeding
\(x^\theta\),
\[
 |C_{x^\theta}\cap[1,x]|=
 \sum_{\substack{d\le x\ {\rm squarefree}\\
                 q\mid d\ \Rightarrow\ q>x^\theta,\ X_q=1}}
              \mu(d)N_C(x/d).
\]
Each \(d\) uses at most \(\lfloor1/\theta\rfloor\) primes.
The measures
\(\sum_{x^\theta<q\le x}X_q q^{-\alpha}\delta_{\log q/\log x}\)
on \([\theta,1]\) converge weakly almost surely to \(\lambda\,dv/v\).
For each fixed interval this follows from convergence of
\(\sum_q(X_q-p_q)q^{-\alpha}\), which has summable variances,
and Mertens' product estimate for \(\sum 1/q\).
For each fixed number \(j\) of prime variables, repeat primes
contribute \(o(1)\) to the product measure since the largest atom is
at most \(x^{-\alpha\theta}\) and the total masses are bounded.
For \(d\le x^{1-\delta}\), (1) applies to \(x/d\) with uniform relative
error tending to zero for any fixed \(\delta>0\).
For \(d>x^{1-\delta}\) use the uniform bound
\(N_C(t)\ll_\omega t^\alpha(1+\log t)^{\lambda-1}\) for \(t\ge1\).
The normalized summands are then bounded by
\(O_\omega(d^{-\alpha}(\delta+1/\log x)^{\lambda-1})\).
Since \(\lambda>1\), letting \(\delta\downarrow0\) controls this
region by the bounded product-measure masses. This gives the first
limit of (17) with
\[
 W_\lambda(\theta)=
 \sum_{j=0}^{\lfloor1/\theta\rfloor}
 \frac{(-\lambda)^j}{j!}
 \int_{\substack{v_i\ge\theta\\ \sum_{i=1}^j v_i\le1}}
            (1-\textstyle\sum_i v_i)^{\lambda-1}
                           \prod_{i=1}^j\frac{dv_i}{v_i},
 \tag{18}
\]
where the term for \(j=0\) equals 1.
The finite integral formula is continuous in \(\theta>0\), since
only finitely many \(j\) are relevant on a neighborhood of each
positive \(\theta\), and the variable boundaries have measure zero.
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

To see positivity, set \(R(u)=u^{\lambda-1}W_\lambda(1/u)\) for
\(u\ge1\), set \(R(u)=u^{\lambda-1}\) for \(0<u<1\), and extend
by zero on negative arguments. Substituting \(t_i=uv_i\) in (18)
gives
\[
 \widehat R(s)=\int_0^\infty e^{-su}R(u)\,du
    =\Gamma(\lambda)s^{-\lambda}
      \exp\!\left(-\lambda\int_1^\infty e^{-st}\frac{dt}{t}\right)
                          \qquad(s>0).
\]
The interchange follows by integrating the absolute series,
whose transform has the same expression with a positive
sign in the exponential. Differentiating gives
\(-\widehat R'(s)=\widehat R(s)\lambda(1-e^{-s})/s\),
so uniqueness of Laplace transforms yields
\[
                   uR(u)=\lambda\int_{u-1}^u R(t)\,dt .
\]
The continuous function \(R\) is positive on \((0,1]\).
A first positive zero would have a strictly positive integral
on its preceding unit interval, a contradiction.

For each fixed \(t>0\), trap \(x^\theta\) between
\((tx)^{\theta-\delta}\) and \((tx)^{\theta+\delta}\) and use the
first limit and continuity, taking a countable family of exponents
and rational \(t\) for the weak convergence argument.
The normalized counting measure of \(a/x\) for
\(a\in C_{x^\theta}\) therefore has limit
\(A W_\lambda(\theta)\alpha t^{\alpha-1}\,dt\) on \([0,2]\).
Apply the simplex calculation of (4) to the two total counts
through \(x\) and \(2x\), and subtract. Tuples with a repeated
coordinate are negligible by the same bound as in (4).
Since \(k\alpha=1\) and \(k(\lambda-1)=1\), this gives the second
limit of (17).

## Stopped martingales at the exact multiplier

**Theorem.** At \(\gamma=1-1/k\), \(\lambda=1+1/k\), and \(k\ge2\),
the probability of cofiniteness is strictly between zero and one.
For the finite pattern \(\mathcal E_y\) excluding all optional primes
through \(y\), one has, as \(y\to\infty\),
\[
 \Pr(kC\text{ cofinite}\mid\mathcal E_y)
       \ll_k (\log y)^{-(k+1)} .
\]
For any fixed \(0<t<1/(2k)\), the conditional probability of at least
\(c_t x^{1-t}\) gaps in infinitely many dyadic blocks \((x,2x]\)
tends to one as \(y\to\infty\).

### A uniform moment bound for restricted prime factors

Write \(N_Q(X)=|C_Q\cap[1,X]|\), where \(C_Q\) uses the selected
primes through \(Q\) and the deterministic prime 2. For \(X,Q\ge2\),
\[
 \mathbb E[N_Q(X)^k]\ll_k X\log(2X)
                  \exp\!\left(-c_k\,\frac{\log X}{\log Q}\right).
 \tag{19}
\]
To prove this, take \(\eta=\min(\alpha/4,1/\log Q)\) and
\(\sigma=\alpha-\eta\), where \(\alpha=1/k\).
In the multivariate moment series, the local coefficient of every
nonzero vector of prime exponents at \(q\le Q\) is the selection
probability \(p_q\), with \(p_2=1\). At \(q>Q\) the factor is 1.
Use the coefficientwise majorant
\[
                    \prod_{i=1}^k\zeta(s_i+\gamma+\eta)^\lambda
                                              H(\mathbf s).
\]
The local factor of \(H\) has nonnegative coefficients and constant
term 1. At \(q\le Q\), give each linear term coefficient
\((p_q-\lambda q^{-\gamma-\eta})_+\) and each term of total degree
at least two coefficient \(p_q\). Take \(H_q=1\) for \(q>Q\).
Multiplication by the local zeta factors supplies the rest of the
linear coefficient, while their constant term preserves every
coefficient of \(H_q\).

The total mass of the linear corrections at \(\mathbf s=(\sigma,\ldots,\sigma)\)
over odd primes is
\[
 \ll_k \sum_{q\le Q}\frac{q^\eta-1}{q}
 \ll_k \eta \sum_{q\le Q}\frac{\log q}{q}\ll_k 1.
\]
Here \(\eta\log Q\le1\). The terms of degree at least two contribute
\(O_k(q^{-\gamma-2\sigma})=O_k(q^{-1-\alpha/2})\) at each odd \(q\),
and the local factor at 2 is bounded since \(\sigma\ge3\alpha/4\).
Thus the coefficient mass of \(H\), weighted by
\(\prod_i d_i^{-\sigma}\), is uniformly bounded.
The zeta factor in each coordinate has summatory coefficients
\[
 \sum_{b\le v}\frac{\tau_\lambda(b)}{b^{\gamma+\eta}}
            \ll_k v^{\alpha-\eta}(1+\log v)^{\lambda-1}.
\]
This follows from the unshifted estimate used in (8) by partial
summation, uniformly because \(\alpha-\eta\ge3\alpha/4\).
Convolving the majorant and using \(k\alpha=1\) and
\(k(\lambda-1)=1\) gives the claimed bound, since
\(k\eta\gg_k1/\log Q\) for \(Q\ge2\).
[Delange, Theorem III](sources/delange-ikehara-1954.pdf)
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

### Compensating a surviving target

Take integer \(x\), put \(X=2x\), and work under the conditional
law \(\mathbb P_y=\mathbb P(\,\cdot\mid\mathcal E_y)\), for \(y\)
beyond the mandatory primes and large enough that \(p_q\le1/2\)
for \(q>y\). The remaining bits are independent.
Expose the primes \(y<q\le X\) in increasing order. For each
target \(n\in(x,2x]\) define \(D_n(q)\) to be the event that one
of its representations uses \(q\) and otherwise uses previously
selected primes. This means \(q\) divides at least one of the
summands, and every other prime divisor of the summands is
already selected. This predicate is increasing in the preceding
bits and is defined even for a target with an earlier representation.

Each initially unrepresented target starts at weight \(W_n=1\).
The other targets start at zero. At a dangerous prime
(\(D_n(q)\) holds), for a positive unfrozen weight \(w\), first
test whether \(w/(1-p_q)>q^t\). If so, freeze \(w\) permanently,
regardless of the upcoming bit at \(q\).
Otherwise replace \(w\) by zero if \(q\) is selected and by
\(w/(1-p_q)\) if it is excluded. All other weights stay fixed.
The freeze decision is measurable before revealing \(q\), and
each update has conditional mean \(w\).
Thus the weights are martingales. A positive terminal weight
on an unfrozen target certifies a gap, since a first
representation would have triggered a lethal update at its
largest prime divisor.

Let \(Z_x\) be the normalized terminal total weight
\(x^{-1}\sum_n W_n\) and \(S_x\) its normalized frozen part.
The initial total satisfies \(\mathbb E_y Z_x=1-o(1)\), since
only polylogarithmically many targets have a representation using
the finitely many mandatory primes. All terminal unfrozen weights
are bounded by \(X^t\).

### Variance and the change of measure

At \(q\), all active dangerous targets that pass the freeze
test are updated by the same Bernoulli bit. If the sum of their
preceding weights divided by \(x\) is \(B_q\), the increment of
the normalized total has conditional variance
\(p_q B_q^2/(1-p_q)\). These preceding weights are at most \(q^t\).
Let \(T_q\) count ordered tuples with sum in \((x,2x]\) using
previously selected primes and the provisionally selected prime
\(q\), with at least one coordinate divisible by \(q\).
Then \(B_q\le q^t T_q/x\).
For \(F_q(\alpha)=\prod_{\ell\le q}(1-\ell^{-\alpha})^{-X_\ell}\)
with \(q\) forced in, the count of eligible multiples of \(q\)
through \(X\) is the semigroup count through \(X/q\).
Rankin's bound therefore gives
\(T_q\le kXq^{-\alpha}F_q(\alpha)^k\).
For its moment, the odd-prime Euler factors away from the
forced prime are
\[
 1+p_\ell\big((1-\ell^{-\alpha})^{-2k}-1\big)
      =1+2k\lambda/\ell+O_k(\ell^{-1-\alpha})
\]
at all large \(\ell\) in the unconditioned model.
The forced factor is uniformly bounded, and exclusion through \(y\)
can only decrease the product. Thus
\(\mathbb E_y F_q(\alpha)^{2k}\ll_k(\log q)^{2k\lambda}\).
Orthogonality of increments now gives, uniformly in \(x\),
\[
 \operatorname{Var}_y Z_x
     \ll_k\sum_{q>y}
        q^{-1-\alpha+2t}(\log q)^{2k\lambda}.
 \tag{20}
\]
[Rosser–Schoenfeld, Theorems 5–8](sources/rosser_schoenfeld_1962.pdf)

For an initially unrepresented target \(n\), its terminal weight
is a density for a new probability measure, since its mean is 1.
Its conditional density at each exposure step is the current
weight, by the martingale property. Under the new measure a
dangerous update forces exclusion, and a step without an update
has its original probability \(p_q\). In particular this includes
the step at which the weight freezes and all steps thereafter.
Couple this tilted process with the original one by shared
independent uniforms at each prime. The tilted path selects a
subset of the original primes, because each of its conditional
selection probabilities is either \(0\) or \(p_q\).

Put
\[
 H_n(Q)=\sum_{y<q\le Q}\mathbf 1_{D_n(q)}[-\log(1-p_q)] .
\]
On a tilted path before freezing, the accumulated log-weight is
the sum over completed dangerous updates. Freezing at \(q\)
implies \(H_n(q)>t\log q\), since the hazard includes that
prospective update. The definition of \(D_n(q)\) makes \(H_n(Q)\)
increasing in the earlier prime selections, including after
freezing. Write \(\theta_j=2^{-j}\) and \(Q_j=X^{\theta_j}\).
If freezing occurs at a prime in \((Q_{j+1},Q_j]\), then
\(H_n(Q_j)>(t/2)\theta_j\log X\). Monotone coupling bounds the
probability of this event under the tilt by the probability
of that hazard inequality under the original conditional law.

Under the original conditional law \(D_n(q)\) is predictable
and \(q\) is selected independently of the past. Therefore
\(p_q\mathbb P_y(D_n(q))\) is at most the expected number
of representations of \(n\) whose largest prime divisor is \(q\).
As \(-\log(1-p_q)\le2p_q\), summing over \(n\) gives
\[
       \sum_{x<n\le2x}\mathbb E_y H_n(Q)
                                 \le 2\mathbb E_y N_Q(X)^k .
\]
Every representation in this count has just one largest prime
divisor. The expectation of an initially unrepresented
target's terminal weight times its freeze indicator equals its
freeze probability under its tilted measure.
Summing the hazard bounds by Markov's inequality therefore gives
\[
 \mathbb E_y S_x\ll_t
     \sum_{\substack{j\ge0\\Q_j>y}}
        \frac{\mathbb E_y N_{Q_j}(X)^k}{X\theta_j\log X}.
 \tag{21}
\]

### Passing from weight to gaps

The summands in (21) are dominated uniformly in \(x,y\) by
\(O_k(2^j e^{-c_k2^j})\), by (19) and monotone domination under
\(\mathcal E_y\). For fixed \(j,y\), the polynomial-cutoff asymptotic (17)
and uniform integrability give
\[
 \frac{\mathbb E_y N_{Q_j}(X)^k}{X\log X}
       \longrightarrow W_\lambda(\theta_j)^k
                          \mathbb E[A^k\mid\mathcal E_y].
\]
Uniform integrability follows by comparison with \(N_C(X)^k\)
and (8) with \(l=2k\). At \(j=0\), \(Q_0=X\) and
\(W_\lambda(1)=1\). Applying the same moment limit unconditionally
in (19) gives
\[
 W_\lambda(\theta_j)^k \mathbb E[A^k]\ll_k e^{-c_k2^j},
 \qquad \mathbb E[A^k]>0 .
\]
Thus dominated convergence in (21) bounds
\(\limsup_{x\to\infty}\mathbb E_y S_x\) by a constant times
\(\mathbb E[A^k\mid\mathcal E_y]\sum_j2^j W_\lambda(\theta_j)^k\).
The sum is finite, and (14) gives
\(\limsup_{x\to\infty}\mathbb E_y S_x\ll_{k,t}(\log y)^{-(k+1)}\).

For fixed \(y\) and large \(x\), Chebyshev bounds
\(\mathbb P_y(Z_x<1/2)\) by a constant times (20).
Markov gives \(\mathbb P_y(S_x>1/4)\le4\mathbb E_y S_x\).
On the complementary event there are at least
\(x/(4X^t)\) genuine gaps in \((x,2x]\), since each positive
unfrozen weight is at most \(X^t\).
Taking the limsup event along dyadic \(x\) bounds its probability
below by the limsup of the separate block probabilities.
The variance tail in (20) tends to zero faster than the stated
logarithmic bound, for example with \(t=\alpha/4\).
This proves the quantitative conditional estimate and the
assertion about infinitely many blocks. Each event \(\mathcal E_y\)
has positive probability, so infinitely many gaps occur with positive
probability unconditionally. The preceding positive-probability
cofiniteness theorem gives \(0<\Pr(kC\text{ cofinite})<1\).
Kolmogorov's zero-one law then excludes tail measurability, also
up to events of probability zero.

## Finite changes at the critical multiplier

**Corollary.** At \(\gamma=1-1/k,\lambda=1+1/k\), almost every
realization admits a finite addition of optional primes producing
a cofinite \(kC\), and a finite removal of optional primes producing
infinitely many gaps. Forcing all optional primes through \(y\) in
makes the conditional probability of cofiniteness tend to one.
Forcing them out makes it tend to zero.

Let \(B\) be the cofinite event. Its saturation under finite
optional changes is the union of countably many measurable
finite-coordinate translates of \(B\). For every finite prefix,
membership in the saturation can be tested with that prefix reset
to one fixed pattern, so the saturation is a tail event. The
same holds for \(B^c\). Both have positive probability by the
endpoint theorem, and Kolmogorov's zero-one law makes both
saturations have probability one. Monotonicity lets a successful
modification retain just its additions, and a modification
producing infinitely many gaps retain just its removals.

For almost every realization, the map forcing all optional
primes through \(y\) in therefore gives a cofinite sumset for
all sufficiently large \(y\), and the map forcing them all out
eventually gives a noncofinite sumset. Their laws are the
conditional laws on the two finite patterns, giving the claimed
limits. Almost surely these changes can switch the least basis
order between \(k\) and \(k+1\). The counting bound excludes
orders below \(k\), and the subcritical coverage argument gives
order \(k+1\) almost surely. These probability-one assertions
hold under all finite optional changes simultaneously since
each finite change is absolutely continuous with respect to the
original law and there are countably many of them.

The counting asymptotic of the selected primes is
\(\lambda x^\alpha/(\alpha\log x)\) almost surely, by partial
summation of the prime number theorem and the variance bound
by the mean on geometric grids. Finite changes alter that
count by \(O(1)\). For the generated semigroup they preserve
the exponents in (1) and multiply \(A\) by
\(\prod_{q\ \mathrm{added}}(1-q^{-\alpha})^{-1}
 \prod_{q\ \mathrm{removed}}(1-q^{-\alpha})\), by the finite Euler factors.

**Open question.** Determine whether \(kC\) has a natural density at
\(\gamma=1-1/k,\lambda=1\), and determine its law if it exists.
The complement has almost-sure lower density at least
\(\exp(-\kappa_k A_{\omega,1}^k)\).

**Open question.** Obtain a pointwise asymptotic or a fluctuation law
for \(r_k(n)\) in the subcritical model. The fixed-modulus averaged
law (4) measures congruence effects after averaging over targets.
