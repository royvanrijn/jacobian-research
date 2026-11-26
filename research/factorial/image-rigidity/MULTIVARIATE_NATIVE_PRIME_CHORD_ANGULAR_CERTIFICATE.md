# Native-prime angular certificates from upper endpoint blocks

Consider an ordinary two-monomial profile with fixed exponent
columns. A prime already present in the upper factorial array can
certify angular separation without modifying the profile or
computing the lower coefficients. The certificate below uses
only four sums of valuations of factorial blocks.
Its primes may vary between adjacent orders.

For an adjacent pair \(n,N=n+1\) with \(N\ge3\), take
\(p>N+1\). In a row of gap \(d>0\) and smaller exponent \(\ell\),
split the upper rising factors into the first block and remainder,
\[
 \prod_{h=1}^{d}(N\ell+h),\qquad
 \prod_{h=d+1}^{Nd}(N\ell+h).
\]
Sum their valuations over the increasing rows to obtain \(E_+,U\),
respectively, and over the decreasing rows to obtain \(E_-,V\).
Set \(T=E_+-E_-+U-V\). Reverse the two columns if needed
so \(T>0\), exchanging the signs and side labels accordingly.

**Theorem.** Suppose \(\gcd(T,N)=1\). Either of the inequalities
\[
                 T>NU,\qquad
                 NE_-+T>N(2U+V)                              \tag{1}
\]
ensures that the actual ordinary moments of orders \(n,N\)
have disjoint root arguments under independent positive
rescalings. Only the upper factorial blocks enter (1).
If a prime satisfying (1) exists at every \(3\le N<J\),
where \(J\) is at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
then the conclusion holds in every degree.

In any box of independently variable smaller exponents with
fixed nonzero gaps and orientations, and with any additional
fixed ordinary rows, the proportion of profiles satisfying
this all-degree arithmetic criterion tends to one. The
same is true within a fixed Cartesian product of arithmetic
progressions for the variable exponents. If the free
progression parameters range from \(1\) to \(H\), the
failure proportion is \(O(e^{-\sqrt{\log H}/4})\), with a
constant depending on the fixed data and the cutoff.
At least one variable row is required.

## A condition on the upper Newton polygon

Use \(v(p)=1\) for the valuation. Let \(a_k\) be the
coefficient valuations of the upper polynomial and set
\[
     T=a_N-a_0,\qquad
     \gamma_k=N a_k-(N-k)a_0-k a_N\quad(1\le k<N).
\]
Suppose \(\gcd(T,N)=1\) and all \(\gamma_k>0\).
Let \(\beta=\min_{1\le k<N}\gamma_k\). A sufficient condition is
that for some \(k\) coprime to \(N\),
\[
                              \gamma_k<2\beta .               \tag{2}
\]
This condition can hold even when the index attaining
\(\beta\) is not coprime to \(N\).

In the completion \(K_0\) of the maximal unramified
extension of \(\mathbb Q_p\), adjoin \(\pi=p^{1/N}\).
Substituting \(z=\pi^{-T}x\) and dividing by \(p^{a_0}\)
gives a polynomial of the form
\[
                         c_0+c_Nx^N+
                     \sum_{k=1}^{N-1} b_k\pi^{\gamma_k}x^k
\]
with \(p\)-adic unit coefficients \(c_0,c_N,b_k\) in \(K_0\).
The binomial has \(N\) unit roots \(\lambda_i\in K_0\)
with distinct residues. Hensel lifting gives \(N\) roots
\(x_i\) near them, and Taylor expansion gives
\[
       \frac{x_i}{\lambda_i}-1
         =-\frac{\sum_{k=1}^{N-1} b_k\lambda_i^k\pi^{\gamma_k}}
                     {N c_N\lambda_i^N}
                       +O(\pi^{2\beta}).                     \tag{3}
\]
Here the error has valuation at least \(2\beta/N\).
The same expansion holds for \(\log(x_i/\lambda_i)\).
In fact the \(j\)-th term of the logarithm has valuation
at least \((j\beta-Nv_p(j))/N\) for \(j\ge2\). This is
at least \(2\beta/N\), since \(j-2\ge(p-2)v_p(j)\),
\(p-2\ge N\), and \(\beta\ge1\).

For distinct upper roots put \(\omega=\lambda_i/\lambda_j\).
This is a nontrivial \(N\)-th root of unity. Subtracting
the two log expansions gives the logarithm of their quotient
divided by \(\omega\). The coefficient at exponent
\(\gamma_k\) in (3) is a unit times
\(\lambda_i^k-\lambda_j^k\), which is a unit for
\(k\) coprime to \(N\). The integers \(\gamma_k\) have
distinct residues modulo \(N\), by the hypothesis on \(T\).
Under (2) the corresponding coordinate of the logarithm
in the basis \(1,\pi,\ldots,\pi^{N-1}\) is nonzero,
with fractional denominator \(N\). The error cannot cancel
it: the basis terms have distinct fractional valuations,
so each error coordinate has valuation at least \(2\beta/N\).

Now suppose the same quotient were \(\alpha/\alpha'\)
for lower roots. Its valuation is zero, so the lower roots
have a common valuation \(s\). Over \(K_0\) they generate
totally tamely ramified extensions of degrees
\(e,e'\le n<p\), respectively. The
[tame basis calculation](MULTIVARIATE_SINGLE_MOVABLE_ROW_PADIC_PROGRESSION.md#the-lower-ratios-cannot-have-that-displacement)
puts them in \(K_0(p^{1/e})\) and \(K_0(p^{1/e'})\).
Since \(s\) belongs to both value groups, divide both roots
by the same power \(p^s\) from the intersection. Their
resulting units have residue ratio equal to the residue of
\(\omega\). Choose unramified unit lifts of their residues
with exact ratio \(\omega\). The logarithm of
\((\alpha/\alpha')/\omega\) is then the difference of the
principal-unit logarithms in the two tame subfields.

In the power basis of \(p^{1/\operatorname{lcm}(N,e,e')}\),
this difference has support only at fractions belonging
to \((1/e)\mathbb Z\) or \((1/e')\mathbb Z\). Neither
contains a fraction with reduced denominator \(N\).
This contradicts the upper coordinate from (3).
Consequently the two orders have disjoint quotients of
distinct roots. Equal nonreal complex arguments would
give equal quotients after dividing by the conjugate
roots. Equal real arguments are excluded by strict
positivity of the even moment on the real line.
This proves the polygon criterion.

## Checking the polygon from the factorial blocks

Strip the scalar \((N\ell)!\) from each row.
Because \(p>N\), the binomial coefficient at every index
is a unit. All increasing first blocks occur at each
\(k\ge1\), and all decreasing first blocks occur at each
\(k<N\). Write \(u_k,v_k\) for the accumulated remainder
valuations in the two orientations. Then \(u_k\) is
nondecreasing from \(0\) to \(U\), and \(v_k\) is
nonincreasing from \(V\) to \(0\). The endpoint difference
is \(T\) as above and the interior scaled heights are
\[
          \gamma_k=NE_-+(N-k)T+N(u_k-U+v_k).                  \tag{4}
\]
In particular, with \(T>0\), every \(\gamma_k\) is at least
\(NE_-+T-NU\), and \(\gamma_{N-1}\le NE_-+T+NV\).
The second inequality of (1) makes these lower and upper
bounds positive and gives \(\gamma_{N-1}<2\beta\).

For the first inequality of (1), \(\gamma_{N-1}\) is
positive. For \(1\le k<N-1\), monotonicity of the two
remainder sums gives
\[
 \gamma_k-\gamma_{N-1}
   =(N-1-k)T+N(u_k-u_{N-1}+v_k-v_{N-1})
   \ge (N-1-k)T-NU>0 .
\]
Thus \(\beta=\gamma_{N-1}\), and (2) follows because
\(N-1\) is coprime to \(N\). The polygon criterion proves
(1). For \(N=2\) strict real positivity of the quadratic
excludes the negative real root of the linear moment.
The absolute angular theorem handles \(N\ge J\).

## Density among unmodified profiles

Choose one variable row and reverse the columns so it
increases. For degree \(N<J\), write its smaller exponent
as \(A\). If a prime \(p\) sufficiently large for the
fixed gaps and fixed rows divides \(NA+1\) exactly once
and divides no other rising factor at order \(N\), then
\(E_+=1\), \(E_-=U=V=0\), and (1) holds.
Such a prime need not exceed the moving exponents.

To find such native primes for most profiles, let the
progression parameters vary independently between \(1\)
and \(H\). For large \(H\), put
\(R=\sqrt{\log H}\), \(P=e^{R/4}\), and \(\sigma=1-1/R\).
In the designated row, \(NA+1\) is an injective positive
affine integer function of its parameter and is at most
\(C H\) for a fixed \(C\ge1\). The number of its values with
every prime factor at most \(P\) is at most
\[
            (C H)^\sigma
          \prod_{\substack{q\le P\\ q\ {\rm prime}}}
                    (1-q^{-\sigma})^{-1}.
\]
Indeed each such integer at most \(CH\) contributes at least
one after weighting it by \((CH/x)^\sigma\), and the product
sums these weights over all integers with those prime factors.
For large \(H\), \(\sigma\ge3/4\), so the sum of the prime
power terms of order at least two in the logarithm of the
product is uniformly bounded by a constant times
\(\sum_{m\ge2}m^{-3/2}\). The first-power terms have sum
at most \(\int_1^P t^{-\sigma}\,dt=R(e^{1/4}-1)\).
After division by \(H\) the smooth-value bound is therefore
\(O(e^{-(2-e^{1/4})R})=O(1/P)\).
The number of parameters for which \(NA+1\) is divisible
by \(p^2\) for any \(p>P\) is \(O(H/P+\sqrt H)\).
For each such \(p\), which is at most \(O(\sqrt H)\),
the congruence has at most \(H/p^2+1\) solutions.
For sufficiently large \(H\), all primes above \(P\)
avoid the fixed progression coefficient and the fixed
row arguments. Summing gives the bound.

Outside this exceptional set choose a prime \(p>P\)
dividing \(NA+1\) exactly once. At large \(H\) this prime
also exceeds \(Nd\) in the marked row, so it divides no
other factor there. Conditional on the marked parameter
and this prime, each factor in another variable row
\(N\ell+h\), \(1\le h\le Nd_i\), is divisible by \(p\)
for at most \(H/p+1\) values of that row's independent
parameter. The union of all such events has fraction
\(O(1/P+1/H)\). This proves a failure proportion
\(O(1/P+1/\sqrt H)\) at order \(N\). The union over the
fixed set \(3\le N<J\) proves the stated density bound.
No factorization or hypothesis at the lower order is used.
