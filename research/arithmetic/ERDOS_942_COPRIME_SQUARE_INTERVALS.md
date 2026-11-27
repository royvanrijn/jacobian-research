# Erdős #942 and coprime occupants of square intervals

Write \(h(q)=\#\{t:q^2\le t<(q+1)^2,\ t\text{ powerful}\}\).
The left endpoint \(q^2\) is powerful, so the count in the
open interval is \(h(q)-1\).
[Erdős #942](https://www.erdosproblems.com/942), from
[Erdős, printed p. 35](https://www.renyi.hu/~p_erdos/1976-39.pdf#page=11),
asks for matching logarithmic power upper and infinitely often lower
bounds. The upper bound remains open.

The [short power window theorem CPW1](ABC_COPRIME_SHORT_POWER_WINDOWS.md#a-positive-density-primitive-construction)
gives positive density of intervals with any fixed number of pairwise
coprime nonsquare powerful occupants. Its constants depend on that number.
Selecting a larger centre for each successive fixed multiplicity
also gives qualitative unboundedness.
[Narumi–Tachiya, Theorems 1–2 and Lemmas 1–2](https://arxiv.org/html/2611.20656v1)
use simultaneous approximation to obtain ordinary occupancy of order
at least \(\log q/\log\log q\) infinitely often. Their approximation
conditions do not impose pairwise coprimality.

## Coprime occupancy with quadratic logarithmic height

**Theorem (Coprime occupancy with exponent one half).** For every sufficiently
large integer \(m\), there is an open interval \((q^2,(q+1)^2)\) containing
at least \(3m/8\) pairwise coprime nonsquare powerful integers, with
\(q\ge m\) and \(\log q=O(m^2\log m)\). In particular the coprime occupancy
is \(\gg(\log q/\log\log q)^{1/2}\) along an unbounded sequence of \(q\).

Use the same Laurent cancellation as below, with prime shifts of size
\(O(m\log m)\). Interpolation over the small prime completions allows
these shifts to stay small. Once the argument \(B^2\) is fixed, we use
the rank of the evaluation forms in the remaining coefficients to
bound the number of pairs that share a larger prime.

### Cancellation and small-prime interpolation

Choose \(P_0=\lceil C m\log m\rceil>4m+5\) with a sufficiently large
absolute constant \(C\). The elementary prime count below supplies \(m\)
distinct odd primes \(d_1,\ldots,d_m\le P_0\) from one residue class
modulo \(8\). Put \(n=2m\), \(c_j=\binom{-3/2}{j}\) and
\[
 R_0(T)=\prod_i(T-d_i)=\sum_{j=0}^m a_jT^j,\qquad
 S(T)=\sum_{\ell=0}^m s_\ell T^\ell,\qquad
 R_0(T)S(T)=\sum_{k=0}^n r_k T^k .
\]
The coefficients of \(S\) will be rational with power-of-two denominators.
Set
\[
 \lambda_k(Z,d)=\frac{1}{c_{k+1}}\sum_{j=0}^k c_j d^j Z^{k-j},
 \qquad
 P(Z)=\sum_{k=0}^n\frac{r_k}{c_{k+1}} Z^k,\qquad
 L_i(Z)=\sum_{k=0}^n r_k\lambda_k(Z,d_i).
\]
The coefficient of \(Z^{-1}\) in \(P(Z)(1+d_i/Z)^{-3/2}\) is
\(d_iR_0(d_i)S(d_i)=0\), and \(L_i(Z)\) is the polynomial part.
For \(E=2T\frac{d}{dT}+3\) and \(\Phi=E^{-1}(R_0S)\), we also have
\(L_i(0)=\Phi(d_i)\). Indeed \(E\) multiplies \(T^k\) by \(2k+3\),
and \(c_k/c_{k+1}=-1+1/(2k+3)\).

For an odd prime \(p\le P_0\), construct
\(\Phi_p\in\mathbb Z_p[T]\) of degree less than \(2m\) such that
\(E\Phi_p(d_i)=0\) for all \(i\) and \(\Phi_p(d_i)\) is a unit
whenever \(d_i\) is. For the shifts in a nonzero class \(a\bmod p\)
use the series \((1+(T-a)/a)^{-3/2}\) in the quotient by the product
of \((T-d_i)^2\) over that class, taking \(a\) to be an integer
representative. The binomial coefficients are \(p\)-integral. If
there are \(r\) shifts in the class, \((T-a)^{2r}\) is divisible by
\(p\) in the quotient, so this series converges there. Evaluation
modulo each \((T-d_i)^2\) gives its value and first derivative at \(d_i\).
The value is \(1\bmod p\) and the series satisfies \(Ef=0\).
For the zero class use the zero polynomial. The defining products for
distinct classes are comaximal over \(\mathbb Z_p\), since they have
unit resultants. The polynomial Chinese remainder theorem gives
\(\Phi_p\).

At \(p=2\), let \(a\) be an odd representative of the common class
modulo \(8\), and write \(U=(T-a)/8\). In the finite free
\(\mathbb Z_2\)-algebra obtained by quotienting by
\(\prod_i(U-(d_i-a)/8)^2\), use \((1+8U/a)^{-3/2}\).
Its \(j\)-th coefficient has valuation at least \(j\), so the series
converges. Its polynomial representative of degree less than \(2m\)
gives a polynomial \(\Phi_2\) in \(T\) with coefficient valuations at
least \(-6m\). Evaluation modulo the squared linear factors gives
unit values at every \(d_i\) and the same differential equation
\(E\Phi_2(d_i)=0\), by the change of variable. The quotient maps
used in both constructions are continuous on the finite free
coefficient modules, so taking the limits preserves the value and
derivative identities.

For each \(p\le P_0\), \(S_p=E\Phi_p/R_0\) is a polynomial of degree
less than \(m\). Division by the monic \(R_0\) preserves the coefficient
valuation bounds. Put
\[
 e_p=1+\lfloor\log_p(4m+3)\rfloor,\qquad D_0=2^{6m},\qquad
 K=D_0\prod_{p\le P_0}p^{e_p},\qquad X=m^{12}.
\]
Choose \(A_\ell\in[0,K)\) by the integer Chinese remainder theorem
so that
\[
 D_0s_\ell=A_\ell+Kh_\ell,\qquad X\le h_\ell<2X,\qquad h_\ell\in\mathbb Z,
\]
gives \(v_p(S-S_p)\ge e_p\) coefficientwise at every \(p\le P_0\).
For this, approximate the \(p\)-integral coefficient of \(D_0S_p\)
to precision \(p^{e_p+v_p(D_0)}\). Applying \(E^{-1}\) loses at
most \(e_p-1\) powers of \(p\). Thus every \(L_i(0)\) is integral at
these primes, with residue \(\Phi_p(d_i)\).

Require \(B\) to be divisible by \(K_B=\prod_{p\le P_0}p^{b_p}\),
where \(b_2=10m+10\) and \(b_p=1+2\lfloor\log_p(4m+3)\rfloor\)
at odd \(p\). Then \(L_i(B^2)\equiv L_i(0)\bmod p\) at all these
primes. To see this, each term of positive \(Z\)-degree in
\(L_i(Z)\) has valuation at least \(1\) upon substituting \(Z=B^2\).
At odd \(p\) the \(r_k,c_j\) are integral and
\(v_p(c_{k+1})\le 2\lfloor\log_p(4m+3)\rfloor\), by the binomial
formula and the factorial valuation formula. At \(2\) we have
\(v_2(r_k)\ge-6m\), \(v_2(c_j)\ge-2j\) and \(v_2(c_{k+1})\le0\),
also by the binomial formula. Beyond \(P_0\), the denominators in
\(L_i\) are units. Consequently \(L_i(B^2)\) is an integer. Also,
every prime \(p\le P_0\) can divide at most one of
\((B^2+d_i)L_i(B^2)\), namely the one with \(d_i=p\) if such a
shift is present.

### Choosing the square argument

The primorial bound below gives \(\log K,\log K_B=O(m\log m)\).
Indeed the extra powers beyond the primorial contribute at most
\(O(m\log m)\) in total, since only \(p\le4m+3\) contribute them.
The coefficient heights of \(P\) and all \(L_i\) are bounded by a
uniform \(H=\exp(O(m\log m))\) for all the allowed \(h_\ell\).
This follows from \(\sum |a_j|\le\prod_i(1+d_i)\), the bounds
on \(s_\ell\) and \(d_i\), and \(1\le |c_j|\le 2j+1\).
The common leading coefficient has modulus
\(|s_m/c_{n+1}|\ge 1/(2n+3)\), since \(s_m\ge1\).
Hence the Laurent remainder estimate below supplies a threshold
of logarithm \(O(m\log m)\) above which \(B^2+d_i\) is strictly
between consecutive squares, \(P(B^2)\) and all \(L_i(B^2)\)
are nonzero, and
\[
 \big|(B^2+d_i)^{3/2}|L_i(B^2)|-|B^3P(B^2)|\big|<\tfrac12 .
\]
The same threshold can ensure \(|B^3P(B^2)|\ge m+2\), using
leading-term dominance uniformly in \(S\).

Choose \(B\) beyond this threshold, with \(K_B\mid B\) and with
no prime \(P_0<p\le Y=m^3\) dividing any \(B^2+d_i\).
Write \(B=K_B b\). For each of these primes at most \(2m\) classes
of \(b\) are forbidden, so their densities \(v_p\) satisfy
\(v_p<1/2\). The bound \(\pi(x)\ll x/\log x\) and partial summation
give \(\sum_{P_0<p\le Y}1/p=O(1)\). Thus
\(V=\prod_p(1-v_p)\ge \exp(-2\sum_p v_p)=\exp(-O(m))\).
As in the progression estimate below, truncate inclusion-exclusion
at a sufficiently large odd order \(O(m)\). With \(S_0=\sum_p v_p\),
the absolute main-term tail is bounded by
\(\sum_{j>k}S_0^j/j!\) at order \(k\), which is at most \(V/2\)
for \(k\) a sufficiently large odd integer of size \(O(S_0+1)\).
The integer counting error is at most the sum of the forbidden
class counts for the subsets used, bounded by
\(\sum_{j\le k}(2mY)^j=\exp(O(m\log m))\).
An interval of \(b\) of length \(\exp(O(m\log m))\), starting beyond
the prescribed threshold, then contains an admissible \(b\). This
gives \(\log B=O(m\log m)\). The numbers \(B^2+d_i\) are pairwise
coprime. At \(p\le P_0\) only \(d_i=p\) can contribute, and at
\(p>P_0\) the differences of the shifts are nonzero modulo \(p\).

### Rank and large-prime collisions

Fix \(Z=B^2\). The \(m\) forms \(L_i(Z)\) in the coefficients of \(S\)
have rank \(m\) over every \(\mathbb F_p\) with \(p>P_0\).
For a relation with weights \(w_i\), put
\(\mu_k=\sum_i w_i\lambda_k(Z,d_i)\). Vanishing of the relation
on coefficients of \(S\) says \(\sum_{j=0}^m a_j\mu_{j+\ell}=0\)
for \(0\le\ell\le m\). Since \(R_0\) is monic and the shifts
are distinct modulo \(p\), the Vandermonde matrix on the first
\(m\) terms and this recurrence yield
\(\mu_k=\sum_i C_i d_i^k\) for \(0\le k\le2m\).
The identity
\[
 (2k+3)\lambda_k(Z,d_i)+2(k+1)Z\lambda_{k-1}(Z,d_i)
   =-2(k+1)d_i^k,\qquad 1\le k\le2m,
\]
therefore implies
\[
 \sum_i\bigl((2k+3)C_i+2(k+1)(ZC_i/d_i+w_i)\bigr)d_i^k=0
 \qquad(1\le k\le2m).
\]
The sequences \(d_i^k\) and \(k d_i^k\) are independent on these
indices. Indeed, relabeling \(k=j+1\) identifies relations
between them with relations between evaluations and first
derivatives on polynomials of degree less than \(2m\), since the
\(d_i\) are nonzero. These evaluation and derivative functionals
are independent by the Chinese remainder theorem modulo
\(\prod_i(T-d_i)^2\). For each \(i\), subtracting the coefficient
of \(k\) in the last displayed relation from its constant
coefficient gives \(C_i=0\), and then \(w_i=0\).

After the substitution \(D_0s_\ell=A_\ell+Kh_\ell\), clear
denominators by a common factor supported at primes at most \(P_0\),
for instance \(D_0\,4^{2n+1}(2n+3)!\). Denote the resulting
integer affine forms by \(G_i\). Their linear parts still have
rank \(m\) over every \(p>P_0\), and hence over \(\mathbb Q\).
Their coefficient logarithmic heights are \(O(m^2\log m)\)
because \(Z=B^2\) and their degree before evaluation at \(Z\)
was at most \(2m\). Primes \(p>P_0\) divide \(L_i(Z)\) exactly
when they divide \(G_i\).

Choose the \(h_\ell\) independently and uniformly in their
integer ranges. For each pair \(i,j\) and each prime
\(P_0<p\le X\), use a nonsingular two-column minor of the two
linear parts and fix all other coordinates. The probability
that \(p\) divides both \(G_i,G_j\) is at most
\((1/p+1/X)^2\le 4/p^2\). The expected number of pairs with a
collision of this kind is \(O(m^2/P_0)=o(m)\).

For \(p>X\), any collision has some coordinate in which at least
one of the two linear coefficients is nonzero modulo \(p\).
For each such candidate coordinate \(\ell\), fix the other
coordinates, and write the forms as \(u h_\ell+v\) and
\(u' h_\ell+v'\). We can skip coordinates with \(u=u'=0\).
The eliminant \(u'v-uv'\) is a nonzero affine form in the other
coordinates, since the two linear parts are independent over
\(\mathbb Q\). Its probability of taking the value zero is at
most \(1/X\). Otherwise its logarithmic size is \(O(m^2\log m)\),
so it has \(O(m^2)\) prime divisors larger than \(X\). For each
such prime at which at least one of \(u,u'\) is a unit, at most one
choice of \(h_\ell\) gives the collision. Summing over
coordinates gives \(O(m^3/X)\) per pair.

Finally consider a prime \(p>P_0\) dividing \(Z+d_i\) and \(G_j\)
for \(j\ne i\). It must be greater than \(Y\). For \(p\le X\),
rank one of the single linear part gives probability at most
\(2/p\). For each \(i\) the sum of \(1/p\) over such divisors of
\(Z+d_i\) is bounded by \(\log(Z+d_i)/(Y\log Y)\).
For \(p>X\) the probability is at most \(1/X\) and the number
of these divisors is at most \(\log(Z+d_i)/\log X\).
Summing over ordered pairs gives \(O(m^3/Y+m^3/X)\).
Together with pairwise coprimality of the bases and the
small-prime conditions, the expected number of pairs sharing a
prime in \((Z+d_i)L_i(Z)\) is \(o(m)\).

Take a choice with at most \(m/4\) such pairs, and remove at most
one index per pair to retain at least \(3m/4\) pairwise coprime
powerful integers \((Z+d_i)^3L_i(Z)^2\). Each is nonsquare
because \(Z+d_i\) is nonsquare. Their square roots lie within
\(1/2\) of the same real number, so they meet at most two open
unit intervals with integer endpoints. At least \(3m/8\) lie
in one such interval \((q,q+1)\). The threshold ensured \(q\ge m\),
and the coefficient and argument bounds give
\(\log q=O(n\log B+\log H)=O(m^2\log m)\). This proves the
claimed lower bound since \(x/\log x\) is increasing for \(x>e\).
Taking a sufficiently large constant multiple of a prescribed
multiplicity parameter gives the hypotheses and conclusion of the
cubic logarithmic height theorem below.

## Earlier cubic logarithmic height bound

**Theorem (Coprime occupancy with exponent one third).** For each integer
\(m\ge2\), there is an interval \((q^2,(q+1)^2)\) containing at least
\(\lceil m/2\rceil\) pairwise coprime nonsquare powerful integers, where
\[
 m\le q,\qquad \log q\ll m^3\log(m+1)
\]
with an absolute constant. In particular the coprime occupancy is
\(\gg \bigl(\log q/\log\log q\bigr)^{1/3}\) along an unbounded
sequence of \(q\).

In the construction below, small primes are handled by scaling the shifts.
For the intermediate primes, extra coefficients in the cancellation
polynomial allow all the values to be units at a prescribed argument.
This reduces the size of the shifts, which determines how far out the
square-root approximations become accurate.

Let \(n=2m+1\), \(d=2m(n+1)\), \(C=8(m+1)\), and put
\[
 A=4^{2n+1}(2n+3)!,\qquad M=A\prod_{p\le C}p.
\]
Use \(m\) shifts \(d_i=M t_i\), where the \(t_i\) will be primes
larger than \(2d\) and of size \(m^{O(1)}\). Let
\[
 \begin{split}
 S(T)&=T^{m+1}+\sum_{\ell=0}^m s_\ell T^\ell,\qquad
 R(T)=T\prod_i(T-d_i)S(T)=\sum_{k=0}^n r_{k+1}T^{k+1},\\
 P(Z)&=\sum_{k=0}^n p_k Z^k,\qquad
 p_k=c_{n+1}r_{k+1}/c_{k+1},\qquad c_j=\binom{-3/2}{j},\\
 L_i(Z)&=\sum_{k=0}^n p_k\sum_{j=0}^k c_j d_i^j Z^{k-j},
 \qquad F_i(Z)=(Z+d_i)L_i(Z).
 \end{split}                                                     \tag{942.1}
\]
The three polynomials in \(Z\) are monic. Since \(R(d_i)=0\),
the same Laurent coefficient calculation as in the construction below
shows that the \(Z^{-1}\) term of \(P(Z)(1+d_i/Z)^{-3/2}\) vanishes.
Take all \(s_\ell\equiv0\bmod M\). Then all coefficients of \(R\)
below the leading one are multiples of \(M\). The denominators of the
rational factors \(c_{n+1}c_j/c_{k+1}\) in (942.1) divide \(A\),
by the binomial formula for \(c_j\). Each nonleading term in a
coefficient of \(P\) or \(L_i\) has a factor \(r_{k+1}\) with \(k<n\)
or a positive power of \(d_i\). The polynomials are therefore integral
and all their nonleading coefficients are divisible by \(M/A\).
This holds also for \(F_i\). Thus \(B\equiv1\bmod p\) keeps all
\(F_i(B^2)\) nonzero modulo every \(p\le C\).

### Choosing the shifts and the intermediate prime conditions

First work with \(S(T)=T^{m+1}\). For each pair \(i\ne j\) the
resultant of \(F_i,F_j\) is a nonzero polynomial in the shifts.
Indeed, specialize only \(d_i=a\ne0\) and set the other shifts to zero.
With \(\lambda=c_{n+1}/c_n\), this gives
\(P(Z)=Z^{n-1}(Z-\lambda a)\) and \(F_j(Z)=Z^n(Z-\lambda a)\).
The polynomial part formula gives
\[
 L_i(0)=a^n(c_n-\lambda c_{n-1}),\qquad
 L_i(\lambda a)=c_n a^n .
\]
Both values are nonzero since \(c_j/c_{j-1}\) varies with \(j\).
Also \(Z+a\) is nonzero at both arguments since \(\lambda\ne-1\).
The product of the pairwise resultants has total degree \(O(m^4)\)
in the shifts, using a Sylvester matrix of size \(2(n+1)\) and the
bound \(n+1\) on each coefficient degree. Substituting \(d_i=M t_i\)
preserves nonvanishing as a polynomial in the \(t_i\).
The grid lemma below therefore allows a choice of the \(t_i\)
among any \(O(m^4)\) distinct primes above \(2d\), with the implicit
constant chosen larger than the degree bound. There are sufficiently
many such primes of size \(m^{O(1)}\). Indeed
\(\binom{2u}{u}\ge4^u/(2u+1)\), and the exponent of each prime \(p\)
in this binomial coefficient is at most \(\lfloor\log_p(2u)\rfloor\):
each power \(p^j\le2u\) contributes zero or one by the factorial
valuation formula. Comparing logarithms gives
\(\pi(2u)\gg u/\log(2u)\). Choose \(u\) a sufficiently large
fixed power of \(m\) times an absolute constant. Fix the resulting
vector of shifts, for which the pairwise resultants are nonzero at
\(s_0=\cdots=s_m=0\).

Fix now a prime \(C<p\le2d\). The \(d_i\) and the \(c_j\),
\(0\le j\le n+1\), are units modulo \(p\). We will take
\(B\equiv0\bmod p\), so it suffices to make each \(L_i(0)\) nonzero.
Write
\[
 T\prod_j(T-d_j)=\sum_{k=0}^m a_k T^{k+1},\qquad a_m=1 .
\]
For \(0\le\ell\le m\), the coefficient of \(s_\ell\) in \(L_i(0)\)
is
\[
 c_{n+1}d_i^\ell
 \sum_{k=0}^m a_k d_i^k\frac{c_{k+\ell}}{c_{k+\ell+1}}
 =
 c_{n+1}d_i^\ell
 \sum_{k=0}^m\frac{a_k d_i^k}{2(k+\ell)+3}.                 \tag{942.2}
\]
The constant \(-1\) in the quotient
\(c_j/c_{j+1}=-1+1/(2j+3)\) disappears because
\(\sum_k a_k d_i^k=0\).
If all \(m+1\) coefficients (942.2) vanished, the rational function
\[
 X\longmapsto\sum_{k=0}^m\frac{a_k d_i^k}{X+2k+3}
\]
would vanish at \(X=0,2,\ldots,2m\) in \(\mathbb F_p\).
None of these points is a pole because \(p>C\). Its numerator over
the common denominator has degree at most \(m\), so it would be the
zero polynomial. The poles \(-2k-3\) are distinct, and the numerator
evaluated at \(-2m-3\) is \(d_i^m\prod_{k<m}(2k-2m)\ne0\).
This contradiction shows that \(L_i(0)=0\) excludes an affine
hyperplane of coefficient vectors modulo \(p\). The \(m\) hyperplanes
together have at most \(m p^m<p^{m+1}\) points. Choose a vector
outside their union for each intermediate prime.

The prescribed residues of \(\mathbf s\), including \(0\bmod M\),
give a Cartesian grid of progressions with common step
\(K=M\prod_{C<p\le2d}p\), with \(\log K=O(m^2)\) by the
primorial bound below. For our fixed shifts the resultant product is
a nonzero polynomial in \(\mathbf s\), as seen at the zero vector.
Its total degree is \(O(m^3)\) since the coefficients of each
\(F_i\) are affine in \(\mathbf s\). Take the first sufficiently many
nonnegative integers in each progression and apply the grid lemma
again. We obtain \(\log(2+\max_\ell|s_\ell|)=O(m^2)\), preserving all
intermediate prime conditions and all nonzero resultants.

### Choosing the square interval

The chosen shifts satisfy \(\log(2+\max_i d_i)=O(m\log(m+1))\).
It follows from (942.1) that the logarithms of all coefficient
heights of \(P,L_i,F_i\) are \(O(m^2\log(m+1))\). To see this,
the product \(T\prod_i(T-d_i)\) has coefficient sum bounded by
\(\prod_i(1+d_i)\), its product with \(S\) costs at most the
coefficient sum of \(S\), and the polynomial part formula uses
at most \(O(m^2)\) terms with one power \(d_i^j\), \(j\le n\).
The bounds \(1\le|c_j|\le2j+1\) control the other factors.

Let \(R_1\) be the product of the absolute values of the pairwise
resultants of the chosen integral \(F_i\). The Sylvester determinants
give \(\log(2R_1)=O(m^5\log(m+1))\). Enforce the chosen residues
\(B\equiv1\bmod p\) for \(p\le C\) and \(B\equiv0\bmod p\) for
\(C<p\le2d\). Their common step is \(\prod_{p\le2d}p\), of logarithm
\(O(m^2)\). A prime \(p>2d\) outside \(R_1\) cannot divide two
different \(F_i(B^2)\). For the remaining primes \(p>2d\) dividing
\(R_1\), exclude the roots of the monic polynomial
\(\prod_i F_i(B^2)\) of degree \(d\). There are at most \(d\)
forbidden residues per prime. The progression estimate below
uses only the step, the bound \(d\) on the number of excluded classes,
and \(J=1+\lceil\log_2 R_1\rceil\) as an upper bound for the number
of primes. In particular it works for the present prescribed residue
just as for \(B\equiv1\) on the small primes. In its notation \(Y=100(d+1)J=m^{O(1)}\), and only primes
\(2d<p\le Y\) contribute to the small-prime sum \(S\).
The primorial bound gives \(\pi(x)\ll x/\log x\) by counting primes
above \(\sqrt x\) with weight at least \((\log x)/2\).
Partial summation therefore gives
\(\sum_{2d<p\le Y}1/p=O(1)\) because \(\log Y/\log(2d)=O(1)\).
Thus \(S=O(d)\), and the error factor in the progression estimate
has \(\log E=O(d\log(m+1))\). With this bound the estimate supplies
an admissible \(B\) within \(\exp(O(m^2\log(m+1)))\) beyond any
positive threshold.

Use the Laurent remainder estimate below with degree \(n\), shift bound
\(D=\max_i d_i\) and coefficient bound \(H\) for \(P,L_i\).
A threshold \(B\) greater than a sufficiently large constant times
\((n+2)^2 H D^{n+2}\) ensures the common square-root error is less than
one, and its logarithm is \(O(m^2\log(m+1))\).
It also ensures that every \(B^2+d_i\) is strictly between \(B^2\)
and \((B+1)^2\) and that all polynomial values in (942.1) are
positive. We can therefore choose \(B\) with
\(\log B=O(m^2\log(m+1))\).

The integers \((B^2+d_i)^3 L_i(B^2)^2\) are pairwise coprime
since their prime supports are those of the \(F_i(B^2)\). They
are powerful and nonsquare, and their square roots differ from
the common integer \(Q=B^3P(B^2)\) by less than one.
At least half lie in one of the two open intervals of length one
with endpoint \(Q\). In either case its lower endpoint \(q\) is
at least \(m\) by the choice of threshold, and
\(\log q=O(n\log B+\log H)=O(m^3\log(m+1))\).
This proves the theorem. It gives a lower bound for coprime occupancy;
the logarithmic uniform upper bound in Erdős #942 remains open.

## Growing coprime multiplicity

**Theorem (Growing coprime occupancy).**
For each integer \(m\ge2\) there is an interval \((q^2,(q+1)^2)\)
containing at least \(\lceil m/2\rceil\) pairwise coprime nonsquare
powerful integers, with
\[
 m\le q,\qquad \log q\ll m^4
\]
with an absolute constant. Consequently along an unbounded sequence
of \(q\) the coprime occupancy is
\(\gg(\log q)^{1/4}\).
In particular the lower bound \((\log q)^{1/5}\) holds for all
sufficiently large \(q\) in that sequence.

We seek integers of the form \((B^2+d_i)^3L_i(B^2)^2\), with one
common integer square root approximation \(B^3P(B^2)\).
Truncating a binomial expansion alone leaves an error of order
\(B\) in the square root. We choose \(P\) so the first discarded
coefficient vanishes for all \(i\). This reduces the error to
\(O(1/B)\) with the construction data fixed. Resultants then let
us choose \(B\) with all the powerful values pairwise coprime.

### Cancellation at all selected shifts

Let \(c_j=\binom{-3/2}{j}\), so \(c_0=1\) and
\[
 c_j=(-1)^j(2j+1)\binom{2j}{j}/4^j,\qquad
 c_{j+1}/c_j=-(2j+3)/(2j+2).
\]
For a vector of shifts \(\mathbf d=(d_1,\ldots,d_m)\), define
\[
 R(T)=T\prod_i(T-d_i)=\sum_{j=1}^{m+1}r_jT^j,\qquad
 P(Z)=\sum_{k=0}^m p_k Z^k,\qquad p_k=c_{m+1}r_{k+1}/c_{k+1}.
\]
In particular \(P\) is monic and \(p_k\) is homogeneous of degree
\(m-k\) in the shifts. Define the monic polynomial
\[
 L_i(Z)=\sum_{k=0}^m p_k\sum_{j=0}^k c_j d_i^j Z^{k-j},\qquad
 F_i(Z)=(Z+d_i)L_i(Z).
\]
These polynomials depend on \(\mathbf d\). In the Laurent expansion
of \(P(Z)(1+d_i/Z)^{-3/2}\) at infinity, \(L_i\) is the polynomial
part. The coefficient of \(Z^{-1}\) is
\(\sum_k p_k c_{k+1}d_i^{k+1}=c_{m+1}R(d_i)=0\).

For example, with \(m=2\) and \((d_1,d_2)=(24,48)\), this yields
\[
 P=Z^2+84Z+1680,\quad L_1=Z^2+48Z-264,\quad L_2=Z^2+12Z-48 .
\]
At \(B=90725\) the two resulting powerful integers are coprime
and both exceed \((B^3P(B^2)-1)^2\) while lying below
\((B^3P(B^2))^2\). These claims can be checked by evaluating the
three integer polynomials and taking a gcd. The proof for arbitrary
\(m\) follows.

### Making the resultants nonzero

For any pair \(i\ne j\), the resultant in \(Z\) of \(F_i,F_j\)
is a nonzero polynomial in the shifts. To prove this, specialize
\(d_i=a\ne0\) and all other shifts to zero. Put
\(\lambda=c_{m+1}/c_m\). Then
\[
 P(Z)=Z^{m-1}(Z-\lambda a),\qquad
 F_j(Z)=Z^m(Z-\lambda a),
\]
so it suffices to check \(F_i\) at \(0\) and at \(\lambda a\).
The polynomial part formula gives
\[
 L_i(0)=a^m(c_m-\lambda c_{m-1}),\qquad
 L_i(\lambda a)=c_m a^m .
\]
For the second equality the two truncated sums from the two
terms of \(P\) cancel except for the term of index \(m\).
The first quantity is nonzero because \(c_m/c_{m-1}\ne\lambda\).
The factors \(Z+a\) are nonzero at both points because
\(a\ne0\) and \(\lambda\ne-1\). This proves the assertion about
the resultant.

The coefficient of \(Z^k\) in \(F_i\) has degree at most
\(m+1\) in the shifts. Its Sylvester determinant with \(F_j\)
thus has degree at most \(2(m+1)^2\). The product over pairs
has degree at most \(m(m-1)(m+1)^2\). A nonzero polynomial of
total degree \(D\) cannot vanish on a Cartesian grid with more
than \(D\) entries in each coordinate. This follows by expanding
in one coordinate, using the univariate root bound, and inducting
on the number of coordinates. Thus we can choose positive integers
\(t_i\le 1+(m+1)^4\) for which all the resultants in the vector
\(\mathbf t\) are nonzero. In particular the \(t_i\) are distinct.

### Choosing the integer parameter

Put
\[
 A=4^{2m+1}(2m+3)!,\qquad d=2m(m+1),\qquad
 M=A\prod_{p\le 2d}p,\qquad d_i=Mt_i.
\]
For the vector \(\mathbf t\), all coefficients of \(P,L_i,F_i\)
have denominator dividing \(A\). Indeed
\((2k+3)\binom{2k+2}{k+1}\) divides \((2m+3)!\) for
\(0\le k\le m\), and the powers of \(4\) from
\(c_{m+1}/c_{k+1}\) and an additional \(c_j\) divide
\(4^{2m+1}\). This also clears the denominators in the sums.
Scaling \(\mathbf t\) by \(M\) multiplies the coefficient of
\(Z^k\) by \(M^{m-k}\) in \(P,L_i\), and \(M^{m+1-k}\) in
\(F_i\). Thus the polynomials for \(\mathbf d\) have integer
coefficients, and all the coefficients below the leading one are
divisible by every prime at most \(2d\).

Use \(F_i(Z;\mathbf t)\) to denote the unscaled polynomials, and put
\[
 R_0=\prod_{i<j}
 \left|\operatorname{Res}_Z(A F_i(Z;\mathbf t), A F_j(Z;\mathbf t))\right|.
\]
This is a positive integer. For a prime \(p\) not dividing \(M R_0\),
no two of the integers \(F_i(B^2;\mathbf d)\) can be zero modulo
\(p\). Otherwise \(B^2/M\) modulo \(p\) would be a common root of
the corresponding unscaled polynomials, by homogeneity and
\(p\nmid M\). Since \(A\) is also a unit, the prime would divide
\(R_0\).

Every prime dividing \(M\) is at most \(2d\). Impose
\(B\equiv1\bmod \operatorname{rad}(M)\). Then all
\(F_i(B^2;\mathbf d)\) are one modulo the primes dividing \(M\).
For every other prime dividing \(R_0\), exclude the zeros
of the monic polynomial \(\prod_i F_i(B^2;\mathbf d)\) in \(B\).
It has degree \(d\), so at most \(d\) classes are forbidden
modulo any such prime \(p>2d\).

For the unscaled vector, each coefficient of \(R(T)\) has modulus
at most \(\prod_i(1+t_i)\), and \(1\le |c_j|\le2j+1\).
The formulas for \(P,L_i,F_i\) give
\(\log\max(2,H(A F_i(Z;\mathbf t)))=O(m\log(m+1))\),
where \(H\) is maximum absolute coefficient.
The Sylvester determinants have size \(2(m+1)\), so
\(\log R_0=O(m^4\log(m+1))\). This estimate uses the unscaled
polynomials, which keeps \(M\) out of the product of resultants.
Also \(\log M=O(m^2)\). To see the primorial bound, the product of
primes in \((l,2l]\) divides \(\binom{2l}{l}\le4^l\).
Sum the logarithmic bound over intervals between powers of two.

Write \(B=1+\operatorname{rad}(M)b\) and let \(b\) run through
\(T\) consecutive integers. Each \(p\mid R_0\), \(p\nmid M\)
excludes \(b_p\le d\) residue classes of \(b\). Put
\[
 J=1+\lceil\log_2 R_0\rceil,\qquad Y=100(d+1)J .
\]
There are at most \(J\) primes to consider. For those at most \(Y\)
put \(g_p=b_p/p\), \(S=\sum g_p\) and \(V=\prod(1-g_p)\).
We have \(g_p<1/2\), \(S=O(d\log Y)\) by the harmonic sum,
and \(V\ge e^{-2S}\). By the Chinese remainder theorem any
\(j\) prescribed prime exclusions intersect in a fraction
\(\prod g_p\) of classes modulo their product. The same
fraction holds in a progression whose step is coprime to all
those primes. In \(L\) successive terms of such a progression
the count differs from \(L\prod g_p\) by at most \(\prod b_p\).

Use inclusion–exclusion truncated at consecutive even and
odd orders of size at least \(10(S+1)\) and at most \(10(S+1)+3\).
The main term tails are at most
\(\sum_{j\ge10(S+1)}S^j/j!\le V/8\), by \(j!\ge(j/e)^j\)
and \(V\ge e^{-2S}\). The total counting error is at most
\(E=(1+dY)^{O(S+1)}\), with an absolute constant in the exponent,
by summing \(\prod b_p\) up to those orders. Thus the count
passing the primes at most \(Y\) is between \((7/8)VL-E\)
and \((9/8)VL+E\), uniformly over those progressions.
Take \(E\ge1\). Its logarithm is \(O(m^2(\log(m+1))^2)\).

For \(Y<p\le\sqrt T\), apply the upper small prime count to
each forbidden step \(p\) progression with \(L\le T/p+1\).
Summing the main terms costs at most \(2VTdJ/Y+2VdJ\),
and the errors cost at most \(EdJ\). The primes exceeding
\(\sqrt T\) cost at most \(dJ(\sqrt T+1)\), even without
a small prime restriction. Take
\[
 T=\left\lceil\big(100(d+1)(J+1)(E+1)/V\big)^2\right\rceil .
\]
The total exclusions are smaller than the lower small prime
count at \(L=T\), since \(dJ/Y\le1/100\) and
\(V\sqrt T\ge100(d+1)(J+1)(E+1)\). An admissible \(B\)
therefore occurs within \(O(\operatorname{rad}(M)T)\) beyond
any prescribed positive threshold. The logarithm of this
length is \(O(m^2(\log(m+1))^2)\).

### Passing to one open interval

Set \(D=\max_i d_i\) and let \(H\ge1\) be the maximum coefficient
modulus over \(P,L_i\) for \(\mathbf d\). The same coefficient
bounds give \(\log(2D)=O(m^2)\) and
\(\log(2H)=O(m^3)\). For \(Z=B^2\ge2D\) the binomial
series converges absolutely. The contribution to the remainder
from the \(k\)-th term after omitting its polynomial part
and the \(Z^{-1}\) term is at most
\[
 |p_k| Z^k\sum_{j\ge k+2}(2j+1)(d_i/Z)^j
 \ll (k+2)|p_k|D^{k+2}Z^{-2}.
\]
The \(Z^{-1}\) terms sum to zero by the cancellation identity.
Consequently for positive \(B\), as soon as \(L_i(B^2)>0\),
\[
 \left| B^3P(B^2)-(B^2+d_i)^{3/2}L_i(B^2)\right|
 \ll (m+2)^2 H D^{m+2}/B .
\]
Choose a common threshold for \(B\) larger than an absolute
constant times \((m+2)^2 H D^{m+2}\). This also ensures
\(d_i<2B+1\), \(B>m+2\), and positivity of \(P(B^2)\) and
all \(L_i(B^2)\), by monicity. The logarithm of the threshold
is \(O(m^3)\). Apply the preceding progression count beyond this
threshold. Its length bound gives \(\log B=O(m^3)\).

Each \(T_i=(B^2+d_i)^3L_i(B^2)^2\) is powerful. The base \(B^2+d_i\)
lies strictly between \(B^2\) and \((B+1)^2\), so \(T_i\) is
nonsquare. Its prime support equals that of \(F_i(B^2;\mathbf d)\).
Thus the \(T_i>1\) are pairwise coprime and distinct. Their
square roots lie at distance strictly less than one from the
integer \(Q=B^3P(B^2)\), and none equals \(Q\). At least half
lie on the same side of \(Q\). The two possible intervals
have lower square roots \(q=Q-1\) and \(q=Q\).
Since \(P\) is monic and the threshold makes its leading term
dominate, \(q\ge m\). Also \(\log Q=O(m\log B+\log H)=O(m^4)\).
This proves the theorem.

The argument is an elementary polynomial construction with an elementary
coprimality sieve. It gives no logarithmic upper bound on \(h(q)\).
