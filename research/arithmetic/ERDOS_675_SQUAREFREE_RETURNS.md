# Squarefree return times in Erdős #675

For a fixed integer \(r\ge2\), let \(S_r\) be the positive integers
having no prime divisor to exponent at least \(r\). Let \(t_r(N)\)
be the least positive shift preserving the indicator of \(S_r\) at
every position \(1\le a\le N\). Erdős asked about the growth of
\(t_2(N)\), including whether it exceeds \(\exp(N^c)\) for some
\(c>0\).
[Erdős, §5, p. 69](https://www.renyi.hu/~p_erdos/1979-22.pdf)

Ho's lower bound answers the latter question. The October 2026
version gives \(t_2(N)>\exp(N^c)\) for every fixed \(c<25/72\)
and the upper bound \(\log t_2(N)\ll N/\log N\).
The count below gives \(\log t_2(N)\ll\sqrt N\) and applies to
every fixed \(r\). The asymptotic order of the unrestricted
least shift remains open.
[Ho's note, Theorem 1.1 and Proposition 1.3](https://boonsuan.github.io/erdos675_squarefree.pdf)

The [density of all returns](#natural-density-of-all-squarefree-returns)
has a sharp logarithmic asymptotic.
The [prime-divisibility profile](#prime-divisibility-patterns-of-power-free-returns)
describes which prime powers typically divide a return and gives
counts at finite height when the negative sites have bounded witness primes.
The [common divisor of all returns](#the-common-divisor-of-the-returns)
is already realized at height \(\exp(O(\sqrt N))\).

## Returns divisible by the small prime powers

**Theorem (Uniform count of power-free prefix returns).**
Put \(Y=N^{1/r}\), \(M=\prod_{p\le Y}p^r\), and
\(k=|S_r\cap[1,N]|\). For a constant \(C_r>0\), uniformly for
integers \(U\ge\exp(C_rY)\), the number of \(1\le u\le U\) for
which \(Mu\) preserves the prefix is
\[
 (1+o_r(1))\,U V,\qquad V=\prod_{p>Y}(1-k/p^r),
\]
as \(N\to\infty\) with \(r\) fixed. In particular
\(\log t_r(N)\ll_r N^{1/r}\). Among shifts divisible by \(M\),
the logarithm of the least positive return is of order \(N^{1/r}\).

Use a progression \(t=b+Mu\) with \(0\le b<M\), with \(p^r\mid M\)
for every \(p\le Y\), and with \(p^r\mid M\) at every prime dividing
\(M\). Suppose the fixed residues modulo those \(p^r\) preserve all
the positive sites and supply a witness at every negative site.
For each fixed \(B>0\), the argument below gives the progression
count \((1+o_{r,B}(1))U V_M\), where
\[
 V_M=\prod_{p\nmid M}(1-k/p^r),
\]
uniformly if \(\log M\le BY\) and \(U\ge \exp(C_{r,B}Y)\).
For \(b=0\) and the primorial power \(M=\prod_{p\le Y}p^r\), the
fixed residues preserve every positive site and retain all the
original negative witnesses. Then \(V_M=V\).

For each \(p\nmid M\) there are exactly \(k\) bad residue classes of
\(u\bmod p^r\), namely \(-(a+b)M^{-1}\) at the positive offsets
\(a\). The classes are distinct because \(p>Y\) and \(p^r>N\).
The proof controls the truncation error with growing \(N\).
Mirsky's fixed-pattern product is a precedent for this local
calculation.
[Mirsky, *Arithmetical Pattern Problems Relating to Divisibility by \(r\)th Powers*](https://doi.org/10.1112/plms/s2-50.7.497)

Set \(Z=N^2\), let \(\mathcal P\) contain the primes \(Y<p\le Z\) with \(p\nmid M\),
and write \(V_0=\prod_{p\in\mathcal P}(1-q_p)\), \(q_p=k/p^r\).
The power-free count from Möbius inversion is
\(k=N/\zeta(r)+O_r(Y)\). In fact its main term comes from
summing \(\mu(d)N/d^r\) over \(d\le Y\), and the floor errors
and the omitted tail are both \(O_r(Y)\). Thus
\(q_p\le\eta_r<1\) for a fixed \(\eta_r\) and all sufficiently
large \(N\). By the prime number theorem and partial summation,
\(\sum_{p>Y} p^{-r}\ll_r Y^{1-r}/\log Y\). Therefore
\[
 S=\sum_{p\in\mathcal P} q_p\ll_r Y/\log Y,\qquad
 -\log V_0\ll_r Y/\log Y .
\]
Also \(\sum_{p>Z} q_p=O_r(N Z^{1-r})=o(1)\), so
\(V_M=(1+o(1))V_0\).

Truncate inclusion-exclusion at two consecutive orders \(L,L+1\),
where \(L=\lceil A_rY/\log Y\rceil\) for a fixed sufficiently
large \(A_r\). The elementary symmetric sum of degree \(j\)
in the \(q_p\)'s is at most \(S^j/j!\). The total omitted
mass \(\sum_{j>L}S^j/j!\) is \(o(V_0)\). Indeed, choosing
\(L\ge2eS\) bounds this sum by \(2^{-L}\), and increasing
\(A_r\) makes \(2^{-L}/V_0\) tend to zero. For a subset of \(j\)
primes, the Chinese remainder theorem gives precisely \(k^j\)
bad classes modulo the product of the prime powers. Each
class has count equal to its expected interval length up to
an absolute error at most one. Thus the total error in either
truncation is at most
\[
 E=\sum_{j\le L+1}(\#\mathcal P\, k)^j\le\exp(O_r(Y)).
\]
The same absolute error holds if we count only within one
specified residue class modulo an extra modulus coprime to all
finite sieve moduli. The main term is divided by that extra modulus.
In particular this works modulo \(p^r\) for \(p>Z\). Consequently the number of parameters passing the
\(\mathcal P\) tests is \(UV_0(1+o(1))+O(E)\). In the
specified residue class the count is
\((U/p^r)V_0(1+o(1))+O(E)\), with the relative error in the
main term uniform in the extra prime and class.

Take \(W=U^{1/4}\), so \(W>Z\) for sufficiently large \(N\).
For \(Z<p\le W\) with \(p\nmid M\), sum the bound in the \(k\) bad classes
over \(p\). The resulting loss from the \(\mathcal P\) count
is at most
\[
 (1+o(1))UV_0\, k\sum_{p>Z}p^{-r}+O(kWE)=o(UV_0).
\]
Here \(C_{r,B}\) can be chosen large enough in \(U\ge\exp(C_{r,B}Y)\)
to absorb the bound on \(E\) and the reciprocal of \(V_0\).
For \(p>W\) with \(p\nmid M\), use the unsieved residue counts. No prime
\(p>(M(U+1)+N)^{1/r}\) can occur. The loss is at most
\[
 O_r\bigl(kU W^{1-r}+k(M(U+1)+N)^{1/r}\bigr)=o(UV_0).
\]
Since \(\log M\le BY\), increasing \(C_{r,B}\) gives \(M\le U^{1/4}\), so both terms
are bounded by \(O_r(kU^{5/8}+kU^{3/4})\). This estimate and
the previous error estimates hold uniformly for the stated
range of \(U\). Thus the full count is \((1+o(1))UV_0\),
and replacing \(V_0\) by \(V_M\) proves the progression count.
For the primorial power the prime number theorem gives \(\log M=(r+o(1))Y\). At
\(U=\lceil\exp(C_rY)\rceil\) the count is positive, so
\(t_r(N)\le MU\). Any positive multiple of \(M\) is at least
\(M\), giving the order statement for the restricted shifts.

### Size of the local product

If \(\alpha_r=1/\zeta(r)\), then
\[
 -\log V\ \sim\ \frac{Y}{\log Y}
           \int_1^\infty -\log(1-\alpha_r x^{-r})\,dx .
\]
For fixed \(B>1\), the prime number theorem gives this
approximation with the sum restricted to \(Y<p\le BY\) and
the integral restricted to \([1,B]\). This follows first for
step functions from prime counts in each subinterval, then
by uniform approximation of the continuous integrand. Replacing
\(k/N\) by \(\alpha_r\) changes the sum by \(o_r(Y/\log Y)\),
using the bounded derivative in the first argument on
\([0,\eta_r]\) and \(\sum_{p>Y}(p/Y)^{-r}\ll_r Y/\log Y\).
The omitted tail is \(O_r((Y/\log Y)B^{1-r})\), by partial
summation. Letting \(B\) tend to infinity gives the formula.

## Natural density of all squarefree returns

**Theorem (Sharp density of squarefree prefix returns).**
Let \(\delta_N\) be the natural density of positive shifts preserving the
squarefree indicator on \([1,N]\). Then
\[
 -\log\delta_N=\bigl(1+1/\sqrt2+o(1)\bigr)\sqrt N.
\]
The density takes the limit in the shift bound for each fixed \(N\).
The unrestricted least-shift estimate requires control in growing finite
intervals and remains open.

Put \(H=\sqrt N\). The density equals the probability of the return event
when the residues \(t_q\bmod q^2\), indexed by the primes, are independent
and uniform. For a finite prime cutoff \(Z\) this is the Chinese remainder
theorem. For each offset \(a\), the number of \(1\le t\le X\) with some
\(q>Z\) satisfying \(q^2\mid t+a\) is at most
\[
 X\sum_{q>Z}q^{-2}+O(\sqrt{X+N}).
\]
The corresponding probability in the product space is at most
\(\sum_{q>Z}q^{-2}\). Summing over the \(N\) offsets and taking \(X\to\infty\)
before \(Z\to\infty\) proves the assertion. The fixed-pattern formulation
also follows from Mirsky's correlations by inclusion-exclusion.
[Cellarosi and Sinai, §3](https://ems.press/content/serial-article-files/31903)

### Lower bound on the density

Fix \(1/\sqrt2<b<1\). Impose \(t_q=0\) at every \(q\le bH\),
an event of probability \(\exp(-(2b+o(1))H)\) by the prime number theorem.
These primes preserve every positive site and cover each nonsquarefree
site except \(p^2\) for \(bH<p\le H\). There are
\(k=(1-b+o(1))H/\log H\) such sites.

For \(bH<q\le H\), a residue represented by exactly one offset \(a\) can
be used without losing a positive site if \(4\mid a\). Every \(a\) in
\([N-q^2+1,q^2]\) is alone in its class. The multiples of 4 in this
interval provide a fraction bounded below in terms of \(b\) of the
\(q^2\) possible choices of \(t_q=-a\). This step costs
\(\exp(-O_b(H/\log H))\) in probability.

Use the \(m=(2+o(1))H/\log H\) primes in \((2H,4H]\) to cover
the \(k\) remaining sites. For each injection from the sites to these
primes, prescribe \(t_q=-p^2\bmod q^2\) at the assigned prime and
require every unused prime in the interval to avoid all \(N\) sites.
Each assigned residue has probability at least \((4H)^{-2}\). Each
unused prime has allowed probability at least \(3/4\).
Different injections give disjoint events because each such prime can
cover at most one offset and the unused primes cover none. Their union
has probability at least
\[
 m(m-1)\cdots(m-k+1)(4H)^{-2k}(3/4)^m
 \ \ge\ \exp(-(1-b+o(1))H).
\]
At every other prime \(q>H\) it suffices to avoid the positive sites.
Writing \(K=|S_2\cap[1,N]|\), this has probability \(1-K/q^2\).
The product of these probabilities is \(\exp(-O(H/\log H))\)
since \(K/N\to1/\zeta(2)<1\) and
\(\sum_{q>H}q^{-2}\ll1/(H\log H)\). Independence gives
\(\delta_N\ge\exp(-(1+b+o(1))H)\). Let \(b\) decrease to \(1/\sqrt2\).

### Upper bound on the density

First, every return has \(t_q=0\) for \(q\le N^{3/10}\), once \(N\)
is sufficiently large. Here is an elementary squarefree witness
argument. For \(Q=q\) or \(q^2\), every reduced class modulo \(Q\)
contains a squarefree positive integer \(\le C Q^{3/2}\) for an
absolute \(C\). For large odd \(q\), set \(x=C Q^{3/2}\) and
\(L=x^{1/3}\). Prime squares \(\ell^2\) with \(\ell\le L\) exclude
at most
\[
 \frac{x}{Q}\sum_\ell\ell^{-2}+O(x^{1/3})
\]
integers in the class, by the Chinese remainder theorem. The term with
\(\ell=q\) excludes none. For \(\ell>L\) write \(n=v\ell^2\).
There are at most \(x/L^2\) possible \(v\), each coprime to \(q\).
For every such \(v\), the congruence for \(\ell^2\bmod Q\) has at most
two roots, and \(\sqrt x<Q\) for sufficiently large \(q\) with \(C\)
fixed. This gives \(O(x^{1/3})\) further exclusions. Since the sum of
prime reciprocal squares is less than 1, taking \(C\) sufficiently
large leaves a squarefree integer. For the finitely many remaining
moduli the reduced class has positive squarefree density, by the same
finite-cutoff and tail argument as above with the class fixed.
Increasing \(C\) handles them. Thus every nonzero class modulo \(q^2\)
has a squarefree representative of size \(O(q^3)\). For the classes
divisible by \(q\), take \(q\) times the squarefree unit representative
modulo \(q\), giving \(O(q^{5/2})\). If \(t_q\ne0\), the witness in
the class \(-t_q\) contradicts preservation of a positive site.

Fix \(0<\eta<1/\sqrt2\) and \(0<\epsilon<\eta\).
For each prime \(\eta H<p\le H\), use the site \(p^2\), and also
use \(2p^2\) when \(p\le H/\sqrt2\). Call this group of one or two
sites the row for \(p\). A prime \(q\) covers a site \(a\) when
\(t_q=-a\bmod q^2\). For large \(N\), a prime \(q\le N^{3/10}\)
covers none of the chosen sites by the forcing just proved. A larger
\(q\) covers at most one site of the form \(p^2\) and at most one of
the form \(2p^2\). Indeed, for distinct odd primes \(p,p'\) the
factors \(p-p'\) and \(p+p'\) have greatest common divisor 2 in
absolute value. If \(q^2\) divided \(p^2-p'^2\), it would divide one
factor, whose absolute value is at most \(2H<q^2\).
The same argument applies to \(2p^2-2p'^2\) because \(q\) is odd.

Call a prime \(q\ge\epsilon H\) exceptional if
\(q^2\mid p^2-2p'^2\) for two chosen sites from distinct rows.
The primes \(p,p',q\) are pairwise distinct, since divisibility by
\(p^2\) or \(p'^2\) would force equality of the two row primes.
They satisfy
\[
 p^2-2p'^2=j q^2,\qquad 0<|j|\le \epsilon^{-2},\qquad q\le H.
\]
The conic estimate below, applied to the finitely many \(j\),
gives \(o_{\eta,\epsilon}(H/\log H)\) exceptional primes.

For a given return, discard every row in which a site is covered by
a prime \(N^{3/10}<q<\epsilon H\) or by an exceptional prime.
At most \(2\pi(\epsilon H)+o_{\eta,\epsilon}(H/\log H)\) rows are
discarded, by the at-most-two-site bound. Among the remaining rows,
specify the subset with \(t_p=0\). For each of these the condition at
\(p\) costs \(p^{-2}\le(\eta H)^{-2}\). Each site in the other
remaining rows requires a distinct prime \(q\ge\epsilon H\).
Two such sites in the same row could share a prime only if
\(q=p\) and \(t_p=0\). Sites in different rows could share a prime
only if it were exceptional. The covering primes of these sites
are also distinct from the primes prescribed to have \(t_p=0\),
since a zero residue there covers only its own chosen row.
By independence and the union bound over injective assignments,
if there are \(v\) sites requiring separate primes, their probability
factor is at most
\[
 \Bigl(\sum_{q\ge\epsilon H}q^{-2}\Bigr)^v
 \ \le\ \bigl(C_\epsilon/(H\log H)\bigr)^v .
\]
Thus a remaining row of length \(d\in\{1,2\}\) costs at least
\(d\log H-O_{\eta,\epsilon}(1)\) in negative log probability.
There are at most \(3^{\pi(H)}=\exp(o(H))\) choices of the discarded
rows and the subset with zero residue. The initial sum of row lengths
is \((1+1/\sqrt2-2\eta+o(1))H/\log H\). Discarding rows removes at
most \(4\pi(\epsilon H)+o_{\eta,\epsilon}(H/\log H)\) sites. Therefore
\[
 -\log\delta_N
 \ \ge\ \bigl(1+1/\sqrt2-2\eta-4\epsilon-o_{\eta,\epsilon}(1)\bigr)H .
\]
Let \(\epsilon\) and \(\eta\) decrease to zero.

### Prime points on the collision conics

For fixed nonzero integers \(A,B,C\) and fixed \(c>0\), there are
\(O_{A,B,C,c}(H/(\log H)^2)\) solutions to
\(A p^2+B p'^2+Cq^2=0\) in pairwise distinct primes in \([cH,H]\).
The conic is nonsingular. If it has a rational point, projection from
that point identifies it with \(\mathbb P^1\). To see the degree of
the parametrization, move the rational point to \([0:0:1]\) over
\(\mathbb Q\) and write the equation as
\(L(X,Y)Z+Q(X,Y)=0\), with \(L\) linear and \(Q\) quadratic.
Then \([u:v]\mapsto [uL(u,v):vL(u,v):-Q(u,v)]\).
The unique zero of \(L\) maps to the chosen point. It cannot also
be a zero of \(Q\) since the conic contains no line. Undoing the
coordinate change and clearing denominators gives three integer binary
quadratic forms \(F_1,F_2,F_3\) without a common projective zero.
No point on the diagonal conic has two zero coordinates, so
any two of the forms have no common projective zero.

For coprime integers \(u,v\), the greatest common divisor of the
three form values is bounded by a constant. Indeed the univariate
Bézout identities for \(F_1,F_2\) in the charts \(v=1\) and \(u=1\),
homogenized and with denominators cleared, express fixed nonzero
multiples of powers of \(v\) and of \(u\) in the ideal of the forms.
Coprimality of \(u,v\) bounds this common divisor. A prime triple
is already primitive, so its coordinates are the form values divided
by this common divisor, up to a common sign. The forms have no common
real zero on the unit circle, whence their largest absolute value
is \(\gg u^2+v^2\). All the triples in question can therefore be
counted using a box \(|u|,|v|\ll\sqrt H\).

Use \(F_1F_2\) for the sieve. Each form has nonzero discriminant
because each coordinate line meets the diagonal conic at two
distinct geometric points. Away from a fixed finite set of primes
the two forms still have disjoint root sets in the projective line.
The congruence \(F_1F_2\equiv0\bmod\ell\) in the affine plane has
\(1+\rho(\ell)(\ell-1)\) solutions, where
\(\rho(\ell)=2+\chi_1(\ell)+\chi_2(\ell)\) and \(\chi_i\) is the
quadratic character of the form discriminant, possibly principal.
If \(D=H^{1/100}\), every prime \(\ell<D\) outside the finite set
can be used for sifting. Enlarge the exceptional set to omit primes
up to the bound on the common divisor. Neither form value of a
prime triple can vanish at any of the sifting primes, since
\(cH>D\) for large \(H\).

Let \(g(\ell)=(1+\rho(\ell)(\ell-1))/\ell^2\) and extend \(g\)
multiplicatively to squarefree products of sifting primes. For a
modulus \(d\), the Chinese remainder theorem gives exactly
\(g(d)d^2\) residue pairs. If the box contains \(V\ll H\) integer
pairs, the count in these residues is
\(Vg(d)+R_d\), where \(|R_d|\ll \sqrt H\,d+d^2\) by counting
the two coordinates in each class. The Selberg upper bound gives
\[
 \frac{V}{G(D)}+
 O\left(\sum_{a,b<D}|R_{[a,b]}|\right),
 \qquad
 G(D)=\sum_{d<D}\prod_{\ell\mid d}\frac{g(\ell)}{1-g(\ell)}.
\]
The sums use squarefree products of the sifting primes. In
Richert's notation the local function \(\omega(\ell)/\ell\)
is \(g(\ell)\). Omit primes up to 11 so that \(0<g(\ell)\le1/2\).
The error is \(O(\sqrt H D^4+D^6)\).
[Richert, Theorem 9.1, printed pp. 106–107](sources/richert-sieve-methods-1976.pdf)

For \(z=D^\gamma\) with sufficiently small fixed \(\gamma>0\),
the product of \((1-g(\ell))^{-1}\) over the usable \(\ell\le z\)
is \(\gg(\log z)^2\). This follows by taking logarithms and
applying Mertens' theorem in each reduced class of the quadratic
characters' fixed moduli. Nonprincipal characters contribute
\(O(1)\) to the weighted prime reciprocal sum, while principal
ones contribute \(\log\log z+O(1)\).
[Williams, §3, Theorem 1](https://people.math.carleton.ca/~williams/papers/pdf/057.pdf)

In the expanded product, normalize by the total mass and regard
each prime as independently selected with probability \(g(\ell)\).
For the resulting product \(d\), \(\mathbb E\log d=
\sum_{\ell\le z}g(\ell)\log\ell\ll\log z\), by the prime number
theorem and \(g(\ell)\ll1/\ell\). Choose \(\gamma\) so that this
expectation is less than \((\log D)/2\). At least half the mass
then comes from \(d<D\), proving \(G(D)\gg(\log H)^2\).
The Selberg bound gives the claimed estimate for the conic.

## Prime divisibility patterns of power-free returns

Put \(Y=N^{1/r}\) and define
\[
 K_r(u)=\#\{j\in S_r: j\le u^{-r}\}\qquad(0<u\le1).
\]
Thus \(K_r(p/Y)\) counts the sites \(jp^r\) for which the
coefficient \(j\) is itself \(r\)-free. Write
\(\vartheta(x)=\sum_{p\le x}\log p\).

**Theorem Erdos675PowerFreePrimeProfiles.** Fix \(r\ge2\) and
\(0<a<1\). Let \(Z\) be any subset of the primes in \((aY,Y]\).
Let \(\delta_{r,N}(a;Z)\) be the natural density of returns \(t\)
with \(p^r\mid t\) at all \(p\le aY\) and exactly the primes in
\(Z\) above \(aY\) through \(Y\). Uniformly in \(Z\),
\[
 -\log\delta_{r,N}(a;Z)
 =r\vartheta(aY)+\log Y\!
   \left(r|Z|+(r-1)\sum_{\substack{aY<p\le Y\\p\notin Z}}K_r(p/Y)\right)
    +o_{r,a}(Y).                                                   \tag{675.P}
\]
In particular, with \(Z=\varnothing\), the leading coefficient of
\(Y\) is \(ra+(r-1)\int_a^1 K_r(u)\,du\).
Write \(A_Z\) for the displayed right side without its error.
There are constants \(B_{r,a},C_{r,a}>0\) such that uniformly
for \(X\ge\exp(C_{r,a}Y)\), the number of these returns
\(t\le X\) for which each negative site has a witness prime
\(q\le B_{r,a}Y\) is
\[
                     X\exp(-A_Z+o_{r,a}(Y)).
\]
In particular every prescribed profile has a return by
\(\exp(C'_{r,a}Y)\) with a common constant.

**Theorem Erdos675PowerFreeDensityConcentration.** For every fixed
\(r\ge2\), the natural density \(\delta_{r,N}\) of all returns satisfies
\[
             -\log\delta_{r,N}
          = (r-1+2^{-1/r}+o_r(1))Y.                                \tag{675.Q}
\]
At a prime \(p\le Y\), compare the costs \(r\) for \(p^r\mid t\)
and \((r-1)K_r(p/Y)\) for \(p^r\nmid t\). For any fixed
\(a_0>0\) and \(\lambda>0\), the proportion, in the natural-density
limit among returns, that make the more expensive choice at
\(\lambda Y/\log Y\) or more primes in \([a_0Y,Y]\) is at most
\(\exp(-cY)\), with \(c>0\) depending on \(r,a_0,\lambda\).
Ties do not count. For \(r\ge3\) the less expensive choice is
divisibility below \(2^{-1/r}Y\) and nondivisibility above it,
with the boundary assigned by \(K_r\). For \(r=2\) the prices
tie on \(Y/\sqrt3<p\le Y/\sqrt2\).

### Density and collisions between rows

For each fixed \(N\), the density equals a probability for
independent uniform residues \(t_q\bmod q^r\) at the primes.
The CRT establishes this with a finite prime cutoff. For each
site, the tail above \(Q\) has upper density at most
\(\sum_{q>Q}q^{-r}\): up to height \(X\), the floor errors total
at most \(O((X+N)^{1/r})\). Let \(X\) tend to infinity first
and then \(Q\) to infinity. The same argument gives a density
with any prescribed choices at finitely many primes.

Fix \(a>0\). For \(aY<p\le Y\), use the row of sites
\(jp^r\) with \(j\in S_r\) and \(j\le (Y/p)^r\).
There are at most \(J=\lceil a^{-r}\rceil\) sites per row.
All coefficients \(j\) are bounded, and different rows
have disjoint sites for large \(N\). A prime \(q\) covers a
site if \(q^r\) divides its translate.

If \(q>c_r N^{1/(r+1)}\), then \(q\) covers at most one site
with a given coefficient \(j\). For distinct \(p,p'\), the
divisibility \(q^r\mid p^r-p'^r\) would require \(q\nmid pp'\).
Since \(q\nmid r\), the cyclotomic factors of \(X^r-1\) have
no shared root modulo \(q\). Thus one homogeneous cyclotomic
factor of \(p^r-p'^r\) must be divisible by \(q^r\), while every
such factor has absolute value at most \((2Y)^{r-1}\).
The latter is smaller than \(q^r\) for all sufficiently large
\(N\), since \(r^2/(r+1)>r-1\). The factor \(j\) is coprime
to \(q\) in this range. In particular \(q\) covers at most \(J\)
sites in total.

For fixed \(\epsilon>0\), consider a prime \(q\ge\epsilon Y\)
that can cover sites in two different rows. It then satisfies
\[
 jp^r-j'p'^r=kq^r,\qquad 0<|k|\le\epsilon^{-r}.
\]
There are finitely many coefficient triples for fixed
\(r,a,\epsilon\), and \(q\le Y\). The three primes are distinct:
if \(q=p\) or \(q=p'\), the congruence forces the other row
prime to equal it, since the coefficients are bounded.
For \(r=2\), the prime-point estimate on collision conics
above gives \(o(Y/\log Y)\) such primes \(q\). For \(r\ge3\),
the displayed ternary form defines a smooth projective
plane curve and is geometrically irreducible, since any two
positive-degree components would meet and give a singularity.
The prime triples are primitive of height at most \(Y\).
Heath-Brown's bound for an irreducible ternary form of degree
\(r\) gives \(O_{r,\nu}(Y^{2/r+\nu})=o(Y/\log Y)\) of them,
taking \(0<\nu<1-2/r\) and summing over the fixed coefficient
triples. Call these primes exceptional.
[Heath-Brown, Theorem 3, equation (1.9)](https://ora.ox.ac.uk/objects/uuid%3A61776df1-3232-4fdc-87de-f3c01f6bce3f/files/m59be1594382f3e2554c57fc1628d3c42)

### Prescribing the divisibility choices

Fix the zero residues in (675.P), at cost
\(\exp(-r\vartheta(aY))\prod_{p\in Z}p^{-r}\).
Every negative site without a witness at one of those primes
is in a row for a prime \(p\notin Z\). Indeed an offset with
two large prime powers would exceed \(N\), and if its cofactor
\(j\le a^{-r}\) has an \(r\)-th power divisor, that divisor
gives a witness at a prime below \(aY\).

Let \(m=\sum_{p\notin Z,\ aY<p\le Y}K_r(p/Y)\). A prime
\(q\le aY\) with zero residue covers none of these \(m\) sites.
The same holds for a prime in \(Z\). A prime \(q>aY\) can
cover two sites within a row only if \(q=p\) and \(t_p=0\),
which is excluded for these rows. Drop all sites that could
be covered by an exceptional prime \(q>aY\). At most
\(o(Y/\log Y)\) sites are dropped by the bound of \(J\)
per prime. The other sites require distinct covering
primes. Independence and summation over injections bound
the probability, after fixing the required zero residues, by
\[
 \left(\sum_{q>aY}q^{-r}\right)^{m-o(Y/\log Y)},
\]
where the exponent is the actual number of retained sites.
The prime number theorem gives
\(\log\sum_{q>aY}q^{-r}=-(r-1)\log Y+O_{r,a}(\log\log Y)\).
This proves the lower bound on the negative logarithm in (675.P).

For the upper bound on the negative logarithm, at each
\(aY<q\le Y\) outside \(Z\), there is an allowed fraction
bounded below by a positive constant depending only on
\(r,a\) for choosing a nonzero residue. To construct such
residues, for each \(0\le j\le J\) assign a distinct fixed
prime \(\ell_j\), all different from \(q\) once \(N\) is large,
and prescribe \(b=-jq^r\bmod\ell_j^r\).
The CRT leaves a fixed positive fraction of choices
\(1\le b<q^r\) for which every \(b+jq^r\) is non-\(r\)-free.
Take \(t_q=-b\). This class preserves every positive
offset and is nonzero. Across these primes the cost is
\(\exp(-O_{r,a}(Y/\log Y))\).

Assign the \(m\) negative sites to distinct primes from
a fixed interval \((LY,2LY]\), with \(L>2\) large enough
in terms of \(r,a\) to provide at least \(2m\) primes for
large \(N\). Prescribe the covering residue at each assigned
prime and avoid all \(N\) sites at each unused prime in
that interval. A prime there covers at most one site,
so different injections give disjoint events.
Writing \(M\asymp_L Y/\log Y\) for the number of primes,
the probability is at least
\[
 (M/2)^m(2LY)^{-rm}(1-L^{-r})^M
     =\exp(-(r-1)m\log Y+o_{r,a}(Y)).
\]
This estimate also holds for \(m=0\). Every other prime
\(q>Y\) need only avoid the positive sites. Their product
costs \(\exp(-O_r(Y/\log Y))\) by the product estimate used
in the uniform count above. The required zero residues
and these choices prove (675.P), using
\(\sum_{p\in Z}\log p=|Z|\log Y+O_a(\pi(Y))\).
The Riemann sum giving its integral consequence follows
from the prime number theorem at each of the finitely many
jumps of \(K_r\) on \([a,1]\).

For the finite-height lower count, retain all the safe choices
and all injections just used. Take one CRT modulus \(M\) with a
factor \(q^r\) at every prime up to \(Y\) and every prime in
\((LY,2LY]\). Its logarithm is \(O_{r,a}(Y)\). The allowed
classes have fraction at least \(\exp(-A_Z+o(Y))\) modulo
\(M\). Each class certifies the negative sites and preserves
the positives at the fixed primes. The progression version of
the uniform count applies with \(U=\lfloor X/M\rfloor-1\)
once \(C_{r,a}\) is large. Its product over the remaining
primes costs \(\exp(-O_r(Y/\log Y))\). Summing over the
allowed classes proves the lower count.

For the upper count take \(B_{r,a}=2L\). Use the same
exceptional-site removal and injective witness count as in
the density upper bound, but restrict the witness primes to
\((aY,2LY]\). The zero conditions and each specified
injection give one CRT class. Its count through \(X\) is at
most its probability times \(X\) plus one. There are at
most \(\pi(2LY)^m=\exp(O_{r,a}(Y))\) injections, so the
summed floor error is negligible beside
\(X\exp(-A_Z+o(Y))\) by increasing \(C_{r,a}\).
This gives the claimed uniform logarithmic count.

### Small forced prime powers

**Lemma Erdos675PowerFreeForcedPrimePowers.** Fix \(r\ge2\).
Every return satisfies \(q^r\mid t\) for
\(q\le c_r N^{1/(r+1)}\), for a fixed \(c_r>0\).
In particular \(\log t_r(N)\gg_r N^{1/(r+1)}\).

An elementary reduced-class bound supplies the witnesses. For
\(Q=q^e\), \(1\le e\le r\), each reduced class modulo \(Q\)
contains an \(r\)-free integer at most \(C_r Q^{1+1/r}\).
For large \(q\), set \(X=C_r Q^{1+1/r}\) and
\(D=\lfloor X^{1/(r+1)}\rfloor\). Prime powers \(\ell^r\)
with \(\ell\le D\) exclude at most
\[
        (X/Q)\sum_\ell \ell^{-r}+O(D).
\]
Here \(\ell=q\) excludes none. For \(\ell>D\), write
\(n=v\ell^r\) with \(v\le X/D^r\) coprime to \(q\). For
each \(v\) the congruence for \(\ell^r\bmod Q\) has at most
\(r\) roots by Hensel lifting when \(q\nmid r\), and
\(X^{1/r}<Q\) for large \(q\). The large powers therefore
exclude \(O_r(X^{1/(r+1)})\) terms. Since
\(\sum_\ell\ell^{-r}<1\), choosing \(C_r\) sufficiently
large leaves an \(r\)-free representative. For the finitely
many small \(q\), the same CRT cutoff and tail argument
gives positive \(r\)-free density in each reduced class,
so increasing \(C_r\) works. A nonzero class modulo \(q^r\)
of valuation \(v<r\) is handled by multiplying the unit
witness modulo \(q^{r-v}\) by \(q^v\). Its size is
\(O_r(q^{r+1})\) and it is still \(r\)-free. Applied to
\(-t\bmod q^r\), any such witness up to \(N\) would destroy
a positive site. This proves the forced divisibility claim.
The product of the forced prime powers has logarithm
\(r\vartheta(c_rN^{1/(r+1)})\gg_r N^{1/(r+1)}\).

### Unrestricted returns and the costly choices

Fix \(0<\epsilon<a<1\) and use the rows above \(aY\).
The forced small primes do not cover any of their sites.
Discard every row touched by a witness
\(c_rN^{1/(r+1)}<q<\epsilon Y\) or an exceptional prime
\(q\ge\epsilon Y\). For each return this discards at most
\(J\pi(\epsilon Y)+o_{r,a,\epsilon}(Y/\log Y)\) rows.
Specify the discarded rows and the zero choices on the
others in \(\exp(o(Y))\) ways. Once those choices are made,
each retained zero row costs \(p^{-r}\). Every site in a
retained nonzero row needs a separate covering prime
\(q\ge\epsilon Y\), none of which is a prescribed zero
prime. As before summation over the witness assignments
costs at most \((\sum_{q\ge\epsilon Y}q^{-r})^v\) for \(v\)
such sites. Thus the cost per intact row, in units of
\(\log Y\) and up to \(o(1)\), is at least the smaller of
\(r\) and \((r-1)K_r(p/Y)\), plus one at every strict
departure from this minimum. The sum of the smaller costs
over all rows, times \(\log Y\), is
\[
 \left(\int_a^1 \min\{r,(r-1)K_r(u)\}\,du+o_r(1)\right)Y .
\]
The discarded rows lose \(O_{r,a}(\epsilon Y)+o(Y)\) from
this bound. Let \(\epsilon\) tend to zero first and then
\(a\) to zero. The integrand is \(r\) for \(u\le2^{-1/r}\)
and \(r-1\) above it, proving the upper density bound.
For the lower density bound use (675.P) with
\(a=2^{-1/r}\), \(Z=\varnothing\), where \(K_r=1\) above
the cutoff.

For concentration use \(a<a_0\) so small that \(ra<\lambda/8\),
then choose \(\epsilon>0\) small enough that the discard
cost, including one per lost strict departure, is below
\(\lambda Y/8+o(Y)\). If there are at least
\(\lambda Y/\log Y\) strict departures in the test range,
the preceding estimate for that subevent has negative
logarithm at least
\((r-1+2^{-1/r}+\lambda/2+o(1))Y\). Dividing by the density
in (675.Q) proves the asserted exponential suppression.
Both density statements take the limit in the shift bound
with \(N\) fixed. They do not bound the first unrestricted
squarefree return on the scale \(\log t_2(N)\asymp\sqrt N\).


## The second-order return density

**Theorem Erdos675SecondOrderReturnDensity.** Fix \(r\ge2\)
and put \(Y=N^{1/r}\) and \(\beta=2^{-1/r}\). The natural
density \(\delta_{r,N}\) of the power-free prefix returns satisfies
\[
 -\log\delta_{r,N}
  =(r-1+\beta)Y
   +(1-\beta)\frac{Y\log\log Y}{\log Y}
   +O_r\!\left(\frac{Y(\log\log Y)^{1/4}}{\log Y}\right).
                                                               \tag{675.R}
\]
For \(r\ge3\) the error is \(O_r(Y/\log Y)\). For \(r=2\)
and any fixed \(\lambda>0\), the proportion in the natural-density
limit among returns with \(p^2\nmid t\) at at least
\(\lambda Y/\log Y\) primes \(p\le\beta Y\) is at most
\[
   \exp\!\left(-(2\lambda+o(1))Y\log\log Y/\log Y\right).
\]
In particular this resolves the leading-order cost tie for
the primes in \((Y/\sqrt3,Y/\sqrt2]\).

The same expansion holds for the negative logarithm of the
fraction of shifts \(t\le X\) that are returns with a witness
prime \(q\le B_rY\) at each negative site, uniformly for
\(X\ge\exp(C_rY)\) and suitable fixed \(B_r,C_r>0\).
For unrestricted returns the expansion holds uniformly for
\(X\ge\exp(C'_rY^2/\log Y)\), with another fixed constant.

### The cost of distinct witnesses

For \(v=O_r(Y/\log Y)\) specified sites, requiring distinct
witness primes has probability at most
\[
 v!\,e_v((q^{-r})_q)
  \le \exp\bigl(-(r-1)v\log Y-v\log\log Y+O_r(Y/\log Y)\bigr),
                                                               \tag{675.S}
\]
where the sum defining the elementary symmetric function runs
over all primes. Any primes already prescribed for other tests
must be excluded from the witnesses. Allowing them in the purely
numerical sum only increases the bound. Indeed,
\[
 e_v((q^{-r})_q)
 \le Y^{-rv}\prod_q(1+(Y/q)^r),\qquad
 \log\prod_q(1+(Y/q)^r)\ll_r Y/\log Y .
\]
For the product, primes above \(Y\) contribute
\(O_r(Y^r\sum_{q>Y}q^{-r})\). Those at most \(\sqrt Y\)
contribute \(O_r(\sqrt Y\log Y)\). Dyadic intervals between
\(\sqrt Y\) and \(Y\) contribute
\(O_r((Y/\log Y)\sum_{j\ge0}(j+1)2^{-j})\)
by the prime count. The bound \(v!\le v^v\) proves (675.S)
when \(v>0\), and the case \(v=0\) is immediate.
The factor \(v!\) counts injections of labelled sites into each
set of \(v\) primes.

A lower bound for \(\delta_{r,N}\) follows by setting \(t_q=0\)
for \(q\le\beta Y\) and \(t_q\ne0\) for \(\beta Y<q\le Y\),
using the safe nonzero choices and the injections into
\((LY,2LY]\) from the proof of (675.P). Only sites \(p^r\)
for \(\beta Y<p\le Y\) remain to be covered. If their number
is \(m\), that construction has probability at least
\[
 \exp\bigl(-r\vartheta(\beta Y)
    -(r-1)m\log Y-m\log\log Y-O_r(Y/\log Y)\bigr).
\]
This uses \(M\asymp Y/\log Y\) in the lower count of injections
and keeps the \(\log\log Y\) term. With
\(m=\pi(Y)-\pi(\beta Y)\), the prime number theorem with its
classical logarithmic error gives the required upper bound
on the negative logarithm in (675.R).

### Averaging squarefree witness collisions

For \(r=2\), take
\[
             \eta=\frac{(\log\log Y)^{1/4}}{\log Y}.
\]
Among the sites \(p^2,2p^2\le Y^2\) with prime \(p\ge\eta Y\),
the number of primes \(q\ge\eta Y\) that can cover both
\(p^2\) and \(2p'^2\) from different rows is
\[
                \ll\frac{Y(\log\log Y)^{1/4}}{\log^2 Y}.
                                                               \tag{675.T}
\]
Such primes satisfy \(p^2-2p'^2=kq^2\) with
\(1\le |k|\le K=\lceil\eta^{-2}\rceil\). The three primes
are distinct for large \(Y\), since equality with a row prime
would imply equality of the row primes. We count prime triples
over all these \(k\) using the norm form for \(\mathbb Q(\sqrt2)\).

The ring \(\mathbb Z[\sqrt2]\) is norm Euclidean. Rounding
both rational coordinates of a quotient to integers gives
an error of absolute norm at most \(1/2\). Since \(q\nmid 2k\)
and \(p,p'\) are distinct from \(q\), the two conjugate primes
above \(q\) cannot both divide \(p+p'\sqrt2\). Indeed an inert
\(q\), or both primes above a split \(q\), would divide both
integer coordinates. Thus \(q\) splits and we can write
\[
 p+p'\sqrt2 = \alpha\gamma^2,\qquad
 N(\alpha)=k,\quad N(\gamma)=\pm q ,
\]
where \(\gamma\) generates the prime whose square divides
\(p+p'\sqrt2\). Multiply \(\gamma\) by a suitable power of
the unit \(1+\sqrt2\) to arrange that both absolute embeddings
of \(\alpha\) are between fixed positive multiples of
\(\sqrt{|k|}\). Then each embedding of \(\gamma\) is
\(O(\sqrt Y/|k|^{1/4})\). For each ideal generated by
\(\alpha\), the number of such balanced generators is bounded,
since their ratios are algebraic integers with both embeddings
in a fixed bounded range. The number of ideals of norm \(n\)
here is \(a(n)=\sum_{d\mid n}\chi_8(d)\), by factoring rational
prime ideals. Thus there are \(O(a(|k|))\) possible \(\alpha\).

Write \(\alpha=A+B\sqrt2\) and \(\gamma=x+y\sqrt2\).
The three prime forms, with either sign for the last one, are
\[
 A x^2+4Bxy+2A y^2,\qquad Bx^2+2Axy+2By^2,\qquad x^2-2y^2 .
\]
Their discriminants are \(-8k,4k,8\). At a prime
\(\ell\nmid 2k\), the first two forms have no common
projective zero because their coefficients in \(x^2+2y^2,2xy\)
form a matrix of determinant \(k\). The norm identity
\(F_1^2-2F_2^2=kF_3^2\) then excludes common zeros involving
the third form as well. Hence the local density of zeros of
their product is
\[
 g(\ell)=\frac{1+\rho(\ell)(\ell-1)}{\ell^2},
 \qquad \rho(\ell)=3+\chi_{-8k}(\ell)+\chi_{4k}(\ell)+\chi_8(\ell).
\]

Apply the same Selberg upper sieve as for the fixed conics,
now to all three forms in the box
\(|x|,|y|\ll\sqrt Y/|k|^{1/4}\).
Use primes at most \(D=Y^{1/100}\) other than the primes
dividing \(2k\) or at most \(13\). None can divide a prime
coordinate of a counted triple, since those primes are at
least \(\eta Y>D\). The lattice counting error at a modulus
\(d\) is \(O(\sqrt Y\,d+d^2)\). The resulting total error
per \(\alpha\) is \(O(Y^{1/2}D^4+D^6)\), as in the fixed
conic estimate.

For a nonzero integer \(u\), put \(L_u=L(1,\chi)\), where
\(\chi\) is the primitive quadratic character of the square
class of \(u\), and put \(L_u=1\) in the principal case.
Let \(C(k)=\prod_{\ell\mid 2k}(1-1/\ell)^{-6}\).
With \(z=D^\gamma\) for a sufficiently small fixed
\(\gamma>0\), the product of \((1-g(\ell))^{-1}\) over
usable \(\ell\le z\) is
\[
                 \gg (\log Y)^3 L_{-2k} L_k/C(k).
\]
In fact \(\log(1-g(\ell))^{-1}=\rho(\ell)/\ell+O(\ell^{-2})\).
The contribution \(3/\ell\) gives the cubic logarithm.
For any nonprincipal character in the display for \(\rho\),
the sum of \(\chi(\ell)/\ell\) through \(z\) equals
\(\log L(1,\chi)+O(1)\) uniformly: its conductor is
\(O(K)\), a fixed power of \(\log z\), and partial summation
of the Siegel–Walfisz estimate bounds the omitted prime tail.
The fixed \(\chi_8\) contributes a bounded term. A principal
character adds a positive logarithm and may be dropped in
this lower bound. Omitting primes dividing \(2k\) loses at
most \(C(k)\) up to an absolute multiplicative constant.
The normalization of the expanded product as a probability
measure, with \(\mathbb E\log d\ll\log z\), gives the same
lower bound up to a constant for \(G(D)\), using
\(\gamma\) small. The main term of the sieve, summed over
\(\alpha\), is therefore
\[
        \ll\frac{Y}{\log^3 Y}
               \frac{a(|k|)C(k)}{\sqrt{|k|}L_{-2k}L_k}.
\]
[Granville and Soundararajan, Theorem 12.1 on Siegel–Walfisz](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf)

The coefficient in this bound satisfies
\[
 \sum_{1\le |k|\le K}
      \frac{a(|k|)C(k)}{\sqrt{|k|}L_{-2k}L_k}
                   \ll \sqrt K\sqrt{\log(2K)}.                    \tag{675.U}
\]
For the first factor in Cauchy–Schwarz,
\(\sum_{n\le K}a(n)^2\ll K\log(2K)\).
At prime powers, \(a(n)^2\le(a*a)(n)\):
at a split prime this is
\((e+1)^2\le\binom{e+3}{3}\), and at inert and ramified
primes it follows from the local ideal counts.
Also \(\sum_{n\le K}a(n)\ll K\), by summing the bounded
partial sums of \(\chi_8\) in the divisor formula.
Summing the convolution proves the second-moment bound.
For the other factor,
\[
           \sum_{1\le |k|\le K} C(k)^2 L_{-2k}^{-2}L_k^{-2}\ll K .
\]
Hölder's inequality uses fixed sixth moments of
\(C(k),L_{-2k}^{-1},L_k^{-1}\). The moments of \(C(k)\)
are bounded on average by expanding its Euler factors
as a sum over squarefree divisors; at each prime the
expanded coefficient is \(O(1/\ell)\) and occurs with
frequency at most \(1/\ell\). For either inverse \(L\)-value,
Granville and Soundararajan give bounded fixed negative
moments averaged over fundamental discriminants, including
the exceptional characters by their remark after Theorem 2.
For arbitrary \(u\) up to \(2K\) in absolute value, write
\(|u|\) as a square times its squarefree part. At a fixed
square factor \(h^2\), the corresponding fundamental
discriminants are at most \(8K/h^2\) in absolute value,
with bounded multiplicity. Summing over \(h\) proves the
mean bound for \(L_k\) and \(L_{-2k}\); principal cases give
one by convention. Cauchy–Schwarz in dyadic intervals and
summation with the weight \(1/\sqrt{|k|}\) prove (675.U).
[Granville and Soundararajan, *The distribution of values of \(L(1,\chi_d)\)*, Theorem 2 and its following remark](https://dms.umontreal.ca/~andrew/PDF/L1chi2002.pdf)

The summed lattice errors are negligible since
\(K=O((\log Y)^2)\) and \(\sum_{n\le K}a(n)\ll K\).
Using (675.U), the main terms total
\(O(Y\sqrt K\sqrt{\log(2K)}/\log^3Y)\).
For the chosen \(\eta\), this is the bound in (675.T).

### Summing the prime choices

Use the rows for \(p\ge\eta Y\) and just the sites \(p^r\)
and \(2p^r\) when present. The small forced prime powers
cover none of these rows. Each prime above the forcing
threshold covers at most one site of either form, by the
cyclotomic argument. Thus primes below \(\eta Y\) can
touch at most \(2\pi(\eta Y)\) rows. Among larger primes,
discard rows touched by cross-row collisions. For \(r=2\)
use (675.T). For \(r\ge3\) take \(\eta=(\log Y)^{-2}\):
the uniform plane curve bound gives
\(O_r(\eta^{-r}Y^{2/r+\nu})=o(Y/\log^2Y)\)
such primes, with \(0<\nu<1-2/r\).
Omitted and discarded rows thus cost at most
\(O_r(Y(\log\log Y)^{1/4}/\log Y)\) in the negative log
probability bound, or \(O_r(Y/\log Y)\) for \(r\ge3\).

Specify the discarded rows and the primes with
\(p^r\mid t\) on the other rows in \(\exp(O(Y/\log Y))\)
ways. The latter cost \(r\log p\) each. Replacing these
costs by \(r\log Y\) loses at most \(O_r(Y/\log Y)\)
because \(\sum_{p\le Y}\log(Y/p)\ll Y/\log Y\).
On the retained nonzero rows every site needs a separate
witness prime, distinct from the prescribed zero primes.
Thus (675.S) gives a negative log cost, up to the total
error \(O_r(Y/\log Y)\), of \((r-1)\log Y+\log\log Y\)
per site. If \(p\le\beta Y\) there are two sites and
the smaller row cost is \(r\log Y\). Above \(\beta Y\)
there is one site and the smaller cost is
\((r-1)\log Y+\log\log Y\).
Summing and applying the prime number theorem proves
the lower bound on the negative logarithm in (675.R).

For \(r=2\), any retained nonzero row at \(p\le\beta Y\)
costs an extra \(2\log\log Y\). The number of discarded
or omitted rows is \(o(Y/\log Y)\). Requiring at least
\(\lambda Y/\log Y\) such primes therefore adds
\((2\lambda+o(1))Y\log\log Y/\log Y\) to the negative log
bound. Dividing by the total density gives the asserted
proportion.

In particular, write \(M_\beta=\prod_{p\le\beta Y}p^2\)
when \(r=2\). For each fixed \(\varepsilon>0\), the proportion
of returns with
\(\log\gcd(t,M_\beta)\le(\sqrt2-\varepsilon)Y\)
is at most \(\exp(-c_\varepsilon Y\log\log Y/\log Y)\) for some \(c_\varepsilon>0\).
Indeed \(\log M_\beta=(\sqrt2+o(1))Y\), and each omitted
prime power removes at most \(2\log Y\).

For the finite-height lower bound use the profile with cutoff
\(\beta\), retaining all its safe residues and injections into
\((LY,2LY]\) and summing the progression counts as in (675.P).
This supplies the same second term with every negative witness
at most \(B_rY=2LY\). For the upper bound with witness primes
at most \(Q\), specify the discarded rows, the zero choices
on retained rows, and the witness injection for the remaining
sites. Each specification defines one CRT class. Its count
through \(X\) is at most its probability times \(X\) plus one.
There are \(\exp(O_r((Y/\log Y)(1+\log Q)))\) specifications,
by the bound \(O_r(Y/\log Y)\) on the number of labelled
sites. The probability sum is bounded by the same
argument using (675.S). For \(Q=B_rY\), the summed floor error
is absorbed at \(X\ge\exp(C_rY)\) with \(C_r\) large.

For unrestricted returns use \(Q=\exp(c_rY)\), with \(c_r\)
large. The shifts with any site covered by a prime \(q>Q\)
number at most
\[
              O_r\bigl(NXQ^{1-r}+N(X+N)^{1/r}\bigr)
\]
by a union bound over offsets and primes. This is negligible
compared with the lower count. At this \(Q\), the floor error
is absorbed when \(X\ge\exp(C'_rY^2/\log Y)\) with
\(C'_r\) large. Both finite-height conclusions require
intervals above the scale of the first unrestricted squarefree
return, whose order remains open.

## The common divisor of the returns

Let \(\mathcal R_N\) be the set of positive squarefree prefix returns and
let \(G_N(X)=\gcd\{t\in\mathcal R_N:t\le X\}\) whenever this set is
nonempty. Write \(G_N=\gcd\mathcal R_N\), and let
\[
 A_p=\{c\bmod p^2: a+c\not\equiv0\bmod p^2
                    \text{ for every }a\in S_2\cap[1,N]\},
 \qquad j_p=\min_{c\in A_p}\min(2,v_p(c)),
\]
with the truncated valuation of the zero class taken as 2. The zero
class belongs to \(A_p\).

**Theorem (Common divisor of short squarefree prefix returns).**
There is an absolute \(C\) such that for all sufficiently large \(N\),
\[
             G_N(\exp(C\sqrt N))=G_N=\prod_p p^{j_p}.
\]
In particular \(G_N\) is determined by the missing residue classes
of the positive offsets. For each \(p\le \sqrt N\), every class
in \(A_p\) occurs among the returns up to \(\exp(C\sqrt N)\).

Every return lies in \(A_p\) at every prime by definition. Qualitatively,
any finite assignment of classes in these sets extends to a
positive-density family of returns. Complete the assignments by
fixing the zero class at all remaining \(q\le\sqrt N\). For each
negative site still lacking a prime-square divisor, select a
distinct unused prime \(q>\sqrt N\) and fix its residue to cover
that site. It cannot cover another offset. At the rest of the
primes \(q>\sqrt N\), the fraction of residues allowed by the
positive sites is \(1-|S_2\cap[1,N]|/q^2\). These fractions have
a positive product. Independence and the fixed-\(N\) density
argument above give the claimed completion. Prescribing a lift
modulo \(p^3\) of a zero class also works by the same finite CRT
and tail argument. This gives the formula for \(G_N\) for every
fixed \(N\).

For the finite-height assertion write \(H=\sqrt N\) and
\(M_0=\prod_{q\le H}q^2\). Fix \(p\le H\) and \(c\in A_p\).
Prescribe \(t=c\bmod p^2\) and \(t=0\bmod q^2\) at the other
primes \(q\le H\). The positive sites are preserved by these
tests. Any negative site still lacking a witness is a multiple
of \(p^2\). If \(c=0\) there are none. If \(c\ne0\), the
squarefree witness argument in the density upper bound shows
that \(p>N^{3/10}\) for all sufficiently large \(N\). Thus
there are at most \(N^{2/5}\) sites still to cover. The primes
in \((2H,4H]\) number \(\gg H/\log H\) by the prime number
theorem. Choose one distinct prime there per uncovered site,
as in the qualitative argument. The resulting progression
\(b+M u\) has
\[
   \log M\le 2\vartheta(H)+2N^{2/5}\log(4H)=(2+o(1))H.
\]
All the prescribed tests preserve the pattern at the sites they
can affect, and every negative site now has a witness. The
progression version of the uniform count gives a return of
size at most \(\exp(C H)\), with \(C\) independent of \(N,p,c\).

It remains to check the extra prime factors of the finite gcd.
Use the baseline progression \(M_0u\) with
\(U=\lceil\exp(C_0 H)\rceil\), where \(C_0\) is large enough
for the uniform count. For any prime \(p\), apply the finite
Bonferroni upper count within \(u=0\bmod p\), omitting \(p\)
from the finite sieve if present. The remaining moduli are
coprime to \(p\). The expansion has the same uniform absolute
error \(O(E)\), with the main term multiplied by \(1/p\).
The finite product \(V_0\) is unchanged if \(p\) is absent,
and is divided by \(1-k/p^2\) otherwise. Consequently the
proportion of the full baseline returns with \(p\mid u\) is
at most
\[
 \begin{cases}
  (1+o(1))/p+o(1),&p\le H,\\[2pt]
  O(1/p)+o(1),&p>H.
 \end{cases}
\]
The errors are uniform in \(p\), by the same CRT remainder
bound and \(E=o(UV_0)\). In particular for every prime there
is a baseline return with \(p\nmid u\), all at the same height
bound \(\exp(C H)\) after increasing \(C\).
For \(p>H\) this excludes \(p\) from the finite gcd. For
\(p\le H\) it bounds the valuation by 2, and the earlier
prescriptions realize the minimum if \(j_p<2\).

**Corollary (Size limit for a common divisor of returns).**
As \(N\to\infty\),
\[
 \log G_N\ \le\
   \left(2\sqrt{2/\zeta(2)}+o(1)\right)
   \sqrt{\frac{N\log\log N}{\log N}} .
\]
The same bound holds for \(G_N(\exp(C\sqrt N))\) with the
constant in the theorem. This is \(o(\sqrt N)\), so a lower
bound obtained solely by exhibiting a divisor common to every
return cannot give the conjectural square-root scale for the
logarithm of the least return.

To prove the bound, fix a small \(\varepsilon>0\). For every
sufficiently large prime \(p\) there is a reduced class
\(a\bmod p^2\) whose first
\[
    L=\left\lfloor (\zeta(2)-\varepsilon)
                         \frac{\log p}{\log\log p}\right\rfloor
\]
positive terms are all nonsquarefree. To see this, use
indices \(0\le j<L\) for the terms \(a+jp^2\), where
\(1\le a\le p^2\). At all primes \(q\le B\), with \(B\)
fixed, impose \(a=0\bmod q^2\). These cover the indices
divisible by some \(q^2\), since \(p>B\). The number \(m\)
of remaining indices is
\(L\prod_{q\le B}(1-q^{-2})+O_B(1)\) by finite periodicity.
Assign a different prime \(q>B\), \(q\ne p\), to each
such index \(j\), prescribing \(a=-jp^2\bmod q^2\).
Use the \(m\) smallest available primes. Their product
together with the small primes, all squared, gives a
modulus \(D\) with
\[
        \log D=(2+o_B(1))m\log m.
\]
This follows from the prime number theorem, and the largest
prime used is \(O_B(L\log L)\ll p\). Choose \(B\) large
enough that
\((\zeta(2)-\varepsilon)\prod_{q\le B}(1-q^{-2})<1\).
There is then a fixed positive margin in \(\log D<2\log p\),
so \(2D<p^2\) for all sufficiently large \(p\).
Each prescription therefore has at least two consecutive
representatives in its residue class modulo \(D\) within
\([1,p^2]\). They cannot both be divisible by \(p\).
This adaptation of Erdős's squarefree-gap construction supplies
the reduced class. [Erdős, equations (3)–(5), p. 176](https://www.renyi.hu/~p_erdos/1966-14.pdf)

For any fixed \(b>\sqrt{2/\zeta(2)}\), take
\(p>b\sqrt{N\log\log N/\log N}\) and \(p\le H\).
Choosing \(\varepsilon\) small enough in terms of \(b\) gives
\(L p^2\ge N\) uniformly in this interval. Thus no positive site
lies in the constructed reduced class \(a\bmod p^2\).
The unit residue \(-a\) belongs to \(A_p\), implying \(j_p=0\).
Primes beyond \(H\) were already excluded by the finite gcd
argument. Therefore
\(\log G_N\le 2\vartheta(b\sqrt{N\log\log N/\log N})\).
Apply the prime number theorem and let \(b\) decrease to
\(\sqrt{2/\zeta(2)}\).

The [sums-of-two-squares problem](ERDOS_675_TRANSLATION_PATTERNS.md)
has a different prime obstruction and remains open.
