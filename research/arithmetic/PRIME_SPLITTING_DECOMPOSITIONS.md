# Actual prime decompositions with splitting restrictions

Prime indicators are not multiplicative Euler indicators. Restricting a
semigroup count to its prime members therefore needs a prime-correlation
theorem, not just the previous multiplicative-set asymptotic.

Here the summands themselves are prime. Classical inputs give full-scale
leading constants, natural-order counts in power-short windows, and
prime-value decompositions for a class-number-three quadratic form.
The every-target two-prime problem remains a separate obligation.
All results are recorded as classical consequences with explicit
application proofs, without an analytic novelty claim.

The statements and primary-source applications have internal written
review. No external human review, formal verification, numerical onset
or complete small-exception classification is asserted.

## 1. Prime stars and their exact local factors

Fix even \(L\ge4\), and let
\[
 {\cal P}_L=\{p:p\text{ prime},\ p=1\bmod L\}.
\]
For fixed \(r\ge1\), count ordered prime stars
\[
 a_j+b_j+c=N_j\quad(1\le j\le r),\qquad
                 c,a_j,b_j\in{\cal P}_L.                 \tag{P.1}
\]
Repetition is initially permitted. Put \(k=2r+1,d=r+1\), and
\[
 m=\min_jN_j,\qquad
 {\cal V}_r(\boldsymbol N)=\int_0^m\prod_{j=1}^r(N_j-t)\,dt.
\]
For \(p\nmid L\), define
\[
 \beta_{p,r}(\boldsymbol N)=
 \frac{p^r}{(p-1)^{2r+1}}
   \sum_{z\in\mathbb F_p^\times}
             \prod_{j=1}^r
                    (p-2+\mathbf1_{z=N_j\bmod p}).        \tag{P.2}
\]
These are ordinary prime unit-root factors at EVERY prime outside
\(L\), including primes in the allowed splitting classes.

**Theorem PrimeStarAsymptotic.** Fix \(C\ge1\). Uniformly for
\(N_j\in[T,CT]\), \(N_j=3\bmod L\),
\[
 R_{r,L}(\boldsymbol N)=
 \frac{L^r}{\varphi(L)^{2r+1}}\,
 {\cal V}_r(\boldsymbol N)
 \frac{\prod_{p\nmid L}\beta_{p,r}(\boldsymbol N)}
                         {(\log T)^{2r+1}}
      +o\!\left(\frac{T^{r+1}}{(\log T)^{2r+1}}\right).    \tag{P.3}
\]
The product and its inverse are uniformly bounded. Every such large
target has a representation. Other target residue classes have zero
count. Requiring all \(2r+1\) primes to be distinct gives the same main
term.

**Proof.** Use Green–Tao, *Linear equations in primes*, Annals 171
(2010), Main Theorem and Corollary 1.7, p. 1763.
The uniform formulations are (1.7)–(1.8), pp. 1757–1759, and
Theorem 1.8, pp. 1765–1766. They allow target-sized moving affine
constants and convex regions; the printed unweighted conversion is
part of the stated scope.

Take free variables \(c,a_1,\ldots,a_r\), with
\(b_j=N_j-c-a_j\). This system has Cauchy–Schwarz complexity at most
one, independently of \(r\). For \(c\), partition the other forms
into all \(a_j\) and all \(b_j\); neither span contains its normal.
For \(a_j\), use \(\{c,a_i:i\ne j\}\) and all \(b_i\).
For \(b_j\), use all \(a_i\) and \(\{c,b_i:i\ne j\}\).
Again neither normal span contains the omitted normal.
Thus the unconditional low-complexity corollary applies.
Fixed progressions are imposed by substituting the free prime
variables \(Lx+1\); all coefficients remain bounded.

At modulus \(L\), the free parameters have exactly one compatible
state. Each prime progression has logarithmic mean \(L/\varphi(L)\).
Their product times the lattice density \(L^{-d}\) gives
\(L^{k-d}/\varphi(L)^k=L^r/\varphi(L)^k\).
At a good prime all prime variables must be units.
Conditioning on the shared unit \(z\), a private unit pair has
\(p-2+\mathbf1_{z=N_j}\) solutions, giving (P.2).
All factors are positive since \(p\ge3\).
The normals of every pair are independent modulo every prime;
inclusion-exclusion yields \(\beta_{p,r}=1+O_r(p^{-2})\) uniformly.

The real volume is the stated integral.
The source's uniform convex-body result supplies the full comparable
target range. In a direct deweighting proof, first restrict all prime
values to \([\varepsilon T,CT]\), where logarithmic weights are
uniformly \((\log T)^k(1+o(1))\). Strip contributions have
\(O(\varepsilon T^d)\) weighted main volume. Values below \(T^\alpha\)
occupy \(O(T^{d-1+\alpha})\) lattice points; the rest of the strip
is bounded using weights at least \(\alpha^k(\log T)^k\).
Let \(T\) tend to infinity and then \(\varepsilon\) to zero.
Prime powers, if using \(\Lambda\), occupy at most
\(O(T^{d-1/2}(\log T)^{O_r(1)})\) weighted mass and are negligible.
Finally, any equality between two prime variables is a proper affine
hyperplane, with \(O(T^{d-1})\) points. Remove their finite union.
\(\square\)

For \(L=12,r=1\), write \(c_p(N)=p-1\) if \(p\mid N\), and \(-1\)
otherwise. The exact triple formula is
\[
 \#\{p_1+p_2+p_3=N:p_i=1\bmod12\}
 \sim\frac3{32}\frac{N^2}{(\log N)^3}
        \prod_{p\ge5}\left(1-\frac{c_p(N)}{(p-1)^3}\right),
 \qquad N=3\bmod12.                                      \tag{P.4}
\]
For the two-equation shared five-prime system, \(L=12,r=2\),
the coefficient is \(9/64\) times
\[
 {\cal V}_2=NMm-\tfrac12(N+M)m^2+\tfrac13m^3.
\]
Its local factor is
\[
 \frac{p^2}{(p-1)^5}
 \left[(p-1)(p-2)^2+(p-2)
       (\mathbf1_{p\nmid N}+\mathbf1_{p\nmid M})
                    +\mathbf1_{N=M\ne0\bmod p}\right]     \tag{P.5}
\]
at all \(p\ge5\). These constants differ from the semigroup constants.

## 2. Every-target natural-order prime counts in short windows

Fix \(2/3<\theta<1\), \(H=T^\theta\), and \(N_j\in[T,CT]\) with
\(N_j=3\bmod L\). Set
\[
 c_0=\tfrac13\min_jN_j,\qquad t_j=\tfrac12(N_j-c_0).
\]
All these centers are comparable to \(T\).

**Theorem PrimeStarWindows.** The number of ordered solutions to (P.1)
with
\[
 |c-c_0|\le H,\qquad |a_j-t_j|,|b_j-t_j|\le H
\]
satisfies
\[
 R^{\mathrm{win}}_{r,L}(\boldsymbol N;H)
                  \asymp_{L,r,C,\theta}
                     \frac{H^{r+1}}{(\log T)^{2r+1}}.     \tag{P.6}
\]
The bounds are uniform in all admissible target valuations.
They also hold with every prime distinct.
No leading constant in these shrinking windows is asserted.

**Proof: lower bound.** The input is Halupczok, *Goldbach's problem
with primes in arithmetic progressions and in short intervals*,
Theorem 3.1, pp. 339–340, with the singular series (1.1)–(1.2), p. 335.
The previously checked [Hq/Tq parameter applications](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#an-every-target-ternary-baseline)
are used with fixed \(L\).
Choose its \(\varepsilon\) with \(0<\varepsilon<\theta-2/3\);
all size inequalities then hold uniformly at centers comparable to
\(T\). Put \(Y=H/10\).

First we need \(\gg Y/\log T\) actual primes \(c\) in
\((c_0,c_0+Y]\), \(c=1\bmod L\).
The first Halupczok application, at \(X_1=X_2=c_0\), with
\(a_1=2,a_2=1,Q_1=Q_2=L\), supplies
aggregate binary logarithmic weight \(\gg Y^2\) over a target block
of length \(2Y\).
For a fixed candidate prime, its partner varies over an interval of
length \(2Y\). Montgomery–Vaughan, *The large sieve*, Theorem 2,
p. 121, equation (1.12), bounds the number of partner primes by
\(O(Y/\log Y)\). Each pair weighs \(O((\log T)^2)\).
A candidate therefore contributes \(O(Y\log T)\), and there are
\(\gg Y/\log T\) candidates.

For each \(j\), apply Halupczok again with
\[
 X_1=t_j-Y/2,\qquad X_2=t_j,\qquad
 a_1=2,\ a_2=1,\ Q_1=Q_2=L.
\]
The source condition \(X_2\le2X_1-Y\) is \(t_j\ge2Y\);
the other hypotheses follow uniformly from \(\theta>2/3\).
Each candidate remainder \(N_j-c\) lies in
\((2t_j-Y,2t_j)\), except a possible endpoint candidate which we discard.
The map from \(c\) to the remainder is injective.

For any fixed \(A\), only \(O(Y/(\log T)^A)\) remainders in this
block have logarithmic binary count below a fixed multiple of \(Y\).
The singular series has a uniform positive lower bound.
The union of the \(r\) bad sets cannot exhaust the
\(\gg Y/\log T\) candidates. Every remaining \(c\) has
\(\gg Y/(\log T)^2\) private prime pairs for EACH equation.
Choices of private pairs are independent once \(c\) is fixed.
Multiplying gives \(\gg Y^{r+1}/(\log T)^{2r+1}\) stars.
The private primes lie within \(2Y\) of their centers. This proves
the required lower bound.

**Proof: upper bound.** We give the sieve estimate to justify the
matching logarithmic power. Use the \(d=r+1\) free progression
variables in a box of side \(O_L(H)\), and let \(F\) be the product
of the \(k=2r+1\) prime forms. Every counted prime is comparable to
\(T\), so it exceeds \(D=H^\eta\), with \(\eta>0\) chosen small.

For a good prime, let \(\nu(p)\) count roots of \(F=0\) in
\(\mathbb F_p^d\), and put \(\rho(p)=\nu(p)/p^d\).
The \(k\) distinct hyperplanes and independent pair normals give
\[
 k/p-\binom{k}{2}/p^2\le\rho(p)\le k/p.
\]
Their complement is nonempty for every odd prime, by (P.2).
For squarefree \(m\), CRT and box counting give
\[
 \#\{F=0\bmod m\}
  =\operatorname{vol}(\text{box})\rho(m)
    +O_{L,r}\!\left(k^{\omega(m)}
              \sum_{h=0}^{d-1}H^h m^{d-1-h}\right).       \tag{P.7}
\]
Here \(\omega(m)\) in the error denotes the number of prime divisors.

For completeness the Selberg quadratic minimization is algebraic.
For squarefree-supported weights \(\lambda_1=1,\lambda_m=0\) at
\(m>D\), put \(g(p)=\rho(p)/(1-\rho(p))\) and
\(G(D)=\sum_{m\le D}\mu(m)^2g(m)\).
Diagonalizing
\(\sum_{a,b}\lambda_a\lambda_b\rho([a,b])\) through common divisors,
or Cauchy–Schwarz with variables
\(y_h=\sum_{h\mid a}\lambda_a\rho(a)\), gives minimum \(1/G(D)\).
Write \(G_m(x)=\sum_{n\le x,(n,m)=1}\mu(n)^2g(n)\).
The minimizing weights satisfy
\[
 \lambda_m=\mu(m)\prod_{p\mid m}(1-\rho(p))^{-1}
                       \frac{G_m(D/m)}{G(D)},\qquad
 |\lambda_m|\ll_{\epsilon,L,r}m^\epsilon .
\]
The weight bound is uniform in the targets; finitely many small
prime factors are bounded using the positive complement in (P.2).

Uniformly \(G(D)\gg_{L,r}(\log D)^k\).
To see this without a moving-target asymptotic, omit finitely many
primes including those dividing \(L\), and use the lower baseline
\(\rho_0(p)=k/p-\binom{k}{2}/p^2>0\).
Its Euler subset mass through \(z=D^{1/(8k)}\) is
\(\prod_{p\le z}(1+g_0(p))\asymp(\log z)^k\).
In the normalized subset law,
\(\mathbb E\log m=\sum_{p\le z}\rho_0(p)\log p
 \le2k\log z+O(1)\).
Markov leaves a fixed positive proportion of this mass at \(m\le D\).
This proves the asserted lower bound for \(G\).

Insert the optimal weights in the upper inequality
\(\mathbf1_{(F,\prod_{p\le D,p\nmid L}p)=1}
 \le(\sum_{m\mid F,m\le D}\lambda_m)^2\).
The main term is \(O(H^d/G(D))\).
Using (P.7), \( [a,b]\le D^2\), and
\(k^{\omega(m)}\ll_\epsilon m^\epsilon\), the total error is
\[
 O_{\epsilon,L,r}
          (H^{d-1}D^{2+O(\epsilon)}+D^{2d+O(\epsilon)}).
\]
Take \(0<\eta<1/(4d+4)\) and then fixed small \(\epsilon\);
this is power smaller than \(H^d\).
Since \(\log H\asymp\log T\), (P.6) follows.

Finally all equality conditions between two prime forms are proper
hyperplanes, containing \(O(H^{d-1})\) parameter points.
Their union is negligible compared with the lower bound.
\(\square\)

For \(r=1\), this proves the natural order
\(N^{2\theta}/(\log N)^3\) for almost-equal prescribed-AP prime
triples on EVERY large admissible target, strengthening the earlier
existence-only Tq statement. For \(r=2\), it gives
\(H^3/(\log T)^5\) shared five-prime decompositions.
The proof is a classical Halupczok/Brun–Titchmarsh/sieve deduction.

## 3. Prime norm generators beyond class number one

Let \(Q_0(x,y)=x^2+xy+6y^2\), of discriminant \(-23\), and let
\[
 {\cal P}_{-23}=\{p\ne23:p\text{ prime and }p=Q_0(x,y)
                                      \text{ for some }x,y\in\mathbb Z\}.
\]
Primitivity is automatic at a prime value.
The existing [class-number-three argument](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md)
proves \(\operatorname{Cl}(\mathbb Q(\sqrt{-23}))=C_3\).

**Theorem PrincipalNormPrimeTriples.** Every sufficiently large odd
integer is a sum of three pairwise-distinct primes in
\({\cal P}_{-23}\). Write \(\chi(N)=\left(\frac N{23}\right)\)
for the quadratic character. More precisely, let
\[
 B_{23}(N)=
 \begin{cases}
 66,&23\mid N,\\
 55,&\chi(N)=1,\\
 60,&\chi(N)=-1.
 \end{cases}
\]
The ordered count satisfies, uniformly over odd \(N\),
\[
 R_{-23}^{\mathbb P}(N)\sim
 \frac{23B_{23}(N)}{66^3}\frac{N^2}{(\log N)^3}
 \prod_{p\ne2,23}
                  \left(1-\frac{c_p(N)}{(p-1)^3}\right).  \tag{P.8}
\]
There is also an eventual coupled version: every sufficiently large
comparable odd \(N,M\) has
\(a+b+c=N,\ c+d+e=M\), with all five distinct prime values in
\({\cal P}_{-23}\).

**Proof.** Let \(K=\mathbb Q(\sqrt{-23})\) and let \(F\) be its
Hilbert class field. Class field theory gives
\(\operatorname{Gal}(F/K)\cong C_3\), with principal prime ideals
exactly those splitting in \(F\). Complex conjugation acts by
inversion, so \(\operatorname{Gal}(F/\mathbb Q)\cong S_3\) and its
maximal abelian subfield is \(K\).
For \(p\ne23\), a prime value of \(Q_0\) is exactly a rational prime
splitting completely in \(F\). A principal prime ideal has a generator
of norm \(p\), and conversely such a generator makes its ideal principal.
The exact class-field source is Milne, *Class Field Theory*, v4.03,
Theorem 0.4 and following Hilbert-class-field paragraph, with full
Artin formulation in Chapter V, Theorems 3.5–3.6 and Example 3.9.
There is no element of norm two, so every prime in this set is odd.

Use Frei–Tinková, *Linear constellations in primes with arithmetic
restrictions*, arXiv:2408.17441v3, Theorem 1.4, pp. 4–5, and
Remark 1.5. Its Chebotarev-prime theorem is unconditional;
the separate primitive-root theorem has different assumptions.
Kane, *An Asymptotic for the Number of Solutions to Linear Equations in
Prime Numbers from Specified Chebotarev Classes*, Theorem 2 in the
retained arXiv:1009.1368v5 edition, already covers the single linear
triple equation; the later theorem supplies the coupled-system scope.
The Frei–Tinková theorem allows moving \(O(T)\) intercepts and convex bodies.
Prescribe the identity Artin class in \(F/\mathbb Q\).
Its density is \(1/6\), and positivity depends only on the maximal
abelian restriction.

The compatible nonzero states modulo 23 are the eleven quadratic
residues \(R\). Three copies of \(R\) cover \(\mathbb F_{23}\).
Explicitly the private-pair count is zero at zero, five at a nonzero
square and six at a nonsquare. Summing over the third residue gives
\(B_{23}\) as displayed. These formulas follow by expanding the
indicator \((1+\chi(x)-\mathbf1_{x=0})/2\).

The intersection \(F\cap\mathbb Q(\zeta_{23})\) is \(K\), since
\(F\)'s maximal abelian subfield is \(K\). Thus the actual Artin
density in each compatible progression is \(1/66\).
The finite factor at 23 is \(23B_{23}/66^3\); the odd-prime condition
contributes factor two at 2, canceling the simplex factor \(1/2\).
All other primes have the ordinary triple unit factors.
This gives (P.8) and uniform positivity.

For the coupled five-prime system, choose a shared residue in \(R\)
different from \(N,M\bmod23\); at least nine choices remain.
Each private pair then has five or six residue representations.
The ordinary odd-prime factors are positive as in (P.2).
The same theorem proves the comparable-target coupled assertion.
Distinctness removes only lower-dimensional equality loci.
\(\square\)

This is an explicit nonabelian Chebotarev-prime consequence.
It is not obtained from Dirichlet characters or from the earlier
class-group semigroup lower bound alone. No almost-equal power-window
claim for \({\cal P}_{-23}\) is inferred from the fixed-volume theorem.

## 4. The exact remaining binary gap

Prime triples do not settle
\[
             N=p+q,\qquad p,q=1\bmod12.                   \tag{P.9}
\]
The existing [Hq theorem](BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md#a-known-density-one-range-in-the-congruence-family)
gives actual prime pairs for almost all admissible targets when
\(\theta>2/3\), with arbitrarily strong logarithmic exception savings.
It does not eliminate every exception.
For (P.9), the expected leading count on \(N=2\bmod12\) is
\[
 \frac34\frac{N}{(\log N)^2}
       \prod_{p\ge5}
             \left(1+\frac{c_p(N)}{(p-1)^2}\right),        \tag{P.10}
\]
where the plus sign uses \(c_p=p-1\) on divisors and \(-1\) otherwise.
This is a prediction, not a theorem established here.

**Theorem PrimeBasisDichotomy.** For fixed even \(L\ge4\), the
asymptotic basis order of \(\{0\}\cup{\cal P}_L\) is either \(L+1\)
or \(L+2\). It is \(L+1\) if and only if every sufficiently large
\(N=2\bmod L\) is the sum of two primes in \({\cal P}_L\).
In particular, the prime \(1\bmod12\) basis order is 13 or 14;
order 13 is equivalent to eventual restricted binary Goldbach.

**Proof.** Infinitely many composite integers \(N=1\bmod L\) cannot
be covered by at most \(L\) positive prime summands: their number
would have to be exactly one. This proves the lower bound \(L+1\).
Fix one prime \(p_0\in{\cal P}_L\).
For each fixed \(j\ge3\), subtract \((j-3)p_0\) and use the ternary
prime theorem to represent every sufficiently large target
\(N=j\bmod L\) by \(j\) primes.
The \(L\) orders \(3,\ldots,L+2\) cover all residues, giving the upper
bound. With at most \(L+1\) summands the class \(2\bmod L\) admits
only two positive summands; every other residue is already covered
by orders \(3,\ldots,L+1\). This proves the equivalence. \(\square\)

Zero is padding, so this basis order means at most the stated number
of positive prime summands. Without that convention no fixed number
of primes \(1\bmod L\) covers every integer residue.

## 5. Why the semigroup theorem does not isolate primes

The \(C_{12}+C_{12}\) count has logarithmic scale \(N/(\log N)^{3/2}\);
two-prime counts have scale \(N/(\log N)^2\).
A semigroup leading asymptotic does not determine a subslice that
small. Its present logarithmic error also exceeds the predicted
prime subslice.

More structurally, the graph reconstruction requires order
\(\log\log N\) successful banks, contributing distinct prime factors
to both endpoints. A prime, or an endpoint with any fixed bounded
number of prime factors, belongs entirely to the discarded rare-bank
part. Prime indicators also fail the multiplication invariance used
by the transfer. Restricting that proof or differentiating a formal
prime-factor parameter does not establish a prime correlation.

The actual prime theorems above instead use prime-distribution inputs:
Green–Tao for the full-scale AP constellations, Halupczok and
Brun–Titchmarsh for the shrinking windows, Kane for the single principal-form prime equation, and Frei–Tinková
for the coupled nonabelian prime systems. The two-prime all-target assertion
and a pointwise prime-plus-semigroup theorem remain separate research
obligations.

A general prime singular series describes primes as local units.
A zero factor can still allow thin solutions with a variable pinned
to that same small prime; it excludes positive-density bulk tuples.
The explicit \(1\bmod L\) star systems above have no such additional
good-prime obstruction. Their target congruence is exact.

All moduli, fields, star sizes and window exponents are fixed.
The results are classical applications with the stated source scopes,
not a proof of binary Goldbach or a numerical classification of small
exceptions.
