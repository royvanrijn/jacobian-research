# The shared-summand C12 stretch: proved bounds and the remaining joint comparison

This note proves uniform upper bounds and unconditional natural-scale lower counts
for a positive proportion of admissible targets. The pointwise three-factor
asymptotic remains open.
Let
\[
 C=C_{12},\quad f=\mathbf1_C,\quad
 T(N)=\#\{a+b=N,\ b+c=N+12:\ a,b,c\in C\}
     =\sum_{1\le a<N}f(a)f(a+12)f(N-a).
\]
Write \(c=c_{12}\) for the Gamma-normalized marginal constant
\(\sum_{n\le X}f(n)\sim cX/(\log X)^{3/4}\), and
\(\mathcal E=\{p:p\nmid6,\ p\not\equiv1\pmod{12}\}\).

## 1. Exact finite admissibility and the model singular series

**Proposition SharedFiniteCriterion.** The system is admissible in every finite
quotient if and only if \(N\equiv2\pmod{12}\).

Every actual triple has \(N\equiv2\pmod{12}\).
This is also sufficient in every finite quotient. At a good forbidden
prime the free parameter must avoid
\[
             0,\quad -12,\quad N\pmod p.
\]
Every such prime is at least five, so at most three forbidden roots
leave a positive number of choices. They lift while preserving units.
At primes dividing twelve the common state \(a,b,c\equiv1\pmod{12}\)
is exactly compatible with \(N\equiv2\pmod{12}\), and this remains true
at every higher prime power.

For completeness, sufficiency needs only unit images, even at allowed
good primes. Given a modulus m, choose a modulo lcm(m,12) with
a=1 mod12, avoiding the three displayed roots at each prime p>=5
dividing m. Then a, N-a and a+12 are all units modulo m and are
1 modulo gcd(m,12). Every such unit image u is represented by a
single prime: the conditions p=u mod m and p=1 mod12 are compatible
by CRT and form a reduced progression, so Dirichlet supplies a prime
in C. Apply this independently to the three coordinate images.
Thus the equations hold modulo m with all three variables in C.
The actual values representing different moduli need not coincide.

For \(p\in\mathcal E\), \(p\nmid12\). Thus 0 and -12 are distinct,
and the number of forbidden roots is exactly
\[
 r_p(N)=3-\mathbf1_{p\mid N(N+12)}.
\]
The normalized local factor is
\[
 B_p(N)=\frac{1-r_p(N)/p}{(1-1/p)^3}.
\]
Put
\[
 C_3=\prod_{p\in\mathcal E}\frac{p^2(p-3)}{(p-1)^3}>0,\qquad
 J_3(N)=\prod_{\substack{p\mid N(N+12)\\p\in\mathcal E}}
                                      \frac{p-2}{p-3}.
\]
The baseline factors are \(1+O(p^{-2})\), so the product converges.
The exact arithmetic model for the desired leading term is
\[
 \boxed{\quad
 T(N)\ \stackrel{?}{\sim}\
       144c^3 C_3 J_3(N)\frac{N}{(\log N)^{9/4}},
       \qquad N\equiv2\pmod{12}.
       \quad}                                                  \tag{ST.1}
\]
The finite factor is \(12^3/12=144\): one free residue state modulo12
and three marginal progression factors. The real parameter length is N.
Unlike the diagonal \(c=a\), the forms \(a\) and \(a+12\) are distinct,
so the exponent is \(3(1-1/4)=9/4\), not3/2.
There are no exceptional good coefficient primes because every
good forbidden prime is coprime to12.

This does not assert (ST.1). The counterexample N26, with targets26,38,
is a genuine small global exception despite complete local admissibility.

## 2. A uniform pointwise upper bound

**Proposition SharedC12Upper.** Uniformly for \(N\ge3\),
\[
 T(N)\ll \mathbf1_{N\equiv2(12)}
                J_3(N)\frac{N}{(\log N)^{9/4}}.                \tag{ST.2}
\]

For an admissible N, sieve the parameter a=1 mod12 in [1,N).
At \(p\in\mathcal E\) remove the \(r_p(N)\) roots of the cubic
\(a(a+12)(N-a)\). Let
\[
 \rho_N(p)=r_p(N)/p,\quad
 g_N(p)=\rho_N(p)/(1-\rho_N(p)),\quad
 G_N(D)=\sum_{d\le D}\mu^2(d)g_N(d),
\]
with all divisors supported on \(\mathcal E\).

Here is a target-uniform lower bound for the Selberg denominator.
Take \(z=D^{1/100}\). The total Euler subset mass is
\[
 P_N(z)=\prod_{\substack{p\le z\\p\in\mathcal E}}
                    (1-\rho_N(p))^{-1}
 =\frac{\prod_{p\le z,p\in\mathcal E}(1-3/p)^{-1}}
        {\prod_{p\le z,p\in\mathcal E,p\mid N(N+12)}(p-2)/(p-3)}
 \gg\frac{(\log z)^{9/4}}{J_3(N)}.
\]
The prime-progression Mertens estimate gives the last exponent.
In the corresponding independent subset law,
\[
 \mathbb E\log d=\sum_{p\le z,p\in\mathcal E}
                         \rho_N(p)\log p\le6\log z+O(1).
\]
Markov therefore leaves a fixed positive proportion of the mass at
\(d\le D\). Consequently
\[
 G_N(D)\gg (\log D)^{9/4}/J_3(N).
\]

The optimal quadratic Selberg weights have \(\lambda_1=1\) and
\[
 |\lambda_d|\le \prod_{p\mid d}(1-\rho_N(p))^{-1}
                         \ll_\epsilon d^\epsilon
\]
uniformly in N. CRT gives \(r_N(m)=\prod_{p\mid m}r_p(N)\)
roots modulo squarefree m, and the progression count has error
\(O(r_N(m))=O_\epsilon(m^\epsilon)\).
Thus the total quadratic-weight error is \(O_\epsilon(D^{2+O(\epsilon)})\).
Taking \(D=N^{1/10}\) makes it negligible compared with
\(N/(\log N)^{9/4}\). The main term is \(O(N/G_N(D))\).
This proves(ST.2).

## 3. Uniform moments of the target correction

**Lemma SharedTargetMoments.** For each fixed integer k>=1,
\[
 \sum_{\substack{X<N\le2X\\N\equiv2(12)}}J_3(N)^k\ll_k X.       \tag{ST.3}
\]

Let
\[
 u_k(p)=((p-2)/(p-3))^k-1\ll_k1/p.
\]
Then \(J_3(N)^k=\sum_{d\mid N(N+12)}u_k(d)\), over squarefree
\(\mathcal E\)-divisors. All relevant primes are at most2X+12.
For any product d of those primes, the condition
\(d\mid N(N+12)\), together with N2mod12, has \(2^{\omega(d)}\)
roots modulo12d, since p does not divide12. Hence
\[
 \sum_NJ_3(N)^k
 \le \frac X{12}\prod_{\substack{p\le2X+12\\p\in\mathcal E}}
          \left(1+\frac{2u_k(p)}p\right)
       +O\left(\prod_{\substack{p\le2X+12\\p\in\mathcal E}}
          (1+2u_k(p))\right).
\]
The first product is bounded, since its logarithm has a convergent
\(O_k(p^{-2})\) sum. The second is a fixed power of log X by Mertens.
This proves(ST.3). The expansion over all subsets is harmless even
when d is huge: the entire endpoint sum is the displayed finite
Euler product. In particular no unjustified O(X) divisor count occurs.

For k=2, \(u_2(p)=(2p-5)/(p-3)^2\); the endpoint product is
\(O((\log X)^4)\). For k=1 the same argument gives the actual mean
\[
 \sum_{\substack{X<N\le2X\\N\equiv2(12)}}J_3(N)
 \sim\frac X{12}
      \prod_{p\in\mathcal E}\left(1+\frac2{p(p-3)}\right).
\]
The exact Euler consistency check is
\[
 C_3\prod_{p\in\mathcal E}\left(1+\frac2{p(p-3)}\right)
 =\prod_{p\in\mathcal E}\left(1-\frac1{(p-1)^2}\right)=C_2.
\]
This checks the proposed singular series; it does not establish(ST.1).

## 4. Unconditional positive-proportion natural-scale counts

**Theorem SharedC12PositiveProportion.** There are constants c1,c2>0
such that, for every sufficiently large X, at least c2 X integers
\(N\in(X,2X]\), \(N\equiv2\pmod{12}\), satisfy
\[
              T(N)\ge c_1 X/(\log X)^{9/4}.                   \tag{ST.4}
\]
This is a classical deduction, independent of the current ZF/G imports.

The needed primary theorem is Iwaniec--Pomykała, *Sums and differences
of quartic norms*, Mathematika40(1993), Theorem1, printedp.234 ([retained primary PDF](sources/iwaniec-pomykala-quartic-norms-1993.pdf)).
Specialize to Q(zeta12), difference12. Its unconditional count gives
\[
 \#\{1\le a\le X/4:a,a+12\in C\}\gg X/(\log X)^{3/2}.
\]
Their coprimality requirement removes nothing here:
\((a,a+12)\mid12\), while C values are coprime to6.
The shift is fixed, so no growing-target use of their conditional
sum theorem is involved.

The ordinary one-variable marginal gives
\[
 \#\{b\in C:X<b\le3X/2\}\gg X/(\log X)^{3/4}.
\]
Every choice of these a and b has N=a+b in(X,2X] and c=a+12 inC.
Thus
\[
 \sum_{\substack{X<N\le2X\\N\equiv2(12)}}T(N)
                        \gg X^2/(\log X)^{9/4}.               \tag{ST.5}
\]
On the other hand (ST.2)--(ST.3) give
\[
 \sum_N T(N)^2\ll X^3/(\log X)^{9/2}.
\]
Discarding values below a sufficiently small fixed multiple of
\(X/(\log X)^{9/4}\) removes at most half the first moment.
Cauchy--Schwarz on the remaining set then proves(ST.4).

This proves positive-proportion coverage with the expected scale.
It neither gives a leading asymptotic on those targets nor removes
the final exceptions.

## 5. A natural-scale three-factor envelope with pair records

The bank records needed by the fixed-shift graph are those of the pair
a,a+12. The third factor changes the required sieve scale, but does
not require recording a third valuation coordinate.

Write L=log N, y=exp(L^(2/5)), E=exp(sqrt L). Let P be the finite union
of the allowed prime banks, all p=1 mod12 and p<=y, after deleting bank
primes dividing N(N+12). For every a set
\[
                     r=a_P,\qquad s=(a+12)_P.
\]
These records are coprime. Their reference law is the ordinary
two-coordinate law
\[
 \mu_p(0,0)=1-2/p,\qquad
 \mu_p(e,0)=\mu_p(0,e)=(1-1/p)p^{-e}\quad(e\ge1).
\]
Write mu(r,s) for its product over P. This is also the marginal of the
three-coordinate valuation law. At a remaining bank prime the
third affine form has no norm restriction; summing its valuation
restores this two-coordinate law. In the usual growing banks,
V_g=sum_(p in P_g)1/p~1 and min P_g tends to infinity, so exactly one
distinct prime at each of a,a+12 has reference probability
V_g^2 exp(-2V_g)+o(1), tending to exp(-2)>0.

**Lemma SharedBankEnvelope.** The following quantitative envelope holds for
all eight actual/model products.
Assume the positive binary density model D has its globally bounded
smooth extension, so throughout the relevant interval
\[
 0\le D(n)\le C\lambda\,1_{n\equiv1(12)} f(n_S),
 \qquad \lambda\asymp L^{-9/20},
\]
where S is the set of primes at most y and n_S is the S-part.
For each W1,W2,W3 independently equal to f or D, and every event B
in the pair-record space,
\[
 \sum_{\substack{1\le a<N\\(r,s)\in B}}
          W1(a)W2(a+12)W3(N-a)
 \ll \frac{N J_3(N)}{L^{9/4}}\mu(B)
       +O\!\left(N^{1/8}+N e^{-cL^{1/10}}\right).                \tag{ST.6}
\]
Constants are uniform in target valuations and B. The model hypothesis
is the same positive pointwise bound used in the binary envelope.

To prove this, first retain r,s<=E. Exact positive valuations at
p^e||rs impose at most rad(rs)<=E² classes modulo
12*r*rad(r)*s*rad(s), with combined density
\[
 \frac{\beta(r,s)}{12},\qquad
 \beta(r,s)=\prod_{p^e\Vert rs}(1-1/p)p^{-e}.
\]
At every bank prime not dividing rs remove the roots0,-12.
These supply the remaining product in mu(r,s).

Let k_D of the factors be models, and k_F=3-k_D be true indicators.
Take D_s=N^(1/20) and z=D_s^(1/64), so z>y eventually. At forbidden
primes p<=y sieve all three forms; at forbidden primes y<p<=z sieve
only the k_F indicator forms. The bank-prime roots just described
form a disjoint prime set. Write rho(p) for the excluded-root fraction.
Uniformly rho(p)<=3/p and 1-rho(p)>=2/5.

For the Selberg denominator
G(D_s)=sum_(d<=D_s)mu²(d)prod_(p|d)rho(p)/(1-rho(p)), consider the
full independent Euler subset law through z, including bank primes.
Its expected log d is at most
\[
 3\log z+2\log y+O(\log\log z)
       =(3/64)\log D_s+o(\log D_s).
\]
Markov leaves a fixed positive fraction of this mass at d<=D_s.
Prime-progression Mertens and the exact collision factors give
\[
 G(D_s)\gg
 \frac{(\log y)^{9/4}
               (\log z/\log y)^{3k_F/4}}{J_3(N)}
 \prod_{\substack{p\in P\\p\nmid rs}}(1-2/p)^{-1}.              \tag{ST.7}
\]
At a low forbidden prime, a collision changes three roots to two,
with reciprocal-factor loss (p-3)/(p-2). For a subset of k_F forms
at a high prime, any collision loss is no larger. Thus J3 bounds
all these losses uniformly, including arbitrary target valuations.

The optimal weights satisfy |lambda_d|<<d^(.01), and a squarefree
sieving modulus has at most3^omega(d)<<d^(.01) root classes. The
quadratic CRT remainder for one record is therefore
O(E² D_s^(2.04)). Summing at most E² retained records gives
\[
 E^4D_s^{2.04}=N^{.102+o(1)}\le N^{1/8}
\]
eventually. The main term from (ST.7), including the model amplitudes,
has precisely the desired logarithmic scale, since
\[
 \lambda^{k_D}(\log y)^{-9/4}
           (\log z/\log y)^{-3k_F/4}
 \asymp
 L^{-9k_D/20-9/10-9k_F/20}=L^{-9/4}.
\]
The exact-valuation density times the inverse bank product is
mu(r,s)/12. Summation over the retained records in B proves the
main term in (ST.6).

Finally put theta=1/log y. Positivity of the divisor-convolution
identity for n_P^theta gives
\[
 \sum_{n\le2N}n_P^\theta
 \le2N\prod_{p\in P}
                 \frac{1-1/p}{1-p^{\theta-1}}\ll N.
\]
The logarithm of this product is
O(theta sum_(p<=y)log p/p)=O(1). Rankin therefore bounds each
omitted tail r>E or s>E by
O(N exp(-theta log E))=O(N exp(-L^(1/10))).
Each actual/model factor is at most one for large N. This proves(ST.6).

For a nonnegative pair-record test W<=L^C, integration of (ST.6)
gives the same main bound with mu(B) replaced by E_mu W and with
the displayed error multiplied by L^C. The resulting power and
exponential savings remain negligible at the natural scale.
Thus rare-bank estimates can be transferred at the true
three-factor scale without dropping the third factor.
This envelope does not supply the missing joint comparison.

## 6. Tested graph lift and the actual obstruction

The fixed-shift12 graph initially looks preferable to the reflected
N-shift graph. In [Walsh v2](sources/walsh-2610.40439v2.pdf),
printed pp.38–41, write
\[
 n=D_Rw,\quad n'=D_R(w+12D_C),\quad D=D_RD_C.
\]
The raw pattern has D_C=1 and represents w=a,n'/D_R=a+12.
On complete target lists let R(n') be the product of all distinct
bank divisors and let Q0 be the product of the designated raw
unshared target slots. Before permutation,
\(R(n')=D_RQ0\), and the tentative marked endpoint
\[
 g_N(n',q)=\beta(n')f((N+12)R(n')-n'Q0)
\]
has exactly the desired raw third factor f(N-a), by allowed-prime
dilation invariance. It is bounded, and the source exact pairing
allows its second endpoint to depend on marks.

However the actual operator is \(\Pi(\sum_\nu A_\nu)\Pi\), with
independent full list-permutation projections. The right projection
averages Q0 over all target labels. It no longer tracks the actual
edge's unique unshared target prime. The desired Q0 is determined
by the shared source-target intersection, not by the target integer
or target vertex alone after that projection.
Therefore this tentative endpoint does not produce the desired raw
weighted correlation. Omitting Pi is not an allowed source specialization.

There is an exact finite negative control on a genuine raw edge. Take
one bank, M=2, shared prime13, unshared source61 and target73, and
\[
 n=13\cdot61=793,\quad n'=13\cdot73=949,\quad N=194.
\]
Then \(n'-n=12\cdot13\), and all three bank primes are \(1\bmod12\)
and avoid \(12N(N+12)\). Here R(n')=949. The designated unshared
target label gives
\[
 f(206\cdot949-949\cdot73)=f(949\cdot133)=0,
\]
since 133=7*19. Swapping the two target slots gives instead
\[
 f(206\cdot949-949\cdot13)=f(949\cdot193)=1,
\]
since193 is a prime1mod12. Thus Pi changes the proposed raw third
weight from zero to1/2. This is a finite algebraic control on the proposed endpoint identity.
It does not by itself construct an admissible bank family in the full
asymptotic parameter regime, and it is not a counterexample to the
desired eventual three-factor theorem.


One can try coloring each bank T disjoint union S and restricting a
complete physical vertex to one T prime and M-1 S primes. The unique
T prime then is numerical and permutation invariant. Removing edges
whose T primes coincide by diagonal pinching forces T not to be
shared. Pinching itself can be represented by averaging numerical
diagonal-unitary conjugations, so the physical pairing estimate is
stable with a polylogarithmic factor.

The difficult issue is normalization and parameter order. All shared
labels are then in S. In a pattern with t_g unshared slots the
shared-label sum contributes
\[
                    (V_{S,g}/V_g)^{M-t_g}
\]
instead of one. This loss must be retained when reading off the
arithmetic relation. Improved colored endpoint norm bounds can pay
part of it; they do not justify simply deleting it.

In the direct strategy that selects active banks from colored successes,
choosing T mass about V_g/M makes the common shared loss bounded per
bank, but makes the initial successful-bank fraction of order M^-2.
That fraction enters the ideal/probe construction before M.
The [weighted-dilation graph source](sources/weighted-dilation-graphs-sep24-2026.pdf),
§4.4, printed pp.32–33, equations(4.22)–(4.24), has for ell=1 the conditions
\[
 M\ge m(m+1),\qquad
 -\frac{c_0}{6}\log(2/m)>E_{\rm id}+3,
\]
where c0 is the lower active-bank density on the loglog scale.
With c0 of order M^-2, these require logm>>M² and logm<=0.5logM
simultaneously. Thus that particular color choice cannot be closed
merely by choosing a small final error exponent after M.
This is a verified obstruction in the available quantitative proof,
not a theorem ruling out all colored variants or a sharper new estimate.

Retaining all ordinary active banks and averaging colors avoids that
particular active-bank selection, but costs the inverse coloring
probability. With about log L banks its M-dependent per-bank loss is
a power of L, rather than a constant. A reconstruction paying that
loss within the inherited error budget has not been proved.
Neither calculation rules out a better weighted reconstruction.

A second failed shortcut would replace the fixed-shift pair by its
model against an arbitrary bounded third weight. Such an assertion is
false: choosing the weight f(a)f(a+12) makes the actual weighted pair
count full size, while the model pair on that support is smaller by
its pointwise density bound. The arithmetic third weight is special;
its structure must enter a valid new comparison theorem.

## 7. A valid one-sided vector Fourier lemma for a future bridge

This lemma removes a spurious sqrt(Y) loss if a future graph lift
produces the indicated one-sided dependence. It does not repair
the permutation obstruction by itself.

**Lemma SharedVectorFourier.** Let nu1,nu2 be probability weights on dyads of lengths U,V, with
nu1(u)<=L^K/U and nu2(v)<=L^K/V. Let |a_u|<=1.
For arbitrary complex B_v(theta), the elementary integer bilinear
bound gives
\[
 \left|\sum_{u,v}\nu_1(u)\nu_2(v)a_uB_v(\theta)e(\theta uv)\right|^2
 \ll L^{2K}\mathcal B(\theta)
                 \sum_v\nu_2(v)|B_v(\theta)|^2,
\]
whenever s>=1, gcd(r,s)=1, |theta-r/s|<=s^-2, and
\[
 \mathcal B(\theta)=1/s+1/U+\log(2s)/V
                                 +s\log(2s)/(UV).
\]
Indeed apply the geometric-sum operator estimate to
A_u=nu1(u)a_u and C_v=nu2(v)B_v, then use
||A||_2²<=L^K/U and
||C||_2²<=(L^K/V)sum nu2|B_v|².

If \(B_v=\widehat{g_v}\), \(g_v(x)=H_v(x)g(x)\), and |H_v|<=1,
Parseval gives
\[
 \int\sum_v\nu_2(v)|B_v(\theta)|^2\,d\theta\le\|g\|_2^2.
\]
Thus against a centered \(\hat\xi\), minor arcs with
mathcalB<=L^(-2A-2K) cost at most
\(L^{-A}\|\xi\|_2\|g\|_2\), without a sqrt(Y) loss.
On major arcs, the trivial weighted Cauchy bound and a suitable
centered-energy estimate give the same type of control.

The hypothesis is genuinely one-sided: H_v may not depend on u or
on the centered endpoint's selected small prime. The current
three-factor graph lift has not been shown to produce that
representation after all permutation and shared-label operations.
Accordingly this is an available analytic tool for a future joint
bridge, not a proof of(ST.1).

## 8. Precise completion target

A successful pointwise argument must construct a source-compatible
joint graph comparison or a distinguished-slot replacement with:

1. correct raw three-factor weight after both permutation projections;
2. recovery/normalization and endpoint exponents controlled before
   the eventual mark choice;
3. either one-sided vector comparison coefficients as in Section7,
   or a new joint Fourier energy estimate handling their true dependence;
4. a bound \(o(NJ_3(N)/(\log N)^{9/4})\), after all coefficient costs.

The positive envelope, exact finite criterion, singular series,
uniform upper bound and unconditional positive-proportion lower
counts are already supplied above. The all-target lower bound and
pointwise leading asymptotic remain unproved.
