# Proof support for the conditional pointwise binary theorem

This file contains the expanded arguments for
[ConditionalBinarySemigroups](PRIMITIVE_NORM_BINARY_CORRELATION.md).
Its full conclusion is conditional on the explicitly identified ZF and
graph inputs. Individual modules state which inputs they use.
The density, finite-AP and arithmetic arguments do not independently
prove the external graph moment or the zero-free theorem.

The source ledger records exact primary editions and internal review
limits. This is written proof support, without an unrestricted
computational or formal-verification claim.

## D. Density models, principal means, bank envelopes and reconstruction

19 November2026. Expanded proof module. This does NOT
claim the missing all-target graph cancellation theorem has been proved.

### Exact source and present theorem boundary

The actual primary inspected is Miguel N.Walsh, Shifted convolutions of
divisor functions, arXiv:2610.40439v1,29October2026,44pages,
https://arxiv.org/pdf/2610.40439v1. Retained primary bytes:
sources/walsh-2610.40439v1.pdf, SHA256
21efec3645622c74acde7587bff811719b6c19801b07b253559b64f85ed259a7.
Sections2--4 and their exact hypotheses were inspected. Theorem1.1 concerns
positive INTEGER divisor orders and a FIXED nonzero shift. Section1.2.4
mentions sums-of-two-squares indicator correlations and sublinear growing
shifts as possible extensions, not as proved theorems. Neither primitive
prime-support indicators nor the reflected growing-target line are a
literal theorem specialization.

This module gives the necessary indicator density model, its principal
means, the model*model main term, a target-uniform bank-record envelope,
and the exact reconstruction reduction. The core cases include the odd
primitive Gaussian core, and the same argument works for every fixed
abelian prime-support subgroup of positive density. The banked comparison is assembled in the main conditional proof from Modules I,G,P,E; it is not a literal specialization of the integer-divisor theorem.

Walsh's displayed high-frequency prime bound, Lemma5.1, uses his reference
[39], OpenAI's claimed zero-free half-plane Re(s)>7/8 manuscript. The
particular P^(-1/16) prime-sum error uniform to heights x^C must not be
attributed to a classical zero-free region without reconstruction. The
low-frequency argument below needs only polylogarithmic moduli and can
use the classical zero-free region plus Siegel's theorem. The new graph
and high-frequency imports require their own source/hypothesis audit.

### 1 Fixed core, notation and exact branches

Fix even L0, a subgroup H of (Z/L0)^* of index kappa>=2, density
rho=1/kappa, and let B0 be the primes dividingL0. Let C contain1 and
integers supported on primes p not inB0 with p modL0 inH. Write f=1_C.
Let U,V be fixed subsets ofH; these can carry all fixed genus states.
All constants may depend on this fixed data. Put

  L=log X, y=exp(L^(2/5)), S={primes<=y},
  V_S=product_{p<=y}(1-1/p), E=exp(sqrt L),
  E_+(y)=product_{p<=y,p allowed}(1-1/p),
  E_-(y)=product_{p<=y,p notinB0,p forbidden}(1-1/p),
  alpha0=phi(L0)/L0.

For largeX, V_S=alpha0 E_+(y)E_-(y) exactly. Dirichlet characters
eta in Hperp give

  F(s)=sum f(n)n^-s=product_eta L(s,eta)^rho * G(s),

where the correction is absolutely convergent for Re(s)>1/2. At a prime
whose order moduloH is r>=2, its factor is(1-p^(-rs))^(1/r).
Let c_C be the actual one-variable constant, so

  F(s)=(s-1)^(-rho) H0(s), H0(1)=Gamma(rho)c_C>0.

The rough core series is

  F_R(s)=sum_{(n,S)=1}f(n)n^-s
        =F(s)product_{p<=y,p allowed}(1-p^-s)
        =(s-1)^(-rho) H_R(s).                           (W1)

Thus H_R(1)=Gamma(rho)c_C E_+(y). In a disk of radius c/log y,
fixed-order Cauchy/derivative estimates give

  |H_R^(j)(1)|/j! <<_j H_R(1)(C log y)^j.                (W2)

Every twist eta inHperp is IDENTICAL to f on the support. All other
fixed twists have no principal branch at1. This is why a scalar model
without a character projection would be wrong, already for GaussianC4.

### 2 A smooth rough density of arbitrary logarithmic precision

For any sufficiently large fixed integer J define, on
X^a<=t<=X^C with fixed0<a<C,

  M_R,J(t)=sum_{j=0}^J H_R^(j)(1)/j!
             * (log t)^(rho-j-1)/Gamma(rho-j).            (W3)

The leading term is c_C E_+(y)(log t)^(rho-1); every following term
has relative size O_j((log y/log t)^j). Hence

  M_R,J(t)=c_C E_+(y)(log t)^(rho-1)(1+O_J(L^(-3/5)))>0. (W4)

Extend it smoothly and nonnegatively below the bulk range, keeping
M_R,J(t)<< E_+(y)L^(rho-1) on0<t<=X^C. Such an extension does not
change any bulk statement below. IncreasingJ is allowed to obtain any
prescribed fixed logarithmic accuracy; J does not grow withX.

Here are the details of the standard one-variable estimate needed in
place of Walsh's integer-order pole computation. If q<=L^B, package its
fixed characters together withL0 in Q=lcm(q,L0). Let W be a smooth dyadic
weight at T, X^a<=T<=X^C, with fixed derivatives bounded by fixed powers
ofL. Character orthogonality reduces progression sums to

  F_R(s,psi)=product_eta L(s,psi*eta)^rho G_psi(s)
                 product_{p<=y,p allowed}(1-psi(p)p^-s). (W5)

The correction G_psi is uniformly absolutely convergent to the right
of1/2. All prime factors ofQ are<=y. If psi is an inflated member
ofHperp, F_R(s,psi)=F_R(s) pointwise, even thoughQ grows. Otherwise
none of the character factors is principal.

Use smooth Mellin inversion with height T0=exp(L^(1/10)), and shift to
Re(s)=1-c/log y. The classical zero-free region gives a substantially
wider zero-free strip at those heights and polylogarithmic conductors.
For a possible exceptional real character, Siegel's theorem with a
sufficiently small fixed exponent, depending onB, puts its zero farther
than c/log y from1 onceX is large. This is ineffective in the numerical
threshold, which is permitted here. Thus the nonprincipal branches
can be defined analytically on the shifted rectangle. The small Euler
product in(W5) is at most L^O(1) there: p^(c/log y)<=exp(c). Standard
polynomial bounds for the L-functions cost exp(O(L^(1/10))), whereas
T^(s-1) contributes exp(-c L^(3/5)). The smooth Mellin tails can be
made exp(-c L^(1/10)) by taking a sufficiently large fixed derivative
order. This proves an error O_A(T L^-A) for nonprincipal twists.

For the principal branches, Hankel integration arounds=1 and Taylor
expansion(W2) give exactly the smooth density(W3). The relative Taylor
remainder is

  O_J((C log y/log T)^(J+1))=O_J(L^(-3(J+1)/5)),

and the outer Hankel endpoint is exponentially small. For a prescribed
A, chooseJ large enough to absorb all fixed derivative/character-count
powers ofL. There is no need to shift a fractional L-power through zeros
at moduli exp(sqrt L). ONLY q<=L^B is being claimed.

The rough UNWEIGHTED population, with series
zeta(s)product_{p<=y}(1-p^-s), has principal densityV_S by the same
argument, and nonprincipal mean O_A(T L^-A). Consequently, for every
bounded periodic b(t) of periodq<=L^B, the two quantities

  sum_t W(t) b(t) f(t)1_{(t,S)=1},
  sum_t W(t) b(t) 1_{(t,S)=1}
        (sum_{eta inHperp}eta(t)) M_R,J(t)/V_S

have the SAME principal mean and differ by O_A(T L^-A), after choosingJ.
For the first, allHperp twists of the rough core coincide; for the second,
only the principal character of the unweighted rough population remains.
The finite-character projection makes these expressions agree exactly
at the principal terms. Fixed unit-residue filters such asU are inserted
into b; no character branch is omitted.

This proves the centered principal-mean interface needed at low
frequencies, including smooth factors t^(iu), |u|<=L^B, by enlarging
J and the derivative budget. It uses classical one-variable inputs,
not an additive correlation theorem.

### 3 The corrected density model and exact dilation invariance

For U subsetH put

  D_U(n)=kappa*1_{n modL0 inU} f(n_S) M_R,J(n_(S'))/V_S,
  xi_U(n)=f(n)1_{n modL0 inU}-D_U(n).                     (W6)

Then D_U is nonnegative and

  D_U(n)<< L^(-(3/5)(1-rho)), |xi_U(n)|<<1               (W7)

on all relevant scales. The same bound is imposed by the smooth
extension outside the bulk. Conditioning on n_S=A<=E preserves the
centered mean estimate of Section2: the residual state isU/A, and
arbitrary small-period weights may depend onA. The sum over suchA
costs only polylogarithmic harmonic factors. For example
sum_{A supported onS}d_K(A)/A is a fixed power ofL, so the bound also
allows Walsh's divisor-bounded record coefficients. Records A>E are
negligible by the elementary tilted moment and Rankin bound

  #{n<=X^C:n_S>E} << X^C L^O(1) exp(-c L^(1/10)).         (W8)

For exact invariance WITH FIXED STATES, select banks only from primes
p=1 modL0. Their density1/phi(L0)>0 is fixed. There are still
c_L loglogX disjoint small banks and big banks in Walsh's ranges,
with reciprocal masses in[1/2,2], by taking wider fixed loglog windows.
Every suchp is allowed and is<=y, so

  f_U(pn)=f_U(n), D_U(pn)=D_U(n), xi_U(pn)=xi_U(n)        (W9)

EXACTLY, including higher prime powers. Multiplication leaves the rough
part unchanged and preserves the fixed residue state. If U=H one may
instead use all allowed primes; identity-class banks avoid having to
assert unchanged genus states under nonidentity primes.

For a targetN comparable toX, delete bank primes dividingN. Uniformly
inN this changes each reciprocal mass by at most

  O(log N/(P_min log P_min)), P_min>=exp(L^(1/8)),

which is smaller than every negative power ofL. Thus the bank count and
mass bounds survive. This deletion will make all endpoint record pairs
coprime. It does not alter the core or its model outside the selected
bank set.

### 4 Unconditional evaluation of model*model on the positive line

Fix coprime B0-supported coefficientsa,b>0 and fixed subsetsU,V ofH.
Let Gamma_{a,b;U,V}(N) beL0/|H|² times the number of unit pairs
(u,v) inU*V with a*u+b*v=N modL0. Put

  J_H(N)=product_{q|N,q notinB0,q forbidden}(q-1)/(q-2),
  C_H=product_{q notinB0,q forbidden}(1-1/(q-1)²).

For Gaussian/Hquadratic this is theJ,C_in of the leading-constant note.
Use a fixed smooth bulk weightPhi(a*x/N) in epsilon<a*x/N<1-epsilon.
After(W8), log(n_(S'))=logN+O(sqrt L) uniformly in the bulk, so(W4) gives

  M_R,J(n_(S'))/V_S
       =lambda_X(1+O_J(L^-1/2)),
  lambda_X=c_C L^(rho-1)/(alpha0 E_-(y)).                 (W10)

Thus it suffices to count the small-prime sieve set with n,m in the
specified fixed residues and no forbidden prime<=y, with their density
model weights replaced bykappa² lambda_X². The fundamental lemma of
the fixed-dimensional combinatorial sieve is asymptotically exact here:
choose levelD=N^epsilon0, while
u=logD/logy asymp L^(3/5) tends to infinity. CRT gives an errorO(nu(d))
at every squarefree sieving modulusd, uniformly in the target. The
dimension upper condition is fixed2(1-rho), since target-prime mergers
only decrease the root densities. Standard fundamental-lemma upper and
lower weights give relative error exp(-c u logu), and their endpoint
sum is N^epsilon0 L^O(1), negligible. All fixed residue classes may be
summed first; their number is fixed.

The raw bad-place probability on the integer-line parameter is

  #fixed roots/L0=Gamma*alpha0²/kappa².

Indeed x=x0+b*k,v=v0-a*k parameterizes the integer line. ModuloL0 this
map bijects k with all roots of a*x+b*v=N: at every prime at least one
of a,b is a unit, and the local kernel has precisely the image(b,-a).
Thus bad coefficients need not be invertible simultaneously moduloL0.

The real line parameter has lengthN/(ab), and every good forbidden
prime has2 roots, merging to1 exactly whenq|N. Consequently the model
correlation has main coefficient

  (N/(ab))*Gamma*alpha0²/kappa² * kappa² lambda_X²
           * product_{q<=y,q forbidden}(1-(2-1_{q|N})/q)
 = (N/(ab))*Gamma*c_C² L^(2rho-2)
           * product_{q<=y,q forbidden}
                (1-(2-1_{q|N})/q)/(1-1/q)².             (W11)

The last product is C_H(y)J_H,y(N). Its omitted target-prime correction
is uniformly1+O(logN/(y logy)); the nontarget Euler tail isO(1/y).
Both errors are superlogarithmically small. This proves

  sum_{a*x+b*v=N} Phi(a*x/N)D_U(x)D_V(v)
   = (c_C² Gamma C_H+o(1))*J_H(N)*N L^(2rho-2)/(ab)
                         * integral_0^1 Phi(t)dt.        (W12)

All constants and the o(1) are uniform in N's prime valuations. A zero
fixed density gives zero exactly. Endpoint strips may be removed by
the same upper sieve, with lengthO(epsilon N), and thenepsilon->0.

This bypasses Walsh's Proposition3.1 model*actual calculation, whose
fractional-order version at moduli exp(sqrt L) would require a separate
argument around L-function zeros. No such large-modulus fractional
Selberg--Delange assertion is needed for(W12).

### 5 A direct target-uniform bank-record enveloping measure

Let C_banks be the union of the selected identity-class banks, after
deleting primes dividingN. For q in C_banks define a probability on
valuation pairs by

  mu_q(0,0)=1-2/q,
  mu_q(e,0)=mu_q(0,e)=(1-1/q)q^-e, e>=1,
  mu_q(e,f)=0 ife,f>=1.                                 (W13)

It has total mass1. Let mu_C be its independent product over the bank
primes. For any nonnegative record functionF(r,s)<=L^B, and any pair
of weightsw1,w2 from{f_U,D_U,f_V,D_V}, one has the envelope

  sum_{a*x+b*v=N} F(x_C,v_C) w1(x)w2(v)
   << N J_H(N)L^(-2(1-rho)) E_mu F
       +N L^O(1) exp(-cL^(1/10))+O(N^(1/2)).             (W14)

The fixed coefficients and fixed supports affect the constant only.
The statement is for a fixed bulk dyadic interval, or a fixed compact
part of the positive line. All target valuations are uniform.

Proof. First remove r,s>E by(W8), since every original/model weight
is bounded. There are at most E² possible retained pairs. If r,s share
a bank prime the stratum is empty, since no bank prime dividesN. For
coprime r,s the divisibilities r|x,s|v define a single additive-line
progression, with lengthN/(abrs). At a bank prime outside r*s,
exactness of the record excludes2 roots of the residual forms. At a
bank prime dividing r*s it excludes only the remaining variable root,
giving1 root. At every forbidden prime both residual slopes remain
units, because the record primes are allowed. Its roots merge iffq|N.
There is no coefficient-radical loss.

Apply the direct Selberg sieve using primes up tow=z^(1/16), z=N^(1/10).
For largeN this includes every bank prime, sincew greatly exceedsy.
The same Markov truncation as MixedPrimeSquareSieve applies, now with
only prime moduli. The Euler product from the bank conditions is
exactly

  product_{q inC_banks,q notdivrs}(1-2/q)
  *product_{q|rs}(1-1/q).

Together with1/(rs) it is mu_C(r,s). The fixed bad conditions can be
dropped for the upper bound.

For f*f the NONBANK SIEVE DENOMINATOR Euler product has size
(logw)^(2(1-rho))/J_H(N). For D*f, only the second variable excludes
forbidden primes above y; its size is bounded below by a constant
multiple of (logw logy)^(1-rho)/J_H(N). Multiply by the bound
D<<L^(-(3/5)(1-rho)) to recover the SAME L^(-2(1-rho)) factor.
For D*D only primes<=y are forbidden; their sieve denominator has size
(logy)^(2(1-rho))/J_H(N), and the two model bounds again give the
same factor. Target corrections can be bounded by the fullJ_H(N)
in every case. These comparisons are uniform because target mergers
only remove one root and their exact enhancement is(q-1)/(q-2).

The CRT endpoint is O(z^(4+epsilon)) uniformly in r,s, even using the
looser mixed-square budget. Sum it over r,s<=E andF<=L^B; this is
N^(0.45)exp(O(sqrt L))L^B=O(N^(1/2)). The retained main terms sum to
the expectation in(W14). This proves the envelope for all four pair
types and hence for |xi_U| times either actual or model weight.

Unlike Walsh Corollary4.2, this proof does not assume a fixed resultant,
positive slopes on both residual forms, or a fixed shift. It is direct
CRT counting on the reflected positive-sum line.

### 6 Rare-bank control and exact reconstruction

For a bankg putI_g=1_{omega_g(x)=omega_g(v)=1}, counting distinct primes.
Under(W13), the events for disjoint banks are independent and

  E_mu I_g=(product_{p inP_g}(1-2/p))
                 *sum_{p!=q inP_g}1/((p-2)(q-2))
           =V_g² exp(-2V_g)+O(sum_{p inP_g}p^-2).

Since1/2<=V_g<=2, this is at least a fixedbeta>0, for example
beta=e^-4/8 for largeX. There are at leastc_L loglogX banks of each
size. Chernoff therefore makes the probability of fewer thanhalf
beta times this many successful banks O(L^-eta_L) for someeta_L>0.
Use(W14) with that bad-bank indicator to obtain

  sum |xi_U(x)| w2(v) (1-G(x,v))
       << N J_H(N)L^(-2(1-rho)-eta_L)+o(N L^-A)           (W15)

for a suitably defined good-bank indicatorG and fixedA. Records with
prime squares from active banks are negligible by the bounded weights
and sum_p p^-2<<exp(-cL^(1/8)). No divisor fourth moments are needed.

For the remaining good-bank set, exactly the Walsh inclusion-exclusion
identity applies: G=sum_T c_T product_{g inT}I_g, where every activeT
has a fixed positive proportion of small AND big banks, and
sum_T|c_T|<=3^(#banks)=L^O_L(1). Put

  Omega_P(n)=product_{g active}1_{omega_g(n)=1}/V_g.

Since all original/model/centered factors are invariant under these
bank primes, there is NO divisor-weight factor kl. The exact restricted
correlation is simply

  xi_U(x)w2(v) product_{g active}I_g
    =(product_{g active}V_g²)
             xi_U(x)w2(v)Omega_P(x)Omega_P(v).           (W16)

The product and total reconstruction coefficient cost at most another
fixed power ofL. Thus a banked estimate of arbitrary logarithmic
strength, uniform for the growing reflected target and forw2=f_V OR
w2=D_V, implies

  sum_{a*x+b*v=N} Phi(a*x/N)xi_U(x)w2(v)
       =o(N J_H(N)L^(-2(1-rho))).                       (W17)

Together with its symmetric counterpart, the exact identity
f_U f_V-D_U D_V=(f_U-D_U)f_V+D_U(f_V-D_V), and(W12), this proves the
actual unweighted binary-core leading term. This is a REDUCTION,
not a claim that the banked graph estimate has already been supplied.

### 7 Interface to the full conditional proof

Walsh Proposition4.1 is printed for integer divisor weights and fixedh.
To turn this module into a complete pointwise theorem one must prove its
banked comparison for the bounded invariant functions above on

  a*x+b*v=N, x,v>0, N comparable toX,

with the reflected slope, target-dependent prime deletion, fixed
identity-class prime banks and all fixed residue states. The low-frequency
principal-mean issue has the explicit rough model of Sections1--3.
Modules I,G,P,E supply the graph and target-grid transfers under the explicit deep inputs. None of those conclusions is inferred merely from the model/envelope calculation.

The generic rho interface includes H={1} modulo12, rho=1/4, and therefore
identifies a genuine route beyond the existing C12+C12 sieve barrier.
The C12 conclusion additionally uses the banked comparison proved conditionally in the main note.
No historical novelty, external review or completed binary asymptotic is
claimed by this module.


### 8 Mixed fixed subgroups: the same module

For two fixed subgroups H1,H2 with densitiesrho1,rho2, use their separate
rough models and projectionskappa1*1_U,kappa2*1_V. The common identity
prime banks preserve both factors. The ordinary joint record law(W13)
and the reconstruction are unchanged. The sieve dimension is

  (1-rho1)+(1-rho2)=2-rho1-rho2.

At a good prime excluded by only one subgroup, the relative pair factor
is(1-1/p)/(1-1/p)=1, whether or not p dividesN. Only primes excluded by
BOTH sides contribute the convergent singular product and the target
factor(p-1)/(p-2). The finite density is

  Gamma=L0/(|H1||H2|)*#fixed unit roots.

Its raw probability isGamma*alpha0²/(kappa1*kappa2), exactly cancelling
the two projection factors. Thus the model main coefficient is

  c_H1*c_H2*Gamma*C_excluded_both*J_excluded_both(N)
            *N(logN)^(rho1+rho2-2)/(a*b).

For example the D_H1*f_H2 envelope has denominator
(logy)^(1-rho1)(logw)^(1-rho2); its first-model scalar supplies the
remaining L^(-(3/5)(1-rho1)), giving the same total dimension.
If a subgroup is the whole unit group, its model is exact and its
centered remainder vanishes; handle that trivial case separately.
No new mixed-group obstruction appears in the modules proved here.
The all-target banked comparison remains the same explicit missing
analytic obligation. This statement does not assert any joint norm
class-group intersection theorem.

## I. The fixed-AP ideal kernels

### Statement and imported source

This note proves the fixed-arithmetic-progression variant of Theorem 4.1 in
OpenAI, *Weighted dilation graphs, smooth shifted primes and totient fibers*
(24 September 2026), pp. 27–34, using its actual construction, not an assertion
that interval prime sets can be replaced without modification.

Inspected source:
- sources/weighted-dilation-graphs-sep24-2026.pdf
- PDF bytes 802390; SHA-256
  351692885d112277e969901e75fe56eeee0d263d5d5847d0b980fe2c8c7744c0.
- Essential source statements: Definition 3.1, Theorem 2.2,
  equations (4.3)–(4.26), Lemmas 4.2–4.3 and Theorem 4.1.
  Section 3 admits arbitrary pairwise-disjoint finite prime sets satisfying
  its location and harmonic-mass bounds. The extra interval assumption
  appears at the start of Section 4 and is used in the proof of (4.9).

Fix an integer L0. Replace each big prime group by
\[
 {\cal P}_g=\{p\in I_g:p\equiv1\pmod {L_0}\}\setminus{\cal E}_g,
 \qquad I_g\subset[\exp(L^c),\exp(L^d)],
 \qquad L=\log x,\quad T=\log L.
 \tag{AP1}
\]
The intervals are pairwise disjoint, the numbers of groups satisfy the
source's bounds, and
\[
 v_-\le V_g=\sum_{p\in{\cal P}_g}1/p\le v_+,\qquad
 \mu_g(p)=1/(V_gp).
\]
The constants c,d,L0,v_-,v_+ and the group-count constants are fixed.
The deletion sets may be empty. More generally allow deletions satisfying
\[
 \sup_g\sum_{p\in{\cal E}_g,\ p\equiv1(L_0)}1/p
       =O_A(L^{-A})\quad\hbox{for every fixed }A.
 \tag{AP2}
\]
All assertions are uniform over such sets when the constants in (AP2)
are uniform.

**Fixed-AP ideal-kernel theorem.** With these laws, Theorem 4.1 remains
valid, including its norm L^{-E_id}, uniformity in every real additive
frequency and in |zeta|<=L^{C4}, polylogarithmic total coefficient cost,
and the requirement that all defining parameters be chosen independently
of the eventual fixed mark length J. Its kernels require the quotient
character construction below in place of the unmodified source (4.7).

This is a consequence of the inspected proof and its Siegel–Walfisz
input. It does not independently verify the full graph-transference
theorem, any claimed quasi-Riemann hypothesis, or the final binary-core
correlation. Those are separate dependencies.

### 1. Why the original character list must change

Let C_Q consist of primitive Dirichlet characters of conductor at most Q.
When all labels are 1 modulo L0, two different members of C_Q can define
the same function on every label product. For example, for L0=4 the
conductor-1 character and chi_4 both equal 1 there. The unmodified sum
in source (4.7) reproduces the constant function twice rather than once.
For larger L0 its constant-function Gram eigenvalue can exceed two.
Thus the source proof of its (4.9), and its reproduction identity (4.18),
cannot be copied literally.

For chi,chi' in C_Q define
\[
 \chi\sim\chi'
 \quad\Longleftrightarrow\quad
 \operatorname{cond}(\chi\overline{\chi'})\mid L_0,
 \tag{AP3}
\]
where the conductor means the conductor of the primitive character
inducing the quotient on a common multiple of the two conductors.
Equivalently, the quotient is induced from a character modulo L0.
Choose one representative of each equivalence class in C_Q, and call
the resulting set R_Q. In particular #R_Q<=#C_Q<=Q^2, and every
representative still has conductor at most Q.

Here is the exact group-theoretic justification. On a common modulus
M=lcm(L0,cond chi,cond chi'), reduction from (Z/M)* to (Z/L0)* is onto.
Its kernel H_M consists of units a congruent to 1 modulo L0. A character
on (Z/M)* is trivial on H_M precisely when it factors through (Z/L0)*,
which is equivalent to its primitive conductor dividing L0. Therefore
(AP3) is an equivalence relation, and equivalent characters agree
EXACTLY on each bank prime and on every product of bank primes.
For inequivalent representatives,
\[
 \sum_{a\in H_M}(\chi\overline{\chi'})(a)=0.
 \tag{AP4}
\]
No unproved assertion about the distribution of the individual product
integers is used here.

An equivalent implementation keeps all members of C_Q but weights each
one by the reciprocal of its equivalence-class cardinality. The
representative construction is simpler and will be used throughout.

### 2. Conditional character orthogonality on AP cells

For a selected nonempty slot set Lambda using at most m slots per group,
let D be the product of its independently sampled labels, and put
I_j=[jh,(j+1)h), nu_j=P(log D in I_j).
For inequivalent representatives chi,chi' and any fixed A,
\[
 \left|E\,1_{\log D\in I_j}
             \chi(D)\overline{\chi'(D)}\right|\ll_A L^{-A},
 \tag{AP5}
\]
uniformly in j, Lambda, and the permitted banks.

To prove this, hold all but one label fixed. The cell and I_g restrict
the remaining prime to a further interval, possibly empty, still with
p congruent to 1 modulo L0. Write theta=chi bar(chi') and use the modulus
M above; M<=L0 Q^2. Since theta is nontrivial on H_M, (AP4) cancels
the prime-number main term after summing over the reduced classes in
H_M. Source Theorem 2.2, with a larger accuracy exponent to pay for
at most M residue classes, gives
\[
 \sum_{\substack{p\le t\\p\equiv1(L_0)}}\theta(p)
       \ll_K M t(\log t)^{-K}
 \qquad(t\ge y:=\exp(L^c)).
 \tag{AP6}
\]
It applies since M is a fixed power of log t once Q is a fixed power
of L. By partial summation the reciprocal-prime sum over ANY subinterval
of [y,exp(L^d)] is
\[
 \ll_K M(\log y)^{1-K}.
 \tag{AP7}
\]
The bound is absolute and requires no lower bound on the interval
length or on its harmonic mass. Choosing K sufficiently large gives
any prescribed L^{-A}. Remove the exceptional primes with error at
most their harmonic mass (AP2), divide by V_g>=v_-, then average the
other labels. This proves (AP5). The number of the other labels causes
no loss: their laws have total mass one and their character factors
have modulus one.

The argument needs only the arbitrary power savings supplied by the
stated source Theorem 2.2. The fixed-modulus prime number theorem can
supply exponentially small errors when constructing the AP bank
masses. It is not legitimate to assert exponential errors for every
growing character modulus merely from the weaker printed statement;
nor is such a strengthening necessary here. The new modulus differs
from that in the source by at most the fixed factor L0.

For the actual deletion p|N, if log N<=L^C for fixed C, then
\[
 \sum_{\substack{p\mid N\\p\ge \exp(L^c)}}1/p
 \le {\log N\over L^c}\exp(-L^c)
 \ll \exp(-L^c/2).
 \tag{AP8}
\]
This is uniform over all target valuations, so target-prime deletion
satisfies (AP2). The lower harmonic-mass bound must still be ensured
when the banks are first chosen, for example by choosing their
predeletion masses in a fixed interior subinterval of [v_-,v_+].

### 3. Kernel, positivity, operator norm and exact reproduction

On retained cells nu_j>=xi set, with the orientation matching the
source reproduction formula,
\[
 K_\Lambda(D,X)=
 \sum_{j:\nu_j\ge\xi}
 {1_{\log D\in I_j}1_{\log X\in I_j}\over\nu_j}
 \sum_{\eta\in R_Q}\overline{\eta(D)}\eta(X).
 \tag{AP9}
\]
D and X are independent products of the same slot type and law.
Conjugating/transposing the convention gives the equivalent frame
orientation if desired; it must not be replaced by an unconjugated
character sum.

The source's elementary estimates are unchanged:
\[
 |K_\Lambda|\le Q^2\xi^{-1},\quad
 \sup_X E_D|K_\Lambda(D,X)|\le Q^2,\quad
 N_\Lambda\ll_m1+TL^d/h.
 \tag{AP10}
\]
The last bound follows from the possible log-product range and not from
prime distribution.

Within one retained cell, the operator
Tf(X)=E_D K(D,X)f(D) is the frame operator
\[
 Tf=\sum_{\eta\in R_Q}\eta\,\langle\eta,f\rangle
\]
on the conditional probability space. In particular it is positive
semidefinite. Its nonzero eigenvalues are those of the Gram matrix.
For large x every bank prime exceeds Q^2 L0, so its diagonal entries
are exactly one. By (AP5) its off-diagonal entries are
O_A(xi^{-1}L^{-A}). Thus
\[
 \|T\|\le1+O_A(Q^2\xi^{-1}L^{-A})\le2.
 \tag{AP11}
\]
Restoring the original, rather than conditional, measure multiplies
both squared norms by the same nu_j. Different cells are orthogonal
blocks. Hence (AP11) holds on the whole product probability space,
without a cell-count or xi^{-1} operator-norm loss. Its adjoint and
transpose have the same norm. A tuple test need not be a function
of its product: pullback to the tuple space gives the same frame
operator, or equivalently factor through the conditional expectation
onto the product sigma-algebra, a contraction.

Let psi be any Dirichlet character modulo r<=Q, and let psi_0 in C_Q
be its primitive inducing character. Every bank product is a unit
modulo r, so psi(D)=psi_0(D). Let eta in R_Q represent psi_0's class.
Since D and X are 1 modulo L0, psi_0(D)=eta(D) and
psi_0(X)=eta(X) exactly. Its own summand in (AP9) therefore contributes
psi(X), and all the other summands satisfy (AP5). Consequently, for
every retained X,
\[
 E_D K_\Lambda(D,X)\psi(D)
   =\psi(X)+O_A(Q^2\xi^{-1}L^{-A}).
 \tag{AP12}
\]
This is the corrected, pointwise reproduction identity corresponding
to source (4.18). No factor equal to the size of an equivalence class
remains.

### 4. Full comparison lemma, with every source use identified

The AP version of source Lemma 4.3 has EXACTLY the same conclusion:
for arbitrary L2 tests on the entire independent label tuples,
\[
 \left|E f(X_B)g(X_A)\left[
 (X_AX_B)^{i\zeta}e(\alpha X_AX_B)
 -K_A(D_A,X_A)K_B(D_B,X_B)
  (D_AD_B)^{i\zeta}e(\alpha D_AD_B)\right]\right|
 \le L^{-A_0}\|f\|_2\|g\|_2,
 \tag{AP13}
\]
uniformly in all alpha and |zeta|<=L^{C4}. Here are all steps, so that
no additional prime-interval hypothesis is concealed.

1. **Product masses and arbitrary tuple tests.**
   Pairwise-disjoint banks and unique factorization, not a prime
   interval theorem, give source (4.3):
   rho_Lambda(n)<=L^kappa/n, with
   kappa=C0 log(m!)+m log max(1,v_-^{-1}).
   Restricting to AP primes or deleting primes changes none of this.
   Cauchy–Schwarz on each product fiber gives source (4.4) for arbitrary
   tuple tests: the associated integer coefficient sequence has
   squared norm <=L^kappa/U times the tuple squared norm on a dyad.

2. **Integer bilinear estimate.**
   Source Lemma 4.2 explicitly permits coefficients on arbitrary
   subsets of containing integer intervals. Its proof extends them by
   zero and estimates a full integer geometric sum. It therefore
   imposes no equidistribution of AP-prime products. The combined
   probability-space norm is source (4.6),
   \[
    \ll L^\kappa
      (V^{-1}+r^{-1}+\log(2r)/U+r\log(2r)/(UV))^{1/2}.
   \]
   The Mellin factor separates and preserves the two test norms.

3. **Dyads and rational approximation.**
   The selected products lie between exp(L^c) and
   exp(mC0 T L^d). There are at most L^{d+1} dyads for each.
   Dirichlet approximation with denominator cap floor(UV/S) is
   unchanged. For a resulting denominator r>Q, move the two kernel
   integrations onto the tests using (AP11), then apply the same
   integer bilinear bound. Log matching enlarges the dyads only to
   fixed constant multiples. The resulting bound is source (4.14):
   \[
    \ll L^\kappa
       (o(L^{-A})+Q^{-1}+L^{d+1}/S)^{1/2}.
   \]
   No prime theorem occurs in this part.

4. **Small denominators.**
   For r<=Q every product is a unit modulo r. The exact finite-group
   expansion
   \[
     e(uv/r)=\sum_{\psi\bmod r}c_\psi\psi(v),
     \qquad \sum_\psi|c_\psi|\le\sqrt{\phi(r)}\le\sqrt Q
   \]
   is still valid at these units. For v=D_AD_B its character factors
   separately. Replace each inducing primitive character by its
   representative in R_Q: this gives the same value on ALL D and X
   under consideration. If several characters now coincide, grouping
   coefficients cannot increase their total absolute mass. Therefore
   (AP12) reproduces the rational phase with pointwise error
   O_A(Q^{5/2}xi^{-1}L^{-A}), exactly source (4.19).

5. **Logarithmic multiplier and discarded cells.**
   On matched retained cells the Mellin/real-phase multiplier
   M(v)=exp(i zeta v)e(beta exp v) has derivative
   O(L^{C4}+S). The two log products change by less than 2h.
   The integrated absolute kernel bounds (AP10) give error
   O(Q^4 h(S+L^{C4})), exactly source (4.17).
   At discarded cells the comparison term is zero; their union has
   probability <=(N_A+N_B)xi. Cauchy–Schwarz and independence of the
   two test spaces bound the raw contribution by
   ((N_A+N_B)xi)^{1/2}||f||_2||g||_2, exactly source (4.16).
   This argument is valid even for tests concentrated in such a cell.

Thus the ONLY use of interval-prime equidistribution anywhere in the
source proof of Lemma 4.3 is replaced by (AP5). All other estimates
are unchanged, or improved by reducing the character list.

### 5. Parameter order and coefficient cost

The source's explicit noncircular hierarchy (4.20) still suffices.
After fixing m,A0,C4 choose
\[
\begin{split}
 B_0&>A_0+2(d+1)+10,\\
 q_*,s_*&>2(B_0+\kappa+10)+d+1,\\
 h_*&>4q_*+\max(s_*,C4)+B_0+10,\\
 x_*&>2B_0+d+1+h_*+10,\\
 D_*&>B_0+x_*+\tfrac52 q_*+10,
\end{split}
\]
and use Q=L^{q_*}, S=L^{s_*}, h=L^{-h_*}, xi=L^{-x_*}.
Finally choose the Siegel–Walfisz accuracy in (AP5) for this D_*.
The extra fixed factor L0 in the modulus and the errors (AP2) are
paid for at this last choice; they do not feed back into earlier
kernel complexity. None of these quantities involves J. Constants
and the threshold for x may depend on L0 and the fixed group data.

The number of cell-character summands in each kernel is a fixed power
of L. Each character used has conductor <=Q. For P=m|B| probes,
the number of signed comparisons is <=2^P<=L^{m C0 log 2}; and
\[
 \sum_\nu\sup|K_\nu|
 \le1+2^P Q^4\xi^{-2}\le L^{C_1},
 \qquad
 C_1>mC_0\log2+4q_*+2x_*+1.
 \tag{AP14}
\]
The raw pattern is included in the initial 1. Expanding both kernels
still produces only a fixed power of L terms. In particular the
finite-AP correction creates no exponential-in-J coefficient loss.

### 6. Symmetric-probe cancellation

For completeness, source Section 4.4 does not use interval primes,
distribution in residue classes, or the full character list. Under
(AP1) all labels in each group still have the same probability law,
so within-group permutation symmetrization is an orthogonal
projection, and the product-coordinate Efron–Stein decomposition is
unchanged.

Average the last ell input/target slots, and split the big groups into
two blocks. For m designated probes per group, a cardinality-j
component of a symmetric function on M=J+ell coordinates survives
the requirement that all m probes occur and all ell last coordinates
be absent in the fraction
\[
 (j)_m(M-j)_\ell/(M)_{m+\ell}
 \le2\ell^\ell m^{-\ell}=:\theta_m
\]
when M>=(m+ell)(m+ell-1). Tensorizing by group cardinalities gives
the source bounds ||F_full||<=theta_m^{s1/2}||F|| and similarly for H.
Choose m so that theta_m<1 and
-(c0/6)log theta_m>E_id+3.

For each pair of missing-probe sets U,V, free nonempty A subset U and
B subset V and use coefficient
\[
 -(-1)^{|A|+|B|}K_A(D_A,X_A)K_B(D_B,X_B).
\]
The asymmetric dependence leaves precisely the required unshared
mean-zero coordinates absent. After conditioning on the remaining
shared product C, (AP13) applies at alpha=Theta C, because it is
uniform for ALL real alpha. Conditional test norms are integrated by
Cauchy–Schwarz exactly as source (4.25). The subset signs sum to one
when U,V are both nonempty and cancel the raw term. The exceptional
full-probe components are summed before taking absolute values, so
their total costs only 2 L^{-E_id-3}||F||||H||, not a subset count.

There are at most 3^P approximation-error triples. Choose
A0>E_id+m C0 log3+10, as in the source. This proves the norm
L^{-E_id}. Finally choose J>=m and
J+ell>=(m+ell)(m+ell-1). The construction and its complexity were all
fixed beforehand, as required.

### 7. Scope and audit verdict

The finite-AP ideal-kernel transfer is proved under the explicit bank
and harmonic-mass hypotheses (AP1)–(AP2), with the character quotient
(AP3) as a necessary correction. It uses only the source's stated
Siegel–Walfisz input for the new orthogonality calculation. It applies
to Gaussian identity banks L0=4 and to arbitrary fixed identity banks
for finite abelian splitting subgroups. Target-prime deletion is
covered uniformly by (AP8).

This establishes the IDEAL-kernel interface, rather than treating it
as an immediate corollary of the printed interval-prime theorem.
Source Section 3 already allows these arbitrary prime groups, so no
further AP change to its hypotheses is needed once this ideal norm
and the unchanged coefficient-complexity bound are supplied. The
validity of the imported Section 3 transference proof, the
large-shift support argument, the prime-correlation input, and the
evaluation of the principal model remain separately identified
dependencies of a full binary-correlation theorem. No external
novelty, publication, independent human review, or unconditional
resolution of those other dependencies is claimed here.

## G. Large shifts, reflection and endpoint pairing

Proof module, 19November2026. This is a deduction from the exact imported k=1 graph results, not an independent verification of their proofs. It does not assert the complete binary norm asymptotic.

Primary source: *Weighted dilation graphs, smooth shifted primes and totient fibers*, September24,2026. Exact retained PDF sources/weighted-dilation-graphs-sep24-2026.pdf, SHA256351692885d112277e969901e75fe56eeee0d263d5d5847d0b980fe2c8c7744c0. Relevant statements: Definition3.3pp11–12, Lemma3.6pp13–14, Theorem3.5p12/proofcompletionp25, Corollary3.11pp26–27.

### LG1: the root moment has no small-shift dependence after deleting bank divisors

Keep the source prime groups, fixed mark/list parameters, signed comparison coefficients, damping and ideal-norm hypothesis. Let k be ANY nonzero integer coprime to every group prime. Define the physical operator A_k by the same edge displacement n'=n+kD. Then the averaged root moment of Theorem3.5 satisfies the same logarithmic bound as for k=1, uniformly in k.

Proof. The absolute row/column bound Lemma3.4 uses only the prime-list normalizations and damping; its proof explicitly does not use the magnitude of integer positions or their translations. It is uniform in k.

Expand a closed walk of the even length Nmom=2ceil(L^.52). For fixed labels/patterns/permutations its offsets are h_i=k s_i with s_i integers determined by the products D. Closure h_Nmom=0 is equivalent to s_Nmom=0, independent of nonzero k. The root dependence is exactly the nonnegative product(3.8) of divisibility and damping indicators n+k s_i=0modp.

Lemma3.6's CRT/Bonferroni proof is uniform in the numerical values of these offsets. Its modulus bound log Qact=O_J(Nmom*T*L^d)=o(L) and its number of terms do not involve the offsets. CRT interval errors areO(1), independently of their residue positions. Repeated offsets modulo a prime are collected into residues, which needs no bound on their size. Thus its relative replacement by independent Haar residues holds with exactly the source precision for every k.

In that independent model, replace the root residue at every group prime p by k^{-1}nmodp. Since p does not divide k, this bijects all active divisibility and damping conditions with the k=1 conditions. All label/pattern weights and the displacement-phase multipliers are unchanged on each closed path. Therefore every fixed-path residue expectation is EXACTLY the k=1 expectation. Sum with absolute coefficients using Lemma3.4 exactly as sourceLemma3.6 does. Restore the uniform replacement error, choosing its precision after the absolute path bound. The k=1 Theorem3.5 supplies the claimed moment. No extension of the memory argument or moment length is required. QED.

This argument also applies to a negative root interval of the same length: replace n by−n and k by−k, or use the uniform CRT interval count directly.

### LG2: the pairing bound loses no factor from a large k

Assume the same endpoint hypotheses as sourceCorollary3.11, including first-endpoint mark independence and the stated L2 norms. Then its pairing estimate holds uniformly for the k of LG1, irrespective of |k|.

Proof. Put Dmax=maxedge|D|=exp(O_J(T L^d)), and let H=|k|Dmax. SourceCorollary3.11 partitions the first endpoint interval I into intervals of lengthceil(H), restricts the second endpoint to their H-neighborhoods, and estimates the positive Gram matrix by its trace. Its only adverse factor is the square root of the numberO(H) of integer positions in a block, raised to1/R.

The operator A_k preserves each residue classmod|k|. Refine EVERY interval block by those residue classes. A refined block contains at most ceil(H/|k|)+1=O(Dmax+1) positions. The Gram matrix is block diagonal in these classes. The same proof now usesO(Dmax) rather thanO(H) in the trace step. The neighborhoods of distinct residue classes are disjoint; within one residue class the interval neighborhoods have bounded overlap. Thus the sums of endpoint squared norms and root traces are the same bounded multiples as in the source proof.

The resulting scalar factor is

  [O(Dmax)*sup_n|f_n|^2]^(1/(2R))
   <= exp(O_J((T L^d+sqrt(L))/L^.52))=O(1),

because d<.47<.52. SourceTheorem3.5 as extended inLG1 bounds the SUM of the traces over the whole original root interval; no CRT estimate in an individual residue classmodk is needed. The source Holder argument then gives exactly its pairing bound. QED.

### Bank primes dividing k and ideal-kernel scope

For |k|<=x^C, removing primes dividing k from any group changes its reciprocal mass by at mostO(log|k|*exp(−L^a)), hence smaller than every fixed negative power of L. LG1/LG2 themselves require only the resulting prime sets satisfy the source group/ideal-norm hypotheses. They do not by themselves prove the ideal-norm construction for intervals with deleted primes, or for banks restricted to a fixed splitting residue class.

The transference theorem permits arbitrary prime-group sets with the stated mass/range bounds. Module I proves the required ideal norm directly for fixed-AP banks with target-prime deletions, using their actual probability laws and corrected character frame. No unproved perturbative operator embedding is needed here; the source interval-prime construction is not invoked without that repair.

### LG3: the large-shift harmonic comparison needs no extra D-localization

Keep |k|<=C_k*x for a FIXED C_k. The shift may be comparable to x, including k~x and the reflected endpoint setup. Walsh Section7.5's endpoint cutoffs can be repaired directly, without inserting any product-localization and without changing the shared-label removal.

For a positive source n, let suppPsi be contained in [c,C], c>0. On the dyadic interval n in[X,2X), an ORIGINAL edge with nonzero weight Psi(n/(xD)) necessarily satisfies

    D <= 2X/(c*x),
    |n+kD| <= (2+2C_k/c)*X.                           (LG3.1)

Thus insert, BEFORE the Mellin inversion, a target cutoff g(n')=0 outside the fixed interval[-C'X,C'X], with C'>2+2C_k/c. It equals1 at every contributing target. There is no requirement that it contain every possible graph edge leaving the source interval.

After insertion, Mellin-invert Psi exactly as in source(7.17), and keep this cutoff in every pairing. Edges on which the cutoff removed a target had Psi=0 initially; their full Mellin integral is exactly zero. Consequently the identity(7.18) remains exact. The endpoint norm is

    ||g||_2 <= X^(1/2)*L^C,

by the same fixed-order divisor/list moment bound, now on an interval of lengthO(X). The source norm and the clipped mark-independent fourth-moment bound are unchanged. LG2, uniform in the Mellin phase, applies. The Schwartz Mellin tail is bounded by the same absolute row/column bound and the same endpoint norms. SummingO(L) source dyadic intervals proves the signed harmonic comparison to every fixed negative logarithmic power.

For reflection use n<0, m=|n| in[X,2X), Psi(m/(xD))/m, and extend the second endpoint coefficient by0 off its required sign. Then

    |n+kD| <= (2+2C_k/c)*X

again. Mellin inversion is in log(m/x), so the phase is |n/x|^(i*zeta); the graph moment proof is valid on a negative root interval by LG1. All endpoint norms depend only on |n|. This establishes the reflected comparison for k~x with a smooth bulk cutoff away from zero.

Crucially, no new cutoff depending on D has been added. The graph kernels still ignore shared labels exactly as sourceLemma7.1 requires. After obtaining the small graph pairing, perform Walsh(7.20)--(7.21) on the ORIGINAL unlocalized edge sum. Shared labels disappear with precisely the source probability calculation. Its only occurrence of fixed h in this stage was the assertion that reduced endpoints have sizeO(xD_C); under |k|<=C_k*x this remains true, and in the reflected bulk both endpoints are comparable to xD_C. Prime collision errors use at mostO(log x) prime factors and therefore retain their exp-small bounds. The absolute shared-label replacement is unchanged.

Proof assurance: LG3 is a direct modification of Walsh's written proof in Section7.5 pp37--41. It uses the already separated ideal-kernel construction and LG1/LG2; it does not claim those imported inputs were independently proved here. The source's stronger fixed-shift statement that its target cutoff contains EVERY edge leaving the source support is unnecessary. Containing all edges with nonzero harmonic weight before inversion suffices.

The full comparison also uses Module I for AP ideal normalization, Module E for conditional target-grid energy, and Module D for density/envelope reconstruction. There is no additional large-shift support loss in this harmonic comparison.

### Parameter order and endpoint moments

For the bounded core/model endpoint functions used here, the moment exponent in sourceLemma2.3 can be fixed BEFORE the graph parametersm,M. Indeed fix an ordered distinct list with productP. For fixedM the range bounds giveP=x^o(1); the sum over multiples ofP in any endpoint interval of lengthO(X), contained in[-C X,C X], isO(X/P). Every original/model/centered factor is bounded by one fixed constant oncex is sufficiently large in terms of the fixed Hankel depth. Summing over lists gives

    sum_lists (1/P)*product_g V_g^(-M) <=1.

Thus the exact-list endpoint r-moment isO_r(X), with noM-dependent log exponent; the eventual threshold can depend on fixedM and depth. The same applies to signed source intervals. This is the sourceLemma2.3 proof with its divisor weight replaced by a bounded invariant weight.

Consequently choose the reconstruction saving first, then the graph accuracy and comparison familym, and only then the target-energy saving after all resulting polynomial-log coefficient exponents are known. Choose the finite Hankel depth last to deliver that energy saving. Earlier endpoint bounds survive because the depth is fixed and the model remains bounded for all sufficiently largex. No circular choice between graph complexity and depth is needed.

## P. The all-conductor prime-factor estimate under ZF

Conditional analytic module, 19November2026. This proves the nonexceptional-character estimate needed by the proposed Walsh reflection argument. The external fixed zero-free-half-plane theorem is an EXPLICIT HYPOTHESIS, not independently verified here. The principal/gcd-record completion is proved in Module E, the centered low-frequency budget in Module C, and the large-shift transfer in Module G. The latter imports the specified k=1 graph theorem.

### Imported hypothesis ZF and exact source

(ZF) There is a fixed sigma_0<1 such that every primitive Dirichlet L-function has no nontrivial zero with real part>sigma_0, allowing the principal pole at1.

The recent source used by Walsh asserts sigma_0=7/8: OpenAI,*The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re(s)>7/8*,30September2026, Theorem1.1 printedp4. I read its actual theorem and scope, not the full194-page proof. Exact primary PDF retained sources/openai-quasi-rh-sep30-2026.pdf SHA2568fe93046f8cf5ef1ba5969c89addc02d76311adc4ee907509ff9cd96f7ec99e7, URL https://raw.githubusercontent.com/openai/math/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf . The following deductions remain conditional onZF unless that imported result's proof is accepted after its separate audit.

Walsh Lemma5.1 printedp20 (both retained v1 and v2 editions; v2 PDF sources/walsh-2610.40439v2.pdf, SHA2565fb593fdd03c79d90e5895177b895c1ce6afcfae660f2bf2f10312f0ff23bad1) derives the needed smooth prime estimate from this exact input. Its statement restricts the modulus tolog^B x, but its contour proof(5.2) uses only log q. Therefore the same proof permits q<=x^B. This extension is proved below with its uniform parameters; it is not obtained merely by changing the statement.

### PP: smooth prime characters at polynomial-size conductor

Let L=logx, P>=exp(L^a), q<=x^B, with a,B fixed positive. Let w be smooth supported[1/2,3], with fixed logarithmic derivatives bounded by powers ofL. Choose sigma_1 between sigma_0 and1. UnderZF,

  sum_p chi(p)w(p/P)p^(it)
   =delta_chi P^(1+it) integral_0^infinity
                    w(u)u^(it)/log(Pu)du
     +O(P^sigma_1 L^C (L+log(|t|+3))^C),                 (P1)

uniformly for every Dirichlet character moduloq and realt. The constant and exponent depend only on the fixed strip gap, derivative exponents,a,B. In particular for |t|<=x^B the error is P^sigma_1 L^C. This is uniform for imprimitive characters; their deleted Euler factors are included.

Proof. The standard local logarithmic-derivative formula, the functional equation and the zero-free gap give, uniformly on sigma_1<=sigma<=2,

  -L'/L(s,chi)-delta_chi/(s-1)
                    << log(q(|Ims|+3)).

For an imprimitive character the extra finite Euler sum is bounded byO(logq), because sigma_1 is fixed positive. For the principal pole subtract1/(s-1) before estimating.

Apply Mellin inversion with w(u)/log(Pu) to the von-Mangoldt sum and move the contour from2 tosigma_1. The only pole is s=1 for the principal character. On the new line the Mellin transform decays faster than every fixed power in Im(s)+t; its weighted integral against log(q(|Ims|+3)) is at mostL^C(L+log(|t|+3))^C. This gives the stated error, and the residue gives the main. Prime powers contributeO(P^(1/2)L^C), absorbed after increasing sigma_1 if necessary. This is precisely Walsh's proof with logq=O(L), rather thanO(logL). QED.

For banks p=1mod a fixed modulus L0, expand the residue indicator in the finitely many characters modL0. The prime main is zero unless chi is induced by a character of conductor dividingL0. All other characters satisfy the error bound(P1), uniformly. For Gaussian banks p1mod4, the exceptional pair is1,chi_-4. One should refine the ambient modulus to a multiple ofL0 first; this costs only a fixed factor for an ERROR bound.

The weakest useful substitute forZF is the conclusion of(P1) without its principal main: for every A, the normalized prime sum for all non-fixed-conductor characters must beO(L^-A), at bank lengths exp(L^a) and conductors up to polynomialx. Any uniform fixed zero-free half-plane suffices. The classical log(qt) zero-free strip does not give this at q~x, since its width times logP tends tozero. This module does not assert thatZF is necessary, only sufficient.

### PP-sharp: actual interval banks, with all-height Mellin control

The exponential bank error can be proved directly under ZF, without a smoothing-boundary error. Choose sigma_1>max(sigma_0,1/2). Uniformly for q<=x^B and X in[P/2,4P],

    theta_chi(X)=sum_(p<=X)chi(p)logp
       =delta_chi X+O(X^sigma_1 L^C).                  (PS1)

Here is the truncated-Perron detail. Replace X by the half-integer immediately after floor(X), changing a prime-weighted sum by at mostO(logX). Use c=1+1/logX and heightT=X² for the von-Mangoldt Perron integral. The elementary truncation error isO(log²X): in the range n~X use distance at least1/2 and sum(X/(T|n-X|))*logX; outside that range use absolute convergence on c. Shift the finite rectangle to sigma_1. ZF and the fixed strip gap give |-L'/L(s,chi)-delta_chi/(s-1)|<=C log(q(|Ims|+3)) throughout the rectangle. Thus the new vertical integral is O(X^sigma_1 log²(qX)), and each horizontal integral isO(X^c log(qX)/(T logX)). The only crossed singularity is the principal pole, with residueX. Imprimitive Euler terms areO(logq) on this fixed right half-plane. Prime powers costO(X^(1/2)log²X), absorbed in(PS1). The half-integer replacement changes delta_chi X byO(1), also absorbed. This proves(PS1) at every real endpoint, not just generic endpoints.

Partial summation against u^(it)/logu on ANY subintervalI of[P,2P] now gives

    sum_(p inI)chi(p)p^(it)
       =delta_chi integral_I u^(it)/logu du
        +O((1+|t|)P^sigma_1 L^C),                     (PS2)

for EVERY realt. The linear |t| factor is harmless: the product-cutoff Mellin transform in(P3) has

    integral |What_beta(t)|(1+|t|)dt <=L^C

when |beta|Y and the source Mellin phase are polynomial-logarithmic. This follows from fixed-order integration by parts and the stated logarithmic derivative bounds. Thus there is no Mellin truncation or smoothing-boundary error in the exponential NONEXCEPTIONAL estimate.

Restricting p=1modL0 expands into the finite fixed-character family as before. If chi is not induced from a character moduloL0, no constituent has a principal pole, so(PS2) has no main. Removing primes dividing the target costs at mostO(logN) terms on each interval, independently of t, and is absorbed by the displayed error. Consequently the literal sharp-bank assumption(P) of the companion variance lemma holds, with normalized error exp(-cL^a) at |t|<=L^D. Its all-height form grows only linearly in |t| and suffices for the full Mellin integral below.

### Parameters of the target-grid problem

Let Y=x H_free with H_free=x^o(1) and H_free>=exp(L^c), with fixed b<c. A small bank lies between exp(L^a) andexp(L^b). Let

  q <= x L^B,        Y/q >= H_free L^-B >> exp(L^b).

After fixed-period refinement the same inequalities hold up to fixed factors. The targetT is comparable tox; the actual Walsh resonance preimages have denominators dividingT r withr<=L^B, so these q cover them.

Consider an endpoint sequence F(n), supported on |n|asympY, with pointwise boundL^C. It is enough to treat positive n; reflection n->-n only inserts fixed character factors and reverses additive phases. Assume the endpoint has exactly one distinct prime from a selected small bank, its structured label coefficient ignores that small prime, and its invariant remainder is unchanged by multiplying by a bank prime. These are the actual complete-small-bank conditions of Walsh Sections6.4/7.3; for norm indicators use allowed split banks, or p1modL0 banks preserving the whole fixed filter.

Prime-square events may first be removed: their squared population is at mostY L^C sum_{p in bank}p^-2 <<Y L^C exp(-L^a). In residue energy moduloq, the maximum number of terms in a class isO(Y/q); hence their contribution isO(Y^2/q *L^C exp(-L^a)). This pays for all such exceptions on the needed scale. Remove those primes dividingq from the bank. For q dividingT times a polylogarithmic integer, deleting primes dividingT already accomplishes this for large x; their harmonic mass is negligible. The ideal-kernel effect of this deletion is a separate one-step operator normalization issue.

### NP: nonexceptional characters, with exact cofactor-length budget

Fix a gcd profile d|q and put r=q/d. On the surviving sequence n=d m, the extracted prime p is coprime toq and occurs to exponentone. Split its bank into dyadic rangesP<=p<2P. Then, after setting aside the small prime-square error,

  F(d m)=sum_{p in bank,p||m} [invariant/other-bank coefficient]

with the bank normalizerV^-1, bounded above and below. Write m=p v. For p in the chosen dyadic range, v is supported on a fixed enlarged interval of length

  M asymp Y/(dP).

The remaining coefficient gamma_d(v) is independent ofp and satisfies|gamma_d(v)|<=L^C. It includes the prohibition on any other small-bank prime dividingv, all other bank counts, the structured BIG-label coefficients and invariant remainder. The latter's p-invariance is essential. No assertion is made for an arbitrary coefficient depending on the extracted p.

For all nonzero profiles,

  M/r asymp Y/(qP) >>1,                                    (P2)

uniformly in d. Thus no loss appears merely because a target divisor d is large.

Let a smooth product cutoff W_beta have support in a fixed compact positive interval and fixed derivatives bounded by powers ofL. This includes the endpoint's original smooth cutoff, n^(i sigma) with|sigma|<=L, and additive phasee(beta n) when|beta|Y<=L^B. Mellin inversion separates it:

  S_(d,chi)(beta)
    = (1/V) integral What_beta(t) (Y/d)^(-it)
                 A_chi(t) B_(d,chi)(t) dt,

  A_chi(t)=sum_{p in dyadic bank}chi(p)p^(it),
  B_(d,chi)(t)=sum_{v asympM}gamma_d(v)chi(v)v^(it).          (P3)

The signs in t depend on Mellin convention and are immaterial. For each fixedj, the weightedL1 norm integral|What_beta(t)|(1+log(|t|+3))^j dt is at most a fixed power ofL. All these exponents are fixed afterB,C and the desired logarithmic saving.

Use the actual sharp bank directly. For a characterchi modulo r not induced from fixedL0, (PS2) gives

  |A_chi(t)| << (1+|t|)P^sigma_1 L^C.

The product cutoff has weighted Mellin L1 norm integral|What_beta(t)|(1+|t|)dt<=L^C, so this all-height estimate can be integrated directly. No smoothing strip is removed. The factorP^(sigma_1-1) beats every fixed log power becauseP>=exp(L^a).

Orthogonality at EVERY realt gives

  sum_{chi modr}|B_(d,chi)(t)|^2
     <=phi(r)(1+O(M/r))*sum_v|gamma_d(v)|^2
     << (r+M) M L^C <<M^2 L^C,                            (P4)

using(P2). This step does not take a supremum over large-conductor character sums of the remaining cofactor. It averages them exactly.

Minkowski in(P3), the weightedL1 bound for the cutoff and(P4) prove

  sum_{chi modr,nonexceptional}|S_(d,chi)(beta)|^2
     << (Y/d)^2 exp(-c L^a) L^C.                         (P5)

The finite dyadic prime partition costs only another fixed log power. The arbitrary-log conclusion follows by increasingx. There is no Mellin truncation error hidden here: the sharp-bank bound(PS2) grows only linearly in|t| and the smooth cutoff has an integrable (1+|t|)-weighted Mellin transform. One can instead truncate at polynomial height and pay the tail in the same vectorL2 norm using(P4).

### Passing from characters to all additive residues

For every residue profile gcd(n,q)=d, setr=q/d. Its residue-energy contribution is exactly

  (1/phi(r)) sum_{chi modr}|S_(d,chi)(beta)|^2,

whereS includes the restriction(v,r)=1 (automatically enforced bychi). Therefore(P5), summed over d|q, is

  << (Y^2/q) exp(-c L^a) L^C,

since r/phi(r) and sum_{d|q}1/d are bounded by fixed powers oflogq. This is stronger than the needed sparse scaleY^2/(q logY). It is uniform in every target valuation.

The finitely many exceptional fixed-conductor characters are NOT discarded. Their contribution needs the centered density model. The low-frequency recordwise cancellation is supplied by the independent module mentioned at the start. Its required uniformity is as follows: for each retained S-smoothd, its principal/field branch must have size at most(Y/d)L^-A, for any prescribedA after enough Hankel terms. Large-good-prime gcd profiles and records n_S>exp(sqrtL) are disposed of in RESIDUE ENERGY by the same max-class-count times globalL2-population bound, not merely by an L1 error.

For norm indicators the remaining record coefficients are bounded bylog powers. Summing records A=dB gives sum_{A Ssmooth,d|A}1/A=(1/d)V^-1, so the independent recordwise proof indeed has the required1/d scaling. All removed/selected record subsets may depend onq. This is why no SD theorem at conductorq is required for the principal branch.

Subject to that centered-exceptional-branch lemma, the final conclusion is

  E_q(beta):=sum_{a modq}|sum_{n=a modq}F(n)e(beta n)|^2
           <<_A (Y^2/q)L^-A                              (P6)

for everyq in the target range and|beta|Y<=L^B. The model centering and full invariant factor must be kept; (P6) is false for an arbitrary uncentered sparse indicator.

### From relative residue energy to the full target resonance set

LetT asympx. The bad set for the free-product multiplier m(Ttheta) is contained in the union, for1<=r<=L^B, of arcs centered at

  (j+a/r)/T,      0<=j<T, (a,r)=1,

with widthO(L^B/(T H_free)). Each center lies on the grid of modulusTr, and refining by fixedL0 is harmless. BecauseY=x H_free andT~x,

  Y*arcwidth <<L^B.                                      (P7)

Thus only polylogarithmic Mellin/additive variation is needed; the fixed-shift source's t-range up tox is not forced here.

AdditiveParseval on each grid gives exactly

  sum_{a modq}|Fhat(a/q+beta)|^2 = q E_q(beta).

Integrating over the displayed width and applying(P6), then summing the polylogarithmically many r, gives

  integral_{target major set}|Fhat(theta)|^2 dtheta
                       <<_A Y L^-A.                     (P8)

The numberT of arcs has canceled through residue-energy normalization; no supremum-times-number-of-arcs bound was used. Together with the second endpointL2 norm and the small multiplier off these arcs, (P8) supplies the missing comparison-correlation estimate in WalshSection6, conditional onZF and the precise finite-AP graph/centering interfaces. This module does not alone assert a complete unconditional binary theorem.

### Separate logical inputs

1. ZF is indispensable to the present all-character prime-factor estimate; its recent primary proof has NOT been verified here.
2. Restricted banks use the corrected ideal character quotient/frame proved in Module I.
3. The full reflected graph application keeps signed cutoffs; LG3 of the large-shift graph module now proves the endpoint support repair directly, with no extra D-localization or changed shared-label coefficients.
4. The independent exceptional-character record lemma supplies the finite main cancellation. This file proves the nonexceptional/orthogonality part, not a substitute for that lemma.

## E. Target-grid energy, gcd profiles and fixed filters

Proof module, 19 November 2026. Conditional ONLY on the explicitly stated all-conductor prime-character estimate(P) below, plus the classical centered rough-model lemma proved in PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md. This module does NOT prove the new all-conductor zero-free theorem or the large-shift graph comparison. It checks the entire arithmetic/variance interface needed if those independent inputs are available.

The bank factorization is Walsh arXiv:2610.40439v1, actual(6.19)--(6.22), but the second-moment argument here is a target-modulus adaptation, not a statement from that paper.

### 1. Parameters, fixed filters and actual coefficients

Let x tend to infinity, L=log x, Y=x^(1+o(1)), and

    y=exp(L^(2/5)), S={p:p<=y}, E=exp(sqrt L),
    P_min>=exp(c L^(1/8)), P_max<=exp(C L^(1/4)).

Fix an even modulus L0 containing the quadratic field character chi's conductor, the core semigroup C, and U a subset of H=ker chi in the units modulo L0. Replace U by U intersect H if necessary. Banks consist only of primes1modL0, all in S. Every bank has reciprocal-prime mass bounded above and below by fixed positive constants; there are O(log L) banks and their list lengths are bounded by one fixed integer. This permits the source's exact normalized list weights W_t^b, all bounded by L^C. Fix a selected SMALL bank g with t_g=1, and assume b is independent of its prime slot, exactly as in Walsh's Proposition6.1. All bank constants may depend on the fixed data.

Use the J-term density M=M_(J,S) from the centered-density proof in this support. For A=n_S and v=n/A define

    f(n)=1_C(n),
    D_U(n)=2*1_U(n)*f(A)*M(v)/V,
    alpha_U,P(n)=1_U(n)*f(n_(S\P))
         [ f(v) - (1+chi(v))*M(v)/V ].                  (G1)

Here V=product_(p<=y)(1-1/p). Since f(A)>0 implies chi(A)=1 and U subset H, the bracket in(G1), after multiplication by its other factors, is exactly f(v)-2M(v)/V. The displayed version retains the two-branch form useful for centering. Multiplication by any bank prime leaves alpha_U,P invariant because that prime is1modL0 and belongs to S.

Let

    F(n)=psi(n/Y)*n^(i sigma)*alpha_U,P(n)*W_t^b(n),      (G2)

where psi is smooth on one fixed compact subset of (0,infinity), its derivatives are bounded by fixed powers of L, and |sigma|<=L^C. The model's smooth bounded extension at v<x^(1/2) is used; consequently |F(n)|<=L^C. The same is true after deleting the selected bank from the list weight.

Let q be any positive modulus such that

    Y/q >= P_max*L^K                                  (G3)

for every fixed K needed in the application. A sufficient typical application is Y~xD1D2, q<=xL^B, and D1D2>=exp(cL^(1/3)). The result has no restriction on the prime factors or valuations of q.

For a sequence G supported on n~Y put

    E_q(G)=sum_(a modq) |sum_(n=a modq)G(n)|^2.         (G4)

Goal: for every B>0, choose a FIXED J sufficiently large so that

    E_q(F) <<_B (Y^2/q)*L^-B.                          (G5)

All constants may depend on fixed filter/bank data, B and the cutoff seminorm bounds.

### 2. The only new analytic prime input

Assume(P): for q1<=C*Y, every character eta modq1 which is NOT induced by a character modulo L0, every dyadic interval at the selected bank scale, and every smooth bank/dyadic weight w with logarithmic derivative bounds,

    |sum_(p in bank g, p not dividing q)
         w(p/P)*eta(p)*p^(it)|
       << P*L^C*exp(-cL^(1/8)),                        (P)

uniformly for |t|<=L^D for every required fixed D.

An all-conductor zero-free half-plane sigma>sigma_*<1, together with the standard log-derivative estimate a fixed distance to its right, implies(P) by the same Mellin contour calculation as Walsh Lemma5.1. For primes1modL0 expand their indicator into the finitely many characters modulo L0; no constituent is principal precisely for the eta covered by(P). Removing primes dividing q costs at most O(log q) terms, negligible after division by P. Module P supplies this prime-input implication under ZF, including sharp intervals and all required frequencies. This proof does NOT silently import it from a theorem restricted to logarithmic conductors.

Only LOGARITHMIC Mellin frequencies are required for the target-grid corollary in Section7, unlike the fixed-shift high-frequency range in the original Walsh application.

### 3. Modulus refinement and the exact all-residue character identity

Put Q=lcm(q,L0). Refining residue classes and applying Cauchy gives

    E_q(F) <= (Q/q)*E_Q(F) <= L0*E_Q(F).               (G6)

Thus no redundant character-frame normalization is needed for this UPPER ERROR estimate. The leading constant is not calculated using(G6).

The support of F is coprime to L0. For d=(n,Q), only d coprime to L0 occur, and Q_d=Q/d remains divisible by L0. Standard multiplicative-character orthogonality gives the EXACT identity

    E_Q(F)=sum_(d|Q) 1/phi(Q_d) *
        sum_(eta modQ_d)
          |sum_m F(d*m)*eta(m)|^2.                     (G7)

Characters vanish on nonunits, so the inner sum imposes (m,Q_d)=1 automatically. It describes exactly the gcd profile d, including arbitrary prime powers. Every character modulo L0 has a unique induced character modulo Q_d; these O_(L0)(1) characters are called exceptional. All other characters satisfy(P).

### 4. Negligible point sets in the residue-energy norm

For any point set B and sequence G with |G|<=L^C,

    E_Q(G*1_B)
       <= (1+C Y/Q) * sum_(n~Y)|G(n)|^2*1_B(n).        (G8)

This is Cauchy inside each residue class. Since Y/Q>>P_max, the factor is O(Y/Q).

The following sets therefore have total energy
O((Y^2/Q)*L^C*exp(-cL^(1/10))):

(a) n_S>E. The elementary tilted-record moment (Walsh Lemma2.2 with order1, or its direct Rankin proof) bounds the population by Y*L^C*exp(-cL^(1/10)).

(b) n divisible by a prime p|Q with p>y. Its population is at most
Y*sum_(p|Q,p>y)1/p <<Y*L/y.

(c) n divisible by a selected-bank prime p|Q. Its population is at most
Y*sum_(p|Q,p>=P_min)1/p <<Y*L/P_min.

(d) n divisible by p^2 for a selected-bank prime p. Its population is at most
Y*sum_(p>=P_min)1/p^2 <<Y/P_min.

These are absolute upper bounds; no distribution of norm values or target valuations is assumed. If the cutoff's support requires a fixed multiple of Y, the same estimates change only fixed constants.

Let F_dagger be F with(c),(d) removed. It retains EXACT selected-prime factorization. Let F_good additionally remove(a),(b). By(G8),

    E_Q(F-F_good)+E_Q(F-F_dagger)
       <<(Y^2/Q)*L^C*exp(-cL^(1/10)).                  (G9)

Do NOT impose n_S<=E before the nonexceptional prime factorization: that would couple p and m through p*m_S<=E. We use the full F_dagger for nonexceptional characters and F_good only for centering. Projections onto either part of the character identity(G7) are contractions in this norm, so this distinction costs at most the tail(G9).

### 5. Nonexceptional character energy

In a gcd profile d, every remaining selected-bank prime is coprime to Q and occurs once. Exact invariance and the selected list-slot condition give Walsh's factorization, with all OTHER banks retained in one coefficient gamma_d(m), |gamma_d(m)|<=L^C. There is NO need to expand another bank product a' and shorten the remaining sum.

After dyadically localizing p~P and Mellin-separating psi(dpm/Y), each character sum is a polynomially weighted combination of

    [sum_(p in g, p not dividing Q) w(p/P)eta(p)p^(it)]
       * [sum_(m~M) gamma_d(m)eta(m)m^(it)],
    M=Y/(dP), |t|<=L^D.                               (G10)

The discarded Mellin tail is smaller than every prescribed negative logarithmic power in ENERGY: truncate the transform at a sufficiently large fixed power of L and use the trivial character-energy bound(G8). Its L1 norm and the number of dyadic pieces cost only L^C.

CRUCIAL UNIFORM RATIO:

    M/Q_d=Y/(Q P)>>1                                  (G11)

independently of d. This is why large target valuations cause no loss.

Elementary orthogonality, with no multiplicativity hypothesis on gamma_d, gives

    1/phi(Q_d) * sum_(eta modQ_d)
       |sum_(m~M)gamma_d(m)eta(m)m^(it)|^2
      <= (1+C M/Q_d)*sum_(m~M)|gamma_d(m)|^2
      << M^2/Q_d * L^C.                               (G12)

Apply(P) to the first factor of(G10), and Minkowski to the fixed-power-log dyadic/Mellin decomposition. The nonexceptional energy of profile d is

    << (Y^2/(Q d))*L^C*exp(-cL^(1/8)).                 (G13)

Finally sum_(d|Q)1/d <= Q/phi(Q)<<loglog(3Q)<<L. This establishes the same bound for the total nonexceptional energy. This step permits arbitrary multiplicative coefficients in every OTHER bank; only the selected small-bank slot has to decouple.

### 6. Exceptional characters: exact centering after the gcd split

We use F_good. Since it has no prime divisor of Q outside S, its exact gcd d=(n,Q) is determined solely by A=n_S. Thus the profile condition only selects records

    A S-smooth, A<=E, gcd(A,Q)=d.

Write n=A v with (v,S)=1. For an exceptional eta induced by eta0 modL0,

    eta(n/d)=eta0(A/d)*eta0(v)

on this profile. The phase eta0(A/d) has modulus1. Every other record weight has modulus <=L^C, while the residual periodic test is

    b_A(v)=1_U(A v)*eta0(v),

of FIXED period L0. The bracket is exactly the centered two-branch expression from(G1).

Apply the recordwise proof of PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md, equations(W6)--(W8), at length Y/A. Its normalized error can be smaller than every desired fixed logarithmic power by choosing J. Uniformly in d, the remaining record sum gains precisely1/d:

    sum_(A S-smooth, d|A)1/A
      = (1/d)*sum_(B S-smooth)1/B
      = 1/(d V) << L^(2/5)/d.                         (G14)

No d_K(d) loss appears: for NORM INDICATORS the actual record/list weight is bounded by L^C independently of A,d. Consequently, for every A0>0, a fixed sufficiently deep model gives

    |sum_m F_good(d*m)eta(m)| <<_(A0) (Y/d)*L^-A0.      (G15)

Restrictions gcd(A,Q)=d, A<=E only delete records from an already ABSOLUTE recordwise bound. The deleted large-Q primes and large records were already paid in(G9); their omission does not alter the principal model.

After dividing by phi(Q_d), summing the O_(L0)(1) exceptional characters and d, the energy is

    << (Y^2/Q)*L^(-2A0) *
           sum_(d|Q) (Q_d/phi(Q_d))/d
    << (Y^2/Q)*L^(-2A0+C).                            (G16)

Choosing A0 after B and adding(G9),(G13) proves(G5).

This proves all gcd layers, including d close toQ; it does not assume that d is bounded, squarefree, or prime to the target. The retained-record calculation needs only d<=E, which follows automatically on F_good; every other layer was paid as an energy tail.

### 7. The target-grid major-energy corollary

Let q~x times at most a fixed power of L and suppose Y/q>>P_max*L^K. Let |beta|<=L^B/Y. The replacement

    psi(u) -> psi(u)*e(beta*Y*u)

has only fixed polynomial-log derivative costs. Thus(G5) is uniform in beta after choosing J accordingly.

Discrete additive orthogonality says exactly

    sum_(h modq) |Fhat(h/q+beta)|^2
       = q*E_q(F_beta).

Integrating gives, for any prescribed A,

    sum_(h modq) integral_(|beta|<=L^B/Y)
           |Fhat(h/q+beta)|^2 d beta
       << Y*L^-A.                                    (G17)

Any union of target-grid arcs is bounded by the same sum. A polynomial-log family of moduli q=x*r, r<=L^B, is also harmless by first demanding a larger A.

For a multiplier m(k theta) with k~x and free product H'=D1D2, the major arcs from approximating k theta have centers on grids of denominator k*r and physical width O(L^B/(kH')). When Y~xH' and k~x, this width is O(L^B/Y). Therefore ONLY logarithmic Mellin frequencies enter this major-energy estimate. The original fixed-shift Walsh range |t|~x*polylog is not the range required by(G17).

This corollary is an upper error bound and introduces no singular-series normalization. Any leading term must still be evaluated from the full density model.

### 8. Scope

Proved in this module: assuming(P), the actual bounded norm-indicator/list coefficients have arbitrarily small target-grid energy after full two-branch/fixed-period centering, uniformly in every target gcd and prime-power profile. There is no residual target-conductor average gap at this particular interface.

The main proof combines this energy estimate with Modules G/I and D. The remaining deep inputs are the explicitly assumed all-conductor zero-free statement and the imported k=1 physical graph moment. Neither is inferred from(G5).


### 9. Exact extension to an arbitrary fixed prime-support subgroup

The proof above is unchanged for a subgroup H of (Z/L0)^* of any fixed index kappa>=2, with f supported on primes in H and rho=1/kappa, provided the centered rough lemma of PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md is used.

Replace(1+chi(v)) in(G1) by

    sum_(eta in Hperp)eta(v)=kappa*1_(v modL0 inH),

and replace the model's factor2 bykappa. For U subsetH and f(A)>0, the condition A*v inU implies v inH, so D_U=kappa*1_U(n)f(A)M_R,J(v)/V exactly. Its size is L^(-(3/5)(1-rho)), hence bounded uniformly. All record/list coefficients remain bounded by L^C.

The low-frequency lemma cancels every Hperp branch and permits any fixed-L0 periodic b_A. The exceptional target characters are still ALL characters induced modulo L0, not merely Hperp; after the U filter their recordwise periodic tests are covered by that same lemma. The target prime input(P), using banks p1modL0, and every character-orthogonality, gcd-profile, energy-tail and grid-integration estimate are otherwise literally identical.

Thus(G5),(G17) are valid for every fixed positive-density prime-support subgroup under the same all-conductor prime input and the explicitly written generic-H centered model. This is an extension of this PROVED-CONDITIONAL variance module; the independent graph, zero-free, density-model and reconstruction modules still must be assembled and audited.

## C. Expanded centered low-frequency budget

Proof module, 19 November 2026. This is a new module adapting the low-frequency calculation in Walsh, *Shifted convolutions of divisor functions*, arXiv:2610.40439v1, actual printed pp19--23. It is not a deduction of a variable-target binary correlation, graph comparison, or reflected sum theorem. The source's explicit further-results paragraph on printed pp5--6 leaves sums-of-two-squares indicator correlations to a future extension.

Primary PDF: https://arxiv.org/pdf/2610.40439v1 .
Local PDF/text: sources/walsh-2610.40439v1.{pdf,txt}.

### 1. Setup and the indispensable second principal branch

Fix an imaginary quadratic primitive character chi and a fixed even bad modulus L0 divisible by its conductor. Let C be the positive integers all of whose prime factors p satisfy p not dividing L0 and chi(p)=1, and let f=1_C. Put

    L=log x, y=exp(L^(2/5)), S={p:p<=y},
    V=product_(p<=y)(1-1/p), E=exp(sqrt L).

All primes dividing L0 belong to S once x is large. Write

    D_S(s)=sum_((v,S)=1) f(v)/v^s
          =product_(p>y, p allowed)(1-p^-s)^-1.

Near s=1, on the branch positive to its right,

    D_S(1+z)=z^-1/2 G_S(z),                              (W1)

where G_S is analytic and nonzero in a disk |z|<=c/log y. Indeed squaring D_S gives zeta(s)L(s,chi), the usual fixed bad Euler factors, the absolutely convergent inert correction, and the squares of the allowed factors removed below y.

For fixed j,

    g_(S,j)=[z^j]G_S(z),
    G_S(0) asymp_(chi,L0) (log y)^-1/2,
    |g_(S,j)| <<_(j,chi,L0) G_S(0)(log y)^j.             (W2)

The derivative estimate follows either from Cauchy on a disk of radius c/log y or by differentiating the finite removed-prime product; sum_(p<=y)(log p)^j/p <<_j (log y)^j. The ratio G_S(z)/G_S(0) stays bounded above and below on a sufficiently small such disk.

For a character eta modulo a logarithmic-sized modulus, with its conductor's primes in S,

    D_S(s,eta)=sum_((v,S)=1) f(v)eta(v)/v^s.

Its square factors as L(s,eta)L(s,eta chi) times the fixed/absolutely convergent corrections and the removed factors. Consequently the ONLY branches at1 are eta=1 and eta=chi. In both these cases D_S(s,eta)=D_S(s) EXACTLY, because f(v)>0 implies chi(v)=1. All other twists have no branch at1.

Thus a uniform-residue rough density M/V alone is wrong. The model must be

    (1+chi(v)) M(v)/V.                                  (W3)

If the tested modulus does not contain the conductor of chi, enlarge it by that fixed conductor. This costs only a fixed factor.

### 2. A finite Hankel density is enough

For an integer J>=1 define, on t in [x^(1/2),x^2],

    M_(J,S)(t)=sum_(0<=j<J)
        g_(S,j) (log t)^(-j-1/2)/Gamma(1/2-j).            (W4)

This is a density, not a summatory function. Its leading term is positive, and all others are relatively O_J(L^-3/5). Hence it is positive for large x and

    M_(J,S)(t) asymp G_S(0)L^-1/2 asymp L^-7/10,
    M_(J,S)(t)/V asymp L^-3/10.                         (W5)

Its logarithmic derivatives satisfy corresponding fixed polynomial bounds. Extend it below x^(1/2) smoothly and nonnegatively with M/V << L^O(1), as in Walsh's harmless low-range extension. The later low-frequency argument only evaluates it where t=x^(1+o(1)).

The following smooth Selberg--Delange estimate holds uniformly for Q=x^(1+o(1)), smooth psi supported in one fixed compact subset of (0,infinity), and archimedean twists |tau|<=L^D:

    Q^-1 sum_((v,S)=1) f(v)eta(v)psi(v/Q)v^(i tau)
      = delta_(eta in {1,chi}) integral_0^infinity
            psi(u)(Qu)^(i tau)M_(J,S)(Qu)du
        + O_J( G_S(0)L^-1/2 L^(-3J/5) )
        + O_A(L^-A).                                   (W6)

Here the second error may be made smaller than every prescribed fixed power; constants and the eventual threshold can depend on fixed D, the logarithmic modulus bound, and finitely many cutoff seminorms. If ||psi||_infinity is allowed a fixed logarithmic factor, multiply the first error by that same factor.

Proof details for the required uniformity:

* Work in the rectangle sigma>=1-c/log y and |Im s|<=L^M, with M a sufficiently large fixed constant. Classical zero-free regions, together with Siegel's ineffective lower bound for a possible real exceptional zero, give a zero-free rectangle for the finitely logarithmic-sized conductor family once x is sufficiently large. Explicitly choose Siegel epsilon with epsilon*B<2/5 when the conductors are <=L^B. Nonreal zeros are excluded by a width >>1/log L, much larger than L^-2/5. This makes the fractional branches unambiguous. No target-sized conductor estimate is being invoked.

* On this rectangle the removed Euler product is polynomially bounded in L: p^(c/log y)<=exp(c), and its logarithm is O(log log y). The L-function factors have polynomial bounds in the logarithmic height and modulus. Mellin cutoff decay makes the top and bottom errors O_A(L^-A), choosing fixed M and then sufficiently many derivatives. The vertical contour is O(L^C exp(-c L^3/5)).

* For eta=1 or chi, deform around the branch cut z in [-c/log y,0]. Taylor-expand ONLY G_S(z) through degree J-1. Do not Taylor-expand the cutoff's Mellin transform, so the phase u^(i tau) is retained exactly. The Taylor remainder is
  O_J(G_S(0)(log y)^J |z|^J).
  Its Hankel integral is bounded by
  O_J(G_S(0)(log y)^J(log Q)^(-J-1/2)),
  which is precisely the first error in(W6). The cutoff Mellin transform along this real cut is bounded by its L-infinity norm independently of tau.

* For each Taylor monomial, the full Hankel integral equals
  (log t)^(-j-1/2)/Gamma(1/2-j). Extending the finite cut to the full Hankel contour costs exponentially little. Interchanging the compact u-integral gives exactly the integral of(W4).

* For every other eta there is no singularity at1, so the contour shift has only the exponentially small/logarithmically arbitrary error.

A FINITE J chosen after the desired saving therefore suffices. There is no need to assert exact cancellation of an infinite asymptotic expansion.

### 3. Evaluating the model itself

The ordinary S-rough progression estimate is Walsh's Lemma5.2 proof with divisor order a=1 (it is also an immediate Mellin contour shift):

    Q^-1 sum_((v,S)=1) eta(v) phi(v/Q)
       = delta_(eta=1) V integral phi(u)du
         + O_A(L^-A),                                  (W7)

uniformly for all characters of logarithmic-sized, S-supported conductor and smooth cutoffs with fixed polynomial logarithmic derivative bounds.

Apply(W7) separately to eta and eta chi with

    phi(u)=psi(u)(Qu)^(i tau) M_(J,S)(Qu)/V.

For eta=1 exactly one of eta,eta chi is principal; the same is true for eta=chi. For all other eta neither is principal. Hence

    Q^-1 sum_((v,S)=1) eta(v)psi(v/Q)v^(i tau)
           (1+chi(v)) M_(J,S)(v)/V
      = delta_(eta in {1,chi}) integral
           psi(u)(Qu)^(i tau)M_(J,S)(Qu)du
        + O_A(L^-A).                                   (W8)

Subtracting(W8) from(W6) cancels BOTH complete retained branches. The remaining relative SD truncation can be as small a fixed logarithmic power as required.

### 4. Full periodic and smooth-record lemma, including the log budget

Fix K,C,D>=0. For every S-smooth A let b_A be q-periodic with |b_A|<=1, where q<=L^C. Let psi_A have fixed compact support and derivatives bounded by fixed powers of L uniformly in A. Define

    S_A(tau)=sum_((v,S)=1) b_A(v)psi_A(Av/R)(Av)^(i tau)
      [ f(v) - (1+chi(v)) M_(J,S)(v)/V ],                (W9)

where R=x^(1+o(1)). For every prescribed B>0 there is a FIXED J such that

    R^-1 sum_(A S-smooth) d_K(A)|S_A(tau)| <<_B L^-B,   (W10)

uniformly for |tau|<=L^D. (Take d_0=1 if needed; any fixed divisor-bounded record majorant also works.)

To prove this, truncate A<=E. For A>E, the absolute f and model terms are bounded by a fixed logarithmic factor times d_max(K,1)(Av) 1_( (Av)_S>E). Walsh's elementary tilted-moment Lemma2.2 gives O(L^C exp(-cL^1/10)) after division by R. No norm-correlation input is used here.

For A<=E put Q=R/A=x^(1+o(1)). Enlarge q by the fixed conductor of chi and expand b_A into multiplicative characters on the S-rough units. Every normalized coefficient has modulus <=1, so the sum of absolute coefficients is <=phi(lcm(q,cond chi)) << L^C. The phase A^(i tau) has modulus1. Equations(W6)--(W8) bound each normalized twist by

    O_J(L^-7/10 L^(-3J/5)) + O_A(L^-A).

Finally

    sum_(A S-smooth) d_K(A)/A
      = product_(p<=y)(1-1/p)^-K << L^(2K/5).           (W11)

Thus, when the cutoffs have bounded sup norm, the explicit first error in(W10) is

    O_J( L^(C+2K/5-7/10-3J/5) ).                       (W12)

Choose, for example,
    3J/5 > B+C+2K/5+2.
A fixed extra cutoff sup-norm factor is absorbed on the right. Choose the arbitrary contour error saving after B+C+2K/5. This proves(W10) with no unaccounted logarithmic loss.

This is the exact analogue needed in Walsh Section6.3. It supplies arbitrarily strong low-frequency cancellation after centering, rather than the finite fixed log saving of an UNCENTERED inverse norm indicator.

### 5. Split prime banks and the separate high-frequency input

Take every active bank to consist of good split primes. Multiplication by such a prime preserves chi, f, and the S-rough cofactor. The invariant-factor construction therefore remains EXACT:

    D_f(n)=f(n_S)(1+chi(n_(S')))M_(J,S)(n_(S'))/V,
    xi_f=f-D_f,
    alpha_P(n)=f(n_(S\P))
        [ f(n_(S'))-(1+chi(n_(S')))M_(J,S)(n_(S'))/V ].

Then xi_f(n)=f(n_P)alpha_P(n), and alpha_P(pn)=alpha_P(n) for every bank prime p. There is no hidden prime-support restriction on the rough centered model: the model is deliberately allowed on nonsupported rough integers and cancels by its character means.

A split-prime version of Walsh Lemma5.1 follows by writing

    1_(chi(p)=1) eta(p) = (eta(p)+(eta chi)(p))/2

away from the fixed bad primes. There are two possible principal twists, each with half the ordinary principal integral. Both decay rapidly at high Mellin frequency, and the source's contour error is unchanged up to the fixed conductor enlargement. This statement imports the SAME explicit zero-free-region input as Walsh's Lemma5.1; that imported 2026 input has not been reproved in this module. It does not require a new target-sized character modulus if all major denominators stay logarithmic.

### 6. What is and is not closed

Closed here: the fractional pole and the second principal branch are NOT an obstruction to the low-frequency centering. A finite, sufficiently deep Hankel model produces arbitrary fixed logarithmic savings with the full smooth-record and periodic-test budget(W10).

Not closed here: norm-indicator envelope/bank reconstruction at the required N/log N scale; the graph comparison and its uniformity for a reflected large target; every comparison's multiplier/major-denominator range; the density-model binary main term with its full local constant; and independent verification of imported 2026 graph/zero-free results. Walsh's source only states fixed-shift integer divisor correlations, and explicitly leaves norm indicators and large-shift uniformity as further work.

This module therefore cannot by itself promote the desired every-target unweighted primitive norm leading asymptotic.

