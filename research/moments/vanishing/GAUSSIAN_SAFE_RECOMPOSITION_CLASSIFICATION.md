# Safe Gaussian recomposition and the unique safe shift

Let K be a characteristic-zero field and E the normalized unweighted
Gaussian Wick functional on K[X], with nonsingular covariance. For a
polynomial phase P define its ORIGINAL marked periods and images by

    Phi_P[Q](t)=E[Q exp(tP)],
    V_P=Phi_P[K[X]], M_P=K[t] V_P.

A phase is safe if every Phi_P[Q] is polynomial; the eventual cutoff may
depend on the original marker Q. It is pure if Phi_P[1]=1.

**Safe recomposition theorem.** Let f in K[T] be nonconstant. Then

    f(P) is safe
      iff there exists lambda in K with
          f(lambda)=0 and P-lambda safe.            (R1)

The admissible lambda is unique; this does not say that f has only one
algebraic root. It is the same unique safe shift for every nonzero safe
recomposition of a fixed P.

The proof uses the exact
[finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
It assumes neither universal safety implies purity nor primitive vacuum
for every finite image. The vacuum p(t)exp(lambda t) is allowed to have
a nonconstant polynomial factor p.

## 1. A common original differential annihilator for each marker

First suppose P is nonconstant, and put d=deg f. The substitution
K[T]->K[X], T->P, is injective, so K[P] is freely generated over K[f(P)]
by

    1,P,...,P^(d-1).                                (R2)

Indeed the leading coefficient of f is invertible, so its degree-d
relation reduces every higher power to this span. For independence, the
nonzero terms A_j(f(T))T^j have highest degrees d deg A_j+j, which are
distinct modulo d. Their largest degree therefore cannot cancel.

Assume f(P) is safe. Fix ONE original marker Q. For the finitely many
markers QP^j, 0<=j<d, choose a common integer N_Q such that

    E[QP^j f(P)^n]=0 for every n>=N_Q and every j<d.

No common exponent for all original Q is asserted. For every k>=0,
use (R2) to write P^k=sum_(j<d)P^j A_(j,k)(f(P)). It follows that

    E[Q f(P)^(N_Q) P^k]=0 for every k>=0.

Thus, with D=partial_t, EVERY Taylor coefficient vanishes in

    f(D)^(N_Q) Phi_P[Q]=0.                          (R3)

This is stronger than a finite Taylor-prefix check: the free-module
reduction supplies all remaining coefficients.

Let L/K be a finite splitting field of f and let Lambda be its distinct
roots. The elementary constant-coefficient differential equation gives

    ker(f(D)^N) in L[[t]]
       =direct sum_(lambda in Lambda)
          exp(lambda t)L[t]_(degree<N mult_lambda(f)).

For one factor this follows from
(D-lambda)^N h=exp(lambda t)D^N(exp(-lambda t)h);
coprime factors split by the polynomial Bezout identity in D. Consequently
all original P-periods lie in the finite L[t]-module

    E_f=direct sum_(lambda in Lambda)exp(lambda t)L[t].                (R4)

The same module is finite over K[t] by restriction of scalars, because
[L:K] is finite. The actual M_P is a K[t]-submodule of it, so Noetherianity
of K[t] proves that M_P is finite over K[t]. This proves the required
descent without assuming in advance that the selected root is K-rational.

## 2. Finite actual image selects a linear exponential

The retained finite-image theorem applies over K. In brief, the exact
Gaussian source has one-dimensional reduction modulo t. Its surjection
onto a finite torsion-free K[t]-image therefore makes that image free of
rank one. Normalize its generator g by g(0)=1. Derivative-stability gives
g'=a(t)g, a in K[t], hence

    M_P=exp(psi)K[t], psi in tK[t],
    Phi_P[1]=p(t)exp(psi), p in K[t], p(0)=1.        (R5)

There is no assertion p=1. Since the generator belongs to (R4), some
power f(D)^N annihilates exp(psi).

This forces psi to be linear. If v=psi' had positive degree s, then

    B_j=exp(-psi)D^j exp(psi), B_0=1,
    B_(j+1)=B_j'+v B_j

would have degree js and leading term v^j. In a nonzero constant-
coefficient operator of order r, its order-r term then gives the unique
highest degree rs after conjugation and evaluation on one. Such an
operator cannot annihilate exp(psi). Thus psi=lambda t. Evaluating
f(D)^N on this exponential gives f(lambda)^N=0, so

    lambda in K, f(lambda)=0, M_P=exp(lambda t)K[t].                  (R6)

In particular lambda lies in the ORIGINAL coefficient field, even when
f initially required a splitting field.

Equation (R6) implies that exp(-lambda t)Phi_P[Q] is polynomial for
every original Q. Hence P-lambda is safe, proving the forward direction
of (R1) for nonconstant P.

For clarity, the following are equivalent for any lambda in K:

    P-lambda is safe;
    M_P=exp(lambda t)K[t];
    for every Q there is N_Q with
       (D-lambda)^(N_Q)Phi_P[Q]=0.                  (R7)

The last equivalence is the same elementary one-factor differential
equation. For the first-to-second implication, shift by lambda. Its
polynomial actual image is a nonzero derivative-stable ideal of K[t],
and every such ideal is K[t]: repeatedly differentiating a polynomial
of least positive degree produces a nonzero constant. When P is
nonconstant, the retained fixed-original-marker theorem further gives

    V_P=M_P=exp(lambda t)K[t].                      (R8)

The constant-phase exception is recorded below. There cannot be a
single annihilator exponent valid for ALL Q in (R7) when P is nonconstant:
(R8) contains exp(lambda t)t^r for arbitrarily large r. Finiteness over
K[t], rather than generic rank one alone, is essential; a rank-one image
with additional rational-time poles need not satisfy (R7).

## 3. Converse, uniqueness, safe ideals and constants

Suppose S=P-lambda is safe and f(lambda)=0. Then
f(lambda+S)=S h(S) for some h in K[T]. For any original Q choose N_Q
such that E[Q S^k]=0 for all k>=N_Q. Every power

    f(P)^n=S^n h(S)^n

is a finite linear combination of S^k with k>=n. Thus it has zero Q-marked
expectation for n>=N_Q. This proves the converse of (R1), with no
primitivity assumption.

A safe shift is unique. If both P-lambda and P-mu are safe, the two
nonzero polynomials

    p_lambda=exp(-lambda t)Phi_P[1],
    p_mu=exp(-mu t)Phi_P[1]

have ratio exp((mu-lambda)t). A nonzero rational function cannot have
nonzero constant logarithmic derivative: its logarithmic derivative
tends to zero at infinity. Therefore lambda=mu. This is an algebraic
formal-series argument in characteristic zero.

If P is the constant c, then f(P) is safe precisely when f(c)=0.
The unique safe shift is lambda=c. Its images are

    V_P=K exp(ct), M_P=exp(ct)K[t],

so (R8) requires nonconstancy exactly as stated. A nonzero constant f
never yields a safe phase, whereas f identically zero yields the zero
safe phase for EVERY P and supplies no information about a safe shift.

Equivalently, for every fixed P the complete set

    S_P={f in K[T] : f(P) is safe}

is either the zero ideal or the maximal ideal

    S_P=(T-lambda)K[T]                              (R9)

for its unique safe shift lambda. This includes the zero recomposition
and treats constants correctly. It is a classification of safe
univariate recompositions, not a classification of all safe phases.

## 4. Every nonzero zero-constant recomposition of a pure unsafe phase

Let P be pure and unsafe. For EVERY nonzero f in T K[T], the phase
f(P) is again pure and unsafe.

Purity follows directly: each f(P)^n, n>=1, is a finite combination of
positive powers of P. If f(P) were safe, (R1) would give a safe shift
P-lambda. But purity of P gives

    Phi_(P-lambda)[1]=exp(-lambda t).

Safety makes this polynomial, forcing lambda=0. Then P itself would
be safe, a contradiction.

The zero polynomial f is excluded: its recomposition is always safe.
The assertion allows the bad original marker to depend on f. It does
not claim that one prescribed marker exposes infinitely many moments
for every recomposition.

The [utility](gaussian_safe_recomposition.py),
[focused checker](verify_gaussian_safe_recomposition.py) and
[verification receipt](gaussian_safe_recomposition_verification.json)
provide bounded exact controls of the displayed algebra and covered
examples. The unrestricted theorem is the written free-module,
differential-equation and original-image argument above; no finite
moment prefix supplies its all-order scope. Independent replay,
formal verification and external review are not asserted.

## 5. Exact filtered original images at the cubic shear

In this section L may be either the Witt Gaussian functional or the
normalized S2 sphere functional. These assertions concern the specific
shear phase P*=W+ZY-Z^3/2; they do not extend (R1) to arbitrary sphere
phases. Use exponential generating functions for both functionals.
The exact filtered original image at P* is

    V_D={L(A exp(sP*)):deg A<=D}=K[s]_(degree<=3D).                  (R10)

Here and below the degree bound is on the ORIGINAL spatial marker.
For completeness, the inverse Gaussian/sphere-preserving shear sends
Z,Y,W to Z,Y+Z^2,W-ZY-Z^3/2. The inverse image of Z^iY^jW^k has
largest Z-minus-W weight n=i+2j+3k, with unique top term
(-1/2)^k Z^n. Its period therefore has leading coefficient

    (-1/2)^k                              Gaussian;
    (-1/2)^k/(2n+1)!!                     sphere.                  (R11)

For n=3k,3k+1,3k+2 choose respectively W^k,ZW^k,YW^k. These markers
have degree ceil(n/3) and form a triangular family of periods in
1,s,...,s^n. The same weight bound gives the reverse inclusion in
(R10), proving it exactly for both functionals.

Let f(U) be nonzero, f(0)=0, e=ord_0 f>=1, and let a_e!=0 be its
coefficient of U^e. On polynomial periods define the finite linear map

    R_f(s^j)=j! sum_(m>=0) [U^j]f(U)^m t^m/m!.                      (R12)

Expansion of the actual original moments gives

    L(A exp(tf(P*)))=R_f(L(A exp(sP*))).                            (R13)

The factorial j!/m! is required by the two exponential generating
functions. No infinite interchange or analytic convergence is used:
the input period is polynomial, and f^m is a finite polynomial.
Since deg R_f(s^j)<=floor(j/e), while the coefficient of t^m in
R_f(s^(em)) is (em)!a_e^m/m!, triangular elimination gives

    V_D^[f(P*)]=K[t]_(degree<=floor(3D/e)).                         (R14)

In particular f(P*) is pure and safe. The least original marker
degree producing t^m, or any period of exact degree m, is precisely

    D_min(m)=ceil(em/3).                                            (R15)

To construct a lift, use the preceding original marker of weight em;
its recomposed period has exact degree m. Eliminate its lower terms
with the markers of weights 0,e,...,e(m-1). Every marker used has
degree at most ceil(em/3). The lower bound follows directly from
(R14). The zero polynomial f is separate: its phase has filtered
image K and positive-degree target monomials have no lifts.

For the explicit degree-fifteen phase S=P*^2+7P*^5, e=2 and the exact
cutoff is floor(3D/2). The following degree-at-most-two markers give
its complete triangular D=2 image:

| Original marker | Gaussian period | Sphere period |
| --- | --- | --- |
| 1 | 1 | 1 |
| Y | 2t | 2t/15 |
| ZW | 1-6t^2 | 1/3-2t^2/315 |
| W^2 | -4t+30t^3 | -4t/105+2t^3/9009 |

The degree-three marker YW^2 supplies the next triangular period:
-36t^2+420t^4 in the Gaussian case and
-4t^2/1155+4t^4/328185 on the sphere. Thus these five markers span
the complete D=3 image. The further degree-two fixture YW has period
-420t for Gaussian moments and -4t/99 for sphere moments. It detects
the coefficient 7 of f, which does not appear in the four D=2
triangular periods above.

## 6. Sharp finite moment-jet budgets after recomposition

Fix S=f(P*) as above, and let H vanish at the parameter origin and
have original spatial degree at most d. Modulo the (r+1)-st power
of the parameter maximal ideal, at most r perturbation factors occur.
A binomial term with j such factors has marker H^j of degree at most
dj. Equation (R14) kills it whenever n-j>(3d/e)j. Hence all r-jet
pure moments vanish identically beyond

    B_r=floor((1+3d/e)r).                                          (R16)

Together with a certified containment m^nu subset I_N, the finite
intervening checks in the
[moment-jet tail theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#2-a-finite-certificate-for-an-infinite-tail)
therefore extend only through floor((1+3d/e)(nu-1)). That theorem's
local/global hypotheses remain required; a jet calculation alone
does not promote local ideal stabilization to global stabilization.

The MOMENT jet budget (R16) is sharp on the FULL perturbation space
of spatial degree at most d, for every r>=1. Put
n=floor(3dr/e), so n+r=B_r. By (R14) some original marker Q of
degree at most dr has L(QS^n)!=0. That marker space is spanned by
r-fold products B_1...B_r with deg B_i<=d: partition the factors of
each monomial into r groups, allowing constants. Thus some such
product has nonzero marked moment. The r-th derivative in these
directions of L((S+H)^(n+r)) at H=0 is

    (n+r)!/n! L(B_1...B_r S^n)!=0.

Polarization in characteristic zero also supplies a diagonal
perturbation B with L(B^r S^n)!=0. This B may depend on r. Therefore
no smaller uniform moment-jet threshold works on the entire degree-d
perturbation space. As with the earlier cubic threshold, this does
not show that B_r is the least IDEAL-prefix stabilization bound:
lower moments may already imply the last nonzero jet equation.

The degree-fifteen benchmark has e=2. For cubic perturbations its
budget is floor(11r/2), and nu=7 requires checks only through 33.
A single explicit direction realizes sharpness for every r:

    H=epsilon H0, H0=W^3+YW^2.

Set B=floor(11r/2), n=B-r=floor(9r/2), j=2n. Then

    [epsilon^r] L((S+epsilon H0)^B)
       =binom(B,r)j!c_r                         Gaussian;
       =binom(B,r)j!c_r/(2j+1)!!                sphere,             (R17)

where c_r=(-1/2)^(3r) for even r and
c_r=r(-1/2)^(3r-1) for odd r. Both are nonzero. Indeed f(U)^n has
leading term U^(2n), and every further term raises its U-degree by
at least three, exceeding the marker bound 9r. In H0^r the unique
largest circular weight 9r comes from W^(3r); the unique next weight
9r-1 comes from rYW^(3r-1). Their inverse-shear top coefficients are
exactly c_r. This proves (R17) with the normalization (R11).

The utility computes (R12)-(R16) with rational coefficients. The
focused checker also compares small original Wick moments, explicit
Gaussian/sphere fixtures and finite perturbation jets. The unrestricted
filtered and sharpness statements follow from the weight, triangular
and polarization arguments above, not a finite rank or moment test.

## 7. Uniform bad markers for specified Gamma and rational branches

For a general known marked period there is a finite candidate list.
If Phi_P[Q] is not annihilated by any power of f(D), with d=deg f,
then at least one of

    Q,QP,...,QP^(d-1)

is an unsafe marker for f(P). Otherwise safety for this finite list
would give (R3), a contradiction. In particular a genuine finite pole
of the original period meets this differential criterion for every
nonconstant f, because constant-coefficient differential equations have
only exponential-polynomial solutions. This list need not identify one
marker that works for every f. The following stronger conclusion does
identify such a marker for the specified periods below.

For an integer r>=1 define the polynomial Gamma functional by

    G_r[F]=1/Gamma(r) integral_0^infinity F(s)s^(r-1)exp(-s)ds,
    G_r[s^k]=(r)_k.

**Gamma lemma.** If F(s)=L s^q+L beta s^(q-1)+... is any nonconstant
complex polynomial, L!=0, then

    G_r[F^n]=L^n Gamma(qn+r)/Gamma(r)
                       * (exp(beta/q)+o(1)).                       (R18)

In particular the moments are nonzero for every sufficiently large n.
Arbitrary complex coefficients are permitted.

To prove this, let S_n have the positive Gamma distribution of shape
qn+r and rate one. After dividing by L^n Gamma(qn+r)/Gamma(r),
the integral in (R18) is the expectation of

    (F(S_n)/(L S_n^q))^n.

On S_n>=qn/2, the logarithm near one gives uniformly

    n log(F(S_n)/(L S_n^q))
         =beta*n/S_n+O(n/S_n^2).

The exact mean and variance are qn+r, so S_n/n tends in probability
to q. Choose C>0 with |F(s)|<=|L|(s+C)^q for all real s>=0;
comparison of coefficients with the binomial expansion supplies such C.
On this region the absolute integrand is bounded by
(1+C/S_n)^(qn)<=exp(2C). Its expectation therefore tends to
exp(beta/q).

The complementary normalized absolute integral is at most

    1/Gamma(qn+r) integral_0^(qn/2)
        exp(-s)s^(r-1)(s+C)^(qn)ds
      <=exp(C) Pr(S_n<=qn/2+C).

The inequality follows by translating s to s+C and using r>=1.
Markov's inequality applied to exp(-S_n), whose expectation is
2^(-qn-r), bounds the last expression by
exp(qn/2+2C)2^(-qn-r), which tends to zero exponentially.
This proves (R18), including its nonzero limit, without assuming
positivity of F or of its coefficients.

Now consider the entire pure rational branch in the
[affine-isotropic classification](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#statement),
including its arbitrary polynomial ell. Its
[original circular elimination](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-circular-elimination)
and [marked-period calculation](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#unsafe-branch-and-admissible-fibre-directions)
give

    Phi_P[Z]=ct/(1-dt), c,d!=0,
    E[Z]=0, E[ZP^k]=c d^(k-1)k! for k>=1.

Write f(T)=a_q T^q+a_(q-1)T^(q-1)+..., a_q!=0, and suppose
f(0)=0. For every n>=1 the exact original marked moment is

    E[Z f(P)^n]=(c/d) G_1[f(ds)^n].                              (R19)

Indeed f(T)^n has zero constant term, so the displayed original
P-moments determine both sides coefficientwise. Applying (R18) gives

    E[Z f(P)^n]
      =(c/d)(a_q d^q)^n(qn)!
            * (exp(a_(q-1)/(q d a_q))+o(1)).                     (R20)

This is nonzero for every sufficiently large n. Thus the SAME original
marker Z proves unsafety for EVERY nonzero zero-constant f in this
whole rational branch. Long's quartic is its c=d=1,ell=0 case.
The nonvanishing threshold may depend on f.
For a general characteristic-zero coefficient field, embed the finitely
generated field of the displayed phase and f coefficients into C.
Injectivity transfers the resulting nonvanishing back to that field.
When f(0)!=0, purity of P already gives the nonpolynomial unit
Phi_(f(P))[1]=exp(f(0)t).

A second explicit example is the
[pure angular cubic boundary](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#a-pure-unsafe-boundary-with-its-full-original-image).
With S=ZW, its original phase and marker are

    P=-2S+S(Y1+iY2)+Y1^2+Y2^2-2Y1+2iY2,
    B=(3-S)/2-Y1, Phi_P[B]=(1-2t)^(-1).

Therefore E[B f(P)^n]=G_1[f(2s)^n], and (R18) makes B a uniform
bad marker for every nonconstant f. For f(0)=0 these recompositions
are pure and unsafe. These uniform conclusions concern the stated
branches; they do not strengthen Section 4 to a prescribed uniform
marker for every pure unsafe phase.
