# Eventual positivity obstructs the polynomial radial normalizer for a Gaussian

Let gamma be the standard N(0,1) probability on R, let X denote its coordinate, and let Lambda be Gaussian expectation followed by Laurent constant term in t on R[X,t,t^-1]. Fix ANY nonzero real polynomial f. Put T=1+t and

    p=(T/t)*f(TX)^2.

An element g is called **pure** for Lambda if Lambda(g^M)=0 for EVERY positive integer M.

**Theorem 1.** Lambda(p^k)>0 for EVERY sufficiently large integer k. Consequently, for EVERY positive integer a, the moments of both p^a and F_a=p^a(1+p) are strictly positive at every sufficiently large order. Neither family is pure. If deg(f)<=1, Lambda(p^k)>0 at EVERY positive k.

For rational f the exact native moments are rational; their eventual nonvanishing persists after EVERY characteristic-zero scalar extension of the rational Gaussian-moment functional. The theorem concerns this real-polynomial normalizer ansatz. It is not a classification of Gaussian Mathieu kernels or a conclusion about arbitrary complex amplifiers.

## 1. Exact Hermite identity

Use probabilists' Hermite polynomials

    He_k(x)=(-1)^k*exp(x^2/2)*(d/dx)^k exp(-x^2/2),
    exp(zx-z^2/2)=sum_(k>=0) He_k(x) z^k/k!.

For any polynomial h, monomial comparison gives

    [t^k](1+t)^k E h((1+t)X)
        = (1/k!)*E[(d/dX)^k (X^k h(X))].                            (1)

Indeed h(X)=X^j gives binom(k+j,k)E X^j on both sides. Gaussian integration by parts k times, with all boundary terms zero for polynomials times the Gaussian density, gives

    E[(d/dX)^k (X^k h(X))]=E[X^k He_k(X)h(X)].                       (2)

Since p^k=t^-k T^k f(TX)^(2k), applying(1)-(2) to h=f^(2k) proves the EXACT identity

    Lambda(p^k)=(1/k!)*E[X^k He_k(X) f(X)^(2k)].                     (3)

This proof is finite polynomial algebra and Gaussian IBP at each k. It needs no exchange of an infinite parameter series with a changing polynomial degree. In particular Lambda(p)=E[X^2 f(X)^2]>0 for nonzero real f, although that first moment alone would not exclude purity of all higher power families.

## 2. The classical diagonal Mehler bound, with conventions proved

If Z is a standard real Gaussian, the generating function gives

    He_k(x)=E[(x+iZ)^k].                                           (4)

For independent standard Gaussians Z,W and |r|<1, absolute domination by exp(|r|(Z^2+W^2)/2) times a fixed constant justifies summing the exponential series inside expectation. Thus

    sum_(k>=0) He_k(x)He_k(y) r^k/k!
        =E exp(r(x+iZ)(y+iW))
        =(1-r^2)^(-1/2)
           *exp([2rxy-r^2(x^2+y^2)]/[2(1-r^2)]).                   (5)

The last equality is the elementary two-Gaussian integral: integrate W first using its complex linear exponential moment, then integrate Z with quadratic coefficient r^2/2<1/2. This also proves the identity directly. At x=y it becomes

    sum_(j>=0) He_j(x)^2 r^j/j!
        =(1-r^2)^(-1/2)*exp(r x^2/(1+r)),       0<r<1.             (6)

Every summand in(6) is nonnegative. For k>=2 choose r=1-1/k. The elementary inequality (1-1/k)^k>=1/4 and 1-r^2>=1/k give

    |He_k(x)| <= 2*k^(1/4)*sqrt(k!)*exp(x^2/4),                    (7)

for EVERY real x. Here r/(2(1+r))<=1/4. The inequality (1-1/k)^k>=1/4 follows, for example, from monotonicity of u log(1-1/u) for u>1, whose derivative is log(1-1/u)+1/(u-1)>0.

Formula(5) is the classical Mehler identity. It agrees with NIST DLMF18.18.28 after the convention change He_k(x)=2^(-k/2)H_k(x/sqrt(2)); the needed identity is proved above rather than imported as an unexamined estimate.

## 3. Finite Jacobi matrix, root bound and second root moment

The generating function implies

    He_0=1, He_1=x,
    He_(j+1)=x He_j-j He_(j-1).                                    (8)

Let J_k be the real symmetric k-by-k tridiagonal matrix with zero diagonal and off-diagonals sqrt(1),...,sqrt(k-1). Expansion of its characteristic determinant along the last row gives exactly recurrence(8). Therefore its eigenvalues r_1,...,r_k are the roots of He_k, including multiplicity.

Every eigenvalue is real. The maximum absolute row sum bounds its absolute value by2sqrt(k). Conjugation by diag(1,-1,1,-1,...) sends J_k to -J_k, so the roots occur in +/- pairs, with any unpaired root zero. Finally

    sum_j r_j^2=tr(J_k^2)=2 sum_(j=1)^(k-1) j=k(k-1).              (9)

In particular x^k He_k(x)>=0 outside |x|<=2sqrt(k). No limiting zero distribution is required.

For y>2 write a_j=r_j^2/k for one root from each +/- pair. Then0<=a_j<=4 and sum_pairs a_j=(k-1)/2. Since a->log(1-a/y^2) is concave on[0,4], it lies above its endpoint chord:

    log(1-a/y^2)>=(a/4)log(1-4/y^2).

The paired product consequently gives

    He_k(sqrt(k)*y)
       >= k^(k/2)*y^k*(1-4/y^2)^((k-1)/8),       y>2.              (10)

The exponent is exactly(k-1)/8 because(9) is halved when one root from each pair is used. Odd k is included through the unpaired zero root.

## 4. Uniform rescaling of an arbitrary fixed real polynomial

Suppose delta=deg(f)>=2 and c is its nonzero leading coefficient. For real y put

    g_k(y)=f(sqrt(k)*y)/(c*k^(delta/2)).

For every fixed compact interval, g_k(y)->y^delta uniformly. This is immediate by dividing each lower-degree term by k^(delta/2). Hence g_k(y)^2->y^(2delta) uniformly there as well. The statements below use suprema and infima of these uniformly convergent continuous functions, so lower coefficients may affect subexponential factors but cannot change the exponential comparison.

Define

    H_k(y)=|y|*g_k(y)^2*exp(-y^2/4),    |y|<=2.

Its maximum converges to the maximum of |y|^(2delta+1)exp(-y^2/4). For0<y<=2 the logarithmic derivative of this last function is

    (2delta+1)/y-y/2>0.

Thus

    max_(|y|<=2) H_k(y) -> 2^(2delta+1)*exp(-1).                    (11)

Set the negative-region exponential base

    B_delta=2^(2delta+1)*exp(-3/2).                                (12)

## 5. All potentially negative mass has base at most B_delta

Let

    N_k=(1/sqrt(2pi))*integral_(|x|<=2sqrt(k))
                  |x|^k*|He_k(x)|*f(x)^(2k)*exp(-x^2/2)dx.

This upper-bounds the absolute value of EVERY potentially negative contribution to the expectation in(3). By(7), the substitution x=sqrt(k)y and the definition of H_k,

    N_k <= C*k^(1/4)*sqrt(k!)*k^((k+1)/2)
              *|c|^(2k)*k^(delta*k)
              *[max_(|y|<=2)H_k(y)]^k.                            (13)

The length4 of the rescaled interval and the Gaussian normalization are absorbed into a fixed C. An elementary factorial bound is enough: since log is increasing,

    log(k!)<=integral_1^k log(u)du+log(k)
             =k log(k)-k+1+log(k).

Therefore sqrt(k!)<=sqrt(e)*k^(1/2)*k^(k/2)*exp(-k/2). Equations(11)-(13) prove that, for ANY fixed b>B_delta and all sufficiently large k,

    N_k <= C*k^(5/4)*|c|^(2k)*k^((delta+1)k)*b^k.                  (14)

All prefactors are only polynomial in k. A sharp Stirling theorem is unnecessary.

## 6. A positive interval has a strictly larger exponential base

For y>2 define

    A_delta(y)=y^(2delta+2)*exp(-y^2/2)*(1-4/y^2)^(1/8).

At y0=sqrt(6),

    A_delta(y0)=6^(delta+1)*exp(-3)*3^(-1/8),
    log[A_delta(y0)/B_delta]
        =delta log(3/2)+(7/8)log(3)-3/2.                            (15)

This is positive for every delta>=2. For an entirely elementary strict bound, log(3/2)>2/5 follows by differentiating log(u)-2(u-1)/(u+1), and log(3)>1 follows from e<3. Thus the right side of(15) is greater than

    4/5+7/8-3/2=7/40>0.

Choose fixed numbers B_delta<b<d<A_delta(y0). By continuity choose a closed interval I of positive length, lying strictly inside(2,infinity) and containing y0, such that min_I A_delta>d. Uniform convergence from Section4 then implies, for all sufficiently large k,

    min_(y in I) [y^2*g_k(y)^2*exp(-y^2/2)
                         *(1-4/y^2)^(1/8)] >= d.                  (16)

Let P_k be the contribution to the expectation in(3) from x=sqrt(k)y, y in I. All this contribution is positive. In(10) the residual factor(1-4/y^2)^(-1/8) is at least1. Using(10),(16) gives

    P_k >= C_I*sqrt(k)*|c|^(2k)*k^((delta+1)k)*d^k,                 (17)

where C_I=length(I)/sqrt(2pi)>0. Comparing(17) with(14),

    P_k/N_k >= C'*k^(-3/4)*(d/b)^k -> infinity                    (18)

whenever N_k>0; if N_k=0 the positivity conclusion is immediate. Since x^kHe_k(x)f(x)^(2k) is nonnegative everywhere outside the inner interval, the complete expectation in(3) is at least P_k-N_k>0 for every sufficiently large k. Division by k! preserves the sign. This proves the theorem for delta>=2.

The proof does not assume that f has no real zeros. Uniform rescaling makes g_k nonzero on I for large k because its limit y^delta is bounded away from zero there. Every possible negative contribution, including arbitrary behavior near all fixed real roots, was already bounded by N_k.

## 7. Constants, linear polynomials and the two power families

If f(x)=alpha*x+beta is nonzero, finite expansion directly gives

    Lambda(p^k)=sum_(j=0)^k binom(2k,2j)
          *alpha^(2j)*beta^(2k-2j)*binom(k+2j,k)*E[X^(2j)].         (19)

Every term is nonnegative and at least one is positive. This covers every nonzero linear or constant polynomial, including beta=0 and alpha=0, and proves positivity at EVERY positive k in these degrees.

Now fix any positive integer a. Once aM is past the positive tail just proved,

    Lambda((p^a)^M)=Lambda(p^(aM))>0.

Likewise the literal positive-coefficient expansion

    Lambda((p^a(1+p))^M)
         =sum_(h=0)^M binom(M,h)*Lambda(p^(aM+h))                  (20)

is strictly positive for all sufficiently large M, since EVERY exponent in it then belongs to the positive tail. Thus neither p^a nor p^a(1+p) is pure. This is an eventual sign theorem for the displayed Gaussian recipe, not a conclusion from checking finitely many powers.

## 8. Rational scalar extensions and exact boundaries

If f belongs to Q[X], all exact moments in(3) or in the original constant-term expression are rational: odd Gaussian monomial moments are0, while E[X^(2j)]=(2j-1)!! is an integer. Define the Gaussian-moment functional over any characteristic-zero field K by K-linearly extending these rational monomial values, and combine it with Laurent constant term. The eventual positive rational values proved above remain nonzero under the unique injection Q->K. Thus the same rational seed/modifier power families are not pure after EVERY such scalar extension. Positivity itself is asserted only in the real realization.

The real-coefficient hypothesis is essential to this proof because f(x)^(2k) is used as a nonnegative weight. No result about arbitrary complex f, arbitrary Laurent phases, a general Gaussian Mathieu conjecture or a Jacobian-conjecture resolution is asserted. The conclusion is the precise obstruction to obtaining a pure seed or modified power family from this particular polynomial radial normalizer over a Gaussian.

Classical identity scope: NIST DLMF [18.18.28, Mehler formula](https://dlmf.nist.gov/18.18.E28), [18.7.11, Hermite normalization](https://dlmf.nist.gov/18.7.E11), and [Table18.9.1, Hermite recurrence](https://dlmf.nist.gov/18.9.T1) give the classical identities corresponding to(5),(8) after the displayed normalization. Mehler, the recurrence, finite root facts and every bound needed here have also been derived in this proof. No Plancherel-Rotach asymptotic or semicircle-limit theorem is imported.
