# Polynomial global root bounds from moving Grace discs

This note proves a global localization theorem for every root of the actual
ordinary factorial moment polynomial. It does not prove consecutive
coprimality, a uniform local saddle theorem, or a polynomial Strong Factorial
threshold.

**Theorem PolynomialGlobalHermiteStrip.** Let d>=2, N>=1 be integers and
let a>=100*d^2 be real. Define

    A_(N;a,d)(x)=sum_(k=0)^N binom(N,k)*(aN+1)_(dk)*x^(N-k).

For every complex root x, put

    w=x/N^d, delta=d/sqrt(a), Z=(w/a^d+1)/delta.

Then

    |Re Z|<7*d^2/sqrt(a),       |Im Z|^2<17*d.             (P1)

These inequalities hold for every stated finite order, not only in an
asymptotic zero distribution. When a is an integer these are precisely
the actual reversed ordinary factorial moments of X^a(A+B*X^d), with
coefficient ratio x=A/B. No fixed-shape substitute is made.

**Corollary PolynomialThinStrip.** In particular, a>=10^6*d^4 gives |Re Z|<7/1000. At a>=10^18*d^4 one
has |Re Z|<7*10^(-9) and delta<=10^(-9)/d. The height bound still grows
as sqrt(d); a uniform outside-one-packet argument is a separate obligation.

## 1. Exact classical inputs and notation

For n>=1, h=n-1, gamma>0 let

    G_(n,gamma)(z)=sum_(k=0)^n binom(n,k)*(gamma)_k*z^(n-k).

The two elementary inputs are proved in the retained account
GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md, Sections1--2:

1. Every root beta of G lies in |beta|<=gamma+3h and satisfies
   Re beta<=-gamma+h. These follow from its explicit complex-symmetric
   Laguerre Jacobi matrix and numerical-range bounds, including h0.
2. For two exact degree-n polynomials, if every root of the first lies
   in a closed disc D(c,R), R>0, and the second has nonzero constant,
   every root of their unsigned binomial composition has the form
   -alpha*beta with alpha in that disc and beta a root of the second.

The second statement is classical Grace--Szego. We additionally inspected
the actual published statement in Martin Lamprecht, *A continuous link
between the disk and half-plane cases of Grace's theorem*, C. R. Acad.
Sci. Paris I353(2015),11--15, Theorem1.3 on printed13, DOI
10.1016/j.crma.2014.10.017. It allows the open or closed disc used here,
requires the two exact degrees n, and gives exactly the MINUS-product
representation. The retained self-contained polar proof supplies this
mathematical input; the new primary inspection confirms its classical
scope. No real-rootedness, positive Laguerre measure, or simultaneous
choice of beta across different discs is required.

For positive gamma_1,...,gamma_d let C_r be the binomial composition of
the first r factors, and M_r=product_(j=1)^r gamma_j. Write root bounds

    |z|<=M_r*(1+e_r),       Re z<=M_r*(-1+u_r),            (P2)

where e_r,u_r are nonnegative. Suppose h/gamma_j<=epsilon and
d*epsilon<=1/100. We prove bounds independent of the individual gammas:

    e_r<=4*r*epsilon,       u_r<=6*r^2*epsilon, r<=d.      (P3)

## 2. The moving-disc recurrence

The origin-disc composition gives

    e_(r+1)=epsilon_(r+1)+(1+epsilon_(r+1))*e_r,
    epsilon_(r+1)=3h/gamma_(r+1)<=3epsilon.

Thus 1+e_r<=product_(j=1)^r(1+3h/gamma_j)<=(1+3epsilon)^r.
Using binom(r,k)<=r^k and r*3epsilon<1 yields

    e_r<=3r*epsilon/(1-3r*epsilon)
        <=(300/97)*r*epsilon<=4r*epsilon.                 (P4)

For the real-part bound, center the NEW containing disc at
-(r+1)*M_r, rather than at the fixed multiple -2*M_r. From(P2),

    |z+(r+1)M_r|^2/M_r^2
      <=r^2+2e_r+e_r^2+2(r+1)u_r.

The square-root inequality sqrt(r^2+t)<=r+t/(2r), t>=0, gives an
admissible disc radius M_r times

    r+(e_r+e_r^2/2+(r+1)u_r)/r.

Apply closed-disc Grace once with THIS disc and the next factor. Its
root beta obeys Re beta<=-gamma+h, |beta|<=gamma+3h. Dividing the
resulting real-part bound by M_r*gamma and using h/gamma<=epsilon gives

    u_(r+1)<=(4r+1)*epsilon
       +(1+3epsilon)*(e_r/r+e_r^2/(2r)+(1+1/r)u_r).       (P5)

The disc chosen for this application is unrelated to the one used for
the origin bound. Different applications may select different beta.
The amplification of u_r is now1+1/r instead of2.

By(P4) and r*epsilon<=1/100,

    e_r/r+e_r^2/(2r)<=4epsilon+8r*epsilon^2<5epsilon.

The base bound is u_1=h/gamma_1<=epsilon. Assuming u_r<=6r^2epsilon,
(P5) is at most

    epsilon*((4r+1)+(1+3epsilon)*(5+6r(r+1))).

Its difference below6(r+1)^2epsilon, after division by epsilon, is

    2r-3epsilon*(5+6r(r+1)).

At every induction step r+1<=d, so epsilon<=1/[100(r+1)]. Multiplying
this last lower bound by100(r+1) gives

    182r(r+1)-15>0.                                     (P6)

At r=1+s, s>=0, the exact positive polynomial is
349+546s+182s^2. This proves(P3). If epsilon0 all factors have h0,
their compositions are (z+M_r)^n with n1, and(P3) holds directly;
there is no division by zero in the claimed bounds.

## 3. Actual factorial parameters and physical rescaling

Take n=N, h=N-1 and gamma_j=(aN+j)/d. Then h/gamma_j<d/a=:epsilon
and d*epsilon=d^2/a=:eta<=1/100. The elementary multiplication identity

    (aN+1)_(dk)=d^(dk)*product_(j=1)^d ((aN+j)/d)_k

gives EXACTLY

    A_(N;a,d)(x)=d^(dN)*C_d(x/d^d).                       (P7)

For y=(x/d^d)/M_d, (P3) therefore gives

    |y|<=1+e, Re y<=-1+B,
    e=4eta, B=6d*eta.

Moreover

    kappa=d^d*M_d/(a^d*N^d)
         =product_(j=1)^d(1+j/(aN)),
    1<=kappa<100/99,       kappa-1<(100/99)*eta.           (P8)

To prove(P8), sum j/(aN)<=eta and use
product(1+x_j)<=1/(1-sum x_j) for nonnegative x_j with sum<1.
For d>=2 the sum is strictly below eta. Thus w/a^d=kappa*y.
The lower real error is at most

    (kappa-1)+kappa*e<(500/99)*eta<6eta,

and the upper real error is at most

    kappa*B<(600/99)*d*eta<7d*eta.

Since delta=d/sqrt(a), these imply the first inequality in(P1).

The height estimate must NOT silently assume B<1. If 0<=B<=1, then
Re y<=-1+B<=0, so

    |Im y|^2<=(1+e)^2-(1-B)^2<=(e+B)*(2+e).

The difference in the last comparison is B*(e+B)>=0. If B>=1, the
origin bound alone gives

    |Im y|^2<=(1+e)^2<=(e+B)*(2+e),

since B*(2+e)-1>=1+e. Thus the SAME height bound holds in both cases.
For d>=2, e+B=(4+6d)*eta<=8d*eta. It follows that

    |Im Z|^2<=kappa^2/delta^2*(e+B)*(2+e)
       <(100/99)^2*8*(51/25)*d<17d,                      (P9)

because eta=delta^2. This completes the all-root theorem.

## 4. Classical comparison and remaining gap

The older retained fixed-gap bound used a center -2M at EVERY step;
its u recurrence doubles the previous real-part error and its displayed
constants grow exponentially with d. The moving center -(r+1)M removes
that loss with the very same exact Grace theorem and factor hypotheses.

The actual primary Bihun--Calogero arXiv1507.07808v1, Section2,
Proposition2.1 was also inspected: it supplies algebraic identities for
zeros of generic hypergeometric polynomials, not the finite-N strip
claimed here. Those identities are not imported as a localization
theorem. Neither classical comparison justifies a novelty claim.

The new global bound supplies a polynomially small REAL strip once
a is a constant times d^4. The remaining uniform analytic obligations
are a complete outside-one-packet bound for the growing imaginary
range |Im Z|<sqrt(17d), and the uniform-in-d local/endpoint and
actual-adjacent-order comparison with explicit constants. No SFC tail
theorem follows from(P1) alone.

## 5. The exact common-root equivalence on the proved region

**Corollary ActualPolynomialStripNoncollisionEquivalence.** Let integers
d>=2, a>=100*d^2 and N>=2 be fixed, and put n=N-1. Define the SAME
upper-order affine argument

 S_N(Z)=N^d*a^d*(-1+(d/sqrt(a))*Z),
 U(Z)=A_(N;a,d)(S_N(Z)),
 L(Z)=A_(n;a,d)(S_N(Z)).

The lower polynomial keeps its own degree n and actual shape a*n+1.
Let R_ad be the open rectangle

 |Re Z|<7*d^2/sqrt(a), |Im Z|^2<17*d.

Then gcd_Q(A_(N;a,d),A_(n;a,d))=1 if and only if U and L have
no common complex zero in R_ad.

The forward implication is immediate. For the reverse, every upper
root lies in R_ad after the invertible affine change S_N, by the
all-root theorem. Any common root must therefore lie there. No common
complex root is equivalent to rational coprimality. The latter is also
equivalent to common-root exclusion over every characteristic-zero
extension field, by rational Bezout transfer. Failure of rational
coprimality does not require a root in a specified nonclosed field.

Quantifying either property over ALL stated integer a,d,N gives
equivalent universal assertions. Neither is proved by this root
localization theorem. Uniform all-gap analytic control, including the
growing imaginary range and an actual adjacent-order comparison,
remains a separate open obligation.

## Exact reproduction

    python3 research/factorial/image-rigidity/verify_iterated_grace_polynomial_root_strip.py \
      --repo-root . --output /tmp/iterated_grace_root_strip_review.json

The exact controls include common-affine substitution in actual
upper/lower moment polynomials. Finite examples support the algebra;
the reverse implication is the written all-root argument above.

Updated 20 November 2026.
