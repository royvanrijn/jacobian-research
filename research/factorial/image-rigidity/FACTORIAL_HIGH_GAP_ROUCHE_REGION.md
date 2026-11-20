# A uniform high-gap region for Gamma moment simplicity

The following sufficient simplicity criterion holds on an infinite parameter
region. A complementary uniform saddle theorem is linked below.

For integers N>=2,d>=1 and a real gamma>0 write

    A_N(z;gamma)=sum_(k=0)^N binom(N,k)(gamma)_(dk) z^(N-k),
    alpha=d^2*N/[4(gamma+d*N)].

**Proposition.** If N*exp(-alpha)<=1/32, all roots of A_N are simple.
Writing rho=((gamma)_(dN))^(1/N), each root lies in a different disk

    |z/rho-zeta|<1/(4N),   zeta^N=-1.                  (H1)

The result is uniform in all three parameters. It concerns an infinite
parameter region rather than a finite scan.

## Strong convexity and two-term dominance

Set g(x)=log Gamma(gamma+dx)-log Gamma(gamma) for 0<=x<=N.
The actual primary formula DLMF5.15.1 is
psi'(u)=sum_(j>=0)(u+j)^(-2) for u>0. Comparing the decreasing
summands with their integral gives psi'(u)>1/u, hence

    g''(x)>=d^2/(gamma+dN)=m_0.

Apply ordinary convexity to g(x)-m_0*x^2/2, using g(0)=0. For integer k,

    (k/N)g(N)-g(k)>=m_0*k*(N-k)/2
                         >=alpha*min(k,N-k).            (H2)

Normalize B_N(w)=rho^(-N) A_N(rho w;gamma). Its first and last terms
are w^N and1. If S is the sum of the N-1 other positive coefficients,
(H2) and binomial symmetry give

    S<=2 sum_(j>=1) binom(N,j)exp(-alpha*j)
      <=2(exp(N exp(-alpha))-1)
      <=2/31.                                          (H3)

The middle sum can harmlessly include every j=1,...,N. The last inequality
uses exp(x)<=1/(1-x) for 0<=x<1 and the hypothesis x<=1/32.

## Separate disks and simplicity

On |w-zeta|=1/(4N), put w=zeta(1+u). Then |u|=1/(4N), and

    |w^N+1|=|(1+u)^N-1|
      >=1/4-(exp(1/4)-1-1/4)
      >=1/6,

since exp(1/4)<4/3. By(H3) the remaining terms have modulus at most

    S(1+1/(4N))^N <=(2/31)exp(1/4)<8/93<1/6.

Rouche's theorem therefore gives exactly one root, counted with
multiplicity, in each disk. The disks are disjoint: adjacent Nth roots
of -1 have separation 2sin(pi/N)>=4/N, while the sum of their disk
radii is1/(2N). This proves(H1) and simplicity.

## Consecutive factorial consequence

The retained GAP_SIX_RESULTANT_STRICT_POSITIVITY.md formula(1), whose
fixed-shape algebraic derivation is valid for EVERY integer gap d>=1,
identifies the lower-one consecutive resultant at terminal N and Gamma
shape beta with the fixed-shape squarefreeness test for A_N, where

    gamma=beta+N-1.

Thus the proposition proves this consecutive coprimality throughout

    N exp(-d^2*N/[4(beta+N-1+d*N)] ) <=1/32.             (H4)

For the ordinary factorial shape beta=1, gamma=N. A simple sufficient
condition is

    d>=1+4log(32N),                                    (H5)

because d^2/[4(d+1)]>(d-1)/4. Therefore every ordinary binomial
X(A+B X^d), with A,B complex and not both zero, has a nonzero moment
in the pair N-1,N whenever(H5) holds. Monomial degeneracies are immediate.
This includes arbitrarily large, non-prime-power gaps with no arithmetic
restriction on N in the stated region.

## Complementary saddle region and remaining scope

The Rouche argument above covers the explicit high-gap wedge. The
[uniform two-saddle theorem](FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md)
proves simplicity for d>=20 and N>=10^12(d+1)^10, and its
[gap-ten extension](FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md) proves that
same onset for EVERY d>=10. These regions overlap and give ordinary
lower-one consecutive nonvanishing at EVERY terminal for every gapd>=512.

The [complete ordinary lower-one theorem](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md)
now proves consecutive nonvanishing at EVERY gapd>=1 and EVERY terminalN>=2.
It sharpens the uniform tail to N>=max(5000,128(d+1)), overlaps the older
sparse theorem atd44, and completes all171812 remaining full-degree finite
modular tests. The enormous predecessor box above was not enumerated
directly; its validity follows from the complete theorem. The Rouche
argument alone retains its stated broader-Gamma wedge. Lower exponents
greater than one and unrestricted Strong Factorial remain separate.

Primary source actually inspected: https://dlmf.nist.gov/5.15.E1,
DLMF5.15.1, the trigamma series for positive real arguments. The integral
comparison, strong-convexity estimate, normalized coefficient bounds and
Rouche argument are given in full above. The exact resultant bridge is
retained in research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md.
