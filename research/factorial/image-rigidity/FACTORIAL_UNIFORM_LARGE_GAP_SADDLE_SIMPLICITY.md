# Uniform large-gap saddle simplicity and all-order lower-one factorial pairs

20 November 2026. Internal written proof, with exact arithmetic support.
All source equations are finite polynomials or convergent contour integrals.
No all-gap or all-binomial Strong Factorial theorem is asserted.

**Theorem UniformLargeGapSaddleSimplicity.** Let d>=20 and N>=10^12(d+1)^10
be integers. Then

    A_N(z;N)=sum_(k=0)^N binom(N,k)(N)_(dk) z^(N-k)

has only simple complex roots.

**Corollary AllSufficientlyLargeLowerOneGaps.** For every integer d>=512,
every nonzero ordinary factorial binomial X(A+B X^d) has a nonzero moment
in every consecutive pair of positive moments. The coefficients A,B may
be complex, or lie in any characteristic-zero field. The same conclusion
holds at every d>=20 and terminal N>=10^12(d+1)^10.

The proof treats ordinary shape beta=1. It gives no all-positive-Gamma-
shape theorem, no conclusion for general lower exponent greater than one,
and no completion of the remaining small fixed gaps or finite index boxes.
The exact resultant reduction and the complementary high-gap Rouche region
are retained mathematical inputs. The contour estimates below are uniform
in d, with every part of the deformed paths included.

## 1. Normalization, domains and uniform numerical bounds

Put D=d+1, H=D^D exp(-D), W=eH=D^D exp(-d), c=D/2. The Gamma integral gives

    A_N(N^d w;N)=N^(DN)/(N-1)! * I_N(w),
    I_N(w)=int_0^infinity s^(N-1)(w+s^d)^N exp(-Ns) ds.

This identity follows by finite binomial expansion and integrating each
integer Gamma moment; it needs no generating-function convergence. Split
I_N=K_N+J_N at c. Both pieces are polynomials in w. We use the domains

    J: |w|<=3W,            K: |w|>=W/3.                 (L1)

We record numerical inequalities in a form that proves their uniformity.
For d>=20, D=d+1>=21:

    8/3<e<87/32<11/4,   W>D(84/11)^d>D*7^d>400D^2,
    40D^2<2^d, 14D^3<2^d, 160D^2<2^d,
    700D<2^d, 800D<5^d, 8D^2<5^d.                    (L2)

Every polynomial-times-exponential comparison in the second line is checked
at d=20. It continues to hold because D^j/2^d decreases for j<=3:
its successive ratio is at most (22/21)^3/2<1. The corresponding ratio
for D^j/5^d, j<=2, is at most (22/21)^2/5<1. The assertion
7^d>400D is proved at d=20 and its ratio grows by 7D/(D+1)>1.
The bounds for e follow from the factorial series; for the upper bound
sum through degree6 plus the geometric tail1/4410 is less than87/32.
The remaining fixed inequalities used below are

    exp(3/5)<2, exp(9/32)<4/3, exp(3/2)>40/9,
    2500*11^21<16^21,
    exp(1+1/100+1/1000)<25/9.                          (L3)

The first three have rational Taylor lower/upper certificates. The last
follows from e<87/32 and exp(11/1000)<=1000/989. These are exact numerical
inequalities, not floating-point evaluations. The accompanying checker
also verifies each base comparison and its stated monotonicity margin.

## 2. The two stationary points

For J use

    Phi(s,w)=D log s+log(1+w/s^d)-s;

for K_N/w^N use

    Psi(s,w)=log s+log(1+s^d/w)-s.

The stationarity equation for both is

    w(s-1)=s^d(D-s).                                   (L4)

For |w|<=3W, Rouche on |s-D|=1/10 gives exactly one root s_J there.
Indeed, rewriting (L4) as s-D=-w(s-1)/s^d, the right side is at most

    3D^2 exp(-d)*(D/(D-1/10))^d <4D^2*2^(-d)<1/10.

We used (D/(D-1/10))^d<exp(1/10)<10/9. The root is simple, hence
analytic throughout the disk domain in (L1). For |w|>=W/3, Rouche on
|s-1|=1/(100D) gives exactly one root s_K. The other side of s-1=
s^d(D-s)/w is less than4D/W<1/(100D). This root is analytic in1/w,
including1/w=0. The bounds hold on neighborhoods of the closed domains.

At a stationary point minus the second phase derivative is

    a(s)=(s^2-Ds+D)/s^2.

The identities

    a(s)-1/D=(s-D)((D-1)s-D)/(D s^2),
    a(s)-1=-D(s-1)/s^2

give respectively

    Re a(s_J)>4/(5D), |a(s_J)|<6/(5D),
    Re a(s_K)>9/10,   |a(s_K)|<11/10.                  (L5)

For the first bound multiply the displayed error estimate by D and put
x=1/D<=1/21. It is at most
(1/10)(1+x/10)/(1-x/10)^2<1/5. The right side increases with x.
For the second use (1/100)(1-1/2100)^(-2)<1/10.
This explicitly proves the disk bounds for all D, rather than checking
only their centers.

Choose logs by continuation from w=0 for J and from1/w=0 for K.
Let S_J=Re Phi(s_J,w), S_K=Re Psi(s_K,w). Taylor's logarithmic remainder
|log(1+u)-u|<=|u|^2/[2(1-|u|)] gives

    |S_J-log H|<1/100,   |S_K+1|<1/100.                (L6)

For J the two errors are at most3/(500D) and5D*2^(-d): the latter
uses |w/s_J^d|<=4D*2^(-d). Their sum is below1/100 at D=21 and
decreases thereafter. For K they are at most1/(10000D^2) and5/W,
again less than1/100 by (L2). These computations specify the branches
and bound actual stationary heights, not only real envelope maxima.

## 3. Local Gaussian replacement with uniform constants

Take horizontal local radii rho_J=1/2 and rho_K=1/10. We shall use cubic
Taylor remainder bounds B_J=2/D^2 and B_K=1.

For J, on |s-D|<=3/5 put r=|w/s^d| and y=r/(1-r). Then
r<=6D*2^(-d) and y<=7D*2^(-d)<=1/100. The absolute third derivative
of the logarithmic correction is bounded by

    [d(d+1)(d+2)y+3d^2(d+1)y^2+2d^3 y^3]/|s|^3.

Since |s|>97D/100, its cube exceeds9D^3/10. The displayed expression
is at most2y<=14D*2^(-d)<=1/D^2, by (L2). The base term2D/s^3 has
modulus less than3/D^2. Thus the whole third derivative is less than
12/D^2, sufficient for B_J. For K, throughout its radius1/10 disk,
89/100<|s|<1101/1000. Here r=|s^d/w|<=3/(D5^d), y<=4/(D5^d).
The same third derivative bound, replacing d+1,d+2 by d-1,d-2,
is less than8D^2*5^(-d)<1; the base2/s^3 has modulus less than3.
Thus the whole third derivative is less than6, sufficient for B_K.

It follows from (L5) that the real local decay constants can be taken as

    c_J=1/(3D), c_K=1/3,
    c_*<=Re a(s_*)/2-B_*rho_*.                         (L7)

For J this is 2/(5D)-1/D^2>=1/(3D), since D>=15; for K it is
9/20-1/10>=1/3.

For clarity the elementary local error formula is recorded. If v=|s_*|,
the difference between the radius-rho local integral and the full Gaussian
leading term is at most

    e^(NS_*)/N *[B/((v-rho)c^2)
            +1/(v(v-rho)c)+1/(v c rho)].               (L8)

The exponential error is bounded by N B|u|^3 exp(-Nc u^2), whose integral
on the real line is B/(Nc^2). The amplitude difference is bounded by
|u|/[v(v-rho)], and the two completed Gaussian tails together by
exp(-Nc rho^2)/(Nc rho v). These prove (L8), including the completion
of the local integral to the whole real Gaussian. For the J data the
bracket is less than9; for the K data it is less than46. For example
the three J terms are bounded by19/D,4/D,7; the K terms by11,4,31.

The leading terms are

    L_J=e^(N Phi(s_J,w))/s_J *sqrt(2pi/(N a(s_J))),
    L_K=w^N e^(N Psi(s_K,w))/s_K *sqrt(2pi/(N a(s_K))).

Square roots use Re a>0. The absolute leading prefactors are between
2/sqrt(ND) and3/sqrt(ND) for J, and between2/sqrtN and3/sqrtN
for K/w^N. These follow from (L5), the saddle disks, and3<pi<22/7.

## 4. Every part of the upper contour

Deform [c,infinity) to the vertical segment from c to c+i Im s_J,
followed by the horizontal ray. The right vertical segment tends to zero
at infinity. The original integrand s^(N-1)(w+s^d)^N exp(-Ns) is entire,
so zeros of w+s^d on a nonlocal path create no contour singularity.

On the horizontal interval 3D/4<=x<=5D/4 let s=x+i Im s_J.
The correction ratio r<=3D*2^(-d), hence y<=4D*2^(-d). The base
D log s has real curvature at most -(99/100)D/x^2: if |Ims|/x<=1/100,
then Re(s^(-2))>99/(100x^2). The logarithmic correction contributes
at most

    D^2(y+y^2)/x^2<=4y<=16D*2^(-d)<=1/(10D).

Consequently

    Re Phi_ss<-[396/625-1/10]/D<-1/(3D).               (L9)

This gives the advertised uniform real curvature on the ENTIRE central
interval. Stationarity at x=Re s_J implies a height loss at least1/(24D)
outside the radius1/2 local interval. For x>=5D/4, the real phase derivative
is at most D/x-1+y<=-1/5+1/100<-1/10. At the joining point the distance
from Re s_J is at least D/5, so the height loss is at least D/150 and
in particular D/200. Thus the far upper ray costs at most

    10/(DN) *exp(NS_J-ND/200).                         (L10)

It remains to control the left part D/2<=x<=3D/4 and its vertical
endpoint. For |Ims|<=1/10, the absolute exponential base is at most

    exp(1/(50D))*x(3W+x^d)exp(-x).

Relative to H the x^D term is at most exp(-D/32), because
exp(9/32)<4/3. The 3W term is at most exp(-D/4): its maximum is
(3eD/2)exp(-D/2), and (3eD/2)<=exp(D/4) for D>=21.
For this last assertion check D=21 with a rational Taylor lower bound;
D/4-log(3eD/2) has positive derivative1/4-1/D thereafter.
Their sum is at most exp(-D/40), since

    log(1+exp(-7D/32))<=exp(-7D/32)<=D/160.

The latter inequality holds at21 by a quadratic exponential lower bound;
its left side decreases and its right side increases. Combining with
(L6), the left path loses at least D/50 from S_J, because
D/200>=1/100+1/(50D).

The left horizontal length and vertical segment contribute less than
exp(NS_J-ND/50). The off-local central interval contributes less than
(2/3)exp(NS_J-N/(24D)). Together with (L10), all nonlocal pieces cost
at most

    4exp(NS_J-N/(24D)).                                (L11)

This explicitly includes the vertical joining segment; it is not a
stationary-phase estimate that drops endpoint terms.

## 5. Every part of the lower contour

Keep [0,1/2] on the real axis, and deform [1/2,c] to height Im s_K,
retaining vertical segments at1/2 and c. On 1/2<=x<=3/2,
|s|<1501/1000, and

    |s^d/w|<=3/(D5^d),

because W>D(84/11)^d and (1501/1000)(11/84)<1/5.
The base log s has curvature at most -99/(100x^2); the correction has
modulus at most8D*5^(-d)/x^2<=1/(100x^2). Therefore

    Re Psi_ss<=-98/(100x^2)<=-98/225<-1/3.              (L12)

Off the radius1/10 local interval the loss is at least1/600. Its total
amplitude cost is at most2exp(NS_K-N/600).

On [0,1/2] the positive modulus envelope
x^(N-1)(1+3x^d/W)^N exp(-Nx) increases for N>=2. Its integral is
at most (8/25)^N: at the endpoint the base is less than8/25,
using exp(1/2)>8/5. On the vertical segment at1/2, the base is less
than33/100 and the integrated amplitude is less than
2|Im s_K|<=2/(100D)<1. By (L6) each of these
two contributions is bounded by exp(NS_K-N/100). These estimates handle
the apparent1/s amplitude at zero in the phase notation; the original
integrand is integrable there and is bounded directly.

For 3/2<=x<=D/2, the real modulus envelope for K/w^N is

    x exp(-x)+(3/e)(x/D)^D exp(D-x).                   (L13)

The first term is at most27/80, by exp(3/2)>40/9. The second is at
most (9/8)(sqrt(11)/4)^D, which is less than9/400 for D>=21.
Precisely, at21 its square comparison is2500*11^21<16^21; its base
sqrt(11)/4 is less than1, so it decreases thereafter. Thus (L13) is
strictly less than9/25. The imaginary correction from |Ims_K|<=1/(100D)
is at most exp(1/(45000D))<exp(1/10000). By the last inequality in (L3)
and (L6), the resulting loss is at least1/2000. The vertical segment at c
obeys the same bound. Including lengths and amplitudes, the entire
nonlocal lower contribution is at most

    4D exp(NS_K-N/2000).                               (L14)

In deriving (L13) no condition |s^d/w|<1 is imposed on the far interval.
This is a modulus bound for the original polynomial integrand and remains
valid when the phase logarithm cannot be continued along that interval.

## 6. Uniform saddle approximations and their analytic errors

For N>=10^12 D^10, both nonlocal bounds (L11),(L14) are at most
e^(NS_*)/N. Indeed exp(v)>=v^4/24 reduces these assertions to

    N^3>=96*24^4 D^4,   N^3>=96*2000^4 D,

respectively, immediate at the stated polynomial threshold. Combining
with (L8) gives

    J_N=L_J(1+epsilon_J),  K_N=L_K(1+epsilon_K),
    |epsilon_J|,|epsilon_K|<=epsilon_0=100D^2/sqrtN<1/100.                (L15)

The upper actual relative bound is at most5sqrtD/sqrtN; the lower is
at most47/(2sqrtN). The common bound in (L15) is deliberately weaker.
The errors are analytic on neighborhoods of their stated domains because
the original pieces and nonzero leading approximants are analytic there.
For the lower exterior domain the analytic factors are functions of1/w
and the explicit single-valued integer power w^N.

In particular the packet magnitudes satisfy

    H^N exp(-N/100)/sqrt(ND) <|J_N|<4H^N exp(N/100)/sqrt(ND),
    exp(-N-N/100)/sqrtN <|K_N/w^N|<4exp(-N+N/100)/sqrtN.               (L16)

## 7. All zeros lie in the overlap and are simple

The coefficients of K_N(w) are nonnegative. For |w|<=W/2,
|K_N(w)|<=K_N(W/2); (L16) then gives

    |K_N/J_N|<=4sqrtD*2^(-N)exp(N/50)<1.

The coefficients of J_N(w)/w^N as a polynomial in1/w are nonnegative.
For |w|>=2W use its value at2W to obtain

    |J_N/K_N|<=4/sqrtD*2^(-N)exp(N/50)<1.

These inequalities use log2>1/2 and N>=12sqrtD, which follows from the
threshold. Thus every zero of I_N lies in W/2<|w|<2W.

A disk of radius W/12 about a point of this annulus lies in both domains
(L1). Cauchy's estimate and (L15) give

    |epsilon_*'/(1+epsilon_*)|<=13epsilon_0/W<1/W.      (L17)

For s_J a Cauchy disk of radius W/2 inside |w|<3W gives
|s_J'|<=1/(5W). Since a'=D/s^2-2D/s^3, multiplying its absolute
bound by D and putting x=1/D gives
(1-x/10)^(-2)+2x(1-x/10)^(-3)<2 at x<=1/21; both terms
increase with x. Thus on its saddle disk |a'|<=2/D and Rea>=4/(5D),
so the s-logarithmic derivative of s^(-1)a^(-1/2) has modulus less
than3/2. Its w derivative is therefore less than1/W. For s_K, Cauchy
gives |s_K'|<=12/(100DW); here the bound
|a'|<=D[(2100/2099)^2+2(2100/2099)^3]<4D, Rea>=9/10, and the same
s-logarithmic derivative is less than3D. Again its w derivative is
less than1/W. These derivative bounds also follow directly from
stationarity, but Cauchy avoids large powers in their formulas.

Differentiating the stationary phases and combining (L17) yields

    |J_N'/J_N-N/(w+s_J^d)|<2/W,
    |K_N'/K_N-N/(w+s_K^d)|<2/W.                       (L18)

For the second identity include the derivative of w^N: it combines
with the explicit w derivative of Psi to give1/(w+s_K^d).
Now |s_K|^d<2 and W>2, hence |1/(w+s_K^d)|>1/(3W).
Also (D-1/10)^d>D^d/2 and2W<D^d/6 by (L2). Consequently

    |1/(w+s_J^d)|<3/D^d=3D exp(-d)/W<1/(100W).

The two main terms in (L18) are therefore separated by more than
N/(4W), while the total error is less than4/W. Since N>16 they
cannot coincide. At any zero of I_N=J_N+K_N, both pieces are nonzero
by (L15), K_N=-J_N, and

    I_N'=J_N*(J_N'/J_N-K_N'/K_N)!=0.

Every zero is simple. This proves UniformLargeGapSaddleSimplicity.

## 8. Overlap with the high-gap region and exact factorial scope

The retained high-gap theorem proves simplicity whenever

    N exp(-d^2/[4(d+1)])<=1/32,

in particular whenever d>=1+4log(32N). For d>=512 the polynomial
threshold N_0=10^12(d+1)^10 satisfies

    log(32N_0)<(d-1)/4.                                (L19)

At d=512, use10^12<2^40 and513<2^10, so32N_0<2^145.
The inequality e>8/3 and (8/3)^3>2^4 gives log2<3/4, hence
145log2<435/4<511/4. The difference between the two sides of (L19)
has derivative1/4-10/(d+1)>0 thereafter. Thus every terminal N>=2
is covered either by the retained high-gap theorem or by the uniform
large-N theorem just proved.

Formula(1) of GAP_SIX_RESULTANT_STRICT_POSITIVITY.md is an exact
all-gap fixed-shape identity: at ordinary beta=1 the consecutive lower-one
pair at terminal N is coprime if and only if A_N(z;N) is squarefree.
The normalizing factors are nonzero positive integers. Applying it proves
AllSufficientlyLargeLowerOneGaps. Monomial degeneracies A=0 or B=0 have
nonzero factorial moments immediately. The rational resultant remains
nonzero over every characteristic-zero coefficient field.

**Corollary ExactFinitePrefixReduction.** For each FIXED gap d>=20,
ordinary lower-one SFC at every terminal N>=2 is equivalent to its tests
at2<=N<10^12(d+1)^10. The forward implication is restriction; conversely
the finite tests handle those N, and UniformLargeGapSaddleSimplicity
handles the entire remaining tail. Using the all-terminal d>=512 result,
ordinary lower-one SFC throughout the RESTRICTED class d>=20 is therefore
equivalent to the single finite union of boxes

    20<=d<=511,   2<=N<10^12(d+1)^10.                 (L20)

Every test here is universal in the complex coefficients A,B, equivalently
nonvanishing of the corresponding rational resultant. These finite tests
have NOT been executed. This retained equivalence concerns d>=20.

The [complete ordinary lower-one theorem](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md)
now proves the assertion at EVERY positive gap and EVERY terminal. It
sharpens the uniform tail to N>=max(5000,128(d+1)), overlaps the older sparse
theorem atd44, and verifies the actual remaining171812 full-degree finite
modular tests. Thus all obligations of the earlier finite reductions are
resolved without enumerating their enormous boxes directly. The theorem
and contour estimates of this note retain their stated d>=20 domain;
the [gap-ten extension](FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md) remains a
valid intermediate antecedent. Arbitrary lower exponents and unrestricted
multivariate Factorial/Strong Factorial remain separate.

## 9. Evidence and source boundary

The checker verifies the exact rational numerical ledger, monotonic ratio
bounds, phase algebra identities and the overlap arithmetic. It does not
replace the full contour proof by numerical integration or a root scan.
The Gamma integral, Rouche's theorem, Cauchy's estimate and Gaussian integral
are used in their standard forms; the previous fixed-gap contour proof is
the method antecedent. The complementary high-gap theorem and exact
resultant normalization are retained proof inputs. No historical priority,
formal verification or external review is asserted.

Proof support: [exact source ledger](factorial_uniform_large_gap_saddle_source_ledger.json), [rational constants](factorial_uniform_large_gap_saddle_constants.json), [checker](verify_factorial_uniform_large_gap_saddle.py), and [receipt](factorial_uniform_large_gap_saddle_verification.json). Written arguments establish the uniform quantifiers; finite checks certify the stated exact ledger. All assurance remains internal.
