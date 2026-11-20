# Gap-ten extension of the uniform saddle theorem

20 November 2026. Internal written proof, with exact rational certificate support. Assurance remains internal.

**Theorem UniformGapTenSaddleSimplicity.** For every integer d>=10 and N>=10^12(d+1)^10, the polynomial

    A_N(z;N)=sum_(k=0)^N binom(N,k)(N)_(dk) z^(N-k)

has only simple complex roots.

**Corollary OrdinaryLowerOneFiniteBoxEquivalence.** The all-terminal Strong Factorial assertion for all ordinary binomials X(A+B X^d), d>=1, is equivalent to its finite restriction

    10<=d<=511,   2<=N<10^12(d+1)^10.

This is an exact finite reduction, not an execution of these tests. Its forward direction is restriction. Its reverse uses the theorem, the retained all-terminal d>=512 corollary, the retained gap-one discriminant formula (9) for d=1, the ordinary growing-shape prime-power theorem for d=2,3,4,5,7,8,9, and the retained full gap-six theorem. The first two inputs are ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex, formula (9) and Theorem prime-power-lower-one / Corollary prime-gap-lower-one: their shape is gamma=N, not the bounded-Gamma-shape constant-plus-power theorem. Coefficients may be complex or in any characteristic-zero field, by the same rational resultant identity. No general lower exponent or Gamma-shape statement is asserted.

Current followthrough: the [complete ordinary lower-one theorem](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md)
now resolves EVERY gap and terminal. Its practical tail, older sparse overlap
and exhaustive171812-case modular bridge establish all tests in the
predecessor finite box without enumerating that enormous box directly.
The theorem and arguments below retain their stated intermediate scope.

## 1. Imported proof and precise replacements

We import the convergent-integral, contour-deformation, local-error, analytic-packet, Cauchy-separation, and exact-resultant arguments of the retained FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md. The purpose of this note is to prove every bound in those arguments for 10<=d<20 too. Put

    D=d+1>=11, H=D^D exp(-D), W=eH=D^D exp(-d), c=3D/8.

Split the SAME integral I_N=K_N+J_N at c. The domains remain |w|<=3W for J and |w|>=W/3 for K. The stationary points, phases, amplitudes, stationary curvature identity, overlap annulus and N threshold are unchanged. Every modification to the previous estimates is proved below. All inequalities are uniform for D>=11, not sampled over a bounded range.

Use 163/60<e<87/32<11/4 and q=3/8. Then

    exp(-d)<q^d, W>D*4^d>400D^2,
    40D^2 q^d<1, 14D^3 q^d<2, 300D q^d<1,
    800D (2/5)^d<1, 8D^2 (2/5)^d<1.                 (T1)

Check each at d=10. The relevant successive ratios are bounded by (12/11)^3*(3/8)<1 and (12/11)^2*(2/5)<1. For 4^d>400D its reciprocal ratio is (D+1)/(4D)<1. The Euler lower bound is the sum through degree5 of its factorial series (with the positive tail). The upper bound is the retained rational certificate.

## 2. Saddles, heights and local packets

The upper Rouche displacement on |s-D|=1/10 is <4D^2 q^d<1/10. The lower displacement is <4D/W<1/(100D). Thus the same simple analytic saddles exist on neighborhoods of both closed domains. In the curvature estimates replace x=1/21 by x=1/11 and 1/2100 by 1/1100. Direct rational bounds still give

    Re a(s_J)>4/(5D), |a(s_J)|<6/(5D),
    Re a(s_K)>9/10, |a(s_K)|<11/10.

The upper stationary height error is now bounded by 3/(500D)+5D q^d<1/100; the lower error bound is unchanged, with D>=11. The leading amplitude bounds 2..3 hold after replacing 210/209,210/211 by 110/109,110/111 and 2100/2099,2100/2101 by 1100/1099,1100/1101. These rational inequalities are explicitly checked; no old D>=21 prefactor assertion is simply imported.

Keep rho_J=1/2, rho_K=1/10. Replace only the upper cubic-remainder and decay data by

    B_J=3/D^2, c_J=1/(5D).                            (T2)

For |s-D|<=3/5, r=|w/s^d|<=6D q^d, and y=r/(1-r)<=7D q^d<1/100. Also |s|>47D/50. The log correction's third derivative, as explicitly differentiated in the imported proof, divided by y is at most

    [(12/11)+3/100+2/10000]*(50/47)^3 <2.

Here d(d+1)(d+2)<= (12/11)D^3, 3d^2(d+1)<=3D^3, and 2d^3<=2D^3. The base term 2D/|s|^3<3/D^2. Therefore the whole third derivative is <3/D^2+14D q^d<5/D^2<18/D^2, which is more than sufficient for the remainder coefficient B_J. Furthermore

    (2/5)/D-(3/D^2)(1/2) >=1/(5D)

for D>=11. Thus the retained local formula (L8) applies, and its upper bracket is <20 (not the old value9). Indeed its three terms are bounded by 75/(D-3/5), 5D/[(D-1/10)(D-3/5)] and 10D/(D-1/10), a decreasing total whose value at11 is <20.

For K replace each occurrence of 5^(-d) in its local and central correction bounds by (2/5)^d. This follows from W>D4^d and 1101/4000<2/5,1501/4000<2/5. Hence r<=3/[D(5/2)^d], y<=4/[D(5/2)^d]. Its log third derivative is still bounded by 8D^2(2/5)^d<1 and the base third derivative by3. Thus B_K=1,c_K=1/3 and the bracket<46 are unchanged.

## 3. Complete upper contour

On 3D/4<=x<=5D/4, the same ratio bound r<=3D*2^(-d) holds since e>8/3. Put r_0=3D*2^(-d), y_0=r_0/(1-r_0). At d=10,

    D(y_0+y_0^2)=371712/982081 <2/5.                 (T3)

It decreases with d: D*y_0=3D^2*2^(-d)/(1-r_0) and D*y_0^2=9D^3*4^(-d)/(1-r_0)^2 are decreasing by the same successive-ratio test. Keep the curvature factor D/x^2 instead of weakening it early. The real curvature is at most

    -[99/100-D(y_0+y_0^2)]D/x^2
    <=-(59/100)*16/(25D)<-1/(3D).

The imaginary-to-real ratio is at most2/165<1/80. The same base curvature factor99/100 holds, since (1-1/80^2)/(1+1/80^2)^2>99/100. The off-local loss is therefore at least1/(24D). On x>=5D/4 use the sharper r<=3D(3/10)^d, so y<1/100, giving derivative<-1/10 and far boundary loss>=D/200. These bounds continue to hold on the entire far ray.

It remains to bound the expanded left interval [3D/8,3D/4] and its endpoint vertical. For |Im s|<=1/10, the modulus base divided by H is at most exp(1/(25D)) times

    (x/D)^D exp(D-x)+3e*x*exp(-x).                   (T4)

The imaginary exponent follows from D*(1/10)^2/[2(3D/8)^2]=8/(225D)<1/(25D). Split at D/2. On the first interval the two terms are bounded respectively by

    (sqrt(e)/2)^D <=1/8,
    (9eD/8)exp(-3D/8) <=9/16.

On the second interval they are bounded respectively by

    ((3/4)exp(1/4))^D <=2/3,
    (3eD/2)exp(-D/2) <=3/16.

For the high terms the bases are <1, so check D=11 by squaring/fourth-power rational inequalities using e<87/32. For the low terms D exp(-alpha D) decreases for D>=11 (or use the ratio bound (12/11)/(1+alpha)<1); at11 their claims follow from exp(25/8)>22 and exp(9/2)>88. The former follows from (163/60)^3*(9/8)>22, and the latter from (163/60)^9>88^2.

Both sums are <7/8. Combining (T4) with S_J>=log H-1/100 gives a height loss greater than1/10, since

    exp(11/100+1/(25D)) <=1/(1-5/44)=44/39 <8/7.

The integrated amplitude of the entire left horizontal plus vertical path is less than log2+(1/10)/c<1. For example log2<3/4 and (1/10)/c<=4/165<1/4. Thus its contribution is <exp(NS_J-N/10). This is stronger than needed for the retained common nonlocal bound4exp(NS_J-N/(24D)). Together with the unchanged central and far upper bounds, every piece of the upper contour is controlled.

## 4. Complete lower contour

The paths near zero and x=1/2 and the central interval[1/2,3/2] are unchanged. The correction curvature bound becomes

    8D(2/5)^d/x^2<1/(100x^2),

by (T1), preserving Re Psi_ss<=-98/(100x^2)<-1/3. The zero endpoint and left vertical estimates are unchanged and hold with D>=11.

For the larger remaining interval[3/2,3D/8], the real envelope is still

    x exp(-x)+(3/e)(x/D)^D exp(D-x).

The first term is <27/80. The second increases up to x=D, so its maximum on this interval is

    (3/e)*[(3/8)exp(5/8)]^D <9/400.                 (T5)

Its exponential base is less than1. At D=11 the strict bound follows by taking eighth powers from the exact rational inequality

    (9/8)^8*(3/8)^88*(87/32)^55 <(9/400)^8.

Here 3/e<9/8. The left/right sides' ratio is about0.9927, but only the exact rational comparison is used. Hence the real envelope is <9/25 as before. The horizontal imaginary correction is <=exp(1/(45000D)), and the right vertical at c obeys the SAME bound. The retained exp(1+1/100+1/1000)<25/9 then gives the same loss1/2000. Integrated amplitudes are bounded by4D. Consequently all nonlocal lower pieces still cost at most4D exp(NS_K-N/2000).

## 5. Threshold and simplicity

No new asymptotic threshold is needed. At N>=10^12D^10, the old fourth-power exponential estimates control the nonlocal upper/lower tails by e^(NS_*)/N for every D>=11. Combine with the local brackets<20 and<46. The same common analytic relative error

    epsilon_0=100D^2/sqrt(N)<1/100

is valid, with much room. The zero-free inner/outer disk arguments use only these errors, heights and prefactors; they are unchanged. The analytic Cauchy estimates on the overlap annulus W/2<=|w|<=2W and the prefactor logarithmic derivative estimates remain valid with D>=11; every newly limiting rational bound is checked at D11. Finally use

    3/D^d<1/(100W),

which follows from 300D q^d<1 in (T1). Thus the same main logarithmic derivatives differ by more than N/(4W), their combined errors are less than4/W, and N>16. Every zero of I_N, hence of A_N, is simple.

This proves the theorem. The exact physical-shape resultant then proves the consecutive-moment tail and the stated finite-box equivalence. The cutoff d>=512 for every N is the ALREADY retained corollary, not a new numerical optimization in this note.

Proof support: [source ledger](factorial_gap_ten_saddle_source_ledger.json), [rational constants](factorial_gap_ten_saddle_constants.json), [checker](verify_factorial_gap_ten_saddle_extension.py), and [receipt](factorial_gap_ten_saddle_verification.json). Written contour arguments establish the uniform quantifiers; the exact replay certifies their numerical ledger. The finite global box is not executed by this packet.
