# Unitary groups with an orthogonal three-dimensional block: explicit product counterexamples in every dimension

The normalized classical HCIZ integral gives a uniform construction with one polynomial of degree four and the same unit marker in every dimension. Its complete confluent transform also gives all native rational radius moments. The argument uses the actual quotient coordinate S throughout.

**Theorem 1.** Let n>=3. In the standard compact rational model G=SU(n), let H=SO(3) be the first three-dimensional block and N=N_G(H) its FULL algebraic normalizer. For the actual invariant radius square S constructed below, put D=2n-3, T=1+t and

    p=T/t * [(1-T^2 S)*(1-T^2 S/4)]^2,
    F=p^D*(1+p),           H_marker=t^2.

On EACH actual affine quotient G/H and G/N, and after EVERY characteristic-zero extension of the native rational model, p^D and F are pure at every order, while the SAME unit t^2 detects every positive order. For every k>=D the exact marked seed value is

    Lambda(t^2 p^k)=beta_n*(81/16)^k*((2k)!)^2/(4k+1)!,
    beta_n=3*2^(2n-7)*(n-2)!*(n-1)!/(2n-5)!.

The binomial formula(11) gives every marked value of F. For EVERY integer a>=D, both p^a and F_a=p^a(1+p) are pure. A scalar monomial t^ell is unsafe for either family IFF ell>=2. It detects EVERY positive order of p^a IFF2<=ell<=a, and of F_a IFF2<=ell<=a+1. Thus2is the LEAST unsafe integer scalar-monomial exponent. Every F_a is absolutely nonpower even up to a unit and has primitive Laurent support on BOTH actual quotients. No sharp pure-power cutoff D is asserted for general n.

The same theorem holds for every specified imaginaryquadratic E=Qsqrt(-d) and positive rational symmetric form q3 directsum q_r for n>=4, where rank(q3)=3 and rank(q_r)=n-3, with G=SU(E,q), H=SO(q3) and the actual x,S in Section6. For n=3 it holds for every specified positive rational ternary q3 as well. These are actual rational models, not a classification of all rational forms.

Over every formally real field, the compact quotient factor and scalar split Laurent factor are individually Mathieu, but their product fails by F and t^2. For n>=4 the spaces G/H are nonspherical by dimension; the full-normalizer quotients are nonspherical by that test for4<=n<=7. All complete native moments are given in(M) below.


## 1. Actual compact models and radius

Let n>=3, m=n-2, and D=2n-3=2m+1. In the standard compact Q-model G=SU(n) for Q(i), let H=SO(3) act in the first three coordinates. For n>=4 put

    x=i diag(-1,0,...,0,1).

For n=3 use the same diagonal, now entirely in the three-dimensional block. Write z=Ad(g^-1)x and let K be the real antisymmetric upper-left three-by-three block of z. It is exactly the trace-orthogonal projection to h. The actual invariant regular function

    S=-2 tr(K^2)=v^2

is the squared SO(3)-radius. It descends to X=G/H; individual entries of K need not descend. Choose a standard rotation generator J in h with eigenvalues i,-i,0,...,0. The component s=-tr(zJ) has range [-2,2]. The SO(3)-action on h is the full rotation group. Haar orbital measure therefore has uniform conditional scalar component on [-v,v]. Its radius law will be computed completely below.

All normalized functionals are native rational Reynolds traces. Their real realizations are the indicated compact Haar laws, by invariance and uniqueness. The scalar Laurent functional Lambda is normalized trace followed by constant term in t. No finite-cover transport or generic descent of an arbitrary witness is used.

## 2. The complete confluent HCIZ determinant

Let a denote the Laplace variable and define

    E_m(a)=sum_(r=0)^(m-1) a^r/r!,
    C_n=m!(m+1)!/4.

Use the HCIZ formula for the two spectra (-1,0^m,1). In the confluent determinant order the middle derivatives as r=0,...,m-1. Its middle block is diagonal:

    B_rs=delta_rs*r!*a^r.

Each confluent Vandermonde has endpoint factor2 and middle factor product_(r=0)^(m-1) r!. The standard HCIZ numerator factorial is product_(j=1)^(m+1) j!. The two-by-two Schur complement has diagonal e^a-E_m(a) and offdiagonal e^-a-E_m(-a). Cancelling the middle determinant against these factors gives the EXACT identity

    E exp(a s)= C_n/a^D * [ (e^a-E_m(a))^2
                                  -(e^-a-E_m(-a))^2 ].                 (1)

This calculation initially uses a!=0 and distinct-spectrum approximants. Continuity of the compact integral supplies all confluent limits. In particular it is not a formal random-eigenvalue hypothesis. To check normalization directly, R_m(a)=e^a-E_m(a) starts with a^m/m!. The first odd term of R_m(a)^2-R_m(-a)^2 is

    4*a^(2m+1)/(m!*(m+1)!),

so the apparent singularity at a=0 in (1) is removable with value1.

For m=1, (1) is the previously known SU3 scalar law. For m=2 its bracket is

    e^(2a)-2(a+1)e^a+4a+2(1-a)e^-a-e^(-2a),

recovering the independently checked SU4 formula.

## 3. Global differentiability from the exact transform

Expand the bracket in (1). The e^(2a) and e^(-2a) coefficients are constants. The e^a and e^-a coefficients have degree at most m-1. The exponent-zero polynomial is

    E_m(a)^2-E_m(-a)^2.

It is odd; for m>=2 its degree is at most2m-3 because the highest possible degree2m-2 is even and cancels. For m=1 this polynomial is zero.

For every integer j>=1 and real c,

    e^(ca)/a^j = integral_(-infinity)^c e^(as)
                              (c-s)^(j-1)/(j-1)! ds,   a>0.           (2)

Consequently (1) has a complete spline inverse with knots only0,+/-1,+/-2 and polynomial pieces of degree at most D-1. Formula (2) initially gives an exponentially integrable signed density. Equality with the actual compact probability follows by analytic continuation in Re(a)>0 and Fourier uniqueness after multiplication by a fixed exponential weight. Thus this inverse is the actual compact scalar density rho_A; no moment-prefix identification is involved.

For m>=2 every nonzero term in this inverse has denominator exponent at least4: the exponent-zero terms have D-(2m-3)>=4, the e^(+/-a) terms have D-(m-1)=m+2>=4, and the outer terms have exponent D>=5. The truncated powers in (2) are therefore globally C2. In particular rho_A is C2 at0. For m=1 there is no zero knot, and rho_A is polynomial near0. In all cases rho_A is even, hence

    rho_A'(0)=0.                                                     (3)

This proves the needed regularity directly from the complete determinant. It imports neither a Morse-Bott theorem nor a differentiability assertion for an unspecified DH law.

For 1<s<2 only the e^(2a) truncated power survives, so

    rho_A(s)=C_n*(2-s)^(D-1)/(D-1)!.                                 (4)

All its coefficients and all complete moments are rational.

## 4. The complete radius-law properties needed by the circuit

Rotational invariance gives rho_A(s)=integral_[|s|,2] (2v)^-1 dnu(v). A radius atom at v>0 would give a jump of rho_A at that v, and a radius atom at0 would give a scalar atom. These are excluded by the actual continuous scalar density. On (0,2), its distributional derivative therefore gives the entire radius probability as

    dnu(v)=rho(v)dv,    rho(v)=-2v*rho_A'(v).                          (5)

In particular no unrecorded singular continuous radius part remains. The density is nonnegative because it is the actual probability. It is polynomial of degree at most D-1 separately on (0,1) and (1,2). By (3), its polynomial on (0,1) has no degree0 or degree1 term. Its exact outer polynomial is

    rho(v)=A_n*v*(2-v)^(D-2),
    A_n=m!(m+1)!/[2*(2m-1)!].                                       (6)

Thus the outer linear coefficient is A_n*2^(D-2)>0. We need no formula or sign claim for the individual higher homogeneous coefficients of the inner polynomial.

### Complete native rational moments

Let m=n-2, D=2m+1, C_n=m!(m+1)!/4, and S=v^2 as in that derivation. The COMPLETE exact moments are

    lambda(S^j)=2 C_n (2j+1)!/(D+2j)!
          * [2^(D+2j)-2 sum_(r=0)^(m-1) binom(D+2j,r)]               (M)

for every integer j>=0.

To prove (M), write R_m(a)=sum_(r>=m) a^r/r!. The bracket in the complete MGF is R_m(a)^2-R_m(-a)^2. Since D+2j is odd, its coefficient there is

    2 sum_(r=m)^(D+2j-m) 1/[r! (D+2j-r)!].

Multiplying by C_n*(2j)! gives E(s^(2j)). Conditional sphere integration gives E(s^(2j)|v)=v^(2j)/(2j+1). Finally the two omitted binomial tails r<m and r>D+2j-m are equal, giving (M). This is a direct extraction from the complete HCIZ transform, not a bounded moment-prefix extrapolation.

For n=3 formula(M) reduces to

    (4^(j+1)-1)/[(j+1)(2j+3)].

For n=4 it reduces to

    6[2*4^(j+1)-j-3]/[(j+1)(j+2)(2j+3)(2j+5)].


## 5. Complete endpoint cancellation and the same fixed unit

Set T=1+t and choose the SAME polynomial of degree four for every n:

    f(v)=(1-v^2)*(1-v^2/4),
    p=T/t*f(Tv)^2
      =T/t * [ (1-T^2 S)*(1-T^2 S/4) ]^2.                            (7)

The second expression defines an actual element of O(X)[t^+/-1]. On each interval I=(alpha,beta) equal to (0,1) or (1,2), write rho(v)=sum_d c_(I,d) v^d. For d>=1 put

    J_(I,d,k)(T)=integral_(T alpha)^(T beta) w^d f(w)^(2k) dw.

At endpoints1,2 the polynomial f has simple zeros; endpoint0 is fixed. Therefore

    J_(I,d,k)(T)-J_(I,d,k)(1)

is divisible by (T-1)^(2k+1). Substitution w=Tv in the averaged p^k gives t^-k T^(k-d-1) J_(I,d,k)(T). Extracting the coefficient of t^(-ell) proves the full exact identity

    Lambda(t^ell p^k)=sum_(I,d) c_(I,d)
         binom(k-d-1,k-ell)*integral_I v^d f(v)^(2k)dv,               (8)

for k>=1 and ell>=-k, with generalized integer binomial coefficients. Only T near1 is used in the integral interpretation; the resulting expression extends polynomially/formally. The endpoint remainder contributes no constant term precisely in the displayed domain.

For k>=D all upper indices k-d-1 are nonnegative and at most k-2. Thus Lambda(p^k)=0, and also Lambda(t^ell p^k)=0 for every ell<=1 in the domain of (8). Every fixed negative ell is eventually in that domain.

For ell=2 only degree1 survives. Section4 says that term occurs only on (1,2), and (6) supplies its coefficient. Substituting u=(v^2-1)/3 yields

    integral_1^2 v*f(v)^(2k)dv
       =(3/2)*(81/16)^k * ((2k)!)^2/(4k+1)!.

Hence for EVERY k>=D,

    Lambda(t^2 p^k)=beta_n*(81/16)^k*((2k)!)^2/(4k+1)! >0,
    beta_n=3*2^(2n-7)*(n-2)!*(n-1)!/(2n-5)!.                         (9)

In particular beta_3=3 and beta_4=12. The explicit phase

    F=p^D*(1+p),       H_marker=t^2                                 (10)

is pure at EVERY order. Its marked value at EVERY order is exactly

    sum_(h=0)^m0 binom(m0,h)*beta_n*(81/16)^(D*m0+h)
                            *((2D*m0+2h)!)^2/(4D*m0+4h+1)! >0.     (11)

All numbers in (11) are positive rational. The same exponent2 is the least unsafe integer monomial exponent for p^D and F, because every exponent<=1 is eventually zero by (8), while exponent2 survives every order. No sharp pure-power cutoff D is asserted for general n. The exact scalar-monomial profiles for every sufficient pure family a>=D are proved in Section7.

## 6. Full normalizer, direct geometry and rational forms

For n>=4, over Qbar the defining SO3-module splits as its irreducible three-dimensional module plus a trivial module of rank r=n-3. Its full commutant in SL_n is

    C={diag(z I3,W): z^3 det(W)=1}.

This is connected reductive: on GL1 times GL_r its defining character has primitive exponent vector(3,1). Every automorphism of so3 is inner by preservation of its Killing form and nonzero alternating bracket-volume. After adjusting a normalizer element by H, it centralizes h and hence H, so

    N_G(H)=H*C.

For n=3 the commutant is mu3, and the same normalizer conclusion holds with this finite centre. These are geometric equalities of characteristic-zero smooth subgroup schemes, and descend to the actual rational models.

The commutant fixes h pointwise, so trace-orthogonal projection to h is unchanged under its adjoint action. Thus S and the SPECIFIC phase (7),(10) are actual functions on the full-normalizer quotient too. This quotient is affine because its full stabilizer is reductive. Its native normalized trace is the restriction from G and has exactly the same displayed values. This argument descends this invariant phase, not an arbitrary upstairs witness.

On EACH actual quotient the identity coset has S=0 and p=1+t^-1. At t=-1/2 the modifier1+p has derivative-4 and p=-1 is a unit. For every integer a>=D, smooth geometric integrality gives an actual prime divisor of multiplicity1 in F_a=p^a(1+p), excluding a unit times any proper power over every algebraic closure. The same identity evaluation specializes F_a to (1+t^-1)^a*(2+t^-1); its coefficients at adjacent exponents0,-1 are2 and2a+1. Thus the literal Laurent support is primitive. These two geometric conclusions do not depend on HCIZ or the moment law.

For specified rational forms with n>=4, let E=Qsqrt(-d) be any imaginary quadratic field and q=q3 directsum q_r a positive rational symmetric form with q3 of rank3 and q_r of rank r=n-3. Let G=SU(E,q), H=SO(q3) upperblock. Choose rational diagonalizers P3,Pr and set

    x=sqrt(-d)*diag(P3 diag(-1,0,0) P3^-1,
                        Pr diag(0,...,0,1) Pr^-1).

This lies in the actual rational Lie algebra. Define K by real-part upperblock projection, and S=-(2/d)tr(K^2). Conjugating by the rational block diagonalizer and then the positive real square root of the diagonalized form preserves the chosen diagonal x. Dividing by sqrt(d) gives exactly the standard SU(n) radius above. Hence ALL complete moments and formulas(7)--(11) hold in these specified native rational models. For n=3 choose a rational diagonalizer P3 of q3 and x=sqrt(-d)*P3 diag(-1,0,1) P3^-1. The same rational and positive-real diagonal conjugacy commutes with this diagonal x and gives the standard n=3 law. This is the specified ternary construction, rederived here by the same actual conjugacy. There is no claim to classify all rational forms of type A.

The resulting rational pair remains pure/unsafe with the SAME unit after EVERY characteristic-zero extension. Compact factors have onlyzero pure elements over every formally real field by the rational Gram argument; the scalar Laurent factor is Mathieu by scalar DvK. Thus these are genuine two-factor instabilities over every formally real field.

For n>=4, dim(G/H)=n^2-4 exceeds dim(B_SLn)=(n^2+n-2)/2, proving these actual spaces nonspherical. The full-normalizer quotient has dimension6n-13; the same dimension test proves its nonsphericity only for4<=n<=7. No broader nonsphericity assertion is made.

## 7. Sharp scalar-monomial profiles for every sufficient pure family

Set m=n-2, D=2m+1, delta=D-2, C_n=m!(m+1)!/4 and A_n=2C_n/delta!. The following positive decomposition is derived from the complete scalar inverse already proved above. It supplies both directions of the marker profiles without classifying the pure-power threshold below D.

### 7.1. Positive inner-density decomposition from the complete inverse

For0<v<1, only the exponential terms with exponents1 and2 contribute to the truncated-power scalar inverse. Thus the COMPLETE inner scalar density is

    rho_A(v)=C_n[(2-v)^(2m)/(2m)!
                 -2 sum_(r=0)^(m-1) (1-v)^(2m-r)/(r!(2m-r)!)].       (S1)

Differentiate this actual formula and use rho(v)=-2v rho_A'(v). With u=1-v and delta=2m-1,

    rho(v)/(A_n v)
       =(1+u)^delta
           -2 sum_(r=0)^(m-1) binom(delta,r)u^(delta-r)
       =sum_(j=0)^(m-1) binom(delta,j)(u^j-u^(delta-j)).               (S2)

The final equality uses odd delta and pairs the lower and upper halves of the binomial expansion. For every index j in the sum, delta-2j>=1, so

    u^j-u^(delta-j)
       =(1-u)u^j sum_(h=0)^(delta-2j-1) u^h.

Since1-u=v, this gives the ACTUAL positive decomposition

    rho_inner(v)=A_n v^2 sum_(r=0)^(delta-1) b_r(1-v)^r,
    b_r=sum_(j=0)^min(m-1,r,delta-1-r) binom(delta,j) >0.             (S3)

In particular every coefficient b_r is a positive integer, including m=1. The known outer polynomial is

    rho_outer(v)=A_n v(2-v)^delta,       1<v<2.                       (S4)

No sign assertion about the alternating homogeneous coefficients was needed. Examples: the inner polynomials divided by A_n v^2 are1 for n=3, 1+4u+u^2 for n=4, and1+6u+16u^2+6u^3+u^4 for n=5.

### 7.2. Exact homogeneous transform and finite sampling interpretation

Keep the fixed polynomial f(v)=(1-v^2)(1-v^2/4), T=1+t and p=T/t*f(Tv)^2. For k>=D the exact local endpoint formula weights a homogeneous density term v^d by

    C_(d+1)=binom(k-d-1,k-ell).                                     (S5)

Fix integers2<=ell<=k. First consider the outer term(S4). Factor A_n v C_2. For0<=j<=delta,

    C_(2+j)/C_2 = (ell-2)_fallj/(k-2)_fallj.                         (S6)

Here C_2=binom(k-2,k-ell)>0. Since k-2>=delta, draw delta distinct objects uniformly without replacement from a population of k-2 objects, of which ell-2 are marked. Let I_1,...,I_delta indicate marked draws. Their joint inclusion moments are exactly the ratios in(S6), with numerator zero when too many marked draws are requested. The transformed outer polynomial is therefore

    A_n v C_2 * E product_(j=1)^delta (2-v I_j).                    (S7)

Every factor is strictly positive when1<v<2. Thus the outer transformed density is strictly positive there.

For the inner term(S3), ell=2 makes EVERY contribution zero because all its degrees are>=2. Suppose ell>=3. Factor A_n v^2 C_3, where C_3=binom(k-3,k-ell)>0. For each0<=j<=r<=delta-1,

    C_(3+j)/C_3 = (ell-3)_fallj/(k-3)_fallj.                         (S8)

Now draw r objects without replacement from k-3, with ell-3 marked. This is valid because k-3>=delta-1>=r. The transformed inner polynomial is

    A_n v^2 C_3 sum_(r=0)^(delta-1) b_r
                      * E product_(j=1)^r (1-v I_j).                (S9)

Every factor is strictly positive for0<v<1; an empty product is1. Thus(S9) is strictly positive there. This is an exact polynomial identity obtained by expansion and the joint moments, not a probabilistic approximation.

The endpoint integral additionally multiplies(S7),(S9) by f(v)^(2k). This multiplier is strictly positive in the two open intervals. Therefore

    Lambda(t^ell p^k)>0 whenever k>=D and2<=ell<=k.                  (S10)

For ell=2 this agrees with the previously established exact factorial beta value. No finite collection of moments is extrapolated.

### 7.3. Boundaries and the exact table within the endpoint domain

The Laurent polynomial p has minimum exponent at least-1. Hence

    Lambda(t^ell p^k)=0 for ell>k                                  (S11)

directly, before any sampling formula. The ratios above are not invoked outside their valid marked-population range. For ell<=1 and ell>=-k, the exact endpoint formula gives zero: every density degree is>=1, so every nonnegative upper index in(S5) is at most k-2 whereas the lower index k-ell is at least k-1. Every fixed negative ell eventually satisfies ell>=-k.

Thus for EVERY k>=D, throughout the exact domain ell>=-k,

    Lambda(t^ell p^k)>0  iff 2<=ell<=k,
    Lambda(t^ell p^k)=0  otherwise.                                (S12)

No fixed-k formula outside ell>=-k is asserted. This boundary matters for sufficiently negative ell.

### 7.4. Sharp unsafe and every-order profiles for sufficient pure families

Let a be ANY integer a>=D. The previously proved pure tail implies p^a and F_a=p^a(1+p) are pure at EVERY positive order. For either family, a scalar monomial t^ell is unsafe exactly when

    ell>=2.                                                        (S13)

Indeed ell<=1 is eventually zero. For any fixed ell>=2, all sufficiently large exponents k occurring in either family's powers satisfy k>=max(D,ell), so(S10) makes the marked values strictly positive. For F_a the coefficients in its finite binomial expansion are positive.

The every-order seed profile is exactly

    Lambda(t^ell (p^a)^m0)!=0 for EVERY m0>=1
                  iff 2<=ell<=a.                                  (S14)

Sufficiency follows from k=a*m0>=max(D,ell). If ell<=1 eventual zero prevents every-order detection; if ell>a, the first-order value is zero by(S11).

The every-order modified profile is exactly

    Lambda(t^ell F_a^m0)!=0 for EVERY m0>=1
                  iff 2<=ell<=a+1.                                (S15)

For sufficiency expand into exponents k=a*m0+h,0<=h<=m0. Every term with k<ell is zero by(S11); every term with k>=ell is positive by(S10), since k>=a>=D. The top exponent(a+1)*m0 is at least ell, so at least one positive term exists at EVERY order, including the first-order boundary ell=a+1. Necessity uses eventual zero for ell<=1 and the first-order Laurent minimum when ell>a+1.

All marked values in(S10),(S12)--(S15) are native rational numbers. Exact actual rational-model conjugacy and normalized trace restriction give the same identities on BOTH G/H and G/N_G(H), and every characteristic-zero extension preserves their stated zero/nonzero values. Real positivity is used to prove the rational values nonzero; no ordering of the extension field is assumed.

The pure-family condition a>=D is sufficient. The proof asserts NO sharp general pure-power cutoff below D and no classification of arbitrary nonmonomial markers. Actual order-one divisor/support certificates remain the independent identity-S0 argument from Section6, valid for every a>=D.

## 8. Classical input and scope

The moment input is Colin McSwiggen, [A new proof of Harish-Chandra's integral formula](https://arxiv.org/abs/1712.03995v4), exact arXiv1712.03995v4: normalized Equation2 on printed page2, Theorem1 on page4 and its complete Section3 heat-flow proof on pages5--10. The retained [primary PDF](sources/mcswiggen-harish-chandra-1712.03995v4.pdf) has SHA256 `dc50b3c21828d4f3a0099d95b8961ba5002ceca996776eea9c95449d14519a5d`. The source's classical radial-operator/Helgason, heat-uniqueness and steepest-descent inputs are retained; no Section4 localization theorem is used. Confluence, the spline inverse, its differentiability, all radius moments and every endpoint calculation required here are derived in this proof.

For native trace restriction, affine homogeneous quotients, rational Gram positivity and scalar Laurent Mathieu inputs, see [the homogeneous quotient trace theorem](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md), [the reductive homogeneous trace construction](REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md) and [the standard compact product constructions](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md). The [ternary unitary-orthogonal construction](TERNARY_UNITARY_ORTHOGONAL_EXPLICIT_RADIUS_COUNTEREXAMPLES.md) gives its additional dimension-three sharp profiles; its special badpair is not needed as a premise here.

This is a classical HCIZ application producing the displayed native-Q family and its scalar extensions. No CDKW or principal-face hypothesis, unspecified Morse-Bott regularity result, inference from finite moment prefixes, arbitrary-witness descent or new complex-existence classification is used. The all-field badpair and the separately Mathieu formally-real factors have different stated scopes.
