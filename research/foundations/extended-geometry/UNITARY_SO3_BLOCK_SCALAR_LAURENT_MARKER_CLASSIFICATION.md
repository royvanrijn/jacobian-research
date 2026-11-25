# Complete scalar Laurent marker classification for the all-n unitary SO3-block circuit

The [all-dimensional unitary SO3-block construction](UNITARY_SO3_BLOCK_ALL_DIMENSIONS_EXPLICIT_RADIUS_COUNTEREXAMPLES.md) supplies an actual rational radius law and pure phases. Gaussian endpoint scaling separates the asymptotic sizes of scalar monomial markers and yields the complete eventual classification of their finite linear combinations.

The models are G=SU(E,q), where E=Qsqrt(-d) is any specified imaginaryquadratic field. For n=3, q=q3 is any positive rational symmetric ternary form and H=SO(q3). For n>=4, q=q3 directsum q_r with r=n-3, both forms positive rational symmetric, and H=SO(q3) in the upper three-dimensional block. Let N=N_G(H) be the FULL algebraic normalizer, and let X be either actual affine quotient G/H or G/N. The rational S and chosen rational Lie element are those explicitly constructed in the linked theorem; no arbitrary rational form is substituted.

Fix n>=3, m=n-2, D=2m+1=2n-3 and an integer a>=D. Put

    f(v)=(1-v^2)(1-v^2/4), T=1+t,
    p=T/t*f(Tv)^2,        F_a=p^a*(1+p),
    M_ell(k)=Lambda(t^ell p^k).

Here f(Tv)^2 means the actual polynomial [(1-T^2S)(1-T^2S/4)]^2 in the quotient coordinate S=v^2. Let K be ANY characteristic-zero extension of the native rational model. A scalar Laurent marker means an element q(t) of K[t,t^-1], with no X-dependent coefficients.

**Theorem 1.** For either pure family p^a or F_a, the scalar Laurent markers whose marked values are eventually zero are EXACTLY

    t*K[t^-1] = {q : every exponent in its finite support is <=1}.

The zero marker is included. Every other scalar Laurent marker has a nonzero marked value at EVERY sufficiently large order. This holds on BOTH actual quotients G/H and G/N_G(H), for all specified rational q/E data and every characteristic-zero extension. It does not classify every-order detection for nonmonomial markers or markers with X-dependent coefficients.

## 1. The positive leading inner coefficient

Use the complete positive inner decomposition, with delta=2m-1 and A_n=m!(m+1)!/[2(2m-1)!]:

    rho_inner(v)=A_n v^2 sum_(r=0)^(delta-1) b_r(1-v)^r,
    b_r=sum_(j=0)^min(m-1,r,delta-1-r) binom(delta,j)>0.

It follows that rho_inner(v)=c2*v^2+O(v^3), where c2=A_n sum_r b_r. Summing in the opposite order gives

    sum_r b_r=sum_(j=0)^(m-1)(delta-2j)binom(delta,j).

The identity

    (delta-2j)binom(delta,j)
       =delta[binom(delta-1,j)-binom(delta-1,j-1)]

telescopes. Therefore

    sum_r b_r=delta*binom(delta-1,m-1),
    c2=m^2*(m+1)/2 >0.                                               (1)

This includes n=3. The complete outer density remains A_n v(2-v)^delta on1<v<2.

Define the fresh positive asymptotic constant

    kappa_n=c2*sqrt(pi)/4*(2/5)^(3/2).                                (2)

This is not the HCIZ normalization constant C_n. Native exact moments remain rational; an irrational real asymptotic constant will only prove noncancellation of rational sequences.

## 2. Elementary endpoint integral estimate

For any fixed nonnegative integer d put

    I_d(k)=integral_0^1 v^d*f(v)^(2k)dv.

On0<=v<1, f(v)>0, f(0)=1 and f(v)<1 for v>0. At0,

    log f(v)=-5v^2/4+O(v^4).

Choose epsilon>0 small. There is c>0 with f(v)^(2k)<=exp(-c*k*v^2) on0<=v<=epsilon. On epsilon<=v<=1 its absolute value is bounded by eta^(2k) with a fixed eta<1. Thus Gaussian comparison and the substitution u=sqrt(k)*v give

    I_d(k)=O(k^(-(d+1)/2)).                                         (3)

For d=2, the rescaled integrand tends pointwise to u^2 exp(-5u^2/2) and is dominated by an integrable Gaussian polynomial on the small interval, extended by zero. The complementary interval is exponentially negligible. Dominated convergence proves the exact leading asymptotic

    I_2(k) ~ sqrt(pi)/4*(2/5)^(3/2)*k^(-3/2).                       (4)

The Gaussian integral is the elementary one obtained by squaring the Gaussian integral and using polar coordinates. No unspecified Laplace/Watson theorem is required.

## 3. Asymptotics of every fixed monomial marker

For fixed ell and all sufficiently large k, the complete local endpoint formula is

    M_ell(k)=sum_(I,d) c_(I,d)*binom(k-d-1,k-ell)
                              *integral_I v^d f(v)^(2k)dv.           (5)

The intervals are I=(0,1),(1,2), and rho_I(v)=sum_d c_(I,d)v^d. The exact domain ell>=-k eventually contains every fixed ell. For a fixed ell>=3, an inner degree d contributes only when d<=ell-1. In that case

    binom(k-d-1,k-ell)=binom(k-d-1,ell-d-1)
                        ~ k^(ell-d-1)/(ell-d-1)!.

The d=2 term, (1) and (4) supply

    c2*binom(k-3,ell-3)*I_2(k)
               ~ kappa_n/(ell-3)! * k^(ell-9/2).                    (6)

For each d>=3, (3) bounds the contribution by

    O(k^(ell-(3d+3)/2)),

which is smaller than (6) by at least k^(-3/2). All sums have finitely many d, independent of k.

On the outer interval, |f(v)|<=9/16: the quadratic in v^2 has its extremum at v^2=5/2. Put

    r0=81/256<1.

For fixed ell, the outer part of (5) is O(k^(ell-2)*r0^k). It is exponentially smaller than (6). We have proved

    M_ell(k) ~ kappa_n/(ell-3)! * k^(ell-9/2),
                         every fixed integer ell>=3.                (7)

For ell=2 the already exact beta formula gives, for k>=D,

    M_2(k)=beta_n*(81/16)^k * ((2k)!)^2/(4k+1)! >0.                 (8)

The beta integral is at most16^-k, so the useful bound

    0<M_2(k)<=beta_n*r0^k                                          (9)

holds. For every fixed ell<=1, M_ell(k)=0 for every sufficiently large k. These statements use the exact complete laws, not moment-prefix extrapolation.

## 4. Seed powers

Write S_ell(M)=Lambda(t^ell*(p^a)^M)=M_ell(aM). For every fixed ell>=3, (7) yields

    S_ell(M) ~ kappa_n/(ell-3)! * a^(ell-9/2)*M^(ell-9/2).           (10)

In particular successive ell have distinct powers of M. Equation(9) makes S_2(M) exponentially small, but strictly positive, while all fixed ell<=1 are eventually zero.

If q(t)=sum_(ell=u)^L q_ell*t^ell is a complex Laurent polynomial with q_L!=0 and L>=3, linearity and (10) show

    Lambda(q(t)*(p^a)^M)
       ~ q_L*kappa_n/(L-3)! * a^(L-9/2)*M^(L-9/2).                 (11)

Every lower exponent is negligible relative to this nonzero leading term. Therefore the marked value is nonzero for every sufficiently large M. If L=2, all lower terms are eventually zero and the value is exactly q_2*S_2(M) for large M, again nonzero. If all exponents are<=1, the finite sum is eventually zero.

## 5. Modified powers and uniform binomial averaging

Let B_ell(M)=Lambda(t^ell*F_a^M). Its literal binomial expansion is

    B_ell(M)=sum_(h=0)^M binom(M,h)*M_ell(aM+h).                     (12)

For fixed ell>=3 put alpha=ell-9/2 and c_ell=kappa_n/(ell-3)!. Equation(7) means M_ell(k)=c_ell*k^alpha*(1+epsilon_k), with epsilon_k->0. This convergence is automatically uniform over k in[aM,(a+1)M], since the supremum of |epsilon_k| over the entire tail k>=aM tends to zero.

Let Z_M be binomial(M,1/2). Dividing (12) by2^M*c_ell*M^alpha therefore gives

    E[(a+Z_M/M)^alpha]+o(1).                                       (13)

The power function is bounded and continuous on[a,a+1], since a>=D>=3. The elementary variance identity E(Z_M/M)=1/2 and Var(Z_M/M)=1/(4M), together with Chebyshev's inequality, gives convergence in probability to1/2; bounded uniform continuity then implies convergence of (13). Hence

    B_ell(M) ~ kappa_n/(ell-3)! * 2^M
                        *(a+1/2)^(ell-9/2)*M^(ell-9/2).            (14)

There is no need for a central limit theorem, a uniform-in-ell assertion, or an asymptotic in the dimension n.

Equation(9) also gives the direct bound

    0<B_2(M)<=beta_n*[r0^a*(1+r0)]^M.                               (15)

Its exponential base divided by2 is strictly less than1, so it is negligible compared with2^M times ANY fixed real power of M. Every fixed ell<=1 is eventually zero throughout (12), because its minimum exponent aM tends to infinity.

For a complex scalar Laurent marker of highest exponent L>=3, (14),(15) prove

    Lambda(q(t)*F_a^M)
       ~ q_L*kappa_n/(L-3)! * 2^M
                         *(a+1/2)^(L-9/2)*M^(L-9/2).               (16)

It is eventually nonzero. For L=2 positivity of (12) and eventual disappearance of the lower terms give the same nonzero conclusion. Markers supported in exponents<=1 are eventually zero.

## 6. Every characteristic-zero coefficient field

Let q have coefficients in an arbitrary characteristic-zero extension field K. Every monomial sequence M_ell(k), S_ell(M), B_ell(M) used above is the exact native rational sequence already proved on the specified actual model. Thus each marked value is a finite K-linear combination of known rational numbers.

The coefficient field K0=Q(q_u,...,q_L) is finitely generated. It embeds injectively into C: choose algebraically independent complex images of a finite transcendence basis, then extend over the remaining finite algebraic extension inside the algebraically closed field C. Only this finite coefficient field is embedded; no embedding of all K or ordering of K is assumed.

Under such an embedding sigma, the marked values become exactly the complex linear combinations in (11) or (16), with leading coefficient sigma(q_L)!=0. Their eventual nonvanishing in C implies eventual nonvanishing in K0 and K. If q is supported below2, eventual zero is already an exact rational identity and transfers directly. This proves the stated classification over EVERY characteristic-zero extension.

No actual variety or extension field has been replaced by a different real form. The complex embedding is applied solely to finite scalar coefficients in the known rational output sequences.

## 7. Both actual quotients and scope

The actual specified rational q/E models have exactly the same native radius moments. The particular S,p,F_a functions descend to G/N_G(H), and native traces there are restrictions of the same group trace. Thus every sequence and the entire classification above hold on BOTH actual quotients, independently of the special form of a scalar coefficient field.

The classification is eventual, not an every-order statement for arbitrary scalar Laurent polynomials. For example, since M_2(a)>0 and M_3(a)>0, the rational scalar marker

    q(t)=t^3-[M_3(a)/M_2(a)]*t^2

has zero first marked moment for p^a. Its highest exponent is3, so the theorem makes it nonzero at every sufficiently large order. The analogous cancellation can be arranged for F_a. Initial cancellations therefore prevent extending the monomial every-order range by a highest-support rule.

Markers with coefficients in O(X_K), rather than K, remain outside the theorem. Nothing here determines the pure-power cutoff below D. The complete HCIZ/radius law and exact local endpoint identity of the linked all-dimensional construction are the model-specific inputs. Section2 proves the Gaussian endpoint estimate by dominated convergence, and Section5 proves the required binomial averaging from its elementary variance bound.
