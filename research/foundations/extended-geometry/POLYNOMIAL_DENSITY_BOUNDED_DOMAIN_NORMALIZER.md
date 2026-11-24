# Polynomial-density bounded-domain normalizers

A polynomial density gives a finite sum of homogeneous endpoint terms. Only its constant homogeneous part contributes at the distinguished marker. A uniform bound on the other homogeneous moments proves an exact unsafe-monomial cutoff; the proved zero tail also reduces purity of every seed power to finitely many exact moments.

## 1. Conditional functional and geometric hypotheses

Let Omega be a bounded nonempty open subset of R^r with 0 in Omega and r>=1. Let rho be a REAL polynomial of total degree b, nonnegative on Omega, with rho(0)>0. Put Z=integral_Omega rho>0. Let f in Q[u_1,...,u_r] satisfy f(0)=1 and vanish on the ENTIRE topological boundary of Omega. Let A be a Q-algebra with actual elements x_1,...,x_r and a normalized Q-valued linear functional lambda satisfying

    lambda(P(x))=Z^-1 integral_Omega rho(u)P(u)du

for EVERY polynomial P with rational coefficients. Compatible characteristic-zero scalar extensions of this given algebra, functional and coordinates are an additional hypothesis. No rationality of rho's coefficients, of its separate homogeneous moments, or of the unweighted domain volume is assumed.

Write rho=sum_(j=0)^b rho_j with rho_j homogeneous of degree j and rho_0=rho(0)>0. Put

    T=1+t, p=T/t f(Tx)^2, a=r+b,
    F=p^a(1+p), H=t^r, Lambda=lambda tensor CT_t.

These are actual regular Laurent expressions overQ and H is a unit. For nonpower and primitive support assume ADDITIONALLY that Spec(A) is smooth and geometrically integral and x has an ACTUAL common geometric zero. The functional law alone supplies none of those geometric hypotheses.

**Theorem 1.** Under the functional hypotheses, F is pure and the ONE fixed unit marker H survives EVERY positive order, with

    Lambda(F^m)=0,
    Lambda(H F^m)=rho_0/Z sum_(h=0)^m binom(m,h)
                           integral_Omega f(u)^(2(am+h))du>0.

Every displayed marked value is rational by the native functional, hence remains nonzero after every compatible characteristic-zero extension.

For each positive integer c, p^c is pure IF AND ONLY IF Lambda(p^(cm))=0 for the finitely many positive m with cm<a. Equivalently, with B_bad={k in{1,...,a-1}:Lambda(p^k)!=0}, the positive pure exponents are precisely the complement in Z_(>0) of the union of the positive divisors of the elements of B_bad. Every c>=a is pure; that bound need not be necessary.

For F, and for every seed power p^c with c>0, a fixed Laurent monomial t^ell is unsafe IF AND ONLY IF ell>=r. On that side its moments are positive at every sufficiently large order; for ell<r they eventually vanish. For c>=a the unit marker t^r survives EVERY positive order of p^c, as it does for F. No every-order iff for the other markers is asserted.

Under the EXTRA geometric hypotheses, F is not a unit times any proper power over an algebraic closure, and its Laurent support differences have gcd1. Each assertion propagates under the stated compatible scalar extensions within its own hypotheses. This theorem is a conditional normalizer criterion, not an existence theorem for an actual model or a global Mathieu classification.

## 2. Homogeneous Euler zero-flux and the exact mixed identity

The [bounded-domain Fubini lemma](UNIFORM_BOUNDED_DOMAIN_NORMALIZER_CRITERION.md), Section2, says that a polynomial vector field vanishing on the whole boundary of ANY bounded open domain has integral of its divergence zero. Its proof needs only coordinate-line fundamental calculus and Fubini: each open bounded line slice is a countable union of intervals with endpoints on the domain boundary; endpoint values vanish, and bounded derivatives justify absolute summation and integration. No regular boundary, finite perimeter, star-shapedness or finite component count is required.

For k>=1 and 0<=j<=b put

    J_(k,j)(T)=T^(r+j)/Z integral_Omega rho_j(u)f(Tu)^(2k)du.

This is a polynomial in T. Differentiating and using Euler's identity u dot grad(rho_j)=j rho_j gives

    J'_(k,j)(T)=T^(r+j-1)/Z
                   integral_Omega div_u[u rho_j(u)f(Tu)^(2k)]du.

Every coefficient of t^n, n<2k, in u rho_j(u)f((1+t)u)^(2k) is divisible by f(u)^(2k-n), hence vanishes on the whole boundary. Apply the zero-flux lemma coefficientwise. Multiplication by T^(r+j-1) cannot lower the order, so

    J_(k,j)(1+t)-J_(k,j)(1) is divisible by t^(2k+1).

Set nu_(k,j)=Z^-1 integral_Omega rho_j(u)f(u)^(2k)du. These individual numbers need not be rational or nonnegative. Coefficient extraction yields the exact formula

    Lambda(t^ell p^k)
       =sum_(j=0)^b nu_(k,j) binom(k-r-j,k-ell),
       k>=1, integer ell>=-k.                             (2.1)

Indeed the left side extracts coefficient t^(k-ell) in sum_j T^(k-r-j)J_(k,j)(T). The selected degree is at most2k. Even if k-r-j is negative, T=1+t is a formal unit and multiplication cannot lower the endpoint remainder's order. A negative lower binomial index means0. Formula(2.1) is not asserted at ell<-k.

## 3. All-order purity and the distinguished fixed unit

For every k>=a=r+b the upper indices k-r-j in(2.1) are nonnegative. At ell=0 they are smaller than k, so every binomial is zero. At ell=r only j=0 survives, with binomial1. Therefore

    Lambda(p^k)=0,
    Lambda(t^r p^k)=nu_(k,0)
       =rho_0/Z integral_Omega f(u)^(2k)du>0, k>=a.        (3.1)

Strict positivity holds because f(0)=1 on a neighborhood of an interior point. Expanding F^m=sum_(h=0)^m binom(m,h)p^(am+h) proves the stated all-order phase and fixed-marker formulas. The nonzero values are rational because F,H and lambda are rational; no separate unweighted-moment rationality is needed.

## 4. Uniform moment domination and the exact unsafe cutoff

Let M_j=max_(closureOmega)|rho_j| for j>=1; these are finite. Since f^(2k)>=0,

    |nu_(k,j)| <= (M_j/rho_0)nu_(k,0) for EVERY k>=1.      (4.1)

Fix ell>=r and put d=ell-r. For k>=max(a,ell) set n=k-r. Formula(2.1) becomes

    Lambda(t^ell p^k)
       =sum_(j=0)^min(b,d)nu_(k,j)binom(n-j,d-j).          (4.2)

If d=0 this is exactly nu_(k,0)>0. If d>=1, then n>=d>=1 and, for 1<=j<=min(b,d),

    binom(n-j,d-j)/binom(n,d)
      =d(d-1)...(d-j+1)/[n(n-1)...(n-j+1)] <=d/n.

With C_d=sum_(j=1)^min(b,d)M_j/rho_0, equations(4.1)-(4.2) imply

    Lambda(t^ell p^k)
      >=nu_(k,0)binom(n,d)(1-d C_d/n)>0

whenever k>=max(a,ell) and k-r>d C_d. Thus every fixed ell>=r has positive moments for ALL sufficiently large k. Signed homogeneous density pieces cannot cancel persistently at a fixed marker exponent. No global-maximizer or asymptotic-concentration assumption on f is used.

If ell<r is fixed and k>=max(a,-ell), formula(2.1) applies with every upper index nonnegative and smaller than k-ell. All terms vanish. Substituting k=cm proves BOTH unsafe directions for each positive c. For F, once am exceeds the same bound, every term in the binomial expansion is positive when ell>=r and zero when ell<r. This proves its exact unsafe iff too. The distinguished marker t^r retains the stronger EVERY-order statement of Section3; other markers may cancel at initial orders, and no every-order classification is inferred.

## 5. A necessary and sufficient finite pure-power reduction

For c>0, the pure moments of p^c are Lambda(p^(cm)). Section3 proves all of them zero when cm>=a. Hence purity is equivalent to the finitely many exact tests with cm<a. Conversely those tests plus the proved unrestricted tail cover EVERY order. Failure occurs precisely when c divides some k in B_bad, proving the stated positive-divisor formula.

This is a complete finite reduction backed by an unbounded tail proof, not extrapolation from a numerical prefix. It does not assert that an arbitrary bad-order set is realizable. All tests are values of the native rational functional, not separate homogeneous-density moments.

The sufficient bound c>=r+b need not be necessary. On Omega=(-1,1) take rho(u)=1+u/2, b=1, and f(u)=1-u^2. The odd density piece integrates to zero against every coefficient of the even-in-u phase p^k, while the normalized even part is the uniform interval law. The rank-one uniform identity makes p itself pure, although r+b=2. Thus the sufficient bound is not a general iff threshold; the finite test is the exact replacement.

## 6. Actual nonpower, support and scalar fields

At an actual common geometric zero x0,

    p(x0,t)=(1+t)/t, 1+p(x0,t)=2+t^-1.

At t=-1/2 the latter has derivative-4 and p=-1 is a unit. On the smooth actual product it cuts a prime divisor of multiplicity1, while p has valuation0 there. Hence F=p^a(1+p) cannot be a unit times any proper power over the algebraic closure. Its specialization is

    F(x0,t)=(1+t^-1)^a(2+t^-1),

with adjacent coefficients2 and2a+1 at exponents0 and-1. The actual coefficient functions are therefore nonzero and the Laurent support differences have gcd1. Neither conclusion follows from the density law alone.

All actual rational phase/marker moments lie inQ. Positivity overR proves nonzero rational values at the stated orders; the zeros and nonzeros persist after EVERY compatible characteristic-zero scalar extension. This uses only the native functional and Q embedding in a characteristic-zero field, not rational coefficients of rho, rational unweighted volume, or an embedding of the whole field intoR. The generic scalar-transfer argument in [the bounded-domain criterion](UNIFORM_BOUNDED_DOMAIN_NORMALIZER_CRITERION.md), Sections1 and3-5, has precisely these hypotheses. The geometric conclusions separately replay over Kbar under the extra smoothness and actual-common-zero hypotheses.

A rational translation q in Omega with rho(q)>0 reduces to the stated constant-term condition. A rational boundary polynomial nonzero atq can be normalized by its rational value atq. For the geometric conclusions an ACTUALq-fiber is still needed; neither a density nor a positive moment law supplies such a fiber by itself.

The reused mathematical input is the arbitrary-bounded-open Fubini zero-flux lemma and its generic conditional scalar-transfer argument. No polynomial-density law on a specific homogeneous model, new primary theorem, finite experiment, checker or global Mathieu classification is a premise of this conditional written result.
