# Eventual Strong Factorial on every fixed integer cusp support

The complete cusp-neighbourhood and compact-complement arguments give eventual consecutive coprimality on every fixed integer cusp support. The cutoff may depend on the support. The remaining finite set of smaller terminal orders is not determined.

For integers s>=1 set

    a=s^2, d=2s+1, D=a+d=(s+1)^2,
    S_c=s(s+1), w_c=S_c^d*(s+1)/s.

For m>=1 define the monic rational moment polynomial

    A_(m;s)(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k), (E1)

where (t)_j=t(t+1)...(t+j-1), (t)_0=1. Every coefficient is an integer. Its constant term (a*m+1)_(d*m) is nonzero.

For a characteristic-zero field K, let Lambda_K:K[X]->K be the K-linear map Lambda_K(X^j)=j!. If B!=0 and f=X^a(A+B X^d), expansion gives the exact normalization

    Lambda_K(f^m)=B^m*(a*m)!*A_(m;s)(A/B).       (E2)

All integers and factorials here are interpreted through the canonical injection Q->K.

## 1. Main theorem

**Theorem EventualIntegerCuspStrongFactorial.** For every fixed integer s>=1 there exists an integer N_s>=2 such that

    gcd_(Q[x])(A_(N;s),A_(N-1;s))=1            (E3)

for every N>=N_s. The same N_s works for every characteristic-zero field K and every nonzero binomial phase

    f=X^(s^2)*(A+B X^(2s+1)), (A,B)!=(0,0):

for every N>=N_s the two values Lambda_K(f^N),Lambda_K(f^(N-1)) are not both zero. The cutoff may depend on s; it does not depend on K,A,or B. No numerical cutoff or uniform bound in s is supplied.

The proof proceeds through actual complex common roots, then descends algebraically through Q. It never assumes that an arbitrary characteristic-zero field embeds in C.

## 2. Exact analytic inputs and their scope

The following written statements are the premises of the synthesis.

1. **Fixed physical cusp neighbourhood.** CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md, Theorem FixedCuspNeighborhoodConsecutiveExclusion, provides r_s>0 and an unevaluated N_local(s) such that no actual common root x of A_(N;s),A_(N-1;s) has

       |x/N^d-w_c|<r_s

   for N>=N_local(s). This is a full fixed physical disc, including every scaled field within it. Its complete contour proof covers outer angular regions by a single Gaussian packet; in the central cone it uses an exact analytic two-period ratio and Cauchy derivative estimates. The bounded-field core uses the fixed-cusp bounded theorem as a sufficient route; the stronger uniform-in-s bounded theorem is optional. No global root-location premise is imported by this input.

2. **Every compact complement.** CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md, Theorem CompactComplementCuspConsecutiveExclusion, provides for each compact K subset of C minus {w_c} an unevaluated N_comp(s,K) such that no actual common root has x/N^d in K for N>=N_comp. This includes upper and lower half-planes, complete finite-sublevel Stokes detours, the entire real boundary away from the cusp, and zero. Its proof supplies actual complete integer-period estimates and a finite cover; no Stokes jump theorem or assumption of two active saddles everywhere is substituted.

3. **The retained global upper modulus.** GAMMA_PRODUCT_ROOT_MODULUS.tex, theorem thm:product-annulus, asserts for all integers n,d>=1 and real gamma>0 that each root t of

       T_(n,gamma,d)(t)=sum_(k=0)^n binom(n,k)*(gamma)_(d*k)*t^k

   satisfies |1/t|<=(gamma+d(n-1))_d. Only this old global UPPER bound is used. No new lower-annulus improvement or new sharp shape-separation theorem is a premise of the present result.

The actual own-order integral appearing in the first two inputs is

    I_m(w)=integral_0^infinity [S^a(w+S^d)]^m*exp(-m*S)dS
          =(a*m)!*A_(m;s)(m^d*w)/m^(D*m+1).    (E4)

Thus a common actual x means I_N(x/N^d)=I_(N-1)(x/(N-1)^d)=0, with the two different actual Gamma shapes. There is no fixed-shape replacement.

## 3. The compact root capture and finite covering

Apply the retained upper theorem to n=N,gamma=a*N+1. Since A_(N;s)(x)=x^N*T_(N,aN+1,d)(1/x) and x0 is not a root, every complex root of A_(N;s) satisfies

    |x|<= (D*N-d+1)_d <= (D*N)^d,
    |x/N^d|<=R_s:=D^d.                         (E5)

The second inequality follows term by term: the last factor in the rising factorial is D*N, and every earlier positive factor is at most D*N. The upper theorem's shape condition is simply gamma>0, satisfied here. Its additional conditional lower inequality and its strictness refinements are not needed.

Fix s. Take the r_s of Input1 and form the fixed compact set

    K_s={w in C: |w|<=R_s, |w-w_c|>=r_s/2}.    (E6)

It excludes w_c, so Input2 applies. Choose an integer

    N_s>=max(2,N_local(s),N_comp(s,K_s)).

All numbers and compact sets have been fixed before N or the coefficients A,B are varied. Their sizes need not be evaluated.

Suppose N>=N_s and x were a common complex root. By(E5), w=x/N^d belongs to the closed disk |w|<=R_s. If |w-w_c|<r_s, Input1 excludes it. Otherwise w belongs to K_s and Input2 excludes it. These alternatives cover every captured root; the overlap r_s/2<=|w-w_c|<r_s is harmless. Hence the two polynomials have no common complex root.

To make the compactness quantifier explicit, Input2 is obtained by choosing, at each point of K_s, a controlled neighbourhood for one of the regular-cell, Stokes-wall, or real-boundary arguments, shrinking it to allow adjacent radial transport, and taking a finite subcover. Its cutoff is the maximum of finitely many corresponding cutoffs. Thus the argument produces one N_s for the entire coefficient plane, not separate cutoffs for individual coefficient ratios.

A nonconstant common divisor over Q would have a complex root, contradicting this conclusion. Therefore(E3) follows.

## 4. Rational Bezout descent to all characteristic-zero fields

Fix N>=N_s. The Euclidean algorithm in Q[x] supplies U_N,V_N in Q[x] with

    U_N*A_(N;s)+V_N*A_(N-1;s)=1.               (E7)

Let K be any characteristic-zero field. Its unique prime-subfield map Q->K is injective, so(E7) remains a polynomial identity over K. If B!=0, evaluate it at A/B. Both A-polynomial values cannot vanish. Formula(E2), with its nonzero factorial and B-power factors, proves the claimed consecutive-moment exclusion. If B=0 and A!=0, every positive moment is A^m*(a*m)!, hence nonzero. The case A=0,B!=0 is also a nonzero monomial, or follows directly from(E2).

This is a field-theoretic argument after the complex analysis. It uses neither a norm nor an embedding of K into C. The same N_s therefore works over every characteristic-zero field and simultaneously for all nonzero coefficient pairs.

## 5. Genuine reverse equivalence: a finite bad-terminal set

For each fixed s define

    BadTerminal_s={N>=2 integer:
                  gcd_(Q[x])(A_(N;s),A_(N-1;s)) !=1}.        (E8)

**Corollary FiniteBadTerminalEquivalence.** For each fixed s the following are equivalent:

(i) BadTerminal_s is finite.

(ii) There is an N_s>=2 such that(E3) holds for every N>=N_s.

(iii) There is an N_s>=2 such that, over every characteristic-zero field and for every nonzero cusp binomial f of this fixed support, no consecutive pair at terminal order N>=N_s vanishes.

These statements are proved by the theorem. This is a reverse equivalence at the actual eventual-coprimality level, not a proposed inverse to an analytic error estimate or a finite checker.

Proof. If the bad set is empty choose N_s=2; otherwise finiteness permits N_s=1+max(BadTerminal_s). Conversely (ii) confines the bad set to {2,...,N_s-1}. Thus (i) and (ii) are equivalent. Statement(ii) implies(iii) by the rational Bezout proof of Section4. Conversely, if(iii) holds and a terminal N>=N_s had a nonconstant rational gcd, choose a root z of that gcd in an algebraic closure of Q. Over that characteristic-zero field, the nonzero polynomial f=X^a(z+X^d) would have both moments zero by(E2), contradicting(iii). This proves the reverse direction without an embedding premise.

The theorem establishes FINITENESS of BadTerminal_s. It does not establish that this set is empty, identify its elements, or verify its remaining finite range. In particular the full all-order Strong Factorial assertion at a support, which would require BadTerminal_s empty, remains a stronger unresolved statement here. The cutoff is not proved uniform in s, and no numerical value is claimed.

## 6. Interpretation of evidence

The analytic inputs contain complete contour, topology, Gaussian and Cauchy arguments; their finite exact/source checkers support those arguments and authenticate the retained bytes. The synthesis uses a compact cover and rational Bezout descent. No finite sample of moments or modular gcd computations is being promoted to an all-order inference. No new lower-annulus packet is needed for this proof route.

