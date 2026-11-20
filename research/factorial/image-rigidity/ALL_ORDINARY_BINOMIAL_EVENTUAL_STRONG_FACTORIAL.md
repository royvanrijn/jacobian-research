# Eventual Strong Factorial for every fixed ordinary binomial support

Every fixed ordinary one-variable binomial support has only finitely many bad terminal orders. The terminal cutoff is uniform over the coefficients and all characteristic-zero fields at that fixed support. No effective cutoff, bound uniform over different supports, or all-order Strong Factorial theorem is asserted.

For integers0<=a<b put d=b-a>=1. For m>=1 define

    A_(m;a,d)(x)
       =sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k).    (M1)

The rising factorial has convention(t)_0=1. Each A_m is a monic integer polynomial of degree m with positive constant term. Let Lambda_K:K[X]->K be the K-linear map Lambda_K(X^j)=j! for a characteristic-zero field K.

## 1. Main theorem

**Theorem EveryFixedOrdinaryBinomialEventualStrongFactorial.**
For every fixed integer support0<=a<b there is an integer N_(a,b)>=2 such that

    gcd_(Q[x])(A_(N;a,d),A_(N-1;a,d))=1                 (M2)

for every terminal order N>=N_(a,b). Equivalently, the two polynomials have no common complex root at those orders.

The same N_(a,b) works for EVERY characteristic-zero field K and EVERY coefficient pair A,B in K with(A,B)!=(0,0). For the nonzero polynomial f=A X^a+B X^b,

    (Lambda_K(f^N),Lambda_K(f^(N-1)))!=(0,0)             (M3)

at every N>=N_(a,b).

The quantifier is for each fixed support there exists one cutoff for all fields and coefficients. It is not a cutoff chosen separately for each coefficient ratio. The theorem supplies no cutoff formula for all supports and does not identify or eliminate the smaller exceptional terminal orders in general.

## 2. Exact actual-order normalization

Writing f=X^a*(A+B X^d), expansion gives, when B!=0,

    Lambda_K(f^m)
      =B^m*(a*m)!*A_(m;a,d)(A/B).                      (M4)

For every positive integer m the analytic period is

    I_m(w)=integral_0^infinity [S^a*(w+S^d)]^m
                                      *exp(-m*S)dS
          =(a*m)!*A_(m;a,d)(m^d*w)/m^((a+d)*m+1).       (M5)

Consequently a common actual root x at orders N,n=N-1 corresponds to

    I_N(x/N^d)=I_n(x/n^d)=0.

At order m, the actual Gamma shape is am+1. The lower own-order parameter is(N/n)^d times the upper one. Every analytic input below treats this exact coupling.

If one coefficient is zero, f is a nonzero monomial, including the possibility of a nonzero constant. Every positive moment is its nonzero coefficient power times a factorial. Those cases never produce a bad consecutive pair.

## 3. Complete exhaustion of the integer supports

The following cases are disjoint after their order is fixed and exhaust0<=a<b.

1. **a=0.** Only a minimal analytic-tail route is needed, not the all-order constant-plus-power finite prefix. In the notation of ALL_UNIVARIATE_TWO_TERM_SFC.tex the reversed polynomial P_(m,d) is exactly A_(m;0,d). Its Theorem thm:analytic proves simplicity whenever m>=10000 and6<=d<=sqrt(m). For each fixed d>=6 it therefore applies at every terminal N>=max(10000,d^2). The exact identity P'_(N,d)=N*P_(N-1,d) gives actual consecutive coprimality there. For d=2,3,4,5, the all-order prime-power theorem in PRIME_POWER_TWO_TERM_SFC.tex supplies simplicity, since all four exponents are prime powers. The remaining d=1 is covered by Theorem thm:main in UNIFORM_PADE_PARABOLA_SFC.tex at a=0. These statements exhaust the constant-plus-power supports. No92935-case finite prefix, finite-prefix checker or certificate is a premise of this route.

2. **a>=1,d=1.** UNIFORM_PADE_PARABOLA_SFC.tex, Theorem thm:coprime and its transfer section, proves the adjacent-monomial assertion at every order, in fact for every positive Gamma shape. For ordinary shape1 use M=aN and l=a(N-1), with M-l=a>=1. Its polynomial R_m satisfies A_(m;a,1)(x)=(-1)^m R_m(am,-x), so this is the actual pair in(M2), and N_(a,a+1)=2 is sufficient.

3. **a>=1,d>=2 and4a<=(d-1)^2.** SQUARE_ROOT_GAP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md, Theorem SquareRootGapEventualStrongFactorial, gives(M2). It combines the strict subcritical support proof and the critical consecutive-square cusp proof. Its cutoff is uniform over the coefficient plane at this fixed support. It uses only the retained old Gamma upper bound for global root capture.

4. **a>=1,d>=2 and4a>(d-1)^2.** These are the remaining above-critical supports. Their two conjugate simple turning values are handled by the new local-neighborhood theorem, and their regular compact complements by the separate above-critical contour theorem. Section4 gives the complete root-capture and finite-cover argument.

No all-order least-exponent-one result is needed as an extra premise: its supports fall within Cases2,3,or4. The present conclusion remains eventual even when some special supports are independently known at every order.

## 4. Above-critical supports: two local neighborhoods and one compact complement

Fix a,d in Case4 and put D=a+d. The actual stationary map is

    Omega(S)=S^d*(D-S)/(S-a).

Its two nontrivial critical points are a+(d-1)/2 plus or minus i*sqrt(a-(d-1)^2/4). Their critical values form a nonreal conjugate pair w_+,w_-, with w_+ in the upper half-plane. They are distinct and nonzero.

The exact analytic inputs are:

- ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md gives, for every compact K subset of C minus {w_+,w_-}, an unevaluated terminal cutoff beyond which no actual common root has x/N^d in K. It covers every regular upper and lower parameter, Stokes walls, and the complete real boundary including0. Equal-height Stokes neighborhoods use local Morse exit loss; no false phase-separation assertion is imposed there. At a regular equimodular pair, the no-vertical-pair lemma supplies the nonzero real critical separation needed by the actual radial comparison.

- SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md, Theorem FixedSimpleTurnNeighborhoodConsecutiveExclusion, supplies a fixed physical disc about w_+ and an unevaluated cutoff excluding every actual common root in that disc. Conjugation supplies the analogous disc about w_-. This is a full fixed neighborhood, including unbounded scaled fields. Its bounded-field Airy expansion and growing-field exact individual-thimble ratio overlap. The individual thimble periods are interpolated analytically in the order only for the needed ratio estimates; the original positive-cycle identity is used at its actual integer orders.

The d=2 contour issue is included in both inputs. A locally smaller physical wedge has lower angle strictly less than pi/2, contains the selected critical points and zero, has inward normalized flow, and has genuine decaying infinity caps. The nondecaying closed boundary of the full angle pi/2 is not used.

Only the retained OLD upper theorem thm:product-annulus from GAMMA_PRODUCT_ROOT_MODULUS.tex is needed for root capture. Apply it with n=N and gamma=aN+1. Since A_N(x)=x^N*T_(N,aN+1,d)(1/x) and A_N(0)!=0, every root satisfies

    |x|<=((a+d)*N-d+1)_d<=((a+d)*N)^d,
    |x/N^d|<=R_(a,d):=(a+d)^d.                         (M6)

The last factor of the rising factorial is(a+d)N, so the second bound is termwise. Its shape hypothesis gamma>0 is automatic. No new sharp lower-annulus theorem is a premise.

Choose radii r_+,r_->0 supplied by the two local exclusions and reduce them if necessary so the discs are disjoint. Form the FIXED compact set

    K={w:|w|<=R_(a,d), |w-w_+|>=r_+/2,
                           |w-w_-|>=r_-/2}.            (M7)

It excludes both turning values. Apply the compact-complement theorem to K and take the maximum of its cutoff, the two local cutoffs and2. If a common root existed beyond that maximum, its normalized w would lie in the root-containing disk by(M6). A point in either physical turning disc is excluded locally; every remaining point belongs to K and is excluded by the compact complement. These alternatives cover the disk, with harmless overlap.

The compact-complement cutoff itself comes from a finite cover by smaller controlled parameter neighborhoods with collars for the adjacent radial transport. Its maximum is chosen before varying x or A/B. Thus the resulting N_(a,b) is uniform for the full coefficient plane at this fixed support. This proves(M2) in Case4 and completes the support exhaustion.

## 5. Descent from complex analysis to rational coprimality and every field

At each terminal order above the appropriate support cutoff, the analytic proofs or the all-order antecedents exclude a common COMPLEX root. Since the actual polynomial coefficients are rational, this is equivalent to coprimality in Q[x]. The Euclidean algorithm supplies U_N,V_N in Q[x] with

    U_N*A_(N;a,d)+V_N*A_(N-1;a,d)=1.                    (M8)

Every characteristic-zero field contains its canonical injective copy of Q. Transport(M8) to K and evaluate at A/B when B!=0. Both actual polynomial values cannot vanish, and the scalar factors in(M4) are nonzero. The monomial cases were already handled. This proves(M3) with the SAME support cutoff over all fields and coefficients.

No embedding of an arbitrary characteristic-zero field into C, no norm argument, and no complex conjugation on K are used. Complex conjugation was used only for the analytic polynomial family over Q before the rational identity was obtained.

Conversely, if a rational actual pair had a nonconstant gcd, that gcd would have a root z in an algebraic closure of Q. Over that characteristic-zero field the nonzero polynomial X^a*(z+X^d) would have both moments zero by(M4). This proves the reverse from a universal field-and-coefficient statement to rational coprimality. A statement about roots only in one specified non-algebraically-closed field would not suffice.

## 6. Genuine reverse equivalence: finite bad terminal orders

For each fixed integer support0<=a<b define d=b-a and

    BadTerminal_(a,b)
       ={N>=2 integer:
          gcd_(Q[x])(A_(N;a,d),A_(N-1;a,d))!=1}.          (M9)

**Corollary OrdinaryBinomialFiniteBadTerminalEquivalence.**
For every fixed support, the following are equivalent:

(i) BadTerminal_(a,b) is finite.

(ii) The rational actual consecutive pair is coprime for every sufficiently large terminal N.

(iii) A single cutoff works over every characteristic-zero field and every nonzero coefficient pair of that support, beyond which consecutive factorial moments do not both vanish.

(iv) A single cutoff works over C for every nonzero coefficient pair of that support.

Proof. A finite bad set permits cutoff2 if empty, or1+max of that set otherwise. Conversely(ii) confines it to the finite interval below its cutoff. The rational Bezout identity proves(ii) implies(iii), which implies(iv). If(iv) held while a later rational gcd were nonconstant, choose any complex root z of it and use f=X^a*(z+X^d), contradicting(M4). Thus(iv) implies(ii). The algebraic-closure argument of Section5 proves the same reverse without a chosen complex embedding.

The theorem establishes all four equivalent statements at EVERY fixed ordinary binomial support. It does not establish BadTerminal_(a,b) empty, enumerate any remaining finite bad set, or close all-order Strong Factorial for the supports not already known by separate theorems.

## 7. Scope of the conclusion and of the evidence

The result is unrestricted in the FIXED pair of ordinary monomial exponents, and uniform over coefficients and characteristic-zero fields once that support is fixed. The cutoff can vary with the support, is not evaluated, and is not proved effective. There is no assertion for arbitrary polynomials with three or more monomials, for multivariate supports, or for general Gamma shapes.

The written analytic inputs provide complete cycles, end bounds, Gaussian errors, shrinking-neighborhood estimates and actual order transport. Source guards freshly authenticate those arguments and their exact finite controls. Such finite controls are supporting evidence, not formal verification or a finite-prefix inference of eventuality. The old Gamma upper bound is sufficient for this proof; the newer sharp lower-annulus packet is an independent result.

