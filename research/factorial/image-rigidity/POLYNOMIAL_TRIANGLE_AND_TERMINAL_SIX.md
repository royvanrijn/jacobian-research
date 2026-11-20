# Polynomial support triangles and the sixth ordinary binomial terminal

20 November 2026. Internal mathematical proof with exact finite-field replay. No unrestricted Strong Factorial theorem is asserted.

Put L(X^j)=j! and, for integers a>=0,d>=1,m>=1, define the ACTUAL normalized moment polynomial

    P_(m,a,d)(t)=sum_(k=0)^m binom(m,k)(am+1)_(dk) t^k.

For f=X^a(A+B X^d) and A!=0,

    L(f^m)=A^m (am)! P_(m,a,d)(B/A).                  (T1)

The two consecutive shapes are a(N-1)+1 and aN+1. They are not silently replaced by a common shape. Coefficient-zero monomial phases have every positive moment nonzero over a characteristic-zero field.

**Theorem OrdinaryBinomialTerminalSix.** Over every characteristic-zero field, for every a>=0,d>=1 and every nonzero ordinary binomial X^a(A+B X^d), the moments of orders5 and6 cannot both vanish. Together with the retained initial-order results, every terminal N=2,...,6 is therefore closed at every pair of exponents.

For N>=2 let e_N be the least integer satisfying 2^(e_N)>=4N(N-1), and put

    D_N=2(N-1)e_N,
    Triangle_N={(a,d):2<=d<D_N,1<=a<d(2N-3)}.          (T2)

The definition is exact integer arithmetic; e_N=(4N(N-1)-1).bit_length(), with no floating logarithm.

**Theorem PolynomialTriangleReduction.** The following are equivalent:

1. Every nonzero ordinary complex binomial X^a(A+B X^d), a>=0,d>=1, has a nonzero moment in every consecutive positive pair.
2. For EVERY integer N>=7 and EVERY (a,d) in Triangle_N, the two rational polynomials P_(N-1,a,d),P_(N,a,d) are coprime.
3. The full assertion in1 holds over EVERY characteristic-zero coefficient field.

Only a finite exponent triangle occurs at each fixed N, but N remains unbounded. In particular this theorem and the terminal-six computation do not prove1. The coprimality tests in2 are not claimed equivalent to phase validity over one prescribed non-algebraically-closed field.

**Corollary PolynomialTriangleSize.** For fixed N, the maximal gap is D_N-1=O(N logN), the maximal lower exponent is (D_N-1)(2N-3)-1=O(N^2 logN), and the NUMBER of exponent pairs is

    (2N-3)*[(D_N-1)D_N/2-1]-(D_N-2)
       =O(N^3 log^2N).                                (T3)

Thus this replaces the older exponentially large lower-exponent rectangle by a polynomially bounded support triangle. The pair count itself is not asserted to be O(N^2).

## 1. Two retained all-exponent complement theorems

The actual Theorem thm:gapbound of lower_exponent_initial_and_boundary.tex proves that for every N>=2, a>=1 and d>=D_N, all first N ACTUAL moment polynomials M_(j,a,d) are pairwise coprime. Its argument uses the same integer e_N as(T2), the uniform moment-root annuli and monotone Gamma centers. Multiplying the polynomials by nonzero rational factorial scalars does not change their common roots, so it applies to the P_j above.

The actual Corollary cor:linear-lower-exponent of GAMMA_PRODUCT_ROOT_MODULUS.tex proves, for every N>=2,d>=1 and positive Gamma shape beta, that a>=d(2N-3) makes T_(j,aj+beta,d),1<=j<=N, pairwise coprime. At ordinary beta=1 these are exactly P_j. Equality at the threshold is allowed, so the residual triangle has STRICT upper bound a<d(2N-3).

The Gamma-product proof uses its exact-degree circular-domain Grace--Walsh--Szego input, including exterior disks. The cached Wagner primary Theorem4.1, PDFpage11, supports exactly that retained premise. No half-plane-only polarization statement is substituted, and no general claim of simplicity from polarization is made here.

For every fixed N, a potential ordinary consecutive counterexample with a>=1,d>=2 is therefore inside Triangle_N: d>=D_N is excluded by the gap theorem; at d<D_N, a>=d(2N-3) is excluded by the linear theorem. This is a direct logical combination of two all-order complement results, not a new finite moment scan or an asymptotic extrapolation.

## 2. Exact completion of the sixth terminal

At N=6 the definition gives e_6=7,D_6=70 and 2N-3=9. Thus the complete residual set is

    2<=d<70, 1<=a<9d,
    sum_(d=2)^69(9d-1)=21658.                        (T4)

Take the prime p=1000003. Exact trial division through1000 proves primality. Every actual factorial argument needed for m=5 or6 is at most

    max_(Triangle_6)(a+d)*6=(620+69)*6=4134<p.

All normalization factors (am)! are units modulo p, both P_5 and P_6 retain their full degrees5 and6, and their constant coefficients equal1. The primary exact checker computes the full gcd over F_p for EVERY one of the21658 pairs in(T4). Every gcd is1. Its complete case stream agrees with the separately obtained original full replay; the enumeration is not a sample.

At each fixed(a,d), the integer resultant Res(P_5,P_6) therefore reduces to a nonzero element modulo p. Since both full degrees are preserved, it is the actual specialized resultant, proving that this integer is nonzero and that the rational polynomials are coprime. A rational Bezout identity transfers the same fact to every characteristic-zero coefficient field. Equation(T1) proves ordinary consecutive nonvanishing throughout(T4).

The complement theorems in Section1 cover all a>=1,d>=2 outside(T4). For a=0, the retained ALL_UNIVARIATE_TWO_TERM_SFC.tex theorem covers every gap and terminal. For d=1, the retained UNIFORM_PADE_PARABOLA_SFC.tex adjacent-exponent theorem covers every lower exponent and terminal. Hence all exponents at N=6 are covered, proving OrdinaryBinomialTerminalSix. The retained Corollary cor:initial and Theorem thm:ordinary-fifth of lower_exponent_initial_and_boundary.tex already cover N<=5.

A separate complete internal replay uses the different prime1000033 and builds coefficients by binomial/sliding-rising recurrences without the primary factorial tables. It checks all21658 full degree5/6 polynomial Bezout identities and the same domain. Its frozen source and receipt are additional independent finite evidence. Each finite-prime success implies rational coprimality; no reverse claim from rational coprimality to success at a prescribed prime is used.

## 3. Both directions of the remaining all-order reduction

Assertion1 of PolynomialTriangleReduction implies2. If a displayed rational pair were not coprime, it would have a common complex root theta. Since both constant coefficients are1, theta!=0. Taking A=1,B=theta gives a nonzero complex phase with two zero moments, contradicting1.

Conversely assume2. Terminals N<=5 have the retained initial proofs and N=6 has Section2. At each N>=7, coefficient-zero phases are monomials, a=0 and d=1 have their retained all-order proofs, and the two complement theorems cover all positive a,d outside Triangle_N. Inside Triangle_N the actual pair is coprime by2. Equation(T1) proves1 in every case.

If1 holds, every fixed rational pair is coprime; otherwise its common complex root yields the same contradiction. Each pair consequently has a rational Bezout identity, which proves3 over arbitrary characteristic-zero fields without an embedding assumption on that field. Assertion3 implies1 by taking the field C. This proves all three directions with the stated quantifiers.

No completed lower-one all-gap theorem is required by this argument: the residual domain deliberately includes a=1. That separately proved theorem can shrink subsequent searches to a>=2, but it is not an input to the present equivalence or terminal-six closure.

## 4. Integer size and assurance

By the integer cutoff definition in(T2), the least binary exponent satisfies e_N<=log2(4N(N-1))+1=O(logN). Hence D_N=O(NlogN). Summing d(2N-3)-1 over d=2,...,D_N-1 gives(T3) exactly. The lower-exponent maximum and pair-count orders follow immediately.

The checker verifies the entire finite sixth-terminal triangle with exact integer/modular arithmetic, literal actual-shape factorial controls, rational moment normalization and rejected domain/polynomial corruptions. Its integer-cutoff controls do not assert vanishing/nonvanishing at later terminal orders. The written complement arguments supply the universal exponent quantifiers; no all-N conclusion is inferred from terminal6. All assurance is internal proof review and exact replay, not formal verification or external review. No historical priority claim is made.

Proof support: [exact checker](verify_polynomial_triangle_terminal_six.py), [finite inputs](polynomial_triangle_terminal_six_inputs.json), [source ledger](polynomial_triangle_terminal_six_source_ledger.json), and [receipt](polynomial_triangle_terminal_six_verification.json). The separate [full Bezout replay](verify_polynomial_triangle_terminal_six_bezout.py) and [its receipt](polynomial_triangle_terminal_six_bezout_verification.json) provide an alternative exact finite-evidence route.
