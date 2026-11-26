# Eventual Strong Factorial for square-root separated ordinary supports

For integers 1<=a<b with sqrt(b)-sqrt(a)>=1, every sufficiently large consecutive pair of factorial moments of a nonzero polynomial A X^a+B X^b cannot vanish simultaneously. The terminal cutoff depends on the fixed support (a,b), and is uniform over the coefficients and all characteristic-zero fields.

Set d=b-a and define, for every integer m>=1,

    A_(m;a,d)(x)
      =sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k).      (Q1)

Here (t)_j=t(t+1)...(t+j-1), with (t)_0=1. These are monic integer polynomials of degree m, with nonzero constant terms. For a characteristic-zero field K, let Lambda_K:K[X]->K be the K-linear functional Lambda_K(X^j)=j!.

## 1. Main theorem and support condition

**Theorem SquareRootGapEventualStrongFactorial.** For every fixed pair of integers 1<=a<b satisfying

    sqrt(b)-sqrt(a)>=1,                                 (Q2)

there exists an integer N_(a,b)>=2 such that

    gcd_(Q[x])(A_(N;a,d),A_(N-1;a,d))=1                  (Q3)

for every integer N>=N_(a,b), where d=b-a. The same cutoff works for every characteristic-zero field K and every coefficient pair A,B in K with (A,B)!=(0,0): if f=A X^a+B X^b, the two values Lambda_K(f^N) and Lambda_K(f^(N-1)) are not both zero for N>=N_(a,b).

This is an eventual theorem. No effective value of N_(a,b), cutoff uniform over different supports, or conclusion at every smaller terminal order is asserted.

**Lemma SquareRootGapSupportEquivalence.** For integers a>=1,d>=1, condition(Q2) with b=a+d is equivalent to

    4a<=(d-1)^2 and d>=3.                               (Q4)

The equality boundary consists exactly of the consecutive-square supports

    a=s^2, d=2s+1, b=(s+1)^2, s>=1 integer.              (Q5)

Proof. The inequality sqrt(a+d)>=sqrt(a)+1 is equivalent to d>=1+2sqrt(a). Both sides of d-1>=2sqrt(a) are nonnegative, so squaring is reversible. Since a>=1 this forces d>=3. Conversely(Q4) and d>=3 permit taking nonnegative square roots and reversing the same steps. If equality holds, d-1=2sqrt(a); an integer with rational square root has integer square root, so sqrt(a)=s is a positive integer and(Q5) follows. Conversely(Q5) gives equality. In particular the strict and equality regimes in Section2 exhaust every support of the main theorem.

## 2. Two independently proved support regimes

The synthesis uses the following written analytic results.

1. **Strict subcritical regime.** SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md, Theorem StrictSubcriticalEventualConsecutiveCoprimality, proves(Q3) for every fixed integer a>=1,d>=3 with

       4a<(d-1)^2.

   Its proof verifies the actual two distinguished stationary branches in the physical lower wedge, the global real-part divider, the positive-real coexistence of two nondegenerate maxima, and the strict inactivity of each spinodal inflection. It revalidates the complete Morse and finite-sublevel Stokes contours for these parameters, rather than applying a theorem stated only at a cusp. Real boundary neighborhoods, the exact own-order consecutive comparison, and the retained global upper modulus yield a compact cover of the entire normalized coefficient plane.

2. **Critical equality regime.** CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md, Theorem EventualIntegerCuspStrongFactorial, proves(Q3) for every fixed integer cusp(Q5). Its proof combines a fixed physical neighborhood of the quartic cusp with complete compact-complement contours, including Stokes walls and every other real boundary point. It then uses the retained global upper modulus to capture all possible common roots. The local bounded-field theorem is only a core within that physical neighborhood, and no assumption that all roots remain in a bounded scaled cusp field is made.

Both results concern the actual polynomials(Q1), whose Gamma shapes depend on their own orders. Both give a threshold for the whole coefficient plane at the fixed support. Neither conclusion is a statement about a chosen coefficient ratio alone.

For clarity, the exact analytic normalization in both inputs is

    I_m(w)=integral_0^infinity [S^a*(w+S^d)]^m*exp(-m*S)dS
          =(a*m)!*A_(m;a,d)(m^d*w)/m^((a+d)*m+1).        (Q6)

Thus a common actual root x at orders N and n=N-1 gives

    I_N(x/N^d)=I_n(x/n^d)=0.

The lower parameter is (N/n)^d times the upper parameter. This own-order transport is part of both analytic proofs.

Only the retained OLD upper theorem thm:product-annulus in GAMMA_PRODUCT_ROOT_MODULUS.tex is used to contain the roots. Applied with gamma=a*m+1, it gives

    |x|<=((a+d)*m-d+1)_d<=((a+d)*m)^d,
    |x/m^d|<=(a+d)^d.                                  (Q7)

The reciprocal-polynomial relation is A_m(x)=x^m T_(m,am+1,d)(1/x), and the constant term of A_m is nonzero. Hence no zero is lost under reversal. In particular no new sharp lower-annulus or shape-separation result is a logical premise of this synthesis.

## 3. Exhaustion of the allowed supports

Fix a,b satisfying(Q2). Lemma SquareRootGapSupportEquivalence gives integers a>=1,d>=3 satisfying(Q4).

If the inequality is strict, Input1 gives an integer threshold N_(a,b) for which(Q3) holds for every later N. If equality holds, the lemma identifies the unique integer s of(Q5), and Input2 gives the threshold N_s; take N_(a,b)=N_s.

These alternatives are exhaustive. No limit in the support parameters is taken, and no uniform threshold is inferred by combining them. At each fixed support the threshold is chosen before the field or coefficient pair is varied.

To state why these are coefficient-uniform conclusions, both inputs first exclude actual common complex roots on controlled normalized neighborhoods. They shrink each neighborhood to allow the adjacent own-order transport, then take a finite cover of a fixed compact disk supplied by(Q7). In the equality case the fixed cusp neighborhood is combined with a compact complement before that final cover. Each support therefore has one threshold for every captured root, hence every coefficient ratio.

This proves(Q3). The remaining arithmetic step transfers that rational polynomial statement to the original factorial functional.

## 4. Rational Bezout descent and arbitrary characteristic-zero fields

For B!=0, expansion of the original binomial gives exactly

    Lambda_K((X^a*(A+B X^d))^m)
       =sum_(k=0)^m binom(m,k)*A^(m-k)*B^k*(a*m+d*k)!
       =B^m*(a*m)!*A_(m;a,d)(A/B).                      (Q8)

All factorials are nonzero in a characteristic-zero field. If B=0 and A!=0, every moment is the nonzero value A^m*(a*m)!. If A=0 and B!=0, the same is true for the single monomial B X^b, or follows from(Q8).

For each N>=N_(a,b), the Euclidean algorithm applied to(Q3) supplies U_N,V_N in Q[x] such that

    U_N*A_(N;a,d)+V_N*A_(N-1;a,d)=1.                    (Q9)

Every characteristic-zero field has a canonical injective map from Q. Transport(Q9) along this map and evaluate at A/B when B!=0. Both actual polynomial values cannot vanish. Their nonzero scalar factors in(Q8) establish the main moment assertion. Thus the same threshold works for every field and every nonzero coefficient pair. There is no assumption that an arbitrary characteristic-zero field embeds in C.

Conversely, a common nonconstant rational divisor would have a root z in an algebraic closure of Q. Over that characteristic-zero field, the nonzero polynomial X^a*(z+X^d) would have both moments zero by(Q8). This is the correct reverse implication from a universal field-and-coefficient assertion to actual rational coprimality; truth only over one specified non-algebraically-closed field would not give it.

## 5. Genuine eventual equivalence and finite bad terminals

For each fixed integer support a>=1,d>=1 define

    BadTerminal_(a,d)
       ={N>=2 integer:
          gcd_(Q[x])(A_(N;a,d),A_(N-1;a,d))!=1}.          (Q10)

**Corollary FiniteBadTerminalEventualEquivalence.** At every fixed integer support, the following are equivalent:

(i) BadTerminal_(a,d) is finite.

(ii) There is an N_0>=2 such that the rational consecutive polynomials are coprime for every N>=N_0.

(iii) There is an N_0>=2 such that, over every characteristic-zero field and for every nonzero coefficient pair, the original binomial's consecutive moments do not both vanish at terminal orders N>=N_0.

Proof. An empty bad set permits N_0=2; otherwise a finite bad set permits N_0=1+max BadTerminal_(a,d). Conversely(ii) confines that set to {2,...,N_0-1}. The rational Bezout argument proves(ii) implies(iii). The algebraic-closure construction at the end of Section4 proves(iii) implies(ii). No analytic estimate or checker is being reversed.

The theorem proves these equivalent statements for every support satisfying(Q2). It therefore proves finiteness of the full set of bad terminal orders at each such support. It does not identify that finite set, prove that it is empty, or resolve the all-order Strong Factorial assertion at those supports.

## 6. Scope and evidence

The conclusion covers both consecutive-square supports and every strict square-root separated support. It does not assert an eventual theorem for the remaining supports with sqrt(b)-sqrt(a)<1. No uniform cutoff over supports or numerical value for a cutoff is obtained.

The source guards authenticate and freshly replay the analytic input packets. Their finite algebraic controls support the written contour proofs, exact own-order normalizations, support equivalence and field descent. They are not formal verification and do not infer eventuality from a finite list of moment or modular-gcd tests. The logical root-capture premise is the retained old Gamma upper bound; the independent newer sharp annulus is not used.

