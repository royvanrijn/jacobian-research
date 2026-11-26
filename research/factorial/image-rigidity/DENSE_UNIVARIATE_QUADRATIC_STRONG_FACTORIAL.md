# Dense univariate quadratic Strong Factorial windows

For every characteristic-zero field, a quadratic with three nonzero coefficients has a nonzero factorial moment in every positive block of length three. The proof uses a degree-bounded Hensel factorization in a divided-power lattice. It does not settle arbitrary binary quadratics or full SFC(2).

Let Z_2 denote the 2-adic integers. In Q_2[Z] put

    eta_k(Z)=2^k Z^k/k!,  k>=0.

The finite polynomials sum b_k eta_k with b_k in Z_2 form a ring H, since

    eta_i eta_j=binom(i+j,i) eta_(i+j).                 (1)

Reduction modulo2 below is reduction of the coefficients b_k in THIS lattice, not ordinary reduction of the polynomial's coefficients in Z_2[Z]. The distinction matters: v_2(2^k/k!)=s_2(k).

## 1. Highest-bit Hensel factorization

**Lemma 1 (binary divided-power splitting).** Let m be positive and even, and let

    Q(Z)=sum_(k=0)^m b_k eta_k(Z),

where every even-index b_k is a 2-adic unit and every odd-index b_k is even. Then Q factors over Q_2, up to a constant unit, into factors Q_e of distinct degrees e, one for each nonzero binary digit of m. Each factor has divided-power integral coefficients and even odd-index coefficients. Its constant coefficient is a unit and its leading divided-power coefficient is a unit. In ordinary coefficients,

    v_2([Z^k]Q_e)>=s_2(k)+1_(k odd),  0<=k<=e,         (2)

with a unit constant and leading ordinary coefficient of valuation1. Its normalized reciprocal is therefore an Eisenstein polynomial of degree e.

Here is the complete Hensel argument. Let e be the highest power of2 in m and write m=e+S, 0<=S<e. In the divided-power algebra modulo2 let A_e be the subspace spanned by eta_0,...,eta_(e-1). It is a subalgebra: if i,j<e and i+j>=e, binary addition has a carry, hence binom(i+j,i) is even. Also, for j<e,

    eta_e eta_j=eta_(e+j) mod2.                        (3)

Indeed the binary digits of e and j are disjoint. These parity facts follow directly from (1+T)^n=product_(r:n_r=1)(1+T^(2^r)) over F_2, so no binomial-unit hypothesis at unrelated indices is being assumed.

The augmentation ideal of A_e is nilpotent. More explicitly every element with zero constant has square zero: each eta_i^2 has even central binomial coefficient for i>0, and cross terms cancel in characteristic2. Thus an element of A_e with constant1 is invertible.

Write the reduction of Q uniquely as

    Qbar=Qlo+eta_e Qhi,
    Qlo in A_e, deg Qhi<=S.                            (4)

Coordinates of Qhi are exactly b_(e+j) modulo2, by(3). In particular Qhi has constant1, has a unit coefficient at degree S, has coefficient1 at every even j<=S and0 at every odd j.

We seek factors G of degree e and U of degree S, normalizing the leading divided-power coefficient of G to be exactly1. Modulo2 define

    Ubar=Qhi,
    Gbar=Glo+eta_e,  Glo=Qlo/Ubar in A_e.               (5)

All products of two degree<e terms crossing e vanish, so(5) really gives Gbar Ubar=Qbar. Both constants are1. Since Qlo and Ubar have only even indices, their quotient in A_e has only even indices as well: the even subalgebra is closed and its unit inverse has only even indices. Hence Gbar also has no odd-index term. The coefficient at eta_e is1 by construction.

Choose any integral divided-power coefficient lifts G_1,U_1 of(5) with the stated degree bounds and G_1's leading coefficient exactly1. Suppose inductively that

    Q-G_n U_n belongs to 2^n H,  n>=1.

The error E=(Q-G_nU_n)/2^n has degree at most m. Decompose its reduction as Elo+eta_e Ehi with deg Ehi<=S. Seek corrections D of degree<e and V of degree<=S. Linearizing modulo2 gives

    D Ubar+Gbar V=Ebar.

Its high-degree part is V exactly: D Ubar has no high part, and only eta_e V contributes to the high part of Gbar V. Thus the unique corrections are

    V=Ehi,
    D=(Elo-Glo V)/Ubar in A_e.                         (6)

In particular deg D<e and deg V<=S, with integral divided-power representatives. Set

    G_(n+1)=G_n+2^n D,  U_(n+1)=U_n+2^n V.

The product is correct modulo2^(n+1); its quadratic correction has factor2^(2n), which is divisible by2^(n+1) for n>=1. The leading divided-power coefficient of G stays exactly1. The limits exist coefficientwise in Z_2, in FIXED finite degrees e and S, so they are polynomials, and Q=GU exactly. This is not an appeal to normality of H or an unbounded formal-series factorization.

The mod2 coefficients remain those in(5). Consequently both factors have unit constants, integral divided-power coefficients, and even odd-index coefficients. The factor U has a unit at EVERY even index and even coefficients at odd indices, because Ubar=Qhi. Its degree S is even and smaller than m, unless S=0, in which case U is a constant unit. Repeat this same construction on U. This proves the asserted factorization by the binary degrees of m.

For any resulting factor Q_e, equation(2) follows from v_2(2^k/k!)=s_2(k), with the extra1 when its divided-power coefficient is even. Its top index e is a power of2, so the top ordinary coefficient has valuation exactly1. All other positive-index ordinary coefficients are divisible by2, and the constant is a unit. Hence

    g_e(B)=B^e Q_e(1/B)/Q_e(0)                         (7)

is monic Eisenstein of degree e. No ordinary coefficient division by a nonunit has been hidden in the Hensel construction.

## 2. Quadratic factorial moments meet the exact hypotheses

For an odd integer C and an even integer m>=2 put

    F_m(A;C)=Lambda((Ax^2+x+C)^m),  Lambda(x^a)=a!,
    s=s_2(m),  P_m(B)=2^s F_m(B/4;C).

Selecting the quadratic factors gives

    [A^i]F_m=binom(m,i)(2i)! C^(m-i) H_i,
    H_i=sum_(l=0)^(m-i) binom(m-i,l)(2i+1)_l C^(-l).     (8)

Since C is odd, H_i is 2-adically integral. For l>=2 the rising product is even; the l=1 term gives

    H_i=1+(m-i) mod2=1+i mod2.

Thus H_i is a unit when i is even, and is even when i is odd. Direct factorial valuations give

    v_2(binom(m,i)(2i)!)=2i+s_2(m-i)-s_2(m).

It follows that P_m has 2-adic integral coefficients, unit leading coefficient, and

    v_2([B^i]P_m)>=s_2(m-i)+1_(i odd),                 (9)

with equality at every even i. Divide by its leading unit to make P_m monic; this does not change any valuation. Its reciprocal Q_m(Z)=Z^m P_m(1/Z) has ordinary coefficients u_k with

    v_2(u_k)>=s_2(k)+1_(k odd),

and equality at every even k, because m is even. Write u_k=b_k 2^k/k!. Then every b_even is a unit and every b_odd is even. Lemma1 applies.

Consequently P_m has one monic Eisenstein factor g_e of degree e for every nonzero bit e=2^r of m. All roots of g_e have valuation1/e; the degrees are distinct, so these are exactly the separate Newton factors. If

    g_e(B)=B^e+sum_(j=0)^(e-1) a_j B^j,

then the reciprocal coefficient bound(2), unchanged by constant-unit normalization, becomes

    v_2(a_j)>=s_2(e-j)+1_(j odd).                     (10)

This proves the desired inheritance of the full divided-power bound and its odd refinement; it is not inferred from finitely many factorizations.

## 3. The common-pair derivative forces the forbidden coefficient

**Theorem 2 (odd-constant upper-even pair).** For every odd integer C and every even m>=2, the rational polynomials F_(m-1)(A;C) and F_m(A;C) are coprime.

We first derive the exact needed derivative identity. Let M_j=F_j(A;C). Write integrals against e^(-x) over the positive real axis only to prove polynomial identities; each identity then holds algebraically over Q in A,C. For m>=1,

    (2A x+1)^2=4A(Ax^2+x+C)+(1-4AC).

Integrating derivatives of x f^m and x f' f^(m-1), or equivalently comparing their native factorial coefficients, yields

    2A^2 partial_A F_m
      =(2Am-1)F_m+m(1-2AC)F_(m-1)+C^m.                (11)

For a direct verification, write I=Lambda(x f^(m-1)) and D=partial_A F_m. Integration of(f^m e^(-x)) gives F_m=m(2A I+F_(m-1))+C^m. The identity f=Ax^2+x+C also gives I=F_m-A D/m-C F_(m-1). Substitution yields(11) exactly. Its boundary term has sign plus.

Suppose the two polynomials had a common root A in an algebraic closure of Q_2. Since F_m(0;C) is odd by(8), A is nonzero. Put beta=4A and let e=2^r be its Eisenstein factor degree in Section2. Then v_2(beta)=1/e and (11) at the common root gives

    v_2(partial_A F_m)=3-2/e,
    v_2(P_m'(beta))=s+1-2/e.                          (12)

Here multiplying P_m by a leading unit changes no valuation. Let S=m mod e, ell=s_2(S), and h=s-ell-1, the number of nonzero bits of m larger than e. In the product of all the other Eisenstein factors at beta, a factor of degree f<e contributes valuation f/e: its roots have strictly larger valuation1/f. A factor of degree f>e contributes valuation1: each root has smaller valuation1/f. Thus

    v_2(product_(f!=e) g_f(beta))=S/e+h.

Using the simple-root product derivative therefore gives

    v_2(g_e'(beta))=ell+2-(S+2)/e.                    (13)

Each nonzero derivative summand j a_j beta^(j-1), 1<=j<=e, with a_e=1, has valuation

    v_2(j a_j)+(j-1)/e.

Their fractional parts are distinct, since j-1 runs from0 through e-1. Thus there is a unique least term and no cancellation among derivative terms. The fractional part forced by(13) identifies its index uniquely as

    j=e-S-1.

Indeed m and S are even, so0<=S<=e-2; this index is odd and lies in1..e-1, including e=2. At this j, equation(13) forces

    v_2(a_j)=ell+1.                                   (14)

But (10) gives

    v_2(a_j)>=s_2(S+1)+1=ell+2,

because S is even. This contradiction proves Theorem2. The argument uses every binary face, including the lowest three: no prior exclusion of initial bit clusters is needed.
## 4. Even constants: a separate degree-preserved argument

For completeness, when C is an even integer put

    U_j(T)=Lambda((T x^2/2+x+C)^j)/j!.

Formal factorial expansion gives

    sum_(j>=0) U_j(T) z^j
      =e^(Cz) sum_(i>=0) (2i-1)!! T^i z^i/(1-z)^(2i+1). (15)

Every double factorial is odd. For every r>=1, v_2(C^r/r!)>=s_2(r)>=1, interpreting a zero term as infinite valuation. Thus U_j lies in Z_(2)[T], has odd leading coefficient (2j-1)!! and degree j. Its reduction p_j over F_2 has generating series and recurrence

    sum p_j z^j=(1+z)/(1+Tz+z^2),
    p_0=1, p_1=T+1, p_j=T p_(j-1)+p_(j-2).             (16)

A common root of p_j,p_(j+1) propagates backward to p_0=0, impossible. Their resultant is nonzero over F_2. Both degrees are preserved, so the rational resultant of U_j,U_(j+1) is a 2-adic unit, in particular nonzero. A rational Bezout identity persists in every characteristic-zero field. We have proved consecutive coprimality for every j>=1 when C is even. No reduction of a possible coefficient-field root modulo2 is being assumed.

## 5. The dense univariate quadratic conclusion

**Theorem 3 (dense univariate quadratic Strong Factorial).** For every characteristic-zero field K, every f=ax^2+bx+c with abc!=0, and every n>=1, at least one of Lambda(f^n),Lambda(f^(n+1)),Lambda(f^(n+2)) is nonzero.

Scale by b and put A=a/b,C=c/b, so A,C are nonzero. For F_j=Lambda((Ax^2+x+C)^j), the exact factorial recurrence is

    F_j=2Aj(2j-1)F_(j-1)
          +j(j-1)(1-4AC)F_(j-2)+C^j+jC^(j-1), j>=2.  (17)

Here is a direct polynomial-identity derivation. With f=Ax^2+x+C and M_r=Lambda(f^r), let J_r=Lambda(f' f^r). Integration of the derivative of f' f^r e^(-x), with boundary -C^r, and (f')^2=4Af+(1-4AC), gives

    J_r=2A(2r+1)M_r+r(1-4AC)M_(r-1)+C^r.

Integration of the derivative of f^(r+1)e^(-x) gives M_(r+1)=(r+1)J_r+C^(r+1), yielding(17). All expressions are polynomial identities over Q; proving them over real parameters suffices for all characteristic-zero fields.

If a three-zero block starts at n, equation(17) at j=n+2 implies

    C^(n+1)(C+n+2)=0.

Since C is nonzero, C=-q with q=n+2>=3 an integer. If q is even, Section4 excludes even the first consecutive zero pair. If q is odd, m=q-1=n+1 is even and C is odd; Theorem2 excludes F_(m-1)=F_m=0. The relevant coprimeness statements are over Q and therefore hold over K by rational Bezout identities. This contradiction proves Theorem3 without assuming that K embeds in C.

## 6. Scope, provenance and the remaining dimensional problem

Theorem3 concerns exactly the three nonzero coefficients of a univariate quadratic. It does not prove the six-monomial inhomogeneous binary quadratic problem, full SFC(2) (where2 counts variables), or SFC(3). Other univariate supports may be added by separate exact arguments; they are not silently included in this dense statement.

The divided-power splitting proof is local and self-contained. Its leading-unit normalization is in the divided-power lattice, not ordinary polynomial coefficient reduction. Only products crossing the highest-bit boundary vanish modulo2; products within the positive ideal need not vanish. Each individual positive element has square zero, which suffices for the inverse in(5)--(6).

The monotone prime-unit native transfer is separate and is not a premise. No finite factorization diagnostic, numerical resultant scan, odd-prime degree sieve, bounded moment prefix or irreducibility black box is a proof input. No literature-priority or formal-verification claim is made.

## 7. Supplementary native controls

The [portable standard-library checker](verify_dense_quadratic_divided_power_hensel.py) and [retained finite receipt](dense_quadratic_divided_power_hensel_verification.json) stress the exact divided-power multiplication and finite-precision degree-bounded Hensel updates on native and deterministic synthetic inputs. These checks are supplementary; the universal result is the all-range argument in Lemma1 and Theorems2--3. From this note's directory, run a fresh verification with ordinary Python:

```text
python verify_dense_quadratic_divided_power_hensel.py --output fresh_dense_quadratic_verification.json
```

An existing output file is never overwritten and optimized Python is refused. The checker does not infer a universal theorem from tested starts or require PARI, a network connection or external inputs.
