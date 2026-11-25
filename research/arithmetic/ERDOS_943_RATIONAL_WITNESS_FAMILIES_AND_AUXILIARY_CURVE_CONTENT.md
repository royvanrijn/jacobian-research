# Rational witness families and auxiliary-curve content for square-full kernel pairs

We prove uniform counting bounds for an ACTUAL rational-function witness family. This does not assert that all compatible kernel pairs admit a small cover by such families, and does not improve the current prescribed-sum exponent 3/10.

Let T>=2 and c T^5<=N<=C T^5, with fixed positive c,C. The original problem counts ordered distinct squarefree kernels x,y in [T,2T), coprime to N, which admit a reduced positive unit fraction u/v of height at most 4T^2 with

\[
 x^3u^2\equiv y^3v^2\pmod N.
 \tag{1}
\]

A degree-d rational witness family is a pair of homogeneous forms U,V in Z[X,Y] of the SAME degree d, not both zero. A pair belongs to the family only if U(x,y),V(x,y) are nonzero of the same sign, and its actual reduced positive fraction is U(x,y)/V(x,y). The original height and unit conditions are still imposed. The forms may be divided by their common integer coefficient content without changing that ratio. Throughout, use this normalization, so the gcd of ALL coefficients of U and V together is one.

Put

\[
 F(X,Y)=X^3U(X,Y)^2-Y^3V(X,Y)^2,
 \qquad D=2d+3.
 \tag{2}
\]

Every actual family witness satisfies N|F(x,y), even if the pointwise gcd of U(x,y),V(x,y) contains primes dividing N: the reduced congruence actually implies divisibility after multiplication by that gcd squared. No inverse of that gcd is used.

**Theorem 1 (coefficient-uniform family bound).** For any such degree-d forms, even when their coefficients depend arbitrarily on N and T, the number of unit integer pairs x,y in [T,2T) with N|F(x,y) is

\[
 \ll_{d,\varepsilon}N^\varepsilon
     \left(T+\frac{T^2}{N^{1/(2d+3)}}\right).
 \tag{3}
\]

The constants do NOT depend on the coefficients of U,V. In particular any quadratic rational witness family contributes

\[
 \ll_{\varepsilon,c,C} T^{9/7+\varepsilon}
 \tag{4}
\]

to the original pair count. This upper bound does not require squarefreeness, distinctness, or the original witness-height restriction.

**Theorem 2 (fixed-family improvement).** If U,V are FIXED independently of N,T, then

\[
 \#\{x,y\in[T,2T): (xy,N)=1,\ N\mid F(x,y)\}
 \ll_{U,V,\varepsilon} N^\varepsilon
       \left(T+\frac{T^2}{N^{1/(2d+1)}}\right).
 \tag{5}
\]

Thus EVERY fixed quadratic family contributes O_(U,V,epsilon)(T^(1+epsilon)). This includes but is not limited to the particular degree-two Pade family. The coefficient dependence in (5) is essential to this proof; it must not be silently removed.

**Theorem 3 (prime modulus and actual squarefree pairs).** If N=p is prime and p>8T^2, there are at most 2d+3 reduced kernel ratios in ANY degree-d family, coefficient-uniformly. With the additional restriction gcd(x,y)<=G, the number of corresponding ordered kernel pairs is O_d(G).

If the family U,V is fixed, then its ACTUAL distinct squarefree pairs satisfy the stronger bound

\[
 \ll_{U,V} T\,p^{-1/(2d+3)}.
 \tag{6}
\]

For a fixed quadratic family at p comparable to T^5 this is O_(U,V)(T^(2/7)). Points where U=V=0 are invalid witnesses and are excluded from this assertion. The bounds for a family do not bound the union of an uncontrolled number of families.

## 1. Elementary prime-power root clusters

**Lemma 4.** Let f in Z[Z] have degree at most D>=1 and at least one coefficient coprime to a prime p. For e>=1 put k=ceil(e/D). The solutions of f(z)=0 modulo p^e lie in at most D residue classes modulo p^k.

**Proof.** If there were D+1 distinct such classes, choose integer representatives z_0,...,z_D that are roots modulo p^e. Each pair difference has p-valuation at most k-1. Lagrange interpolation over Q gives

\[
 f(Z)=\sum_{i=0}^{D} f(z_i)
       \frac{\prod_{j\ne i}(Z-z_j)}
            {\prod_{j\ne i}(z_i-z_j)}.
\]

Each denominator has p-valuation at most D(k-1)<e. Its numerator polynomial has integer coefficients, and f(z_i) is divisible by p^e. Every coefficient of every term on the right therefore has positive p-valuation, hence so does every coefficient of f. This contradicts the hypothesis. This proves the lemma for arbitrary repeated roots and arbitrary coefficient sizes, with no discriminant or Hensel hypothesis.

For a primitive integer polynomial f and arbitrary N, apply the lemma at each prime power p^e exactly dividing N. CRT covers all its roots modulo N by at most D^(omega(N)) classes modulo

\[
 Q_D(N)=\prod_{p^e\Vert N}p^{\lceil e/D\rceil}
       \ge N^{1/D}.
 \tag{7}
\]

Only membership in these coarser classes is asserted; not every element in one is a root modulo N. For fixed D, the number of classes is O_(D,epsilon)(N^epsilon): choose a fixed prime threshold beyond which D<=p^epsilon, and absorb the finitely many smaller primes in the constant. This is the elementary polynomial-root lemma used below. Its classical antecedent is described in Section 10.

## 2. Why the cusp polynomial is primitive at EVERY prime

For each prime p, the dehomogenizations U(Z,1),V(Z,1) are not both zero polynomials modulo p, by the joint coefficient-content normalization. If

\[
 Z^3\overline U(Z,1)^2=\overline V(Z,1)^2
 \quad\hbox{in }\mathbf F_p[Z],
\]

then either both sides are zero, forcing both polynomials zero, or their Z-valuations have opposite parity: 3+2 ord_Z(U) is odd and 2 ord_Z(V) is even. Both possibilities are impossible. This argument also works at p=2.

Consequently f(Z)=F(Z,1) is not identically zero modulo ANY prime. It is a primitive nonzero integer polynomial of degree at most D. In particular no coefficient-height or exceptional-prime loss is needed to apply Lemma 4. Notice that primitivity here follows from this specific square/cube shape; it is not a general assertion about arbitrary resultants.

Because y is a unit modulo N, N|F(x,y) implies f(x/y)=0 modulo N. Also any common kernel gcd g=gcd(x,y) is a unit modulo N. Homogeneity gives

\[
 F(x,y)=g^D F(x/g,y/g),
\]

so the reduced kernel ratio satisfies the same root condition.

## 3. Rational reconstruction and the coefficient-uniform count

**Lemma 5 (elementary packing).** In one unit residue class a/b=r modulo Q, the number of reduced positive fractions a/b with a,b<=Z is O(1+Z^2/Q).

Indeed distinct fractions have nonzero cross-difference divisible by Q. Their vector angles in (0,pi/2) differ by at least Q/(2Z^2), since

\[
 \sin|\theta-\theta'|
  =\frac{|ab'-a'b|}{\|(a,b)\|\,\|(a',b')\|}
  \ge\frac{Q}{2Z^2}.
\]

If this lower bound exceeds one, two fractions cannot occur. Otherwise spacing in an interval of length pi/2 gives the displayed count. The unit-denominator condition justifies the congruence; primitivity and positivity identify a unique numerator/denominator pair for each rational number.

Apply (7) to f=F(Z,1). In one of its at most D^(omega(N)) coarse classes, write x=ga,y=gb with (a,b)=1. At fixed g<=2T the primitive height is at most 2T/g. Lemma 5 bounds the number of ratios by O(1+T^2/(g^2 Q_D)). Summing over g gives

\[
 O\left(T+\frac{T^2}{Q_D}\sum_{g>=1}g^{-2}\right)
  =O\left(T+\frac{T^2}{Q_D}\right).
\]

The bounds in Theorem 1 follow, including the coefficient-uniform exponent 9/7 when d=2 and N is comparable to T^5. This is a count of kernel pairs, not a count of completions for each fixed pair. All original witness filters only decrease it.

## 4. Fixed-family contact multiplicity and the sharper clusters

At any NONZERO finite geometric root of f(Z)=F(Z,1), its characteristic-zero multiplicity is at most

\[
 M=2d+1.
 \tag{8}
\]

Here common factors of U,V are allowed for the upper bound. To prove it, first remove their polynomial gcd C over Qbar, of degree c. The remaining polynomials A,B have degrees at most d-c. At a nonzero root of Z^3 A^2-B^2, A and B are both nonzero. Choose the formal branch gamma(Z)=Z^(3/2) at that root. If h=B/A-gamma vanishes to order m, the nonzero polynomial

\[
 2Z(AB'-BA')-3AB=A^2(2Z h'-3h)
 \tag{9}
\]

vanishes to order m-1 and has degree at most 2(d-c). It is not identically zero: otherwise (B/A)^2/Z^3 would be constant, contradicting odd Z-valuation for a rational square. Thus m<=2(d-c)+1. Restoring C^2 adds at most 2c, proving (8). If one of U,V is identically zero, every nonzero finite root instead has multiplicity at most 2d, so (8) still holds.

For the fixed integer polynomial f, remove any power of Z and factor over Z into distinct primitive irreducible factors:

\[
 f(Z)=\pm Z^h\prod_i f_i(Z)^{m_i},\qquad m_i<=M.
\]

There is no nonunit content, by Section 2. Put g(Z)=product_i f_i(Z). It is squarefree. Choose integer polynomials A_0,B_0 and a fixed nonzero integer Delta with

\[
 A_0g+B_0g'=\Delta.
 \tag{10}
\]

They exist by the Euclidean algorithm over Q followed by clearing denominators. If g is constant, unit roots are absent, and the assertion is immediate.

Fix p^e exactly dividing N, set j=ceil(e/M), and s=v_p(Delta). A unit solution of f(z)=0 modulo p^e satisfies g(z)=0 modulo p^j, because
v_p(f(z))=sum m_i v_p(f_i(z))<=M v_p(g(z)).

If j<=2s, use the single class modulo 1. If j>2s, equation (10) shows v_p(g'(z))<=s. Group these solutions first by their classes modulo p^(s+1). For two solutions z,z' in one such class, put a=v_p(z-z')>=s+1 and s'=v_p(g'(z))<=s. Taylor expansion with integer coefficients shows

\[
 v_p(g(z')-g(z))=a+s',
\]

because all terms of order at least two have valuation at least 2a>a+s'. Both values are divisible by p^j, so a>=j-s. Thus this group is contained in ONE class modulo p^(j-s).

For good primes p not dividing Delta, the number of initial classes is at most deg(g): these are roots modulo p of the nonzero squarefree reduction g. For the finitely many bad primes, the crude bound p^(s+1) is a fixed constant depending only on f. Hence CRT covers the unit roots by O_f((deg g)^(omega(N))) classes modulo a number Q satisfying

\[
 Q\ge \frac{N^{1/M}}{|\Delta|^2}.
 \tag{11}
\]

Indeed j>2s gives local exponent j-s>=e/M-2s, and j<=2s permits exponent zero with the same inequality. This proves the claimed bound for Q including all bad prime powers; no good-prime condition on N is suppressed.

Now use Lemma 5 and the same g=gcd(x,y) summation as in Section 3. Theorem 2 follows. The fixed constant Delta, and therefore the constant in (11), depends on the family. When coefficients vary with N, this improvement cannot replace the coefficient-uniform bound (7).

## 5. Actual squarefree witnesses and prime moduli

For an actual distinct squarefree kernel pair, its reduced kernel coordinates a=x/g,b=y/g are coprime and squarefree. If F(a,b)=0 while U(a,b),V(a,b) are nonzero, then (a/b)^3 is a rational square. Prime-valuation parity forces a/b to be a rational square, and coprime squarefreeness forces a=b=1. This would mean x=y, excluded. Points where U=V=0 are not witness fractions and are also excluded. Therefore actual distinct squarefree witnesses always have

\[
 F(a,b)\ne0.
 \tag{12}
\]

If the family is fixed, its homogeneous coefficient bound gives

\[
 |F(a,b)|\le C_{U,V}\max(a,b)^D.
\]

Since N divides this nonzero integer, every such primitive ratio has height at least c_(U,V) N^(1/D). Consequently the number of common-kernel multipliers g for that ratio is at most C_(U,V) T N^(-1/D). Applying the fixed-family cluster proof but summing only over these g also gives the optional all-N refinement

\[
 \ll_{U,V,\varepsilon}N^\varepsilon
  \left(\frac{T}{N^{1/(2d+3)}}
             +\frac{T^2}{N^{1/(2d+1)}}\right)
 \tag{13}
\]

for ACTUAL distinct squarefree witnesses. This refinement is not asserted for invalid common-zero points or for coefficient-varying families.

For prime N=p, primitivity of f shows that there are at most D roots modulo p, regardless of repeated roots or coefficient sizes. Reduced kernel fractions of height at most 2T inject modulo p when p>8T^2, by their cross-difference. Thus there are at most D primitive kernel ratios. With gcd(x,y)<=G each has at most G admissible multipliers, proving the coefficient-uniform O_d(G) assertion. For fixed coefficients, use the sharper multiplier bound above to obtain (6). This proves Theorem 3, including the T^(2/7) quadratic-family specialization.

**Corollary 6 (separated coarse modulus, including cube-free moduli).** Let D=2d+3 and use Q_D(N) from (7). Whenever

\[
 Q_D(N)>8T^2,
\]

every degree-d family, with arbitrary coefficient dependence, has at most D^(omega(N)) primitive kernel ratios. Its gcd(x,y)<=G portion therefore has O_(d,epsilon)(N^epsilon G) pairs. For a FIXED family its ACTUAL distinct squarefree pairs satisfy

\[
 \ll_{U,V,\varepsilon}N^\varepsilon T N^{-1/(2d+3)}.
\]

Indeed Q_D divides N, so the original unit denominators remain units modulo Q_D. The cross-difference of two distinct reduced fractions of height at most 2T has absolute value less than 8T^2; it cannot be a nonzero multiple of Q_D. Hence at most one primitive ratio occurs in each of the D^(omega(N)) coarse root classes. The fixed-family height lower bound just proved bounds its common multipliers by C_(U,V)T N^(-1/(2d+3)).

Since Q_D(N)>=rad(N), the condition rad(N)>8T^2 is sufficient. In particular EVERY cube-free N comparable to T^5 satisfies it for sufficiently large T: rad(N)>=sqrt(N) is then much larger than T^2. Thus every fixed quadratic family contributes O_(U,V,epsilon)(T^(2/7+epsilon)) throughout this cube-free subclass, not just at primes. The elementary prime route in Theorem 3 remains valid independently of the general CRT root-cluster argument.

On this separated-coarse-modulus subclass, an independently proved ACTUAL cover of the gcd<=T^rho portion by J bounded-degree families would give O_(d,epsilon)(J N^epsilon T^rho) pairs. Together with the elementary O(T^(2-rho)) large-gcd tail, J<=T^(2-rho-eta) would give a saving min(rho,eta) for that subclass. No such cover is asserted, and this subclass conclusion is not a bound uniform over all N.

## 6. The actual missing cover, and its quantitative cost

Suppose one proves that every pair in a specified portion of the full K_N(T) has an ACTUAL reduced witness u/v=U_j(x,y)/V_j(x,y) for at least one of J families of degree at most d. Families may depend on N,T, but each must be normalized as above; actual pointwise nonzero, positivity, height and unit conditions must hold. Theorem 1 then gives

\[
 K_{\rm covered}\ll_{d,\varepsilon}J N^\varepsilon
         \left(T+\frac{T^2}{N^{1/(2d+3)}}\right).
 \tag{14}
\]

At d=2,N comparable to T^5, an independently established cover with J<=T^(5/7-eta) would imply a positive power saving T^(2-eta+epsilon) for that portion. This is a conditional implication, not an assertion that such a cover exists. If only one fixed finite collection independent of N,T is used, Theorem 2 gives O(T^(1+epsilon)) instead.

In the prime case, on the small-gcd portion gcd(x,y)<=T^rho, Theorem 3 gives K_covered=O_d(J T^rho). An independently established cover with J<=T^(2-rho-eta), together with the elementary O(T^(2-rho)) large-gcd tail, would imply a full prime-modulus saving with exponent min(rho,eta). No composite-modulus conclusion follows merely by proving the prime case.

A list formed by assigning one arbitrary constant witness ratio to each unknown pair is not a proved economical cover; its size is already the unknown frequency. Likewise an implicit determinant-method invocation must furnish the actual rational-function representation and the number of families, rather than only a curve containing reductions modulo N. No such full-N cover has been proved in this note. The small-gcd, long-height target and the global exponent 3/10 remain open.


## 7. Actual auxiliary curves and the true-cusp component

Let N be comparable to T^5 and let x,y be distinct squarefree kernels in [T,2T), with all x,y,u,v units modulo N, gcd(u,v)=1 and u,v<=4T^2, satisfying x^3u^2=y^3v^2 modulo N. Put

\[
 s=x/y,\qquad t=v/u,\qquad F(s,t)=t^2-s^3.
 \tag{15}
\]

These rational coordinates have unit denominators modulo N. They satisfy F(s,t)=0 modulo N, but F(s,t) is NONZERO over Q: otherwise (x/y)^3 would be a rational square, forcing x=y by squarefreeness.

Suppose an actual NONZERO integer auxiliary polynomial G(s,t), of bidegrees at most (a,b), satisfies G(s,t)=0 over Q at this actual witness point. A congruence-only relation is not sufficient for the factor-removal step below.

**Theorem 4 (auxiliary curves and exact effective modulus).** Under the preceding actual-witness hypotheses, remove every exact F=t^2-s^3 factor from G. Its nonzero remainder is A(s)+tB(s). If c is the joint coefficient content of A,B, then R=A^2-s^3B^2 has content EXACTLY c^2 and degree at most D=2a+3b. The actual point forces the primitive norm R/c^2 to vanish modulo N_G=N/gcd(N,c^2). For D>=1 a nonconstant primitive norm gives O_(D,epsilon)(N^epsilon[T+T^2/N_G^(1/D)]) covered kernel pairs. A constant primitive norm has no points when N_G>1 and gives only the trivial bound when N_G=1. The theorem supplies the conditional cover consequences in Section 9; it does not assert the existence of a cover.

Divide G by as many exact factors F=t^2-s^3 over Q[s,t] as possible. Monicity in t ensures these divisions are in Z[s,t]. At the actual rational point F is nonzero, so the remaining polynomial G_0 still vanishes there. Its bidegrees do not exceed those of G. Multiplicativity of each separate variable degree over the other polynomial coefficient domain proves this assertion. The remaining G_0 is not divisible by F. Scaling G_0 to its primitive integer representative does not change its rational zero set.

Divide G_0 by F as a polynomial in t. Its unique remainder is

\[
 A(s)+tB(s),\qquad A,B\in\mathbf Z[s].
 \tag{16}
\]

They are not both zero, because all exact F factors were removed. Reduction of each monomial s^i t^j gives

\[
 \deg A\le a+3\lfloor b/2\rfloor,\qquad
 \deg B\le a+3\lfloor(b-1)/2\rfloor,
 \tag{17}
\]

with the second bound omitted if b=0 and B=0.

Define the norm polynomial

\[
 R(s)=A(s)^2-s^3B(s)^2.
 \tag{18}
\]

It is nonzero, by the odd s-valuation obstruction to a rational square if A,B are both nonzero; the one-zero cases are immediate. Also

\[
 \deg R\le2a+3b=:D.
 \tag{19}
\]

The same expression is the resultant Res_t(t^2-s^3,G_0), since multiplication by A+tB in the free rank-two quotient algebra has matrix [[A,s^3B],[B,A]]. No general elimination theorem is needed for the counting argument.

At the actual point, F and G_0 vanish modulo N. Thus A(s)+tB(s)=0 modulo N and t^2=s^3 modulo N, so R(s)=0 modulo N. This uses only unit denominators; no division by an auxiliary leading coefficient is made.

## 8. Exact remainder content and modulus loss

Let c be the positive gcd of ALL coefficients of A and B together. Write A=cA_1,B=cB_1. Then

\[
 R=c^2 R_1,\qquad R_1=A_1^2-s^3B_1^2,
 \tag{20}
\]

and R_1 is primitive at EVERY prime. Indeed if R_1 were the zero polynomial modulo p, either A_1,B_1 would both vanish modulo p or the s-valuations of their two square terms would have opposite parity. Both are impossible. Therefore

\[
 \operatorname{content}(R)=c^2
 \tag{21}
\]

EXACTLY. In particular primitivity of the original G does not imply primitivity of its norm.

The exact effective modulus forced by this argument is

\[
 N_G=\frac{N}{\gcd(N,c^2)}.
 \tag{22}
\]

For every prime exponent, N|c^2 R_1(s) implies N_G|R_1(s), after multiplying the unit denominator by a sufficiently high power. This is a consequence of valuations, not an invalid inversion of c modulo N. No claim that N itself divides R_1 is made when c shares primes with N.

For D>=1 and nonconstant R_1, Lemmas 4--5 above applied to R_1 give the coefficient-uniform bound

\[
 \#\{\hbox{unit kernel pairs covered by this actual curve}\}
 \ll_{D,\varepsilon}N^\varepsilon
       \left(T+\frac{T^2}{N_G^{1/D}}\right).
 \tag{23}
\]

Here D can be the upper bound 2a+3b, and no coefficient-height bound is needed. If R_1 is the primitive constant +1 or -1, there are no points when N_G>1; when N_G=1 the trivial O(T^2) bound applies. The same trivial bound also covers D=0. Thus a curve with complete modulus loss gives no useful savings merely by being primitive.

A simple algebraic example is G=F+N: its remainder is A=N,B=0, so c=N, R=N^2 and N_G=1. There is a stronger example on the actual witness locus. Write its nonzero F(s,t)=K/L in lowest terms. The denominator L is coprime to N, and N divides K. The primitive integer polynomial

\[
 G=L(t^2-s^3)-K
\]

vanishes at this actual point; gcd(K,L)=1 makes it primitive. Its remainder is -K, so again N_G=1. Hence even an actual bounded-bidegree curve with coefficients of polynomial size can lose the ENTIRE target modulus. An auxiliary-cover theorem must control this remainder content explicitly.

## 9. A precise sufficient cover target for the long complement

Fix the exact long complement from the short-witness theorem: pairs having some witness of height at most 4T^2 but NO witness of height at most T^(5/4-delta), for one fixed delta in (0,1/4]. The short part is already O(T^(2-2delta+epsilon)).

Suppose, as an additional hypothesis, that each pair in this long complement has a chosen actual witness point on at least one of J integer auxiliary curves G_j=0 over Q. Remove exact F factors as in Section 7. Assume all reduced curves have bidegrees bounded by fixed (a,b), and after the remainder calculation each has

\[
 N_{G_j}\ge T^\nu
\]

for a fixed nu>0. For nonconstant norm polynomials, summing (23) gives

\[
 K_{\rm long}\ll_{D,\varepsilon}
       J T^\varepsilon\left(T+T^{2-\nu/D}\right),
 \qquad D=2a+3b.
 \tag{24}
\]

Constant primitive norm polynomials with effective modulus greater than one have no points and may be discarded. If an independently proved cover satisfies J<=T^kappa with

\[
 \kappa<\min(1,\nu/D),
\]

then the long complement has a power saving with exponent

\[
 \eta=\min(1-\kappa,\nu/D-\kappa)>0.
 \tag{25}
\]

The exact disjoint short/long partition then gives the full kernel-pair bound with saving min(eta,2delta). This is an implication under actual-cover AND effective-modulus hypotheses, not an assertion of either. It retains the original uniform composite-N target.

For quadratic rational witness families, G=U(s)t-V(s), so (a,b)=(2,1) and D=7. The joint coefficient normalization is exactly c=1, reproducing the coefficient-uniform T^(9/7+epsilon) family bound. General auxiliary curves need not have c=1, as Section 8 demonstrates.

For a prime N=p, a curve with p not dividing c has a nonzero primitive norm modulo p and at most D primitive kernel ratios when p>8T^2. Thus its small-gcd contribution with gcd(x,y)<=T^rho is O_D(T^rho). A cover of that prime-case portion by J such curves would give O_D(J T^rho), together with the elementary large-gcd tail. This is a prime-modulus consequence only; it is not silently promoted to all N.

The unresolved mathematical work is now explicit: construct an economical collection of ACTUAL Q-auxiliary curves for the long pairs and prove that their remainders retain a sufficiently large effective modulus. Counting completions on one fixed cubic slice, asserting that an auxiliary polynomial exists, or checking its primitivity alone does not establish those two requirements. No new full-K estimate follows from the present identities without them.

## 10. Classical antecedents and exact source scope

Lemma 4 is classical. C. L. Stewart, *On the number of solutions of polynomial congruences and Thue equations*, Journal of the American Mathematical Society 4 (1991), 793--835, Theorem 2 on printed page 798, gives the stronger prime-power root-class structure for a p-primitive polynomial of nonzero discriminant. Its stated exponent bound implies classes modulo p^(ceil(e/D)). For repeated-root polynomials, adding c p^e to the constant coefficient for an integer c outside finitely many critical values preserves the congruence, degree and p-primitivity and makes the discriminant nonzero. The local interpolation proof above proves exactly the needed statement directly. No root-clustering novelty is claimed.

Actual author primary:
https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf .
The angular packing step is the existing elementary rational-reconstruction argument in [the dual kernel-graph proof](ERDOS_943_DUAL_KERNEL_GRAPHS.md). The polynomial multiplicity, primitive cusp shape and family counts are written above in full. No Blomer--Schobel fixed-difference estimate, fixed cubic-slice completion theorem, character-sum estimate, finite experiment or unproved auxiliary cover is a premise.

The short/long partition used in Section 9 is [Short witnesses for compatible square-full kernels](ERDOS_943_SHORT_WITNESS_KERNEL_PAIRS.md), Corollary 5 and Theorem 7. It counts the exact minimum-height complement, not all long witness tuples. Its proved short saving is the only input needed to transfer the conditional long-cover saving to the full pair count.

## 11. A universal cap for ACTUAL witness multipliers

The bounds for the larger raw congruence sets in Theorems 1--2 remain unchanged. The following refinements apply only to pairs having an ACTUAL witness with all original distinct-squarefree, positive, reduced, unit and height conditions. Keep the fixed balanced constants c,C from the introduction.

**Proposition 7 (actual-witness multiplier cap).** If an actual witness has height at most H<=4T^2, then

\[
 g=\gcd(x,y)\le G_H:=\left(\frac{16T^3H^2}{N}\right)^{1/3}
                \ll_c T^{2/3}.
 \tag{26}
\]

Indeed, write x=ga,y=gb. Since x,y are units modulo N, so is g, and cancellation gives N|a^3u^2-b^3v^2. This difference is nonzero by the squarefree parity argument in Section 5: equality would force the coprime squarefree a,b both to equal one, hence x=y. Its absolute value is at most16(T/g)^3H^2. A nonzero multiple of N has magnitude at least N, proving (26). If G_H<1 no such pair exists. This is the universal height/gcd relation already stated in [the kernel-pair discussion](POWERFUL_SUM_GRAPHS.md); no diagonal or other special family is a hypothesis.

The actual hypotheses matter. With U=V=1, the raw congruence F=x^3-y^3=0 includes order T diagonal unit pairs at a large prime N, so its T term cannot be replaced by T^(2/3). If squarefreeness is removed, x=9h,y=16h for h in [T/9,T/8) and the fixed reduced witness u/v=64/27 satisfy exact integer equality x^3u^2=y^3v^2, with g=h of order T. At primes comparable to T^5 they satisfy the unit and full-height conditions. Thus even distinct positive reduced witnesses alone do not imply (26).

**Corollary 8 (actual-only truncated counts).** In the generic multiplier summation of Section 3, restrict g to g<=G_H. One coarse class now contributes

\[
 \sum_{g\le G_H}O\left(1+\frac{T^2}{g^2Q}\right)
                       \ll G_H+\frac{T^2}{Q}.
\]

Consequently an actual degree-d rational witness family, with arbitrary coefficient dependence, contributes

\[
 \ll_{d,\epsilon,c}N^\epsilon
       \left(T^{2/3}+\frac{T^2}{N^{1/(2d+3)}}\right).
 \tag{27}
\]

Here Section 2 supplies the primitive cusp polynomial and Section 1 supplies Q>=N^(1/(2d+3)); the original witness filters are retained. For a fixed family, intersect (26) with the separate actual height bound in Section 5 and use its fixed-family root cover from Section 4. This gives

\[
 \ll_{U,V,\epsilon,c}N^\epsilon
 \left(\min\{T^{2/3},T/N^{1/(2d+3)}\}
                         +\frac{T^2}{N^{1/(2d+1)}}\right),
 \tag{28}
\]

with the same coefficient-dependent constants as before.

For one actual auxiliary curve, use its post-cusp remainder norm and EXACT effective modulus M=N_G from Sections 7--8. For a nonconstant primitive norm with any legitimate positive degree bound D, the same summation gives

\[
 \#\text{actual pairs on the curve}\ll_{D,\epsilon,c}
 N^\epsilon\left(T^{2/3}+\frac{T^2}{M^{1/D}}\right).
 \tag{29}
\]

This counts distinct pairs, not witness points. The original actual witness supplies (26), even if the auxiliary norm happens to vanish over Q at that kernel ratio. If the primitive norm is constant, it is +/-1: M>1 gives no actual points, whereas M=1 has only the trivial O(T^2) bound. A zero norm is excluded by the valid factor-removal and nonzero norm arguments in Section 7. D=0 is never inserted into1/D.

If Q_D(M)>8T^2, Corollary 6 instead gives at most D^omega(M) primitive ratios, each with O_c(T^(2/3)) common multipliers. Thus the whole actual collection is O_(D,epsilon,c)(N^epsilon*T^(2/3)). At a retained prime M=N=p>8T^2 the direct degree-root argument gives O_(D,c)(T^(2/3)), with no large-gcd tail. For rational families this also applies at every sufficiently large cube-free balanced N. For auxiliary curves the condition is on N_G after content loss, never merely on the original N.

For J actual curves with fixed norm-degree bound D and every N_G>=T^nu, (29) gives J*T^epsilon*(T^(2/3)+T^(2-nu/D)). Hence J<=T^kappa and kappa<min(4/3,nu/D) suffice for a saving min(4/3-kappa,nu/D-kappa) on the covered pairs. On the exact long complement, the full saving is additionally limited by the proved short exponent2delta. This is a sufficient unproved-cover route. In particular the general quadratic family exponent9/7 is unchanged, since its second term still dominates.

## 12. Primitive ratios, exact witness classes and the pair count

For any subset S of the original full-height kernel pairs, write K(S)=|S| and let R(S) be its number of ordered primitive ratios a/b, represented by unique positive coprime pairs. Reciprocal ratios are distinct. The empty set has K=R=0.

**Proposition 9 (pair--ratio and witness-class identities).** For every such subset, and in fact for any subset of ordered positive integer pairs in [T,2T),

\[
 R(S)\le K(S)\le G R(S)+O(T^2/G)\quad(G\ge1),
 \qquad K(S)\ll T\sqrt{R(S)}.
 \tag{30}
\]

There are at most floor(G) common multipliers at a fixed primitive ratio with gcd at most G. The remaining pairs are bounded by sum_(g>G)(2T/g)^2=O(T^2/G). For1<=R<=T^2 take G=T/sqrt(R); for R>T^2 use K=O(T^2). This proves (30) including its endpoint cases. For actual full-height pairs, (26) separately gives K(S)<<_c T^(2/3)R(S).

Witness membership depends only on the ratio: cancellation of the unit g^3 shows that the entire set of original positive reduced unit fractions satisfying the congruence is the same for all admissible multipliers. Define h_N(a,b) to be the least integer witness height, or infinity if there is no witness. For a primitive squarefree ratio, admissible multipliers are exactly the squarefree integers g with

\[
 \frac{T}{\min(a,b)}\le g<\frac{2T}{\max(a,b)},
              \qquad(g,ab)=1,\qquad(gab,N)=1.
 \tag{31}
\]

Let M_(T,N)(a,b) count them. On every complete class defined by a condition on h_N(a,b), the exact pair count is sum_(a/b) M_(T,N)(a,b). This equality is not asserted for an arbitrary incomplete subset S.

## 13. Content-aware capacity and its forward bound

An actual cover of S chooses, for every pair in S, at least one ORIGINAL full-height witness at a rational zero of a NONZERO integer polynomial G_j(s,t). It need not cover every witness. Each curve undergoes exactly the actual-Q-zero factor removal and primitive normalization in Sections 7--8, giving a primitive nonzero norm f_j, its content parameter c_j, and

\[
 M_j=N/\gcd(N,c_j^2).
\]

For each nonconstant norm choose its OWN positive integer degree bound D_j>=deg(f_j); it may be the actual norm degree. Do not force an unrelated larger global degree onto a curve whose true degree is smaller. Put Q_D(1)=1 and use Q_D(M) from (7) for M>1.

Define its capacity by

\[
 {\cal A}_T(G_j)=
 \begin{cases}
  1+T^2,& M_j=1,\\
  0,& M_j>1\text{ and }f_j\text{ is constant},\\
  D_j^{\omega(M_j)}(1+T^2/Q_{D_j}(M_j)),
     & M_j>1\text{ and }f_j\text{ is nonconstant}.
 \end{cases}
 \tag{32}
\]

A constant primitive norm at M_j>1 has no actual points and may be discarded. Complete content loss M_j=1 MUST be charged order T^2, including constant norms; normalization of G alone cannot turn it into a free cover. The zero-norm alternative is not legitimate. Empty sets may use the empty cover with capacity zero.

**Theorem 10 (capacity and the all-modulus vertical reverse).** For every finite actual cover with total capacity A=sum_j A_T(G_j),

\[
 R(S)\ll A,\qquad
 K(S)\ll T\sqrt A,\qquad
 K(S)\ll_c T^{2/3}A.
 \tag{33}
\]

Conversely, in any permitted class containing the verticals with their own norm degree2, there is an actual cover with A<<_(epsilon,c,C)T^epsilon R(S), as proved in Section 14.

Indeed Lemma 4 and CRT cover the norm roots by at most D_j^omega(M_j) classes modulo Q_(D_j)(M_j). Lemma 5 at primitive height2T bounds each class by O(1+T^2/Q). All denominator units persist because M_j divides N. The constant and M_j=1 branches follow directly as in Section 8. Summing primitive-ratio bounds gives R(S)<<A, and Proposition 9 gives (33). The packing constant is absolute; all degree dependence in (32) is explicit.

There is also the direct bound obtained by summing (29) when degrees are uniformly bounded. It may be sharper than (33). These are complementary upper bounds, not an inference that the actual cover exists or has small capacity.

## 14. Reverse vertical covers, individual degrees and class guards

For the reverse direction the permitted class must contain ALL vertical lines below with their legitimate norm-degree assignment2. Fixed bidegree caps alpha>=1,beta>=0 suffice if no additional restriction excludes those lines. A preassigned smaller class forbidding s-dependence has only the forward implication.

For every primitive ratio a/b in R(S), take

\[
 G_{a,b}(s,t)=bs-a.
\]

This integer curve contains EVERY actual witness of its ratio. Its coefficients are coprime; no true-cusp factor occurs; its remainder has A(s)=bs-a,B(s)=0 and c=1. Therefore

\[
 f_{a,b}=(bs-a)^2,\qquad D_{a,b}=2,\qquad M_{a,b}=N,
\]

and at EVERY balanced integer modulus

\[
 Q_2(N)\ge\sqrt N\ge\sqrt c\,T^{5/2},\qquad
 {\cal A}_T(G_{a,b})
 \le2^{\omega(N)}(1+c^{-1/2}T^{-1/2})
 \ll_{\epsilon,c,C}T^\epsilon.
 \tag{34}
\]

The last step is the elementary subpower bound from Section 1 with epsilon reallocation using N comparable to T^5. The degree TWO is essential to this capacity calculation; a loose unrelated common upper degree can destroy it. No prime-modulus assumption or new coefficient-field point is used.

Thus the R(S) verticals give an ACTUAL cover with

\[
 A\ll_{\epsilon,c,C}T^\epsilon R(S).
 \tag{35}
\]

If B(S) is the infimum of total capacities over any permitted class including these verticals, (33)--(35) give R(S)<<B(S)<<_epsilon T^epsilon R(S). No assertion that the infimum is attained is needed: (35) gives an explicit finite cover for all quantitative reverse implications.

For comparison, at N=p>8T^2 consider a fixed permitted class in which EVERY assigned nonconstant norm bound satisfies D_j<=D, with D>=2 and verticals admitted with D_j=2. Only curves retaining p are counted. Such a curve covers at most D primitive ratios by the direct prime argument in Section 5; the vertical cover gives the reverse inequality. Its minimum retained curve NUMBER C(S) therefore satisfies

\[
                         R(S)/D\le C(S)\le R(S).
 \tag{36}
\]

This comparison concerns ratios and curve number. Equation (30), rather than a same-exponent identity, relates them to the kernel-pair count. Total capacity is the formulation valid for ALL moduli and all content losses.

## 15. Uniform power-saving directions and the exact long complement

A power-saving estimate with exponent eta means C_epsilon*T^(2-eta+epsilon) for EVERY epsilon>0, uniformly for the balanced N,T with fixed parameters in constants. From (30), a K bound with exponent eta gives the SAME exponent for R; an R bound gives K exponent eta/2. Equations (33)--(35) similarly show:

* a K bound gives existence of an actual vertical cover with capacity exponent eta;
* a capacity bound with exponent eta gives K exponent eta/2;
* an R bound gives a vertical capacity bound with the same exponent.

All epsilon factors can be reallocated, since each assertion holds for every positive epsilon. In particular existence of SOME positive saving is equivalent for the pair count, primitive-ratio count and actual-cover total capacity in every permitted class containing the verticals. There is no same-exponent pair/capacity equivalence.

Fix delta in(0,1/4], H_0=floor(T^(5/4-delta)) and H_full=floor(4T^2). The full, short and exact long classes are defined by h_N<=H_full, h_N<=H_0 and H_0<h_N<=H_full respectively. A pair with both short and long witness tuples belongs only to the short class. Proposition 9 gives BOTH exact disjoint set partitions

\[
 P_{\rm full}=P_{\rm short}\sqcup P_{\rm long},\qquad
 R_{\rm full}=R_{\rm short}\sqcup R_{\rm long}.
 \tag{37}
\]

The established short theorem in Section 10 gives K_short<<T^(2-2delta+epsilon). Thus an actual cover of the exact long set with capacity exponent eta implies

\[
 K_{\rm long}\ll T^{2-\eta/2+\epsilon},\qquad
 K_{\rm full}\ll
               T^{2-\min(\eta/2,2\delta)+\epsilon}.
 \tag{38}
\]

Conversely a full-pair saving eta gives a long vertical cover with the SAME capacity exponent, because R_long<=K_full. This proves both directions of existence of a positive saving for full K and for the total capacity of a cover of the exact long complement. It is not a theorem about arbitrary long witness tuples.

A useful sufficient numerical gate follows when J<=T^kappa curves have uniformly bounded assigned norms D_j<=D and every M_j>=T^nu. Constant norms have no assigned points and may be discarded. Root capacities and Q_Dj(M_j)>=M_j^(1/D_j)>=T^(nu/D) give

\[
 A\ll_\epsilon J T^\epsilon(1+T^{2-\nu/D}),\qquad
 \eta_K=\min\{1-\kappa/2,(\nu/D-\kappa)/2\}>0
 \quad\text{if }\kappa<\min\{2,\nu/D\}.
 \tag{39}
\]

Here eta_K is the saving on the COVERED pair set; for an exact-long cover its full consequence is further limited by2delta. This broadens the numerical J gate; the stronger direct bound (29) is still available when useful. At separated coarse moduli with uniformly bounded D_j all curve capacities are subpower, so any actual J<=T^(2-eta) cover suffices for covered K exponent eta/2.

For completeness, if only the g<=G part of the exact long class is covered with capacity A_G, its pair contribution is O(G A_G). Adding the unconditional large-gcd tail gives K_long<<G A_G+T^2/G. If G=T^rho, 0<rho<=1, and A_G<<T^(kappa+epsilon), then the full saving is

\[
            \min\{2-\kappa-\rho,\rho,2\delta\}>0
                \quad\text{when }\kappa+\rho<2.
 \tag{40}
\]

Small-gcd restrictions retain the exact witness-class partitions. No bound on all long witness tuples is substituted.

## 16. A new OPEN capacity predicate equivalent to the original target

**Corollary 11 (positive-saving equivalences).** Fix one delta in(0,1/4], for example delta=1/8. The following existential uniform statements are equivalent:

1. The original full distinct-kernel-pair count K_N(T) has SOME positive power saving below T^2.
2. Its primitive ordered ratio count has SOME positive power saving below T^2.
3. Full kernel pairs admit actual covers with SOME positive-saving total capacity, in a fixed bounded-degree class containing all verticals with their own norm degree2.
4. The exact NO-short/SOME-full complement admits such actual covers with SOME positive-saving total capacity.

Statement1 implies2 with the same exponent, and the vertical construction gives3 or4 with the same exponent. The reverse from3 loses at most a factor2 in the saving exponent; the reverse from4 also takes the minimum with2delta. Thus these are genuine equivalences of the EXISTENCE of a positive saving, with the quantitative directions specified above. No estimate in this list has been proved merely by proving the equivalence.

**Open problem 12 (positive-saving exact-long capacity).** Prove that there exist eta>0 and one fixed bounded-bidegree permitted class containing the verticals such that, for every epsilon>0 and all balanced N,T, the exact long complement at delta=1/8 has an actual finite cover with

\[
                  \sum_j{\cal A}_T(G_j)
                         \le C_\epsilon T^{2-\eta+\epsilon}.
 \tag{41}
\]

Every curve must be NONZERO and vanish OVER Q at a chosen original witness point; its own norm degree and exact post-content modulus enter its capacity. The curve coefficients and the chosen finite cover may depend on N,T (and, in this existence statement, epsilon). Empty sets may use the empty cover. There is no new coefficient-field point, algorithm or auxiliary-cover construction assumed.

This OPEN predicate is equivalent to the original positive-saving target by Corollary 11. The older strict J/retained-modulus obligation in Section 9 remains a merely SUFFICIENT unproved route; it is neither asserted necessary nor silently replaced by a proved cover. Vertical covers establish the reverse implication only after a pair or ratio saving is GIVEN. Choosing curves from the unknown ratio set is not an independent proof of that saving.

The new actual-count and capacity statements do not alter the raw congruence bounds in Theorems 1--2, the content-loss example in Section 8, or the original open frequency targets. No positive-saving capacity estimate, economical full cover, full kernel-pair saving, full six-cycle target, or improvement below the proved prescribed-sum exponent3/10 is established here. All root, packing and norm facts reused in these sections were proved earlier in this same source; the exact short-cutoff input remains the one identified in Section 10.

## 17. A coefficient-uniform structured unit-root refinement

The square/cube norm admits a stronger unit-root bound than the general degree bound. The ordinary polynomial degree and the unit-root contact order play different roles.

**Theorem 13 (uniform structured unit-root cover).** Fix d>=0, set M=2d+1 and D=2d+3, and let U,V in Z[s] have degrees at most d and joint coefficient content1. Either polynomial may be zero, but not both. Put f(s)=s^3U(s)^2-V(s)^2. There are constants A_d,B_d>0, depending only on d, such that for EVERY N>=1 the unit roots of f modulo N occupy at most B_d D^omega(N) classes modulo a divisor Q_d(N) of N with

\[
                         Q_d(N)\ge A_d^{-1}N^{1/M}.
\tag{42}
\]

The divisor can be chosen independently of U,V. At every prime p>D, at most D classes modulo p^ceil(e/M) suffice for unit roots modulo p^e. For the finitely many p<=D there are constants C_(d,p)>=0, independent of the coefficients, and one may take

\[
 Q_d(N)=\prod_{p^e\parallel N}
 p^{\max\{0,\lceil e/M\rceil-\lceil C_{d,p}/M\rceil\}},
 \qquad C_{d,p}=0\quad(p>D).
\tag{43}
\]

The fixed small-prime losses are absorbed in A_d,B_d. In particular the number of coarse classes is O_(d,epsilon)(N^epsilon), uniformly in the coefficients. Class membership remains only a necessary root condition.

The unit restriction is essential: U=s^d,V=0 gives f=s^(2d+3), whose roots near0 need the larger degree exponent. No assertion concerning nonunit roots is substituted for the earlier general root lemma.

## 18. Primitivity and contact, including common roots

The polynomial f is nonzero and primitive at every prime. Over any field an identity s^3U^2=V^2 with U,V nonzero would equate an odd s-valuation with an even one. Thus f identically zero modulo p would force both U,V identically zero modulo p, contrary to joint coefficient primitivity.

We prove that every nonzero root of f over a characteristic-zero field has multiplicity at most M. The same assertion holds over a field of characteristic p>D.

If one of U,V is identically zero, the multiplicity of any nonzero root is at most2d. Otherwise fix c!=0. If the orders of U,V at c differ, f has order twice their smaller order, again at most2d. If the orders agree and equal r, write U=(s-c)^rU_1,V=(s-c)^rV_1. Both new polynomials are nonzero at c and have degree at most d'=d-r. Write m=ord_c(s^3U_1^2-V_1^2). If m=0, the order of f is2r<=2d and the proof is complete. Hence assume m>=1 and put

\[
 W=3U_1V_1+2s(U_1'V_1-U_1V_1').
\]

This polynomial has degree at most2d' and is not identically zero. Indeed, if a,b are the lowest nonzero s-exponents of U_1,V_1, the lowest coefficient in W is their coefficient product times3+2a-2b. This nonzero odd integer has absolute value at most2d'+3<=D, so it is nonzero also in characteristic p>D.

With f_1=s^3U_1^2-V_1^2, direct differentiation gives

\[
                  V_1f_1'-2V_1'f_1=s^2U_1W.
\tag{44}
\]

In characteristic0 the left side has order m-1 at c. In characteristic p>D this remains true, since m<=deg(f_1)<=2d'+3<p. Therefore m-1<=deg(W)<=2d', and the original order is2r+m<=2d+1. This proves the contact bound without assuming that U,V are coprime as polynomials or that f has nonzero discriminant.

## 19. Monic interpolation and the good-prime cover

We use an elementary remainder observation. Let g in Z_p[X], choose z_0,...,z_m in Z_p distinct modulo p^k, and suppose every g(z_i) is divisible by p^e. Division by the monic polynomial P(X)=product_i(X-z_i) gives an integral remainder R of degree at most m. Its Lagrange formula shows that every coefficient has valuation at least

\[
                          v_p(\operatorname{coeff}R)
                                \ge e-m(k-1).
\tag{45}
\]

Each interpolation denominator is a product of m differences, each with valuation at most k-1; every numerator polynomial is integral. This argument applies to the REMAINDER, even when deg(g)>m.

Let p>D and let a be a unit root of the nonzero reduction fbar, with multiplicity m_a<=M. If its lifts modulo p^e occupied more than m_a classes modulo p^ceil(e/m_a), choose representatives in m_a+1 such classes and apply (45). The remainder is zero modulo p. All representatives reduce to a, so the monic divisor reduces to (X-a)^(m_a+1). The division identity would force this power to divide fbar, a contradiction.

Thus this residue disk needs at most m_a classes modulo p^ceil(e/m_a), and hence at most m_a classes modulo the coarser p^ceil(e/M). The sum of all unit-root multiplicities in fbar is at most deg(fbar)<=D. This gives the good-prime assertion, including reductions of lower degree, zero reductions of one polynomial, and characteristic-zero repeated roots.

## 20. Uniform small-prime jets and classes

Fix one prime p<=D. The set of coefficient vectors of U,V in Z_p of degrees<=d with at least one unit coefficient, together with c in Z_p^*, is compact. For j=0,...,M let tau_j(U,V,c) be the coefficient of Z^j in f(c+Z), its jth Hasse derivative. These are integral continuous polynomial functions.

They cannot all vanish at any point of this compact set: the characteristic-zero contact bound just proved would be violated. Consequently there is an integer C=C_(d,p)>=0 for which

\[
                \min_{0\le j\le M}v_p(\tau_j(U,V,c))\le C
\tag{46}
\]

uniformly on the ENTIRE coefficient space. To see this directly, a sequence violating every proposed bound has a convergent subsequence. Its limit remains primitive with c a unit, and continuity would make all M+1 jets zero. Compactness of Z_p^n follows by selecting constant residue subsequences modulo successive powers of p. Thus no coefficient-dependent discriminant constant is hidden in C.

Set r=C+1. If e<=C, the exponent in (43) is0. Otherwise put k=ceil((e-C)/M) and b=max(0,ceil(e/M)-ceil(C/M)), so b<=k. When k<=r, at most p^r classes modulo p^k cover all units and can be coarsened to exponent b.

When k>r, partition the units into disks modulo p^r. Suppose one disk contained roots in M+1 distinct classes modulo p^k. Their monic remainder R in (45), now with m=M, has every coefficient divisible by p^(C+1), since

\[
                         e-M(k-1)>C.
\tag{47}
\]

Choose an integral unit center c of the disk. Every term of the jth Hasse derivative of P at c, for j<=M, contains M+1-j factors c-z_i divisible by p^r. These derivatives all have valuation greater than C. The integral identity f=PQ+R and the Hasse product rule would therefore make every jet in (46) divisible by p^(C+1), a contradiction.

There are fewer than p^r unit disks, so at most M p^r classes modulo p^k, and therefore modulo p^b, suffice. This is a constant depending only on d,p. The proof uses neither a preparation theorem nor a family-dependent bad-prime bound.

## 21. CRT and the classical root-congruence antecedent

Combining the local classes by CRT proves Theorem13. Indeed, (43) satisfies (42) with

\[
                      A_d=\prod_{p\le D}p^{\lceil C_{d,p}/M\rceil}.
\]

The finitely many small-prime class counts contribute B_d; all other counts contribute at most D^omega(N). The elementary subpower estimate from Section1 completes the uniform assertion. The convention for N=1 is one class modulo1. This includes d=0, where M=1.

The congruence-class mechanism is classical. Stewart's Theorem2, printed798 with proof798--803, gives p-adic root disks and the ordinary-degree coarse cover; its stated discriminant hypothesis and the comparison are discussed in Section10. Here the structured contact bound, all small-prime constants, and the coefficient-uniform cover were proved locally in Sections18--20. Stewart is credited as an antecedent, not asserted to state this specific family refinement. No historical-priority conclusion is required.

## 22. Coefficient-uniform rational-family consequences

**Corollary 14 (uniform raw and actual family bounds).** Let homogeneous integer U_h,V_h of the same degree d>=0 be jointly coefficient-primitive. The raw count is for unit pairs x,y in[T,2T) satisfying N|x^3U_h(x,y)^2-y^3V_h(x,y)^2; the actual count retains the original positive reduced unit-witness, distinct-squarefree and height conditions. Theorem13 gives

\[
 \begin{aligned}
 \#\{\hbox{raw unit pairs}\}
    &\ll_{d,\epsilon}N^\epsilon
                  \bigl(T+T^2/N^{1/(2d+1)}\bigr),\\
 \#\{\hbox{actual full-height pairs}\}
    &\ll_{d,\epsilon,c,C}N^\epsilon
                  \bigl(T^{2/3}+T^2/N^{1/(2d+1)}\bigr).
 \end{aligned}
\tag{48}
\]

Here the raw set means the larger unit-pair congruence set, while the second line retains an ORIGINAL positive reduced unit witness of height<=4T^2, squarefree distinct kernels, and the balance cT^5<=N<=CT^5. Coefficients may vary arbitrarily with N,T.

To prove this, dehomogenize at the unit denominator to obtain f(s)=s^3U_h(s,1)^2-V_h(s,1)^2. The actual-witness congruence implication was proved in Section1 by multiplying pointwise gcds, never by inverting them. Theorem13 applies to s=x/y, a unit modulo N. Angular packing and the common-gcd sum give O(G_max+T^2/Q) per coarse class. Take G_max=2T for raw unit pairs and the independent Proposition7 cap G_max=O(T^(2/3)) for actual pairs. Sum the subpower number of classes and use (42).

In particular, for N comparable to T^5,

\[
 \begin{array}{ll}
 d=2:&\#\hbox{raw or actual family pairs}
                      \ll_\epsilon T^{1+\epsilon},\\
 d=0\hbox{ or }1:&\#\hbox{actual family pairs}
                      \ll_\epsilon T^{2/3+\epsilon}.
 \end{array}
\tag{49}
\]

For d=0, the second term in (48) is T^2/N; for d=1 it is of order T^(1/3). These are coefficient-uniform bounds, not fixed-family asymptotics.

A GIVEN actual cover by J<=T^kappa quadratic families therefore contributes O_epsilon(T^(1+kappa+epsilon)) pairs, a saving for kappa<1. A GIVEN cover by degree-at-most1 actual families contributes O_epsilon(T^(2/3+kappa+epsilon)), a saving for kappa<4/3. On the exact long set, the final full-count saving is additionally limited by2delta using the short estimate and exact partition. Neither assertion constructs a cover.

The nonzero integer-height estimate in Section4 still uses the ORDINARY degree2d+3. It is not replaced by2d+1: contact controls unit congruence clusters, not the archimedean size of F(a,b). All content-loss and degree gates for arbitrary auxiliary curves remain as stated earlier. There is no assertion here that an arbitrary curve norm has the special balanced-degree family shape required by Theorem13.

## 23. An actual one-ratio sharpness construction

**Proposition 15 (sharp actual multiplier exponent).** There are infinitely many balanced integer moduli N and scales T, with a subset of at least T^(2/3-o(1)) original distinct squarefree unit kernel pairs in [T,2T), such that:
- every pair has the SAME primitive ratio;
- one positive reduced unit witness of height<=4T^2 works for the whole subset;
- for every fixed delta in(0,1/4], the entire subset is in the exact NO-short/SOME-full complement;
- one primitive retained-content vertical curve, and one jointly primitive constant rational-witness family, cover the subset.

Thus the universal actual per-ratio, per-vertical and constant-family exponent2/3 cannot be decreased by a fixed positive amount over all balanced integer moduli. No prime-modulus or full-ratio-count assertion is made.

Choose arbitrarily large B with B=1 modulo20 and both B,B+1 squarefree; their existence is proved in Section24. Put

\[
 a=B,\quad b=B+1,\quad T=B^3,\quad
 v=B^6+1,\quad u=2B^6+3=2v+1,\qquad
 N=a^3u^2-b^3v^2.
\tag{50}
\]

Then N~3B^15=3T^5. It is positive for B>=2 because u/v>2 and (b/a)^3<4. Moreover (u,v)=1, (a,v)=1 and (b,u)=(B+1,5)=1, the last using B+1=2 modulo5. Together with (a,b)=1, reducing N modulo a,b,u,v gives

\[
                 (abuv,N)=1,\qquad
                 \max(u,v)=2B^6+3\le4T^2.
\tag{51}
\]

Let I_B=[B^2,2B^3/(B+1)) and Q=abN. For every squarefree g in I_B with (g,Q)=1, set x=ga,y=gb. These are distinct squarefree kernels, their gcd is exactly g, both are in [T,2T), and both are units modulo N. The same positive reduced witness works because

\[
 |I_B|=B^2\frac{B-1}{B+1}\sim B^2,\qquad
 x/y=B/(B+1),\qquad
 x^3u^2-y^3v^2=g^3N.
\tag{52}
\]

The subset has exactly ONE primitive ratio; this is not a statement about the full kernel-pair set at that modulus.

## 24. Uniform squarefree counts in the construction

First, there are infinitely many admissible B. In [X,2X] restrict B to1 modulo20. The primes2 and5 cannot make B or B+1 nonsquarefree. For every other prime p, the two exclusions p^2|B or p^2|B+1 remove at most2X/(20p^2)+O(1) terms by CRT. Only p<=sqrt(2X+1) matter, and

\[
 2\sum_{p\ge3,\ p\ne5}p^{-2}
 \le2\sum_{n=3}^{\infty}n^{-2}
 \le2\left(\frac19+\int_3^\infty t^{-2}\,dt\right)
 =\frac89<1.
\tag{53}
\]

The total O(sqrt X) error is smaller than the remaining positive proportion of the progression. Hence arbitrarily large B satisfy the required two squarefree conditions.

We also need squarefree multipliers coprime to the varying Q. For an interval I contained in[1,X], squarefree inversion and coprimality inclusion-exclusion give

\[
 \begin{aligned}
 S(I,Q)&=\#\{n\in I:\mu^2(n)=1,(n,Q)=1\}\\
 &=|I|\frac{\varphi(Q)}Q\prod_{p\nmid Q}(1-p^{-2})
       +O\bigl(\sqrt X\,2^{\omega(Q)}+|I|/\sqrt X\bigr).
 \end{aligned}
\tag{54}
\]

Indeed, write n=d^2k, keep only (d,Q)=1, and count coprime k in the resulting interval as its length times phi(Q)/Q plus O(2^omega(Q)). There are at most sqrt X choices of d. Extending the absolutely convergent sum of mu(d)/d^2 to infinity costs at most O(|I|/sqrt X), uniformly in Q.

For I=I_B take X=2B^2. Here Q~B^17, so2^omega(Q)=B^o(1) by the elementary subpower estimate. Furthermore phi(Q)/Q>=2^-omega(Q), and the remaining Euler product is at least1/2, since the product over all integers n>=2 of1-n^-2 is1/2. The main term in (54) is therefore at least a constant times B^2/2^omega(Q), while the error is O(B2^omega(Q)). Their ratio tends to infinity. Thus

\[
                    S(I_B,Q)\ge B^{2-o(1)}
                                      =T^{2/3-o(1)}.
\tag{55}
\]

No distribution theorem for the prime factors of Q is assumed.

## 25. All witnesses, the exact long complement, and the scope of sharpness

For ANY positive primitive unit witness r/s at the primitive kernels a,b from (50), the integer a^3r^2-b^3s^2 is a NONZERO multiple of N. Nonvanishing follows from the distinct coprime squarefree kernels: equality would make (a/b)^3 a rational square, forcing a=b=1. If H=max(r,s), positivity gives

\[
                 N\le|a^3r^2-b^3s^2|\le b^3H^2,
 \qquad H\ge\sqrt{N/b^3}\asymp B^6=T^2.
\tag{56}
\]

Cancelling g^3 modulo N is valid for every constructed pair. Its complete witness set and minimum height are therefore those of the primitive ratio. The specified u/v has height<=4T^2, whereas EVERY witness has height bounded below by a positive constant times T^2. This exceeds floor(T^(5/4-delta)) for every fixed delta in(0,1/4] and all sufficiently large B. The whole subset is consequently in the EXACT long complement, not merely a set with some long witness tuples.

All its actual points s=x/y,t=v/u lie on the primitive vertical G(s,t)=(B+1)s-B. Its remainder is A=(B+1)s-B,B_remainder=0, its joint content is1, its primitive norm is [(B+1)s-B]^2 of OWN degree2, and its effective modulus is N. Its capacity is

\[
 2^{\omega(N)}\left(1+T^2/Q_2(N)\right)=T^{o(1)},
 \qquad Q_2(N)\ge\sqrt N.
\tag{57}
\]

The same subset is covered by U=u,V=v, a degree0 rational witness family with joint coefficient content1. These coefficients VARY with the modulus and scale. Thus (55) matches, up to subpower factors, the actual degree0 and degree-at-most1 upper exponent in (49), as well as the universal actual multiplier cap and per-vertical bound. This one-ratio construction does not assert sharpness of the quadratic-family exponent1; that separate conclusion is proved next.

Primality of the values N in (50) has not been established. The construction has only one ratio within its specified subset, not necessarily in the entire full set. Its T^(2/3-o(1)) pairs are far fewer than T^(2-o(1)), so it gives no counterexample to a full-K saving, and is consistent with K<<T*sqrt(R). All original squarefree, unit, interval, positivity, reduction, full-height and exact-long conditions have been checked directly.

The new contact and sharpness results leave the sufficient cover obligation in Section9, Open problem12 and the original full frequency targets open. They establish stronger bounds for specified witness families, not an economical cover of all kernels or an improvement below the prescribed-sum exponent3/10.

## 26. Matching quadratic counts in one fixed family

**Proposition 16 (two separate sharpness conclusions).** The quadratic actual-family exponent1 in (49) is optimal: one FIXED jointly primitive quadratic family has at least cT original admissible pairs for infinitely many N at T=N^(1/5). Separately, for every d>=0 the unit-root modulus exponent1/(2d+1) in Theorem13 cannot be increased while keeping a subpower number of covering classes, even for fixed coefficient examples. The second assertion is about root covers, not representability by actual bounded-height kernel fractions.

The quadratic construction is the lower example in [POWERFUL_SUM_GRAPHS.md](POWERFUL_SUM_GRAPHS.md); the verification below makes its fixed-family quantifier and all actual-witness gates explicit. Define

\[
 U(x,y)=-x^2+14xy+35y^2,\qquad
 V(x,y)=35x^2+14xy-y^2.
\]

They are jointly coefficient-primitive and direct expansion gives

\[
 x^3U(x,y)^2-y^3V(x,y)^2
             =(x-y)^5(x^2-23xy+y^2).
\tag{58}
\]

Fix r=21000 and q=r^2-23r+1=440517001. This positive integer is coprime to210. Also (r,q)=1 because q=1 modulo r, and (r-1,q)=1 because q=-21 modulo r-1 and (r-1,21)=1. For m>=5 set

\[
              L=2^m,\qquad N=qL^5/256,\qquad T=N^{1/5}.
\tag{59}
\]

These are positive integer moduli. Since q>256*16^5, we have T>16L. Choose y in[T,33T/32) satisfying y odd and (r-1)y=L modulo q, and put x=y+L. This is one progression modulo the FIXED integer2q, with T/(64q)+O(1) terms. Both kernels lie in[T,35T/32), are distinct and coprime, and x=r y modulo q. Since L,r,r-1 are units modulo q, both kernels are units modulo q; they are also odd.

At primes dividing2q neither kernel is divisible by p. At every other prime p, each condition p^2|y or p^2|y+L excludes at most T/(64q p^2)+O(1) terms of the progression. All relevant p are at most sqrt(35T/32). The excluded total is at most T/(64q) times2 sum_(p>=3)p^-2 plus O(sqrt T). As in (53), 2 sum_(p>=3)p^-2<=2 sum_(n>=3)n^-2<=8/9<1. As q is fixed, at least cT choices survive for all sufficiently large m, with c>0 independent of m. Their kernels are squarefree.

For these choices 1<x/y<17/16, and U,V are positive. In terms of L,

\[
        U=48y^2+12Ly-L^2,\qquad
        V=48y^2+84Ly+35L^2.
\]

Since L is divisible by32 and y is odd, both are16 modulo32. Write h=gcd(U,V)=16h_o, with h_o odd. No prime p|q divides either U or V. Indeed x,y are units modulo p, and (58) with x=r y modulo q implies that divisibility of either U or V by p forces divisibility of both. Then p divides U-V=36(y^2-x^2), so x/y=1 or-1 modulo p. The respective values of U/y^2 are48 and20, impossible because p|q is not2,3,5 or7. Thus (h_o,q)=1 as well.

Set u=U/h,v=V/h. These are positive, reduced, odd and coprime to q, hence units modulo N. Since x>=y and both forms are positive, U,V<=49x^2. Therefore

\[
 \max(u,v)\le\frac{49(35T/32)^2}{16}
                =\frac{60025}{16384}T^2<4T^2.
\tag{60}
\]

The remaining divisibility after this NONunit cancellation must also be checked. Write F=x^2-23xy+y^2. We know q|F. Also h_o^2 divides the left side of (58), and hence L^5 F. Since h_o is odd and L is a power of two, h_o^2|F. With (h_o,q)=1 we obtain qh_o^2|F. Dividing (58) by h^2=256h_o^2 now gives the integer identity

\[
 x^3u^2-y^3v^2
       =\frac{qL^5}{256}\,\frac{F}{qh_o^2}.
\tag{61}
\]

Thus N divides the required difference. This proves at least cT original squarefree, distinct, COPRIME unit kernel pairs in[T,2T), with positive reduced unit witnesses of the required height. All use the same fixed quadratic U,V, independently of N,T. It follows that a uniform O_epsilon(T^(1-delta+epsilon)) family bound with any delta>0 is impossible, by choosing epsilon<delta. No assertion that these particular pairs are in the exact long complement is needed or made.

## 27. A diagonal sequence attains the contact order

For d>=1 put M=2d+1 and define integer polynomials U_d,V_d by the odd/even decomposition

\[
 (w-1)^M(w^2+Mw+1)=V_d(w^2)-w^3U_d(w^2).
\tag{62}
\]

The coefficients of w and w^(M+1) on the left vanish. Its odd part is therefore divisible by w^3 and defines U_d of degree d, while its even part defines V_d of degree at most d. The leading coefficient of U_d and constant coefficient of V_d are both-1, so the pair is jointly primitive. This is the diagonal sequence in [ALL_DEGREE_DIAGONAL_PADE_KERNEL_RESTRICTION.md](ALL_DEGREE_DIAGONAL_PADE_KERNEL_RESTRICTION.md); only the identity is needed here.

Multiplying (62) by the same expression at-w gives directly

\[
 \begin{aligned}
 f_d(s)&=s^3U_d(s)^2-V_d(s)^2=(s-1)^M R_d(s),\\
 R_d(s)&=s^2-(M^2-2)s+1.
 \end{aligned}
\tag{63}
\]

In particular R_d(1)=4-M^2!=0. Also the left side of (62) is0 at w=1 and2^M(M-2) at w=-1. Hence

\[
                U_d(1)=V_d(1)=2^{M-1}(M-2)\ne0.
\tag{64}
\]

The nonzero unit root1 consequently has multiplicity EXACTLY M, with no common zero of U_d,V_d there. For d=0 take U_0=V_0=-1 and f_0=s^3-1.

Choose any fixed prime p>M+2. The integers R_d(1),U_d(1),V_d(1) in (63)--(64) are p-adic units. Thus on s=1 modulo p the exact valuation is

\[
                         v_p(f_d(s))=M v_p(s-1).
\tag{65}
\]

This uses a fixed coefficient family and a good prime; small-prime losses are irrelevant to the example.

## 28. Optimality with a subpower number of divisor classes

Fix d>=1 and the prime p from Section27, and let N=p^e with e increasing. In the disk1 modulo p, the unit roots modulo p^e are precisely the ball s=1 modulo p^a_e, where a_e=ceil(e/M). For every divisor p^k of N, that ball meets EXACTLY

\[
                           p^{\max\{0,k-a_e\}}
\tag{66}
\]

classes modulo p^k. Indeed the ball has p^(e-a_e) residues modulo p^e; fixing a finer residue modulo p^k leaves p^(e-k) when k>a_e.

Suppose a replacement for Theorem13 guaranteed a divisor Q>=c_d N^theta for theta>1/M and O_(d,epsilon)(N^epsilon) root classes for EVERY epsilon>0. Writing Q=p^k gives k>=theta e+O_d(1), so (66) forces at least

\[
                      N^{\theta-1/M}/O_d(1)
\]

classes. Choosing epsilon<theta-1/M contradicts that bound for large e. The same counting obstruction holds if different classes are assigned different divisor moduli, all at least p^k: each class contains at most p^(e-k) residues modulo N. For d=0, exponent1 is already the maximum possible because Q divides N.

This proves the second part of Proposition16. It asserts sharpness of the UNIT-root covering exponent. It does not assert that every root in these prime-power balls has an actual squarefree-kernel fraction of the specified height, or that every higher-degree pair bound in (48) is sharp. The quadratic pair-count sharpness in Section26 and the exact-long one-ratio sharpness in Sections23--25 are separate constructions, with their respective filters proved directly. None yields a nearquadratic lower bound for the full pair count or an economical full cover.


## 29. The genuine fixed-degree union
Fix d>=1, c,C>0 and N,T with T>=2 and cT^5<=N<=CT^5. Count actual pairs x,y which are distinct squarefree units modulo N in[T,2T), and admit a positive reduced unit witness u/v of height at most4T^2 and x^3u^2=y^3v^2 modulo N.

Let E_d(N,T) restrict to those actual pairs for which the witness comes from a jointly coefficient-primitive homogeneous degree-d integer pair(U,V), polynomially coprime after dehomogenization, whose norm

    f(s)=s^3U(s)^2-V(s)^2

has a NONZERO characteristic-zero root of multiplicity M=2d+1. The coefficients and family may be chosen separately for every pair.

**Theorem 17 (genuine fixed-degree maximal-contact unions).** With arbitrary added epsilon,

    E_1 << T^(3/2+epsilon),
    E_2 << T^(7/5+epsilon),
    E_3 << T^(4/3+epsilon),
    E_4 << T^(11/9+epsilon),
    E_5 << T^(109/99+epsilon),
    E_6 << T^(37/39+epsilon).                               (67)

For EVERY fixed d>=7, E_d(N,T)=0 for all sufficiently large T, uniformly in N in the balanced range. The threshold and implied constants depend on d,c,C. Theorem21 below removes any advance bound on the degree; the individual fixed-degree estimates are retained here.

## 30. The fixed canonical sequence and its elementary properties
Recall the integer forms from Section27. For d>=1 put M=2d+1 and define U_d,V_d in Z[s] by

    H(w)=(w-1)^M(w^2+Mw+1)=V_d(w^2)-w^3U_d(w^2).           (68)

The w and w^(M+1) coefficients vanish. Thus deg U_d=d, deg V_d<=d, and

    leadingcoeff(U_d)=-1,  V_d(0)=-1.

The same notation denotes their degree-d homogeneous forms.

Multiplying(68) by its version at-w gives

    s^3 U_d(s)^2-V_d(s)^2
       =(s-1)^M [s^2-(M^2-2)s+1].                           (69)

Also

    U_d(1)=V_d(1)=2^(M-1)(M-2) !=0.                         (70)

The pair is coprime over Qbar, including projective endpoints. Indeed a common finite root beta would be a root of(69). Beta=1 is impossible by(70). Beta=0 is impossible because V_d(0)=-1. The residual quadratic has distinct roots, since its discriminant is M^2(M^2-4)!=0; a common root of U_d,V_d there would have norm multiplicity at least2, contradicting its simple multiplicity in(69). At infinity the leading coefficient of U_d is nonzero. Consequently the homogeneous forms have no common projective zero.

## 31. Classification over every rational centre
Suppose a coprime degree-at-most-d pair U,V has a nonzero norm root alpha of multiplicity M. Since deg f<=2d+3 and d>=1,

    M [Q(alpha):Q] <=2d+3<2M,

so alpha is rational. Coprimality makes U(alpha),V(alpha) nonzero. Their norm equality implies

    alpha=(V(alpha)/(alpha U(alpha)))^2>0.

Write alpha=(R/S)^2, with R,S positive coprime integers, and change the sign of V if needed. Define

    A_(R,S)(a,b)=S^3 U_d(S^2a,R^2b),
    B_(R,S)(a,b)=R^3 V_d(S^2a,R^2b).                        (71)

These are integer forms of degree d. Their joint content is1 because the a^d coefficient of A is-S^(2d+3), and the b^d coefficient of B is-R^(2d+3). They remain polynomially coprime. At alpha their norm has contact exactly M by(69), and their ratio B/A equals(R/S)^3 by(70).

Choose the local branch phi(s)=s^(3/2) at alpha with value(R/S)^3. Both V-phi U and B-phi A vanish to order at least M. Their crossnumerator

    V(s)A(s,1)-B(s,1)U(s)

has degree at most2d and order at least2d+1, so is identically zero. As A,B are coprime and of maximal degree d, every pair of degrees at most d representing the same rational function is a common scalar multiple. Joint integer primitivity makes that scalar+1 or-1. A separate sign on V is harmless for absolute witness heights and squared congruences. This proves that(71) lists ALL genuine families in E_d.

This argument FAILS for standalone d=0: a simple root of s^3U^2-V^2 need not be rational. The union over arbitrary constant witness families is not being estimated by this classification.

## 32. Fixed-degree uniform evaluation gcd
Write x=g a,y=g b with g=gcd(x,y). Then a,b are positive coprime squarefree in[Y,2Y), Y=T/g. Put X=S^2a,Z=R^2b and evaluate A=S^3 U_d(X,Z), B=R^3V_d(X,Z). Let h=gcd(A,B)>0.

If p|R, then p does not divide S. Write alpha_p=v_p(a), beta_p=v_p(b), each0 or1. We have v_p(X)=alpha_p, while v_p(Z)=2v_p(R)+beta_p>alpha_p. The leading-X^d term of U_d has strictly smaller valuation than every other term. Hence

    v_p(A)=d alpha_p.

Every term of V_d has valuation at least d alpha_p, and its prefactor R^3 makes v_p(B) strictly larger. Thus v_p(h)=d v_p(a). Similarly v_p(h)=d v_p(b) for p|S. These statements include small primes and arbitrary powers in R,S.

At p not dividing RS, X,Z are primitive at p. The homogeneous pair U_d,V_d has no common projective zero, so there are fixed integer Bezout identities with right sides C_d X^L and C_d Z^L for some nonzero integer C_d and integer L, depending only on d. One elementary construction applies the univariate Euclidean algorithm over Q, clears denominators, homogenizes, and repeats for reciprocal polynomials. At least one of X,Z is a p-unit, so

    v_p(h)<=v_p(C_d)  for p not dividing RS.

After enlarging C_d to a positive integer, this gives the UNIFORM bound

    h <= C_d gcd(a,rad R)^d gcd(b,rad S)^d
      <= C_d R^d S^d.                                      (72)

It is crucial that a,b are squarefree before deriving(72). With arbitrary nonsquarefree primitive kernels their valuations at primes of R,S would not be bounded by1.

In the actual unit-witness set one also has(RS,N)=1: if p|R and p|N, a,b are units and v_p(A)=0<v_p(B), so the reduced witness cannot have two unit coordinates. The case p|S is symmetric. Discarding this extra centre restriction only enlarges the upper bound.

## 33. Archimedean centre-height restriction
Let B_c=max(R,S). There is c_d>0, depending only on d, such that for a,b in[Y,2Y),

    max(|A|,|B|)>=c_d B_c^(2d+3)Y^d.                        (73)

Here is a direct compactness proof with its boundary checked. Divide R,S by B_c and a,b by Y. The resulting domain is compact: a/Y,b/Y in[1,2], R/B_c,S/B_c in[0,1], with their maximum1. The normalized evaluated pair is continuous. If both normalized R,S are positive, simultaneous vanishing is excluded by the projective coprimality of U_d,V_d. At R=0,S=1, the first coordinate is U_d(a,0)=-a^d, nonzero; at S=0,R=1 the second is V_d(0,b)=-b^d. Therefore the maximum has a strictly positive minimum on the entire compact domain. Homogeneity gives(73).

Divide(73) by(72). The reduced witness height H obeys

    H>=c'_d B_c^3Y^d.

Since H<=4T^2, we obtain

    B_c^3 <<_d T^(2-d)g^d.                                 (74)

For fixed g the number of possible positive primitive centres is consequently

    J_g <<_d T^((4-2d)/3)g^(2d/3).                          (75)

There is no hidden plus-one loss in(75): if the upper height for B_c is less than1 there are NO integer centres; otherwise the number of pairs R,S is bounded by a constant times that height squared.

The independent actual height cap is g<<_(c) T^(2/3). It follows because a^3u^2-b^3v^2 is a nonzero multiple of N, of absolute value at most128T^7/g^3. Nonvanishing uses distinct squarefree primitive kernels; canceling g^3 modulo N uses g being a unit.

For d>2, (74) and B_c>=1 require g>>_d T^(1-2/d). When d>=7 is FIXED, this contradicts g<<T^(2/3) for all sufficiently large T. This proves the asserted emptiness. The d-dependent constants prevent any inference of an unbounded-degree threshold from this argument alone.

## 34. Uniform root packing and the six exponents
For a fixed d and centre, the coefficient-uniform structured unit-root theorem gives O_(d,epsilon)(N^epsilon) root classes modulo a coefficient-independent divisor Q>=c_d N^(1/(2d+1)). For primitive kernel fractions of height Y=T/g, angular packing yields at most

    O_(d,epsilon)(N^epsilon[1+T^(2-5/(2d+1))/g^2])            (76)

possible pairs a,b. Actual reduced-witness congruences imply this necessary norm congruence by MULTIPLYING their pointwise gcds; none is inverted modulo N.

For d=1, multiply(76) by(75), sum g<=G, and use the unconditional tail T^2/G. This gives

    E_1 << N^epsilon[T^(2/3)G^(5/3)+T]+T^2/G.

Taking G=floor(T^(1/2)) proves exponent3/2.

For d=2 the same calculation gives

    E_2 << N^epsilon[G^(7/3)+T G^(1/3)]+T^2/G.

Taking G=floor(T^(3/5)) proves exponent7/5.

For3<=d<=6 use the full actual cap g<=G_max with G_max=O(T^(2/3)). The two sums from(75)--(76) are

    N^epsilon T^((4-2d)/3)
       sum_(g<=G_max) [g^(2d/3)+T^(2-5/(2d+1))g^(2d/3-2)].

Both summation exponents are at least0 when d>=3, so they are bounded by the corresponding endpoint powers times G_max. The resulting T exponents are

    A_d=2-2d/9,
    B_d=8/3-2d/9-5/(2d+1).                                 (77)

Taking their maximum gives the exact table

    d    A_d       B_d           max
    3    4/3       9/7           4/3
    4    10/9      11/9          11/9
    5    8/9       109/99        109/99
    6    2/3       37/39         37/39.

After renaming epsilon using N comparable to T^5, these are(67). The lower restriction on g from(74) may be discarded for this upper bound. Coefficient uniformity is indispensable here because the centre varies with g,N,T and the pair.

## 35. Rational cube witnesses and prime-power exclusion

**Theorem 18 (the cube-witness class).** Let C_N(T) count original actual pairs admitting a reduced positive rational-cube witness. Uniformly over balanced integer moduli, C_N(T)<<_(epsilon,c,C)T^(2/3+epsilon). For all sufficiently large balanced PRIME-POWER moduli, C_N(T)=0. The threshold may depend on the fixed balance constants.

The following argument is independent of the rational-centre classification, so it also applies to the cube case left by polynomial cancellation.

If u/v=(R/S)^3 in lowest terms, the witness-height condition gives max(R,S)<=4^(1/3)T^(2/3). Define the positive rational number

    z=xR^2/(yS^2)=A/B in lowest terms.

Its reduced numerator and denominator are units modulo N, and both are at most

    Z=2*4^(2/3)T^(7/3).

Moreover z^3=1 modulo N. There are at most3^omega(N) unit cube-root residues. For p!=3 this follows by the degree3 root bound and unique lifting, since the derivative is a unit at a cube root (also for p=2). For p=3 a root reduces to1 modulo3. If a=v_3(z-1)>=1, direct expansion gives v_3(z^3-1)=a+1, so at most3 residues modulo3^e can occur. The e=1 case is immediate. CRT proves the bound.

For sufficiently large T, 2Z^2<N because Z^2 is of order T^(14/3), whereas N is bounded below by a constant times T^5. Two DISTINCT reduced positive fractions of height at most Z cannot lie in the same unit residue modulo N: their nonzero crossdeterminant is divisible by N but has absolute value at most2Z^2. Thus only3^omega(N)=T^o(1) rational values of z occur.

For each fixed z, its prime-valuation PARITIES determine the product ab for the primitive squarefree kernel ratio x/y=a/b. They do NOT determine which primes go into a and which into b. Every prime with odd v_p(z) is assigned to exactly one of a,b; primes with even valuation cannot occur in either. Therefore there are at most

    2^omega(A B)=T^o(1)

possible primitive kernel ratios per z. This is the necessary orientation factor; no uniqueness is asserted.

For each such ratio the original full-height actual gcd cap gives at most O(T^(2/3)) multipliers g. Indeed the nonzero integer a^3u^2-b^3v^2 is divisible by N and bounded by128T^7/g^3. Combining the subpower bounds proves the all-modulus estimate in Theorem18. No primality or squarefree-modulus assumption was used.

### No cube witnesses at large prime powers

Suppose N=p^e and keep z=A/B as above, with(A,B)=1. The ratio z is not1: otherwise xR^2=yS^2, and parity at every prime would make the distinct squarefree kernels equal. Thus A-B!=0 and

    p^e | (A-B)(A^2+AB+B^2).

The gcd of the two factors divides3, since reducing the second modulo A-B gives3B^2 and(A-B,B)=1.

If p!=3, the entire p^e divides one of the two factors. For sufficiently large T, however,

    0<|A-B|<=Z<N,
    0<A^2+AB+B^2<=3Z^2<N,

which is impossible. If p=3, the congruence first forces A=B modulo3. Since A,B are3-adic units, expanding at A=B shows

    v_3(A^2+AB+B^2)=1.

Hence3^(e-1) divides A-B. But0<|A-B|<=Z<N/3 for sufficiently large T, another contradiction. This proves the prime-power assertion. It is consistent with composite-modulus lower constructions; it is not a statement about all noncube witnesses.

## 36. Sharp cube pairs in the exact long complement

**Proposition 19 (sharp cube-witness exponent).** Infinitely many composite moduli N at the literal scale T=N^(1/5) admit a subset of at least T^(2/3-o(1)) original actual cube-witness pairs with one primitive kernel ratio. Every pair in the subset belongs to the exact NO-short/SOME-full complement for every fixed delta in(0,1/4]. Thus the exponent2/3 in Theorem18 is optimal over balanced integer moduli. No prime-modulus or full-ratio-count assertion is made.

### Parameters, positivity and units

Take arbitrarily large B=1 modulo20 for which B and B+1 are squarefree. Infinitude follows from the elementary progression argument in Section24. Put

    a=B, b=B+1, S=B^2+1, R=2B^2+3=2S+1,
    u=R^3, v=S^3,
    N=a^3R^6-b^3S^6, T=N^(1/5).

The leading term of N is63B^15. Hence N~63B^15, T~63^(1/5)B^3, and N is a positive integer for all sufficiently large B. In fact R/S>2 and (b/a)^3<4 for B>=2, so N>0 there.

The exact gcd identities are

    (a,b)=1, (R,S)=1, (a,S)=1,
    (b,R)=(B+1,5)=1.

The last equality uses B+1=2 modulo5. Reducing N modulo each of a,b,R,S proves (abRS,N)=1. For example N modulo R is -b^3S^6, a unit modulo R, and the other three cases are identical applications of the displayed cross-gcd identities. Consequently u,v are positive, coprime and units modulo N. Their quotient is the rational cube (R/S)^3 in lowest terms.

The height gate holds with room to spare:

    max(u,v)=R^3~8B^6,
    4T^2~4*63^(2/5)B^6>8B^6.

Thus max(u,v)<=4T^2 for all sufficiently large admissible B. The constants here are absolute.

This modulus is explicitly composite. If A=aR^2 and C=bS^2, then

    N=A^3-C^3=(A-C)(A^2+AC+C^2),
    A-C~3B^5, A^2+AC+C^2~21B^10.

Both positive factors exceed1 for all sufficiently large B. No assertion of prime-modulus sharpness is made.

### Original admissibility conditions

Let Q=abN and use the interval

    I_B=[T/a,2T/b).

It has length T(B-1)/(B(B+1))~63^(1/5)B^2 and lies in[1,6B^2] for all sufficiently large B. For every squarefree g in this interval with (g,Q)=1 put x=ga,y=gb. Then x,y are squarefree, distinct and units modulo N; their gcd is exactly g, and

    T<=x<y<2T,
    x/y=B/(B+1),
    x^3u^2-y^3v^2=g^3N.

Thus the same original positive reduced unit witness u/v=(R/S)^3 of height<=4T^2 works for every such pair. Distinct multipliers give distinct ordered pairs. The subset has one primitive ratio; the entire set at N need not have only one.

### The elementary lower count

Section24 proves that arbitrarily large B=1 modulo20 have B,B+1 squarefree, and equation(54) supplies the squarefree-coprime interval sieve uniformly in Q. Here Q=abN is of order B^17, the multiplier interval has length of order B^2 and lies in[1,6B^2]. Its main term is at least cB^2/2^omega(Q), while the error is O(B*2^omega(Q)). Since2^omega(Q)=B^o(1), their ratio tends to infinity. There are therefore at least B^(2-o(1))=T^(2/3-o(1)) admissible multipliers. Their gcds have order B^2, equivalently T^(2/3). The constants in the interval and Q estimates are absolute; no new uniformity in an independent parameter is used.

### Every witness is long

Let r/s be ANY positive primitive unit witness for a,b at this N, and H=max(r,s). The integer a^3r^2-b^3s^2 is a nonzero multiple of N. Nonvanishing follows from distinct coprime squarefree a,b: equality would make a/b a rational square, forcing a=b=1. Positivity gives

    N<=|a^3r^2-b^3s^2|<=b^3 H^2,
    H>=sqrt(N/b^3)~sqrt(63)B^6~constant*T^2.

For each constructed pair the multiplier g is a unit modulo N, so cancelling g^3 shows that its ENTIRE witness set equals the primitive ratio's witness set. This lower bound therefore excludes ALL witnesses of height<=floor(T^(5/4-delta)) for every fixed delta in(0,1/4] and all sufficiently large B. The explicit u/v supplies SOME full-height witness. Every constructed pair lies in the exact NO-short/SOME-full complement.

Together with Theorem18, this lower construction establishes the sharp power exponent2/3 for that class. The upper proof is separate: this construction alone proves the stated lower and exact-long membership, without a full-K estimate or an economical family cover.

## 37. Polynomial common factors and bounded presentation degrees
Now fix an ORIGINAL presentation degree d>=1, allow U,V of degrees at most d with common polynomial factors, and suppose the norm has a nonzero root alpha of multiplicity2d+1. The original degree/Galois argument already forces alpha to be rational.

The orders of U,V at alpha must agree, say r: otherwise the norm order is twice their smaller order and at most2d. Divide both by(s-alpha)^r. Put j=d-r. The residual polynomials are nonzero at alpha and their norm has contact exactly

    2d+1-2r=2j+1.

At an evaluated point defining an actual nonzero witness, every canceled common factor is nonzero. Cancellation over Q therefore preserves its reduced witness fraction; no modular inverse is used.

If j>=1 they are globally polynomially coprime. Any further common polynomial factor, which cannot vanish at alpha, would reduce their degree bound below j without reducing the contact, contradicting the characteristic-zero contact bound2(degree)+1. Their actual maximal degree is j for the same reason. Thus their rational witness belongs to the genuine degree-j class counted in(67), or is eventually excluded when fixed j>=7.

If j=0, the residual pair is constant. As alpha is ALREADY rational and its residual norm vanishes there, alpha is a rational square and the absolute constant witness ratio is a rational cube. This is the only constant conclusion needed here; it is not a statement about all standalone degree0 families.

The residual cube-witness class is covered by Theorem18, which was proved without using this classification.

**Corollary 20 (bounded presentation degrees).** For a fixed finite range of ORIGINAL presentation degrees1,...,D, take the union over the possible reduced degrees j and the cube case. The largest exponent is3/2. The union is O_(D,epsilon,c,C)(T^(3/2+epsilon)). A large unreduced degree can still encode a lower-degree family by common-factor padding, so only genuinely reduced degrees>=7 are eventually empty.

## 38. Uniform cutoff with no advance degree bound

**Theorem 21 (unbounded exact-maximal-contact presentations).** There is a threshold T_0(c), independent of every degree, centre, family coefficient and balanced modulus, beyond which every genuine actual maximal-contact witness of degree d>=1 has d<=6. Consequently the union over ALL genuine degrees d>=1 is O_(epsilon,c,C)(T^(3/2+epsilon)); the genuine union over ALL d>=2 is O_(epsilon,c,C)(T^(7/5+epsilon)). Allowing ALL positive presentation degrees without bound and permitting polynomial common factors still gives O_(epsilon,c,C)(T^(3/2+epsilon)). Standalone degree0 is excluded.

The classification in Section31 is valid for each integer d without any quantitative constant. To make the cutoff uniform in d, we now replace the fixed-degree resultant bound in Section32 by an evaluation-dependent bound and replace the compactness minimum in Section33 by an explicit exponential estimate.

### A degree-uniform evaluation-gcd bound

For each d>=1 put M=2d+1. Write H_M(w) for H(w) in(68), and define

    Q_M(X,Z)=X^2-(M^2-2)XZ+Z^2.

Write x=g a,y=g b with g=gcd(x,y). Then a,b are positive coprime squarefree in[Y,2Y), Y=T/g. Set

    X=S^2a, Z=R^2b, B_c=max(R,S),
    A=S^3U_d(X,Z), B=R^3V_d(X,Z), h=gcd(A,B),
    D_0=gcd(a,rad R)gcd(b,rad S).

We prove the EVALUATION-DEPENDENT bound

    h <= 2^M |Q_M(X,Z)| D_0^d.                              (78)

This is the step which removes uncontrolled degree-dependent resultant constants from the large-degree argument.

At a prime p|R, S is a unit and v_p(X)=v_p(a) is0 or1, whereas v_p(Z)>=2. The leading-X^d term of U_d has strictly smaller valuation than every other term. Thus v_p(A)=d v_p(a). Every term of V_d has valuation at least d v_p(a), and the R^3 prefactor raises it strictly. Hence v_p(h)=d v_p(a). Similarly v_p(h)=d v_p(b) for p|S. These contributions are exactly absorbed by D_0^d, at all primes including2.

Now take p not dividing RS. If p divides a or b, one of U_d(X,Z),V_d(X,Z) is a unit by its endpoint coefficient-1, so p does not divide h. We therefore need only X,Z p-adic units. Extend the normalized p-valuation to a field containing w with w^2=X/Z. This w is a unit. Homogeneity and the R,S unit factors show that both U_d(w^2),V_d(w^2) have valuation at least e=v_p(h). Thus both H_M(w),H_M(-w) have valuation at least e.

For odd p, at least one of w-1,w+1 is a unit. If w-1 is a unit, the first factorization in(68) gives e<=v(w^2+Mw+1). If w+1 is a unit, use H_M(-w) to get e<=v(w^2-Mw+1). Either bound is at most

    v((w^2+Mw+1)(w^2-Mw+1))=v_p(Q_M(X,Z)),

because Z is a p-unit and both factors are integral. THERE IS NO EXTRA FACTOR M in this estimate.

For p=2 outside RS, X/Z reduces to1, so w reduces to1 in residue characteristic2. Both w^2+Mw+1 and w^2-Mw+1 are units because M is odd. Therefore

    e<=M min(v(w-1),v(w+1))<=M,

since their difference is2. This contributes at most2^M. Combining every prime proves(78), in fact as a divisibility bound up to the sign of Q_M.

The residual Q_M(X,Z) is nonzero for positive rational X/Z. Otherwise its discriminant M^2(M^2-4) would be a rational square. But M^2-4 lies strictly between(M-1)^2 and M^2 for M>=3. Hence its absolute value in(78) is a positive integer.

For comparison, the older coefficient-independent fixed-degree cancellation bound is h<=C_d D_0^d, with C_d depending only on d; the explicit old bound for primitive arguments is2^M M(M-2)^(M/2). That bound alone grows too rapidly to extrapolate a fixed-degree argument. Bound(78) uses the actual quadratic residual instead, and is uniform without such an extrapolation.

### An explicit real lower bound

For every t>0 and odd M>=3,

    max(|H_M(t)|,|H_M(-t)|)
           >=2^(-M) max(1,t)^(M+2).                         (79)

It suffices by reciprocity to treat0<t<=1: the absolute values at1/t are t^(-M-2) times those at t. If t<=1/2, H_M(t) has magnitude at least(1-t)^M>=2^-M. If1/2<=t<=1, then

    Mt-t^2-1>=3t-t^2-1>=1/4,

so |H_M(-t)|>=(1/4)(1+t)^M>=2^-M. This proves(79).

Put t=sqrt(X/Z)=(S/R)sqrt(a/b), and q=(a/b)^(3/2). The EXACT normalization is

    R^3 Z^d H_M(t)=B-qA,
    R^3 Z^d H_M(-t)=B+qA.

Since a/b lies between1/2 and2, q<=2sqrt2 and1+q<4. Also

    R max(1,t)=max(R,S sqrt(a/b))>=B_c/sqrt2,
    R^3 Z^d=R^(2d+3)b^d.

Using b>=Y and(79) gives

    max(|A|,|B|)>=2^(-3d-5) B_c^(2d+3)Y^d.                 (80)

Thus all archimedean constants are explicit and exponential in d; none comes from an unquantified minimum as d varies.

### A universal degree cutoff

The simple size bounds are

    |Q_M(X,Z)|<=4M^2 B_c^4Y^2,
    D_0<=min(RS,ab)<=min(B_c^2,4Y^2).

Combine these with(78),(80). The REDUCED witness height satisfies

    H=max(|A|,|B|)/h
      >= B_c^(2d-1)Y^(d-2)
          /[2^(5d+8) M^2 min(B_c^2,4Y^2)^d].              (81)

For d>=1 the centre-dependent factor obeys

    B_c^(2d-1)/min(B_c^2,4Y^2)^d >=1/(2Y).

If B_c<=2Y the left side is1/B_c; if B_c>=2Y it is increasing in B_c and its value at2Y is1/(2Y). Consequently

    H >= L_d(Y):=Y^(d-3)/[2^(5d+9)(2d+1)^2].               (82)

This bound is uniform in d, R,S, a,b.

The actual conditions independently give

    g<=(128/c)^(1/3)T^(2/3),
    Y>=alpha_c T^(1/3),  alpha_c=(c/128)^(1/3).              (83)

Indeed a^3u^2-b^3v^2 is a NONZERO multiple of N, of absolute value at most128T^7/g^3. Nonvanishing follows from distinct coprime squarefree kernels, and cancellation modulo N uses the original unit conditions.

For Y>=128, the right side of(82) is increasing in every integer d>=1, because

    L_(d+1)(Y)/L_d(Y)
       =(Y/32)((2d+1)/(2d+3))^2 >=36/25>1.

Thus if d>=10,

    H>=L_10(Y)>=alpha_c^7 T^(7/3)/(2^59*21^2).

For all sufficiently large T depending only on c this exceeds4T^2. Therefore ALL d>=10 are excluded SIMULTANEOUSLY, with no advance bound on d or centre heights.

It remains to exclude7,8,9. Here only fixed-degree constants are needed. Since U_d,V_d have no common projective zero, fixed integer Bezout identities with right sides C_d X^L,C_d Z^L give h<=C_d D_0^d. Together with(80) and D_0<=B_c^2,

    H>=c_d B_c^3Y^d>=c_dY^d.

For EACH d=7,8,9, equation(83) then gives a constant times T^(d/3), which exceeds4T^2 eventually. Taking the maximum of these THREE fixed thresholds and the uniform d>=10 threshold produces a single T_0(c), independent of every allowed degree. This proves the uniform cutoff in Theorem21.

### The unbounded unions

Theorem17 supplies the six explicit genuine-degree estimates with exponents3/2,7/5,4/3,11/9,109/99,37/39. After the uniform cutoff, only this FIXED finite set of six degrees remains. Their largest exponent is3/2, and after excluding degree1 it is7/5. Thus their implied constants can be absorbed into constants depending only on epsilon and the fixed balance constants.

For any original positive presentation degree, the algebraic reduction in Section37 either gives a genuine residual degree j>=1 or a cube witness. The reduction uses no degree-dependent quantitative constant: the original degree/Galois argument forces its centre rational, and the characteristic-zero contact bound excludes all remaining common factors. For T>T_0(c), Theorem21's uniform cutoff therefore leaves j<=6, regardless of how large the original padded presentation degree was. Theorem18 handles the residual cube class with exponent2/3. Their union gives the asserted unbounded-presentation estimate. For bounded T, the trivial O(T^2) pair count is absorbed by enlarging the implied constant using the fixed T_0(c); this creates no presentation-degree dependence.

The argument applies only to exact MAXIMAL characteristic-zero contact. It does not assert that high modular contact is exact, supply an economical cover of unrestricted witnesses, or yield a bound for the full K_N(T). All actual distinct-squarefree, positivity, reduction, unit, height and balance gates were retained. Arbitrary standalone constant witness families remain outside this theorem.

## 39. Projective rational reconstruction with audited cancellation

**Proposition 22 (retained-modulus projective reconstruction).** The congruence and height gates below imply equality of the reduced projective points over Q. The pointwise evaluation gcds are charged explicitly; none is inverted modulo the original modulus.

Let m>=1. Let nonzero integer pairs(A,B),(A_0,B_0) satisfy

    A=lambda A_0 modm,  B=lambda B_0 modm,  (lambda,m)=1.

Let h=gcd(A,B)>0 and h_0=gcd(A_0,B_0)>0. Put a=A/h,b=B/h and a_0=A_0/h_0,b_0=B_0/h_0, allowing one reference coordinate to be zero. Then

    m_ret=m/gcd(m,h_0) divides a b_0-a_0 b.                  (84)

This sharper loss, rather than a silent inversion of h or h_0, follows prime by prime. For p^e||m write ell_0=min(v_p(A_0),v_p(B_0)). If ell_0>=e, (84) imposes no p-condition. If ell_0<e, the unit lambda implies min(v_p(A),v_p(B))=ell_0. Moreover

    AB_0-A_0B=(A-lambda A_0)B_0-(B-lambda B_0)A_0

has valuation at least e+ell_0. Dividing by hh_0 loses exactly2ell_0 at p and leaves at least e-ell_0, proving(84). No congruence for the reduced fractions modulo the full m is asserted.

If max(|a|,|b|)<=H and max(|a_0|,|b_0|)<=H_0, then

    |a b_0-a_0 b|<=2H H_0.

Consequently

    m_ret>2H H_0                                             (85)

forces the reduced projective points to be EQUAL over Q. If the actual point has positive finite ratio but the reference has zero, infinity or negative ratio, (85) excludes that actual point. This also handles zero reference evaluations without dividing by zero.

The bound is sharp in its p-adic loss: taking 0<=ell<e and A_0=p^ell a_0,B_0=p^ell b_0, A=A_0+p^e c,B=B_0 with b_0 a p-unit gives reduced crossdeterminant exactly of valuation e-ell for suitable c. Thus a coefficient-congruence-only claim retaining full m is generally false.

## 40. Explicit quadratic and linear coefficient tubes

**Theorem 23 (actual coefficient-tube reconstruction).** Under the original actual witness conditions, the sufficient retained-modulus gates below identify the moving witness with an exact quadratic or linear maximal-contact reference. The moving coefficients need no archimedean bound, and the reference may depend on the pair.

### Quadratic references

Let original kernels x,y be positive distinct squarefree units modulo N in[T,2T), and write x=g a,y=g b, so a,b are primitive squarefree and in[Y,2Y), Y=T/g. Let an actual reduced positive unit witness of height H<=4T^2 be obtained by evaluating a jointly primitive homogeneous quadratic pair(U,V) at(a,b).

Suppose, for some divisor m|N, its coefficient vector is congruent modulo m, up to a unit scalar lambda, to one of the exact reference pairs

    A_0=S^3 U_2(S^2a,R^2b),
    B_0=+/- R^3 V_2(S^2a,R^2b),
    U_2(X,Y)=-X^2+14XY+35Y^2,
    V_2(X,Y)=35X^2+14XY-Y^2,

with positive coprime R,S. The scalar must be a UNIT modulo m. This coefficient condition implies the evaluated pair congruences used in Section39, but does not itself justify canceling their pointwise gcds. All projective comparisons use the signed evaluations; equality also identifies their positive absolute witness fractions.

The actual unit filter forces(RS,m)=1. For if p|R and p|m, then a,b,S are p-units, A_0 is a p-unit, and B_0 is divisible by p. The same is true of the moving evaluations modulo p, so their positive reduced witness cannot have both entries units modulo N. The p|S case is symmetric.

The reference evaluation gcd then satisfies

    gcd(m,h_0)<=17280.                                      (86)

The exact identities

    (497Y-245X)U_2+(115Y-7X)V_2=17280Y^3,
    (115X-7Y)U_2+(497X-245Y)V_2=17280X^3

show, away from RS, that v_p(h_0)<=v_p(17280), because a,b are primitive. The prefactors R^3,S^3 are p-units there. All primes of m are away from RS. This uses a FIXED integer loss and is uniform in R,S.

Let B_c=max(R,S). The reference reduced height has the elementary upper bound

    H_0<=200 B_c^7 Y^2,                                     (87)

obtained by summing the absolute coefficients of the two quadratic forms. It follows from (84)--(87) that the sufficient condition

    m/17280>1600 B_c^7 T^4/g^2                              (88)

forces the moving actual witness to equal the exact reference witness over Q. A sharper pointwise test is simply m/gcd(m,h_0)>8T^2 H_0; it may cover references much larger than the crude bound(87).

For m=N and N>=cT^5, one explicit sufficient centre range is

    B_c^7 < [c/27648000] T g^2.                             (89)

The constants are deliberately conservative; the point is the retained modulus and the powers, not an optimized endpoint.

Thus an arbitrary moving coefficient family congruent to a low-height exact quadratic reference modulo N contributes NO NEW actual witnesses beyond that exact reference class. Adding enormous multiples of N to its coefficients cannot evade the low-height reconstruction argument. The moving coefficients need no archimedean height bound here.

### Linear references

For the exact linear references

    A_0=S^3(5R^2 b-S^2 a),
    B_0=+/-R^3(5S^2 a-R^2 b),

the same proof gives(RS,m)=1. The identities

    (5Y-X)+5(5X-Y)=24X,
    5(5Y-X)+(5X-Y)=24Y

give gcd(m,h_0)<=24. The reference reduced height is at most12 B_c^5Y. Therefore

    m/24>96 B_c^5 T^3/g                                    (90)

forces equality of moving and reference actual witnesses. With m=N>=cT^5 a sufficient centre range is

    B_c^5 < [c/2304] T^2 g.                                 (91)

Again, this is a quantitative coefficient-tube gate, not a statement that arbitrary high-contact modular roots have such rational centres.

## 41. Restricted tube unions and the remaining covering obligation

**Corollary 24 (tube-union bounds).**

The proved restricted-union bounds now apply automatically to every actual pair meeting these tube hypotheses: all quadratic-reference cases satisfying(88) lie in the exact genuine degree2 union of Theorem17 with count O(T^(7/5+epsilon)); all linear-reference cases satisfying(90) lie in the genuine degree1 union of Theorem17 with count O(T^(3/2+epsilon)). One may choose the reference and moving family separately for each pair; reconstruction identifies the ACTUAL witness value, so no sum over moving coefficient vectors is paid.

This is not a bound for the number of all near-contact auxiliary curves. The missing covering theorem would need to show that occupied moving families either satisfy an audited tube condition like(85), or have a positive-saving total primitive-ratio capacity by some other quantitative transversality argument. Merely having nonzero characteristic-zero discriminant does not do this. Merely having a coefficient congruence modulo m does not do this without the unit scalar, the evaluation-gcd loss and the reduced-height inequality.

When the available congruence is only a raw crossdeterminant AB_0-A_0B=0 modm rather than the stronger vector congruence up to a unit, the always valid weaker retained modulus is m/gcd(m,hh_0). It too can be used in(85), but its potentially larger loss must be charged. No pointwise common factor is inverted in this note.


## 42. A fixed-centre contact-four pencil and an exact full-count reverse

The following pencil shows that a fixed nonmaximal contact locus already contains the full remaining witness problem, up to two fixed maximal-contact families.

### The pencil

Define the homogeneous degree2 pairs by dehomogenization

    U_1(s)=4s+20,          V_1(s)=15s^2+10s-1,
    U_2(s)=s^2-2s+25,     V_2(s)=10s^2+16s-2.

Thus, for example, U_1(a,b)=4ab+20b^2. For a nonzero integer parameter pair(lambda,mu) set

    U=lambda U_1+mu U_2,  V=lambda V_1+mu V_2.

Direct expansion gives

    U_1 V_2-U_2 V_1=-15(s-1)^4.                             (92)

The norm has the exact factorization

    s^3U^2-V^2=(s-1)^4 R_(lambda,mu)(s),                    (93)

where

    R=16lambda^2 s-lambda^2
       +8lambda mu s^2+56lambda mu s-4lambda mu
       +mu^2 s^3+48mu^2 s-4mu^2,
    R(1)=15(lambda+mu)(lambda+3mu).                          (94)

These identities are polynomial equalities over Z, not statements inferred from a moment or numerical prefix.

A finite common root of U,V would, by(92), have to be1. Since U(1)=V(1)=24(lambda+mu), this happens precisely when lambda+mu=0. There is no additional projective common root at infinity: the a^2 coefficients of the two homogeneous forms are mu and15lambda+10mu, which cannot both vanish for a nonzero parameter pair. Hence the pair is polynomially and projectively coprime whenever lambda+mu!=0.

If also lambda+3mu!=0, equations(93)--(94) prove that its norm has EXACT contact4 at the fixed nonzero centre1. This includes mu=0: U then has dehomogenized degree1, but V has degree2 and the pair remains genuine of maximal degree2.

If(lambda,mu) is primitive, the joint integer coefficient content of U,V is1. Conversely that content is exactly gcd(lambda,mu): it divides the U coefficient mu and the V constant-lambda-2mu, hence divides both parameters; every coefficient is an integer combination of them.

### The two exceptional parameters

The two exceptional points in the parameter line are distinct. Up to scalar they are

    [lambda:mu]=[1:-1],  [3:-1].

At the first,

    U=(s-1)(5-s),  V=(s-1)(5s-1).

The common factor cancels, leaving the fixed linear diagonal Padé witness, whose reduced norm has contact3 at1. Its unreduced quadratic norm has contact5.

At the second,

    U=-s^2+14s+35,  V=35s^2+14s-1,

which is the fixed genuine quadratic diagonal Padé family of contact5. There are no further exceptional parameter values, by(94) and the coprimality argument.

### Every actual witness is represented over Q

Let x,y be an actual pair: distinct positive squarefree kernels in[T,2T), both units modulo N, with positive reduced unit witness u/v, u,v<=4T^2 and x^3u^2=y^3v^2 modulo N. Write x=g a,y=g b, g=gcd(x,y). Then a,b are primitive distinct positive squarefree units, and cancelling g^3 modulo N is valid.

Evaluate the basis forms at(a,b), and define the INTEGER parameters

    lambda=u V_2(a,b)-v U_2(a,b),
    mu=v U_1(a,b)-u V_1(a,b).                               (95)

By the homogeneous version of(92),

    Delta=U_1(a,b)V_2(a,b)-U_2(a,b)V_1(a,b)
          =-15(a-b)^4 !=0.

Therefore(95) gives the EXACT integer identities

    U(a,b)=Delta u,  V(a,b)=Delta v.                        (96)

In particular the parameter pair is nonzero. It is the unique projective parameter yielding this projective witness value at(a,b), since the evaluation matrix has nonzero determinant Delta.

Let k=gcd(lambda,mu) and divide to obtain a primitive parameter pair. Then k divides both Delta u and Delta v, so k|Delta because(u,v)=1. The evaluations of the normalized pencil pair are therefore

    (Delta/k)(u,v),

which reduce, after a common sign, to the original positive unit fraction. No inverse of Delta or k modulo N was used. They may be nonunits modulo N; exact equality OVER Q and integer divisibility are what certify the actual witness.

Thus EVERY actual witness, and hence every full-K kernel pair, belongs either to a genuine exact-contact4 member of this ONE fixed-centre pencil or to one of the two fixed exceptional rational witness families.

### Parameter heights and actual curves

Let Y=T/g. In the box a,b in[Y,2Y), elementary coefficient bounds give

    |U_1(a,b)|<=96Y^2, |V_1(a,b)|<=104Y^2,
    |U_2(a,b)|<=112Y^2,|V_2(a,b)|<=112Y^2.

Writing H=max(u,v), formula(95) yields

    max(|lambda|,|mu|)<=224 H Y^2
                         <=896 T^4/g^2.                    (97)

Dividing by k only decreases this height. This is an existence and height statement; it does NOT bound the number or frequency of the occupied primitive parameter values.

For a nonexceptional parameter the actual curve

    G(s,t)=V(s)-tU(s)=0,  s=a/b, t=v/u,

has fixed bidegree at most(2,1) and contains its assigned point OVER Q. With primitive parameters its joint remainder coefficient content is1, and its nonzero primitive cusp norm has ordinary degree at most7. Thus no hidden modulus loss is introduced by this t-linear cover. This does not mean its old coarse capacity is small, nor that the cover uses few members. Assigning one member per actual witness is always possible but gives no economy by itself.

### Same-exponent equivalence for savings

Let K_N(T) be the original full actual pair count. Let P_N(T) count those kernel pairs admitting at least one original actual witness from a NONexceptional pencil parameter, with all the original filters retained. Let E_N(T) count those admitting a witness from either fixed exceptional rational family.

**Theorem 25 (fixed-pencil/full-count equivalence).** For every fixed eta in(0,1], a uniform bound P_N(T)<<_epsilon T^(2-eta+epsilon) is equivalent to the corresponding uniform bound for the full K_N(T). Neither bound is asserted here.

As sets,

    K_set = P_set union E_set.

The two sets may overlap if a pair has different witnesses; no disjointness is asserted. Consequently

    P_N(T)<=K_N(T)<=P_N(T)+E_N(T).                          (98)

Corollary14, applied to the reduced fixed exceptional families, gives

    E_N(T)<<_epsilon T^(1+epsilon),                         (99)

since the fixed reduced linear family contributes at most T^(2/3+epsilon), and the fixed quadratic family contributes at most T^(1+epsilon). This uses the actual squarefree/unit/full-height conditions, not a raw congruence count. In the linear exceptional parameter the common factor s-1 is nonzero at the actual a/b, because a!=b, so its removal is legitimate over Q.

There is also an independent fixed-family proof of(99). For each reduced exceptional pair, the original actual congruence implies its necessary norm congruence by multiplying the pointwise evaluation gcd, as in Section1. Theorem2, applied to these TWO FIXED families with degrees1 and2 and M=3 and5, gives respectively O_(epsilon)(N^epsilon[T+T^2/N^(1/3)]) and O_(epsilon)(N^epsilon[T+T^2/N^(1/5)]). Under N comparable toT^5, both are O(T^(1+epsilon)) after reallocating epsilon. Their fixed coefficient-dependent constants are absolute and independent of N,T. This proves the same total exception bound without the new coefficient-uniform root theorem or the actual-height multiplier cap. The common factor s-1 in the first exceptional presentation was already removed over Q at the actual a!=b.

For any FIXED eta with0<eta<=1, equations(98)--(99) prove the exact same-exponent equivalence

    K_N(T)<<_epsilon T^(2-eta+epsilon) uniformly in balanced N
       IFF
    P_N(T)<<_epsilon T^(2-eta+epsilon) uniformly in balanced N. (100)

The forward direction is inclusion. The reverse adds(99), which is no larger than the desired exponent because2-eta>=1. The soft-epsilon quantifier is preserved by choosing the estimates with a smaller epsilon if needed.

For eta>1 this equivalence is not claimed: the fixed quadratic exceptional family already attains a constant times T pairs along infinitely many actual moduli. The natural full-K saving range is0<eta<=1.

### Scope of the reverse

The complete hard problem can therefore be concentrated in a one-parameter fixed-centre genuine quadratic pencil of contact4, with two quantitatively negligible fixed maximal-contact exceptions at the natural saving scales. This is a precise reverse bridge, not an estimate for that pencil. The parameter need not be bounded independently of the actual point; equation(97) is much too weak to supply an economical cover.

In particular, a codimension argument which removes the exact contact5 locus cannot by itself prove a full saving: the submaximal contact4 locus above can already interpolate EVERY actual witness. The stronger all-rational-centre maximal-contact union estimates remain valid, but they do not control this pencil automatically. A positive-saving parameter-frequency, primitive-ratio or total-capacity bound is still required. No full-K lower counterexample, nearquadratic example or new prescribed-sum exponent is asserted.
