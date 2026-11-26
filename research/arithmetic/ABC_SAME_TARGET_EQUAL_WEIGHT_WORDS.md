# Exact completion for equal positive weights on a common target

A mixed word with two equal positive weights and one negative
weight has subpower exact fibres on distinct rows sharing a
complete target. This holds for every fixed choice of the
weights and without an input radical cap. In the balanced
unequal-weight case, the remaining exceptional output has a
specific common divisor of the opposite summands. Its support
rather than its numerical height controls the cost of completion.

## Equal positive weights and arbitrary negative magnitude

For every fixed pair of positive integers \(h,k\) and
\(\varepsilon>0\), consider three pairwise distinct positive
primitive rows \(A_i+B_i=C\) of a common target \(C\le X\).
Put \(x=A_1/C,\ y=A_2/C,\ z=A_3/C\). Prescribe positive
rational numbers \(U,V\) with
\[
 x^h y^h z^{-k}=U,\qquad
 (1-x)^h(1-y)^h(1-z)^{-k}=V.                         \tag{CW1}
\]
The number of ordered triples is \(O_{h,k,\varepsilon}(X^\varepsilon)\)
uniformly in \(U,V\). Reversing the three signs gives the
same conclusion.

First divide the weights by \(g=\gcd(h,k)\) and take the
positive \(g\)-th roots of \(U,V\). If these roots are not
rational there is no completion. Thus assume that
\(\gcd(h,k)=1\). The target can be recovered in subpower
many ways by the common-target recovery in (CT3). Indeed
\(C^{|2h-k|}\) divides a fixed numerator or denominator
gcd when \(2h\ne k\). When \(2h=k\) the equality \(U=V\)
forces the three rows equal, and otherwise \(C\) divides
the reduced numerator of \(U-V\). A nonempty fibre
bounds the heights of both outputs by \(X^{2h+k}\).
[Common-target recovery and mixed-sign fibres (CT3)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#target-recovery-for-arbitrary-mixed-three-slot-weights)

If \(k=1\) and \(h=1\), (CT1) applies. If \(k=1\) and
\(h\ge2\), fix \(C\) and put \(P=A_1 A_2,\ Q=B_1 B_2\).
The exact relation is
\[
                     P^h/U+Q^h/V=C^{2h}.
\]
This is a binary homogeneous power equation at a nonzero
level with polynomial-height integer data after clearing
denominators. The uniform fixed-level bound counts at
most \(X^\varepsilon\) pairs up to an arbitrary adjustment
of \(\varepsilon\). For degree two use the uniform quadratic
bound. For higher fixed degree first divide by the common
divisor of \(P,Q\), whose \(h\)-th power divides the level,
then apply the coefficient-uniform Thue–Mahler bound to the
primitive pair. Each exact \(P\) has divisor-many ordered
factorizations, and the formula for \(z\) recovers the third
row. The same argument works for unequal positive weights
\(a,b\) with negative magnitude \(k\) if \(k\) divides
both and \(\gcd(a,b)/k\ge2\). After taking the \(k\)-th roots
use the products \(P=A_1^{a/\gcd(a,b)}A_2^{b/\gcd(a,b)}\)
and \(Q=B_1^{a/\gcd(a,b)}B_2^{b/\gcd(a,b)}\), at power
\(h=\gcd(a,b)/k\). Knowing \(P\) still allows divisor-many
choices of \(A_1,A_2\).
[Split factorization for weights \((1,1,-1)\) (CT1)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#subpower-completion-for-the-basic-mixed-sign-word)
· [Uniform binary fixed-level recovery](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#proof-of-recovery-and-enumeration)

Suppose \(k\ge2\). Now \(k\nmid h\), and the cap-free case
of (CT3) applies when \(U,V\ne1\). In fact each free letter
of one positive row is a \(k\)-th power times a factor on
a fixed prime support. If \(2h\ne k\), the negative or
positive valuation at a prime of \(C\) makes both \(U,V\)
different from one. The remaining primitive weight pattern
is therefore \((1,1,-2)\) with either \(U=1\) or \(V=1\).
These two possibilities are symmetric.

Take \(U=1\) and \(V=p/q\ne1\) in lowest terms. Then
\[
  xy=z^2,\qquad
  (x-y)^2=(1-z)^2(1-V)\big((1+z)^2-V(1-z)^2\big).
\]
The instance \(U=V=1\) has already been excluded. Let
\(D=A_1-A_2\) and \(B_3=C-A_3\). Clearing denominators gives
\[
 q^2 C^2 D^2
     =(q-p)B_3^2\big(q(C+A_3)^2-p(C-A_3)^2\big).
\]
Since \(B_3\) is coprime to \(C\), this implies
\(J=qD/B_3\in\mathbb Z\). Completing the square yields
\[
 \big((q-p)A_3+C(q+p)\big)^2-(CJ)^2=4pq C^2 .
\]
The right-hand side is a fixed nonzero integer of polynomial
height. Its signed divisor pairs determine \(A_3,J\) because
\(q-p\ne0\). They give \(D\). The exact sum of the first two
summands follows from \(x+y=1+z^2-V(1-z)^2\), so all three
rows follow at divisor cost. This completes (CW1), including
the exceptional zero-total-weight output.

For a balanced chart, the fixed exact record may be replaced
by the two-coordinate modular word of (FP1), with a loss at
most \(O(n^\eta X^\varepsilon)\) for any fixed \(\eta>0\) at
sufficiently small fixed budgets. Use (FP1) to list possible
exact pairs and then (CW1) for each. In a longer fixed word
one may fix the other slots and complete a marked, pairwise
distinct common-target triple with the indicated weights.
These statements also apply to primitive quotients of exact
least repairs, since the common multiplier leaves the two
fraction coordinates unchanged.
[Finite signed probe recovery (FP1)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#finite-probes-by-möbius-inversion)

## Exceptional balanced unequal weights and a common divisor

Let \(a,b>0\) be fixed integers, \(d=a+b\), and prescribe an
exact word with weights \((a,b,-d)\) on pairwise distinct
primitive rows with a common target \(C\le X\). Suppose one
of \(U,V\) is one. After interchanging the summands, assume
\(U=1\) and \(V\ne1\). Recover \(C\) at divisor cost as
above. Let \(P_0\) consist of primes dividing \(C\,abd\)
and the reduced numerators and denominators of \(V\) and
\(V-1\). The prime product of \(P_0\) is of polynomial
height in \(X\) whenever the fibre is nonempty. Let \(G\)
be the part of \(\gcd(B_1,B_2)\) supported outside \(P_0\).

**CW2.** For a prescribed \(G\) the exact word has only
\(O_{a,b,\varepsilon}(X^\varepsilon)\) such completions.
Consequently the restriction \(\operatorname{rad}(G)\le Z\)
has \(O_{a,b,\varepsilon}(Z X^\varepsilon)\) completions,
uniformly for \(Z\ge1\). If \(V=1\) instead, use the gcd
of the other two summands. An existing \(\rho/4\) fibre
bound under the original radical caps is available
independently of the gcd cost.

To prove this, first note that no prime outside \(P_0\)
divides both \(A_1,A_2\). Such a prime divides \(A_3\)
by the equality \(A_1^a A_2^b=A_3^d\). At that prime
all \(B_i/C\) would be units congruent to one, forcing
\(V\equiv1\), contrary to the choice of \(P_0\). Therefore
every prime on \(A_1\) outside \(P_0\) has valuation
divisible by \(r=d/\gcd(a,d)\ge2\).

At a prime outside \(P_0\), the equality
\(B_1^a B_2^b=V B_3^d\) gives
\(a e_1+b e_2=d e_3\), where \(e_i=v_p(B_i)\).
If the prime occurs on \(B_1\) but not on \(B_2\), then
\(r\mid e_1\). If it occurs on both, all three exponents
are positive and equal. For if \(e_1<e_2\), the weighted
mean \(e_3\) lies strictly between them. The identity
\[
        (1-B_1/C)^a(1-B_2/C)^b=(1-B_3/C)^d
\]
has a nonzero linear term of valuation \(e_1\) on the
left modulo \(p^{e_1+1}\) and no such term on the right.
This is a contradiction because \(p\nmid aC\). The case
\(e_2<e_1\) uses \(p\nmid bC\) in the same way.

It follows that \(A_1=\alpha R^r\) and \(B_1=\beta G S^r\)
where \(\alpha,\beta\) are supported on \(P_0\). Fix their
values at subpower cost and use
\(\alpha R^r+\beta G S^r=C\). The same uniform fixed-level
binary bound recovers the first row. Once it is known,
the derivative argument in (SP1) permits at most two
completions of the other two distinct rows. In
the proportional case where that derivative vanishes
identically, the two rows coincide and are excluded.
The divisor enumeration of \(C\) and the enumeration of
the support and full powers in \(G\) cost
\(O(Z X^\varepsilon)\) in total, with the inner allowances
adjusted. This proves (CW2).
[Mixed-sign completion and derivative count (SP1)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#mixed-sign-completion-using-the-three-radical-caps)

In a common-target fibre of size at least \(n^u\), under
individual radical caps \(n^\rho\) on the free letters,
most ordered triples have either free-role pair overlap
at most \(n^\tau\) whenever \(\tau>\rho-2u\) by a fixed
margin. For a fixed preceding row the support of that
role is a prescribed pool. A later row with larger overlap
has outside cap at most \(n^{\rho-\tau}\), and (FL1) gives
\(O(n^{(\rho-\tau)/2+\eta}X^\varepsilon)\) choices.
Sum over positions and choose \(\eta\) and the inner
height loss below the margin. This applies when
\(0<\tau\le\rho\); larger thresholds require no deletion.
In particular near the largest balanced target multiplicity,
(CW2) gives arbitrarily small normalized exact and modular
word fibres on the typical triples, even for the
exceptional output \(U=1\) or \(V=1\). In this application
the exact-pair list of (FP1) is formed first, then (CW2)
is applied to each exceptional pair, while the other
pairs use the cap-free branch of (CT3).
[Fixed-letter counts with prescribed prime pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## A rational map with the same input and output denominator

**CW3.** For every fixed integer \(a\ge1\), a prescribed
\((a,1,-1)\) exact word has \(O_{a,\varepsilon}(X^\varepsilon)\)
pairwise distinct common-target representations, with no
radical hypothesis. When \(a=1\) use (CT1).

Let \(a\ge2\). The target is recovered by the same
denominator gcd formula as in (CT1), now giving \(C^a\).
Eliminating \(z\) gives
\[
 y=\frac{U(V-(1-x)^a)}{V x^a-U(1-x)^a}.
\]
If both polynomials vanish at \(x\), then \(U=x^a\) and
\(V=(1-x)^a\), so \(z=y\) and that point is discarded.
If just the denominator vanishes there is no solution.

Cancel their polynomial gcd and clear denominators to
obtain coprime integer polynomials \(P,Q\) of degree
at most \(a\) and polynomial coefficient height in \(X\).
This follows from Gauss's lemma and the coefficient
bound for factors of a fixed-degree integer polynomial.
In detail, the product has a positive minimum sup norm
on pairs of normalized coefficient vectors of bounded
degree. Each nonzero integral cofactor before
normalization has sup norm at least one.

Suppose \(Q\) is nonconstant. Homogenize \(P,Q\) to their
common maximum degree, obtaining forms \(P_h,Q_h\)
with nonzero homogeneous resultant \(R\). There is no
common root at infinity since at least one polynomial
has that maximum degree. The two homogeneous resultant
identities give \(R\) times a fixed power of each input
coordinate in the ideal \((P_h,Q_h)\). Hence
\(\gcd(P_h(A_1,C),Q_h(A_1,C))\mid R\), since
\(\gcd(A_1,C)=1\). The reduced denominator of \(y\)
is exactly \(C\). Therefore \(Q_h(A_1,C)\) is \(C\)
times a signed divisor of \(R\). For each divisor
the resulting polynomial in \(A_1\) is nonconstant
and has at most \(a\) roots. The resultant has
polynomial height, so the count is subpower, and
the displayed map and the first word equation
recover the remaining rows.

In the polynomial case, the original denominator has
degree \(a\) except when \(a\) is even and \(U=V\), when
it has degree \(a-1\). In the first case a polynomial
quotient would be constant, with value in \((0,1)\) in
a nonempty fibre. But
\[
        Vx^a y+U(1-x)^a(1-y)=UV
\]
then has a nonzero linear coefficient at \(x=0\).
In the second case put \(U=V=t\). A polynomial
quotient must be affine with endpoint values \(1-t,t\).
Comparing the coefficient of \(x\) in
\[
              x^a y+(1-x)^a(1-y)=t
\]
gives \(t=1/(a+2)\). For even \(a\ge4\), its coefficient
of \(x^2\) is then \(-ta(a+1)/2\), again a contradiction.
The remaining possibility \(a=2,t=1/4\) forces \(C^2\)
to divide \(4\) by target recovery, so it contributes
only boundedly many rows. This completes (CW3).
[Common-target denominator recovery (CT1)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#subpower-completion-for-the-basic-mixed-sign-word)

## A flex for the smallest unequal balanced word

**CW4.** The exact word with weights \((1,2,-3)\) has
\(O_\varepsilon(X^\varepsilon)\) pairwise distinct
common-target representations, without radical caps.
The exceptional branch \(U=1,V\ne1\) in fact satisfies
this bound even when the three denominators vary
independently through \(X\).

On a common target, \(U=V=1\) cannot occur for distinct
rows by (CT3), and the generic \(U,V\ne1\) case is its
binary fixed-level branch. For the remaining branch
interchange summands if necessary and take \(U=1,V\ne1\).
Set
\[
       S=y/z,\qquad T=(1-y)/(1-z).
\]
Then \(x/z=S^{-2}\) and \((1-x)/(1-z)=V/T^2\).
Eliminate \(z\) using the two identities
\(Sz+T(1-z)=1\) and \(S^{-2}z+(V/T^2)(1-z)=1\).
Division by \(S-1\), which is nonzero because \(y\ne z\),
gives
\[
 (V-T^2)S^2+(T-1)T^2 S+(T-1)T^2=0.
\]
Here \(T\ne1\) by distinctness. If \(V=T^2\), the equation
requires \(S=-1\), which is impossible. Define
\(\Delta=2(V-T^2)S+(T-1)T^2\). Its square equals
\(T^2(T-1)(T^3+3T^2-4V)\). With
\(s=(T-1)^{-1}\) and \(Z=\Delta s^2/T\) this is
\[
                     Z^2=(3s+1)^2+4(1-V)s^3 .
\]
The cubic is smooth. The discriminant of the polynomial
on the right is \(432 V(1-V)\ne0\). Its point \((0,1)\)
is a flex with tangent \(Z=1+3s\), as substitution
makes the intersection cubic a nonzero multiple of
\(s^3\). With the point at infinity as origin the flex
is a nonzero rational point of exact order three.

Every coordinate formed from \(x,y,z\) and the word
value has height \(X^{O(1)}\). Put \(A=4(1-V)\) and
replace \((s,Z)\) by \((As,AZ)\) to make the cubic
monic. A rational Weierstrass scaling clears coefficient
denominators, with polynomial height for all coefficients
and point coordinates. Apply Dujella's uniform height
count over \(\mathbb Q\) for curves with rational
\(3\)-torsion, taking the height threshold to be a
sufficiently large fixed power of \(X\). It gives
\(O_\varepsilon(X^\varepsilon)\) rational points.
Recover \(T,\Delta,S\) from \(s,Z\), then
\(z=(1-T)/(S-T)\), \(y=Sz\) and \(x=z/S^2\).
The denominator \(S-T\) cannot vanish since it would
force \(S=T=1\). This gives one primitive row per
fraction and proves the assertion.
[Dujella, Corollary 1.3](https://arxiv.org/html/2312.03655v1#S1)
· [Generic unequal-weight completion (CT3)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#target-recovery-for-arbitrary-mixed-three-slot-weights)

These completion estimates control the multiplicity of
occupied signed words. An upper bound on how many such
words occur with both original summand radical caps is
a separate requirement for a global target-tail saving.
