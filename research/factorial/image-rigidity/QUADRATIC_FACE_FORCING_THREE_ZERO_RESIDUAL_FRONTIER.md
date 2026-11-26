# Quadratic face-forcing three-zero residual frontier

The proof uses two precise local exclusions and gives an exact arithmetic
criterion equivalent to the remaining universal three-window. It does not
assert that this criterion holds or prove a full binary quadratic theorem.

Let K have characteristic zero, define
h_j(r,z)=Lambda_1((r x^2+x+z)^j), Lambda_1(x^d)=d!, and
T_j=h_j+j h_(j-1), j>=1. Assume z!=0.

**Theorem 1 (necessary residual family).** If
T_k=T_(k+1)=T_(k+2)=0 for k>=1, and m=k+1, then for an integer q>=1
with q congruent to1 or2 modulo4 one necessarily has

    m=q(4q+3)^2,
    z=-8q(q+1)(2q+1) +/- sqrt(-3q(q+1)).                 (1)

Every embedding of either displayed quadratic value into Q_2bar has
v_2(z)=1/2. In particular k>=48. This is a necessary condition only;
no actual common root r in this residual family is asserted.

## 1. Exact inputs and algebraicity

The native factorial recurrence, or its pivot-free elimination, gives
for every three-zero block

    Pi_m(z)=2(z+m)^3+3(z+m)^2+m=0.                       (2)

The only inputs for the necessary residual classification are these two
proved local theorems:

* The [negative-branch exclusion](QUADRATIC_FACE_FORCING_NEGATIVE_TWO_ADIC_BRANCH.md)
  rules out every triple with v_2(z)=-1. Its only imported ingredient is the
  all-unit divided-power factorization.
* The [positive-branch exclusion](QUADRATIC_FACE_FORCING_POSITIVE_VALUATION_BRANCHES.md)
  rules out every triple with v_2(z)>=1. It includes the two exact native
  resultant/Bezout certificates required after its all-index argument.

For clarity the cubic constraint is itself a local polynomial identity.
The factorial derivative identity gives

 h_j=2rj(2j-1)h_(j-1)+j(j-1)(1-4rz)h_(j-2)+z^(j-1)(z+j).

Put D_i=2i-1+2z. Direct elimination gives the universal syzygy

 z^(i-1)Pi_i(z)=D_i T_(i+1)
   +2(i+1)[1-r(2i+1)D_i]T_i
   -i(i+1)(D_i+2)(1-4rz)T_(i-1), i>=2.                 (2a)

One obtains it by writing
T_i-i(1-4rz)T_(i-1)=2riD_i h_(i-1)+z^(i-1)(z+i),
combining consecutive equations, and using
h_i=T_i-i h_(i-1). No coefficient-dependent pivot is divided.
Thus a triple with middle i gives(2), since z!=0.

Equation(2) makes z algebraic over Q. As a polynomial in r, T_k has degree k
and leading coefficient (2k)!, independent of z. Therefore r is algebraic
over Q(z) as well. Every alleged solution lies in a number field inside K,
and all its number-field embeddings can be extended into Q_2bar. There is
no assumption that the original K embeds as a whole into C or Q_2bar.

Putting u=z+m and w=2u, equation(2) is

    f_m(u)=2u^3+3u^2+m=0,
    w^3+3w^2+4m=0.                                     (3)

Modulo2 the second polynomial is w^2(w+1). It has exactly one unit root
in Q_2bar, which is the unique Hensel lift of the simple root1 and belongs
to Q_2. Its corresponding z has valuation -1. The other two z roots have
valuation v_2(m(m+1))/2: the z polynomial has lower Newton vertices
(0,v_2(m(m+1))),(2,0),(3,1).

## 2. Galois amplification of the negative-branch theorem

If Pi_m is irreducible over Q, every one of its algebraic roots is
Q-conjugate to its negative-valuation root in Q_2bar. An embedding taking
z to that root extends from Q(z) to Q(z,r), preserving all three T equations.
The negative-branch theorem gives a contradiction. Thus Pi_m must be
reducible over Q.

A reducible cubic over Q has a rational root. Its discriminant is
-108m(m+1)<0 for m>=2, so it has exactly one real root. Consequently the
rational root is unique, and its remaining quadratic factor is irreducible
over Q. By the rational-root theorem that rational u is either an integer
or half an odd integer. Since f_m is positive at -1 and0, with both stationary
values positive, its real root lies below -1.

For an integer root u=-h, h>=2, equation(3) is equivalent to

    m=h^2(2h-3),
    z=-h(h-1)(2h-1).                                   (4)

The displayed rational z has valuation at least1, so the positive-branch
theorem excludes it. The negative-valuation root is in the irreducible
quadratic factor, since the rational root is integral. The two quadratic
roots are conjugate over Q, so the same embedding argument excludes both.
There is no solution in this case.

For a nonintegral rational root put u=-b/2 with b odd and b>=3. Then
m=b^2(b-3)/4. Integrality of m forces b congruent to3 modulo4. Since m>=2,
write b=4q+3 with q>=1, giving

    m=q(4q+3)^2,
    u_0=-(4q+3)/2.                                     (5)

This rational root is precisely the negative-valuation root and is excluded.
The remaining two u values solve

    u^2-2q u+q(4q+3)=0.

They are u=q +/- sqrt(-3q(q+1)); subtracting m gives exactly(1).

## 3. The exact residual valuation

Let t=v_2(q(q+1))>=1. In(1) the rational center has valuation3+t,
whereas the square-root term has valuation t/2. These are unequal, so the
valuation of either sum is exactly t/2; no cancellation assumption is made.
If t>=2, the positive-branch theorem excludes the quadratic factor. Thus
t=1, equivalent to q congruent to1 or2 modulo4. This proves the theorem.

The first two residual middle indices are m49 and242, corresponding to
forcing starts k48 and241. Examples of later eligible q are5 and6. This
is an algebraic classification of necessary possibilities, not a search
or a claim that any displayed family is realizable.

## 4. Residual starting indices for symmetric five-zero blocks

The only additional input in this section is the polynomial identity in
[the full symmetric quadratic proof](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md),
Section2, Equation10. It applies to
f=a(x^2+y^2)+bxy+e(x+y)+c even when b=0:

**Corollary 2 (necessary starting indices).** If ec!=0 and five consecutive
factorial moments of this displayed symmetric phase vanish starting at
n>=1, then n has the form(6) below; in particular n>=46.

The imported identity is

    M_(j+1)=(2a+b)(j+1)^2 M_j
       +j(j+1)(e^2-(2a+b)c)M_(j-1)+H_(j+1)+(j+1)e H_j,

where H is the quadratic face sequence. If a hypothetical FIVE-zero M
block starts at n, with ec!=0, it gives the three normalized T zeros
at n+2,n+3,n+4. Thus its middle m is n+3, and this note forces

    n=q(4q+3)^2-3, q>=1, q=1 or2 mod4.                  (6)

In particular no such five-zero block starts at n<=45. The first possible
start under this necessary reduction is46. This conclusion applies to the
displayed symmetric family and remains only a necessary residual condition;
it does not prove the actual support-five b0 assertion at the starts in(6),
and says nothing about general nonsymmetric binary quadratics or SFC(2)/SFC(3).

## 5. An exact remaining coprimeness criterion

Define the following two universal assertions, neither of which is proved
by the preceding necessary reduction.

**Assertion U.** For every characteristic-zero field K, every r,z in K
with z!=0, and every k>=1, the three values T_k,T_(k+1),T_(k+2) are not
all zero.

**Assertion R.** For every q>=1 congruent to1 or2 modulo4, put
m=q(4q+3)^2 and let L_q=Q(omega_q), where omega_q^2=-3q(q+1).
This is a quadratic number field. Put
z_q=-8q(q+1)(2q+1)+omega_q. Then

 gcd_(L_q[r])(T_(m-1)(r,z_q),T_m(r,z_q))=1.              (7)

The minus sign in front of omega_q gives the conjugate criterion and
therefore the same assertion. Both polynomials have nonzero constant
r-leading coefficients, so coprimeness is equivalent to having no common
finite root in an algebraic closure of L_q.

**Proposition 3 (reversible remaining criterion).** Assertions U and R
are equivalent. This is a proved equivalence context, not a proof of
either universal assertion.

For U implies R, suppose the pair in(7) has a common root. The exact
quadratic expression gives Pi_m(z_q)=0. Also D_m=2m-1+2z_q is nonzero:
z_q is a nonrational quadratic element, whereas D_m=0 would make it
rational. Identity(2a), with the pair zero, now gives
D_m T_(m+1)=0, hence T_(m+1)=0. We obtain a three-zero block starting
m-1>=1 in a characteristic-zero field, contradicting U. Equivalently,
at the 2-adic embeddings D_m is a unit because v_2(z_q)=1/2.

For R implies U, any hypothetical triple is classified by Theorem1.
Its m and z identify a Q-embedding of L_q into its coefficient field,
with omega_q sent to one of the two square roots. The common adjacent
pair contradicts the Bezout identity of(7) after that embedding; the
conjugate criterion handles the other sign. Thus the triple is impossible.
No new universal member is seeded merely by proving these implications.

## 6. Conditional connection to the symmetric support-five target

**Corollary 4 (one-way conditional five-window).** If either equivalent
Assertion U or Assertion R holds, then for every characteristic-zero
field, every a,e,c with aec!=0, every b (including b=0), and every n>=1,
f=a(x^2+y^2)+bxy+e(x+y)+c has a nonzero factorial moment in
n,...,n+4. In particular the b=0 family has its actual support-five
Strong Factorial window under this remaining hypothesis.

Indeed five zero interior moments would, by the same imported symmetric
recurrence of Section4 and phase scaling by e, give three zero normalized
T values at n+2,n+3,n+4. Assertion U excludes them. Proposition3 allows
Assertion R as an equivalent sufficient hypothesis. The implication is
one-way: no converse from this symmetric support-five target to U or R
is asserted, and the target is not proved by this corollary alone.

The global necessary family, the starting-index exclusion, and the proved
implications above leave U and R open. General nonsymmetric binary
quadratics and the full dimension-two or dimension-three Strong Factorial
conjectures are outside these conclusions. All proof inputs and nonzero
pivots have been stated explicitly; finite diagnostics are not premises.
