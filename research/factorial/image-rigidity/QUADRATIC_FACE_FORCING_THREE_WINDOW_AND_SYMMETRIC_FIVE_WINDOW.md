# Direct backward-valuation exclusion of the half-valuation forcing branch

Let Lambda(X^d)=d!, h_n(r,z)=Lambda((rX^2+X+z)^n), and
T_n(r,z)=h_n+n h_(n-1), n>=1. Normalize v_2(2)=1.

**Theorem 1 (unconditional local half exclusion).** There are no r,z in
Q_2bar, with v_2(z)=1/2, and no integer k>=1 for which

                    T_k=T_(k+1)=T_(k+2)=0.                       (1)

The proof uses only native factorial coefficients, a polynomial three-term
identity, and exact nonarchimedean valuation comparisons. It does not use
an anti-coefficient hypothesis, a formal completion, a Newton-factor transfer,
or divided-power Hensel factorization.

## 1. Native coefficient bounds with only v(z)=1/2

Put B=1/z and define the ordinary normalized native reciprocal

 Q_n(Z)=(4Z/z)^n T_n(z/(4Z),z)/(2n)!, n>=1,  Q_0=1.             (2)

It has constant1. Literal multinomial expansion, combining the two terms
of T_n, gives

 [Z^j]Q_n=(2^j/j!) b_j^(n),
 b_j^(n)=[1+sum_(l=1)^j 2(i+l)binom(j,l)(2i+1)_(l-1)B^l]
                     /D_(n,j),
 i=n-j, D_(n,j)=prod_(a=0)^(j-1)(2n-2a-1), 0<=j<=n.             (3)

The rising products are literal integer products; D_(n,j) is an odd rational
unit. To see the numerator, the first h_n coefficient has relative term
binom(j,l)(2i+1)_l B^l. The n h_(n-1) contribution is
l binom(j,l)(2i+1)_(l-1)B^l. Their sum is exactly the displayed correction.

**Lemma 2 (the only native coefficient facts needed).** For every n and j
in(3), b_j^(n) is integral. It is a unit if 4 divides j, or if j=n.
Consequently

 v([Z^j]Q_n)>=s_2(j),
 v([Z^j]Q_n)=s_2(j) if 4|j or j=n.                            (4)

Here s_2 denotes binary digit sum. These conclusions need no congruence
condition on z^2/2 and no field-conjugation hypothesis.

Indeed, an interval of l-1 integers has factorial valuation at least
v_2((l-1)!), so a correction term in(3) has valuation at least
l/2-s_2(l-1). For l>=5 this is at least1. The l=1 and l=3 bounds are
at least1/2. At l=4 its valuation is at least
v_2(i+4)+v_2(i+1)>=1. Thus only l=2 might have valuation0. That term is
2(i+2)(2i+1)binom(j,2)B^2. If 4|j, the binomial coefficient is even;
if j=n then i=0 and i+2 is even. In either case this term too has positive
valuation. The numerator is then1 plus terms of positive valuation and is
a unit. Otherwise all corrections are still integral. Finally
v_2(2^j/j!)=s_2(j), proving(4). The elementary inequality for l>=5 follows
from s_2(l-1)<=(l-2)/2 for l>=5, with the integer/half-integer interpretation;
small l=5,6,7,8 can also be checked directly and the doubling induction then
proves the rest.

For p=1 modulo4, every proper highest-bit prefix of p is divisible by4,
and its final prefix p has unit coefficient by Lemma2. The lower Newton
polygon of Q_p is therefore the standard binary polygon: a face of width d
and slope1/d for each binary bit d of p. In detail, its vertices are
(P_j,j), where P_j is the sum of the j largest bits of p. Between consecutive
vertices, an interior index is P_j+u with 0<u<d; its digit sum is
j+s_2(u)>=j+1>j+u/d. This proves both the lower bound and the exact vertices.
For p+1=2 modulo4, the same argument applies: all proper highest-bit prefixes
are divisible by4, and the final prefix has unit coefficient. Its binary
bits equal those of p except that the final bit1 is replaced by bit2.

## 2. Selecting the common higher radius and excluding r=0

A hypothetical triple has middle index M=k+1 and satisfies the exact cubic
constraint

             Pi_M(z)=2(z+M)^3+3(z+M)^2+M=0.                     (5)

A derivation is supplied in Section4. Its expanded z^2 coefficient is odd,
the z^3 and z coefficients at v(z)=1/2 have valuation at least5/2, and its
constant has valuation v_2(M(M+1)). Thus(5) forces v_2(M(M+1))=1:
M is1 or2 modulo4. Select

               p=M if M is odd, p=M-1 if M is even.

Then p=1 modulo4 and both T_p,T_(p+1) are zero. Lemma2 at the top coefficient
shows T_p(0,z)!=0, so r!=0. Hence

                         xi=z/(4r)

is defined, and Q_n(xi)=0 whenever T_n(r,z)=0. Section1 shows that every
root of Q_p has valuation -1/d for a binary bit d of p. The simultaneous
zero of Q_(p+1) excludes d=1, because p+1 has bit2 in place of bit1 and no
bit1. Thus there is a common higher bit e=2^rho>=4 of p and p+1 such that

                         v(xi)=-1/e.                           (6)

In particular p>=5; the possible p=1 case has already been excluded by the
disjoint degree-one and degree-two radii.
Write

 p=D+e+S, 0<S<e, S=1 modulo4,
 H=s_2(D), ell=s_2(S), nu=H-D/e.                                (7)

Here D is the sum of the bits of p strictly above e, so D/e is even.
No assumption is made that e is the highest bit of p.

## 3. The precise native Gauss bound at this radius

**Lemma 3.** At the xi in(6),

 v(Q_n(xi))>=nu if D+e<=n<=D+e+S,
 v(Q_(D+e-1)(xi))=nu.                                         (8)

No assertion of a full standard binary polygon for arbitrary n is needed.
To prove it, put A=D/e (an even integer) and write any coefficient index
j=ea+t, 0<=t<e. By(4), its term has valuation at least

 f(a)+s_2(t)-t/e,  f(a)=s_2(a)-a.

For t>0 the last two terms are strictly positive. Also
f(a+1)-f(a)=-v_2(a+1)<=0. For n<=D+e+S<D+2e, one has a<=A+1, and
f(A)=f(A+1)=nu because A is even. This proves the lower bound.
For n=D+e-1, a<=A. If A>0, f(A-1)=f(A)+v_2(A)>f(A), so every a<A has
strictly larger value; if A=0 there is no such a. Thus the unique minimizing
index is j=D. It is divisible by4, has exact coefficient valuation s_2(D)
by(4), and contributes precisely nu. The minimum cannot cancel, proving
the equality in(8).

## 4. Exact normalized recurrence and all its pivots

The elementary factorial derivative identity gives

 h_n=2rn(2n-1)h_(n-1)
       +n(n-1)(1-4rz)h_(n-2)+z^(n-1)(z+n).                    (9)

For a derivation, apply Lambda(R')=Lambda(R)-R(0) to
phi' phi^(n-1) and phi^n, with phi=rX^2+X+z and
(phi')^2=4r phi+1-4rz. Combining consecutive equations gives the universal
polynomial syzygy, for i>=2 and d_i=2i-1+2z,

 z^(i-1)Pi_i(z)=d_i T_(i+1)
 +2(i+1)[1-r(2i+1)d_i]T_i
 -i(i+1)(d_i+2)(1-4rz)T_(i-1).                               (10)

One can verify the elimination without any pivot division by first writing
T_i-i(1-4rz)T_(i-1)=2ri d_i h_(i-1)+z^(i-1)(z+i),
and then combining the two adjacent indices using h_i=T_i-i h_(i-1).
At the hypothesized triple(10) proves(5), since z!=0.

Set delta_n=2n-1+2z and epsilon_n=4^n/(2n)!. Solving(10) with i=n-1,
then applying normalization(2), gives for n>=3 the EXACT recurrence

 Q_n(Z)=A_n(Z)Q_(n-1)(Z)+C_n(Z)Q_(n-2)(Z)+R_n(Z),
 A_n=1-4Z/[z(2n-1)delta_(n-1)],
 C_n=4(Z^2/z^2-Z)delta_n/[(2n-1)(2n-3)delta_(n-1)],
 R_n=epsilon_n Z^n Pi_(n-1)(z)/[z^2 delta_(n-1)].              (11)

All delta_n are units, being odd integers plus2z of positive valuation.
At xi, the nonconstant term of A_n has valuation3/2-1/e>0, so

 v(A_n(xi))=0,   v(C_n(xi))=b:=1-2/e>0.                       (12)

The latter uses the unique least term xi^2/z^2 in xi^2/z^2-xi,
whose valuation is -1-2/e. Thus no C_n(xi) is zero. The other denominators
are z!=0 and odd rational integers. Finally

 v(epsilon_n)=s_2(n),
 v(R_n(xi))=s_2(n)-n/e+v(Pi_(n-1)(z))-1,                     (13)

with the value infinity if the numerator is zero. The cubic has integer
coefficients and z is integral, so even the weak bound
v(R_n(xi))>=s_2(n)-n/e-1 will suffice in the induction. In fact each Pi_j
has valuation at least1, but that extra unit is not needed there.

## 5. The two initial values in BOTH middle parities

The exact difference is

                   Pi_(i+1)-Pi_i=6(z+i+1)^2.                    (14)

If M=p is odd, Pi_p=0 implies
v(Pi_(p-1))=v(Pi_(p-2))=1: the first uses z+p a unit, and the second subtracts
6(z+p-1)^2, which has valuation2 because v(z+p-1)=1/2.
The initial zeros are Q_p(xi)=Q_(p-1)(xi)=0. At n=p, equations(11)--(13)
give

 v(Q_(p-2)(xi))=nu+ell-1-(S-2)/e.                              (15)

At n=p-1, the residual has value nu+ell-1-(S-1)/e, which is strictly
smaller than the A_(p-1)Q_(p-2) term by1/e. Therefore

 v(Q_(p-3)(xi))=nu+ell-2-(S-3)/e.                              (16)

If M=p+1 is even, the initial zeros instead are
Q_p(xi)=Q_(p+1)(xi)=0. Equation(14) gives

 v(Pi_p)=2,  v(Pi_(p-1))=v(Pi_(p-2))=1.

The first uses z+p+1 of valuation1/2; the next two differences each have
exactly one unit squared term. Applying(11) at n=p+1 first gives

 v(Q_(p-1)(xi))=nu+ell-(S-1)/e.                               (17)

At n=p the residual has value nu+ell-S/e, strictly smaller than the
A_p Q_(p-1) term by1/e, proving exactly(15). At n=p-1 the residual has
value nu+ell-1-(S-1)/e. This is below Q_(p-1) by1 and below
A_(p-1)Q_(p-2) by1/e, proving exactly(16). Thus(15)--(16) hold in both
parities without imposing a zero at the wrong preceding index.
All comparisons are strict; no unanalysed equal-valuation cancellation is used.

If S=1, then ell=1 and(15) is nu+1/e, while p-2=D+e-1. This contradicts
the exact equality in(8). Hence from now on S>=5 and ell>=2.
In this case all recurrence indices just used are >=3.

## 6. Exact backward induction and its stopping index

For d>=0 define the nonnegative digit excess

 E_d=s_2(p-d)-s_2(p)+d.

For 0<=d<p,
E_(d+1)-E_d=v_2(p-d)>=0. Since p=1 modulo4 and p>=5,
E_2=v_2(p-1)>=2. Thus E_d>=2 for every d>=2 in the range below.

We claim that the initial value(16) propagates as

 v(Q_(p-j)(xi))=nu+ell-j+1-(S-2j+3)/e,
                         3<=j<=ell+1.                         (18)

The case j=3 is(16). To prove the next values, put n=p-d with
2<=d<=ell-1. Assume the two previously known values. At d=2,
(15) is strictly above(16) by1-1/e. At every later d, consecutive
values in(18) differ by b=1-2/e>0. Hence A_n Q_(n-1) is below Q_n.

Compare the residual with Q_(n-1). Even using only the WEAK bound following
(13), its valuation exceeds the predicted Q_(n-1) value by at least

                         E_d-1-(d-1)/e>0.                     (19)

For clarity, s_2(p)=H+1+ell, n=p-d, and
s_2(n)=s_2(p)-d+E_d; substitution gives(19) directly. It is strictly
positive because E_d>=2 and d-1<e (indeed d<=ell-1 and ell<=S<e).
Thus A_n Q_(n-1) is the unique least numerator term in

 Q_(n-2)=[Q_n-A_n Q_(n-1)-R_n]/C_n.

Division subtracts exactly b from its valuation, proving(18) at j=d+2.
If ell=2, no induction step is required: the stopping value is already(16).
This completes(18) with every possible zero residual included.

Because S=1 modulo4 with ell set bits, S>=1+4(ell-1)=4ell-3.
For ell>=2 this gives S>2ell-1 and S>=ell+1. At j=ell+1, equation(18)
therefore becomes

 v(Q_(p-ell-1)(xi))=nu-(S-2ell+1)/e<nu.                        (20)

But p-ell-1>=D+e and <=D+e+S. This is in the first range of(8), which
requires its valuation to be at least nu. The contradiction proves
Theorem1. The proof is finite in each hypothetical index but uniform in
all indices; no bounded experiment is extrapolated.

## 7. Pair coprimeness, global three-windows, and the symmetric application

**Corollary 4 (critical half pair).** If M>=2, Pi_M(z)=0 and v(z)=1/2,
then T_(M-1)(r,z) and T_M(r,z) are coprime as polynomials in r over Q_2(z).
Indeed any common algebraic root would, by(10), also satisfy T_(M+1)=0:
d_M is a unit. That would contradict Theorem1. This is a statement about
all roots in the algebraic closure, not just a particular root packet.

**Corollary 5 (universal quadratic-face three-window).** For every
characteristic-zero field, every r,z with z!=0, and every k>=1,
T_k,T_(k+1),T_(k+2) are not all zero.

By the cubic constraint (5), a hypothetical triple at middle index M
makes z a root of

 Pi_M(z)=2z^3+(6M+3)z^2+6M(M+1)z+M(M+1)(2M+1).

Let s=v_2(M(M+1))>=1. The four coefficient valuations are
s,s+1,0,1 in increasing index order, so the lower Newton polygon
has vertices (0,s),(2,0),(3,1). Thus v(z) is -1 or s/2.
The [negative branch theorem](QUADRATIC_FACE_FORCING_NEGATIVE_TWO_ADIC_BRANCH.md)
excludes -1, the
[positive branch theorem](QUADRATIC_FACE_FORCING_POSITIVE_VALUATION_BRANCHES.md)
excludes s/2>=1, and Theorem1 excludes the remaining value 1/2.
For the field transfer, Pi_M makes z algebraic over Q. The polynomial
T_k(r,z) has nonzero leading coefficient (2k)! in r, so a common root
r is also algebraic. Embed their number field into an algebraic closure
of Q_2 and apply the three branches.

The [earlier parametrized residual family](QUADRATIC_FACE_FORCING_THREE_ZERO_RESIDUAL_FRONTIER.md)
still has a direct pair-coprimeness proof. For M=q(4q+3)^2 and
z=-8q(q+1)(2q+1)+/-omega with omega^2=-3q(q+1), one has

 Pi_M(z)=(2(z+M)+4q+3)[(z+M)^2-2q(z+M)+q(4q+3)]=0.

For q=1 or2 modulo4, v(omega)=1/2 and the rational center has
valuation4. Corollary4 proves the pair-coprimeness assertion for
every such q. The universal three-window above uses the cubic
valuation alternatives alone.

Finally use the canonical full symmetric interior recurrence, Section2
Equation10 of [the symmetric quadratic proof](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md),
SHA256 c29933ff9bfee00d55703e049cde90c7ceb1cc577ada03a34dce4b9413530c00:

 M_(j+1)=(2a+b)(j+1)^2 M_j+j(j+1)(e^2-(2a+b)c)M_(j-1)
                          +H_(j+1)+(j+1)eH_j,

where H_j=Lambda_1((a t^2+e t+c)^j).

For f=a(x^2+y^2)+bxy+e(x+y)+c with ec!=0, a five-zero block starting at
n>=1 forces the three normalized face zeros T_(n+2),T_(n+3),T_(n+4).
Corollary5 excludes it. Thus the displayed symmetric family has no five
consecutive zero factorial moments, with no division by a,b,2a+b, a Hessian,
or e^2-(2a+b)c. In particular when b=0 and aec!=0 its actual five monomials
satisfy their Strong Factorial window. The full six-term symmetric family
has this stronger five-zero exclusion as well.

This does not settle a general nonsymmetric binary quadratic, full SFC(2),
or SFC(3), and proves no converse from a closed forcing assertion to the
unused monic anti-coefficient gate. The older conditional gate and formal
series investigations remain separate; no such hypothesis is seeded here.

## 8. Exact proof inputs and scope

Sections1--6 and the critical-half pair statement use factorial
coefficients and valuations directly. Corollary5 adds Theorem1 of the
[negative branch](QUADRATIC_FACE_FORCING_NEGATIVE_TWO_ADIC_BRANCH.md)
and Proposition1 of the
[positive branch](QUADRATIC_FACE_FORCING_POSITIVE_VALUATION_BRANCHES.md).
The negative branch uses the all-unit divided-power factorization
stated in its proof. The symmetric application additionally uses
the pinned polynomial interior recurrence from the
[symmetric quadratic proof](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md).
The monic anti-coefficient gate and the rational classification of
the residual cubic are separate questions.
