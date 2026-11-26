# No three quadratic face-forcing zeros on the negative 2-adic branch

Write v=v_2, normalized by v(2)=1, on an algebraic closure of Q_2.
For r,z define
h_j(r,z)=Lambda_1((r x^2+x+z)^j),
t_j(r,z)=h_j+j h_(j-1), j>=1,
where Lambda_1(x^d)=d!.

**Theorem 1 (negative branch).** If v(z)=-1, then for every k>=1
the three polynomials t_k,t_(k+1),t_(k+2) have no common zero r
in the algebraic closure of Q_2. No assumption about the parity
of the starting index is imposed.

This is only one valuation branch of the possible triple-zero forcing
problem. It does not by itself prove the universal three-window or
a new full binary quadratic Strong Factorial theorem.

The one imported local ingredient is the [all-unit divided-power
factorization](UNIVARIATE_QUADRATIC_ALL_SUPPORTS_STRONG_FACTORIAL.md), Section 1,
digest 9254fb9345762aa5a67e36dd3c01a87e5f44821edcde28bd7f197d68e1c0cfa4.
Its exact consequences used below are recorded in Section 3.

## 1. Triple zeros and the cubic constraint

The factorial identity Lambda(P')=Lambda(P)-P(0) gives, for j>=2,
h_j=2rj(2j-1)h_(j-1)
    +j(j-1)(1-4rz)h_(j-2)+z^(j-1)(z+j).                         (1)
For example, apply that identity to
phi' phi^(j-1), where phi=r x^2+x+z and
(phi')^2=4r phi+(1-4rz), and then to phi^j.
Only division by nonzero rational integers is involved.

Suppose a triple zero starts at k, and put m=k+1>=2.
For i=m,m+1 the adjacent zeros imply
h_i=-i h_(i-1), h_(i-2)=-h_(i-1)/(i-1).
Substitution in (1) gives
2ri(2i-1+2z)h_(i-1)+z^(i-1)(z+i)=0.                              (2)

Multiply the equation at i=m by (m+1)(2m+1+2z), the one at i=m+1
by 2m-1+2z, and add, using h_m=-m h_(m-1).
The h terms cancel. Dividing only by the nonzero z^(m-1) yields
Pi_m(z)=2(z+m)^3+3(z+m)^2+m=0.                                   (3)

Set q=z+m and w=2q. Since v(z)=-1, w is a unit and (3) gives
w^3+3w^2+4m=0.
Modulo 2 this polynomial is w^2(w+1), whose unit root 1 is simple.
Hensel lifting gives a unique unit root, lying in Z_2. Consequently
w and z=(w-2m)/2 actually belong to Q_2. This coefficient-field
fact is necessary below: the factor coefficients have INTEGER
valuations, even though the common root r may lie in an extension.

Also,
(2q-1)(q+1)^2=-(m+1).
Here v(q+1)=-1. If
D=2m-1+2z=2q-1 and t=v(m+1),
then
v(D)=t+2.                                                       (4)
In particular D is nonzero. Further v(z+m)=-1.

## 2. A sharper exact derivative identity

Differentiating with respect to r, the elementary Stein identity gives
2r^2 partial_r h_j
 =(2rj-1)h_j+j(1-2rz)h_(j-1)+z^j.                               (5)

For completeness, if U=Lambda(x phi^(j-1)), then
h_j-z^j=j(2rU+h_(j-1)), while
r partial_r h_j=j(h_j-U-z h_(j-1)).
Elimination of U is (5).

There is the following universal polynomial identity, before imposing
any zero hypotheses:
r^2 D partial_r t_m-z^(m-1)(z+m)^2
 =[rmD-(m+z)]t_m
  +m[(m+z)-rz(D+2)]t_(m-1), D=2m-1+2z.                          (6)

Indeed use (5) at m and m-1 and write
h_m=t_m-m h_(m-1),
h_(m-2)=(t_(m-1)-h_(m-1))/(m-1).
The result is
2r^2 partial_r t_m
 =(2rm-1)t_m+m(1-2rz)t_(m-1)
   -2rm h_(m-1)+z^(m-1)(z+m).
Equation (1) similarly gives
2rmD h_(m-1)+z^(m-1)(z+m)
 =t_m-m(1-4rz)t_(m-1).
Combining these two expressions proves (6).

At the common adjacent zeros, (6) becomes
r^2 D partial_r t_m=z^(m-1)(z+m)^2.                              (7)
The right side is nonzero, so r and the derivative are nonzero.
No additional root equation was used to obtain (6).

## 3. All-unit normalization and its coefficient inheritance

Put B=1/z, so B in Q_2 has v(B)=1, and write r=zA. Set
Q_m(A)=t_m(zA,z)/z^m
 =Lambda((A x^2+B x+1)^m)
   +mB Lambda((A x^2+B x+1)^(m-1)).

The coefficient of A^i in the first term, relative to the
constant-plus-square coefficient binom(m,i)(2i)!, equals
sum_(l=0)^(m-i) binom(m-i,l)(2i+1)_l B^l,
where (u)_l denotes a rising factorial. It is in 1+2Z_2.
For i<m the second term's relative coefficient is
B(m-i) sum_(l=0)^(m-1-i) binom(m-1-i,l)(2i+1)_l B^l,
which belongs to 2Z_2. It is zero for i=m.
Thus every coefficient of Q_m has exactly the baseline valuation
v([A^i]Q_m)=2i-s_2(m)+s_2(m-i).

Put s=s_2(m), and define
P_m(V)=2^s Q_m(V/4).
Then
v([V^i]P_m)=s_2(m-i), 0<=i<=m.                                  (8)
The leading coefficient is a unit. Divide by it to make P_m monic.

The reciprocal of this monic polynomial has all-unit divided-power
coordinates for eta_j(Z)=2^j Z^j/j!. The cited all-unit splitting,
valid for arbitrary m, therefore factors P_m into one monic
Eisenstein factor g_e for each binary degree e=2^r occurring in m.
The factors have distinct root valuations 1/e.
Moreover, writing
g_e(V)=sum_(j=0)^e a_j V^j, a_e=1,
the divided-power integrality of the corresponding unreversed
factor gives the precise inheritance bound
v(a_j)>=s_2(e-j), 0<=j<e,                                        (9)
with v(a_0)=1. To see this directly, the coefficient in question
is a unit multiple of d_(e-j)*2^(e-j)/(e-j)! with
d_(e-j) in Z_2. No odd-coefficient refinement is needed.

## 4. The forced derivative contradicts every possible binary factor

At the hypothesized common root set beta=4r/z.
By (8)--(9), beta is a root of exactly one g_e, with v(beta)=1/e.
For the unnormalized P_m, the derivative scaling in (7) is
P_m'(beta)=2^(s-2) z^(1-m) partial_r t_m(r,z).
Use v(r)=-1+v(beta/4)=-3+1/e, (4), and v(z+m)=-1 in (7).
The result is exactly
v(P_m'(beta))=s-2/e-t, t=v(m+1).                                 (10)

Write S=m mod e, ell=s_2(S), and let h be the number of set bits
of m above e. Thus s=ell+h+1. Every lower-degree factor g_d,
d<e, has value of valuation d/e at beta; every higher-degree
factor has value of valuation 1. Hence the other factors contribute
S/e+h, and
v(g_e'(beta))=ell+1-t-(S+2)/e.                                   (11)

The coefficient valuations of g_e are integers because g_e lies
in Q_2[V]. The nonzero derivative summands
j a_j beta^(j-1), 1<=j<=e,
have pairwise distinct fractional parts (j-1)/e.
Thus no minimum can cancel: (11) must be the valuation of the
unique summand with that fractional part.

First suppose S<e-1. The forced index is
j_0=e-S-1, between 1 and e-1. Comparing its valuation with (11)
gives
v(a_(j_0))=ell-t-v(j_0).                                         (12)
If S is even, then t=0 and j_0 is odd, so (12) is ell,
strictly below the required bound s_2(S+1)=ell+1.
If S is odd, then
t=v(S+1)=v(j_0)>=1,
and (12) is ell-2t, strictly below
s_2(S+1)=ell+1-t.
Both contradict (9).

Finally suppose S=e-1. This includes e=1. Writing e=2^r,
one has ell=r. The forced derivative index is now the leading
index e, whose actual valuation is
v(e beta^(e-1))=r+1-1/e.
But (11) is r-t-1/e, smaller by 1+t. This is impossible.

Every possible factor is excluded. This proves Theorem 1.

## 5. Scope and role

For the original forcing H_j+j e H_(j-1), e!=0, division of the
phase by e reduces to the r,z normalization here. Theorem 1 excludes
a hypothetical triple at any embedding for which v_2(c/e)=-1.
It does not exclude the positive-valuation roots of (3).

The proof above is all-index. Bounded native controls are supplementary;
no finite-to-universal inference is used. The positive-valuation exclusion
and the exact residual family are separate assertions, proved in the
[positive-branch proof](QUADRATIC_FACE_FORCING_POSITIVE_VALUATION_BRANCHES.md)
and the [residual-frontier proof](QUADRATIC_FACE_FORCING_THREE_ZERO_RESIDUAL_FRONTIER.md).
No universal three-window or full binary quadratic theorem is concluded here.
