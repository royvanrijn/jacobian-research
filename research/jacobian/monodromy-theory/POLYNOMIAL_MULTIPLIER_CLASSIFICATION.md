# Polynomial multipliers are finite: completing the k(w)[z] sweep chart

Work geometrically over algebraically closed characteristic-zero constants.
Use S=X(w,z)+gamma(1,w,z), det JS=gamma^2 and nonzero polynomial
normalization lambda*[[-rho,1],[-rho^2,rho]].

**KMFI44.** If rho=h(w) is a nonconstant polynomial, every polynomial
wave solution lambda gives a **finite** sweep S.
There is no nonfinite class left under these hypotheses.
Write d=deg h>=1, s=deg_z lambda and
top_z(lambda)=(alpha*w+beta)z^s.
Put epsilon=1 if alpha!=0 and epsilon=0 otherwise.
The exact degree is

\[
 N=\begin{cases}
 (s+1)((s+2)(d+1)-1),&\epsilon=0,\\
 (s+2)((s+1)(d+1)+1),&\epsilon=1.
 \end{cases} \tag{1}
\]

The wave solutions are parametrized by two integration constants
at each z coefficient, with the finite recursion below.

**KMFI45.** Apart from the following explicit exception, the full
group is natural A_N. The exception occurs precisely when
lambda=alpha*(w-c), alpha!=0, and h is even about c.
Write h(c+v)=a(v^2), deg a=m>=1.
The degree is 4m+4 and the group in its actual action is

- SL2(F3) on the eight nonzero vectors of F3^2, if m=1;
- V_n:A_n on 2n signed letters, n=2m+2, if m>=2,
  where V_n is the even binary sign module.

For m>=2 this extends the KMFI29 full-group argument to arbitrary
coefficients of a. The degree-eight exception is a finite sweep field action.
[KMFI60](BRANCH_CENTERED_COMPLETION.md) now completes this exact
field on ordinary affine4 over Q(sqrt(-3)).
In that exceptional case the arithmetic group over a coefficient
field K is SL2(F3) or GL2(F3), according as -3 is or is not a square
in K, and the exact constants are K(sqrt(-3)).
The displayed centered degree-eight presentation has no unramified
origin and cannot supply its positive-weight Gao stage.

**KMFI46.** Combining KMFI33,38,43--45 completes the possible
geometric actions for every chart rho in k(w)[z]:

- natural S_n, n>=1;
- natural A_n, n>=3, n!=4;
- SL2(F3) on eight nonzero vectors;
- V_n:A_n on 2n signed letters, every even n>=6.

Every action in this list occurs as a sweep field.
Neither specified PSL2(11) action, natural A4, nor any other proper
primitive action outside the listed alternating/symmetric range
occurs through this entire chart. The field-pair exclusions survive
independent stabilization and arbitrary birational recharts.
This is a constructor classification, not a complete classification
of polynomial Keller maps.

Natural S2 and A3 are excluded as actual Keller actions by KMFI1.
The listed S_n for n>=3, A_n for n>=5 and signed family for even
n>=6 already have actual ordinary-affine-space Keller witnesses
elsewhere in the theory. KMFI60 now supplies the actual degree-eight
SL2(F3) realization by affine branch centering.
Its cycle-gcd invariant is4<8, excluding every independently
stabilized curve-induced Keller route for that action.
Its four-pair quotient A4 is not an affine-space Keller realization
of natural A4.

[Formula module](polynomial_wave.py) · [checker](verify_polynomial_chart.py) ·
[model](polynomial_chart_model.json) · [receipt](polynomial_chart_receipt.json).

## 1. All polynomial wave solutions

With h independent of z, the wave equation is

\[
 \lambda_{ww}+2h\lambda_{wz}+h^2\lambda_{zz}
                              +2h'\lambda_z=0.
\]

Write lambda=sum_b f_b(w)z^b. Its top coefficient satisfies
f_s''=0, so f_s=alpha*w+beta.
Choose arbitrary constants alpha_b,beta_b at each lower level and
integrate twice

\[
 f_b''=-2(b+1)(h f_{b+1})'-(b+1)(b+2)h^2 f_{b+2}. \tag{2}
\]

This produces every polynomial solution; there is no denominator
condition or unresolved equation.
The polynomial X is recovered from the closed derivative pairs
in KMFI41, uniquely up to translation.

Set j=d+1, n=s*j+epsilon, a=n+j-1.
For source weights (gamma,w,z)=(a,1,j), the output weights are
(a,a+1,a+j).
Every lower integration constant in (2) has smaller weight than n;
a lower coefficient of h does too.
After scaling its nonzero leading coefficient, the leading model
is h=w^(j-1), lambda=w^n L(z/w^j), with L monic of degree s.
The weight bounds for X follow directly from its derivative pairs.

## 2. The leading model is finite

Put t=z/w^j, A=1-j*t.
The leading wave equation is

\[
 A^2 L''+B L'+n(n-1)L=0,
\]
\[
 B=2(n-1)+(-2jn+j^2+j)t,\qquad B(1/j)=j-1\ne0.
\]

Thus L has no root at1/j and all its other roots are simple:
a multiple root at an ordinary point contradicts uniqueness of
the linear differential equation; at1/j the lowest term B L'
cannot cancel.
Euler differentiation gives

\[
 a X_1=A^2L'+(nA+j-1)L.
\]

As in KMFI42, a common zero of the leading V=X2-wX1 and
W=X3-zX1 at w!=0 would give a curve over the origin.
The Jacobian identity then forces gamma=0 and X=0 along it.
If L=0, simplicity and A!=0 make X1 nonzero.
If L!=0, applying the normalization matrix to the weighted ray
direction forces A=0, and the displayed Euler identity again
makes X1 nonzero. This rules out such a common zero.

At w=0, epsilon0 gives X1=0 and
X2=z^(s+1)/(s+1); hence an origin preimage has z=0.
For epsilon1, X1=-z^(s+1)/(s+1),
X3=-z^(s+2)/(s+2).
The first target equation sets gamma=z^(s+1)/(s+1), and the
third is z^(s+2)/((s+1)(s+2)), again forcing z=0.
Then gamma=0. The leading model therefore has only the origin
over the origin.

All weights are positive. Its source graded ring is finite over
the subring generated by the leading target coordinates.
Lifting a finite homogeneous generating set and inducting on weight
makes the original polynomial source ring finite over k[S].
Thus every lower-weight deformation is also finite.
This proves the finiteness assertion for every polynomial h and
every solution of (2), not merely the monomial controls.

For a finite polynomial map with independent leading coordinates,
the ratio of weighted Hilbert growth coefficients gives the degree:
a(a+1)(a+j)/(a*1*j)=(a+1)(a+j)/j.
Substitution of a=(s+1)j+epsilon-1 gives (1).
This also shows why the apparent nonfinite one-parameter frontier
was empty rather than a missing family of actions.

## 3. The alternating case and its precise exception

Rank two follows from KMFI38: rank-one normalization image would
be a constant line, forcing h constant.
Now KMFI28 applies because finiteness has been proved.
It gives either delta1 with full A_N, or the unique delta2 normal
form, after translating the root of lambda and scaling:

\[
 q=w^2,\quad h=a(q),\quad \lambda=w,\quad
 X=(\sigma,b(q),d_0(q)-\sigma^2/2),
\]
\[
 \sigma=a(q)w-z,\quad b'=q a',\quad d_0'=a b'.
 \tag{3}
\]

Thus delta2 requires s0 and a nonconstant affine lambda, and h is
even about its root. Conversely those hypotheses give (3).
The coprime degrees deg b=m+1 and deg d0=2m+1 show
k(b,d0)=k(q), so delta is exactly2.

For example h=w^2,lambda=w has the binary degree-eight group,
while h=w^2+w,lambda=w has A8 at the same degree. Finiteness and
degree persist under that deformation; the generic group need not.

For m>=2, put n=2m+2.
The eliminant and derivative-square identity of KMFI29 apply to
every polynomial a of degree m. Coprime degrees make the block
group A_n. The sign kernel contains the central flip from infinity,
and the binary-module and attributed Carmichael-presentation proof
forces the whole V_n. That argument requires n>=6 but no special
interpolation coefficients of a.
The exact arithmetic square classes likewise remain
(-1)^m(2m+1) and -(2m+1), over its coefficient field,
as in the resultant/product calculation of KMFI29.

## 4. The exceptional degree-eight action

For m1, translations and scaling reduce the field pair to

\[
 Q(q)=q^4-12Aq^2+24Bq-12A^2=0,\qquad q=w^2. \tag{4}
\]

The full w-curve, over the algebraic closure of k(A), has two
finite branch permutations of type3^2*1^2 and pole inertia(6,2).
It embeds in the even signed group V4:A4 of order96.
The exact finite enumeration retained by the checker constructs
all96 permutations, all32 elements of type3^2*1^2, and all
ordered transitive pairs with product type(6,2).
There are192 pairs; each generates a group of order24.
Their eight distinct subgroups are conjugate in the full signed
permutation group of order384.
An explicit permutation conjugacy identifies one with SL2(F3)
on its eight nonzero vectors. This is an exhaustive finite
action proof, not a specialization guess.

The curve group is normal in the full plane-field group over k(A,B).
The sign kernel of the latter is either the central flip or all V4:
every noncentral vector has weight2 and its A4 orbit spans V4.
Hence the plane group has order24 or96.
The order96 alternative is impossible. Its commutator subgroup
contains V4 and its abelianization is C3, so it has no abelian
quotient of order4 and no normal subgroup of index4.
The full geometric group is therefore the same SL2(F3) of order24.

Its normalizer in S8 is GL2(F3), of order48.
Here is also a size proof: the centralizer has order2, while
Aut(SL2(F3)) injects into Aut(Q8), of order24, by restriction to
the unique Sylow2 subgroup; an automorphism trivial there cannot
change an order3 generator by the central involution.
The explicit GL2(F3) action supplies the matching lower bound.
The block quartic discriminant has square class -3.
Thus the arithmetic group is order24 or48 exactly as stated,
and its sole possible constant extension is K(sqrt(-3)).

## 5. Whole-chart spectrum and the remaining degree-eight gate

A multiplier independent of z but rational in w must be polynomial
by KMFI33. Nonconstant polynomial h is now finite by KMFI44.
Constant h has rank at most one and full natural symmetric group
by KMFI38. A z-dependent multiplier in k(w)[z] is classified
by KMFI43. These exhaust the chart.

All natural S_n occur as rank-one pencils.
Taking lambda1 and h of degree d gives A_(2d+1), every odd n>=3.
Taking lambda=w and a degree-d h not even about zero gives
A_(2d+4), every even n>=6.
No formula gives natural A4.
The even h exception supplies precisely the signed actions and
the degree-eight SL2(F3) action above.
This proves the exact list in KMFI46.

The degree-eight action has cycle types
1^8,2^4,3^2*1^2,4^2,(6,2), so its cycle-gcd lcm is4.
KMFI10 therefore excludes every independent curve-induced Keller
realization of this action, without excluding other field pairs.

For the centered presentation h=w^2,lambda=w,
b=w^4/2,d0=w^6/3, the origin equations give

d0-h*b+b^2/(2w^2)=-w^6/24.

An origin preimage has w=0 and gamma=0.
The centered positive-weight stage requiring gamma!=0 cannot
complete this presentation.
Translations or a different affine presentation are not excluded.
KMFI60 resolves the binary action by a translated branch-centered
presentation. Natural A4 and the absolute PSL2(11) actions remain open. A block quotient is not an affine
Keller factor, and this finite sweep is not itself a Keller map.

All universal finiteness, degree and group arguments are written.
Exact controls are not formal verification or independent whole-proof review.
