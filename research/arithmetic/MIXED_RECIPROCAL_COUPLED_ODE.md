# Mixed reciprocal critical pairs: coupled equations and rational classification

The remaining reciprocal degree-26 cases admit exact coupled differential
equations. With a square reciprocal parameter they give seven-variable
rational systems, and the sign split (2,10) reduces further to four
variables. With a nonsquare parameter the conjugate (6,6) case also gives
seven rational variables, with descent built into the equations.

The same reduction completely classifies degree ten over Q: the known
normalized rational reciprocal model is unique. Mixed signs are impossible
globally, including when the reciprocal parameter is nonsquare. This
degree-ten conclusion uses a projective modulo-11 certificate, with no
integrality or good-reduction assumption on a hypothetical rational model.

No rational degree-26 model or global nonexistence theorem is proved here.
The systems below have explicit nondegeneracy exclusions. Their description
as finite polynomial systems does not assert that their rational-point sets
or their complex solution sets are finite.

## Normalization

Use the setting of the
[reciprocal obstruction note](PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md).
Put b=2h+1, n=4h+2, L=4h+1, and take monic N in Q[u] with
N(1)=N'(1)=0. Its critical polynomial is

\[
                     2uN'(u)-N(u)=(u-1)S_N(u).
\]

The residual roots are paired by u->lambda/u, lambda in Q*, with no
fixed roots. At those roots N^2/u is nonzero and takes precisely h distinct
values, twice each. Put rho^2=lambda, K=Q(rho), r=1/rho, and

\[
 F(v)=\rho^{-b}N(\rho v),\qquad
 \Sigma(v)=\frac{S_N(\rho v)}{L\rho^{2h}},\qquad w=v+v^{-1}.
\]

Then Sigma is monic reciprocal of degree 2h, with constant one, and

\[
                   2vF'-F=L(v-r)\Sigma(v).                  \tag{1}
\]

Write its two aggregate sign factors as

\[
 \Sigma_+(v)=v^pP(w),\qquad \Sigma_-(v)=v^qQ(w),\qquad p+q=h.
\tag{2}
\]

Here P,Q are monic polynomials of degrees p,q over K. A plus root satisfies
F(v)+vF(1/v)=0, and a minus root satisfies F(v)-vF(1/v)=0.
These are aggregate sign factors and need not be irreducible.

Every monic reciprocal even-degree polynomial with constant one has the
form used in (2): induction on its degree removes its highest symmetric
term v^j+v^-j. Thus s(w)=P(w)Q(w) is monic of degree h and
Sigma(v)=v^h s(w).

The rationality of lambda is automatic if such a reciprocal pairing exists
geometrically. The equal-value pairs are unique and hence Galois invariant.
Every paired product is lambda, so Galois fixes lambda. Thus a rational
model cannot evade the classification by taking lambda in an algebraic
extension. The same argument works over any base field with a separable
critical-root algebra.

## Theorem 1. Uniform coupled polynomial equations

Define C_j,D_j by

\[
 C_0=D_0=1,\quad C_1=w-1,\quad D_1=w+1,\qquad
 Z_{j+1}=wZ_j-Z_{j-1}.
\]

Put

\[
 \alpha_h=D_h-rD_{h-1},\qquad
 \beta_h=C_h-rC_{h-1}.
\]

For F satisfying (1), the expressions

\[
 X(w)=\frac{F(v)+vF(1/v)}{v+1},\qquad
 Y(w)=\frac{F(v)-vF(1/v)}{v-1}
\tag{3}
\]

are monic polynomials of degree 2h. They are exactly the polynomial
solutions of

\[
 \boxed{X+2(w+2)X'=L\alpha_h\,PQ,}
 \qquad
 \boxed{Y+2(w-2)Y'=L\beta_h\,PQ.}                              \tag{4}
\]

Both differential operators are invertible on polynomials in characteristic
zero. For a right side expanded as sum a_j(w+2)^j, the first inverse
divides a_j by 2j+1. The second inverse uses powers of w-2 in the same
way. Consequently (4) determines X,Y without integration constants,
root choices, or hidden Laurent-polynomial conditions.

Proof. The recurrence gives the identities

\[
 C_j(w)=\frac{v^{j+1}+v^{-j}}{v+1},\qquad
 D_j(w)=\frac{v^{j+1}-v^{-j}}{v-1}.
\]

For F=sum f_jv^j this makes (3) polynomial:

\[
 X=f_0+\sum_{j=1}^{2h+1}f_jC_{j-1},\qquad
 Y=-f_0+\sum_{j=1}^{2h+1}f_jD_{j-1}.
\]

Its degree and leading coefficient follow from f_(2h+1)=1. Conversely

\[
               F(v)=\frac{(v+1)X(w)+(v-1)Y(w)}2.
\tag{5}
\]

Direct differentiation gives

\[
 2(2vF'-F)=(v-1)\{X+2(w+2)X'\}
                       +(v+1)\{Y+2(w-2)Y'\}.
\]

The recurrence also gives

\[
 (v-1)\alpha_h+(v+1)\beta_h=2(v-r)v^h.
\]

Comparing in the basis 1,v of K(v)/K(w) proves (4). Each operator is
diagonal with nonzero eigenvalues 2j+1 in the indicated shifted-power
basis. This proves uniqueness. In particular, constructing F coefficient
by coefficient from (1) and constructing X,Y from (4) give exactly (3);
the converse does not presume that the Laurent expression (5) is
polynomial. It follows from that uniqueness. ∎

The coefficient reconstruction is explicit. If
Sigma=sum_(j=0)^(2h) s_jv^j, with s_-1=s_(2h+1)=0, then

\[
                  f_j=L\frac{s_{j-1}-rs_j}{2j-1}
                            \quad(0\le j\le2h+1).             \tag{6}
\]

Only nonzero rational integer constants are divided by in (4) and (6).

## Theorem 2. Exact equations for the square-parameter mixed cases

Suppose rho is rational. Choose monic P,Q of degrees p,q, p+q=h,
and r in Q*. Construct Sigma, F, X,Y as above. The exact equations are

\[
                 \operatorname{rem}_P X=0,\qquad
                 \operatorname{rem}_Q Y=0,\qquad F(r)=0.
\tag{7}
\]

Their unknowns are the h coefficients of P,Q and r. The two remainders
supply p+q=h scalar equations, and the last supplies one more.

They are equivalent to the required normalized rational reciprocal
models on the following explicit open set:

* s=PQ is squarefree and s(2)s(-2) is nonzero;
* s(r+r^-1) is nonzero;
* gcd(P,Y)=gcd(Q,X)=1;
* the h squared critical values described below are distinct and nonzero.

The last condition is algebraic and computable. In B=K[w]/s define m by
the Chinese remainder prescriptions

\[
 m\equiv (w-2)Y^2/4\pmod P,\qquad
 m\equiv (w+2)X^2/4\pmod Q.
\tag{8}
\]

The characteristic polynomial of multiplication by m must have nonzero
constant coefficient and nonzero discriminant. No assumption that P
or Q is irreducible is used.

Proof. At a root of Sigma_+, fixed-point-freeness excludes v=-1, so its
sign condition is exactly P|X. Likewise the other sign condition is Q|Y.
Since s is squarefree and avoids 2,-2, the roots of Sigma are simple
and avoid 1,-1. Equation (4) then improves the divisibilities to

\[
                            X=P^2U,\qquad Y=Q^2V,
\tag{9}
\]

where U,V are monic of degrees 2q,2p. For example at a simple P-root,
X=0 and (4) gives X'=0 because w+2 is nonzero.

At a plus root, (5) gives F=(v-1)Y/2; at a minus root it gives
F=(v+1)X/2. The gcd conditions are exactly nonvanishing of F on Sigma.
The corresponding values F(v)^2/v are (8). A squarefree characteristic
polynomial of degree h says that m separates the h geometric points of B.
This is precisely the last open condition.

Finally F(r)=0 and (1) give F'(r)=0, while

\[
                         2rF''(r)=L\Sigma(r)\ne0.
\]

Thus the zero has exactly order two. Setting
N(u)=rho^b F(u/rho), rho=1/r, gives monic rational N with
N(1)=N'(1)=0. The critical-point and critical-value conditions follow
from the four open conditions. Every original model gives exactly these
data by (2)-(3), proving both directions. ∎

The coupled equations after (9) can also be written

\[
 4(w+2)P'U+P\{2(w+2)U'+U\}=L\alpha_h Q,
\]
\[
 4(w-2)Q'V+Q\{2(w-2)V'+V\}=L\beta_h P.                         \tag{10}
\]

They are useful polynomial divisibility equations. A norm identity gives
an additional Pell-type form:

\[
       (w+2)X^2-(w-2)Y^2=4F(v)F(1/v).
\tag{11}
\]

The right side is a polynomial in w of degree 2h+1 with leading
coefficient 4Lr. The double zero at r implies that it is divisible by
(w-r-r^-1)^2, including r=1 or -1. Consequently for monic Z of degree
2h-1,

\[
 (w+2)P^4U^2-(w-2)Q^4V^2
                  =4Lr(w-r-r^{-1})^2 Z(w).
\tag{12}
\]

This is an identity of Pell form, with a variable right side; no
classification of its solutions by ordinary polynomial Pell theory is
claimed.

For degree 26, the previously proved uniform-sign obstruction excludes
p=0 or q=0. Up to changing the sign of rho, the square-parameter cases
are (p,q)=(1,5),(2,4),(3,3). Equations (7) give seven variables and
seven equations in each case, together with the displayed open conditions.

## Theorem 3. Built-in rational descent for the nonsquare case

Let lambda in Q* be nonsquare, K=Q(rho), rho^2=lambda, and gamma its
conjugation. A mixed normalized rational model must have h even and

\[
 p=q=h/2,\qquad Q(w)=(-1)^pP^\gamma(-w).
\tag{13}
\]

Conversely choose monic P of degree p over K, define Q by (13), and
put r=1/rho. Then equations (7) reduce exactly to P|X and F(r)=0.
The former gives p equations over K; the latter is a single equation
over Q after division by rho. Reconstruction gives N in Q[u].

Writing the p coefficients of P as a_j+b_jrho gives h rational
coefficients and the additional rational parameter lambda. Thus there
are h+1 rational unknowns and h+1 rational equations, with lambda
nonsquare and the same nondegeneracy conditions. Degree 26 has p=3,
so this is a seven-variable rational system for the nonsquare (6,6)
case.

Proof. Conjugation of rho swaps the two original signs at the same
u, so their degrees are equal. Each degree is even; hence h is even.
After the coordinate change u=rho v, conjugation changes v to -v.
Monicity gives exactly (13).

Conversely (13) gives s^gamma(w)=s(-w). Equation (1), or its unique
coefficient solution (6), then implies

\[
 F^\gamma(v)=-F(-v),\qquad X^\gamma(w)=Y(-w).
\]

Conjugating P|X therefore gives Q|Y. Moreover

\[
                F(r)^\gamma=F^\gamma(-r)=-F(r).
\]

Its rational component is identically zero. The equation F(r)=0 is
exactly its one remaining component. Finally

\[
 \bigl(\rho^bF(u/\rho)\bigr)^\gamma
                 =(-\rho)^b\{-F(u/\rho)\}
                 =\rho^bF(u/\rho),
\]

because b is odd. This proves descent and both directions of the
equivalence. Coefficients can be reduced using rho^2=lambda and
denominators cleared with lambda nonzero; no further descent
conditions are omitted. ∎

This statement concerns a quadratic field, so lambda nonsquare is an
arithmetic restriction, rather than a Zariski-open condition. Neither
squareness nor integrality at three is silently assumed.

## Theorem 4. Four variables for the degree-26 sign split (2,10)

More generally let p=1,q=h-1. Write P=w-c and V=w^2+aw+b in (9).
The second equation (10) becomes

\[
 4(w-2)VQ'+\{V+2(w-2)V'\}Q
                         =L(C_h-rC_{h-1})(w-c).               \tag{14}
\]

There is a unique monic degree-q Q satisfying all coefficients of
degrees 2 through q+2 of (14). Construct it from highest degree
downwards. The coefficient of w^(j+2) contributed by an unknown
q_jw^j is (4j+5)q_j. Hence for j=q-1,...,0, successively divide the
current error coefficient by -(4j+5). All divisors are nonzero
constants. The top coefficient already agrees, since 4q+5=L.

After this recursion the error has degree at most one. Set its two
coefficients to zero and append X(c)=0 and F(r)=0. These four equations
in a,b,c,r are equivalent to all square-parameter models of this sign
split, on the open set in Theorem 2.

Proof. Equation (14) is linear in Q, and its leading coefficient on
each monomial is the stated nonzero constant, proving the recursion
and uniqueness without a generic-invertibility assumption. If the
two remaining coefficients vanish, Y=Q^2V solves the second equation
(4), hence is its unique polynomial solution. Thus Q|Y is already
built in; X(c)=0 and F(r)=0 are exactly the remaining equations (7).
The converse follows by factoring the monic Y as in (9). ∎

The complete integer polynomials for degree 26 are retained in
[degree26_mixed_2_10_system.json](degree26_mixed_2_10_system.json).
The divisions for q=5 are 21,17,13,9,5. They introduce no variable
denominators. The four stored equations have total degrees 6,6,11,12:
the third has been divided by its automatic factor c+2, which is
nonzero on the required open set; the fourth is F(r)/(2r).
Clearing rational constant denominators makes these integer
polynomials with exactly the same characteristic-zero zero set.

No saturation or complete characteristic-zero elimination of this
degree-26 system has been proved here. In particular its equation
count is not a finiteness proof.

## Theorem 5. Global uniqueness in normalized rational degree ten

The monic rational degree-ten reciprocal model from
[RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)
is the unique model satisfying the present normalization and
nondegeneracy conditions.

First, mixed signs are impossible even if lambda is nonsquare.

Proof. Here h=2 and necessarily p=q=1. Write P=w-c,Q=w-d over K.
On the good open set c,d avoid 2,-2 and r is nonzero. Put

\[
 p_0=\frac{35X(c)}{2(c+2)},\qquad
 q_0=-\frac{35Y(d)}{2(d-2)},\qquad
 z_0=-\frac{35F(r)}{2r}.
\]

These are polynomials, and all three vanish. Let
sigma=c+d, t=cd, delta=c-d. Direct substitution gives

\[
 z_0=21tr^2-9\sigma r^3-105\sigma r+5r^4+42r^2-315,
\]
\[
 p_0-q_0=A+D\delta,\qquad p_0+q_0=C+B\delta,
\]

where

\[
\begin{split}
 A&=-5\sigma^3+24t\sigma+9r\sigma^2-60rt-12\sigma+240r,\\
 B&=-5\sigma^2+14t+9r\sigma-42,\\
 C&=15\sigma^2-84t+36r\sigma+112,\qquad
 D=15(\sigma-6r).
\end{split}
\]

Since r is nonzero, z_0=0 determines

\[
             t=\frac{\sigma r(9r^2+105)-5r^4-42r^2+315}{21r^2}.
\]

Substitute this in

\[
 E_1=AB-DC,\qquad
 E_2=A^2-D^2(\sigma^2-4t),\qquad
 E_3=C^2-B^2(\sigma^2-4t),
\]

and clear nonzero constant and r denominators. These three polynomials
vanish. No division by D or B has occurred, so the cases sigma=6r
and B=0 have not been discarded.

The three pairwise Sylvester resultants in sigma have degrees 68,80,96
in r. Their greatest common divisor is a nonzero constant times

\[
                               r^{24}E_{18}(r^2).
\]

The exact nineteen coefficients of E18 and all three intermediate
polynomials are retained in the
[receipt](mixed_reciprocal_coupled_ode_verification.json).
Its leading coefficient is 20480, and its values modulo 11 at
0,1,...,10 are

\[
                         (7,8,4,3,2,2,10,7,6,3,4).
\]

The leading coefficient is 9 modulo 11. Thus its homogenization has
no point in P1(F11). It has no rational root with any denominator.
But r^2=1/lambda is rational in both the square and nonsquare cases,
a contradiction.

The checker also recomputes a complete lexicographic Groebner basis
of p_0,q_0,z_0 in variables c,d,r. It is zero-dimensional, has three
basis polynomials, and its final polynomial is E18(r^2)/20480.
The basis is retained in full. The independent resultant necessity
argument supplies a second route to that same obstruction.

For uniform signs, rho is rational by the earlier conjugation argument.
The earlier exact classification gives

\[
 F(r)=\frac{2r(2r-15)(4r^3-6r^2-9r-15)}{25}.
\]

The cubic has no projective root modulo seven: its finite values are
(6,2,3,5,4,3,5) and its leading coefficient is four. Thus the only nonzero rational solution of F(r)=0 is r=15/2. Formula (6) determines the
numerator uniquely, giving exactly the retained known model. That
model satisfies every open condition by its existing exact
certificate, so uniqueness is an existence-and-classification
statement for this ansatz. ∎

This is a global arithmetic statement. It does not assume integral
coefficients at any prime, and it does not classify all rational
alternating pencils outside reciprocal critical pairing.

## Verification and remaining degree-26 problem

The [checker](verify_mixed_reciprocal_coupled_ode.py) and
[receipt](mixed_reciprocal_coupled_ode_verification.json) verify both
coupled equations symbolically in h=1,2,3,6, the polynomial inverses,
the norm identity, generic quadratic descent in degree 26, the exact
four-variable recursion, and the complete degree-ten elimination
and projective obstructions.

The all-degree equivalences have written proofs above. The finite
checks supplement them. The degree-26 equation file is a complete
rational search reduction for its stated split and open conditions;
it contains no asserted rational point or global nonexistence result.

The remaining square splits (4,8) and (6,6), and nonsquare conjugate
(6,6), are given without omitted descent conditions by Theorems 2-3.
A rational point satisfying all open conditions would give a rational
degree-26 paired pencil. Its full alternating monodromy and literal
square discriminant would then follow from the retained passport
theorems, since 25 is a rational square.

The [source ledger](mixed_reciprocal_coupled_ode_source_ledger.json)
pins the exact internal proof inputs and the new evidence. No novelty
priority claim, formal verification, independent external replay or
external human review is made.
