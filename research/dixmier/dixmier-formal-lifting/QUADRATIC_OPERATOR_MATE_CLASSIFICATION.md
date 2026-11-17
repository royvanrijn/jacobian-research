# Quadratic scalar operators and every possible Weyl mate

14 November 2026. **QSO2.** Let F be a one-variable algebraic function
field over a characteristic-zero constant field, with a nonzero
derivation delta. In the scalar differential-operator ring F[D;delta],
suppose ord(P)=2 and
\[
                         [P,Q]=\kappa,\qquad
                    0\ne\kappa\in\ker(\delta).                 \tag{1}
\]
Then Q has order one or positive even order. More precisely,
\[
             Q=A(P)+Q_1,\qquad A\in\ker(\delta)[t],
                         \quad\operatorname{ord}(Q_1)=1.        \tag{2}
\]
After a finite algebraic coefficient extension and an affine change of
the differential generator making P=L=E^2+u, every mate has the form
\[
          Q=A(L)+aE+b,\qquad a,b\text{ constant},\quad
                    a\ne0,\quad u'=-\kappa/a.                  \tag{3}
\]
Thus no odd partner order at least three is possible, on any algebraic
coefficient curve or finite cover, in any coefficient degree.

There is a stronger parameter consequence. If P and Q are regular over
a differential valuation ring, P has order two with unit leading
coefficient, and their reductions have orders two and any odd number
at least three, then (1) is impossible regardless of Q's generic order.
Applied to a spatially bounded full Long lift with exact R, this excludes
every finite partner while S keeps momentum order two. A bounded lift
in that class must increase S's order; this note does not exclude that
increase or prove a bounded normalization fixing R.

The proof covers the whole first Painleve hierarchy with arbitrary
constant lower times. It uses formal differential algebra and Puiseux
series, not an extrapolation from computations in several orders.
The order-three case is included in this theorem.

## Normalization and the universal differential polynomials

Extend the constants to an algebraic closure C and choose a component
of the coefficient field after constant extension. The constants of a
nonzero derivation on a one-variable function field over C are exactly C.
Derivations extend uniquely to algebraic field extensions.

Writing P=a_2 D^2+... and Q=b_m D^m+..., the coefficient of order m+1
in the commutator is 2a_2 delta(b_m)-m b_m delta(a_2). Adjoin
rho=sqrt(a_2), set E=rho D and replace its derivation by epsilon=rho delta.
Then P is monic, and Q's leading coefficient b_m/rho^m is constant.
The coefficient-preserving Ore change E -> E-ell_1/2 removes the linear
coefficient ell_1 of P. Hence we may work with L=E^2+u and primes
meaning epsilon. No assumption about an algebraic primitive of 1 is used.

Use the differential-polynomial algebra C[u,u',u'',...] with derivative
partial. Define R_j by B(s)=sum_{j>=0} R_j s^j, R_0=1/2, and
\[
       4B^2-4usB^2-2sBB''+s(B')^2=1.
\]
At each positive s-degree the new coefficient is 4R_j plus known
terms, so this defines differential polynomials uniquely. Differentiation
and cancellation of the unit 2B give the Lenard identity
\[
       R_{j+1}'=\tfrac14R_j'''+uR_j'+\tfrac12u'R_j.              \tag{4}
\]
For example R_1=u/4 and R_2=(u''+3u^2)/16. Each R_j has weight 2j
for wt(u^{(i)})=i+2, with
\[
 R_j=4^{-j}u^{(2j-2)}+\cdots,\qquad
 [u^j]R_j=\beta_j=\frac{\binom{2j}{j}}{2\,4^j}\ne0
                 \quad(j\ge1).                               \tag{5}
\]
The first coefficient follows inductively from (4); setting all
derivatives to zero gives B=1/(2sqrt(1-us)) and the second coefficient.

Define the finite monic odd operators
\[
 B_{2j+1}=\sum_{\ell=0}^j(2R_\ell E-R_\ell')L^{j-\ell}.
\]
Direct Ore multiplication gives
\[
 [L,2RE-R']=4R'L-(R'''+4uR'+2u'R).
\]
Using (4), the sum telescopes:
\[
                         [L,B_{2j+1}]=-4R_{j+1}'.               \tag{6}
\]

We will also use two formal-algebra facts and justify them here.
In the pseudodifferential algebra, powers are bounded above in E and
multiplication is E^n f=sum_{i>=0} binom(n,i)f^{(i)} E^{n-i}.
Each individual coefficient is a finite sum. The unique monic square
root of L has leading term E and differential-polynomial coefficients.
Then
\[
                  B_{2j+1}=(L^{j+1/2})_+,\qquad
                  \operatorname{res}_E L^{j+1/2}=2R_{j+1}.      \tag{7}
\]
For the first equality, the positive part has monic order 2j+1,
homogeneous weight 2j+1 and scalar commutator with L: the commutator
of L with its negative part has order at most zero. Our telescoping
operator has the same properties. If their difference were nonzero,
its leading coefficient would have derivative zero, whereas its
strictly positive weight forbids a nonzero constant. Inducting on
order proves equality. In degree zero the same argument applies.
For the residue identity, the scalar term of [L,c E^{-1}+...] is 2c';
compare with (6). Homogeneity removes the integration constant.

The hierarchy vector fields X_j=4partial R_{j+1}, characterized by
delta_j L=[B_{2j+1},L], commute. Indeed delta_i commutes with partial,
and differentiating the unique monic fractional powers gives
delta_i L^{j+1/2}=[B_{2i+1},L^{j+1/2}]. Projecting the identity that
two powers of L commute gives
delta_i B_{2j+1}-delta_j B_{2i+1}-[B_{2i+1},B_{2j+1}]=0.
The Jacobi identity now gives [delta_i,delta_j]L=0.

Finally, the Frechet derivative DR_j is formally selfadjoint.
The formal trace Tr(A) is res_E(A) modulo total derivatives; direct
binomial multiplication shows Tr([A,B])=0. Put M=L^{1/2}. Cyclicity
applied to variations of M^{2j+1} and L=M^2 yields
\[
 \delta\operatorname{Tr}L^{j+1/2}
   =(j+\tfrac12)\operatorname{Tr}(L^{j-1/2}\delta u)
   =(2j+1)\int R_j\,\delta u.
\]
The two second variations commute; integration by parts therefore
identifies DR_j with its formal adjoint. These are universal identities
of differential polynomials and remain valid at every specialization.

The normalization here is half the residue normalization in
[Takasaki, Sections 2--3](https://arxiv.org/html/nlin/0610073v2).
The [2026 revised correspondence paper, Appendix A](https://arxiv.org/html/2508.18433v2#A1)
also gives the pseudodifferential and string-equation framework.
The arguments above fix the constants and signs used in this note;
no analytic Painleve-property theorem is imported.

## Every odd mate supplies a forced Painleve hierarchy

For any finite Q with [L,Q] scalar, its highest positive-order
coefficient is constant, by the next coefficient of its commutator.
Subtract that constant times L^{m/2} at even order m, or B_m at odd
order m. Both subtractions retain a scalar commutator. Repetition
ends in a scalar coefficient whose derivative must vanish. Hence
\[
 Q=A(L)+\sum_{j=0}^g c_j B_{2j+1},\qquad c_j\in C,              \tag{8}
\]
where a nonzero commutator requires at least one odd summand.
If the remaining highest odd order is 2g+1, take c_g nonzero.
Equations (1) and (6) imply
\[
       \tau=-\frac4\kappa\sum_{j=0}^g c_j R_{j+1}(u),
                       \qquad \tau'=1.                        \tag{9}
\]
Thus the relation itself creates an algebraic time coordinate. The
coefficient field is finite over C(tau). For g>=1, an algebraic u
would solve
\[
      \sum_{j=0}^g c_jR_{j+1}(u)=\alpha\tau+\beta,\qquad
                   c_g\ne0,\quad\alpha\ne0.                    \tag{10}
\]
Allowing beta and arbitrary lower constant c_j will be useful below.

## Integer indicial exponents at every hierarchy pole

Write z=tau-tau_0. A pole u~a z^p with p<0 in (10) must have p=-2.
For p<-2 the pure power u^{g+1} is uniquely most singular. For
-2<p<0 the highest linear derivative u^{(2g)} is uniquely most
singular. To see both assertions in every order, a weight-2j monomial
with m factors of u and its derivatives has exponent
m(p+2)-2j. Its extreme is m=j or m=1 respectively.

At p=-2 put
\[
 R_j(az^{-2})=C_j(a)z^{-2j},\qquad
 C_j(a)=\beta_j\prod_{\ell=0}^{j-1}(a+\ell(\ell+1)).
\]
This follows from (4), which gives
C_{j+1}=(2j+1)(a+j(j+1))C_j/(2j+2), starting with C_0=1/2.
The leading coefficient of (10) forces
\[
                         a=-k(k+1),\qquad1\le k\le g.           \tag{11}
\]

Define d_j(a,rho) by linearizing at this pole:
\[
 DR_j[az^{-2}](z^{\rho-2})
                         =d_j(a,\rho)z^{\rho-2j}.
\]
Here d_1=1/4. Linearizing (4) gives the exact recurrence
\[
\begin{split}
(\rho-2j-2)d_{j+1}
 &= (\rho-2j-1)
       \left(\frac{(\rho-2j)(\rho-2j-2)}4+a\right)d_j\\
 &\quad+\frac{\rho-4j-2}{2}C_j(a).                             \tag{12}
\end{split}
\]
In particular the last numerator is rho-4j-2, not rho-2j-2;
the variation of both u R_j' and u' R_j/2 contributes.

We prove all the needed roots are integers without assuming a
closed resonance formula. At u_k=-k(k+1)z^{-2}, X_k=0, while for
0<=m<k the lower field X_m is a nonzero constant times z^{-2m-3}.
Commutativity of X_k and X_m gives
\(\partial(DR_{k+1}[u_k]\,X_m)=0\). By homogeneity the expression being
differentiated has strictly negative z-degree, so it is zero.
Consequently d_{k+1} has the k roots -1,-3,...,-(2k-1).

Set T=DR_{k+1}[u_k], an operator of order 2k with leading coefficient
4^{-(k+1)}. Homogeneity gives T=z^{-2k}p(z partial), with
p(s)=d_{k+1}(a,s+2). Formal selfadjointness and
(z partial)^*=-z partial-1 imply
p(s)=p(2k-1-s). The roots of d_{k+1} are therefore paired by
rho -> 2k+3-rho. The k reflected roots and k negative roots are
distinct and exhaust its degree 2k. Hence
\[
 d_{k+1}(-k(k+1),\rho)
   =4^{-(k+1)}
       \prod_{m=0}^{k-1}(\rho+2m+1)(\rho-(2k+4+2m)).             \tag{13}
\]
For j>=k+1, C_j(-k(k+1))=0, so (12) becomes
\[
 d_{j+1}
   =\frac{(\rho-2j-1)(\rho-2j-2k-2)(\rho-2j+2k)}
                 {4(\rho-2j-2)}\,d_j.                         \tag{14}
\]
Every displayed zero or possible pole is at an integer. Thus
d_{g+1}(-k(k+1),rho) is nonzero for every noninteger rho.
No convergence or full analytic Painleve property has been assumed.

## Algebraic solutions cannot exist

Suppose u is algebraic over C(tau) and satisfies (10), with g>=1.
Use a Puiseux expansion at a finite tau_0. At a pole, (11) gives
the leading term. If its first nonintegral exponent is rho-2,
then rho>0 is noninteger. At the first fractional order in (10),
the highest hierarchy term contributes
c_g b d_{g+1}(-k(k+1),rho)z^{rho-2g-2}, with b nonzero.
Earlier integer terms cannot cancel it, lower hierarchy terms
have larger exponents, and terms using more than one fractional
correction also have larger exponents. Equations (13)--(14)
make its coefficient nonzero, a contradiction.

At a finite point where u is finite, a first noninteger exponent mu>0
contributes through 4^{-(g+1)}u^{(2g)}. The falling factorial
mu(mu-1)...(mu-2g+1) is nonzero. All nonlinear derivative monomials
and all lower hierarchy terms have strictly larger first fractional
exponents. This is again impossible.

Every branch at every finite point therefore has integral exponents.
The minimal function field C(tau,u) is unramified at every finite
tau, since its generating element u has all its branches in C((z)).
If its degree is d, its genus is G and it has r points over infinity,
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) gives
2G-2=-d-r. This forces d=r=1 and G=0. Thus u is rational.

A rational u with positive integer pole degree m at infinity has
u^{g+1} as the unique term of highest degree, namely (g+1)m>1.
Indeed for a weight-2j monomial its degree is m times its number
of factors minus its total derivative order, strictly smaller than
jm unless it is the pure power. Lower hierarchy terms also have
smaller degree. If u is bounded at infinity, every term on the left
of (10) is bounded. Neither case can equal alpha*tau+beta with
alpha nonzero. This proves the all-hierarchy algebraic exclusion.

In (8), therefore, only c_0 can be nonzero. Comparing [L,aE+b]
directly gives a'=b'=0 and -a u'=kappa, proving (3).
For the original P, subtract its constant-leading even powers
successively; an odd order at least three is impossible by what
was just proved. A scalar remainder cannot have nonzero scalar
commutator with an order-two operator. The remainder has order one,
which proves (2) over the original constant field, without needing
descent of the normalizing cover.

Even orders genuinely occur. In C(v)[D], [D,v]=1, take
P=D^2+v and Q=-D+A(P). Their commutator is one. Choosing A=P^m
gives order 2m for any m>=1. Consequently the parity distinction
and the order-one remainder are essential.

## A valued order-drop obstruction with unrestricted partner order

Let O be a discrete valuation ring in a one-variable algebraic
differential function field F over characteristic-zero constants.
Assume delta(O) is contained in O, so reduction defines differential
operators. Suppose P,Q belong to O[D;delta], satisfy (1), and:
P has generic order exactly two, its leading coefficient a_2 is a
valuation unit, and the reduction of Q has odd order at least three.
Then no such pair exists. The reduction of kappa may be zero.

Here is the argument that covers every generic partner order at once.
If Q has even order 2j, comparison at order 2j+1 gives
\[
          c=b_{2j}/a_2^j\in\ker(\delta)\cap O.
\]
The unit hypothesis ensures c is regular. Subtract c P^j.
The commutator and coefficient regularity persist and the order
strictly drops. Generic odd orders at least three are excluded by
the theorem. The process must end at order one, so
\[
 Q=A(P)+Q_1,\qquad A\in(\ker(\delta)\cap O)[t],
                     Q_1\in O[D;\delta],\quad\operatorname{ord}Q_1=1.
                                                                    \tag{15}
\]
After reduction, a nonconstant surviving A(P) has even order, because
the reduction of P still has order two and a nonzero leading
coefficient. The remaining summand has order at most one and cannot
cancel that leading even term. If no nonconstant A survives, the
reduction has order at most one. Either outcome contradicts the
assumed odd order at least three.

No estimate for ord(Q) is used. No square root of a unit in O is
needed: the normalization cover belongs only to the characteristic-zero
theorem; the regular even-power subtraction occurs in the original O.

## Scope and verification

This is the algebraic-curve finite-mate theorem used by the
[all-degree scalar classification](UNIVERSAL_P_POWER_OBSTRUCTION.md).
The hierarchy, adjoint, Puiseux and valuation arguments above supply the
universal proof. Complementary finite Ore controls are retained with that
classification; no finite prefix substitutes for the all-order argument.
No formal verification, independent whole-proof replay or external review is claimed.
