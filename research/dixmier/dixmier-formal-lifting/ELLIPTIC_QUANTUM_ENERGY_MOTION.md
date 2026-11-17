# Elliptic cohomology forces the quantum energy to move

14 November 2026. **QEM1.** A rational formal quantization of the
normalized Long three-isogeny cannot keep its cubic quantum energy
commuting with the natural target cubic, even through the first
three commutator orders. This holds with arbitrary rational
coefficient poles and no common degree or differential-order bound.

Let the rational Moyal field be k(x,y)[[h]], with [y,x]=h, and put
\[
 b=y^2-x^3,\quad
 u=4y^2/x^2-3x,\quad v=8y^3/x^3-9y.                           \tag{1}
\]
Then {v,u}=27 and v^2-u^3=-27b. Suppose U=u+O(h), V=v+O(h)
satisfy [V,U]=27h, and let
\[
                         E=V\star V-U\star U\star U .
\]
The theorem is
\[
                  2\le v_h([E,b]_\star)\le3.                 \tag{2}
\]
For a Weyl-even pair, the valuation is exactly three.
In particular E is never a formal rational function of b in
quantum functional calculus. Exact energy preservation and every
univariate rational energy reparametrization are excluded.

The obstruction is the nonzero elliptic de Rham class
\[
             -{5\over72b^2}\left[{x\,dx\over y}\right].        \tag{3}
\]
It annihilates every rational correction of the required type;
it is not a residue on a fixed support cutoff.

A separate positive construction gives a flat two-parameter Ore
deformation of QIM1's curve algebra, with a central reduction equal
to a localized Weyl algebra. Thus the objects can be connected
algebraically. The isogeny map must change its energy transversely;
the rigid exact-energy map does not extend.

This is a restriction on a distinguished energy and its isogeny,
not a nonlifting theorem for Long or a DC2 solution.

## Rational Moyal conventions and the elliptic forms

Use
\[
 f\star g=\sum_{n\ge0}{h^n\over2^n n!}\Lambda^n(f,g),\qquad
 \Lambda^n(f,g)=\sum_{i=0}^n(-1)^i{n\choose i}
 f_{y^{n-i}x^i}g_{x^{n-i}y^i}.                               \tag{4}
\]
Each coefficient is rational. Any series with nonzero rational
leading term has a star inverse. We work in this completed rational
division algebra, not only with finite-order scalar operators.

The generic level E_b: y^2=x^3+b is a smooth projective elliptic
curve after adjoining its point at infinity. On it the classes
\[
               \omega_0=dx/y,\qquad \omega_1=x\,dx/y           \tag{5}
\]
are linearly independent modulo exact rational differentials.
The first is nonzero and holomorphic. The second has only a pole
of order two at infinity. An exact primitive of a nonzero
combination with omega_1 would have precisely one simple pole,
giving a degree-one map of an elliptic curve to the projective
line, which is impossible. If only omega_0 remains, a rational
primitive could have no poles and hence would be constant.

The Hamiltonian vector field of b is
D_b=2y*partial_x+3x^2*partial_y. Its rational constants are
exactly k(b): the generic fiber is geometrically integral and
one-dimensional, and a nonzero derivation on its function field
has only algebraic constants over the base.

This is the relative-cohomology mechanism familiar from
[Garay--van Straten, Section 2](https://arxiv.org/html/0802.1647)
and the related
[Sharygin--Talalaev obstruction framework](https://arxiv.org/html/1210.2840).
The calculation and the specific nonzero class here are supplied
below, not imported from those papers.

## Arbitrary first-order corrections can be normalized

Assume, toward (2), that [E,b] is divisible by h^4. Writing the
first two energy coefficients as e_1,e_2, equation (4) gives
{e_1,b}={e_2,b}=0. Therefore
\[
                  E=-27b+h g_1(b)+h^2 g_2(b)+O(h^3),          \tag{6}
\]
with rational g_i. No assumption about higher coefficients is
needed.

Write U=u+hU_1+..., V=v+hV_1+.... Since the rational coordinate
change (x,y)->(u,v) has constant Jacobian 27, the coefficient pair
defines a unique rational vector field X_1 with X_1(u)=U_1,
X_1(v)=V_1. The commutation and energy equations at first order say
\[
             \operatorname{div}X_1=0,\qquad
                         X_1(b)=-g_1(b)/27.                  \tag{7}
\]
Set H=2x*partial_x+3y*partial_y, so H(b)=6b and div(H)=5.
Every field with the second property in (7) is
X_1=-g_1(b)H/(162b)+fD_b. Its divergence equation restricts to
\[
 df={6g_1'(b)-g_1(b)/b\over324}\,\omega_0
                     \quad\hbox{on the generic fiber}.
\]
Nonexactness of omega_0 forces 6b g_1'=g_1. A nonzero rational
solution would have order 1/6 at b=0, impossible for a rational
function. Thus g_1=0, and X_1=f(b)D_b.

There is a formal rational star automorphism fixing b exactly
whose first vector field is f(b)D_b. Here is a construction,
including the possible logarithmic primitives.

After a constant algebraic extension, rational integration writes
f(b)=F'(b)+sum_beta c_beta/(b-beta), with F rational.
The inner derivation ad_star(F(b))/h has the first vector field
F'(b)D_b; rational F(b) here means quantum functional calculus.
For A=b-beta, the automorphism Ad_star(A) is identity modulo h,
fixes b, and
\[
 {1\over h}\log(\operatorname{Ad}_\star A)
     ={1\over h}\sum_{n\ge1}{(-1)^{n+1}\over n}
                   (\operatorname{Ad}_\star A-1)^n            \tag{8}
\]
is a continuous rational-coefficient derivation with leading field
D_b/(b-beta). The logarithm of an identity-mod-h automorphism is a
derivation: modulo h^N its binomial powers are polynomial in an
auxiliary exponent, are multiplicative for every integer exponent,
hence for every exponent; differentiate at zero.
The series in (8) is divisible by h and converges coefficientwise.

Take the indicated linear combination of these derivations and
exponentiate minus h times it. This fixes b and removes U_1,V_1.
It preserves (6) through order h^2, so we may assume U_1=V_1=0.

This normalization makes no claim of finite support or bounded
inverse. Allowing arbitrary formal rational gauges strengthens
the ensuing impossibility statement.

## The second-order obstruction in every rational degree

Define the fixed Moyal defects of the classical pair:
\[
 \begin{split}
 A&={1\over24}\Lambda^3(v,u)={96y^2\over x^8},\\
 Q&={\Lambda^2(v,v)-u\Lambda^2(u,u)-\Lambda^2(u^2,u)\over8}\\
  &=-{6(3x^6+8x^3y^2-8y^4)\over x^8}.
 \end{split}                                                   \tag{9}
\]
Thus A is the h^3 commutator defect and Q the h^2 cubic-energy
defect. For U=u+h^2U_2+..., V=v+h^2V_2+..., let X(u)=U_2,
X(v)=V_2. The two equations become
\[
       \operatorname{div}X=-A/27,\qquad
                            X(b)=(Q-g_2(b))/27.               \tag{10}
\]
Put X_0=(Q-g_2(b))H/(162b). Then X-X_0=fD_b, and solvability
requires a rational solution of
\[
 D_b f=F_0+{6g_2'(b)-g_2(b)/b\over162},\qquad
 F_0=-A/27-\operatorname{div}\bigl(QH/(162b)\bigr).
\]
On y^2=x^3+b this fixed right side is
\[
             F_0=-{56b^2+56bx^3+15x^6\over27b x^8}.
\]
There is the exact rational-differential identity
\[
 {F_0\,dx\over2y}
   =d\left({4y\over27x^7}+{y\over18b x^4}
                         +{5y\over36b^2x}\right)
                         -{5\over72b^2}\,\omega_1.             \tag{11}
\]
It follows either by direct differentiation or successively from
d(y/x^n)=((3/2-n)x^(2-n)-nb*x^(-n-1))*dx/y for n=7,4,1.

Consequently the necessary exactness condition is
\[
 -{5\over72b^2}[\omega_1]
          +{6g_2'(b)-g_2(b)/b\over324}[\omega_0]=0.             \tag{12}
\]
The two classes are independent and the first coefficient is
nonzero. No rational f exists, whatever its poles or degrees.
This contradicts (6), proves [E,b] is not divisible by h^4, and
establishes (2). The lower bound two follows directly from
E=-27b+O(h). For a Weyl-even pair, E is even and its commutator
with b is odd, so its valuation must be three.

The class in (11) also cannot acquire a rational primitive
after a finite algebraic cover of the elliptic fiber: trace a
hypothetical primitive back and divide by the cover degree.
This statement concerns the normalized second-order equation.
The full arbitrary-first-order theorem above is stated for
rational coefficients; no blanket theorem on all algebraic
coefficient covers is claimed.

## Formal lifting survives by moving the energy

There is no contradiction with [QFL1](PROOF.md). A rational
Laurent-position, polynomial-momentum formal pair exists, but
its energy must move. For example take
\[
 f=-{32y^3\over27x^8},\quad
 U_2=f u_y=-{256y^4\over27x^{10}},\quad
 V_2=f v_y=-{256y^5\over9x^{11}}+{32y^3\over3x^8}.              \tag{13}
\]
The vector field f*partial_y has divergence -A/27, so these
coefficients cancel the h^3 commutator defect. Its energy
coefficient is
\[
        G_2=Q-54yf=-18/x^2-48y^2/x^5+112y^4/x^8,
\]
and {G_2,b} is nonzero. The resulting energy commutator already
has the forced h^3 term.

The same unrestricted polynomial homotopy as QFL1 extends
this even jet to all orders in
k[x,x^{-1},y][[h]]. Integration in y remains polynomial in y
and introduces no new position divisors. No common bound on
the Laurent pole orders or momentum degrees is supplied.

## A flat Ore bridge between the two kinds of algebra

The curve model and the Weyl algebra can themselves be linked.
Let S_epsilon be the localized Weyl algebra with [y,x]=epsilon.
For H=y^2-x^3 define the polynomial derivation
delta_epsilon=-ad(H)/(2epsilon). The division notation means
coefficientwise division: every commutator with H is divisible
by epsilon. In particular delta(x)=-y, delta(y)=-3x^2/2.

Form the Ore algebra
\[
 {\cal H}_{\epsilon,\tau}=S_\epsilon[z;\tau\delta_\epsilon],
 \qquad C=2\epsilon z+\tau H\ \hbox{central}.                  \tag{14}
\]
It has a free PBW basis x^i y^j z^k, i in Z and j,k>=0,
over k[epsilon,tau]. It is flat over those parameters.
At epsilon=0 this is exactly the total relative curve Ore
algebra in QIM1. For epsilon invertible it is
S_epsilon[C], by solving for z. On the diagonal epsilon=tau=h,
the element chi=2z+H is central and
\[
                       {\cal H}_{h,h}/(\chi)
                              \simeq S_h.                    \tag{15}
\]
This is a finite algebraic central reduction with noncentral
cubic energy. The same inner-derivation construction works
for any polynomial Hamiltonian; the cubic makes the generator
relations particularly simple. Compare the general
[Ore-extension semiclassical framework](https://arxiv.org/html/1212.4063).

The flat bridge is a construction of algebras, not of an
isogeny homomorphism or of a deformation of QIM1's free module.
A deformation retaining the exact identification H'=-27H,
or even H' a univariate formal rational function of H after
reduction, is excluded by (2). A more flexible map may
deform both the invariant derivation and the energy. It
remains open. Preservation of a central character alone
must not be silently identified with exact energy preservation.

## Application and canonical transport

In QEI1's exact-R Long chart let r=R and use the shifted generator
p=v_pos*a with [p,v_pos]=h. Put
x=8r*v_pos, y=8r*(3r*p-2*v_pos). Then [y,x]=192r^3 h.
The classical source energy coordinates are
u=16-12rS and v=432r^2T+72rS-64. They are the pair (1)
and have commutator [v,u]=27*(192r^3h). Rescale h by the
nonzero scalar 192r^3 over k=K(r).

Therefore every polynomial-coefficient exact-R formal Long
lift, even without a spatial bound, has a normalized quantum
energy E whose commutator with this natural cubic b has
valuation two or three; a lift even in this chart has valuation three.
The classical relation E=-27b cannot persist inside the
natural energy centralizer. Nonpositional deformations of R
are outside this chart and are not excluded.

Source formal canonical automorphisms fixing their cubic
energy leave E unchanged. Target identity-mod-h formal
canonical automorphisms fixing b transport [E,b]
isometrically. The obstruction therefore survives those
energy-preserving changes even when their degrees are
unbounded. Arbitrary changes moving the distinguished
energy are not covered.

The next bounded-realization question is whether the required
transverse energy motion can occur with one common spatial
bound and all four admitted images. A first nonzero
cohomology class alone does not answer that question.

## Verification boundary

[Exact controls](check_energy_motion.py) verify the isogeny
normalization, both Moyal defects, the universal vector-field
linearization, identity (11), the unconstrained formal escape
jet and the Ore central-reduction identities. The nonexactness,
rational logarithmic normalization, all-degree impossibility,
canonical transport and PBW flatness are written proofs.
No finite support census or increasing-order campaign is used.

No bounded full Long lift, all-representative nonlifting theorem,
DC2 solution, independent theorem replay, formal verification,
external review or novelty priority is claimed.
