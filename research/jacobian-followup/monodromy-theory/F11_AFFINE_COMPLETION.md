# The F11 affine-completion problem: excluded target models and an admissible replacement

This continues [the symmetry theory](THEORY.md) and the [F11 audit](F11.md).
The previous turn is progress: its proof files, four ledger entries and
exact replay are present. The new work distinguishes a forbidden affine
model of a cover from a forbidden monodromy action. The full action
spectrum and absolute PSL2(11) realization remain open.

All fields in the geometric statements are algebraically closed of
characteristic zero unless stated otherwise. N,Q,D=N-Q,A,B,C,E,R and
lambda,k,c_0 have the meanings in F11.md.

## 1. A rational curve cover cannot retain a polynomial base parameter

**KMFI5, pole form.** Let R_0=k[x_1,...,x_m], let r=a/b in Frac(R_0)
be nonconstant with gcd(a,b)=1, and let f in k(T) be nonconstant.
If f(r) belongs to R_0, then f has just one pole on P^1. After a Mobius
change of r, both the changed r and the changed outer function are
polynomials. This statement does not require a Keller map.

Here is a direct proof. Express f=P_h/Q_h in coprime binary homogeneous
forms of the same degree. Since k is algebraically closed, their factors
are linear forms. The polynomials P_h(a,b),Q_h(a,b) are coprime: a
common prime factor would divide independent linear combinations of
a,b, and hence both a,b. Polynomiality of their quotient therefore
makes Q_h(a,b) a nonzero constant. If f has two distinct poles, Q_h
has two independent linear factors. Both corresponding linear
combinations of a,b must be constant, forcing a,b constant. This is
a contradiction. With only one pole, choose that linear factor as
the denominator in the new Mobius coordinate. Its evaluated
denominator is constant, proving the final assertion.

F11's beta has six distinct poles. Consequently there is no nonconstant
r in the rational function field of any affine space with beta(r)
polynomial. This excludes every polynomial completion that keeps the
old base parameter t regular, in any number of added variables and
with arbitrary nonlinear feedback. It does not exclude rational
target parameters with poles.

A related classical coordinate statement is worth recording. If a
nonconstant polynomial f has gcd of its partial derivatives equal to
one, k(f) is algebraically closed in k(x_1,...,x_m). Indeed, the
generative-polynomial criterion reduces failure to f=g(h), with
h polynomial and deg g>1; then g'(h) is a common nonconstant factor
of all partial derivatives. The criterion and field formulation are
[Arzhantsev--Petravchuk, Theorem 1 and Lemma 3](https://arxiv.org/html/math/0608157v2).
A primary-source snapshot is retained in
[sources/closed-polynomials-manifest.json](sources/closed-polynomials-manifest.json).
This use of the classical criterion is an attributed input, not a
new closed-polynomial theorem.

For a Keller map F, every nonzero affine-linear combination of its
coordinates has a unimodular gradient. More generally, if a target
polynomial t has no critical divisor, neither does t(F): the invertible
Jacobian matrix identifies the two gradient ideals, and a source
divisor cannot map into a target subset of codimension at least two.
Hence k(t(F)) is algebraically closed in the source rational field.
In particular no nontrivial finite curve cover can be pulled back
along a polynomial target coordinate and retained as an algebraic
subfield of the generic source. This obstruction is about an embedded
cover field, not only its permutation group.

## 2. Explicit polynomial cone and target-divisor Jacobian lift

**KMFI6.** There is a polynomial map on ordinary affine two-space,

\[
 \Phi(x,w)=(U,V)=(wN(x),wQ(x)),
 \qquad
 \det D\Phi=wJ(x),\quad
 J=N'Q-NQ'=\lambda B^4C^4R .
 \tag{1}
\]

It has generic degree twelve and geometric natural PSL2(11) monodromy.
This follows from U/V=beta(x) and w=V/Q(x); the target field is
k(beta,V) and the source is k(x,V). Adding the independent variable V
preserves the original extension and Galois closure. Over the displayed
arithmetic field K=Q(zeta_10), the arithmetic group is PGL2(11), with
the same quadratic constant extension specified in F11.md.

Define the polynomial

\[
 H=\frac{NDQ}{J}=\frac{c_0k}{\lambda}ABCER .
 \tag{2}
\]

Then a polynomial map on ordinary affine three-space is

\[
 \Psi(x,w,z)=(wN,wQ,z\,w^2H),\qquad
 \boxed{\det D\Psi=U(U-V)V.}
 \tag{3}
\]

Its generic degree and group are still twelve and natural PSL2(11):
the last input is recovered by z=Z/(w^2H). Formula (3) follows by
expanding the block determinant, which is w^3 JH=w^3NDQ.
The first two coordinate degrees are thirteen; the third has degree
sixteen. The nonconstant determinant is explicit, so neither map is
Keller.

For U V(U-V) nonzero the full inverse scheme is the degree-twelve
squarefree root algebra of V N(X)-U Q(X), followed by the displayed
unique reconstructions of w and z. Its leading coefficient is a
nonzero constant multiple of V-U. All source factors in (2) are units
in that root algebra: a zero would force U=0, V=0 or U=V. Thus no
generic inverse sheet was lost in passing to these polynomial models.

This construction reaches polynomial affine-space algebraization of
the exact cover, together with a Jacobian equal to a target divisor.
It does not make the determinant constant.

## 3. The cone field is excluded with this affine target

The cone's divisor V=0 has t=U/V of pole order one. Its inertia is the
F11 involution of cycle type 2^6, with no fixed sheet. KMFI1 therefore
excludes *every* Keller realization of this embedded field extension
over this same affine target, including arbitrary source recharts and
identity stabilizations. The auxiliary variable in (3) does not alter
this inertia.

Polynomial target automorphisms cannot repair this: they send V=0
to another affine prime divisor with the same pole order and inertia.
This is stronger than the assertion that the displayed polynomial map
has a nonconstant Jacobian. A genuinely different affine target model
is necessary. The target-divisor equality in (3) is not itself evidence
that such a target model exists.

One first blowdown is completely explicit:

\[
 (U,V)\longmapsto (S,T)=(UV,V).
\]

It is birational, with inverse U=S/T,V=T, and t=S/T^2. The polynomial
source map becomes (w^2NQ,wQ). Its new generic divisor T=0 has pole
order two, killing the order-two inertia. The old derangement divisor
is contracted to (S,T)=(0,0). The first obstruction is removed.

But S is squarefree, and Section 4 excludes every Keller completion
of this new field over A^2_(S,T) as well. This shows why relocating
the derangement branch alone does not solve affine completion.

## 4. Three homogeneous descent equations

**KMFI7, necessary target form.** Suppose a Keller map on ordinary
affine space has, inside its generic source field, an F11 root r with

\[
 \beta(r)=\frac{P(y)}{Q_t(y)},\qquad
 P,Q_t\in k[y_1,\ldots,y_m],\quad \gcd(P,Q_t)=1 .
\]

Write r=a/b with a,b coprime source polynomials. Homogenize the
univariate N,Q to binary degree twelve, and A,B,C,E,R to their own
degrees. Since deg D=11,

\[
 N_h=c_0 A_hB_h^5,\qquad Q_h=kR_h^2,\qquad
 N_h-Q_h=b C_h^5E_h .
\]

The two reduced fractions agree, so a single constant mu in k* gives

\[
 \begin{aligned}
 P(F)&=\mu c_0 A_h(a,b)B_h(a,b)^5,\\
 Q_t(F)&=\mu k R_h(a,b)^2,\\
 (P-Q_t)(F)&=\mu b\,C_h(a,b)^5E_h(a,b).
 \end{aligned}
 \tag{4}
\]

To justify that the multiplier is constant, both fractions are reduced
in the source polynomial UFD. The binary forms N_h,Q_h stay coprime
after substituting a,b, as in Section 1. The target pullbacks P(F),
Q_t(F) are also coprime: any common source divisor would, by
quasi-finiteness, dominate a common target divisor of P,Q_t.

Equations (4), together with etaleness, have three consequences.

* Q_t is a nonconstant square up to a scalar. Its pullback is a
  square; primitive divisor pullback from KMFI1 forces all target
  valuations to be even. Constancy would contradict Section 1.
* P has a prime-divisor multiplicity which is a positive multiple of
  five. The two distinct roots of B give independent binary linear
  forms. At most one evaluates to a constant, so B_h(a,b) is
  nonconstant. A source prime dividing it has multiplicity 5e in
  P(F), since A_h,B_h are coprime. Its image is a target divisor,
  with the same valuation because F is etale.
* P-Q_t also has a prime multiplicity divisible by five, using the
  two distinct roots of C and the third equation in (4).

In particular, neither numerator can be squarefree. The local
fixed-sheet inertia test alone did not impose this requirement.

There is a uniform lower-degree consequence:

\[
 \boxed{\max(\deg P,\deg Q_t)\ge6.}
 \tag{5}
\]

If both degrees were at most five, then Q_t, being a nonconstant
square, would have degree at most four. The other conclusions force
P=cL^5 and P-Q_t=dM^5, where L,M are nonconstant affine-linear
polynomials. Thus Q_t=cL^5-dM^5. If L,M are proportional this is a
nonzero scalar times a fifth power, not a square. Otherwise it is a
product of five distinct nonproportional affine-linear factors; at
most one is constant. It again cannot be a square. This proves (5).
The bound is necessary, not asserted sharp.

The hypotheses are important. This is an obstruction to rational
reparameterizations of the *specified F11 subcover*. It neither
bounds the degree of arbitrary Keller maps nor excludes the abstract
natural PSL2(11) action.

## 5. A rational target model that passes these new tests

The preceding exclusions do not exhaust birational target changes.
There is a concrete replacement with rational source field and only
admissible affine-divisor inertia. On A^3_(a,c,z) set

\[
 \begin{aligned}
 P&=a^5\bigl(a^5+c^5(z+2)\bigr),\\
 Q_t&=(a^5+c^5)^2,\\
 P-Q_t&=c^5(a^5z-c^5),\qquad t=P/Q_t .
 \end{aligned}
 \tag{6}
\]

These fractions are reduced. Define B_t=a^5+c^5(z+2) and
D_t=a^5z-c^5. They are distinct irreducible polynomials, linear in z.
Their coefficients have gcd one, proving irreducibility. Also

\[
 k(a,c,z)=k(t,a,c),\qquad
 z=\frac{(a^5+c^5)^2t-a^5(a^5+2c^5)}{a^5c^5}.
 \tag{7}
\]

Thus the base is still purely transcendental over k(t). Pulling back
F11 preserves the full geometric PSL2(11) group. Its generic source
field is k(r,a,c), which is rational. The model is an actual specified
field extension, not a permutation passport without a cover.

Divisorial branch analysis is complete. At a=0, t has order five;
at c=0, t-1 has order five. Their original order-five inertia is killed.
Each component of a^5+c^5=0 is a pole of order two, killing the original
involution inertia. Only B_t=0 and D_t=0 remain branch divisors,
each with inertia 5^2 1^2 and fixed sheets. This uses the characteristic-zero
tame formula e_new=e/gcd(e,v_D(t-b)) after strict henselization.
Away from t=0,1,infinity, finite etale base change creates no branch.

The rational map realizing the extension can be written on the source
(r,a,c) as the two unchanged outputs a,c and the last output in (7)
with t=beta(r). Its rational Jacobian, with input order (a,c,r), is

\[
 \frac{(a^5+c^5)^2}{a^5c^5}\,\beta'(r).
 \tag{8}
\]

It has poles and zeros. It is not a polynomial Keller map. The
rational parameter in (6) has numerator degree eleven and denominator
degree ten; it passes the square-denominator and fifth-multiplicity
gates. No degree-minimality assertion is made.

This model passes the local tests but is excluded by the stronger
field obstruction below. It is retained as an explicit control:
rational source and target fields, exact prescribed group, square
denominator and admissible affine-divisor inertia are still insufficient.

## 6. A global obstruction to stabilization of a curve cover

**KMFI8, general form.** Let k(r)/k(t) be a finite rational curve
extension of degree n. Suppose one fiber of the corresponding map
P^1_r -> P^1_t has ramification indices e_1,...,e_s, and put

\[
 e=\gcd(e_1,\ldots,e_s)>1,\qquad d=n/e .
\]

If the extension, after adjoining the same independent variables on
both sides, admits any polynomial Keller model on ordinary affine
space, then necessarily

\[
 \gcd(e_i,d)\ \text{divides}\ e_i/e
 \quad\text{for every }i.
 \tag{9}
\]

Equivalently, every such fiber must satisfy

\[
 \boxed{\gcd(e,n/e)=1.}
 \tag{10}
\]

This includes uniform fibers e^s. To see the equivalence, put
f_i=e_i/e, so gcd(f_1,...,f_s)=1. If gcd(e,d)=1, then
gcd(e_i,d)=gcd(f_i,d) divides f_i. If a prime divides both e and d,
some f_i is not divisible by it, contradicting (9).

This is an obstruction to the embedded field pair, preserved under
arbitrary birational identifications of source and target and further
identity stabilization. It does not assume that the offending fiber
stays on an affine target divisor.

To prove it, move the chosen target value to infinity, and choose a
source coordinate r whose infinity is outside the selected fiber.
Write the rational cover in reduced degree-n homogeneous form, with

\[
 t=N(r)/Q(r),\quad Q(r)=cR(r)^e,\quad
 \deg R=d,\quad
 \operatorname{ord}_{r_i}R=e_i/e .
\]

The original source infinity can be chosen outside finitely many points,
so this coordinate choice is always possible over k. Suppose a Keller
model exists for

\[
 K=k(t,u_1,\ldots,u_j)\subset L=k(r,u_1,\ldots,u_j).
\]

In its actual polynomial coordinates, write t=P(y)/Q_t(y) reduced and
r=a(x)/b(x) reduced. Exactly the reduced-fraction and divisor-pullback
argument in Section 4 gives

\[
 Q_t(F)=\mu R_h(a,b)^e.
\]

KMFI1's power-injectivity implies Q_t is an e-th power up to a scalar.
Taking e-th roots over algebraically closed k, some S in K* satisfies

\[
 S=\nu R_h(a,b)=\nu b^dR(r),\qquad
 \boxed{S/R(r)\in L^{*d}.}
 \tag{11}
\]

Now use the valuation v_i of L at r=r_i, keeping every u_l
independent and of valuation zero. Its restriction to K has value
group e_i Z: a rational function of t and the independent u_l has
valuation a multiple of e_i, since t has pole order e_i.
Applying v_i to (11) therefore yields integers m,h with

\[
 e_i m-dh=e_i/e .
\]

This congruence forces (9). The valuation may lie on the boundary of
the putative target affine model; it remains a valuation of the same
field pair and cannot be discarded by a birational rechart.
This proves the obstruction.

For F11, the involution fiber is 2^6, so n=12,e=2,d=6. Equation (11)
would give S/R(r) a sixth power, while v_i(S) is even and v_i(R)=1.
An odd integer cannot be divisible by six. Thus **no stabilization of
the F11 curve-cover field pair has an absolute Keller realization in
any birational affine source and target coordinates**.

In particular, the locally admissible extension (6)--(7) is impossible
as a Keller extension, even after further identity stabilization.
Its source field is k(r,a,c) and target field k(t,a,c), exactly the
pair used in the argument. Locally killing the pole inertia on the
affine target did not change the obstructing valuation at the old
curve boundary. This is the additional obstruction missing from the
local inertia sieve.

## 7. Every natural degree-twelve curve input is excluded

The obstruction is not specific to the displayed F11 formula.
Every nonidentity element in the natural PSL2(11) action has index
6,8 or 10: order two gives 6, orders three and five give 8,
and orders six and eleven give 10. A connected genus-zero cover of
degree twelve satisfies

\[
 \sum_b\operatorname{ind}(\sigma_b)=22 .
\]

There can be at most three branch values, since each index is at
least six. The only index multisets are (6,6,10) and (6,8,8), so
there is always an involution branch of type 2^6. Equation (10)
excludes stabilization of every such cover. In fact the first
multiset has two involutions generating a dihedral group, and the
(2,3,3) and (2,3,5) triangle quotients have orders at most 12 and 60.
Thus full PSL2(11) leaves the familiar (2,5,5) passport. The checker
independently exhausts the finite matrix-group order-three pair
cases; passport agreement alone is never used to identify a group.

A positive-genus curve input cannot help by ordinary stabilization:
a rational function field dominating that curve would give a
dominant rational map from affine space to the curve. Restrict to
a general line where the map is nonconstant, extend to P^1 and apply
Riemann--Hurwitz; the curve must have genus zero.

Consequently **no curve cover at all in the natural degree-twelve
PSL2(11) action can be converted into an absolute Keller map merely
by independent affine stabilization and birational source/target
recharts**. This closes that whole curve-first route, not just one
ansatz or one chosen target boundary.

The action-level question remains open. A realization could arise
from a higher-dimensional field extension not isomorphic to any
such stabilized curve pair. This proof does not rule that out.
The exceptional degree-eleven actions also survive this particular
test: their totally ramified eleven-cycle fiber has d=1, and their
other fibers contain unramified sheets, so no common-index
congruence above fails.

The next action is therefore to work with a genuinely different
higher-dimensional natural-action cover, or to return to a
degree-eleven action and its unresolved polynomial-completion
equations. Repeating stabilization or birational completion of the
natural F11 field cannot solve the absolute problem.

## Reproduction and assurance

Run

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify_affine_completion.py

The checker expands the two- and three-variable Jacobians over the exact
cyclotomic coefficient field, verifies the target-divisor identity,
the binary homogenization and boundary coprimality, tests the first
blowdown, and verifies (6)--(8), the target-factor irreducibility and
the tame branch-index ledger. It checks the global valuation congruence,
the complete genus-zero index budget and the finite-group exclusions.
Small counterexamples check why rational
parameters, characteristic zero and the Keller condition matter.
Universal exclusions have the written proofs above. Finite controls do
not constitute formal verification, a coefficient search, or an
independent whole-proof review.
