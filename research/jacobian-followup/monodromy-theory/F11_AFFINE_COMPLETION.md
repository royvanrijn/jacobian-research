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

The next affine-completion task is now specific: realize the extension
in (6)--(7) on ordinary affine source and target spaces with full
constant Jacobian, or prove that extension impossible. An approach
retaining a,c as literal independent source coordinates and using one
polynomial remaining output cannot work: its Jacobian is the derivative
of that output in r, so constancy forces degree one in r. Source
coordinates must change essentially, or additional genuinely coupled
coordinates must be used. All degree, group and affine-space
conditions still have to be proved. Passing the necessary divisor
tests in (6) does not supply those proofs.

## Reproduction and assurance

Run

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify_affine_completion.py

The checker expands the two- and three-variable Jacobians over the exact
cyclotomic coefficient field, verifies the target-divisor identity,
the binary homogenization and boundary coprimality, tests the first
blowdown, and verifies (6)--(8), the target-factor irreducibility and
the tame branch-index ledger. Small counterexamples check why rational
parameters, characteristic zero and the Keller condition matter.
Universal exclusions have the written proofs above. Finite controls do
not constitute formal verification, a coefficient search, or an
independent whole-proof review.
