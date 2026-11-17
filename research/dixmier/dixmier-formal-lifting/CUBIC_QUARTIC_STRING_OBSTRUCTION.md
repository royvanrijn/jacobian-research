# Cubic and quartic scalar string operators over algebraic curves

**QSO34.** Let C be an algebraically closed characteristic-zero field,
let F/C be an algebraic function field of one variable, and let delta be
a nonzero C-derivation of F. There are no finite operators P,Q in
F[D;delta] with
\[
 \operatorname{ord}P=3,\qquad \operatorname{ord}Q=4,\qquad
                  [P,Q]=\kappa\in C^\times.              \tag{1}
\]
The same exclusion holds over arbitrary characteristic-zero constant
fields after extending constants and choosing a component, and after
every finite algebraic coefficient cover.

This extends the retained quadratic obstruction to the particular
order pair (3,4). It does not exclude order (3,5), arbitrary higher
partner orders for a cubic operator, or general finite Long lifts.
No novelty priority or independent external review is claimed.

## Complete normalization of an arbitrary order-(3,4) pair

Adjoin a cube root rho of the leading coefficient of P and put
E=rho*D. Its derivation epsilon=rho*delta is nonzero with the same
constant field. P is now monic in E. A coefficient shift of E removes
its quadratic term, giving
\[
                       P=E^3+uE+v.                       \tag{2}
\]
These are finite Ore changes over an algebraic function field; no
exponential gauge function or algebraic primitive has been assumed.

The leading commutator equation makes Q's leading coefficient a
nonzero constant. Divide Q by that constant, changing kappa but
keeping it nonzero. Its cubic coefficient is then constant and can
be removed by subtracting a constant multiple of P. Write
Q=E^4+b_2E^2+b_1E+b_0 and let primes mean epsilon.
The coefficients of E^4,E^3,E^2 in [P,Q] give, successively,
\[
 \begin{split}
 b_2&=4u/3+c,\\
 b_1&=2u'/3+4v/3+c_1,\\
 b_0&=2u''/9+2u^2/9+2v'/3+2cu/3+c_0,
 \end{split}                                             \tag{3}
\]
with c,c_1,c_0 in C. Remove c_0 from Q. Set
U=-2u/3 and initially V=2v/3-u'/3. Replacing
V by V+c_1/2 amounts to shifting P by the scalar 3c_1/4.
All these operations preserve the orders and the constant commutator.
The resulting pair is exactly
\[
 \begin{split}
 P&=E^3-\tfrac32UE-\tfrac34U'+\tfrac32V,\\
 Q&=E^4+(c-2U)E^2+2(V-U')E+
          V'+\tfrac12U^2-\tfrac56U''-cU .
 \end{split}                                             \tag{4}
\]
Thus every lower coefficient and constant flow has been retained or
removed by an explicit permitted operation.

Define the differential polynomials
\[
 \begin{split}
 I_1&=\tfrac12V''-\tfrac32UV+\tfrac32cV,\\
 I_2&=\tfrac1{12}U^{(4)}-\tfrac34UU''-\tfrac38(U')^2
       +\tfrac32V^2+\tfrac12U^3-\tfrac c4(3U^2-U'').
 \end{split}                                             \tag{5}
\]
Direct finite Ore multiplication gives the universal identity
\[
                  [P,Q]=-2I_1'E-I_1''-I_2'.              \tag{6}
\]
Consequently (1) forces I_1=-gamma in C and I_2'=-kappa.
In particular
\[
                 x=-I_2/\kappa\in F,\qquad \epsilon x=1. \tag{7}
\]
This constructs the required time primitive inside F. It is nonconstant,
so F is finite over C(x), and epsilon is the extension of d/dx.
The normalized equations are
\[
 \begin{cases}
 \tfrac12V''-\tfrac32UV+\tfrac32cV+\gamma=0,\\
 \tfrac1{12}U^{(4)}-\tfrac34UU''-\tfrac38(U')^2+
       \tfrac32V^2+\tfrac12U^3-\tfrac c4(3U^2-U'')+\kappa x=0.
 \end{cases}                                             \tag{8}
\]

For comparison, [Hayford, arXiv:2405.03260v2](https://arxiv.org/html/2405.03260v2#A2)
gives the normalized operators and equations in Appendix A.2,
(A.8), (A.12), Proposition A.1 and (A.13), with c=5t_5/3 and
kappa=1. The normalization above treats an arbitrary algebraic
coefficient pair and supplies x rather than assuming it.
The exclusion of algebraic solutions below is a written argument;
it is not imported from Hayford's isomonodromy results or from an
irreducibility conjecture.

## All finite-pole dominant balances

Fix a finite value x_0 and put t=x-x_0. At every place of an
algebraic coefficient field above t=0, U,V have Puiseux Laurent
expansions in t^{1/e}. This follows by completing the normalized
curve at that place: t=z^e times a unit, whose e-th root is a
formal unit in characteristic zero.

Suppose U or V has a pole. If V has exponent beta<0, its first
equation has terms of exponents beta-2, alpha+beta and beta,
where alpha is the leading exponent of U. The nonzero leading
coefficient beta(beta-1) of V'' forces alpha=-2: either inequality
would leave a unique most singular term. Thus a pole of V cannot
occur with U bounded.

If U has a pole and V is bounded, the second equation compares
alpha-4, 2alpha-2 and 3alpha. For alpha<-2 the cubic term is
uniquely most singular; for -2<alpha<0 the fourth derivative is.
Hence again alpha=-2. These arguments include an identically
zero V or U wherever the equation permits it.

Write U=u*t^{-2}+..., with u nonzero. The V^2 term in the second
equation rules out every exponent beta<-3. Write V=v*t^{-3}+...,
allowing v=0 when V has smaller pole order or is bounded.
At weight six the second equation gives
\[
                 \tfrac12u(u-2)(u-10)+\tfrac32v^2=0.      \tag{9}
\]
If v is nonzero, the first equation gives u=4. Therefore all
possible principal coefficient pairs are exactly
\[
                    (u,v)=(2,0),(10,0),(4,4),(4,-4).       \tag{10}
\]
In particular unequal pole orders and zero principal V have
not been discarded.

## The first fractional exponent cannot occur

Give U and V the pole weights two and three:
\[
 U=\sum_{r\in(1/e)\mathbb Z_{\ge0}}u_r t^{r-2},\qquad
 V=\sum_{r\in(1/e)\mathbb Z_{\ge0}}v_r t^{r-3}.
\]
Let r>0 be the first noninteger index for which either coefficient
is nonzero. Earlier coefficients have integer indices. Taking
the coefficients of t^{r-5},t^{r-6} in (8) gives
\[
 M_{u,v}(r)\binom{u_r}{v_r}=0,\quad
 M_{u,v}(r)=
 \begin{pmatrix}
 -3v/2&((r-3)(r-4)-3u)/2\\
 A_u(r)&3v
 \end{pmatrix},                                          \tag{11}
\]
where
\[
 A_u(r)=\frac{(r-2)(r-3)(r-4)(r-5)}{12}
 -\frac{3u}{4}((r-2)(r-3)+6)
 +\frac{3u}{2}(r-2)+\frac{3u^2}{2}.                       \tag{12}
\]
Lower-weight c terms and the polynomial forcing only involve
earlier indices or integer indices. A nonlinear term containing
the first fractional coefficient contributes at that same weight
only through the displayed principal coefficients.

Up to sign, the determinant in (11) factors as follows:

| Principal pair | Minus the determinant |
|---|---|
| (2,0) | (r-8)(r-6)(r-4)(r-3)(r-1)(r+1)/24 |
| (10,0) | (r-12)(r-9)(r-8)(r+1)(r+2)(r+5)/24 |
| (4,+/-4) | (r-9)(r-8)(r-4)(r-3)(r+1)(r+2)/24 |

Every root is an integer. The matrix is thus invertible at the
noninteger r, a contradiction. All finite-pole expansions of U,V
have integer exponents. This argument only excludes fractional
terms in an algebraic solution already assumed to exist; it does
not need to prove resonance compatibility or the full analytic
Painleve property.

At a finite point where U,V are bounded, the first equation makes
V'' bounded. If the first fractional exponent of U were below four,
U^{(4)} would be uniquely most singular in the second equation:
UU'' and (U')^2 have strictly larger exponents, and V^2 is bounded.
Thus U''' and V' are bounded as well. The regular first-order
system for (U,U',U'',U''',V,V') has unique formal Taylor solution
for those finite initial values. Equivalently its first fractional
coefficient would occur one exponent earlier on the differentiated
left side than on the polynomial right side, and must vanish.
Hence there is no finite branching at these points either.

## Unramified algebraic functions are rational here

Let J=C(x,U,V). The local result makes its normalized projective
curve unramified over every finite point of the x-line: all U,V
expansions lie in C((t)), so their generated completed field has
ramification index one. Let d=[J:C(x)], let g be the genus, and
let b>=1 be the number of points over infinity. Characteristic
zero makes all ramification tame. The
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
therefore gives
\[
                   2g-2=-2d+(d-b)=-d-b.                 \tag{13}
\]
Since g>=0, this forces d=b=1. Thus U,V are rational functions
of x. Only this elementary curve covering consequence is used,
not an irreducibility theorem for the string equation.

## Infinity excludes every rational solution

Suppose U has positive integer degree A at infinity. If V is
nonzero with integer degree B, the first equation has its UV
term strictly above the degrees of V and V''. It follows that
A+B<=0, since otherwise that term is uniquely leading.
Thus B<=-A; if V=0 the same conclusion about its negligible
contribution is immediate.

The second equation now has its U^3 term of degree 3A>=3.
Every other degree is smaller: those of UU'' and (U')^2 are
at most 2A-2, the cU^2 degree is 2A, V^2 has negative degree,
and the forcing has degree one. This is impossible.

If U is bounded at infinity, all its differential-polynomial
terms are bounded. If V has positive integer degree B, V^2 has
degree 2B>=2 and is uniquely leading. If V is bounded, the
nonzero term kappa*x is uniquely leading. Both are impossible.
The identically zero cases are included. This proves (1).

## A consequence for exact-R Long repairs

In QPH1's Laurent chart, h-regular finite exact-R images S,T
have classical a-orders two and three, with leading coefficients
-3r and r/2 over K(r). QSO2 already excludes generic order two
for S. The present theorem excludes generic orders
\[
                         (\operatorname{ord}_a S,
                           \operatorname{ord}_a T)=(3,4). \tag{14}
\]
It also implies a useful minimal-order consequence together
with QSO2: if generic S has order three, a finite T with the
prescribed Long reduction must have generic order at least five.

Indeed generic order of T is at least three by h regularity.
If it were three, leading commutator comparison would make
a=T_3/S_3 a differential constant. T-aS has order at most two
and the same nonzero constant commutator with S. Order two is
excluded by QSO2, applied with that residual as its quadratic
operator and S as its odd-order mate. Order zero cannot have
constant nonzero commutator with a cubic operator, by its leading
commutator coefficient. If the residual has order one, comparison
of its a^2 coefficient gives T_2=a*S_2. The right side has a
negative h valuation: T_3 has valuation zero, S_3 has positive
valuation, and S_2 has valuation zero. The left side is integral,
a contradiction. Order four is excluded by (1), leaving at least
five. These are generic operator orders, not spatial-degree bounds.

No arbitrary higher-order cubic partner, deformed-R scalar lift,
matrix/projective carrier, finite normalization fixing R, bounded
Long lift, or DC1/DC2 solution is constructed or excluded here.

## Evidence

The standalone checker verifies (3), the full universal identity
(6), all principal balances and the exact resonance matrix/factors,
using independent finite normal-order Ore multiplication and
symbolic differentiation. A deliberately corrupted coefficient
is rejected. It validates the pinned primary PDF and writes a
compact fresh receipt. The all-place Puiseux exclusion, curve
argument, infinity comparison, arbitrary-pair normalization and
Long order consequence are written proofs, not finite scans.

Run from the repository root:

```sh
.venv/bin/python research/dixmier/dixmier-formal-lifting/check_cubic_quartic_string.py
```
