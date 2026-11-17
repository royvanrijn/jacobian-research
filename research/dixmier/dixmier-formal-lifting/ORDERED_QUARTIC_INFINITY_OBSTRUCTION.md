# Every scalar ordering of the quartic square face has no finite mate

Let k have characteristic zero, eta,g,s in k^times and mu,u,w in k.
Work in the finite scalar differential algebra
\[
 \mathcal D=k(V)[A;\eta\delta],\qquad
                         \delta=V^{-1}\partial_V.
\]
Put H=V^2 A and
\[
 \begin{split}
 L&=s\eta^2V^{-7}\bigl((H/\eta)^2-u(H/\eta)+w\bigr)\\
  &=sV^{-3}A^2+(2-u)s\eta V^{-5}A+w s\eta^2V^{-7},\\
 P&=L^2+gA^2+\mu V.                                      \tag{1}
 \end{split}
\]

**QOI1.** No finite T over k(V), or any finite algebraic coefficient
cover, satisfies [P,T]=kappa for a nonzero scalar kappa. There is no
partner-order restriction and no exceptional value of u,w. The rational
coefficient centralizer is exactly
\[
                         \operatorname{Cent}_{\mathcal D}(P)=k[P].
                                                                  \tag{2}
\]
In particular the ordering choices u=5 and u^2-4w=16 do not supply
mates. The case mu=0 is also covered. The assumption g!=0 is necessary
for the polynomial-centralizer assertion: when g=mu=0, L commutes with
P=L^2 but is not a polynomial in P.

The universal obstruction is a growing positive-V coefficient of the
commutator, independent of the ordering parameters. It is distinct
from an obstruction obtained from a finite prefix of negative Laurent
residues. It also closes the globally admitted ordered-square family
described below.

## Every rational finite candidate has a finite even-flow normal form

The leading coefficient of P is p_4=s^2 V^-6. If T has order M and
[P,T] has order at most two, its coefficient at A^(M+3) gives
\[
 4p_4\delta(t_M)-M t_M\delta(p_4)=0,
                 \qquad t_M'/t_M=-3M/(2V).                  \tag{3}
\]
A nonzero rational function has integral valuation at V=0, so M
must be even. For M=2n its leading coefficient is a constant multiple
of s^n V^-3n. An order-zero remainder has constant coefficient.

Let R be the unique formal pseudodifferential square root R^2=P
whose leading term is sV^-3 A^2. Coefficient recursion divides only
by 2sV^-3. Thus each coefficient of R is a finite Laurent polynomial
in V with scalar coefficients rational in the nonzero parameter s.
This is a formal device; it is not asserted to be an actual operator.
For n>=0 set
\[
                           B_{2n}=(R^n)_+.
\]
The negative part of R^n has order at most minus one, so
[P,B_(2n)] has order at most two. Therefore (3) applies after each
leading subtraction. Finite repetition gives, for every actual T
whose commutator has order at most two,
\[
                            T=\sum_n c_{2n}B_{2n},
                                     \qquad c_{2n}\in k.   \tag{4}
\]
Every constant flow is retained. If n is even,
B_(2n)=P^(n/2) and its commutator vanishes. Only the odd n require
exclusion. No hypothesis on the rational position degree of T is used.

## A weighted symbol controls all positive fractional powers

Give V weight one and A weight three. Every nontrivial Ore contraction
lowers weight by five, since [A,V]=eta/V. The weights of P are
6,1,-4,... and its leading commutative weighted symbol is
\[
                       p_*=s^2V^{-6}\xi^4+g\xi^2.          \tag{5}
\]
All dependence on u,w and the term mu*V has lower weight. Recursion
for R and multiplication show that R^n has weights 3n-5l, l>=0.
Its coefficient at A^j therefore has V degree at most 3n-3j.
For B_(2n)=sum_j b_j(V)A^j, its top weighted symbol is
\[
             (p_*^{n/2})_+
       =\left(s^nV^{-3n}\xi^{2n}
          \sum_{a\ge0}{n/2\choose a}
                        (g s^{-2}V^6\xi^{-2})^a\right)_+ .
                                                                  \tag{6}
\]
This contains only even powers of xi. Accordingly
\[
 \deg_V b_j\le 3n-3j,\qquad
 \deg_V b_j\le 3n-3j-5\quad\hbox{if j is odd}.             \tag{7}
\]
The scalar coefficient has the sharper expansion
\[
 b_0=\alpha_n g^n s^{-n}V^{3n}
                     +\text{terms of V degree at most }3n-5,
              \qquad\alpha_n={n/2\choose n}.               \tag{8}
\]
For every odd n>=1 this binomial coefficient is nonzero.

These bounds can also be checked directly during the square-root
recursion: the highest V-degree term at a specified A order is solved
by division by sV^-3, and every derivative, lower ordering term or
mu insertion drops the weight by a positive multiple of five.
Thus the symbol computation in (6) controls all finite coefficients;
it does not exchange an unbounded spatial series with evaluation.

## The uncancellable coefficient for every odd n

For B=B_(2n), the coefficient of A in [gA^2,B] is
\[
                    2g\eta\delta(b_0)+g\eta^2\delta^2(b_1).
                                                                  \tag{9}
\]
Equation (8) gives its term of degree 3n-2 as
\[
              6n\eta\alpha_n g^{n+1}s^{-n}V^{3n-2}.        \tag{10}
\]
All remaining b_0 terms have degree at most 3n-7 after one delta,
and the odd bound for b_1 gives at most 3n-12 after two deltas.

The other parts of P cannot contribute at that degree. In [mu*V,B],
an A term uses b_j with j>=2 and exactly j-1 derivatives of V.
Its degree is at most
\[
                    (3n-3j)+(3-2j)=3n+3-5j\le3n-7.       \tag{11}
\]
Write L^2=sum_(i=0)^4 p_i(V) A^i. Each nonzero p_i is a scalar
multiple of V^(2i-14). A nontrivial commutator contribution at A
uses i+j-1 derivatives, with i+j-1>=1. Its degree is at most
\[
            (2i-14)+(3n-3j)-2(i+j-1)
                                      =3n-12-5j\le3n-12. \tag{12}
\]
Terms with no derivative cancel in the commutator.
It follows that, for every odd n,
\[
 [V^{3n-2}A]\,[P,B_{2n}]
                       =6n\eta{n/2\choose n}g^{n+1}s^{-n}\ne0,
                                                                  \tag{13}
\]
and no A coefficient in that commutator has a greater V degree.

In a finite sum (4), choose the largest odd n with c_(2n) nonzero.
Every smaller odd flow has a strictly smaller maximum V degree in
its A coefficient, while even n commute. Thus (13) cannot cancel.
The commutator cannot be scalar, including the scalar zero.
This proves (2) and excludes every nonzero scalar mate.

For example n=1 gives 3eta*g^2/s at V*A, and n=3 gives
-9eta*g^4/(8s^3) at V^7*A. The latter already excludes every
sixth-order mate for all u,w; a search of exceptional orderings is
unnecessary once this coefficient is retained.

## Arbitrary algebraic covers and fixed scalar matrix padding

For a finite algebraic extension E/k(V), derivations extend uniquely.
Normalized field trace commutes with delta. Applied to the finitely
many normal coefficients of T, it is a bimodule retraction
E[A;eta delta] to k(V)[A;eta delta]. It sends [P,T]=kappa to the
same scalar equation. A cover mate would therefore give the rational
mate just excluded. This proves the stated arbitrary finite-cover
conclusion without claiming that the odd-order parity argument itself
holds on every cover.

Normalized matrix trace gives the same exclusion for the fixed scalar
image P*I_d in every finite matrix size. Neither trace argument treats
a genuinely nonscalar first image, nor asserts polynomial centralizers
after arbitrary algebraic coefficient extension.

## A strict quartic infinity cone

The same proof has a useful more general form. Let
P=sum_(i=0)^4 p_i(V)A^i have arbitrary rational coefficients over k,
with the following asymptotics at infinity:
\[
 \begin{gathered}
 p_4=s^2V^{-6}+O(V^{-7}),\qquad
 p_2=g+O(V^{-1}),\qquad sg\ne0,\\
 \deg_V p_3\le-4,\qquad \deg_Vp_1\le2,\qquad
                                      \deg_Vp_0\le5 .       \tag{14}
 \end{gathered}
\]
A constant square-root extension for the leading coefficient s is
permitted. Then P has no finite scalar mate, on any finite algebraic
coefficient cover, and its rational centralizer is k[P].

Indeed the leading-coefficient equation (3), now with rational p_4,
still gives deg_V(t_M)=M*deg_V(p_4)/4=-3M/2, excluding odd M.
If p_4 is not a square up to a nonzero scalar in k(V), a rational
leading coefficient at M=4j+2 is impossible as well: its square
would make p_4 a scalar multiple of a rational square. Only orders
divisible by four remain, and finite subtraction of scalar powers
of P proves both assertions immediately.

If p_4 is a square up to a scalar, extend constants to choose a
rational square root. Formal R=P^(1/2) and all B_(2n) then have
rational coefficients. The weighted symbol is still (5). The
degree bounds become
\[
 \deg_V b_j\le3n-3j,\qquad
 \deg_V b_j\le3n-3j-1\quad\hbox{for odd j},                 \tag{15}
\]
with b_0's leading coefficient unchanged from (8). Integer infinity
degrees make this gap of one sufficient. Its remainder contributes
at most 3n-3 to (9), while the b_1 term contributes at most 3n-8.

Write P=s^2V^-6 A^4+gA^2+P_<, where every coefficient of P_<
has weighted degree at most five. An A coefficient in [P_<,B_(2n)]
with i+j-1>=1 has degree at most
\[
       (5-3i)+(3n-3j)-2(i+j-1)
                   =3n+7-5(i+j)\le3n-3 .                  \tag{16}
\]
The top quartic term contributes at most 3n-12 by the same count
as (12). Hence (13) remains the unique greatest-degree A term.
The even-flow normal form, finite descending maximum and trace
argument finish the proof. Descending the polynomial centralizer
back from a constant extension uses highest A coefficients, so
its scalar coefficients are already in k.

## The globally admitted ordered-square family

Use the exact global centralizer coordinates from
[the global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md):
\[
 \begin{split}
 R&=X(2-3XQ),\\
 A_0&=Q+3XQ^2-XB/3=a-v,\\
 C&=R A_0-2/3,\\
 S_b&=-3R A_0^2+(18XQ-8)A_0+9Q-3hRX
                                     =4a-3ra^2+2v.
 \end{split}                                               \tag{17}
\]
All these are actual finite global operators commuting with R.
For scalar ordering constants u,w, put
\[
 W_{u,w}=X^3C^2+(2-u)hR X^5C+w h^2R^2X^7,\qquad
 S_{u,w}=S_b+\epsilon h^2W_{u,w}^2,
                                  \qquad\epsilon\ne0.       \tag{18}
\]
These operators are globally polynomial. Every S_(u,w) has Long's
classical second image and generic differential order four.

Center z=a-2/(3r), put h=tau^5, v=tau^2 V, z=tau A, and subtract
the scalar 4/(3r) from S before dividing by tau^2. Since
C=r*tau*(A-tau V), direct finite multiplication gives
\[
 \tau^4 W_{u,w}
     =r^2\bigl[V^{-3}(A-\tau V)^2
                  +(2-u)V^{-5}(A-\tau V)+wV^{-7}\bigr].    \tag{19}
\]
This is integral with critical reduction
r^2[V^-3 A^2+(2-u)V^-5 A+wV^-7]. Consequently the critical
reduction of the rescaled S is (1), with
\[
                  g=-3r,\qquad \mu=2,\qquad \eta=1,
                                    \qquad s^2=\epsilon r^4.
                                                                  \tag{20}
\]
A constant square-root extension for s is harmless.

The strict infinity cone in fact applies to the full globally
admitted S directly, before any parameter rescaling. In the original
relative variable v and momentum a, the three coefficients of W are
\[
 W=r^2v^{-3}a^2+
                 (-2r^2v^{-2}+O(v^{-3}))a+
                              r^2v^{-1}+O(v^{-2}).          \tag{21}
\]
Consequently the coefficients of S have asymptotics
\[
 \begin{split}
 p_4&=\epsilon h^2r^4v^{-6},\qquad p_3=O(v^{-5}),\\
 p_2&=-3r+O(v^{-4}),\qquad p_1=4+O(v^{-3}),\qquad
 p_0=2v+O(v^{-2}).
 \end{split}                                               \tag{22}
\]
Over the generic constant field containing r and h these satisfy
(14), with g=-3r and s^2=epsilon*h^2*r^4. Thus an actual generic
finite rational-position mate commuting with R is impossible in
every order, without a pole-removal or centralizer-finiteness gate.
Trace also excludes finite algebraic coefficient-cover mates.

Thus none of the globally admitted family (18), for any scalar
ordering u,w, can supply the missing finite conjugate, in any
partner order. The same conclusion applies to regular scalar-series
orderings through their critical values, when global admission is
retained. This closes this ordered-square construction; different
negative Laurent critical operators, nonpositional first images and
genuinely nonscalar carriers remain open.

[Exact controls](check_ordered_quartic_infinity.py), with
[their receipt](ordered_quartic_infinity_controls.json), verify declared
fractional positive parts and the decisive coefficients with the ordering
parameters retained, the degree bounds, exceptional-ordering controls,
the commuting g=mu=0 boundary and the finite global rescaling. The
universal weighted-symbol and maximum-degree proof above supplies all
orders. No finite prefix is used as an unrestricted theorem. No bounded
full Long lift, DC1/DC2 solution, formal verification, independent
whole-proof replay, external review or novelty priority is claimed.
