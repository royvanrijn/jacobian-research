# Odd-numerator affine critical trace reduces to an explicit curve

The canonical affine-trace slice of a paired odd Laurent pencil has a
uniform linear polynomial operator. Its monic kernel is generically
unique, and the remaining normalized zero-value condition is one
accessory equation. In degrees ten and twenty-six, the only rational
slope where that triangular construction can fail is minus two;
that slope is impossible for a good model.

Consequently every rational degree-twenty-six model in this slice is
equivalent to a rational point on one explicit degree-twelve equation
over Q(t), subject to exact nonvanishing conditions. The good open is
nonempty: an exact point modulo29 lifts to Q29 and to a number field
of degree at most twelve. No rational Q point or global obstruction
to such points is established.

The scalar trace recurrence has classical Chebyshev background. With
`U_-1=0`, it is exactly `gamma_n(t)=U_(n-1)(t/2)` for `n>=0`, and
`k_j(t)=V_(2j)(t/2)`. Compare the Laurent identities in
[NIST DLMF (18.5.2)--(18.5.3)](https://dlmf.nist.gov/18.5) and
[the U recurrence, Table 18.9.1](https://dlmf.nist.gov/18.9).
For example `k_j(2 cos theta)=cos((2j+1/2)theta)/cos(theta/2)`,
with removable endpoint values. This identifies the recurrence and its
special values as classical; it does not identify the operator (12),
its accessory curve, or the rational roots of the combined diagonals (15)
with an existing classification theorem. Those remain separate arguments.

## 1. Normalization and the trace-square identity

Work in characteristic zero. Let h>=2,m=2h+1,L=4h+1, and consider

\[
 f(z)=N(z^2)/z,\qquad N\ \hbox{monic of degree }m,
 \qquad N(1)=N'(1)=0.
\]

A good model has

\[
        2uN'(u)-N(u)=L(u-1)S(u),                            \tag{1}
\]

with S monic squarefree of degree2h, S(0)S(1)!=0.
On its roots, M=N(u)^2/u has h distinct nonzero values, each twice.
For a paired pair u,v put

\[
 \rho=-N(u)N(v)/M,\qquad \theta=u+v,\qquad
 a=\frac{N(u)-N(v)}{u-v}.
\]

The [canonical quadratic chart](PAIRED_CRITICAL_TRACE_CHART.md) gives

\[
 uv=\rho^2,\quad N(u)=a(u-\rho),\quad
             \theta-2\rho=M/a^2.                           \tag{2}
\]

Here a may be nonconstant in the value algebra. We impose only

\[
                            \theta=s+t\rho.                \tag{3}
\]

If two pairs had the same rho, they would have the same theta
and the same two roots, contradicting squarefreeness of S.
Thus rho separates the pairs. Its monic polynomial q(R) has degreeh,
is squarefree, and q(0)!=0. The algebra B=k[R]/q is the value algebra.

Conversely start with a monic q of degreeh and scalars s,t. Define

\[
 S(u)=\operatorname{Res}_R
       \bigl(q(R),R^2-tuR+u^2-su\bigr),\qquad
 N(u)=\sum_j\frac{[u^j]\,L(u-1)S(u)}{2j-1}u^j.             \tag{4}
\]

Both are monic of their indicated degrees and satisfy(1).
In particular S(0)=q(0)^2 and N(0)=Lq(0)^2.
Equation N(1)=0 also implies N'(1)=0 by(1).

For formal roots u,v of u^2-(s+tR)u+R^2 put

\[
 C_N(R)=\frac{N(u)-N(v)}{u-v},\quad T_N(R)=N(u)+N(v),
 \quad D=T_N-\bigl(s+(t-2)R\bigr)C_N.                       \tag{5}
\]

These are polynomials, including at coincident formal roots.
Let gamma_0=0,gamma_1=1,gamma_(j+1)=t*gamma_j-gamma_(j-1), and write

\[
 k_j(t)=\gamma_{2j+1}(t)-\gamma_{2j}(t),\quad
                   k=k_h,\qquad K=L+kR.                    \tag{6}
\]

**Lemma 1.** Every good affine-trace model satisfies

\[
                              D=2Kq^2.                    \tag{7}
\]

Indeed D vanishes at each q-root by(2). Differentiation at that
root, using N'(u)=N(u)/(2u), uv=R^2 and u+v=s+tR, gives

\[
 T_N'=a(t-2)/2,\qquad
             C_N'=-a(t-2)/(2(s+(t-2)R)).
\]

No derivative of a or inverse of u' or v' is used.
It follows that D' also vanishes. Thus q^2 divides D.
The constant D(0)=2N(0)=2Lq(0)^2 fixes the constant quotient.
For the monomial N(u)=u^j, direct trace recurrence gives

\[
 D_j=2R^j\left[
 \gamma_j(t+s/R)-\gamma_{j-1}(t+s/R)\right]\quad(j\ge1).
 \tag{8}
\]

The leading coefficient of D is therefore2k, proving(7).

## 2. The uniform linear operator and its converse

For a polynomial q(R), reduce

\[
                    q(tu-R)=A_q(R)+B_q(R)u
             \pmod {u^2-(s+tR)u+R^2}.                      \tag{9}
\]

For the monomials q=R^j the remainders satisfy

\[
 A_0=1,\ B_0=0,\qquad
 A_{j+1}=-RA_j-tR^2B_j,\qquad
 B_{j+1}=tA_j+\bigl(ts+(t^2-1)R\bigr)B_j.                  \tag{10}
\]

Put e=s+(t+2)R and

\[
 J(q)=A_q\bigl((t+2)R(R-1)-2s\bigr)
 +R^2B_q\bigl(ts+(t+2)(t-1)R-t-2\bigr),                   \tag{11}
\]
\[
 {\cal L}_{h,s,t}(q)=LJ(q)
   -\bigl((t+2)RK-2eL\bigr)q-4eRKq'.                      \tag{12}
\]

This is linear in q. It is divisible by R^2; define
T_(h,s,t)(q)={\cal L}(q)/R^2.

**Theorem 2.** If t!=-2, then for every monic q of degreeh,
with S,N reconstructed by(4),

\[
                  D=2Kq^2\quad\Longleftrightarrow\quad
                              T_{h,s,t}(q)=0.             \tag{13}
\]

**Proof.** The useful identity, valid without good-open assumptions, is

\[
 H(E):=eR E'+\bigl((t+2)R/2-e\bigr)E=q{\cal L}(q),
              \qquad E=D-2Kq^2.                           \tag{14}
\]

Here is a derivation that also fixes its coefficients.
Choose formal z,w with z^2=u,w^2=v,zw=-R.
Then (z-w)^2=e, and the quadratic trace calculation gives

\[
             f(z)-f(w)=\frac{z-w}{2R}D.
\]

The other root, in R, of R^2-tuR+u^2-su is tu-R.
Consequently S(u)=q(R)q(tu-R), and similarly for v.
Differentiate the displayed f-difference with respect to R and
substitute f'(z)=L(z^2-1)S(z^2)/z^2.

To simplify, reduce (1-1/u)q(tu-R) into the basis1,u:

\[
 A_q-B_q-A_q(s+tR)/R^2
                        +\bigl(B_q+A_q/R^2\bigr)u.
\]

Taking its trace and divided difference gives precisely

\[
                             H(D)=LqJ(q).
\]

Differentiating2Kq^2 gives(14). The calculation on the generic
quadratic algebra proves a polynomial identity, so it remains
valid at every specialization of the polynomial coefficients.

If E is a nonzero polynomial of degree d, the leading coefficient
of H(E) is (t+2)(d-1/2) times that of E.
It cannot vanish in characteristic zero when t!=-2.
Thus H(E)=0 implies E=0. Conversely E=0 implies q{\cal L}(q)=0,
and q is nonzero. This proves(13).

For divisibility by R^2, recurrence(10) gives
A_q(0)=q(0), A_q'(0)=-q'(0). Substituting these into(11)--(12)
shows both {\cal L}(q)(0) and its derivative at zero vanish.
This also proves that the divided operator is polynomial. ∎

## 3. Triangular recurrence and the accessory curve

For0<=j<=h the operator T(R^j) has degree at mostj. Its diagonal is

\[
 \delta_j(t)=(t+2)\bigl(Lk_j(t)-(4j+1)k_h(t)\bigr).          \tag{15}
\]

For j=h this vanishes identically. For j<h it is a nonzero
polynomial: k_h is monic of degree2h while k_j has degree2j.

For completeness, the leading terms in(10), for j>=1, are
A_j=-gamma_(2j-1)R^j+... and B_j=gamma_(2j)R^(j-1)+....
They give(15). Homogeneity in s,R shows T(R^j) is a sum of
terms of total degreesj andj-1; hence its coefficient at R^i
has degree at mostj-i in s.

Let c_(ij)=[R^i]T(R^j). On the open where every delta_j(t),
j<h, is nonzero, the unique monic kernel is

\[
 q_h=1,\qquad
 q_i=-\delta_i(t)^{-1}\sum_{j=i+1}^h c_{ij}(s,t)q_j,
                      \quad i=h-1,\ldots,0.                \tag{16}
\]

All rows of the matrix are accounted for: rowh is zero and
the other h rows are solved in(16).
The coefficients q_i lie in Q(t)[s], of degree at mosth-i.
Define the accessory

\[
                           \Phi_h(s,t)=N(1).               \tag{17}
\]

It has degree at most2h in s. It is not identically zero.
At(s,t)=(0,2), the diagonals are16(h-j), and

\[
 q=(R+L)^h,\quad S=(u+L)^{2h},\quad N=(u+L)^{2h+1},
 \quad \Phi_h(0,2)=(4h+2)^{2h+1}\ne0.                       \tag{18}
\]

These closure data follow directly from(4),(7) and kernel uniqueness.
They are degenerate and do not establish good-open nonemptiness.

At exceptional slopes, equation T(q)=0 remains the exact linear
kernel problem; the unique-kernel formula(16) is simply unavailable.
No all-h assertion that such rational exceptions are absent is made.

## 4. The good open, with no missing descent condition

For a monic q and a solution of(13), work in B=k[R]/q and set

\[
 a=C_N(R)\bmod q,\quad d=s+(t-2)R,\quad
 \Delta=(s+tR)^2-4R^2,\quad M=a^2d.
\]

Let Q_M(Y) be the characteristic polynomial of multiplication by M
on B. Impose

\[
 q(0)\operatorname{Disc}(q)\operatorname{Disc}(S)S(1)\ne0,
 \tag{19}
\]
\[
 \operatorname{Norm}_{B/k}(a)\,
 \operatorname{Norm}_{B/k}(d)\,
 \operatorname{Norm}_{B/k}(\Delta)\,
 \operatorname{Norm}_{B/k}(M)\,
                      \operatorname{Disc}(Q_M)\ne0.        \tag{20}
\]

Some conditions are redundant; they are retained as explicit tests.

**Theorem 3.** For t!=-2, equations(13),(17) with Phi_h=0 and conditions(19)--(20)
are equivalent to good normalized affine-trace models.
All constructions and their inverse are rational over the base field.

**Proof.** Necessity follows from the canonical chart, Lemma1 and
Theorem2. Conversely, write N(u)=B_N+C_Nu in the quadratic
algebra. At each q-root, (7) gives D=2(B_N+RC_N)=0.
Thus N(u)=a(u-R), and its defining quadratic gives N(u)^2/u=a^2d=M.
The discriminant and norm conditions make the 2h critical roots
distinct, nonzero and away from1, with nonzero values.
Disc(Q_M)!=0 makes those h values distinct.
Phi_h=0 and(1) imply N(1)=N'(1)=0; differentiating(1) gives
2N''(1)=LS(1)!=0. The two critical points above u1 supply
the zero branch level.

The prescribed R is the canonical scale, since
N(u)N(v)=a^2(2R^2-R(s+tR))=-RM.
Moreover Disc(Q_M)!=0 makes M generate B, while R generates B by
definition. The rational basis reconstruction in the
[canonical trace chart](PAIRED_CRITICAL_TRACE_CHART.md), Theorem2,
recovers all these coefficients without labeling roots.
When h>=2,1,R are linearly independent in B, so s,t are uniquely
recovered as well. This proves the converse and descent statement. ∎

Every good model has simple paired finite critical points and
poles of orders4h+1 and1. For h>=2 its geometric monodromy is A_(4h+2):
the long cycle forces primitivity and finite inertia supplies a
double transposition. The exact group import is
[Jones, *Characteristically simple Beauville groups, I: cartesian powers of alternating groups*, Proposition6.1(3),p9](sources/jones-beauville-1304.5444v1.pdf),
valid in degrees at least nine. The resultant/leading-coefficient
formula of the
[discriminant bridge](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md)
then gives, with pencil parameter w,

\[
 \operatorname{Disc}_z(N(z^2)-wz)
              =L^L\bigl[wQ_M(w^2)\bigr]^2.                 \tag{21}
\]

When L is a square in the base field, arithmetic and geometric
groups coincide and the extension is regular. This includes
h2,n10 and h6,n26 over Q.

## 5. Complete rational-slope handling in degrees ten and twenty-six

**Lemma 4.** For h>=2, slope t=-2 is impossible for a good model,
over any characteristic-zero field.

**Proof.** Write E=D-2Kq^2 and let q_(h-1) be the subleading
coefficient of q. From the norm and Euler reconstruction,

\[
 [u^{m-2}]S=-hs+tq_{h-1},\qquad
 [u^{m-1}]N=\frac{L}{L-2}(-hs+tq_{h-1}-1).
\]

Using(8) and

\[
 \gamma_n(-2)=(-1)^{n-1}n,\qquad
 \gamma_n'(-2)=(-1)^n n(n^2-1)/6
\]

in the coefficient of R^(2h) at t=-2 gives

\[
             [R^{2h}]E=-\frac43 Lh(h-1)s.                  \tag{22}
\]

The identities for gamma follow by induction from its recurrence.
In(22) all coefficients of q cancel: the monic term contributes
-L*m*(m-1)*s/3, the Euler subleading term contributes
L*(m-1)*s+4L*q_(h-1)+2L, and the target square cancels the
last two terms. Thus(7) forces s0. But then the quadratic is
(u+R)^2, contradicting the distinct critical roots. ∎

**Proposition 5.** For h2 and h6, the only rational zero of any
diagonal(15) is t=-2.

For a proof using integer arithmetic only, put
D_j=Lk_j-(4j+1)k_h. If t=A/B in lowest terms has a denominator
prime ell, monicity gives

\[
 v_\ell(4j+1)-v_\ell(L)=2(h-j)v_\ell(B).                    \tag{23}
\]

For h6 the possible4j+1 are1,5,9,13,17,21, and the right-hand
steps are12,10,8,6,4,2. No positive valuation gap on the left
is a positive multiple of its step. For h2 the values are1,5
and the steps4,2, with the same conclusion. Thus every rational
root is integral. An integer root divides D_j(0).
Those constants are24,-30,16,-38,8,-46 for h6 and8,-14 for h2.
The retained certificate evaluates every signed divisor by the
integer gamma recurrence; exactly minus two vanishes in each row.
This proves the proposition without polynomial factorization. ∎

Consequently the good rational affine-trace locus in degree26
is exactly the good rational points of the single equation
Phi_6(s,t)=0 built by(16)--(17). No rational exceptional slope
or denominator stratum has been omitted. The retained input gives
every coefficient of q and Phi_6 exactly; its degree in s is twelve.
For h2 the accessory has degree four in s.

The known rational degree-ten model is a negative chart control,
not a positive one: its canonical rho is the scalar2/15 while
its two critical-pair traces differ, so no affine relation(3) exists.
For h2 any nonconstant canonical rho separates the two pairs,
and theta is automatically affine in it. The general h2 curve
therefore describes that entire nonconstant-scale open.

## 6. The degree26 good open is nonempty

At t2, the recurrence and accessory have coefficients integral
at29. Modulo29, s22 gives

\[
 q(R)=R^6+11R^5+R^4+22R^3+4R^2+18R+24,
\]
\[
 N(u)=u^{13}+13u^{12}+28u^{11}+22u^{10}+4u^9+11u^8
       +7u^7+10u^6+11u^5+18u^4+13u^3+15u^2+5u+16.
\]

The exact control gives

| Quantity | Value modulo29 |
| --- | ---: |
| q(0),Disc(q),Disc(S),S(1) | 24,23,1,8 |
| Norm(a),Norm(d),Norm(Delta) | 19,25,9 |
| Norm(M),Disc(Q_M) | 6,7 |
| Partial_s Phi_6 | 13 |

The secant is27+4R+23R^2+26R^3+26R^4+25R^5, hence is nonconstant.
The paired product R^2 is also nonconstant in this rank-six algebra.

The simple-root Hensel lemma, in
[Milne, *Algebraic Number Theory*3.00, Proposition7.31,p113](sources/milne-ant300.pdf),
applied to Phi_6(s,2) yields a root
alpha in Q29 congruent22, and all good conditions remain units.
Because this is a nonzero rational polynomial of degree twelve,
alpha is algebraic over Q of degree at most twelve.
Thus the good open has a number-field point in Q29, as well as
a local point over Q29. Equation(21) and the group argument give
full regular A26 for the resulting covers.
No inference that alpha or any other good parameter is rational
over Q is made.

## Evidence and remaining question

The [checker](verify_odd_affine_critical_trace_curves.py) rebuilds all
h2/h6 curve coefficients in Q(t)[s], verifies the full triangular
kernel, and checks the universal Euler/operator identity with generic
q in bounded degrees. It replays the complete rational-slope
certificate and every good-unit test above, including the entire
discriminant identity over F29(w). The written proof supplies
uniformity; finite controls do not prove it by sampling.
The [explicit curve input](odd_affine_critical_trace_curve_inputs.json),
[receipt](odd_affine_critical_trace_curves_verification.json), and
[source ledger](odd_affine_critical_trace_curves_source_ledger.json)
retain the exact equations and scope.

The unresolved arithmetic question is whether the indicated good
open of Phi_6(s,t)=0 has a Q-rational point. This note neither
constructs one nor excludes all of them. It treats N(z^2)/z with
odd-degree N, not arbitrary degree26 polynomial pencils.
No rational parametrization, novelty claim, formal verification,
external human review or polynomial Keller realization of its
generic alternating extension is asserted.

