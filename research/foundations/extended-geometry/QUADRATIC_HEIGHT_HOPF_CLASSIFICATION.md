# Complete quadratic-height Hopf classification

Date: 19 November 2026.

This note classifies all regular winding-one phases with quadratic
height, including a linear height term. The narrower ordinary-ideal
presentation remains in [HOPF_LIFT_CLASSIFICATION.md Section8](HOPF_LIFT_CLASSIFICATION.md#8-complete-complex-polynomial-quadratic-height-classification).
The proof below is a classification over characteristic-zero fields;
it makes no new ideal-theoretic or scheme-reducedness claim. Historical
priority, external review and formal verification are not asserted.

The only unsafe branch has even least balanced degree, with one fixed
marker surviving at every moment order. In particular this entire
class supplies no primitive quintic or septic Image counterexample.
Arbitrary quintics, septics and higher-height phases remain outside
this classification.

Work over any field K of characteristic zero, in
`S=K[x,y,t]/(t²+2xy-1)`, with

\[
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]

The integral is the rational polynomial functional, not an ordered-field
positivity premise. Call p pure when every positive pure moment vanishes,
and safe when `Lambda(g p^m)=0` eventually for every fixed `g in S`.
Let

\[
p=x^{-1}(C(x)+tB(x)+t^2D(x)),\quad C,B,D\in K[x],
\quad C(0)=c_0.
\]

Regularity in S is equivalent to `B(0)=0,D(0)=-c0`: reduce the numerator
modulo x, where the relation is `t²=1`. Under this premise the literal
polynomial representative is

\[
p=2c_0y+\sum_{i\ge1}c_i x^{i-1}
       +t\sum_{i\ge1}b_i x^{i-1}
       +t^2\sum_{i\ge1}d_i x^{i-1}.
\]

## Main theorem

**Theorem QuadraticHeightHopfClassification.** If `c0=0`, the phase is
pure exactly when C,B,D are all divisible by x²; this branch is safe
for every fixed marker, as proved in Section7.

For the remaining case `c0!=0`, put `A=C-B²/(4D)` in K(x), a unit
formal series at zero. Then
`Lambda(p^m)=0` for all `m>=1` if and only if

\[
D=-\frac{A^3}{(A-xA')^2}.                                      (1)
\]

When C,B,D are polynomials, every solution falls into exactly one branch.

1. **Constant center.** `A=a in K*`, `D=-a`, and
   `C=a-B²/(4a)`, with arbitrary `B in xK[x]`. This branch is safe:
   for every fixed `g in S`, `Lambda(g p^m)=0` for all sufficiently
   large m. Its least polynomial degree is odd.

2. **Nonconstant rational center.** There are unique `n>=1` and
   `c,lambda in K*` such that, with `U=1+lambda*x` and
   `S_n(z)=((1+z)^n-1)/z`,

   \[
   u=\frac{x}{A}=c(1-U^{-n}),\qquad
   A=\frac{U^n}{c\lambda S_n(\lambda x)},\qquad
   D=-\frac{S_n(\lambda x)U^{n+2}}{c n^2\lambda}.                (2)
   \]

   The polynomial B is constrained precisely by `B(0)=0` and
   polynomiality of `C=A+B²/(4D)`; equivalently

   \[
   D\mid B^2-\frac{4U^{2n+2}}{c^2n^2\lambda^2}.                (3)
   \]

   Every such member is unsafe with the SAME fixed marker x:

   \[
   \Lambda(xp^m)
    =\kappa_m\frac{(1/n)_m}{\lambda c^m m!}\ne0
       \quad(m\ge1),\qquad
   \kappa_m=\frac{4^m(m!)^2}{(2m+1)!}.                         (4)
   \]

   Here `(a)_m=a(a+1)...(a+m-1)` is the rising factorial.

   Its least polynomial degree is EVEN. More precisely, if
   `b=deg B` (with `deg0=-infinity`), it is

   \[
   d=\begin{cases}
       2n+2,&b\le2n+2,\\
       2b-2n-2,&b\ge2n+3.
      \end{cases}                                             (5)
   \]

Consequently there is no radial-primitive unsafe odd-degree5 or7 phase,
or any unsafe phase of odd least degree, in this entire winding-one
quadratic-height class. This does not classify phases with cubic or
higher height or larger negative angular winding. Section7 treats
`C(0)=0` separately.

## 1. The exact formal purity criterion

Choose the unique formal square root `X=sqrt(-D/A)` with `X(0)=1`,
and put `Y=-B/(2AX)`. Then `Y(0)=0`, and

\[
C=A(1-Y^2),\quad B=-2AXY,\quad D=-AX^2,
\quad C+tB+t^2D=A(1-(Xt+Y)^2).
\]

Let `J_m` be the odd polynomial primitive of `(1-v²)^m`. The interval
average of `(C+tB+t²D)^m` is

\[
\frac{A^m}{2X}\big(J_m(Y+X)-J_m(Y-X)\big).
\]

Because `J_m(v)-J_m(1)` is divisible by `(v-1)^(m+1)`, and similarly
at -1, this is congruent to `kappa_m A^m/X` modulo `x^(m+1)`.
The interval constant is `kappa0=1` and obeys
`kappa_m=2m*kappa_(m-1)/(2m+1)` by integration by parts, giving(4).
Thus the mth pure moment vanishes exactly when

\[
[x^m]A^m L=0,\qquad L=1/X.
\]

For a given unit A and `L(0)=1`, these conditions determine every next
coefficient of L uniquely: its coefficient in the mth equation is
`A(0)^m`, a unit. The unique solution is

\[
L=1-xA'/A,
\]

since `[x^m](A^m-x(A^m)'/m)=0`. This proves (1) in both directions.
No finite-prefix inference occurs.

## 2. Rational pole bookkeeping gives the complete parameter list

Set `u=x/A=P/Q` in lowest terms. It has a simple zero at zero and
`Q(0)!=0`. Equation(1) becomes

\[
D=-\frac{u}{x(u')^2}.                                        (6)
\]

Work temporarily over an algebraic closure. At a nonzero finite zero
of u of multiplicity e, the valuation of D is `2-e`. Thus e<=2.
If e=2, D is regular and nonzero there, whereas A=x/u has a double
pole. Polynomial B cannot make `A+B²/(4D)` regular, contradicting
polynomial C. Hence P is squarefree. At every finite point where u
is nonzero and finite, u' cannot vanish, or D would have a pole.

Write `V=P'Q-PQ'`. At a root of Q of multiplicity e, V has multiplicity
e-1; P and Q are coprime. V has no zeros at P's simple roots and no
other zeros. If `l` is the number of distinct roots of Q, then

\[
\deg V=\deg Q-l.                                             (7)
\]

If `s=deg P` differs from `t=deg Q`, the leading coefficient of V
is `(s-t)P_s Q_t`, giving `deg V=s+t-1`. Equation(7) forces
`s+l=1`, hence `s=1,l=0`. Therefore Q is constant, u is linear,
and A is constant: branch1.

Otherwise `s=t=n>=1`. Let `c=P_n/Q_n` and `k=deg(P-cQ)>=0`.
The expansion at infinity gives `deg V=k+n-1`. Together with(7),
this says `k=1-l`. Since `l>=1`, it follows that `l=1,k=0`.
Consequently Q has one root of multiplicity n and P-cQ is a nonzero
constant. The simple zero `u(0)=0` now gives exactly(2).

These parameters descend to K: the single root of Q is determined by
the ratio of its leading two coefficients, `c=P_n/Q_n`, and n is its
degree. Substitution into(6) gives the stated D. Finally
`AD=-U^(2n+2)/(c² n² lambda²)`, so polynomiality of C is precisely(3).
This is an elementary polynomial-degree argument; it does not import
a Belyi-map or Riemann--Hurwitz classification.

## 3. The fixed marker survives at every exponent

Multiplication by x means extracting coefficient `x^(m-1)`. The same
endpoint congruence used in Section1 gives

\[
\Lambda(xp^m)=\kappa_m[x^{m-1}]A^m(1-xA'/A)
             =\frac{\kappa_m}{m}[x^{m-1}]A^m.                 (8)
\]

For completeness the needed Lagrange coefficient identity follows from
formal residues. If `v(z)` is the inverse of `u=x/A`, then

\[
[z^m]v(z)=\operatorname{Res}_x x u^{-m-1}u'\,dx
         =\frac1m\operatorname{Res}_x u^{-m}\,dx
         =\frac1m[x^{m-1}]A^m.
\]

The middle equality is the vanishing residue of `(x u^-m)'`. For the
change of variable, a monomial `u^j u'` is a derivative and has zero
residue when `j!=-1`; for `j=-1`, `u'/u=1/x+regular` has residue1.
Linearity then proves the needed substitution coefficientwise. Here

\[
v(z)=\lambda^{-1}\big((1-z/c)^{-1/n}-1\big),
\]

whose binomial series proves(4). Every displayed factor is nonzero in
characteristic zero. The result does not depend on B.

## 4. Constant centers are safe by a polynomial measure-preserving shear

Write `B=x b(x)` and `a=A`. Then

\[
p=2ay+b(x)t-\frac{x b(x)^2}{4a}.
\]

The locally nilpotent derivation `delta=x*d/dt-t*d/dy` preserves
`t²+2xy`; multiplying it by any `h(x) in K[x]` remains locally
nilpotent. Its exponential sends

\[
t\mapsto t+xh,\qquad y\mapsto y-th-xh^2/2.
\]

Choose `h=-b/(2a)`. Then `p=2a exp(h delta)(y)`.
This automorphism preserves Lambda: on the localized chart its
infinitesimal action is `x h(x) d/dt`, whose interval integral is the
boundary term

\[
\operatorname{CT}_x\frac{x h(x)}2(f(x,0,1)-f(x,0,-1))=0.
\]

Both boundary evaluations are polynomials in x; the additional x kills
their constant terms. Exponentiating is a finite sum on each polynomial.
For any fixed g, move the inverse automorphism onto g. Its chart
expansion is a finite Laurent polynomial in x, whereas `y^m` contributes
the strictly negative factor `x^-m`. For sufficiently large m their
product has no x-constant term. This proves safety for every fixed marker.

## 5. Exact least-degree parity and the Image normalization

For c0!=0, the least polynomial degree of a regular quadratic-height phase is

\[
\max(1,\deg C-1,\deg B,\deg D+1).
\]

Indeed the displayed regular representative has that degree. For degree
at least2 its top part lies in `K[x,t]` and is nonzero, so it is not
divisible by the leading quadric `t²+2xy`; therefore the sphere relation
cannot reduce that degree. In degree1 the nonzero term `2c0 y` gives
the same conclusion.

In branch1, a nonzero B of degree b>=1 gives degree `2b-1`; B=0 gives
degree1. In branch2, D has degree `2n+1`, A has degree1 at infinity,
and `C=A+B²/(4D)` has degree at most `max(1,2b-2n-1)`. If
`b>=2n+3`, the latter degree is exact and strictly dominates the other
terms. This proves(5), including all possible leading cancellations.

Use the repository bilinears

\[
R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad W=2\xi_2z_1,
\quad T=\xi_1z_1-\xi_2z_2,
\quad T^2+2ZW=R^2.
\]

For least degree d, homogenize the regular p to `F_d=R^d p(Z/R,W/R,T/R)`.
It is radial-primitive because its top sphere-filtration term is nonzero
modulo `T²+2ZW`. The exact two-pair factorial normalization gives

\[
\mathcal E_2(F_d^m)=0,\qquad
\mathcal E_2(ZF_d^m)=(dm+2)!\,
   \kappa_m\frac{(1/n)_m}{\lambda c^m m!}\ne0
\]

in branch2. To transfer safety to every fixed original polynomial marker,
split it into bihomogeneous components Q of bidegree (u,v). If u>v its
contraction against F_d^m is zero by degree. If u<=v, the output is
homogeneous in z of degree v-u. For every multi-index alpha of size v-u,

\[
\partial_z^\alpha\mathcal E_2(QF_d^m)
       =\mathcal E_2(\xi^\alpha QF_d^m).
\]

The marker on the right is balanced, and its contraction eventually
vanishes by angular safety and the homogeneous factorial identity.
There are finitely many alpha, so all top derivatives of the output
eventually vanish and hence so does that homogeneous output.
Taking the maximum cutoff over the finitely many bihomogeneous
components handles an arbitrary original marker. This applies to both
the constant-center and the zero-center safe branches.

## 6. An explicit primitive degree-six witness

Take `n=2,c=1/4,lambda=2`, and

\[
C=1+3x,\quad B=2x(1+2x)^2,\quad
D=-(1+x)(1+2x)^4.
\]

Then `A=(1+2x)^2/(1+x)`, the phase has least degree6, all pure moments
vanish, and the first marked moments from(4) are `2/3,8/5,32/7,...`.
These finite values are controls of the unrestricted formulas, not their
proof. The next research target for primitive odd5/7 must leave this
whole quadratic-height winding-one class.

**Corollary QuadraticHeightDegreeSix.** An explicit homogeneous lift is

\[
\begin{split}
F_6={}&2WR^5+3R^6+2T(R+2Z)^2R^3\\
 &-T^2(9R^4+32ZR^3+56Z^2R^2+48Z^3R+16Z^4).
\end{split}
\]

It has balanced degree6, and for every m>=1,

\[
\mathcal E_2(F_6^m)=0,\qquad
\mathcal E_2(ZF_6^m)=(6m+2)!\frac{2^{2m-1}}{2m+1}\ne0.
\]

Thus every positive power lies in the two-pair Image kernel while one
fixed marker survives at every exponent. It is radial-primitive:
at `(xi1,xi2,z1,z2)=(1,1,1,-1)`, R=0 but `F6=-64`.
It is not a proper power even up to a scalar over field extensions.
Indeed the original-coordinate specialization
`(xi1,xi2,z1,z2)=(1,1,1,s)` is exactly

\[
-4s(11s^2+8s+1)(4s^3-3s^2-4s-1),
\]

which has a simple zero at s=0, with derivative4. A nonzero scalar times
a kth power, k>=2, cannot have a simple specialization zero.

## 7. The zero-center branch is strictly positive-weight and safe

Now allow `C(0)=0`. Regularity also forces `B(0)=D(0)=0`, so p has
only nonnegative angular powers. Its weight-zero part is

\[
h(t)=c+bt+dt^2,\quad(c,b,d)=(c_1,b_1,d_1).
\]

For each m the pure moment is exactly the interval moment of `h^m`;
positive angular weights do not contribute its constant term. The first
three moments force h=0 over every characteristic-zero field. Indeed,
the first gives `c=-d/3`, the second then gives

\[
b^2/3+4d^2/45=0.
\]

Under those two equations the third is

\[
-256d^3/4725=0.
\]

Thus `d=b=c=0`. Equivalently C,B,D are all divisible by x², making
every angular term of p strictly positive. Conversely this condition
gives all pure moments zero and eventual mixed vanishing for every
fixed polynomial marker, by its bounded negative angular support.

Together with Sections1--5 this proves the complete classification for
ALL regular polynomial phases `x^-1(C+tB+t²D)`: the only unsafe branch
is the rational-center family(2)--(3), its fixed marker x survives at
every exponent, and its least balanced degree is always even.


## Sources and finite verification

The angular/factorial normalization is the retained two-pair identity
from [TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
and its self-contained use in
[ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md](ODD_PRIMITIVE_TWO_PAIR_IMAGE_WITNESSES.md).
The narrower B=0 formal classification and ordinary-ideal result are
[HOPF_LIFT_CLASSIFICATION.md Section8](HOPF_LIFT_CLASSIFICATION.md#8-complete-complex-polynomial-quadratic-height-classification).
The endpoint congruence, rational-function pole/degree classification,
formal-residue coefficient identity, shear invariance and zero-center
calculation needed here are proved explicitly above. The known even
primitive Image examples are antecedents, not a novelty claim for even
degree6.

The [source ledger](quadratic_height_hopf_classification_source_ledger.json)
records these scopes and exact retained-source pins. The standard-library
[checker](verify_quadratic_height_hopf_classification.py) tests separate
raw interval expansions against the all-order marker formula for several
parameter and shear choices; checks the zero-center reductions and
corrupted coefficients; and expands the ORIGINAL four-variable F6 to
compare pure/marked factorial contractions through m8. It also checks
the radial value and simple specialization zero.

The [receipt](quadratic_height_hopf_classification_verification.json)
records bounded controls of these identities. The unrestricted theorem,
all-order assertions and polynomial safety follow from the written
proof and are not inferred from a vanishing finite prefix.

```bash
python3 research/foundations/extended-geometry/verify_quadratic_height_hopf_classification.py \
  --write-receipt research/foundations/extended-geometry/quadratic_height_hopf_classification_verification.json
```
