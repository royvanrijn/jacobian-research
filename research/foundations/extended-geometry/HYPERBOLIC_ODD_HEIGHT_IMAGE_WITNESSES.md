# Primitive Image witnesses in prime balanced degrees

Date: 19 November 2026.

All-order identities below are proved by formal algebra. Finite
computations are separate controls; no historical-priority claim is made. The odd cubic subfamily
in Sections1--5 supplies the explicit prime balanced-degree13 example.
Section6 extends it to all odd starting height powers and gives, for
example, prime balanced degrees23,41 and59.

## 1. Family and all-order statement

Let K be any characteristic-zero field, `n>=1` an integer and `lambda in K*`.
Set

\[
S_k(z)=\frac{(1+z)^k-1}{z}
      =\sum_{j=0}^{k-1}\binom{k}{j+1}z^j,\quad U=1+\lambda x,
\]

\[
C=U^{n+1},\qquad X=\frac{U S_{2n}(\lambda x)}{2n},\qquad
D=-CX^2=-\frac{U^{n+3}S_{2n}(\lambda x)^2}{4n^2}.             (H1)
\]

On the sphere `t²+2xy=1`, the phase

\[
p=\frac{t}{x}(C+t^2D)
  =\frac{tC}{x}(1-t^2X^2)                                  (H2)
\]

is polynomial: X(0)=1, and
`(1-t²X²)/x=2y-t²(X²-1)/x` has no denominator. Its least degree is
`d=5n+3`, since D has degree5n+1 and the top term of p is
`-lambda^(5n+1)t³x^(5n)/(4n²)`, nonzero modulo the leading quadric.

Use the standard two-pair bilinears

\[
R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad W=2\xi_2z_1,
\quad T=\xi_1z_1-\xi_2z_2,\quad T^2+2ZW=R^2.
\]

Define the homogeneous degree2n polynomial

\[
Y_h=(R+\lambda Z)\sum_{j=0}^{2n-1}\binom{2n}{j+1}
             \lambda^j Z^jR^{2n-1-j}.                       (H3)
\]

Its value at Z=0 is `2n R^(2n)`, so the quotient in

\[
B_h=8n^2 W R^{4n}
       -T^2\frac{Y_h^2-4n^2R^{4n}}{Z}                       (H4)
\]

is a literal polynomial, not a localization. Finally put

\[
\boxed{F_{n,\lambda}=T(R+\lambda Z)^{n+1}B_h,\quad d=5n+3.}  (H5)
\]

This is `4n²` times the balanced-degree d homogenization of p.

**Theorem HyperbolicOddCubicImage.** For every allowed n,lambda and every
positive m,

\[
\mathcal E_2(F_{n,\lambda}^{m})=0.
\]

For the ONE fixed original marker Z, odd mixed moments vanish and every
even mixed moment is nonzero. More precisely, for m=2r>=2,

\[
\boxed{\mathcal E_2(ZF_{n,\lambda}^{2r})
  =(2rd+2)!\,(4n^2)^{2r}\,\kappa_{2r}\,\mu_{n,\lambda,r}\ne0,} (H6)
\]

where

\[
\kappa_m=\int_0^1 t^m(1-t^2)^m\,dt
       =\frac{2^m m!}{\prod_{j=0}^m(m+2j+1)}\quad(m\text{ even}),
\]

\[
\mu_{n,\lambda,r}
 =\frac{(-n^2\lambda^2)^r}{\lambda}
    \frac{(1/(2n))_r(1-1/(2n))_r}{(1/2)_r r!}.              (H7)
\]

Here `(a)_r=a(a+1)...(a+r-1)`. In addition, `R` does not divide F, and F
is not a proper power even up to a nonzero scalar over field extensions.
Thus this is a family of radial-primitive, non-proper-power failures of
the two-pair Image Mathieu property in balanced degree5n+3.

In particular n=2 gives PRIME balanced degree13, and n=4 gives PRIME
balanced degree23. The ordinary total degrees are twice the balanced
degrees. No assertion that13 is the minimum prime degree is made;
degrees5,7,11 and other primes are not settled by this construction.

## 2. Endpoint averaging and the formal inverse coordinate

Write

\[
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]

Let `A=C/X,L=1/X,u=x/A`, and let `v(u)` be the inverse formal series
of u(x). For the odd cubic-height phase, substitution of `Xt` for height
and order m+1 endpoint contact give, for even m,

\[
\Lambda(p^m)=\kappa_m[x^m]A^mL,
\qquad \Lambda(xp^m)=\kappa_m[x^{m-1}]A^mL.                 (H8)
\]

The full odd-cubic endpoint proof is retained in
[ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md); here it is also immediate
from a polynomial primitive of `s^m(1-s²)^m`, whose derivative has order m
zeros at ±1. Odd m moments vanish by t-parity. Integration by parts of
`t^(a+1)(1-t²)^m`, followed by iteration with a fixed and then a=m,
gives the rational product for kappa in(H7).

Formal residue substitution gives

\[
[x^m]A^mL=[u^m]\Phi(u),\qquad
[x^{m-1}]A^mL=[u^m]v(u)\Phi(u),
\quad \Phi(u)=\frac{v'(u)}{C(v(u))}.                        (H9)
\]

For example the second identity is the substitution in
`Res_x L u(x)^(-m) dx`; the first uses `dx/x`. The residue change is
valid over any characteristic-zero field and follows on monomials from
zero residues of derivatives and residue1 of `u'/u`.

The explicit family has

\[
A=\frac{2nU^n}{S_{2n}(\lambda x)},\quad
u=\frac{U^n-U^{-n}}{2n\lambda},\quad
V(x)=\int_0^x\frac{ds}{C(s)}=\frac{1-U^{-n}}{n\lambda}.       (H10)
\]

The square root with constant term1 satisfies

\[
V(v(u))=u+\frac{1-\sqrt{1+n^2\lambda^2u^2}}{n\lambda}.
\]

Consequently

\[
\Phi(u)=1-\frac{n\lambda u}{\sqrt{1+n^2\lambda^2u^2}}.       (H11)
\]

It has no positive even powers. Equations(H8)--(H9) therefore prove
ALL pure moments zero. This is a constructive example where u is
rational in x but Phi is algebraic and not rational in u; requiring
rational Phi would omit it.

## 3. A self-contained nonvanishing recurrence for the fixed marker

Set `z=n lambda u`, `s=arcsinh(z)` as a formal series at zero, and
`q=(n-1)/n`. Formally `sinh(s)=z`, `cosh(s)=sqrt(1+z²)` and

\[
U(v(u))=\exp(s/n),\quad
v(u)=\lambda^{-1}(\exp(s/n)-1),\quad
\Phi(u)=1-z/\sqrt{1+z^2}.
\]

Taking even parts and using the elementary addition identity for cosh,

\[
\frac{v(u)\Phi(u)+v(-u)\Phi(-u)}2
 =\frac1\lambda\left(
 \frac{\cosh(q\operatorname{arcsinh}z)}{\sqrt{1+z^2}}-1\right). (H12)
\]

This identity works also for n=1; no division by n-1 is used.

Put `w(z)=cosh(q arcsinh z)/sqrt(1+z²)`. Direct formal differentiation
gives

\[
(1+z^2)w''+3zw'+(1-q^2)w=0,
\qquad w(0)=1,\quad w'(0)=0.                               (H13)
\]

One way to check it is to start from `g=cosh(q arcsinh z)`, for which
`(1+z²)g''+zg'-q²g=0`, and substitute `g=sqrt(1+z²)w`.
Write `w=sum_(r>=0)a_r z^(2r)`. Its exact recurrence is

\[
a_0=1,\qquad
a_{r+1}=-\frac{(2r+1)^2-q^2}{(2r+2)(2r+1)}a_r.              (H14)
\]

Every factor is nonzero: `0<=q<1` as a rational number. Therefore
all a_r are nonzero in characteristic zero, and multiplication of the
recurrence yields

\[
a_r=(-1)^r\frac{(1/(2n))_r(1-1/(2n))_r}{(1/2)_r r!}.
\]

The coefficient of u^(2r) in(H12) is exactly(H7). Equations(H8)--(H9)
now prove the normalized all-order marked moments, with no asymptotic
or finite-prefix inference. The product is an identity of rational
numbers and powers of lambda, so the nonvanishing holds in every
characteristic-zero field. Hypergeometric theory is not a proof input.

## 4. Polynomial lift, scaling, radial primitivity and proper powers

The polynomial `Y_h` has homogeneous degree2n and
`Y_h(0,Z)=lambda^(2n)Z^(2n)`. Consequently B_h has degree4n+1, and F
has degree5n+3. Setting R=1 gives `Y_h=2nX`, so(H4) becomes

\[
B_h|_{R=1}=4n^2\left(2y-t^2\frac{X^2-1}{x}\right).
\]

Thus F dehomogenizes to EXACTLY `4n² p`. The original contraction of
a balanced homogeneous polynomial H of degreeD is
`E2(H)=(D+1)! Lambda(H|R=1)`. Applied to F^m and ZF^m, this supplies
the factors `(dm+1)!`, `(dm+2)!` and the crucial scaling `(4n²)^m`.
It proves(H6) in the original four variables.

At the original point `(xi1,xi2,z1,z2)=(1,1,1,-1)`, one has
`R=0,Z=-1,W=2,T=2`, and

\[
F_{n,\lambda}=8(-1)^{n+1}\lambda^{5n+1}\ne0.                (H15)
\]

This proves radial primitivity for ALL allowed lambda.

The bilinear T is absolutely irreducible: its coefficient matrix has
rank2, whereas a nontrivial bilinear factorization would have rank1.
Its polynomial quotient is therefore an integral domain. After dividing
F by T and reducing modulo T, the remaining product is

\[
8n^2 W R^{4n}(R+\lambda Z)^{n+1}\pmod T.                    (H16)
\]

Neither R, W nor R+lambda Z is divisible by T. For R and R+lambda Z
this follows from their nonassociation with the irreducible T; W is a
product of two linear forms. Thus(H16) is nonzero in the prime T
quotient for every lambda, including lambda=-2. The exact T multiplicity
of F is1. Hence F cannot be a nonzero scalar times a proper power,
even after extending the coefficient field.

## 5. An explicit prime balanced-degree13 member

**Corollary PrimeDegreeThirteenImage.** The a=1,n=2 member below is an
explicit radial-primitive, non-proper-power two-pair Image witness of
prime balanced degree13, with fixed original marker Z.

For n=2,lambda=1,

\[
Y_h=4R^4+10R^3Z+10R^2Z^2+5RZ^3+Z^4,
\]

\[
\begin{split}
Q_h={}&80R^7+180ZR^6+240Z^2R^5+208Z^3R^4\\
 &+120Z^4R^3+45Z^5R^2+10Z^6R+Z^7.
\end{split}
\]

The integer-coefficient witness is

\[
\boxed{F_{13}=T(R+Z)^3(32WR^8-T^2Q_h).}                    (H17)
\]

It has `F13=-8` at the R=0 point above, and `(F13/T)=442368` at
`(xi1,xi2,z1,z2)=(1,1,1,1)`. Its normalized ORIGINAL marked values at
m=2 and4 are

\[
\frac{\mathcal E_2(ZF_{13}^2)}{28!}=-1024/35,
\qquad
\frac{\mathcal E_2(ZF_{13}^4)}{54!}=1048576/429.
\]

These are controls of the all-order formula, not its proof. This member
shows in particular that a global obstruction to prime balanced degrees
would be false. It leaves the smaller prime degrees unsettled outside
the previously classified restricted classes.

## 6. All odd starting height powers

Let a>=1 be ODD and keep n,lambda,U,X,A,L from(H1)--(H10). Define

\[
C_a=U^{n+1}X^{a-1},\quad D_a=-U^{n+1}X^{a+1},\quad
p_a=\frac{t^a}{x}(C_a+t^2D_a)
    =\frac A x (Xt)^a(1-(Xt)^2).
\]

This is polynomial by the same endpoint divisibility. Its exact least
balanced degree is

\[
d=(2a+3)n+a+2,\qquad 2d-1=(2n+1)(2a+3).                    (H18)
\]

The cleared original lift is

\[
\boxed{F_{a,n,\lambda}
 =T^a(R+\lambda Z)^{n+1}Y_h^{a-1}B_h,}                     (H19)
\]

using(H3)--(H4). It dehomogenizes to `(2n)^(a+1) p_a`.

**Theorem HyperbolicOddHeightImage.** For all odd a>=1, n>=1 and
lambda!=0 over every characteristic-zero field,

\[
\mathcal E_2(F_{a,n,\lambda}^m)=0\quad(m\ge1).
\]

The same fixed marker Z has zero mixed moments at odd m and, for m=2r,

\[
\boxed{\mathcal E_2(ZF_{a,n,\lambda}^{2r})
 =(2rd+2)!\,(2n)^{(a+1)2r}\,
      \kappa_{a,2r}\,\mu_{n,\lambda,r}\ne0,}                (H20)
\]

where mu is(H7) and

\[
\kappa_{a,m}=\int_0^1 v^{am}(1-v^2)^m\,dv
 =\frac{2^m m!}{\prod_{j=0}^m(am+2j+1)}\quad(m\text{ even}).
\]

Every F is radial-primitive and is not a proper power up to a nonzero
scalar over any coefficient-field extension. Thus(H19) supplies the
broader degree set(H18), including prime41 at(a,n)=(3,4) and prime59 at
(a,n)=(5,4), in addition to prime13 and23 from a=1.

For an odd integer d, this family realizes d if and only if
`2d-1=u*v` with integers `u,v>=5` and `u=v=1 (mod 4)`.
Indeed, (H18) supplies `u=2n+1` and `v=2a+3`; odd d forces n even,
while odd a gives the asserted congruence for v. Conversely, such factors
give `n=(u-1)/2>=2` and the odd integer `a=(v-3)/2>=1`, and (H18)
then gives exactly d. This is an equivalence for the displayed family;
it is not a necessary condition for arbitrary primitive Image witnesses.

**All-order moments and degree.** The height profile is now the odd
polynomial `v^a(1-v²)`. Its mth power still has order m endpoint zeros
at ±1. The endpoint/residue calculation is identical, with kappa replaced
by kappa_(a,m), while A,L,u,Phi and the even marked coefficient mu are
UNCHANGED. This proves(H20); odd a ensures odd m parity vanishing.
The top coefficient of D_a is nonzero, with degree `(2a+3)n+1`, so p_a
has the exact degree(H18). In(H19), `Y_h=2nX` after R=1 and
`B_h=4n²(1-t²X²)/x`, giving the stated clearing factor.

**Radial primitivity.** At the same original R=0 point,

\[
F_{a,n,\lambda}
 =(-1)^{n+1}2^{a+2}\lambda^{(2a+3)n+1}\ne0.                 (H21)
\]

**Proper powers, including noncoprime parameters.** Specialize original
coordinates to `(xi1,xi2,z1,z2)=(1,0,1,z)`. Then R=T=1,W=0. Write
`U=1+lambda*z`, `S(U)=1+U+...+U^(2n-1)` and
`Y(U)=U+...+U^(2n)=U S(U)`. The specialized polynomial is

\[
-\lambda U^{n+a}S(U)^{a-1}
                  \frac{Y(U)^2-4n^2}{U-1}.                 (H22)
\]

The displayed quotient is polynomial, since Y(1)=2n. This formula
includes the factor lambda from division by Z; it is nonzero by hypothesis.

For v>1,

\[
Y(-v)=\frac{v(v^{2n}-1)}{v+1}.
\]

Both factors `v/(v+1)` and `v^(2n)-1` are positive and strictly increasing.
Thus Y(-v) increases strictly from0 to infinity, with positive derivative.
Equivalently, the equation `Y(-v)=2n` has the unique simple root v>1
of

\[
P_n(v)=v^{2n+1}-(2n+1)v-2n.
\]

One has P_n(1)=-4n<0, P_n(2)>=0 and
`P_n'(v)=(2n+1)(v^(2n)-1)>0` for v>1. The inequality at2 follows by
induction from `2^(2n+1)>=6n+2`, starting at n=1.

At the resulting algebraic number U0=-v<-1, neither U0, U0-1 nor
S(U0)=2n/U0 is zero, whereas `Y²-4n²` has a simple zero because
Y(U0)=2n and Y'(U0)!=0. Hence(H22) has a simple zero. This excludes
every proper-power factorization, regardless of gcd(a,n).

The real argument only establishes a nonzero algebraic root and its
nonzero derivative for a polynomial over Q. Transport that root into
an algebraic closure of any characteristic-zero coefficient extension;
the same equalities and nonvanishings persist. Since lambda!=0, passing
from U to the original z is an invertible affine substitution and
preserves simplicity. A scalar times a kth power for k>=2 cannot have
this simple specialization zero. This proves the assertion over every
allowed field and extension.

For reference, the exact T-adic order is a. Modulo the prime T, the
quotient by T^a is
`8n² W R^(4n)(R+lambda Z)^(n+1)Y_h^(a-1)`, which is nonzero. Indeed
R and Z are algebraically independent on the open part Z!=0 of T=0:
take `xi1=1,z2=Z,xi2=R/(2Z),z1=R/2`. Thus the nonzero polynomial
Y_h(R,Z) does not vanish identically modulo T. The simple-zero argument,
rather than T-order alone, excludes powers when a>1.

## 7. Sources, antecedents and verification scope

Classical primary formulas were inspected before writing a special-function
normalization claim:

* NIST DLMF4.37.16, <https://dlmf.nist.gov/4.37.E16>, gives the logarithmic
  inverse-hyperbolic parametrization at the branch used here.
* NIST DLMF15.4.13, <https://dlmf.nist.gov/15.4.E13>, is the exact classical
  hypergeometric representation of w with parameter `a=1/(2n)`.
  The recurrence(H13)--(H14) is proved directly above instead of importing
  that representation.

The actual statements of Long's arXiv2607.19012v1 Theorem2.1,
Proposition3.1 and Theorem4.2 were inspected:
<https://arxiv.org/html/2607.19012v1>. They supply elementary circuit
counterexamples and the SU(2) antecedent; no claim that Mathieu/Image
failure itself is new is made. The retained all-even and odd-composite
primitive families, the quadratic-height classification, and the
odd-cubic formal descent are also antecedents.

The author's published general Mathieu-property classification,
<https://hexagonmath.org/pdf/2611.00556v1>, was checked for its scope;
its quotient constructions and classification are not a prime-degree
claim. Targeted searches did not identify this exact balanced5n+3
hyperbolic family, but a search is not an exhaustive priority audit.
No worldwide novelty, external human review or formal assistant
verification is asserted.

`verify_hyperbolic_odd_height_image.py` independently checks raw
interval expansions for54(a,n,lambda) cases through m14,72 ODE/product
coefficient controls, and ORIGINAL four-variable contractions in11 cases,
including F13 through m6, the prime41/59 lifts, noncoprime a=n=3 and
the lambda=-2 control. It checks the lift's
radial and exact T factors, rejects omitted top coefficients and an
unscaled marked formula, and corrupts an angular coefficient.
`hyperbolic_odd_height_image_verification.json` records the passed
bounded tests. The written proof supplies every a,n,lambda and every m;
finite success is not used to establish those quantifiers.


The [source ledger](hyperbolic_odd_height_image_source_ledger.json) pins
the exact retained normalization and endpoint inputs and records the
classical source scopes. The [checker](verify_hyperbolic_odd_height_image.py)
and [receipt](hyperbolic_odd_height_image_verification.json) bind the
literal lift and bounded controls to these proof bytes.

```bash
python3 research/foundations/extended-geometry/verify_hyperbolic_odd_height_image.py \
  --write-receipt research/foundations/extended-geometry/hyperbolic_odd_height_image_verification.json
```
