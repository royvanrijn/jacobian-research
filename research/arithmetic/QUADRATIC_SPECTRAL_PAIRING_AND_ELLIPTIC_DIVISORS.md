# Quadratic spectral pairing and the elliptic divisor gate

The smooth affine-scale, quadratic-trace slice of the paired-critical
chart has an
elliptic spectral curve, an exact polynomial pairing criterion, and an
unavoidable accessory divisor. In degree twelve that divisor has degree
ten on the elliptic curve. The critical double zeros therefore do not
by themselves give a torsion relation for its marked point.

This note gives an implicit equivalence, not a rational parametrization
or a rational degree-twelve existence theorem. Full alternating monodromy
also requires a separate group argument. It concerns the even
numerator, odd Laurent locus
\[
 f(z)=N(z^2)/z,\qquad \deg N=2h,\qquad \deg f=4h.
\]
For degree twelve, h=3. The scale restriction is M=kappa(R-eta), with R the canonical scale.
For degree twelve, theta being quadratic in a primitive R is automatic;
the affine dependence of M on R is the additional slice condition.
The full nonreciprocal chart is larger. The singular spectral strata
are separated below.

The divisor-class bookkeeping in Section 3 is an application of classical
elliptic-curve theory. In [Milne, *Elliptic Curves* (2006)](https://www.jmilne.org/math/Books/ectext6.pdf),
Chapter I, Proposition 4.10 and the following group-law construction
(printed pp.35--36) identify a smooth genus-one curve with its degree-zero
Picard group once an origin is chosen. Here this applies over an algebraic
closure to the smooth Weierstrass curve (11), with origin O, and justifies
adding the points of the principal divisors (12)--(15). The resulting
identity (14) gives no torsion condition on B. The finite-field nontorsion
control in Section 4 uses the same book's Chapter II, Corollary 4.2
(printed p.63): prime-to-p multiplication is bijective on the reduction
kernel, hence prime-to-p torsion injects under good reduction. These
general mechanisms are classical; the spectral transformation, accessory
polynomial, and exact pairing criterion require their local derivations.

## 1. An exact trace criterion

Work over a characteristic-zero field k. Fix h>=1, L=4h-1, and
\[
 q(R)=R^h+\cdots+q_0,\quad
 \theta(R)=t_0+t_1R+t_2R^2,\quad
 \eta,\kappa\in k,\quad \kappa\ne0.
\]
Require
\[
 q_0q(\eta)\operatorname{Disc}(q)\operatorname{Disc}(S)\ne0,
 \qquad
 S(u)=\operatorname{Res}_R
       \bigl(q(R),u^2-\theta(R)u+R^2\bigr).                 \tag{1}
\]
The polynomial S is monic of degree 2h with S(0)=q_0^2. Define, uniquely,
\[
 N(u)=L\sum_{j=0}^{2h}\frac{[u^j]S}{2j-1}u^j.
 \qquad 2uN'(u)-N(u)=LS(u).                              \tag{2}
\]
In the quadratic algebra over k[R], write
\[
 N(u)=C(R)u+B_N(R)\pmod{u^2-\theta u+R^2}.
\]
Thus, for the two formal roots u,v,
\[
 C=\frac{N(u)-N(v)}{u-v},\quad T=N(u)+N(v)=\theta C+2B_N,
 \quad d=\theta-2R,\quad D=T-dC=2(B_N+RC).                \tag{3}
\]

**Theorem 1.** Subject to the opens (1), the two equations
\[
 \boxed{\quad q\mid D,\qquad
   (\theta-2R)C^2\equiv\kappa(R-\eta)\pmod q\quad}          \tag{4}
\]
are equivalent to a good paired-critical pencil whose canonical scale
is R and whose squared critical value is
\[
                         M=\kappa(R-\eta).               \tag{5}
\]
The secant element is a=C mod q and is automatically a unit. Moreover
q^2 divides D, so replacing the first condition of (4) by q^2|D is
equivalent.

Here “good” has the meaning of the
[canonical paired-critical chart](PAIRED_CRITICAL_TRACE_CHART.md):
the 2h nonzero critical u values are distinct, their nonzero squared
critical values have h distinct levels with two u values at each, and
all 4h finite critical z values are simple. Each of the 2h finite
branch values of f has exactly two critical points.

**Proof.** Over k[R]/q put
\[
 A=(k[R]/q)[u]/(u^2-\theta u+R^2).
\]
Condition Disc(S)!=0 says that all 2h u values are distinct. In
particular A is etale, u is primitive and the two u values above every
R root are different. The factors R, d and theta+2R are units:
R is a unit by q_0!=0, and
\[
                  (u-v)^2=d(\theta+2R)
\]
is a unit by separability. The first equation (4) says
N(u)=a(u-R). Consequently
\[
            N(u)^2/u=a^2(\theta-2R)=\kappa(R-\eta).
\]
The right side is a unit by q(eta)!=0. Hence a and N(u) are units.
The affine function (5) separates all q roots. Equations (1)-(2)
therefore give exactly the good paired passport. The canonical
inverse in the cited chart recovers a and R, including the sign of R.

Conversely, suppose a good paired pencil has canonical scale R,
squared value (5), and trace represented by theta. The element R is
primitive because M is, so its monic polynomial is q. Its critical
polynomial and Euler numerator are precisely (1)-(2). The canonical
identity N(u)=a(u-R) gives (4).

For the double divisibility, differentiate the globally defined trace
polynomials with respect to R at a q root. The Euler equation gives
N'(u)=N(u)/(2u) at both critical u values. Using uv=R^2 and
u+v=theta, without differentiating an assumed function a(R), gives
\[
                  T'=\tfrac12 a d',\qquad
                  C'=-\frac{a d'}{2d}.                   \tag{6}
\]
It follows that D'=T'-d'C-dC'=0 there. Since q is squarefree,
q^2|D. These calculations divide only by units already established. ∎

For h=3, fixing eta and recovering kappa from the linear coefficient,
(4) gives five coefficient equations in the six coefficients of q and
theta: three coefficients of B_N+RC modulo q, and the two conditions that
\[
 (\theta-2R)C^2\bmod q=\kappa(R-\eta)
\]
has degree one and the specified zero. Prescribing kappa as well adds
its linear-coefficient equation. Fixing eta=1 gives the
normalized affine-scale slice. In general eta is retained explicitly;
no rational square root of eta is presumed.

## 2. The accessory polynomial and the natural norm equation

Assume t_2!=0 and put m=2h. The following are exact degree statements,
not generic ones:
\[
 \deg D=2m-1=4h-1,\qquad
 \operatorname{LC}(D)=2t_2^{m-1},\qquad D(0)=-2Lq_0^2.
\]
Indeed, reduction of u^j gives a u-coefficient of degree 2j-2 and
leading coefficient t_2^{j-1}. Its constant coefficient has degree at
most 2j-2. In D=2(B_N+RC), the monic leading term of N supplies the
unique term of degree 2m-1.

Thus every good point of (4) satisfies
\[
 \boxed{\quad D=q^2U,\quad
       \deg U=2h-1,\quad U(0)=-2L,\quad
       \operatorname{LC}(U)=2t_2^{2h-1}.\quad}             \tag{7}
\]
In degree twelve the accessory U has degree five, constant -22, and
leading coefficient 2t_2^5. It cannot be dropped.

There is also an exact obstruction to treating the natural norm
identity as a constant-right-side polynomial Pell equation. Put
\[
 J=N(u)N(v),\quad
 K=J+\kappa R(R-\eta).
\]
At every q root the critical relation gives
\[
             J=-\kappa R(R-\eta),\qquad J'=J/R.
\]
Hence
\[
             q\mid K,\qquad K'(r)=\kappa r\ne0
                    \quad\text{for every }q(r)=0.        \tag{8}
\]
Each forced q root of K is simple, not double. For h>=3, J cannot be
a constant polynomial: its values modulo q equal the nonconstant
quadratic -kappa R(R-eta), which cannot become constant modulo a
monic squarefree polynomial of degree at least three. The natural
Pell-shaped identity is
\[
 T^2-(\theta^2-4R^2)C^2=4J,                              \tag{9}
\]
with this nonconstant norm. This does not exclude a different,
independently derived Pell-Abel reformulation.

## 3. The smooth elliptic spectral curve

Now assume t_0t_2!=0. Set
\[
 \alpha=t_0t_2,\quad b=t_1^2/4-1-\alpha,\quad
                 \alpha(b^2-4\alpha)\ne0.                \tag{10}
\]
The smooth projective model of
\[
                      u^2-\theta(R)u+R^2=0
\]
is the elliptic curve
\[
 E:\quad Z^2=X^3+bX^2+\alpha X,
 \qquad X=t_2u,\quad
 Z=t_2\bigl((1-t_2u)R-t_1u/2\bigr).                      \tag{11}
\]
Its distinguished points are O at infinity,
\[
                  T_2=(0,0),\qquad B=(1,t_1/2).
\]
The first point has order two. The second is a rational marked point,
with no torsion order imposed merely by the spectral curve shape.
The inverse rational functions are
\[
 u=X/t_2,\qquad
 R=\frac{Z+(t_1/2)X}{t_2(1-X)}.
\]
Direct divisors, including t_1=0, are
\[
 \boxed{\quad
 \operatorname{div}(u)=2T_2-2O,\qquad
 \operatorname{div}(R)=T_2+(B+T_2)-O-B.\quad}              \tag{12}
\]
For example B+T_2=(alpha,-t_1 alpha/2). The involution fixing u is
sigma(P)=-P; the involution fixing R is tau(P)=B-P.

Let A_q be the zero divisor of q(R), of degree 2h. It is reduced on
the good locus (1). Let Z_U be the zero divisor of U(R), with
multiplicities, of degree 4h-2. Equations (7) and (12) imply
\[
 [A_q-2hO]=h[B-O],\qquad
 [Z_U-(4h-2)O]=(2h-1)[B-O].                              \tag{13}
\]
The zero-pole divisor of D=q^2U is therefore
\[
 2A_q+Z_U-(4h-1)O-(4h-1)B.
\]
Its degree-zero class is the identity
\[
             2hB+(2h-1)B-(4h-1)B=0.                    \tag{14}
\]
For degree twelve this is 6B+5B-11B=0, not 11B=0.

There is a second useful exact identity:
\[
 S(u)=(1-t_2u)^h\,q(R)q(\sigma R),\qquad
 \operatorname{div}(S(u))=A_q+\sigma A_q-4hO.             \tag{15}
\]
The divisor classes contributed by A_q and sigma A_q cancel.
This remains valid when a zero lies at a u-ramification point, with
the corresponding divisor multiplicity.

Consequently the displayed forced divisor relations do not themselves
impose torsion on B. An additional argument would be needed to show
that the full coefficient conditions (4) impose a further torsion
relation. Neither the presence of rational two-torsion nor this
divisor count establishes such a relation. We make no claim that all
solutions of (4) have nontorsion B.

## 4. The square-root cover and its real gate

Adjoining x=sqrt(u) to E gives the smooth projective genus-one curve
\[
 C:\quad x^2+y^2=t_0-t_1xy+t_2x^2y^2,\qquad
                     R=-xy,\quad u=x^2.                 \tag{16}
\]
It is a connected unramified double cover of E. Indeed (12) has
even valuations, and u is not a square over an algebraic closure:
a square root would make T_2-O a principal divisor, contrary to
T_2 being a nonzero point of E. Riemann-Hurwitz gives genus one.

Over k the cover need not have a rational point. It is then a
genus-one covering torsor, not an elliptic curve with a chosen
rational origin. Rational torsion classification cannot simply be
applied to C.

Over the reals, on the smooth stratum (10), its solubility is exact:

* If t_0>0 or t_2>0, C has real points.
* If t_0,t_2<0, C has real points if and only if
  \[
                   |t_1|>2+2\sqrt{\alpha};
  \]
  equivalently b>0 and b^2>4alpha.

For the first assertion use x=0, or a point at infinity with
y^2=1/t_2. For the second, put c=-t_2>0 and q=xy. The least possible
value of x^2+y^2+t_1xy+c x^2y^2 at a fixed q is
2|q|+t_1q+cq^2. Minimizing over q gives the threshold
|t_1|>=2+2sqrt(alpha). Equality is singular and is excluded by (10).
No real points at infinity occur when t_2<0, so this covers the
projective curve as well.

For an actual good point of (4), there is also an exact test for
**all** finite critical z values to be real:
\[
 \boxed{\quad q\text{ splits over }\mathbb R,\quad
           \theta(r)>2|r|\text{ at every root }r\text{ of }q.\quad} \tag{17}
\]
The two u roots have product r^2>0 and sum theta(r); they are distinct
positive reals precisely under this strict inequality. Conversely,
positive critical u roots give real canonical r and the inequality.
Condition Disc(S)!=0 already separates u values from different pairs.
Real solubility of C alone is only necessary, and is much weaker
than (17).

Two exact controls illustrate the arithmetic distinctions, without
claiming new degree-twelve pencils:

1. t_0=t_2=-1,t_1=1 satisfies (10), but (16) becomes
   (x+y/2)^2+3y^2/4+x^2y^2+1=0. It has no real point.
2. t_0=1,t_1=2,t_2=2 gives
   E: Z^2=X^3-2X^2+2X and B=(1,1), which is nontorsion over Q.
   Indeed #E(F_3)=6, while B modulo 5 has order four; both primes
   have good reduction (the discriminant is -256). If B were torsion,
   its order would be divisible by four, and its prime-to-three
   torsion would inject into E(F_3), a contradiction. The required
   injection is elementary: for a good integral Weierstrass model,
   the kernel of reduction is the formal group; multiplication by
   an integer prime to p has unit linear coefficient, so it cannot
   kill a nonzero element of that kernel.

The second control shows only that the marked point of the ambient
spectral family need not be torsion. It is not asserted to satisfy (4).

## 5. Scope and exact controls

For h=3 the trace theta has degree at most two uniquely modulo q.
The affine-trace case t_2=0 is excluded over the reals by the
[uniform affine-trace theorem](UNIFORM_REAL_AFFINE_TRACE_OBSTRUCTION.md).
The strata t_0=0 or b^2=4alpha have singular spectral normalization
and are not covered by the elliptic arguments (10)-(16). They are
not declared empty here.

The accompanying checker verifies the Weierstrass transformation,
recurrences and exact degrees, the critical-jet identities, norm
derivative, constant and leading accessory terms, real threshold
identity, finite-field nontorsion control, and an actual good rational
degree-eight boundary model extended to a smooth quadratic spectral
curve. The latter checks the criterion and its converse in a case
where a solution is already known; it is not a degree-twelve solution.
The uniform statements rely on the written proof, not bounded tests.

Rational points of (4), and any additional divisor constraint capable
of controlling them, remain the actual degree-twelve arithmetic
problem. The real gate and accessory divisor prevent two unsupported
shortcuts; they do not solve that problem.

