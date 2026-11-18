# A uniform real obstruction to affine critical trace

For every degree 4h>=12, a good real paired odd Laurent pencil cannot
have its critical-pair trace affine in its canonical scale. The secant
element may vary, and the critical points need not be real.

This excludes the affine-trace boundary of the spectral curve for
nonreciprocal pencils. It also gives a uniform real obstruction to a
scalar secant with affine scale, strengthening the
[degree-twelve rational arithmetic obstruction](DEGREE_TWELVE_SCALAR_AFFINE_OBSTRUCTION.md).
It does not exclude every singular spectral cubic or every
nonreciprocal model.

## Theorem

Let h>=3, m=2h, L=2m-1, and let N in R[u] be monic of degree m.
Assume R_N=2uN'-N has m distinct nonzero roots, and M=N(u)^2/u
takes h distinct nonzero values on them, each exactly twice.
For each pair u,v write nu=N(u), nu'=N(v), and define

\[
                   \rho=-\nu\nu'/M,\qquad \theta=u+v.
\tag{1}
\]

These are the canonical scale and trace of the
[paired critical chart](PAIRED_CRITICAL_TRACE_CHART.md).

It is impossible that

\[
                         \theta=s+t\rho,\qquad s,t\in\mathbb R.
\tag{2}
\]

All coefficients and scalars in the statement are real. No assumption
on the locations of the critical points is made. A target scaling
removes the monic restriction without changing the condition (2).

### Corollary: scalar secant and affine scale

In every degree 4h>=12, no good real model has both a scalar secant
a and an affine scale rho=r_0+cM with real scalar coefficients.
Thus every rational such model is excluded as well.

Indeed, the canonical identities are

\[
       N(u)=a(u-\rho),\qquad
                    \theta=2\rho+M/a^2.
\tag{3}
\]

If a is scalar and c!=0, then

\[
             \theta=\left(2+\frac1{a^2c}\right)\rho
                                  -\frac{r_0}{a^2c},
\]

contradicting the theorem. The case c=0 is the previously proved
both-scalar obstruction, whose short proof is recalled below.
The unrestricted nonconstant-secant, nonlinear-trace locus is not
excluded.

## 1. The scale separates the pairs

The canonical definitions give

\[
                         uv=\rho^2,\qquad
                    N(v)=-\rho N(u)/u.                        \tag{4}
\]

All u,v,N(u),M,rho are nonzero. The secant
a=(N(u)-N(v))/(u-v) is nonzero: otherwise equality of the nonzero
N-values and their squared critical values would force u=v.
Equation (3) follows directly from (4), as does

\[
                            d:=\theta-2\rho=M/a^2\ne0.
\]

If two pairs had the same rho, (2) would give them the same theta.
Their u-coordinates would be the same two roots of
u^2-theta*u+rho^2, contradicting distinctness of all critical roots.
Thus rho separates the h pairs.

Let q(R) be the monic degree-h polynomial on these scale values.
It has distinct roots and q(0)!=0. Conjugation preserves the
canonical pairing, so q has real coefficients. The monic critical
polynomial S=R_N/L is the norm of

\[
                     u^2-(s+tR)u+R^2
\]

over R[R]/q. Consequently

\[
                   S(0)=q(0)^2,\qquad
                         N(0)=-Lq(0)^2.                      \tag{5}
\]

## 2. A polynomial square forced by the critical equations

For formal roots u,v of that quadratic, define polynomials in R

\[
 C(R)=\frac{N(u)-N(v)}{u-v},\quad T(R)=N(u)+N(v),\quad
 D(R)=T(R)-d(R)C(R),\quad d=s+(t-2)R .
\tag{6}
\]

Their polynomial definitions extend across coincident formal roots.
In particular deg C<=m-1 and deg D<=m.

At a scale value w, put a=C(w). The critical equations give
N'(u)=N(u)/(2u), and similarly for v. Differentiating locally in
R, where the actual pair roots are distinct, gives

\[
                   T'(w)=\frac{a(t-2)}2,\qquad
                 C'(w)=-\frac{a(t-2)}{2d(w)}.                 \tag{7}
\]

Here no derivative of a is introduced: only its value at the
critical pair is substituted. No inverse of u' or v' is used.
For example, writing delta=u-v, one uses uv=R^2 and

\[
 (\theta+2R)(\delta^2)'-2t\delta^2
                            =4R(\theta+2R)(t-2).
\]

Equations (3), (6) and (7) give D(w)=D'(w)=0.
Since q is squarefree, q^2 divides D. At R=0 the polynomial trace
formula gives D(0)=2N(0), including s=0. From (5),

\[
                            D=-2Lq^2.                        \tag{8}
\]

This also proves deg D=m, so there is no unhandled leading-degree drop.

Let gamma_0=0, gamma_1=1 and gamma_j=t*gamma_(j-1)-gamma_(j-2).
For N(u)=u^j the leading coefficient of D is
2(gamma_j-gamma_(j-1)); this follows from the quadratic trace
recurrences. The leading coefficient of (8) therefore says

\[
                     \gamma_m-\gamma_{m-1}+L=0.              \tag{9}
\]

## 3. Real positivity forces the slope minus two

Let T_0(x)=2, T_1(x)=x and T_(j+1)=x*T_j-T_(j-1), and put

\[
 B_h(x)=\frac{1+\sum_{j=1}^{m-1}T_j(x)-(2m-1)}{x-2}.
\tag{10}
\]

It is an integer polynomial. The trace recurrence, telescoped using
T_j(t)=gamma_(j+1)-gamma_(j-1), gives for even m

\[
             \gamma_m-\gamma_{m-1}+L=(t+2)B_h(-t).            \tag{11}
\]

We include the positivity argument for this witness, also used in the
[earlier scalar-secant note](SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md).

For x>2 choose a positive r!=1 with x=r+r^-1.
The numerator in (10) equals r^(-(m-1)) times

\[
               \sum_{j=0}^{2m-2}r^j-(2m-1)r^{m-1}>0
\]

by strict arithmetic-geometric mean. The denominator is positive.
For x<-2 choose r<0. The displayed sum is positive: its geometric
sum is positive and its subtracted central term is negative.
Since m-1 is odd, the numerator in (10) is negative, as is x-2.
For -2<=x<2 choose |r|=1,r!=1. Then
1+sum T_j(x)=1+2*sum Re(r^j)<2m-1; again both signs are negative.
Finally

\[
                         B_h(2)=\sum_{j=1}^{m-1}j^2>0.
\]

Thus B_h is strictly positive on R. Equations (9)-(11) force t=-2.

If s=0, the critical quadratic is (u+rho)^2, contradicting distinct
critical roots. If s!=0, rescale u,rho by s and replace N(u) by
s^-m N(su). This uses real coefficients even when s<0, and preserves
all good conditions. Hence we may assume

\[
                            \theta=1-2\rho .                 \tag{12}
\]

## 4. A reflection identity and a leading-coefficient contradiction

Write f(z)=N(z^2)/z and

\[
                            Q(z)=q(z(1-z)).
\]

The critical quadratic is (u+rho)^2=u. For each critical u select
z=u+rho, so z^2=u, rho=z(1-z), and the paired root is (1-z)^2.
The m roots of Q are distinct; none equals 0 or 1.
The critical norm factorization gives

\[
 S(z^2)=Q(z)Q(-z),\qquad
                    f'(z)=L\,Q(z)Q(-z)/z^2.                 \tag{13}
\]

Equation (4) gives f(1-z)=-f(z) at each selected critical root.
Both derivatives vanish there. Therefore the polynomial

\[
              (1-z)N(z^2)+zN((1-z)^2)
\]

has m distinct double zeros. Its degree is 2m and its leading
coefficient is -L. Since Q^2 has degree 2m and leading coefficient
one, this proves the rational identity

\[
                 f(z)+f(1-z)=-L\frac{Q(z)^2}{z(1-z)}.         \tag{14}
\]

Differentiate (14), substitute (13), and cancel the nonzero
polynomial Q. Since Q(-z)=Q(z+1), the result is

\[
 \mathcal LQ:=(1-z)^2Q(z+1)-z^2Q(z-1)
          +2z(1-z)Q'(z)-(1-2z)Q(z)=0.                       \tag{15}
\]

For any polynomial of degree n, this operator has degree at most n-1.
If its leading coefficient is b, its coefficient at z^(n-1) is

\[
                           b\,\frac{n(n-2)(n-4)}3.           \tag{16}
\]

For clarity, Taylor expansion shows that the only contributions
at that degree are Q', -zQ'', and z^2Q'''/3. All higher derivative
terms have degree at most n-3. Formula (16) follows directly.

Here n=deg Q=m>=6 and b=(-1)^h. The coefficient in (16) is nonzero,
contradicting (15). This proves the theorem.

For degree twelve the contradiction is especially short. If
Q=q_0+q_1*z(1-z)+q_2*z^2(1-z)^2+z^3(1-z)^3, then

\[
                    \mathcal LQ=8(1-2z)z^2(1-z)^2\ne0.
\tag{17}
\]

The degree-eight threshold is sharp for the trace hypothesis:
the known rational degree-eight model has theta=-1/6-2rho.
Its normalized Q has degree four, and (16) vanishes there.

## 5. The both-scalar endpoint in the corollary

For completeness, suppose a and rho are nonzero base scalars.
At all m critical roots, N-a(u-rho)=0, hence

\[
 S=N-au+a\rho
      =u^m+\frac{a}{2m-2}u+\frac{a\rho}{2m}.
\]

The paired norm gives S(0)=rho^m, so a=2m*rho^(m-1).
Reciprocity under u->rho^2/u then forces the coefficient of u^(m-1)
to be m*rho/(m-1), whereas it is zero when m>=4.
This is the both-scalar argument of the earlier note and proves
the c=0 case needed here.

## Evidence and scope

The [checker](verify_uniform_real_affine_trace_obstruction.py) verifies
the critical-jet polynomial identity, the slope/witness identity and
the operator coefficient formula in bounded degrees. It also verifies
the exact reflection identities on the rational degree-eight boundary
control. The uniform proof is the written argument above.

The [receipt](uniform_real_affine_trace_verification.json) and
[source ledger](uniform_real_affine_trace_source_ledger.json) retain
the actual exact execution and completed internal reviews.
No claim is made for odd numerator degree, arbitrary nonaffine trace,
every singular spectral cubic, rational degree-twenty-six existence,
novelty, formal verification or external review.

