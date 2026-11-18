# A rational degree-twelve obstruction for scalar secants and affine scales

No good rational degree-twelve pencil in the paired critical trace chart
has a scalar secant element and an affine reciprocal scale. This includes
nonreciprocal pairings and nonreal critical points. It completes the
degree-twelve affine-scale slice left open by the
[uniform endpoint obstructions](SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md);
the unrestricted trace chart remains open.

The proof has a uniform algebraic reduction and a finite exact certificate
in degree twelve. The certificate consists of arithmetic identities,
not a numerical search or a local-good-reduction argument. No rational
denominator restriction is imposed.

## Theorem

Let N in Q[u] be monic of degree six. Assume R_N=2uN'-N has six
distinct nonzero roots, and that M=N(u)^2/u takes three distinct
nonzero values on them, each exactly twice. In the
[paired critical algebra](PAIRED_CRITICAL_TRACE_CHART.md), let a,rho
be the canonical units with

\[
                           N(u)=a(u-\rho).
\]

It is impossible that a is a base scalar and rho is affine in M:

\[
              a\in\mathbb Q^*,\qquad \rho=r_0+c_{\rm old}M,
              \qquad r_0,c_{\rm old}\in\mathbb Q.                 \tag{1}
\]

Consequently no all-real rational degree-twelve model lies in this
slice. The theorem does not exclude a nonconstant a or a general
quadratic rho in the three-dimensional paired-value algebra.

Dividing a numerator by its nonzero leading coefficient scales a and M,
and preserves the affine form of rho. Thus the monic normalization
does not cost a rational extension.

The statement uses the coefficient field's scalars. If the canonical
elements of a rational good model become scalar/affine after extension
to C, their coefficients in (1) are already rational: the canonical
elements are defined over Q, and 1,M are linearly independent, making
the affine coefficients unique.

## 1. A uniform rational normalization

First suppose a*r_0*c_old is nonzero. No extraction of roots is needed.
Set

\[
 U=u/r_0,\quad P(U)=N(r_0U)/(ar_0),\quad
 W=M/(a^2r_0),\quad c=a^2c_{\rm old}.
\tag{2}
\]

Then the normalized secant is one and the normalized scale is

\[
                         \rho=1+cW.
\]

Here and below u,M denote these normalized U,W coordinates. The new
polynomial P need not be monic. If N has degree m=2h, its leading
coefficient becomes lambda=r_0^(m-1)/a, a nonzero rational number.
The Euler critical equation scales by 1/(ar_0), so all good conditions
and pairing identities survive. This normalization works over the
original characteristic-zero field.

Conversely, any rational good normalized P with leading coefficient
lambda gives a monic rational numerator N=P/lambda: take a=1/lambda,
r_0=1 and c_old=c*lambda^2. Thus the nonmonic normalized problem loses
no rational models.

In a pair of critical coordinates u,v with common value M, write

\[
 \rho=1+cM,\quad s=u+v=2+(2c+1)M,\quad uv=\rho^2,\quad
                 \delta^2=(u-v)^2=M(M+4\rho).                 \tag{3}
\]

The units and distinct-root hypotheses of the trace chart justify
these identities and division by u,v,rho,M,delta at every pair.

Equivalently, x=u/P(u) parametrizes the critical conic:

\[
 u=\frac{x^2}{x^2-x-c},\quad
 P(u)=\frac{x}{x^2-x-c},\quad
 M=\frac1{x^2-x-c},\qquad x\longleftrightarrow1-x.
\tag{4}
\]

Equation (4) parametrizes the conic; the polynomial and critical
conditions are additional constraints.

## 2. A uniform square and Hermite-interpolation condition

The following derivation applies in every degree m=2h. On the quadratic
algebra given by (3), form polynomials in M

\[
 C(M)=\frac{P(u)-P(v)}{u-v},\quad
 T(M)=P(u)+P(v),\quad D(M)=T(M)-MC(M).                         \tag{5}
\]

They have degrees at most m-1,m,m respectively. Let q be the monic
degree-h polynomial whose roots are the distinct paired values.
At each root w of q, the secant equation gives

\[
                 C(w)=1,\qquad T(w)=w.                       \tag{6}
\]

The critical equations give P'(u)=P(u)/(2u) and its v analogue.
Differentiate locally in M; the two roots u,v vary analytically
because delta is nonzero. Since uv=rho^2,

\[
 T'(w)=\frac12\left(s'-\rho\left(\frac{u'}u+
                                      \frac{v'}v\right)\right)
                         =\frac12.                          \tag{7}
\]

For E=P(u)-P(v)-(u-v), its derivative at the pair is

\[
 E'=-\frac12\left(\delta'
             +\rho\left(\frac{u'}u-\frac{v'}v\right)\right)
                         =-\frac{M+4\rho}{2\delta}.
\]

Indeed, delta^2=M(M+4rho), s'=2c+1 and

\[
 (M+4\rho)(\delta^2)'-2s'\delta^2
                                      =4\rho(M+4\rho).
\]

Since E=delta(C-1), this proves

\[
                    C'(w)=-\frac1{2w}.                      \tag{8}
\]

No inverse of u' or v' is used, so stationary projections of the conic
cause no exception. Equations (6)-(8) imply D(w)=D'(w)=0.
The q-roots are distinct and deg D<=2h, hence

\[
 D=\kappa q^2,\qquad
 q\mid C-1,\qquad q\mid C+2MC'.                              \tag{9}
\]

The scalar kappa is nonzero. At M=-1/c, rho=0 and the quadratic
roots are 0,M, so (5) gives

\[
                         D(-1/c)=2P(0)\ne0.                 \tag{10}
\]

Nonvanishing of P(0) follows from the nonzero critical-root hypothesis:
the Euler polynomial has constant coefficient -P(0).

This is a necessary reduction, not an assertion that arbitrary solutions
of (9) satisfy all critical and good-open conditions.

## 3. The generic degree-twelve equations

For degree twelve, m=6 and h=3. Write

\[
                q=M^3+q_2M^2+q_1M+q_0.
\]

The basis polynomials in (5) are determined without root choices:

\[
 C_0=0,\ C_1=1,\quad T_0=2,\ T_1=s,\qquad
 C_j=sC_{j-1}-\rho^2C_{j-2},\
 T_j=sT_{j-1}-\rho^2T_{j-2},\quad D_j=T_j-MC_j.
\tag{11}
\]

For P=sum n_j u^j, these give C=sum n_j C_j and D=sum n_j D_j.
The degree-j coefficients of D_j are

\[
\begin{array}{c|l}
 j&\delta_j(c)\\ \hline
0&2\\
1&2c\\
2&2c(c+1)\\
3&2c(c^2+3c+1)\\
4&2c(c^3+6c^2+5c+1)\\
5&2c(c+1)(c^3+9c^2+6c+1)\\
6&2c(c^5+15c^4+35c^3+28c^2+9c+1).
\end{array}                                                   \tag{12}
\]

The only rational zeros among these coefficients are 0 and -1:
the remaining monic factors have constant one and no root at 1 or -1,
so the rational-root theorem suffices.

When their product is nonzero, the triangular system
sum n_j D_j=q^2 uniquely determines a polynomial P_0 and its secant
polynomial C_0^*. Since kappa!=0, (9) for the actual P means
P=kappa P_0 and C=kappa C_0^*. Necessary conditions are therefore

\[
 \begin{split}
 [M]\operatorname{rem}_q C_0^*&=0,\qquad
 [M^2]\operatorname{rem}_q C_0^*=0,\\
 \operatorname{rem}_q(C_0^*+2M(C_0^*)')&=0.
 \end{split}                                                 \tag{13}
\]

These are five equations in q_0,q_1,q_2,c, after clearing only
c-dependent denominators. Their complete exact coefficients are retained
in the [reduced system](degree12_scalar_affine_reduced_system.json).
The good locus also has q_0*Disc(q)!=0. Adjoin z with

\[
                         zq_0\operatorname{Disc}(q)-1=0.     \tag{14}
\]

## 4. An explicit arithmetic identity and its specializations

The [certificate bundle](degree12_scalar_affine_arithmetic_certificates.json)
contains a 42-node arithmetic circuit over Q(c) expressing 1 in the ideal
generated by the five cleared equations (13) and (14).
Every node is an input, zero, a scalar multiple, a difference of two
monomial multiples, or an explicitly stored polynomial combination of
previous nodes. The checker reconstructs the final value 1 using exact
addition and multiplication; it invokes no Groebner algorithm.

All 54 distinct denominators of the circuit's 102 stored field
coefficients have rational zeros only in

\[
                       \mathcal E=\{0,-1/4,-1/3,-1/2,-1\}.   \tag{15}
\]

This root-set assertion has its own
[pure arithmetic certificate](degree12_scalar_affine_denominator_certificates.json).
For every circuit denominator and every triangular coefficient (12),
the checker multiplies the supplied factors and verifies the exact
polynomial product, including multiplicities. Linear factors give
precisely the rational exceptions in (15). Each of the nineteen
remaining nonlinear factors has a prime p<=41 such that its integer
leading coefficient is nonzero modulo p and its value at every point
of F_p is nonzero. All these finite evaluations are replayed.

Here is the elementary denominator argument. If a primitive integer
polynomial f has these two properties and a/b is a rational root in
lowest terms, then p cannot divide b: after multiplying by b^deg(f),
the leading term would be nonzero modulo p and all other terms zero.
When p does not divide b, reduction gives a forbidden root in F_p.
Thus f has no rational root, with arbitrary rational denominators
covered. The factor-product and finite-evaluation verification uses no
factorization algorithm or irreducibility claim.

For rational c outside E, every circuit coefficient specializes, as do
the triangular equations; its final identity is still 1. Thus there
is no solution there. This specialization argument has no p-adic
integrality or good-reduction premise.

At each of c=-1/4,-1/3,-1/2,-1 use the direct trace chart instead of
triangular inversion. Put

\[
 S(u)=\operatorname{Res}_M
   (q(M),\,u^2-(2+(2c+1)M)u+(1+cM)^2),
\]
\[
 N_S(u)=11\sum_{j=0}^6\frac{[u^j]S}{2j-1}u^j .
\tag{16}
\]

For an actual good model S is its monic critical polynomial and
P=lambda N_S, by the Euler equation. The necessary six equations
are the two quadratic-remainder coefficients, each reduced modulo q,
of

\[
                   \lambda N_S(u)-(u-1-cM)=0.               \tag{17}
\]

The necessary unit product is

\[
 \lambda q_0\operatorname{Disc}(q)
       \operatorname{Norm}_{B/\mathbb Q}(1+cM)
       \operatorname{Norm}_{B/\mathbb Q}(M(M+4+4cM)),
 \qquad B=\mathbb Q[M]/q .
\tag{18}
\]

Adjoin an inverse for (18). For each of the four fixed parameters,
the bundle contains an arithmetic circuit expressing 1 in this
direct necessary ideal. Their respective node counts are
148,94,79,94 for -1/4,-1/3,-1/2,-1. All coefficients are rational
numbers, so no further parameter exceptions occur.
The checker independently reconstructs (16)-(18), binds the circuit
inputs, and replays all four identities.

This excludes every nonzero rational c in the case a*r_0*c_old!=0.
If r_0=0 and c_old!=0, the frozen uniform real proportional-scale
theorem applies to rational coefficients. If c_old=0, the frozen
both-scalar proposition applies over C. If both vanish, rho would
not be a unit. All cases in (1) are now excluded. This proves the theorem.

## Evidence and limits

Run

~~~sh
.venv/bin/python research/arithmetic/verify_degree_twelve_scalar_affine_obstruction.py
~~~

The [checker](verify_degree_twelve_scalar_affine_obstruction.py) rederives
the trace basis, triangular inverse and five generic equations, checks
the rational exceptional parameters by exact factor products and
projective root-free prime witnesses, reconstructs the four direct
systems, and replays every arithmetic circuit. The
[receipt](degree_twelve_scalar_affine_verification.json) records the
actual result and source hashes; the
[source ledger](degree_twelve_scalar_affine_source_ledger.json) records
proof imports and completed internal reviews.

No unrestricted rational degree-twelve nonreciprocal obstruction or
uniform all-degree affine-scale classification follows. Real-algebraic
models from Hurwitz theory may occupy other chart loci. No rational
degree-twenty-six construction, novelty claim, formal verification,
independent external replay or external human review is asserted.

