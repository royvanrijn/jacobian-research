# Complete z-dependent rational sweep chart

Work over algebraically closed characteristic-zero constants.
Consider the polynomial tautological sweep

S(gamma,w,z)=X(w,z)+gamma(1,w,z), det JS=gamma^2,

whose nonzero normalization matrix is
lambda*[[-rho,1],[-rho^2,rho]], with lambda polynomial and
rho in k(w)[z] depending on z. The full geometric group of
this entire chart is determined below. The assertion is about
these field pairs, not every Keller map or every rational multiplier.

**KMFI41 (complete integrability parameters).**
Up to translations/scaling, choose integers r>=2,s>=0 and a root

\[
 C(r,s)=r(r-1)+2\kappa(r-1)(s+1)
                 +\kappa^2(s+1)(s+2)=0.
\]

Put rho=(kappa*z+R(w))/w. Here R is any Laurent polynomial.
If its lowest exponent m is negative, the exact admissibility
condition is

\[
 r+(s+2)m\ge2. \tag{1}
\]

There is no further condition when R is polynomial.
For each such parameter choice there is exactly one polynomial
lambda with leading z coefficient w^r, satisfying the wave equation
and polynomiality of the normalization matrix.
It is given by the finite coefficient recursion below.
The associated polynomial X is unique up to target translation.
Conversely every z-dependent multiplier in this chart is obtained
in this way. A scalar leading coefficient is restored by scaling.

**KMFI42 (degree and complete affine boundary).**
Remove the affine part of R by the KMFI40 coordinate isomorphism.
If no term remains, the homogeneous KMFI34 case applies.
Otherwise let j be the largest remaining exponent. The generic
degree is exactly

\[
 \boxed{N=(s+2)\bigl(r+(s+1)j\bigr).} \tag{2}
\]

It is composite, including when j is negative.
After centering X(0,z)=0, the sweep is finite over Y2!=0.
All omitted affine normalization branches lie over the plane
E={Y2=0}. At its generic point put
g=gcd(s+1,r-1), e=(s+1)/g.
Its full inertia is

\[
 e^g\,1^{\,N-(s+1)}. \tag{3}
\]

The ramification image D=X(A2) is distinct from E, is birationally
parametrized by X, and has full inertia a single three-cycle.
Thus the two additional ramified branches over D required by the
previous exceptional PSL2(11) gate cannot occur in this chart.

**KMFI43 (full geometric groups).**
For non-affine R, the group in its actual degree-N action is

\[
 \boxed{G=S_N\ \text{if r is even and s is odd; otherwise }G=A_N.}
 \tag{4}
\]

For affine R it is natural A_(r+s+1) in odd degree or S_(r+s+1)
in even degree, by KMFI34.
Thus neither specified PSL2(11) action occurs through this entire
z-dependent rational chart, with any allowed positive z degree
in lambda or any allowed higher pole, whether S is finite or not.
This is a field-pair exclusion, hence survives independent stabilization
and arbitrary birational recharts of the same pair.

The full arithmetic constants of the new positive-z-degree families
are not determined here. When the geometric group is S_N the
arithmetic group is already S_N and the constants are unchanged;
when it is A_N, an arithmetic enlargement is at most quadratic.
KMFI35--37 retain exact constants in their earlier s=0 slice.

[Formula module](rational_wave.py) · [exact checker](verify_rational_chart.py) ·
[model](rational_chart_model.json) · [receipt](rational_chart_receipt.json).

The full inverse problem, absolute PSL2(11) occurrence and other
direction fields remain open. The polynomial-multiplier case is
now completed by [KMFI44--46](POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md).
No novelty priority, formal verification or independent whole-proof
review is asserted.

## 1. The triangular wave recursion

KMFI33 gives the displayed single-pole form, top coefficient w^r,
and the residue equation. Its discriminant is negative rational,
so kappa is nonreal and is not an integer.
Write lambda=sum_(b=0)^s f_b(w)z^b, f_s=w^r.
Let C(a,b) denote the same quadratic with r,s replaced by a,b.
For b<s, define

\[
 Q_b=2(b+1)w^2(Rf_{b+1}/w)'
       +2\kappa(b+1)(b+2)Rf_{b+1}
       +(b+1)(b+2)R^2f_{b+2}. \tag{5}
\]

If Q_b=sum_a q_a w^a, then

\[
 f_b=-\sum_a {q_a\over C(a,b)}w^a. \tag{6}
\]

This is just the z^b coefficient of
lambda_ww+2(rho*lambda)_wz+(rho^2*lambda)_zz=0 after multiplication
by w^2. For every nonnegative integer a and b<s, C(a,b)!=0:
at a=0,1 this follows from nonreal kappa; for a>=2 it follows
from the unique integer-pair argument of KMFI34.
Thus the polynomial solution, if it exists, is unique.

If m and j are respectively the least and greatest exponents of R,
the exponents in f_b lie between r+(s-b)m and r+(s-b)j.
This follows by induction in (5). Terms outside the nonnegative
range must vanish if a polynomial solution exists.

For negative m, take the lowest m-weight component, with
weight(w)=1,weight(z)=m. It has the form
w^(r+s*m)*L(z/w^m), deg L=s, leading coefficient one.
Its constant coefficient is nonzero.
Indeed if its least z power b0 were at least2, the R_m^2 term
of the wave equation at power b0-2 could not cancel.
If b0=1, its power-zero cross term has nonzero coefficient
2(a+m-1+2kappa), again impossible.
Hence the lowest w exponent in f0 is exactly r+s*m.
Polynomiality of rho^2*lambda requires this exponent to be
at least2(1-m), giving (1).

Conversely (1) makes every exponent generated in (6) at least
r+s*m>=2(1-m). Thus all f_b are polynomial and the squared
denominator of rho divides lambda. Every division in (6) is
nonzero and the recursion solves the wave equation.
For polynomial R the same argument uses the bound r>=2.
This proves the exact necessity and sufficiency in KMFI41.

## 2. Reconstructing the polynomial surface

Put A=-(lambda_w+(rho*lambda)_z) and
B=(rho*lambda)_w+(rho^2*lambda)_z.
They are polynomials. The wave equation is exactly A_w=B_z,
so they integrate to polynomial X1 with X1_z=A,X1_w=B.
Then integrate

X2_z=wX1_z+lambda, X2_w=wX1_w-rho*lambda,

X3_z=zX1_z+rho*lambda, X3_w=zX1_w-rho^2*lambda.

Their mixed partials agree by the displayed X1 identities.
All entries are polynomial, so polynomial primitives exist in
characteristic zero. Direct substitution gives the required
nilpotent normalization matrix and det JS=gamma^2.
The kernel of N JH=0 on polynomial three-vectors consists only
of constants, proving uniqueness up to target translation.

Center those constants to make X(0,z)=0, using KMFI40.
Set V=X2-wX1, W=X3-zX1. Eliminating gamma gives the exact fiber:

\[
 A_Y=V+wY_1-Y_2=0,\qquad B_Y=W+zY_1-Y_3=0. \tag{7}
\]

Their z degrees are s+1 and s+2, with leading coefficients

\[
 {w^r\over s+1},\qquad
 b\,w^{r-1},\quad
 b={r+2\kappa(s+1)\over(s+1)(s+2)}\ne0.
 \tag{8}
\]

The same letter b in (8) is a scalar, not the recursion index.

## 3. Nonzero leading resultant: the weighted model

After removing affine R, j>=2 or j<=-1.
The negative-pole bound implies
r+(s+2)j>=2 when j<0.
Put M=r+(s+1)j, H=r+(s+2)j-1.
Then M>=3 and H>=1.
The highest j-weight components of V,W have weights M,H.

A nonzero coefficient of R_j can be absorbed by scaling z
and the remaining scalar parameters for the following root argument.
Normalize R_j=1.
Write the corresponding leading lambda as
w^n L(t), n=r+s*j, t=z/w^j, with L monic of degree s.
Its wave equation is

\[
 A(t)^2L''+B(t)L'+C(n,0)L=0,\quad A=1+(\kappa-j)t.
\]

At t0=-1/(kappa-j),

\[
 B(t_0)=-{j(j-1)\over\kappa-j}\ne0.
\]

A zero of L at t0 would give an uncanceled lowest-order term
B(t0)L'. At any other point the ordinary differential equation
forbids a repeated zero. Thus L has simple roots and L(t0)!=0.

Let a=M-1 be the weight of X1 in this leading model.
Euler differentiation of the X1 identities gives

\[
 aX_1=A^2L'+E L,\quad
 E=nA+(j-1)+2\kappa+\kappa(2\kappa-j-1)t,
\]
\[
 E(t_0)=-{j(j-1)\over\kappa-j}\ne0.
\]

The leading surface has no simultaneous zero at w!=0.
If L=0, simplicity and A!=0 make X1!=0.
If L!=0 and X=0 along its weighted ray, applying N JX to
the ray direction forces A=0; then E(t0)L!=0 again makes X1!=0.

The leading V and W cannot have a common root t.
Otherwise S(-X1,w,t*w^j)=0 along a nonconstant curve.
Its derivative makes det JS vanish, so gamma=-X1=0 and then
X=0, contrary to the preceding paragraph.
Thus the leading weighted resultant is nonzero.

## 4. Exact resultant degree and finiteness away from one plane

Let Rcal(w)=Res_z(A_Y,B_Y). The weighted resultant rule gives

\[
 \deg_w Rcal=(2s+3)r+(s+1)(s+2)j-(s+1),
 \tag{9}
\]

with nonzero constant leading coefficient independent of Y.
Indeed substitute z=w^j*t; the resultant scaling is
w^(j(s+1)(s+2)), and the two equation weights are M,H.

At w=0, take generic Y2!=0. The s+1 roots of A_Y satisfy

z^(s+1)~(s+1)Y2/w^r.

All lower coefficient terms are smaller. For a negative minimum
exponent m, this follows from m+r/(s+1)>0, a consequence of (1);
for polynomial R it is immediate.
At these roots the leading term of B_Y is
b*w^(r-1)*z^(s+2), which dominates Y1*z and every lower term.
The product formula for the resultant gives

\[
 Rcal=w^{(s+1)(r-1)}Q(w),\qquad
 Q(0)=(-1)^{s+2}b^{s+1}Y_2^{s+2}. \tag{10}
\]

There is no further common w factor.
The leading coefficient of Q is constant and nonzero, and
(9)--(10) give deg Q=N from (2).
For w!=0 the two z-leading coefficients in (8) are units.
Thus no other resultant factor is extraneous: its roots count
the full generic fiber with multiplicity.
Generic fibers are reduced because det JS is not identically zero,
so their total length is the field degree N.
The projection to w need not be primitive; resultant multiplicities
still count every common z solution. This proves (2).

On Y2!=0, (10) makes w integral and invertible.
Equation A_Y then makes z integral; gamma=Y1-X1 is integral too.
Consequently S is finite over that open set.
The finite normalization has no missing affine branch outside E.

For j>=2, M>=2(s+2), so N has two factors greater than one.
For j<=-1, M>=2-j>=3, again proving compositeness.

## 5. The complete generic normalization above E

At a generic point of E, Y1,Y3 are nonzero transcendental residues.
All source points with w!=0 are present and unramified, since
the ramification image D is a different divisor.
There is also the present unramified point
w=0,gamma=Y1,z=Y3/Y1.

Every omitted branch has w=0,z=infinity.
If a negative power m occurs in R, a branch with z=O(w^m)
would have X tending to zero by (1), so gamma tends to Y1 and
the third target coordinate diverges. Such a branch is impossible.
Hence z grows faster than R and only the top z terms matter.

Put T=w^(r-1)z^(s+1), beta0=1/(s+1), and use scalar b from (8).
The two equations in (7) give

\[
 Y_2=w(Y_1+\beta_0T)+\text{smaller terms},\qquad
 Y_3=z(Y_1+bT)+\text{smaller terms}.
\]

They force T to tend to -Y1/b.
The ratio Y2/w tends to Y1*(1-beta0/b), which is nonzero
because b!=beta0. Thus w has the same order as Y2.
The nondegenerate binomial local equation is

\[
 w^{r-1}(1/z)^{-(s+1)}=-Y_1/b+\text{higher terms}.
\]

Newton--Hensel lifting gives exactly g=gcd(r-1,s+1) branches,
each of index e=(s+1)/g over E.
Their total sheet count is s+1.
All other sheets are unramified, proving (3).
The inequality in (1) makes every discarded Laurent term
strictly higher order in these expansions.

KMFI40 proves X birational onto its surface D.
Also X2_z has nonzero top coefficient
(1-r-kappa(s+1))*w^r*z^s, so D is not the plane E.
Finiteness on Y2!=0 and the Jacobian identity therefore make
the *full* inertia over D a single three-cycle.
There are no omitted branches above D. This closes the boundary
coincidence left open by KMFI40.

## 6. Full group, not just compatible local permutations

Let G be the monodromy group and N0 the normal subgroup generated
by the conjugates of this single three-cycle.
The Galois quotient by N0 is unramified off E:
all other finite branch inertia is killed.
Purity makes it etale over A3 minus E=G_m times A2, whose finite
geometric fundamental-group quotients are cyclic, generated
by the E-inertia image.

The N0-orbits form a transitive cyclic quotient action.
E-inertia has a fixed sheet (the present unramified point above),
so it fixes an N0-orbit. A cyclic transitive action generated
by an element with a fixed orbit has just one orbit.
Thus N0 is transitive.
The three-cycle support lemma makes N0=A_N.
Hence A_N<=G<=S_N.

The sign of (3) is (-1)^((s+1)-g).
It is odd exactly when r is even and s is odd.
If it is odd, G=S_N.
Otherwise every codimension-one inertia is even, and purity and
affine simple connectedness make the sign cover trivial, so G=A_N.
This proves KMFI43 for all parameters, not a bounded group census.

The excluded PSL2(11) groups remain open through other field pairs:
genuinely two-variable rational charts outside k(w)[z], other
direction fields and higher-dimensional constructions. KMFI44
proves every nonconstant polynomial multiplier in this chart finite.
