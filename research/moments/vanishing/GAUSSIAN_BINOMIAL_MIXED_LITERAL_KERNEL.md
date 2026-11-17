# Exact original kernels and literal images for W(1+Z)^m Y

Let K be any characteristic-zero field. Use one Witt pair Z,W with
E[Z^a W^b]=delta_(a,b)a! and an independent ordinary Gaussian Y. For an
integer m>=1 put

    P_m=W(1+Z)^m Y, Phi_m[A]=E[A exp(tP_m)],
    V_m=Phi_m[K[Z,W,Y]], M_m=K[t] V_m.

Here V_m uses fixed original polynomial markers; M_m also allows
polynomial time coefficients. Write V_(m,d) for the image of original
markers of total spatial degree at most d. All image and kernel claims
below concern this actual Wick functional, not an alternate trace or an
upper bound from a differential equation.

For every m>=2 the complete original kernel is

    ker Phi_m=delta_1(K[Z,W,Y])+delta_2(K[Z,W,Y]),       (1)

where

    delta_1 A=(ZW-Y^2)A-W A_W+Y A_Y,
    delta_2 A=[W-(m-1)ZW+(m-1)]A
                            -(1+Z)A_Z+mW A_W.        (2)

The right side is a sum of K-vector spaces of operator images. It is
not asserted to be an ideal. Every marker has a unique image represented
by a finite linear combination of the canonical monomials

    1; Z^a,W^b (a,b>=1); Y; YZ^a,YW^b (a,b>=1).       (3)

Their images are independent. Reduction never raises spatial degree, so

    dim_K V_(m,0)=1,
    dim_K V_(m,d)=4d, d>=1, m>=2.                    (4)

In particular the filtered kernel has dimension binom(d+3,3)-4d for
d>=1. Equations (1)-(3) give an exact terminating original-kernel
decision, rather than a moment-prefix criterion.

## 1. The two tangent Gaussian Stein operators

Consider the polynomial vector fields

    T_1=W partial_W-Y partial_Y,
    T_2=(1+Z)partial_Z-mW partial_W.

Both annihilate P_m, identically for every m. With the Witt covariance,
the Gaussian divergence of a vector field X is

    delta(X)=W X_Z+Z X_W+Y X_Y
                                      -partial_Z X_Z-partial_W X_W-partial_Y X_Y.

Thus delta_i A=delta(A T_i) are exactly (2). Gaussian integration by
parts, coefficient by coefficient in t, gives

    Phi_m[delta_i A]=t Phi_m[A T_i(P_m)]=0.           (5)

This proves inclusion of the right side of (1) in the original kernel
without localization or a choice of an auxiliary functional.

## 2. Degree reduction and its exact filtration

Assume m>=2. The top-degree multiplication parts of delta_2 and delta_1
are respectively -(m-1)ZW and ZW-Y^2. More precisely the operators

    L_(ZW)=-delta_2/(m-1),
    L_(Y2)=-delta_1-delta_2/(m-1)

have top multiplication terms ZW and Y^2. Their remaining terms send
any degree-e polynomial to degree at most e+1. To remove a monomial
divisible by ZW or Y^2, subtract the corresponding L applied to its
quotient. All new terms have strictly smaller total degree than that
monomial. Apply this to the largest remaining noncanonical degree.
Finite descent gives, for every degree-d marker A,

    A=N+delta_1 B_1+delta_2 B_2,

where N is a linear combination of (3), deg N<=d, and one can take
deg B_i<=d-2. The construction uses operator images at each step; it
does not multiply an already obtained relation as if the image sum
were an ideal.

At degree at most d>=1, the even-Y monomials in (3) number 2d+1 and
the odd-Y monomials number 2d-1. This proves the upper bound 4d.
The following coefficient calculation proves independence and hence
uniqueness, exactness of (1), and (4).

## 3. Exact original coefficients and rational coefficient classes

For all a,b,j,n>=0, direct Wick contraction gives

    [t^n]Phi_m[Z^a W^b Y^j]
      =(n+b)!/n! binom(mn,n+b-a) E[Y^(n+j)].          (6)

The binomial coefficient is zero outside its ordinary finite range.
This identity holds for every exponent, without a truncation bound.
For q>=0 define the nonzero base sequences

    a_q=(2q-1)!! binom(2mq,2q),
    b_q=(2q+1)!! binom(m(2q+1),2q+1),               (7)

using (-1)!!=1. Write (x)_fall_a=x(x-1)...(x-a+1) and
(x)_rise_a=x(x+1)...(x+a-1), with empty products one.
The canonical even periods have only even powers, with ratios

    [t^(2q)]Phi_m[Z^a]/a_q
       =(2q)_fall_a/(2(m-1)q+1)_rise_a,
    [t^(2q)]Phi_m[W^b]/a_q
       =(2(m-1)q)_fall_b.                           (8)

The canonical Y periods have only odd powers, with ratios

    [t^(2q+1)]Phi_m[YZ^a]/b_q
       =(2q+1)_fall_a/(2(m-1)q+m)_rise_a,
    [t^(2q+1)]Phi_m[YW^b]/b_q
       =((m-1)(2q+1))_fall_b.                       (9)

In (9), the denominator is precisely ((m-1)n+1)_rise_a at
n=2q+1; its initial term is 2(m-1)q+m. This fixes the odd indexing.

For even Z^a, a>=1, the poles are the distinct simple points

    q=-j/[2(m-1)], 1<=j<=a.

None cancels: all numerator zeros are nonnegative. Each new a introduces
a nonzero pole absent at all earlier a. The W^b ratios form a polynomial
basis with degrees b. Thus the even canonical periods are independent:
a rational function that vanishes at every q>=0 is zero; its successive
newest-pole residues kill all Z coefficients, and then polynomial
independence kills all W coefficients and the constant coefficient.

For odd YZ^a, a>=1, the simple poles are

    q=-1/2-j/[2(m-1)], 1<=j<=a.

The numerator zeros are all >=-1/2, so again none cancels. The identical
pole argument, followed by the polynomial independence of the YW^b
ratios, proves odd independence. Even and odd powers cannot cancel.
This proves the entire canonical family independent over K; the proof
is algebraic in any characteristic-zero field.

Consequently the exact degree-d literal image, m>=2,d>=1, consists of
the series

    sum_(q>=0) a_q R(q)t^(2q)
         +t sum_(q>=0) b_q S(q)t^(2q),              (10)

where

    R(q)=p(q)+sum_(j=1)^d r_j/(q+j/[2(m-1)]), deg p<=d,
    S(q)=s(q)+sum_(j=1)^(d-1) s_j/(q+1/2+j/[2(m-1)]), deg s<=d-1.

The coefficients in these rational functions are arbitrary elements
of K. To see that the pole classes are arbitrary, note that (8) and
(9) have equal numerator and denominator degrees, so their polynomial
parts are constant, and their residue matrices are triangular with
nonzero diagonal. The polynomial bases already contain constants.
This proves both directions of (10), not only inclusion. The full
literal image is the union of these explicitly given filtered spaces.

## 4. No nonzero literal period is polynomial or analytic

A nonzero rational function R(q) or S(q) has at most finitely many
zeros among nonnegative integers. Therefore (10) immediately gives

    V_m intersection K[t]={0}.                     (11)

This conclusion uses no complex analysis or field embedding.
Over C there is the stronger statement that every nonzero element of
V_m has radius of convergence zero at t=0. Indeed finite factorial
products in (7) show, for m>=2,

    lim_(q->infinity) (a_(q+1)/a_q)/q
      =lim_(q->infinity) (b_(q+1)/b_q)/q
      =2*m^(2m)/(m-1)^(2m-2)>0.                   (12)

For example, the binomial quotient in a_(q+1)/a_q is

    product_(i=1)^(2m)(2mq+i)
     /[(2q+1)(2q+2) product_(i=1)^(2m-2)(2(m-1)q+i)],

and the additional Gaussian quotient is 2q+1. The odd formula has the
same limit by shifting 2q to 2q+1. Multiplication by any nonzero rational
R(q) changes the coefficient ratio by R(q+1)/R(q)->1. Hence a nonzero
even or odd subsequence has factorial growth and radius zero as a
series in t^2. At least one such subsequence exists in every nonzero
period. Thus

    V_m intersection C{t}={0},                     (13)

where C{t} denotes convergent germs. Over other characteristic-zero
fields the coefficient statement (10) and polynomial exclusion (11)
remain the intrinsic algebraic assertions. Complex divergence can be
applied after embeddings of the finitely generated coefficient field.

Nevertheless 1 belongs to M_m. Directly from (6),

    Phi_m[1-t m(1+Z)^(m-1)Y]=1.                    (14)

For n>=1 the cancellation is the elementary identity
m binom(mn-1,n-1)=binom(mn,n). This time-dependent amplitude is allowed
in M_m. Equations (11) and (14) rigorously distinguish V_m from M_m.

## 5. The exceptional m=1 family

Equations (1),(2),(5),(6),(7) still hold when m=1, but the filtered
normal form changes. Now

    delta_2 A=WA-(1+Z)A_Z+W A_W

eliminates W with a strict degree drop. The combination
-delta_1 A+delta_2(ZA) has top multiplication term Y^2 A and all
remaining terms of degree at most deg A+1, so it eliminates Y^2.
Every marker therefore reduces without increasing degree to

    Z^a, YZ^a, a>=0.

Their coefficient ratios are respectively

    binom(2q,a), binom(2q+1,a),

which are polynomial bases of degrees a. They are independent by the
same even/odd argument. Thus (1) is again the complete kernel, while

    dim_K V_(1,d)=2d+1, d>=0.

The exact image formula (10) now uses arbitrary polynomials R,S of
degrees at most d,d-1, respectively, with no finite poles. For d=0 the
odd part is absent. Here a_q=(2q-1)!! and b_q=(2q+1)!!, whose consecutive
ratios divided by q tend to 2. The polynomial and analytic exclusions
(11),(13), and the distinction (14), hold unchanged.

## 6. Relation to the faithful rank-four baseline

For m=2, the companion [rank-four proof](GAUSSIAN_QUARTIC_BASELINE_RANK_FOUR.md)
describes M_2 and its generic
differential module, using basis (1,f,f',k). The present result gives
the complete literal original map V_2 and its filtered spatial kernel,
which are different questions. In particular the complete fixed-marker
kernel is generated as a vector-space sum by (2), while original
degree<=d images have dimension 4d even though their C(t)-rank is four.
The rank-four utility and this spatial normal form can be cross-checked
without replacing either all-order proof by finite series tests.

The companion [bounded checker](verify_gaussian_binomial_mixed_kernel.py)
verifies operator reconstruction,
coefficient ratios, degree counts, and selected direct Wick coefficients.
The unrestricted kernel, image and divergence statements above follow
from the exact descent and rational-pole arguments, not those finite
controls. This is internal proof work; no independent replay, formal
verification or external review is claimed.

The exact implementation is [gaussian_binomial_mixed_kernel.py](gaussian_binomial_mixed_kernel.py).
Its current replay is recorded in the [verification receipt](gaussian_binomial_mixed_kernel_verification.json).
