# Odd-angular D0 finite-moment obstruction

Let K be any field of characteristic zero and let
A=K[x,y,t]/(t^2+2xy-1). Define the normalized spherical functional by
Lambda(F)=CT_x I_t F(x,(1-t^2)/(2x),t), where I_t(t^j)=0 for odd j and
I_t(t^j)=1/(j+1) for even j.

Consider exactly the normalized family
\[
 f=a_1y+a_3y^3+b_1yt^2+y^3t^2+\lambda yt^4
       +xt^2(q_0+t^2).                                      \tag{1}
\]
The coefficients of y^3t^2 and xt^4 are fixed to 1. This note makes no
assertion about charts in which either of those coefficients is zero.

**Theorem 1 (Finite-moment obstruction).**

There is no choice
(a_1,a_3,b_1,q_0,lambda) in K^5 for which
\[
 \Lambda(f^2)=\Lambda(f^4)=\cdots=\Lambda(f^{14})=0.           \tag{2}
\]
Consequently no member of (1) is pure. This is a finite-prefix statement,
not a claim about eventual purity, the full normalized D0 chart, arbitrary
quintics, or the global degree-five problem.

The obstruction is given by exact rational polynomial identities and two
finite cases in lambda. It does not import an ordinary-prime theorem or
the prime-face classification.

## 1. Exact input moments and the two pivots

Write M_m=Lambda(f^m). The angular involution (x,y,t)->(-x,-y,t)
takes f to -f and preserves Lambda. Thus every odd M_m is identically
zero, at all orders. This fact is separate from the finite calculation.

The checker reconstructs every moment directly from (1), by finite
multinomial expansion and the exact monomial formula
\[
 \Lambda(x^i y^j t^k)=
 \begin{cases}
 \displaystyle\frac{2^r r!(2s)!(r+s)!}
 {s!(2r+2s+1)!},&i=j=r,\ k=2s,\\
 0,&\text{otherwise}.
 \end{cases}
\]
This follows by substituting \(xy=(1-t^2)/2\) and integrating the
resulting finite binomial expansion. No stored raw moments are needed.

Set q=q0, b=b1, and
\[
 P_1=7q+3,\qquad
 P_3=1615q^3+2261q^2+1197q+231.                              \tag{3}
\]
The exact M2 is
\[
 \frac{2}{3465}
 (231a_1q+99a_1+99bq+55b+55\lambda q+35\lambda).             \tag{4}
\]
Off P1=0 it determines a1, with denominator 33 P1. After this
substitution M4 is affine in a3, with coefficient 8 P3/4849845.
Off P1P3=0 it determines a3, whose saved reduced denominator is
1089 P1^2 P3. These two explicit rational substitutions are the
"generic_substitutions" entry of the compressed certificate.

For m=6,8,10,12,14, write the fully substituted exact moment as N_m/D_m.
The checker constructs these numerator and denominator polynomials. Each
D_m is a nonzero rational constant times powers of P1 and P3.
The respective numerator term counts are 62,140,212,366,500; their
degrees in b are 3,4,5,6,7. They are not samples or numerical fits.

The resultant Res_q(P1,P3)=2112. Thus the generic chart and the two
exceptional algebraic branches below cover every characteristic-zero
parameter point, without losing any pivot value.

## 2. The generic chart with lambda nonzero

Set h=b/lambda and z=lambda^2. If
epsilon_m is the parity of m/2, the exact numerator identities are
\[
 N_m(h\lambda,q,\lambda)
   =\lambda^{\epsilon_m} G_m(h,q,\lambda^2).                 \tag{5}
\]
The checker reconstructs the five G_m from the native moment numerators
by (5). Hence (2) implies
G6=G8=G10=G12=G14=0 on this chart.

The compressed certificate supplies the following two polynomial
identities, with every rational-function denominator already cleared:
\[
 \sum_{m=6,8,10,12,14} C_m(h,q,z)G_m(h,q,z)=D_A(q),          \tag{6}
\]
and
\[
 \sum_{m=6,8,10,14} E_m(h,q,z)G_m(h,q,z)=D_B(q).             \tag{7}
\]
Their full cleared cofactors are supplied as exact rational polynomials.
The checker verifies each complete polynomial identity exactly over QQ. In particular (6)--(7) are valid at every q; they have no
unrecorded rational-function denominator.

D_A and D_B have degrees 214 and 393. Exactly,
\[
 D_A=P_1^7P_3^3 F_A,\qquad D_B=P_1^7P_3^2 F_B,              \tag{8}
\]
where F_A,F_B are nonzero rational polynomials of degrees 198 and 380.
They are relatively prime over QQ. The package gives a small rigorous
coprimeness certificate as follows.

The checker clears rational coefficient denominators and primitive
integer contents in F_A,F_B. Both resulting integer polynomials retain
their full degrees modulo the prime p=32003. The supplied arrays u,v satisfy
u*F_A+v*F_B=1 modulo p, verified by ordinary integer convolution.
This proves QQ coprimeness: a positive-degree common QQ factor could be
chosen primitive over Z by Gauss's lemma; the nonvanishing leading
coefficients of both inputs modulo p force that factor to retain
positive degree modulo p, contradicting the displayed modular Bezout
identity. This is not an inference from a modular multivariate unit
ideal to a rational unit ideal.

If all G_m vanished, equations (6)--(8), with P1P3 nonzero, would imply
F_A(q)=F_B(q)=0. QQ coprimeness rules this out in every characteristic-zero
extension. No irreducibility claim or large q-resultant is needed.

## 3. The generic chart with lambda zero

This case is proved separately, without the change b/lambda and without
an all-order prime-face input.

First suppose b is nonzero and put z=b^2. The exact polynomials obtained
from N6/b, N8, and N10/b are respectively
\[
 L=A(q)z+B(q),\quad
 Q=C(q)z^2+D(q)z+E(q),\quad
 Q'=C'(q)z^2+D'(q)z+E'(q).                                 \tag{9}
\]
For either quadratic the saved exact identity is
\[
 -(ACz+AD-BC)L+A^2Q=A^2E-ABD+B^2C.                         \tag{10}
\]
The two right sides are the corresponding z-resultants. Their QQ
univariate extended-gcd certificate has gcd equal to a nonzero rational
multiple of P1^4 P3^2. The complete rational coefficient arrays and
ordinary Fraction-convolution verification are saved in
the certificate entry "lambda0_univariate_Bezout". They cannot both
vanish on the generic chart. Formula (10) requires no division by A,
so a possible vanishing leading coefficient introduces no missing case.

If b=0 as well, the scalar q-polynomials from N8 and N12 have gcd a
nonzero rational multiple of P1^4. A second explicit QQ extended-gcd
identity in the same certificate rules out their simultaneous
vanishing whenever P1 is nonzero.

Thus the entire lambda=0 generic chart is excluded.

## 4. All exceptional pivot values

If P1=0, then q=-3/7, and (4) is exactly
16(11b+10lambda)/24255. Thus b=-10lambda/11.
M4 then solves a3 with a nonzero constant denominator, leaving a1
and lambda free. Substituting these identities into M6,...,M14 gives
five polynomials I1,...,I5 in QQ[a1,lambda]. The five literal saved
cofactors satisfy
\[
 \sum_{i=1}^5 U_i I_i=1.                                   \tag{11}
\]
The checker verifies the complete identity by exact polynomial
multiplication.

If P3=0, P1 is invertible modulo P3, so M2 determines a1 by its saved
polynomial representative in QQ[q,b,lambda]/(P3). Keep a3 free.
The six reduced polynomials from M4,...,M14, followed by P3 itself,
give seven polynomials J1,...,J7 in QQ[a3,b,lambda,q], with J7=P3.
Six supplied cofactors are checked modulo P3; exact division of
\(1-\sum_{i=1}^6 V_iJ_i\) by P3 recovers the seventh. Thus
\[
 \sum_{i=1}^7 V_i J_i=1.                                   \tag{12}
\]
This single polynomial identity covers all three geometric roots of P3;
it does not choose one embedding or divide by an expression that might
vanish there. It too was independently multiplied exactly.

Equations (11)--(12) exclude the exceptional branches, including
lambda=0, and complete the finite-moment obstruction.

## 5. Verification and limits

The accompanying checker
[verify_odd_angular_D0_finite_moment_certificate.py](certificates/verify_odd_angular_D0_finite_moment_certificate.py)
takes the compressed
[exact certificate input](certificates/odd-angular-d0-finite-moment-certificate.json.gz). It uses only
Python standard-library exact arithmetic and python-flint. It requires
no Singular runtime, original experiment directory, stored raw moments,
network access, or ordinary-prime theorem.

From the repository root, run:

```sh
python research/foundations/extended-geometry/certificates/verify_odd_angular_D0_finite_moment_certificate.py
```

The retained [verification record](certificates/odd-angular-d0-finite-moment-verification.json)
records a completed exact replay and its checker/input digests.

The checker reconstructs native M1 through M14, verifies the two pivot
identities and their complete exceptional substitutions, constructs all
generic reduced numerators, and checks every cleared polynomial
identity. It verifies the prime, primitive contents, full leading
degrees, and modular Bezout convolution in the QQ-coprimeness argument.
It also checks both rational lambda=0 extended-gcd identities and both
exceptional constant-one identities, including recovery of the seventh
P3 cofactor. Deliberate coefficient, cofactor, native-moment, and modular
certificate corruptions must be rejected.

This is a finite chart certificate, rather than one expanded global
five-variable Bezout expression. Every chart covers all its
characteristic-zero points. The modular multivariate diagnostic and
any incomplete Groebner computation are not proof inputs.

The conclusion concerns exactly (1) and moments through14. It does not
establish eventual-purity exclusion, a result for the full D0 chart, or
a global degree-five classification.
