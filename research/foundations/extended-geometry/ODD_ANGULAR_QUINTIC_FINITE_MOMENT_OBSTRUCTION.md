# Sharp finite-moment obstruction for odd-angular quintics

Let K be any field of characteristic zero, and put
\[
 A_K=K[x,y,t]/(t^2+2xy-1).
\]
The normalized spherical functional used here is
\[
 \Lambda(F)=\operatorname{CT}_x I_t
 F\left(x,\frac{1-t^2}{2x},t\right),\qquad
 I_t(t^{2j})=\frac1{2j+1},\quad I_t(t^{2j+1})=0.
\]
Consider the normalized eight-monomial family
\[
 f=a_1y+a_3y^3+a_5y^5+b_1yt^2+y^3t^2+\lambda yt^4
       +xt^2(q+t^2).                                      \tag{1}
\]
All six displayed parameters are arbitrary elements of K. The coefficients
of y^3t^2 and xt^4 are fixed to1.

**Theorem 1 (Finite-moment obstruction).**

For every parameter point of (1), at least one of
\[
 \Lambda(f^2),\Lambda(f^4),\ldots,\Lambda(f^{14})             \tag{2}
\]
is nonzero. In particular no member of (1) is pure.

Here pure means that every positive moment vanishes. The theorem is a
finite necessary-moment obstruction, with a complete rational certificate.
It does not concern eventual purity. The uniform cutoff14 is sharp as
stated in Theorem3.

**Corollary 2 (Nonzero coefficient normalization).**

The same finite14 obstruction holds for
\[
 a_1y+a_3y^3+a_5y^5+b_1yt^2+C y^3t^2+Ayt^4
       +xt^2(q+B t^2),\qquad BC\ne0,                       \tag{3}
\]
over every characteristic-zero field.

Choose a fourth root rho of C/B in an algebraic extension of K, and set
mu=1/(rho B). The actual quadric automorphism
x->rho*x, y->rho^(-1)*y, t->t preserves Lambda. Multiplication by mu
therefore transforms every mth moment by the nonzero factor mu^m.
The coefficients of xt^4 and y^3t^2 both become1, and every other term
retains its form in (1). Apply Theorem1 over that extension and use the
injectivity of the original field map. No ordering, complex embedding,
or choice of a fourth root in K itself is needed.

There is also a necessary condition in the full even-t stripped space.
Write
\[
 f=P(y,t)+xt^2(q+Bt^2),\qquad \deg P\leq5,
\]
where \(P\) is even in \(t\), and let \(C=[y^3t^2]P\).
If \(BC\ne0\) and the first14 moments vanish, then the angular-even part
\[
 \frac{P(y,t)+P(-y,t)}2
   =p_{00}+p_{02}t^2+p_{04}t^4+p_{20}y^2
       +p_{22}y^2t^2+p_{40}y^4
\]
is nonzero. Indeed, if these six coefficients were all zero, the phase
would be exactly (3), contradicting Corollary2. This is a finite-prefix
necessary condition; it does not assert that a nonzero angular-even
part is sufficient for vanishing.

Neither statement includes B=0 or C=0. The family omits the angular-even
terms of the full even-t stripped space. Thus these results do not
classify that larger space or resolve the global degree-five problem.

## 1. Native moment polynomials and the three pivots

For nonnegative integers i,j,k the exact native monomial formula is zero
unless i=j=r and k=2s; in that case
\[
 \Lambda(x^r y^r t^{2s})
   =\frac{2^r r!(2s)!(r+s)!}{s!(2r+2s+1)!}.                \tag{4}
\]
Substitute xy=(1-t^2)/2 and integrate the finite binomial expansion to
obtain (4). The checker reconstructs the complete moment polynomials
M_m=Lambda(f^m), for m1 through14, by literal multinomial expansion of
the eight terms of (1). No stored raw moments are used.

The angular antipode (x,y,t)->(-x,-y,t) preserves Lambda and takes f to-f.
Consequently every odd moment vanishes at all orders. The even native
moment term counts through14 are6,26,76,189,404,799,1452.

Write
\[
\begin{split}
 P_1&=7q+3,\\
 P_3&=1615q^3+2261q^2+1197q+231,\\
 P_5&=930465q^5+2225025q^4+2314026q^3\\
    &\hspace{13mm}+1285570q^2+376805q+46189.
\end{split}                                               \tag{5}
\]
The native coefficients of a1 in M2, a3 in M4, and a5 in M6 are
\[
 \frac{2P_1}{105},\qquad
 \frac{8P_3}{4849845},\qquad
 \frac{16P_5}{300840735195}.                               \tag{6}
\]
These pivot coefficients do not involve the earlier unknown parameters,
by angular balance. In particular the last is also
\[
 \frac6{32} I_t[(1-t^2)^5t^{10}(q+t^2)^5].
\]
The full M2 is
\[
 \frac2{3465}
 (231a_1q+99a_1+99b_1q+55b_1+55\lambda q+35\lambda).        \tag{7}
\]

The three polynomials in (5) are squarefree and pairwise coprime.
Their ordered resultants are2112,4756480, and
971106570043517241917440 for the pairs(1,3),(1,5),(3,5).
For example Res(P1,P3)=7^3 P3(-3/7)=2112. The checker verifies the
separation by exact polynomial inverse identities as well.

On P1P3P5 nonzero, M2=M4=M6=0 uniquely solve a1,a3,a5. The supplied
rational solutions have monic denominators proportional to
\[
 P_1,\qquad P_1^2P_3,\qquad P_1^3P_3P_5,                  \tag{8}
\]
of degrees1,5,11. Their numerators have4,21,62 terms.
After substitution, let N_m/D_m be the exact fraction for m8,10,12,14,
with N_m a primitive integer polynomial, normalized to positive leading
coefficient. Every D_m has only the pivot factors. The respective
numerator term counts are185,272,526,700.

The checker verifies every rational identity by exact polynomial
cross multiplication. All coefficients are universal rational numbers,
so these identities apply to arbitrary characteristic-zero fields,
including transcendental coefficient fields.

## 2. Generic q and lambda nonzero: iterated resultants

Assume P1P3P5 is nonzero and lambda is nonzero. Put
\[
 h=b_1/\lambda,\qquad z=\lambda^2.
\]
If epsilon_m is the parity of m/2, the exact numerator identities are
\[
 N_m(h\lambda,q,\lambda)
       =\lambda^{\epsilon_m}G_m(h,q,\lambda^2).             \tag{9}
\]
Thus vanishing of the native prefix implies G8=G10=G12=G14=0.
No condition z+2 nonzero is imposed.

The small exact z-resultants have the following contents:
\[
 \operatorname{Res}_z(G_8,G_m)
      =c_m P_1^{u_m}P_3^{v_m}P_5^{w_m}R_m(h,q),           \tag{10}
\]
where R_m is primitive over the integers.

| m | (u_m,v_m,w_m) | c_m | degrees(h,q) of R_m |
|---|---|---|---|
|10|(8,4,1)|34726138722|(10,59)|
|12|(12,6,3)|73213257172147488|(12,80)|
|14|(12,6,3)|73213257172147488|(14,84)|

All constants and removed pivot factors are units on the present
characteristic-zero chart. Hence an actual common zero of the G_m
is a common zero of the three R_m. This implication also holds when a
leading coefficient in z drops degree: the full Sylvester determinant
still vanishes at a common finite root. No leading coefficient is
silently divided out.

Next form Res_h(R10,R12) and Res_h(R10,R14). Their primitive q factors,
after removing respectively P1^120 and P1^140, have degrees1388 and1526.
They are coprime over QQ. The retained coefficient arrays and the exact
CRT verification below certify the complete two resultant identities;
the comparison is not an evaluation at a few sample q values.

To certify coprimeness, take the primitive integer forms of the two
residual factors. Both full leading degrees survive modulo the prime
32003, and the supplied polynomial arrays u,v satisfy
\[
 uF_{1388}+vF_{1526}=1\pmod {32003}.                        \tag{11}
\]
The checker verifies primality and the complete coefficient convolution.
If a positive-degree common rational factor existed, Gauss's lemma
would give a primitive integer common factor. Since neither input loses
its leading degree, that factor would retain positive degree modulo
the prime, contradicting (11). Thus there is no common root over any
characteristic-zero extension.

An actual common zero would force both residual q factors to vanish,
which is impossible. This excludes the generic q, lambda-nonzero chart.

### Exact verification of the large resultants

For integer polynomials A(h,q),B(h,q), put
d_a=deg_h A, d_b=deg_h B, and let ||A||_1 denote the sum of the absolute
values of ALL its integer coefficients. The Sylvester determinant has
q degree at most
\[
 D=d_b\deg_q A+d_a\deg_q B.
\]
Every coefficient of that determinant has absolute value at most
\[
 H=\|A\|_1^{d_b}\|B\|_1^{d_a}.                             \tag{12}
\]
Indeed the polynomial coefficient l1 norm is submultiplicative.
The determinant's permutation products are a subset of all row-wise
entry products; the product of row l1 norms therefore bounds its l1
norm, and consequently each coefficient. This proves (12) without a
factorial or an analytic approximation.

The candidate integer polynomial is its supplied primitive polynomial
times an integer scalar. That scalar is determined, and checked, by one
exact integer h-resultant at a q value preserving both h degrees and
having nonzero candidate value. This scalar evaluation alone is not
the certificate.

For each of many distinct primes, primality is verified by complete
trial division. The checker computes the full modular h-resultant,
requires that both input h degrees survive, and compares EVERY
coefficient through degree D with the candidate. If C is the largest
absolute claimed coefficient and the prime product is greater than
H+C, every coefficient of the integer difference is divisible by that
product and has smaller absolute value; hence it is zero.

The two degree bounds are1508 and1666. The true-coefficient bounds have
7601 and8459 bits;246 and273 verified31-bit primes yield products with
7626 and8463 bits, respectively, exceeding the required exact bounds.
The final P1 divisions and the primitive contents are checked
literally. This is an exact Chinese-remainder argument, not a
probabilistic polynomial identity test.

## 3. Generic q and lambda zero

Assume P1P3P5 is nonzero and lambda=0.

If b1 is nonzero, set z=b1^2 and divide exactly the odd-b1 numerator
parity factor at order10. The resulting F8,F10,F12 have z degrees2,2,3.
Literal Sylvester-cofactor identities give
\[
 U_1F_8+V_1F_{10}=S_{10}(q),\qquad
 U_2F_8+V_2F_{12}=S_{12}(q).                              \tag{13}
\]
No leading coefficient is inverted. The scalar degrees are84 and125.
After removal of their common pivot factor
P1^8 P3^4 P5, the primitive integer factors, of degrees59 and100,
have a complete modular Bezout certificate at32003 with both leading
degrees preserved. The same Gauss argument proves QQ coprimeness.

If b1=0 as well, the reduced scalar M8 and M12 numerators have pivot
factor P1^4P5. Their remaining primitive factors have the same type of
certified QQ coprimeness. Both cases are therefore impossible on the
generic q chart. Only the finite prefix through12 is used here.

These are finite necessary-moment arguments. No all-order prime-face
exclusion or ordinary-reduction assertion is used at lambda=0.

## 4. All exceptional q fibers

The three pivot polynomials are separated, so their exceptional
families consist of1,3,5 distinct geometric q values. The quotient
identities below are actual QQ polynomial congruences and therefore
hold at every root, over any characteristic-zero field.

At P1=0, q=-3/7. Equation(7) gives b1=-10lambda/11.
M4 uniquely solves a3 and M6 uniquely solves a5, retaining a1,lambda.
Four supplied rational-polynomial cofactors of the exact native residual
M8,M10,M12,M14 have sum a nonzero rational constant.

At P3=0, solve a1 modulo P3 and keep a3 and the M4 compatibility.
The P5 coefficient is a unit modulo P3, so M6 solves a5.
Five supplied cofactors of the native residual M4,M8,M10,M12,M14
have sum a nonzero rational constant modulo P3.
This excludes ALL three roots for arbitrary lambda.

At P5=0, P1 and P3 are units. Solve a1,a3 modulo P5, keep a5 and the
M6 compatibility. For lambda=0, five supplied cofactors of the native
residual M6,M8,M10,M12,M14 have sum a nonzero rational constant modulo P5.

In the compressed input, cofactor vectors have been cleared of their
common rational contents, and the corresponding nonzero rational
normalizing constants are retained explicitly. The checker verifies the
complete native products. No old certificate for a5=0 is imported.

It remains to treat P5=0 with lambda nonzero.

## 5. P5 and lambda nonzero: every h-pivot is retained

Keep the native P5 residuals from Section4 and set
\[
 h=b_1/\lambda,\qquad z=\lambda^2,\qquad w=a_5/\lambda.
\]
An actual parameter point has z nonzero. Dividing precisely the native
lambda parity factors gives G6,G8,G10,G12,G14 over QQ[q]/(P5), with
\[
 G_6=A(h)z+B(h),\qquad
 G_8=zC(h)w+E_2(h)z^2+E_1(h)z+E_0(h).                    \tag{14}
\]
The degrees of A,B,C in h are3,1,1. Put
\[
 R_8=E_2B^2-E_1AB+E_0A^2,\qquad P=ABC.
\]
Two explicit pairs of polynomial cofactors satisfy
\[
 U_AA+V_AB=c_A,\qquad U_CC+V_CR_8=c_C\pmod {P_5},          \tag{15}
\]
where c_A,c_C are nonzero rational constants.

If A=0, G6=0 forces B=0, contradicting the first identity.
Thus A is nonzero; since z is nonzero, G6=0 also makes B nonzero.
If C=0, G8=0 and z=-B/A force R8=0, contradicting the second identity.
Consequently ALL possible vanishing h-pivots have been excluded, and
\[
 z=-B/A,\qquad w=R_8/(ABC)                                \tag{16}
\]
are valid at every remaining actual point.

Substitute (16) into G10,G12,G14 and multiply by exactly the known
denominator powers. The resulting h-polynomials H10,H12,H14 have
degrees12,19,20. Their complete identities are checked by reconstruction
from the native moments, not by trusting substituted sample values.

Three supplied cofactors and a polynomial J of h degree3 satisfy
\[
 U_{10}H_{10}+U_{12}H_{12}+U_{14}H_{14}=c_HJ,\qquad
 QJ=(ABC)^3\pmod {P_5},                                   \tag{17}
\]
where c_H is a nonzero rational constant.
If all three H_m vanished, J would vanish, and (17) would contradict
ABC nonzero. This excludes the entire P5, lambda-nonzero chart.

The argument includes h=0 and a5=0, and imposes no z+2 localization.
It does not presume a generic endomorphism ring or select only one
embedding of QQ[q]/(P5). Combining Sections2--5 proves Theorem1.

## 6. Sharpness of the order14 cutoff

**Theorem 3 (Sharpness over characteristic-zero fields).**

There is a member of (1) over Q_101 with a5 nonzero whose first13 positive
moments vanish and whose14th moment is nonzero. There is also such a
member over some number field. Thus14 is the least uniform
nonvanishing cutoff over characteristic-zero coefficient fields,
already for genuine degree-five members of the normalized family.
No rational or real coefficient realization is asserted.

All native coefficients of M2,M4,...,M14 are integral at101. At the
point, in the exact order(a1,a3,a5,b1,q,lambda),
\[
 (23,44,98,12,89,94)\pmod {101},
\]
the complete native computation gives
\[
 M_2=M_4=\cdots=M_{12}=0,\qquad M_{14}=56\pmod {101},        \tag{18}
\]
and
\[
 \det\frac{\partial(M_2,M_4,M_6,M_8,M_{10},M_{12})}
              {\partial(a_1,a_3,a_5,b_1,q,\lambda)}
       =73\pmod {101}.                                   \tag{19}
\]
The checker reconstructs these polynomials and their full6-by6
Jacobian directly. In particular (19) is not a determinant of only
generic triangular equations.

Here is the lifting argument. Let F be the six-polynomial vector in
(19), viewed over Z_101, and start at the displayed integer point v0.
Its Jacobian is invertible over Z_101. Inductively use
\[
 v_{n+1}=v_n-J_F(v_n)^{-1}F(v_n).
\]
If F(v_n) is in101^(2^n) Z_101^6, the correction has the same divisibility.
The exact polynomial Taylor identity has remainder quadratic in the
correction, so F(v_(n+1)) is in101^(2^(n+1)) Z_101^6.
The Jacobian remains a unit modulo101. Thus the sequence converges to
a Z_101 point v with F(v)=0 and the same reduction as v0.
Polynomial Taylor expansion here uses integral binomial coefficients;
no division by factorials is required.

The odd moments vanish identically by angular parity. Equations(18)
and the unchanged residue give M14(v) a101-adic unit; a5(v) is a unit
as well. This proves the Q_101 assertion. It is a smooth-point lifting
argument, not an inference from a bare proper modular ideal.

For the number-field assertion, consider the finite-type Q-algebra
\[
 \mathcal B=
 \mathbb Q[a_1,a_3,a_5,b_1,q,\lambda,u]/
 (M_2,M_4,M_6,M_8,M_{10},M_{12},\,u a_5M_{14}-1).
\]
The lifted point defines a map to Q_101, so this algebra is nonzero.
It has a maximal ideal, and its residue field is finite over Q.
For completeness, a field finitely generated as a Q-algebra is algebraic:
otherwise choose a nonempty transcendence basis, clear the finitely many
denominators in the monic equations for the remaining generators, and
make the field integral over Q[t1,...,tr,1/d]. A subring over which a
field is integral is itself a field, whereas that localization is not
a field: some t1-c, c in Q, does not divide d and remains a nonunit.
This contradiction proves algebraicity, and finite generation gives
finiteness. The resulting number-field point retains a5*M14 nonzero
and gives the required moment behavior.

The witness is not pure, since its14th moment is nonzero. Sharpness
therefore does not supply a degree-five counterexample to purity.

## 7. Portable verification and scope

The accompanying
[checker](certificates/verify_odd_angular_quintic_finite14_certificate.py)
uses the compressed
[exact coefficient input](certificates/odd-angular-quintic-finite14-certificate.json.gz).
From the repository root, run

    python research/foundations/extended-geometry/certificates/verify_odd_angular_quintic_finite14_certificate.py

It needs only Python standard-library exact arithmetic and python-flint.
There is no Singular or SymPy runtime, network access, original
experiment directory, or stored raw moment input.

The checker reconstructs native moments and pivots; verifies every
rational substitution, q-content division, complete CRT comparison and
Gauss coefficient identity; checks all exceptional native residual
cofactors and every P5 h-pivot; and verifies the sharpness point and
full native Jacobian. Deliberate alterations of native data, cofactors,
modular identities, CRT candidates, insufficient modulus products and
sharpness data must be rejected. The retained
[verification record](certificates/odd-angular-quintic-finite14-verification.json)
contains the completed replay and exact code/input digests.

All assurance flags remain separate from the mathematical theorem.
The proof is an exact finite chart certificate, not a formalization.
It excludes the displayed odd-angular family with BC nonzero and
establishes its sharp finite14 cutoff. The full even-t space, other
stripped families, eventual purity, and the global degree-five problem
are outside its conclusion.
