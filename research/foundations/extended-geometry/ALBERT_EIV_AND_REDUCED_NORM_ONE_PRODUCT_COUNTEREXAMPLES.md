# Compact EIV and reduced Albert norm-one products

The same weighted Peirce calculation gives an explicit beta functional on two different Albert homogeneous spaces. In the compact unitary model it is the squared modulus of a diagonal coordinate. On a norm-one hypersurface it is the product of a diagonal coordinate and the corresponding adjugate coordinate. A polynomial cumulative distribution then pulls back the interval--Laurent construction to an explicit phase with one fixed unit marker at every order.

## 1. Two distinct models and their counterexamples

Let J_Q=Herm_3(O_Q) be the positive rational Albert algebra. Write e for its unit, N for its cubic norm, X# for its adjugate, tr for its Jordan trace, and B(X,Y)=tr(X o Y). Let E=Q(i). The involution tau on J_E conjugates ONLY this scalar i; octonion conjugation, used in defining Hermitian matrices, is a different E-linear operation.

Let G_c be the E-linear cubic-norm and Hermitian-form preserving group on J_E, where h(X,Y)=B(X,tau Y), and let H_c be its FULL stabilizer of e. Define the actual rational quotient

X_c=G_c/H_c.

It is the compact rational E6/F4 model, a smooth geometrically integral affine variety of dimension26. For C=g e, set w=C_11 and

D_c=w tau(w), v_c=1-(1-D_c)^8.                         (1.1)

Theorem 1. For every characteristic-zero field K, on the scalar extension X_c,K times the SPLIT multiplicative group, put

F_c=(1-t^-1)(1-v_c+v_c t), H=t^-1.

For the normalized invariant trace Lambda=lambda_X tensor CT_t and every m>=1,

Lambda(F_c^m)=0,
Lambda(H F_c^m)=(-1)^(m-1)/(m+1) !=0.                  (1.2)

The phase is not a unit times any proper power after algebraic closure. Its nonzero Laurent support is not contained in a translate of dZ for any integer d>=2. If K is formally real, each individual factor X_c,K and split Gm is Mathieu, while the product is not.

The second theorem has a different field and model scope. Let K be ANY characteristic-zero field and let J=Herm_3(O,gamma) be ANY reduced Albert algebra overK, with nonzero diagonal gamma parameters and its given K-defined frame e1,e2,e3. Division Albert algebras without such a frame are excluded. Let

Y_J={C in J:N(C)=1}.

Theorem 2. The actual norm-one hypersurface Y_J is smooth and geometrically integral and is homogeneous under the reduced structure group Str_0(J), of type E6, with full F4 stabilizer Aut(J). Put

D_J=C_11(C#)_11, v_J=1-(1-D_J)^8,
F_J=(1-t^-1)(1-v_J+v_J t), H=t^-1.                    (1.3)

Its normalized invariant trace satisfies (1.2), with F_J in place of F_c, at EVERY positive order. The phase is geometrically nonpower even up to a unit and has the same radial support property. Thus the product trace kernel is not Mathieu. These identities persist under every further characteristic-zero extension ofK.

No factor-Mathieu assertion is made for the generally noncompact norm-one forms. For example, the positive real Albert norm-one hypersurface is the noncompact E6(-26)/F4 model, whereas X_c is compact E6/F4. Agreement after complexification does not identify those real or rational varieties. Nor is a phase chosen only after splitting silently descended to a division Albert algebra.

In both theorems the diagonal functional is

lambda(Q(D))=integral_0^1 8(1-u)^7 Q(u)du             (1.4)

for every polynomial Q, interpreted as rational polynomial integration. For the compact real model this is the actual Beta(1,8) probability pushforward. Its relation to classical Shilov--Hardy moments is explained in Section9; no novelty is asserted for that moment law or for the interval--Laurent seed.

## 2. Albert identities and geometric splitting

For the cubic Jordan algebra write L_A(C)=A o C. We use the standard identities

dN_C(Y)=B(C#,Y), C# o C=N(C)e,
(C#)#=N(C)C, B(A o C,Y)=B(C,A o Y).                   (2.1)

They are polynomial identities of the Albert algebra. Its trace form is nondegenerate and associative. The Cayley--Dickson construction supplies the positive rational model; over a general reduced algebra the same identities include the gamma-scaled Peirce forms.

The geometric splitting used below can be seen explicitly. In characteristic different from2, write an octonion composition algebra in Cayley--Dickson form with nonzero doubling parameters alpha,beta,delta. This presentation is obtained by successively choosing nonisotropic orthogonal generators for a quadratic subalgebra, a quaternion subalgebra, and its orthogonal complement. Over an algebraic closure, rescaling the three doubling generators by square roots changes each nonzero parameter to the fixed standard parameter. Thus the octonion algebra becomes the standard split algebra. If Gamma=diag(gamma_i), choose a scalar diagonal S with S^2=Gamma. The map C->S C S^-1 takes Gamma-Hermitian matrices to ordinary Hermitian matrices and preserves the Jordan product. It also preserves the normalized cubic norm. Hence every reduced Albert algebra used here becomes the same standard Albert algebra after algebraic closure. No classification of division Albert algebras is needed.

A norm preserver a satisfies, by the chain rule and nondegeneracy of B,

(a C)#=a^-T C#,                                      (2.2)

where T means transpose relative to B. Norm preservers are stable under inverse transpose: by (2.1)-(2.2),

N(a^-T C#)=N((aC)#)=N(aC)^2=N(C)^2=N(C#).

The map C->C# is dominant. Indeed its differential at e is Y->tr(Y)e-Y, with eigenvalues2 on the scalar line and-1 on the trace-zero subspace. The last polynomial identity therefore implies N(a^-T Y)=N(Y) for all Y.

The full stabilizer of e in the cubic-norm-preserver group is the Jordan automorphism group. Polarize N at e: its linear coefficient gives tr(Y), and its quadratic coefficient gives ((tr Y)^2-B(Y,Y))/2. Polarizing the cubic coefficient then recovers tr(X o(Y o Z)), and nondegeneracy of B recovers the product. Thus a norm preserver fixing e preserves the Jordan product, and conversely a Jordan automorphism preserves N and fixes e. This is the algebraic equivalence of Yokota Lemma2.2.4.

## 3. Full group geometry and the actual varieties

The required exceptional-group identification is the explicit Albert model of Yokota [Y]. The full complex cubic-norm-preserver group is connected, simply connected, and of type E6; its full stabilizer of e is F4. Theorem3.14.1 gives the full group statement, and Section3.7 identifies the stabilizer, rather than only its Lie algebra. All coefficient groups here are algebraic groups in characteristic zero, so these identifications apply to their smooth group schemes and descend to the stated forms.

The action is transitive on the ENTIRE complex hypersurface N=1. Here is the extra step connecting the primary compact diagonalization to that algebraic statement. Yokota Proposition3.8.2 transforms any complex Albert C to a diagonal form by a compact E6 element, with the last two diagonal entries nonnegative real. If N(C)=1, they are nonzero and the product of all three entries is1. Hence the first entry is positive real too. Write these entries lambda_i and choose a_i=lambda_i^-1/2>0, with product_i a_i=1. The linear map multiplying diagonal coordinate i by a_i^2 and octonion entry ij by a_i a_j preserves N: every diagonal-product, diagonal-times-octonion-norm, and triple-octonion term is multiplied by (a1 a2 a3)^2. It maps the diagonal C to e.

Thus the full connected E6 group acts transitively on N=1. Transitivity and the homogeneous-space isomorphism for this rational model can be tested after the faithfully flat extension Qbar toC. By Section2 every reducedJ has the same geometric model. Its full stabilizer is Aut(J), so the orbit map identifies Str_0(J)/Aut(J) with Y_J overK by descent. Euler's identity dN_C(C)=3 on N=1 proves smoothness directly. The orbit of the connected group is irreducible; thus Y_J is geometrically integral and has dimension26. It is affine as the displayed hypersurface.

For the compact rational group, after extending to Qbar the algebra E splits into two copies. An E-linear transformation becomes (a,b), and the Hermitian equation says b=a^-T. The norm condition for b follows from the one for a by Section2. Consequently G_c,Qbar is exactly the same full E6 norm-preserver group. Fixing e gives a Jordan automorphism a, for which a^-T=a; the pair is (a,a). Descent identifies H_c with Aut(J_Q), the full positive rational F4 form. The quotient X_c is smooth and geometrically integral, and its geometric base change is the affine norm-one homogeneous space already proved above. Affineness descends under field extension, so X_c itself is affine. This also explains why no extra Matsushima input is needed here.

Its actual orbit coordinate C=g e descends to X_c. A closed-embedding assertion about this map is unnecessary. In particular the real set {N(C)=1,h(C,C)=3} is not being used as the Q coordinate scheme: those equations alone have the wrong dimension after complexification. Yokota Theorem3.8.3 identifies that REAL set with the compact orbit, which is consistent with the actual quotient used here.

Both models now satisfy the hypotheses of the same-field normalized Reynolds construction [M]. For their connected reductive acting groups, invariants in the coordinate ring are constants. The unique normalized projection to constants is Lie-invariant and commutes with every field extension; tensoring with Laurent constant term gives the product trace. These are precisely the generic trace assertions used below. In the compact model Haar integration on the dense compact real orbit is that trace, by invariance and uniqueness.

## 4. Weighted Jordan Lie coordinates

For a Jordan algebra every commutator [L_A,L_H] is a derivation. The algebraic proof is included. Full polarization of the Jordan identity gives

((a o b)o d)o c+((b o c)o d)o a+((c o a)o d)o b
 =(a o b)o(d o c)+(b o c)o(d o a)+(c o a)o(d o b).

Subtract its instance with triple(b,c,z), middle a from the instance with triple(a,c,z), middle b, or take the negative. The paired product terms cancel by commutativity; the remaining equality is

[L_a,L_b](c o z)=([L_a,L_b]c)o z+c o([L_a,L_b]z).

Thus these are genuine Lie Aut(J) operators. They annihilate e and preserve the trace and norm. Associativity of B makes every L_A selfadjoint, hence makes its commutator skew-adjoint. There is also the direct formula

dN_C(L_A C)=B(C#,A o C)=B(C# o C,A)=tr(A)N(C).         (4.1)

In particular a trace-zero L_A is a Lie operator of the norm-preserver group. In the compact realification, iL_A is skew-Hermitian and also norm-preserving. This agrees with the full compact Lie decomposition der(J) direct-sum iL(J_0) proved in Yokota Theorem3.2.4.

Fix a frame e1,e2,e3 and diagonalize the nondegenerate scalar forms on the two Peirce spaces J_12,J_13, each of dimension8. Choose K-bases H_ia, i=2,3, a=0,...,7, such that

H_ia o H_ib=0 for a!=b,
H_ia^2=d_ia(e1+ei), d_ia in K^times.                  (4.2)

Diagonalization of a nondegenerate symmetric form in characteristic zero requires no extension ofK. The d_ia may be nonsquares or indefinite. For the positive rational compact model use the standard orthonormal octonion units, so d_ia=1.

Expand a general C with diagonal coordinates w_i and Peirce coordinates z_ia. Put w=w_1, A=e1-ei, H=H_ia and Delta_ia=2[L_A,L_H]. A two-coordinate computation gives

L_H w=d_ia z_ia, L_H z_ia=(w+w_i)/2,
Delta_ia w=2d_ia z_ia, Delta_ia z_ia=w_i-w.           (4.3)

For example L_A has zero coefficient on its own J_1i component, while L_H takes that component to d_ia times the e1-coordinate. Other Peirce components contribute neither of the two specified coordinates. No associative identity for octonion matrices is assumed.

If a second coordinate C^flat transforms contragrediently under the trace form, a Lie matrix V acts on it as -V^T C^flat. Since Delta^T=-Delta and L_H^T=L_H, the operator acting as Delta+2L_H on C acts as Delta-2L_H on C^flat. Therefore, writing w^flat,z_ia^flat for its coordinates,

E_ia w=4d_ia z_ia, E_ia w^flat=0,
E_ia z_ia^flat=-2w^flat.                              (4.4)

These weighted coefficients are important: the d_ia occur both here and in the e1-coordinate of the norm identity, so no inverse Gram term can be discarded.

## 5. Instantiating the dual coordinate in each model

On Y_J the coordinate C^flat=C# transforms contragrediently by (2.2), and

C o C^flat=e.

Here E_ia=Delta_ia+2L_Hia is an actual Lie operator overK, with (4.4). No scalar i is introduced.

On the compact model unitarity gives g^-T=tau g tau. Adjugate covariance and e#=e give

C#=(g e)#=g^-T e=tau(g e)=tau C.

Thus C^flat=tau C transforms in the required dual representation, and C o C^flat=e. The real compact Lie operator B_ia=2iL_Hia acts by2iL_H on C and by-2iL_H on tau C. After extending the trace toE, the operator

E_ia=Delta_ia-i B_ia

acts by Delta+2L_H on C and Delta-2L_H on tau C. Consequently it has exactly (4.4), with d_ia=1. In particular E(tau w)=0 and E(tau z_ia)=-2tau w. It is essential here that tau conjugates scalar i; octonion conjugation is a separate operation.

In either model the e1-coordinate of C o C^flat=e is

D+sum_{i,a} d_ia z_ia z_ia^flat=1,
D=w w^flat.                                         (5.1)

In the compact model this is an ordinary sum of complex squared moduli, and0<=D<=1 on its real compact orbit. In the norm-one model it is a polynomial identity overK, without a positivity assertion. All the operators used above induce vector fields on the actual homogeneous coordinate algebra, so its invariant trace annihilates their derivatives.

## 6. One all-order recurrence for both models

Let mu_j=lambda(D^j), with mu_0=1. For each j>=1 and each of the16 Peirce basis elements apply Lie invariance to the regular polynomial

w^j (w^flat)^(j-1) z_ia^flat.

There are no negative coordinate powers, including at j=1. By (4.4),

0=4j d_ia lambda(D^(j-1)z_ia z_ia^flat)-2mu_j.

Sum these16 identities and use (5.1). This gives

0=4j(mu_(j-1)-mu_j)-32mu_j,
(j+8)mu_j=jmu_(j-1).                                 (6.1)

The possibly indefinite Peirce factors d_ia already occur in both equations and cancel in the sum exactly. Since characteristic is zero, induction gives, for EVERY j>=0,

mu_j=j!/(9)_j=8! j!/(j+8)!.

Elementary beta integration now proves (1.4) for every polynomial, over the same field. This is not an inference from a finite prefix of moments. For the compact real model, D lies in[0,1], so equality of all polynomial moments and polynomial density in continuous functions identifies its actual probability law as Beta(1,8). For arbitrary norm-one forms it remains a rational polynomial-functional identity and is not presented as a positive real measure.

## 7. Polynomial CDF, Laurent identity, and direct geometric nonpower

Let v=1-(1-D)^8. Its derivative with respect toD is8(1-D)^7 and its endpoints are0,1. Formula (1.4) and the substitution x=1-(1-u)^8 give

lambda(v^n)=1/(n+1) for EVERY n>=0.                  (7.1)

For F=(1-t^-1)(1-v+vt), polynomial integration in the uniform variable v gives, at EVERY m>=1,

lambda_X(F^m)
 =((t-1)/t)^m integral_0^1(1+u(t-1))^m du
 =(t-1)^(m-1)(t^(m+1)-1)/((m+1)t^m).                (7.2)

The last expression is Laurent polynomial. Its two parts have exponent ranges1,...,m and -m,...,-1. Its constant term is therefore zero, while its coefficient of t^1 is(-1)^(m-1)/(m+1). This proves the exact pure-phase and fixed-unit-marker identities in both theorems. Every coefficient is defined over the stated base field, and normalized trace commutes with all further characteristic-zero extensions; the marked rational numbers remain nonzero.

Nonpower uses the actual product, independently of the moment recurrence. For either smooth geometrically integral variety M=X_c or Y_J, the divisor M_Kbar times{t=1} is prime. The first factor1-t^-1 of F has valuation1 there and the second factor restricts identically to1. Hence ord_(t=1)(F)=1. Units have valuation0 and proper r-th powers have valuation divisible by r>=2. Thus F is geometrically not a unit times any proper power.

There is also a direct support check. Expand

F=v t+(1-2v)-(1-v)t^-1.                               (7.3)

At the identity point C=e, both models have D=v=1. Thus the coefficient of t^1 is nonzero and the coefficient of t^0 is nonzero. A translate of dZ, d>=2, cannot contain both0 and1. Therefore F is radially primitive in the sense stated in the theorems. This support argument is independent of both the moment recurrence and the divisor valuation.

## 8. Compact positivity and formally real factor instability

For the compact rational model G_c(R) is compact and its real group points are Zariski dense. The quotient map is dominant, so a nonzero real regular function a has nonzero pullback to the compact group. It is nonzero on an open set of positive Haar mass. Uniqueness of the normalized trace then gives lambda_R(a^2)>0. This uses only the dense compact orbit, not a claim that every real point of an arbitrary form belongs to it.

The generic rational Gram argument extends this strict-square property to every formally real field. For rational linearly independent a_1,...,a_s the rational matrix (lambda_Q(a_i a_j)) is positive definite overR. Rational LDL elimination has positive rational diagonal pivots. Over an ordered formally real extensionK it remains positive definite, being a sum of positive rational multiples of squares after an invertible rational linear change. Every nonzero a in O(X_c,K) lies in theK-span of some such rational independent set, so lambda_K(a^2)!=0. A pure phase must have square trace zero; hence only0 is pure and the factor X_c,K is Mathieu. No whole-field embedding intoR is used.

The scalar split torus factor is Mathieu over every characteristic-zero field by the ordinary Laurent constant-term theorem [DvK] and the same-phase coefficient-field argument. Together with (7.2), this proves genuine two-factor instability for every formally real K in Theorem1. Theorem2 asserts product failure and nonpower but imports no compact positivity for its generally noncompact norm-one forms.

## 9. The compact beta law as a classical Shilov corollary

The compact Beta(1,8) law is consistent with, and follows independently from, the classical Shilov--Hardy norm formula. Here is the exact specialization, so that the scope of the attribution is explicit.

The Hermitification of the Euclidean Albert algebra has complex dimension27 and Jordan rank3. A primitive idempotent c=e1 is a rank-one tripotent, with B(c,c)=1. Let S be the full Shilov boundary, with its invariant probability. Its determinant-one part is the compact E6/F4 orbit. Equivalently S=U1 E6 e: in the Jordan spectral decomposition a maximal tripotent has all three singular values equal to1. Scalar multiplication by a unit complex number can set its determinant to1; Yokota's real compact orbit theorem then places it in E6 e. Conversely scalar multiples of that orbit are maximal tripotents. This is a comparison of compact real probability models, not an identification of their rational varieties.

The normalized probability on S is therefore the pushforward of product Haar on U1 times E6. The function D=w tau(w) is invariant under the extra scalar circle. Its moments on S consequently agree with those on the determinant-one quotient. No arbitrary phase or purity assertion is transported by this averaging.

In the notation of Upmeier [U], the Fischer kernel for homogeneous degree-j polynomials evaluated at c is w(z)^j/j!. By Lemma4.2 and its full proof, the irreducible kernel E^lambda(.,c) vanishes unless the partition lambda has length at most rank(c)=1. In degreej only lambda=(j,0,0) survives. Thus w^j belongs to that single polynomial summand, and its Fischer norm squared is j!.

The Shilov Hardy/Fischer norm relation displayed in the proof of [U,Lemma4.1] is

(p|p)_F=(dim/rank)_lambda (p|p)_S.

Here dim/rank=27/3=9 and (9)_(j,0,0)=(9)_j. It follows that

integral_S |w|^(2j)=j!/(9)_j,

which, by scalar-circle invariance, is exactly the compact diagonal moment sequence in (6.1). On[0,1] these moments determine the Beta(1,8) law. This classical comparison concerns that compact coordinate law. It does not claim that the whole explicit product theorem, the rational-form construction, or the general reduced-Albert norm-one statement occurs in [U]. The all-order proof in Sections4-7 remains complete if this comparison is omitted.

## 10. References and scope of the inputs

[Y] Ichiro Yokota, Exceptional Lie groups, arXiv:0902.0431v1 (2009), [exact primary PDF](sources/yokota-exceptional-lie-groups-0902.0431v1.pdf). The model inputs are Section2.1; Lemmas2.2.2-2.2.4 and their proofs, printed34-36; the full complex F4 automorphism group in Section2.13, including Lemma2.13.1 and Theorem2.13.2 with complete proofs, printed66-67; the compact cubic/Hermitian model and Lie decomposition in Sections3.1-3.2, printed68-70; full F4 stabilizer in Section3.7, printed82-83; complex diagonalization Proposition3.8.2 and the compact real orbit Theorem3.8.3, printed83-85; and the full complex E6 group Theorem3.14.1, printed104-105. Section3 supplies explicitly the extra diagonal rescaling that extends compact diagonalization to the entire norm-one hypersurface. The real compact set description is not used as an equation-scheme identification.

[U] Harald Upmeier, Stratified Hilbert Modules on Bounded Symmetric Domains, Complex Analysis and Operator Theory17,article74(2023), DOI[10.1007/s11785-023-01377-1](https://doi.org/10.1007/s11785-023-01377-1), [exact primary PDF](sources/upmeier-stratified-hilbert-2023.pdf). Lemma4.1 with its full proof and Lemma4.2 with its full proof, PDFpages18-20, supply only the compact Shilov moment comparison of Section9. They are not premises of the direct Lie-recurrence proofs; they also supply the optional independent compact-moment route in Section9.

[M] J. S. Milne, Algebraic Groups (2017), [author PDF](sources/milne-algebraic-groups-2017.pdf), Sections12l and22.42-22.43: the generic same-field normalized invariant trace and its scalar compatibility. The finite rational Gram argument is the universal one also recorded in [the classical compact-product proof](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3 and13.1. Only these generic lemmas are reused; no classical-family product result is a hypothesis here.

[DvK] J. J. Duistermaat and W. van der Kallen, Constant terms in powers of a Laurent polynomial, Indagationes Mathematicae9(1998),221-231, [primary PDF](sources/dvk-powers-1998.pdf). Its scalar Laurent theorem is used only for the individual split-torus factor in Theorem1; no coordinate-algebra-valued version is invoked.

[Long] C. D. Long, Counterexamples to the xz-Conjecture and the Mathieu Conjecture for SU(2), arXiv:2607.19012v1 (2026), Theorem2.1 and Proposition3.1. The polynomial-CDF pullback here uses the known interval--Laurent mechanism. Its complete relevant moment calculation is reproduced in (7.2), so no imported product-existence premise or priority claim is required.
