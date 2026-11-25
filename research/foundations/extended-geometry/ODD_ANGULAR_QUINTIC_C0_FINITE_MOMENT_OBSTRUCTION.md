# Sharp finite-twelve obstruction for the C=0 odd-angular quintic

Let \(K\) be any field of characteristic zero and put
\[
 {\cal A}_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(F)=\operatorname{CT}_x I_t
       F\!\left(x,\frac{1-t^2}{2x},t\right).
\]
Here \(I_t(t^{2j})=1/(2j+1)\) and \(I_t(t^{2j+1})=0\); no analytic
interpretation or ordering of \(K\) is required. A phase is pure if
all its positive power moments vanish.

Consider the normalized seven-monomial phase
\[
 f=a_1y+a_3y^3+y^5+b_1yt^2+\lambda yt^4
                      +xt^2(q+t^2),\qquad M_m=\Lambda(f^m).        \tag{1}
\]
Its \(y^3t^2\) coefficient is zero and its \(y^5,xt^4\) coefficients
are one. The five parameters are in the order
\((a_1,a_3,b_1,q,\lambda)\).

**Theorem 1 (Finite-twelve obstruction).**
For every \(K\) and every parameter point in (1), at least one of
\[
                  M_2,M_4,M_6,M_8,M_{10},M_{12}
\]
is nonzero. In particular no member is pure.

**Corollary 2 (Nonzero coefficient normalization).**
The same finite-twelve obstruction holds for
\[
 a_1y+a_3y^3+D y^5+b_1yt^2+Ayt^4+xt^2(q+Bt^2),
                         \qquad BD\ne0.                            \tag{2}
\]
Choose \(\rho^6=D/B\) in an algebraic extension of \(K\) and
\(\mu=1/(\rho B)\). The actual quadric automorphism
\(x\mapsto\rho x,\ y\mapsto\rho^{-1}y,\ t\mapsto t\)
preserves \(\Lambda\), as is immediate from its Laurent constant-term
definition. Multiplication by \(\mu\) therefore scales the \(m\)-th
moment by the nonzero factor \(\mu^m\). The transformed coefficients of
\(xt^4,y^5\) are both one. Apply Theorem1 over the extension and use
faithfulness of the field map. No sixth root in \(K\) itself is needed.

**Theorem 3 (Sharpness across characteristic-zero fields).**
There is a genuine quintic of the form (1) over \(\mathbb Q_{103}\)
whose first eleven positive moments vanish and whose twelfth moment
does not vanish. There is also such a phase over a number field.
Consequently twelve is the least single uniform cutoff valid for (1)
over all characteristic-zero coefficient fields.

This is not a claim that twelve is optimal over each individual field,
nor an assertion of a rational or real realization. These statements
do not classify the full even-angular-addition space, prove an
eventual-purity obstruction, or settle a global degree-five problem.

## 1. Fresh native moments and the two pivots

For nonnegative integers \(i,j,k\), the native monomial trace vanishes
unless \(i=j=r,\ k=2s\). In that case
\[
 \Lambda(x^ry^rt^{2s})
       =\frac{2^r r!(2s)!(r+s)!}{s!(2r+2s+1)!}.                     \tag{3}
\]
Substitution \(xy=(1-t^2)/2\) and the finite binomial integral prove
this formula. Applying the multinomial theorem literally to the seven
terms in (1) reconstructs \(M_1,\ldots,M_{12}\) over \(\mathbb Q\).
No moment polynomial from a family containing \(y^3t^2\) is used.

More explicitly, let \(n_1,n_3,n_5,n_b,n_\lambda\) count the five
negative-angular terms. Angular balance forces
\[
 n_1+2n_3+3n_5+n_b+n_\lambda=m/2.
\]
Their total negative winding is
\(r=n_1+3n_3+5n_5+n_b+n_\lambda\).
Split the \(r\) positive-angular factors between \(qxt^2\) and
\(xt^4\), and apply (3) with the resulting even \(t\)-exponent.
This finite enumeration is the checker's native moment generator.
The antipode \((x,y,t)\mapsto(-x,-y,t)\) preserves \(\Lambda\) and
takes \(f\) to \(-f\), so all odd moments vanish identically.

Put
\[
 \begin{split}
 P_1&=7q+3,\\
 P_3&=1615q^3+2261q^2+1197q+231,\\
 P_5&=930465q^5+2225025q^4+2314026q^3\\
    &\hspace{12mm}+1285570q^2+376805q+46189.
 \end{split}                                                       \tag{4}
\]
The native coefficients of \(a_1\) in \(M_2\) and of \(a_3\) in
\(M_4\) are \(2P_1/105\) and \(8P_3/4849845\).
The exact native equation \(M_2=0\) is
\[
       231a_1q+99a_1+99b_1q+55b_1
                           +55\lambda q+35\lambda=0.               \tag{5}
\]
Off \(P_1P_3=0\), solve \(a_1\) and then \(a_3\) uniquely.
The retained rational functions have denominator factors only
\(P_1,P_3\); their native substitutions are checked by exact
cross-multiplication.

After this substitution write \(M_m=N_m/D_m\), with \(N_m\)
primitive integral and \(D_m\) supported only on \(P_1,P_3\).
The displayed \(P_5\) is not an eliminated-variable pivot in this
seven-term family.

## 2. Generic nonzero lambda: exact small resultants

Suppose \(\lambda\ne0\), and set
\[
                    h=b_1/\lambda,\qquad z=\lambda^3.
\]
For \(m=6,8,10,12\), each exponent sum of \(b_1,\lambda\) in \(N_m\)
is congruent to \(m/2\) modulo3. Thus
\[
 N_m(h\lambda,q,\lambda)
        =\lambda^{\epsilon_m}G_m(h,q,\lambda^3),
              \qquad \epsilon_m=(m/2)\bmod3.                       \tag{6}
\]
The exponent map is injective and divides only a power of the
assumed nonzero \(\lambda\). The respective degrees in \((h,q,z)\) are
\[
                 (3,11,1),\ (4,16,1),\
                 (5,18,1),\ (6,25,2).
\]

For \(m=8,10,12\), compute the full polynomial resultants
\(\operatorname{Res}_z(G_6,G_m)\).
Their complete \(q\)-contents, up to nonzero rational constants, are
\[
                    P_1^3P_3,\quad P_1^3P_3,\quad P_1^6P_3^2.
\]
Dividing these known units on the two-pivot chart and clearing only
constant contents gives primitive integral polynomials \(R_8,R_{10},R_{12}\)
of \((h,q)\)-degrees
\[
                         (4,19),\quad(5,21),\quad(6,31).
\]
The checker recomputes the exact \(z\)-resultants and their complete
contents, and verifies these polynomial identities over \(\mathbb Q\).

Two further exact resultants have the factorizations
\[
 \begin{split}
 \operatorname{Res}_h(R_8,R_{10})
       &=c_1P_1^{20}P_3^8P_5^3 F_{120},\\
 \operatorname{Res}_h(R_8,R_{12})
       &=c_2P_1^{24}P_3^{12}P_5^6 F_{148},
 \end{split}                                                       \tag{7}
\]
where \(c_1,c_2\in\mathbb Q^\times\) and the subscripts are the
degrees of the primitive integral residual polynomials.
The full resultant degrees are179 and238. The checker freshly
recomputes both over \(\mathbb Q\), verifies all divisions and
the residual coefficient vectors, and checks the nonzero rational
scalars by exact polynomial comparison. No CRT reconstruction is needed.

Both residual leading coefficients are nonzero modulo the verified
prime32003. The retained complete univariate coefficient vectors
give
\[
               uF_{120}+vF_{148}=1\pmod {32003}.                    \tag{8}
\]
Gauss's lemma therefore proves their coprimeness over \(\mathbb Q\):
a positive-degree primitive rational common factor, and its
quotients in the primitive integer polynomials, would retain nonzero
leading coefficients modulo32003 and contradict (8).
This is a univariate good-reduction argument with full degrees,
not an inference from a multivariate modular unit ideal.

At any actual common zero of the four \(G_m\), all three first
resultants and both outer resultants vanish. The full Sylvester
identity has this necessary vanishing property even after a leading
coefficient drops; none is inverted here. If
\(\lambda P_1P_3P_5\ne0\), (7) forces both coprime residuals to vanish,
a contradiction. The generic nonzero-lambda sector is excluded
through moment12.

## 3. Every zero-lambda point off P1P3

Set \(\lambda=0\). If \(b_1\ne0\), put \(z=b_1^3\) and divide
\(N_m\) only by its explicit nonzero factor
\(b_1^{(m/2)\bmod3}\). Clearing rational constant contents gives
\(F_6,F_8,F_{10},F_{12}\). The first three are linear in \(z\).
Write \(F_6=az+b,\ F_8=cz+d,\ F_{10}=ez+g\), with coefficients
in \(\mathbb Z[q]\). Then
\[
       -cF_6+aF_8=ad-cb=:S_8,\qquad
       -eF_6+aF_{10}=ag-eb=:S_{10}.                                 \tag{9}
\]
The retained exact rational coefficients satisfy
\[
                    uS_8+vS_{10}=\operatorname{monic}(P_1^3P_3).
                                                                        \tag{10}
\]
All leading-coefficient-zero cases remain: no coefficient of these
linear polynomials has been divided out. On \(P_1P_3\ne0\) the right
side is nonzero, excluding this case.

For \(b_1=0\), the two scalar equations \(F_6(0,q)=F_{12}(0,q)=0\)
are necessary. A second complete rational coefficient identity gives
\[
      u_0F_6(0,q)+v_0F_{12}(0,q)
                         =\operatorname{monic}(P_1^3P_3).          \tag{11}
\]
This case is impossible as well. The checker multiplies the full
coefficient vectors in (9)--(11). No \(P_5\) factor is inverted:
this proves the entire zero-lambda two-pivot sector, including all
P5 roots, using moments through12.

## 4. All P1 and P3 fibers

At \(q=-3/7\), the native \(M_2=0\) gives
\(b_1=-10\lambda/11\). Native \(M_4=0\), with a nonzero constant
pivot, gives
\[
 a_3=-\frac{38437}{264}a_1^2
       +\frac{16625}{484}a_1\lambda
       -\frac{26215}{10648}\lambda^2.                              \tag{12}
\]
Let \(E_6,E_8,E_{10},E_{12}\) be the remaining native polynomials
after these substitutions. The retained four rational polynomial
cofactors satisfy
\[
             C_6E_6+C_8E_8+C_{10}E_{10}+C_{12}E_{12}=1
                    \quad\text{in }\mathbb Q[a_1,\lambda].         \tag{13}
\]
Hence all P1 points are excluded, with every \(\lambda\) retained.

Modulo the full polynomial \(P_3\), \(P_1\) is a unit. Solve \(a_1\)
from native \(M_2\), keep \(a_3,b_1,\lambda\) free, and retain
native \(M_4\) as a compatibility equation. In particular the zero
\(a_3\) pivot is not inverted. The supplied five cofactors satisfy
\[
       C_4E_4+C_6E_6+C_8E_8+C_{10}E_{10}+C_{12}E_{12}
                       =1\pmod {P_3}.                             \tag{14}
\]
This is a full identity in \(\mathbb Q[a_3,b_1,\lambda,q]/(P_3)\),
verified by direct polynomial multiplication. It specializes at
every geometric P3 root; no choice of an embedding or field
irreducibility assumption occurs. All P3 points are excluded.

## 5. The remaining P5 points

The polynomials \(P_1,P_3\) are each coprime to \(P_5\), as exact
univariate Euclid verifies. Therefore the two native substitutions
are valid at every P5 root. Zero lambda was already excluded in
Section3. For nonzero lambda, reduce the exact scaled equations
of (6) modulo the full \(P_5\). They satisfy
\[
                         G_6=zA(h,q).
\]
Since \(z=\lambda^3\ne0\), an actual zero requires \(A=0\).
The retained integer polynomial cofactors give the directly checked
identity
\[
          I_AA+I_8G_8+I_{10}G_{10}
                        =\kappa\pmod {P_5},
                    \qquad \kappa\in\mathbb Z\setminus\{0\}.        \tag{15}
\]
Its coefficients have only rational constant denominators before
clearing. No \(h\)-pivot, \(z\)-shift, leading coefficient or branch
discriminant is inverted. Thus \(G_6=G_8=G_{10}=0\) is impossible,
at every P5 root and every \(h\), through native moment10.

Every characteristic-zero parameter point lies on P1, P3, P5,
or their complement. On the complement lambda is zero or nonzero.
Sections2--5 cover all these cases and prove Theorem1. Corollary2
was proved above by the actual coefficient normalization.

## 6. A smooth native prefix-eleven point

The checker reconstructs native moments and their derivatives directly
from (3). Every coefficient through \(M_{12}\) is103-integral.
At the point
\[
             (a_1,a_3,b_1,q,\lambda)=(83,26,73,3,29)
                                  \pmod {103},
\]
it verifies
\[
                M_2=M_4=M_6=M_8=M_{10}=0,\qquad M_{12}=11.
                                                                         \tag{16}
\]
The full native Jacobian, with this parameter order, is
\[
 \begin{pmatrix}
24&0&45&93&56\\
91&78&67&68&81\\
95&59&58&50&91\\
88&86&88&41&46\\
48&36&92&10&39
 \end{pmatrix},
                       \qquad\det=94\pmod {103}.                   \tag{17}
\]
This is not merely the Jacobian of generic reduced equations.
The checker also verifies one literal Newton correction modulo103 squared.

Put \(F=(M_2,M_4,M_6,M_8,M_{10})\) and start at the integer lift \(v_0\)
of (16). Its Jacobian is invertible over \(\mathbb Z_{103}\).
Define
\[
                     v_{n+1}=v_n-J_F(v_n)^{-1}F(v_n).
\]
Inductively \(v_{n+1}-v_n\) is divisible by \(103^{2^n}\).
Polynomial Taylor expansion has a remainder quadratic in this
correction and makes \(F(v_{n+1})\) divisible by \(103^{2^{n+1}}\).
No factorial division is used. The Jacobian remains invertible
modulo103; the sequence is Cauchy, so completeness gives a
\(\mathbb Z_{103}\) zero \(v\) of all five equations with the
same residue. The twelfth moment remains a unit with residue11.
Parity gives all odd moments zero, so its first eleven moments vanish.

This is a smooth-point lifting proof. A bare proper modular ideal
would not justify a characteristic-zero point.

## 7. A number-field witness and optimality

Let \(W=\lambda P_1P_3P_5M_{12}\). The point in Section6 has \(W\ne0\);
the pivots are units at its residue. Thus
\[
 {\cal B}=\mathbb Q[a_1,a_3,b_1,q,\lambda,u]/
            (M_2,M_4,M_6,M_8,M_{10},uW-1)                            \tag{18}
\]
is a nonzero finite-type rational algebra, since it maps to
\(\mathbb Q_{103}\).

For completeness, a nonzero finite-type rational algebra has a maximal
residue field finite over \(\mathbb Q\). To see the needed field lemma,
suppose a field finitely generated as a rational algebra had a nonempty
transcendence basis among its generators. Clear the finitely many
denominators in the equations of the remaining algebraic generators.
The field is then integral over \(\mathbb Q[t_1,\ldots,t_r,1/d]\).
A subring over which a field is integral must itself be a field:
the integral equation of the inverse of a nonzero element expresses
that inverse in the subring. But choose a rational \(c\) such that
\(t_1-c\) does not divide \(d\); it stays a nonunit in that localization,
a contradiction. The residue field is consequently algebraic and
finitely generated over \(\mathbb Q\), hence finite.

A maximal residue of (18) therefore gives some number-field parameter
point with first eleven moments zero and \(M_{12}\ne0\). The coefficient
of \(y^5\) is one, whose top Laurent monomial has angular winding-5
and cannot cancel against another term of (1); the phase is a genuine
quintic. No explicit coordinates, rationality, reality, or algebraicity
of the particular103-adic lift is asserted.

Together with Theorem1, either characteristic-zero witness proves
the uniform optimality statement in Theorem3.

## 8. Portable exact evidence and limits

The retained artifacts are the
[portable checker](certificates/verify_odd_angular_quintic_C0_finite12_certificate.py),
[compressed coefficient input](certificates/odd-angular-quintic-C0-finite12-certificate.json.gz),
and [completed verification receipt](certificates/odd-angular-quintic-C0-finite12-verification.json).
From the repository root the replay command is:

    python research/foundations/extended-geometry/certificates/verify_odd_angular_quintic_C0_finite12_certificate.py --output NEW_RECEIPT.json

The checker reconstructs all native moments through12 rather than
loading raw moment polynomials. It verifies both triangular rational
solutions and their complete denominator support, the small exact
resultants and all contents/divisions, the complete modular univariate
Bezout identity with preserved leading degrees, both zero-lambda
identities, all exceptional native cofactor products, and the full
native smooth point/Jacobian. Sixteen input-mutation controls and four
additional cofactor corruption checks reject changed arithmetic.
An exact compressed-input SHA256 guard rejects changed input bytes;
an existing receipt is not overwritten.

Runtime needs Python's standard library and python-FLINT, without
Singular, SymPy, network access, private-directory paths, raw-moment
input, or a CRT computation. Exact finite-field calculations enter
only the univariate good-reduction argument and the native smooth
point; their characteristic-zero interpretations are proved above.

The coefficient certificate is a finite proof for this specific
seven-term chart and its \(BD\ne0\) normalization. It neither asserts
eventual moment vanishing nor excludes phases with extra angular-even
terms or arbitrary quintic support. Ordinary-prime results, prime-face
parameter rigidity and a global purity-recognition algorithm are not
premises.
