# Factorized homogeneous fields and an actual nonfinite A11 completion

The one-point axis case of KMFI120 does not force cyclic monodromy.
Here we classify an entire homogeneous architecture inside that case,
close it to exceptional11, and give a nonfinite full-rank A11 input
with a fully specified ordinary-A6 Keller completion.

The extra hypothesis is an actual polynomial factorization
\[
                  H=g-fp=F(p)Mp,                          \tag{1}
\]
where F is homogeneous of positive degree n and M is constant.
We assume r>=2, P=Dg-p df nilpotent, and generic rankDX=r.
No diagonalizability or monomial form is assumed in (1).
General homogeneous H without this factorization, nonhomogeneous
fields and the absolute exceptional11 action remain open.

## 1. Rigidity of the factorization

**KMFI121.** Over algebraically closed characteristic-zero constants,
M is invertible and has r distinct eigenvalues a_i.
There is a nonzero constant c, with no a_i=c, such that after a
constant linear source change and matching normalized target change,
\[
 F=\kappa\prod_{i=1}^r p_i^{b_i},\quad b_i\in\mathbb Z_{>0},
 \quad n=\sum_i b_i\ge r,\quad
 f=-cF,\quad g_i=(a_i-c)Fp_i,                               \tag{2}
\]
and exactly the constant matrix condition
\[
 \det(TI-[\operatorname{diag}(a_i)+a b^t])=(T-c)^r          \tag{3}
\]
remains. Conversely (2)--(3), with distinct nonzero a_i different
from c, produce a normalized polynomial input full rank on F!=0.
Its actual axis fiber is just the distinguished point.

First, DH=P-fI, so detDH=(-f)^r. On the other hand
DH=M(FI+p dF), giving detDH=(n+1)det(M)F^r.
Full rank excludes f0, so M is invertible and f=-cF, where
\[
                        c^r=(n+1)\det M\ne0.
\]
A rational function whose r-th power is constant is constant.

Next c is not an eigenvalue of M. Let A=M-cI.
On F!=0, a vector killed by A and dF is killed by DX.
If dimkerA>=2 such a vector exists. If dimkerA1, write
adjA=u v. Nilpotency gives
0=detP=F^(r-1)(dF u)(vMp).
Here vMp=c vp is not identically zero, so dF u0, again
contradicting full rank. Thus A is invertible.

It remains to prove that M has distinct eigenvalues.
The determinant lemma gives the resolvent identity
\[
 \frac{dF(TI-M)^{-1}Mp}{F}
           =1-\frac{(T-c)^r}{\det(TI-M)}.                  \tag{4}
\]
Its expansion at infinity makes every dF M^j p/F constant.
Thus F is an eigenfunction for the linear derivation associated
to every polynomial in M. If that matrix is nilpotent, the
derivation is locally nilpotent on the polynomial ring, so
its eigenvalue on nonzero F must be zero.

Apply this to the nilpotent parts of each generalized eigenspace
projector, which are polynomials in M. All higher-order poles
in the left side of (4) disappear: only the simple projector
terms remain. Since no eigenvalue equals c, the right side
has a pole at each eigenvalue of order equal to its algebraic
multiplicity. Every multiplicity is therefore1.
In eigen-coordinates, each p_i partial_i F is a constant times F.
A polynomial with this property is one monomial.
Its nonnegative integer exponents satisfy
\[
 b_i=-\frac{(a_i-c)^r}
              {a_i\prod_{j\ne i}(a_i-a_j)}.                \tag{5}
\]
None is zero, proving positivity and (2).
Equation(3) is the same determinant identity.

Conversely, on the dense torus P is similar to
F[diag(a_i)+a b^t-cI], hence nilpotent by (3).
The invertible diagonal matrix diag(a_i-c) receives a rank-one
perturbation, so P has one-dimensional kernel. A kernel vector
annihilated by dF would be zero. KMFI116 proves full rank;
polynomial continuation proves normalization everywhere.
Finally T_1=(I+FM)p can have no nonzero zero:
F would be nonzero, all p_i nonzero, and all distinct a_i
would have to equal -1/F. KMFI118 makes the fiber reduced,
so its whole ideal is the origin ideal.

## 2. Exact one-variable field and its critical values

Let h=gcd(b_i), m_i=b_i/h, W=product p_i^(m_i),
A=Y0 and C=product Y_i^(m_i). Then F=kappa W^h and
\[
 C=q_A(W)=W\prod_i(A+a_i\kappa W^h)^{m_i},\qquad
                       \deg_W q_A=n+1=N,                 \tag{6}
\]
\[
 p_i=Y_i/(A+a_i\kappa W^h),\qquad
                       \gamma=A+c\kappa W^h.             \tag{7}
\]
Since (m_i) is primitive, C and r-1 Laurent monomials form
birational coordinates on the target torus. Thus (6) is the
exact field pair with r-1 common spectators, not merely a slice.
The generic degree is N.

Differentiating (6) using (3) gives
\[
 q_A'(W)=(A+c\kappa W^h)^r
                     \prod_i(A+a_i\kappa W^h)^{m_i-1}.    \tag{8}
\]
The h roots of A+c*kappa*W^h have distinct nonzero critical
values: q_A(zeta W)=zeta q_A(W). Each therefore supplies one
single(r+1)-cycle. The remaining finite critical values are0.
At0 the branch cycles have m_i repeated h times; the zero W0
is unramified. These distinctions matter: two critical points
over0 do not supply separate transpositions.

If N is prime and n-r>=3, primitivity and
[Jones Corollary1.3](https://arxiv.org/pdf/1209.5169) give A_N
inside the relative group. The single(r+1)-cycle then gives S_N
for odd r. For even r all branch signs are even:
n is even, and the index at0 is n-hr. The relative group is A_N.
The discriminant below controls enlargement in the parameter A.

## 3. The complete degree11 obstruction inside (1)

**KMFI122.** Every degree11 input satisfying (1) has n10.
Its geometric action is S11 for odd r, A11 for even2<=r<=8,
and regularC11 forr10. Thus none gives exceptionalPSL2(11).
Rank1, outside the reconstruction argument above, has S11 by
the already-proved normalized pencil theorem.

For2<=r<=7 this follows from (8) and the preceding prime-degree
argument. Forr8 the9-cycle fixes2 points; Jones Theorem1.2
allows a proper primitive group only in degreeq+1 for a prime
powerq. Here q10 is impossible. Forr9, positivity of b_i and
sum b_i10 force one exponent2 and eight1, hence h1 and a
single branch transposition at0. Primitivity forces S11.
Forr10 all b_i1 and
q_A'=(A+c*kappa*W)^10, giving a binomial cyclic degree11 field.
It is absolutely forbidden for Keller completion by KMFI2.

There is also an exact arithmetic statement for even2<=r<=8,
over any coefficient field k where (1) is defined, even before
diagonalizing M. Put Z=F(p), B=F(Y_vec), and
\[
 Q_A(Z)=Z\prod_i(A+a_iZ)^{b_i}.
\]
This polynomial is defined over k: it is the coefficient expression
F((AI+ZM)p)/F(p), multiplied by Z.
Equation(7) also gives L=K(Z) and Q_A(Z)=B.
Let L0 be the constant leading coefficient of Q_A and
v_A=Q_A(-A/c). Then
\[
 \operatorname{Disc}_Z(Q_A(Z)-B)=
 -11\left[11^5 L0^5
       B^{(10-r)/2}(v_A-B)^{r/2}\right]^2.                 \tag{9}
\]
The geometric action is therefore A11 over the whole target,
not enlarged by moving A. Arithmetic is A11 iff -11 is a
square in k, otherwise S11, with exact constants k(sqrt(-11)).
For odd r, geometric S11 leaves no arithmetic enlargement.
The regular r10 case instead has its binomial cyclotomic constants.

These are exclusions of the entire factorization(1), including
non-diagonal constant M, repeated-eigenvalue attempts and arbitrary
homogeneous F. They do not make a general homogeneous H factor
as in(1). Equivalently, remove the greatest common polynomial
factor of the homogeneous coordinates of H. If the remaining
tuple has degree1, its matrix is invertible by full rank and
the entire result applies. A homogeneous exceptional11 candidate
must therefore have primitive projective coordinate degree>=2.
This is invariant under constant linear changes in these coordinates.

## 4. An explicit nonfinite full-rank A11 input

**KMFI123.** Let a be a root of the irreducible quartic
\[
                 R(a)=330a^4-165a^3+55a^2-11a+1.           \tag{10}
\]
Set
\[
 c1=2-4a,\quad c2=12a^2-6a+2,\quad
 c3=-30a^3+15a^2-5a+1,
\quad
 C=\begin{pmatrix}0&0&c3\\1&0&-c2\\0&1&c1\end{pmatrix},
\quad M=\operatorname{diag}(a,C).
\]
For p=(p1,p2,p3,p4), put
\[
 F=p1^7\det(p2I+p3C+p4C^2),\qquad
                         X=(-F,F(M-I)p).                  \tag{11}
\]
This has rational coefficients in Q(a), full rankDX4,
homogeneous matrix degree10, generic degree11 and geometric A11.
It is nonfinite and its axis fiber is exactly one point.
Its arithmetic group is S11, with constants Q(a,sqrt(-11)).

Here are exact coefficient and group certificates.
Let
\[
 Q_A(Z)=Z(A+aZ)^7
              (A^3+c1A^2Z+c2AZ^2+c3Z^3).
\]
Modulo R(a),
\[
 Q_A'(Z)=(A+Z)^4(A+aZ)^6.                                 \tag{12}
\]
The cubic factor is separable, has three distinct nonzero roots
in its reciprocal eigenvalue form, and avoids a and1.
For a direct reason, at A1 integrate (1+Z)^4(1+aZ)^6.
Its value at -1/a is -R(a)/(2310a^5).
Thus it has exactly a sevenfold zero there and a simple zero
at0. Its other three zeros are simple; the critical value at
-1 is nonzero, since otherwise roots of orders7 and5 would
exceed degree11. This proves every required separation.
The eigen-coordinate exponents are(7,1,1,1).
KMFI122 gives A11, and (9) gives the arithmetic test.
The quartic has Galois group S4, so Q(a) has no quadratic
subfield and cannot contain sqrt(-11).

Use p*=(1,1,0,0) with target X* and A0=-1.
The reduced polynomial
\[
 \Theta(Z)=\frac{Q_{-1}(Z)-Q_{-1}(1)}{(Z-1)^5}             \tag{13}
\]
has degree6, coprime to Z-1 and det(-I+ZM).
Choose a root beta and recover a by the linear subresultant
of R(a) and Theta(Z). Their resultant is an explicitly retained
irreducible degree24 polynomial over Q. Thus a belongs to Q(beta)
and [Q(beta):Q]=24.
Then
\[
 \gamma0=\beta-1,\qquad
                  p^0=(-I+\beta M)^{-1}X^*_{\rm vec}      \tag{14}
\]
is an exact regular second preimage of X*.
The cofactor equation proves F(p^0)=beta and the first target
equation; every denominator is a unit by the coprimality above.
No numerical marked root is used.

## 5. The actual ordinary-A6 completion and its normalization

Center at p*, use the regular Jordan frame T over Q(a), and
let the centered X be
\[
 X_0^c(w)=f(p*+Tw)-f(p*),
\]
\[
 X_{\rm vec}^c(w)=T^{-1}\{g(p*+Tw)-p* f(p*+Tw)
                         -(g(p*)-p*f(p*))\}.
\]
Set a=(4,5,6,7), e=(1,5,6,7,8), and
\[
 E_i(\gamma,u)=
 [X^c(\gamma^{a_j}u_j)+\gamma(1,\gamma^{a_j}u_j)]_i
                                                    /\gamma^{e_i}.
\]
The point(14) gives the exact zero theta0 of E.
Let L=DE(theta0)^-1.
The model contains every entry of theta0,T,T^-1,L, the finite
curve h(z) moduloz^8 and the lower-unitriangular A(z) from
KMFI103. All are serialized in the24-term power basis of beta.
Put
\[
 \phi(z,v)=(z,h(z)+A(z)v),\qquad
                    \Psi=\theta0+L\phi.
\]
The full map on ordinary coordinates(x,y0,...,y4) is
\[
 \mathcal F=(x\gamma,\ E_0(\Psi)/x,\ E_1(\Psi)/x^5,
                         E_2(\Psi)/x^6,E_3(\Psi)/x^7,
                         E_4(\Psi)/x^8),                  \tag{15}
\]
where (z,v_i)=(xy0,x^{e_i}y_i), and gamma is the first
component of Psi. These formulas and the finite coefficient
model specify every coordinate; no further equation must be solved.

The source automorphism has an explicit polynomial inverse.
The weighted first-stage and finite curve/flag checks prove
global polynomiality, including x0. KMFI104 proves the whole
constant Jacobian beta-1 and reconstructs the exact stabilized
field pair. Hence (15) has generic degree11 and geometric A11.
The degree24 polynomial remains irreducible over Q(sqrt(-11)),
so sqrt(-11) is not in Q(beta). Its arithmetic group is S11
and exact closure constants are Q(beta,sqrt(-11)).
The final six-variable determinant is not claimed to be expanded.

This is an actual nonfinite full-rank input and affine completion,
rather than a passed local group test. A11 was already realized;
no new abstract action or minimum-dimension claim is made.
It also shows that a one-point axis fiber need not have cyclic
generic group.

The geometric normalization boundary can be counted exactly.
Use split eigen-coordinates and exponents(7,1,1,1).
The sweep is finite away from the four target hyperplanes Yi0:
Q_A is monic up to a constant and (7) has unit denominators
when F(Y_vec)!=0. Over the visible branch Q_A(-A)=F(Y_vec),
the full inertia is one5-cycle and all unramified tails are present.

At the exponent7 target hyperplane, Z0 gives a retained
unramified residue-degree-one prime.
The root Z=-A/a has its ramification killed by the seventh-power
base change; its unramified residue extension has degree7 and
is retained. It is not seven separate residue-degree-one primes:
its Kummer unit contains each other target coordinate to exponent1.
The other three simple roots give three omitted unramified primes.

At each of the three exponent1 target hyperplanes, the root Z0
and that hyperplane's simple root are retained unramified primes.
The exponent7 root gives one omitted prime of index7 and residue
degree1; the other two simple roots are omitted unramified.
There are no other boundary divisors. Thus the original
normalization has precisely12 omitted primes: three index7
and nine unramified. KMFI108 makes its class group free of rank12.

For the specified higher-weight completion(15), the original
single5-cycle becomes identity on C0. All11 primes there split;
two are retained and nine unramified primes are omitted.
The visible index5 prime is also omitted after completion.
Consequently the geometric finite normalization has
\[
 \operatorname{Cl}(Z_{\mathcal F})=\mathbb Z^{22},\qquad
 [K_{Z_{\mathcal F}}]=(6,6,6,4,0,\ldots,0),                 \tag{16}
\]
with18 zero entries, in the omitted-prime basis.
This is after algebraic closure of Q(beta), not a claim about
its class group over the degree24 field itself.
It follows from actual boundary counts and the proved class
theorems, not from treating the axis sheet deficit as a divisor count.

## 6. Evidence and remaining scope

[Constructor](factorized_homogeneous.py) ·
[checker](verify_factorized_homogeneous.py) ·
[model with every finite A6 completion coefficient](factorized_homogeneous_model.json) ·
[receipt](factorized_homogeneous_receipt.json) ·
[primary group input](sources/factorized-homogeneous-manifest.json).

Exact checks cover the quartic/cubic identities and separations,
the degree24 marked point, inverse reconstruction, Jordan frame,
entire finite curve and flag residuals, actual polynomial source
inverse, discriminant square class, boundary factors and all
serialized coefficient data. Formal verification, independent
whole-proof replay and external review are not claimed.

The entire homogeneous class(1), including its one-point axis
case, cannot supply exceptionalPSL2(11). The completion is A11.
A general homogeneous H need not have a scalar-times-linear
factorization; genuinely unfactored homogeneous and nonhomogeneous
inputs, the absolute exceptional11 action and the full faithful
transitive Keller spectrum remain open.
