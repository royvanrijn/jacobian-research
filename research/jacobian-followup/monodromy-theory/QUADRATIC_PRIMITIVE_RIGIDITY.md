# Every primitive quadratic rank4 direction is excluded to full rank

The two remaining arbitrary-matrix quadratic families are now closed.
Two polynomial lemmas, one from a highest coefficient of the trace and
one from the leading normal block along a zero stratum, give the result.
Thus a homogeneous rank4 exceptional11 input must have primitive
direction degree at least3. This is a constructor-level obstruction;
the absolute exceptional11 action remains open.

## 1. A translation coefficient forbids a polynomial density

Work over a characteristic-zero field, with polynomial vector fields
on affine space. Let u be a nonzero linear form and a,b constant vectors,
with b nonzero, u(b)=0 and u(a)!=0. Suppose
\[
                         \partial_b V=u\,a.
\]
Let R be any polynomial independent of the b direction.

**KMFI155.** No nonzero polynomial F satisfies
\[
                           \operatorname{div}(FV)=R F.    \tag{1}
\]
Neither F nor V is assumed homogeneous, and no degree bound is used.

Choose constant linear coordinates retaining u, with
b=partial_z and a=c*partial_u, c=u(a)!=0.
This is possible because a,b are independent; choose z(a)=0,z(b)=1
and the other coordinates invariant under both directions.
Then V=c*u*z*e_u+V0, where V0 is independent of z, and
divV=c*z+divV0.
For F=J_m*z^m+lower powers, J_m nonzero, the coefficient of z^(m+1)
in(1) is
\[
                              c(u\partial_u J_m+J_m)=0.
\]
Every monomial has nonnegative u degree, so no polynomial can have
u*partial_u eigenvalue -1. This proves the lemma.

For a normalized homogeneous input H=FV of primitive quadratic degree,
h=lambda*u and the trace identity is
\[
                  \operatorname{div}(FV)=r\lambda uF.     \tag{2}
\]
Each canonical form V=A*K of KMFI143 satisfies
partial_z V=u*A*e_z. Consequently normalization alone forces
\[
                                  A03=0                  \tag{3}
\]
in every one of the four families, regardless of whether F depends on z.
This gives a shorter route through the previously closed single-square
family; no cone or scalar-factor argument is needed for(3).

The polynomial qualification is essential. For V=(u*z,0) in coordinates
(u,z), F=1/u and R=0 solve(1), with u(a)=1.
This rational negative-degree control is not a polynomial multiplier.

## 2. A zero stratum cannot have an invertible normal derivative

Use linear coordinates p=(q,t), with q=(q1,...,qc), and let Z={q=0}.
Suppose polynomial V,h satisfy V|Z=0,h|Z=0, and
\[
                B(t)=D_q V_q(0,t)
\]
is invertible over the rational function field k(t).

**KMFI156.** For every nonzero polynomial F, the rational matrix
\[
                     N=DV+V\,d\log F-hI                  \tag{4}
\]
cannot be nilpotent identically. No rank or homogeneity hypothesis
on the full input is required.

Let m>=0 be the smallest total q degree occurring in F and write
\[
 F(\epsilon q,t)=\epsilon^m(J(q,t)+O(\epsilon)),
\]
where J is nonzero and homogeneous of degree m in q.
Since V vanishes on Z, V_q(epsilon*q,t)=epsilon*Bq+O(epsilon^2)
and V_t is likewise O(epsilon).
Conjugate N(epsilon*q,t) by T=diag(epsilon*I_c,I).
The rational limit exists over k(q,t), and has block form
\[
 \lim_{\epsilon\to0}T^{-1}N(\epsilon q,t)T
       =\begin{pmatrix}B+(Bq)d_q\log J&*\\0&0\end{pmatrix}. \tag{5}
\]
In detail, the tangential-normal block acquires an extra epsilon,
while the tangential-tangential block tends to zero; the
normal-tangential block is finite and need not vanish.
The h term tends to zero by h|Z=0.

If N were nilpotent, so would be(5) and its normal block.
But the determinant lemma and Euler's identity give
\[
        \det(B+(Bq)d_q\log J)
                   =\det B(1+d_q\log J\cdot q)
                   =(m+1)\det B\ne0.                     \tag{6}
\]
This is the contradiction. The leading multiplicity m is any
nonnegative integer, so the obstruction has no vanishing-order bound.
All limits are rational identities in epsilon, not an analytic
properness assumption.

Invertibility of the normal block matters: the linear nilpotent field
V=(q2,0,0), F=1,h=0 vanishes on q1=q2=0 and has nilpotent DV.
Polynomiality matters too: in codimension1, V=(q,0),F=1/q,h=0
makes(4) zero, with the forbidden order m=-1.

## 3. Both single-shear families are closed

Work over algebraically closed characteristic-zero constants in rank4.
Let V=A*K_D with
\[
                  K_D=(u^2,ux,uy,uz+xy),\quad A\in GL_4.
\]
Suppose F gives a normalized full-rank homogeneous input.
KMFI150 says u does not divide F, and KMFI155 gives A03=0.
Hence V_u=uL with L linear. At the generic point of u0, dlogF is regular
and row u of N=DV+V*dlogF-lambda*u*I is
(L|u0,0,0,0). Nilpotency forces L|u0=0.
Thus V_u=a*u^2, with a!=0 since detDV!=0.
This is the matched primitive case excluded by KMFI141.

Every arbitrary output matrix and every homogeneous scalar degree
in the product-shear family is therefore excluded to full rank.
Exactly the same argument applies to K_C=(u^2,ux,uy,uz+x^2), giving
a shorter proof of KMFI154 without its nonlinear-divisor geometry.

## 4. The mixed chain is closed

Now let
\[
                  V=A*K_B,\quad
                  K_B=(u^2,ux,uy+x^2,uz+xy),\quad A\in GL_4.
\]
Again u does not divide F and A03=0.
If A02=0, the same regular critical-boundary row gives a matched square,
already excluded. Thus put d=A02!=0.

The linear plane Z={u=x=0} is a zero stratum of V and h=lambda*u.
Its normal derivative in coordinates(u,x) is
\[
 B(y,z)=\begin{pmatrix}
 A02*y+A03*z&A03*y\\
 A12*y+A13*z&A13*y
 \end{pmatrix},\qquad
 \det B=(A02*A13-A03*A12)y^2.                              \tag{7}
\]
After A03=0, this is d*A13*y^2.
KMFI156 therefore forces
\[
                                  A13=0.                 \tag{8}
\]

It remains to contradict the trace at A03=A13=0,d!=0.
Give y,z total degree1 and u,x degree0, and let F_m be the nonzero
highest homogeneous part of F for that grading.
All four coordinates of V have y,z degree at most1.
The only terms raising that degree in V(F) are
\[
                    u y(d\,\partial_u+A12\,\partial_x)F_m.
\]
The degree1 part of divV is d*y; the h term has degree0.
Thus the degree m+1 trace identity is exactly
\[
       u y(d\,\partial_u+A12\,\partial_x)F_m+d yF_m=0.      \tag{9}
\]
Change x to x-(A12/d)u, retaining u,y,z.
Then(9) gives u*partial_u F_m=-F_m, impossible for a polynomial.
Equivalently apply KMFI155 to the leading vector field
W=u*y*(d,A12,0,0), b=e_y, R=0.
This excludes every remaining parameter and scalar degree in K_B.

## 5. Complete quadratic primitive exclusion and precise survivors

**KMFI157.** No homogeneous primitive quadratic direction V in rank4
admits a nonzero homogeneous polynomial multiplier giving a normalized
input of generic rankDX4.

KMFI125 forces detDV=delta*u^4 with delta!=0.
KMFI142--143 give the four exhaustive primitive families under allowed
constant conjugacy. KMFI152 closes K_A, and Sections3--4 above close K_C,K_D,K_B.
The older KMFI154 single-square proof is retained but not needed for
this shorter closure route.
All arbitrary output matrices and every homogeneous scalar degree
are covered. The scalar-factor degree bound3 is not used in the new
two-family proof.

If a determinant-compatible homogeneous quadratic presentation V is
not primitive, it has form V=u*M*p by KMFI143.
Then H=F*V=G*M*p with G=uF. Every full-rank normalized survivor is
exactly among the KMFI121 data with G divisible by u, and conversely
such data give this presentation by F=G/u.
The whole generic field is the explicit one-variable polynomial
q_A(W) of KMFI121, with its common independent spectators.
No unproved all-degree group classification is asserted from that
field formula. In degree11 and rank4, KMFI122 gives geometric A11,
with its stated arithmetic constant-field test.
The regular C11 and odd-rank S11 cases belong to other ranks.

Consequently a homogeneous rank4 exceptional11 normalized input
must have primitive direction degree>=3. The arbitrary cubic-and-higher
primitive directions and genuinely nonhomogeneous fields remain open.
This does not exclude exceptional11 for unrelated Keller architectures.

The full-rank boundary is real. At lambda=s+2,
A=diag(1,lambda,lambda,lambda), F=u^s with s>=1 give normalized
rank3 mixed-chain inputs and rank2 product-shear inputs, all with
whole sweep Jacobian gamma^4. Exact examples in both families are
retained. The earlier quadratic-field single-square rank3 example
also remains valid. These are not new full-rank inputs or action claims.

## 6. Evidence and the full action objective

[Helper](quadratic_primitive.py) ·
[checker](verify_quadratic_primitive.py) ·
[model](quadratic_primitive_model.json) ·
[receipt](quadratic_primitive_receipt.json).

The unbounded polynomial lemmas and complete family proof are
supplemented by whole coordinate-frame identities, normal-block
limits with tangential terms retained, determinants for arbitrary
leading order, all-parameter mixed-chain normal and trace identities,
and actual normalized lower-rank controls with entire sweep Jacobians.
Formal verification, independent whole-proof replay and external
review are not claimed.

The previous goal turn was progress: KMFI153--154 closed the single-square
family. This turn closes the two remaining families and hence the full
primitive quadratic rank4 direction class.
The full faithful-transitive Keller monodromy spectrum, actual
exceptional11 occurrence or an absolute obstruction, higher primitive
directions and nonhomogeneous inputs remain unresolved.
The user's complete action-level objective remains active.
