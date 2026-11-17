# The entire arbitrary-matrix triangular cubic chain is excluded to full rank

The cubic family left open by KMFI161 is now closed for every invertible
output matrix and every homogeneous scalar degree. A last-coordinate
feedback identity supplies a general multiplicity test. In this family,
it forces the output matrix to be lower triangular, after which the
all-degree triangular rigidity theorem applies.

This closes one full cubic architecture. Arbitrary cubic directions,
nonhomogeneous inputs and the absolute exceptional11 action remain open.

## 1. A last-coordinate feedback multiplicity test

Use coordinates p=(q,z) in rank r>=2 and let
\[
                V=(v(q),\,c(q)z+R(q)),\qquad h=h(q).
\]
Let F be a nonzero polynomial, eta=dlogF, A=DV, and suppose
N=A+V*eta-hI is nilpotent. Assume the Krylov matrix
(V,AV,...,A^(r-1)V) is invertible and the reconstructed
X=(-Fh,F(V-hp)) has generic rank r.
For homogeneous normalized full-rank inputs the cyclicity hypothesis
is supplied by KMFI125.

Put D=cI-Dv and
\[
 W=\det(D)(cz+R)+(z\,dc+dR)\operatorname{adj}(D)v=Kz+L,
\]
so
\[
                 K=c\det D+dc\,\operatorname{adj}(D)v.     \tag{1}
\]

**KMFI162.** Necessarily c-h is nonzero, K is nonzero, and for a positive
integer nu=deg_zF,
\[
                      \nu K=-(c-h)^r.                    \tag{2}
\]
Over k(q), F is a nonzero scalar multiple of
(z+L/K)^nu. This is an exact necessary polynomial-multiplicity test;
it is not a sufficient normalization or generic-group criterion.

At T=c the adjugate vector adj(TI-A)V has all head components zero
and last component W. Cyclicity makes W nonzero: since c is a root of
chi_A, adj(cI-A)=Q(A), with Q(T)=chi_A(T)/(T-c), a monic polynomial
of degree r-1. If Q(A)V were zero, the Krylov columns would be dependent.
This remains valid at repeated eigenvalues; no derivative of chi
or simple-root assumption is used.

The feedback determinant lemma at T=c now gives
\[
                          -\eta_z W=(c-h)^r.              \tag{3}
\]
If c=h identically, then eta_z=0. Both F and X are independent of z,
contradicting full rank. Thus c-h is nonzero and eta_z is nonzero.
If K=0, equation(3) makes eta_z a nonzero constant in z.
A polynomial logarithmic derivative tends to zero at infinity, so
that case is impossible.
For K!=0, its 1/z coefficient at infinity is deg_zF, a positive integer.
Comparing(3) gives(2), and then integrating the one-variable polynomial
identity gives the stated single-root form over k(q).
The head-dependent scalar can be rational; no global polynomial
factorization without denominator clearing is asserted.

An actual normalized rank2 full-rank control is
V=(u,omega*z), F=u*z, h=1+omega, omega^2-omega+1=0.
Here K=omega*(omega-1)=-1 and nu=1, so(2) holds.
Its entire sweep Jacobian is gamma^2. It is an existing scalar-times-linear
input, not a new occurrence action.

## 2. The last quadratic-in-y coupling vanishes

Consider the primitive cubic direction
\[
 K_3=(u^3,u^2x,u^2y+ux^2,u^2z+xy^2),\qquad V=A K_3,\quad A\in GL_4.
\]
A hypothetical normalized full-rank homogeneous multiplier satisfies
the KMFI161 restrictions A02=A03=A13=0 and h=lambda*u^2.
Write the remaining matrix as
\[
 A=\begin{pmatrix}
 a&b&0&0\\d&e&f&0\\g&h&i&j\\k&l&m&n
 \end{pmatrix}.                                          \tag{4}
\]
The letter h in this matrix is a constant entry, distinct from the
polynomial feedback shift lambda*u^2.

First j=0. Suppose j!=0 and make the allowed simultaneous constant
source/output change Z=z-(n/j)y.
The y component has leading term j*x*y^2.
The Z component has no quadratic y term and no Z term, and its
coefficient of u^2*y is r=m-n*i/j.
The x component's coefficient of u^2*y is f.
The divergence is
\[
 (3a+e+i+n)u^2+2(b+f)ux+2jxy.
\]
For F=J(u,x,Z)y^t+lower y powers, J nonzero, the coefficient of y^(t+1)
in div(FV)=4lambda*u^2*F is
\[
           u^2(f\partial_x+r\partial_Z)J=-j(t+2)xJ.        \tag{5}
\]
If the constant derivative on the left is zero, this is already
impossible. Otherwise regard u as a coefficient in its rational
function field and straighten that derivative to partial_w.
A constant derivative of a nonzero polynomial has lower w degree,
while multiplication by a nonzero polynomial cofactor has degree
at least the original degree. Equation(5) is impossible.
This covers f=0, r=0 and every parameter coincidence, and proves j=0.

## 3. Last-coordinate polynomiality removes the head coupling

Now the last coordinate is V_z=n*u^2*z+R(u,x,y).
Invertibility of A gives n!=0. Set
\[
 M_3=\begin{pmatrix}a&b&0\\d&e&f\\g&h&i\end{pmatrix},\qquad
 v=u\,M_3(u^2,ux,uy+x^2)^t,\quad c=n u^2.
\]
The coefficient K in(1) has the entire identity
\[
 K=n\,u^8\det(nI-M_3)
       +2fn^2(a-n)u^7x+bfn^2u^6x^2.                     \tag{6}
\]
All last-row constants k,l,m remain arbitrary; the identity is independent
of them, not a specialization omitting their effects on W.

By KMFI162, c-h=(n-lambda)u^2 is nonzero and
K=-(n-lambda)^4*u^8/nu for an integer nu>0.
The two non-u^8 coefficients of(6) force
\[
                              f(a-n)=0,\qquad bf=0.
\]
If f!=0, then a=n,b=0. But M3 now has first row(n,0,0),
so det(nI-M3)=0, and(6) makes K=0 identically.
That is the separately rejected branch of KMFI162.
Consequently
\[
                                      f=0.               \tag{7}
\]
Cases lambda=n or K=0 are not divided away; their full-rank rejection
was proved in Section1.

## 4. A weighted leading trace removes the remaining unmatched entry

With f=j=0, give y weight1,z weight2 and u,x weight0.
Let F_top be the nonzero part of F of highest such weight.
The leading trace equation is
\[
 \left\{u^2D_B+i u^2y\partial_y+
 (n u^2z+nxy^2)\partial_z+
 [(3a+e+i+n-4\lambda)u^2+2bux]\right\}F_{\rm top}=0,
                                                               \tag{8}
\]
where D_B=(a*u+b*x)partial_u+(d*u+e*x)partial_x.
All lower-weight terms of the vector field drop out exactly.

Let the highest z monomial of F_top be J(u,x)y^s*z^t, J nonzero.
The xy^2*partial_z term lowers z degree. Its top z coefficient in(8) is
\[
                     uD_BJ+(\kappa u+2bx)J=0             \tag{9}
\]
for a constant kappa.
Write J=u^r*J0 with r>=0 and u not dividing J0.
The u0 coefficient of(9), after cancelling u^r, is
\[
                       (r+2)b\,x\,J0(0,x)=0.
\]
Characteristic zero and J0(0,x)!=0 force b=0.
This uses polynomiality at u0 and no upper bound on any degree.

## 5. Complete family obstruction

**KMFI163.** For every A in GL4 and every nonzero homogeneous polynomial F,
the direction A*K3 above admits no normalized input of generic rankDX4.

The preceding steps force A to be lower triangular.
All its diagonal entries a,e,i,n are nonzero.
Rescale V by 1/a and F by a, preserving H=FV and the normalized input.
Its matched first component is u^3. On u1 the projective coordinates are
\[
 x\mapsto(d+ex)/a,\qquad
 y\mapsto(g+hx+iy+ix^2)/a,
\]
\[
 z\mapsto(k+lx+my+nz+mx^2+nxy^2)/a.
\]
This is triangular with constant diagonal.
KMFI129 therefore forces every curvature term to be affine.
But its y component has nonzero quadratic coefficient i/a.
The contradiction closes the entire arbitrary-output-matrix family,
including all repeated or resonant parameter cases and every scalar
degree. No generic-parameter specialization is used.

The exact normalized cubic rank3 controls of KMFI161 remain valid:
A=diag(1,lambda,lambda,lambda), lambda=s+3,F=u^s, s>=0.
They have entire sweep Jacobian gamma^4 but X is independent of z.
Their lambda=n case illustrates why the full-rank condition in KMFI162
and the final obstruction cannot be omitted.

## 6. Evidence and remaining scope

[Helper](cubic_chain.py) · [checker](verify_cubic_chain.py) ·
[model](cubic_chain_model.json) · [receipt](cubic_chain_receipt.json).

The proof is unbounded in the scalar degree and every entry of A.
Whole adjugate and coefficient identities retain the last-row parameters,
the j!=0 coordinate change, repeated-eigenvalue cyclicity identities,
all last-coordinate multiplicity branches and the weighted leading trace.
The rank2 full-rank and cubic rank3 controls retain both sides of the
hypotheses. Formal verification, independent whole-proof replay and
external review are not claimed.

The previous goal turn was progress: KMFI158--161 established the stronger
normal/translation tests and cubic matrix restrictions. This turn closes
that entire arbitrary-matrix cubic family and supplies the broader
last-coordinate multiplicity test.
It does not classify arbitrary cubic directions, other quadratic shapes
of h0 in detDV=delta*h0^4, higher primitive degrees, nonhomogeneous fields
or the full action spectrum. An actual exceptional11 realization or
absolute action obstruction is still missing. The full goal remains active.
