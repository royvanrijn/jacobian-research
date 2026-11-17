# The two remaining quadratic projective normal forms

The remaining rank4 homogeneous quadratic problem with matched V0=u^2
has an exhaustive finite-dimensional reduction. After the already-closed
one-dimensional nonlinear-image case, only two explicit nonlinear
normal forms remain. The arbitrary linear matrix and constant terms
must be retained: triangularizing a normalized quadratic part does not
triangularize the whole projective map.

The reduction imports one precisely attributed linear-triangularization
theorem. It is not an absolute Keller action criterion or a realization.

## 1. Exhaustion of the remaining directions

Work over algebraically closed characteristic-zero constants.
Let V be homogeneous quadratic in p=(u,x,y,z), V0=u^2.
Assume a nonzero homogeneous scalar F gives a normalized input H=FV
with generic full rankDX=4.

**KMFI133.** After a constant homogeneous change preserving u and
the corresponding matching source/target direction change, either
the direction is affine projectively and is classified by KMFI121,
or its projective dehomogenization has exactly one of the forms
\[
              v=M\{(x,y,z)+H(x,y,z)\}+b,\quad\det M\ne0,
\]
\[
              H_A=(0,x^2,y^2),\qquad H_B=(0,x^2,xy).       \tag{1}
\]
M is an arbitrary constant invertible matrix and b an arbitrary constant
vector at this determinant-compatible stage. Every one-dimensional
quadratic-image direction was already excluded by KMFI132.
There is no remaining three-dimensional quadratic-image case.

Conversely each displayed M,b and either H gives a projective quadratic
map with constant nonzero determinant. Its homogeneous direction
V=(u^2,u^2*v(x/u,y/u,z/u)) is primitive quadratic and passes the
Jacobian-power shape gate. This converse does not assert a normalized
polynomial multiplier, full rank, a desired generic group or a completion.

Here is the proof and the exact imported result.
The pure-square first component and KMFI125 force h=lambda*u,
detDv constant nonzero and f=-lambda*u*F.
Write v=Mx+b+Q, where Q is its homogeneous quadratic part.
Since detDv is constant, M=Dv(0) is invertible.
Set H=M^-1 Q. Then
\[
                     \det(I+DH)=1.
\]
Scaling x by an independent t gives det(I+tDH)=1.
Thus DH is nilpotent and its generic rank is at most2.

We use exactly Theorem3.2 of Michiel de Bondt,
[Quadratic polynomial maps with Jacobian rank two](https://arxiv.org/pdf/1601.00579v4):
a nilpotent quadratic Jacobian of generic rank at most2 is similar
over the constant field to a triangular matrix.
The [retained source manifest](sources/two-image-projective-manifest.json)
pins the primary PDF. No inverse-Galois or monodromy assertion is imported.

If T supplies the similarity, the homogeneous quadratic map
T^-1 H(Tx) has triangular Jacobian with zero diagonal.
Integration in characteristic zero and homogeneity give
\[
                      (0,\ a x^2,\ bx^2+cxy+dy^2).       \tag{2}
\]
The nonlinear image has dimension at most2.
If it is at most1, KMFI132 excludes a nonzero quadratic part under full rank.
The zero quadratic case is KMFI121.
Otherwise a!=0 and (c,d)!=(0,0).

For d!=0 put k=c/(2d), B=(b-c^2/(4d))/a and use
\[
 X=x,\qquad Y=(y+kx)/a,\qquad Z=(z-By)/(da^2).
\]
Conjugating(2) by this invertible constant map gives H_A.
For d=0,c!=0 use
\[
 X=x,\qquad Y=y/a,\qquad Z=(z-(b/a)y)/(ac).
\]
It gives H_B. These formulas require no sampling or degree bound.
They give an exhaustive reduction, not uniqueness of the parameters.

Crucially, if the full map was v=M(x+H)+b, its conjugate is
\[
 T^{-1}v(Tx)=(T^{-1}MT)\{x+T^{-1}H(Tx)\}+T^{-1}b.
\]
The new M remains arbitrary. Replacing it by I, or assuming that it
preserves the triangular flag of H, would discard the open cases.

## 2. What the reduction makes finite, and what it leaves open

For either H_A or H_B and arbitrary invertible M,b, the remaining
scalar problem is already degree-free by KMFI125:
h=lambda*u, the Krylov matrix must be invertible, the unique feedback
form must have nonnegative integer logarithmic residues, and the
separate cyclic coframe must be nonzero. If s=degF and Delta=detM,
\[
                         s=\lambda^4/\Delta-2.
\]
This is not a bound on lambda or an assertion that every parameter works.
No multiplier-degree search is needed for a fixed M,b.

The primary theorem triangularizes H=M^-1 Q, not v itself.
Thus the general two forms(1) still include genuinely nontriangular
directions. Their scalar consistency, monodromy and affine completion
are separate gates. In particular this reduction supplies neither
a PSL2(11) input nor an absolute obstruction.

Directions without a matched pure-square component, higher primitive
degrees, nonhomogeneous fields and the full faithful-transitive
Keller action spectrum remain open.

## 3. A coupled companion family meets both normal forms and is excluded

**KMFI134.** For arbitrary constants beta,gamma,tau,sigma and delta!=0,
consider the projective quadratic map
\[
       v=(x^2+y,\ z-\sigma x+\tau y+\beta xy+\gamma y^2,\ \delta x).
                                                               \tag{3}
\]
No nonzero homogeneous polynomial F makes the homogeneous quadratic
direction V=(u^2,u^2*v(x/u,y/u,z/u)) normalized with full rank4.
The conclusion includes every parameter specialization and every
scalar degree. It also applies after matching constant homogeneous
conjugacies; it does not permit an arbitrary unrelated M or b in(1).

Its determinant is detDv=delta, so the exact gate is h=lambda*u.
Let A=DV|_(u=1), w=V|_(u=1), C=(w,Aw,A^2w,A^3w) and
D=detC/delta. If C is singular, KMFI125 already excludes full rank.
Otherwise let eta be the unique characteristic-feedback form for
det(TI-A-w*eta)=(T-lambda)^4.

Write Nx=D*eta_x, Ny=D*eta_y and use the polynomial numerator identity
\[
 D^2(\partial_x\eta_y-\partial_y\eta_x)
   =D(\partial_xNy-\partial_yNx)-Ny\,\partial_xD+Nx\,\partial_yD
   =:\mathcal R.                                         \tag{4}
\]
Its coefficients are defined over the constant parameter ring with
delta inverted. No division by beta or gamma is made.
The exact universal coefficient certificates are
\[
                    [x^4z^5]\mathcal R=-8\beta^2,\qquad
                    [x^2y^6z^3]\mathcal R=32\gamma^4.      \tag{5}
\]
These are identities in all parameters, not sampled values.
A polynomial multiplier requires eta=dlogF, hence a closed form.
Equation(5) forces beta=gamma=0 in characteristic zero.
The explicit polynomial numerator calculation avoids a large rational
simplification, while retaining the entire identity(4).

For beta=gamma=0, normalization is impossible even before full rank.
Let m>=0 be the degree of Fbar=F(1,x,y,z) in x, with nonzero
leading coefficient J(y,z). The trace gate is
\[
  (v-(x,y,z))(Fbar)=
           [4\lambda-2-\tau-\deg F-2x]Fbar.
\]
Comparing x^(m+1) gives
\[
               (-\sigma\partial_y+\delta\partial_z)J=-(m+2)J.
\]
The left operator is a constant directional derivative, locally
nilpotent on polynomials, and cannot have the nonzero eigenvalue -(m+2).
This is a contradiction. Together with the singular-C branch and(5),
it proves the theorem in all parameters and scalar degrees.

This family reaches both remaining nonlinear shapes. Its linear part is
\[
 M=\begin{pmatrix}0&1&0\\-\sigma&\tau&1\\\delta&0&0\end{pmatrix},
\]
and M^-1 times its quadratic part is
\[
                         (0,x^2,-\tau x^2+\beta xy+\gamma y^2).
\]
If gamma!=0 it reduces to H_A; if gamma=0,beta!=0 it reduces to H_B.
It is a genuinely coupled two-image obstruction, not a repetition of
the one-dimensional-image theorem.

## 4. Completely split scalar multipliers cannot leave primitive degree one

**KMFI135.** Let r>=2, p=(u,x), d>=1 and
\[
                    V=(u^d,u^d v(x/u)),
\]
where v is a polynomial automorphism and its degree is at most d.
If a homogeneous F which splits completely into linear forms gives a
normalized full-rank input H=FV, then v is affine.
Consequently V=u^(d-1)*Mp and KMFI121 classifies every survivor.

This is unbounded in rank, direction degree, scalar degree and number
of linear factors. It covers both arbitrary-M,b families(1).
Every full-rank primitive quadratic survivor in those families must
have a nonlinear irreducible multiplier factor; in particular degF>=2.
A quadratic F would have to be irreducible.

The proof uses automorphy, not only constant Jacobian.
The latter gives h=lambda*u^(d-1).
Write
\[
                F=u^e\prod_i l_i(u,x)^{b_i},\qquad b_i>0,
\]
where the l_i are pairwise distinct, not proportional to u, and linear.
On u=1 put q_i=l_i(1,x), Fbar=F(1,x).
The trace identity makes W(Fbar)/Fbar polynomial for W=v-x.
Each separate divisor residue therefore gives q_i divides W(q_i).

Choose an affine coordinate change with first coordinate q_i.
The first coordinate of the conjugated polynomial automorphism v is
W(q_i)+q_i, hence divisible by q_i. A coordinate polynomial of a
polynomial automorphism is irreducible, so
\[
                         q_i(v)=c_i q_i
\]
in affine-coordinate notation, with a nonzero constant c_i.
Precisely the left side denotes dq_i*v+q_i(0), including the affine
constant; it is not the general composition of a nonlinear polynomial.
Thus the displayed equality holds for these affine linear forms.

Let k be the dimension of the span of the constant covectors dq_i.
If k<r-1, let E be the nonzero constant tangent space cut out by
du=0 and all dl_i=0 in the homogeneous source.
The identities for q_i homogenize to
\[
                   l_i(V)=c_i u^{d-1}l_i(p).
\]
They make DV preserve E. Also eta=dlogF annihilates E, so
N=DV+V*eta-hI preserves E and is nilpotent there.
It has a nonzero kernel vector in E over the source function field.
But df=-F(dh+h*eta) annihilates E as well.
The full-rank kernel test KMFI116 is contradicted.
Hence k=r-1. Select that many independent affine linear coordinates q_i.
Their identities determine every component of v as an affine polynomial.

For(1) the automorphism assertion is explicit, without a general
Jacobian-conjecture implication:
\[
 (x,y+x^2,z+y^2)^{-1}=(x,y-x^2,z-(y-x^2)^2),
\]
\[
 (x,y+x^2,z+xy)^{-1}=(x,y-x^2,z-x(y-x^2)).
\]
Multiplication by M and addition of b preserve automorphy.

The conclusion does not exclude multipliers with nonlinear irreducible
factors. Those factors are a necessary remaining feature, not a
construction or a proof of a desired generic action.

## 5. Evidence and scope

[Helper](two_image_projective.py) ·
[checker](verify_two_image_projective.py) ·
[model](two_image_projective_model.json) ·
[receipt](two_image_projective_receipt.json).

The source-imported triangularization is applied only to the homogeneous
quadratic H=M^-1Q. The two normal-form basis changes, arbitrary-M
retention, degree-free scalar gate and coupled companion obstruction
have separate written arguments and exact controls.

The entire arbitrary-M,b families(1) are not excluded or realized here.
Their generic permutation groups are not computed merely by this
reduction. In particular no exceptional11 input, absolute exceptional11
obstruction or complete faithful-transitive classification is supplied.
The next homogeneous quadratic work is confined to the remaining M,b
parameters of these two families, unless it leaves the matched pure-square
hypothesis. Nonhomogeneous and higher primitive directions remain open.

Formal verification, independent whole-proof replay and external review
are not claimed. The goal remains active.
