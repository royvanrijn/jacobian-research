# Branch-centered completion and the binary degree-eight action

Work geometrically over algebraically closed characteristic-zero constants.
A normalized polynomial tautological sweep here means
\[
 S(\gamma,w,z)=X(w,z)+\gamma(1,w,z),\qquad \det DS=\gamma^2,
\]
with nonzero polynomial nilpotent normalization
\[
 P=\begin{pmatrix}-\mu&\lambda\\-\nu&\mu\end{pmatrix},
 \quad\lambda\nu=\mu^2,
\]
and the reconstruction identities of KMFI41. Assume rank DX=2.
Let L=k(gamma,w,z),K=k(S), and let (G,H) be its geometric
normal-closure group and sheet stabilizer.

**KMFI59 (complete completion criterion for this class).**
The following are equivalent:

1. H normally generates G.
2. For the irreducible branch image D=closure(X(A2)), with reduced
   equation d(Y)=0, write d(S)=gamma^3*h. Then h is nonconstant.
3. Some value X(p*) with nonzero nilpotent normalization has a
   preimage (gamma0,p0) under S with gamma0!=0.
4. The field pair K subset L, after adjoining one common independent
   variable C, has a polynomial Keller model on ordinary affine4.

The criterion applies to nonfinite as well as finite sweeps.
The model is constructed by affine branch centering and one
quadratic triangular correction. Over an original nonclosed
coefficient field, its constants may require a finite algebraic
extension. No arithmetic descent to the original field is asserted.

For a finite sweep of degree N and branch degree
delta=[k(w,z):k(X)], the criterion is exactly N>3delta.
By KMFI32, delta<=2, so **every finite sweep in this class with
degree greater than6 admits this completion**.
This is a necessary-and-sufficient result for the specified
normalized rank-two field pairs, not for arbitrary Keller actions.

**KMFI60 (actual binary degree-eight occurrence).**
The faithful action of SL2(F3) on its eight nonzero vectors has an
explicit polynomial Keller realization on ordinary affine4 over
Q(sqrt(-3)). Both geometric and arithmetic generic groups over
that coefficient field are SL2(F3); its normal closure has no
additional constant extension. The complete formulas are below.

**KMFI61 (atomic and not independently curve-induced).**
This map is absolutely and stably atomic under polynomial Keller
composition although it has a unique pair block. Its generic deck
group is C2 and has no nonidentity regular extension to the source.
Its cycle-gcd lcm is4<8, so it lies outside every independently
stabilized curve-cover Keller presentation. The A4 block quotient
does not realize natural A4 as a Keller factor.

[Formula helper](branch_centered.py) ·
[checker](verify_branch_centered.py) ·
[expanded polynomial coordinates](binary_keller_coefficients.json) ·
[model](branch_centered_model.json) · [receipt](branch_centered_receipt.json).
The earlier centered-origin failure is retained in
[KMFI45--46](POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md); branch translation
closes that presentation-specific gap.

## 1. Why a normally generating stabilizer gives a regular tail

Since rank DX=2, D is an irreducible hypersurface of the affine
target. The generic point of gamma0 maps generically finitely to D.
KMFI39's local different calculation gives ramification index3:
the Jacobian has order2, and the residue extension is separable.
Consequently d(S) has **exactly** gamma-adic order3, so
\[
                      d(S)=\gamma^3 h,\quad \gamma\nmid h.   \tag{1}
\]
Away from gamma0 the map S is etale. Every other prime factor
of d(S) therefore has multiplicity1 and maps dominantly to D.
Thus h, up to a constant, is squarefree.

If h is constant, gamma^3 belongs to K. But gamma does not belong
to K: its order along gamma0 is1, whereas every element of K has
order divisible by3 at that valuation. Hence K(gamma)/K is a
nontrivial cyclic Galois extension of degree3, since the geometric
constants contain the cube roots of unity. It is an intermediate
subfield of L. Equivalently H lies in a proper normal subgroup
of G and cannot normally generate G. This proves1 implies2.

If h is nonconstant, take a generic point of a prime factor of h.
It has gamma!=0 and is etale over D. Its image contains a dense
open of D. The image under X of points where rank DX=2 and P!=0
also contains a dense open of D. Their intersection supplies p*
and a regular second preimage with gamma0!=0.
This proves2 implies3. All choices can be made over a finite
algebraic extension of the original coefficient field.

For finite S, the full degree above the generic point of D is
3delta+u, where u is the degree of the unramified tail.
Thus h nonconstant is equivalent to N>3delta.
There is no omitted-normalization assumption in the nonfinite
argument: it uses an actual etale prime of the polynomial source.

## 2. Universal local centering and completion

Choose an affine source change p_old=T p+p* and the matching
target change
\[
 (Y1,Yvec)\mapsto
 \big(Y1-X1(p*),\
 T^{-1}(Yvec-p*Y1-[Xvec(p*)-p*X1(p*)])\big).
\]
The new polynomial X has X(0)=0 and the same tautological direction
(1,w,z). Choose the first column of T along the image of the
nonzero nilpotent P(p*). The normalization at0 is then
\[
                 \lambda(0)\ne0,\quad \mu(0)=\nu(0)=0.
\]
Since nu=mu^2/lambda near0, also dnu(0)=0.

The reconstruction derivatives give
\[
 X2_w=wX1_w-\mu,\quad
 X3_w=zX1_w-\nu,\quad X3_z=zX1_z+\mu.
\]
At0 this implies
\[
 X2_w=0,\quad X3_w=X3_z=X3_{ww}=0.                           \tag{2}
\]
Give w,z weights1,2. Equation(2), together with X(0)=0, is the
complete forbidden-monomial list for weighted orders
X1>=1,X2>=2,X3>=3. Therefore
\[
 E_i(\gamma,u,v)=
 S_i(\gamma,\gamma u,\gamma^2v)/\gamma^i,\quad i=1,2,3         \tag{3}
\]
are polynomials, without any global homogeneity requirement.

The regular second preimage of0 gives theta0=(gamma0,u0,v0)
with E(theta0)=0,gamma0!=0. At this zero the determinant chain is
\[
                 \det DE(\theta0)=\gamma0^{-1}.
\]
Let L0=DE(theta0)^-1, let c1 be its first column, and put
f20=-(D_c1^2 E3)(theta0)/2. The polynomial isomorphism
\[
 \theta=\theta0+
 L0(v1,v2,v3+f20\,v1^2)^t                                  \tag{4}
\]
has constant Jacobian det L0=gamma0.
For source coordinates (x,y,z_src,t_src) set
\[
 v1=xy,\ v2=x^2z_{\rm src},\ v3=x^3t_{\rm src},\
 C=x\gamma,\qquad F=(C,E1/x,E2/x^2,E3/x^3).                  \tag{5}
\]
The weights below1,2,3 in E after(4) vanish: constants vanish at
theta0, linear terms are (v1,v2,v3), and f20 cancels the sole
remaining term v1^2 in E3. Thus all four entries in(5) are
polynomials on ordinary affine4. The full determinant factors are
\[
 x^6\det L0\cdot\gamma\cdot\gamma^2\gamma^3/C^6
                            =\det L0=\gamma0\ne0.          \tag{6}
\]
This identity holds globally by polynomiality.

The target equations S_i=C^i F_(i+1), together with
x=C/gamma,u=w/gamma,v=z/gamma^2 and the inverse of(4), recover
the generic source and target fields in both directions.
Thus(5) has exactly K(C) subset L(C). This proves3 implies4.
Finally4 implies1 by the general Keller normal-generation theorem
KMFI1, which is an action invariant and survives independent
stabilization. All four conditions are equivalent.

No necessity claim here relies on keeping a chosen branch divisor
affine under an arbitrary rechart: the necessity is KMFI1 and the
cyclic intermediate field in(1), both invariant field statements.

## 3. Explicit binary input and its regular tail

Take the finite sweep
\[
 \sigma=w^3-z,\quad
 X=(\sigma,w^4/2,w^6/3-\sigma^2/2),\qquad S=X+\gamma(1,w,z).
                                                               \tag{7}
\]
Its normalization is lambda=w,rho=w^2.
KMFI45 proves its exact geometric action SL2(F3) on8, and over
Q its arithmetic normal closure has group GL2(F3) with constants
Q(sqrt(-3)). These full-group results are retained and reused.

Write A=Y2,B=Y3+Y1^2/2. Its quartic block equation is
\[
 Q(q)=(A-q^2/2)^2+2q^2(A-q^2/2)+2q^4/3-2qB.
\]
The branch image has equation d(Y)=9B^2-8A^3, and direct pullback
gives
\[
                   d(S)=\gamma^3(9\gamma+4w^3)/4.          \tag{8}
\]
The unramified tail is present.

Choose the branch point p*=(1,0). Then X(p*)=(1,1/2,-1/6).
At this target A=1/2,B=1/3 and
\[
                       Q(q)=-(q-1)^3(q+3)/12.
\]
Let r^2=-3. A regular preimage is
\[
             w=r,\quad z=-1-5r/3,\quad \gamma0=4r/3.
\]
Substitution in(7) gives X(p*) exactly, and gamma0!=0.

Use old w=1+w,old z=w+z, with target rows
\[
 Y1'=Y1-1,\quad Y2'=Y2-Y1+1/2,\quad
 Y3'=Y3-Y2+Y1-1/3.
\]
These are affine isomorphisms. They give
\[
\begin{split}
 X1'={}&w^3+3w^2+2w-z,\\
 X2'={}&w^4/2+w^3+z,\\
 X3'={}&-w^6/6-w^5-2w^4+w^3z-4w^3/3
                                      +3w^2z+2wz-z^2/2.
\end{split}                                                   \tag{9}
\]
The new normalization is lambda=1+w,rho=w^2+2w, so rho(0)=0.
Its regular marked coordinates in(3) are
\[
            \theta0=(4r/3,\ (r+3)/4,\ r/2).
\]

## 4. A fully specified ordinary-affine binary map

For this example the polynomials(3) are
\[
\begin{split}
 E1={}&1+2u+\gamma(3u^2-v)+\gamma^2u^3,\\
 E2={}&u+v+\gamma u^3+\gamma^2u^4/2,\\
 E3={}&v+2uv-4u^3/3+
       \gamma(-2u^4+3u^2v-v^2/2)\\
     &\hspace{13mm}+\gamma^2(-u^5+u^3v)-\gamma^3u^6/6.
\end{split}
\]
Use(4) with the exact constants
\[
 L0=\begin{pmatrix}
 -8&-32r/3&56/3\\
 3/2-r&2r+4&5r/2-7/2\\
 r/4+15/4&5r&-10
 \end{pmatrix},\qquad f20=10r-18.                            \tag{10}
\]
The field is Q(r),r^2=-3. The checker verifies
E(theta0)=0, DE(theta0)*L0=I,det L0=4r/3 and the second derivative
defining f20. Equations(5),(9),(10) are the requested polynomial map;
they contain no unsolved choices or coefficient equations.

The expanded coordinates have respectively5,88,137,381 monomials,
and total degrees5,19,22,33. Exact sparse division verifies every
quotient in(5), including x0. An independent computation reconstructs
all coordinates from their serialized rational coefficient pairs and
computes the **entire expanded polynomial Jacobian**, obtaining
exactly4r/3. This is a second global verification of(6); pointwise
Jacobian controls on and off x0 are only supplemental.

The affine transformations and the two-sided field reconstruction
prove degree8 and the exact eight-letter action. Since the coefficient
field already contains sqrt(-3), KMFI45's arithmetic normalizer
calculation gives arithmetic group SL2(F3) as well, with no further
constant extension. This proves KMFI60.

## 5. Blocks, atomicity, generic deck and the curve obstruction

In the nonzero-vector action let H fix a vector. Then |H|=3.
There is exactly one proper intermediate subgroup
\[
             H<K_0<G,\quad |K_0|=6,\quad [K_0:H]=2.
\]
Indeed SL2(F3) has a unique involution, the central -I.
Any order6 subgroup containing H must therefore be the cyclic
group H times<-I>. There is no order12 subgroup because the
abelianization is C3 and hence has no C2 quotient.
The subgroup interval and normalizer are also enumerated exactly
in the checker: N_G(H)=K_0.
Thus the unique proper block is the partition into four pairs,
and N_G(H)/H=C2 is the generic deck group.

A nontrivial polynomial Keller factorization would have a proper
intermediate field. Its only possibility would make the first
factor degree2, a Galois Keller extension, excluded by KMFI1.
Birational Keller factors are automorphisms as used in KMFI30.
The argument survives algebraically closed constant extensions,
identity stabilizations and polynomial left-right equivalence.
This proves absolute and stable atomicity.

In the original coordinates(7), the generic involution is
\[
 (\gamma,w,z)\mapsto(-\gamma,-w,z-2w^3-2\gamma).
\]
It fixes S. With C fixed it sends the actual source coordinate
x=C/gamma to -x. The polynomial gamma in the completed source
restricts to the nonzero constant gamma0 on x0.
A regular extension of the involution would preserve x0 and pull
that restriction back to gamma0, contrary to gamma mapping to
-gamma. Thus the nonidentity generic deck map is not regular,
also in agreement with KMFI47.

Finally the complete permutation census gives cycle-gcd lcm
C(G,Omega)=4<8. KMFI9--10 therefore exclude every independently
stabilized curve-cover Keller presentation of this action, including
birational recharts. The new actual map is another affirmative
instance beyond the curve-generated spectrum.
The intermediate four-sheet A4 quotient is not an affine-space
Keller factor and does not decide absolute natural A4 occurrence.
This proves KMFI61.

## 6. Scope of the progress

All field actions in the completed k(w)[z] chart now have actual
ordinary-affine-space action witnesses except S2 and A3, which
KMFI1 excludes. Its binary exception is no longer an unresolved
absolute-occurrence case. The chart still excludes both PSL2(11)
actions and natural A4 as field inputs.

KMFI59 proves an exact sufficient-and-necessary completion criterion
for normalized rank-two tautological polynomial sweep field pairs.
It does not prove that every permissible action has such a field
presentation. Absolute PSL2(11), natural A4, other construction
classes and the full faithful-transitive action spectrum remain open.

The proofs are written with exact polynomial, finite-group and
expanded-map controls. No formal verification, independent
whole-proof review or claim of novelty priority is asserted.
