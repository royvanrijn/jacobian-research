# Quadratic-tilt monodromy, exceptional boundaries and canonical normal forms

14 November 2026. Four uniform statements extend the [receiver theory](../universal-keller-receiver/PROOF.md). The monodromy lemma is a special case of a published theorem of S. D. Cohen, not a new general Galois-group result. Its application removes the characteristic-zero restriction on arbitrary compiler profiles. We also determine their exceptional boundary fibres and a normal form that excludes every marked equivariant construction linear in its moving variables.

## 1. A quadratic-tilt lemma in every odd characteristic

**UAQTILT1.** Let K be any field of characteristic other than two, let h in K[S] have degree N>=3 and h'(0)!=0, and let B,C be independent indeterminates. Then
\[
 h(S)+BS^2-C
 \tag{1}
\]
has arithmetic and geometric Galois group S_N over K(B,C). Here geometric means after extending K to an algebraic closure; B and C remain variables. There is no restriction involving N, N-1, or the characteristic dividing N.

Cohen's 1980 Theorem 1 proves this immediately. In odd characteristic p, the functions 1,S^2,h(S) are linearly independent over K(S^p): the nonzero S coefficient of h gives a nonzero residue-one component, distinct from residues zero and two. His 1981 correction explicitly preserves Theorem 1. We give a direct proof to expose exactly where odd characteristic matters.

Extend K to an algebraic closure. The polynomial (1) is irreducible and separable: it is linear in the independent C, and its S derivative has nonzero constant term h'(0).

For double transitivity, divide the two-root equation by the diagonal:
\[
 A(X,Y)+B(X+Y)=0,\qquad
 A(X,Y)=\frac{h(X)-h(Y)}{X-Y}.                         \tag{2}
\]
The irreducible polynomial X+Y does not divide A, because
\[
 A(X,-X)=\frac{h(X)-h(-X)}{2X}
\]
has nonzero constant term h'(0). Thus (2), viewed as a polynomial linear in B, is primitive and irreducible in K[X,Y,B]. It stays irreducible over K(B). Localizing away from X=Y gives the integral incidence of two distinct roots of the generic fibre; its projection to C is dominant. Equivalently the point stabilizer is transitive on the other roots. The monodromy is doubly transitive.

There is also a transposition. A critical point r is nonzero and satisfies
\[
 B(r)=-\frac{h'(r)}{2r},\qquad
 C(r)=h(r)-\frac{rh'(r)}2.                           \tag{3}
\]
The rational function B(r) has a simple pole at zero, so B'(r) is not identically zero in any odd characteristic. The extension K(r)/K(B(r)) is finite separable. Differentiation with respect to B preserves its intermediate field K(B(r),C(r)); the identity
\[
 \frac{dC}{dB}=r^2                                  \tag{4}
\]
therefore places r^2 in that field. But B(r) is not in K(r^2), since every pole at zero of a function of r^2 has even order. As K(r)/K(r^2) is quadratic, this proves
\[
 K(B(r),C(r))=K(r).                                 \tag{5}
\]
The critical-value curve is birationally parametrized by r.

For generic B, the critical points are simple, because B(r) is separable; their values are distinct, because (5) is birational. At such a point
\[
 h''(r)+2B=-2rB'(r)\ne0.
\]
The finite branch has exactly one double root and all other roots simple. Since two is invertible, its inertia supplies a transposition. The inertia lies in the generic monodromy even if the ramification at infinity is wild. A doubly transitive group containing a transposition is S_N: conjugation produces every transposition. This proves the assertion over algebraically closed K, hence over the original K.

This proof avoids a characteristic-zero argument that finite inertia generates the group. That argument cannot simply be reused in characteristic p, where nontrivial finite etale covers of A1 exist.

The hypotheses matter. Without a linear term, a polynomial such as S^4+B*S^2-C factors through S^2. In characteristic two, S^8+B*S^2+S-C is linearized up to its constant term; its eight roots form an affine three-dimensional F2 space. Its group is contained in AGL_3(F2), of order 1344, strictly smaller than S8.

## 2. One atomicity theorem for every coefficient section

Use the normalized compiler
\[
 t=1+xy,\qquad q=t^2z+y^2(1+3t),\qquad \Pi=tq,
\]
\[
 B_0=y+3xq+\sum_{j=4}^N j a_jt^2x^{j-2}q^j,\qquad
 C_0=x(5-3t)-x^3z-\sum_{j=4}^N(j-2)a_j(xq)^j,
\]
with output (\Pi,-B_0/2,C_0). Substitute a_j=f_j(Pi), where the profiles lie in k(P). UAPROFILE1 already proves the rational Jacobian identity and the exact criterion for polynomial extension on the original affine source.

**UAALLATOM1.** For every field k of characteristic other than two, every such rational profile map with nonzero top profile f_N has generic degree N and arithmetic and geometric monodromy S_N. When the map extends polynomially, it is absolutely and stably atomic under polynomial composition. This includes the exceptional negative profiles in odd characteristic. The seed-free map has the same assertions with degree three.

Indeed over K=k(pi) the complete generic fibre has equation
\[
 E(S)=h(S)+bS^2-c/2,\qquad
 h(S)=S+\pi S^3+\sum_{j=4}^N f_j(\pi)\pi^jS^j.        \tag{6}
\]
Here h'(0)=1 and deg h=N. The source reconstruction in URKC1 identifies the full function-field extension with this polynomial. UAQTILT1 applies, including when p divides N. Applying it after algebraically closing k(pi) gives an S_N subgroup of geometric monodromy after algebraically closing k; therefore that geometric group and the arithmetic group are S_N.

For clarity, the atomicity step is characteristic-free. A polynomial factorization of a Keller map has Keller factors by the chain rule and the fact that polynomial-ring units are constants. A primitive generic permutation action admits no proper intermediate field, so one factor has generic degree one. An etale birational self-map of affine space is an open immersion by Zariski's Main Theorem. Its complement has no divisor, since such a divisor would restrict to a nonconstant unit on the source; normal Hartogs extension and affineness then make the immersion an isomorphism. Thus the degree-one factor is an automorphism. Constant-field extensions preserve the geometric action, and identity stabilizations adjoin purely transcendental variables without changing it.

The same assertion holds with any number of promoted identity parameters u: allow f_j in k(u)(P), require polynomial landing on the actual whole affine source, and keep pi,b,c independent target coordinates. Apply the argument over k(u,pi). The determinant is block triangular with the profile determinant one. This gives a single monodromy-and-atomicity proof for the full, compressed, trace, sextic and double-trace coefficient families, and for further polynomial sections. It does not prove field-versality of an arbitrary section, classify stable equivalence, or improve target-dimension bounds.

## 3. All exceptional pi-zero fibres

Let k have odd characteristic p. The polynomial-landing classification says exactly
\[
 f_j(P)=g_j(P)+u_jP^{-1}+v_jP^{-2},\qquad g_j\in k[P],
 \tag{7}
\]
with u_j=v_j=0 unless j=2 mod p. Coefficients with j outside 4,...,N are absent. Define the separable polynomial
\[
 \phi(Y)=Y-\sum_{j=2\bmod p}v_j(-Y)^{j-2},
 \qquad \phi'(Y)=1,\qquad d=\deg\phi.                 \tag{8}
\]
Thus d=1 if every v_j vanishes; otherwise d is a positive multiple of p, at most N-2.

**UAMODBOUND1.** The complete fibre over (pi,b,c)=(0,b,c) is canonically, with the displayed coordinates,
\[
 \boxed{
 \operatorname{Spec}\bigl(k[S]/(S+bS^2-c/2)\bigr)[1/(1+2bS)]
 \ \sqcup\
 \operatorname{Spec}\bigl(k[Y,Y^{-1}]/(\phi(Y)-b)\bigr).
 }                                                       \tag{9}
\]
This is a scheme-theoretic whole-fibre formula, including all residue-field extensions, not a point count. Its geometric rank is
\[
 \boxed{
 \begin{cases}
 d,& b=0\ \text{or}\ 1+2bc=0,\\
 d+2,& b\ne0\ \text{and}\ 1+2bc\ne0.
 \end{cases}}
 \tag{10}
\]
It extends UABOUND1, whose regular-profile case has d=1.

To prove it, use the existing comaximality of t and q. On q=0, t is a unit and all negative-profile terms in B_0 vanish: their q exponents are j-1 or j-2, both positive. The original quadratic reconstruction applies:
\[
 D=1+2bS,\qquad y=-2b,\quad x=S/D,\quad
 z=-y^2(D^2+3D).
 \tag{11}
\]
This gives precisely the first component in (9), with ranks one if b=0, zero if b!=0 and 1+2bc=0, and two otherwise.

On t=0 one has y invertible, x=-1/y and q=y^2. The regular and simple-pole terms in B_0 vanish except for y+3xq=-2y. The double-pole terms give
\[
 -B_0/2=y-\sum v_j(-y)^{j-2}=\phi(y).
\]
The C_0 equation uniquely determines z:
\[
 z=y^3\left(c+\frac5y+
            \sum_{j=4}^N(j-2)g_j(0)(-y)^j\right).      \tag{12}
\]
Negative-profile C terms vanish because j-2=0 in k. This proves the second component and its inverse reconstruction. It has rank d for b!=0. For b=0, Y=0 is one simple root of phi, deleted by the localization, so the rank is d-1. Adding the quadratic component proves (10).

The exceptional boundary is controlled only by the double-pole coefficients v_j. The simple-pole coefficients do not change its abstract fibre, and the regular coefficients only affect (12).

### The precise limit of boundary escape

For N>=4 and nonzero top profile, a pi-zero fibre has the full generic rank N if and only if
\[
 N=2\bmod p,\quad v_N\ne0,\quad
 b\ne0,\quad 1+2bc\ne0.                              \tag{13}
\]
In that event it is the product of a rank-two etale algebra and a rank-(N-2) etale algebra. Hence **no pi-zero fibre is a degree-N field extension**, even though its rank can now equal N. A nontrivial field fibre at pi=0 can occur only when b!=0, 1+2bc=0 and phi(Y)-b is irreducible; its degree d is a multiple of p and is at most N-2.

Thus the old rank-three boundary bound is false for exceptional profiles, while their inability to realize a full-degree field at that boundary still has an exact proof. This does not impose a dimension bound on another compiler.

### A uniform example in every odd characteristic

For any odd prime p set N=p+2 and f_N(P)=-P^(-2), with all other profiles zero. The polynomial map is simply
\[
 F_p(x,y,z)=
 \left(tq,\ -y/2-3xq/2+(xq)^p,\
              x(5-3t)-x^3z\right).                  \tag{14}
\]
The additional p-th power has zero derivative, so its Jacobian is one. UAALLATOM1 gives degree p+2, full symmetric monodromy and absolute stable atomicity over every field of characteristic p.

Here phi(Y)=Y-Y^p. Over F_p the exact boundary fibres are
\[
 \begin{array}{c|c}
 (\pi,b,c)&\text{complete fibre algebra}\\ \hline
 (0,1,0)& F_p\times F_p\times F_{p^p}\\
 (0,1,-1/2)& F_{p^p}\\
 (0,0,c)& F_p^{\,\times p}.
 \end{array}                                         \tag{15}
\]
The Artin--Schreier polynomial Y^p-Y+1 is irreducible: its roots are alpha+a for a in F_p, and Frobenius translates alpha by -1, giving one orbit of length p. The quadratic component is respectively split rank two, empty, or rank one. Formula (15) supplies all-prime controls of the distinction between boundary rank and full-degree field realization. This uses familiar Frobenius perturbations; no novelty claim for that general mechanism is made.

## 4. Canonical constructions must be nonlinear in the moving variables

This statement is independent of the quadratic-gauge compiler.

**UANORMAL1.** Let G be a finite group, let k be a field in which |G| is invertible, and let F:V->W be a G-equivariant polynomial Keller map between linear representations. Write
\[
 V=U\oplus T,\quad T=V^G,\quad U^G=0.
\]
Using the equivariant isomorphism DF(0) identifies W with V and the corresponding fixed and moving summands. In these coordinates write F=(F_U,F_T), and put f(t)=F_T(0,t). Then
\[
 A(t)=D_UF_U(0,t)
\]
is a polynomial G-equivariant matrix with nonzero constant determinant. It has a polynomial inverse. After a G-equivariant polynomial automorphism of the source one has the global normal form
\[
 \boxed{F(u,t)=\bigl(u+R(u,t),\ f(t)+Q(u,t)\bigr),\qquad
               R,Q\in (u)^2,\quad \det Df\in k^\times.} \tag{16}
\]
Here ideal membership means every monomial has moving-variable degree at least two; coefficient degrees in t are unrestricted.

The source normalization leaves the fixed target marked; it requires no inverse for the fixed-space Keller map f.

The fixed/moving splitting exists by averaging; characteristic zero is sufficient but not necessary. Equivariance gives F_U(0,t)=0, and there is no invariant linear map U->T, so the moving-linear part of F_T is zero. At u=0 the full derivative has diagonal blocks A(t),Df(t). Its constant determinant yields
\[
 \det A(t)\det Df(t)=\det DF\in k^\times.
\]
Both factors are polynomials, hence both are nonzero constants. Thus A(t)^{-1} is polynomial and G-equivariant. Precompose by
\[
 (u,t)\longmapsto(A(t)^{-1}u,t)
\]
to get (16).

In particular, suppose the moving component of the original F is affine linear in u, with arbitrary polynomial dependence on t. It must be exactly A(t)u. Since A(t) is everywhere invertible, the entire inverse image of W^G is V^G, scheme-theoretically. Every complete fibre over a G-fixed target therefore has trivial G action.

Consequently **no such architecture can have a marked fibre equal to the natural N-point S_N set for any N>=2**, regardless of the number of extra trivial directions. Even arbitrary nonlinear invariant terms in F_T cannot repair an affine-linear moving component. Adding identity coordinates on trivial representations preserves this obstruction.

For V equal to the permutation representation plus r trivial lines, this strengthens the earlier fixed-space test UAEQUIV1. That test requires r+1>=2, and r=1 would imply a plane Keller counterexample. Even at r>=2, a candidate must have genuinely nonlinear terms in its nontrivial representation component. After normalizing its first normal derivative, those terms begin in moving degree two. This supplies a specific search gate for a receiver which forgets the generator line.

The normal-form assertion applies whenever |G| is invertible, to linear actions and a marked fixed target, not to every nonlinear action or every stack receiver. For S_N this includes characteristic zero and characteristic p>N; the earlier JC(2) consequence specifically uses characteristic zero. It does not prove that all nonlinear candidates are impossible, and supplies no global lower bound for unmarked field-versal Keller maps.

## 5. What closes, and what remains

UAQTILT1 and UAALLATOM1 replace many family-specific monodromy arguments by one quadratic-tilt rule, extending arbitrary rational profiles and their polynomial maps to every characteristic other than two. UAMODBOUND1 completes the exceptional boundary geometry. UANORMAL1 rules out all marked constructions affine linear in their moving variables, even after arbitrarily many trivial stabilizations.

These results do not improve the characteristic-zero target-dimension bounds, produce a receiver forgetting the generator line, solve JC(2), or prove stable equivalence of different compiler slices. The live next gate for canonical construction is a nonlinear equivariant polynomial realization which satisfies all subgroup-fixed restrictions and the normal form (16). For compression, another affine completion or another compiler is still required.

The [replay](verify.py) checks representative exact identities and whole quotient-algebra reconstructions, with [source pins](sources.json). The general monodromy, all-prime and normal-form results are written proofs. Formal verification, independent whole-proof replay, external review and a literature-priority claim remain absent.
