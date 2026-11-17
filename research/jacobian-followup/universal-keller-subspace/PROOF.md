# Scalar-affine Keller maps over their actual base field

14 November 2026. This extends [UA1CANCLASS1 and UAAFFSUB1](../universal-keller-symmetry/PROOF.md), without assuming that the base contains the full degree-q finite field. The two results are UASUBDEFORM1 (composition-preserving deformation, universal frames and factorization) and UASCALARDIM1 (canonical dimensions and fixed framed dimensions over small finite bases). These are characteristic-p theorems. They do not close the general characteristic-zero descent problem.

## 1. Setup and the two different dimension questions

Let p be prime, m prime to p, and
\[
 d=\operatorname{ord}_m(p),\quad B=p^d,\quad e=ds,\quad q=B^s,
 \qquad s\ge1.                                      \tag{1}
\]
When m=1, use d=1. Suppose k has characteristic p and contains a specified primitive m-th root a; then it contains F_B. Fix the constant permutation group
\[
 G=({\mathbb F}_B^s,+)\rtimes\langle a\rangle,
 \qquad (b,a^j)z=b+a^jz,\quad |G|=qm.               \tag{2}
\]
Its labels are the q vectors. An F_B-linear identification with F_q identifies this permutation group with G_(q,m) from UAAFFSUB1; that identification does not put F_q inside k.

A **canonical affine-bundle realization** means a G-equivariant polynomial Keller map between affine representations, with a fixed k-rational target and its complete marked fibre identified with these labels. It descends to affine bundles over BG, and it may vary after trivializing a twisted bundle. A **fixed framed universal map** is one polynomial map on a fixed affine space, with a specified G torsor over a target open framing the entire rank-q cover, realizing every G torsor over every extension L/k. Frames matter, even when distinct G torsors induce isomorphic unstructured algebras. We do not identify these two requirements.

Atomicity means every polynomial Keller factorization has a polynomial automorphism as a factor; absolute stable atomicity retains this property after extension of constants and adjoining identity coordinates. This is equivalent to excluding two factors of generic degree greater than one. Indeed, in every characteristic a degree-one Keller self-map is etale and birational, hence an open immersion into affine space by Zariski Main. A divisor in its complement would give a nonconstant unit on its affine-space source. If the complement has codimension at least two, normal extension of regular functions and affineness make the immersion an isomorphism. Thus degree one already means a unit. The retained primitive-monodromy argument is characteristic-free with this observation: a nonunit factorization gives a proper intermediate field and a nontrivial block; adjoining independent identity variables preserves the Galois group and its label action. The converse from imprimitivity to polynomial decomposability is not assumed for arbitrary maps.

## 2. A deformation that respects composition

Write an F_B-linear additive polynomial as
\[
 A(X)=\sum_{j\ge0}\alpha_jX^{B^j}.
\]
Since m divides B-1, all the following exponents are nonnegative integers. Define
\[
 \boxed{\mathcal D_m A(u,x)
       =\sum_{j\ge0}\alpha_j u^{(B^j-1)/m}x^{B^j}.} \tag{3}
\]
This operation preserves addition, identity and composition:
\[
 \mathcal D_m(A\circ C)(u,x)
 =\mathcal D_m A\bigl(u,\mathcal D_m C(u,x)\bigr).   \tag{4}
\]
Indeed, each coefficient in the composite has the same u exponent on both sides, because
\[
 {B^i-1\over m}+B^i{B^j-1\over m}
      ={B^{i+j}-1\over m}.                         \tag{5}
\]
Evaluation at u=1 recovers A, so no polynomial data are lost. If A'=1, the whole plane map (u,x)->(u,D_m A(u,x)) has Jacobian one. Equation (4) carries every normalized additive-polynomial factorization to a global polynomial Keller factorization; it does not merely produce a factorization after inverting u.

Let W be a finite F_B-subspace of k. Put
\[
 P_W(X)=\prod_{w\in W}(X-w),\quad
 c_W=P_W'(X)\ne0,\quad A_W=P_W/c_W.                \tag{6}
\]
The standard subspace-polynomial argument shows P_W is F_B-linear additive. For completeness, adjoin one F_B-line at a time: if W=U+F_B b, then
\[
 P_W=P_U^B-P_U(b)^{B-1}P_U.                        \tag{7}
\]
This proves additivity, B-power support, separability and the assertion about the constant derivative by induction. Thus A_W'=1 and ker(A_W)=W. No extraction of p-th roots from coefficients is used, so imperfect k is allowed.

For any subspace U contained in W, set T=A_U(W). Both normalized polynomials in
\[
 \boxed{A_W=A_T\circ A_U}                          \tag{8}
\]
have derivative one, degree |W| and kernel W; hence they are equal. The image T is an F_B-subspace of k, and |T|=|W|/|U|. A full F_B flag in W consequently gives factors of degree B, all defined over k. This is classical subspace-polynomial algebra; (3)--(4) provide its whole-affine-space Keller realization.

## 3. UASUBDEFORM1: arbitrary scalar-affine translation models

Choose nonzero finite F_B-subspaces W_1,...,W_r of k with
\[
 \sum_i\dim_{F_B}W_i=s,
 \qquad W=\prod_i W_i\simeq F_B^s.                 \tag{9}
\]
The product is an abstract direct product in k^r; the W_i need not be disjoint subsets of k. Define
\[
 f_i(u,x)=\mathcal D_m A_{W_i}(u,x),\qquad
 \boxed{F_{\mathbf W,m}(u,x_1,\ldots,x_r)
       =(u,f_1(u,x_1),\ldots,f_r(u,x_r)).}          \tag{10}
\]
This map has Jacobian one on the **whole** A^(r+1). Off u=0 it is finite etale of degree q. On u=0 it is the identity in the x coordinates, so the complete geometric fibre there has rank one. The nonproper target locus is exactly u=0.

**UASUBDEFORM1.** The map (10) is one fixed framed universal Keller map for (2), over every extension L/k, including finite fields, imperfect fields and disconnected torsors. Both arithmetic and geometric generic permutation groups are exactly G. It is absolutely stably atomic if and only if s=1. If s>1, every F_B flag gives an explicit global polynomial factorization; a full flag gives exactly s factors of degree B, each absolutely stably atomic. No chain of factors with generic degrees greater than one can have more than s factors. The theorem classifies these displayed families, not every Keller map with imprimitive G.

### Actual torsor and all-field versality

On u invertible, use the frame scheme
\[
 \lambda^m=1/u,\qquad f_i(u,\mu_i)=v_i\quad(1\le i\le r).
                                                               \tag{11}
\]
Its right action, using affine maps z->az+b, is
\[
 (\lambda,\boldsymbol\mu)\cdot(a,\boldsymbol b)
      =(a\lambda,\boldsymbol\mu+\lambda\boldsymbol b).
                                                               \tag{12}
\]
If u*lambda^m=1, equation (3) gives the exact identity
\[
 f_i(u,\lambda z)=\lambda A_{W_i}(z).              \tag{13}
\]
Thus the full fibre consists of the q distinct vectors mu+lambda*w. The frame scheme is finite etale of rank mq and (12) acts simply transitively on its geometric points. It is a specified G torsor; the label cover is its quotient by the scalar stabilizer. No splitting-field-only substitute is used.

Given a G cocycle (a_gamma,b_gamma) for any L/k, multiplicative Hilbert90 gives nonzero lambda in L_sep with gamma(lambda)=a_gamma*lambda. The vectors lambda*b_gamma form an additive cocycle, so additive Hilbert90 gives mu with gamma(mu)-mu=lambda*b_gamma, coordinate by coordinate. Then
\[
 u=\lambda^{-m}\in L^*,\qquad
 v_i=f_i(u,\mu_i)\in L.                            \tag{14}
\]
The frame has the prescribed cocycle, up to its usual choice of base frame. This proves versality of actual torsors for all L, rather than only specializations of connected G extensions over infinite fields. For m=1 the lambda equation simply fixes lambda=1/u.

The whole frame scheme over the open target can also be written
\[
 \operatorname{Spec}k[\lambda,\lambda^{-1},z_1,\ldots,z_r],
 \quad u=\lambda^{-m},\quad v_i=\lambda A_{W_i}(z_i).
                                                               \tag{15}
\]
It is geometrically integral, so the generic torsor is connected and has full group G over k and over kbar. Faithfulness of the natural label action gives exactly G as permutation monodromy. Equivalently, a nonzero root difference in any coordinate recovers lambda up to a known nonzero member of W_i. This confirms the generic claim without an essential-dimension or density shortcut.

The leading coefficient of f_i is a unit when u is inverted. Its monic equation gives a free basis of powers less than |W_i|; the product basis proves the degree-q assertion. On u=0 each f_i=x_i. A finite etale map over a neighbourhood of that divisor would have locally constant rank, contrary to ranks q and one, proving the asserted nonproper locus.

### The finite normalization and every omitted branch

There is an explicit finite normalization of the frame field over the whole target, not just the punctured chart. Set t=lambda^(-1) and z_i=t*mu_i, and use
\[
 \boxed{\widetilde T=
 \operatorname{Spec}k[t,\boldsymbol z,\boldsymbol v]/
       (A_{W_i}(z_i)-t v_i)_i,\qquad u=t^m.}        \tag{15a}
\]
The first map to A^(r+1)_(t,v) is finite etale of degree q, since each A_(W_i) has constant derivative one. The second map (t,v)->(t^m,v) is finite of degree m and tame along t=0. The algebra is t-torsion-free and becomes the domain k[t,t^(-1),z_i] after inverting t; hence it is a domain. It is regular by the etale description, so (15a) is the finite normalization in the generic frame field.

The action extends polynomially as
\[
 (t,\boldsymbol z)\cdot(a,\boldsymbol b)
      =(a^{-1}t,a^{-1}(\boldsymbol z+\boldsymbol b)). \tag{15b}
\]
Over the generic point of u=0, there are exactly q divisors indexed by w in W, with z=w and residue field k(v_i). Each has ramification index m and inertia a conjugate of C_m; for m=1 the inertia is trivial. This accounts for every affine boundary divisor of the frame normalization, and asserts nothing about target infinity.

Quotient by the scalar stabilizer H=C_m to obtain the finite normalization of the degree-q label field. Its boundary branches correspond to the H orbits on W. There is one branch of index one, indexed by zero, and (q-1)/m branches of index m, indexed by the nonzero scalar orbits. Thus
\[
 \#\{\text{normalization branches over }u=0\}
       =1+(q-1)/m,\qquad
       1+m\,{q-1\over m}=q.                       \tag{15c}
\]
The indices are cycle lengths of inertia on the labels; here every residue field is already k(v_i). The original affine source retains precisely the zero branch. Indeed its coordinates are x_i=z_i/t. Near z=0, write A_(W_i)(z_i)=z_i*U_i(z_i), with U_i(0)=1: then x_i=v_i/U_i(z_i) is regular. At any w!=0 some z_i is a unit and z_i/t has a pole, excluding that branch from the source. The standard normalization factorization of this quasi-finite etale morphism is an open immersion into the finite label normalization, so this pole test identifies all omitted branches.

For m>1 the unique unramified branch survives and all the ramified branches are removed. For m=1 all q normalization branches are unramified, but the affine source still retains only one. In this family the exact geometry therefore distinguishes tame ramification from pure omission; rank one alone does not distinguish them.

### Atomic factors and the precise stopping rule

Blocks containing zero for the affine permutation group (2) are the additive subgroups of W invariant under its scalar complement. Invariance under a is exactly F_B-linearity, since F_p[a]=F_B. Thus every such block has size B^t, and the action is primitive exactly when s=1. This is a statement about the actual natural permutation action, not a claim that all representations of the group have this property.

For s=1, primitive geometric monodromy proves absolute stable atomicity. For s>1, first factor each A_(W_i) along a full flag using (8); then deform using (4). Apply those factors successively in the relevant coordinate, passing u and the other coordinates through. This gives s global degree-B Keller factors. Each factor has geometric group F_B semidirect C_m acting primitively on B labels, and hence is absolutely stably atomic, also with the other coordinates passed through. For m=1, B=p and the same conclusion uses prime-degree regular monodromy.

Conversely, any factor chain yields intermediate fields and therefore nested blocks in the degree-q geometric cover. All nontrivial relative degrees are powers B^t with t>=1; their product is B^s. There can be at most s factors of degree greater than one. This bounds factor-chain length without claiming uniqueness of factors or inferring polynomial factorization from imprimitivity alone. The exhibited flags supply the polynomial factors.

In particular every canonical marked line model from UA1CANCLASS1, whenever it exists over k, upgrades to a fixed universal plane map over that **same** k. The F_q-constant hypothesis of UAAFFSUB1 was sufficient, not necessary. For m>1 the fixed framed minimum is exactly two: dimension one is impossible by UAFINITE1's tame cyclic quotient obstruction, while (10) with r=1 attains two. For m=1 the single finite map A_W itself realizes every translation torsor, and dimension one suffices.

## 4. UASCALARDIM1: finite-base dimension and the extra parameter

For k=F_Q with Q=p^f containing a, put
\[
 r=\lceil e/f\rceil.                              \tag{16}
\]
Since d divides f, partition s into r positive integers at most f/d and choose W_i in F_Q of those F_B dimensions. Thus (10) exists with exactly r+1 coordinates, even when F_q is not contained in F_Q.

**UASCALARDIM1.** For the constant scalar-affine group (2):

1. Its canonical marked affine-bundle dimension is exactly r over F_Q, and exactly one over every infinite k containing a.
2. For m=1 its fixed framed universal Keller dimension is r over F_Q and one over such infinite k, as in UAFINCOUNT1.
3. For m>1 over any infinite k containing a the fixed framed minimum is two. Over F_Q,
\[
 \boxed{\max\{2,\lceil\log_Q(q+m-1)\rceil\}
   \ \le\ \operatorname{ktdim}^{\rm framed}_{F_Q}(G)
   \ \le\ r+1.}                                  \tag{17}
\]
In particular the exact value is r+1 if f divides e, and is two whenever r=1. For r>1 and f not dividing e, the general minimum remains between r and r+1; (17) records the exact point-budget lower bound, without deciding that gap.
4. Among maps that preserve a source/target coordinate, even after polynomial coordinate changes defined over F_Q, the fixed framed minimum is always exactly r+1 when m>1. This extra bound uses the preserved-coordinate hypothesis. It is not a lower bound for all maps in the unresolved cases.

### Canonical dimension

On A^r let G act by x_i->a*x_i+b_i on the source and by v_i->a*v_i on the target. The map
\[
 \Psi_{\mathbf W}(x_1,\ldots,x_r)
       =(A_{W_1}(x_1),\ldots,A_{W_r}(x_r))          \tag{18}
\]
is equivariant, has Jacobian one and complete zero fibre W, so it descends to the required bundles over BG. Over F_Q its split fibre has q rational points, forcing Q^D>=q in ambient dimension D. This proves D>=r and hence exactness. Over infinite k containing F_B there is an s-dimensional F_B subspace in k, so r=1 is possible; a nontrivial marked fibre excludes dimension zero.

Equation (18) chooses subspaces and a fixed identification of the abstract label module once over the base. It does not choose a basis, origin or normal form functorially for each input torsor. Its existence as a G-equivariant map is what gives the canonical bundle construction.

### Exact weighted point budget for fixed frames

The conjugacy classes in G are: the identity; (q-1)/m classes of nonzero translations; and one class for each of the m-1 nonidentity scalars. A nonzero translation has no fixed label. Every nonidentity scalar, including any translation part, has exactly one fixed label, because 1-a^j is invertible on W. Consequently
\[
 \sum_{[g]}|W^g|=q+m-1,\qquad
 c(G)=(q-1)/m+m\le q+m-1.                          \tag{19}
\]
The UAFINCOUNT1 disjoint-whole-fibre argument now gives Q^D>=q+m-1. Distinct framed torsors require distinct targets even if their underlying label algebras coincide. The independent lower bound D>=2 for m>1 follows because a one-variable Keller map is finite and cannot carry a varying prime-to-p cyclic quotient (UAFINITE1). Together with (10), these prove (17).

If e=rf, then q=Q^r and m>1, so q+m-1>Q^r forces D>=r+1. This matches the construction. If r=1, the independent univariate obstruction gives the matching bound two. If e is not divisible by f and r>1, the point budget can fit within Q^r; no missing lower bound is silently supplied by (10).

If F preserves a coordinate, the target of the split framed torsor has its q rational source points in one coordinate hyperplane. That hyperplane has Q^(D-1) rational points, so Q^(D-1)>=q and D>=r+1. Polynomial coordinate changes over F_Q preserve this argument. Formula (10) attains it. This identifies exactly why simply deleting the retained parameter cannot solve the remaining unrestricted dimension gap.

### An all-rank family over the smallest scalar field

For every s>=1, over F_B and m>1 dividing B-1 with ord_m(p)=d, take W_i=F_B. Then
\[
 G=(F_B^s,+)\rtimes C_m,\qquad
 F(u,x_1,\ldots,x_s)
   =(u,x_1-u^{(B-1)/m}x_1^B,\ldots,
         x_s-u^{(B-1)/m}x_s^B).                   \tag{20}
\]
The canonical marked dimension is exactly s, while the fixed framed universal dimension is exactly s+1. For s>1 this fixed map has a full factorization into s absolutely stably atomic degree-B factors. For example over F4, the natural degree-16 group (C2)^4 semidirect C3 has canonical dimension two and fixed framed minimum three, attained by (u,x-u*x^4,y-u*y^4). F16 is not contained in F4.

## 5. Scope and remaining gates

- All conclusions retain characteristic p, the constant group (2), its natural labels and the specified scalar root in k. A missing root of unity, nonconstant group scheme, nonscalar complement or characteristic-zero source is outside this theorem.
- Canonical minima above have no atomicity requirement. In a line the marked map necessarily has regular translation monodromy and decomposes when e>1. The deformed universal family has the different stopping rule s=1. A fieldwise normal form is not a new canonical gluing theorem.
- For m>1 every displayed universal map has rank-one boundary u=0. For s>1 the group is imprimitive, and the actual flag factors establish non-atomicity of this family. Other absolutely atomic framed universal maps for such groups are not excluded.
- Atomic canonical lifts over weaker bases are not assigned the old two-coordinate bound that assumed F_q constants. The general matrix construction UALANG1 remains available over F_B, but a sharper bound would require a separate argument.
- The unresolved finite-base interval in (17) asks for a map that can use every source coordinate for a split fibre while still varying the tame quotient. A coordinate-preserving ansatz cannot attain its smaller endpoint.
- General characteristic-zero canonical S_N receivers, their nonlinear fixed geometry, and the characteristic-two quartic three-versus-four fixed-map gap remain open. No new result here decides them.
