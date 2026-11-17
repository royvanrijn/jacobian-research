# Quadratic trace obstructions and universal Keller maps in dimension N-3

14 November 2026. This note gives one exact arithmetic criterion in every rank N>=7, rather than a sequence of rank-specific coefficient searches. It constructs a determinant-one map D_N on A^(N-3), with symmetric monodromy, and characterizes exactly which finite etale algebras it realizes. A binary-block trace test makes the same map field-versal under explicit base-field conditions.

In particular, over every characteristic-zero field containing sqrt(-14), the rank-seven map is field-versal on A4. Duncan's essential-dimension theorem then gives the first exact minimum obtained here:

\[
 \boxed{\operatorname{ktdim}_k(S_7)=4
              \quad(\operatorname{char}k=0,\ -14\in k^{\times2}).}
 \tag{1}
\]

A bound on quadratic forms over k alone is not enough to assert field-versality over all extensions of k. The note keeps those quantifiers separate. It does not assert an N-3 bound over every characteristic-zero field.

## 1. The retained compiler and the new sections

Use the normalized polynomial map F_a=(Pi,-B/2,C) from [UACOMP1](../universal-keller-compression/PROOF.md), over Z[1/2]:
\[
 t=1+xy,\qquad q=t^2z+y^2(1+3t),\qquad \Pi=tq,
\]
\[
 B=y+3xq+\sum_{j=4}^N j a_jt^2x^{j-2}q^j,\qquad
 C=x(5-3t)-x^3z-\sum_{j=4}^N(j-2)a_j(xq)^j.
 \tag{2}
\]
Its determinant is one. For pi nonzero the complete fibre algebra is
\[
 (k[S]/E)[1/E'],\quad
 E=S+bS^2+\pi S^3+\sum_{j=4}^N a_j\pi^jS^j-c/2.
 \tag{3}
\]
Covariance under a'_j=alpha^(j+1)a_j is
\[
 F_{a'}(\alpha x,\alpha^{-1}y,\alpha^{-2}z)
       =(\alpha^{-2}\Pi,\alpha^{-1}b,\alpha C).
 \tag{4}
\]

Define D_N by promoting its displayed seed parameters as identity coordinates.

For N>=8 set
\[
 a_4=a_5=w,\qquad a_{N-1}=a_{N-2}=0,
\]
with a6,...,aN other than those two zero coefficients free. This gives N-6 parameters and three moving coordinates.

For N=7 use a different pair:
\[
 \boxed{a_4=h^3,\quad a_5=a_6=0,\quad a_7=h^5,\qquad
 D_7(h,x,y,z)=(h,F_a(x,y,z)).}                         \tag{5}
\]
The surviving weights five and eight are coprime. If A=a4 and B=a7 are nonzero, the exact root-free normalization is
\[
 \alpha=B^3/A^5,\qquad h=B^5/A^8,\qquad
 \alpha^5A=h^3,\quad\alpha^8B=h^5.                    \tag{6}
\]
Its Laurent inverse is A=alpha^(-5)h^3, B=alpha^(-8)h^5. Thus (5) is a polynomial extension of the quotient section, exactly as in the sextic weights-five-and-seven construction.

Every D_N is a fixed polynomial self-map of A^(N-3) over Z[1/2] with determinant one. These formulas define maps over every odd-characteristic field, regardless of whether that field satisfies the realization criterion below.

### The general weight rule

For positive weights w1,...,wr on a coefficient torus, a root-free normalization of every generic coefficient tuple by its scaling parameter exists exactly when gcd(w1,...,wr)=1. With that hypothesis there is an integral unimodular matrix M=[w | B]. Adding a sufficiently large nonnegative multiple of the first column w to each other column makes every entry of B nonnegative without changing the determinant. Hence
\[
 a_i=\tau^{w_i}\prod_{\ell=1}^{r-1}h_\ell^{B_{i\ell}}
\]
is a torus isomorphism, and tau=1 extends to a polynomial monomial section on A^(r-1). This needs a primitive weight vector, not a coprime pair.

For example weights (6,10,15) have no coprime pair but admit
\[
 (a_5,a_9,a_{14})=(h_1,h_2,h_1h_2),
 \qquad \tau=a_5a_9/a_{14}.
\]
The matrix with columns (6,10,15), (1,0,1), (0,1,1) has determinant minus one.

For g=gcd(w_i)>1 the action has generic group-scheme stabilizer mu_g. After an integral torus-coordinate change its moving coordinate transforms by weight g. Normalizing an arbitrary such coordinate requires a g-th root, excluded generically by its order-one valuation. This concerns transport of every coefficient tuple to a chosen section, not existence of a coarse quotient or a chosen representative of each geometric orbit. Target-dependent scaling is a different operation: it does not automatically give affine source and target spaces or a polynomial Keller self-map in the smaller dimension.

## 2. The exact trace-quadric criterion

Let A/k be a finite etale algebra of rank N, with k infinite and characteristic other than two. Set
\[
 V_A=\ker\operatorname{Tr}_{A/k},\qquad
 Q_A(\theta)=\operatorname{Tr}_{A/k}(\theta^2)\quad(\theta\in V_A).
 \tag{7}
\]
The perfect trace pairing makes V_A have dimension N-1.

Its restricted radical is
\[
 \operatorname{rad}(Q_A)=V_A\cap k1
 =\begin{cases}0,&N\ne0\text{ in }k,\\ k1,&N=0\text{ in }k.\end{cases}
 \tag{8}
\]
Indeed V_A=1^perp and V_A^perp=k1. The reduced nondegenerate trace form has dimension N-1 or N-2 respectively. In particular these ranks are at least five for N>=7.

**UADOUBLE1 (exact realization theorem).** For every N>=7 and every infinite field k of characteristic other than two, the following are equivalent:

1. A is a complete fibre of D_N at some k-rational target, anywhere in the affine target.
2. Q_A has a zero outside its radical.
3. A has a primitive generator theta satisfying Tr(theta)=Tr(theta^2)=0 and the additional good coefficient conditions below.

For N=7 those good conditions are d1*d3*d4 nonzero. For N>=8 they are d1*d3*d4*d5 nonzero, where P_theta(T)=sum d_j T^j is the monic characteristic polynomial. Moreover D_N has degree N, geometric and arithmetic monodromy S_N and absolute stable atomicity over every field of characteristic other than two.

### Why one smooth isotropic point supplies a good generator

Over an algebraic closure, the equations in (7) are
\[
 \sum r_i=0,\qquad \sum r_i^2=0.                       \tag{9}
\]
The affine cone X_(1,2) they define is geometrically integral: after eliminating the first equation, its quadratic form has rank N-1 or N-2, hence at least five. Its singular locus is its radical, consisting only of zero or the scalar line. Every primitive point is outside this locus.

The elementary symmetric functions e3,...,eN are algebraically independent on its coefficient quotient. A direct reason is that every monic polynomial with e1=e2=0 has roots over the algebraic closure, so the root-to-coefficient map is surjective on geometric points. In particular each good coefficient function is nonzero on X_(1,2). They correspond to e_(N-j), with indices at least three in the stated ranks. The discriminant is nonzero there as well: T^N+T+c is separable for a generic c in every characteristic, since its derivative has constant term one.

Their product therefore cuts out a nonempty open in X_(1,2), defined over k. If Q_A has a smooth k-point, its projective quadric is k-rational by projection from that point, including when the quadratic form has a radical. Its smooth affine cone is rational as well. Over an infinite field, the good open consequently has a k-point. This proves (2)=>(3). It is the all-degree upgrade from a merely isotropic trace vector to a primitive generator satisfying every compiler denominator condition.

### The realizing target and its converse

Newton's identity, with two invertible, gives
\[
 d_{N-1}=d_{N-2}=0
       \quad\Longleftrightarrow\quad
 \operatorname{Tr}\theta=\operatorname{Tr}(\theta^2)=0. \tag{10}
\]

For N>=8 apply the root-free formulas in UATRACE1 to a good theta, retaining both coefficient zeros. They give a target of D_N and the entire inverse polynomial
\[
 E(S)=\frac{\alpha}{d_1}P_\theta(S/\alpha),\qquad
 \alpha=\frac{d_3d_4}{d_1d_5}.
 \tag{11}
\]
The squarefree degree-N algebra in (3) is A.

For N=7 a monic presentation has the form
\[
 P_\theta(T)=T^7+d_4T^4+d_3T^3+d_2T^2+d_1T+d_0.
\]
Before scaling, a4=d4*d1^3/d3^4 and a7=d1^6/d3^7. The normalization and target are
\[
 \alpha=\frac{d_1^3}{d_3d_4^5},\quad
 h=\frac{d_1^6}{d_3^3d_4^8},\quad
 \pi=\frac{d_3^3d_4^{10}}{d_1^7},\quad
 b=\frac{d_2d_3d_4^5}{d_1^4},\quad
 c=-\frac{2d_0d_1^2}{d_3d_4^5}.                       \tag{12}
\]
Again E=(alpha/d1)*P_theta(S/alpha). These prove (3)=>(1) without extracting roots.

Conversely, a rank-N complete fibre cannot have pi=0, because the complete pi-zero fibre has rank at most three by UABOUND1. A degree drop or derivative localization at a repeated root cannot give rank N. Thus a realizing target has a squarefree degree-N polynomial E in (3). Its monic form has its last two nonleading coefficients zero, so its primitive root satisfies (10). A primitive element cannot belong to the radical k1. This proves (1)=>(2), including targets outside the chosen good chart.

### Generic monodromy and atomicity

The geometrically integral ordered-root variety X_(1,2) has generically free S_N action. Its root function field is the splitting field of its generic polynomial. The coefficient field has splitting degree at most N!, while the S_N permutation action supplies N! automorphisms. Thus the coefficient field is the invariant field and the generic monodromy is S_N in every odd characteristic.

Its N-2 nonleading coefficients other than the first two are algebraically independent. Passing to inverse polynomials with linear coefficient one is a birational change of these coefficient coordinates on d1*d3 nonzero. The sections (6) and a4=a5 identify the uncompressed cover with the D_N cover after adjoining one independent scaling variable. Purely transcendental extension preserves the group and root action, proving the degree and arithmetic and geometric groups claimed.

The characteristic-free primitive-monodromy argument retained in UATRACE1 gives absolute stable atomicity. The degree-one-factor step uses Zariski's Main Theorem, units and normal Hartogs extension, not a characteristic-zero ramification assertion.

## 3. Arithmetic rules, with the quantifiers exposed

Let U(k) denote the supremum of dimensions of all anisotropic nondegenerate quadratic forms over k; this is not the torsion-only convention sometimes used over real fields.

If
\[
 U(k)<\begin{cases}N-1,&N\ne0\text{ in }k,\\
                   N-2,&N=0\text{ in }k,\end{cases}   \tag{13}
\]
then D_N realizes every rank-N finite etale algebra over k. In particular U(k)<=4 suffices for every N>=7. This includes nonarchimedean local fields of characteristic other than two and totally imaginary number fields, by the classical local quadratic-form bound and Hasse--Minkowski.

These are fixed-field realization statements. The bound U(k)<=4 need not hold over extension fields, so it does not by itself establish the programme's field-versal invariant ktdim_k(S_N). For that invariant the statement must hold for every extension L/k.

For a number field k the criterion is especially sharp:
\[
 \boxed{A\text{ is a fibre of }D_N
 \Longleftrightarrow
 A\otimes_{k,v}\mathbb R\ne\mathbb R^N
 \text{ at every real place }v.}                    \tag{14}
\]
Indeed if A_v=R^r x C^s, the trace form on A_v has signature (r+s,s). Its trace-zero hyperplane has signature (r+s-1,s), since the unit has positive norm N and its orthogonal complement is V_A. It is isotropic exactly when s>0 for these ranks. At finite places its dimension N-1>=6 forces isotropy; at complex places this is automatic. Hasse--Minkowski proves (14).

Over any formally real field, the split algebra k^N is excluded because a nonzero sum of squares cannot be zero. Over R, D_N realizes exactly the rank-N finite etale algebras with at least one complex factor. This recovers the old real-split obstruction and proves its sufficiency counterpart.

### A finite binary-block test that does imply field-versality

Write N=n1+...+nr in binary, with distinct powers of two n_i=2^(m_i). Suppose k has characteristic other than two and contains constants c_i such that
\[
 \sum_i n_i c_i=0,\qquad \sum_i n_i c_i^2=0,
 \qquad(c_1,\ldots,c_r)\text{ is not constant}.        \tag{15}
\]

**UABINARY1.** Under (15), D_N realizes every rank-N finite etale algebra over every infinite extension L/k, for N>=7.

Here is the quantifier-preserving argument. An S_N torsor has, after an odd-degree field extension L'/L, a reduction to a fixed Sylow two-subgroup of S_N. Its binary root blocks have sizes n_i; hence A tensor L' has a product decomposition into finite etale algebras of those ranks. The element with scalar c_i on block i satisfies both trace equations and is outside the scalar line.

More explicitly, the reduction space has odd total degree [S_N:P_2]. It is finite etale, so at least one residue field has odd degree. This avoids assuming that every Galois group already has exactly the binary block orbits.

If N is invertible in L, Q_A is nondegenerate. If the characteristic divides N, use its nondegenerate quotient V_A/(L1). In both cases (15) supplies a nonzero isotropic vector over L'. Springer's odd-degree theorem supplies one over L. UADOUBLE1 then supplies a good primitive generator and a complete fibre. The constants c_i lie in k, so this works uniformly over all extensions L/k. This is the classical Sylow/Springer method applied to the Keller compiler, not a new quadratic-form theorem.

When N is invertible in k, condition (15) is also necessary. This imports Brassil--Reichstein, Theorem1.3 and Corollary8.1(c)--(d), at the prime two: their equivalence tests the general rank-N extension, and the generator conclusion applies when N is invertible. Combining it with UADOUBLE1 gives an exact field-versality criterion for this particular polynomial map.

In characteristic zero the criterion has especially simple consequences. If the binary expansion of N has at most two ones, D_N is never field-versal, even over C: the binary trace hyperplane is zero-dimensional or carries a nonzero one-variable quadratic form. If the base field is quadratically closed and the expansion has at least three ones, D_N is field-versal.

For three binary blocks n1,n2,n3, the restricted binary form has determinant n1*n2*N/n3. Thus in characteristic zero its isotropy is equivalent to
\[
 -Nn_1n_2n_3\in k^{\times2}.                         \tag{15a}
\]
In particular D7 is field-versal exactly when -14 is a square. This is a criterion for D_N, not a lower bound excluding every different N-3-dimensional map.

For reference, the hypotheses of Edens--Reichstein's Lemma4.2 imply (15): in odd characteristic p dividing N, binary digit count r>=5 suffices; r=4 suffices if k contains F_(p^2). Their constant-vector construction can take one c_i=0 and is nonconstant. Thus those base fields give field-versal Keller maps in dimension N-3. Their separate N-4 essential-dimension bound does not by itself give an N-4-dimensional affine Keller map: transporting the surrounding polynomial map is another requirement.

### Exact rank-seven minimum

For N=7 the binary blocks are 4,2,1. If delta^2=-14, take
\[
 (c_1,c_2,c_3)=(3,-4+\delta,-4-2\delta).
\]
Directly
\[
 4c_1+2c_2+c_3=0,\qquad
 4c_1^2+2c_2^2+c_3^2=0.                              \tag{16}
\]
In characteristic zero this vector is nonconstant. Hence D7 on A4 is field-versal over every characteristic-zero k containing delta, including Q(sqrt(-14)) and C. Its ordinary coordinate degrees are (1,7,49,47).

**UAK7MIN1.** Duncan proved ed_k(S7)=4 over every characteristic-zero field. The retained compiler-independent bound ed_k(S7)<=ktdim_k(S7), together with the field-versal map (5), proves (1). This computes a global minimum for the programme's defined target-dimension invariant on these base fields. It does not compute a coordinate-degree minimum, assert the same minimum over Q or R, or remove the extension-field quantifier.

## 4. Canonical quadratic lift spaces and the remaining defect

Let T_N be the squarefree finite target chart of D_N with pi and the chosen section coefficients nonzero. For N>=8 also require aN nonzero; for N7 use h nonzero. Let B_N=O_(T_N)[S]/E.

For a fixed A define H_A as the good-generator open in the projective trace quadric P(Q_A=0) inside P(V_A). Its dimension is N-3, and every geometric point is smooth. The normalizer alpha in (11) or (12) has degree minus one under scaling, so beta=alpha(theta)*theta depends only on the generator line. The exact same inverse construction as in UALIFT1 gives
\[
 \boxed{H_A\simeq
 \operatorname{Isom}_{T_N\text{-alg}}(B_N,A\otimes O_{T_N}).}
 \tag{17}
\]
For N7 the canonical normalization condition is d1^3=d3*d4^5; for N>=8 it is d3*d4=d1*d5. These identities follow back from E, proving both directions as functors. Thus H_A->T_N is finite etale of degree N!, a torsor under the group scheme Aut(A).

The quadratic lift space is k-rational exactly when it has a k-point. The forward implication uses a birational affine chart over an infinite field; the reverse uses projection from a smooth isotropic point. In the realizable case H_A(k) and its realizing targets are Zariski dense. Over number fields the rational quadric also gives weak approximation for compatible local targets, by lifting to H_A and approximating on a rational chart. No assertion about arbitrary integral targets or strong approximation is made.

This construction is functorial for finite etale families and base change. Its smooth geometric fibres supply etale-local choices of generators. A point in a residue-field quadric alone is not asserted to lift Zariski locally: the earlier open-projective-bundle lifting argument had no quadratic equation.

**UADOUBLEMOD1.** In characteristic zero the displayed degree-N three-dimensional compiler slices form a pairwise stably inequivalent family of dimension N-6, N>=7. Every realizable fixed A occurs densely among those parameters.

For N7 the parameter space is Gm_h. The invariant a7^5/a4^8=h, and UASPARSE1 separates its points. For N>=8 the parameter space is Gm^2 x A^(N-8), with nonzero w,aN, and its dimension is N-6. The equal nonzero a4,a5 force any residual scaling to be one, so UASPARSE1 again separates all parameters, including their coefficient-zero boundaries. Absolute stable atomicity of each admissible slice is already AMS1.

The map H_A->T_N is finite etale surjective. The target-forgetting projection T_N->Y_N is smooth surjective of relative dimension three. Thus H_A->Y_N has the same property. If A is realizable, its rational lift space has dense k-points, giving the claimed dense parameter image.

Consequently double trace normalization still does not make this compiler canonical, even on the base fields where it is field-versal and rank seven achieves the global target-dimension minimum. Its surviving stable-class defect is N-6, positive already at rank seven. This is not an obstruction to every other compiler.

## 5. A uniform limit of linear normalization

There is only one proper natural homogeneous linear condition on a split generator that can contain primitive generators. Over any field and N>=3, every S_N-invariant linear subspace W of k^N is one of
\[
 0,\quad k(1,\ldots,1),\quad
 \{\sum r_i=0\},\quad k^N.                            \tag{18}
\]
If W contains a nonconstant vector v, a transposition difference is a nonzero multiple of e_i-e_j. All such differences lie in W by permutation invariance and span the trace hyperplane; W is that hyperplane or all k^N. If no such v exists, W lies in the scalar line. This proof also covers characteristics dividing N, when the scalar line lies inside the hyperplane.

Thus a functorial base-change-compatible linear generator subspace of codimension at least two cannot supply primitive generators. Further compression must use a nonlinear constraint or a different receiver. The trace quadratic criterion above is the first such constraint, and shows exactly why arithmetic enters.

The weight-lattice rule, trace-quadric criterion, binary-block sufficient test and stable-boundary fingerprint play different roles. They respectively decide root-free scaling, realizability for an individual algebra, universality over all field extensions, and whether choices yield the same stable map. None can replace the others.

## 6. Literature inputs, assurance and next gates

The [literature and gap map](LITERATURE_AND_GAPS.md) records the primary sources used here. Duncan's lower bound, Springer's odd-degree theorem, local quadratic-form isotropy and Hasse--Minkowski are imported results. The new contributions here are their explicit Keller-map consequences, the exact full-fibre criterion with denominator control, the rank-seven polynomial section, and the resulting exact target-dimension calculation. Priority and independent review are unclaimed.

The [checker](verify.py) retains bounded exact controls of the polynomial Jacobians, generic section identities, complete fibre reconstruction, both trace equations in odd characteristics including characteristics dividing N, the rank-seven binary isotropic vector, quadratic radicals and failure controls. It does not prove universal monodromy, all-field isotropy, the literature theorems or stable classification by finite samples.

The overarching gaps are now more sharply separated:

- Over arbitrary characteristic-zero base fields, the universal N-2 construction still supplies the common bound. D_N is available in N-3 coordinates with the exact trace-form obstruction above.
- When N is invertible, field-versality of D_N is exactly (15); its negative cases now require another construction. When the characteristic divides N, (15) remains sufficient with the nonconstant-vector condition, and a complete converse on the smooth locus is not asserted. A fixed-field u-invariant bound is insufficient in either case.
- Rank-five A3, rank-six A3, and N-4 affine Keller constructions remain open here.
- A canonical quotient receiver should retain the full algebra automorphisms and must carry effective polynomial landing; taking a root/configuration quotient does not by itself supply affine source and target spaces or a natural section.
- The positive stable-class defect persists in the doubly normalized family. Constructing a different presentation-invariant compiler remains separate from reducing universal dimension.

The active goal remains broader than these results; no claim is made that the programme's global open problems have all been closed.
