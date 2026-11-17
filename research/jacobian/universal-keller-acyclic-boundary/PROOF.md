# Canonical cyclic receivers and acyclic degree obstructions

14 November 2026. The [fixed-fibre degree budget](../universal-keller-fixed-fibres/PROOF.md) leaves generic degree at least N+2 for canonical characteristic-zero transitive fibres. We strengthen this to N+3 in every even symmetric rank, and in the affine-action category in every rank. The abstract acyclic-surface obstruction itself is attainable: published Dubouloz--Palka deformations give explicit missing fibres, whose whole omission sets we calculate. These are not whole affine-space constructions.

Work geometrically over an algebraic closure. Generic degree always means geometric function-field degree. Acyclicity means the stated etale cohomology of a point, not recognition as affine space.

## 1. UAACYCLICGAL1: Galois subextensions come from actual etale covers

Let f:S->T be a separated dominant generically finite etale morphism of smooth connected varieties. Assume
\[
 \Gamma(S,\mathcal O_S)^*=k^*,\qquad \operatorname{Pic}(T)\text{ is torsion}.
 \tag{1}
\]
This general statement works in every characteristic. Then:

1. T minus f(S) has codimension at least two.
2. Every finite Galois intermediate field k(T)<=E<=k(S) has normalization C->T finite etale; f lifts to C.
3. If H1(T,Z/ell)=0, no such nontrivial intermediate Galois group can have a quotient C_ell. In particular it cannot be a nontrivial ell-group.

For (1), f(S) is open. If an irreducible divisor D misses it generically, it misses it entirely. Since T is regular, D is Cartier; a positive multiple mD is principal because Pic(T) is torsion. Its defining rational function has no pole on normal T, so is regular and nonconstant. Pullback is a nowhere-zero function on S, hence constant, contradicting dominance. Thus no divisor is omitted.

For (2), at the generic point of every target divisor choose an actual source point above it. The associated valuation of k(S)/k(T) is unramified, including separable residue extension, because f is etale. Its restriction to E is unramified. Since E/k(T) is Galois, conjugacy gives the same conclusion at every place above that divisor. Normalize T in E. This is finite, normal and etale at all codimension-one points; Zariski--Nagata purity, [Stacks0BMB](https://stacks.math.columbia.edu/tag/0BMB), makes it etale everywhere. The lift S->C is regular by integral closure and normality of S. Assertion(3) follows because a C_ell quotient would give a nonzero class in H1(T,Z/ell).

Equivalently, let G be the generic transitive monodromy and P a point stabilizer. The maximal Galois intermediate extension has group
\[
                 G/\langle\!\langle P\rangle\!\rangle.
 \tag{2}
\]
It must be a quotient of the actual etale fundamental group of T. This extends the affine-space normal-generation rule KMFI1 to targets with nontrivial fundamental group, rather than assuming that group away. The classical Galois Keller case is recovered when T is characteristic-zero affine space.

### Automatic hypotheses for characteristic-zero acyclic affine varieties

If Z is smooth affine and rational ell-adic acyclic in characteristic zero, then its units are constants and Pic(Z) is finite. If Z is mod-ell acyclic, that finite group has order prime to ell. This extends the integral supplement in UAACYCLICRIG1 to the rational hypothesis needed here.

Indeed, compactify smoothly. Units modulo constants form a finite free group through their boundary divisors. Its ell-adic completion injects into H1(Z,Z_ell(1)) by Kummer, so rational acyclicity kills its free rank. Restriction H1(Zbar,Q_ell)->H1(Z,Q_ell) is injective; hence Pic^0(Zbar)=0. Neron--Severi finite generation makes Pic(Zbar) finitely generated, and restriction makes Pic(Z) a quotient. Its ell-adic completion injects by Kummer into H2(Z,Z_ell(1)); rational vanishing kills its free part. Thus Pic(Z) is finite. Mod-ell H2=0 further kills its ell-primary part. These use finite-type cohomology and the retained Kummer and Theorem-of-the-Base inputs, not an affine-plane classification.

### Quadratic consequences

Every separable degree-two extension is Galois. Thus an etale map of generic degree two satisfying(1) has finite etale double normalization C->T. If T has no etale double cover, such a map cannot exist, whether surjective or not. For characteristic-zero acyclic affine T, Kummer identifies double covers with Pic(T)[2], since its units are constants.

Consequently a degree-two map between smooth affine mod-2 acyclic varieties is impossible in characteristic zero, in every dimension. For odd ell, a quadratic map between mod-ell acyclic varieties instead requires a nonzero class
\[
 0\ne L\in\operatorname{Pic}(T)[2],\qquad f^*L=0.
 \tag{3}
\]
The pullback vanishes because S lifts to the double cover. This is a concrete obstruction that rational acyclicity alone does not see.

## 2. UAACYCLICATOM1: primitive monodromy gives stable atomicity on acyclic varieties

Let Z be a smooth geometrically connected affine variety with geometrically constant units, torsion canonical bundle, and geometrically mod-ell or rational ell-adic acyclic cohomology. Any etale endomorphism with primitive geometric monodromy is absolutely stably atomic **among regular self-maps of this same variety**: after every extension of constants and every identity stabilization Z times A^r, a factorization into regular self-maps has an automorphism factor. In characteristic zero all the hypotheses besides acyclicity are automatic by Section1.

Choose m with K_Z^m trivial and a nowhere-zero generator omega. For a dominant regular self-map g, write g^*omega=J_g omega using the m-th power of its differential determinant. Chain rule gives J_(h g)=(J_h composed with g)J_g. If h g is etale, the product is a unit. Both factors are units; geometrically constant units and dominance make J_h constant too. Hence both h and g are etale everywhere. Primitive monodromy means that the generic field extension has no proper intermediate field, so one factor is birational. UAACYCLICRIG1 makes that factor an isomorphism. Products with affine spaces preserve every hypothesis. Geometric monodromy and its primitivity are unchanged by constant extension or adjoining independent identity coordinates. In particular every geometric prime-degree etale endomorphism satisfies the criterion, because a transitive action of prime degree is primitive. No converse for imprimitive monodromy is asserted.

The receiver restriction matters. Allowing different intermediate varieties can produce a non-isomorphic birational open-embedding factor; the constructions below explicitly factor that way. This is not atomicity in an unspecified category of all varieties, nor an assertion that Z is affine space.

## 3. UAACYCLICMISS1: published deformations attain exact finite omission sets

For k>=2 define smooth surfaces
\[
 S_k=\{u+u^2v=w^k\},\qquad B_k=\{x^k y=z^k-1\}.
\]
Choose a nontrivial k-th root epsilon. Dubouloz--Palka Corollary5.2 gives an open embedding
\[
 j(u,v,w)=(w,vR(-uv),1-(\epsilon-1)uv),\qquad
 R(t)={ (1+(\epsilon-1)t)^k-1\over t(t-1)},
 \tag{4}
\]
and a finite etale quotient of degree k,
\[
                 \pi(x,y,z)=(x^k,y,xz).
 \tag{5}
\]
The omitted lines of j are y=0,z=epsilon and x=0,z=zeta for zeta^k=1, zeta not1 or epsilon, as computed in the earlier fixed-fibre proof.

The paper's equation(5.9), Proposition5.5 and Example5.7 give the shear method. For any polynomial Q define the polynomial automorphism of B_k
\[
 \alpha_Q(x,y,z)=\left(x,
 y+{(z+x^kQ(x^k))^k-z^k\over x^k},
 z+x^kQ(x^k)\right).
 \tag{6}
\]
The fraction is polynomial and the inverse is alpha_(-Q). Therefore
\[
                       \eta_{k,Q}=\pi\alpha_Q j
 \tag{7}
\]
is an etale endomorphism of S_k of generic degree k. The construction and its published deformation framework are credited to Dubouloz--Palka; the exact image calculation below is our deduction, with no priority certification.

Put
\[
 h(U)=\epsilon+UQ(U),\qquad b(U)={h(U)^k-1\over U}.
 \tag{8}
\]
Both are polynomial, and b(0)=k epsilon^(k-1) Q(0). The moved horizontal omitted line is parametrized by
\[
 (x,y,z)=(x,b(x^k),h(x^k)).
 \tag{9}
\]
The other k-2 omitted lines remain the same whole x=0 lines. We obtain the exact whole image
\[
 \boxed{\eta_{k,Q}(S_k)=S_k\setminus
       \{(a,-1/a,0):h(a)=0\}.}
 \tag{10}
\]
Here h(0)=epsilon!=0, so all a are nonzero and the omitted set is finite.

To prove this and count whole fibres, use the finite etale cover pi. At U!=0 its k points have x^k=U, z=W/x and y=V. If h(U)!=0, at most one lies on(9); precisely one does when V=b(U). If h(U)=0, all k lie on(9) exactly when W=0, equivalently V=-1/U; otherwise none do. At U=0 the sheets have x=0,z^k=1. The sheet z=1 is always retained. The sheet z=epsilon is retained unless V=b(0), and the other sheets are removed. Thus the complete geometric ranks are
\[
\begin{array}{c|c}
 U=0,\ V=b(0)&1\\
 U=0,\ V\ne b(0)&2\\
 U\ne0,\ V=b(U),\ h(U)=0&0\\
 U\ne0,\ V=b(U),\ h(U)\ne0&k-1\\
 U\ne0,\ V\ne b(U)&k.
\end{array}                                                   \tag{11}
\]
Etaleness makes all fibres reduced. This is a classification of the entire geometric image, not a search for rational solutions.

Given any finite subset A of nonzero scalars, take
\[
 h(U)=\epsilon\prod_{a\in A}(1-U/a),\qquad Q(U)=(h(U)-\epsilon)/U.
 \tag{12}
\]
This prescribes exactly the omitted set in(10), with generic degree and source surface fixed. For A empty the map is the earlier surjective undeformed example; every nonconstant h instead omits a target.

### A rational quadratic example with an empty fibre

For k=2, epsilon=-1 and Q=1, put z_0=1+2uv. Then
\[
 \eta(u,v,w)=\left(w^2,\ 4v+2z_0+w^2,\ w(z_0+w^2)\right)
 \tag{13}
\]
is defined over Q and its image is exactly S_2 minus (1,-1,0). Indeed a hypothetical preimage has w^2=1, z_0=-1 and v=0, contradicting z_0=1+2uv. This map is etale on the whole surface, not merely on the finite chart.

### Sharp degree thresholds and stable inequivalence

The surfaces S_k are Q-homology planes with fundamental group C_k. Their H1 is Z/k and higher reduced integral homology vanishes: a smooth affine complex surface has a finite CW model of dimension at most two, whose H2 is free, and rational acyclicity kills that free group. Hence S_k is mod-ell acyclic exactly when ell does not divide k.

Together with UAACYCLICGAL1 and UAACYCLICRIG1, equation(13) and the k=3 instance of(7) prove the exact geometric minimum generic degree for a non-surjective etale map between smooth affine mod-ell acyclic varieties of any fixed dimension d>=2:
\[
 \boxed{2\text{ if ell is odd},\qquad3\text{ if ell=2}.}        \tag{14}
\]
For the cubic witness work over an algebraically closed characteristic-zero field, or one containing epsilon; no cubic construction over every smaller constant field is claimed. Products with A^(d-2) give every d>=2. In dimensions zero and one no such map exists.

At fixed prime k, equations(7),(12) give infinitely many absolutely stably atomic self-maps of one fixed surface, by UAACYCLICATOM1. Taking A={1,...,m}, their omission sets have m points. After r identity coordinates the complement has exactly m connected components, each A^r. Target automorphisms preserve that number, so different m remain inequivalent under polynomial left-right equivalence after every identity stabilization. This supplies a simple stable inequivalence certificate, independent of the paper's stronger parameter-by-parameter surface inequivalence proof. The number of omitted targets is unbounded at fixed generic degree.

None of these stabilized sources is affine space: the geometric fundamental group C_k persists under identity stabilization. They solve the abstract acyclic fixed-surface existence gate, not its realization as the required whole fixed schemes inside an equivariant affine-space Keller construction.

## 4. UAACYCLICRECV1: canonical cyclic receivers on forms of one acyclic surface

Assume the characteristic-zero base contains a primitive k-th root epsilon. Every map eta_(k,Q) above is equivariant for the actual source and target actions
\[
 (u,v,w)\longmapsto(u,v,\zeta w),\qquad \zeta\in\mu_k.
 \tag{15}
\]
The first two outputs are invariant, and the third has weight one. This is a regular whole-surface action, not the generic deck action of the finite cover. Choose c in the base with c!=0 and h(c^k)!=0, once for the construction. The fixed target
\[
                       y_c=(c^k,-c^{-k},0)
 \tag{16}
\]
has complete split fibre w^k=c^k, a regular mu_k set. The projection w/c identifies it equivariantly with mu_k. For Q=1, c=2 works for every k in characteristic zero.

Let a mu_k torsor over a base B correspond to a line bundle L with a specified isomorphism phi:L^k->O_B. Its associated surface is the affine hypersurface bundle
\[
 \mathcal S_L=\{(u,v,w)\in\mathcal O_B^2\oplus L:
                         u+u^2v=\phi(w^k)\}.
 \tag{17}
\]
Put s=phi(w^k), z_0=1-(epsilon-1)uv and
\[
 A_k(z,s)={ (z+sQ(s))^k-z^k\over s}.
\]
The quotient is polynomial. The formula
\[
 (u,v,w)\longmapsto
 \left(s,\ vR(-uv)+A_k(z_0,s),\ w(z_0+sQ(s))\right)
 \tag{18}
\]
glues on L and defines a relative etale endomorphism of S_L, with generic degree k on every geometric fibre. It becomes the published deformation after an etale trivialization of the torsor; hence its etaleness and the complete fibre identification descend effectively. The section(16) has fibre
\[
 \{w\in L:\phi(w^k)=c^k\}\ \simeq\
 \{\ell\in L:\phi(\ell^k)=1\},\qquad w\mapsto w/c,
 \tag{19}
\]
which is exactly the input torsor. Thus(17)--(19) construct a canonical marked lift from B_mu_k into etale self-maps of forms of one smooth rationally acyclic affine surface. For prime k it is absolutely stably atomic by UAACYCLICATOM1. For k=2 this is an explicit characteristic-zero BS2 construction in the enlarged receiver.

The geometry explains why the affine-space no-go theorem UAFIXDEFECT1 does not apply. The whole source and target mu_k fixed schemes are
\[
       \{w=0\}=\{u=0\}\ \sqcup\ \{uv=-1\}
                         \simeq\mathbb A^1\sqcup\mathbb G_m.
 \tag{20}
\]
Every fixed source point maps into the A1 component, since the first output w^k is zero. The marked target lies in the other component. The ambient surface is rationally acyclic but not mod-p acyclic for primes p dividing k. Rational acyclicity alone therefore cannot replace the tame Smith hypothesis in the canonical obstruction. These source surfaces, including their twists, never become affine spaces after identity stabilization.

### A complete quadratic formula over an arbitrary characteristic-zero field

For the torsor K[sqrt(a)], write the twisted surface as
\[
                   S_{2,a}: a(u+u^2v)=z^2,
                     \qquad a\in K^*.
\]
For Q=1 the descended endomorphism is
\[
 (u,v,z)\longmapsto\left(
 z^2/a,\ 4v+2(1+2uv)+z^2/a,\ z(1+2uv+z^2/a)\right).
 \tag{21}
\]
At the fixed target (4,-1/4,0), its whole fibre is
\[
               u=-8/3,\qquad v=15/16,\qquad z^2=4a.
 \tag{22}
\]
The coordinate z/2 identifies it with K[T]/(T^2-a), including the split quadratic algebra. Changing the Kummer presentation is an isomorphism of (L,phi); the line-bundle construction, not a chosen square root, supplies the functorial object. This is a map on a hypersurface, not a claimed Keller self-map of A3.

### The coarse affine quotient is never Keller

The categorical quotient S_k/mu_k is A2 with coordinates(u,v), because w^k=u+u^2v. The induced polynomial map is
\[
 \bar\eta(u,v)=\bigl(s,\ vR(-uv)+A_k(z_0,s)\bigr),
                 \quad s=u+u^2v.
\]
Its exact Jacobian is
\[
 \boxed{\det D\bar\eta=k(1-\epsilon)
                         (z_0+sQ(s))^{k-1}.}                 \tag{23}
\]
For a short proof, use the nowhere-zero residue form omega=du wedge dw/u^2 on its dense chart. Direct differentiation gives eta^*omega=k(1-epsilon)omega. The quotient map q:S_k->A2 satisfies q^*(du wedge dv)=k w^(k-1)omega. Apply this identity to q eta=bar_eta q and cancel on the dense unramified chart; equality of regular functions extends(23) everywhere.

The factor z_0+sQ(s) equals1 on u=0 and epsilon on uv=-1, so it is nonconstant for every Q. Thus the quotient is never a Keller map. Identity stabilization and polynomial left-right equivalence preserve that failure. For the quadratic Q=1 example it is
\[
 (u,v)\longmapsto(u+u^2v,\ 2+u+(u+2)^2v),
\]
with determinant4(1+2uv+u+u^2v). Etaleness of the quotient-stack map does not imply etaleness of its coarse map: source stabilizers are trivial on the marked fibre, while the target stabilizer is all mu_k. Simply forgetting the cyclic cover cannot turn this receiver into an affine-space Keller construction.

### Fixed universality and exact dimension in this receiver

One fixed eta_(k,Q) on S_k also realizes every mu_k torsor over every extension of the base field. Represent the torsor by x^k=a. Since each such characteristic-zero field is infinite, choose b!=0 with h(a b^k)!=0. At the target (a b^k,-1/(a b^k),0), equation(11) gives the complete degree-k fibre x^k=a b^k, which is the required torsor. This is a field-versal fixed map, distinct from the canonical family on twisted surfaces.

Within the category of geometrically rationally acyclic smooth affine sources and targets, dimension two is minimal for these nontrivial cyclic marked constructions: a geometrically rationally acyclic smooth affine curve in characteristic zero is A1 geometrically, and an etale map between such curves has degree one. The same lower bound applies to fixed field-versal maps in this category. For prime k, dimension two is also attained with absolute stable atomicity. For composite k no atomicity assertion is made for eta_(k,Q).

This settles the cyclic receiver in the specified enlarged category, not canonical marked affine-space descent for S_N in general. The [torus and symmetry follow-up](../universal-keller-torus-receivers/PROOF.md) rules out a full natural S_N action, N>=3, on the same unstabilized pseudo-plane: its finite automorphism groups are cyclic. Projection-preserving identity stabilizations fail too, and cyclic generic degree N>=4 rules out every stabilized natural S_N lift. Cubic stabilizations mixing the surface projection and other ambient geometries remain separate. Torus receivers instead give actual natural S3 and S4 constructions in their stated category.

## 5. UAPARITYDEF1: stronger canonical degree costs

Retain the notation and hypotheses of UAFIXDEFECT1: whole affine-space equivariant Keller map, complete transitive marked rank N at a fixed target, generic degree M and delta=M-N.

In characteristic zero, if a 2-subgroup H has no fixed label, its whole fixed varieties are smooth affine mod-2 acyclic. Their fixed restriction has empty marked fibre. UAACYCLICGAL1 excludes degree two and UAACYCLICRIG1 excludes degree one; the degree budget therefore gives delta>=3. For natural S_N with N even, take N/2 disjoint transpositions. Thus arbitrary polynomial symmetric actions satisfy
\[
 \boxed{M\ge N+3\quad(N\text{ even}),\qquad
        M\ge N+2\quad(N\text{ odd}).}                         \tag{24}
\]
The odd bound is the retained UAFIXDEFECT1 bound. If its equality M=N+2 occurs, every deranging tame prime subgroup has a target fixed variety with nonzero Pic[2] killed by the fixed restriction, as in(3). The deformed S_2 maps prove that this abstract surface gate is attainable; their occurrence as actual whole affine-space fixed loci with simultaneous symmetry is not proved.

For affine source and target actions in characteristic zero, take the full group. Averaging gives a source fixed point, both whole fixed schemes are affine spaces, and a transitive fibre of rank N>1 has no full-group fixed point. A degree-two restriction is impossible because its normalization would be a finite etale double cover of affine space. Hence
\[
       \boxed{M\ge N+3\text{ for every nontrivial transitive affine
                       canonical construction in characteristic zero}.} \tag{25}
\]
More generally, over every algebraically closed field of characteristic different from two, no affine-space Keller map has generic degree two: apply Section1 with Pic=0, then Kummer H1(A^t,mu2)=0. Therefore the tame cyclic derangements from UAFULLAFF1 improve the affine natural S_N surplus to three in all characteristics other than two, outside the sole full-degree exception (N,char)=(3,3). In characteristic two retain the earlier bound and its two affine exceptions; no extra assertion about degree-two wild restrictions is made.

These are necessary generic-degree bounds, not constructions or new unmarked dimension minima. Atomicity of an ambient higher-degree canonical map remains separate from the atomicity of a fixed restriction or of a pseudo-plane example. Full-degree stabilization remains excluded, and OP-KDESC remains open for effective affine-space boundary realizations and different receivers.

## 6. Evidence and remaining gates

The all-rank arguments are the divisor/unit proof, purity, the pluricanonical chain rule, the cover image calculation and the fixed-locus degree budget. Exact controls check small all-k instances, complete fibre ideals including empty fibres, the deformation automorphisms, explicit parity derangements and omission/component invariants. They do not formalize purity, compactification, cohomology, or global equivariant gluing. Independent review and publication-priority checks remain pending.

The immediate remaining questions are whether these nontrivial finite fundamental groups and Picard torsion can occur in the required actual tame fixed loci, whether the restrictions can be assembled compatibly for all subgroups, and whether an entire affine-space source with the increased generic degree can be proved. The generic degree threshold for the abstract acyclic map is settled; the canonical ambient realization is not.
