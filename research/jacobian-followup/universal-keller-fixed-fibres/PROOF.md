# Fixed-locus degree defects and the full-degree canonical obstruction

14 November2026. The [Smith extension](../universal-keller-smith/PROOF.md) gave ambient dimension bounds and an acyclic-surface gate. The degree budget below is stronger: a nontrivial transitive canonical fibre cannot have the full generic degree in characteristic zero, in any ambient dimension, and a one-degree surplus is also impossible. Consequently identity stabilization never makes any retained degree-N universal map into a coherent natural marked S_N lift. Higher-degree boundary fibres and other receiver categories remain separate.

## 1. Definitions and a normalization budget

Let F:X=A^D -> Y=A^D be a polynomial Keller map, with geometric generic degree M. A finite group G acts by polynomial automorphisms on source and target, and F is equivariant. At a fixed geometric target y, its complete fibre is a finite-etale G set Omega of rank N. Throughout, the source and target are the actual whole affine spaces, not finite-chart opens. Work geometrically over an algebraic closure.

For H<=G assume the **whole** X^H and Y^H are nonempty smooth connected varieties. The restricted map F_H is etale. To see this even for wild H, lift an H-fixed infinitesimal point along the etale map F. Uniqueness makes its lift H fixed, so the fixed-scheme restriction is formally etale; finite presentation makes it etale. The target contains y and both fixed varieties have the same dimension, so F_H is dominant and generically finite. Put
\[
 r_H=|\Omega^H|,\qquad d_H=\operatorname{gdeg}(F_H),\qquad
 \delta=M-N.
\]

**Fixed-locus degree budget.** Under these hypotheses, in every characteristic,
\[
 \boxed{r_H\le d_H\le r_H+\delta.}                 \tag{1}
\]
The budget also holds for any equivariant separated dominant etale morphism between smooth connected varieties; affine space is only needed for the later Smith conclusions. No flatness of the finite normalization at y is assumed. This matters in high ambient dimension.

### Strict-local proof, including the omitted branches

Normalize the normal affine target in k(X), obtaining a finite normal map B->Y. Zariski Main embeds X as an open in B. Pass to the strict henselization A of O_(Y,y), a normal local domain. The N etale source points split off N rank-one factors:
\[
       B_A=\left(\prod_{\omega\in\Omega} A\right)\times R,
       \qquad\operatorname{rank}_A R=M-N.          \tag{2}
\]
Each factor through an actual source point is wholly in the open X: an open of a local finite branch containing its closed point is the whole branch. It is finite etale over the strictly henselian base and therefore is a copy of A. Every remaining factor has positive generic rank. The strict-local normality and torsion-free rank assertions follow by etale base change of normalization; they do not require B to be locally free.

H acts on A and permutes the N copies as on Omega. The H-fixed base is the corresponding strict-local piece of Y^H. It is a regular local domain; equivariant etale lifting identifies this fixed piece with the fixed-scheme base change. Over its geometric generic point, precisely r_H of the actual N copies are H fixed, and all those copies lie in X_H.

The remaining finite algebra R has at most delta geometric points over any point of Spec A. Indeed localize at that point and strictly henselize again; each nonzero local factor has positive generic rank, and their ranks sum to delta. This familiar point-count bound for finite algebras over a normal base is distinct from a bound on special-fibre scheme length, which may exceed generic rank. Therefore at most delta additional H-fixed source points can appear. The geometric generic fibre of F_H has d_H points, proving (1).

In particular, if delta=0, every H-fixed target in this finite neighbourhood has r_H=d_H>=1. Equivalently one may restrict the finite etale chart at a full-size fibre to Y^H: its H-fixed subcover has constant rank on the connected open Y^H, and dominance prevents rank zero.

## 2. UAACYCLICRIG1: birational etale maps between acyclic varieties

Let ell be invertible in the field. If S and T are smooth geometrically connected varieties of the same dimension and both have the mod-ell etale cohomology of a point, then every birational etale morphism
\[
                         S\longrightarrow T
\]
is an isomorphism. Affineness, equality of the varieties and a Picard-group hypothesis are unnecessary.

By Zariski Main, a birational etale morphism into normal T is an open immersion j. If its complement Z is nonempty, compactly supported localization gives
\[
       H_c^i(S)\longrightarrow H_c^i(T)
                    \longrightarrow H_c^i(Z)\longrightarrow H_c^{i+1}(S).
                                                               \tag{3}
\]
Poincare duality makes the first two varieties' compactly supported cohomology zero except in degree2d, where each is one copy of Z/ell(-d). The top-degree map j_! is an isomorphism: by duality it is the dual of the identity restriction on H0. Thus every H_c^i(Z) vanishes. But a nonempty finite-type variety Z of dimension m has nonzero H_c^(2m), with a fundamental class for each top-dimensional irreducible component. Indeed remove a lower-dimensional closed subset to obtain disjoint smooth dense opens in the top components; compact-support localization and the dimension bound identify their nonzero top groups with H_c^(2m)(Z). This contradiction proves the assertion. The argument is valid in positive characteristic as well; it does not use finite-etale Euler multiplicativity there.

The same argument works with Q_ell coefficients: rational ell-adic acyclicity already suffices for the birational rigidity assertion. The mod-ell hypothesis is only needed for the additional integral constraints below. Consequently a non-surjective etale map between smooth connected mod-ell or Q_ell acyclic varieties has generic degree at least two. This is not a properness theorem for maps of larger degree.

### Extra characteristic-zero geometric constraints

For such a smooth affine variety S in characteristic zero, its units are constants and Pic(S) is finite of order prime to ell. Here is a standard compactification proof, useful when checking actual fixed surfaces.

Choose a smooth projective compactification with boundary divisors. Units modulo constants inject via their boundary divisors into a finite free abelian group. Kummer theory and H1(S,mu_ell)=0 make this free group ell-divisible, hence zero. Restriction H1(Sbar,mu_ell)->H1(S,mu_ell) is injective: every connected finite etale cover of the normal compactification stays connected over the dense open. Thus Pic^0(Sbar) has no ell-torsion and has dimension zero. The Neron-Severi theorem makes Pic(Sbar), and hence its quotient Pic(S), finitely generated. Finally Pic(S)/ell injects into H2(S,mu_ell)=0. A finitely generated ell-divisible abelian group is finite of order prime to ell. The primary Kummer and Neron-Severi inputs are pinned; no affine-plane recognition follows from them.

## 3. UAFIXDEFECT1: every transitive canonical fibre needs a genuine degree surplus

For a tame ell-subgroup H, Smith theory supplies exactly the smooth connected acyclic fixed varieties needed above. Combining (1) with UAACYCLICRIG1 gives
\[
 \boxed{|\Omega^H|=0\quad\Longrightarrow\quad
              M-N\ge d_H\ge2,\qquad \ell\ne\operatorname{char}k.} \tag{4}
\]
Thus a fixed-point-free tame prime-power symmetry forbids both a full-size fibre and a one-degree surplus. This concerns **regular source and target symmetries**, not merely generic monodromy.

In characteristic zero, every nontrivial finite transitive G set has a prime-power derangement by Fein--Kantor--Schacher, Theorem1. Hence
\[
 \boxed{M\ge N+2\quad\text{for every such canonical transitive fibre},
         \qquad N>1.}                              \tag{5}
\]
The general finite-group conclusion imports that CFSG-dependent published theorem. For the natural S_N set no CFSG is needed: choose a prime ell dividing N and the product of N/ell disjoint ell-cycles. Its cyclic subgroup has no fixed label. Therefore (5) holds for natural S_N in every rank and every ambient dimension by an elementary group argument.

In characteristic p, the same conclusion holds whenever the label action has a prime-power derangement of order prime to p. In particular for natural S_N, a canonical full-degree construction with arbitrary polynomial actions can exist only when N is a power of p. This is necessary, not sufficient.

### What the first possible degree must contain

If M=N+2, equation(4) forces d_H=2 for every tame prime-power derangement subgroup H. Thus the first surviving degree requires simultaneous degree-two etale restrictions with empty marked fibres between the actual acyclic fixed varieties. Their generic quadratic involutions cannot extend to regular deck automorphisms in characteristic zero: etale deck freeness would give a finite etale quotient of degree two, whereas compact-support Euler characteristic would read1=2*chi(quotient). This uses characteristic-zero Euler multiplicativity only. In particular generic quadratic Galois symmetry is not enough; one must examine its extension across the nonfinite source boundary.

### Consequences for stabilization and the existing universal maps

A full-degree coherent natural marked lift from BG to polynomial affine-space Keller maps is impossible in characteristic zero for any nontrivial finite transitive label set. The same holds for any representative of the polynomial left-right class of any identity stabilization: those operations preserve M and N.

Every retained fixed universal atomic symmetric map UAM1, UACOMP1, UATRACE1, UASEXT1 and the double-trace maps has generic degree N. **No number of identity coordinates, and no new nonlinear polynomial actions on a stabilized representative, can make it a canonical natural marked S_N construction.** The previous ambient lower bounds remain valid for possible higher-degree boundary constructions, but paying those coordinate costs cannot repair the degree-N category.

The full-degree normal-monodromy and automatic-atomicity theorem UACANATOM1 remains a valid conditional statement. Its characteristic-zero case now has no realizations. It cannot provide automatic atomicity for the remaining characteristic-zero boundary candidates M>=N+2; their generic monodromy and atomicity require separate analysis.

A surviving characteristic-zero construction must increase generic degree and realize the rank-N algebra at a genuinely nonfinite boundary target, or change the receiver category. Equation(5) measures generic-degree surplus, not the number of distinct omitted normalization points: several generic sheets may merge at one omitted ramified point. No attainment of N+2 is claimed.

## 4. UAFULLAFF1: exact full-degree affine symmetric exceptions

For affine source and target actions, any subgroup with a nonempty source fixed locus has affine-space fixed loci on both sides. The fixed-locus budget (1) applies even to wild subgroups; no averaging is needed once the affine fixed locus exists. A birational Keller map of affine spaces is an automorphism in every characteristic. Therefore a fixed-point-free label subgroup with a source fixed point again forces M-N>=2.

In particular a full-degree nontrivial transitive affine-bundle construction has **no source G-fixed origin**, in any characteristic: take H=G otherwise. This is a general necessity, not just a feature of the Frobenius examples. It persists under vector-representation stabilization.

For natural S_N, the full-degree affine category is nonempty exactly in these three cases:
\[
 \boxed{(N,\operatorname{char}k)=(2,2),\ (3,3),\ (4,2).}         \tag{6}
\]
Their canonical affine minima are1,1,2, respectively, over every field of that characteristic. Matching maps are x-x^2, x-x^3 and
\[
                    (x,y)\longmapsto(x-x^2,y-y^2).
\]
The last has the full natural S4 action through AGL2(F2). Its plane model decomposes; the established canonical absolute-atomicity interval3..6 is unchanged. The prime-degree line maps are absolutely stably atomic.

### Exhaustion by tame cyclic derangements

A cyclic affine subgroup whose order is invertible has a source fixed point by averaging. A prime-to-characteristic derangement therefore excludes M=N and M=N+1.

Outside (6), such a symmetric derangement always exists:

- If the characteristic does not divide N, an N-cycle works, including characteristic zero.
- If p is odd, p divides N and N>=5, use disjoint cycles of lengths2 and N-2. Both lengths are prime to p.
- If p=2, N is even and N>=6, use lengths3 and N-3, both odd.
- In ranks2,3,4 outside the displayed exceptions use respectively a2-cycle, a3-cycle or two disjoint2-cycles.

This proves the exhaustion without classifying higher-dimensional affine representations. In fact every affine canonical natural symmetric realization outside (6), whether full degree or a boundary model, must have generic degree at least N+2. No higher-degree boundary realization is supplied.

### The remaining modular full-degree gate is global disconnectedness

Consider arbitrary polynomial actions in characteristic p, N=p^a>=5, with full generic degree N. A prime-to-p cyclic derangement still exists: lengths2,N-2 for odd p; lengths3,N-3 for p=2 (then N>=8). Its order has at least two distinct prime factors, so it is not a tame prime-power group.

Let C be this cyclic subgroup. Haution's Theorem3.3.1 and its Lefschetz proof give nonemptiness and Euler characteristic one for X^C. Both X^C and Y^C are smooth because C is tame. If Y^C were connected, every nonempty source component would map etale and dominantly to that smooth irreducible target fixed variety. A full-degree finite chart at y has no C-fixed sheet; constant rank on its intersection with Y^C gives a contradiction. Thus
\[
 \boxed{Y^C\text{ must be disconnected, and its component through }y
             \text{ receives no point of }X^C.}                 \tag{7}
\]
The target C action must in particular be non-linearizable. This is a precise global-geometry requirement on the remaining modular power ranks, not a construction or an assertion that such disconnected fixed loci exist in the required category. Prime-power Smith connectedness cannot be extended to this composite cyclic group without proof.

## 5. Published acyclic-surface maps: nonproper does not mean missing a target

Dubouloz--Palka exhibit nonproper etale endomorphisms on pseudo-planes. They do not automatically solve the required empty-fibre surface gate. Their degree-two Example5.3 is a particularly clear control:
\[
 S_2=\{u+u^2v=w^2\},\qquad
 (u,v,w)\longmapsto(w^2,4v,w(1+2uv)).               \tag{8}
\]
Its generic quadratic deck involution is
\[
             (u,v,w)\dashmapsto(-u-1/v,\ v,\ -w).
\]
It preserves the source equation and the map, but has a genuine pole along v=0, where the whole fibre rank drops to one. This explicitly illustrates the distinction between generic Galois symmetry and a regular source action. It is surjective, despite being nonproper. More generally the published Corollary5.2 gives, for k>=2 and a nontrivial k-th root epsilon, the open embedding of
\[
 S_k=\{u+u^2v=w^k\}\quad\text{into}\quad
 T_k=\{x^k y=z^k-1\}
\]
with
\[
 t=-uv,\quad R(t)={ (1+(\epsilon-1)t)^k-1\over t(t-1)},\qquad
 j(u,v,w)=(w,vR(t),1+(\epsilon-1)t).
\]
The numerator vanishes at0 and1, so R is polynomial. The finite etale quotient is pi(x,y,z)=(x^k,y,xz), of degree k. We determine the whole image and fibre ranks directly, rather than infer them from nonproperness.

The complement of j consists of k-1 disjoint affine lines:
\[
 \{y=0,z=\epsilon\}\ \cup\!
 \bigcup_{\zeta^k=1,\ \zeta\notin\{1,\epsilon\}}\{x=0,z=\zeta\}.
                                                               \tag{9}
\]
For R(t)y!=0 the inverse is v=y/R(t), u=-tR(t)/y, w=x. At roots of R, a preimage exists exactly when x!=0; for y=0 and R!=0 it exists exactly on z=1. These cases prove (9).

Consequently eta=pi j is etale of generic degree k and its complete geometric fibre ranks at (U,V,W) in S_k are
\[
 \boxed{\begin{array}{c|c}
 UV\ne0&k\\
 V=0,\ U\ne0&k-1\\
 U=0,\ V\ne0&2\\
 U=V=0&1.
 \end{array}}                                                \tag{10}
\]
For U!=0 the pi fibre has k roots x; for U=0 it has x=0 and z running over the k-th roots of1. Remove exactly the points on (9). This proves all four ranks and, in particular, surjectivity. Fibres are reduced by etaleness, so these are whole fibre ranks, not just found rational roots.

The surfaces S_k are Q-homology planes with H1=Z/k and higher reduced integral homology zero; for every ell prime to k they are mod-ell acyclic. Thus higher-degree nonproper etale maps of acyclic surfaces really do exist, while (8)--(10) still have no empty fibre. Rank loss one is also possible here. There is no contradiction with (4): these surfaces are not the whole affine ambient source with the stipulated regular derangement symmetry, and generic deck permutations are not such symmetry.

The calculation in this section concerns the undeformed maps pi j. The published shear family in Proposition5.5 and Example5.7 gives actual missing fibres, as the [acyclic boundary extension](../universal-keller-acyclic-boundary/PROOF.md) now proves by a complete image calculation. That extension also raises the arbitrary-action even-rank and all affine characteristic-zero degree bounds to N+3. The abstract missing-fibre surface gate is therefore attainable; ambient fixed-locus realization and effective whole affine-space gluing remain open.

## 6. Remaining work and assurance

The characteristic-zero full-degree coherent affine-space receiver category is closed by an obstruction, including arbitrary nonlinear polynomial actions and all identity stabilizations. The broader OP-KDESC remains open for genuinely higher-degree boundary models, correspondence-valued receivers, or other source geometry. Such models must preserve the complete marked fibre, prove the entire affine source and effective gluing, and establish generic monodromy/atomicity independently.

The fixed-locus rank budget and compact-support rigidity are written geometric proofs. Finite controls check actual permutation derangements, the three modular matching actions, exact local rank budgets and the complete pseudo-plane fibre formulas; they do not formalize Smith theory, finite normalization, compactly supported cohomology or CFSG. The general finite-group version uses the original FKS theorem, visually inspected at printed page41. No novelty priority, independent whole-proof review or formal verification is asserted.
