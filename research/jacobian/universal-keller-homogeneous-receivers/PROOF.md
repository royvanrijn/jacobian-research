# Universal fibres, homogeneous receivers and the exact isotropy cost

14 November 2026. There are two all-degree extensions. Full-degree field-versality itself forces geometric symmetric monodromy, so atomicity need not be checked again for each compiler. Separately, finite homogeneous receiver pairs exist in every rank, and their optimal number of geometric stabilizer components can be computed exactly. This does not identify their different source and target varieties, or construct a whole affine-space canonical Keller self-map.

All varieties below are of finite type over a field. Geometric monodromy means the faithful permutation group of the geometric generic fibre. A canonical rank-N realization means compatible S_N actions, a fixed target mark, and its complete natural N-point fibre; twisting this split object must realize arbitrary rank-N finite etale algebras, including disconnected algebras. A fixed universal map means one map whose complete fibres realize all such algebras over every extension field. These are different requirements.

For strong atomicity we allow dominant regular factors through arbitrary normal affine integral varieties. For a self-map we also discuss the smaller category of regular self-map factors on the same stabilized variety. Identity stabilization means product with the identity of affine space. Every theorem states its finiteness and characteristic hypotheses.

## 1. UAUNIVMON1: two fibre types force geometric S_N

Let f:X->Y be separated, dominant and etale between geometrically integral geometrically normal varieties over k, of generic degree N>=2. Suppose:

1. one complete fibre at a k-point is the split rank-N algebra;
2. after some extension K/k, one complete fibre has rank N and its Galois closure acts as the full natural S_N.

Then the geometric generic monodromy of f is S_N. In particular this holds for every full-degree fixed universal rank-N map, over any field and in every characteristic.

### The full-degree finite chart

A complete etale fibre with N geometric points lies in the finite etale locus of f on the target. Write the Zariski-main factorization X open in Xbar finite over Y, with Xbar the normalization in k(X). Over the strict henselization at that target, each of the N etale points splits off a copy of the normal local base from the finite algebra of Xbar. These copies already exhaust generic rank N. A remaining factor has generic rank zero and vanishes because the normalization algebra is torsion-free over the base. Thus the boundary has no point above this target. Its finite image is closed, so removing it gives the asserted chart.

Consequently there is a nonempty open U of Y on which f is finite etale and which contains every complete rank-N fibre. This is the finite-chart step used in [UAPSEUDOSYM1](../universal-keller-torus-receivers/PROOF.md), without assuming a symmetry action. Its inverse image remains geometrically integral.

There is a simpler general lemma here, with no symmetric-group classification. A completely split fibre at a k-point identifies arithmetic and geometric monodromy for ANY finite permutation cover of a geometrically connected k-variety. Indeed this point gives a section of the arithmetic-geometric fundamental-group sequence. Its Galois action on every fibre label is trivial, so its image in the faithful monodromy permutation group is trivial. Every arithmetic element is a product of a geometric element and an element of this section. Hence
\[
                    G_{\rm arithmetic}=G_{\rm geometric}. \tag{0}
\]
Equivalently the connected generic Galois closure is regular over k.

Apply this lemma on U using the split fibre. The second fibre assumption forces the arithmetic group to be S_N: its Galois action is a specialization subgroup. Equation(0) gives geometric S_N immediately. This avoids the longer normal-subgroup/sign argument and works for arbitrary finite permutation groups whenever a corresponding arithmetic upper bound and full-group specialization are available.

There is no inverse-Galois assumption in the universality corollary. Over K=k(t1,...,tN)^{S_N}=k(e1,...,eN), the extension k(t1,...,tN)/K has group S_N in every characteristic, and its point-stabilizer fixed field is a separable rank-N algebra with the required Galois closure. The variables are independent, so the generic roots are distinct. Full-degree universality supplies this fibre and the split k-fibre. If universality was only required over infinite extensions, apply the argument after k(u); geometric monodromy is unchanged.

### Atomicity becomes a consequence of universality

If f is finite and X,Y are affine, [UAFINATOM1](../universal-keller-torus-receivers/PROOF.md) gives absolute stable atomicity in the strong normal-affine category: S_N is primitive, all dominant affine factors are finite etale, and a finite birational factor to a normal variety is an isomorphism.

For a possibly nonfinite self-map of a smooth affine acyclic variety satisfying the constant-unit and torsion-canonical hypotheses of [UAACYCLICATOM1](../universal-keller-acyclic-boundary/PROOF.md), that theorem instead gives absolute stable atomicity among regular self-map factors of the same stabilized variety. In particular every full-degree fixed universal polynomial Keller self-map is absolutely stably atomic, independently of its compiler or a separate profile-monodromy calculation. The hypotheses concern full rank N equal to generic degree. Realizing rank-N algebras only on boundary fibres of a map of degree M>N gives no such conclusion about its degree-M monodromy.

Extending constants and adjoining independent identity variables preserve the geometric Galois group of the given function-field extension. This does not assert that every new finite cover in positive characteristic descends.

## 2. UAHOMOCOST1: monodromy measures disconnected isotropy

Work in characteristic zero. Let the geometric target Y be a homogeneous space L/H of a connected algebraic group. Let C be the image of geometric pi1(L) in geometric pi1(Y). Then
\[
  1\longrightarrow C\longrightarrow\pi_1(Y)
       \longrightarrow\pi_0(H)\longrightarrow1,\qquad
                         C\ \hbox{central}.                \tag{1}
\]
Here pi1 is etale pi1 and C denotes its closed image.

For completeness, spread the finite data to a field embedded in C. For connected complex Lie L, take its universal covering group Ltilde->L; its discrete kernel is central. If Htilde is the inverse image of H, the homotopy exact sequence identifies pi1(L/H) with pi0(Htilde). The kernel of pi0(Htilde)->pi0(H) is the image of the central covering kernel. This gives the central extension topologically; finite quotients and Riemann existence give (1). Characteristic-zero invariance under algebraically closed field extension transfers it to every geometric base. This is the standard homogeneous-space covering argument. [Brion--Szamuely](LITERATURE_AND_GAPS.md) supplies the published connected-stabilizer result; we retain the central-kernel argument explicitly for the disconnected case.

For a connected finite etale cover of Y with monodromy M, the image Z of C is therefore central in M and M/Z is a quotient of pi0(H). If the action of M is faithful and primitive, every nontrivial central subgroup is transitive and regular: its orbits are blocks, and a central element fixing one point fixes all points. Its centralizer in Sym(N) is that same regular abelian group. Hence M=Z in this case, and primitivity forces M=C_p of prime degree. Thus
\[
 \boxed{\text{primitive monodromy other than regular }C_p
          \ \Longrightarrow\ M\text{ is a quotient of }\pi_0(H).} \tag{2}
\]
In particular center-free M is always such a quotient, without a primitivity assumption.

Apply this to a full-degree canonical natural S_N fibre. Equivariance at its fixed target normalizes the geometric monodromy, so M is transitive and normal in S_N. For N>=5, M=A_N or S_N and both are center-free. Consequently
\[
                  N!/2\ \mid\ |\pi_0(H)|.                  \tag{3}
\]
For a fixed universal full-degree cover, section1 makes M=S_N. When N>=3 this group is center-free, and
\[
                    N!\ \mid\ |\pi_0(H)|.                 \tag{4}
\]
These costs concern geometric stabilizer components in a homogeneous receiver, not dimension, essential dimension or the number of coordinates.

If H is connected, pi1(Y) is abelian. Connected covers therefore have abelian regular geometric monodromy. Natural S_N can normalize such a group only for N=2,3,4, with C2,C3,V4 respectively. The [torus constructions](../universal-keller-torus-receivers/PROOF.md) attain all three in characteristic zero. Every quartic connected-isotropy realization has imprimitive V4 monodromy and therefore factors through normal affine intermediate covers when the target is affine. It cannot be strongly absolutely atomic. An atomic canonical quartic finite affine receiver must instead have M=A4 or S4, so
\[
                         12\mid|\pi_0(H)|.                \tag{5}
\]
No additional ambient group dimension or affine identity factor can remove these bounds: pi1(Y times A^r)=pi1(Y) in characteristic zero. The assertions are for finite covers, or for a finite normalization genuinely etale on the entire stated homogeneous target. A nonfinite map can have a branched normalization and is outside this argument.

In positive characteristic the published prime-to-p theorem controls quotients of pi1 whose GROUP ORDER is prime to p. Merely assuming p does not divide the fibre degree is insufficient: a degree-prime-to-p cover can have a Galois closure divisible by p. We make no all-characteristic extension of (1)--(5).

## 3. UAQUOTRECV1: one fixed finite receiver for every transitive group

Let Gamma<=S_N be any faithful transitive finite permutation group, over any field k. Set m=N+1 and embed Gamma in L=SL_m by
\[
        \rho(\sigma)=\operatorname{diag}(P_\sigma,
                                             \det P_\sigma^{-1}). \tag{6}
\]
The embedding is faithful even in characteristic two; a nonidentity permutation matrix never becomes the identity. The determinant is one over Z. Let J be a point stabilizer and form right-coset quotients
\[
             X=L/J,\qquad Y=L/\Gamma,\qquad
                         f:X\longrightarrow Y.            \tag{7}
\]
Right multiplication by a constant finite subgroup on L is free. Finite-group quotient theory gives smooth geometrically integral affine quotients and a finite etale cover (7) of degree [Gamma:J]=N, in every characteristic, even when the characteristic divides |Gamma|. No Reynolds averaging or tameness is needed. The dimension of both varieties is (N+1)^2-1; this is a uniform construction, not a dimension minimum.

The connected Gamma-torsor L->Y has group Gamma and remains connected after algebraic closure. Its induced cover (7) has geometric monodromy exactly the prescribed faithful action. Intermediate connected covers correspond exactly to J<=J'<=Gamma and are L/J'. In particular f is absolutely stably atomic through normal affine intermediates exactly when Gamma is primitive on its labels. For the forward direction use UAFINATOM1; for the converse any nontrivial block gives the actual finite etale factorization L/J->L/J'->L/Gamma. This is an exact equivalence in this finite receiver category, not a claimed converse for arbitrary nonfinite Keller maps.

Left Gamma multiplication commutes with every right quotient. It makes f equivariant, fixes y0=eGamma and induces the natural action Gamma/J on the whole fibre. If P is a right Gamma-torsor over any k-scheme S, twist:
\[
       X_P=P\times^\Gamma X,\quad Y_P=P\times^\Gamma Y,
                     \quad f_P^{-1}(y_P)=P/J.              \tag{8}
\]
Affineness, finite etaleness and the complete marked fibre descend effectively. Thus (8) is a canonical receiver for every framed Gamma-torsor, with no primitive generator, chosen basis or root line. Taking Gamma=S_N realizes every rank-N finite etale algebra, including disconnected ones, over arbitrary schemes.

There is also ONE FIXED universal map (7), over every field including finite fields. Indeed E=P times^Gamma L is a right SL_m torsor. Hilbert90 and the determinant exact sequence give H1(K,SL_m)=1 over every extension field K/k. A choice of a point of E identifies E/J->E/Gamma with the fixed map L/J->L/Gamma. Under this identification the canonical mark y_P becomes a K-point of Y whose complete fibre is P/J. Choosing that point of E is not canonical; the construction (8) made no such choice.

For Gamma=S_N the monodromy is S_N and the cover is strongly absolutely stably atomic in every rank N>=2, including composite N and the quartic rank. Its target has exactly N! geometric isotropy components and attains (4). Source and target are different quotient varieties; section5 proves they cannot be made the same by affine identity stabilization in characteristic zero.

The quotient-torsor and Hilbert90 construction is classical versal-torsor theory, credited in the literature note. Its use here is an exact all-rank receiver and factorization model, not a publication-priority claim.

## 4. UAALTRECV1: canonical symmetry costs half as much in high rank

For N>=3 retain L=SL_(N+1), now set H=A_N and J=A_(N-1), with J fixing label0. Choose the explicit transposition tau=(1 2), which fixes label0 and normalizes H and J. Write epsilon:S_N->C2 for permutation parity as an abstract constant-group homomorphism, including in characteristic two. Define
\[
 \alpha_\sigma(gJ)=\rho(\sigma)g\rho(\tau)^{\epsilon(\sigma)}J,
 \qquad
 \alpha_\sigma(gH)=\rho(\sigma)g\rho(\tau)^{\epsilon(\sigma)}H. \tag{9}
\]
The two right translations descend because tau normalizes J,H. Left and right multiplication commute and tau^2=1, so alpha is an S_N action.

The target mark eH is fixed: sigma*tau^epsilon is even. The entire fibre H/J has N points, identified with aJ->a(0). Formula(9) acts on it by a(0)->sigma(a(0)), since tau fixes0. Therefore
\[
                       L/A_{N-1}\longrightarrow L/A_N      \tag{10}
\]
with (9) is a canonical natural S_N receiver in EVERY characteristic. Twisting by the frame torsor of any rank-N etale algebra gives exactly that algebra as the marked complete fibre.

The geometric monodromy is A_N in its natural action. It is primitive: C3 in rank3, and 2-transitive for every N>=4 because A_(N-1) acts transitively on the remaining N-1 labels. Thus (10) is strongly absolutely stably atomic. For N=4 it supplies a canonical atomic finite quartic receiver, whereas every torus quartic receiver fails absolute atomicity.

Every twisted target is geometrically homogeneous. More precisely alpha_sigma intertwines left L multiplication with Int(rho(sigma)), so twisting L by that conjugation action gives a smooth connected group L_P acting transitively on the twisted source and target. Its marked-target stabilizer is the corresponding finite form of A_N. Therefore (10) attains the N!/2 isotropy lower bound for all N>=5, and the12 bound for an atomic quartic receiver, over every characteristic-zero base. No canonical trivialization of the group or quotient varieties is asserted.

One fixed member of this family is NOT symmetric-field-versal. Its geometric monodromy A_N is proper in S_N, contradicting UAUNIVMON1 if it were universal. More precisely the sign resolvent of this finite cover is geometrically trivial, hence a constant quadratic etale torsor over the field of definition. The canonical marked fibre identifies that constant torsor with the input algebra's sign resolvent. All other rational-point fibres of this one cover have that same sign resolvent, after every constant extension. This uses the etale sign torsor in characteristic two as well; an ordinary discriminant squareclass is not substituted there.

### Exact characteristic-zero component minima

Consider finite covers between smooth geometrically integral affine varieties, with homogeneous geometric target L/H. Minimize only the number of geometric components of H; allow different source and target. The canonical column permits effective input-wise twists as above.

| Rank | Canonical natural S_N receiver | Canonical and strongly absolutely atomic | One fixed full-degree symmetric-universal cover |
|---|---:|---:|---:|
|2|1|1|1|
|3|1|1|6|
|4|1|12|24|
|N>=5|N!/2|N!/2|N!|

Lower bounds follow from section2, with UAUNIVMON1 for the fixed column. Tori attain the1 entries. The alternating quotient attains12 and N!/2; the symmetric quotient attains6,24 and N!. These are exact minima for ISOTROPY COMPONENT COUNTS in this stated category. They are not affine-space dimensions, group dimensions or minima for an unrestricted stable-class receiver. Identity affine stabilization cannot improve them.

## 5. UAFINPI1RIG1: finite fundamental groups prohibit self-identification

Let Z be a geometrically integral normal variety whose geometric etale pi1 is finite. Every finite etale endomorphism of Z has degree one. After algebraic closure a connected finite etale cover identifies its pi1 with an index-d subgroup of the target pi1. For a self-map the two groups have the same finite order, so d=1; a finite degree-one etale cover is an isomorphism. This unstabilized statement holds in every characteristic.

In characteristic zero it also holds on Z times A^r for every r, because the geometric etale fundamental group is unchanged. Thus finite coefficient promotion on such a cylinder cannot produce a nontrivial self-cover.

For L=SL_m in characteristic zero, geometric pi1(L)=1: over C it deformation-retracts to the simply connected SU_m, and Riemann existence plus characteristic-zero base change gives the algebraic statement. Therefore for a finite subgroup F,
\[
                         \pi_1(L/F)=F.                    \tag{11}
\]
Equivalently L->L/F is a simply connected finite Galois cover. Notice that SL is essential; pi1(GL_m) is instead procyclic in characteristic zero. Replacing SL by GL in (11) would be wrong.

In the symmetric receiver (7) the source and target pi1 have orders (N-1)! and N!. In the alternating receiver (10) their orders are (N-1)!/2 and N!/2. These orders differ. They remain different under arbitrary affine identity products, so the two varieties cannot become isomorphic after such stabilization, including after extending constants. The argument also rules out a finite nontrivial self-map on either quotient cylinder. The all-rank finite receiver construction consequently cannot be promoted to the original affine-space self-map category by choosing coordinates or adding identity variables.

The same rule shows that characteristic-zero finite-pi1 examples such as the cyclic pseudo-planes necessarily use NONFINITE etale self-maps in degree greater than one. Their successful nonfinite constructions do not contradict (11), while a proposed finite completion on the same cylinder would.

## 6. What this closes and what remains

Full-degree universality gives geometric S_N and stable atomicity without a compiler-specific monodromy calculation. A different finite canonical receiver exists in every rank, with precise subgroup factorizations, effective scheme-level gluing, and sharp geometric isotropy costs. Connected-isotropy geometry cannot support high-rank canonical symmetric covers; disconnected isotropy realizes exactly the missing nonabelian monodromy. The alternating family shows that canonical marked symmetry need not be the geometric monodromy, and its frozen sign explains why a smaller canonical family is not one fixed universal map.

The original characteristic-zero whole affine-space canonical gate remains a genuinely higher-degree boundary construction, with the retained parity costs, actual whole fixed schemes, effective gluing and a separate degree-M atomicity proof. Quotient receiver pairs do not close it: their source and target have different finite fundamental groups and cannot become the same affine-space geometry by stabilization. Likewise no new unmarked dimension minimum, positive-characteristic general homogeneous no-go theorem, or projection-mixing cubic pseudo-plane action has been supplied.

The proof is written mathematics with exact finite group, matrix, quotient-fibre and descent-action controls in verify.py. Those controls do not formalize etale descent, all-rank group theory, topology, homogeneous-space pi1 or the imported atomicity lemmas. Independent full-proof replay, formal verification and external review remain pending. Classical constructions and inputs are explicitly credited; no novelty priority is certified.
