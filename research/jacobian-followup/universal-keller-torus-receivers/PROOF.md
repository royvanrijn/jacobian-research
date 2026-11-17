# Canonical torus receivers and stable atomicity

14 November 2026. The [cyclic acyclic-surface receiver](../universal-keller-acyclic-boundary/PROOF.md) does not immediately extend to symmetric groups. Tori give a different complete classification: natural S2, S3 and S4 receivers exist, with exact dimensions1,1,2, whereas S_N for N>=5 is excluded in every torus dimension. Fixed torus universality is stricter than canonical descent. We also prove that affine identity stabilization cannot create a new factorization of a finite etale torus endomorphism in characteristic zero.

The receiver here allows finite etale maps between torus torsors, with a marked target section and the complete specified fibre. Source and target need not be the same form of a torus unless explicitly stated. This is not a whole affine-space Keller receiver. All finite-factor statements require dominant factors and normal affine intermediate varieties. Geometric degree is denoted N.

## 1. UAFINATOM1: finite maps force finite factors

Let f:X->Y be finite dominant between integral affine varieties and suppose f=h g with g:X->Z and h:Z->Y dominant regular morphisms, Z affine integral. Then both factors are finite. On coordinate rings write A<=B<=C, with C finite over the Noetherian ring A. The subring B is an A-submodule of C and hence finite over A. Every element of C is integral over A, hence over B; finite generation makes C finite over B.

If f is finite etale and X,Y,Z are normal, both factors are finite etale. Indeed B is the normalization of A in its fraction field, an intermediate extension of the finite etale extension. Strict-etale localization splits the latter as copies of the normal local base; normal intermediate algebras are products of the same base. Thus B/A is etale. A morphism between finite etale A-schemes is etale, giving C/B.

When X and Y are geometrically integral and geometrically normal, primitive geometric monodromy makes f absolutely stably atomic among dominant factorizations through normal affine varieties. Every identity stabilization remains finite etale, and both generic factors would give a proper intermediate field. A degree-one finite factor is an isomorphism by normality. Prime geometric degree is a special case, with no acyclicity or trivial-canonical-bundle assumption. This extends the prime/primitive criterion to different normal affine intermediate varieties in the finite category; it does not apply to the nonfinite pseudo-plane maps.

Normality cannot be dropped. In characteristic zero, the square map on Gm factors as
\[
 s\longmapsto(t,b)=(s^2,s(s^2-1))\longmapsto t,
 \quad Z=\{b^2=t(t-1)^2,\ t\ne0\}.                         \tag{1}
\]
The first factor is finite birational but not an isomorphism: it identifies s=1 and s=-1 at the node(1,0). The second has degree two. Their composite is finite etale of prime degree, so arbitrary nonnormal intermediates destroy the asserted atomicity.

## 2. UATORUSSTABLE1: affine stabilization preserves torus factorization

Let T be a torus over a characteristic-zero field. Every finite etale endomorphism G of T times A^r has the form
\[
              G=(\phi\times\operatorname{id})\circ H,        \tag{2}
\]
where phi:T->T is an isogeny followed by translation, and H is an automorphism preserving projection to T. In particular deg G=deg phi.

Over a separable closure, units in the coordinate ring of T times A^r are constants times characters of T. Thus its torus output is independent of the affine coordinates and descends to phi. Dominance makes its character matrix full rank, hence phi a finite etale torus map. If G(t,x)=(phi(t),g(t,x)), take H(t,x)=(t,g(t,x)). It is finite by the preceding ring argument and etale because both G and phi times id are etale. On every geometric torus fibre it is a finite etale endomorphism of A^r in characteristic zero, hence an automorphism by simple connectedness of affine space. Its degree is one, so H is an isomorphism. Global finiteness is essential; this does not prove that arbitrary nonproper affine Keller maps are automorphisms.

Let M=X^*(T), with its actual Galois action, and let A be the injective finite-index endomorphism of M induced by a torus isogeny. Then
\[
 \boxed{\text{the isogeny is stably atomic over the base field}
 \iff A\ne BC\text{ with }B,C\in\operatorname{End}_{\rm Gal}(M),
                   \ |\det B|,|\det C|>1.}                 \tag{3}
\]
A stabilized regular factorization has finite etale factors by UAFINATOM1; (2) makes their degrees equal to their character-lattice indices. Hence it induces the prohibited integral factorization. Conversely an integral factorization gives a torus factorization already before stabilization. Translation factors are automorphisms and cause no change.

After algebraic closure Smith normal form shows that a nontrivial isogeny is absolutely stably atomic exactly when its degree is prime. This last statement holds in every characteristic in which the isogeny is etale. Over the original field, composite-degree stable atoms can exist.

### A composite stable atom that is not absolute

Take a nonsplit separable quadratic extension E/k in characteristic zero and T=Res_(E/k) Gm. Galois swaps a basis of M=Z2. Its equivariant endomorphisms are
\[
 \begin{pmatrix}a&b\\b&a\end{pmatrix},\qquad
                   \det=a^2-b^2.                            \tag{4}
\]
This determinant is never plus or minus2, since the two squares modulo4 have difference0,1 or3. Multiplication by2 on T has degree4, so any nontrivial factorization would need two degree-two factors. Formula(3) therefore makes it atomic after every affine identity stabilization over k. Over E it becomes (x,y)->(x2,y2) and decomposes into two quadratic coordinate maps, so it is not absolutely atomic.

It also factors over k through a different torus: the stable index-two lattice generated by(1,1),(1,-1) has diagonal Galois action, and corresponds to Gm times the norm-one torus of E/k. This explicitly separates the same-variety self-map category from arbitrary normal intermediate tori.

### Unit lattices cannot be repaired by affine stabilization

If T1 times A^r is isomorphic to T2 times A^s, geometric units modulo constants give a Galois-equivariant isomorphism X^*(T1)=X^*(T2). Thus the tori are already isomorphic; their ranks agree and then r=s. More intrinsically, over a separable closure the subalgebra generated by all units is exactly the coordinate ring of the torus factor, so an isomorphism of cylinders restricts to a torus isomorphism and descends. An isogeny between nonisomorphic torus forms is not silently converted into a self-map by identity stabilization. The unit-lattice observation works in every characteristic; the structural factor theorem(2) uses characteristic zero.

## 3. UATORUSCLASS1: a complete permutation criterion for torus receivers

Let G act faithfully and transitively on a finite set Omega of size N>1, over a field of characteristic p>=0. A canonical marked realization by an etale map between torus torsors exists in some dimension if and only if
\[
 \boxed{p\nmid N\quad\text{and}\quad
 \text{there is an abelian regular subgroup }A\le\operatorname{Sym}(\Omega)
                  \text{ normalized by }G.}                  \tag{5}
\]
For characteristic zero the first condition is automatic. Regular means free and transitive; A need not be contained in G. A construction has dimension at most N-1.

Every regular map between split tori is monomials times constants, by the description of units. Etaleness makes its integer exponent matrix invertible in the field. It is therefore finite etale, of degree the absolute determinant, prime to the characteristic. After translations it is an isogeny, whose geometric deck group is an abelian regular group on every fibre. The compatible label action normalizes this group. Thus all marked fibre ranks equal generic degree; tori have no smaller boundary-fibre escape. These facts prove necessity in(5).

For sufficiency, identify Omega with A using one point once in the construction, so that G embeds into Hol(A)=A semidirect Aut(A). Put m=exponent(A), let Delta be the Galois group of k(mu_m)/k, and let
\[
 B=X^*(A)=\operatorname{Hom}(A,k_s^*),\qquad
 M=\mathbb Z[B\setminus\{0\}],\qquad
 L=\ker(M\longrightarrow B).                               \tag{6}
\]
The map sends the basis element indexed by b to b. Aut(A) and Delta act on B with the usual commuting contragredient and cyclotomic actions, so the finite-index inclusion L->M is equivariant. Duality for groups of multiplicative type gives a finite etale isogeny
\[
          1\longrightarrow A\longrightarrow T_M
                 \mathop{\longrightarrow}^{\Phi}T_L
                    \longrightarrow1.                       \tag{7}
\]
Its kernel is the constant group A, not a mistakenly untwisted mu_m: the cyclotomic action is already included in B. The two tori have rank N-1. Aut(A) acts equivariantly; Hol(A) acts on the source by kernel translations and these automorphisms, and on the target through Aut(A). The target identity has exactly the required regular A fibre. Restrict to G and twist by an arbitrary G torsor P. Effective etale descent gives a map of torus torsors with marked fibre P times^G Omega, exactly the input rank-N algebra.

M is a permutation lattice for the combined action. After any fieldwise G twist the source's linear torus remains quasi-trivial, so Hilbert90 trivializes its torsor over the field. This is not a canonical choice of source origin. Source and target tori can also have different character lattices; the stabilization rule above prevents assuming them isomorphic.

Criterion(5) is the abelian-type regular-subgroup condition in Greither--Pareigis Hopf--Galois theory, Theorem3.1. That published classification is credited; the torus receiver and descent above are the present application. Do not strengthen it to requiring A contained in G. For example the regular S3 action on six points is generated inside Hol(C6) by x->x+2 and x->1-x, and normalizes C6 without containing that regular cyclic subgroup.

## 4. UATORUSSYM1: exact symmetric torus classification

For the natural S_N action, A in(5) is an abelian transitive normal subgroup of S_N. The possibilities are exactly C2 at N=2, C3 at N=3 and V4 at N=4. There are none for N>=5. Consequently:

| Rank | Allowed characteristics | Minimum torus dimension | Absolute stable atomicity |
|---|---|---:|---|
|2|not2|1|yes|
|3|not3|1|yes|
|4|not2|2|no, in every torus dimension|
|N>=5|none|no realization|not applicable|

The exclusions survive arbitrary torus dimension and every affine identity stabilization of the fixed cover: generic abelian monodromy and full fibre rank are unchanged, and the finite-chart normalizer argument still applies. No larger generic-degree torus boundary construction exists, since the maps are finite.

### Quadratics

Use squaring on Gm, source C2 translation by -1, target action trivial and marked target1. Twisting gives the full quadratic torsor as fibre. Degree two gives strong finite atomicity by UAFINATOM1.

### Cubics over every field of characteristic different from three

Put E0=k[theta]/(theta2+theta+1), a quadratic etale algebra, and let T0 be its norm-one torus. In coordinates a+b theta,
\[
 T_0:\ a^2-ab+b^2=1,\qquad
 [3]_{T_0}(a,b)=(a^3-3ab^2+b^3,\ 3ab(a-b)).                         \tag{8}
\]
The element theta=(0,1) has norm1 and order3. Translation tau(a,b)=(-b,a-b) and inversion iota(a,b)=(a-b,-b) generate S3. On the target tau acts trivially and iota by inversion. Cubing is equivariant and its entire fibre at(1,0) is
\[
                    (1,0),\quad(0,1),\quad(-1,-1).          \tag{9}
\]
These are rational over every such base field, with the natural S3 action. Thus twisting by the S3 torsor of frames of any cubic algebra produces a canonical marked finite etale degree-three map. This is defined over Q in characteristic zero, without requiring a primitive cube root in the field. The norm-one model also works in characteristic two. The relative logarithmic derivative is multiplication by3, proving etaleness.

The target torus is the sign twist of T0 by the cubic discriminant torsor. In characteristic not2 it is the norm-one torus of the quadratic algebra with square class -3 times the discriminant. The source is a torsor under that same torus, killed by3. Hilbert90 identifies H1 of any one-dimensional norm-one torus with K*/norm(E*), a group killed by2. Thus the source torsor is trivial over every field, though no trivialization is selected by the construction. The source S3 action has no fixed point, since tau is a nontrivial torus translation; vector-representation identity stabilization cannot create a global equivariant origin.

### Quartics over every field of characteristic different from two

Let
\[
          T=\{(x,y,z)\in\mathbb G_m^3:xyz=1\}.               \tag{10}
\]
Coordinatewise squaring has degree4. Its four kernel points are the sign triples with product1, forming V4. Translation by V4 and coordinate permutations S3 give V4 semidirect S3=S4 on the source. The target kills V4 and permutes the coordinates; its identity is marked. The full fibre is the natural four-point S4 set. Twist by the frame torsor of any quartic algebra.

The target is the norm-one torus of its cubic resolvent algebra R. The source torsor is killed by2, whereas H1(K,Res^1_(R/K)Gm)=K*/norm(R*) is killed by3, since scalar norms are cubes. It is therefore trivial over every field. Again no equivariant origin is chosen or available. Over an algebraic closure squaring in two torus coordinates decomposes into the two quadratic coordinate isogenies, so no such quartic torus map is absolutely atomic.

Dimension one suffices for ranks2,3 and dimension zero cannot. A one-dimensional torus has only cyclic finite etale subgroups, so cannot have the required V4 kernel in rank4; dimension two is minimal. These are dimensions in the torus receiver, not unmarked affine-space Keller minima.

## 5. UATORUSFIX1: canonical families do not give fixed symmetric universality

A fixed torus map cannot realize every cubic or quartic algebra over every extension field. All its geometric fibres have its generic degree, so using a larger-degree map does not help.

A degree-three torus cover is a torsor under one fixed finite etale group scheme A of order3. Its fibres' Galois actions are affine transformations of A. In characteristic two, discriminant here means the etale quadratic sign-resolvent torsor, not the ordinary polynomial discriminant square class. Translations are even three-cycles; hence the discriminant of every cubic fibre is the fixed C2 character describing the Galois action on A. It cannot vary with the target.

For a quartic cover the geometric kernel is either C4 or V4. In the cyclic case all fibre Galois groups lie in Hol(C4), a group of order8, so a fibre of type K times a connected separable cubic field is impossible. In the V4 case its fixed three nonzero kernel elements give the fixed cubic resolvent: translations act trivially on the three pair partitions, while the linear Galois action is fixed. It cannot vary with the target.

The every-extension quantifier makes these unconditional obstructions even over a small base. Over K=k(t), take a quadratic etale algebra D_t ramified in t and not descended from the constant field: use sqrt(t) in characteristic not2, or X2+X+t in characteristic two. Cubic K times D_t has variable discriminant. Quartic K2 times D_t has cubic resolvent K times D_t. In the cyclic quartic-kernel case use a separable cubic extension over k(t): X3-t when characteristic is not3, and X3-X-t in characteristic3. This produces the forbidden fibre types without an inverse-Galois assumption.

For N>=5 a fixed torus cover cannot even have a fibre with full natural S_N Galois action, by(5). Such algebras exist over an extension of every field: take independent root variables and the field of symmetric rational functions. Thus one fixed field-versal natural symmetric torus cover exists exactly for rank2 in characteristic not2, where x->x2 on Gm works. Canonical torus families nevertheless exist also for ranks3 and4. Fieldwise triviality of the source torsor does not turn the input-dependent target torus into one fixed torus, and unit lattices prevent an identity-stabilization repair.

In characteristic zero the fixed-cover obstruction also survives arbitrary finite coefficient promotion on a cylinder T times A^r. By UATORUSSTABLE1 every finite etale endomorphism there is a torus isogeny times identity after a source automorphism. Its full fibre algebras are therefore still those of one fixed torus isogeny, so it is not symmetric-universal in any rank N>=3. A viable different construction must drop global finiteness or change this source geometry.

There are fixed framed universal covers for abelian groups. Taking G=A in(6)--(7), T_M is quasi-trivial and H1(K,T_M)=0 over every extension. The connecting map T_L(K)->H1(K,A) is therefore onto. This is a fixed finite etale cover between tori, not automatically a self-map of one torus form. For C3 there is a one-dimensional fixed self-map over every field of characteristic not3: [3]:T0->T0 from(8). Its kernel is constant C3, and H1(K,T0) has exponent2, so all C3 torsors occur. This does not realize cubic algebras with varying S3 discriminant.

For quartic torus covers with V4 kernel, arithmetic monodromy is V4 semidirect Gamma, where Gamma is the Galois action on the cubic resolvent. If the resolvent is connected, Gamma contains a three-cycle; the arithmetic group A4 or S4 is primitive. Such a map is atomic over that field in the strong normal-affine category but decomposes after extending constants. This is an explicit canonical comparison showing why arithmetic primitivity cannot replace geometric primitivity for absolute atomicity.

For the constructed minimal rank-two quartic family, this gives a complete characteristic-zero distinction after choosing any fieldwise source origin (which only changes the map by translation):

| Cubic resolvent R | Stable atomicity among self-maps over k | Atomicity through arbitrary normal affine intermediates over k | Absolute stable atomicity |
|---|---|---|---|
| k times k times k | no | no | no |
| k times a nonsplit quadratic field | yes | no | no |
| a cubic field | yes | yes | no |

The middle row is exactly the restriction-of-scalars example from Section2: the norm-one torus of k times E is Res_(E/k)Gm. In the first row both quadratic coordinate factors already exist over k. In the last row arithmetic primitivity gives the stronger atomicity statement. Thus the table is proved without classifying all imprimitive maps or presuming that a generic intermediate field is the same ambient variety. It concerns this minimal quartic family, not every higher-dimensional quartic torus cover.

## 6. UAPSEUDOSYM1: the original pseudo-plane cannot carry full symmetric symmetry

In characteristic zero every finite subgroup of Aut(S_k), k>=2, is cyclic. Dubouloz--Palka equation(5.12) gives every automorphism in the form
\[
          u\longmapsto\lambda^k u,\qquad
          w\longmapsto\lambda w+u^2Q(u).                    \tag{11}
\]
The multiplier lambda is a group character. Its kernel is the additive group of polynomials Q, which has no nontrivial finite subgroup in characteristic zero. Thus finite subgroups inject into the multiplicative group and are cyclic. The argument is geometric; finite action data over any characteristic-zero field spread to C, where the cited classification applies.

It follows that no etale map on this same unstabilized surface can have a complete natural marked S_N fibre for N>=3 with actual compatible S_N actions: the faithful label action would force a faithful noncyclic finite source subgroup. This is independent of the chosen shear or generic degree.

Nor can identity stabilization repair this while retaining the original surface projection. If a finite source action on S_k times A^r descends to an action on S_k, its label action on a fibre of F times id factors through that cyclic image: all points in the fibre have the same affine coordinate. A faithful natural S_N label action is therefore impossible for N>=3, allowing arbitrary fibrewise polynomial actions in the added coordinates. A successful cubic stabilization must genuinely mix the surface projection.

A stronger monodromy rule closes every stabilized cyclic-degree-N route for N>=4, even with arbitrary mixing. The full-degree finite-chart normalizer argument underlying UACANATOM1 applies to any connected normal source and target, not just affine spaces: a complete N-point fibre at a fixed target exhausts normalization degree, supplies a group-invariant finite etale chart, and the label action normalizes generic monodromy. A cyclic regular C_N subgroup is not normal in natural S_N for N>=4. Since identity stabilization and left-right equivalence preserve geometric monodromy and degree, no representative of the retained eta_(N,Q) class can acquire such a full natural S_N lift. Rank3 is the only nonquadratic cyclic case that passes this necessary monodromy test, and its projection-mixing stabilization remains open.

The new torus constructions change the ambient geometry. They do not supply the missing S3 action on S3 or solve general whole affine-space canonical descent. For N>=5 a surviving different receiver needs nonabelian geometric monodromy in full marked degree, or a genuinely smaller boundary fibre in a higher-degree map outside the finite torus category.

## 7. Evidence and scope

The six claims have written structural proofs: finite factor closure, unit lattices and finite torus-cylinder structure, the equivariant lattice construction, explicit low-rank torsor twisting, frozen resolvents, and the pseudo-plane automorphism/normalizer gates. Exact controls certify the cubic map and all its group relations over integer polynomial rings, complete marked fibres, the quartic natural S4 action, lattice indices and invariant kernels, arithmetic versus geometric permutation actions, and the singular-intermediate countercontrol. They do not formalize Galois descent, torus duality, simple connectedness or Hopf--Galois theory. Published inputs are credited; no independent review, formal verification or novelty priority is asserted.

The original whole affine-space objective remains open in its surviving higher-degree boundary and other receiver categories. The torus category itself is classified, and the old unstabilized pseudo-plane symmetric-extension route is closed. Nothing here changes the retained unmarked affine-space dimension bounds or characteristic-two quartic atomicity interval.
