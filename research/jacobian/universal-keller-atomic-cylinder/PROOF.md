# Atomicity categories and reduction of rigid cylinders to affine Keller data

15 November2026. Written all-degree arguments and narrow exact controls. No independent whole-proof replay, formal verification, external review or publication-priority certificate. OP-KDESC remains open. Classical inputs and retained statements are pinned in [sources.json](sources.json).

There are two distinct conclusions. A nonfinite map is never atomic through all normal affine intermediate varieties. Among self-map factorizations, however, atomic nonfinite maps do exist, and on every geometric torus-quotient cylinder their base map must be an automorphism. Canonical marked data then restrict to a genuine affine fibre, without the earlier metabelian automorphism hypothesis. A dilation construction also gives infinitely many inequivalent atomic nonconstant families on one fixed cylinder; it supplies no new canonical symmetric action.

## Conventions and the classical self-map lemma

A variety is separated and of finite type over a field. All maps considered in factorizations are dominant regular morphisms. Self-atomic means that a factorization into self-maps of the same variety has an automorphism factor. Normal-affine-atomic allows a different normal integral affine intermediate variety and requires an isomorphism factor. Absolute stable versions require the assertion after every field extension and every affine identity stabilization. These categories must be named; they give different answers for nonfinite maps.

A birational etale endomorphism of a normal variety is an automorphism, in every characteristic. Zariski Main makes it an open immersion. Ax's injective-endomorphism theorem makes this self-open-immersion surjective, hence an isomorphism. This is a classical corollary, not a new theorem of Ax type. No acyclicity, affineness, Picard condition or compactification is required. Ax's printed page1 theorem and Corollary1 were read and visually inspected. Between DIFFERENT varieties a birational etale map need not be an isomorphism; UAACYCLICRIG1 retains its separate cohomological role there.

## 1. UANORMAFF1: the exact normal-affine category boundary

Let f:X->Y be a separated dominant etale map of generic degree M>1 between geometrically integral geometrically normal affine varieties. Over an algebraic closure, f is normal-affine-atomic if and only if it is finite and its geometric monodromy is primitive. The same condition characterizes absolute stable normal-affine atomicity over the original field.

If f is nonfinite, let C be the normalization of Y in k(X). Varieties over a field have finite normalization in a finite field extension, so C->Y is finite and C is normal affine. Zariski Main identifies X with an open subset of C. Thus

    X --j--> C --pi--> Y                                    (1.1)

is a factorization through a permitted normal affine intermediate. The first factor is not an isomorphism, since otherwise f would be finite. The second has degree M>1, so is not an isomorphism either. This excludes atomicity. The finite completion may be ramified or singular: the assertion allows dominant REGULAR factors and does not silently restrict both factors to be etale.

If f is finite etale, UAFINATOM1 proves that all normal affine factors are finite etale. A proper block system gives a nontrivial intermediate finite etale cover, while primitive monodromy gives no proper intermediate field. These are the two directions. Geometric monodromy is preserved by extending algebraically closed constants and by adjoining independent identity variables, giving the absolute stable formulation. Conversely absolute atomicity already tests the algebraic closure.

Combining this category boundary with UAINFRASELF1 and UAPPDEGREE1 gives an exact extension of the degree spectrum. Even if finiteness is omitted at the outset and arbitrary affine identity cylinders are allowed, the nontrivial degrees of absolute stable normal-affine-atomic etale endomorphisms of smooth affine geometric torus-quotient cylinders are exactly the prime powers. Finiteness follows from(1.1); the intrinsic-lattice obstruction survives cylinders; existence already holds without adding an affine cylinder. This is an existence spectrum as the source varies. It is not a degree classification for self-atomic nonfinite maps.

## 2. UASELFATOM1: remove acyclicity from the primitive self-map criterion

Let X be smooth geometrically integral, with geometrically constant units and geometrically torsion canonical bundle. An etale endomorphism with primitive geometric monodromy is absolutely stably self-atomic. No acyclicity or affineness hypothesis is needed, and the argument works in every characteristic.

Work over algebraically closed constants and trivialize K_X^m, with generator omega. For a dominant regular endomorphism a, the m-th power of its determinant defines a regular function J_a by a^*omega=J_a omega. For f=b a the chain rule says

    J_f=(J_b composed with a) J_a.                          (2.1)

Since f is etale, J_f is a unit. Both factors on the right are units. Constant units and dominance force both J_a and J_b to be nonzero constants. Thus both regular factors are etale everywhere. Primitivity makes one generic field extension have degree one. The classical self-map lemma above makes that factor an automorphism. This replaces the acyclic birational-rigidity step used in UAACYCLICATOM1; all its retained applications remain valid.

Geometrically constant units, torsion canonical bundle, smoothness and geometric integrality persist under affine identity stabilization. The generic extension and its primitive geometric monodromy persist under that stabilization and extension of constants. Apply the same argument there. In particular geometric prime degree suffices, but the theorem is not restricted to prime degree. No converse for imprimitive monodromy is asserted, and the constant-unit hypothesis is not dropped without a replacement.

## 3. UACYLDICH1: atomic maps cannot spend degree in both base and fibre

Let k have characteristic zero, let Z be a smooth geometrically integral affine geometric torus quotient, and put X=Z times A^r. By UA1RIGIDCYL1 every regular endomorphism projects to Z. For an etale endomorphism write

    F(z,v)=(b(z),g(z,v))=(b times id) H(z,v),
    H(z,v)=(z,g(z,v)).                                      (3.1)

The block differential makes b and H etale; UATORUSDOM1 makes dominant b finite etale. Put D=deg(b), h=deg(H) and M=Dh. Then

    H has degree1 iff H is an automorphism iff F is finite.  (3.2)

The first equivalence is the classical self-map lemma. If H is an automorphism, F is finite. Conversely a finite F makes H finite by UAFINATOM1, and each geometric affine fibre is a connected finite etale cover of A^r. Characteristic-zero affine-space simple connectivity makes those vertical maps automorphisms, giving h=1, as already proved in UA1RIGIDCYL1.

If M>1 and F is self-atomic, exactly one of the following holds:

- h=1: F is a finite base self-cover followed by a source automorphism;
- D=1: b is an automorphism, and all nontrivial degree is in the vertical Keller family.

Indeed D,h>1 would make(3.1) an explicit nontrivial factorization into self-maps. In particular every NONFINITE self-atomic F has D=1. This uses no assumption on Aut(Z), on the label group, or on the base degree being prime. It applies to every absolute stable self-atom and to every identity stabilization. It is a necessary dichotomy, not a converse: D=1 alone does not establish atomicity.

### A sufficient criterion for the vertical case

Suppose b is an automorphism and the geometric monodromy of the generic vertical Keller map over an algebraic closure of k(Z) is primitive. Then F is absolutely stably self-atomic. Here the base parameters are held constant for that relative geometric monodromy test. This condition is sufficient, and is stronger than just checking a few closed fibres.

Consider ANY regular self-factorization F=Q P. Both factors project to dominant base endomorphisms c,d. Since d c=b, UATORUSDOM1 and multiplication of finite degrees force c,d to be automorphisms. Restricting to a geometric affine fibre, the composition is a Keller map between A^r's. The usual polynomial determinant product is a nonzero constant; both factors' determinants are units in a polynomial ring, hence nonzero constants. Thus P and Q are etale on every geometric fibre and, with their etale base maps, everywhere. After using c to identify the intermediate generic fibre, their generic vertical field degrees multiply to h. Relative geometric primitivity forces one to have degree one. Its total base degree is also one, so the self-map lemma makes it an automorphism. This argument does not assume that a pre-existing factor is etale.

The same proof works after every constant extension and identity stabilization: the relative primitive group persists and the base remains the same geometric torus quotient. This handles the nonconstant units of X by its rigid base, rather than applying Section2 with a false constant-unit hypothesis. No necessity of relative primitivity is claimed. A special vertical degree can drop below h, so one closed fibre cannot certify this condition or global atomicity.

## 4. UACYLCAN1: remove the metabelian hypothesis from atomic canonical reduction

Retain Z,X and characteristic zero. Let F be NONFINITE and self-atomic, with compatible polynomial actions of a finite label group G on source and target, a G-fixed target point y, and complete transitive marked fibre Omega of N>1 points. By Section3 the base b is an automorphism. Let y0 be the projection of y and put z0=b^-1(y0). Equivariance shows that y0 and z0 are fixed by their respective projected G actions. Every point of Omega has the same projection z0. Therefore restriction gives

    g_z0:A^r -> A^r                                         (4.1)

with actual polynomial G actions, fixed target and EXACTLY the original complete N-point fibre. Its Jacobian is a nonzero constant, so a linear target normalization gives a genuine polynomial Keller datum. Its generic degree h0 satisfies

    N <= h0 <= M.                                          (4.2)

The upper bound is the strict-henselian point budget from UAFIXDEFECT1, applied at a geometric generic point of this target affine fibre; it is not an assertion that special generic degrees remain constant.

Thus the conclusion of UAVERTRED1 now holds on EVERY geometric torus-quotient base for nonfinite self-atomic data. No metabelian Aut(Z) assumption is needed. For arbitrary transitive labels UASYLOWDEF1 gives

    M >= h0 >= N+2+indicator(2 divides N).                   (4.3)

For natural S_N labels the ordinary arbitrary-action affine bound UANONLINEARCOST1 applies to r itself:

    r>=N+1   if N is a prime power;
    r>=N     if N is a sum of two positive powers of one prime
                 and is not already a prime power;
    r>=N-1   otherwise.                                    (4.4)

The total dimension dim(Z)+r cannot pay this affine-fibre cost. For natural S_N with N>=5, the finite alternative is already excluded by UAINFRASELF1, so nonfiniteness need not be separately assumed for a self-atomic candidate. For example every atomic canonical quartic attempt on ANY nonfinite torus-quotient cylinder requires r>=5 and M>=7; the old W example is no longer the only base controlled. No attainment is claimed.

Restriction preserves the marked symmetry, not an atomicity certificate on(4.1): global atomicity does NOT imply atomicity of every special vertical fibre. Nor does this claim a functorial coordinate choice on all input twists. These points distinguish a genuine reduction from a construction.

### Relative restriction glues, but does not choose affine coordinates

The restriction also makes sense without choices in families. Let B be a reduced finite-type scheme over a characteristic-zero field, and let an etale-locally trivial family with fibre Z times A^r carry a relative etale self-map, compatible G actions, a fixed target section and a complete finite etale marked fibre. Assume every geometric fibre map is nonfinite and self-atomic. All transition automorphisms project to Z by A1-rigidity, so there is an intrinsic Z-bundle. The induced relative base map is etale by the block differential and is an isomorphism on every geometric fibre by Section3. It is therefore an isomorphism. Pulling back the projected marked section through this isomorphism gives the unique source base section z0.

Restrict source and target to their two sections. They are etale-locally A^r bundles, and their relative etale map, G actions and complete marked fibre glue automatically and commute with base change. Source and target bundles need not be naturally identified or have chosen affine coordinates. This gives a receiver PAIR of affine-space forms, not a fixed polynomial self-map of A^r over B. It does not settle triviality of those forms, effective affine-space landing or the original canonical self-receiver. It shows precisely which part of a proposed cylindrical family already belongs to the remaining affine problem.

## 5. UACYLDEF1: general nonconstant atomic families and a stable invariant

Let k have characteristic zero and let f:A^r->A^r be a polynomial Keller map of degree M>1 with primitive geometric monodromy. Normalize f(0)=0 and Df(0)=id by affine target changes, and write its homogeneous components f=v+sum_(j>=2) f_j. For any regular parameter a define the polynomial family

    f_a(v)=v+sum_(j>=2) a^(j-1) f_j(v).                     (5.1)

Where a is nonzero this is a^-1 f(a v); where a=0 it is the identity. Its determinant is identically1, by the identity on the dense open a!=0 and polynomial continuation. This is the standard dilation deformation; no new deformation technique is claimed.

On the FIXED smooth affine variety Gm_t times A^r put

    a_m(t)=product_(j=1,...,m)(t-j),
    F_m(t,v)=(t,f_(a_m(t))(v)),  m>=1.                     (5.2)

The generic vertical map is conjugate to f and has the same primitive relative geometric monodromy. Section3 makes every F_m absolutely stably self-atomic. It preserves

    Omega=dlog(t) wedge dv1 wedge ... wedge dvr.            (5.3)

It has generic degree M, but the vertical map has degree1 exactly at the m distinct nonzero base points t=1,...,m and degree M elsewhere. Consequently F_m is nonfinite. Section1 shows why it is not normal-affine-atomic despite its self-atomicity.

The maps F_m are pairwise inequivalent under arbitrary regular source and target automorphisms of this same cylinder, after every common affine identity stabilization. Indeed all such automorphisms project to automorphisms of Gm, t->c*t or c/t. A left-right equivalence between two maps with identity base forces inverse base changes. Vertical automorphisms preserve each fibre degree. Thus the number m of base points where the vertical degree drops is invariant. Identity stabilization preserves both this degree function and the A1-rigid base. The proof does not assume that the left-right automorphisms were initially declared projection-preserving.

The retained UAM1 determinant-one cubic supplies a concrete seed with M=3. The checker reuses its source_map function, normalizes its first derivative, verifies the full parameter-dependent Jacobian, and checks a collision and degree-drop controls. The written argument proves the assertion for EVERY primitive Keller seed in any dimension and generic degree, without a degree search.

### Exact degree spectra on one fixed fourfold

The all-degree statement can be made on a single geometry X=Gm times A^3. The retained UAALLATOM1 provides a determinant-one three-variable primitive seed in EVERY degree N>=3: take the cubic compiler for N=3, and for N>=4 fix a_N=1 and every other promoted coefficient to zero. For N>=4 its complete generic inverse polynomial is

    S+pi*S^3+pi^N*S^N+b*S^2-c/2,                           (5.4)

of degree N, with geometric S_N by UAQTILT1; in degree3 omit the extra pi^N*S^N term. No coefficient is promoted as an extra source coordinate. Normalize the first derivative and apply(5.1)-(5.2). Thus for every N>=3 there are infinitely many stably inequivalent, volume-preserving, absolutely stably self-atomic NONFINITE endomorphisms of this SAME fourfold, defined over Q.

There is no nonfinite etale endomorphism of degree1 or2 on X. Degree1 is the classical self-map lemma. For degree2 the base/vertical degrees are either (2,1), which is finite by Section3, or (1,2). The latter would give a quadratic polynomial Keller self-map of A^3 over the algebraic closure of k(t), impossible by the Galois/quadratic case retained in UAACYCLICGAL1. Hence the exact nonfinite etale degree spectrum of X is {3,4,5,...}, and every degree in it has the atomic volume-preserving representatives just constructed.

For completeness every finite etale endomorphism of X is a source automorphism followed by (c*t^d,v), with d a nonzero integer. Its geometric monodromy is the regular cyclic group C_|d|. It is normal-affine-atomic exactly when |d| is prime. It is self-atomic exactly in the same prime degrees: a composite |d| factors into actual power self-maps, while prime degree has primitive finite monodromy. A finite volume-preserving prime-degree model is (t,v1,v2,v3)->(t^p,v1/p,v2,v3).

| Maps on the fixed X=Gm times A^3 | Exact nontrivial atomic degrees |
|---|---|
| Finite maps, either named atomicity category | All primes |
| Nonfinite self-maps, absolute stable self-atomicity | All integers at least3 |
| Nonfinite maps, normal-affine atomicity | None |

In particular absolute stable self-atoms with a global preserved volume exist in EVERY degree N>=2 on this one fourfold. This differs from UAPPDEGREE1's strong normal-affine spectrum on varying torus-quotient sources. The entries are existence spectra; they do not claim that every map of an allowed degree is atomic or universal.

The construction proves that nonfinite coefficient variation can give new atomic cylinder maps, whereas finite coefficient variation is removed by UA1RIGIDCYL1. It cannot supply a canonical S_N action absent from the seed. Its exceptional identity fibres have rank one, so they do not create arbitrary smaller transitive marked ranks. The Gm unit also persists after affine identity stabilization, excluding affine-space identification.

## 6. What the extensions do and do not close

The changes close five general gates: the exact normal-affine atomicity category boundary; the unnecessary acyclicity hypothesis for the constant-unit pluricanonical self-map argument; the base/vertical atomic degree dichotomy; the metabelian hypothesis for nonfinite atomic canonical cylinder reduction; and the existence of nonconstant stably inequivalent atomic dilation families at fixed degree and source.

The core higher-degree canonical affine-space problem is still open. A new torus quotient cannot solve it by moving degree into both base and fibre while retaining self-atomicity. Effective simultaneous fixed-locus landing, compatible marked actions, source/target identification over torsor bases, and ambient atomicity for smaller marked rank still require proof. Self-atomicity must not be upgraded to normal-affine atomicity for a nonfinite map. The new families do not solve JC(2), HC(4), dimension minima or prescribed monodromy realization. Classical corollaries and repository deductions are distinguished; no publication-priority conclusion follows from the literature search.
