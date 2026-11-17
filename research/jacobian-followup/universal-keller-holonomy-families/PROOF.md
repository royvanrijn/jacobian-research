# Holonomy families: all-degree atomicity, torsor universality and canonical obstructions

15 November 2026. Written proofs with small exact controls, not formal verification or independent review. Classical torus-cover, norm-cohomology and Hopf--Galois inputs are credited below and pinned in `sources.json`. Publication priority has not been established. OP-KDESC remains open.

## 1. Conventions and the family

Unless stated otherwise work geometrically in characteristic zero. A marked transitive datum means a dominant etale map, compatible source and target actions of a finite group G, a fixed target point, and its COMPLETE finite reduced fibre Omega with a specified transitive G action. For natural symmetric data the action is the faithful natural S_N action. No equality between marked rank and generic degree is presumed for a nonfinite map.

Let A={(v1,v2,v3) in Gm^3: v1*v2*v3=1}. Let T0 be the cyclotomic quadratic norm-one torus from the retained UAQUARTSELF1 construction, with its rational point theta of order3. For m>=1 put

    T_m=T0_x x Gm_z x A_1 x ... x A_m,
    Z_m=T_m/F,  F=S3=<r,s>,
    r(x,z;a_1,...,a_m)=(theta*x,z; r*a_1,...,r*a_m),
    s(x,z;a_1,...,a_m)=(x^-1,-z; s*a_1,...,s*a_m).             (1.1)

On every A_j, r cycles the three coordinates and s swaps the first two. The F action is free: nonidentity rotations translate x by theta or theta^2; reflections change the nonzero z to -z. Thus Z_m is a smooth affine geometric torus quotient of dimension2+2m. Its canonical translation-free torus cover is T_m: the linear S3 action on an A factor is faithful. Z_1 is the previously constructed universal quartic W. Geometrically T0=Gm, with the displayed split x coordinate. The character module is sign plus trivial plus m copies of the integral two-dimensional standard module.

There are two distinct atomicity questions. **Self-atomic** permits factorizations by endomorphisms of the same stabilized source variety. **Strong atomicity** permits different smooth normal affine intermediate varieties. In the finite setting UAFINATOM1 identifies the latter with primitive geometric monodromy. The distinction matters at degree9 below.

## 2. UABLOCKRED1: the general projection-block theorem

Let Z be smooth geometrically integral affine and A1-rigid in characteristic zero. By UA1RIGIDCYL1 any etale cylinder self-map and any cylinder automorphism projects to Z. Write

    G(z,v)=(f0(z),g(z,v))=(f0 x id) H(z,v),
    H(z,v)=(z,g(z,v)),  D=deg(f0), h=deg(H), M=D*h.             (2.1)

Suppose G has marked transitive label set Omega of size N. Partition Omega by equal projections to Z. This is a G-invariant block system. The b distinct projected points form one G orbit, each block has the same size L, and N=bL. If z is one projected point, its stabilizer G_z acts transitively on the L labels in that block. The invariant affine source fibre over z maps etale to the affine target fibre over y0=f0(z). Its complete marked fibre is exactly that block.

For L>1 the whole-Sylow theorem UASYLOWDEF1 applies to G_z, even if its action on the block has a kernel. Let h_z be this vertical polynomial map's generic degree. The normalization point budget of UAFIXDEFECT1 gives h_z<=h; similarly b<=D. Therefore

    M >= D*(L+2+indicator(2 divides L))
      >= N+b*(2+indicator(2 divides L)),       L>1.            (2.2)

No kernel of G_z is quotiented out of its polynomial action. An abstract faithful permutation quotient on the block need not act on the whole affine fibre. For L=1 this argument only gives M>=D>=b=N.

For a faithful primitive label action there are just two cases:

- b=1: the complete datum restricts to one genuine canonical affine-fibre datum on A^r, with the original G action;
- L=1: projection is injective on Omega, so any element acting trivially on Z acts trivially on Omega. Hence G embeds into Aut(Z).

Consequently, whenever the prescribed primitive G does not embed into Aut(Z), canonical data on Z x A^r exist only if they already exist on A^r. Over an algebraically closed field the converse takes the product with id_Z. This extends UAVERTRED1 from metabelian bases to **every** A1-rigid base satisfying the precise nonembedding condition. For split data the argument descends over a characteristic-zero field if Z has a rational point.

If G in(2.1) is finite, UA1RIGIDCYL1 makes H an automorphism. Every complete fibre then projects injectively, so **every faithful transitive label group embeds into Aut(Z)**, without assuming primitivity. Adding finite coefficient coupling does not weaken this requirement.

## 3. UAHOLCOMM1: a reusable automorphism criterion from holonomy

Let Z=T/F be a characteristic-zero geometric torus quotient, with its canonical translation-free torus cover. Write L=X_*(T) and

    E=End_(QF)(L tensor Q).

By UATORUSDOM1 every automorphism lifts and Aut(Z)=N_Aut(T)(F)/F. Conjugation on the deck group defines

    Aut(Z) -> Out(F),
    kernel = C_Aut(T)(F)/Z(F).                                (3.1)

The equality follows because an inner conjugation can be removed by multiplying the lift by a deck transformation. The remaining ambiguity is exactly the centre of F. No assumption that Out(F) is trivial or that F is centreless is needed here.

The affine torus centralizer C has an abelian translation kernel T^(F_lin), and its linear image lies in E^*. If E is commutative, this linear image is abelian. Thus C and the kernel in(3.1) are metabelian. In particular:

1. If Out(F) is solvable of derived length d, Aut(Z) is solvable of derived length at most d+2.
2. If Out(F)=1, Aut(Z) is metabelian.
3. For N>=5, an embedding S_N->Aut(Z) forces an embedding S_N->Out(F). Indeed, the kernel of its outer map is a metabelian normal subgroup of S_N, hence trivial: every nontrivial normal subgroup of S_N contains the non-solvable A_N.
4. For N=4 the outer image of any embedded S4 must be C2, S3 or S4. The only possible metabelian normal kernels are A4,V4 or1. In particular odd-order Out(F) excludes S4.

Combining with Section2 gives affine-fibre reduction for natural S_N data whenever these embedding tests fail, including arbitrary nonfinite coefficient coupling. Finite data are excluded. A useful sufficient hypothesis for commutative E is a sum of pairwise nonisomorphic absolutely irreducible rational constituents. Multiplicity-free over Q alone is insufficient if a constituent has noncommutative endomorphism division algebra.

For the family(1.1),

    E=Q x Q x Mat_m(Q).                                      (3.2)

At m=1 it is commutative, recovering the quartic metabelian obstruction UAQUARTAUT1. At m>=2 the argument no longer gives metabelianity. At m=3 the base actually admits a faithful S4 action: the standard integral representation of S4 on the rank-three augmentation lattice acts by monomial mixing of the three A copies, fixes x,z, and commutes with F. It descends faithfully since the centralizer intersects F only in its trivial centre. Section6 proves that even this genuine S4 subgroup does **not** supply the missing canonical quartic cover. Automorphism-group size is a necessary test, not a sufficient construction.

## 4. UAS3END1: complete endomorphism and stable self-factor classification

Every dominant geometric endomorphism of Z_m is finite etale, by UATORUSDOM1. It has a unique lift to T_m commuting with F. Indeed the induced automorphism of holonomy F is inner, so a deck adjustment makes the lift commute with F; centrelessness gives uniqueness. Conversely every dominant F-commuting torus map descends.

All such lifts are exactly

    x -> epsilon*x^d,
    z -> c*z^e,
    a_j -> zeta_j * product_i a_i^(M_ji),                     (4.1)

where

    d!=0, d congruent1 mod3;
    e!=0 odd;
    M is an integral m by m matrix with det(M)!=0;
    epsilon in mu2, c in Gm, zeta_j in diagonal mu3 subset A_j.

Powers and products of A points are coordinatewise. To prove completeness, torus units make every lift monomial times a constant. The rational sign and trivial summands cannot mix with the standard ones. The rational standard centralizer on m copies is Mat_m(Q); integrality gives Mat_m(Z), because End_(ZF)(L_A)=Z. Commutation with r gives d=1 mod3. Commutation with s gives odd e and epsilon^2=1. The A translations fixed by F are precisely the diagonal cube roots. These conditions are sufficient as well.

The degree and automorphism conditions are

    deg(f)=|d*e|*|det(M)|^2;
    f is an automorphism <=> d=1, |e|=1, |det(M)|=1.           (4.2)

There is a complete self-factorization invariant

    length(f)=Omega(|d|)+Omega(|e|)+Omega(|det(M)|),            (4.3)

where Omega counts prime factors with multiplicity, Omega(1)=0. Composition multiplies d,e and det(M), hence adds this length. The monoid of integers d=1 mod3 is freely generated by the signed primes q for q=1 mod3 and -q for q=2 mod3 (q!=3). Odd e factors into signed odd primes. Smith normal form factors every matrix M into integral matrices of prime determinant, separated by unimodular matrices. All translation constants in(4.1) can be put into automorphism factors. Consequently every noninvertible map factors into exactly length(f) self-atoms, and

    absolute stable self-atomic <=> length(f)=1.              (4.4)

For stabilization, any dominant factorization of a finite map on Z_m x A^r has finite factors by UAFINATOM1. Section1 of the retained rigid-base proof removes their vertical parts by source automorphisms. Their projected base maps compose to the original base map. The additive length gives (4.4) on every affine identity cylinder, over every extension of constants. No arbitrary different intermediate geometry is being identified with the source here.

### Strong atomicity and the exceptional prime3

Put K=kernel of the torus isogeny underlying(4.1). Translations do not affect it. The composite torus cover has group K semidirect F. It is free: a kernel translation cannot cancel the rotational theta on x since 3 does not divide d, or the reflection -1 on z since e is odd. Its faithful permutation monodromy is K semidirect image(F on K), acting affinely on K.

If this action is primitive, K must be elementary abelian, and the F module K must be irreducible. In Smith coordinates,

    K=mu_|d| x mu_|e| x product_i (mu_|s_i|)^2,

where s_i are the nonunit Smith factors of M. As an F_p module, this is the direct sum of sign and trivial one-dimensional pieces and copies of the two-dimensional standard S3 module. That standard module is irreducible exactly when p!=3: for p=3 its diagonal line in {(u,v,w):u+v+w=0} is invariant; for p!=3 there is no common invariant line for a rotation and a transposition (also in characteristic2).

It follows that strong absolute stable atomicity holds **exactly** in the following cases:

- |d| is prime and |e|=|det(M)|=1;
- |e| is an odd prime and d=1, |det(M)|=1;
- |det(M)| is a prime p!=3 and d=1, |e|=1.                (4.5)

The third case has degree p^2. An irreducible standard copy cannot be combined with another nonzero summand while remaining irreducible; Smith changes act on multiplicities and commute with F. These observations prove necessity as well as sufficiency, using UAFINATOM1 primitive-monodromy equivalence.

For det(M)=3 with d=e=1, the map is absolutely stably self-atomic but **not** strongly atomic. In particular the map cubing one A factor, and scaling z by any nonzero constant, has degree9. The diagonal subgroup

    C={(zeta,zeta,zeta) in A_1[3]} congruent C3

is F-invariant, so the map factors through the smooth affine variety

    Z'_m=(T_m/C)/F

by two finite etale degree-three maps. Z'_m cannot be isomorphic, or affine-stably isomorphic, to Z_m: otherwise that factorization would contradict (4.4). This is a concrete all-m boundary between self-atomicity and strong atomicity, rather than a failure of the degree test.

## 5. UAS3VERS1: an all-degree norm-torsor family with an exact obstruction at3

Let n>=2. Work over a field k with characteristic not dividing6n and, for the constant-group formulation, mu_n contained in k. Let

    G_n=A_1[n] semidirect S3 congruent (Z/n)^2 semidirect S3,

with S3 acting by the standard permutation action on product-one triples. Define

    a_m=2 if m is odd, 1 if m is even;
    b_m=1 if m is odd, 2 if m is even;   a_m*b_m=2,

and let Phi_(n,m) on Z_m be induced by

    (x,z;a_1,a_2,...,a_m)
       -> (x,z/n^b_m; a_1^n,a_2,...,a_m).                    (5.1)

This is a finite etale self-map of degree n^2 with geometric monodromy G_n in its natural action on (Z/n)^2. The larger free G_n action on T_m is the S3 action from(1.1) and n-torsion translations on A_1. Quotienting first by the kernel identifies the target with Z_m. Thus the rank-n^2 cover is associated to the finite G_n torsor T_m->T_m/G_n, and

    Phi_(a,m) o Phi_(b,m)=Phi_(ab,m).                         (5.2)

The following is the exact torsor-versality boundary of this family:

    every G_n torsor is realized over every extension field
    of the stated base <=> gcd(n,3)=1.                       (5.3)

For every stated base containing mu_n, necessity is witnessed over its rational-function extension when3 divides n; an individual field need not itself have an obstructed torsor.

For sufficiency let P be any G_n torsor over K and Q its S3 quotient. Let R/K be the corresponding cubic etale algebra. The twist of T_m splits into the following factors.

- The T0 factor is an affine torsor under the sign-twisted quadratic norm-one torus. Its class is killed by3 because cubing removes the rotational theta translation. H1 of that norm torus has exponent dividing2, so the class is zero.
- The z factor is a multiplicative torsor with sign translation, trivial by Hilbert90.
- The first A factor is a torsor under R^1_(R/K)Gm. Its class is killed by n because the nth-power map removes the A[n] translations. The group H1=K^*/Norm(R^*) has exponent dividing3. Coprimality of n and3 makes this class zero.
- Every other A factor is a linear twist of a torus and has its identity point; no translation class occurs there.

Hence the entire twist (T_m)_P has a K-point. For a free finite G action a point of the P-twist is equivalent to realizing P as a fibre of the G torsor. The associated S3-coset cover therefore realizes the full rank-n^2 algebra P times^(G_n) (G_n/S3). This includes finite fields and disconnected algebras. For n=2 the group is S4 and all quartic etale algebras occur, recovering UAQUARTSELF1. For n>2 we claim precisely these G_n-structured algebras, not all algebras of rank n^2.

For failure when3 divides n take the extension of the given base

    K=k(s,t),  R=K(s^(1/3)).

Since mu_n is contained in k and3 divides n, mu3 is contained in k.

This is a cyclic cubic extension; use its C3 torsor as an S3 torsor Q. At the t-adic valuation the extension is unramified of residue degree3, so valuations of norms are divisible by3. Thus t is not a norm, while t^3=Norm(t), and [t] is a nonzero order-three element of H1(K,R^1_(R/K)Gm). Since n kills this class, the Kummer exact sequence supplies a torsor under the n-torsion of this norm torus mapping to [t]. Combining it with Q in the split semidirect product gives a G_n torsor P. The first A factor of (T_m)_P has class [t] and no K-point. Thus P is not realized by the fixed G_n torsor. This proves (5.3) as a statement about specified G_n torsors. It does not assert that every unlabelled associated algebra distinguishes all torsors or all normalizer twists.

### Volume and atomicity

Let omega_0 be an invariant differential on T0, and omega_j an invariant two-form on A_j. The form

    Omega_m=z^a_m * omega_0 wedge dlog(z)
                    wedge omega_1 wedge ... wedge omega_m   (5.4)

descends to a nowhere-zero top form on Z_m. A transposition contributes sign -1 on omega_0 and on each omega_j, and sign (-1)^a_m on z^a_m; the product is1. Rotations preserve everything. Pullback by(5.1) multiplies the A_1 form by n^2 and the density by n^(-a_m*b_m)=n^-2. Therefore Phi_(n,m) preserves Omega_m exactly.

In every stated characteristic the explicit geometric monodromy gives

    Phi_(n,m) is strongly absolutely stably atomic
                                      <=> n is prime and n!=3. (5.5)

Primitive monodromy and finite-factor closure apply after every extension of constants and affine identity stabilization. In characteristic zero Section4 additionally gives absolute stable self-atomicity exactly when n is prime.

For p prime different from3, standard S3 on F_p^2 is irreducible, giving primitive monodromy F_p^2 semidirect S3 and strong atomicity. Composite n has an actual same-variety factorization from(5.2). For n=3 the explicit different-intermediate degree3-by-degree3 factorization in Section4 applies, although in characteristic zero the map remains absolutely stably self-atomic. No positive-characteristic all-endomorphism classification is inferred from the explicit family.

If mu_n is absent, the same twisting proof for gcd(n,3)=1 works with the finite etale group scheme A[n] semidirect S3, retaining its genuine cyclotomic action. It is not thereby a proof of versality for an arbitrarily identified constant group. We impose mu_n for the constant-group assertion and leave any sharper constant-field necessity outside this claim.

## 6. UAS3CAN1: repeated standard factors still cannot make a canonical quartic

For **every** m>=1, no finite etale endomorphism of Z_m, or of Z_m x A^r for any r with arbitrary finite coefficient coupling, carries canonical natural S4 data. This extends UAQUARTAUT1 far beyond the metabelian m=1 case. In particular Z_3 has a faithful S4 automorphism subgroup by Section3, yet no finite canonical quartic map.

It suffices to work on Z_m, by the finite cylinder theorem UA1RIGIDCYL1. Any such map has degree4 because its complete finite marked fibre has rank4. Formula(4.2) and the conditions on d,e leave exactly:

    d=4, |e|=1, |det(M)|=1;          kernel K=C4;
    d=1, |e|=1, |det(M)|=2;          kernel K=V4.             (6.1)

For the first case, the label S4 action would normalize the intrinsic regular abelian group K=C4, impossible since Hol(C4) has order8. To justify the normalization assertion generally, the regular group is the image of the intrinsic translation lattice Lambda in the geometric monodromy (UAINFRASELF1). A compatible target automorphism at the fixed target induces an automorphism of the based fundamental group, preserving characteristic Lambda. Equivariance intertwines its action with the source permutation on the complete fibre. Hence that permutation normalizes the image of Lambda. For torus quotients this image is the same K obtained from the canonical torus lift. This proves the asserted equivariant strengthening, without presuming K lies in the prescribed label group.

In the V4 case, Smith coordinates identify K with the standard S3 module F2^2: the deck holonomy F acts faithfully as GL2(F2). Its centralizer in Aut(K) is trivial. Every source or target S4 automorphism lifts **uniquely** to an automorphism of T_m centralizing F. This follows from Out(S3)=1 and Z(S3)=1 exactly as in Section4. The lifts form genuine S4 actions. The F-equivariant lift f of the endomorphism satisfies

    f alpha_g = beta_g f                                      (6.2)

for these source and target lifts: any deck discrepancy between the two maps centralizes F and is consequently trivial.

Choose a 3-cycle g in S4 and any lift t in T_m of the marked target y. Since y is fixed, beta_g(t)=a*t for some a in F. The lift beta_g has order3 and commutes with F, so a^3=1. Every F-centralizing torus automorphism has x output epsilon*x, by(4.1) with degree1. Its x multiplier is in mu2, hence beta_g fixes x. A nonidentity order-three element of F multiplies x by theta or theta^2 and cannot fix its nonzero value. It follows that a=1 and beta_g(t)=t.

The linear part of alpha_g restricts on K to an automorphism commuting with its faithful F action, so it is the identity on K. By(6.2), alpha_g preserves the K torsor f^-1(t) and acts on it by a translation. Its order divides3, whereas every nontrivial translation of V4 has order2. Thus alpha_g acts trivially on f^-1(t). This torsor identifies with the entire four-point quotient fibre on Z_m, and a natural 3-cycle cannot act trivially there. Contradiction.

The argument uses both the unique sign direction carrying the free order-three translation and the faithful mod2 standard kernel. Repeated standard representations enlarge the centralizer, but do not change that obstruction. We do not exclude other torus quotients with repeated sign directions, different holonomy, or different free translation classes; nor nonfinite coupled cylinders at arbitrary higher generic degree.

## 7. Scope and surviving requirements

This package gives an all-degree positive family with exact norm-cohomology versality and atomicity boundaries, plus general block and holonomy rules. A degree-nine absolute stable self-atom that factors through another smooth geometry illustrates why strong atomicity cannot be inferred from same-variety atomicity. A genuine S4 subgroup in an enlarged base illustrates why faithful automorphisms do not by themselves construct a canonical marked cover.

Still open are the original presentation-free affine descent problem, effective simultaneous fixed-locus landing and gluing, natural quartic canonical self-receivers on different geometries, nonfinite coupled maps on bases where the prescribed primitive group does embed, modular p-power cases and the A1-rich pseudo-plane cubic mixing gate. A smaller marked rank within a larger generic degree continues to require a separate ambient monodromy and atomicity proof.
