# Finite self-maps, torus quotients and a canonical volume lift

14 November 2026. The previous [homogeneous receiver pairs](../universal-keller-homogeneous-receivers/PROOF.md) cannot become self-maps by affine identity stabilization. Here we close a larger finite homogeneous self-map route, construct fixed universal cubic and quartic self-maps on other smooth affine geometries, and give an all-degree volume lift on the canonical line bundle. Whole affine-space canonical landing remains separate.

Strong atomicity allows dominant regular factors through arbitrary normal affine integral intermediates. Absolute stable atomicity allows constant extensions and affine identity stabilization. Canonical marked symmetry means actual compatible S_N actions and a fixed target with the complete natural N-point fibre; fieldwise occurrence of all algebras is a different requirement.

The homogeneous and general torus-quotient obstructions, and arbitrary finite coefficient coupling on cylinders, use characteristic zero. The cubic and quartic constructions work outside characteristics2,3. The surface endomorphism classification works outside characteristic2 with its stated etaleness test. The canonical-line lift works in every characteristic.

## 1. UAHOMOSELF1: finite commutators force abelian self-covers

Let Pi be a profinite group with finite closed commutator subgroup D. Every injective continuous endomorphism phi with finite-index image H satisfies
\[
                 \phi(D)=[H,H]=D,\quad D\subset H.           \tag{1}
\]
Indeed phi(D) is contained in D and has its finite order. Thus H is normal and Pi/H is abelian. A finite etale self-map supplies exactly such a finite-index fundamental-group embedding. Therefore finite geometric commutator subgroup forces abelian regular geometric monodromy, in every characteristic under this hypothesis.

Every geometrically homogeneous variety in characteristic zero satisfies it. UAHOMOCOST1 makes pi1(L/H) central-by-finite. Schur's theorem makes its commutator finite. The elementary transfer proof is short: choose finitely many central coset representatives. Every commutator is one of their finitely many commutators, so D is finitely generated. D intersection the centre C has finite index and is finitely generated abelian. Transfer to C sends c in C to c^n for n=[Pi:C] and kills D. Hence D intersection C has finite exponent, is finite, and D is finite. In a profinite group this finite abstract derived group is closed. This is classical Schur theory.

For finite endomorphisms of a characteristic-zero homogeneous affine variety, and every affine identity cylinder:

- fixed full-degree symmetric universality is impossible for N>=3, by UAUNIVMON1;
- full-degree canonical natural S_N symmetry is impossible for N>=5;
- a full-degree canonical quartic has imprimitive V4 monodromy and is not strongly absolutely atomic.

The last two statements use the full-degree normalizer rule and the normal abelian regular subgroup classification. They concern SELF-MAPS, whereas disconnected-isotropy homogeneous PAIRS remain available. Global finiteness is essential.

An abelian regular geometric cover descends to a torsor under a finite etale commutative constant-field group scheme. If one rational fibre is completely split, that group scheme is constant. Every rational-point fibre is then L^q for ONE abelian Galois field L, with q*[L:k]=N: its Galois action is translation by a homomorphism to the finite abelian label group, and all translation orbits are isomorphic. A mixed algebra k^(N-2) times a quadratic field over k(t) excludes symmetric universality for N>=3 without an inverse-Galois assumption.

## 2. UAINFRASELF1: torus quotients have a sharp symmetric ceiling

Suppose geometrically Z=T/F, where T is a torus and a finite group F acts freely by torus-VARIETY automorphisms; translations are allowed. In characteristic zero every finite etale self-cover of Z has an abelian regular NORMAL subgroup in geometric monodromy. The whole monodromy need not be abelian.

Over C, units make every automorphism of T a monomial integral map followed by translation. Its universal cover C^r gives an affine deck group Gamma with finite linear image. The kernel Lambda of that linear image is a free abelian translation lattice of rank r, containing the original torus lattice with finite index and spanning C^r.

Lambda is intrinsic: it is the unique maximal abelian normal subgroup. If g belongs to an abelian normal subgroup A and has linear part M, then for lambda in Lambda the commutator [g,lambda] is translation by (M-I)lambda and lies in A. Commuting g with this translation gives (M-I)^2 lambda=0. Finite-order M is diagonalizable in characteristic zero, so M=I. Thus A is contained in Lambda. The same proof applies to any finite-index subgroup Gamma', with translation lattice Lambda intersection Gamma'.

A finite self-cover identifies Gamma with Gamma'. It sends their intrinsic lattices isomorphically, so Gamma' has the full linear holonomy image. Consequently
\[
                              \Gamma=\Lambda\Gamma'.       \tag{2}
\]
The image of Lambda on Gamma/Gamma' is abelian, normal and transitive, hence regular. Riemann existence and characteristic-zero spreading transfer the finite permutation statement to every characteristic-zero base. The self-map itself was never assumed affine or equivariant.

For primitive monodromy the regular normal group is elementary abelian C_p^e: otherwise its characteristic subgroups give blocks. In the affine decomposition C_p^e semidirect the point stabilizer, primitivity is exactly irreducibility on F_p^e. Strong absolute atomicity for a finite affine cover is equivalent to primitivity: every intermediate block field normalizes to a finite etale affine cover, and UAFINATOM1 closes all dominant affine factors. Thus an atomic cover here has prime-power degree; prime-power degree alone is insufficient.

UAUNIVMON1 excludes fixed symmetric-universal degree N>=5. Full-degree canonical S_N is excluded too: its monodromy A_N or S_N has no abelian regular normal subgroup. Affine identity stabilization has the same characteristic-zero pi1, so the same intrinsic-lattice proof applies. This does not assert that a cylinder is itself a torus.

The cubic and quartic constructions below attain symmetric ranks3 and4. Thus "virtually abelian fundamental group" cannot replace finite commutator in section1. A finite free quotient of a homogeneous torus need not itself be homogeneous.

## 3. UANORMSELF1: a fixed universal cubic volume self-map

Outside characteristic2 put
\[
 V=\{u\ne0,\ t^2=u(s^2-4)\}.
\]
It is the smooth affine free quotient of Gm_x times Gm_y by
\[
 h(x,y)=(-x,y^{-1}),\quad u=x^2,\ s=y+y^{-1},\
                         t=x(y-y^{-1}).                    \tag{3}
\]
Adjoining x^2=u recovers y=(s+t/x)/2. This proves a finite etale double torsor globally, including t=0. Invariant Laurent units give O(V)^*=k^*u^Z. This nonconstant unit persists under affine identity stabilization.

Define C0=2,C1=s,Q0=0,Q1=1, with recurrence R_(n+1)=sR_n-R_(n-1). Then
\[
 C_n(y+y^{-1})=y^n+y^{-n},\quad
 Q_n(y+y^{-1})=(y^n-y^{-n})/(y-y^{-1}),
\]
\[
                   C_n^2-4=(s^2-4)Q_n^2,\quad C_n'=nQ_n.   \tag{4}
\]
The descent of (x,y)->(x,y^n) is
\[
                    F_n(u,s,t)=(u,C_n(s),tQ_n(s)).          \tag{5}
\]
It is finite of degree n, etale when n is invertible, and F_a F_b=F_ab. The double-torsor square is cartesian, so Chebyshev critical points cause no loss of etaleness.

The form Omega=du wedge ds/t extends globally, using du wedge dt/(u*s) off s=0. The two charts cover because s,t never vanish simultaneously. Its torus pullback is 2 dx wedge dy/y. Therefore F_n^*Omega=n Omega, and
\[
        G_n=(u/n^2,C_n(s),tQ_n(s)/n),\quad G_n^*Omega=Omega, \tag{6}
\]
with G_a G_b=G_ab. For n>=3 geometric monodromy is natural dihedral D_(2n), supplied by rotations y->zeta_n*y and reflection h on the torus normal closure. For n=2 the faithful group is C2. Prime n gives primitive strong absolute stable atomicity; composite n has the actual factorization(5)--(6).

ONE fixed G3 realizes every cubic etale algebra over every extension of every field of characteristic different from2,3:
\[
                   G3=(u/9,s^3-3s,t(s^2-1)/3).            \tag{7}
\]
It suffices to prove this for F3. For a cubic FIELD, Marques--Ward Corollary1.3 gives a generator r with r^3=a or r^3-3r=b. In the pure case r+r^{-1} is still a generator and has the second form with b=a+a^{-1}. A field generator has b!=2,-2. At target (b^2-4,b,b^2-4) the complete fibre is k[r]/(r^3-3r-b): r^2-1 is invertible and determines t uniquely.

The remaining cubic algebras are k times a quadratic etale algebra. The EXACT complete fibre at (u,2,0) is
\[
                           k\times k[z]/(z^2+3u).         \tag{8}
\]
The s=2 branch has t=0; the s=-1 branch has t^2=-3u. Taking u=-d/3 gives k times k[sqrt(d)], including the split algebra. This also proves finite-field universality: the t-coordinate separates the critical branches that a repeated univariate polynomial alone would miss.

More generally for odd invertible n, F_n realizes every degree-n algebra induced from constant dihedral D_(2n)-torsors exactly when the cyclotomic action on mu_n has image in {+1,-1}. Necessity follows from a split fibre: its norm-torus n-torsion must be constant. For sufficiency, fixing u gives the norm-one torus T_u of k[sqrt(u)]. Its n-torsion character is cyclotomic times the quadratic character of u. Choose u to match the input dihedral sign twist. H1(k,T_u)=k^*/Norm has exponent2, so has no n-torsion for odd n; the Kummer boundary is surjective. For n=3 the cyclotomic condition is automatic. This concerns induced algebras, without an asserted global choice of dihedral frame.

## 4. UANORMCYL1: complete endomorphisms and stable canonical distinction

Outside characteristic2 EVERY dominant regular endomorphism of V is, after choosing m>0,
\[
 f_{b,d,m,\eta}=(b^2u^d,\eta C_m(s),
                    b\eta u^{(d-1)/2}tQ_m(s)),             \tag{9}
\]
where b in k^*, d is nonzero ODD, eta=+1 or -1. Its degree is |d|m, and it is etale exactly when the characteristic divides neither d nor m.

Units force u->c*u^d, with d!=0 by dominance. Over a separable closure choose b^2=c and lift through(3) by x'=b*x^d. The other output is a Laurent unit a*x^e*y^m. If d is even, h-invariance forces m=0, contradicting dominance. If d is odd, h-equivariance forces e=0 and a^2=1. The t-output in(9) then forces b in k. Replacing (b,m) by (-b,-m) gives the same quotient map, so m>0 removes that ambiguity. The torus lift proves finiteness and the etaleness test, not just necessity.

For m>=3 geometric monodromy is C_|d| times D_(2m) in product action on |d|m labels. For m=1,2 it is abelian regular C_|d| times C_m. The torus normal closure has rotations mu_|d| times mu_m; h centralizes the first and inverts the second. Its reflection core is trivial exactly for m>=3. Every composite degree factors on V: remove the output automorphism, separate x->x^d and y->y^m, then factor the exponents. Thus every etale map(9) is strongly absolutely stably atomic exactly when its degree is prime.

In characteristic zero EVERY finite etale endomorphism of V times A^r is
\[
                            (f\times id)\circ H,           \tag{10}
\]
with f from(9) and H an automorphism preserving V projection. Indeed units determine u; after the double-torus pullback every torus output is a Laurent unit independent of the affine coordinates. Its V-output is already(9). The remaining H is finite by UAFINATOM1 and etale. Its geometric vertical fibres are finite etale self-maps of A^r, hence isomorphisms in characteristic zero. This is the finite vertical argument of UATORUSSTABLE1 on a NONHOMOGENEOUS torus quotient. Positive-characteristic finite vertical maps are outside(10).

This gives the exact finite-cylinder classification in characteristic zero, up to source and target automorphisms:

| Rank | Fixed symmetric universality | Full natural canonical S_N |
|---|---|---|
|2|d=+/-1,m=2|yes|
|3|d=+/-1,m=3|exactly d=+/-3,m=1, with mu3 in k|
|N>=4|none|none|

Fixed universality requires geometric S_N by UAUNIVMON1. The displayed groups allow only C2 and S3; rank4 gives D8, and higher S_N has no abelian regular normal subgroup. Quadratic universality follows by restricting u=1 to a split torus; cubic universality is section3. Normality in S_N excludes canonical rank4 and every higher rank as well.

For the cubic map with d=+/-1,m=3, no arbitrary affine identity stabilization can become canonical. An order-three TARGET action fixing the mark must fix u: its action on units is u->c*u and u never vanishes. Equivariance forces the SOURCE C3 to fix u too. Over kbar(u), geometric units make its automorphism on the norm torus times A^r descend to that torus independently of affine coordinates. The generic quadratic norm torus has no rational3-torsion: its Galois action on mu3 is inversion. Its order-three automorphism is therefore trivial. The source action preserves V projection and fixes every label in the complete stabilized fibre, a contradiction.

For d=+/-3,m=1 the three source u-values are distinct. If mu3 is absent from k, every source order-three unit action fixes u and cannot permute them, after any stabilization. If zeta3 is in k, then
\[
                            E3=(u^3,s,u*t)                 \tag{11}
\]
has the required S3 action, generated on the source by
\[
 (u,s,t)\mapsto(\zeta_3^2u,s,\zeta_3t),\quad
 (u,s,t)\mapsto(u^{-1},s,t/u).
\]
On the target the rotation is trivial and the reflection is the same. The mark (1,2,0) has full fibre u^3=1,s=2,t=0. Rank2 uses the free involution (u,s,t)->(u,-s,-t).

No one finite cubic map in these cylinders is BOTH fixed universal and canonical, although different prime-degree atoms attain each property.

Formula(9) pulls Omega back by b*d*m*u^((d-1)/2). Every global nonvanishing volume form is c*u^j*Omega. Preservation forces d=1, since the exponent is (d-1)(j+1/2). Hence the canonical E3 category preserves no global volume on V. Section6 removes this obstruction with ONE extra affine line. Over a characteristic-zero field containing mu3, one added line is exactly necessary and sufficient for a finite canonical volume cubic within this cylinder category. This is a cost for this geometry, not a global coordinate minimum.

## 5. UAQUARTSELF1: a fixed universal atomic quartic volume self-map

Over any field of characteristic different from2,3 let T0 be the cyclotomic norm-one torus a^2-a*b+b^2=1. The element theta with theta^2+theta+1=0 is its rational point of order3, with translation and inversion as in UATORUSSYM1. Set
\[
 B=T0\times Gm_z,\qquad
 A=\{(v1,v2,v3)\in Gm^3:v1v2v3=1\}.
\]
Let S3=<r,s> act on B by
\[
                   r(x,z)=(\theta x,z),\quad
                   s(x,z)=(x^{-1},-z),                    \tag{12}
\]
and on A by a corresponding3-cycle and transposition of coordinates. Its diagonal action on the FOUR-dimensional torus T=B times A is free: rotations translate T0 freely, and reflections change nonzero z to -z. Thus W=T/S3 is smooth geometrically integral affine.

The free base is also the REGULAR GALOIS CLOSURE of the cubic surface map, so the two constructions use one mechanism. Write (a0,b0) for T0 coordinates. The quotient of B by the reflection s is exactly V via
\[
               u=-z^2/3,\quad s0=2a0-b0,\quad t=z*b0.
\]
Indeed s0^2+3b0^2=4, and adjoining z^2=-3u recovers b0=t/z and a0=(s0+b0)/2. Cubing T0 sends s0 to C3(s0) and b0 to b0*Q3(s0). It kills the constant C3 translation subgroup, so B/S3 is the same V and the subgroup cover B/C2->B/S3 is exactly F3. The quartic construction reuses this free cubic-frame torus as its base.

The map id_B times[2]_A descends to a finite etale SELF-MAP of W of degree4. Normalize it by z->z/2, obtaining Phi induced by
\[
                (x,z;v1,v2,v3)\mapsto
                         (x,z/2;v1^2,v2^2,v3^2).           \tag{13}
\]
The signs A[2]=V4 and S3 permutations form S4. This affine S4 action on T is free: nontrivial S3 components are free on B, and nontrivial V4 translations are free on A. Quotienting first by V4 identifies T/V4 with T through squaring A, compatibly with S3. Therefore T/S4 is the SAME target W, and the unnormalized map is the subgroup cover T/S3->T/S4. Its full geometric group is natural S4, hence Phi is strongly absolutely stably atomic despite degree4.

### Every quartic algebra is a complete fibre

For an S4-torsor P over any extension field K, it suffices to show that its twist T_P has a K-point. Then the free torsor quotient T->T/S4 is field-versal and the S3-quotient fibre is precisely the input quartic algebra. This includes finite fields.

Let Q=P/V4 be the S3 torsor, and R its cubic resolvent. Twisting the affine T0 action in(12) gives a torsor under the sign-twist of T0. Its class is killed by3: cubing kills the translation by theta and maps to the linear sign action. But H1 of this quadratic norm-one torus has exponent2, so its class is zero. The z factor is a Gm torsor, trivial by Hilbert90. The Q-twist of B therefore has a point.

The P-twist of A is a torsor under the norm-one torus T_R of the cubic resolvent. Its class is killed by2 because squaring maps to the linear S3 twist. On the other hand
\[
                      H^1(K,T_R)=K^*/Norm_{R/K}(R^*)
\]
has exponent dividing3, also for disconnected R. This class is zero too. Twisting the product gives the product of these two twists; thus T_P(K) is nonempty and every quartic algebra is a COMPLETE fibre on this ONE fixed Phi.

The coprime torsor cancellations are the classical norm-cohomology inputs already used in UATORUSSYM1. The FREE S3 base packages the varying resolvent forms into one fixed SELF-MAP. Fieldwise points do not choose natural sections over general schemes, and we do not infer a canonical full S4 action on Phi.

### Preserved global volume and exact rank ceiling

Choose a nonzero invariant1-form omega0 on T0 and omegaA=dlog(v1) wedge dlog(v2) on A. The form
\[
                   z^2\,\omega0\wedge dlog(z)\wedge\omegaA \tag{14}
\]
is S3-invariant: a transposition changes both omega0 and omegaA by -1, and fixes z^2. It descends to a nonvanishing volume on W. Under(13), squaring A contributes4 and z->z/2 contributes1/4 through z^2. Thus Phi preserves that volume exactly.

The nonconstant invariant unit z^2 descends to W and persists under affine identity stabilization. W is not ordinary affine four-space and cannot become affine space this way.

Together with section2 this is sharp: fixed universal finite self-maps on characteristic-zero geometric torus quotients attain symmetric ranks2,3,4 and exclude every N>=5. Ranks2,3 occur on V and rank4 on W. These are attainable ranks, not dimension minima or a classification of canonical quartic self-maps.

## 6. UAKANLIFT1: a canonical volume lift for every etale receiver

Let f:X->Y be etale between smooth d-dimensional varieties, in ANY characteristic. On their canonical covector line bundles there is a natural lift
\[
               \widehat f:\operatorname{Tot}(K_X)
                                      \longrightarrow\operatorname{Tot}(K_Y), \tag{15}
\]
sending xi at x to the unique covolume at f(x) whose df_x-pullback is xi. The total space of K_X has a tautological d-form tau_X. In local coordinates xi=w dx1 wedge...wedge dxd,
\[
 \tau_X=w\,dx1\wedge\cdots\wedge dxd,\quad
       \Xi_X=d\tau_X=dw\wedge dx1\wedge\cdots\wedge dxd.    \tag{16}
\]
Thus Xi_X is global and nowhere zero, including on the zero section. The definition gives hat f^*tau_Y=tau_X and hat f^*Xi_Y=Xi_X. There is no factorial denominator or characteristic condition.

The isomorphism f^*K_Y=K_X makes(15) a BASE CHANGE of f. Hence it is etale, finite when f is finite, with the same degree and geometric monodromy. The latter is just independent-transcendental extension in the target line coordinate. If X,Y are affine their line total spaces are affine, so UAFINATOM1 preserves the strong absolute stable atomicity criterion.

This is a functor: hat(f g)=hat f hat g. Every compatible group action lifts naturally. A target mark becomes (y,0), and its COMPLETE fibre is the original full fibre on the zero section. Fixed universality and canonical marked gluing are preserved, without a choice of primitive generator or volume. Relative canonical bundles and relative exterior derivative give the same statement over arbitrary smooth scheme families.

If X=Y and f^*Omega=j_f Omega for a chosen volume, j_f is a regular unit and the formula is
\[
                         (x,w)\mapsto(f(x),w/j_f(x)).      \tag{17}
\]
Compensation changes the map by a regular source automorphism; preservation belongs to the compensated map. It is source-automorphism equivalent to f times identity. For E3 this becomes
\[
                      (u,s,t,w)\mapsto(u^3,s,u*t,w/(3u)).  \tag{18}
\]
The derivative lifts of its S3 actions preserve the canonical volume, and the zero target mark has its entire original three-point fibre. This proves the sharp one-line volume extension in section4.

This is standard canonical-bundle/tautological-form geometry used as a receiver functor, without a novelty claim. It supplies volume preservation after receiver geometry is available. Units on a line total space equal the units on its base; characteristic-zero pi1 is also unchanged. Thus the unit obstructions of V,W and the source/target finite-pi1 distinction of the old homogeneous PAIRS survive. For a polynomial Keller self-map already having j=1, the lift is just identity stabilization and cannot repair its canonical obstruction.

## 7. Resulting boundary

Finite homogeneous SELF-MAPS cannot be symmetric-universal from rank3 onward. The larger finite torus-quotient class attains exactly ranks2,3,4, with affine-type primitive monodromy. We now have fixed universal volume self-maps in cubic and quartic rank, including a composite-degree absolutely atomic quartic. The cubic model has complete finite-cylinder and canonical-volume classifications. Every smooth etale receiver has a degree-preserving canonical volume lift in all characteristics.

Higher-rank finite self-receivers must escape geometric torus quotients. The original whole affine-space canonical gate still needs its genuinely higher-degree nonfinite boundary, actual simultaneous fixed loci and effective source landing. A canonical full S4 action on a suitable quartic self-receiver remains open here. No new coordinate minimum or publication-priority assertion is made.

The exact controls check norms, Chebyshev identities, full critical fibres, finite-field splitting types, actual monodromy and volume identities. They do not formalize the all-degree group-topology argument, endomorphism classification, etale descent or arithmetic cohomology. Independent whole-proof review, formal verification and external review remain pending.
