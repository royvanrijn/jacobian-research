# Canonical all-rank affine Keller fibres by invariant filtering

15 November2026. Written universal proofs with narrow exact algebra controls. Independent whole-proof replay, formal verification of the new construction, external review and publication priority remain unclaimed. OP-KDESC is narrowed substantially: canonical marked polynomial SELF-MAPS of vector bundles with whole affine-space fibres now exist in every rank. Atomicity and canonical trivialization to one fixed affine space over an arbitrary base are separate requirements and are not supplied.

The construction uses only the finite etale algebra, its unit, multiplication and trace. It chooses no primitive generator, generator line, ordering of roots, or point on an input torsor. Its generic degree is intentionally larger than the marked rank. It has an explicit nontrivial factorization; canonical existence must not be reported as canonical atomic existence.

## 1. A three-dimensional seed with empty, singleton and split triple fibres

Over R0=Z[1/6] put

    t=1+xy, q=t^2*z+y^2*(1+3t),
    K(x,y,z)=(t*q, -(y+3x*q)/2, x*(5-3t)-x^3*z).           (1.1)

This is the determinant-one normalization of the retained F1/UAM1 cubic. Its coordinate degrees are7,6,4 and its generic degree is3. Direct polynomial differentiation gives det(DK)=1. The original seed is published and independently formally verified; these facts do not formally verify the new receiver below.

The exact complete SCHEME fibres used here are

    K^-1(mu)=empty,       mu=(1/3,-1,2/3),
    K^-1(0)={0},
    K^-1(rho)={p0,p1,p2}, rho=(-1/4,0,0),                 (1.2)

where

    p0=(0,0,-1/4), p1=(1,-3/2,13/2),
    p2=(-1,3/2,13/2).                                     (1.3)

These identities hold over R0, including after arbitrary base change, not merely on rational points over Q.

Here is an algebraic verification that also records the necessary denominators. Write the target as(a,b,c). On a invertible, t and q are units. For S=x/t the exact source identities give

    E(S)=S+b*S^2+a*S^3-c/2=0,   t*E'(S)=1.                (1.4)

The cleared identities are

    2*x*t^2+2*b*x^2*t+2*a*x^3-c*t^3=0,
    t^2+2*b*x*t+3*a*x^2=t.

Conversely their standard cubic reconstruction gives x=S/E', t=1/E' and then y,z; this is the retained UAM1 inverse chart. At mu, 3E=(S-1)^3 and E'=(S-1)^2. Nilpotence contradicts tE'=1, so that fibre is empty over every R0-algebra. At rho, E=S-S^3/4 has the three distinct roots0,2,-2, with E' respectively1,-2,-2; reconstruction gives(1.3). Equivalently its fibre ideal is

    (x^3-x, y+3*x/2, z+1/4-27*x^2/4).                    (1.5)

For the zero fibre use the global comaximality identity

    x^2*q-1=t*(x^2*t*z+3*t^2-5*t+1).

Since tq=0, the fibre splits into t=0 and q=0 parts. At t=0, xy=-1, q=y^2 and K2=y, contradicting K2=0. At q=0, K2=-y/2 forces y=0, then t=1, q=z=0 and K3=2x=0. This proves the singleton scheme assertion. The generic degree3 follows from the independent target coefficients in(1.4) and reconstruction; it persists in every characteristic different from2,3.

## 2. UAFILTER1: remove an arbitrary invariant part of a finite fibre

Let k have characteristic zero, G be finite, and F:V=A^d->W=A^d be a G-equivariant polynomial Keller map, with a G-fixed target y and complete finite etale fibre E. Let E0 be any G-stable clopen subscheme of E defined over k. There is a G-equivariant Keller map on the WHOLE affine spaces V times A^3 -> W times A^3 whose complete fibre over(y,0) is canonically E0. Its generic degree is3 times deg(F).

Indeed the characteristic idempotent1_(E0) in O(E) lifts to a polynomial on V. Average a lift over G to obtain a G-invariant polynomial P with values1 on E0 and0 on its complement. This is ordinary finite-set interpolation and the Reynolds operator; it uses no claim that a general free action has a natural point. Define

    Fhat(v,u)=(F(v), K(u)-(1-P(v))*mu).                    (2.1)

Let G act trivially on the three new coordinates on both sides. Invariance makes(2.1) equivariant. Its differential is block triangular with diagonal DF and DK, so its determinant is1. On the first-coordinate fibre E, the desired components have second equation K(u)=0, while the others have K(u)=mu. Equation(1.2) proves the COMPLETE scheme fibre assertion. No open is deleted from the source.

There is an actual factorization

    Fhat=(F times id) composed with
          [(v,u)->(v,K(u)-(1-P(v))*mu)].                    (2.2)

The latter factor has degree3: over k(V) it is a target translation of K. The first has degree deg(F), so the total degree is3 deg(F). If deg(F)>1, these are two noninvertible regular factors. Thus the filter supplies existence and exact fibre selection, not atomicity.

## 3. UACANFIB1: a canonical construction for every finite etale algebra

Fix N>=2 and work over

    R_N=Z[1/(2N)!].                                       (3.1)

In particular all characteristic-zero bases and all fields of characteristic p>2N are permitted. Let B be ANY R_N-scheme, including a nonreduced one, and let A be a finite etale O_B-algebra of rank N. Its underlying vector bundle is part of the input, not a chosen free presentation. Put

    V_A=A^3 direct_sum O_B^3,
    L(X)=(3X^2-X)/2,
    C_A(X)=Tr_(A/O_B)(L(X)),
    lambda_N(T)=-T*product_(j=2,...,2N)(T-j)/(2N-1)!.       (3.2)

Thus lambda_N(j)=1 if j=1 and0 for every other integer0<=j<=2N. Let K_A denote(1.1) evaluated using the unit and multiplication of A. Define the polynomial vector-bundle SELF-MAP

    F_A(X,Y,Z;u)=
      (K_A(X,Y,Z), K(u)-(1-lambda_N(C_A(X)))*mu).           (3.3)

Its marked target is the section

    y_A=(rho*1_A;0).                                      (3.4)

This is a presentation-free canonical marked Keller construction of relative dimension3N+3. It commutes with arbitrary base change and every isomorphism of finite etale algebras. Source and target are the SAME specified vector bundle. All geometric fibres are actual whole affine spaces, and the relative Jacobian determinant is1.

### Whole marked fibre and absence of hidden choices

Etale-locally on B, A=O_B^N and K_A is the product of N copies of K. By(1.2), the whole first-coordinate fibre over rho*1_A is the split finite scheme {p0,p1,p2}^N. The first coordinates of p0,p1,p2 are0,1,-1, and L has respective values0,1,2. A configuration with b entries p1 and c entries p2 therefore has

    C_A=b+2c,  0<=b+2c<=2N.

The selector equals1 exactly when b=1,c=0. Thus(3.3) retains exactly the configurations with one p1 and every other entry p0; it removes every other configuration using the EMPTY seed fibre, while the retained extra fibre is the singleton0. These N configurations form the natural S_N-set. This is an identity of schemes on the split cover: both seed fibre statements hold universally over R_N. Faithfully flat descent proves

    F_A^-1(y_A) canonically isomorphic to Spec_B(A).        (3.5)

The identification can also be written intrinsically. The primitive rank-one idempotent e in A over a test B-scheme gives the fibre point

    (X,Y,Z;u)=(e,-3e/2,-1_A/4+27e/4;0).                   (3.6)

Here e^2=e, L(e)=e and Tr(e)=1. Since the integers up to N are separated over R_N, this trace condition selects precisely rank-one idempotents. Their finite etale scheme is Spec_B(A): over a split algebra they select one component, and this description descends canonically. Formula(3.6) chooses no idempotent over B. It identifies the entire idempotent SCHEME with the input cover. A connected input field need not have any such point over its own ground field.

Every operation in(3.2)-(3.4) is intrinsic to the algebra and commutes with base change. This proves effective polynomial gluing over all the stated bases, not just fieldwise point existence or a reduced-base argument. Classically, this is restriction of scalars along Spec(A)->B; the relevant sheaf, base-change and etaleness properties are in Stacks Section97.11. The explicit split calculation proves the needed representability here as a vector bundle.

### Actual inertia and a receiver over B S_N

Over the split algebra, S_N permutes the N three-coordinate blocks and fixes the last three coordinates. It acts LINEARLY and faithfully, the map is equivariant, y_A is fixed, and its full fibre is the natural faithful N-point set. Associated vector bundles and formula(3.3) therefore give a section over B S_N into the stack of marked polynomial etale vector-bundle self-maps just specified. The full S_N inertia is retained.

This receiver forgets the primitive generator line. It does not map through the line-decorated receiver excluded by UAINERT1, and its generic degree is much greater than N, so it respects UAFIXDEFECT1. Unlike the UARECV1 presentation atlas it requires no good generator or admissible translation from the input algebra. It is a different compiler, not a stabilization or left-right modification of the old degree-N map.

The bundle V_A need not have a canonical trivialization on an arbitrary B. This theorem constructs a canonical SELF-MAP of that bundle, not a natural choice of a basis or a single fixed map on a trivial A_B^(3N+3) realizing all scheme families. In particular it does not contradict UAVERSLIMIT1. Over a field the source and target are affine spaces of the stated dimension; choosing coordinates is an isomorphism, not additional algebra-presentation data in the functor.

### An invariant global frame is already obstructed

A natural frame fixed by every input automorphism is impossible for this receiver, even after adding any number of trivial identity lines. At the split algebra O_B^N the input automorphism group S_N acts faithfully on the A^3 summand. A frame natural under those automorphisms would identify that representation, plus trivial summands, with a framed constant vector space on which all input automorphisms act as identity. Naturality would make the faithful action trivial, a contradiction. The volume line does not remove the faithful summand. This is the elementary inertia obstruction to an invariant coordinate frame, not an obstruction to the map-valued bundle functor just constructed or to choosing coordinates noncanonically over a field.

## 4. UACANFIX1: all subgroup-fixed restrictions are realized simultaneously

Work with the split S_N model in characteristic zero, or characteristic p>2N. Let H<=S_N have s orbits of sizes n1,...,ns on the labels and r singleton orbits. Its WHOLE source and target fixed schemes are the affine spaces

    (V_A)^H=A^(3s+3).

The actual restricted map is(3.3) with one triple for each orbit and trace input

    C_H=sum_(j=1,...,s) n_j*L(x_j).                        (4.1)

It is Keller with determinant1 and generic degree

    d_H=3^(s+1).                                          (4.2)

Its complete marked fibre has exactly r points, naturally the H-fixed part of the N labels. Indeed(4.1) equals1 on the product seed fibre exactly when one SINGLETON orbit has label p1 and all other orbits have label p0. Every non-singleton p1 or any p2 contributes at least2. This proves the statement for every subgroup at once, without enumerating them or separately gluing fixed loci.

In particular every deranging subgroup has an actual empty marked fixed fibre between whole affine spaces. A transitive H has s=1, giving one explicit six-dimensional Keller fixed restriction of generic degree9 and empty marked fibre. The simultaneous geometry requested by UAFIXDEFECT1 and UASYLOWDEF1 is therefore attainable in this higher-degree, non-atomic construction. No exotic fixed surface or nonlinearizable group action is needed. The previous lower bounds remain valid and make no assertion of nonexistence at this large surplus.

## 5. UACANVOL1: dimensions, degree, canonical volume and exact failure of atomicity

Every geometric fibre of F_A has generic degree

    M=3^(N+1),                                            (5.1)

since the N seed factors have degree3^N and the last factor has degree3. In the split coordinates its ordinary coordinate degree is exactly4N: lambda_N(C_A) has that degree, while K has coordinate degrees7,6,4. The trace quadratic is nonzero and the selector leading coefficient is a unit over R_N. No optimality of these bounds is claimed.

### A canonical invariant volume with one more line

The determinant-one relative map on V_A need not come with a chosen trivialization of det(V_A)^*. Apply UAKANLIFT1's WHOLE canonical-line functor. Since the relative canonical line of a vector bundle is pulled back from its determinant dual, the resulting source and target are the same vector bundle

    V_A direct_sum det(V_A)^*,    rank3N+4.                (5.2)

Its determinant is canonically trivial by evaluation. The canonical lifted map preserves the tautological volume, retains the exact fibre(3.5) at the zero covolume section, and keeps generic and coordinate degrees. In split coordinates the extra S_N line is the sign representation: three permutation copies contribute sign^3 and the dual determinant contributes sign, giving total determinant1. This is a whole affine-line volume lift, not a punctured line with a new unit. All source fibres remain affine spaces.

### The remaining atomicity gate is real

Formula(3.3) factors, over EVERY permitted B and functorially in A, as

    H_A:(P,u)->(P,K(u)-(1-lambda_N(C_A(X)))*mu),
    R_A:(P,v)->(K_A(P),v),
    F_A=R_A H_A.                                          (5.3)

Both are polynomial self-maps of V_A with determinant1, and have geometric degrees3 and3^N. They are equivariant under all input automorphisms. Thus F_A is not self-atomic, even in the equivariant category. The volume lift remains decomposable because both factors have their canonical volume lifts. This factorization persists under every extension of constants, identity stabilization and polynomial left-right equivalence; none of those operations can repair atomicity of this example.

Over geometric constants, the maximum number of nonautomorphism factors is exactly N+1. The N separate seed blocks and the last translated seed exhibit that length; all have geometric prime degree3 and are absolutely stably self-atomic by UASELFATOM1. Conversely regular factors of a polynomial Keller map are Keller, every degree-one factor is an automorphism, and the nontrivial factor degrees are powers of3 whose exponents sum to N+1. This proves the upper bound. The same geometric argument applies after stabilization. For nonsplit A the individual ordered block factorization need not descend, while the two canonical factors(5.3) always do. No full wreath-product monodromy assertion is inferred from the degree tower.

The generic monodromy is consequently imprimitive; no automatic S_M or atomicity follows from the much smaller marked rank N. An atomic canonical construction must use a different coupling that destroys this actual intermediate extension while preserving its prescribed complete fibre. Merely composing more filters, adding identity variables or changing coordinates cannot achieve that.

## 6. Precisely which open requirements are now closed

For the newly specified polynomial vector-bundle receiver, canonical marked realization in every rank, whole affine-space geometric source fibres, actual natural S_N symmetry, identification of source and target, and effective gluing over all R_N-schemes are constructed. The all-subgroup fixed restrictions exist simultaneously, and canonical volume costs one explicitly identified line. These were separate existence gates; they are not left as an unexecuted plan.

Still open are a canonical ATOMIC version, substantially smaller generic degree and dimension, optimality, and questions about existential freeness of the underlying bundle on prescribed bases. A natural frame fixed by all input automorphisms is already obstructed as in Section3, and a single fixed quasi-compact all-scheme universal map retains the published UAVERSLIMIT1 obstruction. Neither should be relabelled an unresolved positive construction gate. OP-KDESC is therefore kept open for its strengthened atomic/degree-controlled and receiver-comparison requirements, while its earlier unqualified non-atomic polynomial-bundle existence gap must no longer be reported as open. The old degree-N compiler remains obstructed: this construction changes generic degree to3^(N+1) and does not transport its stable class.

No assertion about JC(2), HC(4), full generic S_M monodromy, atomicity of this receiver, or publication priority is made. The all-rank proofs are the exact seed fibres, invariant interpolation, trace selection and ordinary descent; small exact controls are not a substitute for them.
