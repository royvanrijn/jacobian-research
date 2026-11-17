# Versal isogeny self-covers: arbitrary free bases, exact cohomology and the scheme boundary

15 November2026. Written arguments with narrow exact controls; no formal verification, independent whole-proof replay or priority certificate. The original affine-space canonical descent problem remains open. References and retained inputs are pinned in `sources.json`.

## Conventions

A torsor is **field-versal** here if every torsor for the specified group over every extension field, including finite fields, is a pullback at a field-valued point. This is the all-fields point-existence version of weak versality. It does not assert a natural choice of that point, density on every invariant open, or versality over general affine schemes. Duncan--Reichstein's twist criterion is the classical framework. Source/target maps below are regular self-maps of one fixed smooth affine variety, but that variety has nonconstant units and is not affine space.

Finite dominant affine factors are finite (UAFINATOM1). Strong absolute stable atomicity permits different normal affine intermediate varieties and is equivalent in this finite setting to primitive geometric monodromy. Absolute stable **self**-atomicity allows only endomorphisms of the same identity-stabilized variety. These two properties must be kept distinct, as in UAS3END1.

## 1. UAFREEBASE1: a free field-versal torus base for every finite group

Let F be a finite constant group over k, with char(k) not dividing |F|. There is a k-torus B, with an affine monomial F action, such that:

- the action is free on the entire torus;
- for every F torsor Q over every extension K/k, the twist B_Q has a K-point;
- the linear torus underlying B_Q is quasi-trivial. In particular B_Q is a rational torsor under a quasi-trivial torus; over infinite K its points are dense.

The construction is elementary induced-character geometry, not an asserted new classification of versal actions. Choose one subgroup H of order ell for each conjugacy class of prime-order subgroups of F. Put

    E_H=k(zeta_ell),  C_H=Res_(E_H/k) Gm.

A chosen generator of H acts on C_H by multiplication by the rational point zeta_ell of C_H(k), of exact order ell. Form the multiplicatively coinduced affine H-variety

    B_H=Coind_H^F(C_H),  B=product_H B_H.                     (1.1)

After choosing coset representatives, B_H is a product of [F:H] copies of C_H. If g*t_j=t_i*h, its action sends the jth coordinate to the ith coordinate multiplied by the chosen character of h. The linear part only permutes coordinates. This is an algebraic action over k because zeta_ell is a k-point of the restriction-of-scalars torus; no root of unity is assumed to belong to k itself.

If h has prime order and is conjugate into H, it fixes a coset coordinate and multiplies that coordinate by a nonidentity root of unity. Thus it has no fixed point on B_H. A fixed point of any nonidentity g would also be fixed by a prime-order power of g. This proves freeness of F on B.

For any Q, the affine twist (B_H)_Q is a torsor under the twist of the linear permutation torus. Arithmetic Galois permutes the embeddings of E_H and F permutes the cosets, so the resulting character lattice is still a permutation lattice. The torus is quasi-trivial; Hilbert90 gives H1=0 over any field. Therefore every (B_H)_Q, and hence B_Q, has a point. Its noncanonical trivialization also proves rationality. The dimension is explicitly

    dim(B)=sum_H [F:H]*[k(zeta_ell):k].                       (1.2)

No minimum dimension is claimed. B is an invariant torus open in the corresponding induced linear representation, but its global freeness and finite-field point assertion have been proved separately rather than inferred from generic freeness.

For cyclic F=C_m one can use the smaller base B=Res_(k(zeta_m)/k)Gm with the generator acting by multiplication by zeta_m. It is free, and all its affine twists are Hilbert90-trivial. This works for every char(k) not dividing m.

## 2. UAISOCRIT1: the exact cohomological criterion for a fixed self-cover

Let S be a split k-torus with an action of F by group automorphisms. Let phi:S->S be an F-equivariant isogeny of degree D, with char(k) not dividing D|F|. Put A=ker(phi) and let the finite etale group scheme G=A semidirect F act on B x S: A translates S and F acts diagonally. The action is free because F is free on B and nonidentity elements of A translate S freely. Set

    X=(B x S)/F.

The map id_B x phi descends to a finite etale self-map f:X->X of degree D. Quotienting B x S first by A identifies the target (B x S)/G with the same X. The G torsor B x S -> X has f as its associated F-coset cover.

For an F torsor Q/K, write S_Q for the twisted torus and phi_Q its endomorphism. Then

    the fixed G torsor is field-versal
       <=> ker[phi_Q:H1(K,S_Q)->H1(K,S_Q)]=0
           for every K/k and every F torsor Q/K.             (2.1)

Proof. A G torsor P above Q twists S into a torsor under S_Q. The possible classes are exactly the image of H1(K,A_Q) in H1(K,S_Q), hence exactly ker(phi_Q), by the cohomology sequence of the isogeny. Every such A_Q cocycle combines with Q to give a G cocycle because G is the split semidirect product. The B factor always has a point by Section1. Thus the P twist of B x S has a point precisely when its S_Q torsor class vanishes. The usual free-action twist criterion proves both directions of(2.1), including necessity for specified G torsors. Associated unlabelled algebras need not determine their G structures uniquely.

A useful sufficient condition is

    gcd(D,|F|)=1.                                            (2.2)

Indeed Q splits S_Q, and restriction/corestriction kills H1(K,S_Q) by |F|. The integral adjugate isogeny psi, commuting with F, satisfies psi*phi=[D]; therefore ker(phi_Q) is killed by D. Bezout gives zero when(2.2) holds. More sharply, |F| can be replaced by any uniform annihilator c of all the groups H1(K,S_Q). This is why the cubic norm example works at D=4 although2 divides |S3|.

Geometric monodromy is A(kbar) semidirect image(F in Aut(A(kbar))), acting affinely and transitively on A(kbar). The self-cover is strongly absolutely stably atomic exactly when A(kbar) is elementary abelian of prime-power order and the resulting F module is irreducible. This is the primitive-affine criterion UAINFRASELF1, not a supposition that the abelian regular group must lie in the prescribed label group.

When A is not constant, (2.1) is a statement about its genuine finite etale group scheme. For a constant-group formulation impose the necessary splitting of A, for example mu_n in k when phi=[n]. Field-versal torsors yield every associated structured rank-D etale algebra. They are not thereby universal for all rank-D algebras.

The canonical covolume-line functor UAKANLIFT1 turns f into a volume-preserving finite etale self-map of Tot(K_X), with the same degree, complete fibres and monodromy; strong atomicity and field versality persist. It adds one dimension without an arbitrary root choice. No new volume obstruction needs to be solved.

For every nontrivial cover D>1 a nonconstant invariant unit excludes affine space. If F is nontrivial, norms of an induced torus coordinate, followed by the product over its finite F orbit, give a nonconstant invariant monomial unit from B. If F is trivial, B may have rank zero, but S has positive rank because D>1; use a torus coordinate of S instead. The unit descends to X and pulls back to every canonical-line lift and affine identity cylinder. This unit obstruction is not asserted for the degenerate point/isomorphism case.

## 3. UANORMGEN1: the exact norm criterion for every transitive group

Let F act transitively on r>=2 letters, and let

    S=ker[product:Gm^r->Gm],  phi=[n],  n>=2.

Use the base of Section1. Assume char(k) does not divide n|F| and mu_n is in k for the constant G=S[n] semidirect F statement. Then the fixed G torsor from Section2 is field-versal exactly when

    gcd(n,r)=1.                                             (3.1)

For sufficiency, an F torsor Q defines a rank-r etale algebra R/K. Its norm-one torus satisfies

    H1(K,S_Q)=K*/Norm_(R/K)(R*),

of exponent dividing r. The kernel in(2.1) is its n-torsion and vanishes under(3.1). This recovers UAS3VERS1 without requiring F=S3 or the special two-dimensional free base.

For necessity let ell divide gcd(n,r), and let P be a Sylow ell-subgroup of F. By the elementary stabilizer argument in UASYLOWDEF1, P fixes no letter. All its orbit lengths on the r letters are positive powers of ell greater than1. Let d be their greatest common divisor; d is an ell-power and d>=ell.

Choose a generic P Galois extension L/K0, for example L=k(x_g:g in P), K0=L^P with P permuting the independent variables regularly. Induce its P torsor to an F torsor Q. The associated etale algebra is a product of fields L_i/K0 whose degrees are the P orbit lengths. Pass to K=K0(t). At the t-adic valuation each constant extension L_i(t)/K has residue degree [L_i:K0], so the valuations of norms from the product are exactly dZ. Scalar t and Bezout among those degrees also show that t^d is a norm. Consequently

    [t] has order d in K*/Norm(R*),
    [t^(d/ell)] is nonzero of order ell.                     (3.2)

This class is killed by n. Kummer exactness lifts it to H1(K,S_Q[n]), and hence to a G torsor whose S_Q twist has no point. Thus the fixed G torsor is not field-versal. The obstruction works over an extension of every stated base. It concerns specified G torsors, not an asserted one-to-one correspondence between those torsors and their unlabelled associated algebras.

For F=S_r and n=p prime not dividing r, the augmentation F_p module of dimension r-1 is irreducible. A nonzero vector in it cannot be constant since p does not divide r; a transposition difference produces e_i-e_j, and conjugating produces every such difference. Thus these covers are strongly absolutely stably atomic of degree p^(r-1). In particular r odd, p=2 gives every degree 2^(2j), j>=1, over characteristic-zero fields without an extra cyclotomic assumption. These are universality statements for the specified affine groups, not full symmetric groups on p^(r-1) letters.

## 4. UAPRIMEPOWER1: strong atomic fixed self-covers in every odd prime-power degree

Let e>=1 and take F=(C2)^e semidirect S_e, acting on S=Gm^e by coordinate inversions and permutations. Use the free base of Section1. The twist S_Q for an F torsor Q is a restriction of scalars of a quadratic norm-one torus: more explicitly, Q determines an etale algebra R/K of degree e and a quadratic etale algebra E/R, and

    H1(K,S_Q)=R*/Norm_(E/R)(E*),                              (4.1)

with exponent dividing2. This follows from the description of signed permutation lattices and Hilbert90, componentwise if R is disconnected. Hence phi=[n] is field-versal for every odd n. When n=p is an odd prime, the F action on F_p^e is irreducible: a coordinate sign change isolates a nonzero coordinate by division by2, and permutations generate all coordinate vectors. The action is faithful, so monodromy is

    F_p^e semidirect ((C2)^e semidirect S_e).                 (4.2)

Consequently every odd prime-power degree p^e has a fixed strongly absolutely stably atomic self-cover in this construction. Assume char(k) does not divide p|F|, and mu_p in k for the constant-group version. The genuine finite group scheme formulation works without that last assumption. In particular all these constructions work in characteristic zero. A canonical volume-preserving lift costs one dimension by UAKANLIFT1; no minimum dimension is claimed.

For completeness, even n is not field-versal for the full specified group in this architecture. Choose Q induced by the C2 that inverts just one coordinate. Over K=k(s,t), the corresponding quadratic extension K(sqrt(s))/K has norm quotient containing the nonzero order-two class of t, detected by the t-adic valuation. Its Kummer lift obstructs(2.1). This is the exact odd/even torsor boundary, including e=1; it does not forbid an unlabelled degree-two algebra from being realized by a different cover or different auxiliary torsor.

Together with Section3, the positive families now include every odd prime-power degree, and every binary degree 2^(2j). The remaining odd binary exponents are not excluded by these proofs. Section5 supplies the first additional one, degree8, by an integral isogeny that is not multiplication by an integer.

## 5. UACYCIDEAL1: cyclotomic ideals distinguish self-atoms, strong atoms and stable geometry

Let m>=2, O=Z[zeta_m], and let S be the split torus with cocharacter lattice O. Let F=C_m act on S by multiplication by zeta_m. Use the cyclic base B=Res_(k(zeta_m)/k)Gm, with F acting by translation by zeta_m, and put X_m=(B x S)/F. For nonzero a in O the endomorphism multiplication by a defines a self-isogeny phi_a and hence a fixed self-cover f_a of X_m, of degree |Norm(a)|. Work in characteristic zero for the full atomicity and stable-isomorphism assertions below.

For nonunit a:

    f_a is absolutely stably self-atomic
                           <=> a is irreducible in O;        (5.1)
    f_a is strongly absolutely stably atomic
                           <=> (a) is a maximal ideal of O.  (5.2)

Proof of(5.1). Over kbar the B lattice is a sum of trivial F modules, whereas O tensor Q is the nontrivial irreducible cyclotomic rational representation. A dominant map of X_m lifts to its canonical torus cover and intertwines holonomy by some automorphism zeta_m->zeta_m^u. Its B block is an integral matrix, while its S block is necessarily a* sigma_u for some a in O: any rational intertwiner is multiplication by an element of Q(zeta_m) followed by sigma_u, and integrality on O forces that element into O. No off-diagonal block between B and S is possible. Translations affect neither these blocks nor degrees.

In any self-factorization of f_a, the B blocks multiply to the identity and are unimodular. The S blocks therefore express a, up to a unit, as a product a1*sigma_u(a2). A factor is invertible exactly when its a_i is a unit. Thus irreducibility of a excludes every nontrivial self-factorization; a reducible a gives one by composing its multiplication maps. Finite-factor closure and UA1RIGIDCYL1 transfer this argument to every affine identity stabilization, over every extension of constants. This proves absolute stable self-atomicity, rather than merely a restricted factorization within commuting lifts.

For(5.2), the kernel's finite F module is controlled by O/aO, up to the usual dual identification for torus characters. Its invariant additive subgroups are exactly O-submodules, because O is generated by Z and zeta_m. The affine monodromy is primitive precisely when this module is simple, which is equivalent to (a) being maximal. UAFINATOM1 gives the strong statement after every constant extension and identity stabilization. Thus irreducible nonprime elements, when O is not a unique factorization domain, supply further self-atoms that fail strong atomicity. No particular class-number computation is assumed here.

There is an exact geometric version of the ideal obstruction. For any nonzero ideal I of O let S_I have cocharacter lattice I, and X_I=(B x S_I)/F. Then

    X_I is geometrically affine-stably isomorphic to X_O
                         <=> I is principal.                (5.3)

A stabilized isomorphism descends to an isomorphism of the A1-rigid bases by UA1RIGIDCYL1. It lifts to the characteristic torus covers, allowing a holonomy automorphism sigma_u. The trivial B and nontrivial cyclotomic blocks cannot mix. An integral semilinear isomorphism O->I is multiplication by c after sigma_u, so I=cO. Conversely such a principal generator gives an isomorphism. This also proves that a nonprincipal prime ideal gives a genuine degree-Norm(I) cover between different quotient geometries which cannot be turned into a self-cover on X_O by affine identity stabilization. It is a limitation of this fixed cyclotomic lattice architecture, not a universal obstruction for every possible holonomy lattice.

If a generates a prime ideal above p with p not dividing m, the residue field is F_(p^f), f=ord_m(p). The cover is strongly atomic of degree p^f, and is field-versal for ker(phi_a) semidirect C_m by(2.2), since gcd(p^f,m)=1. The constant-group formulation requires mu_p in k; the group scheme formulation does not. Principality of a prime ideal, rather than merely existence of its finite residue field, is the additional self-isogeny gate.

### An explicit degree-eight strong atom

Take m=7 and

    a=zeta_7+zeta_7^2+zeta_7^4,
    a_bar=zeta_7^3+zeta_7^5+zeta_7^6.

Direct reduction modulo Phi_7 gives a+a_bar=-1 and a*a_bar=2. The multiplication matrix has determinant8 and Smith factors1,1,1,2,2,2. Moreover

    O/(a)=F2[z]/(z^3+z+1) congruent F8.                      (5.4)

The polynomial is irreducible and z has order7. Hence geometric monodromy is the full affine group AGL1(F8)=F8 semidirect C7, and the degree-eight map is strongly absolutely stably atomic. Its specified constant-group torsor is field-versal over every k with char(k) different from2,7: mu2 is already in k, and the cyclic base handles zeta7 by restriction of scalars. The dimension of X_7 is 6+[k(zeta7):k], at most12 and equal7 when zeta7 belongs to k. No minimality is asserted.

Invariant top logarithmic forms on B and S descend because the cyclotomic determinant Norm(zeta7)=1. Relative to their product, f_a has constant Jacobian8. Thus its canonical-line volume lift is explicitly

    (x,w)->(f_a(x),w/8)

on X_7 x A1, of dimension at most13, or8 over a field containing zeta7. It retains the entire eight-point fibres, field versality and strong absolute stable atomicity. This is not universality for arbitrary octic algebras, whose full symmetric monodromy is excluded in torus-quotient self-geometry by UAINFRASELF1.

There is also one uniform model over Z[1/14]. Use the degree-six finite etale algebra E=Z[1/14][z]/Phi_7(z) and B=Res_(E/Z[1/14])Gm, rather than the smallest cyclotomic field for an individual k. The integral matrices above, translation by z in B and the free C7 quotient define the same degree-eight construction in every permitted fibre. Its relative source dimension is12, and the volume lift has dimension13. Over every field of characteristic different from2,7 the fibre is field-versal for the constant AGL1(F8) torsor and strongly absolutely stably atomic. Hilbert90 applies to the finite etale algebra E tensor k even when it is a product. This uniform source does not change dimension after extending constants; the dimension7 model when zeta7 belongs to k is a separate smaller choice of B. No all-affine specialization assertion is inferred from this integral model.

## 6. UAVERSLIMIT1: field universality, scheme universality and natural choice are different gates

The following obstruction is an application of a published theorem, not a new general algebraic-group theorem. First--Florence--Rosengarten, arXiv2401.04458v2, Theorem1.1 and Theorems6.3--6.4, prove that an affine algebraic group admits a torsor over a quasi-compact scheme specializing to all torsors over finite-type regular affine schemes exactly when it is unipotent. Over a perfect base, smooth affine tests suffice. If G embeds in GL_s and the receiver base is quasi-projective of dimension d, a missed regular affine test already has dimension at most s(d+1)+2; over a perfect field it can be smooth.

In characteristic zero every nontrivial finite group is non-unipotent. Therefore no fixed torsor over any quasi-compact scheme can realize all of its torsors over all smooth affine bases. In particular a fixed rank-N finite etale cover over a quasi-compact scheme cannot realize all rank-N etale algebras on all smooth affine bases: use its S_N frame torsor. For a quasi-projective receiver of dimension d the standard faithful S_N representation of dimension N-1 gives a missed smooth affine base of dimension at most

    (N-1)(d+1)+2,  N>=2.                                    (6.1)

The published obstruction includes no requirement of a Keller map, a finite ambient parameter space over a prescribed degree, or linear source actions. It applies even to field-versal covers from Sections1--5. The quasi-compact *scheme* hypothesis matters: it does not rule out stack-valued receivers retaining inertia, dimension-dependent approximations or unbounded stable constructions.

There is a useful exact invariance statement. For any torsor E->Y and any p:Y'->Y with a section, the pullback torsor has exactly the same image in H1(U,G) as E, for **every** test scheme U. Projection gives one inclusion and the section gives the reverse. Thus identity stabilization and canonical covolume-line lifting preserve, rather than improve, which scheme families a fixed torsor realizes. For a base change without a section the image can only shrink. The same applies to the frame torsor of a fixed finite etale cover.

Nor can fieldwise existence select points naturally with compatible torsor isomorphisms. A section of the classifying map Y->BG would send the inertia group G of BG to the trivial inertia of the scheme Y, and back to G as the identity, impossible for nontrivial G. Equivalently, for a free G-variety, a natural choice of a point on every torsor twist would already give a G-fixed point on the split twist. Our free torus base has none. This obstruction to point selection is distinct from canonical map-valued descent, which can retain nontrivial automorphisms in its receiver stack; OP-KDESC is not closed by it.

### What the newer literature changes

Scavia, arXiv2607.25118v1, Theorem1.2, supplies a faithful action of (C7 semidirect C3) x C2 on a degree-two Del Pezzo surface for which every Sylow restriction is very versal, while the full action is not versal. Therefore successful Sylow-versality checks cannot certify a general positive versality claim. The use of a Sylow subgroup in Section3 is a direct construction of a negative norm obstruction; it does not use the disproved converse. UASYLOWDEF1's fixed-fibre degree bounds are likewise unaffected.

Kresch--Tschinkel, arXiv2606.10965v1, Theorem2.5, gives a valid alternative in their algebraically closed characteristic-zero setting: the listed properties of all twists, including rational points, can be tested on one generic torsor arising from a genuinely versal generically free action. That is a geometric generic-twist test, not a bounded sample or an automatic finite-field principle. Duncan--Reichstein's Theorem1.1 distinguishes point existence, density, unirationality and stable linearizability. We use the point-existence assertion precisely; none of the new constructions asserts functorial choices or scheme universality from a finite list of field examples.

No source used here is being described as AI-assisted without an author statement establishing that fact.

## 7. Remaining requirements

The general construction removes the need to invent a special free torus base for each group. The exact isogeny cohomology test and uniform annihilators replace degree-by-degree experiments. Norm-one tori have a sharp all-group coprimality boundary, and cyclotomic principal ideals add strong primitive self-covers beyond integer power maps. These constructions remain outside affine space and are field-versal for specified groups rather than full symmetric universality in arbitrary rank.

The original higher-degree nonfinite affine canonical boundary, effective simultaneous fixed-locus landing and gluing, smaller marked-rank ambient atomicity, modular cases and A1-rich pseudo-plane mixing remain open. Other odd binary exponents need suitable integral isogenies or different holonomy geometry; nonprincipal ideal obstructions must not be confused with a proof of impossibility in every architecture. The all-affine scheme obstruction explains why field-versal covers alone do not supply the sought global canonical receiver.
