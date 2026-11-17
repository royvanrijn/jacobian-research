# Ruling descent, product-leaf topology and sheet-loss budgets

15 November 2026. Three written general claims with selected exact controls. Independent whole-proof replay, formal verification, external review and novelty remain unverified. OP-KDESC, JC(2), nonlinear integration and canonical absolute stable atomicity remain open.

The main extension is broader than the explicit Danielewski family: ANY affine-line ruling defined over the retained field on a proper native intermediate forces a geometric affine-plane atlas. Its target restrictions would disprove JC(2). Geometric ruling existence alone does not supply descent to that retained field.

Use the actual UAINTERSURF1 diagram of a whole retained-signal factorization F=E composed G in characteristic zero:

    U --g--> V --e--> A2, K0=k(u,a),
    U=Spec K0[S,Y,r^(-1)], r=1-SY.

For G noninvertible, UAUNITFACT1 gives constant GEOMETRIC units on V and an arithmetic UFD A=O(V). It also has a specified K0-point, carried from (u,0,0,a). UALOGVOL1 gives negative geometric log dimension. For native F, UAGEOPIC1 gives nonzero geometric Picard group and UAPLANECHART1 gives the whole carried curve Gamma=A1 with embedded target image. These are actual whole-factor inputs, not conditions to be inferred from an abstract surface model.

## 1. UAARITHRULING1: a descended ruling produces a plane atlas

Suppose V, in the noninvertible G branch, admits a surjective morphism

    phi:V->C

DEFINED OVER K0, where C is a smooth geometrically integral curve and the generic fibre is A1 over K0(C). An affine-line form in characteristic zero is also split: its smooth projective completion has genus zero and its unique geometric boundary point is rational, giving P1 minus a rational point. Thus the usual geometric-generic A1 hypothesis gives the same generic fibre. We do not assume that a geometric ruling descends.

### The base is A1, including the exclusion of complete-type rulings

First, pullback on Picard groups is injective:

    phi*:Pic(C)->Pic(V) is injective.                 (1.1)

Represent a class by a divisor D. If phi*D=div(h), restrict h to the generic affine-line fibre. There it has no zeros or poles, hence is a unit of K0(C)[w] and belongs to K0(C)*. Thus h is a base rational function. Surjectivity of phi makes pullback of base divisors injective (every closed-point fibre is a nonempty divisor), so D=div_C(h). This proves (1.1).

Arithmetic factoriality gives Pic(V)=0, hence Pic(C)=0. The dominant rational map U->C makes the smooth projective completion of C geometrically rational. Its K0-point, the image of the specified point of V, identifies that completion with P1_K0. One can justify genus zero directly by restricting a dominant rational map from the rational plane to a general rational curve and using Riemann-Hurwitz.

C cannot be the complete P1, whose Picard group is Z. If the geometric boundary of C in P1 had at least two points, C would have a nonconstant geometric unit, pulling back to one on V. Consequently the boundary has exactly one geometric point; in characteristic zero it is a K0-rational point. Therefore

    C is A1_K0.                                      (1.2)

This excludes descended complete-type rulings. It says nothing about a complete-type ruling existing only after a scalar extension.

### Every arithmetic closed fibre is integral; every geometric fibre is reduced

Write phi=f:V->A1 and K0(T)=K0(f). Surjectivity and normality give

    A intersection K0(f)=K0[f].                       (1.3)

Indeed, a rational function of T with a pole at a closed base point retains a pole on at least one component of its nonempty fibre.

Let p(T) be ANY irreducible polynomial in K0[T]. Factor p(f) in the UFD A:

    p(f)=c*product h_j^(n_j), c in K0*.

Each prime h_j is vertical. On the generic A1 fibre it is a unit, so h_j belongs to K0(f)*. Equation (1.3) then puts every h_j in K0[f]. The displayed equality is therefore a factorization of the irreducible p(T) inside K0[T]; it has one factor with exponent one. Thus

    A/(p(f)) is a DOMAIN for every closed base point.  (1.4)

Over its residue field this fibre is integral. That field has characteristic zero, so the fibre is geometrically reduced. Geometric irreducibility is NOT asserted: the arithmetic prime may split into several geometric components. Those components form one orbit under the Galois group of the residue field.

This is a general factorial-fibration argument, not a newly claimed classification theorem. It is stronger here than merely excluding a common multiplicity. In particular no geometric fibre component of a descended ruling has multiplicity greater than one.

### A geometric plane atlas and the precise JC(2) consequence

We explicitly import two classical surface facts, restated in Dubouloz2020, Definition4 and Proposition3: degenerate fibres of a smooth affine-line-fibred surface are disjoint unions of affine lines with multiplicities, and an affine-line fibration with all fibres irreducible and reduced is a Zariski locally trivial A1-bundle.

Apply them after algebraic closure. Our fibres are reduced, so every component is a disjoint A1. There are finitely many reducible fibres. In EACH choose one component and delete all other components. The remaining open W is affine: the deleted union is an effective Cartier divisor on the smooth affine surface, and its complement is affine locally on the original affine target, hence affine. Its morphism to A1 has all fibres A1, irreducible and reduced, so it is an A1-bundle. Over affine A1 every such bundle is trivial (Pic(A1)=0 and H1(A1,O)=0). Hence

    W is A2.                                         (1.5)

Varying the selected components gives a FINITE geometric A2 atlas covering V. If there were no geometrically reducible fibre, V itself would be A2.

For native F, a proper factor has noninvertible e. UAPLANECHART1 now makes its restriction to EVERY chart (1.5) a nonautomorphic plane Keller map. Therefore

    proper native factor + a ruling over K0
          implies a JC(2) counterexample.             (1.6)

This is an IMPLICATION, not a constructed counterexample or a proof of JC(2). Assuming JC(2), every affine-line ruling of a proper native intermediate must fail to descend to K0.

There is also an unconditional curve requirement. If the whole Gamma meets at most one component in EACH reducible fibre, select those components in (1.5). The resulting plane chart contains all of Gamma, contradicting UAPLANECHART1. Thus at least one fibre contains TWO different components hit by Gamma. The ruling restricted to Gamma cannot be constant or degree one; its polynomial degree is at least two. This generalizes the earlier single-bad-fibre Danielewski curve test to every descended ruling, with no bound on fibre count, component count or polynomial degrees.

### Why the descent hypothesis cannot be dropped

UALOGVOL1 supplies geometric negative log dimension, and hence geometric affine-line rulings via the standard surface theorem. Dubouloz--Kishimoto2014 explicitly distinguishes ruling existence after a FINITE FIELD EXTENSION from existence over the original field. We do not identify them. Arithmetic factoriality can fail after that extension, so (1.4) cannot silently be repeated there.

This locates a concrete remaining gap: either realize the compatible JC(2)-counterexample atlas forced by a descended ruling, or control a genuine failure of ruling descent together with its full whole-polynomial data. It is not enough to cite geometric ruledness, a rational point or arithmetic UFD as a descent proof.

### A Galois-invariant ruling already descends

There is an important positive descent case in the SAME published source. Dubouloz--Kishimoto2014, Section1.3.1 and Proposition4, proves that a ruling invariant under a finite Galois extension's action, up to base-curve automorphisms, descends to the original field. The base automorphisms are unique because the ruling is surjective, so they satisfy the cocycle identity; the quasi-projective base curve and affine morphism descend effectively. The generic A1 form is trivial in characteristic zero.

In particular, a UNIQUE geometric ruling up to base automorphism descends after choosing a finite Galois field of definition. Applying (1.6),

    a proper native factor with a Galois-invariant ruling,
    in particular a unique geometric ruling,
        would disprove JC(2).                         (1.7)

Consequently, CONDITIONAL on JC(2), every geometric ruling of a proper native intermediate must belong to a NONTRIVIAL Galois orbit of inequivalent rulings; at least two rulings with different general fibres are required. This is a published descent criterion applied to our arithmetic-fibre theorem, not a general proof of ruling descent. Neither the existence of such an orbit nor its compatibility with the whole factor is constructed.

## 2. UAPRODUCTTOPO1: finite stages and fundamental groups of product leaves

This part needs neither a native curve nor two residual dimensions. Suppose a whole Keller map has a retained target tuple p whose GEOMETRIC generic source fibre is

    X_(N,h)=U^N times A^h, N>=1, h>=0.                 (2.1)

The actual generic intermediate diagram exists by base change for every whole factor, with arbitrary coordinate and identity mixing. No assertion is made that an arbitrary changed-signal canonical receiver satisfies (2.1).

### The source has Euler characteristic one and abelian fundamental group

Use complex models. Compactly supported Euler characteristic gives

    chi_c(U)=chi_c(C^2)-chi_c({SY=1})=1-0=1,
    chi_c(X_(N,h))=1.

For pi1(U), cover U by A={S!=0} and B={abs(SY)<epsilon}, 0<epsilon<1. A is C* times C* in coordinates (S,r). B is contractible by first shrinking S to0 and then Y to0. Their intersection is C* times a disk, and its fundamental group maps to the S-loop in pi1(A)=Z^2. Van Kampen therefore gives

    pi1(U)=Z, generated by the r-loop,
    pi1(X_(N,h))=Z^N.                                 (2.2)

These are classical topological calculations; Hatcher Theorem1.20 is the imported van Kampen theorem.

If g:X_(N,h)->V is FINITE etale of degree d, finite-cover multiplicativity gives 1=d*chi_c(V), forcing d=1. The finite birational map to normal V is an isomorphism. For an actual whole factor, UAALLATOM1 then makes G an automorphism. Conversely an invertible G gives an invertible, hence finite, g. Similarly a finite e:V->A^(2N+h) is trivial in characteristic zero, since complex affine space is simply connected. Thus

    G invertible iff its generic g is finite;
    E invertible iff its generic e is finite.         (2.3)

Consequently BOTH generic arrows of a proper product-leaf factor are nonfinite. The proof works geometrically in every characteristic-zero field by descending the finite morphism data and degree to finitely many coefficients and using a complex embedding. It is not a positive-characteristic finite-cover theorem.

### A finite-index abelian subgroup whose index DIVIDES the degree

Let g:X_(N,h)->V be dominant and generically finite of degree d, with V smooth and geometrically integral. Over a dense open V0 it is a connected finite etale cover; its source X0 is open in X_(N,h). Deleting complex subvarieties of positive codimension induces SURJECTIONS on pi1, since loops can be perturbed off their real-codimension-at-least-two strata.

Write B=pi1(V0), K=image pi1(X0), and Q=pi1(V). Then [B:K]=d. The image H=g_*pi1(X_(N,h)) is the image of K under the surjection B->Q. The inverse image of H contains K, so the index formula gives

    [Q:H] DIVIDES d.                                  (2.4)

By (2.2), H is abelian and generated by at most N elements. Hence every actual intermediate has VIRTUALLY ABELIAN geometric fundamental group of rank at most N, with an explicit abelian subgroup of index dividing the whole source-prefix degree. For one block it is VIRTUALLY CYCLIC. The same finite-cover consequence holds for geometric etale pi1 by comparison; no topological pi1 is assigned to an arbitrary noncomplex field.

Finite cover theory and the subgroup index formula are classical (Hatcher Proposition1.32). Cofinite image does not make g finite, and (2.4) does not make g injective on fundamental groups.

### An exact sheet-loss cocycle

For any dominant quasi-finite etale f:X->Y of degree d define

    n_f(y)=number of geometric points in f^(-1)(y),
    h_f=d-n_f, D(f)=integral h_f d chi_c.

Etaleness implies 0<=n_f<=d, and n_f is constructible. Finite algebraic stratification makes the restrictions finite covers, so additivity and finite-cover multiplicativity prove

    D(f)=d*chi_c(Y)-chi_c(X).                           (2.5)

For f:X->Y of degree a and j:Y->Z of degree b, actual fibre counts give

    h_(j composed f)=a*h_j+j_!(h_f),
    D(j composed f)=a*D(j)+D(f).                       (2.6)

Thus a product-leaf diagram with target affine space has

    D(g)=a*chi_c(V)-1, D(e)=b-chi_c(V),
    D(e composed g)=a*b-1.                            (2.7)

This is standard constructible Euler calculus, explicitly derived here; no novelty is claimed for Fubini or the cocycle identity. Although h_f is pointwise NONNEGATIVE, its Euler integral need not be: a positive constant on a curve of negative Euler characteristic has negative integral. Do not infer separate positive losses at both stages from (2.7).

In dimension two, nonfiniteness also forces a missing DIVISOR in the normal finite completion of the source. Zariski Main gives X open in a finite normal affine Z over Y. If the boundary had only finitely many points, normal Hartogs extension would identify their coordinate rings, and affineness would force X=Z. Thus a proper surface stage cannot be repaired merely by finitely many added source points. This boundary is different from the at-most-finitely-many points missing from its TARGET IMAGE.

The old all-degree J_d(S,Y)=(S,Y*(1+r+...+r^(d-1))) controls are sharp: they have d sheets for S!=0 and ONE sheet on the entire S=0 line. Therefore h_(J_d)=(d-1)*1_(S=0) and D(J_d)=d-1. Their nonconstant unit index and whole-lift poles still exclude them as whole Keller factors by UAUNITFACT1; the Euler identity alone does not construct one.

## 3. UAORBFIB1: geometric multiple-fibre patterns and an even-degree gate

This statement applies to every V reached from (2.1) by the degree-d map g, even when no ruling descends. Let

    phi:V->A1

be a surjective morphism with connected geometric generic fibre. For each multiple fibre let m_i>=2 be the gcd of its component multiplicities. We explicitly import the classical orbifold fundamental-group quotient, restated in Cogolludo-Agustin--Elduque, arXiv:2306.07359v3, Definition2.8 and Remark2.11:

    pi1(V) surjects onto C_(m1) * ... * C_(ms).        (3.1)

The cited statement is about connected generic fibres and the maximal orbifold structure; it is NOT a classification of our surfaces or a theorem about every component having multiplicity m_i. Surjectivity of phi is an actual hypothesis: omitted values add free generators and are not silently counted as multiple fibres.

By (2.4) every quotient of pi1(V) is virtually abelian. A free product of s nontrivial finite cyclic groups is virtually abelian only if

    s<=1, OR s=2 and m1=m2=2.                         (3.2)

For completeness, map the free product to Q=product C_(mi). Its kernel acts freely on the usual free-product tree; the finite quotient graph has |Q| central vertices, sum |Q|/m_i other vertices, and s|Q| edges. The kernel is free of rank

    1+|Q|*(s-1-sum 1/m_i).

Except for a single finite factor or two order-two factors, that rank is at least2. A virtually abelian group cannot contain such a free group. This is elementary free-product/graph theory, not a bounded census.

In the two-double case the quotient is infinite dihedral. The image of the finite-index abelian H of (2.4) is an infinite abelian finite-index subgroup of that quotient. It contains only rotations: an abelian subgroup containing a reflection is finite. Its index is therefore even. This index divides [pi1(V):H], which divides d. Hence

    two double fibres force 2 divides deg(g).         (3.3)

In particular an ODD whole source-prefix degree excludes every geometric connected-generic fibration over A1 with two multiple fibres. More than two multiple fibres, or two of any other multiplicities, are excluded in ALL degrees and all product-leaf ranks.

For the one-block proper class, UAARITHRULING1 is stronger when an A1 ruling DOES descend: every geometric component is then reduced, so there are no multiple fibres at all. Equations (3.2)-(3.3) retain content for rulings which exist only after a field extension, where arithmetic factoriality is no longer available.

## 4. Remaining research gate and assurance

The main geometric alternatives are now explicit. A proper native factor with a ruling over K0 would produce a geometric plane atlas of actual JC(2) counterexamples; its whole carried line must enter two components of some split fibre. A candidate escaping that implication needs a PROVED failure of ruling descent, while retaining the negative-log, unit, Picard, exact-volume, cofinite-arrow and whole-source constraints.

Product-leaf topology extends a finite-stage and fundamental-group obstruction beyond two residual dimensions, but only under the stated actual product generic-fibre hypothesis. It does not establish that hypothesis for the mutually changing blocks of the canonical trace receiver, or replace their separate marked-scheme and gluing requirements.

Both proper arrows are genuinely nonfinite; missing source divisors and the sheet-loss identities must be realized, not replaced by finite covers or formal neighborhoods. The orbifold gate excludes whole multiplicity patterns and the odd-degree two-double case. It does not construct a surviving intermediate.

All imported topology, surface-fibration and orbifold statements are identified in the literature note. The deductions require independent proof review and a fuller priority audit. Complete canonical Spec(A), actual nonlinear polynomial integration, effective global gluing, improved costs, JC(2) and absolute stable atomicity remain open.
