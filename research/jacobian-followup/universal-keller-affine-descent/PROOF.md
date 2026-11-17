# Descent of affine-type rulings on actual Keller intermediates

15 November 2026. Two written general results. The new descent statement is unconditional in its stated proper-native branch. Its exclusion of that branch is conditional on JC(2). Complete-type Galois-permuted rulings, OP-KDESC and absolute stable atomicity remain open. Selected symbolic controls are not independent whole-proof replay, formal verification, external review or novelty evidence.

Use the actual retained-signal diagram of UAINTERSURF1:

    F=E composed G, K0=k(u,a), characteristic zero,
    U=Spec K0[S,Y,(1-SY)^(-1)] --g--> V --e--> A2_K0.

Let A=O(V), B=A tensor_K0 K0bar. V is smooth, affine and geometrically integral, and its canonical bundle is trivial: omega=e*(db wedge dc) is nowhere zero and exact. The whole section supplies P in V(K0). If G is not invertible, UAUNITFACT1 gives constant geometric units and arithmetic factoriality. A proper native factor means BOTH whole G and E are noninvertible; UARELPLANE1 then excludes V_geometric=A2. The same conclusions apply to another retained-signal map with its required embedded-curve witness. Identity stabilization and coordinate mixing are unrestricted.

The two results below use these established hypotheses. They do not add a new candidate-surface checklist.

## 1. UAMLFACTOR1: trivial geometric ML is exactly the target-automorphism branch

The Makar-Limanov invariant ML(B) is the intersection of the kernels of all locally nilpotent K0bar-derivations of B; if there are no nonzero such derivations it is B.

For EVERY actual retained-signal diagram, even without a native witness,

    ML(B)=K0bar  implies  V_geometric=A2.             (1.1)

For the native/witness class this gives the exact equivalence

    E is an automorphism  iff  ML(B)=K0bar.           (1.2)

We explicitly import Daigle, Affine surfaces with trivial Makar-Limanov invariant, Theorem3.2, arXiv:0705.0334v1. For a smooth affine surface with ML equal to the ground field, trivial canonical bundle is equivalent to a presentation

    B=K0bar[u,v,z]/(u*v-P(z)), P nonconstant.         (1.3)

This is a published classification theorem, not a classification established by the selected controls. Smoothness makes P square-free. If deg(P)=1, (1.3) is the affine plane. If deg(P)>=2, UADANVOL1 proves that its canonical volume has a nonzero period and that no etale map from this surface to A2 exists. The actual e is such a map, so this case is impossible. This proves (1.1).

For native maps UARELPLANE1 turns the resulting geometric plane into whole E invertibility. Conversely whole E invertibility makes V the plane, whose two translation derivations have kernel intersection equal to the constants. This proves (1.2).

In particular every proper native intermediate has NONCONSTANT geometric ML. We do not classify surfaces with nonconstant ML, or infer that two complete-type rulings make ML trivial.

## 2. UAAFFDESC1: every affine-type geometric ruling actually descends

For a proper native intermediate, ANY geometric A1-ruling over an AFFINE smooth curve is unique up to a base automorphism, descends to K0, and admits a nonzero locally nilpotent derivation defined over K0.

The descent is constructive at the level of its coordinate:

    choose f over a finite Galois extension L/K0;
    f0=f-f(P);
    choose c in L with h=sum_sigma sigma(c*f0) nonzero;
    then h belongs to A and gives the descended A1-ruling. (2.1)

There is no assumption that a chosen geometric ruling was already Galois invariant. We prove that invariance from (1.1).

### 2.1 Affine bases and a volume-defined additive action

First work geometrically and let phi:V_geometric->C be a surjective A1-ruling with C affine. The dominant rational map from U to C makes the projective completion of C rational. Constant geometric units exclude two or more boundary points; affineness supplies at least one. Thus C=A1. Write phi=f.

For a surjective ruling on a normal surface,

    B intersection K0bar(f)=K0bar[f].                 (2.2)

A base rational function with a finite pole has a pole on a component of that nonempty fibre. This proves (2.2).

Define D_f by contraction with the actual nowhere-zero volume:

    i_(D_f) omega = df.                              (2.3)

The volume identifies the tangent and cotangent bundles, so D_f is a global regular derivation. It kills f. On the generic fibre B tensor_K0bar[f] K0bar(f)=K0bar(f)[w], the volume has the form

    omega=q(f)*df wedge dw, q(f) in K0bar(f)*.

Indeed its coefficient is a nowhere-zero regular function on that affine line and hence is a unit in K0bar(f)[w]. Consequently

    D_f=-q(f)^(-1)*partial_w

on this localization. It is locally nilpotent there and therefore on B, since every element of B has finite w-degree and D_f preserves B. Its kernel is precisely K0bar[f] by (2.2).

This use of the A1 generic fibre is essential. A regular etale lift of a plane translation is not automatically locally nilpotent.

### 2.2 Nonconstant ML forces uniqueness among AFFINE-type rulings

Suppose two affine-type rulings f1,f2 have inequivalent general fibres. Each supplies the locally nilpotent D_fi in (2.3). Let Fi=K0bar(fi).

The fields Fi are relatively algebraically closed in Frac(B), because Frac(B)=Fi(wi). If the two kernels had a nonconstant common element h, both Fi would be algebraic over K0bar(h). Relative algebraic closedness would force F1=F2. Equation (2.2) would then make K0bar[f1]=K0bar[f2], so f2=alpha*f1+beta and the rulings would be equivalent.

Thus two inequivalent affine-type rulings have kernel intersection equal to K0bar. ML(B), contained in this intersection, would be K0bar. UAMLFACTOR1 excludes that possibility for a proper native factor. Hence there is exactly ONE affine-type ruling class whenever one exists.

Galois conjugation preserves affine type. Therefore every conjugate of our ruling belongs to that SAME class. This closes the Galois-permuted AFFINE-type escape, rather than presuming it absent.

### 2.3 An explicit descended coordinate

Choose a finite Galois extension L/K0 over which f, its A1-fibration and a coordinate on its affine base are defined. For each sigma in Gal(L/K0), uniqueness gives

    sigma(f)=a_sigma*f+b_sigma, a_sigma in L*, b_sigma in L.

The coefficients are in L: both sides belong to A tensor L and 1,f are linearly independent over L. Center f at the specified rational point P. Since sigma(P)=P,

    f0=f-f(P), sigma(f0)=a_sigma*f0.                 (2.4)

The semilinear one-dimensional space L*f0 is therefore stable. For c in L form the weighted trace in (2.1). It is Galois fixed and lies in A, and

    h=(sum_sigma a_sigma*sigma(c))*f0.

Some c makes the coefficient nonzero: otherwise the distinct field automorphisms sigma would satisfy a nontrivial L-linear relation, contradicting their linear independence. This is the usual one-dimensional Galois descent/Hilbert90 argument, written here as an explicit trace formula.

Over L, h is a nonzero scalar multiple of f0, so its ruling is exactly the original one after changing the base coordinate. Surjectivity descends faithfully flatly. The generic fibre is a characteristic-zero form of A1 and hence is A1; equivalently use the published ruling-descent criterion of Dubouloz--Kishimoto2014, Section1.3.1. Thus h:V->A1_K0 is the requested ruling over the RETAINED field.

Finally define D_h over K0 by i_(D_h)omega=dh. After scalar extension its generic-fibre calculation is (2.3), so D_h is nonzero and locally nilpotent on A. This constructs an additive action over K0 as well as the descended ruling.

No arithmetic UFD hypothesis is transported to L. Factoriality is used only AFTER descent, in UAARITHRULING1.

### 2.4 What this closes and what it does not

UAARITHRULING1 now applies to h. Every arithmetic closed fibre is integral, every geometric fibre is reduced, and there is a finite geometric plane atlas. For a proper native factor every target restriction to a chart would be a JC(2) counterexample. Therefore, unconditionally,

    proper native factor + ANY geometric affine-type A1-ruling
        implies a JC(2) counterexample.              (2.5)

The standard quotient theorem for additive actions on normal affine surfaces identifies every nontrivial additive action with an affine-type A1-ruling over a normal affine curve. This fact is explicitly recalled by Dubouloz2015, introduction, arXiv:1507.05802v1. It and (2.1)-(2.3) show that, for a proper native intermediate, the following existence statements are equivalent:

    geometric affine-type A1-ruling;
    retained-field A1-ruling of any type;
    nonzero geometric locally nilpotent derivation;
    nonzero retained-field locally nilpotent derivation. (2.6)

For the second implication UAARITHRULING1 makes the descended base affine A1. For the third-to-first implication we import the standard surface quotient theorem just stated.

Assuming JC(2), every proper native intermediate must consequently be GEOMETRICALLY ADDITIVELY RIGID, ML(B)=B. Its geometric A1-rulings supplied by UALOGVOL1 must ALL have complete base P1. No complete-type ruling can be Galois invariant: its published descent would give a retained-field ruling, which UAARITHRULING1 proves cannot have complete base. This last complete-type non-invariance is unconditional for the actual noninvertible-source branch.

Thus the remaining descent question concerns COMPLETE-type rulings only, if JC(2) holds. The affine-type descent problem has been solved in this actual diagram. This is not a proof of descent for the remaining complete-type orbit, nor a solution of JC(2) or stable atomicity.

## 3. Why averaging distinct rulings does not solve the remaining descent

The weighted trace above works because uniqueness first puts all conjugates in ONE one-dimensional space. Averaging arbitrary ruling coordinates can destroy the generic A1 fibre even after centering at a rational point.

Here is an exact control over Q(sqrt(2)), with d=2 and s^2=d:

    f_plus=x+s*(y+s*x^2)^2,
    f_minus=x-s*(y-s*x^2)^2=sigma(f_plus).

Each is a polynomial coordinate: (x,y)->(x,y+s*x^2)->(x+s*(y+s*x^2)^2,y+s*x^2) is a composition of triangular automorphisms, and conjugation gives the other one. Both vanish at P=(0,0). But

    (f_plus+f_minus)/2=x+2*d*x^2*y=H.

On the generic H-fibre, x is invertible because

    H=x*(1+2*d*x*y);
    y=(H-x)/(2*d*x^2).

Hence its coordinate ring is Q(H)[x,x^(-1)], so its generic fibre is Gm, not A1. This is an exact symbolic counterexample to the proposed trace shortcut, not a bounded miss. The plane is not claimed to be a proper intermediate; it is a control for the averaging operation itself.

For a complete-type ruling V_geometric->P1 there is not even a global affine coordinate f to center. Pullback of O(1) is nontrivial by the generic-A1 Picard injection in UAARITHRULING1. Taking the norm of its finite Galois orbit gives an invariant Picard class, which UAGEOPIC1 makes zero, but this only trivializes a tensor product. It does NOT identify the conjugate ruling classes or supply an A1-fibration. The two P1 projections on (P1 times P1) minus the diagonal already have opposite nonzero Picard classes; the ordinary affine quadric realizes this geometric pattern and UADANVOL1 excludes it as an actual target. Thus a zero orbit sum is not a ruling-descent proof.

Published rigid smooth affine cubic surfaces with complete-type rulings (Dubouloz2013 and2015) also prevent replacing “two AFFINE-type rulings” by “two arbitrary rulings” in Section2. These published surfaces are not asserted to satisfy our whole-factor hypotheses or to be surviving intermediates.

The next descent work must identify conjugate COMPLETE-type fibrations, or exclude their actual whole realization using a new argument. Adding further necessary conditions or assuming that a trace retains A1 fibres does not do that.
