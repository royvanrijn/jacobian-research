# Joint selectors, Bertini and the auxiliary first-integral ring

15 November 2026. Written general proofs with selected exact controls. No independent whole-proof replay, formal verification, external review or publication-priority certificate. OP-KDESC remains open. This extends the [source-only selector theorem](../universal-keller-selector-rigidity/PROOF.md) to joint polynomial selectors beta(P,U).

Work over a characteristic-zero field k. Use N>=1 retained three-coordinate compiler blocks P_i, any finite higher degrees j>=4, and arbitrary polynomial profiles g_ij(U,sigma_1,...,sigma_N). The first output sigma_i is independent of all coefficients. All maps and factors are polynomial on the entire stated affine space. Source-dependent profiles beyond these retained signals, higher-rank return and positive characteristic remain outside.

## 1. The lifted derivation and all constant Jacobians

Let Phi(P,U) denote the full profiled compiler. UAPARALLEL1 gives determinant one when U is retained as identity coordinates. Let G:A^m->A^m be ANY polynomial determinant-one auxiliary map, m>=2, and let mu be a nonzero constant vector. Put

    Delta=(DG)^(-1)*mu,
    A=D_P Phi,
    W=A^-1*(D_U Phi)*Delta,
    D=Delta-W on k[P,U].                                (1.1)

These are polynomial vector fields, since both inverse Jacobian matrices are adjugates. The fields Delta and W act on U and P respectively. We have D(Phi)=0 and D(G)=mu. For ANY joint polynomial selector beta(P,U), form

    T_beta(P,U)=(Phi(P,U),G(U)-beta(P,U)*mu).

The exact Schur identity is

    det DT_beta=1-D(beta).                              (1.2)

Thus a Keller member with determinant c in k^times satisfies D(beta)=1-c. Unlike the source-only case, c need not equal one.

Choose a linear output change M in SL_m(k) with M*mu=e_m. Write M*G=(R,t), where R has m-1 coordinates and t is the last one. Then

    Delta R=0, Delta t=1, D t=1.

For c nonzero, (1.2) is equivalent to

    beta=(1-c)*t+h,    h in ker D.                      (1.3)

No local nilpotence is inferred from this polynomial slice.

Write E_aux=k(U), F_aux=k(R), L=ker(Delta on E_aux), and C=ker(Delta on k[U]). UAGEOMFEED1 identifies L as the relative algebraic closure of F_aux in E_aux: their transcendence-degree difference is one, and Delta t=1. Moreover

    Frac(C)=L.                                         (1.4)

To see (1.4), an element a in L is algebraic over k(R). Clearing the denominators of its monic minimal polynomial makes q(R)*a integral over k[R] for some nonzero q. The normal ring k[U] then contains q*a. It is a Delta-constant, and q is too. This does not assert C is finitely generated, polynomial, or integral over k[R].

Call a compiler block active if some Delta g_ij is nonzero; let I be that set. Process just these blocks, retaining the others and all U coordinates, as in UAPARALLEL1. Denote the resulting 3N compiler coordinates by

    Z=(Phi_i for i in I, P_l for l outside I).

The map H_I:(P,U)->(Z,U) is a polynomial Keller self-map. Its source variables Z,U are algebraically independent. All original signals are polynomially recoverable in Z. Since W preserves every signal and its inactive block components vanish, D fixes Z and acts on U as Delta.

## 2. The precise Bertini field step

We use the classical generic-field Bertini lemma, Stacks Lemma37.32.1. Its application can be stated as follows. If K/F is a geometrically integral finitely generated field extension and u,v in K are algebraically independent over F, then for an independent variable b,

    K(b) is geometrically integral over F(b,v+b*u).      (2.1)

Here is the hypothesis check from the cited lemma. Use two independent variables x_1,x_2 and g_1=u,g_2=v, taking the constant combination delta=u. The four elements x_1,x_2,x_1*u+x_2*v,u are algebraically independent, since x_2 is nonzero and u,v are independent. The lemma applies. Set b=x_1/x_2 and divide the universal linear combination by x_2; removing the extra pure transcendental variable x_2 gives (2.1). Characteristic zero supplies geometric reducedness. This is a generic-field statement, not an assertion about every numerical line.

For a geometrically integral curve field K/F, let S be independent and

    f(S)=S+sigma*S^3+sum_j a_j*S^j,

with a_j in K and at least one coefficient nonconstant along a nonzero derivation of K fixing F. Then S^2 and f are algebraically independent over F in K(S). Indeed the derivation extended to fix S annihilates F(S) and does not annihilate f. In characteristic zero f cannot be algebraic over F(S). Apply (2.1) to the surface field K(S), with u=S^2 and v=f. It follows that the generic inverse equation

    S+b*S^2+sigma*S^3+sum_j a_j*S^j-c/2=0               (2.2)

defines a geometrically integral curve over F(b,c). Both independent target parameters b,c are essential to this application. The argument requires no bound on coefficient degrees or pole orders on the auxiliary curve.

## 3. UAJOINTKER1: the full joint first-integral ring

For the general auxiliary G and active set I above, the exact rational and polynomial kernels are

    ker(D on k(P,U))=L(Z),                              (3.1)
    ker(D on k[P,U])=C[Z].                             (3.2)

The latter is a polynomial extension in 3N independent coordinates over the FULL auxiliary constant ring C, not a claim that C itself is a polynomial ring. No output-recovery or finite-generation hypothesis on C is needed.

### Geometric integrality after each active block

Begin with

    F0=L(sigma_i:i in I, P_l:l outside I),
    K0=k(U)(sigma_i:i in I, P_l:l outside I).

The auxiliary curve field k(U)/L is geometrically integral by relative algebraic closedness and characteristic zero, so K0/F0 is too. The horizontal derivation restricts to K0, fixes F0 and has D t=1.

Adjoin active blocks one at a time. The retained inverse coordinates from UAM1 are S_i=x_i/(1+x_i*y_i) and r_i=1/(1+x_i*y_i). If b_i is the second compiler output, then

    r_i=1+2*b_i*S_i+3*sigma_i*S_i^2
                     +sum_j j*g_ij*sigma_i^j*S_i^(j-1).

Thus (S_i,b_i) are independent rational coordinates for the two source variables over the already retained sigma_i. Its third output c_i obeys

    c_i/2=S_i+b_i*S_i^2+sigma_i*S_i^3
                        +sum_j g_ij*sigma_i^j*S_i^j.    (3.3)

On the preceding curve field the horizontal derivation acts on every g_ij(U,sigma) as Delta. Activity says the derivative of the last coefficient sum is a nonzero polynomial in S_i. Extend this preceding derivation temporarily to fix the new independent S_i; it proves the algebraic independence required in Section2. The Bertini step now shows the enlarged source field is geometrically integral over the enlarged target field obtained by adjoining the independent pair b_i,c_i.

The actual horizontal derivation fixes these new outputs, with D S_i=-sum_j(Delta g_ij)*sigma_i^j*S_i^j/r_i, and still has D t=1. It restricts to each field in this induction; profiles in other blocks depend on signals and U, not on unadjoined leaf variables. Repeating Section2 handles every active block, even when different coefficient profiles depend on each other's retained signals.

The resulting source field is the original k(P,U), by the inverse reconstruction, and it is geometrically integral over L(Z), with transcendence-degree difference one. A nonzero derivation of a one-variable function field in characteristic zero has precisely the relative algebraic closure of the base as constants. Therefore its constants are L(Z). This proves (3.1). Inactive blocks must be left in Z as full source coordinates; treating their inverse pencils as active would invalidate the algebraic-independence hypothesis.

### Polynomial saturation of the intermediate affine space

An elementary regularity fact applies to H_I. For any dominant polynomial Keller map H between affine spaces, a rational function on its target which pulls back to a polynomial is itself polynomial. To prove this, take a nonconstant irreducible reduced denominator q. Dominance and constant affine units make q composed with H a nonconstant nonunit. A height-one prime above it contracts to (q) by flat going-down, since H is etale. The coprime numerator has valuation zero there, so its pullback has a pole, a contradiction.

Now a polynomial D-constant belongs by (3.1) to L(Z), which lies inside the intermediate function field k(Z,U). Saturation for H_I makes it a polynomial in Z,U. On this polynomial ring D is exactly Delta acting coefficientwise in U. Each coefficient is therefore in C. This proves (3.2), and conversely every element of C[Z] is visibly annihilated. This regularity argument does not algebraize an arbitrary imprimitive intermediate field: H_I is already an explicit polynomial model.

## 4. UAJOINTSEL1: classify every polynomial joint selector

Combining (1.2)-(1.3) with UAJOINTKER1 gives the complete classification, for any c in k^times:

    det DT_beta=c
      iff beta=(1-c)*t+h(Z,C), h in C[Z].               (4.1)

In particular determinant-one joint selectors are EXACTLY the polynomials in processed active outputs, untouched inactive source coordinates and auxiliary polynomial first integrals. This closes the beta(P,U) gate, conditional only on the retained polynomial-profile and one-constant-direction architecture. It does not require computing C for an arbitrary auxiliary block.

Equivalently, after the explicit active-block map H_I, every allowed selector descends polynomially and becomes an affine function of the auxiliary time t, with coefficients in the auxiliary constant ring. This general theorem does not use the coefficient-leading-term or reciprocal sign argument of UASELKERNEL1. That source-only theorem remains a sharper description of the intersection with k[P].

## 5. UAJOINTFACT1: polynomial and stably polynomial auxiliary constants force factors

Assume the FULL auxiliary ring C is a polynomial algebra

    C=k[eta_1,...,eta_(m-1)].

Since every R coordinate is in C, write R=r(eta), and define H_aux(U)=(eta(U),t(U)). The identity M*G=(r(eta),t) and the determinant chain rule give

    det Dr(eta(U))*det DH_aux(U)=1.

Both factors are polynomials. Each is a nonzero constant, since their product is a unit in k[U]; dominance of eta also makes det Dr constant on its entire affine source. Thus r:A^(m-1)->A^(m-1) and H_aux:A^m->A^m are actual polynomial Keller self-maps. Rescaling one eta generator can normalize both determinants to one. If their generic degrees are a and b respectively, then

    deg G=a*b.                                         (5.1)

Take ANY joint-selector Keller member from (4.1). It has polynomial stages in the following order:

1. Process active compiler blocks, giving Z and retaining U.
2. Apply H_aux, followed by the polynomial automorphism

       (Z,eta,t) -> (Z,eta,c*t-h(Z,eta)).

3. Process inactive compiler blocks, using their coefficient expressions in eta and the retained signals.
4. Apply r on eta, retaining the last auxiliary coordinate, then the fixed output change M^-1.

The selector is evaluated before the inactive source coordinates are lost. The inactive profiles belong to C[sigma] coefficientwise, so they are recoverable from eta throughout Stage3. UAPARALLEL1 supplies every compiler stage on the whole affine source. All intermediate spaces are the SAME affine space of dimension3N+m, with different coordinate names. No invariant quotient variety is substituted for an affine self-map.

The stage degrees, omitting empty compiler stages and automorphisms, are

    product_active d_i, b, product_inactive d_i, a.

Hence deg T=deg G*product_i d_i. If deg G>1 and N>=1, at least one auxiliary stage and at least one compiler stage are noninvertible. Every such T is self-decomposable. Extra joint dependence of beta, any number of polynomial controls, constant directions and arbitrary degrees do not change this conclusion.

For the retained cubic, C is the polynomial linear output-quotient ring in EVERY nonzero constant direction by UACUBICFEED1. Here a=1,b=3. Splitting the compiler stages block by block gives the same explicit maximal atomic chain of length N+1 and degree multiset {3,d_1,...,d_N}, now for EVERY joint selector in the class and every nonzero constant Jacobian c. UAALLATOM1 proves the displayed compiler stages atomic; the auxiliary stage is polynomial left-right equivalent to the stabilized cubic. All maximal chains are not classified.

For N>=2, the two or more compiler stages already force decomposability even when G is an automorphism. At the complementary endpoint N=1 with G an automorphism, every auxiliary stage is an automorphism, so T is polynomial left-right equivalent to one stabilized compiler profile and is absolutely stably self-atomic by UAALLATOM1. This is a positive boundary within the known stable atomic classes, not a new canonical construction.

The polynomial-kernel hypothesis is strictly broader than output recovery. For G(B,V)=(K(B),K(V)) and return direction only in the second cubic, its full constant ring is k[B_1,B_2,B_3,R_1(V),R_2(V)]. It is polynomial, but B coordinates are not polynomials in the original auxiliary outputs K(B),R(V). The construction above first exposes these constants, processes the compiler and only then applies the transverse cubic K(B). Both auxiliary degrees are3. This is an actual family beyond the previous output-recovery case, not a new atomic example.

### Stable polynomiality and extension of constants

If C becomes polynomial after adjoining finitely many independent variables, stabilize G by the corresponding affine identity map, with return direction zero in those added coordinates. Its constant ring is exactly the enlarged C. Applying the proved factorization shows T is NOT stably self-atomic even when C itself is not known polynomial.

More generally field extension commutes with the kernel of the k-linear map Delta, by flatness. If for ANY K/k and any s>=0 the ring (C tensor K)[z_1,...,z_s] is polynomial, T_K times id_(A^s) decomposes. Thus absolute stable atomicity in this architecture requires the auxiliary constant ring to remain non-stably-polynomial under every field extension. This is only a necessary condition; no such ring or atomic receiver is constructed here.

### Locally nilpotent auxiliary directions

If Delta is locally nilpotent, its slice Delta t=1 gives k[U]=C[t] by the elementary finite Taylor formula. For f in k[U], the finite sum pi(f)=sum_j(-t)^j*Delta^j(f)/j! is constant, and f=sum_j t^j*pi(Delta^j f)/j!. Hence C is at least one-stably polynomial, giving stable decomposability immediately.

There is also an unstabilized factorization in this case. After H_I, the change fixing Z,C and sending t to c*t-h(Z,C) is a polynomial automorphism of k[Z,U]=C[Z,t]. It fixes every inactive profile. The remaining map is therefore the compiler with inactive profiles, followed by G, composed with that automorphism. Together with the active compiler stages, this gives actual nontrivial compiler and auxiliary factors whenever deg G>1. Local nilpotence must be proved; a polynomial slice alone does not suffice.

## 6. Natural marked fibres exclude the all-active symmetric case

Suppose N>1, S_N permutes compiler blocks, fixes the auxiliary source, the profiles are covariant, beta is invariant, and the target mark is fixed. The active set is then all blocks or none. In the all-active case, (4.1) writes beta as (1-c)*t+h(Phi,C). At the marked target Phi is fixed, so the auxiliary equation depends on U alone. For each such U all block targets and coefficient profiles specialize identically. Its compiler fibre is a product T_U^N.

If the complete marked fibre is nonempty, choose any of its geometric points and repeat one block preimage in all N slots, keeping the same U. The formula for beta still holds at the fixed Phi target, so this diagonal point lies in the complete fibre and is fixed by S_N. The natural transitive N-label fibre has no S_N-fixed point for N>1. Thus the all-active case cannot retain that complete canonical fibre, for ANY auxiliary G; polynomiality of C was not used.

Consequently a natural marked one-direction candidate with fixed auxiliary symmetry must be all-inactive. Its joint selector is affine in t with remaining coefficients in C and the full compiler sources. UAJOINTFACT1 excludes it whenever C is polynomial or stably polynomial (over the relevant constants), and excludes locally nilpotent directions even without choosing polynomial constant generators. UACANFIB1 is correctly in the all-inactive, non-atomic class. No full marked fibre for a different auxiliary geometry is asserted.

## 7. What remains

Joint selectors no longer supply an unclassified gate in the polynomial retained-signal, one-constant-return architecture. Their entire kernel and affine-slice classification are (3.2)-(4.1). The new factorization covers every auxiliary map with polynomial or stably polynomial constant ring, including the cubic and products whose invariants were not recoverable from the original outputs.

Remaining departures include genuinely higher-rank return, profiles beyond the retained-signal polynomial class with a new whole-source Jacobian/landing proof, another compiler, or a non-locally-nilpotent auxiliary direction with a constant ring that does not become stably polynomial. In the natural fixed-auxiliary symmetric case all compiler coefficients must then be first integrals. None of these conditions proves existence or absolute atomicity. Canonical absolute atomic realization, improved dimensions/degrees and stronger receiver requirements remain open. The general proofs above are not replaced by the selected exact controls, and novelty and independent review remain unverified.
