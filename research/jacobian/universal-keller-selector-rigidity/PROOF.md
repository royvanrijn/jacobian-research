# Complete selector rigidity for polynomial feedback in one auxiliary direction

15 November 2026. Written all-degree proofs with selected exact controls. No independent whole-proof replay, formal verification, external review or publication-priority certificate. OP-KDESC remains open. This closes the different-source-selector gate from [rank-one feedback](../universal-keller-rank-one-feedback/PROOF.md); it does not construct an absolutely atomic canonical receiver.

Work over a characteristic-zero field k. The return selector beta is an arbitrary polynomial in the COMPILER source blocks P, independent of the auxiliary source U. Coefficient profiles may be arbitrary polynomials in U and all retained signals. Auxiliary-dependent selectors beta(P,U), higher-rank return, rational profiles and other compilers are outside the classification.

## 1. Leaf coordinates and the exact response

For block i set t_i=1+x_i*y_i, q_i=t_i^2*z_i+y_i^2*(1+3*t_i), and sigma_i=t_i*q_i as in UAM1. On the function field introduce

    S_i=x_i/t_i,    r_i=1/t_i.

Then

    x_i=S_i/r_i,    y_i=(1-r_i)/S_i,
    z_i=sigma_i*r_i^3-y_i^2*(r_i^2+3*r_i).

Consequently k(P)=k(sigma_1,...,sigma_N)(S_1,r_1,...,S_N,r_N), with all displayed variables algebraically independent over the indicated base field. These coordinates are used only to prove function-field identities; the maps and factors constructed later are polynomial in the original coordinates on the whole source.

Let g_ij(U,sigma) be any finite family of polynomial higher coefficient profiles, j>=4. Write

    A_i(U,S_i)=sum_j g_ij(U,sigma)*sigma_i^j*S_i^j,
    J_i=S_i^2*partial_(S_i)^2 A_i-S_i*partial_(S_i) A_i.

Let Delta be ANY polynomial derivation of k[U], extended to fix all signals, and put H_i=Delta A_i. All coefficient fields below may include the signals; Delta need not be locally nilpotent, and no group action for Delta is assumed. The response field W=(D_P Phi)^(-1)*V*(Delta g) from UARANKONE1 preserves every sigma. In leaf coordinates its block components are

    W_i(S_i)=H_i/r_i,
    W_i(r_i)=H_i/S_i-partial_(S_i)H_i
       +H_i*(-1/S_i+3*sigma_i*S_i+J_i/S_i^2)/r_i.       (1.1)

Indeed, with all profiles frozen at their values, the two remaining outputs satisfy

    Bout_i=(r_i-1-3*sigma_i*S_i^2-partial_(S_i)A_i)/(2*S_i),
    Cout_i=S_i*(r_i+1)-sigma_i*S_i^3
                         -(S_i*partial_(S_i)A_i-2*A_i).

Their Jacobian in (S_i,r_i) is -r_i/S_i. The coefficient response is

    (-partial_(S_i)H_i/(2*S_i), -(S_i*partial_(S_i)H_i-2*H_i)).

Inverting that two-by-two matrix gives (1.1). Although profiles depend on every signal, d sigma*W=0 makes all profile-signal derivatives vanish along W; UAPARALLEL1 therefore justifies the frozen-block computation. This is not a diagonal-profile assumption.

## 2. A reciprocal-derivation parity lemma

Let F be any characteristic-zero field, and let S_i,r_i, i in a nonempty finite set T, be independent variables over F. Choose nonzero f_i in F(S_i). For

    V=sum_(i in T) f_i(S_i)/r_i * partial_(r_i),

every rational first integral of V is invariant under EACH independent sign change r_i -> -r_i, fixing the S variables and all other r variables.

Here is the exact constant-field proof. Enlarge the coefficient field by all S_i and any untouched variables, calling it F0. Set t=r_1^2/(2*f_1) and

    d_i=r_i^2/(2*f_i)-t for i other than 1;    d_1=0.

Then V(t)=1, V(d_i)=0. The field E=F0(r_i:i in T) is

    F0(d_i:i!=1)(t)( sqrt(2*f_i*(t+d_i)):i in T ).       (2.1)

The t+d_i are distinct independent linear factors. Every nonempty product of their radicands has odd valuation at at least one of these factors, even after algebraic closure of the coefficient field F0(d). Their square classes are independent, so (2.1) has degree 2^|T| both arithmetically and geometrically. Thus F0(d) is relatively algebraically closed in E. A nonzero derivation in a one-variable function field in characteristic zero has precisely that relative algebraic closure as its constant field: a transcendental constant would make the whole derivation zero on a finite separable extension. It follows that

    ker(V on E)=F0(d).

This field is fixed by every sign change. No claim that V integrates to a complete additive action is used. Distinct-factor valuations, not a bounded group enumeration, prove the statement for every |T|.

## 3. UASELKERNEL1: the exact common first-integral ring

Expand W in the algebraically independent U variables:

    W=sum_gamma U^gamma*W_gamma(P).

Call a block ACTIVE when some Delta g_ij is nonzero. Let I be the active set and I^c its complement. Then the complete rational and polynomial common kernels are

    intersection_gamma ker(W_gamma on k(P))
       =k(sigma_i:i in I, P_l:l in I^c),                (3.1)

    intersection_gamma ker(W_gamma on k[P])
       =k[sigma_i:i in I, P_l:l in I^c].                (3.2)

If I is empty this is the whole source field/ring. Any number of polynomial auxiliary controls, every block count and every finite collection of degrees is allowed. There is no degree bound on a proposed rational or polynomial first integral.

### Highest auxiliary term and sign separation

Put F=k(sigma_1,...,sigma_N) and work in F(S_i,r_i). Choose positive integer weights on U which distinguish the finitely many distinct exponent vectors occurring in the A_i,H_i and their products. For every active block let D_i be the leading weighted degree of A_i and E_i that of H_i. Since Delta annihilates U-independent coefficients, D_i>0, whereas E_i>=0. The operator A -> S^2*A''-S*A' cannot kill a nonzero polynomial with only powers S^j, j>=4. Thus J_i has the same leading auxiliary monomial as A_i.

Let M=max_(i in I)(D_i+E_i). It is STRICTLY greater than every E_i. Therefore a leading auxiliary coefficient of W has only the quadratic product term from (1.1), and has the form

    V=sum_(i in T) f_i(S_i)/r_i * partial_(r_i),         (3.3)

for a nonempty subset T of the active blocks. Here

    f_i=(S_i^2*c_i''-S_i*c_i')*h_i/S_i^2,

where c_i,h_i are the leading U coefficients of A_i,H_i respectively. Both are nonzero polynomials of S order at least four. The weight choice ensures that blocks attaining the same leading weight have the same leading exponent vector; no other terms can cancel that coefficient. All remaining source variables may be treated as coefficients in the parity lemma.

Let beta be a rational common first integral. It annihilates V, so the parity lemma makes it invariant under every sign change r_i -> -r_i for i in T. In (1.1), the part of W_i odd under this sign change is

    Y_i(H_i)=H_i/r_i * partial_(S_i)
              +(H_i/S_i-H_i')*partial_(r_i).

Every other block and the remaining part of W_i are even. Applying the sign change to W(beta)=0 and subtracting therefore gives Y_i(H_i)(beta)=0 for every U. Taking its leading U coefficient gives Y_i(h_i)(beta)=0.

Both V and Y_i(h_i) annihilate beta, so their commutator does too. It has support only in block i. Put S=S_i,r=r_i,h=h_i,f=f_i and l=h/S-h'. The two block vectors are

    Y=(h/r,l),
    [V,Y]=(-h*f/r^3, (l*f-h*f')/r^2).

Their determinant is

    h*(2*l*f-h*f')/r^3.                                 (3.4)

It is nonzero for EVERY c,h allowed above. Indeed, with Jc=S^2*c''-S*c', set

    P0=h^3*Jc/S^4.

This is a nonzero polynomial of strictly positive S order (at least12), hence nonconstant. Since f=Jc*h/S^2,

    2*l*f-h*f'=-h*f*P0'/P0 !=0                         (3.5)

in characteristic zero. Equations (3.4)-(3.5) force beta independent of both S_i and r_i. Remove every block in T and repeat. There are finitely many active blocks, so the induction removes all of them. This proves that beta belongs to the right side of (3.1). Conversely each signal and each inactive source coordinate is annihilated by W, giving equality.

### Polynomial regularity

The polynomial projection

    P -> (sigma_i for i in I, P_l for l in I^c)

is dominant and smooth: it is a product of retained first-output projections (each the cubic Keller map followed by a linear projection) and identities. Its relative dimension is 2*|I|. A reduced rational function on its affine-space target cannot pull back regularly with a nonconstant denominator. Such a denominator pulls back to a nonconstant nonunit; a height-one prime in its zero scheme contracts to the original divisor by flat going-down. Its coprime numerator does not vanish there. This gives a pole, a contradiction. Thus intersection with k[P] gives exactly the polynomial ring (3.2). No rational invariant is silently assumed to descend polynomially.

## 4. UAALLSELECT1: classify EVERY compiler-source selector

Take ANY polynomial auxiliary Keller map G:A^m->A^m with determinant1 and a nonzero constant vector mu. Let Delta=(DG)^(-1)*mu. For the arbitrary polynomial profiles above and ANY beta(P) in k[P], consider

    T(P,U)=(Phi(P;g(U,sigma(P))), G(U)-beta(P)*mu).

The following are equivalent:

1. T is Keller, with any nonzero constant determinant.
2. Its determinant is identically one.
3. beta is a polynomial in the retained signals of ACTIVE blocks and the complete source coordinates of INACTIVE blocks:

       beta in k[sigma_i:i in I, P_l:l in I^c].          (4.1)

The Schur identity of UARANKONE1 says det DT=1+W(beta). At P=0 the response vanishes, so any constant is one. UASELKERNEL1 then proves the classification. This removes the b(C) restriction of UARANKONE1 and closes its common-kernel gate for selectors independent of U. It does not assert Delta g=0 for all allowed selectors; the earlier signal-only counterexample remains valid.

In particular if every block is active, every Keller selector is signal-only. If a selector depends on every block beyond its signal, all coefficient profiles must be Delta first integrals. An arbitrary polynomial source selector cannot simultaneously separate feedback-active blocks and preserve the Keller property.

## 5. UAONECHANNEL1: active blocks, auxiliary map, inactive blocks

Assume additionally that the auxiliary first-integral ring is output-recoverable:

    ker(Delta on k[U])=k[l_1(G),...,l_(m-1)(G)],         (5.1)

for independent linear forms annihilating mu. UAGEOMFEED1 gives the exact geometric generic-line connectedness criterion for (5.1). The retained auxiliary cubic satisfies (5.1) in EVERY nonzero constant direction by UACUBICFEED1.

Every Keller T from Section 4 has actual polynomial Keller stages in this order:

1. Process all active compiler blocks, leaving U and inactive source blocks unchanged.
2. Apply the auxiliary G and subtract beta*mu.
3. Process all inactive compiler blocks, leaving the already processed blocks and auxiliary outputs unchanged.

To construct Stage1, use UAPARALLEL1 partial processing with U as promoted identity coordinates. Every original signal is available after any subset has been processed: it is the first coordinate of a processed block or the original signal polynomial of an untouched block.

By (4.1), beta can be recovered after Stage1 from those active first coordinates and the untouched inactive source coordinates. Thus Stage2 is a genuine polynomial map. It is a target translation of G times identities and has degree deg(G). No inverse compiler is used.

For an inactive block Delta g_ij=0 coefficientwise in the independent signal variables. Under (5.1) write g_ij=h_ij(l(G),sigma). If Vout=G(U)-beta*mu, then l(G(U))=l(Vout), so every inactive profile is recoverable after Stage2. UAPARALLEL1 now constructs Stage3. These identities prove equality with T on the entire affine source.

The three generic stage degrees are

    product_(i in I) d_i,   deg(G),   product_(i notin I) d_i,

where empty block stages are omitted. Hence deg T=deg(G)*product_i d_i. If N>=1 and deg(G)>1, the auxiliary stage and at least one compiler stage are nontrivial. Every such T is decomposable, including all constant extensions and affine identity stabilizations.

For the cubic, process each block separately within Stages1 and3. This gives N+1 noninvertible polynomial Keller stages with degrees {3,d_1,...,d_N}. Each compiler stage is an allowed single-block profile with all other coordinates retained as identities. UAALLATOM1 makes it absolutely stably self-atomic, and the auxiliary stage is polynomial left-right equivalent to the cubic times identities. Thus this particular chain is maximal and atomic, of length N+1. We do NOT classify all maximal decomposition chains of T.

### The canonical symmetry boundary

If S_N permutes the compiler blocks, fixes the auxiliary source, the profiles are covariant and the selector is invariant, the active set is invariant. For N>1 it is therefore either empty or all blocks. When all are active the factors are compiler first, auxiliary second; when none are active they occur in reverse order. These are compatible with the given covariance. Thus no polynomial compiler-source selector, of any degree, can make this one-direction cubic architecture an atomic canonical replacement. The non-atomic canonical UACANFIB1 construction remains valid.

This is a factorization obstruction even before checking a proposed marked fibre. It does not claim that arbitrary choices of profiles and selector retain that fibre, or that the original canonical fibre can survive every forward-order member.

For another auxiliary G WITHOUT (5.1), Section 4 remains valid. If every block is active and deg(G)>1, the forward two-factor construction already applies, without (5.1). Consequently an S_N-covariant noninvertible auxiliary-only candidate can escape these proved polynomial factors only when every block is inactive and some polynomial first integrals are not output-recoverable. It faces the retained generalized closed-subring/geometric-line gate. No such auxiliary map or atomic construction is supplied.

## 6. What has closed, and what remains

The proposed separating common invariant beyond active signals does NOT exist in this class. Its exact rational and polynomial common kernels are now computed, including mixed active/inactive blocks. Replacing the old selector b(C) by a more complicated polynomial in P cannot remove the cubic factor: the proper order of processing changes, but polynomial factors remain.

A new atomic construction must leave a stated hypothesis: for example genuinely higher-rank return, selector dependence mixing P with U and its own Jacobian analysis, profiles outside polynomial retained-signal dependence with a separate whole-source landing proof, an auxiliary map with non-output-recoverable constants in the inactive case, or a different compiler. These are necessary departures, not constructions. The original objective, canonical absolute atomicity, quantitative improvements and stronger receiver requirements remain open. No numerical rank bound or finite first-integral ansatz is used to infer this theorem.
