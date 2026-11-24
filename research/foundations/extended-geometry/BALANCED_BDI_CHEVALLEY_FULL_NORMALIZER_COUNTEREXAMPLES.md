# Balanced real Grassmann full-normalizer product counterexamples

## 1. Statement and exact input

Let k>=2, let q0=sum_(i=1)^(2k) z_i^2 over Q, and let G0=SO(q0). Set S0=diag(I_k,-I_k). The Q-involution theta0(g)=S0 g S0 is defined even when det S0=-1. Put

    H=G0^theta0=S(O_k times O_k).

Let X be the actual rank-k self-adjoint projector model from the canonical classical proof, with base point P0=diag(I_k,0). Its affine quotient by P maps to I-P will be denoted Y. The assertion proved here is an explicit compatible Q-identification: the compact simply connected Chevalley D_k model with its Chevalley involution has adjoint full-fixed quotient exactly Y. For k=2 use D2=A1 times A1, and for k=3 use the D3=A3 root system. The k=1 torus boundary is excluded.

The mathematical construction imported afterward is [the character-domain theorem](SPLIT_COMPACT_SYMMETRIC_CHARACTER_DOMAIN_PRODUCT_COUNTEREXAMPLES.md), SHA256 `9c6bc0b8002820ed7dd46e456c42e718abcd04b787d3134f8a298a6b11f2a185`, including its semisimple-product corollary. The projector scheme geometry comes from [the classical-model proof](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), SHA256 `77e79160d3cc1478cabd9eea55fdc4c0374ac937312157a59fb90c99d047b9de`, Section2.1. This is a Q-model identification, not a real-form label or an assertion for arbitrary quadratic forms.

**Theorem 1.** For each k>=2, the standard positive rank-k real projector variety in dimension 2k, modulo its complement involution, is the actual adjoint full-fixed compact Chevalley D_k quotient constructed below. Over every characteristic-zero field it has the specified rank-k character-domain phase F and fixed unit marker t^k, with all unmarked powers of F having trace zero and a nonzero marked value at every positive order, with absolute unit-nonpower, primitive Laurent support and the exact seed-power/monomial-marker cutoffs in Section 5. Over every formally real field its two individual factors are Mathieu and their product is not. The D2 semisimple case and D3=A3 identification are included; k1 is excluded.

## 2. Split matrix model and explicit last-root signs

Let V_sp=Q^k plus Q^k with coordinates (x,y) and split quadratic form

    q_sp(x,y)=sum_i x_i y_i.

Write J=[[0,I_k],[I_k,0]], so its polar matrix is J and J squared=I. Let G_sp,vec=SO(J). Its split torus is diag(t_1,...,t_k,t_1^-1,...,t_k^-1). Denote matrix units by E_ab, and put H_i=E_ii-E_(k+i,k+i).

For i<j choose positive and negative root matrices

    X_(e_i-e_j)=E_ij-E_(k+j,k+i),
    X_(-e_i+e_j)=E_ji-E_(k+i,k+j),
    X_(e_i+e_j)=E_(i,k+j)-E_(j,k+i),
    X_(-e_i-e_j)=E_(k+j,i)-E_(k+i,j).

Every matrix X satisfies X^T J+JX=0, its square is zero, and X_negative=X_positive^T. Direct matrix-unit multiplication gives

    [X_(e_i-e_j),X_(-e_i+e_j)]=H_i-H_j,
    [X_(e_i+e_j),X_(-e_i-e_j)]=H_i+H_j.

Thus x_alpha(u)=I+u X_alpha is the usual pinned root subgroup on the vector representation. In particular for the simple roots alpha_i=e_i-e_(i+1), 1<=i<=k-1, and alpha_k=e_(k-1)+e_k, the LAST root has the explicit signs

    X_alpha_k=E_(k-1,2k)-E_(k,2k-1),
    X_-alpha_k=E_(2k,k-1)-E_(2k-1,k),
    [X_alpha_k,X_-alpha_k]=H_(k-1)+H_k.

The first k-1 simple pairs have coroot H_i-H_(i+1). These formulas include k2 and k3 with their stated semisimple/isomorphic root systems.

The automorphism theta_sp(g)=g^-T equals JgJ on SO(J), by the defining orthogonal equation. On the chosen pinning it satisfies

    theta_sp(x_alpha(u))=x_-alpha(-u),
    dtheta_sp(X_alpha)=-X_-alpha, dtheta_sp(H_i)=-H_i.

Indeed inverse transpose of I+uX is I-uX^T. This checks the Chevalley convention including the last simple root, rather than deducing it from an unspecified transpose convention. Lift this pinned involution to the simply connected split Spin cover; it is precisely the Chevalley involution used in the canonical construction.

## 3. Descent of the standard module over Q(i)

Let E=Q(i), with scalar conjugation sigma. Descend the split vector space by the semilinear involution

    rho(v)=J sigma(v).

Its fixed vectors are EXACTLY

    (x,y)=(a+ib,a-ib), a,b in Q^k.

The displayed vectors are a Q-basis after allowing all a,b, and scalar extension to E is the original 2k-dimensional vector space. Restriction of q_sp to this fixed Q-space is

    sum_i(a_i+ib_i)(a_i-ib_i)=sum_i(a_i^2+b_i^2)=q0(a,b).

No square root of2 or unverified similarity of rational quadratic forms occurs. Equivariance gives descent on the orthogonal group by

    g maps to J sigma(g) J = theta_sp(sigma(g)).

Consequently the descended vector group is exactly SO(q0) over Q as an algebraic group. This is verified on the vector space and form, so is not merely a coincidence of compact real points.

The linear involution J commutes with rho. On fixed vectors it sends (a+ib,a-ib) to (a-ib,a+ib), hence acts as (a,b) maps to (a,-b). Therefore the descended Chevalley involution on SO(q0) is EXACTLY theta0, conjugation by S0. The descended simply connected group maps onto SO(q0); by uniqueness of the simply connected cover of a connected semisimple Q-group it is its standard Spin cover, with the lifted compatible involution. This applies also to the semisimple D2 group.

## 4. Full adjoint fixed group adds exactly the complement

For k>=2 the center of G0 is {+I,-I}, so G_ad=PSO(q0). Consider the FULL fixed subgroup H_ad=(G_ad)^theta_ad, and its inverse image N in G0. After algebraic closure every projective orthogonal class lifts. The condition is exactly

    g S0 g^-1 S0^-1 in {+I,-I},

or equivalently g S0 g^-1=+S0 or -S0. The plus case preserves the two k-dimensional eigenspaces and is H. The minus case interchanges them. The rational matrix

    R=[[0,I_k],[-I_k,0]]

has determinant1, is orthogonal for q0, interchanges them, and has R squared=-I in H. Thus

    N=H union RH, N/H=constant Z/2 over Q.

All these equations define smooth group schemes in characteristic zero, so the geometric equality descends to Q and does not ignore disconnected fixed components. Its right coset action on G0/H is literally P maps to I-P, because R P0 R^-1=I-P0. It is fixed-point-free on the idempotent variety: P=I-P would require P=I/2, contrary to P squared=P. The invariant quotient is therefore finite etale, smooth geometrically integral affine, and

    Y=(G0/H)/Z2=G0/N=PSO(q0)/H_ad.

The last equality follows by quotienting the central kernel {+-I}, which is contained in H. This is exactly the adjoint full-fixed quotient of the descended simply connected Chevalley group in Section3, a central quotient explicitly allowed by the canonical character-domain theorem.

For completeness, N is also the full normalizer of the connected block subgroup H0=SO_k times SO_k. Work over an algebraic closure. Its Lie algebra h is the plus eigenspace of theta0 on so_(2k); the minus eigenspace p is its Killing orthogonal complement, since the Killing form is nondegenerate for semisimple so_(2k), including k2. A normalizer preserves both eigenspaces, so its adjoint action commutes with theta0. Hence g S0 g^-1 S0^-1 lies in the kernel of Ad, namely {+-I}, giving g in N. Conversely every block-preserving or block-interchanging g normalizes H0. Smooth geometric equality descends. It also normalizes H, and any normalizer of H normalizes H0, so N is the full normalizer of H as well. No accidental low-rank extra permutation is omitted.

## 5. Consequences with their exact scope

The installed character-domain construction on the adjoint full-fixed quotient now applies to this ACTUAL standard balanced BDI complement quotient over Q and every characteristic-zero extension. In its notation, use its specified Coxeter parameter qstar, its rational character tuple Phi and adjoint discriminant d. With T=1+t and y=Phi(c)-qstar, the upstairs recipe is

    f(y)=d(qstar+y)/d(qstar),
    p=T/t f(Ty)^2, F=p^k(1+p), H_marker=t^k.

The theorem proves that these actual coefficients descend to the adjoint full-fixed quotient identified with Y above; individual character coordinates are not presumed to descend. For mu_j=lambda_Q(f(y)^(2j))>0 rational, the exact every-order values are

    Lambda(F^m)=0,
    Lambda(t^k F^m)=sum_(h=0)^m binom(m,h)mu_(km+h)>0, m>=1.

Absolute unit-nonpower and primitive Laurent support are also part of that same actual-quotient theorem. For positive integer a, p^a is pure IFF a>=k. Given a>=k, t^ell is unsafe for p^a IFF ell>=k, and survives EVERY positive order IFF k<=ell<=a. For F the corresponding conditions are ell>=k and k<=ell<=k+1. These are properties of this displayed character-domain phase, not a universal lower bound across all possible phases.

Over every formally real extension, the compact quotient factor and split Gm are separately Mathieu, whereas their product is not. Both trace restriction and scalar transfer are the canonical same-model statements; no failure is descended from an arbitrary upstairs phase.

This closes the previously conditional standard-Q identification for odd balanced k>=3. Even balanced k also admits the separate symmetric polynomial-CDF construction from the Beta(k/2,k/2) diagonal coordinate; it can have a simpler seed, without contradicting the rank cutoff for a different character-domain seed. The excluded k1 is a torus boundary, for which this semisimple argument and its product-failure conclusion are not asserted.

## 6. Attribution and assurance

The matrix realization, Chevalley involution and quadratic descent are classical constructions, written here with exact coordinates to match the actual rational model. The counterexample is an application of the already proved canonical character-domain theorem, with its retained polar/Steinberg and endpoint attributions. No novelty, external review, formal verification, independent replay or worldwide priority is claimed. The mathematical argument is complete and uses the stated exact model scopes.
