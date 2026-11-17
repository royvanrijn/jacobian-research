# Arbitrary polynomial feedback through one auxiliary direction

15 November 2026. Written general proofs with selected exact controls. These extend the scalar classification in [feedback rigidity](../universal-keller-feedback-rigidity/PROOF.md). Independent whole-proof verification, formal verification, external review and publication priority are not established. OP-KDESC remains open.

Work over a characteristic-zero field k. Every feedback classification here fixes the return scalar beta=b(C), where C is the retained quadratic trace below and b is nonconstant. Arbitrary other source functions beta(P) require a new response analysis and are not classified. All maps below are polynomial on the WHOLE stated affine source. Self-factorizations have the same affine space as intermediate source and target. Absolute stable atomicity includes extension of constants as well as affine identity stabilization.

## 1. UARESPINJ1: the coefficient-response functional is injective

Use the retained normalized compiler on N>=1 blocks P_i=(x_i,y_i,z_i):

    t_i=1+x_i*y_i,
    q_i=t_i^2*z_i+y_i^2*(1+3*t_i),
    sigma_i=t_i*q_i,
    Phi_i=(sigma_i,
      -y_i/2-3*x_i*q_i/2
        -sum_(j>=4) (j/2)*a_ij*t_i^2*x_i^(j-2)*q_i^j,
      x_i*(5-3*t_i)-x_i^3*z_i
        -sum_(j>=4) (j-2)*a_ij*(x_i*q_i)^j).             (1.1)

There are finitely many active indices, with no bound on their degrees. Let U be any additional independent coordinates. Substitute arbitrary profiles

    a_ij=g_ij(U,sigma_1,...,sigma_N).

Freeze U when taking A=D_P Phi. UAPARALLEL1 gives det A=1. Write V for the matrix of independent coefficient-response columns, evaluated at g; its columns have zero first output in each block. For arbitrary polynomial coefficient directions h_ij(U,sigma), set

    W_h=A^-1*V*h.

Let ell(x)=(3*x^2-x)/2, C(P)=sum_i ell(x_i), and beta(P)=b(C(P)) for ANY nonconstant polynomial b in k[T]. Then the k[U]-linear map

    h |-> d_P beta * W_h                                (1.2)

is INJECTIVE on all finite polynomial coefficient directions. In particular its value is zero exactly when every h_ij is zero. The profiles and directions may have unrelated dependence on U and all retained signals. This is stronger than checking directions of the form partial_a g for one scalar a.

### Proof: retain all signals and isolate one block

Every signal differential satisfies d sigma_i*A^-1*V=0. Thus W_h preserves ALL original signals, even when the profile dependence has cycles. Since d_P g is a linear combination of those differentials, the signal-profile part of A makes no contribution along W_h. Each block response can therefore be computed using the compiler with its coefficient values frozen.

Choose a block i with some h_ij nonzero. Put x_l=y_l=0,z_l=pi_l in every other block l. Their signals pi_l remain algebraically independent. Their coefficient-response columns vanish because the higher terms contain x_l to power at least two. Their source responses are therefore zero. In the active block set y_i=0, x=x_i, z=z_i; its signal is z. Localize at x and put w=xz, so z=w/x. All other signal arguments and U remain independent coefficient variables.

Suppress i and write

    B=sum_j h_j(U,w/x,pi_other)*w^j,
    C0=sum_j j*(j-2)*g_j(U,w/x,pi_other)*w^j.

The exact first-coordinate response is

    (W_h)_x=sum_j (j+1)*h_j(U,w/x,pi_other)*w^j
             -3*x*w*B-x^-1*C0*B.                        (1.3)

Here is a direct derivation, rather than an extrapolation from examples. The inverse equation and reconstruction are

    E(S)=S+Bout*S^2+sigma*S^3
           +sum_j a_j*sigma^j*S^j-Cout/2=0,
    x=S/E'(S),    y=(1-E'(S))/S.                         (1.4)

At y=0, S=x and E'=1. The target variation Vh gives W_h(S)=B. With coefficients frozen and signals constant,

    W_h(E')=(3*w+C0/x^2)*B
                 -x^-1*sum_j j*h_j*w^j.

Differentiating x=S/E' gives (1.3). Thus no derivative of a signal-dependent profile was dropped without justification.

B is nonzero. Indeed a monomial of active signal degree e in h_j gives x^-e*w^(j+e); the pair of exponents uniquely recovers (e,j). There is no cancellation between distinct coefficient directions. If -m is its largest x power, m>=0, the middle term of (1.3) has largest x power 1-m. The first term has largest power at most -m and the last at most -1-m, because C0 has only nonpositive x powers. The leading term cannot cancel. Finally d beta on this slice multiplies this response by b'(ell(x))*ell'(x), which is nonzero in characteristic zero. This proves injectivity.

Extra independent identity coordinates can be included in U; the argument is unchanged. Negative powers of the signals and positive-characteristic zero derivatives are outside this theorem.

## 2. UARANKONE1: exact classification for an arbitrary auxiliary Keller map

Let G:A^m->A^m be ANY polynomial map with det DG=1, m>=2, and let mu be a nonzero constant vector in k^m. Define the polynomial derivation

    Delta=(DG)^-1*mu on k[U].                            (2.1)

For the arbitrary profiles g of Section 1, form

    T_g(P,U)=(Phi(P;g(U,sigma(P))), G(U)-beta(P)*mu).
                                                                    (2.2)

Then the following are equivalent:

1. T_g is Keller, with any nonzero constant Jacobian.
2. Its Jacobian is identically one.
3. Delta(g_ij)=0 for EVERY effective coefficient profile, holding all signals constant.

No restriction is placed on the number of independent source controls through which the coefficients were presented. The assertion concerns the resulting coefficients g, not redundant controls which might disappear under substitution.

### Proof

Write D=DG and B_U=D_U Phi=V*D_U g. The Schur complement gives

    det DT_g=1+d beta*A^-1*V*(Delta g).                  (2.3)

All matrices and inverses used here are polynomial. The second term is exactly the injective functional (1.2) evaluated on h=Delta g. Thus det DT_g=1 exactly when every coefficient is a Delta first integral. At P=0 all response columns V vanish, so det DT_g=1 there for every U. Any constant determinant must consequently equal one. This proves all three equivalences.

The feedback is rank one in the direction FROM the compiler TO the auxiliary block: it is one scalar beta times a constant vector mu. Merely adding more controls in the other direction does not change this rank. The full small-matrix Schur equation from the earlier note allows more possibilities only after leaving this rank-one architecture, or after finding effective polynomial first integrals not recoverable from the auxiliary outputs.

## 3. UAGEOMFEED1: the exact geometric gate for output recovery

Choose independent linear forms l_1,...,l_(m-1) annihilating mu and a linear form tau with tau(mu)=1. Put

    R=(l_1(G),...,l_(m-1)(G)),
    B=k[R], F=k(R), E=k(U),
    L=the relative algebraic closure of F in E.

The quotient map pi=R:A^m->A^(m-1) is smooth of relative dimension one: it is G followed by a smooth linear projection. Delta(R)=0 and Delta(tau(G))=1. Then

    ker(Delta on E)=L,                                  (3.1)

and the following are equivalent:

(a) The geometric generic fibre of pi is connected.
(b) That fibre is integral.
(c) F is algebraically closed in E.
(d) ker(Delta on k[U])=B.

The fibre is smooth, so its geometric irreducible components are disjoint. This explains the equivalence of connectedness and integrality here; no singular-fibre assertion is made.

### Rational constants and regular descent

If a Delta-constant f were transcendental over F, then E would be finite algebraic over F(f), since trdeg_F E=1. In characteristic zero the zero derivation on F(f) has only the zero extension to E, contradicting Delta(tau(G))=1. Conversely every element algebraic over F is annihilated by Delta, by differentiating its separable minimal polynomial. This proves (3.1). The classical geometric-irreducibility/relative-algebraic-closedness criterion gives (a)-(c); smoothness supplies geometric reducedness.

If L=F, a polynomial first integral is a reduced fraction p(R)/q(R). A nonconstant irreducible factor q of its denominator pulls back to a nonconstant nonunit polynomial: dominance excludes constant pullback, and the only units on A^m are constants. Pick a prime divisor in its zero scheme. Smoothness and relative dimension one force this divisor to dominate the corresponding divisor of A^(m-1), rather than a codimension-two subset. Equivalently use flat going-down and the dimension formula at its generic point. Its valuation of p is zero, giving a pole in p/q. Hence a regular pullback has constant denominator, proving (d).

Conversely suppose L is larger than F. For a in L outside F, clear the finitely many denominators of its monic minimal polynomial: some nonzero q in B makes q*a integral over B. Since E=k(U) and k[U] is integrally closed, q*a belongs to k[U]. It is a Delta-constant outside F and hence outside B. This disproves (d). Thus the criterion is an equivalence, not only a sufficient test.

### Polynomial reverse factors and the published closedness boundary

Under these equivalent conditions, Section 2 writes every allowed profile as

    g_ij=h_ij(R,sigma).

Set Vout=G(U)-beta(P)*mu. Then R(U)=l(Vout), and

    H(P,U)=(P,G(U)-beta(P)*mu),
    Q(P,Vout)=(Phi(P;h(l(Vout),sigma(P))),Vout),
    T_g=Q composed with H.                              (3.2)

These are actual whole-source polynomial Keller self-maps. Their generic degrees are respectively deg(G) and product_i d_i>=3^N, where d_i are the effective compiler degrees. In particular if deg(G)>1 this is a nontrivial self-factorization. If G is an automorphism, H has degree one, so this displayed two-factor formula alone does not assert non-atomicity; for N>1 the output-controlled block factors of Q still give it.

A change of auxiliary map which evades THIS output-recovery obstruction while retaining (2.2) must therefore have disconnected geometric generic inverse images of its mu-directed lines. Equivalently B must fail to be algebraically closed in k[U]. Its Jacobian minors generate the unit ideal, because the rows of D(pi) extend to the invertible matrix DG after a linear output change. Thus such an escape would violate the corresponding algebraic-closedness assertion JC_sub(m-1,m;k), in the notation of the generalized Jacobian conjecture discussed by Jedrzejewicz and Zielinski. It is not asserted to disprove the ordinary plane Jacobian conjecture: failure of a stronger generalized assertion does not have that logical consequence.

For m=2 the one-generator closedness theorem JC_sub(1,2;k) is classical and proved (Ayad; as recorded in Section 1 of that primary paper). Consequently every plane auxiliary Keller map satisfies the output-recovery condition. Polynomial automorphism auxiliaries in every dimension also satisfy it, since pi is a coordinate projection transported by an automorphism. Neither observation supplies a new noninvertible auxiliary example.

## 4. UACUBICFEED1: every constant cubic direction and every polynomial control

Take the retained auxiliary cubic K. For EVERY nonzero constant direction mu in k^3, not only the particular filter direction, its quotient pi=l(K) has rational geometric generic fibre. Indeed in output coordinates a=K1,b=K2,c=K3 the inverse equation is

    S+b*S^2+a*S^3-c/2=0.

Write (a,b,c)=s(R)+tau*mu for a linear section s. The equation is linear in tau, with nonzero coefficient

    mu_1*S^3+mu_2*S^2-mu_3/2.                            (4.1)

Hence tau is rational in R,S. Reconstruction gives

    k(U)=k(R,S).

This is a rational one-variable extension of k(R), so the geometric generic fibre is integral. UAGEOMFEED1 gives

    ker Delta on k[U]=k[l_1(K),l_2(K)].                  (4.2)

Combining Sections 1-3, for ANY polynomial profiles g_ij in ALL auxiliary coordinates and all signals, ANY block count N>=1, ANY finite active degrees, and ANY nonconstant scalar b(C), the map (2.2) is Keller exactly when every profile is polynomial in the two invariant auxiliary outputs and the retained signals. It then has the reverse polynomial factorization (3.2) with degrees3 and product_i d_i. It is never self-atomic. Constant extensions and identity stabilizations preserve these explicit factors. Independent controls cannot atomicify it, and changing only the constant output direction cannot do so either.

For the original canonical filter direction mu=(1/3,-1,2/3), take R1=K2+3*K1 and R2=K3-2*K1. Then (4.2) agrees with UACUBICKER1. When b(C)=1-lambda_N(C), N>=2, the profiles are S_N covariant and

    h_ij(0,0,rho_1,...,rho_1)=0,

the exact complete marked fibre of UACANFIB1 is preserved. Scheme-theoretically the marked output ideal imposes R=0 and all sigma_i=rho_1, hence all coefficient profiles vanish in that quotient. The old trace-filter equations are recovered exactly. Every such canonical member still has the reverse factors. No assertion is made that arbitrary different directions retain the old marked fibre.

## 5. The exact gate for a different polynomial source selector

For arbitrary beta in k[P], keep the same profiles and auxiliary block and put

    W(P,U)=A^-1*V*(Delta g).

The Schur identity (2.3) still holds, without assuming beta=b(C):

    det DT_g=1+W(beta).                                  (5.1)

The response at P=0 still vanishes, so every constant determinant is one. Expand the polynomial vector field in the independent auxiliary coordinates,

    W(P,U)=sum_gamma U^gamma*W_gamma(P).

This is a finite sum of polynomial k-derivations of k[P]. Therefore the EXACT condition for a polynomial source selector to keep the map Keller is

    beta belongs to intersection_gamma ker(W_gamma).    (5.2)

All retained signals belong to this common kernel, because d sigma*W=0. These automatic invariants do not solve the problem. If beta=b0(sigma), there are FORWARD factors

    H0(P,U)=(Phi(P;g(U,sigma)),U),
    Q0(Pout,U)=(Pout,G(U)-b0(first_signals(Pout))*mu),
    T_g=Q0 composed with H0.                             (5.3)

Both are polynomial Keller self-maps. Their degrees are product_i d_i and deg(G); for the auxiliary cubic both are nontrivial. This is also an exact counter-boundary to extending the b(C) classification to all selectors: signal-only beta can make T_g Keller even when Delta g is nonzero, but the map then factors in the other order.

For the cubic, a possible atomic candidate with a different selector must have W nonzero (otherwise Sections 3-4 give reverse factors), and must find an element of the common kernel (5.2) outside k[sigma]. A canonical filtering selector must additionally distinguish the required configurations at their shared diagonal signal value and give the prescribed selector values there. A function only of sigma is constant on that entire unfiltered fibre, so it cannot perform this selection. These necessary conditions do not assert that such a separating invariant exists, or that it would imply absolute atomicity.

This replaces an unstructured search for another beta by a common-first-integral problem with exact interpolation conditions. It does not classify all polynomial selectors, compute the common kernel in general, or construct a new canonical atomic map.

## 6. Consequences and remaining work

The earlier suggested multiple-control route requires a correction: any number of controls into the compiler are now covered when the auxiliary return is the single fixed direction beta(P)*mu. Multiple controls by themselves are insufficient. Changing the cubic return direction is insufficient too.

A candidate using this compiler must instead change the return geometry (for example a different source selector outside b(C), with a new response proof, or at least two independent return functions/directions), permit profiles outside the polynomial retained-signal class with a separate whole-source landing proof, or leave the compiler. Replacing the auxiliary block while keeping rank-one return faces the exact generic-line/closed-subring gate in Section 3. A formal nilpotent feedback matrix still is not a realizable polynomial map or an atomicity certificate. These necessary departures are not constructions, and a particular rank-two ansatz can still factor for other reasons.

The canonical non-atomic bundle theorem and the arithmetic Weil atomicity theorem remain valid. This package neither gives an absolutely atomic canonical replacement nor improves the canonical degree/dimension bounds. Its general proofs are the arguments above; exact controls verify selected identities and actual maps, not every geometric or all-degree statement. No publication novelty or AI-assisted external breakthrough is claimed.
