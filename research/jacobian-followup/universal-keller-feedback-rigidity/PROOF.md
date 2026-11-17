# Output-controlled products, arithmetic atoms and the complete scalar-feedback obstruction

15 November2026. Written all-degree proofs with narrow exact controls. No independent whole-proof replay, formal verification, external review or publication-priority certificate. OP-KDESC remains open for atomic canonical realization. This package closes two large coupling classes and gives a base-sensitive atomicity theorem; it does not produce an absolutely atomic canonical replacement.

Work over a characteristic-zero field k. Self-atomicity concerns dominant regular self-map factors of the SAME affine space; k-stable means after every affine identity stabilization over k. Absolute requires every extension of constants too. These are distinct from normal-affine atomicity through different intermediate varieties.

## 1. UAPARALLEL1: all output-signal-controlled products have polynomial factors

For block i with source P_i=(x_i,y_i,z_i), put

    t_i=1+x_i*y_i, q_i=t_i^2*z_i+y_i^2*(1+3*t_i),
    sigma_i=t_i*q_i.

Use the normalized compiler

    Phi_i=(sigma_i,
      -y_i/2-3*x_i*q_i/2
        -sum_(j>=4) (j/2)*a_ij*t_i^2*x_i^(j-2)*q_i^j,
      x_i*(5-3*t_i)-x_i^3*z_i
        -sum_(j>=4) (j-2)*a_ij*(x_i*q_i)^j).             (1.1)

Choose ANY finitely many polynomial profiles

    a_ij=f_ij(sigma_1,...,sigma_r,Z),                      (1.2)

where Z are optional promoted identity coordinates. Set F=(Phi_1,...,Phi_r,Z). Let d_i be the largest j with nonzero f_ij, or3 when all higher profiles vanish. Then:

1. F is polynomial on the whole source with determinant1 and generic degree product d_i.
2. Every subset of blocks gives an actual polynomial Keller factorization.
3. Arithmetic and geometric monodromy are product_i S_(d_i), in the PRODUCT action on product_i{1,...,d_i}.
4. All embedded intermediate fields form the Boolean lattice on the r blocks. Every node has a polynomial model. Every maximal polynomial decomposition has length r and degree multiset {d_1,...,d_r}, including after arbitrary constant extensions and identity stabilization.

Thus arbitrary polynomial cross-coupling through the retained first outputs does not make an absolutely atomic map when r>1. This is stronger than an obstruction based only on imprimitive monodromy: here the actual polynomial factors are constructed.

### Jacobian and the partial-processing formula

For fixed coefficient parameters, each compiler Jacobian is1. Its coefficient-response columns V_i have first output entry zero, while the first row of its source Jacobian A_i is d sigma_i. Hence

    d sigma_i * A_i^-1 * V_i=0.                            (1.3)

For the full profile substitution the determinant lemma applies to the block diagonal A and the coefficient response. Every entry of the small contraction matrix is zero by(1.3), so det DF=1. This extends UAPROFILE1's one-block argument to all cross profiles simultaneously.

For a subset S, define P_S by processing blocks in S with(1.1)-(1.2) and leaving all other blocks unchanged. To process a remaining block after this step, recover each original signal by

    sigma_i^S=current first coordinate    if i in S;
              sigma_i(current block)     if i not in S.    (1.4)

Use f_ij(sigma^S,Z) for the remaining blocks. This defines a polynomial Q_S, with F=Q_S P_S. Both factors are Keller by the same determinant argument, and their degrees are products of the d_i on the processed and remaining subsets. Any cycle in the profile-dependence graph is immaterial: the original signal of a processed block is still available as its first output.

### Monodromy, Boolean fields and maximal polynomial chains

Over K=k(sigma_1,...,sigma_r,Z), the inverse equations are

    E_i(S_i)=S_i+B_i*S_i^2+sigma_i*S_i^3
             +sum_j f_ij(sigma,Z)*sigma_i^j*S_i^j-C_i/2.  (1.5)

Each has the independent parameter pair B_i,C_i and nonzero linear coefficient. UAQTILT1 and UAALLATOM1 give S_(d_i), geometrically over K. The regular splitting fields in disjoint independent parameter pairs are linearly disjoint: their geometrically integral Galois covers have geometrically integral products. Their compositum therefore has group product S_(d_i), and reconstruction identifies its product root action with the original cover. This also proves the degree assertion. No full wreath product is inferred from a composition.

A point stabilizer is H=product S_(d_i-1). Any overgroup J between H and G projects in coordinate i either to the maximal point stabilizer or to S_(d_i). In the latter case J intersect S_(d_i) is normal in the full projection and contains the point stabilizer. Its normal closure is S_(d_i), since d_i>=3. Thus J contains that whole coordinate factor. This proves the Boolean interval, without classifying all subgroups of G.

The models P_S realize every node, and(1.4) realizes every nested-subset transition. KRP1's affine-open rigidity applies to every embedded field already supplied with such a model: any other compatible polynomial model differs by an intermediate automorphism. Thus a polynomial chain skipping a Boolean rank can be refined, and every maximal chain has r factors. No assertion is made that all order-dependent factors have the same stable class as a fixed coefficient-free seed, or that the state-dependent factors are literally commuting endomorphisms. Those stronger separated-product statements are not needed. The construction and interval persist under constants and identity variables.

### A simultaneous canonical-fibre obstruction in this family

If r=N>1, profiles are S_N covariant and Z is fixed at the mark, an S_N-fixed target has all block coordinates equal. Its signal vector is diagonal, so all coefficient profiles specialize to the same constants in every block. Its complete fibre is therefore T^N for one scalar finite etale fibre T. If nonempty, it contains the diagonal T of S_N-fixed points. It cannot be the nontrivial transitive natural N-point S_N set. This closes the canonical complete-fibre route within these output-controlled tensor products, in addition to their absolute-atomicity failure. It does not exclude source-dependent feedback outside this family.

## 2. UAWEILATOM1: exact base-field atomicity of restriction of scalars

Let f_d:A^3->A^3 be the fixed normalized constant-profile compiler of degree d>=3: the cubic for d3, and a_d=1 with other higher coefficients zero for d>=4. It is defined over Q and has geometric and arithmetic S_d by UAALLATOM1. Let A/k be ANY finite etale algebra of rank N and write

    A=A_1 times ... times A_c,    each A_i a field.

Evaluating f_d with the algebra operations gives a polynomial Keller self-map W_(A,d) of the k-vector space A^3. Then:

- Its degree is d^N and geometric monodromy is S_d^N in product action.
- If L is the splitting field of A and Q=Gal(L/k) in its action on the N embeddings, arithmetic monodromy is S_d^N semidirect Q in product action; the constant field of its normal closure is exactly L.
- Its arithmetic intermediate-field lattice is Boolean on the c field components. Every maximal k-polynomial decomposition has length c and degree multiset {d^[A_i:k]}. All component factors are k-stably self-atomic.
- Consequently W_(A,d) is k-stably self-atomic exactly when A is a FIELD. After any extension K/k the exact length becomes the number c(A tensor K) of its field components. Absolute stable atomicity holds exactly when N=1.

This produces nonfinite affine-space k-stable atoms in degrees d^N for every finite separable field input, including many non-prime-power degrees. They become decomposable over splitting constants. It does not construct an absolutely atomic map or assert the complete natural rank-N marked fibre for W_(A,d).

### Arithmetic closure and the subgroup interval

Over L the map is a separated product of N copies. Its regular disjoint splitting fields have group S_d^N and compositum M with constants L. Present each splitting field by a full ordered root list. Q acts on L and permutes the block index, keeping each ordered-root position fixed; this gives a semidirect action on the generic splitting-field presentation. The field of one selected root per block has stabilizer

    (S_(d-1))^N semidirect Q.                              (2.1)

Its fixed field is the original k(source) by ordinary Galois descent. The action on d^N root tuples is faithful, so the core of(2.1) is trivial: M is precisely the normal closure, not a larger auxiliary Galois extension. This proves the exact arithmetic and constant-field assertions.

An overgroup of(2.1) meets S_d^N in a Boolean subgroup from Section1, necessarily stable under Q. The subset of expanded coordinate factors must therefore be a union of Q orbits, which are exactly the field components A_i. Conversely every such union gives an overgroup. All are realized polynomially by processing the corresponding algebra idempotent components. KRP1 again makes every maximal polynomial decomposition field-maximal. When c=1 the arithmetic action is primitive, so the polynomial determinant chain rule and the birational self-map lemma in UASELFATOM1 give k-stable atomicity. This is the familiar primitive product-action wreath phenomenon; no new general permutation-group theorem is claimed.

After extension to K the same proof applies to A tensor K. Over an algebraic closure it has N components. In particular connectedness over k must not be confused with geometric connectedness, nor arithmetic primitivity with absolute atomicity. At the split cubic target, the first-stage fibre is the algebra-valued configuration scheme {p0,p1,p2}^N, usually not split over k even though p0,p1,p2 themselves are rational. This respects the split-fibre monodromy criterion UAUNIVMON1.

## 3. UACUBICKER1: the cubic has exactly two polynomial feedback invariants

Use the cubic K of UAFILTER1, with source U=(u,v,w) and mu=(1/3,-1,2/3). Its Jacobian is1, so

    D=DK,  delta=D^-1*mu                                  (3.1)

is a POLYNOMIAL vector field, viewed as a k-derivation. Set

    R1=K2+3*K1,    R2=K3-2*K1.

Then

    ker(delta on k(u,v,w))=k(R1,R2),
    ker(delta on k[u,v,w])=k[R1,R2].                       (3.2)

Indeed delta K=mu gives delta R1=delta R2=0. Put a=K1, b=K2, c=K3, and S=u/(1+uv) in the generic inverse equation. Replacing b=R1-3a and c=R2+2a gives

    S+R1*S^2-R2/2+a*(S^3-3*S^2-1)=0,
    a=(R2/2-S-R1*S^2)/(S^3-3*S^2-1).                     (3.3)

Reconstruction proves k(u,v,w)=k(R1,R2,S). Since delta a=1/3, delta is nonzero in this one-variable rational function field over k(R1,R2). Characteristic zero makes its rational constants exactly k(R1,R2).

For the polynomial assertion, the map (R1,R2):A^3->A^2 is smooth: it is K followed by a rank-two linear projection. Thus it is flat, by the classical smooth-flat theorem. A rational function p(R1,R2)/q(R1,R2) in reduced form cannot pull back regularly unless q is constant. Otherwise q(R1,R2) is a nonconstant nonunit on A^3. A component of its zero divisor dominates the corresponding base divisor: smoothness gives relative dimension one and an open image, so a divisor cannot be supported only over a base codimension-two subset. At its generic point p is nonzero and the pullback has a pole. This contradiction proves the polynomial kernel identity. This argument does not assume that a polynomial rational first integral automatically descends regularly.

Neither local nilpotence nor completeness of delta is asserted. In positive characteristic its kernel contains p-th powers;(3.2) is not exported there.

## 4. UASCALARFEED1: every scalar feedback that remains Keller is still decomposable

This classification permits arbitrary block number N>=2, arbitrary finite coefficient degrees, arbitrary polynomial dependence on ALL retained signals, and arbitrary polynomial dependence on one common scalar parameter. Use profiles

    a_ij=f_ij(a,sigma_1,...,sigma_N),                      (4.1)

and assume genuine a-dependence: at least one partial_a f_ij is nonzero. Let Phi(P;a) be the full output-controlled compiler from Section1, with a frozen as parameter. Write X_i for the first source coordinate in block i and set

    C(P)=sum_i (3*X_i^2-X_i)/2,
    beta(P)=1-lambda_N(C(P)),
    T_alpha(P,U)=(Phi(P;alpha(U)), K(U)-beta(P)*mu),        (4.2)

where alpha is ANY polynomial in the three auxiliary source coordinates. This looks like feedback in both directions. If its Jacobian determinant is any nonzero constant, that constant must be1: all parameter-response columns vanish at P=0, where the block determinant is1 for every U. Thus the following classifies ALL Keller members, not just one normalization. The exact theorem is

    det D T_alpha=1
      iff alpha is in k[R1,R2].                           (4.3)

Every such Keller map has an actual reverse-order polynomial self-factorization. In particular NONE is self-atomic, absolutely stably atomic, or an atomic replacement for UACANFIB1. This is an all-degree statement, not a finite ansatz computation.

### Schur complement and the nonzero loop response

Put A=D_P Phi at frozen a. Section1 gives det A=1 even with all the signal profiles. Set

    L(P,a)=d beta * A^-1 * partial_a Phi.

Both A^-1 and D^-1 are polynomial. The block determinant lemma gives the exact polynomial identity

    det D T_alpha=1+L(P,alpha(U))*delta(alpha)(U).          (4.4)

The key point is that L is NOT the zero polynomial for any genuinely nonconstant profile family(4.1).

Choose an active block i and put every other block at X_j=Y_j=0, Z_j=pi_j, with pi_j independent. All their parameter-response columns vanish because every higher term has X_j power at least two. Their contributions to L vanish. In the active block set Y_i=0, x=X_i, z=Z_i and w=xz; localize at x and regard x,w as independent. The other signals are independent constants, and the active signal is pi_i=z=w/x. Suppress the index i and let

    B(w)=sum_j partial_a f_j(a,w/x)*w^j,
    C0(w)=sum_j j*(j-2)*f_j(a,w/x)*w^j.

All other signal arguments remain present as independent parameters. Exact differentiation of the inverse equation, using d sigma*A^-1*Phi_a=0, gives the first source response

    (A^-1*Phi_a)_x=
       sum_j (j+1)*partial_a f_j(a,w/x)*w^j
       -3*x*w*B(w)-x^-1*C0(w)*B(w).                      (4.5)

Because the profiles are polynomial in the active signal, they have only NONPOSITIVE powers of x after pi_i=w/x. Let -m be the largest x power in the nonzero B. The middle term of(4.5) has largest power1-m; the first has power at most-m, and the last at most-1-m. Its nonzero leading coefficient cannot cancel. Thus the response is nonzero. Multiplication by -lambda_N'(L(x))*L'(x), with L(x)=(3x^2-x)/2, remains nonzero in characteristic zero. This proves the assertion about the loop response L(P,a).

If alpha is nonconstant, substitution a=alpha(U) is injective in the polynomial ring with the independent P variables. Therefore L(P,alpha(U)) is nonzero. Equation(4.4) makes Jacobian one equivalent to delta(alpha)=0. UACUBICKER1 then gives alpha=h(R1,R2). Constant alpha already has that form. This proves(4.3), including constants.

### Every allowed control is recoverable from the other output

Let V=K(U)-beta(P)*mu. Since R1,R2 are linear combinations annihilating mu,

    R1(U)=V2+3*V1,    R2(U)=V3-2*V1.

Hence if alpha=h(R1,R2), define

    H(P,U)=(P,K(U)-beta(P)*mu),
    Q(P,V)=(Phi(P;h(V2+3V1,V3-2V1)),V).                 (4.6)

Then T_alpha=Q H. Both are polynomial Keller self-maps, the first of degree3 and the second of degree product d_i, where d_i are the effective output-controlled block degrees after substitution. Every d_i>=3. Thus both factors are noninvertible. This conclusion does not rely solely on a field-theoretic monodromy test.

If the profiles are S_N covariant, h(0,0)=a0 and

    f_ij(a0,rho_1,...,rho_1)=0  for every i,j,

the complete natural marked rank-N fibre from UACANFIB1 is retained: the target auxiliary coordinate V=0 forces a=a0, and first outputs force the diagonal signal vector, so all higher coefficients vanish on that fibre. These are genuine canonical feedback families, but every Keller member merely reverses the old factor order. A nonzero new profile generally increases generic degree; it does not remove the nontrivial intermediate map.

## 5. What remains outside the closed coupling classes

A general two-block feedback with polynomial Keller families F(x,a),G(y,b), a=alpha(y), b=beta(x), has

    det J =det(I_r-M),
    M=D_y alpha*(D_y G)^-1*D_b G*
      D_x beta*(D_x F)^-1*D_a F.                          (5.1)

This is a small r-by-r matrix, where r is the number of independent a controls; all inverses are polynomial when the frozen determinants are1. In one control it is the scalar product above. With two or more controls det(I-M)=1 admits cancellations and nilpotent matrices with nonzero entries; a formal matrix solution need not arise from actual whole-source Keller families.

The next atomic construction must leave the output-recoverable block class and the complete scalar-feedback class(4.1)-(4.2): for example multiple genuinely independent controls satisfying(5.1), a different auxiliary geometry, or a different compiler. No realization or primitive monodromy certificate for such a multi-control candidate is claimed. A shared cyclic graph by itself is not such a certificate. Positive-characteristic feedback with zero derivatives is outside these characteristic-zero kernel arguments.

The canonical non-atomic construction remains valid. The new arithmetic Weil atoms do not upgrade it: they realize the input cover only as a selected subconfiguration before filtering, and lose atomicity upon splitting constants when N>1. Full marked fibre, arithmetic atomicity and absolute atomicity remain separate gates. All general proofs are written above; controls verify selected exact identities and finite actions, not the universal group/regularity arguments or publication priority.
