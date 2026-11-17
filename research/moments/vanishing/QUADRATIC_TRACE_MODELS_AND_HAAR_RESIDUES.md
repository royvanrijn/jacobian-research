# Quadratic trace models and exact Haar regularity in every degree

15 November 2026. Internal written proofs; no novelty, independent
specialist replay, external review or formal verification is claimed.
The [checker](verify_quadratic_trace_models.py),
[frozen polynomial images](quadratic_trace_polynomial_images.json),
[receipt](quadratic_trace_verification.json), and
[source audit](quadratic_trace_source_audit.json) delimit the evidence.
The universal torsor, topology and trace arguments are written proofs;
the checker supplies exact algebraic and nonzero moment controls.

The central result is an **equivalent all-degree model** for the entire
coordinate algebras of both remaining exceptional quotients. Each is
an explicitly constrained subalgebra of a rank-two quadratic extension
of C[SL_2]. Its Haar functional is exactly ordinary SL_2 Haar after
normalized quadratic trace. Derivative divisibility describes precisely
which expressions and markers are regular in the original quotient.
There is no coefficient-degree or support restriction.

This resolves the transfer and regularity problem for this model. It
does not decide whether its constrained trace kernel is MZ. Forgetting
the constraints would give a different algebra, already containing
known SL_2 counterexamples. Both exceptional MZ questions remain OPEN.

## 1. The additive quotient and every higher derivative image

Let Gamma be any finite subgroup of SL_2, and put

    A_Gamma=(C[x,y,z,w]/(xw-yz-1))^(right Gamma),
    D=z partial_x+w partial_y,
    E=x partial_z+y partial_w,
    B_Gamma=C[z,w]^Gamma, m=B_Gamma,positive.

Left D and E commute with right Gamma. Localizing at z or w gives
a rational additive slice x/z or y/w.
Their difference is1/(zw), and intersecting their kernel rings proves

    ker(D|A_Gamma)=B_Gamma.                             (1)

Write B_Gamma,n for ordinary bottom-row degree n. Complete reducibility
for left SL_2 and its elementary ladder formulas give, for every k>=0,

    D^k(A_Gamma) intersect B_Gamma
       = direct_sum_(n>=k) B_Gamma,n.                  (2)

Indeed D^k a in ker D implies a in ker D^(k+1). In each irreducible
V_n the bottom line can be reached in k steps exactly when k<=n.
Conversely, for a homogeneous q in B_Gamma,n,

    D^k E^k q = k! n!/(n-k)! *q,  0<=k<=n.             (3)

This proves every order in (2), rather than checking a finite prefix.
For k=1 it says that the plinth ideal is exactly m, so the additive
quotient has base

    U_Gamma=Spec B_Gamma minus {m},
    X_Gamma=SL_2/Gamma -> U_Gamma.                       (4)

It is a principal additive-group torsor. This can be seen by descending
SL_2->C^2 minus0 through the free finite Gamma action; a nonidentity
finite determinant-one matrix has no eigenvalue1. Alternatively, if
q=Da belongs to ker D, the slice a/q trivializes the action on q!=0.
The usual finite Taylor expansion proves A_Gamma,q=B_Gamma,q[a/q].
For the general local-slice terminology see
[Dufresne--Kraft, Definition5.1 and Proposition5.2](https://arxiv.org/pdf/1402.4299v1).
The standard punctured-plane SL_2 torsor is also discussed in
[Dubouloz--Heden--Kishimoto, Example2.1](https://arxiv.org/pdf/1707.08768v1).
No extension-classification theorem from that paper is needed here.

For later use, let f,g be homogeneous invariants of degrees n,m whose
only common zero is the origin. Define

    a_f=(x f_z+y f_w)/n, a_g=(x g_z+y g_w)/m,
    c={f,g}/(nm), {f,g}=f_z g_w-f_w g_z.

Then

    Da_f=f, Da_g=g, a_f*g-a_g*f=c.                      (5)

The last equality follows from Euler's identities and xw-yz=1.
The local slices s_f=a_f/f and s_g=a_g/g consequently have transition

    s_f-s_g=c/(fg).                                    (6)

Thus the full regular algebra is the exact chart intersection

    A_Gamma=B_Gamma,f[s_f] intersect B_Gamma,g[s_g],
    s_g=s_f-c/(fg),                                    (7)

inside its common function field. All multivariable regularizations
must pass this condition; merely clearing one displayed denominator
does not replace it.

## 2. The two complete Klein base algebras

For2O take

    f=z^8+14z^4w^4+w^8, n=8,
    g=z^2w^2(z^4-w^4)^2, m=12.

For2I take

    f=zw(z^10+11z^5w^5-w^10), n=12,
    g=z^20-228z^15w^5+494z^10w^10+228z^5w^15+w^20, m=20.

With c normalized as in (5), direct exact differentiation gives

    c^2=K_Gamma(f,g),
    K_O(F,G)=F^3 G/36-3G^3,
    K_I(F,G)=12F^5+G^3/144.                            (8)

The checker verifies f,g under exact finite-group generators, their
coprimality and absence of a common projective zero, and (8).
The degrees of c are18 and30. The second icosahedral form is the
normalized Hessian of f. Jacobian and Hessian covariance give the
invariance of c and g as well.

These are **full** invariant rings:

    B_Gamma=C[F,G,C]/(C^2-K_Gamma(F,G)).                 (9)

Here is a complete generation argument. The nonzero Jacobian proves
f,g algebraically independent. Scalar multiplication (z,w)->i(z,w)
fixes f,g and negates c, so1,c are independent over C(f,g). Hence
the displayed subring is free of rank2 over C[f,g], with Hilbert series

    (1+t^(n+m-2))/((1-t^n)(1-t^m)).                     (10)

This equals the exact all-degree Molien series already certified in
the [two-model reduction](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md).
The graded inclusion therefore has equal dimension in every degree
and is equality. This uses a rational character-series identity, not
a finite list of checked dimensions. In both cases nm=2|Gamma|.

## 3. A regular birational map from a quadratic cover

Use a separate copy of SL_2 with coordinate ring

    A_*=C[X,Y,F,G]/(XG-YF-1), delta=F partial_X+G partial_Y.

Form the affine quadratic cover

    R_Gamma=A_*[eta]/(eta^2-K_Gamma(F,G)),
    Z_Gamma=Spec R_Gamma,
    p:Z_Gamma -> SL_2.                                 (11)

The projection to (F,G,eta) makes Z_Gamma an additive torsor over the
same U_Gamma as (4). Its local slices S_f=X/F and S_g=Y/G have
transition1/(FG). In particular Z_Gamma is smooth, and delta is
locally nilpotent with delta eta=0. The cover p is finite flat of
degree2, branched at eta=0; it is not étale.

Multiplying each fibre coordinate by eta changes the transition
1/(FG) into eta/(FG), exactly (6). Therefore there is a **regular**
birational morphism over U_Gamma

    Psi:Z_Gamma -> X_Gamma,
    Psi^*f=F, Psi^*g=G, Psi^*c=eta,
    Psi^*a_f=eta X, Psi^*a_g=eta Y.                    (12)

On F!=0 it is s_f=eta S_f, and similarly on G!=0. These local formulas
glue even at eta=0. The map is an isomorphism when eta!=0 and contracts
the additive fibre when eta=0. No claim that it is a finite cover or
surjective is made. Its pullback on coordinate rings is injective.

Both projections to U_Gamma are locally trivial bundles with
contractible complex-affine-line fibre. They are homotopy equivalences.
Since Psi is over the identity of U_Gamma, Psi itself is a homotopy
equivalence. This supplies the topological information needed below;
birationality alone would not supply it.

For comparison, on eta!=0 the inverse of Psi followed by p gives the
covariant rational map

    X_Gamma -->> SL_2,
    [[a_f/c,a_g/c],[f,g]].                              (13)

Its determinant is1. It is generically of degree2, with poles at c=0.
Equation (12), not an unsupported change of contour through (13), is
the regular map used in the proof.

## 4. The all-degree quadratic Haar identity

Let h_Gamma denote normalized Haar on X_Gamma and h_* on the target
SL_2. Every Psi^*a has a unique expression a_0+eta a_1 with a_i in A_*.
Then

    h_Gamma(a)=h_*(a_0)
       =h_*((1/2)Tr_(R_Gamma/A_*) Psi^*a).             (14)

This is the essential transfer theorem. Here is its proof, including
the branch locus.

Use the same invariant algebraic volume form Omega_0 on each SL_2,
written on the z!=0 slice chart as d(x/z) wedge dz wedge dw. Let
Omega_Gamma be the descended form on X_Gamma; its compact period is
nonzero. On the f-chart it is

    Omega_Gamma=ds_f wedge df wedge dg/(nm c).

This follows from df wedge dg=nm c dz wedge dw and the fact that
s_f-x/z depends only on the bottom row. Using (12),

    Psi^*Omega_Gamma=(1/(nm)) p^*Omega_0.               (15)

The calculation is immediate on eta!=0; both sides are regular forms
on all of Z_Gamma, so it holds globally. In particular p^*[Omega_0]
is nonzero: Psi is a homotopy equivalence and [Omega_Gamma]!=0.

The deck involution tau sends eta to-eta and fixes A_*. Its action
on H^3(Z_Gamma;C) is the identity. To see this, its action on U_Gamma
is precisely the one induced by scalar multiplication (z,w)->i(z,w):
n,m are divisible by4 and n+m-2 is2 modulo4. The complex scalar
action has a path from1 to i and is a homotopy to the identity.
The projection Z_Gamma->U_Gamma is a homotopy equivalence and
intertwines the two involutions, proving the assertion.

Consequently every anti-invariant closed3-form on Z_Gamma has zero
cohomology class. In particular

    [eta a_1 p^*Omega_0]=0.

On each reductive homogeneous space, integration on the compact core
identifies [a Omega]=h(a)[Omega] in its one-dimensional H^3, as in the
[intrinsic Reynolds theorem](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md).
Apply this on X_Gamma and on the target SL_2, pull back through (15),
and discard the preceding anti-invariant class. The nonzero common
class p^*[Omega_0] then gives (14). This proof does not assert that
Z_Gamma is reductive homogeneous, or invoke an étale trace theorem
for the ramified map p.

For P,Q regular downstairs, write their images as p_0+eta p_1 and
q_0+eta q_1. Put

    E_m=sum_j binom(m,2j) p_0^(m-2j) K^j p_1^(2j),
    O_m=sum_j binom(m,2j+1) p_0^(m-2j-1) K^j p_1^(2j+1).

The exact marked identities in every degree are

    h_Gamma(P^m)=h_*(E_m),
    h_Gamma(QP^m)=h_*(q_0 E_m+K q_1 O_m).              (16)

The sums run over nonnegative exponents. These are identities for
every m, not extrapolations of the checker's finite controls.

As a nonzero check on both exceptional models, the checker expresses
the **entire** highest form f(x,y) as a regular even polynomial in
R_Gamma (20 terms for2O,40 for2I in the retained normal form). It
checks all8 or12 derivative conditions below and obtains

    h_Gamma(f(x,y) f(z,w))=8/15 on2O, 25/1092 on2I

both directly and through (14). A separate trivial-group quadratic
model verifies the normalization on monomials. These checks support,
but do not replace, the written cohomological proof.

## 5. Exact finite regularity conditions for each expression

The image of Psi^* is precisely

    S_Gamma={a in R_Gamma:delta^j a in eta^j R_Gamma
                                for every j>=1}.       (17)

Because delta is locally nilpotent, only finitely many j must be
checked for each a. This is an exact all-degree description of the
original algebra, not a new ambient algebra in which markers may be
chosen freely.

To prove (17), work on F!=0, where R_Gamma,F=B_Gamma,F[S_f]
and delta=partial_(S_f). The conditions imply that each coefficient
of S_f^j is divisible by eta^j: take the constant term of the jth
derivative. Thus the polynomial belongs to B_Gamma,F[eta S_f].
The converse is immediate by differentiation. The same statement
on G!=0 gives B_Gamma,G[eta S_g]. These are exactly the two regular
charts of X_Gamma in (7), so they glue to its coordinate algebra.
The opens F!=0 and G!=0 cover Z_Gamma; divisibility can therefore
be checked locally or globally.

For a=a_0+eta a_1, (17) is equivalently the pair of tests inside A_*:

    delta^j a_0 divisible by K^ceil(j/2),
    delta^j a_1 divisible by K^floor(j/2), every j>=1.   (18)

Both the phase and **every marker** in (16) must satisfy (18).
The complete exceptional problem is therefore equivalent to deciding
MZ for normalized quadratic Haar trace on the explicit algebra
S_Gamma, for the two polynomials K in (8). No group averaging or
unbounded search through right-invariant representations is required
to state or test membership in this model.

The ambient R_Gamma contains A_* and its known bad SL_2 phases;
that does not make them members of S_Gamma. For example neither the
old Hopf nor normalizer circuit, nor any nonconstant unary rational
reparameterization of either, passes (17) for these covariant maps.
Indeed a regular rational function of either circuit in R_Gamma is
already in A_*: R_Gamma is integral over the normal ring A_* and
intersects Frac(A_*) there. The constant-unit argument then makes
the reparameterization polynomial.

On the target additive chart X=Fs,Y=Gs-1/F, the leading terms of the
circuits are

    P_H=-16 F^4G^4 s^6+lower s-powers,
    P_N=4096 F^12G^12 s^20+lower s-powers.              (19)

The cusp component F^3-108G^2 of K_O=0, or G^3+1728F^5 of K_I=0,
has F,G both nonzero generically. Thus each circuit restricts to a
nonconstant polynomial in s on that component. For any nonconstant
polynomial R, delta R(P)=partial_s R(P) remains nonzero there. It
is not divisible by K, already violating the j=1 even test in (18).
This closes direct unary pullback through the new rational covariant
maps; it does not close multivariable expressions satisfying (18).

## 6. The poles are topologically essential, and the residue is explicit

The rational-slice obstruction has a period counterpart that applies
to the entire algebra. Let s be any rational slice Ds=1 on X_Gamma,
and let V be its regular locus. The invariant volume satisfies
Lie_D Omega_Gamma=0. Cartan's formula gives on V

    d(s i_D Omega_Gamma)=Omega_Gamma.                  (20)

Since [Omega_Gamma] is the nonzero generator of H^3(X_Gamma;C),
the inclusion H_3(V;C)->H_3(X_Gamma;C) is zero. No compact cycle
representing the Haar class can avoid the poles of s. This is why
merely deforming a contour away from a rational slice's poles cannot
supply the missing change-of-variables theorem. Equation (20) also
holds on any affine reductive homogeneous space with orientable compact
core,an invariant volume,a volume-preserving unipotent derivation and a rational slice;
the same nonzero-top-period conclusion follows in its dimension.

There is nevertheless a direct exact residue description. On the
ambient z*w!=0 chart put r=xw, so

    x=r/w, y=(r-1)/z,
    A[1/(zw)]=C[r,z^(+-1),w^(+-1)].

A finite Laurent expression sum_(i,j) a_(ij)(r) z^i w^j comes from
a regular function on SL_2 exactly when

    r^max(0,-j) (r-1)^max(0,-i) divides a_(ij)(r)       (21)

for every pair i,j. The four coordinate generators prove necessity
by their two-torus weights, and the same generators reconstruct each
term satisfying the divisibility, proving sufficiency.

For every regular f, including every invariant phase and marker,

    h(f)=CT_(z,w) integral_0^1
              f(r/w,(r-1)/z,z,w) dr.                  (22)

On x^a y^b z^c w^d, the constant term vanishes unless c=b,d=a;
in that case the integral is

    integral_0^1 r^a(r-1)^b dr
       =(-1)^b a! b!/(a+b+1)!,

the normalized SU_2 Haar monomial rule. This proves (22) for every
regular function, not just the finite monomial controls.

Equivalently, (22) is the double residue of the additive-torsor
transgression. Integrate f in the z-chart slice s=x/z and in the
w-chart slice t=y/w with zero integration constant at each local
slice origin. Their primitives differ by

    integral_0^(1/(zw)) f(zs,ws-1/z,z,w) ds.

Its z^(-1)w^(-1) coefficient is (22). The interval-and-constant-term
formula is this algebraic transgression; it is not a closed compact
Haar cycle in the pole-free slice chart, which (20) forbids.

## 7. Filtration control and the remaining mathematical question

There is one further all-degree safeguard. For F_k A_Gamma=ker D^(k+1),
the leading map a->D^k a/k! identifies

    gr_D A_Gamma = direct_sum_(k>=0) B_Gamma,>=k *T^k. (23)

Multiplication follows from Leibniz, and (2) proves surjectivity in
every degree. Homogeneous generators of B_Gamma and their finite
raising ladders generate A_Gamma and this associated graded algebra.
The Rees algebra is finitely generated and torsion-free over C[lambda],
hence flat, with general fibre A_Gamma and special fibre (23).

The induced special-fibre functional is an augmentation: its kernel
is an ideal, and hence MZ. This does **not** imply MZ upstairs. For
Gamma trivial the actual SL_2 Haar kernel has the known Hopf failure,
while its special-fibre kernel is that augmentation ideal. Thus this
natural flat degeneration does not justify positive MZ transfer from
the leading algebra. The lower filtration terms carrying the torsor
and residue cannot be discarded.

The open question is now an explicit all-degree one: do (16) and
(18), for K_O or K_I, admit a pure-null admissible phase and an
admissible marker with infinitely many nonzero marked moments, or
do the derivative constraints force eventual safety? The unrestricted
quadratic ambient algebra, a rational coordinate change, a leading
degeneration, or a finite moment prefix does not answer it.

No such exceptional phase/marker pair or general positive theorem is
proved here. The new result is the exact global algebra and functional
model, with its branch and regularity constraints retained. It provides
a concrete all-degree target for the remaining construction or safety
proof and does not propose a new support or degree sweep.
