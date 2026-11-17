# Ramified branch flattening and global Haar periods

15 November 2026. Internal written proofs with exact symbolic controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_ramified_branch_flattening.py),
[receipt](ramified_branch_verification.json) and
[source audit](ramified_branch_source_audit.json).

Every regular exceptional section can be made arbitrarily flat along
the quadratic branch by an explicit regular nonlinear map. The full
Jacobian-weighted bad tail survives. However, normalized Haar is not
continuous for the branch-adic topology, or for the topology of any
nonzero proper spatial ideal on any connected reductive homogeneous
space. The formal limit loses an essential pole divisor. These facts
explain exactly why arbitrarily accurate branch jets do not supply
unweighted purity. Both exceptional MZ decisions remain **OPEN**.

## 1. A general obstruction to spatial-adic Haar extension

Let X=G/H be a connected affine complex reductive homogeneous space,
A=C[X] and h its normalized Haar functional. Its bilinear pairing

    (a,b) -> h(ab)                                           (1)

is nondegenerate. Choose compatible compact forms K,L. On K/L the
complex conjugate of a regular a extends to a regular a^dagger in A:
on K it is a coefficient of the dual representation, and right
L-invariance extends to right H-invariance by Zariski density.
Since K/L is Zariski dense in X, for0!=a,

    h(a*a^dagger)=integral_(K/L) |a|^2>0.                    (2)

For every nonzero ideal I of A and every N>=1, therefore,

    h(I^N)=C.                                               (3)

Indeed A is a domain, so choose0!=a in I^N. Nondegeneracy supplies b
with h(ab)!=0, and scalar rescaling supplies any value. In particular
ker h contains no nonzero ideal.

For any nonzero proper I there is **no continuous extension** of h
to the I-adic completion of A, with C given its usual topology. Choose
a_N in I^N with h(a_N)=1. Then a_N tends to0 in that completion, but
its asserted period does not. The same obstruction holds with C
discrete. This statement is independent of whether ker h is MZ: it
also applies to tori, where the MZ property is known to hold.
For example, on C[z,z^(-1)] take I=(z-1) and
a_N=(z-1)^N z^(-N). Its constant term is1 for every N, although
a_N belongs to I^N.

The ideal here is an ideal of functions on X. This does not invalidate
finite coefficient-family moment tests, degree-bounded calculations,
or adic manipulations accompanied by a separate global period
comparison. It forbids recovering the full normalized Haar functional
by continuous passage to a spatial formal neighbourhood.

## 2. An explicit semigroup of regular ramified ambient maps

Use the [quadratic exceptional model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md):

    R=C[X,Y,F,G,eta]/(XG-YF-1,eta^2-K),
    K_O=F^3 G/36-3G^3, K_I=12F^5+G^3/144,
    delta=F partial_X+G partial_Y,
    lambda_O=-1/3, lambda_I=144,
    E=lambda*eta^2*X^3, r=XG.

The following J is a polynomial, with FJ=E-1:

    J_O=Y*(r^2+r+1)-X^3 F^2 G/108,
    J_I=Y*(r^2+r+1)+1728X^3 F^4.                           (4)

For every integer k>=0 let S_k=1+E+...+E^(k-1), with S_0=0, and
define a regular map theta_k:Z->Z over the common punctured base by

    theta_k^*X=X E^k,
    theta_k^*Y=Y+XG J S_k,
    theta_k^*(F,G,eta)=(F,G,eta).                           (5)

The determinant identity follows directly from FJ=E-1:

    (X E^k)G-(Y+XG J S_k)F=1.

On the generic additive fibre, X maps to a nonzero constant times
X^(3k+1). Hence the generic field degree is3k+1, and theta_k is not
an automorphism when k>0.
They contract the additive fibres at eta=0, so they are not finite
for k>0. Also

    theta_k^*E=E^(3k+1),
    theta_k theta_l=theta_(k+l+3kl).                        (6)

The indices3k+1 multiply. These maps are equivariant for the inherited
original torus: E has weight0 because3m=2(n+m-2) in BOTH models.

Every theta_k induces the identity on H^3(Z;C), since it is over the
identity of the punctured base and both additive torsor projections
are homotopy equivalences. Its generic degree3k+1 must not be used
instead of that cohomology scalar1. With the nowhere vanishing
canonical form Omega_Z, direct differentiation gives

    theta_k^*Omega_Z=j_k Omega_Z,
    j_k=(3k+1) E^k.                                        (7)

Consequently for all a in R, with the two exact period functionals
T(a_0+eta*a_1)=h_*(a_1), L(a_0+eta*a_1)=h_*(a_0),

    T(j_k theta_k^*a)=T(a),
    L(j_k theta_k^*a)=L(a).                                (8)

The second identity uses theta_k^*eta=eta. This is weighted transfer,
not ordinary L-preservation. Thus it does not contradict the proved
rigidity of Haar-preserving **automorphisms** of R.

## 3. Exact scalar moments close the direct ambient Hopf route

For all integers a,b>=0,

    L(E^a r^b)=1/(3a+b+1).                                 (9)

Indeed expand (lambda K)^a X^(3a+b) G^b. Every term involving F
has zero ordinary Haar moment, since it has no matching Y factor.
The unique term containing no F is X^(3a+b)G^(3a+b), with coefficient1
in both models. Its Haar value is the stated fraction.

In particular L(theta_k^*r)=1/(3k+2), while L(r)=1/2. For the old
target Hopf phase P_H, inherited torus weights and the
[height criterion](HOPF_HEIGHT_CRITERION_AND_REAL_OBSTRUCTION.md) give

    L(theta_k^*P_H)
      =1-3*L((2r E^k-1)^2)
      =-2k(6k+1)/((3k+2)(2k+1)).                           (10)

This is nonzero for EVERY k>=1. Thus this entire ramified ambient
postcomposition family destroys the old phase's purity already at
the first moment. Its Jacobian-weighted pure identity remains valid
by (8), and cannot be substituted for (10).

These maps also fail exact admissibility preservation. The admissible
element eta X maps to lambda^k eta^(2k+1) X^(3k+1). Its derivative of
order3k+1 is

    lambda^k (3k+1)! eta^(2k+1) F^(3k+1),

which is not divisible by eta^(3k+1) for k>0: F is a unit at the
branch divisor and eta has order1 there. Thus the automorphism
admissibility theorem cannot be extended to these weighted maps.
None of these observations decides ordinary purity after pulling
back through an original regular section; its functional is different.

## 4. Arbitrarily flat regular sections in the full original algebra

Let sigma_b be ANY of the [regular Bezout sections](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md),
with coordinates U_b=U+Fb,V_b=V+Gb,C and

    H_b=H+Db, W_b=C H_b,
    h_Gamma(W_b sigma_b^*a)=L(a).

Pull back E,J,S_k through sigma_b, denoting the results e_b,j_b,s_(b,k).
The composed section theta_k sigma_b has parameter

    b[k]=b+U_b*j_b*s_(b,k),
    U_(b[k])=U_b e_b^k, V_(b[k])=V_b+G U_b j_b s_(b,k).     (11)

These are actual regular invariant polynomials for all k. Their
heights and normalized Jacobian weights satisfy the exact identities

    tau_(b[k])+1=(tau_b+1) e_b^k,
    W_(b[k])=(3k+1) e_b^k W_b.                             (12)

Because e_b=lambda C^2 U_b^3, (12) gives

    tau_(b[k])+1 in C^(2k) A_Gamma,
    W_(b[k]) in C^(2k+1) A_Gamma,
    h_Gamma(W_(b[k]))=1.                                  (13)

Thus every prescribed finite branch order can be attained by a
regular section, starting from ANY section. If b originally has
weight L, so does b[k]. No degree bound or enumerated search is used.

All these sections retain the exact normalized weighted transfer and
the fixed positive Hopf marked tail supplied by the Bezout theorem.
Their **unweighted** pure moments are not given by (8). In particular
the nonzero AMBIENT value in (10) is a W_b-weighted source moment,
not an asserted ordinary source moment. No exceptional counterexample
or purity theorem for the new regular sections is claimed.

## 5. The formal limit and the missing global divisor

In the C-adic completion of A_Gamma, e_b tends topologically to0
under its powers. Also F is a unit there: the identity

    1=F*(-j_b)+e_b                                         (14)

gives the geometric-series inverse of F. Formula (11) therefore
has the formal limit

    b[infinity]=b+U_b*j_b/(1-e_b)=-U/F,
    U_infinity=0, V_infinity=-1/F, tau_infinity=-1.           (15)

The parameter limit is independent of the initial b. Its Jacobian
weight is0. Nevertheless every finite weight in (13) has Haar
average1. This is an explicit witness for the general noncontinuity
theorem, not a failed numerical limit.

The rational b[infinity] is genuinely not regular on all of X:
if U/F were regular, (UG-VF)=1 would make F a unit. But F has a
nonempty zero divisor. The formal neighbourhood inverts F precisely
because F=0 and C=0 are disjoint on X: in the Klein base their common
zero would force G=0, the omitted origin.

The lost divisor is also topologically essential. The rational slice
a_f/F has D(a_f/F)=1 and is regular on F!=0. The established identity

    d((a_f/F) i_D Omega_Gamma)=Omega_Gamma

holds there. Since the Haar top class is nonzero, no representative
of its compact homology class can avoid F=0. The formal collapse
therefore removes a divisor that the actual global period must meet.

This does not prove that finite branch jets can never decide purity
under additional hypotheses. It proves the precise non-transfer:
arbitrary branch flatness and the resulting formal collapse do not
give a continuous extension of the actual Haar functional. Any
unweighted purity proof must still control the global period, including
the missing F=0 boundary. The complex-height, other homogeneous-weight,
mixed-section and full exceptional MZ questions remain open.
