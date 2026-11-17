# Profile-factor rigidity and the canonical plane completion

15 November 2026. Three written general claims and selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This extends [native profile orbits](../universal-keller-native-orbits/PROOF.md) using a classical plane line-injectivity theorem, rather than a further coefficient search.

Work over a characteristic-zero field k, with retained identity parameters u. All assertions about a relative block fix u. Put

    t=1+xy, q=t^2*z+y^2*(1+3t), sigma=t*q, X=x*q,
    K=(sigma,B0,C0),
    B0=-y/2-3*x*q/2, C0=x*(5-3t)-x^3*z.

A native map is K plus a finite polynomial source list times

    v_j=(0,-(j/2)*sigma^2*X^(j-2),-(j-2)*X^j), j>=4.

A retained profile has these coefficients in k[u,sigma]. UANATIVE2GEN1 reduces every native list to two unique polynomial effective controls. Here "degree" for a field extension is generic covering degree; degree in r and total polynomial degree are separate.

## 1. Published inputs and an attribution refinement

**HC4MYGL0, imported earlier elsewhere in this repository**, is Gwozdziewicz's [Injectivity on one line (1993)](https://arxiv.org/html/alg-geom/9305008v1), Theorem1.1: a plane Keller map over an algebraically closed characteristic-zero field that is injective on one affine line is an automorphism. Its Lemma2.3 covers a map fixing a line pointwise. The paper explicitly invokes Abhyankar--Moh rectification of a polynomially embedded affine line and Newton-polygon similarity. We use these established results as imports, not as new proofs of JC(2). Over a general characteristic-zero field, pass to an algebraic closure and descend the unique inverse. The primary HTML and the retained import record are pinned.

The characteristic-zero algebraically closed saturation assertion underlying UAQFSAT1 was already stated publicly in **Anton Trushin, [Contracted divisors and Degree-Two Maps (2026)](https://arxiv.org/html/2605.26390v1), Proposition3.1.1(2)**. That proposition assumes only no contracted divisors, which is weaker than quasi-finiteness. Its coprimality proof is the same elementary dimension mechanism. Our prior note already disclaimed novelty, but should be read with this explicit earlier attribution. The arbitrary-field quasi-finite proof in UAQFSAT1 remains valid; novelty of that scope is not asserted. We import no unrestricted degree-two classification from this preprint.

In particular every retained profile H with identities is etale and satisfies

    k[u,P] intersection k(u,H)=k[u,H].                   (1.1)

This is polynomial rational descent, not invertibility of H. UAPROFILE1 and UAALLATOM1 supply its whole-source polynomial landing and absolute stable atomicity, respectively.

## 2. UAPROFFACT1: every retained-profile first factor is already a coordinate explanation

Let H=(sigma,Bh,Ch) be ANY retained profile, and T=(sigma,Bout,Cout) ANY whole polynomial native Keller map, on the same source with identities u. No degree restriction or affine-r assumption is made on T. Then the following are equivalent:

1. Bout,Cout belong to the RATIONAL field k(u,H).
2. There is a polynomial target map E with T=E composed H, fixing u.
3. There is a UNIQUE h in k[u,a,c] such that

       T=(sigma,Bh+sigma^5*Ch^3*h(u,sigma,Ch),Ch).        (2.1)

Thus E is automatically the explicit target automorphism from UANATIVEORBIT1. It is impossible in this native class for a nontrivial polynomial target factor to follow any retained-profile H. Every recovered map has exactly H's generic degree and inherits its absolute stable atomicity. This is a statement about the specified literal profile first factor; arbitrary intermediate maps and source changes moving the signal field are not classified.

### Proof without a plane inverse conjecture

Assume1. Saturation (1.1), using UAQFSAT1, gives a polynomial E on the whole target, with T=E composed H. Dominance makes E retain u and a=sigma. The chain rule and det DH=det DT=1 give det DE=1. On x=0 all native maps have exactly

    (sigma,B,C)=(z+4*y^2,-y/2,0).

This source plane covers every target point (a,b,0), including all retained u. Therefore E fixes the entire target plane c=0 POINTWISE.

Over K0=k(u,a), its last two components form a polynomial plane Keller map fixing the line c=0 pointwise. HC4MYGL0 makes this generic plane map an automorphism. Coordinate irreducibility then forces E_c=gamma*c and E_b=b+f(c), with gamma in K0*. Both H and T have C-C0 in (x^4), while Ch=x*L with L modulo x equal2. Hence gamma=1.

Consequently T differs from H only by f(sigma,Ch). Since E was already a whole polynomial map, f is in k[u,a,c]. The exact target-coordinate landing ideal of UANATIVEORBIT1 is (a^5*c^3); it gives (2.1). Conversely (2.1) is an explicit polynomial factorization and rational recovery. Uniqueness follows from dominance.

This is the needed extra step in the earlier UACUBICDESC1 conditional analysis: **the native branch "target factor of degree greater than one" is empty**, because its line-fixed generic plane map is invertible. The earlier conditional factorization statement was valid; the stronger result removes the unresolved alternative.

### Consequences in every profile degree

- Adding rational recovery through SOME retained profile to the four equivalent conditions of UANATIVEORBIT1 gives another equivalent condition. In particular Cout genuinely nonlinear in r rules out EVERY literal retained-profile first factor, not just the original cubic.
- Two distinct retained profiles have incomparable literal output fields: an inclusion k(u,Hg) subset k(u,Hf) would put Hg in (2.1). Both have affine-r B outputs, whereas a nonzero permitted shear has c-degree at least3. The shear vanishes, hence Hg=Hf. This is literal field inclusion, not inequivalence under arbitrary source coordinates.
- For the original cubic, UACUBICDESC1's two nonconstant basis remainders vanish exactly for the cubic target-shear orbit (2.1), necessarily of generic degree3. The cubic test remains exact for arbitrary polynomial candidates; the empty native degree-greater-than-one branch is now closed.
- The UASELFFORMAL1 tangent has a C output with nonzero r dependence relative to the fixed coefficient S. If a finite polynomial native family integrated that tangent, its generic member could not lie in the classified profile orbit. The new theorem rules out ANY retained-profile first factor for that hypothetical generic member. No such finite integration is constructed, and other polynomial factors remain open.
- Inserting these recovered coordinate changes into a full canonical receiver preserves its existing factors. Single-block factor rigidity does not assert that mutually changing blocks are retained identities.

No coefficient-degree bound, arbitrary Keller-map inverse theorem, or assumption that a rational chart is the whole source occurs in this argument.

## 3. UARPLANE1: the nonnegative-r sector has a genuine polynomial plane completion

Pass only for proof to the generic signal field K0=k(u,a), a=sigma. The generic source leaf has coordinate ring

    R=K0[x,y,(1+xy)^(-1)].

Indeed a is nonzero in K0, so t is invertible and z is reconstructed from a. Put

    S=x/t, Y=y, r=1-S*Y=1/t.
    x=S/r, y=Y,
    z=a*r^3-r*(r+3)*Y^2.                               (3.1)

Then R=K0[S,Y,r^(-1)]. Moreover

    R intersection K0[S,S^(-1),r]=K0[S,Y].              (3.2)

To prove it, K0[S,S^(-1),r]=K0[S,S^(-1),Y]. The polynomials S and r=1-SY are comaximal: 1=r+SY. Thus D(S) and D(r) cover the ENTIRE affine (S,Y) plane. A rational function regular on both is a polynomial; alternatively cancel the coprime source denominators in the UFD K0[S,Y]. This is an exact ring identity.

Therefore a whole-source output has NO negative r powers in its finite Laurent expression precisely when its generic leaf map extends polynomially across r=0 to A2_(S,Y). No negative-r sector is silently included. The new boundary is the smooth irreducible hyperbola

    D=V(1-SY), isomorphic to Gm.

For T native and polynomial in r, write its plane completion F=(B,C). The chart Jacobian is B_S*C_r-B_r*C_S=-r/S, hence

    det D_(S,Y) F = r=1-SY,
    F(0,Y)=(-Y/2,0).                                   (3.3)

The unique two effective controls also lie in K0[S,Y]: the native response matrix depends on S,a alone, so absence of negative r powers of both outputs implies the same for the controls, and (3.2) applies. This describes an actual plane completion of the generic leaf, not an identification of the original leaf or whole threefold with affine plane.

The base completion is

    B0=-Y/2-3*a*S/2,
    C0=2*S-S^2*Y-a*S^3.                                (3.4)

Source jets become B-B0 in (S^2), C-C0 in (S^4). Writing B-B0=S^2*U and C-C0=S^4*V, polynomial membership in the generic native module is exactly

    V-U in (S),
    A_eff=(5*V-6*U)/(2*a^4),
    B_eff=2*(U-V)/(a^5*S).                              (3.5)

The parameter a is a unit in K0. Equation (3.5) is not whole-source landing over a=0. For example the generic retained profile a4=1/a is polynomial on this plane and has determinant (3.3), but its source correction -2*X^4/sigma has an uncancelled pole on t=0. Generic completion alone cannot construct a whole polynomial native map.

### Leading coefficients without a degree sweep

The finite r-degree reduction from UANATIVEORBIT1 remains valid. If m,n are positive degrees and g=gcd(m,n), their leading coefficients p,q satisfy

    n*p'/p=m*q'/q.

At every valuation v of K0(S), this gives n*v(p)=m*v(q). Therefore

    m/g divides every v(p), n/g divides every v(q).     (3.6)

After algebraic scalar extension this includes all finite roots and poles and infinity. A reduced non-dividing pair requires m/g,n/g>=2, so its nonconstant leading coefficients must be proper powers in their rational-function divisors. In particular a simple root or pole in p forces m|n and another reduction; likewise for q. Constant leading coefficients have zero divisor and provide no such exclusion. This is a necessary all-degree gate, not an existence or atomicity criterion.

### What the boundary does and does not say

There is only one critical divisor of the completion, namely D, with Jacobian multiplicity one. If it is not contracted, the associated transverse ramification index is2: at generic smooth source and image points, the separable tangential derivative is a unit and a local target transverse parameter pulls back as a unit times r^e, making the Jacobian order e-1. If D is contracted, this transverse argument does not apply. Trushin's contraction/branching distinction is relevant context; neither alternative is excluded here.

In particular HC4MYGL0 does NOT make F invertible: its Jacobian is nonconstant. The base completion already has collisions while mapping the line S=0 injectively.

## 4. UAPLANEATOM1: irreducible Jacobian plus an embedded line forces plane indecomposability

Here is a general characteristic-zero plane statement, independent of native formulas.

Let F:A2->A2 be polynomial, with Jacobian a nonzero constant times an irreducible nonconstant polynomial. Suppose the restriction of F to SOME affine line L is a CLOSED EMBEDDING. Then in every polynomial self-map factorization

    F=E composed G,

at least one factor is an automorphism over k. If the Jacobian factor is GEOMETRICALLY irreducible, the statement holds after EVERY characteristic-zero field extension. This is indecomposability in the polynomial PLANE category (absolute in the geometric case), not identity-stable atomicity or a general theorem about Keller maps in higher dimensions.

### Proof

Both factors are dominant. The chain rule factors the irreducible Jacobian as

    J_F=J_G*(J_E composed G).

If J_G is constant, G is Keller and is injective on L because F is. HC4MYGL0 makes G an automorphism.

Otherwise irreducibility forces J_E composed G constant, and dominance forces J_E constant. Thus E is Keller. Let s parametrize L. Since F|L is a closed embedding,

    k[s]=k[F1|L,F2|L] subset k[G1|L,G2|L] subset k[s].

Equality makes G|L a closed embedding as well. Over an algebraic closure, the Abhyankar--Moh theorem rectifies its image G(L) by a polynomial plane automorphism A. The map E composed A is then Keller and injective on the straight line A^(-1)(G(L)), because E is injective on G(L). HC4MYGL0 makes E composed A, hence E, an automorphism. Descent gives the statement over the original field. This proves the theorem over k without assuming JC(2). When the Jacobian factor is geometrically irreducible, repeat the argument after any field extension.

An arbitrary finite composition consequently has only one nonautomorphism stage. The nonconstant Jacobian prevents F itself from being an automorphism.

### Application and its precise full-map boundary

Every native polynomial-r completion in (3.3) satisfies the theorem. It is therefore indecomposable through polynomial plane self-maps in ALL source-coefficient and r degrees, including any non-dividing degree survivors. This does not imply that the whole native Keller map is absolutely stably atomic: its actual generic source leaf is the OPEN complement of D, not this plane, and intermediate coordinates may have poles on D or move the retained signal.

There is a useful whole-map factor exclusion. Suppose a native polynomial-r T has a whole polynomial factorization

    T=E composed G,
    G=(u,sigma,Gb,Gc), E=(u,a,Eb,Ec),

so both stages retain the displayed signal and all identities. If Gb,Gc also have no negative r powers on the generic leaf, G extends to a polynomial plane map Gbar by (3.2). Both whole factors are Keller, by the chain rule for T. Over K0 the target pair Ebar is Keller and

    F=Ebar composed Gbar, J_Gbar=constant*r.

The plane theorem forces Ebar to be an automorphism, since Gbar is not. Thus E is birational on the whole affine space and Keller, hence an automorphism by the birational-etale self-map criterion used in UAALLATOM1.

**Any nontrivial signal-retaining factorization must therefore use intermediate coordinates with a genuine pole on D**, or depart from the stated retained-coordinate architecture. This pole condition concerns the intermediate map, not the final polynomial-r outputs. We do not assume that a general plane Keller target cannot hide such a pole; that would exceed the imported line theorem. Identity stabilization mixing the retained coordinates is outside this exclusion.

## 5. What is now left

The open polynomial-r native problem is an exact plane-map problem with determinant 1-SY, the native module congruence, and separate whole-source landing at special signals. Its completion is already plane-indecomposable. Any new full-map decomposition respecting the signal must cross the boundary through a pole, while arbitrary intermediate signals and unrestricted source changes remain separate.

A genuinely nonlinear C output excludes every literal retained-profile first factor automatically by UAPROFFACT1. The scalar tangent would escape all such factors if it algebraized, but its finite polynomial integration remains unconstructed. The full canonical Spec(A), all thickened marked restrictions, effective gluing, absolute stable atomicity and improved costs remain open requirements. No further coefficient-degree or formal-order campaign is scheduled in a closed sector.

The exact controls verify formulas and sharp boundary examples; the line-injectivity, embedded-line rectification and general factor arguments require independent proof review. Publication priority of these compiler applications is unverified. The saturation statement has explicit earlier public precedent.
