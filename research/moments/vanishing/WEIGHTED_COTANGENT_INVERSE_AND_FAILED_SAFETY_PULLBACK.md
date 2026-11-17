# Weighted cotangent inverses and failure of forward Gaussian safety

15 November 2026. Internal written proofs. Independent replay, external review and formal verification remain pending. Finite exact controls are not the proof of the infinite tails. No novelty or improved minimum-dimension claim is made.

This resolves OP-GMC-ENDOMORPHISM-SAFETY from the [preceding note](GAUSSIAN_KELLER_ENDOMORPHISMS_AND_ONE_WAY_TRANSFER.md) negatively. The counterexample uses the already known weighted six-dimensional cotangent map, credits its existing map/moment-map sources, and proves an additional mixed-moment identity through its formal inverse. The general result covers arbitrary degrees and weighted Keller maps, and also gives a direct generic Gaussian test for any normalized Keller inverse.

## 1. Weighted cotangent Gaussian evaluation is formal inverse substitution

Let K have characteristic zero. Let F in K[x_1,...,x_n]^n have F(0)=0, JF(0)=I, and det JF=1. Write A=JF. Suppose, for a diagonal matrix W=diag(w_i) with each w_i a nonzero element of K,

    A(x) W x=W F(x).                                        (1)

Thus F respects the indicated nonzero Euler weights. Let D=W^(-1), let g be the unique formal inverse of F at0, and give (x,p) the split Wick pairings

    L(x_i p_j)=D_i delta_ij,  L(x_i x_j)=L(p_i p_j)=0.

Set Y=A(x)^(-T)p. For every polynomial amplitude Q(x) and constant direction a,

    L(Q(x) exp(s a^T Y))=Q(g(s D a)).                        (2)

This is an identity in K[[s]], also valid over an exact extension field containing a. Every requested coefficient is finite. The denominator is exactly1; it has not been discarded.

Proof. Put h_a(x)=D A(x)^(-1)a. Rescale the momentum variables by W so their pairings with x are1. The retained [formal Gaussian--Lagrange lemma GL1](../../foundations/extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md), including its nonzero-constant-term completion, gives

    L(Q exp(s a^T Y))=Q(b(s))/det(I-s Jh_a(b(s))),
    b=s D A(b)^(-1)a.                                      (3)

This classical Lagrange--Good mechanism is credited to its sources; see [Abdesselam's account](https://arxiv.org/html/math/0208174v1) and the retained local GL1 proof for the precise formal residue convention. It is also the architecture used in Long's [retained source](sources/long_gaussian_de4082c.tex). No new Gaussian inversion formula is claimed as a classical mechanism.

Multiplying the fixed-point equation by A(b)W gives

    A(b)Wb=s a, hence W F(b)=s a,

so b(s)=g(sDa). To compute its Jacobian factor, differentiate (1) in direction x_k:

    (partial_k A)Wx=W A e_k-A W e_k.

Since partial_k A^(-1)=-A^(-1)(partial_k A)A^(-1), at the fixed point the k-th column of I-s Jh_a is

    e_k+D A^(-1)(partial_k A)Wb=D A^(-1)W A e_k.

Therefore

    I-s Jh_a(b)=D A(b)^(-1)W A(b),
    det(I-s Jh_a(b))=1.                                    (4)

Equations (3)--(4) prove (2). In particular Q=1 gives L(exp(s a^TY))=1: every such scalar phase is pure-null. The formal branch may move away from0 in s; the completed-ring argument in GL1 is essential and avoids an invalid x-adic change of variables.

## 2. A generic dual phase detects the polynomial inverse in all degrees

Take a_1,...,a_n algebraically independent over K and regard a as exact parameters. For

    P_a=a^T A(x)^(-T)p                                      (5)

over K(a), the following are equivalent:

1. F has a polynomial inverse over K.
2. P_a is safe against every polynomial multiplier.
3. Each of the n fixed position markers x_i has an eventually zero mixed tail.

Indeed (2) gives

    sum_m L(x_i P_a^m) s^m/m!=g_i(sDa).

A univariate specialization at algebraically independent a is polynomial in s exactly when g_i is a polynomial: each nonzero homogeneous component of g_i evaluates to a nonzero polynomial in a. Thus3 implies1. A polynomial formal inverse is a two-sided polynomial inverse because both formal composition identities then are polynomial identities. If1 holds, the cotangent map is a polynomial Gaussian automorphism and transports the safe isotropic target linear phase a^Tp, proving2 by the earlier exact automorphism theorem. Finally2 implies3 directly.

If F is noninvertible, some FIXED coordinate marker x_i has an infinite nonzero tail for the generic phase (5). This is an inverse criterion, not a claim that every special direction detects noninvertibility. Special lines can have polynomial inverse restrictions. If K is a finitely generated coefficient field inside C, one may choose actual complex a algebraically independent over K; the same conclusion then gives a complex counterexample. No general inverse recognition algorithm or search for the best direction is executed here.

## 2a. Complete safe-direction criterion, including all source markers

For a FIXED direction a, the phase P_a in (5) is safe if and only if the individual inverse restriction b(s)=g(sDa) is polynomial. This strengthens the generic criterion without asserting that every special direction sees a failure.

Necessity follows from the n coordinate markers in (2). For sufficiency introduce a momentum source rho and rescale p to Wp. The same completed formal-residue proof gives

    L(Q(x) exp(sP_a+rho^T Wp))=Q(b(s,rho))/j_s(b(s,rho)),
    b=rho+s h_a(b), j_s(x)=det(I-sJh_a(x)).                  (13)

At rho=0, b=b(s) and j_s(b)=1 by (4). Every derivative of b in rho is polynomial in s when b(s) is polynomial: first implicit derivatives use (I-sJh_a(b))^(-1), whose determinant is1 and whose adjugate is polynomial; higher derivatives follow inductively by differentiating the polynomial implicit equation. Derivatives of the reciprocal determinant have the same property. Differentiating (13) finitely many times in rho produces every polynomial momentum marker. Thus every source polynomial has a polynomial exponential moment, which is exactly eventual mixed vanishing.

There is a uniform finite direction certificate over characteristic-zero fields. If D0>=1 bounds deg F, put B0=D0^(n-1). Every polynomial inverse restriction b(s) has degree at most B0. For a!=0 choose n-1 independent target linear forms annihilating Da. Their pullbacks by F cut out a smooth pure one-dimensional preimage of the target line, because JF is invertible. A generic affine hyperplane cuts this curve in finitely many isolated points. If b(s) has degree e, a generic hyperplane composed with b is a separable degree-e polynomial; its e roots give e distinct points because F(b(s))=sDa. The affine Bezout isolated-point bound applied to the n equations of degrees at most D0,...,D0,1 gives e<=D0^(n-1). For a=0 the restriction is zero. Only this classical upper bound is imported; see Mondal, *How Many Zeroes?*, TheoremI.2, printed p.2, [author's preview](sources/mondal_how_many_zeroes_2021_preview.pdf). No equality case or positive-dimensional multiplicity assertion is used.

Compute the formal branch through degree B0 and call the resulting polynomial b_B0(s). Then

    F(b_B0(s))=sDa as a polynomial identity                  (14)

is equivalent to safety of the phase for that field direction. If (14) holds, formal uniqueness identifies b_B0 with b; conversely the degree bound makes (14) necessary. Coefficients of b_B0 are polynomial functions of a, so (14) defines the exact Zariski-closed FIELD safe-direction locus, with a theoretical finite algorithm and an implemented guarded small-input constructor. It is an inverse-section certificate, not a claimed short pure-moment prefix.

A common sufficient linear cutoff also follows directly. Let H0 bound deg_x h_a uniformly, for example H0=(n-1)(D0-1). Define

    C_p=(2n-1)(1+B0 max(H0-1,0)), C=max(B0,C_p).

For a safe direction, a monomial marker of position degree alpha and momentum degree beta has exponential-moment degree at most B0 alpha+C_p beta. To check this bound, write (13) as N/j_s^r and apply each rho derivative through the vector field (adj(I-sJh_a)/j_s) partial_x. One derivative changes the numerator to a sum of j_s*adj*partial N and -r*N*adj*partial j_s, and increases the denominator power by2. Numerator s-degree increases by at most2n-1 and x-degree by at most(2n-1)max(H0-1,0). After beta derivatives evaluate x=b(s), of degree<=B0, and j_s(b)=1. This gives the displayed bound. Thus m>C deg Q suffices for every safe field direction, without general Gaussian field stratification.

The field-exact degree bound is not an all-ring converse. Nilpotent direction scaling can truncate a genuinely infinite inverse series at an arbitrarily large nilpotent index. Identity (14) remains an ordinary SUFFICIENT certificate whenever it holds over a Q-algebra; it is not claimed to represent every safe nonreduced point. No large degree49 branch calculation for the six-variable example is run: its nonpolynomial branch is already explicit in (8).

## 3. An explicit direction in the known six-dimensional map

Use the normalized foundational F1 map f from the preceding note:

    u=1+xy,
    f1=x-3x^2y/2-x^3z/2,
    f2=y+3xu^2z+3xy^2(4+3xy),
    f3=u^3z+y^2u(4+3xy),
    W=diag(1,-1,-2), D=diag(1,-1,-1/2).                      (6)

Its inverse Jacobian A^(-1) is polynomial because det A=1. In the six split variables (x,y,z,p1,p2,p3), define the EXPLICIT polynomial phase and linear marker

    P=(1,-1,0) A(x)^(-T)p,
    Q=x.                                                   (7)

The phase has degree12. Formula (7), together with the displayed three small polynomials f and det Jf=1, specifies a polynomial by its adjugate; no rational inverse of f is used. The [receipt](weighted_cotangent_inverse_verification.json) includes its three coefficient polynomials in p and the exact fixed linear conversion to standard Gaussian coordinates.

For a=(1,-1,0), the target line in (2) is sDa=(s,s,0). Its formal inverse is explicit. Let

    U(s)=(1-2s^2)^(-1/2),
    x(s)=(U(s)-1)/s, y(s)=s,
    z(s)=-s^2(1+3U(s))/U(s)^2.                              (8)

All three are regular formal series with zero constant term. Substituting1+x(s)y(s)=U(s) makes f3=0 and f2=s. On f3=0 the first component simplifies to

    f1=x(2+xy)/(2(1+xy)^2).

The identity (1-2s^2)U(s)^2=1 then gives f1=s. Thus (8) is the unique formal inverse branch at0 on that line.

Apply (2) with Q=x:

    L(exp(sP))=1,
    L(x exp(sP))=((1-2s^2)^(-1/2)-1)/s.                     (9)

For every k>=1, coefficient extraction yields

    L(x P^(2k-1))=(2k-1)! binomial(2k,k)/2^k !=0,
    L(x P^(2k))=0.                                        (10)

This is a FIXED linear witness with infinitely many nonzero mixed moments. The universal tail follows from the written formal inverse and binomial series identities, not from checking finitely many moments.

The target phase P0=p1-p2 is an isotropic linear form, hence pure-null and safe with cutoff m>deg Q_target. The known cotangent map K_f=(f,A^(-T)p) preserves the full Gaussian state by GMCNONINV620261115, yet P=P0 composed K_f fails safety by (10). Therefore forward safety under moment-preserving Keller endomorphisms is FALSE, already for a safe linear target phase and a linear source witness. This resolves OP-GMC-ENDOMORPHISM-SAFETY negatively.

For ordinary independent real variance-one inputs use the fixed complex linear change

    x_i=U_i+iV_i, p_i=(U_i-iV_i)/(2w_i).

The source marker becomes U1+iV1, and the safe target phase becomes

    P0=(U1-iV1+U2-iV2)/2.

Conjugating the same map by this fixed change gives the counterexample to forward safety in ordinary Gaussian coordinates. The all-moment-preserving map is the previously known cotangent collision; its foundational map/collision are Alpoge's, the cotangent lift is credited to Mayner and SW1, and its moment-map identity to Kirtchakov's pinned TheoremC(i), as documented in [the prior source audit](gaussian_keller_endomorphism_source_audit.json). We claim no new map, independent formal assurance, novelty or improved smallest GMC dimension/degree. The existing ternary quartic failure remains the protected smaller witness.

## 4. Direct generic Gaussian phase for an arbitrary normalized Keller map

The weighted criterion is not restricted to Keller maps that initially have an Euler grading. For any normalized Keller F in n variables, use the preceding scaling family F_t=t^(-1)F(tx), A_t=J_xF_t. Its position map

    G(x,t)=(F_t(x),t)

is a normalized Keller map with nonzero Euler weights (1,...,1,-1). Its signed cotangent lift is exactly the previously constructed H_F in dimension2n+2. Define

    R_a=z+(a+partial_t F_t)^T A_t^(-T)y.                    (11)

This is a polynomial phase, with coefficients in K(a) for generic parameters a. For every polynomial position amplitude Q(x,t), formula (2) on G gives

    L(Q(x,t) exp(s R_a))=Q(s^(-1) g(a s^2),s),              (12)

where g is the formal inverse of F. The source Wick pairings are L(x_i y_j)=delta_ij and L(tz)=1. In the signed cotangent notation the final momentum is eta=-z; the direction is (a,-1), and its covariance inverse weight is also-1. These signs give t=s and the plus sign in (11).

In particular every positive pure moment of R_a is zero. Over K(a), full mixed safety of R_a, eventual vanishing for just the n fixed source markers x_i, and polynomial invertibility of F are equivalent. If F is noninvertible, a fixed x_i gives a nonpolynomial series s^(-1)g_i(a s^2). If F is invertible, H_F is a polynomial Gaussian automorphism and transports the safe isotropic target linear phase a dot y+z.

For any FIXED a, Section2a applied to G shows more: full safety of R_a is equivalent to polynomiality of s^(-1)g(a s^2), hence to polynomiality of the inverse restriction g(at). The finite inverse-section certificate F(g_B(t))=ta uses the bound B=(deg F)^(n-1) from the same affine-line argument, which never required Euler weights. Therefore this whole phase family for an arbitrary normalized Keller map has an exact closed FIELD safe-parameter locus. The ordinary equations remain only a sufficient certificate over nonreduced algebras.

Thus every noninvertible normalized Keller map in n variables yields a generic pure-null Gaussian failure in dimension2n+2 with a fixed linear source witness. This is a direct all-degree construction, without cubic homogenization or dimension-optimality assertion. Its specified-input phase formula and covariance/sign checks are implemented. No large expansion of the eight-variable foundational suspension is needed: the explicit weighted six-variable direction (7) already resolves the transfer question.

The generic coefficient field is a rational function field over the finite coefficient field of F. In the exactly presented effective-field setting, the retained theoretical fixed-point mixed-safety algorithm could decide the resulting phase and recover a witness. That general solver is not implemented or executed here; the formulas, inverse-equivalence proof and explicit six-variable witness do not depend on it.

## 5. Verification and remaining work

The guarded inverse-section constructor implements (14) for supplied small inputs and refuses larger configured orders without a classification. The [checker](verify_weighted_cotangent_inverse.py) verifies the Euler derivative/Jacobian-cancellation identities, the algebraic inverse branch relation in a quadratic function field, exact Wick moments using finite componentwise coefficient bounds, covariance/sign/transpose mutations, and a small universal-constructor example. The wrong transpose already changes the second pure moment to-3; omitting the covariance sign changes the third marked value from9 to-9. These are transcription controls only. The [source audit](weighted_cotangent_inverse_source_audit.json) binds proof, receipt, source, changed records and protected evidence.

The earlier one-way identities remain valid. Their previously missing converse now has an explicit counterexample; this is an actual negative closure, not a lack of proof. Full automorphism transfer, pure-prefix94 and all ordinary-ring distinctions remain unchanged. For a fixed noninvertible map, classify the phases whose safety does survive and identify their source-marker defects; that broader classification remains separate. Independent specialist review, full moment-map necessity without an inverse, general formal-symbol eventuality-to-linearity, actual higher-rank modules and general solver implementations remain open.
