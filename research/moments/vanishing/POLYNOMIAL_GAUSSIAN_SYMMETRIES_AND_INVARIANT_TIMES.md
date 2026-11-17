# Full polynomial Gaussian symmetries and invariant-time positive families

15 November 2026. Written internal deductions; independent replay, external review and formal verification remain pending. Exact finite controls are not universal proofs. No novelty or complete safe-phase classification is asserted.

This completes the moment-symmetry criterion left outside the scope of [the earlier Gaussian shear note](GAUSSIAN_SHEAR_NORMAL_FORM_AND_TRANSFERS.md). Quadratic-form preservation alone had not been used there as an unproved general invariance rule. Below it is proved sufficient for polynomial automorphisms over reduced coefficient rings, and the precise correction over arbitrary Q-algebras is derived. All invertibility assumptions are explicit: no Jacobian-conjecture inference is made.

## 1. Setup and the separating Gaussian pairing over coefficient rings

Let A be any commutative Q-algebra, R=A[X_1,...,X_n], and B an invertible symmetric n-by-n matrix over A. Define

    S(X)=X^T B^(-1)X/2,
    L_B(exp(s dot X))=exp(s^T B s/2).

The second equation defines each polynomial moment by finite coefficient extraction. Standard Gaussian expectation is B=I. It satisfies

    L_B(partial_i F)=L_B((B^(-1)X)_i F).                       (1)

It is the unique A-linear functional with L_B(1)=1 and (1), since multiplication by B converts (1) to the recursive relations for every X_i F.

Its bilinear polynomial pairing is separating even when A has nilpotents. Indeed the formal identity

    L_B(exp(t^T B^(-1)X-t^T B^(-1)t/2) F(X))=L_B(F(X+t))

expanded as sum_alpha H_alpha t^alpha/alpha!, defines polynomials H_alpha and yields L_B(H_alpha F)=L_B(partial^alpha F). For a nonzero F, choose a nonzero coefficient c_alpha of maximal total degree. Then partial^alpha F=alpha! c_alpha, and its pairing is nonzero because alpha! is a unit in A. Thus L_B(FG)=0 for all G implies F=0. No positivity, embedding of A in C, finite-dimensional Frobenius assertion or Noetherianity is needed.

## 2. Complete criterion for a supplied polynomial automorphism

Let Phi=(phi_1,...,phi_n) be an A-polynomial automorphism, with polynomial inverse. Write j=det Jac(Phi), a polynomial unit, and u=j(0), a unit in A. Every nonconstant coefficient of a polynomial unit is nilpotent: reduce modulo each prime of A, where units in a polynomial ring over a domain are constant, and intersect all primes. There are finitely many such coefficients; if their individual nilpotence indices are a_i, their ideal has zero power 1+sum_i(a_i-1), even if A is not Noetherian. Consequently

    ell=log(j/u)=sum_(k>=1) (-1)^(k+1)(j/u-1)^k/k             (2)

is a finite, canonical polynomial with ell(0)=0. Its derivatives satisfy d ell=j^(-1) d j.

**Criterion.** Phi preserves every polynomial Gaussian moment,

    L_B(F composed Phi)=L_B(F) for every F in R,               (3)

if and only if the polynomial

    C_Phi=S composed Phi-S-ell                               (4)

is constant in X. The constant may be nonzero and nilpotent; requiring it to be zero would incorrectly reject genuine symmetries.

Proof. Put D_i=sum_k v_(ki) partial_k with column v_i the i-th column of (Jac Phi)^(-1). Then D_i(phi_j)=delta_ij and the D_i form an R-basis of derivations. The polynomial cofactor identity

    div(j v_i)=0,   hence div(v_i)=-D_i ell                   (5)

is the Piola identity: expand each cofactor as an alternating sum of products of partial derivatives of all phi_j with j!=i. Differentiating and pairing terms with the two exchanged differentiation indices cancels every term by commutativity of mixed partials. This proof is valid over A; dividing by the unit j gives the second identity.

Polynomial Gaussian integration by parts gives, for every G,

    L_B(D_i G)=L_B((D_i S-div(v_i)) G)
              =L_B((D_i S+D_i ell)G).                        (6)

If (3) holds, transport (1) through Phi, using its inverse to let G range over every polynomial. This gives

    L_B(D_i G)=L_B(((B^(-1)Phi)_i)G).

The separating pairing in Section1 and D_i(S composed Phi)=(B^(-1)Phi)_i imply D_i C_Phi=0 for every i. Since the D_i are a derivation basis, all ordinary partials of C_Phi vanish. A Q-algebra polynomial with this property is constant.

Conversely, if (4) is constant, equation (6) gives the transported integration-by-parts relations. Therefore F->L_B(F composed Phi) is normalized and satisfies (1), so uniqueness from Section1 proves (3). This proves the equivalence in every dimension and degree over every Q-algebra.

An equivalent infinitesimal test for a polynomial derivation V=sum v_i partial_i is

    L_B(VF)=0 for every F  iff  div(v)=V(S).                  (7)

This follows directly from (6) without a change of coordinates and the separating pairing. If V is locally nilpotent, exp(V) is a polynomial Gaussian symmetry under (7); all positive powers of V have zero expectation. More generally exp(tV) preserves the functional for formal t whenever interpreted coefficientwise; a polynomial automorphism requires separate local nilpotence or nilpotent time. The implemented map certificate asks for and checks a polynomial inverse; it does not solve invertibility from a Jacobian determinant.

## 3. Reduced rings: exactly the quadratic-form automorphisms

If A is reduced, every polynomial unit j is constant, so ell=0. More generally this conclusion applies to any A-polynomial automorphism whose Jacobian determinant is constant. In this case (3) holds if and only if

    S composed Phi=S.                                       (8)

Indeed Section2 first says S composed Phi-S is constant. Taking its gradient at0 gives (Jac Phi)(0)^T B^(-1)Phi(0)=0. The matrices are invertible, so Phi(0)=0. Evaluation at0 makes the constant zero. The converse follows from Section2. This proves, in particular, that over every characteristic-zero field the full polynomial Gaussian moment symmetry group is exactly Aut(S), with no tameness assumption. These are polynomial automorphisms preserving the quadratic polynomial, rather than only linear orthogonal transformations.

For general nonreduced A, neither quadratic-form preservation alone nor constant volume alone suffices. The full logarithmic criterion (4) is the boundary. If (8) holds, (3) is equivalent to j being constant: ell must be constant and ell(0)=0, and finite log/exp on the nilpotent polynomial ideal are inverse.

## 4. Exact positive and negative nilpotent controls

Use a circular pair Z,W with L(ZW)=1 and S=ZW (plus any independent real spectators).

Over A=Q[e]/(e^2), the map

    Z -> Z-e Z^2,
    W -> W+e(ZW-1)                                          (9)

has inverse obtained by replacing e with -e. It has j=1-eZ, ell=-eZ, and S composed Phi-S=-eZ. Thus C_Phi=0 and every Gaussian moment is preserved, despite failure to preserve S. It is exp(eV) for V=-Z^2 partial_Z+(ZW-1)partial_W, with div V=V(S)=-Z. The field derivation V is not locally nilpotent; the nilpotent time e makes the exponential finite.

A second true symmetry shows why the constant in (4) is allowed. Set

    V=(1-Z^2-ZW)partial_Z+(-1+ZW+W^2)partial_W

and use Phi=exp(eV) modulo e^3. Then div V=V(S)=W-Z. Its coordinate values are id+eV(id)+e^2 V^2(id)/2. Here Phi(0)=(e,-e), j(0)=1-e^2, and C_Phi=-e^2, a NONZERO constant in Q[e]/(e^3). Criterion (4) accepts it correctly.

For a negative control add a real T and S=ZW+T^2/2. Over e^2=0 set

    Z -> Z,   T -> T+eZT,   W -> W-eT^2.                    (10)

This is an automorphism with inverse e->-e and it preserves S exactly. But j=1+eZ, ell=eZ, and C_Phi=-eZ is nonconstant. In fact L(W composed Phi)=-e!=0. Thus quadratic-form preservation without constant volume is insufficient on nonreduced rings. An ordinary translation has constant Jacobian but fails (8), giving the separate volume-only boundary.

## 5. Invariant-time null rotations, including Nagata type

In three split coordinates let q=2ZW+T^2, D=Z partial_T-T partial_W. Then Dq=DZ=0 and div D=D(S)=0 for S=q/2. For arbitrary f in A[Z,q], the derivation fD is locally nilpotent, since D is locally nilpotent and Df=0. Its exponential is

    Phi_f:
      Z -> Z,
      T -> T+Zf,
      W -> W-Tf-Zf^2/2.                                    (11)

Its inverse is Phi_(-f), q is fixed, and j=1. One can also verify the determinant directly: for a general polynomial f, it is 1+Z f_T-T f_W=1+Df, and the displayed f has Df=0. Thus (11) preserves all Gaussian moments over every Q-algebra.

Taking f=q gives, after x=2W,y=T,z=Z, the classical Nagata formula. The formula, the invariant-time construction and its inverse are established prior work: see Huarcaya--Palacios, [arXiv2401.12523v2](https://arxiv.org/html/2401.12523v2), Proposition2 and its proof, with the determinant in Proposition1. We rederive the required identities and add the Gaussian implications; no novelty or wildness theorem is imported. The source's wider wildness/classification results are not dependencies here.

There is a stronger cutoff than the coordinate degree bound. If f has Z-degree at most j when regarded as a polynomial in the independent variables Z,q, then

    P=W composed Phi_f

has the Gaussian cutoff m>(2j+1)deg Q, independent of the q-degree of f. To prove it, use exact symmetry to replace Q by Q composed Phi_(-f) while the phase becomes W. Give Z,W,T angular charges1,-1,0, so q has charge0. The inverse coordinates have maximum charges1,j+1,2j+1 respectively. Every degree-delta marker therefore has maximum charge at most(2j+1)delta. A nonzero expectation with W^m requires charge m in the marker. This proves the cutoff with no degree bound on q.

The coefficient is sharp for the entire monomial family f=Z^j q^k, j,k>=0. Put m=(2j+1)d. For Q=W^d the sole maximum-charge term in Q composed Phi_(-f) is (-1/2)^d Z^m q^(2kd). Hence

    L(W^d P^m)=(-1/2)^d m! 2^(2kd) (m+3/2)_(2kd) !=0.       (12)

Here (a)_b is the rising factorial. The radial identity used is L((ZW)^m q^s)=m!2^s(m+3/2)_s; it follows by extracting t^s in L((ZW)^m exp(tq))=m!(1-2t)^(-m-3/2), or directly by the standard Gaussian radial moment integral. It is a rational Wick identity, so holds over every Q-algebra. In particular f=q^k gives phases of degree4k+1 with the SAME sharp coefficient1 for all k. Arbitrarily large degree in this invariant does not force a growing cutoff.

These maps and phases extend the old f(Z)-only shear family. They do not classify all Aut(q) or all safe phases, and no assertion that each resulting PHASE lies outside every older phase orbit is needed.

## 6. Full orbit transfer and the differential boundary

For every symmetry from Section2, with inverse Psi,

    L_B((P composed Phi)^m)=L_B(P^m),
    L_B(Q (P composed Phi)^m)=L_B((Q composed Psi)P^m).        (13)

Thus all pure sequences, qualitative mixed safety/failure and actual fixed witness sequences transfer exactly. A common bound C deg Q transfers with coefficient C deg(Psi); use maximum coordinate degree and an upper bound if coefficients cancel. Exact filtered marker spaces transport bijectively; the same ordinary degree cutoff need not be preserved. In shared parameter coefficient rings satisfying the symmetry and inverse identities, pure sequences and their full ordinary coefficient ideals agree literally, including nilpotents. No equality of same-degree mixed ideals is asserted.

Over C, enlarge the ternary cubic orbit to ALL polynomial Aut(q) maps. Every nonzero pure phase in this orbit is safe; first94 pure tests suffice and its least universal pure prefix still equals the cubic least value in[20,94], since all pure sequences agree and the identity map is allowed. The ordinary all-ring gate is not94: for arbitrary coefficient rings use full ordinary purity of the base cubic and its retained sufficient ring slope, then multiply the slope by the inverse degree. Likewise the harmonic algebras in [the previous proof](HARMONIC_POLYNOMIAL_ALGEBRAS_AND_GAUSSIAN_CUTOFFS.md) generate further safe orbit classes with transported cutoffs. The retained quartic failure and all its Gaussian-symmetry transforms remain failures, so these positive orbit classes do not exhaust all pure-null phases.

Gaussian invariance does not imply preservation of the ordinary Laplacian. In the three split coordinates Delta=2 partial_Z partial_W+partial_T^2, the Nagata choice f=q gives

    Delta(W composed Phi_q)=-10T-14Zq !=0,

although Delta W=0 and every positive pure Gaussian moment of W composed Phi_q vanishes. Even the first-power differential hypothesis can fail. Pullback transports the differential operator to one with variable coefficients; it does not prove a constant-coefficient GVC statement. The direct harmonic-core GVC theorem remains separately valid.

## 7. Implementation and audit boundary

The [checker](verify_polynomial_gaussian_symmetries.py) implements exact symmetry certificates for supplied maps AND supplied inverses over Q(i) and Q(i)[e]/(e^L), checking inverse identities, the Jacobian unit, its finite normalized logarithm and constancy of (4). It supports fixed invertible rational covariance matrices, including split pairs. It does not infer invertibility from a constant Jacobian, solve general inverse recognition, classify every quadratic-form automorphism, or decide whether an arbitrary phase belongs to a safe orbit. The nilpotent corrections and the nonzero constant control are essential.

The [source audit](polynomial_gaussian_symmetry_source_audit.json) binds this proof, source snapshots, [receipt](polynomial_gaussian_symmetry_verification.json), checker and retained frontiers. The theorem is proved by Gaussian IBP, a separating pairing, the cofactor identity and uniqueness; exact finite maps/moments only check transcription. Cubic94, least intervals, formal-symbol eventuality-to-linearity and the other protected endpoints remain unchanged. Independent review and a general safe normal-form classification remain open.
