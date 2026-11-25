# The even D=0 quintic face: balanced Prym multiplication and ordinary primes

Let K have characteristic zero, A_K=K[x,y,t]/(t^2+2xy-1), and

    Lambda(f)=CT_x (1/2) integral[-1,1]
                 f(x,(1-t^2)/(2x),t) dt.

The integral is the rational polynomial rule. Eventual purity means Lambda(f^m)=0 for all sufficiently large positive integers m.

Consider the even stripped space

    f=P(y,t)+x*t^2*(q0+B*t^2),
    P=sum_(j+r<=5,r even) p_jr*y^j*t^r,
    A=p_14, C=p_32, D=p_50.                              (1)

All lower coefficients are arbitrary. The main new result concerns D=0 with BC!=0. It uses one number-field specialization and does not require the parameter-algebraicity theorem, a universal moment cutoff, a generic endomorphism assumption, or an ordinary-prime conjecture.

**Theorem 1.** For EVERY characteristic-zero K and EVERY phase(1) with D=0 and BC!=0, infinitely many positive moments are nonzero. In particular no such phase is eventually pure.

The subsequent boundary consequence uses the separately established exact Newton necessary conditions; it is not silently included as an assumption of Theorem1.

## 1. A balanced rank-two CM action gives quaternion multiplication

Here is the Hodge-theoretic fact needed for all special members.

**Lemma 2.** Let A be a complex polarized abelian fourfold. Suppose E is a quartic CM field with maximal real subfield F, and E acts faithfully through End^0(A). Suppose Rosati restricts to complex conjugation on E and the E-signature is(1,1,1,1). Then End^0(A) contains a quaternion F-algebra D containing E. This assertion requires neither simplicity of A nor equality End^0(A)=E. The algebra D is split at both real places of F.

Proof. Write V=H_1(A,Q), so V is two-dimensional over E. The E-compatible polarization gives a nondegenerate E/F-Hermitian form h on V, preserved by the Hodge circle. One construction starts from the alternating Q-valued polarization psi. The trace pairing on E gives the unique E-valued sesquilinear form phi with

    Tr_(E/Q)(e*phi(v,w))=psi(e*v,w) for every e in E.

The Rosati relation psi(e*v,w)=psi(v,bar(e)*w) shows that phi is linear in v, conjugate-linear in w and skew-Hermitian. Dividing by a nonzero purely imaginary element of E turns it into h. The polarized Hodge condition shows that h has opposite signs on the two Hodge-circle eigenspaces at each real embedding of F. Each such eigenspace has dimension one by the assumed signature. Hence h has signature(1,1) at both real places.

The Hodge circle acts E-linearly and preserves h. At every embedding of E, it has one eigenvalue z and one eigenvalue z^-1 for z in the unit circle. Its E-determinant is therefore1. Thus the Hodge circle lies in the real points of Res_(F/Q) SU(h).

Choose an E-basis diagonalizing h as diag(a,b), with a,b in F nonzero. Define the Q-linear, E-antilinear map

    j(v1,v2)=(b*bar(v2),-a*bar(v1)).                      (2)

It satisfies

    j*e=bar(e)*j,       j^2=-a*b.

It commutes with SU(h). For completeness, a matrix of determinant1 preserving diag(a,b) has the form

    M=[[r,s],[-(a/b)*bar(s),bar(r)]].

Direct multiplication gives j*M=M*j, regarding j as its displayed matrix followed by coefficient conjugation. Therefore j commutes with the Hodge circle and is a rational Hodge endomorphism of V.

A rational Hodge endomorphism of H_1(A,Q) is an element of End^0(A): clear denominators on H_1(A,Z), use commutation with the complex structure to obtain an analytic homomorphism of the complex torus, and use the algebraicity/full-faithfulness theorem for complex abelian varieties. This is the precise standard input stated in Milne, Abelian Varieties v2.0, Theorem2.9; no higher-codimension Hodge conjecture is involved.

The algebra E+E*j has dimension four over F and the quaternion relations in(2). Its action is faithful because the corresponding quaternion algebra is simple and the action is unital. Since a*b is negative at both real places, j^2=-ab is positive there. Thus this quaternion algebra is split over both real completions. QED.

This argument proves an EMBEDDED quaternion algebra. Additional endomorphisms, isogeny decompositions, special CM members and the case D=M_2(F) are all allowed.

## 2. A rational imaginary quadratic field in this quaternion algebra

**Lemma 3.** Let F be a real quadratic field and D any quaternion F-algebra split at its real places. Then there is an imaginary quadratic field K0/Q with an F-algebra embedding FK0->D. If D acts on the fourfold in Lemma2, the induced K0-action has signature(2,2).

Proof of existence. Only finitely many finite places v of F ramify in D. At each underlying rational prime p choose a rational p-adic squareclass d_p which remains nonsquare in EVERY F_v at a ramified v over p.

If p splits in F, its completions are Q_p, so any nonsquare class works at both. If p does not split, F_v is quadratic over Q_p. Choose a quadratic extension Q_p(sqrt(d_p)) different from F_v; then d_p cannot become a square in F_v. Such a choice exists also at p=2: a p-adic field Q_p has distinct ramified and unramified quadratic extensions, and a fixed quadratic F_v cannot contain both. These cases exhaust a quadratic field F; there is no place of larger local degree to handle.

Each chosen squareclass contains an open p-adic congruence class of integers after choosing valuation0 or1 representatives. The Chinese remainder theorem gives an integer d in all these classes; subtracting a sufficiently large multiple of their joint modulus makes d<0. Put K0=Q(sqrt(d)). The field K0 is imaginary and disjoint from the real F, so L=FK0 is a quadratic FIELD over F. By construction no ramified finite place of D splits in L. No real place is ramified here; in any case L is complex at both real places.

The exact classical input is Voight, Quaternion Algebras, Proposition14.6.7: a quadratic extension L/F embeds in a quaternion F-algebra if and only if no ramified place of that algebra splits in L. It supplies the embedding L->D. This is a field embedding, not a statement about integral orders or a chosen maximal order.

Proof of signature. The F-action decomposes Lie(A) into two complex spaces of dimension2, because each real embedding of F has two E-embeddings, each of Lie multiplicity1. On either space D tensor_(F,sigma) C is M_2(C), and its unital representation in dimension2 is the standard module. An element sqrt(d) of K0 has reduced characteristic polynomial X^2-d over F; after each complex embedding it has the two distinct eigenvalues sqrt(d),-sqrt(d), once each. The two F-places together give multiplicities2,2 for the two embeddings of K0. QED.

Only the actual real-quadratic center is used. No claim about an arbitrary higher-degree center or prescribed imaginary quadratic field is made.

## 3. The exact ordinary-prime theorem and geometric descent

The primary ordinary input is Fite, Ordinary primes for some varieties with extra endomorphisms, arXiv2005.10185v3, THEOREM4: every abelian fourfold over a number field whose GEOMETRIC End^0 contains an imaginary quadratic field acting with signature(2,2) has a nonzero-density set of ordinary primes. It assumes neither geometric simplicity nor equality of the displayed field with End^0. Theorem2 of that source is a different, half-ordinary statement and is not substituted here.

Consequently every number-field fourfold satisfying Lemma2 has positive-density ordinary primes. To justify the word geometric, endomorphisms obtained by the complex-torus argument are defined over Qbar when A is. Indeed an algebraic complex homomorphism descends to some finite-type Qbar coefficient variety. After shrinking, it is a family of homomorphisms A×T->A. Choose a Qbar point of the integral base T. The difference from that specialized homomorphism vanishes on A times that point and on the zero section; the rigidity theorem makes the family constant. This is the ordinary homomorphism base-extension invariance for abelian varieties, and follows from the same complete-source rigidity input (Milne Theorem1.1). Each of finitely many resulting endomorphisms is then defined over a finite number-field extension.

Thus the quaternion and imaginary quadratic subfield of Lemmas2--3 lie in the GEOMETRIC endomorphism algebra required by Fite. If necessary enlarge the number field first to define the selected factors, endomorphisms and deck roots of unity. Fite's theorem still applies over that enlarged field. No generic/full-endomorphism equality is used.

## 4. One specialization preserves all tail identities

Suppose, for a contradiction to Theorem1, that a phase(1) with D=0 and BC!=0 is eventually pure over some K. Its finitely many coefficients generate a finite-type Q-domain R inside K. Invert BC in R and, when A^2+2BC is nonzero, invert that discriminant as well. The result is still a nonzero finite-type Q-algebra. A maximal ideal has a number-field residue field. The resulting single specialization preserves the chosen nonzero coefficients and discriminant branch.

Every sufficiently late moment is an exact rational polynomial identity in these finitely many coefficients, zero in R. The SAME specialization preserves all of them simultaneously. There is no prime-dependent or moment-dependent choice of specialization. Therefore it is enough to exclude eventual purity with all coefficients in a number field.

Over a finite extension choose rho^4=C/B and multiply the actual angular transform x->rho*x,y->rho^-1*y by mu=1/(rho*B). This fixes B=C=1 and replaces A by

    lambda=A/(B*rho^2),  lambda^2=A^2/(BC).               (3)

The support, D=0 condition and moment vanishing are preserved. No algebraicity assertion about the original unspecialized lambda is needed. After this normalization, enlarge the number field F to contain zeta_5 and all fixed data used below.

## 5. The exact leading residue and its repeated-root branch

On the Laurent chart the ordinary weight W=h+i of x^i*t^h is at most5. With D=0 its equality face is

    x*t^4 * [1-(lambda/2)z-z^2/8],  z=t^2/x^2.           (4)

For y^j*t^r the highest Laurent term is (-1/2)^j*x^-j*t^(2j+r), with W=j+r; its lower terms have smaller W. Thus lower coefficients of P and q0 do not contribute to this face.

At p=10k+1, m=2k, angular balance selects z-degree k. The maximum angular-constant t height is10k=p-1. For the number-field coefficients, omit the finite primes where a coefficient is nonintegral. Multiplying the exact rational moment by p is integral at all remaining such primes. Only integration of t^(p-1) has a denominator divisible by p, so reduction gives

    p*Lambda(f^(2k)) modmathfrak p
       =H_p(lambda),
    H_p(lambda)=[z^k](1-(lambda/2)z-z^2/8)^(2k).          (5)

This is an exact identity; all lower heights have prime-to-p denominators. Any possible odd height integrates to zero. For all sufficiently large primes in this progression, eventual purity would force H_p(lambda)=0.

If lambda^2=-2, the polynomial in(4) is (1-lambda*z/4)^2. Hence

    H_p(lambda)=binom(4k,k)*(-lambda/4)^k !=0 modp,        (6)

since 4k<p and lambda is a unit outside finitely many primes. There are infinitely many primes p=1 mod10. This already contradicts eventual purity in the repeated-root case. No ordinary-prime theorem is used for this branch.

It remains to treat lambda^2+2!=0, including lambda=0. Put

    Rtilde(z)=z^2-(lambda/2)z-1/8.

Its two roots are distinct and nonzero.

## 6. The genus-eight cover and the balanced Prym

Take the smooth projective normalization of

    C_lambda: w^5=u^2*Rtilde(u^2).                       (7)

The six branch points have raw exponents(2,4,1,1,1,1): exponent2 at u=0, exponent4 at infinity in the homogeneous degree10 completion, and four simple nonzero finite roots. Each exponent is prime to5, so the cover is connected and Riemann--Hurwitz gives genus8.

For deck characters i=1,2,3,4, normalization gives the H^1(O) dimensions

    2i-1-floor(2i/5)-floor(4i/5)=(1,2,2,3).             (8)

Here the eigensheaf is O(-2i+sum floor(i*a/5)); the four simple exponents add no floors for i<5. This is the elementary integral-generator formula: near a branch of exponent a, divide w^i by the local parameter to floor(i*a/5).

The involution iota:u->-u,w->w commutes with the order-five deck action. Its quotient is the normalization of

    Cbar_lambda: w^5=z*Rtilde(z),  z=u^2,

with raw exponents(1,1,1,2), genus4 and character dimensions

    i-1-floor(2i/5)=(0,1,1,2).                           (9)

Let P_lambda be the image of1-iota on Jac(C_lambda), up to the usual isogeny splitting. Its dimension is4 and subtracting(9) from(8) gives E=Q(zeta_5) signature

    (1,1,1,1).                                         (10)

The restriction of the Jacobian polarization to this abelian subvariety is a polarization. The deck automorphisms preserve it, so the Rosati involution on E is zeta_5->zeta_5^-1. The E-action is faithful since its nontrivial order-five deck characters occur; Phi_5 is irreducible over Q. All hypotheses of Lemma2 hold at EVERY lambda with distinct nonzero roots. In particular P_lambda has the quaternion algebra and the signature(2,2) rational imaginary quadratic subfield of Lemmas2--3, even when it is nonsimple or has extra endomorphisms. Fite's Theorem4 supplies positive-density ordinary primes for P_lambda over F (or after a fixed further enlargement).

## 7. The selected Frobenius line lies in the Prym

The character-one eigensheaf for C_lambda is exactly O(-2); none of the raw branch exponents in(7) reaches5. A Cech generator of its one-dimensional H^1 is w/u. The involution iota acts on it as minus1. Thus this ENTIRE selected line belongs to P_lambda, not to the quotient factor.

At a sufficiently large good absolute-degree-one prime mathfrak p of F, its rational characteristic satisfies p=1 mod5 because F contains zeta_5. Since p is odd, write p=10k+1. Frobenius on w/u is computed by w^p=w*(u^2*Rtilde(u^2))^((p-1)/5). Its selected coefficient is

    [u^(p-1)] (u^2*Rtilde(u^2))^((p-1)/5)
      =[z^(3k)] Rtilde(z)^(2k)
      =[z^k](1-(lambda/2)z-z^2/8)^(2k)
      =H_p(lambda).                                    (11)

The middle equality is polynomial reversal, not an identification of different coefficients. In Cech terms this is the coefficient of u^-1 after applying pth power to the generator; the other powers are coboundaries. Good-prime smooth normalization is secured by excluding the fixed finite branch/discriminant and isogeny primes.

If H_p(lambda)=0, Frobenius kills this nonzero invariant line in H^1(O_(P_lambda)). Therefore P_lambda is nonordinary. Indeed an ordinary abelian variety has bijective semilinear Frobenius on H^1(O); the p-rank is its stable rank. This standard criterion also follows from the Artin--Schreier sequence and is preserved under isogeny away from the fixed splitting degree.

But(5) makes this zero-line condition hold at EVERY sufficiently large absolute-degree-one prime of F under the hypothetical eventual purity. Prime ideals of larger absolute residue degree have density zero. Thus ordinary primes for P_lambda would have density zero, contradicting the positive density given by Fite's theorem. This proves Theorem1 in the distinct-root branch. Together with(6) it proves all of Theorem1. Equivalently, every phase in the stated chart has nonzero moments at infinitely many orders. The proof is independent of PFP parameter algebraicity and global purity-recognition effectiveness.

## 8. The complete mixed-even D=0 boundary consequence

Use the existing canonical Newton result SNF-EVEN-TWO-NEGATIVE-TOPCOEFFICIENTS-NECESSARY, in STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md. It states exactly: an even eventually pure phase in the stripped24-space, of actual degree5 and with BOTH positive and negative angular terms, has B=[x*t^4]f!=0 and at least two nonzero coefficients among

    A=[y*t^4]P, C=[y^3*t^2]P, D=[y^5]P.

Its proof also excludes the v2 even case and the v4 degreeP<=4 case; a single negative degree-five leader is excluded by an exposed two-term face. This is a declared old premise for the following consequence, not an unproved classification of generic four-root covers.

**Corollary 4.** An even stripped eventually pure phase of actual degree5 with both angular signs MUST have D=[y^5]P!=0. In particular the entire genuine MIXED even D=0 quintic sector is excluded. Positive-only and negative-only pure safe possibilities are not excluded.

Proof. If D=0, the old necessary condition forces both A and C nonzero as well as B nonzero. This is a subcase of Theorem1, a contradiction. QED.

For the precisely normalized PFP even D=0,BC!=0 pure-parameter chart, the exact geometric PURE parameter set of the effective elimination theorem is consequently EMPTY. Its squarefree output polynomial is1, by this analytical proof. This does not mean that the large elimination algorithm has been executed, nor does it assert an exact parameter-set conclusion for an unbounded union of unrelated eventual thresholds. Here each phase in the particular chart has been excluded directly.

No conclusion about the full even BD!=0 leading face, its Prym sixfold, arbitrary stripped quintics outside this sector, or global primitive degree-five existence is asserted.

## 9. Primary inputs and honest scope

1. Fite, arXiv2005.10185v3: actual Theorem4 and its complete supplied Section8 proof were read; the earlier Proposition14 and representation-theoretic inputs remain classical imported results unless separately audited. The theorem covers ALL fourfolds with the stated geometric imaginary quadratic action/signature, not just full-End or simple members.
2. Voight, Quaternion Algebras v1.0.5, Proposition14.6.7: the actual statement and complete supplied proof on printed224 were read. Its ABHN/global classification and local splitting facts are imported; no whole-book proof audit is claimed. The embedding target is a field, not an integral order.
3. Milne, Abelian Varieties v2.0: Theorem2.9 is the complex-torus algebraicity/full-faithfulness input; Theorem1.1 is rigidity. The explicit Hodge-circle/quaternion argument is supplied above; no general Hodge conjecture is imported.
4. The normalized cyclic-cover eigensheaf and Cech calculations in(8)--(11), specialization and prime isolation are proved here. The old Newton node is used ONLY in Corollary4. PFP/globalU are not premises of the main obstruction.
