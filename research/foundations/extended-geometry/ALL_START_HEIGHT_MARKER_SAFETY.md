# All-marker safety for pure two-term starting-height phases

Date: 19 November 2026.

Work over any characteristic-zero field K. Put

    S=K[x,y,t]/(t²+2xy-1),
    Lambda(f)=CT_x (1/2) integral_(-1)^1
                    f(x,(1-t²)/(2x),t) dt.

The integral is its exact rational polynomial functional. Let a>=0 be
an INTEGER and consider a regular polynomial phase

    p=t^a(C(x)+t²D(x))/x,       C,D in K[x].              (S1)

Regularity is C(0)+D(0)=0. Purity means Lambda(p^m)=0 for EVERY m>=1.
Safety means that for every fixed polynomial marker G in S,
Lambda(Gp^m)=0 for all sufficiently large m, with the threshold allowed
to depend on G.

**Theorem AllStartHeightMarkerSafety.** Assume p is pure. The following
are equivalent:

1. p is safe for every fixed polynomial sphere marker.
2. Each of the TWO fixed markers x,x² has eventual mixed vanishing.
3. Either C(0)=D(0)=0 and x² divides both C,D, or
   C=c0,D=-c0 for some nonzero c0 in K.

For EVEN a, marker x alone suffices in statement2. Every other pure
phase in this class has one fixed marker x^j, 1<=j<=2, nonzero at
infinitely many positive orders. For odd a those orders are even.

The same dichotomy holds for every fixed ORIGINAL polynomial marker
of a standard balanced two-pair lift F: an unsafe phase has an original
witness among Z,Z². No claim is made that every mixed order is
nonzero, that a finite initial moment segment decides eventual safety,
or that this theorem classifies which nonconstant C,D are pure.

This is a complete safety classification CONDITIONAL ON ALL-ORDER
PURITY within(S1), including a=0. It does not classify phases with more
than the two height coefficients C+t²D.

## 1. The exact endpoint criterion

At a unit center scale p by C(0)^-1, so C(0)=-D(0)=1. Scaling by a
nonzero field element changes no vanishing or eventual-safety claim.
Set

 X=sqrt(-D/C), X(0)=1, A=C/X^a, u=x/A=xX^a/C=x+O(x²).

Let v(u)=u+O(u²) be the formal inverse. For every j>=0 and every m>=1
with am even, formal endpoint contact and residue substitution give

 Lambda(x^j p^m)=kappa_(a,m) [u^m] v(u)^j Phi(u),
 Phi(u)=T_a(v(u)), T_a(x)=X^(a-1)/(C u'),             (S2)

 kappa_(a,m)=integral_0^1 s^(am)(1-s²)^m ds
            =2^m m!/product_(q=0)^m(am+2q+1) !=0.

The factors are nonzero rational numbers, hence nonzero over K. If am
is odd the marked moment is zero by t parity. To recall the endpoint
argument, substitute s=Xt: a primitive of s^(am)(1-s²)^m has contact
order m+1 at both endpoints. Replacing endpoints +-X by +-1 therefore
changes no coefficient of x-order <=m. The factor x^j has nonnegative
order. The resulting coefficient [x^m]x^j A^m/X becomes the coefficient
in(S2) under the invertible formal substitution u=x/A.

For odd a, purity is exactly

    Phi(u)+Phi(-u)=2.                                  (S3)

For even a, kappa is nonzero at EVERY positive m, so purity is exactly

    Phi(u)=1.                                         (S4)

These are formal all-order identities, not finite-prefix tests.

## 2. Rationality needed by the argument

For odd a>=1, the retained OddStartRationalCenter theorem gives
u,X in K(x) over the ORIGINAL field K. Its hypotheses are exactly
unit center and purity in(S1); no polynomial-descent or rational-return
hypothesis is added here.

For even a, u=x(-D/C)^(a/2)/C is already rational. Equation(S4), after
substituting u(x), gives X^(a-1)=C u'. If a=0 this directly makes X
rational; if a>=2 divide by (X²)^((a-2)/2). Thus X is rational for all
even a as well. In every case T_a belongs to K(x). Write n=deg(u)>=1.

## 3. TwoMarkerRationalReturnRigidity

Let a be odd and suppose BOTH markers x,x² eventually vanish. Put
x+=v(u), x-=v(-u) and write their weights as T+=T_a(x+), T-=T_a(x-).
Purity gives T++T-=2. The even marked series are

 P(u²)=x+ T+ +x- T-,
 Q(u²)=x+² T+ +x-² T-.

By(S2), eventual vanishing is equivalent to P,Q being polynomials in
u². Indeed these are twice the even coefficients of the corresponding
marked series and every even kappa is nonzero. In particular P(0)=0.
A finite formal series is an actual polynomial identity in the formal
field; no convergence assumption is used.

Write x=x+, y=x-. The elementary two-point recurrence gives

 Q=(x+y)P-2xy,
 y=(Q-Px)/(P-2x).                                    (S5)

Here P,Q belong to K(u(x)²), a subfield of K(x). The denominator is
nonzero: P has x-order at least two whereas -2x has order one.
Thus y=j(x) is rational. The formal identities u(j(x))=-u(x),
j(0)=0,j'(0)=-1 and inverse-branch uniqueness imply j(j(x))=x.
Rational-map degrees multiply, so j is a nontrivial Möbius involution:

 j(x)=-x/(1+beta x), beta in K.

Put t0=xj(x), b0=x+j(x). Then b0=-beta t0, and(S5) gives

 Q=-(beta P+2)t0,
 t0=-Q/(beta P+2) in F=K(u²).                        (S6)

The denominator has constant2 and is nonzero. Therefore x satisfies

 x²+beta t0*x+t0=0 over F.

If n=deg u, the rational-map degree formula gives [K(x):F]=2n.
For completeness, if u=P0/Q0 with coprime P0,Q0, the polynomial
P0(X)²-wQ0(X)² is irreducible over K(w): it is primitive and linear
in w over K[X], and Gauss's lemma applies. Its X-degree is2n. This
proves the claimed field degree after w=u². The quadratic equation
forces2n<=2, hence n=1.

This argument constructs rational return only AFTER the two eventual
marker hypotheses. It assumes neither rational return nor canonical
polynomial descent for a general phase. In particular n>1 already
gives an unsafe ORIGINAL marker among x,x².

## 4. The degree-one odd-start case needs at most two markers

Assume a is odd and n=1. Normalization gives

    u=x/U, U=1+lambda x,
    C=U X^a, D=-U X^(a+2), T_a=U/X.                   (S8)

Polynomial D forces rational X to have no finite pole. Away from U=0
this is immediate; at its simple zero, a negative order h would give
ord(D)=1+(a+2)h<0. Hence X is a polynomial.

If lambda=0, equation(S3) becomes 1/X(x)+1/X(-x)=2. Therefore

    (2X(x)-1)(2X(-x)-1)=1.

Units of K[x] are constant. Together with X(0)=1 this forces X=1,
so C=1,D=-1 and p=2t^a y, the safe phase.

Let lambda!=0 and assume marker x is eventually annihilated. Put
v+=u/(1-lambda u), v-=-u/(1+lambda u) and H(u)=Phi(u). By(S3),
S_1=(v+-v-)H+2v-. Solving gives

 H=1-lambda u +(1-lambda²u²)S_1/(2u).                 (S9)

S_1 is an even polynomial with zero constant term, so(S9) makes H a
polynomial. A zero of X away from U=0 would make H=U/X have a pole
at a finite u-value, impossible. Over an algebraic closure all X-roots
therefore equal -1/lambda; normalization gives X=U^r. Polynomiality of
H forces r>=1. Now H=(1-lambda u)^(r-1), and(S3) forces r=1 or r=2:
for r-1>=2 the coefficient of u² is the nonzero number
2 binomial(r-1,2)lambda².

For r=1,

    S_1=2lambda u²/(1-lambda²u²),

which is not polynomial, contradicting eventual vanishing of x.
For r=2, S_1=0 but

    S_2=2u²/(1-lambda²u²),

so marker x² is unsafe. This completes the odd-start unit case.

This argument classifies safety, not all rational-coordinate fibers:
when marker x is unsafe, no assertion X=U^r is made.

## 5. Even starts need only marker x

For even a, equation(S4) makes the marked coefficient in(S2) equal
[u^m]v(u) for marker x, at EVERY positive m. If these eventually vanish,
v(u) is a polynomial. The identity u(v(u))=u is then an identity of
rational functions. Rational-map degrees multiply, so both maps have
degree1. The normalization v(0)=0,v'(0)=1 makes v=u and u=x.
Then C=X^a, and(S4) gives 1/X=1. Thus C=1,D=-1. The converse is the
one-sided safe phase. This includes a=0.

## 6. The zero center, including a=0

If C(0)=D(0)=0, all angular weights are nonnegative. Their weight-zero
part is c t^a+d t^(a+2), where c=C'(0),d=D'(0). The second and fourth
moments depend only on this core. If d=0, the second moment forces
c=0. If d!=0, put z=c/d and clear positive rational denominators in

 P2=z²/(2a+1)+2z/(2a+3)+1/(2a+5),
 P4=sum_(j=0)^4 binomial(4,j)z^(4-j)/(4a+2j+1).

For A2=(2a+1)(2a+3)(2a+5)P2 and
A4=product_(j=0)^4(4a+2j+1)P4, the exact Sylvester determinant is

 det Sylvester(A2,A4)=4096 P(a),
 P(a)=13312a^8+132096a^7+544768a^6+1209600a^5
      +1565824a^4+1195296a^3+518064a²+115200a+10125.

The polynomial has positive coefficients, so it is a nonzero integer
at every integer a>=0, including a=0. It remains nonzero in K.
Therefore P2 and P4 have no common root, a contradiction. Hence c=d=0
and x² divides C,D. Conversely that divisibility makes p strictly
positive in angular weight, unless p=0. Every fixed marker has finite
weight support, proving purity and eventual all-marker safety.

The checker verifies the determinant as a POLYNOMIAL identity in a,
not merely at finitely many starting heights.

## 7. Safety and unsafe markers in the original Image contraction

The constant unit phase 2c0 t^a y has weight-1; the zero-center phase
has strictly positive weights. Both are safe for every sphere marker.
For an original balanced lift F of degree d, split a fixed original
marker into bidegrees(A,B). If A>B, its contraction is zero. Otherwise
the output z-degree is B-A, independent of m, and output weights lie
in[-(B-A),0]. The contraction is equivariant under
xi1->s xi1,z1->s^-1 z1. Input weights shift by -m in the constant-unit
case and tend to +infinity in the zero-center case, so eventually lie
outside that fixed output interval. This proves safety for EVERY
original polynomial marker; it does not assert that all nonzero input
weights are annihilated.

Conversely the exact radial normalization for an ORIGINAL marker is

    E2(Z^j F^m)=(dm+j+1)! Lambda(x^j p^m),             (S10)

with any actual scalar clearing of p carried to its mth power. Every
factorial is nonzero. Thus the unsafe markers above transfer to the
SAME fixed original markers Z,Z². This proves the announced
original all-marker dichotomy.

## 8. The second monomial marker is necessary at every odd height

For every odd a>=1 and lambda!=0 put U=1+lambda x and take

 C=U^(2a+1), D=-U^(2a+5), X=U², u=x/U.

Then Phi(u)=1-lambda u. Every positive pure moment is zero, and

 v(u)Phi(u)=u,
 v(u)²Phi(u)=u²/(1-lambda u).

Thus marker x has ZERO mixed moments at EVERY positive order, while
marker x² is nonzero at every positive EVEN order2r, with value

 kappa_(a,2r)lambda^(2r-2) !=0.                       (S11)

The degree-two bound among powers of x cannot be replaced by degree
one. This says that the test with x alone fails; it does not assert
that every other single polynomial marker of a different form fails.
The a=1 degree9 member is a retained antecedent, not a new priority
claim.

An explicit original balanced lift, with R=xi1*z1+xi2*z2,
Z=xi1*z2, W=2xi2*z1, T=xi1*z1-xi2*z2 and U_h=R+lambda Z, is

 B=2WR^4-T²*(U_h^4-R^4)/Z,
 F_a=T^a U_h^(2a+1)B,       d=3a+6.                  (S12)

The quotient in B is a literal polynomial. Dehomogenizing R=1 gives
exactly the phase above, with no scalar clearing. The balanced degree
is d and ordinary total degree2d. Equation(S10) gives

 E2(Z F_a^m)=0 for all m>=1,
 E2(Z² F_a^(2r))=(2rd+3)! kappa_(a,2r)lambda^(2r-2) !=0.

These witnesses are radial-primitive and not proper powers, even up
to scalar over field extensions. At(xi1,xi2,z1,z2)=(1,1,1,-1),
R=0,T=2,Z=-1,W=2 and F_a=-2^(a+2)lambda^(2a+5)!=0.
The bilinears T and U_h are absolutely irreducible and nonassociated
(their coefficient matrices have rank2). Modulo T, B=2WR^4!=0,
so F_a has exact T-order a. Modulo U_h,

 B=R^4(2W+T²/Z)=R^6/Z=lambda^6 Z^5 !=0,

using T²+2ZW=R². Thus its exact U_h-order is2a+1. These orders are
coprime, excluding a scalar times any proper power. This proves the
stated original sharpness family without asserting a minimum degree
outside the specified class.

## 9. Scope and evidence

The source ledger pins the retained odd-start rationality theorem,
formal endpoint criterion, all-degree radial normalization and repaired
original-marker weight argument. The two-point recurrence, rational-return construction and field-degree
argument are supplied explicitly, so no normal-closure or polynomial-
descent theorem is imported for this new step.

A classical method antecedent is F. Pakovich, *On rational functions
orthogonal to all powers of a given rational function on a curve*,
arXiv:0910.2105v1 (2009), section5, Theorems5.1--5.2, printed pages21--22.
Theorem5.1 concerns complex rational P,Q and curve integrals
integral_gamma P^i Q^j Q': rational generating functions for deg(P)
consecutive j imply a common right factor and a closed image curve
with the specified pole-side property. Its proof uses inverse branches
and a Vandermonde determinant. Theorem5.2 treats eventual vanishing
for all i>=i0 and j>=j0 and obtains the corresponding composition and
pole-side criterion. Those curve, consecutive-marker and all-j
hypotheses differ from the normalized two-branch mass2 identity here.
Neither theorem is logically imported for our two-marker criterion.
The explicit recurrence and Mobius involution above avoid the
normal-closure/three-marker route entirely. Weighted inverse-branch
and composition arguments have this classical antecedent; no novelty
of those methods is asserted.

The standard-library checker verifies the two-point recurrence and rational
identities, the all-a Sylvester polynomial identity, exact negative
controls, raw interval moments across both start parities,
the two degree-one branches, and original-contraction normalization
controls. These controls supplement the unrestricted proof. They do
not establish purity from a finite prefix, or constitute formal proof
assistant or external human verification.

The [source ledger](all_start_height_marker_safety_source_ledger.json) pins the exact retained inputs and distinguishes the [classical method antecedent](sources/pakovich_rational_moments_0910_2105v1.pdf), Pakovich arXiv0910.2105v1 Section5, Theorems5.1--5.2, from the premises of this theorem. Those theorems concern complex rational curve integrals with different double-moment hypotheses; their inverse-branch/composition method informed the investigation, while the two-point recurrence above supplies the present result directly. The actual statements and adjacent proofs were inspected. No claim that moment-composition or inverse-branch methods are new is made.

The [checker](verify_all_start_height_marker_safety.py), [exact inputs](all_start_height_marker_safety_inputs.json) and [receipt](all_start_height_marker_safety_verification.json) bind finite algebraic controls to this written proof.

```bash
python3 research/foundations/extended-geometry/verify_all_start_height_marker_safety.py \
  --repo-root . --output research/foundations/extended-geometry/all_start_height_marker_safety_verification.json
```
