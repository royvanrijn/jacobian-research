# Full filtered Gaussian safety without a purity premise

16 November2026. Internal all-degree/all-ring extension of the sealed
[filtered critical-action theorem](FILTERED_GAUSSIAN_CRITICAL_ACTION_AND_EXACT_ALL_RING_GATES.md).
No earlier claim or packet is changed. Independent replay,formal verification
and external review remain false;priority unassessed. The
[finite gate](filtered_gaussian_full_safe_gate.py),
[checker](verify_filtered_gaussian_full_safe_gate.py),
[certificate](filtered_gaussian_full_safe_gate_certificate.json) and
[receipt](filtered_gaussian_full_safe_gate_verification.json) retain controls.

## 1. Entire supplied architecture and retained exact source formulas

Let Q0 be a fixed nonsingular rational Gaussian quadratic form of weight2
for supplied NONNEGATIVE rational coordinate weights. Let P have weighted
degree AT MOST2,with arbitrary polynomial coefficients in weight-zero
variables. In the canonical split frame,

    Q0=Z.W+(1/2)T^T B0 T,
    P=W.H(Z)+F(Z,T), 0<w(T_i)<2.                       (1)

Constant and lower-weight phase terms are included. This is the ENTIRE
filtered architecture already proved,not the assertion that every Gaussian
phase has such a frame. The split covariance becomes standard after the
usual coefficient extension by i;no nonzero real pure phase is asserted.

The previous finite block/Morse construction supplies

    g=tH(g), J=det(I-tH_Z(g)), sigma_g=sigma(t,g),
    psi(t,Z)=nu(t,Z)/sigma(t,Z)^r, psi(0,Z)=0,
    f=Phi[1]=exp(psi_g)/(J sqrt(sigma_g)),
    Phi[Z^alpha]=f*g^alpha.                            (2)

Here sigma,nu are POLYNOMIALS,sigma(0,Z)=1,nu(0,Z)=0,r>=0;
all inverse/square-root branches are normalized at zero. These data come
from a finite polynomial change in the ORIGINAL potential with Jacobian1.
Every original residual and W-momentum period is f times a finite expression
in g,J^-1,sigma_g^-1 and the retained critical/Morse/source derivatives.
Only the exact formulas from the old proof are used below.

## 2. Safety forces purity over every characteristic-zero FIELD

Without any purity premise,the exact field classification is

    P safe iff g in K[t]^n, J=1, sigma_g=1, psi_g=0.     (3)

In particular EVERY safe phase in this architecture has ALL positive pure
moments zero. This removes the old purity premise from its field safe
normal-form theorem. It does not establish safety implying purity outside
this architecture or its Gaussian moment-preserving polynomial symmetry orbit.

Necessity has three short steps. If P is safe,f is a nonzero polynomial,
with f(0)=1. Each f*g_i=Phi[Z_i] is polynomial too,so g_i is RATIONAL.
All powers of the ORIGINAL coordinate give f*g_i^m in K[t] for every m.
At any finite pole of g_i,this product has valuation

    ord(f)+m ord(g_i),

which becomes negative. Hence g_i has no finite pole and is a POLYNOMIAL.
The argument may be made over an algebraic closure and descends to K.
It is essential to use all powers,not only the linear coordinate periods.

Now J and sigma_g are polynomials,and psi_g is rational. Equation(2)
makes exp(psi_g)=f*J*sqrt(sigma_g) algebraic. On the proper algebraic
curve of the square-root extension,the logarithmic differential of an
algebraic function has at most simple poles. At any pole of nonconstant
psi_g,its differential has pole order at least two. Thus psi_g is constant,
and its zero initial value makes it0. Equivalently this is the retained
elementary algebraic exponential obstruction,not a new transcendence import.

It follows that f^2 J^2 sigma_g=1 in K[t]. Every factor is a polynomial
unit. Their constant coefficients are1,so f=J=sigma_g=1. This proves
necessity in(3),including purity. Sufficiency is the full original-amplitude
formula in(2):all denominators are1,the seed is1,and g is polynomial.
There is no need to compute an actual-image kernel or an ambient Hom rank.

This also applies to nonnegative-weight homogeneous phases,with psi=0.
It includes the polynomial-seed phases of the previous turn. Their
nonconstant polynomial seed is consistent with(3) because an original
coordinate g has a pole and its sufficiently high powers reveal it.
Their PURE POWER phases have weight greater than2 and escape(1).
The dimension-only marker-test obstruction is unchanged.

## 3. A finite coefficient ideal for the full ALL-RING safety functor

Let R be ANY commutative Q-algebra. Define

    D_i=max(0,deg_Z H_i), B=product_i(D_i+1),
    q=g mod t^(B+1),
    J_q=det(I-tH_Z(q)), sigma_q=sigma(t,q), nu_q=nu(t,q).

For no zero positions the empty product is1 and q is empty. The old
Bezout graph argument proves that any polynomial field fixed point has
degree at most B. Retain the following FINITE ORDINARY coefficient ideal:

    I=(coeff_t(q-tH(q)), coeff_t(J_q-1),
          coeff_t(sigma_q-1), coeff_t(nu_q)).            (4)

Over fields,(3) says exactly I=0. Over EVERY R,the stronger result is

    P safe iff I is NILPOTENT,
    polynomial-Hom normalization = Gamma_I(R)
       ={c : I^L c=0 for some L}.                      (5)

No purity premise occurs. I is a finite,possibly nonradical presentation
of the safety topology;it is not asserted equal to an ordinary pure ideal,
an individual cutoff ideal or the radical universal safe ideal.

For necessity,safety descends to every prime residue field. The field
degree bound gives g=q there,and(3) makes every generator in(4) zero.
Thus the finitely generated ideal I is nil,hence nilpotent,without
Noetherianity or a bounded nilradical hypothesis.

For constructive sufficiency,assume I^L=0. Modulo I,the polynomial q is an
exact fixed point and its Jacobian determinant is1. Its polynomial inverse
Jacobian lifts through I by the finite geometric series. Apply polynomial
Newton correction to F(z)=z-tH(z),starting at q. The residual lies in
I R[t]^n and squares its ideal power on each step;ceil(log_2 L) steps
suffice. Inverting each updated polynomial Jacobian remains finite because
its determinant is1 modulo I. All corrections have positive time order.
The resulting polynomial p solves p=tH(p),p(0)=0,so formal uniqueness
identifies p=g. This is the same nilpotent Newton mechanism as the old proof.

Since g=q modulo I,J_g and sigma_g are polynomial units congruent to1
modulo I,and nu_g has coefficients in I. Therefore psi_g is a POLYNOMIAL
with nilpotent coefficient ideal. Its exponential is a FINITE polynomial;
the normalized inverse square root of sigma_g and J_g^-1 are finite too.
Every original residual/Morse/source coefficient in(2) is now polynomial.
This proves all-marker safety directly,and its full actual image is R[t]
because the polynomial unit f has a polynomial inverse.

For normalization,apply safety over R/Ann(c). Every original series
becomes polynomial in that quotient precisely when its high coefficients
are killed by c. The finite gate(4) commutes with that coefficientwise
base change. Criterion(5) there is exactly I^L c=0,proving the Gamma formula.

Safe phases over R need not be pure:their positive pure coefficients are
nilpotent defects. More precisely,f is a polynomial unit congruent to1
modulo I,so the ideal of its positive coefficients is nilpotent. The old
finite logarithmic identity is still needed if actual purity over R is
requested;one cannot replace it with the field equations.

## 4. Why the new determinant and action blocks are necessary

The old residual-only finite criterion REQUIRED full purity. The example
H=e+Z over Q[e]/(e^2) has residual ideal(e),but g=et/(1-t) is not
polynomial. In(4),J_q=1-t contributes the coefficient-1,so the NEW full
ideal is the unit ideal and correctly rejects safety. No previous gate
or failed example has been erased.

Conversely,H=1+eZ with no residual F gives I=(e). Over e^4=0,its exact
polynomial fixed point is t+et^2+e^2t^3+e^3t^4,exceeding the field
bound B=2. Nilpotent lifting,not truncation at B,is essential.

For the genuine three-variable family

    P=W(1+eZ)+e(T^2-1)/2, E[ZW]=E[T^2]=1,

the full gate is(e). Over e^2=0,it is safe and nonpure:

    Phi[1]=1+et, Phi[Z]=t+2et^2,
    Phi[W]=0, Phi[T^2]=Phi[ZW]=1+2et.

The focused checker verifies these by25 ORIGINAL Wick contractions,not
just by reading a normal-form flag. Eight gate models,a two-coordinate
fixed point,two polynomial Newton steps,five pole-valuation controls and
four input rejections transcribe the finite construction. These controls
do not establish the universal quantifiers or provide an ideal oracle
for arbitrary abstract coefficient rings. In finitely presented rational
rings,the existing primary/saturation backends can consume(4).

The new constructor also normalizes each critical block residual by exact
polynomial expansion BEFORE testing which variables occur. The retained
older routine could reject a legitimate coefficient such as -(1+Z)XY
because an unexpanded cancelling expression still listed X among its free
symbols. The new pure-unsafe model exercises that case. The old sealed
code/receipt and its written theorem remain unchanged;this is a symbolic
presentation fix in the new constructor,not a change to the block proof.

## 5. Published Stein operators: useful identities and the exact boundary

[Azmoodeh--Gasbarra--Gaunt,On algebraic Stein operators for Gaussian
polynomials](https://arxiv.org/html/1912.04605v4),published Bernoulli29(1)
(2023),350-376,DOI10.3150/22-BEJ1460,provides an algebraic operator-chain
and bounded null-control construction. We read Definition1.1,Section2's
IBP setup,Section3.1 Definitions3.1-3.2/Propositions3.1-3.2 and their proofs,
and Section3.3's control correspondence. Its distinction between algebraic
and general Stein operators is retained. This is useful provenance for
Gaussian IBP/operator searches;no scalar Stein operator is substituted
for every ORIGINAL marker image here. The source does not resolve our
general safety-to-purity question. Its probabilistic/conditional-expectation
statements are not imported into arbitrary complex or nonreduced inputs.

The proof(3)-(5) uses the older exact filtered Gaussian source construction,
the polynomial pole argument and elementary nilpotent descent. Classical
Morse/Gaussian/Lagrange/Newton credit remains unchanged. Sabbah/Maisonobe
retain convergent rational-slope credit;nonconvergent formal and all-ring
extension priority remains separately unassessed. No historical novelty
claim is made for the present extension.

## 6. Remaining general obligation

This closes the full supplied filtered class WITHOUT a purity premise,
over fields and all coefficient rings,including its Gaussian-preserving
polynomial symmetry orbit. The all-degree original-marker-filtration
obstruction from the previous turn remains valid outside this class.

Still OPEN: field safety implying purity for phases admitting no such
frame;exhaustive safe phase normal forms retaining their source filtration;
and an ORIGINAL regular-transcendental Gaussian realization or a general
obstruction. The field proof needs a seed times an algebraic fixed-point
field and all coordinate powers. The Kloosterman/Airy families show that
this is not a universal normal form. No finite dimension/degree campaign,
terminal UNKNOWN restart,delegation,outreach,assurance promotion or broad
goal completion is undertaken.
