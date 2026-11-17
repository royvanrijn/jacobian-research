# Boundary images of arbitrary Gaussian invariant cores

Let K have characteristic zero. Use an original Witt Gaussian pair Z,W
with E[ZW]=1 and an independent ordinary Gaussian Y with E[Y^2]=1.
Put S=ZW and let P(S,Y) be ANY polynomial, with boundary G(Y)=P(0,Y).
There is no degree, parity or radial-degree restriction. Write

    Phi_P[A]=E[A exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

The exact boundary identity is, for every polynomial A(S,Y),

    Phi_G[A(0,Y)]
       =Phi_P[A-partial_S A-t P_S A].                (B1)

The right side is a polynomial-TIME original amplitude. In particular
the entire ordinary-Gaussian boundary image embeds as

    M_G subset M_P.                                  (B2)

If the full original M_P is finite over K[t], then

    G(Y)=c+bY,
    M_P=M_G=exp(ct+b^2t^2/2)K[t].                    (B3)

The original boundary amplitude 1-tP_S gives this exponential generator.
The distinguished original vacuum may still be p(t) times that generator;
no general assertion p=1 is made. In particular (B3) is a NECESSARY
finite-image condition. An affine boundary alone is not sufficient.

## 1. The original endpoint identity, with all markers retained

For the factorial radial moments E_S[S^n]=n!, direct monomial comparison
gives

    E_S[F-partial_S F]=F(0).                         (B4)

For a positive monomial degree the factorial terms cancel; for the
constant monomial the value is one. Apply (B4), coefficientwise, to
F=A(S,Y)exp(tP(S,Y)), and then apply the original Y Gaussian functional.
At every time order these are finite polynomial operations. The chain
rule gives exactly (B1), without a convergence assumption or a choice
of a different scalar solution.

The general identity is the endpoint term of the retained
[radial Gamma Stein formula](GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md#1-the-actual-radial-scalar-source-retains-its-endpoint).
Its affine-radial special case is the q=1-tH boundary amplitude in
[the affine radial source](GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md#1-exact-gamma-elimination-preserves-the-original-endpoint).
The conclusion here uses the endpoint to embed the FULL original boundary
image for arbitrary P. It does not replace the relative source by an
absolute quotient or declare the endpoint term zero.

For a fixed ordinary marker F(Y), take A=F in (B1):

    Phi_G[F]=Phi_P[F]-t Phi_P[P_S F].                (B5)

Both spatial markers on the right are original polynomials. Extending
K[t]-linearly proves (B2), including every boundary marker and every
polynomial-time combination. It is not a claim that the literal fixed-
marker space V_G is contained in V_P. The displayed transport itself
has time degree at most one.

There is also exact parameter compatibility. Set B_P=1-partial_S-tP_S
on polynomial-time amplitudes. Then

    [partial_t+P,B_P]=0.                              (B6)

Indeed the derivative -P_S of B_P cancels [P,-partial_S]=P_S.
Evaluation at S=0 sends multiplication by P to multiplication by G.
Thus the full boundary identity retains the original parameter
connection as well as every scalar value.

If A has original spatial degree d_A after S=ZW and P has degree d_P,
then partial_S lowers degree by two when nonzero. Hence the amplitude
in (B1) has spatial degree at most

    d_A+max(0,d_P-2),                                (B7)

and time degree at most one. If P_S=0 its time degree is zero. The same
bound with d_A=deg F applies to (B5). No inverse or rational spatial
marker is introduced.

## 2. A finite exponential image has no proper differential submodule

We use the retained
[general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image):
every nonzero finite original Gaussian image has the unique normalized
form

    M_P=exp(psi)K[t], psi in tK[t],
    Phi_P[1]=p exp(psi), p in K[t],p(0)=1.           (B8)

Let N be any nonzero K[t]-submodule of exp(psi)K[t] stable under time
differentiation. The set I=exp(-psi)N is a nonzero ideal of K[t]. It
is stable under ordinary differentiation: if q is in I, then
exp(psi)(q'+psi'q) is in N, as is exp(psi)psi'q, so q' is in I.
If q generates I, then q divides q'. In characteristic zero a nonconstant
q has derivative of strictly smaller degree. Thus I=K[t], and

    N=exp(psi)K[t].                                  (B9)

Original Gaussian images are derivative-stable because
partial_t Phi_H[F]=Phi_H[HF], with polynomial-time coefficients handled
by the product rule. The nonzero boundary image M_G in (B2) therefore
qualifies for (B9): its vacuum has constant coefficient one. Consequently

    M_P finite ==> M_G=M_P.                         (B10)

This argument proves finiteness of M_G as well. It does not assume that
the distinguished boundary or ambient vacuum generates either module.

The same identity (B1), inclusion (B2), and finite equality (B10) hold
with Y replaced by ANY independent Gaussian vector of nonsingular
covariance. No classification of that vector boundary is used. The next
explicit affine and zero-boundary conclusions concern ONE ordinary
variable; a nonzero isotropic linear vector phase can itself be safe.

## 3. The explicit one-ordinary-variable consequence

The [one-variable finite-image theorem](GAUSSIAN_W_AFFINE_ARBITRARY_FIBRE_FINITE_IMAGES.md#one-ordinary-gaussian-variable-finite-image-means-affine-phase)
states that M_G is finite exactly when G=c+bY, and then

    Phi_G[1]=exp(psi_G), psi_G=ct+b^2t^2/2,
    M_G=exp(psi_G)K[t].                               (B11)

Its proof retains the original absolute ordinary Gaussian source. It
excludes degree at least two by the stated leading-coefficient comparison;
it does not discard a Gamma endpoint or apply a quadratic-only theorem
to a cubic or higher polynomial.

Applying (B10)-(B11) proves (B3). Equivalently, inclusion of the two
normalized finite exponential generators gives
exp(psi_G-psi)=q in K[t]. Differentiating q'= (psi_G'-psi')q and comparing
polynomial degrees forces q=1 and psi_G=psi. Both proofs preserve the
same embedded original images, rather than only their generic ranks.

The actual amplitude 1-tP_S has period Phi_G[1]=exp(psi_G), so it
explicitly realizes the generator. Nevertheless (B8) still allows

    Phi_P[1]=p exp(psi_G),
    Phi_P[P_S]=((p-1)/t)exp(psi_G).                  (B12)

Here (p-1)/t is a polynomial. Neither (B1) nor the existence of this
generator amplitude forces p=1. General finite-vacuum primitivity and
safety-to-purity for arbitrary mixed-parity invariant cores are not
resolved by this boundary reduction.

For example P=S has G=0, but

    Phi_P[S^m]=m!/(1-t)^(m+1).

These pole orders are unbounded, so M_P is not finite over K[t]. Its
boundary amplitude 1-t still has period one. This proves that even
the zero boundary is not a sufficient finite-image condition.

## 4. Safety and strict one-sided additions

Safety means every original fixed-marker period is a polynomial. Then
M_P is a nonzero submodule of K[t], hence finite, and (B3) applies.
Its exponential generator is a polynomial only when ct+b^2t^2/2=0.
Since K is a field and Y has variance one, this forces c=b=0. Conversely
if M_P is finite and G=0, (B3) makes M_P=K[t], so every original marker
is safe. Thus the exact criterion is

    P is safe iff M_P is finite AND P(0,Y)=0.        (B13)

In particular every safe such invariant core is divisible by S. This
is only a necessary algebraic condition on its coefficients; the
finiteness condition in (B13) remains essential.

Now add a polynomial R with strictly positive Z-minus-W weights, or
strictly negative weights after reversing the pair. The
[angular-core full-image theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md)
gives M_(P+R)=M_P and preserves the vacuum. All finite-image, exponent
and safety assertions above therefore transfer with the SAME boundary
G. The explicit identities (B1) and (B5) also remain valid with Phi_(P+R)
on their right: every amplitude displayed there has weight zero, so
all terms containing R vanish by angular balance. The degree cost is
still controlled by P, independently of the degree of R.

No fixed spatial original-marker transport, full image-kernel formula,
finite-if-boundary-affine converse, universal primitivity theorem or
general mixed-parity purity-to-safety theorem is asserted.

## 5. Several independent radial variables

Let S_i=Z_iW_i for n independent Witt pairs, let Y be an independent
Gaussian vector, and let P(S_1,...,S_n,Y) be arbitrary. Repeated endpoint
evaluation gives

    M_(P(0,...,0,Y)) subset M_P.                     (B14)

There is an exact finite original amplitude for it. Set

    B_i=1-partial_Si-tP_Si.

The B_i commute: the twisted derivatives partial_Si+tP_Si commute
because mixed partial derivatives of P agree. Iterating (B4) yields

    Phi_(P(0,...,0,Y))[F(Y)]
       =Phi_P[B_1 ... B_n F(Y)].                    (B15)

Each endpoint is retained, and product coefficients contain all required
mixed derivatives. The amplitude has time degree at most n. This is
not a product of arbitrary fixed-marker substitutions. At any intermediate
step the remaining phase and its derivatives are restricted to the
endpoint already taken, exactly as the formula requires.

If M_P is finite, (B9) upgrades (B14) to equality, without classifying
the Gaussian-vector boundary or assuming its vacuum is primitive.
For ONE ordinary Y the boundary must again be c+bY, and the common finite
image is exp(ct+b^2t^2/2)K[t]. In that scalar case (B15) with F=1
explicitly realizes its generator. It does not establish that the
original full vacuum is that generator. The finite exponent has degree
at most two, independently of radial degrees or the number of pairs.

## Evidence and replay

Run

    python verify_gaussian_invariant_core_boundary_reduction.py \
      --output gaussian_invariant_core_boundary_reduction_verification.json

from this directory. The
[checker](verify_gaussian_invariant_core_boundary_reduction.py) compares
bounded ORIGINAL Wick contractions for nonlinear mixed-parity P and
S,Y markers, checks the parameter commutator and degree cost, retains
strict remainders, and checks explicit boundary generators and the
nonfinite zero-boundary example. Its
[receipt](gaussian_invariant_core_boundary_reduction_verification.json)
pins every linked mathematical import. The unrestricted identity and
module consequences are the written proofs above. No independent
replay, formal verification, external review or priority claim is made.
