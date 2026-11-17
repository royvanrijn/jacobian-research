# Null suspension, polynomial coordinates and exact literal Gaussian images

Let K have characteristic zero, let E be a normalized Gaussian Wick
functional of nonsingular covariance on X, and let P in K[X]. Introduce
an INDEPENDENT Witt pair Z,W with

    E[Z^a W^b]=delta_(a,b) a!,
    Ptilde(X,Z,W)=P(X)+Z.                             (N1)

For a phase Q write Phi_Q[A]=E[A exp(tQ)], V_Q=Phi_Q[K[variables]],
and M_Q=K[t]V_Q. Original markers are spatial polynomials independent
of t. The null suspension has the exact embedded identities

    Phi_Ptilde[1]=Phi_P[1],
    V_Ptilde=M_Ptilde=M_P.                           (N2)

Moreover Ptilde is a polynomial coordinate and partial_Z Ptilde=1.
Consequently restricting universal Gaussian safety-to-purity or finite
vacuum-primitivity to polynomial coordinate phases does not narrow either
problem when two independent null variables may be added. The new
coordinate change is not asserted to preserve the Gaussian functional.

## 1. The exact original map and its degree filtration

For every original A(X) and nonnegative a,b, direct Wick contraction gives

    Phi_Ptilde[A Z^a W^b]=0, b<a;
    Phi_Ptilde[A Z^a W^b]
       =b!/(b-a)! t^(b-a) Phi_P[A], b>=a.             (N3)

Indeed the only surviving auxiliary term of exp(tZ) has degree b-a,
and its factorial contraction is b!. This is an equality of formal
series: at every time order only finitely many polynomial contractions
are used. In particular Phi_Ptilde[A W^k]=t^k Phi_P[A]. Every original
marker is a finite sum of the separated markers in (N3), proving (N2)
including the literal equality. The distinction is material: V_P itself
need not equal M_P.

Let V_P,D be the literal image of markers of total spatial degree at most
D, and put V_P,D=0 for D<0. Then, for every D>=0,

    V_Ptilde,D=sum_(k=0)^D t^k V_P,D-k.               (N4)

For the forward inclusion a surviving term in (N3) has k=b-a and costs
a+b>=k auxiliary degrees. For the reverse inclusion use A(X)W^k, whose
degree is deg A+k. Thus (N4) retains the ORIGINAL filtered map and is not
an assertion only about a rationally localized module.

Writing F_P,D=K[t]V_P,D gives the further exact identity

    F_Ptilde,D=F_P,D.                                 (N5)

Every summand on the right of (N4) belongs to F_P,D, while k=0 gives
the reverse inclusion after taking K[t]-spans. The degree cost needed
to turn a polynomial-time amplitude into a single original marker is
therefore explicitly bounded by its time degree in (N4).

The submersion Stein lift is especially simple here. The polynomial
operator TQ=WQ-partial_Z Q satisfies

    Phi_Ptilde[TQ]=t Phi_Ptilde[Q].                   (N6)

This follows by Gaussian integration by parts in Z. On markers independent
of Z it is just multiplication by W. It is the constant-gradient instance
of the [general polynomial submersion lift](GAUSSIAN_ANGULAR_CORE_LITERAL_REGULAR_SECTOR.md#1-a-polynomial-submersion-makes-the-literal-image-time-closed).

## 2. Polynomial-coordinate reduction preserves the actual obstruction

The map fixing X,W and sending Z to Z+P(X) is a polynomial algebra
automorphism, with inverse Z to Z-P(X). Thus Ptilde is literally one
coordinate of an ambient polynomial coordinate system. This is an
ordinary polynomial automorphism; Gaussian preservation is unnecessary
for (N3), and generally fails for this automorphism.

The same embedded M and same original vacuum in (N2) imply equivalences
between P and Ptilde for each of the following properties:

* safety, meaning every original marked period is polynomial;
* purity, meaning the original vacuum is one;
* finiteness of M over K[t];
* primitivity, meaning M=Phi[1]K[t];
* any specified scalar condition on Phi[1], including p(t)exp(psi(t)).

Their faithful generic images and distinguished vacuum also agree.
Hence each of the following universal assertions is equivalent to its
restriction to polynomial coordinate phases:

    every safe Gaussian phase is pure;
    every finite actual image has primitive original vacuum;
    every finite actual image with a specified scalar unit condition
        has primitive original vacuum.               (N7)

Here the quantification permits all Gaussian dimensions and the
nonsingular covariances obtained by adding one Witt block. In one
direction apply the restricted assertion to Ptilde; in the other take
a subclass. In particular unimodularity of the polynomial gradient
alone cannot reduce these universal problems. No fixed-dimension or
real-positive-covariance reduction is asserted.

Over C the Witt block may be expressed using independent standard
ordinary variables U,V by Z=U+iV and W=(U-iV)/2. Thus the reduction
also stays within standard complex Gaussian coordinates after this
linear change. The covariance remains nonsingular. Over a field without
i, the displayed split rational Witt covariance already gives (N1).

The polynomial critical set of Ptilde is empty, but Gaussian stationary
actions are not thereby removed. In the split auxiliary coordinates
their phase is tP(X)+tZ-ZW plus the old Gaussian quadratic form. Its
auxiliary critical equations are W=t and Z=0, and its auxiliary critical
action is zero. Thus the old Gaussian stationary equations and actions
are preserved exactly.

## 3. Null-variable elimination and all safe recompositions

There is a useful extension of (N2). Let G(X,Z) be ANY polynomial,
independent of W, and put G0(X)=G(X,0). Then

    Phi_G[1]=Phi_G0[1], M_G=M_G0.                    (N8)

For completeness define polynomials H_k(t,X,Z) recursively by

    H_0=1, H_(k+1)=partial_Z H_k+t G_Z H_k.

Thus partial_Z^k exp(tG)=H_k exp(tG). Wick contraction now gives,
with k=b-a>=0,

    Phi_G[A Z^a W^b]
       =b!/k! Phi_G0[A H_k(t,X,0)],                 (N9)

and gives zero when b<a. The right side is a polynomial-time amplitude
in the original X variables, hence belongs to M_G0. Conversely the old
markers A(X), with auxiliary marker one, give exactly Phi_G0[A]. This
proves both inclusions of (N8) without replacing the actual source by a
different scalar solution. It does not assert V_G=M_G for arbitrary G.

This is also a direct instance of the
[angular-core full-image theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md):
give Z,W weights 1,-1 and all old variables weight zero, so that G-G0
is strictly positive. Equations (N3)-(N5) add an exact literal degree
calculation for the linear suspension, rather than a second account
of that general module theorem.

Applying (N8) to G(X,Z)=g(P(X)+Z), for EVERY g in K[T], gives

    Phi_(g(Ptilde))[1]=Phi_(g(P))[1],
    M_(g(Ptilde))=M_(g(P)).                          (N10)

In particular the entire safe-recomposition ideals are identical:

    S_Ptilde=S_P,
    S_P={g in K[T]:g(P) is safe}.                    (N11)

This uses equality of all original marked modules, not a scalar-purity
test. The [safe recomposition classification](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
then assigns the same zero ideal or the same unique safe shift to both.

## 4. Exact controls and the boundary of a scalar-only hypothesis

For P=0, one has V_P=K but V_Ptilde=K[t], and (N4) gives exactly
V_Ptilde,D=K[t] of degree at most D. This includes the constant-phase
exception to literal V=M and shows why the suspension can enlarge V.
For P=S^2 in an old Witt pair, (N2) preserves its rank-two nonfinite
global image, while making it literal, consistently with the
[explicit angular example](GAUSSIAN_ANGULAR_CORE_LITERAL_REGULAR_SECTOR.md).

Scalar polynomiality plus unimodular gradient does not force purity.
The [arbitrary polynomial-seed construction](POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md#1-arbitrary-polynomial-unit-periods-in-seven-variables)
already has partial_A P=1: its only A term is A. For every p(0)=1,
its M=1 instance has original vacuum p and full image K[t,p^-1]. A
nonconstant p therefore gives an explicit polynomial submersion with
nonpure polynomial vacuum, and with an UNSAFE, nonfinite full image.
For example p=1-t has the original bad marker (1+Y)^2, of period
(1-t)^-1. No safe or finite-image counterexample is claimed. An
independent exponential-cancelling block from the
[reversible vacuum construction](SAFE_SEED_DIVISIBILITY_AND_PHASE_MARKER_FILTRATION.md#4-an-exact-auxiliary-gaussian-exponential-cancellation)
can replace p by p exp(psi) without removing the constant A derivative.

The companion checker
`verify_gaussian_null_suspension_and_coordinate_reduction.py` checks
bounded original Wick identities, the exact filtered linear map,
the lift, coordinate inverse, and nonlinear elimination/recomposition
identities. Its receipt pins the proof and every linked mathematical import.
The all-degree statements are proved above. No independent replay,
external review, formal verification or priority claim is made.
