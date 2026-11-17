# Actual scalar bridge and finite-image exclusion at a simple origin

Internal written proof. The scalar bridge below is independent of an
original-source localization isomorphism. First work over C; the final
exclusion descends to every characteristic-zero field. Let

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

with constant a and polynomial other coefficients, and suppose

    ell=zL, h=zH0, c=zC0, L(0)!=0.

All Gaussian coefficients use the original Witt pair Z,W and an
independent ordinary Gaussian Y. Put f_P(t)=E[exp(tP)].

There are two logically distinct assertions.

1. If the full actual image is finite, its explicit origin period is
   its primitive exponential generator: necessarily a=0 and
   psi=k(0)t+u(0)^2t^2/2. This alone does NOT prove primitive vacuum.
2. A separate explicitly constructed horizontal functional on the
   localized curve agrees with the ACTUAL scalar vacuum and all its
   parameter-covariant derivatives. With deg ell>=2 and every root of
   gcd(ell,c^2) of multiplicity<=2, the retained positive-volume
   comparison excludes any identity f_P=p exp(psi) with nonzero
   polynomial p and psi in tC[t]. Hence no finite actual image exists
   in this simple-origin volume range.

The second assertion does not transport the original full marker map
through localization. In particular its value on the Gysin marker can
differ from the original value. The proof below constructs the needed
cyclic rational quotient directly and checks every parameter iterate.

## 1. What finite-image classification and the actual kernel already force

The exact actual homotopy with

    F=L(z)W^2+[H0(z)+C0(z)Y]W,
    Jker=1-t partial_W F

gives, for every W-free polynomial A,

    Phi_P[Jker A]=E_Y[A(0,Y) exp(tQ0(Y))],
    Q0=aY^2+u(0)Y+k(0).

Here Jker is a polynomial-time original amplitude. In particular

    N(t)=Phi_P[Jker]
      =d^(-1/2) exp(R(t)),
    d=1-2at, R=tk(0)+t^2u(0)^2/(2d), N(0)=1.

If the FULL actual image is finite, the
[general finite-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives M=exp(psi)C[t], psi in tC[t]. Therefore N=exp(psi)q with polynomial
q, q(0)=1. Squaring gives

    q^2 d=exp(2(R-psi)).

The exponential of a rational function is rational only when that
function is constant: differentiate logarithmically and compare pole
orders, then compare at infinity. Thus R=psi and q^2d=1. A polynomial
product equal to one gives q=d=1. Consequently a=0 and the claimed
psi follows. The explicit origin period N itself supplies the normalized
generator of any finite actual image. This scalar assertion makes no claim about
the rank or injectivity of the original localization kernel. The vacuum
can still have form f_P=exp(psi)p; no assertion p=1 follows from this
kernel argument alone.

## 2. An exact actual scalar Wick-residue identity

Use d=1-2at and set

    Psi0(t,z)=t k(z)+t^2u(z)^2/(2d), D=exp(Psi0),
    A=t ell+t^2c^2/(2d)=t z Ltilde,
    Ltilde=L+t z C0^2/(2d),
    B=t h+t^2c u/d, z-B=zF,
    F=1-tH0-t^2C0u/d,
    H=2t Ltilde/F^2.

Exact completion of the ordinary Y Gaussian gives

    f_P=d^(-1/2) E_ZW[D exp(AW^2+BW)].

Expanding W and applying the Witt rule proves, coefficientwise,

    f_P=d^(-1/2) sum_(j>=0) (2j)!/j!
                         Res_z D A^j/(z-B)^(2j+1)
       =d^(-1/2) sum_(j>=0)(2j-1)!!
                         Res_z (D/F)H^j/z^(j+1).                    (1)

For example the expansion of (z-B)^(-2j-1) contributes
binom(2j+k,k), giving exactly (2j+k)!/(j!k!), the original W-Wick
coefficient. The cancellation uses only the displayed powers of z.
All denominators are expanded t-adically BEFORE residue; each fixed
t coefficient is a finite Laurent polynomial in z.

Let z=z(t,S) solve z=H(t,z)S^2. Then

    z=2L(0)tS^2+O(t^2), J=1-H_z(t,z)S^2=1+O(t).

The formal inverse exists uniquely, with polynomial S coefficients.
The elementary residue change w=z-yH(t,z), a t-adically near-identity
map, gives

    sum_j y^j Res_z b(z)H(z)^j/z^(j+1)
      =b(z(t,y))/(1-yH_z(t,z(t,y))).

Indeed its left side is Res_z b(z)/(z-yH(z))dz, and changing to w
extracts the constant coefficient of b(z(w)) dz/dw. The change of
formal residue is coefficientwise valid; invariance also follows by
differentiating the near-identity change and taking the residue of a
total z derivative. Apply y=S^2 and ordinary Gaussian contraction to
obtain the exact ACTUAL scalar identity

    f_P=d^(-1/2) E_S[D/(FJ)]_(z=z(t,S)).                             (2)

Only ordinary polynomial Gaussian contraction occurs in (2).

## 3. A full localized functional with a fixed Laurent prescription

Over K=C(t), take the smooth affine curve algebra

    Acurve=K[z,U,U^-1]/(U^2-q),
    q=c^2+2d ell/t=2dA/t^2,
    phi=Psi0-(z-B)^2/(4A)=Psi0-z/(2H),
    eta=dz/(tU).

The phase is regular in Acurve because A is a nonzero scalar multiple
of U^2. A disconnected normalization causes no problem. Since U is
invertible, every relative one-form is uniquely f dz. Fixed-z and
fixed-t derivatives extend by

    delta_t U=q_t/(2U), delta_z U=q_z/(2U),

and preserve Acurve. Thus EVERY finite parameter-covariant iterate
of every allowed form remains in this same algebra, not a larger
uncontrolled localization.

On the selected formal root define

    U(t,S)=sqrt(d) F(t,z)H(t,z)S/t=2L(0)S+O(t),
    Psi(t,S)=Psi0(t,z(t,S)).

Its square is q, and sqrt(d) has constant coefficient one. Therefore
the full curve algebra and its forms map to

    B=C[S,S^-1]((t)), respectively B dS.

In particular U^-1 has leading coefficient (2L(0)S)^-1. No inverse
of a nonmonomial Laurent polynomial in S is introduced. Each element
has a finite lower t bound and finite Laurent coefficients. Also
Psi in tC[S^2][[t]], so exp(Psi) acts on B coefficientwise.

Define a FIXED linear Laurent Gaussian G by

    G(S^(2m+1))=0 for every integer m,
    G(S^(2m))=(2m-1)!! for m>=0,
    G(S^(-2n))=(-1)^n/(2n-1)!! for n>=1.

It obeys G[(partial_S-S)b]=0 for every b in B, and commutes with
partial_t. This is a coefficientwise Laurent prescription; no pole
prescription is reset as t varies and no integral at S=0 is asserted.
For a pulled-back form omega=b(t,S)dS set

    Lambda(omega)=(1/2)G[exp(Psi)b].                               (3)

Its value on any spatial twisted exact form is zero, by the preceding
Stein identity. Thus it descends to the FULL localized curve cohomology.

## 4. The covariant chain rule proves horizontality for ALL iterates

Write hats for evaluation at z(t,S),U(t,S). For an arbitrary source
form omega=f dz put z_t=partial_t z and z_S=partial_S z. Evaluation
has the exact derivation rules

    partial_t fhat=(delta_t f)^hat+(delta_z f)^hat z_t,
    partial_S fhat=(delta_z f)^hat z_S.

Furthermore phi(t,z(t,S))=-S^2/2+Psi(t,S). The chain rule gives

    partial_t(fhat z_S)+Psi_t fhat z_S
      -[(delta_t f+phi_t f)^hat]z_S
       =[partial_S-S+Psi_S](fhat z_t).                            (4)

Every term belongs to B; the identities for U derivatives above
include f involving arbitrary positive or negative U powers. Multiply
by exp(Psi) and apply G. The right side vanishes as
G[(partial_S-S)(exp(Psi)fhat z_t)]. Hence

    partial_t Lambda(omega)=Lambda(nabla_t omega).                 (5)

The proof applies to every allowed form, and nabla_t preserves its
source algebra. Induction therefore proves (5) for EVERY finite
parameter-covariant iterate, rather than only a finite list or just
the first scalar value.

Finally z_S=2HS/J, so eta pulls back to
2/[sqrt(d)FJ]dS. Equations (2),(3) yield

    Lambda(eta)=f_P(t).                                           (6)

This scalar equality is exact. It does NOT identify original W-marked
periods with the alternate localized values. For example Jker has
localized class zero but original period N!=0. No such identification
is used in (4)-(6).

The formal-Hom uniqueness theorem for the ORIGINAL source is not
contradicted. The alternate functional takes general original markers
to C((t)), possibly with marker-dependent poles, rather than always
to C[[t]]. In the exact model ell=z(1+z), h=c=a=u=k=0, its W-period
is 1/(2t), whereas the original W-period is zero. Its unit period and
all unit parameter derivatives still agree exactly with the actual
ones. This is precisely why the cyclic argument below is legitimate
without an all-original-marker transport claim.

## 5. The genuine rational cyclic quotient under a finite-image hypothesis

Suppose f_P=exp(psi)p with psi in tC[t], p nonzero polynomial. Its
constant coefficient makes p(0)=1. Shift the curve phase to phi-psi,
and define Lambda'=exp(-psi)Lambda. This is horizontal for the shifted
parameter connection by (5). The form eta/p is in the same K-module
and has Lambda'(eta/p)=1.

Let Cseed be its differential K-span. For every n>=0,

    Lambda'(nabla_t^n(eta/p))=partial_t^n(1).

Thus Lambda' on Cseed takes values in K, not merely C((t)), and
defines a nonzero horizontal K-linear quotient Cseed->K. Every source
relation is respected because Lambda' was constructed on full curve
cohomology. This is the exact cyclic quotient required by the retained
comparison theorem. It does not require the original full source to
survive localization or its functional to descend there.

Assume now deg ell>=2 and every root of gcd(ell,c^2) has multiplicity
at most two. The retained cover-volume proof gives positive order
for eta on EVERY component of EVERY sufficiently refined marked model;
its finite/infinity inequalities concern q and eta alone. Hence they
remain valid after a spatially constant polynomial phase shift psi.
The generic horizontal poles of eta are logarithmic or holomorphic.
Dividing by polynomial p increases its infinity order by deg p>=0.
The phase is spatially nonconstant: its leading t^-1 coefficient is
-z/(4L), nonconstant since L(0)!=0. Scalar subtraction cannot remove it.

The retained full compact-piece, annulus, descending-tail and period
comparison therefore excludes exactly the quotient just constructed.
This contradicts f_P=exp(psi)p. In particular no finite actual image
exists in the specified simple-origin volume range. The comparison
dependencies and pending external-review status are exactly those
already declared in the canonical double-cover theorem.

The construction uses no other-root source-cleanliness assumptions.
If those are retained for the original Gysin discussion they cause no
problem, but they are unnecessary for this independent scalar bridge.
Over characteristic-zero coefficient fields the exponential-polynomial
identity and the finitely many polynomial coefficients descend to a
finitely generated field admitting an embedding in C; nonzero leading
coefficients and the gcd condition are preserved. Thus the exclusion
descends as well.

## Proof dependencies and exact controls

The original homotopy and source context are retained in the
[unclean-origin section](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#the-actual-kernel-at-an-unclean-simple-origin).
The geometric input is exactly the retained
[finite-centre volume proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#positive-volume-order-over-finite-centres)
and [complete infinity/period comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison).
The [focused checker](verify_gaussian_simple_origin_alternate_scalar.py) and
[receipt](gaussian_simple_origin_alternate_scalar_verification.json) replay
exact rational-model covariant identities, actual scalar coefficients and
the original-marker discrepancy. These bounded controls do not prove the
imported comparison theorem or unrestricted all-order scope. The written
argument supplies the latter. Independent replay, formal verification and
external mathematical review are not claimed.
