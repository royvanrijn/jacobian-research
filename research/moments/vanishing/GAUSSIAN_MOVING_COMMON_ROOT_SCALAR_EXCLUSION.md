# Moving common-root scalar bridge and residual-source exclusion

Internal written proof. Consider

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

with constant a and polynomial other coefficients. Put K=C(t), and

    d=1-2at, Psi0=tk+t^2u^2/(2d), D=exp(Psi0),
    A=t ell+t^2c^2/(2d), B=t h+t^2cu/d, p=z-B,
    q=2dA/t^2=c^2+2dell/t,
    phi=Psi0-p^2/(4A), eta=dz/(tU), U^2=q.

The source of the bridge is the ORIGINAL rational curve algebra

    Acurve=K[z,U,U^-1]/(U^2-q).

Its phase is regular since A is a nonzero scalar multiple of U^2.

## 1. The precise moving-root hypothesis

Assume A and p have a common algebraic germ beta(t) in tC[[t]],
such that A_z(t,beta) is nonzero. Thus beta is a simple generic
A-root, although its specialization at t=0 can be a multiple raw
ell-root. Define the translated quantities by z=beta+x:

    A(t,beta+x)=x A1(t,x),
    p(t,beta+x)=x F(t,x),
    H(t,x)=2A1/F^2.

Because A is t-divisible, A1 belongs to tC[x][[t]]. Because p=z+O(t)
and beta(0)=0, F=1+O(t) in C[x][[t]]. Let

    A1(t,0)=a0 t^nu+O(t^(nu+1)), a0!=0, nu>=1.

The case raw ord_0 ell=1 has nu=1 and a0=ell'(0). The construction
below also allows nu>1. These statements use the selected embedding
of the finite algebraic field K(beta) in C((t)); no rationality of
beta over K is assumed.

## 2. Translation preserves the ACTUAL scalar residue

Exact ordinary-Y elimination followed by original circular Wick
contraction gives

    f_P=d^(-1/2) sum_(j>=0)(2j)!/j!
                   Res_z D A^j/p^(2j+1).                           (1)

Every denominator is expanded t-adically. Each fixed t coefficient
is a finite Laurent polynomial in z, and only finitely many j occur.
Translation z=x+beta(t) preserves its formal residue coefficientwise.
Indeed for each integer k the expansion of (x+beta)^k has x-residue
one if k=-1 and zero otherwise. For negative k its further terms have
strictly smaller x exponents; for nonnegative k it is polynomial.
Since beta is t-divisible, this substitution is t-adically locally
finite. No contour or branch outside the original selected residue
has been added.

Consequently (1), after the exact cancellation of x powers, becomes

    f_P=d^(-1/2) sum_(j>=0)(2j-1)!!
                   Res_x (D(t,beta+x)/F) H^j/x^(j+1).              (2)

The implicit equation x=H(t,x)S^2 has a unique t-adic solution with
polynomial S coefficients. Since H(t,0)=2a0t^nu+..., while H_x is
t-divisible,

    x(t,S)=2a0 t^nu S^2+O(t^(nu+1)),
    J=1-H_x(t,x(t,S))S^2=1+O(t).

The same near-identity formal residue change as in the fixed-origin
proof gives the exact actual scalar formula

    f_P=d^(-1/2) E_S[D(t,beta+x)/(FJ)].                              (3)

The expectation in (3) is ordinary Gaussian contraction of polynomial
S coefficients, so this is an identity of the original Wick series.

## 3. The full original curve domain still maps to a fixed Laurent ring

Set z(t,S)=beta(t)+x(t,S) and

    U(t,S)=sqrt(d) F(t,x)H(t,x)S/t.

Its square is q(t,z(t,S)), and

    U(t,S)=2a0 t^(nu-1)S(1+O(t)).                                  (4)

Thus U^-1 is a Laurent-t monomial times S^-1 times a formal unit.
The ENTIRE original curve algebra, every regular one-form f dz, and
every finite parameter-covariant iterate evaluate in the fixed rings

    C[S,S^-1]((t)), C[S,S^-1]((t))dS.

Here the derivative lift is the intrinsic original one:
delta_t U=q_t/(2U), delta_z U=q_z/(2U). Both preserve Acurve.
No translated curve with formal coefficients is used as a global
comparison model. Beta is used only to construct this formal map
from the original K-curve at t=0.

Writing Psi(t,S)=Psi0(t,z(t,S)), one has

    phi(t,z(t,S))=-S^2/2+Psi(t,S),
    eta -> 2/[sqrt(d)FJ]dS,
    Psi in tC[S^2][[t]].

Use the FIXED Laurent Gaussian G with odd moments zero and
G(S^(2m))=(2m-1)!! for m>=0,
G(S^(-2n))=(-1)^n/(2n-1)!! for n>=1. It kills partial_S-S on the
entire Laurent ring and commutes with partial_t. Define

    Lambda(b dS)=(1/2)G[exp(Psi)b].

This kills all spatial twisted exact forms and gives Lambda(eta)=f_P
by (3). For every original form omega=f dz its moving-coordinate
chain-rule defect is EXACTLY

    partial_t(fhat z_S)+Psi_t fhat z_S
       -[(delta_t f+phi_t f)^hat]z_S
        =[partial_S-S+Psi_S](fhat z_t),                            (5)

where z_t=beta'(t)+x_t. The beta' term is retained. Every term belongs
to the same fixed Laurent ring. Hence Lambda is horizontal on full
curve cohomology and therefore on every parameter-covariant iterate.

The original W/Gysin marker map need not agree with Lambda and may
have different regularity at t=0. This does not affect the exact
scalar equality or its covariant derivatives.

## 4. Finite-image exclusion under the retained volume hypotheses

Suppose deg ell>=2 and every root of gcd(ell,c^2) has multiplicity
at most two. Assume also that phi is spatially nonconstant, as verified
for the residual branch in Section 5. If the actual vacuum were exp(psi)p, with psi in tC[t]
and nonzero polynomial p, gauge the ORIGINAL curve phase to phi-psi
and put Lambda'=exp(-psi)Lambda. On the differential K-span of eta/p,

    Lambda'(nabla_t^n(eta/p))=partial_t^n(1), n>=0.

This defines a nonzero K-valued horizontal cyclic quotient. It is
not inferred from an original-source localization isomorphism.

The cover and normalized volume are still the original q and eta.
The retained positive-volume inequalities apply to every algebraic
component independently of the scalar phase shift; multiplying eta
by 1/p increases its infinity order by deg p. All horizontal volume
poles remain logarithmic or holomorphic. The retained complete
comparison therefore excludes this cyclic unit quotient.

The phase is spatially nonconstant in the residual case below: its
leading t^-1 coefficient is -z^2/(4ell). If this were constant, ell
would be proportional to z^2, which cannot have residual e with a
simple origin and C(0)!=0, since ord_0g is even. Thus the ordinary-
parameter exponential comparison applies. A finite actual image
always supplies such an exponential-polynomial vacuum, so it is
excluded. The comparison dependencies and review limits are exactly
the retained double-cover volume theorem's, without a new geometric
or stationary-phase import.

## 5. Algebraic verification on the residual Rmov=0 branch

Let

    g=gcd(ell,c^2), e=ell/g, C=c^2/g,
    f=4ae-C, J0=-2e+tf=-2de-tC,
    Rmov=C[zf-2eh]+4e^2cu.

Assume e=kappa z with kappa!=0, C(0)!=0, and Rmov is identically
zero. This is the entire nonconstant-e residual branch: reducing
Rmov modulo e gives -zC^2, and gcd(e,C)=1 forces e to divide z.
Thus a nonconstant e must be a nonzero scalar multiple of z.
Write e1=kappa and C(0)=C0. The z coefficient of Rmov gives

    h(0)=-C0/(2e1)!=0.

The polynomial 2de+tC has a unique algebraic root beta in tC[[t]],
with beta=t h(0)+O(t^2), by the formal implicit-function theorem.
Since g(beta) is nonzero in C((t)), beta is a simple generic root
of A=t g(2de+tC)/(2d). At J0(beta)=0 the exact source relation is

    t^2 Rmov(beta)=-4d e(beta)^2 p(beta).

Here e(beta) is nonzero; hence p(beta)=0. This proves the common-root
hypothesis in Section 1, to ALL orders, without formally inverting e
at zero.

If gamma=ord_0g and g_gamma is its leading coefficient, then

    A1(t,0)=e1 g_gamma h(0)^gamma t^(gamma+1)+O(t^(gamma+2)).

Thus nu=gamma+1, with nonzero leading coefficient. Since c^2=gC
and C0!=0, gamma is even. The bridge itself permits arbitrary gamma;
the stated volume theorem restricts gamma to zero or two when it
is used at the origin. In particular it covers both

    ord_0ell=1, c(0)!=0;
    ord_0ell=3, ord_0c=1.

In the first case ell1=ell'(0)!=0 and the relation is precisely

    2ell1 h(0)+c(0)^2=0,

as proposed. In the second case U has leading nonzero multiple of
t^2S, which remains invertible in the fixed Laurent ring.

All algebraic identities are defined over the coefficient field;
the formal beta coefficients are determined there recursively.
Finite-coefficient embedding into C gives characteristic-zero descent
of the finite-image exclusion. The residual constant-e branch, or
branches whose relevant roots escape t=0, are not covered here.

## Examples and verification scope

The degree-one model ell=z,c=1,h=-1/2,a=u=k=0 has beta=-t/2 and
is pure and safe by strict angular support. It satisfies the bridge,
but fails the deg ell>=2 volume premise, so is not excluded.

A nonrational-beta test is

    ell=z(1+z^2), c=1+z^2, h=-(1+z^2)/2, a=u=k=0,
    beta=(sqrt(1-t^2)-1)/t.

Here g=C=1+z^2, e=z and all volume multiplicities are one. The
high-valuation control ell=z^3,c=z,h=-1/2 has gamma=2, beta=-t/2.
Finite controls should compare original Wick moments with (3), verify
translation for Laurent monomials, and check (5) with beta' retained.
The unrestricted conclusions above use the exact formal residue,
domain and covariant identities, not any bounded moment prefix.

## Dependencies and verification

The fixed-coordinate covariant argument and its exact original-marker
scope are given in the [simple-origin scalar bridge](GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md).
The geometric input is the retained [finite-centre volume proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#positive-volume-order-over-finite-centres)
and [complete infinity comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison).
Finite actual images supply exponential-polynomial vacua by the
[general original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
The [focused checker](verify_gaussian_moving_common_root_scalar.py) and
[receipt](gaussian_moving_common_root_scalar_verification.json) replay
bounded exact original Wick, translation, algebraic-beta and moving
covariant identities. The unrestricted scalar bridge and exclusion
are the written argument above. Imported comparison dependencies retain
their existing review limits. Independent replay, formal verification
and external mathematical review are not claimed.
