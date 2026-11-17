# Actual unit divergence and finite images for a one-root mixed coefficient

Internal written proof, 17 November 2026. This is a negative finite-actual-
image result obtained from the ORIGINAL Wick unit, followed by a complete
one-distinct-root finite-image criterion. It does not use an alternate
localized trace or infer an infinite claim from moment prefixes. Internal
audits checked the original source, selected inverse sheet, fixed critical
value and uniform contour/convolution argument. Independent mathematical
review, formal verification and external review are not claimed; priority
is unassessed.

The source identity below is derived directly, and is also the one-base
specialization of the [exact multibase Wick proof](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#1-a-self-contained-multibase-wick-contraction).
The algebraic finite-image normalization and the distinction between V
and M are retained in the [original-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
The core divergence proof in Sections 1-4 uses only the displayed Wick
identities and elementary complex implicit inversion, Cauchy contours and
coefficient bounds; it does not import a stationary-phase comparison.
The complete family criterion in Section 7 additionally uses the
[W-affine H5 finite-image criterion](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#7-complete-finite-image-criterion-in-the-nonconstant-c-h5-family),
its [unshared-root estimates](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#4-unshared-roots-of-arbitrary-multiplicity),
and the [exact mixed zero-root theorem](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#9-a-separate-mixed-source-extension-with-zero-root),
with the dependencies and assurance boundaries stated there.

Let m>=2, alpha beta!=0, c(z)=(alpha+beta z)^m, and

    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

with arbitrary complex polynomial h,u,k and constant a. Then the actual
unit series Phi_P[1] is divergent at t=0. In particular the whole actual
image M cannot be finite over C[t], and P is not pure or safe.
In fact the two ORIGINAL periods Phi_P[1] and Phi_P[P] are linearly
independent over C(t), so the actual generic marked-image rank is at
least two.

There is no common-root condition on h,c. Thus this includes genuinely
non-H5 branches such as m=4 and ord_(z=1)h=2. The case m=1 has a
different unperturbed singularity (at spatial infinity) and is not
proved by this argument; its no-finite-image conclusion already follows
from the separate squarefree/H5 theorem.

Equivalently, over C the result covers every c=kappa(z-alpha0)^m with
kappa!=0, alpha0!=0 and m>=2: choose beta^m=kappa and alpha=-beta alpha0.
No m-th root in the original coefficient field is assumed for the
characteristic-zero consequences in Section 6.

## 1. Exact original source kernel, including all markers

Set d=1-2at and choose sqrt(d) with constant coefficient one. Exact
Gaussian Y completion gives

    Phi_P[1]=d^(-1/2) E_ZW[D(t,Z)exp(A(t,Z)W^2+B(t,Z)W)],
    A=t^2 c^2/(2d),
    B=t h+t^2 c u/d,
    D=exp(tk+t^2u^2/(2d)).

Using exp(AW^2)=E_S exp(t c S W/sqrt(d)), then circular Lagrange
contraction, let z=z(t,X) be the analytic germ satisfying

    z=B(t,z)+c(z)X/sqrt(d),   z(0,0)=0,
    J=1-B_z(t,z)-c'(z)X/sqrt(d).

Define

    K(t,X)=D(t,z(t,X))/[sqrt(d)J(t,X)].

The exact actual-unit identity is

    Phi_P[1](t)=E_S[K(t,tS)].                        (D1)

The expectation is formal: expand K in X, substitute X=tS, then
contract Gaussian powers. Every t coefficient involves finitely many
polynomial S coefficients. The Gaussian square completion and circular
Lagrange identity concern the ORIGINAL Wick functional throughout.

For completeness this kernel retains all original markers. Introduce
sources vW+xiY. Replace B,D by

    B_xi=t h+t c(tu+xi)/d,
    D_xi=exp(tk+(tu+xi)^2/(2d)),

and let z=v+B_xi(t,z)+t c(z)S/sqrt(d). Then

    Phi_P[Z^r W^s Y^j]
      =partial_v^s partial_xi^j
        E_S[z^r D_xi/(sqrt(d)(1-(B_xi)_z-tc'S/sqrt(d)))] at v=xi=0. (D2)

Every fixed source derivative is coefficientwise finite. Equation (D2)
is an exact original-marker formula, not a claim that its full image
has already been classified. Only its unit is needed for divergence.

## 2. The unperturbed inverse has one selected simple critical point

At t=0, the selected inverse satisfies z=c(z)X. Put

    gamma=beta alpha^(m-1),
    rho=(m-1)^(m-1)/m^m,
    zstar=alpha/[(m-1)beta],
    Xstar=rho/gamma.

The inverse function G_0(z)=z/c(z) has exactly one finite nonzero
critical point zstar, with G_0''(zstar)!=0. The other ramification
values are at zero or infinity and do not affect the selected germ
z(0,X) at X=0: this is the unique finite root there, so local
monodromy around zero fixes it. In the coordinate v=1+beta z/alpha,
its equation is

    v=1+gamma X v^m.

Along the scaled real interval gamma X from zero to rho, the selected
root runs from 1 to m/(m-1). Hence this selected root meets zstar at
Xstar; the simple critical singularity belongs to the actual branch.

Formal Lagrange inversion gives the exact coefficients

    K(0,X)=sum_(n>=0) b_n gamma^n X^n,
    b_n=binom(mn,n).                                (D3)

Indeed b_n gamma^n=[z^n]c(z)^n. In particular no coefficient vanishes.
On a disk cut along the outward radial slit from Xstar and with zero
temporarily removed, the fundamental group is generated by a loop
around zero. The preceding monodromy observation therefore makes the
selected branch single-valued, and its local implicit-function germ
fills zero. It remains bounded on the closure of a slightly smaller
slit disk, with the two banks kept as distinct boundary copies: use
polynomial-root bounds away from zero, the local finite germ at zero,
and the finite square-root limit at Xstar. A nonzero finite X cannot
send its root to infinity.

## 3. Fixing the perturbed critical value and its selected sheet

For small t the inverse function is

    G_t(z)=sqrt(d)[z-B(t,z)]/c(z).

The simple critical point and its nonzero value persist analytically:
zstar(t), Xstar(t), with the displayed values at zero. Normalize

    T=t Xstar(0)/Xstar(t).

This analytic parameter change has derivative one and inverse t=t(T).
Put eta(t)=Xstar(t)/Xstar(0) and

    Ktilde(T,X)=K(t(T),eta(t(T))X).

The distinguished critical value is now fixed at Xstar(0), and (D1)
becomes exactly

    Phi_P[1](t(T))=E_S[Ktilde(T,TS)].                (D4)

There is no S-dependent parameter change.

Here is the uniform selected-branch justification. Work with the
RESCALED equation z=B(t,z)+c(z)eta(t)X/sqrt(d), t=t(T), whose selected
critical value is fixed at Xstar(0). Remove a small disk
around Xstar from a slightly smaller closed slit disk. The unperturbed
selected root graph on this compact set is bounded, avoids the root
of c, and its defining z derivative is nonzero. Uniform implicit
inversion in finitely many isolated root neighborhoods gives a joint
analytic branch for all small complex T; the local branches glue by
uniqueness near the unperturbed graph. Extra roots or critical points
created by arbitrary h,u,k at unbounded z cannot enter these bounded
isolated neighborhoods, even if their critical values approach zero.
Inside the removed disk the simple quadratic critical point gives a
joint analytic square-root inverse after the critical-value normalization.

Thus Ktilde has, uniformly for small T, a fixed-slit continuation and
local form

    P(T,X)+Q(T,X)(1-X/Xstar(0))^(-1/2),              (D5)

with P,Q jointly analytic and bounded near the critical point. The
amplitude D/sqrt(d) is an analytic unit on the bounded graph, and J
has a simple zero in the square-root coordinate, so no stronger
singularity occurs.

## 4. Uniform coefficient bounds and Gaussian divergence

Choose the outer radius just beyond |Xstar|, so the entire short slit
lies in the fixed quadratic-critical-point disk. Deform the X coefficient
contour to that outer contour and the two slit banks; the tiny circle
around the critical point contributes O(epsilon^(1/2)) and vanishes.
The jump bound in (D5) gives

    |[X^n]Ktilde(T,X)|
      <=M |Xstar|^(-n) integral_0^delta
                      s^(-1/2)(1+s)^(-n-1)ds+O(R^(-n))
      <=M' |Xstar|^(-n)(n+1)^(-1/2).

This is uniform on a fixed small closed T disk. The exact factorial
ratio for b_n=binom(mn,n) gives

    b_n/b_(n-1)=rho^(-1)(1-1/(2n)+O(n^-2)).

Therefore b_n is comparable to rho^(-n)(n+1)^(-1/2): after multiplying
by rho^n sqrt(n), its successive ratios are 1+O(n^-2), whose positive
product converges. Cauchy's estimate in T consequently gives

    |[T^r X^n]Ktilde(T,X)|<=M R0^(-r)|gamma|^n b_n. (D6)

Let

    L_q=(2q-1)!! |gamma|^(2q) b_(2q),   L_0=1.

Its successive ratio satisfies

    (L_q/L_(q-1))/q ->2 |gamma|^2 rho^(-2)>0.        (D7)

In the coefficient of T^(2q) in (D4), odd powers of S vanish. The
term with X power 2q and T power zero contributes exactly

    (2q-1)!! gamma^(2q)b_(2q).

Every other surviving term has X power 2j with j<q. By (D6), their
total absolute contribution is at most

    M sum_(j<q) L_j (R0^-2)^(q-j)=O(L_q/q).         (D8)

For the last estimate split j>=q/2 and j<q/2. On the first range
factorial ratios give L_j/L_q<=(C/q)^(q-j), so the sum is geometric
of order 1/q. On the second, eventual monotonicity and domination of
the finite initial segment bound L_j by L_floor(q/2), whose ratio
to L_q is at most (C/q)^floor(q/2). This defeats the fixed exponential
majorant even after summing q terms.

Hence the coefficient in (D4), divided by its displayed main term,
tends to one. Its magnitude grows faster than every exponential by
(D7). The actual unit is therefore not analytic at t=0, because an
invertible analytic parameter change preserves analyticity.

If M were finite, the retained exact source/PID theorem would give
Phi_P[1]=exp(psi)p with polynomial psi,p; this is analytic. The
contradiction proves the finite-image obstruction.

There is also an explicit original-marker rank obstruction. Write
f=Phi_P[1], so Phi_P[P]=f'. If f'/f were rational, its formal expansion
would be regular at zero because f(0)=1. Its Taylor series is therefore
convergent, and the first-order equation f'=(f'/f)f with initial value
one has a convergent analytic solution. Uniqueness of its formal solution
contradicts the proved divergence. Thus f'/f is nonrational, and f,f'
are C(t)-linearly independent. This identifies a fixed original marker
P witnessing rank at least two, rather than merely failing finiteness
of an unspecified lattice.

## 5. The proposed m=4, ord h=2 prototype

For h=(z-1)^2, c=(z-1)^4, a=k=0, u=-1/4, ordinary Gaussian shifting
separates the exact scalar factor exp(t^2/32). The original normalized
unit is

    exp(-t^2/32)Phi_P[1]
       =sum_(n>=0)t^n[z^n]E_Y[(h+c(Y-t/4))^n].      (D9)

Its first coefficients are exactly

    1-2t+35t^2-(775/2)t^3+... .

These are algebraic controls only. The unrestricted proof is Sections
1-4, with zstar=-1/3 and Xstar=-27/256. It proves the actual unit and
its normalized version divergent, so this example is neither a finite
image nor a safety-to-purity counterexample. The quartic zero-volume
chart at parameter infinity is consistent with this separate original
source obstruction at t=0.

## 6. General scale and characteristic-zero descent

Over C, let c=kappa(z-alpha0)^m with kappa!=0, alpha0!=0 and m>=2.
Choose beta^m=kappa and put alpha=-beta alpha0. The entire preceding
proof applies to c=(alpha+beta z)^m. The choice of m-th root affects
only that coordinate description, not the actual source or conclusion.

Now let the phase have coefficients in ANY characteristic-zero field K,
with the same displayed factorization of c over K. Its actual original
periods are defined by the same rational Wick contractions. If the full
actual image were finite over K[t], the general source/PID theorem
would give f=exp(psi)p with polynomial psi,p. Adjoin to Q the finitely
many phase, psi and p coefficients and embed that finitely generated
field in C. The exact formal unit identity is preserved. Its embedded
right side is analytic at zero, contradicting the proved complex
divergence. Thus no such finite image exists over K.

Likewise suppose f'/f were rational over K(t). Include the finitely
many coefficients of this rational function in the embedded field.
Its regularity at zero is forced by f(0)=1, and the embedded first-order
equation would again make f analytic. Therefore

    f'/f is not in K(t),
    rank_(K(t)) span{Phi_P[1],Phi_P[P]}=2.           (D10)

This is an explicit ORIGINAL-marker obstruction: the fixed marker P
itself witnesses actual generic rank at least two. No localization or
replacement of the marked map is involved.

## 7. Complete finite-image classification with one distinct root

Let K be a characteristic-zero field and

    c(z)=kappa(z-alpha0)^m,
    kappa!=0, m>=1,
    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),

where a is constant and h,u,k are arbitrary polynomials. Then the WHOLE
actual image satisfies

    M is finite over K[t]
      iff alpha0=0, m>=2, a=0, h(0)=h'(0)=0.        (D11)

There is no restriction on u or k. In every finite case, writing
u0=u(0), k0=k(0), the exact primitive image is

    psi=k0t+u0^2t^2/2,
    f=exp(psi), V=M=exp(psi)K[t].                   (D12)

Here V is the literal K-linear image of all fixed original polynomial
markers. The phase is necessarily nonconstant because c!=0, so there
is no constant-phase exception. Safety is equivalent to (D11) together
with u0=k0=0; every such safe phase is pure.

For alpha0!=0 and m>=2, Sections 1-6 prove necessity by divergence.
For alpha0!=0 and m=1, c is squarefree and the imported H5 criterion
excludes finite images. For alpha0=0 with h(0)!=0, the sole c-root is
unshared, and the imported unshared-root argument again excludes them.
For alpha0=h(0)=0, the exact mixed zero-root theorem says precisely
a=h'(0)=c'(0)=0. Since c'(0)!=0 for m=1 and is zero for m>=2, this
gives exactly (D11). That theorem also proves sufficiency for ALL
original markers, the unit formula in (D12), and its safety condition.
The equality V=M follows from the general finite fixed-marker theorem.

The imported H5 rejections are stated over C, and extend to the field
form of (D11) by finite-coefficient descent. To see that the WHOLE
finite-image hypothesis is preserved, write M=exp(psi)K[t] and embed
the field generated by the phase and psi coefficients. Every normalized
original MONOMIAL period has coefficients in that same field, because
they are obtained from Wick coefficients and exp(-psi). Its tail
vanishes over K and hence over the smaller field. After embedding,
every such normalized period is still polynomial, so the complex
image is a submodule of exp(psi)C[t] and is finite. This contradicts
the imported complex obstruction. No infinite collection of new
marker coefficients needs to be adjoined. The zero-root sufficiency
is already algebraic over every characteristic-zero field.

## 8. Verification and exact limits

The [checker](verify_gaussian_w_affine_single_root_divergence.py) and
[receipt](gaussian_w_affine_single_root_divergence_verification.json)
retain finite exact controls of Gaussian completion, critical geometry,
baseline Lagrange coefficients, independent original Wick coefficients
for the prototype, and its factorial ratio. Replay from the repository
root with

    .venv/bin/python research/moments/vanishing/verify_gaussian_w_affine_single_root_divergence.py --output research/moments/vanishing/gaussian_w_affine_single_root_divergence_verification.json

Those bounded controls do NOT certify the unrestricted selected-sheet
continuation, contour estimate, or all-marker classification. The
written proofs and explicitly linked imports supply those conclusions.
The actual generic image in the divergent cases is proved to have rank
at least two; its entire higher-rank image and original kernel are not
classified here. General field safety implies purity remains open.
