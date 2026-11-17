# Finite original images with an arbitrary W-affine polynomial fibre

Let K have characteristic zero. Use the normalized Gaussian Wick
functional with E[ZW]=E[Y^2]=1 and all other basic pairings zero. Write
Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]], and M_P=K[t]V_P. Safety
means that EVERY fixed original polynomial marker has polynomial period.

This is an internal complete written extension of the retained finite-image
proof. Its curve comparison, full stationary-ball theorem, and general
finite-image theorem retain their existing assurance. Bounded exact checks
are supplied separately; they do not prove the comparison imports.

## Criterion

For c NONCONSTANT and arbitrary h,c,G over a characteristic-zero field,

    P=W(h(Z)+c(Z)Y)+G(Z,Y)

has finite full actual image if and only if

    h,c belong to Z^2 K[Z],
    G(0,Y)=uY+k is affine.                          (C1)

In every finite case the vacuum is primitive and, since P is nonconstant,

    psi=kt+u^2t^2/2,
    Phi_P[1]=exp(psi), V_P=M_P=exp(psi)K[t].         (C2)

Safety is equivalent to (C1) plus u=k=0. In particular every safe phase
in this entire architecture is pure. No converse from a scalar unit alone,
or general purity-to-safety classification for arbitrary G, is asserted.

## 1. The logarithmic original marker outside the common origin

Let

    P=W(h(Z)+c(Z)Y)+G(Z,Y),

where G is ANY complex polynomial, c is nonconstant, and
(h(0),c(0))!=(0,0). Then M is not finite. In particular every safe member
of this enlarged architecture must lie in the common-zero-origin branch.

Put s=(z/t-h)/c, phi=-s^2/2+tG(z,s), eta=dz/(tc). The original source is

    zeta=t[h(zeta)+c(zeta)Y],
    J=1-t[h'(zeta)+c'(zeta)Y],
    eta -> dY/J,
    phi(t,zeta)=-Y^2/2+tG(zeta,Y).

The exact circular Lagrange identity gives

    Phi[A(Z,Y)]=E_Y[A(zeta,Y)exp(tG(zeta,Y))/J]

for every original polynomial A independent of W. The spatial and parameter
Stein identities give the full rational-curve connection exactly as in the
retained W-affine source proof; the degree of G in Y is irrelevant.

If c(0)!=0, c(zeta) is a t-adic unit. If c(0)=0,h(0)!=0, then
zeta=t h(0)+O(t^2) and c(zeta) has leading coefficient
c_m h(0)^m t^m independent of Y. Thus all required c-localized amplitudes
have Laurent-t expansions with polynomial Y coefficients, and their
original Gaussian evaluation remains well-defined and horizontal.

Take monic rho=gcd(c,c') and c_red=c/rho. Then

    rho eta=r dz/c_red, r=1/t,

with fixed logarithmic form omega=dz/c_red. The reduced rational phase
has fixed spatial poles only at roots(c) and infinity, regardless of G.
If c(0)!=0, Phi[rho] at t=0 equals rho(0)!=0. Otherwise write
rho=rho_s z^s+..., s=ord0(c)-1; the exact first term is

    Phi[rho]=rho_s h(0)^s t^s+O(t^(s+1))!=0.

The fibre exponential and J^-1 start with one, so higher Y degree in G
cannot alter this term. Suppose M_P were finite. The general finite-image theorem gives
M_P=exp(psi)K[t]. Hence the nonzero ORIGINAL period Phi[rho] equals
exp(psi)p_rho for a nonzero polynomial p_rho. Replace the rational phase
by phi-psi and use seed (r/p_rho)dz/c_red. Its exact original value is one,
and its scalar infinity order is 1+deg p_rho>0. The retained
[fixed-puncture logarithmic seed bound](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound)
excludes that unit seed quotient on the complete marked curve. Its
hypotheses allow arbitrary rational phase with these fixed poles,
including the spatially constant shift psi. It keeps finite components,
annulus gluing, regular punctures and descending tails, as in Section14
of the [constant-quadratic-fibre proof](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#14-a-fixed-original-logarithmic-marker-removes-every-root-order-restriction).
This proves the claim. All-field descent adjoins the finite phase, gcd, generator and
marked-polynomial coefficients, then preserves the original formal identity
and this nonzero leading term under embedding in C.

## 2. The common origin

Assume now h(0)=c(0)=0.

### Exact original cancellation into the one-variable source

Put eta=h'(0), gamma=c'(0) and

    J0(t,Y)=1-t(eta+gamma Y), G0(Y)=G(0,Y).

The unique formal circular root zeta=t[h(zeta)+c(zeta)Y] is now
IDENTICALLY zero. The exact original Lagrange identity consequently gives

    Phi_P[A(Y)]=E_Y[A(Y)exp(tG0(Y))/J0(t,Y)].

The expression on the right is a valid coefficientwise identity:
J0 has constant coefficient one, and every t coefficient is polynomial
in Y. For every original polynomial A(Y), extend Phi K[t]-linearly to
the finite time-dependent amplitude J0 A. Then

    Phi_P[J0(t,Y)A(Y)]=E_Y[A(Y)exp(tG0(Y))].          (C3)

Thus the WHOLE one-variable Gaussian image M_G0 is contained in M_P.
Explicitly (C3) is a combination of the fixed original markers A and YA,
with coefficients 1-t eta and -t gamma. It does not assume that J0 is
itself a fixed spatial marker.

If M_P is finite, its K[t]-submodule M_G0 is finite by Noetherianity.
The following one-variable fact therefore applies.

### One ordinary Gaussian variable: finite image means affine phase

For any g in K[Y], the full ordinary Gaussian image is finite over K[t]
if and only if deg g<=1. If g=uY+k, the exact Gaussian formula and
polynomial marked Wick formulas give

    M_g=exp(kt+u^2t^2/2)K[t].

For necessity first work over C. Suppose deg g=D>=2 and M_g is finite.
The general actual-image theorem gives

    E_Y exp(tg(Y))=p(t)exp(psi(t)),
    p in C[t], p(0)=1, psi in tC[t].

Use the polynomial spatial phase and normalized volume

    phi=-Y^2/2+t g(Y)-psi(t), h(t)dY=dY/p(t).         (C4)

The original Gaussian Stein identities define a horizontal functional on
the absolute polynomial exponential de Rham system with this phase and
volume. Its unit seed has value one. Equivalently the scalar-normalized
functional is exp(-psi)Phi_g/p, with the connection term h'/h retained.

The exact imported statement is Formal stationary balls, theorem
thm:stationaryballs in

    unpublished/gaussian-boundary/polynomial_stationary_balls.tex,

also retained as the Leading-coefficient obstruction in the
constant-mixed-source section of GAUSSIAN_AFFINE_ISOTROPIC_CHART.md.
For a polynomial phase of degree D>=2 with leading coefficient a_D and
volume scalar h, it excludes a rational horizontal unit-seed quotient
whenever D v_infinity(h)>v_infinity(a_D). Its full formal frame includes
colliding, multiple and escaping critical points.

In (C4), the actual spatial degree is D and

    a_D=t g_D                   if D>2,
    a_2=t g_2-1/2              if D=2.

Thus v_infinity(a_D)=-1, while v_infinity(h)=deg p>=0. The strict
criterion holds, contradicting the unit seed quotient. The spatially
constant shift psi does not affect a_D. This proves necessity for
EVERY D>=2, including the quadratic case. It uses the retained formal
leading-coefficient theorem, not a divergent-moment prefix.

The same conclusion holds over every characteristic-zero field: a
hypothetical finite M_g provides p,psi over that field. Embed the finitely
generated field of g,p,psi coefficients in C; its exact scalar identity
and nonzero highest coefficient survive. Sufficiency is algebraic.

### The remaining zero-angular core and all-marker sufficiency

By (C3), finite M_P forces G0=uY+k. Write S=ZW. Then

    P=Pzero+Pplus,
    Pzero=S(eta+gamma Y)+uY+k,

where EVERY monomial of Pplus has strictly positive Z-minus-W weight.
Indeed the remaining h,c terms have Z order at least two and one W,
and G-G0 is divisible by Z with no W. Arbitrary Y degree in G creates
no exception to this angular decomposition.

This is exactly the scalar-ordinary specialization of the retained
weak-one-sided quadratic-core theorem, Section10 of

    research/moments/vanishing/GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md.

Its necessity gives eta=gamma=0: the constant matrix A is zero, the
radial-square coefficient is zero, and the scalar null-chain condition
gamma^2=0 forces gamma=0 over a field. These are precisely the remaining
conditions in (C1). This application retains arbitrary positive-weight
Pplus, not just a quadratic G.

For completeness, sufficiency and primitive vacuum follow directly.
Under (C1), P=uY+k+Pplus with Pplus strictly positive in angular weight.
For an original marker Q of degree D, angular balance truncates the
positive-part expansion after at most D factors. Gaussian translation
in Y yields the exact formula

    Phi_P[Q]=exp(kt+u^2t^2/2)
       *sum_(n=0)^D t^n/n!
          E[Q(Z,W,Y+ut)Pplus(Z,W,Y+ut)^n].           (C5)

Every summand is polynomial in t. With Q=1 all n>=1 terms vanish by
positive angular weight, leaving precisely exp(psi). Thus all original
marked periods lie in exp(psi)K[t] and the vacuum itself generates that
module. The nonconstant-phase fixed-marker theorem gives V_P=M_P.
If u=k=0, (C5) is polynomial and the unit is one. Conversely safety makes
exp(psi) polynomial, forcing u=k=0 in a characteristic-zero field.

Combining this source argument with Section1 proves the unrestricted
classification (C1)-(C2). The fixed-logarithmic comparison and formal
stationary-ball imports retain their existing internal assurance status;
no additional external review or formal verification is claimed.

## 3. Marker cutoffs and finite perturbation jets

Let b be the maximum Y degree of Pplus in (C5), with b=0 when Pplus=0.
For a marker monomial Z^i W^j Y^alpha, only n<=j-i can contribute.
The translated amplitude has t degree at most alpha+bn; the nth
summand of (C5) therefore has degree at most alpha+(b+1)n. Consequently

    deg_t(exp(-psi)Phi_P[Q]) <= (b+1) deg Q.         (C6)

When u=0, there is no translation, so the sharper bound is deg Q.
In particular EVERY safe member of this architecture satisfies

    E[Q P^m]=0 whenever m>deg Q.                    (C7)

The general [mixed-cutoff jet theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md)
now gives, for any polynomial or formal coefficient perturbation of
spatial degree at most d at such a safe base point,

    I_infinity + m^(r+1) = I_((d+1)r) + m^(r+1).    (C8)

If an actual local nilpotence certificate m^nu subset I_N is available,
only the intervening moments through (d+1)(nu-1) require membership
checks. This is a sufficient bound and assumes the stated certificate;
it does not assert any particular nilpotence index or least ideal prefix.

The degree bound (C6) can be attained when u!=0. For example

    P=2Y+3+W(Z^2+Z^3Y)+ZY^4+Z^2Y^3

has exp(-psi)Phi_P[W]=16t^5+24t^3+5t and
exp(-psi)Phi_P[W^2]=256t^10+1792t^8+3456t^6+1840t^4+159t^2.
These exact formulae follow from the finite sum (C5). They describe a
cutoff, not surjectivity onto every intermediate time degree.

## Evidence and boundary

The proof uses the full actual image to normalize the nonzero logarithmic
marker. It makes no unrestricted claim that an exponential-polynomial
unit, or a pure unit, alone implies safety. The assumption that c is
nonconstant is essential to this argument: dz/c_red is then logarithmic
at infinity. Constant c can admit finite nonquadratic Gaussian shears
with a variable ordinary-fibre coefficient.

The retained [exact checker](verify_gaussian_extended_fibre_finite_images.py)
verifies original J0 cancellation for three markers through moment six,
including cubic G(0,Y), and compares six original marked periods with
(C5). It includes a control that fails if J0 is omitted. The receipt is
[gaussian_extended_fibre_finite_images_verification.json](gaussian_extended_fibre_finite_images_verification.json).
These are internal bounded exact controls. Independent replay, formal
verification and external review are not claimed. General Gaussian
safety-to-purity and finite-vacuum primitivity remain open.
