# Variable radial quadratic Gaussian fibres: finite images and primitive vacua

Let K be a characteristic-zero field. Use one Witt pair Z,W and an
independent standard Gaussian Y, with E[Z^iW^j]=delta_(i,j)i! and
E[Y^(2j)]=(2j-1)!!. Put S=ZW and consider

    P=P0+R, P0=g(S)+B(S)Y+a(S)Y^2,                 (V1)

where g,B,a are arbitrary polynomials. Every monomial of R has strictly
positive Z-minus-W weight; strictly negative support is included by
reversing the weights. R=0 is allowed, and no degree bound is imposed
on R, including its Y degree.

For ORIGINAL polynomial markers write Phi[Q]=E[Qexp(tP)],
V=Phi[K[Z,W,Y]], and M=K[t]V. The following are equivalent:

1. The actual unit is p(t)exp(psi(t)), with p polynomial, p(0)=1 and
   psi in tK[t].
2. M is finite over K[t].
3. a=0 and g,B are constant.

In this branch, writing g=c and B=b,

    p=1, psi=ct+b^2t^2/2, M=exp(psi)K[t].           (V2)

The original vacuum is primitive. For nonconstant P its literal image
V equals M by the retained original fixed-marker theorem. If P=c is
constant, V=Kexp(ct), while M=exp(ct)K[t].

Throughout (V1), purity and all-original-marker safety are equivalent.
Both require zero core: a=g=B=0, so P=R. The safe recompositions in a
separate indeterminate U form (U-c) precisely when a=B=0 and g=c
constant, and form the zero ideal otherwise.

The proof uses the original Gamma endpoint, a full anti-invariant
relative functional on a moving double cover, and complete marked-curve
comparison. For affine nonconstant a the volume has a double pole at
infinity; its irregular tail is treated explicitly, not by a logarithmic-
pole theorem. These are internal written proofs with bounded exact
controls. Independent replay, formal verification and external review
are not claimed. Arbitrary invariant cores of higher Y degree are not
classified by the main theorem.

## 1. Angular reduction and the constant-a import

The weight-zero core has exactly the original radial moments
E_S[S^m]=m!. Every term of the original vacuum containing R has positive
weight, so Phi_P[1]=Phi_P0[1]. In fact

    M_P=M_P0                                           (V3)

as actual embedded modules, by the
[angular-core theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md).
Necessity may therefore be proved on P0. The
[general finite-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
supplies the scalar premise pexp(psi) whenever M is finite; primitivity
of that unit is proved here, not assumed. The same retained theorem
gives V=M for a nonconstant phase once finiteness is established.
Sufficiency for every original marker is given in Section 6 below.

When a is constant, the full scalar criterion, finite-image criterion
and primitive vacuum in the statement are exactly the
[retained radial theorem](GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md).
It allows arbitrary g,B and retains the relative Gamma endpoint. It
remains to exclude the scalar premise for every nonconstant a.

## 2. Exact full anti-relative source on the normalized cover

Assume a is nonconstant, work first over C and set r=1/t. Define

    d=1-2t a(S), q=a(S)-r/2, V^2=q,
    c(t)^2=-r/2, V=c(t)sqrt(d),
    U=tg(S)+t^2B(S)^2/(2d),
    phi=-S+g(S)/r-B(S)^2/(4rV^2)=-S+U,
    eta=c(t)dS/V.                                    (V4)

The selected sqrt(d) has constant t-coefficient one. Work on the harmless
finite square-root time field F=C(t,c), embedded using this branch in
Laurent-sqrt(t) series. Rational unit quotients and the comparison below
are taken on that finite time cover. The curve is connected
and generically smooth: a(S)-r/2 is squarefree over C(r), since a
common root with a' would identify transcendental r with a fixed
critical value; its simple zeros also prevent it from being a square.

Coefficientwise ordinary Gaussian completion is EXACT:

    Phi_P0[1]=E_Gamma[d^(-1/2)exp(U)].                (V5)

Here E_Gamma[S^m]=m!, with endpoint S=0. This is not an analytic-integral
assumption, and that endpoint will not be discarded.

Let A=F[S,d^-1]. Its elements have Laurent expansions on this square-root time cover,
with finite principal part and polynomial S coefficients, since d is
a t-adic unit.
The anti-invariant one-forms on the open cover V!=0 are A(S)dS/V,
and its anti-invariant functions are F(S)/V with F in A.
The relative endpoint divisor is the TWO sections

    D0={S=0,V=±sqrt(a(0)-r/2)}.

They are distinct generically. A function F(S)/V vanishes on D0 exactly
when F(0)=0; evaluation is legitimate because d(t,0) is a nonzero
rational time scalar.

Define the full anti-relative functional by

    Lambda(A(S)dS/V)
        =c^-1 E_Gamma[A(S)d^(-1/2)exp(U)].            (V6)

All coefficients use finitely many original factorial moments. It is
well-defined on the entire rational ring A, not only on the seed.
In particular Lambda(eta)=Phi_P0[1].

For a relative function F/V with F(0)=0, direct differentiation gives

    (d_S+d_S phi)(F/V)
       =[F'+(-1+U_S-(partial_S d)/(2d))F]dS/V.

Its image under Lambda is

    c^-1 E_Gamma[(partial_S-1)
                         (F d^(-1/2)exp(U))]=0,      (V7)

because the EXACT Gamma identity is
E_Gamma[H'-H]=-H(0). Thus the endpoint condition is essential.
The formula proves descent to the full anti-relative de Rham quotient.

Time differentiation in (V6) is coefficientwise and horizontal. The
derivative is taken at FIXED S; on curve coefficients,

    (V_t/V)|S=-1/(2td), c'/c=-1/(2t),
    nabla_t(A/V dS)
       =[A_t+(U_t+1/(2td))A]dS/V.                   (V8)

It preserves F(0)=0 in relative degree zero. On amplitudes multiplying
eta this is partial_t+U_t+a/d, exactly the derivative of (V5).
Consequently all parameter-covariant iterates are retained.

If the actual unit is pexp(psi), gauge phi to phi-psi and normalize
eta/p. The gauged functional sends this seed to one and all its
positive covariant iterates to zero. Its restriction to their
differential span is a genuine rational unit quotient. Other marker
values need not be rational, which causes no difficulty.

No cleanliness or source-localization isomorphism is assumed. When
B!=0 the generic branch points V=0 are phase-irregular: B cannot vanish
at a generic root of a(S)-r/2. When B=0 the phase is regular there and
the marked points may remain; alternatively their Gysin residue
dV/V is invariant under V->-V, so its anti-invariant part is zero.
The relative endpoint anti-line is different and MUST remain.

## 3. All component orders for degree(a)>=2

Let L=deg a>=2 and normalize v(r)=1. The scalar c has order 1/2.
At a nonconstant finite S residue, or a finite centre where a is a
unit, eta has order at least 1/2; a further displacement only improves
it. At a root alpha of a of multiplicity m, put N=v(S-alpha)>0.

* If mN<1, then v(V)=mN/2 and
      v(eta)>=1/2+(1-m/2)N
        = (1-mN)/2+(mN)/m
        >=min(1/2,1/m)>0.
* If mN>1, then v(V)=1/2 and v(eta)>=N>0.
* At mN=1 without extra cancellation, v(eta)>=1/m.
* With cancellation, v(V)>1/2 and v(a')=(m-1)/m.
  The EXACT relative identity 2VdV=a'dS gives
      eta=2c dV/a',
      v(eta)>=1/2+v(V)-(m-1)/m>1/m.                 (V9)

At infinity v(S)=-N<0, a dominates r and

    v(eta)>=1/2+(L/2-1)N>=1/2.                      (V10)

These are ALL possible valuation centres and remain valid under every
phase-dependent refinement and scalar gauge. Division by p adds
deg p to the orders. They therefore apply after marking all phase
fibres 0,1,infinity, phase poles, branch sections and endpoint sections.

The generic horizontal form is regular at a simple branch point.
At infinity it is logarithmic when L=2 and holomorphic for L>=3.
Indeed for L=2m the coordinate 1/S gives dS/V proportional to
x^(m-2)dx; for L=2m+1 the coordinate S^-1/2 gives
x^(2m-2)dx. Thus every generic seed pole in this case is logarithmic
or holomorphic.

## 4. Affine a: the actual top coefficient and every infinity scale

Now a=a1S+a0, with a1!=0. Use the global cover coordinate V:

    beta=(r/2-a0)/a1, S=beta+V^2/a1,
    eta=(2c/a1)dV,
    phi=-S+g(S)/r-B(S)^2/(4rV^2).                  (V11)

Both relative endpoints are bounded at r=0. In these coordinates
time transport is the FIXED-S lift

    partial_t|S=partial_t|V+(r^2/(4V))partial_V.     (V12)

This lift preserves S and the endpoint condition. Using only the
partial derivative at fixed V would give the wrong relative connection.

The positive-power part of phi is an EVEN polynomial in V. Every
positive-power coefficient has valuation at least -1. Its ACTUAL
highest degree D is even, at least two, and its coefficient C_D has
valuation -1 or zero. Here is the complete cancellation check.

Let m=deg B. If m>=2, the possible leading competition is between
g of degree 2m-1 and B^2/V^2, at degree V^(4m-2). Cancellation requires

    deg g=2m-1, g_(2m-1)=b_m^2/(4a1).

If g has larger degree its leading term is alone; if smaller, the
B term is alone. In the cancelling case the coefficient of
V^(4m-4) is exactly

    C/r-b_m^2/(8a1^(2m)),
    C=g_(2m-2)/a1^(2m-2)
          -b_m b_(m-1)/(2a1^(2m-1))
          +b_m^2 a0/(4a1^(2m)).                     (V13)

The nonzero constant cannot cancel. As 4m-4>=4, the term -S does
not affect it. If m<=1 and deg g>=2, the g-leading term is alone.
Otherwise the V^2 coefficient is

    -1/a1+[g_1/a1-b_1^2/(4a1^2)]/r,                (V14)

with missing coefficients zero. This has nonzero constant term.
Thus infinity remains phase-irregular for EVERY parameter.

The only negative V-power is V^-2, with coefficient
-B(beta)^2/(4r) of valuation at least -1. All other scalar terms,
including -psi, are spatial constants.

Consider v(V)=-M<0. If M<1/2, the unzoomed volume has positive order
1/2-M, improved by any further zoom. At M=1/2, a nonconstant leading
phase face has positive pole action R0>=1 because the highest actual
positive-power coefficient has valuation at most zero. A finite-phase
unit-centre ball must cancel that nonconstant polynomial face at a
root of finite multiplicity mu, so its additional scale epsilon
satisfies epsilon>=R0/mu>0. Its volume has order epsilon>0.
A scalar face dominating alone has no finite-phase unit-centre ball.

If M>1/2, the highest positive monomial is strictly below every lower
positive monomial: the power gaps are at least two and the worst
coefficient disadvantage is one, giving strict gap at least 2M-1.
The V^-2 term has order at least -1+2M>0. Only a scalar gauge can tie
the highest monomial. The resulting nonzero unit-centre roots are
simple. With R0=DM-v(C_D)>=DM, exact differentiation of the phase
therefore forces epsilon>=R0. Hence

    v(eta)>=1/2-M+epsilon
             >=1/2+(D-1)M-v(C_D)>0.                 (V15)

Finite V centres, including V=0 and both endpoints, have volume order
at least 1/2 before any shrinking displacement. Thus all finite-phase
components on every refined model have strictly positive seed order.

The form dV has a DOUBLE pole at spatial infinity. Its descending tail
is handled directly: choose the DESCENDING coordinate T=-C_D V^D,
so the leading exponential is exp(-T). Then, with a local unit-root choice,

    eta=(2c/(a1D))(-C_D)^(-1/D) T^(1/D-1)dT,

whose scalar order is 1/2-v(C_D)/D>=1/2. Exact inverse-phase coordinates
have the same leading factor. Every fixed additional algebraic power
is integrable against descending exp(-T). The M=1/2 argument treats
the crossover where lower positive monomials can still compete.
This is an irregular-tail argument, not a logarithmic-pole bound.

## 5. Complete relative comparison and scalar contradiction

Use the complete proper marked-curve comparison retained in the
[H21 fixed-puncture section](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound):
its three underlying comparison facts keep full component cohomology,
finite/finite annulus complexes, finite/pole tails and all gluing maps.
The moving-cover valuations here supply its new hypotheses; a
fixed-puncture-only lemma is not being applied verbatim.

Mark the relative endpoint divisor D0 throughout. At an ordinary time,
the full relative rapid-decay pairing is perfect. Its additional
anti-invariant endpoint line is represented by a path from one endpoint
to a descending sector minus its deck image. It is not discarded and
no splitting of the relative extension is assumed. This is the same
explicit endpoint-frame/cone construction as in the
[retained relative proof](GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md#2-an-endpoint-sensitive-leading-coefficient-bound).
The finite group anti-invariant summand commutes with the relative cone.

After finite base change and resolution, relative endpoint sections
meet separated smooth marked components. If their gauged action has
a nonzero polar part, their formal exponential character contributes
no trivial quotient. If it is regular, they are allowed boundaries
of finite component chains. The seed is regular at these generic
endpoints, and its positive component factor bounds paths ending
there as well as closed paths. Endpoint collisions with branch or
phase-pole sections are resolved; (V9)-(V10), or all finite-V bounds
and (V15), apply to every component introduced.

At a finite/finite node both vertical orders are positive, so its
annulus integral is a positive power times at most logarithmic length.
At a finite/pole node the positive finite-side factor survives
descending exponential integration, which absorbs each fixed extra
pole order. At V=0 in the affine case the form is regular; its
finite/pole estimate uses that finite-side factor. The double pole
at infinity is instead covered by the explicit tail in Section 4.
Algebraic phase/inverse branches have only algebraic growth, dominated
by the deep exponential cutoff. No endpoint, neck, crossover or
puncture contribution is omitted.

The precise quantitative path input is
[analytic_core.tex](../../../unpublished/gaussian-boundary/analytic_core.tex),
Lemma lem:curvelength. The full ordinary-point passage from a cyclic
horizontal quotient to a period is
[seed_period_comparison.tex](../../../unpublished/gaussian-boundary/seed_period_comparison.tex),
Lemma lem:seedfunctional, with the endpoint cone/frame just described.
The retained fixed-high-pole tail argument is also explicit in
[the variable-W-square proof](GAUSSIAN_VARIABLE_W_SQUARE_FINITE_IMAGES.md#3-complete-neck-and-high-pole-tail-estimates).
The [primary manifest](../../../unpublished/gaussian-boundary/sources/laurent_primary_manifest.json)
identifies the three formal-nearby-cycle, finite/pole-node and ordinary
pairing inputs. These imports retain their existing review and
primary-source scope.

The phase is spatially nonconstant. For B!=0 it has a generic branch
pole; for B=0, -S+g(S)/r cannot be constant. In the affine case Section 4
also proves permanent irregularity at infinity. After removing finitely
many time poles the indicated affine source is smooth; spatial
nonconstancy makes the time and phase functions algebraically
independent. Thus the relevant exponential curve comparison applies. The gauged actual functional
(V6)-(V8) would supply a rational cyclic unit quotient; at an ordinary
parameter its extension is represented by a full relative period.
Every zero-exponential seed period just shown decays at r=0, whereas
that quotient assigns the seed the constant one. This contradiction
excludes EVERY nonconstant a, with arbitrary g,B,p,psi.

All arguments so far are over C. A hypothetical scalar identity over
an arbitrary characteristic-zero K uses only finitely many coefficients
of a,B,g,p,psi. Their finitely generated field embeds in C, preserving
nonconstancy of a and the full formal identity. The complex
contradiction proves necessity over K. Sufficiency is algebraic below.

## 6. Every original marker, safety and recomposition

In the remaining branch P0=c+bY, let gamma be the minimum strict
angular weight of a nonzero R after choosing its positive orientation.
For every original marker Q of degree D, angular balance and ordinary
Gaussian translation give EXACTLY

    Phi_P[Q]=exp(ct+b^2t^2/2)
        sum_(n=0)^floor(D/gamma) t^n/n!
                 E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^n].        (V16)

For R=0 only n=0 is used. Every summand is polynomial in t. Thus
M is contained in exp(psi)K[t], and the original unit exp(psi) proves
equality. This proves sufficiency and primitivity without replacing
original markers by an enlarged presentation.

If b=0 the normalized degree is at most floor(D/gamma). In particular
every safe member has original cutoff m>deg Q (C=1), and degree-d
perturbation r-jets vanish beyond time degree (d+1)r by the
[retained mixed-cutoff theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md).
If b!=0 and R!=0, writing M_Y=deg_Y R gives the valid bound
D+(1+M_Y)floor(D/gamma); for R=0 the bound is D.

Polynomiality of the original unit, purity and safety all force
c=b=0 by (V2), so the core is zero. Conversely zero core is safe
by (V16). The safe-recomposition ideal follows from the
[general recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md):
a safe shift exists exactly when b=0 in the finite branch, and its
unique value is c. No general safety-to-purity premise is used.

## 7. Full finite images through original degree four

The [invariant-core boundary theorem](GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md)
is a FULL-IMAGE statement: finite M for any invariant core P0(S,Y)
forces the ordinary boundary G(Y)=P0(0,Y) to be affine. Combining it
with this theorem classifies every invariant core of original
Gaussian spatial degree at most four, with arbitrary one-sided R:

    M_(P0+R) finite iff P0=c+bY.                    (V17)

Indeed after the boundary elimination a quartic core has exactly

    P0=g2 S^2+g1 S+c+(b1 S+b)Y+a1 S Y^2.

It is quadratic in Y, so the main theorem forces
g2=g1=b1=a1=0. Finite vacua are primitive with exponent ct+b^2t^2/2,
and safety holds exactly for P0=0. The literal constant-phase
exception remains as above.

This corollary does NOT infer a scalar-only quartic criterion or
pure-quartic-implies-zero. A pure unit alone does not provide the
all-marker boundary reduction. Within the finite branch, purity
of course forces c=b=0.

After the same boundary reduction, degree-five support not
automatically covered by the
[S-affine](GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md)
or quadratic-Y architectures
has form

    P0=S H3(Y)+S^2(alpha Y+beta)+bY+c,

with a nonzero cubic coefficient in H3 and (alpha,beta)!=(0,0).
This records an analytic frontier, not a classification of every such
parameter choice. Other parity-defined subfamilies retain their
separate established results.

## 8. Bounded internal support

Run

    python verify_gaussian_variable_radial_quadratic_fibre_finite_images.py \
      --output gaussian_variable_radial_quadratic_fibre_finite_images_verification.json

from this directory. The
[checker](verify_gaussian_variable_radial_quadratic_fibre_finite_images.py)
retains exact original Wick/Y-completion controls, full-source conjugation,
endpoint and moving-lift identities, symbolic all-multiplicity volume
algebra, the affine top-coefficient cancellation, and finite-branch
original markers. Its
[receipt](gaussian_variable_radial_quadratic_fibre_finite_images_verification.json)
pins these files and the retained proof dependencies. Bounded moment
controls do not prove the unrestricted curve comparison or the theorem.
