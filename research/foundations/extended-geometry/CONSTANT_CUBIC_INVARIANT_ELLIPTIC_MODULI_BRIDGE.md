# Constant cubic invariants reduce to the fixed-shape classifications

Work over a characteristic-zero field K in the Hopf sphere algebra
A=K[x,y,t]/(t^2+2xy-1). Let f=P/x be regular, with
P=C+Bt+Dt^2+Et^3 and C,B,D,E in K[x], E not identically zero.
Regularity means C0+D0=B0+E0=0. Set

    Jc=3BE-D^2,
    I=27CE^2-9BDE+2D^3,
    Disc=D^2 B^2-4EB^3-4D^3 C-27E^2 C^2+18EDBC.

Call the cubic invariant constant when either I=Jc=0 identically,
or the rational projective pair [I^2:Jc^3] is a constant point of
P^1(K). The zero pair is stated separately and is not a projective point.

**Proposition ConstantInvariantFixedShapeBridge.** If E0!=0, the
invariant is constant if and only if there exist s,h in K[[x]],
s0=-1,h0=2, and a CONSTANT lambda in K such that

    P=E h^3 Phi_lambda((t-s)/h),
    Phi_lambda(v)=v(v-1)(v-lambda).

If E0=0 and the invariant is constant, then all four center
coefficients C0,B0,D0,E0 vanish.

**Corollary ConstantInvariantLowDegreeSafety.** Every pure regular
mixed cubic with E not identically zero, constant invariant, and least
degree5 or7 is safe. This uses the fixed NONCENTRAL shape classification
when lambda!=1/2, and the CENTRAL invariant-zero low-degree theorem
when lambda=1/2. The zero-center branch is safe in every degree.
This does not exclude varying invariant, or arbitrary higher height.
The hypothesis E not identically zero keeps the lower-height case
outside this theorem.

## 1. The discriminant identity and the nonunit boundary

Direct polynomial algebra gives

    I^2+4Jc^3=-27 E^2 Disc.                            (B1)

Suppose E0=0 and D0!=0. Regularity gives B0=0 and C0=-D0;
therefore Jc0=-D0^2, I0=2D0^3. The constant projective pair must
have I^2/Jc^3=-4, so (B1) implies E^2 Disc=0. Since E is nonzero
and K[x] is a domain, Disc=0. But Disc(0)=4D0^4!=0, a
contradiction. Hence D0=0 and regularity gives C0=B0=0.
Every pure zero-center mixed cubic is safe: it belongs to K[x,t],
and the first-four interval-moment certificate for cubics forces its
x-constant term to vanish, leaving strictly positive angular weight.
No classical analytic interval theorem is needed for this height bound.

## 2. A formal affine trivialization over K itself

Assume e0=E0!=0. Depress the cubic using

    c=-D/(3E),  p=Jc/(3E^2),  q=I/(27E^3),
    P/E=(t-c)^3+p(t-c)+q.

Here c,p,q are formal series over K. Their constant terms cannot
have p0=q0=0: P0/E0=(t^2-1)(t+D0/E0) has the two distinct roots
-1 and1 and is not a cube. In particular the identically zero pair
I=Jc=0 cannot occur at this unit-center boundary.

If p0!=0, let rho=(p/p0)^(1/2), with rho0=1. Such a formal
square root exists uniquely in K[[x]] in characteristic zero.
Constant [I^2:Jc^3] implies q^2/p^3=q0^2/p0^3.
When q0=0 this makes q=0; otherwise q/(q0 rho^3) is a formal
series of square1 and constant1, hence is1. Thus in either case

    p=p0 rho^2,  q=q0 rho^3.                          (B2)

If p0=0, then q0!=0, so the constant projective pair forces
Jc identically zero. Set rho=(q/q0)^(1/3), rho0=1; again(B2)
holds. All formal coefficients belong to K, without extending constants.

Put r0=-D0/E0, lambda=(r0+1)/2, and c0=r0/3. The fixed
monic depressed cubic is the translate of
(t+1)(t-1)(t-r0). Consequently the formulas

    s=c+rho*(-1-c0),  h=2rho

give s0=-1,h0=2 and exactly P=Eh^3 Phi_lambda((t-s)/h).
Conversely an affine presentation of this form gives
p=h^2 p_lambda, q=h^3 q_lambda; its projective invariant is constant.
This proves both directions, including repeated-root lambda0,1
and the p_lambda=0 exceptional complex shapes.

## 3. Safety consequence and the elliptic connection

If lambda=1/2, q_lambda=0, hence I=0 and the central theorem applies.
For every other lambda, the noncentral fixed-shape theorem gives the
unsafe degree spectrum6,9,12,...; a pure degree5 or7 must therefore
be safe. Section1 handles E0=0. This proves the corollary.

When Disc is nonzero as a polynomial, the smooth projective model of
w^2=P(x,t), regarded over K(x), is an elliptic curve. Explicitly set
X=Et and Y=Ew, then translate X by -D/3. Its short Weierstrass
coefficients are Jc/3 and I/27, so its classical invariant is

    j=6912 Jc^3/(4Jc^3+I^2).                           (B3)

Thus the remaining nonsingular degree5/7 cubic search must have
NONCONSTANT elliptic j. This is a reduction of the remaining frontier,
not a theorem excluding all varying-j phases. If Disc=0 and (I,Jc) is not the zero pair, the projective ratio is -4.
The zero-pair/triple-root case is covered separately by the definition
and the nonunit-boundary argument.
This Weierstrass change is defined over K(x) itself. The formal affine
trivialization in Section2 was also proved directly over K[[x]].

For classical context, Andrew Sutherland's MIT18.783 Fall2025 Lecture13,
Definition13.11 and Theorems13.13--13.14, printed4--6, give the usual
short-Weierstrass j and its algebraic-closure isomorphism criterion:
https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec13.pdf
The actual displayed statements and proof were read. This is context;
(B1)--(B2) and the affine equivalence are proved directly above.
No historical-priority or external-review claim is made.

## 4. Mathematical source boundary and reproduction

The noncentral implication uses only the actual fixed-shape classification
in NONCENTRAL_CONSTANT_AFFINE_CUBIC_SHAPE_CLASSIFICATION.md, whose
hypothesis is lambda!=1/2. The midpoint implication uses only the central
I=0 theorem in INVARIANT_ZERO_MIXED_CUBIC_LOW_DEGREE_CLASSIFICATION.md.
The zero-center implication uses only the direct first-four-moment
certificate in Section2 of
MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md; its analytic
interval-rigidity discussion and conditional off-center theorem are not
premises here. No purity-to-constant-invariant implication is asserted.

The zero-pair exception is genuine: P=x^3 t^3 gives the regular pure safe
phase x^2 t^3 of least degree5, with I=Jc=Disc=0. Its invariant is covered
by the separate zero-pair definition, not a projective point of ratio -4.
For Disc!=0, the remaining UNSAFE pure degree5/7 search has nonconstant j.
Constant j is not asserted to trivialize an elliptic curve over K(x);
quadratic or higher twists are retained.

Run the portable exact/source support with

    python3 research/foundations/extended-geometry/verify_constant_cubic_invariant_moduli.py --repo-root . --output /tmp/cubic-invariant-moduli.json

The checker proves displayed polynomial identities by exact expansion and
checks literal finite formal-series and exceptional-shape controls. It
authenticates the three narrow mathematical source roles above. The
written proof supplies the all-field formal existence, necessity and
safety arguments; finite controls do not prove these universal statements.
The MIT citation is classical context only and is not a mathematical
source prerequisite. No formal-verification or independent-replay
assurance is claimed.

## 5. A necessary central branch for prime-degree unsafe phases

**Corollary ConstantInvariantPrimeUnsafeCentral.** If a pure UNSAFE
phase in this theorem has PRIME least degree, then E0 is a unit,
C0=D0=0 and I=0. Its generic elliptic curve has j=1728.

Indeed the nonunit-center branch is safe. In the unit branch, every
noncentral unsafe degree is one of6,9,12,..., all composite, by the
noncentral classification. Thus lambda=1/2. The normalized third root
is r0=2lambda-1=0, so D0=C0=0, and q_lambda=0 gives I=0.
Regularity gives B0=-E0, hence Jc(0)=-3E0^2 is nonzero. Therefore
Disc is not identically zero by(B1), and(B3) gives j=1728.
This is a necessary condition and supplies no unsafe phase, no minimum
degree and no classification of the remaining central branch.
