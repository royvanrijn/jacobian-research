# Hopf quadric automorphisms, exact orbit degrees and power descent

Work over a characteristic-zero field K. All assertions below are algebraic and
remain valid after every coefficient-field extension. Put

    A=K[x,y,t]/(t^2+2xy-1),    J=x partial_t-t partial_y,

and give A its ordinary polynomial degree filtration. Write deg_A(f) for the
least total degree of a polynomial representative of a nonzero f. The normalized
moment functional is

    Lambda(f)=CT_x (1/2) integral_(-1)^1
                    f(x,(1-t^2)/(2x),t) dt.                  (O1)

The integral is formal coefficientwise polynomial integration. A phase f is
pure if Lambda(f^m)=0 for every m>=1. It is unsafe if there is ONE fixed g in A
for which Lambda(g f^m) is nonzero at infinitely many positive integers m.

The retained Hopf notes already prove the moment-preserving shear and the
original contraction normalization. They are antecedents, not new invariance
claims here. This note proves the general degree and proper-power statements
needed to transport an arbitrary primitive nonpower unsafe phase along that
shear orbit. It asserts no new degree5 or7 witness.

## 1. The degree filtration has an absolutely integral associated graded ring

Put q=t^2+2xy. The associated graded ring is

    gr A = K[x,y,t]/(q).                                      (O2)

Indeed the highest homogeneous term of H(q-1), for nonzero H, is q times
the highest term of H. Thus the initial ideal of the principal ideal (q-1)
is exactly (q).

The quadratic form q has rank3. A product of two homogeneous linear forms has
rank at most2, so q is irreducible even after algebraic closure. Hence (O2) is
an integral domain over every coefficient extension. In particular, for all
nonzero f,g,

    deg_A(fg)=deg_A(f)+deg_A(g).                               (O3)

The highest nonzero symbols multiply to a nonzero symbol. Consequently

    deg_A(g^r)=r deg_A(g)                                    (r>=1),

and the only units of A are nonzero field constants. This eliminates a
possible hidden nonconstant-unit qualification from proper-power descent.

The ring A itself embeds in K[x,x^(-1),t] by y=(1-t^2)/(2x). To justify
injectivity, q-1 is a primitive linear polynomial in y over K[x,t]:
its coefficients 2x and t^2-1 are coprime. It is irreducible by Gauss's
lemma, and its quotient is a domain. Localization at the nonzero x therefore
is injective and gives the displayed Laurent polynomial ring.

## 2. Local nilpotence, the kernel and the last derivative

On the polynomial generators,

    Jx=0,  Jt=x,  Jy=-t,  J^2y=-x,  J^3y=0.

Thus J is locally nilpotent and descends to A. Each summand of J has
ordinary degree change0, so J preserves the filtration:

    deg_A(Jf)<=deg_A(f)            whenever Jf!=0.             (O4)

On the localization it is x partial_t. Its kernel there is K[x,x^(-1)].
Moreover

    A intersection K[x,x^(-1)] = K[x].

For if a Laurent polynomial f in A had a pole of order r>0 at x=0,
then x^r f=P(x) with P(0)!=0. Reducing this equality modulo x in the
nonzero quotient A/(x)=K[y,t]/(t^2-1) gives 0=P(0), a contradiction.
It follows that

    ker_A J=K[x].                                             (O5)

Let f!=0, set d=deg_A(f), and define

    nu=max{j>=0:J^j f!=0},    P(x)=J^nu f,    c=deg P.          (O6)

Then P is a nonzero polynomial in x and c<=d. When nu=0, f lies in K[x].
The degree-propagation result below concerns nu>=1.

## 3. A completely explicit eventual degree bound

For an integer k>=0 let

    phi_k=exp(x^k J).

The exponential is finite on each polynomial because J(x)=0. Explicitly,

    phi_k(x)=x,
    phi_k(t)=t+x^(k+1),
    phi_k(y)=y-tx^k-x^(2k+1)/2.                               (O7)

Its inverse is exp(-x^k J). For nu>=1, choose

    k>=d-c+1.                                                (O8)

Then

    phi_k(f)=sum_(j=0)^nu x^(kj) J^j f/j!.                    (O9)

For j<nu the corresponding term has degree at most d+kj. On the other
hand the last term has degree c+k nu by (O3). The difference is

    c+k nu-(d+kj)=c-d+k(nu-j)>=c-d+k>=1.

Thus the last term supplies the UNIQUE highest homogeneous term:

    (lc(P)/nu!) x^(nu k+c).

Since q does not divide this nonzero pure x power, its class in (O2) is
nonzero. Therefore the exact least degree is

    deg_A(phi_k f)=nu k+c.                                    (O10)

This proves a whole arithmetic progression, not a finite range of examples.
If desired, the sufficient threshold can be sharpened by writing
d_j=deg_A(J^j f), and taking any integer k>=0 with

    k>max_(0<=j<nu) (d_j-c)/(nu-j).                           (O11)

No such sharpening is needed for the theorem.

## 4. Balanced original forms and the quadric algebra

In the original polynomial ring S=K[xi1,xi2,z1,z2], use

    R=xi1*z1+xi2*z2,  Z=xi1*z2,
    W=2xi2*z1,        T=xi1*z1-xi2*z2.

Then T^2+2ZW=R^2. A form is balanced of degree e if it has bidegree
(e,e) in the two xi and the two z variables.

The balanced subring of S is

    C=K[R,Z,W,T]/(T^2+2ZW-R^2),                              (O12)

graded with all four bilinears of degree1. Here is an elementary verification,
without an unstated invariant-theory theorem. Every balanced monomial is a
product of bilinears xi_i*z_j. With a=xi1*z1, b=xi1*z2,
c0=xi2*z1, e0=xi2*z2, their only relation is a e0-b c0.
After localizing a, set xi1=1,z1=a,z2=b,xi2=c0/a; this shows that a,b,c0
are algebraically independent there. The quotient by a e0-b c0 is an
integral domain, so localization loses no possible further relation.
Changing coordinates by a=(R+T)/2,e0=(R-T)/2,b=Z,c0=W/2 gives (O12).

Dehomogenization R=1 identifies C/(R-1) with A. For every D>=0 it gives a
vector-space isomorphism

    C_D -> A_(<=D).                                         (O13)

Surjectivity follows by homogenizing a monomial x^i y^j t^l of degree<=D
as R^(D-i-j-l) Z^i W^j T^l. For injectivity, if a homogeneous H of degree D
is (R-1)Q, take the least and greatest nonzero homogeneous degrees of Q,
say a0 and b0. The least degree of (R-1)Q is a0 and the greatest is b0+1;
both extreme terms are nonzero since C is a domain. Such a product cannot
have only one homogeneous degree because a0<=b0.

Write H_D(f) for this unique homogenization when deg_A(f)<=D. Formula
(O13) also gives

    H_D(f) H_E(g)=H_(D+E)(fg).                               (O14)

A balanced form F in C_D is divisible by R in the ORIGINAL ring S iff it
is divisible by R in C: if F=RH, its quotient H is automatically balanced
of bidegree(D-1,D-1), by separating the bihomogeneous components. Hence

    R does not divide H_D(f)  iff  deg_A(f)=D.                (O15)

Such a lift is called radial-primitive. This identifies primitive original
balanced forms with least-degree homogenizations, over every extension.

## 5. Scalar proper powers descend and ascend exactly

Let F=H_D(f) be radial-primitive, so D=deg_A(f). For every extension L/K
and every integer r>=2,

    F=alpha H^r in L[xi1,xi2,z1,z2], alpha in L*,
       iff
    f=alpha g^r in A_L.                                     (O16)

For the forward direction, first prove that H must itself be balanced.
Decompose H by total xi degree. If a0 and b0 are its least and greatest
nonzero xi degrees, the extreme xi degrees of H^r are r a0 and r b0:
their extreme products are nonzero in the polynomial domain. Since F has
the single xi degree D, one has a0=b0=D/r. Repeat this argument using
total z degree. Therefore H is bihomogeneous of degree(D/r,D/r), not just
invariant up to an unspecified character. It lies in the balanced subring
C_L by (O12). Its R=1 restriction g belongs to A_L and gives f=alpha g^r.
This argument proves descent in the FULL original ring.

Conversely, if f=alpha g^r in A_L, put e=deg_(A_L)(g). The absolute
graded-domain argument (O3) gives D=r e. Equation(O14) then yields exactly

    F=alpha H_e(g)^r.

No extra R power is present, because F has least degree. These arguments
also handle scalar factors and require no unmentioned assertion about
units in S/(R-1) or characters of a torus. Both grading-extreme arguments
are explicit, and the units in A were already checked in Section1.

In particular, if the primitive F is not a scalar proper power after
ANY field extension, then f has the same property in A after every
extension. Since each phi_k is an automorphism defined over K,
phi_k f also has that property. Its primitive least-degree original lift
therefore is not a scalar proper power over any extension.

## 6. All pure and marked moments are transported

For any h(x) in K[x] let J_h=h(x)J. On the localized chart J_h acts as
xh(x) partial_t. Formal integration gives

    Lambda(J_h a)
       =CT_x [xh(x)/2 * (a(x,0,1)-a(x,0,-1))]=0            (a in A).

The endpoint evaluations are POLYNOMIALS in x, so the additional x
annihilates their constant coefficients. Applying this identity to every
iterate and summing its finite exponential proves

    Lambda(exp(J_h)a)=Lambda(a)                            (a in A). (O17)

Consequently, with f_k=phi_k(f) and g_k=phi_k(g),

    Lambda(f_k^m)=Lambda(f^m),
    Lambda(g_k f_k^m)=Lambda(g f^m)                         (m>=0). (O18)

The marker g_k is one fixed polynomial for each k; it does not vary with m.
If g belongs to K[x], then g_k=g, so the same marker works throughout
the orbit.

Let D_k=nu k+c, k satisfying(O8), and let

    F_k=H_(D_k)(f_k),    e_k=deg_A(g_k),    Q_k=H_(e_k)(g_k).

The original balanced contraction formula is

    E2(B)=(b+1)! Lambda(B|_(R=1))       for B in C_b.

Thus

    E2(F_k^m)=(mD_k+1)! Lambda(f^m),
    E2(Q_k F_k^m)=(mD_k+e_k+1)! Lambda(g f^m).               (O19)

All factorials are nonzero in characteristic zero. If f is pure and g
survives infinitely many marked moments, every F_k is an actual original
Image failure with the fixed original marker Q_k.

The unique top term in Section3 gives an explicit original radial test:

    F_k(1,0,0,1)=lc(P)/nu! !=0.                             (O20)

At that point R=T=W=0 and Z=1. This supplies radial primitivity without a
generic-point or algebraically closed field assumption.

## 7. An arbitrary original marker can first be balanced

This section is needed only when the starting data are an ORIGINAL
balanced Image failure rather than a specified sphere pair f,g.

Let a nonzero F have bidegree(D,D), with E2(F^m)=0 for all positive m.
Suppose a fixed original polynomial Q has E2(QF^m)!=0 for infinitely many m.
Decompose Q into finitely many bihomogeneous components. At least one
component Q0, of bidegree(u,v), has a nonzero contraction for infinitely
many m. If u>v that contraction is identically zero, so u<=v. Its output

    P_m=E2(Q0 F^m)

is homogeneous in z of the FIXED degree L=v-u, independent of m. Among its
finitely many degree-L coefficients choose a multi-index beta with
|beta|=L whose coefficient is nonzero infinitely often. Characteristic zero
ensures

    E2(xi^beta Q0 F^m)=partial_z^beta P_m

is nonzero at all those orders. The one fixed marker xi^beta Q0 is balanced
of bidegree(v,v). Its sphere restriction g belongs to A, and the original
normalization gives Lambda(g f^m)!=0 infinitely often.

No unbalanced polynomial is being silently identified with a function on
the quadric. This finite-component reduction may change the infinite set of
surviving orders; it does not promise every even order from completely
arbitrary original data. If the starting balanced marker already survives
every even order, equation(O18) preserves exactly that set.

Conversely (O19), together with the elementary characteristic-zero
Image-kernel identity, converts a sphere unsafe pair into a genuine
original fixed-marker Image failure. Thus safety/unsafety of balanced
original phases and their sphere restrictions agree for arbitrary fixed
markers.

## 8. The general unsafe-orbit theorem

**Theorem HopfDerivationOrbitDegreePropagation.** Let F be a nonzero
radial-primitive balanced original Image witness of degree d over a
characteristic-zero field K. Assume it is not a nonzero scalar times a
proper power over any coefficient-field extension. Let f be its R=1
restriction. Define nu,P,c by(O6). Then nu>=1, and for every integer

    k>=d-c+1

the explicit least-degree homogenization of exp(x^k J)f is a radial-primitive,
absolutely nonproper-power original Image witness of balanced degree

    nu k+c.

It has one fixed original marker for each k. If the initial sphere marker
belongs to K[x], that sphere marker is unchanged throughout the orbit;
its least-degree original lift is therefore unchanged too. All scalar
contraction normalizations are given by(O19).

To prove nu>=1, suppose otherwise. Then f lies in K[x] by(O5), and purity
forces f(0)=0. For nonzero f in xK[x], every fixed g in A has a finite
Laurent expansion in x on the chart of(O1). The least x exponent of g f^m
eventually is positive, so all sufficiently large marked moments vanish.
Such f is safe, contradicting the witness hypothesis. The zero phase is
excluded by assumption. All other conclusions have been proved above.

This is an EVENTUAL arithmetic-progression theorem along each specified
orbit. It neither asserts coverage of every degree by an arbitrary orbit
nor supplies a witness at the globally unresolved degrees5 or7.

## 9. Checks against the retained examples

For the unscaled quartic sphere phase

    f=(1+x)/(2x) * (1-t^2(1+x)^2),

one has nu=2 and J^2 f=-x(1+x)^3. Hence d=c=4, the sufficient threshold is
k>=1, and the exact orbit degree is 2k+4. This explains the even-degree
pattern through the last derivative; it is not a new even-degree witness
claim.

For the cleared a=1 hyperbolic seed P_n=4n^2 p_n, write

    D_n^*=-(1+x)^(n+3) [((1+x)^(2n)-1)/x]^2.

The coefficient of t^3 in its localized expression is D_n^*/x. Therefore

    J^3 P_n=6x^2 D_n^*,    nu=3,
    d=c=5n+3,    lc(J^3 P_n)/3!=-1.

The theorem gives EXACT degree 3k+5n+3 for every k>=1 and the radial value
-1, recovering the already proved shear transport. Its fixed sphere
marker is x, so the original marker remains Z. The separate degree-coverage
arithmetic of that construction is not enlarged here.

## 10. Sources and scope

The raw antecedents actually read are:

* QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md Section4: finite-exponential
  shear and formal boundary proof of Lambda invariance.
* TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md, the full coefficient space,
  equations(2)--(4) and the balanced-marker reduction(6): quadric
  coordinates, original contraction normalization and marker balancing.
* HOPF_LIFT_CLASSIFICATION.md Section7, Proposition7.3: associated-graded
  least degree and the full original fixed-marker interpretation.
* HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md Sections1--4: the all-n a=1
  seed and its exact all-order marker formula.

Classical LND context was checked in Daniel Daigle's author-hosted
[Polynomials f(X,Y,Z) of low LND-degree](https://aix1.uottawa.ca/~ddaigle/articles/lndeg.pdf),
introduction, printed page1, and Section3.1(b), printed page4: derivation
height is an additive degree function, and a local slice gives the
corresponding localized polynomial coordinate. Gene Freudenburg's
[Local slice constructions in k[X,Y,Z]](https://doi.org/10.18910/4752),
Osaka Journal of Mathematics34 (1997), Proposition1 on printed page757
([publisher PDF](https://ir.library.osaka-u.ac.jp/repo/ouka/all/4752/ojm34_04_02.pdf)),
assumes an integral characteristic-zero domain and Ds=1 and identifies
the ring with (ker D)[s]. On A_x the applicable slice for J is t/x,
since J(t/x)=1. These are classical context, not new logical premises:
the specialized kernel, localization and all power/degree arguments are
proved above. Derivation height nu_D is distinct from the ordinary
sphere filtration degree deg_A used in(O3).

The general kernel, explicit degree threshold, balanced-ring identification
and all-extension scalar-power descent are proved directly in this note.
The written proof and its bounded computational controls have separate
scopes. No external review, worldwide priority or formal verification is
asserted.


## 11. Every unsafe phase has derivation height at least two

The following strengthening uses only purity at orders1 and2.

**Lemma HeightOneSafety.** Every pure f with J^2 f=0 is safe. In particular
an unsafe phase in Theorem HopfDerivationOrbitDegreePropagation has nu>=2.

On the Laurent chart J=x partial_t, so write f=a(x)t+b(x), with Laurent
polynomials a,b. If their common minimal pole order were N>0, write

    x^N f=A(x)t+B(x),

where A,B are polynomials and at least one of A(0),B(0) is nonzero.
Reducing in A/(x)=K[y,t]/(t^2-1), whose 1,t are linearly independent over
K[y], forces both constants to vanish, a contradiction. Thus a,b lie in
K[x]. Let a0=a(0),b0=b(0). The first two pure moments are

    Lambda(f)=b0,
    Lambda(f^2)=b0^2+a0^2/3.

Their vanishing gives b0=a0=0 over EVERY characteristic-zero field; no
ordered-field positivity is used. Hence f lies in xK[x,t]. Multiplication
of f^m by a fixed marker with bounded negative Laurent x support eventually
has strictly positive x support, proving safety.

## 12. Derivation height at most two is the FULL quadratic-height class

Suppose J^3 f=0. Its Laurent t degree is at most2, so for a minimal common
pole order N>=0,

    x^N f=A2(x)t^2+A1(x)t+A0(x).

If N>0, reduction modulo x gives

    A1(0)=0,    A0(0)=-A2(0).

Minimality then forces a=A2(0)!=0. Write
A2=a+xB2, A1=xB1, A0=-a+xB0. The quadric relation gives

    x^(N-1)f=-2ay+B2(x)t^2+B1(x)t+B0(x).

If N>=2, reduce this equality modulo x again. The coefficient of y in
the free K[y]-module K[y,t]/(t^2-1) is -2a, whereas the remaining terms
give only a constant plus a multiple of t. This is impossible. Thus N<=1.

It follows that EVERY f in ker J^3 is of the exact regular form

    f=x^(-1)(C(x)+t B(x)+t^2 D(x)),
    B(0)=0,    D(0)=-C(0).                                 (O21)

Conversely every displayed regular form belongs to ker J^3, since its
Laurent t degree is at most2. Thus this is an exhaustive identification,
not an assumption on a selected subclass of ker J^3.

The actual retained theorem QuadraticHeightHopfClassification in
QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md classifies ALL such regular phases,
including the zero-center case and a nonzero linear-height term. Its
zero-center pure branch and its constant-center branch are safe; its
only unsafe branch has EVEN least degree, with the exact parity formula(5).
Therefore its precise hypotheses apply to every pure unsafe f with nu<=2.

**Corollary OddUnsafeHeightBound.** A pure unsafe phase of ODD least degree
has nu>=3. In particular every primitive degree5 or7 Image witness, if one
exists, has nu>=3 for the chosen J. HeightOneSafety excludes nu<=1
without a classification input; the nu=2 exclusion at odd degree uses the
retained full quadratic-height theorem explicitly.

The statement holds for every LINEAR orthogonal conjugate of J, also after
field extension. To make this qualification algebraic, let

    Delta_q=2 partial_x partial_y+partial_t^2.

For a homogeneous polynomial H of ordinary degree2r, the formal moment
formula is

    Lambda(H)= [Delta_q^r H](0)
                    /(2^r r! (2r+1)!!),                   (O22)

and it is zero at odd homogeneous degree. Direct monomial evaluation gives
the retained sphere moments. A linear map preserving q preserves Delta_q,
so (O22) proves Lambda invariance over every characteristic-zero field.
Such a map and its inverse preserve the degree filtration, hence its least
degree. If sigma is such a map, then
J_sigma=sigma J sigma^(-1) has nu_(J_sigma)(f)=nu_J(sigma^(-1)f).
The phase sigma^(-1)f remains pure, unsafe and of the same odd least degree.
Applying the previous corollary proves nu_(J_sigma)(f)>=3.

This does NOT assert the same height bound for arbitrary nonlinear
conjugates, which can change the least degree.

## 13. The first remaining mixed cubic-height box

A useful general pole observation makes the next finite-dimensional target
precise. Reduction by the monic relation

    xy=(1-t^2)/2

gives a unique normal form

    A0(t)+sum_(i>=1) x^i Ai(t)+sum_(j>=1) y^j Bj(t).          (O23)

For example choose a degree-compatible monomial order with leading term xy;
there is one Groebner generator, so precisely the monomials with no xy
factor survive. Substituting y=(1-t^2)/(2x), the coefficient at the most
negative angular power x^(-N), N>0, is

    2^(-N)(1-t^2)^N BN(t).

It has t degree at least2N and cannot cancel against any other x power.
Hence a regular phase of Laurent t degree nu has pole order at most
floor(nu/2).

In particular ker J^4, whose Laurent t degree is at most3, is EXACTLY

    f=x^(-1)(C(x)+t B(x)+t^2 D(x)+t^3 E(x)),
    C(0)+D(0)=0,    B(0)+E(0)=0.                           (O24)

Every such expression is regular: write c0=C(0),b0=B(0), and its literal
normal form is

    2c0 y+2b0 yt
      +(C-c0)/x+t(B-b0)/x
      +t^2(D+c0)/x+t^3(E+b0)/x.                            (O25)

All quotients in (O25) are polynomial. Conversely (O23) and the pole
bound give (O24). Its derivation height is EXACTLY3 iff E is nonzero.

With deg0=-infinity, every nonzero regular phase in(O24) has exact least
degree

    max(deg C-1, deg B, deg D+1, deg E+2).                  (O26)

The normal form (O25) contains no xy monomial. A nonzero top homogeneous
part of a reduced polynomial cannot be divisible by q, since q has
leading monomial xy under the same degree-compatible order. Thus no hidden
sphere reduction can lower the displayed maximum. This also treats the
constant coefficients c0,b0, whose y and yt terms have degrees1 and2.

For d>=3, the COMPLETE least-degree-at-most-d box therefore is

    deg C<=d+1, deg B<=d, deg D<=d-1, deg E<=d-2,
    C0+D0=0, B0+E0=0.                                     (O27)

It has 4d independent coefficients: in degrees5 and7 the dimensions are20
and28. Exact least degree d means that at least one of
C_(d+1), B_d, D_(d-1), E_(d-2) is nonzero; exact height3 additionally
requires E!=0.

This is a concrete FIRST remaining low-height target: classify the pure
phases in the two boxes d=5,7 and prove safety or produce an unsafe point.
It is not a claim that every hypothetical odd witness can be put in
height3 by an orthogonal coordinate change. Phases whose height is at
least4 in every chosen direction remain outside this proposed target.

The old odd-cubic class is only the subspace C=D=0 in(O24). The retained
odd-cubic low-degree and degree7 safety results do NOT classify the mixed
box. Similarly THREE_HEIGHT_QUINTIC_IMAGE_SAFETY.md treats the ODD profile
t(C+t^2D+t^4E)/x, not the mixed height0,1,2,3 box(O24). The quadratic
classification covers E=0, and the odd-only results cover C=D=0; neither
is an exhaustion of the remaining mixed coefficients.

These reductions add no new witness or degree5/7 exclusion. They identify
an exact next coefficient domain, with the old classification input and
the still-open mixed-cubic obligation separated.


## 14. Every polynomial quadric automorphism preserves Lambda

The preceding orbit transport is not restricted to exponentials. In fact

    Lambda(sigma(f))=Lambda(f)
       for EVERY sigma in Aut_K(A) and every f in A.         (O28)

Here is a purely algebraic proof, valid after all characteristic-zero
field extensions and requiring no classification of the automorphism group.

Classical geometric context is Kutzschebauch--Leuenberger,
[Lie algebra generated by locally nilpotent derivations on Danielewski
surfaces](https://arxiv.org/pdf/1311.1075v2), version2 dated24 November2014.
Section2 on printed page3 records the standard volume form
(dx/x) wedge dz; Proposition2.1 on pages3--4 gives complex H^2 dimension
deg(p)-1. The proof of Theorem3.9 on printed page11 explicitly uses the
nonzero volume class for the quadratic p=z^2-1. These are complex
geometric antecedents only: no comparison theorem or analytic cohomology
is imported below. The paper's Theorem2.7 requires deg(p)>=3, so its
automorphism-group description is NOT used for this quadric.

Define J_y=y partial_t-t partial_x and H=x partial_x-y partial_y. The
module of algebraic two-forms on the smooth affine quadric is free:

    Omega_A^2=A omega,
    omega=2y dx wedge dt-t dx wedge dy.                     (O29)

The polynomial expression is global. Differentiating t^2+2xy=1 gives

    x omega=dx wedge dt,
    t omega=-dx wedge dy,
    y omega=-dy wedge dt.                                  (O30)

These identities also show that omega generates every two-form, since the
three elementary wedge products are multiples of omega. It has no
annihilator: after localizing x, omega=(dx wedge dt)/x is a free generator,
and localization A->A_x is injective. Thus(O29) is a free rank-one
description. Smoothness follows directly because t,x,y cannot all vanish
on t^2+2xy=1.

For every f in A, direct differentiation and(O30) give

    d(f dx)=-(Jf) omega,
    d(f dy)=(J_y f) omega,
    d(f dt)=(Hf) omega.                                    (O31)

Since dx,dy,dt generate Omega_A^1, the coefficient space of exact
two-forms is precisely

    V=J(A)+J_y(A)+H(A),    d Omega_A^1=V omega.              (O32)

It is a K-vector space, not asserted an ideal. Formula(O17) gives
Lambda(Jf)=0. Angular weights wt(x,y,t)=(1,-1,0) give Lambda(Hf)=0.
The moment formula is symmetric under exchanging x and y, so the same
argument, or the linear orthogonal invariance(O22), gives Lambda(J_y f)=0.
Consequently V is contained in ker Lambda.

Conversely the weight decomposition and characteristic zero imply that
every nonzero angular-weight part of A belongs to H(A). The normal
form(O23) shows that the zero-weight part is K[t]. Hence A/V is spanned
by the classes of powers of t. For every j>=0,

    J(y t^j)=j/2 t^(j-1)-(j+2)/2 t^(j+1),                  (O33)

where the first term is zero for j=0. Thus odd powers of t vanish
modulo V, and induction gives

    [t^(2r)]=[1]/(2r+1)              in A/V, r>=0.

It follows that dim_K(A/V)<=1. Since Lambda annihilates V and
Lambda(1)=1, this dimension is exactly1 and

    V=ker Lambda.                                          (O34)

Now let sigma be any K-algebra automorphism of A. Its action on
differential forms preserves exact forms and sends the free generator
omega to u omega for a unit u in A. The units calculation in Section1
gives u in K*. Therefore(O32) implies sigma(V)=V: for v in V,
sigma(v)u omega is exact, and multiplication by the nonzero constant
u does not change V; use sigma^(-1) for equality. Apply this to
f-Lambda(f)*1, whose class lies in V. Since sigma fixes K and1,
formula(O34) proves(O28).

This argument does NOT require pointwise preservation of omega and does
not assert it. A constant nonzero multiplier is enough. The proof also
does not use a topological realization of the sphere or an external
algebraic de Rham comparison theorem.

**Theorem AllQuadricAutomorphismMomentTransport.** Every polynomial
automorphism sigma of A preserves purity, safety and unsafety. More
precisely,

    Lambda((sigma f)^m)=Lambda(f^m),
    Lambda((sigma g)(sigma f)^m)=Lambda(g f^m)       (m>=0).  (O35)

For safety use the inverse image sigma^(-1)g of an arbitrary fixed
marker. For unsafety the one fixed marker g is transported to sigma g.
The surviving set of orders is unchanged for a specified pair f,g.

If F is a radial-primitive, absolutely nonproper-power balanced original
Image witness, let f be its sphere restriction and put
D_sigma=deg_A(sigma f). The primitive least-degree homogenization
H_(D_sigma)(sigma f) is again an absolutely nonproper-power original
Image witness, with a fixed original marker obtained by homogenizing
sigma g. This uses the exact original-ring power descent(O16), the
marker balancing in Section7 and the nonzero factorial normalization(O19).

The theorem applies to ALL polynomial automorphisms, while the explicit
degree progression nu k+c and its threshold remain assertions for the
specific family phi_k=exp(x^k J). Ordinary degree need not be preserved
by an arbitrary automorphism, and this theorem by itself neither
constructs nor excludes degree5 or7 witnesses.

