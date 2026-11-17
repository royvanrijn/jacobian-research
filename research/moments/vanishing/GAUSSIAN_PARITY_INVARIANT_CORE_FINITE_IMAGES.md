# Parity invariant Gaussian cores: scalar rigidity and finite actual images

Use one Witt Gaussian pair Z,W with E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1; all other basic pairings vanish. Work over a
characteristic-zero field K and put S=ZW, T=Y^2/2. All markers and
moments are ORIGINAL unweighted Gaussian Wick moments:

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

Consider either invariant core

    P0=F(S,T),                    F arbitrary;       (P1a)
    P0=c+Y B(S,T),                B arbitrary.       (P1b)

Add an arbitrary polynomial remainder R whose monomials all have strictly
positive Z-minus-W weight, or all strictly negative weight. There is no
W-degree or Y-degree bound on R. The case R=0 is included.

The scalar and full-image classification is

    Phi_(P0+R)[1]=p exp(psi), p in K[t],p(0)=1,psi in tK[t]
      iff M_(P0+R) is finite over K[t]
      iff P0=c+bY.                                  (P2)

In (P1a) this means F=c is constant; in (P1b) it means B=b is constant.
Every finite case has primitive vacuum

    p=1, psi=ct+b^2t^2/2,
    M_(P0+R)=exp(psi)K[t].                            (P3)

For nonconstant P0+R the literal fixed-marker image V equals M. For the
constant phase c, the exception is V=K exp(ct) and M=exp(ct)K[t].
In BOTH complete families, purity and safety are equivalent and require
zero invariant core, so the phase is just the strict remainder R.

This is an internal written proof extending a specific radial shape in
Long's binary-Factorial strategy. It does not infer weighted Gamma
rigidity from ordinary FC(2). The exact kernel and continuation are
proved below. The unrefereed Long/Long--Le Blanc manuscript inputs,
retained internal audits and transitive published inputs keep their
existing review limits. External review and formal verification are not
claimed. General mixed-parity invariant cores are outside this theorem.

## 1. Exact statements and algebraic specialization

The first scalar theorem is: for an EVEN-Y invariant core f(S,Y^2/2),
a vacuum p(t)exp(psi(t)), with p polynomial, p(0)=1 and psi in tK[t],
exists if and only if f is constant. In that case p=1 and psi=f t.
Consequently full actual M is finite exactly for constant cores, and
its vacuum is primitive. Arbitrary strict one-sided angular remainders
have the same classification, by original angular balance.

A second theorem handles a CONSTANT-PLUS-ODD core c+Y b(S,Y^2/2): its
vacuum has that scalar form, equivalently its full M is finite, exactly
when b is constant. The primitive unit is exp(ct+b^2t^2/2). Arbitrary
strict one-sided remainders are included. General mixed parity is not
claimed. Proof of this second theorem is Section 7.

Let L be the product-Gamma functional

    L(S^i T^j)=i! (1/2)_j.

The first theorem in particular gives STRONG tail rigidity:

    L(f^n)=0 for all sufficiently large n  ==> f=0, (G1)

for every f in K[S,T], K characteristic zero. In the original Gaussian
with E[ZW]=E[Y^2]=1, S=ZW and T=Y^2/2 realize this exact functional.
The substitution K[S,T] -> K[Z,W,Y] is injective.

For the tail assertion alone, it suffices to exclude nonzero f with
every positive moment zero. If a
tail starting at N vanished, the nonzero polynomial f^N would have every
positive moment zero. This reduction preserves the same Gamma functional;
it does not change the shape or insert an n-dependent marker.

Here is the exact specialization used for the stronger scalar theorem.
Fix the finite coefficient supports of the hypothetical nonconstant f,
p and psi. Impose the COUNTABLY MANY coefficient equations in

    L(exp(tf))=p(t)exp(psi(t)), p(0)=1, psi(0)=0.

Each is a polynomial equation with rational coefficients in these finitely
many unknown coefficients: the Gamma moments are rational. In the nonflat
case preserve total radial degree D and positive angular degree d by
linear upper-degree constraints and an auxiliary inverse of a_d=[x^d]A.
It follows also that f_D!=0. In the flat case first divide f by its
nonzero flat leading constant, rescaling time in p and psi, so A=1;
this affine linear condition preserves D. These constraints are imposed
on ORIGINAL coefficients of f in S,T, so the structure (G2) is preserved.

The given complex point makes this ideal proper. Noetherianity reduces
the ideal to finitely many generators, and the Nullstellensatz over Qbar
gives an algebraic point satisfying EVERY original coefficient equation.
The p/psi degrees may drop; they are still finite polynomials and their
vacuum is entire. No preservation of monodromy, critical values, or
finer degrees of A_j is needed. All subsequent EIT/Lindemann arguments
are applied only to this ALGEBRAIC witness.

## 2. Uniform angular measure and the new radial shape

Take independent S~Gamma(1,1), T~Gamma(1/2,1), and put

    S=r(1-x^2), T=r x^2, 0<=x<=1.

The Jacobian is 2rx and the joint density becomes

    (2/sqrt(pi)) exp(-r) r^(1/2) dr dx.

Thus x is UNIFORM and r~Gamma(beta,1), beta=3/2. If f has total degree
D>=1, write

    P(x,r)=f(r(1-x^2),rx^2)
          =A(x)r^D+A_1(x)r^(D-1)+...+A_D(x).        (G2)

Each A_j is an EVEN polynomial of degree at most 2(D-j). Restriction of
homogeneous forms to S+T=1 is injective, so A is nonzero. Its angular
degree d may be as large as 2D. The cases A nonconstant and A constant
are treated separately.

Define rationally normalized moments and their germ

    M_n=L(f^n)/(3/2)_(Dn),
    H(w)=w^-1 sum_(n>=0)M_n w^-n.

M_0=1. Under all-positive-moment vanishing H=1/w. Under the stronger
scalar identity write a_n=[t^n]p exp(psi). The entire function satisfies
|a_n|<=C_R R^(-n) for EVERY R>0, and

    M_n=n! a_n/(3/2)_(Dn), (3/2)_(Dn)>=n!.

Thus limsup |M_n|^(1/n)=0, and H is single-valued HOLOMORPHIC on C*.
It may have an essential singularity at zero. The leaf proof below needs
only holomorphy at every nonzero leaf and single-valued continuation
away from zero; it does not require H to equal 1/w or be rational.

The exact kernel is

    H(w)=integral_0^1 C(w,x)dx,
    C(w,x)=w^-1 integral_0^infinity exp(-r)r^(1/2)
                       E_(D,3/2)(P(x,r)/w)dr.       (G3)

The factor Gamma(3/2) cancels between the original density and the rising-
factorial normalization. No irrational constant is hidden in M_n.

## 3. The new kernel interface, including global continuation

The new kernel argument is given explicitly here. All assertions below
hold uniformly on compact subsets of V={A!=0}.

The entire root filter is

    E_(D,3/2)(X^D)=D^-1 sum_(j=0)^(D-1)E_(1,3/2)(omega^jX),
    E_(1,3/2)(z)=(2/sqrt(pi))integral_0^1 exp(z(1-v^2))dv.

This holds also for D=1, without a correction term. Its bound by a
constant times exp(|X|) gives the original kernel convergence for large
|w|. At radial infinity set

    Phi_x(r)=(P/A)^(1/D)=r+p1+O(r^-1), p1=A_1/(DA),
    chi(x,Y)=Y-p1+O(Y^-1),
    F(x,Y)=exp(Y-chi)(chi/Y)^(1/2)chi_Y
          =sum_(m>=0)F_m(x)Y^-m, F_0=exp(p1).       (G4)

The square root is the exterior unit with initial value one. These are
convergent exterior Laurent series; |F_m|<=C R0^m on compact x-sets.

Split the radial integral at a uniformly exterior radius. In the initial
small-rho convergence domain, Y=Phi_x(r)=r+O(1) permits deformation of
the exterior tail to the positive Y-ray: all exponentials decay in the
intervening right half-sector, and connecting contours at infinity vanish.
Finite connectors are holomorphic by local primitives. This deformation
is used only initially; the formulas below construct its continuation.
Each root-filter summand is a sum of

    2/(D sqrt(pi)w) F_m(x) Gamma(3/2-m) I_m(rho),
    rho^D=A/w,
    I_m(rho)=integral_0^1(1-rho+rho v^2)^(m-3/2)dv, (G5)

and an entire remainder. Explicitly, with s_m=3/2-m,

    J_m(eta)=integral_R^infinity e^(-eta Y)Y^(s_m-1)dY
            =Gamma(s_m)eta^(-s_m)+E_m(eta),
    E_m(eta)=-R^s_m sum_(k>=0)(-eta R)^k/[k!(s_m+k)].

The equality holds first at positive shape and extends by elementary
integration by parts and Gamma(s+1)=s Gamma(s). Since all s_m+k are
nonzero half-integers, |s_m+k|>=1/2, giving the EXACT all-m bound

    |E_m(eta)|<=2 R^(3/2-m) exp(LR), |eta|<=L.

For R>R0 the remainder sum therefore converges normally on every compact
rho-set, and is entire in rho. The remaining functions satisfy

    I_0=1/(1-rho), I_1=atanh(sqrt(rho))/sqrt(rho),
    I_m=[1+(2m-3)(1-rho)I_(m-1)]/[2(m-1)], m>=2.   (G6)

Their only possible finite continuation exceptions are rho=0,1.
The recurrence follows by integrating the derivative of
v(1-rho+rho v^2)^(m-3/2). It gives at most exponential growth in m on
every compact continued sheet. For m>=1,
|Gamma(3/2-m)|=pi/Gamma(m-1/2), so the Gamma factor decays factorially.
Thus the full kernel has global continuation with no finite
w exception except 0 and A(x). This is a normal-convergence proof, not a
sectorial Mittag--Leffler asymptotic assumption.

The m=0 root sum is EXACTLY

    exp(p1)/(w-A).

For m>=1 the logarithmic coefficient is

    B(w,x,rho)=1/(D w sqrt(rho))
       sum_(m>=1)(-1)^m F_m(x)(1-rho)^(m-1)/(m-1)!. (G7)

Indeed the m=0 identity is D^(-1)sum_j(1-omega^j rho)^(-1)
=(1-rho^D)^(-1). For m>=1 the recurrence writes

    I_m=P_m(1-rho)+c_m(1-rho)^(m-1)I_1,
    c_m=Gamma(m-1/2)/(sqrt(pi)Gamma(m)).

Near rho=1 on sqrt(rho)=+1,
I_1=[log(1+sqrt(rho))-(1/2)log(1-rho)]/sqrt(rho).
Gamma reflection then gives exactly the sign and factorial in (G7).
The factorially convergent series is holomorphic on
sigma^(2D)=A(x)/w, rho=sigma^2, with wA!=0. Locally near w=A!=0 the
remaining term is jointly holomorphic.

At A(x)=0 the ORIGINAL branch fills by a normally convergent integral.
On a compact nonzero w-set, shrink the x-neighborhood so
|P(x,r)/w|^(1/D)<=r/2+C(1+r)^(1-1/D). Its integrand has the uniform
integrable bound

    C r^(1/2) exp(-r/2+C'(1+r)^(1-1/D)).

For a fixed compact w-continuation path away from zero, shrink that
neighborhood until |A(x)/w|^(1/D)<1/2 along the WHOLE path. The bound
is then uniform along it, and uniqueness identifies the normally
convergent original integral with that fixed continued germ exactly
where angular deformation uses the fill. This is the comparison in the
pinned Long lemma `lem:normality`. It does NOT assert removability of
EVERY monodromy branch, or cover paths shrinking around A(x) as A(x)->0.

Since the continued H is holomorphic at EVERY nonzero leaf, the local
pole-plus-log extraction from the retained Long proof gives at each such
target leaf the FULL corrected identity (the pinned Long
`lem:local-extraction`, whose hypothesis is meromorphy of that chosen
H germ, not global rationality)

    sum_i n_i[ exp(p1(x_i(w)))x_i'(w)
               +integral_(c_i)^(x_i(w))B(w,x,rho)dx ]=0. (G8)

Discarding the integral correction would be invalid. Continue (G8)
with every lower endpoint, history path and root sheet retained, as in
the pinned Long `thm:transport`. Holomorphic single-valued H on C* is
sufficient for this continued zero identity. On a
base disk choose s(w)^(2D)=w; then sigma=tau/s(w), A(x)=tau^(2D), makes
this cover a product with that disk. Fixed compact history integrals are
holomorphic in w. At a nonzero critical value lambda, a short terminal
piece has a joint holomorphic primitive and contributes only NONNEGATIVE
Puiseux powers.

At every moving upper endpoint A(x_i(w))/w is identically one. Its
continuous lift has sigma=+1, the starting physical singular sheet.
Lower endpoint sheets are never reset. This condition prevents an
incorrect sign change in the pole residue. Thus every NEGATIVE critical
coefficient of a continued (G8) comes solely from

    exp(A_1(x_i(w))/(D w))/A'(x_i(w)).              (G9)

The fractional radial weight creates no new negative history exponent.

## 4. Degree audit of the retained nonflat argument

The primary source is Long's
[pinned binary-Factorial manuscript](../../factorial/image-rigidity/sources/long-fc2.tex),
commit
[8b1eb46](https://github.com/long-mathematics/factorial-conjecture-two-variables/tree/8b1eb46ff530122ed25f588b78bf9fdc98cb9577).
Its nonflat theorem, exact-kernel, corrected-leaf, transport, rational
path-module, critical coefficient and endpoint arguments are recorded in
the retained
[complete internal audit](../../factorial/image-rigidity/LONG_EIT_AND_FC2_AUDIT.tex),
section "The nonflat FC(2) proof". Exact primary identities and transitive
reference locators are pinned in the
[source ledger](../../factorial/image-rigidity/source_ledger.json).

In that source deg A_j<=D-j occurs in the original binary radial setup
and its specialization constraints. The kernel proof only requires
polynomial A_j, while the exterior coefficients need only be holomorphic
on V. The monodromy degree is d=deg A, independent of radial D; no step
of the rational path-module, inverse-derivative trace, or critical
arithmetic proof uses d<=D or deg A_1<=D-1. The optional EIT-free flat
shortcut needs a small degree and is NOT used here. Our degrees
2(D-j) are therefore compatible with the entire retained nonflat
argument, now that (G3)-(G9) supply its exact kernel interface.

For clarity, the remaining logical steps are explicit.

Let G_A be the polynomial monodromy of A, acting on its d sheets. The
rational signed path module M generated by the nonzero leaf vectors
contains the top rational component U_d. The representation-theoretic input is
[Pakovich--Muzychuk, Theorem3.1 and Section4.1/Proposition4.1](../../factorial/image-rigidity/sources/pakovich-muzychuk-07104085v2.html).
The retained Long argument checks the actual endpoint averages and
planar path containment; the irreducibility is rational, not absolute
complex irreducibility.

The unweighted inverse-derivative vector u_i=1/A'(x_i(w)) belongs to
U_d tensor C. For every proper polynomial decomposition A=P o Q,
the sum on each Q-block is zero:

    sum_(Q(c)=xi) 1/A'(c)
      =P'(xi)^(-1)sum_(Q(c)=xi)1/Q'(c)=0.

The last equality is the coefficient of x^-1 in 1/(Q(x)-xi), whose
order at infinity is at least two. These block equations are enforced
simultaneously. This is why retaining the ORIGINAL unweighted angular
dx is essential; an arbitrary density would change this vector.

If lambda!=0 is a critical value, let k be the largest local degree
above it. The coefficient beta=(beta_i) of u at exponent 1/k-1 is
nonzero, algebraic and lies in U_d tensor Qbar. By (G8)-(G9), for every
monodromy continuation the weighted vector

    beta_i exp(A_1(c_i)/(D lambda))

annihilates U_d. No smaller local degree or positive-order variation of
the exponential factor contributes at this strongest negative exponent.
Since U_d is rational, conjugate(beta) also belongs to it. Pairing gives

    sum_(beta_i!=0) |beta_i|^2
                       exp(A_1(c_i)/(D lambda))=0.       (G10)

Group equal algebraic exponents. Every coefficient is a strictly positive
real algebraic number. Lindemann--Weierstrass contradicts (G10).
Thus all critical values of A are zero.

Counting root and critical-root multiplicities then gives
A=a(x-c)^d. In OUR problem A is even and nonconstant. Hence d is even
and c=0, so A(0)=0 and A(1)=a!=0. At the nonzero endpoint leaf1, the
shrinking correction integral vanishes and its pole density is the
single nonzero value exp(A_1(1)/(D a))/A'(1). The other endpoint lies
at target0, which is excluded from this nonzero leaf. Thus the endpoint
cannot cancel, contradicting (G8). No equal-endpoint comparison or
extra amplitude theorem is needed in this shape-specific endpoint case.
This excludes every nonflat algebraic data set for which H is
single-valued holomorphic on C*. In particular it excludes the specialized
scalar p exp(psi) witness. The endpoint proof has used the original
exp(p1) residue and sigma=+1 upper sheet, so no fractional-sheet sign
or extra weighted density can cancel that residue.

## 5. The flat top and its Gamma-three-halves bound

If A is constant, normalize A=1 and use algebraic specialization as
above. For a multinomial index k=(k_1,...,k_D), put
K=sum k_j and J=sum j k_j. Its normalized radial coefficient is

    C_(n,k)=(n!/(n-K)!)
       Gamma(Dn-J+3/2)/Gamma(Dn+3/2), K<=n,

and zero if K>n. One has 0<=C_(n,k)<=1. Indeed J>=K and J<=DK, and
the smallest K denominator factors beginning at Dn-J+3/2 dominate
n-K+1,...,n. All additional denominator factors are at least one.
For each fixed k,

    C_(n,k)~D^(-J)n^(K-J).

The summable majorant product_j ||A_j||^k_j/k_j! gives uniform dominated
convergence. Only J=K survives, so

    lim_(n->infinity) M_n
      =integral_0^1 exp(A_1(x)/D)dx.                (G11)

The exponent is a polynomial with algebraic coefficients and the
amplitude is EXACTLY dx. The pinned
[Long--Le Blanc EIT I single-integral corollary](../../factorial/image-rigidity/sources/long-eit-pinned.tex),
source label `cor:single`, gives this value nonzero; a constant A_1 is
immediate. Its exact scope is a polynomial exponent over Qbar, distinct
algebraic endpoints and the UNWEIGHTED differential dx. Here the endpoints
are 0 and 1 and the parameter is 1. No arbitrary-amplitude extension is
used. The retained
[EIT dependency audit](../../factorial/image-rigidity/LONG_EIT_AND_FC2_AUDIT.tex)
records its E-function lifting, Fourier and published transcendence
inputs; the [Beukers primary mirror](../../factorial/image-rigidity/sources/beukers-math0405549.html)
is one of those transitive inputs.
For the algebraically specialized scalar p exp(psi) witness, root
growth of M_n is zero. This contradicts the NONZERO limit (G11).
The claim of a nonzero flat limit is made ONLY after specialization;
complex exponential integrals at transcendental coefficients can vanish.
This closes the flat case and proves the scalar criterion. Constants
f=c have vacuum exp(ct); polynomial-exponential rigidity gives p=1,
psi=ct. The STRONG tail result (G1) follows either by the f^N reduction
in Section 1 or by applying the scalar theorem directly to the polynomial
vacuum of a phase with an eventual zero moment tail.

## 6. Complete even-core finite images and one-sided remainders

For P=f(S,Y^2/2)+R with arbitrary strict one-sided angular R, its ACTUAL
vacuum equals that of f. Thus its vacuum is p exp(psi) exactly when
f=c is constant, and then p=1,psi=ct. If M_P is finite it supplies that
scalar identity by the [general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image).
Conversely for constant core, angular balance gives for every ORIGINAL
marker Q of degree D

    Phi_P[Q]=exp(ct) sum_(n=0)^floor(D/gamma)
                          t^n/n! E[Q R^n],

where gamma is the minimum strict angular weight after reversing a
negative remainder. For R=0 retain only n=0. Hence the full original
module is finite, primitive and equals exp(ct)K[t]. For nonconstant P
its literal V equals M by the retained original-image theorem; for the
constant phase c the literal image is only K exp(ct).

Polynomial unit, purity, safety, and zero core are equivalent in this
WHOLE even-core-plus-one-sided family. Safe members have original cutoff
m>floor(deg Q/gamma). For a separate recomposition indeterminate U, the
[general safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
gives S_P=(U-c) when f=c is constant, and S_P=(0) otherwise. This tests every
polynomial recomposition even when it leaves the original architecture.

## 7. Constant-plus-odd cores: a separate entire-growth argument

Now let P0=c+Y b(S,T), T=Y^2/2, with b arbitrary. Suppose its ACTUAL
vacuum is p(t)exp(psi(t)). Then the entire function

    exp(-ct)Phi_P0[1]=sum_(m>=0)h_m t^m

is the vacuum of Y b. Put Q(S,T)=2T b(S,T)^2. Its raw moments satisfy

    L(Q^n)=(2n)! h_(2n).

If deg b=m>=1, Q has total degree D_Q=2m+1>=3, and

    M_n(Q)=(2n)! h_(2n)/(3/2)_(D_Q n).

Every denominator factor dominates the corresponding factor of (2n)!,
so |M_n(Q)|<=|h_(2n)|. Entire growth gives root growth zero, hence H_Q
is single-valued holomorphic on C*.

Specialize the ORIGINAL coefficients c,b,p,psi and their scalar formal
identity, not a putative scalar identity for Q. Preserve m by an inverse
of a nonzero leading b coefficient. Q then has algebraic coefficients,
keeps D_Q, and its angular top is

    A_Q(x)=2x^2 b_m(1-x^2,x^2)^2.

It is nonzero by injectivity of homogeneous restriction and vanishes at
x=0. Thus it is necessarily NONFLAT. The nonflat kernel/monodromy
obstruction in Sections 3-4 contradicts holomorphic H_Q. No flat EIT
argument is used in this odd-core necessity proof.

Therefore b is constant. Conversely for constant b and arbitrary strict
one-sided R, ordinary Gaussian translation and angular truncation give

    Phi_(P0+R)[Q]=exp(ct+b^2t^2/2)
      *sum_(n=0)^floor(D/gamma) t^n/n!
                    E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^n].

Hence scalar polynomial-exponential vacuum, full finite image and
constant b are equivalent, with primitive exponent ct+b^2t^2/2.
Purity and safety both require c=b=0, so the core vanishes. For R!=0
and b!=0 the normalized cutoff is bounded by
D+(1+deg_Y R)floor(D/gamma); for b=0 it is floor(D/gamma). If R=0,
retain only n=0 in the displayed sum and use normalized degree at most D
(at most zero if also b=0), without assigning a degree to the zero R. Safe recompositions exist exactly when
b=0, and then S_P=(U-c); otherwise S_P=(0).

The odd core Y itself is finite with vacuum exp(t^2/2), so polynomial-
unit rigidity must not be confused with absence of every finite image.
GENERAL mixed-parity cores F(S,T)+Y B(S,T), with F nonconstant and B
nonzero, remain outside these two results. Multiplying a safe phase by
its Y-reflection is not a valid shortcut: safety is not an ordinary
ideal property (Z and W are safe, while ZW is not).

## 8. What is old and what is new in this reduction

The one-dimensional eventual-moment theorem for interval polynomials
and its Gegenbauer extensions are in
[Francois--Pakovich--Yomdin--Zhao, Bull.Sci.Math.135(2011)](https://www.math.bgu.ac.il/~pakovich/Publications/FPYZ.pdf),
Corollary4.1 and Proposition4.2. They already exclude nonzero
Gaussian-HOMOGENEOUS axial invariant cores by sphere reduction. No novelty
is claimed for that consequence.

Ordinary binary FC(2), the nonflat corrected-leaf/monodromy strategy,
and the EIT I input are credited to Long and Long--Le Blanc, with the
published representation and transcendence inputs identified in the
retained audit. Ordinary FC(2) alone does not imply the weighted result.
The additional work here is the exact Gamma-three-halves kernel and
global/history continuation, its degree and endpoint audit, and the
flat Gamma-three-halves normalization. Those steps address arbitrary nonhomogeneous Gamma(1,1/2) phases,
including a finite polynomial-exponential scalar gauge. They yield full
finite-image/primitivity classifications for even invariant cores and,
by the separate squared-odd growth argument, constant-plus-odd invariant
cores. No general mixed-parity theorem is claimed. Independent/external mathematical review and
formal verification remain unclaimed; no priority claim is made.


## 9. Original-image bounds, evidence and scope

The exact finite expansions in Sections 6-7 retain every ORIGINAL marker.
The stronger all-degree identity M_(P0+R)=M_P0, with its polynomial-time
kernel transports, is in the
[angular-core full-image theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md).
For a nonzero strict remainder let gamma be its least absolute angular
weight after reversing negative support. Safe members have cutoff
m>floor(deg Q/gamma), independently of the degree of R. For R=0 only
the zero exponent contributes. In the finite odd branch with b!=0,
ordinary translation gives normalized time degree at most
D+(1+deg_Y R)floor(D/gamma) for original marker degree D; with R=0
use the ordinary degree bound D instead.

At a safe member P=R, the
[mixed-cutoff finite-jet and certified-tail theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#1-the-finite-jet-theorem)
gives moment r-jets through floor((1+d/gamma)r) for perturbations of
original spatial degree at most d. An actual certificate for the
coefficient maximal ideal m, with m^nu contained in a chosen prefix ideal, leaves only intervening moments
through floor((1+d/gamma)(nu-1)) to check. For R=0 take gamma=infinity.
These are sufficient jet/tail budgets at a known safe phase, not claims
of a least ideal-generating prefix or a computed nilpotence index.

The [focused exact checker](verify_gaussian_parity_invariant_core_finite_images.py)
and [receipt](gaussian_parity_invariant_core_finite_images_verification.json)
retain the stated bounded kernel identities, original Gamma/Gaussian
fixtures and all-marker finite-branch examples. They do not certify
normal convergence, analytic continuation, the unrestricted monodromy
argument, EIT I or the all-degree theorem by extrapolation. Those are the
written arguments with the pinned dependencies above. All scalar identities
are actual original vacuum identities; no alternate W-marker map is used.

The complete even-core and constant-plus-odd-core results extend beyond
the retained [constant-fibre radial theorem](GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md).
Previously settled mixed-parity subclasses remain valid, but no arbitrary
F(S,T)+Y B(S,T) theorem is inferred. General Gaussian safety-to-purity,
full finite-vacuum primitivity outside the proved families, and generic
ranks of arbitrary mixed-parity invariant cores remain open here.
No independent replay, formal verification, external mathematical review,
publication or priority is claimed.


## 10. A coefficient-nilpotence consequence, without an effective prefix

Let A be ANY commutative Q-algebra and f in A[S,T]. The strong field
Gamma(1,1/2) tail theorem gives

    L_A(f^n)=0 eventually
      iff every coefficient of f is nilpotent in A
      iff f is nilpotent in A[S,T].                 (P4)

For necessity, reduce modulo any prime ideal and pass to its fraction
field; (G1) makes every coefficient zero there. Thus each coefficient
belongs to the nilradical. There are only finitely many coefficients.
If their nilpotence exponents are e_i, the coefficient ideal has power
1+sum_i(e_i-1) equal to zero by the pigeonhole principle. This proves
polynomial nilpotence and the converse. Nilpotence also gives eventual
vanishing against EVERY original Gaussian marker. No Noetherian
hypothesis on A is needed for (P4).

For a fixed finite coefficient family f=sum_i u_i m_i(S,T) in DISTINCT
monomials m_i over a field, let I_infinity be the ideal of all positive pure moments.
The field theorem gives radical(I_infinity)=(u_1,...,u_N). Noetherianity
then supplies some finite moment prefix I_n=I_infinity and some index
nu with (u_1,...,u_N)^nu contained in I_n. Hence the CORE coefficient
slice has finite colength. No effective value of n or nu is proved here.
With independent strict-remainder parameters, the vacuum ideal is the
extension of this core ideal; finite colength is not asserted in those
additional free directions.

This does NOT promote field safety-to-purity to rings. For epsilon^2=0
and epsilon!=0, the constant phase epsilon is safe and has vacuum
1+epsilon t, which is not one. An effective jet/tail conclusion still
requires an actual nilpotence certificate and its intervening checks.
