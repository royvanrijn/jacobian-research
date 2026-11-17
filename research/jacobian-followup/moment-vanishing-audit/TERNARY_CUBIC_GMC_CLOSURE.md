# Ternary cubic GMC: repeated lines complete the degree-three proof

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

Current lower-bound refinement,15 November2026: the [odd-invariant proof](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md)
raises the least global pure-prefix interval to[20,1855] and proves
homogeneous lower14, by algebraic existence without displayed coefficients.
Earlier lower10 statements remain valid explicit-witness evidence.

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

15 November 2026. Internal written proof with exact auxiliary algebra.
No independent specialist replay, external review or formal verification.
The checker does not certify the period lemmas or the global case cover.

## Results

**G3REPEATED20261115.** Every ternary cubic with nonreduced leading
form L^2*M or L^3 is GMC-safe. If all positive pure moments vanish,
it has a Gaussian W-linear frame.

**G3CUBICALL20261115.** For any complex polynomial P of total degree
at most3 in three independent standard real Gaussian coordinates,
all positive pure moments vanish if and only if BOTH:

    there exists v!=0 with B(v,v)=0, grad F3(v)=0, F2(v)=0;
    E(P^j)=0 for j=1,...,888.                       (1)

Every such P satisfies E(QP^m)=0 for m>3*deg Q, for every polynomial Q.
The coefficient3 is sharp. This closes OP-GMC3-CUBIC as an internal
written theorem. Equation(1) is a finite algebraic certificate with a
frame condition, not a claim that888 pure moments alone suffice for
arbitrary cubics. No explicit global pure-only prefix is obtained.

**G3DEGREEMIN20261115.** Degree4 is the least possible degree of a
three-variable GMC failure, and four Gaussian variables are the least
possible dimension of a cubic GMC failure. This combines the new
degree-three theorem with Long's explicit examples, not a new claim
of global five-term optimality.

The broader moment/vanishing programme is not complete: external review,
explicit global pure prefixes, global term-count optimality, coefficient
schemes and lower-power GVC remain open. The general GMC is still false
in dimension3 once degree4 is allowed, and cubic failures exist in
dimension4. No Jacobian/GVC dimension transfer is silently inferred.

## 1. Inputs and the one-variable exclusion used here

The [squarefree-leading theorem](LINE_CONIC_AFFINE_EXCLUSION.md) already
handles every reduced leading cubic, and forces a Gaussian W-linear
frame on any pure-null member. It remains to treat a repeated linear
factor L. The [exact frame criterion and complete W-linear theorem](RANK_ONE_WT2_CLASSIFICATION.md)
give GMC safety whenever some nonzero v satisfies the first line of(1).

For the complementary cases below, exact quadratic elimination over
C(t) produces a polynomial phase Phi(t,x), degree3 or4 in x, with a
scalar Gaussian determinant prefactor. The relevant one-variable
exclusions are the shrinking-period and two-monomial escape arguments
proved in the [rank-one WT^2 note](RANK_ONE_WT2_CLASSIFICATION.md).
We specify their uses so that no missing parameter degeneration is
replaced by a generic Morse assertion:

* If Phi=t^p*(V(x)+O(1/t)), p>0, and V is nonconstant of the full
  actual x-degree, zero-exponential scalar periods shrink. All critical
  multiplicities are retained in the inverse-phase cut construction.
* If the leading moving term has degree n<d, and the retained degree-d
  term has lower t-power, its d-n escaping critical dimensions have
  scale x=t^((p-r)/(d-n))*zeta. Their action has nonzero leading
  coefficient (1-n/d)*v*zeta^n; the n-1 bounded dimensions shrink
  whenever their exponential action is zero.
* A family tending coefficientwise to a fixed degree-d phase has a
  full fixed-phase block. A scalar prefactor t^-beta, beta>0, makes
  every zero-exponential seed value decay.

All our determinant prefactors have nonnegative beta. The degree drops
below are only(d,n,p,r)=(4,1/2/3,1,0),(4,2,2,1),
(3,2,2/3,1),(3,1,2,1). Bounded plus escaping dimensions total d-1.
The retained construction is the complete one-variable cut/local-block
argument, not a count of Hessians; coalescing points are allowed. For
the contour terminology see [DLMF2.4](https://dlmf.nist.gov/2.4).
That page is not used as a substitute for the parameter-uniform bounds
and full-rank hypotheses checked in the retained internal proof.

The rational seed-submodule bridge applies to each exact eliminated
module. If all pure moments vanished, it would supply a nonzero constant
normalized period. The above complete one-variable exclusions preclude
it. Seed cyclicity and convergence of the original Wick series are not
assumed. Nonzero exponential factors are retained as full local blocks;
we do not discard them merely because a scalar leading term cancels.

## 2. A nonisotropic repeated linear form

Take L=u with E(u^2)=1, and independent standard Gaussian Y=(Y1,Y2).
Every repeated-line cubic has

    P=Y^T*A*Y+B(u)^T*Y+g(u),
    B(u)=p*u^2+q*u+r, deg g<=3, A=A^T.           (2)

For L^2*M with M independent of L, p!=0. For L^3, p=0 and g has
nonzero cubic coefficient. The original Gaussian factorizes in u,Y.
Put C(t)=(I-2tA)^-1,D(t)=det(I-2tA). Exact quadratic elimination gives

    D(t)^(-1/2)*integral exp(Phi(t,u)) du,
    Phi=t*g(u)-u^2/2+(t^2/2)*B(u)^T*C(t)*B(u).   (3)

This identity is also coefficientwise Gaussian integration near t=0,
since D(0)=1. At infinity the prefactor has order t^-beta with
beta=1,1/2,0 according as det A!=0, rank A=1 with tr A!=0, or A
is nilpotent. It never supplies growth that could cancel shrinking.

### 2a. The double-line case p!=0

The quartic coefficient is exactly

    K4(t)=t^2*[p.p-2t*(tr(A)*p.p-p.A.p)]/(2D(t)). (4)

If it vanishes identically, p.p=p.A.p=0. The direction v=(0,p) is
nonzero and isotropic, lies in the singular plane u=0, and satisfies
F2(v)=p.A.p=0. Thus the W-linear frame criterion applies immediately,
without any assumptions on q,r,g. This handles the actual possible
pure-null branch, rather than trying to pure-exclude it.

If K4 is nonzero, the following exhaustive matrix split excludes pure zeros:

| A | Effective phase at infinity | Exclusion |
|---|---|---|
| det A!=0, p.A^-1.p!=0 | t*V4(u)+O(1), deg V4=4 | Full-degree shrinking, beta=1. |
| det A!=0, p.A^-1.p=0 | t*V(u)+W4(u)+O(1/t), deg V<=3, W4 quartic nonzero | Degree-drop escapes if V nonconstant; otherwise fixed quartic. beta=1. |
| rank1, eigenvalues lambda!=0,0; p has kernel component | Leading t^2 times the square of a quadratic | Full-degree shrinking, beta=1/2. |
| Same rank1, p has no kernel component | Nonzero quartic of order t; a kernel linear source may add a leading quadratic of order t^2 | Full-degree shrinking or(4,2,2,1) escape model. A constant kernel source adds only a scalar t^2 action. |
| A!=0,A^2=0 | Leading t^3 quartic | Full-degree shrinking, beta=0. |
| A=0 | Leading t^2 quartic | Full-degree shrinking, beta=0. |

Here are the potentially delicate boundaries. For invertible A,
C(t)=-(2t)^-1 A^-1+O(t^-2), so every coefficient of Phi is O(t).
If its leading quartic cancels, (4) shows the remaining quartic tends
to p.p/(8det A)!=0; otherwise K4 would be identically zero. Hence all
degree-drop actions have order greater than t and cannot be canceled
by a scalar O(t) term.

A rank1 nonnilpotent complex symmetric2-by-2 matrix is orthogonally
diagonalizable as diag(lambda,0). If the kernel component of p is0,
the active component is nonzero and supplies a nonzero order-t
quartic. A nonconstant kernel source is linear, supplying a quadratic
at order t^2; its escapes have nonzero order-t^3 action, above every
scalar source of order at most t^2. A constant kernel source either
adds a nonzero scalar t^2 action or vanishes.

Finally a nonzero symmetric nilpotent A is eta*n*n^T with n.n=0.
If p.A.p=0 then p is proportional to n, so p.p=0 and (4) was the
frame branch. Otherwise p.A.p!=0 supplies the order-t^3 quartic.
For A=0, nonzero(4) means p.p!=0. Thus every no-frame case in the
table has the stated nonzero quartic degree or the explicit degree drop.

### 2b. The triple-line case p=0

Now Phi has exact cubic coefficient t*[u^3]g!=0. Invertible A gives
a full-degree order-t cubic, possibly with a scalar higher action
absent here, so zero-exponential periods shrink. In rank1 semisimple
coordinates, a nonconstant kernel source adds a leading quadratic
of order t^2 and yields the(3,2,2,1) model; a constant nonzero kernel
source adds scalar t^2, while the moving cubic retains full degree.
A=0 likewise gives a leading quadratic or linear of order t^2,
or a scalar t^2 plus a full-degree moving cubic.

For nonzero nilpotent A=eta*n*n^T, a nonzero n.q gives a leading
quadratic at order t^3, with nonzero order-t^7 escaping action.
If n.q=0, then q is proportional to n. If n.r!=0 there is a nonzero
scalar t^3 term, and the largest remaining moving term is linear
of order t^2 or the original cubic of order t. Its possible order-t^(5/2)
actions cannot cancel the scalar t^3. If n.r=0, B.B=B.A.B=0 and only
the full-degree moving cubic remains. This exhausts the nilpotent
boundary. Every zero-exponential bounded period has a strictly positive
shrinking power, and all escaping actions are nonzero or dominated by
a higher nonzero scalar action. Thus every nonisotropic triple-line
member is pure-nonvanishing.

## 3. An isotropic repeated linear form

Use circular coordinates L=Z, E(ZW)=E(T^2)=1. The unique isotropic
direction in the singular plane Z=0 is the W direction. If its W^2
coefficient lambda in F2 is0, the exact W-linear frame theorem applies.
Assume lambda!=0. Gaussian null rotations fixing Z put the leading
form into exactly one of

    I: Z^2*W+aZ^3; II: Z^2*T; III: Z^3.           (5)

For a double line, the coefficient of W in M distinguishes I from II;
in II the Z component of M is removed by a null rotation. Multiplying
P by a nonzero scalar normalizes the displayed coefficient. All lower
quadratic, linear and constant terms remain arbitrary.

Write P=lambda*W^2+W*H(Z,T)+G(Z,T), and G2 for arbitrary degree<=2
terms. In I, H=Z^2+bZ+cT+d,G=aZ^3+G2; in II and III,
H=bZ+cT+d and G=Z^2*T+G2 or Z^3+G2 respectively. Completing W in
S_t=tP-ZW-T^2/2 gives a quadratic Gaussian factor of order t^-1/2 and

    Psi=t*(G-H^2/(4lambda))+Z*H/(2lambda)
              -T^2/2-Z^2/(4lambda*t).            (6)

This is an exact generic-t twisted-module calculation over C(t),
not a conditional t-adic average in the null variable W at t=0.
The rational seed quotient is inherited from the original nonsingular
Gaussian module before this elimination.

Put A=[T^2]G2-c^2/(4lambda). The T quadratic coefficient of (6) is
C_T=A*t-1/2. Complete T as well. The total determinant prefactor has
order t^-1 if A!=0, or t^-1/2 if A=0, always nongrowing.

For I, the exact quartic coefficient of the resulting Z phase is

    K4=-t/(4lambda)-c^2*t^2/(16lambda^2*C_T).      (7)

If A!=0 and [T^2]G2!=0 it has nonzero order-t limit. If A!=0 but
[T^2]G2=0, necessarily c!=0 and

    K4=-t/[2(c^2*t+2lambda)] -> -1/(2c^2)!=0.

All other coefficients are O(t), so this is exactly the fixed quartic
degree-drop family from Section1, with prefactor t^-1. If A=0,c!=0,
the leading quartic is c^2*t^2/(8lambda^2). If A=c=0, the quartic
coefficient -t/(4lambda) remains and the T linear source can contribute
a leading order-t^2 quadratic; this is the(4,2,2,1) model, or a scalar
t^2 plus a full-degree quartic. Every alternative excludes a constant.

For II, the T linear source has nonzero quadratic Z coefficient t.
Its elimination gives quartic coefficient -t^2/(4C_T), of order t
when A!=0 and exactly t^2/2 when A=0. The degree is stable and full
in either case. For III, the cubic coefficient is exactly t. If
A!=0 all moving coefficients are O(t). If A=0 a nonconstant T linear
source gives an order-t^2 quadratic and the(3,2,2,1) escape model;
a constant source instead adds a scalar t^2, or vanishes. Thus III
is also excluded. In all three types the one-variable quotient has
the full actual rank degree-1 and the determinant prefactor cannot
offset shrinking. This proves pure exclusion whenever lambda!=0.

Combining Sections2/3, every pure-null repeated-line cubic must have
the required Gaussian W-linear frame. No actual possible pure-null
branch was lost in the pure exclusions.

## 4. Full cubic closure, finite certificate and sharp mixed bound

For degree3, combine the repeated-line result with the complete
squarefree-leading theorem. For degree<=2 use the retained quadratic
classification. Pure nullity therefore forces the frame condition
in(1). Conversely a frame puts P in the complete sixteen-coefficient
W-linear family; its first888 vanishing moments force all-order pure
nullity and mixed vanishing for m>3*deg Q. This proves (1) and GMC.

For sharpness take P=W+ZT-Z^3/2. It is the Gaussian-preserving
polynomial shear of W, so all pure moments vanish. The exact source
identity is

    E exp(tP+beta*W)=exp[-beta*t*(t+beta)^2/2].

Consequently E(W^D*P^(3D))=(3D)!*(-1/2)^D!=0 for every D>=1.
Only the beta*t^3 term can attain this highest t-degree. Thus the
strict cutoff m>3*deg Q is sharp on the full ternary cubic class.

The finite certificate in(1) uses a projective algebraic frame test
in addition to888 pure moments. It does not identify a global initial
pure-moment prefix. The known sharp cusp witness still proves a lower
bound10 for such a prefix; its upper numerical value remains open.

## 5. Minimal degree and dimension of known failures

[Long's primary paper](https://arxiv.org/abs/2607.18186v1) supplies

    P3=(1+Z)*(W-(2+Z)*T^2/2), Q3=Z,
    P4=(1+Z2)*(W1*(1-Z1)+W2), Q4=Z2.

Their degrees/dimensions are4/3 and3/4. In both cases the all-order
coefficient identity is E(A(Z)P^m)=m!*[z^m]A(z)(1+z)^(m-1), using
Z2 in the second case. Circular contraction and the central-binomial
series prove it; it is not an extrapolation from the checker. Explicitly,
formal circular contraction of exp(t*(1+Z)*W)*f(Z) evaluates f at
z=t/(1-t) with factor1/(1-t). In P3, averaging T first gives
f(Z)=A(Z)/sqrt(1+t*(1+Z)*(2+Z)); at that fixed point the square-root
factor is1-t, canceling the circular factor. In P4, averaging the
first circular pair gives1/(1+t*(1+Z2)); its value at the same fixed
point also is1-t. Thus in both cases the source series is exactly
A(t/(1-t)), whose coefficient of t^m is the displayed coefficient
for m>=1. These are formal identities, so analytic convergence is
unnecessary. Taking
A=1 gives zero pure moments, and A=z gives mixed m! for every m>=1.
The new cubic theorem excludes all smaller degrees in three variables,
and by ignoring extra coordinates it excludes cubic failures in dimensions
one and two as well. These examples therefore attain both asserted
degree/dimension minima. The all-dimension quadratic theorem also
excludes every degree<=2 failure, so degree3 is the least failing
degree when dimension is unrestricted. No global term-count minimum
is proved beyond the retained rank-one ansatz.

## Evidence and outstanding work

The [checker](verify_cubic_gmc_closure.py) and
[receipt](cubic_gmc_closure_verification.json) verify exact resolvent
and elimination identities, degree-drop balances/counts, sharp shear
regressions and the two known witness regressions. The
[source audit](cubic_gmc_closure_source_audit.json) identifies the internal
period lemmas and external example source. Written analytic exclusions,
the exhaustive case cover and all-order deductions are not machine proofs.

Independent specialist replay, external review, a numerical global
pure-only gate, global term optimality, coefficient schemes and lower-power
GVC remain open. Closing the cubic objective does not complete the
broader requested research programme.
