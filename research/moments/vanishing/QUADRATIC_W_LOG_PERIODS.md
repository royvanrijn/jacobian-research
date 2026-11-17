# Nondegenerate quadratic W coefficients and logarithmic periods

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

14 November 2026. Internal written proof; exact algebraic regressions are
separate from its analytic argument. No external review, formal verification
or priority is asserted.

The rank-one cubic boundary left by this argument is now complete in the
[all W-linear cubic theorem](RANK_ONE_WT2_CLASSIFICATION.md). This note
retains its all-degree nondegenerate exclusions and finite bounds.

## Statements

Let

    h(Z,T)=a+bZ+cT+dZ^2+eZT+fT^2,
    delta=e^2-4df,
    P=W*h(Z,T)+G(Z,T),       E(ZW)=E(T^2)=1.

**Theorem.** If delta!=0, no polynomial G, of any degree, makes all positive
pure Gaussian moments of P zero.

**Cubic corollaries.** If deg G<=3:

1. At least one of the first 888 moments is nonzero when delta!=0.
   In the f=0,e!=0 face, 447 suffice.
2. The complete fifteen-coefficient family with f=0 satisfies GMC.
   Its first 447 zero moments imply all-order pure vanishing and GMC safety:
   e!=0 is excluded, while e=0 is the preceding completed
   [fourteen-coefficient family](MIXED_PERIOD_AT_INFINITY.md), with gate 320.
3. Within the full sixteen-coefficient W-linear cubic family, a pure-null member not
   classified by this theorem alone must have f!=0 and delta=0. A linear null
   rotation then reduces it to

       P=W*(a+bZ+cT+fT^2)+G(Z,T),    f!=0, deg G<=3. (A)

These are local-family results. Higher W powers remain outside the family,
and no exhaustive ternary circular-frame reduction is supplied.

## The logarithmic-period estimate

**Lemma.** Let Phi(t,x) be a Laurent polynomial in x with coefficients in
C(t). On a generic sector at t=infinity, the zero-exponential part of a
flat period

    integral exp(Phi(t,x))*dx/x

is O(1+log|t|). Consequently no normalized period

    (1/t)*integral exp(Phi(t,x))*dx/x              (1)

is a nonzero constant. Here zero-exponential means that the positive-power
Puiseux part of the critical or regular-end action is zero; a bounded
constant action is included. Coincident actions are grouped. The estimate
includes multiple saddles and regular punctures at 0 or infinity.

The logarithmic differential is essential. Replacing dx/x by a form with
a higher-order pole can produce algebraic growth and does not satisfy this
lemma merely by the same argument.

### Proof of the estimate, including boundary pieces

All coefficients, critical points and critical values have convergent
Puiseux expansions after a finite cover at infinity. A nonzero exponential
action has a nonzero positive-power leading term; choose a small sector
avoiding its finitely many limiting Stokes directions. The usual
one-variable inverse-phase cut construction separates these exponential
blocks. Here is a contour bound for the zero block, without a nonuniform
Gaussian approximation near a puncture.

For fixed t use the finite branched covering x -> Phi(t,x) away from its
critical values and finite values at regular ends. Cuts from these values
to the left half-plane lift to rapid-decay paths. Adjacent inverse branches
at a critical point of multiplicity m give m local differences. At a
regular end include the small circular detours about the removed point
and the adjacent branch differences when that value is critical there.
These paths and detours span the relative homology: cut the value plane
to a simply connected domain, lift its sheets, and glue their boundaries.
This description keeps the regular-end classes as well as finite critical
classes. It does not assume a Morse phase or discard a class when a
critical point approaches 0 or infinity.

Use piecewise horizontal/vertical cuts with a bounded number of pieces.
They can avoid the other finitely many critical values, and for bounded
starting values remain in Re(w)<=C, with C independent of large t.
Coincident values use the local branch star and are grouped first.
Separating unbounded action blocks is the standard exponential
steepest-descent decomposition of these lifted cuts; cancellations within
a common action do not create a new action. A zero block is represented
by the bounded-action cuts and regular-end detours. Any exponentially
small additions from other blocks do not change its power/logarithm bound.

Choose M so that all nonzero algebraic critical locations, all coefficient
transition scales and all needed end detours lie in

    |t|^(-M) <= |x| <= |t|^M.                     (2)

Such M exists because there are only finitely many rational coefficients
and algebraic branches. On this annulus every lifted horizontal or vertical
piece lies on a real algebraic curve of uniformly bounded degree D:
Im Phi=constant or Re Phi=constant, after multiplying by a power of |x|^2.
Its length inside a disk of radius R is at most pi*D*R. Indeed average
the total variation of orthogonal projections over directions; each line
meets a degree-D real algebraic curve at most D times outside exceptional
line components. The one-dimensional coarea formula gives this bound;
line components satisfy it as well.

In a dyadic annulus r<=|x|<=2r, therefore,

    integral |dx|/|x| <= length(curve in |x|<=2r)/r
                       <=2*pi*D.

There are O(log|t|) annuli in (2). Small end circles have bounded logarithmic
length. On the bounded-action representatives |exp Phi|<=exp C, so their
contribution is O(1+log|t|).

Outside (2), an irregular end is dominated by its leading Laurent monomial.
On its lifted descending tail, w=Phi is a coordinate with
|d log x/dw|<=C'/|w|. Enlarge M so this is uniform beyond all transition
scales. Integrating exp(Re w)*|dw|/|w| gives a bounded tail. At a regular
end close the path by a small circle on which Phi is uniformly close to
its bounded end action; its logarithmic length is bounded. This also
covers a saddle coalescing with a regular puncture: use the small detour,
not the nonuniform value of 1/x at the saddle.

Thus the zero-action scalar periods are at most logarithmic, including
all boundary classes. Nonzero actions cannot supply the constant block
of (1). The zero block of (1) tends to zero. If Phi is independent of x,
the quotient has only the residue class and a period is a constant times
exp(Phi(t))/t; it cannot be a nonzero constant because a rational function
has no derivative 1/t. This proves the lemma.

For context, [NIST DLMF 2.4(vi)](https://dlmf.nist.gov/2.4#vi) treats
coalescing saddles and poles as a separate uniform-asymptotic problem.
The logarithmic contour bound above avoids assuming the ordinary
nondegenerate saddle formula remains uniform there. The analytic input
for comparison of horizontal functionals with actual rapid-decay periods
is again [Hien--Roucairol 0704.1739v1, Theorems 1.1/3.7 and 2.4](https://arxiv.org/html/0704.1739v1).
Neither source is claimed to state our Gaussian exclusion or the degree
gates below.

## The smooth-conic model when f!=0

Fix a square root eta of delta, eta!=0, and put

    rho=(e+eta)/(2f),    sigma=(e-eta)/(2f),
    u=(rho*c-b+1/t)/eta,
    v=(b-1/t-sigma*c)/eta,
    k=u*v-a/f.

Then rho-sigma=eta/f and, on the Lagrange constraint z=t*h(z,T),

    U=T+rho*z+u,    V=T+sigma*z+v,
    h(z,T)-z/t=f*(U*V-k).

The coefficient of 1/t^2 in k is -1/eta^2, so k is not identically zero.
Outside finitely many t values parametrize the conic by U=x,V=k/x:

    z=(f/eta)*(x-k/x-u+v),
    T=x-u-rho*z,
    Phi(t,x)=t*G(z,T)-T^2/2.                      (3)

This is a Laurent polynomial for every polynomial G. Direct differentiation
on the conic gives

    dz/h_T = dx/(eta*x),
    dT/(1-t*h_Z)=dx/(t*eta*x).                   (4)

Hence the normalized period is (1/(t*eta))*integral exp(Phi)*dx/x.
The factor 1/t is retained, just as in the previous polynomial-phase proof.

The formal fixed point g=t*h(g,T) gives

    x=T+rho*g+u=1/(eta*t)+O(1).

All Laurent powers of x therefore have well-defined formal Laurent series
whose coefficients are polynomials in T. This checks that localization at
x is legitimate on the actual formal branch; it is not merely an analytic
parametrization unrelated to that branch.

## The WZT face f=0,e!=0

Set z0=-c/e, x=z-z0 and A(z)=a+bz+dz^2. The inverse relation is

    T=F(t,z)=(z/t-A(z))/(e*x),
    Phi=t*G(z,F)-F^2/2,
    dT/(1-t*h_Z)=dx/(t*e*x).                     (5)

Again Phi is Laurent in x. If c!=0, x=g+c/e is a formal unit. If c=0,a!=0,
x=g=a*t+O(t^2), still a Laurent unit with polynomial T coefficients.
These are the smooth-conic cases; the numerator at x=0 is identically
zero only in the following separate case.

### The formal branch on the removed component: a=c=f=0

Here h=Z*(b+dZ+eT), g=0 identically, and GL1 gives directly

    I(t)=E_T[exp(t*r(T))/(1-t*(b+eT))],
    r(T)=G(0,T).                                 (6)

Thus all other G terms and d disappear from the pure generating series.
Do not invert g on this branch. Instead put

    T0=(1/t-b)/e,    x=T-T0,
    Phi=t*r(x+T0)-(x+T0)^2/2.

Now x=-1/(e*t)+O(1) is a valid formal Laurent unit and

    dT/(1-t*(b+eT))=-dx/(t*e*x).                 (7)

This polynomial phase on C* has exactly the same logarithmic period form.
It closes the component that a blind localization at Z would have lost.

## From the formal Wick series to a constant period

This step is common to (3), (5), (7). Write kappa=eta,e,-e respectively,
D=1-t*h_Z(g,T) (or 1-t*(b+eT) in (6)), and B=t*G(g,T) (or t*r(T)).
The formal coordinate x=X(t,T) satisfies

    X_T=t*kappa*X/D.                              (8)

Define

    J(q)=E_T[q(X)*exp(B)/D],       q in C(t)[x,x^(-1)].

All its coefficients are valid formal Laurent series by the unit checks
above, and J(1)=E(exp(tP)). On the constraint Phi(t,X)=B-T^2/2.
The logarithmic twisted quotient and normalized tangent connection are

    H=C(t)[x,x^(-1)]/Dlog(C(t)[x,x^(-1)]),
    Dlog(q)=x*q_x+x*Phi_x*q,
    nabla[q]=[q_t+(Phi_t-1/t)*q].                  (9)

These are vector-space quotients and connections. Gaussian integration by
parts and (8) give J(Dlog q)=0. For the derivative identity put w=1/D and
R=X_t/X_T. Equation (8) implies

    w_t+w/t = (w*R)_T.

The difference between differentiating J at fixed T and applying nabla is
the Gaussian integral of the T-derivative of w*R*q(X)*exp(Phi(t,X)).
It is zero coefficientwise, so J is horizontal. No convergence assumption
on the original Wick series is made.

The quotient is finite dimensional. If the phase has positive Laurent
order p and negative order q, both nonzero, its basis is
x^(-q),...,x^(p-1), of dimension p+q. Reduce high powers by the top term
of x*Phi_x and low powers by the bottom term. A high reduction never
crosses below -q; a low reduction never crosses above p-1. For independence,
a nonzero exact expression with support inside that interval would require
its preimage to have largest exponent at most -1 and smallest exponent at
least zero, which is impossible. With one irregular end, use
1,x,...,x^(p-1) or the corresponding inverse powers. At a regular end the
nonzero exponent in x*q_x reduces the remaining outside powers. A constant
phase leaves the basis 1. In every case [1] is nonzero.

If all positive pure moments vanished, J(1)=1. Exactly as in the preceding
[seed-submodule argument](MIXED_PERIOD_AT_INFINITY.md), its restriction to
the submodule generated by [1],nabla[1],... is rational and horizontal,
with values 1,0,... . It extends locally to a horizontal functional on H.
For nonconstant Phi the Hien--Roucairol comparison applies to C* x C* and
represents that functional by a flat logarithmic period. The lemma excludes
its constant normalized period. For a constant-in-x phase the rank-one
connection would require Phi_t=1/t, impossible for rational Phi. Thus
all cases contradict pure vanishing and prove the theorem.

## Explicit cubic prefixes: 447 and 888

In (3), z and T have Laurent orders at most one in x, with coefficient
t-powers between -2 and 0. Therefore, when deg G<=3,

    x-orders of Phi are in [-3,3],
    t^5*Phi has polynomial t-coefficients of degree <=6.

In (5) and (7), coefficient t-powers are instead between -1 and 0 before
forming the phase, and

    t^2*Phi has polynomial t-coefficients of degree <=3.

Let m=5 or 2 respectively and L=m+1=6 or 3. Choose the actual Laurent
orders and the quotient basis from (9); its rank r is at most six. The
coefficients of x*Phi_x are A_j/t^m with deg A_j<=L, and Phi_t-1/t has
denominator t^(m+1) and numerator degree at most L.

For two irregular ends of orders p,q, reducing a connection column needs
at most p top cancellations or q bottom cancellations. A common denominator
is

    H=t^(m+1)*A_p^p*A_(-q)^q,
    deg H<=B=L*(r+1),     deg N_ij<=B,
    connection matrix=N/H.                       (10)

At a regular end omit the absent A factor; the same bound holds, including
the constant phase after harmless overestimation. Derivative columns then
have N_j/H^j with deg N_j<=jB. At minimal seed rank s<=r, the polynomial
cofactor argument of the preceding note gives scalar coefficient degrees
at most B*s*(s+1)/2. Its zeroth coefficient is nonzero, since a constant
normalized period is impossible. Hence a sufficient nonzero prefix is

    B*r*(r+1)/2+r.

For r<=6 this is at most 888 when L=6 and 447 when L=3. The bound uses
the actual Laurent orders at every coefficient point and does not promote
a generic matrix rank to special parameters. It involves no expansion
through moment 888 and no restarted standard-basis calculation.

## All-degree finite exclusions within the nondegenerate family

The same estimate applies when deg G<=D for any integer D>=3; no new
contour theorem or higher-degree search is needed. In the general conic
coordinates, the Laurent orders are at most D at each end and

    t^(2D-1)*Phi has coefficient degree <=2D.

In the f=0,e!=0 and trapped-coordinate models, instead,

    t^(D-1)*Phi has coefficient degree <=D.

The Gaussian term has t-pole order at most four or two respectively and
fits these estimates because D>=3. Thus r<=2D and the connection bound
is B<=L*(2D+1), with L=2D generally or D on the WZT face. A sufficient
nonzero initial-moment prefix is therefore

    N_general(D)=2D^2*(2D+1)^2+2D,
    N_WZT(D)=D^2*(2D+1)^2+2D.                     (11)

At D=3 these recover 888 and447. At D=4 they give2600 and1304.
These are finite exclusions for delta!=0, not all-degree GMC theorems
on the degenerate-discriminant boundary: that boundary already contains
known quartic failures. No optimality is claimed.

## The remaining W-linear cubic family

When delta=0 and f!=0 use the Gaussian null rotation
T -> T+lambda*Z, W -> W-lambda*T-lambda^2*Z/2 with lambda=-e/(2f).
It changes the quadratic coefficients to f*T^2, since e+2f*lambda=0 and
d+e*lambda+f*lambda^2=0. The extra non-W part still has degree at most
three. This proves reduction (A). When f=0, delta=0 implies e=0, the
already completed fourteen-coefficient family. Thus (A), higher W powers,
and an unproved global frame reduction were the gates left by this argument.
The rank-one case is now classified by the all W-linear theorem linked above;
higher W powers and the global frame criterion remain open. The
remaining family contains P=W*T^2, whose pure moments vanish by circular
weight and whose mixed moments eventually vanish by the same weight count.
A blanket pure-null exclusion of (A) would therefore be false; its goal is
classification and GMC safety.

The [checker](verify_quadratic_w_log_periods.py) verifies conic and Jacobian
identities, the formal coordinate unit conditions, logarithmic quotient
reductions, universal degree implications, representative connection
matrices, the finite-bound arithmetic and direct Wick regressions. The
logarithmic growth lemma and period comparison remain written analytic
mathematics, not computationally certified contour estimates.
