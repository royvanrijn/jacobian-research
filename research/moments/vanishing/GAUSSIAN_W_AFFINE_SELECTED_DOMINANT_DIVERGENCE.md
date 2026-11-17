# Actual unit divergence from a selected dominant inverse singularity

This written theorem extends the exact original-source divergence proof
from a single nonzero root to selected dominant inverse singularities,
including explicit and generic multiroot coefficient families. No
moment-prefix inference is used.

Let c(0)!=0 and deg c=m>=2. Let a be constant and h,u,k arbitrary complex
polynomials, and put

    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z).

The selected inverse z0(X) is the analytic germ

    z0=X c(z0), z0(0)=0,
    G0(z)=z/c(z), K0(X)=1/[1-X c'(z0(X))].

## 1. A selected dominant singularity criterion

Assume there is rho!=0 such that the selected z0 admits a bounded
analytic continuation to a disk |X|<R, R>|rho|, cut along the outward
radial segment from rho, with its two banks treated separately, and:

* the only selected singularity in that slit closure is rho;
* z0 tends to a finite zstar at rho, with c(zstar)!=0;
* G0'(zstar)=0 and G0''(zstar)!=0.

Then the ORIGINAL Gaussian unit f=E exp(tP) is divergent at t=0.
Consequently its full actual image is not finite over C[t], and P is
neither pure nor safe. In fact f and its fixed original marked period
f'=Phi[P] are C(t)-independent.

Here and below the slit can be shortened by reducing R toward |rho|.
No condition is imposed on critical values belonging only to other
inverse sheets inside the disk.

Proof. Use exactly the original completion and circular Lagrange source:

    d=1-2at, B=th+t^2cu/d, D=exp(tk+t^2u^2/(2d)),
    z=B(t,z)+c(z)X/sqrt(d),
    J=1-B_z-c'X/sqrt(d),
    K(t,X)=D(t,z)/[sqrt(d) J],
    f(t)=E_S K(t,tS).                               (S1)

Every fixed t coefficient in (S1) is a finite Wick contraction. These are
the actual original periods; no localized alternate trace is introduced.
The unrestricted source proof is Section 1, (D1)-(D2), of
[the canonical original-source proof](GAUSSIAN_W_AFFINE_SINGLE_ROOT_DIVERGENCE.md).

The simple critical point of G_t(z)=sqrt(d)(z-B)/c persists analytically,
with critical value rho(t), rho(0)=rho. Put

    T=t rho/rho(t), eta(t)=rho(t)/rho,
    Ktilde(T,X)=K(t(T),eta(t(T))X).

Then t(T) is analytic with derivative one and

    f(t(T))=E_S Ktilde(T,TS).                        (S2)

The distinguished singularity is fixed at rho. On the compact selected
graph outside a small critical disk, the defining z derivative is a
unit. Finitely many isolated root neighborhoods and the implicit-function
theorem give a joint analytic branch for small complex T. They glue by
uniqueness near the unperturbed graph. Inside the critical disk, the
quadratic critical point gives a joint square-root inverse, now with
fixed branch value. This constructs the same bounded fixed slit uniformly
in T. Extra critical points created at unbounded z by arbitrary h,u,k
cannot enter these bounded isolated neighborhoods even when their values
approach zero. Likewise an interior other-sheet critical value is harmless
because the selected defining derivative there is nonzero.

Near rho the kernel therefore has a uniform local form

    P(T,X)+Q(T,X)(1-X/rho)^(-1/2),                   (S3)

with P,Q jointly analytic, Q(0,rho)!=0. Indeed the amplitude D/sqrt(d)
is an analytic unit and J has a simple zero in the square-root coordinate.

Write b_n=binom(2n,n)/4^n. Choose the outer radius so close to |rho| that the entire short slit lies
inside the quadratic-critical disk. In deforming the coefficient contour,
the tiny circle about rho contributes O(epsilon^(1/2)) and vanishes.
The slit-contour estimate gives, uniformly on
a fixed closed T disk,

    |[X^n]Ktilde(T,X)|<=C |rho|^(-n)(n+1)^(-1/2).

Subtract Q(0,rho)(1-X/rho)^(-1/2) at T=0. Its remaining slit jump is
O(s^(1/2)), so the same contour gives

    a_n:=[X^n]K0
       =Q(0,rho) rho^(-n)b_n
          +O(|rho|^(-n)(n+1)^(-3/2)).              (S4)

The analytic part of the contour is exponentially smaller. Thus
a_n=Q(0,rho) rho^(-n)b_n(1+O(1/n)). In particular the selected unique
dominant singularity forbids eventual cancellation of its even coefficients.
Cauchy's estimate in T gives

    |[T^r X^n]Ktilde|<=C R0^(-r)|rho|^(-n)(n+1)^(-1/2).  (S5)

Define L_q=(2q-1)!! |rho|^(-2q)(2q+1)^(-1/2). Its successive ratio
divided by q tends to 2|rho|^(-2)>0. The coefficient of T^(2q) in (S2)
has main contribution (2q-1)!! a_(2q). Every other surviving term has
X degree 2j<2q and, by (S5), its total absolute contribution is at most

    C sum_(j<q) L_j (R0^-2)^(q-j)=O(L_q/q).         (S6)

For completeness, split j>=q/2 and j<q/2. The first range has
L_j/L_q<=(C1/q)^(q-j), giving a geometric O(1/q) bound. On the second,
eventual monotonicity bounds all L_j by a constant times L_floor(q/2);
its ratio to L_q is at most (C2/q)^floor(q/2), which defeats any fixed
exponential factor and the q summands.

The main term in (S4) has magnitude comparable to L_q, so it cannot be
cancelled by (S6). The coefficients of f(t(T)) grow faster than every
exponential. An invertible analytic parameter change preserves analyticity,
hence f diverges.

A finite actual image would give f=p exp(psi) with polynomial p,psi,
which is analytic. Purity and safety likewise require an analytic unit.
Finally, if f'/f were rational, it would be regular at zero since f is a
formal unit. The regular analytic first-order equation for f would make
f convergent, a contradiction. This proves the rank-at-least-two claim.

## 2. A finite critical-value criterion without positivity

The selected-dominant hypothesis follows from the following explicit
algebraic and absolute-value conditions:

    every genuine finite critical point of G0 is simple;
    their critical values have pairwise distinct absolute values.   (S7)

A genuine finite critical point means a root zeta of c-zc' with
c(zeta)!=0; roots at poles of G0 are excluded. In particular zeta!=0
and its value zeta/c(zeta) is nonzero. Simple means
G0''(zeta)!=0.

To prove (S7) sufficient, the selected inverse has positive finite radius
of convergence. Positivity follows from the implicit-function theorem at
(0,0). If it were entire, then for |X|>=1 it would solve

    c(z)-z/X=0,

whose leading coefficient is fixed nonzero and whose other coefficients
are uniformly bounded. Polynomial-root bounds make every such root
uniformly bounded. The selected entire inverse would therefore be a
bounded entire function, hence constant by Liouville, contradicting
z0'(0)=c(0)!=0.

At a finite nonzero X, the same local polynomial-root bounds exclude an
infinite spatial limit. Every selected singularity is consequently a
genuine finite critical value. This can be seen directly from the
implicit-function theorem: on a neighborhood with no critical selected
root all finite roots remain isolated analytic branches. A first selected
singularity must therefore have J=0; since c!=0 there and its critical
point is simple, its local inverse is a square root.

At a boundary value, finitely many isolated neighborhoods of the bounded
polynomial roots select its limit: on the connected interior sector the
branch must stay in one neighborhood. At regular roots the implicit
function theorem applies, and at a simple critical root the local quadratic
model supplies its finite limit.
At the selected convergence radius, (S7) allows just one critical value.
It is a simple singularity of the actual selected sheet. At every other
boundary point the selected root extends regularly. Compactness of the
closed selected graph off a critical disk, followed by the local
square-root model in that disk, yields the bounded slit collar required
in Section 1. Choose the outer collar radius before the next larger
global critical-value modulus. Local boundary charts glue to the interior
by uniqueness. Interior critical values on other sheets do not obstruct
this graph. The selected sheet is already analytic throughout its open
convergence disk, including zero; no assumption about transitivity of the
other sheets' monodromy is needed.

For fixed degree, (S7) supplies a Euclidean-open dense sufficient family
on the squarefree/simple-critical locus. Here is an explicit density
justification. This locus is nonempty: c=1+z^m has simple roots and
c-zc'=1-(m-1)z^m has m simple roots disjoint from c. Its complement
is algebraic and proper. Locally label the m distinct critical points
z_i. Their critical values X_i vary according to

    delta X_i=-z_i delta c(z_i)/c(z_i)^2.

The omitted delta z_i term vanishes because G0'(z_i)=0. Polynomial
interpolation of degree at most m-1 permits arbitrary independent complex
values delta c(z_i), so the critical-value map is locally a submersion.
The finitely many real equations |X_i|=|X_j| therefore have empty interior;
their complement is open dense. Thus (S7) is a verifiable coefficient
criterion rather than a claim requiring the selected radius to be located.

## 3. An explicit positive-coefficient family with multiple roots

A second sufficient family is

    c(z)=C A(z/zeta), C zeta!=0,

where A has nonnegative real coefficients, A(0)>0, deg A>=2, and

    gcd{j>=1 : [z^j] A(z)>0}=1.                     (S8)

In particular this includes every product of at least two factors
(1+lambda_j z), lambda_j>0, with arbitrary repetitions or distinct
lambda_j, followed by the displayed complex scalings. No common-root
condition is placed on h; this includes unresolved-equality local orders
whenever the chosen c has the corresponding even-multiplicity root.

It suffices to prove the geometry for c=A. There is a unique tau>0 with

    tau A'(tau)=A(tau), rho=tau/A(tau)>0.

Indeed A-zA'=A0-sum_(j>=2)(j-1)A_j z^j decreases strictly from A0 to
minus infinity on the positive axis, with derivative -z A''(z)<0.
The real inverse runs from zero to tau on 0<=X<=rho.

The selected formal inverse has nonnegative coefficients. Its iteration
z_(n+1)=X A(z_n), z_0=0, stabilizes coefficientwise. For each real
0<=x<rho the iterates are bounded by the real inverse below tau, so
every partial coefficient sum is bounded and the selected power series
converges absolutely at x. As x rises to rho its value tends to tau.
Monotone convergence yields sum z_n rho^n=tau, and the Weierstrass
M-test gives uniform convergence on |X|<=rho. In particular

    |z0(X)|<=tau on |X|<=rho.

Inside the disk, |X A'(z0)|<rho A'(tau)=1, so J is a unit. If a boundary
point has J=0, then

    1<=rho A'(|z0|)<=rho A'(tau)=1,

forcing |z0|=tau. Equality in z0=X A(z0) and the triangle inequality
then forces every nonzero summand A_j z0^j to have the same phase as
the positive constant term. Writing z0=tau exp(i theta), one gets
exp(i j theta)=1 for all positive-support j. The gcd condition gives
theta=0, and hence z0=tau and X=rho. This is the unique selected
boundary critical point. Moreover

    G0''(tau)=-tau A''(tau)/A(tau)^2!=0.

Regular continuation across the remaining compact boundary and the local
quadratic model at rho give the required bounded slit collar. The
leading singular coefficient is explicitly positive:

    K0(X)~kappa (1-X/rho)^(-1/2),
    kappa=sqrt(A(tau)/(2 tau^2 A''(tau)))>0.          (S9)

For c(z)=C A(z/zeta), replace rho by zeta rho/C; the kernel's leading
coefficient kappa is unchanged. Section 1 therefore proves original-unit
divergence for the whole family (S8), with arbitrary h,u,k,a.

The selected-dominant geometry also persists in a complex coefficient
neighborhood of any fixed such A. Keep its compact slit graph and
normalize the continuing simple critical value. The joint implicit-function
and quadratic-critical arguments of Section 1 apply with the coefficients
of c as extra parameters. Therefore this construction supplies open
multiroot families beyond the real positive slice.

A concrete multiple-root equality example is c=(1+z)^4(1+2z), with
h=(1+z)^2 and arbitrary a,u,k. At alpha=-1 it has (m,ell)=(4,2),
while c also has the distinct root -1/2. Here

    c-zc'=(1+z)^3(1-3z-8z^2),
    tau=(sqrt(41)-3)/16>0.

The positive-family proof gives original-unit divergence in this example,
without asserting a uniform positive-volume proof at the equality face.

The positive square-root geometry is classical analytic combinatorics.
For context see Flajolet and Sedgewick, Analytic Combinatorics,
[23 October 2006 draft](https://algo.inria.fr/flajolet/Publications/book061023.pdf),
printed pp. 442-446 (physical PDF pp. 452-456), Definition VII.4,
Theorem VII.3, Lemma VII.3 and Notes 15-16. Their smooth implicit-function
schema and aperiodicity analysis contain the same type of selected
square-root and unique-dominance geometry. This is context rather than
an indispensable imported theorem: all geometry needed here was proved
above. The original Wick amplitude, moving critical-value normalization
and Gaussian-coefficient divergence argument are the additional steps
of the present statement.
## 4. The concrete bottleneck beyond these criteria

If several selected singularities have the same smallest modulus, the
baseline coefficient is a sum

    sum_i kappa_i rho_i^(-n)b_n + lower-order terms.

Its even subsequence can have genuine phase cancellation. Moreover each
rho_i moves independently after inserting arbitrary h,u,k; a single
analytic change of t fixes only one critical value. The uniform coefficient
majorant at the EXACT main exponential rate can then fail on every fixed
complex t disk. Merely knowing that K0's even part is not a polynomial,
or taking longer moment prefixes, does not repair that loss.

The present theorem avoids this concrete obstruction by unique selected
dominance, not by assuming that all finite critical values belong to the
selected sheet. It settles an open dense multiroot coefficient family
and the explicit positive aperiodic families above; it does not claim
all c, nor all equality ell=m/2 cases.

All formal algebraic exclusions also descend to characteristic-zero
coefficient fields when a complex embedding of their finite coefficient
field satisfies one of the displayed analytic hypotheses. The moduli
and positivity criteria themselves are embedding-dependent conditions.


## Verification scope

The [focused exact checker](verify_gaussian_w_affine_selected_dominant_divergence.py)
verifies the multiroot quadratic critical values, the explicit equality
example's polynomial factor and critical point, and original Wick versus
completed coefficients through order four for one polynomial perturbation.
It also checks generic completion identities and coefficient-ratio limits,
and rejects zero constant term, degree below two and period-two support
from its positive-representative test. The
[internal verification receipt](gaussian_w_affine_selected_dominant_divergence_verification.json)
records the exact proof/checker bytes and bounded scope.

These controls do not certify selected-sheet geometry for arbitrary
inputs or infer any universal conclusion from a moment prefix. The
unrestricted continuation, contour, generic-density and divergence claims
are proved above. No independent replay, external review or formal
verification is asserted.
