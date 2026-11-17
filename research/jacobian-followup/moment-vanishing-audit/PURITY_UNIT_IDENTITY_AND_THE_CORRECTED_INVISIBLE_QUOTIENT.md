# Purity, the unit identity and the corrected invisible quotient

17 November 2026. Written all-degree purity consequence. Both FULL exceptional
Mathieu--Zhao cases remain **OPEN**. Independent replay, formal verification
and external review are not claimed; review remains deferred.

The actual corrected obstruction of EVERY eventually pure phase lies in the
constant annihilator of the SPECIFIC original unit's boundary residue vector.
No eigenvalue-one semisimplicity of the fibre is needed. This extends the
previous semisimple-unit detector to the full canonical all-genus obstruction.
The proof uses the exact unit identity, a global differential quotient and
the already proved local quotient theorem; it does not choose a generic
cyclic vector or assume the unit generates the marked module.

## 1. Use the exact physical identity in the full source

Let X=SL2/Gamma, in particular either exceptional quotient. Fix a nonconstant
original regular phase f and normalized original volume Omega, with h(1)=1.
Let z denote the phase parameter, K=C(z), and put

    B=H^3_dR(X_K minus {f=z}), u=[Omega],
    a=[Omega/(z-f)].

The physical periods are period(u)=1 and period(a)=C_1(z). If f is
all-positive pure, compact convergence gives the EXACT germ

    C_1(z)=1/z near infinity.                         (1)

Its analytic continuation as a period germ is the same rational function.
This assertion uses the unbounded component/physical germ; it does not
identify independently defined Cauchy integrals on every component of the
complement of the compact image. For eventual purity replace 1/z throughout
by the actual finite Laurent polynomial r(z)=C_1(z).

Use the retained canonical quotient

    E=B/C_f, W4 E=L=K*u,
    A=W5 E, P=A/L, T=E/A.                            (2)

C_f is the globally semisimple compact-middle system with NO trivial
constituent. P is the pure genus system; T is finite Tate. These are the
actual original systems, with all branches retained. Define the SPECIFIC
unit-residual cyclic submodule

    w=bar(a)-u/z in E,
    W_unit=D_K*w,
    Q_unit=E/W_unit.                                 (3)

For eventual purity use w=bar(a)-r(z)u. The quotient is a global differential
quotient, not necessarily Hodge and not necessarily a scalar-function module.
The top Tate part of w is independent of which rational r is used.

Under purity, Q_unit RETAINS THE UNIT. To prove this without an actual-kernel
assumption, let tilde(W)=D_K*(a-u/z) in B. Its actual physical periods are
zero by (1), including all derivatives and continuations. If u belonged to
W_unit, some tilde(w) would satisfy tilde(w)=u+c with c in C_f. Its period
would give period(c)=-1. But the period image of the globally semisimple
nontrivial C_f has no trivial subquotient and cannot contain literal 1.
This contradiction proves

    W_unit intersect L=0.                            (4)

The same proof works with r(z) for eventual purity. Thus (3) is the universal
global quotient imposing bar(a)=u/z as a cohomological relation, and the
actual boundary period quotient factors through it under purity. This
does not assert that unit retention of (3) produces a pure physical period.

## 2. The genus part cannot hide a local unit in any global quotient

We need the following consequence of the completed intrinsic local proof:

    Every GLOBAL quotient A -> A' retaining L has
    a LOCAL unit retraction at every puncture.         (5)

Here A=W5 E, and no global source-unit splitting is assumed. For completeness,
apply the lifted-genus-socle construction from the
[full local criterion](FULL_INTRINSIC_LOCAL_TATE_CRITERION.md) to A. It gives
a pure weight-five nontrivial semisimple subvariation C, with Abar=A/C having
no nontrivial simple submodule. Any global quotient of A retaining the unit
can first be quotiented by the image of C without losing the unit. Its
remaining kernel Kbar injects into the pure semisimple genus quotient, so
Kbar is trivial and lies in F=Abar^global.

The mixed fixed-part theorem makes F a constant MHS, with W4 F=L. Strict
limit inclusion gives F intersect M4 Abar=L. The source A itself has a
local unit retraction: its connecting map from ker N_P, of weights at most5,
to L(-1), of weight6, is zero. Removing C reflects that local unit test by
M4 C=image N_C. If a quotient then had u=Nx+k with k in Kbar, both u,Nx are
in M4, so k belongs to F intersect M4=L. Since the kernel avoids L, k=0,
contrary to the source local retraction. Finite averaging restores the
original local eigenvalue-one sector. This proves (5), including complex
global quotient kernels. It is the same admissible argument with top
Tate grade zero, not a statement for arbitrary LOCAL quotients.

## 3. Compute the precise unit-invisible Tate space

Let

    rho(z)=projection_T bar(a),
    tau=D_K*rho subset T.                            (6)

Projection of W_unit onto T is exactly tau. Put A'=image(A->Q_unit). By
(4), A' is a global quotient of A retaining the unit, so (5) applies. There
is an exact underlying differential-module sequence

    0 -> A' -> Q_unit -> T/tau -> 0.                 (7)

The finite-monodromy quotient T/tau is the top unit-invisible space. Its
dual is described particularly simply on a finite labelled cover trivializing
T. Use a CONSTANT flat graph basis and the one finite algebraic field
L_alg/C(z) of the original residue construction. Write rho=(rho_1,...,rho_r).
Then

    K_unit={c in C^r:sum_j c_j rho_j(z)=0 in L_alg}
           =(T/tau)^dual as a flat space.            (8)

This is the actual original unit vector. A generic vector or a parameter-
dependent change of flat frame is not a replacement. Its dimension is
computed by the rank of the matrix of the first r covariant derivatives
of rho; on the labelled flat cover these are ordinary derivatives. Thus
the full Wronskian is nonzero exactly when K_unit=0. Constant dependence,
not pointwise dependence at one phase value, is intended in (8).

## 4. The CORRECTED class factors through the invisible space, with no SS gap

Fix phase zero, or any chosen puncture when a rational unit identity is
given there through (3). Let I denote the ORIGINAL finite local monodromy
on T. The completed local Tate normal form gives

    (E/L)^local-invariants -> T^I onto,
    beta_E:T^I -> L,                                 (9)

where beta_E is the canonical local connecting row and is the ENTIRE
corrected intrinsic obstruction. It kills invariant genus classes. It is
equivalent to the normalized Neron component nu_0(sigma_f) under the
retained Betti/Tate comparison. None of this assumes semisimple genus
or unreduced middle-fibre monodromy.

There is a UNIQUE row beta_blind on (T/tau)^I such that

    beta_E=beta_blind composed with [T^I -> (T/tau)^I]. (10)

To prove it, the image A' in (7) has a local unit retraction by (5).
Hence the connecting map for 0->L->Q_unit->Q_unit/L->0 kills the invariants
of A'/L. Its invariant domain maps ONTO (T/tau)^I: lift an invariant class
first to T^I, possible because I is finite, then use (9), then project
to Q_unit/L. The connecting map therefore factors uniquely through
(T/tau)^I. Naturality gives (10). The projection in (10) is onto, so

    beta_E=0 iff beta_blind=0.                       (11)

Equivalently, Q_unit retains the FULL local corrected obstruction. This
also follows from the global unit-retaining quotient theorem. The argument
above additionally identifies its exact Tate carrier in (7). Neither
Q_unit nor its quotient in (7) is asserted to embed in scalar Cauchy functions.

In the normalized toric coordinates of the actual boundary motive, (10)
proves for EVERY eventually pure phase

    nu_0(sigma_f) in K_unit intersect
                        (X_*(T_G) tensor Q)^I,       (12)

after complexification in (8) and the standard dual Tate identification.
This is the ACTUAL corrected vector, not an arbitrarily projected Haar
variation. The proof avoids a literal formula for all ordinary log
coefficients in the nonsemisimple full source; no such formula is claimed.

Consequently

    K_unit^I=0 implies full all-marker safety
               for every eventually pure f.          (13)

More strongly it suffices that K_unit contain no nonzero RATIONAL
I-invariant vector. A nonzero complex kernel alone may contain no admissible
rational component. A nonzero kernel never by itself constructs a failure.

For a practical unit-specific test, project rho to T^I using the finite
averaging projector and use its coordinates in a constant basis of T^I.
Their derivative Wronskian being nonzero is equivalent to K_unit^I=0.
For dim T^I=0 the empty determinant is1. This extends the earlier
semisimple-fibre Wronskian implication to ALL genus and ALL fibre monodromy.

## 5. A geometric support bound for what the unit can see

Choose a labelled relative SNC compactification over the regular phase
base. At an actual node let u,v be boundary coordinates and z the phase
coordinate. The unit fibre form eta_1=Res(Omega/(z-f)) has double Laurent
residue rho_e. If Omega is holomorphic along either adjacent horizontal
boundary divisor, eta_1 has no negative power of that boundary coordinate.
Therefore

    rho_e=0 whenever either incident component
                    is volume-regular.              (14)

Higher poles along the other component do not change this conclusion.
The relative SNC coordinates and generic smoothness are essential; this
is a statement on the regular generic fibre, not the entire special fibre
without its additional components. All residues still satisfy incidence
balance, as proved in the original-residue theorem.

Let Gamma_pol be the subgraph consisting of polar boundary vertices and
edges with BOTH ends polar. Include other vertices as isolated if desired.
Then in the constant graph frame

    rho(z) in H_1(Gamma_pol;C) tensor L_alg,
    tau subset H_1(Gamma_pol;C).                     (15)

Derivatives do not create coordinates on the identically zero edges. Hence
K_unit contains the annihilator of H_1(Gamma_pol), canonically identified
with the cohomology pulled back from the graph contracting every connected
polar component. In particular

    dim K_unit >= b1(Gamma_C)-b1(Gamma_pol).          (16)

For an actual component take rational ORIGINAL inertia invariants of this
subspace; a full graph-dimension bound alone does not guarantee a nonzero
admissible local direction. The graph contraction here computes an
annihilator only; it is not descent of the corrected section or a blowdown.

Every loop and parallel attachment edge in that quotient remains. If
Gamma_pol is a forest, the top unit residue vector is identically zero
and this detector leaves the whole corrected torus space invisible.
It does not follow that the actual obstruction is nonzero or that other
consequences of purity supply no constraints.

For a pure exceptional phase, the retained strict-leading/capacity theorem
places all possible finite original-volume poles above its unique dominant
left root. Thus the visible cycle subspace (15) lies in the polar part of
that region. This localizes the UNIT DETECTOR, not the corrected section.
It is not permission to contract the volume-regular components. The
balancing correction and the original volume remain in the actual class.

## 6. Calibrate against the actual examples, not generic cyclicity

There is a stronger exactness statement for an actual eigenphase. Suppose
delta is a regular vector field on X, L_delta Omega=0 and delta(f)=c*f,
c!=0. Then d(i_delta Omega)=0. Since X is affine and H^2_dR(X)=0, this
regular two-form is exact. On F_z the relative unit form is, up to the
fixed residue sign, (i_delta Omega)/(c*z), so its ENTIRE H^2 class is zero.
The full Gysin sequence puts a in the ambient line. Stokes applied to

    d(f^m i_delta Omega)=m*c*f^m Omega

gives purity. Its actual unit period fixes the ambient coefficient, hence

    a=u/z in the FULL B, rho=0, W_unit=0.             (17)

No completeness of delta is required for (17). This is a unit-exactness
statement and supplies no all-marker safety proof for arbitrary such fields.

* The exceptional equivariant Bezout controls with pairings 2,-2 have
  proper ORIGINAL-volume models over the finite phase line. On a generic
  compact fibre eta_1 extends holomorphically, so rho=0. Their entire
  corrected Tate obstruction is invisible to this unit residue detector.
  Those controls are proved NONPURE; their scalar unit transform is not
  asserted to be 1/z. This explains the stated diagnosis rather than
  producing an exceptional pure failure.
* On both known PURE SAFE weight-one exceptional biform families, (17)
  applies to the actual volume-preserving torus field. Thus rho is ZERO,
  not merely bounded in rank by1. The ENTIRE graph space is invisible;
  the weight-one action trivializes the family over C*, so original
  inertia is trivial as well. Its rational invariant invisible dimension
  equals b1(Gamma_C), which is unbounded, while nu_0=0. Thus even the
  admissible rational invisible space is nonempty on the actual exceptional
  pure locus; it is not itself a counterexample.
* For the nonexceptional pure Long phase, the full source has rank2 and
  is Kummer, with no compact middle. Its actual period map is injective
  by the rank-two unit argument. Since C_1=1/z, a=u/z in that source:
  W_unit=0 and Q_unit=E. The actual nonzero obstruction is wholly invisible
  to the unit. This is a genuine pure control outside the two target quotients.

These examples keep three separate notions apart: a nonzero invisible
space; a nonzero corrected class in it; and an actual exceptional pure
phase realizing that class. Only the last would settle the cases negatively.

## 7. Relation to the one-dimensional moment literature and the remaining task

[Pakovich--Muzychuk, PLMS 99 (2009), 633--657](https://doi.org/10.1112/plms/pdp010)
reduce their interval polynomial moment problem to specific algebraic
branch relations and a monodromy-invariant subspace; see Theorem2.2 and
Section2.2 in the primary paper. Its polynomial decomposition classification
and its tree model are NOT imported for these surface fibres. The relevant
parallel is to use the invariant span of the actual data and its annihilator.
Equations (3)--(12) supply that step here using the existing geometric
admissibility and corrected-source machinery, rather than assuming its
transfer from dimension one.

On the surviving mixed exceptional families, the missing task is now
sharper: determine the actual rational I-invariant part of (8), then show
the balanced ruling/Neron vector has zero projection there, or construct
an actual pure phase with a nonzero vector there and a fixed original
bad marker. No proof that this space is zero on those unclassified families
is supplied. The previous family-transfer, capacity, automorphism and
primitive constraints remain inputs. No new degree campaign, ordinary
resolvent algebraicity, Gaussian-EGF implication or external review is
claimed. Both FULL exceptional Mathieu--Zhao cases remain OPEN.
