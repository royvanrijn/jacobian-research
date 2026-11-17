# Canonical finite Ga trace and loss of the original unit

17 November 2026. Written geometric obstruction to a proposed detector,
uniform in phase degree. Both FULL exceptional Mathieu--Zhao cases remain
**OPEN**. External review is deferred; independent replay and formal
verification are not claimed.

Every possible exceptional failure has Ga degree at least five on 2O
and seven on 2I. On the canonical finite Ga chart of EVERY such phase,
the original unit fibre form has trace ZERO. Moreover that chart kills
the original ambient Haar line before any trace is applied. Thus this
projection supplies no unit-retaining bridge to the intrinsic splitting
calculation. The trace-free cover data and the deleted divisor must both
be retained; no localization of the corrected Neron vector follows.

## 1. A canonical finite chart for every Ga-nonconstant phase

Let A=C[SL2]^Gamma, Gamma=2O or2I, and fix

    D=z partial_x+w partial_y,
    B=ker D=C[z,w]^Gamma.

Use the full torsor and higher-plinth theorem from
[the quadratic trace model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md).
For f in A with positive Ga degree d, write

    D^(d+1) f=0,       ell=D^d f/d! !=0.

Then ell is an ORIGINAL regular element of B and belongs to its positive
degree ideal: D^d A intersect B has only binary degrees at least d.
Consequently Spec B_ell lies in the smooth punctured Klein base U_Gamma.
Define on X_ell, X=SL2/Gamma,

    v=D^(d-1) f/(d! ell),       Dv=1.                    (1)

The locally nilpotent Taylor expansion gives

    A_ell=B_ell[v],       f=P(v),       P in B_ell[v],
    deg_v P=d,           leading coefficient ell.        (2)

Indeed a slice identifies the localized algebra with the polynomial
algebra over its kernel by subtracting successive highest derivatives.
The inverse of ell is D-invariant, so the kernel is precisely B_ell.
Formula (1) also makes the v^(d-1) coefficient zero, though this is not
needed below. No arbitrary rational leading coefficient is substituted.

For a generic phase value s the map

    pi_s:F_s intersect X_ell -> Spec B_ell

is FINITE FLAT of degree d, with algebra

    B_ell[s,v]/(P(v)-s).

This follows by dividing P-s by the invertible ell to make it monic.
The statement is over the phase base, and its base change at a generic s
has the same rank. Ramification is retained. If a phase cover separates
geometric fibre components, trace still sums ALL of them.

## 2. The actual original volume and an exact trace identity

Let omega_B be the descended bottom symplectic form, with normalization
chosen so that

    i_D Omega=pi^*omega_B,
    Omega=dv wedge pi^*omega_B on X_ell.                (3)

Upstairs omega_B is a nonzero constant multiple of dz wedge dw. Its
primitive is the same multiple of (z dw-w dz)/2, invariant under right
Gamma because Gamma is contained in SL2. Thus there is a REGULAR form
lambda_B on the entire punctured base with

    d lambda_B=omega_B.                                (4)

On the unramified locus of pi_s the original marked fibre form is

    eta_Q = Q(v) omega_B/P'(v),

with convention df wedge eta_Q=Q Omega. For any ORIGINAL regular marker
Q, regard its restriction as a polynomial in v over B_ell. Let R_Q be
its Euclidean remainder modulo P(v)-s, of v-degree less than d. Then

    Tr_Ga(eta_Q)
      = [v^(d-1)]R_Q / ell * omega_B.                  (5)

This is the rational trace of the actual differential form, not the
trace of the coefficient Q alone. To prove (5), let v_1,...,v_d be the
distinct roots over a splitting field at the generic point. Lagrange
interpolation of R_Q gives

    sum_i R_Q(v_i)/P'(v_i)=[v^(d-1)]R_Q/ell.

It is an identity in the original function field, so no branch is
discarded. The right side is regular on the finite base chart: remainder
division only uses ell^-1. In particular its identity extends across
the branch values wherever a regular trace of top forms is used. No
trace operation through the DELETED ell=0 divisor is asserted.

For the original UNIT, (5) gives the exact result

    Tr_Ga(eta_1)=0 if d>=2;
    Tr_Ga(eta_1)=omega_B/ell if d=1.                    (6)

The zero in (6) is cancellation among the actual Ga branches. It does
NOT say eta_1 is exact upstairs or that its intrinsic boundary residue
vector rho vanishes. The finite Ga trace and the specific unit-annihilator
K_Q of the intrinsic boundary theorem are different constructions.

## 3. A nonzero marked trace can still be entirely exact

The first obvious regular marker with nonzero trace is Q=Df=P'(v).
Equation (5) gives

    Tr_Ga(eta_(Df))=d omega_B.                          (7)

But on the ENTIRE smooth original generic fibre, including its ell=0
part,

    eta_(Df)=(i_D Omega)|F_s=d(lambda_B|F_s).            (8)

The contraction identity df wedge i_D Omega=(Df)Omega proves (8);
the base primitive in (4) is global. Thus this marked class is exactly
zero in fibre cohomology, not a detector of the corrected obstruction.
Independently, Haar invariance under D gives for every m>=0

    h((Df) f^m)=h(D(f^(m+1)))/(m+1)=0.                 (9)

These are moments with ONE fixed original regular marker, and no purity
premise. The ordinary marked Cauchy germ is identically zero. Therefore
a nonzero rational trace such as (7) cannot be promoted to an actual
failure witness or a nonzero corrected cycle pairing.

## 4. The finite chart loses the original ambient unit

Equations (3)--(4) also give the ORIGINAL primitive on X_ell

    Omega=d(v i_D Omega).                             (10)

The form in parentheses is regular on X_ell. Consequently

    H^3_dR(X) -> H^3_dR(X_ell) sends [Omega] to zero.   (11)

The same is true of the original unit in the complement source
H^3_dR(X minus {f=s}) after restriction to its ell!=0 chart. It is
already killed by restriction, before the finite trace or any graph
contraction is considered. In particular this is NOT a global
unit-retaining source quotient of the sort required by the intrinsic
local reflection theorem.

The localization exact sequence does prove the ambient support statement

    [Omega] lies in image[H^3_{ell=0}(X)->H^3(X)].      (12)

Formula (10) supplies its transgression on the complement of that support.
This is support of the ambient unit, not support of the balanced ruling
line bundle or corrected Neron component. No implication between those
supports has been proved here.

Every component over ell=0 must be retained in a proper phase model,
including degree-drop, vertical and additional resolution components.
The Ga-leading divisor ell=0 also must not be identified with the unique
dominant LEFT leading root or with an ambient resolution tree. Those
are separate data. The earlier fixed-divisor and raw-ruling calculations
retain their exact scopes.

## 4a. A relative ambient interface that retains the unit and every marker

The lost unit can be retained in an explicit relative model. For ANY
original R in A, let H_R(v) be its polynomial Ga antiderivative over
B_ell, normalized by H_R(0)=0. Define

    beta_R=H_R(v) omega_B,       d beta_R=R Omega.

The equality follows because all differentials in the two-dimensional
base wedge to zero with omega_B. Every beta_R is regular on the SAME
X_ell. In the de Rham cone for the pair (X,X_ell), the closed cocycle

    (R Omega,beta_R) in H^3_dR(X,X_ell)

maps to [R Omega]=h(R)u in H^3_dR(X). In particular the R=1 cocycle
maps to the nonzero original Haar line. This is a unit-retaining RELATIVE
AMBIENT interface, not an identification with the full intrinsic
phase-complement extension. It keeps all localization data on ell=0;
replacing it by its X_ell entry alone again loses the unit.

This construction also retains the ORDINARY marked generating series.
For ONE fixed original Q put, formally near s=infinity,

    H_Q(s,v)=sum_(m>=0) H_(Q f^m)(v)/s^(m+1),
    partial_v H_Q(s,v)=Q(v)/(s-P(v)), H_Q(s,0)=0.        (12a)

Let v_i be ALL the roots of P(v)-s in a finite splitting field over
Frac(B)(s). Divide Q(v)/(s-P(v)) into a polynomial W_Q(v;s) plus partial
fractions. At the generic point (12a) is exactly

    integral_0^v W_Q(t;s)dt
      -sum_i Q(v_i)/P'(v_i) log(1-v/v_i).              (12b)

The logarithms are chosen by their primitive vanishing at v=0. Differentiation
and that normalization prove the formula. For Q=1 the polynomial part
is zero. The trace identity (6) says only that the SUM of its logarithmic
coefficients is zero for d>=2: the common log(v) at Ga infinity cancels.
The individual root logarithms need not vanish or become algebraic.

The field here has the TWO-dimensional Klein base as well as s. Formula
(12b) is a Ga primitive before the remaining relative localization or
Haar pairing. It is not a formula asserting that the actual ordinary
Haar resolvent C_Q(s) is algebraic or is a single logarithm. No exponential
generating function is used. All branches and the original fixed Q remain.

There is an actual exceptional calibration. For the minimal invariant
bottom form F of degree n=8 or12, put

    f=a_F=(x F_z+y F_w)/n,       Df=F.

It has the strict left weight 2-n, so it is pure SAFE. Its Ga degree is
one, v=f/F, and its ordinary Ga unit primitive is exactly

    -(1/F) log(1-f/s).

As a function of the Ga coordinate this has a nonzero logarithmic residue
and is not algebraic. Yet its ACTUAL normalized unit Haar resolvent is
C_1(s)=1/s. The remaining relative Haar pairing removes these logarithms.
Thus even on the original exceptional pure locus, algebraicity of a Ga
primitive and algebraicity of the actual ordinary resolvent differ.

The missing multiplicative step is visible already in the exact identity

    H_(QR)(v)=Q(v)H_R(v)-integral_0^v (DQ)(t)H_R(t)dt.  (12c)

Integration by parts proves (12c), including its zero value at v=0.
Thus unit moment vanishing in the relative ambient model does not
automatically give marked vanishing: multiplying its primitive produces
the displayed correction. Computing the corrected boundary coupling
requires this correction and the ell=0 data, not just the traceless
branch sum. Neither vanishing of that coupling nor an actual exceptional
pure failure is established by (12a)--(12c).

## 5. All possible exceptional failures fall within the traceless range

This elementary filtration bound makes (6) relevant to the entire
remaining problem. The smallest nontrivial Peter--Weyl left highest
weight is n0=8 on2O and n0=12 on2I. In a left module of highest weight
n, D-degree at most d means only ladder positions j<=d occur. Their
left torus weights are

    2j-n.

If d<=n0/2, every nonconstant component has nonpositive torus weight;
constants have weight zero. The zero-weight part f0 has precisely the
same pure moments as f, since negative weights cannot contribute to a
weight-zero moment. The completed torus theorem gives

    eventual purity of f implies f0=0.

Thus f has only strictly negative weights and is safe for ALL original
markers, uniformly when 2m>deg_matrix Q. This proves that any genuine
failure must have

    d>=5 on2O,       d>=7 on2I.                        (13)

The argument applies after EVERY algebraic automorphism and every left
conjugation of D, because normalized Haar, purity and safety are preserved.
These are Ga degrees, not full matrix degrees or support counts. This
safe sector recovers part of the already closed one-sided/nullcone
families; no new mixed safe normal form is claimed from this bound.

The threshold for merely encountering mixed PURE phases cannot be raised
by a nonpurity argument. In the minimal left module take the binary vector

    Y^(n/2-1) (X+Y)^(n/2+1),       n=8 or12.

Its top-row ladder degree is n/2+1, its weights range from -n to2, and
its root at X+Y=0 has multiplicity n/2+1>n/2. It is in the left nullcone,
so its ORIGINAL invariant matrix-coefficient phase is pure SAFE. This
is a left interchange of the retained sharp sparse controls, not a
finite moment test or an actual failure at the bound (13).

## 6. Consequence for the full geometric bridge

Every candidate failure therefore enters a canonical finite projection
where its unit form is traceless and its ambient unit has been removed.
This proves a uniform obstruction to detecting the corrected class by
that projected unit alone. Neither a nonzero marked trace nor an ambient
support statement repairs the missing unit-retaining comparison.

The next geometric calculation must keep the trace-free cover part AND
the ell=0 localization data, with the original volume and all components,
then compute the actual corrected cycle coupling. We have not shown
that this coupling vanishes under purity or constructed a pure exceptional
phase realizing a nonzero coupling. The full unit-residue monodromy
criterion remains available on the actual intrinsic boundary; (6) does
not imply its residue vector is zero.

The capacity, automorphism, primitive-fibre and family-transfer constraints
remain inputs. No finite-degree campaign, extra decidability milestone,
ordinary/EGF identification or external review is introduced. Both FULL
exceptional Mathieu--Zhao cases remain OPEN.
