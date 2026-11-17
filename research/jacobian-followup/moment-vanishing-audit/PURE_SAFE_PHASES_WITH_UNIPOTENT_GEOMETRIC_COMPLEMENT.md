# Pure safe phases with unbounded unipotent geometric complement monodromy

16 November 2026. Internal written existence/construction theorem with a
named classical semistable spectral-sequence input. Independent replay,
formal verification and external review are not claimed; review remains
deferred. Both FULL exceptional Mathieu--Zhao cases remain **OPEN**.

On BOTH exceptional quotients there are original regular invariant phases
with just TWO POSITIVE left weights, pure and safe in every marked order,
yet with arbitrarily large nontrivial unipotent eigenvalue-one blocks in
the WHOLE geometric middle-fibre and complement cohomology at phase value0.
Every actual original marked Cauchy resolvent is nevertheless a finite
Laurent polynomial. This distinguishes the full geometric source from its
selected actual Haar-period image and from the intrinsic SOURCE-UNIT test.

The retained sufficient condition of eigenvalue-one semisimplicity on the
whole geometric complement remains sufficient. It cannot be a consequence
of purity on the ENTIRE pure locus. These already-safe examples do NOT
refute such a condition explicitly restricted to a further unresolved
candidate locus. No actual pure MZ failure is constructed.

## 1. The all-degree biform input and the new two-weight phase

Retain the completed [invariant biform construction](UNBOUNDED_SAFE_BOUNDARY_HOMOLOGY_ON_EXCEPTIONAL_QUOTIENTS.md).
Let G=Gamma/{+-I}, and put

    Gamma     |G|=N     exponent e
    2O         24          12
    2I         60          30.

For EVERY k divisible by e with k>=2N, a nonempty open set in the
invariant biform space

    V=H0(P1 times P1,O(k+1,k-1))^G

defines smooth irreducible curves, with the exact quotient genus

    g=k(k-2)/N>0.                                    (1)

The linearization is the retained binary-group one. Its compatible first
jets, nonfree-point behavior and genus calculation are in the input proof.
We do not assume invariant sections are basepoint-free on the projective
surface; they are NOT, and the pencil resolution below retains this fact.

Choose two independent invariant sections a,b giving a generic pencil.
A general pencil member is in the stated smooth open set. Their ratio
rho=a/b is a rational function on the projective quotient

    Z=(P1 times P1)/G.

This surface is rational by the retained ruling/diagonal-section argument.
Resolve Z and the base ideal of the pencil to obtain a smooth projective
RATIONAL surface W and a morphism

    p:W->P1, rho as its rational affine coordinate.   (2)

The general fibre is the smooth curve C_H/G and has genus g. Base-point
blowups and resolution do not change its function field or smooth genus.
A general member is irreducible/reduced; hence the morphism has connected
smooth general fibre. After replacing a by a-rho0*b for a regular value
rho0, we may and do assume

    C=p^-1(0) is smooth, connected, genus g.           (3)

This replacement stays within the same invariant biform space. Regular
values and the retained smooth-section open have nonempty intersection.
This handles quotient fixed points and infinitely near pencil base points;
it is not a claim that a,b were basepoint-free initially.

Now restrict a,b to det=1 on the ORIGINAL X=SL2/Gamma and define

    f=a+b^2.                                         (4)

Both a,b have original left weight2. Thus f has weights2 and4 and is a
regular invariant phase. Its filtered matrix degree is4k: b is not
divisible by det, so its top square cannot lose degree modulo det-1.
For every original regular marker q,

    h(f^m)=0, m>=1,
    h(q f^m)=0 whenever 2m>deg_matrix(q).              (5)

Every term of f^m has weight at least2m, and the marker's weights lie
between minus and plus its degree. Therefore (5) is an exact all-order
weight proof of PURE SAFETY, not a moment-prefix inference.

## 2. The original generic fibre has an exact double-cover field

The effective left torus T/{+-I} is G_m; b has effective weight ONE.
Its rational quotient field is C(Z), since removing the diagonal does
not change that field. The retained weight-one torsor trivialization gives

    C(X)=C(Z)(v), v=b,
    a=rho*v, f=v^2+rho*v.                            (6)

Consequently the generic original fibre f=s has the function field

    C(W)(v)/(v^2+rho*v-s).

Put w=2v+rho. This gives the EXACT birational model

    w^2=rho^2+4s.                                    (7)

All of this is on the original quotient: no nonregular phase is inserted
and no change of period functional is made. The rational coordinates are
used only to describe a projective model of its actual generic fibre.

Equation (7) is geometrically connected for s!=0 sufficiently near0:
the two branch values in rho are distinct regular values of p, and their
fibres have multiplicity1. The quadratic function is therefore not a
square even after algebraic closure of the parameter field. In particular
f is primitive on the original quotient. Primitivity is not needed for
purity or for the countercontrol below.

## 3. A proper semistable model with exactly two rational components

Let L=p^*O_P1(1) on W and let A,B be the pullbacks of the two base
coordinate sections, with A/B=rho. They are BASEPOINT-FREE on W.
They need not equal the original sections before resolving their pencil.
Over a small disk in s, form the finite double cover of W times the disk

    Y_s: omega^2=A^2+4s B^2, omega in L.              (8)

It is projective over the disk and has (7)'s function field. For s!=0
small its branch divisor consists of two disjoint smooth fibres of p,
so Y_s is a smooth projective surface. At s=0,

    Y_0=W_+ union_C W_-, W_+=W_-=W,
    omega=+A or omega=-A.                             (9)

Along C, B is nonzero by basepoint-freeness. With local coordinates
u=omega/B and r=A/B, the family equation is

    (u-r)(u+r)=4s.                                   (10)

The total space is smooth and the central fibre is REDUCED simple normal
crossing with TWO components and NO triple intersections. Away from C
one of the omega derivatives is already nonzero. Thus no parameter cover
or hidden multiplicity is needed for this projective semistable family.
Original local monodromy on H2(Y_s) is unipotent.

## 4. The genus forces rank 2g monodromy in proper H2

Use the classical semistable weight spectral sequence, with its monodromy
weight identification and E2 degeneration. A precise primary statement
and the E1 differential description were read in FULL in Sukjoo Lee,
[Section2.2,Theorems2.6-2.7](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/semistable-degenerations-of-calabiyau-manifolds-and-mirror-pw-conjectures/9A7E93052B2CB9C09A02C48DACF2620D),
Forum Math.Sigma12(2024),e94. This is the classical background input,
credited there to Peters--Steenbrink; no mirror-symmetry or Calabi--Yau
hypothesis/result is imported. Our family is proper projective semistable,
which is more restrictive than needed for that background setup.

For two components, the only columns are p=-1,0,1:

    E1^(-1,q)=H^(q-2)(C)(-1),
    E1^(0,q)=H^q(W_+) direct_sum H^q(W_-),
    E1^(1,q)=H^q(C).                                 (11)

The differential at H1(C) is the alternating restriction from H1 of the
components. Both components are rational, so H1(W_+)=H1(W_-)=0. Therefore

    Gr^M_1 H2_lim(Y_s)=H1(C), dim=2g,
    Gr^M_3 H2_lim(Y_s)=H1(C)(-1).                    (12)

Here H2 and H1 denote COHOMOLOGY, not homology. M is the monodromy weight
filtration centred at2. Its nilpotent map N_Y identifies grade3 with
grade1 after the usual twist. There are no grades0 or4 because there are
no triple intersections. Hence

    N_Y^2=0, rank N_Y=2g.                             (13)

This is a property of the actual projective family (8), not an arbitrary
abstract variation with a prescribed nilpotent matrix.

## 5. Why the unipotent blocks survive on the original affine fibre

Two standard geometric facts are needed, with their hypotheses retained.

First, birational smooth projective surfaces have H2 differing only by
Tate classes from point blowups. For families, take a common resolution
over C(s) and spread out over a common punctured algebraic parameter base.
Closed blowup centres can be geometrically labelled after a finite cover;
their exceptional classes have FINITE permutation monodromy. Pullback and
the orthogonal polarization decomposition therefore preserve every
nontrivial unipotent Jordan block at original eigenvalue1. This can be
applied over the original parameter field: no eigenvalue is silently
changed to1 by taking a cover. Thus ANY smooth projective compactification
Z_s of the original affine fibre F_s={f=s} has rank2g nilpotent part in
its generalized ORIGINAL eigenvalue1 subspace of H2(Z_s).

Second, for the geometric restriction

    H2(Z_s)->H2(F_s),

the kernel is generated by the classes of boundary divisors. It is an
image of H0 of their irreducible components, with twist(-1), so is a PURE
TATE subvariation with finite monodromy. A polarization of the proper
pure H2 variation restricts nondegenerately to this Hodge subvariation;
its orthogonal complement is flat and gives a direct sum. Use a genuine
polarization of the FULL pure variation, via a relative ample Lefschetz
decomposition, rather than assuming the raw intersection pairing is
nondegenerate on boundary classes. The kernel therefore cannot remove
the nontrivial unipotent blocks of (13).
Consequently, on the actual ORIGINAL affine fibre,

    rank N_(F,1) >= 2g.                              (14)

The argument concerns Hodge/geometric restriction kernels and proper
polarizations. It does not assert the same reflection for arbitrary
non-Hodge actual-period quotients, which can and do kill these blocks.
There is no claim that every compact cohomology class injects into the
affine fibre, or that its full monodromy is unipotent on other weights.

## 6. The whole complement is nonsemisimple while actual periods are rational

The geometric Gysin sequence and H4(X)=0 give a monodromy-equivariant
surjection

    H3(U_s)->H2(F_s)(-1), U_s=X minus F_s.             (15)

The ambient unit is still retained as in the
[full Haar weight/source theorem](HODGE_WEIGHTS_OF_HAAR_VARIATION_AND_BOUNDARY_GRAPHS.md).
Generalized eigenvalue1 decomposition and its unipotent logarithm commute
with this surjection. Equation (14) therefore implies

    rank N_(U,1) >= 2g=2k(k-2)/N.                    (16)

The rank grows WITHOUT BOUND on BOTH original exceptional quotients as k
runs through the allowed multiples. The cohomology in (16) is the WHOLE
geometric complement H3, with actual phase f, not the minimal connection
of selected Haar functions.

In contrast, the SAME f has ALL its actual ordinary marked periods

    C_q(s)=sum h(q f^m)s^(-m-1)

finite Laurent polynomials by (5). Their analytic continuations have
TRIVIAL monodromy. In particular they are algebraic, their integer-log
coefficients vanish, and the intrinsic reduced boundary SOURCE UNIT
splits by the retained safety theorem. This statement does not identify
the entire reduced boundary source with its unit, or claim the full
geometric source has semisimple monodromy.

Thus the source can have arbitrarily large unipotent blocks that do not
obstruct the actual Haar marker packet. A proof cannot replace control of
the source UNIT extension by semisimplicity of every source class.

## 7. Exact consequence for the whole-quotient proof strategy

The recorded condition that the WHOLE complement monodromy be semisimple
at eigenvalue1 remains a correct SUFFICIENT condition for a pure phase.
This construction proves it is NOT NECESSARY and is FALSE as a universal
consequence of purity on the ENTIRE exceptional pure locus, even among
primitive two-positive-weight phases.

The phases here are ALREADY SAFE and in the original left nullcone. The
construction does NOT refute that sufficient condition when explicitly
restricted to some further unresolved candidate locus, and does not claim
new closure of such a locus. No result forces every original pure phase
into these controls or into the automorphism saturation of a closed class.

The viable full target remains the ACTUAL all-original-marker ordinary
resolvents, their integer logarithms at phase0, or the intrinsic reduced
source-unit splitting condition. A proof of the stronger source condition
on ALL remaining candidate phases would also suffice but is UNPROVED.
The [ordinary/exponential normalization audit](ORDINARY_RESOLVENTS_AND_HOMOGENEOUS_GAUSSIAN_CONVERGENCE.md)
remains in force: Gaussian exponential algebraicity is a different object.

This is an all-degree countercontrol to a blanket hypothesis, NOT an MZ
counterexample. Both full exceptional MZ cases remain OPEN. All capacity,
automorphism, primitive, family-transfer results/goals and the fixed
original bad-marker requirement are preserved. External review remains
deferred; no finite-degree/marked-prefix campaign, new general decision
procedure or external review/outreach is introduced.
