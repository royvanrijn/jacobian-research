# Semistable log periods and all-degree single-isotypic safety

16 November 2026. Internal written proofs with precise classical imports
and selected exact controls. Independent specialist replay, external
review and formal verification remain pending. No novelty is claimed.
See the [checker](verify_semistable_isotypic_safety.py),
[receipt](semistable_isotypic_verification.json) and
[source audit](semistable_isotypic_source_audit.json).

On BOTH exceptional quotients, EVERY single Peter--Weyl isotypic
component, with all its right-invariant multiplicities, now has pure
locus exactly the left-SL2 nullcone, in EVERY degree. Every such pure
phase is marked-safe on the full original algebra, with the sufficient
uniform field cutoff 2m>deg_matrix(Q). This includes both full V12
modules; their earlier sparse-support restrictions are unnecessary
for the field pure-locus/safety conclusion.

A full failing phase must therefore mix at least two different
nonconstant Peter--Weyl degrees. Its leading rank-one boundary form
must have a common LEFT factor of multiplicity at least half its
degree. Real seeds in one isotypic component admit no nonzero
polynomial pure phase at all. Consequently the full thin-phase
degree gates improve to36 on2O and60 on2I, and their fixed bad-marker
bounds improve to D-12 and D-20. A general2I failure has degree>=20.

These are all-degree closures and necessary conditions. Both FULL
exceptional MZ decisions remain **OPEN** on mixed-isotypic phases.
An arbitrary matrix-homogeneous representative is not a single
Peter--Weyl component: determinant powers can hide lower components.

## 1. Projective holomorphic-volume periods in arbitrary dimension

Let Y be smooth of complex dimension N>=2, f:Y->A^1 a nonconstant
PROJECTIVE morphism, Omega a holomorphic algebraic N-form on Y,
and M a compact oriented real N-cycle with integral_M Omega!=0.
No isolated-critical-point assumption is needed. For every
nonconstant p in C[t], the moments integral_M p(f)^m Omega are
not eventually zero.

Here is the extension of the preceding
[surface argument](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md).
At each finite bad base value s0, apply characteristic-zero
semistable reduction over a curve. The exact source is
[KKMSD, ChapterII, theorem on printed53--54](https://page-one.springer.com/pdf/preview/10.1007/BFb0070320),
whose complete statement is in this public two-page preview.
After a finite base change, locally s-s0=t^e, and a projective
modification, the resulting W is smooth and its central fibre
is a REDUCED simple-normal-crossing divisor. It dominates the
original base change, so Omega_W is the actual regular pullback
of Omega. Projectivity is preserved. Generic components and the
finite Stein factor can be split first and treated separately.
The construction is algebraic and its local analytic parameters
are convergent, not arbitrary formal series.

In a semistable chart t=x_1...x_r, write

    Omega_W=J dx_1 wedge ... wedge dx_N.

Its relative form Omega_W/dt is, up to sign,

    J (dx_2/x_2) wedge ... wedge (dx_r/x_r)
                   wedge dx_(r+1) wedge ... wedge dx_N.   (1)

It is a holomorphic section of the top RELATIVE LOGARITHMIC
de Rham sheaf, not a meromorphic lattice with arbitrary powers
of t inserted. The relevant classical input is
[Steenbrink, Theorem2.7 and Section2 setup](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/steenbrink2.pdf),
printed533--535: relative log de Rham cohomology is locally free,
with a logarithmic Gauss--Manin connection of nilpotent residue
after the reduced base change. His Section2 treats projective
families in arbitrary dimension, independently of the paper's
later isolated-singularity application. The setup and theorem
were checked on the original scanned pages.

The class of(1) is a holomorphic section of that log-cohomology
bundle. Nilpotent residue makes its flat periods polynomially
bounded in log(t) on a fixed sector. Explicitly a logarithmic
connection with nilpotent residue R has a holomorphic invertible
gauge times exp(R log t); the latter is a finite polynomial in
log t. Thus all these continued periods are O((1+|log t|)^K)
for some finite K. Returning to the ORIGINAL relative form
gives the factor 1/(e*t^(e-1)), hence a bound

    O(|s-s0|^(-1+1/e)(1+|log(s-s0)|)^K),                  (2)

strictly integrable at every finite value, for every branch.

We retain the compact-cycle part of the argument as well.
Cut the regular base into contractible sectors and arcs.
Closedness of Omega and its zero restriction as an N-form to
each complex (N-1)-fibre express the nonvertical part of a
compact period through the one-form part of the Gauss--Manin
complex, paired with flat middle fibre cycles. In a finite
semistable cover, logarithmic de Rham comparison supplies the
same description across the central strata. In topological
terms subdivide the compact chain into the finite plumbing
cells over the base arcs and vertical cells. The vertical
pieces have zero N-form integral. The remaining Cauchy terms
are finite sums of integrals

    integral_arc g(t)dt/(s-f_base(t))

and the residue additions from tubes, with g from the above
log-relative cohomology. This does not assume that a higher-
dimensional Milnor fibre is a bouquet of middle spheres.
Relative boundary contributions on regular pieces are kept
fixed and are holomorphic in s. The finite base change is
returned by chain trace, summing generic sheets with their
multiplicities and dividing by the degree; a transverse full
pullback of a compact cycle has no new boundary. Projective
modifications preserve these integrals of the original form.

The local Cauchy estimate already proved in the preceding
notes, together with(2), makes every finite continuation
singularity of A_f(s)=integral_M Omega/(s-f) subpolar. The
polynomial root-trace identity preserves this property for
p composed f, with exponents delta/k-1>-1. A rational germ
therefore has no finite poles and tends0 in its original
infinity germ; it is zero. Eventual vanishing would make it
a polynomial in1/s with nonzero zeroth-mass coefficient, a
contradiction.

The semistable and log-cohomology theorems are classical
imports. The compact-cycle gluing, finite trace and this
moment application are INTERNAL WRITTEN arguments. They
are not certified by the algebraic controls below. This
section asserts PROJECTIVE maps; it does not assert that
every affine map automatically has such a volume model.

## 2. Adjunction improves the boundary threshold to (N-1)/D

Retain barQ_N, its smooth boundary B, and F=H/u^D of
minimal degree D. Let Z be the zero divisor of the leading
section H|B. The sufficient condition is now

    (B, ((N-1)/D)*Z) is klt,
    equivalently lct_B(Z)>(N-1)/D.                        (3)

Here klt means that all divisorial log discrepancies of
the pair are strictly positive; coefficients along its
own divisor components must be<1 as well. It is not merely
smoothness or reducedness of Z.

Set c=(N-1)/D. By the klt/plt form of inversion of
adjunction,

    (barQ_N, B+c*div(H)) is plt near B.                    (4)

This is exactly the smooth Cartier case of
[Kollar, Theorem7.5.1](https://arxiv.org/pdf/alg-geom/9601026v2),
printed/PDF43. The coefficients are rational and H does
not contain B, so its hypotheses apply. A finite-boundary
divisor E over the graph has a=ord_E u>0, E!=B, and
ord_E H>=D*a. The plt inequality is

    A(E)>a+c*ord_E H>=N*a.

Consequently the ORIGINAL volume order

    A(E)-1-N*a>=0                                        (5)

by integrality. Resolve the projective graph, as an
isomorphism on the affine quadric, and restrict over A^1.
Every finite boundary volume pole is removed, including
vertical exceptional divisors. The map is projective and
the volume is holomorphic on the smooth resolution.
Section1 excludes every nonconstant polynomial in F.

This improves the previous sufficient ambient threshold
N/D: the coefficient1 boundary in(4) contributes the
missing a. It uses the full plt inequality, not just
the unweighted multiplicity lower bound. It is sufficient,
not necessary, for an actual holomorphic volume model.
Non-klt leading data do not certify purity or failure.

For SL2 the dimension is3 and B=P1_left times P1_right,
so EVERY nonzero eventually pure original exceptional
phase of degree D must satisfy

    lct_B(Z)<=2/D.                                       (6)

There is no reality or thin-image assumption in(6).

## 3. At the critical coefficient, only common left factors matter

Let s(u,v) be an exceptional invariant leading section
of O(D,D) on B. Expand in ANY basis f_i(v) of the right
binary invariants:

    s=sum A_i(u)f_i(v),  A=gcd_i A_i,
    s=A(u)*G(u,v).                                       (7)

A is homogeneous and G has NO left fibre as a component.
Let k_max be the maximum root multiplicity of A, or0
if A is constant. This is independent of the chosen
right invariant basis. A common root defines a fibre
{u=u0} times P1_right. These are called left fibres here.

For every such fibre, G restricts to a NONZERO invariant
right binary form of degree D. Right projective orbits
have size>=6 for O and>=12 for I. Thus every root of
each restriction has multiplicity at most D/r, where
r=6 or12 respectively. Put c=2/D. On each such fibre
the restricted pair with coefficient c is klt, because
c*m<=2/r<1.

Off the vertical components of A, restriction to that
smooth fibre and inversion of adjunction show local klt
for (B,cZ). At a component L of A of multiplicity k,
the curve pair (L,c*G|L) is klt. Inversion of adjunction
now gives (B,L+c*G) plt near L. If c*k<1, decreasing
the L coefficient from1 to c*k makes the pair klt.
The different left fibres are disjoint, so this covers
all points. Conversely c*k>=1 already prevents klt
at the generic point of L. Therefore the exact test
at this particular coefficient is

    (B,(2/D)Z) klt  iff  k_max<D/2.                       (8)

No general log-canonical-threshold formula is needed.
Combining(6) and(8), every nonzero pure phase, in EVERY
original degree, has

    a common left factor with root multiplicity>=D/2.     (9)

It follows that leading singularities without such a
half-degree fibre cannot support a pure phase, however
singular their other curve components may be. This is
a stronger condition than merely singular leading divisor.
It still does not imply that a phase satisfying(9) is pure
or failing. Algebraic automorphisms preserve Haar, so the
same necessary test applies anew to every automorphic
image, using that image's own degree and leading class.

## 4. Every single Peter--Weyl isotypic component is safe

Let

    W_D=V_D tensor (V_D^*)^Gamma subset A_Gamma.           (10)

This is the ENTIRE degree-D left isotypic component,
including all right-invariant multiplicities. The leading
map W_D->H0(B,O(D,D))^Gamma is an equivariant isomorphism.
One can see this from the Peter--Weyl decomposition and
the associated graded determinant ring. Thus(7) is also
the tuple of left binary forms representing P in W_D.

Consider P in C+W_D with all positive moments zero.
Its constant term is0 by its first moment. If P!=0,
(9) gives k_max>=D/2. Move that common root by a left
SL2 transformation to the coordinate origin. All the
remaining left weights then have one sign, except
possibly weight0. If k_max=D/2, the corresponding
cocharacter has a nonzero limit P0 consisting exactly
of the zero-weight terms. At least one such coefficient
is nonzero by the definition of the gcd multiplicity.
The limit is an ORIGINAL torus-invariant regular function.
Haar moments are invariant polynomials on the finite
module, so purity passes to P0. This contradicts
[the complete torus-invariant closure](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md).
This limit stays inside the actual left-SL2 action on
the fixed component. No matrix scaling, special-fibre
augmentation or assumed Haar-adic continuity is used.

Hence

    Nonzero P pure in C+W_D iff its constant is0 and
    the tuple in(7) has a common root of multiplicity>D/2.
                                                               (11)

The reverse implication follows immediately from its
strict one-sided weights. This is equivalently the
left-SL2 nullcone, by the retained one-parameter criterion.
There is no bound on D or on dim (V_D^*)^Gamma.

With k>D/2 the minimum absolute weight is delta=2k-D.
For any ORIGINAL marker Q,

    m*delta>deg_matrix(Q) => h_Gamma(QP^m)=0.              (12)

The binary centre forces D even, so delta>=2. Therefore
2m>deg_matrix(Q) is a sufficient uniform cutoff across
ALL these isotypic degrees. Left conjugation preserves
the original degree filtration. This closes both full
V12 modules over C, as well as all later single isotypic
components, not only their sparse sectors.

For an exact algebraic input known to be in C+W_D,
this is a finite decision procedure without moment
prefixes: check its constant, substitute a rank-one
matrix in its leading form, take the gcd of the left
coefficient forms, and compute its greatest projective
root multiplicity, including infinity. Condition(11)
decides purity, and(12) supplies its safety certificate.
Verify single-isotypic membership first; the same gcd
on the leading part of a mixed input is only necessary.

The ordinary pure-moment ideal on a fixed W_D need not
equal the ideal generated by all positive invariants.
What is now proved is equality of their RADICALS, both
defining the nullcone. Thus the previous V12 invariant-
quartic deficit remains a valid ordinary-ideal result;
its previously unresolved radical question is closed.
No numerical sufficient pure prefix for V12 or general
D is computed. Noetherianity gives finite-prefix
existence, not a value or a new implemented campaign.

The same uniform marked cutoff holds over ANY REDUCED
commutative C-algebra. At every prime, apply the field
statement to the finite set of phase and marker
coefficients using the fixed Qbar form, as in the prior
base-change proof. Each late mixed moment vanishes in
every residue field, hence is0 in the reduced ring.
Ordinary nilpotent thickenings are not settled merely
by replacing the pure ideal with its radical. The
previous all-ring torus-invariant result is retained;
the new reduced-ring statement is not silently copied
to arbitrary nonreduced isotypic families.

Since the first two nonconstant degrees on2I are12
and20, every full function of matrix degree<=19 lies
in C+W12 and is safe when pure. Thus any FULL2I MZ
failure has matrix degree>=20. On2O the first mixed
space C+W8+W12 still starts in degree12; its mixed
problem is not closed. More generally a full failure
must involve at least TWO DISTINCT nonconstant
Peter--Weyl degrees. An arbitrary homogeneous matrix
representative is insufficient to identify(10): det^j
can homogenize a mixture without removing its lower
isotypic components.

## 5. Real single-isotypic seeds exclude every polynomial phase

Let h=c+H be a nonconstant regular real-valued seed
on the compact core, with H in a single W_D. The
dagger real structure preserves this decomposition,
so c is real and H is self-adjoint. On the rank-one
boundary, dagger acts by

    (u,v)->(j(u),j(v)),  j(u0,u1)=(-conj(u1),conj(u0)).

The projective involution j has NO fixed point:
j(u)=lambda*u would imply -u=|lambda|^2*u. Thus the
vertical roots of the leading form occur in distinct
antipodal pairs with equal multiplicity.

If k_max<D/2, Sections1--3 supply a proper model with
holomorphic original volume for h. Every nonconstant
complex polynomial p(h) is excluded from eventual
purity. If k_max>=D/2, antipodal pairing and total
left degree D force exactly two roots, each of
multiplicity D/2, exhausting the whole left degree.
After a left SU2 rotation the left tuple is

    (u0*u1)^(D/2) tensor v_right.

Its leading class is weight0. Since H is in ONE
isotypic component and the leading map is equivariant
and injective, H itself is torus-fixed. Then h and
every p(h) are torus-fixed, and the complete restricted
closure excludes nonzero pure p(h) again. The case
k_max>D/2 is impossible for a self-adjoint seed.

Therefore every real seed in C+W_D has NO nonzero
polynomial pure phase, in all degrees. This argument
uses reality; it does not infer that a general complex
single-isotypic seed has the same antipodal pairing.

For an arbitrary nonzero thin pure phase P=p(h), take
the normalized regular real seed supplied by the
curve-normalization theorem. The seed must now mix
at least two nonconstant isotypic degrees. Its minimum
original degree is consequently>=12 on2O or>=20 on2I.
The prior general polynomial-path gate deg p>=3 then
gives the improved full thin-phase matrix gates

    deg_matrix(P)>=36 on2O, >=60 on2I.                    (13)

The associated-graded domain still gives a fixed bad
seed-power marker of degree<=deg(P)-deg(h). Hence at
least one marker in ANY basis of the respective
filtered space of degree

    deg(P)-12 on2O, deg(P)-20 on2I                       (14)

has infinitely many nonzero marked moments. These are
necessary degree and marker gates, not existence
certificates or finite tail cutoffs. The earlier
24/36 bounds remain weaker historical consequences;
original torus-invariant thin existence is already
closed in every degree.

Both FULL exceptional problems remain open on mixed
isotypic phases, including full non-T thin phases
and image-interior marked safety. Retain the finite-
stabilizer condition, actual Haar graph constraints,
original regularity and the independent symmetrized
height gate. No expanding degree/support/density/
moment-prefix campaign has been performed.

In a genuinely mixed phase, lower isotypic terms can
have opposite weights even when its leading tuple
has the common factor in(9). The above one-parameter
limit then need not exist for the full phase. Removing
those lower terms by matrix-degree degeneration would
lose precisely the global period information retained
here, so no full MZ conclusion follows from that step.
