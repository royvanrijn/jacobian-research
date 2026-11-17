# Boundary volumes and closure of the original torus-invariant subproblems

15 November 2026. Internal written proofs using precise classical imports
and selected exact controls. Independent specialist replay, external
review and formal verification remain pending. No novelty is claimed.
See the [checker](verify_torus_volume_closure.py),
[receipt](torus_volume_verification.json) and
[source audit](torus_volume_source_audit.json).

For BOTH exceptional groups, the entire original left-torus-invariant
algebra has **no nonzero eventually pure phase**, in every degree:

    rad(ker(h_Gamma) intersect A_Gamma^T) = {0},
    A_Gamma=C[SL2/Gamma], Gamma=2O or2I.                    (1)

Here the radical means eventual vanishing of all sufficiently high
powers under the Haar functional. Consequently these restricted
Mathieu--Zhao subproblems are closed, including their image-interior
phases and every polynomial reparametrization. This is not a reduction
of every full phase to a torus-invariant one. Both FULL exceptional
MZ decisions remain **OPEN**.

The proof removes the smooth-leading restriction from the preceding
[quadric argument](PROPER_QUADRIC_MODELS_AND_LEADING_DISCRIMINANTS.md)
on the sphere models. Divisor valuations control every finite boundary
volume pole. The polyhedral orbit sizes make the required inequality
strict for every invariant leading form. A proper surface-period
argument, including nonisolated critical fibres, then excludes purity.

## 1. Proper surface periods do not need isolated critical points

Let Y be a smooth complex algebraic surface, f:Y->A^1 proper and
nonconstant, Omega a holomorphic algebraic2-form, and M a compact
oriented real2-cycle with V=integral_M Omega!=0. Then for EVERY
nonconstant complex polynomial p the moments

    integral_M p(f)^m Omega                                (2)

are not eventually zero. In particular only0 is pure in C[f]. This
extends the surface case of the previous proper-period theorem:
critical curves and multiple special fibres are now permitted.
It does not remove that theorem's isolated-critical-point hypothesis
in higher dimensions.

We verify the additional local input. Near any finite base value s0,
absorb the finite Stein-factor extension and split generic components
by a finite base change. The generic fibres are smooth projective
curves. By semistable reduction, after a further finite base change

    s-s0=t^e, e>=1,

they admit a proper regular model W whose central fibre is reduced
and nodal. The exact imports are
[Stacks0CDN](https://stacks.math.columbia.edu/tag/0CDN), all genera,
and [Stacks0CDE](https://stacks.math.columbia.edu/tag/0CDE), which
resolves nodal total-space singularities while retaining nodal fibres.
Over C the resulting local smooth charts are t=u and t=uv. An analytic
unit in the finite base map is absorbed into the parameter. The
disconnected generic-fibre case is treated component by component.
Apply reduction over the algebraic local DVR C[s]_(s-s0) and spread
the finite-presentation models to a neighbourhood of the resulting
algebraic base curve. Their analytic charts are therefore convergent;
no analytic growth conclusion is extracted from an arbitrary formal model.

The pullback of Omega is holomorphic on W, even if a chosen semistable
model does not dominate the original one. The generic-fibre
identification gives a rational map W->Y over the base. Properness
of Y over the base extends this map at every codimension-one point
of W by the valuative criterion. Its pullback form is therefore
holomorphic away from finitely many points; the coefficients of a
2-form extend across those points by Hartogs on the smooth surface W.
Equivalently use a common resolution and birational invariance of
regular top forms. This step retains the original volume form; an
unrelated semistable fibre differential cannot be substituted.

At a node, write the pulled-back form as J(u,v)du wedge dv with J
holomorphic. Its relative form for t=uv is

    eta_t = -J(u,v)du/u,  dt wedge eta_t=Omega_W.           (3)

On a fixed sector, a flat fibre cycle can be represented by finitely
many compact regular portions and plumbing-annulus portions. J stays
bounded. Regular portions have bounded integrals; annular circle and
radial portions are bounded by C(1+|log t|). Each fixed continuation
branch has such a representation. Thus every continued fibre period
of Omega_W/dt has at most logarithmic growth. Since ds=e*t^(e-1)dt,
the ORIGINAL relative form has periods bounded by

    C |s-s0|^(-1+1/e)(1+|log(s-s0)|),                      (4)

which is strictly integrable. No bound on e is needed.

Here is why this controls the compact Cauchy period, not just its
fibre residues. On the regular nodal model W, cut a small base disk
along arcs to its finitely many nodes. As in the preceding
proper-period proof, a compact cycle decomposes into a regular part,
thimbles and tubes; their Cauchy terms are integrals of the periods
in(3) over those arcs and their residue additions. Local relative
pieces have the same logarithmic annulus bound. The original local
Cauchy period is recovered by finite-base-change trace: sum the
pullback chains over all generic sheets and divide by their total
degree. A compact chain can first be put transverse to the finite
branch locus; its full pullback has no new boundary. Proper
modifications preserve the period, since the forms extend and the
exceptional vertical pieces have zero2-form integral. This argument
also applies to relative chains whose regular boundary terms are
kept fixed. It justifies the local trace on a disk; no global base
cover parametrized by A^1 is being assumed.

The already proved Cauchy estimate and polynomial root trace send
these logarithmic terms through s-s0=t^e to a power strictly greater
than-1, with finite logarithms. Hence EVERY continuation branch of

    A_f(s)=integral_M Omega/(s-f)

has only subpolar finite singularities. The same is true for A_p(s)
with f replaced by p(f), either by another root trace or by applying
the surface argument to the proper composition p composed f.
If A_p is rational it has no finite pole and tends0 in its original
infinity germ, so is zero. Eventual vanishing in(2) would make A_p a
finite polynomial in1/s with nonzero coefficient V at1/s, a
contradiction.

Semistable existence and the nodal regular model are classical
imports. The extension of the original volume, the sector bounds,
the compact-chain trace and this MZ application are internally
written arguments, not machine or external certification.

## 2. A valuative certificate for finite boundary volume regularity

Let X=barQ_N be the smooth projective completion of the affine
quadric Q_N={q=1}, let B={u=0} be its boundary, and let F have minimal
polynomial degree D. Locally write

    F=H/u^D,  Omega=u^-N J du_1...du_N,

where H is the homogenized numerator in a local trivialization and
J is a holomorphic unit. Resolve the projective graph of F, keeping
the resolution an isomorphism on Q_N, and restrict over the finite
base A^1. This produces a smooth variety Y proper over A^1. The
pulled-back Omega is a meromorphic top form on Y.

For any divisorial valuation E over X with a_E=ord_E(u)>0, let
A_X(E)=1+ord_E(K_(resolution)/X) be its log discrepancy. Its volume
order is exactly

    ord_E(Omega)=A_X(E)-1-N*a_E.                           (5)

If E is a boundary divisor remaining over a FINITE value of F, then

    ord_E(H)>=D*a_E.                                      (6)

Such a divisor is centred on the base locus B intersect{H=0}.
Suppose at every base point P one has

    lct_P(H)>N/D.                                         (7)

The valuative definition of log canonical threshold gives

    A_X(E)>=lct_P(H)*ord_E(H)>N*a_E.

All orders are integers, so(5)>=0. Thus Omega has no pole at ANY
finite boundary divisor of this graph resolution. It is already
holomorphic inside Q_N, hence is holomorphic on Y. This includes
vertical exceptional divisors over special finite values, not just
horizontal generic-fibre punctures.

The classical threshold facts used here are the log-discrepancy
definition and the multiplicity lower bound

    lct_P(H)>=1/mult_P(H),                                (8)

in [Mustata, IMPANGA notes, definition(1) and Property1.18](https://www.impan.pl/~pragacz/mustata1.pdf).
For a base point whose leading restriction H|B has multiplicity m,
the total local multiplicity of H is at most m. Therefore the
elementary sufficient condition

    D>N*m at every leading base point                     (9)

implies(7). Lower-degree coefficients may decrease the total
multiplicity and improve this certificate; none are discarded.
The strict inequality is sufficient, not a sharp necessary bound.
Failure of(7) or(9) is not a purity or failure certificate.

For N=2, the surface theorem in Section1 applies to the resulting
proper model regardless of its critical curves. Thus(7), or(9),
closes the whole polynomial subalgebra C[F] at nonzero real-sphere
volume. F can be complex and have image with interior. In higher
dimensions this section proves the volume extension, not a new
nonisolated-period theorem; the earlier period theorem still needs
its stated critical-point hypothesis there.

## 3. Polyhedral orbits make the certificate automatic on the sphere

Take N=2, Q_2={a^2+b^2+c^2=1}. The infinity conic B is P^1, and
the leading restriction of a minimal-degree nonconstant F is a
nonzero section of degree2D on this P^1. If F is invariant under
the proper octahedral or icosahedral rotation group, its zero
divisor is invariant too. This remains true independent of a
linearization character; we only use the zero divisor and its
multiplicities.

A point stabilizer in a finite subgroup of PSU2 is cyclic: under
the Riemann-sphere identification it consists of rotations about
the corresponding axis. The maximum rotation orders are4 in the
octahedral group of order24, and5 in the icosahedral group of
order60. These orders are also in the retained exact projective
element distributions. Consequently EVERY orbit has size at least

    6 for O, 12 for I.

If a root has multiplicity m, all points of its orbit have the
same multiplicity. The total degree2D therefore gives

    2D>=6m for O, 2D>=12m for I.

In particular

    D>=3m>2m for O, D>=6m>2m for I.                        (10)

More explicitly, the proof of(5)--(8) gives volume order at least
a_E-1 for O and4*a_E-1 for I at every finite boundary divisor.
Both are nonnegative since a_E>=1.

This is exactly(9) in dimension2, in EVERY degree, for EVERY
invariant F and every choice of its lower terms. Hence every
nonconstant invariant sphere function has a proper model carrying
its original holomorphic volume form. By Section1 it is not even
eventually pure. Nonzero constants are immediate; zero is the
only eventually pure element of C[Q_2]^O or C[Q_2]^I.

For comparison, a cyclic rotation group has fixed points on this
P^1, so(10) does not hold. Its invariant isotropic heights can
have one-sided weights and be pure. The tetrahedral minimum orbit
size4 gives only D>=2m; the strict argument does not settle its
equality cases. No unproved extension to those groups is made.

## 4. The exact original algebra and all automorphic images

For the original left diagonal torus T in SL2 and right Gamma,
the actions commute, so

    A_Gamma^T=(C[SL2]^T)^Gamma=C[Q_2]^(barGamma).           (11)

Explicit coordinates are a=xz-yw, b=i(xz+yw), c=xw+yz,
as in the preceding notes. On SU2 their law is the uniform
sphere law. The binary centre acts trivially; barGamma is O or I.
Thus Section3 proves(1) on the complete original algebras, not
only on a sparse support, low degree or selected real seed.
Complex phases whose sphere image has interior are included.
The Haar kernel in either algebra is a Mathieu subspace because
its eventual radical is zero. Any phase in this restricted algebra
with all positive moments zero is0, so all original markers also
vanish trivially. No finite universal moment prefix is computed.

Every left maximal torus is conjugate to T, and the normalized
Reynolds functional is left-SL2 invariant. Hence(1) holds for
every such torus. Moreover the retained
[intrinsic Haar theorem](INTRINSIC_REYNOLDS_AND_FINITE_QUOTIENT_DESCENT.md)
preserves Haar under every algebraic automorphism of the full
homogeneous space. Therefore NO nonzero pure phase lies in any
automorphic image of any A_Gamma^T. These images can contain
non-T-invariant functions. This is a larger exclusion, not an
assertion that their union is the full algebra.

The previous minimal-seed, cubic/Jordan, exact density and
hypergeometric arguments remain valid evidence and useful boundary
examples. Their former T-invariant existence gates are now
subsumed by(1), rather than left as active open questions. The
particular abstract quartic graph law remains a valid positive
model; no law on that graph can arise from a pure ORIGINAL
T-invariant exceptional phase. Its possible realization by a
full non-T-invariant three-dimensional seed remains open.
The separate symmetrized Hopf-height law is not a pure-phase
condition and is not decided by(1).

## 5. Coefficient rings, nilpotents and bounded-degree uniformity

Let R be ANY commutative C-algebra and put B=A_Gamma^T,
B_R=B tensor_C R. The exact conclusion after base change is

    {P in B_R: h_R(P^m)=0 for all sufficiently large m}
       = {nilpotent elements of B_R}.                     (12)

Indeed the group, a filtered basis of B and Haar moments have a
Qbar form. Write P in finitely many elements of the resulting
complex vector-space basis. At each prime of R, the field generated
by these finitely many residue coefficients over Qbar embeds into
C fixing Qbar. A nonzero residue phase would then contradict
Section4, with every moment identity preserved by that embedding.
Equivalently one may use extension of the zero-locus statement
for the bounded coefficient ideals below. Thus every coefficient
of P belongs to every prime, and is nilpotent. The finitely many
coefficients generate a nilpotent ideal, so P itself is nilpotent.
The converse follows from P^L=0. This is an eventual-radical
statement; a nilpotent need not have its FIRST moment zero.

Consequently the restricted Haar kernel is a Mathieu submodule
over every such R: a phase with all positive moments zero is
nilpotent, so all marked moments eventually vanish as elements,
even for markers in the FULL original algebra A_Gamma tensor R.

There is stronger fixed-degree uniformity. Choose an algebraic
basis b_1,...,b_s of a bounded-degree filtered space in B, form
the universal phase P=sum c_i b_i, and let I be its **ordinary**
pure coefficient ideal generated by h(P^m), m>=1. Equation(1)
and the Nullstellensatz give

    radical(I)=(c_1,...,c_s),
    (c_1,...,c_s)^L subset I                              (13)

for some finite L depending only on this filtered space. Each
of the finitely many generators in the inclusion has a finite
ideal-membership certificate over Qbar, so there is a finite prefix N
with the same inclusion into I_N=(h(P),...,h(P^N)). Thus over
EVERY R, a phase in this space whose first N moments vanish
satisfies P^L=0. In particular every fully pure bounded-degree
family has an elementwise all-marker cutoff L independent of
R and of the marker degree. No values of L or N are computed.

There is a terminating theoretical certificate algorithm:
enumerate N and L, compute the exact moment polynomials and
test (c)^L subset I_N by Groebner ideal membership over Qbar.
Existence in(13) guarantees termination. This is not an
implemented prefix campaign, and supplies no general producer
for the FULL exceptional MZ problems. It preserves the ordinary
ideals; it does not replace I by its radical or assert I reduced.

The base-change boundary is visible already over C[epsilon]/epsilon^2.
For a nonzero mean-zero harmonic seed h_O, epsilon*h_O is nonzero
and pure, but its square is0. The constant epsilon is nilpotent
and eventually pure, while its first mean is epsilon!=0. Thus
the complex conclusion 'only0 pure' cannot simply be copied
to nonreduced coefficient rings.

## 6. Every full failing phase has finite left stabilizer

For a full original regular phase P, its left-SL2 stabilizer is
an algebraic subgroup, since P lies in a finite-dimensional
rational left module. In characteristic0 it is finite exactly
when its Lie algebra is zero.

Suppose P is nonzero and pure and its Lie stabilizer contains
X!=0 in sl2. If X is semisimple, conjugate it to a diagonal
element. The zero derivative condition makes P fixed by its
whole torus, contradicting Section4. If X is nilpotent, its
exponential unipotent subgroup fixes P. After left conjugation
this is the standard row shear with kernel C[z,w]^Gamma.
Purity removes its constant term. All its remaining torus
weights have one strict sign, with minimum absolute weight
n_min>=8 on2O or>=12 on2I. For ANY original marker Q,

    m*n_min > deg_matrix(Q) => h_Gamma(Q P^m)=0.           (14)

This is the retained one-sided-weight safety argument, without
a homogeneity assumption on P. Conjugation preserves Haar and
the matrix filtration. Every nonzero trace-zero2x2 matrix is
either semisimple or nilpotent, so these cases exhaust X.
Therefore every nonzero pure phase with positive-dimensional
left stabilizer is safe; it cannot witness a failure.

Equivalently every possible full failing phase must have finite
left stabilizer. Its three original derivatives

    D=z*partial_x+w*partial_y,
    E=x*partial_z+y*partial_w,
    H=x*partial_x+y*partial_y-z*partial_z-w*partial_w

applied to P must be linearly independent over C. This is a
necessary condition, not an implemented general classifier.
The same statement holds for every automorphism-conjugated
left action, by intrinsic Haar and safety transport. A
finite-stabilizer phase need not be pure or failing: the
self-adjoint invariant sums f_Gamma(x,y)+f_Gamma(z,w) already
have independent derivatives and strictly positive second
moments. Conversely many known safe phases can have finite
stabilizers; none are removed from the safe class.
The degree bound(14) is for linear left conjugation. A nonlinear
automorphism transports safety without asserting that identical
original-degree bound.

Both FULL exceptional MZ questions therefore remain open on
the remaining singular-leading, finite-stabilizer architecture,
outside all automorphic torus-fixed subalgebras. Retain full
non-T-invariant thin phases, their conditional fixed bad
markers, the full image-interior marked-safety problem, actual
Haar graph constraints and the independent symmetrized height
gate. No expanding degree/support/density/prefix campaign is
proposed.
