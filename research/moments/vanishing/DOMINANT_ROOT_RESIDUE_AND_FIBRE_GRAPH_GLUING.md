# Dominant-root residue and fibre-graph gluing

16 November 2026. Internal written calculation. Both FULL exceptional
Mathieu--Zhao cases remain **OPEN**. External review remains deferred.
No degree campaign, new holonomic run or counterexample is claimed.

The unique dominant left-root direction localizes possible poles of the
ORIGINAL volume. The first blowup of that line has an explicit, NONZERO
residue class. Thus this region is not cohomologically empty, despite its
rational geometry. The missing step is to transport its actual class into
the corrected intrinsic boundary extension, then prove zero cycle periods.
This note computes the residue and gives the precise graph gluing equation;
it does NOT assert that this transport or its vanishing has been proved.

## 1. The first blowup exposes a symplectic quadric

Work on PGL2, since the binary centre acts trivially on every original
exceptional phase and marker. Write projective matrices as [x:y:z:w],
Delta=xw-yz, with determinant boundary B. By a left SL2 change, put the
unique strict leading root at the ruling line ell={x=y=0}. This preserves
the original invariant volume, up to its fixed normalization. The retained
[capacity theorem](CRITICAL_SECTION_CAPACITIES_AND_SHARP_SPARSE_SAFETY.md),
Section 2, proves pole localization for an eventually pure phase; it does
not prove localization of the full intrinsic extension.

Blow up ell in P3. On z=1 put

    u=w, r=Delta=xu-y, x=r v, y=r(vu-1).             (1)

The last equations give the blowup chart complementary to the strict
transform of B on the exceptional divisor E. Fix the scalar of the
original volume by writing it on z=1 as

    Omega=dx wedge dy wedge du / Delta^2.

Direct substitution, with NO phase-dependent Jacobian replacement, gives

    Omega=dlog r wedge dv wedge du,
    Res_E Omega=omega=dv wedge du.                    (2)

The coefficient is 1 in this convention. A different original Haar
normalization multiplies both forms by the same nonzero constant.

Geometrically E=P(N_(ell/P3))=P1 times P1: one coordinate is the bottom
row line, the other the first-order top row line. Its intersection with
the strict determinant boundary consists of equal row lines. Hence

    E intersect B'=diagonal,
    E^o=E minus B'=(P1 times P1) minus diagonal.       (3)

In particular E^o is the same affine quadric that
appears in the retained
[global torus bundle](GLOBAL_TORUS_BUNDLE_AND_UNIT_INVISIBLE_INFINITY_POLES.md).
The blowup is an isomorphism on the ORIGINAL PGL2 open. It is not a new
phase or an auxiliary Gaussian lift.

On w=1 use

    u'=1/u, r'=r/u^2, v'=u^2 v-u.                    (4)

In this chart the volume has the opposite coordinate sign and its
residue is -dv' wedge du'=omega. The two affine primitives are

    alpha=v du, alpha'=-v' du',
    d alpha=d alpha'=omega,
    alpha-alpha'=du/u.                               (5)

The transition term in (5) is essential. Local exactness of the residue
does not give a global primitive on E^o.

There is also an exact ambient cohomological localization. Put

    Y=Bl_ell(P3) minus the strict transform of B.

Then Y is the TOTAL LINE BUNDLE Hom(ell_bottom,ell_top) on E^o=Q2,
including its zero section E^o, and X=PGL2 is its nonzero part.
Indeed a point is the distinct pair of bottom row and limiting top row
lines, together with an arbitrary scalar map between them. Zero scalar
is the exceptional divisor. Points with bottom row zero belong to the
removed strict determinant boundary. This is the one-sided completion
of the retained global torus bundle, without a rational frame.

The projection Q2->P1 to its bottom row line has affine-line fibres.
It is locally trivial by (4), so is a homotopy equivalence. A line
bundle retracts to its zero section. Thus H^3(Y)=H^4(Y)=0, whereas
H^2(E^o)=Q(-1). The Gysin sequence for E^o in Y gives the isomorphism

    Res: H^3(PGL2;Q) -> H^2(E^o;Q)(-1).              (5a)

The usual de Rham residue of Omega is (2), with
the conventional 2*pi*i factor in Betti comparison. Consequently the
WHOLE ambient Haar line, not merely a selected local differential,
is represented by this residue. This assertion is BEFORE imposing f.
It does not identify the canonical source of a finite-phase family
with the source of this nonproper, generally non-phase-compatible Y.

## 2. Its class is nonzero and exceptional symmetry preserves it

One can verify nonvanishing without a general comparison theorem. The
antipodal row-line locus is a compact real sphere in E^o. In the chart
above it has

    v=bar(u)/(1+|u|^2),
    omega=dbar(u) wedge du/(1+|u|^2)^2.

With the complex-plane orientation its integral is

    integral_(sphere) omega
      =2i integral_(R2) dx dy/(1+x^2+y^2)^2
      =2*pi*i.                                       (6)

The two-chart transition (4) makes this sphere and form smooth at
infinity. Thus omega has no global smooth, and hence no algebraic,
primitive. Equivalently (5) records a nonzero first Chern class.

Right projective SL2 acts diagonally on the two row-line factors and
preserves the original Omega. It preserves ell, its blowup and E, so
functoriality of residue proves that it preserves omega as well.
In particular averaging over O or I does not annihilate this class.
Remove the finite fixed-point set of the projective exceptional action
on E^o. Removal of points from a smooth complex surface preserves H2
cohomology: the relative groups in degrees 2 and 3 vanish by the local
four-dimensional ball calculation. The quotient on this free open is
finite etale, and omega descends. Its descended class is nonzero, since
a primitive downstairs would pull back to one upstairs. No assertion
about the singular quotient's de Rham complex is needed.

In fact (5a) descends to the ORIGINAL exceptional quotient. All fixed
points of the projective exceptional action on Y lie in E^o, since
the action on X=PGL2 is free. They are finite in number: every
nonidentity projective element fixes just two points on P1, hence
only two ordered distinct pairs in Q2. Remove their union P from Y
and E^o. In a smooth complex threefold removing points changes neither
H^3 nor H^4, by the same local-ball argument in real dimension six.
The finite free quotient (Y minus P)/G therefore has H^3=H^4=0,
using rational transfer. Its smooth divisor (E^o minus P)/G has
H^2=Q(-1); transfer identifies this with the invariant H^2 upstairs,
which is generated by omega. Gysin gives, for G=O or I,

    H^3(SL2/Gamma;Q)
       ~= H^2((E^o minus P)/G;Q)(-1).                (6a)

Thus the entire ORIGINAL ambient Haar line has this residue
representation. The deleted points are boundary points, so no part
of the original open was removed. This still concerns an ambient
partial completion, not the proper finite-phase source.

This nonvanishing is an INPUT to the desired transport calculation,
not evidence of a failure. Pure safe phases already have nonzero
finite-phase volume residues and zero intrinsic corrected class.

## 3. The first exceptional divisor need not be over finite phase values

For a Peter--Weyl summand of matrix degree 2d, write f_d=H_d/Delta^d.
Separate H_d by its degree j in the top row (x,y). Substitution in (1)
is the EXACT finite expansion

    f_d=sum_(j=0)^(2d) r^(j-d) A_(d,j)(u,v),
    A_(d,j)=H_(d,j)(v,vu-1,1,u).                     (7)

The exponent j-d is the effective left charge. For fixed j the leading
v coefficient is H_(d,j)(1,u,1,u); in the irreducible PW summand it is
the corresponding right binary coefficient. Distinct degrees at a
fixed charge have distinct possible leading v degrees d+charge.
Consequently the nonzero summand of largest such degree prevents
cancellation of a nonzero charge coefficient on the generic E.

For an unresolved mixed phase in the leading-root gauge, a negative
charge is present. Formula (7) therefore has a genuine negative lowest
r exponent: the generic E maps to phase infinity. Finite f fibres meet
further resolved centres, including zeros of the leading polar
coefficient. The nonzero class (2) must be followed through THOSE
modifications; it is not already a class on the finite-fibre graph.

The top PW tuple has j>d by strict leading purity, but lower PW terms
can have j<d. The leading-root theorem does not remove them. Retain the
entire opposite-side capacity bound <=k-2, all allowed automorphisms,
the six/eight minimum failure support and the primitive constraints.
No bound on finitely many degrees completes this step.

## 4. What the graph equation would have to compute

The intrinsic local criterion is

    eventually pure f is SAFE for all original markers
      iff xi_(boundary,0)=0.                          (8)

Use the canonical reduced local Tate invariant domain in the
[full local theorem](FULL_INTRINSIC_LOCAL_TATE_CRITERION.md). Do not
replace it by an arbitrary unreduced fibre graph in a nonsemisimple
degeneration. Boundary genus supplies no additional independent failure
gate. The equivalent corrected component is that of the ACTUAL
geometric 1-motive, not an arbitrarily chosen torus section.

Here is the exact cycle calculation WHEN an actual Picard gluing
presentation has been supplied. On a splitting cover s=u0^e choose
the actual identity-component correction b and integer n as in the
[toric correction theorem](ACTUAL_TORUS_CORRECTIONS_AND_RAMIFIED_COMPONENT_VALUATIONS.md).
Trivialize the corrected component line bundles on the normalizations
of the boundary curves. Let g_a be their actual gluing scalars at the
oriented graph edges a, and put

    w_a=ord_(u0)(g_a)/(e*n).                          (9)

For every surviving actual character represented by a graph cycle
z=(z_a), its corrected valuation is

    <w,z>=sum_a w_a z_a.                             (10)

Changing vertex trivializations adds a coboundary to w and leaves (10)
unchanged. A supplied realization of the actual torus as a subquotient
of the graph Picard torus must retain its character lattice and the
map to the canonical reduced domain. This note does NOT identify an
arbitrary collection of node scalars with that actual presentation.

In a full graph-torus presentation, zero on every cycle is equivalent
to the existence of vertex potentials p_v with

    w_a=p_(head a)-p_(tail a).                        (11)

For an actual smaller character domain, (11) is sufficient and may
be stronger than necessary. Thus the desired mechanism is a compatible
primitive/gluing construction which makes the ACTUAL corrected
cochain a coboundary, or at least annihilates every surviving character.
The primitives (5) alone do not do this: their transition is nonzero,
and the finite-phase transport and correction b are still missing.

The unit two-form residues r_1 on a generic fibre are different data.
Under the extra semisimple-fibre hypothesis the
[actual log formula](ALGEBRAIC_HAAR_LOG_COEFFICIENTS_AND_ORIGINAL_RESIDUES.md)
says B_(Q,0)=constant*<c_0,r_Q>. Purity gives <c_0,r_1>=0;
it does not give c_0=0. The present computation does not confuse r_1,
the first-blowup residue omega, and the corrected extension cochain w.

## 5. A relative forest, rather than an ambient tree

Suppose, AS A SEPARATE PROVED INPUT, actual corrected gluing scalars
away from a designated bad subgraph are units. Then w=0 on the good
edges. Contract each connected component of the good subgraph to a
vertex, retaining all bad edges, including loops and parallel edges.
Call the resulting graph K. The cochain descends to K. It is a
coboundary on the original graph iff it is a coboundary on K: any
potential with zero differences on good edges is constant on each
good component, and conversely such potentials pull back.

Therefore if K is a forest, EVERY supported w is a coboundary. If K
has a cycle, some supported w has a nonzero cycle period. This is an
exact linear-algebra statement; it supplies no unproved localization
of the actual w. Even a bad path which is itself a tree can create
a loop when both endpoints meet one connected good component.

For a direct illustration of the ambient/fibre distinction, take
ambient A3 with coordinates x,y,z and the SNC divisors x=0 and y=0.
Their ambient graph has two vertices and one edge. For

    f=x+y+z^d, d>=2, s!=0,

the two boundary curves in F_s are connected affine lines and meet
transversely at the d distinct roots of z^d=s. The fibre graph has
two vertices, d parallel edges and first Betti number d-1. Resolving
a proper completion does not erase these cycles: the distinct finite
nodes are retained or their edges subdivided. This is a geometric
control, not an exceptional phase or an MZ counterexample.

## 6. The remaining proof, without a new classification detour

The focused obligation now has concrete input (2)--(7): track the
original nonzero residue and its transition (5) through the finite-phase
resolution, construct the actual Picard/NS gluing and abelian correction,
and prove that (10) vanishes on the canonical surviving domain for
EVERY exceptional eventually pure phase. If it does not, one must
produce an ACTUAL such pure phase and a fixed original bad marker.

Two unproved bridges must stay explicit: pole localization does not
yet localize the intrinsic extension, and the canonical residue class
has not been shown to induce a coboundary after the required transport.
The ambient residue isomorphism (5a) fixes the starting Haar class; it
does not supply either of these phase-dependent bridges.
The original all-degree pure-to-safe obligation, family-transfer result,
capacity and automorphism inputs remain intact. No ordinary-resolvent
algebraicity follows from this residue calculation, and no exponential
generating-function theorem supplies it. Both FULL cases remain OPEN.
