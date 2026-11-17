# Hodge weights of Haar variation and boundary graphs

16 November 2026. Internal written proof with classical mixed Hodge inputs
pinned below. Independent replay, formal verification and external review
are not claimed. Both COMPLETE exceptional Mathieu--Zhao cases remain
**OPEN**. Review remains deferred.

When nearby middle-fibre monodromy is semisimple at eigenvalue 1, the
actual Haar variation is confined to the lowest homological weight of
that fibre. This is dual to the top-weight, Tate part of H2 cohomology,
determined by the boundary dual graph. It has no component in the genus
graded piece contributed by boundary curves. This strengthens the earlier
intersection-radical constraint and applies before any purity assumption.

The smaller space still cannot be eliminated uniformly on the entire
pure locus: actual one-weight PURE SAFE exceptional phases already have
arbitrarily large invariant Tate boundary spaces. The required universal
result remains vanishing of the ACTUAL coupling, not of its possible
target space.

See the [checker](verify_haar_tate_variation.py),
[receipt](haar_tate_variation_verification.json) and
[source audit](haar_tate_variation_source_audit.json).

## 1. The weight of the ambient top class

Let X=SL2/Gamma for ANY finite Gamma in SL2(C). Its rational cohomology
in degrees 2 and 4 is zero and H3(X;Q) is one-dimensional. In mixed Hodge
structures, more precisely,

    H3(X;Q)=Q(-2), of weight 4.                         (1)

To check the weight rather than just the dimension, compactify SL2 as the
smooth quadric Q3={xw-yz=t^2} in P4. Its boundary is the smooth quadric
Q2=P1 times P1. The residue/Gysin sequence identifies H3(SL2) with the
kernel of

    H2(Q2)(-1) -> H4(Q3).                              (2)

Both sides of (2) are pure Tate of weight 4. The two ruling classes on
Q2 map nontrivially to the same one-dimensional H4(Q3), so the kernel is
one copy of Q(-2). Finite etale pullback from X is injective on cohomology,
commutes with mixed Hodge structures, and is an isomorphism in degree 3
by the rational homology of the compact cores. This proves (1).

For a nonconstant regular f and a nearby regular value s, put
F_s={f=s} and U_s=X minus F_s. Retain the FULL bifurcation set, including
finite bifurcations caused by infinity. The retained integral tube sequence,
dualized over Q, is the EXACT mixed Hodge sequence

    0 -> Q(-2) -> V_s=H3(U_s;Q)
                       -> K_s=H2(F_s;Q)(-1) -> 0.     (3)

Here H2 denotes cohomology; homology below has a subscript written as H_2.
The twist in (3) is essential. The logarithmic residue sequence gives the
cohomological Gysin sequence as a sequence of mixed Hodge structures.
Its dual map is the ordinary tube map, with the complex normal orientation
and period 2*pi*i. Thus the dual of (3) is

    0 -> H_2(F_s;Q)(1) --Tube--> H_3(U_s;Q)
                                      -> Q(2) -> 0.    (4)

The possible weights of H2(F_s) are 2,3,4. Consequently those of V_s are
4,5,6, and H_3(U_s)=W_(-4)H_3(U_s). The fibre need not be connected;
the bounds and sequence hold componentwise and with the same ambient line.

## 2. The mixed Hodge input and a necessary parity step

The algebraic family of complements gives a graded-polarizable geometric
variation of mixed Hodge structures. After a finite base change at the
chosen puncture its unipotent monodromy logarithm N has a RELATIVE
monodromy filtration M=M(N,W). The filtration W here is the ORIGINAL
nearby-fibre weight filtration. M is not silently substituted for W.
These are the geometric limit/admissibility results audited in Sections
2.5--2.8 and 3.1 of Brosnan--El Zein, together with the logarithmic mixed
Hodge construction. No tameness, isolated singularity or properness of f
is required.

For graded polarizability, the logarithmic weight construction expresses
the graded pieces as subquotients of the polarized variations from proper
smooth boundary strata. Their kernels and quotients remain polarizable.
The invariant polarization below is this pure-variation polarization;
it is not an assumed nondegenerate intersection form on the open fibre.

Suppose T_F on H_2(F_s) is semisimple on its generalized eigenvalue-1
space. Write V_1 for the generalized eigenvalue-1 part of the ORIGINAL
complement monodromy. Exactness and spectral projectors applied to (3)
show that, on V_1,

    image(N) subset Q(-2),    N^2=0.                   (5)

The base-change multiplicity scales N by a nonzero integer and has no
effect on these assertions. Other original monodromy eigenvalues may
have arbitrary unipotent parts allowed by geometry.

It would be WRONG to infer immediately from (5) that N is zero on gr_W.
The constant line and some quotient classes can both have weight 4.
Here graded polarizability supplies the missing step. The weight-4 pure
variation has a nondegenerate SYMMETRIC flat polarization. Its generalized
eigenvalue-1 part is nondegenerate for that pairing: different reciprocal
eigenvalues pair, so the other eigenvalues are orthogonal to this part.
The logarithm N is skew for the pairing, while (5) bounds its rank by 1.
A skew endomorphism for a nondegenerate symmetric form cannot have rank 1.
Indeed if Nx=ell(x)b, skewness with x=x0 and ell(x0)!=0 gives
<b,x0>=0; skewness with x0 and arbitrary y then gives <b,y>=0 for all y,
forcing b=0. Therefore

    N=0 on EVERY gr_W of V_1.                          (6)

In weights 5 and 6 this already follows from (5). The parity argument
was needed specifically in weight 4; it cannot be omitted or replaced
by a rank-one extension argument for an odd-weight alternating pairing.

The restriction of the relative filtration to V_1 is legitimate even
though the spectral projector need NOT preserve the moving Hodge
filtration. The finite-order semisimple part of the original monodromy
preserves W and commutes with N. By uniqueness it preserves M. Its
spectral projector therefore preserves both W and M, and the defining
relative-filtration isomorphisms restrict to its image. No claim that
this projector is a Hodge morphism at each nearby fibre is needed.

When (6) holds, those defining isomorphisms say that M on gr_W^k is
concentrated in degree k. Thus M=W: the two finite filtrations have no
off-diagonal double graded pieces. By the relative-filtration axiom,

    N W_k(V_1) subset W_(k-2)(V_1).                    (7)

This is a consequence of admissibility and (6), not a general rule for
monodromy acting on the nearby weight filtration.

## 3. The actual Haar variation is a Tate boundary class

Continue the ACTUAL compact Haar cycle from a physical large parameter,
as in the retained [cycle theorem](HAAR_CYCLE_VARIATION_AND_FIBRE_BOUNDARY_LINKING.md).
For any such continuation define the integral d_f by

    (T_U-I)gamma_s=Tube(d_f),
    d_(f,1)=pi_1(T_F)d_f.                              (8)

No purity hypothesis is made yet. Dualizing (7), or applying the relative
filtration to the dual of (3), gives

    d_(f,1) in W_(-4)H_2(F_s;Q) intersect ker(T_F-I).   (9)

In detail, the dual N lowers weight by two as well. The projected gamma_s
lies in W_(-4)H_3(U_s), so its motion lies in W_(-6). On the dual
eigenvalue-1 extension N^2=0, so N differs from T_U-I only by the nonzero
base-change factor. Strictness of the tube map in (4) identifies
W_(-6) of its image with W_(-4) of H_2(F_s), twisted by (1). This proves
the weight assertion. The kernel assertion follows from semisimplicity
of T_F at 1. It is the same ACTUAL lift as in (8), not a chosen fixed lift.

Let F_s=Y minus D be any smooth projective SNC compactification of the
nearby surface. Its boundary dual complex is a graph Delta(D), with a
vertex for each irreducible curve and an edge for each intersection point.
The top-weight identification, with its twist retained, is

    Gr_W^4 H2(F_s;Q)=H_1(Delta(D);Q)(-2),
    W_(-4)H_2(F_s;Q)=H^1(Delta(D);Q)(2).               (10)

One can see the indices directly in the logarithmic weight spectral
sequence E1^(-p,q)=H^(q-2p)(D^(p))(-p). In total degree 2 and weight 4,
its only relevant differential is

    Q^(edges)(-2) -> Q^(vertices)(-2),                 (11)

the signed incidence map. There are no triple intersections on an SNC
surface. Its kernel is H_1 of the graph. Higher differentials vanish by
the distinct pure weights. Payne Theorem 4.4 gives the same canonical
top-weight identification. The homological space in (10) is the DUAL
graph cohomology, not an unexplained identification of homology with
cohomology of a graph.

Boundary genera enter instead through H1 of the individual curves,
in weight 3 of H2(F_s). Statement (9) kills that GRADED component of
the variation. It does not choose a canonical splitting of all boundary
homology, and does not assert that an arbitrary torus representative
has a canonically separate genus and graph part.

**Tate-variation theorem.** For any finite Gamma and any nonconstant
regular f, semisimple eigenvalue-1 middle-fibre monodromy implies (9),
with (10). For an EVENTUALLY PURE f, all-original-marker safety follows
if additionally

    (W_(-4)H_2(F_s;Q))^(T_F)=0.                        (12)

In particular it follows if the boundary graph is a forest. Purity is
used only in this last conversion of zero eigenvalue-1 variation to
ALL original marked tails. Semisimplicity at 1 remains essential to
this stated theorem; it is not proved for all exceptional pure phases.

The criterion also applies to the primitive branch theorem BEFORE
primitive purity. If (12) and semisimplicity hold at every induced
primitive branch over H^-1(0), an eventually pure H(f) is safe.
If they hold at EVERY finite primitive puncture, they kill the
eigenvalue-1 logarithm on the full primitive complement system. The
[all-composite classification](PRIMITIVE_HAAR_VARIATION_AND_ALL_POLYNOMIAL_COMPOSITES.md)
then applies in every H degree: the eventually pure composites are
exactly H(c)=0 for a unique safe centre c, or none if no safe centre exists.
There is no inference that f-c itself is pure.

## 4. Even the invariant Tate target is unbounded on actual safe phases

The earlier [safe boundary-excess family](UNBOUNDED_SAFE_BOUNDARY_HOMOLOGY_ON_EXCEPTIONAL_QUOTIENTS.md)
does not disappear after this weight refinement. Retain its all-degree
existence proof: N=24, e=12 for 2O and N=60, e=30 for 2I. For EVERY
e|k with k>=2N a nonempty open set of invariant biforms of bidegree
(k+1,k-1) gives actual weight-2 phases P on X, pure and all-marker safe.
Fibre and complement monodromy at 0 are identity and the actual d_P=0.

In the projective quotient Z=(P1 times P1)/G, its two distinguished
boundary curves are the image C of the smooth biform curve and the
image A of the diagonal. Upstairs they meet transversely in 2k points.
The construction avoids ALL nonfree diagonal points, so G acts freely
on these intersections. Consequently C and A meet transversely at

    r=2k/N distinct smooth quotient points.             (13)

Resolve the remaining boundary and quotient singularities, retaining
these already smooth transverse intersections. The SNC boundary graph
contains the two strict-transform vertices with their r parallel edges.
Additional vertices and edges cannot kill cycles in a graph; there are
no 2-cells. Equivalently, blowing up such a node only subdivides its
edge. Thus its first Betti number is at least r-1. Equation (10) gives
the NEW all-degree bound

    dim (W_(-4)H_2(F_s;Q))^(T_F) >= 2k/N-1.           (14)

It is invariant because T_F=I. This grows without bound on BOTH actual
PURE SAFE exceptional loci. The earlier quadratic genus bound
2k(k-2)/N for the larger radical remains valid; (14) concerns the
different Tate graded piece and follows from intersections, not genus.
For k=48 on 2O or k=120 on 2I, (14) is at least 3 in each case. These
are formula evaluations, not expanded phases or least-degree assertions.

Therefore universal zero invariant Tate space on the ENTIRE pure locus
is also false. The new sufficient theorem is useful but not exhaustive.
The family is already in the closed one-weight sector; this is no new
MZ sector closure, and does not refute a space-vanishing hypothesis
explicitly restricted to some further unresolved candidate locus.

## 5. What is narrowed, and what still requires a whole-quotient proof

Every full exceptional failure must have either nontrivial unipotent
middle-fibre monodromy at eigenvalue 1 or an ACTUAL nonzero invariant
Tate boundary-graph coupling (9). In the semisimple case, arbitrary
positive-genus boundary contributions cannot by themselves supply it.
This statement is all-degree and retains every original marker.

The smaller possible target is nevertheless unbounded even on actual
safe phases by (14). A proof must control its coupling to the physical
Haar cycle using right-2O/2I regularity or equally restrictive geometry.
Neither rational-homology-sphere topology, boundary genus, graph rank,
semisimplicity alone, nor a different lift proves that vanishing.
Both COMPLETE exceptional MZ cases remain OPEN. No general Betti/MHS
producer, all-ring topological promotion, uniform cutoff, degree/support
campaign or external-review gate is introduced.

Classical inputs are Deligne's functorial logarithmic MHS, residue/twist,
weight spectral sequence and strictness; geometric relative monodromy
and graded polarizability as in Brosnan--El Zein; and Payne's top-weight
boundary theorem. The rank-one parity argument, its application to the
full Haar extension, the resulting actual-cycle weight restriction and
the safe-family graph bound are written above. Exact finite controls
do not machine-verify those universal geometric arguments.
