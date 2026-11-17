# Radical plinth rigidity and automatic quotient recovery

Canonical claims: HC4MYKD0, HC4MYKD1, HC4MYKD2. Mathematical authority:
[MATH_STATUS.json](../../MATH_STATUS.json). This removes the source
codimension and one-reached-value assumptions from JV0 on affine four-space.
The remaining embedded kernel and radical plinth hypotheses are substantive;
no new action on S_rho or K2 is constructed.

## 1. The global presentation needs no extra source codimension hypothesis (KD0)

Let R be the coordinate ring of a smooth irreducible factorial affine
complex fourfold X. Let delta be a nonzero locally nilpotent derivation
with actual embedded kernel

    A=ker delta=B[h],  B=C[b1,b2],

where b1,b2,h are algebraically independent. Suppose, as an equality
of ideals rather than merely of zero sets,

    pl(delta)=I A,
    I=(F,G) subset B the radical ideal of N>=1 points b_i.    (1)

Reduced finite point ideals in a polynomial plane have such a generating
pair: move the points to (0,i) by JV2's actual interpolation, then use
the pair (X,product(z-i)). No globally Keller generators are assumed.

The kernel of an LND in characteristic zero is factorially closed in
its source domain. Indeed the delta-degree of a nonzero product is
the sum of the two delta-degrees: the highest nonzero Leibniz term
has nonzero binomial coefficient and nonzero factors. If the product
is invariant, both degrees are zero. Consequently F and G, coprime
in A, are coprime in R. A common prime divisor in the UFD R would,
along with its complementary factors, lie in A and contradict that
coprimality. Thus V(IR) has codimension at least two automatically.
In particular

    R=R_F intersect R_G                                (2)

in the common fraction field. Normality and the absence of any common
height-one zero justify this intersection.

Since F,G are in the plinth, choose global p,q with

    delta(p)=G, delta(q)=-F, a=Fp+Gq in A.

The slice theorem gives R_F=A_F[q], R_G=A_G[p]. Their transition is
a/(FG), with the ordering of JV/JU. Equation (2) gives the exact
annihilator identity

    pl(delta)=ann_A[a/(FG)]=IA:a.                       (3)

For clarity, a coefficient f kills the Cech class precisely when
local primitives with delta-derivative f glue; (2) makes the glued
function global. The injective map A/IA into A_FG/(A_F+A_G), proved
in JU by the regular-sequence colon calculation, then gives the
second equality. Since A/IA=product_i C[h], (1)--(3) imply

    a_i(h):=a(b_i,h) is a nonzero polynomial for every i. (4)

Repeat the presentation argument in JV0 before its one-point-image
assumption is used. It gives the entire source, not merely an open part:

    R = S := B[h,P,Q]/(FP+GQ-a(b,h)).                   (5)

Here are the details needed for the enlarged source scope. The equation
is primitive and irreducible, being linear in P with coprime coefficient
data. The slice charts identify S_F,S_G with R_F,R_G. Outside V(F,G),
S is smooth. At b_i and a root h=r of a_i, a singular point can occur
only if a_i'(r)=0; then the two base derivative equations determine
P,Q uniquely, since J_b(F,G)(b_i)!=0. Thus the singular locus of S
is finite. The four-dimensional hypersurface domain is Cohen--Macaulay
and regular in codimension one, hence normal. By (4), V(F,G) in S
is a finite union of two-dimensional planes, or empty. Normality
therefore gives S=S_F intersect S_G. Together with (2), this proves
(5). These are the same normality and intersection imports already
checked in JV0; polynomiality of the source was not used here beyond
factoriality, smoothness and dimension four.

Smoothness of X now forces every a_i squarefree, including the
possibility that a_i is a nonzero constant. Write

    d_i=deg a_i>=0,   M=sum_i d_i,
    E=V(IR)=disjoint union_(i,r:a_i(r)=0) A2_(P,Q).     (6)

Projection pi:X -> A2_b is smooth. Off the b_i it is a Zariski locally
trivial A2-bundle, using h and a slice as coordinates. At E the
nonzero a_i'(r) lets h be solved locally, so the two base parameters
remain independent. Each component of E has trivial complex normal
bundle of rank two, identified by pi with the normal space at b_i.

The quotient map to Spec A reaches exactly the d_i points
(b_i,r), a_i(r)=0, on the i-th plinth line. At each reached point
its fibre is A2; elsewhere on that line it is empty. No restriction
on d_i has yet been inferred from radicality of the plinth alone.

## 2. Integral homology forces exactly one reached value (KD1)

Under Section 1's hypotheses, let n_0 be the number of indices with
d_i=0. Then the entire positive integral homology of X is

    H3(X,Z)=Z^(n_0),
    H4(X,Z)=Z^(sum_i max(d_i-1,0)),
    Hk(X,Z)=0 for k>0, k not 3 or 4.                  (7)

In particular these groups have no torsion. Here homology is of
complex points with the Euclidean topology.

Put U=A2_b minus {b_i}. The bundle pi:X minus E -> U has contractible
fibre C2, so it is a homotopy equivalence. Its local trivializations
are actual polynomial charts, and a topological fibre bundle with
contractible fibre over this paracompact CW-type base has this
homotopy equivalence. Therefore X minus E has reduced homology
Z^N in degree 3 and zero otherwise, identified by pi with the
oriented small three-spheres around the b_i.

The rank-two complex normal framing in (6), tubular-neighbourhood
excision and the oriented Thom isomorphism give

    Hk(X,X minus E;Z)=H_(k-4)(E;Z),

which is Z^M at k=4 and zero at other degrees. The long exact sequence
of this pair reduces to

    0 -> H4(X) -> Z^M --T--> Z^N -> H3(X) -> 0.        (8)

It is essential to identify T, rather than just compare ranks. A
generator associated with E_(i,r) is a small normal four-disk at
one point of that plane. Its boundary projects with degree +1 onto
the small sphere about b_i: the map of complex normal spaces is an
isomorphism. It does not meet any other b_j. Thus

    T(e_(i,r))=e_i.                                  (9)

Naturality of the relative boundary and of the normal identification
makes (9) independent of local charts. The block for i is the row
of d_i ones; it has image Z e_i if d_i>0 and zero if d_i=0. Its
kernel has basis e_(i,r)-e_(i,r0) for r!=r0. Equations (8)--(9)
prove (7) for every N and every list of degrees, without a finite
degree experiment or a cancellation of contributions from different
plinth lines.

The smooth-pair purity theorem and its transverse naturality are
explicitly supplied by [Asok--Dubouloz--Ostvaer, Theorem 2.3.1](https://arxiv.org/html/2112.08241v2).
The full primary cache was retained in JZ; the naturality clause
and its transversality condition were rechecked here. Both pairs
(E in X) and ({b_i} in A2) are smooth, E is the full inverse image,
and the map on normal bundles is an isomorphism as verified above.
On complex points the calculation is precisely the usual Thom
calculation just given. It needs neither a global section of pi nor
properness. The source HTML's reversed ambient/complement names
in one displayed purity formula are not used; the smooth-pair
statement, naturality clause and proof fix the conventional notation.

It follows that

    H3(X,Z)=H4(X,Z)=0  iff  d_i=1 for every i.          (10)

An Euler characteristic test alone is insufficient. In (5) take

    B=C[x,z], F=x, G=z(z-1),
    a(z,h)=1-z+z h(h-1).

Then a_0=1 and a_1=h(h-1), so (d_0,d_1)=(0,2). The presentation
is smooth by the same derivative test. Its Euler characteristic is
1+M-N=1, but (7) gives H3=H4=Z. A missing line and an extra reached
value cannot cancel in integral homology. This is an exact negative
control within the smooth presentations, not an affine-space example.

When (10) holds, write a_i=lambda_i(h-c_i), lambda_i!=0. Interpolate
lambda,d in B and correct P,Q by invariant polynomials to turn (5)
into

    R=B[h,P,Q]/(FP+GQ-lambda(b)h-d(b)).                (11)

The row (F,G,-lambda) is unimodular. JV0's Quillen--Suslin completion
and explicit inverse matrix now construct globally

    R=B[U,V], h=d0(b)+U f(b)+V g(b), (f,g)=I.          (12)

This step uses the retained projective-module theorem, not a new
cancellation assumption. Delta becomes a B-linear translation.
Its coefficient ideal is I: compare its image ideals through the
invertible completion in (11), then contract to B. Its primitive
linear invariant generates B[h], forcing h to be affine-linear
in U,V with a constant nonzero scale, as in JV0. All coordinates
and the polynomial inverse are global.

Thus a smooth factorial affine fourfold with H3=H4=0 and the exact
kernel/plinth data (1) is itself affine four-space, with the simple
cotangent presentation (12). This is a recognition theorem under
those action hypotheses; arbitrary acyclic fourfolds are not classified.
For the actual source X=A4, (10) is automatic. Both JV0's source
codimension and its one-reached-point assumptions can be omitted.

## 3. The remaining HC4 action gate has only kernel and plinth data (KD2)

For any polynomial h on A4, the following conditions are equivalent:

1. There are actual full polynomial coordinates (b1,b2,U,V), with
   a polynomial inverse, in which h has the simple cotangent form
   (12), with nonempty finite radical coefficient zero scheme.
2. There is one nonzero h-preserving LND delta with actual embedded
   kernel C[b1,b2,h], the three generators algebraically independent,
   and exact plinth ideal I C[b1,b2,h], where I in C[b1,b2] is the
   radical ideal of a nonempty finite set of points.

The implication 2 -> 1 is Sections 1--2 applied to the actual A4
source. No source codimension, quotient image, fixed-plane, generic
freeness, separately supplied slices or affine-socle hypothesis is
an additional input. The plinth equality itself supplies the two
global local slices; acyclicity supplies their one-root residue form.

For 1 -> 2, take delta=g partial_U-f partial_V. Coprimality gives
ker delta=B[h] by the height-one translation-kernel argument in KC.
Its exact plinth is (f,g)B[h]. The proof in JU0 uses a formal local
(f,g) chart at each reduced point and the universal raising/lowering
weight identity; it needs only simple zeros there, not a globally
constant coefficient Jacobian. It therefore proves this extension
as well. The reverse inclusion follows from delta(U)=g,delta(V)=-f.
Its fixed planes, their contracted values, generic freeness and
the full two-slice presentation are consequences.

Every such h has exactly as many finite Morse critical points as
there are points of V(I); one occurs over each coefficient zero,
and there are no others. Its critical Hessian determinant is the
square of the corresponding ordered coefficient Jacobian. Hence
the number of plinth lines need not be specified separately when
the native critical count is already known.

Suppose now h has N>=2 finite Morse points, a common nonzero critical
Hessian determinant, and critical-value multiset c. By KB0, condition
2 is equivalent to an actual source equivalence h~R_(N,c). By KB1,
its full-source unit-Hessian existence is therefore exactly

    E_N AND condition 2,                              (13)

where E_N is existence of one actual N-point fibre of a plane Keller
map. E_N is independent of this pure alignment criterion.

In particular S_rho~R_rho, for every rho including collisions, is
equivalent to condition 2; its point count then must be three. For
K2~h1 it must be two. Neither action is constructed here. The
universal affine frontier remains the two endpoint inversion
requirements, absence of E2, and not E3 OR no rho satisfies condition
2. A nonempty radical plinth is essential; the unit-ideal branch is
the separate coordinate endpoint. The singleton inversion endpoint
also remains distinct from native quadratic unit existence.

The exact kernel splitting B[h] and the radical plinth ideal in that
same embedded kernel have not been inferred from freeness, from the
radical of the plinth alone, or from collective flexibility. The new
criterion recovers missing-image geometry after those algebraic data
are supplied. It does not classify all LNDs, their kernels or tori.

## Sources and checks

The [JV recovery proof](TWO_SLICE_RECOVERY_AND_ORIENTED_COEFFICIENT_MODELS.md)
is retained unchanged; its extra assumptions are removed above, not
silently deleted from that proof. JU, KB and KC supply the exact
plinth, normal-form, Keller-gate and collective-action boundaries.

The current primary [Shakhmatov cylinder/plinth note](https://arxiv.org/html/2605.00138v1)
was read at Theorem 2.6, Corollary 2.7 and the relevant counterexamples.
It identifies the union of principal cylinders for one action with
the complement of its plinth zero locus. This gives context for the
remaining action gate; no recognition theorem is imported from it.
It does not make the exact plinth radical or supply the kernel splitting.
The nonnegative-log-Kodaira descent hypothesis in the inspected
[Dubouloz torsor criterion](https://arxiv.org/html/1610.01409v1)
does not provide a canonical punctured-plane quotient here.

[The checker](verify_radical_plinth_rigidity.py) controls the universal
presentation/smoothness algebra, normal map and integral boundary
matrix, the Euler-balanced negative example and the exact gate logic.
The all-degree recognition and homology proofs are in Sections 1--3,
not consequences of a finite list of examples. [The context](sources/radical-plinth-rigidity-context.json)
pins retained evidence and the precise primary naturality import;
[the receipt](receipt-radical-plinth-rigidity.json) records the controls.
Independent replay, formal verification and external review are not
claimed. Source alignment, E2/E3, the original endpoint inversions
and unrestricted HC4/JC2 remain open; the original goal is not complete.
