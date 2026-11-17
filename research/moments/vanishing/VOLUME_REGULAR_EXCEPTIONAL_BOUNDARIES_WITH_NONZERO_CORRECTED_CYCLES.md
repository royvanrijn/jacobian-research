# Volume-regular exceptional boundaries with nonzero corrected cycles

17 November 2026. Written deduction combining the retained ORIGINAL-volume
and intrinsic transfer theorems with the new ruling presentation. Both FULL
exceptional Mathieu--Zhao cases remain **OPEN**. Review remains deferred.
Independent replay, formal verification and external review are not claimed.

Volume regularity alone does NOT make corrected cycles vanish, even on the
actual exceptional quotients and even with a unique strict dominant leading
root. The existing equivariant Bezout phases give a uniform obstruction.
Below an actual rational character has corrected coupling exactly 2 at zero
and -2 at minus one, while the original volume is holomorphic on a proper
model over the whole finite phase line. These phases are ALREADY PROVED
NONPURE; this is not an exceptional MZ counterexample.

## 1. Keep the actual original family and its fixed marker

On X=SL2/Gamma, Gamma=2O or 2I, retain the regular Bezout maps

    Phi_b=[[U+F_0 b,V+G_0 b],[F_0,G_0]],
    (n,m)=(8,12) or (12,20).

Here F_0,G_0 are the coprime bottom-row Klein forms of degrees n,m.
For this note restrict b to the equivariant weight n+m branch, including b=0.
Then Y_b=V+G_0 b has weight n, and F_0 has weight -n. The actual phases are

    g_b=(1+Y_b)(1+F_0+Y_b F_0)
        =1+Y_b+F_0+2Y_b F_0+Y_b^2 F_0.              (1)

Every function in (1) is regular and right-Gamma invariant. The retained
primitive integral H^3 action of Phi_b is exactly 2, not its generic field
degree. Its normalized Jacobian is 2W_b. For the SAME fixed original marker

    C_b=W_b(1+Y_b),
    h_Gamma(C_b R(g_b))=integral_0^1 R(-t)dt          (2)

for every polynomial R. In particular its actual ordinary Cauchy function
is log((s+1)/s), whose derivative has residues -1,+1 at 0,-1. Neither the
Jacobian W_b nor the original marker is discarded.

These statements are the retained
[cohomological transfer and uniform-marker theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md),
not a new calculation of finitely many moments.

## 2. The original volume is holomorphic on the full finite-phase model

The charged decomposition of (1) is

    g_b=p+A+B,
    p=1+2Y_bF_0, A=Y_b(1+Y_bF_0), B=F_0,
    weights 0,n,-n.                                 (3)

Both charged functions are nonzero. The retained
[constant-volume correspondence theorem](CONSTANT_VOLUME_MODELS_AND_ACTUAL_MARKER_DEFECTS.md),
Sections 1--3, proves that (3) has a smooth proper ORIGINAL model over A1
carrying the SAME holomorphic Omega. It transfers the auxiliary model by
a finite function-field correspondence with a NONZERO CONSTANT ratio of
top forms; it does not drop the nonconstant Jacobian of Phi_b. The valuative
conclusion is

    ord_v(Omega)>=0 for EVERY divisorial valuation
                  for which g_b has finite value.    (4)

Choose a smooth projective resolved pencil W->P1 dominating that model,
the fixed projective quotient and the root blowup used in the ruling
construction. Such further modifications preserve (4). Put D=W minus X.
For a regular finite s, every component of the genuine boundary C_s=D
intersect Y_s is contained in a volume-regular horizontal boundary divisor.
All relevant finite-value exceptional divisors remain volume-regular as
well. Ramified base change pulls a regular top form to a regular top form;
it does not justify omitting any component or phase multiplicity.

Thus on this full finite-phase model EVERY genuine generic boundary graph
component is volume-good, in the original-form sense. This is stronger than
merely choosing a raw ruling representative trivial off a root. It is not
the assertion that every marked Q Omega extends: (2) already detects the
failure of that assertion in the actual global period image.

## 3. These phases also have a unique strict dominant leading root

This extra fact does not require purity. Use the original matrix-degree
filtration. Its associated graded ring is C[x,y,z,w]/(det), a domain; finite
Gamma invariants retain a domain. Hence highest degrees multiply without
loss. Let d_Y>0 be the matrix degree of Y_b. The unique highest-degree term
in (1) is Y_b^2 F_0, of degree

    2D=2d_Y+n.

It has standard left weight n. In the highest Peter--Weyl summand its left
weight space has dimension one, so its common left factor is proportional to

    X^(D+n/2) Y^(D-n/2).                            (5)

The multiplicity D+n/2 is strictly larger than D; the other is smaller.
This gives a unique strict dominant left-root direction for every b in
this branch. The right coefficient remains exceptional invariant. Thus
exceptional invariance, strict leading-root instability and even (4)
together are still insufficient for corrected-cycle vanishing.

The all-degree capacity and automorphism restrictions on PURE candidates
are not contradicted. These phases have the already excluded three-weight
form (3), so they do not enter the surviving exceptional pure locus.

## 4. An explicit actual character has corrected component 2

On the target SL2 retain Long's phase

    F=(1+y)(1+z+yz).

Its canonical boundary 1-motive has torus G_m and Kummer coordinate

    sigma_F=s/(s+1),                                (6)

up to a constant factor and the already recorded rational isogeny
normalization. Choose the rational target character so its normalized
orders at 0,-1 are respectively 1,-1. Constant factors have zero order.

The canonical complement pullback, AFTER removing compact-middle systems,
sends the ambient unit by 2. Dualizing gives the actual morphism of geometric
1-motives in the rational isogeny category

    M_(g_b) -> M_F,
    lattice map 2, semiabelian map h_b,
    h_b(sigma_(g_b))=[2]sigma_F.                     (7)

This is the retained intrinsic transfer theorem before purity. Restrict h_b
to the ACTUAL generic torus and pull back the character chosen in (6).
It gives a rational intrinsic character chi_b. Neron functoriality yields

    <nu_0(sigma_(g_b)),chi_b>=2,
    <nu_(-1)(sigma_(g_b)),chi_b>=-2.                 (8)

Here chi_b is allowed in X^*(T_(g_b)) tensor Q. Clearing a denominator gives
an honest character with the correspondingly multiplied NONZERO integer
value. Equation (8) does not assert that chi_b is primitive integral, or
that every coordinate of the entire component vector is known.

This calculation retains the actual abelian correction automatically:
a character sends an identity-component lift to a unit, so that correction
does not change its order. It uses the generic torus, not a larger torus
from abelian degeneration. It does not infer an intrinsic class merely
from a marked logarithm or a raw divisor residue; (7) supplies the actual
geometric comparison.

## 5. Feed this obstruction directly into the ruling calculation

Use ANY compatible full marked semistable boundary presentation of the
proper-volume model in Section 2. The
[ruling construction](RULING_CLASS_BALANCING_AND_THE_ORIGINAL_HAAR_SECTION.md)
gives a genuine line bundle L_K and actual quotient psi:Pic^0(C_K)->G_(g_b),
with psi[L_K]=[a]sigma_(g_b), a>0 the retained lattice multiplier. For local
parameter s=u^e at zero, put

    delta_(v,j)=deg(L|E_(v,j)), M_v d_v=delta_v,
    z_b=psi^*chi_b in H_1(Gamma_C;Q).

The previously proved corrected raw-degree formula and (8) give exactly

    (1/(a e)) sum_edges (z_b)_edge
        (d_tailAttachment-d_headAttachment)=2.       (9)

At minus one, use s+1=u^e and the corresponding models; the same normalized
formula is -2. This is an actual test character/cycle for the Mathieu
extension calculation, for every equivariant b on BOTH original quotients.
No finite degree collection, unspecified cut graph or scalar branch norm
replaces it. All components of the model used in (9) are volume-regular
at finite phase values; some actual corrected attachment period is nonzero.

Therefore the implication

    original Omega holomorphic on a proposed full good boundary
       => all its corrected cycle periods vanish               (10)

is FALSE before purity, even with exceptional invariance and a strict
dominant leading root. Contracting every connected component of Gamma_C
would leave a graph without edges, but violates the exact bounded-descent
condition: (9) is a nonzero obstruction on the very lattice being discarded.
An ambient resolution tree would not alter (9).

## 6. What this does and does not change for full MZ

The retained all-Long theorem has ALREADY proved that every nonzero
polynomial H(g_b) in this equivariant branch is NOT eventually pure.
Do not reopen a b or H search here. The fixed marker (2) and nonzero
cycle (8) do not construct an exceptional pure failure.

The result sharpens the missing proof obligation. One cannot justify
good-component contraction from volume regularity plus exceptional
symmetry and strict leading-root localization alone. A successful full
argument must use an ADDITIONAL purity consequence that controls the
balanced NS/Picard class on the actual good cycle lattice, while retaining
all components and node attachments. The existing capacity, minimum-support
and automorphism exclusions are useful inputs and already exclude these
controls; whether they and the remaining purity identities force (9)'s
analogue to be zero throughout the surviving pure locus remains OPEN.

All prior goals, original volume, markers, branches, primitive constraints
and completed family transfer are preserved. Both FULL exceptional MZ
cases remain OPEN; external review remains deferred.
