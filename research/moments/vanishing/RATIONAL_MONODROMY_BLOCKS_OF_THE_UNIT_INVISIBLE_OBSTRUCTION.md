# Rational monodromy blocks of the unit-invisible obstruction

17 November 2026. Written conditional safety theorem in all degrees, with
all genus and fibre monodromy retained. Both FULL exceptional Mathieu--Zhao
cases remain **OPEN**. External review remains deferred. No independent
replay, formal verification or universal machine certificate is claimed.

The specific unit detects either all or none of a simple rational global
boundary-monodromy block. Consequently, a nonzero actual unit residue vector
and rational irreducibility force safety under purity. The geometric group
in this statement is the phase-dependent monodromy of the ACTUAL generic
fibre boundary, not the exceptional quotient group itself. No surviving
mixed exceptional family is asserted to satisfy the new hypotheses here.

## 1. Retain the original unit and the corrected obstruction

Use the notation and theorem of
[the exact unit identity and invisible quotient](PURITY_UNIT_IDENTITY_AND_THE_CORRECTED_INVISIBLE_QUOTIENT.md).
For an original regular phase f on X=SL2/Gamma, write

    B=H^3_dR(X_K minus {f=z}), K=C(z),
    u=[Omega], a=[Omega/(z-f)], E=B/C_f,
    A=W5 E, T=E/A.

The original normalized volume and every boundary component are retained.
The top Tate system T has finite global monodromy. On one finite labelled
cover, its actual unit vector

    rho(z)=projection_T a

has algebraic coordinates: the double Laurent residues of the original
unit fibre form. The frame is CONSTANT and flat. An algebraic or rational
function vanishing below means identical vanishing, not evaluation at one z.

Let V_Q be the rational dual Tate fibre, with the standard Tate comparison
with rational generic torus cocharacters. Define

    phi:V_Q -> L_alg,       phi(v)=<v,rho(z)>,
    K_Q=ker(phi).

Let G be the finite image of GLOBAL monodromy on V_Q, and I its ORIGINAL
local inertia image at phase zero. The previous theorem proves

    eventual purity implies nu_0(sigma_f) in K_Q^I.       (1)

Here nu_0 is the actual corrected normalized Neron component, including
the abelian identity correction. No semisimplicity of the middle fibre is
assumed. Its vanishing is equivalent to all-original-marker safety under
eventual purity, by the completed intrinsic criterion. All-positive purity
uses the physical germ C_1(z)=1/z near infinity; eventual purity uses the
actual finite Laurent germ instead. Neither is a statement about separately
chosen Cauchy integrals on every component outside the compact image.

## 2. The rational invisible kernel is a global monodromy submodule

**Lemma.** K_Q is a QG-submodule of V_Q.

**Proof.** If <v,rho>=0, continue this identical zero along any loop in the
regular phase base. The original de Rham class a, and hence rho as a
section of the algebraic bundle, is globally defined. Continuation of its
pairing with the flat vector v is therefore the pairing with the transported
vector gv. Thus <gv,rho>=0. Inverse conventions for the monodromy action
give the same conclusion. Rationality of the local system preserves V_Q.
The action factors through G since this is the actual finite Tate system.
One may equivalently use a common finite Galois labelled cover and its
semilinear action on the residue functions. This proves the lemma. QED.

In particular, the rational kernel need not be a rational form of the
entire complex kernel. Complex invisible directions that contain no
nonzero rational vector cannot carry the normalized rational class (1).

**Unit-specific irreducibility theorem.** Suppose V_Q is nonzero and simple
as a rational G-module. If the ACTUAL rho is not identically zero, every
eventually pure f with these data is safe for every original regular marker.

**Proof.** Simplicity and the lemma imply K_Q=0 or V_Q. The second case
would force all coordinates of rho to be zero, since rational basis vectors
span the complexified dual Tate fibre. Nonzero rho excludes it. Equation
(1) gives nu_0=0, and the intrinsic criterion gives safety. QED.

If V_Q=0, safety under purity follows directly from (1). The theorem
requires rational irreducibility, not absolute irreducibility and not a
full complex Wronskian. It does not replace rho by a generic cyclic vector:
nonvanishing must concern this original unit. In particular, on an actual
irreducible block, a pure failure would require rho IDENTICALLY ZERO.

## 3. Multiplicities and the precise invisible blocks

Maschke's argument applies over Q: averaging a projection onto a submodule
over the finite group makes it G-equivariant. Hence V_Q and K_Q are
semisimple. Choose representatives S_i of the simple rational constituents,
and put

    M_i=Hom_QG(S_i,V_Q),
    M_i^blind={F in M_i : <F(s),rho>=0 for EVERY s in S_i}.

Then the exact description is

    K_Q = direct_sum_i sum_(F in M_i^blind) image(F).     (2)

Indeed every such image lies in K_Q. Conversely semisimplicity expresses
K_Q as a sum of simple submodules, whose embeddings into V_Q belong to
the indicated spaces. This also proves that M_i^blind is precisely
Hom_QG(S_i,K_Q). For the division algebra D_i=End_QG(S_i), these are
modules under precomposition. Let m_i^blind be their dimension over D_i,
equivalently the multiplicity of S_i in K_Q. Then

    dim_Q K_Q^I = sum_i m_i^blind dim_Q S_i^I,
    dim_Q S_i^I = (1/|I|) sum_(h in I) tr(h|S_i).         (3)

The second formula is the trace of the averaging projector. It uses the
original inertia, not invariants after a ramified cover trivializes it.

If V_Q is multiplicity-free, (2) says K_Q is a direct sum of ENTIRE simple
blocks on which the actual unit pairing is identically zero. One nonzero
pairing on a block excludes that whole block from K_Q. With repeated
isomorphic constituents, checking a nonzero central isotypic projection
does not suffice: a diagonal copy can be invisible even when other copies
are observed. The definition of M_i^blind retains exactly these relations.

Thus a pure failure must live in an actual unobserved rational block, or
in an unobserved multiplicity combination, that also survives original I.
Equations (2)--(3) identify a carrier, not an existence theorem for failure.

## 4. Separate attachment blindness from blindness within polar cycles

Choose a common Galois-compatible relative SNC model of the generic fibre.
Write Gamma_C for its ENTIRE boundary graph and Gamma_pol for the subgraph
of polar vertices and edges whose two ends are polar for the original
volume. Polar orders are preserved by continuation, so this is a
G-invariant subgraph. The unit support theorem gives

    rho in H_1(Gamma_pol;Q) tensor L_alg.

Use V_Q=H^1(Gamma_C;Q) with the preceding dual Tate identification. Collapse
each connected component of Gamma_pol to its OWN vertex, retaining every
other edge, loop and parallel attachment. Call the resulting graph Q_pol.
For these finite graphs there is an exact sequence

    0 -> H_1(Gamma_pol;Q) -> H_1(Gamma_C;Q)
      -> H_1(Q_pol;Q) -> 0.

For example, an edge cycle in Q_pol lifts by joining its endpoints inside
the connected polar pieces; the ambiguity is exactly polar cycles. This
also proves exactness for disconnected graphs component by component.
Dualizing and restricting the pairing defines

    K_pol,Q={w in H^1(Gamma_pol;Q):<w,rho>=0}

and gives the sharper exact sequence of rational G-modules

    0 -> H^1(Q_pol;Q) -> K_Q -> K_pol,Q -> 0.             (4)

To check surjectivity, extend any rational functional w on polar homology
to all graph homology. Since rho is supported on polar homology, its
pairing with any such extension is still zero. The kernel is exactly the
annihilator of polar homology, namely H^1(Q_pol;Q).

Finite averaging makes (4) split as rational representations, though not
canonically, and makes its I-invariant sequence exact. In particular,

    dim K_Q^I = dim H^1(Q_pol;Q)^I + dim K_pol,Q^I.       (5)

The first term is invisible independently of the nonzero residue values;
the second measures cancellations within the actual polar cycles. This
graph contraction computes a linear annihilator ONLY. It neither descends
the corrected section nor contracts a surface or a special-fibre model.
The attachment correction, special components and original volume remain.

If V_Q is irreducible, its dual graph homology is irreducible too. Therefore
the rational polar homology submodule is either zero or all graph homology.
Nonzero rho forces the latter. It does not follow conversely that full
polar homology forces nonzero rho; double Laurent coefficients can vanish.

The unique dominant root and capacity results constrain where polar
components can occur for pure exceptional phases. They do not remove the
first term of (4), nor prove the balanced ruling component vanishes in it.
An ambient resolution tree is not Gamma_C or Q_pol.

## 5. An all-degree graph criterion, conditional on the ACTUAL model

Suppose the entire generic boundary graph has two vertices joined by d>=2
parallel edges, and its finite GLOBAL label monodromy is doubly transitive
on the d edges. If its actual rho is nonzero, every eventually pure phase
with this model is safe, whatever the genera of its boundary components.

To prove it, orient all edges from one vertex to the other. Graph homology
is the augmentation module

    { (c_1,...,c_d) in Q^d : sum c_i=0 }.

If monodromy exchanges the two vertices, its action on oriented edges is
the permutation action times the same one-dimensional sign character;
this twist does not affect irreducibility. For a doubly transitive edge
action the commuting algebra of Q^d has dimension two: a commuting matrix
is constant on each orbit of ordered edge pairs, the diagonal and its
complement. Transitivity gives precisely one trivial constituent. If the
augmentation module were reducible, semisimplicity would give projections
onto at least two further summands, making that commuting algebra have
dimension at least three. This contradiction proves rational simplicity.
Duality and the possible sign twist preserve it. Section2 now applies.

This is a uniform geometric sufficient condition, not an assertion that
the surviving exceptional families have this graph or monodromy. Two
components alone are insufficient; every parallel edge is relevant.

## 6. Actual calibration of observed and blind blocks

For the original nonexceptional phase P=(xw)^2 on SL2, write z=v^2. The two
generic torus branches have unit residues, in the retained orientation,

    rho=(-1/(2v),+1/(2v)).

Global C2 monodromy swaps the branches. Thus V_Q=Q^2 has a trivial diagonal
block and a sign block. The actual unit sees the sign block and kills the
entire diagonal block. Original inertia at zero fixes that diagonal. The
actual corrected Neron vector is nonzero diagonal, as also follows by
recomposing the retained xw component calculation through z=(xw)^2.
This uses the corrected-vector convention; it does not identify its sign
with the physical projected Haar-cycle convention. The original marker xw
detects what the unit integer-log coefficient misses:

    C_(xw)^P(z)=(1/2) log(z/(z-1)).

This is an ORDINARY marked resolvent and one FIXED original marker. The
phase is NONPURE, h(P^m)=1/(2m+1), and is not exceptional-invariant. It
calibrates (2)--(3), not the exceptional conclusion.

On the already known pure SAFE exceptional weight-one families, rho=0
and the entire rational graph space is invisible, with unbounded dimension,
while the actual corrected vector is zero. The proper-volume exceptional
Bezout controls with corrected pairings 2,-2 also have rho=0 but are proved
NONPURE. The nonexceptional pure Long example has rho=0 and a nonzero
corrected vector. None contradicts the nonzero-rho irreducibility criterion.

## 7. Source relation and the unsolved exceptional implication

[Gavrilov--Pakovich, CMH89 (2014),125--155](https://doi.org/10.4171/CMH/314)
uses monodromy-invariant spaces for one-dimensional moment problems.
Its Theorem2.2 and the following remark discuss full symmetric and doubly
transitive monodromy; Section3 uses the invariant span of an actual cycle.
These are useful comparisons. Their Riemann-surface vanishing conclusion,
polynomial classification and inverse-branch model are NOT transferred to
these Haar surface fibres. The finite group arguments above are proved
directly and combined with the previously established all-genus corrected
unit-annihilator theorem.

For either exceptional quotient, the next obligation is to determine the
actual G-module, its original inertia, the observed multiplicity spaces
and the balanced ruling component on the surviving mixed pure locus.
In an irreducible sector only identically zero rho can conceal a pure
failure; in a multiplicity-free sector only entire unobserved blocks can.
Equations (4)--(5) isolate attachment cycles separately. Their corrected
coupling still needs a purity argument or a genuine pure producer.

No such calculation has been supplied for every surviving mixed phase.
The capacity, automorphism, primitive and full family-transfer constraints
are retained as inputs. No all-ordinary-resolvent algebraicity, Gaussian
EGF transfer, pure exceptional counterexample or full exceptional closure
is claimed. Both FULL cases remain OPEN; no external review was sought.
