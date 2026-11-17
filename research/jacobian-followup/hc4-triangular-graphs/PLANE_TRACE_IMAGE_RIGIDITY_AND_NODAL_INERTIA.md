# Plane trace-image rigidity and nodal inertia

Claims HC4MYKP0/1/2. All complex plane Keller maps and all degrees.
The results constrain the actual whole-source trace image and ordinary
nodes of the TOTAL nonproperness divisor. They do not prove that divisor
empty, exclude cusps or tangencies, or decide HC4/JC2.

## 1. Setup and the single new standard finiteness input

Let K:A2_X -> A2_Y be Keller, R=C[x,y], B=C[P,Q], F=Frac(B), L=Frac(R),
d=[L:F]. Let S normalize B in L. Retained KMFI94 gives a canonical open
X in Z=Spec(S), with complement a union of ramified and possibly additional
unramified prime divisors. PFN2 identifies their target images with the
nonproperness locus. Write this reduced curve as D=V(f), with f=1 if empty.
Retained KO1 gives the exact height-one pole support of source power traces.

Use ordinary FIELD trace and set

    M=Tr_(L/F)(R) subset B_f, A=B_f, W=D_B.

The commuting Nambu derivations dual to P,Q give the source its W-action.
Field trace commutes with these derivations, so M is a W-submodule of A.
Multiplication by B is included in W. Trace(1)=d makes B a submodule of M.

The only newly imported proved statement is finite W-generation of B_f,
from Mustata's author-hosted D-module notes, Theorem 4.15. Its proof bounds
the Bernstein growth of all finitely generated submodules to deduce
finite generation. We use the actual polynomial localization in two
variables; no generic connection, semisimplicity of the source, or
automatic trace retraction is assumed. The inspected proof and retained
byte pin are recorded in [the context](sources/plane-trace-core-context.json).
The univariate boundary-core construction in retained MVBOUNDARYSPACE is
credited as a precursor; its one-variable rank-one localization theorem
must not be transferred to two variables. The node below explicitly
refutes that transfer.

## 2. The trace image is full at every generic missing-curve point

At a height-one prime p not dividing f, M_p=B_p. At p dividing f, some
source polynomial r has a pole at an omitted divisor over p. We claim

    M_p=F.

An elementary valuation proof avoids a trace-surjectivity assumption.
If the B_p-submodule M_p of F were not F, all its elements would have
valuations bounded below by -N for some integer N>=0. In a valued
splitting field let q>0 be the number of negative-valuation conjugates of r,
with multiplicity; let their valuation sum be -a<0. For every positive
integer h, the first q power sums of the negative conjugates of r^h have
valuation >=-N: subtract the integral conjugates from Tr(r^(h*j)). Newton
identities then give lower bound -q*N for their product. Its actual
valuation is -h*a, a contradiction for large h. This also works after
scaling the extended valuation. Thus M_p contains arbitrarily deep target
poles and equals F. No bounded experiment proves this assertion.

Consequently C=A/M vanishes at every height-zero or height-one target
prime. Since W is Noetherian and A is finite over W, C is W-finite. Its
support is a FINITE set of closed points: for finite W-generators c_i,
the support equals the union of the closed supports of the cyclic
B-modules B*c_i. Derivatives cannot enlarge a generator's support, since
localization commutes with the W-action. The absence of height-one support
therefore gives codimension at least two. This argument uses W-finiteness,
not false B-finiteness of the localization or of M.

## 3. A top de Rham obstruction, and the exact trace core: HC4MYKP0

Place the affine de Rham complex of a W-module in degrees 0,1,2. Its last
cohomology is N/(partial_P*N+partial_Q*N), so H2_dR is RIGHT EXACT for
every W-module, without any holonomicity hypothesis on that module.
The de Rham complex of R with the Nambu derivations is the ordinary
source complex: dP,dQ are an R-basis of one-forms because the Jacobian is
a scalar unit. Every polynomial two-form on A2 has a polynomial primitive.
Thus H2_dR(R)=0 and, by the surjection R -> M,

    H2_dR(M)=0.

For a closed target point p, let Delta_p=W/(W*(P-P(p))+W*(Q-Q(p))).
It is simple and H2_dR(Delta_p)=C. Explicitly Delta_p=C[a,b], with
partial_P,partial_Q acting by a,b and the centered coordinates acting by
-partial_a,-partial_b. Derivatives reduce any nonzero polynomial to a
constant, and a,b then generate all of Delta_p. Hence M has NO nonzero
quotient Delta_p.

Every W-finite module supported at finitely many points is a finite
direct sum of the Delta_p. One can see this directly, without importing
semisimplicity for general holonomic modules: split the point supports,
use local nilpotence of the centered coordinates, and apply the usual
Taylor/commutator decomposition to get C[partial_P,partial_Q] tensor V_p,
where V_p is the joint coordinate kernel. Finitely many W-generators make
V_p finite dimensional. This is the point case of the equivalence already
used in MVBOUNDARYSPACE; it is not a claim about modules with curve support.

Define a submodule depending ONLY on the target curve:

    J_f = intersection ker(phi), over ALL W-maps phi:B_f -> Delta_p
          and ALL closed target points p.

Then the exact equality is

    Tr_(L/F)(R)=J_f.

Proof. The quotient C=A/M is a finite direct sum of point modules, so the
universal intersection J_f is contained in M. Conversely any phi restricts
to zero on M: a nonzero restriction would surject onto the simple Delta_p,
contradicting H2_dR(M)=0 and right exactness. Thus M is contained in J_f.
This also shows that only finitely many independent point-quotient maps
enter the intersection for an actual Keller curve. The map A -> C is the
maximal quotient of A supported at points. Point-supported quotients at
smooth points of D vanish: in formal coordinates f=u times a unit, the
localization is generated by powers of u^-1, killed by partial_v, whereas
partial_v acts injectively on Delta_p. Therefore C is supported among the
singular points of D.

This is a plane-specific rigidity statement: higher-dimensional codimension
two support need not be a collection of points. It is NOT the univariate
assertion M=B_f. It says the complete embedded linear differential trace
image is fixed by the target nonproperness curve, independent of degree,
monodromy, and the choice of a hypothetical source realization. It supplies
no multiplication on a finite cover or on a source algebra.

## 4. Exact nodal residue model: HC4MYKP1

Suppose p is an ordinary node of the TOTAL reduced curve D. After
completion choose parameters u,v with f=unit*u*v. Put O=C[[u,v]]. Then

    M tensor_B O = O[u^-1]+O[v^-1] subset O[(uv)^-1],
    (A/M) tensor_B O = H2_(u,v)(O) = Delta_p.

Here tensoring includes localization at p; completion is flat. Delta_p
has its natural formal-coefficient action through finite jets, and its
tensor with the completed base is again Delta_p. No inverse-limit
completion of the point module is being substituted for this tensor.

To verify the assertion explicitly, quotient O[(uv)^-1] by the displayed
sum. Its basis consists of u^-i*v^-j, i,j>=1. It is generated as a
differential module by u^-1*v^-1, and is the simple point module. There
are no point-quotient maps from either single-branch localization: the
other parameter's derivative kills u^-k or v^-k and is injective on the
point module. Thus EVERY point-quotient map of O[(uv)^-1] factors through
this one copy of Delta_p. The Cech quotient defines a global W-map
B_f -> Delta_p by formal expansion, and it is nonzero (1/f has a nonzero
u^-1*v^-1 leading class). Its restriction to M must vanish by Section 3.
The finite-point quotient C therefore has exactly one Delta_p at this
node, giving the claimed equality after flat base change.

In particular the Laurent expansion of EVERY whole-source field trace has
NO term with both u and v exponents negative. This excludes all mixed
principal parts, not just one scalar residue. The restriction is stable
under every target differential operator and target polynomial multiplier.
Neither O[(uv)^-1] nor a mere list of pole orders is the trace image here.
The trace image is not a ring: it contains u^-1 and v^-1 after base change,
but not their product. Global field trace is linear, not multiplicative.

## 5. Overlapping inertia at a node is impossible: HC4MYKP2

Let sigma,tau be the two geometric meridian permutations on the d sheets
around the node in Section 4. They commute. Indeed the whole finite-cover
branch locus is contained in D, so retained KMFI96's strict-local tame
Kummer/Abhyankar description applies to these two normal-crossing
parameters. This is a statement about the two actual meridians in one
common sheet labeling, not arbitrary conjugate class representatives.

We claim their moved-sheet sets are DISJOINT:

    supp(sigma) intersect supp(tau) is empty,
    |supp(sigma)|+|supp(tau)| <= d.

Suppose a common sheet were moved. In its joint orbit the abelian local
group acts transitively, so each generator either fixes the whole orbit
or fixes none of it. Both are nontrivial on this orbit. In the corresponding
complete local normalization factor, EVERY prime over u=0 and EVERY prime
over v=0 is therefore ramified. All are omitted from the actual source,
which is etale over the target. There are no other boundary prime divisors
over this formal neighborhood, since the TOTAL nonproperness curve there
is u*v=0. Retained KMFI94's pure divisorial complement description makes
the source in this factor exactly its localization at u*v.

Base changing the affine source ring to O splits the normalization factors.
The element which is 1/(u*v) on this source factor and zero on the others
belongs to R tensor_B O. Its field trace is a/(u*v), where a>=1 is the
factor's field degree, a scalar unit in C. This has a nonzero mixed Cech
class, contradicting Section 4. Thus no such joint orbit exists.

This proof retains source multiplication and the actual open complement:
the forbidden element is in the SOURCE after flat base change. An
arbitrary abstract symmetric trace matrix does not provide it. It also
uses the WHOLE nonproperness divisor being nodal at p; an extra branch,
a cusp, tangency, or unresolved singularity invalidates this local model.

An immediate group-theoretic screen is now available in EVERY degree.
At a node with nontrivial inertia on both branches, their moved-sheet
counts must sum to at most d. If every nonidentity inertia candidate in
an action moves more than d/2 sheets, such ramified nodes are impossible.
For each exceptional eleven-sheet PSL2(11) action, retained KMFI95 gives
at most3 fixed sheets for every nonidentity element (the cyclic marks
for orders2,3,5,6,11 are3,2,1,0,0). Each moves at least8 sheets, so two
nontrivial inertia branches cannot meet in an ordinary node of the total
nonproperness curve of a PLANE Keller realization. In particular a V4
crossing orbit cannot be rescued merely by its one common fixed sheet.
This strengthens KMFI96's all-dimensional image bound in the plane; it
does NOT exclude the exceptional actions themselves or import the plane
restriction to their higher-dimensional realizations.

Disjoint transpositions on four sheets remain allowed, as do two disjoint
3-cycles on six sheets. Thus the retained quartic (2+2) endpoint packets
and degree-six collision packets are not eliminated. A node may still be
omitted from the target image if there is no common fixed sheet; the result
does not prove surjectivity. It forbids overlapping moved orbits even when
other sheets DO reach the target point.

## 6. Effect on the HC4 programme

Through retained HCTANG2, these are constraints on the plane branch of
HC4. The trace-core theorem is a genuinely two-variable source result;
it is not asserted directly for arbitrary four-dimensional Keller maps.
It shows why further LINEAR differential refinements of the complete trace
image cannot distinguish two hypothetical plane sources with the same
target nonproperness curve: their embedded image is already J_f. A new
global obstruction must control actual multiplication/source-open geometry,
a marked observable such as the canonical energy and its powers, or the
global curve-complement monodromy. None is supplied automatically here.

The [checker](verify_plane_trace_core.py) verifies finite transcription
controls for point-module actions, top cohomology, Laurent Cech projections,
actual local Kummer trace factors, and the group support inequality. It
does not enumerate new candidate groups or rerun any registered checker.
The universal proofs are written above. No characteristic-zero plane
counterexample, independent replay, formal verification, external review,
or full HC4/JC2 solution is claimed.
