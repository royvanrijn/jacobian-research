# Source-module layers and the trace-zero kernel

Claims HC4MYKQ0/1/2. Complex plane Keller maps, all degrees throughout.
This note classifies the semisimple layers of the ACTUAL source module.
It does not construct that source, its algebra multiplication, or an HC4
counterexample, and does not prove a universal inverse.

## 1. Attribution audit and retained geometry

The omitted-boundary structure needed here was ALREADY imported in
[PF2BATT1](../jc2-boundary-attachment/PROOF.md), from Borisov's published
Stein-boundary theorem. The affine omitted primes are disjoint, each has
unibranch A1 normalization, and every singular point of the affine finite
normalization is a cyclic quotient singularity on that boundary. Each
boundary prime has at most one such point and is smooth elsewhere.

In particular the nodal moved-support restriction in KP2 was already
implied by this stronger retained boundary theorem: two ramified local
branches cannot meet at a single omitted point. The KP trace-core proof
is an INDEPENDENT DIFFERENTIAL DERIVATION of that exclusion. KP0's exact
embedded trace image and KP1's exact Cech image remain the additional
trace statements; KP2 should not be advertised as a newly excluded
geometric architecture beyond PF2BATT1. All original proof/checker/receipt
bytes are preserved, and the authoritative scope/dependencies record this
attribution correction.

Write K:X=A2 -> Y=A2, B=C[P,Q], R=C[x,y], Z=Spec(S) its full finite
normalization, pi:Z -> Y, j:X -> Z, and E=Z-X. Let E_1,...,E_b be its
omitted primes. Let D_1,...,D_r be their DISTINCT target images, and f_i
the degree of E_i^nu -> D_j^nu. Put d=[C(x,y):C(P,Q)]. The source R is
a module over W=D_Y through the commuting Nambu derivations.

Two further standard inputs are screened in [the context](sources/source-module-layers-context.json):
regular Riemann-Hilbert and its direct-image compatibility; and the finite
(small) direct-image identity for intersection complexes. Finite monodromy
is semisimple by averaging over its finite group. We also use Borisov's
Theorem3.5: each normalized boundary map is linear or, after affine
coordinate changes, t -> t^(f_i). This is a newly used part of the already
cached primary; PF2BATT1 had used only the weaker unramified consequence.
We do NOT infer endpoint-curvette hypotheses for all cyclic quotients from
a one-dicritical local theorem, or conclude that every normalization is smooth.

## 2. The boundary localization sequence: HC4MYKQ0

Use perverse normalization IC_Z=C_Z[2]. This equality holds because Z is
a rational homology manifold: at a cyclic quotient point its link is a
finite free orientation-preserving quotient of S3, hence a rational
homology sphere. The same equality is automatic at smooth points.
Each E_i is homeomorphic to its normalization A1, by proper bijective
normalization and unibranchness, so IC_(E_i)=C_(E_i)[1].

There is rational local cohomological purity along E:

    R Gamma_E(C_Z) = direct_sum_i C_(E_i)[-2].

Here is a direct duality check, including singular points. For the closed
inclusion i:E -> Z, Verdier duality gives i^! D_Z = D_E i^*. The oriented
rational homology manifold Z has D_Z(C_Z)=C_Z[4], while the disjoint
topological affine lines have D_E(C_E)=C_E[2]. Therefore
i^! C_Z=D_E(C_E)[-4]=C_E[-2], proving the displayed purity. This uses
neither a smooth algebraic surface nor local flatness of an embedded knot.

The localization triangle, shifted by2, therefore gives

    0 -> IC_Z -> Rj_* C_X[2] -> direct_sum_i IC_(E_i) -> 0

in the perverse category. One may check perversity directly from the
triangle; after finite push to Y the middle term is the regular holonomic
module R, since K is affine etale and its D-module direct image is its
ordinary source ring with the Nambu action. The other two terms are
perverse, and the triangle gives the same short exact sequence there.

Consequently there is a canonical exact sequence of W-modules

    0 -> I -> R -> T -> 0,
    I = IC extension of the d-sheet permutation local system,
    T = direct_sum_i (pi|E_i)_* IC_(E_i).

Both I and T are SEMISIMPLE. For I, finite pushforward commutes with IC
extension and the generic permutation representation has finite monodromy.
For each T_i the normalized curve map is a finite power map, whose cyclic
permutation representation splits into its f_i rank-one characters.
After the finite normalization map into D_j, these remain simple IC modules.
Thus length(T)=sum_i f_i, even when the target curve is singular.

Moreover I is exactly Soc(R), the sum of its simple submodules. R is
B-torsion-free, so it contains no simple submodule supported on a curve or
a point. Every full-support simple submodule is the IC extension of a
generic irreducible summand and hence lies in I. Thus a noninvertible plane
Keller source has EXACTLY TWO semisimple layers, with no extra point layer.
Its socle is not assumed to be the finite algebra S: S need not be stable
under the target derivations. Nor is semisimplicity of the quotient a
splitting of the sequence; a nonzero curve-supported direct summand would
contradict B-torsion-freeness.

## 3. Exact trace-zero structure: HC4MYKQ1

The generic permutation representation contains its trivial summand once,
so I=B direct_sum I_0, with I_0 the semisimple IC extension of the
augmentation representation, of generic rank d-1. Let

    M=Tr(R)=J_f, N=ker(Tr:R -> M), C=M/B.

The trace on I is the finite-cover IC trace and takes values in B; it
kills I_0 and is nonzero on B. The induced map T -> C is onto. The
generic curve-supported part of C along each D_j has rank one, because
M is the full meromorphic localization at its generic point. Since C is
a quotient of semisimple T, it follows that

    C = direct_sum_j IC_(D_j).

In particular length(M)=1+r. This statement concerns the WHOLE trace image,
not a marked list of energy moments. It uses the actual retained boundary
geometry, not a general branch-number formula for every arbitrary curve.

Put T_0=ker(T -> C). The exact trace diagram gives

    0 -> I_0 -> N -> T_0 -> 0.

Every T_i has one trivial character IC_(D_j) and f_i-1 nontrivial cyclic
characters. The curve trace is nonzero on its trivial character and kills
the others. Thus T_0 consists of all nontrivial characters plus b_j-1
trivial IC copies for each D_j, where b_j is the number of boundary primes
over D_j. Its length is sum_i f_i-r. Residue choices can rescale the
nonzero row traces; no unjustified universal ramification-index coefficient
is assigned to them here.

As before, Soc(N)=I_0. Therefore the exact all-degree criterion is

    N is semisimple
      iff sum_i f_i = r
      iff EACH target nonproperness component has exactly ONE omitted
          prime over it,of residue degree ONE.

This does NOT say R is semisimple. R is semisimple iff E is empty, which
is the inverse case. In particular semisimplicity of N is automatic in
many of the hardest remaining one-boundary packets and gives no separate
contradiction to their existence.

For clarity, the cohomological budgets introduce no hidden obstruction.
Use unshifted affine de Rham degrees0,1,2. Since E_i^nu=A1, finite direct
image gives H1(T)=C^b and the other groups zero. Since R is the source
polynomial ring, its only cohomology is H0(R)=C. Hence

    H0(I_0)=H1(I_0)=0, H2(I_0)=C^b,
    H1(T_0)=C^(b-r), H0(T_0)=H2(T_0)=0,
    H0(N)=H1(N)=0, H2(N)=C^r.

The connecting map H1(T) -> H2(I) is an isomorphism, the boundary Gysin
map. These are complex vector-space statements, not integral unimodularity
claims at quotient singularities. Numerically they reproduce boundary and
curve budgets; they are not a new contradiction or a degree bound.

## 4. The one-boundary doubly transitive remainder: HC4MYKQ2

Assume there is one omitted prime, its residue degree is one, and the
generic sheet action is doubly transitive. The augmentation representation
is irreducible over C. It follows that N=I_0 is SIMPLE and

    length(N)=1, length(M)=2, length(R)=3,
    Soc(R)=B direct_sum N, R/Soc(R)=IC_D.

R is INDECOMPOSABLE. In a decomposition, simplicity of its top would make
one summand semisimple and contained in B direct_sum N. A B summand is
forbidden by KO2's unique-normalized-trace retraction criterion under
noninvertibility. An N summand is forbidden by H2(N)=C and H2(R)=0.
This treats both possible simple summands of that socle.

The residual extension is precise. From 0 -> B -> M -> IC_D -> 0 and
H0(N)=H1(N)=0, the contravariant Ext sequence gives

    Ext1_W(M,N) = Ext1_W(IC_D,N).

The actual source extension 0 -> N -> R -> M -> 0 has connecting map
H1(M)=C -> H2(N)=C an isomorphism. Thus the remaining object is a
particular nonzero extension with a nonzero cohomological connecting map,
together with compatible COMMUTATIVE POLYNOMIAL SOURCE MULTIPLICATION.
No theorem rules out this extension, and no arbitrary extension is
asserted to carry that multiplication or to be a polynomial ring.

This applies in every degree with these hypotheses, including the natural
S_d and A_d doubly transitive actions when applicable. The surviving
one-boundary quartic S4 and degree-six A6/S6 packets would therefore have
an already SIMPLE trace-zero kernel. Its simplicity or semisimplicity
cannot be treated as an unsupplied extra recognition certificate for those
packets. The module extension and multiplication remain the real gate.

## 5. Boundaries, imports, and verification

All results use the actual plane normalization and the FULL Borisov
Stein-boundary input already retained in PF2BATT1. They do not apply the
surface argument to the three-dimensional Cox base-change order in QLD4.
They do not prove all cyclic endpoint lifts smooth, assume every omitted
prime ramified, or discard allowed unramified primes. Universal HC4 use
is through the retained Tang plane reduction.

The three new proved imports are the regular Riemann-Hilbert package,
finite/small pushforward of IC, and Borisov's normalized power-map statement.
The point-quotient trace proof, earlier finite controls, and source bytes
remain untouched. The [checker](verify_source_module_layers.py) supplies
exact Kummer, semisimple projection, rank and low-layer controls; it does
not verify the topological/category theorems or construct the required
extension. Independent replay, formal verification and external review
remain absent. Full HC4/JC2 remains open.
