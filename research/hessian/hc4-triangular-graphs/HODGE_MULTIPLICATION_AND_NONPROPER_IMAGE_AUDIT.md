# Hodge multiplication and nonproper direct image: route audit

16 November 2026. ROUTE AUDIT ONLY. No new HC4 claim is registered.
The question is whether purity, polarizations, or the lowest Hodge piece
can supply the missing ordinary trace pairing. They do not supply it by
the proposed formal arguments below. This does not rule out every possible
use of Hodge theory. Full HC4/JC2 remains OPEN.

## 1. First retain the multiplication theorem already proved

For an ACTUAL complex plane Keller source, use the notation of
[KQ0](SOURCE_MODULE_LAYERS_AND_TRACE_ZERO_KERNEL.md) and
[KU2](CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md):

    B=C[P,Q], R=C[x,y], S=full finite normalization,
    Lambda=trace-dual of S, I=Soc_(D_B)(R), A_can=B[Lambda].

The retained all-degree theorem already says

    D_B*S=I,       image(Sym^2_B I -> R)=A_can.

A_can fills ALL omitted unramified primes. It need not equal R. In
particular the socle I is generally NOT an algebra. Reproving its generic
quadratic span or its semisimplicity would not add a missing constraint.
KR2 already distinguishes holonomic self-duality from an ordinary
B-valued field-trace pairing. The calculation below tests whether adding
Hodge-theoretic purity repairs precisely that distinction.

## 2. Pure characters still multiply to a polar integer character

Work transversely at a generic omitted ramified divisor, after the same
faithfully flat tame splitting used in KU2. Write

    O=C[[u]], K=C((u)), S_e=C[[t]], t^e=u, e>=2.

An extra smooth coordinate can be included for a surface. In this factor,
the retained IC module and codifferent have the elementary forms

    I_e=O*1 plus direct_sum_(a=1)^(e-1) K*t^a,
    Lambda_e=t^(1-e)*S_e.

For every 1<=a<e the two elements t^(-a),t^(-(e-a)) belong to Lambda_e
and to the nontrivial IC character summands. Their product and trace are

    t^(-a)*t^(-(e-a))=u^(-1),
    Tr_(C((t))/K)(u^(-1))=e/u.

Each nontrivial Kummer character comes from a finite unitary local system
of Hodge type(0,0); its minimal extension is pure. Purity of the two
factors therefore does NOT make their ordinary product remain in the
regular integer-character lattice O. This is the same codimension-one
defect already present in KR2, now screened against Hodge purity.

There is also a direct tensor calculation, without any Hodge convention.
For M_a=K*t^(-a), multiplication induces an isomorphism of differential
modules

    M_a tensor_O M_(e-a) -> K,       v tensor w -> v*w.

The result is the localized trivial connection, rather than its minimal
extension O. Thus one cannot replace the ordinary O-module tensor product
in the actual multiplication map by the tensor product of generic local
systems followed by minimal extension. Extension across the divisor and
ordinary multiplication have not commuted. A theorem about tensor products
of variations on the complement does not establish that commutation.

This local model is NOT an actual plane Keller counterexample. In an
actual source its appearance is conditional on an omitted ramified prime;
it is precisely what a global plane-specific argument would have to
exclude. No residue-degree-one assumption on the original prime is used.
Unramified omitted primes and finite-point issues remain as in KU2.

## 3. The lowest Hodge piece allows the same product

The primary source screen is Henry Dakin,
[*Weight filtration and generating level*, arXiv:2503.14216v2](https://arxiv.org/html/2503.14216v2),
Theorem3.24 in Section3.3. Its smooth one-variable specialization gives
F_0 M(u^(-alpha))=O*u^(-alpha) for 0<alpha<1; these nontrivial
characters are already minimal extensions. Hence both factors above
occur in their respective lowest Hodge pieces. Their product still
has a pole. The theorem statement and the pertinent explicit proof were
inspected; its general V-filtration machinery was NOT independently
reproduced. This is a source screen, not a new imported HC4 theorem.

Even the trivial-character localization exhibits the indexing issue.
Mustaţă--Popa,
[*Hodge ideals*, arXiv:1605.08088v4](https://arxiv.org/html/1605.08088v4),
Sections3,4,6 and Proposition8.2, give, for a smooth divisor u=0,

    F_k O(*D)=u^(-(k+1))*O,       k>=0.

Consequently its actual ordinary products satisfy

    F_p O(*D)*F_q O(*D)=F_(p+q+1) O(*D),

including F_0*F_0=F_1, not F_0. This equality follows just by multiplying
the displayed principal modules. A degree-zero multiplicative filtration
on source functions cannot be assumed to be this filtration. Their
projective strict direct-image statement also retains its projectivity
hypothesis. The relevant formulas, Proposition8.2's proof, and the
localization/direct-image discussion were read; the full Hodge-module
foundations were NOT independently audited.

No global identity F_0 I=Lambda, no complete weight-filtration description
of R, and no new Hodge generation-level theorem is needed or adopted here.
The transverse countercalculation suffices to reject the proposed formal
purity-to-product inference. It does not settle a plane-specific global
Hodge identity or supply ordinary trace integrality.

## 4. Relative dimension zero does not transport a trivial filtration

Let j:G_m -> A1 be the open immersion, a smooth etale morphism of relative
dimension zero. Its source structure sheaf has the trivial source Hodge
filtration F_0 O_(G_m)=O_(G_m). The unfiltered direct image is
O_(A1)(*0). Its canonical Hodge filtration is the one in Section3:

    F_0(j_* O_(G_m))=u^(-1)*C[u].

Simply pushing forward the source degree-zero piece instead yields

    j_*(F_0 O_(G_m))=C[u,u^(-1)],

which is not even a finite C[u]-module. These two modules differ, for
example on u^(-2). Thus an etale relative-dimension-zero calculation of
the UNFILTERED direct image does not identify its canonical Hodge
filtration with the naive pushforward filtration. The localization at
the missing boundary is essential. This example has a non-plane source;
it refutes the unrestricted formal inference, not JC2.

For the actual map A2 -> A2 one may still identify the unfiltered
D-module direct image with R carrying its Nambu derivations, as in KQ0.
It is not valid to append 'therefore every source polynomial is in the
canonical F_0 of the target direct image'. That additional assertion
needs its own argument through the nonproper boundary. Nor does purity
of the source constant Hodge module make its nonproper direct image pure;
projective direct-image theorems have different hypotheses.

## 5. What changes in the research direction

Retire automatic purity-to-ordinary-product closure, polarization-to-
regular-field-trace, and naive nonproper filtered-pushforward shortcuts.
The local obstruction occurs in codimension one, so normal Hartogs
extension cannot remove it. KU2 already assembles the quadratic image
globally; another additive IC decomposition or a generic character tensor
decomposition is not the missing algebra theorem.

Revisit Hodge theory only with an independently justified constraint that
uses the ACTUAL polynomial plane source to exclude these polar products,
or with a theorem whose hypotheses on the nonproper direct image have
been established. No existence of such a theorem is asserted here.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION:

    ACTUAL polynomial plane source
       ==> Tr(Lambda*Lambda) subset B.               OPEN

Full marked-energy integrality is also OPEN. No new mathematical checker,
registered replay, primary cache, theorem registration, independent replay,
formal verification, or external review is supplied by this audit.
External review and outreach remain DEFERRED. Prior theorem and checker
bytes are preserved. Published generic-degree-four plane exclusion remains
retained; this audit neither creates an intermediate plane nor reconstructs
a missing artifact.
