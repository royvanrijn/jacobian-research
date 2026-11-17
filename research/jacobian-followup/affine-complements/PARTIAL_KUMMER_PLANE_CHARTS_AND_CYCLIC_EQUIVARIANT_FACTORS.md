# Partial Kummer plane charts remove the mixed equivariant alternative

PF2RULED142. 16 November 2026. ALL complex polynomial plane Keller degrees.
The degree descent is unconditional on minimality; its cyclic-factor
conclusion is for a hypothetical counterexample of MINIMAL GENERIC DEGREE.
The full Wright/HC4/JC2 and finite trace-dual implications remain OPEN.

The new step is to allow a finite base cover which is ramified on some
components, then retain an open plane on which the map is etale. PF130/139
required an entire finite etale universal cover. That stronger requirement
left a mixed fibre unnecessarily as a live minimal-counterexample case.
The original finite normalization, its boundary and trace-dual are retained;
the new plane defines a DIFFERENT comparison Keller map.

## 1. A component multiplicity gives an etale plane chart

Let V be a smooth connected AFFINE complex surface with an affine
A1-fibration rho:V->A1_u. Assume every fibre other than possibly the one
at a has a component of multiplicity ONE. Choose ANY component C of the
fibre at a, of multiplicity n>=1. Then there exists a dominant ETALE
morphism

    j:A2 -> V, generic degree n.                           (1)

It need not be finite, surjective, or a chart in a finite etale cover
of all of V. No volume, fundamental-group or Keller assumption is needed
for this geometric lemma.

Use the retained affine-line-fibration structure in HC4MYML0/PF130:
fibres have finitely many disjoint smooth A1 components, and an affine
A1-fibration with every fibre irreducible and reduced is a Zariski
A1-bundle. The precise primary reference is Dubouloz2002.12838v1,
Definition1, Proposition3 and the paragraph after Definition4. Its stated
proof-sketch/adopted-foundation limitations are unchanged. The primary
statements were rechecked; no stronger classification is imported.

Normalize V times_(A1_u) A1_t, where

    u=a+t^n,

to obtain a finite map q:W->V. The normalization is finite because the
rings are finitely generated over C. It is connected and has degree n:
the generic fibre of rho is A1 over C(u), so C(V)=C(u)(s), and adjoining
t with t^n=u-a has degree n even after adjoining s. The induced generic
fibre over C(t) is again A1. For n=1 take W=V directly.

We check smoothness and etaleness at EVERY point of the special fibre,
rather than invoking purity for the entire q. Near a component B of
rho^-1(a) of multiplicity mu, a local equation is

    u-a=eta*z^mu, eta a unit.

After an etale localization adjoining an n-th root of eta, the root
equation is t0^n=z^mu. Put g=gcd(n,mu). Its normalization has g branches
with local parameters v and a transverse parameter along B, and

    z=v^(n/g), t0=zeta*v^(mu/g).                         (2)

It is smooth. Normalization commutes with this etale localization, so
W is smooth near the whole fibre. Away from t=0 it is an etale base
change of smooth V, hence W is smooth everywhere.

The reduced inverse image of B is a finite ETALE cover of B=A1 of
degree g. Indeed (2) has no tangential branching at any point; unit-root
localizations are etale. Every finite etale cover of complex A1 is
trivial. Thus it consists of g disjoint A1s, each mapping isomorphically
to B, with multiplicity mu/g in the t-fibre and normal index n/g for q.
In particular, above the chosen mu=n component C there are n disjoint
REDUCED A1s, and q is etale along each ENTIRE line. Above any B with
n not dividing mu, q is ramified; none of those components is retained.

Choose one of the reduced lines above C. For every other degenerate
fibre of W->A1_t, choose one multiplicity-one A1: the cover is unramified
on the base there, and the corresponding V-fibre has such a component
by hypothesis. Delete ALL other components of these finitely many
fibres, including every other component at t=0. Call the result O.
The deleted curves are an effective Cartier divisor in smooth W, so
O->W is an affine open immersion. Since W is affine, O is affine.
The morphism O->A1_t is affine and flat (torsion-free over the base
local DVR), with unchanged generic A1 and EVERY closed fibre one
reduced A1. Proposition3 gives an A1-bundle. Pic(A1)=H1(A1,O)=0 makes
this bundle A2. The restriction q|O is etale: off t=0 q already is
etale, and the retained special line has normal index one by (2).
An open restriction does not change its generic degree n. This proves (1).

The deletion is essential. For a fibre 2*C2+3*C3 and n=2, q is ramified
over C3. The proof NEVER calls W->V finite etale and NEVER identifies W
with the universal cover. It deletes the whole inverse image of C3.

## 2. An explicit check on the retained mixed surface

PF130's existing control is

    S={X^2*Z=(Y^2-X)*((Y-1)^3-X)}.

It is smooth with unit volume omega=dX wedge dY/X^2, interpreted by
its regular residue expression at X=0. Its fibre at zero is
2*A1 disjoint union 3*A1. The following NEW map is polynomial:

    X=t^2,
    Y=t+t^3*w,
    Z=(2*w+t^2*w^2)*((t+t^3*w-1)^3-t^2).                (3)

Substitution gives the defining equation, and on t!=0

    j^*omega=2*dt wedge dw.                              (4)

Both sides are globally regular, so (4) holds everywhere. The map is
etale on all A2, has generic degree TWO (t^2=X and w=(Y-t)/t^3), and
maps t=0 isomorphically onto the multiplicity-two line, Z=-2*w.
It omits the multiplicity-three line. This is an actual etale plane
map to the retained surface, not an etale target map S->A2.

There is also an all-parameter version, without a coefficient search.
For any r,s>=2 put S_(r,s)={X^2*Z=(Y^r-X)*((Y-1)^s-X)}. It is smooth
by the same derivative test at X=0 and has the same unit residue volume.
Set X=t^r, Y=t+t^(r+1)*w, and

    Z=[((1+t^r*w)^r-1)/t^r]*((Y-1)^s-t^r).              (5)

The bracket is polynomial by the binomial theorem. Equation(5) gives
an everywhere etale plane map of generic degree r, with volume pullback
r*dt wedge dw. Interchanging the roots gives one of degree s.
No exactness or finite plane target is asserted for all S_(r,s).
The existing exact primitive of S=S_(2,3) is retained, not recomputed.

## 3. Apply the chart to EVERY replica-equivariant factor

Let B=C[P,Q] subset R=C[x,y], Jac(P,Q)=1, F=Frac B, L=Frac R,
d=[L:F]. Retain an ACTUAL proper field factor

    F subsetneq E subsetneq L, delta(E) subset E,
    delta=p(y)*partial_x, p nonzero,
    m=[L:E]>1, b=[E:F]>1, d=b*m.                         (6)

HC4MYMR0, with its original source/classification limits, supplies

    V=Spec(R intersection E), X=A2 -> V SURJECTIVE ETALE,
    V -> A2_target ETALE, u=h(y), deg h=m,
    E=C(h,s), s=x/p(y)+beta(y), beta RATIONAL,
    p=h'*r(h), r a nonzero POLYNOMIAL.                   (7)

Its whole source-line isomorphisms give: every root of h-a is assigned
to a V-fibre component with EXACTLY its root multiplicity; every
component is reached. No beta specialization or entire-etale-locus
replacement is used.

Recall PF139's elementary all-degree count. A fibre with q roots and
all multiplicities at least two contributes m-q>=m/2 zeros to h'.
Two such fibres contradict deg h'=m-1. There is therefore at most ONE
V-fibre with no component of multiplicity one.

If none exists, Section1 with n=1 gives an etale plane map A2->V,
and composing with V->A2_target gives an ACTUAL polynomial plane
Keller map of degree b=d/m<d.

Otherwise choose in the unique such fibre a root of MINIMUM multiplicity
n, and its actual V-component. Section1 gives an actual comparison
Keller map of degree

    d_new=b*n=d*n/m.                                    (8)

If h is NOT a shifted pure power, that fibre has at least two roots,
so n<=floor(m/2). Thus in EVERY non-power case there is a smaller
comparison counterexample with

    1<d_new<=d/2.                                       (9)

These maps are noninvertible because b>1. The composition is polynomial
since its source is A2; etaleness makes its Jacobian a nowhere-zero
polynomial, hence a nonzero constant. No entire-cover surjectivity of
the new chart, or equality with the original source, is needed.

## 4. Minimal-degree consequence: all these factors are actual cyclic decks

For a hypothetical counterexample of MINIMAL GENERIC DEGREE, (9) is
impossible. Thus EVERY proper factor (6), for EVERY nonzero source LND
in EVERY polynomial coordinate frame, has

    h=A*(y-c)^m+a, A!=0,
    p=(y-c)^(m-1)*r0((y-c)^m), r0 a nonzero polynomial.  (10)

Retained MR0/PF139's rational-slice argument now gives an ACTUAL
Galois cyclic extension, without assuming p is monomial. Multiply the
slice s by r(h)=p/h' inside E, and rescale, to write

    E=C((y-c)^m,z), L=C(y-c,z),
    z=x/(y-c)^(m-1)+gamma(y), gamma rational.             (11)

Rotation of y-c, fixing z, gives m relative automorphisms. They fix F,
so E is the fixed field of a subgroup C_m of the actual generic rational
deck group D. In particular m divides |D| and deg p+1. MP0 already
proves D cyclic for actual plane Keller maps. Consequently there is at
most one such E for each divisor m of |D|, across all source replicas.
This classifies these factors, not full monodromy or arbitrary surface pi1.

For MR0's finite source-derivative field E_delta=Frac(C[delta^j P,
delta^j Q]), (10)--(11) apply whenever E_delta is proper in L. E_delta=F
is excluded by retained KX0. In minimal generic degree with TRIVIAL D,

    E_delta=L for EVERY nonzero source LND,              (12)

not only bare translations, monomial or square-free replicas. With
nontrivial D, every proper E_delta is an actual cyclic deck fixed field.
This is a FIELD statement: no raw-algebra equality, finite B-module
landing, or stability of the full S or Lambda follows.

## 5. What is closed, and what remains

PF139's mixed normalized fibres, INCLUDING every coprime mixed profile
with all parts at least two, are no longer possible equivariant factors
of a MINIMAL-degree counterexample. The old square-free restriction is
unnecessary for the cyclic-factor conclusion. A 2+N fibre gives a
comparison map of degree 2*b instead of (N+2)*b, independently of N
and of the rational slice. In HC4MYMU0's actual degree-five top factor,
the comparison degree is 2*b instead of 5*b. Its finite-order defects
and marked trace-pole statements remain valid conditional statements
for arbitrary original degree; they are no longer a live MINIMAL-degree
mixed-family obstruction to solve.

The general chart lemma also applies to other actual affine ruled opens
with at most one fibre lacking a unit component. It does NOT apply to
arbitrary complete-base rulings or to several such bad fibres without
an appropriate base cover. The partial cover is not the original finite
normalization and supplies no new residual-defect equality for it.

Source-equivariance itself is NOT automatic. Cyclic/trivial-deck
finite-target realization, non-equivariant factors, original ramified
boundary, full trace-dual multiplication/marked-energy integrality and
full Wright/HC4/JC2 remain OPEN. All original primes are retained.

The full MR0, PF130, PF139, ML0 and MU0 proofs were read. The primary
Dubouloz Definition1/Proposition3/fibre paragraph were rechecked, retaining
their already adopted assurance limits. The local normalization (2),
whole-component deletion and degree computation are written universal
proofs. New cheap exact controls test their arithmetic, the explicit
mixed plane chart, and the scope boundaries, not a parameter census or
surface classification. No predecessor mathematical checker is replayed.
No new external theorem, novelty, formal verification, independent whole
proof, external review or full-goal completion is claimed.
