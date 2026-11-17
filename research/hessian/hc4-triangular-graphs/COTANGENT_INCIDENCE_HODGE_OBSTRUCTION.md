# Cotangent incidence sources: an all-degree Hodge obstruction

Date: 2026-11-16. Claims HC4MYLI0 and HC4MYLI1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the full marked-energy alternative is also OPEN. External review/outreach
remains DEFERRED. The preceding goal turn was PROGRESS (LH0/1).

This note tests a concrete counterexample route, rather than assuming
trace regularity. It closes the ENTIRE primitive-beta moment-chart source
question left open in the historical natural-torus reduction note. A
general incidence calculation applies in every dimension and degree.
It does not exclude every primitive chart, construct a Hessian potential,
or resolve HC4. In particular the other primitive gamma chart passes this
Hodge test and is not declared polynomial affine space.

## 1. LI0: the general all-degree incidence formula

Let F=(F_1,...,F_n):A^n->A^n be any complex polynomial Keller map, n>=2,
of generic field degree d. Choose an affine target coordinate frame and
write the last coordinate as F_n. For every c in C set

    X_c={F_n(x)+sum_(i=1)^(n-1) q_i*F_i(x)=c}
                         subset A^n_x times A^(n-1)_q,
    Z={F_1=...=F_(n-1)=0} subset A^n,
    Z_c={F_1=...=F_(n-1)=0,F_n=c}.

The target incidence hypersurface eliminates its last base coordinate and
is A^(2n-2). The induced map X_c->A^(2n-2) is its base change under F,
hence ETALE. In particular X_c is smooth affine of dimension2n-2.
One may also check smoothness directly: its source gradient is
DF^T*(q_1,...,q_(n-1),1), which cannot vanish.

Z is a smooth affine curve, possibly DISCONNECTED or EMPTY. Each Z_c is
a finite reduced set, possibly empty; write N_c=#Z_c. Projection of X_c
to A^n has fibre A^(n-2) off Z, no fibre over Z minus Z_c, and fibre
A^(n-1) over Z_c. Off Z this is a Zariski locally trivial affine bundle:
on F_i!=0 solve for q_i. Thus in the Grothendieck ring of varieties

    [X_c]=L^(2n-2)-L^(n-2)*[Z]+N_c*L^(n-1).          (1)

Here L=[A1]. NO CANCELLATION of L in that Grothendieck ring is used.
It can have zero divisors. We instead apply the compactly supported
Hodge--Deligne polynomial E, valued in the integral domain Z[u,v], with
E(A1)=uv, denoted L again in the formulas below.

Suppose Z has k connected components; let g be the SUM of their completion
genera and s the TOTAL number of punctures in their smooth projective
completions. Each nonempty affine component has at least one puncture,
so s>=k. For the empty curve put k=g=s=0. The standard smooth-curve Hodge
formula is E(Z)=k*uv-g*(u+v)+k-s. Consequently

    E(X_c)=L^(2n-2)+(N_c-k)*L^(n-1)
                   +g*L^(n-2)*(u+v)+(s-k)*L^(n-2).   (2)

These different bidegrees cannot cancel. Thus E(X_c)=L^(2n-2) iff
g=0,s=k,N_c=k. The first two conditions say that EVERY component of Z
is A1. On each such component F_n is a nowhere-critical regular function,
because F is etale. It is a polynomial on A1 with nowhere-zero derivative,
hence LINEAR and an isomorphism to the target axis. This automatically
gives N_c=k for EVERY c.

We obtain an exact all-degree dichotomy:

    E(X_c)=E(A^(2n-2)) for SOME c
      iff Z is a disjoint union of affine lines
      iff E(X_c)=E(A^(2n-2)) for EVERY c.             (3)

The empty union is allowed. If even ONE component of Z is not A1, then
for EVERY c and EVERY r>=0,

    X_c times A^r is NOT isomorphic to A^(2n-2+r).    (4)

Indeed its Hodge difference from affine space remains nonzero after
multiplication by (uv)^r. Thus polynomial changes of coordinates and
ordinary affine stabilization cannot repair that source. Equality of
Hodge polynomials in the other branch is NECESSARY, not sufficient, for
an affine-space isomorphism.

### Generic affine target lines exclude the incidence descent in all degrees

If d>1, for a GENERAL affine target line the preimage Z is geometrically
irreducible and the restricted map Z->A1 has generic degree d. Hence Z
cannot be A1: an etale polynomial map A1->A1 has degree1. Equations(3)--(4)
then exclude every level and every ordinary affine stabilization.

We explicitly IMPORT the Bertini irreducibility input used here: Stacks
Lemma37.32.3, tag0G4C, including its hypotheses, read with its complete
proof. Apply it successively to the linear system spanned by1 and the
remaining target affine coordinates. The base locus is empty. At a
source point, two independent affine target equations through its image
pull back, by etaleness, to a nonempty codimension-two intersection;
the section1 has empty zero divisor. This verifies the lemma's required
three-section condition at every stage of dimension at least2. Iterated
general hyperplanes give the stated irreducible preimage of a general
line. Choose the line also meeting the dense open where F is finite
etale of rank d, so the restricted generic degree really is d.
This is a precise classical input, not independently reproved by the
finite checker. There is NO assertion that EVERY special target line
has irreducible preimage.

The statement concerns incidence source recognition. An etale incidence
map alone is not a gradient map, and a polynomial affine-space source
alone would not identify its symplectic form with standard Darboux
coordinates in dimensions greater than two.

### The sharper plane consequence

For n=2, the target axis is a divisor, so retained KMFI1 guarantees Z is
NONEMPTY. If E(X_c)=L^2, (3) supplies an embedded A1 component of Z on
which F is injective. By the classical Abhyankar--Moh line embedding
theorem, followed by the already retained Gwozdziewicz one-line criterion,
F is a polynomial automorphism. This exact closed-component consequence
is already retained as HC4MYGR0 in CRITICAL_COMPONENTS.md, including its
rectification source record. Conversely for an
automorphism every X_c is A2 by target coordinates. Therefore

    F plane Keller inverse
      iff E(X_c)=L^2 for ANY ONE level c
      iff X_c is A2 for ANY ONE level c.              (5)

If F were noninvertible, NO level of this one-helper incidence family
could be made affine two-space, even stably. The classical line inputs
are credited; this does not prove that the equality in (5) holds
automatically. It is an obstruction to another source-reconstruction
shortcut, not a new unconditional JC2 theorem.

## 2. LI1: the historical primitive-beta source is never affine four-space

Use the KNOWN Alpoge foundational map F=(A,B,C), already retained in
KO2/KW1, with u=xy:

    A=(1+u)^3*z+y^2*(1+u)*(4+3u),
    B=y+3*x*(1+u)^2*z+3*x*y^2*(4+3u),
    C=2*x-3*x^2*y-x^3*z.

The historical note NATURAL_TORUS_HAMILTONIAN_REDUCTION.md gives the
primitive beta!=0 chart and gauge beta=1. Its source at moment level c is

    X_c={B+2*alpha*A-gamma*C=c} subset A5.             (6)

The target is standard A4 after eliminating B. Equation(6) is exactly
the incidence family of Section1 with coefficient pair(2A,-C), last
coordinate B. Its unresolved polynomial source question in that old
note is now CLOSED NEGATIVELY at EVERY c.

Compute Z={A=C=0}. At x=0 it is

    L0={(0,t,-4t^2):t in C}, B=t,

an affine line. On x!=0 put u=xy,v=x^2*z. Then

    C=x*(2-3u-v),
    x^2*A=(1+u)*((1+u)^2*v+u^2*(4+3u)).

C=0 gives v=2-3u, reducing x^2*A to (1+u)*(u+2).
The other two disjoint closed components are therefore

    L1={(t,-1/t,5/t^2):t in C*}, B=2/t,
    L2={(t,-2/t,8/t^2):t in C*}, B=-2/t.

They are closed in A3 because xy=-1 or -2 already forces x invertible.
This is the COMPLETE set decomposition, not a selected open subset:

    Z=A1 disjoint_union Gm disjoint_union Gm,
    [Z]=3L-2,
    N_c=1 when c=0, and N_c=3 when c!=0.             (7)

Substituting into (1) yields

    [X_0]=L^4-2L^2+2L,
    [X_c]=L^4+2L                 (c!=0).             (8)

The same formulas hold for E(X_c;u,v), with L=uv. The L-term has
coefficient2 in BOTH cases, so no level has the affine-space Hodge
polynomial. In particular

    X_c times A^r != A^(4+r), EVERY c and EVERY r>=0. (9)

Euler characteristic alone would miss the zero level:
chi(X_0)=1, whereas chi(X_c)=3 for c!=0. The Hodge polynomial proves
the uniform exclusion, including the exceptional zero level. This closes
the entire beta-chart source question without coefficient/degree searches,
quantization, or a source cancellation conjecture.

The old archived note and its checker are retained BYTE-FOR-BYTE; they
are historical evidence, not rerun. The known Jacobian/collision facts
are reused rather than replayed. This result does not alter the older
rank-two symplectic construction or any Dixmier status.

## 3. A different primitive chart passes the test: no overstatement

The target cotangent gamma coordinate has primitive weight+1. Gauging
gamma=1 gives another source family

    Y_c={2*alpha*A+beta*B-C=c}.

Its relevant axis preimage is Z'={A=B=0}. On x=0 the only point is the
origin. On x!=0 define H=(1+u)^2*v+u^2*(4+3u). Then
x^2*A=(1+u)H and x*B=u+3H. Thus B=0 gives H=-u/3,
and A=0 forces u=0 or -1. The case u=-1 has H=1 and B!=0,
so is impossible; u=0 forces v=0. Consequently

    Z'={(x,0,0):x in C}=A1, C=2x.

Every level point has one preimage on this axis, and (1) gives

    [Y_c]=L^4, E(Y_c)=L^4 for EVERY c.              (10)

This is a sharp boundary of the invariant, NOT a proof Y_c is A4.
No polynomial trivialization or compatible Darboux coordinates are
constructed. The present theorem does not discard all primitive momentum
choices, nonlinear charts, or every dimension-reduction architecture.
It also does not claim that passing the Hodge test gives an HC4 potential.

## 4. Provenance and assurance

One explicit proved external input is adopted: Stacks Lemma37.32.3 for
generic irreducibility. Its full proof was read; it is not independently
reproved or formalized. Hodge--Deligne additivity, affine-bundle
multiplicativity and the smooth-curve formula are retained standard
foundations already used in the archived oriented-cubic motivic theorem
and KT0. The classical Abhyankar--Moh and Gwozdziewicz line inputs are
credited as retained foundations, not claimed new classifications.

The finite checker controls the new component identities, their complete
case elimination, fibre formulas and Hodge defects, and binds the retained
proof/source bytes. It does not rerun the old map Jacobian, collision or
natural-torus checkers. The general all-n/all-degree/all-c/all-stabilization
proof is written above, not inferred from sampled dimensions or levels.
No formal verification, independent whole-proof replay or external review
is claimed. Full HC4/JC2, actual trace-pairing regularity and full marked
energy integrality remain OPEN. No HC4 counterexample is constructed.
