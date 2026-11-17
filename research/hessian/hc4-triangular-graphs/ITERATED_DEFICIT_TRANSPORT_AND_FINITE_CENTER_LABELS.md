# Exact iteration of missing sheets and finite-boundary valuations

2026-11-16. HC4MYMH0. Actual complex polynomial plane Keller maps, ALL
generic degrees and ALL positive iterate numbers. No novelty priority
is claimed for composition of nonproperness sets or Euler pushforward.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda*Lambda) subset B, or full marked-energy integrality. OPEN.
External review/outreach remains DEFERRED.

MG0 locates every omitted finite-normalization prime outside the tight
dynamical tree. Here we retain its finite target center, its valuation
multiplicity, and the missing sheets at EVERY finite target point.
The result is an exact transport identity, not a vanishing theorem.

## 1. The deficit is the length of the entire omitted fiber

Let F:A2->A2 be Keller of generic degree d. For n>=1 put

    r_n(y)=#(F^-n(y)), delta_n(y)=d^n-r_n(y).

All these fibers are finite and reduced. The full finite normalization
pi_n:N_n->A2 for F^n is finite flat of degree d^n: N_n is a normal
surface, hence Cohen-Macaulay, and the regular target gives flatness
as in retained KU1. The ACTUAL source A2 is its Zariski-Main open.

Every retained fiber point is etale and has length one. Consequently

    delta_n(y)=sum_(z in (N_n-A2) intersect pi_n^-1(y))
                   length(O_(pi_n^-1(y),z)) >=0.          (1)

This includes ramified AND unramified omitted points, all residue
copies, special points and intersections of target images. There is
no replacement of scheme length by the number of omitted points.

The integer function delta_n is constructible. Its support is exactly
the nonproperness set A(F^n): it is the finite image of the omitted
closed boundary. A fiber has no omitted point precisely when it has
a finite neighborhood wholly in the source open.

## 2. An exact formula at every point for every iterate

For an integer constructible function phi on A2 define

    T phi(y)=sum_(F(x)=y) phi(x).

The sum is finite; no properness of F is assumed. In particular
T^n 1=r_n. Counting the first predecessor of each preimage gives

    delta_(n+1)=d^n delta_1+T delta_n,
    delta_n=sum_(j=0)^(n-1) d^(n-1-j) T^j delta_1.         (2)

For clarity, T delta_n=d^n T1-T^(n+1)1 and
d^n delta_1=d^(n+1)1-d^n T1; adding proves the recurrence.
Induction proves the second identity. These are identities of WHOLE
constructible functions, including the finite collision corrections.
They do not suppose that forward image curves are distinct.

Since every summand is nonnegative, (1)-(2) give the exact set identity

    A(F^n)=union_(j=0)^(n-1) F^j(A(F)).                    (3)

Images here are actual finite-point images. Their closedness in this
setting also follows from (1)-(2); alternatively use the finite
polynomial parametrizations of the nonproperness components.

## 3. Euler integration reproduces the budget, without amplification

Use compactly supported Euler characteristic. For constructible phi,

    integral T phi dchi_c = integral phi dchi_c.           (4)

This is ordinary Euler pushforward/Fubini. It does not require a
proper map: stratify into algebraic constructible pieces over which
the finite fiber counts are locally constant and use finite-cover
multiplicativity and additivity of chi_c. Equivalently use the usual
compactly supported constructible pushforward.

Both source and target have chi_c(A2)=1. Thus

    integral delta_n dchi_c=d^n-1,
    integral delta_1 dchi_c=d-1.

Integrating the full function identity (2) gives exactly

    (d-1)*sum_(j=0)^(n-1) d^(n-1-j)=d^n-1.               (5)

This is valid even with coincident forward components, new singular
image points and special omitted unramified points. Counting curve
components while dropping point corrections would not compute this
integral. Positivity of delta_n pointwise does NOT make its Euler
integral the sum of positive curve contributions.

Thus the iterated Euler budget itself adds no inequality to the
original budget: it transports it by (2). This does not assert that
every other geometric constraint on iterates is automatic.

The dropped point terms can be arbitrarily large even for a simple
normalization. For ALL m>=1, the known type of rational nodal curve

    gamma_m(t)=(t^2,t^(2m+1)-t),
    D_m={y^2=x(x^m-1)^2}

has normalization A1, with exactly m ordinary nodes at x^m=1,y=0.
The derivative of gamma_m never vanishes. Away from these nodes it
is one-to-one; each node has the two parameters t and -t. Their two
tangent vectors (2t,2m),(-2t,2m) have determinant 8mt!=0.
There are no other singularities: the implicit equation and its
derivatives reduce them to x^m=1,y=0; at x=0 the x-derivative is -1.
The equation is irreducible because x(x^m-1)^2 is not a square in C(x).
Consequently

    (gamma_m)_*1=1_(D_m)+sum_(nodes)1_node,
    chi_c(D_m)=1-m, integral (gamma_m)_*1 dchi_c=1.

This is an ALL-m curve control, not an ambient plane Keller map or
a hypothetical counterexample realization. It shows explicitly why
using only the image curve's Euler characteristic loses the budget.

## 4. Normalized A1 image curves are transported with degree ONE

We need an elementary curve fact. Let C be any closed irreducible
source curve with normalization A1. For a polynomial etale F, F(C)
is a closed curve with normalization A1, and the induced map of
normalizations is an ISOMORPHISM. The restriction on the singular
curves themselves need not be injective on points.

Indeed, if gamma:A1->C is its normalization, F composed with gamma
is a nonconstant polynomial parametrization. It is finite onto its
image: one nonconstant coordinate f(t) already makes t integral over
C[f], hence over the intermediate image coordinate ring. The image
is closed. Its normalized projective completion is dominated by P1,
so has genus zero. Its affine normalization has exactly one point at
infinity, since that point's inverse image can only be the single
infinity of the source P1. It is therefore A1.

The induced normalization map is a polynomial h:A1->A1. At every
finite parameter F is a local isomorphism of ambient surfaces, and
hence an isomorphism on the chosen curve-branch normalization germ.
So h' vanishes nowhere. A nonconstant complex polynomial with this
property has degree one. This proves the fact, with no global
curve rectification or polynomial inverse of F assumed.

The full nonproperness components have A1 normalization by the
retained boundary theorem PF2BATT1/KQ0, so this applies to EVERY
forward image F^j(C) of every such component.

In particular let D_n be the effective curve cycle whose coefficient
on a component of A(F^n) is the GENERIC value of delta_n there.
The curve-cycle pushforward counts normalization degrees; those
degrees are one in the present situation. Taking generic curve
coefficients in (2) gives

    D_n=sum_(j=0)^(n-1) d^(n-1-j) F_*^j D_1.             (6)

Coincident image components have their coefficients ADDED.
Equation(6) is not a substitute for (2) at singular/intersection
points. The point data remain in the full constructible function.

## 5. Positive-discrepancy labels survive after the first finite image

Keep any original omitted prime E in N_1-A2, its UNNORMALIZED
valuation v_E=ord_E, its image curve C, and its ramification index
e_E>=1. For valuations here use the honest restriction/pullback rule

    F_* v(H)=v(H composed with F),

even when the resulting valuation has center in the finite affine
plane. Do NOT replace it by the zero valuation when the normalized
escape degree is zero.

The definition of the ramification index at the target DVR gives

    F_* v_E=e_E ord_C.

Put C_j=F^j(C). The map F is etale at the generic point of each C_j;
the transverse ramification index of C_j over C_(j+1) is one.
Therefore, for ALL n>=1,

    (F^n)_* v_E=e_E ord_(C_(n-1)).                        (7)

The residue degree of E over C is not being set to one; it is
additional field information not represented by the scalar multiple
in (7). The degree-one statement of Section4 concerns F between the
successive finite target curves, not E->C.

For a prime divisor valuation v on a smooth surface use the
unnormalized log-discrepancy coefficient
a(v)=1+ord_v(source volume), extended homogeneously to multiples.
MG0 gives a(v_E)=e_E. A finite-plane curve has a(ord_C_j)=1.
Thus (7) preserves the positive coefficient e_E forever after the
first finite image. No iteration drives it to zero or into the
tight tree.

This does NOT bound the ramification indices of ALL omitted primes
of N_n. New omitted valuations can spend several steps at infinity
before their first finite image. They are not the fixed v_E followed
in (7), and they are retained in (1)-(2).

## 6. Stable image is not the vanishing of the deficit

Primary context: Peretz, Nguyen Van Chau, Gutierrez and Campbell,
*Iterated Images and the Plane Jacobian Conjecture*, inspected
[arXiv:math/0411245v2](https://arxiv.org/html/math/0411245v2).
Its complete pertinent Sections1-6, including the proofs of the
stabilization theorem and the equivalence in Theorem2, were read.
The [publisher record](https://www.aimsciences.org/article/doi/10.3934/dcds.2006.16.455)
identifies the published 2006 paper. No PDF was cached and no new
external theorem is used to prove (1)-(7).

That paper studies eventual equality of the image SETS. The missing
image of F^n is

    A2-F^n(A2)={y:r_n(y)=0}={y:delta_n(y)=d^n}.

It is not the positive-deficit locus {delta_n>0}. Stabilization of
image sets says nothing by itself about deficits on curves with
retained sheets. The paper's own Theorem2 identifies the additional
surjective-etale-on-a-cofinite-set injectivity statement as JC2,
rather than proving it automatically.

Do not infer from a missing point that its whole forward orbit
remains missing from every iterate, or assume it is periodic. Other
preimages can restore membership in an iterate's image.

## 7. Remaining gap

MH0 retains exactly the finite-center data that the normalized
tight-tree dynamics omits. Its deficit transport and Euler identity
are consistent in every degree and at every iterate. The positive
ramification coefficient is transported, not annihilated.

No contradiction to an actual nonzero deficit, no polynomial plane
counterexample, no trace-dual multiplication landing and no full
marked-energy integrality is obtained. Do not start an iterate-degree,
component-count or normalization-parameter census to prove a stronger
inequality by omitting the collision terms.

The checker supplies finite-function and symbolic local
transcription controls; (1)-(7) are written all-degree proofs, not
extrapolations. All previous proof/checker/context/receipt bytes are
sealed. No registered mathematical checker is rerun. Independent
replay, formal verification and external review remain false.
Full HC4/JC2 remains OPEN.
