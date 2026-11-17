# Replica-equivariant factors descend in degree or require one mixed fibre

PF2RULED139. 16 November 2026. ALL actual complex polynomial plane
Keller degrees and polynomial source-coordinate frames. This combines
the NEW shared HC4MYMR0 replica theorem with PF2RULED130's retained
universal-cover plane-chart construction. It gives an actual smaller
plane Keller map, or a necessary mixed-fibre profile; it does not close
the full Wright/HC4/JC2 or finite trace-dual implication.

## 1. The actual equivariant source image and polynomial fibre data

Let B=C[P,Q] subset R=C[x,y], Jac(P,Q)=1, F=Frac(B), L=Frac(R),
d=[L:F]. Let F subsetneq E subsetneq L be an embedded intermediate
field stable under delta=p(y)*partial_x, p nonzero, k=deg p. Use the
ACTUAL source image

    V=Spec(R intersection E), m=[L:E]>1, b=[E:F]>1.            (1)

We deliberately exclude E=F here: retained KX0 already settles descent
to the target field. HC4MYMR0 gives V smooth AFFINE and etale over the
original target, a SURJECTIVE etale X=A2->V, and an affine A1-ruling
u=h(y) with

    deg h=m, E=C(h(y),s), s=x/p(y)+beta(y),
    p=h'*r(h), r a nonzero POLYNOMIAL,
    k+1=m*(deg r+1).                                         (2)

These are shared retained theorems, credited rather than re-proved.
The polynomial factorization and its finite-pole proof in MR0 Section3
are essential: dividing p inside an unknown E is not allowed.

For a in C write

    h(y)-a=c*product_(i=1)^q(y-c_i)^n_i,
    g_a=gcd(n_1,...,n_q).                                    (3)

MR0's surjectivity and etaleness say that EVERY fibre component of V
is an A1 reached isomorphically by at least one WHOLE source line.
Its multiplicity is precisely n_i for every line assigned to it.
Every root occurs and every component is reached. Consequently

    gcd of V-fibre multiplicities = g_a;
    some component has multiplicity g_a
      iff some n_i=g_a.                                    (4)

Repeated roots may map to the same component; (4) does not count them
as distinct components. This uses MR0's exact positive source-line
counts, not a guessed fibre passport or an intermediate plane model.

## 2. Only one fibre can have no simple polynomial root

For a fibre of (3), the total zero order it contributes to h' is

    sum_i(n_i-1)=m-q.                                     (5)

Different fibres contribute at disjoint points. If all n_i>=2, then
q<=m/2, so (5) is at least m/2. Two such fibres would contribute at
least m zeros to a derivative of degree m-1, impossible. Thus

    AT MOST ONE h-fibre has all root multiplicities>=2.     (6)

This is an all-degree derivative count, not a classification of
polynomial covers. In particular at most one g_a>1. The GENERAL
affine A1-ruling fundamental-group presentation reused in PF130
Section3 now gives

    pi1(V)=C_g,
    g=the unique g_a>1 if one exists, and g=1 otherwise.    (7)

Here V is the source-image open, and g is its fundamental-group order.
It is NOT automatically the order of a subgroup of the original generic
deck group. Universal-cover automorphisms need not lift to L.
The imported general ruling presentation and fibre structure retain
their original primary-source/proof-sketch assurance limits.

## 3. A genuine degree descent when every fibre has a gcd component

Suppose EVERY fibre in (3) has some n_i=g_a. By (4) the actual
V-ruling has a component of its WHOLE gcd multiplicity in every fibre.
Apply precisely the chart construction in PF130 Section3: normalize
the degree-g Kummer base cover at the one whole multiple fibre, with
g=1 meaning no cover. Connectedness of the generic A1 and codimension-one
unramifiedness plus purity give the UNIVERSAL finite etale cover

    W->V of degree g.

Over the multiple fibre, each component of multiplicity g*a_j lifts
to g components of multiplicity a_j. The hypothesis supplies a_j=1;
every other fibre already has a unit component. Choose one unit A1
in each exceptional W-fibre and delete the other components. The
retained affine-line bundle input makes the remaining open A2.

This construction uses only the actual smooth affine ruled V and its
proved group (7). It does NOT require V to be the ENTIRE etale locus
of the finite intermediate normalization. We do NOT transfer PF130
Section2's residual-defect equalities, its identification of U_N, or
finiteness of the original relative map to this selected open.

Compose the resulting plane chart with W->V and V->A2_target. This is
an ACTUAL polynomial plane Keller map, of generic degree

    d_new=b*g=d*g/m>1.                                    (8)

If g<m this is a strictly smaller-degree counterexample. If g=m,
(3) at the unique whole multiple fibre has only one root: its degree
is m, all multiplicities are divisible by m, and at least one exists.
Thus

    h(y)=A*(y-c)^m+a, A!=0.                               (9)

In this case L/E is GALOIS CYCLIC of order m, even if p itself is not
monomial. Indeed multiply s in (2) by r(h)=p/h', which belongs to E,
and rescale to obtain t=x/(y-c)^(m-1)+b_0(y), b_0 rational. Then
E=C((y-c)^m,t), L=C(y-c,t), and rotating y-c gives C_m. Its elements
fix F and are actual generic rational decks. This reuses MR0's slice
argument AFTER (9); it does not assume p is monomial.

Therefore in a hypothetical counterexample of MINIMAL GENERIC DEGREE,
every proper equivariant factor satisfying the gcd-component condition
is the cyclic alternative (9). More generally (8) gives an explicit
degree-reduction mechanism without assuming minimality.

## 4. The remaining mixed alternative in minimal degree

Suppose the original counterexample has minimal generic degree and
(9) does not hold. Section3 forces some fibre a_* to have NO root
of its gcd multiplicity. Let g=g_(a_*). Then

    n_i/g>=2 for EVERY root of h-a_*, gcd_i(n_i/g)=1.       (10)

It is wholly critical and is the unique such fibre by (6); all other
fibres have a SIMPLE root, hence a unit V-component. The value g in
(10) agrees with (7), including g=1. Over the algebraically closed C,
absorb the leading coefficient and write

    h-a_*=H(y)^g,
    ell=deg H=m/g=sum_i(n_i/g)>=5.                        (11)

The inequality is elementary and UNBOUNDED. There must be at least
two parts in (10), since a single part would have gcd equal to itself.
Two or more integers at least2 with gcd1 have sum at least5; equality
is precisely the partition2+3. Therefore

    m>=5*g, k+1=m*(deg r+1)>=5*g.                        (12)

Some normalized part in (10) is at least3: if all were2 their gcd
would be2. At its corresponding root, h' has order n_i-1>=3g-1.
Since p=h'*r(h) is a polynomial product, we also have

    maximum finite root multiplicity of p >= 3g-1.        (13)

This is a root-multiplicity constraint, not merely a total-degree bound.

The universal W->V need not have a plane chart in this case. We do
NOT prove absence of every chart from the fibre partition alone, or
construct the necessary profile as an actual Keller map. A different
proved plane chart would still give (8) and contradict minimality.

This closes ALL Morse polynomial h of degree>=3 as equivariant factors
of a minimal counterexample: a Morse critical value has one double
root and m-2 simple roots, while ordinary values are reduced. Here
Morse means simple critical points with DISTINCT critical values.
For m=2 the polynomial is already of form (9). More generally EVERY
non-power h whose every fibre has a root of its gcd multiplicity is
removed in all degrees, not merely a bounded list of derivative degrees.

As a direct corollary, replica degree k<=3 permits only the cyclic
alternative in minimal generic degree. The normalized mixed case starts
at ell=5, with partition2+3. This is a threshold consequence of (12),
not a proposal for an increasing low-degree coefficient campaign.

## 5. An all-degree square-free replica exclusion

Suppose p is SQUARE-FREE, with arbitrary degree. Equation(13) rules out
the mixed alternative for every g>=1. Equivalently h' divides p, so all
roots of h' are simple. Every h-fibre then has root multiplicities only1
or2: if a simple root occurs it equals the gcd1, and otherwise every
root has multiplicity2 and equals the gcd2. Thus the chart condition
holds in every fibre even when critical VALUES collide; a Morse
hypothesis is unnecessary here.

Minimality forces (9). At its centre c, p=h'*r(h) has a root of order
at least m-1. Square-freeness gives m<=2, and properness gives

    m=2, L/E an ACTUAL quadratic generic deck extension,
    p(y)=(y-c)*q((y-c)^2), q a nonzero polynomial.         (14)

The harmless constants from (9) are absorbed into q. In particular
deg p is ODD and p is odd about c: p(c+t)=-p(c-t).
We claim necessity, not existence of this quadratic field for such p.

Consequently in a hypothetical MINIMAL-degree counterexample:

* EVERY square-free replica p of EVEN degree preserves NO proper
  intermediate field, in ANY polynomial source-coordinate frame.
* A square-free odd-degree replica can preserve a proper intermediate
  only through (14). If the full actual generic deck group has odd
  order (including trivial), it preserves NO proper intermediate field
  for ANY square-free replica, of ANY degree.

This excludes entire replica/factor classes, not Keller maps merely
because such a source derivation exists. Existence of a suitable
equivariant intermediate is a hypothesis. The first mixed threshold
p has a repeated root, exactly as (13) requires.

For MR0's finite source-derivative field E_delta=Frac(J_delta), the same
minimal-degree hypotheses therefore give [L:E_delta]<=2 for EVERY
square-free p, and E_delta=L for even deg p or odd-order generic deck.
The possibility E_delta=F is already excluded by retained target-field
LND descent. These remain FIELD conclusions; no raw derivative-algebra
normality, finite B-lattice, full trace regularity or inverse is inferred.

## 6. What remains and the assurance boundary

The surviving mixed fibre is an ACTUAL SOURCE-IMAGE requirement coupling
the polynomial h, exact component multiplicities, and the replica
factorization p=h'*r(h). It is stronger than an arbitrary exact-volume
surface control. PF130's x-dependent exact mixed(2,3) surface and PF138's
all-r free-additive controls have NO such supplied polynomial source
cover and etale target; they do not realize (1)--(2) here.

MP0's cyclic-DECK theorem does not eliminate a non-Galois L/E with
polynomial relative monodromy. No noncyclic deck survivor is reopened.
The source image V need not be the entire finite-intermediate etale
locus; the original full finite normalization, all ramified primes,
residual discriminant and finite trace-dual remain unchanged.

The next useful question is incompatibility of the mixed profile
(10)--(12) with the ACTUAL finite target order or a chart in its
universal cover. It is not further canonical-volume residue counting:
PF138 already completes that test when its whole-atlas hypotheses hold.
The cyclic alternative, trivial-deck cases, non-equivariant factors
and the full trace-dual/marked-energy implication remain OPEN.

No new external theorem or classification is imported. The complete
MR0 and PF130 proofs are read and credited, with their original source
limits. Small exact derivative/partition/degree controls support this
written all-degree argument; no predecessor mathematical checker is
executed. No novelty, formal verification, independent whole-proof
replay, external review or full-goal closure is asserted.
