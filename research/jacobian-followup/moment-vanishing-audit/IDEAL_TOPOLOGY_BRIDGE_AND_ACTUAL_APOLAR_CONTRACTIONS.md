# Fixed ideal topologies and the actual apolar contraction

16 November 2026. Internal written argument and small exact transcription
checks. No independent replay, formal verification or external review.
The full actual-gap implication and both complete exceptional MZ questions
remain **OPEN**. This is an adoption check for a proposed universal proof
mechanism, not another positive subclass or normal-form programme.

## 1. The published theorem and its actual hypothesis

Irena Swanson, *Linear equivalence of ideal topologies*, Math. Z.234
(2000),755-775, [DOI](https://doi.org/10.1007/s002090050007),
[author PDF](https://www.math.purdue.edu/~iswanso/topology.pdf).
The cached author manuscript has21 physical pages. Its text encoding is
garbled; physical pages1,2,11,12,13 were rendered and visually read.
The precise general statement and its reduction proof are on pages12-13.
This audit does not claim to reprove the deep valuation inputs to Theorem3.1.

Main Theorem3.3, as stated in the author manuscript, concerns a
**commutative Noetherian local ring R and two fixed ideals I,J**.
Among its equivalent conditions are:

* the filtration `I^n:J^infinity` defines the same topology as the powers
  `I^n`;
* some integer k satisfies `I^(kn):J^infinity subset I^n` for every n>=1.

Here saturation means the union of `I^n:J^a` over a>=0. The prime
symbolic-power statement in the abstract is a familiar consequence.
Finite generation of the symbolic or saturated Rees algebra is **not**
a hypothesis. It is therefore incorrect to reject this lead merely
because a symbolic Rees algebra has not been proved finitely generated.

The missing adoption step is an EXACT model of the actual deficit
comparison by the appropriate powers/saturations of a FIXED ideal.
Nesting alone does not provide that step. The following control excludes
the most immediate proposed model, even after contraction to the two
variables for the fixed operator and the passive suspension derivative.
No Swanson theorem is applied to apolar heights in this note.

## 2. One fixed pure safe pair, with a genuinely passive suspension

Work over any characteristic-zero field K. Start with

    P=x^2+v,   T=partial_x^2 partial_v,

and adjoin a variable u ignored by T. Put

    F=P+u=x^2+v+u,   eta=partial_u.

The full ideals Ann(F^m) are nested: differentiating F^(m+1) by u
recovers (m+1)F^m, and every constant-coefficient operator commutes
with partial_u. Their contraction below is therefore nested as well.

The pair is PURE and SAFE with exact least real mixed cutoff1.
The operator T^m has order3m and QF^m has degree at most degQ+2m,
so m>degQ forces vanishing. At the boundary Q=v^m,

    T^m(v^m F^m)=(2m)! m! != 0.                         (1)

Its actual height is h_m=floor(m/2): first applying partial_v^k
leaves a nonzero scalar times F^(m-k), whose x-degree is2(m-k).
Applying partial_x^(2k) is nonzero exactly when2k<=m. Thus the
deficit is ceil(m/2), with a perfectly positive gap. No failure of
individual eventuality-to-linearity is present here.

Let A=K[[L,E]]. Define its actual annihilator contractions

    K_m={a(L,E): a(T,eta) F^m=0}.                       (2)

These are ideals, since all participating operators commute. Formal
evaluation is locally finite. The map L->X^2 V,E->U into the ambient
formal symbol ring is injective. Consequently (2) is a literal pullback
of the actual embedded apolar ideal, not an abstract height sequence.

## 3. Exact all-index product obstruction

Define the Hermite polynomials by

    H_q(X,E)=q! sum_(j=0)^floor(q/2)
               (-1)^j X^(q-2j) E^j / ((q-2j)! j!).     (3)

For f=x^2+y the translation identity gives

    exp(t partial_x-t^2 partial_y) f^m=(f+2xt)^m.

Taking the coefficient of t^q proves

    H_q(partial_x,partial_y) f^(q-1)=0.                 (4)

For the square use two independent parameters:

    exp((t+s)partial_x-(t^2+s^2)partial_y) f^q
      =(f+2x(t+s)+2ts)^q.

Its coefficient at t^q s^q is2^q: total degree2q requires selecting
2ts from every one of the q factors. Therefore

    H_q(partial_x,partial_y)^2 f^q=2^q(q!)^2 != 0.      (5)

For n>=1 now define a polynomial in the SAME two fixed variables L,E:

    R_n(L,E)=(2n)! sum_(j=0)^n
                 (-1)^j L^(n-j) E^(2j)
                   / ((2n-2j)! j!).                    (6)

On every F^m, partial_v and partial_u agree; this equality survives
all constant-coefficient differentiations. Thus (6) acts exactly as

    eta^n H_(2n)(partial_x,eta).

Equations(4)-(5) and differentiation by u give the exact identities

    R_n(T,eta) F^(3n-1)=0,
    R_n(T,eta)^2 F^(4n)=2^(2n)(4n)!(2n)! != 0.         (7)

Hence R_n belongs to K_(3n-1), but R_n^2 does not belong to K_(4n).
For ANY fixed nonnegative integer c, choose n with2n>=c+2. Then
6n-2-c>=4n. Nesting and (7) imply

    K_(3n-1) K_(3n-1) is NOT contained in K_(6n-2-c).   (8)

In particular there is **no fixed additive index error c** for the
product containment K_a K_b subset K_(a+b-c), even for this ONE fixed
actual nested PURE SAFE pair and this fixed two-parameter contraction.

The small checker tests (3)-(7), the literal passive suspension,
several exact heights and the sharp marker. The proof of EVERY n and
EVERY proposed c is (3)-(8), not the finite replay.

## 4. What this excludes, and the full target it leaves

Powers of a fixed ideal are multiplicative filtrations. So are symbolic
powers and saturated powers `I^m:J^infinity`: if J^a r is in I^m and
J^b s is in I^n, then J^(a+b)rs is in I^(m+n). Pullback along any
ring map preserves this product property.

Consequently the literal family K_m cannot be identified, index for
index, with the contraction of such a filtration in any commutative
Noetherian ring. A fixed additive shift of the indices cannot repair
the problem, because it would imply a product law with a fixed
additive error, contradicted by (8).

This does **not** exclude a different exact ideal model of the deficit
comparison, a comparison with genuinely different linear constants,
or another use of Swanson's theorem. In particular it does not assert
that all actual apolar filtrations lack useful linear comparisons.
The very pair above already has exact linear heights and cutoff1.

R_n varies with n. Thus (8) does **not** refute bounded-error
quasi-subadditivity for the powers of ONE fixed symbol, which remains
a separate unproved possibility. Nor does it construct an actual
eventual-but-nonlinear pair. An arbitrary nested or multiplicative
abstract filtration remains insufficient evidence for that purpose.

The owned full-theorem obligation stays unchanged: for EVERY actual
pure polynomial-phase/formal-symbol pair prove b_m->infinity implies
inf_m b_m/m>0, or realize a genuine counterexample. If this published
topology theorem is used, its fixed ideal and the exact correspondence
to actual deficits must first be supplied. The [scalar-equation
countercontrol](UNIVERSAL_REVERSED_TIME_EQUATION_COUNTERCONTROL.md)
does not address that different possible bridge.

The exceptional quotient programme already has the actual finite
four-variable radial linear cutoffs. Extra formal-symbol finiteness
does not prove universal intrinsic source splitting for
its actual invariant period module. Both full exceptional cases and
the full actual gap remain OPEN; review and outreach stay deferred.

See [checker](verify_apolar_topology_bridge.py),
[receipt](apolar_topology_bridge_verification.json) and
[source audit](apolar_topology_bridge_source_audit.json).
