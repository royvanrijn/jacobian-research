# Other generic fibres do not supply the retained jump blocks

The fresh/historic panel does not obtain its large jumps by importing the
rational generic sections of the seven tested family presentations through
isomorphisms of elliptic fibres. This is an exact exclusion of a specific
simultaneous block construction, using only equations.

There are **23 rational presentations in 112 target/family cells**. Every
one is the existing presentation or a globally base-equivalent coordinate
change. Each cell has one remaining irreducible parameter polynomial of
degree 23 or 24. In **all 112 remaining fields**, the isomorphism requires a
genuine quadratic extension. Its sign involution negates every transported
source section. Consequently the span of those sections and all their
conjugates has **zero rational non-torsion subspace**.

There is a uniform reason for the sign obstruction: its norm squareclass
is a fixed nonsquare resultant of the source family. Thus a new block by
this route requires an **additional factorization of the parameter
correspondence at t**. An intact residual preimage field cannot merely
switch its isomorphism class from nonsquare to square. The certified
irreducibility on the entire panel rules out that required event here.

Thus this proposed source of a large simultaneous rational block fails
on the highs themselves, as well as the controls. It cannot explain the
[strict classes outside the inherited global pool](CLASS_CREATION_REQUIRES_NEW_UNRAMIFIED_COVERS.md).
It does not exclude other unramified covers, other families, multisections
on the K3, or additional points in the same number fields.

## The proposed specialization event

Fix a target E_t: y^2=x^3+a(t)x+b(t), and another source family
E'_s: y^2=x^3+A(s)x+B(s). All target a,b here are nonzero and nonsingular.
An isomorphism of these short equations requires

\[
 H_t(s)=b(t)^2 A(s)^3-a(t)^3 B(s)^2=0.
\tag{1}
\]

At a smooth preimage, put

\[
 r_t(s)=\frac{b(t)A(s)}{a(t)B(s)}.
\]

Equation (1) implies r_t(s)^2 A(s)=a(t) and r_t(s)^3 B(s)=b(t).
The isomorphism is

\[
 (x,y)\longmapsto\bigl(r_t(s)x,\;u^3 y\bigr),
 \qquad u^2=r_t(s).
\tag{2}
\]

For an irreducible parameter factor q(s) of H_t, write
F_q=Q[s]/(q). This is the **parameter preimage field**, not the cubic
two-division field studied in the strict-class construction.

A precise necessary entry condition for a rational block made by tracing
source sections from this component is

\[
\boxed{r_t(s)\in F_q^{*2}.}
\tag{3}
\]

Equivalently there must be z(s)∈Q[s], deg z<deg q, such that

\[
 q(s)\mid a(t)B(s)z(s)^2-b(t)A(s).
\tag{4}
\]

This is an explicit splitting condition on a finite cover at t, with no
exceptional point supplied. If it holds, (2) transports source sections
over F_q and their traces are rational points. Those traces may still
vanish or be dependent/inherited. Neither (3) nor a factorization of H_t
by itself proves a rank jump.

The experiment asks whether the successful t values satisfy this gate
on any component beyond their existing rational presentation. They do not
within the frozen family set.

## Why a nonsquare excludes even combinations of conjugate sections

Suppose q is irreducible and r is not square in F_q. Then
M_q=F_q(u), u^2=r, is a quadratic field extension. For any section
P(s)∈E'(Q(s)), its transported specialization Q∈E_t(M_q) satisfies

\[
 \sigma Q=-Q,\qquad
 \operatorname{Tr}_{M_q/F_q}Q=O,
 \operatorname{Tr}_{M_q/\mathbb Q}Q=O,
\]

where sigma fixes s and sends u to −u. This remains true at a pole of a
section's coordinate functions, using specialization on the smooth proper
elliptic curve. No point search or height calculation is involved.

Take a finite Galois closure containing all these points, for all source
families under consideration, and tensor their generated group with Q.
Each generator has zero Galois average. Averaging is a projection onto
the invariant subspace in characteristic zero. Therefore the entire span
has zero invariant subspace. Relations between conjugates, or mixing
different preimage fields and different sections, cannot produce an
invariant non-torsion direction when every generator has zero average.

This is stronger than reporting a zero trace for a displayed basis. It
applies to **every rational generic section of the tested source family**
at that component, and to their conjugates. It does not apply to arbitrary
points of E_t(M_q) that are not obtained by this construction.

## Frozen panel and rational preimages

The [first protocol](SECOND_GENERIC_PRESENTATION_PROTOCOL.json) fixes the
sixteen existing target equations and these seven source presentations:
074d9, 103b2, 11952, a1-fibration-01, a1-fibration-05,
historic-lineage-074d9 and curve398-p16875. They are precisely the source
models already used for the generic-capacity work, not a newly selected
atlas. Arithmetic inputs contain their coefficients and the existing
own-family parameter as a positive control; no points or rank labels.

For each target/source pair the worker forms (1), checks infinity in the
degree-(8,12) chart, and removes the known control root when applicable.
It either gives a root-free reduction modulo a prime, or factors the
remaining rational polynomial. Every rational preimage is checked against
the square condition in (2). All sixteen own-family controls are recovered.
All target workers complete within their separate 15-second caps.

The seven additional rational presentations are exhausted by two identities:

- 074d9 and historic-lineage-074d9 are globally base-equivalent. This
  accounts for both compact 074d9 targets and historic curves 356 and 385.
- a1-fibration-01 and curve398-p16875 are globally base-equivalent. This
  accounts for the fresh MW16 high, its low control, and curve398.

The [independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_second_generic_presentation_verification_v1.json)
reconstructs the two rational Möbius maps from three matching parameters
and verifies the complete polynomial identities

\[
 (ct+d)^8 A_{\rm right}\!\left(\frac{at+b}{ct+d}\right)
   =\lambda^4 A_{\rm left}(t),
\]
\[
 (ct+d)^{12} B_{\rm right}\!\left(\frac{at+b}{ct+d}\right)
   =\lambda^6 B_{\rm left}(t),\qquad \lambda\in\mathbb Q^*.
\]

The matrices and scalings are stored exactly. Three matched points alone
would not prove these identities. No equality of specialized section
bases is assumed: sections transported through a global family isomorphism
already belong to the original family-level global pool L. Their images
cannot provide the large outside-pool block. Division or saturation can
change mod-two representatives, but cannot increase the free rank of their
rational span.

## The nonrational components and their sign obstruction

After removing those rational roots, every remaining parameter polynomial
is irreducible over Q: **23 cells have degree 23 and 89 have degree 24**.
Thus there are no hidden quadratic or other lower-degree parameter factors
in these residual polynomials.

The separate [sign-trace protocol](J_PREIMAGE_SIGN_TRACE_PROTOCOL.json)
checks the same 112 cells with a 30-second total worker cap and primes
at most 251. It constructs no number fields. At each full-degree squarefree
reduction it records the irreducible factor degrees. A proper rational
factor degree must be a subset sum at every such prime; their intersection
is empty in every cell.

For the sign certificate, let h be one irreducible factor modulo p. The
worker verifies that r is a unit and a nonsquare in F_p[s]/(h). This proves
it is nonsquare in F_q: if r had a square root in F_q, that root would be
integral and a unit at this unramified prime, and would reduce to a square
root modulo h. The monic, full-degree, squarefree reduction ensures the
local order is maximal, so no polynomial-index ambiguity is hidden here.

All 112 cells pass both tests; the largest prime used is 239. The
isomorphism fields therefore have degrees **46 or 48**, with the sign
involution required by the trace proof. All of the tested nonrational
source-section constructions have rational span zero.

The [CAS-free verifier](verify_j_preimage_sign_trace.py) checks **4871**
finite irreducible factors by Rabin tests, exact products, all subset-degree
exclusions, and the nonsquare Euler characters. It checks the scaling
residues by multiplication instead of using the worker's inversion routine.
The rational-preimage replay independently uses direct finite-root
enumeration instead of the worker's polynomial gcd test.

## The sign obstruction comes from a fixed family resultant

The local certificates have a uniform algebraic explanation. Suppose
deg A=8, deg B=12, H=b^2 A^3−a^3 B^2 has degree24 and leading coefficient
L, and R=Res(A,B)≠0. These hypotheses hold in every retained cell.
Resultant multiplicativity and evaluation at the roots of A and B give

\[
 \operatorname{Res}(H,A)=a^{24}R^2,\qquad
 \operatorname{Res}(H,B)=b^{24}R^3.
\]

For a polynomial C, the norm of C(s) in Q[s]/(H) equals
Res(H,C)/L^{deg C}. Therefore the isomorphism scaling satisfies

\[
\boxed{
 N_{\mathbb Q[s]/(H)\,/\,\mathbb Q}
 \left(\frac{bA(s)}{aB(s)}\right)
 =\left(\frac ba\right)^{24}
   \frac{a^{24}R^2/L^8}{b^{24}R^3/L^{12}}
 =\frac{L^4}{R}.
}
\tag{5}
\]

Its rational squareclass is **[R]**, independent of the target a,b.
The exact [resultants](../../artifacts/generated-results/elliptic-curves/rank_jump_resultant_norm_obstruction_v1.json)
are positive and nonsquare for all seven source presentations. This is
certified by integer square-root bounds, with no integer factorization.
A second [verifier](verify_resultant_norm_obstruction.py) recomputes all
seven through 20-by-20 integer Sylvester determinants and fraction-free
elimination, independently of the worker's Sage resultant operation.
It also verifies that the two pairs of globally equivalent presentations
have the same resultant squareclass.

Every removed rational factor in this panel has square scaling class.
Dividing its norm contribution out of (5) preserves the nonsquare norm.
The single remaining irreducible field therefore has nonsquare scaling
class: a square element would have square norm. This reproduces all 112
finite-place sign obstructions from just seven fixed source resultants
and the parameter-polynomial irreducibility certificates.

This sharpens the specialization criterion. On the full-degree open
locus of (5), after removing the inherited rational isomorphism components,
**a further factorization is necessary** for a new rational trace block.
Only then can a component acquire square scaling while another component
carries the nonsquare total norm. Factorization is not sufficient: its
component scaling must satisfy (3), and the resulting traces must be
nonzero and independent modulo the generic subgroup.

This is a concrete cover-splitting requirement rather than an invariant
being proposed to correlate with rank. It fails on the successful fibres
themselves, so it excludes this construction as their explanation. The
fixed resultant is not a rank predictor.

## Paired outcomes and what has been learned

Retained labels are joined only after the arithmetic. The full
[comparison artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_j_preimage_sign_trace_comparison_v1.json)
preserves all sixteen rows.

| Comparison | Retained gains | Rational presentations up to base equivalence | Rational span from nonrational source components |
|---|---:|---:|---:|
| 103b2, 3726/881 versus −1049/2296 | +10 / observed 0 | 1 / 1 | 0 / 0 |
| 11952, −2448/11 versus −1171/1683 | +10 / observed 0 | 1 / 1 | 0 / 0 |
| 11952, 110314/102227 versus 130349/28916 | +10 / observed 0 | 1 / 1 | 0 / 0 |
| MW16-01, −1867/270 versus −3187/3697 | +11 / observed +1 | 1 / 1 | 0 / 0 |
| Historic 356, 385, 398 | +12 / +12 / +14 | 1 / 1 / 1 | 0 / 0 / 0 |

The low labels remain censored. This calculation gives no upper bound on
their full ranks. It proves the proposed constructor contributes no new
outside-pool block even on the successful fibres; this is not a population
correlation or a claim that their full incidence structures coincide.

The lesson is to distinguish two possible mechanisms. A second rational
high-rank presentation would supply soluble sections immediately, but it
does not occur here beyond global duplicates. Algebraic presentations
also do not help automatically: their isomorphism signs can eliminate
the whole rational trace span. In these examples that elimination is exact.

The strongest remaining target is still a point-independent construction
of the specialized cubic field's additional unramified, S-split classes.
Ramified family covers or K3 multisections remain possible sources. The
current result excludes transporting rational generic sections through
these fibre-isomorphism correspondences; it does not exclude those other
constructions. The condition from t to new incidence, and its subsequent
simultaneous rational solubility, remain missing.

No positive prospective feature or scoring change follows. The result is
an **incidence exclusion for a specified rational-block constructor**.
Condition (3) concerns splitting of its auxiliary isomorphism cover;
it must not be promoted to evidence for a new Selmer or Mordell–Weil class
without nonzero, independent trace points. Nothing here measures search
visibility or modifies Agent1's search.

Replay:

```sh
timeout 30 sage -python elliptic-curves/rank-jump/verify_second_generic_presentation.py check
timeout 60 python3 elliptic-curves/rank-jump/verify_j_preimage_sign_trace.py check
timeout 30 python3 elliptic-curves/rank-jump/verify_resultant_norm_obstruction.py check
```

All scripts, inputs and outputs are new rank-jump-specific files. Only
the preceding rank-jump note receives a navigation link. No search protocol,
candidate population, worker limit, frozen search output or mathematical
status entry is changed.
