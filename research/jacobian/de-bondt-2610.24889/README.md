# De Bondt 2610.24889v1: exact audit and constructive follow-up

Date: 13 November 2026. Scope: the user's requested examination of
[Michiel de Bondt's preprint](https://arxiv.org/abs/2610.24889), with bounded
attempts to contribute. No new dimension record is established here.

The concrete 12-, 18- and 90-variable constructions pass the checks below.
The follow-up adds an explicit conjugacy with our earlier map, a direct
integer matrix certificate, a simple lifting proof, and nine scoped
graph-deletion obstructions. General minimality remains unknown.

The exact reproduction and graph obstruction are registered as BCRDB1 and
BCRDBGRAPH1 in the live [mathematical-status authority](../../MATH_STATUS.json).
Existing archived claim scopes have not been edited or reactivated. The repository-wide status
audit already failed before this work because claim `KDSQ6` names the absent
consumer `artifacts/generated-results/README.md`. Missing historical evidence
was not reconstructed to bypass that check.

## Provenance and reproduction

The [arXiv v1 source](https://arxiv.org/src/2610.24889v1) is retained as
[`source-v1.tar`](source-v1.tar), SHA-256
`0e48fed1f0804ad263c24e7f4b146eb8d0ff618bf22a8f975907829b5cfe49ba`.
Its extracted manuscript and ancillary worksheets are byte-preserved in
`source/`. No Maxima worksheet was executed. Python implementations transcribe
and check the relevant identities independently of Maxima.

From the repository root:

```sh
.venv/bin/python research/jacobian/de-bondt-2610.24889/run_audit.py
# Optional repetition of the specified finite discovery attempts:
.venv/bin/python research/jacobian/de-bondt-2610.24889/run_audit.py --bounded-search
```

Each worker has a 60-second cap. The graph audit takes approximately 15 seconds
in this environment; the remaining exact checks are short. The optional
grouping/basis trials take approximately 11 seconds. SymPy and python-flint
are used by the producers and graph RREF calculation. The independent checker
uses only the Python standard library and the archived Fraction sparse engine.

[`receipt.json`](receipt.json) records input and result digests.
[`independent.json`](independent.json) records the independent replay boundary.
No formal verification, external review, or independent replay of the full
Vanishing Conjecture differential-operator argument is claimed.

## What was checked

| Object | Exact verification |
|---|---|
| Section 4.6 cubic map | Reconstruct extension and elementary target map; degrees `3,3,3,3,3,2,2,2,2,2,2,2`; cubic-output rank 5; rational collision |
| Its Jacobian | Independent full 12-by-12 sparse polynomial determinant equals 1; 438 memoized determinant states |
| Homogeneous companion | Explicit 18-variable map, homogeneous cubic nonlinear part, collision, and written Schur-complement determinant identity |
| Section 4.5 and our K12 | Twelve exact polynomial conjugacy identities with an explicit invertible monomial matrix |
| Power decomposition | Reconstruct 104 ungrouped cubes and 90 grouped cubes; independently check all 18 coefficient identities, rank 18 of B and rank 90 of the selected cubes |
| Cubic-linear lift | Explicit integer 90-by-90 matrix, two distinct rational points, equal images, and determinant proof by exact factorization and Sylvester's identity |

See [`verify.py`](verify.py), [`equivalence.py`](equivalence.py),
[`polarize.py`](polarize.py), [`druzkowski.py`](druzkowski.py), and
[`check_independent.py`](check_independent.py).

## An exact identification of the earlier construction

Let K be our [archived twelve-variable map](../../foundations/verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md),
and J the six-cubic map in the paper's section 4.5. Define

\[
Tx=(x_1,2x_2,2x_3,4x_4,4x_6,2x_5,4x_{10},2x_8,
12x_{11},4x_{12},2x_9,-x_7/2).
\]

Exact substitution gives

\[
K(Tx)=T J(x),\qquad \det T=-49152.
\]

Thus these two six-cubic maps are conjugate by a permutation and rational
coordinate scalings. This sharpens the qualitative identification into a
replayable certificate. It does not identify the new five-cubic section 4.6
map with K. The archived K is reconstructed from its pinned F13 source;
[`equivalence.json`](equivalence.json) retains the transformation.

## A direct integer cubic-linear lift

This is a simple specialization of the pairing idea, not a claim to a new
general reduction theorem. It avoids square matrix completion and the need
for a right inverse over the integers.

Write the five-cubic map as K(x)=x+Q(x)+C(x), with C supported in its first
five output coordinates. In coordinates (x,w,tau), let

\[
H=(\tau Q(x)-\tau^2\iota w,\ C_1(x),\ldots,C_5(x),\ 0),
\qquad G_\lambda=I+\lambda H.
\]

Here iota includes five coordinates into twelve. The Schur complement of the
five companion coordinates gives

\[
\det DG_\lambda
=\det\bigl(I+\lambda\tau DQ(x)+\lambda^2\tau^2 DC(x)\bigr)
=\det DK(\lambda\tau x)=1.
\]

No division by lambda or tau is needed in this polynomial identity. This
also proves that JH is nilpotent. For any nonzero rational lambda and any
original collision point p, the lift

\[
\widehat p_\lambda=(p,-\lambda C(p),1/\lambda)
\]

satisfies

\[
G_\lambda(\widehat p_\lambda)=(K(p),0,1/\lambda).
\]

This provides rational collisions for every such lambda directly; no square
root scaling or general specialization argument is required.

For the integral variant in worksheet `section4d.mac`, let L have the 90
linear forms as rows. Multiplying the first 17 coefficient rows by six gives
an integer matrix B such that

\[
B(Lz)^{\circ3}=6H(z)+e_{18}z_{18}^3.
\]

The list contains the pure form z18. If its index is i, then
L_i=e18^T and row 18 of B is e_i^T. Let L0 be L with row i replaced by zero,
and set

\[
A=L_0B\in\mathbb Z^{90\times90},\qquad
F(y)=y+(Ay)^{\circ3}.
\]

The other entries in column i of B need not vanish: t18 cubed is shared by
the polarization formulas. Set b=B e_i-e18; its last entry is zero. Then

\[
BF(y)=S(G_6(By)),\qquad S(z)=z-bz_{18}^3.
\]

S is an elementary target shear with determinant one. This shear must be
retained; omitting it would give an incorrect projection identity, although
the collision and determinant conclusions below survive. Write Gtilde=S G6.
Sylvester's determinant identity therefore gives

\[
\det DF(y)
=\det\bigl(I_{18}+3B\operatorname{diag}((L_0By)^2)L_0\bigr)
=\det D\widetilde G(By)=1.
\]

B has rank 18 over Q, so choose any rational right inverse C0. Given two
distinct collision points p,q of G6, define

\[
u=C_0p,\quad v_0=C_0q,\quad
v=F(u)-(Av_0)^{\circ3}.
\]

For collision points with equal last coordinate, S preserves their common
image, so Gtilde(p)=Gtilde(q). One has Bv=q and Av=Av0, hence F(v)=F(u);
u and v are distinct because
their B-images differ. This proves noninjectivity and supplies explicit
coordinates without an integral right inverse. The saved certificate uses
the displayed rational G6 collisions with tau=1/6.

[`druzkowski90.json`](druzkowski90.json) contains A, B, both rational points
and their common image. Exact multiplication gives the additional structural
information

\[
\operatorname{rank} A=18,\quad
\operatorname{rank} A^2=2,\quad
\operatorname{rank} A^3=1,\quad A^4=0.
\]

These ranks concern this particular certificate. Nilpotence of A alone is
not being used as a substitute for the Jacobian determinant proof.

## A new bounded obstruction for the five-cubic map

Apply the repository's old graph-deletion/target-completion method to the
new map, rather than repeating its earlier calculations on K.

For each literal triangular component j=4,...,12, impose K_j=0 and delete
that coordinate. Test a one-stage elementary correction of a retained bad
component using every monomial of target degree at most three in the other
ten retained raw outputs. This includes all linear mixing and has 285
columns. The test is an exact rational linear system on terms of ordinary
degree greater than three.

| Deleted j | Obstructed component i | Coefficient rank | Augmented rank |
|---:|---:|---:|---:|
| 4 | 2 | 276 | 277 |
| 5 | 3 | 276 | 277 |
| 6 | 1 | 275 | 276 |
| 7 | 1 | 277 | 278 |
| 8 | 2 | 276 | 277 |
| 9 | 3 | 275 | 276 |
| 10 | 2 | 276 | 277 |
| 11 | 3 | 275 | 276 |
| 12 | 2 | 275 | 276 |

Thus none of these nine literal deletions admits the specified correction
stage. The calculation uses exact rational RREF, not an inference from a
finite-field miss. [`graph_audit.json`](graph_audit.json) retains both
degree-two and degree-three checks.

This does not exclude nonliteral graph coordinates, source changes,
sequential nonlinear target stages, or a different starting construction.
It proves neither an 18-variable homogeneous lower bound nor a
90-variable cubic-linear lower bound.

## Bounded improvement attempts

[`bounded_search.json`](bounded_search.json) records:

- all 512 subsets of the nine published regrouping operations, with the
  fixed seven-cube polarization: minimum 90;
- sixteen specified invertible triangular rational basis changes, with
  five seeded greedy polarization restarts each: minimum 95;
- [`polarization.json`](polarization.json) also records thirty seeded
  restarts allowing the four-cube and signed seven-cube identities on
  the fully grouped expression: best 90.

The grouping enumeration is exhaustive only for those nine binary choices.
The greedy trials do not exhaust decompositions. No improved dimension
was found; no minimality conclusion follows.

## Corrections and limits in the written preprint

1. **Origin-fixing assumption.** Theorem 2.2 should explicitly require
   F(0)=0 if “linear part is the identity” means its degree-one homogeneous
   part. For F(x)=x+1, n=1, d=2, its stated dimension formula gives m=1,
   so there are no auxiliary coordinates. Its preimage identity requires
   q=tp. At p=0 this would force t^-1 H(0)=t, impossible for a fixed
   polynomial H. The worksheet's normalization check already excludes
   this example. The concrete Keller constructions are origin-fixing.
   An assertion of exact degree d also needs a nontrivial nonlinear part;
   the identity map is a degenerate exception.

2. **Inverse missing in the normalization formula.** In the proof of
   Theorem 3.2(3), the displayed Phi=A(X+t(H-H(0))) has derivative A^2 at
   zero, not the asserted A^-1 A. Replacing the left A by A^-1 fixes this
   step. For the nilpotent linear example H=(x2,0), its printed formula
   has derivative [[1,2t],[0,1]], directly demonstrating the mismatch.
   This is a local proof correction, not a counterexample to the concrete
   constructions or a verdict on the full arbitrary-ring theorem.

3. **Independence is not global optimality.** Ninety independent chosen
   cubes, each used by B, prevent deletion from that fixed list. They do
   not exclude a different list with fewer powers. Phrases suggesting an
   absolute structural limit should be restricted to the specified
   construction or replaced by an upper-bound statement.

The scalar-extension and arbitrary-ring claims, all intermediate-degree
tracks, and the full Vanishing Conjecture operator argument have not been
completely audited here. No messages have been sent to the author.

## What would be worth trying next

The most direct open improvement target is a different simultaneous cube
decomposition of the same 18-output tensor, with fewer than the89 forms now certified in the [shared-basis follow-up](../cubic-reduction/POWER89.md).
The new rank-five cubic map is also a better starting point than the archived
rank-six map for nonliteral source-coordinate/target-completion work. The
nine literal deletions above should not be rerun unchanged. An improvement
requires an exact reconstruction and collision/determinant certificate;
neither a greedy miss nor a flattening bound settles the remaining cases.


## Coordinate-support follow-up

The [36-term integral conjugate](../cubic-reduction/PROOF12.md), KCUB3,
simplifies the displayed support of the existing12-variable map from40
to36 nonlinear occurrences. Two elementary source/target conjugations
align repeated carrier combinations, and diagonal conjugation clears
denominators. Its determinant and integer collision are independently
replayed. This changes neither the12/18/90 dimension bounds nor any
minimality conclusion.

## Shared-cube89-variable follow-up

The [exact89-form basis](../cubic-reduction/POWER89.md), KCUB5, improves
the audited cubic-linear construction from90 to89 variables. It replaces
the fixed cube list with a different shared basis, allows temporary
support growth during discovery, and retains the necessary padding
target shear. Its integer matrix, determinant factorization and complete
rational collision have an independent exact replay. The original
90-variable theorem and fixed-list independence result remain valid;
neither is a global minimum statement.

For the36-term coordinate frame, [KCUB4](../cubic-reduction/GRAPH_DELETION_OBSTRUCTION.md)
also strengthens the fixed-graph target-repair test to arbitrary polynomial
degree by explicit finite dual certificates. General source and target
sequences remain outside these construction-specific obstructions.
