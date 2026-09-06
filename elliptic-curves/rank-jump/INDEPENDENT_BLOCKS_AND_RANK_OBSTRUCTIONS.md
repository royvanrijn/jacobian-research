# Adding independent blocks and finding efficient rank obstructions

This extends [J1–J6](RANK_JUMP_MECHANISM_THEOREMS.md) with four theorems
that turn a proposed construction into either a certified contribution or
a useful exclusion. The general statements have written proofs below;
the code verifies exact finite applications, not universal formal proofs.

The intervening [triple panel](LOW_DEGREE_TRIPLE_RELATIONS_CONCENTRATE_ON_THE_PLUS4_CONTROL.md)
completed the former first search priority: all 11 degree-six/eight
incidences occur at the +4 control. The
[coverage audit](RESEARCH_PIVOT_AFTER_THE_CARRIER_PANEL.md) excludes the
displayed native dictionary as an explanation of the whole retained
quotient at 62 of 69 positive published-R17 addresses. These results
motivate the present direction: certify independent additions, and target
the part of a Cassels–Tate calculation that can change a rank bound.

The current operational rules are [SEARCH_THEOREM_GATES_V2.json](SEARCH_THEOREM_GATES_V2.json).
The original gate file and J1–J6 note remain immutable replay inputs;
their proposed triple task is now completed, not a current instruction.
No new search, descent, Selmer upper bound, or original-curve rank result
is asserted here.

## J7. A lower certificate that does not mistake saturation for rank

Let A=E(K) for a number field K, let M⊂A have certified rank r, and put
L=M+<P₁,…,Pₙ>. Let t_ell=dim_Fell A[ell], and let

\[
\phi:A\longrightarrow\mathbf F_\ell^N
\]

be a **verified homomorphism**, for example concatenated good-reduction
Kummer characters, with ell prime. Write a=dim φ(M), b=dim φ(L).
Suppose a matrix R of verified rational relations among the Pᵢ modulo
M⊗Q has rank c over Q.

**Theorem J7 (finite-signature and relation sandwich).** If
q=rank L−r, then

\[
\boxed{\max(0,b-r-t_\ell)\le q\le n-c.} \tag{J7}
\]

Equivalently, the observed signature increment b−a yields only

\[
q\ge \max(0,(b-a)-(r+t_\ell-a)).
\]

Matching bounds give the exact rank of this **specified subgroup quotient**,
without a full Mordell–Weil basis or a Selmer computation.

**Proof.** The homomorphism on L factors through L/ell L, of dimension
rank L+dim L[ell]≤r+q+t_ell. Hence b≤r+q+t_ell. The relation upper bound
is J2. Subtract a for the second form. ∎

The defect budget r+t_ell−a is nonnegative by the same argument on M.
It accounts for both undiscovered finite signatures and possible
ell-saturation/torsion effects. If a=r+t_ell, no correction is needed.

The correction is necessary. In A=Z, M=2Z and P=1, reduction modulo two
gives a=0,b=1,r=1,t₂=0. The signature increment is one but q=0.
In A=Z⊕Z/2, a torsion point can increase the signature span without
increasing rational rank. Conversely P=2e₂ modulo M=<e₁> in A=Z²
has zero new mod-two signature but genuine quotient rank one.
These are finitely generated group counterexamples to logical shortcuts,
not new elliptic-curve examples.

**Search use.** A positive result must contain exact point membership,
the named generic subgroup and its rank, a genuine finite homomorphism,
and a torsion-dimension proof or upper bound. A table of arbitrary bits
does not satisfy J7. Never add ranks measured over different prime fields.
A failed lower-bound test leaves q UNKNOWN unless other evidence closes it.

For each retained quartet, concatenate the existing verified Kummer
characters of the 17 generic points and the four native points. Their
ranks are a=17 and b=20, with E(Q)[2]=0. A single exact relation gives
n−c=3. Thus J7 independently closes q=3 using only these 21 points and
one relation; the full rank-24/rank-25 witness basis is unnecessary for
this quotient certificate. It remains useful for measuring the unexplained
four/five directions.

## J8. Combining constructions: intersections are the missing term

Work in (E_b(K)⊗Q)/(M⊗Q). Let A and B now denote the rational spans of
two finite point blocks (the letter A in J7 denoted the whole group).

**Theorem J8 (block addition and diminishing increments).**

\[
\boxed{\dim(A+B)=\dim A+\dim B-\dim(A\cap B).} \tag{J8}
\]

For finite label sets I⊂J and another block B,

\[
\dim(A_I+B)-\dim A_I\ \ge\
\dim(A_J+B)-\dim A_J.
\]

In particular block ranks can be added exactly when their intersection
is zero. For several blocks the correct rank is the sum of their
successive certified increments, not the sum of their individual ranks.

**Proof.** The kernel of A⊕B→A+B is {(v,−v):v∈A∩B}.
Rank–nullity proves J8. The increment from B is
dim B−dim(B∩A_I); enlarging A_I enlarges this intersection. ∎

These are vector-space identities applied to **specialized** points.
Disjoint cover labels or distinct generic quadratic characters do not
prove a zero specialized intersection. A rational branch change
Pᵢ↦wᵢ−Pᵢ preserves its quotient line.

Exact coordinates in an independently certified witness basis compute
every term of J8. Alternatively apply J7 with the *already accumulated
subgroup* as M. Its new defect budget must be recomputed: saturation
certified for an earlier M need not hold after adjoining another point.

A practical certificate for full complementarity is a combined finite
signature rank r+q_A+q_B+t_ell (with the assumptions of J7), together
with exact individual ranks q_A,q_B. It proves q_union≥q_A+q_B,
and subadditivity gives equality. One can instead use a smaller
torsion-corrected bound if torsion in the combined subgroup is known.

**Search use.** Freeze the candidate blocks before consulting exceptional
oracle points, then prove their increments after specialization.
A block contributing zero is redundant relative to the current span;
another representative might still improve search coordinates, but it
adds no rank. The replay tests all three pair-versus-pair partitions in
each quartet, with exact specialized overlap dimensions.

## J9. A radical partner certifies two obstruction dimensions

Let S be a finite-dimensional F₂-space with an alternating bilinear form
λ. In the arithmetic application S=Sel₂(E/K) and λ is the actual
Cassels–Tate pairing. Let V⊂S have a certified independent basis of
size d. Suppose its complete restricted matrix A has rank c, so
R=rad(λ|V) has dimension d−c.

For any certified classes x₁,…,xₑ∈S, form the rectangular matrix

\[
C_R=(\lambda(r_i,x_j)),\quad a=\operatorname{rank}C_R,
\]

where rᵢ is a basis of R. The xⱼ need not be known independent.

**Theorem J9 (radical-partner amplification).**

\[
\boxed{\operatorname{rank}\lambda|_{V+\langle x_j\rangle}
       \ge c+2a.} \tag{J9}
\]

Furthermore the xⱼ span at least a dimensions modulo V, and

\[
\dim\bigl(V\cap\operatorname{rad}(\lambda|_S)\bigr)\le d-c-a.
\]

The lower bound c+2a requires no pairings between the xⱼ themselves.
In particular, one nonzero pairing with a restricted radical vector
certifies a new independent Selmer class modulo V and a rank increase
of at least two in the obstruction pairing.

**Proof.** Split V=H⊕R with λ|H nondegenerate, dim H=c.
An invertible a×a minor of C_R selects independent radical combinations
r₁,…,r_a and a partners x₁,…,x_a. A combination of these partners lying
in V would pair trivially with R, so their images modulo V are independent.
Subtract suitable elements of H from the partners to make them orthogonal
to H; this leaves their pairings with R unchanged. On their span with
r₁,…,r_a the pairing has block matrix

\[
\begin{pmatrix}0&B\\B^{\mathsf T}&D\end{pmatrix},
\qquad \det B\ne0.
\]

Its kernel is zero: the first block row forces the partner coordinates
zero, then the second forces the radical coordinates zero. This works
for every alternating D. Together with H it gives c+2a nondegenerate
dimensions. Finally a vector in V orthogonal to all of S must lie in R
and in the kernel of the pairing with the partners, of dimension d−c−a. ∎

The distinction from arbitrary rectangular pairing rank is essential.
Only pairings against the **restricted radical** receive this factor-two
amplification over the already certified c. A partner lying in V has zero
such column. A class built from a known rational point pairs to zero
with all S and is equally uninformative as an obstruction partner.

A deterministic useful query is therefore: reduce existing CT data to
its radical, construct a new genuine Selmer class, and test its column
against that radical first. Retain a nonzero column even before computing
new-new entries. If all columns vanish, the lower bound stays c; additional
obstructions can still live among the new classes themselves or in
uncomputed classes. This is not a theorem that nonzero columns can always
be found or that the remaining radical is rationally soluble.

The retained replay applies J9 to every principal prefix of all seven
fixed-cubic matrices, comparing the predicted lower bound against the
independently computed rank of the whole retained matrix. This is a
retrospective algebra check on certified arithmetic matrices, not a new
class construction or a measurement of prospective savings.

## J10. Stop an exclusion calculation when the minor is sufficient

Assume a certified unconditional upper bound U≥dim Sel₂(E/K), and set
t₂=dim E(K)[2]. If a valid CT certificate, including one from J9, proves
pairing rank at least c, then

\[
\boxed{\operatorname{rank}E(K)\le U-t_2-c.} \tag{J10}
\]

**Proof.** Rational Kummer classes have dimension rank E(K)+t₂ and lie
in the full CT radical, whose dimension is at most dim Sel₂−c.
Subtract t₂ and use U. The Kummer image and its vanishing in CT are
the standard facts in
[Fisher–Schaefer–Stoll, Section 2](https://www.dpmms.cam.ac.uk/~taf1000/papers/casselspairing.pdf).
No finiteness assumption on Sha is required. ∎

To exclude a target rank R*, it suffices to reach the smallest nonnegative
even integer c* with

\[
c^*\ge U-t_2-R^*+1.
\]

Equality U−t₂−c=R* does **not** exclude rank R*. Once c≥c*, additional
pairing entries are unnecessary for this target. A certified independent
point lower bound equal to U−t₂−c proves exact rank instead.

Without U the routine must return UNKNOWN. Neither dim V, a lower
Selmer span, nor a relative dimension involving an unbounded unknown
excess can be substituted for U. The six existing fixed-cubic controls
have bounds 1+ε,1+ε,3+ε,1+ε,2+ε,2+ε, with the same UNKNOWN
ε=dim Sel₂(E₀)−20. Replaying their matrix ranks does not turn these
expressions into numerical rank upper bounds.

**Search use.** First obtain a valid upper envelope or a mathematically
new feasible way to certify one. Then choose the smallest CT minor that
can exclude the desired threshold, using J9 to prioritize partners.
Do not restart a previously stalled class-group descent unchanged merely
because the stopping theorem is now explicit. The immediate arithmetic
objective remains one uncensored same-family negative control or a
certified small rational/Sha switch, with its transfer limits stated.

## Current evidence and operational priorities

The new compact certificate independently recomputes:

- both 17+4-point finite-signature ranks, the exact quotient relations,
  the rational two-torsion exclusions, and all six two-block overlap tables;
- seven fixed-cubic alternating ranks and all their principal-prefix
  radical-partner lower bounds, without using new-new entries for the bound;
- the latest pair-plus-triple relation counts at 3/8, including
  7+11−13=5 overlapping relation dimensions and six new constraints;
- UNKNOWN handling when a complete Selmer upper envelope is absent.

The u=−1 arithmetic input additionally has its full retained 153-entry
Fisher-formula certificate replayed by the existing Sage verifier.
The other six matrices are inherited arithmetic evidence with exact
matrix replays here; they are not claimed as fresh local-symbol calculations.

The current priority order is: (1) identify independently reachable quotient
directions beyond the deficient native dictionaries; (2) seek a new feasible
source of a certified Selmer upper bound; (3) use radical partners to reach
an exclusion threshold; (4) study a controlled rational/Sha switch;
(5) only then evaluate a candidate rule on separately frozen arithmetic
controls. The completed degree-six/eight panel is a negative regression.
The original J1–J6 mathematical statements remain in force.

## Replay

~~~sh
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_block_rank_theory.py
python3 elliptic-curves/rank-jump/verify_block_rank_theory.py check
sage -python elliptic-curves/cas/verify_fixed_cubic_cassels_tate.sage --check
sage -python elliptic-curves/rank-jump/verify_low_degree_triple_panel.py check
python3 scripts/render_status.py --check
~~~

The implementation is [block_rank_theory.py](block_rank_theory.py), and
the immutable application certificate is
[rank_jump_block_rank_theory_v1.json](../../artifacts/generated-results/elliptic-curves/rank_jump_block_rank_theory_v1.json).
New tests exhaust all alternating matrices through size five and all
eight completions of a three-partner missing block in selected examples;
these check the implementation and the claimed missing-entry boundary.
The general proofs are J7–J10 above. No formal proof-assistant verification,
external review, full new rank bound, or new point-search result is claimed.
