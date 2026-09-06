# Governing cochains change obstruction while Selmer incidence stays fixed

For an elliptic curve over Q with S3 two-torsion action, suitable positive
prime twists can preserve **every local Kummer image and the full 2-Selmer
group**, yet prescribe an independently specified block's Cassels--Tate
pairing. This is an existence result and a second-descent condition, not
an explicit twist or a rational-point construction.

We derive the prime-twist statement below, verify its finite S3 hypotheses
for three retained production controls, and independently check the
smallest nontrivial governing group. No twist prime, point, parameter,
class group, or new CT entry is computed.

## External input and applicability

Morgan's [arXiv:2309.02374v2](https://arxiv.org/pdf/2309.02374v2),
Proposition 3.3, gives a Frobenius formula for CT variation under locally
trivial twists ramified at fixed-point-free two-torsion primes. Proposition
4.13 and Corollary 4.18 describe the governing extension and permit arbitrary
pairing values. Theorem 5.8 preserves the full Selmer group while prescribing
an admissible form. Section 2 resets the setting to arbitrary principally
polarized abelian varieties; Section 5.2 requires its stated module
hypotheses, not the introduction's dimension restriction for Kummer varieties.
These pairing results do not require Sha finiteness.

This removes the full-rational-2-torsion restriction of the older route
considered in our CT notes. It does not apply a twist theorem to the
non-twist fixed-cubic pencil. The S3 specialization and prime refinement
below are our deductions from these inputs.

Write V=E[2]=F2^2 and G=GL2(F2)=S3. The action is simple since G acts
transitively on the three nonzero vectors. Its commuting endomorphisms
are zero and identity; a three-cycle has no nonzero fixed vector.
Also H1(G,V)=0: averaging kills positive-degree C3 cohomology, V^C3=0,
and inflation--restriction applies to the normal C3 in S3. The elliptic
polarization has rational symmetric representative O, so c=0. The
theorem's inflated Selmer intersection is therefore zero, and its
admissible forms are alternating.

The [finite gate](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_cochain_gate_v1.json)
enumerates six invertible matrices, sixteen endomorphisms and all 4^6
one-cochains. Exactly four are cocycles and they are exactly the four
coboundaries. It rechecks G=S3 on the frozen A1/MW16-05 307/206 and
published R17 -2300/843 and -1561/3133 models.

## The extra arithmetic datum

Fix independent Selmer classes a_1,...,a_k and representing V-valued
cocycles. Local isotropy and Brauer reciprocity make their global pairwise
cup products zero. Choose cochains gamma_ij satisfying

\[
d\gamma_{ij}=a_i\cup a_j.
\]

At a Galois element sigma acting without fixed vectors, set

\[
P_{i,\sigma}=(\sigma-1)^{-1}a_i(\sigma),\qquad
\psi_{ij}(\sigma)=e_2(P_{i,\sigma},a_j(\sigma))+
\gamma_{ij}(\sigma)\in\mathbf F_2.
\tag{1}
\]

Pairings are written additively as bits. The cochain term is additional
information beyond the ordinary class cocycles. Let Sigma contain 2,
infinity, bad reduction and governing ramification. For a twist trivial
on Sigma, ramified elsewhere only at fixed-point-free primes, its CT
matrix on the block changes by

\[
\Delta_{ij}=\sum_{p\text{ newly ramified}}\psi_{ij}(\mathrm{Frob}_p).
\tag{2}
\]

Thus vanishing on the block requires Delta to equal its original CT
matrix. Full annihilation also requires zero cross-rows against the
remaining Selmer classes. This is a **solubility-obstruction condition**;
neither version guarantees rational points.

For two classes, Gamma=(V x V) semidirect S3 has order 96. Its cup cocycle
is omega((a,b,g),(a',b',h))=e_2(a,gb'). Adjoin a central bit z with
(x,z)(y,z')=(xy,z+z'+omega(x,y)). The gate verifies the cocycle equation
on all 96^3=884,736 triples and (1) under 12,288 conjugations. Of 64
fixed-point-free lifts, 32 have each governing value. Changing only z
reverses that value while leaving both class values and g unchanged.

The [independent tuple-group replay](../../artifacts/generated-results/elliptic-curves/rank_jump_governing_prime_compatibility_v1.json)
reconstructs multiplication, finds commutator subgroup order 96 in the
order-192 group, and checks that all fixed-point-free lifts lie in it.
This is finite algebra, not an explicit production number field.

## Prime-twist refinement over Q

Put m=binomial(k,2). The governing group is a central extension

\[
1\longrightarrow\mathbf F_2^m\longrightarrow\mathcal E
\longrightarrow V^k\rtimes S_3\longrightarrow1,
\qquad |\mathcal E|=6\cdot2^{2k+m}.
\]

Its abelianization is C2. Indeed, commuting translations supported on
distinct coordinates i,j produces exactly their central ij bit
e_2(v_i,v_j); suitable vectors produce each central basis vector. After
killing these commutators, commutation with a three-cycle spans V^k
because g-1 is invertible. The remaining abelianization is S3^ab=C2.
For k=1 the same argument omits the central part.

Let F be the actual governing field and choose

\[
M=8\prod_{\ell\in\Sigma,\ \ell\text{ odd finite}}\ell.
\]

The maximal abelian subfield of F is its cubic discriminant field, whose
conductor divides M. Hence F intersects Q(zeta_M) in exactly that
quadratic field. Every element projecting to a three-cycle is trivial
on the intersection, independently of its central bits. Chebotarev in
the compositum gives infinitely many positive primes p outside Sigma with

\[
p\equiv1\pmod M,\qquad
\mathrm{Frob}_p|_{E[2]}\text{ a three-cycle},\qquad
(\psi_{ij}(\mathrm{Frob}_p))_{i<j}=D
\tag{3}
\]

for any prescribed vector D. The character of Q(sqrt(p)) is trivial at
Sigma and ramifies only at p. At p the entire local H1(Q_p,E[2]) is zero;
at other good unramified places both local point images are unramified.
Thus **all local Kummer images agree**, the full Selmer group agrees,
and (2) gives Delta=D.

No prime satisfying (3) is computed here. This condition concerns a
**twist parameter**, not t in the original MW17/MW16 family.

## A second-descent baseline

Fix the block, cochains and M. Conditional on the first two conditions
in (3), all 2^m governing vectors have equal Chebotarev density: for each
fixed three-cycle and translation part, the central bits translate
bijectively onto the vector. The congruence constraint removes none of
these central choices, by the preceding abelianization argument.

The density of any prescribed restricted CT matrix is therefore 2^-m.
For the three existing strict rational blocks, the original matrix is
zero. The resulting theorem-derived baseline is:

| Original fibre | Fixed block dimension k | Pairing slots | Conditional prime density of zero restricted CT |
|---|---:|---:|---:|
| A1/MW16-05 307/206 | 10 | 45 | 2^-45 |
| R17 -2300/843 | 8 | 28 | 2^-28 |
| R17 -1561/3133 | 6 | 15 | 2^-15 |

These are not measured search frequencies or a rank distribution. The
denominators reflect the chosen block sizes and do not discriminate the
original gains. The classes are labelled retrospective rational controls.

For a full Selmer space of dimension s, a fixed k-dimensional subspace
lies in the full CT radical precisely when binomial(k,2)+k(s-k)
independent alternating-matrix entries vanish. With governing data for
that full space, the corresponding conditional density is
2^(-binomial(k,2)-k(s-k)). This separates internal isotropy from full
annihilation; it remains a second-descent statement.

## Consequences and missing implications

An alternating form on an s-dimensional space can have rank
2 floor(s/2). Apply (3) to a full Selmer basis with this target form.
There are therefore positive prime twists with the identical full
Selmer group and

\[
\operatorname{rank}E^{(p)}(\mathbf Q)\le s\bmod2\le1.
\]

This is unconditional: E^(p)(Q)[2]=0 and rational Kummer classes lie in
the full CT radical. Odd s gives at most one, not exact rank one.
An unknown full Selmer dimension prevents an explicit construction here
but does not prevent this existential conclusion.

Likewise, prescribe a nonsingular form on any of the even-dimensional
strict blocks in the table and extend it to the full space. Some positive
prime twists then make **every nonzero class in that same block Sha**,
while every local Kummer condition and the entire Selmer group remain
unchanged. This is stronger than our old -1 contrasts, where boundary
dimensions also contracted.

Conversely, zero CT only gives second-descent lifting. Higher-divisible
Sha can survive, even if Sha is finite. Further descent or rational
witnesses remain necessary. Twisting also need not preserve the original
rational generic subgroup: these existential low-rank curves are not
ordinary fibres of the original R17/A1 family. None has been computed.

The resulting mechanism to investigate is **degeneration of a governing
obstruction matrix while incidence stays fixed**. The precise arithmetic
level is now identified: cochains and their Frobenius values. We have not
proved that such degeneration forces a large rational jump, or identified
a special low-degree carrier for the unexplained production directions.

## Next falsifiable computation

Construct one actual gamma_ij for an independently fixed pair of
production Selmer classes. Represent its degree-192 governing extension
by a tower or equivalent cocycle data and certify its ramification. This
two-class task does not require a full Selmer basis. Known rational
classes are permitted only as retrospective controls.

Success must identify the extra central value arithmetically and make
(2) evaluable without a point search. Reconstructing only the order-96
individual-class field, or fitting gamma to a known CT entry, fails this
endpoint. A later twist test must meet the local-square and inertness
hypotheses; the existing -1 twists do not. No new field construction or
twist run is claimed in this turn.

## Replay and provenance

```bash
python3 elliptic-curves/rank-jump/governing_cochain_gate.py check
python3 elliptic-curves/rank-jump/verify_governing_prime_compatibility.py check
```

Both finite computations completed in under a second per command. The
versioned external PDF has SHA256
`ec4d757f37167d4c42a43580ddd3c208e6386eb6e16d15cc45fa02d7128835e1`.
It is retained only in the ignored local research directory. The finite
certificates bind existing production class and Galois evidence. The
theorem application, general-k group argument and Chebotarev deduction
are written proofs, not outputs of a descent solver. Concurrent theorem
navigation and all active search work are left untouched.
