# Additive actions on the canonical etale completion

2026-11-16. Claims HC4MYLM0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The preceding goal turn made
PROGRESS by retiring the fixed gamma-map HC4 conversion queue (LK0/1).

We now return to the finite canonical source algebra. The results hold
in EVERY field degree. They constrain the maximal etale completion's
additive actions and the GENERIC RATIONAL DECK group, not the monodromy
group. They do not prove trace regularity or eliminate all unramified
boundary. In particular, known higher-power Danielewski controls prevent
an exact-volume-only affine-plane recognition argument.

## 1. Retained canonical completion and the prior classification argument

Let K:X=A2->Y=A2 be an actual complex polynomial Keller map, with
R=C[x,y], B=C[P,Q], L=Frac(R), F=Frac(B). Write S for the full finite
normalization over B, Lambda for its trace-dual and

    A=B[Lambda], U=Spec(A).

KU0/1 identifies U with the FULL etale locus of Spec(S)->Y. It is
smooth affine and contains X openly. U minus X is the disjoint union
of exactly b0 omitted UNRAMIFIED primes, all smooth A1s. Retained
LA0/1 proves

    A*=C*, Pic(U)=Z^b0, pi1(U)=1,
    H2(U,Q)=Q(-1)^b0.

The etale target map supplies the nowhere-zero EXACT volume

    omega=dP wedge dQ=d(P*dQ).

The finite generic rational deck group

    G=Aut_F(L)=N_M(H)/H

acts REGULARLY and FREELY on U. Here M is the full normal-closure
monodromy group and H its sheet stabilizer. G is not M. LA0 gives
|G| dividing both d=[L:F] and 1+b0; LA1 gives b0>=2|G|-1 when G is
noncyclic. None of these retained proofs/checkers is rerun.

We reuse the PREVIOUSLY IMPORTED Daigle classification from
UAMLFACTOR1: Theorem3.2 of *Affine surfaces with trivial Makar-Limanov
invariant*, arXiv:0705.0334v1, published J.Algebra319(2008),3100-3111.
For a smooth affine complex surface with ML=C, trivial canonical bundle
forces a presentation uv=p(z). Smoothness makes p square-free.
Its complete statement and proof passage were re-read in the retained
primary HTML. This is an existing proved import, not a new classification
proved by our checker. Over C the cited proof invokes Bandman--Makar-Limanov.

UADANVOL1 already excludes uv=p(z) of degree>=2 as an etale source over
A2: its unique-up-to-scalar nonvanishing volume has a NONZERO period,
whereas the target volume is exact. Degree1 gives A2. The combination
of classification and periods is already the proof of UAMLFACTOR1;
we credit and REUSE it, rather than presenting it as a new theorem.

## 2. LM0: trivial ML recognizes exactly the absence of added unramified primes

For the ACTUAL canonical completion above,

    ML(A)=C  iff  b0=0  iff  U=X=A2.                 (1)

Indeed ML=C, the trivial canonical bundle and the retained
classification/period argument give U=A2. Then Pic(U)=0, so the
retained boundary Picard basis forces b0=0. Conversely b0=0 gives U=X
and its two translations have kernel intersection C.

This does NOT prove that U has trivial ML. It also does NOT force the
Keller inverse: b0=0 permits omitted RAMIFIED primes outside U, and
A=R can hold in a still-unresolved nonautomorphic plane source.

### Proper completions have at most one additive ruling

If b0>0 and A has any nonzero locally nilpotent derivation D, then
EVERY nonzero locally nilpotent derivation has the SAME kernel C[f].
Thus ML(A)=C[f], of transcendence degree1. If no such D exists, then
ML(A)=A and U is additively rigid. These are the only two possibilities
for a proper actual U.

Here are the standard surface-kernel details, matching the argument
already retained in UAAFFDESC1. The quotient of a nontrivial additive
action on a normal affine surface is an affine normal curve, with
generic fibre A1. The quotient theorem is already explicitly imported
there. Since Frac(A)=C(x,y), that curve is rational; since its units
embed in A*=C*, its smooth affine model is A1. Hence ker(D)=C[f] and
Frac(A)=C(f)(t).

If two kernels C[f1],C[f2] shared a nonconstant h, their fraction fields
would both be algebraic over C(h). Each C[fi] is algebraically closed
in A and each C(fi) is relatively algebraically closed in Frac(A),
by the generic A1 fibre (or the characteristic-zero derivation
constant-field argument). Both fields therefore equal the relative
algebraic closure of C(h) in Frac(A). The kernels themselves are
A intersect C(fi), so they are equal. Distinct kernels consequently
have intersection C, which would imply ML(A)=C, contradicting (1).

This concerns additive actions, equivalently affine-base A1 rulings.
It does NOT identify distinct rulings with COMPLETE base P1 or assert
that every rational affine surface has an additive action.

### A usable test on the finite canonical algebra

If two polynomial source LNDs with DISTINCT invariant fields both
preserve the embedded algebra A=B[Lambda], then b0=0. Their restrictions
remain locally nilpotent on A; the fraction fields are unchanged, so
their kernels there are distinct and (1) applies.

This is a condition on the actual finite canonical generators, not on
a guessed abstract surface. There is no automatic preservation of A by
source translations. Conversely an action on U need NOT preserve its
chosen open plane X; restricting a complete field to an arbitrary open
subset need not remain complete. We do not infer a source LND in that
direction without preserving X. Also A-preservation is weaker than
preserving the target field F; they must not be confused.

## 3. LM1: noncyclic generic decks force additive rigidity

For EVERY actual complex plane Keller source, in EVERY degree,

    U admits a nonzero LND  ==>  G=Aut_F(L) is CYCLIC. (2)

Equivalently,

    G NONCYCLIC  ==>  ML(A)=A,
                       no nonzero additive action on U.           (3)

If b0=0, LA0 already gives |G| dividing1, so G is trivial. Suppose
b0>0 and U has an LND. By LM0 there is a UNIQUE additive ruling and
ML(A)=C[f]. Every automorphism of U preserves this intrinsic algebra.
Thus G acts on A1_f by affine transformations. We prove that action
is FAITHFUL; merely knowing that an action on the base exists would
not be enough.

Suppose 1!=sigma in G acts trivially on f. On the generic fibre it
acts as an automorphism of C(f)[t], hence

    sigma(t)=a*t+b, a in C(f)*, b in C(f).

The restriction is nontrivial, since the localization has fraction
field L. If a=1, finite order in characteristic zero gives b=0 and
sigma=1, contradiction. If a!=1, it has the fixed point t=b/(1-a)
on the WHOLE generic affine line. This is a point of the generic
fibre of the fixed locus of sigma in U, so that fixed locus is
nonempty; after clearing finitely many denominators it has a closed
complex point. That contradicts the FREE deck action from LA0.
Therefore G embeds in Aut(A1).

Every finite group of affine transformations of C has a common fixed
point: average the finite orbit of any point. Translating that point
to0 embeds the group in C*, whose finite subgroups are cyclic. This
proves (2), including the faithfulness step, without a finite-group
classification or an integral boundary-permutation assumption.

Consequently every NONCYCLIC subgroup of the generic deck group
requires the rigid branch. LA1's previous cost b0>=2|G|-1 is retained;
the new restriction adds that NO source polynomial LND can preserve
A in that branch. It does NOT exclude the noncyclic group itself.

In particular noncyclic MONODROMY does not imply (3). A transitive
monodromy action can have self-normalizing H and G=N_M(H)/H=1.
The natural S3 action is a familiar control of that distinction;
it is not asserted to be a possible plane degree3 Keller map.

## 4. Sharp cyclic and exact-volume controls: do not overextend the result

Reuse LA1's KNOWN surfaces, for EVERY k>=2:

    W_k={x^2*y=z^k-1},
    (x,y,z)->(epsilon*x,epsilon^-2*y,epsilon*z), epsilon^k=1.

They are smooth affine, simply connected, contain an actual whole
A2 chart with k-1 omitted A1s, and have scalar units, Pic=Z^(k-1)
and an EXACT nonvanishing volume. Their free CYCLIC action and all
these properties are retained evidence, not recomputed here. No
etale target pair W_k->A2 is supplied.

There is an explicit additive action compatible with the new restriction:

    D(x)=0, D(z)=x^2, D(y)=k*z^(k-1).

It preserves x^2*y-z^k+1. It is LOCALLY NILPOTENT: D^2(z)=0 and

    D^j(y)=k!/(k-j)! * z^(k-j)*x^(2j-2), 1<=j<=k,
    D^(k+1)(y)=0.

It is nonzero everywhere: for x!=0 use D(z); at x=0 the equation
gives z^k=1 and D(y)!=0. Its kernel is C[x], since after localizing
x it is x^2 partial_z and a rational function of x regular on every
nonempty x-fibre has no finite pole. The generic cyclic action on
the ruling is x->epsilon*x and is faithful, exactly as in (2).

Thus proper exact-volume plane enlargements CAN have one additive
action and a free nontrivial cyclic symmetry. The false inference
"exact volume plus a whole plane chart forces A2" is specifically
excluded by the retained family. The absence of a compatible etale
target remains essential. No all-Ga-rigidity theorem for proper U,
automatic trivial ML, automatic cyclic generic deck, or JC2 proof
is asserted.

## 5. Assurance and remaining target

All classification, quotient and period inputs above were ALREADY
retained in UAMLFACTOR1/UAAFFDESC1/UADANVOL1. The Daigle source was
re-read to check that no Gizatullin, hypersurface or free-Ga hypothesis
had been silently added or dropped. No NEW external theorem is adopted.
The broad period idea first fails on the retained higher-power family;
only the ML-trivial branch is excluded, as the classification permits.

The checker verifies the NEW additive derivation/nilpotence identities,
the affine fixed-point algebra, generic-deck/monodromy distinction and
retained proof/source byte pins. It does not replay the earlier source
chart, exact primitive, free cyclic action, classification or completion
checkers. Universal assertions are written proofs, not finite searches.
Independent replay, formal verification and external review remain false.

Full HC4/JC2 remains OPEN. Even b0=0 only identifies the polynomial
source with its maximal etale completion; it does not make the finite
normalization etale or force Tr(Lambda*Lambda) subset B. Noncyclic decks
with a rigid completion, trivial decks, and the ramified-boundary gap
remain open. This is an all-degree structural restriction, not a
replacement success condition for the full objective.
