# Exceptional PSL2(11): tame residue covers and a global parity obstruction

Work over an algebraically closed characteristic-zero field k.
Suppose a polynomial Keller self-map F:A^d_X -> A^d_Y has exceptional
degree11 geometric action PSL2(11)/A5. Write K=k(Y), L=k(X), and M/K
for its generic closure. Both conjugacy classes of A5 stabilizers are
covered below.

This is an absolute constraint on any such Keller map, in every
dimension. It is not a construction-specific test. It excludes
certain arrangements of residue covers, but does not exclude the
exceptional action itself.

[Exact group helper](psl11_residues.py) ·
[checker](verify_psl11_residues.py) · [model](psl11_residue_model.json) ·
[receipt](psl11_residue_receipt.json) ·
[retained primary inputs](sources/psl11-residues-manifest.json).

## Statements

**KMFI81 (prime-to-degree Brauer visibility).** Let a characteristic-zero
polynomial Keller map have generic degree n, and let l be a prime
not dividing n. Every nonzero alpha in Br(K)[l] has a nonzero residue
after restriction to L at some prime divisor of the ordinary affine
source. That source prime maps to an affine target prime, has
ramification index1, and its residue is the restriction of the
target residue. In particular a nonzero class cannot lose all of
its ramification merely because the corresponding finite-normalization
divisors are absent from the affine source.
Over a nonclosed characteristic-zero coefficient field, the same assertion
holds for classes not coming from Br(k)[l]; constant classes must be excluded.

**KMFI82 (complete tame decomposition and surviving-sheet table).**
For the exceptional11 action, nontrivial affine inertia can have
orders only2,3,5. In the following table I is inertia, D is the finite
decomposition group, and the listed degrees are the D/I-orbit sizes
on the I-fixed sheets. The actual affine source retains a nonempty
subset of these orbits; each retained orbit is one source prime
of that residue degree. C_r means cyclic order r, and D12 means
dihedral order12.

| I | D | Available unramified residue degrees | Ramified finite-normalization primes (e,f) |
|---|---|---|---|
| C2 | C2 | 1,1,1 | (2,1) four times |
| C2 | V4 | 1,2 | (2,1) twice; (2,2) |
| C2 | C6 | 3 | (2,1); (2,3) |
| C2 | D12 | 3 | (2,1); (2,3) |
| C3 | C3 | 1,1 | (3,1) three times |
| C3 | C6 | 2 | (3,1); (3,2) |
| C5 | C5 | 1 | (5,1) twice |

Every ramified prime in the right column is absent from the affine
source. Unretained unramified orbits may also be absent.
In particular the rows with a single available orbit force a unique
source prime of degree3,2 or1 as listed. In these single-orbit rows,
a reduced equation of the target prime has irreducible pullback under F.
A decomposition group need
not fix one sheet: only its inertia does. Order6 is forbidden as
affine inertia but permitted in the decomposition groups in this table.

**KMFI83 (exact Schur Brauer residues and the central lifting gate).**
Let beta in Br(K)[2] be the inflation of the nonsplit central extension
\[
           1\to C2\to SL2(11)\to PSL2(11)\to1.                \tag{1}
\]
Its affine target residue is nonzero exactly at the involution
divisors with D=V4 or D=D12. It is the sign character of the residue
action on the three I-fixed sheets, equivalently the discriminant character of
their full rank3 etale residue algebra, including omitted tails.
In all other rows it is zero.
On the affine source this residue remains nonzero precisely on
a retained degree1 orbit in the V4 row, or the retained degree3
orbit in the D12 row. The degree2 orbit in the V4 row kills it.

Consequently beta=0, equivalently the generic PSL2(11) representation
has a proper SL2(11) Galois lift, if and only if every affine
involution decomposition group is cyclic (C2 or C6).
This is a field lift, not an affine Keller realization of SL2(11).
The lifted A5 stabilizer is the perfect binary group of order120,
so the cyclic etale source-quotient construction from KMFI65 does not apply.

**KMFI84 (global obstruction to retaining only the quadratic tails).**
If any affine involution divisor has noncyclic decomposition, then
some such divisor must have an actual affine source prime of odd
residue degree. Equivalently, at least one of the following occurs:

- a D12 involution divisor, with its mandatory degree3 source prime;
- a V4 involution divisor retaining its degree1 orbit.

Thus a nonempty collection of V4 involution divisors cannot all retain
only their degree2 orbits when no D12 divisor occurs.
The individual quadratic-orbit choice is locally allowed. The
simultaneous choice everywhere is impossible for an ordinary-affine
Keller source of odd degree11.

Cyclic decomposition groups throughout are still allowed by these
results. The exceptional degree11 action and the full spectrum
therefore remain open.

## 1. The general visibility lemma

Restriction and corestriction on Brauer groups satisfy
\[
       \operatorname{cores}_{L/K}\operatorname{res}_{L/K}
                        =[L:K]=n.
\]
Since l does not divide n, restriction on l-torsion is injective.
Hence alpha_L is nonzero. If it had zero residues at every prime
of A^d_X, Brauer purity would extend it to Br(A^d_X)[l]=0,
a contradiction.

A source prime maps to an affine target prime: F is polynomial,
etale and quasi-finite, so codimension is preserved.
The ramification index is1. Functoriality of Brauer residues
therefore identifies its residue with restriction of the target
one. This proves KMFI81 without assuming F finite.

Over a nonclosed field, source purity instead puts an unramified
alpha_L in Br(A^d_X)[l]=Br(k)[l]. Subtract that same constant class
from alpha in Br(K)[l]. Its restriction is zero, so injectivity
makes the difference zero. Thus only constant Brauer classes can
lose all affine-source residues. Without this qualification, an
identity Keller map over a field with nonzero Br(k) is a counterexample.

The cited inputs are restriction/corestriction (Section3.8),
residue functoriality (Proposition1.4.6), purity (Corollary3.7.3)
and affine vanishing (Theorem5.1.1) in the retained 2019
Colliot-Thelene--Skorobogatov draft. The residue edge used below is
Theorem1.4.1 and Definition1.4.2. Tame inertia and its conjugation
are Stacks Lemmas15.114.5 and15.114.8. Exact locators and source
hashes are in the retained-input manifest.

## 2. Tame centrality and all decomposition groups

The cycle types of the exceptional action are
\[
\begin{array}{c|c}
\text{element order}&\text{cycle type}\\\hline
2&2^4\,1^3\\
3&3^3\,1^2\\
5&5^2\,1\\
6&6\,3\,2\\
11&11.
\end{array}
\]
KMFI1 excludes the last two as affine inertia because they fix
no sheet. It also ensures that at least one fixed-sheet orbit
is present above each affine target divisor.

Inertia is cyclic in characteristic zero. Its conjugation by the
decomposition group is the action on the corresponding roots of
unity in the residue field. All these roots lie in k, so the action
is trivial: D is contained in C_G(I), not just its normalizer.
The centralizers of generators of orders2,3,5 have orders12,6,5,
respectively. They are D12,C6,C5.

The overgroups of C2 inside D12 correspond to the subgroups of
D12/C2=S3, giving C2,V4,C6,D12. The overgroups of C3 inside C6
are C3,C6; the order5 case has only C5.
Their actions on the fixed sheets give the left half of the table.

For an orbit O of D on all eleven sheets, the common I-orbit
size is its ramification index e, and f=|O|/e.
This produces the right half of the table. The finite normalization
of the target has one prime for each such D-orbit.
Every actual source prime has e=1, and is the same generic prime
as one of these unramified orbits by Zariski's main theorem.
The affine source is an open of the finite normalization; it may
omit additional unramified primes. The table does not assume they
are all retained.

The checker constructs both A5 coset actions directly in the
660-element matrix group. A5 subgroups contain a (2,3,5) generating
pair, and all involutions of G are conjugate. Enumerating pairs with
one fixed involution and then conjugating exhausts the A5 subgroups:
there are two disjoint conjugacy classes of eleven subgroups.
All centralizers, their overgroups and both actions are checked.

## 3. The local residue as a commutator character

Only inertia2 requires a new calculation. Odd inertia orders3,5
have unique prime-to-2 complements in their central preimages,
as in KMFI65, giving zero residue. Trivial inertia also gives
zero residue.

At an involution i, its centralizer is D12. Choose a matrix lift A
of i in SL2(11), and a lift B of d in C_G(i). Their commutator
is in {I,-I}. Define
\[
           \chi_i(d)=[A,B]\in C2.
\]
It is independent of lift choices, is a homomorphism and is
trivial on I=<i>. Thus it is a character of D/I.

Here is why it is exactly the Brauer residue, rather than merely
a finite-group diagnostic. Over the henselian local field at the
divisor, the tame inertia is procyclic of cohomological dimension1.
The split residue sequence identifies the residue of a central
2-cocycle with its Hochschild--Serre edge in
H1(k(D),Hom(I_tame,C2)).
As all roots of unity are fixed, inertia and a chosen residue
section commute. For a normalized cocycle c, the edge evaluated
on an inertia generator tau and residue element sigma is
c(tau,sigma)-c(sigma,tau), precisely the commutator of their lifts.
After the finite quotient it is chi_i(d).

The lift A has order4, but this causes no additional residue term.
Procyclic inertia itself can lift to order4; its H2 vanishes.
A purely cyclic C4->C2 obstruction inflates locally to zero when
fourth roots of unity are present. The mixed commutator is the
residue term that remains.

In the full involution centralizer, chi has kernel the cyclic
torus C6. It coincides with the sign of the permutation on the
three i-fixed sheets. Therefore it is nontrivial on V4/C2 and
D12/C2, and zero on C2/C2 and C6/C2.
The residue Galois group is exactly D/I, so a nontrivial such
character gives a nonzero field residue.

Equivalently, form the etale rank3 residue algebra of all unramified
sheets in the finite normalization over the target prime. Its
discriminant character is the permutation sign on its embeddings,
so it is precisely chi. It is essential to include unramified tails
that the actual source may omit; otherwise the lifting obstruction
could be incorrectly erased.

For a retained source orbit and a point s in it, the source residue
field corresponds to the stabilizer D_s/I. Pullback restricts chi
to this subgroup. In the V4 row the fixed orbit of size1 retains
the character, while the size2 orbit kills it. In the D12 row
S3 acts on three letters; a point stabilizer C2 still has nontrivial
sign. This proves the exact target and source residue assertions.

## 4. The global lifting criterion

If all involution decomposition groups are cyclic, every affine
target residue of beta is zero. Purity and Br(A^d_Y)[2]=0 give
beta=0. Conversely any V4 or D12 divisor gives a nonzero residue,
so beta is nonzero. This proves the equivalence with the table.

The vanishing of beta is the continuous central embedding-problem
condition for(1). A weak lift has image surjecting onto PSL2(11).
If it missed the prime central kernel, it would be a group section,
contradicting nonsplitting of(1). Hence every such weak lift is proper.

This argument does not remove order4 lifts of involution inertia.
Nor does the preimage of A5 have a nontrivial cyclic quotient:
both lifted stabilizers have order120 and are perfect.
Thus a proper SL2(11) lift here is not the forbidden etale source
cover from KMFI65, and it is not a Keller realization.

## 5. The global source obstruction

By the preceding section beta is nonzero whenever at least one
V4 or D12 involution divisor occurs. Since [L:K]=11 is odd,
KMFI81 at l=2 forces a nonzero residue of beta_L at some affine
source prime.

The complete local table shows exactly where this can happen:
the degree1 orbit at a V4 divisor or the degree3 orbit at a D12
divisor. If every V4 divisor retained only its quadratic orbit
and there were no D12 divisor, beta_L would be unramified at
every affine source prime and hence zero. Corestriction would
then give 11*beta=0, contrary to beta being nonzero2-torsion.
This proves KMFI84, including nonfinite maps and all dimensions.

For a local control, the V4 central preimage is Q8. Its residue
character can be represented by the quaternion symbol (s,pi)
over k(s)((pi)). The quadratic residue extension adjoining sqrt(s)
kills the residue. An odd-degree residue extension cannot kill a
nonzero quadratic character. Thus the hidden quadratic orbit is
locally legitimate; the contradiction is genuinely global.
The odd degree is essential to the injectivity step.

These results neither assume that every decomposition group fixes
a sheet nor that a source retains all unramified tails.
They do not infer nonoccurrence from a formal allowed action.
The general geometric proofs, cited inputs and exact complete
finite group checks have separate roles; no formal verification
or independent whole-proof review is claimed.
