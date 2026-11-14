# Xns121: adopted methods and exact finite-matrix replay

The recent complete rational-point claim is by Santiago Arango-Piñeros, Filip
Najman and David Zureick-Brown,
[Rational points on the non-split Cartan modular curve of level121](https://hexagonmath.org/2611.00231v1),
November10,2026. Its authors explicitly say their proofs have not yet been fully
checked by humans. The [source repository](https://github.com/F-Najman/Xns-121-)
is pinned to root tree8426c0ae44265f06805020d72042542fee265fe3 in this import;
the seven finite data files used below are pinned by exact Git blob identity.

The manuscript and code have been read for their argument structure.
Neither the full claimed theorem nor its global Selmer computations is adopted
as an established dependency of our EC programme.

## What the argument actually does

The genus511 modular curve maps to an elliptic quotient of rank1. Rational
points upstairs satisfy additional Kummer/integrality conditions on that
quotient. A six-prime Mordell-Weil sieve restricts the quotient multiplier to
seven residue classes modulo4356. Each residue class leads to a genus5 curve
over a degree15 number field. An isogeny induced by1-zeta11 on its Jacobian is
then used for Selmer-group Chabauty.

The important ingredients for transfer are separate:

1. a low-genus quotient carries the global parameter information;
2. global admissible Kummer classes and local conditions are intersected by
   exact linear algebra, before expensive point or basis searches;
3. a complete global Selmer bound localizes injectively at a useful place;
4. the formal direction forced by a hypothetical new point lies outside that
   localized space; repeated isogeny division converts this into an exclusion.

The formal lifting lemma requires both localization injectivity and equality
of local and global isogeny-kernel torsion. Completeness of the global Selmer
bound is essential. A restricted span of known classes or a matrix of sampled
local points cannot supply that premise.

The manuscript exploits cyclotomic automorphisms and controlled cuspidal
reduction of curves Y^2=X^11+b^2. Our generic K3/third-cover curves have no such
automorphism supplied. Its special1-zeta11 formulas cannot simply be copied.
An analogous2-isogeny or Jacobian2-covering strategy requires its own proof.

## Independent finite arithmetic

[verify_xns121_localization.py](../cas/verify_xns121_localization.py) uses only
Python's standard library. It parses data-only Magma assignments, checks exact
Git blob identities and reimplements Gaussian elimination over F11. It imports
no Magma, Sage or upstream executable code.

For the four even residue classes(-2,0,2,4), it verifies at each of the two
relevant places:

- the supplied localized Selmer rows have dimensions5,5,4,4, respectively,
  retaining the supplied global row rank;
- the torsion line is in their span;
- the formal line is outside that span;
- the intersection of the Selmer span with formal-plus-torsion is exactly
  the torsion line.

For each odd class(-1,1,3), the one-dimensional supplied Selmer localization
and the formal direction are nonzero and independent at both relevant places.

All14 class/place checks pass. The
[replay certificate](../../artifacts/generated-results/elliptic-curves/xns121_localization_replay_20261113_v1.json)
records source hashes and the checked dimensions. This proves finite matrix
statements only. We have not independently constructed the complete global
Selmer spaces, identified their arithmetic local bases, justified p-adic
precision or replayed the formal lifting proof. Therefore this reproduction is
not a proof that Xns121 has exactly seven rational points.

## Applied to our new objects

The [1152 pilot](../../elkies-k3/RECENT_PAPERS_X1152_TRANSFER_2026-11-13.md)
implements a first exact local gate for the15 paired quadratic-cover bases:
primitive projective parameter classes modulo2^3,3^2,5^2,7^2,11^2 must make
both radicands squares. Denominators are cleared by squares, and only even
common valuations are removed. The rational constant squareclass is preserved.

None of these15 pairs is excluded. This is evidence to change the next action:
do not merely add good finite primes to this sieve. A smooth genus-one curve
over a finite field has a point; good-prime local existence will usually carry
no global obstruction. Global covering/Selmer-class compatibility is the
appropriate next gate.

The positive-rank genus-one parameter base is intended to have infinitely many
rational points. The finite-point target appears only after extra conditions,
for example a third independent quadratic cover. Three quadratic radicands with
six distinct geometric branch points give a degree8 multiquadratic cover of
genus5: Riemann-Hurwitz gives2g-2=8*(-2)+6*4=8. This cover maps to the three
paired genus-one quotients. Their complete global covering data could constrain
possible rational points without building a large model of the entire curve.

Before using a quotient Mordell-Weil sieve globally, certify the quotient
group/index information or an appropriate complete Selmer substitute. A sieve
on multiples of one known point only controls that cyclic subgroup. A rank
upper bound or a saturated basis must not be inferred from a point search.

The Hill paired parameter curve is explicitly available. A bounded PARI
ellrank(0) diagnostic on its minimal model returned the interval[1,3]. This is
a diagnostic, not a certified exact rank or a complete Mordell-Weil basis.
Its retained nontorsion point remains the established lower-bound evidence.
No full-rank-dependent sieve has been run on this curve.

The November14 [paired-Jacobian transfer](../../elkies-k3/X1152_PAIRED_JACOBIANS_AND_FIRST_DESCENT_2026-11-14.md)
supplies the previously missing global gate for the first generic pair[0,1]:
complete PARI2-descent gives rank2 and Sha[2]=0, while an independent replay
certifies solubility at every local place. This proves an infinite rational
paired base, with the complete global computation explicitly retained as a
trusted CAS dependency. It does not establish the Xns121 global hypotheses,
a saturated quotient basis or a third-cover finite-point result.

## Replay and remaining verification

    python3 research/elliptic-curves/cas/verify_xns121_localization.py

The finite matrices are under[data/external/xns121](../data/external/xns121/).
Our existing SubspaceDescent/LocalKummerBasis infrastructure can carry restricted
known-class calculations, but a future global exclusion must explicitly supply
a complete global bound and verified local/formal hypotheses. This import does
not change its completeness boundary.
