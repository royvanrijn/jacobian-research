# F7: the sharp root-multiplicity bound five

**GAOF7MULT1:** every root of Gao F7's degree-twelve eliminant has
multiplicity at most five. The bound is sharp at Y=0:
R=t^5*(t^7-2t^4+2t-1), with seven other simple roots.

Combined with [GAOF7CUBE1](../gao-f7-three-point-certificate/PROOF.md),
this bound left four lists. The subsequent
[GAOF7FIVE1 proof](../gao-f7-quintuple/PROOF.md) classifies all118 normalized
fivefold targets, each with seven-point fibres. Only(4,3,2) and(3,2,2,2),
global three-point existence and JC2 remain **OPEN**.

[Proof and exact scope](HIGH_ROOT_PROOF.md),
[2426-check receipt](high-root-receipt.json),
[separate verifier](verify_high_root.py),
[certificate, about739KB](high-root-certificate.json),
[coefficient input](high-root-bivariate-input.json),
[exact sharpness witness](fivefold-witness.json),
[tested dependencies](requirements-high-root.txt).

The replay takes about29seconds. It reconstructs all multiple roots,
checks three polynomial resultants at1157 exact integer points with two
determinant engines, verifies a full623x623 integer Sylvester determinant,
and replays both exceptional-case unit identities in two integer engines.
There is no independent end-to-end review, formal proof or established
novelty priority.

The earlier [four-variable seven-two derivation](derive.py),
[equations](seven-two-input.json) and [transport receipt](derivation-receipt.json)
are retained as a superseded route. Both its rational lift and the later
bivariate direct lift hit their20second caps; neither is treated as proof.
The general linear-family coefficient of w remains
384*((z-u)*(z+3r)+2*(v-h)); any future division still requires a zero case.
