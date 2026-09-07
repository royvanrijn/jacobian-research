# Collision defects are pairwise, but span every native character

The fixed collision support does give a smaller arithmetic description of
the simultaneous lift. It does **not** collapse the successful quartets
to one shared local defect. Every pair of native factors supports a
certified local valuation defect, in both successful systems and the
obstructed control. Thus the pair graph is complete in all three cases.

This is a **solubility** calculation on fixed, retrospectively selected
cover systems. No exceptional coordinates enter its local witnesses.
Neither the collision graph nor its dimension is a Mordell–Weil rank
predictor.

## A finite squareclass target with pair restrictions

Use the [previous carrier and fixed-prime criterion](COLLISION_PRIMES_CONTROL_THE_REMAINING_LIFT.md):

\[
H:\ y^2=\prod_{i=1}^4 f_i(t),\qquad
C:\ u_i^2=f_i(t),\qquad
S=\{p:p\mid\prod_{i<j}\operatorname{Res}(f_i,f_j)\}.
\]

For a nonbranch local product point, write its valuation defect as

\[
\epsilon_p(t)=\big(v_p(f_1(t)),\ldots,v_p(f_4(t))\big)\bmod2.
\]

This vector has even weight. At a primitive projective parameter modulo
p, its nonzero entries can occur only among the factors that vanish
there together. Consequently each common projective root gives an
upper bound: the even-weight masks supported on its vanishing factors.
Take the union over roots, together with zero. This is a union of local
possibilities, not generally a vector space.

The new certificate computes every common projective root over Fp for
every p in S, including infinity. No prime can be removed just because
its collision has no Fp-rational branch parameter: all collision primes
in these systems have at least one such parameter.

Nevertheless, the allowed masks reduce the global squareclass target
substantially. A rational product point in either successful system has
all four native values positive and even valuations outside S. Its four
global squareclasses are therefore determined by the masks at S.

| System | Support primes | Necessary choices at each prime | Upper bound on global squareclass tuples |
|---|---:|---|---:|
| +7, 08234-003 | 18 | Two at every prime | 2^18 = 262,144 |
| +8, 08234-009 | 23 | Three at 2 and 7; two at each other prime | 3^2 · 2^21 = 18,874,368 |

These are upper bounds, not counts of realized rational classes. No
twist population is enumerated. The three affine descent coordinates
are **squareclass-valued coordinates**, not three bits: before collision
restrictions, their prime support alone would allow 2^54 or 2^69 tuples.
The remaining global problem is which of the bounded tuples the curve
H(Q) actually realizes, particularly the all-zero tuple.

For the control, the corresponding valuation-mask count is
3^2 · 4^2 · 2^21 = 301,989,888. This omits its separate real sign
possibilities and is not a bound on its full signed squareclass image.
The control is an obstructed quartet, not an observed low-rank fibre.

## Exact local image at a simple pair collision

Call a prime ordinary here if p≥5, it occurs to exponent one in exactly
one pair resultant, there is one affine collision root r, and both
colliding factors have nonzero derivative at r. These conditions are
checked, not inferred from the label.

Suppose the pair is i,j. Write t=r+ps and set

\[
a_k=f_k(r)/p\pmod p,\qquad b_k=f'_k(r)\pmod p
\quad(k=i,j).
\]

The two linear functions a_i+b_i s and a_j+b_j s have distinct roots;
the determinant certifying this is retained and independently recomputed.
The other two factors are units. Put

\[
G(s)=(a_i+b_i s)(a_j+b_j s)
       \prod_{k\ne i,j}f_k(r)\in\mathbf F_p[s].
\]

If G(s) is a nonzero square, both colliding values have valuation exactly
one and the product has valuation two with square unit. Hensel's lemma
gives a Qp-point on H. It does not lift to C over Qp, since the two
individual odd valuations forbid square roots. Its mask is e_i+e_j.

Let L be the leading coefficient of G and let χ be the quadratic
character. Exactly

\[
\frac{p-2-\chi(L)}2
\]

residues s give this outcome. To see this, translate the distinct roots
to obtain a scalar multiple of u(u−1). The character sum of u(u−1) is
−1: completing the square reduces the point count to
v^2−w^2=1, which has p−1 solutions for odd p. There are two zero values
of G, so separating its positive and negative character counts gives
the formula. In particular the count is positive for p≥5.

The computation additionally retains explicit witnesses, using at most
64 residues for each one. Separate witnesses give zero valuation at all
four factors and a square product unit. Thus the **exact** local parity
image at each certified ordinary prime is

\[
\epsilon_p(H(\mathbf Q_p)_{\rm nonbranch})=\{0,e_i+e_j\}.
\]

Zero valuation does not itself imply a native local lift: nonsquare
units may remain. The equality above concerns valuation parity only.

## The paired result rejects a common one-character explanation

| System | Ordinary primes with both local witnesses | Realized pair masks | Span in the native even-parity space |
|---|---:|---|---:|
| +7 quartet | 13 | AB, AC, AD, BC, BD, CD | 3 |
| +8 quartet | 17 | AB, AC, AD, BC, BD, CD | 3 |
| Obstructed quartet | 18 | AB, AC, AD, BC, BD, CD | 3 |

For example, the +7 system realizes AC at 163, BD at 937 and AB at
4,706,591. These three masks are independent. The full output records
all six pairs and all prime witnesses. Every requested witness was
found inside the declared bound; the largest s used was nine.

The number three here refers to the span after forgetting which prime
carried the defect. The separate prime coordinates are distinct: a
defect at one prime cannot cancel a defect at another in Q*/Q*^2.
Using the retained rational product point only as an explicitly
retrospective reference at all other places, either local witness can
be chosen independently at each of the 13 or 17 primes. This gives
2^13 or 2^17 different valuation patterns of adelic product points.
It does **not** show those patterns arise from H(Q).

Equivalently, the Chinese remainder theorem can combine the finitely
many local residue conditions on t. That supplies no rational y solving
the product equation. The loss of that implication is precisely the
global rational-point problem, not a shortage of ways to satisfy these
individual local conditions.

At the remaining primes the certificate gives only necessary masks.
It does not assert that each upper-bound mask is locally realized there.
In particular the control's known obstruction at 23 remains separate
from its eighteen ordinary-prime witnesses.

## Consequence for the mechanism search

The equation-defined condition is sharper: the affine target is a
finite collection of pair-supported squareclass tuples. But all three
systems have the same complete graph of realizable pair defects.
That graph cannot explain why the retained two quartets are rationally
soluble while the control is obstructed. A smaller global image or a
global construction could still synchronize the classes; local span
three does not exclude either possibility.

The next high-value target is therefore the **global image of H(Q) in
this affine squareclass set**, via the degree-eight isogeny descent or
an explicit rational intersection component. A Jacobian rank bound or
a local graph alone does not establish that the embedded curve meets
the zero class. No incidence or visibility claim follows from this
local audit.

## Reproduction and failure record

The [protocol](COLLISION_DEFECT_PROTOCOL.json) fixes three existing
systems, one worker and sixty seconds per system. The first implementation
called Sage's root routine on a constant polynomial over a large finite
field; all three workers failed with an NTL degree error. The original
script and `rank_jump_collision_defect_v1.json` retain that failed run.
The corrected [producer](collision_defect_v2.py) returns no roots for
constant nonzero polynomials. All three corrected workers completed.

The immutable successful result is
`artifacts/generated-results/elliptic-curves/rank_jump_collision_defect_v2.json`.
Its [independent verifier](verify_collision_defect.py) uses Python integer
arithmetic, its own polynomial Euclidean algorithm, the quadratic
discriminant root count, and exact valuation and Euler-criterion checks.
It neither uses Sage nor searches for witnesses.

```sh
python3 elliptic-curves/rank-jump/verify_collision_defect.py check
```

The verification certificate is
`artifacts/generated-results/elliptic-curves/rank_jump_collision_defect_verification_v1.json`.
All work is confined to rank-jump analysis files and immutable outputs.
