# Two bounded arithmetic transfer tests

These tests evaluate proposed external inputs on concrete current models.
They prove no improved abc exponent, semigroup intensity law, or two-copy
energy bound.

## A determinant-ratio cell in a saturated abc fibre

[rStation Paper 125/v1, Lemma 2](https://rstation.ai/problems/riemann_hypothesis/papers/125/v1)
counts primitive rows in the specific half-open strip
`R_A={(a,b): A≤a<min(2A,Q), Q≤b<Q+a, gcd(a,b)=1}`.
For ordered distinct rows `r=(a,b), s=(c,e)`, it sets
`D=|ae-bc|`, chooses `zc+we=1`, and takes `h_D=za+wb (mod D)`.
For `d|D` and one fixed unit residue `h mod d`, the number of
matching row pairs is at most `2(D+8A²)/d`. Its proof uses primitivity
in the determinant lattice and the half-open second-coordinate strip.
The paper's conditional RH application is not used here.

Use the [primitive two-divisor fibre](ABC_DIFFUSE_FACTOR_GRAPH.md#primitive-two-divisor-fibres-attaining-the-unsaved-count):
`x=pu`, `y=qv`, `U=x³(x+2y)`, `V=y³(2x+y)`, and
`T=(x+y)³(x−y)`, with `U=V+T`. The two full divisors are
`p³|U`, `q³|V`; each distinct actual triple has weight one.
The exact dictionary for the tested subfibre is

| rStation datum | abc datum |
|---|---|
| `A=N`, `Q=3N`, `r=(a,b)` | `r=(y,x)=(qv,pu)`, `N≤y<2N`, `3N≤x<3N+y`; the last condition selects the actual half-open strip |
| `D=|ae−bc|` | for two triples, `D=pq|vu'−uv'|` |
| `d|D` | `d=pq`; this is a row-level divisor, not the full moduli `p³,q³` of the abc terms |
| `h_D mod d` | `h_D≡v/v' (mod p)` and `h_D≡u/u' (mod q)`, when rows are primitive |
| `W(r,s)` | no depth weight occurs in this abc fibre; each triple is counted once |

The last ratio is **not fixed** by the `p,q` record. The deterministic
[finite checker](verify_abc_determinant_ratio_fibre.py) uses `N=10^10`,
`p=10000019`, `q=10000169`, `x≡1`, `y≡12 (mod 72)`,
`gcd(u,v)=1`, and the strip cut above. Python 3.12 with SymPy 1.14
reported 166 rows, 27,390 ordered nonzero-determinant pairs, 8,776
distinct determinant values and 27,390 distinct `(D,h mod d)` cells.
There were 130,302 ordered collisions at the same `D` and **zero** at
the same `(D,h)` in this finite box. The largest actual cell was one;
the largest displayed Lemma 2 upper bound was 24,819,411. These are
finite diagnostics, not uniform estimates.

The loss in the proposed transfer is exact: for the actual unit weights,
`sum_h M(D,h)=N_D`, where `M(D,h)` is the cell multiplicity. To use the
lemma's `1/d` one must either fix `h` through additional abc data or
average with a compensating character/denominator factor. The fixed
`p,q` record supplies neither. Summing all eligible residues can cost
up to `φ(d)`, cancelling the proposed gain. The next viable test is an
existing abc joint-determinant configuration with a *genuine* ratio
lock; its source and weight cost must be included before claiming any
saving.

## Pouly's stopped expansion on conditional `k=2`

[Pouly, Hexagon 2611.00447v1, Lemma 4.1 and Proposition 4.2](sources/pouly-2611.00447v1.pdf)
give a polynomial identity and a signed-measure remainder for any finite
ordered family of dependent indicators and a hereditary retained set
`D`. The ordered boundary consists of first excluded prefixes,
`I∉D` with `I\{max I}∈D`, not all inclusion-minimal excluded sets.
The bound is `||Law(Z)−T_D||_1≤sum_{I∈∂D}2^|I| E∏_{i∈I}X_i`.
The paper's further Poisson rates require arithmetic large-order tail
estimates not imported here.

At a fixed target `n`, expose the selected primes through `Y`. For a
representation event `i=(a,b)`, `a+b=n`, let `S_i` be the union of
prime supports of `a,b` above `Y`. If a factor at most `Y` is excluded,
the event is zero. Otherwise, at critical `k=2` selection probabilities
`p^{-1/2}`,

`E[∏_{i∈I}X_i | F_Y] = ∏_{p∈∪_{i∈I}S_i} p^{-1/2}`.

This union, rather than the product of separate event probabilities,
retains every shared unexposed prime. The
[finite checker](verify_semigroup_k2_stopped_expansion.py) takes `n=100`,
`Y=5`, primes `2,3,5` exposed and selected, and the eight displayed
representations in its output. The overlap graph has one cluster at
prime 7 and one at prime 11; the retained hereditary family allows
independent events and at most one interacting cluster. There are 192
retained subsets and 48 ordered boundary prefixes. Exhaustive
enumeration of the seven unexposed Bernoulli choices gives actual
`L¹` error `0.39317163`; the boundary bound is `2.35195599`.
The direct overlapping-pair mass is `0.55697228`. On this finite test
the boundary estimate does not supply a stronger usable upper bound.

The same calculation separates the two-copy energy. Conditional on
`F_Y`, independent completions have
`E[Z^(1)Z^(2)|F_Y]=(E[Z|F_Y])²`; here that value is `2.30010759`,
while the one-copy second moment is `4.16127308`. Within each
representation, `S_i` is already a union of the two summands' prime
supports. A one-copy stopped expansion cannot replace the sharp
[two-copy energy upper bound](SHIFTED_LOCAL_CORRELATIONS_AND_SEMIGROUP_ENERGY.md#sharp-second-moment-route-to-the-random-prime-intensity-law).
For an asymptotic transfer, a hereditary retained class must make its
full boundary sum `o(1)` (or beat the existing overlap target) uniformly
over the relevant targets and cutoffs, while preserving the exposed
local factors and the two independent completions. This remains open.
