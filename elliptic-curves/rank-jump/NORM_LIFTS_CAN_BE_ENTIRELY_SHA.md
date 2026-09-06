# A complete unramified norm-lifting block can consist entirely of Sha

The two strict classes with explicit unramified norm lifts are **both
non-rational on the original small curve**. They generate its entire
`Sha[2]`. On its `-1` twist, the same two classes have explicit independent
rational representatives. This closes the rational/Sha distinction in the
small control and gives a counterexample to promoting an unramified
norm-lifting block to a Mordell–Weil block.

| Fixed small quotient | Exact rank | Full 2-Selmer dimension | Strict rational dimension | `dim Sha[2]` |
|---|---:|---:|---:|---:|
| `E0: y^2=x^3-11x^2-14x-1` | **1** | 3 | **0** | **2** |
| `E+: y^2=x^3+11x^2-14x+1 = E0^(-1)` | **3** | 3 | **2** | **0** |

These are exact ranks of two fixed arithmetic controls, not new high-rank
curves or rank upper bounds on production specializations.

## Fixed block and bounded descent

Keep `K=Q(theta)`, `f(theta)=theta^3-11theta^2-14theta-1=0`, and
`S={2,163,infinity}`. The preceding
[nonscalar calculation](NONSCALAR_CUP_BLOCK_AND_SELF_GLUING.md) fixes
\[
U=\langle\beta_0,\beta_1\rangle,\qquad
\beta_0=\theta^2-10\theta+1,\quad
\beta_1=\theta^2-13\theta+12.
\]
The classes were obtained from cubic arithmetic, before the rational-point
test. They form the complete strict S-class character space. Both belong
to both elliptic 2-Selmer groups: their local classes at `S` are zero,
and elsewhere their classes are unramified and the curves have good
reduction. Thus the experiment holds strict incidence fixed.

The [protocol](SMALL_QUOTIENT_SOLUBILITY_PROTOCOL.json) permits one
deterministic PARI `ellrank` computation with effort zero and a 30-second
cap per curve. It supplies only the immediately verified points
`(-1,1)` on `E0` and `(0,1)` on `E+`. The
[descent output](../../artifacts/generated-results/elliptic-curves/rank_jump_small_quotient_descents_v1.json)
is `[1,1,2,[(-1,1)]]` and
`[3,3,0,[(-1,5),(0,1),(2,5)]]`, respectively. The third field reports
`dim Sha[2]/2 Sha[4]`; it is not treated as a point count.

The proof below re-establishes the rank bounds from retained arithmetic,
exact point witnesses and the scalar cup identity. It does not use these
new descent upper bounds as its proof.

## Independent full-Selmer bound

Let `V_S` be the norm-square cubic classes with even ideal valuations
outside `S`. The retained certified class group is `C2 x C2`.
There is one prime above 2, since `f` is irreducible modulo 2. There is
one prime above 163, since
\[
f(Z+58)=Z^3+163Z^2+8802Z+157295
\]
is Eisenstein at 163. These primes are principal: `(2)` is the first;
the second has class killed by three, hence trivial in the class group
of exponent two.

Multiplication by 2 and 163 removes arbitrary valuation parities at
these two primes. Every resulting everywhere-even class is a product
of a unit and the two half-ideal classes supplied by `beta0,beta1`.
Those half ideals form a class-group basis: their retained Artin matrix
is invertible. Unit squareclasses have dimension three; their norm map
has rank one, detected by `-1`. The two half-ideal representatives have
norm `625`, a square. The rational factors 2 and 163 have independent
non-square norm valuations, so the norm-square condition removes both.
Consequently
\[
\dim V_S=(3-1)+2=4.
\]

The norm-one units `-1-theta` and `1+12theta-theta^2` have real sign
vectors `011` and `101`. They span the even-sign plane, while both
strict classes are positive. For either elliptic curve the real Kummer
image is a one-dimensional line in that plane. Hence the real condition
alone gives `dim Sel2 <= 3`. Conversely `U` and the supplied point with
nontrivial real Kummer image give three independent Selmer classes on
each curve. Thus **both full 2-Selmer dimensions are exactly three**.

In fact the supplied points identify the full spaces in the common cubic
norm kernel:
\[
\mathrm{Sel}_2(E_0)=U\oplus\langle[-1-\theta]\rangle,\qquad
\mathrm{Sel}_2(E_+)=U\oplus\langle[\theta]\rangle.
\]
Their intersection is exactly `U`, since the two nonzero real sign lines
are distinct. This identifies the inherited block and the changing local
boundary separately.

The two lines in the common real sign plane differ between `E0` and `E+`.
Their full Selmer groups need not be identical; their common strict block
is the same. No uncomputed bad-place condition is replaced by a heuristic:
the lower bound already reaches the real-condition upper bound.

## Rational witnesses and the obstruction on the original curve

Use Kummer labels `x-theta` on `E0` and `x+theta` on `E+`.
Exact good-prime characters at 37, whose cubic roots are `25,27,33`,
together with one real character, give point fingerprints

| Curve / points | Four-bit masks | Rank |
|---|---|---:|
| `E0: (-1,1)` | `8` | 1 |
| `E+: (-1,5), (0,1), (2,5)` | `13,8,6` | 3 |

These prove independence; no numerical heights or analytic rank are used.
The twist's three point classes therefore fill its full 2-Selmer group.
In particular `Sha(E+)[2]=0` and its CT pairing is zero.

There are small explicit representatives of the strict basis:
\[
Q_0=(6,23)=(-1,5)+(0,1),\qquad Q_1=(2,5).
\]
Their Kummer identities are
\[
\beta_0=(6+\theta)
 \left(\frac{4\theta^2-45\theta-16}{23}\right)^2,
\qquad
\beta_1=(2+\theta)(\beta_1/5)^2.
\]
The [point/block certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_small_quotient_block_v1.json)
retains the exact combinations and square roots.

On `U`, the independently verified scalar cup matrix is
`H=[[0,1],[1,0]]`. The strict twist-comparison identity gives
\[
\operatorname{CT}_{E_0}|_U+
\operatorname{CT}_{E_+}|_U=H.
\]
Thus `CT_E0|U=H`. Every nonzero element of `U` is detected by another
element, so none is rational on `E0`. They inject into `Sha(E0)[2]`.
The full Selmer dimension three now bounds the rank by one; `(-1,1)`
proves the matching lower bound. This gives exact rank one and
`Sha(E0)[2] ~= (Z/2)^2`, with `U` mapping isomorphically onto it.

This argument uses the general strict scalar comparison proved in
[the scalar cup note](INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md), without
assuming that the original classes are rational. Rationality on the
twist was proved above.

## The paired two-cover equations

For `z=u+v*theta+w*theta^2`, write
\[
\beta_i z^2=Q_{i,0}+Q_{i,1}\theta+Q_{i,2}\theta^2.
\]
The two covers for each class differ by a sign:
\[
\begin{array}{ll}
E_0:&Q_{i,2}=0,\quad Q_{i,1}+h^2=0,\\
E_+:&Q_{i,2}=0,\quad Q_{i,1}-h^2=0.
\end{array}
\]
In monomial order `(u^2,uv,uw,v^2,vw,w^2)`, the explicit forms are

| Class | `Q1` coefficients | `Q2` coefficients | Rational projective point on the twist cover |
|---|---|---|---|
| `beta0` | `(-10,30,30,15,730,4240)` | `(1,2,52,26,602,3676)` | `(2,11/5,-1/5,1)` |
| `beta1` | `(-13,52,-54,-27,108,242)` | `(1,-4,8,4,34,241)` | `(2/5,1/5,0,1)` |

The map is `x=Q0/h^2`, `y=25*N(z)/h^3`. It sends the displayed points
to `Q0=(6,23)` and `Q1=(2,5)` respectively. All forms and point checks are
in [the cover certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_small_quotient_covers_v1.json).
Both original covers are everywhere locally soluble but have no rational
point, by their nonzero Sha classes. There is no omitted rational point
at `h=0`: then `beta_i*z^2` would be rational; taking norms would make
that rational number a square, contradicting the nontrivial class `beta_i`.

## What this settles, and what it does not

1. **Solubility:** this is an exact simultaneous switch of two globally
   rational elliptic classes, with the same strict incidence space and
   an independently computed obstruction. The soluble side is established
   by explicit points, not by a zero CT radical alone.
2. **Disproved sufficiency:** the entire explicit unramified norm-lifting
   block for `gamma-=-(1+theta)` consists of Sha classes on `E0`.
   Passing that lifting test is therefore strictly weaker than rational
   solubility, even when the auxiliary construction is an elliptic
   self-gluing. The local Jacobian lifting conditions are not determined
   by this computation and are not silently declared satisfied.
3. **Weak Galois explanation:** the smaller normal closure and the larger
   unramified norm image occur in the sign case whose second elliptic
   quotient has rank one. The other quotient has rank three. Neither
   smaller descent fields nor more norm lifts is a monotone rank predictor.
4. **Missing production bridge:** the nine nonscalar production CT bits
   and their non-strict local corrections remain uncomputed independently.
   The production cubics are `S3`, whereas this control is cyclic. Its
   self-identification and ordinary-unit norm-sign criterion cannot be
   transplanted to those fibres.
5. **Next mechanism target:** isolate a production-available condition
   that annihilates the relevant full CT obstruction and also supplies
   rational representatives of several independent quotient classes.
   A zero restricted map or an auxiliary norm solution is now excluded
   as a sufficient replacement for that second requirement.

The [independent verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_small_quotient_block_verification_v1.json)
uses rational arithmetic, prime and real characters, the retained certified
class group, and the established cup identity. It uses no new descent,
class-group computation, norm solver, or analytic estimate. No production
search file, policy, population, or status entry was modified.

```sh
sage -python elliptic-curves/rank-jump/small_quotient_block.py --check
python3 elliptic-curves/rank-jump/verify_small_quotient_block.py --check
python3 elliptic-curves/rank-jump/small_quotient_covers.py --check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_small_quotient_block.py
```
