# A known Q80 gain realizes exactly one common halving class

The first retained cheapest rational bisection, `alternate-orbit-11ae6`,
gives **exact generic rank18** on a rational quadratic base. It cannot
supply a second direction at any height. This is a positive control for the
[branch specialization mechanism](Q80_BRANCH_TRACE_SPECIALIZATION_GATE_2026-09-14.md),
not a solution of the [two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).

Write its literal cover as w²=q(t), with coefficients in ascending order:

```
98872578977556376101742883789238871426743106890776868117310060729
80148528835876681718319913092201091573928018290138078125694109004
16124196532885861882702787740901506922195833708883364381645924036
```

The [checker](scripts/verify_q80_one_cover_halving_control.py) freezes the
first record of the existing cheapest1024 bank. It does not search that
bank or reconstruct any bisections. It checks that q is irreducible, and
checks an exact rational point inherited from the first retained V4 pair.
Thus this smooth conic has infinitely many rational points. The existing
lifted section Q has nonzero coefficient of w in its abscissa; its equation
is rechecked by exact polynomial arithmetic in Q[t,w]/(w²-q).

Let alpha be the unique closed quadratic branch value, and let M be the
full saturated17-section group of the direct11952 Q80 model. At each simple
root of q modulo a good prime p, reduction gives a necessary test for

```
K_alpha = ker(M/2M -> E_alpha(Q(alpha))/2E_alpha(Q(alpha))).
```

The checker uses ascending primes5..997 and stops at matrix rank16. It
reduces the actual17 sections, enumerates each entire finite elliptic group,
constructs the subgroup of doubles, and labels the quotient by cosets.
It independently computes the rational-root Kummer characters and checks
agreement of row spaces at every residue specialization. The combined
matrix has rank16. Its kernel in the supplied section basis is

```
{0,65538} = span(e_1+e_16),   with zero-based indices.
```

This matches the retained `section_basis_w` trace parity (-e_1+e_16).
The different `direct_alternate_w` lattice coordinates must not be used
as coordinates in this matrix. Finite reduction alone proves only that
K_alpha is contained in this line; the known gain proves equality below.

At p=23, q has the simple root t=0, its leading coefficient is a unit,
and the reduced parent cubic has no root. The fibre is smooth there.
A rational2-torsion point over Q(alpha) would reduce to a nonzero2-torsion
point over F23 at this degree-one place, which is impossible. This also
certifies characteristic-zero smoothness of the conjugate branch fibres.
Consequently the actual branch image of every anti-invariant is zero.

The retained Q80 unramified Kummer theorem and
[exact trace-norm rank formula](Q80_GOOD_BRANCH_CODE_AND_TRACE_NORMS_2026-09-14.md)
therefore give

```
rank gain = dim(Tr(E(Q(C)))/2M) <= dim K_alpha <= 1.
```

On this base change all fibres are irreducible and chi=4. The height
formula excludes nonzero torsion. Since the verified Q differs from its
conjugate, Q-sigma(Q) is nonzero and nontorsion, giving gain at least1.
Hence the gain, realized trace dimension and dim K_alpha are all exactly1.
This conclusion has no height bound and does not require good branch
reduction at131. It does not bound ranks of individual rational fibres.

The [certificate](../artifacts/generated-results/elkies-k3-q80-one-cover-halving-control-v1/result.json)
contains all residue matrices, both quotient computations, input hashes,
branch coefficients and the rational base point. The bounded computation
took under0.1 seconds with20 CPU seconds and1GiB allowed. Replay with:

```
python3 research/elkies-k3/scripts/verify_q80_one_cover_halving_control.py --record /tmp/q80-one-cover-replay.json
```

Parent rank/saturation, the global unramified theorem and the geometric
rank argument are inherited or written mathematics, not formal verification.
The dual finite calculations share one harness; no separate independent
replay of the whole theorem is claimed. Only this fixed cover is closed.
For a new two-gain candidate with branch-field2-torsion absent, a common
halving kernel of dimension at least2 is necessary, and its classes must
still be realized by global trace torsors on the same cover.
