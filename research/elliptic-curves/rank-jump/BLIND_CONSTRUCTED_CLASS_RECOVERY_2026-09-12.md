# Blind recovery of both constructed strict classes

**Both frozen covers were solved without the known point oracle.** On MW16-05
at `t=3/17`, a separate worker received only their equations and rational norm
invariants. It recovered rational points on both. Independent exact arithmetic
checks the entire transport, and a fresh finite quotient certificate proves
that these two points add two independent directions to the generic sixteen.

Ordinary V3, independently started from that same generic subgroup, recovered
six extra directions and also supplies exact lifts of both fixed classes.
This closes the point-oracle gap on this control. It does not identify a
parameter condition forcing the block, establish a new rank record, or prove
that the complete class-construction route is faster than V3.

## Inputs and the successful route

The earlier [class-block report](CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md)
already established two additional strict classes and ordinary ideal-class
directions, then verified rational lifts using known points. Its
[adaptive principal-relation constructor](TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md)
is the required prior implementation for future collectors. Its successful
class construction and the later sparse two-field relation pilots have
different targets; neither computes a full class-group upper bound.

This experiment reused compaction columns **6 and 7** unchanged. The
[frozen protocol](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/protocol.json)
binds their [equations](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/covers.json)
and the [generic subgroup](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/generic.json).
No relation collection or previous search was restarted. The cover worker had
a fresh context and a single mathematical-input allowlist: the cover file.
It received no exceptional coordinates, field square roots, point combinations,
or V3 data. This was an instruction and dataflow boundary on a shared filesystem,
not an operating-system access sandbox. The fibre and classes were user-selected
positive controls; their selection is retrospective.

The established [Cremona–Fisher–Stoll algorithms](https://www.dpmms.cam.ac.uk/~taf1000/papers/minred-234.pdf)
are implemented by Magma's documented
[genus-one minimisation and reduction](https://magma.maths.usyd.edu.au/magma/handbook/text/1592).
The worker used degree-four `Minimise` and `Reduce`, parametrised the surviving
ternary conic, then used degree-two `Minimise` and `Reduce`. PARI searched the
resulting binary quartics. All rational transformations were retained.

| Frozen column | Original representative bits | Reduced degree-four bits | Quartic parameter | Primitive reduced degree-four point |
|---|---:|---:|---|---|
| 6 | 2464 | 39 | `1/2` | `[614050,-4,-1,0]` |
| 7 | 1766 | 44 | `3/10` | `[14345125,-4,1,-6]` |

Both quartic searches succeeded at abscissa-height bound **1000**. Before that,
degree-four `PointsQI` searches at bounds1000 and100000 returned no points.
The displayed degree-four representatives have heights614050 and14345125,
outside those bounds. These are witness heights, not globally minimal heights.
The misses did not justify treating the covers as insoluble.

Degree-two minimisation with `CrossTerms=false` reported positive level at2;
this permitted feature of that model convention was not used as a local
obstruction. The degree-four positive-level sets were empty. Independent replay
certifies equivalence and the rational witnesses, without asserting a separate
proof of global minimality of the returned models.

## Exact verification

The [witness bundle](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/cover-worker-witness.json)
contains the original primitive `[u,v,w,s]` coordinates and complete4-by4
coordinate and2-by2 equation matrices. Two independent Fraction-only
implementations verify:

- every coefficient in the degree-four transformation;
- the conic parametrisation, discriminant quartic and binary change of variables;
- each returned square and its transport to the original cover;
- `beta*xi^2=4*x-theta`, the norm map and the original elliptic equation.

The resulting original Weierstrass abscissas are

```
column 6: 393529041691378317987472481 / 23565934144
column 7: 6645456506656578661752762977 / 514456169536
```

Their exact ordinates are in the
[independent cover replay](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/independent-cover-replay.json).
The [separate quartic replay](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/independent-quartic-replay.json)
binds the intermediate search hits. A fresh
[finite independence certificate](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/blind-lift-independence.json)
proves rank18 for the generic sixteen plus these two points. This confirms
independence directly from the recovered points; no historical point group is
needed. Five tamper regressions reject changed maps, points, classes and roots.

## Comparison with V3

The [comparison record](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/results.json)
separates setup, maps, search, replay and failed attempts. The maintained MW16
upper-shell adapter used the existing V3 selector, exact CVP, bounded maps,
immediate rebuild after a gain, and a100-chart cap at height125000. The generic
points were reconstructed from the sixteen Kummer representatives. Their exact
atlas transport, including sign changes in the height Gram, was checked first.
The baseline's artifact-read guard rejected external mathematical data.

| Measured component | Result | Recorded seconds |
|---|---|---:|
| Blind cover arm, charged arithmetic calls including network and failed preparations | Both fixed covers recovered | 18.401 |
| Its two final quartic searches | Both succeed at height1000 | 0.0247 |
| V3 preparation plus time to its second new direction | Generic16 to18 in2 charts | 38.011 |
| V3 preparation plus time to its sixth new direction | Generic16 to22 in6 charts | 103.962 |
| Full V3 preparation and100-chart exposure, including the initial failed preparation | All100 charts complete; final lower bound22 | 239.149 |

The two V3 timing milestones include the0.902-second failed preparation. The
full search phase took214.119 seconds, with60.129 seconds in point-search
backends. Its original and corrected replays together took19.677 seconds;
the successful replay checks every retained square, map, adaptive gain source
and final finite rank proof. It does not independently regenerate the full
V3 selector/CVP schedule or rerun point boxes.

After V3 sealed, a separate evaluator matched the fixed classes to combinations
of its **newly recovered** points. Exact field square roots and independent
rational group-law replay certify
[lifts of both covers](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/v3-fixed-cover-lifts.json).
Both combinations require its rank22 prefix. This is an exact Kummer-class
statement inside the certified subgroup, not merely character agreement.
Neither this evaluation nor its output entered either search arm.

These measurements do **not** establish end-to-end superiority. The cover arm
starts after the earlier class constructor has finished; that upstream cost is
excluded. Reduction ran on remote Magma2.29-10, whereas V3 used local
Sage10.9/PARI2.17.3. Algorithm implementation and documentation work are excluded
from the arithmetic timings: the worker's retained artifact interval spans
about20 minutes. A larger population and a common accounting of construction,
preparation and search would be needed for a performance claim.

## Retained failures and reproduction

The first verbose Magma export exceeded its normal output limit; a compact
export retained the same final transformation. An initial conic matrix
orientation and a subsequent variable-scope error failed before any quartic
search. The corrected identities passed before search. V3's first preparation
hit its external-cache guard and was replaced by explicit in-memory quotients.
Its first replay passed the arithmetic checks but compared tuple fields to JSON
lists incorrectly. The corrected replayer changes that serialization comparison;
the successful search and its frozen adapter are preserved unchanged.

The [evidence archive](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/evidence.zip)
and [manifest](../../artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/evidence-manifest.json)
retain inputs, programs, responses, misses, failed versions, maps, searches and
certificates. Restore only into an empty replay checkout, rooted at `research/`.
Do not overwrite a live tree or reconstruct any unrelated missing artifacts.

From the repository root, use fresh output filenames:

```sh
python3 research/elliptic-curves/rank-jump/verify_constructed_class_blind.py --witness research/artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/cover-worker-witness.json --output /tmp/blind-cover-replay.json
python3 research/elliptic-curves/rank-jump/replay_constructed_class_quartic.py --witness research/artifacts/local/elliptic-curves/constructed-class-blind-cover-v1/cover-worker-witness.json --output /tmp/blind-quartic-replay.json
sage -python research/elliptic-curves/rank-jump/certify_constructed_class_blind.py --output /tmp/blind-independence.json
python3 research/elliptic-curves/rank-jump/evaluate_constructed_class_v3.py verify --source research/artifacts/generated-results/elliptic-curves/constructed_class_blind_recovery_v1/v3-fixed-cover-lifts.json --output /tmp/blind-v3-class-replay.json
python3 -m unittest discover -s research/elliptic-curves/rank-jump -p test_constructed_class_blind.py
```

The reusable preparation lesson is concrete: freeze explicit classes, retain
the source-cover map, minimise and reduce, and exploit an available conic to
search a degree-two model. The remaining mathematical gap is the parameter
condition creating the additional arithmetic block; the remaining algorithmic
question is its full construction-and-recovery cost against V3.
