# Rational j-incidence does not add sections to the first 32 curves

The pinned cohort of the first 32 certified new curves was checked against the six compact
R17 presentations, five compact MW16 presentations, and published R17 model.
The **384 curve/presentation pairs** have exact outcomes: **336 have no
rational j-preimage**, and 48 have completely certified rational preimages.
There are no unresolved residual factors or preimages at infinity.

Every curve occurs in its original family. The sixteen curves belonging to
the published R17/`08234` pair each occur in both presentations. Those are
the same generic family and the same integral generic point subgroup:

\[
s=-26t-50,\qquad
A_{08234}(s)=26^4 A_{\rm pub}(t),\qquad
B_{08234}(s)=26^6 B_{\rm pub}(t).
\]

After dividing the compact section coordinates by \(26^2,26^3\), all seventeen
sections are exact integral combinations of the published sections. The
retained matrix has determinant **−1**. Every group identity was checked over
\(\mathbf Q(t)\), then replayed without numerical heights. Thus these duplicate
presentations add no directions to the existing specialized subgroups.

This closes one specific possible source of extra points on the first 32
curves. It excludes neither other families nor rational points outside the
retained generic subgroups, and it is not an upper bound on any curve's rank.

## Exact projective exclusions

For each family the rational function

\[
j(t)=\frac{6912A(t)^3}{4A(t)^3+27B(t)^2}
\]

is reduced to a coprime integral numerator/denominator pair, homogenized to
degree 24. A usable reduction prime requires no common zero over the algebraic
closure of the residue field, including infinity. The checker verifies
polynomial coprimality modulo that prime and that both leading coefficients
do not vanish together. The resulting map is a morphism of unchanged degree
on the projective line. Therefore any rational preimage must reduce to a
preimage in the finite projective line.

For each of the 336 excluded pairs, a retained good prime at most 251 gives a
complete projective image that omits the target j-value. This excludes **all
rational parameters**, rather than just a height box. Bad map reductions
are not used for exclusion. Target j-values with a denominator divisible by
the prime are checked against the image of infinity in the target line.

Surviving equations are factored over Q. Every linear factor gives an exact
rational preimage; each residual factor has its own good prime at most 997
where its homogenization has no projective root. Exact multiplication checks
the entire factorization, and the degree test handles parameter infinity.
The Sage-free replay recomputes all finite images, polynomial identities,
rational roots, residual exclusions and all 384 pair outcomes.

## Calibration corrections retained with the proof

Sixteen generic section words needed a numerical proposal before exact
verification. The initial adapter mistakenly passed projective triples to
PARI's height matrix routine; the corrected adapter supplied affine pairs.
At the initially chosen parameter t=1, the displayed seventeen sections
only gave finite mod-2 image rank 15 through 251, and the numerical Gram inverse
did not yield useful relation proposals. This is not an exact rank 15 claim.

The successful proposal uses the existing new curve at t=33/119, whose exact
certificate verifies the independence of its first seventeen generic points.
Its height matrix is explicitly computed at 320-bit precision. All proposed
words are then accepted only by exact function-field identities and a
unimodular matrix; the numerical calculation is unnecessary for replay.

An additional API audit found that Cypari2's `ellheightmatrix` defaults to
64-bit precision even after setting PARI's displayed decimal precision. The
new optional `prospective_half_lattice_v2.sage` explicitly requests 384 bits.
On each of the three certified rank 25 point sets, the old and new computations
gave identical Gram entries after the existing multiplication by 10^6 and
integer rounding. The largest observed entry difference was below
1.85×10^−17. These are finite numerical comparisons, not rigorous height
intervals. Exact point/rank proofs do not depend on either precision setting.

## Evidence and replay

The [portable manifest](../../artifacts/generated-results/elliptic-curves/compact_cross_family_incidence_evidence_v1.json)
and accompanying ZIP include the input models, exact certificates, replay
sources, all successful logs, and both failed proposal adapters.

```sh
python3 elliptic-curves/cas/replay_compact_cross_family_incidence.py \
  --input artifacts/generated-results/elliptic-curves/compact_cross_family_j_incidence_v1.json \
  --output artifacts/local/elliptic-curves/cross-family-replay-new.json
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/audit_compact_published_r17_transport_v3.sage --check \
  artifacts/generated-results/elliptic-curves/compact_published_r17_generic_transport_v1.json
```

The j-incidence discovery and replay had 300/180-second allowances and 1.5-GiB
memory limits. Each generic transport attempt and replay had 180 seconds;
the separate three-curve numerical precision comparison had 120 seconds.
