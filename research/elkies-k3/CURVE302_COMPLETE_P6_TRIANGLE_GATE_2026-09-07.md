# Complete P6 triangle gate for MW14–15 parents of 302

Authority: `EC-K3-CURVE302-COMPLETE-P6-TRIANGLE-GATE`.

On the pinned determinant-948 K3, consider every triangle
`D=O+P6+Q` with `height(Q)=6` and `<P6,Q>=3` in the original
[11952 presentation](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
**None gives a generic-rank-at-least-14 Jacobian parent of 302.**

There are 880 configurations modulo `Q -> P6-Q`. Of these, 786 have
arithmetic generic rank at most 13 if Jacobian. The other 94 comprise
88 exact MW14 pencils and the [six completed MW15 pencils](CURVE302_SIX_MW15_TRIANGLES_2026-09-07.md);
every one excludes 302 at all rational parameters. The 88 MW14 pencils
are pairwise inequivalent. This completes one entire fixed-anchor
triangle class, rather than just the former two-term-word dictionary.

The scope is essential: other anchors, nontriangle degree-three divisors,
other fibrations and other surfaces remain open. The 786 lower-rank
configurations are **not** excluded as possible parents of 302. No
generic-rank-at-least-14 parent of 302 has been recovered.

## Complete lattice coverage

The full arithmetic NS has rank 19 and matrix
`N=diag([[-2,1],[1,0]],-G)`, where G is the certified rank-17 height matrix.
An old section word w has class `S(w)=(1,w*G*w/2,w)`.

Exact rational LDL enumeration gives all 53,290 signed height-six words
and all 2,626 signed height-four words. The enumeration uses exact
rational costs and conservative integer bounds at every branch; it is
independent of PARI. It visits 412,168 nodes for height six, below the
declared five-million-node limit. Exactly 1,760 height-six words satisfy
`<P6,Q>=3`. The involution `Q -> P6-Q` pairs these into 880 configurations.
On the surface, the automorphism `R -> P6-R` transports the corresponding
triangles and pencils, so one representative covers both.

For every height-four old section R,

`D.R = 6 - <P6+Q,R>`.

For 786 representatives, two such sections R1,R2 are vertical. The five
classes `O,P6,Q,R1,R2` are independent: each record supplies a nonzero
five-column determinant. Their span contains D, so their images modulo
D provide four independent fibre-component directions. If D defines a
Jacobian fibration, Shioda–Tate bounds its arithmetic MW rank by
`19-2-4=13`. This argument does not require a complete root census for
those 786 cases.

The remaining 94 representatives match the frozen 88 MW14 and six MW15
protocols exactly, allowing `Q -> P6-Q`. No word or coefficient-support
restriction is imposed on the complete height-six shell.

The [gate certificate](../artifacts/generated-results/elkies-k3-curve302-anchor6-triangle-gate-v1.json)
records all 880 representatives, their rank witnesses or protocol
references, enumeration counts and a digest of the complete height-six
shell. The checker regenerates that shell from G.

## The 88 MW14 equations and inverse tests

The [word roster](../artifacts/generated-results/elkies-k3-curve302-anchor6-mw14-roster-v1.json)
was frozen before testing 302. The [equation protocol](../artifacts/generated-results/elkies-k3-curve302-anchor6-mw14-protocol-v1.json)
reconstructs every rational pencil and its degree-one rational zero map.
For each frame, PARI and independent rational LDL enumeration agree on
eight signed roots. Their primitive root lattice is `A2+A1`. Thus the
exact arithmetic generic rank is 14, torsion is trivial and the full
abstract MW determinant is `948/6=158`. Full abstract height matrices
are supplied; Weierstrass-coordinate MW bases are not constructed.

The frozen screen uses primes
`1013,1021,1009,1031,1033,1039,1049,1051`, one worker and a 600-second
invocation cap, saving a checkpoint after each model. It performs 256
generic modular conversion attempts. Nine fail a required conversion
guard and remain UNKNOWN at that prime. They contribute no exclusion.
The [preserved screen](../artifacts/generated-results/elkies-k3-curve302-anchor6-mw14-v1.json)
excludes 86 pencils and leaves rows 4 and 33 unresolved.

Both survivors are then converted exactly over Q(s). Explicit rational
coordinate maps in both directions, a plane cubic, Jacobian Weierstrass
coefficients and the j-map are supplied. Each j-map has degrees 24/21.
Each primitive integral comparison with j302 has degree 24 and no
projective root modulo 61. This closes both inverse equations without
a rational-parameter height bound. These final obstructions use exact
rational polynomials and do not depend on surface good reduction at 61.

At 1013 the 88 normalized binary branch discriminants have degree 43
and are pairwise different. They prove the 88 elliptic fibrations
inequivalent even over the algebraic closure, while all remain on X948.
One separately normalized smooth fibre per pencil agrees with its
generic modular j-map. Those 88 checks validate a separate computational
method at the sampled fibres; they are not an interpolation proof.

## Certificate and replay

The [final certificate](../artifacts/generated-results/elkies-k3-curve302-anchor6-triangles-complete-v1.json)
combines the complete lattice gate, all 88 MW14 inverse results and the
existing six-MW15 exclusion. The earlier six-MW15 proof is retained as
a dependency, not silently replaced.

```sh
sage -python elkies-k3/scripts/verify_curve302_anchor6_triangles.sage
```

The checker rebuilds all 88 frames and rational pencils, recomputes the
modular screens, independently enumerates the complete fixed-anchor
shell, verifies both exact rational inverse witnesses, and reconstructs
the 88 branch signatures and independent fibre checks. The existing
six-MW15 certificate and its input hashes are checked as a dependency.
SageMath 10.9, one worker, 600 seconds. Build also verifies the saved JSON
round trip. Detailed logs and per-model checkpoints remain under
`artifacts/local/elkies-k3/`.

MW14, MW15 and MW16 remain eligible construction targets. This theorem
provides a stopping rule for the fixed P6 height-six triangle route;
it does not complete the parent-recovery objective.
