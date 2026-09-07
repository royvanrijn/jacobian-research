# Umbral complement-orbit pilot (2026-09-02)

## Result

The proposed experiment gives a useful negative control before an all-Niemeier
search.  The six requested embeddings have the same Niemeier provenance and
the same effective umbral stabilizer type:

\[
 X=A_7^2D_5^2,\qquad G^X\cong \operatorname{Dih}_4,
 \qquad \operatorname{im}(\operatorname{Stab}_N(K)\to G^X)
 =\{1A,2A\}.
\]

In every case the induced complement action is exactly
\(\{+I_M,-I_M\}\), with class \(2A\) acting as \(-I_M\).  Consequently it
acts trivially on \(M/2M\).  All rational bisection cosets are singleton
orbits, so this action does **not** distinguish the 41,421 bisections of
`NS0032-F011` from the 39,120 of published R17.  The larger count is a genuine
lattice statistic, but it is not explained by different umbral orbit types in
this pilot.

This is an exact finite computation for stabilizers, norm-four vectors, and
degree-two cosets.  Degree-three results below are exact only inside the stated
deterministic samples.  No equality with an umbral module, and no new
geometric or arithmetic theorem, is claimed.

## Why all six have the same ambient

The current foundry catalogue was deliberately a finite
`one-root-control-shell-v1` experiment inside \(N(2A_7+2D_5)\), seeded by the
published R17/Q80 construction.  It was designed to compare auxiliary
embeddings in one controlled ambient, not to classify embeddings across all
23 rooted Niemeier lattices.  Thus the shared \(X\) is a search-design
restriction, not evidence that this Niemeier lattice is uniquely relevant.

For the shorthand `NS0024`, this pilot uses `NS0024-F005`, the frame selected
by the existing foundry multisection work.

## Exact and sampled counts

Here \(N_4\) is the number of signed norm-four vectors.  The \(d=2\) column is
the complete rational-bisection coset count.  The final two columns give the
size of the invariant \(M/3M\) sample and its number of qualifying rational
cosets.

| construction | \(\det M\) | \(N_4\) | unoriented norm-4 orbits | rational \(d=2\) cosets/orbits | sampled \(d=3\) | rational in sample / orbits |
|---|---:|---:|---:|---:|---:|---:|
| published R17 (`NS0001-F001`) | 948 | 2,622 | 1,311 | 39,120 / 39,120 | 513 | 64 / 32 |
| Q80 (`NS0001-F002`) | 948 | 2,626 | 1,313 | 39,147 / 39,147 | 513 | 62 / 31 |
| `NS0024-F005` | 950 | 2,632 | 1,316 | 39,012 / 39,012 | 513 | 80 / 40 |
| `NS0032-F011` | 1,124 | 2,430 | 1,215 | **41,421 / 41,421** | 513 | 70 / 35 |
| `NS0028-F005` | 1,132 | 2,394 | 1,197 | 41,376 / 41,376 | 513 | 78 / 39 |
| `NS0033-F026` | 1,088 | 2,446 | 1,223 | 40,912 / 40,912 | 513 | 56 / 28 |

For every row:

- the literal stabilizer in the chosen chamber-preserving eight-element
  section of \(G^X\) has order 1;
- the full ambient stabilizer \(\operatorname{Stab}_{\operatorname{Aut}(N)}(K)\)
  has order 3,840;
- its image in \(G^X\) has order 2 with classes \(1A,2A\);
- its induced image in \(\operatorname{Aut}(M)\) is \(\{+I,-I\}\);
- signed norm-four vectors form two-element \(\{v,-v\}\) orbits;
- unoriented norm-four pairs and all rational \(d=2\) cosets are fixed;
- qualifying nonzero sampled \(d=3\) cosets form two-element
  \(\{c,-c\}\) orbits.

The distinction between the first two stabilizers matters.  An outer class
can stabilize \(K\) only after multiplication by a Weyl element, so merely
intersecting the stored chamber-preserving section with the literal
stabilizer misses class \(2A\).  The full computation enumerates compatible
\(\operatorname{Aut}(K)\times\operatorname{Aut}(M)\) pairs and tests whether
they extend across the primitive gluing \(K\oplus M\subset N\).

Across these six preselected examples, \(N_4\) and the rational-bisection
count are strongly anticorrelated (Pearson \(-0.993\), Spearman \(-0.886\)).
With \(n=6\), one ambient, and a selected sample, this is only an exploratory
diagnostic.  It does suggest that the richer bisection spectrum is not simply
"more short vectors": the richer complements have fewer norm-four vectors.

## Twining convention and comparison

There is no canonical scalar \(\operatorname{tr}(g\mid v)\) attached to one
lattice vector.  The concrete permutation-character version is the
fixed-point theta series

\[
 \Theta^{\mathrm{fix}}_{M,g}(q)
 =\sum_{\substack{v\in M\\gv=v}}q^{(v,v)/2}.
\]

For the only nonidentity class visible here, \(2A=-I_M\), the coefficient at
\(q^2\) is zero, while the identity coefficient is \(N_4\).  On the rational
subset of \(M/2M\), however, the permutation character is

\[
 \chi_{d=2}(1A)=\chi_{d=2}(2A)=N_2,
\]

because \(-c=c\pmod {2M}\).  On the sampled qualifying subset of \(M/3M\),
\(2A\) has no fixed points.

The standard lambency-eight umbral tables give

| datum | \(1A\) | \(2A\) | \(2BC\) | \(4A\) |
|---|---:|---:|---:|---:|
| \(H_{g,1}\), discriminant 31 | 2 | 2 | -2 | 2 |
| \(H_{g,2}\), discriminant 28 | 4 | -4 | 0 | 0 |

These are comparison controls, not matched gradings: no canonical map from
those mock-modular coefficients to the norm-four, \(d=2\), or sampled
\(d=3\) shells has been established.  The group, conjugacy-class labels, and
trace controls are from Cheng--Duncan--Harvey,
[Umbral Moonshine and the Niemeier Lattices](https://arxiv.org/abs/1307.5793),
especially Tables 2, 18, and 38--44.

## Actionable next experiment

The present data say where to broaden the experiment.  Enumerate auxiliary
embedding orbits across all rooted Niemeier ambients, compute their full
ambient stabilizer images, and prioritize images containing \(2B\), \(2C\),
or \(4A\)-type component permutations rather than only scalar \(-I\).  Those
actions can remain nontrivial modulo 2 and therefore can produce meaningful
orbit-resolved bisection-coset characters.  For each such orbit, record the
full fixed-point distribution on the rational subset of \(M/2M\) before using
any total-count correlation as a search heuristic.

Only that cross-\(X\), cross-orbit dataset can test the proposed claim that
distinguished umbral orbit types systematically predict richer complements.

## Reproduction

Run:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/analyze_lattice_foundry_umbral_orbits.sage \
  --d3-orbit-seeds 512 --pari-stack-gb 4
```

This writes
`artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json`.
Each target runs in an isolated child process so that PARI shell-enumeration
memory is released between frames.  The artifact pins the SHA-256 hashes of
all three input catalogues and stores the ambient and induced action matrices,
fixed counts, and orbit histograms.
