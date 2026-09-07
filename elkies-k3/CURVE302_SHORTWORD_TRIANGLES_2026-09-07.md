# Ten inequivalent MW14 triangle fibrations exclude302

Authority: `EC-K3-CURVE302-TEN-SHORTWORD-TRIANGLES`.

The frozen short-word dictionary on the determinant948 K3 supplies ten
explicit rational triangle pencils of exact generic arithmetic rank14.
**All ten exclude302 at every rational parameter.** Their j-maps have
pairwise different branch divisors, proving the fibrations inequivalent
even over Qbar. They remain fibrations on one K3 surface.

This completes this finite construction-and-inverse experiment. It does
not classify all degree-three fibrations, all section configurations, or
all parents of302. The alternative-parent objective remains open.

The [two original triangle proofs](CURVE302_TRIANGLE_MW14_2026-09-07.md)
remain necessary construction and transport controls. This note is the
active summary of the larger dictionary and its completed inverse tests.

## Frozen dictionary and exact pencils

Use the seventeen one-based sections of the certified
[11952 presentation](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
The dictionary consists of all289 basis words with support at most two,
coefficients ±1 and first nonzero coefficient positive. It contains65
height-six words and67 height-four words.

For every pair of height-six words P,Q with pairing ±3, choose the sign
of Q so the pairing is3. Then `O,P,Q` meet pairwise once. Normalize the
three unoriented edge vectors `P,Q,P-Q` to remove translations and signs,
retaining the first representative of each of127 edge keys. Admit that
representative if one of the signed height-four words supplies a section
Z with `(O+P+Q).Z=1`. The exact roster contains123 admitted configurations;
four representatives are outside this frozen zero-section rule.
This edge normalization alone is not a complete fibration deduplication.

The complete root census selects ten with generic rank14. For all ten,
PARI and independent rational LDL enumeration agree on eight signed
roots; their root lattice is primitive `A2+A1`. The full rational NS rank
is19. Hence every selected pencil has exact generic arithmetic rank14,
zero torsion and full abstract MW determinant `948/6=158`.

The [protocol](../artifacts/generated-results/elkies-k3-curve302-shortword-triangles-protocol-v1.json)
records all123 rows and all ten exact rational pencils before the new
target tests. It contains their fibre and zero classes, frame matrices,
root vectors, rational-function pencil coefficients, degree-one zero maps
and implicit trigonal equations. The pole-cancellation proof from the
original triangles is checked for every selected word overQ. In particular,
the three finite intersection parameters are distinct, and the new base
map has no vertical poles. Both old pencils are recovered exactly as
transport controls. Row0 chooses another valid zero in the new dictionary;
the original P7 zero map is checked separately, and its j-map is unchanged.

## Completed inverse tests

| Row | Triangle `O+P+Q` | Exclusion witness |
|---:|---|---:|
| 0 | `O+P4+P11` | 1013, retained control |
| 12 | `O+P6+P8` | 1013, retained control |
| 17 | `O+P6+(P16-P17)` | 139, exact QQ inverse polynomial |
| 30 | `O+P10+P12` | 1009 |
| 33 | `O+P10+(P4-P1)` | 1013 |
| 56 | `O+P16+(P6-P1)` | 1021 |
| 84 | `O+(P3-P13)+(P14-P17)` | 1013 |
| 98 | `O+(P6-P14)+(P6-P15)` | 1013 |
| 108 | `O+(P9-P11)+(P17-P16)` | 1013 |
| 120 | `O+(P12-P13)+(P12-P14)` | 1021 |

The first stage freezes primes1013,1021,1009, at most64 sample parameters
per prime,300 seconds per probe, one worker and1200 seconds total.
Fourteen supervised probes all finish successfully, in585.7 seconds total.
Each uses52 genus-one fibre normalizations,49 interpolation samples and
three additional checks. Exact reduction of the QQ pencil and its zero
map is checked before each probe. Degree24 is retained throughout.

The preserved [first-stage record](../artifacts/generated-results/elkies-k3-curve302-shortword-triangles-v1.json)
therefore contains **nine exclusions and one UNKNOWN**, rather than
silently counting row17 as excluded. Row17 survives at all three primes:
its necessary residue sets have sizes1,2,3 respectively. These residues
are not rational parameters or evidence of an actual302 fibre.

## Exact generic conversion closes row17

The [generic conversion certificate](../artifacts/generated-results/elkies-k3-curve302-triangle-generic-conversion-v1.json)
constructs a plane cubic over Q(s), verifies both birational directions,
and supplies its explicit Jacobian Weierstrass coefficients and j-map.
It uses no specialized-fibre normalization or interpolation.

After making the trigonal equation monic in x, write

\[
f=x^3+a(u)x^2+b(u)x+c(u),\qquad
\deg_u(a,b,c)=(4,8,12).
\]

Its discriminant is `q(u)^2*h(u)`, with degrees9 and6 respectively; h is
squarefree and coprime to q. The repeated root modulo q is

\[
r=\frac{9c-ab}{2(a^2-3b)}\pmod q.
\]

Set `alpha=a+r mod q`, `beta=-r(a+2r) mod q`. Polynomial row reduction
of the candidate basis

\[
1,\quad x,\quad \frac{x^2+\alpha x+\beta}{q}
\]

with weights `(0,4,8)` finds a new coordinate
`v=g0(u)+g1(u)*x+g2(u)*x^2`. Its characteristic polynomial has the form

\[
v^3+A_1(u)v^2+A_2(u)v+A_3(u),\qquad \deg A_i\leq i.
\]

Homogenization is a plane cubic. Crucially, the proof does not rely on
trusting this basis-finding procedure: the certificate gives g0,g1,g2 and
an explicit inverse `x=h0(u)+h1(u)*v+h2(u)*v^2`. Companion-matrix identities
verify both the cubic relation and the inverse exactly over Q(s)(u).
The new cubic discriminant is squarefree of u-degree6; the squared basis
determinant also verifies its relation to the old discriminant.

The standard ternary-cubic invariants give an explicit Jacobian
`Y^2=X^3+A(s)*X+B(s)`, as documented in
[Sage's Weierstrass-form construction](https://doc.sagemath.org/html/en/reference/schemes/sage/schemes/toric/weierstrass.html).
The source pencil already has a Q(s)-point, so its smooth genus-one fibre
is isomorphic to this Jacobian. No Weierstrass-coordinate MW basis is
asserted here. The computed j-map has numerator/denominator degrees24/21.

The primitive integral polynomial comparing this exact j-map with j302
has degree24 and maximum coefficient size15,315 bits. Its reduction at139
has no finite root and a nonzero leading coefficient. Thus it has no
projective rational root. **This closes row17 at every rational parameter.**
Because the QQ inverse polynomial itself is known, this final witness
does not require a separate surface-good-reduction claim at139.

## Independent check of all modular j-maps and deduplication

The generic plane-cubic method is also run over Fp(s) for every saved
map, including the two controls. All sixteen entire rational j-maps agree
exactly with the earlier fibre-normalization method, covering832 saved
sample values. This provides a different computational method, not merely
another interpolation of the same values. It shares the certified
rational-pencil input and does not independently reconstruct the K3.

The [branch certificate](../artifacts/generated-results/elkies-k3-curve302-ten-triangle-branches-v1.json)
computes the binary discriminant of `n_h-J*d_h` for all ten maps at1013.
All ten normalized degree43 polynomials differ. Their common factors are
`J^16*(J-1728)^12`; the monic residual degree15 polynomials distinguish them.
The degree46 discriminant bound permits47 interpolation values plus three
checks. The nonaffine PGL2 control and homogeneous leading-term-loss
convention are inherited from the original branch-separation checker.
Thus no two of these ten fibrations are duplicate presentations, even
over Qbar. This does not imply ten different K3 surfaces.

## Replay and retained evidence

```sh
sage -python elkies-k3/scripts/verify_curve302_shortword_triangles.sage
```

This composite replay rebuilds the dictionary, rank checks and QQ pencils,
checks the saved modular inverse polynomials, verifies the stored QQ
birational conversion, independently reconstructs all sixteen generic
modular j-maps, and checks all ten branch signatures. It performs no new
point or parameter search. The software is SageMath10.9.

Detailed logs and per-sample checkpoints are retained under
`artifacts/local/elkies-k3/curve302-shortword-triangles-v1/` and its parent
directory. The original first-stage outcome remains9 exclusions/1 unknown;
the generic conversion certificate supplies the tenth exclusion. Initial
interface failures, the timed-out general symbolic-normalization attempt,
and the successful direct algebraic pilots remain local evidence. No
failed or timed-out computation contributes an exclusion.

The next parent candidate must come from outside this completed dictionary
or a different construction pattern. The full302 parent objective remains
unchanged and unresolved.
