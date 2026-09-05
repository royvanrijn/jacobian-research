# Regulator obstruction: the remaining four products have arithmetic rank zero

Each of the four remaining product twists has
`E^(d)(QQ(u))={O}` and arithmetic rank zero. No target is
regulator-compatible with rank one, so the explicit section-solving queue
is empty.

| Product | Good primes used | Regulator squareclasses | Classification |
|---|---|---|---|
| `11ee2:0c36e` | 131, 137 | 1006, 182 | rank 0 |
| `0c10b:17a1a` | 131, 137 | 670, 78 | rank 0 |
| `0f82c:025be` | 137, 151 | 2, 6 | rank 0 |
| `11ae6:0f82c` | 131, 137 | 1194, 1 | rank 0 |

Only one additional Frobenius calculation was needed: `0f82c:025be` at
151. Its stored polynomial at 131 has analytic rank **two**, so it cannot
provide a rank-one regulator constraint. All other rows reuse the stored
131/137 polynomials.

Together with the [previous `19bad:083ad` exclusion](R17_PRODUCT_19BAD_083AD_ARITHMETIC_RANK_ZERO_2026-09-05.md)
and the [twelve geometric rank-zero products](R17_ALTERNATE_Q80_ALL17_PRODUCT_TWIST_CLASSIFICATION_2026-09-04.md),
this closes the arithmetic product-character worklist for all seventeen
selected bases. The five former Tate survivors still have geometric rank
`UNKNOWN` in `[0,2]`. The other twelve have geometric rank zero.

## Exact local and leading-term data

Use the same normalization as the single-target proof:

\[
L_p(T)=\det(1-T\operatorname{Frob}_p),\qquad
L_p^*=\left.\frac{L_p(T)}{1-pT}\right|_{T=1/p},\qquad
C_p=\frac{p^3 L_p^*}{\prod_v c_v}.
\]

The eight rows below all have analytic rank exactly one and geometric
fibres `4I0*+24I1`, with smooth infinity and `chi=4`.

| Product | `p` | `L_p^*` | `prod c_v` | `C_p` | Squareclass |
|---|---:|---:|---:|---:|---:|
| `11ee2:0c36e` | 131 | `65929216/2248091` | 1 | 65929216 | 1006 |
| `11ee2:0c36e` | 137 | `584450048/2571353` | 4 | 146112512 | 182 |
| `0c10b:17a1a` | 131 | `175636480/2248091` | 4 | 43909120 | 670 |
| `0c10b:17a1a` | 137 | `5111808/2571353` | 4 | 1277952 | 78 |
| `0f82c:025be` | 137 | `1179648/2571353` | 4 | 294912 | 2 |
| `0f82c:025be` | 151 | `14155776/3442951` | 1 | 14155776 | 6 |
| `11ae6:0f82c` | 131 | `78249984/2248091` | 16 | 4890624 | 1194 |
| `11ae6:0f82c` | 137 | `12845056/2571353` | 4 | 3211264 | 1 |

At each closed branch place, the checker refactors the original residual
two-division cubic and computes `c_v=1+#roots`. Each closed place contributes
once to the Tamagawa product, irrespective of its residue degree. All other
local factors are one. These factorizations also agree with the boundary
permutation factors in the stored complete Frobenius certificates.

The exact nonsquare ratios give short valuation witnesses:

| Product | Ratio | Odd valuation |
|---|---:|---:|
| `11ee2:0c36e` | `C_131/C_137=2012/4459` | `v_7=-3` |
| `0c10b:17a1a` | `C_131/C_137=1340/39` | `v_3=-1` |
| `0f82c:025be` | `C_137/C_151=1/48` | `v_3=-1` |
| `11ae6:0f82c` | `C_131/C_137=597/392` | `v_2=-3` |

The rank-two row at 131 for `0f82c:025be` is retained as
`NOT_APPLICABLE`, with no rank-one regulator squareclass attached. A
hypothetical characteristic-zero section would only force positive rank
there, not rank equality; moreover a rank-two regulator does not determine
an individual section's height squareclass.

## Proof

Apply the [height-preserving specialization argument](R17_PRODUCT_19BAD_083AD_ARITHMETIC_RANK_ZERO_2026-09-05.md#why-this-is-an-unconditional-obstruction)
to each pair of usable primes. Its hypotheses are the same for all four
targets: the discriminant and quartic retain full degree, are squarefree
and coprime, and the good surface reductions preserve `U+4D4`.
Specialization preserves the intersection form and projection away from
this geometric trivial lattice, hence preserves the positive Shioda height
of a nontorsion section, including sections with denominators. The
torsion-free `48I1` double cover rules out any nonzero torsion section.

If a nonzero `P` over `QQ(u)` existed, its reduction at each usable prime
would be nontorsion and force algebraic rank equal to analytic rank one.
Rank equality implies refined BSD and finite square-order `Sha` for these
elliptic curves over finite function fields. This is the theorem in
[Ulmer, Theorems 6.2.6 and 6.3.1 and Proposition 6.3.3](https://math.stanford.edu/~conrad/BSDseminar/refs/Ulmer.pdf);
the rational origin makes all indices and periods in the square-order
formula equal to one.

With heights divided by `log p`, BSD gives

\[
C_p=\operatorname{Reg}(E_p)\,
      \frac{\#\Sha(E_p)}{\#E_p(\mathbf F_p(u))_{\rm tors}^{\,2}}.
\]

The specialization of `P` is an integer multiple of the rank-one
generator modulo torsion. Consequently `h(P)/C_p` is a nonzero rational
square at each usable prime. Every ratio in the preceding table is
nonsquare, a contradiction. This proves arithmetic rank zero at every
height, not merely absence in a bounded polynomial box.

Thus the anti-invariant group, integral character glue, point-Kummer image,
and Tate quotient vanish for each of these four twists. All rational
height-eight, height-ten, and higher section boxes are empty. This does
not compute the full two-Selmer groups or assert that `Sha[2]` vanishes.

## Certificate, controls, and replay

The [aggregate certificate](../artifacts/generated-results/elkies-k3-r17-product-regulator-sweep-v1.json)
pins the model, shortlist, prior height gate, all nine Frobenius inputs,
their audited inputs, and the retained control files. Its
`explicit_section_solving_queue` is `[]`.

The [generalized checker](scripts/certify_r17_product_regulator_sweep.sage)
shares its rank-one local/BSD calculation with the original single-target
checker. The original certificate replays unchanged. Independent Magma
`AnalyticInformation(E,L)` controls agree with all eight `C_p` values;
they check local/BSD normalization using the supplied certified
L-polynomials. Exact jobs and successful raw XML, plus the new toric
input/output and log, are preserved in
[`r17-product-regulator-sweep-controls/`](../artifacts/generated-results/r17-product-regulator-sweep-controls/).

```bash
sage -python elkies-k3/scripts/certify_r17_product_regulator_sweep.sage --check
sage -python elkies-k3/scripts/certify_r17_product_regulator_sweep.sage --self-test
sage -python elkies-k3/scripts/certify_r17_product_19bad_083ad_rank_zero.sage --check
```

These are local, cheap replays. The self-test checks square rescaling,
nonsquare valuation witnesses, exclusion of rank-two data from the
rank-one comparison, and rejection of a corrupted moment.
`--export-magma` regenerates the eight control jobs without making network
requests. Omitting `--check` regenerates the aggregate certificate.

The only new cohomological calculation is reproducible with:

```bash
bash elkies-k3/scripts/run_r17_product_toric_frobenius_extra_prime.sh \
  'alternate-orbit-0f82c:alternate-orbit-025be' 151
```

Its independent first-two-moment audit and complete degree-28 certificate
both pass. The toric calculation took about 261 seconds; routine regulator
replay does not rerun it. No section solver or Selmer search was launched.

<!-- status-consumer: EC-K3-R17-PRODUCT-REGULATOR-OBSTRUCTION-SWEEP f86dead53d55babe -->
