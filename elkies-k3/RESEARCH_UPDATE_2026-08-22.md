# Research update — 2026-08-22

## Q80 secondary route: low-q corridor found

The canonical Q80 six-step rootless certificate in
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md)
remains unchanged and remains the only Q80 route certified to MW17/rootless.

A new secondary corridor has now been found and checked exactly through MW8:

```text
E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5
 --q6--> D7+D4/MW6
 --q4--> A6+A4/MW7
 --q4--> A6+A3/MW8.
```

The new q6 is a non-shortest root correction in the **same MW coset** as the
old pinned third q12 move. After chamber reduction it has old-fibre degree two
and the same marked horizontal section, replacing the marked-trisection gate
by a marked-chord/binary-quartic gate.

At CM24 the q6 has an exact characteristic-zero equation over
`QQ(sqrt(-6))`. The exact new base is

```text
V = (z_Q - 162*s + (39/4)*s*(W+27/2)) / (W+27/2)^2,
s^2=-6,
```

and the child fibres are

```text
I4* + I2* + 2 I2 + 2 I1
ADE = D8+D6+2A1
MW = 2.
```

The preferred next q4 is Weyl orbit 424. Its CM24 marked section specializes
to the exact 2-torsion point

```text
T = (1944*s*V^3 + 12150*V^2 - 4401*s*V - 3036, 0),
```

so its new base becomes

```text
U = Y / ((X-T_x)*(V-13*s/18)).
```

The resulting CM24 child is exceptionally clean:

```text
2 I8 + 8 I1
ADE = 2A7
MW = 4.
```

Exact Weyl-orbit scans show the same local pattern in both checked corridor
frames:

- q2: all loops;
- q3: all loops;
- q4 `(1,4)`: all loops;
- q4 `(2,2)`: first productive shell, dropping root rank by one.

From `A6+A4/MW7`, the retained rank-9 continuation is orbit 1222:

```text
A6+A4/MW7 --q4 (2,2)--> A6+A3/MW8.
```

At CM24 this specializes to `2A6+3A1`; its marked section remains non-torsion
with `P.O=1`, height `25/8`, and a vertical correction supported on one A7
fibre with only three nonzero coefficients. The same CM24 root signature also
occurs for the old pinned q12 child, but equality of the two elliptic
fibrations has **not** been proved.

Full details, exact vectors, equations, scan counts, and open questions are in
[`Q80_LOW_Q_ALTERNATE_2026-08-22.md`](Q80_LOW_Q_ALTERNATE_2026-08-22.md).
Machine-readable records and replay scripts are:

```text
data/fibrations/kumar_q80_lowq_alternate_prefix.tsv
data/fibrations/kumar_q80_lowq_weyl_summary.tsv
scripts/verify_q80_lowq_alternate_prefix.sage
scripts/verify_q80_lowq_cm24_equations.sage
```

## Reusable compiler update

The Q80 third-q12 module work also identified that the exact neighbour
compiler's finite-quotient stack was unnecessarily restricted to `QQ`.
A compatibility layer now makes the local quotient/module intersection
machinery field-generic while preserving the historical `QQ` default:

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage
```

This supports exact number-field module intersections such as the Q80 CM24
compositum `QQ(sqrt(-6),sqrt(-3))`.

## Current priority

The live fallback branch is now

```text
D7+D5/MW5
 --q6 escape--> D7+D4/MW6
 --q4 orbit424--> A6+A4/MW7
 --q4 orbit1222--> A6+A3/MW8.
```

Next work should first determine whether this q4 root-peeling corridor
continues cheaply from `A6+A3/MW8`, and whether the orbit-1222 CM24 child is
actually the same `2A6+3A1` fibration reached by the old q12 route.

The H3 route remains the primary source-polarization route. Q80 remains a
secondary route until this alternate corridor is continued to rootless and
its generic characteristic-zero marked sections are algebraized.
