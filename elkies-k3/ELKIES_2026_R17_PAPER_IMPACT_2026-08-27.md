# Elkies 2026 rank-17 paper: impact on the H3 -> R17 reconstruction

Date: 2026-08-27

Source: Noam D. Elkies, *An elliptic K3 surface X/Q(t) with Mordell-Weil rank 17, I: Formulas for X and base changes of ranks 18 and 19*, arXiv:2608.25406v1, submitted 2026-08-26.

https://arxiv.org/abs/2608.25406

## Executive conclusion

This paper changes the role of the final R17 lift.

Before the paper, the repository had independently identified `rank17_gram.txt` as the likely Elkies essential lattice and was trying to reconstruct an equation all the way from H3 through a marked neighbor chain.  Elkies now publishes the exact rootless rank-17 Weierstrass model and an explicit height Gram for 17 integral sections.

Therefore the final equation is no longer an unknown object that must be *discovered* by the neighbor chain.  It is now an independent, authoritative endpoint oracle.  The remaining reconstruction problem should be split into:

1. construct enough of the H3-side neighbor chain to prove that the selected final fibration is the pinned R17 `U`-embedding;
2. identify the final base coordinate / Weierstrass scaling with the published model;
3. use the published equation as the canonical endpoint instead of re-reconstructing its huge coefficients from scratch.

This does **not** make the intermediate q8/orbit376 proof gaps disappear.  A random semistable `4A1` child is still not automatically the marked q8/orbit376 child.  But it substantially lowers the cost and risk of the final q12 closeout.

## What the paper now gives exactly

Elkies publishes

```text
y^2 = x^3 - 27*S(t)*x + (27/4)*T(t)
```

with `deg(S)=8`, `deg(T)=12`, and gives one complete integral point `P1=(x1,y1)`.  He also publishes the x-coordinates of 16 further integral sections and a 17x17 height Gram matrix of determinant

```text
948.
```

The paper states that the Mordell-Weil lattice has exactly `1311` +/- pairs of norm-4 vectors.  These are precisely the two strongest fingerprints already used in this repository to identify `rank17_gram.txt`: determinant 948 and 1311 pairs of minimal height-4 sections.

A new verifier on this branch,

```text
elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
```

checks the published equation, rootless discriminant, Gram determinant/minimal shell and then asks Sage/PARI `qfisom` for an **integral isometry** between the published Gram and `data/lattice/rank17_gram.txt`.  A PASS therefore upgrades the old high-confidence identification to an explicit basis transformation.

## The most useful new computational facts

### 1. Rootless R17 has a trivial height formula

Because the published fibration has no reducible fibres, Elkies proves that if

```text
x = N/D
```

with homogeneous degrees `d` and `d-4`, then the canonical height is exactly `d`.

Thus for the final R17 model:

```text
height(P) = 4 + 2*(P.O).
```

There are no local correction terms to resolve.  This means final-section recognition can be done from denominator degree alone.

### 2. Integral-section pairings are just gcd degrees

For two integral sections `P=(x,y)` and `P'=(x',y')`, Elkies gives

```text
<P,P'> = 2 - deg gcd(x-x', y-y')
```

with the projective contribution at infinity included.

This is much cheaper than repeated function-field group-law height calculations.  Once a q12-derived endpoint is mapped to the published coordinates, a large part of its MW marking can be certified by polynomial gcds.

### 3. The published basis was deliberately chosen for cheap chord arithmetic

For `P2,...,P16`, Elkies chose

```text
<P1,Pi> = -2.
```

Consequently

```text
mi = (yi-y1)/(xi-x1)
```

is only a quadratic polynomial.  The paper publishes all fifteen `mi`, with coefficients dramatically smaller than the original section coordinates.

This is directly relevant to our compiler philosophy: **chord data can be far smaller than section coordinates**.  It supports the q8/o376 post-collision approach, although the current branch's chord sign convention must first be corrected as documented in `Q8O376_RR_BRANCH_AUDIT_2026-08-26.md`.

It also suggests a new possible final-edge strategy: compile the inverse q12 neighbor from the published R17 side using the 17 integral sections / cheap quadratic chords, and meet the H3-side 4A1 model in the middle.

## Recommended change to q12/orbit5867 closeout

The current promoted lattice route is

```text
... -> 4A1/MW13 --q12/o5867--> pinned rootless R17.
```

Previously we expected the final q12 compiler to manufacture a huge exact R17 equation and then prove that its essential lattice was pinned R17.

That is now unnecessary.

### Preferred endpoint procedure

After obtaining a physically marked P1229-zero `4A1` equation:

1. compile only enough of q12/o5867 to obtain the rootless child and the induced base function `U`;
2. compute the child's `j(U)`;
3. match it to the published

   ```text
   j_pub(t) = 1728 * 4*A_pub(t)^3 / (4*A_pub(t)^3 + 27*B_pub(t)^2)
   ```

   by solving for a rational Mobius transformation

   ```text
   t = (a*U+b)/(c*U+d);
   ```

4. after the base transformation is known, solve only for the standard Weierstrass scaling/twist needed to identify the child with Elkies's published model;
5. certify that this endpoint is the pinned R17 `U`-embedding using the already-certified NS transport plus the new published-Gram -> pinned-Gram integral isometry.

This should be much cheaper than reconstructing the final ~800 digits coefficient-by-coefficient from modular q12 calculations.

### Strong normalization/check oracles

The paper publishes four high-rank specializations:

```text
rank >= 25: t = -2/377
rank >= 26: t = -308/251
rank >= 27: t = 2456/135
rank >= 28: t = -9529/5471
```

These are useful independent checks after finding the Mobius base transformation.  In particular the historical rank-28 fibre supplies a very distinctive endpoint normalization fingerprint.

## A second strategy worth testing: reverse q12 from published R17

Because all 17 published generators are integral and many pair with `P1` by `-2`, the published endpoint is exceptionally friendly for explicit divisor/chord calculations.

The q12 lattice certificate already tells us, in pinned R17 coordinates, which old `4A1` fibre and zero must be recovered when running the final edge backwards.  Once the new published-to-pinned basis isometry is available, transport those classes into Elkies's published section basis.

Then try to realize the inverse q12 fibre directly using:

- integral sections of height 4;
- their degree-2 chord slopes `mi`;
- small sums/differences of the published sections;
- divisor-first RR rather than recovering a large forward q12 horizontal.

If the inverse old-fibre class decomposes into a small set of published integral sections/chords, this could make the final edge dramatically easier than the forward q12 compiler.  It also gives a strong meet-in-the-middle test: the inverse construction should land on a `4A1/MW13` model birationally equivalent to the H3-side q8/o376 child.

## Proposition 8 / quadratic sections

The paper's final section gives another useful lattice criterion, although it is not needed for the immediate H3 -> R17 lift.  For the rootless R17 fibration, a divisor associated with a trace `tau` gives a rational quadratic section precisely when there is no MW point `P` with

```text
h(tau - 2P) = 6.
```

This characterizes quadratic base changes raising the rank to at least 18.  Elkies reports 39120 mod-2 MW cosets with no norm-6 representative and exhibits explicit quadratic covers, including one pair whose compositum yields rank at least 19 over a positive-rank elliptic base.

This is a natural follow-up once the R17 endpoint is fully attached: the repository can reproduce and then systematically search these mod-2 cosets using the already-developed lattice tooling.

## Interaction with the q8/o376 branch audit

The paper reinforces one idea in this branch and weakens one motivation:

- **reinforces:** using low-degree chord functions instead of explicit high-degree sections is exactly the kind of compression Elkies himself uses for the R17 basis;
- **weakens:** there is no reason to reconstruct the final R17 coefficients through CRT merely to discover them; they are now published.

The blocking q8/o376 audit items remain unchanged:

1. fix the chord sign convention;
2. prove the exact marked q8 fibre class, not merely the `4A1` ADE fingerprint;
3. prove the complete `12 -> 4 -> 2` RR rank;
4. enforce artifact/script provenance;
5. attach enough child marking to run or meet the final q12 edge.

## Immediate action

Run:

```bash
~/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage
```

If `PASS_EXACT_PUBLISHED_R17_IS_PINNED_R17`, retain the emitted integral basis transformation as a new endpoint certificate.  Then prioritize a **published-endpoint q12 meet-in-the-middle** over blind reconstruction of the huge final equation.
