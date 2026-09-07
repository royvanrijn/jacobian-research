# Attack the MW3 node directly

The selected node has:

    ADE = A10 + A2 + A1^2
    MW rank = 3
    torsion = 0
    |disc(NS)| = 948
    disc(root lattice) = 11*3*2^2 = 132
    det(MW height lattice) = 948/132 = 79/11

Shimada's classification lists A10 + A2 + 2A1 with trivial torsion only.

## Step 1: recover the exact MW height Gram

    sage elkies-k3/scripts/recover_mw3_height_lattice.sage

This projects the full frame lattice through the root lattice, includes finite
glue exactly, and returns the rational 3x3 MW height Gram.

## Step 2: enumerate possible Kodaira realizations

    sage elkies-k3/scripts/describe_mw3_fiber_models.sage

The simplest/semistable realization is:

    I11 + I3 + I2 + I2 + 6 I1

Normalize reducible fibers to infinity, 0, 1, lambda. Then for short
Weierstrass coefficients A(t), B(t):

    deg A <= 8
    deg B <= 12
    Delta = -16(4A^3+27B^2)

and I11 at infinity forces deg Delta = 13. In the semistable case:

    Delta = c*t^3*(t-1)^2*(t-lambda)^2*R6(t).

The reducible-fiber configuration alone is a 4-dimensional family. The three
MW sections (with the exact height Gram from step 1) cut out the target
1-dimensional rank-19 Shimura family.

Next reconstruction step:
- derive component labels of a reduced MW basis at I11/I3/I2/I2 from the
  discriminant/glue data;
- translate those component labels into local valuations for x(t),y(t);
- solve the resulting three-section Weierstrass system.
