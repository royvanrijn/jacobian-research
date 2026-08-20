# Three MW3 reconstruction branches

Run all three initial jobs:

    bash elkies-k3/scripts/start_three_mw3_branches.sh

Logs:

    artifacts/local/elkies-k3/three-branches/a10.log
    artifacts/local/elkies-k3/three-branches/e6.log
    artifacts/local/elkies-k3/three-branches/a6.log

## A10 branch

ADE:

    A10 + A2 + A1^2

Current preferred fibers:

    I11 + I3 + I2 + I2 + 6 I1

This branch starts the existing local-Tate P1 modular probe.

## E6 branch

ADE:

    E6 + A3^2 + A1^2

First model to try:

    IV* + I4 + I4 + I2 + I2 + 4 I1

This is attractive because putting IV* at infinity directly forces:

    deg A <= 5
    deg B <= 8

so it avoids the long I11 discriminant-cancellation staircase.

The initial job computes its exact MW height lattice/glue.

## A6/A4 branch

ADE:

    A6 + A4 + A1^4

First semistable model:

    I7 + I5 + 4 I2 + 4 I1

The initial job computes its exact MW height lattice/glue.

After these finish, use the two new MW height Grams to enumerate component
labels exactly as was done for the A10 branch. Then build local-Tate section
scaffolds for each branch.
