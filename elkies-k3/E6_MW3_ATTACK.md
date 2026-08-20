# E6 MW3 direct attack

Node:

    ADE = E6 + A3^2 + A1^2
    MW rank = 3

Exact reduced MW height Gram:

    (1/12) * [[23,-10,-8],[-10,23,1],[-8,1,23]]

Preferred Kodaira model:

    IV* + I4 + I4 + I2 + I2 + 4 I1

The IV* fiber is placed at infinity. In short Weierstrass form this directly
forces:

    deg A <= 5
    deg B <= 8

which is substantially simpler than the A10/I11 branch.

Run both preparation jobs:

    bash elkies-k3/scripts/start_e6_attack.sh

or individually:

    sage elkies-k3/scripts/enumerate_e6_component_labels.sage
    sage elkies-k3/scripts/build_e6_mw3_fiber_scaffold.sage

The component labels use the correction tables:
- IV*: nonzero diagonal 4/3; distinct nonzero pair 2/3
- I_n: standard A_(n-1) inverse-Cartan corrections.

If these collapse to a small symmetry orbit, build section ansatzes next.
