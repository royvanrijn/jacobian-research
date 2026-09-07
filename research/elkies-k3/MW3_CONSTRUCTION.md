# MW3 explicit-construction scaffold

Target fibration:

- ADE = A10 + A2 + A1^2
- preferred Kodaira realization = I11 + I3 + I2 + I2 + 6 I1
- MW rank = 3
- torsion = 0

Exact reduced MW height Gram:

    (1/66) * [[79,17,-1],[17,106,19],[-1,19,259]]

Canonical component labels (I11,I3,I2,I2):

    P1 = (2,1,0,1),  P1.O = 0
    P2 = (6,2,1,1),  P2.O = 1
    P3 = (10,2,0,1), P3.O = 1

Mutual intersections:

    P1.P2 = 1
    P1.P3 = 2
    P2.P3 = 2

## Stage A: fiber family

Run:

    sage elkies-k3/scripts/build_mw3_fiber_scaffold.sage

Short Weierstrass normalization:

    y^2 = x^3 + A(t)x + B(t)
    deg A <= 8
    deg B <= 12
    A8 = -3
    B12 = 2

I11 at infinity is imposed by deg Delta <= 13.
I3 is at 0, I2 fibers at 1 and lambda.

This produces 21 variables / 17 equations, so expected dimension 4.

## Stage B: add sections

Run:

    sage elkies-k3/scripts/build_mw3_section_system.sage --stage p1
    sage elkies-k3/scripts/build_mw3_section_system.sage --stage p12
    sage elkies-k3/scripts/build_mw3_section_system.sage --stage all

P1.O=0:

    x1=X1(t), deg X1<=4
    y1=Y1(t), deg Y1<=6

P2.O=P3.O=1:

    z_i=t-r_i
    x_i=X_i/z_i^2, deg X_i<=6
    y_i=Y_i/z_i^3, deg Y_i<=9

After clearing denominators:

    Y_i^2 = X_i^3 + A X_i z_i^4 + B z_i^6.

## Stage C: modular probe

Start with P1 only:

    python3 elkies-k3/scripts/run_mw3_construction_probe.py       --stage p1 --p 101 --threads 8 --timeout 300

Do not jump directly to the full system if P1 is already computationally hard.

The component labels are deliberately kept as a later local-Tate filter instead
of being expanded into the first global polynomial system.
