# Elkies rank-17 K3 model

Expected Weierstrass representation:

    y^2 + a1(t) x y + a3(t) y
      = x^3 + a2(t) x^2 + a4(t) x + a6(t)

JSON coefficient arrays are stored low-degree first:

    [c0,c1,c2,...]

for

    c0 + c1*t + c2*t^2 + ...

Once recovered, store the canonical model as:

    elkies-rank17.json

and generic MW sections separately.
