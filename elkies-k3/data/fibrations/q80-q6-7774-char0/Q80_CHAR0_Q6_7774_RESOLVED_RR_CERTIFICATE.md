# Q80 q6_7774 — exact characteristic-zero resolved RR certificate

Status: **PASS_EXACT_Q6_7774_RESOLVED_RR**

Field: `QQ(sqrt(-3))`.

## Selected old fibres

- I7: `9/76*j - 135/76` (reduces to `V=6` at `j -> 17 mod 73`)
- I2: `-15471/5668*j + 18765/5668` (reduces to `V=5`)

## Resolved Riemann–Roch space

Ambient marked-chord space:

    (1, V, V^2, m),  m=(y+P3y)/(x-P3x)

The selected I7/A6 support was resolved torically in characteristic zero.
Its component restrictions collapse to one exact connected quotient row:

    (1, 9/76*j - 135/76, -1215/2888*j + 8991/2888, -5147030574667159309232932999492283490/575440030160433757457430217*j + 36482341338745219526761278033990297165/1150880060320867514914860434)

with orientation `direct`, support `(1, 4, 5, 6)` and quotient coordinate

    c7 = -5147030574667159309232932999492283490/575440030160433757457430217*j + 36482341338745219526761278033990297165/1150880060320867514914860434.

The I2/A1 exceptional restriction is the intrinsic nodal value

    c2 = 52910815220669667719436421144608397354905/3200608605869260877638343770033*j - 135440088161091608174329024836644263908375/3200608605869260877638343770033.

The stacked exact condition matrix has rank `2` on ambient dimension `4`,
hence kernel dimension and certified `h0(D)` are both `2`.

The exact kernel is the pencil

    d(V) = -35*V^2 + (-2460465/26923*j + 1445850/26923)*V - 39505725/215384*j + 37022265/215384
    a(V)+m,  a(V) = (277610270613495594051820346318069949105/171662846892070449329932340524*j + 2031453203570192491242235076958635251845/171662846892070449329932340524)*V + 3575994339198526543019597053284230977575/343325693784140898659864681048*j - 3469004253354616873638551659297641079935/343325693784140898659864681048

and therefore

    T = (a(V)+m)/d(V).

Modulo 73 this is the pinned 7774 parameter after the certified parent gauge.

## Child

Binary-quartic invariants from this same resolved pencil give an exact
characteristic-zero Jacobian with fibres

    I6 + 2 I5 + 2 I2 + 4 I1.

Its reduction is the pinned GF(73) q6_7774 equation (constant Weierstrass
scale marker `16`).

Next: propagate the exact child through q4_1938.
