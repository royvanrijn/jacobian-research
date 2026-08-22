# H3 q=8 current frontier

Status: 2026-08-22, after the source-side module-intersection and
normalization audit.

Detailed ledger: [`H3_Q8_MODULE_INTERSECTION_2026-08-22.md`](H3_Q8_MODULE_INTERSECTION_2026-08-22.md).

## Objective

Execute the exact second H3 neighbor

```text
H3 E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4
```

at equation level over characteristic zero.

The generic q8 fibre has degree 18 on the old fibration and exact basis

```text
B = 1,m,...,m^9, x,x*m,...,x*m^7,
m=(y-y(P1))/(x-x(P1)).
```

The literal source-nef vertical difference is

```text
(-11,0,2,3,4,6,5,5,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0).
```

The global `-11F` coefficient must be applied exactly once.  Do not inherit
and then double-count the q6 representative's `-F` factor through a ninth
power.

## Trusted local data

### Smooth h-collisions

At the four `O.(-P1)` fibres:

```text
p=y(P1)/x(P1),
q=(m-p)/h,
X=h^2*x,
```

and the exact local lattice is

```text
<1,q,...,q^9,X,Xq,...,Xq^7>.
```

`q` is local, not global.  The finite factor

```text
den(p)=h*u^2*d0,   deg(d0)=10
```

kills the fake eleven-dimensional global q-frame.  At both `43` and `59` the
associated 450-row transition matrix has rank 11 on the 11 candidates.

### E7 target

The resolved q8 E7 restriction degrees are `(0,1,0,0,2,0,1)`, and the exact
exceptional-cycle comparison is

```text
c8 = 9*c6 + (2,5,6,4,6,3,5).
```

This comparison is about the exceptional class; the common fibre
representative must be normalized separately.

After stripping the q6 helper fibre factor, the first reduced E7 saturation
step starts from `t^-1 B`.  The exceptional rows vanish at this layer and the
exact affine E7 identity

```text
y^2=x^3,  m=y/x,  x=m^2
```

gives ten rows of rank ten.  At both `43` and `59` the eight-dimensional
kernel is exactly

```text
x*m^j - m^(j+2),  j=0,...,7.
```

Hence use the adapted coordinate

```text
z = x-m^2.
```

**This is the current E7 frontier.**  Reduced E7 saturation step 2 is next.
The earlier absolute t-exponents obtained from the q6^9 helper normalization
are diagnostic only and are not the final q8 E7 lattice.

### E8 target

In actual chart order `(B1,B2,B3,B4,N3,N40,N4B,N4inf)`, the q8 exceptional
cycle is

```text
(-2,-4,-6,-10,-4,-7,-5,-8).
```

In integral II* coordinates

```text
u=1/t, X=u^4*x, Y=u^6*y,
Q=u^2*m,
```

the singular integral twist has exact complete ideal

```text
I=(u^2,X,Y),  quotient basis=(1,u).
```

Inside `C=(1,Q,...,Q^9,X,XQ,...,XQ^7)`, with `s=Y(P1)/X(P1)`, the full
preimage of this **singular-ring** ideal has basis

```text
u^2,
Q^b-s^b      (1<=b<=9),
X*Q^b        (0<=b<=7).
```

Do not call this the resolved E8 pushforward lattice yet.  Resolve/saturate it
against the actual II* blow-up valuations before global assembly.  The
previous helper-normalized determinant `342` is not a global q8 invariant.

## Retractions / do-not-repeat list

Do not revive these without new structural input:

- uniform `extra_e7_pole` enlargement as a route to the pencil;
- the local q/X collision frame as a global base frame;
- the child-side fractional finite/infinity shortcut whose 95-column test
  produced the 48-dimensional base-twist ladder;
- the inference `P0 smooth at IV* => identity component`;
- the provisional E7 determinant target `-98`;
- the claim that helper-normalized E7 saturation exponents were final;
- the claim that the singular integral E8 preimage plus helper `u^9` was the
  resolved saturated E8 lattice.

## Next exact work

1. Continue **reduced E7 saturation** from the mixed `(m^b,z*m^b)` frame,
   dividing by `t` and computing exact/two-prime obstruction kernels until a
   reverse-saturation step has full rank.
2. Perform the analogous **resolved E8 saturation** using the actual blow-up
   valuation atlas and the exact q8 E8 target cycle.
3. Express both reduced endpoint lattices in the common generic frame `B`,
   with no inherited q6 fibre powers.
4. Intersect them with the exact h-adic collision lattice using its local
   transition matrix.  Keep `q` local; the degree-10 `d0` factor is the guard
   against accidental globalization.
5. Apply the single global `-11F` twist and compute the global section space.
   The expected target is `h0(D)=2`.
6. From that two-dimensional kernel form the q8 pencil, eliminate it to a
   genus-one quartic/Jacobian, and verify the expected `D13/MW4` child before
   continuing the H3 lattice path.

## Final checksum, not a derivation tool

Only after both reduced resolved endpoint lattices are complete and expressed
in one common frame should determinant bookkeeping be used.  Since the q8
fibre class is isotropic on a K3,

```text
chi(O(D))=2,
rank(pi_*O(D))=18,
deg(pi_*O(D))=-16.
```

Use `-16` as a final vector-bundle checksum, not as a substitute for any local
saturation calculation.
