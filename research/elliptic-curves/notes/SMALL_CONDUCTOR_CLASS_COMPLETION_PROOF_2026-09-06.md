# MW16 at 3/17: a sixteen-generator certificate under GRH

This is the proof behind the [current status and replay index](SMALL_CONDUCTOR_CLASS_TARGET_2026-09-06.md).
`MATH_STATUS.json` remains the mathematical authority. The calculation concerns
only `a1-fibration-05` at `3/17`, curve `new-20260905-36`.

## Statement and assumption

For the cubic field

```text
K = Q[z]/(z^3 + z^2 - 2919231625641258502793755607986240*z
          + 45440201616242830029801770634418828098464819545088),
disc(K) = 128900477062442043600727490102612931938219670661531295245188752203875468,
```

the certificate proves `dim_F2 Cl(K)/2 = 16`, conditional on GRH for the
nontrivial quadratic characters of its ordinary ideal class group. It actually
uses GRH only for such characters trivial on the sixteen specified anchors.
The known independent points and the previously proved local and parity bounds
then imply **rank E(Q) = dim_F2 Sel_2(E/Q) = 22 under that assumption**.
An unconditional upper bound is not asserted.

## Exact algebraic input

The lower-bound and character certificates construct sixteen independent
everywhere-unramified quadratic extensions. Their character values on sixteen
listed prime ideals form an invertible matrix over F2. These sixteen ideal
classes therefore span an unconditional sixteen-dimensional subspace
`H` of `Cl(K)/2`.

The retained rational-prime, norm and point-derived principal parity relations
are exact relations in `Cl(K)/2`. Outside coordinates are retained until they
cancel. After the final residual wave the formal quotient on the earlier
generating base still has dimension **18**. This formal presentation is not
assumed to be complete.

For each prime ideal of norm below `T=50000`, reduce its coordinate using the
audited principal relations. A prime is called known only when its resulting
coordinate lies entirely in the anchor span. An uncancelled outside coordinate
means unknown. A second full-coordinate row-space calculation verifies every
claimed representation and independently reproduces the known/unknown
classification. The inherited maximal-order audit includes every prime ideal
above every rational prime up to 400000, so the prime-power sum is complete.

## Quadratic-character exclusion

Put `L=log T`, and let `F(x)=max(1-x/L,0)` on the nonnegative real axis.
Its even extension has Fourier transform
`L*(sin(tL/2)/(tL/2))^2`, including its limit `L` at zero, hence nonnegative.
Write

```text
S_chi = sum_P sum_(m>=1) log(NP) F(m log(NP)) chi(P)^m / (NP)^(m/2),
I = integral_0^infinity (1-F(x))/(2 sinh(x/2)) dx,
J = integral_0^infinity F(x)/(2 cosh(x/2)) dx,
C = log(disc(K)) - 3*(gamma+log(8*pi)) + 3*I - 3*J.
```

For a nontrivial quadratic ordinary ideal-class character, the explicit formula
has no pole contribution and gives

```text
sum_rho Phi(rho) = C - 2*S_chi.
```

Under the stated GRH, every term on the left is nonnegative. Thus
`2*S_chi-C <= 0`. This is the character formula in
[Belabas–Diaz y Diaz–Friedman (2008), equation (3), pp. 1188–1189](https://doi.org/10.1090/S0025-5718-07-02003-0),
also available from the [author's institutional repository](https://repositorio.uchile.cl/handle/2250/154648).
The ordinary class character is unramified at finite and real places, so its
archimedean terms agree with the trivial character's terms. No assumption about
the zeros of the trivial character's L-function is needed for this argument.

If H were proper in `Cl(K)/2`, finite F2-linear duality would give a nontrivial
quadratic character chi trivial on H. It equals +1 at every known prime. At an
unknown prime its even powers are +1 and its odd powers are at least -1.
Since all weights in S are nonnegative, writing U for the weighted sum over
unknown prime ideals and odd powers gives

```text
S_chi >= S_1 - 2*U,
2*S_chi - C >= (2*S_1-C) - 4*U.
```

The certificate encloses the right side in a strictly positive interval.
This contradicts the explicit formula and proves **H = Cl(K)/2**. In particular
the sixteen independent anchors generate, and `g=16` under GRH. The two missing
formal norm relations are not needed, and are not manufactured by this proof.

## Interval arithmetic and independent check

The first calculation uses 256-bit MPFI intervals and 4096 terms of the
exponential series for I and J. For `lambda=k+1/2`, the omitted absolute tail is
bounded by

```text
exp(-(N+1/2)*L)/((N+1/2)^2*(1-exp(-L))).
```

This bound encloses the positive and alternating series tails. The earlier
triangle checker provides `2*S_1-C`; the new checker subtracts exactly the
worst-case odd-power penalty `4*U`. At T=50000 the corrected margin is about
**17.161534309987**, safely positive.

A separate 160-bit calculation directly accumulates the worst signed sum and
uses the weaker archimedean inequality

```text
C < log(disc(K)) - 3*(gamma+log(8*pi)) - 3*pi/2
    + (3*pi^2/2 + 12*Catalan)/log(T).
```

It follows by replacing `1-F(x)` by `x/L` in the positive archimedean kernels;
it uses no truncated exponential series. Its lower margin must also be positive.
Both calculations enumerate every prime power below T, including powers of
ramified and non-degree-one prime ideals.

## Curve consequence and reuse

The independently replayed arithmetic and parity implication is
`22 <= rank E(Q) <= dim Sel_2(E/Q) <= g+7`, with even 2-Selmer dimension.
Putting `g=16` gives an upper bound 23, then 22 by parity. This proves exact
rank 22 under the class-character GRH assumption. It does not produce a
twenty-third point or establish unconditional algebraic rank parity.

The reusable step is the **sufficient test for a proposed class-group span**:
verify principal relations, classify prime ideals whose classes are in that
span, assign worst-case signs to the rest, and test the explicit-formula margin.
This can stop relation collection before a full formal presentation is found.
The function `conservative_character_margin` accepts a general field's degree,
signature, discriminant and certified prime data. Each application still needs
its own exact arithmetic, span certificates and elliptic Selmer bounds. A
nonpositive margin is inconclusive; no universal fast or unconditional
rank-proving algorithm is claimed.

The [completion certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_class_completion_v1.json)
records the exact intervals, all known representations, all unknown prime
ideals and source hashes. Run the bounded replay through the current index.
