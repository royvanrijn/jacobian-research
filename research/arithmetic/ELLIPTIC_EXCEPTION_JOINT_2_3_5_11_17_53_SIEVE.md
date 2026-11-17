# A correlated fifty-three-adic refinement of the elliptic sieve

Adding the test at fifty-three to `2,3,5,11,17` reduces the exact
finite-place coefficient to
`2081912199403/9893175318390`, approximately `0.2104392303`.
It also reduces the certified quadratic class's coefficient from
`1/8` to `13361/114520`. The new residue condition is correlated
with the earlier tests at eleven and seventeen; those correlations
are retained explicitly.

## 1. The finite image and the needed saturation

Retain the curve and points of
[the five-prime sieve](ELLIPTIC_EXCEPTION_JOINT_2_3_5_11_17_SIEVE.md):
\[
 E:y^2=x^3-10656x^2+27091730x,\qquad
 X=x-3552,\quad t=-X/y,\quad F(x)=128x-25865,
\]
\[
 P=(373248/49,59459616/343),\quad T=(0,0),\quad R=P+T,
 \qquad S=(92179201/6400,553221611649/512000).
\]
Write `Gamma=E(Q)`, `Gamma_0=<R,S,T>`, and `q_0=2709173`.
No displayed free basis of `Gamma` is assumed.

**Lemma J531 (finite image, five-saturation and new formal image).**
The index `[Gamma:Gamma_0]` is prime to `2,3,5`. At fifty-three,
\[
 E(\mathbb F_{53})\simeq C_2^2\times C_{15},\qquad
 \bar R=(22,11),\quad\bar S=(44,8),\quad\bar T=(0,0).
                                                               \tag{1.1}
\]
Put `A=(30,0)`, `B_3=(13,45)` and `B_5=(36,27)`. One has
\[
 [15]\bar R=A,\quad[15]\bar S=[15]\bar T=\bar T,
\]
\[
 [20]\bar R=B_3,\quad[20]\bar S=[20]\bar T=O,
 \qquad
 [12]\bar R=[12]\bar S=B_5,\quad[12]\bar T=O.       \tag{1.2}
\]
Here `A,T` generate the two-primary group, and `B_3,B_5` have
orders three and five. The actual formal image is `E_1(53)`.

**Proof.** Prime-to-two index is ETF1; prime-to-three index is the
finite saturation certificate ETC3a. Its second point is
`Q=(3442688/529,304670016/12167)` and satisfies `S=-P-Q`, so
`<R,S,T>=<P,Q,T>`. For five use the two additive
characters on `Gamma` from the five-prime proof: the normalized
formal digit at five and the order-five reduction coordinate at
eleven. On `R,S` their matrix is
\[
       \begin{pmatrix}2&1\\1&1\end{pmatrix}\pmod5,
                    \qquad\det=1.                \tag{1.3}
\]
Both characters kill `5 Gamma` and the order-two torsion.
Consequently `R,S` are independent in `Gamma/5 Gamma`, which
has dimension two over `F_5` because the rank is two and there is
no five-torsion. They form a basis of that quotient. Equivalently
the determinant of their free-lattice coefficient matrix is a
five-adic unit, so `5` does not divide `[Gamma:Gamma_0]`.
This is a finite saturation argument, not an assertion of a full
Mordell--Weil basis.

Complete finite-field enumeration gives sixty points at fifty-three;
addition verifies (1.1)--(1.2). In particular `A+T=(26,0)` is the
third nonzero two-torsion point. The points in (1.2) generate a
subgroup of order `4*3*5=60`, proving the full finite image.
Exact rational multiplication also gives
\[
       v_{53}(t([60]P))=1,\qquad
       t([60]P)=2491=53\cdot47\pmod{53^2}.         \tag{1.4}
\]
The formal filtration therefore gives the actual image `E_1(53)`.
\(\square\)

In the primary coefficient coordinates `mR+nS+cT`, `c in {0,1}`,
the new two-, three- and five-primary zero conditions are thus
\[
       mA+(n+c)T=O,\qquad m=0\pmod3,
                         \qquad m+n=0\pmod5.     \tag{1.5}
\]
Prime-to-`2,3,5` index makes these relations valid for the actual
primary completion of `Gamma`, not just for the known subgroup.
The new fifty-three-primary formal coordinate is independent of all
earlier tests: their local finite quotients use only the primes
`2,3,5`, and their formal coordinates use `2,3,5,11,17`.
Chinese remainder approximation separates the new `53` factor.

The local group inputs are Milne, *Elliptic Curves* (2006),
Chapter II, Theorem 4.1 and Lemma 5.9, printed pp.62--63 and 67:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The existing two-/three-saturation inputs are ETF1 and ETC3a in
[the two/five-adic classification](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md)
and [the complete three-adic classification](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).

## 2. Pole probability after the earlier five tests

Write `S_0={2,3,5,11,17}`. All conditional probabilities below
are in the actual rational closure and are conditioned on the
specified global class and solubility at every prime in `S_0`.

**Lemma J532 (the correlated pole fractions).** The conditional
probability of reduction to `O` at fifty-three is
\[
 \rho_a=
 \begin{cases}
 1/(615\cdot533),&a=1,\\
 0,&a=2\text{ or }5q_0,\\
 1/15,&a=10q_0.
 \end{cases}                                                   \tag{2.1}
\]
There are no finite zero residues for `F` at fifty-three. Given
reduction to `O`, the soluble formal fraction is `1/2863`.

**Proof of the two-primary conditions.** The old coordinate and
class conventions of J2171 and J235172 give
\[
 \alpha=(5q_0)^{m\bmod2}(10q_0)^c,
 \qquad \bar P_3=cA_3+(n\bmod2)B_3',
\]
where `A_3,B_3'` are the old two-torsion residue basis at three.
The normalized formal coordinate at two has parity `m+n`.
For class one, `m` is even and `c=0`; the old two-adic soluble
test forces `n` even, since the normalized formal coordinate is
divisible by four. For the certified class `10q_0`, one has
`m` even and `c=1`; the old three-adic test forces `n` odd,
landing in its always-soluble `X=2` disc. Thus the new two-primary
condition in (1.5) holds in both cases. In the other two classes
`m` is odd, so it cannot hold: `A` and `T` are independent.
This proves their zero pole probabilities without presuming
independence of the finite two-primary coordinates.

**The class-one three-primary factor.** On the formal group at
three, choose the normalized additive coordinate `z_3` used in
J2171. Its coefficient on `R` is divisible by three and that on
`S` is a unit. The old order-nine coordinate at seventeen is
`w_3=m+n mod9`, uniform and independent of the entire `z_3`.
Class one lies in the `O` disc at three; solubility there forces
formal depth at least three. The actual formal closure starts
at depth two, so `z_3=0 mod3`, and hence `n=0 mod3`.
The new order-three zero condition `m=0 mod3` is therefore
equivalent to `w_3=0 mod3`.

In class one the two-primary coordinate at seventeen is zero.
Exactly the three values `w_3=0,3,6` are branch residues, each
with old soluble mass `1/307`; the other six values are fully
soluble. Consequently
\[
 \mathbb P(m=0\bmod3\mid\text{old tests})
            =\frac{(1/3)(1/307)}{2/3+(1/3)(1/307)}
            =\frac1{615}.                         \tag{2.2}
\]
The old three-adic unit-cube condition does not change this ratio,
because `w_3` is independent of the entire formal coordinate.

**The class-one five-primary factor.** The new zero condition
`m+n=0 mod5` is precisely the old order-five pole coordinate at
eleven. By J23511171 this coordinate is uniform and independent
of the entire formal coordinate at five. The finite part at five
depends on the two-primary data, which we may fix first. Whether
that state puts the point in the five-adic formal disc or outside
it, the test at five does not change the uniform eleven residue
coordinate. The old eleven test retains mass `1/133` at its zero
coordinate and mass one at its other four values. Thus
\[
 \mathbb P(m+n=0\bmod5\mid\text{old tests})
                =\frac{(1/5)(1/133)}{4/5+(1/5)(1/133)}
                =\frac1{533}.                     \tag{2.3}
\]
After fixing the old two-primary state, the old three-/seventeen-
primary data and the five-/eleven-primary data are product images.
The conditional ratios (2.2) and (2.3) are constant in that state.
They therefore remain independent after all old tests and after
integrating over the allowed two-primary states. Their product
proves the class-one row of (2.1).

**The certified class.** Here the old conditions reduce to the
two-primary restrictions `m` even, `n` odd, `c=1` just described.
The old two-adic translated component is always soluble, the
three-adic `X=2` disc is always soluble in the actual closure,
the five-adic formal disc is absent, and the tests at eleven and
seventeen are automatic. Thus the three- and five-primary
coordinates remain unrestricted and independent. The zero
conditions `m=0 mod3`, `m+n=0 mod5` have respective probabilities
`1/3` and `1/5`. This gives `1/15`.

Finally `(51730/53)=-1`, so `F` has no finite zero residues.
By ELV1 only the pole can obstruct solubility. Its depth-one
formal image from (1.4) has soluble mass
\[
       \sum_{k\ge1}(1-53^{-1})53^{1-3k}
                    =\frac1{53^2+53+1}=\frac1{2863}. \tag{2.4}
\]
The new formal coordinate is independent of every old condition
and of reduction at fifty-three. This completes the proof.
\(\square\)

## 3. Exact six-prime coefficients and Hasse bounds

**Theorem J533 (the six-prime count).** Let
`S=S_0 union {53}`. The admissible full field targets of class
`a` passing every test in `S` satisfy
\[
 M_{S;a}(B)=c_a C_E\log B+O_E(\sqrt{\log B}),
\]
with whole-slice normalized coefficients

| Class | Coefficient `c_a` |
| --- | --- |
| `1` | `1251298964/1648862553065` |
| `2` | `205/2456` |
| `5q_0` | `533/55860` |
| `10q_0` | `13361/114520` |

In particular
\[
 \boxed{M_S(B)=\frac{2081912199403}{9893175318390}C_E\log B
                                +O_E(\sqrt{\log B}),}           \tag{3.1}
\]
\[
 \boxed{A_E(B)\le\frac{2081912199403}{9893175318390}C_E\log B
                                +O_E(\sqrt{\log B}),}           \tag{3.2}
\]
and for the certified quadratic field
\[
 \boxed{A_{K_0}(B)\le\frac{13361}{114520}C_E\log B
                                +O_E(\sqrt{\log B}).}           \tag{3.3}
\]

**Proof.** Multiply the coefficients of J23511173 by
`1-rho_a(1-1/2863)` from J532. The class-one retention is
`312824741/312825695`; the certified-class retention is
`13361/14315`. The other two factors are one. Exact rational
arithmetic gives the table and (3.1). ESQ3 supplies the fixed-set
square-root logarithmic error. Every all-places Hasse target
passes these six tests, giving (3.2)--(3.3). \(\square\)

The exact coefficients count a fixed finite local prefix. Neither
the smaller certified-class bound nor five-saturation proves
Hasse infinitude, determines the intrinsic all-primes Haar limit
or controls the arithmetic tail in EFA3.

## 4. Exact controls and scope

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_elliptic_exception_joint_2_3_5_11_17_53.py --write-receipt
```

The checker enumerates the complete reduction group at fifty-three,
checks the primary projections and group generation, replays the
five-prime controls, verifies the saturation character matrices,
and computes (1.4) using exact rational arithmetic. Its finite
weighted tables check (2.2)--(2.3), including both possible finite
five-adic states, and all six-prime rational coefficients.
The receipt pins the source and exact inputs. These controls
support the written actual-image and conditional-probability
arguments; they do not replace their unrestricted formal-group
or counting inputs. No full displayed basis, formal verification,
external review or all-primes asymptotic is claimed.
