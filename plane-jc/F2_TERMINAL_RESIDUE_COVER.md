# The degree-six terminal residue cover in F2 `(75,125)`

## Result and claim boundary

The uniquely normalized F2 terminal block does more than fix five endpoint
coefficients.  It determines an exact toroidal source-to-target boundary row.
After extracting the target ray `(5,2)`, the terminal source divisor has
transverse index one and residue map

\[
h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}. \tag{1}
\]

This map has degree six, three branch values, branch passport

\[
(5,1),\qquad(3,3),\qquad(3,1,1,1), \tag{2}
\]

and monodromy group `A_6`.  Its branch cycles satisfy the complete global
meridian relation.

This supplies genuine target-side data that were absent from the earlier
contact-only handoff.  It does **not** yet exclude `(75,125)`: the extracted
target component lies in the toroidal target boundary, and the remaining
source boundary, spectator orbits, affine sheets, and target-transfer gluing
are not yet classified.

The exact checker is
[`cas/verify_f2_terminal_residue_cover.py`](cas/verify_f2_terminal_residue_cover.py).

## 1. Terminal block

Put

\[
s=X^{17}y^5.
\]

The forced terminal type-I block is

\[
P=X^4y(1+s),
\]

\[
Q=-X\left(1+3s+\frac95s^2\right). \tag{3}
\]

Direct differentiation gives

\[
[P,Q]_{X,y}=X^4. \tag{4}
\]

On the torus `Xy != 0`, form the target character

\[
r=\frac{P^5}{(-Q)^3}.
\]

Equation (3) gives

\[
r=h(s)
=\frac{s(1+s)^5}
       {(1+3s+\frac95s^2)^3}
=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}. \tag{5}
\]

Conversely, once a root `s` of `h(s)=r` is chosen,

\[
X=\frac{-Q}{1+3s+\frac95s^2},
\]

\[
y=\frac{P}{X^4(1+s)}. \tag{6}
\]

Substituting (6) gives `X^17 y^5=s`.  Therefore

\[
[k(X,y):k(P,Q)]=\deg h=6 \tag{7}
\]

for the terminal torus block.  The numerator and denominator in (1) are
coprime degree-six polynomials.

## 2. Source and target toroidal rays

In the Laurent coordinates

\[
t=Xy,\qquad z=y^{-1},
\]

the terminal supports are

\[
P:t^4z^3+t^{21}z^{15},
\]

\[
-Q:tz+3t^{18}z^{13}+\frac95t^{35}z^{25}. \tag{8}
\]

Both support segments have direction `(17,12)`.  The primitive normal

\[
\nu=(12,-17) \tag{9}
\]

has constant values

\[
\nu(P)=-3,
\qquad
\nu(Q)=-5. \tag{10}
\]

Thus the source toric divisor has pole orders `(3,5)` for `(P,Q)`.
Near the `Q`-dominant target-infinity point put

\[
a=(-Q)^{-1},
\qquad
b=P/(-Q). \tag{11}
\]

Their pullback orders along the source divisor are

\[
(\nu(a),\nu(b))=(5,2). \tag{12}
\]

Hence the required target extraction is the primitive ray `(5,2)`.  Its
regular neighboring rays may be chosen as `(3,1)` and `(2,1)`, since

\[
\det\begin{pmatrix}3&5\\1&2\end{pmatrix}=1,
\qquad
\det\begin{pmatrix}5&2\\2&1\end{pmatrix}=1.
\]

On the chart adjacent to `(3,1)`, a transverse uniformizer and residue
coordinate are

\[
\pi=\frac{b^3}{a},
\qquad
\eta=\frac{a^2}{b^5}. \tag{13}
\]

Their source orders are

\[
\nu(\pi)=3\cdot2-5=1,
\qquad
\nu(\eta)=2\cdot5-5\cdot2=0. \tag{14}
\]

Therefore the source divisor maps to the extracted target divisor with
transverse index

\[
\boxed{e=1}. \tag{15}
\]

Moreover

\[
\eta^{-1}=\frac{b^5}{a^2}
          =\frac{P^5}{(-Q)^3}
          =h(s). \tag{16}
\]

The residue degree is consequently

\[
\boxed{f=6}. \tag{17}
\]

This is a certified toroidal boundary row `(e,f)=(1,6)`, not an unsupported
promotion of an edge-root contact multiplicity.

## 3. Branch passport

The derivative of (1) collapses to

\[
h'(s)=
\frac{625(s+1)^4}{(9s^2+15s+5)^4}. \tag{18}
\]

The three branch fibers are:

- over `0`: `s=0` with index one and `s=-1` with index five;
- over `infinity`: the two roots of `9s^2+15s+5`, each with index three;
- over `125/729`: `s=infinity` with index three and the three simple roots of
  `135s^3+405s^2+396s+125`.

The quadratic discriminant is `45`; the cubic discriminant is `-98415`, so
all displayed finite roots are distinct.  The total different is

\[
4+2+2+2=10=2\cdot6-2, \tag{19}
\]

which verifies Riemann--Hurwitz.

## 4. Monodromy and meridian relation

Choose a branch cycle over `0` of type `(5,1)`.  An exhaustive permutation
calculation in `S_6` finds five compatible cycles of type `(3,3)` for the
second branch value such that the inverse product has type `(3,1,1,1)`.
They form one centralizer orbit.  Every resulting transitive triple generates
all `360` even permutations:

\[
\boxed{G=A_6}. \tag{20}
\]

For every triple

\[
\sigma_0\sigma_\infty\sigma_{125/729}=1, \tag{21}
\]

so (21) is the actual global meridian relation of the terminal residue cover,
not an abstract endpoint matching.

## 5. Consequences for the F2 programme

Combined with
[`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md), the terminal
chain now supplies:

1. one source Kummer orbit with a completely known principal Newton block;
2. the source terminal ray `(12,-17)`;
3. the target extraction ray `(5,2)`;
4. transverse index `e=1`;
5. residue degree `f=6`;
6. the complete residue ramification passport and `A_6` monodromy;
7. an exact global meridian factorization.

The next gap is global rather than local.  One must attach this row to the
original `A^2` completion, classify the simple `R` spectator orbits, and in
the double-root row determine whether the two identical `A_6` packets land on
the same or distinct target boundary components.  Only after that gluing can
the class-group, unit, canonical, and finite-normalization ledgers be run.

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
```

Expected final marker:

```text
F2_TERMINAL_RESIDUE_COVER_PASS
```
