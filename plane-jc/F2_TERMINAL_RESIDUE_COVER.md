# The degree-six terminal residue cover in F2 `(75,125)`

## Result and claim boundary

The uniform extension to every F2 parameter is proved in
[`F2_MODIFIED_LAURENT_FAMILY.md`](F2_MODIFIED_LAURENT_FAMILY.md): the terminal
residue has degree `2r`, passport
`(2r-1,1)|(r,r)|(3,1^(2r-3))`, and geometric monodromy `A_(2r)`.  This note
retains the fuller toroidal and global-ledger analysis of its `r=3`
specialization.

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

and geometric monodromy group `A_6`.  Its branch cycles satisfy the complete
global meridian relation.  The natural degree-six `A_6` action is
four-transitive and primitive, so the residue cover is geometrically
indecomposable; its target-fixed deck group is trivial.  Over
`Q(r)` the arithmetic monodromy is instead `S_6`, with quadratic constant
field `Q(sqrt(5))` in the Galois closure.  After rescaling it is a Belyi map
whose regular geometric closure has signature `(5,3,3)` and genus `25`.
The two target toric nodes have three preimages in the source-divisor
interior, giving three forced attachment points in any global boundary
realization.

This supplies genuine target-side data that were absent from the earlier
contact-only handoff.  It does **not** yet exclude `(75,125)`: the extracted
target component lies in the toroidal target boundary, and the remaining
source boundary, spectator orbits, affine sheets, and target-transfer gluing
are not yet classified.

None of the derived consequences below proves that these global attachments
exist or are impossible.  In particular they do not exclude `(75,125)`,
improve the conditional degree frontier `125`, make the conditional 14- and
22-equation modified systems exhaustive, or prove `JC(2)`.

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

In characteristic zero the transverse contribution to the different is
`e-1`, hence it is zero on this row.  All nontrivial ramification recorded
here is residue ramification, not ramification in the normal direction.

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
all displayed finite roots are distinct.  The four nonzero residue-different
coefficients are

\[
(4,2,2,2), \tag{19}
\]

and the total different is

\[
4+2+2+2=10=2\cdot6-2, \tag{20}
\]

which verifies Riemann--Hurwitz.

After the target rescaling

\[
\beta(s)=\frac{729}{125}h(s),
\]

the branch values are `0`, `infinity`, and `1`.  Thus `beta` is a degree-six
Belyi map over `Q` with the same passport.

There is also a forced endpoint-incidence pattern.  The toric endpoints of
the source divisor are `s=0,infinity`, while the target toric endpoints are
`h=0,infinity`.  Directly,

\[
h(0)=0,
\qquad
h(\infty)=\frac{125}{729}. \tag{21}
\]

The other point over the first target endpoint is the interior point `s=-1`
with index five.  The two points over the second target endpoint are the
distinct roots of `9s^2+15s+5`, both interior and both of index three.
Consequently the two target nodes have exactly three preimages in the
interior of the terminal source divisor.  In any resolved morphism of
boundary pairs, each is a forced attachment point for another source-boundary
branch.  Conversely, the source endpoint `s=infinity` maps with tangential
index three to the smooth target point `125/729`.

## 4. Monodromy and meridian relation

Choose a branch cycle over `0` of type `(5,1)`.  An exhaustive permutation
calculation in `S_6` finds five compatible cycles of type `(3,3)` for the
second branch value such that the inverse product has type `(3,1,1,1)`.
They form one centralizer orbit.  Every resulting transitive triple generates
all `360` even permutations:

\[
\boxed{G_{\mathrm{geom}}=A_6}. \tag{22}
\]

The natural action of this group on the six sheets has no nontrivial block
of size two or three and is transitive on ordered tuples of up to four
distinct sheets.  It is therefore primitive and four-transitive.
Equivalently, the degree-six residue cover has no nontrivial intermediate
cover and cannot factor as covers of degrees two and three (in either order).
The centralizer of this action in `S_6` is the identity, so

\[
\operatorname{Aut}_{\mathbf P^1_r}(\mathbf P^1_s)=1. \tag{23}
\]

Thus two copies of this residue cover, if global gluing places them over the
same target component, have at most one target-fixed identification.

For every triple

\[
\sigma_0\sigma_\infty\sigma_{125/729}=1, \tag{24}
\]

so (24) is the actual global meridian relation of the terminal residue cover,
not an abstract endpoint matching.

There is a useful arithmetic/geometric distinction.  Over `Q(r)` the cover
is defined by

\[
\Phi_r(s)=125s(s+1)^5-r(9s^2+15s+5)^3, \tag{25}
\]

and exact elimination gives

\[
\operatorname{disc}_s(\Phi_r)
=5^{17}r^4(729r-125)^2. \tag{26}
\]

Its squareclass is `5`.  The geometric group in (22) is therefore enlarged
by an odd permutation over `Q(r)`, and the only possible overgroup is

\[
G_{\mathrm{arith},\mathbf Q(r)}=S_6. \tag{27}
\]

The fixed quadratic field in the Galois closure is the constant extension
`Q(sqrt(5))`; after adjoining `sqrt(5)`, the arithmetic and geometric groups
both equal `A_6`.

The geometric Galois closure is the regular `A_6` cover with inertia orders
`(5,3,3)`.  Its Riemann--Hurwitz calculation is

\[
2g-2
=360\left(-2+\left(1-\frac15\right)
              +\left(1-\frac13\right)
              +\left(1-\frac13\right)\right)
=48,
\]

so its genus is

\[
\boxed{g=25}. \tag{28}
\]

## 5. Processed implications for the F2 programme

Combined with
[`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md), the terminal
chain now supplies:

1. one source Kummer orbit with a completely known principal Newton block;
2. the source terminal ray `(12,-17)`;
3. the target extraction ray `(5,2)`;
4. transverse index `e=1`;
5. residue degree `f=6`;
6. the complete residue ramification passport and geometric `A_6` monodromy;
7. an exact global meridian factorization;
8. a four-transitive, primitive, geometrically indecomposable residue packet
   with trivial target-fixed deck group;
9. zero transverse different contribution, with residue-different packet
   `(4,2,2,2)` and total degree `10`;
10. a three-point Belyi normalization whose regular `A_6` closure has triangle
    signature `(5,3,3)` and genus `25`;
11. arithmetic `S_6` over `Q(r)`, geometric `A_6`, and quadratic constant
    field `Q(sqrt(5))`;
12. three forced interior attachment points over the target toric nodes, while
    the source endpoint `s=infinity` maps to the smooth third branch value;
13. a parameter-free residue map: equation (1) contains no coefficient of
    the quadratic cofactor `R`, so the same packet occurs uniformly on every
    admissible principal-chain stratum.

Together with the Kummer-orbit theorem this leaves the following finite case
split for the first global ledger:

- if `R` is squarefree and nonzero at zero, there is one principal `A_6`
  packet; its two simple-root Kummer orbits are spectators;
- if `R` has a nonzero double root, there are two identical principal `A_6`
  packets, and they either land on the same target boundary component or on
  two distinct target components.

### 5.1 Global degree ledger

Let

\[
d=[k(x,y):k(P,Q)]
\]

for any global realization of this F2 row.  The fundamental equality over
the target valuation extracted by `(5,2)` is

\[
d=\sum_{D\mid T}e(D/T)f(D/T). \tag{29}
\]

The certified row contributes `1*6`, so every such realization satisfies

\[
\boxed{d\ge6}. \tag{30}
\]

On the nonzero double-root stratum, the two principal chains give two
distinct source valuations.  If both lie over this same target valuation,
their contributions add and

\[
\boxed{d\ge12}. \tag{31}
\]

If they land on distinct target components, the two fundamental equalities
are separate and only the floor `d>=6` follows.  This distinction is exactly
why the same-target versus distinct-target split cannot be collapsed.

There is no affine-sheet `+1` in (30) or (31).  Equations (10) show that both
`P` and `Q` have poles on `T`; the valuation is centered on the target
boundary at infinity.  The affine-companion theorem for a nonproperness
curve applies to target curves inside `A^2`, not to this boundary valuation.
In particular, this calculation does **not** prove the tempting bounds
`d>=7` or `d>=13`.

### 5.2 Forced attachment and different placement

The endpoint calculation after (21) already fixes part of the global graph
before any spectator classification.  The point `s=-1` over `h=0` and the two
denominator roots over `h=infinity` lie in the interior of the source
divisor.  Because both target values are boundary nodes, a resolved boundary
morphism must carry another source-boundary branch through each of these
three distinct points.  Their residue-different contributions are `4,2,2`,
of total degree eight.

The remaining contribution `2` occurs at the source toric endpoint
`s=infinity`, which maps to the smooth target value `125/729` with index
three.  The three simple interior points in that fiber contribute no
different.  Thus the global gluing ledger must distinguish

\[
\underbrace{(4,2,2)}_{\text{target-node attachments}}
\quad\text{from}\quad
\underbrace{(2)}_{\text{source endpoint over a smooth target point}}.
\]

### 5.3 Purity and global monodromy

Because `T` is centered at infinity and `e=1`, this row supplies no
codimension-one branch component of the canonical finite normalization over
`A^2`.  Purity and the absence of nontrivial connected finite etale covers of
`A^2` therefore force any global
noninvertible realization to contain some additional missing-boundary row
with `e>1` over an affine nonproperness curve.  That separate row has its own
positive affine companion, but its fundamental equality is over a different
target valuation and cannot be added to (29).

The surface different coefficient on the terminal row is zero.  The residue
map nevertheless has the four-point different packet `(4,2,2,2)`.  Any
global boundary graph must place those four marked ramification points, with
total residue different ten, among its nodes and smooth marked points; it
may not reuse ten as a transverse surface-different coefficient.

Finally, in the geometric Galois closure of a global realization, the
decomposition group at `T` has the residue `A_6` Galois group as a quotient.
Consequently the global geometric monodromy is nonsolvable and its order is
divisible by `360`; equivalently, `A_6` is a nonabelian simple composition
factor.  If equality holds in (30), the global monodromy is a
transitive subgroup of `S_6` of order divisible by `360`, hence is `A_6` or
`S_6`.  This is a global group-theoretic restriction, not an exclusion.

The five compatible branch-cycle triples are one centralizer orbit and all
generate `A_6`.  Thus the residue-side monodromy input is finite and already
classified; there is no smaller `2`-by-`3` factorization to analyze.  Since
the local deck group is trivial, a same-target identification of the two
double-root packets, if it exists, is unique.

The next gap is global rather than local.  One must attach this row to the
original `A^2` completion, classify the simple `R` spectator orbits, and in
the double-root row determine whether the two identical `A_6` packets land on
the same or distinct target boundary components.  Only after that gluing can
the class-group, unit, canonical, and finite-normalization ledgers be run.
None of (30)--(31) changes the coordinate-degree frontier `125`.

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
```

Expected final marker:

```text
F2_TERMINAL_RESIDUE_COVER_PASS
```
