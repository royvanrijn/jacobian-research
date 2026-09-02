# E6 rank-sum-four ordered incidence and rational quotient — 2026-09-02

## Outcome

The systematic E6 incidence search reaches geometric rank four on a
one-dimensional ordered incidence curve.  The genus-zero **unordered**
quotient is rational over `QQ`, but the ordering needed to define the four
individual sections is a genus-one double cover.  Consequently the quotient
does not give a rank-four family over `QQ(k)`.

More precisely, the quotient has the smooth rational point

```text
(S,M)=(2,16)
```

and the exact normalization

```text
S=-(k^2+1)*(k^4+2*k^2+13)/(4*k*(k^2+3)),
M=-2*(k^2+1)*(k^2+3)/k^3.
```

The point corresponds to `k=-1`.  It lies on the diagonal boundary `v=w=1`,
but is smooth on the quotient and therefore suffices to identify its
normalization with `P1_QQ`.  The ordered lift is

```text
r^2=k^4+6*k^2+13.
```

This curve has genus one and is birational to the Cremona curve `52a2`, of
rank zero and rational torsion `Z/2`.  Its only rational points are the two
points above `k=infinity`, so there is no nondegenerate affine rational
incidence.  Rank four holds over the function field of this genus-one curve;
the descended surface over `QQ(k)` has arithmetic Mordell--Weil rank two.

Put

```text
z_v=v^2+2,       m_v=v*(v^2+3),
z_w=w^2+2,       m_w=w*(w^2+3),
a=2*(v^2+v*w+w^2+3)/(v+w),
c=-2*(v^2*w^2+2*v^2-v*w+2*w^2+6)/(v+w).
```

Then

```text
E: y^2=x^3+(a*u-3)*x+(u^2+c*u-2)
```

has the two sections

```text
P=(z_v,u+m_v),       Q=(z_w,u+m_w).
```

The generic rational-surface fibres are `IV*+4I1`, with the rational simple
fibre at `u=0`.  Its node is `x=-1`.  Branch the quadratic cover at this
`I1` and at the smooth value `u=M`:

```text
d(u)=u*(u-M),       u=M/(1-t^2).
```

For `x=v,w`, introduce a line slope `L=ell*u-x` through the negative of the
corresponding marked section.  The residual intersection discriminant is

```text
delta=L^4-6*(x^2+2)*L^2-8*(u+x*(x^2+3))*L
      -3*(x^2+2)^2-4*(a*u-3).
```

It has the required squareclass precisely when

```text
delta=u*(u-M)*(ell^2*u+q)^2,
q=-2*ell*x+M*ell^2/2.
```

After coefficient comparison, the complete node-collision linear-chord
conditions for the first section are

```text
(v+w)*ell^2*M*(ell*M-4*v)^2 = 32*(w^2+3),
ell*(3*ell^2*M^2-8*v*ell*M-16*v^2-48) = 32,
```

and the second pair is obtained by swapping `v,w` and using a second slope.
Thus the construction is an incidence curve, not a height box.

## Incidence components

Eliminating the two slopes and passing to the unordered coordinates

```text
S=v+w,       P=v*w
```

gives two irreducible plane components after eliminating `P`.  Their exact
geometric genera are

```text
0 and 2.
```

The genus-zero factor is

```text
-10368*S^3*M^5 +331776*S^4*M^2 +497664*S^3*M^3
+42624*S^2*M^4 +13728*S*M^5 +2197*M^6
-2359296*S^3*M -245760*S^2*M^2 -335872*S*M^3
-68352*M^4 +4194304*S^2 +3145728*S*M +589824*M^2 = 0.
```

The second factor, the eliminated `P` equations, and both slope equations are
stored in the generated certificate.  At the rational singular origin of the
genus-zero plane quotient, the first tangent is `8*S+3*M=0`; after that
blow-up the next tangent cone is

```text
16384*(256*X^2+9*M^2).
```

Those two branches above the origin are defined over `QQ(i)`, but this does
not obstruct rationality: `(S,M)=(2,16)` is a different smooth rational point.
Writing `P=v*w`, exact elimination and the normalization above give

```text
P=(-3*k^8-8*k^6+2*k^4+112*k^2-39)/(16*k^2*(k^2+3)),

S^2-4*P = h(k)^2*(k^4+6*k^2+13),
h(k)=(k^2-1)*(k^2+7)/(4*k*(k^2+3)).
```

Thus the quotient parameterization does not by itself order `v,w`.  On the
ordered cover `r^2=k^4+6*k^2+13`, all incidence data become explicit:

```text
v=(S+h*r)/2,                 w=(S-h*r)/2,

ell_v=k^2*(k^4+4*k^2+11+(1-k^2)*r)
      /(4*(k^2+1)*(k^2+3)),
ell_w=k^2*(k^4+4*k^2+11-(1-k^2)*r)
      /(4*(k^2+1)*(k^2+3)).
```

Literal substitution verifies both node-collision equations.  The quotient
surface itself descends to `QQ(k)`: substituting the symmetric functions gives

```text
a=-(k^8+7*k^6+29*k^4+37*k^2+22)
   /(2*k*(k^2+1)*(k^2+3)),

c=(9*k^12+62*k^10+243*k^8+612*k^6+1071*k^4+446*k^2+117)
   /(32*k^3*(k^2+1)*(k^2+3)).
```

The ordered curve is not a conic.  The birational maps to

```text
E_0: Y^2=X^3-64*X-192
```

are

```text
X=2*(r+k^2)+2,              Y=4*k*(r+k^2+3),
k=Y/(2*(X+4)),              r=(X-2)/2-k^2.
```

Sage's proved Mordell--Weil calculation identifies `E_0` as `52a2` and gives
`E_0(QQ)={O,(-4,0)}`.  Both points are at infinity in the ordered incidence
chart.  Hence a nondegenerate affine point needs an extension of degree at
least two.  Degree two is attained, for example at `k=2` over `QQ(sqrt(53))`:

```text
S=-185/56,                  M=-35/4,
v=(-185+33*sqrt(53))/112,   w=(-185-33*sqrt(53))/112,
ell_v=(43-3*sqrt(53))/35,   ell_w=(43+3*sqrt(53))/35.
```

This is minimal in degree, but it does not rationalize the family: geometric
genus one is unchanged by constant field extension.

The half-collision branches at the nodal fibre were also separated.  Their
small resultant component reduces to

```text
w=-v-r,       r^2=v^2+3,
```

and simultaneous half--half collision forces either duplicate marked
sections or the excluded divisor `v+w=0`.  The remaining half resultants are
higher-genus components.  This is why the D6 exact correspondence obstruction
was retained as the parallel control; it does not replace the successful E6
node--node component.

## Four independent directions on the ordered cover

Let `sqrt(d)` denote the chosen square root on the cover.  For either marked
section define

```text
x(R)=(L^2-z+(ell^2*u+q)*sqrt(d))/2,
y(R)=-(u+m)+L*(x(R)-z).
```

The discriminant identity proves that `R` is a section after base change.
Its conjugate satisfies `R+R'=P`, so

```text
T=R-R'=2*R-P
```

is anti-invariant.  Applying this to `P,Q` gives `T1,T2`.

At the complete good-reduction incidence point

```text
p=11,       (v,w,M,ell1,ell2)=(2,7,4,3,4),
```

clearing the component groups by multiplication by six and reading exact
pole growth gives

```text
Gram(T1,T2) = [22/3  4/3]
              [ 4/3 22/3],       determinant 52.
```

The invariant block is `A2*(2)`:

```text
Gram(P,Q) = [ 4/3 -2/3]
            [-2/3  4/3].
```

The quadratic-base-change deck characters make the two blocks orthogonal.
Hence the rank over the ordered incidence field is at least `2+2=4`.  The K3
fibres are

```text
2IV* + I2 + 6I1,
```

so the root rank is thirteen.  The marked family is nonconstant; therefore
generic geometric Picard rank is nineteen and Shioda--Tate makes the displayed
rank exactly four over the ordered field.

The parameter involution `r -> -r` swaps `P,Q` and swaps `T1,T2`.  Since these
four directions exhaust the geometric Mordell--Weil space, its fixed subspace
is exactly

```text
span(P+Q, T1+T2).
```

Thus the descended family over `QQ(k)` has exact arithmetic rank two, split
as one invariant and one anti-invariant direction for the K3 base-change deck
involution.  It is not the desired rank-four `QQ(k)` source.

## Saturation and NS determinant

The character basis is not saturated.  The construction itself supplies

```text
2*R1=P+T1,       2*R2=Q+T2.
```

In the basis `(P,Q,R1,R2)` the height matrix is

```text
[ 4/3 -2/3  2/3 -1/3]
[-2/3  4/3 -1/3  2/3]
[ 2/3 -1/3 13/6  1/6]
[-1/3  2/3  1/6 13/6],
```

with determinant `13/3`.  The root determinant of `2E6+A1` is `18`, hence

```text
abs(det NS)=18*(13/3)=78.
```

This is squarefree, so no further finite-index saturation is possible.  The
component groups allow only `2`- and `3`-primary generic torsion, while the
good control fibre at `p=11,t=2` has order `17`.  Specialization therefore
eliminates generic torsion.

The good-reduction section incidences also recover the full integral marking.
In the basis

```text
O,F,E6a_1,...,E6a_6,E6b_1,...,E6b_6,A1_1,P,Q,R1,R2
```

the section profiles, up to swapping the two `IV*` fibres and the `E6`
diagram involution, are

```text
P,Q:   (E6a_1,E6b_1,I2_identity),
R1,R2: (E6a_1,E6b_identity,A1_1).
```

The only nonzero intersections between distinct displayed sections are

```text
P.R2=Q.R1=1.
```

These two intersections occur at a smooth good-reduction fibre; the remaining
common Weierstrass zeros lie over the resolved `IV*` fibres.  The resulting
rank-19 integral Gram has Smith invariants

```text
1,...,1,78
```

and determinant `78`.  Splitting the marked `U=<F,O+F>` is unimodular and
gives the stored positive frame
[`data/lattice/e6_rank4_det78_frame.txt`](data/lattice/e6_rank4_det78_frame.txt).

## Rootless-MW17 feasibility

For a rootless fibration on this Picard-rank-19 K3, the MW frame would be even,
positive definite, rank seventeen, minimum at least four, and determinant
`78`.  Its required Hermite invariant is

```text
4/78^(1/17) = 3.0957102762...
```

while Blichfeldt gives

```text
gamma_17 <= 3.2821242499....
```

Thus the determinant screen **passes**.  Unlike the earlier determinant-24
E6 family, rootless MW17 is not ruled out by this numerical bound.

## Complete low-degree rootless search

<!-- status-consumer: EC-K3-E6-RANK4-ROOTLESS-Q2Q4-CENSUS 2351738f44774cfe -->

Write the integral split as `NS=U+M(-1)`, with `F=e` and `O=f-e`.  Proposition
C2.2 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
shows that the first zero-neutral old-degree-`q` layer is

```text
D=q*e+q*f-w,       w.M.w=2*q^2.
```

The complete Weyl-dominant census through `q=4` is:

| old degree | norm | dominant orbits | primitive | minimum signed child roots | rootless |
|---:|---:|---:|---:|---:|---:|
| 2 | 8 | 280 | 277 | 74 | 0 |
| 3 | 18 | 6,242 | 6,239 | 52 | 0 |
| 4 | 32 | 73,601 | 73,321 | 36 | 0 |

For each of the `79,837` primitive classes, the checker splits a primitive
marked `U` and asks PARI for the exact number of norm-two vectors in the
rank-17 child frame.  Every child remains rootful.  Dynkin-label enumeration
plus an exact rank-four closest-vector reconstruction proves completeness of
each shell modulo the full old `2E6+A1` Weyl group.  The physical
fixed-component degree distributions are stored as an independent chamber
audit.

This is an exact bounded negative result, not a global obstruction.  Degree
at least five and multi-step routes remain open.  A separate seeded
`p`-neighbour discovery run visited `63,521` frames and reached root rank ten
but not zero; because that run samples rather than exhausts the genus, it is
routing evidence only and is not part of the proof claim.

## Global frame-genus obstruction

<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL bd12c183aa886b15 -->

The later equation-free argument in Theorem H5 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
closes that bounded-search gap.  The determinant-78 rank-seven Nishiyama
auxiliary contains a primitive `A3+A2+A1`.  A complete residual-Weyl anchor
enumeration over all 23 rooted Niemeier lattices finds 1,591 primitive anchor
representatives, and every residual root system has rank at least 14.  A
rootless complement would force the final auxiliary vector to have positive
Dynkin label on every residual simple root, requiring projected norm at least
`14/2=7`; its exact available norm is only `13/4`.  This contradiction proves
that every rank-17 lattice in the target genus is rootful.

Thus degree at least five and multi-step routes are no longer open at the
`O(NS)`/J2 lattice level: no rootless MW17 frame exists for this saturated
Neron--Severi lattice.  This does not classify all frame isometry classes or
elliptic equations.

## Replay and status boundary

Run

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_linear_chord_incidence.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6_rank4_rootless_low_degree_search.sage

sage -python \
  elkies-k3/scripts/classify_e6_rank4_det78_niemeier_frames.sage \
  --rootless-obstruction
```

The generated certificate is
[`../artifacts/generated-results/elkies-k3-e6-rank4-linear-chord-incidence-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-linear-chord-incidence-v1.json).
The complete low-degree search artifact is
[`../artifacts/generated-results/elkies-k3-e6-rank4-rootless-low-degree-search-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-rootless-low-degree-search-v1.json).
The genus-wide obstruction artifact is
[`../artifacts/generated-results/elkies-k3-e6-rank4-det78-rootless-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-e6-rank4-det78-rootless-obstruction-v1.json).

Proved here: the exact incidence decomposition; the displayed `QQ(k)`
parameterization of the unordered genus-zero quotient; the genus-one ordered
normalization `52a2`, its complete rational point set, and a quadratic
nondegenerate point; exact arithmetic rank two over `QQ(k)` and geometric rank
four over the ordered incidence field; generic geometric `rho=19`, saturated
determinant `78`, and passage of the necessary rootless-MW17 determinant
screen; the full geometric integral NS marking; and absence of a rootless
child in the complete zero-neutral old-degree two, three, and four shells.
The cited Theorem H5 additionally proves genus-wide rootless nonexistence at
J2 level.

Not proved here: a rank-four family over `QQ(k)`, which this incidence cannot
supply; a classification of all J2 frame classes; or a J1 classification
modulo automorphisms of the K3 surface.

The E6 rational-surface normal form and quadratic-base-change rank splitting
are compatible with Kimura's construction; the incidence component, sections,
saturation, and determinant computation are the exact calculations recorded
here.  See Y. Kimura, [*F-theory models on K3 surfaces with various
Mordell--Weil ranks*](https://arxiv.org/abs/1802.05195), and M. Schuett--T.
Shioda, [*Elliptic surfaces*](https://projecteuclid.org/journals/algebraic-geometry/volume-3/issue-4/Elliptic-surfaces/10.14231/AG-2016-020.full), for the background.
