# Quartic weighted geometry program

This note follows the proof order for the exact quartic-sheet model before any
claim is generalized to arbitrary seeds.

Put

\[
u=1+3xy,\qquad \gamma=1-4xy-x^2z,
\]

and define

\[
G=\left(
\frac{2u+u^2-3u^4\gamma^2}{x^2},
\frac{1+u-2u^3\gamma^2}{x},
x\gamma
\right).
\]

Both apparent quotients cancel. Exact expansion gives coordinate degrees
`(12,11,4)`, determinant `-6`, and

\[
G(1,0,0)=G(-1,0,2)=(0,0,1).
\]

## Implementation record

- [x] Generalized weighted-seed builder.
- [x] Exact integer quartic map, polynomiality, determinant, degrees, and
  collision.
- [x] `C!=0` root reconstruction and pole classification.
- [x] Discriminant and repeated-root normalization.
- [x] Cusp/node classification.
- [x] Direct `C=0` fiber calculation.
- [x] Explicit escaping sequences for both anticipated components of `S_{G}`.
- [x] Converse properness outside those components.
- [x] Radical certificate for the singular locus.
- [x] Exact image and fiber-cardinality theorem.
- [x] Connectedness and geometric/arithmetic monodromy `S_4`.
- [x] Universal inverse, discriminant-normalization, and `S_n` monodromy
  theorem for arbitrary characteristic-zero seeds.
- [x] Generalization of the image and nonproperness statements to the
  canonical family.
- [x] Boundary and image theorem for one additional simple primitive zero.
- [x] Omitted-value classification for several or repeated extra zeros.
- [x] Boundary saturation theorem for repeated extra zeros.
- [x] Good-reduction finite-field Chebotarev law from universal `S_n`
  monodromy.
- [x] Generic nodal-cuspidal discriminant theorem and generic surjectivity in
  inverse degree at least five.

## Quartic inverse on `C!=0`

For target `(A,B,C)`, the inverse equation is

\[
E(W)=W^2-W^4-2BCW+AC^2=0.
\]

Its derivative satisfies

\[
E'(W)=-2\gamma,
\]

and a simple root reconstructs by

\[
\gamma=BC-W+2W^3,\quad x=C/\gamma,\quad u=W/\gamma,
\]

\[
xy=(u-1)/3,\qquad x^2z=1-4xy-\gamma.
\]

Thus, on `C!=0`, finite source points are in bijection with the simple roots
of `E`; repeated roots are exactly the reconstruction poles. The control target
`(1,0,1)` gives the squarefree polynomial `-(W^4-W^2-1)`, proving generic
degree four.

Run:

```bash
.venv/bin/python scripts/verify_quartic_weighted_map.py
```

The broader seed scan in `WEIGHTED_SEED_SCAN.md` remains exploratory. Its role
is diagnostic; the completed quartic and all-degree theorems use the exact
certificates listed above rather than extrapolation from the scan.

## Discriminant normalization and singularities

Put `s=BC` and `t=AC^2`. Up to the nonzero discriminant factor `-16`, the
quartic discriminant is

\[
\Delta(s,t)=27s^4-36s^2t-s^2+16t^3+8t^2+t.
\]

Its repeated-root normalization is

\[
s=r-2r^3,\qquad t=r^2-3r^4.
\]

The parameter is integral through

\[
r^2-3sr+2t=0,
\]

and is recovered generically from

\[
r={s(12t-1)\over18s^2-4t-1}.
\]

The discriminant curve has exactly three singular points:

- a node `(s,t)=(0,-1/4)`, with normalization preimages
  `r=+/-1/sqrt(2)` and quartic `-(W^2-1/2)^2`;
- two cusps `(s,t)=(+/-sqrt(6)/9,1/12)`, each with one normalization
  preimage `r=+/-1/sqrt(6)`, a triple root, and one simple root.

The node is the first omitted-value candidate on `C!=0`, since both of its
root pairs are repeated. Each cusp retains the source point belonging to its
remaining simple root.

Run:

```bash
.venv/bin/python scripts/verify_quartic_discriminant.py
```

## Monodromy of the inverse cover

The polynomial `E` is irreducible over `Q(s,t)`: viewed in
`Q[s,W][t]`, it is linear and monic in `t`, so any factor independent of `t`
would divide its unit leading coefficient. Thus the four-sheet cover over the
discriminant complement is connected and its geometric monodromy is a
transitive subgroup of `S_4`.

The normalization gives exact representatives for all three relevant local
branch types:

- at the smooth point `r=1`, `(s,t)=(-1,-2)`, the factorization is
  `-(W-1)^2(W^2+2W+2)`; a transverse `t`-slice splits the double root and gives
  a transposition;
- at either cusp `r=+/-1/sqrt(6)`, the local equation in `W=r+h` and
  `t=t_0+epsilon` starts with `epsilon+c h^3`, so the three roots at the triple
  root give a 3-cycle;
- at the node, `s=0`, `t=-1/4+epsilon` gives
  `-(W^2-1/2)^2+epsilon`, so the two double roots split simultaneously and
  give a double transposition.

If `G<=S_4` is the geometric monodromy group, transitivity implies `4` divides
`|G|`, while the cusp 3-cycle implies `3` divides `|G|`. Hence `|G|` is `12`
or `24`. The smooth-point transposition is odd, excluding the unique
index-two subgroup `A_4`; therefore

\[
G=S_4.
\]

The arithmetic monodromy over `Q(s,t)` contains the geometric group and is
also `S_4`.

Run:

```bash
.venv/bin/python scripts/verify_quartic_monodromy.py
```

## Direct fibers on `C=0`

Because `G_3=x gamma`, the affine source over the target plane `C=0` is the
union of `x=0` and `gamma=0`.

On `x=0`,

\[
G(0,y,z)=(3(29y^2+2z),y,0),
\]

so every `(A,B,0)` has the unique point

\[
(0,B,(A-87B^2)/6).
\]

On `gamma=0`, writing `v=xy` gives

\[
A={3(v+1)(3v+1)\over x^2},\qquad
B={3v+2\over x}.
\]

Eliminating `v` yields

\[
(B^2-A)x^2=1.
\]

Thus the `C=0` fibers have three affine points when `A!=B^2` and one when
`A=B^2`. The specialized inverse quartic

\[
E(W)=W^2(1-W^2)
\]

has a double root `0`, the finite boundary root `+1`, and the additional root
`-1`. The latter is not an affine point and is the source of the anticipated
plane component in the nonproperness set.

Run:

```bash
.venv/bin/python scripts/verify_quartic_c0_fibers.py
```

## Both nonproperness components: forward inclusion

Substitution in the quartic discriminant gives

\[
\Delta(BC,AC^2)=C^2Q_4(A,B,C),
\]

where

\[
Q_4=A-B^2+C^2(27B^4-36AB^2+8A^2)+16A^3C^4.
\]

The factor `C^2` is geometrically real rather than a denominator artifact.
The two anticipated components are therefore

\[
V(C)\quad\text{and}\quad V(Q_4).
\]

For every target `(A,B,0)`, the exact root jet

\[
W=-1-BC+\frac{3B^2-A}{2}C^2
\]

(with `C` tending to zero) produces moving targets converging to `(A,B,0)` and
reconstructed sources with `C^2z -> 2`. This proves `V(C) subset S_G`.

For `C!=0`, every point of `V(Q_4)` is represented by

\[
BC=r-2r^3,\qquad AC^2=r^2-3r^4.
\]

Perturbing the root `r` while adjusting `A` keeps the inverse equation exact;
`gamma` tends to zero and `x=C/gamma` diverges. A separate scaled path
`r=BC` verifies the intersection `C=0`, `A=B^2`. Hence
`V(Q_4) subset S_G` as well.

Run:

```bash
.venv/bin/python scripts/verify_quartic_nonproperness_paths.py
```

## Converse properness

Outside `V(C) union V(Q_4)`, the inverse quartic has nonzero discriminant, so
all four roots are finite and simple. Its projective homogenization has value
`-1` at `[W:S]=[1:0]`, excluding an unexamined root-at-infinity chart. The
reconstruction denominators are products of `C` and `E'(W)`, both units on
this open set. Consequently every local root branch reconstructs to a bounded
source branch, proving

\[
S_G=V(C)\cup V(Q_4).
\]

Run:

```bash
.venv/bin/python scripts/verify_quartic_properness_converse.py
```

## Singular locus of the nonproperness hypersurface

For `N=CQ_4`, the radical of the Jacobian ideal is the intersection of three
reduced strata:

\[
(C,A-B^2)
\cap(B,4AC^2+1)
\cap(12AC^2-1,27B^2C^2-2).
\]

They are respectively the intersection of the two nonproperness components,
the lift of the discriminant node, and the two lifts of the discriminant cusps.
The executable certificate checks containment of the Jacobian ideal in this
intersection and verifies that the square of every displayed radical generator
lies in the Jacobian ideal.

Run:

```bash
.venv/bin/python scripts/verify_quartic_singular_locus.py
```

## Exact image and fibers

On the discriminant normalization, the inverse quartic factors as

\[
E(W)=-(W-r)^2\bigl(W^2+2rW+3r^2-1\bigr).
\]

The residual quadratic meets the repeated root when `r^2=1/6` (the cusps) and
itself becomes double when `r^2=1/2` (the node). Consequently the complete
affine fiber table is:

| Target stratum | Fiber size |
|---|---:|
| `C!=0`, `Q_4!=0` | 4 |
| smooth part of `Q_4=0` | 2 |
| either cusp lift | 1 |
| node lift `B=0`, `4AC^2+1=0` | 0 |
| `C=0`, `A!=B^2` | 3 |
| `C=0`, `A=B^2` | 1 |

Thus

\[
G(\mathbb C^3)=\mathbb C^3\setminus
V(B,4AC^2+1),
\]

while

\[
S_G=V(C)\cup V(Q_4).
\]

Run:

```bash
.venv/bin/python scripts/verify_quartic_image.py
```

## Independent central-algebra audit

A dependency-free sparse-polynomial implementation reconstructs the
coordinate divisions, obtains `det DG=-6`, checks both collision points and
the inverse equation, and independently recovers

\[
\operatorname{disc}_W(E)=-16C^2Q_4.
\]

It also verifies

\[
E_{r-2r^3,\,r^2-3r^4}(W)
=-(W-r)^2(W^2+2rW+3r^2-1),
\]

the `C=0` relation `(B^2-A)x^2=1`, and the omitted-node factorization.  Run

```bash
python3 scripts/audit_c14_independent.py
```

The former standalone audit narrative is retained in
[archive/geometry-support](../../archive/geometry-support/README.md).
