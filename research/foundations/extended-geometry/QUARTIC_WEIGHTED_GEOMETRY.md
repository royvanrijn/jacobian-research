# Quartic weighted geometry program

This note follows the proof order for the exact quartic-sheet model before any
claim is generalized to arbitrary seeds.

The map below is the point `alpha=-1/2` on the normalized quartic weighted
seed line

\[
 H_\alpha(W)=
 \alpha W^4-(1+2\alpha)W^3+(1+\alpha)W^2
 =W^2(W-1)(\alpha W-\alpha-1),                       \tag{0.1}
\]

where `H_alpha'(1)=-1`.  Exact degree and weighted admissibility remove
`alpha=0,1`, respectively.  Thus the presentation space is the curve

\[
 \mathcal W_4=\operatorname{Spec}
 k[\alpha,(\alpha(\alpha-1))^{-1}].                  \tag{0.2}
\]

For `alpha!=-1`, zero is a double root and the fourth primitive root is
`rho=1+1/alpha`.  At `alpha=-1` it joins the zero cluster, giving the triple-
zero seed.  After dividing the universal factor `C^2`, the leading
target-boundary trace is, up to a unit,

\[
 (\alpha+1)^2\bigl(4(\alpha+1)A-B^2\bigr).           \tag{0.3}
\]

It is reduced for every split seed.  At `alpha=-1` the displayed leading
trace vanishes, the exact `C`-order rises from two to three, and the next
trace is a cube; its boundary-contact index is three.  In particular no
quartic weighted seed has the six-fold contact of degree-four cancellation.

The scoped exhaustion question, the proof under straightened-suspension
hypotheses, and the remaining escape hatches are stated in
[Degree-four marked-root classification](DEGREE_FOUR_MARKED_ROOT_CLASSIFICATION.md).

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

## Quartic inverse on `C!=0`

For target `(A,B,C)`, the inverse equation is

\[
E(W)=W^2-W^4-2BCW+AC^2=0.
\]

This is twice the normalized pencil attached to `H_(-1/2)`; multiplying an
inverse equation by a nonzero scalar does not change its marked-root cover.

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

## Rank-two degree-drop completion

The seed \(H=(W^2-W^4)/2\) is also the specialization
\((\kappa,\tau)=(-5,0)\) of the normalized degree-five seed surface.  Thus it
is the natural classical viability test on the degree-drop divisor.  The
uniform rank-two shear specializes to

\[
s_2=-\frac{261}{28}.
\]

Substitution before parameter-field simplification gives a polynomial
Hamiltonian and an exact polynomial symplectic map

\[
(R,T,D,S):\mathbb A^4\longrightarrow\mathbb A^4.
\]

The checker verifies all six Poisson brackets, the polynomial canonical
adapted-coordinate change, and polynomial left--right equivalence to
\(G\times\operatorname{id}_{\mathbb A^1}\).  Hence the completed map has
generic degree four.  The two stored quartic points transport to two distinct
points over \((R,T,D,S)=(2,0,0,0)\), so noninvertibility survives the
degree-drop specialization.

The fiber differential orders are

\[
\deg_Z(S,T)=(4,3),
\]

one step below the generic degree-five orders \((5,4)\).  In adapted
coordinates the four output degrees are \((3,15,30,18)\), and after
substitution into standard canonical coordinates they are
\((3,18,37,22)\).

This proves **classical viability only**.  It does not construct an
endomorphism of \(A_2\).

### Rebuilt \((4,3)\) restricted deformation complex

For this specialization the exact sparse symbols have

\[
 (\deg_B S,\operatorname{ord}_Z S)=(22,4),\qquad
 (\deg_B T,\operatorname{ord}_Z T)=(18,3),
\]

where \(\deg_B(X^iQ^jZ^k)=i+j+3k\).  Rebuilding the normal-ordered
correction spaces with the inherited rule

\[
 \deg_B(F_n)\leq \deg_B(F)-2n,\qquad
 \operatorname{ord}_Z(F_n)\leq \operatorname{ord}_Z(F)-n
\]

changes the previous boundary diagnosis.  The full parity-preserving
\(\hbar^3\) equation has \(615\) columns, rank \(592\), and a
23-dimensional affine solution space.  At \(\hbar^5\), however, only the
120-dimensional \(S_4\) space remains; \(T_4=0\) in this filtration.  The
allowed current corrections have rank \(115\).  After adjoining all 299
constant, linear, and quadratic coefficient vectors obtained by varying the
complete \(\hbar^3\) affine family, the span has rank \(143\), while adjoining
the constant fifth-order defect raises the rank to \(144\).

An exact six-term dual functional supported on

\[
 X^{12},\ X^{13}Q,\ X^{14}Z,\ X^{14}Q^2,\
 X^{15}QZ,\ X^{15}Q^3
\]

annihilates the full rank-143 span and pairs to one with the defect.  It is
therefore a gauge-invariant dual cocycle for this restricted
parity-preserving complex, and proves that **every** lift in the rebuilt
affine \(\hbar^3\) family is obstructed at \(\hbar^5\).

The unrestricted first-correction calculation is not empty.  Its
\(600+324\) columns have rank \(890\) and nullity \(34\).  The complete
admissible target-Hamiltonian gauge has rank \(14\), leaving a
20-dimensional quotient.  Projection of its 210 Maurer--Cartan quadrics to
the next cokernel has rank \(21\).  Five coordinate axes survive this
quadratic projection, but all five give inconsistent coupled
\(\hbar^2/\hbar^3\) systems.

The bounded low-support closure test is also exact.  The quadrics vanish on
the coordinate \(\mathbf P^4\) with coordinates

\[
 (x_0,x_1,x_7,x_8,x_{17}).
\]

Outside its ten coordinate lines, there are exactly nine rational
exact-support-two directions and no algebraic ones.  Every one of those nine
directions fails the coupled \(\hbar^2/\hbar^3\) equations.  A uniform linear
relaxation of the genuine third-order equations reduces the whole
\(\mathbf P^4\) to the necessary residual plane spanned by

\[
 x_7+2x_0,\qquad
 x_8+\frac{28}{9}x_0,\qquad
 x_{17}+\frac{824}{81}x_0.
\]

The three displayed basis directions fail individually.  Reapplying the
uniform relaxation does not shrink their projective plane.  The genuine
determinantal compatibility calculation nevertheless collapses exactly:
after quotienting by the fixed third-order correction image, all 23
lower-lift columns span one fixed six-dimensional space, and the remaining
obstruction is

\[
 \frac{(21a+28b+64c)^3}{21^3}.
\]

Coefficientwise, the projected kernel couplings, linear right side, and
cubic right side are all divisible by \(21a+28b+64c\).  Consequently the
exact nonzero-scale \(\hbar^3\) locus on the residual plane is the rational
projective line

\[
 21a+28b+64c=0,
\]

with basis \((4,-3,0)\), \((0,16,-7)\).  This is a genuine odd-correction
resonance, not a numerical sample.

The complete resonance line is obstructed at fourth order.  Parameterize it
by

\[
 e_0+t e_1,\qquad e_0=(4,-3,0),\quad e_1=(0,16,-7).
\]

Over \(\mathbb Q(t)\), the joint second/third correction family has dimension
38.  After adjoining all 779 constant, linear, and quadratic coefficient
vectors from that complete lower-lift family to the 120 allowed fourth
corrections, the span has rank \(143\); adjoining the fourth-order defect
raises it to \(144\).  The resulting six-term dual cocycle has sole
denominator factor \(t\).  Exact audits at \(t=0\) and at projective infinity
both give the same rank jump.  Hence every point of the projective resonance
line is obstructed at \(\hbar^4\).

This closes the entire bounded low-support sector: the nine isolated
support-two directions fail at \(\hbar^3\), and the only coordinate
\(\mathbf P^4\) component reduces to the resonance line and fails at
\(\hbar^4\).  Possible higher-support components of the 20-variable
quadratic scheme remain separate.

These are exact statements about the displayed symbols, filtration, and
normal ordering.  They do not rule out a different ordering, larger
correction spaces, nonstandard polarization, Hamiltonian reduction, or a
different classical map, and hence do not prove \(DC_2\).  The rational
six-term cocycle and all ranks are recorded in
[`../artifacts/generated-results/quartic_degree_drop_quantization.json`](../artifacts/generated-results/quartic_degree_drop_quantization.json).

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

Apply the
[generic Morse-slice lemma](../verified/UNIVERSAL_SYMMETRIC_MONODROMY.md#the-generic-morse-slice-lemma)
to the quartic `H`.  Some vertical line `s=sigma` is the cover of a Morse
quartic, hence has monodromy `S_4`; its monodromy embeds in that of the
two-parameter pencil.  Therefore the geometric monodromy of `E` over
`bar(Q)(s,t)` is `S_4`, and the arithmetic monodromy over `Q(s,t)` is also
`S_4`.

The following exact local branch types are useful consistency checks and
describe the special discriminant strata, but they are no longer separate
inputs to the monodromy proof:

- at the smooth point `r=1`, `(s,t)=(-1,-2)`, the factorization is
  `-(W-1)^2(W^2+2W+2)`; a transverse `t`-slice splits the double root and gives
  a transposition;
- at either cusp `r=+/-1/sqrt(6)`, the local equation in `W=r+h` and
  `t=t_0+epsilon` starts with `epsilon+c h^3`, so the three roots at the triple
  root give a 3-cycle;
- at the node, `s=0`, `t=-1/4+epsilon` gives
  `-(W^2-1/2)^2+epsilon`, so the two double roots split simultaneously and
  give a double transposition.

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
python3 scripts/audit_quartic_independent.py
```

The former standalone audit narrative is retained in
[archive/geometry-support](../archive/geometry-support/README.md).
