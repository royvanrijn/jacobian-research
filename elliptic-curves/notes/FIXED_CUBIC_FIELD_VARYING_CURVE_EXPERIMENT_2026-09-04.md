# Fixed cubic field, varying elliptic curve

Status: **exact class-group-free local-intersection pilot complete; point
realization open**.

This is a mechanism experiment, not a replacement for the MW16/MW17
production lanes.  It fixes a cubic 2-division field and a certified space of
global Kummer classes, varies the elliptic curve, and asks which combinations
remain everywhere locally admissible.

## Fixed-field family

Let

\[
 E_0:y^2=f(x),\qquad f(x)=x^3+Ax+B,
 \qquad K=\mathbf Q(\theta),\quad f(\theta)=0.
\]

For

\[
 \alpha_u=\theta+u\theta^2
\]

the norm polynomial is

\[
\begin{aligned}
F_u(x)={}&x^3+2Au x^2+(A+3Bu+A^2u^2)x\\
          &+B+ABu^2-B^2u^3.
\end{aligned}
\]

Exact expansion gives

\[
 \operatorname{disc}(F_u)
 =\operatorname{disc}(f)(1+Au^2+Bu^3)^2.
\]

Write \(D_u=1+Au^2+Bu^3\).  The experiment excludes \(D_u=0\).  When
\(D_u\ne0\), the inverse change of cubic generator is explicit:

\[
 \theta=\frac{-2Bu^2+(1-Au^2)\alpha_u-u\alpha_u^2}{D_u}.
\]

Thus the roots are not merely roots in isomorphic cubic fields: the formula
fixes the labelled Galois-equivariant identification
\(E_0[2]\simeq E_u[2]\).

## Anchor and fixed global span

The pilot uses the pinned Fermigier--Mestre rank-at-least-20 specialization at
family parameter \(28917/20\).  Its certified short model is

\[
y^2=x^3+Ax+B
\]

with

```text
A = -5750886029903523759416717668139307
B = 167347710468055045100164888198438918505621536951206.
```

The repository already certifies that the twenty pinned points have
20-dimensional image in \(E_0(\mathbf Q)/2E_0(\mathbf Q)\).  Since \(f\) is
irreducible, \(E_0(\mathbf Q)[2]=0\), and their Kummer representatives

\[
 \beta_i=x(P_i)-\theta
\]

form a 20-dimensional subspace

\[
 W\subset H^1(\mathbf Q,E_0[2])
       \simeq H^1(\mathbf Q,E_u[2]).
\]

Every norm identity
\(N_{K/\mathbf Q}(\beta_i)=y(P_i)^2\) is replayed exactly.  This establishes
the global cohomology classes; it does not assert that they lie in the local
Kummer image for a new \(E_u\).

## Full-span local computation

The first bounded policy is all integers \(|u|\le2\).  For each \(u\), the
runner checks:

1. every finite prime dividing \(2\operatorname{disc}(f)D_u\), including all
   newly bad primes of \(E_u\);
2. the real place;
3. the restriction of all twenty \(\beta_i\) to the exact local squareclass
   space of \(K\otimes\mathbf Q_v\);
4. an explicit local-point Kummer basis for \(E_u(\mathbf Q_v)/2E_u(\mathbf
   Q_v)\).

At odd \(p\), the required local-image dimension is
\(\dim E_u(\mathbf Q_p)[2]\), and at \(p=2\) it is one larger.  The number of
local factors of the fixed cubic gives the torsion dimension.  Local points
are searched on both the raw and a \(p\)-minimal model, but a search result is
accepted only when exact local squareclass independence reaches this known
dimension.  Consequently the bounded point enumeration is a witness finder,
not a completeness assumption.  A failure to reach the dimension suppresses
the entire \(W_u\) claim.

Outside the declared support, both integral cubics have good reduction at odd
primes and the fixed classes are unramified, so the good-prime local condition
is automatic.  No class group or global Selmer basis is computed.

The local quotient rows are concatenated and one linear map

\[
 W\longrightarrow\bigoplus_v
 H^1(\mathbf Q_v,E_u[2])/
 \delta_v(E_u(\mathbf Q_v)/2E_u(\mathbf Q_v))
\]

is formed.  The reported space is its kernel.  In particular, the code does
not discard basis elements one at a time: combinations of individually
inadmissible \(\beta_i\) can and do remain eligible for the kernel.

## First result

All local Kummer images completed, and the five curves have pairwise distinct
rational \(j\)-invariants.

| \(u\) | checked finite places | newly bad primes | local-condition rank | \(\dim W_u\) |
|---:|---:|---:|---:|---:|
| \(-2\) | 17 | 5 | 7 | 13 |
| \(-1\) | 14 | 2 | 2 | 18 |
| \(0\) | 12 | 0 | 0 | 20 |
| \(1\) | 15 | 3 | 7 | 13 |
| \(2\) | 16 | 4 | 7 | 13 |

The \(u=0\) positive control recovers the whole 20-dimensional space.  The
sharpest new-curve observation is \(u=-1\): eighteen independent fixed global
classes remain everywhere locally admissible.  This is exactly the signal the
experiment was designed to isolate.

It is not a rank lower bound.  The exact sequence still allows surviving
classes to be represented by nontrivial elements of \(\Sha(E_u)[2]\).

## Point-realization covers

For a surviving span class, use a representative

\[
 \beta=\prod_{i\in I}\beta_i=b_0+b_1\theta+b_2\theta^2
\]

and write

\[
 \gamma=a+b\theta+c\theta^2.
\]

Reducing \(\beta\gamma^2\) modulo
\(\theta^3+A\theta+B\), the projective covering in variables
\((a:b:c:d)\) is the intersection of the two quadrics

\[
 [\theta](\beta\gamma^2)+d^2=0,
 \qquad
 [\theta^2](\beta\gamma^2)+u d^2=0.
\]

For \(d\ne0\), a rational point gives

\[
 x=\frac{[1](\beta\gamma^2)}{d^2},
 \qquad
 y=\frac{s_\beta N(\gamma)}{d^3},
 \qquad s_\beta^2=N(\beta).
\]

The artifact stores an exact basis of every \(W_u\), including the anchor
basis masks, the reduced power-basis coefficients of \(\beta\), and
\(s_\beta\).  These rows are ready for cover minimization and rational-point
search.  The next bounded phase should start with the 18 basis covers at
\(u=-1\), then test rank growth only after mapping any cover point back to
\(E_{-1}\) and certifying independence.

## Reproduction and claim boundary

```sh
sage -python elliptic-curves/cas/run_fixed_cubic_field_curve_family.sage --check
python3 -m unittest elliptic-curves/tests/test_fixed_cubic_field_curve_family.py
```

The canonical output is
[`fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json`](../../artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json).

The output certifies the family identities, fixed Kummer span, complete local
intersections, newly bad-prime coverage, and explicit cover inputs.  It does
not compute the full class group, the full 2-Selmer group of any \(E_u\), a
rational point on a new curve, a new rank lower bound, or an exact rank.
