# Determinant 800: the full marked curve is X_0(200)

The row K3-643eed54b4d256fc has literal transcendental lattice

\[
T=U(2)\oplus\langle200\rangle=2T_0,\qquad
T_0=U\oplus\langle100\rangle.
\]

There is no rank-19 K3 over Q with this actual transcendental lattice and
its full geometric NS rationally marked. Its full marked period curve is
the canonical Q-curve X_0(200), of genus 19, and the degree-four forgetful
map to X_0(50) excludes rational noncuspidal points.

This is the last of the historical 21 coarse-low-genus diagnostic rows to
be excluded after the two real-place filters. It is not an exhaustion of
the different-NS foundry: 721 rows of the historical queue remain unresolved.

## The full integral group

Put N=50 and use coordinates (e,h,f), so the quadratic form of T_0 is
ef+Nh². Its even Clifford order is the split Eichler order
R^0(N), with basis I, N E12, E11, E21. This follows directly from the
Clifford products of the two isotropic generators and h. It is conjugate
to R_0(N) by J=[[0,-1],[1,0]]. Scaling T_0 by two does not change its
integral isometry group, though it does change the discriminant kernel.

The trace-zero representation sends (e,h,f) to
[[Nh,Ne],[f,-Nh]]. For g=[[a,b],[c,d]] of determinant delta, conjugation
acts by

\[
A(g)=\delta^{-1}\begin{pmatrix}
a^2&-2ab&-b^2/N\\
-ac&ad+bc&bd/N\\
-Nc^2&2Ncd&d^2
\end{pmatrix}.
\]

All isometries of T_0 induce automorphisms of its Clifford order. By
Skolem–Noether, their projective representatives normalize this order.
Locally at p^n exactly dividing N, the order is the intersection of the
two maximal orders at the endpoints of a length-n segment in the
Bruhat–Tits tree. A normalizer either preserves or exchanges the endpoints.
At other primes it preserves the single maximal order. Thus the positive
rational projective normalizer has the four exact-divisor cosets
Gamma_0(50) W_d for d=1,2,25,50. Representatives, in R_0 coordinates, are

\[
W_2=\begin{pmatrix}2&1\\50&26\end{pmatrix},\quad
W_{25}=\begin{pmatrix}25&1\\600&25\end{pmatrix},\quad
W_{50}=\begin{pmatrix}0&-1\\50&0\end{pmatrix}.
\]

Their determinants are respectively 2,25,50. The projective representation
of O also permits a common sign on A(g); both signs must be checked.
The replay calculates the literal action on A_T, whose dual-generator
denominators are (2,200,2). The Gamma_0(50) image has order four.
Of all eight signed cosets, only +Gamma_0(50) contains the identity
discriminant action. No reflection or additional Atkin–Lehner coset has
been omitted.

Inside Gamma_0(50), the kernel is exactly

\[
\Gamma=\Gamma_0(50)\cap\Gamma(2).
\]

The replay checks every generator of the conjugate of Gamma_0(200)
inside Gamma_0(50); all act trivially. This subgroup has index four,
matching the discriminant image order, so containment is equality.
Conjugation by diag(1/2,1) sends Gamma to Gamma_0(200).

## The Q-model, including determinant units

The complex group alone would not determine the arithmetic descent.
At 2 write an Eichler unit in R^0 coordinates as
g=[[a,50b],[c,d]], with a,d odd and delta odd. The displayed action shows
that A(g) is stable on A_(T,2) precisely when b and c are even. The
checks involving the h-dual vector give this condition; with it all the
other entries of (A(g)-I)G^-1 are 2-integral.

No -A(g) is stable: the off-diagonal conditions again force b,c even,
whereas the h coefficient of A(g) is then 1 modulo eight, not -1.
The W_2 coset, with either sign, exchanges the two isotropic directions
modulo two, so cannot be stable. The replay exhausts all 1024 unit
residues modulo eight, in each signed coset, and obtains 256 stable
residues in the positive unit coset and none in the other three. Their
determinants cover all four odd residues. The rational action formulas
show that these tests depend only on those residues, so this is a full
local calculation, not a finite-depth approximation to an unknown kernel.

In particular no determinant-minus-one orthogonal isometry is stable
locally at 2. Hence the orientation character of a full rational NS
marking is trivial. It is one global quadratic character, not a choice
of unrelated local signs; compatibility of the transcendental cohomology
representations gives the same character at every prime.

At 5, all Eichler units are stable, while W_25 has multiplier -1 on
the cyclic order-25 part of A_T and is excluded. At other primes the
level is maximal. After the same rational conjugation, the stable GL
level is therefore exactly the usual Gamma_0(200) level, including all
determinant units. This identifies the canonical model as X_0(200)
over Q. There is no unexamined quadratic twist or orientation quotient.

The passage from a full rational NS marking to this canonical orthogonal
Shimura model uses the same period-map inputs as the
[corrected NS0031 proof](NS0031_FULL_STABLE_MARKING_OBSTRUCTION_2026-09-14.md):
[Dolgachev, Proposition 3.3](https://arxiv.org/abs/alg-geom/9502005)
for the stable group and
[Rizov, Section 3.9 and Theorem 3.16](https://arxiv.org/abs/math/0508018)
for the canonical arithmetic period map. Here the local 2-adic calculation
excludes the orientation ambiguity that required a separate quotient
for NS0031.

## Rational-point obstruction

The rational cyclic-isogeny classification allows degrees 1 through 19,
21,25,27,37,43,67,163; 50 is absent. See the theorem recalled in
[Banwait, Cyclic isogenies of elliptic curves over fixed quadratic fields](https://arxiv.org/abs/2206.08891).
Thus X_0(50)(Q) has no noncuspidal point. The canonical degree-four map
X_0(200) to X_0(50) preserves cusps and noncuspidal periods, so neither
does X_0(200). In particular no rational non-CM marked period exists.

This excludes only the stated full rank-19 marking over Q. It does not
exclude geometric realizations, partial markings, or larger fields.

## Replay and assurance boundary

```sh
sage -python research/elkies-k3/scripts/certify_det800_marking.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det800-full-marking-v1.json)
contains the exact finite discriminant image, signed-coset decisions,
kernel generator checks, genera, and local determinant-unit calculation.
The normalizer completeness and canonical period-map arguments are
explicit written theorem inputs. Independent implementation, formal
verification and external review are not claimed. No equation or
specialization search was run. The preflight retains the failed use of a
global integrality assertion on a 2-adic unit; the corrected calculation
tests 2-integrality instead.
