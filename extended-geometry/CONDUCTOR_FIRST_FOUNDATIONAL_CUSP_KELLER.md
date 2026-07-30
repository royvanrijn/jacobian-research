# A conductor-first cusp realization of the foundational Keller map

> **Status.** Starting from the cusp conductor
> \(k[u^2,u^3]\subset k[u]\), one recovers the cubic tangent incidence,
> its finite marked-root algebra, the reconstruction pole, and the
> distributed determinant ledger of the foundational weighted Keller map.
> All polynomiality equations and the determinant identity hold
> simultaneously.  Thus the requested existence criterion is met: the
> scaled foundational degree-three map is a Keller family whose selected
> repeated-root discriminant has nontrivial cusp-conductor gluing.  This is
> not a new construction mechanism; the conductor-first data reconstructs
> the known weighted tangent suspension.

This reverses the usual presentation.  The cusp conductor is the input and
the seed polynomial is reconstructed from it.

## 1. Begin with the conductor

Let

\[
A_{\rm cusp}=k[u^2,u^3]\subset B=k[u].                     \tag{1.1}
\]

Its normalization quotient has basis \(u\), and its conductor in \(B\) is

\[
\mathfrak c=u^2k[u].                                       \tag{1.2}
\]

Choose normalized target coordinates

\[
S=-3u^2,\qquad V=-2u^3.                                   \tag{1.3}
\]

They satisfy

\[
4S^3+27V^2=0,                                             \tag{1.4}
\]

and generate \(A_{\rm cusp}\).  Thus (1.3) is the prescribed finite
normalization map, with nontrivial conductor at \(u=0\).

Set

\[
W=\frac13+u,\qquad
s=\frac13+S,\qquad
t=\frac1{27}+\frac S3+V.                                  \tag{1.5}
\]

The unique cubic tangent primitive in this normalization is

\[
H(W)=W^2(1-W).                                             \tag{1.6}
\]

Indeed, on (1.5),

\[
s=H'(W),\qquad t=WH'(W)-H(W).                             \tag{1.7}
\]

Consequently the conductor-first parametrization is exactly the critical
restriction of the tangent incidence.

## 2. The finite marked-root algebra descends

Consider the cubic inverse equation

\[
E(W;s,t)=H(W)-sW+t.                                       \tag{2.1}
\]

It is monic up to the constant leading coefficient \(-1\), so it defines a
finite degree-three marked-root algebra over \(k[s,t]\).  Its discriminant
is

\[
\Delta(s,t)
=-4s^3+s^2+18st-27t^2-4t.                                \tag{2.2}
\]

After (1.5),

\[
\Delta=-4S^3-27V^2.                                      \tag{2.3}
\]

Thus the reduced discriminant divisor is precisely the cusp (1.4), and its
coordinate algebra descends through the prescribed conductor.  At the
conductor point,

\[
E\left(W;\frac13,\frac1{27}\right)
=-\frac{(3W-1)^3}{27},                                    \tag{2.4}
\]

so all three marked roots meet at \(W=1/3\).

## 3. Restore the transverse critical coordinate

Introduce \(\gamma\) and define

\[
\begin{aligned}
s&=H'(W)+\gamma
  =2W-3W^2+\gamma,\\
t&=W\bigl(H'(W)+\gamma\bigr)-H(W)
  =W\gamma+W^2-2W^3.                                     \tag{3.1}
\end{aligned}
\]

The plane-core Jacobian is

\[
\det\frac{\partial(s,t)}{\partial(W,\gamma)}=-\gamma.
                                                                  \tag{3.2}
\]

Hence \(\gamma=0\) is the selected critical divisor, and its finite
birational image is the conductor-glued cusp of Section 1.

The inverse derivative is

\[
E_W=H'(W)-s=-\gamma.                                      \tag{3.3}
\]

The required reconstruction coordinate will therefore carry a pole on
\(\gamma=0\), including above the cusp conductor point.

## 4. Solve polynomiality and the determinant ledger together

On affine source space put

\[
\begin{aligned}
\gamma_0&=1-\frac32xy-\frac12x^2z,\\
W_0&=(1+xy)\gamma_0,\\
C_0&=x\gamma_0.                                           \tag{4.1}
\end{aligned}
\]

Substitute \((W,\gamma,C)=(W_0,\gamma_0,C_0)\) into (3.1), and define
the target coordinates \(A,B,C\) by

\[
BC=s,\qquad AC^2=t,\qquad C=C_0.                          \tag{4.2}
\]

The apparent divisions in (4.2) cancel exactly.  The resulting polynomial
map is

\[
G=(A,B,C)=\left(F_1,\frac{F_2}{2},\frac{F_3}{2}\right),    \tag{4.3}
\]

where

\[
\begin{aligned}
F_1&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
F_2&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
F_3&=2x-3x^2y-x^3z.                                      \tag{4.4}
\end{aligned}
\]

Thus polynomiality is solved, rather than assumed.

For the source and target vertical charts

\[
\rho=(W_0,\gamma_0,C_0),\qquad
\mu(A,B,C)=(BC,AC^2,C),                                   \tag{4.5}
\]

the determinants are

\[
J_\rho=-\frac12x^3\gamma_0^2,\qquad
J_{\rm core}=-\gamma_0,\qquad
J_\mu=-C^3.                                               \tag{4.6}
\]

Since \(C=x\gamma_0\),

\[
J_{\rm core}J_\rho
=\frac12C^3
=(J_\mu\circ G)\left(-\frac12\right).                     \tag{4.7}
\]

The chain rule therefore gives

\[
\boxed{\det DG=-\frac12.}                                 \tag{4.8}
\]

This is the simultaneous determinant-ledger and polynomiality solution.

## 5. Reconstruction pole and noninjectivity

Equations (3.3) and (4.1) give

\[
x=\frac C\gamma=-\frac C{E_W}.                            \tag{5.1}
\]

On the \(C\ne0\) part of the repeated-root discriminant this has a genuine
pole.  The cusp-conductor point (2.4) lies on that same divisor, so the
required reconstruction pole meets the nontrivial conductor gluing.

The map is not merely Keller.  The three source points

\[
(0,0,-1/4),\qquad
(1,-3/2,13/2),\qquad
(-1,3/2,13/2)                                             \tag{5.2}
\]

all map under \(G\) to

\[
(-1/4,0,0).                                               \tag{5.3}
\]

Thus (4.3) is an everywhere-etale, noninjective polynomial map with the
declared conductor-glued marked-root discriminant.

## 6. What the success means

The existence half of the requested criterion is now literal:

\[
\boxed{
\text{normalization }\mathbb A^1
+\text{ cusp conductor}
+\text{ finite marked roots}
+\text{ reconstruction pole}
+\text{ polynomial Keller ledger}.}                       \tag{6.1}
\]

However, it does not yield a new stable class or construction mechanism.
The distributed source/target ledger in (4.6) is precisely the weighted
tangent suspension already underlying the foundational map.  This also
explains why the
[separated one-chart obstruction](CONDUCTOR_FIRST_ONE_CHART_OBSTRUCTION.md)
does not apply: the source is not obtained by simply localizing the
normalized conductor line.  The pole cancellation is distributed between
two vertical charts, and the final source remains \(\mathbb A^3\).

A genuinely new conductor mechanism must therefore have conductor gluing
that changes the vertical ledger or its polynomiality equations, rather
than merely decorating the discriminant of a weighted tangent seed.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_conductor_first_foundational_cusp_keller.py
```

The command writes
`artifacts/generated-results/conductor_first_foundational_cusp_keller.json`.
