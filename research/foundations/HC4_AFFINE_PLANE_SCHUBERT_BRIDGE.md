# Affine-plane and second-order flatness obstruction for HC4

## Status

This note continues `HC4-MR` on the regular `[4]` stratum.  It does **not**
assume that the moving Jordan frame is affine parallel.

The combined relative-nilpotent result is registered as `HC4MR1` in
`MATH_STATUS.json`.  The labels `HC4RSD77--80` are retained as local proof-map
identifiers for the final packet.

**Correction, 5 September 2026.** The
[motion-frame audit](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md) shows that frozen
normalization controls `pq/a^2`, not `pq`. The extra equation `d(pq)=0`
used below is not established by the geometric hypotheses. `HC4RSD79`
therefore remains conditional as a local augmented-system statement.
The positive maximal-motion sign is now excluded by differentiated branch
identities; the negative sign is now excluded globally by
[HC4MRA2](HC4_NEGATIVE_MOTION_POLYNOMIAL_OBSTRUCTION.md). These two replacement
arguments restore the full `HC4MR1` reduction.

The middle Jordan distribution is an affine-plane foliation.  Its first-order
Grassmann motion is not automatically Schubert: one extra scalar survives.
The old augmented flatness prolongation excludes motion only if
`d(pq)=0` is supplied as an additional hypothesis.
The lower-motion split is closed by a
degree-one incidence argument for affine-hyperplane foliations.

> **Theorem HC4RSD77 — affine-plane middle foliation.**  Let
> \[
> S=\operatorname{Hess}\psi,\qquad
> T=\operatorname{Hess}A,\qquad
> N=S^{-1}T
> \]
> be a polynomial regular nilpotent `[4]` pencil with `det S` a nonzero
> constant, and assume the complete Jordan flag is Frobenius.  On the generic
> constant-rank locus put
> \[
> E_2=\ker N^2=\operatorname{im}N^2.
> \]
> Then `E_2` is autoparallel for the ambient flat affine connection.  Hence
> every connected leaf of `E_2` is an open subset of an affine two-plane.

> **Theorem HC4RSD78 — exact first-order Grassmann normal form.**  In an
> `S`-adapted Jordan frame, the two transverse derivatives of the direction
> plane `E2` have the form
> \[
> A_3=\begin{pmatrix}0&-(a+q)/2\\0&0\end{pmatrix},
> \qquad
> A_4=\begin{pmatrix}a&r\\0&q\end{pmatrix}.
> \tag{0.1}
> \]
> Thus first order alone does not force the rank-one Schubert condition
> `q=0`.

> **Conditional theorem HC4RSD79 — augmented maximal-motion obstruction.** On the
> smooth Gauss-rank-two locus of the linearly-independent regular `[4]`
> packet, assume additionally that `d(pq)=0` in the adapted frame.
> The zero-curvature
> prolongation of the full Hessian/Codazzi, Frobenius, quasi-translation, and
> unit-volume system is then inconsistent.  Consequently the final regular
> `[4]` packet has no such locus satisfying that extra hypothesis.

> **Theorem HC4RSD80 — lower-motion hyperplane-pencil closure.**  In the
> linearly-independent regular `[4]` packet, the projective source-kernel line
> cannot have rank at most one.  Rank zero is the fixed-kernel branch.  In
> rank one, flatness gives the exact split `pr=0`: the component `r=0` has a
> constant middle plane, while `p=0` has an affine-hyperplane foliation whose
> leaf hyperplanes form a pencil.  The pencil either gives a constant linear
> invariant or forces the middle plane to be constant.  Both alternatives
> return to the already closed linearly-dependent packet.

`HC4RSD80` addresses the lower-motion alternatives. The corrected
maximal-motion closure uses HC4MRA1 and HC4MRA2; the latter requires global
polynomiality on an affine leaf. The complete `HC4-MR` reduction follows
from these replacement arguments and the earlier proof map.

## 1. Why `E_2` is affine

Choose an `S`-adapted Jordan frame

\[
Ne_1=0,\quad Ne_2=e_1,\quad Ne_3=e_2,\quad Ne_4=e_3
\]

with anti-diagonal `S`.  Let

\[
\nabla_{e_i}e_j=\Gamma^k_{ij}e_k.
\]

The exact linear Codazzi system for both `S` and `T=SN`, together with
Frobenius of `E_2=<e_1,e_2>`, forces

\[
\Gamma^3_{ij}=\Gamma^4_{ij}=0
\qquad(i,j\in\{1,2\}).
\tag{1.1}
\]

Equation (1.1) is the vanishing of the affine second fundamental form of each
leaf.  Thus `E2` is autoparallel and its leaves are affine two-planes.

The same system gives

\[
\Gamma^3_{4,1}=\Gamma^2_{3,1}.
\tag{1.2}
\]

## 2. First-order Grassmann motion

For `L=E2`,

\[
T_L\operatorname{Gr}(2,4)=\operatorname{Hom}(L,V/L).
\]

Put

\[
a=\Gamma^3_{4,1}=\Gamma^2_{3,1},\qquad
r=\Gamma^3_{4,2},\qquad
q=\Gamma^4_{4,2}.
\tag{2.1}
\]

Reading every remaining coefficient from the exact `47`-rank linear system
gives (0.1).  In particular, `q` is not forced to vanish.  For example,
`a=q=1,r=0` gives

\[
A_3=\begin{pmatrix}0&-1\\0&0\end{pmatrix},\qquad
A_4=I_2.
\tag{2.2}
\]

This is an exact formal first-order jet with rank-two projective motion, but
it is not yet required to satisfy flatness at the next order.  Hence it is not
a polynomial `HC4` counterexample.

## 3. Flatness forces two non-Schubert signs

Parameterize the `17`-dimensional first-order solution space and introduce all
`68` directional derivatives of those parameters.  For the row convention

\[
(\Gamma_i)_{jk}=\Gamma^k_{ij},
\]

the flatness equations are

\[
e_i\Gamma_j-e_j\Gamma_i
+\Gamma_j\Gamma_i-\Gamma_i\Gamma_j
-\sum_s(\Gamma^s_{ij}-\Gamma^s_{ji})\Gamma_s=0.
\tag{3.1}
\]

There are `96` scalar equations.  Eliminating the `68` derivative variables
has constant linear rank `48` and leaves four quadratic compatibility
equations, generated up to nonzero scalar multiples by

\[
a(p-q)=0,
\tag{3.2}
\]

\[
4pa-3a^2-4aq+3q^2=0,
\tag{3.3}
\]

where

\[
p=\Gamma^4_{3,3}.
\]

On the maximal-motion open `a!=0`, equations (3.2)--(3.3) give

\[
p=q,\qquad q^2=a^2.
\tag{3.4}
\]

Thus flatness does not send the branch to the Schubert value `q=0`.  It forces
the two non-Schubert signs `q=+a` and `q=-a`.

## 4. The conditional constant-motion calculation

Let `theta^1,...,theta^4` be the dual moving coframe.  The HC4-selected
Gauss-kernel line is represented by

\[
\ell=S e_1=\theta^4.
\]

Because the components of `S` are fixed in the adapted frame,

\[
\nabla_{e_i}\ell=-\sum_j\Gamma^4_{ij}\theta^j.
\tag{4.1}
\]

On the transverse directions `(e3,e4)` and the projective target directions
`(theta^3,theta^2)`, the first-order solution therefore gives

\[
D[\ell]=
\begin{pmatrix}
-p&-s\\0&-q
\end{pmatrix},
\qquad
\det D[\ell]=pq,
\tag{4.2}
\]

where `s=Gamma^4_{4,3}` is irrelevant to the determinant.  The maximal-motion
identity of `HC4RSD72` applies in a differently normalized frozen frame.
The exact transition sends this determinant to `pq/a^2`. It does not
preserve the normalized matrix of `N`, so the sign-only centralizer
argument cannot make its factor constant. The
[transport audit](HC4_MOTION_FRAME_TRANSPORT_AUDIT.md) gives the full matrices.

The following equation must therefore be treated as an **extra hypothesis**:

\[
pq\in K^*,\qquad e_i(pq)=0\quad(1\le i\le4).
\tag{4.3}
\]

Add the four equations `d(pq)=0` before eliminating the derivative variables
from (3.1).  The derivative rank becomes `50`, and one new compatibility
generator is

\[
p(2pa-aq+3q^2)=0.
\tag{4.4}
\]

On `a!=0`, equations (3.4) make `p=q!=0`.  Equation (4.4) gives

\[
a=-3q.
\tag{4.5}
\]

Together with `q^2=a^2`, this gives

\[
8q^2=0,
\]

contradicting characteristic zero and `q!=0` under that extra hypothesis.
Equivalently, the three
polynomials (3.2), (3.3), and (4.4), saturated by `a`, have Gröbner basis
`[1]`. This proves the conditional version of `HC4RSD79`.

## 5. Lower-motion closure

The previously proposed Schubert split is not the correct route:

1. first-order geometry allows both `q=0` and `q!=0`;
2. flatness plus maximal motion forces `q=+/-a`, so `q=0` is impossible;
3. independently established constancy of `pq` would eliminate both signs.

The negative sign survives the local transport audit but is excluded by the
global polynomial-leaf proof in HC4MRA2. The fixed and linearly-dependent packets retain
their reductions to `HC2` or the exact `JC2` cotangent packet.

It remains to close projective source-kernel motion of rank at most one.  The
exact projective derivative of the line `[e1]`, with columns indexed by
`e1,...,e4` and rows by `e2,e3,e4`, is

\[
D[e_1]=
\begin{pmatrix}
0&0&a&b\\
0&0&0&a\\
0&0&0&0
\end{pmatrix}.
\tag{5.1}
\]

If its rank is zero, then `a=b=0`; flatness gives `q=0`, and `[e1]` is a
constant affine line.  This is the fixed-kernel theorem `HC4RSD65`.

Assume henceforth that its rank is one.  Then

\[
a=q=0,\qquad b\ne0.
\tag{5.2}
\]

After adjoining `d(a)=d(q)=0`, the curvature compatibility ideal is generated
by

\[
pr=0.
\tag{5.3}
\]

Because the generic coordinate ring is a domain, (5.3) gives two components.

### 5.1 The component `r=0`

Under (5.2), the two transverse derivatives of `E2` are

\[
A_3=0,\qquad
A_4=\begin{pmatrix}0&r\\0&0\end{pmatrix}.
\tag{5.4}
\]

The derivatives along `E2` already preserve `E2` by `HC4RSD77`.  Thus `r=0`
makes `E2` parallel in every ambient direction.  It is a constant two-plane,
so the associated source-kernel field `e1` has two constant linear invariants.
This is inside the linearly-dependent branch closed by `HC4RSD70`.

### 5.2 The component `p=0`

The exact affine second fundamental form of

\[
E_3=\ker N^3=\langle e_1,e_2,e_3\rangle
\]

along `E3` is

\[
\operatorname{II}_{E_3}=
\begin{pmatrix}0&0&0\\0&0&0\\0&0&p\end{pmatrix}e_4.
\tag{5.5}
\]

Hence `p=0` makes every generic `E3` leaf an open subset of an affine
hyperplane.  We use the following elementary global lemma.

> **Affine-hyperplane pencil lemma.**  Let a codimension-one algebraic
> foliation of a dense open subset of `A4` have generic leaves which are open
> subsets of affine hyperplanes.  Then the projective closures of the generic
> leaves form a line in the dual projective space `(P4)*`.

Indeed, the affine tangent hyperplane through a point is determined
algebraically by the distribution, so it defines a rational leaf-hyperplane
map

\[
h:\mathbb P^4\dashrightarrow C\subset (\mathbb P^4)^*.
\]

The image cannot be a point because one affine hyperplane cannot contain a
dense open subset of `A4`; hence its closure `C` is an irreducible curve.  Let

\[
I_C=\{(x,H):x\in H,\ H\in C\}.
\]

The closure of the graph of `h` is a dimension-four irreducible subvariety of
the irreducible dimension-four incidence variety `I_C`; hence it equals
`I_C`.  The projection `I_C -> P4` is therefore birational.  On the other
hand, its fiber over a general point `x` is the hyperplane section

\[
C\cap x^\perp,
\]

whose length is `deg C`.  Therefore `deg C=1`, proving the lemma.  This
argument is unchanged after extending the constant field to its algebraic
closure; the resulting line and its span descend because the foliation is
defined over `K`.

Apply the lemma to `E3`.  Write an affine hyperplane in the pencil as
`H(t)` and let `lambda(t)` denote its linear direction conormal.  Either this
direction is constant, or, after a rational parameter `t` on the dual line,

\[
\lambda(t)=\lambda_0+t\lambda_1.
\tag{5.6}
\]

If the direction is constant, then `E3` is a constant three-plane direction,
so `e1` has a constant linear invariant and `HC4RSD70` applies.

Suppose the direction moves, so `lambda0,lambda1` are linearly independent.
Equations (5.1)--(5.2) and the middle-plane
matrices show that `[e1]` and `E2` are constant along every `E3` leaf and vary
only transversely.  Choose a representative `v(t)` of `[e1]`.  Since `b!=0`,

\[
E_2(t)=\langle v(t),v'(t)\rangle.
\tag{5.7}
\]

Both vectors lie in `E3(t)=ker lambda(t)`, so

\[
\lambda v=0,\qquad \lambda v'=0.
\tag{5.8}
\]

Differentiate the first identity and use the second:

\[
0=(\lambda v)'=\lambda_1v+\lambda v'
\quad\Longrightarrow\quad
\lambda_1v=0.
\]

Equation (5.6) then also gives `lambda_0v=0`.  Hence `v(t)` takes values in a
fixed space, and differentiating once more gives

\[
v(t),v'(t)\in
W:=\ker\lambda_0\cap\ker\lambda_1.
\]

The fixed space `W` has dimension two, and `b!=0` makes `v,v'` independent.
Thus `E2(t)=W` is constant, contradicting `r!=0`.  Therefore the component
`p=0,r!=0` is empty.  This proves `HC4RSD80`.

Combining rank zero, the `r=0` component, and the contradiction on
`p=0,r!=0` closes every lower-motion boundary of the final regular `[4]`
packet.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_affine_plane_bridge.py
.venv/bin/python scripts/verify_hc4_affine_plane_prolongation.py
```

The first checker proves `HC4RSD77--78`.  The second constructs all curvature
equations, performs both derivative eliminations, verifies the residual frame
gauge, checks the saturated unit ideal under the extra `d(pq)=0` assumption,
and certifies the
complete local flag tensors and the split `pr=0` used in `HC4RSD80`.  The
degree-one incidence argument in Section 5.2 is the global proof step and is
not replaced by a bounded computation. The corrected transport, valid
positive-sign closure, and compatible negative-sign jet are replayed by
`scripts/verify_hc4_motion_frame_transport.py`. The final global obstruction
has its own written proof and
`scripts/verify_hc4_negative_motion_polynomial_obstruction.py`.
