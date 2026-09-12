# Sparse quartic--sextic collision carriers in the Meng `HC_4` chart

## Status

This note proves a bounded characteristic-zero obstruction in the normalized
Meng descent:

> No homogeneous sextic correction supported on at most four monomials can
> retain the Meng--Yang collision and give constant nonzero Hessian
> determinant.

The logical dependency is the short branch

\[
 \texttt{HC5T1}\longrightarrow\texttt{HC4MS6}.
\]

It is parallel to the sparse-quartic/full-cubic branch
`HC5T1 -> HC4MQ1 -> HC4MCK`.  The present theorem is sextic-only: it does
not combine a sextic with quartic or cubic corrections, and it does not
exclude sextics supported on five or more monomials.

The later mixed theorem `HC4MQS6` excludes zero-gradient sextics supported
on at most four monomials over the 234 quartic principal parts selected by
`HC4MQ1`.  Its short dependency branch is

\[
 \texttt{HC4MQ1},\texttt{HC4MS6}\longrightarrow\texttt{HC4MQS6}.
\]

The joint theorem `HC4JQS4` removes the old sextic-free quartic selection:
it treats every genuinely mixed quartic--sextic collision support of total
size at most four.  Together with the pure quartic and pure sextic
theorems, it gives the combined statement

> No correction \(h_4+h_6\) of total monomial support at most four can
> retain the Meng--Yang collision and give constant nonzero Hessian
> determinant.

There is no cubic term in this statement.  The short dependency branch is

\[
 \texttt{HC4MQ1},\texttt{HC4MS6}\longrightarrow\texttt{HC4JQS4}.
\]

The dense theorem `HC4DCK` removes every support bound when the quartic and
sextic have a constant common Hessian-kernel direction in the precise sense
below.  The determinant layers force that condition automatically when the
sextic Hessian has generic rank three.  Its proof is structural rather than
enumerative:

\[
 \texttt{HC5T1}+\text{Gordan--Noether}+\texttt{HC}_{\le3}
 +\text{the degree-five plane theorem}
 \longrightarrow\texttt{HC4DCK}.
\]

## 1. Collision normalization

Use

\[
 \psi_0=2yr+4xs,\qquad \det\operatorname{Hess}(\psi_0)=64
\]

and the antipodal Meng points

\[
 p=\left(1,-\frac32,6,\frac{81}{8}\right),\qquad -p.
\]

For a homogeneous sextic

\[
 h_6=\sum_e c_e w^e,
\]

its gradient is odd.  The gradients of \(\psi_0+h_6\) collide at
\(p\) and \(-p\) exactly when

\[
 \nabla h_6(p)=-H_0p.
\]

Put \(d_e=c_ep^e\) and multiply the \(i\)-th gradient equation by
\(p_i\).  The collision equations become

\[
 \sum_e d_e e
 =
 \left(-\frac{81}{2},18,18,-\frac{81}{2}\right).
\]

There are 84 degree-six monomials in four variables.  The checker solves
this four-row exponent system before expanding a Hessian determinant.

## 2. Exact support census

Every support of size at most four is enumerated.  The collision census is:

\[
\begin{array}{c|r|r|r|r}
\text{support size}&\text{inconsistent}&\text{isolated}
&\text{lines}&\text{planes}\\ \hline
1&84&0&0&0\\
2&3480&6&0&0\\
3&92936&2344&4&0\\
4&198450&1723488&7562&1
\end{array}
\]

The isolated census is performed modulo \(1000003\).  It is also a
characteristic-zero census.  After doubling the target, every coefficient
minor is bounded by \(6^4\), while Hadamard's inequality bounds every
augmented \(4\times4\) minor by

\[
 6^3\sqrt{2\cdot81^2+2\cdot36^2}<1000003.
\]

Thus no nonzero rank minor vanishes on reduction.  The coordinates of
\(p\) are nonzero modulo the certificate prime as well.

## 3. Principal-part cancellation

If

\[
 \det\operatorname{Hess}(\psi_0+h_6)
\]

is constant, its degree-sixteen term forces

\[
 \det\operatorname{Hess}(h_6)=0.
\]

Five modular leading-Hessian evaluations reduce the 1,725,838 isolated
collision points to 748 candidates.  Full determinant evaluation rejects
727 at the first spatial point and the remaining 21 at the second.  A
rational constant-determinant point would reduce to zero at every one of
these evaluations, so this excludes every isolated characteristic-zero
point.

All 7,566 collision lines are then reconstructed over \(\mathbb Q\).
Exact rational gcds of the principal evaluation polynomials leave only two
lines, with monomial supports

\[
\begin{aligned}
 &(0,2,2,2),(0,3,3,0),(1,2,2,1),(2,2,2,0),\\
 &(2,0,2,2),(2,1,1,2),(2,2,0,2),(3,0,0,3).
\end{aligned}
\]

The first family has the factor \((yr)^2\); its principal gcd is
\(t+81/8\).  The second has the factor \((xs)^2\); its principal gcd is
\(t-9/2\).  Exact symbolic coefficient extraction gives determinant
degrees \(16,12,8,4\), and the degree-twelve coefficient gcd is already one
on both principal roots.

The unique collision plane is supported on

\[
 (yr)^3,\quad xs(yr)^2,\quad (xs)^2yr,\quad (xs)^3.
\]

It is the binary-cubic radial plane in \(yr\) and \(xs\).  The complete
degree-sixteen spatial-coefficient ideal has Gröbner basis \(1\) over
\(\mathbb Q\), so it has no principal survivor.

Combining the isolated, line, and plane calculations proves the stated
support-at-most-four theorem.

## 4. Mixed quartic--sextic theorem

The following is theorem `HC4MQS6`.  It is not part of the sextic-only
statement `HC4MS6`.

Among the 234 quartic principal parts from `HC4MQ1`, exactly four have zero
determinant-degree-two signature.  This condition is immutable when a
sextic is added without a cubic.  For each of those four quartics, enumerate
sextic supports of size at most four subject to

\[
 \nabla h_6(p)=0.
\]

The homogeneous exponent kernel has 976 three-support lines, 205,494
four-support lines, and 519 four-support planes.  Principal Hessian
cancellation leaves 121,146 lines.

For each quartic, 52,686 lines have an affine unit gcd over
\(\mathbb F_{1000003}\) together with a nonzero cubic scale coefficient.
The affine and projective parameter charts are therefore both empty, which
promotes these lines to characteristic zero.  The remaining 68,460
lower-scale-degree lines have unit exact rational gcds at the same two
spatial evaluations.

For the planes, the full determinant at three exact spatial points gives a
three-generator ideal in two sextic parameters.  Direct affine-matrix
determinant expansion and Singular Gröbner bases prove that all
\(519\cdot4=2076\) ideals are unit ideals over \(\mathbb Q\).

Consequently no zero-gradient sextic supported on at most four monomials
repairs any of the 234 quartic principal parts from `HC4MQ1`.

This does not classify the general mixed quartic--sextic chart.  A sextic
changes determinant degrees sixteen through four and can in principle
rescue a quartic that failed the old sextic-free degree-eight principal
screen.  Such quartics are outside the set of 234 treated here.
The later theorem `HC4JQS4` excludes all such mixed cases of combined
support at most four; the unrestricted statement in this paragraph
remains open only beyond that bound.

Replay the complete characteristic-zero certificate with:

```bash
.venv/bin/python scripts/verify_hc4_meng_mixed_quartic_sextic.py
```

## 5. Joint total-support-four theorem

The following is theorem `HC4JQS4`.  Unlike `HC4MQS6`, its genuinely mixed
part does not start from the 234 quartics surviving the old sextic-free
principal screen.

Write

\[
 h=h_4+h_6
\]

with at least one quartic and at least one sextic monomial.  Since both
degrees are even, the collision equation at \(p\) and \(-p\) is the single
four-row system

\[
 \sum_{|e|=4}d_e e+\sum_{|e|=6}d_e e
 =
 \left(-\frac{81}{2},18,18,-\frac{81}{2}\right),
 \qquad d_e=c_ep^e.
\]

The checker exhausts all 6,133,820 genuinely mixed supports of total size
two, three, or four.  Their collision census is

\[
\begin{array}{c|r|r|r|r}
\text{total support}&\text{inconsistent}&\text{isolated}
&\text{lines}&\text{planes}\\ \hline
2&2930&10&0&0\\
3&165514&6446&30&0\\
4&695358&5219228&44270&34
\end{array}
\]

The same Hadamard bounds used above show that reduction modulo \(1000003\)
preserves every collision rank.  At each spatial evaluation the checker
interpolates

\[
 \det\!\left(H_0+zH_4+z^2H_6\right)-64
\]

at \(z=0,\ldots,8\), and inspects its coefficients from \(z^8\) down to
\(z\).  Thus the spatial determinant layers are tested in descending
degrees sixteen through two.  The rejection counts for the 5,225,684
isolated solutions are

\[
\begin{array}{c|rrrrrrrr}
\text{determinant degree}&16&14&12&10&8&6&4&2\\ \hline
\text{rejected}&3839383&655804&641880&42559&44656&12&1390&0.
\end{array}
\]

For every positive-dimensional collision family, the support system is
reconstructed over \(\mathbb Q\).  Exact evaluation polynomials give unit
gcds for all 44,300 lines: 517 are rejected by the first spatial
evaluation and 43,783 by the second.  Exact bivariate Gröbner bases give
the unit ideal for all 34 planes after the third spatial evaluation.
Therefore no genuinely mixed correction of combined support at most four
has constant Hessian determinant.

The pure quartic boundary is `HC4MQ1`, and the pure sextic boundary is
`HC4MS6`.  Adding those two already-proved cases gives the combined
total-support-at-most-four statement above.  Replay the new genuinely
mixed certificate with

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_joint_quartic_sextic_total_support_four.py
```

## 6. Dense common-kernel theorem

Let

\[
 \psi=q_2+h_4+h_6
\]

be an even four-variable potential over a characteristic-zero field, with
\(\operatorname{Hess}(q_2)=H_0\) nonsingular, \(h_4\) homogeneous quartic,
and \(h_6\) homogeneous sextic.  The following is theorem `HC4DCK`:

> Suppose there is a nonzero constant direction \(v\) with
> \(D_vh_6=0\) and \(D_v^2h_4=0\).  Then
> \(\det\operatorname{Hess}(\psi)\) cannot be a nonzero constant while
> \(\nabla\psi\) has an antipodal collision.  In particular, the conclusion
> holds whenever \(\operatorname{Hess}(h_6)\) has generic rank three.

There is no support restriction on either homogeneous layer.

### 6.1 The sextic kernel and the quartic

Under \(w\mapsto\lambda w\), put \(z=\lambda^2\).  Constancy of the
determinant gives

\[
 \det(H_0+zH_4+z^2H_6)=\det(H_0).
\]

The coefficient of \(z^8\) says
\(\det\operatorname{Hess}(h_6)=0\).  By the characteristic-zero
Gordan--Noether theorem in four variables, \(h_6\) has a constant Hessian
kernel direction.  Choose coordinates \((t,u_1,u_2,u_3)\) along a direction
which also satisfies \(D_t^2h_4=0\).  This is the common-kernel hypothesis.
If the sextic Hessian has generic rank three, its kernel is one-dimensional
and the coefficient of \(z^7\) is

\[
 \det\operatorname{Hess}_u(h_6)\,\partial_t^2h_4.
\]

so the common-kernel hypothesis is automatic.  In every case covered by
the theorem,

\[
 h_4=t\,a_3(u)+b_4(u).
\]

Write

\[
 q_2=\frac{\kappa}{2}t^2+t\ell(u)+q(u).
\]

If \(\kappa\ne0\), the critical equation in \(t\) has a polynomial
solution.  Its Schur complement is a three-variable constant-Hessian
potential and carries the antipodal collision.  This contradicts
\(\mathrm{HC}_3\).

### 6.2 The isotropic bordered case

It remains to take \(\kappa=0\) and write

\[
 \psi=t\,s(u)+\phi(u),\qquad s=\ell+a_3.
\]

Put \(C=\operatorname{Hess}(a_3)\) and
\(v=\ell+\nabla a_3\).  The coefficient of \(t^2\) in the Hessian
determinant is

\[
 -v^{\mathsf T}\operatorname{adj}(C)v.
\]

Its degree-six part is

\[
 -\frac32a_3\det C.
\]

Thus either \(a_3=0\), or \(\det C=0\).  Gordan--Noether in three
variables gives a constant kernel direction \(m\) for \(a_3\).  The
degree-two part
\(\ell^{\mathsf T}\operatorname{adj}(C)\ell=0\) lets \(m\) be chosen in
\(\ker\ell\).  Hence, after a linear change in \(u\),

\[
 s=x+a(x,y)
\]

is independent of \(m\).

Let

\[
 R=(\nabla s)^{\mathsf T}
   \operatorname{adj}(\operatorname{Hess}_{x,y}a)\nabla s.
\]

The coefficient of \(t\) factors exactly as

\[
 -\phi_{mm}R.
\]

There are two cases.

If \(\phi_{mm}=0\), then

\[
 \phi=m\,g(x,y)+h(x,y)
\]

and direct block expansion gives

\[
 \det\operatorname{Hess}(\psi)
 =
 \operatorname{Jac}(s,g)^2.
\]

Here \(s\) and \(g\) are odd polynomials of degrees at most three and five.
The degree-at-most-100 plane Jacobian theorem makes \((s,g)\) a polynomial
automorphism.  An antipodal critical pair would give a nonzero common zero
of \((s,g)\), which is impossible.

If \(R=0\), write

\[
 a=\alpha x^3+\beta x^2y+\gamma xy^2+\delta y^3.
\]

The coefficients of \(x\), \(y\), and then \(x^3\) in \(R\) give

\[
 2\gamma=0,\qquad6\delta=0,\qquad-4\beta^2=0.
\]

Thus \(a=\alpha x^3\).  If \(\alpha\ne0\), the vector
\(\nabla s=(1+3\alpha x^2,0)\) vanishes over the algebraic closure, forcing
the full Hessian determinant to vanish there.  Hence \(\alpha=0\) and
\(s=x\).  Expansion along the \(t,x\) hyperbolic block gives

\[
 \det\operatorname{Hess}(\psi)
 =
 -\det\operatorname{Hess}_{y,m}(\phi).
\]

At \(x=0\), \(\mathrm{HC}_2\) makes the remaining two-variable gradient an
automorphism, so the only critical point is the origin.  This also excludes
the antipodal collision.

The exact determinant identities and the three-coefficient binary
elimination are checked by:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_rank_three_sextic_reduction.py
```

The external structural inputs are the
[Gordan--Noether classification](https://arxiv.org/abs/1501.05168),
the known Hessian conjecture in dimensions at most three, and
[Moh's plane degree bound](https://www.math.purdue.edu/~ttm/jacobian.pdf).

## 7. Closure of the even chart and what remains

The later theorem `HC4E46`, proved in
[`HC4_SOURCE_DUAL_BIGRADING.md`](HC4_SOURCE_DUAL_BIGRADING.md), closes the
entire support-free chart

\[
 \psi=q_2+h_4+h_6
\]

with no cubic correction.  Thus `HC4JQS4`, `HC4MQS6`, and the dense
rank-three theorem above are retained as exact checkpoints and proof
components, not as the boundary of the current result.

The pure homogeneous quartic chart, including arbitrary nondegenerate
quadratic renormalizations, is `HC4HQ1` in
[`HC4_MENG_SPARSE_QUARTIC_AUDIT.md`](HC4_MENG_SPARSE_QUARTIC_AUDIT.md).
The full dense cubic--quartic chart is `HC4CQ1` in
[`HC4_MENG_DENSE_CUBIC_QUARTIC.md`](HC4_MENG_DENSE_CUBIC_QUARTIC.md).

The rank-at-most-two sextic locus first becomes a binary or ternary
vanishing-Hessian problem over a polynomial base.  The source/dual
reorganization identifies the dual-linear part with \(JC(2)\), and
weighted Hessian faces synchronize the remaining binary kernels to one
rational projective cone.  Bihomogeneity then leaves only
\(c(X^{\mathsf T}MU)^2\) as a moving rank-two cone.  If \(M\) has rank
one, `HC4DCK` applies; if it has rank two, the dual-degree-four part of the
spatial \(z^4\) face is

\[
 48c^4\det(M)^2(X^{\mathsf T}MU)^4,
\]

so it cannot cancel.  The one-base rank-one boundary synchronizes to a
constant direction by scaling and is also `HC4DCK`.  Rank zero is
`HC4HQ1`.

The earlier obstruction model

\[
 h_4=(xt+ym)^2
\]

has no constant binary kernel direction, but its full Hessian contributes
\(48(xt+ym)^4\).  It is exactly the moving cone killed by the new
determinant face, not a constant-Hessian collision candidate.

The coordinate-chart problem with simultaneous cubic and sextic
interaction and homogeneous support \(\{2,3,4,6\}\) is now closed.  The
rank-three part of
\(q_2+h_3+h_4+h_6\) is `HC4T31` in
[`HC4_MENG_TRIPLE_RANK_THREE.md`](HC4_MENG_TRIPLE_RANK_THREE.md), and
the rank-two part is `HC4T21` in
[`HC4_MENG_TRIPLE_RANK_TWO.md`](HC4_MENG_TRIPLE_RANK_TWO.md).  The
rank-one part and the exhaustion theorem `HC4TC1` are in
[`HC4_MENG_TRIPLE_RANK_ONE.md`](HC4_MENG_TRIPLE_RANK_ONE.md); rank zero
is `HC4CQ1`.  The pure quartic chart is `HC4HQ1`, and the even
quartic--sextic chart is `HC4E46`.  Quintic and higher homogeneous layers,
and non-coordinate coisotropic embeddings, remain separate routes.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_sextic_obstruction.py
```

The exact replay takes about one minute on the reference development
machine and requires only the pinned Python environment.
