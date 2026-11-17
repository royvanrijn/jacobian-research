# Degree-four equivariant Keller maps in three variables

## 1. The classification

Let \(k\) be a characteristic-zero field.  In displayed affine coordinates,
let a nontrivial linear \(\mathbb G_m\)-action on \(\mathbb A^3\) have
integer weights

\[
 w=(w_1,w_2,w_3)\ne(0,0,0).
\]

The source and target actions may initially be written separately.  The
invertible linear part of a Keller map identifies their weight multisets.
After a target permutation, equivariant linear normalization, and (in the
zero-weight space) a target translation, the problem is therefore

\[
 F_i=x_i+H_i,\qquad
 \operatorname{wt}(H_i)=w_i,\qquad
 2\le \deg H_i\le4,\qquad
 \det JF=1.                                           \tag{1.1}
\]

The computation below proves:

> **Theorem.** Every degree-\(\le4\) Keller map
> \(F:\mathbb A^3_k\to\mathbb A^3_k\) equivariant for a nontrivial linear
> \(\mathbb G_m\)-action is a polynomial automorphism.  It is tame.

Thus the bounded equivariant search produces no coordinate-degree-four
counterexample.  Any degree-four counterexample in three variables must
have no nontrivial linear torus symmetry in its degree-four coordinates
(or acquire such a symmetry only after a nonlinear change that does not
preserve the degree bound).

This is a classification for a **linear action in the displayed
coordinates**.  Although \(\mathbb G_m\)-actions on complex affine
three-space are linearizable, a nonlinear linearizing conjugacy can raise
coordinate degree and is not covered by the bounded monomial census.

## 2. Why the weight search is finite

For a nonlinear monomial \(x^\alpha\) in \(H_i\), equivariance is exactly

\[
 (\alpha-e_i)\cdot w=0,\qquad 2\le|\alpha|\le4.       \tag{2.1}
\]

There are 93 labelled relations and 55 distinct relation vectors.  Their
hyperplane arrangement has three kinds of strata.

1. Two independent active relations determine a primitive weight ray.
   Modulo coordinate permutation and simultaneous reversal \(w\mapsto-w\),
   there are 92 distinct rank-two support types.
2. Active relations spanning one line give 11 rank-one support types.
3. With no active relation, \(F\) is linear.

These counts include repeated and zero weights.  They are generated from
the bounded relations, rather than from an arbitrary cutoff on the entries
of \(w\).

The exact census is implemented in
[`classify_degree_four_equivariant.py`](../scripts/classify_degree_four_equivariant.py).

## 3. The dependency-graph cutoff

Put a directed cross-edge \(j\to i\) when an allowed monomial in \(H_i\)
uses \(x_j\), with \(j\ne i\).  Self-dependence is deliberately not drawn.

> **Acyclic lemma.** If this cross-dependency graph is acyclic, every Keller
> point of the support is triangular and tame.

Indeed, topologically order the vertices.  The Jacobian matrix is triangular
apart from possible diagonal self-dependence, so

\[
 1=\det JF=\prod_i\left(1+\frac{\partial H_i}{\partial x_i}\right).
\]

Every factor is a unit of \(k[x_1,x_2,x_3]\), hence constant.  Since \(H_i\)
has no linear term, \(\partial H_i/\partial x_i=0\).  The map is triangular
in the chosen order.

Of the 92 rank-two supports, 36 are acyclic immediately and 56 contain a
cycle.  Of the 11 rank-one supports, eight are acyclic and three contain a
cycle.  The exact coefficient ideals of the latter three have radical equal
to the full nonlinear coefficient ideal, so only the linear map survives.

## 4. Mixed nonzero signatures

After removing zero weights, every cyclic rank-two support is mixed-sign.
There are 51 such support types.  For each one the checker:

1. constructs the complete normalized support;
2. expands \(\det JF-1\) and forms its exact coefficient ideal over
   \(\mathbb Q\);
3. computes the radical and all minimal associated primes in Singular;
4. tests the surviving support graph on every prime.

The 51 radicals have 78 minimal components: 25 signatures have one
component, 25 have two, and the foundational signature has three.  For 44
of the 51 signatures every component has acyclic surviving support and is
therefore tame by the lemma.

Exactly seven signatures retain a support cycle on at least one radical
component:

| primitive weight, up to permutation and sign | nonlinear coefficients | radical components | outcome |
|---|---:|---:|---|
| \((1,-4,-2)\) | 5 | 1 | two shears |
| \((1,-3,-2)\) | 5 | 1 | two shears |
| \((1,-2,-2)\) | 8 | 1 | binary quadratic shear |
| \((1,-2,-1)\) | 9 | 3 | three foundational tame families |
| \((1,-2,1)\) | 11 | 1 | binary cubic shear |
| \((1,-2,3)\) | 5 | 1 | two shears |
| \((1,-1,-1)\) | 8 | 1 | binary quadratic shear |

The first two one-prime radicals give, respectively,

\[
\begin{aligned}
(x,\ y+a z^2,\ z+b x^2(y+a z^2)),\\
(x,\ y+a xz^2,\ z+b x(y+a xz^2)).
\end{aligned}                                       \tag{4.1}
\]

For \((1,-2,3)\) the family is

\[
 (x+a y(z+b x^3),\ y,\ z+b x^3).                    \tag{4.2}
\]

These are displayed compositions of two elementary shears.

For \((1,-2,-2)\) and \((1,-1,-1)\), the repeated-weight pair has the form

\[
 (x,\ (y,z)+s(x)\,Q(y,z)),                           \tag{4.3}
\]

where \(s=x^2\) and \(s=x\), respectively, and the radical relations say
that the binary quadratic vector \(Q\) is invariant under its own
translation:

\[
 Q\bigl((y,z)+sQ(y,z)\bigr)=Q(y,z).
\]

The inverse is obtained by replacing \(+sQ\) with \(-sQ\).  The
\((1,-2,1)\) component has the identical form on the equal-weight
\((x,z)\)-pair, with parameter \(y\) and a binary cubic vector.  The checker
reduces both inverse compositions modulo the exact radical ideals.

The remaining weight \((1,-2,-1)\) is a permutation of the foundational
\((1,-1,-2)\) signature.  Its three-prime radical and three explicit
two-shear inverses are proved in
[`DEGREE_FOUR_FOUNDATIONAL_WEIGHT_EXCLUSION.md`](../extended-geometry/DEGREE_FOUR_FOUNDATIONAL_WEIGHT_EXCLUSION.md).

## 5. Zero weights

There are three qualitatively different parabolic cases.

### Two zero weights

For \(w=(0,0,1)\), equivariance gives

\[
 F=(A(x,y),B(x,y),zC(x,y)),\qquad
 \det JF=C\,\operatorname{Jac}(A,B).
\]

The product is a nonzero constant, so \(C\) is constant and \((A,B)\) is a
plane Keller map of degree at most four.  Moh's degree-\(\le100\) theorem
makes it a plane automorphism, and Jung--van der Kulk makes the lifted
three-variable map tame.

### One zero weight, same-sign pair

For \(w=(0,p,q)\) with \(pq>0\), the zero-weight coordinate is affine in
itself.  If \(p\ne q\), ordering the two nonzero weights makes the other two
coordinates triangular over \(k[x]\).  If \(p=q\), they form a matrix in
\(\mathrm{GL}_2(k[x])\).  Both cases are tame.

### One zero weight, opposite-sign pair

After permutation and scaling, the only degree-four ratios with a nonlinear
invariant are

\[
 (0,1,-q),\qquad q=1,2,3.
\]

Put \(t=y^qz\).  Every normalized equivariant map is

\[
 F=(A(x,t),\,yB(x,t),\,zC(x,t)),\qquad
 T=tB^qC.
\]

A direct change-of-variables calculation gives the quotient Keller equation

\[
 \operatorname{Jac}_{x,t}(A,T)=B^{q-1}.              \tag{5.1}
\]

On \(t=0\), the original determinant identity is

\[
 A_x(x,0)B(x,0)C(x,0)=1.
\]

Thus every pure-\(x\) nonlinear coefficient vanishes.  Substitution into
(5.1) leaves a small exact ideal.  Its radicals give

\[
\begin{array}{c|c}
q&F\\ \hline
1&(x+a\,yz+b\,(yz)^2,\ y,\ z),\\
2&(x+a\,y^2z,\ y,\ z),\\
3&(x+a\,y^3z,\ y,\ z).
\end{array}                                         \tag{5.2}
\]

All are elementary shears.

## 6. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_degree_four_equivariant_classification.py
```

The checker verifies the arrangement counts, all three cyclic rank-one
radicals, the 51 mixed nonzero radical decompositions, the componentwise
graph cutoff, the six new exceptional radicals and inverses, the existing
foundational three-prime certificate, and the three quotient ideals in
(5.2).

The calculation is set-theoretic over the coefficient schemes: replacing a
Keller ideal by its radical loses nilpotent and embedded structure but loses
no characteristic-zero Keller map.  No numerical solving is used.

## 7. Literature boundary

- T. Shaska, *Graded Keller maps and the Jacobian Conjecture*,
  arXiv:2607.20210 (2026), <https://arxiv.org/abs/2607.20210>.
  It proves the all-one-sign theorem in every dimension, the two-variable
  equivariant theorem for every signature, and the quotient-Jacobian
  framework.  It does not perform this bounded three-variable coefficient
  classification.
- T. T. Moh, *On the Jacobian conjecture and the configurations of roots*,
  J. Reine Angew. Math. **340** (1983), 140--212,
  <https://doi.org/10.1515/crll.1983.340.140>.
  Its plane degree bound is used only for the \((0,0,1)\) cylinder case.
- A. van den Essen, *Polynomial Automorphisms and the Jacobian Conjecture*,
  Progress in Mathematics 190, Birkhäuser (2000), for the standard
  normalization and tame plane-automorphism background.

