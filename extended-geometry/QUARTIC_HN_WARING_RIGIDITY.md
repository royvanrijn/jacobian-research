# Quartic HN Waring rigidity and the first viable ordinary-Laplacian architecture

**Research note — 4 August 2026**

## Status

This note records a theorem-sized reduction for the search for a
six-variable homogeneous quartic Hessian-nilpotent polynomial whose gradient
map is noninjective.

The principal conclusions are:

> **There is no essential six-variable quartic HN polynomial of Waring rank
> at most eight. Every essential rank-nine quartic HN polynomial has a
> polynomially invertible gradient map. Consequently, every essential
> six-variable quartic HN counterexample has Waring rank at least ten.**

The proof has four layers:

1. coordinate-free trace and Gram identities;
2. an exact codimension-two Gale calculation for rank eight;
3. a complete one-latitude classification for the rank-nine stratum with one
   relation-free Waring term; and
4. a complementary-minor/Cauchy--Binet obstruction for the fully supported
   rank-nine Gale matroid.

The finite calculations are over `QQ` and `QQ(lambda)`, not bounded
coefficient searches. The rank-nine top-determinant step uses the standard
matroid base-packing theorem together with finite exact checks of the two
exceptional multiplicity patterns.

I have not located these exact Waring-rank statements in the HN literature
checked so far. That is not a priority or novelty claim; expert literature
review is still required.

---

## 1. Waring--Gram form of a quartic Hessian

Work over an algebraically closed characteristic-zero field. Let

\[
 P(z)=\sum_{i=1}^r c_i\ell_i(z)^4
\]

be a minimal Waring decomposition, with every \(c_i\ne0\). Absorb fourth
roots of the \(c_i\) into the linear forms. Thus we may write

\[
 P(z)=\sum_{i=1}^r\ell_i(z)^4.
\]

Let \(v_i=\nabla\ell_i\), put

\[
 V=(v_1\ \cdots\ v_r),\qquad
 G=V^{\mathsf T}V,\qquad
 l=V^{\mathsf T}z,
\]

and let

\[
 D(l)=\operatorname{diag}(l_1^2,\ldots,l_r^2).
\]

Then

\[
 \operatorname{Hess}P=12VD(l)V^{\mathsf T}.
 \tag{1.1}
\]

Consequently,

\[
 \operatorname{tr}(\operatorname{Hess}P)
 =12\sum_i g_{ii}l_i^2,
 \tag{1.2}
\]

and

\[
 \operatorname{tr}((\operatorname{Hess}P)^2)
 =144\left(
 \sum_i g_{ii}^2l_i^4+
 2\sum_{i<j}g_{ij}^2l_i^2l_j^2
 \right).
 \tag{1.3}
\]

If \(P\) is HN, both expressions vanish identically.

The space of possible value vectors \(l\) is

\[
 L=\operatorname{im}V^{\mathsf T}
   =(\ker V)^\perp\subseteq k^r.
 \tag{1.4}
\]

This turns the first two nilpotence traces into low-degree relations on one
linear subspace.

---

## 2. Two one-relation lemmas

Let

\[
 R=k[x_1,\ldots,x_r]/(u_1x_1+\cdots+u_rx_r),
\]

where the support \(S=\{i:u_i\ne0\}\) has at least three elements.

### Lemma 2.1 — diagonal-square rigidity

If

\[
 \sum_i a_ix_i^2=0\quad\text{in }R,
\]

then every \(a_i=0\).

#### Proof

The quadratic is divisible by \(U=\sum u_ix_i\), so it equals
\(U\sum b_ix_i\). Comparing the coefficient of \(x_ix_j\), \(i\ne j\),
gives

\[
 u_ib_j+u_jb_i=0.
\]

If \(i\notin S\), comparison with any \(j\in S\) gives \(b_i=0\). For
three distinct \(i,j,k\in S\), the pair equations give

\[
 \frac{b_i}{u_i}
 =-\frac{b_j}{u_j}
 =\frac{b_k}{u_k}
 =-\frac{b_i}{u_i}.
\]

Characteristic zero forces all these ratios to vanish. Hence \(b=0\), and
so \(a=0\). \(\square\)

### Lemma 2.2 — square-pair rigidity

If

\[
 \sum_{i<j}A_{ij}x_i^2x_j^2=0\quad\text{in }R,
 \tag{2.1}
\]

then every \(A_{ij}=0\).

#### Proof

For three distinct \(i,j,k\in S\), set every other variable to zero and use

\[
 x_k=-\frac{u_ix_i+u_jx_j}{u_k}.
\]

In the restricted quartic, the coefficients of \(x_i^4\) and \(x_j^4\)
first give \(A_{ik}=A_{jk}=0\); the coefficient of \(x_i^2x_j^2\) then gives
\(A_{ij}=0\). Thus every edge internal to \(S\) vanishes.

If \(p,q\notin S\), the vector supported only on \(p,q\) lies in the
hyperplane, so \(A_{pq}=0\).

Finally let \(p\notin S\). For distinct \(i,j\in S\), use the vector supported
on \(p,i,j\), with \(x_j=-(u_i/u_j)x_i\). Since \(A_{ij}=0\),

\[
 \frac{A_{pi}}{u_i^2}+
 \frac{A_{pj}}{u_j^2}=0.
\]

Three indices in \(S\) make these signs inconsistent unless every
\(A_{pi}=0\). \(\square\)

---

## 3. One-relation Waring rigidity

### Theorem 3.1

Let

\[
 P=\sum_{i=1}^r\ell_i^4
\]

be a minimal Waring decomposition. Assume the \(\ell_i\) have at most one
independent linear relation. If \(P\) is HN, then

\[
 \langle v_i,v_j\rangle=0
 \qquad\text{for all }i,j.
 \tag{3.1}
\]

In particular,

\[
 (\operatorname{Hess}P)^2=0,
\]

and \(z-\nabla P\) is the quasi-translation with polynomial inverse
\(z+\nabla P\).

#### Proof

If the \(\ell_i\) are independent, (1.2) and (1.3) give (3.1)
coefficientwise.

Otherwise the unique relation in a minimal Waring decomposition has support
at least three: support one would give a zero form, and support two would
make two forms proportional. Apply Lemma 2.1 to (1.2), then Lemma 2.2 to
(1.3). Hence \(G=0\).

Equation (1.1) now gives

\[
 (\operatorname{Hess}P)^2
 =144VDG\,DV^{\mathsf T}=0.
\]

Euler's identity for a quartic gives

\[
 (\operatorname{Hess}P)z=3\nabla P.
\]

Therefore

\[
 (\operatorname{Hess}P)\nabla P=0.
\]

The derivation \(\nabla P\cdot\nabla\) fixes \(\nabla P\), so its exponential
is \(z\mapsto z+t\nabla P\). The inverse assertion follows. \(\square\)

### Corollary 3.2

An essential quartic HN polynomial in \(n\) variables which is not a
quasi-translation has Waring rank at least \(n+2\).

The next section removes equality when \(n=6\).

---

## 4. The rank-eight trace dichotomy

Assume now that \(P\) is essential in six variables and has a minimal
eight-term decomposition. Then

\[
 \operatorname{rank}V=6,\qquad \dim\ker V=2.
\]

Let

\[
 a_i=g_{ii}.
\]

Equation (1.2) says that the diagonal quadratic
\(\sum a_ix_i^2\) vanishes on the six-plane \(L\subset k^8\).

Let \(T=\{i:a_i\ne0\}\), \(m=|T|\), and let \(s_T\) be the dimension of the
span of the coordinate forms indexed by \(T\). Since the quadratic on
\(k^T\) is nondegenerate and its restriction to that span is zero,

\[
 s_T\le\left\lfloor\frac m2\right\rfloor.
 \tag{4.1}
\]

There are only two global relations, so

\[
 s_T\ge m-2.
 \tag{4.2}
\]

Thus \(m\le4\). Minimality excludes \(m=1,2,3\): the corresponding forms
would span at most \(0,1,1\) dimensions as appropriate, forcing a zero or
proportional Waring term. Hence

\[
 m=0\quad\text{or}\quad m=4.
 \tag{4.3}
\]

If \(m=4\), then equality holds throughout (4.1)--(4.2): the four forms in
\(T\) span a two-plane and the complete relation space is supported on
\(T\). The other four forms are free coordinates on \(L\). Restricting
(1.3) to those four free coordinates forces their four coefficient vectors
to be pairwise orthogonal and isotropic. They are linearly independent, so
they would span a four-dimensional totally isotropic subspace of a
nondegenerate six-dimensional quadratic space. The Witt index is three.
Contradiction.

Therefore

\[
 \boxed{g_{ii}=0\quad(1\le i\le8).}
 \tag{4.4}
\]

Every Waring vector in a rank-eight candidate is isotropic.

Equation (1.3) reduces to

\[
 \sum_{i<j}g_{ij}^2l_i^2l_j^2=0\quad\text{on }L.
 \tag{4.5}
\]

---

## 5. Codimension-two square-pair classification

Choose a \(2\times8\) Gale matrix \(K\) whose row space is the relation space
of the \(\ell_i\). Its columns consist of:

- \(z\) zero columns;
- nonzero columns in \(m\) projective directions, with multiplicities
  \(\nu_1\ge\cdots\ge\nu_m\).

A nonzero relation among the Waring forms has minimum support

\[
 (8-z)-\nu_1.
 \tag{5.1}
\]

Minimality requires this number to be at least three.

### 5.1 Five or more projective Gale directions

Suppose \(m\ge5\).

Edges joining two coordinates in the same projective Gale class vanish by
evaluating (4.5) on the two-coordinate vector in \(\ker K\).

Choose one representative from five distinct classes. For a triple
\(i,j,k\), let

\[
 d_{ij}=\det(K_i,K_j).
\]

The circuit vector supported on the triple gives, after dividing by the
nonzero determinant squares,

\[
 B_{ij}+B_{ik}+B_{jk}=0,
 \qquad
 B_{ij}=\frac{A_{ij}}{d_{ij}^2}.
 \tag{5.2}
\]

On four vertices these equations say that opposite edges are equal. A fifth
vertex makes every edge equal; a triangle then gives \(3B=0\). Thus all
edges between nonzero Gale classes vanish.

If \(p\) is a zero Gale column, the coefficient of \(x_p^2\) is a diagonal
quadratic on the nonzero-column kernel. Any such quadratic has support at
most four by the same isotropic-dimension argument as in Section 4. If two
distinct Gale directions remain outside its support, no relation is
supported there and the coordinate forms are independent. If exactly one
direction remains, the supported forms have at most one relation, and
Lemma 2.1 applies. Hence this diagonal quadratic is zero. Edges joining two
zero columns vanish directly.

Therefore (4.5) has only the zero coefficient vector when \(m\ge5\).

### 5.2 At most four projective directions

It remains to consider

\[
 z\le3,\qquad 2\le m\le4,\qquad
 (8-z)-\nu_1\ge3.
 \tag{5.3}
\]

There are exactly 25 multiplicity patterns.

Normalize the directions to

\[
 (0,\infty),\qquad
 (0,1,\infty),\qquad
 (0,1,\infty,\lambda)
\]

for \(m=2,3,4\), respectively. For each pattern, parameterize
\(L=\ker K\) by six variables and form the coefficient matrix

\[
 \mathcal M_{z,\nu}:
 k^{28}\longrightarrow\operatorname{Sym}^4(k^6),
 \qquad
 (A_{ij})\longmapsto
 \sum_{i<j}A_{ij}l_i^2l_j^2.
 \tag{5.4}
\]

The checker exhibits one exact \(28\times28\) minor in every case.

| \(z\) | multiplicities \(\nu\) | certified determinant |
|---:|---|---:|
| 0 | (5,3) | \(-4096\) |
| 0 | (5,2,1) | \(4096\) |
| 0 | (4,4) | \(4096\) |
| 0 | (4,3,1) | \(-4096\) |
| 0 | (4,2,2) | \(-8192\) |
| 0 | (3,3,2) | \(-8192\) |
| 1 | (4,3) | \(-4096\) |
| 1 | (4,2,1) | \(-4096\) |
| 1 | (3,3,1) | \(-4096\) |
| 1 | (3,2,2) | \(-8192\) |
| 2 | (3,3) | \(-4096\) |
| 2 | (3,2,1) | \(-4096\) |
| 2 | (2,2,2) | \(-8192\) |
| 3 | (2,2,1) | \(4096\) |
| 0 | (5,1,1,1) | \(4096\lambda^6\) |
| 0 | (4,2,1,1) | \(-8192\lambda^4\) |
| 0 | (3,3,1,1) | \(-8192\lambda^4\) |
| 0 | (3,2,2,1) | \(8192\) |
| 0 | (2,2,2,2) | \(-8192\) |
| 1 | (4,1,1,1) | \(-4096\lambda^6\) |
| 1 | (3,2,1,1) | \(-8192\lambda^4\) |
| 1 | (2,2,2,1) | \(-8192\) |
| 2 | (3,1,1,1) | \(-4096\lambda^6\) |
| 2 | (2,2,1,1) | \(-8192\lambda^4\) |
| 3 | (2,1,1,1) | \(4096\lambda^6\) |

The four-direction normalization requires \(\lambda\ne0,1\), so every listed
minor is nonzero. Thus (5.4) is injective in all 25 cases.

The only remaining Gale pattern has four zero columns and four distinct
nonzero directions. Its square-pair kernel has dimension five. But the four
zero-column Waring vectors are free, independent, isotropic, and pairwise
orthogonal, again contradicting Witt index three.

---

## 6. Main consequence

### Theorem 6.1 — rank-eight Waring obstruction

There is no essential six-variable quartic HN polynomial of Waring rank at
most eight.

#### Proof

Ranks at most seven have relation nullity at most one, so Theorem 3.1 applies.

At rank eight, Section 4 first makes every Waring vector isotropic. Sections
5.1--5.2 show that the square-pair identity (4.5) has no nonzero coefficient
vector compatible with an essential minimal decomposition. Hence
\(g_{ij}=0\) for every \(i,j\), so the six-dimensional span of the \(v_i\)
would be totally isotropic. Contradiction. \(\square\)

### Corollary 6.2

Any six-variable homogeneous quartic HN counterexample to the
ordinary-Laplacian vanishing conjecture has

\[
 \boxed{\operatorname{WaringRank}(P)\ge9.}
\]

Since any minimum-dimensional counterexample is essential, this is a direct
search lower bound for the six-variable moonshot.

This is not merely a bounded coefficient search: every coefficient and every
linear form in Waring ranks at most eight is covered.

---

## 7. A second exact obstruction: three relation-free terms at rank nine

Rank nine is the first possible Waring rank. Its relation nullity is three.

A natural first architecture leaves three Waring terms outside every linear
relation. In Gale language, the \(3\times9\) relation matrix has three zero
columns and the remaining six columns have rank three.

This complete architecture is impossible, even when the six-term internal
Gram block is arbitrary.

### Proposition 7.1 — no \(3+6\) relation-support split

Let

\[
 P=e_1^4+e_2^4+e_3^4+\sum_{i=1}^6 b_i^4
 \tag{7.1}
\]

be a minimal essential rank-nine Waring decomposition, where the three
\(e_p\) occur in no linear relation and the six \(b_i\) span a
three-dimensional complementary space. Then \(P\) is not HN.

#### Proof

The value coordinates \(e_1,e_2,e_3\) are free on the Waring-value
six-plane. Restricting the first two trace identities to the locus
\(b_1=\cdots=b_6=0\) shows that their coefficient vectors are isotropic and
pairwise orthogonal. Hence they span a maximal isotropic three-plane \(E\).

Let \(W\) be the span of the six remaining coefficient vectors. Essentiality
and the relation count give

\[
 k^6=E\oplus W,\qquad \dim E=\dim W=3.
\]

The pairing \(E\times W\to k\) is nondegenerate: if \(w\in W\) is orthogonal
to \(E\), then \(w\in E^\perp=E\), hence \(w=0\). Choose bases so that the
cross-pairing is the identity. Write the six block vectors as columns
\(c_i\in k^3\). Their Waring values are \(c_i^{\mathsf T}t\), where
\(t=(t_1,t_2,t_3)\) are independent coordinates on \(W\), and

\[
 \langle e_p,c_i\rangle=(c_i)_p.
\]

The coefficient of the free monomial \(e_p^2\) in
\(\operatorname{tr}((\operatorname{Hess}P)^2)\) is

\[
 288\sum_{i=1}^6(c_i)_p^2(c_i^{\mathsf T}t)^2.
 \tag{7.2}
\]

Put

\[
 F_W(t)=\sum_{i=1}^6(c_i^{\mathsf T}t)^4.
\]

Equation (7.2) is \(\frac1{12}\partial_{t_p}^2F_W\). HN therefore forces

\[
 \partial_{t_1}^2F_W=
 \partial_{t_2}^2F_W=
 \partial_{t_3}^2F_W=0.
\]

A quartic in three variables annihilated by all three pure second
derivatives would have exponent at most one in every variable, but no
degree-four monomial has that property. Thus \(F_W=0\).

The six block fourth powers therefore cancel identically, contradicting
minimality of the nine-term Waring decomposition. \(\square\)

### Corollary 7.2

A rank-nine six-variable HN candidate has at most two zero Gale columns.
Equivalently, at least seven of its nine Waring terms must participate in
the three-dimensional relation space. Proposition 7.3 below removes the
two-zero-column case as well.

If exactly two Gale columns are zero, let \(e_1,e_2\) be their free
values and let \(l_1,\ldots,l_7\) be the active values. Write \(A\) for the
\(2\times7\) free-to-active Gram block and put

\[
 R=A\operatorname{diag}(l_1^2,\ldots,l_7^2)A^{\mathsf T}.
\]

The second trace gives \(R_{11}=R_{22}=0\). In the fourth trace, the
coefficient of \(e_1^2e_2^2\) is a nonzero scalar multiple of
\(R_{12}^2\): the two quadratic-in-\(e\) contributions are the top-left and
bottom-right alternating blocks, and every term using the internal active
Gram block has lower \(e\)-degree. Hence \(R_{12}=0\) as well.

Equivalently, if \(t_1,t_2\) are the active coordinates dual to the two free
isotropic vectors, then

\[
 \partial_{t_1}^2F_W=
 \partial_{t_1}\partial_{t_2}F_W=
 \partial_{t_2}^2F_W=0.
\]

The seven-term active quartic therefore has the sharper exact form

\[
 \boxed{
 F_W=
 t_1C_3(t_3,t_4)
 +t_2D_3(t_3,t_4)
 +E_4(t_3,t_4).
 }
 \tag{7.3}
\]

Thus the first relation-support stratum not removed by the present argument
is not merely separately linear in two latitude variables: it is jointly
affine-linear in the entire two-plane.


In fact this remaining two-latitude stratum is also empty.

### Proposition 7.3 — no two relation-free terms at rank nine

A minimal essential rank-nine HN Waring decomposition cannot have two zero
Gale columns.

#### Proof

Keep the notation above. Let \(c_1,\ldots,c_7\in k^4\) be the coefficient
vectors of the seven active linear forms, and let

\[
 \Phi:k^7\longrightarrow\operatorname{Sym}^2(k^4),
 \qquad
 (r_i)\longmapsto\sum_i r_ic_ic_i^{\mathsf T}.
\]

Let \(a,b\in k^7\) be the two rows of the free-to-active Gram block \(A\).
The second and fourth trace calculations just made say

\[
 a^{\circ2},\quad a\circ b,\quad b^{\circ2}\in\ker\Phi.
 \tag{7.4}
\]

The \(c_i\) span four dimensions and are pairwise nonproportional. Choose
four independent ones. Their symmetric squares are independent; in that
basis they are four distinct diagonal matrix units. Any fifth
nonproportional vector has at least two nonzero coordinates, so its square
has an off-diagonal entry and leaves that four-space. Therefore

\[
 \operatorname{rank}\Phi\ge5,\qquad \dim\ker\Phi\le2.
 \tag{7.5}
\]

The rows \(a,b\) are independent: otherwise one nonzero vector in the free
two-plane would be orthogonal to the free plane and to the active span,
contradicting nondegeneracy of the ambient pairing. Their coordinatewise
symmetric square has dimension at least two. By (7.4)--(7.5), it has
dimension exactly two.

Hence a nonzero binary quadratic vanishes on every nonzero projective column

\[
 [a_i:b_i]\in\mathbb P^1.
\]

Over the algebraic closure it has at most two roots. Since \(A\) has rank
two, both roots occur. Change the basis of the row space of \(A\) so that
the two roots are the coordinate points. Because the whole coordinatewise
symmetric square of the row space lies in \(\ker\Phi\), the transformed
rows \(a',b'\) still satisfy

\[
 \sum_{i\in I}(a'_i)^2c_ic_i^{\mathsf T}=0,\qquad
 \sum_{j\in J}(b'_j)^2c_jc_j^{\mathsf T}=0,
 \tag{7.6}
\]

where \(I,J\) are disjoint, every displayed coefficient is nonzero, and
both sets are nonempty. Columns of \(A\) equal to zero belong to neither
set.

For a relation

\[
 \sum_{i\in I}\alpha_ic_ic_i^{\mathsf T}=0
 \qquad(\alpha_i\ne0),
\]

the span of the \(c_i\), \(i\in I\), is totally isotropic for a
nondegenerate diagonal form on \(k^I\). Its dimension is therefore at most
\(\lfloor |I|/2\rfloor\). Pairwise nonproportionality rules out
\(|I|=1,2,3\). Thus \(|I|\ge4\), and similarly \(|J|\ge4\).

But \(I\) and \(J\) are disjoint subsets of seven indices. Contradiction.
\(\square\)

### Corollary 7.4 — first live rank-nine stratum

Every rank-nine six-variable HN Waring candidate has at most one zero Gale
column.

With exactly one zero Gale column, the eight active terms span five
dimensions. If \(t\) denotes the active coordinate dual to the free
isotropic vector, the second trace forces their active quartic to satisfy

\[
 \partial_t^2F_W=0,
\]

hence

\[
 \boxed{F_W=t\,C_3(x_1,x_2,x_3,x_4)+E_4(x_1,x_2,x_3,x_4).}
 \tag{7.7}
\]

This one-latitude, eight-channel form is the first Waring relation-support
architecture not excluded by the relation-count argument alone.

### Proposition 7.5 — canonical one-latitude block

Assume exactly one Gale column is zero. Let its Waring covector be denoted
by \(s\), and let \(W\) be the five-dimensional span of the other eight
covectors. The first trace makes \(s\) isotropic. The functional
\(w\mapsto\langle s,w\rangle\) on \(W\) is nonzero, since otherwise \(s\)
would be in the radical of the ambient quadratic form.

Choose \(t\in W\) with \(\langle s,t\rangle=1\), and put
\(U=W\cap s^\perp\). The induced form on \(U\) is nondegenerate. After
replacing \(t\) by \(t+u\) for a suitable \(u\in U\), the Gram matrix on
\(ks\oplus kt\oplus U\) is

\[
 J=\begin{pmatrix}
 0&1&0\\
 1&\beta&0\\
 0&0&Q
 \end{pmatrix},
 \tag{7.8}
\]

where \(Q\) is a nondegenerate symmetric \(4\)-by-\(4\) matrix. The second
trace says that the active eight-term quartic is affine-linear in the
\(t\)-coordinate. Hence

\[
 \boxed{P=s^4+tC(x)+E(x)},
 \qquad x=(x_1,x_2,x_3,x_4),
 \tag{7.9}
\]

with \(C\) cubic and \(E\) quartic.

Put

\[
 u=\nabla C,\qquad v=Qu,
 \qquad B=Q\bigl(t\operatorname{Hess}C+\operatorname{Hess}E\bigr),
 \qquad a=12s^2.
\]

The endomorphism corresponding to the Hessian of \(P\) is

\[
 A=J\operatorname{Hess}P
 =\begin{pmatrix}
 0&0&u^{\mathsf T}\\
 a&0&\beta u^{\mathsf T}\\
 0&v&B
 \end{pmatrix}.
 \tag{7.10}
\]

A direct Schur-complement calculation gives

\[
 \boxed{
 \det(\lambda I-A)
 =\lambda^2\det(\lambda I-B)
 -(\beta\lambda+a)
 u^{\mathsf T}\operatorname{adj}(\lambda I-B)v.
 }
 \tag{7.11}
\]

Since \(a=12s^2\) is an independent coefficient, \(A\) is nilpotent if and
only if

\[
 \det(\lambda I-B)=\lambda^4
 \tag{7.12}
\]

and

\[
 u^{\mathsf T}B^jv=0,
 \qquad 0\le j\le3.
 \tag{7.13}
\]

Thus the first rank-nine problem is a four-dimensional nilpotent Hessian
pencil plus four explicit transport moments. The parameter \(\beta\) drops
out completely.

### Proposition 7.6 — the cubic profile is binary isotropic

The \(j=0\) equation in (7.13) is

\[
 u^{\mathsf T}Qu=0.
 \tag{7.14}
\]

Consequently the gradient image of \(C\) lies in the isotropic quadric of
\(Q\). Differentiating (7.14) makes the Hessian of \(C\) singular. The
classical Hesse theorem in at most four variables, applied successively to
the essential variable space of \(C\), shows that \(C\) depends on at most
two linear forms. If it depends essentially on two, its gradient map is
dominant on their two-plane, so (7.14) makes that plane totally isotropic.
The one-variable case is the same statement with an isotropic line.

After an orthogonal change in the \(x\)-space, therefore,

\[
 \boxed{C=C(x_1,x_2)},
 \qquad
 Q=\begin{pmatrix}0&I_2\\I_2&0\end{pmatrix},
 \tag{7.15}
\]

where \(x_1,x_2\) span a maximal isotropic plane and the remaining variables
are denoted \(y_1,y_2\). Up to the \(\operatorname{GL}_2\)-action preserving
this split pairing, a nonzero binary cubic has exactly three root types:

\[
 x_1^3,
 \qquad x_1^2x_2,
 \qquad x_1x_2(x_1-x_2).
 \tag{7.16}
\]

### Theorem 7.7 — nontriple cubic profiles are polynomial automorphisms

For each of the last two normal forms in (7.16), put

\[
 M=Q\operatorname{Hess}E,
 \qquad N=Q\operatorname{Hess}C.
\]

The three linear necessary conditions

\[
 \operatorname{tr}M=0,
 \qquad
 \operatorname{tr}(MN)=0,
 \qquad
 u^{\mathsf T}Mv=0
 \tag{7.17}
\]

have an exact ten-dimensional solution space. The checker computes this
space over \(\mathbb Q\); in invariant notation it is precisely

\[
 \boxed{
 E=h_4(x_1,x_2)+y_1H_1(x_1,x_2)+y_2H_2(x_1,x_2),
 \qquad
 \partial_{x_1}H_1+\partial_{x_2}H_2=0.
 }
 \tag{7.18}
\]

No term quadratic in \(y_1,y_2\) survives. Write \(H=(H_1,H_2)\). For
(7.18), the four-dimensional pencil in (7.12) is block lower triangular:

\[
 B=\begin{pmatrix}
 JH&0\\
 *&(JH)^{\mathsf T}
 \end{pmatrix}.
 \tag{7.19}
\]

All transport moments (7.13) vanish automatically. Hence the six-variable
potential is HN exactly when \(JH\) is nilpotent.

The divergence equation in (7.18) writes

\[
 H=(\partial_{x_2}\psi,-\partial_{x_1}\psi)
\]

for a binary quartic \(\psi\). Nilpotence of \(JH\) is equivalent to
\(\det\operatorname{Hess}\psi=0\), so the binary Hesse theorem gives
\(\psi=\gamma L^4\). Therefore

\[
 H=\kappa(b,-a)(ax_1+bx_2)^3,
 \qquad (JH)^2=0.
 \tag{7.20}
\]

The associated Keller map can be inverted triangularly. In the coordinates
of (7.8), it is

\[
\begin{aligned}
 s'&=s-C(x),\\
 t'&=t-4s^3-\beta C(x),\\
 x'&=x-H(x),\\
 y'&=(I-(JH)^{\mathsf T})y-t\nabla C(x)-\nabla h_4(x).
\end{aligned}
 \tag{7.21}
\]

First recover \(x=x'+H(x')\), then \(s\), then \(t\), and finally

\[
 y=(I+(JH)^{\mathsf T})
 \bigl(y'+t\nabla C(x)+\nabla h_4(x)\bigr).
 \tag{7.22}
\]

Thus neither the double-root nor the squarefree cubic profile can yield an
ordinary-Laplacian counterexample.

### Proposition 7.8 — the triple-root branch splits into two channels

For \(C=x_1^3\), the linear system (7.17) has dimension sixteen. It writes

\[
 E=y_1A_3(x_1,x_2,y_2)+B_4(x_1,x_2,y_2),
 \qquad
 \partial_{x_1}A_3+
 \partial_{x_2}\partial_{y_2}B_4=0.
 \tag{7.23}
\]

The cubic coefficient is

\[
\begin{aligned}
 A_3={}&p_{13}x_1^3
 +(p_2x_2^3+p_6x_1x_2^2+p_{10}x_1^2x_2)\\
 &+(p_1y_2^3+p_5x_1y_2^2+p_9x_1^2y_2).
\end{aligned}
 \tag{7.24}
\]

In the checker coordinates, the coefficient of \(t\) in
\(\operatorname{tr}(B^3)\) is exactly the edge ideal of \(K_{3,3}\):

\[
 (p_1,p_5,p_9)(p_2,p_6,p_{10})=0.
 \tag{7.25}
\]

At every field-valued point, therefore, one of the two triples vanishes. The
first removes the \(y_2\)-dependent part of \(A_3\), while the second removes
the \(x_2\)-dependent part. Swapping \(x_2\) and \(y_2\) preserves the split
quadratic form and exchanges the two branches. Hence the only remaining
one-zero-Gale search may be normalized to

\[
 \boxed{C=x_1^3,
 \qquad A_3=A_3(x_1,x_2).}
 \tag{7.26}
\]

The remaining trace equations in fact close this last one-latitude branch.

### Theorem 7.9 — no relation-free Waring term at rank nine

A minimal essential rank-nine HN Waring decomposition which has exactly one
zero Gale column always gives a polynomial automorphism.

#### Proof

The nontriple cubic profiles are Theorem 7.7. It remains to use (7.26).
Keep the checker parameters of (7.24), and use the branch

\[
 p_1=p_5=p_9=0,
 \tag{7.27}
\]

the other branch being obtained by \(x_2\leftrightarrow y_2\).

If \(p_0\ne0\), exact coefficients of the second trace give

\[
 p_2=p_6=p_{10}=p_7=p_{11}=p_3=p_{13}=0.
 \tag{7.28}
\]

Thus \(A_3=0\), and the remaining polynomial is triangular in the paired
variables.

If \(p_3\ne0\), the coefficient identities

\[
 p_2p_4=0,
 \qquad p_2p_8+3p_3p_4=0,
 \qquad 8p_0p_3+p_2p_4=0,
 \qquad p_2^2p_8=0
 \tag{7.29}
\]

first give \(p_0=p_4=0\). If \(p_2=0\), then
\(p_6^2=3p_{10}p_2=0\), while

\[
 p_{10}p_8=0,
 \\qquad
 p_{10}^2+3p_{12}p_2-2p_{13}p_6+4p_3p_8=0
\]

forces \(p_8=0\). If \(p_2\ne0\), (7.29) already gives \(p_8=0\).
Hence \(E\) is linear in both \(y_1,y_2\), and the argument of Theorem 7.7
applies.

We may therefore assume

\[
 p_0=p_3=0.
 \tag{7.30}
\]

The exact monomial trace coefficients include

\[
 p_2p_4=p_2p_8=p_6p_4=p_6p_8
 =p_{10}^2p_4=p_{10}^2p_8=0.
 \tag{7.31}
\]

If one of \(p_2,p_6,p_{10}\) is nonzero, then \(p_4=p_8=0\); again \(E\)
is linear in \(y_1,y_2\) and Theorem 7.7 applies.

In the last case

\[
 p_2=p_6=p_{10}=0.
\]

The remaining necessary equations include

\[
 2p_{11}p_8+9p_{13}^2=0
 \tag{7.32}
\]

and

\[
 8p_{11}^2p_8^2+108p_{11}p_{13}^2p_8+81p_{13}^4=0.
 \tag{7.33}
\]

Substituting (7.32) into (7.33) gives
\(-243p_{13}^4=0\), hence \(p_{13}=0\) and \(A_3=0\). Moreover

\[
 p_{11}p_8=p_7p_8=p_{11}p_4=p_4p_7=0.
 \tag{7.34}
\]

These equations split \(B_4\) into a polynomial in \(y_2\) plus a term
linear in \(x_2\), or a polynomial in \(x_2\) plus a term linear in
\(y_2\). The \((x_2,y_2)\)-map is triangular. Since \(x_1\) is fixed when
\(A_3=0\), the complete six-variable map is triangular as well.

Thus every field-valued HN point in the one-zero-Gale stratum is a
polynomial automorphism. \(\square\)

### Corollary 7.10 — the first live rank-nine architecture

Every essential rank-nine quartic HN counterexample must have **no zero Gale
columns**. All nine Waring terms must participate in the three-dimensional
relation space.

### Proposition 7.11 — every fully supported rank-nine covector is isotropic

Let \(T=\{i:g_{ii}\ne0\}\), put \(m=|T|\), and let \(s_T\) be the dimension
of the span of those Waring covectors. The first trace gives, as in Section
4,

\[
 m-3\le s_T\le\left\lfloor\frac m2\right\rfloor.
 \tag{7.35}
\]

Minimality removes \(m\le3\). Hence \(m\in\{4,5,6\}\).

For \(m=5\), one has \(s_T=2\), so the three-dimensional relation space is
entirely supported on \(T\). The four complementary Gale columns are zero.
For \(m=6\), one similarly has \(s_T=3\), and the three complementary Gale
columns are zero. Both contradict Corollary 7.10.

It remains to take \(m=4\). Then \(s_T=2\), and exactly two independent
relations are supported on \(T\). Choose them as the first two rows of the
Gale matrix. On the four-dimensional source subspace annihilating the
\(T\)-span, the other five Waring values have one relation. Since none of
their Gale columns is zero, that relation has full support. Restricting the
second trace to this subspace and applying Lemma 2.2 gives

\[
 g_{ij}=0
 \qquad(i,j\notin T).
\]

Those five complementary covectors are isotropic and pairwise orthogonal.
Their images modulo the \(T\)-span already have dimension four, so their
actual span is a totally isotropic space of dimension at least four. This
contradicts the Witt index three.

Therefore \(T=\varnothing\):

\[
 \boxed{g_{11}=\cdots=g_{99}=0.}
 \tag{7.36}
\]

Thus the sole live rank-nine architecture is a fully supported configuration
of nine isotropic Waring covectors.

## 8. The rank-nine top-determinant obstruction

The fully supported rank-nine architecture of Corollary 7.10 is also
impossible. The decisive invariant is the determinant of the Hessian at one
specially chosen source point, rather than another closed-walk trace.

### Lemma 8.1 — cyclic-complement determinant gate

Let \(V\) be an \(n\times(n+c)\) matrix of rank \(n\), and let \(K\) be a
\(c\times(n+c)\) Gale matrix: its row space is \(\ker V\). Put

\[
 H(l)=V\operatorname{diag}(l_1^2,\ldots,l_{n+c}^2)V^{\mathsf T},
 \qquad l\in\ker K.
 \tag{8.1}
\]

Let \(T\) be a basis of the column matroid \(M(K)\), and put
\(S=E\setminus T\), so \(|S|=n\). If the restriction \(M(K)|S\) has no
coloops, then there is an \(l\in\ker K\) with

\[
 \operatorname{supp}l=S
 \tag{8.2}
\]

and

\[
 \boxed{
 \det H(l)=\det(V_S)^2\prod_{i\in S}l_i^2\ne0.
 }
 \tag{8.3}
\]

#### Proof

A finite set is cyclic precisely when its matroid restriction has no
coloops. Thus every coordinate of \(\ker K_S\) is nonzero on some vector.
Over an infinite field, the union of the finitely many coordinate
hyperplanes cannot cover \(\ker K_S\); choose a full-support vector and
extend it by zero on \(T\). Since

\[
 \operatorname{im}V^{\mathsf T}=\ker K,
\]

this vector is the Waring-value vector at some source point.

Complementary maximal minors of a Gale pair vanish simultaneously, so
\(K_T\) invertible implies \(V_S\) invertible. Cauchy--Binet gives

\[
 \det H(l)=
 \sum_{|I|=n}\det(V_I)^2\prod_{i\in I}l_i^2.
\]

Only \(I=S\) survives, proving (8.3). \(\square\)

For an HN quartic, every Hessian specialization is nilpotent and therefore
singular. Hence any Gale matroid admitting such a basis \(T\) is excluded
immediately.

### Theorem 8.2 — no fully supported rank-nine architecture

Let \(P\) be an essential six-variable quartic with a minimal rank-nine
Waring decomposition. If its \(3\times9\) Gale matrix has no zero column,
then \(P\) is not HN.

#### Proof

Let \(M\) be the rank-three column matroid of the Gale matrix. Minimality of
the Waring decomposition says that every nonzero linear relation has support
at least three. Projectively, every line therefore contains at most six of
the nine columns, counted with multiplicity. In particular, a projective
point has multiplicity at most five.

Let \(m\) be the maximum projective-point multiplicity.

**Case 1: \(m\le3\).** For every flat \(X\) of \(M\),

\[
 |X|\le3\operatorname{rk}X:
\]

the rank-one case is the multiplicity bound, the rank-two case is the
six-on-a-line bound, and the rank-three case is \(|E|=9\). The matroid
base-packing theorem therefore partitions \(E\) into three disjoint bases.
Take one as \(T\); the complement is the union of two bases and has no
coloops. Lemma 8.1 applies.

**Case 2: \(m=4\).** Let four copies lie at the projective point \(p\).
Each line through \(p\) contains at most two of the five outside columns.
Their direction-multiplicity partition is therefore one of

\[
 (2,2,1),\qquad(2,1,1,1),\qquad(1,1,1,1,1).
\]

The five outside columns can be split into a three-set and a two-set such
that each set uses distinct directions from \(p\). Let \(T\) consist of one
copy of \(p\) and the two-set; let \(S\) contain the other three copies and
the three-set. Then \(T\) is a basis. Removing any element of \(S\) still
leaves a spanning set: the repeated copies handle an element at \(p\), and
\(p\) together with either two remaining outside directions handles an
outside element. Thus \(S\) has no coloops, and Lemma 8.1 applies.

**Case 3: \(m=5\).** Five columns lie at \(p\), while the four outside
columns lie on distinct lines through \(p\). Proposition 7.11 makes all nine
Waring covectors isotropic. Restrict the source to the subspace on which the
four outside Waring values vanish. The remaining five values have exactly
one full-support relation: the restrictions of the Gale rows to the five
repeated columns span the line generated by \((1,1,1,1,1)\), after scaling.
The restricted value space therefore has dimension four.

The second trace, together with Lemma 2.2, forces the five corresponding
Waring covectors to be pairwise orthogonal. Their restrictions span four
dimensions, so their actual span has dimension at least four. This is a
four-dimensional totally isotropic subspace of a nondegenerate
six-dimensional quadratic space, contradicting Witt index three.

All cases are impossible. \(\square\)

### Theorem 8.3 — rank-nine HN maps are automorphisms

Let \(P\) be an essential six-variable quartic HN polynomial of Waring rank
nine. Then the map

\[
 z\longmapsto z-\nabla P(z)
\]

is a polynomial automorphism.

#### Proof

Three or two zero Gale columns are excluded by Propositions 7.1 and 7.3.
Theorem 7.9 proves polynomial invertibility when there is exactly one zero
Gale column. Theorem 8.2 excludes the fully supported case. \(\square\)

### Corollary 8.4 — counterexample Waring rank starts at ten

Every essential six-variable homogeneous quartic HN counterexample satisfies

\[
 \boxed{\operatorname{WaringRank}(P)\ge10.}
 \tag{8.4}
\]

This is a counterexample lower bound, not a statement that every HN quartic
has Waring rank at least ten: the rank-nine one-latitude strata may contain
HN polynomials, but their gradient maps are polynomial automorphisms.

---

## 9. Rank-ten parallel-class obstructions

Lemma 8.1 removes every rank-ten Gale matroid having a basis with cyclic
six-element complement.  Exact characteristic-zero examples show that this
matroid condition can fail when projective Gale points have high
multiplicity.  The HN traces eliminate the first two exceptional mechanisms
completely.

### Lemma 9.1 — a large Gale parallel class is totally isotropic

Let \(R\) be a projective parallel class of \(q\ge3\) columns in the
rank-four Gale matrix of an essential rank-ten Waring decomposition.  Then
all Waring covectors indexed by \(R\) are isotropic and pairwise orthogonal.
Consequently

\[
 q\le4.
 \tag{9.1}
\]

#### Proof

After scaling the columns of \(R\), the value vectors supported on \(R\)
form the full-support hyperplane

\[
 l_1+\cdots+l_q=0.
\]

Restrict the first two Hessian traces to this value subspace.  Lemmas 2.1 and
2.2 force

\[
 g_{ii}=g_{ij}=0\qquad(i,j\in R).
\]

The restricted Waring values span \(q-1\) dimensions, so the actual Waring
covectors indexed by \(R\) span at least \(q-1\) dimensions.  Their span is
totally isotropic, while the Witt index of the ambient nondegenerate
six-dimensional quadratic space is three.  Thus \(q-1\le3\). \(\square\)

### Proposition 9.2 — multiplicity four is impossible

A rank-ten six-variable essential quartic HN Waring decomposition has no
four-element projective Gale parallel class.

#### Proof

Let \(R\) be such a class and let \(O\) be its six-element complement.
Lemma 9.1 shows that the four Waring covectors in \(R\) span exactly a
three-dimensional totally isotropic space \(E\).

Gale duality gives

\[
 \operatorname{rk}_V(O)
 =|O|-\operatorname{rk}K+
   \operatorname{rk}_K(R)
 =6-4+1=3.
 \tag{9.2}
\]

Let \(W\) be this three-dimensional outside span.  Essentiality gives

\[
 k^6=E\oplus W.
\]

After one invertible source-coordinate change, the quartic splits as

\[
 P=P_E(x_1,x_2,x_3)+P_W(x_4,x_5,x_6),
 \tag{9.3}
\]

where both summands are essential ternary quartics: \(P_E\) has four Waring
terms spanning \(E\), and \(P_W\) has six terms spanning \(W\).  Hence the
Hessian determinant factors, up to a nonzero square from the coordinate
change,

\[
 \det\operatorname{Hess}P
 =c\,
   \det\operatorname{Hess}P_E\,
   \det\operatorname{Hess}P_W.
 \tag{9.4}
\]

The classical Hesse theorem is valid for ternary forms in characteristic
zero: a homogeneous ternary form with identically zero Hessian determinant
is a cone.  Both summands in (9.3) are essential, so both determinant factors
are nonzero polynomials.  Equation (9.4) contradicts Hessian nilpotence.
\(\square\)

### Proposition 9.3 — two triple classes are impossible

A rank-ten six-variable essential quartic HN Waring decomposition cannot
have two distinct three-element projective Gale parallel classes.

#### Proof

Let the classes be \(A\) and \(B\).  Lemma 9.1 makes each internal Waring
span totally isotropic.  Restrict the Waring-value space by setting the other
four values to zero.  Since the two Gale directions are distinct, the
remaining value space is the product of two full-support three-term
hyperplanes.

The cross part of the second trace is

\[
 \sum_{i\in A,\,j\in B}
 g_{ij}^2l_i^2m_j^2=0.
 \tag{9.5}
\]

The three squares on either full-support hyperplane are linearly independent
by Lemma 2.1.  Their nine tensor products are therefore independent, and
(9.5) forces every cross pairing to vanish.  The six Waring covectors in
\(A\cup B\) consequently span a totally isotropic space of dimension at most
three.

Let \(C\) be the four-element complement.  The dual rank formula gives

\[
 \operatorname{rk}_V(A\cup B)
 =6-4+\operatorname{rk}_K(C)
 =2+\operatorname{rk}_K(C).
 \tag{9.6}
\]

Thus \(\operatorname{rk}_K(C)\le1\).  But the two distinct class directions
span rank two, so

\[
 \operatorname{rk}K
 \le \operatorname{rk}_K(A\cup B)+\operatorname{rk}_K(C)
 \le2+1<4,
\]

contradicting that the Gale matrix has rank four. \(\square\)

### Corollary 9.4 — the sole remaining rank-ten matroid gate

In any rank-ten HN candidate:

1. every projective Gale class has multiplicity at most three; and
2. at most one class has multiplicity three.

The remaining purely combinatorial target is now:

> **Rank-ten cyclic-complement lemma.**  Let \(M\) be a characteristic-zero
> representable rank-four matroid on ten elements.  Assume every cocircuit
> has size at least three, every parallel class has size at most two except
> possibly one class of size three.  Then \(M\) has a basis whose
> six-element complement is cyclic.

This lemma would combine with Lemma 8.1 and Propositions 9.2--9.3 to prove:

\[
 \boxed{
 \text{every essential six-variable quartic HN counterexample has
 Waring rank at least eleven.}
 }
 \tag{9.7}
\]

Exact reconnaissance supports the lemma:

- every sampled rational configuration under its hypotheses has a
  cyclic-complement basis;
- the characteristic-zero configurations without such a basis found by the
  search have either one class of multiplicity four or two classes of
  multiplicity three;
- binary and ternary finite-field exceptions exist, so a proof cannot simply
  be an unrestricted abstract-matroid statement.

This is the current dead end.  It is a small representable-matroid theorem,
not a remaining polynomial-system computation.  The HN side of every
observed exceptional architecture is already closed above.

---

## 10. The exact search representation after the matroid gate

If the cyclic-complement lemma is proved, the first possible Waring rank is
eleven.  Without it, a rank-ten search should only inspect Gale matroids
satisfying Corollary 9.4 and failing the determinant gate.  For each such
\(4\times10\) Gale matrix \(K\), impose

\[
 T_j(l)=\operatorname{tr}
 \left(
   \bigl(\operatorname{diag}(l_1^2,\ldots,l_{10}^2)G\bigr)^j
 \right)=0,
 \qquad l\in\ker K,
 \quad 1\le j\le6,
 \tag{10.1}
\]

on symmetric Gram matrices \(G\) of rank six with \(\ker G\) equal to the
Gale row space.  The nonradial collision is tested only after the trace
scheme survives.

Thus neither route returns to the generic 126-dimensional quartic
coefficient space.

---

## 11. Reproduction

Run

```bash
.venv/bin/python scripts/verify_quartic_hn_waring_rigidity.py
.venv/bin/python scripts/verify_quartic_hn_rank9_one_zero.py
.venv/bin/python scripts/verify_quartic_hn_rank9_top_determinant.py
.venv/bin/python scripts/verify_quartic_hn_rank10_parallel_obstructions.py
```

Expected final markers:

```text
QUARTIC_HN_WARING_PATTERN_COUNT=25
QUARTIC_HN_WARING_EXCEPTIONAL_KERNEL_DIMENSION=5
QUARTIC_HN_WARING_RANK8_GATE_PASS
QUARTIC_HN_WARING_RANK9_LINEAR_PASS profile=triple rank=19 dimension=16
QUARTIC_HN_WARING_RANK9_LINEAR_PASS profile=double rank=25 dimension=10
QUARTIC_HN_WARING_RANK9_LINEAR_PASS profile=squarefree rank=25 dimension=10
QUARTIC_HN_WARING_RANK9_TRIPLE_K33_PASS
QUARTIC_HN_WARING_RANK9_LINEAR_GATE_PASS
QUARTIC_HN_WARING_RANK9_ONE_ZERO_GATE_PASS
QUARTIC_HN_RANK9_COMPLEMENTARY_MINOR_PASS
QUARTIC_HN_RANK9_ONE_TERM_CAUCHY_BINET_PASS
QUARTIC_HN_WARING_RANK9_TOP_DETERMINANT_GATE_PASS
QUARTIC_HN_RANK10_CLASS4_BLOCK_FACTOR_PASS
QUARTIC_HN_RANK10_TWO_TRIPLE_TENSOR_RANK
QUARTIC_HN_WARING_RANK10_PARALLEL_GATE_PASS
```

The first script prints every exact minor determinant and verifies the
rank-nine linear normal forms. The second independently replays the finite
trace coefficients used in Theorem 7.9.

---

## 12. Scope and next step

This note does **not** construct an ordinary-Laplacian counterexample or
prove that six variables are impossible. It establishes the following
complete structural results:

1. every Waring rank at most eight is impossible for an essential quartic HN
   polynomial;
2. every essential rank-nine quartic HN polynomial has a polynomially
   invertible gradient map;
3. therefore every essential six-variable quartic HN counterexample has
   Waring rank at least ten.

The next theorem-directed problem is precisely the rank-ten
cyclic-complement lemma in Corollary 9.4.  The determinant gate removes the
generic Gale matroid, while the HN traces remove every observed
characteristic-zero obstruction family.  I do not currently have a clean
proof of the residual matroid statement, and the small-characteristic
counterexamples show that representability in characteristic zero must be
used essentially.  This is therefore the current mathematical dead end,
rather than a computational timeout or an unexamined coefficient boundary.
