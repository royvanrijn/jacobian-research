# The absolute \(D_5\) affine-modification frontier

## 1. Outcome

Work over a characteristic-zero field \(k\) containing \(\sqrt5\).  For

\[
 P=P_5(a,u)=a^5-5a^3u+5au^2,\qquad
 J=J_5(a,u)=\partial_aP,
\]

put

\[
 \alpha=\frac{3+\sqrt5}{2},\qquad
 \beta=\frac{3-\sqrt5}{2},
\]

and

\[
 C=a^2-4u,\qquad
 R_+=a^2-\alpha u,\qquad
 R_-=a^2-\beta u.
\]

Then

\[
 \frac J5=R_+R_-,
 \qquad
 P^2-4u^5=CR_+^2R_-^2.                         \tag{1.1}
\]

Thus the determinant-one derivative-unit chart is

\[
 (a,u,z)\longmapsto
 \left(u,P,-\frac zJ\right).                       \tag{1.2}
\]

This note carries out the requested precomputation before a coefficient
search.  It gives three conclusions.

1. The two ramification colors \(R_+\) and \(R_-\) are individually easy:
   filling either one by a Danielewski equation gives affine three-space.
   Filling both at once produces a factorial fourfold, but it has an
   isolated quadratic-cone singularity and is therefore not affine
   four-space.
2. The one-equation fill \(xy=R_+R_-\) is singular and has class group
   \(\mathbb Z\).  Adding a free coordinate preserves both defects.  Hence
   neither the total-product fill nor the separated two-color fill meets
   the affine-space requirement.
3. There is an all-degree obstruction to the first genuinely coupled
   coefficient ansatz.  If \(u\) is retained and all other outputs are
   affine-linear in any number of new coordinates, with arbitrary
   polynomial \(a,u\)-dependent masks, constant Jacobian forces the target
   field to recover \(a\).  The generic degree is then one, not five.

No absolute Keller map is constructed here.  The surviving search must
both modify the \(u\)-coordinate and contain nonlinear dependence on the
new coordinates.  Constructing such a map would give a degree-five
counterexample to the Jacobian conjecture, so the remaining gate should
not be treated as a routine pole clearing.

Over \(\mathbb Q\), the displayed Dickson row has discriminant square class
\(5\), so its arithmetic group is not the natural
\(D_5\subset A_5\).  The conclusions below concern the geometric
split-field cover.  An arithmetic \(D_5\) construction must either descend
from the cosine field with its monodromy checked or begin with a different
rational \(D_5\) polynomial.

## 2. Boundary, units, and divisor lattice

Let

\[
 B^\circ=\operatorname{Spec}
 k[u,v,(v^2-4u^5)^{-1}]
\]

and let \(X^\circ\) be its inverse image under
\((a,u)\mapsto(u,P)\).  Equation (1.1) gives

\[
 X^\circ=\operatorname{Spec}
 k[a,u,(CR_+R_-)^{-1}].
\]

Both rings are localizations of polynomial UFDs, so

\[
 \operatorname{Cl}(B^\circ)=
 \operatorname{Cl}(X^\circ)=0.
\]

Modulo constants, their unit lattices are

\[
 \mathcal O(B^\circ)^\times/k^\times
   =\mathbb Z[v^2-4u^5],
\]

\[
 \mathcal O(X^\circ)^\times/k^\times
   =\mathbb Z[C]\oplus\mathbb Z[R_+]\oplus\mathbb Z[R_-].
\]

In the ordered source basis \(([C],[R_+],[R_-])\),

\[
 \pi^*[v^2-4u^5]=(1,2,2),\qquad
 [J]=(0,1,1).                                      \tag{2.1}
\]

The primitive derivative vector is not an integral multiple of the sole
target boundary vector.  In particular, a target-only monomial Cox
correction cannot supply \(J^{-1}\): the unramified color \(C\) is the
parity obstruction.  A successful modification must distinguish the two
ramification colors on the source and nevertheless descend to the single
target branch divisor.

The curves \(R_+=0\) and \(R_-=0\) are smooth, but they meet at
\((a,u)=(0,0)\) with intersection multiplicity two.  Indeed, substituting
\(u=a^2/\alpha\) into \(R_-\) gives

\[
 R_-\equiv
 \left(1-\frac{\beta}{\alpha}\right)a^2
 \pmod{R_+}.                                      \tag{2.2}
\]

This tangency is exactly where the naive simultaneous Cox fills become
singular.

## 3. The product Danielewski fill

The smallest total-product fill is

\[
 T_{\rm prod}
 =\{xy=R_+R_-\}\subset\mathbb A^4_{a,u,x,y}.        \tag{3.1}
\]

It has one singular point, the origin.  The partial derivatives of its
defining equation vanish there, and solving

\[
 x=y=\partial_a(R_+R_-)=\partial_u(R_+R_-)=0
\]

gives only \(a=u=0\).

The ring of (3.1) is normal: it is a hypersurface and its singular locus
has codimension three.  Localizing at \(x\) gives a UFD.  The two height-one
primes above \(x=0\) are

\[
 \mathfrak p_\pm=(x,R_\pm),
\]

and

\[
 \operatorname{div}(x)=\mathfrak p_++\mathfrak p_-.
\]

Nagata's theorem therefore gives

\[
 \operatorname{Cl}(T_{\rm prod})
 \simeq
 \frac{\mathbb Z[\mathfrak p_+]\oplus
       \mathbb Z[\mathfrak p_-]}
      {\mathbb Z(\mathfrak p_++\mathfrak p_-)}
 \simeq\mathbb Z.                                  \tag{3.2}
\]

The positive grading

\[
 \deg a=\deg x=\deg y=1,\qquad \deg u=2
\]

shows

\[
 \mathcal O(T_{\rm prod})^\times=k^\times.          \tag{3.3}
\]

Consequently \(T_{\rm prod}\) is neither factorial nor smooth, hence is
not affine three-space.  Polynomial stabilization does not repair this:
for every \(m\ge0\), \(T_{\rm prod}\times\mathbb A^m\) remains singular,
has the same unit group, and has class group \(\mathbb Z\).

## 4. Separated two-color Cox fill

Filling one color alone is harmless:

\[
 T_\pm=\{x_\pm y_\pm=R_\pm\}.
\]

Since \(R_\pm\) is linear in \(u\), eliminating \(u\) gives

\[
 T_\pm\simeq\mathbb A^3_{a,x_\pm,y_\pm}.            \tag{4.1}
\]

Thus each single-color fill is smooth, factorial, and has only constant
units.  It cannot clear \(J\), because the other ramification valuation
remains with coefficient one in (2.1).

The separated simultaneous fill is

\[
 T_{\rm sep}=
 \left\{
 \begin{aligned}
 x_+y_+&=R_+,\\
 x_-y_-&=R_-
 \end{aligned}
 \right\}.                                         \tag{4.2}
\]

Eliminate \(u=(a^2-x_+y_+)/\alpha\).  Its coordinate ring becomes

\[
 k[a,x_+,y_+,x_-,y_-]/
 \left(
 x_-y_-
 -\frac{\beta}{\alpha}x_+y_+
 -\left(1-\frac{\beta}{\alpha}\right)a^2
 \right).                                          \tag{4.3}
\]

The quadratic form in (4.3) is nondegenerate.  Hence (4.2) is the
four-dimensional affine quadric cone with unique singular point at the
origin.  Its projectivization is a smooth quadric threefold, whose Picard
group is generated by the hyperplane class.  The standard cone class-group
sequence gives

\[
 \operatorname{Cl}(T_{\rm sep})=0,\qquad
 \mathcal O(T_{\rm sep})^\times=k^\times.           \tag{4.4}
\]

So the separated Cox ring is factorial, but factoriality is not enough:
the vertex excludes \(T_{\rm sep}\simeq\mathbb A^4\).  The failure is caused
by the tangency (2.2), not by an unaccounted divisor class.

Adjoining the unramified color \(C\) in the same separated-product manner
also fails before any coefficient search.  At the common origin, the
three relation gradients for \(C,R_+,R_-\) are all proportional to \(du\);
the resulting complete intersection is singular.

## 5. Affine-linear mask theorem

The preceding Cox models suggest feeding new variables into the cover
outputs.  The entire affine-linear version of that idea can be excluded.

### Theorem 5.1

Let \(K\) be a characteristic-zero field, let \(r\ge1\), and put
\(R=K[a]\).  Take

\[
 h\in R^{r+1},\qquad
 M\in\operatorname{Mat}_{r+1,r}(R),
\]

and define

\[
 F(a,z)=h(a)+M(a)z,\qquad z=(z_1,\ldots,z_r)^T.     \tag{5.1}
\]

If

\[
 \det\frac{\partial F}{\partial(a,z_1,\ldots,z_r)}
 \in K^\times,                                     \tag{5.2}
\]

then

\[
 K(F_0,\ldots,F_r)=K(a,z_1,\ldots,z_r).             \tag{5.3}
\]

In particular, (5.1) has generic degree one.

### Proof

Condition (5.2) implies that \(M\) has rank \(r\).  Let
\(p=(p_0,\ldots,p_r)\) be its signed maximal-minor row, so \(pM=0\).
Expansion along the \(a\)-derivative column gives

\[
 \det DF
 =p\,h'+\sum_{j=1}^r z_j\,p\,M'_j.                 \tag{5.4}
\]

Since (5.4) is a nonzero constant,

\[
 pM'=0,\qquad ph'\in K^\times.                     \tag{5.5}
\]

Write \(p=dq\), where \(d\in R\) is the gcd of the entries and \(q\) is
primitive.  Then \(qM=qM'=0\).  Differentiating \(qM=0\) gives \(q'M=0\).
The left kernel of \(M\) over \(K(a)\) is one-dimensional, so \(q'\) is
a rational multiple of \(q\).  All projective ratios of the entries of
\(q\) are therefore constant in \(a\).  Primitivity forces

\[
 q\in K^{r+1}.                                     \tag{5.6}
\]

Now \(d(qh')\in K^\times\).  Both factors lie in \(K[a]\), hence
\(d\in K^\times\) and

\[
 qh=\lambda a+\mu,\qquad
 \lambda\in K^\times,\quad\mu\in K.                \tag{5.7}
\]

Because \(qM=0\), equations (5.1) and (5.7) recover

\[
 a=\frac{qF-\mu}{\lambda}.
\]

After recovering \(a\), any nonzero maximal minor of \(M\) recovers
\(z\) rationally.  This proves (5.3). \(\square\)

### Corollary 5.2

Consider a polynomial map

\[
 (a,u,z_1,\ldots,z_r)
 \longmapsto
 \bigl(u,F_0,\ldots,F_r\bigr),                     \tag{5.8}
\]

where the \(F_i\) are affine-linear in the \(z_j\), with arbitrary
coefficients in \(k[a,u]\).  If (5.8) has nonzero constant Jacobian, then
its generic degree is one.

Indeed, apply Theorem 5.1 over \(K=k(u)\).  Therefore no two- or
three-new-coordinate modification of the form (5.8) can preserve the
degree-five \(D_5\) function-field extension.  This is stronger than the
block-triangular obstruction: all the affine-linear masks may depend on
\((a,u)\) and may enter every output except the retained \(u\).

## 6. Fibre and Jacobian gates for the surviving search

The precomputation leaves a sharply smaller search space.

1. **Both ramification colors must participate.**  A one-color fill leaves
   the other coefficient of \([J]=(0,1,1)\) uncancelled.
2. **The target branch alone is insufficient.**  Its pullback vector
   \((1,2,2)\) cannot produce \((0,1,1)\) integrally.
3. **The standard product and separated Cox spaces are not affine
   spaces.**  The former is nonfactorial and singular; the latter is
   factorial but singular.
4. **The \(u\)-coordinate must feed back.**  Retaining \(u\) and using
   affine-linear auxiliary masks forces degree one by Corollary 5.2.
5. **The auxiliary dependence must be nonlinear.**  Merely increasing the
   number of affine-linear masks does not evade the theorem.
6. **Degree and complete fibres must be checked by elimination.**  On
   \(v^2-4u^5\ne0\), the final auxiliary system must contribute exactly one
   point above each of the five roots of \(P(A,u)-v\), with no localization
   and no extra auxiliary degree.  Polynomiality alone does not guarantee
   this.

A minimal surviving ansatz therefore has the schematic form

\[
\begin{aligned}
 U&=u+\text{terms nonlinear in }z_1,z_2,\\
 V&=P(a,u)+\text{coupled terms},\\
 Z_1&=\text{coupled terms},\\
 Z_2&=\text{coupled terms},
\end{aligned}                                      \tag{6.1}
\]

with at least one genuinely nonlinear auxiliary term and with an explicit
elimination certificate recovering the original quintic incidence field.
The vertex calculation shows where its linear part must differ from the
naive Cox relations: the conormal rows at the common tangent point must be
made independent without deleting or merging either ramification color.

This is the first coefficient search worth running.  A search that fixes
\(U=u\), uses only affine-linear masks, or starts from a freely stabilized
product fill is already closed by the results above.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_absolute_dihedral_d5_modification_frontier.py
```

The checker verifies the \(D_5\) factorization and divisor vectors, the
order-two intersection of the ramification colors, the unique singular
vertices of the product and separated fills, the nondegenerate separated
quadric, the singular conormal rank of the three-color fill, and the
maximal-minor determinant ledger behind Theorem 5.1 for two and three new
coordinates.
