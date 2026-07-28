# The \(2+2\) source/dual bigrading in \(HC_4\)

## 1. Status and scope

This note proves theorem `HC4SDW`.  It reorganizes the remaining even
Meng--Yang four-variable chart by a source/dual splitting and proves the
weighted-face lemma needed to interpret its rank-two leading cone.  It
establishes two exact structural statements:

1. the part which is at most linear in the two dual variables is precisely
   a cotangent lift of a plane Keller map; and
2. successive vanishing weighted Hessian faces of a rank-one binary
   residual block determine one rational projective cone, rather than
   unrelated null directions.

The second statement is the **successive-cone synchronization lemma**.
It does not by itself say that the synchronized projective cone is
constant or a polynomial coordinate.  The quartic \((xt+ym)^2\) shows
that such a strengthening is false.  Section 8 uses the degree-four
bigrading and the next full determinant face to exclude that moving cone
inside the even quartic--sextic chart.  This proves that chart, not
unrestricted \(HC_4\).

Throughout, the coefficient field has characteristic zero.

## 2. The source/dual ledger

Use source variables \(X=(x,y)\) and dual variables \(U=(t,m)\).  The
quadratic Meng form

\[
 4xs+2yr
\]

becomes \(tx+my\) after setting \(t=4s\) and \(m=2r\).  Thus the even
quartic--sextic chart has the canonical form

\[
 \Psi(X,U)=tx+my+
 \sum_{a+b=4}h_{a,b}(X,U)+
 \sum_{a+b=6}k_{a,b}(X,U),
\tag{2.1}
\]

where the indices record source degree \(a\) and dual degree \(b\).
The two homogeneous ledgers are

\[
\begin{array}{c|rrrrr}
(a,b)&(4,0)&(3,1)&(2,2)&(1,3)&(0,4)\\ \hline
\dim h_{a,b}&5&8&9&8&5
\end{array}
\]

and

\[
\begin{array}{c|rrrrrrr}
(a,b)&(6,0)&(5,1)&(4,2)&(3,3)&(2,4)&(1,5)&(0,6)\\ \hline
\dim k_{a,b}&7&12&15&16&15&12&7 .
\end{array}
\]

They total \(35\) quartic and \(84\) sextic monomials.  This is a
reindexing of the whole chart, not a support restriction.

## 3. The dual-linear locus is \(JC(2)\)

Assume that every nonlinear term in (2.1) has dual degree at most one.
Then, uniquely,

\[
 \Psi=tF(x,y)+mG(x,y)+H(x,y),
\tag{3.1}
\]

where

\[
\begin{aligned}
 F&=x+f_3+f_5,\\
 G&=y+g_3+g_5,
\end{aligned}
\qquad
 H=h_4+h_6
\]

with the subscripts denoting source degree.  In the variable order
\((t,m,x,y)\),

\[
 \operatorname{Hess}\Psi=
 \begin{pmatrix}
 0&J(F,G)\\
 J(F,G)^{\mathsf T}&*
 \end{pmatrix}.
\]

The block determinant identity gives

\[
 \det\operatorname{Hess}\Psi
 =\det J(F,G)^2.
\tag{3.2}
\]

The source-only summand \(H\) does not occur in (3.2).  Hence a nonzero
constant Hessian determinant is equivalent to the plane map
\((F,G)\) being Keller.

There is also an exact collision statement.  A critical point of (3.1)
satisfies \(F=G=0\).  Conversely, at any common zero \(X\), the remaining
two critical equations are

\[
 J(F,G)(X)^{\mathsf T}U+\nabla H(X)=0.
\]

They determine a unique \(U\), because (3.2) makes the plane Jacobian a
nonzero constant.  Since \(F,G\) are odd and \(H\) is even, a nonzero
solution occurs with its negative and gives an antipodal critical pair.
The origin is a simple zero because the linear part of \((F,G)\) is
\((x,y)\).  Thus the antipodal collision obstruction on this stratum is
exactly a second zero of an odd plane Keller map, not a new four-variable
Hessian phenomenon.

In the present chart \(\deg(F,G)\le5\).  Moh's plane theorem through
degree \(100\) excludes such a second zero.  We will refer to (3.1) as
the **dual-linear \(JC(2)\) locus**.  It has \(13\) quartic and \(19\)
sextic nonlinear coefficients before imposing the collision and
determinant equations.

## 4. Weighted Hessian faces

Let \(z=(z_1,\ldots,z_n)\), let
\(\rho=(\rho_1,\ldots,\rho_n)\) be an integral weight, and put

\[
 D_\rho(\lambda)=
 \operatorname{diag}(\lambda^{\rho_1},\ldots,\lambda^{\rho_n}).
\]

If \(p_d\) is a polynomial face of \(\rho\)-weight \(d\), then the
derivative shifts disappear after diagonal conjugation:

\[
 D_\rho(\lambda)
 \operatorname{Hess}(p_d)(\lambda^\rho z)
 D_\rho(\lambda)
 =\lambda^d\operatorname{Hess}(p_d)(z).
\tag{4.1}
\]

Consequently, for \(p=\sum_d p_d\),

\[
 K_\rho(\lambda):=
 D_\rho(\lambda)\operatorname{Hess}(p)(\lambda^\rho z)
 D_\rho(\lambda)
 =\sum_d\lambda^d\operatorname{Hess}(p_d)(z).
\tag{4.2}
\]

This is the correct Hessian face pencil.  Taking the initial form of
each uncorrected Hessian entry separately would introduce different
derivative shifts and is not invariant.

If \(\det\operatorname{Hess}(p)=c\), then

\[
 \det K_\rho(\lambda)
 =c\,\lambda^{2\sum_i\rho_i}.
\tag{4.3}
\]

Thus every determinant face except the single exponent on the right
vanishes.  For the source/dual weight
\(\rho=(0,0,1,1)\), formula (4.2) is simply

\[
 D\operatorname{Hess}\Psi(X,\lambda U)D
 =\sum_b\lambda^b\operatorname{Hess}(\Psi_{\bullet,b})(X,U),
\quad
 D=\operatorname{diag}(1,1,\lambda,\lambda)
\]

when the variables are ordered \((x,y,t,m)\), and (4.3) is
\(\det K_\rho=c\lambda^4\).

## 5. Successive-cone synchronization

The relevant local algebra is only two-dimensional.

### Lemma 5.1 (successive-cone synchronization)

Let \(A\) be a domain, let \(\epsilon\) define a filtration, and let

\[
 S(\epsilon)=
 \begin{pmatrix}
 a(\epsilon)&b(\epsilon)\\
 b(\epsilon)&c(\epsilon)
 \end{pmatrix}
 \in\operatorname{Sym}_2(A[[\epsilon]]).
\]

Suppose the first nonzero face of \(S\) has rank one.  On either principal
chart on which a diagonal pivot of that face is a unit, factor out its
common power of \(\epsilon\) and call the resulting upper-left entry
\(a\).  Then:

1. if \(\det S=0\), there is a unique series \(r=b/a\) such that
   \[
    S=a
    \begin{pmatrix}1&r\\r&r^2\end{pmatrix};
   \tag{5.1}
   \]
2. its kernel is the single projective line generated by
   \((-r,1)\);
3. if only the first \(N+1\) determinant faces vanish, the same
   conclusions hold modulo \(\epsilon^{N+1}\).

In particular, successive weighted Hessian cones cannot choose
independent null lines.  They are the successive coefficients of one
projective line.  On the other pivot chart the reciprocal slope gives
the same line, so the construction glues over the fraction field of
\(A\).

#### Proof

After the indicated normalization, \(a\) is a unit.  Put \(r=b/a\).
The Schur identity is

\[
 \det S=ac-b^2=a(c-ar^2).
\tag{5.2}
\]

Therefore \(\det S=0\) is equivalent to \(c=ar^2\), proving (5.1);
the displayed kernel follows immediately.

For the face-by-face assertion write

\[
 a=\sum_{i\ge0}a_i\epsilon^i,\quad
 b=\sum_{i\ge0}b_i\epsilon^i,\quad
 r=\sum_{i\ge0}r_i\epsilon^i.
\]

The equation \(b=ar\) determines the next slope coefficient uniquely:

\[
 a_0r_n=b_n-\sum_{i=1}^{n}a_i r_{n-i}.
\tag{5.3}
\]

Once \(r_0,\ldots,r_n\) are known, the next vanishing determinant face,
using (5.2), forces

\[
 c_n=\sum_{i+j+k=n}a_i r_jr_k.
\tag{5.4}
\]

Equations (5.3)--(5.4) prove synchronization successively and also prove
the truncated assertion. \(\square\)

### Corollary 5.2 (weighted Hessian form)

Suppose a weighted Hessian pencil has a nondegenerate quotient block.
Eliminate that block over its principal open and let \(S(\epsilon)\) be
the \(2\times2\) Schur complement on the residual kernel.  Every
successive vanishing determinant face synchronizes \(S\) by Lemma 5.1.
If the residual determinant vanishes identically, the resulting null
line is rational over the polynomial base.

Indeed, if the quotient block is \(Q\), block elimination gives

\[
 \det K_\rho=\det(Q)\det(S).
\]

Since \(\det(Q)\) is a unit on the chosen chart, the vanishing determinant
faces are exactly the corresponding Schur-complement faces.

This is an exact statement over the localized coefficient ring.  A
denominator-free polynomial direction requires a separate
algebraization argument.

## 6. Application to the rank-two sextic face

Write

\[
 \Psi=q_2+h_4+h_6
\]

as in the dense common-kernel audit.  Suppose
\(\operatorname{Hess}(h_6)\) has generic rank two and let \(W\) denote
its two-dimensional kernel over the relevant function field.  The
degree-twelve spatial determinant face is

\[
 \det\!\left(
 \overline{\operatorname{Hess}(h_6)}\bigm|V/W
 \right)
 \det\!\left(\operatorname{Hess}(h_4)\bigm|W\right).
\tag{6.1}
\]

The first factor is nonzero.  Hence

\[
 \det\!\left(\operatorname{Hess}(h_4)\bigm|W\right)=0.
\tag{6.2}
\]

On a chart which trivializes the quotient and kernel, choose any integral
weight and form the corrected Hessian pencil (4.2).  After the
nondegenerate sextic quotient is eliminated, (6.2) is a binary residual
symmetric-form identity.  Lemma 5.1 says that all its successive faces
define one rational projective null line on each pivot chart.  Thus the
rank-two residue is not a collection of unrelated facewise kernels: it
is a single synchronized rational cone

\[
 [a(X):b(X)]\in\mathbf P^1(\operatorname{Frac}K[X]).
\tag{6.3}
\]

There are two already-excluded subloci:

1. **constant cone:** \([a:b]\) is constant; this is excluded by the
   dense common-kernel theorem `HC4DCK`;
2. **dual-linear locus:** all nonlinear terms have dual degree at most
   one; this is the \(JC(2)\) locus of Section 3 and is excluded in the
   present degree range by Moh's theorem;

These subloci can meet.  Their complement, and the only part which
survives at the synchronization stage, is the **nonlinear moving-cone
locus**: the synchronized line (6.3) is nonconstant and some term has
dual degree at least two.  Section 8 excludes it using the next full
determinant face.

## 7. Why synchronization is the sharp conclusion

For

\[
 h_4=(xt+ym)^2
\]

the dual Hessian is

\[
 2\begin{pmatrix}x^2&xy\\xy&y^2\end{pmatrix}
\]

and its kernel is \((-y,x)\).  Give \(y\) weight one and the other
variables weight zero.  The weighted face series is

\[
 2x^2
 \begin{pmatrix}
 1&\epsilon y/x\\
 \epsilon y/x&\epsilon^2y^2/x^2
 \end{pmatrix}.
\tag{7.1}
\]

This is exactly (5.1) with \(r=\epsilon y/x\).  The middle face does not
choose an independent cone: its Schur square produces the final face.
But the synchronized direction is nonconstant and has a denominator on
the \(x\)-chart.

Therefore synchronization alone calls for an algebraization or exclusion
theorem for the nonlinear rational map (6.3), compatible with Hessian
integrability and the transported collision.  The degree-four
algebraization and exclusion are supplied next.

## 8. Cone-degree rigidity closes the moving locus

The required algebraization is forced in degree four.  This gives theorem
`HC4E46`.

### Lemma 8.1 (primitive cone-degree rigidity)

Let \(A=K[x,y]\), \(F=K(x,y)\), and let \(h_4\in A[t,m]\) be homogeneous
of total degree four.  Suppose

\[
 \det\operatorname{Hess}_{t,m}(h_4)=0.
\tag{8.1}
\]

Modulo a polynomial which is at most linear in \(t,m\), exactly one of the
following holds:

1. \(h_4\) has a nonzero constant second-derivative kernel direction; or
2. \[
    h_4=c\,(X^{\mathsf T}MU)^2
    +\text{terms of dual degree at most one},
   \tag{8.2}
   \]
   where \(c\in K^\times\), \(X=(x,y)\), \(U=(t,m)\), and
   \(M\in\operatorname{GL}_2(K)\).

#### Proof

Apply the two-variable singular-Hessian classification over \(F\).
Up to terms linear in \(U\),

\[
 h_4=g(a(X)t+b(X)m)
\tag{8.3}
\]

for \(a,b\in F\).  The synchronization lemma shows equivalently that the
dual-degree faces in (8.3) use one projective cone \([a:b]\).

Let \(d\ge2\) be the highest nonlinear dual degree.  Represent
\([a:b]=[p:q]\) by coprime homogeneous polynomials of the same degree
\(\delta\).  Such a representation exists because the highest face is
homogeneous in \(X\): scaling \(X\) preserves its projective binary
power, hence preserves \([a:b]\), so this rational map descends to
\(\mathbf P^1_X\).  The highest face is

\[
 \gamma(X)(p(X)t+q(X)m)^d.
\tag{8.4}
\]

Because \(p,q\) are coprime and both endpoint coefficients of (8.4) are
polynomial, every denominator of \(\gamma\) divides both \(p^d\) and
\(q^d\).  Hence \(\gamma\in A\).  Bihomogeneity gives

\[
 \deg\gamma=4-d-d\delta\ge0.
\tag{8.5}
\]

If \(\delta=0\), the cone is constant; a constant direction transverse
to \(pt+qm\) annihilates the second directional derivative of (8.3).
This gives the first case.  If the cone is moving, then \(\delta\ge1\).
Formula (8.5) has the unique solution

\[
 (d,\delta,\deg\gamma)=(2,1,0).
\]

Thus the moving term is the square of a bilinear form.  If \(p,q\) are
linearly dependent, its projective dual direction is constant.  Otherwise
they are the rows of an invertible matrix \(M\), giving (8.2).
\(\square\)

The classification input used here is only the characteristic-zero
two-variable singular-Hessian theorem in
[de Bondt's low-dimensional classification](https://arxiv.org/abs/1501.05168).
The denominator step is the displayed elementary UFD argument; no
formal-series algebraization is being assumed.

### Lemma 8.2 (uncancellable bilinear-square face)

Assume the moving alternative (8.2), and suppose the sextic is independent
of \(U\).  In

\[
 \det(H_0+zH_4+z^2H_6),
\tag{8.6}
\]

the dual-degree-four part of the coefficient of \(z^4\) is

\[
 48c^4\det(M)^2(X^{\mathsf T}MU)^4.
\tag{8.7}
\]

In particular it is nonzero.

#### Proof

For \(L=X^{\mathsf T}MU\), direct differentiation gives

\[
 \det\operatorname{Hess}(cL^2)
 =48c^4\det(M)^2L^4.
\tag{8.8}
\]

Every term of \(h_4-cL^2\) has dual degree at most one.  Its \(UU\),
\(UX\), and \(XX\) Hessian blocks therefore have dual degrees
\(-\infty,0,\le1\), respectively, whereas the corresponding blocks of
\operatorname{Hess}(cL^2)\) have degrees \(0,1,2\).
Replacing any entry of the determinant (8.8) by a lower quartic entry
strictly lowers dual degree.

The sextic Hessian occupies only the \(XX\) block.  A \(z^4\) determinant
term involving it contains two fewer quartic entries; the block row and
column constraints prevent those entries from attaining dual degree four.
The quadratic Hessian has dual degree zero.  Therefore nothing else in
(8.6) contributes to the indicated bidegree, proving (8.7).
\(\square\)

### Theorem 8.3 (even quartic--sextic closure)

Let

\[
 \Psi=q_2+h_4+h_6
\]

be an even four-variable potential over a characteristic-zero field, with
\(\operatorname{Hess}(q_2)\) nonsingular and \(h_4,h_6\) homogeneous of
degrees four and six.  If
\(\det\operatorname{Hess}\Psi\) is a nonzero constant, then
\(\nabla\Psi\) has no antipodal collision.

#### Proof

Extend scalars to an algebraic closure if necessary; the determinant
identity and a collision both persist.  Classify by the generic rank of
\(\operatorname{Hess}(h_6)\).

- **Rank three.**  The degree-fourteen spatial face supplies the common
  quartic kernel direction, so `HC4DCK` applies.
- **Rank two.**  Iterating the low-dimensional Gordan--Noether theorem
  makes the two-dimensional sextic kernel \(W\) a constant plane.  Choose
  it as the dual plane \(U\); then \(h_6=h_6(X)\).  The degree-twelve face
  is a nonzero quotient determinant times
  \(\det\operatorname{Hess}_U(h_4)\), so (8.1) holds.  Lemma 8.1 gives
  either a constant common kernel, excluded by `HC4DCK`, or (8.2).
  Lemma 8.2 excludes (8.2), because every nonconstant spatial determinant
  coefficient must vanish.
- **Rank one.**  The sextic depends on one linear form and has a constant
  three-plane kernel \(W\).  The degree-ten face makes the ternary
  Hessian of \(h_4\) on \(W\) singular.  The characteristic-zero
  three-variable singular-Hessian normal form, applied over the
  one-variable quotient field, supplies a nonzero rational direction
  \(v(s)\) for which \(D_{v(s)}^2h_4=0\) identically in the
  \(W\)-variables: after a linear change the Hessian is zero below the
  anti-diagonal, so its final diagonal entry is zero.  Let
  \(Z_s\subset\mathbf P(W)\) be the full scheme of such directions.
  Homogeneity gives
  \[
   \operatorname{Hess}_W(h_4)(\lambda s,\lambda U)
   =\lambda^2\operatorname{Hess}_W(h_4)(s,U).
  \]
  Since the identity defining \(Z_s\) holds for every \(U\), replacing
  \(U\) by \(\lambda U\) shows \(Z_{\lambda s}=Z_s\).  The quotient of
  the nonzero one-dimensional base by scaling is a point.  Thus a
  direction in the generic nonempty \(Z_s\) can be chosen constant; the
  polynomial identity extends across \(s=0\).  Again `HC4DCK` applies.
- **Rank zero.**  A homogeneous sextic with zero Hessian is zero.  The
  remaining quadratic-plus-quartic chart is `HC4HQ1`.

These cases exhaust the sextic Hessian ranks. \(\square\)

The theorem closes the even quartic--sextic chart with no cubic correction.
For simultaneous cubic--quartic--sextic corrections, `HC4T31` in
[`HC4_MENG_TRIPLE_RANK_THREE.md`](HC4_MENG_TRIPLE_RANK_THREE.md) closes
the rank-three sextic-Hessian stratum.  Sextic Hessian rank at most two
still has half-integral spatial layers, and the binary residual form is no
longer isolated by (6.1).

The rank-one normal-form input is the determinant-zero case used in
[de Bondt's three-variable Hessian theorem](https://arxiv.org/abs/1203.6605).
The homogeneous rank filtration is the low-dimensional
Gordan--Noether classification, also reviewed in the first source above.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_source_dual_bigrading.py
```

The checker verifies the \(35+84\) bidegree ledger, weighted Hessian
covariance for every quartic and sextic monomial under four weights, the
cotangent determinant (3.2), the synchronization recursion through order
eight, and the rotating-cone identity (7.1).

The all-order proof is the Schur identity (5.2); the finite symbolic
recursion is retained as an independent regression.

Replay the closure coefficient with:

```bash
.venv/bin/python scripts/verify_hc4_even_quartic_sextic_closure.py
```

The checker verifies the cone-degree ledger, the relative invariant
(8.8), and (8.7) with a generic quadratic Hessian, every compatible
dual-linear/source-only quartic Hessian block, and a generic source-only
sextic Hessian block present.  It also verifies the rank-one kernel
Hessian scaling identity on all \(35\) quartic monomials.
