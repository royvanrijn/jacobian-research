# The balanced cubic two-variable GVC theorem

## 1. Statement

Let \(k\) be a characteristic-zero field.  Let \(\Lambda\) be a homogeneous
constant-coefficient differential operator of order three in \(x,y\), and
let \(P\in k[x,y]\) be homogeneous of degree three.

> **Theorem 1.1.** If
> \[
>  \Lambda^m(P^m)=0,\qquad 1\leq m\leq4,
>  \tag{1.1}
> \]
> then all pure contractions vanish and, for every \(Q\in k[x,y]\),
> \[
>  \Lambda^m(QP^m)=0
>  \qquad\text{whenever }m>\deg Q.
>  \tag{1.2}
> \]

Thus the balanced homogeneous cubic stratum satisfies GVC in two
variables.  Under the Segre embedding

\[
 f(\zeta,z)=\sigma_\Lambda(\zeta)P(z),
 \tag{1.3}
\]

the first four SIC moments force \(f\) into the pair-linear one-sided
nullcone.  Consequently every fixed SIC multiplier, including multipliers
depending on both dual and coordinate variables, also contracts to zero
for all sufficiently large powers.

This theorem treats homogeneous order and polynomial degree exactly three;
it did not by itself prove unrestricted \(\operatorname{GVC}(2)\), which is
now supplied by the independent
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md).  It does not prove
\(\operatorname{SIC}(2)\).

## 2. Covariance and the three symbol orbits

Write the binary cubic symbol of the operator as

\[
 A(u,v)=\sigma_\Lambda(u,v).
\]

The scalar \(\Lambda^m(P^m)\) is the natural apolar pairing of \(A^m\) and
\(P^m\).  It is unchanged by the simultaneous contragredient
\(\operatorname{GL}_2\)-action, up to a nonzero scalar that does not affect
vanishing.  We may extend scalars to an algebraic closure: a polynomial
identity proved there descends to \(k\).

The zero symbol is trivial.  Every nonzero binary cubic has one of three
root-multiplicity types.  After a linear change and rescaling, take

\[
 A=u^3,\qquad A=u^2v,\qquad\text{or}\qquad A=uv(u+v).
 \tag{2.1}
\]

Write

\[
 P=ax^3+bx^2y+cxy^2+dy^3.
 \tag{2.2}
\]

The proof now consists of three elementary calculations.

## 3. Triple-root symbol

Let \(A=u^3\), so \(\Lambda=\partial_x^3\).  The first moment is

\[
 \Lambda(P)=6a.
 \tag{3.1}
\]

Hence (1.1) gives \(a=0\), and \(P\) has \(x\)-degree at most two.  For a
fixed \(Q\),

\[
 \deg_x(QP^m)\leq \deg Q+2m<3m
\]

when \(m>\deg Q\).  Therefore
\(\partial_x^{3m}(QP^m)=0\), proving (1.2) in this orbit.  The same degree
inequality with \(Q=1\) gives every remaining pure moment.

## 4. Double-root symbol

Let \(A=u^2v\), so \(\Lambda=\partial_x^2\partial_y\).  Direct contraction
gives

\[
\begin{aligned}
 \Lambda(P)&=2b,\\
 \Lambda^2(P^2)&=48(2ac+b^2),\\
 \Lambda^3(P^3)&=4320(3a^2d+6abc+b^3).
\end{aligned}
\tag{4.1}
\]

The first three vanishings imply

\[
 b=0,\qquad ac=0,\qquad a^2d=0.
\tag{4.2}
\]

If \(a=0\), then

\[
 P=y^2(cx+dy)
\tag{4.3}
\]

has \(x\)-degree at most one.  Thus

\[
 \deg_x(QP^m)\leq\deg Q+m<2m
\]

for \(m>\deg Q\), and the factor \(\partial_x^{2m}\) kills the product.

If \(a\ne0\), equations (4.2) give \(c=d=0\), so \(P=ax^3\).  Then
\(\deg_y(QP^m)=\deg_y Q<m\), and the factor \(\partial_y^m\) kills the
product.  This proves (1.2) in the double-root orbit.

## 5. Squarefree symbol

Let

\[
 A=uv(u+v),\qquad
 \Lambda=\partial_x\partial_y(\partial_x+\partial_y).
\tag{5.1}
\]

The first moment is

\[
 \Lambda(P)=2(b+c),
\]

so put \(c=-b\).  After removing nonzero integer factors, the next three
moments are

\[
\begin{aligned}
 R_2={}&-2ab+3ad-b^2+2bd,\\
 R_3={}&ad(a+d),\\
 R_4={}&14a^2b^2-42a^2bd+45a^2d^2
          +12ab^3-40ab^2d+42abd^2\\
      &\quad+3b^4-12b^3d+14b^2d^2.
\end{aligned}
\tag{5.2}
\]

More precisely, the scalar moments are \(48R_2\), \(12960R_3\), and
\(414720R_4\).  The equation \(R_3=0\) gives three branches.

### Branch \(a=0\)

Here

\[
 R_2=b(2d-b),\qquad
 R_4=b^2(3b^2-12bd+14d^2).
\]

The choice \(b=2d\) gives \(R_4=8d^4\), so the common zero set has \(b=0\).
It consists of \(P=dy^3\).

### Branch \(d=0\)

Here

\[
 R_2=-b(2a+b),\qquad
 R_4=b^2(14a^2+12ab+3b^2).
\]

The choice \(b=-2a\) gives \(R_4=8a^4\), so the common zero set has \(b=0\).
It consists of \(P=ax^3\).

### Branch \(d=-a\)

Here

\[
\begin{aligned}
 R_2&=-(a+b)(3a+b),\\
 R_4&=(3a+b)^2(5a^2+6ab+3b^2).
\end{aligned}
\]

The choice \(b=-a\) gives \(R_4=8a^4\), while \(b=-3a\) gives

\[
 P=a(x-y)^3.
\]

Consequently the full common zero set of the first four moments in the
squarefree orbit is exactly

\[
 P=ax^3,\qquad P=dy^3,\qquad P=a(x-y)^3.
\tag{5.3}
\]

Each form in (5.3) is annihilated by one of the three commuting linear
factors of \(\Lambda\):

\[
 \partial_y(x^3)=0,\qquad
 \partial_x(y^3)=0,\qquad
 (\partial_x+\partial_y)(x-y)^3=0.
\tag{5.4}
\]

Let \(D\) be the corresponding factor.  Since \(D(P)=0\),

\[
 D^m(QP^m)=P^mD^mQ=0
\]

for \(m>\deg Q\).  The remaining commuting factors of \(\Lambda^m\) do not
change zero.  This proves (1.2) in the squarefree orbit and completes the
proof of Theorem 1.1.

## 6. Nullcone and frontier consequences

The normal forms above are visibly one-sided:

- for \(u^3\), the polynomial side omits the opposed extreme weight;
- for \(u^2v\), equations (4.3) and \(P=ax^3\) have a strict weight
  separator; and
- for the squarefree symbol, a linear factor of the operator annihilates
  the pure cube on the polynomial side.

Thus

\[
 \Sigma_3\cap V(\mu_1,\mu_2,\mu_3,\mu_4)
 =\Sigma_3\cap N_3,
\tag{6.1}
\]

where \(\Sigma_3\) is the rank-one Segre cone and \(N_3\) is the two-pair
one-sided nullcone.  In particular, the extra semistable component of the
first thirteen full bidegree-\((3,3)\) moment equations does not meet the
Segre cone.

Combined with the complete bidegree-\((2,2)\) SIC theorem, this initially
moved the first open balanced homogeneous separable GVC stratum to degree
four.  The later
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) now proves the GVC
conclusion in every balanced homogeneous degree.  The present theorem
remains stronger in giving a four-moment premise, exact Segre-nullcone
equality, and the cutoff \(m>\deg Q\).  Nonhomogeneous pairs and genuinely
nonseparable SIC forms remain open.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_two_variable_cubic_gvc.py
```

The checker constructs the apolar moments independently from coefficient
expansion, verifies (4.1) and (5.2), checks every squarefree branch
factorization, and verifies the three annihilating directions in (5.4).
The degree arguments proving the all-order conclusion are the written
proof above.
