# Independent-parity quartic obstruction below the GVC(3) witness

## 1. The complete independent-linear repair

Work over a characteristic-zero field. On the affine sphere quadric

\[
xy+t^2=1,
\]

write Long's three-variable sphere polynomial as

\[
L=E+O,
\]

with antipodally even and odd parts

\[
E=xy-2t^2-x^2t^2,
\qquad
O=y-3xt^2.
\tag{1.1}
\]

The earlier common-linear parity repair multiplied both odd terms by one
linear form. The complete independent-linear version instead takes

\[
H_1=a_1x+b_1y+c_1t,
\qquad
H_3=a_3x+b_3y+c_3t,
\]

and

\[
\boxed{
F=\alpha E+yH_1-3xt^2H_3.
}
\tag{1.2}
\]

It is the sphere restriction of the homogeneous quartic

\[
\widetilde F=
\alpha\bigl(\rho(xy-2t^2)-x^2t^2\bigr)
+\rho yH_1-3xt^2H_3,
\qquad \rho=xy+t^2.
\tag{1.3}
\]

Thus (1.2) allows independent linear profiles on the degree-one and
degree-three odd pieces of Long's polynomial. It is larger than the previous
common-factor repair, but it is not the full fifteen-dimensional space of
ternary quartics.

Let

\[
\mu_m(F)=\int_{S^2}F^m\,d\sigma.
\tag{1.4}
\]

For homogeneous quartics these are, up to nonzero radial constants, exactly
\(\Delta^{2m}(\widetilde F^m)\), where
\(\Delta=4\partial_x\partial_y+\partial_t^2\).

## 2. The first moment and projective chart cover

Direct phase extraction and the height integral give

\[
\boxed{
\mu_1=\frac{2}{15}(5a_1-3b_3).
}
\tag{2.1}
\]

Hence every pure-zero point satisfies

\[
b_3=\frac53a_1.
\tag{2.2}
\]

The remaining projective parameter space is covered by three cases:

1. \(\alpha\ne0\), normalized by \(\alpha=1\);
2. \(\alpha=0,a_1\ne0\), normalized by \(a_1=1\);
3. \(\alpha=a_1=0\), hence also \(b_3=0\).

This cover is exhaustive, including projective infinity.

## 3. Exact characteristic-zero elimination on the two open charts

On \(\alpha=1\), substitute (2.2) and use variables

\[
(a_1,b_1,c_1,a_3,c_3).
\]

The primitive numerators of \(\mu_2,\ldots,\mu_6\) generate the unit ideal
over \(\mathbb Q\). An exact monic Buchberger replay reaches the constant
one after 119 critical pairs; the basis has 73 elements immediately before
the unit appears.

On \(\alpha=0,a_1=1\), the primitive numerators of
\(\mu_2,\ldots,\mu_5\), in variables

\[
(b_1,c_1,a_3,c_3),
\]

also generate the unit ideal over \(\mathbb Q\). The exact replay reaches
one after 26 critical pairs, with 22 basis elements immediately before the
unit.

No modular inference is used in either statement. The checker implements
ordinary Buchberger reduction over \(\mathbb Q\), uses only the safe product
criterion, and processes every remaining critical pair.

## 4. The deep boundary

It remains to put

\[
\alpha=a_1=b_3=0.
\]

Then

\[
F=b_1y^2+c_1yt-3a_3x^2t^2-3c_3xt^3.
\tag{4.1}
\]

The next three moments are

\[
\mu_2=-\frac4{35}(4a_3b_1+3c_1c_3),
\tag{4.2}
\]

\[
\mu_3=-\frac8{385}(11a_3c_1^2-15b_1c_3^2),
\tag{4.3}
\]

and

\[
\mu_4=\frac{48}{5005}
\left(
48a_3^2b_1^2+120a_3b_1c_1c_3+35c_1^2c_3^2
\right).
\tag{4.4}
\]

Equation (4.2) gives

\[
a_3b_1=-\frac34c_1c_3.
\tag{4.5}
\]

Substitution into the parenthesis in (4.4) gives exactly

\[
-28c_1^2c_3^2.
\tag{4.6}
\]

Thus

\[
c_1c_3=0,
\qquad
a_3b_1=0.
\tag{4.7}
\]

Equation (4.3) now rules out every crossed choice. Consequently either

\[
b_1=c_1=0
\tag{4.8}
\]

or

\[
a_3=c_3=0.
\tag{4.9}
\]

Both alternatives plainly satisfy every pure moment: in (4.8) all phase
weights are positive, while in (4.9) they are all negative.

## 5. Complete radical and all-order termination

Let \(I_6=(\mu_1,\ldots,\mu_6)\) in the seven parameter coordinates. The
three-chart argument and the two explicit one-sided components give

\[
\boxed{
\sqrt{I_6}
=(\alpha,a_1,b_3,b_1a_3,b_1c_3,c_1a_3,c_1c_3).
}
\tag{5.1}
\]

Equivalently,

\[
\sqrt{I_6}
=(\alpha,a_1,b_3,b_1,c_1)
\cap
(\alpha,a_1,b_3,a_3,c_3).
\tag{5.2}
\]

The components are terminal for every fixed multiplier. On the first,

\[
F=-3xt^2(a_3x+c_3t),
\]

so every monomial has positive phase weight. On the second,

\[
F=y(b_1y+c_1t),
\]

so every monomial has negative phase weight. The phase weight is preserved
by \(\Delta\). If a fixed multiplier has ordinary degree \(q\), then the
output of \(\Delta^{2m}(Q\widetilde F^m)\) has degree at most \(q\), while
its phase weight has absolute value at least \(m-q\). It is therefore zero
once

\[
m>2q.
\tag{5.3}
\]

Hence no component in (5.2) is a GVC counterexample.

> **Theorem 5.1 — independent-linear parity obstruction.** The complete
> family (1.2) contains no homogeneous quartic GVC(3) counterexample. Its
> first six pure moments cut out exactly the two one-sided terminal
> components (5.2).

This strengthens the common-linear obstruction. It does not classify
repairs using a genuine harmonic cubic multiplier on the degree-one piece,
or arbitrary homogeneous ternary quartics.

## 6. Next frontier

The original endpoint-contact family cannot begin below \(\Delta^6\), and
now neither the common nor independent linear homogenizations produce a
quartic witness for \(\Delta^2\). The next natural parity repair is degree
six:

\[
F=\alpha E+H_3O,
\tag{6.1}
\]

with a complete homogeneous cubic profile \(H_3\). The radial-linear part
\(H_3=\rho H_1\) is already closed by Theorem 5.1; only the harmonic-cubic
profile is genuinely new.

## 7. Reproduction

Run

```bash
python3 scripts/verify_gvc3_independent_parity_quartic.py
```

The checker derives the moments from (1.2), performs both exact
characteristic-zero Buchberger calculations, verifies the boundary
factorization (4.5)--(4.6), and writes the radical and component data to the
generated artifact.
