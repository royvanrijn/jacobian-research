# The minimal Hasse-failing geometric degree is five

This note sharpens the degree-eight construction in
[HASSE_PRINCIPLE_KELLER_FIBER.md](HASSE_PRINCIPLE_KELLER_FIBER.md).
It combines the weighted-suspension theorem with the classical minimal
intersective polynomial of Berend and Bilu.

## The theorem

There is a polynomial Keller map

\[
 F:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q}
\]

of geometric degree five and a rational target \(q\) such that

\[
 F^{-1}(q)(\mathbb Q)=\varnothing,
\qquad
 F^{-1}(q)(\mathbb R)\ne\varnothing,
\qquad
 F^{-1}(q)(\mathbb Q_p)\ne\varnothing
\]

for every prime \(p\). No finite etale scheme of degree at most four can
have this local-global behavior. Consequently

\[
 \boxed{d_{\mathrm{HP}}=5.}
\]

Here \(d_{\mathrm{HP}}\) is the least geometric degree of a Keller map having
a complete regular fiber which is everywhere locally soluble but has no
rational point.

## The minimal intersective polynomial

Start with

\[
 f(X)=(X^3-19)(X^2+X+1).
\]

It has no rational root. For a prime \(p\ne3\):

- if \(p\equiv1\pmod3\), then \(X^2+X+1\) has a simple root modulo \(p\);
- if \(p\equiv2\pmod3\), then cubing is a bijection of
  \(\mathbb F_p^\times\), so \(X^3-19\) has a simple root;
- \(p=2\) is included in the second case.

At \(p=3\), the unit \(19\) lies in \(1+9\mathbb Z_3\), which is the image
of the cube map on \(1+3\mathbb Z_3\). Thus \(f\) has a root in every
\(\mathbb Q_p\). Its cubic factor has one real root.

Berend and Bilu prove both this example and the impossibility of degree
less than five in
[Polynomials with roots modulo every integer](https://doi.org/10.1090/S0002-9939-96-03210-8).
The lower bound is purely arithmetic, so it applies to every complete finite
etale Keller fiber.

## Tangent normalization by an affine change

Put \(X=10-27W\) and divide by \(27\). The resulting primitive polynomial is

\[
\begin{aligned}
P(W)={}&-531441W^5+1003833W^4-758889W^3\\
      &\quad+286497W^2-53901W+4033.
\end{aligned}
\]

The affine change preserves the finite etale scheme over every
characteristic-zero field. Direct calculation gives

\[
 P(1)-P(0)=P'(0)=-53901.
\]

Remove the affine part:

\[
 H(W)=P(W)-P(0)-P'(0)W
 =-531441W^5+1003833W^4-758889W^3+286497W^2.
\]

Then

\[
 H(0)=H'(0)=H(1)=0,\qquad
 c=-H'(1)=345546,
\]

and

\[
 \frac{H''(1)}c=-\frac{586}{79}\ne-2.
\]

Thus \(H\) is an admissible weighted seed. Its source coefficient is

\[
 a_0=-\frac{1+H''(1)/c}{2+H''(1)/c}=-\frac{507}{428}.
\]

## The Keller map and complete fiber

Write \(h=H'\) and set

\[
 v=xy,\qquad S=x^2z,\qquad
 \gamma=1-\frac{507}{428}v+S,\qquad
 W=(1+v)\gamma.
\]

Define

\[
\begin{aligned}
F_1&=\frac1{x^2}
 \left(1+v+\frac{Wh(W)-H(W)}{c\gamma^2}\right),\\
F_2&=\frac1x\left(c+\frac{h(W)}\gamma\right),\\
F_3&=x\gamma.
\end{aligned}
\]

The weighted polynomiality theorem makes these polynomial coordinates.
Their total degrees are \(17,16,4\), and

\[
 \det DF=345546.
\]

Take

\[
 q=\left(\frac{4033}{345546},\,53901,\,1\right).
\]

The inverse pencil is

\[
 H(W)-BCW+cAC^2,
\]

which equals \(P(W)\) at \(q\). Since \(P\) is squarefree and \(C=1\), every
root is in the regular reconstruction chart and no boundary point can occur.
Therefore

\[
 \boxed{F^{-1}(q)\simeq\operatorname{Spec}\mathbb Q[W]/(P).}
\]

The local and global point assertions follow from the corresponding
assertions for \(f\).

## Verification

Run

```bash
.venv/bin/python scripts/verify_minimal_hasse_keller_fiber.py
```

The checker verifies the local residue covering, the \(3\)-adic lifts,
absence of rational roots, squarefreeness, tangent normalization,
weighted admissibility, polynomial expansion, constant Jacobian,
suspension identities, target conversion, and quotient-ring reconstruction.
