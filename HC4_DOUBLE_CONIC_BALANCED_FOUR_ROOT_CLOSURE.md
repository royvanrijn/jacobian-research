# Double-conic balanced four-root closure

## Status

This note proves `HC4NHM18`.  It removes the finite exceptional-cross-ratio
locus left by `HC4NHM15` in the binary-decic partition ((3,3,2,2)).  Thus
every four-support decic is excluded from the clean double-conic target.
Forms with at least five support points remain open, and no Schur solution or
HC4 counterexample is constructed here.

The exact replay is

```bash
.venv/bin/python scripts/verify_hc4_double_conic_balanced_four_root_closure.py
```

## 1. Normal form and residual line

Put

\[
 q=xz-y^2,
 \qquad
 f=s^3t^3(s-t)^2(s-\lambda t)^2,
 \qquad \lambda\ne0,1,
\tag{1.1}
\]

and write

\[
 h_5=\operatorname{lift}(f)+qG_3
\tag{1.2}
\]

with (G_3) a general ternary cubic.  Use the monomial ordering

\[
G_3=g_9x^3+g_8x^2y+g_7x^2z+\cdots
       +g_4xz^2+g_1yz^2+g_0z^3,
\tag{1.3}
\]

as fixed by the checker.  Let (D=\det\operatorname{Hess}(h_5)), and define
the prospective residual-line coefficients by

\[
 A=[x^5z^4]D,
 \qquad B=[x^4yz^4]D,
 \qquad C=[x^4z^5]D.
\tag{1.4}
\]

Set

\[
 R=D-q^4(Ax+By+Cz).
\tag{1.5}
\]

The clean target (D=q^4\ell), (ell\ne0), is equivalent to (R=0)
and ((A,B,C)\ne(0,0,0)).

## 2. The two endpoint chains

Five coefficients at the (x)-endpoint give

\[
\begin{aligned}
 [x^9]R&=32g_9^3,\\
 [x^7yz]R\big|_{g_9=0}&=12(g_8+1)^3,\\
 [x^6yz^2]R\big|_{g_9=0,\,g_8=-1}
   &=-144(g_7-2\lambda-2)^2,\\
 A\big|_{g_9=0,\,g_8=-1,\,g_7=2\lambda+2}&=0.
\end{aligned}
\tag{2.1}
\]

Over characteristic zero, (R=0) therefore forces

\[
 g_9=0,
 \qquad g_8=-1,
 \qquad g_7=2\lambda+2,
 \qquad A=0.
\tag{2.2}
\]

At the opposite endpoint,

\[
\begin{aligned}
 [z^9]R&=32g_0^3,\\
 [xyz^7]R\big|_{g_0=0}&=12(g_1+\lambda^2)^3,\\
 [x^2yz^6]R\big|_{g_0=0,\,g_1=-\lambda^2}
   &=-144\lambda^2(g_4-2\lambda^2-2\lambda)^2,\\
 C\big|_{g_0=0,\,g_1=-\lambda^2,
                 \,g_4=2\lambda^2+2\lambda}&=0.
\end{aligned}
\tag{2.3}
\]

Since (lambda\ne0), these force

\[
 g_0=0,
 \qquad g_1=-\lambda^2,
 \qquad g_4=2\lambda^2+2\lambda,
 \qquad C=0.
\tag{2.4}
\]

No division by a possible exceptional polynomial occurs in either chain.

## 3. The middle contradiction

After (2.2) and (2.4), introduce

\[
 u=\lambda^2+4\lambda+g_5+1,
 \qquad
 v=g_6-2\lambda-2.
\tag{3.1}
\]

The residual line and two further normal-layer coefficients reduce exactly
to

\[
 \boxed{
 B=16u^3,
 \qquad
 [x^4y^5]R=-16(3u-v^2),
 \qquad
 [x^4y^4z]R=-16v(2u-v^2).
 }
\tag{3.2}
\]

Because (A=C=0), a nonzero residual line would require (B\ne0), hence
(u\ne0).  But (R=0) first gives (v^2=3u).  The last equation then
becomes (uv=0), so (v=0); substituting back gives (u=0), a
contradiction.

For regression, the adjacent coefficients are

\[
 [x^4y^3z^2]R=-32u(u-v^2),
 \qquad
 [x^4y^2z^3]R=32u^2v.
\tag{3.3}
\]

## 4. Conclusion

> **Theorem `HC4NHM18` -- Balanced four-root closure.**  Let
> (q=xz-y^2), and suppose that the restriction of a ternary quintic to
> (q=0) has binary root partition ((3,3,2,2)).  Then no nonzero linear
> form (ell) satisfies
> \[
> \det\operatorname{Hess}(h_5)=q^4\ell.
> \]
> Hence the complete ((3,3,2,2)) row is empty, including every special
> cross-ratio fiber.

Together with `HC4NHM15`, this excludes every binary decic supported on at
most four points from the clean double-conic packet.  The next double-conic
target is the family of decics supported on at least five points.  Only a
survivor of those normal-layer equations should proceed to the Schur system.
The invariant continuation and the necessary residual-line saturation are
formulated in
[`HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md`](HC4_DOUBLE_CONIC_INVARIANT_SATURATION_GATE.md).
