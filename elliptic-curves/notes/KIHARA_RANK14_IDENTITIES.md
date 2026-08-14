# Algebra behind Kihara's rank-at-least-14 family

## Scope

This note separates the exact polynomial identities in Kihara's construction
from the later independence argument.  It proves why the twelve root sections
and the three additional printed abscissas lie on the quartic.  It does **not**
give a new rank lower bound: the rank-at-least-14 status still comes from the
independent specialization certificate described in
[`families/kihara_rank14.json`](../families/kihara_rank14.json).

The formulas are checked by
[`cas/derive_kihara_rank14_identities.py`](../cas/derive_kihara_rank14_identities.py)
over exact rational-function fields.

## The intrinsic six-center condition

Let

\[
A(X)=\prod_{i=1}^6(X-a_i),\qquad
F_u(X)=A(X-u)A(X+u)
       =\prod_{i=1}^6\bigl(X-(a_i+u)\bigr)\bigl(X-(a_i-u)\bigr).
\]

There is a unique monic sextic (G_u(X)) for which (G_u^2-F_u) has
degree at most five: recursively match the coefficients in degrees 12 down
through 6.  Put (r_u=G_u^2-F_u).  If

\[
A(X)=X^6+c_5X^5+c_4X^4+c_3X^3+c_2X^2+c_1X+c_0,
\]

then exact coefficient comparison gives

\[
[X^5]r_u=-u^2\left(
24c_1-8c_2c_5-12c_3c_4+7c_3c_5^2+8c_4^2c_5
-6c_4c_5^3+c_5^5\right).
\]

Translate the centers by their mean, so (z_i=a_i-\bar a) and

\[
B(Z)=\prod_i(Z-z_i)
=Z^6+c_4Z^4+c_3Z^3+c_2Z^2+c_1Z+c_0.
\]

The obstruction becomes

\[
[X^5]r_u=-12u^2(2c_1-c_3c_4).
\]

Thus, on the nondegenerate locus (u\ne0), the paired degree-twelve
construction has quartic remainder if and only if

\[
\boxed{2c_1=c_3c_4}.
\]

Writing (e_j=e_j(z_1,\ldots,z_6)), this is the translation-invariant
relation

\[
\boxed{2e_5=e_2e_3},
\]

because (c_4=e_2, c_3=-e_3, c_1=-e_5).  This is the essential algebraic
relation among the six centers.  Neither (c_2) nor (c_0) enters the
degree-five obstruction.

For Kihara's displayed (a_i(p,q)), define

\[
\begin{aligned}
\alpha={}&16p^8+56p^7q+116p^6q^2+110p^5q^3+86p^4q^4\\
&+41p^3q^5+82p^2q^6+24pq^7+12q^8,\\
\beta={}&p(p-q)(p+2q)(2p-q)(2p+5q)(4p+q)\\
&\quad\cdot(p^2+2q^2)(2p^2+2pq+3q^2)(4p^2+4pq+3q^2).
\end{aligned}
\]

Direct expansion of the mean-centered polynomial gives

\[
c_4=-\frac{2\alpha}{3},\qquad
c_3=\frac{2\beta}{27},\qquad
c_1=-\frac{2\alpha\beta}{81}
     =\frac{c_3c_4}{2}.
\]

This short factorization is the reason the degree-five remainder disappears.
At each of the twelve roots (b_i=a_i\mathbin\pm u), (F_u(b_i)=0), hence

\[
r_u(b_i)=G_u(b_i)^2.
\]

That accounts for the first twelve sections without any further square test.

## The thirteenth section is already generic in (p,q,u)

Set

\[
\begin{aligned}
D&=2p^2+2pq+3q^2,\\
L&=2p^2+4pq+5q^2,\\
M&=8p^6+28p^5q+58p^4q^2+69p^3q^3
   +76p^2q^4+40pq^5+22q^6.
\end{aligned}
\]

Kihara's thirteenth abscissa is

\[
x_{13}=\frac{Lu+M}{D}.
\]

Write (H=H_0+H_1u+H_2u^2+H_3u^3), where

\[
\begin{aligned}
H_0={}&q(p+q)(4p^2+5q^2)(2p^2+pq+8q^2)D^2\\
&\quad\cdot(2p^2+3pq+2q^2)L(3p^2+2pq+2q^2),\\
H_1={}&2(64p^{12}+576p^{11}q+2656p^{10}q^2+8192p^9q^3
 +19120p^8q^4\\
&\quad+35264p^7q^5+52339p^6q^6+62738p^5q^7
 +60769p^4q^8\\
&\quad+45996p^3q^9+26455p^2q^{10}+10260pq^{11}+2352q^{12}),\\
H_2={}&4L(8p^6+40p^5q+84p^4q^2+125p^3q^3
 +157p^2q^4+105pq^5+48q^6),\\
H_3={}&24q(p+q)(2p^2+3pq+4q^2).
\end{aligned}
\]

The exact identity in (mathbf Q(p,q,u)) is

\[
\boxed{
r_u(x_{13})=
\left(
\frac{2puq^2(2p-q)(p+q)^2(2p^2+pq+2q^2)H}{D^2}
\right)^2.}
\]

Consequently (P_{13}) does **not** require Kihara's later one-parameter
base change.  It is extra structure of the particular two-parameter center
formulas, not a formal consequence of (2e_5=e_2e_3) alone.

## The last two squares use the one-parameter base change

Now impose

\[
\begin{aligned}
p&=t^2(8+3t^2),\\
q&=-6(t^2+2)(t^2+4),\\
u&=\frac{4(t^2+2)U(t)V(t)}{t},
\end{aligned}
\]

with

\[
\begin{aligned}
U(t)&=9t^8+150t^6+928t^4+2400t^2+2304,\\
V(t)&=18t^8+201t^6+860t^4+1632t^2+1152.
\end{aligned}
\]

Substitution of the printed (x_{14}(t)) and (x_{15}(t)) into (r_t)
factors as exact squares.  To state the result compactly, put

\[
\begin{gathered}
A=t^2+3,\quad B=3t^2+4,\quad C=3t^2+8,\quad D=t^2+4,\\
E=t^2+2,\quad J=3t^4+28t^2+48,\\
W=27t^8+339t^6+1580t^4+3168t^2+2304,
\qquad \kappa=2359296.
\end{gathered}
\]

Let (K_{14}) and (K_{15}) be the primitive degree-58 and degree-44
integer factors produced by the exact factorization.  Their ascending-
coefficient SHA-256 values (newline separated, with a final newline) are

```text
K14  b6b5bb73a584cacddcbe1e45d0fb7e839f1489ffaffa3f2db14e606d0784326a
K15  bedb18dc27516b36336fc7690319a05e2b32b1052ca4068c5cd0412860fcd585
```

The verifier's `--full-polynomials` option prints every coefficient.  With
those canonical primitive factors, it derives

\[
\boxed{
r_t(x_{14})=
\left(
\frac{\kappa ABCD^2E^4J^2UV^3K_{14}}{tW^2}
\right)^2,}
\]

and

\[
\boxed{
r_t(x_{15})=
\left(
\frac{\kappa ABCD^2E^3J^2UV^3K_{15}}{t}
\right)^2.}
\]

These two identities hold in (mathbf Q(t)).  The displayed base change is
a sufficient rational parametrization making them squares; the computation
does not claim it is the unique or necessary parametrization.

## Essential relations versus coordinate artifacts

- **Intrinsic and necessary for the quartic degree drop:** after centering,
  (2e_5=e_2e_3) (equivalently (2c_1=c_3c_4)) for (u\ne0).
- **Construction-specific sufficient structure:** Kihara's full formulas
  (a_i(p,q)) imply the generic (P_{13}) identity.  They describe a much
  smaller locus than the single intrinsic hypersurface.
- **Base-change-specific sufficient structure:** the displayed (p(t),q(t),u(t))
  forces the two additional (P_{14},P_{15}) squares.  This is not a
  classification of all ways to obtain them.
- **Artifacts of coordinates:** setting (a_1=0) is a translation choice;
  scaling ((p,q,u,x,y)\mapsto(\lambda p,\lambda q,\lambda^4u,
  \lambda^4x,\lambda^{24}y)) changes normalization only; permuting the six
  centers or replacing (u) by (-u) preserves the root set.
- **Artifact of the rank proof:** choosing (P_{15}) as the group-law origin
  is useful for the fourteen-point independence certificate but plays no role
  in any square identity above.

## Exact replay

With SymPy 1.14 installed:

```bash
python3 elliptic-curves/cas/derive_kihara_rank14_identities.py
python3 elliptic-curves/cas/derive_kihara_rank14_identities.py --full-polynomials
```

The first command is compact; the second also emits the full degree-58 and
degree-44 factors and exact rational-function ordinates.  Both reconstruct
the degree-twelve product and its square approximant from the stored primary
formulas rather than trusting pre-expanded quartic coefficients.
