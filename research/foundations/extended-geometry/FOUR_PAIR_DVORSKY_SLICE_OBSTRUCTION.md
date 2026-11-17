# A full-symbol four-pair obstruction in the Dvorsky slice

## 1. Statement and scope

Put

\[
 P=(t+a+b+d)(ad+bt)
\]

and let

\[
\begin{aligned}
 R={}&A\partial_a^2+B\partial_a\partial_b
      +C\partial_a\partial_d+D\partial_b^2\\
    &+E\partial_b\partial_d+F\partial_d^2,\qquad
 \Lambda=\partial_tR.                                      \tag{1.1}
\end{aligned}
\]

Thus \(R\) is the general ternary quadratic constant-coefficient operator;
there is no coefficient box.

> **Theorem 1.1.**  Over a characteristic-zero field, if \(R\ne0\) and
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),                           \tag{1.2}
> \]
> then, up to scalar, \(R\) is one of
> \[
>  \boxed{
>  \partial_a^2,\quad \partial_d^2,\quad
>  (\partial_a-\partial_b)^2,\quad
>  (\partial_b-\partial_d)^2.}                              \tag{1.3}
> \]
> The first eight identities in (1.2) already force this conclusion
> set-theoretically.

For \(R=0\), one has \(f=0\), so the Mathieu conclusion is immediate.

Every operator in (1.3) satisfies a strict one-sided weight condition.
Consequently, for

\[
 \lambda(\zeta)
 =\zeta_tR(\zeta_a,\zeta_b,\zeta_d),\qquad
 f=\lambda(\zeta)P(z),                                      \tag{1.4}
\]

one has

\[
 \mathcal E(gf^m)=0
\]

for every fixed \(g\in\mathbb C[\zeta,z]\) and all sufficiently large \(m\).
Hence no counterexample to the Mathieu property of

\[
 \mathcal M_4
 =\sum_{x\in\{t,a,b,d\}}(\partial_x-\zeta_x)
   \mathbb C[\zeta,z]                                      \tag{1.5}
\]

occurs anywhere in the full six-parameter symbol slice (1.1) with this
fixed \(P\).

This is not a proof of \(\operatorname{SIC}(4)\).  It does not treat other
cubic \(P\), quadratic symbols involving \(\partial_t\), nonfactorized
symbols, or nonseparable \(f\).  It upgrades one natural slice of the
bounded search in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](DVORSKY_GVC5_COUNTEREXAMPLE.md)
from an integral-box experiment to an exact characteristic-zero theorem.

## 2. The eight moment equations

For \(m\geq1\), define

\[
 c_m(A,B,C,D,E,F)=\Lambda^m(P^m).                           \tag{2.1}
\]

Both \(\lambda\) and \(P\) are homogeneous of degree three, so \(c_m\) is a
scalar homogeneous polynomial of degree \(m\) in the six coefficients.
The exact checker constructs it without symbolic differentiation.  If

\[
 (k_A,k_B,k_C,k_D,k_E,k_F)\in\mathbb N^6,\qquad
 \sum k_\bullet=m,
\]

then the corresponding term of \(\lambda^m\) has derivative exponent

\[
 k_A(1,2,0,0)+k_B(1,1,1,0)+k_C(1,1,0,1)
 +k_D(1,0,2,0)+k_E(1,0,1,1)+k_F(1,0,0,2),                 \tag{2.2}
\]

in the variable order \((t,a,b,d)\).  Reading the matching coefficient of
\(P^m\) and multiplying by the derivative factorial gives \(c_m\) exactly.

Let

\[
 H=(c_1,\ldots,c_8)\subset\mathbb Q[A,B,C,D,E,F].            \tag{2.3}
\]

The four claimed projective points all meet the chart
\(A+D+F\ne0\).  Normalize by

\[
 A+D+F=1
\]

and put

\[
 I=H+(A+D+F-1).                                             \tag{2.4}
\]

The normalized representatives of (1.3) are

\[
\begin{array}{c|rrrrrr}
&A&B&C&D&E&F\\ \hline
\partial_a^2&1&0&0&0&0&0\\
\partial_d^2&0&0&0&0&0&1\\
(\partial_a-\partial_b)^2&1/2&-1&0&1/2&0&0\\
(\partial_b-\partial_d)^2&0&0&0&1/2&-1&1/2.
\end{array}                                                 \tag{2.5}
\]

Their reduced ideal \(J\) has Gröbner basis

\[
\begin{aligned}
J=\langle&
 C,\ B+2D+E,\ A+D+F-1,\ 4F^2-E-4F,\\
&2EF-E,\ 4DF+E,\ E^2+E,\ 2DE-E,\ 2D^2-D
\rangle.                                                    \tag{2.6}
\end{aligned}
\]

Exact ordered reduction gives

\[
 I\subseteq J.                                              \tag{2.7}
\]

Conversely, the eighth power of each generator in (2.6) lies in \(I\);
for \(A+D+F-1\) the first power suffices, and for \(2D^2-D\) the fourth
power suffices.  Therefore

\[
 J\subseteq\sqrt I.
\]

Since \(J\) is the intersection of the four distinct rational maximal
ideals in (2.5), it is radical.  Combining this with (2.7) gives

\[
 \boxed{\sqrt I=J.}                                        \tag{2.8}
\]

As auxiliary scheme data,

\[
 \dim_{\mathbb Q}\mathbb Q[A,\ldots,F]/I=64,\qquad
 \dim_{\mathbb Q}\mathbb Q[A,\ldots,F]/J=4.                 \tag{2.9}
\]

The nilpotent thickness in \(I\) is why a literal comparison of ideals
would be wrong.

It remains to justify that the normalization misses no projective
component.  The homogeneous ideal

\[
 H+(A+D+F)
\]

is zero-dimensional.  A homogeneous zero-dimensional affine cone is
supported only at the origin.  Thus every nonzero solution of \(H\) has
\(A+D+F\ne0\), and (2.8) proves Theorem 1.1.

## 3. All-order Mathieu vanishing

We use one elementary weight lemma.  Suppose that, after an invertible
pair-linear change of variables, every monomial of \(\lambda(\zeta)\) has
weight \(u\), every monomial of \(P(z)\) has weight at most \(v<u\), and
the coordinate weights are nonnegative.  For a fixed polynomial \(g\),
the \(\zeta\)-weight of every monomial of \(g\lambda^mP^m\) eventually
exceeds its \(z\)-weight.  But

\[
 \mathcal E(\zeta^\alpha z^\beta)\ne0
 \quad\Longrightarrow\quad \alpha\leq\beta
\]

coordinatewise, which is impossible under nonnegative weights.  Hence
\(\mathcal E(gf^m)=0\) for all sufficiently large \(m\).

For \(R=\partial_a^2\), assign weight one to \(t,a\) and zero to \(b,d\).
The symbol \(\zeta_t\zeta_a^2\) has weight three, while every monomial of
\(P\) has weight at most two.  The case \(R=\partial_d^2\) is identical
with weight one on \(t,d\).

For \(R=(\partial_a-\partial_b)^2\), put

\[
 a=x,\qquad b=y-x.
\]

Then \(\partial_x=\partial_a-\partial_b\) and

\[
 P=(t+y+d)\bigl(x(d-t)+yt\bigr),                            \tag{3.1}
\]

which has \(x\)-degree at most one, whereas the transformed symbol has
\(\zeta_x\)-degree two.  Weight only \(x\).  Similarly, for
\(R=(\partial_b-\partial_d)^2\), put

\[
 b=x,\qquad d=y-x.
\]

Now

\[
 P=(t+a+y)\bigl(x(t-a)+ay\bigr),                            \tag{3.2}
\]

and the same degree gap applies.  Pair-linear changes preserve the
contraction identity and \(\mathcal M_4\), so all four cases satisfy the
Mathieu conclusion for arbitrary fixed \(g\).

This proves the theorem.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_four_pair_dvorsky_slice_obstruction.py
```

The checker constructs \(c_1,\ldots,c_8\) by exact sparse arithmetic,
uses Singular over \(\mathbb Q\) for the two ideal containments and chart
complement, and verifies the four strict weight gaps.
