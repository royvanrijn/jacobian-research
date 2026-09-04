# Finite certificates and an effective cutoff for binary GVC

## 1. Statement and scope

Let $k$ be a field of characteristic zero, let

\[
 \Lambda=\lambda(\partial_x,\partial_y),\qquad
 0\ne\lambda\in k[X,Y],\qquad 0\ne P\in k[x,y],
\]

and assume \(\lambda(0,0)=0\). Write
\(r=\min\{|\alpha|:\alpha\in\operatorname{supp}\lambda\}\),
\(R=\deg\lambda\), and \(d=\deg P\).

> **Theorem 1.1 (binary support certificate).** The identities
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1)
>  \tag{1.1}
> \]
> hold if and only if there are a linear coordinate change over **$k$**
> and positive integers $u,v\leq d+R$ such that, in those coordinates,
> \[
>  g:=\min_{\alpha\in\operatorname{supp}\lambda}(u\alpha_x+v\alpha_y)
>       -\max_{\beta\in\operatorname{supp}P}(u\beta_x+v\beta_y)>0.
>  \tag{1.2}
> \]
> If $d\geq r$, the Hall direction is unique and defined over $k$.
> Taking it as the first coordinate direction, with any complement over
> $k$, permits $u>v>0$.

For $d<r$, no coordinate change is needed: use $u=v=1$.
The theorem strengthens the conclusion extracted from the existing
[Hall-envelope proof](BINARY_GVC_ENVELOPE_CLOSURE.md); its necessity uses
that proof's Hall localization and shifted-ray endpoint theorem. This is
an internal written theorem, without external review or a complete Lean
proof. The exact computations below verify certificates and regressions;
they do not prove the universal necessity direction.

> **Corollary 1.2 (degree-only mixed cutoff).** Under (1.1), every nonzero
> $Q\in k[x,y]$ satisfies
> \[
>  \Lambda^m(QP^m)=0\qquad\bigl(m>(d+R)\deg Q\bigr).
>  \tag{1.3}
> \]
> More precisely, a certificate (1.2) gives vanishing for
> \[
>  m\geq\left\lfloor\frac{K_Q}{g}\right\rfloor+1,
>  \qquad K_Q=\max_{\gamma\in\operatorname{supp}Q}(u\gamma_x+v\gamma_y),
>  \tag{1.4}
> \]
> with $Q$ expressed in the same coordinates.

If \(\lambda=0\), $P=0$, or $Q=0$, the corresponding contractions
vanish directly. A nonzero constant term of \(\lambda\), with $P\ne0$,
already makes \(\Lambda(P)\ne0\), by its highest ordinary-degree part.

## 2. Hall localization, uniqueness, and descent

Suppose $d\geq r$ and (1.1). The homogeneous part of degree $m(d-r)$
gives

\[
 \lambda_r(\partial)^m(P_d^m)=0\qquad(m\geq1).
 \tag{2.1}
\]

Over an algebraic closure factor
\(\lambda_r=\prod_{i=1}^rD_{v_i}\) and
\(P_d=\prod_{j=1}^dL_j\), absorbing nonzero scalar factors. For a generic
translation $z$, coefficient extraction identifies the constant term of

\[
 \left(\frac{P_d(z+\sum_i t_i v_i)}{t_1\cdots t_r}\right)^m
\]

with \(\lambda_r(\partial)^m(P_d^m)(z)/(m!)^r\).
The Duistermaat--van der Kallen theorem therefore excludes
\((1,\ldots,1)\) from the Newton polytope of the numerator. A matching
from the derivative copies $i$ to distinct factors $j$ with
\(L_j(v_i)\ne0\) would put that vector in the Minkowski sum of the factor
polytopes, which equals the product polytope. Hence Hall's theorem gives
a deficient subset of derivative copies.

A subset containing two nonparallel directions sees every one of the
$d$ factors. Since \(d\geq r\geq |I|\), such a subset cannot be
deficient. Thus one direction of full multiplicity $e$ annihilates at
least \(d-e+1\) factors of $P_d$. This proves the Hall lemma also when
$d=r$; the earlier lemma stated only the range $d>r$.

Call a direction with these multiplicities a **Hall candidate**. There
cannot be two distinct candidates. If their symbol multiplicities are
$e_1,e_2$ and their annihilating linear factors in $P_d$ have
multiplicities $h_1,h_2$, then

\[
 e_1+e_2\leq r,\qquad h_1+h_2\leq d,
 \qquad
 h_1+h_2\geq 2d-e_1-e_2+2\geq2d-r+2\geq d+2,
\]

a contradiction. This uniqueness holds whenever $d\geq r$, independently
of whether the full pair satisfies (1.1).

The set of Hall candidates is invariant under every automorphism of
\(\overline k/k\). A unique candidate is therefore fixed. Its projective
coordinate is algebraic over $k$, and characteristic zero makes that
extension separable, so the coordinate belongs to $k$. The direction
thus descends to $k$.

Choose it as the first coordinate vector and any independent vector over
$k$ as the second. In these coordinates

\[
 \lambda_r=X^eC_{r-e}(X,Y),\quad C_{r-e}(0,1)\ne0,
 \qquad P_d=y^{d-e+1}D_{e-1}(x,y).
 \tag{2.2}
\]

In particular the minimum $X$-exponent of \(\lambda_r\) is $e$,
and the maximum $x$-exponent of $P_d$ is some $t<e$.

## 3. A strict separator in that same frame

Use the envelopes

\[
 L(s)=\min_{\alpha\in\operatorname{supp}\lambda}(s\alpha_x+\alpha_y),
 \qquad
 U(s)=\max_{\beta\in\operatorname{supp}P}(s\beta_x+\beta_y),
 \qquad \Delta(s)=U(s)-L(s).
\]

If $d=r$, then \(\Delta(1)=0\), and immediately to the right of $1$
its slope is $t-e<0$. Thus \(\Delta(s)<0\) for some rational $s>1$.

If $d>r$, [the envelope proof, §§4--5](BINARY_GVC_ENVELOPE_CLOSURE.md)
gives a first rational $s_*>1$ with \(\Delta(s_*)=0\). Throughout the
preceding positive component, the complete operator face lies strictly to
the right of the complete polynomial face in the $x$-coordinate: the
shifted-ray theorem forbids overlap or reversal before the first zero.
At $s_*$, the zero-output case of the same theorem makes the equality
faces disjoint. Both faces contain their left-hand limiting active
exponents. Those exponents are ordered with the operator to the right,
so the two complete equality faces must have that same order. In particular

\[
 \min_{\alpha\in\operatorname{supp}A_{s_*}}\alpha_x
 >\max_{\beta\in\operatorname{supp}B_{s_*}}\beta_x.
\]

By finiteness of the supports, immediately to the right of $s_*$ the
right derivative of \(\Delta\) is precisely the right side's maximum
minus the left side's minimum, and is strictly negative. Hence again
\(\Delta(s)<0\) for some rational $s>1$. All coordinates used here
remain in $k$; only the necessity proof passes through scalar extension.

Conversely, (1.2) kills each monomial selection in
\(\Lambda^m(P^m)\): its total differential weight exceeds the input
weight by at least $mg>0$. A derivative that acts nontrivially cannot
increase this weight. No cancellation argument is needed for this direction.

## 4. A small integral weight

It remains to bound the weights, rather than merely clear an unspecified
denominator. Each support pair \(\alpha,\beta\) imposes the strict
inequality

\[
 as+b>0,\qquad a=\alpha_x-\beta_x,\quad b=\alpha_y-\beta_y.
\]

Together with $s>1$, the feasible slopes form a nonempty open interval
\((\ell,h)\), possibly with $h=+\infty$. An inequality with $a=0$
must already have $b>0$. Write \(\ell=p/q\) in lowest terms. Either
\(\ell=1\), or a lower-bound inequality with $a>0$ gives

\[
 \ell=\frac{\beta_y-\alpha_y}{\alpha_x-\beta_x}.
\]

Thus \(1\leq p\leq d\) and \(1\leq q\leq R\), including
\(\ell=1\) since $d\geq r\geq1$. If the upper endpoint is finite,
write \(h=p'/q'\) in lowest terms. Since $h>\ell\geq1$, it arises
from $a<0$ and satisfies

\[
 h=\frac{\alpha_y-\beta_y}{\beta_x-\alpha_x},
 \qquad 1\leq p'\leq R,\quad1\leq q'\leq d.
\]

The mediant
\(s=(p+p')/(q+q')\) lies strictly between distinct endpoints, so
\(w=(u,v)=(p+p',q+q')\) has both entries at most $d+R$.
If $h=+\infty$, take \(w=(p+q,q)\), corresponding to
$s=\ell+1$, with the same bound. Dividing out a common factor is optional.
The gap $g$ is a positive integer, hence at least $1$.

For $d<r$, the weight $(1,1)$ has gap $r-d>0$ and satisfies the
same bound. This finishes Theorem 1.1. A linear coordinate change preserves
total degree, so
\(K_Q\leq\max(u,v)\deg Q\leq(d+R)\deg Q\).
Every selected mixed term has differential weight greater than its input
weight once $mg>K_Q$. This proves both cutoff formulas.

## 5. A finite rational-arithmetic decision procedure

For rational inputs the proof gives this terminating procedure:

1. Handle zero inputs, a nonzero constant symbol, and $d<r$ directly.
2. Factor \(\lambda_r(T,1)\) over \(\mathbb Q\), including the possible
   direction at infinity. Test its rational linear factors for the Hall
   multiplicity condition (2.2). There is at most one candidate.
3. If there is none, reject the pure premise. A nonrational candidate
   would have at least two distinct conjugates, contradicting uniqueness.
4. In the unique candidate's rational frame, solve the finitely many
   strict linear inequalities above. An empty interval rejects the pure
   premise. Otherwise the mediant gives a positive certificate and the
   multiplier cutoff.

For a finite root $a$, the implementation uses

\[
 M=\begin{pmatrix}1&0\\-a&1\end{pmatrix},\qquad
 P'=P(Mz)=P(x,y-ax),\qquad
 \lambda'=\lambda(M^{-T}(X,Y))=\lambda(X+aY,Y).
\]

The direction at infinity uses the coordinate swap. This inverse-transpose
transport is essential: the symbol and polynomial cannot be substituted by
the same matrix. All coefficients and interval endpoints use exact rational
arithmetic. Nonlinear irreducible factors are retained in the output's
factorization record.

The implementation is
[`jcsearch/binary_gvc_certificate.py`](../jcsearch/binary_gvc_certificate.py).
Its negative answers mean **the all-powers premise fails**, not that GVC
fails. Negative completeness uses the written theorem; the routine does
not supply a uniform bound on the first nonzero pure power and does not
decide by checking a prefix of powers. Positive certificates are checked
independently by symbolic coordinate substitution and support comparison.

Replay with:

```bash
.venv/bin/python scripts/verify_binary_gvc_finite_certificate.py
```

The pinned [artifact](../artifacts/generated-results/binary-gvc-finite-certificate-v1.json)
contains 12 rational fixtures, all 2025 ordered pairs of one- or two-element
supports in positive degrees through three, six instances of the family
below, and rejected weight/coordinate-transport mutations. These finite
regressions check the implementation, not the universal theorem. The default
replay is read-only and requires byte-identical output; `--write` deliberately
regenerates it. The verifier also accepts `--input path.json`, where `symbol`,
`polynomial`, and optional `multiplier` are lists of `[x-exponent, y-exponent,
"rational coefficient"]` rows. Floating-point coefficients are rejected.

## 6. Degree dependence and a sharp family cutoff

Let $d,q\geq1$, and set

\[
 \Lambda=\partial_x+\partial_y^{d+1},\qquad P=y^d,\qquad Q=x^q.
\]

Every pure contraction vanishes: any \(\partial_x\) kills $P^m$,
while the remaining derivative has $y$-order $(d+1)m>dm$.
In a mixed contraction a term with $j$ copies of \(\partial_x\)
can survive only if

\[
 j\leq q,\qquad(d+1)(m-j)\leq dm,
 \quad\text{equivalently }m\leq(d+1)j.
\]

Consequently every $m>(d+1)q$ vanishes. At $m_0=(d+1)q$, only
$j=q$ survives and gives the nonzero constant

\[
 \Lambda^{m_0}(QP^{m_0})
 =\binom{m_0}{q}\,q!\,(dm_0)!.
\]

The separator $(u,v)=(d+1,1)$ has gap $1$, and (1.4) gives exactly
$m_0+1$. Thus the certificate cutoff can be sharp, and no bound
$m>C\deg Q$ with an absolute degree-independent $C$ works for all
binary pairs. This does **not** claim optimality of the general coefficient
$d+R$ in (1.3).

## 7. Proof dependencies and formal boundary

The necessity direction depends on the translated
Duistermaat--van der Kallen/Hall argument and on the shifted-ray and
no-reversal arguments of `GVC2UFT` and `GVC2ENV`. Uniqueness, ground-field
descent, the equal-degree extension, the mediant bound, and the exact
algorithm are the additional deductions recorded here as `GVC2SC`.

[`formal/gvc/GVC/BinaryReduction.lean`](../formal/gvc/GVC/BinaryReduction.lean)
already proves the finite-support implication from an ordered common face
to a strict positive weight and eventual mixed vanishing. It does not
currently prove Hall-direction descent, the $d+R$ bound, or completeness
of this rational decision procedure. No new result here changes the formal
status of the unrestricted binary theorem. No literature-priority claim is
made for the extracted classification or bound.
