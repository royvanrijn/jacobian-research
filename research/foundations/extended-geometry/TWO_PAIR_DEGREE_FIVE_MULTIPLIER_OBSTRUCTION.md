# The degree-five bilinear-multiplier obstruction

## 1. Statement

Use the two-pair quartic witness
\[
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right)
\]
from
[the two-pair counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md),
and let
\[
 L=aR+bZ+cW+eT
 \tag{1.1}
\]
be an arbitrary bilinear.  Thus \(LF\) has balanced degree five.

> **Theorem 1.1.** Over a characteristic-zero field, if
> \[
>  \mathcal E_2((LF)^m)=0\qquad(1\leq m\leq4),
>  \tag{1.2}
> \]
> then \(L=aR\).  Consequently every all-order pure-moment-zero point in
> the bilinear-multiplier family is the radial propagation \(aRF\).

In particular, multiplying the quartic seed by a noninvariant bilinear
cannot produce an \(R\)-primitive failure of \(\mathrm{MN}_5\).  This is a
theorem about the four-dimensional multiplier family (1.1), not a
classification of the full thirty-six-dimensional space \(V_5\).

## 2. The first three equations

Exact contraction gives
\[
\begin{aligned}
 \frac{\mu_1(LF)}{24}
 &=10b-3c,\\
 \frac{\mu_2(LF)}{190080}
 &=56ab+4ac+28b^2-4bc+3c^2+4e^2,\\
 \frac{\mu_3(LF)}{1045094400}
 &=3432a^2b+6864ab^2-156ac^2+1144b^3\\
 &\quad+1716b^2c+130bc^2-49c^3-260ce^2.
\end{aligned}
\tag{2.1}
\]
The first equation gives \(c=10b/3\).  After this substitution, the next
two normalized polynomials become
\[
 \frac43(52ab+36b^2+3e^2)
 \tag{2.2}
\]
and
\[
 \frac{8b}{27}
 (11583a^2+17316ab+21916b^2-2925e^2).
 \tag{2.3}
\]

If \(b=0\), then \(c=0\), and (2.2) gives \(e=0\).  This is exactly the
radial line \(L=aR\).

## 3. Excluding the nonradial branch at moment four

Suppose \(b\ne0\).  By homogeneity set \(b=1\), write \(u=e/b\), and solve
(2.2):
\[
 a=-\frac{36+3u^2}{52}.
 \tag{3.1}
\]
Then (2.3) is a nonzero scalar multiple of
\[
 q(u)=8019u^4-623736u^2+3219760.
 \tag{3.2}
\]
On the same substitution, the fourth moment is a nonzero scalar multiple
of
\[
 p(u)=136323u^6-5359284u^4-174020976u^2-802761152.
 \tag{3.3}
\]
Exact Euclidean reduction over \(\mathbb Q\) gives
\[
 \gcd_{\mathbb Q[u]}(p,q)=1.
 \tag{3.4}
\]
Equivalently,
\[
 \operatorname{Res}_u(p,q)
 =
 97842802725670657299880334741299150646484103418606547061702656
 \ne0.
 \tag{3.5}
\]
Thus the nonradial branch has no point over an algebraic closure, proving
Theorem 1.1.

## 4. What this blocks and what remains

The theorem closes the most direct degree-five operations on the quartic
seed:

1. multiplication by any balanced bilinear;
2. every linear combination of the invariant, phase, and height
   bilinears \(R,Z,W,T\);
3. any proposed degree-one angular factor before homogenization.

It does not cover a Cartan projection followed by correction in lower
summands, a sum of several independently raised profiles, a nonlinear
differential intertwiner, or a general point of \(V_5\).  Those operations
do not have the form \(LF\).

The conclusion is equivariant.  The contraction-preserving
\(\mathrm{PGL}_2\)-action fixes \(R\) and acts invertibly on the four
bilinears.  Hence for every nonzero scalar \(s\) and every
\(g\in\mathrm{PGL}_2\), the first four moments of
\[
 L\,(s\,gF)
\]
vanish only when \(L\) is proportional to \(R\).  In particular, this
holds uniformly on the nondegenerate two-parameter family \(F_{a,b}\)
from the
[local-moduli calculation](TWO_PAIR_COUNTEREXAMPLE_LOCAL_MODULI.md),
whose members are scaled diagonal translates of \(F\).

The next primitive search should therefore begin with either a genuine
two-profile sum in \(V_5\) or an intrinsic trace-free ansatz, rather than
another multiplicative lift of \(F\).

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_two_pair_degree_five_multiplier_obstruction.py
```

The checker constructs \(F\) and a symbolic \(L\), derives the first four
contractions exactly, verifies (2.1)--(3.3), and checks both the polynomial
gcd and the displayed resultant.  The calculation is an exact
characteristic-zero certificate, not a bounded numerical search.
