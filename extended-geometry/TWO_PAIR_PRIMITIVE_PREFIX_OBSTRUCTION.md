# Primitive sharp finite prefixes in every balanced degree

## 1. The family

Let \(F\) be the quartic two-pair counterexample, let
\[
 A_d=R^{d-4}F\in V_d,
 \qquad d\geq4,
 \tag{1.1}
\]
and put
\[
 G_{d,\lambda}=A_d+\lambda Z^d.
 \tag{1.2}
\]

> **Theorem 1.1.** For every \(d\geq4\),
> \[
>  \mathcal E_2(G_{d,\lambda}^m)=0
>  \qquad(1\leq m\leq d),
>  \tag{1.3}
> \]
> while
> \[
>  \boxed{
>  \mathcal E_2(G_{d,\lambda}^{d+1})
>  =(d+1)\lambda\,
>  (d(d+1)+1)!\,
>  \frac{d!}{(2d+1)!!}.}
>  \tag{1.4}
> \]
> Hence for \(\lambda\ne0\), \(G_{d,\lambda}\) is an
> \(R\)-primitive point whose first \(d\) pure moments vanish and whose
> next pure moment is nonzero.

This is a sharp finite-prefix obstruction, not an Image-Mathieu
counterexample: the all-order premise fails exactly at the displayed next
moment.

The same argument closes the whole positive-phase triangular sector.
Let
\[
 H=\sum_{j=1}^d c_jR^{d-j}Z^j
 \tag{1.5}
\]
be nonzero, and let \(s\) be the least index with \(c_s\ne0\).

> **Theorem 1.2.** For \(B_d=R^{d-4}F+H\),
> \[
>  \mathcal E_2(B_d^m)=0\qquad(1\leq m\leq s),
>  \tag{1.6}
> \]
> and
> \[
>  \boxed{
>  \mathcal E_2(B_d^{s+1})
>  =(s+1)c_s\,
>  (d(s+1)+1)!\,
>  \frac{s!}{(2s+1)!!}\ne0.}
>  \tag{1.7}
> \]

Consequently no nonzero correction in (1.5) preserves all pure moments.
Among this entire \(d\)-dimensional Borel sector, \(H=\lambda Z^d\)
uniquely maximizes the zero prefix length.

The first genuinely two-sided degree-five correction is also rigid.  Put
\[
 C_{a,b}=RF+aZ^5+bW^5.
 \tag{1.8}
\]

> **Theorem 1.3.** If
> \[
>  \mathcal E_2(C_{a,b}^m)=0\qquad(m=2,4,6),
>  \tag{1.9}
> \]
> then \(a=b=0\).

The three exact branch formulas are
\[
\begin{aligned}
 \mathcal E_2(C_{a,b}^2)
 &=921600\,ab,\\
 \mathcal E_2(C_{0,b}^4)
 &=-12009388769280000\,b,\\
 \mathcal E_2(C_{a,0}^6)
 &=569547266090245736292679680000000\,a.
\end{aligned}
\tag{1.10}
\]
The second moment first gives \(ab=0\).  The fourth and sixth formulas
then close the two branches.  Thus adding the two extreme weights of the
primitive decic quotient does not produce a degree-five continuation.

Odd height produces a longer false prefix.  Define
\[
 J_{d,\lambda}=R^{d-4}F+\lambda Z^{d-1}T.
 \tag{1.11}
\]

> **Theorem 1.4.** For every \(d\geq4\),
> \[
>  \mathcal E_2(J_{d,\lambda}^m)=0
>  \qquad(1\leq m<2d),
>  \tag{1.12}
> \]
> whereas
> \[
>  \boxed{
>  \mathcal E_2(J_{d,\lambda}^{2d})
>  =\binom{2d}{2}\lambda^2(2d^2+1)!
>  \frac{(2d-2)!}{(4d-1)!!}.}
>  \tag{1.13}
> \]

For \(\lambda\ne0\), this is again \(R\)-primitive.  Thus any uniform
consecutive-moment exclusion of primitive points in \(V_d\) must reach
order at least \(2d\), not merely \(d+1\).

The same parity argument closes every opposite odd-height monomial pair.
Let
\[
 r=d-s\geq1\quad\text{be odd},\qquad
 U=Z^sT^r,\quad V=W^sT^r,
 \tag{1.14}
\]
and set
\[
 K_{a,b}=R^{d-4}F+aU+bV.
 \tag{1.15}
\]

> **Theorem 1.5.** If all pure moments of \(K_{a,b}\) vanish, then
> \(a=b=0\).

The exact branch certificates are
\[
\begin{aligned}
 \mathcal E_2(K_{a,b}^2)
 &=2ab(2d+1)!
 \frac{s!(2r-1)!!}{(2d+1)!!},\\
 \mathcal E_2(K_{a,0}^{2s+2})
 &=\binom{2s+2}{2}a^2(d(2s+2)+1)!
 \frac{(2s)!(2r-1)!!}{(4s+2r+1)!!},\\
 \mathcal E_2(K_{0,b}^{s+2})
 &=\binom{s+2}{2}b^2(d(s+2)+1)!
 \frac{(-1)^s(2s)!(2d-1)!!}
 {2^s(4s+2d+1)!!}.
\end{aligned}
\tag{1.16}
\]
All displayed scalar factors are nonzero in characteristic zero.  The
first equation gives \(ab=0\), and the other two close the branches.
For \(d=5,s=4,r=1\), this is the interior pair
\(Z^4T,W^4T\), detected at moments \(2,10,6\).

Even height is detected earlier.  Retain (1.14)--(1.15), now with
\(r=d-s\geq0\) even, and assume \(s\geq3\).  Put
\[
 n=\lceil s/2\rceil,\qquad
 \gamma_s=
 \begin{cases}
  1,&s\ \text{even},\\
  3n,&s\ \text{odd}.
 \end{cases}
 \tag{1.17}
\]

> **Theorem 1.6.** If all pure moments of \(K_{a,b}\) vanish, then
> \(a=b=0\).

Here the three exact certificates are
\[
\begin{aligned}
 \mathcal E_2(K_{a,b}^2)
 &=2ab(2d+1)!
 \frac{s!(2r-1)!!}{(2d+1)!!},\\
 \mathcal E_2(K_{a,0}^{s+1})
 &=(s+1)a(d(s+1)+1)!
 \frac{s!(r-1)!!}{(2s+r+1)!!},\\
 \mathcal E_2(K_{0,b}^{n+1})
 &=(n+1)b(d(n+1)+1)!
 \frac{\gamma_s(-1)^ns!(r+2n-1)!!}
 {2^n(r+2n+2s+1)!!}.
\end{aligned}
\tag{1.18}
\]
We use the convention \((-1)!!=1\).  The assumption \(s\geq3\) makes
the first line free of terms linear in \(a,b\).  It gives \(ab=0\), and
the other two nonzero formulas close the branches.  In particular this
subsumes the extreme pair \(Z^d,W^d\).

There is one even-height interior pair left in degree five after Theorems
1.5--1.6:
\[
 L_{a,b}=RF+aZT^4+bWT^4.
 \tag{1.19}
\]

> **Theorem 1.7.** If
> \[
>  \mathcal E_2(L_{a,b}^m)=0\qquad(m=2,3,4),
>  \tag{1.20}
> \]
> then \(a=b=0\).

After removing nonzero scalar factors, the three moments are
\[
\begin{aligned}
 q_2={}&70ab+198a-165b,\\
 q_3={}&490ab+858a+21b^2-845b,\\
 q_4={}&8580a^2b^2+35112a^2b+45220a^2\\
 &-36036ab^2+45220ab+100776a\\
 &+52269b^2-142120b.
\end{aligned}
\tag{1.21}
\]
An exact lexicographic basis for \((q_2,q_3)\) is
\[
 528a-21b^2-310b,\qquad
 b(735b^2+12929b-12870).
\tag{1.22}
\]
The normal remainder of \(q_4\) is
\[
 \frac{4b(155631189b+54155455)}{94325}.
\tag{1.23}
\]
The resultant of the quadratic and linear factors in (1.22)--(1.23) is
\[
 -418538718730248905250\ne0.
\tag{1.24}
\]
Thus \(b=0\), and then (1.22) gives \(a=0\).  Theorems 1.5--1.7 together
exclude every opposite monomial pair
\[
 Z^sT^{5-s},\quad W^sT^{5-s}\qquad(1\leq s\leq5)
\tag{1.25}
\]
as an all-moment-preserving correction of \(RF\).

## 2. Phase proof

On the Hopf sphere write
\[
 A_d=S^dp,\qquad Z^d=S^dx^d.
 \tag{2.1}
\]
The phase support of \(p\) has minimum weight \(-1\).  In a term containing
\(k\geq1\) copies of \(x^d\) and \(m-k\) copies of \(p\), the minimum
phase weight is
\[
 dk-(m-k)=(d+1)k-m.
 \tag{2.2}
\]
For \(m\leq d\), this is positive.  Every term involving \(x^d\) therefore
has zero phase constant term, while the term \(p^m\) has zero expectation
by the quartic all-order identity.  This proves (1.3).

At \(m=d+1\), only the term with \(k=1\) can have phase weight zero.
The full lower-jet identity for the quartic seed gives
\[
 \mathbb E_U(x^dp^d)=\frac{d!}{(2d+1)!!}.
 \tag{2.3}
\]
The radial degree is \(d(d+1)\), so
\[
 \mathbb E(S^{d(d+1)})=(d(d+1)+1)!.
 \tag{2.4}
\]
Multiplying (2.3)--(2.4) by the binomial coefficient \(d+1\) and by
\(\lambda\) proves (1.4).

For Theorem 1.2, every angular monomial of \(H\) has phase weight at least
\(s\).  Replacing \(d\) by \(s\) in the phase estimate shows that all mixed
terms vanish through order \(s\).  At order \(s+1\), only one copy of the
lowest term \(c_sx^s\) can contribute.  The lower-jet identity
\[
 \mathbb E_U(x^sp^s)=\frac{s!}{(2s+1)!!}
 \tag{2.5}
\]
and the radial moment
\[
 \mathbb E(S^{d(s+1)})=(d(s+1)+1)!
 \tag{2.6}
\]
give (1.7).  Higher phase terms and products of two or more terms from
\(H\) have strictly positive phase.

For Theorem 1.4, the angular correction is \(x^{d-1}t\).  Terms containing
an odd number of corrections integrate to zero by \(t\)-parity.  A term
with an even number \(k\geq2\) of corrections has minimum phase
\[
 (d-1)k-(m-k)=dk-m.
 \tag{2.7}
\]
This is positive for \(m<2d\).  At \(m=2d\), only \(k=2\) contributes, and
the phase boundary forces every one of the remaining \(2d-2\) copies of
\(p\) to use
\[
 [x^{-1}]p=\frac{1-t^2}{2}.
 \tag{2.8}
\]
Consequently
\[
\mathbb E_U\left(x^{2d-2}t^2p^{2d-2}\right)
=2^{-(2d-2)}\int_0^1t^2(1-t^2)^{2d-2}\,dt
=\frac{(2d-2)!}{(4d-1)!!}.
\tag{2.9}
\]
The radial moment is \((2d^2+1)!\), and choosing the two corrections
contributes \(\binom{2d}{2}\lambda^2\).  This proves (1.13).

For Theorem 1.5, the mixed product satisfies
\[
 UV=(xy)^st^{2r}
 =2^{-s}(1-t^2)^st^{2r},
\]
which gives the first line of (1.16).  On either single branch, odd powers
of the correction vanish by height parity.  On the positive branch, two
copies first reach phase zero at order \(2s+2\); the minimum coefficient
of \(p\) gives
\[
 2^{-2s}\int_0^1t^{2r}(1-t^2)^{2s}\,dt
 =\frac{(2s)!(2r-1)!!}{(4s+2r+1)!!}.
\tag{2.10}
\]
On the negative branch, two copies first meet the maximum phase of \(p\)
at order \(s+2\).  Since
\[
 y=\frac{1-t^2}{2x},
 \qquad [x^2]p=-\frac{t^2}{2},
\]
the boundary integral is
\[
 \frac{(-1)^s(2s)!(2d-1)!!}
 {2^s(4s+2d+1)!!}.
\tag{2.11}
\]
This proves the remaining two lines of (1.16).

For Theorem 1.6, the first line of (1.18) is the same mixed-product
calculation.  Since \(s\geq3\), neither \(A_dU\) nor \(A_dV\) has phase
zero.  On the positive branch, one correction first meets the minimum
phase of \(p\) at order \(s+1\), and
\[
 2^{-s}\int_0^1t^r(1-t^2)^s\,dt
 =\frac{s!(r-1)!!}{(2s+r+1)!!}.
\tag{2.12}
\]
On the negative branch, one correction first meets the maximum phase at
order \(n+1\).  If \(s=2n\), the boundary term is
\((-t^2/2)^n\).  If \(s=2n-1\), it is
\[
 n\left(-\frac{3t^2}{2}\right)
 \left(-\frac{t^2}{2}\right)^{n-1}.
\]
Multiplication by \(y^st^r=2^{-s}x^{-s}(1-t^2)^st^r\) and beta
integration give the last line of (1.18).  Terms with two or more
same-sign corrections cannot reach phase zero at either first-detection
order.

For \(d>4\), setting \(R=0\) kills \(A_d\), whereas \(Z^d\) survives at
\[
 (\xi _1,\xi _2,z_1,z_2)=(1,1,1,-1).
\]
Thus \(G_{d,\lambda}\) is not divisible by \(R\) when \(\lambda\ne0\).
For \(d=4\), its top phase term likewise shows primitivity.

## 3. Consequences

1. No moment cutoff at or below \(d\) can distinguish the all-moment-zero
   locus from primitive points in \(V_d\).
2. In degree five, the first five moments do not exclude primitive points,
   even though the bilinear-multiplier family \(LF\) is already forced
   radial by its first four moments.
3. Any proposed consecutive cutoff for the primitive quotient must include
   a moment of order at least \(2d\).
4. Positive-phase triangular corrections cannot turn radial propagation
   into an all-order primitive witness; a successful correction must mix
   both phase signs or introduce a different height profile.
5. No opposite monomial pair in degree five preserves all pure moments.
   Any successful two-sided correction must mix several weights or use
   independent nonmonomial height profiles.
6. Odd height can conceal failure twice as long as a pure highest-weight
   perturbation, so parity must be incorporated into search cutoffs.
7. In every degree, every opposite odd-height pair and every opposite
   even-height pair of phase at least three are excluded.  Only the
   low phases \(s=1,2\) remain open in the uniform even-height argument;
   Theorem 1.7 separately closes the phase-one case in degree five.

The family also explains why adding a highest-weight correction to radial
propagation looks successful in short exact searches: the failure is
delayed linearly with the balanced degree.

## Reproduction

Run

```bash
python3 scripts/verify_two_pair_primitive_prefix_obstruction.py
```

The dependency-free checker performs direct sparse contractions for
\(4\leq d\leq8\) through the first nonzero moment and verifies (1.4).  It
also verifies (1.6)--(1.7) for \(4\leq d\leq7\) and every
\(1\leq s\leq d\), and checks the three exact branch formulas (1.10).
It additionally checks (1.12)--(1.13) for \(4\leq d\leq7\).  The
checker verifies the three formulas (1.16) in the same degree range, the
even-height formulas (1.18) for all applicable pairs in that range, and
the exact degree-five certificate (1.21)--(1.24).  The phase-support
argument above proves the all-degree statements.
