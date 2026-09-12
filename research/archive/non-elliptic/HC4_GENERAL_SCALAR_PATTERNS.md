# All-degree patterns in the scalar reverse-Schur branch of `HC4`

## Status

The degree-eight and degree-nine computations reveal that the root-partition
census is not the natural parameterization. The same local mechanisms admit
all-degree formulas, and the relevant complexity is the **transverse excess**
rather than the total leading degree.

This note records four exact all-degree results and a research program for the
remaining scalar branch.

> **Theorem HC4RSD44 — resonance square identity and transverse-excess
> compression.** Let `f` be a degree-`d` binary leading form in the
> synchronized scalar reverse-Schur packet, let a root have multiplicity
> `m<d`, and let the first transverse coefficient `g` have degree `e<=d-2`.
> If `n=ord_r(g)`, then
> \[
> C_{d,m,e,n}=(dn-em)^2+m(d-m)(d-2e-1).             \tag{0.1}
> \]
> Hence proper-root resonances occur only when `2e+1>d`. Away from resonance,
> write
> \[
> f=A^2B,
> \]
> where `B` is the squarefree product of roots of odd multiplicity and put
> `o=deg B`. Then
> \[
> g=ABH,\qquad W:=\deg(AB)=\frac{d+o}{2},
> \]
> and, with `h=e-W`,
> \[
> q=BC,\qquad \deg H=h,\qquad \deg C=2h.           \tag{0.2}
> \]
> Thus the first Schur problem is controlled primarily by `h`, not by `d`.

> **Theorem HC4RSD45 — all-even minimal-excess closure.** If every root
> multiplicity is even and `h=0`, then
> \[
> f=A^2,\qquad g=aA,\qquad q=\frac{a^2}{d}.
> \]
> For
> \[
> c=A^2+azA+\frac{a^2}{2d}z^2
> \]
> one has
> \[
> [z^4]J(c)=\frac{a^6}{d^2}\det\operatorname{Hess}A.
> \]
> Hence `A` is a power of a linear form and the packet is a fixed cylinder.

> **Theorem HC4RSD46 — sparse-odd-support minimal-excess closure.** Let
> `h=0`, let `o=deg B>0`, and put
> \[
> k:=\deg A=\frac{d-o}{2}.
> \]
> If
> \[
> d>3o\qquad\text{equivalently}\qquad k>o,
> \tag{0.3}
> \]
> then every such scalar packet is a fixed cylinder. There is no degree
> bound. Indeed the complete same-weight face has no `z^3` term and is
> exactly
> \[
> c=A^2B+azAB+\frac{\kappa a^2}{2}z^2B.             \tag{0.4}
> \]
> If `\kappa\ne0` and `o>=2`, its highest coefficient is
> \[
> [z^6]J(c)=
> -\frac{2(o+2)}{o-1}a^8\left(\frac\kappa2\right)^4
> B^2\det\operatorname{Hess}B.                     \tag{0.5}
> \]
> Since `B` is squarefree, this is impossible. If `o=1`, normalize `B=y`.
> The `z^5` coefficient first forces `A_{xx}=0`; homogeneity then gives
> `A=y^{k-1}(ux+vy)`. If `u\ne0`, the two multiplicities are `(2k-1,2)` and
> the next coefficient is a nonzero scalar multiple of
> \[
> \frac{16(k-1)^3(k+1)}{(2k-1)(2k+1)}u^2,
> \]
> so `u=0` and the top is again a pure power. Finally, if `\kappa=0`, then
> \[
> [z^2]J(c)=
> -\frac{k+o+1}{k+o-1}a^4(AB)^2
> \det\operatorname{Hess}(AB),                     \tag{0.6}
> \]
> which forces `AB` to be a power of a linear form and again gives a fixed
> cylinder.

> **Theorem HC4RSD47 — Pell sieve for zero interpolation value.** In an
> `h=0` packet with at least one odd-multiplicity root, suppose the common
> interpolation value is `\kappa=0`. Then the odd-root count `o` must satisfy
> \[
> 8o(o+1)=s^2.
> \]
> Equivalently, writing `u=2o+1` and `s=2v`,
> \[
> u^2-2v^2=1.                                      \tag{0.7}
> \]
> Hence the only possible odd-root counts are the Pell sequence
> \[
> o=1,8,49,288,1681,\ldots .                       \tag{0.8}
> \]
> For such an `o`, every odd root multiplicity `m` must satisfy
> \[
> \frac md=
> \frac{3o+2\pm2v}{(o+2)^2},                       \tag{0.9}
> \]
> while every even root multiplicity, if present, must satisfy
> \[
> \frac md=\frac{4(o+1)}{(o+2)^2}.                 \tag{0.10}
> \]
> Thus the zero-value exceptional set is an explicitly arithmetic family,
> not an unconstrained root-partition locus.

These statements do **not** yet close all `h=0`: when `d<=3o`, complete
same-weight `z^3,z^4,...` tails can occur and must be retained.

## 1. Resonance as a square

The old expression

\[
 C_{d,m,e,n}
 =d^2m+d^2n^2-2demn-2dem-dm^2-dm
  +e^2m^2+2em^2+m^2
\]

completes to (0.1). Its discriminant as a quadratic in `n` is

\[
 \Delta_n=-4d^2m(d-m)(d-2e-1).
\]

The resonance equation is therefore

\[
 (dn-em)^2=m(d-m)(2e+1-d).                          \tag{1.1}
\]

For the minimal nonresonant order `n=ceil(m/2)`, the root interpolation value
is

\[
 \kappa_{d,e}(m)
 =\frac{2e+1-d}{d}
  -\frac{(d\lceil m/2\rceil-em)^2}{dm(d-m)}.        \tag{1.2}
\]

For two distinct even multiplicities,

\[
 \kappa(m_1)-\kappa(m_2)
 =-\frac{(d-2e)^2(m_1-m_2)}
         {4(d-m_1)(d-m_2)}.                         \tag{1.3}
\]

For two distinct odd multiplicities, equality requires

\[
 d(m_1+m_2-d)
 +m_1m_2(d-2e)(d-2e+2)=0.                          \tag{1.4}
\]

## 2. Transverse excess

If `o` roots have odd multiplicity, then

\[
 W=\frac{d+o}{2}.
\]

The first transverse coefficient vanishes when `e<W`; otherwise set
`h=e-W`. Then

\[
 g=ABH,\quad \deg H=h,
\qquad
 q=BC,\quad \deg C=2h.                              \tag{2.1}
\]

At every distinct root `r`,

\[
 C(r)=\kappa_{d,e}(m_r)H(r)^2.                     \tag{2.2}
\]

The `j`-th same-weight `z` tail has binary degree

\[
 \deg R_j=je-(j-1)d
 =d+\frac j2(o-d)+jh.                               \tag{2.3}
\]

For `h=0`, writing `d=2k+o`, this becomes

\[
 \deg R_j=o-(j-2)k.                                 \tag{2.4}
\]

This is why the region `k>o` is especially clean: every `j>=3` tail is
absent, so (0.4) is the **complete** face rather than a truncation.

## 3. Universal highest-coefficient identities in the tail-free `h=0` region

Put `s=\kappa/2`. For

\[
 c=B(A^2+azA+sa^2z^2)
\]

a direct symbolic differentiation gives

\[
 [z^6]J(c)=2Ba^8s^4
 \left(
 2B\det H_B-3\nabla B^{\mathsf T}\operatorname{adj}(H_B)\nabla B
 \right).                                          \tag{3.1}
\]

For a homogeneous binary form of degree `o>1`, Euler gives

\[
 \nabla B^{\mathsf T}\operatorname{adj}(H_B)\nabla B
 =\frac{o}{o-1}B\det H_B,
\]

and (3.1) becomes (0.5). A homogeneous binary form of zero Hessian is a power
of a linear form, incompatible with squarefree `B` unless `o=1`.

When `s=0`, a different coefficient survives. Put `P=AB`, of degree `k+o`.
The exact identity is (0.6). Thus even the resonant zero-value case reduces
to the binary zero-Hessian theorem.

## 4. Pell arithmetic of the zero-value case

At `h=0` one has `e=(d+o)/2`. For an odd multiplicity `m`, the equation
`\kappa(m)=0` becomes, with `x=m/d`,

\[
 (o+2)^2x^2-(6o+4)x+1=0.                            \tag{4.1}
\]

Its discriminant is `32o(o+1)`. Rational root multiplicities therefore
require `8o(o+1)` to be a square. Setting `u=2o+1` converts this to the Pell
equation (0.7), yielding (0.8). The two odd-root ratios are (0.9).
For even `m`, `\kappa(m)=0` is linear and gives (0.10).

The first two Pell families explain exact packets already seen computationally:

* `o=1`: degrees are multiples of `9`, with the basic multiplicity ratio
  `(8,1)`;
* `o=8`: the basic normalized ratios are `25`, seven copies of `1`, and an
  even multiplicity `18`, summing to `50`.

The next Pell value is `o=49`, where several integral count combinations are
possible. These are now generated arithmetically rather than discovered by a
degree census.

## 5. What remains

The scalar strategy is now naturally two-dimensional:

1. **minimal excess with tails:** `h=0` and `d<=3o`;
2. **first positive excess:** `h=1`, a weighted Veronese interpolation problem.

The first should be attacked by classifying the finite possible tail lengths
from (2.4), not by total degree. In particular:

* `d=3o` means `k=o` and permits only a scalar `z^3` tail;
* `2o<d<3o` still permits only `z^3`, of degree `o-k<o`;
* smaller `d/o` ratios permit longer tails but force many odd roots, where
  the existing squarefree and one-double-root all-degree theorems already
  remove large open strata.

For `h=1`, write

\[
 H=uX+vY,
 \qquad C=c_0X^2+c_1XY+c_2Y^2.
\]

Then (2.2) is a weighted Veronese interpolation problem. A determinant
criterion in terms of root cross-ratios and the multiplicity weights
`\kappa(m_r)` should replace the old chart-by-chart Gröbner calculations.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_general_scalar_patterns.py
```

The checker proves symbolically the resonance square identity, parity
coincidence formulas, transverse-excess bookkeeping, the all-even terminal
coefficient, both tail-free `h=0` highest-coefficient identities, and the Pell
reduction. No fixed total degree is used.
