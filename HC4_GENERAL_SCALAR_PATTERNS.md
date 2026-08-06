# All-degree patterns in the scalar reverse-Schur branch of `HC4`

## Status

The degree-eight and degree-nine computations reveal that the root-partition
census is not the natural parameterization.  The same local mechanisms admit
all-degree formulas, and the relevant complexity is the **transverse excess**
rather than the total leading degree.

This note records two exact all-degree results and a research program for the
remaining scalar branch.

> **Theorem HC4RSD44 — resonance square identity and transverse-excess
> compression.**  Let `f` be a degree-`d` binary leading form in the
> synchronized scalar reverse-Schur packet, let a root have multiplicity
> `m<d`, and let the first transverse coefficient `g` have degree `e<=d-2`.
> If `n=ord_r(g)`, then the root-valuation resonance polynomial of
> `HC4RSD27` satisfies
>
> \[
> C_{d,m,e,n}
> =(dn-em)^2+m(d-m)(d-2e-1).                       \tag{0.1}
> \]
>
> Hence its discriminant as a quadratic in `n` is
>
> \[
> \Delta_n=-4d^2m(d-m)(d-2e-1).                   \tag{0.2}
> \]
>
> In particular no proper-root resonance is possible when
> `2e+1<=d`; resonances can occur only in the upper half of the transverse
> degrees.
>
> Away from resonance, write
>
> \[
> f=A^2B,
> \]
>
> where `B` is the squarefree product of the roots having odd multiplicity,
> and put `o=deg B`.  The valuation bound gives
>
> \[
> g=ABH,
> \qquad
> W:=\deg(AB)=\sum_r\left\lceil\frac{m_r}{2}\right\rceil
>            =\frac{d+o}{2}.                       \tag{0.3}
> \]
>
> If `h=e-W`, then every polynomial Schur solution has
>
> \[
> q=BC,
> \qquad \deg H=h,
> \qquad \deg C=2h.                                \tag{0.4}
> \]
>
> Thus the first Schur problem is a degree-`h` / degree-`2h` weighted-square
> interpolation problem.  Its algebraic complexity is controlled primarily
> by `h`, not by `d`.

> **Theorem HC4RSD45 — all-even minimal-excess closure.**  In HC4RSD44,
> suppose every root multiplicity of `f` is even and `h=0`.  Then `d=2k`,
> `f=A^2`, and after scaling the unique nonzero first transverse direction is
>
> \[
> g=aA,
> \qquad q=\frac{a^2}{d}.                           \tag{0.5}
> \]
>
> For the complete face
>
> \[
> c=A^2+azA+\frac{a^2}{2d}z^2,                     \tag{0.6}
> \]
>
> the full bordered polynomial satisfies
>
> \[
> [z^4]J(c)
> =\frac{a^6}{d^2}\det\operatorname{Hess}_{x,y}A. \tag{0.7}
> \]
>
> Therefore a nonzero transverse direction forces
> `det Hess(A)=0`.  Over characteristic zero a nonzero homogeneous binary
> form with zero Hessian is a power of a linear form.  Hence `A=ell^k`, the
> packet is a fixed cylinder, and it reduces to `HC2` or the exact `JC2`
> cotangent endpoint.  This closes the all-even `h=0` family in every even
> degree at once.

## 1. Why the resonance identity matters

The old expression

\[
 C_{d,m,e,n}
 =d^2m+d^2n^2-2demn-2dem-dm^2-dm
  +e^2m^2+2em^2+m^2
\]

looked degree-specific.  Completing the square gives (0.1) immediately.
The resonance equation is therefore the Diophantine conic

\[
 (dn-em)^2=m(d-m)(2e+1-d).                          \tag{1.1}
\]

This has several useful consequences.

* The entire lower half `2e+1<=d` is resonance-free.
* Exceptional local directions can be generated arithmetically before any
  Groebner calculation is attempted.
* The degree-nine `(m,e,n)=(3,5,1)` event is simply the first small solution
  of (1.1), not an isolated accident.

For the minimal nonresonant order `n=ceil(m/2)`, the root interpolation value
can also be written as

\[
 \kappa_{d,e}(m)
 =\frac{2e+1-d}{d}
  -\frac{(d\lceil m/2\rceil-em)^2}{dm(d-m)}.        \tag{1.2}
\]

Separating parity gives especially rigid coincidence formulas.  For two
even multiplicities `m1 != m2`,

\[
 \kappa(m_1)-\kappa(m_2)
 =-\frac{(d-2e)^2(m_1-m_2)}
         {4(d-m_1)(d-m_2)}.                         \tag{1.3}
\]

Thus distinct even multiplicities never have the same local value unless
`e=d/2`.  For two odd multiplicities,

\[
 \kappa(m_1)=\kappa(m_2)
\]

with `m1 != m2` only if

\[
 d(m_1+m_2-d)
 +m_1m_2(d-2e)(d-2e+2)=0.                          \tag{1.4}
\]

This explains the repeated phenomenon in degree nine where almost every
partition dies simply because two roots demand different constants.

## 2. Transverse excess instead of total degree

Let `o` be the number of odd-multiplicity roots.  Since

\[
 W=\frac{d+o}{2},
\]

the first transverse coefficient vanishes immediately when `e<W`.
Otherwise set

\[
 h=e-W.
\]

Then

\[
 g=ABH,\qquad \deg H=h,
\]

and the Schur coefficient has degree

\[
 \deg q=2e-d=o+2h.
\]

The parity factor `B` accounts for exactly the `o` forced roots, leaving

\[
 q=BC,\qquad \deg C=2h.
\]

At every distinct root `r`, the first local equation becomes

\[
 C(r)=\kappa_{d,e}(m_r)H(r)^2.                     \tag{2.1}
\]

This reorganizes the previous calculations:

* `h=0`: `H` and `C` are constants.  A nonzero packet requires all relevant
  root values `kappa` to coincide.  This is why parity/value comparisons
  killed most degree-eight and degree-nine partitions immediately.
* `h=1`: `H` is linear and `C` quadratic.  Exceptional sets are projective
  cross-ratio loci; the old Groebner systems were solving this same universal
  interpolation geometry repeatedly.
* `h=2`: `H` is quadratic and `C` quartic.  The first genuinely larger
  interpolation problem appears here, independently of how large `d` is.

The `j`-th same-weight `z` tail has binary degree

\[
 \deg R_j
 =je-(j-1)d
 =d+\frac j2(o-d)+jh.                              \tag{2.2}
\]

So even the size of the **complete** face is determined by `(d,o,h)` and is
usually small when `h` is small.

## 3. The all-even `h=0` family

When `o=0` and `h=0`, one has `e=d/2`.  Formula (1.3) degenerates exactly as
expected: all even roots demand the same value.  The Schur equation therefore
has the universal ray

\[
 f=A^2,\qquad g=aA,\qquad q=a^2/d.
\]

The important point is that this apparent survivor is killed without looking
at the root partition of `A`.  Differentiating (0.6) abstractly gives

\[
 [z^4]J(c)=\frac{a^6}{d^2}
 (A_{xx}A_{yy}-A_{xy}^2).
\]

Thus the old four-double-root obstruction was a special case of a general
binary-Hessian obstruction.  No cross-ratio census is necessary.

## 4. Research program suggested by the compression

The next useful target is **not degree ten**.  It is the classification by
small transverse excess.

### `h=0`

The remaining problem is purely arithmetic/projective: classify multiplicity
multisets for which all required values in (1.2) coincide, then derive a
partition-free complete-face obstruction analogous to HC4RSD45.  Equations
(1.3)--(1.4) already make this much smaller than a partition census.

### `h=1`

Write

\[
 H=uX+vY,
 \qquad C=c_0X^2+c_1XY+c_2Y^2.
\]

Equation (2.1) is a weighted Veronese interpolation problem: the evaluations
of `C` must equal prescribed weights times the square of one linear form.
The natural invariant is the cross-ratio of roots together with the finite
set of multiplicity weights `kappa`.  A useful goal is a determinant criterion
for solvability which replaces all chart-by-chart Groebner calculations.

### Resonant roots

Equation (1.1) should be treated first as a number-theoretic sieve.  For fixed
transverse defect `d-e`, its solutions may admit a finite parametrization in
`m,n`; this could isolate all possible resonance families before any
polynomial geometry enters.

## 5. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_general_scalar_patterns.py
```

The checker proves symbolically the square identity, its discriminant, the
parity coincidence formulas, the transverse-excess degree bookkeeping, and
the all-even terminal coefficient (0.7).  No fixed total degree is used.
