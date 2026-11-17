# Minimal-excess tail descent in the scalar `HC4` packet

## Status

This note continues `HC4RSD44--47` from
[`HC4_GENERAL_SCALAR_PATTERNS.md`](HC4_GENERAL_SCALAR_PATTERNS.md).
It replaces another large root-partition range by identities that are
independent of the total leading degree.

Let the synchronized scalar face have transverse excess `h=0`. Write

\[
 f=A^2B,
 \qquad o=\deg B,
 \qquad k=\deg A=\frac{d-o}{2},
 \qquad e=k+o,
\]

where `B` is the squarefree odd-multiplicity factor. The `j`-th same-weight
coefficient has binary degree

\[
 \deg R_j=d-jk=o-(j-2)k.                            \tag{0.1}
\]

> **Theorem HC4RSD48 — universal highest-tail and mixed-tail identities.**
> Let a weighted face satisfying the bordered equation have highest nonzero
> `z`-tail
> \[
> \frac{z^j}{j!}R(x,y),\qquad j\ge3,
> \]
> with homogeneous `R` of degree `r>0`. Ignoring the nonzero factorial scale,
> the absolute highest `z` coefficient is
> \[
> jR\left(jR\det H_R-(j+1)\nabla R^{\mathsf T}
> \operatorname{adj}(H_R)\nabla R\right).
> \tag{0.2}
> \]
> By Euler homogeneity this equals
> \[
> -\frac{j(j+r)}{r-1}R^2\det H_R                  \tag{0.3}
> \]
> for `r>1`. Hence the highest tail is a power of a linear form; the case
> `r=1` already has that form.
>
> Normalize `R=x^r`. Let `S z^i/i!`, `i<j`, be the next nonzero lower tail.
> The highest coefficient linear in `S` is, again up to the nonzero factorial
> scale,
> \[
> -jr(r+j)x^{3r-2}S_{yy}.                          \tag{0.4}
> \]
> Remarkably, the index `i` cancels completely. Thus every next tail satisfies
> `S_{yy}=0`, so
> \[
> S=x^{\deg S-1}(\alpha x+\beta y).                \tag{0.5}
> \]

> **Theorem HC4RSD49 — complete minimal-excess closure for `d>=2o`.**
> Every `h=0` synchronized scalar packet satisfying
> \[
> d\ge2o                                                   \tag{0.6}
> \]
> has a fixed ruling and reduces to `HC2` or the exact `JC2` cotangent
> endpoint. There is no total-degree bound.

The proof of HC4RSD49 splits only by the number of same-weight tails, not by
root partitions.

## 1. Tail-free region `d>3o`

This is `HC4RSD46`: `k>o`, so (0.1) permits no `j>=3` tail. The complete face
is

\[
 c=A^2B+azAB+\frac{\kappa a^2}{2}z^2B.             \tag{1.1}
\]

For `\kappa\ne0`, its `z^6` coefficient forces `det Hess B=0`; squarefreeness
then leaves only `o=1`, where the next two coefficients force the top to be a
pure power. For `\kappa=0`, the `z^2` coefficient forces
`det Hess(AB)=0`, again giving a fixed cylinder.

## 2. One-tail region `2o<d<=3o`

Now (0.1) permits at most

\[
 \frac{z^3}{6}R,
 \qquad r:=\deg R=o-k\ge0.                          \tag{2.1}
\]

At `d=3o`, `r=0` and `R` is scalar. If `R\ne0`, exact differentiation gives

\[
 [z^8]J(c)=\text{const}\cdot \kappa^2R^2\det H_B
 \tag{2.2}
\]

when `\kappa\ne0`, and

\[
 [z^6]J(c)=\text{const}\cdot R^2\det H_{AB}
 \tag{2.3}
\]

when `\kappa=0`. Thus `R=0` unless the packet is already a cylinder, and
(1.1) finishes the row.

For `2o<d<3o`, one has `r>0`. HC4RSD48 makes

\[
 R=x^r.
\]

If `\kappa\ne0`, the next exact coefficient is

\[
 [z^9]J(c)=
 \text{const}\cdot \kappa\,r^2(7-3r)x^{3r-2}B_{yy}. \tag{2.4}
\]

The integer factor `7-3r` is never zero. Hence `B_{yy}=0`, so

\[
 B=x^{o-1}(\alpha x+\beta y),
\]

contradicting squarefreeness because this strip has `o>=3`.

If `\kappa=0`, the corresponding coefficient is

\[
 [z^8]J(c)=
 \text{const}\cdot r^2(7-3r)x^{3r-2}(AB)_{yy}.     \tag{2.5}
\]

Then `AB` has at most two distinct projective roots, impossible because the
squarefree divisor `B` has degree `o>=3`. Thus the `z^3` tail vanishes and the
tail-free identities finish the packet.

## 3. Boundary `d=2o`

Here `k=o/2`; admissibility `e<=d-2` forces `o>=4`. The complete additional
tails are

\[
 \frac{z^3}{6}R_k+\frac{z^4}{24}t.                 \tag{3.1}
\]

If `t=0`, the one-tail argument above works verbatim. Suppose `t\ne0`. The
highest coefficient is

\[
 [z^{12}]J(c)=\frac{t^2}{1296}\det H_R.            \tag{3.2}
\]

Thus `R=x^k`. Since `k>=2`, the next coefficient is

\[
 [z^{11}]J(c)=
 \text{const}\cdot \kappa\,t^2k(k-1)x^{k-2}B_{yy}
 \tag{3.3}
\]

if `\kappa\ne0`, and

\[
 [z^{10}]J(c)=
 \text{const}\cdot t^2k(k-1)x^{k-2}(AB)_{yy}      \tag{3.4}
\]

if `\kappa=0`. Again squarefreeness of the degree-`o>=4` factor `B` makes
either alternative impossible. Hence `t=0`, reducing to the one-tail case.

Sections 1--3 prove HC4RSD49.

## 4. Why this looks inductive

Equation (0.4) is the important new feature. It does not know how far the next
tail lies below the highest one. Every top tail first becomes a pure power,
and then every first descendant becomes at most linear in the transverse
projective coordinate. The remaining range `d<2o` merely allows a longer tail
chain:

\[
 R_3,R_4,R_5,\ldots .
\]

The next target should therefore be a **tail-chain recurrence theorem** rather
than another `d`-census: compute the first coefficient where a quadratic
combination of an already constrained descendant can meet the linear term of
the next descendant, and show that the recurrence forces all tails to share a
single ruling.

## 5. Geometric interpretation

The equation used throughout is the vanishing determinant of the bordered
Hessian, also known as the Universal Field Equation. In three variables its
regular level surfaces are developable. The tail descent above is therefore
consistent with the classical developable-surface picture: the highest
weighted direction produces a ruling and subsequent coefficients are forced
toward the same ruling. This suggests that a classification of polynomial or
weighted-homogeneous solutions of the Universal Field Equation may provide a
conceptual proof of the full scalar packet.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_minimal_excess_tail_descent.py
```

The checker verifies (0.2)--(0.4), the one-tail coefficients (2.2)--(2.5), and
the boundary coefficients (3.2)--(3.4) symbolically with abstract jets. No
root partition and no fixed total degree is used.
