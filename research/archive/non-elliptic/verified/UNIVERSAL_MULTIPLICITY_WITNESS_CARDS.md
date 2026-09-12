# Universal multiplicity witness cards in degrees four, five, and six

The universal multiplicity theorem is not merely asymptotic.  This note gives
three small connected fields and, for each field, three exact Keller-map
presentations with pairwise distinct stable invariants.  The degree-four card
uses weighted maps; the degree-five and degree-six cards use quadratic
gauges.

All maps below are determinant one and have geometric degree equal to the
degree of the displayed field.  Every target fiber is complete.

## 1. The two map templates

### Weighted template

For a normalized weighted seed `H` with

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1,
\]

put

\[
 \kappa=H''(1),\qquad
 a_0=-\frac{1+\kappa}{2+\kappa}.
\]

With

\[
 v=xy,\quad S=x^2z,\quad
 \gamma=1+a_0v+S,\quad W=(1+v)\gamma,\quad C=x\gamma,
\]

the weighted map `F_H=(A,B,C)` is determined by

\[
 BC=H'(W)+\gamma,\qquad
 AC^2=W(H'(W)+\gamma)-H(W).                            \tag{1.1}
\]

The weighted polynomiality theorem says that the displayed quotients are
polynomials, and `det DF_H=1`.  Its inverse equation at `(A,B,C)` is

\[
 E(W)=H(W)-BCW+AC^2.                                  \tag{1.2}
\]

### Quadratic-gauge template

For

\[
 G(S)=g_1S+\cdots+g_NS^N,\qquad g_1g_3\ne0,
\]

put

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t).
\]

The determinant-one map `F_G=(Pi,U,C)` is

\[
\begin{aligned}
\Pi={}&tq,\\
U={}&-\frac12\left(
y+3\frac{g_3}{g_1}xq
+2\frac{g_2}{g_1}tq
+\sum_{k=4}^Nk\frac{g_k}{g_1}t^2x^{k-2}q^k
\right),\\
C={}&x(5-3t)-\frac{g_3}{g_1}x^3z
-\sum_{k=4}^N(k-2)\frac{g_k}{g_1}(xq)^k .
\end{aligned}                                         \tag{1.3}
\]

At a target `(1,0,c)`, its inverse equation is

\[
 E(S)=G(S)-\frac{g_1}{2}c.                            \tag{1.4}
\]

Thus every row below specifies an actual polynomial map by substituting its
coefficient vector into (1.3).

## 2. Degree four

Let

\[
 P_4(T)=T^4-3T^2-1.
\]

It is irreducible modulo `7`, and

\[
 \operatorname{Disc}(P_4)=-2704.
\]

Hence

\[
 A_4=\mathbb Q[T]/(P_4)
\]

is a connected quartic finite etale algebra.  Its trace-zero generator has
`Tr(T^2)=6`.  For each row choose `(e,u)` satisfying

\[
 e^2+2u^2=3,
\]

put `d=e-2u`, and define

\[
 P(W)=-\frac{P_4(u+dW)}{2d^3e},\qquad
 H(W)=P(W)-P(0)-P'(0)W.
\]

The exact rows are:

| `e` | `u` | `d` | `alpha=u/e-1/2` | weighted seed `H_alpha` | target `(A,B,C)` |
|---:|---:|---:|---:|---|---|
| `5/3` | `-1/3` | `7/3` | `-7/10` | `-W^2(W-1)(7W+3)/10` | `(107/3430,5/49,1)` |
| `5/3` | `1/3` | `1` | `-3/10` | `-W^2(W-1)(3W+7)/10` | `(107/270,-5/9,1)` |
| `-1` | `-1` | `1` | `1/2` | `W^2(W-1)(W-3)/2` | `(-3/2,-1,1)` |

Here

\[
 H_\alpha(W)=W^2(W-1)(\alpha W-\alpha-1).
\]

The corresponding weighted parameters `a_0` are respectively

\[
 -\frac{22}{17},\qquad-\frac{18}{13},\qquad-2.
\]

Equation (1.2) at each displayed target is exactly the corresponding
nonzero scalar multiple of `P_4(u+dW)`.  Thus all three complete fibers are
`\operatorname{Spec}A_4`.

All three seeds lie in the exact-double, Hessian-clean, boundary-clean
weighted locus.  Their distinct normalized parameters

\[
 -\frac7{10},\qquad-\frac3{10},\qquad\frac12
\]

are separated by weighted selected-root Torelli, so the three maps are
pairwise stably inequivalent.

## 3. Degree five

Let

\[
 P_5(T)=T^5+T^3+1.
\]

It is irreducible modulo `2`, with discriminant `3233`.  Put

\[
 A_5=\mathbb Q[T]/(P_5).
\]

For `s=1,2,3`, use

\[
 G_s(S)=P_5(s+S)-P_5(s).
\]

| `s` | `(g_1,g_2,g_3,g_4,g_5)` | target `(Pi,U,C)` | `I=a_5^5/(a_3a_4^6)` |
|---:|---|---|---:|
| `1` | `(8,13,11,5,1)` | `(1,0,-3/4)` | `64/171875` |
| `2` | `(92,86,41,10,1)` | `(1,0,-41/46)` | `529/2562500` |
| `3` | `(432,279,91,15,1)` | `(1,0,-271/216)` | `256/1421875` |

In every row, equation (1.4) is `P_5(s+S)`.  Translation identifies its
quotient algebra with `A_5`.  Every coefficient needed for the clean torus
is nonzero, and the three values of the complete quintic stable invariant
`I` are distinct.  The three determinant-one maps are therefore pairwise
stably inequivalent and have the same connected complete fiber.

## 4. Degree six

Let

\[
 P_6(T)=T^6+T^4+1.
\]

It is irreducible modulo `7`, with discriminant `-61504`.  Put

\[
 A_6=\mathbb Q[T]/(P_6).
\]

Again set `G_s(S)=P_6(s+S)-P_6(s)`.

| `s` | `(g_1,g_2,g_3,g_4,g_5,g_6)` | target `(Pi,U,C)` | `J_6=a_4a_6/a_5^2` |
|---:|---|---|---:|
| `1` | `(10,21,24,16,6,1)` | `(1,0,-3/5)` | `4/9` |
| `2` | `(224,264,168,61,12,1)` | `(1,0,-81/112)` | `61/144` |
| `3` | `(1566,1269,552,136,18,1)` | `(1,0,-811/783)` | `34/81` |

The inverse equations are `P_6(s+S)`, so every complete fiber is
`\operatorname{Spec}A_6`.  The displayed coefficient vectors lie in the
clean torus and the three `J_6` values are distinct.  Hence the three maps
are pairwise stably inequivalent.

This degree is the first one where the uniform top-three invariant appears:

\[
 J_6(s)=\frac5{12}+\frac1{36s^2}.
\]

The table checks it without suppressing any lower coefficient or exceptional
translation.

## 5. Exact regression

Run

```bash
.venv/bin/python scripts/verify_universal_multiplicity_witness_cards.py
```

The checker proves the three modular irreducibility certificates,
discriminants, trace-chord equations, weighted normal forms, exact targets,
translated inverse identities, clean-torus conditions, and pairwise
distinct stable invariants.
