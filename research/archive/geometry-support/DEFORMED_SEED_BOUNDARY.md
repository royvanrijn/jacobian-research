# Deformed-seed boundary and one-extra-root theorem

This note extends the canonical-family image calculation to admissible
primitives with additional zeros. The first part gives the uniform boundary
mechanism. The second completely solves the family with one additional simple
zero.

## Boundary-clean primitives

Let `H` be an admissible characteristic-zero primitive of degree `n>=3`, with

\[
H(W)=W^m(W-1)R(W),\qquad m\ge2,
\]

where `R(0)R(1)!=0` and every zero of `R` is simple. Write `p=H'`, retain the
weighted normalization constant `c!=0`, and assume `p'(1)/c!=-2` as in the
general construction.

For a target `(A,B,C)`, the inverse polynomial is

\[
E(W)=H(W)-BCW+cAC^2.
\]

On `C!=0`, its simple roots are exactly the affine source points and its
repeated roots are exactly the reconstruction poles.

### Direct `C=0` charts

The chart `x=0` remains triangular and maps bijectively onto the entire target
plane `C=0`.

The only other finite chart has `gamma=0`, hence `W=0`. If

\[
H(W)=h_2W^2+O(W^3),\qquad k=h_2/c\ne0,
\]

then on this chart

\[
A={u+ku^2\over x^2},\qquad
B={c(1+2ku)\over x}.
\]

Eliminating `u` gives

\[
\left({B^2\over c^2}-4kA\right)x^2=1.
\]

Thus a double zero at `W=0` contributes two points off this conic and none on
it. Together with the `x=0` point, the `C=0` fiber sizes are `3` and `1`.

If `m>=3`, the higher terms vanish on `gamma=0` and

\[
A={u\over x^2},\qquad B={c\over x}.
\]

This contributes one point when `B!=0` and none when `B=0`, giving total
fiber sizes `2` and `1`.

Every root of `R` is a projective inverse branch but not a finite `C=0` source
point. Indeed, for a simple extra root `rho`, implicit continuation gives a
root `W(C)->rho`, while

\[
\gamma\longrightarrow-{H'(\rho)\over c}\ne0,
\qquad x={C\over\gamma}\longrightarrow0.
\]

Since the bounded `x=0` chart has `W->1`, a branch with `rho!=1` must escape in
`y` or `z`. Extra primitive zeros therefore create genuine boundary branches;
they are not denominator-lost affine points.

### General image and nonproperness schema

Let

\[
D_H(s,t)=\operatorname{Disc}_W(H(W)-sW+t).
\]

Under the boundary-clean hypotheses, the pullback has exact order `m`:

\[
D_H(BC,cAC^2)=C^mQ_H(A,B,C).
\]

For `m=2`, the trace `Q_H(A,B,0)` is a nonzero scalar multiple of the conic
above. For `m>=3`, it is a nonzero scalar multiple of `B^m`.

Let `Omega_H` be the finite set of `(s,t)` for which `H(W)-sW+t` has no simple
root. Since every `C=0` target is supplied by the `x=0` chart,

\[
G_H(\mathbb C^3)=\mathbb C^3\setminus
\{(A,B,C):C\ne0,\ (BC,cAC^2)\in\Omega_H\}.
\]

Except for the minimal cubic case (`n=3`, `m=2`, `R` constant), the direct
`C=0` fiber contains fewer than the generic `n` sheets. Étaleness prevents
distinct bounded branches from merging, so at least one branch escapes over
every target on `C=0`. Consequently, for every boundary-clean inverse degree
`n>=4`,

\[
S_{G_H}=V(C)\cup V(Q_H).
\]

This statement explicitly accounts for both the discriminant saturation and
all primitive-zero boundary branches.

## One additional simple zero

Every normalized primitive with a double zero at `0`, the distinguished zero
at `1`, and one additional simple zero `rho` has the form

\[
H_\rho(W)={c\over1-\rho}W^2(1-W)(W-\rho),
\]

where

\[
\rho\notin\{0,1,2\}.
\]

The first two exclusions preserve the stated zero profile. The last is exactly
the forbidden construction value `p'(1)/c=-2`. Here

\[
{p'(1)\over c}=-2{2\rho-3\over\rho-1},
\qquad
a=-{3\rho-5\over2(\rho-2)},
\qquad
k={\rho\over\rho-1}.
\]

The `C=0` fiber therefore has three points off

\[
{B^2\over c^2}-4{\rho\over\rho-1}A=0
\]

and one point on it. The additional root `rho` supplies an escaping fourth
branch everywhere on the target plane, so `V(C)` is a nonproperness component
even off this conic.

### Exact omitted locus

A quartic with no simple root must be a square of a quadratic. Coefficient
comparison gives exactly one such inverse value:

\[
s_\rho={c(1-\rho^2)\over8},
\qquad
t_\rho=-{c(1-\rho)^3\over64},
\]

and the exact factorization

\[
H_\rho(W)-s_\rho W+t_\rho
=-{c\over1-\rho}
\left(
W^2-{1+\rho\over2}W-{(\rho-1)^2\over8}
\right)^2.
\]

It is generically a double-double node. The quadratic itself becomes double
exactly when

\[
3\rho^2-2\rho+3=0,
\]

giving a quadruple-root degeneration.

Thus the exact affine image is

\[
G_\rho(\mathbb C^3)=\mathbb C^3\setminus
V\left(
BC-{c(1-\rho^2)\over8},
cAC^2+{c(1-\rho)^3\over64}
\right),
\]

and, writing

\[
D_{H_\rho}(BC,cAC^2)=C^2Q_\rho(A,B,C),
\]

the nonproperness set is

\[
S_{G_\rho}=V(C)\cup V(Q_\rho).
\]

The earlier exceptional quartic is the specialization `rho=-1`, up to its
chosen coordinate normalization.

Run:

```bash
.venv/bin/python scripts/verify_deformed_seed_boundary.py
```

The verifier checks the symbolic `rho`-family, the direct boundary conic,
exact `C^2` saturation, the omitted square factorization, and four rational
deformations of the actual polynomial map.

The omitted-value problem for several or repeated extra zeros is solved by the
exact multiplicity-partition classifier in `OMITTED_VALUE_CLASSIFICATION.md`.
Repeated-root saturation and its additional boundary multiplicities are proved
in `REPEATED_ROOT_BOUNDARY.md`.
