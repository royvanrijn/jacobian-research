# Exact computational exclusion of the remaining \((72,108)\) configurations

**Status:** exact characteristic-zero computation completed; external mathematical review is still required before treating this as an established result.

**Date:** 2026-07-21

## 1. Scope of the result

Guccione–Guccione–Horruitiner–Valqui, Proposition 4.3 of arXiv:2204.14178, show that a counterexample in their remaining \((8,28)\) case—and hence in the unresolved degree pair \((72,108)\)—would produce Laurent polynomials \(P,Q\in K[x,x^{-1},y]\) satisfying

\[
[P,Q]=x^2
\]

with one of two pairs of Newton polygons:

\[
\begin{aligned}
\text{Case 1: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16),(0,8)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24),(0,12)\};
\end{aligned}
\]

or

\[
\begin{aligned}
\text{Case 2: }N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24)\}.
\end{aligned}
\]

The computations in this package give the unit ideal in characteristic zero for the complete coefficient systems of **both** cases. Subject to checking that Proposition 4.3 has been transcribed and normalized exactly, this excludes the remaining \((72,108)\) configuration.

This does **not** prove the two-dimensional Jacobian conjecture. Combined with the enumeration in the cited paper, it would raise the known lower bound for the maximum degree of a planar counterexample from \(108\) to \(125\).

## 2. Exact Laurent reduction

Set

\[
t=xy^2,\qquad z=y^{-1}.
\]

Then \([t,z]_{x,y}=-1\) and \(x^2=t^2z^4\). The upper bands of either polygon therefore have the form

\[
P=A(t)z^2+B(t)z+C(t),\qquad
Q=D(t)z^3+E(t)z^2+F(t)z+G(t).
\]

The coefficient identities of \([P,Q]=x^2\) are

\[
2AD'-3A'D=t^2, \tag{J4}
\]
\[
2(AE'-A'E)+(BD'-3B'D)=0, \tag{J3}
\]
\[
(2AF'-A'F)+(BE'-2B'E)-3C'D=0, \tag{J2}
\]
\[
2AG'+(BF'-B'F)-2C'E=0, \tag{J1}
\]
\[
BG'-C'F=0. \tag{J0}
\]

`verify_laurent_reduction.py` checks the chain rule, all five identities, and the lattice-point counts. `verify_case1_reduction.py` independently checks all bands in Case 1.

## 3. Normalization

Let the coefficients at the three nonzero vertices \((1,0),(8,14),(8,16)\) of \(P\) be \(a_1,a_8,c_8\). For

\[
\widetilde P=\rho P(\lambda x,\mu y),\qquad
\widetilde Q=\sigma Q(\lambda x,\mu y),
\]

we have

\[
[\widetilde P,\widetilde Q]=\rho\sigma\lambda^3\mu\,x^2.
\]

Because the ground field is algebraically closed and all five scaling parameters involved are nonzero, one can choose \(\lambda,\mu,\rho,\sigma\) so that

\[
a_1=a_8=c_8=1
\]

and the right-hand side remains \(x^2\). Explicitly, choose

\[
\mu^2=\frac{a_8}{c_8},\qquad
\lambda^7=\frac{a_1}{a_8}\mu^{-14},\qquad
\rho=(a_1\lambda)^{-1},\qquad
\sigma=(\rho\lambda^3\mu)^{-1}.
\]

Constants can be subtracted from \(P,Q\) without changing the bracket.

## 4. The exact first block

Write

\[
A=t+a_2t^2+\cdots+a_7t^7+t^8,
\qquad
D=d_2t^2+\cdots+d_{12}t^{12}.
\]

Equation (J4) solves all eleven coefficients of \(D\), and leaves six compatibility equations in \(a_2,\dots,a_7\). An exact rational Gröbner/FGLM calculation gives a lexicographic basis consisting of

- one irreducible degree-35 polynomial \(H(a_7)\), and
- five equations linear in \(a_2,\dots,a_6\).

Thus the first-block quotient is the number field

\[
K=\mathbb Q[u]/(H(u)),\qquad [K:\mathbb Q]=35.
\]

`verify_firstblock_exact.py` checks:

1. all eleven triangular substitutions for \(D\);
2. all six original compatibility equations vanish modulo the lexicographic basis;
3. \(H\) is irreducible over \(\mathbb Q\).

No finite-field assumption remains after this point.

## 5. Case 2: exact unit ideal

Equations (J3) and (J2) are solved by exact Gaussian elimination over \(K\). They leave three free parameters. Equation (J1), followed by (J0), gives 25 residual equations. Four of them already generate the unit ideal:

\[
I_2=(R_0,R_1,R_7,R_9)=K[r,s,h].
\]

Singular's exact algebraic-number-field routine returns

```text
INPUT_SIZE=4
SIZE=1
J[1]=1
UNIT
```

in `case2_compact4_exact.out`.

## 6. Case 1: exhaustive exact collapse

The extra negative \(z\)-bands are introduced in complete slices. Solving the first five new bracket layers gives 13 compatibility equations in six parameters

\[
r,s,h,u_1,u_2,u_3.
\]

A fresh regeneration of the checkpoint reproduced the original data byte-for-byte.

### 6.1 Exhaustive split in \(s\)

The first compatibility equation factors over \(K\) as

\[
\kappa(s-c)^2(s+c)^2,
\qquad \kappa,c\in K^\times.
\]

The second compatibility equation has the same two factors. Hence every solution lies in exactly one of the exhaustive branches

\[
s=c\quad\text{or}\quad s=-c.
\]

The factorization is checked both by Singular and by direct multiplication in `exact_core.K`.

### 6.2 Forced affine relation for \(r\)

On either branch, exact equations \(E_2,E_6\) satisfy

\[
E_6-\lambda E_2=S L^2,
\qquad S\in K^\times,
\]

where

\[
\lambda=\frac{13}{9}\quad\text{or}\quad -\frac{13}{9},
\]

and

\[
L=r+\alpha h+\frac13u_2+\gamma,
\qquad
\alpha=-\frac23\quad\text{or}\quad \frac23,
\qquad \gamma\in K.
\]

Therefore \(L=0\), with no division by a polynomial parameter. Substitution removes \(r\). Three remaining equations are exact constant multiples of another equation and are discarded without changing the zero set.

### 6.3 Invertible affine change in \(h\)

An affine translation of \(h\) makes the coefficient of \(u_3\) in the simplest remaining equation equal to a nonzero field scalar times the new \(h\). This is an invertible coordinate change—not a case assumption. Each branch is reduced to seven equations in

\[
h,u_1,u_2,u_3.
\]

### 6.4 Exact standard bases

For each branch, Singular returns the exact standard basis \(\{1\}\):

```text
# branch s=c
INPUT=7
SIZE=1
J[1]=1
UNIT

# branch s=-c
INPUT=7
SIZE=1
J[1]=1
UNIT
```

The fresh runs are stored in

- `case1_branch1_after_w_nfmod.out`, and
- `case1_branch2_after_w_nfmod.out`.

## 7. Why these are characteristic-zero certificates

The final rings are polynomial rings over the algebraic number field \(K=\mathbb Q[u]/(H)\). The routine `nfmodStd` uses modular images only as an algorithmic acceleration. It reconstructs a basis over characteristic zero and then performs an exact final test: the original generators must reduce to zero by the reconstructed basis, and the reconstructed basis is itself checked to be a standard basis. The returned element is literally \(1\) in the characteristic-zero ring.

Thus the final outputs are not merely repeated finite-field obstructions.

## 8. Conclusion and appropriate claim

The exact computation establishes:

> There are no Laurent-polynomial pairs satisfying Proposition 4.3, Case (1) or Case (2), with bracket \([P,Q]=x^2\) and the stated Newton polygons.

Consequently, **assuming the cited reduction is applied exactly and no transcription/normalization condition has been missed**, the last \((72,108)\) configuration is excluded. Together with Theorem 2.1 of arXiv:2204.14178, this would imply:

\[
\boxed{\text{Every planar Jacobian counterexample has }\max\{\deg P,\deg Q\}\ge 125.}
\]

It does not settle JC2. Before public announcement, the reduction, normalization, and computer algebra certificate should be independently checked by specialists and ideally by a second CAS.

## 9. Core reproduction commands

```bash
python verify_laurent_reduction.py
python verify_case1_reduction.py
python verify_firstblock_exact.py

python case2_exact_generate.py
Singular -q case2_compact4_exact.sing

python case1_cascade_analysis.py
python case1_cascade_w.py
Singular -q case1_branch1_after_w_nfmod.sing
Singular -q case1_branch2_after_w_nfmod.sing
```

A complete regeneration of Case 1 from the lattice slices is available through `case1_descent_checkpoint.py`; remove or rename `case1_checkpoint.pkl` first to force a fresh rebuild.
