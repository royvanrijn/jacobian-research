# Dessin-first closure of the no-vertical `(72,108)` Laurent branch

## Result and scope

Inside the two corrected Proposition-4.3 Laurent systems certified by
`PJ43A1`, the no-vertical-edge system has no characteristic-zero solution.
This is an alternative exact closure of that one Laurent branch.  It depends
on the complete intrinsic first-block graph certified by `PJ72Q1`; it is not a
stand-alone proof of the published degree reduction, the classification of
the two Laurent polygons, or the planar Jacobian conjecture.

The implementation follows the finite route suggested by the ramification
geometry:

1. enumerate the five degree-21 dessins;
2. reconstruct all five conjugate Belyi maps in one exact quintic model;
3. compile `(B,E)`, `(C,F)`, and `G` as kernels and cokernels of explicit
   linear maps over that quintic field;
4. decide only the resulting five-parameter terminal compatibility ideal.

No Gröbner basis of the original bivariate coefficient system is formed.

## 1. The five dessins

For

\[
\beta(X)=\frac{D(X)^2}{A(X)^3},\qquad
2AD'-3A'D=X^2,
\]

cancel the common powers `A=Xa`, `D=X^2d` and rescale the target so that the
third branch value is one.  The passport is

\[
(2^{10},1),\qquad (3^7),\qquad (17,1,1,1,1).
\]

Fixing the 17-cycle and solving the residual matching problem gives `85`
labelled center sets.  Quotienting by the centralizer leaves exactly

\[
(0,3,7,11),\ (0,3,7,12),\ (0,3,8,11),\
(0,3,8,13),\ (0,3,9,13).
\]

The checker emits a branch-cycle triple for every representative, checks the
product relation and transitivity, and computes trivial automorphism group in
all five cases.

## 2. Exact reconstruction and arithmetic action

Use the `mu_7`-quotient normalization from `PJ72Q1`:

\[
A=X+x_2X^2+x_3X^3+x_4X^4+x_5X^5+x_6X^6+qX^7+qX^8.
\]

The complete first-block graph presents `x_2,...,x_6` over

\[
K=\mathbf Q[q]/(H_5(q)),
\]

where the exact coefficients of `H_5` and of the graph are pinned in
`plane-jc/cas/firstblock_mu7_quotient_lex_basis.txt`.  The triangular
Wronskian recurrence then reconstructs all coefficients of `D` in `K`.
The generated certificate records every coefficient of `A` and `D` in the
basis `1,q,...,q^4`.

Put

\[
c=\frac{d_{10}^2}{a_7^3},\qquad
\widetilde\beta=\frac{X d^2}{c a^3}.
\]

Exact Euclidean arithmetic over `K` proves

\[
\deg(Xd^2-ca^3)=4,
\]

and proves that `a`, `d`, and `Xd^2-ca^3` are squarefree and pairwise
coprime.  Hence the three fibers have exactly the displayed passport; the
degree-four numerator gives multiplicity `17` at infinity.

The quintic `H_5` is irreducible over `Q`, has signature `(1,2)`, and its
Galois closure has group `S_5`.  Thus its five embeddings form one transitive
Galois orbit of normalized Belyi maps.  Since the five dessins have trivial
automorphism group, the five exact conjugates match the five combinatorial
types.  The arithmetic shortcut does **not** discard a type: the simple point
`X=0`, the unique index-17 point `X=infinity`, and the polynomial
normalization are all defined over the corresponding quintic field and are
preserved under conjugation.  In particular, demanding that the distinguished
point be rational over `Q` would be an unjustified stronger condition.

## 3. The three linear maps

Use coefficient bases compatible with the Newton windows and the common
target basis `1,X,...,X^19`.  The exact maps and their ranks are:

| stage | source | linear map | shape/rank | kernel |
|---|---|---|---|---|
| `(B,E)` | `B:X^1..X^8`, `E:X^2..X^12` | `-2AE'-BD'+3DB'+2EA'` | `20 x 19`, rank `17` | dimension `2` |
| `(C,F)` | `C:X^0..X^8`, `F:X^1..X^12` | `3DC'-2AF'+A'F` | `20 x 21`, rank `18` | dimension `3` |
| `G'` | `G:X^1..X^12` | `-2AG'` | `20 x 12`, rank `12` | zero |

The omitted constant of `G` is the one-dimensional target-translation
kernel.  The checker also verifies closed-form complete kernel bases.  For the
first map they are the density-preserving modes

\[
\phi=(A/X)^2-1,\qquad \phi=X,
\]

inserted in

\[
B=\phi(A'-A/X)-\tfrac12\phi'A,
\quad
E=\phi(D'-3D/(2X))-\tfrac34\phi'D.
\]

For the homogeneous second map the complete modes are

\[
f=D/X^2,\qquad f=1,\qquad f=X,
\]

inserted in

\[
C=fA'-(2f'/3+4f/(3X))A,
\quad
F=fD'-(f'+2f/X)D.
\]

The inhomogeneous `(C,F)` forcing always has zero cokernel class.  After
introducing two first-stage and three second-stage parameters, the `G'` map
leaves seven nonzero cokernel equations.  The last Wronskian equation
`-BG'+FC'=0` contributes eighteen nonzero coefficient equations.  Altogether
there are only `25` sparse equations, of total degrees three and four, in five
parameters.  Exact matrices, cokernel-functional bases, analytic kernels, and
terminal equations are hash-pinned by the artifact.

## 4. Terminal open and conclusion

The degree condition `deg(B)=8` is a principal open.  In the RREF parameter
basis the leading coefficient is a unit of `K` times `p_1`, so introduce `z`
and impose

\[
zB_8-1=0.
\]

Singular computes, over `K`,

\[
\bigl(\text{seven }G'\text{ cokernel equations},\
      \text{eighteen }J_0\text{ coefficients},\ zB_8-1\bigr)=(1).
\]

This is a Gröbner calculation only in the terminal five deformation
parameters.  It proves that none of the five Belyi top maps extends through
the four lower Wronskian equations with the required Newton degree.  Therefore
the no-vertical-edge Laurent branch is empty.

## 5. Reproduction

Generate the exact certificate:

```bash
.venv/bin/python plane-jc/cas/jc2_degree108_belyi_deformations.py
```

Replay the pinned certificate, including the terminal Singular computation:

```bash
.venv/bin/python scripts/verify_jc2_degree108_belyi_deformations.py
```

For a quick replay of the dessins, coefficient reconstruction, arithmetic
audit, Belyi factorization, and all linear algebra while trusting the pinned
terminal unit-ideal result:

```bash
.venv/bin/python scripts/verify_jc2_degree108_belyi_deformations.py --quick
```

The generated result is
`artifacts/generated-results/jc2_degree108_belyi_deformations.json`.
The software assumptions are `.python-version`, `requirements.txt`, SymPy for
the `S_5` audit, and Singular 4 for the terminal characteristic-zero standard
basis.
