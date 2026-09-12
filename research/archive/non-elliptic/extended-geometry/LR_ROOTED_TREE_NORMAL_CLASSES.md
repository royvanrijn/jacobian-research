# Rooted-tree normal classes for the formal LR target lift

This note connects the formal target lift in
[complexity-filtered contact](COMPLEXITY_FILTERED_CONTACT.md) with the
rooted-tree model of a pre-Lie algebra.  For the degree-five map \(F_2\), it
also gives an exact all-order family of nonzero rooted-tree residues in the
third summand of the saturated LR normal module.

The result has a deliberate boundary.  It proves nonvanishing of one
individual tree class in every arity.  The subsequent
[balanced mixed BCH calculation](LR_MIXED_BCH_NORMAL_CLASSES.md) proves that
the tree sum does not cancel in one linear-in-\(X\), opposite-weight sector.
That sector still vanishes on a valid lower-target-jet locus.  Thus neither
result alone is the universal filtered obstruction required by `OP-CCDM`.

The generated certificate is
[`lr_rooted_tree_normal_classes.json`](../artifacts/generated-results/lr_rooted_tree_normal_classes.json).

The standard combinatorial input is the Chapoton--Livernet description of
the free pre-Lie algebra by labelled rooted trees
([arXiv:math/0002069](https://arxiv.org/abs/math/0002069)).  Tree expansions
of Magnus and continuous BCH series are developed, in a closely related
planar/dendriform model, by Ebrahimi-Fard and Manchon
([arXiv:1203.2878](https://arxiv.org/abs/1203.2878)).  The new content here is
the specialization of that dictionary to the nonlinear LR target lift and
its exact saturated normal residues.

## 1. Elementary differentials and the target-lift defect

On polynomial vector fields use the map-composition pre-Lie product

\[
 P\mathbin{\triangleright}Q=(DP)Q.
\]

Let

\[
 \ell_F(Y)=(DF)^{-1}(Y\circ F)
\]

be the differential of the formal target lift.  For a rooted tree \(\tau\)
whose vertex \(v\) is decorated by a target field \(Y_v\), let
\(E_Y(\tau)\) be its target elementary differential and let
\(E_X(\tau)\) be the source elementary differential obtained after replacing
every \(Y_v\) by \(\ell_F(Y_v)\).  Recursively, if the root has decoration
\(Y\) and children \(\tau_1,\ldots,\tau_k\), then

\[
\begin{aligned}
 E_Y(\tau)
   &=D^kY[E_Y(\tau_1),\ldots,E_Y(\tau_k)],\\
 E_X(\tau)
   &=D^k\ell_F(Y)[E_X(\tau_1),\ldots,E_X(\tau_k)].
\end{aligned}
\]

Define the rooted-tree target-lift defect

\[
 \boxed{\Delta_F(\tau)=E_X(\tau)-\ell_F(E_Y(\tau)).}       \tag{1.1}
\]

For the one-edge tree with root \(Y_1\) and child \(Y_2\), differentiation of
\(DF\,\ell_F(Y_1)=Y_1\circ F\) gives

\[
\begin{aligned}
\Delta_F(Y_1(Y_2))
={}&(D\ell_F(Y_1))\ell_F(Y_2)
       -\ell_F((DY_1)Y_2)\\
={}&-(DF)^{-1}D^2F[\ell_F(Y_1),\ell_F(Y_2)].
                                                               \tag{1.2}
\end{aligned}
\]

The right side is symmetric.  Therefore (1.2) is exactly the second
fundamental form \(\operatorname{II}_F(Y_1,Y_2)\) of Proposition 6.4 in the
complexity-filtered note, even before symmetrizing the directed pre-Lie
defect.

Formula (1.1) is the higher rooted-tree continuation of that tensor.  It is
canonical once the ordinary affine connections on source and target are
fixed.  The compiler implements (1.1) for trees decorated by the three
constant target directions.

## 2. The invariant-coordinate compiler

For \(F_2\), write

\[
 F_2=(x^{-2}a(u,\gamma),x^{-1}b(u,\gamma),x\gamma),
\qquad
 \operatorname{wt}(A,B,C)=(-2,-1,1).
\]

The constant fields

\[
 A:=\partial_A,\qquad B:=\partial_B,\qquad C:=\partial_C
\]

have field weights \(2,1,-1\), respectively.  A homogeneous source field of
weight \(p\) is stored as

\[
 x^p(q_r\partial_r+q_u\partial_u+q_\gamma\partial_\gamma),
\qquad r=\log x,
                                                               \tag{2.1}
\]

with \(q_r,q_u,q_\gamma\in\mathbf Q[u,\gamma]\).  The logarithmic
differential matrix \(J\) from the
[torus module](TORUS_FILTERED_LR_MODULE.md) gives

\[
 \ell_{F_2}(A)=x^2J^{-1}e_A,\qquad
 \ell_{F_2}(B)=xJ^{-1}e_B,\qquad
 \ell_{F_2}(C)=x^{-1}J^{-1}e_C                         \tag{2.2}
\]

in the coordinates (2.1).

The calculation must still use the original affine connection in \(x,y,z\);
pre-Lie elementary differentials are not invariant under nonlinear coordinate
changes.  The compiler retains the forced powers of \(x\) and uses

\[
\begin{aligned}
\partial_x(x^kf)
 &=x^{k-1}\left(kf+(u-1)f_u+
   \left(2\gamma-2+\frac87(u-1)\right)f_\gamma\right),\\
\partial_y(x^kf)
 &=x^{k+1}\left(f_u-\frac87f_\gamma\right),\\
\partial_z(x^kf)&=x^{k+2}f_\gamma.                    \tag{2.3}
\end{aligned}
\]

Thus arbitrary branching uses exact higher affine derivatives while all
coefficients remain bivariate polynomials.

For a weight-zero defect, transport by \(J\) and reduce in

\[
 N_R=
 R/(a,b^2)\oplus R/(b,a\gamma)\oplus R/(\gamma),
\qquad R=\mathbf Q[u,\gamma].                         \tag{2.4}
\]

The third row of \(J\) is \((\gamma,0,1)\).  Consequently the third normal
residue is simply

\[
 \rho_\tau(u)=q_\gamma(\Delta_{F_2}(\tau))\big|_{\gamma=0}
 \in\mathbf Q[u].                                    \tag{2.5}
\]

The checker first compiles the one-edge tree \(B(C)\) and recovers

\[
 \rho_{B(C)}
 =-\frac{30}{7}
 (4896u^5-25092u^4+15232u^3-1887u^2+126u-21),        \tag{2.6}
\]

with associated-graded symbol \(-146880u^5/7\).  This independently connects
the tree convention to the already verified quadratic tensor.

## 3. An all-order ladder family

For a label \(R\), write \(R(\tau)\) for grafting a new unary root decorated
by \(R\) above \(\tau\).  Define

\[
 \tau_2=B(C),\qquad \tau_3=A(C(C)),\qquad
 \boxed{\tau_{n+2}=B(C(\tau_n)).}                    \tag{3.1}
\]

Both seeds have total weight zero, and the added pair has weight
\(1-1=0\).  Hence \(\tau_n\) is a weight-zero rooted ladder for every
\(n\ge2\).

Because all decorations are constant target fields,
\(E_Y(\tau_n)=0\) for \(n\ge2\).  Let \(M_A,M_B,M_C\) be the exact
\(3\times3\) polynomial matrices that graft the corresponding source lift
above a homogeneous child.  If \(q_n\) is the logarithmic coefficient vector
of \(\Delta_{F_2}(\tau_n)\), then

\[
 q_{n+2}=Tq_n,\qquad T=M_BM_C.                       \tag{3.2}
\]

Specialize the clean third normal summand at
\((u,\gamma)=(1/6,0)\).  The transfer matrix becomes

\[
 N=
 \begin{pmatrix}
 55165/378&17510/63&8995/144\\
 31195/567&1445/378&1445/144\\
 0&170/7&605/756
 \end{pmatrix}.                                      \tag{3.3}
\]

Its characteristic polynomial is

\[
 \chi_N(\lambda)=
 \lambda^3-\frac{113825}{756}\lambda^2
 -\frac{157250425}{10584}\lambda
 -\frac{15477069875}{428652}.                        \tag{3.4}
\]

Put

\[
 s_k=\rho_{\tau_{2+2k}}(1/6),\qquad
 t_k=\rho_{\tau_{3+2k}}(1/6).
\]

Cayley--Hamilton gives both sequences the recurrence

\[
 r_{k+3}=
 \frac{113825}{756}r_{k+2}
 \frac{157250425}{10584}r_{k+1}
 \frac{15477069875}{428652}r_k.                      \tag{3.5}
\]

All three coefficients are positive.  The exact initial signs are

\[
\begin{aligned}
s_0&=\frac{170}{63}>0,\\
s_1&=-\frac{2022575}{11907}<0,\quad
s_2=-\frac{256808972125}{9001692}<0,\quad
s_3=-\frac{45742985609279375}{6805279152}<0,
\end{aligned}
\]

and

\[
\begin{aligned}
t_0&=-\frac{516875}{9072}<0,\quad
t_1=-\frac{65504665}{979776}<0,\\
t_2&=-\frac{455884893661325}{5184974592}<0.
\end{aligned}
\]

Equation (3.5) preserves strict negativity.  Therefore every \(s_k\) and
every \(t_k\) is nonzero.

### Proposition 3.1

For every \(n\ge2\), the rooted-tree defect
\(\Delta_{F_2}(\tau_n)\) has nonzero image in the third saturated normal
summand \(R/(\gamma)\).  In particular, its weighted associated-graded
normal symbol is nonzero.

### Proof

The preceding recurrence proves
\(\rho_{\tau_n}(1/6)\ne0\), so
\(\rho_{\tau_n}\) is a nonzero polynomial in \(\mathbf Q[u]\).
The filtration \(\deg u=2\) is separated, and the ideal \((\gamma)\) is
homogeneous.  Hence the highest weighted part of \(\rho_{\tau_n}\) is
nonzero in
\(\operatorname{gr}(R/(\gamma))\).  QED

The certificate records explicit symbols through order twelve as regression
data.  Those finite expansions illustrate the family; the
Cayley--Hamilton argument, not the cutoff at twelve, proves the all-order
statement.

## 4. What this says about BCH, and what it does not say

The free pre-Lie algebra on decorated rooted trees maps to polynomial vector
fields by elementary differentials.  Antisymmetrizing the pre-Lie product
gives the map-composition Lie bracket, so every BCH coefficient is a rational
linear combination of the same decorated trees.  This makes (1.1) a concrete
compiler from BCH combinatorics to the filtered LR normal module.

However, Proposition 3.1 alone is not an all-order LR obstruction.

1. The trees \(\tau_n\) contain only target decorations.  Since
   \(\ell_F\) is a Lie homomorphism, target-only **Lie words** lie in the
   lifted target Lie algebra.  Their individual pre-Lie trees can have
   nonzero coordinate-normal defects while cancelling in the
   antisymmetrized Lie combination.
2. A BCH coefficient is a sum of trees.  Nonvanishing of one summand does not
   exclude cancellation with other trees having the same associated-graded
   face.
3. The universal LR obstruction must hold over the lower-jet solution scheme
   in (6.13)--(6.14) of the complexity-filtered note.  The present calculation
   fixes the constant target directions and does not quotient higher tree
   classes by nonlinear combinations generated at lower orders.

The next decisive experiment was therefore narrower than another blind tree
search:

1. add one distinguished weight-zero deformation label \(X\);
2. compile the actual mixed BCH coefficients, not their individual trees;
3. apply the third-summand functional before expanding;
4. prove that one leading tree or one transfer-eigenspace has a nonzero total
   BCH coefficient after all same-face cancellations; and
5. check that this functional descends over the allowed lower target jets.

The first three steps are now carried out in
[the mixed BCH note](LR_MIXED_BCH_NORMAL_CLASSES.md).  The balanced
linear-in-\(X\) coefficient is nonzero in every odd order.  The fourth step
fails for that sector alone because its coefficient is proportional to
\(s^kt^k\) in the two optional opposite-weight target amplitudes.  The
remaining task is therefore a stratified lower-jet dichotomy, not BCH tree
cancellation.

## 5. Reproduction

Generate the exact certificate and the order-two through order-twelve
regression table with

```bash
.venv/bin/python scripts/compile_lr_rooted_tree_classes.py --max-order 12
```

Replay the all-order matrix identity and sign induction using only the Python
standard library with

```bash
python3 scripts/audit_lr_rooted_tree_normal_classes.py
```

The first command uses the pinned SymPy environment.  The second reads the
generated rational matrix and seeds, verifies the \(3\times3\)
Cayley--Hamilton identity with `Fraction` arithmetic, and checks the signs
needed for induction.
