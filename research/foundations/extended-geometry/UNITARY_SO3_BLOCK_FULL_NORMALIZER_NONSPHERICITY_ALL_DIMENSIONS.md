# Full SO3-block normalizer quotients are nonspherical in every dimension at least four

Let k be algebraically closed of characteristic zero. Let V=U directsum W, with dim(U)=3 and dim(W)=r=n-3>=1. Fix a nondegenerate quadratic form q on U. Embed H=SO(U,q) in G=SL(V), acting trivially on W, and let N=N_G(H) be its actual full normalizer. Write theta for the highest root of SL_n.

**Theorem.** One has

\[
\dim L(2\theta)^N=2.
\]

In particular G/N is nonspherical for every n>=4. Two explicit full-N-fixed vectors are given in Section 4. The geometric conclusion applies to the specified rational unitary forms of the SO3-block construction after scalar extension to an algebraic closure. No descent of these representation vectors to arbitrary rational forms is asserted.

## 1. The actual full normalizer

The three-dimensional standard H-module U is irreducible and nontrivial, whereas H acts trivially on W. Thus W=V^H and U is the span of all (h-1)V. Both are preserved by every normalizer element. An element of N consequently has blocks (A,B) on U,W.

The space of H-invariant bilinear forms on U is one dimensional, spanned by q: the form q identifies U with U*, and Schur's lemma applies to the irreducible U. Therefore A^*q is a nonzero scalar multiple of q, so A belongs to GO(U,q). Over k, choose a scalar z whose square is this multiplier. The isometry z^{-1}A has determinant either1 or-1. Changing z to -z changes that determinant in odd dimension3, so choose z with z^{-1}A in SO(U,q). Hence

\[
N=H\,C,\qquad
 C=\{\operatorname{diag}(z I_U,B):z^3\det B=1\}.
\tag{1}
\]

Conversely each displayed element normalizes H, so (1) is equality with the full normalizer. The factors H,C commute and intersect trivially: a scalar in SO_3 has z^2=z^3=1, hence z=1.

For r>=1, C is a connected reductive group. It is the kernel of the primitive central character (z,B)↦z^3 det(B) of G_m times GL_r; explicitly its variety is isomorphic to G_m times SL_r by writing B=diag(z^{-3},1,...,1)B_0. Projection C→GL(W) is surjective over k, since the required cube root of det(B)^{-1} exists. These facts include r=1.

## 2. Four distinct irreducible summands

The symmetric square decomposes as

\[
\operatorname{Sym}^2V=
 k q^{-1}\ \oplus\ \operatorname{Sym}^2_0U\ \oplus\
 (U\otimes W)\ \oplus\ \operatorname{Sym}^2W.
\tag{2}
\]

The first two summands have C-character z^2. The SO_3 module U has dimension3, and its traceless symmetric square is irreducible of dimension5. One can see the latter on the simply connected SL_2 cover: U=Sym^2(k^2), and the highest vector in its symmetric square generates Sym^4(k^2), of dimension5; the invariant quadratic tensor supplies the remaining dimension1.

The module Sym^2W is irreducible for GL(W), including when r=1, and stays irreducible on C by the surjective projection. The C-action on U tensor W is the irreducible twist z times W, while H acts irreducibly on its U factor. Since the factors commute, U tensor W is irreducible under N.

These four N-modules are pairwise nonisomorphic. Their SO_3 types distinguish all but the two SO_3-trivial summands. The central one-parameter subgroup of C given by

\[
z=s^r,\qquad B=s^{-3}I_W
\]

has weights respectively 2r, r-3 and -6 on Sym^2U, U tensor W and Sym^2W. These three weights are distinct for every r>=1, so the two trivial SO_3 types cannot be confused even at r=1. Thus Schur's lemma gives

\[
\dim\operatorname{End}_N(\operatorname{Sym}^2V)=4,
\qquad \dim\operatorname{End}_N(V)=2.
\tag{3}
\]

## 3. The actual contraction kernel is L(2theta)

Use the natural symmetric pairing

\[
\langle\phi\psi,uv\rangle=
 \tfrac12\{\phi(u)\psi(v)+\phi(v)\psi(u)\}.
\]

Identify End(Sym^2V) with Sym^2V tensor Sym^2V*. Define the G-equivariant partial contraction

\[
\begin{split}
c(uv\otimes\phi\psi)=\tfrac14\{&
 \phi(u)v\otimes\psi+\phi(v)u\otimes\psi\\
 &+\psi(u)v\otimes\phi+\psi(v)u\otimes\phi\}.
\end{split}
\tag{4}
\]

In particular c(u^2 tensor phi^2)=phi(u)u tensor phi. Every rank-one operator with nonzero trace is in the image up to a nonzero scalar. Such operators span End(V): if phi(u)=0, choose psi(u)!=0 and write u tensor phi=u tensor(phi+psi)-u tensor psi. Therefore c is surjective.

For the positive roots of an upper triangular Borel, e_1^2 tensor(e_n*)^2 is a highest-weight vector of weight2theta. It belongs to ker(c). Characteristic-zero complete reducibility and highest-weight theory put a copy of L(2theta) in that kernel.

Let d=n(n+1)/2. Surjectivity gives

\[
\dim\ker(c)=d^2-n^2=\frac{n^2(n-1)(n+3)}4.
\tag{5}
\]

Weyl's dimension formula [E, Proposition26.8] gives exactly this dimension for 2theta. In type A_(n-1) the weight has tuple (2,0,...,0,-2), equivalently the GL_n partition (4,2,...,2,0). The nontrivial factors in the Weyl product are

\[
 \prod_{j=2}^{n-1}\frac{j+1}{j-1}=\frac{n(n-1)}2,
 \quad \frac{n+3}{n-1},\quad
 \prod_{i=2}^{n-1}\frac{n-i+2}{n-i}=\frac{n(n-1)}2.
\]

Their product is n^2(n-1)(n+3)/4. Hence the included irreducible summand exhausts the kernel:

\[
\ker(c)=L(2\theta).
\tag{6}
\]

This identification uses an actual contraction and the displayed dimension calculation; it does not require a tensor-product decomposition table.

## 4. An explicit basis of the full-normalizer invariants

Choose a q-orthonormal basis e_1,e_2,e_3 of U and any basis w_1,...,w_r of W. Put

\[
 Q_U=\sum_i e_i^2,\qquad
 \tau_U=\sum_i(e_i^*)^2,
\]

extending tau_U by zero on the other blocks of Sym^2V. Then tau_U(Q_U)=3. Define actual endomorphisms

\[
 P_0=\tfrac13Q_U\otimes\tau_U,\quad
 P_5=P_{\operatorname{Sym}^2U}-P_0,\quad
 P_{\mathrm{cross}}=P_{U\otimes W},\quad
 P_W=P_{\operatorname{Sym}^2W}.
\tag{7}
\]

Here the block projectors refer to the decomposition of Sym^2V induced by V=U directsum W. The formula for P_0 is independent of the chosen q-orthonormal basis. It projects onto the quadratic tensor using contraction with q. Each P commutes with the full N. In particular the scalar z multiplies Q_U and tau_U inversely, so it fixes P_0. The four projectors form a basis of the invariant endomorphism space in (3).

The pairing convention makes the dual of v_i v_j for i<j equal to 2v_i*v_j*. Formula (4) consequently gives

\[
 c(v_i^2\otimes(v_i^*)^2)=v_i\otimes v_i^*,
\quad
 c(v_i v_j\otimes2v_i^*v_j^*)=
 \tfrac12(v_i\otimes v_i^*+v_j\otimes v_j^*).
\]

Summing the diagonal and off-diagonal basis tensors on a d-dimensional block yields c(Id_{Sym^2(k^d)})=(d+1)Id/2. Summing mixed tensors yields one half of each of its two coordinate projections per pair (i,alpha). In Q_U tensor tau_U the terms with different indices contract to zero. Thus

| Endomorphism | Id_U coefficient in its contraction | Id_W coefficient |
| --- | ---: | ---: |
| P_0 | 1/3 | 0 |
| P_5 | 5/3 | 0 |
| P_cross | r/2 | 3/2 |
| P_W | 0 | (r+1)/2 |

The images of P_0 and P_W are independent. The restriction of contraction to the four-dimensional invariant space has rank2 and kernel of dimension2, with the explicit basis

\[
\boxed{\quad 5P_0-P_5,\qquad
 P_{\mathrm{cross}}-\tfrac{3r}{2}P_0-
                       \tfrac{3}{r+1}P_W.\quad}
\tag{8}
\]

The first has a nonzero P_5 component and no mixed component, whereas the second has mixed component one. They are independent also at r=1. Equations (6)--(8) prove the theorem. These calculations are identities for every r>=1, rather than conclusions from a finite sample.

For comparison, the same dimension follows from the exact sequence of N-modules

\[
0\longrightarrow L(2\theta)\longrightarrow
 \operatorname{End}(\operatorname{Sym}^2V)
 \mathrel{\mathop{\longrightarrow}^{c}}
 \operatorname{End}(V)\longrightarrow0:
\]

reductivity of N in characteristic zero makes taking invariants exact, and (3) gives4-2=2. The explicit calculation above supplies the two vectors directly.

## 5. The geometric nonsphericity consequence

Let v_1,v_2 be independent vectors in L(2theta)^N, and choose a nonzero Borel eigenvector phi in the dual irreducible module. The matrix coefficients

\[
 f_i(gN)=\phi(gv_i)
\]

are regular functions on the homogeneous quotient: on G they are regular and right-N-invariant, so descend through G→G/N. They are linearly independent. Indeed, if a linear combination vanishes on G, the span of the G-translates of phi annihilates the same combination of the v_i; that span is the whole dual irreducible representation.

Both functions have the same left Borel character, with the inverse convention if translation is by b^{-1}. Their ratio is therefore a nonconstant B-invariant rational function on G/N. An open Borel orbit would be dense, and every B-invariant rational function would be constant on it. This contradiction proves nonsphericity.

Only the implication from these explicit invariant vectors to nonsphericity is needed. No classification of spherical subgroups or converse multiplicity-free theorem is used.

## 6. Geometric forms, boundary and classical inputs

For every specified positive q_3 directsum q_r unitary model of the [all-dimensional radius construction](UNITARY_SO3_BLOCK_ALL_DIMENSIONS_EXPLICIT_RADIUS_COUNTEREXAMPLES.md), scalar extension to an algebraic closure identifies the pair and its actual full normalizer with the pair above. Thus those normalizer quotients are geometrically nonspherical for all n>=4, including the range n>=8 beyond the earlier dimension test. This is a statement about geometric quotients, not new descent of the tensors in (8) to arbitrary rational forms.

For n=3 the proof does not apply: W has dimension0 and the four-block argument no longer holds. The quotient by SO_3 times mu_3 is not declared nonspherical by this certificate. No new Mathieu-existence theorem or field descent of arbitrary bad phases is inferred. The radius construction and its moment, trace and nonpower assertions retain their independent proofs.

The imported representation foundations are characteristic-zero complete reducibility, Schur's lemma and highest-weight theory. In the standard matrix model all tensor operators and weight data can be constructed over Q. Their dimensions and the ranks of the displayed linear maps are preserved in characteristic zero; equivalently the usual algebraic highest-weight framework gives the same L(2theta) over each algebraically closed characteristic-zero field. The theorem is therefore not restricted to complex coefficients.

[E] Pavel Etingof, MIT18.745, *Lie Groups and Lie Algebras I*, Fall2020, [Lecture26: The Weyl character formula](https://ocw.mit.edu/courses/18-745-lie-groups-and-lie-algebras-i-fall-2020/mit18_745_f20_lec26.pdf), Proposition26.8 and its complete Section26.5 derivation, printed142--143. The retained [exact PDF](sources/etingof-weyl-character-lecture26-2020.pdf) has SHA256 efbc773b03a439d9398d1605224d1fa4ffdf0a36c9e5acddecddd5f5d85dce02. The lecture's preceding complete-reducibility/highest-weight framework is a classical input; no audit of the entire course is asserted. This source is distinct from the previously retained 2023 Chevalley-restriction lecture.
