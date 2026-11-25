# Orthogonal block normalizers: the spherical boundary and a safe coordinate algebra

Work over an algebraically closed field k of characteristic zero. Let n>=2, 1<=p<=n, V=U directsum W with dim(U)=p and dim(W)=n-p, and let H=SO(U,q) act trivially on W. Put G=SL(V), N=N_G(H), the actual full scheme normalizer, and X=G/N. By choosing a q-orthonormal basis of U, take q=I_p. Fix the upper triangular Borel B_+ of SL_n.

**Theorem 1 (exact spherical boundary).** X is spherical if and only if p=1 or p=n.

**Theorem 2 (explicit coordinates and Borel obstruction).** For every proper block 2<=p<n, there is an explicit dominant regular map

\[
 (C,\Delta):X\longrightarrow\mathbb A^2
\]

whose coordinates have B_+ characters chi and chi^2, where chi(b)=b_nn/b_11. The rational function Delta/C^2 is a nonconstant B_+-invariant function. The polynomial subalgebra k[C,Delta] satisfies the exact purity/safety criterion of Theorem 3 below. In particular the rational ratio does not itself supply a regular unsafe phase.

The proof uses classical homogeneous-quotient and same-weight-ratio methods; Section 6 gives their precise roles.

## 1. Full normalizers, including the reducible SO2 standard module

Assume 2<=p<=n. The fixed subspace V^H is W, while the span of all (h-1)V is U. To see the latter assertions without irreducibility, use the SO_p elements which negate two selected coordinate vectors and fix the others. Every U-coordinate is changed by one such element, and 2 is invertible. These descriptions are preserved under normalizers, so every element of N is block diagonal on U,W.

The space of H-invariant quadratic forms on U is exactly kq, also at p=2. Indeed let S be the symmetric matrix of an invariant form. Differentiating the action of J_ij=E_ij-E_ji in so_p gives

\[
 [S,J_{ij}]=0\qquad(i<j).
\]

The diagonal entries of this equation give S_ij=0, and its (i,j) entry gives S_ii=S_jj. Hence S is scalar. When p=2 this uses only J_12 and proves the claim even though the standard SO2 module splits into two characters over k. It is uniqueness in Sym^2(U*), not irreducibility of U, that is needed.

Consequently, writing A_0 for the block on U, every normalizer element satisfies

\[
 A_0^T A_0=\lambda I_p\quad\text{for some}\quad\lambda\in k^*.
\]

Conversely a similitude conjugates SO_p to itself. Thus for r=n-p>0,

\[
\boxed{N=\{\operatorname{diag}(A_0,B_0):
 A_0\in\operatorname{GO}_p,
 B_0\in\operatorname{GL}_r,
 \det(A_0)\det(B_0)=1\}.}
\tag{1}
\]

For p=n, the formula is N=GO_n intersection SL_n.

These are full scheme identities, rather than only identity-component formulas. The argument is functorial over k-algebras: the fixed submodule and the span of (h-1)V are recovered using the same constant diagonal rotations with 2 invertible; the linear equations for an invariant symmetric matrix have the same scalar solutions over every k-algebra. The multiplier is a unit because A_0 is invertible. The converse conjugation check is valid over every such algebra. In particular all reflection components in GO2 are retained, with their determinants compensated in the W block when r>0. No rational coset representative or inner-automorphism theorem is assumed.

For p=1, H=SO1 is the trivial group, so N=G. Formula (1) is not applied to that case.

## 2. Actual regular semi-invariants on the quotient

Suppose 2<=p<n and put P=diag(I_p,0). For g in SL_n define

\[
 A(g)=\sum_{i=1}^{p}g_{ni}^2,\qquad
 B(g)=\sum_{i=1}^{p}(g^{-1})_{i1}^2,\qquad
 C(g)=(gPg^{-1})_{n1},\qquad \Delta(g)=A(g)B(g).
\tag{2}
\]

All are regular functions on SL_n, since matrix inversion is regular there. A and B need not individually descend to X.

Under right multiplication by the full normalizer element in (1), A is multiplied by lambda and B by lambda^{-1}; C is unchanged because that element commutes with P. Thus C and Delta are right-N-invariant and descend as actual regular functions on G/N. The homogeneous quotient is smooth, and right-N-invariant regular functions descend through it [Brion, Theorem 1.16]. The full group in (1), including its reflection components, is used here.

For b in the upper triangular Borel, its last row is b_nn times the last coordinate row, and the first column of b^{-1} is b_11^{-1} times the first coordinate column. Hence

\[
 C(bg)=\chi(b)C(g),\qquad
 \Delta(bg)=\chi(b)^2\Delta(g),\qquad
 \chi(b)=b_{nn}/b_{11}.
\tag{3}
\]

The same equalities with inverse characters hold under the usual inverse-left-translation action on functions. In either convention

\[
 R=\frac{\Delta}{C^2}\in k(X)^{B_+}.
\tag{4}
\]

In the natural variables (g,g^{-1}), C has bidegree(1,1) and Delta has bidegree(2,2). They are not asserted to have uniformly bounded ordinary polynomial degree in the entries of g alone: replacing inverse entries by cofactors gives degrees n and 2n, respectively.

## 3. A two-parameter slice and the nonconstant ratio

Use the actual lower-unipotent matrices

\[
 g_{x,y}=I+xE_{n1}+yE_{n2},\qquad
 g_{x,y}^{-1}=I-xE_{n1}-yE_{n2}.
\]

They are available exactly because p>=2 and n>p. The square of xE_n1+yE_n2 is zero. Direct substitution in (2) gives

\[
 C(g_{x,y})=x,\qquad
 A(g_{x,y})=x^2+y^2,\qquad B(g_{x,y})=1,
 \qquad\Delta(g_{x,y})=x^2+y^2.
\tag{5}
\]

The induced map (x,y)↦(x,x^2+y^2) is dominant (indeed surjective on geometric points). Thus C and Delta are algebraically independent, k[C,Delta] is a polynomial algebra in two variables, and the map in the theorem is dominant.

Taking x=1 gives R(g_{1,y})=1+y^2, so R is nonconstant. If B_+ had an open orbit in X, every B_+-invariant rational function would be constant on that dense orbit and hence constant in k(X). This contradicts (4)--(5). Therefore X is nonspherical for every proper block 2<=p<n, with no dimension bound and no representation-table or Weyl-formula premise.

The ratio R itself is not regular on all X. If it were, its pullback along the regular map (x,y)↦g_{x,y}N would be regular on the affine plane; but by (5) that pullback is (x^2+y^2)/x^2, which has a pole at x=0,y=1. This is why the rational Borel obstruction is not automatically a regular bad phase.

## 4. The two spherical boundary cases

If p=1, N=G and X is a point, hence spherical.

If p=n>=2, the multiplier on N=GO_n intersection SL_n satisfies lambda^n=1, since (det A_0)^2=lambda^n and det A_0=1. Its kernel is SO_n, so

\[
 \dim N=n(n-1)/2,\qquad
 \operatorname{Lie}N=\mathfrak{so}_n.
\]

At the identity coset the Borel stabilizer has Lie algebra b_+ intersection so_n=0: an upper triangular skew-symmetric matrix is zero. Therefore its dimension is zero, and

\[
 \dim(B_+\cdot eN)=\dim B_+
   =n(n+1)/2-1
   =\dim(SL_n/N).
\]

A group orbit is locally closed and has dimension equal to the group dimension minus the stabilizer dimension [Brion, Proposition 1.11]. Thus an orbit of full dimension in the irreducible homogeneous variety is open. Hence X is spherical. This treats the actual full normalizer, not only SO_n. It includes the original n=p=3 exception and p=n=2. Together with Section 3 it proves both directions of the exact criterion.

For any characteristic-zero form of this pair that becomes the displayed pair over an algebraic closure, the same result is a geometric statement after that scalar extension. This does not assert new rational points, rational Borels, rational tensor bases or classification of arbitrary forms.

## 5. An exact obstruction to seeds confined to the two-coordinate algebra

Assume 2<=p<n. Let \mathcal R=k[X], and suppose Lambda:\mathcal R→k is any normalized left-G-invariant linear functional. Put S=k[C,Delta].

**Theorem 3 (exact conditional purity/safety criterion).** For every f in S, the following three conditions are equivalent: every positive pure moment vanishes; f(eN)=0; and for every fixed marker h in \mathcal R the marked moments Lambda(hf^m) vanish for all sufficiently large m. The assertion is conditional on the given functional and does not infer its existence from nonsphericity.

Consider the diagonal one-parameter subgroup

\[
 \gamma(s)=\operatorname{diag}(s^{-1},1,\ldots,1,s).
\]

Equation (3) gives

\[
 C(\gamma(s)g)=s^2 C(g),\qquad
 \Delta(\gamma(s)g)=s^4\Delta(g).
\]

Thus S has strictly positive weights on both generators, with degree-zero subspace k. Every h in \mathcal R has finite integer weight support for this torus: its pullback to SL_n is a regular polynomial in matrix entries and inverse entries, each of which has a definite integer weight.

Left invariance of Lambda kills every nonzero weight. For f in S write c=f(0,0)=f(eN). The weight-zero term of f^m is exactly c^m, so

\[
 \Lambda(f^m)=c^m\qquad(m>=1).
\tag{6}
\]

If c=0, every weight of f^m is at least2m. If the minimum weight of a fixed marker h is j_min, then all weights of hf^m are positive for 2m+j_min>0, and hence Lambda(hf^m)=0 for all such m. If c!=0, the marker h=1 detects every order by (6). We obtain the exact within-subalgebra equivalence

\[
\boxed{
 \forall m>=1:\Lambda(f^m)=0
 \quad\Longleftrightarrow\quad f(eN)=0
 \quad\Longleftrightarrow\quad
 \forall h\in\mathcal R:\Lambda(hf^m)=0\text{ for }m\gg_h0.
}
\tag{7}
\]

Therefore every pure phase obtainable as a polynomial in C and Delta is all-marker safe in the full ambient coordinate ring. This does not prove a global Mathieu property of X, nor exclude phases involving other coordinates, opposite torus weights or Laurent variables. It states exactly why the regular functions producing the nonsphericity ratio do not, by themselves, give the desired unsafe seed.

## 6. Classical sources and scope

[Brion] Michel Brion, *Introduction to actions of algebraic groups*, Les cours du CIRM 1(1) (2010), 1--22, [doi:10.5802/ccirm.1](https://doi.org/10.5802/ccirm.1); [retained primary](sources/brion-actions-of-algebraic-groups-2010.pdf).

Proposition 1.11 and its proof on printed4 supply the usual orbit-dimension and local-closedness facts used in Section 4. Theorem 1.16 and its proof on printed6 give smooth homogeneous quotients and descent of invariant regular functions, as used in Section 2. These standard algebraic arguments apply over algebraically closed characteristic-zero fields; no embedding of the whole coefficient field into C is assumed. The normalizer and coordinate computations in Sections1--5 are explicit in this note.

For comparison, Example2.11.4 on printed18 identifies the nonzero determinant hypersurfaces in symmetric matrices with SL_n/SO_n and gives their open principal-minor Borel cell. Lemma2.12 on that page gives the classical same-weight-ratio characterization of affine spherical varieties. The present proof needs only the elementary direction that a nonconstant invariant rational function rules out a dense open Borel orbit; it does not import a spherical-subgroup classification or a converse Mathieu assertion.

The theorem about S concerns only the stated polynomial coordinate algebra and the given normalized functional. It does not prove a global Mathieu property of X or exclude phases involving other coordinates or opposite weights. The ratio Delta/C^2 has a genuine pole on the displayed two-parameter slice and is not declared a regular phase.
