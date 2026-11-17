# Polynomial position connections preserve finite-support admission

Let \(k\) have characteristic zero and
\[
 A_h=k[[h]]\langle x_1,\ldots,x_n,p_1,\ldots,p_n\rangle,\qquad
 [p_i,x_j]=h\delta_{ij},\quad [x_i,x_j]=[p_i,p_j]=0.
\]
Here \(A_h\) consists of finite PBW sums with scalar-series coefficients.
Its \(h\)-completion \(\widehat A_h\) additionally permits unbounded
spatial degrees across parameter orders.

**Theorem.** Let \(M_i\in M_r(k[[h]][x_1,\ldots,x_n])\) be polynomial
position matrices with finite common spatial support and
\[
       \partial_iM_j-\partial_jM_i+[M_i,M_j]=0.            \tag{1}
\]
Put \(L=\max(1,\deg M_1,\ldots,\deg M_n)\). The assignment
\[
 \iota(x_i)=x_iI_r,\qquad
 \iota(p_i)=p_iI_r+hM_i,\qquad \iota(h)=hI_r               \tag{2}
\]
is an injective homomorphism \(A_h\to M_r(A_h)\), with a continuous
extension to the completions. For
\(F=\sum_{m\geq0}h^mf_m\in\widehat A_h\), \(f_m\in A_0\) in PBW order,
\[
 \sup_m\deg_{\rm PBW} f_m<\infty
       \quad\Longleftrightarrow\quad
 \sup_m\deg_{\rm PBW}[h^m]\iota(F)<\infty.                 \tag{3}
\]
The degree of a matrix is the maximum of its entry degrees.
More precisely, a bound \(D\) on either side supplies a bound \(LD\)
on the other side. This concerns composition with (2), not arbitrary
matrix actions with the same classical connection data.

## Actual source admission

All images in (2) are actual finite polynomial matrices.
The position relations and mixed relations follow immediately.
For the momentum relation, ordinary matrix multiplication and the
source Weyl rule give
\[
 [\iota(p_i),\iota(p_j)]
 =h^2\bigl(\partial_iM_j-\partial_jM_i+[M_i,M_j]\bigr)=0.
\]
Thus the homomorphism exists over the uncompleted finite-PBW ring.
It reduces modulo \(h\) to scalar amplification of the identity map.

The completion extension is well-defined: all entries of every
\(\iota(f_m)\) are regular in \(h\), so only finitely many terms
contribute at any given parameter order. If \(F\ne0\), its first
nonzero parameter coefficient maps modulo \(h\) to that same
coefficient times \(I_r\), which is nonzero. This proves injectivity,
also on the completion. No rational matrix gauge or scalarization
theorem is a premise.

## A filtration proves both directions of boundedness

Assign weight one to each position \(x_i\), weight \(L\) to each
momentum \(p_i\), and weight zero to \(h\). Write this positive
weighted PBW filtration as \(\mathcal F_D\). Formula (2) preserves it,
since every \(M_i\) has position degree at most \(L\).

For every fixed \(D\) there is the sharper equivalence
\[
 \bigl(f_m\in\mathcal F_D\ \hbox{for all }m\bigr)
 \quad\Longleftrightarrow\quad
 \bigl([h^m]\iota(F)\in M_r(\mathcal F_D)\
                                      \hbox{for all }m\bigr).   \tag{4}
\]
The forward implication follows from filteredness. For the converse,
reduce the image modulo \(h\) to obtain \(f_0I_r\), hence
\(f_0\in\mathcal F_D\). Inductively assume the same for \(f_0,\ldots,f_{m-1}\).
Subtract
\(\sum_{j<m}h^j\iota(f_j)\) from \(\iota(F)\).
The remainder still has every coefficient in \(M_r(\mathcal F_D)\).
Divide by \(h^m\) and reduce modulo \(h\). The result is \(f_mI_r\),
so \(f_m\in\mathcal F_D\). This proves (4).

Every PBW monomial satisfies
\[
             \deg_{\rm PBW}\leq\deg_{\mathcal F}
                                  \leq L\deg_{\rm PBW}.
\]
Using these inequalities in (4) proves both quantitative directions
of (3). This is an induction on parameter order, not an inference
from a finite sample or from a formal inverse to (2).
It also handles any fixed finite parameter pole after multiplying
the series by a central power of \(h\). The connection matrices in
the hypothesis themselves remain parameter-regular.

In particular a scalar formal Long tuple cannot acquire a uniform
spatial bound merely by applying (2). If its compositional image
lies in \(M_r(A_h)\), the original scalar tuple already lies in \(A_h\).
There is no claim that every formal quantization of the same
projective connection is compositional.

## An actual nonconstant global rank-two source datum

Take two positions \(X,Q\), momenta \(P,Z\), and \([P,X]=[Z,Q]=h\).
Set
\[
 t=XQ,\qquad M(t)=
 \begin{pmatrix}0&1\\t&0\end{pmatrix},\qquad
 M_X=QM(t),\quad M_Q=XM(t).                              \tag{5}
\]
The position matrices commute, and
\[
       \partial_X(XM(t))=\partial_Q(QM(t))
                            =M(t)+XQ\,M'(t).
\]
Thus (1) holds exactly. Formula (2) becomes
\[
 \begin{split}
 \iota(X)&=XI_2,\qquad \iota(Q)=QI_2,\\
 \iota(P)&=PI_2+h
                  \begin{pmatrix}0&Q\\XQ^2&0\end{pmatrix},\\
 \iota(Z)&=ZI_2+h
                  \begin{pmatrix}0&X\\X^2Q&0\end{pmatrix}.       \tag{6}
 \end{split}
\]
All six Weyl relations hold as actual finite polynomial matrix
identities over the original global source algebra.
Here \(L=3\), so (3) has factor-three bounds in both directions.

The underlying connection is \(d+M(t)\,dt\). It has no invariant
rational line over the POSITION field \(k(X,Q)\). Indeed, hold \(Q\)
constant and write \(k(X,Q)=k(Q)(t)\). Invariance under
\(\partial_X+QM(t)\) implies invariance under
\(\partial_t+M(t)\). For a line generated by \((f,g)^t\),
\(f=0\) is impossible from the first component. For \(f\ne0\),
\(w=g/f\) satisfies \(w'-w^2=-t\). Thus \(u=-w\) would solve
\[
                               u'+u^2=t.                \tag{7}
\]
There is no rational solution over any characteristic-zero constant
extension. If the degree at infinity of nonzero \(u\) is an integer
\(a\geq1\), its square has degree \(2a\), while its derivative has
degree \(a-1\), so their sum cannot have degree one. If \(a\leq0\),
both terms have degree at most zero. The zero function also fails.

After algebraically closing constants, a flat constant-matrix
position connection has commuting constant matrices and a common
eigenline. A rational position gauge transporting (5) to such a
connection would give the line just excluded. Therefore this datum
is not a rational POSITION gauge of a constant connection.
No irreducibility or gauge-exclusion claim over the full Weyl
division ring is made.

## What changes, and what remains missing for Long

For any scalar formal quantization \(F_i\) of the canonical full
Long graph, the compositional matrices \(\iota(F_i)\) have classical
residue \(f_iI_2\), hence classical trace \(2f_i\), a polynomial in
the original source variables. This differs from the fixed shifted
rank-five model's trace \(20/(3r)\). The obstruction for that fixed
generic rank-five representation does not exclude this different
rank-two graph datum. Formula (6) supplies actual global SOURCE
admission, while the formal target action may still have unbounded
support.

This datum does not itself correct the target relations. If chosen
finite scalar representatives \(G_i\) have nonzero curvature
\(C_{ij}=[G_i,G_j]-hJ_{ij}\), then
\[
 [\iota(G_i),\iota(G_j)]-hJ_{ij}I_2=\iota(C_{ij})\ne0,
\]
by injectivity. Applying (6) to uncorrected Long representatives
therefore preserves their nonzero curvature.
Applying it to a corrected scalar formal lift preserves the
unboundedness issue by (3).

A genuinely finite noncompositional matrix/graph target action,
even with this classical projective datum, remains outside the
theorem. Source-finite projectivity, finite target matrices, any
needed scalar realization, and a bounded Long lift remain separate
obligations. The theorem does not classify all nonscalar carriers
and gives no DC1/DC2 conclusion.

The [exact checker](check_polynomial_connection_admission.py) verifies
the concrete polynomial Weyl relations (6), the flatness identity
and finite filtration controls, and rejects an incorrect-sign
connection. Its [receipt](polynomial_connection_admission_controls.json)
pins the proof and checker. The universal boundedness result is the
written filtration induction above; finite controls do not prove
its universal scope. Formal verification, independent whole-proof
replay and external review are not asserted.
