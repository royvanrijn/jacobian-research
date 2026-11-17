# Nonlinear polynomial reparametrizations of the first image have no global conjugate

Let \(K\) have characteristic zero, let \(L=K((h))\), and use the
finite global PBW algebra
\[
 A=A_2(L),\qquad [P,X]=[Z,Q]=h.
\]
Put
\[
 \begin{split}
 r&=X(2-3XQ),\qquad w=1-3XQ,\\
 B&=3X^2P+2wZ,\\
 D_0&=(1+3XQ)P/2-3Q^2Z.
 \end{split}                                                   \tag{1}
\]
Here \(B\) is the unshifted vector-field operator, without the
position term \(-9Q^2\).

**Theorem.** For \(f\in L[t]\), there is a finite global
\(D\in A\) satisfying
\[
                              [D,f(r)]=h                     \tag{2}
\]
if and only if \(f(t)=\alpha t+\beta\), with
\(\alpha\in L^\times\). In this case all solutions are
\[
                D=\alpha^{-1}D_0+C,\qquad C\in\operatorname{Cent}_A(r).
                                                               \tag{3}
\]
The partner has no differential-order or coefficient-shape bound.
Parameter poles in its scalar coefficients are allowed. The exclusion
is global: rational position coefficients do supply a conjugate and
are not being excluded.

## The global polynomial derivation frame

The ordered coefficient rows of \(D_0,B\) in \(P,Z\) have determinant
\[
 \det\begin{pmatrix}(1+3XQ)/2&-3Q^2\\3X^2&2w\end{pmatrix}=1.
\]
Their exact polynomial inverses are
\[
 P=2wD_0+3Q^2B,\qquad
 Z=-3X^2D_0+(1+3XQ)B/2.                                    \tag{4}
\]
No rational position denominator or infinite series is used. Direct
first-order vector-field calculations give
\[
 [D_0,r]=h,\qquad [B,r]=0,\qquad
 [D_0,B]=\frac92hQB.                                        \tag{5}
\]
On positions write
\[
 \delta_0=\frac{1+3XQ}{2}\partial_X-3Q^2\partial_Q,
 \qquad
 \Delta=3X^2\partial_X+2w\partial_Q.
\]
Then \(\delta_0r=1\), \(\Delta r=0\), and
\([\delta_0,\Delta]=(9/2)Q\Delta\); these identities prove (5)
as actual differential operators.

Every global finite operator has a unique finite expansion
\[
                       \sum_{i,j}a_{ij}(X,Q)B^iD_0^j,
                                  \qquad a_{ij}\in L[X,Q].      \tag{6}
\]
Indeed, in differential-order associated graded, the symbols of
\(D_0,B\) are an invertible polynomial linear change of the two
momentum symbols, by the determinant-one calculation. Their monomials
are a basis in every homogeneous momentum degree. Subtracting a lift
of the top degree proves spanning by finite descending induction;
the same symbol independence proves uniqueness.

Consequently the actual global centralizer and Ore extension are
\[
 \mathcal C=\operatorname{Cent}_A(r)
       =L[X,Q][B;h\Delta],\qquad
                  A=\mathcal C[D_0;\operatorname{ad}D_0].       \tag{7}
\]
For the first equality, commuting (6) with \(r\) differentiates its
\(D_0\)-powers and gives \(jhD_0^{j-1}\). Uniqueness and
characteristic zero force all \(j>0\) coefficients to vanish.
For the second, (5) and the polynomial action of \(\delta_0\) on
positions show that \(\operatorname{ad}D_0\) preserves
\(\mathcal C\). Formula (6) supplies both spanning and absence of
additional relations. This is the global normal form already
recorded in [QEP1](EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md),
with the proof recalled to fix its coefficient ring.

The element \(r\), and hence every \(f(r)\), is central in
\(\mathcal C\). Also \(r\) is transcendental over \(L\): it is a
nonconstant polynomial, and \(\deg g(r)=3\deg g\) for every
nonzero nonconstant \(g\in L[t]\).

## The conjugate cannot hide in a higher transverse order

Suppose \(f\) is nonconstant and write a hypothetical partner uniquely
as
\[
                       D=\sum_{j=0}^M a_jD_0^j,
                                     \qquad a_j\in\mathcal C.
\]
Because \([D_0,f(r)]=hf'(r)\), the highest term of the commutator
for \(M>0\) is
\[
 [D,f(r)]
       =Mh\,a_Mf'(r)D_0^{M-1}
                           +\text{terms of smaller }D_0\text{-order}.
                                                               \tag{8}
\]
The coefficient is nonzero: \(A\) is a domain, \(a_M\ne0\), and
\(f'(r)\ne0\) in characteristic zero. Thus \(M\geq2\) is
incompatible with (2). Order zero commutes with \(f(r)\), so
\(M=1\). Equation (2) becomes
\[
                              a_1f'(r)=1.                     \tag{9}
\]
The coefficients \(a_j\) have been allowed every finite \(B\)-order;
that freedom does not make (9) possible for a nonlinear \(f\).
In \(\mathcal C=L[X,Q][B;h\Delta]\), multiplication by the nonzero
central position polynomial \(f'(r)\) preserves \(B\)-order.
Hence (9) forces \(a_1\in L[X,Q]\). The only units in that
polynomial ring are \(L^\times\), so \(f'(r)\) must be constant.
Transcendence of \(r\) then forces \(f\) affine. Conversely (5)
gives (3), and (7) shows that it lists all solutions. A constant
\(f\) has zero commutator with every operator. This proves the
theorem in all finite partner orders.

## A formal inverse does not give a global PBW inverse

For example
\[
                         R_h=r+h r^4
\]
is an actual finite polynomial first image with residue \(r\), but
it has no finite global conjugate. Its derivative
\(1+4hr^3\) is a unit in the h-completion, where
\[
 D_h=(1+4hr^3)^{-1}D_0
       =\sum_{n\geq0}(-4h)^nr^{3n}D_0,
                         \qquad [D_h,R_h]=h.                  \tag{10}
\]
The spatial degrees of these coefficients grow. After rational
position localization (10) is a finite first-order operator. More
generally \(f'(r)^{-1}D_0\) is a finite rational-position conjugate
for every nonconstant \(f\). This positive control is why global
admission, rather than ordinary finite differential order alone,
is essential to the theorem.

The cap-23 odd first-jet direction that changes the first image by a
multiple of \(r^4\) is not removed by this result. What is excluded
is its **exact continuation entirely within \(L[r]\)**, with a
finite global conjugate. A bounded continuation can have higher
position or momentum corrections not lying in \(L[r]\); those
are outside this proof. No evenization, grading normalization,
all-odd-jet obstruction or unrestricted deformed-first-image
classification is asserted.

[Exact frame controls](check_polynomial_first_image.py) and their
[retained receipt](polynomial_first_image_controls.json) check the
polynomial frame, inverse, vector-field commutators, and affine and
rational-localized controls. The global normal-form and arbitrary
partner-order arguments are written proofs. No DC1/DC2 result,
formal verification or external review is claimed.
