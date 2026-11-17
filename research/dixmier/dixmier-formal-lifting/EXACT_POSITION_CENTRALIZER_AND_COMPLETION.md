# Exact positional Long images and automatic polynomial centralizers

**QEP1.** Keep Long's first image \(R=X(2-3XQ)\) exactly positional
in the finite PBW algebra over \(K[[h]]\), with
\([P,X]=[Z,Q]=h\) and \(K\) of characteristic zero.
Let an actual finite \(S\), with Long's second classical image,
commute with \(R\). If an actual generic finite \(T\) satisfies
\[
                    [R,T]=0,\qquad [S,T]=h,                 \tag{1}
\]
then, writing \(L=K((h))\),
\[
 \operatorname{Cent}_{A_2(L)}(R,S)=L[R,S],\qquad
 \operatorname{Cent}_{A_{2,K[[h]]}(h)}(R,S)=K[[h]][R,S].     \tag{2}
\]
The separate centralizer-finiteness gate in QFC1 is therefore automatic
in this exact-\(R\) category as soon as the \(S\)-mate exists.

If an actual generic individual conjugate of \(R\), commuting with
\(S\), also exists, both conjugates can be made integral, matched to
Long's classical tuple, and flattened by a finite polynomial curvature
correction. This is a conditional bounded-lift theorem. No such actual
individual conjugates or bounded normalization fixing \(R\) are supplied.

## An equality of finite localized quantum algebras

Use [QPH1's exact quantum chart, equations (1)--(2)](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md).
With \(r=R,v=X^{-1}\) and its operators \(a,d,F,G\),
\[
 A_2(L)[X^{-1}]
 =L[r,v,v^{-1}][a;h v^{-1}\partial_v][d;h\partial_r],
 \quad [d,r]=h,\quad [a,v]=h/v.                            \tag{3}
\]
All other chart-generator brackets vanish. In particular this is an
equality of finite Ore algebras. Explicit inverse expressions are
\[
 \begin{split}
 X&=v^{-1},& Q&=(2v-rv^2)/3,\\
 Z&=-3v^{-2}(d-G),&
 P&=-v^3(a-F)+2(rv-1)(d-G).
 \end{split}                                             \tag{4}
\]
The indicated ordering is exact; both chart directions use Laurent
powers of \(v\), polynomial dependence on \(r,a,d\), and no \(r\)
denominators.

In the \(d\)-PBW order commutation with \(r\) is \(h\partial_d\).
Thus \(S,T\), and every element commuting with \(R\), belong to
\[
                \mathcal O=L[r,v,v^{-1}][a;h v^{-1}\partial_v].
                                                               \tag{5}
\]
The classical reduction of \(S\) is
\[
                         s=4a-3ra^2+2v.                  \tag{6}
\]
Its generic \(a\)-order is positive, including when \(h\)-small
higher-order coefficients increase that order.

## A global finite Ore normal form for all admitted fixed-R images

There is also an exact normal form before localizing. Put
\[
 B=3X^2P+2(1-3XQ)Z,\qquad
 D_0=(1+3XQ)P/2-3Q^2Z .                                  \tag{15}
\]
Their coefficient matrix in \(P,Z\) has determinant \(-1\), and
\[
 [B,R]=0,\quad [D_0,R]=h,\quad
 [D_0,B]=\tfrac92hQB,
\]
with exact inverses
\[
 P=2(1-3XQ)D_0+3Q^2B,\qquad
 Z=-3X^2D_0+(1+3XQ)B/2.                                  \tag{16}
\]
These are polynomial vector fields, so these first-order commutators
are verified by their ordinary derivation brackets. The same frame
appears classically in QCR1; the assertions here concern the actual
finite quantum algebra.

Every finite operator has a unique finite expansion
\[
                    \sum_{i,j}f_{ij}(X,Q)B^iD_0^j.         \tag{17}
\]
Indeed differential-order associated graded is the polynomial algebra
in the two momentum symbols. The invertible coefficient matrix in
(15) makes the symbols of \(B,D_0\) another polynomial basis of that
algebra. Lifting and subtracting its top homogeneous piece proves
spanning by descending order; the same symbol independence proves
uniqueness. This works over \(K[[h]]\) as well as \(K((h))\).
Commutation with R in (17) is precisely
\[
 h^{-1}[R,-]=-\partial_{D_0}.
\]
Consequently the complete actual global centralizer of R is
\[
 \boxed{\quad
 \operatorname{Cent}(R)=K[[h]][X,Q][B;h\Delta],\quad
 \Delta=3X^2\partial_X+2(1-3XQ)\partial_Q .\quad}           \tag{18}
\]
The generic version replaces \(K[[h]]\) by \(K((h))\). Formula (18)
classifies all globally admitted finite S and T commuting with R;
it is not a condition checked only in a rational coefficient chart.
It also gives a finite right inverse for the normalized commutator
with R, by integrating (17) in \(D_0\) with a minus sign. It supplies
no right inverse for commutation with an unknown S.

The exceptional R-fiber clarifies what this normal form does not
exclude. Since \(R=X(2-3XQ)\), its two factors generate comaximal
ideals. Both ideals are preserved by \(\Delta\). Thus the quotient
of (18) by R is exactly the product
\[
 A_1(K[[h]];h)\ \times\
 K[[h]][v,v^{-1}][b;h\partial_v].                         \tag{19}
\]
On the first component \(X=0\), use position Q and momentum B/2.
On the second, \(Q=2/(3X)\), use \(v=X^{-1}\), \(b=-B/3\).
Over \(K((h))(r)\), the generic fiber R=r is just the second Laurent
Weyl algebra: r nonzero makes X invertible. Both component algebras
in (19) admit Weyl pairs, so their existence alone is no obstruction
to the required global individual mates. The remaining search is for
finite elements of the single algebra (18) with the prescribed Long
reduction and the required spectral translation relations.

## A finite-order centralizer lemma

For a characteristic-zero constant field \(C\) and \(h\in C^\times\),
put \(\mathcal D=C(v)[a;h\delta]\), \(\delta=v^{-1}\partial_v\).
The centralizer \(Z(U)\) of every positive-order \(U\in\mathcal D\)
is finite as a \(C[U]\)-module.

Indeed let \(m=\operatorname{ord}_aU>0\). In every residue class modulo
\(m\) attained by a nonzero member of \(Z(U)\), choose a member \(B_i\)
of least order \(n_i\). The class zero has \(B_0=1\).
For \(V\in Z(U)\) of order \(n=n_i+km\), \(k\ge0\), leading
coefficient comparison in \([U,V]=0\) gives
\[
           m u_m\delta(v_n)-n v_n\delta(u_m)=0.             \tag{7}
\]
The leading coefficient of \(U^k B_i\) satisfies the same equation.
The ratio of these leading coefficients has derivative zero, so
belongs to \(C\). Subtracting that multiple of \(U^k B_i\) lowers
the order. Repetition terminates. The at most \(m\) chosen \(B_i\)
therefore span \(Z(U)\) over \(C[U]\).
No commutativity theorem, mate-order bound, monicity, or leading
coefficient regularity at \(h=0\) was assumed in this argument.

## The rational-coefficient classical centralizer

Let \(k=K(r)\), \(g=4a-3ra^2\). The substitution \(v=(s-g)/2\)
identifies \(k(v,a)=k(s,a)\), where
\[
                         \{s,-\}=-\frac4{s-g}\partial_a
                       \quad\hbox{with \(s\) held fixed}.
\]
Its fraction-field kernel is \(k(s)\), and
\[
                \operatorname{Cent}_{\{\,,\,\},k(v)[a]}(s)=k[s].
                                                               \tag{8}
\]
For the intersection assertion, \(s\) is a nonconstant polynomial
in \(a\). Coprime \(p(t),q(t)\in k[t]\) stay coprime after evaluation
at \(s\), by their Bezout identity over \(k\). A nonconstant \(q(s)\)
has positive \(a\)-degree and cannot divide \(p(s)\) in \(k(v)[a]\).
Thus a rational function of \(s\) polynomial in \(a\) is polynomial
in \(s\).

This is the centralizer in the larger rational coefficient algebra,
not merely the polynomial global centralizer.

## Local polynomiality with the \(h\)-adic argument made explicit

QET1's finite parameter-pole descent makes the mate \(T\) integral:
it subtracts finite polynomials in \(R,S\), preserving (1), using
QSC1's classical joint-centralizer identity \(K[R,S]\).
See [QET1, Remove parameter poles](ETALE_CENTRALIZERS_AND_GALOIS_COMPLETION.md).

Embed \(L(r)\) faithfully into \(C=K(r)((h))\) and work in
\(\mathcal D=C(v)[a;h\delta]\).
Let \(\mathcal E\) be the finite operators with nonnegative
coefficientwise \(h\)-Gauss valuation. It is a \(K(r)[[h]]\)-algebra,
is \(h\)-torsion-free and \(h\)-separated, and
\[
 \mathcal E/h\mathcal E=K(r)(v)[a],\qquad
                       \mathcal E[h^{-1}]=\mathcal D.      \tag{9}
\]
Normalizing the numerator and denominator of each rational
coefficient in \(v\) proves these assertions. Infinite \(a\)-order
or a spatial completion is not admitted.

The classical centralizer (8) gives QFC1's reduction hypothesis.
The order-semigroup lemma supplies generic module finiteness.
The [finite centralizer theorem QFC1](FINITE_CENTRALIZERS_AND_COMPLETION.md)
now yields
\[
                              Z_{\mathcal D}(S)=C[S].      \tag{10}
\]
For clarity, the commutativity and constant-field conclusions used
here do not require an imported scalar centralizer theorem.
The unique successive polynomial expansion supplied by (8)
embeds the integral centralizer in the commutative completed
polynomial algebra. The mate acts there coefficientwise as
\(-h^{-1}\operatorname{ad}(T)=\partial_s\). Its generic constants
are exactly \(C\). QFC1's Fitting argument makes the finite generic
centralizer etale over \(C[S]\).

In this one-variable application the final covering input can be
limited to the classical curve theorem: over algebraically closed
characteristic-zero constants, a finite connected etale cover of
the affine line has degree one. Its smooth projective completion,
of degree \(e\), genus \(g\), and with \(b\) points above infinity,
satisfies \(2g-2=-e-b\) by Riemann--Hurwitz. Hence \(e=b=1\).
Over general constants a finite etale cover of the affine line
is a constant extension; the horizontal-constant calculation
excludes that extension. The covering step remains an explicit
external geometric input; see QFC1's reference, and
[Riemann--Hurwitz, Stacks 0C1B](https://stacks.math.columbia.edu/tag/0C1B).
No independent verification of those external theorems is claimed.

Descend (10) to \(L(r)\). For an element with coefficients in
\(L(r)(v)\), write the finite expression \(\sum c_jS^j\), \(c_j\in C\).
The distinct powers have distinct \(a\)-orders. Its highest
coefficient ratio puts the highest \(c_j\) in
\(C\cap L(r)(v)=L(r)\), because \(v\) is transcendental over \(C\)
and a rational function over \(L(r)\) constant in \(v\) lies in \(L(r)\).
Subtract and repeat. Therefore
\[
       Z_{L(r)(v)[a;h\delta]}(S)=L(r)[S].                  \tag{11}
\]

## Specialization at every \(r\)-divisor removes denominators

A global joint-centralizer element \(b\) belongs both to (5) and
to \(L(r)[S]\). Write \(b=p(r,S)/q(r)\), with
\(p\in L[r,t]\), \(q\ne0\), and no common irreducible factor
between \(q\) and all numerator coefficients.

If a nonconstant irreducible \(\pi(r)\) divided \(q\), reduction
of \(q(r)b=p(r,S)\) in (5) modulo \(\pi\) would give
\(p(\bar r,\bar S)=0\) over \(L_\pi=L[r]/(\pi)\).
But \(T\) specializes in the same algebra and
\[
                            [\bar S,\bar T]=h .
\]
Thus \(\bar S\) is transcendental over \(L_\pi\).
For positive \(a\)-order this follows from order additivity.
At order zero, \(\bar S\in L_\pi[v,v^{-1}]\); an element there
algebraic over \(L_\pi\) is scalar and has zero commutator.
Both cases exclude algebraicity.
Every coefficient of \(p\) is therefore divisible by \(\pi\),
a contradiction. This proves \(q\) constant and (2) generically.

This argument permits generic operator order to drop at every
special divisor. Preservation of the mate relation supplies the
needed transcendence.

For an integral \(b=\sum c_{ij}R^iS^j\), \(c_{ij}\in L\), a negative
minimum coefficient valuation would, after clearing that pole
and reducing modulo \(h\), give a nonzero polynomial relation
between the independent classical \(R,S\). Thus all \(c_{ij}\)
are integral, proving the second equality in (2).

## An order-drop control for divisor specialization

The specialization argument does not require the leading operator
coefficient to stay nonzero. In the same Laurent chart the exact pair
\[
             S_*=v+ra,\qquad
             T_*=-ra^2/2-va-h/(2v)
\]
satisfies \([S_*,T_*]=h\). At \(r=0\) its first member drops from
\(a\)-order one to zero but is the transcendental coefficient \(v\);
the same specialized mate still has commutator \(h\).
This is a control for the saturation proof, not a Long deformation.

## Finite curvature completion and the remaining construction

Assume an actual generic \(D\) also satisfies
\[
                          [D,R]=h,\qquad [D,S]=0.          \tag{12}
\]
QET1 pole descent makes \(D,T\) integral. By finite polynomials
in \(R,S\), their reductions can be matched to Long's prescribed
classical individual conjugates \(d_0,t_0\): each difference
commutes classically with both positions, so belongs to
\(K[R,S]\).

Jacobi and (2) give \(C=h^{-1}[D,T]\in K[[h]][R,S]\).
The finite polynomial
\[
                         b(R,S)=-\int C(R,S)\,dR
\]
satisfies
\[
                   [D,T+b]=0,\quad [S,T+b]=h,\quad
                   [R,T+b]=0.                            \tag{13}
\]
Indeed \([D,-]/h=\partial_R\) on the centralizer. The prescribed
classical tuple is flat, so \(C\) and \(b\) are \(h\)-divisible;
the correction preserves that entire tuple.

All images are actual finite PBW elements and thus share one
finite spatial bound. For Long's noninvertible classical frame,
QFL2 would then give a nonsurjective complex Weyl endomorphism.
This proves a conditional construction, not existence of the
inputs (1) and (12).

In the exact-\(R\) scalar category the remaining gate consists
of the two actual universal spectral translation classes.
There is no additional centralizer-finiteness or curvature
obligation once these vanish. One individual conjugate does
not here imply the other. A formal, division-ring, or
selected-fiber mate cannot replace an actual finite mate.

Genuinely matrix-valued positions, general projective carriers,
and finite normalization of deformed first images remain outside
the theorem. No DC1/DC2 solution, actual bounded Long lift,
formal verification or external review is claimed.

