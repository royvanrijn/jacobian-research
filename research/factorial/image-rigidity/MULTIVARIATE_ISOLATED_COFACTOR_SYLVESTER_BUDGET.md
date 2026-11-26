# Cofactor Smith data for isolated common factors

Let \(F,G\in\mathbb Z[z]\) be primitive of degrees \(n,N\) with
nonzero endpoint coefficients. Suppose their greatest common divisor
over \(\mathbb Q\) has primitive integer form \(H\) of degree \(d\ge1\).
Write \(F=HA\), \(G=HB\), \(a=n-d\), and \(b=N-d\). Assume \(a,b\ge1\).
The cofactors are primitive and coprime. The application to an
isolated nonreal common quadratic has \(d=2\).

**Theorem.** The Smith entries of the Sylvester matrix of \(F,G\) are
the Smith entries of the Sylvester matrix of \(A,B\), with \(d\) extra
entries equal to one and \(d\) extra entries equal to zero. In
particular, if \(\Delta\) is the positive greatest common divisor of
the \((n+N-d)\)-minors of the first matrix, then
\[
                        \Delta=|\operatorname{Res}(A,B)|.       \tag{1}
\]
Neither \(\Delta\) nor the rest of the Smith data depends on \(H\)
when the cofactors are fixed and \(H\) remains primitive of degree
\(d\).

At a prime \(p\), let \(s_1,\ldots,s_a\) and \(t_1,\ldots,t_b\) be the
cofactor root valuations. Let the corresponding roots be
\(\alpha_i,\beta_j\) in a common algebraic closure of \(\mathbb Q_p\).
The exact valuation of (1) is
\[
 \begin{split}
 v_p(\Delta)
   ={}&\sum_{\substack{i,j\\s_it_j>0}}\min(|s_i|,|t_j|)\\
     &+\sum_{i,j}\bigl(v_p(\alpha_i-\beta_j)-\min(s_i,t_j)\bigr).
                                                               \tag{2}
 \end{split}
\]
The second sum is nonnegative and is supported on pairs with equal
valuations. The first sum is therefore a lower bound for \(v_p(\Delta)\)
obtained from the two parent Newton polygons and a proposed allocation
of the \(d\) common root valuations. It charges all overlaps of strictly
positive or strictly negative cofactor valuations, including unequal
valuations.

Let \(M\) denote Mahler measure. The global budget is
\[
       \Delta\, M(H)^{a+b}
               \ \le\ 2^{ab}M(F)^b M(G)^a.                     \tag{3}
\]
For an irreducible nonreal quadratic
\(H=h_2z^2+h_1z+h_0\), one has
\(M(H)=\max(|h_2|,|h_0|)\). If \(D\) is the reduced denominator of
the reciprocal trace, then \(D\le |h_2h_0|\le M(H)^2\). Thus (3)
also gives
\[
       \sum_p v_p(\Delta)\log p + \frac{a+b}{2}\log D
          \ \le\ b\log M(F)+a\log M(G)+ab\log2.                \tag{4}
\]
One can use local lower bounds from (2) in (4) without computing the
actual common factor. The tested root allocation must be the same in
these lower bounds and in any bound on \(D\). A proof at a diffuse
reference still needs a uniform estimate on the right side after this
debit.

## The integral Sylvester comparison

Use binary homogeneous forms of the declared degrees. In the
Sylvester map the last step is multiplication by \(H\), from the
lattice of forms of degree \(n+N-d-1\) to that of degree
\(n+N-1\). Its image is saturated. Indeed a form divisible by \(H\)
over \(\mathbb Q\) has an integral quotient by primitivity and Gauss's
lemma. Alternatively, after dehomogenizing, the product of a
primitive integer polynomial and a primitive rational quotient
numerator is primitive, so no denominator survives. The homogeneous
degree bound is preserved when rehomogenizing. The quotient lattice
is free of rank \(d\). The other map is
\[
 (U,V)\longmapsto AU+BV,\qquad
     \deg U=b+d-1,\quad \deg V=a+d-1,
\]
where the displayed degrees describe the homogeneous spaces. More
generally let \(\ell\ge0\) replace \(d\) in the multiplier and output
degrees. For \(\ell=0\) this is the ordinary cofactor Sylvester map.

Fix a prime \(p\). Over \(\mathbb Z_p\), make a finite unramified
extension of the coefficient ring with residue field large enough
that the nonzero reductions of \(A,B\) are simultaneously nonzero at
some projective point. Make that point infinity by an invertible
linear change of the homogeneous variables. The change of variables
induces invertible operations in every homogeneous degree. In this
chart, both leading coefficients are units and we may divide them
out. The use of a projective chart is necessary because the parent
leading coefficients need not be units. A change of affine chart
also appears in a different setting in
[Lindy and Noferini, §4.3](https://arxiv.org/html/2512.08550v1#S4.SS3).

In the new chart we are working over a discrete valuation ring
\(\mathcal O\), with monic ordinary polynomials \(A,B\).
For each \(\ell\) the cokernel of the cofactor map identifies with
\(\mathcal O[z]/(A,B)\). There is a surjection onto this quotient
because every class has an \(A\)-remainder of degree below \(a\).
Conversely, a polynomial \(C\) of degree at most \(a+b+\ell-1\) in
the ideal has a representation \(C=AU+BV\) with \(\deg V<a\), by
subtracting an \(A\)-multiple from \(V\). The resulting \(U\) has
degree at most \(b+\ell-1\), by monic division and the bound on
\(\deg(BV)\). This is exactly the required degree constraint.

The \(\ell=d\) cofactor map therefore has the same nonunit Smith
entries as the \(\ell=0\) map, with \(d\) extra unit entries
accounting for the larger target dimension. Both cofactor maps
have full row rank. Unramified extension preserves all valuations
of the Smith entries, so this identifies them primewise over
\(\mathbb Z\). The saturated multiplication by \(H\) adds a free
cokernel of rank \(d\). This gives the asserted Smith form and
(1).

## The valuation and height bounds

The Gauss norm identity gives
\(v_p(\operatorname{lc} A)=\sum_i(-s_i)_+\) for primitive \(A\),
and the analogous formula for \(B\). The root-product formula
applied to (1) therefore writes \(v_p(\Delta)\) as the sum over
\(i,j\) of
\[
        (-s_i)_++(-t_j)_++v_p(\alpha_i-\beta_j).
\]
On subtracting the nonnegative valuation excess in (2) from
each summand, the remaining expression is
\((-s_i)_++(-t_j)_++\min(s_i,t_j)\). It is zero if the two
valuations have opposite signs or either is zero, and
\(\min(|s_i|,|t_j|)\) otherwise. This proves (2). Further local
information can be obtained from [truncated Newton faces](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).

At the complex place,
\(|\alpha_i-\beta_j|\le
2\max(1,|\alpha_i|)\max(1,|\beta_j|)\). The same root-product
formula now gives
\(\Delta\le 2^{ab}M(A)^b M(B)^a\).
Multiplicativity of Mahler measure proves (3). For a nonreal
quadratic its roots have the same absolute value, whose square
is \(h_0/h_2\). Its reciprocal trace is \(h_1^2/(h_0h_2)-2\).
These identities prove (4).

## Positive parents with fixed Smith data and growing denominators

For every \(n\ge5\), put \(k=n-2\), and define
\[
 A=(z+1)^k+1,\qquad B=(z+1)A+1,\qquad
 H_t=t z^2-(t-1)z+t\quad(t\in\mathbb Z,\ t\ge2).
\]
Then \(F_t=H_tA\) and \(G_t=H_tB\) are primitive squarefree
polynomials of degrees \(n,n+1\) with strictly positive
coefficients. Their only common factor is the irreducible
quadratic \(H_t\).

To check positivity, set \(r=(t-1)/t\). For \(j\) in the
interior of the coefficient range and \(0\le j-1\le k\),
\[
 \frac{\binom{k}{j}+\binom{k}{j-2}}{\binom{k}{j-1}}
  = (k+1)\left(\frac1j+\frac1{k-j+2}\right)-2
  \ \ge\ 2-\frac4{k+2}>1>r.
\]
Thus \((H_t/t)(z+1)^k\) has positive coefficients.
The extra term for \(F_t/t\) is \(H_t/t\), which only
reduces the coefficient of \(z\). Its final value at that
index is \(k-2r>0\). For \(G_t/t\) the main term has
exponent \(k+1\) in place of \(k\) and the extra polynomial
is \((H_t/t)(z+2)\). The only possible negative extra
coefficient is at \(z\), where the final value is
\(k+2-3r>0\).

The negative discriminant makes \(H_t\) irreducible.
For a root \(z\) of \(H_t\) one has \(|z|=1\) and
\(\sqrt2<|z+1|<\sqrt3\), so neither \(A\) nor \(B\)
vanishes there. Indeed \(A\) would require \(|z+1|=1\),
and \(B\) would require \(|z+1|^{k+1}=|z+2|\), which is
impossible under these bounds. The polynomial \(A\) is
squarefree. A repeated root of \(B\), in the variable
\(u=z+1\), would satisfy \(u^{k+1}+u+1=0\) and
\((k+1)u^k+1=0\). The two equations force
\(u=-(k+1)/k\) and \(|u|<1\), a contradiction.
The identity \(B=1\) at the roots of the monic polynomial
\(A\) gives \(\operatorname{Res}(A,B)=1\).

The Smith form of every parent matrix is consequently the
identity of size \(n+N-2\) with two zero entries. The
reciprocal trace is \(((t-1)/t)^2-2\), of reduced denominator
\(t^2\) and limit \(-1\). The coefficient arrays after
division by \(t\) converge while these denominators grow.
For any fixed modulus \(M_0\), restricting to \(t\equiv1
\pmod{M_0}\) also fixes both integer parent arrays modulo
\(M_0\). The quadratic fields still need not be bounded:
for any new odd prime \(p\nmid M_0\), impose in addition
\(t\equiv p-1\pmod{p^2}\). The discriminant
\(-(3t-1)(t+1)\) then has valuation one at \(p\), which
ramifies in the quadratic field.

This example concerns positive parent arrays. It makes no
claim that they are linked factorial arrays. It shows that
even the full Sylvester Smith form, bounded coefficient
shapes and fixed finite residue data cannot by themselves
give a uniform height bound on the isolated quadratic.


The same obstruction can be placed inside the narrow angle range
\((\pi/n,2\pi/(n+2))\) for every \(n\ge5\).
Choose \(\theta\in(\pi/n,\pi/(n-1))\) so that
\(r=2\cos\theta=u/v\) is rational, with positive coprime \(u,v\),
and put \(H_*(z)=z^2-rz+1\).
The series for \(1/H_*\) through degree \(n-2\) has coefficients
\(\sin((j+1)\theta)/\sin\theta>0\).
Multiplying this truncation by \(H_*\) gives
\[
       1-\frac{\sin(n\theta)}{\sin\theta}z^{n-1}
         +\frac{\sin((n-1)\theta)}{\sin\theta}z^n,
\]
whose two nonconstant coefficients are positive. Give each
coefficient at \(z,\ldots,z^{n-2}\) a sufficiently small positive
value by adjusting the truncated coefficients successively.
The two highest coefficients stay positive by continuity.
The feasible cofactor set is open. Choose a rational squarefree
cofactor coprime to \(H_*\) in this set, then clear denominators
and divide out its content to obtain \(A\).

Take \(B=L(z+1)A+1\) with \(L\) a sufficiently large positive
integer. Then \(H_*B\) has positive coefficients, and \(B\) is
primitive because \(B(-1)=1\). The cofactors are coprime since
\(B\) takes value one at all roots of \(A\). There are only
finitely many exceptional \(L\) for squarefreeness and for
coprimality of \(B\) with \(H_*\). For squarefreeness this follows
by comparing \((z+1)A\) with its derivative at a repeated root
of \(B\). For the other assertion compare values at the two
roots of \(H_*\).

Now use \(H_t=vt(z^2+1)-(ut-1)z\) for positive integers \(t\)
divisible by \(v\) and tending to infinity. These quadratics
are primitive since \(\gcd(vt,ut-1)=1\). The factors \(H_t\)
have nonreal roots converging to those of \(H_*\), so eventually
the parents are squarefree and strictly positive. Their Smith
form is fixed by (1) and the full Smith assertion. The trace
denominator is exactly \((vt)^2\), and the roots stay in the
indicated angle range.
