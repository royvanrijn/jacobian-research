# Regularity along the cube divisor in sextic Veronese pencils

**Theorem SexticKellerThreeCubeVeroneseExclusion.**
There is no characteristic-zero three-variable Keller map with
highest tuple \((P^2,Pt^3,t^6)\) in any of the following three
frames:
\[
\begin{array}{c|c|c}
   P&t&\text{source weights}\\ \hline
   x^2z+y^3 & x &(1,1,1)\\
   x^2z     & y &(1,1,1)\\
   x^2y+z^2 & x &(2,2,3).
\end{array}
\]

The argument tracks the two summands \(tS,t^2S\) of the polynomial
constant ring through the same binomial regularity calculation,
then uses the small source weights of the remaining tangent
parameters to finish the descent.

**Theorem SexticKellerCubeVeroneseDivisorRegularity.**
In those three frames, rescale the components by their weights,
write \(A=F_{1,\epsilon}^{1/2}\) for the formal root with leading
term \(P\), and put \(B=F_{2,\epsilon}/A\).
Fix \(A\) as a generic parameter on a formal chart at \(t=0\),
with a free fibre variable and, if needed, algebraic extension
of the scalar field \(k(A)\).
There is no Keller completion if, after subtracting a series in
\(\epsilon,A\) independent of the fibre variables and with zero
constant in \(\epsilon\), the minimum-weight part of \(B\)
on the chart is
\[
                   t^3+\epsilon^a t^b v
\]
with \(v\) transcendental over the scalar field and
\((a,b)\in\{(1,0),(1,1),(2,0)\}\) in any of the three frames,
or \((a,b)\in\{(2,1),(4,0)\}\) in the weighted frame. The
weights of \(t,\epsilon\) in this comparison are \(a,3-b\).

**Theorem SexticKellerCubeVeroneseAlgebraicRegularity.**
Use the same notation, with no subtraction from \(B\).
There is no Keller completion if the minimum-weight part of
\(B\) is \(t^3+\epsilon^a t^b v\) with
\(v^p\in k(A)^\times\), \([k(A)(v):k(A)]=p\), and
\((a,b,p)\) one of
\[
 \begin{array}{c|c}
 \text{ordinary weights}&(1,0,3),\ (1,1,3),\ (2,0,3)\\
 (2,2,3)&(1,1,2),\ (3,0,2).
 \end{array}
\]

**Theorem SexticKellerOrdinarySplitCubeExclusion.**
There is no characteristic-zero three-variable Keller map with
highest ordinary tuple \((P^2,Pt^3,t^6)\) for \(P=x^2z,\ t=y\).

## The constant ring and the reference tuples

Work over an algebraic closure \(k\). Let \(d\) be the weight
of \(t\), so that \(P,Q=t^3\) both have weight \(3d\) and
the leading components have weight \(6d\). The rational
kernel of \(\det J(P,t,\cdot)\) is \(k(P,t)\) in all three
frames. In the first and third frames fix \(t=x\ne0\) and
solve for \(z\) or \(y\), respectively, and in the second fix
\(t=y\) and solve for \(z\) on \(x\ne0\). In each case one
free coordinate remains. The pair \((P,t)\) is onto over
\(k\), including at \(t=0\), where \(P\) restricts to \(y^3\),
\(x^2z\), or \(z^2\), respectively. Polynomial contraction
therefore gives the kernel
\[
             T=k[P,t]=S\oplus tS\oplus t^2S,\qquad
             S=k[P,Q].
\]
Multiplication of the derivation by \(3t^2\) does not change
this kernel. The leading adjugate calculation is the one used
for the [closed-kernel Veronese descent](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md#cubic-veronese-tuples-with-a-closed-constant-ring).
See also [Wang–Yan, Section 3.4](sources/wang-yan-2610.36495v1.pdf)
for the related quartic computation.

Write \(n=(Q^2,-2PQ,P^2)\). For a polynomial reference tuple
factoring through two source polynomials, with the given
weight-\(6d\) tuple, the first positive error \(E_j\) of weight
\(j<6d\) has \(n\cdot E_j\in T\). This follows from the
determinant coefficient \(2\det J(P,Q,n\cdot E_j)\) of
weight \(12d+j-\sum w_i>0\), which must vanish. If the normal
contraction is zero, the tangent error has form
\[
                     (2P U, Q U+P V,2Q V)
 \tag{1}
\]
with \(U,V\) of weight \(j-3d\), by coprimeness of \(P,Q\).
Adding \(U,V\) to the two parameter polynomials of the reference
absorbs (1) at weight \(j\).

When the normal contraction belongs to \(S\), its
weight \(6d+j\) forces \(j=3d\) unless the contraction is
zero. At \(j=3d\) it is a cubic in \(P,Q\), the contraction
with \(n\) of a linear vector in \(P,Q\). Consequently one may
descend with references \(G(\widetilde P,\widetilde Q)\), where
the quadratic part of \(G\) is \((p^2,pq,q^2)\) and its other
terms are affine in \(p,q\). If the descent could continue at
every positive weight, absorption of constants would give
an actual two-parameter factorization and zero determinant.
There must therefore be a first error with
\(n\cdot E_j\in T\setminus S\).

For a reference factoring through \(\widetilde P,t\) one
can instead absorb an arbitrary normal contraction in \(T\).
Indeed every monomial \(P^i t^l\) of weight at least \(6d\)
satisfies \(i+\lfloor l/3\rfloor\ge2\), so it belongs to
\((P^2,PQ,Q^2)T\). Absorb the normal error by replacing
the required monomials of \(P\) by those of \(\widetilde P\)
in the reference. The remaining error is tangent as in (1).
To preserve the two-parameter reference it suffices that
\(V\in k[t]\). More explicitly one keeps the reference of the
form \(\nu(\widetilde P,\widetilde Q)+L(\widetilde P,t)\),
where \(\nu=(p^2,pq,q^2)\), \(L\) has lower weight, and
\(\widetilde Q\in k[t]\). The first change in the reference
from adding \(U,V\) to \(\widetilde P,\widetilde Q\) is (1).

## Binomial regularity on a transverse edge

Let \(K\subset L\) be characteristic-zero fields and
\(v\in L\) transcendental over \(K\). Suppose \(B'\) is
regular in \(t,\epsilon\) over \(L\), with \(B'(t,0)=t^3\),
and its minimum-weight part is \(t^3+\epsilon^a t^b v\)
in weights \(a,3-b\), where \(a>0\) and \(0\le b<3\) are
integers. In the \(t\)-Laurent series ring truncated modulo
\(\epsilon^N\), use the formal root \((B')^{1/3}=t+O(\epsilon)\).
Let
\[
 R(s,\epsilon)=\sum_{i=0}^{N-1}\epsilon^i R_i(s),
 \qquad R_i(s)\in K((s)).
\]
If \(R((B')^{1/3},\epsilon)\) is regular in \(t\) modulo
\(\epsilon^N\), then every monomial \(\epsilon^i s^l\) of
\(R\) with \((3-b)i+a l<(3-b)(N-a)\) has
\(l\in3\mathbb Z_{\ge0}\).

In fact, suppose there is an offending monomial, and take
the least offending weight \(W\). Subtract all terms of
weight at most \(W\) with nonnegative exponents divisible
by three. Their substitutions are regular powers of \(B'\).
Other terms of weight above \(W\) cannot contribute at
weight \(W\), because the binomial expansion of \(B'/t^3\)
has no negative weights. The weight-\(W\) terms from the
remaining monomials come just from \(t^3+\epsilon^a t^b v\).
For an offending coefficient \(c_{il}\ne0\) of
\(\epsilon^i s^l\) of weight \(W\), choose the first index
\(j\ge i\), with \(j\equiv i\pmod a\), for which
\((3-b)j>W\). The weight bound guarantees \(j<N\).
At that index put \(n_0=(j-i)/a\). The coefficient of
\(\epsilon^j t^{(W-(3-b)j)/a}\) contains the nonzero term
\[
                   c_{il}\binom{l/3}{n_0}v^{n_0}.
\]
The binomial is nonzero because \(l/3\) is either a
negative integer or a noninteger. All other terms of
weight \(W\) at that index contribute distinct powers
of \(v\), because their \(\epsilon\)-indices differ.
They cannot cancel it over \(K\). The exponent of \(t\)
is negative, contradicting regularity.

## Applying the regularity calculation to the descent

Under rescaling the determinant is a nonzero scalar times
\(\epsilon^N\), where \(N=18d-\sum w_i\), thus \(15\)
for ordinary weights and \(29\) for the weighted pencil.
A chart fixing \(A\) exists at the generic point of \(t=0\):
start with \(P=y^3,x^2z,z^2\) there in the three cases
and solve formally for \(y,z,z\), respectively. In the
first and third cases fix a root of the corresponding
scalar equation over \(k(A)\); in the second use \(x\ne0\).
The formal derivative of \(P\) in the solved variable
is a unit on each chart.

Assume the edge in the theorem occurs for \(B'=B-\beta\),
where \(\beta\) is the series subtracted in the theorem.
For contradiction take the first non-\(S\) normal error
of the descent through references \(G(\widetilde P,\widetilde Q)\)
above. Put \(r=6d-j\) for its loss. In the rescaled
reference, solve its first two equations for the two
parameters as formal functions of \(A,B\).
At \(\epsilon=0\) the equations are \(p^2=A^2,\ pq=AB\)
and the derivative determinant in \(p,q\) is \(2A^2\).
Thus the reference third component is a formal function
\(\Phi_\epsilon(A,B)\) with coefficients polynomial in
\(B,A,A^{-1}\). The actual difference
\[
                         C=F_{3,\epsilon}-\Phi_\epsilon(A,B)
\]
is regular in \(t\) on the chart. It begins at order
\(\epsilon^r\) with coefficient \((n\cdot E_j)/P^2\),
by the first variation of \(F_2^2/F_1\) at the top tuple.

Through order \(\epsilon^{N-1}\), \(C\) is a formal series
in \(A\) and \((B')^{1/3}\), with scalar Laurent series
coefficients in the latter on the fixed-\(A\) chart.
One way to see this first with \(B\) in place of \(B'\)
is to subtract any coefficient already matched by a
rational function of \(A,B^{1/3}\).
At the next coefficient its differential wedged with
\(dP,dt\) vanishes by the Keller determinant, since
the formal differentials of \(A,B^{1/3}\) span the
same plane as those of \(F_{1,\epsilon},F_{2,\epsilon}\)
over the generic field. The coefficient therefore
belongs to \(k(P,t)\) and can be lifted by replacing
\(P,t\) with \(A,B^{1/3}\).
This gives the series modulo \(\epsilon^N\).
Now substitute \(B=B'+\beta\) and expand in \(\epsilon\)
with \(B'\) invertible in the Laurent chart.
This leaves the first coefficient at order \(r\) unchanged,
and yields the claimed \(R((B')^{1/3},\epsilon)\).
The coefficient field can be extended algebraically
to contain the scalar functions used in the chart
and in \(\beta\).

Since \(n\cdot E_j\) is in \(T\setminus S\),
some monomial in its quotient by \(P^2\) has \(t\)-exponent
\(l\) not divisible by three. It has
\(0<l\le(12d-r)/d\).
Its weight \((3-b)r+a l\) is less than \((3-b)(N-a)\)
in each applicable pair \((a,b)\), since
\(1\le r\le6d-1\) and \((d,N)=(1,15)\) or \((2,29)\).
This contradicts binomial regularity.

The same weight comparison also gives the algebraic test
when \(\beta=0\). In that case the coefficients of \(R\)
are over \(K=k(A)\). Fix a weight \(W\) after subtraction of
the safe terms \(s^{3j}\), \(j\ge0\), at weights up to \(W\).
If a least offending weight is negative, its first
\(\epsilon\)-index already contributes a negative \(t\)-power
with nothing to cancel it. Thus it suffices to test
\[
                   0\le W\le (3-b)r+a(12-r/d),
 \tag{2}
\]
where \(r\) is the loss of the first non-\(S\) normal error.
The normalized series starts at \(\epsilon^r\). Also
\(r\) must be divisible by \(d\) but not by \(3d\), because
the normal weight \(12d-r\) belongs to \(T\setminus S\).
For each residue class \(\rho\bmod ap\), take as columns
the indices \(i\) in \([r,N-1]\) congruent to \(\rho\)
for which \(l_i=(W-(3-b)i)/a\) is an integer outside
\(3\mathbb Z_{\ge0}\). Take rows \(n\) in \([0,N-1]\) with
\((3-b)n>W\) and \(n\equiv\rho\pmod a\).
The coefficient of \(\epsilon^n\) at weight \(W\) supplies
the row with entries
\[
              \binom{(W-(3-b)i)/(3a)}{(n-i)/a},
 \tag{3}
\]
where an index below zero gives entry zero.
Indeed the factor \(v^{(n-i)/a}\) splits as
\(v^{(n-\rho)/a}(v^p)^{(\rho-i)/(ap)}\).
For a fixed row the distinct residue classes modulo \(ap\)
give independent powers of \(v\) modulo \(p\), and the
second factor amounts to nonzero column scaling over \(K\).

The [exact rank checker](verify_sextic_cube_veronese_regularity.py)
enumerates (2) and row-reduces (3) over the rational numbers.
In the order of the five cases in the theorem it checks
\(228,192,334,176,384\) nonempty matrices, all of full
column rank. Thus no offending coefficient occurs through
the bound (2). In particular the coefficient
\((n\cdot E_j)/P^2\) could not contain any power of \(t\)
outside the \(S\)-summand. This contradiction proves
algebraic regularity.

## First corrections in ordinary weights

In the ordinary pencils, the first normal correction
is represented by \(t^2 M(P,Q)\) with \(M\) linear,
as the only degree-five monomials in \(T\) are \(t^5,t^2P\).
The remaining degree-five error is (1).
On a chart fixing \(A\), the order-one coefficient of
\(B-t^3\) agrees modulo \(t^2\) with \(V\) restricted
to the leading \(P\)-level.
This follows by expanding \(F_{2,\epsilon}/A\):
the contribution involving \(F_{1,\epsilon}\)
is divisible by \(t^3\), and the first reference normal
correction to \(F_{2,\epsilon}\) is divisible by \(t^2\).
The change from fixing \(P\) to fixing \(A\) does not
affect the order-one coefficient, because \(B_0=t^3\)
depends only on \(t\).

## Complete descent for \(P=x^2z,\ t=y\)

At the fixed-\(P\) level, \(x\) is free and \(z=P/x^2\).
A quadratic form in \(x,z\) independent of the free
variable is zero: its three terms have \(x\)-powers
\(2,-1,-4\). Thus binomial regularity with
\((a,b)=(1,0)\) first forces \(V\bmod y=0\).
A nonzero linear combination of \(x,z\) also varies
with \(x\), so \((a,b)=(1,1)\) then forces the
coefficient of \(y\) in \(V\bmod y^2\) to vanish.
Consequently \(V\) is a scalar multiple of \(y^2\),
and the first correction is absorbed in a reference
factoring through \(\widetilde P,y\) as above.

At weight four, absorb the normal contraction into
a vector of weight four in \(T\), whose entries are
linear combinations of \(y^4,yP\), both divisible by \(y\).
The remaining tangent corrections \(U,V\) are linear.
In the \(A\)-chart the coefficient \(B_2(0)\) is
\(V|_{y=0}\): the reference second component vanishes
on \(y=0\) through this order, the first normal
correction was divisible by \(y^2\), and the new normal
correction is divisible by \(y\). With weights
\(2,3\) on \(t,\epsilon\), the first corrected
coefficient of \(B\) has weight greater than six
including its factor of \(\epsilon\), since it is
divisible by \(y^2\). Hence a nonzero \(V|_{y=0}\)
gives the edge \((a,b)=(2,0)\) with transcendental
coefficient. This is impossible, so \(V\in k[y]\).
The second correction is absorbed into the same
two-parameter reference.

At weight three both tangent parameters are constant,
and at weights two and one the tangent part is zero.
At every remaining weight the normal contraction can
be absorbed in \(T\), and finally one absorbs the
constant error. This factors the original map through
\(\widetilde P,y\), contradicting its nonzero determinant.

## Complete descent for \(P=x^2z+y^3,\ t=x\)

On a fixed-\(P\) chart at \(t=0\),
\(y=y_0+O(t^2)\), \(y_0^3=P\), and \(z\) is free.
A \(yz\) or \(z^2\) term of the first \(V\) makes
\(B_1(0)\) transcendental, and is excluded.
If a \(y^2\) term remains, \(B_1(0)\) is a nonzero
multiple of \(y_0^2\), of degree three over \(k(P)\)
with cube in \(k(P)\). Apply the algebraic edge
\((1,0,3)\). Thus \(B_1(0)=0\).
An \(xz\) term of \(V\) now gives a transcendental
coefficient of \(t\) in \(B_1\), excluded by \((1,1)\).
An \(xy\) term on its own gives a nonzero multiple
of \(y_0\) instead and is excluded by \((1,1,3)\).
Here the derivative in \(t\) of \(y\) at fixed \(P\)
vanishes at \(t=0\), so no extra term enters any of
these coefficients. The first \(V\) is therefore
a scalar multiple of \(t^2\).

At weight four the same normal absorption as for
\(P=x^2z,\ t=y\) gives \(B_2(0)=V|_{t=0}\)
for the new linear tangent parameter \(V\).
A nonzero \(z\) coefficient would give a
transcendental edge \((2,0)\).
With that excluded, a nonzero \(y\) coefficient
gives the algebraic edge \((2,0,3)\).
The first correction \(B_1\) is divisible by \(t^2\),
so neither of these edges is altered by order one.
Thus the new \(V\) is a multiple of \(t\), and the
same two-parameter descent finishes at weights
three, two and one.

## Complete descent for \(P=x^2y+z^2,\ t=x\)

At weight eleven the normal contraction vanishes,
because \(T\) has only even weights. The tangent
parameter \(V\) in (1) is a linear combination of
\(xz,yz\). On the leading fixed-\(P\) chart
\(z=z_0+O(t^2)\), \(z_0^2=P\), and \(y\) is free.
The \(yz\) term would give a transcendental value
of \(B_1(0)\), excluded by \((1,0)\).
A nonzero \(xz\) term then gives a coefficient
of \(t\) in \(B_1\) of degree two over \(k(P)\),
excluded by \((1,1,2)\). Therefore \(V=0\).
A source change \(z\mapsto z+\lambda x+\mu y\)
removes \(U\), which is a linear combination
of \(xz,yz\) as well. Thus \(F_{11}=0\) for the
whole tuple and \(B_1=0\).

At weight ten a normal correction has entries in
\(\langle t^5,t^2P\rangle\). The tangent \(V\) has
weight four and is a binary quadratic in \(x,y\).
As at the first correction in the ordinary pencils,
\(B_2\) on the fixed-\(P\) level agrees with \(V\)
modulo \(t^2\). A \(y^2\) term in \(V\) gives a
transcendental edge \((2,0)\), so is excluded.
To test an \(xy\) term, consider the next weight
nine. Since \(F_{11}=0\), its determinant equation
uses only the homogeneous top tuple and \(F_9\)
at loss three. Its normal contraction has odd
weight and vanishes. Hence \(F_9\) itself is
tangent, with its parameter \(V_3=c z\) in
weight three. Expanding the quotient gives
\(B_3=V_3\) at this order. In this calculation fixing
\(A\) changes the source starting at order two
and does not change \(B_3\), since \(B_1=0\)
and \(B_0=t^3\).

Now subtract \(\epsilon^3 c z_0\) from \(B\) on
the fixed-\(A\) chart. An \(xy\) term in the
weight-four tangent parameter gives the edge
\((2,1)\), with coefficient a nonzero multiple of
the free variable \(y\). Indeed \(B_1=0\),
\(B_2\) is divisible by \(t\), and the constant
term of \(B_3\) has just been subtracted. This
edge is impossible. Thus the weight-four tangent
parameter is a multiple of \(t^2\) and can be
absorbed into a reference factoring through
\(\widetilde P,t\). In particular \(B_2\) is
divisible by \(t^2\).

Return to the weight-three tangent parameter
\(V_3=c z\). If \(c\ne0\), \(B\) itself has the
algebraic edge \((3,0,2)\), since its order-two
coefficient is divisible by \(t^2\) and its
order-three constant in \(t\) is \(c z_0\).
This excludes \(c\ne0\); thus \(B_3=0\) and the
weight-nine correction is absorbed by changing
only \(\widetilde P\).

At weight eight the new normal vector has entries
in \(\langle t^4,tP\rangle\), both divisible by \(t\).
The new tangent \(V\) is linear in \(x,y\), and
\(B_4(0)=V|_{t=0}\). To see this last identity one
can evaluate the second component at \(t=0\):
the reference vanishes there through loss four,
and the first three terms of \(B-t^3\) vanish at
\(t=0\), so fixing \(A\) has no effect on \(B_4(0)\).
A \(y\) term in \(V\) now gives the edge \((4,0)\)
with transcendental coefficient. All earlier terms
have weight above the edge or vanish, since
\(B_1=B_3=0\) and \(t^2\mid B_2\).
Thus \(V\) is a multiple of \(t\) and is absorbed.

At weight seven there are no tangent parameters
of weight one, and at weight six the tangent
parameters are constant. At all smaller positive
weights the tangent part is zero. Absorbing all
remaining normal corrections in \(T\) and the error
of weight zero yields a factorization through two
source polynomials, again contradicting the Keller
condition.


**Corollary SexticKellerOrdinaryCubicPencilsAndFourthPowerReduction.**
In the independent sextic leading-pencil classification, the
case of a closed pencil of ordinary degree three is impossible.
In the line-pencil case with quadratic outer map and binary quartic
common factor, a triangular source change gives a map of ordinary
degree at most six and sextic leading rank at most one.

The [closed-kernel and binary quartic reductions](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md)
leave just the two ordinary pencils and the \((2,2,3)\) pencil
of SexticKellerThreeCubeVeroneseExclusion. Combining their
statements with that exclusion gives the two conclusions.
