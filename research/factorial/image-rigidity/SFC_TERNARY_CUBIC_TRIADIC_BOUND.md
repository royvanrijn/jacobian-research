# A uniform three-adic restriction on the shifted ternary cubic gate

For \(j\ge1\), set
\[
 F_j(A,B,C)=[z^j]\exp\!\left(j(Az+Bz^2/2+Cz^3/6)\right).
\]
These are the [noncollinear ternary limit forms](SFC_HOMOGENEOUS_TERNARY_SHIFTED_GATES.md#noncollinear-offsets).
Normalize \(v=v_3\) by \(v(3)=1\). If \(C\ne0\) and the three forms
at orders \(m,m+1,m+2\) vanish, let \(n\) be the unique multiple of
three in the window. Then \(A,B\ne0\) and, at each extension of the
three-adic valuation to the coefficient field,
\[
 0<v(A)-\tfrac12v(B)<\tfrac13,\qquad
 v(C)-\tfrac32v(B)\ge\tfrac12v(n).                            \tag{1}
\]
In particular the centered case \(A=0\) is excluded at every start, even
over an arbitrary characteristic-zero field. Rational triples are also
excluded, since \(v(A^2/B)\) would have to lie strictly between zero
and \(2/3\). The centered assertion applies over arbitrary
characteristic-zero fields by the weak Nullstellensatz applied to the
rational system with \(A=0\) and an auxiliary inverse for \(C\).
The inequalities leave other cubic parameters open.

## The lower bound and the cubic coefficient

Use weighted scaling \((A,B,C)\mapsto (\lambda A,\lambda^2B,\lambda^3C)\).
It preserves the zeros. Take \(C=1\) in the algebraic closure.
Write \(P_j(X,Y,Z)=j![z^j]\exp(Xz+Yz^2/2+Zz^3/6)\).
Its coefficients count labelled partitions into singletons, pairs,
and triples. The three-cycle on three marked labels gives
\[
 P_{j+3}(X,Y,Z)\equiv(X^3+Z)P_j(X,Y,Z)\pmod3 .
\]
Indeed, a fixed partition leaves the three labels as singletons or as
a single triple. If three blocks are permuted, each contains just one
of the marked labels and no fixed label, and if a block contains the
three marked labels it is already full. This is also the specialization
to maximum block size three of [Carlitz's congruence as stated by Gessel, (5.5)](https://people.brandeis.edu/~gessel/homepage/papers/combinatorial-congruences.pdf).

At \(n=3k\), divide \(n!F_n\) by \(n^k\). The pure triple term has
three-adic valuation zero, as \(v(n!)=k+v(k!)\). All other terms
are integral partition coefficients times powers
\((n^{2/3}A)^a(n^{1/3}B)^b\) with \(a+b>0\). They therefore
have positive valuation if \(v(A)>-2v(n)/3\) and
\(v(B)>-v(n)/3\). In particular, vanishing at \(n\) implies
\(\min(v(A),v(B)/2,0)<0\).

Scale so that this minimum becomes zero. The rescaled \(C\) now has
positive valuation. At the order \(j\equiv1\bmod3\), the partition
congruence for \(P_j(jA,jB,jC)\) forces the rescaled \(A\) to
have positive valuation. Thus, in the normalization \(C=1\),
\(B\ne0\) and \(v(A)>v(B)/2\), allowing \(v(A)=+\infty\).
If \(v(B)>-v(n)/3\), then \(v(A)>-2v(n)/3\) too, which was
excluded by the pure triple term. Undoing the weighted scaling proves
the lower bound for \(v(A)-v(B)/2\) and the bound for \(C\) in (1).

## Reduction by replacing singletons

Normalize \(B=1\), write \(\alpha=v(A)\), \(\gamma=v(C)\), and
suppose \(\alpha\ge1/3\). The case \(A=0\) is included by taking
\(\alpha=+\infty\). At start one the equations \(F_1=F_2=0\)
force \(A=B=0\), so assume \(m\ge2\).
At any order \(j\), let \(e=v(j)\). For \(a+2b+3c=j\), the
valuation of the corresponding term in \(F_j\) is
\[
 a(\alpha+e)+b e+c(\gamma+e-1)-v(a!)-v(b!)-v(c!).            \tag{2}
\]
Replace \(2r\) singletons by \(r\) pairs, where \(r=\lfloor
a/2\rfloor\). The excess valuation of the original term is at least
\[
 2r\alpha+re-v(a!)+v(r!).
\]
For \(a\ge2\), this is positive except possibly at unit orders
with \(a=3,9\). To check the bound, count multiples of \(3^i\)
in \((r,2r+1]\). With \(L=\lfloor\log_3(2r+1)\rfloor\), their
total for \(r\ge18\) is at most
\((r/2)(1-3^{-L})+L<2r/3\), since \(L\le r/6\).
At \(1\le r\le17\) the same count is
\[
 1,1,1,3,3,3,4,4,4,5,5,5,8,8,8,9,9 ,
\]
and for \(a=2,8\) the count ends at \(2r\) and is respectively
zero and one. This proves \(v(a!)-v(r!)<2r/3\) except for
\(a=3,9\). When \(e>0\), both exceptions are strictly lowered
by the same replacement since \(v(a!)\le r\).

We now show that \(\gamma>1\) is impossible. If an even order
\(j\) in the window is prime to three, compare with the pure pair
term. The number of removed pairs is \(R=(a+3c)/2\), so the
pair factorial contributes at least \(v(R!)\). For \(a=0,1\)
this offsets \(v(c!)\). For \(a=3,9\) the parity implies \(c\ge1\).
At \(a=3\) the excess valuation is at least
\(3\alpha-1+c(\gamma-1)>0\). At \(a=9\), \(R-c\ge5\)
gives \(v(R!)\ge v(c!)+1\), so the excess is at least
\(9\alpha-3+c(\gamma-1)>0\). Thus the pure pair term is
uniquely least.

Otherwise there are two odd unit orders around the even multiple
of three. At the one with \(j\equiv1\bmod3\), the replacement
for \(a=3,9\) increases the pair count to \(b'\equiv0\bmod3\).
The removed strip from \(b'!\) then contributes at least one
for length one and two for length four. The two exceptional
terms are therefore strictly above their replacements by (2).
At fixed \(a\in\{0,1\}\), the least admissible number of triples
is \(c_*=1-a\). Increasing it by \(2r\), \(r\ge1\), removes
\(3r\) pairs. This raises (2) since \(\gamma>1\) and
\(v((c_*+2r)!/c_*!)\le v((3r)!)\).
For finite \(\alpha\) the only possible tie at this \(j\) requires
\[
 \gamma=\alpha+1-v((j-1)/2)\le\alpha .
\]
This excludes \(\gamma>1\) for \(\alpha\le1\). If \(\alpha>1\)
then pair replacement lowers all terms with \(a\ge2\), including
at the odd order \(j'\equiv2\bmod3\). The same comparison there
requires \(\gamma=\alpha+1\), incompatible with the first tie.
When \(A=0\) only the \(a=0\) branch exists and its minimum is
unique. Thus \(\gamma\le1\).

For \(\gamma\le1\), replace three or nine singletons by one or
three triples, respectively. The respective excess valuations
in (2) are at least \(3\alpha+2e-\gamma\) and
\(9\alpha+6e-3\gamma\). The latter bound uses that the three new
consecutive triple indices contribute at least one to their factorial
valuation. The replacement is strict unless \(\alpha=1/3\),
\(\gamma=1\), and \(e=0\). At \(j\equiv0,2\bmod3\), every
\(a=1\) term has \(b\) a three-adic unit. Replacing its singleton
and a pair by a triple is strict, with excess in (2) equal to
\(\alpha+e+1-\gamma+v(c+1)-v(b)>0\).

## Opposite bounds on the centered Newton slopes

For \(j\ge2\), choose \(b_0\in\{0,1,2\}\) with \(2b_0\equiv
j\bmod3\), and set \(c_0=(j-2b_0)/3\), \(K=\lfloor c_0/2\rfloor\),
and \(s=c_0-2K\). After a nonzero monomial factor, the centered
polynomial is a polynomial in \(W=B^3/C^2\), of degree \(K\).
For \(b=b_0+3k,\ c=c_0-2k\), its coefficient valuation relative
to its constant coefficient is
\[
 U_k=(1+v(j))k+v\bigl(c_0!/(c_0-2k)!\bigr)-v(k!).
\]
This follows by taking the consecutive coefficient ratio and using
\(v((b_0+3k)!/b_0!)=k+v(k!)\).

At \(3\mid j\), all slopes of the lower coefficient polygon are
at least two. Indeed \(v(c_0!/(c_0-2k)!)\ge v((2k)!)\ge v(k!)\).
Equality in \(U_k\ge2k\) at \(k>0\) requires \(v(j)=1\), \(k=1\),
and \(c_0\equiv2\bmod3\), because \(v((2k)!/k!)>0\) when \(k\ge2\).
In particular slope two can occur only under these necessary conditions.

At \(3\nmid j\), all slopes are at most two. For \(1\le r\le K\),
\[
 U_K-U_{K-r}
 =r+v\bigl((s+2r)!/s!\bigr)-v\bigl(K!/(K-r)!\bigr)\le2r .
\]
Here \(v((s+2r)!)\le r\). Equality requires \(s=1\),
\(2r+1\) a power of three, and a length-\(r\) factorial strip
consisting of units. The strip condition gives \(r\le2\), so
\(r=1\) and \(3\nmid K\). In particular slope two requires
\(c_0\) odd and \(3\nmid K\).

A sum can vanish only if its least valuation is attained at least
twice. When only centered terms can minimize (2), their polygon
must therefore have slope \(-v(W)=2\gamma\). At the multiple of
three only centered terms can minimize even when \(\alpha=1/3\).
The lower slope bound and \(\gamma\le1\) force \(\gamma=1\),
slope two, and \(n/3\equiv2\bmod3\).

First let \(\alpha>1/3\). At the order equal to two modulo
three, only centered terms can minimize (2), so it too must
have slope two. Write \(n=3l\). If the other order is \(3l-1\),
its \(c_0=l-1\) must be odd, but then \(K=(l-2)/2\) is
divisible by three. The only remaining placement starts at
\(3l\) and contains \(3l+2\), whose \(c_0=l\) is odd.
The remaining \(3l+1\) is an even unit order. There the
\(a=1\) terms are strictly above pure pairs by the even-order
comparison at \(\gamma=1\), and its centered \(c_0=l-1\) is
even. Its centered polygon has no slope two.

Finally let \(\alpha=1/3\). At any even unit order, the
least centered term is unique at \(c=0\), since the centered
\(c_0\) is even and all slopes are strictly less than two.
An exceptional replacement of three or nine singletons by
triples cannot produce this term. Terms with \(a=1\) are
strictly above pure pairs as before. Thus there can be no
even unit order. The multiple of three \(n=3l\) must be the
even middle order, and \(j=n+1\) is an odd unit order.
Its centered \(c_0=l-1\) is odd, but \(K=(l-2)/2\) is
divisible by three. Its uniquely least centered term
therefore has \(c=1\). The \(a=1,c=0\) term lies strictly
below it, since \(v((j-1)/2)=v(3l/2)=1\), and the excess
valuation over the \(a=0,c=1\) term is \(\alpha-1\).
Within the \(a=1\) branch, increasing \(c\) from zero to
\(2r\) removes a strip of \(3r\) pairs. That strip has
valuation at least \(r>v((2r)!)\), so those terms are higher.
The exceptional \(a=3,9\) terms are at least as high as
centered terms, and all other replacements were strict.
This leaves a unique least term and excludes the boundary.
