# Shallow dyadic windows for two prime-power anchors

Let \(K\) carry a real valuation \(v\) extending \(v_2(2)=1\).
For a scalar \(\tau\ne0\) and a sequence \(U_m\) in \(K\)
such that every \(U_m\) is congruent to an odd rational integer
modulo \(8\mathcal O_K\), put
\[
                  S_n(z)=\sum_{m=0}^n(n)_{\underline m}\,m!\,U_m z^m .
\]
For every \(n\ge1\) and every \(v(\tau)<31/16\), at least one of
\(S_n(\tau^{-1}),S_{n+1}(\tau^{-1})\) is nonzero.
There is no restriction on the ramification or residue field
within this valuation range. The proof gives a finite exact
replay [of the digit and derivative tables](verify_sfc_dyadic_shallow_windows.py).

In particular this applies to \(F=t+x^{2^b}y^{2^c}\), since
on writing \(\omega=2^b+2^c-2\) the quotient
\[
       U_m=\frac{(2^bm)!(2^cm)!}{2^{\omega m}(m!)^2}
\]
is an odd rational integer and
\(\Lambda(F^n)/t^n=S_n(2^\omega/t)\).
Combine it with the
[deep two-window theorem for two prime-power anchors](SFC_HIGHER_DEGREE_CONSTANT_SHIFT.md#5-two-window-reduction-for-two-prime-power-anchors),
which applies for \(v(t)-\omega\ge2\lambda=2\) at \(p=2\).
A common zero must therefore have
\[
                         31/16\le v(t)-\omega<2.                       \tag{1}
\]
Consequently every two-term phase \(t+x^{2^b}y^{2^c}\) over
a number field with a place over \(2\) of ramification index
at most fifteen has Strong Factorial at every order, for
arbitrary \(b,c\ge0\). For complex \(t\), an exact common
zero is necessarily algebraic since the moments are nonzero
polynomials in \(t\) with integer coefficients. A common root
would have ramification index at least sixteen at every place
over \(2\), and degree at least sixteen over \(\mathbb Q\).
The same two-window conclusion in the shallow scalar range
holds for other phases whose normalized source moments give
the stated \(U_m\bmod 8\) condition.

The proof uses Kummer's count of borrows and Legendre's
factorial valuation as in
[Granville, Section 2](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node3-an.shtml).
The fixed number of low digits comes from the scalar slope gap.
At the finitely many possible common leading residues, the next
derivative terms separate the windows. Initial residues alone
leave a common factor in these cases.

## 1. A first-block reduction for factorial kernels

The reduction works at any prime. For an integer \(r\ge1\),
consider \((n)_{\underline m}(m!)^{r-1}U_m/\tau^m\)
with \(U_m\) units and \(0\le m\le n\). Write
\(\lambda=1/(p-1)\) and \(\gamma=r\lambda-v(\tau)>0\).
The valuation of the term is
\[
                 w_n(m)=\gamma m-r\lambda s_p(m)+B(n,m),
\]
where \(B(n,m)=v_p\binom nm\) counts the borrows in
the subtraction of \(m\) from \(n\).
Let \(P=p^H\), and write \(m=P M+b\), \(0\le b<P\).
The low digits and the borrow into position \(H\) agree
for subtraction of \(m\) and of \(b\). In the remaining
positions, subtraction of \(b\) uses zero input digits,
so its outgoing borrow never exceeds the borrow from \(m\).
This follows recursively because outgoing borrow is
nondecreasing in the subtracted digit and incoming borrow.
Thus
\[
                   w_n(PM+b)-w_n(b)\ge(\gamma P-r\lambda)M.          \tag{2}
\]
Here \(b\le n\) whenever \(m\le n\), and \(s_p(M)\le M\).
If \(\gamma P>r\lambda\), all minimizing indices lie in the
first block, uniformly over every \(n\).

There is a corresponding bound for differentiation by
\(\theta=z\,d/dz\). For a gauge \(g\) with \(0\le g<P\),
the coefficient multiplier of \(\theta-g\) is \(m-g\).
If \(b\ne g\), then \(v_p(m-g)=v_p(b-g)\) and (2) carries over
to that derivative term. If \(b=g\) and \(M>0\), its
derivative valuation is at least
\[
                      w_n(g)+(\gamma P-r\lambda)M+H.                 \tag{3}
\]
This treats the higher block even if the derivative term
vanishes at its counterpart \(g\). These estimates allow
finite precision beyond the leading residue without a size
bound on the ordinary order.

## 2. The finite leading-residue comparison

Return to \(p=2,r=2\). Put \(\beta=v(\tau)\).
For \(\beta<1\), the even orders suffice: for even \(n\),
\(v_2((n)_{\underline m})\ge\lceil m/2\rceil\) and
\(v_2(m!)\ge\lfloor m/2\rfloor\), so the constant term
is the unique minimum. First handle
\(1\le\beta<15/8\).

Write \(n=P N+r_0\), with \(0\le r_0<P=2^H\), and
let \(s=v_2(N)\) when \(N\ne0\). For \(m<P\), using
\(\gamma=2-\beta\), the valuation is
\[
 w_n(m)=\gamma m+c(m),\quad
 c(m)=
 \begin{cases}
   s_2(r_0-m)-s_2(r_0)-s_2(m), &m\le r_0,\\
   s_2(P+r_0-m)-s_2(r_0)-s_2(m)+s-1,&m>r_0 .
 \end{cases}                                                       \tag{4}
\]
The second case uses \(s_2(N-1)-s_2(N)=s-1\) and is
absent when \(N=0\). To determine minima, cap \(s\) at
\(2H\) and use this cap for \(N=0\) as well. Indeed the
upper borrow carries incurred when \(m>r_0\) give
\(B(n,m)\ge s\). Thus if \(s\ge2H\), the second case
has valuation at least \(\gamma m-2s_2(m)+s>0\) and
cannot be a minimizing term, since \(m=0\) has valuation
zero. The same argument discards the fictitious second
case for \(N=0\) after capping. In the derivative table
below, terms with \(m>r_0\) at capped \(s\) are also
above all layers used.

For adjacent orders with lower residue less than \(P-1\),
both have the same capped \(s\). At the wrap from \(P-1\) to
zero, the possible capped pairs of valuations of the high
parts are \((0,j)\) or \((j,0)\), \(1\le j\le2H\),
including the case of a zero high part. Use \(H=3\) for
\(1\le\beta<7/4\), and \(H=4\) for
\(7/4\le\beta<15/8\); in both bands \(\gamma P>2\).

A common zero needs at least two minimizing indices in
each order, since a unique minimum gives nonvanishing.
Thus its \(\gamma\) is an intersection of two distinct
lines in (4) for each order. There are finitely many
rational candidates. At a candidate, if
\(\beta=a/q\) in lowest terms, adjoin \(\pi\) with
\(\pi^q=2\), put \(Z=\pi^a/\tau\) and multiply the
polynomial by a power of \(\pi\) to make its coefficient
valuation minimum zero. Each leading coefficient has
residue one over the residue field because the \(U_m\)
and the stripped integer factorial units are odd.
Any common zero must therefore give a common nonzero
root of the two residual polynomials with monomials
at the minimizing indices. Removing initial powers of
\(Z\), the coprimality test can be done in \(\mathbb F_2[Z]\).

The complete list of failing initial residue tests is as
follows. The order residue is that of the upper order.

| \(\beta\) | \(P\) | Upper order residues |
|---|---:|---|
| \(3/2\) | 8 | \(7\) |
| \(7/4\) | 16 | \(0,11,12,14,15\) |
| \(11/6\) | 16 | \(3,4\) |

The attached replay only performs the finite comparison
in (4) and checks the derivative layers used next. It
enumerates the capped adjacent classes, forms all line
intersections in each band, and uses the binary polynomial
Euclidean algorithm on the two leading sets. Thus the
range check is over ordinary starts of every size and
valued extensions of arbitrary ramification, since
possible slopes are obtained from intersections of lines
rather than a bounded ramification grid.

## 3. Separation at the three surviving slopes

For an upper order \(n\ge2\), direct comparison of
falling factorials gives
\[
                        S_{n-1}(z)=S_n(z)-\frac{\theta S_n(z)}n .
\]
A common zero is therefore a zero of both \(S_n\) and
\((\theta-g)S_n\), with any integer \(g\).
Choose \(g\) to be the first minimizing index of \(S_n\).
Let \(\widetilde S\) and \(\widetilde D\) denote the
respective polynomials after substituting \(z=\pi^{-a} Z\)
and normalizing each minimum coefficient valuation to zero.
The coefficient at its layer \(j\) has the form
\(\pi^j u\), where \(u\) is congruent to an odd rational
integer modulo \(4\). A layer table lists exponents with
the indicated coefficient valuation \(j/q\).
For the derivative apply (4) to the moment weights and
add \(v_2(m-g)\). Equations (2)–(3) put the higher blocks
above the displayed layers. All the relevant layers are
in the following table, whose contents are checked by
the replay. An empty set means the layer has no term.

| \(\beta\), upper residue | \(\widetilde S\) layers | \(\widetilde D\) layers |
|---|---|---|
| \(3/2,7\) | \(0:(3,7),\ 1:\varnothing\) | \(0:(2,6),\ 1:(1,5,7)\) |
| \(7/4,0\) | \(0:(7,15)\) | \(0:(0)\) |
| \(7/4,11\) | \(0:(3,7,11),\ 1:\varnothing\) | \(0:(2,6,10),\ 1:(7)\) |
| \(7/4,12\) | \(0:(3,7,11)\) | \(0:(4,12)\) |
| \(7/4,14\) | \(0:(6,14),\ 1,2:\varnothing\) | \(0:(5,13),\ 1:\varnothing,\ 2:(3,7,11)\) |
| \(7/4,15\) | \(0:(7,15),\ 1,2,3:\varnothing,\ 4:(3,11),\ 5:\varnothing\) | \(0:(6,14),\ 1,2:\varnothing,\ 3:(5,13),\ 4:(2,10),\ 5:(3,11,15)\) |
| \(11/6,3\) | \(0:(3,15),\ 1:\varnothing\) | \(0:(2,14),\ 1:(15)\) |
| \(11/6,4\) | \(0:(3,15)\) | \(0:(4)\) |

For the upper residues zero and four in the rows with
a single derivative minimum, nonvanishing of the derivative
follows immediately. At \(7/4\), upper residue twelve,
the derivative residue vanishes only at \(Z=1\) among
nonzero \(Z\) in characteristic two, whereas the leading
moment residue \(Z^3(1+Z^4+Z^8)\) is nonzero there.

In each of the other rows, \(Z\widetilde D\) and
\(\widetilde S\) have the same layer-zero exponents.
Their differences at those exponents have valuation at least
one, because both unit coefficients are odd modulo \(2\).
For \(3/2\), upper residue seven, a common root must reduce
to \(Z=1\), and the three layer-one monomials in
\(Z\widetilde D-\widetilde S\) have sum one there.
For \(7/4\), upper residue eleven, and \(11/6\), upper
residue three, the layer-one difference is a single monomial.
For \(7/4\), upper residue fourteen, the three layer-two
terms of the difference also sum to one at the only possible
residue \(Z=1\). In each case the nonzero layer occurs below
valuation one and forces nonvanishing of the difference.

It remains to check \(\beta=7/4\), upper residue fifteen.
If \(\widetilde S(Z)=0\) and \(Z\) is a unit, then
\(1+Z^8\) belongs to \(2\mathcal O_K\): the moment has
no terms of layers one through three and both leading units
are odd modulo \(2\). At such a \(Z\), every pair of odd
unit coefficients attached to exponents eight apart has
sum of valuation at least one. In
\(Z\widetilde D-\widetilde S\), the layer-three contribution
is such a pair, as are both contributions from the listed
layer-four terms. These therefore have valuations at least
\(7/4\) and \(2\), respectively. The differences of the
layer-zero units at exponents seven and fifteen give
\(2 Z^7(u+u' Z^8)\), with \(u,u'\) each congruent to a
rational integer modulo two by the hypothesis on \(U_m\).
Their combined valuation is either exactly one or at
least two. If it is one, it dominates the difference.
Otherwise the layer-five terms from \(Z\widetilde D\)
dominate. Those three terms reduce to one at the only
possible residue \(Z=1\), so their sum has valuation
exactly \(5/4\). This rules out a common zero in the last
case and proves (1).


## 4. The next four slopes

For \(15/8\le\beta<31/16\) take \(H=5\) in (2) and (4).
The replay verifies that common initial residues occur only
at the following upper-order residues modulo \(32\).

| \(\beta\) | Upper order residues |
|---|---|
| \(15/8\) | \(0,1,22,23,24,25,27,28,29,30,31\) |
| \(53/28\) | \(3,4\) |
| \(23/12\) | \(6,7,8,9\) |
| \(27/14\) | \(3,4\) |

Use the same notation for normalized moments and derivative
gauges as in Section 3. For \(\beta=53/28\) and \(27/14\),
upper residue four has a unique derivative minimum.
At upper residue three, subtract \(\widetilde S\) from
\(Z\widetilde D\) to align layers zero at exponents \(3,31\).
The first higher term of the difference is a single monomial
at exponent \(32\), in layer three for \(53/28\) and layer
one for \(27/14\). Both layers are below valuation one,
so their minima rule out a common zero.

For \(\beta=23/12\), upper residue nine has a unique
derivative minimum. At upper residues six and eight the
difference \(Z\widetilde D-\widetilde S\) has aligned layers
zero and a single monomial in layer two. At upper residue
seven, the failing classes have odd high part in (4), so
the upper order has residue \(39\) modulo \(64\). Use
\(H=6\) for the layer check. The two normalized polynomials
\(\widetilde S,Z\widetilde D\) have the same exponents
\((7,31)\) at layer zero, \((23)\) at layer four, and
\((15,39)\) at layer eight. There are no other terms
through layer thirteen except the derivative terms
\((6,30)\) at layer eleven and \((32)\) at layer thirteen.
At a unit zero of \(\widetilde S\), one has
\(v(1+Z^{24})\ge1/3\) by the absence of layers one through
three. Hence the derivative pair at layer eleven has
valuation at least \(11/12+1/3\). In the difference of
the two polynomials, the aligned layer-zero coefficients
give \(2Z^7(a+bZ^{24})\) with \(a,b\) rational integers
modulo \(2\). This has valuation either one or at least
\(1+1/3\). The aligned layers four and eight cancel to
valuation at least \(1+4/12\) because their corresponding
coefficients are odd. So the layer-thirteen singleton
either follows a dominating layer-zero difference or is
itself the unique contribution of lowest valuation.

For \(\beta=15/8\), the derivative has a unique minimum
at upper residues one and, when the high part is even,
zero. At upper residue twenty-five the derivative
residue has nonzero roots only at \(Z=1\); the moment
residue is \(Z^7(1+Z^8+Z^{16})\), which is nonzero there.
For the remaining rows, the replay verifies the
following comparisons. The shifts shown are operations
on the normalized derivative before comparison.

| Upper residue | Comparison |
|---|---|
| \(0\), odd high part | \(Z\widetilde D-\widetilde S\) aligns layer zero and leaves a singleton in layer two. |
| \(22\) | The same difference aligns layer zero and leaves a singleton in layer two. |
| \(24\) | The same difference leaves the layer-two pair \((9,25)\); its residue is nonzero at a root of \(1+Z^8+Z^{16}\). |
| \(27\) | The same difference aligns layers zero and four, leaving the layer-five triple \((8,16,24)\); its residue is nonzero at \(Z=1\). |
| \(28\) | \(Z^{-1}\widetilde D\) has leading pair \((11,27)\) and no terms in layers one through five. |
| \(29\) | \(Z^2\widetilde D\) aligns layer zero of \(\widetilde S\), with pair \((14,30)\) in layer one and triple \((9,17,25)\) in layer four. |

In the aligned differences in the first four rows the
coefficient differences at aligned exponents have valuation
at least one. For residue twenty-eight, a derivative zero
forces \(v(1+Z^{16})\ge6/8\). The leading and layer-one
pairs in the moment then have valuation at least \(6/8\),
so its layer-four triple \((7,15,23)\) dominates and is
nonzero at residue \(Z=1\). For residue twenty-nine the
same bound on \(1+Z^{16}\) follows from the moment,
which has no terms in layers one through five.
The shifted derivative pair at layer one is pushed to
valuation at least \(7/8\), and its layer-four triple
gives nonvanishing.

The remaining residues at \(15/8\) are twenty-three,
thirty and thirty-one.

- At twenty-three the common leading exponents of
  \(\widetilde S,Z\widetilde D\) are \(7,15,23\).
  Writing \(P(Z)=1+Z^8+Z^{16}\), a moment zero forces
  \(P(Z)\in2\mathcal O_K\) since it has no terms in layers
  one through seven. The difference of the layer-zero
  coefficients is \(2 Z^7\) times a polynomial in \(Z^8\)
  of degree at most two with rational mod-two coefficients.
  If this polynomial has nonzero residue at the root,
  the difference has valuation one. Otherwise it is a
  multiple of \(1+Y+Y^2\) modulo \(2\) in its variable
  \(Y\), so its evaluation has valuation at least one
  and the layer-zero difference has valuation at least two.
  The layer-seven triple of the shifted derivative
  evaluates in \(Z^6(2\mathcal O_K)\), hence contributes
  at valuation at least \(15/8\).
  The possible layer-eight term has the same exponent on
  both sides and cancels to valuation at least two.
  The layer-nine singleton of the shifted derivative thus
  cannot cancel.
- At thirty the common leading exponents are \(14,30\).
  The moment has one pair at layer seven, separated by
  sixteen, and no earlier higher layers. At its zero,
  \(1+Z^{16}\in2\mathcal O_K\). To see this, reduce the
  equation modulo \(2\); the moment is \(1+Z^{16}\)
  times \(Z^{14}+\pi^7 Z^{13}\), a unit. The layer-zero
  difference of \(Z\widetilde D-\widetilde S\) is then
  either of valuation one or at least two, using rational
  mod-two coefficients after division by two.
  Layers six, seven and eight consist of pairs separated
  by sixteen with odd coefficients and gain at least one
  in valuation. The layer-ten triple \((8,16,24)\) of the
  shifted derivative is nonzero at residue \(Z=1\).
  Its valuation \(10/8\) cannot cancel.
- At thirty-one, use \(H=6\) to check the additional
  layers; the lifted residues are \(31,63\).
  The common leading exponents are \(15,31\).
  Every nonempty moment layer from eight through sixteen
  is a pair separated by sixteen, with no layers one
  through seven. At a moment zero \(1+Z^{16}\in2\mathcal O_K\).
  All higher moment layers then have valuation at least
  two, so \(Z^{16}\) is congruent to an odd rational
  integer modulo \(4\), using the leading-unit quotient.
  It follows that \(v(Z^8-1)\ge1/2\), by factoring
  \(Z^{16}-1\).
  In \(Z\widetilde D-\widetilde S\), the layer-seven pair
  either evaluates at valuation \(15/8\) or has valuation
  at least \(23/8\): its coefficients are odd rational
  modulo \(4\) and \(Z^{16}\) is rational modulo \(4\).
  The layer-zero difference is \(2Z^{15}\) times a linear
  polynomial in \(Z^{16}\), whose coefficients are
  rational modulo \(4\) by the mod-eight unit hypothesis.
  It has valuation one, two, or at least three.
  The aligned layer-eight difference is \(4Z^7\) times a
  linear polynomial in \(Z^{16}\) with rational
  coefficients modulo \(2\), so it has valuation two
  or at least three. If both contributions have
  valuation two, their sum has valuation at least
  \(5/2\), since after factoring \(4Z^7\) they are
  congruent to \(1+Z^8\) modulo \(2\).
  All other layers below seventeen are empty or consist
  of pairs separated by sixteen at layer twelve or
  higher, hence have valuation at least \(5/2\).
  Therefore either one of the earlier nonzero valuations
  \(1,15/8,2\) dominates, or the layer-seventeen triple
  \((8,24,32)\) dominates. The latter reduces to one at
  residue \(Z=1\).

This completes the separation up to \(31/16\). The proof
keeps the actual two-term window at arbitrarily large
prime-power exponents. Above the proved shallow range the
number of candidate digits grows as the slope approaches
two. Uniform control of subsequent common residues, or a
different argument for the remaining interval in (1), is
still needed over arbitrary ramified coefficient fields.
