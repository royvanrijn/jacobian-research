# Prime-block transfer for composite even-quadratic windows

Let \(Q=ax^2+bxy+cy^2\), \(f=Q+k\), and \(M_j=\Lambda_2(f^j)\),
where \(\Lambda_2(x^iy^j)=i!j!\). Put \(t=b/2\) and \(\Delta=ac-t^2\).
For \(m\ge2\), define
\[
 R_m(X)=\sum_{\ell=0}^m\binom m\ell(\ell!)^2X^{m-\ell},
 \qquad E_m(X,Z)=\sum_{\ell=0}^m\binom m\ell(2\ell)!Z^\ell X^{m-\ell}.
\]
Take a prime \(p>2m+6\) with \(p\nmid\operatorname{Res}(R_{m-1},R_m)\).
This integer is nonzero by the
[constant-plus-product two-window](CONSTANT_PLUS_PRODUCT_OF_TWO_VARIABLES_ALL_ORDERS.md).

**Theorem.** Suppose \(abck\ne0\) and four consecutive moments
vanish in a block containing both \(mp,mp+1\).
At a place over \(p\), take \(v(p)=1\) and let bars denote residue.
After exchange and scaling, one can arrange \(a=1\) and
\[
 v(k)=2/p,\quad E_m(\bar K,1)=0,\quad K=k^p/p^2.                \tag{1}
\]
The residual profiles satisfy both conditions
\[
\begin{split}
 &\bar c\in\{0,1\}\ \text{or}\ \bar t\in\{1,\bar c\}\
                        \text{or}\ E_m(\bar K,\bar c^p)=0,\\
 &\bar\Delta=0\ \text{or}\ \Lambda_2((\bar Q^{[p]}+\bar K)^{m-1})=0,
 \qquad Q^{[p]}=a^p x^2+b^p xy+c^p y^2.                       \tag{2}
\end{split}
\]
There are at most \(m(m+4)(m+1)\) triples \((\bar c,\bar t,\bar K)\)
in this normalization. When \(\bar\Delta\ne0\), the additional
order-\(m\) conditions (5) below hold. If \(mp\) is penultimate
and \(\bar c=0\), then also \(\bar t=0\). The triples in
Frobenius coordinates \((\bar c^p,\bar t^p,\bar K)\) lift
to a fixed list of characteristic-zero templates depending
only on \(m\). For fixed \(m\) the theorem
applies at infinitely many composite windows. Further lifting is
needed to decide if the residue profiles survive. If instead
\(mp\) is the last index, the two charts specified below
are still excluded.

Over an arbitrary characteristic-zero field, use the weak
Nullstellensatz with the nonzero conditions to realize the
equations over \(\overline{\mathbb Q}\), then take a place over \(p\).

## The homogeneous blocks

Use integral coefficients throughout the local calculations.
Congruence modulo \(p\) means valuation at least one,
while the bars denote reduction modulo the maximal ideal.
Put \(h=(p-1)/2\) and \(H_i=\Lambda_2(Q^i)\).
For \(i=\ell p+r\le mp+3\), \(0\le r<p\),
\(v(H_i)\ge2\ell\) because \(v(i!)=\ell\) and
[the primitive form is integral](BINARY_QUADRATIC_PRIMITIVE_HOMOGENEOUS_RESULTANT.md#primitivity).
For \(r\ge h\) the coefficientwise congruence is
\[
             p^{-2\ell}H_{\ell p+r}\equiv
                  H_\ell(Q^{[p]})H_r\pmod p .                \tag{3}
\]
A monomial of total degree \(2\ell p+2r\) has its exponent
quotients on division by \(p\) sum to at least \(2\ell\),
since the two remainders sum to at most \(2p-2\).
Terms with larger quotient sum vanish after division by
\(p^{2\ell}\). For the surviving terms use
\(Q^{\ell p+r}\equiv(Q^{[p]}(x^p,y^p))^\ell Q^r\bmod p\).
Monomials of \(Q^r\) with an exponent at least \(p\)
increase the quotient sum and vanish. The quotient unit
\((up+v)!/p^u\equiv(-1)^u u!v!\pmod p\) for \(up+v<p^2\),
\(0\le v<p\), follows by grouping complete \(p\)-blocks.
The signs cancel at total quotient \(2\ell\).

For \(0\le r<h\), put
\(L=(Q(x,y)^p-Q^{[p]}(x^p,y^p))/p\) and define \(J_r^x,J_r^y\)
by the respective sums
\[
 \sum_{\substack{0\le u,v<p\\u+v=p+2r}}u!v!
     [x^{p+u}y^v]LQ^r,\qquad
 \sum_{\substack{0\le u,v<p\\u+v=p+2r}}u!v!
     [x^uy^{p+v}]LQ^r .
\]
The quotient sum can now be \(2\ell-1\). There the term
with no \(L\) is absent since each exponent of \(Q^r\)
is below \(p\). The first \(L\)-term is
\(\ell p(Q^{[p]}(x^p,y^p))^{\ell-1}LQ^r\), and higher
\(L\)-terms vanish. The quotient sum from \(LQ^r\) must be one,
giving the two sums above. The total quotient sum is odd,
giving the sign in
\[
 p^{-2\ell}H_{\ell p+r}\equiv H_\ell(Q^{[p]})H_r
  -\ell\Lambda_2((J_r^x x+J_r^y y)(Q^{[p]})^{\ell-1})
                                                          \pmod p.
\]
At \(\ell=0\) omit the correction. The bound \(2(mp+3)<p^2\)
ensures the factorial unit formula applies throughout.

## The interior reduction at the transition

Suppose \(v(k)=2/p\). Set \(F=\Lambda_2((Q^{[p]}+K)^{m-1})\)
and for \(0\le d\le3\) put
\[
 \begin{split}
 U_d&=\frac{m}{p}\sum_{r=h+1}^{p-1}
       \frac{-d!}{r!(p+d-r)!}\,k^{p+d-r}H_r,\\
 G_d&=\sum_{\ell=0}^m\binom m\ell K^{m-\ell}
                         \frac{H_{\ell p+d}}{p^{2\ell}} .
 \end{split}
\]
Writing \(o_v(0)\) for valuation strictly above zero gives
\[
                     M_{mp+d}/p^{2m}=F U_d+G_d+o_v(0).      \tag{4}
\]
For \(i=\ell p+r\), the valuation of the binomial coefficient
is one when \(r>d\) and zero otherwise, by the factorial
valuation formula. When \(r>d\), dividing the binomial by \(p\)
gives the residue \(-m\binom{m-1}\ell d!/[r!(p+d-r)!]\)
by the factorial unit formula. For \(r=d\) the residue is
\(\binom m\ell\). After division of the moment by \(p^{2m}\),
the classes \(r<d\) and \(d<r\le h\) have strictly positive
valuation. The class \(r=d\) gives \(G_d\).
For \(r>h\) apply (3) and the binomial residue. The discrepancies
gain at least one, and the least possible valuation before
the gain is \(-1+2(d+1)/p>-1\). The terms that remain sum to \(FU_d\).

When \(\bar\Delta\ne0\), \(H_{p-1}\) is a unit by the
[central coefficient calculation](BINARY_EVEN_QUADRATIC_PRIME_TERMINAL_FOUR_WINDOW.md#reduction-at-the-terminal-prime).
Only \(r=p-1\) attains \(v(U_d)=-1+2(d+1)/p\).
As \(G_0\) is integral, \(M_{mp}=0\) gives \(v(F)\ge1-2/p\).
For \(d>0\) this puts \(FU_d\) strictly above valuation zero.
Thus another zero \(M_{mp+d}=0\) implies
\[
 0=\overline{G_d}
  =\overline{H_d\Lambda_2((Q^{[p]}+K)^m)
    -m\Lambda_2((J_d^x x+J_d^y y)(Q^{[p]}+K)^{m-1})}.       \tag{5}
\]
The last equality sums the low-remainder congruence with
the binomial weights.

## Two excluded integral residue charts

Make all coefficients integral with some unit.
At \(mp\), Frobenius gives \(M_{mp}\equiv k^{mp}\)
in the residue field, so \(s=v(k)>0\). Consider either
two nonunit axes, or a unit axis scaled to \(a=1\)
with \(\bar c=0,\ \bar b=2\). The primitive form
\(H_i/(i!)^2\) specializes to \(\gamma^i\) at the
respective reductions, with \(\gamma=b\) or \(4\).
The first identity is the moment of \(bxy\).
The second follows from the
[single-face recurrence](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#2-the-planes-b2a-and-b2c)
with zero constant and zero opposite axis. Hence
\(v(H_{\ell p+r})=2\ell\) throughout both charts.

If \(s<2/p\) the least term in each \(r\) class uses
\(\ell=0\). At \(mp\) only \(i=0,p-1\) can be least.
Their valuations are \(mp s\) and \(1+(mp-p+1)s\),
so a zero forces \(s=1/(p-1)\). At \(mp+1\) and this
strength \(i=1\) is uniquely least, at \(mp s\).
Every class \(r>1\) has valuation at least \(mp s+s\).
If \(s>2/p\) each \(r\) class instead uses its largest
\(\ell\). For \(2/p<s<1\), the \(i=mp-1\) term of
\(M_{mp}\) is uniquely least, at \(2m-1+s\). For
\(s\ge1\), the top term of \(M_{mp+1}\) is uniquely
least, at \(2m\), since the classes \(r>1\) are at least
\(2m-1+2s\). Thus the strength must be \(s=2/p\).

Here \(\bar\Delta\ne0\) in both charts, so (4) gives
\(\bar F=\overline{G_1}=0\). But the factorial unit
formula and the two specializations give
\[
      \overline{p^{-2\ell}H_{\ell p+r}}
           =(\ell!r!)^2\bar\gamma^{\ell p+r}.
\]
Consequently \(\bar K/\bar\gamma^p\) would be a common
root of \(R_{m-1}\) and \(R_m\) in the residue field.
The hypothesis on \(p\) excludes both charts.

The same two charts are excluded if the block ends at \(mp\).
The moment at \(mp-1\) is available in place of \(M_{mp+1}\).
In \(M_{mp-1}\) all binomial coefficients are \(p\)-units.
For \(s<2/p\) the term at \(i=p-1\) is uniquely least by
first minimizing over \(\ell\) and then maximizing \(r\).
For \(s>2/p\) the same argument gives the unique minimum
at \(i=mp-1\). At \(s=2/p\), (3) and the unit binomial
formula give
\[
   p^{-2(m-1)}M_{mp-1}
       = F\sum_{r=h}^{p-1}\binom{p-1}r H_r k^{p-1-r}
                       + O_v(1),
\]
where \(O_v(1)\) has valuation at least one. The terms at
\(r<h\) have valuation at least \(2(p-h)/p>1\), while
every discrepancy in the remaining congruences gains at
least one. The sum is a unit in these charts by its
\(r=p-1\) term. Therefore \(v(F)\ge1\), so (4) now gives
\(\overline{G_0}=0\). As above this would give a common
root of \(R_{m-1},R_m\), using \(r=0\) in the chart
specialization.

## Weighted faces

Normalize \(a=1\) by a unit axis and put
\(B_j(Z)=\Lambda_1((Zx^2+k)^j)\). By the
[exact coupled recurrence](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation),
the weighted values \((c-t)B_j(c)+(1-t)B_j(1)\) vanish
for the last two indices, of the form \(mp+d\) with
\(0\le d\le3\). At least one of the two has \(d>0\). Equal axes and exactly zero weights
are covered by the
[single-face theorem](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md).

At unit \(Z\), the \(i=\ell p+r\) term of \(B_{mp+d}\) has valuation
\[
 2\ell+\lfloor2r/p\rfloor+\mathbf 1_{r>d}
                       +((m-\ell)p+d-r)s .
\]
The unique minimum is at \(i=d\) for \(s<2/p\) and
at \(i=mp+d\) for \(s>2/p\). Among \(r\le d\), first
minimize over \(\ell\) and then maximize \(r\).
For \(r>d\) the optimal \(\ell\) is respectively
\(0\) or \(m-1\). Its excess valuation over the claimed
minimum is respectively \(1+\lfloor2r/p\rfloor-(r-d)s\)
or \(-1+\lfloor2r/p\rfloor+(p+d-r)s\), both positive
in the corresponding regime. At equality only \(r=d\)
attains valuation \(2m\).

For nonunit \(c\), use one boundary index with \(d>0\).
The weight \(1-t\) is a unit by the chart exclusions. In either strict regime the positive-index
minimum at \(Z=1\) lies strictly below all terms at \(Z=c\).
For unit \(c\), the quotient of the adjacent boundary
polynomials, after stripping a scalar, reduces to \(Z\)
on the unit arguments. It separates distinct residues.
On a residue disk its difference quotient is a unit,
since corrections to \(Z\) have coefficients of positive
valuation and denominators that are units at both arguments.
The weighted cancellations would equate the quotient at
\(1,c\). Thus \(s=2/p\) is necessary.

At \(i=\ell p+d\) the factorial unit formula now yields
\[
        \overline{p^{-2m}B_{mp+d}(Z)}
                     =(2d)!Z^d E_m(\bar K,Z^p).
\]
If \(\bar c=0\), use an index with \(d>0\) to obtain
\(E_m(\bar K,1)=0\) from the unit axis. If \(d=0\)
is also a boundary index, its reduced equation now
reads \(-\bar t\bar K^m=0\), giving \(\bar t=0\). For distinct unit residues of the
axes the two adjacent equations give a nonsingular
system for the two weighted contributions. Hence every
axis with unit weight gives a root of \(E_m(\bar K,Z^p)\),
and at least one weight is a unit. If \(\bar c=1\) and
\(\bar t\ne1\), the weight sum is a unit and again
\(E_m(\bar K,1)=0\). If \(\bar c=\bar t=1\) and
\(E_m(\bar K,1)\ne0\), the common residual factor
cancels in the quotient on that disk, leaving a unit
scalar times \(Z\). The same injectivity contradicts
equality at distinct arguments.

Exchange and scale when needed to put the axis giving
a root at \(a=1\). This proves (1), and the weighted
system gives the first line of (2). Formula (4) gives
the second line. The monic polynomial \(E_m(X,1)\)
has at most \(m\) roots. For a fixed \(\bar K\), either
\(\bar c\) belongs to the set of at most \(m+2\) fixed
values from the first line of (2), or \(\bar t=1\) or
\(\bar t=\bar c\). Frobenius is injective on the residue
field. For fixed \(\bar c\), the rank-one equation gives
at most two \(\bar t\), while the lower moment has
degree \(m-1\) in \(\bar t^p\) with unit leading coefficient
from the pure \(xy\) term. On each of the planes
\(\bar t=1,\bar t=\bar c\) the same count gives at most
\(m+1\) values of \(\bar c\). The leading quadratic
directions there are \(y^2\) and \(y^2+2xy\), respectively,
whose moments at order \(m-1\) are units. This gives
the stated bound.

## A finite list over characteristic zero

Use variables \(C,T,X\) for the three Frobenius coordinates.
Replace \(\bar c^p,\bar t^p,\bar K\) by \(C,T,X\) in
(1)–(2), in particular replacing the second line of (2) by
\[
 C-T^2=0\quad\text{or}\quad
 \Lambda_2((x^2+2Txy+Cy^2+X)^{m-1})=0.
\]
These are fixed integer equations. Each branch first solves
the monic equation \(E_m(X,1)=0\). Next either \(C\)
is chosen from \(0,1\) and the roots of \(E_m(X,C)\)
or one substitutes \(T=1\) or \(T=C\).
The last equation, of degree at most \(m-1\) or two,
has the unit leading coefficients used in the count.
After inverting \((2m)!\), all the branch equations
can thus be made monic in the successive unknowns.
The characteristic-zero lists in these branches have
total size at most \(m(m+4)(m+1)\).

Every residual candidate is obtained by reducing a member
of this list at a place over \(p\). Indeed, a monic
polynomial with integral coefficients splits into integral
roots over \(\overline{\mathbb Q}_p\), and reducing that
factorization gives all its residual roots. Apply this
first to \(X\), then to each successive monic equation
with the earlier lifted values. The lifted values are
algebraic over \(\mathbb Q\) because the equations and
earlier values are algebraic. This realizes the claimed
finite list independently of \(p\). It is a lift of the
leading equations only.

## Ramification of the coefficient field

Take \(L=\mathbb Q(a/k,b/k,c/k)\) for an algebraic solution.
At every place of \(L\) over an eligible \(p\) with \(mp,mp+1\)
in the block, the renormalized constant in (1) is
either \(k/a\) or \(k/c\), hence lies in \(L\).
Its valuation \(2/p\) implies that \(p\) divides the
ramification index at that place. Consequently
\(p\mid[L:\mathbb Q]\), by the sum of ramification indices
times residue degrees. This holds simultaneously for
each eligible prime from any of the first three indices.
In particular a field of degree less than \(p\), or any
number field with an unramified place above \(p\), already
satisfies the four-window for these indices.

## A positive density of full windows over fields of bounded degree

Let \(\mathcal A\) be the positive starts \(n\) at which all
even four-term quadratics with
\([\mathbb Q(a/k,b/k,c/k):\mathbb Q]\le n^{2/3}\)
satisfy the four-window. Then
\[
       \liminf_{Y\to\infty}
         \frac{\#(\mathcal A\cap[1,Y])}{Y}\ \ge\ \log(3/2).
\]
In particular, this gives the same density bound for any
fixed number field. It uses the start itself as the index
\(mp\), leaving its next index available in the block.

Fix \(2/3<\alpha<1\). An integer \(n\le Y\) has at most
one prime divisor \(p>Y^\alpha\). The number with such a
divisor is therefore
\[
                \sum_{Y^\alpha<p\le Y}\lfloor Y/p\rfloor
                       =Y\log(1/\alpha)+o(Y)
\]
by [Mertens' prime harmonic formula](https://www.degruyterbrill.com/document/doi/10.1515/crll.1874.78.46/html).
That formula also gives
\(\pi(Y)=o(Y)\), by bounding the primes above \(Y^\theta\)
by \(Y\) times their reciprocal sum and letting the fixed
\(\theta<1\) increase to one. Thus discarding \(m=1\)
costs only \(o(Y)\).

For \(n=mp\) in this set with \(m\ge2\), one has
\(m\le Y^{1-\alpha}\), and \(p>2m+6\) for large \(Y\).
The exceptional primes have negligible total contribution.
Indeed \(\binom mj(j!)^2\le m^{2j}\), so the monic
\(R_m\) has all roots at radius less than \(2m^2\) by
leading term dominance. On that disk
\(\lvert R_{m-1}\rvert\le 2(2m^2)^{m-1}\), by the
same coefficient bound and the geometric series.
The monic resultant formula gives
\(\log|\operatorname{Res}(R_{m-1},R_m)|=O(m^2\log(m+1))\).
The number of distinct exceptional \(p>Y^\alpha\) at
a given \(m\) is at most this logarithm divided by
\(\alpha\log Y\). Summing to \(Y^{1-\alpha}\) gives
\(O(Y^{3(1-\alpha)})=o(Y)\).

For every remaining \(n\) a putative counterexample has
coefficient field degree divisible by \(p>Y^\alpha>n^{2/3}\).
Hence the lower density is at least \(\log(1/\alpha)\).
Letting \(\alpha\) decrease to \(2/3\) proves the bound.

For the prime harmonic estimate, write
\(\Lambda_{\mathbb Z}(p^j)=\log p\) for the von Mangoldt
function. Each dyadic prime block has log sum \(O(x)\)
from its central binomial coefficient, so
\(\sum_{q\le x}\Lambda_{\mathbb Z}(q)=O(x)\), the higher
prime powers contributing \(O(\sqrt{x}\log x)\).
The identity
\(\sum_{q\le x}\Lambda_{\mathbb Z}(q)\lfloor x/q\rfloor
=\log(\lfloor x\rfloor!)\) then yields
\(\sum_{q\le x}\Lambda_{\mathbb Z}(q)/q=\log x+O(1)\).
The higher prime powers in this last sum are absolutely
summable. Partial summation of the prime terms against
\(1/\log t\) gives
\(\sum_{p\le x}1/p=\log\log x+B+O(1/\log x)\):
the bounded error in \(\sum_{p\le t}(\log p)/p-\log t\)
is integrable against \(1/(t\log^2 t)\).

The [exact replay](verify_binary_even_quadratic_composite_prime_block_transfer.py)
checks the block identities and the boundary reduction
with exact factorials and ramified valuations.
