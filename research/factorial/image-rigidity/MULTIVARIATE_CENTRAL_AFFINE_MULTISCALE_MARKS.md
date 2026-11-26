# Central factorial marks for independently repeated reflected rows

Fix an ordinary two-monomial core and finitely many possible appended
unequal row types. In type \(i\) the smaller exponent is
\[
                  \ell_i(T)=c_iT^{e_i}+a_i
\]
with \(c_i,e_i\) positive integers and \(a_i\) a nonnegative integer.
The gap \(d_i\) is a fixed positive integer. Either orientation is
allowed, and the nonnegative integer multiplicities are independent.
Call a row type present if its multiplicity is positive. Call a
present type paired if its reflection, with the same smaller exponent
function and gap, is also present. The two multiplicities may differ.
A degree block consists of the types with the same \(e_i\).
Assume that distinct degrees in the template have distinct
\(2\)-adic valuations. This last condition can be replaced by a
finite coprimality test between the literal ramp polynomials of
different degree blocks at the relevant pairs of orders.

Use a computable
[absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
\(J\ge5\). There is a computable \(T_0\), depending on the fixed core
and the template alone, with the following property. At each residual
adjacent pair \(n,N=n+1\) with \(5\le N<J\), either of the following
conditions proves physical coprimality for \(T\ge T_0\).

- Some occupied nonlinear block has a paired type maximizing
  \[
         (m a_i+m d_i/2)/c_i
  \]
  among the present types in its maximal \(v_p(c_i)\) tier. Here
  \(m\) is the even order of the pair and \(p\) is any prime dividing
  \(m/2\). In particular, an occupied nonlinear block
  in which all present types are paired suffices at every pair.
- In the occupied affine block, all present slopes are prime to
  \(nN\). At each of the two orders \(j=n,N\) a paired type maximizes
  \[
                       (j a_i+d_i\lceil j/2\rceil)/c_i
  \]
  among the present types of that block. The two maximizing types
  can differ.

If the condition holds at each residual pair, the original moments
satisfy SFC(2) in every degree. The affinely growing case uses two
different moving integer values. The nonlinear case uses one of them
at the even order. Both conclusions hold uniformly over all
independent admissible counts, without a signed balance condition.
For example, one can choose any number of degree blocks at scales
\(T^2,T^4,T^8,\ldots\), with arbitrary slopes and fixed gaps. If one
occupied block contains only reflected pairs in the sense of presence,
the counts of its two orientations and all other counts can vary
independently. The angular regularity of the fixed core is irrelevant.

## Exact central coefficient minima

Fix the orders \(n,N\). Strip the scalar \((j\ell_i(T))!\) in order
\(j\). The remaining factors from type \(i\) come from
\[
             j c_i T^{e_i}+j a_i+h,\qquad 1\le h\le j d_i .
\]
An increasing row contributes the initial segment of length \(d_i k\)
at coefficient \(k\), and a decreasing row contributes the initial
segment of length \(d_i(j-k)\). Take a finite prime set \(S\)
containing \(2\), the prime divisors of the fixed core's stripped
factors and of \(N!\), the prime divisors of the leading coefficients
of all literal factors, and the prime divisors of the nonzero integer
resultants of pairs of nonproportional literal factors. All
nonproportional factors used below have nonzero resultant whenever
one of them is the selected mark. Factors of the same degree block
are proportional exactly when their normalized constant offsets
agree. For different degree blocks, a common root \(t\) would satisfy
\(t^{e_i}=-A\) and \(t^{e_l}=-B\) for \(A,B>0\). Put
\(g=\gcd(e_i,e_l)\). Raising these two equations to powers \(e_l/g\)
and \(e_i/g\) gives values of opposite signs. This proves the
separation of different degree blocks when \(v_2(e_i)\ne v_2(e_l)\).

Suppose \(F(T)=j c_i T^{e_i}+\eta\) is a literal factor and
an outside prime divides \(F(T)\). Only factors proportional to
\(F\) can be nonunits at that prime. In particular all the companion
order's stripped coefficient ordinates are zero when it has no
proportional factor. Suppose also that the marked order's ordinates
vanish precisely on the middle indices
\[
                           \lfloor j/2\rfloor,\ldots,\lceil j/2\rceil .
\]
This gives horizontal Newton polygon width at most one. Indeed a
horizontal lower edge has supporting line of slope zero, and its
endpoint ordinates are therefore global minima. Off the middle
indices all ordinates are positive. The companion contributes only
slope zero. Any ordinary common factor would have degree at least
two by even real positivity, and its slopes, with multiplicities,
must belong to both parents. This contradicts the width bound.
Here the Newton polygon product rule is used in the same form as
in the
[exact affine cone theorem](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md#all-but-one-mark-gives-an-exact-slope-constraint).
[Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)

We now produce the indicated middle ordinates. At order \(j\)
select a paired type \(i\) maximizing
\((j a_i+d_i\lceil j/2\rceil)/c_i\) in the types allowed to share a
mark. Mark an offset
\[
                     \eta=j a_i+d_i\lceil j/2\rceil+s,
                     \qquad 1\le s\le d_i .
\]
The offset stays within that type's literal full ramp for \(j\ge4\).
Every allowed type has both central coefficient arguments strictly
below the mark in normalized offset. Thus no proportional factor
occurs at the central indices. For \(k<\lfloor j/2\rfloor\) the
decreasing row of type \(i\) contributes the mark, since
\(j-k\ge\lceil j/2\rceil+1\). For \(k>\lceil j/2\rceil\) the
increasing row contributes it. This proves exactly the central
minimum assertion, with no lower bound on the two positive counts.

For the affine clause, all slopes in its block are units modulo
\(nN\). Choose \(s\) such that \(j\nmid\eta\). If \(d_i\ge2\)
one of \(s=1,2\) works. If \(d_i=1\), the residue
\(\lceil j/2\rceil+1\) is nonzero modulo \(j\) for \(j\ge4\).
Proportionality to an affine factor at the other order would give
\[
                   j c_i(j' a_l+h)=j'c_l\eta ,
\]
which is impossible modulo \(j\). The central-maximum comparison
used the whole affine block, so every proportional marked-order
factor is permitted. The two orders give distinct linear marks.
By the effective two-variable \(S\)-unit bound as used in the
[exact affine cone theorem](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md#all-but-one-mark-gives-an-exact-slope-constraint),
two distinct positive affine marks cannot both take \(S\)-smooth
values beyond a computable common cutoff. At least one outside
prime proves coprimality.

For the nonlinear clause, write \(j=m\) for the even order.
If \(4\mid m\) one may use \(p=2\), and otherwise any prime
dividing the odd integer \(m/2>1\) works. Each such \(p\) divides
\(m/2\), so \(\eta=m a_i+m d_i/2+1\) is a \(p\)-unit.
Proportionality within the marked order requires
\(c_i(m a_l+h)=c_l\eta\). Since \(i\) has maximal \(v_p(c_i)\)
among the present types in the block, any such \(l\) belongs to
the same maximal tier. This is the set used in the central
comparison. A proportional factor at the odd order \(j'\) would
require \(m c_i(j'a_l+h)=j'c_l\eta\). The left side has valuation
at least \(v_p(m)+v_p(c_i)\) and the right side at most \(v_p(c_i)\);
here \(p\nmid j'\). This excludes all proportional companion
factors. The mark \(F(T)=m c_iT^{e_i}+\eta\) has \(e_i\ge2\)
distinct roots, because \(\eta>0\) and its derivative can vanish
only at zero. The effective \(S\)-part bound used in the
[polynomial native-mark theorem](MULTIVARIATE_POLYNOMIAL_HEIGHT_NATIVE_MARKS.md#exact-prime-isolation-and-the-retained-charge)
supplies an outside prime at every sufficiently large integer
\(T\). These effective bounds can be taken over the finite menu
of orders and literal factors, independently of which maximizing
paired types are present. The
[low-order theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
treats terminals through four and the absolute cutoff handles
orders beyond the finite menu. Rational coprimality gives
coprimality over all characteristic-zero fields.

## Three affine projections for unit gaps

Consider just the affine block, with every gap one, and suppose the
slopes \(c_i\) are all congruent to the same unit modulo \(J!\).
For \(5\le N<J\), let \(o\) be the odd order of \(n,N\).
Write \(b_i\) for the signed count of type \(i\), positive for
increasing and negative for decreasing, and associate the line
\(L_i(x)=(a_i+x)/c_i\) to that type. Beyond a computable cutoff,
a collision would require exact signed conservation on every fiber
of all three projections
\[
                             x=1/n,\qquad 1/N,\qquad 2/o.
\]
In particular every occupied line has an opposite-orientation
occupied line meeting it at each of these three distinct points
unless it already has its reflection present. In that case the
central theorem excludes a collision. Distinct affine lines meet
at most once. It follows that a template with at most two types
in at least one orientation satisfies SFC(2) in every degree above a
computable height. More generally one can check whether the
three finite set partitions admit any nonzero signed vector with
the required orientation signs. This extends the
[first-mark circulation criterion](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md#a-circulation-condition-for-unit-gaps).

For the proof use the marks \(j c_i T+j a_i+h\) with
\((j,h)=(n,1),(N,1),(o,2)\). A proportional factor at the other
order is impossible by the preceding modular comparison.
A same-order factor \(j c_lT+j a_l+h'\) can be proportional
only if \(h'c_i\equiv h c_l\bmod j\). All slopes are equal units
modulo \(j\); since \(1\le h'\le j\), one has \(h'=h\).
The mark classes within an order are exactly the fibers of
\(L_i(h/j)\). The classes across these three projections are
disjoint: across orders the companion comparison applies, and
within the odd order the factor residues differ.

Take a good prime for a mark and divide the coefficient ordinates
by the valuation of its primitive linear polynomial. The constant
multipliers of all literal factors are units at this prime. Write
\(U,V\) for the total increasing and decreasing counts in the fiber.
The marked array is the sum of \(U1_{\{k\ge h\}}\) and
\(V1_{\{k\le j-h\}}\). If \(U\ne V\), it attains its minimum
on the \(h\) indices at one endpoint, because \(j-h\ge h\).
It then has horizontal width \(h-1\le1\), against a flat
companion polygon. The Newton comparison forces \(U=V\).
At most one linear mark across all three families can lack
a good prime past the effective \(S\)-unit threshold. A complete
family gives total signed count zero. If a fiber in another family
lacks a good prime, subtract its other signed sums from this same
total, which forces its signed count to be zero as well.
This proves conservation at all three projections.

There are templates with a directed cycle in the first two
partitions that the third projection excludes. Fix a residual
pair \(n,N\), put \(L=6J!\), and take a rectangular array of
integers
\[
 X_e=X+eL,\qquad Y_f=X-1-2fL,\qquad
 R_{ef}=X_e-Y_f=1+(e+2f)L ,\qquad e,f\in\{0,1\}.
\]
The four \(R_{ef}\) are pairwise coprime: any common prime divides
the difference of two, whose primes divide \(L\), whereas all four
integers are one modulo \(L\). Each \(R_{ef}\equiv1\bmod J!\).
Choose \(X\) by the Chinese remainder theorem so that
\(X_e\equiv N R_{ef}\bmod nN R_{ef}\) for the four cells, and
choose it large enough for all resulting \(a_{ef}\) to be
nonnegative. The congruences agree modulo the common factor
\(nN\), as \(L\) is a multiple of \(nN\). Set
\[
 W=nN\prod_{e,f}R_{ef},\qquad
 c_{ef}=\frac{W}{nN R_{ef}},\qquad
 a_{ef}=\frac{X_e}{nN R_{ef}}-\frac1n .
\]
These are admissible integer affine types with
\(L_{ef}(1/n)=X_e/W\) and \(L_{ef}(1/N)=Y_f/W\).
Taking the diagonal types increasing and the off-diagonal
types decreasing gives a directed four-cycle. Independent counts
on this template nonetheless satisfy the all-degree result above
because there are only two types in each orientation and no
reflected pair of identical affine functions.

A convenient source for the two smoothness inputs is
[Bérczes, Evertse and Győry, Theorem 2.1](https://math.leidenuniv.nl/reports/files/2008-16.pdf)
for effective \(S\)-unit finiteness, and
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/pdf/1708.08290)
for effective bounds on the \(S\)-part of a polynomial value
when the polynomial has at least two distinct complex roots.
The latter applies to the nonlinear marks themselves, even
when they are reducible over the rationals.

## Varying row types in a growing height band

There is a computable constant \(\kappa_J>0\) with an additional
template-free consequence. Fix any ordinary core and append
any finite multiset of unequal shared rows with both exponents in
\([H,H+W]\), where \(H,W\) are positive integers and
\[
           W\le\kappa_J
                 \frac{(\log_2 H)(\log_3 H)}{\log_4 H}.
\]
Use the iterated logarithms truncated below by one. There is no bound
on the number of row types or rows, no relation between their heights,
and no equal-multiplicity requirement. Suppose that, at each residual
pair \(n,N\) and each of its two orders \(j\), at least one reflected
pair is present whose larger central factorial argument
\[
                    j\min(a,b)+|b-a|\lceil j/2\rceil
\]
is maximal among the appended types. In particular this holds if
each present type is paired. Then the original ordinary moments
satisfy SFC(2) in every degree once \(H\) exceeds a computable
threshold, uniformly over all these choices.

For a pair \(n,N\), set \(B=\max(3,J,JM_C,J^2 W)\), with \(M_C\)
the largest core exponent or zero for the empty core. At each order
\(j=n,N\) take a maximizing pair and a mark \(x_j\) between one and
its gap above its larger central argument, with \(j\nmid x_j\).
The same residue choice used for affine slopes works. Each mark
lies between \(jH+1\) and \(j(H+W)\). If \(p>B\) divides \(x_j\),
every marked-order term divisible by \(p\) in the tail's stripped
ramps is exactly \(x_j\), since all literal factors in that order
differ by at most \(jW\). The central values are then zero and
the remaining values are positive, using both orientations of the
selected pair. In the companion order \(q\), a tail ramp factor
\(F\) divisible by \(p\) would make \(jF-qx_j\) divisible by \(p\).
Its absolute value is at most \(jqW\), and it is nonzero by
\(j\nmid x_j\) and \(\gcd(j,q)=1\). The choice of \(B\) also
makes the stripped core factors and binomial coefficients units.
The middle-width test proves coprimality in this case.

A collision would therefore make both \(x_n,x_N\) \(B\)-smooth.
The positive integers \(U=Nx_n\) and \(V=nx_N\) would be distinct
\(B\)-smooth integers in \([nNH,nN(H+W)]\), by the mark residues
and \(B\ge J\). Their distance is at most \(J^2W\). Define
\(\Theta(x)=(\log_2 x)(\log_3 x)/\log_4 x\), which is increasing
for sufficiently large \(x\), and use a monotone adjustment on
the initial bounded range. Let
\[
       y(x)=\max(3,J,JM_C,J^2\kappa_J\Theta(x)).
\]
For sufficiently large \(H\), both \(U,V\) lie in the sequence
of positive integers with largest prime at most \(y\) of that
integer. By
[Stewart, Theorem 1 and (6)](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/on_sequences-2.pdf),
consecutive members \(u<v\) of that sequence have
\[
         v-u>\frac{u}{(\log u)^{\delta(c y(v))}},
         \qquad \delta(x)=\exp\!\left(\frac{x\log_2 x}{\log x}\right)
\]
for an effective constant \(c>0\) and sufficiently large \(u\).
Choose \(\kappa_J>0\) so that \(cJ^2\kappa_J<1/4\). With this
choice, the logarithm of \(\delta(c y(v))\) is at most
\(\frac12\log_2 v\) for all sufficiently large \(v\). Indeed
the ratio of that logarithm to \(\log_2 v\) tends to
\(cJ^2\kappa_J\), by the definition of \(\Theta\). These
elementary logarithmic inequalities have computable
thresholds by bounding each iterated logarithm in the ratio.
For \(u,v\) between \(U,V\) in either ordering, one has
\(u\ge H\) and \(v\le 2J^2H\) for large \(H\). Thus the
displayed lower bound on the gap exceeds \(H^{1/2}\), since
\[
          \delta(cy(v))\log_2 u
             \le (\log(2J^2H))^{1/2}\log_2(2J^2H)
                           =o(\log H).
\]
There would be consecutive members of the sequence between
the smaller and larger of \(U,V\). Their gap cannot exceed
\(J^2W\), which is less than \(H^{1/2}\) beyond another
computable threshold. This contradiction closes the residual
pairs and hence all degrees. It permits growing bands and
mixed gaps in the paired case, complementing the
[widening-band first-block theorem](MULTIVARIATE_NARROW_HEIGHT_BAND_NATIVE_MARKS.md).
