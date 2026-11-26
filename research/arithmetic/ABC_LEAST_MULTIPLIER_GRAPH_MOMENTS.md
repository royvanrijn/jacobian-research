# Graph moments of actual least-repaired summands

Use genuine positive primitive rows \(a+b=c\) in
\(X/2<c\le X\), with all three original caps.
Write \(n=X^L\), \(M=\lceil K/L\rceil\),
\(Z=Z_M(c)\), and \(U=aZ\).
The height limit is taken with \(K,L\) fixed.
Use the balanced target cap and parameter choices of LM1.
Write \(\gamma=1/6-\Delta\) for its exponent at a
prescribed nonzero weighted difference. Its proof also
gives that exponent for any prescribed nonzero signed
weighted sum. The primitive coefficient form of SHF1
allows either sign, and the supported-factor argument
for LM1 works without change. Each exact \(Z\) has at most
\(n^{1/K}X^\varepsilon\) genuine target preimages.
[Fixed-shift least multiplier count (LM1)](ABC_LEAST_MULTIPLIER_SHIFT_COUNTS.md#a-strict-half-bound-at-a-prescribed-weighted-multiplier-shift)

## Collective differences across all original summand fibres

**LG1.** Let \({\cal E}\) have \(N\) actual rows and \(A\)
complete original \(a\)-values. Let \({\cal D}\) be the
union of the nonzero exact differences \(U_1-U_2\) over
ordered pairs with \(a_1=a_2\). Then
\[
 |{\cal D}|\, n^\gamma X^\varepsilon
            \ \gg_\varepsilon\
       N^2/A-O_\varepsilon(N n^{1/K}X^\varepsilon),
\]
where the two height allowances can be chosen separately.
For a fixed nonzero difference \(F\), the common original
\(a\) divides \(F\). Choose it at divisor cost and apply
LM1 to \(Z_1-Z_2=F/a\). Both original rows follow from
their targets and \(a\). Thus \(F\) has at most
\(O(n^\gamma X^\varepsilon)\) preimages across the whole
family. On the other hand, writing \(r_a\) for the actual
original degree, the total number of pairs is
\(\sum_a r_a^2\ge N^2/A\). Equal \(U\) at fixed \(a\) means
equal \(Z\), with at most \(n^{1/K}X^\varepsilon\) original
targets per first row.

In particular at row exponent \(d\) and complete
\(a\)-dictionary exponent \(\alpha\), if
\(d-\alpha>1/K\) by a fixed margin then the actual
nonzero difference catalogue has exponent at least
\(2d-\alpha-\gamma\) up to the inner allowances.
At balanced saturation with \(\alpha=1/6+o(1)\) this gives
\(|{\cal U}-{\cal U}|\ge n^{1/3+\Delta-o(1)}\) on the
entire repaired \(a\)-dictionary \({\cal U}\), whose
cardinality is at most \(N\). The global union is used
rather than the size of a separately chosen fibre.

For any fixed \(k\ge3\), the same method works for
nonzero signed words \(\sum_{i=1}^k w_i U_i\) in a
common original \(a\), with prescribed nonzero integral
weights of polynomial height. Fix \(a\mid F\), then all
but two positions with weights of the same sign. Such
positions always exist. Their residual weighted sum of
positive multipliers is nonzero, so LM1 applies. If
\(n^u\) up to allowances is an original upper bound on
the degree of a complete \(a\), the nonzero word
multiplicity is at most \(n^{(k-2)u+\gamma}X^\varepsilon\).
The zero word has at most \(A n^{(k-2)u+\gamma}X^\varepsilon\)
preimages by fixing \(a\) directly. Thus whenever
\(k d-(k-1)\alpha>\alpha+(k-2)u+\gamma\) by a fixed
margin, the nonzero word dictionary exponent is at least
\(k d-(k-1)\alpha-(k-2)u-\gamma\), up to allowances.
For positive weights no subtraction for zero is needed.

At occurrence there is a different one-record statement.
For each ordered pair in the full family, put
\(F=\operatorname{lcm}(a_1,a_2)(Z_1-Z_2)\).
When \(F\ne0\), both \(a_i\) divide \(|F|\). Their choices
cost height-subpower and their least common multiple
makes the multiplier shift exact. Thus this scalar also has
at most \(O(n^\gamma X^\varepsilon)\) preimages. In a
balanced near-full family, the pairs with equal \(Z\)
have count \(O(N n^{1/K+1/6+\tau/2+\eta}X^\varepsilon)\)
by reconstructing the other target and using FL1.
For sufficiently large \(K\) and small tolerances this
is power-small relative to \(N^2\). Hence in every
balanced near-full family, including occurrence, these
actual scalar outputs have exponent at least
\(1/2+\Delta-o(1)\).
[Original fixed-target fibre bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Determinants on genuine repaired grids

**LG2.** Form the bipartite graph with complete original
\(a\)-labels on one side and actual multiplier labels \(Z\)
on the other. For each \(Z\) choose one genuine target \(c\)
with maximum incidence in that fibre. This costs at most
\(n^{1/K}X^\varepsilon\) in the number of rows and makes
the graph simple. At each retained edge the repaired
summand is \(U=aZ\) and the other original summand is
\(c-a\).

For an ordered complete \(k\)-by-\(k\) grid with distinct
labels on each side, let \(f_k(x)=\sum_{r=0}^{k-1}x^r\).
The determinant of its actual repaired evaluations is
\[
 \det\bigl(f_k(a_i Z_j)\bigr)_{i,j=1}^k
     =\prod_{i<j}(a_j-a_i)\prod_{i<j}(Z_j-Z_i).
\]
Indeed the matrix is the product of two \(k\)-by-\(k\)
Vandermonde matrices. For a fixed nonzero determinant, all
signed offsets from \(a_1\) and from \(Z_1\) are divisors
of the determinant. At fixed \(k\) they have subpower
enumeration. Choose one actual \(a\)-anchor. One nonzero
fixed \(Z\)-offset, by LM1, then determines the two \(Z\)
values at exponent \(\gamma\). The remaining labels follow
from the offsets. Thus each determinant has at most
\(O(A n^\gamma X^\varepsilon)\) preimages among the ordered
complete grids.

Write \(d,\alpha,\beta\) for the retained edge and part
exponents, with the actual target selection cost already
subtracted in \(d\). The repeated Jensen bound for ordered
homomorphisms of \(K_{k,k}\) has exponent at least
\[
                     h=k^2d-(k^2-k)(\alpha+\beta).
\]
Indeed apply convexity first to the number of common right
neighbors of \(k\) left vertices, then to the right degrees.
Repeated labels account for at most
\(O_k(A^{k-1}|{\cal Z}|^k+A^k|{\cal Z}|^{k-1})\) grids.
When \(h>\max((k-1)\alpha+k\beta,k\alpha+(k-1)\beta)\)
by a fixed margin, the nonzero determinant catalogue
has exponent at least \(h-\alpha-\gamma\), up to the
height allowances. At balanced minimum dictionaries
\(\alpha=\beta=1/6+o(1)\), for any fixed \(k\) and for
\(K\) large enough to absorb the target selection loss,
this lower bound approaches \((k-1)/3+\Delta\).
This is a restriction on the exact arithmetic output
of the original complete grids with their actual repairs.

## A polynomial graph cannot model cheap repair transport

**LG3.** Let \(F(C,Z)\) be a nonzero polynomial of total
degree at most \(D\ge1\) with rational coefficients of
polynomial height. Then uniformly over \(F\), the number
of integer targets \(c\le X\) with \(F(c,Z_M(c))=0\) is
\[
                   O_\varepsilon(X^{(D+1)/M+\varepsilon}).
\]
Clear denominators and prescribe the coefficient prime
pool. At a prime outside this pool with original exponent
\(0<e<M\), the corresponding \(Z\)-exponent is \(M-e\).
A vanishing sum needs two monomials with the same minimal
valuation. Distinct monomials can tie only when their
\(C\)-orders and \(Z\)-orders change in opposite directions.
Write their positive differences as \(r,s\) respectively.
The tie gives \(e/M=s/(r+s)\), with \(1\le r,s\le D\).
In particular \(e\ge M/(D+1)\). The remaining nondeficient
primes already have \(e\ge M\). Thus the outside support
radical is at most \(X^{(D+1)/M}\), and uniform fixed-support
enumeration finishes. In fact the admissible \(e/M\) lie
in the finite set of ratios from the monomial support
of \(F\). For an affine line with nonzero intercept the
constant term is uniquely minimal at every deficient prime
outside the coefficient pool. The count improves to
\(O_\varepsilon(X^{1/M+\varepsilon})\) on these lines.

The degree-\(D\) exponent is attained at the level of the
integer repair map when \(D+1\mid M\). Take
\(F(C,Z)=C^D-Z\), \(c=R^{M/(D+1)}\) for squarefree \(R\),
and then \(Z=c^D\). This construction is about the target
map. In a genuine near-full abc family, FL1 requires the
actual target dictionary to have exponent at least \(1/6-o(1)\).
A cover of its repair graph by polynomial curves of fixed
degree \(D\) therefore has catalogue price at least
\(1/6-(D+1)/K-o(1)\). This excludes cheap exact polynomial
transport of the complete original target dictionary.
[Original fixed-value count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
