# Truncated Newton Smith widths without prime factorization

The [gcd-free residual sieve](MULTIVARIATE_GCD_FREE_NEWTON_RESIDUAL_SIEVE.md)
compares the first Newton residual polynomials of two rational
coefficient arrays at all primes in their coefficient support.
The same coprime bases also support arbitrary precision comparisons.
A Sylvester matrix on a truncated Newton chart is a rational matrix
with monomial row and column scalings. This form is preserved by
Gaussian elimination, so an integral Smith pivot can be chosen at
all primes of a composite modulus using rational arithmetic and gcd
splits. The resulting bounds apply to the original factorial profile,
with arbitrary row sizes, multiplicities and orientations.

## A common Newton chart

Use the coefficients \(c_{j,k}\) of the two constant-normalized
moments at adjacent orders \(j=n,N\). For one coprime base \(B\),
let \(e_{j,k}\) be its exponent in \(c_{j,k}\), and let \(a/b\)
be a common lower polygon slope for these raw ordinates, in
lowest terms with \(b>0\). Set
\[
       m_j=\min_k (b e_{j,k}-a k),\qquad
       h_{j,k}=b e_{j,k}-a k-m_j .
\]
Choose an integer precision \(L\ge1\) and put
\(l_j=\min\{k:h_{j,k}<L\}\).
In an algebra over \(\mathbb Q\) with \(T^b=B\), form
\[
       A_j(X)=
       \sum_{k:\,h_{j,k}<L}
                         c_{j,k}T^{-a k-m_j}X^{k-l_j}.        \tag{1}
\]
The coefficients are \(T^{h_{j,k}}\) times rational units at
every \(p\mid B\). At such a prime choose any root \(T\) of
\(T^b=B\) in \(\overline{\mathbb Q}_p\) and use the finite
extension \(K=\mathbb Q_p(T)\). With \(v=v_p(B)\), write
\(\omega=b v_p/v\), so \(\omega(T)=1\).
No irreducibility or ramification restriction on \(T^b-B\)
is needed.

Compute the Smith entry valuations \(\lambda_i\) of the
Sylvester matrix of \(A_n,A_N\) using \(\omega\).
Every rational common factor of the original pair has at most
\[
                          \#\{i:\lambda_i\ge L\}             \tag{2}
\]
roots of valuation \(-va/b\), counted with multiplicity.
The precision ideal here is \((T^L)\) in the valuation
ring \(\mathcal O\) of \(K\). If either \(A_j\) is constant,
the bound is zero.

Indeed, in both orders substitute \(z=T^{-a}X\) and divide
by \(T^{m_j}\). The factors on the \(r\) roots of the proposed
common factor with valuation \(-va/b\) become the same
monic divisor \(H(X)\) in \(\mathcal O[X]\), with all roots
units. To see that its coefficients are in \(K\), the multiset
of roots of the specified valuation is invariant under the
absolute Galois group of \(K\). Monic division gives divisibility
over \(\mathcal O\). Reduce modulo \(T^L\), discard terms of
height at least \(L\), and divide by \(X^{l_j}\), which is
invertible modulo \(H\) because its constant coefficient
is a unit. Each \(A_j\) is nonzero in the residue field,
by construction, so reducing first modulo the maximal ideal
gives \(r\le \min(\deg A_n,\deg A_N)\).
The truncated Sylvester cokernel over \(\mathcal O/(T^L)\)
surjects onto \((\mathcal O/(T^L))[X]/(H)\), which is free
of rank \(r\). If \(T^L\) has uniformizer valuation \(e\),
multiply by the \((e-1)\)-st power of a uniformizer and
compare dimensions over the residue field to obtain (2).
This is the local truncated-face bound at the chart's
common normalized precision. [Physical comparisons of truncated Newton faces](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)

## Simultaneous Smith pivoting

Index a row of the Sylvester matrix by output degree \(i\)
and a column by the order \(j\) and the input multiplier
degree \(u\). A nonzero entry comes from \(k=i-u+l_j\)
and is
\[
      T^{r_i} q_{i,(j,u)} T^{s_{j,u}},\qquad
      r_i=-a i,\quad s_{j,u}=a(u-l_j)-m_j,\quad
      q_{i,(j,u)}=c_{j,k}.                                   \tag{3}
\]
For a zero entry put \(q_{i,(j,u)}=0\).
This follows from the coefficient indexing of (1), even if
there are gaps between retained coefficients. For any two
nonzero entries their ratio to the \(b\)-th power is
rational. Its reduced denominator records exactly the primes
where the numerator entry has smaller valuation. Therefore
one can compare two candidate pivots simultaneously over
a modulus \(M\) made of whole prime powers of \(B\): for
denominator \(D\), compute
\[
                \gcd(M,D^{\,1+\lfloor\log_2 M\rfloor}).
\]
On its support the numerator entry is strictly smaller.
At all other primes of \(M\) it is no smaller. Both parts
retain whole prime powers. A scan of the entries, splitting
as needed, selects an entry of minimum valuation for
every prime of each resulting block. Ties can be resolved
by keeping the current pivot.

Permute the minimum pivot to the first position. Subtracting
row and column multiples clears its row and column using
integral operations at every prime of the block. The
Schur complement has entries of the same form (3), with
each coefficient replaced by the ordinary rational Schur
complement of the \(q\)-matrix. The remaining valuations
are at least that of the pivot, because each original entry
has at least that valuation, as does each correction term
formed by dividing a product of two entries by the pivot.
Repeating this procedure yields a diagonal matrix whose
valuations are nondecreasing at every prime of the block.
Thus the pivot valuations are exactly the Smith entry
valuations. A zero remaining matrix gives infinite
valuations for the remaining entries.

Every pivot is stored via its rational \(b\)-th power
\(R=B^{r_i+s_{j,u}}q_{i,(j,u)}^b\). Its contribution to
(2) is equivalent to \(v_p(R)\ge v_p(B^L)\), which is
decided by the same saturated gcd on the denominator of
\(R/B^L\). This includes different values of \(v_p(B)\)
at different primes of one block. Thus the entire vector
of threshold tests, and therefore the exact cap (2), is
obtained on a coprime partition of \(B\) without prime
factorization. The number of pivoting steps per block is
bounded by the matrix size, and there are at most as many
final blocks as prime divisors of \(B\).

There is an equivalent expression using rational minors.
For row and column sets \(I,J\) of size \(t\), define
\[
       R_{I,J}=B^{\sum_{i\in I}r_i+\sum_{(j,u)\in J}s_{j,u}}
                                  \det(q[I,J])^b .
\]
It is the \(b\)-th power of the corresponding chart minor.
If \(G_t\) is the gcd of the reduced numerators of the
\(R_{I,J}\), with \(G_0=1\), then while \(G_t\ne0\),
\[
             \sum_{i=1}^t \lambda_i=\frac{v_p(G_t)}{v_p(B)},
             \qquad
             \lambda_t\ge L
               \quad\Longleftrightarrow\quad
             v_p(G_t)\ge v_p(B^L G_{t-1}).
\]
If \(G_t=0\), \(\lambda_t\) and all later Smith valuations are infinite.
This formula also proves independence of the chosen root
of \(T^b=B\). All rational minors have nonnegative valuation
at the primes under consideration after the indicated
scaling, since (1) is integral there.

## Degree and projective consequences

One can refine the coprime blocks across several slopes
and several precisions, and at each slope take the minimum
of its first residual width and the bounds (2). Increasing
\(L\) generally changes the truncated matrix, so these
separate caps can be combined without assuming monotonicity.
The reduced slope denominator is still a divisor of the
multiplicity at each prime. Hence these narrower widths
can be used in the same block allocation, perfect-power,
Jacobi and projective denominator tests as the first
residual widths. Supplemental moduli with horizontal
polygons also admit the method by taking \(b=1\), zero
raw ordinates, and a modulus supported on unit-coefficient
primes. The primes outside the coefficient support still
contribute zero to projective denominators.

The \(p=17\) ordinary profile in the
[exact precision example](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md#an-exact-precision-example)
gives strict improvement over first residual widths.
At orders \(3,4\) the adjacent valuation differences are
\((0,2,7)\) and \((0,3,4,7)\), so \(v_{17}(B)=1\)
in the coprime base containing \(17\). The two shared
slopes \(0,7\) each have first residual width one, and
these widths permit a quadratic. At raw slope \(7\) and
\(L=2\), both polynomials (1) are linear and have unit
leading coefficients. The ratio of their normalized
last-block ratios is \((3/4)^{16}\), whose difference
from one has valuation one. The Smith valuations are
\(0,1\), so (2) deletes this slope. At most one root
can remain at slope zero, excluding a common rational
factor of positive even degree.

These are adaptive tests on the original factorial
profiles. If every required pair below the absolute
eventual cutoff is excluded, they certify all degrees
of that profile. A proof ensuring that this happens
for every unrestricted profile is still open.
