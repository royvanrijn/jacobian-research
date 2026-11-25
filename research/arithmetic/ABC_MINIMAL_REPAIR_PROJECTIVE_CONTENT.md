# Primitive row content and congruences after one-target repair

Let \(a_i+b_i=c_i\) be positive primitive triples. Fix an order \(M\)
and write
\[
 Z_i=\prod_{p\mid c_i}p^{(M-v_p(c_i))_+},\qquad
 (A_i,B_i,N_i)=Z_i(a_i,b_i,c_i).
\]
The repaired target \(N_i\) is \(M\)-full. The common content of its
row is precisely \(Z_i\). This content can absorb all the apparent
new divisibility in homogeneous determinant calculations. The
statements below concern exact descent of given congruences and a
specific Taylor divisor, rather than bounds for every use of repair.

## Exact descent of multihomogeneous congruences

Let \(P\) be any integral polynomial in the rows that is homogeneous
of degree \(d_i\ge0\) separately in row \(i\). If an integer \(D\)
divides \(P((A_i,B_i,N_i)_i)\), the corresponding assertion on the
primitive rows is exactly
\[
 \frac{D}{\gcd(D,\prod_i Z_i^{d_i})}
 \ \bigm|\  P((a_i,b_i,c_i)_i).
\]
This is immediate from the equality
\(P((A_i,B_i,N_i)_i)=(\prod_i Z_i^{d_i})P((a_i,b_i,c_i)_i)\),
and includes the case where the polynomial value vanishes.

For example, the common repaired target divisor gives
\(\gcd(N_i,N_j)\mid A_iB_j-A_jB_i\), since this minor is also
\(A_iN_j-A_jN_i\). The descended modulus is
\[
 \frac{\gcd(N_i,N_j)}{\gcd(\gcd(N_i,N_j),Z_iZ_j)}.
\]
At a common target prime let \(e_i=v_p(c_i), e_j=v_p(c_j)\),
\(z_i=v_p(Z_i),z_j=v_p(Z_j)\). Its descended exponent is
\[
 [\min(e_i+z_i,e_j+z_j)-z_i-z_j]_+
 =[\min(e_i-z_j,e_j-z_i)]_+
 \le \min(e_i,e_j).
\]
The bound on the right was already present on the primitive rows.
If \(e_i,e_j\le M/2\), the descended exponent is zero, although
\(p^M\) divides both repaired targets. The single-row congruence
\(N_i\mid A_i+B_i\) likewise descends exactly to \(c_i\mid a_i+b_i\).

## The uncentred Taylor divisor can already follow from content alone

Let the square matrix of evaluations on \(k\) repaired rows have
columns indexed by distinct monomials \(A^\alpha B^\beta\). For total
degree \(d\), let \(k_d\) be the number of such columns. Write
\[
 S=\sum_d d k_d,\qquad K=\sum_d \binom{k_d}{2}.
\]
Fix a prime \(p\) in all the target supports and suppose that
\(0<e_i=v_p(c_i)\le M/2\) on every row. Thus
\(v_p(N_i)=M\), but \(v_p(A_i),v_p(B_i)\ge M/2\).
An ordinary uncentred Taylor expansion from
\(B_i=-A_i\bmod p^M\) supplies divisibility of its determinant by
\(p^{MK}\): in each degree block, using the same Taylor order twice
gives proportional columns, so at least
\(0,\ldots,k_d-1\) are needed.

In fact that same divisor follows with no use of
\(A_i+B_i=N_i\): by row content each individual determinant term is
divisible by \(p^{\lceil MS/2\rceil}\). Since there are only \(d+1\)
distinct ordinary monomials of total degree \(d\),
\[
 S=\sum_d d k_d\ge\sum_d k_d(k_d-1)=2K.
\]
Thus for this entire class of evaluation matrices, even with mixed
total degrees, the naked common-slope Taylor divisor is implied by
row content. If all columns have the same degree \(d\), row-content
cancellation is exact: divide the \(i\)-th row by \(Z_i^d\) and the
remaining matrix is the primitive evaluation matrix. Accounting for
extra Taylor divisibility from the primitive congruences is then a
calculation at the original orders \(e_i\). The claim does not
disallow methods using additional arithmetic restrictions on the
primitive rows.

In particular this applies at the primes of bounded scaled original
target degree \(e_i L\le D_0\) whenever the repair is chosen with
\(ML\ge2D_0\). Large formal target fullness supplied by this
repair cannot by itself replace a congruence estimate at the
original scaled orders.

## Sharpness on actual primitive triples

Fix a finite set of distinct ordinary monomial columns as above, let
\(D=\max(\alpha+\beta)\), and take \(1\le e<M\).
For every prime \(p>1+M\max(D,1)\) there exist genuine positive
primitive triples with targets
\[
 c_i=p^e q_i^M,\qquad \gcd(p,q_i)=1,\qquad
 \gcd(q_i,q_j)=1\quad(i\ne j)
\]
whose repaired evaluation determinant has valuation exactly
\[
                  (M-e)S+eK .
\]
The \(a_i,b_i,c_i\) can all be taken comparable in height. Thus
any universal extra divisibility beyond the common row content in
this setting has to use the original order \(e\), or further side
conditions. This sharpness example does not prescribe small radicals
on \(a_i,b_i\).

For each row write \(b_i=-a_i+p^e q_i^M\) and consider the primitive
evaluation determinant first. In each degree \(d\), its lowest
Taylor order uses \(0,\ldots,k_d-1\), and the block coefficient
determinant is that of
\(\binom{\beta_j}{r}\) with \(0\le r<k_d\), up to signs.
This is the Vandermonde in the distinct \(\beta_j\) divided by
\(\prod_{r=0}^{k_d-1}r!\), hence a \(p\)-adic unit. After factoring
\(p^{eK}\) the lowest term is a unit times the evaluation determinant
of the distinct monomials \(a^{d-r}z^r\),
\(0\le r<k_d\), with \(z=q_i^M\bmod p\). All terms of higher
total Taylor order vanish modulo \(p\).

The grid
\[
        \mathbf F_p^\times\times(\mathbf F_p^\times)^M
\]
has more than \(D\) points in each coordinate, because the second
factor has size \((p-1)/\gcd(p-1,M)\). On this grid the lowest-order
monomials are linearly independent by interpolation in each
coordinate. Choose the required number of independent evaluation rows.
Choose distinct positive pairwise coprime \(b_i'\) in the
prescribed nonzero root classes modulo \(p\), recursively imposing
\(b_i'\equiv1\bmod\prod_{j<i}b_j'\). Take a common positive
multiple \(Q\) of \(p\) and of every nonzero difference
\(b_i'-b_j'\) and set \(q_i=Qv+b_i'\) for large positive \(v\).
Any common prime would divide \(Q\) and both \(b_i',b_j'\),
so the \(q_i\) are pairwise coprime. They are in the required
classes and \(q_i/q_j\to1\). For each \(i\) choose
\(a_i\equiv1\bmod q_i\) and its prescribed nonzero class modulo \(p\)
with \(c_i/3\le a_i\le2c_i/3\). This is possible since the
progression gap \(p q_i\) is \(o(c_i)\) as \(v\to\infty\). Then \(b_i=c_i-a_i>0\) and
\(\gcd(a_i,c_i)=1\), so the triple is primitive. Its least repair
factor is exactly \(p^{M-e}\). Consequently the repaired determinant
is the primitive determinant multiplied by
\(p^{(M-e)S}\), giving the asserted exact valuation.
