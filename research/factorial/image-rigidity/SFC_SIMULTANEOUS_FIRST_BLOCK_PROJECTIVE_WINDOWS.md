# First-block factors and projective separation for ordinary moment windows

Let \(a_1,\ldots,a_s\) be distinct nonnegative integer exponent columns,
\(s\ge2\), and fix a maximum order \(J\ge s\). Suppose the multisets of
entries of the columns are equal. In particular the pure factorial
coefficients \(\prod_i(ja_{i,\nu})!\) are equal as \(\nu\) varies, at
each order \(j\). The proof only needs this pure factorial equality.

Omit rows constant across the columns. Select a nonempty set of rows
in which the entries take just two values, \(\ell_i\) and
\(M_i=\ell_i+d_i\), with \(d_i>0\). Put
\[
 M=\max_{\rm selected}M_i,\quad X=\max_{\rm selected}d_i,\quad
 e_i=M-M_i,\quad
 F=\max_i(e_i+d_i),\quad G=\max_i(e_i+(1-1/J)d_i),
\]
where the last two maxima are also over the selected rows.
Fix an integer \(B\ge0\) with \(e_i\le BX\). Suppose one can choose
\(G<T\le F\) with \(T-G\ge\gamma X>0\). Call a selected row persistent
if \(e_i+d_i\ge T\). On these persistent rows, any two distinct
columns must differ in at least \(2R-1\) entries, for some integer
\(R\ge1\).

For a retained row put \(M_i=\max_\nu a_{i,\nu}\) and
\(d_i=M_i-\min_\nu a_{i,\nu}\). Define \(\mathcal L\) by summing
\(d_i\log(JM_i)\) over the retained rows, omitting rows with
\(JM_i\le HX\) and rows whose entries are all in either
\([M-G,M]\) or \([M-(B+1)X,M-F]\). Here \(H=J(B+1)+s\).

**Theorem.** Suppose
\[
 \begin{gathered}
 X\ge\max\{4,2/\gamma\},\qquad
 M\ge\max\{64,\,2(B+1)X,\,X^{3/2}\},\qquad
 \log X\ge144H/\gamma,\\
                     \mathcal L\le\frac{R\gamma X\log M}{64J^2}.
 \end{gathered}                                               \tag{1}
\]
Then no nonzero projective coefficient direction has vanishing
factorial moments throughout any window of \(s\) consecutive
positive orders at most \(J\). The conclusion holds over every
characteristic-zero field for the original columns.

For example, one may use selected blocks each consisting of \(s\)
rows with the same lower exponent \(\ell\) and gap \(d\).
Take the common selected upper exponent \(M=X^2\).
In such a block the \(i\)-th row has exponent \(\ell+d\) in column
\(i\) and \(\ell\) elsewhere. All the pure factorial coefficients
of this block agree. Arbitrary fixed retained columns with the
same exponent multiset may be included. At a common selected
upper exponent, a sufficiently large gap in one block gives a
window even when the gaps of the other selected blocks vary.
The retained fixed factors eventually fall below the prime cutoff.
The conclusion concerns windows bounded by the prescribed \(J\).

## Mixed coefficients have an extra marked factor

For each \(u\in[G,T)\cap\mathbb Z\), put \(y=M-u\).
At every degree \(1\le j\le J\), strip the common factorial
\((j\min_\nu a_{i,\nu})!\) from each row of the form.
In a selected row with \(u<e_i+d_i\), the factor \(jy\)
lies inside its first block of length \(d_i\) in the stripped
ramp. In fact the difference from the smaller endpoint at order
\(j\) is \(j(e_i+d_i-u)\), which lies in \((0,d_i]\)
by \(u\ge G\). In the other selected rows the factor has
already fallen below or at the stripped ramp's lower endpoint.
This is the same mechanism as the
[binary simultaneous first-block comparison](MULTIVARIATE_SIMULTANEOUS_FIRST_BLOCK_ANGULAR_SEPARATION.md).

In an active row \(i\), a coefficient indexed by
\(\alpha\in\mathbb Z_{\ge0}^s\), \(\sum_\nu\alpha_\nu=j\),
contains \(jy\) if at least one column with positive
\(\alpha_\nu\) has the higher exponent in row \(i\).
Let \(V_\nu\) be the set of active rows where column \(\nu\)
has its higher exponent. The number of copies of the factor is
therefore \(\left|\bigcup_{\alpha_\nu>0}V_\nu\right|\).
Put \(r=\min_\nu |V_\nu|\). Distinct columns differ in at least
\(2R-1\) persistent rows, and all persistent rows are active.
For two distinct indices \(v,w\),
\[
 |V_v\cup V_w|-\min(|V_v|,|V_w|)
   =\max(|V_v\setminus V_w|,|V_w\setminus V_v|)\ge R .
\]
Thus each mixed coefficient has at least \(r+R\) copies of \(jy\).

Take a prime \(p>HX\) dividing \(y\) and put \(E=v_p(y)\).
The selected factors other than \(jy\) are units at \(p\).
Indeed every stripped factor in their full ramps has signed
difference \(j(u-e_i)-h\) from \(jy\), with
\(0\le h<jd_i\), and the nonzero differences have magnitude
at most \(J(B+1)X\). The omitted retained rows are also units.
The smaller ramps lie below the prime cutoff. A retained row
in the upper permitted interval starts strictly above \(jy\),
and a row in the lower one ends strictly below \(jy\). The
same offset bound applies in either interval.
All multinomial coefficients are units since \(j<p\).

Write \(D_j\) for the product of the full stripped order-\(j\)
ramps of the retained rows charged to \(\mathcal L\).
The common pure coefficient \(b_j\) has valuation at most
\(rE+v_p(D_j)\) by the pure column attaining the minimum \(r\).
Every mixed coefficient has valuation at least \((r+R)E\).
If \(RE>\max_{1\le j\le J}v_p(D_j)\), dividing each stripped
form by its \(b_j\) gives simultaneous reductions
\[
                             \sum_{\nu=1}^s z_\nu^j .
\]

## A prime for all the orders and the residue argument

Such a prime exists by the
[progression product comparison](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root).
The \(y\)'s form a step-one progression of length
\(K\ge\gamma X/2\) with \(K\le X\), since \(F-G\le X/J\).
Every \(y>M/2\) has logarithm greater than
\((5/6)\log M\) by \(M\ge64\). If no prime succeeds, the
\(R\)-th power of the product of the large-prime parts of all
marks divides \(\prod_{1\le j\le J}(D_j)_{>HX}\).
This uses that each prime above \(HX\) divides at most one mark.
The latter product has logarithm at most \(J^2\mathcal L\)
by the order-\(j\) ramp lengths.

Removing one maximal-valuation mark at each smaller
prime leaves a small-prime part dividing \((K-1)!\).
The elementary bound on the number of primes at most \(HX\)
and (1) give at most \(K/8\) deletions.
Since \(\log K\le(2/3)\log M\), the logarithm of the
remaining product is both greater than
\((7/8)(5/6)K\log M\) and at most
\[
                 \log((K-1)!)+J^2\mathcal L/R
                               \le(2/3+1/32)K\log M .
\]
This is impossible.

A common projective zero in characteristic zero would give one
over \(\overline{\mathbb Q}\) because the equations have rational
coefficients. Take a place above the successful prime of a finite
field of definition and rescale the coordinates to be integral
there with at least one unit. Their reduction would give a nonzero
point where the displayed pure power sums vanish at \(s\) consecutive
orders. Let \(t\le s\) be the number of distinct nonzero coordinate
values \(\lambda_1,\ldots,\lambda_t\), with multiplicities \(m_i\).
Using the first \(t\) orders of the window, starting at \(n\), gives
the Vandermonde system with matrix \((\lambda_i^k)_{0\le k<t,1\le i\le t}\)
and unknowns \(m_i\lambda_i^n\). The matrix is invertible and the
unknowns are all nonzero because \(p>s\). This contradiction proves
the projective statement.
