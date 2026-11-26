# Factorial quotient transfer by borrowing across exponent digits

A row of positive ordinary exponents obtained by borrowing across
a base-\(p\) digit gap transfers a quotient of two shifted factorial
rows to the actual moment coefficients. The comparison preserves
the valuation of each row contribution, up to the exact stated
scalar and diagonal factors. It works at arbitrary prescribed
precision below any finite moment-order bound.

A consequence is an all-degree closure property for ordinary
two-monomial profiles. One can append arbitrarily many balanced
cycles of borrowed rows to any all-degree regular core. This
applies to fully shared profiles and allows arbitrarily many
unequal rows, with both signs of gap or with gaps all of one sign.
The proof uses the absolute eventual theorem only above its
existing cutoff.

## The borrowed factorial identity

Fix a prime \(p\), a positive integer \(J\), and a desired precision
\(p^a\), \(a\ge1\). Let \(A,B\) be two integer rows with positive
entries in every column. For a nonnegative count vector \(m\) with
\(1\le\sum_\nu m_\nu<J\), put \(x=\sum_\nu A_\nu m_\nu\) and
\(y=\sum_\nu B_\nu m_\nu\). Thus \(x,y>0\). Take integer scales \(e>s\ge0\)
and assume
\[
  s+1\ge a,\qquad
  e-\lfloor \log_p(J\max_\nu B_\nu)\rfloor\ge a.           \tag{1}
\]
In particular \(p^e>y\), and \(p^e A-B\) is an ordinary row of
positive exponents. Define
\[
 F_p(M)=\prod_{\substack{1\le h\le M\\p\nmid h}}h,\qquad
 \Lambda_b=p^{1+p+\cdots+p^{b-1}}\prod_{h=1}^b F_p(p^h),
 \qquad \Lambda_0=1 .
\]
Successively removing \(p\)-divisible factors gives
\[
 (p^e x)!=(p^s x)!
   (\Lambda_e/\Lambda_s)^x\prod_{h=s+1}^e U_h(x),\qquad
 U_h(x)=F_p(p^h x)/F_p(p^h)^x .                           \tag{2}
\]
This is the same inflation identity used in
[higher-precision row inflation and packing](SFC_PADIC_ROW_INFLATION_HIGH_PRECISION.md).
Every \(h\) in (2) satisfies \(h\ge a\). The unit residues in
each successive interval of length \(p^h\) have the same
product modulo \(p^a\), so \(U_h(x)\equiv1\bmod p^a\).
This complete-block argument works at every prime.

The last \(y\) factors in \((p^e x)!\) have product
\[
 p^e x\,(-1)^{y-1}(y-1)!
       \prod_{h=1}^{y-1}(1-p^e x/h).
\]
All factors in the final product are \(1\bmod p^a\) by (1),
including when \(p\mid h\). Dividing (2) by this product
proves the exact comparison
\[
 (p^e x-y)!
   =-p^{s-e}(\Lambda_e/\Lambda_s)^x(-1)^y
       \frac{(p^s x-1)!}{(y-1)!}\,U(x,y),\qquad
                    U(x,y)\in1+p^a\mathbb Z_p .            \tag{3}
\]
The product over \(h\) is empty at \(y=1\). Equation (3)
does not require either factorial in the quotient to be a
unit. In fact it follows by exact division in \(\mathbb Q\)
and gives the exact valuation after scalar and diagonal
normalization.

In a moment polynomial, the factor
\((\Lambda_e/\Lambda_s)^x(-1)^y\) is the monomial with
exponent vector \(m\) evaluated at the entries
\((\Lambda_e/\Lambda_s)^{A_\nu}(-1)^{B_\nu}\).
One inverse rational diagonal scaling of the coefficient
variables therefore removes it at every order. The scalar
\(-p^{s-e}\) can also be removed. This identifies the
virtual contribution of the borrowed row with
\((p^s x-1)!/(y-1)!\), up to a quotient that is a unit
and one modulo \(p^a\). Separate borrowed rows can be
combined by multiplying these diagonal normalizations.
If the virtual array has denominators at \(p\), compare the
arrays after dividing both by the same minimal coefficient
valuation.

## Balanced cycles of ordinary rows

Take a finite directed multigraph with equal indegree and
outdegree at every vertex, counting repetitions. Assign a
positive exponent row \(A_v\) and an integer \(s_v\ge a-1\)
to vertex \(v\). At edge \(v\to w\), append the ordinary row
\[
                         p^{e_{vw}} A_v-p^{s_w}A_w        \tag{4}
\]
with \(e_{vw}>s_v\) and
\(e_{vw}-\lfloor\log_p(Jp^{s_w}\max_\nu A_{w,\nu})\rfloor
 \ge a\). In particular all rows in (4) have positive
entries. With \(x_v=\sum_\nu A_{v,\nu}m_\nu\), each edge
has virtual contribution
\((p^{s_v}x_v-1)!/(p^{s_w}x_w-1)!\).
Their product is one, by balance. After the joint
diagonal and scalar normalization, every coefficient of
the augmented moment is the corresponding core coefficient
times a unit in \(1+p^a\mathbb Z_p\), for all orders
\(1\le j<J\). The number of vertices or edges places no
additional demand on the precision.

The same calculation for an unbalanced graph leaves the
explicit quotient with each vertex factorial raised to
its outdegree minus indegree. It can therefore compare
ordinary profiles with prescribed virtual shifted
factorial quotients at the specified orders.

## Application to every moment order

Let the core consist of two distinct monomials and
suppose its consecutive moment polynomials are coprime in
every positive pair. Let \(J\) be the
[effective absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md),
so every ordinary profile is regular at terminal \(N\ge J\).
Choose any prime \(p\). For the finitely many core pairs
below \(J\), divide each parent by the \(p\)-power of its
minimum coefficient valuation. The resulting
\(\mathbb Z_p\) polynomials have nonzero resultants.
Choose \(a\) larger than all of their resultant valuations.

Now append any balanced collection of (4), using the given
\(p,J,a\). The normalized augmented arrays agree modulo
\(p^a\) with the corresponding normalized core arrays.
The nominal degrees used in the resultant are the same,
even if a leading coefficient is a nonunit. Its Sylvester
determinant is an integer polynomial in those coefficients,
so the augmented resultant is congruent to the nonzero core
resultant at precision greater than its valuation. Hence
each augmented pair below \(J\) is coprime. For \(N\ge J\)
the absolute theorem applies to the actual ordinary
augmented support. The diagonal and scalar changes are
rational and invertible, so coprimality holds over every
characteristic-zero field.

For example, start with the ordinary single shared row
\((1,2)\), which is regular by the
[unit-gap Gamma comparison](MULTIVARIATE_SHARED_UNIT_GAP_GAMMA_NUMERICAL_RANGE.md).
Take a two-edge cycle on \(A=(1,3)\) and \(B=(4,5)\).
Its two rows are \(p^e A-p^{s_B}B\) and
\(p^f B-p^{s_A}A\). For large enough \(e,f\), both
gaps are strictly positive, and their sizes can be
chosen independently. Arbitrarily many independently
chosen cycles yield arbitrary dimension in a fully shared
ordinary family. Taking one of the vertex gaps negative
instead gives both signs.

This transfer is conditional on the regularity of the
comparison core or virtual quotient arrays. The large
digit gaps enforce the comparison in (3), and an arbitrary
ordinary row need not admit such a representation at the
selected precision. Conversely, knowing coprimality of
an augmented profile without a valuation bound for its
resultants does not prove the corresponding congruent
core coprime.

## Exact checks

The [standard-library checker](verify_multivariate_borrow_cycle_row_transfer.py)
tests the exact rational units in (3) at several primes,
including high precision at \(2\) and \(3\), and checks
both balanced cancellation in a three-vertex cycle and
the predicted quotient for an open path. These are
transcription controls for the written uniform identities.

    python3 research/factorial/image-rigidity/verify_multivariate_borrow_cycle_row_transfer.py
