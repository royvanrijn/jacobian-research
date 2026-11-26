# Multi-degree matched tanh moments do not determine native prime certificates

Fix a finite collection of nonnegative integer monomials in the
positive Gauss tanh values \(T_j\) at any finite collection of
orders \(j\ge2\). Independently fix a nonempty finite set of target
upper orders \(m\ge3\). There are two ordinary profiles with the
same specified moment sums, the same number of shared rows, and
arbitrarily small shared curvature, with the following property.
For each target \(m\), a prime \(p_m\) gives an actual native prime
certificate of angular separation for the first profile at
\((m-1,m)\), whereas the degree-\(m\) Newton polygon of the second
profile is horizontal at the same prime.

Both profiles consist of a private unit pair and increasing
shared unit-gap rows. Their row multiplicities are fixed along a
common integer progression. Equality holds for the signed as well
as unsigned moment sums since the orientations coincide. Taking
all targets \(3\le m<J\) makes the first profile native-certified
in every degree using the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
for the remaining orders. The second profile may have other
certifying primes, and the construction gives no angular collision.

This is a strengthening of the
[two-order positive-moment obstruction](MULTIVARIATE_FINITE_MATCHED_POSITIVE_MOMENTS_NATIVE_PRIME_OBSTRUCTION.md).
It uses finite differences along the poles corresponding
to the observed orders, rather than a two-dimensional
array of interpolation nodes.

## Exact identities by finite differences

Select any reference \(n\ge2\) and set \(N=n+1\), \(L_0=2nN\).
Parametrize a unit-gap row by integer coordinates \(u,v\)
with \(D=v-u>0\) and an integer parameter \(B\) through
\[
              a(B)=\frac{B-u}{L_0 D}-\frac{N}{2n}.            \tag{1}
\]
The
[exact Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
for \(d=1\) gives \(T_j=(j-1)/(2ja+j+1)\). Define
\[
                 \lambda_j=N-\frac{nN}{j},\qquad
                 w_j=u+\lambda_j D .
\]
Then \(T_j\) is \(D/(B-w_j)\) times the nonzero constant
\(L_0(j-1)/(2j)\). The different \(j\) have different
\(\lambda_j\).

Let \(\mathcal J\) be the finite set of observed orders, adding
one order with zero exponent if needed. For each \(j\in\mathcal J\)
let \(k_j\) be one plus the largest exponent of \(T_j\) among
the specified monomials. Choose positive integers \(h_j\)
with \(\lambda_j h_j\) integral. Form a grid indexed by
\(\mathbf t=(t_j)_{j\in\mathcal J}\), \(0\le t_j\le k_j\), with
\[
     u_{\mathbf t}=-\sum_j \lambda_j h_jt_j,\qquad
     D_{\mathbf t}=D_0+\sum_j h_jt_j,\qquad
     v_{\mathbf t}=u_{\mathbf t}+D_{\mathbf t}.
\]
Give this cell the signed weight
\[
                    W_{\mathbf t}
                     =\prod_j(-1)^{t_j}\binom{k_j}{t_j}.
\]
The two profiles use the positive and negative cells,
with multiplicities \(|W_{\mathbf t}|\) respectively.
Their counts agree since the total signed weight is zero,
and the corner \(\mathbf t=0\) has positive weight one.

For a specified monomial with powers \(r_j\), its value
at a cell is a constant times
\[
               \frac{D_{\mathbf t}^{\sum_j r_j}}
                    {\prod_j(B-w_j)^{r_j}}.                  \tag{2}
\]
To inspect the signed sum at a pole corresponding to order
\(j\), fix all coordinates except \(t_j\).
The pole \(w_j=u_{\mathbf t}+\lambda_jD_{\mathbf t}\)
is independent of \(t_j\). At each cell the difference from
any other order's pole is
\(w_i-w_j=(\lambda_i-\lambda_j)D_{\mathbf t}\).
The principal part coefficient of \((B-w_j)^{-r_j+h}\)
for \(0\le h<r_j\) is therefore a constant times
\(D_{\mathbf t}^{r_j-h}\). The constant is independent of
the cell: the remaining denominator factors in (2) are
expanded at \(w_j\), and each successive term costs one
power of the pole difference, hence of \(D_{\mathbf t}\).
This is a polynomial in \(t_j\) of degree at most \(r_j\).
The finite difference of order \(k_j>r_j\) with the weights
\((-1)^{t_j}\binom{k_j}{t_j}\) annihilates it.
This cancels every principal part, even if several groups
of cells happen to have the same pole. For a monomial of
positive total degree the resulting rational function of
\(B\) tends to zero at infinity, so it is identically zero.
For total degree zero the signed counts already sum to
zero. This proves all prescribed matched moment equalities.

## Ensuring integral rows and isolated native markers

For a target \(m\), set
\[
          \lambda_m^*=N^2-\frac{L_0}{m}>N,\qquad
          g_{m,\mathbf t}=u_{\mathbf t}+\lambda_m^*D_{\mathbf t}.
\]
Formula (1) gives
\[
           m a_{\mathbf t}(B)+h
                =\frac{m}{L_0D_{\mathbf t}}
                      \left(B-g_{m,\mathbf t}
                                +\frac{L_0(h-1)}{m}D_{\mathbf t}\right).
                                                               \tag{3}
\]
Let \(C\) be the number of cells. Choose the steps \(h_j\)
in an arbitrary fixed order, each divisible by \(L_0\),
by the denominator of \(\lambda_j\) with a further factor
\(L_0\), and by every prime at most \(C\).
Choose each step sufficiently large that its own increment
exceeds the sum of maximal preceding increments, both for
the sums \(\sum t_j h_j\) and for every target sum
\(\sum t_j(\lambda_m^*-\lambda_j)h_j\).
All the weights \(\lambda_m^*-\lambda_j\) in these latter
sums are positive. Hence each family of sums is distinct
across the cells. The coordinates \(u_{\mathbf t}\) are
integers divisible by \(L_0\).

There is a choice \(D_0\equiv1\bmod L_0\) for which the
positive integers \(D_{\mathbf t}\) are pairwise coprime.
For each prime dividing a difference of two of the distinct
sums \(\sum t_jh_j\), avoid all negatives of those sums
modulo the prime. At a prime at most \(C\) all the sums
are zero, so this excludes just one residue. At a larger
prime there are fewer than the prime number of residues
to avoid. At primes dividing \(L_0\) the required residue
one works. The Chinese remainder theorem gives arbitrarily
large such \(D_0\).

For a fixed target \(m\), every rational number
\[
 g_{m,0}-g_{m,\mathbf t}
                 +\frac{L_0(h-1)}{m}D_{\mathbf t},
 \qquad
 1\le h\le m,\quad (\mathbf t,h)\ne(0,1),                     \tag{4}
\]
can be made nonzero. For \(h=1\) this follows from
distinctness of the corresponding target sums. For \(h>1\)
the last term grows with \(D_0\), whereas the first
difference is independent of \(D_0\). Choose \(D_0\)
sufficiently large for all the finitely many conditions.
Then choose distinct primes \(p_m>m+1\) avoiding
\(L_0\), every \(D_{\mathbf t}\), and the numerators and
denominators of (4) for their target order.

The conditions
\[
 \begin{split}
 B&\equiv N^2 \pmod {L_0},&
 B&\equiv u_{\mathbf t}\pmod {D_{\mathbf t}}\quad
                                \text{for every cell},\\
 B&\equiv g_{m,0}+p_m \pmod {p_m^2}\quad
                                \text{for every target }m
 \end{split}
\]
are compatible by the Chinese remainder theorem.
Rational residues at the primes are interpreted using
their unit denominators. The first two conditions give
\(B-u_{\mathbf t}-N^2D_{\mathbf t}\equiv0
\bmod L_0 D_{\mathbf t}\), so (1) is integral.
On sufficiently large positive members of the progression
all row bases are positive and grow linearly with \(B\).
The fixed row multiplicities give shared curvature
tending to zero.

At \(p_m\), equation (3) gives valuation one for
\(m a_0+1\), and valuation zero for every other sampled
factor, by (4). The positive profile has this single
marked row with multiplicity one. After subtracting the
row scalars, its degree-\(m\) coefficient valuations are
\((0,1,\ldots,1)\), because the private pair cancels the
binomial coefficient up to \(m!\), a unit at \(p_m\).
This is the first case of the
[native prime chord certificate](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md).
The negative profile has no marked factor, so all its
normalized upper coefficients are units at that prime.
