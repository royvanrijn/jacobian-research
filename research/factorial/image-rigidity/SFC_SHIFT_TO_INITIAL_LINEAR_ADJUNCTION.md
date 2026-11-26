# Exact reduction of shifted Strong Factorial windows to initial windows

Fix an ordinary support of size \(s\) and let \(f\) use all its
monomials over a characteristic-zero field. Write
\(a_j=\Lambda(f^j)/j!\) and \(a_0=1\).
Fix \(k\ge0\). Over an algebraic closure the following are
equivalent:

- \(a_k\ne0\) and \(a_{k+1}=\cdots=a_{k+s}=0\).
- There are \(k\) nonzero coefficients \(d_i\) of fresh independent
  linear monomials such that
  \(g=f+\sum_{i=1}^k d_iY_i\) has its first \(k+s\) moments zero.

For \(k=0\) no new monomials are needed. The full multiset of the
\(d_i\) in the second statement is determined by \(f\). More
generally, if a length-\(s\) window at positive start \(n\)
vanishes, at most \(n-1\) extra linear monomials give an initial
zero window of their actual enlarged support size.
Thus general Strong Factorial in arbitrary finite dimension
is equivalent to its initial-window version in arbitrary finite
dimension. For a fixed support the shifted problem is exactly
the initial-window problem for that support with independent
linear monomials adjoined.

## A triangular ideal over the symmetric coefficient ring

Let \(h(t)=1+h_1t+\cdots+h_kt^k\) with independent
formal coefficients. For independent linear adjunction the
exponential moment series satisfy
\[
 A(t)=h(t)B(t),\quad
 A(t)=\sum_{j\ge0}a_jt^j,\quad B(t)=\sum_{j\ge0}b_jt^j,\qquad
 h(t)=\prod_{i=1}^k(1-d_it).
 \tag{1}
\]
Here \(b_j=\Lambda(g^j)/j!\), and the formula follows from
\(\Lambda(Y_i^j)=j!\). The forward adjunction identity is
used by [Rocks, Theorem 3.12 and Corollary 3.13](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf).
A reverse lift by one independent linear term is also in the
[restricted quadratic adjunction theorem](QUADRATIC_STRONG_FACTORIAL_LINEAR_ADJUNCTION_BRIDGE.md).
The symmetric formulation gives the full reverse and the exact
corresponding ideal.

In the universal ring of phase coefficients and \(h_1,\ldots,h_k\)
the ideal equality is
\[
 (b_1,\ldots,b_{k+s})
   =(a_1-h_1,\ldots,a_k-h_k,\ a_{k+1},\ldots,a_{k+s}).
 \tag{2}
\]
It follows coefficientwise from (1): for \(j\le k\),
\(a_j-h_j=b_j+\sum_{i=1}^{j-1}h_i b_{j-i}\), and for \(j>k\)
the corresponding leading term is \(a_j=b_j+\sum_{i=1}^k h_i b_{j-i}\).
These are triangular with diagonal one, which proves both ideal
inclusions over any coefficient \(\mathbb Q\)-algebra.
Eliminating \(h_1,\ldots,h_k\) therefore recovers the shifted
moment ideal, and localizing at \(h_k\) corresponds to localizing
at \(a_k\).

If \(a_k\ne0\), choose \(h(t)=\sum_{i=0}^k a_i t^i\) and factor
over an algebraic closure. Its constant term one and exact
degree \(k\) give \(k\) nonzero \(d_i\), including multiplicities.
Equation (2) gives the equivalence and uniqueness of their multiset.
For a window starting at \(n\) with \(a_{n-1}\) possibly zero,
take \(h(t)=\sum_{i=0}^{n-1}a_i t^i\) and use its actual degree
\(k'\le n-1\). Then \(B-1=(A-h)/h\) vanishes through \(n+s-1\),
which is at least the actual support size \(k'+s\).

The equivalence uses arbitrary finite dimension because the added
linear terms require independent new variables. It holds over any
one algebraically closed characteristic-zero field and for the
universal assertion quantified over all characteristic-zero fields.
For the ordinary binomial problem the
[absolute eventual cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
bounds the number of linear auxiliaries needed, although the
ordinary exponent profiles remain unrestricted.
