# Dense two-letter graphs with sparse repaired additive relations

Two complete radical caps and dense incidence between their letters
need not force additive structure in the associated least-repaired
values. This remains true for two completely diffuse valuation
profiles. In the example, almost every value of the third letter
has outside radical close to the product of the two allowed
radical budgets. A general entropy calculation identifies the
same obstruction when the third radical test is replaced by a
bound on the number of allowed labels.

## Two diffuse complete caps

Fix \(K>4\) and \(B_0>0\). There is a sequence of fixed budgets
\(\lambda_m\downarrow0\) such that, at each fixed \(m\) and all
sufficiently large heights \(X\), there are distinct-value
dictionaries \(\mathcal A,\mathcal C\) with
\[
 \begin{gathered}
 n=X^{\lambda_m},\qquad
 |\mathcal A|=|\mathcal C|=n^{1/3+o_X(1)},\\
 \mathcal A\subset[X/8,X/4],\qquad
 \mathcal C\subset[3X/4,X],\qquad
 \operatorname{rad}(a),\operatorname{rad}(c)\le n^{1/3}.
 \end{gathered}
\]
Every \((a,c)\in\mathcal A\times\mathcal C\) gives a primitive row
\(a+b=c\). Both selected letters have full scaled valuation
measures tending to \((1/3)\operatorname{Unif}[2,4]\), where
\(\operatorname{Unif}\) is a probability measure. The height
limit is taken at each fixed budget before the outer limit
\(m\to\infty\).

Given any prescribed set \(S\) of primes of product at most
\(X^{B_0}\), write \(\operatorname{rad}_{\overline S}\) for the
radical restricted to primes outside \(S\). Uniformly in \(S\),
all but an \(O_{m,B_0}(n^{-1/m})\) fraction of the rows have
\[
              \operatorname{rad}_{\overline S}(b)
                         > n^{\,2/3-4/m}.
 \tag{TA1}
\]
The complete dictionary of \(b\) has size at least
\(n^{2/3-1/m-o_X(1)}\). For
\(M=\lceil K/\lambda_m\rceil\), all least-repaired values of
\(a,c\) are exact \(M\)-th powers. If \(\mathcal D\) is their
united dictionary, then for each fixed \(h\ge2\) and
\(\eta,\varepsilon>0\), for sufficiently large \(m\) the
number of equal sums of two ordered \(h\)-tuples of
\(\mathcal D\) off the permutation diagonal is at most
\[
 O_{h,K,\eta,\lambda_m,\varepsilon}
     \bigl(|\mathcal D|^{h-1}n^\eta X^\varepsilon\bigr).
 \tag{TA2}
\]
Thus even a full Cartesian projection on two tall diffuse
letters with both radical caps can coexist with nearly
diagonal energy for their repaired heights.

### Construction and the third radical

For each \(m\ge7\), choose an integer \(t_m\) divisible by
\(m-1\), with \(t_m\) strictly increasing and sufficiently large
as specified below. Set \(\lambda_m=1/t_m\) and
\[
   e_j=2t_m+\frac{2t_mj}{m-1}\qquad(0\le j<m),
   \qquad E_0=\sum_j e_j=3mt_m,\qquad
   U=X^{1/E_0}.
\]
Take \(m\) disjoint prime pools in
\([8^{-1/E_0}U,4^{-1/E_0}U]\) for the first letter and \(m\)
in \([(3/4)^{1/E_0}U,U]\) for the last letter. The intervals
are disjoint. By the prime number theorem, all \(2m\) pools
can have a common size \(n^{1/(3m)+o_X(1)}\). These are fixed proportional
intervals in \(U\) at each fixed budget. A letter uses one
prime per pool, to the power \(e_j\) in slot \(j\).
The pools and unique factorization recover every chosen
root tuple from its complete letter. Their products have
the specified heights. The prime support of a letter has
radical at most \(U^m=n^{1/3}\), and the two roles use
disjoint primes, so every row is primitive. The individual
root masses tend to \(1/(3m)\) as \(X\to\infty\), and the
scaled degrees \(e_j/t_m\) form an equally spaced grid
on \([2,4]\).

At a fixed \(b>0\), the root tuples satisfy the nondegenerate
relation
\[
             \prod_j q_j^{e_j}-\prod_j p_j^{e_j}-b=0.
\]
There are \(2m\) pairwise coprime variables of size at most
\(U=n^{1/(3m)}\), and the exponent spread in each variable
is at least \(2t_m\). The bound (EB1) is uniform in \(b\).
Choose \(t_m\) so large that
\[
            \frac{C_{2m,3}}{3m}
             (2t_m)^{-1/(2(2m-1))}\le\frac1m.
\]
It gives \(O_{m,\lambda_m,\varepsilon}(n^{1/m}X^\varepsilon)\)
rows for each \(b\), implying the stated size of the third
dictionary. In particular the budget can be chosen before
\(X\) and the growing number of degrees causes no
interchange of limits.
[Nondegenerate monomial count (EB1)](ABC_NONDEGENERATE_POWER_COLLISIONS.md#exclusive-boundary-release-for-nondegenerate-power-collisions)

For any \(n^z\) cap on the outside radical of \(b\), number its
outside support products first. After adjoining the prescribed
pool, each support gives \(X^{o_X(1)}\) complete values \(b\le X\),
uniformly over its polynomial-height product. Thus the count
under the cap is \(O(n^{z+1/m}X^\varepsilon)\), uniformly in
the prescribed \(S\) and \(z\ge0\).
At \(z=2/3-4/m\), the ratio to the full family is at most
\(n^{-1/m}\) for sufficiently large \(X\) at fixed \(m\),
choosing the arbitrary height-exponent allowance smaller
than \(\lambda_m/m\). This proves (TA1).
[Fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

All \(e_j\) are at most \(4t_m<M\), so the least-repaired
heights are \(\operatorname{rad}(a)^M\) and
\(\operatorname{rad}(c)^M\). In the shifted-sums and energy
proof (HP3)–(HP4), the tails for these values are identically
one. The same proof therefore gives (TA2) with no tail
charge. The number of membership patterns for an energy
of fixed order \(h\) is bounded in terms of \(h\),
independently of \(m\). This allows the remaining
release loss to be made smaller than any fixed \(\eta\)
as \(m\to\infty\).
[Shifted sums and true additive diagonals](ABC_REPAIRED_HEIGHT_POLYNOMIAL_RELATIONS.md#shifted-sums-and-true-additive-diagonals)

## Arbitrary fixed orders on the two capped letters

A separate version allows any fixed integers \(p,q\ge2\).
Use all pairs of prime roots \(x,y\) with
\(x^p\in[X/8,X/4]\) and \(y^q\in[3X/4,X]\).
The prime intervals are disjoint when \(p=q\), and their
scales separate for large \(X\) when \(p\ne q\).
The resulting primitive rows have
\(\operatorname{rad}(x^p)\le X^{1/p}\) and
\(\operatorname{rad}(y^q)\le X^{1/q}\). There are
\(X^{1/p+1/q+o(1)}\) of them. The two repaired dictionaries
are exact powers of any common order \(M\ge\max(p,q)\).

For fixed \(b>0\), the curve \(y^q=x^p+b\) is absolutely
irreducible by Eisenstein at an arbitrary simple root of
\(x^p+b\) over \(\mathbb C\). Its degree is \(\max(p,q)\)
and its root coordinates are bounded by
\(X^{1/\min(p,q)}\). The coefficient-uniform
Bombieri–Pila bound gives \(O_{p,q,\varepsilon}
(X^{1/(pq)+\varepsilon})\) pairs at this level. Thus the
number of rows with \(\operatorname{rad}_{\overline S}(b)
\le X^z\) is
\[
                   O_{p,q,B_0,\varepsilon}
                         (X^{z+1/(pq)+\varepsilon})
 \tag{TA3}
\]
for any prescribed polynomial-height prime pool.
In particular \(p=m,q=m+1\) are coprime orders with
\(\lambda=3/m\). Their two normalized radical caps
approach \(1/3\), whereas the generic lower threshold
for the normalized third outside radical approaches
\(2/3\), since \(1/p+1/q-1/(pq)=2/(m+1)\).
[Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)
· [Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)

## A ceiling for third-dictionary cardinality in entropy bounds

Take any finite family of distinct rows and its uniform
law \(\mathbb P\) on \((A,C)\). The row is recovered from
this pair, so \(H(\mathbb P)\) is the logarithm of its row
count. Let \(\mathbb Q=\mathbb P_A\otimes\mathbb P_C\)
and \(\mathcal B\) be a dictionary containing every
\(C-A\) on the family. On the event
\(E=\{C-A\in\mathcal B\}\), \(\mathbb P(E)=1\).
The log-sum inequality, or conditioning \(\mathbb Q\) on
\(E\), gives
\[
 I_{\mathbb P}(A;C)=D(\mathbb P\Vert\mathbb Q)
                    \ge-\log\mathbb Q(E).
\]
If \(H(A)+H(C)\le (d_A+d_C)\log n\),
\(|\mathcal B|\le n^\beta\), and the independent
product-measure convolution has pointwise bound
\(\mathbb Q(C-A=z)\le n^{-d_A-d_C+\kappa}\), this argument gives
only
\[
                          H(\mathbb P)\le(\beta+\kappa)\log n.
 \tag{TA4}
\]
The diffuse construction realizes the ceiling at the critical
exponent on actual primitive rows with the two complete
radical caps. Let the two sides have size \(D=n^{1/3+o_X(1)}\).
Delete the pairs failing (TA1) for a fixed prescribed pool.
Now greedily choose edges with distinct \(a,c,b\), one at a
time. If \(u\) vertices remain on each side, at most
\(O(D^2n^{-1/m})\) pairs are excluded by the initial
deletion. At most \(O(Dn^{2/m})\) further pairs use an
already chosen value of \(b\), by the fixed-\(b\) bound
with height-exponent allowance at most \(\lambda_m/m\).
The process continues while \(u^2\) exceeds the sum of
those bounds. For fixed \(m\ge7\) this leaves \(o_X(D)\)
vertices, so the resulting matching has size
\(D'=(1-o_X(1))D\).

The restriction has exactly \(D'\) values on every letter
and retains the valuation profiles and (TA2) on the two
capped roles. All its third values satisfy (TA1).
For its uniform row law, \(I(A;C)=\log D'\) and the
dictionary of its used differences has independent
product mass at least \(1/D'\), from the matching edges.
Each difference has product mass at most
\(O(n^{1/m}X^\varepsilon/(D')^2)\), by the fixed-\(b\)
bound on the original graph. Consequently as \(m\to\infty\)
the matching attains the entropy ceiling in normalized
exponents with \(\beta\to1/3\) and \(\kappa\to0\),
even though its third radical mass exceeds that budget.

At \(\beta=1/3\) even the ideal pointwise exponent
\(\kappa=0\) in (TA4) reaches only \(1/3\). Improving
it through this argument requires a further reduction
of product mass on the actual radical-constrained
dictionary, or an independent entropy constraint.
The diffuse construction above demonstrates why the
corresponding third-letter hypothesis cannot be removed
when transferring original additive relations to
least-repaired heights.
