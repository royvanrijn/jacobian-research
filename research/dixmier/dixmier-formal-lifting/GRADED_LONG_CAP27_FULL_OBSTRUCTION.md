# Full arbitrary-parity obstruction for the canonically graded Long tuple at spatial cap 27

**Theorem GradedLongCap27FullExclusion.** Let \(k\) be any characteristic-zero field and \(f=(R,T,D,S)\) the canonical global polynomial Long tuple. Use scalar Moyal symbols in the variables \((X,Q,P,Z)\) with \(\{P,X\}=\{Z,Q\}=1\) and weights
\[
\operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
\operatorname{wt}(R,T,D,S)=(-1,2,4,1).
\]
There is no weight-homogeneous lift of \(f\) of common total spatial degree at most \(27\) satisfying the six Weyl commutators through \(h^7\). Odd and even coefficients in all four images are admitted. Thus the complete graded scalar degree-\(27\) scheme has no field-valued points, and a full lift in this grading requires common spatial degree at least \(28\).

The input tuple is the literal [canonical Long coefficient expansion](../dixmier-audit/verify_long_pbw.py). The [cap-\(26\) result](GRADED_LONG_CAP26_FULL_JET_OBSTRUCTION.md) and the [cap-\(27\) even computation](GRADED_LONG_CAP27_EVEN_JET_OBSTRUCTION.md) also establish that \(27\) is the minimum for jets through \(h^5\). The theorem concerns the displayed grading and fixed classical tuple. The nongraded problem, degree \(28\) and higher, nonscalar constructions, and the rank-two Dixmier conjecture remain open. A positive point of the complete scheme would meet the common-bound requirement in [bounded realization](BOUNDED_REALIZATION.md).

## The complete coefficient scheme

Write \(F_i=\sum_m h^m F_{m,i}\), with \(F_{0,i}=f_i\), and let \(w_i\) be the weight of \(f_i\). The full basis at order \(m\) consists of \(X^aQ^bP^cZ^d\) with
\[
a=b+4c+2d+3m-w_i\geq0,\qquad
2b+5c+3d+3m\leq27+w_i.
\]
In image order \((R,T,D,S)\) the largest possible \(m\) are \(8,9,10,9\), respectively. The table counts every positive-order coefficient.

| \(m\) | \(R\) | \(T\) | \(D\) | \(S\) | Total | Rational rank | Nullity |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 127 | 172 | 206 | 156 | 661 | 647 | 14 |
| 2 | 91 | 127 | 156 | 114 | 488 | 480 | 8 |
| 3 | 62 | 91 | 114 | 80 | 347 | 345 | 2 |
| 4 | 40 | 62 | 80 | 54 | 236 | 235 | 1 |
| 5 | 24 | 40 | 54 | 34 | 152 | 152 | 0 |
| 6 | 13 | 24 | 34 | 20 | 91 | 91 | 0 |
| 7 | 6 | 13 | 20 | 10 | 49 | 49 | 0 |
| 8 | 2 | 6 | 10 | 4 | 22 | 22 | 0 |
| 9 | 0 | 2 | 4 | 1 | 7 | 7 | 0 |
| 10 | 0 | 0 | 1 | 0 | 1 | 1 | 0 |

The total is \(2054\). For every \(i<j\), the commutator coefficient at \(h^N\) is
\[
\sum_{\substack{a+b+\ell=N\\\ell\geq1\ {\rm odd}}}
\frac{\Pi^\ell(F_{a,i},F_{b,j})}{2^{\ell-1}\ell!},
\]
where \(\Pi^\ell\) is the \(\ell\)-fold contraction with the constant Poisson tensor. A nonzero \((P,X)\) contraction consumes five units from the combined \(2b+5c+3d\) budget of the two monomials. A \((Z,Q)\) contraction consumes \(3+2=5\) units. Consequently \(3(a+b)+5\ell\leq54+w_i+w_j\leq60\), so \(a+b+\ell\leq19\). This bound is sharp before imposing the equations. The order-eight monomial \(X^{26}P\) in \(T\) and the order-ten monomial \(X^{26}\) in \(D\) have nonzero Poisson bracket. The coefficient polynomial at \(h^{19}\) contains the product of their two independent coefficients. The full scheme imposes the commutators through \(h^{19}\).

## Sequential rational elimination

The checker verifies the classical brackets at \(h\). At each subsequent order \(h^{m+1}\) the new tuple \(U=F_m\) occurs only through
\[
(A_mU)_{ij}=\{U_i,f_j\}+\{f_i,U_j\}.
\]
The right side is the negative of the known terms in the contraction formula. The table gives the ranks of these exact rational maps. Reduction modulo \(1000003\) selects pivot rows and columns. The checker inverts that square submatrix over \(\mathbb Q\), then verifies the complete kernel by multiplying each basis vector against the full matrix and checking its free coordinates. Solving the pivot equations and projecting away their contribution gives a residual in the full row space. A tuple extends to this order exactly when all residual entries vanish. When they do, the free kernel coordinates give every extension. Taking a rational row basis of the residual polynomials preserves their common zeros.

There are \(25\) free kernel coordinates across all ten orders. The [cap-preserving source and target shears](GRADED_LONG_CAP26_FULL_JET_OBSTRUCTION.md) remove the \(XQ\) and \(RS\) Hamiltonian directions at order one. Explicitly, these are \(P\mapsto P-\lambda hQ,\ Z\mapsto Z-\lambda hX\) on the source and \(T\mapsto T+\mu hR,\ D\mapsto D-\mu hS\) on the target. Substitution preserves the cap because these transformations are linear in spatial variables. They preserve the weights and commutators. The checker completes the resulting nine old directions by three new kernel directions. It uses every free coordinate at the subsequent orders.

## Exact elimination of the residual equations

The following equations denote exact rational linear combinations of residual polynomials at the indicated order, after the preceding substitutions and reductions. For polynomial ideal reductions the replay computes a Gröbner basis over \(\mathbb Q\), verifies a transformation expressing each basis element in the original generators, and verifies that each original generator reduces to zero. The substitutions deduced from powers vanishing hold in every field of characteristic zero.

**Orders \(h^2\) and \(h^3\).** After the two shears the first-order kernel coordinates are \(a,b,c,d,e_{\rm old},f,\eta,j,k,x_0,x_1,x_2\). The \(h^3\) residual contains \(x_0^2,x_1^2,x_2^2\). With these coordinates zero, the next equations include
\[
 (2j+k)^2=0,\qquad d^2-250df+35625f^2=0.
\]
Write \(j=3z,\ k=-6z,\ e_{\rm old}=e+2z\). The checker verifies that the full \(h^3\) residual vanishes modulo the displayed conic. It uses an order-two particular together with the full eight-dimensional homogeneous space, with coordinates \(u,v,\gamma,\delta,\chi,\psi,\kappa,\sigma\).

**Order \(h^4\).** Modulo the conic the residual contains \(f^2(d+375f)\) and \(e^2z\). If \(f\ne0\), substitution of \(d=-375f\) into the conic gives \(270000f^2=0\). Hence \(d=f=0\) for every point, and \(ez=0\). All other \(h^4\) equations under these substitutions are kept. The two free order-three coordinates are \(\omega,\tau\). The \(h^4\) ideal contains \(z(\kappa-162z^2)\).

**Order \(h^5\).** The residual modulo the preceding ideal contains
\[
 \sigma^2,\qquad \kappa(\kappa-162z^2),\qquad
 -27\delta\kappa+4374\delta z^2+11\kappa^2,\qquad e^3\eta .
\]
Combining the preceding \(z(\kappa-162z^2)=0\) with the second equation gives \((\kappa-162z^2)^2=0\). The third equation then gives \(\kappa=z=0\). Thus \(\sigma=\kappa=z=0\) and \(e\eta=0\). The full residual also supplies
\[
 3\gamma+10\chi+9\psi-576c^2+3888ce=0.
\]
Eliminate \(\gamma\) by this relation. Keep the entire one-dimensional order-four kernel, denoted by \(w\). There are now thirteen coordinates
\[
 a,b,c,e,\eta,u,v,\delta,\chi,\psi,\omega,\tau,w.
\]

**Order \(h^6\).** The map \(A_5\) is injective. Use the \(h^5\) ideal to reduce its right side, take the unique pivot solution, and add all residual equations to the lower-order equations. The exact basis of the resulting ideal, up to nonzero rational factors, is
\[
\begin{gathered}
e,\ c+36\eta,\ 3a+2176\eta,\ \psi\tau,\ \chi\tau,\ \delta\tau,\
\psi^2-3888\eta\tau,\ \chi\psi,\ \eta\psi,\ \eta\tau^2,\ \eta^2\tau,\\
\delta(5\chi-10077696\eta^2),\qquad
\eta^2(5\chi-10077696\eta^2).
\end{gathered}
\]

**Order \(h^7\).** Reduce the full right side for \(A_6\) by this exact ideal and use all 91 coefficient columns. The rational row span of the remaining residual polynomials contains \(1\). Thus no point extends through \(h^7\). Every solution of the complete scheme would restrict to such a point, so the higher commutator coefficients are unnecessary for the exclusion.

## Replay

Run the [exact checker](verify_graded_long_cap27_full.py) with Python, python-flint, and Singular. For example, from the repository root use:

~~~sh
python research/dixmier/dixmier-formal-lifting/verify_graded_long_cap27_full.py
~~~

The [verification receipt](graded_long_cap27_full_verification.json) records the input hashes, coefficient counts, ranks, reduced \(h^6\) basis, and the \(h^7\) residual test. The checker reconstructs the matrices and nonlinear defects. It also compares the symbolic Moyal arithmetic with the retained rational implementation at a fixed rational specialization at orders \(h^3\) and \(h^7\). Assurance is internal and computational, with exact rational arithmetic.
