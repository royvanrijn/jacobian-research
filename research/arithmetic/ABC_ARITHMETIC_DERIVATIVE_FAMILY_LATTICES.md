# Simultaneous arithmetic derivatives on primitive abc families

A collection of primitive rows imposes a joint lattice condition on its actual prime valuations and its additive equations. This gives a way to combine a common powered modulus with valuation ranks. The estimate is uniform in the number of rows, and can be applied directly to the primitive quotients of least-repaired triples.

## The family inequality

Let \(\mathcal F\) be a finite family of positive primitive triples \(a_i+b_i=c_i\) with \(X/2<c_i\le X\). Choose any subset \(S\) of the primes occurring in the family. For a positive integer \(T\) write \(V_T=(v_p(T))_{p\in S}\), and set
\[
 L_T(\boldsymbol z)=\sum_{p\in S} v_p(T)z_p/p .
\]
Define the rational subspaces
\[
 K=\{\boldsymbol z:a_iL_{a_i}+b_iL_{b_i}-c_iL_{c_i}=0
                         \text{ on every row}\},
 \qquad
 Z=\{\boldsymbol z:L_{a_i}=L_{b_i}=L_{c_i}
                         \text{ on every row}\}.
\]
The \(L_T\) in these equalities are evaluated at \(\boldsymbol z\). Let \(A\) be a matrix of a maximal independent subset of the original additive rows, with integer entries \(a_i v_p(a_i)/p+b_i v_p(b_i)/p-c_i v_p(c_i)/p\). Take \(r=\operatorname{rank}A\), \(u=\dim(K/Z)\), and let \(\Delta(A)\) be the positive gcd of the \(r\)-minors of \(A\). At \(r=0\) the determinant and the gcd below are both one.

Put \(R=\max_i\operatorname{rad}(a_i b_i c_i)\) and \(H=X\log2/(2R\log X)\). Assume \(H\ge1\). A prime \(p\in S\) is eligible if \(p\le H\) and some occurrence of \(p\) in the family has valuation not divisible by \(p\). If \(I\) is any set of eligible primes whose coordinate functionals are linearly independent on \(Z\), then
\[
 H^u\prod_{p\in I}p
      \ \le\ \frac{\sqrt{\det(AA^{\mathsf T})}}{\Delta(A)}
      \ \le\ \frac{\prod_{i\ {\rm in}\ A}\|A_i\|_2}{\Delta(A)}.
 \tag{AD1}
\]
In particular, for any fixed \(0<\lambda<1/2\) and \(R\le X^\lambda\), all occurring primes greater than \(\log_2 X\) are eligible for sufficiently large \(X\). Each row norm is at most \(X\log X/\log2\). All constants in these displays are independent of the family size.

## Proof of the size and divisibility bounds

For integer \(\boldsymbol z\), the map \(D(T)=TL_T(\boldsymbol z)\) has integer values and obeys the Leibniz identity on positive integers. It is additive on each row when \(\boldsymbol z\in K\). If \(\boldsymbol z\in K\setminus Z\), take a row on which the three logarithmic derivatives are unequal. Its Wronskian
\[
 W=aD(b)-bD(a)=aD(c)-cD(a)=cD(b)-bD(c)
\]
is nonzero, by additivity. For each of \(T=a,b,c\), \(T/\operatorname{rad}(T)\) divides \(D(T)\) and \(T\). The three quantities \(T/\operatorname{rad}(T)\) are coprime, so \(abc/\operatorname{rad}(abc)\) divides \(W\). On the other hand,
\[
 |W|=ab|L_b-L_a|
        \le ab\,\|\boldsymbol z\|_\infty
             \sum_{p\mid ab}v_p(ab)/p
        \le ab\,\|\boldsymbol z\|_\infty\,\log X/\log2 .
\]
It follows that \(\|\boldsymbol z\|_\infty\ge H\). This uses the primitive letters even if the original input is a repaired multiple.

If \(\boldsymbol z\in Z\cap\mathbb Z^S\), the three \(L_T\) on any row have pairwise coprime denominator supports. Their common value is therefore integral. If \(p\in S\) occurs in one letter at a valuation coprime to \(p\), the term \(v_p(T)z_p/p\) is integral at \(p\), since all other summands of \(L_T\) are integral there. Thus \(p\mid z_p\).

Apply the product form of the Bombieri–Vaaler Siegel lemma to \(A\). It gives a rational basis of \(K\) consisting of integer vectors with the product of their supremum norms at most \(\sqrt{\det(AA^{\mathsf T})}/\Delta(A)\). At least \(u\) of the vectors are outside \(Z\), so their norms are at least \(H\). Order \(I=\{p_1>\cdots>p_k\}\), and inside \(Z\) impose \(z_{p_1}=\cdots=z_{p_j}=0\). The resulting subspace has codimension \(u+j\) in \(K\), so at least \(u+j\) basis vectors lie outside it. Each of these has norm at least \(p_j\), either by the \(H\) bound or by a nonzero coordinate divisible by one of \(p_1,\ldots,p_j\). Ordering the basis norms from largest to smallest now gives the first bound of (AD1). Hadamard gives the second. For the stated row bound, split the coordinates by letter, which makes the three subvectors disjoint. The sum of the absolute coefficients on each letter is at most \(X\log X/(2\log2)\), since \(p^{-1}\le1/2\) and \(\sum_{p\mid T}v_p(T)\le\log X/\log2\). The Euclidean norm of the row is at most \(\sqrt3\) times this bound. [Pastén, Theorems 2.5 and 3.3](https://arxiv.org/pdf/2106.16165v3)
· [Bombieri and Vaaler, Theorem 2](https://link.springer.com/article/10.1007/BF01393823)

## Minor divisibility from all three letters

For a prime \(p\), an integer \(m\ge1\), and the particular rows chosen in \(A\), let \(J_{p,m}\) be the rows with \(p^m\) dividing one letter. In each of these rows form the difference of the two other valuation vectors \(V_T\), in either order. Let \(s_{p,m}\) be the rational rank of these difference vectors, and put
\[
 \ell_p=\sum_{m\ge1}\left(|J_{p,m}|-s_{p,m}\right).
\]
Then
\[
                    v_p(\Delta(A))\ge(\ell_p-\mathbf 1_{p\in S})_+ .
 \tag{AD2}
\]
The omitted letter can vary by row and prime. The inequality does not require that the prime be a selected coordinate.

To prove it, multiply each column indexed by \(q\in S\) by \(q\), and call the resulting matrix \(A'\). All maximal minors are multiplied by the product of the column labels they select. Hence \(v_p(\Delta(A'))\le v_p(\Delta(A))+\mathbf 1_{p\in S}\). Modulo \(p^m\), a row of \(A'\) indexed by \(J_{p,m}\) is a unit multiple of its stated difference vector. For example, if \(p^m\mid c_i\), then \(a_i\) is a unit and \(b_i\equiv-a_i\pmod{p^m}\), so the row is \(a_i(V_{a_i}-V_{b_i})\) modulo \(p^m\). The other two cases follow by replacing \((a_i,b_i)\) by the two \(p\)-free letters and using the equation with its signs.

The rational rank of those difference vectors is also their rank over \(\mathbb Q_p\). Add the \(r-|J_{p,m}|\) remaining rows of \(A'\) to their span and intersect with \(\mathbb Z_p^S\). This gives a saturated direct summand of rank at most \(s_{p,m}+r-|J_{p,m}|\). The image generated by the rows of \(A'\) modulo \(p^m\) is contained in the reduction of that summand. A submodule of a free \(\mathbb Z/p^m\mathbb Z\)-module of rank \(d\) is generated by at most \(d\) elements. By Smith normal form the number of invariant factors of \(A'\) with \(p\)-adic order less than \(m\) is therefore at most \(s_{p,m}+r-|J_{p,m}|\). Summing the complementary inequality over all \(m\) bounds the sum of the orders of the invariant factors, which is \(v_p(\Delta(A'))\). This proves (AD2).

The loss of one \(p\)-adic digit in (AD2) is sharp on actual rows. Take distinct primes \(p,q>3\) and a positive integer \(m\) with \(p\nmid m\) and \(q<2p^m\). For \(i=2,3\) use the primitive row
\[
              q+(ip^m-q)=ip^m
\]
and the coordinate set \(S=\{p,q\}\). The two normalized derivative rows are \((-im p^{m-1},1)\), so the unique maximal minor has absolute value \(m p^{m-1}\). The two rows lie in the shell \(X=3p^m\). At every depth up to \(m\), the pair of valuation differences on the two \(p\)-free letters has rank one, giving \(\ell_p=m\) and \(v_p(\Delta(A))=m-1\) exactly. This verifies sharpness of the normalization for unbounded shared powered moduli. It does not impose a cap on the radical of \(ip^m-q\).

For a basis row \(i\), let \(p_{\min,i}\) be the smallest selected prime present on that row. The row is nonzero, so this is defined. Its Euclidean norm is at most \(3X\log X/(p_{\min,i}\log2)\), by summing the absolute coefficients. Combining (AD1) and (AD2) gives the following normalized inequality for any fixed \(0<\lambda<1/2\) with \(R\le X^\lambda\) and \(X\) sufficiently large:
\[
 (1-\lambda)u
 +\sum_{p\in I}\log_X p
 +\sum_{i\ {\rm in}\ A}\log_X p_{\min,i}
 +\sum_p(\ell_p-\mathbf 1_{p\in S})_+\log_X p
 \ \le\ r+O\!\left((r+u)\frac{\log\log X}{\log X}\right).
 \tag{AD2'}
\]
The implied constant is absolute. Here the height limit can be taken at any fixed budget in the stated range and the family size is unrestricted.

For example, suppose fixed pairwise coprime \(D_A,D_B,D_C\) divide the corresponding complete letters throughout the family. For a color \(T\), let \(s_T\) be the rank of the valuation differences on the other two letters, using the rows selected in \(A\). Then (AD2) shows that \(\Delta(A)\) is divisible by
\[
 \frac{\prod_T D_T^{\,r-s_T}}
      {\displaystyle\prod_{T:r>s_T}\ \prod_{\substack{p\in S\\p\mid D_T}}p}.
 \tag{AD3}
\]
The denominator has each selected prime at most once by coprimality of the \(D_T\). Formula (AD2) can be stronger because it uses all primes and orders present in the rows.

Let \(Z_A\) be the equal-derivative space defined using just the rows selected in \(A\), and let \(u_A=\dim(K/Z_A)\). On \(K\) the other two logarithmic derivatives on a basis row are equal if and only if all three are equal. Thus \(u_A\le s_T\) for each color. The full family gives the bound \(u_A\le u\). For any two distinct colors \(T,U\), their two difference systems together define \(Z_A\) over the full space \(\mathbb Q^S\), so \(s_T+s_U\ge r+u_A\).

If the family contains two distinct ordered primitive rows, \(D_A D_B D_C\le X^2\). Indeed all three moduli divide their nonzero projective minor \(a_i c_j-a_j c_i\), using the equivalent expression \(b_j c_i-b_i c_j\) for \(D_B\). Primitivity implies that equality of ratios \(a_i/c_i=a_j/c_j\) gives identical rows. Also \(D_T\le X\) individually. Put \(d_T=r-s_T\), \(h_T=\log_X D_T\), \(K_0=r-u_A\). We have \(0\le h_T\le1\), \(\sum_T h_T\le2\), and \(d_T+d_U\le K_0\) for distinct colors. These inequalities imply
\[
                   \sum_T h_T d_T\le K_0 .
 \tag{AD4}
\]
If all \(d_T\le K_0/2\), use \(\sum h_T\le2\). Otherwise the other two \(d_U\) are at most \(K_0-d_T\), and
\(h_Td_T+(2-h_T)(K_0-d_T)\le K_0\).
Thus the three common moduli alone can certify at most \(X^{r-u_A}\) in the numerator of (AD3). Any use of a larger normalized minor divisor needs the individual depth sets in (AD2) or further information. In (AD1), the prime coordinate flags, the Gram determinant and the extra equality constraints on redundant additive rows can all contribute separately.

## Sparse supports and the remaining input

Suppose that the selected prime sets on distinct rows are disjoint, with each of the three letters on each row still represented. For a family with \(k\) rows one has \(r=u=k\). In fact the additive space has codimension one per row and the equal-logarithmic-derivative subspace has codimension two per row, because the three linear forms have disjoint nonempty supports. Both dimensions add across rows. The Gram determinant and the maximal-minor gcd split:
\[
 \sqrt{\det(AA^{\mathsf T})}=\prod_i\|A_i\|_2,\qquad
 \Delta(A)=\prod_i\gcd\{A_{i,p}:p\in S,\ A_{i,p}\ne0\}.
\]
Every nonzero maximal minor uses one coordinate from each row, and all such selections occur. Each row gcd is at most \((\log_2 X)^2\), by taking one represented prime in each of two letters. The corresponding coefficients have the form \((T/p)v_p(T)\), and the two factors \(T/p\) are coprime, so their gcd divides the product of the two valuation exponents. Thus the maximal-minor normalization in (AD1) is at most a power of \(\log X\) per row in this situation.

If all primes outside a fixed discarded pool occur in at most one row, (AD2) gives no gain at those primes, as \(J_{p,m}\) is a singleton with a nonzero difference vector whenever it is nonempty. In the example where each row has exactly one prime per letter selected, the possible independent sets \(I\) select at most one prime per row. More generally, if every selected prime is eligible and all three letter supports are nonempty, a full coordinate basis on the equal-derivative space omits exactly two selected primes per row, and the omitted primes must be from distinct letters. This follows by checking the two equality equations: the two columns of omitted coordinates must be independent, and their two-vectors are proportional exactly for positions on the same letter. The entire inequality is then bounded by products of per-row conditions.

Actual primitive rows can realize this sparse pattern outside a fixed prime pool. Fix distinct positive integers \(h_1,\ldots,h_k\) and take
\[
        a_i=t+h_i,\qquad b_i=2(t+h_i)+1,\qquad
        c_i=3(t+h_i)+1 .
\]
They satisfy additivity and \(\gcd(a_i,b_i)=1\). The \(3k\) linear polynomials have distinct rational roots. Put every prime divisor of their pairwise resultants into the discarded pool, and select all remaining primes. Each selected prime occurs in at most one of the \(3k\) letter values, since a common prime would divide the corresponding resultant. All three letters are represented for all but \(O((1+\log T)^{O_k(1)})\) positive parameters \(t\le T\): a fixed pool allows at most a product of logarithmic factors many integers in \([1,O_k(T)]\), and each linear polynomial is injective in \(t\). The rows lie in one dyadic shell for sufficiently large \(t\). This example concerns sparsity and rank with actual additive equations. It has no asserted small radical budget.

The known single-triple Wronskian argument supplies the lower bound for every derivative outside \(Z\). The exact single-row Wronskian image was calculated by Laniewski–Müller. The additional conditions here use the occurrence of the same prime coordinate on distinct actual rows and the \(p\)-adic rank of the other two valuations. A counting application still needs a rank or minor-divisibility estimate for an extracted subfamily with all three radical budgets. The quotient rank \(u\) can be smaller than \(r\). The bound (AD4) describes the capacity of fixed simultaneous common moduli in terms of the equality rank on basis rows. The family inequality supplies a necessary condition on such an extraction, with no inferred unrestricted exponent saving. [Laniewski–Müller, Section 4](https://arxiv.org/html/2610.26490v1#S4)
