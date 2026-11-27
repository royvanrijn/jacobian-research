# Primitive 4-full triples in Erdős problem 939

The case \(r=4\) asks whether \(a+b=c\) is possible for positive 4-full integers with \(\gcd(a,b)=1\). An integer is 4-full if every prime in its factorization has exponent at least \(4\), so \(1\) is included. The equation makes the three terms pairwise coprime. One example would settle this existence question. The problem page records examples at \(r=5\) and an infinite construction for each \(r\ge6\), where coprimality is joint. The \(r=4\) existence question remains open. [Erdős, p. 33](https://www.renyi.hu/~p_erdos/1976-39.pdf) · [Erdős Problems #939](https://www.erdosproblems.com/939)

The Campana point is \((a:c)\) on \((\mathbb P^1,\frac34([0]+[1]+[\infty]))\) and has height \(c\). The [powerful-number manuscript](../../unpublished/powerful-number-representations/main.tex) gives the bound \(N_4(X)\ll_\varepsilon X^{2/5+\varepsilon}\). Its canonical factors and rational reconstruction also give an explicit search, described below. [Browning–Verzobio, §1](https://arxiv.org/html/2608.24512v1#S1)

## Exhaustive search in a finite height range

For a maximum-summand scan through \(B\), generate each 4-full value through \(2B\) by recursively choosing increasing primes and exponents at least four. The prime bound \(\lfloor(2B)^{1/4}\rfloor\) follows from \(p^4\mid n\). The sorted prefix through \(B\) is the summand list. The checker buckets this prefix by residue modulo \(16\cdot81\) and visits each eligible unordered pair once. For \(p=2,3\), a pair is excluded if both summands are divisible by \(p\), or if their sum is divisible by \(p\) but fails divisibility by \(p^4\). An absent target residue or unset target bit also certifies absence. Every remaining pair is tested for exact target-set membership and then for coprimality.

At \(B=10^{22}\), the scan generated 2,671,046 possible summands and 3,191,130 possible targets. Of 3,567,244,702,581 unordered pairs, 798,439,545,154 survived the residue tests. The checker found 318,157 exact 4-full sums among those survivors and zero primitive triples. Thus any solution has \(\max(a,b)>10^{22}\). A second generator based on the canonical factorization reproduced both input lists.

[Manji-Jung's search through \(c\le9\cdot10^{18}\)](https://github.com/joaquimmj/erdos-939) · [Earlier computational report](https://www.erdosproblemaday.com/report/939) · [Full search checker](verify_erdos_939_bounded.cpp)

## ABC triples

Put \(R=\operatorname{rad}(abc)\). For any solution,
\[
R^4\mid abc,\qquad abc\le \frac{c^3}{4},\qquad
q(a,b,c)=\frac{\log c}{\log R}>\frac43.
\]
Thus every solution occurs among the \(abc\)-triples with \(R<c\). ABC@home describes its first archive as complete for \(c<10^{18}\). Its second archive contains additional triples below \(2^{63}\). The exact scan reads every row and uses the identity \(b=c-a\). It constructs the 4-full membership set through the archive height using the canonical factorization below. [ABC@home archives](https://pub.math.leidenuniv.nl/~smitbde/abc/index.php?set=1)

| Archive | Rows | 4-full targets | 4-full target and at least one 4-full summand | Primitive 4-full triples |
| --- | ---: | ---: | ---: | ---: |
| \(c<10^{18}\) | 14,482,065 | 1,269,560 | 100,171 | 0 |
| Additional rows with \(10^{18}\le c<2^{63}\) | 9,345,651 | 630,346 | 32,620 | 0 |

High-quality \(abc\)-triples have been examined in connection with this problem since at least Browkin and Brzeziński's 1994 table. The necessary threshold \(4/3\) is below the usual \(1.4\) cutoff for curated high-quality tables. [Browkin–Brzeziński, Remark 6](https://artefacts-discovery.researcher.life/full_text/DA-2/b4/b454fe563255359f879217a088e015db/full_text/8415a403048c5756f392b86da2303c78.pdf) · [de Smit's good triples](https://pub.math.leidenuniv.nl/~smitbde/abc/?set=2)

The [curated-table checker](verify_erdos_939_curated.py) independently disqualifies all 241 entries in de Smit's good-quality table by finding an exact prime valuation below four. Across the good, merit and unbeaten tables, the checker reconstructs 539 distinct triples and certifies 529 exclusions. Of the 371 entries with \(c\ge2^{63}\), it certifies 361 exclusions. The other ten entries have a hidden factor in the source and reported quality at most \(1.0883\). Their exclusions are supported by the tabulated quality and lack an independent valuation or radical certificate here. [Good table mirror](https://numberdb.org/files/T62/abc-quality.html) · [Merit table archive](https://web.archive.org/web/20251211064253id_/https://pub.math.leidenuniv.nl/~smitbde/abc/?set=3) · [Unbeaten table archive](https://web.archive.org/web/20251211063930id_/https://pub.math.leidenuniv.nl/~smitbde/abc/?set=4)

## Cohn's cubic tangent construction

Cohn's tangent step on a diagonal cubic gives a map on primitive signed triples. Write \(u+v+w=0\), where \(u,v,w\) are nonzero and pairwise coprime. Put \((d_1,d_2,d_3)=(v-w,w-u,u-v)\) and define
\[
T(u,v,w)=\frac{(u d_1^3,v d_2^3,w d_3^3)}{G},
\qquad
G=\begin{cases}
27&u\equiv v\equiv w\pmod 3,\\
1&\text{otherwise}.
\end{cases}
\]
The raw sum is zero by substitution of \(w=-u-v\). The map commutes with permutations up to a common sign. The only primitive nonzero seeds with a zero difference are the permutations and sign changes of \((1,1,-2)\), whose image has a zero coordinate. [Cohn, p. 439](https://pub.math.leidenuniv.nl/~smitbde/edu/bachsem/3-powerful_abc-triple.pdf)

**Proposition (Odd valuations and the first 4-full iterate).**

For any seed with three nonzero differences, \(T\) is primitive and preserves each positive odd-prime valuation on its original coordinate. The even coordinate has valuation at least \(5\) in the image. The image is 4-full exactly when all odd primes dividing the seed have exponent at least \(4\), every \(p\ge5\) dividing a difference divides it at least twice, and, if \(3\nmid uvw\), the unique difference divisible by \(9\) is divisible by \(27\). If the seed is outside the 4-full locus, the first iterate is its only possible first visit to that locus.

To prove this, the six linear forms \(u,v,w,d_1,d_2,d_3\) in \((u,v)\) have pairwise determinants in \(\{\pm1,\pm2,\pm3\}\). At \(p\ge5\) at most one form vanishes modulo \(p\), so an old valuation survives and a new valuation is \(3v_p(d_i)\). Exactly one of \(u,v,w\) is even, say \(u\). Then \(d_1=u+2v\) is even and the other two differences are odd. If \(v_2(u)=1\), then \(v_2(d_1)\ge2\), giving image valuation at least \(7\). If \(v_2(u)\ge2\), then \(v_2(d_1)=1\), giving image valuation \(v_2(u)+3\). In particular an even valuation of \(4\) in a solution rules out a Cohn preimage.

If \(3\) divides one coordinate, the other two are opposite units modulo \(3\), so all differences are units and the old valuation survives. Otherwise the three coordinates are congruent to the same unit \(s\bmod3\). Write them as \(s+3u_i\). The identity gives \(\sum_i u_i=-s\), so the three \(u_i\) are neither all equal nor all distinct modulo \(3\). Exactly one difference is divisible by \(9\), and the other two have valuation \(1\). Dividing by \(27\) leaves a power of \(3\) in just the coordinate of the exceptional difference and the new exponent is \(3(v_3(d_i)-1)\). These calculations give the content \(G\) and pairwise coprimality after division. They also show that after the first step all defects are odd-prime valuations below \(4\), which persist at every subsequent step in the domain.

**Corollary (Height bound for 4-full tangent images).**

If a primitive positive 4-full triple comes from \(T\), its largest term exceeds \(10^{72}/27\). The complete ABC@home archive below \(10^{18}\) contains no eligible primitive preimage. The additional archive below \(2^{63}\) contains none among its listed rows.

In fact, normalize a nondegenerate preimage to \((a,b,-c)\) with \(a<b\), \(a+b=c\). Then \(c>2\) and its odd prime valuations must all be at least \(4\). For \(R=\operatorname{rad}(abc)\), this gives \(R^4\le8abc\le2c^3<c^4\), since the exponent at \(2\) is at least \(1\). Such a seed must appear in the complete \(R<c\) archive if \(c<10^{18}\). The [seed checker](verify_erdos_939_cohn.py) tested all 14,482,065 rows in that archive and all 9,345,651 rows in the additional archive. It found no row whose three odd parts are 4-full. The height in the image is \(b(a+c)^3/G\), since the second signed coordinate is negative and the other two are positive. For \(t=b/c\in[1/2,1)\), the function \(t(2-t)^3\) is decreasing and exceeds \(1\). Thus the height in the image exceeds \(c^4/27\). [ABC@home archives](https://pub.math.leidenuniv.nl/~smitbde/abc/index.php?set=1)

**Corollary (Tangent seeds with bounded prime support).**

If every prime dividing a primitive positive triple is at most \(53\), or if its radical is at most \(10^7\), the triple has an odd prime valuation below \(4\) unless it is \((1,1,2)\). In either domain no triple is 4-full and no seed yields a 4-full tangent iterate.

For these two domains, von Känel and Matschke's Theorem B supplies the completeness of the normalized \(S\)-unit tables. The checker verified each equation, gcd, radical and supported factorization by exact division. The tables contain 1,663,449 and 432,408 rows, with largest terms 7,058,844,079,092,540,383 and 8,251,953,408 respectively. In each, \((1,1,2)\) was the only row with all odd prime valuations at least four. Its tangent step has a zero coordinate. The completeness theorem is an external input to the replay. [Theorem B](https://arxiv.org/pdf/1605.06079#page=8) · [Authors' tables 1 and 2](https://bmatschke.github.io/solving-classical-diophantine-equations/#s-unit-equation)

## Fixed-coefficient quartics

The canonical factorization of a 4-full integer is \(n=x^4d_5^5d_6^6d_7^7\), with \(d_5,d_6,d_7\) squarefree and pairwise coprime. For a prime of exponent \(4k+r\) with \(1\le r\le3\), put it in \(d_{4+r}\) and put \(p^{k-1}\) in \(x\). At residue \(r=0\), put \(p^k\) in \(x\). This recovers the factors uniquely and permits \(x\) to share primes with its own coefficient \(D=d_5^5d_6^6d_7^7\).

For fixed canonical coefficients \(A,B,C\), search the primitive equation \(Ax^4+By^4=Cz^4\). The coefficients have disjoint prime supports across terms. A local sieve checks the primes \(2,3,5,7,11,13,17,29,41\) and the coefficient primes. At a prime \(p\), at most one term can vanish modulo \(p\), and a nonunit term vanishes modulo \(p^4\). For each choice of vanishing term, the unit ratio of the other two has to be a fourth power. For three unit terms, normalize the first fourth power to \(1\) and check the remaining unit fourth powers. The checker uses modulus \(16\) at \(2\) and \(p\) at the odd primes.

For a fixed \(z\) put \(M=Cz^4\), \(X=\lfloor((M-1)/A)^{1/4}\rfloor\) and \(Y=\lfloor((M-1)/B)^{1/4}\rfloor\), with the bounds evaluated by integer division and integer roots. Since \(\gcd(Axy,M)=1\) in a primitive solution,
\[
 (x y^{-1})^4\equiv -BA^{-1}\pmod M.
\]
Factor \(M\) from \(C\) and \(z\). Find the unit fourth roots by successive modular square roots at each prime power and combine them by the Chinese remainder theorem. For a root \(r\) and a candidate pair, \(ry-x=kM\) for an integer \(k\). Since \(M\ge16\) for a nontrivial 4-full target and \(2XY<M\), one has
\[
 \left|\frac{r}{M}-\frac{k}{y}\right|=\frac{x}{My}<\frac{1}{2y^2}.
\]
Primitivity gives \(\gcd(x,y)=1\), hence \(\gcd(k,y)=1\). Legendre's continued-fraction criterion makes \(k/y\) a convergent. There is at most one reduced fraction \(x/y\) per root in the indicated box, since \(M\) divides any difference determinant and the determinant has absolute value less than \(M\). The Euclidean algorithm recovers the numerator and denominator, and the exact equation and gcd verify each candidate.

With \(A,B,C\le10^6\) and \(z\le150{,}000\), there are 31 possible coefficients and 2,431 primitive coefficient triples with the summands unordered. The sieve leaves 49 triples. The search tests 3,335,151 coefficient–target choices with nonempty coordinate bounds, of which 418,514 admit modular fourth roots. No primitive solution occurs in this domain. Its height ceiling is \(10^6(150{,}000)^4\). [Catalog and quartic checker](verify_erdos_939.py)

## Squares among the three terms

A 4-full square is the square of a square-full integer. The [square-case checker](verify_erdos_939_squares.cpp) uses this observation in three finite domains. It found zero primitive solutions in each.

| Squares required | Height \(c\) through | Search |
| --- | ---: | --- |
| All three terms | \(10^{40}\) | Euclid parameters for a primitive Pythagorean triple |
| The target and one summand | \(10^{30}\) | Coprime factors of a difference of squares |
| Both summands | \(10^{34}\) | Gaussian prime factors of the target |

For the first row, the square roots are \(u^2-v^2,2uv,u^2+v^2\), with \(\gcd(u,v)=1\) and opposite parity. The odd parts of \(u,v\) and both \(u+v,u-v\) must be square-full. The odd linear factors are coprime, so checking them separately is necessary. The range \(c\le10^{40}\) gives \(u^2+v^2\le10^{20}\). The search checked 4,167,496,440 pairs after the local residue and height tests. None had coprime parameters with both required square-full linear factors.

For a square summand \(x^2\) and square target \(z^2\), if \(x,z\) have opposite parity, \(z+x,z-x\) are coprime odd 4-full factors of \(b\). If \(x,z\) are odd, \((z+x)/2,(z-x)/2\) are coprime, one is even and divisible by 4, and their odd parts are 4-full. In both branches the checker reconstructs \(x,z\) and tests their square-fullness.

For two square summands \(x^2,y^2\), opposite parity is necessary since \(v_2(x^2+y^2)=1\) when both are odd. Primitivity prevents an inert prime or both conjugate Gaussian primes above a split prime from dividing \(x+iy\). Thus \(c\) has only primes \(p\equiv1\bmod4\). The checker enumerates each such rational prime with exponent at least four and chooses one of its Gaussian factors. Units and complex conjugation allow one representative for the first choice. The run reached 132,709,025 primitive norm representations, with 11,417,025 passing the local tests for square-full coordinates.

## Replay

Run from the repository root. The Python replays used Python 3.12.3, with SymPy 1.14.0 for the fixed-coefficient run. The native checkers used GCC 13.3.0, with OpenMP for the unrestricted search.

~~~sh
curl -fLO https://pub.math.leidenuniv.nl/~smitbde/abc/abctriples_below_1018.gz
curl -fLO https://pub.math.leidenuniv.nl/~smitbde/abc/big_triples.gz
python3 research/arithmetic/verify_erdos_939.py catalog abctriples_below_1018.gz 999999999999999999 14482065
python3 research/arithmetic/verify_erdos_939.py catalog big_triples.gz 9223372036854775807 9345651
.venv/bin/python research/arithmetic/verify_erdos_939.py fixed 1000000 150000
g++ -std=c++20 -O3 -march=native -fopenmp research/arithmetic/verify_erdos_939_bounded.cpp -o /tmp/erdos939-bounded
OMP_NUM_THREADS=12 /tmp/erdos939-bounded 10000000000000000000000 /tmp/erdos939-height22
python3 research/arithmetic/verify_erdos_939.py enumeration 10000000000000000000000
g++ -std=c++20 -O3 research/arithmetic/verify_erdos_939_squares.cpp -o /tmp/erdos939-squares
/tmp/erdos939-squares all-square 10000000000
/tmp/erdos939-squares square-c 1000000000000000
/tmp/erdos939-squares square-summands 100000000000000000
~~~

The archive SHA-256 digests are `8d40b46d6962aa14b2f4723eafcbaa68644a786a5b35db0d706a22ce51fd912c` for `abctriples_below_1018.gz` and `a7e6b9c755a156f44407fb944d425faf688e9a1398ca66015f8ca6afd0979c2b` for `big_triples.gz`. The independently generated native and canonical-factorization lists use sorted unsigned 16-byte little-endian values. Their SHA-256 digests are `f3629b893469607a71274662ada4dae823e4cf55a08bece81da1ffe20b776df5` for the summands and `fc56c9269c915823da7395f12298f515a7086863f7df226253d313f9dce3293b` for the targets.

For the curated replay, save the three linked HTML snapshots as `good.html`, `merit.html` and `unbeaten.html` in one directory. The [curated-table checker](verify_erdos_939_curated.py) pins their SHA-256 digests. Run `python3 research/arithmetic/verify_erdos_939_curated.py DIRECTORY`.

For the tangent seed checks, reuse the two ABC@home archives above and download the authors' two \(S\)-unit tables. The [seed checker](verify_erdos_939_cohn.py) verifies the SHA-256 digests of all four inputs.

~~~sh
curl -fL https://bmatschke.github.io/solving-classical-diophantine-equations/s-unit-equation/solutions_2_3_5_7_11_13_17_19_23_29_31_37_41_43_47_53.txt -o sunits53.txt
curl -fL https://bmatschke.github.io/solving-classical-diophantine-equations/s-unit-equation/solutions_uptoRadical_10000000.txt -o sunits-radical.txt
python3 research/arithmetic/verify_erdos_939_cohn.py abc-small abctriples_below_1018.gz
python3 research/arithmetic/verify_erdos_939_cohn.py abc-additional big_triples.gz
python3 research/arithmetic/verify_erdos_939_cohn.py primes-53 sunits53.txt
python3 research/arithmetic/verify_erdos_939_cohn.py radical-1e7 sunits-radical.txt
~~~
