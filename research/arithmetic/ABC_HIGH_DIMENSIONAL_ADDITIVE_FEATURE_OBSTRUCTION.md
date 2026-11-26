# High-dimensional complete letters with no large small-feature subset

Exact addition, primitivity, a common height shell, and the near-quadratic complete-label edge exponent do not by themselves force a large subset in the fixed-rank affine models used by (CF1). A carry-free digit construction matches the endpoint cardinality scales and has subpower doubling on each actual letter. Long sumsets of every large subset give an obstruction to small maximal feature minors, uniformly in the real coefficients of the putative affine model.

## Maximal feature minors control long sumsets

Let \(U\subset\mathbb Z\) be finite and nonempty. Suppose each \(T\in U\) has an integer feature vector \(z(T)\in\mathbb Z^R\) and \(T=\ell(z(T))\) for one real affine function \(\ell\). Let \(r\) be the affine dimension of these vectors. Define \(\Delta\) to be the maximum absolute value of every \(r\)-by-\(r\) coordinate minor of \(r\) feature differences from one common feature anchor, where the anchor, the other \(r\) vectors, and the coordinate columns range over their possibilities. Set \(\Delta=1\) when \(r=0\). Then for every positive integer \(m\),
\[
                         |mU|\le \Delta(2m+1)^r.
 \tag{HF1}
\]
Here \(mU\) is the set of sums of exactly \(m\) elements of \(U\), with repetitions allowed.

For \(r>0\), choose an anchor, \(r\) other features, and \(r\) coordinate columns attaining \(\Delta\). Form the square integer matrix \(B\) whose columns are their projected differences. The coordinate projection is injective on the affine span, because those differences project to independent vectors. In the basis given by \(B\), each projected feature difference has coordinates in \([-1,1]^r\). Cramer's rule gives this bound, since replacing one column by that projected difference yields another anchored minor of absolute value at most \(\Delta\). Each coordinate vector belongs to \(B^{-1}\mathbb Z^r\). This lattice contains \(\mathbb Z^r\) with quotient of order \(|\det B|=\Delta\), by the map \(x\mapsto Bx\) onto \(\mathbb Z^r/B\mathbb Z^r\). Its intersection with \([-m,m]^r\) contains at most \(\Delta(2m+1)^r\) points, since each translate of \(\mathbb Z^r\) has at most \((2m+1)^r\) points in that cube. Every sum of \(m\) feature vectors corresponds to one such point. Applying the affine formula, with \(m\) times its constant term, bounds the number of sums of \(m\) original values. A zero-dimensional feature space gives a singleton. This proves (HF1) with the same original maximal-minor condition as (CF1). If several features are available for the same value, apply the argument to the original feature set before taking its image, using the original rank and maximal minors. [Original maximal-minor condition (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)

## Exact primitive graphs at the endpoint cardinality scale

Fix \(0<L<1\) and put \(t=6/L-1\). Let \(d\) tend to infinity, choose positive integers \(H=H(d)\) with \(\log H\to\infty\) and \(\log d=o(\log H)\), and set
\[
     N=H^d,\qquad m=\lceil d H^t\rceil,\qquad
     Q=4mH,\qquad M=Q^d .
\]
For \(\boldsymbol x=(x_0,\ldots,x_{d-1})\), write
\(F(\boldsymbol x)=\sum_{j=0}^{d-1}x_j Q^j\). Let
\[
     S=F(\{0,\ldots,H-1\}^d).
\]
For each \(a_0,b_0\in[M,2M)\cap\mathbb Z\), take the actual sets
\[
 A=a_0+S,\qquad B=b_0+S,\qquad C=a_0+b_0+(S+S).
\]
Every pair in \(A\times B\) gives an edge in the addition graph. Digit expansion gives \(|A|=|B|=N\), \(|C|=(2H-1)^d=N^{1+o(1)}\), and \(\max S\le M/(4m)\). For \(X=a_0+b_0+2\max S\), all corresponding targets lie in \((X/2,X]\) and all three letters lie in \([X/5,X]\). The height \(X\asymp M\) satisfies \(\log X=(t+1+o(1))\log N\). Thus \(N=X^{L/6+o(1)}\) at each fixed \(L\).

Choose \(a_0,b_0\) independently uniformly among the \(M\) integers in \([M,2M)\). For a fixed \(s_1,s_2\in S\) and any prime \(p\le3M\), the probability that \(p\) divides both \(a_0+s_1\) and \(b_0+s_2\) is at most \((1/p+1/M)^2\). The probability of a nonprimitive pair is therefore at most
\[
 \sum_{p\le3M}(1/p+1/M)^2
     \le \sum_{j\ge2}j^{-2}+O(\log M/M)<1
\]
for sufficiently large \(M\), with an absolute fixed distance below one. No larger prime can divide either positive summand. Averaging gives shifts for which the graph has \(\gg N^2=X^{L/3+o(1)}\) primitive edges. Passing to the part labels that occur on these edges leaves \(\gg N\) labels in each part, since a fixed letter has at most \(N\) partners from one other part. Every fixed \(k\)-fold sumset on a part has cardinality at most \((2kH)^d=N^{1+o(1)}\), before or after this deletion. In particular each active part has subpower doubling and additive energies of all fixed orders \(2k\) at least \(N^{2k-1-o(1)}\), by Cauchy–Schwarz on the \(k\)-fold representation function.

## Each part excludes large low-rank, small-minor models

Let \(U\) be any nonempty subset of one of these three part sets. Removing the fixed translation places it in \(F(\{0,\ldots,2H-2\}^d)\). Its digit vectors have some affine dimension \(d'\). Projection onto a suitable \(d'\) digit coordinates is injective on their affine span, so
\[
                   |U|\le (2H)^{d'}.
\]
The \(m\)-fold digit sums have individual digits at most \(m(2H-2)<Q\), hence no carries. Select \(d'+1\) affinely independent digit vectors of \(U\). The sums from these vectors with nonnegative multiplicities of total \(m\) are distinct, by affine independence, and their images under \(F\) are distinct by the no-carry condition. Consequently
\[
 \begin{split}
            |mU|&\ge {m+d'\choose d'},\\
            \log |mU|&\ge
             \frac{t\log H}{\log(2H)}\log |U| .
 \end{split}
 \tag{HF2}
\]
For \(d'>0\) the second inequality uses \({m+d'\choose d'}\ge(m/d')^{d'}\) and \(m/d'\ge H^t\). For \(d'=0\) both statements hold with \(|U|=1\).

Combining (HF1) and (HF2) gives a uniform constraint on every affine feature presentation of \(U\) of affine rank \(r\le d\) and maximal anchored minor at most \(X^\theta\), for fixed \(\theta\ge0\):
\[
             \frac{\log |U|}{\log N}
                \le \frac{t+1}{t}\theta+\frac rd+o(1).
 \tag{HF3}
\]
The error is independent of the presentation and of which subset is chosen. In particular if \(r=o(d)\), no such model contains \(N^\sigma\) labels for any fixed \(\sigma>(t+1)\theta/t\). For fixed rank this applies to the small-minor class of (CF1) at any threshold \(\theta<\sigma t/(t+1)\). Since each complete label supports at most \(N\) primitive edges, the same conclusion prevents an \(N^{1+\sigma}\)-edge subfamily with such a model on any one letter. A cover by \(N^{o(1)}\) such models still captures at most \(N^{(t+1)\theta/t+o(1)}\) labels in that part.

A version for ordinary progressions follows without the determinant comparison. For a rank-\(r\) progression of coefficient-box volume \(V\), the \(m\)-fold sumset has at most \(m^r V\) values, by summing the coefficient intervals. In particular any progression with \(r=o(d)\) and \(V\le N^K\) meets a part in at most \(N^{K/t+o(1)}\) values. The maximal-minor proof handles arbitrary real output coefficients and feature sets of arbitrary shape, rather than only a given progression presentation.

The earlier primitive graph in the modular progression analysis uses random generators and tests fixed-rank progressions after reduction modulo a prescribed modulus. The construction here is explicit in its digits. It tests arbitrary integer features through their original maximal minors and it uses long integer sums for the obstruction. Neither construction imposes the small radical caps on the three parts. Thus the example establishes a limitation on structural extraction using the integer additive graph, primitivity, numerical heights, and even all fixed-order additive energies. It does not contradict the arithmetic count (CF1). [Earlier modular progression obstruction](ABC_MODULAR_PROGRESSION_SLICING.md#a-primitive-additive-graph-with-uniform-progression-avoidance)

## Two inverse transfers that use the radical caps

Write \(n=X^L\) and choose one prescribed complete letter. Take \(S\) to be the set of its actual values in a family of positive primitive triples with \(X/2<c\le X\) and \(\operatorname{rad}(abc)<c^L\). Each of the following two conditions gives a strict uniform improvement below the one-third exponent, for all sufficiently small fixed \(L>0\). They impose no condition on the limiting valuation profiles. The threshold in \(L\) and the size of the saving may depend on the fixed parameters of the condition.

To reduce to the exact hypotheses of (CF1), once its fixed feature rank \(R\) is chosen, discard the triples on which one individual radical exceeds \(n^{1/3+\kappa_R}\). This sector has a uniform exponent gap by (GD3): bin the individual normalized radical masses \(l_T\) within a sufficiently small fixed mesh, so that \(\sum l_T\le1+o(1)\) up to mesh error and \(\max l_T\ge1/3+\kappa_R\). The two-letter coefficient \((\sum l_T-\max l_T)/2\) is then below \(1/3\) by a fixed amount. Choosing exact valuation truncation losses and the (GD3) allowance smaller than the gap preserves it. Both conditions below are inherited by the remaining subset of \(S\). [Two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

**A many-fold sumset of bounded normalized size.** Fix \(b,u>0\), set \(k=\lfloor n^u\rfloor\), and suppose
\[
                              |kS|\le n^b .
 \tag{HF4}
\]
For nonempty \(S\), translate by an element of \(S\) so that the translated set contains zero. Its \(k\)-fold sumset has the same size, and this is at most \(k^\gamma |S|\) for all large \(X\) with fixed \(\gamma=2b/u\), since \(k\ge n^{u/2}\). Nguyen–Vu's long-range inverse theorem gives a proper symmetric integer progression of rank \(r=O(\gamma)\) and size \(O_\gamma(k^{-r}|kS|)\) containing the translated set. Thus the rank is bounded in terms of \(b,u\), and the progression size is \(O_{b,u}(n^b)\). Its integer coefficient features have absolute values \(O_{b,u}(n^b)\), and every required maximal minor of rank at most \(R\) is \(O_{b,u}(n^{bR})\). For sufficiently small fixed \(L\) this is within (CF1). This proves the strict saving on (HF4). In particular, in a balanced saturating sequence every subfamily that still saturates has \(|\lfloor n^u\rfloor S|> n^b\) for every fixed \(u,b>0\), in the height-first convention. Bounds on a fixed number of additive moments do not give (HF4). [Nguyen–Vu, Theorem 3.2](https://arxiv.org/pdf/1004.3967v2#page=10)

**Short fixed-rank progressions at many smaller primes.** Fix a positive integer \(r\) and \(0<\vartheta<1/(r+2)\). There are constants \(C_{r,\vartheta},\omega_{r,\vartheta}>0\) such that it suffices to have a scale
\[
                   (\log X)^{C_{r,\vartheta}}\le y\le X^{\omega_{r,\vartheta}}
\]
and at least \(y/(4\log y)\) primes \(p\in[y,2y]\) for which \(S\bmod p\) lies in a GAP of rank at most \(r\) and cardinality at most \(p^\vartheta\). Indeed take \(0<\eta<1/2\) and \(\varepsilon_0>0\) sufficiently small in terms of \(r,\vartheta\) that
\[
              \vartheta < \frac{(1-\vartheta)(1-\eta)}{r+1-\eta}
                                      -\varepsilon_0,
   \qquad \varepsilon_0<(1-\vartheta)/r .
\]
The choice \(\delta=(2y)^{-\vartheta}\) is a uniform lower bound for the density of each nonempty reduction inside its prescribed progression, and is above the density threshold of Croot–Mao–Yip's Theorem 1.11 for sufficiently large \(y\). Its other height hypothesis holds by making \(C_{r,\vartheta}\) large enough. The theorem covers \(S\) by a bounded number, in terms of \(r,\vartheta\), of translates of a rank-at-most-\(r\) integer GAP of cardinality
\[
                     O_{r,\vartheta}\big(y^{\vartheta(1+r/(1-\eta))}|S|\big).
\]
After the balanced restriction above, \(|S|\le n^{1/3+\kappa_R+o_X(1)}\) by radical and fixed-support enumeration. Even if the resulting progression is improper, each side with a nonzero step has length at most the cardinality of the progression, as varying that side alone gives distinct integers. Discard zero steps. The anchored feature minors are bounded by a constant times the \(r\)-th power of the resulting cardinality bound. First take \(\omega_{r,\vartheta}\) sufficiently small and then \(L\) sufficiently small to put all those minors within (CF1). The finite cover preserves the exponent loss from (CF1). This is a gluing criterion across many small moduli. The earlier modular progression slicing argument instead begins from one prescribed modulus and counts the values that land in its progression. No deduction of the many-prime hypothesis from balanced saturation is presently available. [Croot–Mao–Yip, v2, Theorem 1.11](https://arxiv.org/html/2510.08862v2#S1) · [Single-modulus slicing](ABC_MODULAR_PROGRESSION_SLICING.md#fixed-rank-modular-progressions-on-regularized-complete-letters)

## What the quantitative inverse results provide

Gowers–Green–Manners–Tao, Theorem 1.3, gives for a set of doubling at most \(K\) in \(\mathbb Z^D\) a subset with polynomial loss in \(K\) and affine dimension \(O(\log K)\). Directly in the integers the affine dimension is already at most one and the statement has no feature-size conclusion. For a high-dimensional integral lift, the theorem can return an unbounded dimension when \(K=N^{o(1)}\), and does not bound the anchored minors of integer coordinates. Raghavan, v2, Theorem 1.4 gives a convex coset progression cover in an abelian group, with dimension and the logarithmic size and cover overhead \(O(\log(3K)\log\log(3K))\). The extra \(\log\log\log\) factor occurs in v1. The improved bound still allows unbounded dimension and, from the hypothesis \(\log K=o(\log N)\) alone, even permits the overhead estimate to exceed \(\log N\). [Gowers–Green–Manners–Tao, Theorem 1.3](https://arxiv.org/html/2311.05762v2#S1) · [Raghavan, v2, Theorem 1.4 and Question 1.2](https://arxiv.org/html/2512.11217v2#S1)

The recent energy form of Freiman's lemma gives a set of size comparable to the ambient affine dimension whose mixed energy with the original set is \(O(|A||S|)\). It concerns a mixed low-energy test, which can coexist with the near-maximal fixed-order self energies above. The many-fold rank lemma of Harrison–Mudgal–Schmidt assumes a bound \(|mA|\le K|A|\) in a one-dimensional algebraic group and gives a subset of rank \(O(1+\log(4K)/\log m)\). In the carry-free range, (HF2) allows \(\log(|mA|/|A|)\) of order \(\log N\) and \(\log m\) of order \((\log N)/d\). Thus the rank bound obtained from that input scale need not be bounded. It also concerns subgroup rank and does not bound numerical feature minors. [Mudgal–Schildkraut–Zakharov, Theorem 1.1](https://arxiv.org/html/2611.04644#S1) · [Harrison–Mudgal–Schmidt, Lemma 6.1](https://arxiv.org/html/2603.06483#S6)

Shparlinski–Stewart, Theorem 1.2, counts sets of \(s\) primes up to \(H\) that actually occur in primitive solutions to a prescribed \(S\)-unit equation with fixed coefficients. Its unconditional bound for three varying letters uses a fixed short-letter condition and has main \(H\)-exponent \(\lfloor 2s/3\rfloor\) up to logarithmic factors depending on \(s\). Removing that shortness in their stated estimate is conditional on \(abc\). These bounds do not give a varying-product radical budget estimate on balanced complete letters, with unbounded support length and no short-letter gap. [Shparlinski–Stewart, Theorem 1.2](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/countingsunitsr.pdf#page=5)
