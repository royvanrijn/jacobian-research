# Simultaneous angular and ray class tests from private prime packets

At a balanced limiting abc profile, every letter supplies a growing set of pairwise coprime full divisors of positive numerical weight. They give projective residue tests for rational directions and scalar ray tests for norm ideals. A nonzero secant or conjugate difference bounds the sum of the weights of all matching tests. The argument allows the aggregate least common multiple to have superpolynomial height.

## Private packets of uniformly positive numerical weight

Use the height-first convention at small budget \(L>0\), with \(n=X^L\). The positive-scale lemma in (SH1) gives absolute \(0<r_0<R_0\) and \(\gamma>0\) such that every letter of every limiting profile at normalized exponent \(1/3\) has degree mass at least \(\gamma\) in \([r_0,R_0]\). Enlarge to \([r_0/2,2R_0]\) and work at a fixed truncation past the upper end of the corresponding valuation range. Weak convergence gives normalized radical mass at least \(\gamma/2\) there on each letter of sufficiently late boxes. Inner height allowances can be taken smaller than this margin. [Positive scales in every endpoint letter](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)

Fix \(C>0\). Choose the private packets of (GP2) with loss less than \(\gamma/16\) on each letter. Keep only the full powers in the enlarged band. At sufficiently late fixed budgets and then at large \(X\), this gives at least \(\tfrac12\lfloor(\log X)^C\rfloor\) actual rows, with mutually coprime packets on all three roles across the selected rows. Each packet has \(\log M\ge\nu\log X\) for an absolute \(\nu>0\). For any subsequently prescribed \(P\le X^B\) one can delete all primes of \(P\) from these packets at the cost of discarding \(O_{B,L}(1)\) rows. Indeed (GP3) bounds the number of rows for which this deletion costs more than \(\gamma/16\) of normalized radical mass. Each remaining packet still has radical mass at least \(\gamma/4\) in the band, so one can use \(\nu=r_0\gamma/8\) even after the deletion. This applies to a common projective modulus or a polynomial-height field discriminant chosen after seeing the sample. Each such deletion uses one fixed pool for the subsequent comparisons. The resulting packets always permit the full scalar modulus below. The omitted-arm version requires its additional splitting hypothesis. [Growing selections and uniform pool avoidance](ABC_GROWING_PRIVATE_PRIME_PACKETS.md#growing-selections-and-uniform-pool-avoidance)

## Projective tests in an angular interval

Consider distinct primitive directions \((u:v)\), with \(0<u\le U\), \(0<v\le V\), whose angles after scaling by \((U^{-1},V^{-1})\) lie in an interval of length \(A\). Suppose they are in one projective class modulo \(M_0\), allowing \(M_0=1\). For pairwise coprime \(M_i\), all coprime to \(M_0\), let a match of two directions at \(i\) mean equality in \(\mathbb P^1(\mathbb Z/M_i\mathbb Z)\). The sum of \(\log M_i\) over their matches is at most
\[
       H_*=\log\max(1,2UV A/M_0),                              \tag{PSR1}
\]
since their nonzero determinant is divisible by \(M_0\prod_{\rm matches}M_i\) and is at most \(2UV A\) in magnitude. This is the pair bound behind the projective larger sieve with an angular gain. [Zywina, Theorem 3.1 and Lemma 3.3](https://pi.math.cornell.edu/~zywina/papers/Quantitative-HIT.pdf)

For each \(i\) allow a list of \(K_i\) projective classes, and suppose every direction in a query family \(\mathcal Q\) is in its allowed list at \(s\) or more indices. If every \(M_i\ge X^\nu\) and \(k\le s\) is a fixed integer with \(k\nu\log X>H_*\), then
\[
 |\mathcal Q|\binom{s}{k}\ \le\
          \sum_{\substack{I\subseteq\{1,\ldots,m\}\\|I|=k}}
                    \prod_{i\in I}K_i .                         \tag{PSR2}
\]
For each \(k\)-set \(I\), any specified class vector determines at most one direction by (PSR1). Counting matches over all such \(I\) proves (PSR2). In particular, with \(m\) polylogarithmic in \(X\), \(K_i\le X^{\kappa L}\) and \(s\ge k\), the count is at most \(X^{k\kappa L+o_X(1)}\). The lists can be chosen after the private packets are selected. The bound applies to cofactor secants in a fixed affordable record using their actual angular scale and any common projective modulus. [Angular secant estimate (PS1)](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#fixed-record-setup)

## A collision-weight bound for scalar ray tests

Let \(K\) be a number field of degree \(d\ge2\). Let \(M_1,\ldots,M_m\ge2\) be pairwise coprime rational integers, prime to \(2\operatorname{Disc}(K)\). For each \(i\) use one of two ideal moduli. The full modulus is \(\mathfrak m_i=M_i\mathcal O_K\) and has weight \(w_i=\log M_i\). In degree \(d\ge3\), the other option for each \(p^a\parallel M_i\) keeps all primes above \(p\) at exponent \(a\) except for one prescribed degree-one prime above \(p\). This modulus has weight \(w_i=(1-2/d)\log M_i\). Write \(S=\sum_i w_i\).

Use (PR3) with this \(\mathfrak m_i\), denoting its compact group by \(G_i\), the ideal map by \(g_i\), and the image of the centered logarithms by \(h_i(V_\infty)\). Explicitly the quotient defining \(G_i\) identifies \(((\gamma),-L(\gamma))\) with zero if \(\gamma\) reduces to a rational scalar unit at all retained primes, with the same rational unit across the modulus. Here
\[
  L(\gamma)_\sigma=\log|\sigma(\gamma)|-\tfrac1d\log|N\gamma|,
 \qquad C_W=\{t\in V_\infty:\|t\|_\infty\le W\}.
\]
No real sign condition is required. Fix \(W>0\), independent of \(X\). For each rational prime allowed in the ideal dictionaries, require it to be unramified in \(K\) and select at most one degree-one prime of \(K\) over it. The selection is the same for all the tests. [Scalar ray groups and the single-modulus comparison (PR3–PR5)](ABC_PROJECTIVE_RAY_CLASS_NORM_MOMENTS.md#the-compact-group-and-the-diagonal-moment)

For distinct integral ideals \(\mathfrak A,\mathfrak B\) supported on the selected primes with \(N\mathfrak A,N\mathfrak B\le X^H\), and any positive integer \(e\), we have
\[
 \sum_{\substack{i:\,(\mathfrak A\mathfrak B,\mathfrak m_i)=1\\
             e\,g_i(\mathfrak A/\mathfrak B)\in h_i(C_W)}}w_i
       \ \le\ C_{d,W}\bigl(eH\log X+1\bigr).                    \tag{PSR3}
\]
For each matching test there is a generator \(\gamma_i\) of \((\mathfrak A/\mathfrak B)^e\) with \(\|L(\gamma_i)\|_\infty\le W\) and the scalar residue condition. These generators differ by units. Relative to any one of them, the ratio unit has all conjugates bounded by \(\exp(2W)\). Bounded degree and the integer minimal-polynomial coefficient bounds give a uniformly finite list of such units. Thus only \(O_{d,W}(1)\) generators can occur.

Fix one used generator \(\gamma\). The integer \(D_2=(N\mathfrak B)^e\) clears its denominator, and every conjugate of \(D_2\gamma\) has magnitude at most \(\exp(W)X^{eH}\) because
\[
 |\sigma(D_2\gamma)|
 \le \exp(W)\,(N\mathfrak A)^{e/d}(N\mathfrak B)^{e(1-1/d)}.
\]
A rational \(\gamma\) would force \(\mathfrak A=\mathfrak B\). At every allowed rational prime a nonzero valuation of a rational generator would occur at all prime ideals above it, whereas the dictionary permits one degree-one prime and the unramified prime has at least one other prime upstairs. Choose unequal conjugates of \(D_2\gamma\), and let \(\Delta\ne0\) be their difference in a normal closure \(E\) of degree \(d_E\).

For a full modulus all primes of \(E\) over \(p^a\parallel M_i\) contribute at least \(a\) to their respective valuations of \(\Delta\), since both conjugates have the same rational residue. For an omitted-arm modulus the field \(E\) is unramified at \(p\). Indeed the completions of \(K\) and their local conjugates at \(p\) have unramified compositum. For either embedding of \(K\) in the difference, exactly a fraction \(1/d\) of the primes of \(E\) over \(p\) restrict to the omitted degree-one prime. The total count over \(p\) is \(d_E/f_E\), where \(f_E\) is the common residue degree, and the count over that one prime is \(d_E/(df_E)\). Thus at least a fraction \(1-2/d\) of the primes retain both congruences, each to depth \(a\). Over the probes using \(\gamma\) all rational prime supports are disjoint. Therefore
\[
    \exp\!\left(d_E\sum_{i\text{ using }\gamma}w_i\right)
         \le |N_{E/\mathbb Q}\Delta|
         \le \bigl(2\exp(W)X^{eH}\bigr)^{d_E}.
\]
Sum over the bounded number of used generators to obtain (PSR3).

## Averaged moments and a diagonal joint test

Let \(\mathcal D\) be a dictionary of \(B\) integral ideals of norm at most \(X^\alpha\), with the common prime selection above. Fix integers \(k,e\ge1\). For each \(i\) let \(\mathcal D_i\) consist of members prime to \(\mathfrak m_i\), and put \(F_i(\chi)=\sum_{\mathfrak a\in\mathcal D_i}\chi(e\,g_i(\mathfrak a))\). For \(G_i\) use the autocorrelation \(K_i\) of \(\mathbf 1_{h_i(C_{W/2})}\), divided by the Haar measure of \(h_i(C_{W/2})\) in \(G_i\). Then \(K_i(0)=1\), \(0\le K_i\le1\), its support lies in \(h_i(C_W)\), and \(\widehat K_i\) is a nonnegative probability weight by Parseval. The probability statement also holds for a zero-dimensional torus. In the height limit with \(d,W,\alpha,k,e\) fixed,
\[
 \frac1S\sum_i w_i\sum_\chi\widehat K_i(\chi)|F_i(\chi)|^{2k}
     \ \le\ B^k X^{o_X(1)}
       +\frac{C_{d,W}(ek\alpha\log X+1)}{S}B^{2k}.             \tag{PSR4}
\]
Indeed expand each moment over pairs of products of \(k\) ideals. The identical products cost at most \(B^k X^{o_X(1)}\) by the ordinary divisor bound for factorizations of ideals of bounded-degree polynomial norm. For a nonidentical pair use (PSR3) with \(H=k\alpha\). Allowing the full dictionary for this upper bound handles probes at which some products are ineligible.

For any subcollection \(I\) with \(\sum_{i\in I}w_i>C_{d,W}(ek\alpha\log X+1)\), restrict \(\mathcal D\) to ideals prime to all its moduli. Sum the products \(\prod_{i\in I}\chi_i(e\,g_i(\mathfrak a))\) over this joint dictionary. Its \(2k\)-th moment against \(\prod_{i\in I}\widehat K_i\) is exactly the number of pairs of equal ideal products. In the expansion a nonidentical product pair cannot pass all the tests by (PSR3), while every identical pair contributes one. Thus the joint moment is diagonal whenever the chosen probes jointly cross the stated threshold. With the positive-weight packets above, one may use a fixed number of probes for this purpose. No bound on the product of all available moduli is needed.

### Sharpness of the off-diagonal averaging term

Take \(K=\mathbb Q(i)\), full scalar moduli \(M_i\) of comparable logarithmic weight \(\nu\log X\), and \(M_1=3^{\lceil\nu\log X/\log 3\rceil}\). For example, raise distinct odd primes of size \(X^{o_X(1)}\) to their first powers above \(X^\nu\). For fixed \(\beta>\nu\) take the ideals \((a+iM_1)\) with \(a\le X^\beta/2\) positive even and prime to \(3\). The norm is odd and the two coordinates are coprime, so every good prime used is split and only one arm occurs in any one ideal. Independent uniform choices of the allowed arm at each split prime retain \(X^{-o_X(1)}\) of the dictionary, because each norm has \(O(\log X/\log\log X)\) distinct prime factors. The ideals are distinct and their norms are at most \(X^{2\beta}\) for large \(X\). Each generator reduces to the rational unit \(a\) at \(M_1\), hence its ray record is identical in \(G_1\), whose logarithmic torus is trivial. The first probe alone contributes \((w_1/S)B^{2k}\asymp B^{2k}/m\) to (PSR4). The associated triples \((a^2,M_1^2,a^2+M_1^2)\) are primitive. Their target radical is unrestricted in this construction.

The projective and ray tests therefore give a count when the candidate directions admit small residue lists or share suitable scalar ray records at several of the private packets. Deriving that simultaneous constraint from the three original radical caps remains the arithmetic step. In particular the averaged moment and the existence of polylogarithmically many packets alone do not eliminate a power-sized off-diagonal term.
