# Additive energy at the minimum nonatomic \(abc\) endpoint

Consider a limiting sequence of \(abc\) boxes with normalized
count exponent \(d=1/3\) and total nonatomic mass \(\Gamma=1/2\).
Write \(N=X^L\) for the root scale and \(c=1/6\) for the
single-letter capacity exponent. A subset of distinct complete
values is at capacity in the limit when its cardinality is
\(N^{c+o(1)}\), in the same order of limits as the complete-label
box estimates.

## Theorem. Equivalent doubling and energy exclusions

For every such subset \(S\) of a single complete letter,
\[
 E(S)=\#\{(a,b,a',b')\in S^4:a+b=a'+b'\}=o(|S|^3).
\]
The universal assertion over these subsets is equivalent to
the assertion that \(|S+S|/|S|\to\infty\) for every such subset.

For the forward transfer from doubling to energy, if
\(E(S)\ge\eta |S|^3\) along a subsequence for a fixed
\(\eta>0\), the Balog–Szemerédi–Gowers theorem gives
\(U\subseteq S\) of size at least \(c_\eta |S|\) with
\(|U-U|\le K_\eta |U|\). Apply Plünnecke's inequality to
\(U+(-U)\) to obtain a nonempty \(V\subseteq U\) with
\(|V-2U|\le K_\eta^2|V|\). Since \(V-2U\) contains a
translate of \(-2U\), it follows that
\(|U+U|\le K_\eta^2|U|\). This is bounded doubling on
a subset that still has cardinality \(N^{c+o(1)}\), contrary
to the complete-letter obstruction.
[Gowers, Proposition 7.3](https://www.cs.umd.edu/~gasarch/TOPICS/vdw/sz-thm-gowers-proof.pdf) · [Petridis, Theorem 3.1](https://arxiv.org/pdf/1101.3507v3) · [The complete-letter obstruction](ABC_DIFFUSE_FACTOR_GRAPH.md#a-uniform-obstruction-to-bounded-rank-additive-progressions)

For the reverse transfer, Cauchy–Schwarz gives
\(|S|^4\le |S+S|E(S)\). Thus the energy estimate forces
growing doubling for every such \(S\).

## Corollary. Higher moments and additive edge density

For each fixed integer \(m\ge2\), the number \(E_m(S)\) of
ordered equal sums of two \(m\)-tuples in \(S\) is
\(o(|S|^{2m-1})\). Indeed, Young's convolution inequality
gives \(\|1_S^{*m}\|_2\le |S|^{m-2}\|1_S*1_S\|_2\).

Let \(A,B,C\) be subsets at capacity in the limit from the
three complete letters. If \(M\) counts selected triples with
\(a+b=c\) in \(A\times B\times C\), then
\[
 M=o\!\left((|A||B||C|)^{2/3}\right).
\]
By Cauchy–Schwarz,
\(M^2\le |C| E(A,B)\), where
\(E(A,B)=\sum_h r_{A-A}(h)r_{B-B}(h)\le\sqrt{E(A)E(B)}\).
Use the same argument for the other two pairs, reversing signs
as necessary. Multiplying the three bounds gives
\[
 M^3\le (|A||B||C|)^{1/2}(E(A)E(B)E(C))^{1/2}.
\]
The energy estimate on each letter proves the claimed cut.
For comparable part sizes the proportion of pairs that give a
selected triple tends to zero. The resulting bound is
qualitative and permits the limiting normalized count
exponent \(1/3\).
