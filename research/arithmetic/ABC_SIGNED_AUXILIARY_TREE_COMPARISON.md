# Signed comparison across auxiliary prime tests

A finite-field binary-tree estimate gives signed cancellation for paired
histories whose shared slots have a sufficiently connected overlap graph.
Its local mixed-character hypothesis holds with a power saving for the
centered transform of a suitably sparse residue image. Actual occurrence
dictionaries supply such images at auxiliary prime scales above their
cardinality exponent. A transfer still needs the arithmetic realization
of the trees and a signed lower bound from the genuine triples.

[Signed Jacobi sum obstacle](ABC_MODULAR_JACOBI_INCIDENCE.md#scope-at-the-prime-sizes-needed-for-true-occurrence)
· [Induced matchings with all three caps (IM1–IM3)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md)

## The tree estimate and the law of shared slots

Use the Fourier convention
\[
 \widehat F(v)=\mathbf E_{x\in\mathbb F_q}F(x)e_q(-vx),
 \qquad e_q(u)=\exp(2\pi i u/q)
\]
at an odd prime \(q\). The local input in Lemma 6.1 of the inverse
Goldbach manuscript is \(g(0)=0\),
\(\mathbf E_{v\in\mathbb F_q}|g(v)|^2\le1\), and a bound on
\[
 \epsilon_q(g)=
 \max_{\chi,\ a\in\mathbb F_q}
       |\mathbf E_{v\in\mathbb F_q}g(v)\chi(v)e_q(av)|,
\]
where \(\chi\) runs over all multiplicative characters of
\(U=\mathbb F_q^\times\), extended by zero at zero.

The allowed diagrams are full ordered binary trees of fixed depth
\(\ell\), with \(r=2^\ell\) leaves having unit variables \(M_i\).
Write \(M_j\) for the product of the leaves below a node \(j\).
The fixed unit data are a global \(D'\), node frequencies \(s_j\),
an internal triple \(C_{j,+},C_{j,-},u_j\) at each internal node,
and current entries \(X_{j,+},X_{j,-}\) fixed at the root. The
recursion is
\[
 \begin{split}
 H_{j,\pm}&=X_{j,\pm}C_{j,\pm}M_{j,\pm},\\
 p_j&=(s_{j,+}H_{j,-}-s_{j,-}H_{j,+})/(s_ju_j),\\
 y_j&=s_j/(D'H_{j,+}H_{j,-}),\\
 y_{j,\pm}&=s_{j,\pm}/(D'p_ju_jH_{j,\pm}).
 \end{split}
\]
At an internal child \(j,\pm\) pass down the current entries
\((p_j,X_{j,\pm})\). The fixed constants must obey
\(C_{j,\pm,+}C_{j,\pm,-}=u_j C_{j,\pm}\) whenever the child is
internal, so the two definitions of its \(y\) agree. Any vanishing
\(p_j\) makes the diagram value zero. Otherwise its value \(W\)
is a product of \(g(y_i)\) and \(\overline{g(y_i)}\) at the leaves,
with opposite conjugations in each bottom sibling pair.

With independent uniform unit leaf variables, Lemma 6.1 gives
\(\mathbf E|W|^2\le3^r\). For two such diagrams, partition each leaf
set into \(b\) nonempty blocks, pair those blocks, and take uniform
measure on the subgroup where the products of the leaf variables
of paired blocks agree. For \(\ell\ge2\) and \(b\le3r/4\), the
lemma gives
\[
 |\mathbf E W_1\overline{W_2}|
            \ll_\ell \epsilon_q(g)+q^{-1/4}.
 \tag{AT1}
\]
The same \(g\) is used in both diagrams at the given prime. The fixed
unit parameters and the sibling conjugation choices can differ.
[Inverse Goldbach manuscript, Lemma 6.1 and its hypotheses](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/paper.pdf#page=46)

The subgroup law follows for two arrangements of independent uniform
unit slots into leaves. Make each slot an edge of the bipartite graph
of first and second leaves. Leaf variables are incident-edge products.
In each connected component, the products on the two sides agree.
Conversely, prescribe any leaf values with these equalities, set
edges outside a spanning tree of each component to one, and solve
the remaining edge values in order from its terminal vertices. The
last vertex condition holds by the component product equality.
The resulting homomorphism is onto the stated subgroup, so the
pushforward is uniform. This verifies \(b\) as the number of connected
components, as in Section 8 of the manuscript.
[Arithmetic realization of the overlap law, Section 8](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/paper.pdf#page=65)

There is a direct stability version of (AT1). The two marginals of
this subgroup law are uniform on their respective full unit leaf
spaces, since each projection is a surjective homomorphism.
If the actual joint law has density \(w\) with
\(\|w-1\|_\infty\le \rho\) relative to the subgroup law, its signed
correlation differs by at most
\[
 \rho\mathbf E|W_1\overline{W_2}|
   \le \rho(\mathbf E|W_1|^2\,\mathbf E|W_2|^2)^{1/2}
   \le 3^r\rho .
 \tag{AT2}
\]
Merely bounding the density above by a constant supplies the
unsigned majorization and leaves signed cancellation to be proved.

### Different masks at different leaves

The argument of Lemma 6.1 also gives (AT1) when every leaf is assigned
its own mask, all satisfying the same local bounds with
\(\epsilon=\max_i\epsilon_q(g_i)\). The constant can depend on depth.
To see this, first allow a separate mask for each bottom sibling
pair, with opposite conjugations of that pair's mask. The one-tree
\(L^2\) proof takes a product of the separate squared norms, all at
most one. In the quartet step the source already keeps separate
bottom-pair coefficients \(P_\chi\), \(P_\psi\) and second moments
\(R_L,R_R\). Each pair's autocorrelation bounds and the uniform
Mellin bound in equation (6.14) depend only on its own \(L^2\) and
mixed-correlation bounds. In the two-nonprincipal case the
Cauchy–Schwarz bound uses the two separate fourth moments. In the
other cases it uses the corresponding individual Parseval and
maximum bounds. Consequently each quartet has its own coefficient
majorants with the same uniform small supremum and bounded sum.
The cycle projection only multiplies these majorants, so it
allows different ones at different quartets. The second diagram
needs just its one-tree \(L^2\) bound in that projection step.

For two different leaf masks \(f,g\) of a bottom pair, put
\(h_\omega=(f+\omega g)/2\) for \(\omega\in\{1,i,-1,-i\}\).
Then
\[
 f(x)\overline{g(y)}
    =\sum_{\omega\in\{1,i,-1,-i\}}
                        \omega h_\omega(x)\overline{h_\omega(y)}.
\]
Each \(h_\omega\) vanishes at zero, has \(L^2\) norm at most one,
and has mixed-correlation bound at most \(\epsilon\).
Expand the pairs this way and apply the preceding version to
each term. The number of terms depends only on depth. A sign
at any particular leaf can be absorbed into its mask by
conjugation, which preserves the local hypotheses. Thus the
local estimate can keep the three role masks separate.
[The quartet and cycle estimates in the proof of Lemma 6.1](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/paper.pdf#page=47)

## Sparse residue masks give a polynomial local input

**AT3.** Let \(S\subset\mathbb F_q\) have \(1\le s=|S|\le q/2\),
put \(\sigma=s/q\), and use the centered unweighted mask
\[
 F(x)=\frac{\mathbf 1_S(x)-\sigma}{\sqrt{\sigma(1-\sigma)}},
 \qquad g(v)=\sqrt q\,\widehat F(v).
\]
Then \(g(0)=0\), \(\mathbf E|g|^2=1\), and
\[
 \epsilon_q(g)\le
     \max\left\{\sqrt{\frac{q-s}{qs}},
                 \sqrt{\frac{s}{q-s}}\right\}.
 \tag{AT3}
\]

For the principal character, Fourier inversion gives
\(q^{-1/2}F(a)\) as the correlation before taking its absolute
value, since \(g(0)=0\). As \(\sigma\le1/2\), its largest possible
absolute value is \(q^{-1/2}\sqrt{(1-\sigma)/\sigma}\).
For a nonprincipal \(\chi\) the finite-field Gauss sum gives
\[
 \mathbf E_v g(v)\chi(v)e_q(av)
  =\frac{\tau(\chi)}{\sqrt q}\,
                  \mathbf E_x F(x)\overline{\chi(a-x)},\qquad
 \tau(\chi)=\sum_v\chi(v)e_q(v),\quad|\tau(\chi)|=\sqrt q.
\]
The formula includes \(x=a\), where the character sum vanishes.
The constant term in \(F\) contributes zero. Bounding the
remaining sum on \(S\) by \(s\) gives the second term of (AT3).

Take \(n=X^L\) at fixed \(L>0\). Suppose three complete dictionaries
of distinct positive integers at most \(X\) have sizes
\(n^{t+o_X(1)}\), with \(t>0\). For any fixed \(\beta>t\), choose
\(q\) uniformly among primes in \([n^\beta,2n^\beta]\).
Each distinct difference of two labels has at most
\(O_{L,\beta}(1)\) prime divisors in the interval. The prime number
theorem gives expected cardinality loss on reduction at most
\(n^{2t-\beta+o_X(1)}\) for each dictionary, by counting its
unordered congruent pairs. Markov's inequality and a union bound
give \(s=n^{t+o_X(1)}\) simultaneously for the three residue
images outside a proportion \(n^{-(\beta-t)/2+o_X(1)}\) of primes.
For each image as an unweighted mask, (AT3) now yields
\[
 \epsilon_q(g)\le
 n^{-\frac12\min(t,\beta-t)+o_X(1)}.
 \tag{AT4}
\]
This argument applies to the separate residue dictionaries of an
actual induced matching. Its conclusion uses their unweighted
masks and needs no prime-coverage hypothesis. The weighted version
below applies when a transfer produces weights with enough
effective support.

The local estimate also permits complex weights with suitable
effective support. Let \(a\) be a nonzero function supported on
\(S\), let \(\mu=\mathbf E_x a(x)\), and use
\(F=(a-\mu)/(\mathbf E|a|^2-|\mu|^2)^{1/2}\) in place of the
unweighted \(F\). With counting norms for \(a\), the nonprincipal
mixed correlations remain at most \(\sqrt{s/(q-s)}\), while the
principal correlations are at most
\[
 \frac{\|a\|_\infty/\|a\|_{\ell^2}+\sqrt{s}/q}{\sqrt{1-s/q}}.
 \tag{AT3w}
\]
Indeed \(|\mu|\le\sqrt{s}\|a\|_{\ell^2}/q\) by Cauchy–Schwarz.
Hence the normalizing denominator is at least
\(q^{-1/2}\sqrt{1-s/q}\|a\|_{\ell^2}\). The Gauss-sum
formula and the bound \(\|a\|_{\ell^1}\le\sqrt{s}\|a\|_{\ell^2}\)
give the nonprincipal assertion. Fourier inversion and the bound
on \(|\mu|\) give the principal assertion. This supplies a power
saving for a sparse mask when
\(\|a\|_{\ell^2}^2/\|a\|_\infty^2\) grows by a fixed power,
including signed weights with this property.

There is a sharp distinction when the mask is close to half the
field. For the set of nonzero squares in an odd prime field, write
\(\chi_2\) for the quadratic character. Then
\(\mathbf 1_S=(\mathbf 1_{v\ne0}+\chi_2)/2\).
Its Gauss sum yields
\[
 |\mathbf E_v g(v)\chi_2(v)|
    =\frac{q-1}{\sqrt{q^2-1}}\longrightarrow1 .
\]
Thus balance, normalization and an \(L^1\) norm tending to one can
hold together with a mixed-character bias of order one.

## What the positive input would need to retain

For the mask in (AT3), its normalized additive transform has
\[
 \|g\|_4^4=\frac{E(S)-s^4/q}{s^2(1-s/q)^2},
 \qquad
 \|g\|_1\ge\frac{s(1-s/q)}{\sqrt{E(S)-s^4/q}},
 \tag{AT5}
\]
where \(E(S)\) is the number of ordered additive quadruples in
\(S\) in \(\mathbb F_q\), and all norms use probability measure.
For \(v\ne0\) one has \(g(v)=(s(1-s/q))^{-1/2}\sum_{x\in S}
e_q(-vx)\). Expanding its fourth moment and subtracting the
zero-frequency contribution proves the equality. The inequality
follows by interpolating between \(L^1\) and \(L^4\) at \(L^2\).
For a complete dictionary \(T\) of size \(n^{t+o_X(1)}\), its
average modular additive energy over \(q\asymp n^\beta\) is at most
\[
 E_{\mathbb Z}(T)+n^{4t-\beta+o_X(1)}
\]
by the same divisor argument, applied to the nonzero four-term
errors. Choosing one representative per residue injects the
quadruples of \(S=T\bmod q\) into the corresponding modular
quadruples of \(T\). Thus any available upper bound on the genuine
additive energy can be used in (AT5) on most primes, with a small
normalized Markov allowance. The false-quadruple term is below
the diagonal exponent when \(\beta>2t\).

Sparse support alone gives weak lower bounds for \(\|g\|_1\).
For an interval of \(s\) residues with \(1\ll s\le q/2\), the
geometric-series bound yields
\(\|g\|_1\ll\log(2s)/\sqrt{s}\). Indeed, for the nonzero frequencies
the interval sum is bounded by a constant times
\(\min(s,\|v/q\|^{-1})\), whose probability \(L^1\) norm is
\(O(\log(2s))\). Small-radical one-letter dictionaries can have
these images at coarse scales. Fix \(e=3^k\), \(n=X^{3/e}\), and
take \(q\asymp n^\beta\) prime with \(q\equiv2\bmod3\) and
\(\beta<1/3\). Every nonzero residue is the image of a unique
unit root class under \(z\mapsto z^e\), and each root class meets
\([(X/2)^{1/e},X^{1/e}]\) for large \(X\). Choose one root
for each residue in an interval of length \(s=n^{t+o_X(1)}\)
with \(0<t<\beta\). Its power is tall, has radical at most
\(n^{1/3}\), and has least \(M\)-repair at most \(X^{M/e}\)
for \(M\ge e\). This is an obstruction for one-letter hypotheses
in that range. Its exponent \(t\) lies below the balanced complete
matching exponent.

The inverse Goldbach manuscript obtains its input masks from
disjoint reductions of Cartesian summand tails at each prime.
Its Corollary 4.2 gives qualitative \(\epsilon_q=o(1)\) outside a
set of small harmonic mass in log-log prime bands. Its Proposition
5.1 separately obtains a uniform \(L^1\) supply on a set of
positive harmonic mass through coverage of all sufficiently large
primes. In Sections 7–9 the positive statistic comes from that
coverage and the two-summand structure. The arithmetic comparison
introduces independent unit residues at spectator primes only
after treating the other arithmetic restrictions. An induced abc
matching provides genuine diagonal pairs and forbids cross
equalities. It gives neither Cartesian residue disjointness nor
coverage of an ambient prime population.
[Corollary 4.2, Proposition 5.1, and the paired-history construction](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/paper.pdf)

For occurrence, a transfer could use (AT4) without those two
Goldbach hypotheses if the masks it actually creates are the
sparse unweighted images. It would still need an amplitude with
a quantitative lower bound forced by the genuine rows, retaining
the third letter's radical and exact repair restrictions. It
would need integral histories realizing the signed recursion,
with those same restrictions retained in the estimates. The
paired upper bound must justify a subgroup law or its signed
approximation after conditioning the arithmetic and real data,
and control histories outside the connected-overlap condition
and the discarded zero terms. The caps and a prime test alone
provide no construction of the required positive statistic.
[False three-cap matchings at coarse mixed precision (MB1)](ABC_OCCURRENCE_MIXED_BIN_FALSE_MATCHINGS.md)
· [False representatives for support-prime and auxiliary tests (AU1)](ABC_OCCURRENCE_AUXILIARY_UNIT_TESTS.md)

Qualitative decay at each of finitely many fields does not entail
a fixed power of \(n\). Even a fixed contraction factor less than
one at every auxiliary prime would require a number of independent
factors of order \(\log n\) to yield a fixed power by multiplication
alone. A product of distinct primes bounded by \(n^{O(1)}\) has
only \(O(\log n/\log\log n)\) factors. In contrast, (AT4) does
give a power at a single polynomial-size field, once the correct
history measure is available. These statements concern the local
comparison and its hypotheses, with the height limit kept at
fixed budget and fixed depth.
