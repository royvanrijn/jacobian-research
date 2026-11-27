# The degree-fourteen alternating Hurwitz passport has one absolute orbit

There is one absolute full-alternating braid orbit for
\[
             \bigl((2,2,1^{10})^7,(13,1)\bigr).
\]
This is a finite classification, obtained from two complete degree-thirteen
polynomial censuses and the retained uniform folding and attachment
theorems. It does not prove connectivity in every degree. The next
unclassified even degree on this route is sixteen.

The largest census visits 446,265,625 absolute polynomial tuples. Its
independent total bound comes from exact finite-group character arithmetic,
rather than from a presumed enumeration formula or random exploration.

## The two polynomial domains

Put \(P=(0\,1\,\cdots\,12)\), write DT for a double transposition, and use
\[
               Q_i(a,b)=(b,b^{-1}ab).
\]
Consider ordered tuples whose product is \(P\):
\[
 \mathcal X_D=\{(g_1,\ldots,g_6): g_i\text{ is DT}\},
 \qquad
 \mathcal X_T=\{(g_1,\ldots,g_6):
          g_1\text{ is a 3-cycle},\ g_2,\ldots,g_6\text{ are DT}\}.
 \tag{1}
\]
All these tuples are transitive because their product is a full cycle.
Their finite indices sum to twelve, so their source genus is zero.
The first factor in \(\mathcal X_T\) is marked.

Simultaneous conjugation preserving \(P\) is exactly conjugation by
\(C_{S_{13}}(P)=\langle P\rangle\). These cyclic quotients are the
absolute polynomial domains. For \(\mathcal X_D\) use all \(Q_i\);
for \(\mathcal X_T\) use
\[
                 H=\langle Q_1^2,Q_2,Q_3,Q_4,Q_5\rangle.
 \tag{2}
\]
This is the full subgroup fixing the first branch label, as proved in
the [uniform marking reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md),
Lemma 1.

**Theorem 1.** The absolute domain \(\mathcal X_D/\langle P\rangle\)
is one orbit of size \(446265625\). The absolute marked domain
\(\mathcal X_T/\langle P\rangle\) is one \(H\)-orbit of size \(89253125\).
Every tuple in either domain generates \(A_{13}\).

The proof consists of the exact upper bounds and exhaustive orbit
algorithm below. Neither a sampled orbit nor an inequality in the wrong
direction would suffice.

## An independent exact total bound

For a partition \(\lambda\vdash13\), let \(m_d\) be its multiplicities and
\(z_\lambda=\prod_d d^{m_d}m_d!\). For the standard permutation
representation with the trivial summand removed, the character of its
\(k\)-th exterior power is
\[
 \chi_k(\lambda)=[t^k]\,
        \frac{\prod_d(1-(-t)^d)^{m_d}}{1+t},
                      \qquad 0\le k\le12.
 \tag{3}
\]
This follows directly from the eigenvalues of a permutation matrix.
The checker enumerates all 101 partitions of thirteen and verifies, over
the rationals,
\[
 \sum_{\lambda\vdash13}
       \frac{\chi_i(\lambda)\chi_j(\lambda)}{z_\lambda}
                         =\delta_{ij}.
 \tag{4}
\]
Thus these thirteen actual representation characters are distinct
irreducible characters. Here we use ordinary finite-group character
orthogonality: averaging on a Hom space computes the character inner
product, and complete reducibility gives the irreducibility criterion.

Their values on \(P\) are \((-1)^k\). The sum of their squared values is
thirteen, the centralizer order of \(P\). Column orthogonality therefore
shows that every remaining irreducible character vanishes at \(P\).
Consequently only these thirteen terms contribute to the class-product
count.

For completeness, if \(C_1,\ldots,C_r\) are conjugacy classes in a finite
group \(G\), their central class sums act on an irreducible representation
of dimension \(d_\chi\) as the scalars
\(|C_i|\chi(C_i)/d_\chi\). The coefficient of \(P\) in their product is
the regular-representation trace of \(P^{-1}\) times that product, divided
by \(|G|\). Decomposing the regular representation gives
\[
 \#\{g_1\cdots g_r=P:g_i\in C_i\}
 =\frac{\prod_i|C_i|}{|G|}
       \sum_\chi\frac{\chi(P^{-1})\prod_i\chi(C_i)}
                         {d_\chi^{r-1}}.
 \tag{5}
\]
This is also the Frobenius formula stated in
[Jones--Zvonkin, Section 3.2, printed pages 139--140](sources/jones-zvonkin-braid-cacti-2002.pdf).

In our case \(D=3\binom{13}4=2145\), \(T=2\binom{13}3=572\),
\(d_k=\binom{12}k\), and
\[
 \chi_k(D)=[t^k](1+t)^{10}(1-t)^2,\qquad
 \chi_k(T)=[t^k](1+t)^9(1+t^3).
\]
The complete character rows are
\[
\begin{split}
 d_k&=(1,12,66,220,495,792,924,792,495,220,66,12,1),\\
 \chi_k(D)&=(1,8,26,40,15,-48,-84,-48,15,40,26,8,1),\\
 \chi_k(T)&=(1,9,36,85,135,162,168,162,135,85,36,9,1).
\end{split}
\]
Substitution into (5) yields, exactly,
\[
 \frac{D^6}{13!}\sum_{k=0}^{12}
       \frac{(-1)^k\chi_k(D)^6}{d_k^5}=5801453125,
\quad
 \frac{TD^5}{13!}\sum_{k=0}^{12}
       \frac{(-1)^k\chi_k(T)\chi_k(D)^5}{d_k^5}=1160290625.
 \tag{6}
\]
The second count has its 3-cycle in the first position, without an extra
factor of six.

The cyclic quotient is free. A nonidentity power of a thirteen-cycle
has thirteen-element orbits, so cannot preserve the support of a DT.
A tuple fixed by that power would in particular have its DT factors
fixed, which is impossible. Dividing (6) by thirteen gives the two
absolute total bounds in Theorem 1.

These counts agree with the Goulden--Jackson cactus formula, but that
formula is not a proof premise here. Its original 1992 paper prints
\(\sum l(\lambda_i)=n+1\) in Theorems 3.1--3.2 for general numbers of
factors; the cactus argument and the last line of its proof use the
correct \((r-1)n+1\). Jones--Zvonkin Section 1.3 states the latter
planarity convention correctly. The independent character calculation
avoids relying on an unstated correction of that printed hypothesis.

## Collision-free orbit enumeration

The [checker](verify_degree_fourteen_hurwitz_orbit.py) retains both
literal C++ engine sources. It compiles them with assertions enabled
and runs the high-memory cases sequentially.

There are 2145 DTs. Their full permutations fit into thirteen nibbles
of a 64-bit word. A DT itself is reconstructed from its two unordered
support pairs. The engine enumerates all four-element supports and
their three matchings, and precomputes exact conjugation and cyclic
rotation tables.

For the DT domain, the first factor has 165 free cyclic orbits.
Normalize it to the smallest DT identifier in its orbit, and rotate
every other factor by the same unique power of \(P\). Store its
representative number in eight bits and the next four DT identifiers
in twelve bits each. The sixth factor is forced:
\[
                g_6=(g_1g_2g_3g_4g_5)^{-1}P.
\]
Thus the normalized tuple has a 56-bit injective key. A full-permutation
assertion checks that the reconstructed last factor is indeed a DT.

For the marked domain the first factor has 44 free cyclic orbits.
Its representative needs six bits; the next four DTs again need
twelve each. The last factor is forced, giving a 54-bit injective key.
The smaller degree-nine control has a first 3-cycle with a possible
cyclic stabilizer. Its separate encoder minimizes the complete stored
tuple over *every* shift giving the minimal first factor; it does not
assume a unique first-factor shift in that control.

The DT neighbors are
\((a,b)\mapsto(b,bab)\). In the marked first pair they are
\[
 a'=bab,\qquad b'=(a')^{-1}ba',
\]
which is exactly \(Q_1^2\). The other marked neighbors are ordinary DT
moves. All tables use literal permutations and preserve the product.
After each neighbor, normalize by the cyclic action and insert the
complete key into an open-addressed hash table. Comparisons use the
entire key; hash collisions are resolved by exact equality.

The queue starts with a known full-alternating seed. In the DT case use
the cyclic seed from the uniform marking reduction with \(r=6\):
\(g_j=P^{2j}(02)(34)P^{-2j}\), \(0\le j<6\).
For \(e=g_0(g_0g_1)^2\),
\[
                   e=(04)(23),\qquad (eg_2)^2=(064).
\]
The odd-degree support argument gives primitivity, and Jordan's
3-cycle theorem gives \(A_{13}\).
For the marked case start with \((012)\), then append
\((0\,\,2j+1)(1\,\,2j+2)\), \(1\le j\le5\), and relabel its full-cycle
product to \(P\). Its expanded graph is a tree; it is transitive,
primitive by the same support argument, and already contains a
3-cycle. It too generates \(A_{13}\).

The finite queue is exhausted, rather than stopped when its observed
size first reaches the expected bound. Only positive moves are needed:
on a finite set each elementary braid is a bijection, so its inverse
is a positive power. The visited sizes are respectively exactly the
two bounds following (6). Hence each queue contains its entire
absolute domain. Generated groups are preserved by Hurwitz moves,
and cyclic conjugation preserves \(A_{13}\), proving the last
assertion of Theorem 1.

The retained control censuses have sizes

| Polynomial degree | DT absolute orbit | Marked 3-cycle absolute orbit |
|---:|---:|---:|
| 7 | 7, inside a total domain of 56 | 28 |
| 9 | 6561 | 2187 |
| 11 | 1362944 | 340736 |
| 13 | 446265625 | 89253125 |

The degree-seven DT seed belongs to a proper order-168 component.
This deliberate exception prevents treating a single seed orbit as
the full domain without comparing its size with an independent bound.
The checker also counts all 392 degree-seven DT triples directly.

## Transfer to the degree-fourteen passport

**Theorem 2.** The passport
\(\bigl((2,2,1^{10})^7,(13,1)\bigr)\) has one absolute transitive
Hurwitz orbit, and every such tuple has group \(A_{14}\).

**Proof.** Apply the
[uniform cycle-folding theorem](ALTERNATING_HURWITZ_CYCLE_FOLDING.md).
Every original tuple braids to the adjacent doubled-edge form with
the pole letter isolated except in that pair. Coalescing the adjacent
pair gives a degree-thirteen polynomial tuple of one of the two types
(1), with its coalesced branch marked and one outside attachment.

The marked-3-cycle domain is connected by Theorem 1. For the DT domain,
the cyclic seed has a literal stabilizer inducing a six-cycle on
the six branch labels: \((Q_1\cdots Q_5)^{13}\).
The uniform unmarked-to-marked lemma therefore upgrades its unmarked
connectivity to connectivity with the coalesced branch marked.

These assertions can also be used literally at fixed \(P\), as required
by the attachment theorem: the central full twist fixes every branch
label and effects simultaneous conjugation by \(P^{-1}\). It removes
any residual cyclic relabeling in a path between absolute representatives.

Corollary 6 of the folding note gives transitivity on all outside
attachments in both polynomial types. Thus all expansions of each
type lie in one original orbit. Its uniform cross-type braid link
(Corollary 4) joins those two expansion types. The complete
coalescence-and-expansion reduction (Corollary 5) now gives exactly one
absolute orbit of the original passport.

Every coalesced group is \(A_{13}\) by Theorem 1 and fixes the isolated
pole letter. The original group contains it and moves that letter.
The point stabilizer \(A_{13}\) is maximal in the natural \(A_{14}\):
any larger subgroup is transitive and contains a point stabilizer,
so its order is \(14\cdot|A_{13}|=|A_{14}|\).
All original factors are even. The original group is therefore
\(A_{14}\). This also follows from the single-orbit conclusion and
the retained alternating generating tuple. ∎

The theorem classifies a finite Nielsen passport. It supplies no
rational coefficients, rational point or unirational parametrization
of its Hurwitz component, and no braid word identifying an unspecified
algebraic polynomial with a particular tuple. Connectivity in degree
sixteen and uniform connectivity remain open.

For the residual pure-DT polynomial problem on \(N\ge15\) sheets,
the [global compatibility restrictions](TERMINAL_MASS_HURWITZ_DESCENT.md#corollary-11-positive-minima-have-restricted-mass-and-terminal-sizes)
force a hypothetical positive orbit minimum to satisfy \(3m\ge N-4\),
with equal terminal sizes or a singleton at a minimizing factor.
Full alternating monodromy has not yet been shown to contradict these
remaining cases. These restrictions neither extend the finite census above
nor resolve the separate mixed polynomial problem needed for degree sixteen.

## Evidence

The [foreground census record](degree_fourteen_hurwitz_census_inputs.json)
and [receipt](degree_fourteen_hurwitz_verification.json) pin the exact
character arithmetic, literal engine sources, compiler, execution
commands, observed orbit counts and memory scope. The original full
foreground DT census took 343.665 seconds; the marked census took
44.2615 seconds and measured peak RSS 2812372 KiB. The DT engine
allocates an eight-GiB hash table and reserves a four-GiB queue;
its observed RSS during the run was 10502376 KiB, and its peak RSS
was not instrumented in that first run. These are execution records,
not asymptotic complexity claims.

The Python wrapper and engines can repeat the complete calculation.
Its explicit collection mode binds already completed foreground engine
outputs to their exact source and executable hashes while rerunning
the character and smaller controls; it does not claim a second large
census. The default mode reruns both large engines.
The [source ledger](degree_fourteen_hurwitz_source_ledger.json)
distinguishes this finite proof from its uniform imports. No novelty,
external human review, independent external replay or formal
verification is asserted.
