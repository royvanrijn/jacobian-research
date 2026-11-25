# Sharpness of the absolute Jacobi bound on primitive balanced boxes

Fix three finite lists of positive root caps \(p_{ij}\) and positive
integral powers \(e_{ij}\), with \(\sum_jp_{ij}=1\) for each letter
\(i\in\{A,B,C\}\). Put \(W_i=\sum_jp_{ij}e_{ij}\) and assume
\(W_C\ge\max(W_A,W_B)\). Fix \(\beta>2\), set \(Q=N^\beta\),
and fix \(M\ge\max_{i,j}e_{ij}\). All data are fixed before
\(N\to\infty\). There exist factorized dictionaries of sizes
\(D_i=N^{1+o(1)}\) of the form
\[
       U_i=\prod_jx_{ij}^{e_{ij}},\qquad x_{ij}\in A_{ij},
\]
with pairwise disjoint prime pools for the original slots, such that
all tuples \(a,b,c\) are pairwise coprime, have \(c>\max(a,b)\),
and satisfy the exact combined cap
\(\operatorname{rad}(abc)<c^{3/W_C}\) and
\(H_M(c)=\operatorname{rad}(c)^M\). There is no true equality \(a+b=c\).
For a positive proportion of the primes \(q\in[Q,2Q]\), their
absolute high-order Jacobi term has size
\[
 \frac{1}{(q-1)^2}
 \sum_{\substack{\chi,\psi\\
                  \chi\in\mathcal G_A,\ \psi\in\mathcal G_B,\
                  \chi\psi\in\mathcal G_C}}
  |J_q(\chi,\psi)\widehat w_A(\chi)\widehat w_B(\psi)
                         \overline{\widehat w_C(\chi\psi)}|
       \asymp \sqrt q\, (D_A D_B D_C)^{1/2}.
 \tag{AJ1}
\]
Here \(\mathcal G_i\) excludes all characters whose order divides
any power \(e_{ij}\) on letter \(i\), and \(w_i\) counts exact
integer labels on unit residues. The corresponding prime average of the absolute value
of the signed sum is \(O(N^{2-\beta/2+o(1)})\).

Choose \(0<s<1/(100\max_i W_i)\). Before thinning, choose
the prime-root pools \(P_{ij}\) so that
\(x^{1/p_{ij}}/N\) lies in a fixed closed subinterval of
\((1,1+s)\) for \(A,B\) and of \((1+2s,1+3s)\) for \(C\).
Choose disjoint subintervals for slots with an equal cap,
also separating the two letter intervals for \(A,B\).
Different caps give disjoint prime intervals for large \(N\).
Each pool has size \(m_{ij}\asymp N^{p_{ij}}/\log N\)
by the prime number theorem in fixed-ratio intervals.
All resulting products are injective as functions of
the root choices, and all tuples are pairwise coprime.
The radicals of \(A,B\) lie in \((N,(1+s)N)\), and that of
\(C\) lies in \(((1+2s)N,(1+3s)N)\). The numerical height
of each value is between the corresponding range endpoints
raised to \(W_i\), and the minimum repair is the asserted one.
Since
\[
   (1+s)^2(1+3s)<(1+2s)^3,
\]
the exact combined radical budget holds. If both summands
have \(W_i=W_C\), then \(a+b>c\) since
\((1+3s)^{W_C}<2\). If at least one has smaller \(W_i\), then
\(a+b<c\) for all large \(N\): at most one summand has the
same height exponent as \(C\), and this one is uniformly smaller
by a fixed factor. Also \(c>\max(a,b)\) in both cases.

Each full, unthinned letter has at most \(O(N)\) labels of
fixed polynomial height. Any nonzero difference between two
labels is divisible by at most \(O(1)\) primes in \([Q,2Q]\).
Thus the residue map for all three full label populations is
injective at all but \(O(N^2)\) test primes, a vanishing
proportion when \(\beta>2\). These primes and their injectivity
are independent of the thinning.

Now independently keep each prime in each slot with probability
one half, and call the retained pools \(A_{ij}\).
If \(z_1,\ldots,z_m\) are arbitrary unit complex numbers and
\(\xi_j\) are independent Bernoulli variables with parameter
one half, then
\[
            \mathbb E\left|\sum_j\xi_jz_j\right|
                          \gg \sqrt m.                     \tag{AJ2}
\]
Indeed, let \(\xi_j'\) be independent copies.
The sum \(Z=\sum_j(\xi_j-\xi_j')z_j\) has
\(\mathbb E|Z|^2=m/2\) and \(\mathbb E|Z|^4\ll m^2\):
a nonzero term in the fourth-moment expansion uses each
index at least twice. Hölder implies \(\mathbb E|Z|\gg\sqrt m\),
and the triangle inequality gives (AJ2).
The phase choices for a slot in the Fourier transform are
\(z_x=\overline{\chi^{e_{ij}}(x)}\), all of unit modulus.
Factorization and independence across slots and letters show
that the expected product of the three transform magnitudes
is \(\gg \prod_{i,j}\sqrt{m_{ij}}\), uniformly in all mode
choices. There are \((q-1)^2-O(q)\) high-order character pairs
and on all these pairs \(|J_q|=\sqrt q\). The expected
absolute mass in (AJ1) is therefore
\(\gg \sqrt q\prod_{i,j}\sqrt{m_{ij}}\).

On an injective test prime, Parseval gives
\(\sum_\chi|\widehat w_i(\chi)|^2=(q-1)D_i\).
Successive Cauchy inequalities bound the absolute mass in
(AJ1) by \(\sqrt q(D_A D_B D_C)^{1/2}\).
This is at most \(\sqrt q\prod_{i,j}\sqrt{m_{ij}}\) for all
thinnings. Thus a fixed positive proportion of thinnings
reach a fixed positive fraction of this latter upper bound
for each injective prime. The event that some pool loses
more than three quarters of its primes has probability
\(O(\sum_{i,j}\exp(-c m_{ij}))=o(1)\). It follows by averaging
over primes and good thinnings that a deterministic choice
of the pools, with \(D_i\asymp\prod_j m_{ij}\) for every
letter, reaches (AJ1) at a positive proportion of test primes.

Finally a false tuple contributes to at most \(O(1)\)
prime tests, so its weighted modular count averages
\(O(N^{3-\beta+o(1)})\). The principal term has the same
bound and the union of bounded-order terms contributes at
most \(O(N^{2-\beta/2+o(1)})\) by residue collision averaging
and Cauchy. Subtract these terms in (MF1) to obtain the
claimed upper bound by the triangle inequality on the prime average.
[Jacobi decomposition and bounded-order estimate](ABC_MODULAR_JACOBI_INCIDENCE.md#jacobi-sums-and-the-removal-of-bounded-order-modes)


## Exact resummation in one character coordinate

The same construction, possibly with a different deterministic
thinning, gives a lower bound even after performing the inner
character sum exactly. For dictionaries on unit residues \(B,C\)
and a character extended by zero at zero, orthogonality gives
\[
 S_{\chi}=\frac1{q-1}\sum_{\psi}J_q(\chi,\psi)
                  \widehat w_B(\psi)
                  \overline{\widehat w_C(\chi\psi)}
           =\sum_{b,c} w_B(b)w_C(c)\chi(c-b).
 \tag{AJ3}
\]
Indeed, in the Jacobi definition the \(\psi\)-sum forces
\((1-u)c=b\), and the remaining \(\chi\)-phase is
\(\chi(u)\chi(c)=\chi(c-b)\). If \(c=b\) the term vanishes.

Let \(L_q\) consist of the average of
\(|\widehat w_A(\chi) S_\chi|\) over \(\chi\in\mathcal G_A\),
with normalization \(1/(q-1)\). At a positive proportion of
test primes the construction can be thinned so that
\[
                 L_q \gg (D_A D_B D_C)^{1/2}.               \tag{AJ4}
\]
Also this holds if the inner sum over \(\psi\) in \(S_\chi\)
is restricted to the other two good modes, with the same
absolute value on the outside.

To prove this, discard also the \(O(N^2)\) primes dividing
a nonzero difference \(c-b\) of full unthinned labels,
so all the phases in (AJ3) have modulus one. All these
differences are nonzero because \(c>b\).
Set \(k=\sum_i\#\{j\}\).
For a fixed good \(\chi\), the product
\(\widehat w_A(\chi)S_\chi\) is a decoupled multilinear
Bernoulli polynomial \(T_\chi(\xi)\), with one slot
variable from each of the \(k\) independent pools.
Its coefficients indexed by choices of one prime per
slot all have modulus one. Write \(\xi=(1+\sigma)/2\),
where the \(\sigma\)'s are independent Rademacher signs.
The component of full degree \(k\) then shows by
orthogonality that
\(\|T_\chi\|_{L^2(\xi)}
                 \ge 2^{-k}\prod_{i,j}\sqrt{m_{ij}}\).

For every such multilinear Bernoulli polynomial \(T\),
even with arbitrary complex coefficients, one has
\(\|T\|_{L^4(\xi)}\ll_k\|T\|_{L^2(\xi)}\).
Expand by the subset of slots whose signs occur, an
orthogonal decomposition in \(L^2\).
For each decoupled sign component, successively apply
the elementary Rademacher fourth-moment inequality
\(\|\sum_x z_x\sigma_x\|_4\le3^{1/4}
(\sum_x|z_x|^2)^{1/2}\) in each slot, with Minkowski to
move the \(\ell^2\) sum outside \(L^4\). This bounds
its \(L^4\) norm by \(3^{k/4}\) times its \(L^2\) norm.
Summing the at most \(2^k\) components, with Cauchy,
proves the claimed inequality. Interpolation therefore
gives \(\|T\|_2\ll_k\mathbb E|T|\).
In particular
\(\mathbb E_\xi L_q\gg_k\prod_{i,j}\sqrt{m_{ij}}\)
and, again by Minkowski,
\(\|L_q\|_{L^2(\xi)}\ll_k\mathbb E_\xi L_q\).
The latter bound and Cauchy give a fixed positive
probability that \(L_q\) exceeds half its mean, at
each of the retained test primes. Exponentially
unlikely pool compressions cost \(o(1)\) of this
probability. Average over primes to fix a good
deterministic thinning satisfying (AJ4) at a
positive proportion of primes.

To restrict the inner sum in (AJ3) to the other
two good modes, use the triangle inequality on the
outside. The change in (AJ4)'s left side is at most the
absolute sum of the deleted modes after Jacobi
normalization. By (MF2), pointwise on the injective
test primes, this is \(O(N^{2-\beta/2+o(1)})\),
strictly smaller than the right side of (AJ4).
Thus fully resumming one character coordinate still
leaves an absolute sum too large for a power saving
below \(N\). The corollary concerns these two
specific absolute-value operations and leaves signed
one-dimensional estimates open.
