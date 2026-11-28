# Transferring Gaussian stability to balanced voting

Let \(k\ge 2\). A vote profile \(X\) is uniform on \([k]^n\).
Construct \(Y\) by independently retaining each entry with
probability \(\rho\) and replacing it by a fresh uniform label
with probability \(1-\rho\). For a rule \(f:[k]^n\to[k]\), its
stability is \(\mathbb P[f(X)=f(Y)]\). Balance means
\(\mathbb P[f(X)=j]=1/k\) for every \(j\).

At a fixed positive correlation, Gaussian simplex optimality
in every dimension implies plurality optimality as the maximum
voter influence tends to zero. Conversely, any fixed balanced
Gaussian partition has one sequence of deterministic, exactly
balanced, low-influence implementations whose entire transition
matrices converge uniformly for \(0\le\rho\le1\). In particular a
strict Gaussian counterexample at any positive correlation would
give a deterministic voting counterexample at that correlation.

## Gaussian optimality implies the voting inequality

For \(0<\rho<1\) and \(\epsilon>0\), Isaksson and Mossel's
Theorem 7.1 supplies a degree \(D\) and \(\tau>0\) such that if
all degree-\(D\) coordinate influences of the label indicators
of \(f\) are at most \(\tau\), then an ordinary Gaussian
partition in dimension \(n(k-1)\) has the same label masses
and stability within \(\epsilon\) of that of \(f\). The
degree-\(D\) influences are bounded above by the usual
conditional-variance influences. Writing \(X^{(i)}\) for the
profile with vote \(i\) independently resampled, the influence
of the \(j\)-indicator is
\[
 \frac12\mathbb E
 \left(1_{\{f(X)=j\}}-1_{\{f(X^{(i)})=j\}}\right)^2
 \le \frac12\mathbb P[f(X)\ne f(X^{(i)})].
\]
Thus a small uniform bound on the probability that changing
one vote changes the winner suffices. For a balanced rule the
transferred partition is balanced. Assuming the simplex maximizes
at this \(k,\rho\) in every admissible dimension, we obtain
\[
 \mathbb P[f(X)=f(Y)]\le S_\rho(T_k)+\epsilon,               \tag{VT1}
\]
where \(T_k\) is the regular simplex partition, whose stability
is unchanged by adjoining independent Gaussian coordinates.
The proof of (VT1) uses the Gaussian maximal value. It makes
no uniqueness assertion for discrete rules.
[Isaksson–Mossel, Theorem 7.1 and proof of Theorem 1.10](https://arxiv.org/pdf/0903.3362v3)

For an exactly balanced deterministic plurality rule, among
tied leaders choose the candidate whose first appearance in
the voter sequence is earliest. Every leader appears when
\(n\ge1\). The rule commutes with every candidate permutation,
so each candidate wins with probability \(1/k\). A single
resampling changes the winner only when two original counts
differ by at most two. For candidates \(a,b\), condition on
\(M=N_a+N_b\), where \(N_j\) is the count of candidate \(j\).
Then \(M\) is binomial with parameters \(n,2/k\), and
\(N_a-N_b\) conditionally has law \(2\operatorname{Bin}(M,1/2)-M\).
The central binomial bound and a binomial tail bound show
\(\mathbb P(|N_a-N_b|\le2)=O_k(n^{-1/2})\). A union bound
over candidate pairs gives the same order for the maximum
voter influence.

For the value of this rule, use the simplex vectors
\(v_j=\sqrt k(e_j-k^{-1}\mathbf1)\) in the hyperplane
\(\mathbf1^\perp\subset\mathbb R^k\). With a uniform vote
\(U\), the mean is zero and the covariance of \(v_U\) is
the identity on that hyperplane. For a noisy copy \(V\),
the cross covariance is \(\rho I\). The two normalized
vote-vector sums converge jointly to the corresponding
correlated Gaussian pair by the multivariate central limit
theorem. Away from ties, plurality takes the argmax of the
coordinates of each sum. The Gaussian comparison hyperplanes
have measure zero for both marginals. This proves
\[
 \lim_{n\to\infty}\mathbb P[\operatorname{Plur}_{k,n}(X)
                =\operatorname{Plur}_{k,n}(Y)]=S_\rho(T_k).
\]
Together with (VT1), this gives precisely the asymptotic voting
comparison at each positive correlation where balanced Gaussian
optimality is known.
[Isaksson–Mossel, Lemma 7.3](https://arxiv.org/pdf/0903.3362v3)

## A single deterministic implementation for all correlations

Fix a balanced measurable Gaussian partition \(A=(A_1,\dots,A_k)\)
of \(\mathbb R^d\). There are deterministic rules \(f_\ell\),
each exactly balanced and each with finite voter count, such
that their maximum outcome-change influences tend to zero and
\[
 \max_{a,b}\ \sup_{0\le\rho\le1}
 \left|\mathbb P[f_\ell(X)=a,\ f_\ell(Y)=b]
              -\mathbb P[G\in A_a,\ H\in A_b]\right|
                  \longrightarrow 0.                     \tag{VT2}
\]
Here \(G,H\) are standard Gaussian vectors with cross
covariance \(\rho I\). Isaksson–Mossel's reverse transfer,
Theorem 7.4, concerns one correlation and allows a
simplex-valued discrete rule. The construction below uses
deterministic labels and treats the entire transition matrix
simultaneously.
[Isaksson–Mossel, Theorem 7.4](https://arxiv.org/pdf/0903.3362v3)

Approximate \(A\) in Gaussian disagreement measure by a
partition \(B\) defined by finitely many coordinate threshold
cuts. Finite unions of boxes approximate each cell in Gaussian
measure, and resolving overlaps and uncovered points gives
such a partition with arbitrarily small disagreement
\(\eta=\gamma_d(\{A\ne B\})\). Label mismatch at either
argument of a pair changes each transition probability by
at most \(2\eta\), uniformly in the correlation.

Use \(d\) disjoint blocks of \(m\) votes. In each block,
take the normalized sum of the bounded scalar vote
contrast
\[
             b(u)=\frac{k1_{\{u=1\}}-1}{\sqrt{k-1}} .
\]
It has mean zero and variance one under a uniform vote,
and cross covariance \(\rho\) with the noisy copy. Apply
\(B\) to the vector of the \(d\) block sums to obtain a
deterministic preliminary rule \(F_m\). By the paired
central limit theorem, \((F_m(X),F_m(Y))\) converges in law
to \((B(G),B(H))\) at each \(\rho\in[0,1]\). This includes
\(\rho=1\), since the limiting pair is then supported on the
diagonal but each coordinate threshold still has marginal
Gaussian measure zero. The maximum voter influence of
\(F_m\) tends to zero. A resampled vote changes one block
sum by at most \(k/\sqrt{(k-1)m}\), so a change of label
requires that this block sum lie within that distance of
one of its finitely many thresholds. The marginal central
limit theorem and continuity of the Gaussian distribution
make these probabilities tend to zero, uniformly over the
voter positions.

The transition matrices for this fixed \(B\) converge
uniformly in \(\rho\in[0,1]\). For the indicator \(h\) of a
single label or a union of two labels, the expectation
\(\mathbb E[h(X)h(Y)]\) is increasing in \(\rho\). Indeed the
orthogonal-degree expansion for the product noise operator
has nonnegative squared-coefficient weights. The Gaussian
expectation has the same monotonicity from
its Hermite expansion, and it is continuous on \([0,1]\)
by summability of the squared coefficients. Pointwise
convergence of increasing functions to a continuous limit
on this compact interval gives uniform convergence by
squeezing between successive points of a finite mesh.
Finally the mixed entry \(a,b\), for \(a\ne b\), is
one half of the difference between the expectation for
the union and the sum of the expectations for the individual
labels, by symmetry of the two profiles. This proves the
uniform claim.

To repair the masses, move individual profiles from
overrepresented labels to underrepresented labels. Each
target count is the integer \(k^{dm-1}\), so the exact
fraction of changed profiles can be taken to be
\[
                     \delta_m=\frac12\sum_{j=1}^k
               \left|\mathbb P[F_m(X)=j]-1/k\right|.
\]
For fixed \(B\), the marginal central limit theorem gives
\(\limsup_m\delta_m\le\eta\). Relabeling this fraction changes
each transition probability by at most \(2\delta_m\)
uniformly in \(\rho\) and each outcome-change influence by
at most \(2\delta_m\), since the noisy or singly resampled
profile has the same uniform marginal. Taking box
approximations with \(\eta\to0\) and then choosing \(m\)
large for each one proves (VT2). In particular the sum of
the diagonal entries, the agreement probability, converges
uniformly to \(S_\rho(A)\).
