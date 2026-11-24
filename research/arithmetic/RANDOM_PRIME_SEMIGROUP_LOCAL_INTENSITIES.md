# Local intensities and critical gaps for random prime semigroups

Use the [multiplier model](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md)
with \(2\) admitted, independent selection of each odd prime \(p\) with
probability \(\min(1,\lambda p^{-\gamma})\), and unordered representations
\(r_k(n)\) with repetitions. Write \(\alpha=1-\gamma\), \(k\ge2\), and let \(A>0\)
be the almost-sure constant in the semigroup counting asymptotic. All
probability-one assertions below are for fixed model parameters. Haar
expectations refer to independent auxiliary variables with the selected
prime set held fixed.

## Local convolution and the parity obstruction

For a prime \(p\), put \(t_p=p^{-\alpha}\) if selected and \(t_p=0\)
otherwise. Set \(b_p=(pt_p-1)/(p-1)\) and \(\theta_p=pt_p^k\).
For Haar-almost every \(u\in\mathbb Z_p\), with \(v=v_p(u)<\infty\), define
\[
 g_p(u)=1+b_p^k\left((p-1)\sum_{j=0}^{v-1}\theta_p^j-\theta_p^v\right).
 \tag{1}
\]
Use an empty sum at \(v=0\) and \(0^0=1\). Thus for an excluded prime \(g_p\)
is \(1-b_p^k\) on the units and \(1+(p-1)b_p^k\) on \(p\mathbb Z_p\).
For a fixed modulus \(m\), let
\[
 S_m(a)=\prod_{p^f\parallel m}
       \mathbb E\!\left[g_p(U_p)\mid U_p\equiv a\bmod p^f\right],
 \qquad
 B=\frac{\Gamma(1+\alpha)^k}{k!\,\Gamma(1+k\alpha)}A^k,
\]
where the \(U_p\) are independent Haar-uniform variables. Then
\[
 \frac{m}{X^{k\alpha}(\log X)^{k(\lambda-1)}}
        \sum_{\substack{n\le X\\n\equiv a\bmod m}}r_k(n)
 \ \longrightarrow\ B S_m(a) \quad\text{almost surely}.                 \tag{2}
\]

To evaluate the local convolution in [equation (4) of the multiplier
theorem](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md#counting-and-residue-classes),
use the \(p\)-adic variable with valuation \(j\) of probability
\((1-t_p)t_p^j\) and uniform unit part. For an additive character of exact
conductor \(p^\ell\), \(\ell\ge1\), its mean on this variable is
\[
 t_p^\ell-\frac{(1-t_p)t_p^{\ell-1}}{p-1}=b_p t_p^{\ell-1}.
\]
Indeed the unit average is \(1\) at valuations \(j\ge\ell\), is
\(-1/(p-1)\) at \(j=\ell-1\), and vanishes for \(j<\ell-1\) by subtracting
the multiples of \(p\) from the full residue sum. The \(k\)-fold convolution
raises these means to the \(k\)-th power. Inverse Fourier transform on
\(\mathbb Z/p^f\mathbb Z\), with densities relative to uniform measure, gives
partial sums of
\[
 1+\sum_{\ell\ge1} b_p^k t_p^{k(\ell-1)}c_{p^\ell}(u),
 \qquad
 c_{p^\ell}(u)=p^\ell\mathbf1_{p^\ell\mid u}
              -p^{\ell-1}\mathbf1_{p^{\ell-1}\mid u}.
\]
Each term has Haar \(L^1\) norm at most \(2|b_p|^kt_p^{k(\ell-1)}\).
Consequently the series converges in \(L^1\), integrates to \(1\), and its
conditional expectation modulo \(p^f\) is the finite Fourier density.
At \(v_p(u)=v<\infty\) it is the finite sum (1). The Chinese remainder theorem
now gives (2) from the fixed-modulus averaged representation law.

In particular \(S_2(a)=1+(-1)^a b_2^k\) with \(b_2=2^{1-\alpha}-1>0\).
The normalized difference between the even and odd cumulative counts in
(2), without the prefactor \(m\), tends to \(B b_2^k\). The sequence
\(n^{k\alpha-1}(\log n)^{k(\lambda-1)}\), for \(n\ge2\), has asymptotically
equal normalized sums over the two parity classes. This follows from
sum–integral comparison for regularly varying powers of \(n\) and \(\log
n\), or by pairing \(2j\) with \(2j+1\) and using the ratio tending to one.
Thus for every finite real \(c\), even chosen from the realized prime set,
\[
 \liminf_{X\to\infty}
 \frac{\sum_{2\le n\le X}|r_k(n)-c n^{k\alpha-1}(\log n)^{k(\lambda-1)}|}
      {X^{k\alpha}(\log X)^{k(\lambda-1)}}\ \ge B b_2^k.                  \tag{3}
\]
This follows by applying the triangle inequality to the signed sum with
sign \((-1)^n\). In particular a single finite limit for
\(r_k(n)/(n^{k\alpha-1}(\log n)^{k(\lambda-1)})\) is impossible. Formula (3)
applies in particular to the subcritical pointwise representation node
\(\texttt{GAP-ARITH-RANDOM-PRIME-POINTWISE-REPRESENTATIONS}\).

## Stronger almost-sure critical gap bound

Now take \(\lambda=1\) and \(\alpha=1/k\). In this case \(B=\kappa_k A^k\)
with \(\kappa_k=\Gamma(1+1/k)^k/k!\). Fix any modulus \(m\) and residue \(a\).
Almost surely,
\[
 \liminf_{X\to\infty}\frac{m}{X}
       |\{n\le X:n\equiv a\bmod m,\ n\notin kC\}|
       \ \ge\ \exp(-B S_m(a)).                                        \tag{4}
\]

Here is the residue-class extension of the conditional concentration
argument in [Gaps below the logarithmic transition, equations (8)–(10)](RANDOM_PRIME_SEMIGROUP_LOGARITHMIC_TRANSITIONS.md#gaps-below-the-logarithmic-transition).
Use \(Y=\exp((\log X)^b)\) for \(0<b<1\), and let \(\mathcal F_Y\) expose
prime choices through \(Y\). Use the two observables in that argument,
restricting the gap count and representation total to \(n\equiv a\bmod m\)
and normalizing each by \(m/X\). Switching on a prime \(q>Y\) changes each
observable by at most
\[
 \frac{mk}{X}N_{C^+}(2X/q)N_{C^+}(2X)^{k-1},
\]
where \(C^+\) has \(q\) forced in. The restricted counts change
monotonically and every added representation has a member divisible by
\(q\). By Hölder and the uniform forced-prime \(2k\)-th moment bound in
(8), the squared influence has mean \(O_{m,k}(q^{-2\alpha})\). Applying
the variance argument above (9) gives expected conditional variance
\(O_{m,k}(\sum_{q>Y}p_q q^{-2\alpha})=O_{m,k}(Y^{-\alpha})\). Thus along
each fixed geometric grid both observables differ from their conditional
means by \(o(1)\) almost surely, by Chebyshev and Borel–Cantelli. The
conditional representation mean tends to \(B S_m(a)\) by (2).

On that grid the number of targets represented using only exposed primes
is \(o(X)\), by the Rankin bound preceding (9). Within the class under
consideration this is \(o(X/m)\) since \(m\) is fixed. On every other
target, each possible representation has conditional probability at most
\(Y^{-\gamma}\). The event that a particular representation fails is
decreasing in the unexposed Bernoulli variables. Harris' inequality bounds
the probability of a gap below by the product of the failure probabilities,
which is at least \(\exp(-\mu_n/(1-Y^{-\gamma}))\) for
\(\mu_n=\mathbb E[r_k(n)\mid\mathcal F_Y]\).
Jensen averaging over the undiscarded targets yields (4) along the grid.
The unnormalized residue-class gap count is monotone. Comparing with the
preceding grid point and then taking a countable set of grid ratios
decreasing to \(1\) proves (4). Intersect over the countably many classes.

At a selected prime \(\theta_p=1\), so \(g_p=1+b_p^k((p-1)v_p-1)\).
For an excluded prime use the two values following (1). All factors are
strictly positive Haar-almost surely, since \(p=2\) is admitted and
\(|b_p|<1\) for all the remaining cases. Each has mean \(1\), and
\[
 \mathbb E|g_p(U_p)-1|=2(1-1/p)|b_p|^k.
\]
For selected primes \(|b_p|^k\ll_k p^{-1}\), with
\(\sum_{p\ \mathrm{selected}}p^{-1}<\infty\) almost surely since its
expectation over odd primes is \(\sum_p p^{-1-\gamma}<\infty\).
For excluded primes the sum is bounded by \(\sum_{p\ge3}(p-1)^{-k}\).
It follows that \(S=\prod_p g_p(U_p)\) exists and is positive and finite
Haar-almost surely. For fixed selected primes the finite products are
Cauchy in \(L^1\) by independence, telescoping, and the mean-one property.
In particular \(\mathbb E S=1\). Also
\(\mathbb E[S\mid U\bmod m]=S_m(U\bmod m)\) by independence and \(L^1\)
convergence. Along a divisible sequence of moduli generating the full
product sigma field, these conditional expectations converge in \(L^1\)
to \(S\).

Sum (4) over the classes and then pass along such a sequence using the
bounded Lipschitz function \(z\mapsto e^{-Bz}\) on \(z\ge0\). This proves
\[
 \underline d(\mathbb N\setminus kC)
 \ \ge\ \mathbb E_U e^{-B S(U)}
 \ \ge\ e^{-B}\cosh(B b_2^k)
 \ >\ e^{-B}.                                                         \tag{5}
\]
The second inequality is conditional Jensen using \(U\bmod 2\), and the
strict inequality uses \(B,b_2>0\). This strengthens the lower bound in
\(\texttt{GAP-ARITH-RANDOM-PRIME-CRITICAL-NATURAL-DENSITY}\).
Existence of the natural density and its law remain open.

## Transfer from finite-group subset sums

The [pinned Ma–Tang resolution](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/ma-tang-random-subset-sums/PROOF.md)
addresses [Problem 4.1 in v2](https://arxiv.org/html/2602.05768v2#S4).
It concerns the limsup over all group orders and uses \(G=\mathbb F_2^d
\times\mathbb F_3\). Its elementary certificate can be verified as
follows. A surjective binary coordinate map has fibres \(u+K\) for
\(K=\ker A\). Parameterize a dimension-\(r\) kernel by binary columns
\(\ell_j\) so that the \(j\)-th coordinate of its point at \(t\) is
\(\ell_j\cdot t\). For ternary marks \(c_j\), the nonconstant Walsh
coefficient at \(b\ne0\) of the ternary coordinate on \(u+K\) is
\[
 -\tfrac12\sum_{\ell_j=b}c_j(-1)^{u_j}.
\]
This follows by writing each binary bit as \((1-(-1)^{u_j+\ell_j\cdot t})/2\)
in \(\mathbb F_3\). All coefficients vanish whenever every nonzero
frequency bin has zero or at least two nonzero marks: two can be
paired with opposite signs, and an odd number at least three can use
three equal signs and pair the rest. Then one entire binary fibre
hits a single ternary value. The independently sampled column model
has exactly-one probability at most
\((2k/3)e^{-2(k-1)/(3\cdot 2^r)}\) by the union bound. At
\(k=d+r\) and \(2^r\asymp d/(10\log d)\), with \(2^r\le d/(10\log d)\),
this tends to zero. The rank and distinctness conditioning costs tend
to zero by the elementary union bounds \(2^{-d}\) for the random kernel
embedding and \(O(k^2/|G|)\) for duplicate samples. This verifies
the elementary lower-bound mechanism. The bounded-remainder variant
adds a characteristic-three additive-basis input.
[Pinned supplement](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/ma-tang-random-subset-sums/BOUNDED-REMAINDER.md)

The prime-semigroup problem retains residue-dependent fibres, for
which (1)–(5) quantify the effect on gap bounds. A representation of
a fixed target by members of \(C\) does not supply a binary affine
fibre with freely selectable signs. The finite model's cancellation
cannot be used to supply the missing pair-overlap estimates below.

## A sufficient condition for the critical density

Ma–Tang's [arXiv v1, §§3.1–3.3, Proposition 3.1 and Lemma 3.2](https://arxiv.org/html/2602.05768v1#S3)
uses addition of independent uniform variables in a prime field. An \(r\)-tuple
of representation equations with coefficient matrix \(M\) has joint
probability \(p^{-\operatorname{rank} M}\) for a consistent right-hand
side, by counting a coset of \(\ker M\), and zero otherwise. Bonferroni then
turns the resulting factorial-moment estimates into an estimate for the
probability of no representation. In the prime-semigroup model, if
representations require remaining prime supports \(T_i\), the joint
probability is exactly \(\prod_{q\in\cup_i T_i}p_q\). Row rank does not
control this product. For example the supports \(\{q_1,q_2\}\) and
\(\{q_1,q_3\}\) have rank two and their joint probability differs from the
product of the individual probabilities by a factor \(1/p_{q_1}\).
The reusable step is to estimate void probabilities from a local mean and
control of joint events.

For the present critical model take the same \(Y\) on every member of a
countable family of geometric grids with ratios tending to \(1\). Set
\(\mu_n=\mathbb E[r_k(n)\mid\mathcal F_Y]\). Let \(\Delta_n\) sum the
conditional intersection probabilities over distinct unordered pairs of
possible representations of \(n\) with intersecting unexposed prime
supports. Events with zero conditional probability can be dropped.
The following two estimates, if proved almost surely on these grids,
suffice for the gap density to be \(\mathbb E_U e^{-B S(U)}\):
\[
 \frac1X\sum_{n\le X}\Delta_n=o(1),\qquad
 \frac1X\sum_{n\le X}\delta_{\mu_n}
       \ \Longrightarrow\ \operatorname{Law}_U(B S(U)).                 \tag{6}
\]
For the implication, remove the \(o(X)\) targets already represented by
exposed primes. On the other targets every event probability is at most
\(\epsilon=Y^{-\gamma}\) and Harris gives the bound in the proof of (4).
Janson–Warnke's equation (2) at the zero event gives
\(\mathbb P(n\notin kC\mid\mathcal F_Y)
 \le \exp(-\mu_n^2/(\mu_n+2\Delta_n))
 \le \exp(-\mu_n+2\Delta_n)\) when \(\mu_n>0\).
If \(\mu_n=0\) the gap event has probability one.
Therefore the absolute discrepancy from \(e^{-\mu_n}\) is
\(O(\epsilon+\min(1,\Delta_n))\), uniformly in \(\mu_n\ge0\).
For the lower bound use
\(e^{-u}-e^{-u/(1-\epsilon)}\le \epsilon u e^{-u}/(1-\epsilon)\).
For the upper bound use \(e^{2\Delta}-1\ll\Delta\) for \(0\le\Delta\le1\)
and the trivial bound \(1\) otherwise. Averaging, applying (6), using
concentration (10) for the total gap count and interpolating between
grid points gives the implication. The arithmetic estimates in (6) remain
open. The fixed-modulus mean in (2) constrains the putative intensity law
but allows fluctuations within each class.
The Janson and Harris bounds have the product-space hypotheses checked
above. [Janson–Warnke, equations (2)–(3)](sources/janson-warnke-1406.1248v1.pdf)

The new unconditional claims use the retained counting, residue and
forced-prime moment results. Ma–Tang's finite-field threshold is
provenance for the comparison of dependence mechanisms. The statements
here have written proofs. A direct cyclic-convolution check with exact
rational arithmetic matched the Fourier formula on 1008 residues, covering
moduli 2, 4, 16, 3, 27, 25, and 7, orders 2 through 4, and values
\(t=0,1/3,2/5,3/5\). This tests the finite local formula.
