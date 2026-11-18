# Forbidden-prime marginals and dilation laws

The following reductions concern the raw simultaneous Gaussian--Eisenstein
coordinate count. They prove neither its eventual positivity nor a marginal
asymptotic. They show that, in the stated valuation regimes, a bad-prime
marginal is equivalent to regularity of the raw count under multiplication
of the target by \(p^2\). A quantitative lower bound of natural scale is
not needed for this equivalence.

Put \(h(n)=r_G(n)r_E(n)/24\), restricted to \((n,6)=1\). For a fixed
set \(I\subset[0,1]\), let
\[
 R_I(N)=576\sum_{\substack{a+b=N\\a/N\in I}}h(a)h(b).
\]
The unrestricted count uses \(I=[0,1]\). Let \(B_{p,I}(N)\) count the
same coordinates with \(p\mid ab\), and set \(G_{p,I}=R_I-B_{p,I}\).
All targets below are positive \(N\equiv2\pmod{12}\). The restriction
\(I\) stays fixed under dilation; a window defined by a varying relative
width is not included. Fix \(p\equiv5,7\pmod{12}\), and put
\[
 q=p^{-2},\quad k_j=(2j+1)^2,\quad
 K(z)=\sum_{j\ge0}k_jz^j=\frac{1+6z+z^2}{(1-z)^3},\quad
 S_s=\sum_{j=0}^s k_jq^j.
 \tag{1}
\]
Write \(\beta_p(N)\) for the exact local bad proportion from the
[simultaneous-norm density calculation](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md).
At \(v_p(N)=2s+1\), the local integral is
\(I_{p,2s+1}=(1-p^{-1})S_s\), and \(\beta_p(N)=1-S_s^{-1}\).
As \(v_p(N)\to\infty\) through either parity,
\[
 \beta_p(N)\longrightarrow 1-K(q)^{-1}.
 \tag{2}
\]
The error in (2) is \(O_p((s+1)^2q^s)\) at valuations \(2s\) or
\(2s+1\). This follows directly by summing the polynomial local
coefficients against the geometric valuation measures.

## Exact identities on dilation orbits

**Proposition SimNormRenewal.** For \(v_p(M)=1\), define
\(N_s=p^{2s}M\), \(A_s=R_I(N_s)\), \(G_s=G_{p,I}(N_s)\), and
\(B_s=A_s-G_s\). Then, with negative indices zero,
\[
 A_s=\sum_{j=0}^s k_jG_{s-j},
 \tag{3}
\]
\[
 B_s+6B_{s-1}+B_{s-2}
       =9A_{s-1}-2A_{s-2}+A_{s-3}.
 \tag{4}
\]
In particular
\[
 B_{p,I}(p^2M)=9R_I(M)\quad(v_p(M)=1).
 \tag{5}
\]

Indeed \(h(p^{2j})=2j+1\) and \(h(p^{2j+1})=0\). In a supported
pair at odd target valuation the two component valuations are the same
even integer \(2j\). Removing their common \(p^{2j}\) gives a pair of
units, multiplies its weight by \((2j+1)^2\), and preserves \(a/N\).
This is a bijection on each disjoint stratum, proving (3). The generating
series identity \(A(z)=K(z)G(z)\) gives (4).

For later use take instead \(p\nmid M\). For \(s\ge1\),
\[
 A_s=\sum_{j=0}^{s-1}k_jG_{s-j}+H_s,\qquad
 H_s=(2s+1)(A_0+2sU_0),\quad 0\le U_0\le A_0.
 \tag{6}
\]
To see this, let \(u_p\) have the same local factors as \(h\) except
\(u_p(p^{2j})=1\), \(u_p(p^{2j+1})=0\), and put
\(U_0=576\sum_{a+b=M,a/M\in I}u_p(a)u_p(b)\).
One has \(h(p^{2s}t)=h(t)+2s u_p(t)\). Since \(p\nmid M\), at least
one component is a unit, where \(h=u_p\). Thus the sum of the two
ordered mixed correlations is \(A_0+U_0\). Expanding the scaled weight
gives (6), and \(H_s\le k_sA_0\). In particular
\[
 B_{p,I}(p^2M)=3R_I(M)+6U_0\quad(p\nmid M).
 \tag{7}
\]
The two ordered mixed correlations must both be kept when \(I\) is
not reflection-symmetric.

Write
\[
 J(z)=K(z)^{-1}=\frac{(1-z)^3}{1+6z+z^2}
                      =\sum_{j\ge0}\gamma_jz^j.
\]
Its first coefficients are \(1,-9,56,-328\), and thereafter
\(\gamma_j=-6\gamma_{j-1}-\gamma_{j-2}\). Hence
\(|\gamma_j|\le2\cdot7^j\). In odd orbits,
\(G_s=\sum_{j=0}^s\gamma_jA_{s-j}\) exactly. In even orbits (6) gives
\[
 G_s=\sum_{j=0}^s\gamma_jA_{s-j}+O(7^sA_0)\quad(s\ge1).
 \tag{8}
\]
For this calculation set the formal \(G_0=0\), so that
\(A(z)=K(z)G(z)+H(z)\). The quadratic polynomial \(H_s\) has generating
series \(P(z)/(1-z)^3\), with \(\deg P\le2\) and coefficients
\(O(A_0)\); multiplication by \(J\) proves (8).

## Equivalence assuming only positivity

**Theorem SimNormOddDilation.** Suppose \(R_I(N)>0\) for every
sufficiently large admissible target with odd \(v_p(N)\). Then the
following two uniform limits are equivalent:
\[
 \frac{B_{p,I}(N)}{R_I(N)}-\beta_p(N)\longrightarrow0
       \quad(N\to\infty,\ v_p(N)\text{ odd}),
 \tag{9}
\]
\[
 \frac{R_I(N/p^2)}{R_I(N)}
    -p^{-2}\frac{I_{p,r-2}}{I_{p,r}}\longrightarrow0
       \quad(N\to\infty,\ r=v_p(N)\ge3\text{ odd}).
 \tag{10}
\]
There is no bound on the odd valuation. The valuation-one marginal is
identically zero whenever defined. In particular, at valuation three,
(9) is equivalent to
\[
 \frac{R_I(p^2M)}{R_I(M)}\longrightarrow p^2+9
       \quad(M\to\infty,\ v_p(M)=1).
 \tag{11}
\]
The required factors at \(p=5,7\) are \(34,58\), respectively.

**Theorem SimNormDeepDilation.** Suppose there is an integer \(r_0\)
such that \(R_I(N)>0\) for every admissible \(N\) with
\(v_p(N)\ge r_0\). Across both parities the uniform deep limits
\[
 \frac{B_{p,I}(N)}{R_I(N)}-\beta_p(N)\longrightarrow0
 \quad\Longleftrightarrow\quad
 \frac{R_I(N/p^2)}{R_I(N)}\longrightarrow p^{-2}
 \tag{12}
\]
are equivalent as \(v_p(N)\to\infty\). Uniform means that for each
error tolerance one threshold on \(v_p(N)\) works for every base of
the dilation orbit. Eventual positivity on each base separately, with
unbounded starting indices, is not the stated hypothesis.

### Uniform comparison of predecessor masses

In odd orbits put \(a_s=q^sA_s/S_s\); in the deep argument use
\(a_s=q^sA_s\). Common multiplication of both summands by \(p^2\)
shows \(A_s\ge A_{s-1}\), since the local coefficients of \(h\)
increase. Zeros therefore form an initial segment. Put
\(r_s=a_{s-1}/a_s\) when \(a_s>0\), and \(r_s=1\) otherwise.
Then \(a_{s-1}=r_sa_s\) at every stage, and
\(0\le r_s\le C_0=q^{-1}K(q)\), uniformly in the base.

The dilation hypothesis is exactly \(r_s\to1\). The marginal
hypothesis also gives a useful bound. Both (3) and (6) imply
\(B_s\ge9G_{s-1}\), at all odd stages \(s\ge1\) and even stages
\(s\ge2\). Writing \(t_s=G_s/A_s\) at positive stages gives
\[
 \frac{A_{s-1}}{A_s}\le\frac{1-t_s}{9t_{s-1}}.
\]
Under (9), \(t_s=S_s^{-1}+o(1)\); in the deep case it tends to
\(K(q)^{-1}\). In the respective normalizations this implies
\[
 r_s\le\frac{K(q)-1}{9q}+o(1)
 =\frac{9-2q+q^2}{9(1-q)^3}+o(1)<\frac65
\]
for sufficiently far stages. Indeed \(q\le1/25\) and
\((25/24)^3<6/5\). Either proposed asymptotic therefore ensures
\(r_s\le C=6/5\) eventually.

Only uniformly many early steps can fail this bound. In the odd
all-target limit a failed step has target below a fixed threshold;
in the deep limit its index is below a fixed threshold. Combining
these exceptions with \(r_s\le C_0\) gives, at every positive stage,
\[
 \frac{a_{s-j}}{a_s}\le C_*C^j\quad(0\le j\le s),
 \tag{13}
\]
with \(C_*\) independent of the base. For any fixed \(j\), the last
\(j\) steps are eventually all good, giving the sharper bound \(C^j\).
This also covers bases with zero initial count. No upper or lower
estimate in terms of \(N\mathcal T(N)\) has been used.

### Dilation implies the marginal

On odd orbits inversion and
\(\sum_{j=0}^s\gamma_jq^jS_{s-j}=1\) give
\[
 S_s\frac{G_s}{A_s}-1
 =\sum_{j=0}^s\gamma_jq^jS_{s-j}
                    \left(\frac{a_{s-j}}{a_s}-1\right).
\]
For each fixed lag the ratio tends to one under the dilation hypothesis.
Equation (13) controls the tails, because
\(7qC\le42/125<1\). This proves (9) uniformly in odd valuations.
In deep orbits use (8) instead. Its boundary divided by \(A_s\) is
\(O(C_*(7qC)^s)\), so it vanishes uniformly in the base. The same
summable-tail argument gives \(G_s/A_s\to J(q)=K(q)^{-1}\), proving
the deep marginal by (2).

### The marginal implies dilation

On odd orbits write \(e_s=S_sG_s/A_s-1\), with value zero at zero
stages. The direct renewal, after division by \(a_s\), gives
\[
 S_s=\sum_{j=0}^s k_jq^j\frac{a_{s-j}}{a_s}(1+e_{s-j}).
\]
Fixed lags have \(e_{s-j}\to0\), and (13) controls the remaining
polynomial-geometric tail. Thus
\[
 a_s=\sum_{j=1}^s w_{s,j}a_{s-j}+o(a_s),\qquad
 w_{s,j}=\frac{k_jq^j}{S_s-1},\quad\sum_jw_{s,j}=1.
 \tag{14}
\]
In deep orbits substitute \(G_t/A_t=K(q)^{-1}+o(1)\) in (6).
The direct boundary is at most \(C_* k_s(qC)^s=o(1)\) relative
to \(A_s\). This yields (14) with weights
\(w_j=k_jq^j/(K(q)-1)\), truncated at \(s-1\); their missing mass
is exponentially small.

Put \(\Delta_s=1-r_s\). Subtract \(a_{s-1}\) from (14) and divide
by \(a_s\) to obtain
\[
 \Delta_s=-\sum_{j\ge2}w_{s,j}
       \sum_{i=1}^{j-1}\Delta_{s-i}\frac{a_{s-i}}{a_s}+o(1).
 \tag{15}
\]
Only present indices occur. The truncated deep weight has an additional
\(o(1)\) term. The \(\Delta_s\) are bounded by \(1+C_0\).
For each fixed lag use the sharper \(C^i\) bound in (13), and use
its full \(C_*C^i\) bound only for the uniformly summable tail.
The limsup coefficient mass is at most
\[
 \frac1{9q}\sum_{j\ge2} k_jq^j(j-1)C^{j-1}
 =\frac19\sum_{t\ge1}t(2t+3)^2(qC)^t
 <\frac{25}{108}<\frac14.
\]
For the last inequality use
\(t(2t+3)^2\le25\cdot8^{t-1}\) and \(qC\le6/125<1/20\).
If \(D\) denotes the relevant limsup of \(|\Delta_s|\), (15) gives
\(D\le D/4\). Hence \(D=0\), or \(r_s\to1\), proving both
converses. This completes the two equivalence proofs.

## Scope and dependencies

The exact local coefficients and \(\beta_p\) are documented in the
[local-density proof](SIMULTANEOUS_NORM_LOCAL_DENSITIES.md); the raw
coordinate normalization and earlier finite-stratum identities are in
[the pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md).
The new reductions use positivity, multiplicativity and convergent
power-series recurrences. No external analytic theorem is imported here.
They do not prove positivity, either asymptotic, the valuation-zero
marginal, bounded positive even valuation marginals, a varying-ratio
power window, or every-target \(C_{12}+C_{12}\) coverage.
At \(p\equiv11\pmod{12}\), the previously known exact identity
\(B_{p,I}(N)=R_I(N/p^2)\) for \(p^2\mid N\) already supplies the
corresponding reduction; it is not a new theorem of this note.

The exact renewal and the uniform equivalences received separate internal
reconstructions, including a fresh adversarial review of the version
assuming only positivity. The reviews explicitly checked zero prefixes,
uniform bases, even boundaries and the contraction. External mathematical
review, formal verification and historical priority are not established.
