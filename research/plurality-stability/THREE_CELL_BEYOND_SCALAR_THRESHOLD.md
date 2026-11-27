# Direct comparison beyond the first separated scalar obstruction

Use \(d\nu(t)=e^{-t}dt\) and the radial coefficients of
[(DG15)](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md),
writing \(K_\rho=K_+\) and \(J_\rho=J_+\). For a bounded scalar
function \(z\) on \((0,\infty)\) with \(\int z\,d\nu=0\), set
\[
 q_\rho(z)=\langle z,\mathcal K_\rho z\rangle_\nu
                       -2\int J_\rho z^2\,d\nu .
\]
Define the first separated scalar obstruction *above the certified interval* by
\[
 \rho_{\mathrm{s}}=
 \inf\{\,\rho\in[11/25,1/2]:
          q_\rho(z)>0\text{ for some bounded real }z,\ \int z\,d\nu=0\,\}.
 \tag{ST1}
\]
This definition needs no claim about the scalar criterion at
correlations below the certified range. The actual Gaussian comparison
there is already supplied by (DG10).

One has \(11/25<\rho_{\mathrm{s}}<1/2\). Moreover there are
\(\varepsilon,c>0\) such that, for every ordinary balanced planar
partition \(A\) and \(|\rho-\rho_{\mathrm{s}}|<\varepsilon\),
\[
 S_\rho(T)-S_\rho(A)\ge
      c\sum_i\int J_\rho(t)z_i(t)^2\,d\nu(t),\qquad
 z_i(t)=p_i(\sqrt{2t})-\tfrac13 .
 \tag{ST2}
\]
Equality in the simplex comparison forces a rotated labeled simplex
up to null sets. Consequently the direct comparison holds throughout
\(0<\rho<\rho_{\mathrm{s}}+\varepsilon\), and this interval is strictly
longer than the one certified in (DG10). The endpoint and the amount
of extension are implicit, with no new explicit numerical cutoff.

## Weighted compactness for the separated criterion

On any compact \(I\subset(0,1/2)\) we use the following elementary
bounds from (DG15), stated with constants depending on \(I\).
Write \(d_\rho=1-\rho^2\), \(k=\rho^2/d_\rho\), and
\(a=2\rho\sqrt{ts}/d_\rho\). Then
\[
 \begin{split}
 |K_\rho(t,s)|&\le C_I e^{-k(t+s)+a},\\
 c_I e^{-kt}\sqrt t(1+t)^{-3/2}
             &\le J_\rho(t)\le C_I\min(\sqrt t,1).
 \end{split}
 \tag{ST3}
\]
Indeed the angular integral in \(K_\rho/R_0\) has numerator
bounded by \(C e^a\) after cancellation of \(I_0(a)\), and
the remaining \(a e^{-a}\) term satisfies the same bound.
The term \(w(a)\) in \(J_\rho\) is nonnegative: the angular
weight \(1-3v\) integrates to zero on \([0,2/3]\), changes
sign at \(1/3\), and \(e^{a\cos(\pi v)}\) is decreasing.
The positive \(\theta(a)\) term alone gives the lower
bound on \(J_\rho\). After cancellation of \(I_0(a)\)
restrict its radial integration to
\(1/(2(1+t))\le s\le 1/(1+t)\). There \(a\) stays bounded
and \(\int\sqrt s\,ds\) is bounded below by a constant
times \((1+t)^{-3/2}\), while the factor \(\sqrt t\)
comes from \(a\).

For the upper bound, \(w(a)+\theta(a)\le C\min(a,1)\).
For \(a\le1\) use the zero integral of the angular
weight defining \(w\) and Taylor expansion. In general
the weighted angular integral divided by \(I_0(a)\) is
bounded, as is \(a e^{-a}/I_0(a)\).
The radial factor \(R_0(t,s)d\nu(s)\) is a conditional
probability measure, with \(\mathbb E[\sqrt s\mid t]\le
C_I(1+\sqrt t)\) from the conditional Gaussian representation.
This gives the upper bound for \(t\le1\) by the factor
\(a\), and for \(t\ge1\) by boundedness.

On \(H=L^2(\nu)\) put \(w_\rho=J_\rho^{-1/2}\) and let
\(A_\rho\) be the self-adjoint operator with kernel
\(K_\rho(t,s)/\sqrt{J_\rho(t)J_\rho(s)}\). It is convenient
here to leave the kernel uncentered. For \(y=\sqrt{J_\rho}z\),
the zero-mean constraint is \(y\perp w_\rho\) and
\[
             q_\rho(z)=\langle y,(A_\rho-2I)y\rangle_H .
\]
The bound (ST3) gives \(w_\rho\in H\) and continuity of
\(w_\rho\) in \(H\) by dominated convergence. The squared
kernel times \(d\nu(t)d\nu(s)\) is dominated on \(I\) by
\[
 C_I (ts)^{-1/2}(1+t)^{3/2}(1+s)^{3/2}
                      e^{-c_I(t+s)}\,dt\,ds .
\]
In the exponent before domination the numerator is
\(-(t+s-4\rho\sqrt{ts})/d_\rho\), which has a strictly
negative Gaussian rate since \(\sup I<1/2\).
This proves Hilbert--Schmidt compactness and continuity
of \(A_\rho\) on \(I\). Pointwise continuity of \(J_\rho\)
follows by dominated integration in (DG15).

The vectors \(\sqrt{J_\rho}z\) for bounded mean-zero \(z\)
are dense in \(w_\rho^\perp\). Given \(y\) in that
hyperplane, truncate \(z=y/\sqrt{J_\rho}\) to bounded
functions, which converge in weighted \(L^2(J_\rho\,d\nu)\).
Their ordinary means converge to zero by \(w_\rho\in H\).
Subtract these means, preserving convergence since
\(\int J_\rho\,d\nu<\infty\). Consequently semidefiniteness
on the tests in (ST1) is exactly semidefiniteness on
\(w_\rho^\perp\). The orthogonal projection onto
\(w_\rho^\perp\) varies continuously in operator norm.

At \(11/25\), (DG47) gives a strict cap below two for its
centered weighted operator. The additive centering terms
vanish on \(w_\rho^\perp\), so its cap bounds the compressed
uncentered operator used here. Continuity above keeps the
cap strict in a neighborhood. On the other hand (DG5)
gives a bounded mean-zero \(f\) with \(q_{1/2}(f)>0\).
For a fixed bounded test, its form is continuous also
at \(1/2\): use the unweighted kernel and include the
two \(\nu\) factors, whose resulting exponent satisfies
\[
 -[t+s-2\rho\sqrt{ts}]/d_\rho\le-c(t+s)
\]
on a compact neighborhood of \(1/2\). The linear term
\(\int J_\rho f^2\,d\nu\) is continuous as well, using
its defining double integral and the same domination.
Thus (ST1) is nonempty and gives
\(11/25<\rho_{\mathrm{s}}<1/2\). The form is
nonpositive at \(\rho_{\mathrm{s}}\) by continuity on
fixed bounded mean-zero tests, hence also on its
weighted hyperplane.

## A nonzero scalar null profile has a strict angular loss

Suppose bounded balanced circular mass deviations \(z_i\),
\(p_i=1/3+z_i\in[0,1]\), satisfy
\(q_{\rho_{\mathrm{s}}}(z_i)=0\) for every \(i\). If one is
nonzero, the profile has a nonzero loss in the coupled
bound (DG29).

To see this, for each scalar null function \(z_i\),
test \(z_i+\tau v\) with any bounded mean-zero \(v\)
and both signs of \(\tau\). Semidefiniteness and nullity
force the bilinear first derivative to vanish. Thus
\[
      \mathcal K_{\rho_{\mathrm{s}}}z_i
                    -2J_{\rho_{\mathrm{s}}}z_i=\lambda_i
                    \quad\text{a.e.}
 \tag{ST4}
\]
for some constant \(\lambda_i\), by orthogonality to all
bounded mean-zero test functions. The two functions
\(\mathcal K_\rho z_i\) and \(J_\rho\) are real analytic
on \(t>0\). For clarity, in the variable \(x=\sqrt t\)
and on a fixed bounded complex neighborhood of \(x>0\),
the unnormalized kernel or \(J_\rho\) integrand from
(DG15) after canceling \(I_0\) is bounded by a fixed
polynomial in \(\sqrt s\) times \(e^{C\sqrt s}\) times
the decaying radial exponential in \(s\).
Together with \(d\nu(s)\) and bounded \(z_i\) this gives
holomorphy by integration. Since \(J_\rho(t)>0\), (ST4)
gives analytic representatives of the \(p_i\) on the
positive axis. Values in \([0,1]\) extend from almost
everywhere to every positive \(t\). Each coordinate
has positive mean and therefore is nonzero except at
isolated zeros. A nonzero \(z_i\) is nonconstant by its
zero mean. There are positive-measure pairs \((t,s)\)
with both mass vectors interior and unequal.

By (DG29), the excess over simplex stability is at most
\[
 \sum_i q_\rho(z_i)-
 \iint V_\rho(t,s)\,
       \Gamma(p(\sqrt{2t}),p(\sqrt{2s}))\,d\nu(t)d\nu(s),
 \tag{ST5}
\]
where \(V_\rho>0\) for \(t,s>0\).
The function \(\Gamma\) of (DG26) is nonnegative,
continuous for two interior profiles, and positive
there when they are unequal. The strictness follows
from its rank-two Gram numerator: it can vanish only
when the positive sine-length triples are proportional.
These triples are the side lengths of triangles with
angles \(\pi p_i\) and \(\pi q_i\), respectively.
Proportional sides give identical labeled angles.
Thus the integral in (ST5) is strictly positive for
a nonzero scalar null profile. This proves the strict
Gaussian comparison at \(\rho_{\mathrm{s}}\) directly.

## Uniform direct comparison across the threshold

It remains to ensure that the strict angular loss is
relative to the size of the radial defect for nearby
correlations. Write \(N_\rho(z)=\sum_i\int J_\rho z_i^2d\nu\).
Suppose there is no \(\varepsilon,c>0\) giving (ST2).
There would then be \(\rho_j\to\rho_{\mathrm{s}}\),
ordinary balanced planar \(A^{(j)}\) and
\[
 0<N_j:=N_{\rho_j}(z^{(j)}),\qquad
 S_{\rho_j}(T)-S_{\rho_j}(A^{(j)})<N_j/j .
 \tag{ST6}
\]
The case \(N_j=0\) already satisfies the comparison
by (DG2). The \(N_j\) are uniformly bounded by (ST3).

If \(N_j\to0\) on a subsequence, the annular squared
radial deficit hypothesis of (DG48) holds for large \(j\):
\(J_{\rho_j}\) is uniformly bounded below on the fixed
annulus in that theorem. Also (ST3), with \(r=\sqrt{2t}\),
bounds \(N_j\) by a constant times
\[
 \int_0^{R_J}r\sum_i|p_i(r)-1/3|^2d\mu(r)
 +\int_{R_J}^\infty\sum_i|p_i(r)-1/3|^2d\mu(r).
\]
Apply (DG48) on a fixed compact correlation interval
around \(\rho_{\mathrm{s}}\). Its uniform positive
radial deficit contradicts (ST6).

Otherwise pass to a subsequence \(N_j\to N>0\) and
\(y_{j,i}=\sqrt{J_{\rho_j}}z_i^{(j)}\rightharpoonup y_i\)
in \(H\). Weighted mean-zero passes to the limit by
strong continuity of \(w_\rho\). Compactness of
\(A_{\rho_{\mathrm{s}}}\) and operator norm continuity
give convergence of the summed quadratic evaluations.
By (ST5)--(ST6) and nonnegativity of the integral,
\[
 2N\le\sum_i\langle y_i,A_{\rho_{\mathrm{s}}}y_i\rangle
                       \le 2\sum_i\|y_i\|^2\le 2N .
\]
The middle inequality is the weighted semidefiniteness,
and the last is weak lower semicontinuity. Hence \(y_j\)
converges strongly, all components are scalar null, and
at least one is nonzero. On a further subsequence
\(y_{j,i}\to y_i\) almost everywhere. Since
\(J_{\rho_j}(t)\to J_{\rho_{\mathrm{s}}}(t)>0\) for each
\(t>0\), the corresponding \(z_i=y_i/\sqrt{J_{\rho_{\mathrm{s}}}}\)
are bounded feasible mass deviations, with the required
means. The argument of (ST4) makes their mass vectors
interior almost everywhere and nonconstant. For almost
every pair of positive radii the limiting mass vectors
are interior, so \(\Gamma\) is continuous along this
sequence of pairs there. Also \(V_{\rho_j}\) converges
pointwise. Fatou's lemma and strictness in (ST5) give a
strictly positive liminf for its subtracted integral.
But (ST5)--(ST6) and the scalar quadratic convergence
force that integral to tend to zero. This contradiction
proves (ST2).

For correlations from \(11/25\) to \(\rho_{\mathrm{s}}\),
semidefiniteness from (ST1) proves the comparison and
the null-profile argument gives equality rigidity
at every such correlation. Below \(11/25\) the comparison is supplied by (DG10).
The equality conclusion here concerns the new interval
from \(11/25\) through the threshold neighborhood.
As in (DG10), existence and dimension reduction of a global
maximizer extend the comparison and equality to all
ambient dimensions \(d\ge2\).

## Angular loss and the discrete comparison

The same estimates control the actual first angular moment loss.
Use the notation of (DG26)–(DG29), so $u_i(t)$ is the first
circular vector moment of label $i$ at $\sqrt{2t}$ and $a_i(t)$
is its cap $\sin(\pi p_i(\sqrt{2t}))/\pi$. Put
\[
 E_\rho(A)=\iint V_\rho(t,s)\sum_i
   [a_i(t)a_i(s)-u_i(t)\cdot u_i(s)]\,d\nu(t)d\nu(s) .
\]
It is nonnegative. The cosine subtraction giving (DG29)
first bounds the Gaussian deficit from below by
$E_\rho(A)-\sum_i q_\rho(z_i)$, before replacing the exact
loss by $\Gamma$. The weighted operator norm is uniformly
bounded near $\rho_{\mathrm{s}}$, so $\sum_iq_\rho(z_i)\le C N_\rho(z)$.
Together with (ST2), this gives
\[
 S_\rho(T)-S_\rho(A)\ge c'[N_\rho(z)+E_\rho(A)]
 \tag{ST7}
\]
uniformly on a possibly smaller neighborhood of $\rho_{\mathrm{s}}$.

In that neighborhood, every balanced planar partition obeys
\[
 \inf_{O\in O(2)}\sum_i\gamma_2(A_i\mathbin\triangle OT_i)
       \le C [S_\rho(T)-S_\rho(A)]^{1/2} .
 \tag{ST8}
\]
Let $\Delta=S_\rho(T)-S_\rho(A)$ be small and positive.
Use the fixed reference interval $I_0=[1,2]$ in the $t$
coordinate and put $M_i=\int_{I_0}u_i\,d\nu$ and
$C_i=\int_{I_0}a_i\,d\nu$. The radial bound and Lipschitz
continuity of the cap give $C_i=a_*\nu(I_0)+O(\sqrt\Delta)$,
where $a_*=\sqrt3/(2\pi)$. Since $V_\rho$ is bounded below
on $I_0^2$, (ST7) gives $C_i^2-|M_i|^2\le C\Delta$ by
integrating the angular loss on $I_0^2$. The $M_i$ thus have
lengths bounded away from zero, with prescribed lengths
up to $O(\sqrt\Delta)$. They sum to zero. By the law of
cosines their unit directions $v_i=M_i/|M_i|$ agree to
$O(\sqrt\Delta)$ with an orthogonal image of the regular
three directions, after aligning one direction and the cyclic order.

Take $W(t)=\sqrt t\exp[-k_*t-C_0\sqrt t]$, with
$k_*=\max\kappa_\rho<1$ on the compact correlation neighborhood
and with a sufficiently large constant $C_0$. The explicit
kernel gives $V_\rho(t,s)\ge cW(t)$ for all $t>0$ and $s\in I_0$.
Also $J_\rho(t)\ge c W(t)$ by (ST3). The integrand defining
$E_\rho$ is nonnegative. Integrating first in $s\in I_0$
therefore bounds $E_\rho$ from below by
\[
 c\sum_i\int W(t)[C_i a_i(t)-M_i\cdot u_i(t)]d\nu(t)
 \ge c'\sum_i\int W(t)[a_i(t)-v_i\cdot u_i(t)]d\nu(t).
\]

On each circle with $p_i\in[1/6,1/2]$, cosine sorting about
$v_i$ bounds the squared measure discrepancy from the centered
arc of width $p_i$ by $C(a_i-v_i\cdot u_i)$. To see this,
the discarded and added portions relative to the optimal arc
each give the least cosine cost when placed next to its two
endpoints. Since the endpoint sine is bounded below on these
widths, the cost is bounded below by a constant times the
squared discrepancy. For discrepancies larger than a fixed
small constant, retain only a fixed small portion of the
boundary cost. The squared discrepancy on other circles is
bounded by $C|p_i-1/3|^2$. Changing widths to $1/3$ and changing
the directions to the common regular directions adds a squared
cost bounded by $C(|p_i-1/3|^2+\Delta)$ on each circle. Using
(ST7) and $J\ge c W$ bounds the resulting weighted squared
discrepancy by $C\Delta$. Finally $\int W^{-1}d\nu$ is finite
since $k_*<1$. Cauchy–Schwarz removes the weight and proves
(ST8). At $\Delta=0$ use equality rigidity, and for
$\Delta$ bounded away from zero enlarge the constant.

The proof of (ST2) by two cases also gives uniform constants
on the whole compact interval from $11/25$ to
$\rho_{\mathrm{s}}$, since every limiting parameter there has
the scalar semidefiniteness used in (ST4)–(ST6). Combining it
with the threshold neighborhood gives (ST7)–(ST8) uniformly
through any closed upper cutoff chosen inside that neighborhood.

Fix any closed cutoff $r_1$ strictly between
$\rho_{\mathrm{s}}$ and the upper end of that neighborhood.
For each $0<\rho\le r_1$ and $\epsilon>0$, Isaksson–Mossel's
invariance theorem gives $D$ and $\tau>0$ such that any
balanced rule $f:[3]^m\to\Delta_3$ under the uniform law
with all degree-$D$ coordinate influences at most $\tau$
has stability at most the ternary plurality limit plus
$\epsilon$. In their Theorem 7.1 the Gaussian realization is
an ordinary label assignment in dimension $2m$ with the
same masses and stability within $\epsilon$. The present
Gaussian comparison bounds its stability by the simplex,
which equals the plurality limit by their Lemma 7.3. The
perturbation retains each discrete coordinate independently
with probability $\rho$ and resamples uniformly otherwise.
[Isaksson–Mossel](https://arxiv.org/pdf/0903.3362)
