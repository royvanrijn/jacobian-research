# Central Taylor gains in the signed-shift bootstrap

The signed-shift variance bootstrap can retain both a central
contribution to the quartic loss and a separate central contribution
to the positive Taylor slack. Both gains can be evaluated using the
same projected reference and tilted residual estimate.
The directed cover proves the balanced three-cell comparison through
$339/500$.
[Signed-shift bootstrap](THREE_CELL_SIGNED_SHIFT_BOOTSTRAP.md)

## Marked-pair estimate for the positive slack

For $p,q$ independent under $\mu_1$, put $Z=p-p_*$ and
$Z'=q-p_*$, $\ell(p)=\min(1/4,\min_i p_i)$,
$e_0=\mathbb E_{\mu_1}\ell(p)$ and $V=\operatorname{Var}_{\mu_1}(p)$.
Independently of the radial draws, sample two uniforms on $[0,1]$.
Mark each draw if its uniform is at most $4\ell$ of its own profile.
Conditioned on profiles $p,q$, the union of the mark events has
probability $4[\ell(p)+\ell(q)-4\ell(p)\ell(q)]$. For
$x,y\in[0,1/4]$ the expression $x+y-4xy$ is at most both
$x+y$ and $1/4$, since $1/4-(x+y-4xy)=(1-4x)(1-4y)/4$.
Consequently
\[
 m(p+q)=\min(1/4,\min_i(p_i+q_i))
 \ge \tfrac14\Pr(\text{at least one mark}\mid p,q).         \tag{CS1}
\]

The unconditional mark probability is $g=4e_0$. Denote by
$m_G,m_B$ and $V_G,V_B$ the conditional means and variances of
$Z$ given a mark or its absence. When $0<g<1$, the contribution
of the marked union to $\mathbb E|Z+Z'|^2$ is
\[
 2g V_G+2g(1-g)V_B+4g^2|m_G|^2+2g(1-g)|m_G+m_B|^2.
\]
At fixed difference $m_G-m_B$, completing the square in $m_G$
bounds the last two terms below by
$2g^2(1-g)|m_G-m_B|^2/(2-g)$. The first two terms are at least
$2g[gV_G+(1-g)V_B]/(2-g)$. Thus the union captures at least
$2gV/(2-g)$, with the same bound directly when $g=0$ or $g=1$.
By (CS1),
\[
 \mathbb E[m(p+q)|Z+Z'|^2]\ge \frac{e_0}{1-2e_0}\,V .     \tag{CS2}
\]
Apply (UC1) with the undamped first radial feature. In the
normalization of (SV6), (CS2) contributes
$\frac{21\pi^3}{1024}\frac{e_0}{1-2e_0}$ times $V$ from the
positive Taylor slack. The factor increases on $[0,1/4]$.
The negative Taylor slack contributes the larger of (SV9) and
(UC5). Thus, with $M\ge\mathbb E_{\mu_1}|Z|^2$ and
$m\ge|\mathbb E_{\mu_1} Z|$ with $m<\sqrt{2/3}/2$, the new gain
is bounded below by
\[
 \max\left\{\frac43-2\sqrt{2/3}\,m+2m^2-3M,\,
             \frac23-M+2e_0\right\}
          +\frac{21\pi^3}{1024}\frac{e_0}{1-2e_0}.            \tag{CS3}
\]
Here and in the checker, a lower bound on $e_0$ can replace $e_0$.
The positive and negative Taylor terms are separate in the arc
bound. Neither uses the actual first angular deficit. Equation
(CS3) replaces the previous stage gain in its entirety.
[Quartic and positive estimates](THREE_CELL_UNDAMPED_CENTRAL_VARIANCE.md)
[Tilted mean bounds](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md)

## Centrality of the projected trial

Keep the signed-shift stage notation, with $z=Avu+g_0$ and
$h=\Pi_C(Avu)-Avu$. For fixed $v$, put
$a_+=-\min_i v_i$ and $a_-=-\min_i(-v_i)$. If $s=A|u(t)|$,
the exact minimum coordinate of $p_*+\Pi_C(Avu(t))$ is
\[
 \max\{0,\,1/3-a_\pm s\},                                  \tag{CS4}
\]
where the sign is that of $u$. Indeed a feasible reference is
fixed by the projection. An infeasible reference lies in the
same sum-one affine plane as the simplex, and an interior
projection in that plane would have to equal the reference.
Hence an infeasible reference projects onto the boundary, where
its minimum is zero.

For a bin $[A_-,A_+]$, let $\tau=r^2-l_2$ be the remaining
squared $H$-norm budget of $(g_0-h)/A$, and $\Lambda$ the bound
on its $\mu_1$ norm from (SV13). Write $W=d\mu_1/(J'_m\,d\nu)$.
The functional $\ell$ on simplex profiles is $\sqrt{2/3}$-Lipschitz
in the sum-one affine plane. This follows by projecting each label
coordinate functional onto the centered plane and taking the
minimum and then the cap. On every radial set $E$, (CS4) and two
applications of Cauchy–Schwarz give
\[
 \begin{split}
 e_0\ge{}&\int_E\min\{1/4,(1/3-a_\pm A_+|u|)_+\}\,d\mu_1\\
 &-\sqrt{2/3}\,A_+
 \min\left\{\Lambda\sqrt{\mu_1(E)},\,
             \sqrt{\tau\int_E W\,d\mu_1}\right\}.          \tag{CS5}
 \end{split}
\]
On an angular direction bin use the upper bound on $a_+$ or
$a_-$ from the appropriate endpoint of the sorted arc. A panel
with an uncertain sign uses $\sqrt{2/3}$ as an upper bound for
both. This is the opposite enclosure direction from the lower
clipping costs (SS2). In each panel of $x>0$, enclose $|xu(x^2)|$
above and divide by the lower endpoint of $x$ for a bound on
$|u|$. The panel touching zero is omitted. The applicable lower
integrand in the first line of (CS5) is constant $1/4$ below one
threshold in this upper trial magnitude and affine between it
and a second threshold. Sorting the upper trial magnitude for
each sign allows prefix sums for the two cases. For the affine
term its value is nonnegative on every selected panel, so the
lower tilted mass gives a lower integral bound.

The checker uses five sets $E$, choosing all panels of that sign
whose bound in (CS4) is at least $0,1/100,1/30,1/15,1/10$, where
the zero-threshold case includes only panels with nonnegative
unclipped minimum bound. For the error term it sums directed
upper bounds for the tilted mass and for the weighted tilted
mass, with $W(x^2)=4d/(Z_1 Q(x))$. Take the maximum of zero,
the bounds (CS5) from those sets, and the raw-moment bound
$(1/4-2M)_+$ from (UC3). This is the central mass input in (CS3).
The signed-shift clipping estimate and the two tilted second
moment bounds are the same ones used in the preceding cover.

## Verified correlation cover

The [checker](check_three_cell_central_shift_bootstrap.py) and
[full receipt](three_cell_central_shift_bootstrap_verification.json)
cover $[671/1000,339/500]$ in 242 adjacent rational slabs, with
largest full operator norm bound below $1.998819$. The checker
uses adaptive subdivision from width $1/16000$ and has at most
23 gain advances on accepted slabs. All other spectral, radial,
direction and amplitude parameters are as in the signed-shift
cover, including the $10^{-5}$ preliminary radial reserve and
the $1999/1000$ final norm cap. The previously rejected
signed-axis, graded and localized calculations do not enter this
cover.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_central_shift_bootstrap.py --workers 8

The adjoining lower interval and fixed-mass dimension reduction
give regular simplex optimality and equality rigidity for every
balanced three-cell Gaussian partition in dimension $d\ge2$ and
every $0<\rho\le339/500$. On the added interval $\kappa<1$ and
$\gamma\sqrt2<4$. The radial reserve and retained first angular
deficit provide a planar transfer-cost margin with
$w(t)=\sqrt{t}\,e^{-t-4\sqrt{t}}$. The same margin and (CT5) give a
positive unevaluated right neighborhood.
[Planar transfer cost](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md)
[Correlation continuation](THREE_CELL_CORRELATION_CONTACT.md)

For each fixed $0<\rho\le339/500$ and $\varepsilon>0$, there are
$D\in\mathbb N$ and $\tau>0$ such that, uniformly in $n$, every
$F:[3]^n\to\Delta_3$ with uniform mean $\mathbf1/3$ and
$\operatorname{Inf}_i^{\le D}(F_j)\le\tau$ for all $i,j$ has
stability under retain-and-resample noise of correlation $\rho$
at most the infinite-voter plurality limit plus $\varepsilon$.
Theorem 7.1 supplies the balanced Gaussian partition for transfer,
and Lemma 7.3 identifies the limiting plurality stability.
[Isaksson–Mossel](https://arxiv.org/pdf/0903.3362)
