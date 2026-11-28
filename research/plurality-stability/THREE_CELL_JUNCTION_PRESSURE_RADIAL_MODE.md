# Junction location, prices, and the first radial feature

At a pressure-free triple junction at the origin, the degree-zero feature of the sharp scalar radial kernel has zero pairing with all three circular mass deviations. A displacement or an unequal price vector changes this constraint. We give the exact identity at a displaced junction, a bound on its displacement term, and two local constructions with nonzero degree-zero pairing. One construction has zero prices and a displaced genuine junction. The other has a centered genuine junction and unequal prices. Both obey the winning-score rule and strict second-variation stability inside a disk.

## The exact pressure identity

Fix $0<\rho<1$, and write $d=1-\rho^2$, $\kappa=\rho^2/d$, $t=|y|^2/2$, and $d\nu(t)=e^{-t}dt$. On the circle use $d\sigma=d\theta/(2\pi)$. For a balanced ordinary planar partition $A$, write $p_i(t)$ for its circular masses and $z_i(t)=p_i(t)-1/3$. Set
\[
 f_0(t)=e^{-\kappa t},\qquad
 h_i=\int f_0(t)z_i(t)\,d\nu(t),\qquad
 U_i=T_\rho 1_{A_i}.
\]
Thus $f_0$ is the degree-zero feature from the [sharp radial series](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md#consequences-for-the-radial-operator), and
\[
                  h=d(U(0)-\mathbf1/3).
\]
The prices are normalized by $\sum_i\beta_i=0$. At any point $x_0$ where all three adjusted scores $U_i-\beta_i$ agree, their common value is $1/3$. This applies in particular at a genuine triple junction under the winning-score rule and gives $U(x_0)=\mathbf1/3+\beta$.

Write $x_0=R(\cos\theta_0,\sin\theta_0)$, $a(t)=\rho R\sqrt{2t}/d$, and
\[
 m_{i,n}(t)=\int e^{in\theta}
        (1_{A_i}(\sqrt{2t},\theta)-1/3)\,d\sigma(\theta).
\]
Expanding the Gaussian likelihood ratio in the angular variable gives
\[
 \begin{split}
 d e^{\rho^2R^2/(2d)}\beta_i
  =\int e^{-\kappa t}\Big[
       I_0(a(t))z_i(t)
       +2\sum_{n\ge1} I_n(a(t))
          \operatorname{Re}\!\big(e^{-in\theta_0}m_{i,n}(t)\big)
                    \Big]\,d\nu(t).                       \tag{JP1}
 \end{split}
\]
Indeed the conditional density of $Y$ given $X=x_0$, relative to $\gamma_2$, is $d^{-1}e^{-\rho^2R^2/(2d)}e^{-\kappa t}e^{a(t)\cos(\theta-\theta_0)}$. Here $I_n$ are the Fourier coefficients of $e^{a\cos\theta}$. Their sum with the angular multiplicities is $e^a$, so boundedness of the indicators and Gaussian radial decay justify the series integral. When $R=0$, all terms of positive angular degree vanish and (JP1) becomes $h=d\beta$.

There is a global displacement bound,
\[
                |h-d\beta|\le \rho\sqrt{d/\pi}\,|x_0|.    \tag{JP2}
\]
For a unit vector $v$ and $Z\sim\gamma_2$, differentiation of the Gaussian conditional density gives
\[
 \partial_v U_i(x)=\frac{\rho}{\sqrt d}\,
           \mathbb E[(v\cdot Z)1_{A_i}(\rho x+\sqrt d\,Z)].
\]
The three scalar expectations sum to zero. The sum of their positive entries is at most $\mathbb E(v\cdot Z)_+=1/\sqrt{2\pi}$, as is the sum of the magnitudes of their negative entries. The sum of their squares is therefore at most $1/\pi$. This gives $\|\nabla U(x)\|_{\rm op}\le\rho/\sqrt{\pi d}$. Integrating on the segment from $0$ to $x_0$ proves (JP2).

## The consequence for a trial direction

Let $J$ be a positive radial weight with $\int J^{-1}\,d\nu<\infty$, and assume the $z_i$ belong to $\mathcal H=L^2(J\,d\nu)$. Take a unit trial $\phi\in\mathcal H$ with $\int\phi\,d\nu=0$. Write $z_i=\alpha_i\phi+g_i$ by orthogonal projection onto $\phi$, $A=|\alpha|$ and $G^2=\sum_i\|g_i\|_{\mathcal H}^2$. Global balance gives $g_i\perp J^{-1}$ in $\mathcal H$. Define
\[
 \begin{split}
 C_\phi&=|\langle\phi,f_0/J\rangle_{\mathcal H}|,\\
 B_\phi&=\left\|
      P_{(\operatorname{span}\{\phi,J^{-1}\})^\perp}
              (f_0/J)\right\|_{\mathcal H}.
 \end{split}
\]
The two removed functions are orthogonal by the mean-zero condition on $\phi$. Both representers are in $\mathcal H$ since $0<f_0\le1$. Then
\[
 |h|\ge A C_\phi-G B_\phi.
\]
In particular, if a separate scalar estimate shows $G\le H A$, (JP2) gives
\[
           A(C_\phi-H B_\phi)_+
                   \le d|\beta|+\rho\sqrt{d/\pi}|x_0|.     \tag{JP3}
\]
This removes such a trial coefficient at zero price and zero center when $C_\phi>H B_\phi$. With unknown center and prices it bounds a joint obstruction instead.

## Balanced sources with genuine local winner junctions

Fix $\rho$ and the regular simplex partition $T$. Choose a bounded real function $f$ of $t$, supported on a compact interval in $(0,\infty)$, such that
\[
   \int f\,d\nu=0,\qquad
   \ell=\int f_0 f\,d\nu\ne0.
\]
For example, use the difference of the indicators of two separated compact intervals of equal $\nu$-mass. The different integrals against the strictly decreasing $f_0$ give the second property. Fix a nonzero centered vector $v\in\mathbf1^\perp$. Select a bounded annulus in the $t$ coordinate, represented by an interval $E\subset(0,\infty)$ of positive $\nu$-measure and disjoint from the support of $f$. We use it for exact Gaussian mass correction.

For a radius $b>0$ small enough that the ball $B_b$ misses both exterior supports, let $F$ be any ordinary partition of $B_b$. Write
\[
 m_i(F)=\gamma_2(F_i)-\tfrac13\gamma_2(B_b).
\]
Outside $B_b$ take consecutive arcs in the same cyclic order as $T$, with one endpoint fixed at its simplex angle and prescribed masses
\[
        p_i^{\epsilon,F}(t)=\frac13+\epsilon v_i f(t)
                  -m_i(F)\frac{1_E(t)}{\nu(E)}.             \tag{JP4}
\]
Call the full input partition $\mathcal A^{\epsilon,F}$. The masses in (JP4) sum to one and are strictly positive for small $b$ and $|\epsilon|$, uniformly over $F$. Since both exterior supports lie outside $B_b$ and $\int f\,d\nu=0$, the full partition is exactly balanced. On the complement of the two supports the exterior labels coincide with $T$. Moving the arc endpoints by their cumulative mass changes bounds the exterior indicator $L^1(\gamma_2)$ difference between two such inputs by a constant times $\sum_i|m_i(F)-m_i(G)|$, with a constant independent of $b$. This sum is at most $\|F-G\|_{1,b}$, where we use the metric
\[
       \|F-G\|_{1,b}
           =\sum_i\int_{B_b}|1_{F_i}-1_{G_i}|\,d\gamma_2 .
\]
The partitions of $B_b$ modulo null sets are complete in this metric, as the almost everywhere limits of a convergent subsequence of the indicator vectors are again indicator partitions.

Define two maps of this space into itself. In the first, assign the labels inside $B_b$ by the unadjusted scores of $\mathcal A^{\epsilon,F}$. In the second, subtract the values of those scores at the origin before taking their winners. Either map preserves the exterior construction rule (JP4), which will be reapplied to the resulting interior partition. We use these maps to obtain local input-output agreement.

All changes from $T$ are supported in a fixed compact set comprising a bounded origin ball and the two exterior supports. On a fixed small ball the scores and their first derivatives are consequently Lipschitz in the $L^1(\gamma_2)$ change of the inputs. This follows by differentiation under the Mehler integral and boundedness of the kernel and its derivatives on the two compact sets. At the origin the simplex pair gradients are nonzero, and two independent pair gradients span the plane. In fact $\nabla U_{T_i}(0)$ is a common positive multiple of the regular label vectors, by the centered Gaussian first moments of the three sectors. The unadjusted simplex scores have the simplex winners, by reflecting across each pair line in the convolution. For small $b$ and $|\epsilon|$, every score field used in the two maps is close in $C^1$ to the corresponding simplex scores. In particular on $B_b$ each pair contrast has directional derivative bounded below by a fixed positive constant in the direction of its simplex pair gradient.

If two inputs $F,G$ produce different winning labels at a point in $B_b$, one of the pair contrasts for the first input is within $C\|F-G\|_{1,b}$ of zero. For each contrast, disintegrate Lebesgue measure on lines in the fixed direction with derivative bounded below. On each line the near-zero band has length at most $C'\|F-G\|_{1,b}$, and the transverse set has width at most $2b$. Gaussian density is bounded on the ball. Summing the three bands shows the map Lipschitz constant is at most $C''b$. This also shows that the tie sets have measure zero. In the origin-adjusted case the score difference between two inputs vanishes at zero, so the near-zero band has width parameter $Cb\|F-G\|_{1,b}$ and the Lipschitz constant can be replaced by $C''b^2$. Taking $b$ small gives a unique fixed interior partition $F_\epsilon$ for either map. On the disk the resulting balanced ordinary partition agrees with its own strict score winners for the corresponding prices.

At $\epsilon=0$ the fixed point is $T|_{B_b}$. The output of that input changes in metric by $O(b|\epsilon|)$ for the raw map and $O(b^2|\epsilon|)$ for the origin-adjusted map, using the same contrast-band estimates for the exterior perturbation. Dividing by one minus the contraction constant gives those bounds for the fixed points as well. Therefore for either resulting partition, denoted by $A^\epsilon$, its centered score vector is
\[
 U_{A^\epsilon}(0)-\mathbf1/3
             =\frac{\epsilon\,\ell}{d}v+O(b|\epsilon|).
                                                               \tag{JP5}
\]
Indeed the exterior profile with the simplex interior has exactly the displayed leading score. Replacing the interior by its fixed point and correcting balance outside costs $O(b|\epsilon|)$ in $L^1(\gamma_2)$ on the fixed compact support. Choose $b$ small enough that this error is smaller in norm than half the leading coefficient, and then take nonzero $|\epsilon|$ sufficiently small. Equations (JP5) and $h=d(U(0)-\mathbf1/3)$ give $h\ne0$.

For the origin-adjusted map take $\beta=U_{A^\epsilon}(0)-\mathbf1/3$. The adjusted scores tie at zero and their two independent gradients stay independent. Thus all three local winner sets accumulate at zero, making it a genuine triple junction of $A^\epsilon$ with unequal prices. For the raw map the prices are zero. The two score differences have derivative close on $B_b$ to the invertible simplex derivative at zero and value $O(|\epsilon|)$ there. They have a unique common zero near zero for small $|\epsilon|$. For example, subtract their values after multiplication by the inverse simplex derivative to obtain a contraction of a smaller closed ball into itself. The common zero is displaced because the centered score vector at zero is nonzero. Independence of the differences' gradients there gives all three local strict winners, so this point is a genuine triple junction as well.

The local stability form $Q_0$ is strictly positive for every nonzero regular normal speed $\eta$ compactly supported inside $B_b$. Orient $\Sigma_{ij}$ from $i$ to $j$ for $i<j$ and put $\nu_{\rm int}=\gamma_2\,d\mathcal H^1$ on these arcs, using the Gaussian density for $\gamma_2$ in this expression. The local term is $\sum_{i<j}\int_{\Sigma_{ij}}|\nabla(U_i-U_j)|\eta_{ij}^2\,d\nu_{\rm int}$. From this term subtract the Mehler interaction, which pairs $\eta_{ij}(x)(e_i-e_j)$ and $\eta_{kl}(y)(e_k-e_l)$ over every ordered pair of arcs. Here $e_i$ are the coordinate vectors of $\mathbb R^3$. These are the two terms in minus one half of the second variation of the priced stability $\sum_i\langle1_{A_i},T_\rho1_{A_i}\rangle-2\sum_i\beta_i\gamma_2(A_i)$. The local stationary pair levels follow from the local winner rule, and the positive normal score derivatives are bounded below uniformly. Each pair tie set in the ball is a graph in its fixed normal direction, with bounded slope by the same derivative estimates. Hence the total interface length there is $O(b)$. The Mehler kernel is bounded on the compact pair domain, so by Cauchy–Schwarz the absolute interaction is at most $Cb\|\eta\|_{L^2(\nu_{\rm int})}^2$. The local term is at least $c\|\eta\|_{L^2(\nu_{\rm int})}^2$. This proves strict stability for small $b$, in particular for zero first-order volume flux. The winning-score, stationary-level and stability conclusions here concern only the disk.

## Radial approximation by the same construction

Assume now that the radial weight $J(t)$ is bounded on compact subintervals of $(0,\infty)$ and satisfies $J(t)\le C\sqrt t$ near zero. This holds for the usual truncated and full scalar row weights at fixed positive correlation. The centered circular mass profile of either construction approximates $\epsilon v f$ in $L^2(J\,d\nu)$, with relative error $O(\sqrt b)$ for the raw map and $O(b^{3/2})$ for the origin-adjusted map.

To see this, the full input changes from $T$ by $O(|\epsilon|)$ in $L^1(\gamma_2)$ on the fixed compact supports. The scores therefore change by $O(|\epsilon|)$ in $C^1$ on the ball. On a radius $r\le b$, the simplex pair scores divided by $r$ have simple angular zeros with a uniformly positive angular derivative magnitude near the zeros. This follows by linearizing at zero and using nonzero simplex pair gradients. The raw scores change mass on that circle by at most $C\min(1,|\epsilon|/r)$. The adjusted scores change it by $C|\epsilon|$ because the score perturbation vanishes at zero. Squaring and integrating against $J(r^2/2)r e^{-r^2/2}dr$ gives squared errors at most $Cb\epsilon^2$ and $Cb^3\epsilon^2$ respectively. The exterior correction in (JP4) contributes norms $O(b|\epsilon|)$ and $O(b^2|\epsilon|)$ on its fixed annulus.

When $\int J^{-1}\,d\nu<\infty$, bounded compactly supported functions away from zero can approximate any mean-zero trial in $L^2(J\,d\nu)$ while retaining its mean-zero constraint. Truncate first, and correct the mean by a fixed bounded bump of integral one in an annulus. Continuity of the mean and of the $f_0$ projection in this weighted norm preserves any nonzero limiting $f_0$ projection. The constructions thus permit radial perturbations arbitrarily close in direction to such an adverse trial together with either kind of locally stable winner junction. An exclusion at a first tie requires an additional global property.
