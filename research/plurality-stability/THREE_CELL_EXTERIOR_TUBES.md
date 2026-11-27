# Exact exterior comparison in three conical tubes

Let $T$ be the central planar $120^\circ$ simplex. For each outward
interface ray, use orthonormal coordinates $(s,t)$, with $s$ along
the ray and $t$ its signed normal. For a sufficiently small fixed
$\eta>0$ and a sufficiently large fixed $R>0$, let the three
exterior tubes be $s\ge R,\ |t|\le\eta s$ in the respective
coordinates. There is $\rho_0<1$ such that for
$\rho_0<\rho<1$, every ordinary partition $U$ agreeing with $T$
off these tubes has
\[
                  S_\rho(T)-S_\rho(U)\ge0 ,
\]
with equality only for $U=T$ almost everywhere. All three labels
may occur in each tube. The comparison requires neither balance
nor smoothness, and the tubes have no outer radial cutoff.

Write $e=\sqrt{2(1-\rho)}$, $\sigma=\sqrt{1-\rho^2}$, and
let $\phi,\Phi,\gamma_1$ denote the scalar Gaussian density,
distribution function and measure. Throughout, $\sigma$ and $e$
are comparable. Denote Gaussian ordered disagreement by $N=1-S$.

## Binary threshold graphs

Let $H=\{t>0\}$, and consider the binary competitor whose fiber at
$s$ is $t>g(s)$, with $g=0$ on $s<R$. The measurable
thresholds may take infinite values. Set
\[
 p(s)=\Phi(g(s))-\tfrac12,\quad
 V=\int p^2\,d\gamma_1,\quad
 \mathcal B=\mathbb E_\rho
        \min\left(\frac{(p(S)-p(S'))^2}{e},\,|p(S)-p(S')|\right),
\]
where $S,S'$ are standard tangential Gaussian observations of
correlation $\rho$. The correlated normal pair is independent
of the tangential pair.

For real thresholds $u,v$, the [binary graph identity (P12)](THREE_CELL_NEAR_ONE.md#a-uniform-comparison-for-graph-deformations-away-from-the-vertex) expresses the
normal pair disagreement difference as half the sum of the diagonal
changes $D_\rho(u),D_\rho(v)$ plus
$\Pr_\rho(X,Y\in I)$, with $I$ the interval between $u,v$.
The identity extends to infinite thresholds by continuity. The
diagonal change satisfies
\[
 D_\rho(u)\ge-Ce\,(\Phi(u)-\tfrac12)^2 .
\]
For $|u|\le1$, integrate the (P12) derivative bound
$|D_\rho'(u)|\le Ce|u|$ and use
$|\Phi(u)-1/2|\ge c|u|$. For $|u|>1$, the reference normal
disagreement is at most $Ce$, by integrating the normal midpoint
density on the interval between the two endpoints, and
$|\Phi(u)-1/2|\ge c$.

For any interval $I$ of Gaussian measure $d$, one has
\[
 \Pr_\rho(X,Y\in I)\ge c\min(d^2/e,d).
 \tag{CT1}
\]
Keep a subinterval on one side of zero of measure at least $d/2$,
and reflect it if needed. Write it as $I_0=[a,b]$, $a\ge0$.
If $b-a\le\sigma$, the Mehler density relative to
$\gamma_1\otimes\gamma_1$ is at least $c/\sigma$ on
$I_0^2$. Indeed this density is
\[
 \sigma^{-1}\exp\left[
 -\frac{\rho}{2\sigma^2}(x-y)^2+
 \frac{\rho}{2(1+\rho)}(x^2+y^2)\right].
\]
If $b-a\ge\sigma$ and $a\sigma\ge1$, the first block
$[a,a+\sigma]$ contains a constant fraction of $\gamma_1(I_0)$.
For $x\ge a$ use
$\phi(x+\sigma)\le e^{-a\sigma}\phi(x)$ and integrate to see
that $\gamma_1([a,a+\sigma])\ge(1-e^{-a\sigma})\gamma_1([a,\infty))$.
The same near-diagonal density bound on this block gives (CT1).
If $b-a\ge\sigma$ and $a\sigma\le1$, restrict the first
endpoint to $I_0\cap[a,a+2]$, containing a fixed fraction of
$\gamma_1(I_0)$ by the same Gaussian tail ratio with shift two.
For such $x$, $|x-\rho x|\le C\sigma$, and
$I_0\cap[x-\sigma,x+\sigma]$ has length at least $\sigma$.
The conditional Gaussian assigns this intersection probability
at least $c$. This proves (CT1) in the last case. Passing
from $I_0$ back to $I$ changes only the absolute constant.
These facts and (P12) give
\[
 N(\{t>g(s)\})-N(H)\ge c_1\mathcal B-C_1eV.
 \tag{CT2}
\]

The antisymmetric test $\operatorname{sgn}(S-S')$ supplies a
drift estimate. Conditional on $S=s$, its expectation is
$2\Phi((1-\rho)s/\sigma)-1$, which is at least
$c\min(eR,1)$ for $s\ge R$. Thus, by exchanging $S,S'$,
\[
 \begin{aligned}
 c\min(eR,1)V
 &\le \mathbb E_\rho[\operatorname{sgn}(S-S')p(S)^2]\\
 &=\tfrac12\mathbb E_\rho[
       \operatorname{sgn}(S-S')(p(S)^2-p(S')^2)]
 \le C\sqrt{e\mathcal B V}+C\mathcal B .
 \end{aligned}
 \tag{CT3}
\]
For the last inequality, on $|p(S)-p(S')|\le e$ use
Cauchy--Schwarz with the quadratic part of $\mathcal B$;
on its complement use $|p(S)|+|p(S')|\le1$ and the linear
part. At each fixed $R\ge1$ and all small $e$,
(CT3) implies $\mathcal B\ge c\min(R^2,R)eV$: for $V>0$
divide by $eV$ and use $eR\le1$. The $V=0$ case is
immediate. Choose $R$ fixed large enough that (CT2) gives
\[
                   N(\{t>g(s)\})-N(H)\ge c_2 e V .
 \tag{CT4}
\]

## Binary fiber sorting and elimination of the extra label

For a binary, ordinary competitor $A$ that agrees with $H$
on $s<R$, let $A^*$ consist of the increasing half-line
fibers with the same Gaussian masses as $A_s$, and define
$p(s)=1/2-\gamma_1(A_s)$ as above. By Borell's two-set
inequality, for each tangential fiber pair the sorted fibers
increase the probability of joint membership and hence
decrease binary disagreement. This is the one-dimensional case of Borell's two-set comparison.
[Mossel–Neeman, Theorem 1.1](https://arxiv.org/pdf/1210.4126).

Define the conditional sign score
\[
 B_\rho(t)=2\Phi(\rho t/\sigma)-1,\qquad
 L_s=\int |1_{A_s}-1_H|\,|B_\rho(t)|\,d\gamma_1(t),\qquad
 L_s^*=\int |1_{A_s^*}-1_H|\,|B_\rho(t)|\,d\gamma_1(t).
\]
Every wrong binary label has sign opposite to $1_H$.
When the other fiber is the unchanged $H$, the exact
binary disagreement improvement on replacing $A_s$ by
$A_s^*$ is
\[
 2\int(1_{A_s^*}-1_{A_s})T_\rho^{(1)}1_H\,d\gamma_1
     =L_s-L_s^*,
\]
because the fiber masses agree and
$2T_\rho^{(1)}1_H=1+B_\rho$. In particular $L_s\ge L_s^*$.
Put $q_R(s)=\Phi((R-\rho s)/\sigma)=\Pr(S'<R\mid S=s)$.
Using the nonnegative Borell improvement on the other fiber
pairs, counting both orders of the pairs with one unchanged
fiber, and invoking (CT4), gives
\[
         N(A)-N(H)\ge c_2eV+
                  2\int_R^\infty q_R(s)(L_s-L_s^*)\,d\gamma_1(s).
 \tag{CT5}
\]
Also $L_s^*\le C p(s)^2/\sigma$. For $|g(s)|\le1$ this
follows from $|B_\rho(t)|\le C|t|/\sigma$ and
$|g(s)|\le C|p(s)|$. For larger $|g(s)|$ use
$L_s^*\le1$ and $|p(s)|\ge\Phi(1)-1/2$.

Now fix one tube with incident labels $A,B$ and remaining
label $M$. For the binary reference $H$, give the two
halfplanes the corresponding incident labels and treat $M$
as absent. Let $U_i^H$ use the labels of $U$ inside this
tube and those of $H$ elsewhere. The notation $i$ here
indexes the ray. Merge every introduced $M$ in $U_i^H$
into either $A$ or $B$, taking the one of the two global
merges with the larger cut decrease, and call the resulting
binary partition $A_i$. Put \(m_i(s)=\gamma_1\{t:
U_i^H(s,t)=M\}\). The average cut decrease from the two
possible merges is the full ordered interaction between
$M$ and its complement divided by two. In particular,
counting only pairs with the other tangential coordinate
$S'<R$ in its two orders, that average is at least
$\int_R^\infty q_R(s)m_i(s)\,d\gamma_1(s)$. Therefore
\[
 N(U_i^H)-N(H)\ge c_2e V_i+
 2\int_R^\infty q_R(s)(L_{i,s}-L_{i,s}^*)\,d\gamma_1(s)
 +\int_R^\infty q_R(s)m_i(s)\,d\gamma_1(s).
 \tag{CT6}
\]
Here $V_i,L_{i,s},L_{i,s}^*$ are the quantities for the
post-merge $A_i$. The merge changes only points previously
assigned $M$. Thus the incident-to-incident transfers of
$U_i^H$ are still transfers in $A_i$ and their binary
linear regret in a fiber is bounded by $L_{i,s}$.

## Passing from the binary references to the tripod

The two incident halfspaces for ray $i$ agree with the
corresponding cells of $T$ after removing the third wedge.
This wedge is in the tangential region $s'<0$, with
centered interval normal fibers. Accordingly the difference
between the halfspace and tripod pair-score differences is
\[
 T_\rho\big(\operatorname{sgn}(t)1_{T_M}\big)(s,t).
\]
At $t>0$, the conditional integral of $\operatorname{sgn}$
over any centered interval is nonnegative and at most
$B_\rho(t)$. Indeed pair each positive normal coordinate
with its negative and compare the two Gaussian densities
with conditional mean $\rho t$; enlarging the centered
interval adds only nonnegative contributions. There is the
odd reflection formula for $t<0$. Thus the incident-transfer
error in the linear regret is at most
$\int_R^\infty q_0(s)L_{i,s}\,d\gamma_1(s)$, with
$q_0(s)=\Phi(-\rho s/\sigma)$.
For a transfer into $M$, each single score difference
between $T$ and $H$ is bounded by $q_0(s)$, since the
changed reference cell indicators are supported in the
third wedge. The additional regret error is at most
$2\int_R^\infty q_0(s)m_i(s)\,d\gamma_1(s)$.
Uniformly on $s\ge R$,
\[
                q_0(s)\le C e^{-cR^2/\sigma^2}q_R(s).
 \tag{CT7}
\]
If $\rho s<R$, $q_R\ge1/2$ and $\rho s\ge\rho R$.
If $\rho s\ge R$, the two Gaussian tails have thresholds
$(\rho s-R)/\sigma$ and $\rho s/\sigma$, and shifting
the Gaussian tail by $R/\sigma$ bounds the ratio by
$e^{-R^2/(2\sigma^2)}$.

Let $h_i=1_{U_i^H}-1_H$ be the vector label defect with
the absent third label included. This also equals the actual
change of $U$ relative to $T$ on that tube. The individual
interaction $Q(h_i)$ is thus identical in the binary
reference and tripod comparison. To estimate mixed terms,
write the points in two distinct tubes as
$s\tau_i+t n_i$ and $s'\tau_k+t'n_k$. Separation of
the outward rays at $120^\circ$ and $|t|\le\eta s$,
$|t'|\le\eta s'$ give
\[
 |x-y|\ \ge\
 \sqrt{s^2+s'^2+ss'}-\eta(s+s')
 \ \ge\ (1-\sqrt2\eta)\sqrt{s^2+s'^2}.
\]
Choose $\eta$ small and then take $\rho$ close enough to
one so that $\rho(1-\sqrt2\eta)^2/2>a$ for some fixed
$a>1/4$. By the [conjugated kernel formula (P16)](THREE_CELL_NEAR_ONE.md#relative-control-in-the-extreme-gaussian-tail), the kernel is bounded by
$C\sigma^{-2}e^{-a(s^2+s'^2)/\sigma^2}$ on these pairs.

For each tube the normal slice obeys
\[
  \left(\int |h_i(s,t)|\sqrt{\phi(t)}\,dt\right)^2
       \le C(s+\sigma)\big(L_{i,s}+m_i(s)\big).
 \tag{CT8}
\]
For the incident transfers with $|t|\le\sigma$, this follows
from the half-line first-moment bound on each side of zero
and $|B_\rho(t)|\ge c|t|/\sigma$. For the incident transfers
with $|t|\ge\sigma$, use Cauchy--Schwarz over the normal
interval of length $2\eta s$ and the saturated score margin.
For transfers into $M$, use Cauchy--Schwarz over the same
interval with their mass $m_i(s)$. The label vector norm
on a changed point is $\sqrt2$, which alters only the
constant.

Combining (P16) and (CT8) bounds the mixed interaction by
$C\sigma^{-2}$ times the product of the two factors
\[
 \int_R^\infty e^{-a s^2/\sigma^2}
       \sqrt{(s+\sigma)\phi(s)(L_{i,s}+m_i(s))}\,ds .
\]
By Cauchy--Schwarz each factor is at most
\[
 \left(\int_R^\infty q_R(s)(L_{i,s}+m_i(s))\,d\gamma_1(s)\right)^{1/2}
 \left(\int_R^\infty
     \frac{(s+\sigma)e^{-2a s^2/\sigma^2}}{q_R(s)}\,ds\right)^{1/2}.
\]
The last integral is $O(e^{-cR^2/\sigma^2})$ since
\(q_R(s)\ge\Phi(-s/\sigma)\ge c(1+s/\sigma)^{-1}
e^{-s^2/(2\sigma^2)}\) and $2a>1/2$. Polynomial constants,
including $\sigma^{-2}$, can be absorbed by decreasing $c$
and taking $\rho$ close to one. Thus the mixed interaction
in absolute value is bounded by
\[
 C e^{-cR^2/\sigma^2}\sum_i
         \int_R^\infty q_R(s)(L_{i,s}+m_i(s))\,d\gamma_1(s).
 \tag{CT9}
\]

Finally use the exact deficit expansion $2L-Q$ to sum the three
single-tube comparisons. The self interactions are those of
(CT6); the linear regret errors and all mixed interactions
are bounded by (CT7) and (CT9). Write
$L_{i,s}=(L_{i,s}-L_{i,s}^*)+L_{i,s}^*$.
The terms with $L_{i,s}-L_{i,s}^*$ and $m_i(s)$ are
absorbed by their positive summands in (CT6). Those with
$L_{i,s}^*$ are $O(\sigma^{-1}e^{-cR^2/\sigma^2})\sum_i V_i$
and are absorbed by $c_2e\sum_i V_i$.
This leaves a positive multiple of
\[
 e\sum_i V_i+\sum_i\int_R^\infty
               q_R(s)(L_{i,s}-L_{i,s}^*+m_i(s))\,d\gamma_1(s).
\]
If the deficit is zero, then in each tube $p_i=0$,
$m_i=0$ and $L_{i,s}=L_{i,s}^*=0$ almost everywhere,
since $q_R$ is strictly positive. The binary sign score
vanishes only on the normal origin, so $L_{i,s}=0$ implies
its fiber equals the original half-line almost everywhere.
The vanished $m_i$ shows there are no extra-label changes
to undo, proving $U=T$.

## Exterior changes separated from the reference interfaces

There is also a strict comparison for a separate class of exterior
competitors. Fix $\delta>0$ and sufficiently large $R$. Suppose
every changed point lies outside $B_R$, and inside its original
simplex cell at distance at least $\delta |x|$ from its two supporting
lines. Then $U$ strictly loses to $T$ for $\rho$ sufficiently
close to one, unless the labels agree almost everywhere.

Write $m$ for the Gaussian area of the change. At such a point
the conditional Gaussian centered at $\rho x$ almost surely stays
inside the same reference cell up to probability
$C e^{-c_\delta R^2/\sigma^2}$. The simplex regret therefore
obeys $L\ge(1-C e^{-c_\delta R^2/\sigma^2})m$. On the entire
exterior region the conjugated Mehler kernel (P16) bounds
the interaction by

\[
 Q\le 2\lambda_{\rho,R}m,\qquad
 \lambda_{\rho,R}=\rho^{-1}
                     \exp[-(1-\rho)R^2/(2(1+\rho))].
\]

Indeed the second factor of (P16) gives the stated discount
and the Gaussian convolution left over has $L^2$ norm at
most $\rho^{-1}$ in the plane. The sum of squared indicator
changes is two on each changed point. For $R>2$ fixed,
$1-\lambda_{\rho,R}\ge c_R(1-\rho)$ near one, by expanding
its logarithm. This dominates the exponentially small
regret error. The deficit $2L-Q$ is consequently positive
if $m>0$.

## Arbitrary changes beyond a growing exterior radius

There are absolute $K<\infty$ and $\rho_0<1$ such that, for
$\rho_0<\rho<1$ and $R\ge K e^{-1/3}$, every ordinary planar
partition agreeing with the central simplex inside $B_R$ has
strictly smaller stability unless it agrees with the simplex
almost everywhere. No volume or angular support condition is
imposed. If $L$ denotes its simplex regret on the changed set,
the deficit is at least $c\min(e^2R^2,1)L$.

First strengthen the comparison for a single tube supported in
$s\ge r\ge R_0$. For its halfspace reference $H$, write $L_H$
for the linear regret, and use the single global merge and its
$V$ from (CT6). At a fixed $s$, merging the third label decreases
the diagonal normal cut, and one-dimensional fiber sorting
decreases it further. The normal diagonal bound preceding (CT1)
then gives a diagonal deficit of at least $-Ce p(s)^2$. By the
deficit identity this says
\[
  \langle h_s,T_\rho^{(1)}h_s\rangle\le 2L_s^H+Ce p(s)^2,
\]
with scalar products summed across labels. On the tangential
line the Gaussian operator compressed to $s\ge r$ has norm at
most
\[
 \lambda_1(r)=\rho^{-1/2}\exp[-(1-\rho)r^2/(2(1+\rho))],
\]
by the one-dimensional version of (P16). Factor the full
operator into tangential and normal factors and apply that norm
bound to $(T_\rho^{(1)})^{1/2}h_s$ in the normal variable. The
result is $Q(h)\le\lambda_1(r)(2L_H+CeV)$. For $R_0$ fixed
sufficiently large, $1-\lambda_1(r)\ge c\min(e^2r^2,1)$ for
all $r\ge R_0$ and sufficiently small $e$, by the formula for
the logarithm of $\lambda_1$. A fixed small positive weight on
this bound, with the rest of the weight on (CT6), gives
\[
 N(U_i^H)-N(H)\ge c\left(\min(e^2r^2,1)L_H+eV+
       \int_r^\infty q_r(s)(L_s-L_s^*+m_i(s))\,d\gamma_1(s)\right).
 \tag{CT10}
\]
The constants are uniform in $r$, since the drift estimate (CT3)
may use the fixed support threshold $R_0$. Changing the
reference from the halfspace to the tripod costs at most
$C\int q_0(s)(L_s+m_i(s))\,d\gamma_1(s)$ by the regret
comparison preceding (CT7). This is absorbed in the last two
terms of (CT10), since $L_s^*\le Cp(s)^2/\sigma$ and the
ratio in (CT7) is exponentially small uniformly for $r\ge R_0$.
Also $L_H$ and the simplex regret are comparable in this tube.
For incident transfers their difference is at most $q_0(s)$ times
the halfspace regret. For a transfer to the third label the
halfspace regret is at least $1/2$ and the difference is at most
$2q_0(s)$. Thus the tripod deficit of a single tube change is
at least $c\min(e^2r^2,1)$ times its simplex regret. Changes
wholly in the deep exterior of the cells likewise have deficit
at least this factor times their Gaussian mass, using the
two-dimensional compression bound above with radius $r$ and
the saturated score margin there.

To combine these bounds, choose fixed $0<\theta_1<\theta_2$
sufficiently small for the tube comparison, and choose $\theta$
uniformly in $[\theta_1,\theta_2]$. Split the changed set $D$
into four disjoint groups, one for each ray's angular tube of
half-width $\theta$, and the remaining cell interiors. For
each group compare to the partition making only that group's
label changes. In a tube the tangential support starts at
$r=R\cos\theta_2$. Put $d_R=\min(e^2R^2,1)$ and let $S$ be
the enlarged angular seam, where the smallest angle to the
rays lies in $[\theta_1/2,2\theta_2]$. The seam is at distance
at least $c|x|$ from the incident supporting lines, so simplex
regret per changed point is at least a fixed constant there.
The sum of the four localized deficits is therefore at least
$c d_R(L+\gamma_2(D\cap S))$ uniformly in $\theta$.

The full deficit differs from the sum by the negative of the
mixed interactions. Consider first pairs with
$|x-y|\le\delta\min(|x|,|y|)$ for a sufficiently small fixed
$\delta>0$. If such a pair can be assigned to different groups,
both points belong to $S$, and the probability of separation is
at most $C|x-y|/R$. They cannot lie in two distinct tubes, and
separation between a tube and the bulk occurs only when $\theta$
lies between their angles to that ray. Those angles differ by
at most $C|x-y|/\min(|x|,|y|)$. Multiply the conjugated kernel
(P16) by $|x-y|$ and drop its confinement factor. The resulting
convolution kernel has $L^1$ norm $O(\sigma)$. The defect vector
has squared norm two on $D$, so the expected absolute mixed
interaction of the close pairs is at most
$C(\sigma/R)\gamma_2(D\cap S)$.

All other changed pairs have $|x-y|\ge\delta R$. Extract
$e^{-cR^2/\sigma^2}$ from (P16), leaving a positive Gaussian
convolution of width comparable to $\sigma$ with bounded
$L^1$ norm. Its interaction on the absolute defect is at most
$CL$. For each ordered transfer split by the inward distance
$t$ to its pair line. On $t\ge\sigma$ use the convolution norm
and the saturated regret. On $0\le t\le\sigma$ drop the normal
Gaussian factor and bound squared section length by the first
normal moment, as in (P14). The tangential convolution then has
norm $C/\sigma$, and the regret controls the first moment divided
by $\sigma$. Sum the pieces by Cauchy–Schwarz for this positive
convolution operator.

Subtracting the expected absolute mixed interactions gives
\[
 S_\rho(T)-S_\rho(U)\ge
 \bigl(c d_R-C e^{-c'R^2/\sigma^2}\bigr)L
 +\bigl(c d_R-C\sigma/R\bigr)\gamma_2(D\cap S).
 \tag{CT11}
\]
The second coefficient is nonnegative for $R\ge K e^{-1/3}$,
using $\sigma\asymp e$ and $e^2R^3\ge K^3$ when $eR\le1$.
The exponential term in the first coefficient is at most
$c d_R/2$ there for sufficiently small $e$. The regret is
strictly positive for every change of positive Gaussian area
since the simplex score ties lie on lines. This proves the claim.

## Changes away from higher junctions for any cell count

Fix \(k\ge3\), let \(n=k-1\), and let \(b_i\) be the cell moments of
the regular simplex in its active space, with
\(\ell=|b_i-b_j|\) for \(i\ne j\). Use the ordered scores
\(s_{(1)}(x)\ge s_{(2)}(x)\ge\cdots\) of \(s_i(x)=b_i\cdot x\), and set
\(g_j=(s_{(1)}-s_{(j)})/\ell\). There are constants
\(K_k,c_k>0\), \(\rho_k<1\), such that for every
\(\rho_k<\rho<1\) and \(R\ge K_k e^{-1/3}\), every ordinary partition
\(U\) that agrees with \(T\) off \(\{g_3\ge4R\}\) satisfies

\[
 S_\rho(T)-S_\rho(U)\ \ge\ c_k\min(e^2R^2,1)L_\rho(U;T).
 \tag{CT12}
\]

The same statement holds in every \(d\ge n\), with the score gaps
computed in the active projection and the same constants. There is
no volume constraint. In particular, equality forces agreement
almost everywhere. Thus changes can be arbitrary on the part of the
exterior where the third-largest score is sufficiently separated
from the winner. This includes changes into labels other than the
two leading labels. For three cells \(g_3\) is comparable to \(|x|\)
in the active plane, recovering the order of the radius in (CT11).
For higher cell counts, the condition excludes neighborhoods where
three or more scores are close, even at large radius.

### One pair of incident cells

Work first in the active space. For a pair \(i,j\) put
\(x=z+tn_{ij}\), \(n_{ij}=(b_i-b_j)/\ell\), and
\[
 m_l=\frac{b_i+b_j-2b_l}{\sqrt3\ell}\quad(l\ne i,j),
 \qquad G_R=\{z\in n_{ij}^{\perp}: m_l\cdot z\ge R
                                  \text{ for every }l\ne i,j\}.
\]
The \(m_l\) are unit tangential vectors. The incident scores at \(z\)
are equal and exceed score \(l\) by
\(\sqrt3\,\ell(m_l\cdot z)/2\). For general \(t\), the greater
incident score increases by \(\ell|t|/2\) and the other labels'
scores are independent of \(t\). Thus \(T\) is the binary reference
\(H\) on \(G_R\times\mathbb R\), and at any \(z\) the region in the
normal fiber occupied by other reference labels is a centered
interval, possibly empty.

Consider changes supported on \(G_R\times\mathbb R\) against \(H\).
Write \(q_G(z)=\Pr(Z'\notin G_R\mid Z=z)\) for correlated tangential
Gaussians. Merge every introduced nonincident label into one of the
two incident labels, using the better of the two global merges.
Neither merge increases disagreement. Counting pairs with one
endpoint outside \(G_R\), in both orders, gives a gain in cut of at
least \(\int q_G m\,d\gamma_{n-1}\), where \(m(z)\) is the total
nonincident mass in the original changed fiber. Sort the merged
binary fibers by the one-dimensional two-set Borell comparison in
(CT5). Every pair of fibers gains, and pairs with one unchanged
fiber gain \(2\int q_G(L_z-L_z^*)\,d\gamma_{n-1}\) after summing both
orders. Here \(L_z,L_z^*\) are the unsorted and sorted binary regrets.

The disagreement calculation (CT1)–(CT2) for the sorted thresholds
uses only the normal pair. In (CT3), test with the sign of the
difference of the two coordinates in one fixed direction \(m_l\).
On the support \(G_R\) this coordinate is at least \(R\).
The conditional expectation of the sign is the same scalar drift
as in (CT3). Hence (CT4) with a fixed large support threshold
\(R_0\le R\) gives a reserve \(c e V\), where
\(V=\int p(z)^2\,d\gamma_{n-1}\) and \(p\) is the sorted mass change
in the notation of (CT2).

The diagonal bound of (CT10) also holds for each normal fiber with
any number of extra labels: merging and sorting the same fiber
can only decrease its cut, and the sorted diagonal deficit is at
least \(-Ce p(z)^2\). With \(h\) the vector transfer against \(H\)
and \(L_H\) its full binary-reference score loss (all extra reference
labels absent), this gives
\[
 Q(h)\le\lambda_1(R)(2L_H+CeV),\qquad
 \lambda_1(R)=\rho^{-1/2}
              \exp[-(1-\rho)R^2/(2(1+\rho))].
\]
Indeed \(G_R\) is contained in the tangential halfspace
\(m_l\cdot z\ge R\), so the tangential OU compression is bounded
by the scalar compression there, tensored with an OU contraction
in the remaining directions. Apply it after the square root of
the normal OU operator. For \(R\ge R_0\) and \(\rho\) near one,
\(1-\lambda_1(R)\ge c\min(e^2R^2,1)\).
Combining a fixed small weight of this bound with the sorting and
merging gain gives
\[
 S_\rho(H)-S_\rho(U^H)\ \ge\
 c_k\left(\min(e^2R^2,1)L_H+eV+
          \int q_G(L_z-L_z^*+m)\,d\gamma_{n-1}\right).
 \tag{CT13}
\]

To change the reference back to \(T\), let
\(q_0(z)=\Pr(\min_l m_l\cdot Z'<0\mid Z=z)\).
For \(z\in G_R\), (CT7) applied to each scalar projection gives
\[
 q_0(z)\le \sum_l\Phi(-\rho m_l\cdot z/\sigma)
          \le C_k e^{-cR^2/\sigma^2} q_G(z).
\]
The centered normal intervals of nonincident reference labels
make the incident-to-incident regret error at most
\(\int q_0 L_z\,d\gamma_{n-1}\), by the density comparison before
(CT7). Those transfers survive the merge. For transfers into any
extra label the error is at most \(2\int q_0m\,d\gamma_{n-1}\).
Split \(L_z=L_z-L_z^*+L_z^*\) and use \(L_z^*\le C p(z)^2/\sigma\).
The first and extra-label errors are absorbed by the \(q_G\)
reserve in (CT13). The sorted error is absorbed by \(eV\), since
\(\sigma^{-1}e^{-cR^2/\sigma^2}\ll e\) uniformly for \(R\ge R_0\).
The quadratic self-interaction is identical for the two references.
Their regrets are comparable: for incident changes the error
is bounded pointwise by \(q_0\) times binary regret, and for
nonincident changes the binary regret is at least \(1/2\).
This proves a simplex deficit of at least
\(c_k\min(e^2R^2,1)L\) for any single guarded pair.

### Combining different pairs and the interiors

Fix a partition with changed set \(D\subset\{g_3\ge4R\}\).
For \(\theta\in[1,2]\), split its changes into a bulk group
\(g_2>\theta R\) and pair groups \(g_2\le\theta R\) indexed
by the two highest scores. In each pair group both incident
scores exceed every other by at least \(2\ell R\). Its tangential
coordinates lie in \(G_R\), so (CT13) and the reference comparison
apply to the partition making only that group's changes.
Two distinct pair groups have mutual distance at least \(4R\):
take one label in each pair outside the other pair and test their
score difference at the two points.

In the bulk the winner beats every other by more than \(\ell R\),
and \(|x|\ge R\). Conditional Gaussian tails along all pair normals
give loss at least \(1-C_k e^{-cR^2/\sigma^2}\) per changed point.
The conjugated kernel (P16) in dimension \(n\) has prefactor
\((2\pi\sigma^2)^{-n/2}\). Its compression to \(|x|\ge R\)
has norm at most
\(\lambda_n(R)=\rho^{-n/2}\exp[-(1-\rho)R^2/(2(1+\rho))]\),
by the norm of the remaining Gaussian convolution.
The squared label defect is two per changed point. Thus for
\(R\ge R_0(k)\) the bulk deficit also bounds
\(c_k d_R L_{\rm bulk}\), where \(d_R=\min(e^2R^2,1)\).

Let \(D_{\rm seam}=D\cap\{R/2\le g_2\le3R\}\).
On the seam the score loss per change is bounded below by a
positive constant, by the same conditional Gaussian estimate.
The sum of all isolated deficits is therefore at least
\(c_k d_R(L+\gamma_n(D_{\rm seam}))\) for every \(\theta\).
Only mixed interactions remain. For pairs \(x,y\) at distance at
most \(\delta R\) with a sufficiently small fixed \(\delta>0\),
separation by the grouping can occur only between bulk and a
pair group. Both points then belong to the seam. The function
\(g_2\) is Lipschitz with a fixed constant, since the ordered
scores inherit the maximum Lipschitz constant of the linear
scores. With \(\theta\) uniform on \([1,2]\), separation probability
is at most \(C_k|x-y|/R\).
Multiplying the conjugated kernel by \(|x-y|\) leaves a convolution
of norm \(O_k(\sigma)\) after discarding the confinement factor.
Since each defect has squared length two on the seam, the expected
absolute close-pair error is at most
\(C_k(\sigma/R)\gamma_n(D_{\rm seam})\).

For all other pairs extract \(e^{-cR^2/\sigma^2}\) from the
conjugated kernel and bound the remaining widened Gaussian
convolution against the absolute defects. Its interaction is
at most \(C_k L\) as follows. For each ordered transfer split at
normal distance \(\sigma\) to its pair plane. On the thicker
part the convolution norm and the saturated score margin
control the squared Gaussian norm by its loss. On the thinner
part drop the normal Gaussian factor; the remaining tangential
convolution has norm \(C_k/\sigma\).
The squared section length, with the Gaussian square-root
density, is bounded by \(C\) times the first normal moment,
as in (P14). The saturated margin bounds that moment divided
by \(\sigma\) by its loss. Sum the pieces by Cauchy–Schwarz for
the positive convolution. This gives the stated uniform bound
in the active dimension.

As in (CT11), the resulting two errors are bounded by
\(C_k e^{-cR^2/\sigma^2}L+C_k(\sigma/R)\gamma_n(D_{\rm seam})\).
They are absorbed for \(R\ge K_k e^{-1/3}\) and \(\rho\) close
enough to one. Indeed \(\sigma/R=O(e/R)\) and
\(d_R=e^2R^2\) for \(eR\le1\), so \(eR^3\) is the
relevant ratio. If \(eR>1\), \(d_R=1\). This proves (CT12).
The saturated margin makes \(L>0\) whenever the changed set
has positive Gaussian measure.

For the ambient lift use the uncentered interaction \(Q\) in
the exact expansion \(2L-Q\), since the slices need not balance.
Condition on the correlated transverse variables \(w,w'\).
Positive definiteness of the active OU operator gives
\(\langle h_w,T_\rho^{\rm act}h_{w'}\rangle
\le[Q_{\rm act}(h_w)+Q_{\rm act}(h_{w'})]/2\).
Integrating bounds the full interaction by the mean of the
active self-interactions, while the full loss is exactly the
mean of active losses. Each slice has the required active support,
so (CT12) follows with the same constants.

The convolution argument also gives the dimension-\(n\) form of
(P14): for changes supported on \(|x|\ge R\), retain the full
confinement discount in (P16) to obtain
\(Q\le C_k e^{-c e^2R^2}L\). Thus arbitrary changes solely beyond
\(K_k/e\) strictly lose when \(K_k\) is sufficiently large,
even if they lie near higher junctions. The same slicing gives
this estimate for an active radial cutoff in higher ambient
dimension. These arguments leave changes near higher junctions
at smaller radii and their mixed interactions with other changes
as separate obligations. A three-label normal section would
require control across unequal fiber masses; the balanced
planar near-one theorem does not supply that comparison.

### Necessary junction loss for a tie and the scale of the estimate

For arbitrary changes, split the score loss into \(L_{\rm sep}\) on
\(\{g_3\ge4R\}\) and \(L_{\rm jun}\) on its complement. Apply (CT12)
to the separated piece alone. The other piece has
\(Q_{\rm jun}\le C_k L_{\rm jun}\) by the convolution bound above at
\(R=0\). The same bound on \(Q_{\rm sep}\) and positivity of OU
give a mixed term of size at most
\(C_k\sqrt{L_{\rm sep}L_{\rm jun}}\). Consequently
\[
 S_\rho(T)-S_\rho(U)
 \ge c_k d_R L_{\rm sep}
     -C_k\big(L_{\rm jun}+\sqrt{L_{\rm sep}L_{\rm jun}}\big).
\]
Thus any tie or improvement must have
\(L_{\rm jun}\ge c'_k d_R^2 L_{\rm sep}\) when (CT12) applies:
if \(L_{\rm jun}\le L_{\rm sep}\), divide by \(L_{\rm sep}>0\)
and square the resulting bound on the square root of the ratio;
otherwise the conclusion follows from \(d_R\le1\).
This holds also for arbitrary ambient dimension by active slicing.

Even exactly balanced changes along one deep facet can have relative
loss of order at most \(e^2(R+1)^2\). To see this, take two disjoint
smooth tangential bumps on balls of radius \(a=(R+1)^{-1}\),
centered a distance \(3a\) apart along an interior facet ray at
distance \(O_k(R+1)\) from the origin. Choose the positions so all
guards \(m_l\cdot z\) exceed \(6R\) on the supports. Gaussian
densities on the two balls are comparable to their common central
density, since displacement times the center norm is bounded.
Take opposite multiples of a fixed smooth bump to give a nonzero
mean-zero function \(p\) with
\(\|\nabla p\|_2^2\le C_k (R+1)^2\|p\|_2^2\).
Change the binary threshold to
\(t>\Phi^{-1}(1/2+\alpha p(z))\), for \(\alpha\) sufficiently small.
The integral of \(p\) being zero gives exact volume balance, and
all changes remain in the separated region.

At fixed \(\rho,R\), as \(\alpha\to0\) the binary-reference cost
satisfies \(L_H/\alpha^2\to(\rho/\sigma)\|p\|_2^2\).
This follows from \(B_\rho'(0)=2\rho\phi(0)/\sigma\) and the
threshold derivative \(p/\phi(0)\). The normal kernel at zero is
\(\sigma^{-1}\) relative to standard Gaussian measure, so
\(Q/\alpha^2\to(2/\sigma)\langle p,T_\rho^{\rm tan}p\rangle\);
the compact supports and small graph heights justify both limits.
The simplex loss is between a fixed positive multiple of \(L_H\)
and \(L_H\) by the incident reflection comparison, with identical
interaction. Thus the relative simplex deficit has limiting
upper bound \(C_k(1-\rho)(R+1)^2\), by
\(1-\rho^m\le m(1-\rho)\) in the tangential Hermite expansion.
This excludes a radius-independent positive relative margin at
\(R=o(1/e)\), even for smooth separated perturbations.
