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
as separate obligations. The guarded triple comparison below treats three-label normal
sections with unequal fiber masses.

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

## Changes away from four-label junctions

Let \(k\ge4\) and use the score gaps of (CT12). There are
\(K_k,c_k>0\) and \(\rho_k<1\) such that for \(\rho_k<\rho<1\)
and \(R\ge K_k e^{-1/3}\), every ordinary partition \(U\)
agreeing with \(T\) off \(\{g_4\ge8R\}\) satisfies
\[
 S_\rho(T)-S_\rho(U)\ge c_k d_R L_\rho(U;T),
 \qquad d_R=\min(e^2R^2,1).
 \tag{CT14}
\]
There is no volume constraint. For \(k=4\), \(g_4\) is comparable
to the active radius since the score differences span the active
space. Thus (CT14) covers arbitrary changes solely beyond a constant
times \(e^{-1/3}\), including along triple junction rays.

### A guarded triple of incident cells

For a triple \(I\), split the active space into the plane of its
site differences and its orthogonal complement, with coordinates
\(y,z\). Put \(\bar b=\frac13\sum_{i\in I}b_i\) and
\(m_l=(\bar b-b_l)/(\ell\sqrt{2/3})\) for \(l\notin I\).
These are tangential unit vectors. In
\(G_R=\{z:m_l\cdot z\ge R\text{ for all }l\notin I\}\), the
reference \(T\) is the product \(H\) of the balanced tripod with
the tangent space. The incident normal scores have sum zero, and
the excluded sites have zero normal scores.

For changes supported on the guard, write \(B_z\) for the normal
partition, taking all \(k\) labels, and \(L_z\) for its transfer
cost against \(H\).
Put \(m(z)=|\boldsymbol\gamma_2(B_z)-\boldsymbol\gamma_2(H)|\),
\(F_z=S_\rho(H)-S_\rho(B_z)\) on the normal plane,
\(L_H=\int L_z\,d\gamma\), and \(V=\int m^2\,d\gamma\).
Merging all extra labels into one incident label can only decrease
disagreement. The merged mass defect has norm at most \(C_k m(z)\).
By the [quadratic mass bound (R15)](THREE_CELL_NEAR_ONE_EXACT.md#a-uniform-quadratic-penalty-for-unequal-masses),
\(F_z\ge-C_ke\,m(z)^2\), or \(Q_z\le2L_z+C_ke\,m(z)^2\) in
terms of the normal self-interaction. Use the halfspace compression
norm preceding (CT10), after the normal OU square root, in the halfspace \(\{m_l\cdot z\ge R\}\) containing
the guard. For \(F=S_\rho(H)-S_\rho(B)\) this gives
\[
 F\ge c d_R L_H-C_ke V.                                               \tag{CT15}
\]

There is a sharper relation between mass and loss,
\[
 m(z)^2\le C_k(\sigma+m(z))L_z,\qquad \sigma=\sqrt{1-\rho^2}.          \tag{CT16}
\]
Split incident transfers at distance \(\sigma\) from their pair
planes as in (P14). The squared changed mass in the thin region
is at most \(C_k\sigma L_z\), by the first-normal-moment strip
bound and Cauchy–Schwarz in the section variable. Thick changed
mass costs at least a constant per point. The same is true of
transfers into an extra label, since the winning incident OU score
is at least \(1/3\). Therefore \(m^2\le C_k(\sigma L_z+L_z^2)\).
Taking the minimum with \(m^2\) on the right and using
\(\min(m^2,C_kL_z^2)\le\sqrt{C_k}\,m L_z\) proves (CT16).

Let \(Z,W\) have the correlated pair law in the tangent space,
write \(b_z=1_{B_z}\) as a vector of indicators, and put
\(D_{zw}=\|T_\rho^{1/2}(b_z-b_w)\|_2^2\), \(D=\frac12\mathbb E D_{ZW}\).
The operator here acts on the normal plane. Expansion gives
\(F=\int F_z\,d\gamma+D\). Write \(N_z=1-S_\rho(B_z)\). Then
\[
 |m(z)-m(w)|\le\sqrt2\,\gamma_2(B_z\ne B_w)
 \le (D_{zw}+2N_z+2N_w)/\sqrt2 .
\]
The last inequality uses the positive form \(I-T_\rho\) to bound
its energy on \(b_z-b_w\). Conditional on \(Z=z\), the mean of
\(\operatorname{sgn}((Z-W)\cdot m_l)\) is
\(2\Phi((1-\rho)m_l\cdot z/\sigma)-1\), at least
\(c\min(eR,1)\) on \(G_R\). The mass error vanishes off the guard.
Exchange \(Z,W\) to write its expected product with \(m(Z)^3\)
as half the expectation of this sign times \(m(Z)^3-m(W)^3\).
The last display bounds that by
\(C_k\mathbb E[(m(Z)^2+m(W)^2)(D_{ZW}+N_Z+N_W)]\).
Since \(m\) is bounded and \(N_z\le Ce+(F_z)_+\), this is at most
\(C_k(D+\int(F_z)_+d\gamma+eV)\).
Use \(\int(F_z)_-d\gamma\le C_keV\) to conclude
\[
 F\ge c_k\min(eR,1)\int m^3d\gamma-C_keV.                            \tag{CT17}
\]

Average (CT15) and (CT17). If \(eR\le1\), take a large fixed
\(M_k\) and absorb the negative mass term by the cubic term
on \(m\ge M_k/R\). On \(m<M_k/R\), (CT16) bounds that term
by \(C_ke(\sigma+M_k/R)L_H\). The \(1/R\) contribution is
absorbed by \(d_R=e^2R^2\) on increasing \(K_k\) in
\(R\ge K_ke^{-1/3}\), and the \(\sigma\) contribution is
absorbed since \(\sigma=O(e)\). If \(eR>1\), use (CT15) and
\(m(z)^2\le C_k L_z\). Thus \(F\ge c'_k d_R L_H\).

On \(G_R\) the true reference gives the same defect as \(H\).
One also has \(|L_\rho(B;T)-L_H|\le C_k e^{-c R^2/\sigma^2}L_H\).
Conditional on \(z\), an excluded reference label can win at a
correlated point only when some \(m_l\cdot Z'<0\), with probability at
most \(C_k e^{-c R^2/\sigma^2}\). Given \(Z'\), the union of the excluded
winning regions in the normal plane is invariant under the reflection
exchanging two incident labels: it is defined by the maximum of their three
scores being below a common threshold. The two kernels paired by
that reflection have a difference of fixed sign on the winner
side. The lost incident margin is therefore between zero and
the full tripod margin times the stated probability, by normal
and tangential independence. A transfer into an extra label has
\(H\)-margin at least \(1/3\), so its absolute error of at most twice
that probability is also a relative bound. The identity \(2L-Q\)
gives the guarded comparison for \(T\).

### Combining the triples and the separated region

Choose \(\theta\) uniformly from \([4,5]\) and group the defect on
\(\{g_4\ge8R\}\) into \(g_3>\theta R\) and groups indexed by the
three highest labels where \(g_3\le\theta R\). The first part is
covered by (CT12). In the other parts the three incident labels
all beat each excluded label by at least \(3\ell R\), so
\(z\in G_R\). Distinct triples are separated by \(6R\), by
testing a difference of labels in their symmetric difference.
A closer pair separated by the grouping therefore crosses the
single threshold, with probability at most \(C_k |x-x'|/R\)
by Lipschitz continuity of \(g_3\).

Multiplying the conjugated kernel by distance is bounded by
\(C_k\sigma\) times a widened Gaussian. The widened convolution
argument following (CT13) bounds its interaction on the absolute
defect by \(C_k L_\rho(B;T)\). On pairs at distance at least
a constant times \(R\) the same argument gains
\(e^{-cR^2/\sigma^2}\). Averaging the isolated comparisons leaves
mixed error at most \(C_k(\sigma/R+e^{-cR^2/\sigma^2})L_\rho(B;T)\).
If \(eR\le1\) the ratio of the first error to \(d_R\) is
\(O_k(1/(eR^3))\). If \(eR>1\) it tends uniformly to zero.
Choose \(K_k\) sufficiently large. This proves (CT14) in the
active dimension. Its ambient lift uses the same active slicing
as (CT12).

### Necessary core loss for an unrestricted tie

Take \(L_{\rm far}\) on \(\{g_4\ge16R\}\) and \(L_{\rm in}\)
on its complement. Every ordinary partition satisfies
\[
 S_\rho(T)-S_\rho(U)\ge c'_k d_R L_{\rm far}-C_k L_{\rm in}.
 \tag{CT18}
\]
To see this, choose a random threshold \(q\) uniformly on
\([8R,16R]\) and split the defect into transfers with \(g_4\ge q\)
and \(g_4<q\). The isolated exterior comparison controls at least
\(c_k d_R L_{\rm far}\) by (CT14). The isolated interior piece loses
at worst \(-C_k L_{\rm in}\), using the absolute-defect convolution
bound on \(Q\). The probability that a close pair is routed to
different pieces is at most \(C_k|x-x'|/R\), because \(g_4\) is
Lipschitz. The same widened-kernel argument as in (CT14) bounds the
mean mixed error by \(C_k(\sigma/R+e^{-cR^2/\sigma^2})
(L_{\rm in}+L_{\rm far})\). Increasing \(K_k\) absorbs the
\(L_{\rm far}\) contribution into half of the exterior margin.
The remaining error is bounded by \(C'_k L_{\rm in}\).
This proves (CT18), first in the active space and then by slicing.
Thus any tie or improvement must have
\(L_{\rm in}\ge c''_k d_R L_{\rm far}\).
The same two-threshold argument for \(g_3\) gives the analogous
linear-in-\(d_R\) comparison for (CT12), with loss beyond \(8R\)
and on its complement. The sign of the inner deficit remains
the unrestricted near-one obligation.

### Exponential localization for actual maximizers

Fix \(k\ge4\) and let \(A\) be a balanced actual maximizer. Use
Heilman–Tarter to work in dimension \(k-1\), embedding a smaller
active space there if necessary, and align the reference \(T\).
Let \(\ell_A(x)\) denote the pointwise OU transfer cost, \(L=\int\ell_A\,d\gamma\),
and \(L(r)=\int_{\{g_4\ge r\}}\ell_A\,d\gamma\).
There are \(K,c,C>0\) depending only on \(k\) such that, uniformly
near \(\rho=1\), with \(r_0=K e^{-1/3}\), one has
\[
 L(R)\le C L
   \exp\!\left(-c\int_{r_0}^{R}\min(e t^2,t)\,dt\right)
 \qquad (R\ge r_0).                                                 \tag{CT19}
\]
For \(k=4\), \(g_4\) is comparable to the active radius. In particular,
a fixed fraction of the transfer cost is confined inside radius
\(O_k(e^{-1/3})\), and the fraction beyond radius
\(C_{k,p}[e^{-1}\log(1/e)]^{1/3}\) is at most \(e^p\) for any fixed \(p>0\).

The near-one compactness (P2)–(P3) holds in the fixed dimension
\(k-1\). The simplex noise energy divided by \(e\) converges to
\(\sqrt{2/\pi}\) times its Gaussian perimeter. This follows by
directional slicing across its flat facets in (P2); on a fixed
direction the translated Gaussian facet weights converge to their
original weights. The difference quotient integrated in the
midpoint is bounded by \(C_k |z|\exp(C e^2|z|^2)\), by counting
facet crossings over the translated segment and bounding the
translated weights by a wider Gaussian. This bound is integrable
uniformly in small \(e\). Milman–Neeman's Gaussian multi-bubble
theorem at \(k\le (k-1)+1\) therefore gives \(L^1(\gamma)\) convergence
to the simplex orbit for these maximizers. Choose a reference
alignment giving this convergence. Exact stationarity of \(A\) or a
description of its interfaces is unnecessary for the following
fixed-volume repair.

Let \(D=\{A\ne T\}\) and \(D_r=D\cap\{g_4\ge r\}\). The function
\(g_4\) is continuous and piecewise linear, with gradient of length
one on each ordering chamber, so \(g_4\le |x|\) and it is 1-Lipschitz.
The strip argument for (P14) gives
\[
 \gamma(D)^2\le CL,\qquad \gamma(D_r)^2\le C e^{-c r^2}L(r).
\]
For the exterior thin strips the transverse projection is at
radius at least \(\sqrt{r^2-\sigma^2}\), and Cauchy–Schwarz in that
projection gives the Gaussian tail factor. The thick strips use
\(\gamma(D_r)\le\gamma\{|x|\ge r\}\) and the saturated margin.

Truncate the defect on \(\{g_4\ge r\}\) to obtain \(A^{<r}\).
The \(k-1\) compact volume bumps on independent simplex facets
have a uniformly invertible volume matrix for \(A^{<r}\) by
\(L^1\) alignment. For example take independent label differences
from the edges of a tree and choose disjoint bump supports away
from higher junctions. Their supports are inside \(\{g_4<r_0\}\)
when \(e\) is small. The volume adjustment has parameter
\(O(\gamma(D_r))\) by the implicit function theorem, whose Jacobian
control follows from bounded derivatives of the smooth flow
Jacobian. Its cost in disagreement is at most
\(C(\gamma(D)\gamma(D_r)+\gamma(D_r)^2)\le C e^{-c r^2}L\).
To verify the cost, the derivative at \(T\) is zero since incident
OU scores agree on each facet. The fixed-flow first derivative at
\(A^{<r}\) differs by \(O(\gamma(D))\) and its second derivative
is uniformly bounded. The pullback proof in (P15) applies in
dimension \(k-1\): each standardized flow difference is bounded
by a constant times the standardized endpoint displacement by
the Lipschitz flow derivatives. A first derivative is supported
where at least one endpoint lies in a fixed bump ball. Integrating
the displacement polynomial conditional on the changed endpoint
is uniformly bounded, using Gaussian tails if the other endpoint
must enter the ball from far away. This proves the first-derivative
comparison for measurable partitions and the uniform second
derivative bound.

Write \(h=1_A-1_T\), \(h_{<r}=h\,1_{\{g_4<r\}}\) and
\(h_{\ge r}=h-h_{<r}\). Compare the balanced repaired core to
the maximizer using the exact deficit \(2L-Q\), and use (CT14)
on \(h_{\ge r}\) with scale \(r/8\). With
\(d(r)=\min(e^2r^2,1)\), increasing \(K\) gives
\[
 c d(r)L(r)\le
 2|\langle h_{<r},T_\rho h_{\ge r}\rangle|+C e^{-c r^2}L
 \quad(r\ge r_0).                                                  \tag{CT20}
\]

Choose a small fixed \(\epsilon>0\) and define
\[
 \Phi(s)=\epsilon\int_{r_0}^{\max(s,r_0)}\min(e t^2,t)\,dt,
 \qquad W(x)=\exp(\Phi(g_4(x))).
\]
Integrate (CT20) against \(\Phi'(r)e^{\Phi(r)}/d(r)\,dr\) on
\(r\ge r_0\). The left side is \(c\int(W-1)\ell_A\,d\gamma\).
For a pair \(x,y\) set \(u=|x-y|/\sigma\) and
\(v=\sigma(|x|+|y|)\). The integral of the cutoff weight over
the cutoffs that separate the pair is at most
\[
 C\epsilon\,u(1+v)\sqrt{W(x)W(y)}\,\exp(C\epsilon u v).
\]
The cutoff range has length at most \(|x-y|\), the ratio
\(\Phi'/d\) is \(\epsilon\max(1/e,r)\), and the endpoint
oscillation of \(\Phi\circ g_4\) is at most
\(\epsilon (|x|+|y|)|x-y|\).
The conjugated kernel (P16) has an exponential bound
\(\exp(-c u^2-c v^2)\). For small \(\epsilon\) it absorbs the
polynomial and the \(uv\) term, leaving a widened Gaussian with
a spare factor \(\exp(-c'\sigma^2(|x|^2+|y|^2))\).
Use the absolute defects multiplied by
\(\sqrt W\,\exp(-c'\sigma^2|x|^2)\) in the conjugated space.
The interaction is at most \(C\epsilon\int W\ell_A\,d\gamma\).
On thick strips use the convolution norm and the saturated
margin. In a thin section with tangent coordinate \(z\), the
oscillation of \(\log W\) is at most
\(\epsilon\sigma(|z|+\sigma)\) on the interval of width \(\sigma\).
Its exponential is absorbed by the spare factor
\(\exp(-2c'\sigma^2|z|^2)\) in the squared-section estimate.
The first-normal-moment bound in (P14) thus applies with \(W\).

The repair term integrates to \(O(L)\) since \(\Phi(r)\le
\epsilon r^2/2\) and \(r_0\ge K e^{-1/3}\), which absorbs the
factor \(1/e\) in the cutoff weight. These bounds give
integrability of \(W\ell_A\) by \(g_4\le|x|\). Absorb the mixed
term by taking \(\epsilon\) small. Then
\(\int W\ell_A\,d\gamma\le C L\), proving (CT19).
For \(R\ge2r_0\), the integral from \(r_0\) to \(R\) in (CT19)
is at least a constant times \(\min(e R^3,R^2)\).

For the tetrahedron, this reduces an exact near-one proof to a
comparison in a growing ball. Given any fixed \(B>0\), take
\(R=C_{k,B}[e^{-1}\log(1/e)]^{1/3}\) so that \(L(R)\le e^{2B}L\)
and let \(\widehat A\) be the balanced repaired truncation.
The uncentered interaction bound \(Q\le C L\), Cauchy–Schwarz
for OU and the repair estimate give
\[
 |\Delta_\rho(A)-\Delta_\rho(\widehat A)|\le C_B e^B L .
\]
Also \(L_\rho(\widehat A;T)=L-L(R)+O(\sigma^{-1}e^{-cR^2}L)\).
For this last bound pull back the score integrals by the repair
flow. The derivative of the difference of the integrals at
\(A^{<R}\) and \(T\) is at most \(C\gamma(D)/\sigma\) since
the scores have gradient \(O(1/\sigma)\) on the bump supports.
The score cost of flowing \(T\) itself is quadratic in the
flow parameter divided by \(\sigma\): transferred points
lie in thin compact strips and the score vanishes on the
original facets. Using the \(O(\gamma(D_R))\) parameter bound
proves the assertion. Hence a lower bound
\(\Delta_\rho(U)\ge c e^p L_\rho(U;T)\) for balanced defects
supported in these balls, for some fixed \(p\) and the radius
constant corresponding to \(B>p\), would imply exact
tetrahedral optimality and equality rigidity near one.

## Volume density and outgoing transfers of actual maximizers

Fix \(k\ge3\), work in active dimension \(d=k-1\), and let
\(A^\rho\) maximize balanced Gaussian noise stability. Set
\(e=\sqrt{2(1-\rho)}\), \(a^2=(1+\rho)/2\). Align the reference
simplex \(T\) so \(A^\rho\to T\) in Gaussian measure as
\(\rho\uparrow1\), by the compactness step of
[(CT19)](#exponential-localization-for-actual-maximizers).
That step uses [Milman–Neeman, Theorems 1.1 and 1.2](https://arxiv.org/html/1805.10961v3),
with \(k=d+1\) and positive prescribed masses. The simplex
facet approximation also gives \(1-S_\rho(A^\rho)=O_k(e)\).

For every fixed bounded physical set \(K\), the following
statements hold uniformly for \(\rho\) sufficiently close to one.

1. At every adjusted-score winning point \(x\in K\) of label \(i\),
   \[
     |A_i^\rho\cap B_s(x)|\ge c_K s^d
              \qquad(R_0 e\le s\le r_0).                       \tag{WD1}
   \]
   In particular the closed winner sets converge locally in
   distance to their reference simplex cells.
2. Write \(L=L_\rho(A^\rho;T)\) for exact OU score transfer
   cost. In simplex coordinates on \(\mathbf1^\perp\), a point
   \(x\in K\cap T_i\) where label \(j\) is a winner satisfies
   \[
                 x_i-x_j\le C_K(e+L^{1/d}).                   \tag{WD2}
   \]
   Let \(B_K\) be the opposing outgoing interaction in [(BR1)](TWO_BLOCK_SIMPLEX_JUNCTION.md#a-fixed-radius-junction-bound-with-arbitrary-incoming-transfers)
   restricted to transfers in \(K\), and let \(g_3\) be the
   raw largest-to-third coordinate gap. Then
   \[
     e^d B_K\le
       C_K L_\rho\big(A^\rho;T;\ K\cap
            \{g_3\le C_K(e+L^{1/d})+CDe\}\big)
            +C_K e^{-cD^2}L
                    \qquad(D\ge1).                           \tag{WD3}
   \]
   The notation with a final set argument restricts the
   integration defining \(L_\rho\) to that set. For
   \(D\) a sufficiently large multiple of \(\sqrt{\log(1/e)}\),
   only the shrinking higher-junction tubes can supply
   an outgoing expense appreciable relative to \(L\).
### Local comparison after volume repair

At a maximizer, first variation along convex combinations of
fixed-mass fractional label profiles gives linear optimality
at its masses. The finite concave linear-assignment value as
a function of the masses admits a supporting plane there.
Hence the cells choose the largest of
\(T_\rho1_{A_i^\rho}-\beta_i\), with centered prices \(\beta\),
as in the [pointwise winner argument](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality).
For use in that argument, fractional fixed-mass partitions
may be weak-star approximated by ordinary partitions of the
same masses using nonatomic subdivision of small cells.
Ties have measure zero, since distinct adjusted score
differences are nonconstant real analytic. A constant OU
image of an indicator difference would force the difference
itself to be constant, by its Hermite expansion.

For each edge of a label spanning tree, choose two small,
mutually disjoint physical patches in its reference facet
interior. Take all patches separated. A compact normal bump
on one patch from each edge gives an invertible \(k-1\)
column mass derivative at \(T\). The columns are nonzero
multiples of tree incidence vectors. Choose \(r_0>0\) so
small that for any center in \(K\) one can choose one patch
from each edge a uniformly positive distance from the ball
\(B_{r_0}\) about that center. The volume derivative columns
for \(A^\rho\) converge to those for \(T\) by the bulk formula
\(\int_{A_i^\rho}\mathrm{div}(\phi_1 V)\,dx\)
and Gaussian measure convergence. Second mass derivatives
are uniformly bounded by smooth density pullback. Hence an
ordinary change inside \(B_{r_0}\) can be volume-corrected
by these remote flows with parameters
\(O(\gamma(A^\rho\ne U))\). The finite choice of families
makes the inverses uniform. Take \(r_0\) small enough that
every change has mass error in the inversion radius.

The derivative of the disagreement of the original \(A^\rho\)
along these small fixed-support flows is \(O(e)\). This is
the pullback estimate in the [microscopic planar proof](THREE_CELL_NEAR_ONE.md#blow-ups-of-moving-maximizers-minimize-the-euclidean-cut),
and it works in dimension \(d\). In the joint density write
\(w=(x+y)/2,\ z=(y-x)/e\). Whenever the derivative is nonzero,
one endpoint lies in a fixed compact set, so
\(|w|\le C+e|z|/2\). The parameter derivatives of the density
and the Jacobians are bounded by \(C(1+|z|^2)\) times the
density at the flow endpoints. For small fixed parameters
the flow is uniformly bi-Lipschitz with derivative close
to the identity. The difference exponent has a fixed
Gaussian reserve compared to the joint density at noise
length \(2e\). It absorbs both the polynomial and the
midpoint change \(O(1+e|z|+e^2|z|^2)\) in the exponent.
The Hermite inequality
\(1-(1-2e^2)^m\le4[1-(1-e^2/2)^m]\)
bounds original disagreement at \(2e\) by \(O(e)\).
Stationarity of \(N_\rho+2\beta\cdot m\), tested on the
invertible family, therefore gives \(|\beta|=O(e)\).
Here \(N_\rho=1-S_\rho\) and \(m\) is the vector of cell
masses.

The repair effect on \(N_\rho(U)\) differs from its effect
on \(N_\rho(A^\rho)\) only at pairs with one changed endpoint
and one endpoint in the remote flow patches. Their separation
is uniformly positive, so their joint measure is
\(O(\gamma(A^\rho\ne U)e^{-c/e^2})\), by conditioning on
the changed endpoint. Optimality after the repair gives
\[
     N_\rho(A^\rho)\le N_\rho(U)+C_K e\,\gamma(A^\rho\ne U)
                                                               \tag{WD5}
\]
for every ordinary change supported in such a ball. No
interface regularity is required.

### Phase density and the outgoing interaction

In coordinates of length \(e\) centered at \(x\in K\), the
joint density divided by \(e^d\) is
\(W(u,v)=\phi_a(x+e(u+v)/2)G(v-u)\), where \(G=\phi_1\)
is the Euclidean unit Gaussian kernel. It satisfies
\(W\le C G\) everywhere and \(W\ge c_K>0\) on
\(|u|\le r_0/e,\ |v-u|\le1\). Let \(E=A_i^{\rho,e}\cap B_r\),
\(V(r)=|E|\), and write
\(C(r)=\int_{E\times(A_i^{\rho,e}\setminus B_r)} W\).
Recolor \(E\) to any other label \(j\). The ordered cut
difference divided by \(2e^d\) is
\(C(r)-\int_{E\times A_j^{\rho,e}}W\).
Summing the bounds from (WD5) gives
\(P_W(E):=\int_{E\times E^c}W\le k C(r)+CeV(r)\).

A unit-scale smooth mollification \(v\) of \(1_E\) has
\(\|v-1_E\|_1+\|\nabla v\|_1\le C P_W(E)\), using the
unit-displacement lower kernel bound on pairs of which one
belongs to \(E\). If the approximation error is at most
\(V(r)/4\), \(v\) is at least \(1/2\) on a subset of \(E\)
of volume at least \(V(r)/2\). The \(L^1\) Sobolev inequality
then gives \(\|\nabla v\|_1\ge c V(r)^{(d-1)/d}\).
Otherwise the approximation error is at least \(V(r)/4\).
Thus \(P_W(E)\ge c\min(V(r),V(r)^{(d-1)/d})\).
Choose \(r_0\) small enough to absorb \(CeV(r)\) since
\(eV(r)^{1/d}\le C r_0\).

Integrating from \(R/2\) to \(R\) gives
\(V(R)\ge c R\min(V(R/2),V(R/2)^{(d-1)/d})\).
Indeed each pair in \(C(r)\) contributes for an interval
of radii of length at most \(|v-u|\), and
\(\int_{R/2}^R C(r)dr\le C V(R)\) by the Gaussian first
absolute moment. The starting mass \(V(R_0)\) is bounded
below uniformly at a winner point. Its conditional Gaussian
score is at least \(1/(2k)\), since the prices are \(O(e)\)
and the unadjusted scores sum to one. In the scaled coordinates
the conditional Gaussian has bounded density and mean
\(O_K(e)\), and a fixed \(R_0\) captures all but \(1/(4k)\)
of its probability. Iterating the recurrence at dyadic radii,
the ratio \(V(R)/R^d\) has a positive lower bound because
each step replaces the preceding ratio by at least a fixed
positive multiple of its power \((d-1)/d\). Monotonicity
interpolates between dyadic radii. This proves (WD1).
Distance convergence follows by applying it to a hypothetical
winner point a fixed positive distance from \(T_i\), and the
reverse distance follows by measure convergence on balls in
the interior of \(T_i\).

For (WD2), if \(g=x_i-x_j\ge C e\), a ball of radius
\(s=\min(r_0,cg)\) about \(x\) still has gap at least \(g/2\)
from coordinate \(i\) to \(j\) and satisfies \(s\ge R_0e\).
Every actual label-\(j\) point there pays transfer cost
bounded below by a positive constant. In fact if \(b\) is
the reference winner at \(y\), the pair score against \(j\)
is at least \(c\min((y_b-y_j)/e,1)\). Use reflection between
\(b,j\). In the conditional Gaussian at \(\rho y\), restrict
to a size-\(e\) displacement patch inside \(T_b\) with
\(b,j\) gap at least \(ce\). The patch added to \(\rho y\)
stays in \(T_b\), and the reflected kernel ratio is at most
\(\exp(-c(y_b-y_j)/e)\). This proves the saturated bound.
By (WD1) the points paying the bound have volume at least
\(cs^d\), and hence \(L\ge c_K s^d\).
The case \(s=r_0\) is excluded as \(L\to0\) by Gaussian
measure convergence. This proves (WD2).

If two transfers share their old label \(i\), but their new
labels are \(j,\ell\), then at \(x,y\) separated by at most
\(De\), both \(x_i-x_j\) and \(x_i-x_\ell\) are at most
\(C_K(e+L^{1/d})+\sqrt2 De\), by (WD2) and the pair-gap
Lipschitz bound. The same holds at \(y\). Thus both are in
the tube in (WD3). Their interaction is bounded by the
integrated saturated pair-gap weight on this tube in
scaled units, via the thin-section and thick-set estimates
in [(ET4)](EUCLIDEAN_SIMPLEX_JUNCTION_TRANSFER.md).
For \(|x-y|>De\), the original scaled kernel is bounded by
\(Ce^{-cD^2}\) times a widened Gaussian convolution and
the same estimates apply on all transfers in \(K\).
Multiplication by \(e^d\), the saturated score bound and
the lower bound on the Gaussian physical marginal on \(K\)
convert these estimates to (WD3).
