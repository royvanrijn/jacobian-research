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
