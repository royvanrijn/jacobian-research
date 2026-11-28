# Undamped central radial variance

Use the [sharp arc proof (SP1–SP4)](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)
and [its weighted Taylor slack (WV1)](THREE_CELL_WEIGHTED_ANGULAR_VARIANCE.md). For $s=p+q$ set $z=s-(2/3)\mathbf1$ and
$m(s)=\min\{1/4,\min_i s_i\}$. The positive part of the aligned-arc Taylor
slack has the undamped lower bound
\[
 \mathscr P_a(s):=\frac14\left(P_\sharp(a)|z|^2
               -\sum_i z_i^2R_a(z_i)\right)
 \ge \frac{21\pi}{1024}\,a\,m(s)|z|^2,\qquad a\ge0.            \tag{UC1}
\]
In particular it can be used independently of the quartic *negative*
Taylor remainder and the full actual first-angular deficit extracted
before the scalar aligned-arc comparison.

## The positive Taylor estimate

Write $d_0=2/3-m(s)\ge5/12$ and
$P_h=(2R_a(-h)+R_a(h/2))/3$ for $h\in[0,2/3]$.
As in SP2, $P_h$ is nondecreasing from $e^{-a/2}$ to $P_\sharp$.
If $z$ has one negative component $-h$, then $h\le d_0$.
For each positive displacement $\delta\in[0,h]$, the Taylor
remainder minus $P_h\delta^2/2$ is a concave function of $\delta$.
Equalizing the two positive coordinates therefore bounds the
total Taylor ratio by $P_h\le P_{d_0}$.

To quantify this, put $\theta=\pi(2/3-hv)$ in the derivative of
$P_h$. The derivative is at least
\[
 a\pi\int_0^1 v(1-v)\sin\theta\,e^{a\cos\theta}\,dv .
\]
Indeed the combined coefficients of the negative-displacement
and positive-displacement integrals are $4/3$ and $-1/3$,
respectively, and $\sin\theta\ge\sin(\theta/2)$ and
$\cos\theta\ge-\cos(\theta/2)$ on $0\le\theta\le2\pi/3$.
For $h\in[5/12,2/3]$, retain $v\ge1/(6h)$. Here
$0\le\theta\le\pi/2$, so $e^{a\cos\theta}\ge1$ and
$\sin\theta\ge2(2/3-hv)$. The resulting integral is
\[
 A(h)=\frac29-\frac h6-\frac5{324h^2}+\frac{13}{7776h^3}
       \ge\frac{21}{256}.
\]
The rational function is concave on $[5/12,2/3]$ by differentiating
twice, and the smaller endpoint value is $21/256$ at $2/3$.
Thus $P_\sharp-P_{d_0}\ge(21\pi/256)a\,m(s)$.
This proves (UC1) in the one-negative case.

In the two-negative case write $\omega$ for their actual fraction of
the total squared norm. Then $\omega\ge1/3$. The positive-coordinate
estimate in SP2 uses $P_0=R_a(-2/3)$ and $E_0=e^{-a/2}$ and puts
the total ratio at most $(P_0+E_0)/2$. Keeping that positive-coordinate
estimate unchanged, replace $P_0$ by $R_a(-d_0)$ on the two negative
coordinates. This saves $\omega(P_0-R_a(-d_0))$. The derivative of
$R_a(-h)$ is at least $2a\pi A(h)$ on the same interval. Therefore the
total improvement over the sharp coefficient is at least
\[
 \frac23\frac{21\pi}{256}a\,m(s)
       +\left(P_\sharp-\frac{P_0+E_0}{2}\right).
\]
The second term is at least $21a/800$. To verify it, set
$F(a)=(P_0+2P_1)/3-e^{-a/2}$ and let $C=\cos(\pi v)$ under the
probability density $3v/2$ on $[0,2/3]$ and $12(1-v)$ on $[2/3,1]$.
Then $F(a)=\mathbb E e^{aC}-e^{-a/2}$. The cosine coefficients at
$l\ge1$ are
\[
 C_l=-\frac{3\sin(2l\pi/3)}{l\pi}
       +\frac{(27/2)\cos(2l\pi/3)-3/2-12(-1)^l}{(l\pi)^2}.
\]
In particular $\mathbb E C^2=(1+C_2)/2>11/25$,
$\mathbb E(C^3+C^2/2)=(3C_1+C_3+1+C_2)/4>-1/20$, and
$F'(0)=1/2+C_1>21/400$, by the same rational enclosures for
$\pi,\sqrt3$ as in SV4. For $a\in[0,2]$ use the exponential
tangent to obtain $\mathbb E C^2e^{a(C+1/2)}>11/25-a/20>1/4$.
For $a\ge2$, restricting to $v\in[0,1/4]$ gives a lower bound
$(3/128)e^{12/5}>1/4$, using the degree-six exponential sum.
Thus $F''\ge0$, $F(0)=0$, and the claim follows since
$P_\sharp-(P_0+E_0)/2=F/2$. Because $m(s)\le1/4$,
$21/800\ge(21\pi/256)m(s)/3$. The total improvement again
bounds $(21\pi/256)a\,m(s)$ before the factor $1/4$.
The balanced case $z=0$ is immediate.

## Undamped coefficient in the first radial feature

Write $\mu_1=f_1\nu/Z_1$ as in SV1, $u=p-p_*$ and
$M=\mathbb E_{\mu_1}|u|^2$. Since $Q(t,s)a=(\gamma/d)f_1(t)f_1(s)$,
(UC1) gives the undamped radial penalty
\[
 \frac{21\pi^2}{2048}\rho d\,
       \mathbb E_{\mu_1}\min(1/4,\min_i p_i)\,
       \operatorname{Var}_{\mu_1}(p).                         \tag{UC2}
\]
Indeed $\min(1/4,\min_i(p_i+q_i))\ge\min(1/4,\min_i q_i)$,
and conditioned on $q$, $\mathbb E|u(p)+u(q)|^2$ is at least
$\operatorname{Var}_{\mu_1}(p)$. The radial prefactor is
$\gamma Z_1^2/d=\rho\pi d/2$. In particular
\[
 \mathbb E_{\mu_1}\min(1/4,\min_i p_i)\ \ge\ (1/4-2M)_+,       \tag{UC3}
\]
since $\min_i p_i\ge1/3-\sqrt{2/3}|u|$, and
$(\sqrt{2/3}|u|-1/12)_+\le2|u|^2$ follows by completing a square.
In the normalization $\rho d/(2\pi)$ of SV6, the coefficient
added by (UC2)–(UC3) is
$(21\pi^3/1024)(1/4-2M)_+$. For example $M\le7/200$ adds
more than $0.114$ to the quartic-only coefficient. In a stage of
the variance bootstrap already using the original quartic estimate
to supply the entire scalar gain, this counts new positive Taylor
slack on all pairs. If adding to WV3 instead, restrict to profiles
$q$ with $\min q\ge1/4$, so all $q_i\le1/2$. Its dominant-label
products vanish and the positive Taylor gap there is unused.
This yields the undamped analogue of WV11 with $\mu_1$ in place
of $\mu_0$ and coefficient $(21\pi^2/8192)\rho d$ times the central
radial mass and variance. This central loss is available independently of the directed correlation covers.

## A general convex angular kernel

If the exponential kernel is replaced by \(k(v)=g(\cos(\pi v))\)
for a differentiable increasing convex \(g\) on \([-1,1]\), the
same positive Taylor slack is at least
\[
       \frac{7\pi}{512} g'(0)m(s)|s-2p_*|^2 .                 \tag{UC4}
\]
Use the same proof with \(g'\) in place of \(a e^{a\cos\theta}\).
On the retained part of the derivative integral \(\theta\le\pi/2\),
so \(g'(\cos\theta)\ge g'(0)\). In the two-negative case Jensen
still gives \(P_\sharp\ge(P_0+E_0)/2\); using the saved negative
coefficient gives \(2/3\) of the one-negative lower bound.
The improvement in (UC1) from differentiating \(F(a)\) is
specific to the exponential; the other steps do not use its exact
functional form.

## A stronger quartic loss for a mixed profile

The following estimate is independent of the radial kernel and
balance. For any probability law \(\mu\) on radial profiles, draw
independent \(p,q\) under \(\mu\). Put \(b=\mathbb E p\),
\(V=\mathbb E|p-b|^2\), \(P=\mathbb E(1-|p|^2)\), and
\(\ell(p)=\min\{1/4,\min_i p_i\}\) with \(e_0=\mathbb E\ell\).
For \(D=|p-q|^2\), the quartic term obeys
\[
\begin{split}
 \mathbb E[D(1-D/2)] \ \ge\ &(P+2e_0)V
       +\mathbb E[(1-|p|^2+(2-6e_0)\ell(p))|p-b|^2]\\
       &+6|\mathbb E[\ell(p)(p-b)]|^2
 \ \ge\ (P+2e_0)V .                                          \tag{UC5}
\end{split}
\]
In fact (UC5) holds for any \(\ell(p)\) with
\(0\le\ell(p)\le\min_i p_i\). The key identity is
\(1-D/2=(1-|p|^2+1-|q|^2)/2+p\cdot q\). Conditioned on \(q\),
\(\mathbb E_p D=V+|q-b|^2\). The first term in the identity
therefore contributes \(PV+\mathbb E[(1-|p|^2)|p-b|^2]\).
For the last term write \(p=\ell(p)\mathbf1+\tilde p\),
\(q=\ell(q)\mathbf1+\tilde q\) with both tilded vectors nonnegative.
Then
\[
        p\cdot q\ge \ell(p)+\ell(q)-3\ell(p)\ell(q).
\]
By independence,
\[
 \mathbb E[\ell(p)\ell(q)D]
 =2e_0\mathbb E[\ell(p)|p-b|^2]-2|\mathbb E[\ell(p)(p-b)]|^2 .
\]
Substitute these expressions into the quartic term; \(2-6e_0\ge0\)
because \(\ell\le1/3\). This proves (UC5).

In the sharp-arc comparison, the full negative Taylor remainder
has weight at least \(a/\pi^2\) on \(D(1-D/2)\) by SR2 and the
first odd coefficient of \(M-I_0\). For \(\mu=\mu_1\), (UC5) thus
gives normalized variance coefficient \(2/3-M+2e_0\) in the
units \(\rho d/(2\pi)\) of SV6. It can be maximized with the SV9
coefficient \(4/3-2r_*m+2m^2-3M\) at the same stage, where
\(m=|\mathbb E_{\mu_1}(p-p_*)|\). The exact angular deficit and the
positive Taylor reserve (UC2) are separate from both negative
Taylor estimates. Combining them gives the new normalized gain
\[
 \left[
 \max\{2/3-M+2e_0,\;4/3-2r_*m+2m^2-3M\}
       +\frac{21\pi^3}{1024}e_0
 \right]\operatorname{Var}_{\mu_1}(p) .                        \tag{UC6}
\]
The bootstrap replaces its preceding scalar gain by this
coefficient; it does not add it to the preceding scalar gain.

## Lower bounds on the central radial mass using projection

There are useful \(e_0\) lower bounds even when \(M\ge1/8\)
makes (UC3) vanish. In the trial decomposition of SV7, write
\(\bar p=p_*+\alpha\phi\), \(A=|\alpha|\),
\(\hat p=\pi_\Delta(\bar p)\), and \(h=\hat p-\bar p\).
For any \(A\in[A_-,A_+]\) and any measurable radial set \(E\),
\[
 \begin{split}
 e_0\ge{}&
 \int_E\min\left\{\frac14,\left(\frac13-r_*A_+|\phi|\right)_+\right\}
                  d\mu_1\\
 &-r_* A_+\sqrt{R^2-I_2(A_-)}
       \left(\int_E \frac{d\mu_1}{J'_m\,d\nu}\,d\mu_1\right)^{1/2},
                                                               \tag{UC7}
 \end{split}
\]
provided \(G_0/A\le R\), where \(I_2\) is the expense in SV10.
Use \(I_2(0)=0\). A negative residual budget makes the amplitude
bin infeasible and no lower bound is then needed for that bin.

Indeed for centered perturbations \(\min_i p_i\) is
\(r_*\)-Lipschitz in Euclidean norm. The projection does not
decrease the explicit lower bound
\(\min_i\hat p_i\ge(1/3-r_*A|\phi|)_+\): if the parentheses are
positive then the predicted vector is already nonnegative
with sum one, and otherwise the bound is zero. Thus
\(\ell(p)\ge\min\{1/4,(1/3-r_*A_+|\phi|)_+\}-r_*|g-h|\).
Pointwise projection gives
\(\|g-h\|_{J'_m}^2\le G_0^2-\|h\|_{J'_m}^2
\le A^2(R^2-I_2(A_-))\). Integrate over \(E\) and apply
weighted Cauchy–Schwarz to its last term. In directed computation
one can choose \(E\) separately for each amplitude bin as those
radial panels where the explicit lower bound on predicted
centrality exceeds a fixed positive threshold. Formula (UC7)
keeps the local tilt ratio inside the integral, rather than
using its global supremum near the pole. The best lower bound
among those panel selections is valid in that bin.
