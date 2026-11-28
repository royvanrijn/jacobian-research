# A Gaussian smoothing obstruction to periodic interface tails

Fix \(0<\rho<1\) and \(\sigma=\sqrt{1-\rho^2}\). In coordinates
\((s,y)\in\mathbb R\times\mathbb R^{d-1}\), let \(1_A\) be an
ordinary partition that chooses the pointwise largest of the
adjusted Gaussian scores \(T_\rho1_{A_i}-\beta_i\), with fixed
prices \(\beta_i\). In particular this applies to an actual
positive-correlation maximizer with positive prescribed masses.
The following conclusions use the score rule alone.

- If \(A\) approaches a fixed profile periodic in the outward
  coordinate, that profile is independent of the outward
  coordinate. This holds for every positive correlation,
  including reciprocal integers.
- A sequence of outward charts cannot persistently approach
  horizontally nonconstant periodic profiles from a compact
  template family with periods bounded above and below, with
  every label present in every template. In particular, a
  three-label cyclic stripe end with a slowly modulated phase
  of bounded positive frequency is impossible.

The approach is local in translated physical coordinates.
No regularity of the partition or its interfaces is assumed
in these statements.
[Pointwise assignment and fixed-mass prices for maximizers](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality)

## Periodic profiles and horizontal smoothing

Let \(P\) be a horizontally periodic label tuple of period
\(L>0\). For \(b>0\), let \(H_b^s\) denote ordinary Gaussian
convolution of variance \(b^2\) in the \(s\) coordinate.
On the product of the uniform circle of length \(L\) and
\(\gamma_{d-1}\), the operator
\[
                     K_b=H_b^s\otimes T_\rho^y
\]
is self-adjoint and positive definite. Write
\(\mathcal E_b(P)=\langle P,K_bP\rangle_L\), summing over labels,
and let \(m(P)\) be the label mean for this product measure.
These definitions are independent of the choice of period.
For periodic functions of \(s\) a mean over a period equals the
large-interval mean, which does not depend on the period used.

Write \(q(\theta,y)=P(\theta/\omega,y)\) for
\(\omega=2\pi/L\). If \(q_{i,n}\) is the \(n\)-th Fourier
coefficient on the standard \(2\pi\) circle, then
\[
 \mathcal E_b(P)=\sum_i\sum_{n\in\mathbb Z}
   e^{-b^2 n^2\omega^2/2}
       \langle q_{i,n}, T_\rho^y q_{i,n}\rangle_{\gamma_{d-1}} .
                                                               \tag{PT1}
\]
The wrapped Gaussian multiplier comes from its characteristic
function. The Gaussian transverse pairings are real and
nonnegative, and strictly positive for a nonzero coefficient
since all Hermite eigenvalues of \(T_\rho\) are positive.
It follows that
\(\mathcal E_\sigma(P)-\mathcal E_{\sigma/\rho}(P)>0\)
for every horizontally nonconstant profile. The exact
scaling identities for \(P^\uparrow(s,y)=P(s/\rho,y)\) are
\[
             \mathcal E_\sigma(P^\uparrow)
                       =\mathcal E_{\sigma/\rho}(P),
              \qquad m(P^\uparrow)=m(P).                       \tag{PT2}
\]
For example \(H_\sigma^s P^\uparrow(s,y)=
(H_{\sigma/\rho}^sP)(s/\rho,y)\) by changing the Gaussian
displacement. The horizontal mean is invariant under the
change of scale.

## Dissipation along a sequence of scaled charts

For \(t_j\to+\infty\), suppose that for every \(k\ge0\),
after subsequence selection,
\[
           1_A(\rho^k t_j+s,y)\ \longrightarrow\ P_k(s,y)
                   \qquad\text{in }L^1_{\rm loc}.              \tag{PT3}
\]
Local convergence and Gaussian tails imply uniform convergence
of the conditional scores on compact sets. For example, on
each such observation set the Gaussian kernel is bounded
above on a fixed source ball and the mass outside that ball
tends uniformly to zero as the source radius grows. Applying
this to the source chart \(k+1\) and passing the winner
inequalities to the limit gives
\[
    P_k^\uparrow(s,y)=P_k(s/\rho,y)
          \ \text{chooses a winner of }K_\sigma P_{k+1}(s,y)-\beta .
                                                               \tag{PT4}
\]
The same implication holds for a single pair of successive charts.
[Scaling-linked local winner equation (EI6)](THREE_CELL_EXTERIOR_INDEX_AND_DILATION_CHAIN.md)

Suppose all profiles in (PT3) are horizontally periodic, and
every label used on positive measure in \(P_k\) is present
on positive measure in \(P_{k+1}\). Then \(P_k^\uparrow\) has
every period of \(P_{k+1}\). Indeed restrict the winner scores
in (PT4) to the labels present in the source. Any two
distinct present labels have a nonconstant indicator contrast.
Its smoothed contrast is analytic and nonconstant
since \(K_\sigma P(s,y)=T_\rho P(s/\rho,y)\) and
\(T_\rho\) is injective. Subtracting a constant price difference
can therefore cause ties only on null sets. The unique
restricted winner inherits every source period. Restricting to
the present labels here also handles the possibility of two
absent source labels with identical adjusted scores. A source
with just one present label gives the period assertion directly.

Put \(a=P_k^\uparrow,\ b=P_{k+1}\) and integrate over any
horizontal period of \(b\) and the standard transverse
Gaussian measure. Both profiles have that period. The
pointwise winner equation gives
\(\langle a-b,K_\sigma b-\beta\rangle_L\ge0\).
Expanding the square shows
\[
 \begin{split}
 &\mathcal E_\sigma(a)-2\beta\cdot m(a)
              -\mathcal E_\sigma(b)+2\beta\cdot m(b)\\
 &\hspace{18mm}
   =2\langle a-b,K_\sigma b-\beta\rangle_L
                    +\langle a-b,K_\sigma(a-b)\rangle_L
   \ \ge\ \langle a-b,K_\sigma(a-b)\rangle_L .
 \end{split}
\]
Define \(F(P)=\mathcal E_\sigma(P)-2\beta\cdot m(P)\). By (PT2)
this yields
\[
      F(P_k)-F(P_{k+1})
          \ \ge\ \mathcal E_\sigma(P_k)
                         -\mathcal E_{\sigma/\rho}(P_k).
                                                               \tag{PT5}
\]
The right-hand side is nonnegative, and is strictly positive
for a horizontally nonconstant profile. Also
\(-2\max_i\beta_i\le F(P)\le1-2\min_i\beta_i\), by
\(0\le \mathcal E_\sigma(P)\le\|P\|_2^2=1\).
Consequently the losses on the right of (PT5) are summable
along every such chain. A uniformly positive loss in every
chart is impossible.

## A fixed periodic tail

Suppose a fixed profile \(P\) of horizontal period \(L\)
satisfies
\[
   1_A(t+s,y)-P(t+s,y)\longrightarrow0
        \quad\text{in }L^1_{\rm loc}(ds\,dy)
        \text{ as }t\to+\infty.                               \tag{PT6}
\]
Choose integers \(n_j\to\infty\) such that the distance of
\(\rho n_j\) to \(\mathbb Z\) tends to zero, and take
\(t_j=n_jL\). Existence of this sequence follows from recurrence
of a circle rotation, including its finite orbit in the rational
case. The observation and source chart limits are both \(P\).
Here translation continuity in local \(L^1\) applies to the
small source phase shift modulo \(L\). Apply (PT4)–(PT5) to
this one link, which has the same set of present labels
at both ends. Its left side vanishes, so (PT1) forces \(P\)
to be horizontally constant.

In particular a genuinely oscillatory periodic partition
cannot satisfy the one-dimensional Gaussian score winner
equation at a reciprocal integer correlation. The contradiction
uses the decrease of the horizontal Fourier energies when
the winner is pulled back from observation coordinates to
source coordinates.

## Slowly modulated stripes and compact variable templates

A more general assumption is that for all sufficiently large
\(t\), the translated partition is approximated in
\(L^1_{\rm loc}(ds\,dy)\), with error tending to zero as
\(t\to\infty\), by
\[
                  P_t(s,y)=q_t(\omega_t s+\vartheta_t,y).
                                                               \tag{PT7}
\]
Assume \(0<\omega_-\le\omega_t\le\omega_+<\infty\), and that
\(q_t\) belongs to a compact family of ordinary templates
in \(L^1(\mathbb T\times\gamma_{d-1})\). Assume each template
has every label present and is horizontally nonconstant.
There is no assumption that the masses of a label agree
between different templates.

Take any sequence \(t_j\to+\infty\). For each fixed \(k\ge0\),
pass to a subsequence on which the template, frequency and
phase at \(\rho^k t_j\) converge, then use a diagonal
subsequence. The corresponding periodic profiles converge
in \(L^1_{\rm loc}\). To see this when the templates are
merely measurable, first control the difference between two
templates integrated at the same affine phase on a bounded
coordinate interval. The pushforward of that interval to the
circle has bounded density uniformly for the frequencies in
\([\omega_-,\omega_+]\). Next approximate the limiting template
by a continuous one in \(L^1\). Convergence of the affine
phases then proves the claim on each bounded chart. Restricting
the transverse Gaussian to a bounded region permits the
equivalent Lebesgue bounds there. The approximation in (PT7)
therefore gives (PT3) with all labels present in every chart.

Formula (PT1) is continuous jointly in the template and
frequency on this compact range. Indicator \(L^1\)
convergence implies \(L^2\) convergence, the operators are
contractions, and the horizontal multipliers converge
uniformly over integer frequencies. Their high-index tails
are uniformly small since \(\omega_->0\). Thus the horizontal
loss \(\mathcal E_\sigma-\mathcal E_{\sigma/\rho}\) has a
positive minimum over the compact family and the frequency
interval, by its strict positivity in (PT1). This contradicts
summability in (PT5). Hence (PT7) cannot describe a stationary
tail under these hypotheses.

For example, fix a cyclic order of the three labels and
assign each to one nonempty angular interval in a partition
of the circle. Suppose that in a tail the labels are those
of the fixed template evaluated at a phase \(\theta(s)\),
independently of the transverse coordinate there, where
\(\theta\) is \(C^1\) and
\[
 \omega_-\le\theta'(s)\le\omega_+,\qquad
   \sup_{|u|\le R}|\theta'(s+u)-\theta'(s)|\longrightarrow0
                   \quad\text{for every fixed }R .
\]
Then \(\theta(t+u)-\theta(t)-\theta'(t)u\to0\) uniformly
for bounded \(u\). The phase approximation can change a
label only near a finite set of angular endpoints, and
the affine phase has frequency bounded below. The label error
on every bounded chart tends to zero, giving (PT7) and hence
a contradiction. Equal stripe widths are unnecessary.

For either application, exact agreement with the described
profile on \(\{s>S,\ |y|<w(s)\}\) suffices if \(w(s)\to\infty\).
Every translated fixed compact then eventually lies inside
this region, and discrepancies outside it are irrelevant to
(PT3). This covers a fixed opening wedge, and allows widths
growing sublinearly or errors tending to zero in translated
local measure. The argument leaves unbounded or vanishing
frequencies, noncompact template families and intermittent
oscillatory windows unresolved.

## Horizontal dissipation in a fixed strip

The same argument yields a statement for arbitrary stationary
winners, without assuming periodicity or a recurrent template.
For every positive integer \(j\),
\[
 \begin{split}
 &\frac1R\int_0^R\int_{\mathbb R^{d-1}}
       \sum_i|\partial_s^jT_\rho1_{A_i}(s,y)|^2
              \,d\gamma_{d-1}(y)\,ds \longrightarrow0,\\
 &\int_1^\infty\int_{\mathbb R^{d-1}}
       \sum_i|\partial_s^jT_\rho1_{A_i}(s,y)|^2
              \,d\gamma_{d-1}(y)\,\frac{ds}{s}<\infty .
 \end{split}                                                   \tag{PT8}
\]
This holds in every chosen outward direction. The weight is
\(ds/s\) on a fixed transverse strip, with the latter included
by restricting the transverse Gaussian measure. The related
radial estimate (OT5) uses the full-dimensional weight
\(|(s,y)|^{-d}\).
[Full-dimensional radial dissipation, (OT5)](GAUSSIAN_OUTWARD_THRESHOLD_LIMIT.md#energy-along-the-outward-orbit)

Write \(f=1_A\), \(a(s,y)=f(s/\rho,y)\) and
\(w_R=1_{[0,R]}\). Then \(a\) chooses a winner of
\(K_\sigma f-\beta\) on the full product with measure
\(ds\,d\gamma_{d-1}(y)\). Define
\[
 \begin{split}
 \mathcal E_{b,R}(g)&=R^{-1}\langle w_R g,K_b(w_Rg)\rangle,\\
 m_R(g)&=R^{-1}\int w_R g\,ds\,d\gamma_{d-1},\\
 F_R&=\mathcal E_{\sigma,R}(f)-2\beta\cdot m_R(f),\\
 \Delta_R&=\mathcal E_{\sigma,R}(f)-\mathcal E_{\sigma/\rho,R}(f)
                      \ \ge0 .
 \end{split}
\]
The nonnegativity follows by the real-line horizontal Fourier
transform and the positive transverse OU operator. All the
quadratic pairings here are on \(L^2(ds\,d\gamma_{d-1})\)
after applying the cutoffs.

In the winner inequality on \([0,R]\), replace \(K_\sigma f\)
by \(K_\sigma(w_R f)\). The integrated error before division
by \(R\) is at most a constant times
\[
 \int_0^R\int_{\mathbb R\setminus[0,R]}
                 \phi_\sigma(s-t)\,dt\,ds
       \le 2\int_0^\infty h\,\phi_\sigma(h)\,dh .
\]
This uses only boundedness of the label coordinates, the Markov
property in the transverse variables and the Gaussian density
in the first variable. Completing the square with the cutoff
now gives
\[
  \mathcal E_{\sigma,R}(a)-2\beta\cdot m_R(a)
       \ge \mathcal E_{\sigma,R}(f)-2\beta\cdot m_R(f)
                            -C_\rho/R .
\]
Changing both horizontal endpoints in the quadratic term,
\(s=\rho u,t=\rho v\), gives
\(\mathcal E_{\sigma,R}(a)=\mathcal E_{\sigma/\rho,R/\rho}(f)\).
The mass term scales in the same way. Consequently
\[
                 F_{R/\rho}-F_R
                       \ge \Delta_{R/\rho}-C_\rho/R .          \tag{PT9}
\]
As before \(F_R\) is bounded between \(-2\max\beta_i\) and
\(1-2\min\beta_i\). On any geometric grid
\(R_n=R_0\rho^{-n}\), nonnegativity of \(\Delta_R\) and
telescoping give \(\sum_{n\ge1}\Delta_{R_n}<\infty\).

In fact \(\Delta_R\to0\) without selecting a grid. For
\(S\ge R>0\), extending the interval to \(S\) changes either
unnormalized pairing in \(\Delta_R\) by at most \(2(S-R)\).
The new pairs have an endpoint in \([R,S]\) and the kernel has
unit row integrals. The change of normalization contributes
at most \((S-R)/S\), since the old unnormalized energy is at
most \(R\). This gives a uniform Lipschitz bound for \(\Delta_R\)
in \(\log R\). If \(\Delta_{r_l}\) stayed bounded away from
zero for some \(r_l\to\infty\), write \(r_l=t_l\rho^{-n_l}\)
with \(t_l\in[1,\rho^{-1}]\). Pass to \(t_l\to t_*\).
The Lipschitz bound would make \(\Delta_{t_*\rho^{-n_l}}\)
stay bounded away from zero, contrary to summability on that
fixed grid.

To obtain derivative control, for \(j\ge1\) compare multipliers
at horizontal frequency \(\xi\) and transverse Hermite degree
\(m\). The squared derivative multiplier of \(K_\sigma\) is
\(\xi^{2j}e^{-\sigma^2\xi^2}\rho^{2m}\). It is at most
a constant depending on \(\rho,j\) times
\[
         \big(e^{-\sigma^2\xi^2/2}
              -e^{-\sigma^2\xi^2/(2\rho^2)}\big)\rho^m .
\]
This holds near zero since the difference is comparable to
\(\xi^2\), at infinity by the strict Gaussian reserve, and on
the intervening compact region by continuity. Hence
\(R^{-1}\|\partial_s^j K_\sigma(w_R f)\|_2^2\le C_{\rho,j}\Delta_R\).
Removing the cutoff for output \(s\in[0,R]\) costs \(O_{\rho,j,k}(1/R)\)
in the normalized squared bound. Indeed the pointwise absolute
derivative error is bounded by
\(\int_{\mathbb R\setminus[0,R]}|\phi_\sigma^{(j)}(s-t)|\,dt\),
which is at most \(C_{\rho,j}(e^{-c s^2}+e^{-c(R-s)^2})\).
Its square integrates to \(O_{\rho,j}(1)\) uniformly in \(R\).

Put \(W_j(s)=\int\sum_i|\partial_s^jK_\sigma f(s,y)|^2\,d\gamma_{d-1}(y)\).
The previous bounds give
\[
               R^{-1}\int_0^R W_j(s)\,ds
                     \le C_{\rho,j,k}(\Delta_R+1/R).
\]
This proves convergence of that average to zero. On the grid
\(R_n=\rho^{-n}\), the averages are summable. For every \(s\ge1\)
the first grid point not below \(s\) is at most \(s/\rho\), so
\(\sum_{n:R_n\ge s}R_n^{-1}\ge\rho/s\).
Fubini therefore gives \(\int_1^\infty W_j(s)\,ds/s<\infty\).
Finally \(T_\rho f(s,y)=K_\sigma f(\rho s,y)\).
The chain rule and the change \(s\mapsto\rho s\) give both
parts of (PT8).

For a precise consequence about recurrent windows, fix a
nonconstant-in-\(s\) profile \(P\). Then
\(\partial_s K_\sigma P\) is nonzero somewhere. Otherwise
horizontal translations of \(P\) would have identical
smoothed profiles and injectivity would make \(P\) horizontally
constant. Choose fixed \(L<\infty\) and a bounded transverse box
on which its squared derivative integral is positive. By
Gaussian derivative tails and boundedness of the labels,
agreement with \(P\) to a sufficiently small fixed
\(L^1\) error on a sufficiently large fixed observation box
around a translated center \(t\) forces
\(\int_{t-L}^{t+L}W_1(s)\,ds\ge\epsilon\) for some fixed
\(\epsilon>0\). The constants depend on the chosen profile.
The set of such centers beyond \(2L+1\) has finite
\(\int dt/t\). Indeed, integrating the displayed window
energy with \(dt/t\) and using Fubini bounds it by a
constant times \(\int_1^\infty W_1(s)\,ds/s\). This excludes
positive-density returns to a fixed oscillatory profile
with this uniform local accuracy. It allows isolated windows
whose centers are arbitrarily far apart.
