# Localized quartic bounds for the three-cell comparison

The centered covariance and localized quartic estimates below give
sufficient bounds for comparing an adverse radial mass profile to the
simplex. The corrected checker fails on the first slab of the attempted
interval $[2661/4000,2681/4000]$. The previous passing receipt is invalid
because an amplitude-bin counter overwrote the origin weight entering
later amplitude caps. The estimates here supply no new correlation
cutoff. [Corrected checker](check_three_cell_localized_quartic_extension.py)

## The trial's centered covariance slope

For any of the tilted laws $\mu$ used in the comparison, let
$W=d\mu/(J\,d\nu)$ in the physical mixed weighted space. The bounded
positive operator $\mathcal T_\mu$ representing tilted variance has

\[
                   \mathcal T_\mu\phi=W(\phi-\mathbb E_\mu\phi).
\]

For the unit reference $\phi$ in this space, write
$V_\phi=\operatorname{Var}_\mu(\phi)$. Its variance form on
$z=\alpha\phi+g$ has the lower bound

\[
 \frac{\operatorname{Var}_\mu(p)}{A^2}
 \ge (\sqrt{V_\phi}-s_\mu r)_+^2,\qquad
 s_\mu=\frac{\|P_{\phi^\perp}\mathcal T_\mu\phi\|_J}{\sqrt{V_\phi}},
 \qquad r=G/A .                                               \tag{LQ1}
\]

This follows by projecting the centered form onto the unit vector
in the direction of $(\alpha/A)\mathcal T_\mu^{1/2}\phi$.
The orthogonality of $g$ to $\phi$ allows the projection in
$s_\mu$. In place of the universal bound $s_\mu\le\sqrt{B-V_\phi}$
from (SR15), the checker also uses the identity

\[
 \|P_{\phi^\perp}\mathcal T_\mu\phi\|_J^2
       =\mathbb E_\mu[W(\phi-\mathbb E_\mu\phi)^2]-V_\phi^2 . \tag{LQ2}
\]

To enclose this identity, use anchor $x=\sqrt{t_0}$ at slab upper
endpoint $u$. Denote $d_u=1-u^2$, $\beta_u=1/d_u$ and
$\kappa_u=u^2/d_u$. The mixed weighted measure has density
$x^2Q(x)e^{-\beta_u x^2}/(2d_u)$, with $Q$ its row polynomial.
For the anchor trial $\widehat\phi$, each tilted law has density
$c_\lambda x^{n+1}e^{-\lambda x^2}$ with
$c_\lambda=2\lambda^{1+n/2}/\Gamma(1+n/2)$.
Here $n=1,3$ and $\lambda=\beta_u$ for $\mu_1,\mu_3$.
For $\mu_0$ take $n=1$ and enclose $\lambda=(1+\rho)\beta_u$
on the slab. Define

\[
 W_\lambda(x)=2d_uc_\lambda x^{n-1}e^{-(\lambda-\beta_u)x^2}/Q(x),
 \qquad m_\lambda=\mathbb E_\lambda\widehat\phi .
\]

The physical trial is $\phi=\sqrt{d_u/d}\,\widehat\phi$ and the
physical density ratio is $W=(d_u/d)W_\lambda$. Thus variance
scales by $d_u/d$, while the squared norm in (LQ2) before subtracting
variance squared scales by $(d_u/d)^2\le1$. The directed upper
bound on the latter uses the anchor integral

\[
 \int c_\lambda x^{n-1}e^{-\lambda x^2}W_\lambda(x)
         (x\widehat\phi(x^2)-m_\lambda x)^2dx,                 \tag{LQ3}
\]

which is regular at the origin. The checker integrates it and the
mean on $[0,12]$ over 3600 intervals. For the tails, use $Q(x)\ge Q(12)$
and the polynomial bound of the raw tilted moment in the
excess-purity checker. In the square expansion the two Gaussian
rates are $\lambda$ and $\lambda-2\kappa_u$, both positive on the
stated range. The omitted part of the mean is bounded by the square
root of tail mass times tail raw second moment. With that full mean
enclosure, the tail of (LQ3) is at most
$2B[\mathbb E_\lambda(\widehat\phi^2 1_{x>12})+m_\lambda^2\mu_\lambda(x>12)]$.
Subtract a lower physical variance squared from the upper (LQ3)
and divide the square root by a lower $\sqrt{V_\phi}$. The checker
takes the smaller valid slope from (LQ2) and (SR15).

## Localizing the quartic difference on central profiles

Put $D=|p-q|^2$ for independent profiles under $\mu_1$ and put
$\mathcal Q_1=\mathbb E[D(1-D/2)]$. For $\tau>0$ let
$G_\tau=\{|p-p_*|\le\tau\}$, set $g=\mu_1(G_\tau)$ and write
$V_G=\operatorname{Var}_{\mu_1}(p\mid G_\tau)$. If at least one
partner is in $G_\tau$, then $D\le(\sqrt{2/3}+\tau)^2$. If both
are in $G_\tau$, then $D\le4\tau^2$. Write
$C_\tau=1-(\sqrt{2/3}+\tau)^2/2$ and $C_2=1-2\tau^2$.
The paired difference on the union of the two central events has
expectation at least $2g V_1$. Indeed it has the exact decomposition
$2g V_G+2g(1-g)V_B+2g(1-g)|m_G-m_B|^2$ in terms of the conditional
means and variances. On the both-central event its expectation is
$2g^2 V_G$. Hence

\[
           \mathcal Q_1\ge2g C_\tau V_1+2g(C_2-C_\tau)gV_G . \tag{LQ4}
\]

This retains the distinction between pairs with one central partner
and pairs with two. Subtract $(12/25)V_1$ in (LQ4) for the degree-one
variance already absorbed into the scalar form.

To bound $gV_G$ from below, take a residual and amplitude bin and
use the axial-adjusted reference and its orthogonal budget from (FR2)
and (FR3). For $0<q<1$ choose the anchor cells $S$ on which the trial
has certified magnitude at most $q\tau/[A_+(1+\theta_+)]$.
Let $m_S$ be their lower mass, $I_S$ their lower centered trial
energy, and $B_S$ an upper tilted density ratio there. A failure of
central membership in $S$ requires orthogonal residual larger than
$(1-q)\tau$, so

\[
  \mu_1(S\setminus G_\tau)\le b:=B_S(A_+r_\perp/[(1-q)\tau])^2,
  \qquad g\ge g_-=m_S-b .
\]

Let $c_S$ be the mean of $\phi$ on $S\cap G_\tau$. If $g_->0$ then
$|c_S|\le\sqrt{M_\phi/g_-}$ from the upper raw trial moment.
Subtracting the removed cells from the trial centered energy gives

\[
 \int_{S\cap G_\tau}(\phi-c_S)^2d\mu_1
 \ge I_S-b\left[
       \frac{q\tau}{A_+(1+\theta_+)}+\sqrt{M_\phi/g_-}
                        \right]^2.                          \tag{LQ5}
\]

When the right side is positive, say $I_-$, restrict the actual
centered energy to $S\cap G_\tau$ and use the triangle inequality,
optimizing over its centering constant. This gives

\[
             g V_G\ge A^2(\sqrt{I_-}-\sqrt{B_S}\,r_+)_+^2 .  \tag{LQ6}
\]

The error norm here is the original residual $G$ because the trial
in (LQ5) is the original affine trial. The orthogonal residual budget
enters only in certifying central membership. The checker encloses
the selected mass, raw trial moment, trial mean and density ratio
by cumulative interval sums. It uses five fixed radii $\tau$ from
$1/20$ to $3/20$ and five fractions $q$ from $1/5$ to $3/5$.

## Combining the reserves and corrected replay

For each amplitude, residual and direction bin, the checker takes
the maximum of the excess-quartic bound (EP4), the median bound
(PT5)–(PT6) after the same credit, and the localized quartic bound
from (LQ4)–(LQ6). It uses (LQ1)–(LQ3) for the required variances.
The graded positive central reserve adds to each alternative because
the same-dominant-label products vanish on its union and the
localized estimate only uses the negative Taylor remainder.
[Graded reserve](THREE_CELL_GRADED_CENTRAL_EXTENSION.md)

The checker attempts 20 adjacent correlation slabs of width $1/4000$,
128 residual bins, 24 direction bins and 12 amplitude bins when
needed. The finite spectral block is unchanged, with 120-bit Arb,
168 kernel degrees, 54 row terms and 4096 radial panels through
radius 14. The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_localized_quartic_extension.py

On the corrected run, $[2661/4000,1331/2000]$ has an upper gap of
approximately $0.05581$ in the scalar-versus-remainder test. The
amplitude-bin upper endpoint is approximately $0.14157$, and the
residual-ratio upper endpoint is approximately $0.28770$. A positive
gap here leaves the sufficient test unresolved for that slab.
