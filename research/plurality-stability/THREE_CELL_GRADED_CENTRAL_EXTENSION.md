# Graded angular reserves for central circular mass profiles

The comparison strengthens the angular reserve using four central
neighborhoods and the pairs where both profiles are central.

The attempted cover of $[2641/4000,2661/4000]$ used the same overwritten
radial weight as the [signed-direction attempt](THREE_CELL_SIGNED_AXIS_EXTENSION.md).
The [corrected checker](check_three_cell_graded_central_extension.py)
rejects the first slab. The analytical angular reserves below remain
valid on their stated domains.

## One central partner and two central partners

Put $p_*=\mathbf1/3$, $r_*=\sqrt{2/3}$ and use the central predicate

\[
 G_h=\{t:|p(t)-p_*|\le(h-1/3)/r_*\},
 \qquad h\in\{7/20,3/8,2/5,5/12\}.                         \tag{GC1}
\]

On $G_h$ all coordinates lie in $[2/3-h,h]$, because each centered
coordinate is bounded by $r_*|p-p_*|$. In particular the same-dominant-label
products in (WV1) vanish when either partner is in $G_h$.

For the linear cosine kernel, put
$G(v)=(1-\cos(\pi v))/\pi^2$ and

\[
 R(v)=\frac{2[G(2/3+v)-G(2/3)-G'(2/3)v]}{v^2},\qquad
 R(0)=-1/2,
 \qquad C_*=-\frac3{2\pi^2}. 
\]

If one partner is in $G_h$, the displacements of $p(t)+p(s)$ from
$2p_*$ are between $-h$ and $5/6$. If there is one negative
displacement, use the isosceles coefficient $K_1=(2R(-h)+R(h/2))/3$.
The increasing isosceles ratio in (SP2) gives $K_1\ge-1/2$ and
bounds the ratio at each smaller isosceles displacement. Subtracting
the quadratic with coefficient $K_1$ leaves a concave function on
the positive displacement range, where the cosine is at most $-1/2$.
Equalizing the two positive coordinates increases that remainder,
so the isosceles bound applies to the original triple. If there
are two negative displacements, their
squared share is at most one half, and the positive remainder ratio
is at most $-12/25$ by the first-four-fifths split in (WV10). These
two cases bound the ratio by

\[
 K(h,u)=\max\left\{ \frac{2R(-h)+R(h/2)}3,
                           \frac{R(-h)+u}2\right\},
                          \qquad u=-12/25 .                \tag{GC2}
\]

The checker verifies $R(-h)>u$ for each of these cases. The bound
holds also for a zero displacement vector by continuity.
If both partners are in $G_h$, every displacement has magnitude at
most $d_h=2(h-1/3)\le1/6$. The same argument gives $K(d_h,-1/2)$,
since a positive displacement at most $1/6$ stays in the region
where the cosine is at most $-1/2$.

Thus the positive Taylor slack for the linear cosine kernel is at
least $c_1(h)|z+w|^2$ if either profile is central and
$c_2(h)|z+w|^2$ if both are central, with

\[
 c_1(h)=\frac{C_*-K(h,-12/25)}4,\qquad
 c_2(h)=\frac{C_*-K(d_h,-1/2)}4 .                              \tag{GC3}
\]

Directed evaluation of the trigonometric endpoint formulas gives
the following lower bounds. The checker uses the directed values
of (GC3) themselves.

| $h$ | lower $c_1(h)$ | lower $c_2(h)$ |
| --- | ---: | ---: |
| $7/20$ | $0.0411$ | $0.0831$ |
| $3/8$ | $0.0376$ | $0.0770$ |
| $2/5$ | $0.0342$ | $0.0706$ |
| $5/12$ | $0.0319$ | $0.0663$ |

Separation of the linear cosine term lifts these bounds with a
factor $ae^{-a}$ to the positive Taylor slack for the exponential
kernel. As in (WV11), the dominant-label terms charged in either
variance comparison vanish on the whole union of the two central
events. Thus this reserve is available after either the excess
quartic or median estimate. [Sharp remainders](THREE_CELL_SHARP_SCALAR_REMAINDERS.md)
[Weighted variance](THREE_CELL_WEIGHTED_ANGULAR_VARIANCE.md)

## Variance retained by the union

Use the damped measure $\mu_0$ of (SV1), put $X=p-p_*$ and let
$g=\mu_0(G_h)$. For independent $X,Y$, the contribution from the
union of the two central events is at least

\[
                    \frac{2g}{2-g}\operatorname{Var}_{\mu_0}(p).
                                                               \tag{GC4}
\]

Indeed, conditional means $m_G,m_B$ and variances $V_G,V_B$ on
$G_h$ and its complement give the union integral of $|X+Y|^2$ as

\[
 2g V_G+2g(1-g)V_B+4g^2|m_G|^2+2g(1-g)|m_G+m_B|^2.
\]

With $d_m=m_G-m_B$, completing the square in $m_G$ bounds the
last two terms below by $2g^2(1-g)|d_m|^2/(2-g)$. The other two
coefficients dominate the corresponding within-group weights of
$(2g/(2-g))\operatorname{Var}_{\mu_0}(p)$. The endpoint cases
$g=0,1$ follow directly. In addition the both-central part retains
at least $2g\int_{G_h}|z|^2d\mu_0$. It is used with the extra
coefficient $c_2-c_1$ after charging $c_1$ to the whole union.

## Selected trial bounds for the central set

Fix an amplitude bin with upper endpoint $A_+$. Let $\theta_+$ and
$r_\perp$ be the axial and orthogonal residual bounds of (FR3).
For a fraction $q\in\{1/10,\ldots,9/10\}$ let $S$ consist of anchor
intervals whose upper absolute trial value satisfies

\[
                  A_+(1+\theta_+)|\widehat\phi|
                           \le q (h-1/3)/r_* .               \tag{GC5}
\]

Denote the lower mass of $S$ by $m_S$, the lower physical trial
square integral over $S$ by $U_S$, and an upper density ratio on $S$
by $B_S$. These are cumulative bounds over the sorted anchor cells
as in (FR5)–(FR6). The checker includes the scale factor
$(1-u^2)/(1-l^2)$ in $U_S$ to transfer from anchor to physical
trial norms. The axial-adjusted projected reference has distance at
most the left side of (GC5) from the center. If an actual profile in
$S$ is outside $G_h$, its orthogonal residual must therefore have
size more than $(1-q)(h-1/3)/r_*$. In particular

\[
 g\ge g_-:=\left[m_S-B_S\left(
       \frac{A_+r_\perp}{(1-q)(h-1/3)/r_*}\right)^2\right]_+.
\]

The same exceptional mass, multiplied by the trial square upper
bound on $S$ from (GC5), gives

\[
 \int_{S\cap G_h}\phi^2d\mu_0\ge
        U_S-B_S r_\perp^2\left(\frac{q}{(1-q)(1+\theta_+)}\right)^2 .
\]

If the right side is a positive value $u_S$, the original
decomposition $z=\alpha\phi+g$ and the triangle inequality on
$S\cap G_h$ give

\[
        \int_{G_h}|z|^2d\mu_0
              \ge A^2(\sqrt{u_S}-\sqrt{B_S}\,r_+)_+^2 .     \tag{GC6}
\]

Here $r_+$ is the original residual-ratio bound. The orthogonal
residual was used for the central-set test, and the original
residual is used for the energy comparison to the unclipped trial.

For each $h,q$ the checker multiplies the sum of the one-sided
union bound and the both-central gain by the common damped
coefficient $(\gamma/d)Z_0^2$. It takes the best of these finite
choices and zero. The angular gain adds to the maximum of the
remaining quartic and median alternatives from the signed-axis
comparison. No actual first angular moment deficit is charged.

## Attempted directed cover

The attempted cover uses 20 rational slabs of width $1/4000$, 128
residual bins, 24 label-direction bins and 12 amplitude bins.
The spectral block uses the preceding checker’s 120-bit Arb,
168 kernel degrees, 54 row terms and 4096 radial panels through
radius 14. The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_graded_central_extension.py

## Combining graded central angular reserves

Let \(Z,W\) be independent copies of the centered circular-mass profile
under the damped first-feature law \(\mu_0\) of (SV1). For the nested sets
\(G_1\subset\cdots\subset G_m\) from (GC1), write \(U_j\) for the event
that at least one partner is in \(G_j\) and \(B_j\) for the event that
both are. Let \(a_j=c_1(h_j)\), \(b_j=c_2(h_j)\), with zero appended
at \(m+1\). The ordered choices in (GC3) have \(a_j,b_j\) decreasing,
as verified by the directed trigonometric evaluations. Nonnegative
decreasing lower values of their exact coefficients also suffice.

After factoring out the common damped kernel, (GC2)--(GC3) give a
positive Taylor slack of at least \(a_j|Z+W|^2\) on \(U_j\) and
\(b_j|Z+W|^2\) on \(B_j\). The proof also applies when extra thresholds
with those same one- and two-partner bounds are available.

### Simultaneous grading

For each anchor index \(s\) define \(d_j=(b_j-a_s)_+\) including
\(d_{m+1}=0\). There is the pointwise lower bound
\[
 \left(\sum_{j=s}^m(a_j-a_{j+1})1_{U_j}
       +\sum_{j=1}^m(d_j-d_{j+1})1_{B_j}\right)|Z+W|^2 .
 \tag{SG1}
\]
Indeed, the first sum is either zero or \(a_j\) for a threshold whose
union event holds. It is at most \(a_s\). If the second sum is positive,
it equals \(b_k-a_s\) for the first true both-central event, so the
total coefficient is at most \(b_k\). If it is zero, the first sum
alone is available.

For \(s=j\) the expression in (SG1) dominates the earlier
\([a_j1_{U_j}+(b_j-a_j)1_{B_j}]|Z+W|^2\) since its first sum equals
\(a_j\) on \(U_j\) and its second sum is at least \(b_j-a_j\) on
\(B_j\). Outside \(U_j\) every summand in (SG1) is nonnegative.
Each integral in (SG1) therefore uses the best of all its own
(GC4)--(GC6) bounds. In particular the anchor-cell cutoffs in (GC5)
can be chosen independently for each event, and independently for
the union and both-central estimates at the same threshold.

### Exact finite variance test for all grades

Write \(I\) for the index of the first central threshold met by a
profile, with \(I=m+1\) if no threshold is met, and similarly \(J\)
for its independent partner. The preceding one- and two-partner
estimates give the pointwise weight
\[
 H_{ij}=\max\{a_{\min(i,j)},b_{\max(i,j)}\}.
\]
The full weighted pair integral admits a small matrix comparison
that uses the whole weight \(H\). Suppose \(\pi_i=\Pr(I=i)\),
\(D=\operatorname{diag}(\pi)\), and set
\[
 M=\operatorname{diag}(D H\pi)+DHD,\qquad
 C=D-\pi\pi^{\mathsf T}.
\]
For any \(\lambda\ge0\) with
\[
 2(H\pi)_i\ge\lambda\quad(\pi_i>0),
 \qquad
 2M-\lambda C\succeq0,
 \tag{SG2}
\]
one has
\[
 \mathbb E H_{IJ}|Z+W|^2\ge\lambda\operatorname{Var}(Z).
 \tag{SG3}
\]
To verify this, set \(m_i=\mathbb E[Z\mid I=i]\) and
\(v_i=\operatorname{Var}(Z\mid I=i)\) where the stratum has positive
probability. The conditional expected squared norm of \(Z+W\)
in strata \(i,j\) is \(v_i+v_j+|m_i+m_j|^2\). Summation gives
\(2\sum_i\pi_i(H\pi)_i v_i+2m^{\mathsf T}Mm\). The total variance is
\(\sum_i\pi_i v_i+m^{\mathsf T}Cm\). The quadratic forms act separately
on the label coordinates. These two identities prove (SG3).
Condition (SG2) characterizes the best variance coefficient that
can be inferred from the masses alone if the conditional means
and variances are left free. Take zero means and variance in a
single class to test the first condition, or constant variables
within each class to test the second.

Directed lower bounds \(g_i\le\Pr(I\le i)\) suffice in (SG2), after
taking nondecreasing \(g_i\) and \(g_{m+1}=1\). Assign an auxiliary
uniform coordinate \(u\) by sampling uniformly on the actual
cumulative-probability interval of the observed grade. The
unconditional \(u\) is uniform. Reassign the grade to the first index
with \(u\le g_i\), which is at least the original grade
almost surely and has the chosen masses \(\pi_i=g_i-g_{i-1}\).
The matrix weight \(H\) is nonincreasing in both indices. Applying
this downgrade with independent auxiliaries to the two profiles
can only lower the weighted integral. The variance of the
profiles is unchanged, so (SG3) holds with the chosen masses.

The last matrix test and (SG1) apply to the same reserve. When using
both one should take the stronger lower bound.
