# Conditional rigidity from two obtuse isolated ends

Fix \(0<r<1\). Suppose a balanced planar Gaussian noise-stability
maximizer has bounded isolated full-half-plane graph ends for two
distinct pairs, in nonantipodal outward directions at a strictly
obtuse angle. Let \(A_0\) be their common cell and \(W\) the open wedge
between their outward rays from the origin. If \(A_0\) agrees with
\(W\) outside a compact set up to null sets, the partition has a
centered regular simplex exterior. If the common strict winner is
connected, the entire partition is a centered simplex. Neither
conclusion requires any initial description of the third interface.

## The backward slice identity

The two ends are exactly straight and centered, and the label
pressures are equal, by the two-end pressure result. Thus the cells
agree a.e. with the strict winners of the conditional scores
\(U_i=T_r1_{A_i}\). For distinct labels the contrast is a nonzero
real analytic function by injectivity of Gaussian convolution.
[Two ends and zero pressures](THREE_CELL_ONE_STRAIGHT_END_REFLECTION.md#two-isolated-ends-and-a-second-order-tie)

Write \(\theta\) for the outward angle, \(a=\theta/2\),
\(c=\cos a\), \(d=\sin a\), and take \(u_+=(c,d)\), \(u_-=(c,-d)\),
\(n_+=(-d,c)\). The noncommon labels of the upper and lower ends
are \(+\) and \(-\). In the upper frame \(x=s u_++v n_+\), put
\(\sigma^2=1-r^2\), \(\kappa=\sin\theta/(-\cos\theta)\). The lower
forward half-plane excludes \(+\) whenever
\[
 v<t/\kappa-S_-/\sin\theta,\qquad s=-t,
\tag{CE1}
\]
where \(S_-\) is its longitudinal cutoff. The analytic \(0,+\)
contrast vanishes on the whole upper supporting line. Injectivity
of longitudinal Gaussian convolution gives equality of the
variance-\(\sigma^2\) transverse masses of the two labels on a.e.
slice. At sufficiently negative \(s=-t\) the common slice is the
wedge slice \(v<-\kappa t\). Thus for a.e. late \(t\)
\[
 \int 1_{A_+}(-t,v)p_\sigma(v)\,dv=\Phi(-\kappa t/\sigma).
\tag{CE2}
\]
Here \(p_\sigma\) is the centered normal density of variance
\(\sigma^2\) and \(\Phi\) the standard normal CDF. Reflection of the
coordinate frame and interchange of the two noncommon labels give
the analogous lower-end formulas. Comparing the Gaussian tail rate
in (CE2) with the tail allowed by (CE1) gives \(\kappa\ge1\), hence
\(\theta\le3\pi/4\).
[The Gaussian slice identity and angle bound](THREE_CELL_TWO_END_ANGULAR_COMPATIBILITY.md#the-common-phase-and-the-two-slice-constraints)

The connected-common-winner hypothesis implies the wedge hypothesis.
The common winner avoids both supporting lines and meets their
forward overlap in \(W\), so by connectedness it is confined to \(W\).
Conversely each point of \(W\) sufficiently far out is in at least
one defining forward half-plane, since one of its outward projections
is at least \(c\) times its norm. There it belongs to the common cell
a.e. by the two end descriptions.

## Conditional comparison with a half-plane

Let \(B=\{(s',v'):\ell(s',v')>0\}\) with \(\ell=v'+\kappa s'\).
At observation \(p=(-t,v)\) write \(P_p\) for the source law
\(N(rp,\sigma^2 I)\). Take \(S_0\) so (CE2) holds a.e. on
\(s'\le -S_0\). For \(v>0\) and such a transverse slice, subtract
from the tilted density \(p_\sigma(v'-rv)\) a multiple of
\(p_\sigma(v')\) chosen so that the two densities agree at
\(v'=-\kappa s'\). By (CE2) this leaves the integral against
\(1_B-1_{A_+}\) unchanged. The new density difference times this
indicator difference is pointwise nonnegative, by monotonicity
of the likelihood ratio. On \(A_-\cap\{\ell\ge L\}\) it is at least
\((1-e^{-rvL/\sigma^2})p_\sigma(v'-rv)\), for \(L>0\). On the other
slices the original product has its full positive contribution on
this same set and can only be negative where \(v'<\kappa S_0\).
Consequently
\[
 \begin{split}
 U_+(p)&\le P_p(B)+P_p(s'>-S_0),\\
 U_+(p)&\le P_p(B)-\tfrac12 P_p(A_-\cap\{\ell\ge L\})
              +P_p(s'>-S_0,\ v'<\kappa S_0)
        \quad\text{if }rvL\ge\sigma^2\log 2,\\
 P_p(B)&=\Phi\!\left(\frac{r(v-\kappa t)}{\sigma\sqrt{1+\kappa^2}}\right).
 \end{split}
\tag{CE3}
\]

For any fixed \(0<\varepsilon<\kappa\), label \(+\) occurs on
positive transverse measure in \(v\in[(\kappa-\varepsilon)t,
(\kappa+\varepsilon)t]\) on a.e. sufficiently late slice \(s=-t\).
Indeed (CE1) excludes nonpositive \(v\), and the first bound in
(CE3) gives \(U_+\to0\) uniformly for \(0<v\le(\kappa-\varepsilon)t\).
The largest score is at least \(1/3\), so \(+\) is excluded there.
The tail above \((\kappa+\varepsilon)t\) cannot supply (CE2) on its
own. The reflected-coordinate proof gives the lower concentration.

For any \(p\) with \(p\cdot u_+\le0\) and \(p\cdot u_-\le S_-\),
one has \(p\cdot z\le C|z|\) on \(W\): write
\(z=\alpha u_++\beta u_-\) with \(\alpha,\beta\ge0\) and
\(c(\alpha+\beta)\le|z|\). The common cell lies in \(W\) and a
fixed disk up to null sets. Expanding the conditional Gaussian
density at its mean gives constants \(C_1,C_2\) such that
\[
 U_0(p)\le C_1\exp\!\left(-\frac{r^2|p|^2}{2\sigma^2}+C_2|p|\right).
\tag{CE4}
\]
The wedge contribution uses integrability of
\(e^{C'|z|-|z|^2/(2\sigma^2)}\), and the disk contribution uses
\(p\cdot z\le C''|p|\).

## The angle is exactly \(2\pi/3\)

If \(\theta<2\pi/3\), then \(\tan a<\kappa\). On the negative
horizontal bisector the ratio \(v/t\) is \(\tan a\) in both
frames. By (CE3) and (CE4), all three conditional scores
tend to zero there. Their sum is one, a contradiction.

Suppose \(2\pi/3<\theta\le3\pi/4\). Let \(d_+\) be the unit
normalization of \(-u_++\kappa n_+\) and \(d_-\) its horizontal
reflection. These are the directions of concentration for the
noncommon labels. In ordinary polar angles they are \(3a\) and
\(2\pi-3a\), separated by \(\Delta=3\theta-2\pi\in(0,\pi/4]\).
The positive normal of \(\ell\) has angle \(3a-\pi/2\), so
\(\ell(d_-)>0\). Consider any \(p=(-t,v)\) in the upper frame
with \(\kappa t-1\le v\le\kappa t\). Then
\(p=|p|d_++O(1)\). For each small fixed \(\eta>0\) and all large
\(t\), lower-frame concentration provides a point \(q\) where
label \(-\) wins and
\[
 |q-|p|\cos\Delta\,d_-|\le\eta|p|.
\]
Indeed a slice whose longitudinal coordinate is within one of
the target can be chosen from the full-measure set of late
slices, and the relative width of the concentration band can
be chosen sufficiently small in terms of \(\eta\). Since
\(U_-(q)\ge1/3\), a fixed radius \(R\) with
\(P_q(|z-rq|>R)<1/6\) gives \(A_-\cap B(rq,R)\) uniformly
positive Lebesgue measure by the conditional density upper
bound. For small enough \(\eta\), throughout that ball
\[
 \ell(z)\ge C_3|p|,\qquad
 |z-rp|\le r(\sin\Delta+\eta)|p|+O(1).
\]
The first inequality uses \(\ell(d_-)>0\) and \(\cos\Delta>0\).
Since \(\sin\Delta<1\), the density lower bound on the ball gives
\[
 P_p(A_-\cap\{\ell\ge C_3|p|\})
 \ge C_4\exp\!\left(-\xi\frac{r^2|p|^2}{2\sigma^2}\right)
 \qquad\text{for some }\xi<1.
\tag{CE5}
\]

For points of the band not already excluded by (CE1), (CE4)
applies. The early-slice error in the second bound of (CE3) is
\(O(\exp(-r^2|p|^2/(2\sigma^2)+C|p|))\), by independence and
the two one-dimensional Gaussian tail bounds. With
\(L=C_3|p|\), the loss in (CE5) dominates this error and
\(U_0(p)\). Since \(P_p(B)\le1/2\) on the band, it follows that
\(2U_+(p)+U_0(p)<1\), so \(-\) beats \(+\). For
\(v<\kappa t-1\) not excluded by (CE1), one has \(v>0\) for late
\(t\). Now \(P_p(B)\) is uniformly below \(1/2\), and (CE3)
and (CE4) give the same strict exclusion. At \(v=\kappa t\),
either the strict score inequality or the open half-plane
condition in (CE1) holds, so \(+\) is excluded on a nonempty
interval above \(\kappa t\) as well. The remaining slice has
mass strictly smaller than (CE2), a contradiction.

## Exterior and global rigidity

At \(\theta=2\pi/3\), compare the partition with the simplex
\(T\) that has common cell \(W\) and divides the other two labels
by the horizontal axis, with \(+\) on top. Use usual coordinates
\((X,Y)\), longitudinal coordinates \(s_\pm=cX\pm dY\), and
the centered Gaussian measure \(\mu=N(0,\sigma^2 I)\).
Write \(f=1_{A_+}-1_{T_+}\) and \(w(s)=\min(s+S,0)\), with
\(S\) larger than the radius of the exceptional disk. Subtract
the longitudinal slice identities of \(A\) and \(T\) on each of
the two supporting lines, then integrate them against \(w\).
On the support of either weight, the common cells agree. The
upper slice identity then concerns \(f\) and the lower one
concerns \(-f\), so
\[
 \int (w(s_+)-w(s_-)) f\,d\mu=0.
\tag{CE6}
\]
Where a weight is nonzero we have \(Yf\le0\) because \(T\)
uses the upper label above the axis outside the common cell.
The weight difference has the sign of \(Y\), strictly for
\(Y\ne0\) and \(\min(s_+,s_-)<-S\), since \(w\) is increasing
strictly below \(-S\). Hence (CE6) forces agreement there.
The original two forward half-planes have the prescribed
simplex labeling, and the region omitted by these four
half-planes is compact because the two longitudinal coordinates
are independent. This proves simplex exterior rigidity.

If the common strict winner is connected, it is confined to
\(W\) globally. Its Gaussian mass \(1/3\) is the full mass of
\(W\) at this angle. One simplex cell in a maximizer determines
the whole simplex, up to rotation, relabeling and null sets.
[One simplex cell determines the whole partition](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)
