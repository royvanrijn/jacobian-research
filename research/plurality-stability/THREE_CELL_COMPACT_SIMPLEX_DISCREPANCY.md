# Bounded components inside a simplex exterior

Fix an actual balanced planar Gaussian noise-stability maximizer \(A\)
at a scalar correlation in \((0,1)\). Suppose that \(A\) agrees with a
centered regular simplex \(T\) outside a compact set. If \(A\ne T\),
every label has a bounded strict-winner component containing changed
points. Moreover the same fixed partition can be stationary at only
finitely many positive correlations. The latter assertion concerns
fixed geometry. Shapes can vary between correlations.

## Every label must receive a bounded discrepancy

Write \(L_{ij}\) for the supporting line of the reference \(i,j\)
interface. At a correlation \(r\) where the fixed partition is
stationary, the \(i,j\) pressure holds on the exterior ray and so on
all of \(L_{ij}\) by analytic continuation. The reference contrast
vanishes there by reflection. The score of the compact difference
decays on the line as \(|x|\to\infty\), from the Gaussian conditional
density. Hence the pair pressure is zero. In coordinates \((s,v)\)
along and normal to \(L_{ij}\), injectivity of longitudinal Gaussian
convolution gives the slice constraint
\[
 \int_{\mathbb R}(1_{A_i}-1_{A_j})(s,v)
        e^{-v^2/(2(1-r^2))}\,dv=0
                                    \qquad\text{for a.e. }s .
                                                               \tag{CD1}
\]
[Weighted slice constraint (PG10)–(PG11)](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#moment-constraints-on-the-compact-discrepancy)

At the maximizing correlation, a strict-winner component of label
\(i\) avoids both reference pair lines incident to \(i\), since it
ties another score there. If the component is unbounded, it meets
the exterior of the compact exceptional set and must lie in the
same open region of those two lines as \(T_i\). Indeed, strict
winners are open and the ordinary labels agree outside the compact
set almost everywhere. Thus any point of the component outside a
large disk and away from the lines has the reference label \(i\).
The shared pair lines then prevent the component from leaving
\(T_i\). This conclusion concerns the strict-winner representative
of a cell, which agrees with its indicator almost everywhere.
[Geometric restrictions from exterior agreement](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#balanced-cuts-using-only-two-labels)

If \(A\ne T\), no single cell \(A_i\) can agree with \(T_i\) almost
everywhere, by single-simplex-cell rigidity. Exact balance therefore
gives \(\gamma_2(A_i\setminus T_i)>0\) for every \(i\). These incoming
changes lie entirely in bounded strict-winner components, modulo
null sets. In particular any such discrepancy requires at least
three bounded strict-winner components, including at least one
with each label. Thus one connected strict-winner set suffices
to force the simplex globally from exterior agreement. The
argument covers either orientation of the simplex germ at the
origin.
[One simplex cell determines the whole partition](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

When \(A\) ties \(T\), the positive incoming mass is uniform for each
label over compact positive-correlation bands. Otherwise there
would be tied pairs with one cell distance tending to zero on a
compact band. The simplex can be fixed by rotation and relabeling.
Strong maximizer compactness gives a limit maximizer sharing that
cell with \(T\), hence equal to \(T\). The uniform measurable local
simplex gap then excludes the sequence of nonsimplex ties.
[Compact-band single-cell exclusion](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#balanced-cuts-using-only-two-labels)

## A fixed compact discrepancy has finitely many stationary correlations

For each pair \(i,j\) and longitudinal slice \(s\) define
\[
 g_{ij,s}(v)=
  (1_{A_i}-1_{A_j})(s,v)+(1_{A_i}-1_{A_j})(s,-v),
                   \qquad v>0 .
\]
Each of these slice functions is supported in a fixed bounded
interval, since outside the exceptional disk the reference pair
contrast is odd in \(v\). If every such function vanished almost
everywhere in \((s,v)\), the reflection in \(L_{ij}\) would exchange
labels \(i,j\) and fix the third label for every pair. Indeed the
pair contrast takes the distinct values \(1,-1,0\) on the three
labels. The reflections generate the corresponding cyclic
rotation, so all three angular masses agree on almost every
circle. The balanced-circle comparison at the maximizing
correlation would give \(A=T\) up to a rigid motion and relabeling.
The fixed exterior then identifies the whole partition as \(T\).
[Equal angular masses force a maximizing simplex](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#circular-balance-and-compatibility-of-first-angular-moments)

Hence for a nontrivial compact discrepancy there is a pair with a
positive-measure set of \(s\) for which \(g_{ij,s}\ne0\). For almost
every such \(s\) it is a finite step function on its support. To
see this, restrict the three real-analytic adjusted score differences
at the original maximizing correlation to the slice. To avoid
identically zero slices, for each difference use one coefficient
of its Taylor expansion in \(v\) that is a nonzero analytic
function of \(s\). Its zero set is discrete. On every other slice
the three differences have finitely many zeros in the compact
interval, and their signs determine all the labels between zeros.

Suppose infinitely many distinct correlations were stationary for
this fixed \(A\), and choose a countable infinite subfamily.
Select a single generic \(s\) with \(g_{ij,s}\ne0\) satisfying
(CD1) for the whole subfamily by countable intersection. Let
\(m\) be the number of sign changes of its nonzero pieces, in order
of increasing \(v\). Then
\[
   \mathcal L(t)=\int_0^\infty g_{ij,s}(v)e^{-t v^2}\,dv
                                                               \tag{CD2}
\]
has at most \(m\) distinct real zeros. For \(m=0\) this follows
from strict positivity of the weight. For \(m>0\), choose \(c>0\)
at one transition in \(v^2\) between opposite signs, allowing
an intervening gap where \(g_{ij,s}=0\). The derivative of
\(e^{ct}\mathcal L(t)\), after division by \(e^{ct}\), is
the transform of \((c-v^2)g_{ij,s}(v)\), which has exactly
\(m-1\) sign changes. Rolle's theorem and induction prove
the bound. Each correlation in the subfamily supplies a
different zero by (CD1), at \(t=1/(2(1-r^2))\). This contradiction
proves finiteness. After finiteness is known, the same argument
on any slice satisfying (CD1) at all admissible correlations
gives the sign-change bound whenever its even discrepancy is
nonzero.

## Both cyclic orientations occur in the defect

Let \(w_i\) be the unit center of \(T_i\), and let \(e_i\) be the
outward ray direction of the \(i,i+1\) interface, with cyclic indices.
Write \(H_i=1_{A_i}-1_{T_i}\) and \(\mu=\gamma_{(1-r^2)I_2}\) at the
maximizing correlation. A changed point in \(T_i\) is a forward
transfer if its new label is \(i+1\), and a reverse transfer if its
new label is \(i-1\). Any nontrivial compact discrepancy has both
types on sets of positive measure.

For \(p>1\) define \(g_p(s)=s_+^p+2^p s_-^p\), where
\(s_+=\max(s,0)\) and \(s_-=\max(-s,0)\). The slice conditions give
\[
     \int \sum_i g_p(e_i\cdot x)(H_i-H_{i+1})(x)\,d\mu(x)=0.
                                                               \tag{CD3}
\]
At a point in the interior of \(T_i\), put \(s=e_i\cdot x\),
\(t=e_{i-1}\cdot x\), and \(u=e_{i+1}\cdot x=-w_i\cdot x\).
Then \(s+t+u=0\) and \(u<\min(s,t)\). In particular,
\[
  2g_p(s)<g_p(t)+g_p(u),\qquad 2g_p(t)<g_p(s)+g_p(u).
                                                               \tag{CD4}
\]
If \(s,t\ge0\), then \(g_p(u)=2^p(s+t)^p>2\max(g_p(s),g_p(t))\).
Otherwise exchange \(s,t\) as needed so \(s>0>t\). Then
\(s=|t|+|u|\) with \(|t|<|u|\). Strict convexity gives
\(g_p(t)+g_p(u)>2s^p\), and \(g_p(t)<g_p(s)<g_p(u)\)
gives the other inequality.

The integrand in (CD3) vanishes at unchanged points. At a forward
transfer it is \(g_p(t)+g_p(u)-2g_p(s)>0\), and at a reverse transfer
it is \(-[g_p(s)+g_p(u)-2g_p(t)]<0\). Thus the integrals of the
two positive weights over the two orientations agree for every
\(p>1\). Both orientations must occur if either occurs.

They also reach the same largest essential source-center projection
\(w_i\cdot x\), taking the source index \(i\) at each changed point.
Indeed in the interior of \(T_i\) one has
\(g_p(u)=(2w_i\cdot x)^p\) and each incident weight is the \(p\)-th
power of a strictly smaller nonnegative number. On either orientation,
the \(p\)-th root of its positive integral therefore tends to twice
the essential supremum of \(w_i\cdot x\) there. For the lower bound,
restrict to a positive-measure set on which the source projection is
within any fixed distance of its supremum and both incident ratios
are bounded strictly below one. Such a subset exists by a countable
exhaustion of the strict pointwise ratio bounds. The upper bound
follows by bounding the two positive terms by the supremum. Equality
of the two integrals proves equality of the projections.

## A spectator bisector forces the simplex

Fix a reference pair \(i,j\) and write \(k\) for the third label.
Use orthonormal coordinates \((s,v)\) along and normal to its reference
line, with \(v>0\) on the \(T_i\) side. Define the spectator imbalance
\[
 b_k=\gamma_2(A_k\cap\{v>0\})-\gamma_2(A_k\cap\{v<0\}).
\]
If \(b_k=0\), the maximizer is the simplex. Thus a nonsimplex compact
defect has nonzero imbalance for each pair, and each label's defect
has a nonzero net mass on both sides of its sector-center line.

Freeze the spectator and write \(f=1_{A_i}-1_{A_j}\) and
\(m(s,v)=\operatorname{sgn}(v)f(s,v)\). The stability terms depending
on the split of the two labels are exactly
\(\frac12\langle f,T_rf\rangle_{\gamma_2}\), since their union is fixed.
Put \(\sigma^2=1-r^2\), \(a=r/\sigma^2\), and
\(\mu=\gamma_{\sigma^2 I_2}\), with one-dimensional marginal \(\mu_1\).
Factoring the joint Gaussian density gives
\[
 \langle f,T_rf\rangle_{\gamma_2}
     =\sigma^2\iint e^{a(ss'+vv')}f(s,v)f(s',v')\,d\mu\,d\mu.
\]
Subtract the transverse constant term \(1\) from \(e^{avv'}\). Its
quadratic form is nonnegative, since the longitudinal kernel
\(e^{ass'}\) is positive semidefinite, and it vanishes for \(f\) by
(CD1), which sets \(\int f(s,v)d\mu_1(v)=0\) for almost every \(s\).
The remaining quadratic form is the pairing of \(m\) against the
kernel
\[
 K_+(s,v;s',v')=\sigma^2 e^{ass'}(e^{avv'}-1)
                                   \operatorname{sgn}(vv').
\]
This kernel is strictly positive when \(vv'\ne0\).
All forms are absolutely finite by the Gaussian bounds and \(r<1\).

Now orient labels \(i,j\) correctly throughout their union, replacing
\(f\) by \(f_*=\operatorname{sgn}(v)|f|\) and \(m\) by \(|m|\). Write
\(D=\{m=-1\}\) and \(P=\{m=1\}\). The stability of the oriented
partition is at least
\[
      S_r(A)+2\iint_{D\times P} K_+\,d\mu\,d\mu.          \tag{CD5}
\]
Indeed, the strictly positive kernel gains on the two cross terms
between \(D\) and \(P\), and the added transverse constant form of
\(f_*\) is nonnegative. Neither term uses an equality case of a
rearrangement theorem.

If \(b_k=0\), the oriented partition is exactly balanced. Exterior
agreement gives \(P\) positive measure, so maximality and (CD5) imply
that \(D\) is null. The pair is then affinely separated in the original
partition. The affine-separation theorem makes every strict-winner
component unbounded. By the compact discrepancy result above there
can be no changed points, so \(A=T\). Since \(T_k\) itself is bisected,
balance also shows that confining the entire discrepancy of label
\(k\) to one side would give \(b_k=0\).
[Affine separation and unbounded strict winners](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md)

Without a bisector, the oriented label masses are
\(1/3-b_k/2\) and \(1/3+b_k/2\). Repairing them moves label mass
\(|b_k|/2\) between the two labels and costs at most \(|b_k|\) in
stability by the union bound for the two observations. Hence
\[
           |b_k|\ge 2\iint_{D\times P}K_+\,d\mu\,d\mu.    \tag{CD6}
\]
For nonsimplex maximizers with all defects supported in a fixed disk
and correlations in a compact subinterval of \((0,1)\), each \(|b_k|\)
has a uniform positive lower bound. Otherwise fix the exterior
simplex by rotation and relabeling and take a strongly convergent
maximizer subsequence at a limiting correlation in that interval.
The exterior agreement, balance and vanishing imbalance pass to
the limit. The bisector result makes the limit the simplex.
The uniform measurable local simplex comparison then contradicts
a nonsimplex sequence of global maximizers tending to the simplex.
[Positive-correlation compactness](SMALL_POSITIVE_CORRELATION_INTERVAL.md#a-rigorous-continuation-boundary-in-correlation)
[Uniform measurable local comparison](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)

## Repairing the oriented split in the exterior

The oriented split used in (CD5) differs from \(T\) only inside
the same exceptional disk. Its volume correction can be made
at superlinear cost. For fixed \(r,R_0\) and any
\(0<\epsilon<r^2\), a maximizer with discrepancy supported
in \(B(0,R_0)\) satisfies
\[
       2\iint_{D\times P}K_+\,d\mu\,d\mu
               \le C_{r,R_0,\epsilon}|b_k|^{\,1+r^2-\epsilon}.
                                                               \tag{CD8}
\]
The constants can be uniform on compact positive-correlation
bands for a fixed disk and a fixed \(\epsilon<\inf r^2\).

Let \(B\) be the oriented partition and
\(\delta=|b_k|/2>0\). It suffices to handle small \(\delta\);
larger values follow from (CD6). On the common reference ray
use longitudinal position \(s\in[L,L+1/L]\) and normal
distance \(0<v<w\) on the overfull side. Take
\(L=\sigma\sqrt{2\log(1/\delta)}\). This rectangle lies
outside the exceptional disk for small \(\delta\). The
longitudinal Gaussian mass is bounded below by a constant
times \(e^{-L^2/2}/L\). Thus the width giving Gaussian mass
\(\delta\) can be chosen with
\(w\le C\delta L e^{L^2/2}=C L\delta^{r^2}\), in particular
with \(w<1\). Transfer the rectangle to the underfull label.

The reference pair score vanishes along the ray and its
transverse derivative is bounded at fixed \(r\). Its contrast
in the rectangle is therefore \(O(w)\). The score difference
between \(B\) and \(T\), whose source is supported in the
exceptional disk, is bounded there by the conditional
Gaussian probability of that disk times a fixed constant.
Projection along the ray bounds this by
\(C e^{-(rL-R_0)^2/(2\sigma^2)}\). Hence the linear stability
cost of the mass correction is at most
\[
       C\delta\bigl(w+e^{-(rL-R_0)^2/(2\sigma^2)}\bigr)
             \le C_{r,R_0,\epsilon}\delta^{\,1+r^2-\epsilon}.
\]
The quadratic term of the transfer in the stability expansion
is nonnegative by positive definiteness of \(T_r\). The corrected
competitor is balanced, so maximality of \(A\) and (CD5) give
(CD8).

## Projection constraints allow fractional asymmetry

Exact balance at the standard Gaussian and all identities (CD1)
permit a compact fractional discrepancy with three nonzero
spectator imbalances and unequal angular masses. This example
uses fractional assignments and makes no assumption of maximality.

Let \(n_1,n_2,n_3\) be the oriented normals to the three spectator
axes, related by rotations of \(2\pi/3\). Write
\[
 D=\partial_{n_1}\partial_{n_2}\partial_{n_3},\qquad
 P(x)=\prod_i(n_i\cdot x),\qquad
 H_\varphi(x)=e^{|x|^2/(2\sigma^2)} D\varphi(x)
\]
for any smooth compactly supported \(\varphi\). On each axis's
transverse slice, the Gaussian weight in (CD1) removes the
transverse part of the exponential. It leaves the integral of a
normal derivative of a compactly supported function, which is
zero on every slice. Thus every \(H_\varphi\) has zero projections
in all three directions. For a smooth radial function \(h(|x|^2)\),
\[
       D h(|x|^2)=8P(x)h'''(|x|^2),                         \tag{CD7}
\]
since \(n_i\cdot n_j=-1/2\) for \(i\ne j\) and
\(n_1+n_2+n_3=0\).

First take \(\varphi(x)=\psi(|x|^2)\), with \(\psi\ge0\) smooth and
supported in a small positive interval. Then \(H_\varphi\) is
odd under inversion, so its standard Gaussian mass is zero.
Its signed half-mass across \(n_i\cdot x=0\) is the same and
nonzero for each \(i\). To check this, orient the polar angle
from \(n_1\). The polynomial \(P\) equals
\(\frac14|x|^3\cos(3\theta)\), with a nonzero integral of its
angular factor over the half-plane. The remaining radial
factor is a constant times
\[
  \int_0^\infty \psi'''(s)s^{3/2}e^{\alpha s}\,ds,\qquad
                      \alpha=\frac{1}{2\sigma^2}-\frac12>0 .
\]
The third derivative of \(s^{3/2}e^{\alpha s}\) is negative
for sufficiently small positive \(s\). With nonzero \(\psi\)
supported there, integration by parts makes this integral
strictly positive. Denote this odd perturbation by \(H_{\rm o}\).

For a perturbation with nonzero angular mean, take instead
\(\varphi_i=P\psi_i(|x|^2)\), with nonnegative nonzero smooth
\(\psi_1,\psi_2\) supported in disjoint positive intervals.
Both \(H_{\varphi_i}\) are even. Their standard Gaussian masses
are proportional, by integration by parts and (CD7), to
\[
                 -8\alpha^3\int P(x)^2
                         \psi_i(|x|^2)e^{\alpha|x|^2}\,dx<0 .
\]
Take a nonzero linear combination \(H_{\rm e}\) with zero mass.
Its circular mean is nonzero on a positive measure of radii.
Indeed, testing \(D\varphi_i\) against \(|x|^6\) and integrating
by parts via (CD7) gives \(-48\int P^2\psi_i\ne0\).
This is a radial test, so \(H_{\varphi_i}\) has a nonzero
circular mean as well. Their supports are in disjoint annuli,
which precludes cancellation in the combination. The even
perturbation's signed half-masses are zero by inversion.

Choose a disk containing all these supports and use the
baseline tuple equal to \((1/3,1/3,1/3)\) inside it and to
the simplex indicators outside. It is balanced and obeys
(CD1), since each truncated exterior contrast is odd in
the transverse coordinate of its own axis. Add
\(\epsilon(1,1,-2)(H_{\rm o}+H_{\rm e})\) for sufficiently
small nonzero \(\epsilon\). The tuple still lies in the probability
simplex and is balanced. Its spectator imbalances are all
nonzero by the odd component, and its angular masses are
unequal on a positive measure of circles by the even component.
Thus the slice identities alone in the relaxed class do not
provide the bisector or the equal-circle comparison.
