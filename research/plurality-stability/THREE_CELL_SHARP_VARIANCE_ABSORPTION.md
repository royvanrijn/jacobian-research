# Unconditional variance absorption for the sharp arc bound

The sharp aligned-arc scalar comparison admits an explicit radial
variance penalty for every mass profile. It uses the positive Taylor
slack for two profiles favoring the same label and the quartic slack
for profiles that differ. Both can be retained along with the
actual first angular deficit. In particular the penalty can be included
in the linear norm test even when the uncorrected scalar operator has
an adverse direction. The reserved bootstrap below proves the global
three-cell comparison through 0.648 and gives a positive unevaluated
right neighborhood. The argument uses circular masses and exact score
transfer cost in place of a global first-moment loss.

Use the notation of the [sharp arc bound](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)
and its [remainders](THREE_CELL_SHARP_SCALAR_REMAINDERS.md).
For fixed \(\rho\) let
\(d=1-\rho^2,\ \kappa=\rho^2/d,\ \gamma=2\rho/d\) and
\(d\nu(t)=e^{-t}dt\). Write \(a=\gamma\sqrt{ts}\), and use
\(\operatorname{Var}_\mu(p)=\int |p-\int p\,d\mu|^2d\mu\).
Define two normalized radial measures
\[
 \mu_0=\frac{f\,\nu}{Z_0},\quad f(t)=\sqrt t\,e^{-(\kappa+\gamma/2)t},
 \quad Z_0=\Gamma(3/2)(1-\rho)^{3/2},\qquad
 \mu_1=\frac{f_1\,\nu}{Z_1},\quad f_1(t)=\sqrt t\,e^{-\kappa t},
 \quad Z_1=\Gamma(3/2)d^{3/2}.
\]
One may subtract either of the following penalties from the sharp
scalar bound:
\[
 \beta_0\operatorname{Var}_{\mu_0}(p),\qquad
 \beta_0=\frac{\rho(1-\rho)^3}{4\pi d^2},\qquad
 \beta_1\operatorname{Var}_{\mu_1}(p),\qquad
 \beta_1=\frac{6\rho d}{25\pi}.                               \tag{SV1}
\]
Any convex combination is valid. These estimates are for the same
angular remainder, so their sum does not follow from (SV1).

## Retaining the same-label cosine slack

Put \(s=p+q\) for two circular mass profiles. For the linear cosine
kernel, its contribution to the *positive* Taylor slack in the
aligned-arc pairing is
\[
 H(s)=\frac{1}{2\pi^2}\left[\sum_i\cos(\pi s_i)+\frac32
                  -\frac34|s-\mathbf{2}/3|^2\right].
\]
This is nonnegative by (SP2) for the linear kernel. If \(s_i\ge1\),
write \(x=2-s_i\), so \(s_j+s_k=x\in[0,1]\). The cosine
addition identity gives
\(\cos(\pi s_j)+\cos(\pi s_k)\ge1+\cos(\pi x)\), and the square
sum is at most \(8/3-4x+2x^2\). Thus the bracket is at least
\(2(1-\cos(\pi(1-x)))-3(1-x)^2/2\).
The bound \(\sin(\pi u/2)\ge u\) on \([0,1]\) then gives
\[
 H(p+q)\ge \frac5{\pi^2}\sum_i(p_i-1/2)_+(q_i-1/2)_+.
\]
Indeed both positive parts can occur for at most one common label
and \((r+s)^2\ge4rs\) for \(r,s\ge0\).

For the exponential kernel the positive Taylor slack is at least
\(ae^{-a}H\). The sharp Taylor estimate is linear in the
kernel and holds for the increasing convex remainder
\(e^{ac}-ae^{-a}c\), as proved before (SR4).
Independently the negative Taylor slack (SR2) is at least
\((a/\pi^2)|p-q|^2(1-|p-q|^2/2)\).
Use the common factor \(ae^{-a}\) for these two lower bounds.
Including the angular pair density, it is bounded below by
\((\gamma/d)f(t)f(s)\), since \(2\sqrt{ts}\le t+s\).
For \(p,q\) independent under \(\mu_0\) let
\(V=\operatorname{Var}_{\mu_0}(p)\) and \(Q=\int |p|^2d\mu_0\).
Pointwise \(|p|^2\le1/2+\sum_i(p_i-1/2)_+\):
if the largest coordinate is at most \(1/2\), use
\(|p|^2\le\max_i p_i\), and otherwise the right side equals
that maximum. Hence the averaged sum of products of the
positive parts is at least \((Q-1/2)_+^2/3\).
By (SR7) and \(V\le2/3\), the two losses together give
\[
 \pi^{-2}\left[V(1-Q)+\frac53(Q-1/2)_+^2\right]
 \ge\frac{2}{5\pi^2}V,
\]
where \(1-Q+\frac52(Q-1/2)_+^2\ge2/5\).
Multiplying by \((\gamma/d)Z_0^2\) gives a preliminary bound.
These are slacks in the full aligned-arc upper bound. Applying
(DG27)–(DG29) to the exact angular pairing also subtracts the
actual first-mode deficit, with no charge to either Taylor slack.

The coefficient in this estimate can be improved to \(1/(2\pi^2)\).
Set \(u_i=(p_i-1/2)_+\), \(v_i=(q_i-1/2)_+\),
\(u=\sum_i u_i,\ v=\sum_i v_i\), and \(D_*=|p-q|^2\).
One has \(|p|^2\le1/2+2u^2\), by concentrating the mass remaining
after the dominant coordinate. If \(D_*>3/2\), both profiles have
a dominant coordinate, and their dominant labels differ. A profile
without a dominant coordinate has square norm at most \(1/2\).
Two profiles with the same dominant label are in a simplex of
diameter squared \(1/2\). Since \(u^2+v^2-uv\le1/4\) on
\([0,1/2]^2\), the case \(D_*>3/2\) obeys \(D_*-3/2\le2uv\).
Thus
\[
 D_*(D_*-3/2)_+\le4\sum_{i\ne j}u_i v_j .
\]
If \(m_i=\int u_i\,d\mu_0\), independence gives
\(\mathbb E D_*(D_*-3/2)\le 8\sum_i m_i^2\).
With \(\mathbb E D_*=2V\), it follows that
\[
 \mathbb E[D_*(1-D_*/2)]+5\sum_i m_i^2
 \ge \tfrac12 V+\sum_i m_i^2 .
\]
This proves the strengthened \(\beta_0\) in (SV1).
For the two extracted penalties with the common product weight, the
factor \(1/2\) cannot be increased. Give mass \(1-\varepsilon\) to
\((0,1/2,1/2)\) and mass \(\varepsilon\) to \((1,0,0)\).
The quartic contribution is \(V/2\) and \(\sum_i m_i^2=\varepsilon^2/4\).
This sharpness concerns only these two lower bounds on the full
angular slack.

## A stronger penalty on the existing first radial feature

The extra exponential factor in the positive same-label bound
can be removed. For \(s_i=p_i+q_i\ge1\), write
\(s=(2-x,r,x-r)\) after permuting labels, with \(0\le r\le x\le1\).
Let \(B_a(s)\) denote the positive Taylor slack in the sharp
aligned-arc pairing. Subtract its value at \(s'=(2-x,x,0)\).
Twice the difference is
\[
 G_a(x)-G_a(r)-G_a(x-r)-P_\sharp(a)r(x-r)
 =r(x-r)\big(\mathbb E k_a(U+V)-P_\sharp(a)\big),
\]
where \(U,V\) are independent uniform variables of interval lengths
\(r,x-r\), each starting at zero. The formula extends to zero
lengths by continuity. Stretching the two intervals by \(1/x\)
can only decrease the expectation. When their lengths sum to one,
the symmetric convolution whose smaller interval has length
\(b\le1/2\) has density \(v/[b(1-b)]\) for \(0\le v\le b\)
and \(1/(1-b)\) for \(b\le v\le1/2\).
Compared to the equal-length density \(4v\) on this half-interval,
the difference changes sign once from positive to negative and
has integral zero. The symmetrized function
\(\cosh(a\cos(\pi v))\) decreases there, so the smallest expected
\(k_a\) comes from equal lengths. On the full interval the
equal-length density \(4\min(v,1-v)\), minus the density in
(SP4), changes sign from positive to negative at \(4/7\) and
has integral zero. Its expectation of the decreasing \(k_a\)
is therefore at least \(P_\sharp(a)\). This proves \(B_a(s)\ge B_a(s')\).

Write \(s'=(1+y,1-y,0)\) with \(y=1-x\). The formula for \(G_a\)
and equality at \(y=0\) give
\[
 B_a(s')=\frac{y^2}{2}
   \left[P_\sharp(a)-2\int_0^1(1-v)e^{-a\cos(\pi yv)}\,dv\right]
 \ge \frac{y^2}{2}(P_\sharp(a)-N(a)),
 \quad N(a)=2\int_0^1 v k_a(v)\,dv .                            \tag{SV3}
\]
The integrand increases with \(y\in[0,1]\), and reflection
at \(y=1\) gives \(N\).

Here an elementary estimate valid for all \(a\ge0\) is
\[
 P_\sharp(a)-N(a)\ge \frac{39}{200}a
                         \ge \frac{48}{25\pi^2}a .             \tag{SV4}
\]
For details let \(b(v)\) be the signed density of this difference.
It equals \(v\) on \([0,2/3]\) and \(6-8v\) on \([2/3,1]\).
It has zero integral and changes sign from positive to negative
at \(3/4\). Put \(m_j=\int b(v)\cos^j(\pi v)\,dv\).
Folding at \(1/2\) shows \(m_{2j}\le0\) for \(j\ge1\):
the sum density is \(9v-2\) up to \(1/3\) and \(1\) afterwards,
with zero integral. Also \(m_{2j+1}+m_{2j+2}\ge0\) for \(j\ge0\).
After factoring \(\cos^{2j+1}(\pi v)\) the remaining density is
\(2-7v+(9v-2)\cos(\pi v)\) up to \(1/3\) and
\(2v-1+\cos(\pi v)\) afterwards. The former is nonnegative:
write it as \(2(1-c)+v(9c-7)\), and if \(9c<7\) use
\(v\le1/3,\ c=\cos(\pi v)\ge1/2\). The latter follows from
concavity of cosine on \([0,\pi/2]\).

For \(0\le a\le5\), every odd-even pair after degree eight
in the exponential expansion is thus nonnegative. Lower bounds
for \(10^5 m_j\), \(1\le j\le8\), are
\[
 (25330,-17098,23500,-18167,21523,-17633,19881,-16848).
\]
These follow by binomial cosine expansion from
\(\int b(v)\cos(l\pi v)dv=
[9\cos(2l\pi/3)-1-8(-1)^l]/(l\pi)^2\) for \(l>0\) and zero
for \(l=0\), using \(3.1415926<\pi<3.1415927\).
Call the eight listed integers \(m_j^-\).
The lower-bound polynomial for \((P_\sharp-N)/a\) is
\(\sum_{j=1}^8 m_j^- a^{j-1}/(10^5 j!)\).
The floors of \(10^5\) times its smallest Bernstein coefficient
on each successive half-unit interval of \([0,5]\) are
\[
 (21950,20099,19521,19524,20128,21931,25033,29548,35504,42692).
\]
For direct rational verification, expand the polynomial as
\(\sum_{j=0}^7 c_j z^j\) after substituting \(a=(l+z)/2\),
\(l=0,\ldots,9\). Its Bernstein coefficients are
\(\sum_{j=0}^i c_j\binom{i}{j}/\binom{7}{j}\) for \(0\le i\le7\).
Their lower bounds prove (SV4) up to five. For \(a\ge5\),
subtract the threshold value \(e^{-a/\sqrt2}\) from the kernel
in the signed integral. The integrand is then everywhere
nonnegative. The contribution on \([0,1/4]\) is at least
\((e^{7a/10}-e^{-7a/10})/32\), whose ratio to \(a\) increases.
At \(a=5\) the ratio exceeds \(39/200\), as follows by summing
the exponential series for \(e^{7/2}\) up to degree nine.
The last inequality in (SV4) uses the stated rational bound for
\(\pi\).

By (SV3)–(SV4) the positive Taylor slack is at least
\((96/(25\pi^2))a\sum_i u_i(p)u_i(q)\), with
\(u_i(p)=(p_i-1/2)_+\). Use \(\mu_1\) for both draws and
write \(S=\sum_i(\int u_i\,d\mu_1)^2\) and
\(V=\operatorname{Var}_{\mu_1}(p)\).
The earlier quartic argument gives the normalized quartic
term at least \(V/2-4S\), and it is also nonnegative.
The sum of the positive and quartic losses, after removing
the factor \(a/\pi^2\), is at least
\[
 \max\{V/2-(4-96/25)S,\ (96/25)S\}\ge(12/25)V .
\]
For \(S\le V/8\) use the first bound; otherwise use the second.
Multiplying by \((\gamma/d)Z_1^2=\rho\pi d/2\) gives the
second penalty in (SV1). Both are full aligned-arc Taylor
slacks, so the first angular deficit remains available.

## The corrected scalar operator

A penalty \(\beta\operatorname{Var}_{\mu}(p)\), with
\(\mu=f_\mu\nu/Z_\mu\), replaces the row and kernel in (GH3) by
\[
 J^\#=J+\frac{\beta f_\mu}{2Z_\mu},\qquad
 K^\#=K+\frac{\beta}{Z_\mu^2} f_\mu\otimes f_\mu.             \tag{SV2}
\]
Indeed the additional \(K\) term restores the square of the
\(\mu\)-mean after subtracting the second moment from \(2J\).
For \(\mu_1\) only \(b_1\) in (GH5) and \(k_1\) change,
by \(2d\beta/Z_1\) and \(\beta/Z_1^2\), respectively.
For \(\mu_0\) add the one feature \(f\) to the centered Gram
matrix. It has only a different Gaussian exponent in the
radial integrals. In either case the infinite tail bound relative
to the old full row weight \(J\) also applies to \(J^\#\).
A directed norm bound below two gives a uniform radial
margin and keeps the full actual angular deficit. The argument
of [FA12–FA13](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md#fa12-scalar-extension-and-transfer-cost-margin-through-0535)
then converts the two reserves to a planar margin relative to
transfer cost and a positive continuation past the certified endpoint.

## Directed comparison through correlation 0.623

For all \(d\ge2\), balanced ordinary Gaussian three-cell partitions
satisfy the simplex inequality through \(623/1000\).
On the added interval \([29/50,623/1000]\), ties are cylindrical
regular tripods. Use the undamped \(\mu_1\) penalty in (SV1).
The preceding [sharp positive arc comparison](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)
supplies the lower correlations.

For the infinite feature tail, the even folded signed density
in (SP4) for \(M_{2j}-P_{2j}\) has positive mass \(55/512\)
on \([0,1/16]\), where \(\cos^2(\pi v)>24/25\).
The negative part starts at \(2/9\), where the squared cosine
is below \(3/5\). The remaining positive part has squared
cosine at least that of the negative part. As in (GH9),
\[
 J(t)\ge\frac{55}{2048}e^{-\kappa x^2}
        \left(e^{24\kappa x^2/25}-e^{3\kappa x^2/5}\right),
 \qquad x=\sqrt t.
\]
For \(x\ge2\) and \(\kappa\ge1/2\), the factor
\(1-e^{-9\kappa x^2/25}\) is at least \(1/2\).
Hence \(J^\#(t)\ge J(t)\ge(55/4096)e^{-\kappa x^2/25}\)
on this range. For the trace use \(J\ge J_1\) on \(x\le2\)
as in (GH10). On \(x\ge2\), with \(\eta=1+49\kappa/25\),
integrating each Gaussian monomial bounds the outer trace by
\[
 \frac{4096(\gamma/\eta)^{N+1}}
 {55d\eta(1-\gamma/\eta)}.                                    \tag{SV5}
\]
Both geometric ratios in the inner and outer estimates are
verified below one on every certified correlation slab.

The checker encloses the signed centered feature Gram matrix.
Use the positive row polynomial truncated at \(m\) and modified
by (SV2). For its radial moments use the linear interpolation and
second derivative bound (SP7). Integrating each interpolation
basis function against each monomial gives explicit rational
weights, enclosed once in Arb. Beyond the radial cutoff, any
positive term of the row polynomial bounds the denominator
from below, giving the same Gaussian integral estimate as
(GH8) with the power shifted by that term's degree.

Let \(G\) be the exact finite Gram matrix, \(q\) the Rayleigh
quotient of a normalized trial \(u\), and \(h\) a directed upper
bound on \(\sqrt{\|G\|_{\mathrm{HS}}^2-q^2}\).
Since \(G\) is positive and \(q\le\lambda_1(G)\), \(h\) bounds
\(\lambda_2(G)\). With \(q>h\), the top eigenvalue is at most
\[
 q+\frac{\|(G-q)u\|^2}{q-h}.
\]
Indeed the expectation of \((\lambda_1-X)(X-h)\) for the
spectral measure of \(u\) is nonpositive, giving this bound by
expanding it. Floating-point diagonalization picks the trial
only. Its dyadic coordinates, quotient, residual and norm
bounds are then used with directed arithmetic on the entire slab.

The [retained checker](check_three_cell_sharp_undamped_extension.py)
and [receipt](three_cell_sharp_undamped_through_622_verification.json)
cover \([29/50,311/500]\) with 282 directed correlation slabs
obtained by dyadic subdivision from width \(1/1000\).
They use 110-bit Arb, \(m=72,\ N=144\), and 4096 radial
panels on \([0,12]\) with analytic tails. The combined finite
operator norm and infinite trace bound is below \(1999/1000\).
The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_undamped_extension.py --start 29/50 --stop 311/500 --width 1/1000 --modes 144 --rows 72 --panels 4096 --cap 1999/1000

A second [receipt](three_cell_sharp_undamped_through_623_verification.json)
covers the adjoining slab \([311/500,623/1000]\) in 96 directed
intervals. Use \(m=72,\ N=168\), 8192 radial panels and dyadic
subdivision starting at width \(1/8000\). The resulting full
norm bound is below \(19999/10000\).

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_undamped_extension.py --start 311/500 --stop 623/1000 --width 1/8000 --modes 168 --rows 72 --panels 8192 --cap 19999/10000

The norm bound leaves a radial deficit of at least
\(\frac1{10000}\int J^\#|p-\mathbf1/3|^2d\nu\). It also leaves
the entire first angular deficit from (DG29).
The first row term and the angular kernel dominate the two
reference-slice weights used in (FA27), since
\(\kappa\le2/3\) and \(\gamma\sqrt2\le3\) on this interval.
Thus the argument of FA11–FA13 gives a uniform planar
transfer-cost margin here. At the upper endpoint (CT5) extends
the comparison to a positive unevaluated right neighborhood.
The existing maximizer dimension theorem passes both optimality
and rigidity to every \(d\ge2\).

## Increasing the variance gain by feasibility

The undamped variance gain in (SV1) uses no purity assumption.
A larger gain is available on profiles that could have
nonpositive deficit. Indeed (SR7) alone yields
\[
 \lambda(P)\operatorname{Var}_{\mu_1}(p),
 \qquad \lambda(P)=\frac{\rho d}{2\pi}P,\qquad
 P=\int(1-|p|^2)\,d\mu_1 .
 \tag{SV6}
\]
Use this loss in place of the undamped gain once a lower bound
on \(P\) has been obtained. The two penalties are never added.

Here is the directed feasibility test for a gain \(\lambda(P_0)\)
already known to apply to a class of profiles. In the Hilbert
space \(H=L^2(J'_m\,d\nu)\), \(J'_m\) is the modified truncated
row weight. Let \(F\) be synthesis of the centered features
\(\sqrt{k_i}(f_i-\overline f_i)/J'_m\) with
\(\overline f_i=U_i/W\) as in (GH6), and \(G=F^*F\).
For a nonzero trial \(v\) in the finite feature coordinates
set \(u=Fv/\sqrt{v^{\mathsf T}Gv}\). Its exact Rayleigh quotient is
\(a_u=|Gv|^2/(v^{\mathsf T}Gv)\). For the orthogonal projection
\(P_\perp\) off \(u\) one has
\[
 \|P_\perp FF^*u\|^2
   \le\frac{\|G\|_{\rm HS}}{v^{\mathsf T}Gv}|Gv-a_u v|^2,\qquad
 \|P_\perp FF^*P_\perp\|^2
   \le\|G\|_{\rm HS}^2-a_u^2 .
\]
The first bound follows by writing the residual as
\(F(Gv-a_u v)/\sqrt{v^{\mathsf T}Gv}\). The second follows
by expanding the squared Hilbert–Schmidt norm in the two
orthogonal blocks.

Write \(z=\alpha u+g\) with \(g\perp u\) componentwise,
\(A=|\alpha|\), \(G_0=\|g\|_H\). Include the old-row trace
bound as in (SV5). The resulting scalar upper bound is
\[
 q_{\rm pre}(z)\le a_0 A^2+2r_0 A G_0-\delta G_0^2 ,
 \quad a_0,r_0\ge0,\quad\delta>0 .                           \tag{SV7}
\]
In fact a tail bound \(\tau\) against the full modified row
gives negative coefficient \(2-\tau\) against \(J'_m\).
Subtract the trial Rayleigh quotient and the complementary
norm respectively from this coefficient, and bound the
cross term by the projected residual.

Fix a reserve \(\epsilon\ge0\). If
\(q_{\rm pre}\ge-\epsilon(A^2+G_0^2)\) and \(A>0\), then
\(G_0/A\le H_0\), where \(H_0\) is the positive root of
\((\delta-\epsilon)h^2-2r_0h-(a_0+\epsilon)=0\).
If \(A=0\), the same antecedent forces \(G_0=0\), provided
\(\delta>\epsilon\). For the new desired purity \(P_1\) use
\[
 A_0\le\sqrt{\frac{2/3-P_1}{B_0(1+H_0^2)}},\qquad
 B_0=\frac{4d}{Z_1 b'_1},
\]
where \(b'_1\) is the first row coefficient after the gain.
Indeed \(d\mu_1/(J'_m\,d\nu)\le B_0\) by that first term.
If \(A\ge A_0\), the mass simplex bound
\(|p-\mathbf1/3|\le\sqrt{2/3}\) gives, with
\(Q(x)=\sum_{j=1}^m b'_j x^{j-1}\),
\[
 \frac{G_0^2}{A^2}
 \ge\int_0^1\frac{Q(x)e^{-(1+\kappa)x^2}}{2d}
        \bigl(|x u(x^2)|-\sqrt{2/3}\,x/A_0\bigr)_+^2\,dx . \tag{SV8}
\]
The substitution uses
\(J'_m(x^2)d\nu(x^2)=Q(x)x^2 e^{-(1+\kappa)x^2}dx/(2d)\).
It removes the apparent origin pole from the numerical
integral. Enclose the right side from below on rational
subintervals. When it exceeds \(H_0^2\), admissible
profiles have \(A<A_0\), whence
\[
 \int(1-|p|^2)d\mu_1
 \ge \frac23-B_0(A^2+G_0^2)\ge P_1.
\]
Apply (SV6) at the new purity bound. This step can be repeated.
At each stage a profile with
\(q_{\rm pre}<-\epsilon\int J'_m |z|^2d\nu\) already has a
radial margin, while the remaining profiles pass to the next
stronger scalar bound. Every such row dominates the original
first radial term. Thus a positive reserve at the feasibility
tests and a strict final norm bound give a uniform radial
margin, with the entire first angular deficit available.

### A stronger gain from the balanced tilted mean

The quartic loss can exceed the purity coefficient in (SV6).
Write \(z=p-\mathbf1/3\), \(M=\int|z|^2d\mu_1\),
\(m=|\int z\,d\mu_1|\), and \(V=M-m^2\).
The fourth moment expansion (SR7), together with
\(|p-\int p\,d\mu_1|\le\sqrt{2/3}+m\) and
\(\operatorname{tr}(\Sigma^2)\le V^2\), gives the quartic
bound \([2-(\sqrt{2/3}+m)^2-3V]V\).
For balanced profiles and \(29/50\le\rho\le2/3\), one has
\[
 m^2\le \chi V,\qquad
 \chi=\frac{\pi(1-\rho^2)^2}{2\sqrt{1-2\rho^2}}-1\le\frac12.
\]
Indeed \(d\nu/d\mu_1-1\) has squared \(L^2(\mu_1)\)-norm
\(\chi\), by direct integration, and its pairing with the
\(\mu_1\)-centered profile is the negative tilted mean.
The function \(d^2/\sqrt{1-2\rho^2}\) increases with \(\rho^2\)
above \(1/3\); at \(2/3\) the value of \(\chi\) is
\(25\pi/54-1<1/2\). Thus \(m\le\sqrt{M/3}\).
If \(M\le1/6\), \(2\sqrt{2/3}\,m-2m^2\) increases throughout
this \(m\)-range, giving the lower coefficient of \(V\)
\[
 \frac{4-2\sqrt{2M}-7M}{3}.                                   \tag{SV9}
\]
To pass to gain \(\rho d P_1/(2\pi)\) for \(0<P_1<4/3\),
it suffices to bound \(M\) above by
\[
 \min\!\left\{\frac16,
     \left(\frac{\sqrt{120-84P_1}-2\sqrt2}{14}\right)^2\right\}.
\]
Substitute this cap for \(2/3-P_1\) in the expression for
\(A_0\). Gains above \(2/3\) are thus feasible targets for
the same comparison and do not require a change of feature base.

### Endpoint bounds for a fixed scalar gain

The finite Gram matrix also admits a comparison using only the
endpoints of a correlation slab. Fix a gain \(b\,\rho d/(2\pi)\)
with \(b\ge0\), and use \(v=(\rho/\sqrt d)\sqrt t\). If the
truncated row is written \(J'_m=e^{-v^2}vP_m(v)/4\), its
polynomial is independent of \(\rho\). Its \(v^{j-1}\) coefficient
before the gain is \(2^j(M_j-P_j)\Gamma(1+j/2)/j!\); the gain
adds \(b/(\pi\Gamma(3/2))\) to the constant coefficient.
The kernel feature coefficients in these variables are
\(d^{-1}\rho^{-j}a_j\), where
\(a_j=2^{j-1}(M_j+P_j)/j!\) and the gain adds \(2b/\pi^2\) to
\(a_1\).

Let \(C_0(\rho)\) be the Gram matrix of \(v^j e^{-v^2}\),
after subtracting their means in the measure
\[
 \frac{\exp[-(\rho^{-2}-2)v^2]}{P_m(v)}\,dv .
\]
That measure is finite for \(0<\rho<1/\sqrt2\). The exact
centered feature Gram matrix in (GH6), after the gain, is
\(8D(\rho)C_0(\rho)D(\rho)\) with
\(D_{jj}(\rho)=\rho^{-1-j/2}\sqrt{a_j}\).
This follows by writing \(d\nu/J'_m\) in the new radial variable.
The matrices \(C_0(\rho)\) increase in positive semidefinite
order: every one of their quadratic forms is the infimum over
a constant \(c\) of the integral of \((f-c)^2\), with the
displayed density increasing pointwise in \(\rho\).

Hence on \([\rho_-,\rho_+]\) the exact finite norm is bounded
above by the norm of
\(8D(\rho_-)C_0(\rho_+)D(\rho_-)\).
After the positive semidefinite upper bound by \(C_0(\rho_+)\),
use the equality of the nonzero spectra of
\(D C_0 D\) and \(C_0^{1/2}D^2 C_0^{1/2}\).
For entry enclosures of \(C_0(\rho)\), the difference from the
right endpoint is positive semidefinite and its \(ij\) entry
in absolute value is bounded by the geometric mean of the
\(i,i\) and \(j,j\) entries of \(C_0(\rho_+)-C_0(\rho_-)\).
The diagonal feature factors must still be included separately.
These endpoint bounds provide an alternative finite-block
certificate; they are not used in the retained interval runs.

### Directed bootstrap through 0.630

The [bootstrap checker](check_three_cell_sharp_variance_bootstrap_anchor.py)
and [retained cover](three_cell_sharp_variance_bootstrap_anchor_through_625_verification.json)
apply (SV7)–(SV8) on \([623/1000,5/8]\). The preliminary
gain starts at \(P_0=12/25\) and the directed test gives
\(P_1=3/5\), with a \(10^{-5}\) reserve in the antecedent.
There are 18 accepted correlation slabs starting from a cover
of width \(1/8000\), with \(N=144,\ m=72\), 2048 radial panels
on \([0,12]\), and 768 clipping panels on \([0,1]\).
The final norm enclosure is below \(1999/1000\).
The implementation uses a fixed exact dyadic anchor for the
spectral residual in feature coordinates. The residual around
the actual quotient has smaller norm, so the off-diagonal
bound following (SV6) can also use the anchor residual.
Likewise Temple's inequality with that larger residual gives
a bound on \(\|G\|\) and on \(a_u\), improving the
feasibility test.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_variance_bootstrap_anchor.py --start 623/1000 --stop 5/8 --width 1/8000 --modes 144 --rows 72 --panels 2048 --purity 3/5 --cap 1999/1000

The [adjoining cover](three_cell_sharp_variance_bootstrap_anchor_through_630_verification.json)
continues the same reserved test through \(63/100\), with
the same purity threshold \(3/5\). It uses 98 accepted slabs,
\(m=76,\ N=160\), and 4096 radial panels. The final norm bound
is below \(19999/10000\).

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_variance_bootstrap_anchor.py --start 5/8 --stop 63/100 --width 1/8000 --modes 160 --rows 76 --panels 4096 --purity 3/5 --cap 19999/10000

Together with the preceding slabs this proves the unrestricted
balanced comparison for \(0<\rho\le63/100\) in every \(d\ge2\),
with simplex equality rigidity on the added slabs. Since the
bootstrap keeps a positive radial reserve on both branches
and the exact first angular deficit, the same reference-slice
argument gives a positive unevaluated right neighborhood
beyond \(63/100\).

### Hermite and endpoint cover through 253/400

For \([63/100,253/400]\) the same reserved test with
\(P_1=311/500\) uses a [cubic interpolation module](sharper_hermite_radial_quadrature.py)
and evaluates the final Gram norm at the right endpoint with
the diagonal rescaling described above. The preliminary test
still encloses each whole slab.

For \(F(x)=8d\,e^{-\eta x^2}/Q(x)\) on \([l,r]\), put
\(T_j=Q^{(j)}(r)/Q(l)\). The positivity of the coefficients
bounds \(Q^{(j)}(x)/Q(x)\le T_j\) and \(F(x)\le F(l)\).
Define
\[
\begin{split}
 L_1&=2\eta r+T_1,\qquad L_2=2\eta+T_2+T_1^2,\\
 L_3&=T_3+3T_1T_2+2T_1^3,\\
 L_4&=T_4+4T_1T_3+3T_2^2+12T_1^2T_2+6T_1^4 .
\end{split}
\]
These bound the absolute derivatives of \(\log F\) through
order four. So the fourth derivative of \(F\) is bounded by
\(F(l)(L_1^4+6L_1^2 L_2+3L_2^2+4L_1L_3+L_4)\).
The cubic using values and slopes at both endpoints has absolute
error at most \((r-l)^4/384\) times this bound, by the Hermite
interpolation remainder. Integrate its four basis polynomials
and the error bound against each radial monomial. The
[directed checker](check_three_cell_bootstrap_hermite_right.py)
uses 512 radial panels, 768 clipping panels, \(m=80,\ N=168\),
and 110-bit Arb. The [full receipt](three_cell_bootstrap_hermite_right_630_6325_verification.json)
contains 106 adjacent accepted intervals and an improved full
norm bound below \(19999/10000\), with the same \(10^{-5}\)
preliminary reserve.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_bootstrap_hermite_right.py --start 63/100 --stop 253/400 --width 1/2000

The relative margin proof uses
\(w(t)=\sqrt t\exp[-3t/4-4\sqrt t]\) in place of the previous
reference-slice lower weight. The explicit radial and angular
weights dominate it because \(\kappa\le3/4\) and
\(\gamma\sqrt2\le4\) on this interval. The comparison and
equality rigidity therefore continue through \(253/400\)
and on a positive unevaluated right neighborhood.

### Clipping with the radial orthogonality constraint

There is a further lower bound in (SV8) from the orthogonality
of \(g\). Pointwise let \(y=Awu\) with \(w=\alpha/A\), let
\(z_*(t)\) be the Euclidean projection of \(y(t)\) onto
\(\Delta_3-\mathbf1/3\), and put \(g_*=z_*-y\). Projection
Pythagoras gives \(\|g\|_H^2\ge\|g_*\|_H^2+\|g-g_*\|_H^2\).
Since \(g\perp wu\) and \(\|wu\|_H=1\), the second term bounds
\(|\langle g_*,wu\rangle_H|^2\). The projection fixes zero
and is nonexpansive, so \(-g_*(t)\cdot y(t)\ge0\) and is at least
\(|y(t)|(|y(t)|-\sqrt{2/3})\) where \(|y(t)|>\sqrt{2/3}\).
Thus
\[
 \frac{G_0^2}{A^2}\ge I_2(A)+ I_1(A)^2,\qquad
 I_j(A)=\begin{cases}
 \int J'_m |u|(|u|-\sqrt{2/3}/A)_+\,d\nu,&j=1,\\
 \int J'_m (|u|-\sqrt{2/3}/A)_+^2\,d\nu,&j=2 .
 \end{cases}                                                   \tag{SV10}
\]
The full proof uses a projection only to justify the extra term;
both integrals can be lower bounded using the same one-variable
clip panels. A direction-dependent version replaces \(\sqrt{2/3}\)
by the support function of \(\Delta_3-\mathbf1/3\) at
\(\operatorname{sgn}(u(t))w\). Also \(\int g\,d\nu=0\) and
\(\int u\,d\nu=0\), so if \(1/J'_m\in H\), one can add
\(|\int g_*\,d\nu|^2/\int (J'_m)^{-1}d\nu\) to the original
projection inequality before majorizing its terms. These refinements
do not enter the preceding directed cover.

### Balanced second-moment bootstrap through 0.641

The [two-stage checker](check_three_cell_sharp_variance_bootstrap_SR14_hermite.py)
uses the sharper tilted-mean expression in (SV9) before replacing
\(\chi\) by \(1/2\). On \([63/100,641/1000]\), it starts with
coefficient \(12/25\), clips with \(M\le47/300\), and verifies
coefficient \(31/50\) using
\(m\le\sqrt{\chi M/(1+\chi)}\) in the quartic bound.
With that coefficient the next clip gives \(M\le7/150\),
which verifies coefficient \(1\). Each stage uses its own row
weight and the preliminary reserve \(10^{-5}\).
The scalar bound with the last coefficient is below \(1999/1000\).

The [full cover](three_cell_sharp_variance_bootstrap_hermite_through_641_verification.json)
has 349 accepted rational slabs, from initial width \(1/1000\).
The checker uses 110-bit Arb, \(N=192,\ m=96\), 512 cubic
interpolation panels and 768 clipping panels. It checks the
tilted mean on each slab with the exact formula for \(\chi\).
Both preliminary tests and the final Gram norm enclose each
entire slab, with the analytic tails as in (SV5).

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_variance_bootstrap_SR14_hermite.py --start 63/100 --stop 641/1000 --width 1/1000 --modes 192 --rows 96 --panels 512 --cap 1999/1000

Thus the unrestricted balanced Gaussian comparison holds in
every \(d\ge2\) for \(0<\rho\le641/1000\), with the regular
tripod as the only equality case on the added slabs up to
rotation, labels and null sets. The weight
\(\sqrt t\exp[-3t/4-4\sqrt t]\) still works through \(641/1000\).
The radial reserve and exact angular deficit give a uniform
planar transfer-cost margin on the new interval and a positive,
unevaluated right neighborhood beyond \(641/1000\).

### Projection of the tilted mean

For further use in the bootstrap, let
\(w=P_{\rm bal}(f_1/(Z_1 J'_m))\), with \(P_{\rm bal}\) the
weighted orthogonal projection onto the ordinary-mean-zero
subspace. The tilted mean of \(z=\alpha u+g\) then satisfies
\[
 \left|\int z\,d\mu_1\right|
 \le A|\langle u,w\rangle_H|
     +G_0\sqrt{\|w\|_H^2-|\langle u,w\rangle_H|^2}.           \tag{SV11}
\]
This uses \(g\perp u\) and the mass constraint.
The inputs are already in the centered feature calculation:
\(\|w\|_H^2=(T_2-U_1^2/W)/Z_1^2\) and
\(\langle u,w\rangle_H=(Gv)_1/(Z_1\sqrt{k_1}\sqrt{v^{\mathsf T}Gv})\),
with feature index 1 denoting the first nonconstant monomial
degree. Taking the smaller of (SV11) and the bound from \(\chi\)
can improve the quartic coefficient in (SV9). The exact
complementary feature Gram matrix is
\(G-(Gv)(Gv)^{\mathsf T}/(v^{\mathsf T}Gv)\), whose nonzero
spectrum is that of \(P_\perp FF^*P_\perp\). It can replace the
Hilbert–Schmidt bound for the complementary norm in (SV7).

### Projected-mean bootstrap through 129/200

The [projected-mean checker](check_three_cell_sharp_variance_mean_SR19.py)
uses (SV10) in the clipping test and (SV11) for the tilted mean.
It also uses the exact two-dimensional multiplier bound (SR19).
To upper bound \(M_u=\int u^2d\mu_1\), let
\(B_0\ge\sup d\mu_1/(J'_m d\nu)\) as above. Since
\(\|u\|_H=1\), one may subtract
\[
 \int_{0\le\sqrt t\le2}
   \left(B_0-\frac{d\mu_1}{J'_m d\nu}\right)u(t)^2 J'_m(t)\,d\nu(t)
\]
from \(B_0\). The integrand is nonnegative and is bounded below
by directed Riemann sums on 1024 panels. If the resulting upper
bound for \(M_u\) is at most \(B_0/(1+H_0^2)\), (SR19)
improves the bound for \(M/A^2\) from \(B_0(1+H_0^2)\)
to \((\sqrt{M_u}+H_0\sqrt{B_0-M_u})^2\), using the upper
value of \(M_u\). Otherwise use the original coefficient.

On \([641/1000,129/200]\) the stages use \(M\le1/6\) and
\(M\le7/200\) to increase the variance coefficient from \(12/25\)
to \(13/20\), then \(11/10\). Both stages keep the reserve
\(10^{-5}\). The [directed cover](three_cell_sharp_variance_mean_SR19_through_645_verification.json)
uses 180 accepted rational slabs from width \(1/4000\),
\(N=204,\ m=104\), 512 cubic radial panels, 768 clipping
panels and 110-bit Arb. The final norm is below \(1999/1000\).

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_variance_mean_SR19.py

Together these give the unrestricted balanced Gaussian comparison
through \(129/200\), with simplex rigidity on the added intervals.
The same enlarged reference-slice weight and the retained angular
deficit give a uniform score margin, and (CT5) gives a positive
unevaluated right neighborhood of \(129/200\).

### Adaptive projected-mean bootstrap through 81/125

The [adaptive checker](full_angular_directmean_bootstrap.py)
uses the exact complementary Gram matrix following (SV11),
which permits a tighter complementary norm bound in (SV7).
It chooses the amplitude cap by rational bisection using (SV8).
For the input \(M_u\) in (SR19), it directly integrates the
tilted trial square. Writing \(P(x)=\sum_j\sqrt{k_j}v_j x^j\)
and \(c_0=\sum_j\sqrt{k_j}v_j U_j/W\), this square is
\[
 M_u=\frac{32d^2}{Z_1(v^{\mathsf T}Gv)}
   \int_0^\infty
    \frac{\big(e^{-\kappa x^2}P(x)-c_0\big)^2
                       e^{-(1-\kappa)x^2}}{Q(x)^2}\,dx .
\]
The checker encloses the integral with 4096 panels on
\([0,12]\). Beyond 12 it uses
\((a-b)^2\le2a^2+2b^2\) and the lower bound on \(Q\) by one
positive monomial. Expanding \(P^2\) in absolute monomials
then gives explicit Gaussian tails as in (GH8).

Combining this upper bound with (SR19) gives a bound on the
raw second moment of any profile passing the clip. The
checker bounds its mean by the smaller of (SV11) and the
\(\chi\)-bound, and rounds down the resulting directed lower
bound on the quartic coefficient to a rational multiple of
\(1/1000\). It verifies that strict lower bound and that the
coefficient increases before starting the next stage. The
slab is accepted only when the final directed norm is below
\(1999/1000\). The stage estimates retain \(10^{-5}\) of
their weighted squared radial norm in the antecedent.

The [complete cover](three_cell_directmean_bootstrap_645_648_verification.json)
of \([129/200,81/125]\) consists of 172 accepted rational
slabs using two to four gain stages each. It uses
110-bit Arb, 512 cubic radial panels, \(N=192,\ m=96\),
and 768 clip panels. The three adjoining internal replays
reproduce the exact list of accepted intervals. The maximum
directed upper norm bound is \(1.9988760641266756\).

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/full_angular_directmean_bootstrap.py --start 129/200 --stop 81/125 --width 1/4000

Thus the balanced Gaussian comparison holds in every
\(d\ge2\) for \(0<\rho\le81/125\), with simplex equality
rigidity on the added intervals. The \(10^{-5}\) radial reserve
and the full actual first angular deficit give the same
uniform score margin. The reference-slice weight
\(\sqrt t\exp[-3t/4-4\sqrt t]\) still works on the added slab.
Applying (CT5) at \(81/125\) gives a positive unevaluated
right neighborhood.

## Discrete plurality

For any fixed positive correlation in the certified global interval,
balanced three-candidate Plurality Is Stablest follows with the usual
low-degree influence hypothesis. Given \(\varepsilon>0\), there
are an integer \(D\) and \(\tau>0\) such that every
\(f:[3]^n\to\Delta_3\) with uniform mean \(\mathbf1/3\) and
\(\operatorname{Inf}_i^{\le D}(f_j)\le\tau\) for all \(i,j\)
has noise stability at most the infinite-voter plurality limit
plus \(\varepsilon\).
The lower-degree influences sum the Fourier weights of degree at
most \(D\) whose supports contain the specified coordinate.
[Isaksson–Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
give a balanced ordinary Gaussian partition with stability
within \(\varepsilon\) of the discrete rule and identify the
plurality limit with the Gaussian simplex stability.
