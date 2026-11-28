# Circular purity at high three-cell stability

Let \(A\) be an ordinary balanced planar three-cell partition and \(T\)
the regular simplex. For \(63/100\le\rho\le68/100\), the condition
\(S_\rho(A)\ge S_\rho(T)\) forces the mean circular mixing mass to
exceed \(119/500\) under the Gaussian radial law and \(11/80\)
under the first tilted radial law. The conclusion applies to every
balanced planar partition that reaches the simplex value.
In the angular Taylor comparison there is also an independent
positive remainder when two circular profiles favor different
labels and their combined mass in the third label is appreciable.

Use \(t=r^2/2\), \(d\nu(t)=e^{-t}dt\), and write \(p_i(t)\) for
the circular masses. Put \(p_*=(1/3,1/3,1/3)\),
\(L(t)=1-|p(t)|^2\), \(P=\int L\,d\nu\), and, for the tilted law
\(\mu_1\) below, \(P_1=\int L\,d\mu_1\). Thus \(L\in[0,2/3]\) and
larger \(L\) means more mixing. Let \(W_m\) denote the total Gaussian
Hermite weight at degree \(m\).

## Hermite weights in terms of circular mixing

For \(j\ge1\) define the angular moments
\[
 u_i^{(j)}(t)=\int_{A_i(t)}(\cos(j\theta),\sin(j\theta))
                                      \frac{d\theta}{2\pi}.
\]
For any projection direction, the corresponding cosine has the
same value distribution as \(\cos\theta\). Taking the largest values
at fixed mass bounds \( |u_i^{(j)}(t)|\le\sin(\pi p_i(t))/\pi\).
The sine cap in (SR17) of the [sharp scalar remainders](THREE_CELL_SHARP_SCALAR_REMAINDERS.md#rank-one-bounds-for-the-quartic-penalty)
gives \(\sin(\pi p_i)\le4p_i(1-p_i)\). Since
\(\sum_i p_i(1-p_i)=L\) and \(p_i(1-p_i)\le L/2\),
\[
 \sum_i|u_i^{(j)}(t)|^2\le \frac8{\pi^2}L(t)^2.              \tag{PH1}
\]
The scalar cap follows by applying the chord bound
to the convex function \(x\mapsto\cos(\pi\sqrt x/2)\) on
\([0,1]\) at \(x=(2p_i-1)^2\). The bound on the individual
purity term follows from \(L=2\sum_{i<j}p_i p_j\).

The first Gaussian moment uses only \(j=1\), and (PH1) gives
\[
 W_1\le \frac8{\pi^2}\left(\int rL\,d\nu\right)^2.          \tag{PH2}
\]
The degree-two angular part \(W_{2,\mathrm{ang}}\) uses the two
polynomials \(r^2\cos(2\theta)\), \(r^2\sin(2\theta)\), each
of squared Gaussian norm \(4\). Similarly the angular degree-three
space is spanned by the cosine and sine pairs with radial parts
\(r(r^2-4)\) and \(r^3\), of angular orders \(1\) and \(3\)
and squared norms \(8\) and \(24\), respectively. These four
polynomials are mutually orthogonal and orthogonal to all lower
degrees. For example \(r(r^2-4)\cos\theta=x(x^2+y^2-4)\) is the
sum of \(x^3-3x\) and \(x(y^2-1)\), of squared norms \(6\) and
\(2\). The radial orthogonality to degree one also follows from
\(\int r^4d\nu=4\int r^2d\nu\). The other pair are harmonic
homogeneous cubics. Applying (PH1) to the projections gives
\[
 \begin{split}
 W_{2,\mathrm{ang}} &\le \frac2{\pi^2}
                  \left(\int r^2L\,d\nu\right)^2,\\
 W_3 &\le \frac1{\pi^2}\left(\int r|r^2-4|L\,d\nu\right)^2
           +\frac1{3\pi^2}\left(\int r^3L\,d\nu\right)^2 .    \tag{PH3}
 \end{split}
\]
In these integrals \(r=\sqrt{2t}\). Each right side also permits
a one-dimensional upper bound at fixed \(P\), by putting the
mass \(L\le2/3\) on the largest values of its nonnegative radial
weight. The norm bounds apply to arbitrary angular shapes.

## A necessary mixing floor for high stability

Separate each cell indicator into its angular mean and its
angular-mean-zero part. These spaces are Gaussian orthogonal and
invariant under the Ornstein–Uhlenbeck operator. The nonconstant
radial weight is \(2/3-P\), and the total angular weight is \(P\).
The degree-two radial weight is at most
\[
 M_r=\frac{4a^2+(2a-b)^2+b^2}{9}<.217,
 \qquad a=\log(3/2),\quad b=\log 3.                           \tag{PH4}
\]
Indeed it is the squared norm of the three integrals of \(t-1\)
against \(p_i\). To maximize that norm, fix a vector of integrals
and sort its label coefficients. A linear functional of the
integrals then maximizes by filling the three ordered quantile
intervals of the exponential measure, with endpoints \(a,b\).
This uses the upper-tail integral of \(t-1\) for each nonnegative
increment in the sorted coefficients. Cauchy–Schwarz bounds the
original norm by the norm of this quantile assignment. Integrating
over the intervals gives (PH4), including for fractional profiles.
The remaining radial degrees are at least four. Applying the
bound \(\rho^2\) to all angular degrees other than one yields
\[
 S_\rho(A)\le\frac13+\frac23\rho^4
    +(\rho^2-\rho^4)(M_r+P)+(\rho-\rho^2)W_1.                \tag{PH5}
\]
This is a convenient version of (JA9) of the
[joint angular analysis](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md#sharp-circular-mass-penalty-for-the-first-moment).
It applies to every balanced planar partition.

One may also keep all three angular degrees bounded in (PH2)–(PH3).
Using \(\rho^4\) for the remaining angular part gives
\[
 S_\rho(A)\le\frac13+\frac23\rho^4+(\rho^2-\rho^4)
   \min(M_r,2/3-P)+(\rho-\rho^4)W_1
   +(\rho^2-\rho^4)W_{2,\mathrm{ang}}+(\rho^3-\rho^4)W_3 .
\]
The three angular weights in this bound sum to at most \(P\).
At fixed \(P\) and known individual caps, maximizing this linear
expression fills the three weights in decreasing order of their
coefficients, up to the total allowed \(P\).

Suppose \(P\le119/500\). The largest possible value of
\(\int rL\,d\nu\) allows \(L=2/3\) on an outer tail of mass
\(357/1000\). Let the threshold in \(t\) be
\(t_0=\log(1000/357)>1.030\). For the cumulative gamma law \(F\)
of shape \(3/2\) and unit scale, the integral of
\(\sqrt t e^{-t}\) over this tail is
\(\Gamma(3/2)[1-F(t_0)]\). Thus (PH2) gives
\[
 W_1\le \frac{16}{9\pi}[1-F(t_0)]^2<.178.                    \tag{PH6}
\]
The last inequality uses \(F(1.030)>.4399\), verified by the
gamma series below. All decimals in the bounds here denote
rationals. The logarithm bounds follow from
\[
 \log x=2\sum_{n=0}^m\frac{y^{2n+1}}{2n+1}+\mathcal R_m,
 \quad y=(x-1)/(x+1),\quad
 0<\mathcal R_m\le\frac{2 y^{2m+3}}{(2m+3)(1-y^2)}
\]
for \(x>1\), using \(m=9\). This gives
\(.4054<a<.4055\), \(1.0986<b<1.0987\), and the indicated
bound on \(t_0\).

Substitute (PH4) and (PH6) into (PH5) with \(P\le.238\).
The resulting upper bound is
\[
 C(\rho)=\frac13+\frac23\rho^4
                 +.455(\rho^2-\rho^4)+.178(\rho-\rho^2).
\]
On \([.63,.68]\), \(C'(\rho)>.73\). The
[exact simplex stability](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)
has derivative
\[
 S_T'(\rho)=\frac3{2\pi^2}\left[
       \frac{\arccos(-\rho)}{\sqrt{1-\rho^2}}
       +\frac{\arccos(\rho/2)}{2\sqrt{1-\rho^2/4}}\right]<.66
\]
there. For the upper bound use
\(\arccos(-\rho)\le\pi/2+\rho/\sqrt{1-\rho^2}<2.51\),
\(\arccos(\rho/2)\le\pi/2\), and the two square roots greater
than \(.73,.94\). At the right endpoint \(C(.68)<.627716\) and
\(S_T(.68)>.6279\). For the latter use
\(\arcsin(.68)>.74769\) and \(.3469<\arcsin(.34)<.347\) in
\[
 S_T(\rho)=\frac13+\frac3{4\pi^2}
 (\arcsin\rho+\arcsin(\rho/2))(\pi+\arcsin\rho-\arcsin(\rho/2)).
\]
These bounds follow from the binomial series for the derivative
of arcsine, integrating through degree thirteen. The omitted
integrated terms at \(x=.34\) are bounded by
\(x^{15}/[15(1-x^2)]\). The rational enclosures
\(3.14159<\pi<3.14160\) suffice throughout. Since \(S_T-C\)
decreases on this interval and is positive at its right endpoint,
\[
 S_\rho(A)\ge S_T(\rho),\quad .63\le\rho\le.68
               \quad\Longrightarrow\quad P>119/500 .       \tag{PH7}
\]

## Transferring the floor to the first tilted law

Put \(d=1-\rho^2\), \(\kappa=\rho^2/d\), and
\[
 w_\rho(t)=\frac{\sqrt t\,e^{-\kappa t}}{\Gamma(3/2)d^{3/2}},
 \qquad d\mu_1=w_\rho\,d\nu.
\]
For fixed \(t>0\), \(\log w_\rho(t)\) is concave in \(\kappa\),
since its nonaffine term is \(3\log(1+\kappa)/2\). Thus
\(w_\rho(t)\ge\min(w_l(t),w_u(t))\) for \(l=.63,u=.68\).
For any \(\beta>0\) and every profile,
\[
 P_1\ge\beta P-\frac23\int(\beta-\min(w_l,w_u))_+d\nu .      \tag{PH8}
\]
The analogous bound with \(w_\rho\) in place of the minimum is
the density-threshold bound at a fixed correlation. It is optimal
among measurable \(L\in[0,2/3]\) with prescribed \(P\) when
\(\beta\) is chosen so that the low-density set has measure
\(3P/2\). This follows by filling that set with \(L=2/3\).

Take \(\beta=7/8\). The first threshold crossing of \(w_l\) is
in \((.163,.164)\), and the last threshold crossing of \(w_u\)
is in \((1.679,1.680)\). The densities are unimodal, both are
above \(\beta\) on \([.164,1.679]\), and their ratio \(w_u/w_l\)
is decreasing. At \(.164\) the ratio exceeds one, while at
\(1.679\) it is below one. It follows that the respective low
and high tails in the integral in (PH8) use \(w_l\) and \(w_u\).
Writing \(d_l=1-l^2,d_u=1-u^2\), that integral is at most
\[
 \beta(1-e^{-.164}+e^{-1.679})
                -F(.163/d_l) -[1-F(1.680/d_u)] < .105663,
\]
The following rational intervals verify the four
threshold conditions, their density order at the two central
endpoints, and the needed tail integrals.

| Quantity | Lower bound | Upper bound |
|---|---:|---:|
| \(w_l(.163)\) | .8737 | .8738 |
| \(w_l(.164)\) | .8758 | .8759 |
| \(w_u(.164)\) | 1.0067 | 1.0069 |
| \(w_l(1.679)\) | 1.0339 | 1.0341 |
| \(w_u(1.679)\) | .8751 | .8753 |
| \(w_u(1.680)\) | .8746 | .8748 |
| \(e^{-.164}\) | .8487 | .8488 |
| \(e^{-1.679}\) | .1865 | .1866 |
| \(F(.163/d_l)\) | .0900 | .0902 |
| \(F(1.030)\) | .4399 | .4401 |
| \(F(1.680/d_u)\) | .8998 | .9000 |

For direct verification, expand \(e^{-x}\) using Taylor
polynomials of degrees 19 and 20, which are lower and upper
bounds for \(x\ge0\). The gamma enclosures use their termwise
integrals, or equivalently
\[
 F(x)=\frac{2x^{3/2}}{\sqrt\pi}
                \sum_{j\ge0}\frac{(-x)^j}{j!(j+3/2)}.
\]
Square roots in the table can be enclosed to \(10^{-8}\) by
squaring rational endpoints. Combining (PH7) and (PH8) gives
\[
 S_\rho(A)\ge S_T(\rho),\quad .63\le\rho\le.68
               \quad\Longrightarrow\quad P_1>11/80.        \tag{PH9}
\]
In particular the pure-vertex radial profile does not compete on
this interval. The estimate does not require a first tie or a
derivative sign for the fixed-partition deficit.

## A positive Taylor reserve for mixed pairs

Use \(k_a(v)=e^{a\cos(\pi v)}\), and let \(B_a(s)\) denote the
positive Taylor slack in the sharp aligned-arc bound, normalized
as in (SV3) of the [variance absorption](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md#a-stronger-penalty-on-the-existing-first-radial-feature).
For \(s=(2-x,r,x-r)\) with \(0\le r\le x\) the exact
difference from the edge comparison is
\[
 B_a(s)-B_a(2-x,x,0)=\frac{r(x-r)}2
                        [\mathbb E k_a(U+V)-P_\sharp(a)],    \tag{PH10}
\]
where \(U,V\) are independent uniform variables on \([0,r]\)
and \([0,x-r]\). This is the two integrations of \(k_a\) in
(SV3). It extends to a zero-length interval by continuity.
The edge slack is nonnegative by (SP2).

For \(x\le1\) the argument of (SV3) bounds the bracket below
by \(W(a)=\mathbb E k_a(U_0+V_0)-P_\sharp(a)\), with
\(U_0,V_0\) independent uniform on \([0,1/2]\). The signed
density of \(W\) is \(v\) up to \(1/2\), \(4-7v\) from
\(1/2\) to \(2/3\), and \(-2(1-v)\) from \(2/3\) to \(1\).
Its integral against \(\cos(\ell\pi v)\) for \(\ell>0\) is
\[
 \frac{2(-1)^\ell-1+8\cos(\ell\pi/2)-9\cos(2\ell\pi/3)}
      {(\ell\pi)^2},
\]
by subtracting (SP4) from the Fourier coefficient of the
triangular density \(4\min(v,1-v)\). Its zeroth integral is zero.
Write \(m_j\) for its moments against \(\cos^j(\pi v)\). On
folding at \(1/2\), the even density is \(-v\) up to \(1/3\)
and \(8v-3\) thereafter, and the odd density is \(3v\) and
\(3(1-2v)\), respectively. The even density has zero integral
and one negative-to-positive crossing against the decreasing
even powers, so \(m_{2j}\le0\) for \(j\ge1\). The odd density
plus the even density times \(\cos(\pi v)\) is nonnegative.
For \(1/3\le v\le3/8\) this follows by using \(\cos(\pi v)\le1/2\)
in the negative term. Hence \(m_{2j+1}+m_{2j+2}\ge0\).
The Fourier formula gives
\[
 \pi^2(m_1,m_2,m_3,m_4)=(3/2,-5/16,19/24,-53/256).
\]
For \(0<a\le5\), all subsequent odd-even pairs in the exponential
series are nonnegative, since the even term within a pair has
relative factor at most \(a/6\). It follows that
\[
 \frac{\pi^2 W(a)}a\ge \frac32-\frac{5a}{32}
       +\left(\frac{19}{144}-\frac{265}{6144}\right)a^2
       >\frac{143}{100},
\]
where the last bound follows by completing the square. For \(a\ge5\)
the signed density changes sign once at \(4/7\). Subtract
the kernel at the crossing to make its product with the density
nonnegative. Restricting to \(v\le1/4\) then gives
\(W(a)\ge(e^{.7a}-1)/32\). Its ratio to \(a\) is increasing
and already exceeds \(143/(100\pi^2)\) at \(a=5\).
For \(s=p+q\), this proves
\[
 B_a(s)\ge B_a(2-x,x,0)+\frac{143a}{200\pi^2}s_j s_k,
           \qquad s_i\ge1,\quad x=s_j+s_k .                  \tag{PH11}
\]

The reserve persists for \(s_i\ge.995\), with constant \(.04a\)
in place of \(143a/(200\pi^2)\). To verify the extension for
\(1\le x\le1.005\), normalize the length of \(U+V\) to one.
The symmetrized test function
\(k_a(xv)+k_a(x(1-v))\) decreases on \(0\le v\le1/2\).
Indeed write \(\theta=\pi xv\le\psi=\pi x(1-v)\).
If \(\psi\le\pi\), then \(\sin\theta\ge\sin\psi\) by the
sine subtraction formula and \(\cos(\pi x/2)\le0\), and
\(\cos\theta\ge\cos\psi\). If \(\psi>\pi\), its sine is
negative. The signed-density comparison of (SV3) therefore still
minimizes the expectation at equal interval lengths. Put
\(T_0=U_0+V_0\); it suffices to bound
\(W_x(a)=\mathbb E k_a(xT_0)-P_\sharp(a)\).

For \(0<a\le5\), the negative \(x\)-derivative of the expectation
in \(W_x\) is at most
\(a\pi^2 x\mathbb E[T_0^2 e^{a\cos(\pi T_0)}]\).
To see this ignore the nonnegative derivative from \(xT_0>1\),
use \(\sin(\pi xT_0)\le\pi xT_0\), and use that the cosine
decreases between \(T_0\) and \(xT_0\le1\). The expectation in
this bound is convex in \(a\). At zero it is \(7/24\). At five,
bound its integral for \(T_0\le1/2\) by
\(4e^5\int_0^\infty t^3 e^{-20t^2}dt=e^5/200<.745\), using
\(1-\cos(\pi t)\ge4t^2\) up to \(1/2\). The remaining part
is at most \(1/2\). Thus the derivative bound is \(12.5a\),
and \(W_x(a)\ge W(a)-12.5a(x-1)>.08a\).

For \(a\ge5\), fold the density of \(xT_0\) at \(1\) using
the cosine symmetry there. The density on \([0,1]\) is
\(4v/x^2\), \(4(x-v)/x^2\), and \(8(x-1)/x^2\) on the
successive intervals cut by \(x/2\) and \(2-x\). Subtracting
the \(P_\sharp\) density gives signs \(+,-,+\), with crossings
\(v_1\in[.55,.6]\) and \(v_2\ge.8\). These intervals follow
from the linear densities on each piece. At \(x=1\) take
the second crossing to be \(1\). The signed first cosine
moment is positive, since its reduction from \(3/(2\pi^2)\)
at \(x=1\) is at most \(\pi(x-1)/2\). Subtract the affine
secant of \(e^{ac}\) across \(c_1=\cos(\pi v_2)\),
\(c_2=\cos(\pi v_1)\). Its signed integral is nonnegative
by the first-moment bound and the positive secant slope.
The remainder has the sign of the signed density by convexity.
On \(v\le1/4\) the density there is at least \(.92v\),
and the secant is at most \(4e^{-.1a}\), since
\(-.32<c_2<-.1\) and \(c_1<-.8\). Thus
\[
 W_x(a)\ge\frac{.92}{32}(e^{.7a}-4e^{-.1a})>.08a
 \qquad(a\ge 5).
\]
The last ratio increases for \(a\ge5\) and satisfies the
bound at five. By (PH10) this proves the stated extension.

In particular the pair consisting of a pure label and a half
split of the other two labels receives a strictly positive
additional remainder in (PH11). The older same-label extraction
vanishes for that pair. The reserve comes from \(B_a(s)-B_a(2-x,x,0)\),
so it can be added to the edge bound itself, the negative Taylor
slack, and the first angular moment compatibility loss. Uniformly
pure slices have zero mixed reserve, so the separate purity
bound (PH9) is a relevant restriction on any high-stability
candidate. A sign comparison for all the remaining profiles is
an open obligation. The analytic proof gives all domain bounds.
The [rational checker](check_three_cell_purity_of_high_stability.py)
replays the numerical inequalities using integer arithmetic and
rational enclosures. The evidence has internal checking only.
