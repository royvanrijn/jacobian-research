# Four-cell tests of the finite reflection comparison

For balanced four-cell Gaussian voting, the original
[(FR4)](FINITE_REFLECTION_TEST_FOR_MAXIMIZERS.md#a-finite-sufficient-comparison-on-a-positive-correlation-band)
passes on the compact band \([2^{-8771},2^{-8770}]\) with
\(N=2^{8420}\) and reflection tolerance \(2^{-4200}\).
The sign follows analytically over the entire moment body.
A first-moment distance constraint reduces the near-zero comparison
to degree one and gives the explicit sufficient range
\(0<\rho\le2^{-1655}\), uniformly for \(d\ge3\).
This quantifies a portion of the known qualitative initial interval.

At moderate correlations the tail padding can defeat the comparison.
An explicit anisotropic tetrahedral fan with all six interface
adjacencies fails the tests through degree three for \(1/2\le\rho<1\) and
reflection tolerance at most \(1/50\), even with the first-moment
constraint. The strict comparison on unresolved bands remains open.

## The simplex tail can enter the distant set

Write \(T\) for the regular tetrahedral partition, \(P^T\) for its
moment tuple, and \(\tau_N(T)\) for its Hermite weight beyond degree
\(N\). Its exact transposition reflections give \(D_N(P^T)=0\).
For every finite \(N\) and \(0<\rho<1\) the padded value at \(P^T\)
is strictly larger than \(S_\rho(T)\), since the simplex indicators
have positive weight beyond degree \(N+1\). A function with only
finitely many Hermite coefficients is a polynomial, and a
nonconstant Gaussian indicator cannot equal a polynomial almost
everywhere. Thus \(P^T\) prevents FR4 whenever
\(4\tau_N(T)\ge\varepsilon\) and the tail cap admits it. The cap
from the proof of (FR2) admits \(T\) as well: its noise scores
select the simplex by reflection, and its label contrast on
each pair union is already that union's median split along
the pair normal.

Here is a lower bound on the necessary degree. Represent the
simplex by the winner of four independent normals, adding a
common-coordinate Gaussian to the active representation. Fix
\(i\ne j\), put \(t=(Z_i-Z_j)/\sqrt2\), and use the orthogonal
remaining coordinates as tangential coordinates. On the event
that the tied average \((Z_i+Z_j)/2\) exceeds both spectators,
the \(t\)-fiber of \(T_i\) is exactly \(1_{\{t>0\}}\).
That tangential event has measure \(p=\mathbb E[\Phi(Y)^2]\ge1/4\),
where \(Y\) is normal of variance \(1/2\), by symmetry and Jensen.
The squared degree-\((2r+1)\) coefficient of the halfline
indicator in normalized Gaussian Hermites is
\((2\pi)^{-1}\binom{2r}{r}/(4^r(2r+1))\), by integration
by parts. The terms of degree \(>N\) in \(t\) belong to those
of total degree \(>N\). Thus with \(M=\lceil N/2\rceil\ge1\),
\[
 \tau_N(T)\ge 4p\sum_{r\ge M}
          \frac{\binom{2r}{r}}{2\pi\,4^r(2r+1)}
       \ge \frac1{6\pi\sqrt M}.
 \tag{FB1}
\]
For the last inequality use \(2r+1\le3r\),
\(\binom{2r}{r}/4^r\ge 1/(2\sqrt r)\) by induction, and
\(\sum_{r\ge M}r^{-3/2}\ge2/\sqrt M\). The raw threshold therefore
requires \(M>(2/(3\pi\varepsilon))^2\) before the simplex tuple
itself can be discarded.

## A first-moment constraint on the distant set

Use \(n=k-1\), \(F=F_k=W_1(T_k)\), \(B=P_1\) as the row
matrix of first moments, and
\[
 \Gamma(B)=F-\sqrt{F/n}\,\|B\|_* ,
\]
where \(\|\cdot\|_*\) sums the singular values. A singular value
decomposition aligns the regular moment matrix \(B_*\) so that
\(\Gamma=F-\langle B_*,B\rangle_F\). This is the exact linear
score regret against the aligned simplex. The construction
applies to balanced fractional realizers of the moment body.

The top score gap \(G\) of this regular moment simplex satisfies
\(\Pr(G\le t)\le nt\). In the independent normal representation,
the top spacing density at \(s\ge0\) is
\[
 k(k-1)\int\phi(x+s)\phi(x)\Phi(x)^{k-2}\,dx .
\]
Substituting \(y=x+s/2\) shows it decreases on \(s\ge0\).
At zero it equals \(\mu=\mathbb E\max_i Z_i\) by integration
by parts. The regular moment score scale is \(\mu/n\),
so the scaled density is bounded by \(n\).
For \(m=\int(1-f_{T(x)})\,d\gamma\), this gives
\(\Gamma\ge\int_0^\infty(m-nt)_+\,dt=m^2/(2n)\).
It follows that \(d_\gamma(f,1_T)=2m\le\sqrt{8n\Gamma}\).
For ordinary \(A\), each exact simplex reflection has defect
at most \(2d_\gamma(A,T)\) by the triangle inequality.
Consequently a true bad-defect maximizer with \(D(A)\ge\varepsilon\)
must satisfy
\[
             \Gamma(B)\ge\frac{\varepsilon^2}{32n}.
 \tag{FB2}
\]
This necessary constraint removes the simplex tuple and a
neighborhood of its full moment orbit at every degree. The
original FR4 proof remains valid with (FB2) imposed on its
moment body, since any maximizer failing the required
reflection bound must obey it.

## Verification on a positive correlation band

The local comparison proved below is strict for any ordinary
four-cell partition at nonzero orbit distance at most
\(\eta=2^{-650}\) in the active dimension, for
\(0<\rho\le2^{-50}\). Its numerical spanning-tree consequence
uses \(\varepsilon=2^{-4200}\) throughout that same interval.
The tetrahedral first-moment comparison gives an explicit
constant \(K_4^2\le2^{354}\) in
[(G2)](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#gaussian-margin-and-partition-distance).
Indeed (C1) gives \(\alpha_4>2^{-170}\), \(\beta_4>2^{-341}\).
With \(F\in(2^{-4},1)\) this gives \(L_4\le2^{344}\),
\(M_4\le2^5\), and \(K_4^2=8M_4(1+L_4)\le2^{354}\).
The proof of (G2) applies to a balanced fractional realizer
\(f\) of any \(P\) here, using \(1-f_{T(x)}\) as the weight
in the same linear score regret.

Let \(d=d_\gamma(f,1_T)\) for the moment-aligned reference
and \(W_1=\|P_1\|_2^2\). The simplex reflections give
\(D_N(P)\le4d\), by \(\|f-1_T\|_2^2\le d\).
The pointwise estimate \(1-\sum f_i^2\le2(1-f_{T(x)})\)
and the squared triangle estimate on the omitted modes
give \(\tau_N(P)\le2\tau_N(T)+3d\).
So a bad FR4 tuple obeys \(\varepsilon\le16d+8\tau_N(T)\).
The tetrahedron has total Gaussian perimeter \(<6\),
by bounding each of its twelve oriented facets by a whole
central plane. The rotation estimate in (FR2) gives
\(\tau_N(T)\le12/\sqrt{N+1}\).
Take \(N=2^{8420}\), for which \(\tau_N(T)<\varepsilon/16\).
Then \(d\ge\varepsilon/32\), and (G2) gives
\(F-W_1\ge2^{-8764}\). Positivity of the Hermite weights gives
\[
 H_N(P,\rho)+\rho^{N+1}\tau_N(P)
 \le \tfrac14+\rho W_1+\rho^2(3/4-W_1)
 <\tfrac14+\rho F\le S_\rho(T)
 \qquad (0<\rho\le2^{-8770}).
 \tag{FB3}
\]
The strict margin is at least
\(\rho(2^{-8764}-\rho)\), uniform over the entire bad moment body.
This proves the original FR4 on the stated compact band.
The additional maximizer-tail cap was unnecessary for its upper
bound.

For the degree-one version, use the local radius directly:
by (FB2)'s score estimate, any maximizer with
\(\Gamma\le\eta^2/24\) is in the strict local neighborhood.
Thus the finite comparison need only reject moment matrices
with \(\Gamma\ge\eta^2/24\).
Equation (M3) bounds \(\Gamma\le2^{344}(F-W_1)\), so this
set has \(F-W_1>2^{-1649}\). The upper bound in (FB3) at
\(N=1\) gives a strict margin for every \(0<\rho\le 2^{-1655}\).
Dimension reduction reduces any putative maximizer to
the active dimension, giving the stated global conclusion.
Uniqueness includes fractional assignments by the strict
convexity argument used to establish (FR1).

## A four-cell degree-two obstruction with all six adjacencies

Let \(A=A_{1/2}\) and, more generally, let \(A_\lambda\) select the
largest of \(sx+ty+\lambda stz\), \(s,t\in\{\pm1\}\), for
\(\lambda>0\). This is an invertible linear pullback of the tetrahedral
partition. All six interfaces have positive surface measure. For
fixed magnitudes the signs selecting the \(++\) cell are the all-plus
pattern and the pattern with a minus at the least of
\(|x|,|y|,\lambda|z|\). The coordinate sign changes of product one
permute the cells transitively and preserve Gaussian measure.

Write \(\phi_0=\phi(0)\), \(Q=1-\Phi\). The centroid row of
the \(++\) cell is \((x_\lambda,x_\lambda,z_\lambda)\). For a
coordinate given comparison weight ratios \(r,s\) with the other
two coordinates, its value is
\[
 \tfrac14\left(\sqrt{2/\pi}
              -8\int_0^\infty t\phi(t)Q(rt)Q(st)\,dt\right).
\]
Integrating by parts turns this into
\(2\int_0^\infty\phi(t)[r\phi(rt)Q(st)+s\phi(st)Q(rt)]dt\).
By the Gaussian quadrant-angle identity,
\(\int_0^\infty\phi(t)\phi(rt)Q(st)\,dt
=\phi_0(1+r^2)^{-1/2}
 [1/4-\arctan(s/\sqrt{1+r^2})/(2\pi)]\).
Take ratios \(1,2\) for a first coordinate and \(1/2,1/2\) for the
third at \(\lambda=1/2\). Thus
\[
 z_{1/2}=\frac{2\phi_0}{\sqrt5}
      \left(\frac12-\frac{\arctan(1/\sqrt5)}{\pi}\right),
 \qquad
 x_{1/2}=z_{1/2}+\frac{\arctan(1/\sqrt2)}{2\pi^{3/2}}.
\]
Let \(d=x_{1/2}-z_{1/2}\). Elementary rational bounds give
\(z_{1/2}>13/100\) and \(d>11/200\). For the first use
\(447/1000<1/\sqrt5<448/1000\),
\(\arctan(1/\sqrt5)\le(353/375)/\sqrt5<422/1000\),
\(\phi_0>199/500\), and \(\pi>157/50\).
The resulting bound is
\((199/500)(447/500)(287/785)>13/100\).
For the second use
\(\arctan(1/\sqrt2)\ge c/\sqrt2>c(707/1000),\ c=727/840+1/144-1/352\)
and \(\pi<22/7\), \(\sqrt\pi<1773/1000\). Thus
\(d>c(707/1000)/(2(22/7)(1773/1000))>11/200\).
Both arctangent bounds are alternating-series truncations.

All diagonal degree-two coefficients vanish since each magnitude
triple has conditional label probability \(1/4\). The coefficient
of \(xy\) in the \(++\) cell is
\(1/4\) times the expectation of \(|xy|\) where the third weighted
magnitude is least, which equals \(1/(2\pi\sqrt{1+2\lambda^2})\)
by the half-normal truncated first moment
\(\mathbb E[|X|\,1_{|X|>u}]=\sqrt{2/\pi}e^{-u^2/2}\).
The same argument with \(|y|\) least gives the \(xz\) and \(yz\)
coefficients \(1/(2\pi\sqrt{2+\lambda^{-2}})\). Therefore
\(W_2(A_\lambda)=\pi^{-2}\), using
\((1+2\lambda^2)^{-1}+2(2+\lambda^{-2})^{-1}=1\).

For a centroid matrix \(B\) with these four rows and
a transposition \(\sigma\), minimizing the degree-one squared
reflection defect gives
\[
 2\|B\|_F^2-2\,\operatorname{tr}(B^\top\sigma B)
                 +4\lambda_{\min}(B^\top\sigma B).
\]
Indeed \(R_v=I-2vv^\top\), so the maximizing inner product is
the trace minus twice the minimum eigenvalue.
For a pair adjacent in the square of sign indices \((s,t)\),
the matrix \(B^\top\sigma B\) has one diagonal eigenvalue
\(4x_\lambda^2\) and a two-by-two block with zero diagonal and
off-diagonal entries \(\pm4x_\lambda z_\lambda\). The defect is
\(8(x_\lambda-z_\lambda)^2\).
For a square-opposite pair the coordinate swap reflection
preserves the other two rows and the defect is zero.
Any spanning tree uses an adjacent square pair. Thus
\[
 W_1(A)=4(2x_{1/2}^2+z_{1/2}^2)>17/50,\qquad
 D_1(A)=8d^2>1/50,\qquad D_2(A)\ge D_1(A).
\]

For the regular tetrahedron \(\mu=\mathbb E\max_{i\le4}Z_i
=6p/\sqrt\pi\), by integrating the derivative of \(\Phi^3\)
against \(\phi\). Here \(p=\mathbb E\Phi(Y)^2\) as in (FB1),
and \(p=1/4+\arcsin(1/3)/(2\pi)\) by the Gaussian quadrant
formula. For this tetrahedron write \(F=W_1(T)\) and \(G=W_3(T)\).
Its first moments give
\(F=12(1/4+\arcsin(1/3)/(2\pi))^2/\pi<177/500\) by
\(\arcsin(1/3)<17/50\), where \(\sin(17/50)>1/3\)
already follows from the cubic lower bound.
To bound \(G\), write the third-degree coefficients in the sign representation. The unnormalized coefficient
\(\int x_1(x_2^2-1)1_{T_{++}}\,d\gamma\) is \(-u\) with
\(u=2\int_0^\infty t^2\phi(t)^2 Q(t)\,dt\le1/(8\sqrt{\pi})<1/14\).
This follows from the same conditional sign rule and the
truncated second centered moment \(2t\phi(t)\).
If \(b\) is the centroid coordinate of the regular all-plus
cell, conicity and \(\mathbb E R^3=4\mathbb E R\) for a
\(\chi_3\) radius \(R\) give the unnormalized pure degree-three
coefficient \(\int (x_1^3-3x_1)1_{T_{++}}d\gamma=-b+2u\).
The three-index odd coefficient is zero by the sign rule.
Thus \(G=2(b-2u)^2+12u^2\).
Since \(b=(1/4+\arcsin(1/3)/(2\pi))/\sqrt\pi\) lies in
\((1/6,7/40)\), bounding by the larger endpoint value at
\(u=0,1/14\) yields \(G<1/15\).

Set \(\delta=F-W_1(A)<7/500\) and
\(U=3/4-F-\pi^{-2}-G>227/1000\).
The difference between the degree-two padded value for \(A\)
and the degree-three padded upper bound for \(T\) is
\[
 \rho(1-\rho)[\rho^2 U-\delta(1+\rho)].
 \tag{FB4}
\]
For \(\rho\ge1/2\) it is strictly positive by monotonicity of
\(\rho^2/(1+\rho)\) and \(U/4>3\delta/2\).
The tail constraint in FR4 is automatic at \(N=2\). For
every admissible perimeter bound \(C\), the rotation bound
at correlation zero gives \(C\ge(3/2)\sqrt{2/\pi}\).
The FR2 cap is then at least
\((3/4)\sqrt{2/3}/(1-e^{-1})>3/4\).
Thus this tuple defeats FR4 at \(N=2\) for all \(1/2\le\rho<1\)
whenever \(0<\varepsilon\le1/50\). The first-moment condition
(FB2) still admits it, since its raw defect is at least its
degree-one defect. The same witness defeats the degree-one
comparison, which pads all degrees above one by the larger
factor \(\rho^2\).

## The degree-three test on the same cone

The same \(\lambda=1/2\) witness also defeats degree three.
Use comparison weight ratios \(r=w_i/w_j,\ s=w_i/w_l\) for
\(w=(1,1,\lambda)\), where \(i,j,l\) are the three coordinate
indices. The \(++\) cell's unnormalized coefficient
\(B_{ij}=\int x_i (x_j^2-1)1_{A_{++}}\,d\gamma\) equals
\[
 B_{ij}=-2r\int_0^\infty t^2\phi(t)\phi(rt)Q(st)\,dt
       =-\frac{r}{\pi}J_2(1+r^2,s),
\]
by the same conditional sign rule. Here
\[
 J(A,s)=\frac{\sqrt{2\pi}}{\sqrt A}\left(\frac14
                           -\frac{\arctan(s/\sqrt A)}{2\pi}\right),
 \qquad
 J_2(A,s)=-2\partial_AJ(A,s).
\]
The unnormalized pure coefficient is
\(-b_i-\sum_{j\ne i}B_{ij}\) by the \(\chi_3\) radius identity.
The three-odd coefficient vanishes. So the weight of degree
three is
\[
 W_3(A)=4\sum_i \left[
   \frac{(b_i+\sum_{j\ne i} B_{ij})^2}{6}
            +\frac12\sum_{j\ne i}B_{ij}^2\right].
\]
For the regular tetrahedron the all-even degree-four
coefficients vanish by the conditional sign rule. The remaining
unnormalized coefficients have types \((3,1,0)\) and \((1,1,2)\).
They both equal \(-1/(3\pi\sqrt3)\): integrate the two odd
coordinates above the least magnitude \(t\), leaving
\((t^2-1)e^{-t^2}/(2\pi)\) in the half-normal law of the least
coordinate. Including permutations, the factorial normalizations
and the four cells gives \(W_4(T)=10/(27\pi^2)\).

Write \(\delta_1=W_1(T)-W_1(A)\),
\(\delta_3=W_3(A)-W_3(T)\), and
\(R=3/4-W_1(T)-W_2(T)-W_3(T)-W_4(T)\).
The difference of the degree-three padded value at \(A\)
and the degree-four padded bound at \(T\) is
\[
 \rho(1-\rho)\big[\rho^3R-\delta_1(1+\rho+\rho^2)+\delta_3\rho^2\big].
 \tag{FB5}
\]
The [rational checker](check_four_cell_fr4_benchmark.py)
and its [verification record](four_cell_fr4_benchmark_verification.json)
give \(\delta_3>13/2500\), \(R>1/5\), and a bracket greater than \(7/500\) at \(\rho=1/2\).
Because both \(\rho^2/(1+\rho+\rho^2)\) and
\(\rho^3/(1+\rho+\rho^2)\) increase, the difference is positive
throughout \(1/2\le\rho<1\). The tail constraint is automatic
at \(N=3\) by the same rotation-bound lower estimate on \(C\).

## The explicit local exclusion

Let \(A\) be ordinary and balanced in active dimension three, with
\(d_4(A)\le2^{-650}\) and \(0<\rho\le2^{-50}\).
Use the nearest regular moment matrix \(B_*\) as the reference, rotating it to the usual active \(E=\mathbf1^\perp\). Write \(F=F_4\), \(\ell=\sqrt{2F/3}\), \(v=\ell/4\). The regular flux identity and \(\mathbb E\max_{i\le4}Z_i\ge1/\sqrt\pi\) give \(\ell\ge 2\phi_0/3>1/4\), hence \(v>1/16\), \(F>2^{-4}\). Also \(F\le1\) by projection. Write
\[
 h_{ij}=1_{T_i\cap A_j},\quad t=x\cdot n_{ij},\quad
 e_i=1_{A_i}-1_{T_i},\quad L=\sum_{i\ne j}\int U_{\rho,ij}h_{ij}\,d\gamma,\quad s^2=L/\rho.
\]
For the lower score bound use the reflected kernel
[(T3)](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#global-score-margins-and-facet-normal-derivatives).
If \(m=\rho x/\sqrt{1-\rho^2}\in T_i\), restrict it to \(m+z\)
with \(z\in T_i\) and \(z\cdot n_{ij}\ge1/2\). The retained
Gaussian mass in \(z\) is at least \(1/4-\phi_0/2>1/20\).
The kernel ratio is at most \(e^{-m\cdot n_{ij}}\), so concavity
on \([0,1]\) gives \(U_{\rho,ij}/\rho\ge (1/64)\min(t,1)\).
The derivative lower bound
[(B7)](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation)
gives \(\partial_n U_{\rho,ij}/\rho\ge\ell\) on the facet. The standard reflected upper bound gives \(U_{\rho,ij}/\rho\le t\) since \(2\phi_0/\sqrt{1-\rho^2}<1\). Thus
\[
 \sum_{ij}\int\min(t,1)h_{ij}\le64s^2,\qquad
 s^2\le (F-\langle B_A,B_*\rangle)/\ell\le 4\sqrt{d_4(A)}.
\]
The last step uses first-Hermite projection and nearest alignment, so \(\|B_A-B_*\|_F\le\sqrt{d_4(A)}\). Hence \(s\le2^{-160}\) under the assumed radius. The case \(s=0\) is exactly the reference simplex by the positive score margins. Otherwise \(H=(B_A-B_*)/s\) is symmetric in \(E\), by the nearest-moment condition.

Put \(w(x)=(1+|x|^2)e^{|x|^2/16}\). For \(m\le3\),
\[
 \int_{|x|>R}w^2\,d\gamma_m\le 2^{11}e^{-R^2/16}.
\]
Indeed \((1+u)^2\le 256 e^{u/8}\) (the exact maximum of \((1+u)^2e^{-u/8}\) occurs at \(u=15\)), and inserting \(e^{(|x|^2-R^2)/16}\) leaves Gaussian moment \(\mathbb E e^{5|X|^2/16}=(8/3)^{m/2}<5\). For \(t<1\), \(w(z+tn)\le4w(z)\), \(|x|>R\) implies \(|z|>\sqrt{R^2-1}\), and the squared normal mass is at most its first normal moment by one-dimensional filling, since \(2\phi_0<1\). For \(t\ge1\) use Cauchy–Schwarz against \(h_{ij}\). Combining the two bounds and Cauchy–Schwarz over the twelve ordered transfers gives
\[
 \sum_{i\ne j}\int_{|x|>R}w h_{ij}\,d\gamma
 \le 2^{15}s e^{-(R^2-1)/32},\quad R\ge1,\qquad
 \sum_{i\ne j}\int w h_{ij}\,d\gamma\le2^{15}s.
\]
For example each split costs at most \(4\sqrt{2^{11}}+\sqrt{2^{11}}\) times the square root of the corresponding normalized mass-cost, and \(\sqrt{12\cdot64}<2^5\).

Let \(R=64,\ h=2^{-80}\). Keep transfers on \(|x|\le R,\ 0\le t\le h\) whose foot \(z\) belongs to the actual facet. Let \(\theta_{ij}\) be the projected kept measure divided by \(s\,d\nu\), let \(f_{ij}=\theta_{ji}-\theta_{ij}\), and identify the resulting signed labelled interface measures and their moments with \(f\) and \(\mathcal Df\). The following upper bounds apply to the discrepancy between the original normalized signed transfers and the projected ones, before multiplying by \(\|e_j-e_i\|=\sqrt2\):
- Exterior normalized mass and first moment cost at most \(2^{15}e^{-(R^2-1)/32}<2^{-128}\).
- Within \(R\), the normalized mass and first moment cost for \(t>h\) is at most \((1+R)64s/h\le2^{13}s/h\le2^{-67}\).
- For \(t\le h\) and off-facet foot, a spectator score violates a facet inequality by \(\le h\). For the tetrahedron the projected score difference in the tie plane has variance \(3/2\). The union of at most two width-\(h\) strips thus has \(\gamma_{n-1}\)-measure at most \(2h\). Normal filling gives \(\sum_{\mathrm{ordered}}\int \theta_{ij}^2\,d\nu\le (2/s^2)\sum\int_{t<1}t h_{ij}\le128\) even before removing the bad strips. The total normalized mass discarded is at most \(\sqrt{128\cdot24h}<2^6\sqrt h\); include \(1+R\) for mass plus first moment, giving \(2^{13}\sqrt h\).
- For the kept projection, the normalized first moment changes by at most \(2^{15}h\).

It follows conservatively that the Euclidean row mass vector \(c\) of \(f\) has norm \(\|c\|_2\le2^{-24}\), and \(\|\mathcal Df-H\|_F\le2^{-24}\). The reference masses and actual masses both equal \(1/4\). The actual moment matrix error is centered and symmetric; no balance or rotation correction of \(f\) is required.

For the kept transfer, the normal coefficient at the foot is at least \(\ell\) after dividing by \(\rho\). Directional differentiation of the Gaussian smoothing twice shows \(|\partial_{nn}U_{\rho,ij}|/\rho \le 2\rho/(1-\rho^2)\le1\), using \(|1_{T_i}-1_{T_j}|\le1\) and \(\mathbb E|Z^2-1|\le2\). Hence \(U_{\rho,ij}/\rho\ge(\ell-h)t\) there. Normal filling gives \(\sum_{\mathrm{ordered}}\int \theta_{ij}^2\,d\nu\le2/(\ell-h)\). As \(\theta_{ij}\ge0\), summation over unordered faces yields \(\|f\|^2\le2/(\ell-h)\).

Use the proof of [(L2)](SMALL_POSITIVE_CORRELATION_INTERVAL.md#an-interface-estimate-in-the-moment-gauge) with the *symmetric test matrix* \(H\), rather than \(\mathcal Df\). The conditional variance estimate after subtracting face means gives a pairing at most \(\sqrt{3\ell/4}\,\|f\|\|H\|_F\). The conditional face means are \(\alpha(H_{ii}-H_{jj})\), since \(\mathbb E[z\mid F_{ij}]=\alpha P(e_i+e_j)\) and \(\|P(e_i+e_j)\|=1\) for \(k=4\). Here \(|\alpha|\le\sqrt{\phi_0/v}<4\) by Cauchy–Schwarz in the supporting standard Gaussian plane. Thus their pairing is bounded by \(4\|c\|_2\|H\|_F\). With the moment discrepancy,
\[
 \|H\|_F\le\sqrt{\frac{3\ell/2}{\ell-h}}+4\|c\|_2+\|\mathcal Df-H\|_F,
 \qquad \|H\|_F^2<7/4.
\]
This is immediate as \(h/\ell\le2^{-78}\), \(5\cdot2^{-24}<2^{-20}\).

Fix \(r=1/16\). The exact Mehler bound \(M_r(x,y)\le2e^{(|x|^2+|y|^2)/32}\le2w(x)w(y)\) and \(\sum_i\int w |e_i|\le2^{16}s\) imply \(\sum_i\langle e_i,T_r e_i\rangle\le2^{33}s^2\). Since \(\rho^m\le\rho^2r^{m-2}\) for \(m\ge2\), the order-\(\ge2\) noise remainder divided by \(\rho s^2\) is \(\le 2^{41}\rho\le2^{-9}\). The first order term divided by \(\rho s^2\) is \(\|H\|_F^2<7/4\); the constant term vanishes by balance. Therefore \(2L-\sum_i\langle e_i,T_\rho e_i\rangle>0\) for \(s>0\).

## The explicit spanning-tree tolerance

For the balanced reflected repair \(B\) of a single approximate edge of error \(\delta\), use the existing repair \(d(A,B)\le\delta\). The gain from any polarization of \(B\) is at most \(2\rho\sqrt\delta\): by balance, write \(S_\rho(A)-S_\rho(B)=\langle A-B,(T_\rho-\Pi_0)(A+B)\rangle\) and use operator norm \(\rho\), difference norm at most \(\sqrt\delta\), and sum norm at most 2.

Retain the sign-majority set above height \(1/(16\phi_0)\ge1/8\), with mass \(\ge1/16\), and sign-minority set above height \(m/(2\phi_0)\ge m\), with mass \(\ge m/2\), where \(m\) is the minority mass. The positive kernel difference on this product is at least \((\rho m/8)M_\rho\), using \(1-e^{-\rho m/4}\ge \rho m/8\). [Reverse Gaussian hypercontractivity](https://arxiv.org/pdf/1306.2410v2) with \(p=q=1/2\) applies since \((1-p)(1-q)=1/4\ge \rho^2\), giving joint mass at least \(2^{-10}m^2\). The exact polarization gain (RG1) is at least \(2^{-13}\rho m^3\) (discarding extra factors). Therefore \(m\le 2^5\delta^{1/6}\).

Orient the edge by the majority. For \(m\le2^{-18}\), the good side of \(B_i\) above height \(1/8\) has mass at least \(1/8\) whereas the wrong-side first moment has absolute value at most \(\sqrt m\); \(b_i(B)-b_j(B)\) is parallel to the oriented normal and has length \(>1/64\). Individual centroid errors of the repair are at most \(\sqrt\delta\), testing against unit directions. When \(2\sqrt\delta<1/128\), the A-centroid distance is \(>1/128\) and its unit normal changes by at most \(2^8\sqrt\delta\), using the normalized vector difference inequality. Gaussian central halfspace symmetric differences are at most their unit normal difference. Consequently the total wrong-side mass of \(A_i,A_j\) relative to their own pair-centroid direction is \(\le 2m+\delta/2+2^9\sqrt\delta \le 2^{11}\delta^{1/6}\).

Each label transposition can be expressed in \(\le5\) edges of a four-vertex spanning tree (along a path of length at most three). Any permutation requires at most three such transpositions, hence at most 15 edges. Conjugating one tree reflection by such compositions makes every edge error at most \(31\varepsilon\). The preceding bounds apply with \(\delta\le32\varepsilon\), and the same compositions show that all pairwise centroid lengths differ from a fixed one \(\Delta\) by at most \(2\sqrt{15\varepsilon}<8\sqrt\varepsilon\). Every centroid length is at most 1 (Cauchy–Schwarz against a unit direction), so the squared distance matrix differs from the equidistant one entrywise by at most \(32\sqrt\varepsilon\). After double centering the Gram error is at most \(64\sqrt\varepsilon\) in Frobenius norm. Since \(\Delta>1/128\), the aligned regular matrix with distance \(\Delta\) differs from the A centroid matrix by at most \(2^{14}\sqrt\varepsilon\), using scalar square-root comparison against \((\Delta^2/2)P\). The resulting unit pair normal differences are at most \(2^{23}\sqrt\varepsilon\). Summing the wrong-side and angle errors over six pairs, including the factor two for labelled distance, gives \(d_4(A,T)\le2^{29}\varepsilon^{1/6}\). With \(\varepsilon=2^{-4200}\) this is at most \(2^{-671}<2^{-650}\), so the above strict local inequality applies to any non-simplex at that reflection threshold.

## Additional constraints from exposed moment tuples

The cone in (FB4) witnesses a failure of the padded comparison.
Its moment tuple can be excluded from the set of actual maximizers
by the following finite-dimensional constraints.

First, this cone violates the actual score equation on an adjacent
interface. The face between labels \((+,+)\) and \((-,+)\) has
normal \((1,0,\lambda)\). If the cone were a maximizer, the score
difference would be constant on the open face. Analyticity
extends that equality to the entire pair plane. Its value at zero
is zero by conicity and balance. Differentiating tangentially at
zero would then give \(b_i-b_j\) parallel to \((1,0,\lambda)\):
indeed \(\nabla T_\rho(1_{A_i}-1_{A_j})(0)=
(\rho/\sqrt{1-\rho^2})(b_i-b_j)\).
For \(\lambda=1/2\) the moment difference has direction
\((x_{1/2},0,z_{1/2})\), with \(x_{1/2}\le\phi_0/2<1/5\)
and \(z_{1/2}>13/100\). These directions differ.

Moreover this failure excludes an open first-moment neighborhood
of that cone from actual maximizers, uniformly on any compact
positive band. The cone's tuple is exposed by its own linear
scores, whose pointwise winner is unique almost everywhere.
If a sequence of realizers has first moment matrices converging
to the cone's, their regrets against those linear scores converge
to zero. Splitting into an arbitrarily thin score strip and its
complement shows convergence to the cone in Gaussian label
distance. Stability is continuous in that distance and uniformly
Lipschitz in correlation on a compact band, by the Hermite series.
The maximum stability is also continuous in correlation, using
the same uniform Lipschitz bound. A sequence of actual maximizers
approaching the cone tuple would therefore make it optimal at a
limiting band point, contradicting its score equation.

At \(\rho=1/2\) the excluded first-moment neighborhood can be
taken explicitly to have Frobenius radius \(2^{-136}\) about the
cone's matrix \(B_A\) in its stated orientation. To see this,
write \(i=(+,+),j=(-,+)\), put
\(v=(-1/2,0,1)/\sqrt{5/4}\), \(s=2^{-8}\), and let
\(x_\pm=s e_y\pm s v\), where \(e_y=(0,1,0)\).
Both points lie on the common face with margin greater than
\(s\) over the spectator scores.
For \(F=T_\rho(1_{A_i}-1_{A_j})\), the derivative
\(\partial_v F(0)>1/40\) by the moment bounds on \(x_{1/2}\)
and \(z_{1/2}\). Gaussian smoothing gives
\(\|\nabla F\|\le1\) and \(\|\nabla^2 F\|_{\mathrm{op}}\le1\).
Indeed the directional bounds are \((\rho/\sqrt{1-\rho^2})
\mathbb E|Z|\) and \((\rho^2/(1-\rho^2))\mathbb E|Z^2-1|\)
for a standard normal \(Z\).
Taylor's inequality gives \(F(x_+)-F(x_-)>s/25\).

Shift \(x_-\) into label \(i\) and \(x_+\) into label \(j\) by
\(r=s/1024\) along the pair normal, and take balls of radius
\(r/2\) about the shifted points. The balls lie in their assigned
cone cells by the spectator margins and are contained in the
unit ball. Each has Gaussian mass \(>2^{-61}\). The \(F\) values
on the second ball exceed those on the first by \(s/25-3r\)
or more. If \(d_\gamma(A',A)\le2^{-63}\), the two respective
labels of \(A'\) still occur on positive measure in the balls.
On the unit ball, Mehler density relative to the Gaussian
measure in its other variable is at most
\((1-\rho^2)^{-3/2}\exp(1/2)<3\) at \(\rho=1/2\).
Thus replacing \(A\) by \(A'\) changes the score difference there
by at most \(3d_\gamma(A',A)\). The strict score ordering
persists, and a balanced exchange between the two balls violates
the first variation of \(A'\).

For the cone's defining linear scores the Frobenius norm of
the score matrix is \(3\). The top gap has sublevel probability
at most \(8t\) at height \(t\), by the union bound on the
six Gaussian pair-score differences, each of standard deviation
at least \(2\). Integrating this bound as in (FB2) gives
\(d_\gamma(A',A)\le20\sqrt{\|B_{A'}-B_A\|_F}\).
At moment distance at most \(2^{-136}\) this is \(<2^{-63}\).
The exclusion also holds for rotations and relabelings of
the cone tuple.

A second finite constraint can be applied directly. Let
\(W_\rho P=(\rho^jP_j)_{1\le j\le N}\) and let
\(\sigma_{\mathcal K_N}\) denote the exact balanced support
function (FR5). For the moments \(P\) of an actual maximizer,
\[
 2\big(\sigma_{\mathcal K_N}(W_\rho P)
          -\langle W_\rho P,P\rangle\big)
     \ \le\ \rho^{N+1}\tau_N(P).
 \tag{FB6}
\]
For any \(Q\in\mathcal K_N\) its fractional realizers have actual
stability at least \(H_N(Q,\rho)\). Global maximality and (FR3) give
\(H_N(Q,\rho)-H_N(P,\rho)\le\rho^{N+1}\tau_N(P)\). The left
side is at least \(2\langle W_\rho P,Q-P\rangle\), since the
remaining weighted quadratic term is nonnegative. Taking the
supremum proves (FB5). A lower bound from a single realizable
challenger \(Q\) already gives a necessary cut. This allows
constraints from actual maximality to be added to the finite
comparison alongside the first-moment distance.
