# Uniform regular exclusion at the shared quartic transition

Use the fully shared ordinary profiles, curvature \(W,\rho\),
and centered parameters \(w_n,w_{n+1}\) of the
[moving quartic window](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
The varying critical center \(w_*\) there is real. Write
\(\delta=\max\rho-1\).

**Theorem.** On any fixed band \(0<W_-\le W\le W_+<\infty\),
a sequence of consecutive common roots with \(n\to\infty\)
and \(\delta\to0\) must satisfy both
\[
                  w_n-w_*\longrightarrow0,\qquad
                  w_{n+1}-w_*\longrightarrow0 .
\]
In particular, for every finite \(K\) there is one eventual
cutoff on the entire band \(\sqrt n\,|\delta|\le K\).
The total-degree difference is arbitrary and the cutoff
is unevaluated.

## The critical limiting graph

Take a retained and diffuse limit \(G_*\) from the
[compactness theorem](MULTIVARIATE_UNIFORM_SADDLE_MAP_COMPACTNESS.md#retained-and-diffuse-limits).
The curvature converges on \([0,1]\). The limiting
\(\rho_*\) is strictly concave with maximum one, at a
unique interior point \(q_c\), and \(w_*\to w_c=X_*(q_c)\).
Here \(X_*(q)=G_*'(q)+\log(1-q)-\log q\).
The two parameters of a prospective common root have
bounded real parts by the uniform root annulus and
have identical imaginary parts. Work in the upper
half-plane for the original parameter \(t\); the
strip for \(w\) is \(-\pi<\Im w<0\).

The critical vertical line \(\Re q=q_c\) misses all fibers
over this strip. Indeed, adding \(\lfloor\epsilon A\rfloor\)
rows \((A,A+1)\) with \(A\to\infty\) adds diffuse mass
\(\epsilon>0\) to \(G_*\). The limiting curvature at \(q_c\)
becomes \(1+\epsilon q_c(1-q_c)>1\). The finite-profile
[vertical-divider inequality](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md#equal-total-degrees-at-and-above-the-real-cusp)
applies at \(q_c\), since its proof uses the curvature and
the integer gap multiplicities. Pass first to the limit
along the augmented sequence and then let \(\epsilon\)
decrease to zero. This gives \(\Im X_*(q_c+iy)\ge0\)
for \(y>0\).

There are two boundary inverses at each point of
\(\Im w=-\pi\), one in each of the real gaps outside
\([0,1]\), with \(X_*'>0\). The gap ranges are as in the
[limiting supercritical geometry](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#the-limiting-phase-and-its-real-boundary-geometry).
Their local inverses enter the strip on opposite sides
of the divider. Since the interior map is proper of
degree two and the divider cannot be crossed, it has
one simple branch on each side over the whole strip.
Label them \(q_L,q_R\). Put
\[
 f(q,w)=G_*(q)+(q-1)\log(1-q)-q\log q-wq,\qquad
 D(w)=f(q_R,w)-f(q_L,w).
\]
Then \(\Re D'=\Re(q_L-q_R)<0\). At fixed imaginary
part the left branch approaches zero as \(\Re w\to+\infty\)
and the right approaches one as \(\Re w\to-\infty\),
from the endpoint logarithms. The fixed divider gives
a strictly negative bound for \(\Re D'\) at each
horizontal end. Hence there is one equal-modulus point
per height, and these points form a smooth graph \(E\).
The two complementary regions are connected.
The two real gap branches give a unique lower
endpoint with \(D=i\pi\).

At the top, \(X_*(q_c+h)=w_c-c h^3+O(h^4)\) with \(c>0\),
since \(X_*'=(\rho_*-1)/(q(1-q))\) and
\(\rho_*''(q_c)<0\). The local analytic cubic coordinate
shows that, on a small half-disk below \(w_c\), both
upper branches lie in the cusp chart and coalesce
at \(q_c\). On each real radius away from \(w_c\), the
cubic inverses consist of a real branch and a nonreal
upper branch. Their limiting real parts differ by
\(3(|w-w_c|/c)^{1/3}/2+O(|w-w_c|^{2/3})\).
The real inverse continues from the upper half-plane
because its derivative \(X_*'\) is negative.
Integration of \(\Re D'\) horizontally from \(w_c+iy\),
followed by \(y\uparrow0\), gives opposite strict signs
for the limiting \(\Re D\) at \(w_c\pm\epsilon\) for
all sufficiently small \(\epsilon>0\).
Thus \(E\) ends at \(w_c\) and \(D\to0\) there.
On \(E\) the Cauchy--Riemann equations give
\(d(\Im D)/d(\Im w)=|D'|^2/\Re D'<0\).
Consequently \(0<\Im D<\pi\) at every interior point.

## Actual contours and the adjacent comparison

The lower-boundary anchoring and interior finite flow of
the [uniform supercritical contour proof](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#finite-contour-deformation-for-the-limiting-profile)
apply with this limiting graph. In detail, the signs at
real cuts are unchanged and the contour compactness proof
at infinity needs the same strict velocity estimates.
For infinitely many retained rows or positive diffuse mass,
its finite selections give those estimates as stated.
For finitely many retained rows of one sign with zero
diffuse mass, their absolute gaps must still sum to at
least three. A sum of at most two gives
\(\rho_*(q)\le q(1-q)\sum_i d_i^2\le1\);
equality would require a nonzero-gap row to have
\(a_i+d_iq=1\) at an interior \(q\), which is impossible
for distinct positive integer entries. This checks the
only strict hypothesis at infinity that used the previous
supercritical margin. The finite-loop arguments therefore
give both actual Gaussian coefficients on the interior
of \(E\), a single leading nonzero coefficient in each
complementary region, and the two gap saddles near the
lower endpoint. The fixed homotopies transport to the
approximating finite profiles through the same cut-free
compact neighborhoods. Their stationary phases converge
with derivatives, irrespective of the rate at which the
profile moves.

At the upper boundary, away from \(w_c\), the real
limiting entropy phase has a unique nondegenerate maximum.
Indeed its real second derivative is nonpositive and
vanishes only at \(q_c\), and the endpoint derivatives
have opposite infinite limits. Cotangent summation on a
small interval about the maximum gives the one-saddle
exclusion, as in the moving quartic proof. Thus every
limiting parameter of a zero is on the closed graph.

At a common root the two limiting parameters have the
same height and must agree, because the closed graph
has one point per height. If this is an interior point
or the lower endpoint, the
[adjacent Gaussian comparison](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#uniform-comparison-of-consecutive-orders)
applies with \(\Re D'<0\), \(\Delta/n\to0\), and
\(e^D\ne1\). For \(\Delta\ne0\), the real logarithm
of the ratio of the two stationary multipliers has
leading term \(-\Delta\Re D'\), which contradicts the
converging coefficient ratio. For \(\Delta=0\) the ratio
tends to \(e^D\ne1\), with the same contradiction.
These comparisons remain uniform along the chosen
subsequence. The only possible common limit is \(w_c\),
which proves the localization. The fixed-radius moving
quartic theorem then excludes it on \(\sqrt n|\delta|\le K\).
