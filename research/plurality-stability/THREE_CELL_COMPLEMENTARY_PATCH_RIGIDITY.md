# Complementary patches of three-cell maximizers in distinct planes

Fix \(0<\rho<1\). Let \(T\) be a balanced regular simplex
maximizing Gaussian noise stability, and let \(A\) be another
ordinary balanced three-cell maximizer. Suppose their active
planes \(V_T,V_A\) are distinct. If \(B,C\) are also ordinary
balanced maximizers and
\[
                         1_A+1_T=1_B+1_C                     \tag{CP1}
\]
as label-vector indicators almost everywhere, then
\(\{B,C\}=\{A,T\}\) with the given labels, up to null sets.

In particular put \(e=1_A-1_T\). If \(E\) is a measurable cut
with \(\int_E e\,d\gamma=0\) and both \(E\) and \(E^c\) contain
positive-measure disagreement, then
\[
           \sum_i\langle 1_E e_i,T_\rho(1_{E^c}e_i)\rangle>0.  \tag{CP2}
\]
This applies to arbitrarily small nonzero tilts of a maximizing
plane in a larger ambient Gaussian space. The argument leaves cuts
whose two original partitions have the same active plane open.

## Rigidity of the two-term identity

Dimension reduction gives an active plane for every ordinary
three-cell maximizer. The plane is genuinely two-dimensional,
and every individual cell uses both active directions.
We use these results for the three sources apart from \(T\).
The dimension reduction theorem is also what permits cylindrical
extension of actual optima without changing their status.
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)
[Full active dimension for every maximizing cell](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions)

We need two distributional observations in an active plane.

**A simplex wedge.** Let \(W\) be a pointed simplex cell in
\(\mathbb R^2\), with independent vectors \(r_1,r_2\) generating
its two boundary rays. The Fourier transform of \(1_W\), off the
two lines \(\xi\cdot r_1=0\) and \(\xi\cdot r_2=0\), is the
nonvanishing smooth function
\[
                  \frac{|\det(r_1,r_2)|}
                   {(i\xi\cdot r_1)(i\xi\cdot r_2)} .
\]
Indeed in coordinates \(x=t_1r_1+t_2r_2\) on the cone, multiply
by \(e^{-\epsilon(t_1+t_2)}\) and use the product of exponential
integrals on two positive half-lines. The damped indicators
converge in tempered distributions as \(\epsilon\downarrow0\),
and the displayed transforms converge uniformly on compact
sets off the two lines. Hence no product of finitely many
nonzero directional derivatives on the plane can annihilate
\(1_W\). Its Fourier multiplier is nonzero on some open set
away from the two lines, where the transform is nonzero.

**A binary indicator with a vanishing mixed derivative.** A
bounded indicator \(f\) on \(\mathbb R^2\) satisfying
\(D_aD_b f=0\) in distributions for nonzero \(a,b\) is
cylindrical in some direction. If \(a,b\) are parallel, smooth
by a compact mollifier. The smoothed function is affine in
their direction on each full line, and boundedness makes it
constant there. This invariance passes to \(f\).
If \(a,b\) are independent, use them as coordinates. The
smoothed functions satisfy the zero rectangular difference
identity, which passes by local \(L^1\) convergence to almost
every choice of the four rectangle coordinates. Fubini
then gives \(f(s,t)=p(s)+q(t)\) almost everywhere, by fixing two
generic base coordinates in that identity. If \(p\) is not
essentially constant, choose two typical \(s\) values at which
\(p\) differs. The two binary values then differ by a fixed
nonzero amount for almost every \(t\), so both are fixed and
\(q\) is essentially constant. If \(p\) is constant, the
conclusion is immediate.
[Distributional cylinder argument for pointwise-distinct triples](THREE_CELL_LATIN_MAXIMIZER_RIGIDITY.md#a-pointwise-distinct-triple-lies-in-one-active-plane)

The ambient distributional derivatives in (CP1) are legitimate
because Gaussian and Lebesgue null sets agree in any finite
ambient dimension. A derivative of a cylindrical function
acts only through the projection of its direction onto the
active plane. Choose \(z\perp V_A\) with a nonzero projection
onto \(V_T\), which is possible since the two planes differ.
If \(V_T\) differed from both \(V_B,V_C\), we could likewise choose
\(v\perp V_B\) and \(w\perp V_C\) with nonzero projections onto
\(V_T\). Applying \(D_zD_vD_w\) to (CP1) would annihilate the
three other sources and leave the asserted vanishing of
three nonzero directional derivatives of a simplex wedge.
This is impossible by the first observation. Interchange
\(B,C\) if needed so that \(V_B=V_T\).

If \(V_C\ne V_A\), choose \(v\perp V_T\) and \(w\perp V_C\)
with nonzero projections onto \(V_A\). Now \(D_vD_w\) applied
to (CP1) would annihilate every term except \(1_A\).
The second observation contradicts the full active dimension
of its cells. Thus \(V_C=V_A\).

The vector difference \(1_B-1_T=1_A-1_C\) therefore depends
only on \(V_A\cap V_T\). Indeed it is invariant under
translations in both orthogonal complements and hence under
their sum. The sum is the orthogonal complement of the
intersection. For a locally integrable representative this
invariance modulo null sets follows, for example, by first
convolving with compact smooth mollifiers, then passing to the
local \(L^1\) limit. If the intersection is zero-dimensional,
the difference is constant, and balance makes it zero.
If it is a line, each coordinate difference in the \(T\)
plane is a function of just that line coordinate. A value
of \(-1\) on a positive-measure set of line coordinates would
make the corresponding \(T\) cell contain an entire
transverse line up to one-dimensional null sets on almost
every such fiber. Its closed pointed wedge contains no full
affine line. The complement of that wedge on the fiber is
relatively open, so every such fiber has positive length
outside the wedge. The coordinate difference is therefore
nonnegative almost everywhere. Balance makes it zero.
This proves (CP1)'s rigidity.

## Strict cost of a nontrivial cut

For the cut in (CP2) form \(B=T+1_E e\) and \(C=T+1_{E^c}e\)
as indicator vectors. Both are ordinary and balanced. Expansion
of the quadratic functional gives
\[
    S_\rho(A)+S_\rho(T)-S_\rho(B)-S_\rho(C)
           =2\sum_i\langle 1_E e_i,T_\rho(1_{E^c}e_i)\rangle.
                                                               \tag{CP3}
\]
Each patched stability is at most the common maximum. A zero
right-hand side would make both patches actual maximizers and
(CP1) would force a trivial cut on the disagreement set.
This proves strictness.

There is a positive lower bound in (CP2), uniform in the
ambient dimension, when \(\rho\) lies in a fixed compact
positive-correlation band, \(\|P_{V_A}-P_{V_T}\|_{\rm op}\)
is bounded below by a fixed positive number, and each of
\(\gamma(E\cap\{A\ne T\})\), \(\gamma(E^c\cap\{A\ne T\})\)
is bounded below by a fixed positive number. If this failed,
pass to a sequence with correlations converging in the band
and with the right side of (CP3) tending to zero. The two
patches are near-maximizers by (CP3). Apply dimension-free
near-maximizer compactness with the varying correlations to
each patch and to the original \(A\). The limiting partitions
are actual maximizers at the limiting correlation.
[Dimension-free near-maximizer compactness](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)

For clarity, the subsequential limits can all be recovered
in a single fixed Gaussian space. The original two planes
span at most four directions. The two near-maximizing patches
are within \(o(1)\) in label distance of cylindrical copies
of fixed planar limit maximizers, adding at most four more
directions. After padding by unused directions and rotating
this span into \(\mathbb R^8\), pass to convergent frames
for all four planes. This does not need the cuts themselves to
be cylindrical. For the original \(A\), compactness is applied
in its two active coordinates, so its limiting frame comes
from the actual original plane. Fixed planar indicator tuples
vary continuously in Gaussian label distance under convergent
orthogonal frames, by approximating them in \(L^1\) by bounded
continuous functions and using invariance of the marginal
standard Gaussian. The limiting \(A,T\) planes remain distinct
by the projector bound. Their two limiting patches obey
(CP1), because the identity held before the \(o(1)\)
approximations. The lower bounds on disagreement ensure \(B\)
stays at positive label distance from both original sources.
This contradicts the rigidity already proved.

Neither the strict comparison nor its compactness consequence
gives an angle-independent bound as the planes approach one
another.
