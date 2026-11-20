# Quartic Lee–Yang zeros and a quantitative Robin branch bound

This is a research-only model theorem. It does not identify an actual Gamma
period with the model, supply a quartic uniform remainder, or prove a new
Strong Factorial case.

Put

\[
 Z(h,\tau)=\int_{\mathbb R}
 \exp(-u^4/12+\tau u^2/2+hu)\,du,
 \qquad F(x,\tau)=Z(ix,\tau).
\tag{Q1}
\]

Here \(\tau\) is real unless stated otherwise. Both variables of \(Z\) are
entire. Differentiation under the integral is uniform on compact subsets.

## 1. Actual classical theorem and its exact specialization

Charles M. Newman, *Fourier transforms with only real zeros*, Proceedings AMS
61 (1976), 245–251, Theorem 2, printed pp.246–248, classifies the finite even
positive measures \(\rho\) such that
\(\int e^{ht-bt^2}d\rho(t)\) has only imaginary zeros for every \(b>0\).
The sufficient class includes
\(d\rho(t)=K\exp(-\alpha t^4-\beta t^2)dt\), with \(\alpha>0\) and
**arbitrary real** \(\beta\). This is the empty-product, \(m=0\) case of
(1.8). The proof of sufficiency is written on printed pp.247–248.

For any given real \(\tau\), take \(\alpha=1/12\),
\(\beta=-\tau/2-1\), and \(b=1\). This is a finite measure regardless of
the sign of \(\beta\). Consequently all zeros of \(Z(\cdot,\tau)\) are
imaginary, for **every real \(\tau\)**. No finite-size scaling argument is
required for this zero-location statement.

The exact primary PDF inspected is
`quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf`, downloaded from
<https://sites.math.northwestern.edu/~auffing/papers/Newman.pdf>.
This is an actual theorem specialization; it is not a claim of novelty.

## 2. Differential identities and simplicity

Integration of the derivative of the integrand gives

\[
 Z_{hhh}=3\tau Z_h+3hZ,\qquad Z_\tau=\tfrac12Z_{hh}.
\tag{Q2}
\]

Equivalently

\[
 F_{xxx}+3\tau F_x-3xF=0,\qquad F_\tau=-\tfrac12F_{xx}.
\tag{Q3}
\]

Every zero is simple. Indeed, a zero of multiplicity at least three would
give zero initial data for the third-order linear ODE in (Q3), forcing
\(F\equiv0\), whereas \(F(0,\tau)>0\). At a hypothetical double zero
\((x_0,\tau_0)\), put \(c=F_{xx}(x_0,\tau_0)/2\ne0\). For real
\(\eta\to0\), uniformly for bounded \(z\),

\[
 \eta^{-2}F(x_0+\eta z,\tau_0-\eta^2)
      \longrightarrow c(z^2+1).
\tag{Q4}
\]

Rouché on small circles around \(i\) and \(-i\) gives nonreal zeros for
\(\tau_0-\eta^2\), contradicting Section 1. Thus every zero of the model
is imaginary and simple.

The zeros form an infinite sequence \(\pm i x_j(\tau)\), \(x_j>0\).
For completeness, the elementary quartic Young inequality bounds the entire
order of \(F\) by \(4/3\). If it had only finitely many zeros, Hadamard
factorization would make it an exponential of a linear polynomial times a
polynomial. Evenness eliminates the nonconstant exponential. But on the
imaginary axis \(F(-iy,\tau)=Z(y,\tau)\) grows faster than every polynomial,
as follows by integrating over a fixed interval translated to a positive
multiple of \(y^{1/3}\). This contradiction proves infinitude. A local
simple-zero branch is real analytic in real \(\tau\).

## 3. Critical Curie–Weiss limit and strictly decreasing zeros

The monotonicity import used here is precisely Theorem 2 of Qi Hou,
Jianping Jiang and Charles M. Newman, *Motion of Lee–Yang zeros*,
arXiv:2208.00917v1, printed p.5, equation (17): the positive principal zeros
of

\[
 \sum_{\sigma\in\{-1,1\}^n}
       \exp(tS^2+ixS),\qquad S=\sum\sigma_j,
\]

are strictly decreasing as \(t\ge0\) increases **when \(n\) is even**. We
use only even \(n\) throughout the limiting argument. (For odd \(n\), the
last principal zero \(\pi/2\) is constant.) This theorem compares
couplings at the **same finite degree**. We apply it at fixed \(n\) and then
take a locally uniform limit; we do not infer interlacing between different
degrees.

Take \(\beta_n=1+\tau/\sqrt n>0\),
\(t_n=\beta_n/(2n)\), and external field \(h/n^{3/4}\). Gaussian
linearization and a translation of its entire Gaussian contour show that,
up to a positive factor independent of \(h\), the partition function is

\[
 e^{-h^2/(2\beta_n\sqrt n)}
 \int_{\mathbb R}\exp\left\{
 n\left[\log\cosh(n^{-1/4}u)
             -\frac{u^2}{2\beta_n\sqrt n}\right]
                  +\frac{hu}{\beta_n}\right\}du.
\tag{Q5}
\]

It converges to (Q1) locally uniformly in complex \(h\), uniformly when real
\(\tau\) ranges over a fixed compact interval, as even \(n\to\infty\).
Here is a domination argument:
for real \(v\), the function
\(\log\cosh v-v^2/2\) is strictly negative away from zero, equals
\(-v^4/12+O(v^6)\) near zero, and is at most \(-v^2/4\) for sufficiently
large \(|v|\). The small-\(v\) part is bounded by
\(\exp(-c u^4+C_Tu^2+C_H|u|)\). The intermediate compact annulus contributes
\(O(e^{-c n})\) times its length, since its negative \(n\)-term dominates
the \(O(\sqrt n)\) temperature correction and \(O(n^{1/4})\) field term.
The large-\(v\) part has a Gaussian upper bound with vanishing tails. These
three bounds justify locally uniform convergence and all fixed derivatives.

Each fixed positive zero in a compact window is the limit of the
correspondingly ordered finite principal zero, by simplicity and Rouché.
There are no uncounted zeros near the origin since \(F(0,\tau)>0\).
Because \(t_n\) increases with \(\tau\), the limiting positive zero branches
are nonincreasing. Write

\[
 r=\frac{F_{xx}(x_j,\tau)}{F_x(x_j,\tau)},\qquad
 q=-r=-2x_j'(\tau)\ge0.
\]

Equations (Q3) give the exact system

\[
 x'=-q/2,\qquad q'=(3x-3\tau q-q^3)/2.
\tag{Q6}
\]

In fact \(q>0\). If \(q(\tau_0)=0\), (Q6) gives
\(q'(\tau_0)=3x(\tau_0)/2>0\), so \(q<0\) just before \(\tau_0\), a
contradiction. Thus each positive zero decreases strictly.

## 4. A quantitative logarithmic-derivative bound at every zero

At \(h=\pm ix_j(\tau)\), one has

\[
 \left|\frac{Z_{hh}}{Z_h}\right|
 \le \max\left\{\sqrt{12(|\tau|+1)},
                    [12(x_j+1)]^{1/3}\right\}
 \le 4(1+\sqrt{|\tau|}+x_j^{1/3}).
\tag{Q7}
\]

Here is a blow-up barrier proof, avoiding an unproved asymptotic estimate.
Suppose that at \((x_0,\tau_0)\) the positive value \(q_0\) exceeds the
first maximum in (Q7). Run (Q6) backwards, writing \(X(t)=x(\tau_0-t)\),
\(Q(t)=q(\tau_0-t)\). As long as \(t<1\), \(X<x_0+1\), and \(Q\ge q_0\),

\[
 \dot X=Q/2,\qquad
 \dot Q=(Q^3+3(\tau_0-t)Q-3X)/2>Q^3/4.
\]

The latter inequality uses
\(Q^2>12(|\tau_0|+1)\) and \(Q^3>12(x_0+1)\). Moreover
\(dX/dQ<2/Q^2\), so \(X-x_0<2/q_0<1\). Integrating
\(d(Q^{-2})/dt<-1/2\) forces blow-up before \(2/q_0^2<1\).
Thus none of the auxiliary bounds can fail first. But \(X\) stays bounded
and increasing and \(\tau\) stays in a compact interval. At the finite
limiting pair, continuity gives a zero of \(F\), and its simplicity makes
\(F_{xx}/F_x\) bounded there. This contradicts blow-up. The conversion from
\(q\) to the absolute quotient of \(Z\) uses \(F_x=iZ_h\),
\(F_{xx}=-Z_{hh}\); negative zeros follow by evenness.

## 5. A quantitative Robin branch, with an explicit scope restriction

Fix a model zero \(h_0=ix_j(\tau)\), real \(\tau\), and
\(K=1+\sqrt{|\tau|}+|h_0|^{1/3}\). For any complex \(\Lambda\) satisfying
\(|\Lambda|K\le1/100\), the root of

\[
 Z(h,\tau)+\lambda Z_h(h,\tau)=0
\tag{Q8}
\]

starting at \(h_0\) continues along \(\lambda=s\Lambda\), \(0\le s\le1\),
and obeys \(|h-h_0|\le2|\Lambda|\). Along that branch put
\(Q=Z_{hh}/Z_h\). Direct differentiation gives

\[
 \frac{dh}{d\lambda}=-\frac1{1+\lambda Q},\qquad
 \frac{dQ}{d\lambda}
   =-\frac{3\tau-3\lambda h-Q^2}{1+\lambda Q}.
\tag{Q9}
\]

The bounds \(|Q|\le5K\), \(|h-h_0|\le2|\Lambda|\) close by continuation:
the denominator is at least \(19/20\); using (Q7), the numerator is less
than \(28.031K^2\), so \(|dQ/d\lambda|<30K^2\). Hence
\(|Q-Q(0)|<0.3K\), while \(|Q(0)|\le4K\); also
\(|h-h_0|\le(20/19)|\Lambda|<2|\Lambda|\).
No finite continuation obstruction remains on this bounded set.

This statement tracks a branch issuing from a genuine model zero. It does
**not** assert that every root of (Q8), with an unrestricted zero index, is
obtained in this manner under a condition stated only at the perturbed root.
Escaping-root exclusion and quadratic-amplitude control remain separate
obligations for an actual quartic period.

### Complex temperature shifts are not covered by the Lee–Yang assertion

There is nevertheless an explicit local continuation bound from each genuine
zero at a **real** temperature. Fix \((h_0,\tau_0)\) as above and take complex
\(\Lambda,\Sigma\) with
\(|\Lambda|K\le1/1000\), \(|\Sigma|K^2\le1/1000\). Along
\(\lambda=s\Lambda\), \(\tau=\tau_0+s\Sigma\), the zero of
\(G=Z+\lambda Z_h\) issuing from \(h_0\) continues and satisfies

\[
 |h-h_0|\le2|\Lambda|+3K|\Sigma|,\qquad |Q|\le5K.
\tag{Q9a}
\]

On a Robin zero the exact identities, with the partial derivatives of
\(Q=Z_{hh}/Z_h\) taken in the ambient \((h,\tau)\) variables, are

\[
 \frac{G_\tau}{Z_h}
 =\frac{Q+3\lambda\tau-3\lambda^2h}{2},\quad
 Q_h=3\tau-3\lambda h-Q^2,\quad
 Q_\tau=\frac{3h-3\lambda+3\lambda hQ}{2}.
\tag{Q9b}
\]

In the proposed region, \(|h|\le1.005K^3\), \(|\tau|\le1.001K^2\), and
\(|1+\lambda Q|\ge.995\). Thus \(|Q_h|<29K^2\),
\(|Q_\tau|<2K^3\), and implicit differentiation gives
\(|\dot h|<2|\Lambda|+3K|\Sigma|\). It follows that
\(|\dot Q|<58K^2|\Lambda|+89K^3|\Sigma|\), so
\(|Q-Q(0)|<.147K<1K\). This strictly closes the continuation bounds and
proves (Q9a). To check extension at a finite limiting parameter, bounded
\(Q\) rules out \(Z_h=0\): otherwise \(Z=Z_h=Z_{hh}=0\), and (Q2) would
make the entire function of \(h\) identically zero. That is impossible by
uniqueness of the Fourier transform of the nonzero integrable density in
(Q1), also for complex \(\tau\). Thus the implicit-function denominator
remains nonzero. This is a branch bound, not an exhaustive all-index theorem.

For example, on the *coupled* perturbation
\(\tau=\eta h\), \(\lambda=\eta\) near \(\tau_0=0\), the local
implicit-function theorem gives

\[
 h'(0)=-1-\tfrac12h_0Q_0,\qquad
 Z_h+\tfrac h2 Z_{hh}=\tfrac16 Z_{hhhhh}\quad(\tau=0).
\tag{Q9c}
\]

This odd quintic correction is not a pure real-temperature displacement.
In particular, allowing complex field along the actual finite-\(d\) Gamma
unfolding generally makes its temperature parameter complex as well. The
real-temperature zero theorem alone cannot locate that moving zero curve.

## 6. A two-sided compact perturbation strip

For every finite \(T,H>0\), there is a constant \(c(T,H)>0\) such that

\[
 |Z(h,\tau)|\ge c(T,H)|\Re h|
       (|Z|+|Z_h|+|Z_{hh}|),\quad
 |\tau|\le T,\quad |h|\le H.
\tag{Q10}
\]

At every zero, the implicit-function theorem gives
\(Z=(h-i x_j(\tau))g\), with \(g\) nonzero and continuous on a product
neighborhood. Since \(|h-i x_j(\tau)|\ge|\Re h|\), (Q10) follows there
after shrinking the neighborhood. The zero set in the displayed compact
parameter region is compact; finitely many such neighborhoods suffice.
On their complement, \(Z\ne0\); compactness and the nonvanishing of the
joint norm finish the estimate.

Consequently every zero in this compact region of any pointwise perturbation
\(P\) satisfying
\(|P-Z|\le\epsilon(|Z|+|Z_h|+|Z_{hh}|)\) obeys
\(|\Re h|\le\epsilon/c(T,H)\). This converse covers **all** perturbed
zeros in that compact window. The constant has not been evaluated.

A positive absolute derivative lower bound uniform in the zero index is
impossible even at a fixed \(\tau\): \(F_x(x,\tau)\to0\) as
\(|x|\to\infty\) by the Riemann–Lebesgue lemma, and the model has
infinitely many zeros. This does not rule out a useful *relative* vector-norm
estimate after rescaling high zeros.

## 7. Remaining Gamma-cusp bridge

The actual Gamma stationary cusp at \((d-1)^2=4a\) has a real quartic
maximum, but the following are not consequences of the model theorem:

1. Identification of its actual complete cycle with (Q1), including all
   nonlocal pieces and the chosen quartic branches.
2. A uniform quartic remainder in the joint norm
   \(|Z|+|Z_h|+|Z_{hh}|\), with parameter and order ranges stated.
3. A uniform comparison of the adjacent-order field **and temperature**
   parameters and amplitudes, strong enough to dominate the remainder and
   amplitude-induced zero motions. The finite-\(d\) mixed temperature jet
   is nonzero and must be included; it cannot be replaced by real \(\tau\)
   when the field is complex.
4. If the relevant scaled field window is unbounded, an exhaustive
   high-index Robin/vector perturbation theorem; (Q9) alone is a
   branch-existence statement.

Finite-dimensional coefficient convergence and the exact cubic obstruction
in `GAMMA_CURIE_WEISS_CRITICAL_BRIDGE.md` give no replacement for these
obligations. No assertion is made here about a new all-order Gamma family,
a uniform threshold near \(a=d^2/4\), or the full Strong Factorial Conjecture.

## Reproduction and primary reading

Run `python3 verify_quartic_lee_yang_zero_and_robin_bridge.py --output RECEIPT`.
The stdlib checker verifies the universal differential/flow algebra,
backwards-barrier and Robin-budget inequalities, critical CW coefficient
normalization, source pins, and deliberate sign/budget corruptions. It does
not numerically certify Newman's theorem or substitute finite tests for the
analytic proofs above. The source ledger gives exact inspected theorem
locators and PDF hashes. All review/replay assurance flags remain false.
