# Sharp two-mode energy and a finite-amplitude radial comparison

The first and second angular Fourier moments obey a sharp tradeoff at every
circular mass. This improves the factor \(4/3\) in AM7 to the optimal
factor \(1\). It gives a global two-mode matrix bound, including a
shared-radius version, and an explicit mass cube on which both moment
penalties can be retained. On that cube the penalties outweigh the
positive scalar contribution of the bounded radial profile from DG5 at
\(\rho=3/5\), uniformly over an explicit amplitude interval and all
centered label directions.

This is not the all-partition comparison at \(3/5\), or at every positive
correlation. The full direct inequality remains open. The analytic results
and the directed cube certificate have internal checking only; independent
replay, formal verification and external review are unperformed.

## TME1. The exact energy of the first two Fourier moments

For a measurable \(E\subset S^1\), use probability measure
\(d\sigma=d\theta/(2\pi)\), mass \(p=\sigma(E)\), and
\(z_m(E)=\int_E e^{im\theta}d\sigma\). Put
\(a(p)=\sin(\pi p)/\pi\). For every \(r,s\ge0\),

\[
 \max_{\sigma(E)=p}\{r|z_1(E)|^2+s|z_2(E)|^2\}
   =a(p)^2\max\{r+s\cos^2(\pi p),\ s\}.
 \tag{TME1}
\]

The first entry is attained by one arc of mass \(p\), and the second
by two equal antipodal arcs, each of mass \(p/2\). In particular

\[
 |z_2(E)|^2-a(p)^2\cos^2(\pi p)
 \le\sin^2(\pi p)\,[a(p)^2-|z_1(E)|^2].
 \tag{TME2}
\]

For every \(0<p<1\), the coefficient \(\sin^2(\pi p)\) is best
possible, by the two equal antipodal arcs. A single arc is an energy
maximizer exactly when \(r\ge s\sin^2(\pi p)\), apart from the
trivial zero-energy cases. This is a one-set theorem; applying its two
harmonic maxima separately without the tradeoff would be invalid.

Here is a complete reduction. Complementation negates both nonconstant
moments and leaves their squared magnitudes unchanged, so assume
\(0<p\le1/2\). The endpoint masses are trivial. Allow densities
\(0\le f\le1\), \(\int f=p\). This set is weak-star compact and the
energy depends continuously on four moments, so it has a maximizer.
If at least one weight is positive its maximum is positive: choose an
arc for the first harmonic or two antipodal equal arcs for the second.
First variation along every feasible density makes the maximizer a
linear maximizer for the score
\(r\operatorname{Re}(\overline{z_1}e^{i\theta})+
 s\operatorname{Re}(\overline{z_2}e^{2i\theta})\).
This is nonconstant, since the energy is positive. Every level set
has measure zero. Selecting its upper level set of mass \(p\) is
the unique linear maximizer almost everywhere: subtract the threshold
and compare the signs on the selected set and its complement.
Thus the maximizing density is ordinary and consists of at most two
arcs, since a nonconstant trigonometric polynomial of degree two has
at most four roots at any level.

Let their masses be \(p_1,p_2\), with sum \(p\), and let their smaller
center separation be \(\delta\in[\pi p,\pi]\). Set
\(A_j=\sin(\pi p_j)\), \(C_j=A_j\cos(\pi p_j)\), and
\(u=\cos\delta\). Both \(C_j\) are nonnegative because \(p\le1/2\).
Multiplying the energy by \(\pi^2\) gives

\[
 r(A_1^2+A_2^2+2A_1A_2u)
 +s(C_1^2+C_2^2-2C_1C_2+4C_1C_2u^2).
\]

This is convex in \(u\). Its maximum occurs when the arcs touch,
giving one arc, or when they are antipodal. In the antipodal case put
\(x=\pi p\), \(c=\cos x\), and
\(v=\cos(\pi(p_1-p_2))\in[c,1]\). The energy becomes

\[
             r(1+c)(1-v)+s(1-c^2)v^2.
\]

It is again convex. The endpoint \(v=c\) gives the single arc,
and \(v=1\) gives the equal antipodal pair, proving (TME1).
Take \(r=\sin^2(\pi p)\), \(s=1\) to obtain (TME2).
The argument is self-contained; no historical-priority claim is made.

## TME2. Global pair and shared-radius matrix certificates

For a circular three-cell partition let \(u_{m,i}\in\mathbb R^2\)
be its real harmonic moment in degree \(m=1,2\). Each moment triple
closes. On \(H=\mathbf1^\perp\) take arbitrary positive definite
matrices \(X_1,X_2\), extended by zero on \(\mathbf1\), and use
\(\lambda_i(X)=2X_{ii}+X_{jk}\) for distinct indices. Define the
explicit support

\[
 \Phi_p(r,s)=a(p)^2\max\{r+s\cos^2(\pi p),s\},\qquad r,s\ge0.
\]

For two circular partitions of masses \(p,q\), and \(c_1,c_2\ge0\),

\[
\begin{split}
 \sum_{m=1}^2c_m\sum_i u_{m,i}\cdot v_{m,i}
 \le\frac12\sum_i\big[&
 \Phi_{p_i}((c_1\lambda_i(X_1))_+,(c_2\lambda_i(X_2))_+)\\
 &+\Phi_{q_i}((c_1\lambda_i(X_1^{-1}))_+,
                         (c_2\lambda_i(X_2^{-1}))_+)\big].
\end{split}
 \tag{TME3}
\]

This is valid for every mass vector, including the boundary, and can
be minimized over the two matrices. It follows from the squared
matrix norm inequality separately in the two modes. The resulting
traces are \(\sum_i c_m\lambda_i(X_m)|u_{m,i}|^2\); discard negative
coefficients and apply the exact one-set support (TME1). There is no
optimization over measurable angular sets in the displayed upper bound.

A shared-radius version also follows. Write \(d\nu(t)=e^{-t}dt\)
and take any nonnegative integrable weight \(w\). With prescribed
mass profile \(p(t)\), choose measurable positive definite fields
\(X_1(t),X_2(t)\) and put

\[
 C_X=\int w\sum_i\Phi_{p_i(t)}
       ((c_1\lambda_i(X_1(t)))_+,(c_2\lambda_i(X_2(t)))_+)d\nu,
 \quad D_m=\int wX_m^{-1}d\nu.
\]

When these integrals are finite,

\[
 \sup_A\sum_{m=1}^2c_m\sum_i\left|\int w u_{m,i}\,d\nu\right|^2
        \le C_X\max_{m=1,2}\lambda_{\max}(D_m).
 \tag{TME4}
\]

The supremum uses common ordinary slices with their prescribed masses.
Dualize the joint norm of the two integrated moment matrices. Its
unit dual Grams \(B_1,B_2\succeq0\) satisfy
\(\operatorname{tr}(B_1+B_2)=1\). Matrix Cauchy--Schwarz and (TME1)
bound twice the dual pairing by
\(C_X+\sum_m\operatorname{tr}(D_mB_m)\), at most
\(C_X+\max_m\lambda_{\max}(D_m)\).
Scale both fields by the same positive scalar and minimize to obtain
(TME4), including zero cases by limits. The matrices are dual
certificates; the bound does not assert joint realizability of arbitrary
independently chosen harmonic vectors.

For \(P_2(z)=(1+z)^2=3/2+2\cos\theta+\tfrac12\cos2\theta\), set
\(a_{m,i}(t)=\sin(m\pi p_i(t))/(m\pi)\), allowing the signed
second coefficient. Define the common-slice polynomial loss

\[
 \mathfrak D^{(2)}_w(p)=
 2\sum_i\left(\int w a_{1,i}d\nu\right)^2
 +\tfrac12\sum_i\left(\int w a_{2,i}d\nu\right)^2
 -\sup_A\left[2\sum_i\left|\int w u_{1,i}d\nu\right|^2
       +\tfrac12\sum_i\left|\int w u_{2,i}d\nu\right|^2\right].
 \tag{TME5}
\]

It is nonnegative: circular rearrangement applies to the nonnegative
kernel \((1+\cos\theta)^2\), which decreases with circular distance.
The constant Fourier term cancels because masses are fixed. Taking
\(c_1=2,c_2=1/2\) in (TME4), subtracting its upper bound from the
two cap terms, and taking the positive part gives a global explicit
lower bound for (TME5) at every amplitude. The second harmonic must
remain coupled to the first in this argument; its individual arc
moment is not a global norm cap.

## TME3. Adding the second loss to the global radial comparison

Use \(d=1-\rho^2\), \(v=2\rho\sqrt{ts}/d\),
\(Q=d^{-1}e^{-\rho^2(t+s)/d}\), and
\(\eta=\rho/(1-\rho)\). Put

\[
 h_1(t)=\sqrt t e^{-\eta t},\quad h_2(t)=t e^{-\eta t},\qquad
 b_1=2\rho/d^2,\quad b_2=2\rho^2/d^3.
\]

The exact expansion
\(e^{v\cos\theta}=e^{-v}\sum_{k\ge0}
 v^k(1+\cos\theta)^k/k!\) has nonnegative monotone circular
kernels term by term. Consequently each term has a nonnegative
rearrangement deficit. The pointwise inequalities

\[
 Qe^{-v}v\ge b_1h_1(t)h_1(s),\qquad
 Qe^{-v}v^2/2\ge b_2h_2(t)h_2(s)
 \tag{TME6}
\]

follow from \(2\sqrt{ts}\le t+s\). They are scalar kernel
inequalities applied to nonnegative deficits, not operator inequalities.
Writing \(U_\rho(p)=\sum_i q_\rho(p_i-1/3)\), and keeping the
notation \(V_\rho\), \(\mathfrak D_{h_1}\) and any valid first-mode
pair penalty \(\gamma\) from GAC1--GAC7, gives

\[
\begin{split}
 S_\rho(A)-S_\rho(T)\le U_\rho(p)
 &-\iint(V_\rho-b_1h_1\otimes h_1)\gamma(p(t),p(s))d\nu^2\\
 &-b_1\mathfrak D_{h_1}(p)-b_2\mathfrak D^{(2)}_{h_2}(p).
\end{split}
 \tag{TME7}
\]

This holds for every ordinary balanced planar partition and every
\(0<\rho<1\). It strengthens (GAC7) by an additional nonnegative
shared-radius polynomial loss. The two losses come from distinct Taylor
terms, so they are not counted twice. It is not valid to append this
loss to an arbitrary larger cosine extraction without checking the
remainder again.

The scalar envelope in this step is the inspected all-parameter
Corollary A.3 of
[Heilman, arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1),
as used in DG6. The rearrangement input is
[Friedberg--Luttinger, Theorem I, printed pp. 1--2, (2)--(5)](https://digital.library.unt.edu/ark:/67531/metadc1025197/m2/1/high_res_d/4446751.pdf),
with the hypotheses checked in DG6 and GAC1. Integrability follows
from the full nonnegative exponential kernel and its radial probability
integral. Nonnegative summation justifies extracting the first two
terms. No stationarity or regularity of the partition is required.

## TME4. A certified cube retaining both penalties

Suppose \(|p_i-1/3|,|q_i-1/3|\le1/20\).
Let \(K_m(p)\) be the Gram on \(H\) of consecutive-arc moments
in harmonic \(m\), \(m=1,2\), and let
\(F(A,B)=\operatorname{tr}[(A^{1/2}BA^{1/2})^{1/2}]\).
Write \(\Gamma_m(p,q)=\sum_i a_{m,i}(p)a_{m,i}(q)
                              -F(K_m(p),K_m(q))\).
The directed full cube certificate proves

\[
\begin{split}
 c_1\sum_i u_{1,i}\cdot v_{1,i}
  +c_2\sum_i u_{2,i}\cdot v_{2,i}
 &\le c_1F(K_1(p),K_1(q))+c_2F(K_2(p),K_2(q)),
                    &&c_1\ge4c_2\ge0,\\
 \Gamma_1(p,q)&\ge\tfrac1{10}|p-q|^2,
 &\Gamma_2(p,q)&\ge\tfrac3{40}|p-q|^2.
\end{split}
 \tag{TME8}
\]

In particular the \(P_2\) pair deficit is at least
\(19|p-q|^2/80\). This is an explicit finite-amplitude version of
the local two-mode statement AM9, with a stronger one-set tradeoff.

Here is what is certified. Both Grams are positive definite on the
cube. For each mode set
\(X_m=(K_m(q)+\sqrt{\det K_m(p)\det K_m(q)}K_m(p)^{-1})/F_m\).
It satisfies \(X_mK_m(p)X_m=K_m(q)\). The six coefficients of
\(X_2,X_2^{-1}\) are positive, as are those for mode one, and

\[
 \lambda_i(X_1)\ge\tfrac14\sin^2(\pi p_i)\lambda_i(X_2),\qquad
 \lambda_i(X_1^{-1})\ge\tfrac14\sin^2(\pi q_i)\lambda_i(X_2^{-1}).
 \tag{TME9}
\]

In the two trace bounds the second-mode excess is bounded by (TME2)
times the first-mode cap deficit. Equations (TME9) absorb every such
excess, proving the first line of (TME8).

For the coercivities, the cap-Gram fidelity is attained by a rotation
of two consecutive-arc moment triples of the same orientation. Its
first-mode deficit satisfies

\[
 \Gamma_1(p,q)\ge
 \sum_{i<j}
 \frac{w_iw_j}{\sum_\ell w_\ell}
                    (1-\cos(\pi(p_k-q_k))),\qquad
 w_i=\frac{\sin(\pi p_i)\sin(\pi q_i)}{\pi^2},
 \tag{TME10}
\]

where \(k\) is the remaining label. This follows by
subtracting the magnitude of the weighted complex sum from its total
weight, as in AM5. On the cube replace the cosine quotient by
\(\pi^2/2-\pi^4/2400\), compress the resulting diagonal quadratic
form to \(H\), and certify it is at least \(I/10\).
The folded masses \(\widetilde p_i=1-2p_i\) lie in
\([7/30,13/30]\), and
\(K_2(p)=K_1(\widetilde p)/4\).
Apply (TME10) to them with cosine quotient lower bound
\(\pi^2/2-\pi^4/600\). The factors \(1/4\) in the Grams and
\(4\) in \(|\widetilde p-\widetilde q|^2\) cancel. The second
compressed form is at least \(3I/40\).

The [checker](check_three_cell_sharp_two_mode.py) exhausts the exact
four-dimensional root box for \(p_1,p_2,q_1,q_2\) in
\([17/60,23/60]\), clipping the complementary third coordinates to
the same range. It bisects a rational coordinate whenever positivity
of a Gram, a coefficient, an absorption margin, or either compressed
matrix minor has not yet been certified. Only boxes with impossible
third-coordinate ranges are discarded. The complete run used 100-bit
Arb intervals: 66241 attempts, 31528 accepted boxes and 1593 infeasible
boxes. The full binary-tree count is checked. All matrix and sine
expressions enclose entire boxes; no sampled mass pairs prove (TME8).
Rejected controls test excessive second weight, excessive coercivity,
and an incomplete cover.

## TME5. The penalties pay for an adverse radial profile

Put
\(f(t)=e^{-t}+2e^{-16t}-21/34\), so
\(\int f\,d\nu=0\), \(|f|\le81/34\).
For every centered unit vector \(v\in\mathbb R^3\), every
\(|\epsilon|\le1/40\), and every ordinary planar partition with
circular masses

\[
       p_i(t)=1/3+\epsilon v_i f(t),
\]

one has, at \(\rho=3/5\),

\[
                   S_{3/5}(T)-S_{3/5}(A)\ge\epsilon^2/500.
 \tag{TME11}
\]

The profiles are balanced and feasible, and
\(\max_i|p_i-1/3|\le(81/1360)\sqrt{2/3}<1/20\).
The statement allows arbitrary angular shapes, not just arcs.

The scalar contribution is exactly \(\epsilon^2q_{3/5}(f)>0\).
The checker replays the full directed angular integration from DG5--DG8
at 120-bit precision, with 8192 covering intervals, obtaining
\(q_{3/5}(f)\in[0.0096962141463,0.0101457669783]\), rounded outward
here. The first line of the Taylor decomposition, using
\(\Gamma_1\ge|p-q|^2/10\), costs more than
\(0.0092917675494\epsilon^2\). Its coefficient is
\((1/10)\iint V_\rho(t,s)(f(t)-f(s))^2d\nu^2\), equal to
\(8/5\) times the elementary DG31 correction replayed by the checker.
The second line uses
\(\Delta_{P_2}\ge19|p-q|^2/80\) and (TME6).
Its radial integrals are exact finite sums: with
\(\beta=1/(1-\rho)=5/2\), exponents \((0,1,16)\), coefficients
\((-21/34,1,2)\), and \(h_2=t e^{-\eta t}\), put

\[
 W=\beta^{-2},\quad F_1=\sum_j\frac{c_j}{(\beta+\lambda_j)^2},\quad
 F_2=\sum_{j,k}\frac{c_jc_k}{(\beta+\lambda_j+\lambda_k)^2}.
\]

The additional cost is
\((19/40)(2\rho^2/d^3)(WF_2-F_1^2)\epsilon^2
 >0.0032223371417\epsilon^2\).
Thus the full upper coefficient is below \(-0.0023683<-1/500\).
The entire amplitude interval follows from the cube theorem and exact
quadratic dependence, not a finite list of amplitude tests.

The [verification record](three_cell_sharp_two_mode_verification.json)
retains the full cube counts and margins, scalar replay, exact inputs,
source pins and assurance limits. Run from the canonical repository root:

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_sharp_two_mode.py --output .scratch/sharp-two-mode-production.json
```

## TME6. Why shared first moments alone do not close this route

The [equal-pair analysis](THREE_CELL_EQUAL_PAIR_ANGULAR_MOMENTS.md),
EP1--EP4, shows that for profiles
\(p=(x,(1-x)/2,(1-x)/2)\) sufficiently close to equal masses,
one common family of consecutive arcs attains every pairwise optimal
first-moment interaction simultaneously. This holds for every
nonnegative integrable radial interaction weight. Consequently replacing
the exact two-slice first-moment penalty by the exact shared first-moment
optimization gives no further improvement on these profiles. This
comparison uses the strengthened pair penalty in DG29, rather than
its earlier elementary penalty.

Take the bounded mean-zero six-exponential function \(g\) following
DG40 and the unit label direction \((2,-1,-1)/\sqrt6\).
For sufficiently small nonzero \(\epsilon\), these profiles are of
the preceding equal-pair type. The first-moment optimum is the full-cap
fidelity in this neighborhood, and EP4 supplies its common realizer.
It applies also to the positive mass-dependent scalar weight
\(R_0c_a^*(p_t,p_s)\) in DG39. The quadratic coefficient of that
global first-moment upper bound is therefore still
\(E_{99/100}(g)>1.20298\cdot10^{-6}\), certified in DG9.
The label direction has unit norm; the same coefficient follows by
the per-label second derivatives in DG40. The positive coefficient
persists for all sufficiently small nonzero amplitudes.

Thus fully shared first moments cannot establish the all-correlation
sign **within those separated remainder envelopes**. This is not a
counterexample to actual simplex optimality, nor a claim about every
conceivable use of first moments. A stronger remainder comparison or
additional angular information is necessary in this route. The sharp
second-mode support and (TME7) provide such additional information, but
their domination of the complete finite adverse radial space at every
feasible amplitude remains open.
