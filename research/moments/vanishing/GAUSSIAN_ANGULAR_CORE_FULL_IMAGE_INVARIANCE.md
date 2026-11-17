# Angular cores preserve the full actual marked image

Let a characteristic-zero moment functional kill every nonzero angular
weight. Adding arbitrary terms of strictly positive angular weight to a
weight-zero phase preserves its **embedded full actual marked module**.
No safety, purity, finite-generation or fibre-degree hypothesis is needed.
The original fixed-marker image can enlarge, even for a nonconstant core.
Finite polynomial-time amplitude transports give its exact kernel.

## 1. Statement and original Gaussian and sphere functionals

Let (A=\bigoplus_{q\in\mathbb Z}A_q) be a commutative graded algebra over
a characteristic-zero field (K), and let (E:A\to K) be linear with

\[
 E(A_q)=0\quad(q\ne0).
\]

Every element has finite homogeneous support. For (P,Q\in A), define

\[
 \Phi_P[Q]=\sum_{m\ge0} E(QP^m)\frac{t^m}{m!},\qquad
 V_P=\Phi_P[A],\qquad M_P=K[t]V_P\subset K[[t]].       \tag{A1}
\]

The (K[t])-linear extension on (A[t]) has image (M_P). A phase is
safe when every original period is polynomial, and pure when its vacuum
is one. These definitions use the normalized functional (E(1)=1).

Suppose

\[
 P=P_0+R,\qquad P_0\in A_0,\qquad
 R\in\bigoplus_{q\ge\gamma}A_q,\quad\gamma\ge1.       \tag{A2}
\]

Then, as actual submodules of (K[[t]]),

\[
 \boxed{M_P=M_{P_0}},\qquad
 V_{P_0}\subseteq V_P,\qquad
 \Phi_P[1]=\Phi_{P_0}[1].                             \tag{A3}
\]

Thus safety, purity and finite generation are individually equivalent
for (P) and (P_0). Their generic marked ranks also agree. Negative
one-sided support is covered by reversing weights. If (R=0) all
identities below use only the summand (n=0).

For the original Gaussian take Witt pairs (Z_i,W_i) of covariance

\[
 E(Z_iW_j)=\delta_{ij},\quad E(Z_iZ_j)=E(W_iW_j)=0,
 \qquad \operatorname{wt}Z_i=1,\quad
 \operatorname{wt}W_i=-1,
\]

and independent ordinary variables (Y_j) of weight zero. Wick
contraction kills nonzero total angular weight. With several pairs the
core belongs to the whole invariant algebra (A_0); it need not be a
polynomial in (\sum_iZ_iW_i). Other integral weights preserving the
covariance are allowed. The same theorem applies to the original sphere
functional: rotation of each complex Witt pair preserves that functional
and kills nonzero angular weight. The proof uses only this invariance.

## 2. Two finite identities prove actual equality

Write (q_{\min}(Q)=\min\{q:Q_q\ne0\}) for nonzero (Q), and set

\[
 N_Q=\max\left(0,\left\lfloor
                    \frac{-q_{\min}(Q)}\gamma\right\rfloor\right).
\]

For (n>N_Q), every weight of (QR^n) is positive. Since (P_0) and
(P) have nonnegative weights, both corresponding periods vanish. Hence
the formal exponential identities yield the **finite** formulas

\[
 \begin{aligned}
 \Phi_P[Q]&=\sum_{n=0}^{N_Q}\frac{t^n}{n!}
                              \Phi_{P_0}[QR^n],\\
 \Phi_{P_0}[Q]&=\sum_{n=0}^{N_Q}\frac{(-t)^n}{n!}
                              \Phi_P[QR^n].
 \end{aligned}                                      \tag{A4}
\]

For example, multiply (e^{-tR}e^{tP}=e^{tP_0}) by (Q) and apply
(E) coefficientwise. Every coefficient is a finite algebraic sum; the
outer terms with (n>N_Q) vanish identically. No convergent integral or
infinite interchange is involved. Formula (A4) holds for amplitudes in
(A[t]) as well. Its two directions prove the module equality in (A3).

Let (Q_0) denote the weight-zero projection. Then

\[
 \Phi_{P_0}[Q]=\Phi_{P_0}[Q_0]=\Phi_P[Q_0].          \tag{A5}
\]

Every term containing (R) in the last expression has positive weight.
This proves the literal inclusion at the same original marker degree
whenever degree respects homogeneous projection. It also proves the
vacuum identity. Safety means (M_P\subseteq K[t]), so its equivalence
follows from actual equality, without any finite-image theorem.

## 3. Exact kernel transport and the time connection

Let (A_+=\bigoplus_{q>0}A_q), and form the vector-space quotient

\[
 \mathcal B=A[t]/A_+[t],\qquad \pi:A[t]\to\mathcal B.
\]

This is a quotient by a **linear subspace**. It is generally not a
quotient algebra: (Z\in A_+) but (ZW\notin A_+). Both period maps
kill (A_+[t]). Multiplication by (P_0) and (R) preserves that
subspace, so it induces commuting endomorphisms on \(\mathcal B\).

Multiplication (N) by (R) is locally nilpotent there: after enough
applications every weight of a given input is positive. Therefore

\[
 T_+=e^{tN},\qquad T_-=e^{-tN}                       \tag{A6}
\]

are mutually inverse (K[t])-linear automorphisms. Each exponential
is a finite sum on each input, with no uniform marker-independent bound.
The descended maps satisfy

\[
 \overline\Phi_P=\overline\Phi_{P_0}T_+,\qquad
 \ker\overline\Phi_P=T_-\ker\overline\Phi_{P_0}.     \tag{A7}
\]

In particular the exact kernel on original fixed spatial markers is

\[
 \ker(\Phi_P|_A)=
 \{Q\in A:T_+\pi(Q)\in\ker\overline\Phi_{P_0}\}.    \tag{A8}
\]

This is an exact transport to polynomial-time amplitudes. It does not
give a time-independent spatial-marker automorphism or an ideal kernel.

The time operators on \(\mathcal B\) are

\[
 D_P=\partial_t+P_0+R,\qquad D_0=\partial_t+P_0.
\]

Because (\partial_tT_+=NT_+) and the multipliers commute,

\[
 D_0T_+=T_+D_P.                                    \tag{A9}
\]

Both period maps are horizontal for these time operators, and
(T_+\pi(1)=\pi(1)). Spatial derivatives need not preserve (A_+);
no isomorphism of full spatial Stein quotients is asserted.

## 4. Original degree filtration and finite-jet bounds

For standard angular weights, \(\deg Q\le D\) implies
(q_{\min}(Q)\ge-D). Put (s=\deg R), and let (V_{P,D}) be the
image of original markers of degree at most (D). Then

\[
 V_{P_0,D}\subseteq V_{P,D}\subseteq
 \sum_{n=0}^{\lfloor D/\gamma\rfloor}
                         t^nV_{P_0,D+ns}.          \tag{A10}
\]

For a polynomial-time input of time degree (e) and spatial degree
(D), either transport has a representative of time degree at most
(e+\lfloor D/\gamma\rfloor) and spatial degree at most
(D+s\lfloor D/\gamma\rfloor). For generator weights of absolute
value at most \(\Lambda\), replace (D) inside the floor by
(\Lambda D), or retain the exact (N_Q).

Suppose a common scalar exponent \(\psi\) satisfies

\[
 e^{-\psi}\Phi_{P_0}[Q]\in K[t],\qquad
 \deg_t(e^{-\psi}\Phi_{P_0}[Q])\le C_0\deg Q+E_0.
\]

The finite formula (A4) gives the original-marker bound

\[
 \deg_t(e^{-\psi}\Phi_P[Q])\le
 C_0D+E_0+(1+C_0s)\lfloor D/\gamma\rfloor.           \tag{A11}
\]

Here the degree of zero is \(-\infty\). If the core is constant,
(C_0=E_0=0), so the normalized cutoff is simply
(\lfloor D/\gamma\rfloor), independently of \(\deg R\).
If that constant is zero, the original moments are pure and safe.
The [mixed-cutoff finite-jet theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md)
therefore applies with \(C=1/\gamma\): degree-at-most-\(d\)
perturbations through order \(r\) require moments only through
\(\lfloor(1+d/\gamma)r\rfloor\). This is a sufficient bound;
the least ideal-generating prefix can be smaller.

## 5. Literal images can strictly enlarge

First, (P_0=0,R=Z) gives

\[
 V_{P_0}=K,\quad V_P=K[t],\quad M_{P_0}=M_P=K[t],
 \qquad \Phi_Z[W^j]=t^j.
\]

There is also a nonconstant-core example. Put

\[
 S=ZW,\quad P_0=S^2,\quad P=S^2+Z,\qquad
 f=\Phi_{P_0}[1]=\sum_{n\ge0}\frac{(2n)!}{n!}t^n.
\]

Unbalanced original markers have zero core period, while

\[
 [t^n]\Phi_{P_0}[S^j]=\frac{(2n)!}{n!}
                         \prod_{i=1}^j(2n+i).      \tag{A12}
\]

Thus every core literal period has coefficients ((2n)!/n!) times a
polynomial in (n); ordinary independent Gaussian markers only multiply
these by constants. A polynomial period would make that polynomial in
(n) vanish for every sufficiently large integer. Consequently

\[
 V_{P_0}\cap K[t]=\{0\}.                            \tag{A13}
\]

Angular balance and the factorial identity give

\[
 \Phi_P[W]=t\Phi_{P_0}[S],\qquad
 f-1=2t\Phi_{P_0}[S],\qquad
 \Phi_P[1-2W]=1.                                   \tag{A14}
\]

Hence (V_{P_0}\subsetneq V_P), although the actual full modules agree.
The distinction concerns original spatial markers, not a localization.

For a nonconstant Gaussian core with finite full image, the
[retained literal-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
gives (V=M). In that case both nonconstant phases have the same literal
image too. The module in (A12)-(A14) is not finite; its vacuum has
factorially growing coefficients, whereas every finite Gaussian module
has an exponential-polynomial vacuum. The all-degree equality (A3)
does not use this imported theorem.

## 6. Repeated face reductions and remaining classification

For any integral covariance-preserving weight with phase support in the
closed nonnegative half-space, its weight-zero Newton face has exactly
the same full actual module and distinguished unit. A sequence of such
face reductions preserves them at every step. Formula (A7) composes the
finite amplitude transports, retaining an exact original kernel test.
This is an equality of actual embedded images, stronger than an abstract
rank or degeneration statement.

With one Witt pair, the remaining invariant core is an arbitrary
polynomial in (ZW,Y_1,\ldots,Y_s). Its classification remains a genuine
problem. The [radial constant-fibre theorem](GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md)
settles the subfamily (g(ZW)+B(ZW)Y+aY^2), with arbitrary one-sided
remainder. No general invariant-core or Gaussian safety-to-purity
classification is claimed here.

## 7. Evidence and scope

The proof of (A3)-(A11) is algebraic and applies to every invariant graded
functional described above. The literal finite-image observation alone
imports the retained Gaussian classification. The companion
`verify_gaussian_angular_core_full_image.py` checks bounded original Wick
moments, both finite transports, the time-connection identity and the
strict literal-image example using exact arithmetic. Its receipt records
the bounds and source hashes. These are internal checks, not independent
review, formal verification or an all-input computational proof.
