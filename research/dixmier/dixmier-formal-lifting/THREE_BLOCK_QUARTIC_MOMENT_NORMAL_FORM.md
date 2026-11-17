# Three canonical quartic blocks have no essential surviving mate direction

Let k have characteristic zero and eta,s in k^times. Work in
k(V)[A;eta V^{-1}partial_V], and put

\[
 H=V^2A,\qquad c(H)=H(H-2\eta),\qquad L=sV^{-3}A^2.
\]

Fix the ordered bottom coefficient c(H-7eta)c(H), and consider

\[
 \begin{split}
 P={}&s^2V^{-14}c(H-7\eta)c(H)
          +V^{-9}C_3(H)+V^{-4}C_2(H)+\mu V,\\
 C_3(H)&=a_3H^3+a_2H^2+a_1H+a_0,\\
 C_2(H)&=gH^2+e_1H+e_0,\qquad g\mu\ne0.
 \end{split}                                             \tag{1}
\]

The first term is exactly L^2. There is no finite nonzero scalar mate
of P over k(V) or any finite algebraic coefficient cover, in any
partner order. The classification below identifies every initial-moment
kernel direction as an actual square correction, a position potential,
or the already controlled first-order momentum residual.

This is a theorem about the fixed bottom polynomial and three indicated
canonical weights -14,-9,-4, with the position term of weight one.
It is not a classification of arbitrary quartic Laurent operators or
of all actual globally admitted Long images. Polynomial-centralizer
claims are restricted to the branches where a moment certificate proves
them; the final selfadjoint branch has only the no-mate conclusion here.

## The non-square top branch

Give V weight one and A weight three. The top symbol of (1) is

\[
           s^2V^{-6}\xi^4+a_3V^{-3}\xi^3+g\xi^2.          \tag{2}
\]

All lower coefficients have strictly lower weight. In particular its
A^1 coefficient has degree at most -2 and its scalar coefficient at
most one, within the strict bounds of the
[mixed quartic top-symbol theorem](MIXED_QUARTIC_TOP_SYMBOL_OBSTRUCTION.md).
If a_3^2-4s^2g!=0, that theorem proves rational Cent(P)=k[P] and
no finite scalar mate, including on finite algebraic coefficient covers.
Thus only the perfect-square relation

\[
                              a_3^2=4s^2g                \tag{3}
\]

needs further analysis. Since g!=0, its cubic coefficient is nonzero.

## The full cubic initial residue functional

Use the canonical weight wt(V)=1, wt(A)=-2. The middle block in (1)
has gap five from L^2; the upper block has gap ten and the position
term gap fifteen. Set theta=H/eta and c_0(x)=x(x-2), so
L=s eta^2 V^{-7}c_0(theta).
The first middle-block response of the formal square root is
V^{-2}F(theta), where

\[
 c_0(x-2)F(x)+c_0(x)F(x-7)=\frac{C_3(\eta x)}{s\eta^2}.
                                                                  \tag{4}
\]

Each coefficient is solved by division by two. Its expansion begins

\[
 F(x)=\frac{a_3\eta}{2s}x+
       \left(\frac{a_2}{2s}+\frac{15\eta a_3}{4s}\right)
                          +\ell_3x^{-1}+O(x^{-2}),
 \quad \ell_3=\frac{a_1+4\eta a_2+19\eta^2a_3}{2s\eta}.
                                                                  \tag{5}
\]

Matching the coefficients of x^3,x^2,x in (4) proves this formula;
no finite mate ansatz is involved.

If ell_3!=0, the entire branch is closed by the a=-2 moment certificate
from the mixed-top note, not only its first flow. For completeness,
write mathcal L(B)=[x^{-1}]F(x)B(x) and
p_n(x)=product_{j=1}^n c_0(x+7j). Since the right side of (4) is a
polynomial, formal translation-invariant residue gives

\[
 \mathcal L\{c_0(x-2)B(x)+c_0(x+7)B(x+7)\}=0.
\]

The alternating translated response terms of the positive half-power
P^{k+1/2} reduce to E_k=mathcal L(p_{2k}), E_0=ell_3.
The exact polynomial certificate is

\[
\begin{split}
 c_0(x+7(n+1))c_0(x+7(n+2))-r_n^-
   &=c_0(x-2)Q_n^-(x)+c_0(x+7(n+1))Q_n^-(x+7),\\
 Q_n^-(x)&=\frac{x^2}{2}+\frac{7(3n+4)}2x
                         +\frac{7(n+1)(35n+58)}4,\\
 r_n^-&=-\frac{7(n+1)(7n+9)(7n+11)(7n+16)}4.
 \end{split}                                             \tag{6}
\]

Multiplying by p_n and applying the residue functional gives
E_k=r_{2k-2}^- E_{k-1}!=0 for every k>=1.
Only the A^{-1} coefficient of the negative part contributes to A^2
in its commutator with P. It yields the nonzero most polar coefficient

\[
 [P,(P^{k+1/2})_+]_{A^2}
 =(56k+22)E_ks^{2k+2}\eta^{4k+2}V^{-12-14k}
                         +\text{terms of greater degree}.        \tag{7}
\]

Indeed the negative coefficient is
E_k*s^{2k}*eta^{4k+1}*V^{-4-14k}, and differentiation against
p_4=s^2V^{-6} gives (7). A second middle insertion, an upper insertion,
or a position insertion has strictly greater canonical weight.
All rational finite scalar-commutator candidates have even order by
the leading coefficient's V-divisor -6. Leading subtraction therefore
expresses them as a polynomial in P plus finitely many positive odd
powers of P^{1/2}. The largest-index nonzero term (7) cannot cancel.
This proves rational Cent(P)=k[P] and no finite scalar mate when
ell_3!=0. Coefficients descend from any constant square-root extension
by comparing highest A coefficients.

## The exact kernel and its square normal form

Now assume (3) and ell_3=0. Define

\[
 b=\frac{a_3}{2s},\qquad
 d=\frac{a_2}{2s}+\frac{15\eta a_3}{4s},\qquad
                         Q=L+bA+dV^{-2}.                 \tag{8}
\]

An exact finite Ore expansion gives

\[
\begin{split}
 \{L,bA+dV^{-2}\}
   =V^{-9}\bigl[&sbH(H-2\eta)(2H-11\eta)\\
                &+sd(2H^2-8\eta H+8\eta^2)\bigr].
\end{split}                                               \tag{9}
\]

Its cubic, quadratic and linear coefficients equal those of C_3,
because ell_3=0. The remaining constant is

\[
                         \kappa=a_0-4\eta^2a_2-30\eta^3a_3.
                                                                  \tag{10}
\]

The upper square correction is also exact:

\[
 (bA+dV^{-2})^2
       =V^{-4}(bH+d)(bH+d-2\eta b).                      \tag{11}
\]

Condition (3) says g=b^2. Hence after subtracting (11), the upper
residual has degree at most one in H. Put

\[
 \tau=e_1+2\eta b^2-2bd,\qquad
 \gamma=e_0-d^2+2\eta bd.
\]

The complete normal form is

\[
                 P=Q^2+\kappa V^{-9}+\tau V^{-2}A
                                      +\gamma V^{-4}+\mu V.       \tag{12}
\]

This proves the finite-dimensional quotient statement. The cubic
polynomial space has dimension four and ell_3 has one-dimensional
image. Its three-dimensional kernel contains the two-dimensional
anticommutator image (9); the remaining direction is the scalar
position term kappa*V^{-9}. After the upper quadratic coefficient is
fixed by (3), its two residual directions are the momentum term
tau*V^{-2}A and the position term gamma*V^{-4}. There is no additional
essential direction in these two blocks.

## The selfadjoint branch has no mate

Let D=eta V^{-1}partial_V. The coefficient-fixing antiinvolution

\[
 A^\dagger=-A+u(V),\qquad
 u(V)=\frac{D(sV^{-3})-b}{sV^{-3}}
                         =-3\eta V^{-2}-(b/s)V^3         \tag{13}
\]

fixes Q; its scalar coefficient dV^{-2} does not affect this formula.
If tau=0, (12) is selfadjoint. A finite scalar mate would yield a
nonzero skewadjoint mate, which necessarily has odd A-order.
But the leading coefficient of P is s^2V^{-6}. Highest commutator
comparison forces the leading coefficient of a mate of order M to
have V-degree -3M/2, which is integral only when M is even.
This contradiction excludes every finite rational mate.

No rational polynomial-centralizer theorem is inferred in this branch.
The argument only excludes scalar mates. In particular it does not by
itself justify parameter-pole removal for an arbitrary deformation with
this critical reduction.

## The nonselfadjoint upper residual retains the second certificate

Suppose tau!=0. Write

\[
 P_0=Q^2+\kappa V^{-9}+\gamma V^{-4}+\mu V,
                         \quad P=P_0+\tau V^{-2}A.
\]

The complete base P_0 is fixed by (13), and its formal square root
with leading term L is also selfadjoint by uniqueness. Its positive
odd root powers are therefore selfadjoint. Their commutators with
P_0 are skewadjoint and have order at most two. Thus their A^2
coefficients vanish exactly, at every canonical weight. This is
stronger than ignoring the more polar position term kappa*V^{-9}.

The tau insertion has gap ten from L^2. The difference Q-L and the
kappa position term have gap five, the gamma term gap ten, and mu*V
gap fifteen. Every contribution containing tau and another base
correction has gap strictly greater than ten; tau squared has gap
twenty. Hence the most polar A^2 anomaly is precisely the first tau
response around L. No tau-independent A^2 term survives by the exact
adjoint cancellation just proved.

That response and its moment recurrence are the a=3 certificate of
[the nonsymmetric perturbation theorem](TWISTED_ADJOINT_SQUARE_OBSTRUCTION.md):

\[
 R_\tau=\frac{\tau}{s\eta}V^3F_+(\theta),\qquad
 c_0(x+3)F_+(x)+c_0(x)F_+(x-7)=x,
\]

with F_+=(1/2)x^{-1}+O(x^{-2}). The translated residues give
Delta_0=1/2 and

\[
 \Delta_k=-2(7k-6)(7k-5)(7k-4)(14k-3)\Delta_{k-1}\ne0.
                                                                  \tag{14}
\]

This is the proved all-index recurrence, not a finite computation.
For B_{4k+2}=(P^{k+1/2})_+ the decisive coefficient is

\[
 [P,B_{4k+2}]_{A^2}
 =(56k+2)\Delta_k\tau s^{2k+1}\eta^{4k+1}V^{-7-14k}
                       +\text{terms of greater degree}.          \tag{15}
\]

Highest-index subtraction and noncancellation as in the ell_3!=0
branch give rational Cent(P)=k[P] and no scalar mate for tau!=0.
This proves that the previously known first-order residual remains
closed even with all square-completion and position terms in (12).

In every branch, finite algebraic coefficient-cover mates are excluded
by normalized coefficient trace before applying the rational argument.
The centralizer equalities in the two certificate branches are over
the rational coefficient field only.

## What an actual global application must verify

This theorem starts with (1), a critical operator over a generic constant
field that may contain the relative parameter r. A Laurent expression
alone does not prove polynomial admission in the original global
centralizer K[[h]][X,Q][B;h Delta]. A caller must supply an actual finite
global position, verify its exact commutation with R and prescribed
classical reduction, and justify the finite critical normalization and
quantum rescaling that produce (1).

If ell_3!=0 or if ell_3=0,tau!=0, the rational polynomial centralizer
proved above supplies QCL1's finite parameter-pole descent. Under
those caller hypotheses, a generic mate's finitely many parameter
poles are removed by subtracting finite polynomials in the integral
rescaled position; its reduction would contradict the critical
no-mate theorem. Algebraic-cover mates may first be traced down.

If ell_3=tau=0, the critical no-mate statement remains valid, but that
generic descent does not follow from this note alone. It requires an
independent polynomial-centralizer hypothesis, an already integral
rescaled mate, or a direct obstruction for the complete actual operator.
An exact full selfadjoint-square structure with the odd leading divisor
can give such a separate direct obstruction. No missing centralizer or
pole-descent premise is silently supplied here.

The bottom polynomial in (1) is fixed. Other resonant bottom patterns,
additional canonical blocks, and genuinely different lower momentum
structures are not classified. No new surviving global candidate,
bounded Long lift, DC1/DC2 solution, formal verification, independent
whole-proof replay, external review or novelty is claimed. The finite
identities (5),(9),(11) were checked by symbolic expansion; the
unrestricted partner-order conclusions use the written recurrences,
adjoint argument and finite normal-form exhaustion.
