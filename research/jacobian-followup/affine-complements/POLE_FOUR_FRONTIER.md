# Four sections describe the next quartic pole pair

14 November 2026. PF2RULED3 gives a necessary normal form for the first
remaining pair (3,4) after [PF2RULED2](POLE_REDUCTION.md), and excludes one
explicit family of cubic leading data. The subsequent
[complete marked-atlas proof](POLE_FOUR_EXCLUSION.md), PF2RULED4, excludes
the entire (3,4) layer. This note retains the necessary normal form and
the first control certificate.

Use S=F1, B=C0+2f, L=O_B(B)=O_P1(3), V_m=H0(S,O(mB)),
and the required volume divisor E-3B. This argument retains the whole
allowed contact-four boundary family; the control in the last section
examines one specified family inside its coefficient problem.

## From two high-pole functions to four low-pole functions

Normalize the leading symbols of P,Q to tau^3,tau^4 for some nonzero
tau in H0(L). This follows from their leading differential equation:
q^3 is a constant times p^4, so div(p)=3D and div(q)=4D for a degree-three divisor D. Target scalings absorb the constants.

Choose T in V1 with symbol tau and write P=T^3+A2, A2 in V2. Then

\[
 D_3=Q-T^4-\tfrac43 T A_2\in V_3,
\qquad
 dP\wedge dQ=-\tfrac43 A_2\,dT\wedge dA_2+dP\wedge dD_3.
\]

The first term has pole at most six. Vanishing of the pole-seven term in
the second forces the leading symbol of D3 to be a scalar times tau^3.
An affine target change Q -> Q-cP therefore puts D3 in V2; call it D2.

Let a2,d2 in H0(L^2) be the symbols of A2,D2, allowing zero symbols.
The remaining differential is

\[
 -\tfrac43 A_2\,dT\wedge dA_2+
 3T^2\,dT\wedge dD_2+dA_2\wedge dD_2.
\]

Its pole-six coefficient is

\[
 -\tfrac83\tau'a_2^2+\tfrac43\tau a_2a_2'
 +6\tau^2\tau'd_2-3\tau^3d_2'
 =
 -3\tau^5\left(
 \frac{d_2}{\tau^2}-\frac{2a_2^2}{9\tau^4}
 \right)'.
\]

It follows that

\[
 d_2=\frac{2a_2^2}{9\tau^2}+\kappa\tau^2.
\]

All symbols are regular sections. At each zero of tau, the last identity
forces tau^2 to divide a2^2, hence tau divides a2. Write a2=tau*alpha
for alpha in H0(L), and lift alpha to A in V1. Then

\[
 B_1=A_2-TA\in V_1,\qquad
 C=D_2-\tfrac29 A^2-\kappa T^2\in V_1.
\]

Renaming B1 as H avoids confusing this function with the infinity divisor.
Every candidate therefore has, after the target changes above,

\[
 \boxed{
 \begin{aligned}
 P&=T^3+TA+H,\\
 Q&=T^4+\tfrac43T^2A+\tfrac43TH+\tfrac29 A^2+\kappa T^2+C,
 \qquad T,A,H,C\in V_1,\quad\kappa\in\mathbb C.
 \end{aligned}}
 \tag{1}
\]

Since dim V1=5, (1) has21 scalar coefficients before further equations,
nondegeneracy conditions and quotienting redundant descriptions. This is
a necessary parameterization, not a dimension assertion for a moduli space.

## The next leading relation is cubic

In the independent variables T,A,H,C, the coefficients of dP wedge dQ are

\[
\begin{array}{c|l}
 dT\wedge dA & \tfrac49 A^2-\tfrac43TH-2\kappa T^2\\
 dT\wedge dH & -\tfrac43TA-\tfrac43H-2\kappa T\\
 dT\wedge dC & 3T^2+A\\
 dA\wedge dH & -\tfrac49A\\
 dA\wedge dC & T\\
 dH\wedge dC & 1 .
\end{array}
\]

Write tau,alpha,eta,gamma in H0(L) for their leading symbols. Vanishing
of the pole-five term says the derivative of

\[
 3\gamma/\tau-\tfrac43\alpha\eta/\tau^2+
 \tfrac4{27}\alpha^3/\tau^3-2\kappa\alpha/\tau
\]

is zero. Thus, for a scalar lambda,

\[
 \boxed{
 4\alpha^3-36\tau\alpha\eta+81\tau^2\gamma
 -54\kappa\tau^2\alpha=\lambda\tau^3.}
 \tag{2}
\]

Pole four must still vanish. The remaining pole-three coefficient must
define the required smooth section E, with its contact-four intersection
with B. Generic degree four and finiteness must then be established.

The odd-degree conic argument from pole three does not carry over directly.
For example, with kappa=lambda=0 the binary cubics

\[
 \tau=s^3,\quad\alpha=3s^2t,\quad
 \eta=3st^2,\quad\gamma=\tfrac83t^3
\]

satisfy (2) and have no common zero. Thus an obstruction must inspect the
next coefficient or use additional global geometry.

## An exact obstruction to one cubic-symbol family

Normalize S as the blowup of [0:0:1] in P2, and B as the strict transform
of XZ-Y^2=0. In the chart Z=1 set b=x-y^2 and put

\[
 g_1=y/b,\quad g_2=y^2/b,\quad
 g_3=xy/b,\quad g_4=x^2/b.
\]

Together with1 these form a basis of V1. Their global symbols, in the
natural trivialization along the conic, are1,w,w^2,w^3. Consider all
constant lifts and the indicated scalar freedoms of the displayed cubic
data:

\[
 \begin{aligned}
 T&=g_4+t_0,& A&=3g_3+a_0,\\
 H&=3g_2+h_0,&
 C&=\tfrac83g_1+\tfrac23\kappa A+\ell T+c_0.
 \end{aligned}
\]

Substitute these into (1), with arbitrary
t0,a0,h0,c0,kappa,ell in C. All satisfy the pole-five relation (2).
In the regular local coordinates z=b,w=y at the generic point of B,
the coefficient of z^(-4) dz wedge dw is

\[
 \boxed{
 w^7\big[
 (-4h_0-3\ell+8)w^3-6\kappa w^2+12t_0w+4
 \big].}
 \tag{3}
\]

Its w^7 coefficient is4, independently of every lift parameter.
Consequently this entire six-parameter control family cannot meet the
required pole-three differential condition. The retained exact
[calculation](explore_pole4_cubic_control.py) and
[coefficient certificate](pole4_cubic_control.json) establish this small
exclusion.

This does not classify all cubic solutions of (2), including changes in
their root partitions or their position relative to the marked point of
the conic. No such equivalence is assumed.

## The next finite calculation

Work with the21 coefficients in (1), impose (2), and then expand the
pole-four coefficient. Use the common leading cubic and its root
partitions to organize the equations before elimination. The
[714-column linear kernel](pole_frontier.json) remains a complete
independent calibration of the unreduced necessary pole equations.

The [complete exclusion](POLE_FOUR_EXCLUSION.md) now covers every
leading-data family, using exact pole and contact contradictions.
The [active research goal](RESEARCH_GOAL.md) has moved to height-five
pairs and the actual-normalization compatibility problem.

