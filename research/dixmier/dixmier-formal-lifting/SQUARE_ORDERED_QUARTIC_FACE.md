# The square-ordered quartic face has no finite scalar mate

Let k have characteristic zero, eta,g,s in k^times, and mu in k.
In the finite scalar differential algebra

\[
 k(V)[A;\eta V^{-1}\partial_V],\qquad
 L=sV^{-3}A^2,\qquad P=L^2+gA^2+\mu V,                   \tag{1}
\]

the rational-coefficient centralizer is exactly k[P], and no finite
operator over k(V), or any finite algebraic coefficient cover, has
nonzero scalar commutator with P. There is no partner-order bound.
The assumption g!=0 matters: if g=mu=0, L commutes with P=L^2 and
is not polynomial in P.

The proof gives an exact recurrence for the most polar commutator
coefficient of every positive half-integer power. Its nonvanishing
is established by one polynomial identity with a symbolic index.
The [ordered quartic infinity theorem QOI1](ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md)
provides a separate proof and covers additional scalar orderings;
that theorem is not a premise of the residue argument here.

## Every finite candidate has an even-flow normal form

Write a_4=s^2V^{-6}. If a finite rational T has positive order M
and [P,T] has order at most two, its highest commutator coefficient is

\[
 4a_4\delta(t_M)-M t_M\delta(a_4)=0,\qquad
       \frac{t_M'}{t_M}=-\frac{3M}{2V},
                         \quad\delta=V^{-1}\partial_V.  \tag{2}
\]

The order of a rational function at V=0 is integral, so M is even,
and t_M is a constant multiple of s^{M/2}V^{-3M/2}.
An order-zero remainder must be constant, since otherwise its
commutator has order three.

Let R=P^{1/2} be the unique formal pseudodifferential square root
with leading term L. Multiplication uses

\[
 A^n f=\sum_{j\ge0}{n\choose j}\eta^j\delta^j(f)A^{n-j}.
\]

Each individual coefficient is a finite sum. Successive coefficients
of R^2=P are solved by division by 2sV^{-3}, proving existence
and uniqueness. The coefficients are finite Laurent polynomials in V,
polynomial in g,mu, with coefficients rational in nonzero s and eta.
The formal root is a proof device, not an asserted finite operator.
For k>=0 put B_{4k+2}=(P^{k+1/2})_+.
Its commutator with P has order at most two, since the negative part
has order at most minus one. Leading subtraction using (2) gives
for every finite T with [P,T] of order at most two

\[
                 T=q(P)+\sum_{k=0}^N c_kB_{4k+2},
                             \quad q\in k[t],\ c_k\in k. \tag{3}
\]

The sum is finite and every lower constant flow is retained. At orders
divisible by four subtract a power of P; at orders two modulo four
subtract the corresponding B. Odd remainders are excluded by (2),
and the process strictly lowers differential order. There is no
restriction on the rational position degrees of T.

## The decisive coefficient and its residue sequence

Give V weight one and A weight minus two. Ore contractions preserve
this weight. The three parts L^2, gA^2 and mu*V have weights -14,-4,1.
In P^{k+1/2}, each insertion of g raises weight by ten relative to
its g=mu=0 part L^{2k+1}; each insertion of mu raises it by fifteen.
This follows coefficientwise during the square-root recursion and
multiplication, where division uses only sV^{-3}. Projection onto
nonnegative A-powers preserves each monomial's weight.

The zero-insertion commutator vanishes because L^{2k+1} is an actual
operator commuting with L^2. Thus the most polar possible part of
[P,B_{4k+2}] is its first-order-in-g part, of weight -14k-11.
All mu terms and higher powers of g have greater weight.
If the coefficient of A^{-1} in the g-linear part of P^{k+1/2} is

\[
       b_{-1}=D_k g s^{2k-1}\eta^{4k+1}V^{1-14k},
\]

then the A^2 coefficient of minus its commutator with P is

\[
 [P,B_{4k+2}]_{A^2}
 =(56k+2)D_k g s^{2k+1}\eta^{4k+2}V^{-7-14k}
                    +\text{terms of greater V degree}. \tag{4}
\]

Indeed only the A^{-1} term of the negative part can contribute to
A^2. If p_4=s^2V^{-6}, that coefficient is
-eta*(4p_4 delta(b_{-1})+b_{-1}delta(p_4)), giving (4).
The remaining task is to prove D_k!=0 for every k.

Put theta=V partial_V, c(x)=x(x-2). Then
L=s eta^2 V^{-7}c(theta). The first g-derivative of R at g=mu=0 is
s^{-1}V^3 F(theta), where the formal Laurent series F at infinity
is uniquely determined by

\[
 c(x+3)F(x)+c(x)F(x-7)=c(x),\qquad
           F(x)=\tfrac12-\tfrac3{2x}+O(x^{-2}).           \tag{5}
\]

This is the linearized equation L R_g+R_g L=A^2. Its coefficient
recursion divides by two, so it defines a unique series. Differentiating
R^{2k+1}, shifting the coefficient of x^{-1} under translations, gives

\[
 D_k=\sum_{i+j=2k}\mathcal L\bigl(b_i(x)p_j(x)\bigr),
 \quad\mathcal L(B)=[x^{-1}]F(x)B(x),                    \tag{6}
\]

where
b_i(x)=product_{ell=0}^{i-1}c(x+3-7ell) and
p_j(x)=product_{ell=1}^{j}c(x+7ell).
To verify the residue extraction, the unshifted jth summand is

\[
 \Bigl[ x^{-1}\Bigr]
 \left(\prod_{\ell=0}^{2k-j-1}c(x+3-7j-7\ell)\right)
             F(x-7j)\left(\prod_{\ell=0}^{j-1}c(x-7\ell)\right).
\]

Translation x to x+7j gives (6). Coefficients of positive A-powers
are polynomials in theta, while A^{-1}=eta^{-1}V^2(theta+2)^{-1};
therefore its first theta^{-1} coefficient is precisely D_k with
the scalar normalization in b_{-1} above.

Taking the coefficient of x^{-1} in (5) times any polynomial B gives

\[
 \mathcal L\{c(x+3)B(x)+c(x+7)B(x+7)\}=0.               \tag{7}
\]

The polynomial right side has zero residue, and formal residue is
invariant under translation. Taking B=b_{i-1}(x-7)p_j(x) in (7)
shows
mathcal L(b_i p_j)=-mathcal L(b_{i-1}p_{j+1}).
The alternating 2k+1 terms of (6) consequently sum to

\[
                              D_k=\mathcal L(p_{2k}).    \tag{8}
\]

## A polynomial certificate proves nonvanishing in every order

For an indeterminate n define

\[
 \begin{split}
 Q_n(x)&=\frac{x^2}{2}+\frac{21n+23}{2}x
                         +\frac{245n^2+511n+246}{4},\\
 r_n&=-\frac{(7n+2)(7n+4)(7n+6)(7n+11)}4.
 \end{split}
\]

The following identity is polynomial in both x and n:

\[
\begin{split}
 c(x+7(n+1))c(x+7(n+2))-r_n
   =c(x+3)Q_n(x)+c(x+7(n+1))Q_n(x+7).                   \tag{9}
\end{split}
\]

Multiplying by p_n makes the right side exactly the expression
in (7) with B=p_n Q_n, because
c(x+7)p_n(x+7)=c(x+7(n+1))p_n(x).
Thus mathcal L(p_{n+2})=r_n mathcal L(p_n) for every integer n>=0.
Using (8) and D_0=-3/2 yields

\[
 D_k=-2(7k-6)(7k-5)(7k-4)(14k-3)D_{k-1},\quad k\ge1.    \tag{10}
\]

Every factor is nonzero in characteristic zero. Hence D_k!=0
for every k, by induction. This is a universal symbolic certificate,
not an inference from finitely many values.

In any nonempty finite sum in (3), take the largest index with c_k!=0.
Its nonzero term (4) is strictly more polar than the A^2 coefficients
of every lower-index flow. It cannot be canceled. A scalar commutator,
including zero, therefore forces all c_k=0. The remaining q(P)
commutes with P. This proves Cent(P)=k[P] and no nonzero scalar mate
over k(V). Normalized coefficient trace extends mate nonexistence
to every finite algebraic coefficient cover; the centralizer equality
itself is asserted over the rational coefficient field.

## Independent exact low-order controls

The order-two operator is

\[
                    B_2=sV^{-3}A^2+\frac{gV^3}{2s}.      \tag{11}
\]

Writing B_6=sum_{j=0}^6 b_j A^j, its exact coefficients are

\[
\begin{split}
 b_6={}&s^3V^{-9},\\
 b_5={}&-18\eta s^3V^{-11},\\
 b_4={}&159\eta^2s^3V^{-13}+\tfrac32gsV^{-3},\\
 b_3={}&-780\eta^3s^3V^{-15}-\tfrac92\eta gsV^{-5},\\
 b_2={}&1800\eta^4s^3V^{-17}+\tfrac{21}2\eta^2gsV^{-7}
                    +\frac{3g^2V^3}{8s}+\tfrac32\mu sV^{-2},\\
 b_1={}&\tfrac32\eta\mu sV^{-4},\\
 b_0={}&18\eta^4gsV^{-11}-\frac{3\eta^2g^2V^{-1}}{4s}
            -2\eta^2\mu sV^{-6}-\frac{g^3V^9}{16s^3}
            +\frac{3g\mu V^4}{4s}.
\end{split}                                               \tag{12}
\]

Only the square root through A^{-4} is needed to derive (12), because
P has order four. The retained checker rederives those root coefficients,
checks their square to the required orders, multiplies PR, and compares
with (12), before forming either commutator.

The order-two commutator is

\[
\begin{split}
 [P,B_2]={}&-3\eta^2gsV^{-7}A^2\\
 &+\left(12\eta^3gsV^{-9}+\frac{3\eta g^2V}{s}
                                      -2\eta\mu sV^{-4}\right)A\\
 &+36\eta^4gsV^{-11}+\frac{3\eta^2g^2V^{-1}}{2s}
                                      +\eta^2\mu sV^{-6}.
\end{split}                                               \tag{13}
\]

For the sixth-order positive part, the needed commutator coefficient is

\[
 [P,B_6]_{A^2}=
 \frac{9\eta^2g}{8sV^{21}}
 \left(10208\eta^4s^4-12\eta^2gs^2V^{10}
                    -11g^2V^{20}+4\mu s^2V^{15}\right).  \tag{14}
\]

The full coefficients at A^2,A,1 are retained in the checker's output;
all higher coefficients vanish exactly. Formula (14) has the unique
most polar term 11484 eta^6*g*s^3*V^{-21}. It is nonzero for every
allowed eta,g,s and is independent of mu. No scalar multiple of (13)
can cancel it, because (13)'s A^2 coefficient is supported only at
V^{-7}. Thus

\[
                    c_6[P,B_6]+c_2[P,B_2]
                       \text{ scalar}\quad\Longrightarrow\quad
                              c_6=c_2=0.                 \tag{15}
\]

These finite identities independently check the first two nonzero entries
of the all-order residue sequence above: D_0=-3/2 and D_1=198.
They also give the truncated centralizer k+ kP+ kP^2 through order
eight. The unrestricted theorem uses the symbolic recurrence, not
a finite prefix of these computations.

## The globally admitted square ordering is closed

Use the exact global operators from
[the global admission note](GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md):

\[
 \begin{split}
 A_0&=Q+3XQ^2-XB/3=a-v,\qquad C=R A_0-2/3,\\
 S_b&=-3R A_0^2+(18XQ-8)A_0+9Q-3hRX,\\
 S_{\rm sq}&=S_b+\epsilon h^2(X^3C^2)^2,
                                      \qquad\epsilon\ne0.
 \end{split}                                             \tag{16}
\]

These are actual finite members of the global centralizer of exact R,
and S_sq has the prescribed Long reduction. After centering
z=a-2/(3r), C=r(z-v). The finite rescaling h=u^5, v=u^2V, z=uA,
subtracting 4/(3r) and dividing by u^2, has integral critical reduction
(1), with g=-3r, mu=2, eta=1 and s^2=epsilon*r^4. Adjoining a
constant square root of epsilon is harmless for nonexistence.

A hypothetical finite mate can first be traced to rational position
coefficients and scaled to commutator h. Its rescaling has finitely
many poles in the u-Gauss valuation. At a highest pole, the residue
commutes with P and is q(P) for a finite polynomial q. Subtract
the corresponding polar multiple q of the integral rescaled S_sq.
This preserves the commutator and strictly lowers the parameter pole.
After finitely many steps an integral mate would reduce to a nonzero
constant mate of P, already excluded. Each subtraction is finite;
no uniform bound on an unknown partner's order is required.
Thus S_sq has no finite rational-position or finite algebraic-cover
mate of any order, and in particular no globally admitted finite mate.
The same argument excludes finite deformations with this exact
critical face, even if their generic order increases.

## Evidence and scope

[The exact symbolic controls](check_square_ordered_quartic_face.py)
verify (9) with both indices symbolic, rederive the formal F coefficients,
check (6)--(10) at low orders as independent controls, and retain the
full finite B_6 calculation, both commutators, corrupted-coefficient
rejection and commuting control g=mu=0. The all-order conclusion
uses the written residue reduction and induction, not the finite
controls alone. No actual Long mate, bounded lift, DC1/DC2 resolution,
formal verification, independent whole-proof replay or external
review is claimed. Mixed critical faces with a leading cubic term
remain outside this square-only theorem.
