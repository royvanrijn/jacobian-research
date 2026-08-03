# Quintic graph and relative-linear Schur frontiers

## Status

No `HC(4)` counterexample is constructed here. Starting with the explicit
five-variable Meng--Yang potential, this note proves two reductions.

> **Theorem `HC4MYG5N` -- rational quintic-graph normal slice.** Let
> \(r=R(x,y,p,q)\), with \(\deg R\leq5\), and pull back the scaled
> Meng--Yang family \(\Psi_{L,M,N}=LA^2+MA+NB\), where \(LN\ne0\).
> Constant Hessian determinant forces the homogeneous quintic graph part
> \(R_5\) to have a constant direction in the \((y,p,q)\)-space and
>
> \[
> R_5(0,y,p,q)=\zeta y^5.
> \]
>
> Over \(\mathbb Q\), necessarily \(\zeta\ne0\). The kernel direction then
> lies in the \((p,q)\)-plane, and \(R_5\) is on one of the charts
>
> \[
> R_5=\zeta y^5+xT_4(x,y,p-\kappa q),\qquad
> R_5=\zeta y^5+xU_4(x,y,q).
> \]
>
> This is a necessary normal form, not a completed determinant identity.

> **Theorem `HC4MYG5K` -- complete degree-two lower-trace obstruction on
> the first constant-kernel chart.** For the v2 potential, no
> degree-at-most-five graph on the
> \(\partial_q\)-kernel chart whose plane trace has the form
>
> \[
> \kappa y^5+d y^3p+\rho y^2q+h y^3+e yp^2
> +T_{\le2}(y,p,q)
> \]
>
> has constant Hessian determinant. All ten coefficients of
> \(T_{\le2}\), the target constant, and every allowed \(x^2U_3\) repair are
> included. For \(a=[pq]T\), the generic branch gives a three-equation unit
> ideal in \(\mathbb Q[a,\rho]\), with the kernel-denominator chart excluded
> separately. A quartic-normal square forces \([q^2]T=0\), and the remaining
> nonzero \([yq]T\) branch ends in two coprime polynomials in \(\rho\).

> **Theorem `HC4MYR89` -- relative-linear obstruction through degree
> 89.** Complete the polynomial unit pivot \(t=A\), eliminate it, and allow
> an arbitrary base-dependent linear change of the remaining dual variables
> with \(C(x,y)\in\operatorname{SL}_2(\mathbb Q[x,y])\). If every entry of
> \(C\) has total degree at most \(89\), no constant-Hessian member of this
> class has a gradient collision. Its gradient is a polynomial automorphism.

The first theorem uses the projective-gradient cone classification only
after the leading graph square has been isolated. The second closes the
first cone-compatible transverse slice. The third is a genuine polynomial
Schur elimination. The open search is now in broader cubic and quartic traces
on the two quintic cone charts, residual-linear corrections of degree at
least 90, or transformations nonlinear in the residual dual variables.

The complementary exact calculation in
[the graph-obstruction note](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md)
solves the determinant equation on \(x=0\) as a unit-affine equation for the
first normal jet, proves the same vertical-quintic trace rigidity without an
external cone theorem, and excludes one complete sparse quintic trace family
by immutable first-transverse coefficients. Thus the normal-slice theorem
and the transverse exclusion constrain different layers of the same
degree-five problem.

## 1. Five-variable input

Put \(u=1+xy\) and

\[
\begin{aligned}
A&=u^3p+3xu^2q-x^3r,\\
B&=y^2u(4+3xy)p+
   \bigl(y+3xy^2(4+3xy)\bigr)q+(2x-3x^2y)r.
\end{aligned}
\]

For \(\Psi_{L,M,N}=LA^2+MA+NB\), the five-variable Hessian determinant is
\(8LN^4\). The v2 member is \((L,M,N)=(1,13,2)\), and

\[
(1,-3/2,0,0,0),\qquad(-1,3/2,0,0,0)
\]

have common gradient \((0,0,-1/2,0,0)\). Any graph through both points
transfers this collision because the omitted \(r\)-gradient component is
zero there. The affine and degree-at-most-four graph cases are already
excluded in
[`HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md`](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md).

That note also excludes every arbitrary-degree graph whose complete 1-jet
along \(x=0\) has degree at most four, because an \(x^2\)-divisible tail is
invisible to the contradictory quartic plane calculation.

## 2. Complete two-slope graph jet

Let \(D_R\) be the Hessian determinant after substituting
\(r=R(x,y,p,q)\), and restrict to

\[
(x,y,p,q)=(0,t,ct,dt).
\]

Only the \(x^0\)- and \(x^1\)-coefficients of \(R\) enter the Hessian on
this line. The checker retains every such coefficient through degree five,
91 parameters in total, and expands the determinant only through \(t^{10}\).
The leading square is

\[
[t^{10}]D_R
=256N^4\bigl(\partial_qR_5(0,1,c,d)\bigr)^2.
\]

It kills all 15 \(x\)-free quintic coefficients containing \(q\).
The \(t^8\) coefficient is a second square. Solving it and substituting into
the next row gives

\[
\begin{aligned}
[t^7]D_R=-LN^3&
(4s_{0140}c^3+3s_{0230}c^2+2s_{0320}c+s_{0410})\\
&\cdot(2s_{0140}c^4-s_{0230}c^3-4s_{0320}c^2
       -7s_{0410}c-10\zeta).
\end{aligned}
\]

The coefficient ring is a domain. If the second factor vanishes, the first
factor vanishes as well. Hence every branch satisfies

\[
R_5(0,y,p,q)=\zeta y^5.
\]

For the next rows set

\[
\begin{gathered}
\eta=[q^3]R_3,\quad\theta=[pq^2]R_3,\quad
\tau=[p^2q]R_3,\quad\iota=[yq^2]R_3,\\
\sigma=[ypq]R_3,\quad\rho=[y^2q]R_3,\\
a=[p^4]R_4,\quad b=[yp^3]R_4,\quad
e=[y^2p^2]R_4,\quad f=[y^3p]R_4.
\end{gathered}
\]

The exact \(t^6\) coefficient is

\[
N^3\bigl(NS(c,d)^2+L\zeta T(c,d)\bigr),
\]

where

\[
\begin{aligned}
S={}&-8ac^3+(16\tau-6b)c^2+32\theta cd
 +(16\sigma-4e)c+48\eta d^2+32\iota d+16\rho-2f+89,\\
T={}&240\eta d^2+160\theta cd+80\tau c^2
 +160\iota d+80\sigma c+80\rho+492.
\end{aligned}
\]

Extreme square coefficients followed by the \(c^5t^5\) and \(c^3t^5\)
rows force

\[
a=\eta=\theta=\iota=\tau=\sigma=e=b=0.
\]

With \(\beta=16\rho-2f+89\), the remaining \(t^6\) equation is

\[
N\beta^2+L\zeta(80\rho+492)=0.
\]

If \(\zeta=0\), then \(\beta=0\), and the \(ct^5\) coefficient becomes

\[
2LN^3Q(\rho),\qquad Q(\rho)=160\rho^2+1968\rho+6021.
\]

Its discriminant is \(19584=576\cdot34\), so \(Q\) has no rational root.
Thus every rational graph candidate has \(\zeta\ne0\). Over an algebraic
closure the zero-\(\zeta\) axis branch remains over
\(\mathbb Q(\sqrt{34})\); this is not a complex exclusion.

## 3. Projective-gradient cone classification

The pullback potential has leading part

\[
h_{16}=LF^2,\qquad F=x^3R_5,\qquad\deg F=8.
\]

For a homogeneous degree-\(d\) form \(F\), with \(g=\nabla F\) and
\(H=\operatorname{Hess}F\), the rank-one determinant identity and Euler
relations give

\[
\det\operatorname{Hess}(F^2)
=2^4\frac{2d-1}{d-1}F^4\det H.
\]

At \(d=8\), the scalar is \(240/7\). Constancy of the full determinant
therefore forces \(\det\operatorname{Hess}(x^3R_5)=0\).

The low-dimensional Gordan--Noether theorem supplies a constant direction
\(v\) with \(D_v(x^3R_5)=0\). Its \(x\)-component is zero: otherwise every
translation orbit in direction \(v\) meets \(x=0\), where \(x^3R_5\)
vanishes, forcing the whole invariant polynomial to vanish. Thus
\(D_vR_5=0\) for a direction in \((y,p,q)\).

For a rational graph this direction can be chosen over \(\mathbb Q\): the
four first partial derivatives of \(x^3R_5\) are coefficient vectors over
\(\mathbb Q\), so a linear dependence after scalar extension is already a
dependence over \(\mathbb Q\).

On the rational branch \(\zeta\ne0\), restriction to
\(R_5(0,y,p,q)=\zeta y^5\) forces the \(y\)-component of \(v\) to vanish.
Normalizing the two charts of the remaining projective line gives exactly
the two forms in `HC4MYG5N`.

This determines the reduced top support, not the later normal-cone
multiplicities or lower \(X_0\)-adic module.

## 4. Complete degree-two lower trace on the \(\partial_q\) chart

The plane normal-jet identity from
[`HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md`](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md)
becomes especially effective after imposing a constant top-kernel direction.
For the v2 member, consider the trace

\[
 T=\kappa y^5+d y^3p+\rho y^2q+h y^3+e yp^2
   +a pq+\ell q+T_{\le2}(y,p),                       \tag{4.1}
\]

where all six coefficients of the q-free polynomial \(T_{\le2}\), together
with \(\ell\), are arbitrary. Put

\[
 V=16\rho-2d+89,\quad
 P=8\rho^2+99\rho+279,\quad
 Q=160\rho^2+1968\rho+6021,\quad W=2a+3.           \tag{4.2}
\]

The degree-six and degree-five plane faces are

\[
\begin{aligned}
V^2+2\kappa(20\rho+123)&=0,\\
32Va+40a\kappa-8Ve-12V\rho-39V+60\kappa+Q&=0.
\end{aligned}                                      \tag{4.3}
\]

The coefficient of \(y^3q\) in the forced quartic normal jet is
\(-(WV-P)/2\). On the generic \(W\ne0\) chart it vanishes exactly when

\[
 V=P/W,\qquad
 \kappa=-\frac{V^2}{2(20\rho+123)},                 \tag{4.4}
\]

and (4.3) uniquely determines \(e\). The plane equation then uniquely fixes
the complete quartic normal jet \(S=(\mathcal F(T)-C)/64\), for an arbitrary
target constant \(C\).

The divisions in (4.4) lose no generic branch. If \(P=V=0\), the plane
faces force \(\kappa=Q=0\), but
\(\operatorname{Res}(P,Q)=16959456\). Also
\(P(-123/20)=-2727/100\), so \(20\rho+123\) cannot vanish. The separate
boundary \(W=0\) is treated below.

Write \(R=T+xS+x^2U\). The first transverse variation of the determinant
under \(U\) is \(-192U(0,y,p,q)\). Since \(\deg U\le3\), it cannot alter
the coefficients of total degree at least four in \([x]D_R\).

First set \(a=0\). Exact extraction gives

\[
\begin{aligned}
[y^7x]D_R&=\frac{8P}{20\rho+123}
 (656\rho^3+11166\rho^2+70911\rho+168912),\\
[y^4qx]D_R&=48(16\rho^3+192\rho^2+801\rho+54).
\end{aligned}                                      \tag{4.5}
\]

Every coefficient in \(h,T_{\le2}\), \(\ell\), and \(C\) cancels from
(4.5). The two terminal cubics have resultant

\[
-108117004020524928
=-2^7 3^{17}\cdot11\cdot13\cdot53\cdot863\ne0.    \tag{4.6}
\]

Thus the \(a=0\) diagnostic slice has no point over any characteristic-zero
field.

Now retain the new coefficient \(a=[pq]T\). Normalize the three immutable
coefficients \([xy^7]D_R,[xy^5p]D_R,[xy^4q]D_R\) by the nonzero factors in
(4.4), calling the resulting polynomials \(E_7,E_6,E_5\). Their bidegrees in
\((a,\rho)\) are \((1,3),(3,7),(1,3)\). The two smaller equations are

\[
\begin{aligned}
E_7={}&1280a\rho^3+23232a\rho^2+156744a\rho+383022a\\
     &+1968\rho^3+33498\rho^2+212733\rho+506736,\\
E_5={}&64a\rho^2+768a\rho+3114a
       -16\rho^3-192\rho^2-801\rho-54.             \tag{4.7}
\end{aligned}
\]

The exact lexicographic Gröbner basis is

\[
 (E_7,E_6,E_5)=(1)\subset\mathbb Q[a,\rho].        \tag{4.8}
\]

As an elimination cross-check, the two resultants against \(E_7\) have
monic gcd \(32\rho^2+384\rho+1557\); the third equation removes that
apparent projection. On the missing \(W=0\) chart, the kernel equation forces
\(P=0\). The surviving nonzero-\(V\) branch has

\[
 [xy^4q]D_R\bmod P=540(25\rho+141),
 \qquad \operatorname{Res}(P,25\rho+141)=-15552,   \tag{4.9}
\]

while \(V=0\) is already excluded by \(\operatorname{Res}(P,Q)\).
Consequently (4.1) is empty over every characteristic-zero field.

To complete the degree-two lower trace, now adjoin
\(b yq+cq^2\). Before solving for \(V,\kappa,e\), the new degree-five plane
faces and one coefficient of the quartic normal jet are

\[
\begin{aligned}
[y^5]\mathcal F(T)&=128b(4V+5\kappa),\\
[y^4q]\mathcal F(T)&=256c(4V+5\kappa),\\
[y^2q^2]S_4&=256c^2.                               \tag{4.10}
\end{aligned}
\]

The \(\partial_q\)-kernel chart requires \(S_4\) to be q-independent, so
the square in (4.10) forces \(c=0\). If \(b=0\), this is (4.1). If
\(b\ne0\), then \(4V+5\kappa=0\). Together with the first equation of
(4.3), this gives \(V=0\) or

\[
 V=\frac{8(20\rho+123)}5,qquad
 \kappa=-\frac{32(20\rho+123)}{25}.                \tag{4.11}
\]

The \(V=0\) branch would force both \(P=0\) and \(Q=0\), which is impossible
by \(\operatorname{Res}(P,Q)=16959456\). On (4.11), the kernel equation
determines

\[
 a=\frac{40\rho^2+15\rho-1557}{16(20\rho+123)},
\]

and the second equation of (4.3) uniquely determines \(e\). Exact extraction
then gives

\[
\begin{aligned}
[xy^7]D_R&=\frac{192(20\rho+123)}{25}
 (160\rho^2+1968\rho+5841),\\
[xy^5p]D_R&=\frac3{25(20\rho+123)}
 \bigl(3968000\rho^4+98457600\rho^3+896972480\rho^2\\
 &\hspace{38mm}+3553345920\rho+5132433339\bigr).   \tag{4.12}
\end{aligned}
\]

Both coefficients are independent of \(b,\ell,h,T_{\le2}\), and \(C\).
The resultant of the quadratic and quartic factors in (4.12) is

\[
986335129354383654912000
=2^{27}3^8 5^3\cdot11\cdot24223\cdot33629\ne0.    \tag{4.13}
\]

Thus the nonzero-\(b\) branch is empty over every characteristic-zero field
as well. This proves `HC4MYG5K` for the complete
\(T_{\le2}(y,p,q)\). Broader cubic and quartic trace terms on this chart and
the other constant-kernel charts remain open.

## 5. Finite-field scout

For \(L=1,N=2\), the projected axis locus is

\[
2\beta^2+\zeta(80\rho+492)=0.
\]

The checker exhausts this quadric over \(\mathbf F_{101}\) and
\(\mathbf F_{103}\), then separately counts zero-\(\zeta\) points that also
pass \(Q(\rho)=0\). A small rational point on the nonzero-\(\zeta\)
projection is

\[
(\rho,\beta,\zeta)=\left(0,6,-\frac6{41}\right).
\]

This is deliberately recorded only as a **projected jet survivor**. It is
not a graph satisfying the full determinant identity, and no collision
claim is attached. The next search must impose the remaining line layers
and an off-axis normal slice before rational reconstruction.

There is also an exact, stronger near miss in the graph-obstruction note.
Its trace is

\[
 \frac{51}{100}y^5-\frac{47}{10}y^3p-\frac{123}{20}y^2q,
\]

and its uniquely forced normal jet makes the determinant equal to
\(17165601/25\) on all of \(x=0\) while placing both marked collision points
on the graph. It fails at the first transverse term,
\([xy^7]D=22032/125\). Hence rational collision containment and complete
plane flattening are both achievable; transverse compatibility is the
genuine remaining gate.

For the stronger transverse slice of Section 4, the two cubics in (4.5)
have no common root over \(\mathbf F_{101}\) or \(\mathbf F_{103}\).
Exactly one modular common root appears at each nontrivial bad resultant
prime:

\[
2\pmod {11},\quad10\pmod {13},\quad34\pmod {53},
\quad717\pmod {863}.
\]

These are bad-characteristic collapses, not reconstructible rational points.
After adjoining \(a pq\), the three-equation locus has no point over
\(\mathbf F_{103}\). Over \(\mathbf F_{101}\) its only raw point is
\((\rho,a)=(9,4)\), where \(P(9)=0\); saturation by the parameterization
denominators removes it. Thus both good-prime admissible fibers are empty.
Over \(\mathbb Q\), the full ideal is \((1)\), so this enlarged slice has no
candidate to reconstruct. The exceptional nonzero-\(yq\) pair in (4.12)
also has no common root over \(\mathbf F_{101}\) or \(\mathbf F_{103}\).
Its bad-prime collapses are

\[
0\pmod {11},\qquad4365\pmod {24223},\qquad
30101\pmod {33629}.
\]

The full-potential replay recomputes all three generic witnesses and all
three exceptional-branch witnesses directly from the original five-variable
formula at independent rational specializations. Thus the truncated
\(x^3\)-jet extraction is checked by a second exact route; the finite-field
counts remain scouting data rather than the characteristic-zero proof.

## 6. Relative-linear unit-pivot descent

The explicit polynomial coordinate \(t=A\) from
[`HC5_NONLINEAR_TORIC_DESCENT.md`](HC5_NONLINEAR_TORIC_DESCENT.md)
has a unit quadratic pivot. After eliminating it and applying any
\(C(x,y)\in\operatorname{SL}_2(\mathbb Q[x,y])\), the reduced potential is

\[
\psi=f(x,y)+2G_1(x,y)r+2G_2(x,y)s,\qquad G=\gamma C,
\]

with

\[
\gamma=\bigl(-yP(xy),xQ(xy)\bigr),
\]

\[
P(v)=18v^5+81v^4+120v^3+60v^2-1,\qquad Q(v)=(v+1)(v+2).
\]

For every upper-left Hessian block,

\[
\det\operatorname{Hess}\psi=16(\det DG)^2.
\]

Moreover \(\gamma=(y,2x)+O(2)\) and \(\gamma(0)=0\), so derivatives of
\(C\) do not enter \(DG(0)\), and \(\det DG(0)=-2\det C(0)=-2\).

The entries of \(\gamma\) have degrees 11 and 5. If the entries of \(C\)
have degree at most \(d\), then \(\deg G\le d+11\). For \(d\le89\),
Moh's plane theorem through degree 100 makes every constant-Jacobian
\(G\) a polynomial automorphism.

The last two coordinates of \(\nabla\psi\) are \(2G(x,y)\), so they recover
\((x,y)\). With the base fixed, the first two recover \((r,s)\) through
\(2DG(x,y)^{\mathsf T}\). Hence \(\nabla\psi\) is injective and, over
characteristic zero, a polynomial automorphism. This proves `HC4MYR89`.

The result permits constant determinants but excludes collisions. Degree
at least 90 and transformations nonlinear in \(r,s\) remain open.

## 7. Reproduction

Run

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_graph_normal_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_graph_normal_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_q_kernel_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_q_kernel_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_relative_linear_obstruction.py
```

The first checker retains the complete 91-parameter line jet, verifies the
displayed identities, checks both kernel charts, enumerates the projected
finite-field loci, and replays the rational projected point. The second
derives the two-parameter plane parametrization, proves the three-witness
unit ideal and its exceptional-chart resultant, classifies the remaining
degree-two lower trace by a quartic-normal square and a second transverse
resultant, exhausts the finite-field fibers, and independently replays both
branches from the original five-variable potential. The third verifies the block
determinant, degree ledger, and local Jacobian. Gordan--Noether and Moh are
external theorem inputs.

The remaining exact targets are:

1. add the broader cubic and quartic trace terms on the \(\partial_q\) chart
   and run the same first-transverse ideal on the other constant-kernel
   charts;
2. residual-linear corrections of degree at least 90; and
3. unit-pivot changes nonlinear in the two residual dual variables.

Repeating an arbitrary quartic or quintic potential search addresses none
of these loci.

## References

- [Gordan--Noether/singular-Hessian classification](https://arxiv.org/abs/1501.05168).
- T. T. Moh,
  [*On the Jacobian conjecture and the configurations of roots*](https://doi.org/10.1515/crll.1983.340.140).
