# Linear escape of a one-sided isolated Gaussian interface end

Let \(0<r<1\). In orthonormal coordinates about the Gaussian origin,
suppose a stationary pair interface is a finite continuous graph
\(y=h(x)\) for \(x>S_0\). In this half-plane the lower label fills
the region below the graph and the upper label fills the region
above it, up to null sets. The remaining half-plane may contain
other labels. Let \(\beta\in(-1,1)\) be the adjusted pressure
difference of the lower and upper labels. Define the flat offset
\(b\) by
\[
  \psi((1-r)b)=\beta,\qquad
  \psi(u)=2\Phi(u/\sigma)-1,\qquad \sigma=\sqrt{1-r^2}.
\]
If \(h(x)\ge b\) for all sufficiently large \(x\), then either
\(h=b\) on a tail or there is a constant \(c>0\) such that
\(h(x)-b\ge c x\) for all sufficiently large \(x\).
If \(h(x)\le b\) eventually, the alternatives are exact
straightness and \(h(x)-b\le-c x\) on a tail.
Every isolated finite graph end, with or without a sign condition,
satisfies \(|h(x)|\le Cx\) for a constant \(C\) and all large \(x\).

The argument first truncates the height to obtain a bounded positive
quantity with a uniform lower bound. The lower bound then gives
linear growth by a comparison of tail infima. No monotonicity or
upper bound on the height is needed.
[Conditional wall equation and bounded one-sided ends](THREE_CELL_ISOLATED_END_RIGIDITY.md)

## An unconditional linear upper bound

Choose \(C_1>0\) so that \(\psi(-C_1)<\beta<\psi(C_1)\).
For a sufficiently large fixed \(R\), the probability
of \(S_0<S_x\le rx+R\) is as close to one as needed for all
large \(x\), where \(S_x\sim N(rx,\sigma^2)\).
If \(r h(x)\) lay above
\(\sup_{S_0\le s\le rx+R}h(s)+C_1\), the conditional contrast
on this event would be at most \(\psi(-C_1)\).
The contrast is bounded by one on the complement, so its total
expectation would be strictly below \(\beta\).
In the same way, \(r h(x)\) cannot lie below
\(\inf_{S_0\le s\le rx+R}h(s)-C_1\). Here \(S_0\) can be
increased to lie in the continuous graph domain.

Let \(B(x)=\sup_{S_0\le s\le x}|h(s)|\). The two comparisons
imply \(r|h(x)|\le B(rx+R)+C_1\) for all large \(x\).
Because \(B\) is nondecreasing, taking a supremum over
the prefix up to \(x\) gives
\[
                  r B(x)\le B(rx+R)+C_2                  \tag{OE0}
\]
for all large \(x\), increasing \(C_2\) to include the initial
compact interval. With \(a=R/(1-r)\) and
\(H(t)=B(t+a)+C_2/(1-r)\), (OE0) gives \(rH(t)\le H(rt)\).
Iterate until \(r^nt\) first enters a fixed compact interval
in the tail where the comparison holds. The resulting bound
on \(H(t)/t\) proves the upper bound.

## A positive lower bound

Increase \(S_0\) so that \(q=h-b\ge0\) there. With \(d=(1-r)b\),
let \(S_x\) have law \(N(rx,\sigma^2)\), and let \(p_\sigma\)
be the centered normal density of standard deviation \(\sigma\).
The conditional wall equation, with the flat pressure subtracted,
is
\[
  0=\int_{s>S_0}[\psi(d+q(s)-r q(x))-\psi(d)]
             p_\sigma(s-rx)\,ds+E(x),\qquad
  |E(x)|\le 2\Pr(S_x\le S_0).                             \tag{OE1}
\]
Indeed, the conditional vertical probability gives the displayed
integrand on \(s>S_0\), and the other part is a difference of
contrasts bounded in magnitude by two.

Choose \(M>0\) sufficiently small that
\[
  L_-=\min_{[d-rM,d+M]}\psi'>0,\qquad
  L_+=\max_{[d-rM,d+M]}\psi',\qquad
  \frac{L_-}{rL_+}>1 .
\]
Such a choice is possible since \(\psi'(d)>0\) and \(r<1\).
Put \(g=\min(q,M)\) on \(s>S_0\).
When \(q(x)<M\), replacing \(q(s)\) by \(g(s)\) can only
decrease the integral in (OE1). Apply the mean-value formula
to the new integrand. Its slope \(K(x,s)\) lies in
\([L_-,L_+]\), so
\[
 \begin{split}
  rL_+ g(x)&\ge \int_{s>S_0} K(x,s)g(s)p_\sigma(s-rx)\,ds+E(x)\\
            &\ge L_-\int_{s>S_0}g(s)p_\sigma(s-rx)\,ds
                           -2\Pr(S_x\le S_0).
 \end{split}                                             \tag{OE2}
\]
If \(q\) is identically zero on a tail there is nothing to prove.
Otherwise choose a compact interval \([a,a']\subset(S_0,\infty)\)
of positive length on which \(g\) has a positive lower bound.
For \(rx>a'\), its contribution to the integral in (OE2) is
bounded below by a positive constant times \(p_\sigma(a-rx)\).
This term strictly dominates \(2\Pr(S_x\le S_0)\) after
multiplication by \(L_-\) and for all sufficiently large \(x\).
The ratio grows exponentially in \(x\) by the Gaussian tail bound
and the strict inequality \(a>S_0\). In particular \(g(x)>0\)
for all sufficiently large \(x\).

Fix \(0<\ell<r<u<1\). For large \(x\), the interval
\([\ell x,ux]\) is disjoint from \([a,a']\).
When \(q(x)<M\), the contribution from \([a,a']\) absorbs
the error and (OE2) gives
\[
           g(x)\ge\int_{\ell x}^{ux}g(s)p_\sigma(s-rx)\,ds .
                                                               \tag{OE3}
\]
When \(q(x)\ge M\), (OE3) also holds since \(g(x)=M\),
\(g(s)\le M\), and the Gaussian mass is at most one.
The retained Gaussian mass is at least \(1-C e^{-c_1x^2}\)
for constants \(C,c_1>0\).

Choose \(R_0\) large enough that (OE3) holds for \(x>R_0\), the
mass lower bounds there are positive, and \(g\) has a positive
minimum \(m_1\) on \([\ell R_0,R_0]\). Iterate (OE3) with its
Gaussian transition normalized on \([\ell x,ux]\), stopping at
the first point \(x_N\le R_0\). Each transition reduces the
parameter by a factor at most \(u\), so the stop is finite and
\(x_N\in[\ell R_0,R_0]\). The preceding points satisfy
\(x_{N-1-j}>R_0u^{-j}\). Along every path the product of the
retained masses is thus at least
\(\prod_{j\ge0}(1-Ce^{-c_1R_0^2u^{-2j}})>0\).
The resulting bound \(g(x)\ge m_1\) times this product supplies a
constant \(m_*>0\) with \(q(x)\ge m_*\) for all sufficiently
large \(x\).

## From the lower bound to linear growth

Put \(m(x)=\inf_{s\ge x} q(s)\), and fix
\(0<\delta<(1-r)m_*\). Choose \(R>0\) so that
\(p=\Pr(N(0,\sigma^2)\ge-R)\) satisfies
\[
              p[\psi(d+\delta)-\psi(d)]>2(1-p).
\]
For large \(x\), if \(r q(x)<m(rx-R)-\delta\), then on
\(S_x\ge rx-R>S_0\) the conditional pair contrast exceeds
the pressure by at least \(\psi(d+\delta)-\psi(d)>0\).
On the complementary event the difference is at least \(-2\).
This contradicts the wall equation. Hence
\(r q(x)\ge m(rx-R)-\delta\) for all large \(x\).
Taking the infimum over larger \(x\) gives
\[
                     r m(x)\ge m(rx-R)-\delta .
                                                               \tag{OE4}
\]
Set \(a=R/(1-r)\) and \(H(t)=m(t-a)-\delta/(1-r)\). Then
\(H(rt)\le rH(t)\) for all sufficiently large \(t\). Also
\(H(t)\ge m_*-\delta/(1-r)>0\) there. Iterating until the first
entry of \(r^nt\) in a fixed compact interval of this positive
tail bounds \(H(t)/t\) below by a positive constant. Thus
\(q(x)\ge m(x)\ge c x\) for all large \(x\). Reversing the
vertical coordinate interchanges the lower and upper labels,
negates \(b\), and proves the negative version.

## Properness of a centered affine-separated graph

Consider a balanced three-cell maximizer of the centered form with
finite heights in
[the affine-separation reduction](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md#finite-heights-and-centered-separation).
For \(s>0\), let \(h_+(s)=h(s)\) and \(h_-(s)=h(-s)\).
The two outer labels have the same pressure, so the two outward
profiles have the same flat offset \(b\). The axis identity (ASR3)
also gives the same pushforward of \(p_\sigma(s)\,ds\) under
\(h_+\) and \(h_-\) on \(s>0\).

Suppose each profile is eventually on a fixed side of \(b\).
If either is exactly flat on a tail, analyticity on that open
half-line makes it constant throughout. Equality of the height
distributions and continuity then make the other profile constant
as well. Such a partition is invariant under the central reflection
with transposition, and is excluded for a maximizer by the strict
outward slope theorem.
[Strict slope for reflected maximizers](THREE_CELL_REFLECTION_REDUCTION.md#the-translation-equation-fixes-the-sign-of-the-slope)

Otherwise each profile escapes at least linearly in absolute value
by the preceding result. Their signs of escape agree. If one
diverged to \(+\infty\) and the other to \(-\infty\), their ranges
would have different bounded endpoints because their heights are
finite and continuous up to the separating line. Their pushforward
distributions could then not agree. Thus both profiles tend to the
same signed infinity. The preimage of every compact interval of
finite heights is compact in either closed half-line.
