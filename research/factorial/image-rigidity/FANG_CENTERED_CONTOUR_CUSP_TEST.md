# A centered same-contour test at the cubic Gamma cusp

[Fang, Hexagon 2611.01129v1, Section 5, equations (5) and (8)](sources/fang-2611.01129v1.pdf)
uses neighboring contour integrals on one contour, centers an insertion,
and in one regime applies exact integration by parts before dividing by
a small saddle-separation factor. The paper proves eventual interlacing
for its own set-partition polynomials. Its threshold and root geometry
do not supply a factorial cutoff.

Here is a test on one concrete obligation in
[the active-cycle and global-error gap](CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md):
can a common-contour centered signal at the first ordinary cusp be
separated from its complete-contour error? Set `d=3`, `a=1` in the
physical positive-axis moment integral. Its cusp is `S_c=2`, `w_c=16`.
For integer `m≥1` put

`I_m(w)=∫₀^∞ S^m(w+S³)^m e^(−mS) dS`.

At `w=16`, write `Φ(S)=log S+log(16+S³)−S` and
`f(S)=1/(16+S³)`. Differentiation under the integral gives the exact
same-contour identity

`J_m := I'_m(16)/m − I_m(16)/24
     = ∫₀^∞ exp(mΦ(S)) (f(S)−f(2)) dS`.              (1)

No separate estimate of the two terms on the left is needed.
The exact derivative is

`Φ'(S)=−(S−2)³(S+2)/(S(S³+16))`.

Integration by parts on this same real contour also gives
`∫₀^∞e^{mΦ(S)}(SΦ'(S)+1/m)dS=0`, since the endpoint terms
vanish. Thus the insertion in (1) may be changed by any constant
multiple of `SΦ'+1/m`. Here `SΦ'` has a cubic zero at `S=2`, so
this exact cancellation does not remove the linear centered
insertion. It tests the integration-by-parts step without dividing
by a saddle-separation parameter.

Thus `S=2` is the unique global maximum on the positive axis; the
endpoints have exponentially smaller phase. With `t=S−2`, exact Taylor
coefficients are

`Φ(2+t)=log 48−2−t⁴/48+t⁵/80−t⁶/288+O(t⁷)`,

`f(2+t)=1/24−t/48+0·t²+t³/288+O(t⁴)`.

The leading insertion term is odd against the symmetric quartic
model, so its leading integral cancels. Rescale `t=m^(−1/4)u`.
The first even surviving term is
`−u⁶/(3840)m^(−1/2)` inside the rescaled integrand. Standard
dominated Taylor expansion on a shrinking symmetric neighborhood,
with the strict global phase gap outside, yields the scoped asymptotic

`J_m = −(3/320) M₂ e^{m(log 48−2)} m^(−3/4)
       + o(e^{m(log 48−2)}m^(−3/4))`,                 (2)

where `M₂=∫_R u²e^(−u⁴/48)du=(1/2)48^(3/4)Γ(3/4)>0`.
The moment identity `∫u⁶e^(−u⁴/48)du=36M₂` gives the coefficient.
The [exact coefficient checker](verify_fang_cusp_centered_signal.py),
run with Python 3.12 and SymPy 1.14, confirms the factorization and
Taylor data. It is an internal algebra check, not a proof of the
whole-contour estimate.

This test produces a nonzero centered signal on the *positive real*
cycle, but the signal is only order `m^(−1/2)` relative to `I_m(16)`.
The unresolved factorial-zero problem needs a bound on the complete
original cycle after continuation into the relevant complex cut cells,
including its saddle coefficients, other active saddles and connectors.
The existing stationary-cover result does not identify those
coefficients. A qualitative `o(1)` saddle error is insufficient for
the signal in (2). The next test must construct a common contour in
one specified complex small-angle cell and prove a uniform error
`o(m^(−1/2))` relative to its dominant packet, with a nonvanishing
centered coefficient there. Equation (2) does not address the
multivariate finite residual `15≤N<N_ang`.
