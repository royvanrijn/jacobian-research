# Universal Keller-fiber multiplicity

This note answers the universal nonuniqueness question for finite etale
Keller fibers.  The answer is stronger than the proposed three-class lower
bound: over every characteristic-zero field, every algebra of rank at least
three occurs in infinitely many stable classes.

For a characteristic-zero field `K` and a finite etale `K`-algebra `A`, let
`\mathcal R_K(A)` denote the set of stable polynomial left--right classes of
Keller maps having a complete fiber isomorphic to `Spec A`, as in
[the common-fiber note](COMMON_ARITHMETIC_FIBERS.md#1-the-invariant).

## The theorem

> **Universal Keller-fiber multiplicity theorem.**
>
> If `K` is a characteristic-zero field and `A` is a finite etale
> `K`-algebra of rank `N>=3`, then
>    \[
>      \boxed{|\mathcal R_K(A)|=\infty.}
>    \]

In rank three the representatives are determinant-one fiber-invisible cubic
gauge lifts.  In every rank at least four, one common power shift on all
higher decorations gives the representatives and separates them by Fitting
Newton area.  Thus the result concerns arbitrary abstract finite etale
algebras, not merely a Hilbert family or a specially selected polynomial.

The rank-three assertion is proved in the
[fiber-invisible cubic gauge note](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md).
The uniform assertion in all ranks at least four is proved in the
[all-degree power-shifted gauge note](UNIVERSAL_POWER_SHIFTED_GAUGE_MULTIPLICITY.md).
Its first case is the
[power-shifted quartic gauge note](UNIVERSAL_QUARTIC_GAUGE_MULTIPLICITY.md).
The earlier
[quartic trace-chord note](UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md)
gives a smaller-degree weighted family over number fields, and more
generally whenever its trace quadric is isotropic.
Sections 1--4 below retain the earlier translated minimal-diagonal argument
for ranks at least five.  Sections 5--6 record the two low-rank discoveries.
The complete failure-mode review is in the
[adversarial audit](UNIVERSAL_MULTIPLICITY_ADVERSARIAL_AUDIT.md), and
[three connected witness cards](UNIVERSAL_MULTIPLICITY_WITNESS_CARDS.md)
make degrees four, five, and six fully explicit.

## 1. A centered primitive generator

Assume `N>=5`.  Put

\[
 A_0=\ker(\operatorname{Tr}_{A/K}).
\]

The primitive elements form a nonempty Zariski open in `A_0`: start with
any primitive element and subtract its trace divided by `N`.  The trace
pairing is nondegenerate, and its restriction to `A_0=1^\perp` is
nondegenerate because `Tr(1)=N`.  Since `K` is infinite, there is therefore
a primitive `eta in A_0` with

\[
 \operatorname{Tr}(\eta^2)\ne0.                       \tag{1.1}
\]

Write its monic characteristic polynomial as

\[
 P(T)=T^N+c_{N-2}T^{N-2}+\cdots+c_0.
                                                               \tag{1.2}
\]

Newton's identity gives

\[
 \boxed{c_{N-2}=-\frac12\operatorname{Tr}(\eta^2)\ne0.} \tag{1.3}
\]

For `s in K`, translate the generator by `T=s+S` and define

\[
 G_s(S)=P(s+S)-P(s)=\sum_{j=1}^Ng_j(s)S^j,
 \qquad
 g_j(s)=\frac{P^{(j)}(s)}{j!}.                         \tag{1.4}
\]

Each `g_j` is a nonzero polynomial.  Hence outside a finite set of `s`,

\[
 g_1(s)g_3(s)\cdots g_N(s)\ne0,                        \tag{1.5}
\]

which places the normalized seed on the clean quadratic-gauge coefficient
torus.  The coefficient `g_2` is irrelevant because a polynomial target
shear removes it.

## 2. Every translation has the same complete fiber

Let `F_s` be the quadratic-gauge Keller map attached to `G_s`.  At the
target

\[
 y_s=\left(1,0,-\frac{2P(s)}{g_1(s)}\right),           \tag{2.1}
\]

the inverse equation is

\[
 G_s(S)-\frac{g_1(s)}2
       \left(-\frac{2P(s)}{g_1(s)}\right)
 =P(s+S).                                              \tag{2.2}
\]

The root-engineered reconstruction theorem therefore identifies the
complete scheme-theoretic fiber as

\[
 F_s^{-1}(y_s)
 \simeq\operatorname{Spec}K[S]/(P(s+S))
 \simeq\operatorname{Spec}A.                           \tag{2.3}
\]

The displayed maps have determinant `-2`; one fixed target scaling makes
them determinant one without changing the fiber or the stable class.

It remains only to prove that translation moves through infinitely many
stable classes.

## 3. The quintic invariant

For `N=5`, after normalizing by `g_1`, put `a_j=g_j/g_1`.  The exact
[quadratic-gauge stable-moduli theorem](QUADRATIC_GAUGE_STABLE_MODULI.md)
assigns weights

\[
 w_3=(-2,-1),\qquad w_4=(-3,-4),\qquad w_5=(-4,-5).
\]

Their primitive relation is `-w_3-6w_4+5w_5=0`, so

\[
 I(s)=\frac{a_5(s)^5}{a_3(s)a_4(s)^6}
     =\frac{g_5(s)^5g_1(s)^2}{g_3(s)g_4(s)^6}          \tag{3.1}
\]

is a stable invariant.  From (1.2),

\[
 g_3(s)=10s^2+c_3,\qquad g_4(s)=5s,\qquad g_5(s)=1.
\]

Here `c_3=c_{N-2}` is nonzero, while `g_1` has order at most two at
`s=0`.  Thus (3.1) has a pole at zero and is nonconstant.  This is the
detailed argument of the
[quintic multiplicity note](UNIVERSAL_QUINTIC_FIBER_MULTIPLICITY.md).

## 4. One invariant for every degree `N>=6`

For every index `j>=4`, the normalized quadratic-gauge coefficient
`a_j=g_j/g_1` has stable weight

\[
 w_j=(1-j,-j).                                         \tag{4.1}
\]

The top three weights obey the universal second-difference relation

\[
 w_{N-2}+w_N-2w_{N-1}=0.                               \tag{4.2}
\]

Since `N>=6`, all three indices in (4.2) are at least four, so

\[
 J_N(s)
 =\frac{a_{N-2}(s)a_N(s)}{a_{N-1}(s)^2}
 =\frac{g_{N-2}(s)g_N(s)}{g_{N-1}(s)^2}                \tag{4.3}
\]

is a stable invariant on the clean coefficient torus.

Only the two highest terms in (1.2) contribute to these derivative jets:

\[
 g_N(s)=1,\qquad
 g_{N-1}(s)=Ns,\qquad
 g_{N-2}(s)=\binom N2s^2+c_{N-2}.
\]

Consequently

\[
 \boxed{
 J_N(s)=\frac{N-1}{2N}
        +\frac{c_{N-2}}{N^2s^2}.
 }                                                       \tag{4.4}
\]

Equation (1.3) makes this rational function nonconstant.  A nonconstant
rational function on the affine line over an infinite field has infinite
image.  Removing the finite exceptional set (1.5) still leaves infinitely
many values of `J_N`.

In the first uniform degree `N=6`, no lower-coefficient path is hidden by
this top-jet calculation.  For

\[
 P(T)=T^6+c_4T^4+c_3T^3+c_2T^2+c_1T+c_0,
\]

the relevant jets are

\[
\begin{aligned}
g_1&=6s^5+4c_4s^3+3c_3s^2+2c_2s+c_1,\\
g_3&=20s^3+4c_4s+c_3,\qquad
g_4=15s^2+c_4,\qquad g_5=6s,\qquad g_6=1.
\end{aligned}
\]

Every clean-locus failure is therefore among finitely many roots of explicit
nonzero polynomials.  Moreover

\[
 J_6(s)=\frac5{12}+\frac{c_4}{36s^2},
\qquad
 J_6'(s)=-\frac{c_4}{18s^3}.
\]

The sole constant-`J_6` path is `c_4=0`, equivalently
`Tr(eta^2)=0`; the generator choice in Section 1 excludes it.  This closes
all sextic translation and lower-coefficient exceptions.

If two maps in this family were stably polynomially left--right equivalent
over `K`, they would remain so over an algebraic closure.  The exact
stable-moduli theorem would then force their invariant values to agree.
Thus infinitely many values of `I` in rank five, or `J_N` in ranks at least
six, give infinitely many stable classes, all with the common complete
fiber (2.3).

## 5. Rank four: the power-shifted gauge

For diagonal quadratic-gauge maps of degree four the coarse stable
coefficient quotient is a point, so no coefficient invariant can detect
translation.  The weighted construction turns its quartic presentation
condition into the five-variable trace quadric

\[
 \operatorname{Tr}(\eta^2)=2e^2+4u^2.
\]

Over a number field this quadric is locally isotropic at every place:
it is indefinite at every real place, automatic at complex places, and
five-dimensional over nonarchimedean local fields of `u`-invariant four.
Hasse--Minkowski supplies a `K`-point.  Rationality and weighted
selected-root Torelli then produce infinitely many stable classes of small
weighted maps.

Over a general characteristic-zero field the trace-chord quadric can be
anisotropic.  The replacement mechanism keeps the selected quartic inverse
polynomial fixed at `P=1` and changes its lift

\[
 g_4P^4S^4\longmapsto g_4P^{m+4}S^4,\qquad m\ge0.
\]

Every lift is a polynomial determinant-`-2` map of geometric degree four
with the same complete fiber.  On the normalized ramified stratum its
relative Fitting generator has Laurent support

\[
 \{(0,0),(1,2),(m+4,3)\}.
\]

The generated affine lattice has index `2m+5`, which is preserved by stable
polynomial left--right equivalence.  These indices give infinitely many
classes without a rational-point or trace-form hypothesis.

More generally, applying one common extra power `P^m` to every decoration
of degree at least four gives Fitting support

\[
 \{(0,0),(1,2)\}\cup
 \{(j+m,j-1):4\le j\le N\}.
\]

Its Newton polygon has normalized area

\[
 2N-3+(N-2)m.
\]

This strictly increasing stable invariant gives one mechanism for every
rank `N>=4`; the quartic lattice index is its triangular first case.

The weighted, cancellation, and minimal diagonal quadratic-gauge mechanisms
still collapse in rank three, as proved in
[the boundary note](LOW_RANK_MULTIPLICITY_BOUNDARIES.md).  The next section
explains how changing the cubic lift escapes that scoped collapse.

## 6. Rank three: fiber-invisible cubic lifts

For a translated cubic

\[
 G(S)=g_1S+g_2S^2+g_3S^3,
\]

replace its minimal lift by

\[
 G_{P,n}(S)
 =g_1S+g_2PS^2+
   g_3P(1+P^{n-1}-P^2)S^3,\qquad n\ge4.
\]

At `P=1` the extra factor is one, so every map has the same selected
inverse polynomial and the same complete cubic fiber.  The paired slope and
intercept corrections remain polynomial and preserve determinant `-2`.

The degree-drop polynomial

\[
 h_n(P)=1+P^{n-1}-P^2
\]

has exactly `n-1` simple nonzero geometric roots.  Over each root, two
inverse sheets remain affine and one unramified sheet escapes to the
canonical boundary.  Over `P=0`, the two `q=0` sheets and the `t=0` sheet
are all affine.  The remaining boundary image is the irreducible ramified
discriminant.  Hence the complete canonical boundary has exactly `n`
geometric target components.  Stable normalization functoriality preserves
this count, so different values of `n` give different stable classes.

## 7. Exact regressions

Run

```bash
.venv/bin/python scripts/verify_universal_quartic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_quartic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_power_shifted_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_quintic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_higher_degree_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_multiplicity_witness_cards.py
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
```

The first checker verifies the quartic tangent-chord and normalized weighted
seed identities.  The second checks the exceptional quintic weights and
pole.  The third proves symbolically that the top-three weights satisfy
(4.2) for symbolic `N`, derives (4.4), and checks degrees six through twelve
as concrete regressions.  The arithmetic isotropy and stable-separation
steps are written proofs, not bounded searches.
