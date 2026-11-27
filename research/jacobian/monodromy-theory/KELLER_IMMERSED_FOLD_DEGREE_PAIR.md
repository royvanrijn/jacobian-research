# Degree pairs at an immersed six-sheet polar fold

**Proposition SixSheetImmersedFoldDegreePair.** Suppose a six-sheet
complex plane Keller map has a sole ramified target curve \(B\),
with transposition meridian and every affine analytic branch
immersed. Suppose that a prime above \(B\) is an exceptional
\((2,1,5)\) polar fold, with normalized polynomial source
coordinates \(a,b\) and whole chart \(a=b^{-2}+bt\). Then the
total degree of \(B\) in any polynomial target frame is at
least ten and
\[
                 \max(\deg_a P,\deg_a Q)\ge60
\]
in the normalized source frame. At equality an affine target
change gives restriction degrees \((8,10)\) on the fold line
and longitudinal degrees \((48,60)\). With \(u=ab^2\) the two
lowest blocks are nonzero scalar multiples of
\[
                  b^{-24}u^{40}(u-1)^8,\qquad
                  b^{-30}u^{50}(u-1)^{10}.
\]

## The coprime obstruction

The
[immersed sole-branch projection theorem](KELLER_IMMERSED_SOLE_BRANCH_PROJECTION.md#six-sheet-consequence)
gives composite restriction degrees at least eight for every
target coordinate in this setting. We also use the coprime
step of its proof directly. In any target frame with
restriction degrees \(L<M\) on the branch normalization,
if \(\gcd(L,M)=1\), choose the coordinate of degree \(M\)
as \(x\) and a local parameter \(w\) at infinity with
\(x=w^{-M}\). The companion's first pole term has order \(L\).
Its Puiseux exponent has reduced denominator \(M\), so the
knot at infinity has just one characteristic pair \((L,M)\).
In the conventions of Neumann–Norbury the terminal number
is \(d_1=L+M>1\), by §4 and equation (7). Since this is an
immersed rational component with one point at infinity, their
Lemma 5.4 gives an abelian affine complement group. Its
Orevkov invariant would thus be abelian, contrary to the
[actual-source obstruction](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#local-degree-and-the-finite-branch-quotient).
The hypotheses of the one-pair calculation and of Lemma 5.4
are stated in the
[Neumann–Norbury paper](sources/neumann-norbury-orevkov-invariant.pdf).

The total degree of a polynomially parametrized plane curve
with birational normalization equals the maximum restriction
degree of the two target coordinates, by intersecting with
a generic affine line. An affine target change can make the
two restriction degrees distinct by canceling a common largest
term. Thus total degree at most nine would force the pair
\((8,9)\), which is coprime and impossible. At total degree ten
the candidate \((9,10)\) is coprime too, so only \((8,10)\)
remains.

## Equality in the original polynomial source

By the
[whole-fold-chart theorem](../../hessian/hc4-triangular-graphs/EXCEPTIONAL_POLAR_FOLD_WHOLE_CHARTS_AND_BIRATIONAL_IMAGES.md),
the restriction on this fold line is the branch normalization.
The
[coordinatewise weighted fold bound](../../hessian/hc4-triangular-graphs/EXCEPTIONAL_FOLD_WEIGHTED_JACOBIAN_OBSTRUCTION.md#2-applying-the-equation-on-the-whole-fold-chart)
therefore gives the lower bound sixty in the original source
frame, since some restriction degree is at least ten.
Assume equality. The original restriction degrees cannot
exceed ten. By the preceding argument and a constant affine
change we may normalize the restriction degrees of \(P,Q\)
to eight and ten, respectively, keeping both \(a\)-degrees
at most sixty.

In each target coordinate, if the lowest block is \(b^{-m}p(u)\)
for \(m>0\), the weighted fold bound gives
\(\deg p\ge2m\), with equality only when \(p\) has precisely
the two roots \(0,1\). Indeed its companion in the weighted
Jacobian equation has degree two at equality and vanishes
at all the roots of \(p\). The degree-ten restriction therefore
forces the lowest block of \(Q\) to have weight \(-30\) and
coefficient degree sixty. Its coefficient has order ten at
\(u=1\) by the whole-chart restriction rule, so it is
proportional to \(u^{50}(u-1)^{10}\).

For the lowest block \(b^{-m}p(u)\) of \(P\) we have
\(24\le m\le30\). The vanishing bracket of the two lowest
blocks gives multiplicities proportional to \(m:30\).
In particular at \(u=1\) the order of \(p\) is \(m/3\).
Thus \(3\mid m\). For \(m=27,30\), that block would directly
contribute a nonzero \(t^{m/3}\) term on the fold line,
contrary to restriction degree eight. The contribution cannot
cancel with a different weight. It follows that \(m=24\)
and \(p\) is proportional to \(u^{40}(u-1)^8\).

Write \(d=\deg_a P\) and \(A(b),C(b)\) for the leading
\(a\)-coefficients of \(P,Q\), respectively. Then \(d\ge48\),
\(\deg_a Q=60\), and \(\operatorname{ord}_b C=90\).
The coefficient of \(a^{d+59}\) in the constant Jacobian
vanishes and gives \(d A C'-60 A' C=0\), so
\(\operatorname{ord}_b A=3d/2\). If \(d>48\), every term
of \(P\) with \(u\)-degree \(d\) has weight at least \(-23\),
since its weight-\(-24\) coefficient has degree forty-eight.
Hence \(\operatorname{ord}_b A\ge2d-23>3d/2\), impossible.
This proves \(d=48\). The argument imposes the displayed
necessary form and leaves its polynomial admission open.
