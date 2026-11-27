# Paying for adverse radial directions with same-label angular loss

The sharp aligned-arc bound has an unconditional variance reserve
\(17\rho\pi(1-\rho^2)/1000\).  This is more than five times the
first-feature reserve in [SV1](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md).
It combines the quartic difference remainder with a stronger positive
remainder when both circular slices favor the same label.  The reserve
holds at every amplitude and every \(0<\rho<1\).  It can therefore
repair a positive direction of the uncorrected scalar operator without
assuming closeness to the simplex or continuing a family of maximizers.

The corrected operator has a directed norm bound below two throughout
\([29/50,61/100]\). Together with the preceding sharp scalar comparison,
this proves the unrestricted direct comparison for every
\(0<\rho\le61/100\), in every dimension at least two.

The all-correlation comparison still requires a sign proof for the
corrected operator, or further angular information.  This note does not
assert that its corrected scalar criterion holds at every correlation.

## SLV1: The positive Taylor slack

Use \(k_a(v)=e^{a\cos(\pi v)}\), its even second primitive \(G_a\),
and \(P_\sharp,M,I_0\) from
[SP1--SP4](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md).
For \(s\in[0,2]^3\) with \(\sum_i s_i=2\), set

\[
 H_s=\sum_i(s_i-2/3)^2,\qquad
 D_a(s)=\frac{H_s}{2}P_\sharp(a)
             -\sum_iG_a(s_i)+3G_a(2/3).
\]

This is nonnegative by SP2.  Its contribution to the angular pairing
slack is \(D_a(s)/2\).  At the pure same-label sum \(s_*=(2,0,0)\),

\[
 D_*(a):=D_a(s_*)=M(a)+P_\sharp(a)-2I_0(a).                 \tag{SLV1}
\]

Indeed \(G_a(2)=2I_0(a)\), \(2G_a(1)=M(a)\) and
\(P_\sharp=6G_a(1)-9G_a(2/3)\).

The aligned-arc reduction is the same one used in the preceding
all-mass comparisons.  For primary context, Heilman's
[Theorem A.1 and Appendix A.1, arXiv:2608.00333v1](https://arxiv.org/html/2608.00333v1#A1)
give a balanced comparison through \(2/5\) and an all-parameter
aligned-arc estimate.  The theorem and its parameter restriction were
inspected.  The sharp coefficient and the additional remainder used here
are the internal results SP1--SP4 and SR2, rather than a claim that the
primary theorem already has the larger scope.

## SLV2: A cumulative-density comparison for a dominant label

For every \(a\ge0\), if \(s_1=S\ge4/3\), then

\[
                         D_a(s)\ge(S-1)D_*(a).             \tag{SLV2}
\]

Here is an explicit proof, including the part where pointwise density
comparison would be insufficient.  Write \(h=2-S\le2/3\),
\(s_2=r,s_3=t\), so \(r+t=h\), and put \(b=2/3\).
The probability density for \(P_\sharp\) is \(3v\) on \([0,b]\)
and \(6(1-v)\) on \([b,1]\).  Write its cumulative distribution as
\(R(v)=3v^2/2\) and \(R(v)=1-3(1-v)^2\), respectively.
After reflecting the integrals defining \(G_a(u)\) at 1, the
signed density defining \(D_a(s)\) is

\[
 \frac{H_s}{2}R'(v)
 -\sum_i\big[(s_i-v)_++(s_i-2+v)_+\big]+3(b-v)_+.
\]

Its integral is zero.  Let \(F_s\) be its cumulative integral.
For the pure sum, the cumulative integral is

\[
 F_*(v)=\begin{cases}
             v^2/2,&0\le v\le b,\\
             6v-4v^2-2,&b\le v\le1.
           \end{cases}                                    \tag{SLV3}
\]

We claim \(F_s\ge(1-h)F_*\).  For \(v\le b\), twice this
difference is

\[
 (3-5h+3h^2-3rt)v^2
       -(v-r)_+^2-(v-t)_+^2-(v-h)_+^2.
\]

For each fixed \(c\ge0\), \((v-c)_+^2/v^2\) increases with
\(v>0\).  Evaluate the three quotients at \(v=b\), where
\(r,t,h\le b\).  Subtracting their sum from the coefficient gives
\(h-3h^2/2+3rt/2\ge0\).  The endpoint \(v=0\) is immediate.
For \(v\ge b\), the difference is exactly

\[
              (1-v)^2(2h-3h^2+3rt)\ge0.
\]

Integration by parts, with zero endpoint terms, now gives

\[
 D_a(s)-(1-h)D_*(a)
   =a\pi\int_0^1[F_s(v)-(1-h)F_*(v)]
              \sin(\pi v)e^{a\cos(\pi v)}dv\ge0.
\]

This proves SLV2.  More generally the same argument works for an
increasing differentiable function of \(\cos(\pi v)\).

## SLV3: A uniform bound for the pure same-label slack

For all \(a\ge0\),

\[
                         D_*(a)\ge\frac{19}{100}a.         \tag{SLV4}
\]

By SLV3 the function

\[
 F(a)=D_*(a)/a
   =\pi\int_0^1F_*(v)\sin(\pi v)e^{a\cos(\pi v)}dv
\]

extends continuously to zero and is convex on \([0,\infty)\).
Its integration weight is nonnegative.  At \(a_0=3/2\), the
directed certificate proves

\[
 F'(a_0)>0.0004618631,\qquad
 F(a_0)-a_0F'(a_0)>0.1945673542>19/100.
\]

The supporting affine function therefore exceeds \(19/100\)
on the entire half-line, proving SLV4.  This is an infinite-parameter
consequence of convexity and a single finite certificate.

For reproducibility, let \(M_n,P_n,U_n\) be the cosine-power moments
of the probability measures defining \(M,P_\sharp,I_0\), and
\(d_n=M_n+P_n-2U_n\).  Here \(U_n=0\) for odd \(n\) and
\(U_n=2^{-n}\binom n{n/2}\) for even \(n\).  Then

\[
 F(a)=\sum_{n\ge1}\frac{d_n a^{n-1}}{n!},\qquad |d_n|\le4.
\]

The checker evaluates the series and its derivative through \(n=40\)
in 140-bit Arb arithmetic.  At \(a_0\) the omitted absolute values
are bounded by

\[
 \frac{4a_0^{40}}{41![1-a_0/42]},\qquad
 \frac{160a_0^{39}}{41![1-41a_0/(40\cdot42)]},
\]

respectively.  It also verifies \(F(a_0)<1/5\), rejecting that
larger proposed uniform constant at an actual parameter value.

## SLV4: An unconditional angular variance reserve

Let \(p,q\) be independent draws from any probability distribution
on \(\Delta_3\), and write \(Q=\mathbb E|p|^2\) and
\(V=\operatorname{Var}(p)\).  For
\(u_i(p)=(p_i-2/3)_+\), SLV2--SLV4 imply

\[
 \frac12D_a(p+q)\ge
      \frac{171}{200}a\sum_i u_i(p)u_i(q).                 \tag{SLV5}
\]

If one product is positive, it is the only such label and
\(p_i+q_i\ge4/3\).  With \(x=p_i-2/3,y=q_i-2/3\in[0,1/3]\),
\(1/3+x+y\ge9xy\); the bilinear expression has nonnegative
values at all four corners of its rectangle.  These facts prove SLV5.

The elementary purity estimate

\[
 |p|^2\le5/9+\frac43\sum_i u_i(p)
\]

follows by maximizing the square sum with every coordinate at most
\(2/3\), or, for a dominant coordinate \(m\), from
\(m^2+(1-m)^2\le5/9+(4/3)(m-2/3)\).
Consequently
\(\sum_i(\mathbb E u_i)^2\ge(3/16)(Q-5/9)_+^2\).
Also \(V\le2/3\).

The separate quartic difference slack in SR2 is at least
\(a\pi^{-2}|p-q|^2(1-|p-q|^2/2)\).  The expectation of its
last factor is at least \(V(1-Q)\), by SR7 or QE2.
Thus the sum of the two distinct Taylor slacks has expectation at least

\[
 aV\left[\frac{1-Q}{\pi^2}
            +\frac{1539}{6400}(Q-5/9)_+^2\right]
                         \ge\frac{17}{500}aV.             \tag{SLV6}
\]

For the last inequality use \(\pi^2<10\).  If \(Q\le5/9\),
the first term is at least \(2/45\).  Otherwise completing the
square gives the lower bound
\(2/45-16/1539>17/500\).  No balance, small-amplitude or
regularity assumption on this distribution is required.

Now put \(d=1-\rho^2\), \(\kappa=\rho^2/d\),
\(d\nu(t)=e^{-t}dt\), \(f(t)=\sqrt t\,e^{-\kappa t}\),
\(Z=\sqrt\pi d^{3/2}/2\), and \(d\mu=f\,d\nu/Z\).
In the planar polar pairing, \(a=2\rho\sqrt{ts}/d\), and its
remaining radial factor is \(d^{-1}e^{-\kappa(t+s)}\).
The extracted slacks in SLV5 and SR2 both have a constant lower
coefficient times \(a\); hence they integrate against the same
product probability measure \(\mu\otimes\mu\).  Applying SLV6 yields

\[
 S_\rho(A)-S_\rho(T)
 \le \sum_iq_\rho(p_i-1/3)
       -\beta_\rho\operatorname{Var}_{\mu}(p),\qquad
                  \beta_\rho=\frac{17\rho\pi d}{1000}.     \tag{SLV7}
\]

This is a global bound for arbitrary circular masses and angular shapes.
It uses the positive and negative Taylor slacks once each.

## SLV5: The corrected radial operator and its slab comparison

The reserve SLV7 replaces the sharp scalar kernel and weight by

\[
 K^\#=K+\frac{\beta_\rho}{Z^2}f\otimes f,\qquad
 J^\#=J+\frac{\beta_\rho}{2Z}f.                           \tag{SLV8}
\]

Then the right side of SLV7 is exactly
\(\sum_i[\langle z_i,K^\#z_i\rangle-2\int J^\#z_i^2d\nu]\).
In the GH5 positive series only the first row coefficient and the
first kernel-feature coefficient change:
\(b_1\mapsto b_1+2d\beta_\rho/Z\) and
\(k_1\mapsto k_1+\beta_\rho/Z^2\).

For a rational slab \([l,u]\subset[29/50,61/100]\), form these
coefficients at \(u\).  Multiply the \(n\)th row coefficient by
\((l/u)^n\), and call the truncated weight \(J_{\rm mix}\).
Use 32 row terms and 161 kernel terms.  The fixed-feature comparison
of [QE3](THREE_CELL_QUARTIC_EXTENSION_THROUGH_THREE_FIFTHS.md)
applies also to the corrected coefficients: in the scaled
coordinate \(r=\sqrt{t/d}\), both new first coefficients are linear
in \(\rho\), with positive constants.  More explicitly the added
row polynomial coefficient is \(17\rho\sqrt\pi/250\) and the
added kernel-form coefficient is \(68\rho/125\).
Thus the mixed row at \(l\) is a lower bound for the row, and the
anchored kernel at \(u\) is an upper bound throughout the slab.
Indeed the auxiliary row norm in the scaled variable has density
\((d/2)r^2e^{-r^2}Q_l(r)dr\). Writing the norm variable as
\(y=\sqrt{d/2}\,r e^{-r^2/2}\sqrt{Q_l(r)}z\), the original finite
kernel becomes a sum of the rank-one forms with features
\(e^{-r^2/2}r^n/\sqrt{Q_l(r)}\) and coefficients
\(4(2\rho)^n(M_n+P_n)/n!\). The added coefficient just displayed
is also increasing with \(\rho\). This proves the comparison directly.

At the anchor use the centered Gram matrix and positive coefficients
as in GH6.  Its weighted absolute row-sum bound is denoted by
\(\lambda\).  The physical mean-zero hyperplane differs from the
anchored one.  Put

\[
 W=\int J_{\rm mix}^{-1}d\nu,\quad
 w=(WJ_{\rm mix})^{-1/2},\quad P=I-w\otimes w,
 \quad L=\frac{u^2-l^2}{1-u^2}
          \sqrt{\frac{\int t^2/J_{\rm mix}\,d\nu}{W}}.
\]

The transport \(\widehat z(t_0)=\sqrt{d/d_u}z(t)\),
\(t_0=d_u t/d\), turns the physical mean-zero condition into
orthogonality to \(e^{-(d/d_u-1)t_0}/\sqrt{J_{\rm mix}}\).
The inequality \(1-e^{-x}\le x\) therefore bounds
\(|\langle v,w\rangle|\le L\|v\|\).
If \(B\) is the anchored uncentered finite operator,
\(b_w=\|PBw\|\) and \(k_w=\langle w,Bw\rangle\), its norm
on transported mean-zero profiles is at most

\[
                   \lambda+2b_wL+k_wL^2.                  \tag{SLV9}
\]

The quantities \(b_w,k_w\) are computed from the same moments,
as in QE4.  For an inexpensive certified upper bound on the leakage
moment, write \(x=\sqrt{t_0}\), \(\eta=1-\kappa_u>0\), and
\(Q(x)=\sum_{n=1}^{32}b_n^{\rm mix}x^{n-1}\).  Then

\[
 \int t_0^2/J_{\rm mix}\,d\nu
   =8d_u\int_0^\infty\frac{x^4e^{-\eta x^2}}{Q(x)}dx
          \le\frac{3d_u\sqrt\pi}{b_1^{\rm mix}\eta^{5/2}}.
\]

SP5 bounds the omitted positive feature terms against the full
uncorrected row, hence also the full corrected row.  Denote that
tail by \(\tau\).  Therefore

\[
               \lambda+2b_wL+k_wL^2+\tau<2                \tag{SLV10}
\]

is sufficient for strict comparison at every nonuniform mass profile
on the entire slab.  This is a linear operator criterion, valid at
all amplitudes, with the adverse scalar direction already paid for.

## SLV6: Directed radial evidence

The [checker](check_three_cell_same_label_variance.py) uses 140-bit
Arb arithmetic, 4096 radial panels on \([0,12]\), 32 row terms,
and kernel terms of degree 0 through 160.  It uses the GH8 linear
interpolation enclosure with its second-derivative error.  For each
Gaussian moment tail it uses the full positive polynomial at the
cutoff: \(Q(x)\ge Q(12)\) for \(x\ge12\).  In particular, for
\(\eta>(n-1)_+/(2R^2)\),

\[
 8d_u\int_R^\infty\frac{x^ne^{-\eta x^2}}{Q(x)}dx
 \le\frac{4d_u R^{n-1}e^{-\eta R^2}}
              {Q(R)[\eta-(n-1)_+/(2R^2)]}.
\]

All coefficients of \(Q\) are positive.  This improves the older
first-coefficient-only enclosure, which can manufacture large
finite-operator bounds when the kernel truncation order is increased.
The improved enclosure keeps the entire infinite radial domain.

The contiguous directed cover of \([29/50,61/100]\) passes on 40
accepted rational slabs, of widths \(1/1000,1/2000,1/4000,1/8000\).
The largest corrected norm upper bound, including leakage and the
infinite feature tail, is less than \(1.998842<1999/1000\).
The cover starts with the 30 consecutive slabs of width \(1/1000\)
and bisects only when the displayed inequality cannot yet be proved.
Every accepted slab is checked on its full parameter interval; no
sampling interpolation in correlation is used. The separate point
check at \(\rho=61/100\) gives an upper bound below \(1.994114\).
The [saved receipt](three_cell_same_label_variance_verification.json)
pins the source and imported checker helpers. All computations reported
here are internal directed verification.

## SLV7: Unrestricted direct comparison through 61/100

For every ordinary measurable balanced Gaussian three-cell partition
\(A\) in \(\mathbb R^n\), \(n\ge2\), and every
\(0<\rho\le61/100\),

\[
 S_\rho(A)\le S_\rho(T),\qquad
 \widetilde Q_\rho(1_A-1_T)\le2L_\rho(A;T).                \tag{SLV11}
\]

The interval through \(29/50\) is the preceding
[sharp scalar theorem](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md).
On the new interval SLV7--SLV10 and the contiguous directed cover
give strict planar comparison whenever the circular masses are
nonuniform. Uniform circular masses have the simplex as their only
equality case by DG6. Thus equality on the new interval requires
a simplex up to rotation, relabeling and null sets.

Existence and dimension reduction of a global maximizer are the same
imported Heilman--Tarter Lemma 2.3 and Theorem 1.9 used in GH1 and
DG4. Their hypotheses concern ordinary balanced Gaussian partitions
and apply here without alteration. These inputs reduce a maximizing
competitor to a plane; they do not transform an arbitrary competitor
into a planar partition. The planar comparison therefore implies the
first inequality in every stated dimension. The exact expansion B1
then gives the requested second inequality, as in DG4.
This adoption is not an independent replay of those imported results.

There is no continuation assumption in this proof and no smallness
hypothesis on the competitor. The interval \((61/100,1)\) is not
settled by this certificate; the separate qualitative near-one result
does not establish an overlap closing that gap.
