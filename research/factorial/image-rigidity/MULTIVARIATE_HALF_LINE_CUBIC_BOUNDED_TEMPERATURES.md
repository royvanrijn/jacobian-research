# Uniform zero separation for the cubic model on the central temperature interval

Let
\[
 F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy
\]
as in the [moving one-private cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
The notation \(H_T\) denotes the derivative of a simple zero branch.

**Theorem.** For every real \(T\in[-1,1]\), all zeros of
\(F(T,\cdot)\) are simple. Distinct zeros have distinct imaginary
heights. At every upper zero,
\[
                       \Im(2H-T H_T)>0 .
\]
By conjugation the corresponding imaginary part is negative
at every lower zero. This removes the entire interval
\([-1,1]\) from the finite exceptional-temperature set in the
[half-line cubic separation result](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md).
The proof uses explicit exterior estimates and a reproducible
finite interval certificate.

## Uniform exterior derivative bounds

For \(|T|\le1\), all upper zeros lie in \(\Re H>-1/4\).
Indeed at \(\Re H\le-1/4\), the density in the sine transform
is strictly decreasing, since its logarithmic derivative is
\(\Re H-Ty-y^2\le-1/4+y-y^2\le0\).
Pairing adjacent half periods gives a strictly positive imaginary
part. There are no real zeros by positivity of the integral.

For \(\Im H\ge0\), \(\Re H\ge-1/4\), and \(R=|H|\ge8\), write
\(a=T/2\), \(Z=H+a^2\), \(\omega=e^{2\pi i/3}\).
Use the [Airy decomposition](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#an-exact-airy-decomposition)
\[
 F(T,H)=E+S,quad E=\omega F(T\omega^2,\omega H),\quad
 S=2\pi i\omega^2 e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2 Z).
\]
Set \(d=-\Re(\omega H)\). One has \(d\ge R/2\), from
\(d=(\Re H+\sqrt3\Im H)/2\) and the lower bound
\(-1/4\) on the real part. On the real integration ray for
the endpoint term, the real residual exponent is at most
\(\max_{y\ge0}(-y^3/3+y^2/4)=1/48\).
The inequality \(|e^v-1|\le |v|e^{\max(0,\Re v)}\) therefore
gives
\[
 \left|E+\frac1H\right|\le\frac{21}{20}(2d^{-4}+d^{-3}),
 \quad |E_H|\le\frac{21}{20d^2},\quad |E_T|\le\frac{21}{20d^3}.
\]
In particular \(4/(5R)\le|E|\le6/(5R)\),
\(|E_H/E|\le6/R\) and \(|E_T/E|\le11/R^2\).

The Airy argument has principal argument of absolute value
at most \(2\pi/3\), and \(|Z|\ge R-1/4\).
The explicit primary-sector errors for Ai and Ai-prime give
\[
       \left|\omega^2\frac{\operatorname{Ai}'(\omega^2Z)}
                               {\operatorname{Ai}(\omega^2Z)}
                        -\sqrt Z\right|\le\frac3{|Z|}.
\]
Indeed their individual leading relative errors are at most
\((25/24)|Z|^{-3/2}\) and \((49/48)|Z|^{-3/2}\), respectively.
These are the bounds of the
[Airy envelope estimate](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4).
In particular \(S\) is nonzero.

Let \(\Xi\) denote a local logarithm of \(S/E\), and put
\(P=\sqrt Z-a\). The derivative bounds just obtained give
\[
 \begin{split}
 |\Xi_H-P|&\le3/|Z|+6/R,\\
 |\Xi_T+P^2/2|&\le3/(2|Z|)+11/R^2 .
 \end{split}
\]
Since \(\Re P\ge\sqrt{(R-1/2)/2}-1/2\), the horizontal
derivative \(\Re\Xi_H\) is strictly positive for \(R\ge8\).
Thus exterior zeros are simple. In particular at each height
\(\Im H\ge8\) a horizontal segment joining two proposed zeros
would have strictly increasing \(\log|S/E|\), ruling out
the pair.

At an exterior zero the implicit derivative gives
\(H_T=-\Xi_T/\Xi_H\). The estimates above also show that
\(2H\Xi_H+T\Xi_T\) differs from
\(2HP-TP^2/2=P^2(2\sqrt Z+a)\) by less than \(19\).
The latter term has modulus at least \((9/4)^2\cdot5>25\),
because \(\sqrt{|Z|}>11/4\) and \(|a|\le1/2\).
Thus \(2H-T H_T\ne0\) on the exterior.

The sign of its imaginary part also follows from the
moduli of the two terms. At a zero, \(\log|S/E|=0\) confines the argument to
\(\arg H>\pi/4\) when \(R\ge8\). In fact if
\(0\le\arg H\le\pi/4\) one has
\(\Re Z^{3/2}\ge R^{3/2}\cos(3\pi/8)>R^{3/2}/3\).
It follows that
\[
 \Re(2Z^{3/2}/3-aH-2a^3/3)
   >(2/9)R^{3/2}-R/2-1/12>0 .
\]
The prefactor already satisfies
\(\sqrt\pi |Z|^{-1/4}(1-(25/24)|Z|^{-3/2})>6/(5R)\),
so saddle dominance excludes these fields. Thus at a zero
\(\pi/4<\arg H<3\pi/4\) and \(\Im H>R/\sqrt2\).
The derivative bounds give
\[
 |\Xi_H|\ge \sqrt{|Z|}-1/2-3/|Z|-6/R>1,\qquad
 |\Xi_T|\le (\sqrt{|Z|}+1/2)^2/2+3/(2|Z|)+11/R^2<R .
\]
The last upper bound uses
\(R/2+\sqrt{R+1/4}/2+3/4<R\) for \(R\ge8\).
Hence \(2\Im H+\Im(T\Xi_T/\Xi_H)>0\), proving
the exterior imaginary order-derivative bound. The same saddle-dominance argument already
excludes \(\Re H\ge15\), \(0\le\Im H\le10\), since this
region has angle less than \(\pi/4\) and \(R\ge15\).

## Uniform bounded-field certificate

For each \(T\in[-1,1]\), consider the rectangle
\[
             [-1/4,15]\times[0,9+17T/20]
]
in \((\Re H,\Im H)\). Its left, lower, and right sides are
zero-free by the estimates above. The certified upper side
also avoids zeros, and its height is always greater than eight.
The boundary is a continuous deformation of the \(T=0\)
rectangle; hence its winding gives exactly three zeros,
counted with multiplicity, for every \(T\) in the interval.

The
[bounded-temperature verifier](verify_half_line_cubic_bounded_temperatures.py)
does the analytic enclosure and Rouché comparisons described
below. Run from the repository root with python-flint installed.

```text
.venv/bin/python research/factorial/image-rigidity/verify_half_line_cubic_bounded_temperatures.py
```

The replay gives initial winding three with \(1067\)
certified segments, \(14520\) certified upper-boundary cells,
and \(768\) certified root disks. The root disks consist of
three disks for each of \(256\) rational temperature slabs.
At every temperature in a slab, the three disks are disjoint
in their imaginary projections and each contains one zero.
The verifier also certifies throughout every disk
\[
                           4\Im(H)\,|F_H|-1>5.
]
Each single-zero count proves simplicity. Using
\(H_T=1/(2F_H)-T/2\) at a simple zero, the last inequality
implies \(\Im(2H-T H_T)>0\) there.

Here are the enclosure rules. Write
\(m_k=\int_0^\infty y^k e^{-y^3/3}\,dy\); then
\(m_{k+3}=(k+1)m_k\), and the first three moments are
\(3^{-2/3}\Gamma(1/3),3^{-1/3}\Gamma(2/3),1\).
The series for \(F(T,0)\) and \(F_H(T,0)\) have coefficients
\(m_{j+2l}/(2^l l!)\) in powers \((-T)^l\), for
\(j=0,1\). After the first 80 terms the ratios in each of
the three residue classes are less than \(1/2\) for
\(|T|\le1\), by the moment recurrence. Twice the first
omitted block bounds the tail. Now use the field series
of \(F\), with coefficients \(a_k\) from
\(F_{HH}+T F_H-HF=1\). For \(|H|\le20\), the weights
\(|a_k|20^k\) satisfy
\[
 |a_k|20^k\le \frac{20^3}{k(k-1)}|a_{k-3}|20^{k-3}
                       +\frac{20}{k}|a_{k-1}|20^{k-1}.
]
For \(k\ge420\) the sum of the two factors is below
\(1/4\), so successive blocks of three contract at that
rate. This also bounds the derivative tail by summing
the same weights multiplied by \(k/20\).
Point values and derivative enclosures are obtained
with these series in Arb arithmetic.

For interval propagation write \(p=F_H\).
The field evolution satisfies \(p_H=HF+1-Tp\).
The temperature evolution satisfies
\[
 F_T=(Tp-HF-1)/2,\qquad
 p_T=((TH-1)F-(T^2+H)p+T)/2.
]
Thus the \(1\)-norm of \((F,p,1)\) grows with exponential
rate at most \(|H_0|+1\) along temperature segments at
fixed \(H_0\), and at most \(|H_0|+d_H+2\) on a
field segment of length at most \(d_H\) at fixed
\(|T|\le1\). The verifier integrates these inequalities
first along temperature, then along field. To enclose
a top cell it uses the temperature radius \(d_T\)
and field radius \(d_x+(17/20)d_T\), where \(d_x\)
is the horizontal half-width. The center value
exceeds the resulting variation bound, so the
whole cell avoids zero.

On the root slabs the field radius is \(1/64\).
The approximate Newton iterates are used only to
select rational disk centers. With point values and
the two evolution bounds, the verifier encloses
\(F(T,H)-F_H(T_0,H_0)(H-H_0)\) on each disk, for all
temperatures in the slab, by a quantity strictly
less than \(|F_H(T_0,H_0)|/64\). Rouché gives its
one-zero count. Applying the same derivative bound
throughout the disk gives the stated speed test.
Since the enclosing rectangle has count three, the
disks exhaust all zeros with height below the top.

The strict exterior height inequality covers every
zero above the top; mixed-height pairs are separated
by the zero-free top itself. The exterior speed and
the verified bound inside give the stated conclusions.
\(\square\)
