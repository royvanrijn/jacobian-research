# Zero separation for the half-line cubic on the extended temperature interval

Let
\[
 F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy .
\]
At a simple zero, \(H_T\) denotes differentiation of its local zero
branch in temperature.

**Theorem.** For every real \(-2\le T\le2\), the zeros of
\(F(T,\cdot)\) are simple, distinct zeros have distinct imaginary
heights, and \(\Im(2H-T H_T)>0\) at every upper zero.
This extends the [central interval certificate](MULTIVARIATE_HALF_LINE_CUBIC_BOUNDED_TEMPERATURES.md)
to \([-2,2]\).

## Explicit exterior comparison

The decreasing-density sine-transform argument excludes
\(\Re H\le-1\) for \(|T|\le2\): there the derivative of the
log density is at most \(-1+2y-y^2\le0\).
There are no real zeros by positivity. Consider henceforth
\(\Re H\ge-1\), \(\Im H\ge0\), and \(R=|H|\ge12\).
Set \(a=T/2\), \(Z=H+a^2\) and \(\omega=e^{2\pi i/3}\).
The [exact Airy identity](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#an-exact-airy-decomposition)
writes \(F=E+S\) with
\[
 E=\omega F(T\omega^2,\omega H),\qquad
 S=2\pi i\omega^2e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2 Z).
\]
Put \(d=-\Re(\omega H)\ge R/2\). The remaining exponent
on the endpoint ray has real part at most
\(\max_{y\ge0}(-y^3/3+y^2/2)=1/6\). Thus
\[
 |E+1/H|\le\tfrac65(2d^{-4}+2d^{-3}),\qquad
 |E_H|\le\tfrac65d^{-2},\qquad |E_T|\le\tfrac65d^{-3}.
\]
Here the first estimate compares with the linear integral,
using \(|e^v-1|\le |v|e^{\max(0,\Re v)}\).
For \(R\ge12\) these imply
\(5/(6R)\le |E|\le7/(6R)\) and \(|E_H/E|\le6/R\).

The [explicit Airy-sector estimate](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4)
gives
\[
 \left|\omega^2\frac{\operatorname{Ai}'(\omega^2 Z)}
                        {\operatorname{Ai}(\omega^2 Z)}
                      -\sqrt Z\right|\le3/|Z|.
\]
Indeed \(|Z|\ge R-1\) and the Airy argument stays in
\(|\arg(\omega^2Z)|\le2\pi/3\).
For \(\Xi=\log(S/E)\), interpreted locally, we have
\[
 \Re\Xi_H\ge\sqrt{(R-2)/2}-1-3/(R-1)-6/R>0.
\]
At a zero this proves simplicity. On a fixed height
\(\Im H\ge12\), the entire segment between two proposed
zeros is exterior, so strict increase of \(\log|S/E|\)
excludes a shared height.

At an exterior zero one in fact has \(\arg H>\pi/5\).
If \(0\le\arg H\le\pi/5\), then
\[
 \Re(2Z^{3/2}/3-aH-2a^3/3)
      > R^{3/2}/3-R-2/3>0 .
\]
Here \(\cos(3\pi/10)>1/2\), and the Airy prefactor
already exceeds the endpoint upper bound \(7/(6R)\).
This proves saddle dominance in that cone. The same
argument excludes \(\Re H\ge30\), \(0\le\Im H\le18\)
by using the slightly weaker angle bound \(\arg H<\pi/4\)
and \(R\ge30\).

At an exterior zero \(\arg H>\pi/5\) gives
\(d\ge4R/5\). For nonnegative real part this follows
from \(d/R=\cos(\arg H-\pi/3)\); for negative real part
use \(\Re H\ge-1\) and \(R\ge12\).
Consequently \(|E_H/E|\le9/(4R)\), \(|E_T/E|\le3/R^2\).
With \(P=\sqrt Z-a\), the saddle derivatives are
\(P\) and \(-P^2/2\), up to errors at most
\(3/|Z|\) in each variable. The resulting bounds give
\[
 \begin{split}
 |\Xi_H|&\ge\sqrt{R-1}-1-3/(R-1)-9/(4R)>1.83,\\
 2|\Xi_T|/R&\le 1+2\sqrt{R+1}/R+2/R+
                   6/[R(R-1)]+6/R^3<1.82 .
 \end{split}
\]
The first bound increases, and the second decreases,
from \(R=12\). Thus \(|T\Xi_T/\Xi_H|<R\).
Since \(\Im H>R/2\) in the remaining cone
(the negative real part is at most one in modulus),
\(\Im(2H-T H_T)=2\Im H+\Im(T\Xi_T/\Xi_H)>0\).

## Enclosed low-field calculation

The rectangle
\[
                    [-1,30]\times[0,15+43T/40]
\]
has no zeros on its left, lower, or right side by the
exterior and real-axis estimates. Its upper height is
always greater than twelve. Only the part with
\(3\le\Re H\le15\) needs interval verification on the top:
at \(\Re H=3\) the ratio modulus \(|S/E|\) is below one,
and at \(\Re H=15\) it is above one. Then the strict
horizontal derivative excludes both complementary
pieces. To check these endpoint inequalities, at
\(\Re H=3\) one has
\(\arg Z>2\pi/5\) since \(\Im H\ge257/20\) and
\(\Re Z\le4\). The real saddle exponent is at most
\(-13/3\), using \(\cos(3\pi/5)<-3/10\) and
\(|Z|^{3/2}>40\). The remaining factors in \(|S/E|\)
multiply to less than \(48\), while \(e^{13/3}>48\).
At \(\Re H=15\), set \(x=\Re Z\ge15\).
One has
\(2x-|Z|\ge2x-\sqrt{x^2+18^2}>13/2\) and
\(\Re\sqrt Z>19/5\). Since
\(\Re Z^{3/2}=(\Re\sqrt Z)(2x-|Z|)\),
the real saddle exponent is positive. Its prefactor
already dominates the endpoint as above.

The
[extended-temperature verifier](verify_half_line_cubic_extended_temperatures.py)
uses the same Taylor and interval method as the central
certificate. In the temperature series it takes 220 terms
and uses the coefficients \(m_{j+2l}/l!\) in powers
\(-T/2\), for \(j=0,1\). From term 220 onward the ratios
in each residue class are below \(1/2\), by
\(m_{k+3}=(k+1)m_k\). In the field series it uses
\((R,K)=(20,520),(24,800),(36,1050)\) as the
disk radius and first omitted index. In each case
\(R^3/[K(K-1)]+2R/K<1/4\), which is the same
three-term block bound for the field Taylor tail
with \(|T|\le2\). Derivative tails follow by multiplying
each term by its index divided by \(R\).

Point values use Arb with 512-bit precision. The bound on
\(T\)-evolution of the 1-norm of \((F,F_H,1)\) at fixed
\(H_0\) is \(3|H_0|/2+3\). The subsequent field evolution
bound on a segment of length \(d_H\) is \(|H_0|+d_H+3\).
The same two evolution equations as in the central proof
give the value and derivative variations. The verifier
tightens the second field derivative bound by substituting
the resulting first value and derivative variations back
into \(F_{HH}=HF+1-TF_H\). Each such substitution uses only
bounds already valid on the full cell.

The certified initial winding is seven, with 3512 certified
boundary segments. The upper side is checked with 95340
cells in the \((T,\Re H)\) rectangle.
For every one of the \(2048\) rational temperature slabs,
seven Rouché disks of radius \(1/128\) lie inside the
contour and have disjoint imaginary projections. The replay
certifies all 14336 disks. The field centers are selected by numerical Newton steps
and the final inequalities are all interval tests.
Each disk contains exactly one zero by comparison
with the linear field term at its center. The verifier
also checks \(4\Im H\,|F_H|-2>1\) throughout these disks,
which ensures the positive imaginary order derivative
since \(H_T=1/(2F_H)-T/2\).
Together with the constant outer winding, these disks
exhaust the interior. Exterior zeros have heights
above the top and are covered by the strict horizontal
and imaginary derivative bounds, completing the proof.
\(\square\)

Run the certificate from the repository root with python-flint
installed.

```text
.venv/bin/python research/factorial/image-rigidity/verify_half_line_cubic_extended_temperatures.py
```
