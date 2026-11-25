# Effective exclusion at the negative one-private boundary

Use ordinary profiles with exactly one private row \((0,1)\),
all other unequal rows shared, and the centered fields
\(\xi_j=\log(j^\Delta t)+\sum_i d_i\log a_i\) in the upper
parameter half-plane. Write \(\sigma=\sum_i d_i^2/a_i\).

**Theorem.** There are computable \(\epsilon_*>0\) and an integer
\(N_*\) such that, for \(0\le\sigma\le5/4\) and \(N\ge N_*\),
a consecutive common root cannot have
\(\pi-\epsilon_*\le\Im \xi_j\le\pi\).
This is uniform in the shared rows, dimension, and integer
degree difference. Conjugation gives the lower counterpart.

This gives a quantitative version of the negative-real part of
the [uniform one-private regular comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md)
on the band that meets the near-critical endpoint.
The proof uses only two vertical contours and an exact reflection
on the right contour.

## Uniform gap saddle

Use \(r_i=d_i/a_i\), \(B(q)=\sum_i d_i\log(1+r_iq)\), and
\(V=\sum_i d_i^2/\min(a_i,b_i)\).
As in the [effective positive-boundary proof](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_POSITIVE_ENTRY.md#coefficient-ratios-and-the-negative-real-fields),
\(-5/8\le r_i\le5/4\), \(V\le10/3\), and every actual
root satisfies \(|\Re\xi_j|\le13\).
Allow real fields \(|x|\le14\) in the construction below.

Let \(p_R=\min_{r_i<0}(-1/r_i)\), with \(p_R=\infty\)
when there is no negative gap.
On the real right gap \(1<q<p_R\), \(B\ge0,\ B'\ge0\),
and \(Y(q)=\log(q-1)+B(q)\) strictly increases from
\(-\infty\) to \(+\infty\).
At a finite right pole, the corresponding negative gap term
of \(B\) tends to \(+\infty\); in the infinite case the logarithm
does so. Thus there is exactly one real solution \(s=s(x)\)
of \(Y(s)+x=0\). It satisfies
\[
 \eta\le s-1\le e^{14},\qquad
 1+r_i s\ge e^{-38}\quad(r_i<0),\qquad
                     Y'(s)\ge e^{-14},                    \tag{1}
\]
where \(\eta=e^{-24}\) can be replaced by any smaller
fixed rational when bounding analytic neighborhoods.
The upper bound follows from \(B\ge0\).
At \(q=1+\eta\) the logarithm equals \(-24\) and
\(B(q)<5\), by \(\sum |r_i|\le5/4\) and
\(1+r_i q\ge19/64\) throughout \([0,1+\eta]\).
This proves \(s-1\ge\eta\).
Then \(B(s)=-x-\log(s-1)\le38\), which bounds each of
its nonnegative summands and yields the factor bound.
The derivative bound follows from \(Y'=1/(q-1)+B'\).

Set \(\zeta=\xi-i\pi\) near this boundary and define the right
phase on the right gap and its vertical neighborhood by
\[
 P(q,\zeta)=\zeta q+(q-1)\log(q-1)-q+\int_0^q B(u)\,du.
\]
The integrand defining the primitive is taken in the domain of
positive real parts of all \(1+r_i u\); its straight real
value extends analytically there. The right logarithm has real
value on \(q>1\). The saddle \(q_*(\zeta)\) solves
\(Y(q_*)+\zeta=0\), and its stationary value \(D=P(q_*,\zeta)\)
has \(D'(\zeta)=q_*(\zeta)\).

All derivatives of \(Y\) on a sufficiently small disk around
\(s(x)\) have explicit uniform bounds from (1) and
\(\sum |d_i r_i|=\sigma\le5/4\).
Choose a positive rational disk radius small enough that
\(|Y'(q)-Y'(s)|\le e^{-14}/4\).
Rouché continues the saddle for \(\zeta=x+i\beta\) with
\(|\beta|\) uniformly small. Its displacement is at most
\(2e^{14}|\beta|\), by the same derivative bound.
Shrinking the imaginary cap gives
\(\Re q_*>1+\eta/2\) and
\(\Re Y'(q_*)\ge e^{-14}/2\).
These local analytic saddles glue on the field rectangle
because each agrees with the unique real solution there,
and the uniform small disks can be chosen on a slightly
larger real interval than the root cap.

## The two exact vertical integrals

Let \(z=jq\), and use the branch
\(\arg(-t)=\arg t-\pi\) near negative real \(t\) above the axis.
The normalized moment is the integral on the upward left line
minus the integral on the upward right line, with kernel
\(K_j(z)/(2\pi i)\), where
\[
 K_j(z)=
 \frac{j!\,\Gamma(-z)\Gamma(1+z)}{\Gamma(j-z+1)}
 \prod_i\frac{\Gamma(ja_i+d_i z+1)}{\Gamma(ja_i+1)}
 (-t)^z .
\]
Use the left line \(\Re z=-1/2\) and the right line
\(\Re z=j\,\Re q_*(\zeta)\).
The only intervening poles are the coefficient poles \(0,\ldots,j\),
with residues minus their respective summands.
Indeed the ratio \(\Gamma(-z)/\Gamma(j-z+1)\)
is \(1/\prod_{k=0}^j(k-z)\), and all remaining Gamma
arguments have positive real parts in the rectangle.
Horizontal integrals tend to zero at infinite height.
For fixed \(j\) this follows from the same rational ratio,
the elementary gamma absolute ratio bound below for
\(\Gamma(1+z)\), and the bound by real Gamma values for the
remaining factors. Thus these are exact integrals.

For the left line take \(\theta=\Im\xi\) with
\(|\theta-\pi|\le\pi/4\).
The identity \(\Gamma(-z)\Gamma(1+z)=-\pi/\sin(\pi z)\)
gives decay \(O(e^{-\pi|\Im z|})\) there.
The reciprocal falling Gamma factor, relative to its real
value, costs at most \(e^{\pi|\Im z|/2}\), by the Euler
product. Indeed for \(A\ge1\) the reciprocal modulus
logarithm is at most
\[
 \tfrac12\int_0^\infty\log(1+v^2/s^2)\,ds=\pi|v|/2,
\]
using monotonicity to bound the Euler sum.
The shared factors are bounded in modulus by their real
Gamma values, and the explicit real Stirling ratios after
extracting the field scale are bounded uniformly.
The field factor costs at most \(e^{\pi|\Im z|/4}\).
This gives a uniform integrable exponential majorant.

On \(|\Im z|\le C_0\log j\), the exact remaining Gamma
ratios after extracting \((-e^\xi)^z\) differ from one
by \(O((1+|\Im z|)^2/j)\) with an effective constant.
Indeed compare each \(\log\Gamma(A+w+1)\) to its linear
Taylor term \(w\log A\), using
\(|\psi(A+1)-\log A|\le1/A\) and the
trigamma bound \(O(1/A)\) on \(\Re(A+w)\ge A/2\).
In the row sum
\(\sum_i |d_i|/a_i\le\sigma\) and
\(\sum_i d_i^2/a_i=\sigma\).
The conditions \(|w|\le A/2\) and a small total
Taylor defect follow for large \(j\), with explicit
thresholds using \(\log j/j\) and \((\log j)^2/j\).
The same argument treats the falling Gamma ratio.
Take \(C_0\) large enough that the omitted tails
are \(O(j^{-2})\).
The exact upward left integral thus satisfies
\[
                    E_j(\xi)=(1-e^\xi)^{-1}+O(1/j),         \tag{2}
\]
with computable error. The limiting integral is the
[Mellin–Barnes inversion at the private endpoint](MULTIVARIATE_ONE_PRIVATE_SIDE_RICE_EVENTUAL.md#one-private-linear-variable).
The leading coefficient is bounded and bounded away from
zero throughout the prescribed field rectangle.

## The right line and its uniform error

On the right the exact reflected kernel is
\[
 K_j(z)=(-1)^{j+1}j!\,\Gamma(z-j)
 \prod_i\frac{\Gamma(ja_i+d_i z+1)}{\Gamma(ja_i+1)}(-t)^z .
\]
This follows either from reflection or the rational ratio.
All Gamma arguments have positive real parts there.
The Euler product gives, for \(A>0,\ |v|\ge A\),
\[
 \log \frac{|\Gamma(A+iv)|}{\Gamma(A)}
       \le-\frac12\int_0^{|v|}
               \log(1+v^2/(A+s)^2)\,ds
       \le -g_0|v|,\qquad g_0=\tfrac12\log(5/4).
\]
Use this on \(\Gamma(z-j)\), and bound each shared
numerator by its real Gamma value.
For \(|\beta|\le g_0/2\), the field factor on a vertical line
costs at most \(e^{g_0|\Im z|/2}\).
Hence the tails \(|\Im q|\ge Y_0\), for any fixed
\(Y_0>e^{14}+2\), are exponentially small relative to
the reflected kernel at the real point on that line.

On the finite segment of the line, Stirling gives
\[
 K_j(jq)=(-1)^{j+1}A(q)e^{jP(q,\zeta)}(1+O(1/j)),
 \qquad
 A(q)=\frac{2\pi}{\sqrt{q-1}}
           \exp\!\left(\tfrac12\sum_i\log(1+r_iq)\right).
\]
The remainders and fixed derivatives are bounded effectively,
since \(\Re(q-1)\ge\eta/2,\ \Re(1+r_iq)\)
are bounded below and the sum of shared normalized
errors uses \(\sum_i1/a_i\le\sigma\).
Binet's formula provides the uniform complex bounds on
these right-half-plane arguments by the same segment
estimate as in the
[effective cubic field proof](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md#effective-lattice-representation-on-a-fixed-small-rectangle).
The amplitude is bounded and bounded away from zero on
the uniform saddle disks, by the log ratio bounds and
\(\sum_i |r_i|\le5/4\).

At \(\zeta=x\) and on \(q=s(x)+iy\), \(y>0\), the real
vertical phase derivative is at most
\(-\arctan(y/(s-1))\).
All row logarithms have nonnegative imaginary contribution
there. The lower vertical half has the same bound when
its coordinate increases downward.
This gives computable strict loss outside any small fixed
saddle segment through the full finite line.
The loss persists for small \(|\beta|\) on the vertical line
through \(q_*(x+i\beta)\), by the uniform complex phase bounds
on the finite segment, (1), and the saddle displacement.
The central exponent differs from the real-phase exponent
by a computably small amount, so the absolute tail estimate
also gives strict loss relative to the saddle exponent
upon reducing \(|\beta|\) if needed.

Local integration on the shifted vertical line gives the
right contribution, with its residue orientation, as
\[
         R_j(\xi)=\pm\sqrt{j}\,s(\xi)e^{jD(\zeta)}
                             (1+O(1/j)),                   \tag{3}
\]
where \(s(\xi)\) is nonvanishing and analytic, with
\(|s|,|1/s|,|s'|\) uniformly bounded.
It is the leading amplitude \(A(q_*)\) divided by
\(\sqrt{2\pi Y'(q_*)}\), with the positive real choice
on real \(\zeta\).
For the relative \(O(1/j)\) error integrate on a
symmetric vertical segment about the saddle.
Subtract the odd linear amplitude term and the odd cubic
phase correction against the even complex Gaussian; their
integrals vanish. Taylor remainder and the quadratic loss
bound the residual by the integral of
\(C(v^2+jv^4+j^2v^6)e^{-c jv^2}\)
times the leading exponential, plus the \(O(1/j)\) Stirling
error. Gaussian moments give the relative estimate.
All constants used to shrink \(\epsilon_*\) are computable,
and the sign in (3) does not enter the modulus comparison.

## Consecutive modulus separation

At a zero, (2) and (3) give
\[
 j\Re D(\xi_j-i\pi)=-\tfrac12\log j+C(\xi_j)+O(1/j),
       \qquad C(\xi)=\log|(1-e^\xi)^{-1}/s(\xi)| .
\]
The function \(C\) is uniformly bounded and Lipschitz along
horizontal field segments throughout the rectangle.
The exact drift is
\(\xi_N-\xi_n=\Delta\log(N/n)\in\mathbb R\).
The modulus of the change in \(\Re D\) along that segment
is at least \(|\xi_N-\xi_n|\), because \(\Re D'=\Re q_*>1\).
Subtracting the two zero identities gives
\[
         (N-C_1)|\xi_N-\xi_n|
                     \le C_2 (1+\log n)/n
\]
with computable constants, using the lower-order identity
to bound \(|\Re D(\xi_n-i\pi)|\).
If \(\Delta\ne0\), then \(|\xi_N-\xi_n|\ge1/N\),
contradicting this bound above a computable terminal.
If \(\Delta=0\), subtraction gives \(\Re D=O(1/n)\),
which is incompatible with
\(n\Re D=-\tfrac12\log n+O(1)\).
This proves the theorem.
