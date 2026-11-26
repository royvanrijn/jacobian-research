# Mixed multiplicities from the first complex jets

Consider a positive interior physical contact of the
[fixed-terminal reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at adjacent orders \(n,N=n+1\). Assume that the contacting root at
one order \(j_S\) is simple and the root at the other order \(j_M\)
has multiplicity \(m\ge2\). The ordinary retained core is infinite
and has summable curvature. In the individual endpoint scales write
\[
 G_j(x)=e^{-jx/2}F_j(T_*,e^x),\qquad
 G_j(\xi_j+x)=g_jx^{m_j}
           \left(1+\sum_{h\ge1}c_{j,h}x^h\right),
 \quad m_{j_S}=1,\quad m_{j_M}=m .
\]
Here \(\xi_j\) are logarithms of the contacting roots, with common
imaginary part \(\vartheta_*\), and \(g_j\ne0\). Put
\[
 r_S=\min\{h\ge1:\Im c_{j_S,h}\ne0\},\qquad
 r_M=\min\{h\ge1:\Im c_{j_M,h}\ne0\}.
\]
These indices are finite by independence of the distinct real
exponential frequencies, as in the
[first complex jet argument](MULTIVARIATE_SIMPLE_CONTACT_COMPLEX_JET_UNFOLDING.md#the-first-complex-jet-and-how-it-moves).
Suppose \(r_S,r_M\ge2\). The case where one of the indices is one
is covered by the
[angular rank theorem](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md).

**Theorem.** The contact is a limit of simple transverse interior
contacts of finite ordinary cores in the following cases.

1. \(r_S<r_M\).
2. \(r_S=r_M\ge3\).
3. \(r_M<r_S\) and \(r_M\) is odd.
4. \(m\) is odd, \(r_M\) is even, and \(r_M\le r_S\).
5. \(m=r_M=2\) and \(2A_{j_M}/j_M^2\ne A_{j_S}/j_S^2\), where
   \(A_j=\Im c_{j,2}\).

The orders and the total integer degree difference are preserved.
The strengths, biases and relative physical root scales converge
to the reference values, and every strict total curvature cap
persists. Case 1 also applies to a finite retained core.

Thus among these mixed contacts the residual case on infinite
cores has even multiple-root multiplicity and either a strictly
earlier even multiple-root jet or a tie at order two.

## Local angles from a real operator

At a scale \(s\downarrow0\), take a real coefficient multiplier
on the original frequencies whose normalized symbol, on the
dilated coordinate \(X=(x-\xi_j)/s\), tends coefficientwise on
finite jets to a real operator \(Q_j(D)\), \(D=\partial_X\).
The symbols in this proof are normalized to constant one and
linear coefficient zero. The discarded high-order terms have
arbitrarily high order in \(s\) on the original finite frequency
sets, by Taylor expansion of the log multipliers. Put
\[
                           P_{j,k}(X)=Q_j(D)X^k .
\]
Let \(r=\min(r_S,r_M)\) and \(A_j=\Im c_{j,r}\), allowing zero.
At any real simple root \(h\) of \(P_{j,m_j}\), the corresponding
logarithmic root angle has the expansion
\[
 \vartheta_* - s^{r+1} A_j
                       \frac{P_{j,m_j+r}(h)}{P_{j,m_j}'(h)}
                        +o(s^{r+1}).                           \tag{1}
\]
To see this, divide the dilated branch polynomial by \(g_js^{m_j}\).
Use the real Taylor terms through relative order \(r\) with the
actual normalized real symbol. This gives a real polynomial
with a real simple root tending to \(h\), and its derivative there
tends to \(P'_{j,m_j}(h)\). The difference has leading imaginary
term \(is^r A_j P_{j,m_j+r}(X)\), since each earlier normalized
Taylor coefficient is real. The coefficient errors after division
are \(o(s^r)\) by truncating the original symbols sufficiently far.
Simple-root perturbation and multiplication of the displacement
by \(s\) prove (1). A nonreal leading root of \(P_{j,m_j}\) has
imaginary displacement of order \(s\), with sign given by that
root's imaginary part. Real affine shifts of the logarithmic
variable do not affect any of these comparisons.

In particular, if two fixed effective clocks have different counts
of multiple-branch angles strictly below the simple-branch angle,
and each tested angle difference is strict, some fixed increasing
angular rank reverses. Finite truncation and the
[ordinary regularization lemma](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
then give a simple transverse angular contact. In transferring
a strict comparison at scale \(s\), the truncation and regularizing
row perturbations can be chosen smaller than all the tested
angle differences. Root labels for the multiple branch are
needed only by increasing angular rank on its full isolating disk.

## Ordinary operator models

The [midpoint and cloud construction](MULTIVARIATE_PAIRED_MULTIPLE_MIDPOINT_CLOUDS.md)
and its [first-jet variant](MULTIVARIATE_SIMPLE_CONTACT_COMPLEX_JET_UNFOLDING.md#an-ordinary-family-of-scaled-multipliers)
give the following real models:
\[
 \begin{split}
 Q_j(D)&=\exp\bigl(a_j(t)D^2+\beta D^3/j^2-gD^4/j^3\bigr),\\
 a_j(t)&=t/(2j)-\sqrt{q w_R}/(4j^2),\qquad
                     g=(3-q)/12,\qquad 0<q\le3 .               \tag{2}
 \end{split}
\]
Here \(\beta\) is arbitrary and \(w_R>0\) can be as small as desired.
For the endpoint \(q=3\) the reflected cloud cancels the quartic
term exactly in the model. These are iterated limits, first in
the weak scale \(\eta=d/M\) for a fixed midpoint shift \(L\), then
in \(L\to\infty\) with \(s=L^{1/4}\eta\). Thus a strict sign test in
the model is implemented by first fixing \(L\) large enough and
then taking \(\eta\) sufficiently small. The two zeroes of \(a_j(t)\)
are distinct because \(q w_R>0\). Without a midpoint shift,
a single cubic cloud and the common effective clock supply the
model \(\exp(tD^2/(2j)+\beta D^3/j^2)\) for finite or countable cores.

We also need a quintic model. Raise the selected midpoint by \(L\),
use \(s=L^{1/5}\eta\), and set effective strength to
\(T_*+\eta^2(L+tL^{2/5})\). Cancel the shifted cubic term with a
signed unit cloud of coefficient \(-L/(3j^2)\) at scale \(\eta\).
For the reflected cloud take limiting curvature
\(\widetilde w_R=w_R L^{-1/5}\) and midpoint parameter
\[
                   A_R=\sqrt{(3L+12c L^{4/5})/\widetilde w_R}
\]
in the notation of that construction. Any fixed real \(c\) is
admissible for large \(L\). At scale \(\eta\), this cloud adds
quadratic coefficient \(-\widetilde w_R A_R/(4j^2)\) and quartic
coefficient \((3L+12c L^{4/5})/(12j^3)\).
The midpoint log coefficient of degree \(h\ge2\) is
\((-1)^{h+1}L/(h j^{h-1})\) at that scale.
After dilation by \(L^{1/5}\) the remaining coefficients of degree
\(h>5\) tend to zero, and the model becomes
\[
 Q_j(D)=\exp\left((t/(2j)-\sqrt{3w_R}/(4j^2))D^2
                          +cD^4/j^3+D^5/(5j^4)\right).          \tag{3}
\]
This follows using the same analytic estimates on the original
finite frequency sets, with the cloud parameters held fixed
during the first limit for each \(L\). In particular the cubic
cancellation leaves no cubic term in (3). Compensating the two
cloud variances and recovering the bias works as in the paired
midpoint theorem. For each choice of the small cloud curvatures,
the chosen clocks approach the original effective strength.
The center changes by at most their total limiting curvature
in the bias parameter. The subsequent diagonal choice with both
curvatures tending to zero retains the physical interior and
strict curvature margins.

## Earlier or tied simple-branch jets

If \(r_S<r_M\), take negative common heat and a sufficiently small
nonzero cubic coefficient. In the limit where the cubic is zero,
all the multiple model roots are simple and real. They stay real
for sufficiently small cubic. The simple model polynomial is \(X\),
and \([D^{r_S+1}]\exp(aD^2+bD^3)\) is nonzero for a suitable small
\(b\), with fixed \(a<0\). For even power this holds already at
\(b=0\), and for odd power the linear term in \(b\) is nonzero.
Equation (1) gives a nonzero simple-branch shift, while the
multiple-branch shifts at that order all vanish. At a sufficiently
large positive clock the multiple leading model has nonreal roots
on both sides of the real axis, by dilation to pure positive heat.
The simple leading root remains real, so some angular rank reverses.
This proof uses only the single-cloud model.

Suppose \(r_S=r_M=r\ge3\). Use (2) with \(q=3\), and fix \(w_R>0\).
Take the clock where \(a_{j_M}(t)=-h^2\) with \(h>0\) small.
Then \(a_{j_S}(t)\) tends to a nonzero number as \(h\to0\).
For odd \(r\), take \(\beta=0\). For even \(r\), take
\(\beta=b h^3\) for a fixed sufficiently small \(b\ne0\).
Dilating \(X=hY\) in the multiple model leaves fixed negative
heat and a sufficiently small cubic. All its roots are therefore
real and simple. Its ratio in (1) is \(O(h^{r+1})\), by dilation.
At the simple branch, the ratio tends to a nonzero constant for
odd \(r\). For even \(r\) its size is asymptotic to a positive
constant times \(h^3\), from the linear term of the cubic.
Since the even tied case has \(r\ge4\), the multiple shifts are
smaller than the simple shift in either case. Choose \(h\) fixed
sufficiently small and proceed to positive clocks as in Case 1.

## An earlier odd multiple-branch jet

Suppose that \(r=r_M\) is odd and \(r_S>r\). Use the even model
in (2), with \(\beta=0\) and \(0<q<3\). In the multiple
order write \(P_k=e^{aD^2+cD^4}X^k\) with \(c<0\) fixed and
\(a\) tunable. There are only finitely many \(a\) for which \(P_m\)
has a repeated root or shares a root with \(P_{m+r}\).
For the repeated roots this follows by comparison at large
negative \(a\) with pure Hermite heat. For the common roots,
rescale that comparison to \(e^{-D^2/2+\delta D^4}\) with
\(\delta\to0\), \(\delta\ne0\). At a common nonzero root \(x\)
of the original Hermite polynomials of degrees \(m,m+r\),
the normalized fourth derivatives are
\(x^3+(3-2k)x\) for \(k=m,m+r\), by the Hermite differential
equation. Thus the infinitesimal motions under the quartic
are different since \(rx\ne0\). There is no common root zero,
by opposite parity and the nonzero even Hermite constant terms.
This verifies that the resultant in \(a\) is nonzero.

At suitable large positive \(a\), pure positive heat after
dilation gives \(m/2\) negative angular displacements and as
many positive ones if \(m\) is even. If \(m=2u+1\) is odd,
there are \(u\) on each side in the leading model and a real
center zero. The ratio in (1) at that zero is positive,
since \(P'_m(0)\) and \(P_{m+r}(0)\) are positive for sufficiently
large positive \(a\). The simple branch is negligible at this
order. Multiply angle differences by a common sign if needed
so that the odd-degree center at this positive clock counts as
positive. This convention means the sign of a real-root correction
is the sign of \(P_{m+r}/P'_m\).

There is a second clock with a different number of negative
angles. If \(m\) is even, choose \(a\) such that \(P_m(0)<0\).
The constant term is \(m![z^{m/2}]\exp(az+cz^2)\), a
positive multiple of a Hermite polynomial after a real
rescaling of \(a\), so that such an interval exists. Avoid the
finitely many exceptions above. There are then an odd number
\(v\) of positive real roots of \(P_m\), from its sign at
zero and positive leading coefficient. The nonreal roots
contribute \((m-2v)/2\) negative angles. Each real pair
\(\pm x\) has equal corrections in (1), by opposite parity of
\(P_m,P_{m+r}\). It follows that the total negative count
has parity \(m/2-v\), different from the positive-clock count.

If \(m=2u+1\), choose \(a\) with \(P_{m+r}(0)<0\) by the
same constant-term argument, avoiding exceptions. Write \(v\)
for the number of positive real roots of \(P_m\). The sign of
\(P_m'(0)\) is \((-1)^v\), since each pair of real roots gives
a negative contribution to the product and each conjugate
nonreal pair a positive one. Hence the zero root contributes
one negative angle precisely when \(v\) is even. The other
real pairs contribute an even number of negative angles,
and the nonreal pairs contribute \(u-v\). Thus the negative
count has parity opposite to \(u\), the positive-clock count.
Both parity arguments give a strict rank reversal.

## The central angle at odd multiplicity

Let \(m\) be odd, \(r=r_M\) even, and \(r_M\le r_S\).
For \(r\equiv0\bmod4\), use (2) with \(g>0\) arbitrarily
small and a small cubic. Test the multiple order at
\(a=-1/2\) and \(a=1/2\). In the zero-cubic limit, the
center root is zero at both clocks by parity. At the negative
clock the noncentral real-root corrections form opposite
nonzero pairs after choosing small \(g>0\). Indeed the
pure Hermite common-root issue is resolved to first order
by the fourth-derivative computation above at each nonzero
candidate. At the positive clock the noncentral leading
roots are nonreal, evenly distributed on the two sides.
These properties persist for sufficiently small cubic, with
the central angular rank corresponding to the root near zero.

The derivative of the central ratio in (1) with respect to
the coefficient \(b\) of \(D^3\), at pure heat and \(b=0\),
is
\[
 \frac{P_{m+r}'(0)}{P_m'(0)}
           \left(\frac{P_{m+r}'''(0)}{P_{m+r}'(0)}
                        -\frac{P_m'''(0)}{P_m'(0)}\right)
     =\frac{P_{m+r}'(0)}{P_m'(0)}\frac{r}{2a}.
\]
Here all \(P_k=e^{aD^2}X^k\) in the displayed formula.
It follows by differentiating the root and using the Hermite
differential equation. The ratio of the first derivatives
is a positive constant times \(a^{r/2}\). Since \(r/2-1\)
is odd, the central ratios have opposite signs for the two
clocks at small nonzero \(b\). This derivative sign persists
when the quartic is small. First choose the quartic to resolve
the noncentral zeros and then take the cubic sufficiently
small for the central test and the opposite noncentral pairs.
The simple branch contributes only at a higher order than this
central comparison when \(r_S>r\). A tie here has \(r\ge4\) and is
already included in Case 2.

For \(r\equiv2\bmod4\), use (3) with a fixed nonzero quartic
parameter \(c\), at the two clocks \(a=\pm A\) for large fixed
\(A>0\). After rescaling by \(\sqrt A\), the quartic has
size \(O(A^{-2})\) and the positive quintic has size
\(O(A^{-5/2})\). With zero quintic, the same fourth-derivative
argument gives opposite nonzero corrections for the noncentral
real roots at the negative clock, with magnitude at least
a constant times \(A^{-2}\) in the rescaled model if they
were shared at pure heat. The smaller quintic preserves
those signs. At the positive clock the noncentral roots
remain nonreal.

The derivative of the central ratio in the quintic parameter
at pure heat is, by the same root differentiation,
\[
 \frac{P_{m+r}'(0)}{P_m'(0)}
           \left(\frac{P_{m+r}^{(5)}(0)}{P_{m+r}'(0)}
                        -\frac{P_m^{(5)}(0)}{P_m'(0)}\right)
     = \frac{P_{m+r}'(0)}{P_m'(0)}
                       \frac{r(2m+r-4)}{4a^2}.
\]
The derivative identity used here is
\(P_\ell^{(5)}(0)/P_\ell'(0)=(\ell-1)(\ell-3)/(4a^2)\) for odd
\(\ell\), including \(\ell=3\). Now \(r/2-2\) is odd,
so the two signs are opposite. The central ratio vanishes
identically at zero quintic for all even operators, so the
small quartic changes the linear-in-quintic sign by \(o(1)\)
after the \(\sqrt A\) rescaling. It follows that the central
rank reverses for large fixed \(A\). If \(r_S>r\) this proves
the comparison directly. A tie with \(r\ge6\) was already
covered by Case 2. At \(r_S=r=2\), the simple shift in (1)
vanishes because the symbol in (3) has no cubic coefficient.
It therefore cannot mask the central reversal. This completes
Case 4.


## The double root at relative order two

In Case 5 use (2) at \(q=3\), with a fixed negative heat
coefficient \(a=a_{j_M}(t)\) and arbitrary cubic coefficient
\(b=\beta/j_M^2\). The roots of \(P_{j_M,2}\) are
\(h=\pm\sqrt{-2a}\), unaffected by the cubic. At either root,
\[
     \frac{P_{j_M,4}(h)}{P'_{j_M,2}(h)}=2ah+12b.
\]
The simple ratio is \(6\beta/j_S^2\). Thus the difference of
the leading angular shifts is the sum of a common term
proportional to
\(\beta(2A_{j_M}/j_M^2-A_{j_S}/j_S^2)\) and an opposite
bounded term at the two ranks. Choose a fixed \(|\beta|\)
large enough for both differences to have the same strict sign.
At a clock with \(a_{j_M}(t)>0\) the two leading roots are
purely imaginary of opposite signs. Some rank therefore reverses.
