# Effective exclusion for escaping shared quartic fields

Use the ordinary shared profiles, curvature \(W,\rho\), and centered
real quartic coefficients \(q_*,w_*,a_2,g\) of the
[moving quartic comparison](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
In particular \(g<0\). Set
\[
 \delta=\max\rho-1,\qquad \beta_j=w_*-w_j,\qquad
 \ell_0=(-4g)^{-1/3},\quad
 \epsilon_j=|\beta_j|^{1/3},\quad M_j=j\ell_0\epsilon_j^4 .
\]
A nonzero field is understood when using \(\epsilon_j\) in a denominator.

**Theorem.** Given rational \(0<W_-\le W_+\), one can compute
positive rational bounds \(d_0,b_0,\tau_0,M_0\) and an order \(n_0\)
as follows. Suppose \(W_-\le W\le W_+\), \(|\delta|\le d_0\),
\(n\ge n_0\), \(N=n+1\), and both fields have \(|\beta_j|\le b_0\).
A consecutive common root is excluded if at either order
\[
               \frac{|a_2|\ell_0}{\epsilon_j^2}\le\frac{\tau_0}{2},
                 \qquad M_j\ge8 M_0 .
\tag{1}
\]
The conclusion is uniform in the exponent profiles and the integer
degree difference \(\Delta\). The constants are unevaluated.

## Effective finite contour near the center

Use the derivatives and center construction of the
[effective bounded quartic proof](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md#effective-center-and-localization-on-the-real-lattice).
For \(h=q-q_*\), the phase after exact centering is
\[
             f_\beta(h)=\beta h+a_2h^2+gh^4+h^5E(h).
\tag{2}
\]
The bounds for \(-g\), \(E\) on a fixed disk, and \(|a_2|/|\delta|\)
are uniform and computable. Take a fixed rational retained radius
\(r>0\) well inside that disk. By decreasing \(r\), the \(h^5E\) term
and its real derivative are small relative to the quartic term
and its derivative. Decrease the caps for \(|\delta|\) and
\(|\beta|\) so that at \(h=\pm r/2\) the real phase is strictly
negative and its real derivative points inward, with computable
uniform margins. Outside this interval, \(\rho(q)-1<0\) by the
strict real concavity bound on \(\rho\) and \(|q_*-q_m|=O(|\delta|)\).
Thus the derivative signs extend the phase loss through both
endpoints. The ordinary entropy bound on the exact binomial
coefficient and the shared Gamma ratio estimate then bound all
discarded summands by a polynomial in the order times a fixed
exponential loss.

On the retained interval and a fixed complex neighborhood, use the
amplitude \(A_j\) of the
[finite integral reduction](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md#uniform-centers-and-finite-integral-reduction).
It satisfies
\[
     A_j=A_0(1+O(j^{-1})),\qquad
     A_{j+1}/A_j=1+O(j^{-2}).
\tag{3}
\]
The magnitudes of \(A_0\) and its reciprocal have uniform bounds
there. In the logarithm of its row part each term is bounded by a
constant times \(|d_i|/m_i\), which sums to at most \(W_+\).
For the remainder, take a fixed relative complex-order disk about
each positive order. The normalized row arguments remain in a
strict positive-real sector and have moduli bounded below by
a fixed multiple of \(j m_i\). The explicit
[Stirling remainder](https://dlmf.nist.gov/5.11.ii) and
\(\sum_i1/m_i\le W_+\) give a computable \(O(j^{-1})\) bound on
that disk. Cauchy's estimate in order gives \(O(j^{-2})\) for its
derivative. The fixed spatial disks give computable bounds for the
needed derivatives of \(A_j\) and \(f_\beta\) as well.

The cotangent rectangle uses endpoints at neighboring half-grid
points. Choose a rational height \(\nu>0\) so that the vertical
walls stay within the strict endpoint phase loss and the increase
\(C\nu^2+|\Im\beta|\nu\) on the horizontal sides is less than
\(\pi\nu/2\). This increase bound follows by Taylor expansion at
the real point, because the nonlinear coefficients in (2) are
real. Decrease the original field and curvature caps so that
the maximum of the real retained phase is less than \(\pi\nu/2\).
The vertical cotangent is bounded by one, and the cotangent's
error from its constant value on the horizontal sides is bounded
by \(2e^{-2\pi j\nu}/(1-e^{-2\pi j\nu})\). The constant-kernel
integrals deform back to the real interval. Moving the half-grid
ends back to \(\pm r\) costs exponentially little by the same
loss at those endpoints. After dividing the moment by its nonzero
scalar center and common order factor, the exact value is therefore
\[
        I_j(\beta)+{\cal R}_j(\beta),\qquad
        I_j(\beta)=\int_{-r}^r A_j(q_*+h)e^{j f_\beta(h)}\,dh,
        \qquad |{\cal R}_j|\le C e^{-\gamma j},
\tag{4}
\]
with computable \(C,\gamma>0\). The factors polynomial in \(j\)
are absorbed by reducing the positive exponential constant.
The fixed endpoint loss also holds on any of the small vertical
end connectors below, after shrinking the field caps.

## One-field confinement and two horizontal packets

At one order set \(\epsilon=|\beta|^{1/3}\) and \(h=\ell_0\epsilon z\).
Division of (2) by \(\ell_0\epsilon^4\) gives
\[
       G_\zeta(z)+\tau z^2+\epsilon\ell_0^4 z^5 E(\ell_0\epsilon z),
       \qquad G_\zeta(z)=\zeta z-z^4/4,\qquad
       \zeta=\beta/\epsilon^3,\quad \tau=a_2\ell_0/\epsilon^2 .
\tag{5}
\]
In this substitution \(|\zeta|=1\). On each fixed scaled disk the
perturbation and a fixed number of its derivatives are bounded by
computable multiples of \(|\tau|+\epsilon\). In the horizontal tails
the same estimate uses \(\epsilon\ell_0 |z|\le2r\) to absorb
the remainder into a small multiple of \(|z|^4\) by the choice of
\(r\). For \(|\Im z|\le2\) and \(|\Re z|\) beyond a computable
constant, this gives a bound \(-c|\Re z|^4\) on the real phase,
uniformly for bounded \(\zeta,\tau\).

The model at \(\zeta=i\) has saddles
\(z_R=\sqrt3/2+i/2\) and \(z_L=-\sqrt3/2+i/2\).
Their horizontal half-lines have strict quadratic loss and
unique half-line maxima. The real stationary heights are
\(-3/8\), and the split point \(i/2\) has height \(-33/64\).
Also \(\Re(z_R-z_L)=\sqrt3\), and the imaginary part of the
stationary value difference is \(3\sqrt3/4\).
Thus a computably small convex disk of values of \(\zeta\)
about \(i\), and computably small bounds on \(|\tau|,\epsilon\),
preserve the strict margins, using the actual continued saddles.
For example, take disjoint small saddle disks, use Rouché
against the simple linearizations of \(G'_\zeta\) there, and
apply interval subdivision on the bounded remaining portions
of the model half-lines. Taylor's inequality supplies the
quadratic loss locally, and the quartic tail bound supplies
the unbounded parts. Vertical connectors close to the split
point at height \(1/2\) retain the strict loss below both
heights. These inequalities are strict for the unperturbed
model, so the searches produce positive margins. The scaled
phase and its derivatives are uniformly bounded on the fixed
model disks used in the tests.

Choose a small rational cone \(|\Re\beta|\le s|\Im\beta|\)
so that any two fields in this cone of equal positive imaginary
part, scaled by \(|\beta_n|\), and their joining segment
lie in the chosen convex disk about \(i\). This follows from
\(|\Im\beta_n|/|\beta_n|\in[1/\sqrt{1+s^2},1]\).
The cone can also be chosen so that the saddle disks used for
the two orders after their own normalizations refer to the
same physical saddle branches on this disk. Indeed the ratio
of the two lengths \(\epsilon_j\) tends uniformly to one
as \(s\) decreases and the centers are the same for both
orders. Use smaller saddle disks in (5) if necessary. The
negative imaginary cone is treated by conjugation.

Off these cones, use the single horizontal contour through the
simple dominant model saddle. For \(\zeta=e^{i\theta}\),
\(0\le\theta<\pi/2\), it is \(e^{i\theta/3}=u+iv\).
On its horizontal line the real phase, up to a constant, is
\[
            -x^4/4+3v^2x^2/2+x\cos\theta .
\]
The positive half-line wins against the reflected negative
half-line and has exactly one critical point, a strict maximum,
since the cubic derivative has exactly one positive zero and
\(u^2>v^2\). Reflection and conjugation treat the other angles.
The off-cone closed angular set gives a compact family of strict
gaps on bounded parts of the contour away from the model
saddles. Rational interval searches verify these gaps as above.
Choose common perturbation caps for these tests and the
two-saddle cone tests.

Here are the uniform packet estimates following from the gaps.
A shift of the full real interval, or a split at scaled
height \(1/2\) followed by separate horizontal shifts of its
two halves, consists of horizontal saddle paths and the stated
connectors. Each half can be assigned its own outer physical
connector. This is an exact contour replacement for (4).
On a short horizontal segment through a saddle \(z_s\), write
the scaled phase difference at the saddle as
\(\Psi''(z_s)v^2/2+O(|v|^3)\), where \(v=z-z_s\) is real.
The interval can be chosen so the cubic error leaves strict
quadratic loss in the whole line segment from the quadratic
exponent to the actual exponent. Integrating the difference of
the exponentials bounds its absolute integral by a constant
times \(M\int |v|^3 e^{-cMv^2}dv\), with \(M=j\ell_0\epsilon^4\).
Replacing the amplitude by its value at the saddle costs at
most a constant times \(\int |v|e^{-cMv^2}dv\).
The complex Gaussian coefficient is explicit, nonzero and
uniformly bounded above and away from zero. Bounded complements
have a fixed scaled loss below the relevant saddle, and the
unbounded horizontal pieces beyond a fixed large scaled radius
have quartic loss. These estimates give relative \(O(M^{-1/2})\)
errors and computable constants. They are direct Gaussian
estimates on the chosen contours, without a singular saddle limit.

In particular each retained saddle packet has magnitude between
positive computable multiples of
\(\epsilon M^{-1/2}e^{j\Re f_\beta(h_s)}\). It has an absolute
Gaussian majorant on its local segment of the same size.
Its magnitude is bounded below by
\(c\epsilon M^{-1/2}e^{-C M}\).
First shrink the field cap so the rate \(CM/j\) is less than
one quarter of the physical gap in (4) and the end connectors.
As \(M\ge1\) implies \(\epsilon^{-1}\le(j\ell_{0,+})^{1/4}\),
relative physical remainders are exponentially small, with
computable constants. These estimates imply a single-order
confinement statement. For \(|\tau|\le\tau_0\), sufficiently
small \(|\beta|\), sufficiently large \(M\) and \(j\), a root
must lie in the cone. Outside it there is one nonvanishing
dominant packet. All constants can be chosen simultaneously
for both conjugate cones.

If (1) holds at one order of a common root, that order lies
in the cone and \(|\Im\beta_j|\ge|\beta_j|/2\).
The two fields have equal imaginary parts, because their
difference is exactly \(\Delta\log(N/n)\).
The other field has \(|\beta|\ge|\beta_j|/2\), so it also
has \(|\tau|\le\tau_0\) and \(M\ge M_0\) when (1) holds.
It therefore lies in the cone as well.

## Linked comparison at a common root

Use the lower-order \(\epsilon\) and \(M=n\ell_0\epsilon^4\).
Write \(D(\beta)\) for the right minus left stationary value
of the physical centered phase \(f_\beta\). At both roots,
cancellation of the two actual packets in (4) implies
\[
           |\Re D(\beta_j)|\le C/j .
\tag{6}
\]
Indeed the Gaussian amplitudes and relative errors are bounded
above and away from zero after taking out the common
\(\epsilon_j M_j^{-1/2}\) factor, and the remainder is smaller
than either packet. Along the real segment joining the two
fields, analyticity and stationarity give
\[
        D_\beta=h_R-h_L,\qquad
        \Re D_\beta\ge c\epsilon,\qquad
        |\Im D(\beta_n)|\ge c'\epsilon^4 .
\tag{7}
\]
These are the strict normalized margins at \(\zeta=i\)
just verified. Integration of the first bound across the
segment and (6) give \(d=|\Delta|\epsilon\le C_1\), using
\(\log(N/n)\ge1/N\).

Now keep the exact lower-order contours at the upper order.
Their integrand multiplier is \((A_N/A_n)e^{B(h)}\), where
\[
          B(h)=f_{\beta_n}(h)+N\Delta\log(N/n)h .
\tag{8}
\]
On either half the integral at the upper order, divided by the
lower packet, equals \(e^{B(h_s)}(1+E_s)\), with
\[
          |E_s|\le C_2\left(\frac{\epsilon^4}{M}
                        +\frac{d}{\sqrt M}+n^{-2}\right)
                      +C_2 e^{-\gamma' n}
\tag{9}
\]
for computable \(C_2,\gamma'>0\).
To obtain this estimate, subtract \(B(h_s)\) within the ratio.
On the local segment, lower-order stationarity gives
\[
 |B(h)-B(h_s)|\le C\bigl(\epsilon^4 |v|^2+d|v|\bigr).
\]
Insert \(e^{B(h)-B(h_s)}-1\) in the lower integral. The inequality
\(|e^u-1|\le |u|e^{|u|}\), the Gaussian majorant, and the bound
on \(d\) give the first two terms of (9) after division by the
packet. On a bounded complement the same insertion is at most
a constant times \(\epsilon^4+d\) in addition to the fixed
scaled loss. On the horizontal tail the insertion is bounded
by a constant times
\[
     [\epsilon^4(1+|x|^4)+d(1+|x|)]
                 \exp(C[\epsilon^4(1+|x|^4)+d(1+|x|)]).
\]
A sufficiently large fixed cutoff in \(x\), the bound on \(d\),
and the quartic decay absorb the exponential uniformly for
\(M\ge M_0\). The tail and the bounded complement are then at
most \(C(\epsilon^4+d)\) times an exponentially small function
of \(M\) relative to the packet, which fits (9).
The amplitude ratio (3) costs \(O(n^{-2})\).
At the fixed outer physical connectors the multiplier has at
most \(e^{C(1+|\Delta|)}\) growth. Since
\(|\Delta|\le C_1/\epsilon\le C_1(n\ell_{0,+})^{1/4}\),
the fixed physical loss absorbs it. This proves the last term
of (9) and gives the same bound for the compared exact
remainder in (4). In particular these estimates preserve a
factor of the actual increment, rather than subtracting two
independent \(O(M^{-1/2})\) asymptotic estimates.

The two stationary multipliers are bounded above and away from
zero, since \(|B(h_s)|\le C(\epsilon^4+d)\). If \(\Delta\ne0\),
the real difference of their exponents is
\[
   \Re D(\beta_n)+N\Delta\log(N/n)\Re(h_R-h_L).
\]
By (6), (7), and \(d\ge\epsilon\), its absolute value is at least
a positive multiple of \(d\) for \(M\) large. In fact
\(n^{-1}=\ell_0\epsilon^4/M\) bounds the first term by \(C d/M\).
Thus the two multipliers differ by at least \(c_3 d\).
When \(\Delta=0\), the exponent difference is \(D(\beta_n)\).
It has modulus bounded below by \(c'\epsilon^4\) from (7) and
above by \(C\epsilon^4\) from (5). Shrink the field cap so it
has modulus at most \(1/2\). The two multipliers then differ
by at least \(c_4\epsilon^4\).

At two roots, the two lower packets sum to the negative of the
exponentially small lower remainder. The two upper packets sum
to the negative of the upper remainder. Substitute the exact
ratios (9) in the latter equality and use the former.
It bounds the difference of the stationary multipliers by
a constant times the error in (9) plus the two relative
physical remainders. Increase \(M_0\) so this bound is smaller
than the stated separation in both cases. The terms \(n^{-2}\)
and the physical remainders are also dominated effectively,
since \(\epsilon^4\ge M_0/(n\ell_{0,+})\) and in the nonzero
degree case \(d\ge\epsilon\ge\epsilon^4\).
This is a contradiction.

All parameter choices have a finite order. Fix the physical
radius and its strict losses, certify the model contours on
bounded intervals and the quartic tail, then choose small
temperature, angular and field caps. Next choose the saddle
comparison threshold \(M_0\) from the computable bounds in
(6)–(9). Finally increase \(n_0\) to absorb the physical
remainders and any polynomial prefactors. This proves the
one-order criterion (1) with no numerical cutoff assumed.
