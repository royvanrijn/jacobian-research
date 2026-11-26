# Simple stationary contacts unfolded by the first complex root jet

Consider the [fixed-terminal Gaussian reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data)
at adjacent orders \(n,N=n+1\), \(n\ge2\). Fix a positive interior
physical contact at \(T_*>0,\ |\lambda_*|<T_*\), simple in both
orders. Assume the ordinary retained core is infinite and has
summable curvature. Use \(F_j(T,z)\) for the polynomials in their
own endpoint scales. At logarithmic roots \(\xi_j\) with a common
angle \(\vartheta_*\), put
\[
 G_j(x)=e^{-jx/2}F_j(T_*,e^x),\qquad
       \frac{G_j(\xi_j+x)}{G_j'(\xi_j)}=\sum_{h\ge1}c_{j,h}x^h.
\]
In particular \(c_{j,1}=1\). Set
\[
 k_j=\min\{h:\Im c_{j,h}\ne0\},\qquad k=\min(k_n,k_N).
\]

**Theorem.** If \(k\ge4\), the contact is a limit of simple
transverse interior Gaussian contacts of finite ordinary cores.
The integer total degree difference and the orders are fixed.
The strengths, biases and relative physical root scales converge,
and every strict total curvature cap persists. The same method
works for \(k=3\) if
\[
             \frac{\Im c_{n,3}}{n^2}
                              -\frac{\Im c_{N,3}}{N^2}\ne0 .
\]
In particular the \(k\ge4\) case covers simultaneous vanishing of
the two individual angular velocities and of the two individual
cubic angular responses. It permits bounded or unbounded gaps,
and no fixed gap needs to repeat.

## The first complex jet and how it moves

The indices \(k_j\) are finite and at most \(j\).
Indeed, \(G_j(\xi_j+x)/G'_j(\xi_j)\) is a combination of \(j+1\)
distinct real exponential frequencies \(h-j/2\) with nonzero
coefficients. If all its imaginary Taylor coefficients up to
degree \(j\) were zero, the Vandermonde matrix on those frequencies
would force every exponential coefficient to be real. Their
consecutive ratios are a positive number times \(e^{i\vartheta_*}\),
which is nonreal by the nonreal contact assumption.

For a scale \(s\downarrow0\), consider a real Fourier multiplier
on the original frequencies with the following Taylor control.
To a fixed order greater than \(k\) its symbol is a polynomial
in \(sD\) with bounded real coefficients, converging coefficientwise
on the relevant jets to \(Q(D)\), where \(D\) denotes differentiation
after dilation by \(s\). The remainder on the original finite
frequency sets is \(o(s^k)\). Normalize away positive scalar and
real dilation factors so \(Q(0)=1\) and \(Q'(0)=0\), and take
the same vanishing linear term in the approximating symbols.
If \(A_j=\Im c_{j,k}\), the resulting simple root angle is
\[
       \vartheta_*-s^k A_j\,k!\,[D^k]Q(D)+o(s^k).                \tag{1}
\]
For the proof, dilate \(x=sX\) after centering at the reference
log root, and divide by \(sG_j'(\xi_j)\). The resulting expansion
has leading real polynomial \(X\). Its real Taylor part through
degree \(k\) has a simple real zero \(X_s=O(s)\), with derivative
tending to one. Since all \(c_{j,h}\) with \(h<k\) are real,
the first imaginary term is \(i s^{k-1} A_j Q(D)X^k\).
The symbol remainder after division is \(o(s^{k-1})\).
Simple-root perturbation at \(X_s\) therefore gives imaginary
displacement \(-s^{k-1} A_j Q(D)X^k|_{X=0}+o(s^{k-1})\).
Multiplying by \(s\) yields (1). This argument also applies
when \(k_j>k\), in which case the leading displacement is zero.

## An ordinary family of scaled multipliers

Choose a sequence of retained rows with midpoint \(M\) and signed
gap \(d\), oriented by a common column reversal if needed so
\(d>0\). Set \(\eta=d/M\). Their curvatures \(d^2/M\to0\) by
summability. The following construction uses the
[exact midpoint and cloud operators](MULTIVARIATE_PAIRED_MULTIPLE_MIDPOINT_CLOUDS.md).
Take a fixed positive integer \(L\), raise both endpoints of
the selected row by \(L\), and use effective time
\[
                      T=T_*+\eta^2(L+\tau\sqrt L).
\]
Up to affine terms in the log frequency, the row shift contributes
\(jL\log(1+\eta D/j)\) in the unscaled logarithmic coordinate,
with holomorphic remainder \(O_{j,L}(M^{-1})\) in a fixed disk
in the variable \(\eta D/j\). This follows directly by the finite
Gamma quotient with \(jL\) factors. Cauchy estimates show
that every fixed Taylor coefficient of the remainder tends to
zero in that variable.

Add increasing or decreasing unit-gap rows at a shared weak
midpoint. Their compensated contribution, after dilation by
\(\eta\), can prescribe any fixed cubic \(B_L D^3/j^2\).
Indeed take a cloud with limiting curvature \(w_C>0\),
midpoint \(H_C\sim w_C/(6|B_L|\eta^3)\), and
\(\lfloor w_C H_C\rfloor\) rows of orientation
\(-\operatorname{sgn}B_L\). The case \(B_L=0\) needs no cloud.
Choose \(B_L=\beta L^{3/4}-L/3\), where \(\beta\) is real.
Also add a reflected unit-gap cloud with limiting total curvature
\(w_R>0\). For \(0<q<3\), take its midpoint
\(H_R\sim1/(a\eta^2)\) where \(a=\sqrt{q L/w_R}\), and
\(\lfloor w_R H_R/2\rfloor\) reflected pairs. At scale \(\eta\)
its nonvanishing log terms are
\[
                      -\frac{w_Ra}{4j^2}D^2
                             +\frac{w_Ra^2}{12j^3}D^4 .
\]
The midpoints in both clouds can be rounded to half-integers.
In each profile reduce the true Gaussian strength by the exact
midpoint curvatures of the added clouds.

Now work at the dilation \(s=L^{1/4}\eta\). Cancel affine
terms by positive scalars and a real dilation of the variable
in each order. The quadratic contribution of the midpoint shift
cancels the part \(L\eta^2\) of the effective time change.
The cubic cloud cancels the large cubic of that shift and leaves
the chosen \(\beta\). At each fixed degree in the symbol
the terms of degree \(h\ge5\) from the shift are
\(O_{j,h}(L^{1-h/4})\) after dilation. Consequently the iterated
limits \(\eta\to0\), then \(L\to\infty\), give
\[
  Q_j(D)=\exp\left(a_j(\tau)D^2+\frac{\beta}{j^2}D^3
                                         -\frac{g}{j^3}D^4\right),
 \qquad
  a_j(\tau)=\frac{\tau}{2j}-\frac{p}{4j^2},
  \quad p=\sqrt{q w_R},\quad g=\frac{3-q}{12}>0 .                \tag{2}
\]
The formula is for finitely many jets, which suffice in (1).
For clarity, the unit-gap compensated log is analytic in inverse
midpoint at zero, and its inverse-midpoint term of order \(h\)
has frequency degree at most \(h+1\). In the cubic cloud an
inverse-midpoint layer \(h\ge3\) after dilation by \(\eta\) is
\(O(\eta^{2h-4})\), since the row count is \(O(H_C)\).
In the reflected cloud odd frequency powers cancel. The
inverse-midpoint layer \(h\ge4\), of even frequency degree
at most \(2\lfloor(h+1)/2\rfloor\), vanishes after dilation.
The leading quadratic and quartic terms of that cloud are the
displayed ones at layers two and three. Terms of lower frequency
degree tend to zero after the indicated affine removal. Truncating
these analytic expansions to sufficiently high order gives the
required \(o(s^k)\) remainder on the finite frequency grids,
for each fixed \(L,w_R,w_C,\beta,\tau,q\). This justifies applying
(1) before taking the second limit in \(L\).

## Comparing the two first complex jets

In the limiting models, the phase difference \(\arg z_N-\arg z_n\)
divided by \(s^k\)
has limit
\[
            \Phi(\tau,\beta)
                   = k!\bigl(A_n[D^k]Q_n-A_N[D^k]Q_N\bigr).
                                                               \tag{3}
\]
At least one of \(A_n,A_N\) is nonzero.

If \(k=2\ell\ge4\), put \(\beta=0\). By the generating function
of the [probabilists' Hermite polynomials](https://dlmf.nist.gov/18.12.E16),
\[
 [D^{2\ell}]Q_j
   =\frac{(2g/j^3)^{\ell/2}}{\ell!}\operatorname{He}_{\ell}(t_j),
 \qquad
       t_j=\frac{\tau\sqrt j/2-p/(4\sqrt j)}{\sqrt{2g}}.
\]
At \(\tau_0=-p/(2\sqrt{nN})\) the two arguments are equal to
\[
       t=-\frac{p}{4\sqrt{2g}}\left(\frac1{\sqrt n}
                                             +\frac1{\sqrt N}\right).
\]
For each fixed \(w_R>0\), varying \(q\) in \((0,3)\) makes
this common argument range over all negative reals.
If \(A_n n^{-3\ell/2}-A_N N^{-3\ell/2}\ne0\), choose \(t\)
on either side of a negative simple zero of \(\operatorname{He}_\ell\).
Such a zero exists by its real simple zeros and parity for
\(\ell\ge2\). This realizes either sign in (3).
If the weighted difference is zero, the common weight is
nonzero. Choose \(t<0\) where \(\operatorname{He}'_\ell(t)\ne0\).
The derivative in \(\tau\) at the matched arguments is nonzero,
since the two argument derivatives differ by
\((\sqrt N-\sqrt n)/(2\sqrt{2g})\). Varying \(\tau\) slightly
again gives both signs.

For odd \(k=2r+3\ge5\) the expression (3) vanishes at \(\beta=0\),
and its \(\beta\)-derivative there is the same Hermite comparison
with order \(r\ge1\) and weights \(A_j j^{-2-3r/2}\), up to a
common nonzero factor. If the weighted difference is nonzero,
choose equal arguments at a nonzero value of \(\operatorname{He}_r\).
If it is zero, the same derivative-in-\(\tau\) argument makes
the comparison nonzero. Then choosing the sign of small \(\beta\)
realizes both signs. For \(k=3\) the derivative reduces to the
explicit weighted difference in the theorem.

## Transfer to exact ordinary finite cores

If the original angular difference changes sign, ordinary
finite truncation and the
[transverse ordinary regularization](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
already give the desired contacts. Otherwise its zero has
finite even order, by the endpoint-binomial quotient resultant.
Choose parameters in (2) making (3) have the sign opposite to
the original angular difference at nearby strengths. For
each small \(w_R,w_C\), first take a sufficiently large fixed
integer \(L\), then take a sufficiently weak retained row
from the sequence. Formula (1) and the strict limiting sign
give the desired opposite sign at the adjusted effective time.
Choose finite truncation of the retained core, containing the
selected row, with coefficient and center errors small enough
to preserve that sign. The ordinary perturbations tend to the
reference in the compensated coefficients. Thus they also
preserve the original sign at an endpoint of a small fixed
strength interval around \(T_*\). The interval may be shrunk
along the sequence. The regularization lemma with arbitrarily
weak extra ordinary rows makes the resulting sign-changing
contact simple and transverse.

The shift preserves the signed gap of the selected row and
changes its curvature and center by \(o(1)\).
The reflected cloud's center tends to zero. The signed unit
cloud has limiting center correction of magnitude \(w_C/(2nN)\).
At the equal-angle roots the physical bias is recovered from
the exact center and the log ratio of their moduli. In the
limit with fixed cloud curvatures its deviation from
\(\lambda_*\) is at most \(w_C\), while the actual strength
tends to \(T_*-w_R-w_C\). Take \(w_R,w_C\) small enough for the
strict interior condition. Sending them to zero diagonally
gives the specified limiting strength and bias. The total
signed gap stays fixed through the Gaussian residual data,
and the compensated variance and sum of retained curvature
preserve the strict cap.

For a simple even contact with \(k=2\), the two individual
angular velocities are equal and nonzero, since they are
\(-\Im c_{j,2}/j\). The
[earlier add-delete and signed-cubic tests](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md#signed-unfolding)
already reduce that case on an infinite retained core.
Thus possible remaining simple even contacts of an infinite
core have first complex order three with zero weighted cubic
difference, subject to the separate stationary response tests.
The conclusions concern Gaussian reference contacts.
