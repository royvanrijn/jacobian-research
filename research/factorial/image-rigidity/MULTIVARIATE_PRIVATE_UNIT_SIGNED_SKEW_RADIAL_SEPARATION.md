# Cutoff free radial separation from integer drift and matched skew

Consider an ordinary two-monomial profile with one private unit
row in each orientation and arbitrary shared rows. Write \(\Delta\)
for the sum of the signed gaps of the shared rows. Omit equal
rows and suppose \(\Delta\ne0\). Set
\[
 \kappa=\sum_{\rm shared}\frac{d_i^2}{\ell_i+|d_i|+1/2},
        \qquad b=1-(1+\kappa/4)^{-1},
\]
where \(d_i\) is a signed gap and \(\ell_i\) is the smaller
exponent. One may remove any complete mirrored pairs before
computing the two quantities below. Define
\[
 \begin{split}
 L_+&=\sum_{\operatorname{sgn}d_i=\operatorname{sgn}\Delta}
                                  \log(1+|d_i|/\ell_i),\\
 K_m(v)&=\sum_i d_i^2\,m\psi'(1+m(M_i+d_i(v-1/2))),
          \qquad M_i=\ell_i+|d_i|/2,\\
 H&=\sup_{m\ge2,\ 0\le v\le1}
                              |K_m(v)-K_m(1-v)|/2 .
 \end{split}
\]

**Theorem.** The ordinary physical profile satisfies SFC(2)
in every degree if
\[
               |\Delta|-L_+/2\ge\frac1{16},
                   \qquad H\le b^5/8192.                     \tag{1}
\]
The upper \(H\) bound may be replaced by an explicit matched
tanh moment condition below. The proof has no terminal cutoff.
It uses separation of the absolute values of the physical roots.
The total curvature and the number of shared rows may be large.

To evaluate (1) by rational tests, one can bound \(L_+\) by
\(\sum_{\operatorname{sgn}d_i=\operatorname{sgn}\Delta}|d_i|/\ell_i\).
For the rows kept in \(K_m\), put
\[
 \delta=\max |d_i|/M_i,\qquad S=\left|\sum d_i^3/M_i^2\right|,
                              \qquad D=\sum d_i^4/M_i^3 .
\]
If \(\delta\le1\), then \(H\le S/2+5D/12\) by the
exact odd-curvature estimate for the matched symbol.
In the joint tanh coordinates,
\(\delta=2\max v\), \(S=4|\sum\epsilon v^2|\), \(D=8\sum v^3\),
where \(\epsilon\) is the orientation and
\(v=(n/T'-(n-1)/T)^{-1}\); sums count all unmatched slots.
Complete ordinary rows and \(\Delta\) are reconstructed from the
same signed data.
[Odd curvature estimate](MULTIVARIATE_PRIVATE_UNIT_DIFFUSE_SIGNED_TANH_ANGLES.md#even-and-odd-curvature-of-the-residual)
[Exact matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## A symmetric comparison at each individual order

After canceling the private factorials, make both endpoint
coefficients one by positive scalings and write \(c_k\)
for the normalized degree-\(m\) coefficients. If
\(h_k=\min(k,m-k)\) and \(q=1-b\), the trigamma curvature
bound gives both \(c_k\le q^{h_k}\) and \(c_{m-k}\le q^{h_k}\)
by the quadratic Green chord bound.
Let \(\widehat c_k=(c_k+c_{m-k})/2\).

The coefficients \(c_k\) are the moments \(\mathbb E T^k\)
of a positive variable, under the order-\(m\) endpoint tilt,
with \(\mathbb E T^m=1\). The law of \(1/T\) under the
additional tilt by \(T^m\) has moments \(c_{m-k}\).
Mix the two laws equally. The result has moments
\(\widehat c_k\) and is invariant under order-\(m\)
size bias followed by reciprocal inversion.
Thus its midpoint log law is even. The geometric
resolvent formula puts the roots of
\(\widehat P_m(z)=\sum_{k=0}^m\widehat c_k z^k\)
simply on the circle with the usual angle boxes.
[Geometric resolvent and boxes](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md#the-exact-resolvent-and-angle-boxes)

The same-order derivative proof in the
[near mirror chart](MULTIVARIATE_PRIVATE_UNIT_NEAR_MIRROR_QUADRATIC_ANGLES.md#the-mirrored-reference-slope)
uses only this symmetric positive moment law and
\(\widehat c_k\le q^{h_k}\). It gives
\(\lvert \widehat P_m'(\zeta)\rvert\ge mb/8\) at each
circle root \(\zeta\), and disjoint root disks of radius
\(\epsilon_m=b^2/(128m)\) with boundary modulus at least
\(b^3/2048\).

The odd part of the normalized log coefficient has
modulus at most \(Hh_k/2\), by the Green chord estimate
from its signed curvature. The complete mirrored pairs
have zero odd part. Since for real \(a,b\)
\(\lvert e^a-e^b\rvert/(e^a+e^b)=|\tanh((a-b)/2)|\),
\[
        |c_k-\widehat c_k|\le q^{h_k}H h_k/2 .
\]
On the root disks the evaluation error is at most
\(2H/b^2\), because \((1+\epsilon_m)^m<2\) and
\(\sum_k h_k q^{h_k}\le2/b^2\).
Under (1) this is at most \(b^3/4096\).
Rouché pairs each normalized ordinary root with a
symmetric root at distance less than \(\epsilon_m\).
Every normalized log radius is therefore bounded by
\(2\epsilon_m=b^2/(64m)\).
This estimates each degree by its own symmetrization,
without comparing the symmetrizations across degrees.

## The physical endpoint drift

Write \(G_m\) for the positive argument scale from
endpoint normalization, so a physical root \(z\)
has normalized value \(G_m z\).
The private unit rows cancel in \(G_m\).
A shared row with exponents \(a,b\) contributes
\(\int_a^b\psi(1+mx)\,dx\) to \(\log G_m\).
Its order derivative has leading term \((b-a)/m\)
and an error with opposite sign. In magnitude the
error is at most \(\log(1+|b-a|/\min(a,b))/m^2\),
using
\[
              0\le1-u\psi'(1+u)\le \frac1{1+u}<\frac1u .
\]
This follows from the upper and lower integrals of the
decreasing trigamma summands. Mirrored pairs cancel
exactly, so for \(N=n+1\) the bound in the drift direction is
\[
 \operatorname{sgn}(\Delta)\log(G_N/G_n)
        \ge|\Delta|\log(N/n)-L_+/(nN).
\]
Since \(n\ge2\) and \(\log(N/n)>1/N\), (1) makes this
strictly greater than \(1/(16N)\).
On the other hand overlapping physical radii require
\[
 |\log(G_N/G_n)|
        \le\frac{b^2}{64}(1/n+1/N)\le \frac{5}{128 N},
\]
which contradicts \(5/128<1/16\).
The first and second moments are coprime because the
first has a real root and the second is strictly positive
on the real axis. Rational Bézout identities transfer
complex coprimality to every characteristic-zero field.
