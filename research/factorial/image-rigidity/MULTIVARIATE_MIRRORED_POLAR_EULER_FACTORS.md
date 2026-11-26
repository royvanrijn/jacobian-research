# Exact Euler factors for the mirrored polar multiplier

Let \(F_n,F_N\) be the normalized ordinary mirrored moment
polynomials for \(N=n+1\), and let \(E=z\partial_z\).
Write \(H_N=F_N-zF_N'/N\), as in the
[polar-half reduction](MULTIVARIATE_MIRRORED_POLAR_HALF_ANGULAR_REDUCTION.md).
The exact multiplier taking \(F_n\) to \(H_N\) splits into
a symmetric Gamma-shape advance and one decreasing Euler
factor for each upper Gegenbauer factor. The shape advance
stays within the family of simple unit-circle polynomials,
and the subsequent decreasing factors keep all zeros
outside the closed disk. An interlacing obstruction below
shows what these two separate stability facts leave open.

## Factorization and the attached shapes

For one mirrored pair
\((\ell,\ell+d),(\ell+d,\ell)\), where \(\ell\ge0\) and \(d\ge1\)
are integers, define
\[
       \alpha_h=\frac{d}{n\ell+h},\qquad
       \beta_h=\frac{d}{n(\ell+d)+h}.
\]
The exact polar multiplier of the pair at index \(k\) is
\[
 \begin{split}
 M_{\ell,d}(k)
 &=\prod_{h=1}^{\ell}(1+\alpha_h k)
                          \prod_{h=1}^{\ell+d}(1-\beta_h k)\\
 &=\prod_{h=1}^{\ell}\bigl(1+\alpha_h\beta_h k(n-k)\bigr)
                          \prod_{j=1}^d(1-\beta_{\ell+j} k).
 \end{split}                                                    \tag{1}
\]
In the factorial ratio \(e^{g_0(k/n)}\) defining the polar
multiplier, a row starting at \(a\) and having signed gap \(b-a\)
contributes
\[
 \frac{(na+(b-a)k+1)_a}{(na+1)_a}
        =\prod_{h=1}^a\left(1+\frac{(b-a)k}{na+h}\right).
\]
The endpoint argument normalization is identically one by
mirror symmetry. This proves the first line of (1),
and \(\alpha_h-\beta_h=n\alpha_h\beta_h\) proves the second.
Products over mirrored pairs give the full multiplier.

In operator form, the \(h\)-th symmetric factor
\(1+\alpha_h\beta_h E(n-E)\) advances both Gamma shapes
from \(h\) to \(h+1\) at the fixed degree \(n\).
This follows by multiplying each normalized Gamma
coefficient by its matching normalized rising factor.
Thus the symmetric stage, taken successively in \(h\)
in each pair, keeps the polynomial an exact mirrored Gamma
moment with simple unit-circle roots.
When that stage ends, each pair has both shapes \(\ell+1\).
Factoring its rising symbols into residues modulo \(d\)
then gives \(d\) normalized Gegenbauer factors with shapes
\[
                         \lambda_{N,j}=\frac{N\ell+j}{d},
                                \qquad 1\le j\le d.
\]
These are precisely the shapes of the upper degree factors.
This description uses the
[fixed-degree circular composition theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)
for the exact Gamma family.

The remaining \(j\)-th factor in (1) is
\[
                         1-\frac{E}{n+\lambda_{N,j}}.            \tag{2}
\]
It is itself the exact multiplier from the normalized
degree-\(n\) Gegenbauer factor of shape \(\lambda_{N,j}\)
to the polar half of the degree-\(N\) factor of the same
shape. Indeed their binomial-normalized coefficient
ratio is \((\lambda_{N,j}+n-k)/(\lambda_{N,j}+n)\).
Equivalently, if the intermediate reciprocal polynomial
after the symmetric stage is denoted \(G_n\), then the lower
moment array uses the lower shapes, \(G_n\) uses the upper
shapes at the lower degree, and \(H_N\) is obtained by
applying (2) for each shape attached to \(G_n\).
The private case \(\ell=0\) has only the residual stage.

The multiplier in (1) is positive and strictly decreasing
and log-concave on the real interval \(0\le k\le n\).
Positivity and strict log-concavity follow from the
linear factors, whose roots are outside this interval.
At \(k=0\) the log derivative is
\(d(R_n(\ell)-R_n(\ell+d))\), where \(R_n(0)=0\) and
\[
 R_n(t)=\sum_{h=1}^t\frac1{nt+h}
       =\int_0^\infty e^{-nx}
                    \frac{1-e^{-x}}{t(e^{x/t}-1)}\,dx
                         \quad(t\ge1).
\]
The integral follows by summing a geometric series after
scaling the integration variable by \(t\). It is
strictly increasing in \(t\), since
\((x/t)/(e^{x/t}-1)\) increases with \(t\).
Strict log-concavity then gives a negative log derivative
throughout the interval.

## What the residual stability gives

Let \(P\) have degree \(n\), with all zeros outside the open
unit disk and any unit-circle zeros simple.
For \(0<\beta<1/n\), the polynomial \(P-\beta EP\) is
nonvanishing on the closed disk.
For \(|z|\le1\) away from a zero of \(P\), each root
\(\xi\) contributes \(1-(1-z/\xi)^{-1}\) to \(EP/P\),
with real part at most \(1/2\).
Thus the real part of \(1-\beta EP/P\) is positive.
At a simple circle root of \(P\),
\(P-\beta EP=-\beta zP'\ne0\).
Every factor (2) satisfies \(0<\beta<1/n\), so in any
ordering the residual chain has exterior-disk zeros
after its first step. It also retains positive coefficients.

For a single factor the lifted polynomial
\(P-\beta EP+z^{n+1}(P-\beta EP)(1/z)\)
strictly interlaces a reciprocal \(P\) with simple circle
roots. At its lower roots the ratio of polar halves is
the constant \(n\beta>0\), and the
[lower-root polar sign](MULTIVARIATE_MIRRORED_POLAR_HALF_ANGULAR_REDUCTION.md#an-exact-lower-root-test)
applies. Already two factors require more information
about the root configuration.

## A positive sign theorem on one Gegenbauer core

Let \(G\) be the normalized degree-\(n\) Gegenbauer factor
of shape \(\lambda>0\). For \(r=2\) or \(3\) put
\[
 H=\prod_{j=1}^r(1-\beta_j E)G,\qquad
 F(z)=H(z)+z^{n+1}H(1/z),\qquad 0<\beta_j<1/n.
\]
If \(r=2\), \(F\) strictly interlaces \(G\).
For \(r=3\) the same holds when
\[
       A_1+A_2+A_3\ge\lambda+\tfrac12,\qquad
                      A_j=\beta_j^{-1}-n/2 .
\]
The condition holds, for instance, when all three
\(\beta_j=1/(n+\lambda)\).

To prove this, put \(t=\theta/2\),
\(g(\theta)=e^{-in\theta/2}G(e^{i\theta})\),
and \(q=n(n+2\lambda)/4\).
The [Gegenbauer equation](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)
gives \(g''+\lambda\cot t\,g'+qg=0\).
Hence at a zero of \(g\),
\[
 \frac{g''}{g'}=-\lambda\cot t,\qquad
 \frac{g'''}{g'}=\lambda^2\cot^2t+\frac{\lambda}{2}\csc^2t-q .
\]
In centered variables the residual product is
\((\prod_j\beta_j)\prod_j(A_j+i\partial_\theta)\).
Divide the centered lift at a lower zero by
\(2(\prod_j\beta_j)g'\sin t\), which extracts a
positive factor from the ratio to \(g'\).
For two factors the resulting value is
\(A_1+A_2+\lambda\cot^2t>0\).
For three it is
\[
 q-\frac{\lambda}{2}+\sum_{i<j}A_iA_j
    +\left(\lambda\sum_i A_i-\lambda^2-\frac{\lambda}{2}\right)
                                                    \cot^2t>0 .
\]
Here \(q-\lambda/2=n^2/4+(n-1)\lambda/2>0\).
The sign criterion proves both claims.

This specialization is a theorem for one input Gegenbauer
factor under two or three Euler operators. In an ordinary
profile the number \(r=\sum d_i\) of native residual factors
also counts the factors composing \(G_n\) in the
symmetric stage. For \(r\ge2\) and \(n\ge r+2\), that
composition cannot itself equal a single positive-shape
Gegenbauer factor. Indeed equality would give, for
\(k=0,\ldots,n-1\),
\[
       \prod_{j=1}^r\frac{\lambda_j+k}{\lambda_j+n-k-1}
                         =\frac{\lambda+k}{\lambda+n-k-1}.
\]
After clearing denominators, the degree-\((r+1)\)
polynomials agree at \(n>r+1\) points. They cannot
agree identically since their numbers of negative
roots with multiplicity are \(r\) and \(1\), respectively.
The specialization therefore does not give a new
all-degree ordinary profile at \(r=2\) or \(r=3\).
The composition stage would need a replacement for the
single Gegenbauer force identity.

## Two decreasing factors can fail the sign

Take \(n=6\) and
\[
               G(z)=\prod_{c\in\{1/2,\,1,\,101/100\}}
                                      (1+cz+z^2).
\]
Its roots are simple and circular. Its coefficients are
positive and its binomial-normalized coefficients are
strictly log-convex. Indeed the normalized first,
second, and third coefficients are
\(251/600,\ 1003/3000,\ 221/800\), respectively.
They satisfy
\[
    (251/600)^2<1003/3000,\qquad
    (1003/3000)^2<(251/600)(221/800),\qquad
    221/800<1003/3000 ,
\]
and reciprocity covers the remaining inequalities.

Let \(H=(1-E/7)^2 G\) and
\(F(z)=H(z)+z^7H(1/z)\). The residual argument shows
that \(H\) is disk-zero-free with positive coefficients.
The lifted \(F\) has simple unit-circle roots.
Indeed for its centered phase, the degree-six stable half
gives a phase derivative greater than \(7/2-6/2=1/2\),
and the phase goes from \(0\) to \(7\pi\).

Put \(g(\theta)=e^{-3i\theta}G(e^{i\theta})\), which is
\(\prod_c(2\cos\theta+c)\).
More generally, for
\(H=(1-\beta_1 E)(1-\beta_2 E)G\) and \(a_j=1-n\beta_j/2\),
at a root \(\theta\) of \(g\) the centered lift satisfies
\[
 \frac{e^{-i(n+1)\theta/2}F(e^{i\theta})}{2g'(\theta)}
     =(a_1\beta_2+a_2\beta_1)\sin(\theta/2)
               -\beta_1\beta_2\,\frac{g''(\theta)}{g'(\theta)}
                                      \cos(\theta/2).
\]
This follows from \(E=n/2-i\partial_\theta\) on the
centered function. At the root with \(\cos\theta=-101/200\),
\[
       \frac{g''(\theta)}{g'(\theta)}
              =\cot\theta+(400+400/51)\sin\theta>319.
\]
Here \(\sin\theta>4/5\) and \(\cot\theta>-1\).
But with \(\beta_1=\beta_2=1/7\) the positive sign would
require \(g''/g'<8\tan(\theta/2)<16\), since
\(\tan^2(\theta/2)=301/99<4\).
The sign is strictly negative, which rules out strict
interlacing by the lower-root test.

This counterexample uses an exterior-disk residual chain
and positive log-convex coefficients. Its initial \(G\)
is a general reciprocal polynomial, so it does not supply
an ordinary factorial counterexample. For the factorial
problem the exact remaining condition must relate the
residual factors to the attached Gegenbauer shape array,
or retain equivalent angular information about that array.
